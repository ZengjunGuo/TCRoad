"""Contract tests for the official-C15 / CLIMADA-Petals TCR seam."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

import c15_climada_tcr as bridge  # noqa: E402
from c15_climada_tcr import (  # noqa: E402
    C15FixedR0WindProfileProvider,
    C15ProfileDomainError,
    C15WindProfileProvider,
    C15_WIND_MODEL_NAME,
    patched_c15_tcr,
)


def _si_track(
    *,
    ntime: int = 2,
    vmax_ms: float = 50.0,
    rmax_m: float = 25_000.0,
    coriolis_s: float = 5e-5,
) -> xr.Dataset:
    return xr.Dataset(
        data_vars={
            "vmax": ("time", np.full(ntime, vmax_ms)),
            "rad": ("time", np.full(ntime, rmax_m)),
            "cp": ("time", np.full(ntime, coriolis_s)),
        },
        coords={"time": np.arange(ntime)},
        attrs={"latsign": 1.0},
    )


class OfficialC15ProviderContractTest(unittest.TestCase):
    def test_official_anchor_cache_and_public_mu_companion(self):
        provider = C15WindProfileProvider()
        track = _si_track()
        radius_m = np.array(
            [[0.0, 25_000.0, 300_000.0], [0.0, 25_000.0, 300_000.0]]
        )
        active = np.ones_like(radius_m, dtype=bool)

        gradient = provider(track, radius_m, active, cyclostrophic=False)
        # The existing official-PDF test suite fixes this rmaxinput anchor.
        self.assertAlmostEqual(gradient[1, 2], 9.760671, delta=0.003)
        np.testing.assert_array_equal(gradient[0], 0.0)

        companion = provider(track, radius_m, active, cyclostrophic=True)
        f = 5e-5
        rm = 25_000.0
        vmax = 50.0
        query = radius_m[1]
        momentum = query * gradient[1] + 0.5 * f * query**2
        maximum_momentum = rm * vmax + 0.5 * f * rm**2
        expected = np.zeros_like(query)
        expected[1:] = rm * vmax * (momentum[1:] / maximum_momentum) / query[1:]
        np.testing.assert_allclose(companion[1], expected, rtol=1e-13, atol=1e-13)
        np.testing.assert_array_equal(companion[0], 0.0)

        hits, misses, size = provider.cache_info
        self.assertGreaterEqual(hits, 1)
        self.assertEqual((misses, size), (1, 1))

    def test_inactive_radius_is_never_queried_and_active_beyond_r0_errors(self):
        provider = C15WindProfileProvider()
        track = _si_track()
        profile = provider.profile_for(50.0, 25_000.0, 5e-5)
        beyond_r0 = profile.outer_radius_m + 1000.0
        radius_m = np.array([[0.0, beyond_r0], [300_000.0, beyond_r0]])

        # The point outside both 300 km and r0 is inactive, so it is not sent
        # to C15 and remains exactly zero.
        inactive_outside = np.array([[False, False], [True, False]])
        result = provider(track, radius_m, inactive_outside, cyclostrophic=False)
        self.assertEqual(result[1, 1], 0.0)
        self.assertGreater(result[1, 0], 0.0)

        active_outside = np.array([[False, False], [False, True]])
        with self.assertRaises(C15ProfileDomainError):
            provider(track, radius_m, active_outside, cyclostrophic=False)


class OfficialC15FixedR0ProviderContractTest(unittest.TestCase):
    def test_official_r0input_anchor_and_exact_finite_endpoint(self):
        provider = C15FixedR0WindProfileProvider(outer_radius_m=900_000.0)
        profile = provider.profile_for(50.0, 5e-5)

        # Official r0input PDF anchor, independently frozen in the adapter's
        # example-curve tests.
        actual_300km = np.interp(
            300_000.0, profile.radius_m, profile.gradient_wind_ms
        )
        self.assertAlmostEqual(actual_300km, 9.528959, delta=0.003)
        self.assertEqual(profile.outer_radius_m, 900_000.0)
        self.assertEqual(profile.radius_m[-1], 900_000.0)
        self.assertEqual(profile.gradient_wind_ms[-1], 0.0)
        self.assertAlmostEqual(
            profile.radius_m[np.argmax(profile.gradient_wind_ms)],
            profile.radius_max_wind_m,
            delta=1e-8,
        )

    def test_event_r0_is_constant_while_solver_infers_dynamic_rmax(self):
        calls: list[tuple[float, float, float]] = []

        def fake_solver(vmax_ms, r0_m, coriolis_s, parameters):
            del parameters
            calls.append((vmax_ms, r0_m, coriolis_s))
            inferred_rmax_m = 10_000.0 + 500.0 * vmax_ms
            return (
                np.array([0.0, inferred_rmax_m, r0_m - 1000.0]),
                np.array([0.0, vmax_ms, 0.1]),
                inferred_rmax_m,
            )

        provider = C15FixedR0WindProfileProvider(
            outer_radius_m=880_000.0, solver=fake_solver
        )
        weak = provider.profile_for(40.0, 4e-5)
        strong = provider.profile_for(60.0, 4e-5)

        self.assertEqual(weak.outer_radius_m, strong.outer_radius_m)
        self.assertEqual(weak.outer_radius_m, 880_000.0)
        self.assertNotEqual(
            weak.radius_max_wind_m, strong.radius_max_wind_m
        )
        self.assertEqual([call[1] for call in calls], [880_000.0, 880_000.0])

    def test_active_query_beyond_fixed_r0_is_rejected_without_tail(self):
        def fake_solver(vmax_ms, r0_m, coriolis_s, parameters):
            del coriolis_s, parameters
            return (
                np.array([0.0, 25_000.0, r0_m - 1000.0]),
                np.array([0.0, vmax_ms, 0.1]),
                25_000.0,
            )

        provider = C15FixedR0WindProfileProvider(
            outer_radius_m=300_000.0, solver=fake_solver
        )
        track = _si_track()
        radius_m = np.array(
            [[0.0, 0.0], [300_000.0, 300_000.001]], dtype=float
        )
        active = np.ones_like(radius_m, dtype=bool)
        with self.assertRaisesRegex(C15ProfileDomainError, "no|fixed finite|exceeds"):
            provider(track, radius_m, active, cyclostrophic=False)

    def test_official_wrapper_uses_fifth_return_value_as_rmax(self):
        adapter = mock.Mock()
        adapter.ER11E04_nondim_r0input.return_value = (
            np.array([0.0, 1.0]),
            np.array([0.0, 2.0]),
            333.0,  # rmerge: must not be exposed as RMW
            444.0,  # Vmerge: must not be exposed as RMW
            555.0,  # rmax
        )
        with mock.patch.object(
            bridge, "_load_frozen_c15_adapter", return_value=adapter
        ):
            radius_m, wind_ms, rmax_m = bridge._official_r0input_solver(
                50.0, 900_000.0, 5e-5, bridge.C15_R0INPUT_DEFAULTS
            )
        np.testing.assert_array_equal(radius_m, [0.0, 1.0])
        np.testing.assert_array_equal(wind_ms, [0.0, 2.0])
        self.assertEqual(rmax_m, 555.0)

    def test_official_wrapper_converts_system_exit_to_diagnostic_error(self):
        adapter = mock.Mock()
        adapter.ER11E04_nondim_r0input.side_effect = SystemExit(7)
        with mock.patch.object(
            bridge, "_load_frozen_c15_adapter", return_value=adapter
        ):
            with self.assertRaisesRegex(
                RuntimeError, r"r0input aborted.*r0=900000.*exit=7"
            ):
                bridge._official_r0input_solver(
                    50.0, 900_000.0, 5e-5, bridge.C15_R0INPUT_DEFAULTS
                )


class PatchedPetalsContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            from climada_petals.hazard import tc_rainfield
        except ImportError as error:  # pragma: no cover - environment routing
            raise unittest.SkipTest(f"CLIMADA Petals is not installed: {error}")
        cls.tc_rainfield = tc_rainfield

    def test_c15_tcr_dispatch_calls_provider_not_er11_and_restores(self):
        tc_rainfield = self.tc_rainfield
        original = tc_rainfield._windprofile
        self.assertNotIn(C15_WIND_MODEL_NAME, tc_rainfield.MODEL_VANG)

        calls: list[bool] = []

        def fake_provider(si_track, radius_m, active_mask, *, cyclostrophic):
            calls.append(cyclostrophic)
            return np.full_like(radius_m, 7.0, dtype=float)

        track = _si_track(ntime=3)
        distance = np.full((3, 2), 100_000.0)
        distances = {
            "": distance,
            "+": distance + 2000.0,
            "-": distance - 2000.0,
            "dir": np.zeros((3, 2, 2)),
        }
        active = np.ones((3, 2), dtype=bool)

        # If the C15 dispatch leaks to CLIMADA's default path, this sentinel
        # turns the test into an immediate failure.
        with mock.patch.object(
            tc_rainfield,
            "compute_angular_windspeeds",
            side_effect=AssertionError("ER11/default angular-wind path was called"),
        ):
            with patched_c15_tcr(fake_provider):
                model = tc_rainfield.MODEL_VANG[C15_WIND_MODEL_NAME]
                winds = tc_rainfield._horizontal_winds(
                    track, distances, active, model
                )

        self.assertIs(tc_rainfield._windprofile, original)
        self.assertNotIn(C15_WIND_MODEL_NAME, tc_rainfield.MODEL_VANG)
        self.assertIn(False, calls)
        self.assertIn(True, calls)
        np.testing.assert_array_equal(winds["r,t"][0], 0.0)
        np.testing.assert_array_equal(winds["nocoriolis"][0], 0.0)
        np.testing.assert_array_equal(winds["r,t"][1:], 7.0)


if __name__ == "__main__":
    unittest.main()
