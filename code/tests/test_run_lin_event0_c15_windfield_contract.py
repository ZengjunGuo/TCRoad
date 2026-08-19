"""Narrow contract tests for the first-Lin-event C15 wind-field runner."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from run_lin_event0_c15_windfield import (  # noqa: E402
    BACKGROUND_CCW_ROTATION_DEG,
    BACKGROUND_REDUCTION_FACTOR,
    CLIMADA_MAJORITY_HEMISPHERE_RULE,
    climada_majority_hemisphere_sign,
    compute_wind_field,
    lin_chavas_background_wind,
    spherical_distance_and_outward_bearing,
)


class _LinearProfileProvider:
    """Small deterministic profile seam; it does not stand in for production."""

    parameters = {}
    cache_info = (0, 1, 1)
    outer_radius_m = 400_000.0

    def profile_for(self, vmax_ms, coriolis_s):
        del vmax_ms, coriolis_s
        return SimpleNamespace(
            radius_m=np.asarray([0.0, 200_000.0, 400_000.0]),
            gradient_wind_ms=np.asarray([0.0, 20.0, 0.0]),
            outer_radius_m=400_000.0,
            radius_max_wind_m=50_000.0,
        )


def _prepared(*, translation_u: float = 0.0, translation_v: float = 0.0) -> xr.Dataset:
    return xr.Dataset(
        data_vars={
            "lat": ("time", [20.0]),
            "lon": ("time", [0.0]),
            "circular_wind": ("time", [30.0]),
            "radius_max_wind": ("time", [50.0]),
            "translation_u": ("time", [translation_u]),
            "translation_v": ("time", [translation_v]),
        },
        coords={"time": [0]},
        attrs={
            "outer_radius_m": 400_000.0,
            "outer_radius_fixed_for_event_lifetime": 1,
        },
    )


def _track(ntime: int = 1) -> xr.Dataset:
    times = np.asarray(
        [
            np.datetime64("1995-09-17T08:00:00") + np.timedelta64(hour, "h")
            for hour in range(ntime)
        ]
    )
    return xr.Dataset(coords={"time": times})


def _prepared_lats(lats: list[float], *, translation_u: float = 0.0) -> xr.Dataset:
    ntime = len(lats)
    return xr.Dataset(
        data_vars={
            "lat": ("time", np.asarray(lats, dtype=float)),
            "lon": ("time", np.zeros(ntime, dtype=float)),
            "circular_wind": ("time", np.full(ntime, 30.0)),
            "radius_max_wind": ("time", np.full(ntime, 50.0)),
            "translation_u": ("time", np.full(ntime, translation_u)),
            "translation_v": ("time", np.zeros(ntime)),
        },
        coords={"time": np.arange(ntime)},
        attrs={
            "outer_radius_m": 400_000.0,
            "outer_radius_fixed_for_event_lifetime": 1,
        },
    )


class LinEventWindfieldContractTest(unittest.TestCase):
    def test_spherical_cardinal_bearings(self):
        distance, bearing = spherical_distance_and_outward_bearing(
            20.0,
            0.0,
            np.asarray([21.0, 20.0]),
            np.asarray([0.0, 1.0]),
        )
        self.assertAlmostEqual(distance[0] / 1000.0, 111.195, places=3)
        self.assertAlmostEqual(bearing[0], 0.0, places=12)
        self.assertAlmostEqual(np.rad2deg(bearing[1]), 90.171, places=3)

    def test_oblique_bearing_is_point_local_final_not_center_initial(self):
        # This oblique, non-equatorial arc curves enough that its center-local
        # initial and point-local outward bearings differ by about 30 degrees.
        _, outward = spherical_distance_and_outward_bearing(
            35.0,
            -75.0,
            np.asarray([45.0]),
            np.asarray([-30.0]),
        )
        self.assertAlmostEqual(np.rad2deg(outward[0]), 89.604709, places=6)
        center_initial_bearing_deg = 59.677510
        self.assertGreater(
            abs(np.rad2deg(outward[0]) - center_initial_bearing_deg), 29.0
        )

    def test_lin_chavas_background_is_uniform_vector_transform(self):
        background_u, background_v = lin_chavas_background_wind(
            np.asarray([10.0, 10.0]), np.asarray([0.0, 0.0])
        )
        expected_magnitude = BACKGROUND_REDUCTION_FACTOR * 10.0
        expected_angle = np.deg2rad(BACKGROUND_CCW_ROTATION_DEG)
        np.testing.assert_allclose(
            background_u, expected_magnitude * np.cos(expected_angle)
        )
        np.testing.assert_allclose(
            background_v, expected_magnitude * np.sin(expected_angle)
        )
        self.assertAlmostEqual(
            float(np.hypot(background_u[0], background_v[0])),
            expected_magnitude,
            places=12,
        )

    def test_nh_tangent_and_direct_linear_c15_profile(self):
        # North and east points have nearly equal radius.  For a NH cyclone,
        # the vortex points west at the north point and north at the east point.
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", [21.0, 20.0]),
                "lon": ("centroid", [0.0, 1.0]),
            },
            coords={"centroid": [0, 1]},
        )
        output, _ = compute_wind_field(
            _prepared(), _track(), grid, provider=_LinearProfileProvider()
        )
        north_u = float(output["near_surface_wind_u"][0, 0])
        north_v = float(output["near_surface_wind_v"][0, 0])
        east_u = float(output["near_surface_wind_u"][0, 1])
        east_v = float(output["near_surface_wind_v"][0, 1])
        self.assertLess(north_u, 0.0)
        self.assertAlmostEqual(north_v, 0.0, places=5)
        self.assertAlmostEqual(east_u, 0.0, delta=0.04)
        self.assertGreater(east_v, 0.0)

        north_distance, _ = spherical_distance_and_outward_bearing(
            20.0, 0.0, np.asarray([21.0]), np.asarray([0.0])
        )
        expected = np.interp(
            north_distance[0],
            [0.0, 200_000.0, 400_000.0],
            [0.0, 20.0, 0.0],
        )
        self.assertAlmostEqual(-north_u, expected, places=5)

    def test_sh_cyclonic_tangent_and_background_rotation_reverse(self):
        prepared = _prepared(translation_u=10.0)
        prepared["lat"][:] = -20.0
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", [-19.0]),
                "lon": ("centroid", [0.0]),
            },
            coords={"centroid": [0]},
        )
        output, metadata = compute_wind_field(
            prepared, _track(), grid, provider=_LinearProfileProvider()
        )
        # At the point north of an SH center, clockwise cyclonic flow is east.
        self.assertGreater(float(output["near_surface_wind_u"][0, 0]), 0.0)
        self.assertLess(float(output["near_surface_wind_v"][0, 0]), 0.0)
        self.assertEqual(metadata["hemisphere_sign"], -1.0)
        self.assertEqual(metadata["cyclonic_background_rotation_degrees"], -20.0)

        background_u, background_v = lin_chavas_background_wind(
            np.asarray([10.0]),
            np.asarray([0.0]),
            hemisphere_sign=-1.0,
        )
        self.assertGreater(float(background_u[0]), 0.0)
        self.assertLess(float(background_v[0]), 0.0)

    def test_oblique_vortex_has_no_point_local_radial_leakage(self):
        point_lat = np.asarray([21.5])
        point_lon = np.asarray([1.5])
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", point_lat),
                "lon": ("centroid", point_lon),
            },
            coords={"centroid": [0]},
        )
        output, _ = compute_wind_field(
            _prepared(), _track(), grid, provider=_LinearProfileProvider()
        )
        _, outward = spherical_distance_and_outward_bearing(
            20.0, 0.0, point_lat, point_lon
        )
        radial_east = np.sin(outward[0])
        radial_north = np.cos(outward[0])
        wind_u = float(output["near_surface_wind_u"][0, 0])
        wind_v = float(output["near_surface_wind_v"][0, 0])
        radial_leakage = wind_u * radial_east + wind_v * radial_north
        self.assertAlmostEqual(radial_leakage, 0.0, places=6)

    def test_no_averaging_conversion_and_event_maximum(self):
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", [21.0]),
                "lon": ("centroid", [0.0]),
            },
            coords={"centroid": [0]},
        )
        output, metadata = compute_wind_field(
            _prepared(translation_u=10.0),
            _track(),
            grid,
            provider=_LinearProfileProvider(),
        )
        self.assertEqual(output.attrs["ten_minute_sustained_wind_claim"], 0)
        self.assertEqual(output.attrs["wind_averaging_period_conversion_applied"], 0)
        self.assertIn("unspecified", output.attrs["wind_averaging_period"])
        self.assertEqual(
            output.attrs["surface_wind_reduction_factor_applied_to_c15"], 0
        )
        self.assertEqual(output.attrs["inflow_angle_applied_to_c15"], 0)
        self.assertEqual(output.attrs["radial_background_decay_or_taper_applied"], 0)
        self.assertAlmostEqual(
            float(output["event_maximum_near_surface_wind_speed"][0]),
            float(output["near_surface_wind_speed"][0, 0]),
            places=6,
        )
        self.assertTrue(metadata["linear_profile_interpolation"])

    def test_majority_hemisphere_sign_matches_climada_61(self):
        sign, north, south = climada_majority_hemisphere_sign(np.asarray([20.0]))
        self.assertEqual((sign, north, south), (1.0, 1, 0))
        sign, north, south = climada_majority_hemisphere_sign(np.asarray([-20.0]))
        self.assertEqual((sign, north, south), (-1.0, 0, 1))

        # Event 14472 analog: 22 N then 227 S, first node is northern.
        crossing = np.concatenate(
            [np.full(22, 3.80), np.full(227, -7.37)]
        )
        sign, north, south = climada_majority_hemisphere_sign(crossing)
        self.assertEqual((sign, north, south), (-1.0, 22, 227))

        # Tie, equator nodes ignored, all-equator: CLIMADA defaults to N.
        self.assertEqual(
            climada_majority_hemisphere_sign(np.asarray([1.0, -1.0])),
            (1.0, 1, 1),
        )
        self.assertEqual(
            climada_majority_hemisphere_sign(np.asarray([1.0, 0.0, -1.0, -2.0])),
            (-1.0, 1, 2),
        )
        self.assertEqual(
            climada_majority_hemisphere_sign(np.asarray([0.0, 0.0])),
            (1.0, 0, 0),
        )

    def test_equator_crossing_majority_south_produces_sh_wind_field(self):
        # Cheap 14472-shaped analog: first node NH, southern nodes win.
        lats = [3.80, 2.00] + [-7.00] * 5
        prepared = _prepared_lats(lats)
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", [3.80, -6.00]),
                "lon": ("centroid", [0.0, 0.0]),
            },
            coords={"centroid": [0, 1]},
        )
        output, metadata = compute_wind_field(
            prepared, _track(ntime=len(lats)), grid, provider=_LinearProfileProvider()
        )
        self.assertEqual(metadata["hemisphere_sign"], -1.0)
        self.assertEqual(metadata["hemisphere_northern_node_count"], 2)
        self.assertEqual(metadata["hemisphere_southern_node_count"], 5)
        self.assertEqual(
            metadata["hemisphere_rule"], CLIMADA_MAJORITY_HEMISPHERE_RULE
        )
        self.assertEqual(metadata["cyclonic_background_rotation_degrees"], -20.0)
        self.assertTrue(np.all(np.isfinite(output["event_maximum_near_surface_wind_speed"])))
        # Point due north of an SH center: clockwise cyclonic flow is east.
        sh_hour = 2
        self.assertGreater(float(output["near_surface_wind_u"][sh_hour, 1]), 0.0)
        self.assertAlmostEqual(float(output["near_surface_wind_v"][sh_hour, 1]), 0.0, places=5)


if __name__ == "__main__":
    unittest.main()
