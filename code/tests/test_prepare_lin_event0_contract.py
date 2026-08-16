"""Narrow tests for the Lin-event synthetic-environment prepare candidate."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from netCDF4 import Dataset, date2num
import numpy as np


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from prepare_lin_event0_c15_climada import (  # noqa: E402
    KNOTS_TO_MPS,
    MPS_TO_KNOTS,
    derive_c15_r0input_rmw,
    emanuel_v64_qs950,
    emanuel_v64_shear,
    load_fixed_r0_for_event,
    official_utrans_adapted_to_native_hour,
    sample_t600_lin_public_adapter,
)


class EmanuelPublicAdapterContractTest(unittest.TestCase):
    def test_qs900b_known_vector(self):
        t600 = np.array([275.5, 276.0, 277.0])
        circular = np.array([12.0, 20.0, 30.0])
        actual = emanuel_v64_qs950(t600, circular)
        expected = np.array(
            [0.016297854910440654, 0.01703411817375113, 0.01853850366742712]
        )
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-15)

    def test_raingen_shear_units_and_drift(self):
        utrans_knots = np.array([-10.0, 4.0])
        vtrans_knots = np.array([2.0, 8.0])
        u850_ms = np.array([-2.0, 1.0])
        v850_ms = np.array([0.5, -1.0])
        lat = np.array([20.0, 30.0])
        ushear, vshear = emanuel_v64_shear(
            utrans_knots, vtrans_knots, u850_ms, v850_ms, lat, 12.0
        )
        expected_u = 5.0 * (utrans_knots * KNOTS_TO_MPS - u850_ms)
        expected_v = 5.0 * (
            vtrans_knots * KNOTS_TO_MPS
            - 1.5 * np.cos(np.deg2rad(lat))
            - v850_ms
        )
        np.testing.assert_allclose(ushear, expected_u, rtol=0, atol=1e-14)
        np.testing.assert_allclose(vshear, expected_v, rtol=0, atol=1e-14)

    def test_utrans_one_hour_dimensional_adaptation(self):
        native = np.arange(0, 9, dtype=np.int64)
        time = native.astype(float) * 3600.0
        # Constant 0.1 degree/hour eastward and 0.05 degree/hour northward.
        track = {
            "lon_trks": 250.0 + 0.1 * native,
            "lat_trks": 20.0 + 0.05 * native,
            "native_index": native,
            "target_native_index": native[2:7],
            "time": time,
        }
        u_knots, v_knots = official_utrans_adapted_to_native_hour(track)
        expected_v = np.full(5, 3.0)
        expected_u = 6.0 * np.cos(np.deg2rad(track["lat_trks"][2:7]))
        np.testing.assert_allclose(v_knots, expected_v, rtol=0, atol=1e-12)
        np.testing.assert_allclose(u_knots, expected_u, rtol=0, atol=2e-5)
        # Guard against accidentally applying a second m/s-to-knot conversion.
        self.assertTrue(np.all(np.abs(u_knots) < 10 * MPS_TO_KNOTS))

    def test_utrans_official_endpoint_extrapolation_is_finite(self):
        native = np.arange(0, 5, dtype=np.int64)
        track = {
            "lon_trks": 100.0 + native,
            "lat_trks": 20.0 + 0.2 * native,
            "native_index": native,
            "target_native_index": native,
            "time": native.astype(float) * 3600.0,
        }
        u_knots, v_knots = official_utrans_adapted_to_native_hour(track)
        self.assertTrue(np.all(np.isfinite(u_knots)))
        self.assertTrue(np.all(np.isfinite(v_knots)))
        self.assertAlmostEqual(v_knots[0], 12.0, places=12)
        self.assertAlmostEqual(v_knots[-1], 12.0, places=12)

    def test_fixed_r0_provider_derives_every_hourly_rmw(self):
        class Provider:
            outer_radius_m = 800_000.0
            cache_info = (0, 3, 3)

            def profile_for(self, vmax_ms, coriolis_s):
                self.assert_positive = (vmax_ms > 0.0 and coriolis_s > 0.0)
                return type("Profile", (), {"radius_max_wind_m": vmax_ms * 1000.0})()

        intensity = np.asarray([12.0, 20.0, 30.0])
        rmw_km, metadata = derive_c15_r0input_rmw(
            intensity,
            np.asarray([18.0, 22.0, 26.0]),
            800_000.0,
            provider=Provider(),
        )
        np.testing.assert_allclose(rmw_km, intensity, rtol=0.0, atol=0.0)
        self.assertEqual(metadata["c15_input_mode"], "r0input")
        self.assertTrue(metadata["outer_radius_fixed_for_event_lifetime"])
        self.assertFalse(metadata["clipping_applied"])

    def test_catalogue_requires_exact_position_and_event_id(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            catalogue = root / "fixed.nc"
            manifest = root / "fixed.manifest.json"
            with Dataset(catalogue, "w", format="NETCDF4") as output:
                output.createDimension("event", 1)
                output.createVariable("event_position", "i4", ("event",))[:] = [0]
                output.createVariable("event_id", str, ("event",))[:] = np.asarray(
                    ["e0"], dtype=object
                )
                output.createVariable("outer_radius_m", "f8", ("event",))[:] = [800_000.0]
                output.status = "FROZEN_IMMUTABLE"
                output.source_sample_sha256 = "sample-sha"
                output.distribution_contract_sha256 = "distribution-sha"
                output.event_outer_radius_binding_sha256 = "binding-sha"
            import hashlib, json

            digest = hashlib.sha256(catalogue.read_bytes()).hexdigest()
            manifest.write_text(
                json.dumps(
                    {
                        "status": "FROZEN_IMMUTABLE",
                        "source_sample": {"sha256": "sample-sha"},
                        "artifacts": {
                            "fixed_r0_catalogue_netcdf": {
                                "bytes": catalogue.stat().st_size,
                                "sha256": digest,
                            }
                        },
                    }
                )
            )
            radius, _ = load_fixed_r0_for_event(
                catalogue,
                manifest,
                {"event_position": 0, "event_id": "e0"},
                expected_sample_sha256="sample-sha",
            )
            self.assertEqual(radius, 800_000.0)
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                load_fixed_r0_for_event(
                    catalogue,
                    manifest,
                    {"event_position": 0, "event_id": "wrong"},
                    expected_sample_sha256="sample-sha",
                )

    def test_january_boundary_and_periodic_longitude_are_mechanical(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ta.nc"
            with Dataset(path, "w", format="NETCDF4") as source:
                source.createDimension("time", 2)
                source.createDimension("plev", 1)
                source.createDimension("lat", 2)
                source.createDimension("lon", 4)
                time = source.createVariable("time", "f8", ("time",))
                time.units = "days since 1850-01-01 00:00:00"
                time.calendar = "proleptic_gregorian"
                time[:] = date2num(
                    [
                        datetime(1995, 1, 16, 12),
                        datetime(1995, 2, 15),
                    ],
                    time.units,
                    time.calendar,
                )
                plev = source.createVariable("plev", "f8", ("plev",))
                plev.units = "Pa"
                plev[:] = [60_000.0]
                source.createVariable("lat", "f8", ("lat",))[:] = [-10.0, 10.0]
                source.createVariable("lon", "f8", ("lon",))[:] = [0, 90, 180, 270]
                ta = source.createVariable(
                    "ta", "f8", ("time", "plev", "lat", "lon")
                )
                ta.units = "K"
                field = 275.0 + np.cos(np.deg2rad([0, 90, 180, 270]))
                ta[0, 0, :, :] = np.tile(field, (2, 1))
                ta[1, 0, :, :] = np.tile(field + 1.0, (2, 1))

            with patch(
                "prepare_lin_event0_c15_climada.validate_identity",
                return_value={"path": str(path), "bytes": path.stat().st_size, "sha256": "fixture"},
            ):
                result, metadata = sample_t600_lin_public_adapter(
                    path,
                    np.asarray([359.5, -0.5]),
                    np.asarray([0.0, 0.0]),
                    1995,
                    1,
                )
        self.assertTrue(metadata["boundary_month_field_used"])
        self.assertEqual(metadata["time_left_index"], 0)
        self.assertEqual(metadata["time_right_index"], 0)
        self.assertAlmostEqual(float(result[0]), float(result[1]), places=12)


if __name__ == "__main__":
    unittest.main()
