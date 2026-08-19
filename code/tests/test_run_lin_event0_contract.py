"""Contract tests for the first-Lin-event C15--CLIMADA runner."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from run_lin_event0_c15_climada import (  # noqa: E402
    CENTRAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA,
    ENVIRONMENTAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA,
    MPS_TO_KNOTS,
    build_climada_track,
    build_rainfall_output,
)
from run_irene_c15_climada import (  # noqa: E402
    MAX_DISTANCE_EYE_KM,
    build_moving_union_grid,
    shortest_lon_delta_deg,
    unwrap_longitude_deg,
    wrap_lon_deg,
)


class LinEventRunnerContractTest(unittest.TestCase):
    def setUp(self):
        count = 94
        self.prepared = xr.Dataset(
            data_vars={
                "native_index": ("time", np.arange(56, 150, dtype=np.int32)),
                "time_seconds_from_seed": ("time", np.arange(count) * 3600.0),
                "lat": ("time", np.linspace(18.6, 26.2, count)),
                "lon": ("time", np.linspace(273.9, 260.0, count)),
                "circular_wind": ("time", np.linspace(12.0, 28.0, count)),
                "radius_max_wind": ("time", np.linspace(62.0, 101.0, count)),
                "q950": ("time", np.full(count, 0.017)),
                "ushear": ("time", np.full(count, 2.0)),
                "vshear": ("time", np.full(count, 3.0)),
            },
            attrs={
                "event_id": "stream0000-year1995-track000002",
                "outer_radius_m": 1_000_000.0,
                "outer_radius_fixed_for_event_lifetime": 1,
            },
        )
        self.manifest = {
            "event": {
                "threshold_genesis_datetime": "1995-09-17 08:00:00",
                "threshold_genesis_region": "NA",
            }
        }

    def test_track_uses_circular_wind_and_c15_r0input_rmw(self):
        track = build_climada_track(self.prepared, self.manifest)
        np.testing.assert_allclose(
            track["max_sustained_wind"],
            self.prepared["circular_wind"] * MPS_TO_KNOTS,
        )
        np.testing.assert_allclose(
            track["radius_max_wind"], self.prepared["radius_max_wind"] / 1.852
        )
        np.testing.assert_allclose(
            track["central_pressure"], CENTRAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA
        )
        np.testing.assert_allclose(
            track["environmental_pressure"],
            ENVIRONMENTAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA,
        )
        self.assertEqual(track.attrs["schema_pressure_placeholders_unused_by_c15_tcr"], 1)
        self.assertEqual(track.attrs["outer_radius_m"], 1_000_000.0)
        self.assertNotIn("size_predictor_vmax", track.variables)

    def test_hourly_accumulation_semantics(self):
        track = build_climada_track(self.prepared, self.manifest)
        grid = xr.Dataset(
            data_vars={
                "lat": ("centroid", [20.0, 21.0]),
                "lon": ("centroid", [260.0, 261.0]),
            },
            coords={"centroid": [0, 1]},
        )
        rates = np.column_stack(
            [np.arange(1, 95, dtype=float), np.ones(94, dtype=float)]
        )
        output = build_rainfall_output(track, grid, rates, {})
        self.assertEqual(float(output["event_total_rainfall"][1]), 94.0)
        self.assertEqual(float(output["maximum_24h_rainfall"][1]), 24.0)
        self.assertEqual(
            float(output["maximum_24h_rainfall"][0]), float(np.sum(np.arange(71, 95)))
        )

    def test_event_shorter_than_24h_uses_whole_available_window(self):
        prepared = self.prepared.isel(time=slice(0, 7)).copy()
        manifest = dict(self.manifest)
        track = build_climada_track(prepared, manifest)
        grid = xr.Dataset(
            data_vars={"lat": ("centroid", [20.0]), "lon": ("centroid", [260.0])},
            coords={"centroid": [0]},
        )
        output = build_rainfall_output(
            track, grid, np.arange(1, 8, dtype=float)[:, None], {}
        )
        self.assertEqual(float(output["event_total_rainfall"][0]), 28.0)
        self.assertEqual(float(output["maximum_24h_rainfall"][0]), 28.0)
        self.assertEqual(output.attrs["maximum_24h_window_hours_used"], 7)
        self.assertEqual(output.attrs["event_shorter_than_24h"], 1)


class MovingUnionGridPeriodicLongitudeTest(unittest.TestCase):
    def _track(self, lat, lon):
        return xr.Dataset(
            data_vars={
                "lat": ("time", np.asarray(lat, dtype=float)),
                "lon": ("time", np.asarray(lon, dtype=float)),
            },
            coords={"time": np.arange(len(lat))},
        )

    def test_unwrap_crosses_the_date_line_by_the_short_arc(self):
        unwrapped = unwrap_longitude_deg(np.asarray([179.0, -179.0, -178.0]))
        self.assertAlmostEqual(unwrapped[0], 179.0)
        self.assertAlmostEqual(unwrapped[1], 181.0)
        self.assertAlmostEqual(unwrapped[2], 182.0)
        self.assertAlmostEqual(float(shortest_lon_delta_deg(-179.0, 179.0)), 2.0)
        self.assertAlmostEqual(float(wrap_lon_deg(181.0)), -179.0)

    def test_non_crossing_track_still_builds_a_300km_disk(self):
        grid = build_moving_union_grid(self._track([20.0, 20.1], [140.0, 140.2]))
        self.assertGreater(grid["centroid"].size, 0)
        self.assertLessEqual(float(grid["lon"].max()), 180.0)
        self.assertGreaterEqual(float(grid["lon"].min()), -180.0)

    def test_antimeridian_track_is_kept_and_covers_both_sides(self):
        grid = build_moving_union_grid(self._track([10.0, 10.0], [179.5, -179.5]))
        lon = np.asarray(grid["lon"].values, dtype=float)
        self.assertGreater(lon.size, 0)
        self.assertTrue(np.any(lon > 170.0))
        self.assertTrue(np.any(lon < -170.0))
        # Every centroid is within 300 km of at least one hourly eye.
        d0 = 2.0 * 6371.0 * np.arcsin(np.sqrt(np.clip(
            np.cos(np.radians(10.0)) ** 2
            * np.sin(np.radians(shortest_lon_delta_deg(lon, 179.5)) / 2.0) ** 2,
            0.0, 1.0,
        )))
        d1 = 2.0 * 6371.0 * np.arcsin(np.sqrt(np.clip(
            np.cos(np.radians(10.0)) ** 2
            * np.sin(np.radians(shortest_lon_delta_deg(lon, -179.5)) / 2.0) ** 2,
            0.0, 1.0,
        )))
        nearest = np.minimum(d0, d1)
        self.assertTrue(np.all(nearest <= MAX_DISTANCE_EYE_KM + 1.0))


if __name__ == "__main__":
    unittest.main()
