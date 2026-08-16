"""Non-science contract tests for the generic single-event worker."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from run_lin_event_worker import build_compact_hazard_footprint  # noqa: E402


class GenericEventWorkerContractTest(unittest.TestCase):
    def test_compact_footprint_preserves_spatial_fields_but_not_hourly_arrays(self):
        centroid = np.arange(3, dtype=np.int32)
        lat = np.asarray([-20.0, -20.1, -20.2])
        lon = np.asarray([359.9, 0.0, 0.1])
        rainfall = xr.Dataset(
            data_vars={
                "rainfall_rate": (("time", "centroid"), np.ones((2, 3))),
                "event_total_rainfall": ("centroid", [2.0, 4.0, 6.0]),
                "maximum_24h_rainfall": ("centroid", [2.0, 4.0, 6.0]),
                "lat": ("centroid", lat),
                "lon": ("centroid", lon),
            },
            coords={"time": [0, 1], "centroid": centroid},
        )
        wind = xr.Dataset(
            data_vars={
                "near_surface_wind_speed": (
                    ("time", "centroid"),
                    np.ones((2, 3)),
                ),
                "event_maximum_near_surface_wind_speed": (
                    "centroid",
                    [10.0, 20.0, 30.0],
                ),
                "lat": ("centroid", lat),
                "lon": ("centroid", lon),
            },
            coords={"time": [0, 1], "centroid": centroid},
            attrs={"wind_averaging_period": "unspecified in frozen Lin source"},
        )
        event = {
            "event_id": "fixture",
            "event_position": 17,
            "event_weight_climate_fixed_effect_ht_analysis_yr": 0.25,
            "outer_radius_m": 800_000.0,
            "fixed_r0_catalogue_sha256": "catalogue-sha",
            "fixed_r0_distribution_contract_sha256": "distribution-sha",
        }
        compact = build_compact_hazard_footprint(
            rainfall, wind, event, right_censored=True
        )

        self.assertNotIn("time", compact.dims)
        self.assertNotIn("rainfall_rate", compact.variables)
        self.assertNotIn("near_surface_wind_speed", compact.variables)
        self.assertEqual(
            set(compact.data_vars),
            {
                "lat",
                "lon",
                "event_total_rainfall",
                "maximum_24h_rainfall",
                "event_maximum_near_surface_wind_speed",
            },
        )
        np.testing.assert_allclose(compact["lon"], lon)
        np.testing.assert_allclose(compact["event_total_rainfall"], [2, 4, 6])
        np.testing.assert_allclose(
            compact["event_maximum_near_surface_wind_speed"], [10, 20, 30]
        )
        self.assertEqual(compact.attrs["right_censored_at_15_day_limit"], 1)
        self.assertIn(
            "lower bound", compact.attrs["cumulative_rainfall_interpretation"]
        )
        self.assertEqual(compact.attrs["hazard_thresholds_applied"], 0)
        self.assertEqual(compact.attrs["damage_or_loss_model_applied"], 0)
        self.assertEqual(compact.attrs["outer_radius_m"], 800_000.0)
        self.assertEqual(
            compact.attrs["fixed_r0_distribution_contract_sha256"],
            "distribution-sha",
        )


if __name__ == "__main__":
    unittest.main()
