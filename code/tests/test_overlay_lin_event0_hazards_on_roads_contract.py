"""Contract tests for the first-event wind/rain road exposure overlay."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from overlay_lin_event0_hazards_on_roads import (  # noqa: E402
    EXPECTED_EVENT_ID,
    RAIN_VARIABLE,
    WIND_VARIABLE,
    build_road_overlap,
    canonical_longitude,
    load_hazard_contracts,
    match_event_centroids_to_road_centres,
    summarize_by_road_class,
    weighted_empirical_quantile,
    weighted_pearson_correlation,
)


def _road_fixture(path: Path) -> None:
    lat = np.asarray([10.05, 10.15], dtype=np.float32)
    lon = np.asarray([-179.95, -179.85], dtype=np.float32)
    lengths = np.zeros((5, 2, 2), dtype=np.float32)
    lengths[0] = [[1.0, 2.0], [3.0, 4.0]]
    lengths[1] = [[2.0, 1.0], [0.0, 3.0]]
    lengths[2] = [[0.5, 0.5], [0.5, 0.5]]
    lengths[3] = [[0.0, 1.0], [1.0, 0.0]]
    lengths[4] = [[4.0, 3.0], [2.0, 1.0]]
    dataset = xr.Dataset(
        data_vars={
            "road_length_by_class": (("road_class", "lat", "lon"), lengths),
            "road_length": (("lat", "lon"), lengths.sum(axis=0)),
        },
        coords={
            "road_class": np.arange(5, dtype=np.int8),
            "lat": lat,
            "lon": lon,
        },
        attrs={"grid_resolution_degrees": 0.1},
    )
    dataset["road_class"].attrs[
        "flag_meanings"
    ] = "highways primary secondary tertiary local"
    dataset["road_length_by_class"].attrs["units"] = "km"
    dataset["road_length"].attrs["units"] = "km"
    dataset.to_netcdf(path)


def _hazard_fixtures(
    rain_path: Path, wind_path: Path, *, mismatch: bool = False
) -> None:
    # The middle point (10.10, 180.00) is a valid 0.05-degree event point but
    # is not a 0.1-degree road-cell centre and must not be silently aggregated.
    centroid = np.arange(5, dtype=np.int32)
    lat = np.asarray([10.05, 10.05, 10.10, 10.15, 10.15])
    lon = np.asarray([180.05, 180.15, 180.00, 180.05, 180.15])
    rainfall = xr.Dataset(
        data_vars={
            RAIN_VARIABLE: ("centroid", [10.0, 20.0, 999.0, 30.0, 40.0]),
            "lat": ("centroid", lat),
            "lon": ("centroid", lon),
        },
        coords={"centroid": centroid},
        attrs={"event_id": EXPECTED_EVENT_ID},
    )
    rainfall[RAIN_VARIABLE].attrs["units"] = "mm"
    rainfall.to_netcdf(rain_path)

    wind_lat = lat.copy()
    if mismatch:
        wind_lat[0] += 0.01
    wind = xr.Dataset(
        data_vars={
            WIND_VARIABLE: ("centroid", [1.0, 2.0, 99.0, 3.0, 4.0]),
            "lat": ("centroid", wind_lat),
            "lon": ("centroid", lon),
        },
        coords={"centroid": centroid},
        attrs={
            "event_id": EXPECTED_EVENT_ID,
            "ten_minute_sustained_wind_claim": 0,
            "wind_averaging_period_conversion_applied": 0,
            "wind_averaging_period": "unspecified in frozen Lin v_trks source",
        },
    )
    wind[WIND_VARIABLE].attrs["units"] = "m s-1"
    wind.to_netcdf(wind_path)


class LinEventRoadOverlapContractTest(unittest.TestCase):
    def test_periodic_longitude_and_coincident_centre_only_matching(self):
        np.testing.assert_allclose(
            canonical_longitude(np.asarray([180.05, 180.15])),
            [-179.95, -179.85],
            atol=1e-12,
        )
        event_index, lat_index, lon_index = match_event_centroids_to_road_centres(
            np.asarray([10.05, 10.05, 10.10, 10.15, 10.15]),
            np.asarray([180.05, 180.15, 180.00, 180.05, 180.15]),
            np.asarray([10.05, 10.15]),
            np.asarray([-179.95, -179.85]),
        )
        np.testing.assert_array_equal(event_index, [0, 1, 3, 4])
        np.testing.assert_array_equal(lat_index, [0, 0, 1, 1])
        np.testing.assert_array_equal(lon_index, [0, 1, 0, 1])

    def test_hazard_contract_rejects_coordinate_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rain_path = root / "rain.nc"
            wind_path = root / "wind.nc"
            _hazard_fixtures(rain_path, wind_path, mismatch=True)
            with self.assertRaisesRegex(ValueError, "latitude coordinates differ"):
                load_hazard_contracts(rain_path, wind_path)

    def test_overlap_retains_five_lengths_without_threshold_or_loss(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            roads_path = root / "roads.nc"
            rain_path = root / "rain.nc"
            wind_path = root / "wind.nc"
            _road_fixture(roads_path)
            _hazard_fixtures(rain_path, wind_path)
            hazard = load_hazard_contracts(rain_path, wind_path)
            overlap, metadata = build_road_overlap(roads_path, hazard)

        self.assertEqual(overlap.sizes["road_class"], 5)
        self.assertEqual(overlap.sizes["road_cell"], 4)
        self.assertEqual(metadata["fine_event_centroids_not_at_road_centres"], 1)
        np.testing.assert_allclose(overlap[WIND_VARIABLE], [1.0, 2.0, 3.0, 4.0])
        np.testing.assert_allclose(overlap[RAIN_VARIABLE], [10.0, 20.0, 30.0, 40.0])
        self.assertEqual(overlap.attrs["hazard_thresholds_applied"], 0)
        self.assertEqual(overlap.attrs["fragility_or_damage_function_applied"], 0)
        self.assertEqual(overlap.attrs["spatial_interpolation_applied"], 0)
        self.assertIn(
            "not the full C15 r0 footprint", overlap.attrs["evaluation_support"]
        )
        self.assertNotIn(99.0, overlap[WIND_VARIABLE].values)
        self.assertNotIn(999.0, overlap[RAIN_VARIABLE].values)

    def test_weighted_statistics_are_explicit_and_threshold_free(self):
        value = np.asarray([1.0, 2.0, 3.0])
        weight = np.asarray([1.0, 1.0, 8.0])
        self.assertEqual(weighted_empirical_quantile(value, weight, 0.5), 3.0)
        self.assertAlmostEqual(
            weighted_pearson_correlation(value, 10.0 * value, weight), 1.0
        )

        overlap = xr.Dataset(
            data_vars={
                "road_length_by_class": (
                    ("road_class", "road_cell"),
                    np.tile(weight, (5, 1)),
                ),
                WIND_VARIABLE: ("road_cell", value),
                RAIN_VARIABLE: ("road_cell", 10.0 * value),
            },
            coords={"road_class": np.arange(5), "road_cell": np.arange(3)},
            attrs={"event_id": EXPECTED_EVENT_ID},
        )
        summary = summarize_by_road_class(overlap)
        self.assertFalse(summary["definition"]["hazard_thresholds_applied"])
        self.assertFalse(summary["definition"]["damage_or_loss_model_applied"])
        self.assertEqual(len(summary["road_classes"]), 5)
        self.assertEqual(
            summary["road_classes"][0]["road_length_weighted_p50_max_wind_m_s"],
            3.0,
        )


if __name__ == "__main__":
    unittest.main()
