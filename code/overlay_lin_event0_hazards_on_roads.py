#!/usr/bin/env python3
"""Overlay the first Lin-event wind and rain fields on the frozen road grid.

This is deliberately a descriptive exposure overlay, not an impact model.  It
samples the 0.05-degree event fields at coincident 0.1-degree road-cell centres,
retains all five road-length classes, and reports road-length-weighted
descriptive summaries.  It does not apply hazard thresholds, fragility or
damage functions, failure rules, or monetary valuation.

The two longitude conventions are reconciled periodically: event longitudes
may be in 0--360 degrees, while the road grid is in [-180, 180).  A road cell is
included only when its centre coincides with an event centroid after periodic
normalisation (within floating-point tolerance).  There is no spatial
interpolation or extrapolation at the moving-domain edge.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np
import xarray as xr


SCRIPT_VERSION = "1.0.0"
EXPECTED_EVENT_ID = "stream0000-year1995-track000002"
ROAD_GRID_STEP_DEG = 0.1
COORDINATE_TOLERANCE_DEG = 2.0e-5
ROAD_CLASS_IDS = np.arange(5, dtype=np.int8)
ROAD_CLASS_NAMES = ("highways", "primary", "secondary", "tertiary", "local")
RAIN_VARIABLE = "event_total_rainfall"
WIND_VARIABLE = "event_maximum_near_surface_wind_speed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_longitude(longitude_deg: np.ndarray) -> np.ndarray:
    """Map finite longitude values periodically to [-180, 180)."""

    longitude = np.asarray(longitude_deg, dtype=float)
    if not np.all(np.isfinite(longitude)):
        raise ValueError("longitude values must be finite")
    return (longitude + 180.0) % 360.0 - 180.0


def _periodic_absolute_difference(a_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    return np.abs(canonical_longitude(np.asarray(a_deg) - np.asarray(b_deg)))


def _require_1d_coordinate(
    dataset: xr.Dataset, name: str, dimension: str
) -> np.ndarray:
    if name not in dataset.variables or dataset[name].dims != (dimension,):
        raise ValueError(f"{name} must have dimensions ({dimension},)")
    values = np.asarray(dataset[name].values)
    if values.ndim != 1 or not np.all(np.isfinite(values.astype(float))):
        raise ValueError(f"{name} must be a finite one-dimensional coordinate")
    return values


def _validate_regular_road_axis(values: np.ndarray, name: str) -> None:
    numeric = np.asarray(values, dtype=float)
    if numeric.size < 2 or np.any(np.diff(numeric) <= 0.0):
        raise ValueError(f"road {name} coordinate must be strictly increasing")
    if not np.allclose(
        np.diff(numeric), ROAD_GRID_STEP_DEG, rtol=0.0, atol=COORDINATE_TOLERANCE_DEG
    ):
        raise ValueError(f"road {name} coordinate is not a regular 0.1-degree axis")


def load_hazard_contracts(
    rainfall_path: Path,
    windfield_path: Path,
    *,
    expected_event_id: str = EXPECTED_EVENT_ID,
) -> dict[str, Any]:
    """Read only the frozen event-total rain and event-maximum wind variables."""

    rainfall_path = rainfall_path.resolve()
    windfield_path = windfield_path.resolve()
    for path in (rainfall_path, windfield_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    with xr.open_dataset(rainfall_path) as rainfall:
        rain_event_id = str(rainfall.attrs.get("event_id", ""))
        if rain_event_id != expected_event_id:
            raise ValueError(f"unexpected rainfall event_id: {rain_event_id!r}")
        rain_lat = _require_1d_coordinate(rainfall, "lat", "centroid").astype(float)
        rain_lon = _require_1d_coordinate(rainfall, "lon", "centroid").astype(float)
        rain_centroid = _require_1d_coordinate(rainfall, "centroid", "centroid").copy()
        if RAIN_VARIABLE not in rainfall.variables:
            raise ValueError(f"rainfall file lacks {RAIN_VARIABLE}")
        rain = rainfall[RAIN_VARIABLE]
        if rain.dims != ("centroid",) or str(rain.attrs.get("units", "")) != "mm":
            raise ValueError(f"{RAIN_VARIABLE} must be centroid-dimensional in mm")
        rain_values = np.asarray(rain.values, dtype=float)

    with xr.open_dataset(windfield_path) as windfield:
        wind_event_id = str(windfield.attrs.get("event_id", ""))
        if wind_event_id != expected_event_id:
            raise ValueError(f"unexpected wind event_id: {wind_event_id!r}")
        if int(windfield.attrs.get("ten_minute_sustained_wind_claim", -1)) != 0:
            raise ValueError("wind field must not claim an unsupported 10-minute mean")
        if (
            int(windfield.attrs.get("wind_averaging_period_conversion_applied", -1))
            != 0
        ):
            raise ValueError(
                "wind averaging-period conversion is outside this contract"
            )
        averaging_period = str(windfield.attrs.get("wind_averaging_period", ""))
        if "unspecified" not in averaging_period.lower():
            raise ValueError("wind averaging period must remain explicitly unspecified")
        wind_lat = _require_1d_coordinate(windfield, "lat", "centroid").astype(float)
        wind_lon = _require_1d_coordinate(windfield, "lon", "centroid").astype(float)
        wind_centroid = _require_1d_coordinate(windfield, "centroid", "centroid").copy()
        if WIND_VARIABLE not in windfield.variables:
            raise ValueError(f"wind field lacks {WIND_VARIABLE}")
        wind = windfield[WIND_VARIABLE]
        if wind.dims != ("centroid",) or str(wind.attrs.get("units", "")) != "m s-1":
            raise ValueError(f"{WIND_VARIABLE} must be centroid-dimensional in m s-1")
        wind_values = np.asarray(wind.values, dtype=float)

    if not np.array_equal(rain_centroid, wind_centroid):
        raise ValueError("rainfall and wind centroid identifiers differ")
    if not np.allclose(rain_lat, wind_lat, rtol=0.0, atol=1e-10):
        raise ValueError("rainfall and wind latitude coordinates differ")
    if not np.all(_periodic_absolute_difference(rain_lon, wind_lon) <= 1.0e-10):
        raise ValueError("rainfall and wind longitude coordinates differ periodically")
    if not np.all(np.isfinite(rain_values)) or np.any(rain_values < 0.0):
        raise ValueError("event-total rainfall must be finite and non-negative")
    if not np.all(np.isfinite(wind_values)) or np.any(wind_values < 0.0):
        raise ValueError("event-maximum wind must be finite and non-negative")

    return {
        "event_id": expected_event_id,
        "centroid": rain_centroid,
        "lat": rain_lat,
        "lon": rain_lon,
        "event_total_rainfall_mm": rain_values,
        "event_maximum_wind_m_s": wind_values,
        "wind_averaging_period": averaging_period,
    }


def match_event_centroids_to_road_centres(
    event_lat: np.ndarray,
    event_lon: np.ndarray,
    road_lat: np.ndarray,
    road_lon: np.ndarray,
    *,
    tolerance_deg: float = COORDINATE_TOLERANCE_DEG,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return one-to-one coincident event/road indices, sorted by road cell.

    Fine-grid event points that lie between road-cell centres are intentionally
    not used.  This implements centre-point nearest-neighbour sampling without
    interpolation or domain-edge extrapolation.
    """

    event_latitude = np.asarray(event_lat, dtype=float)
    event_longitude = canonical_longitude(np.asarray(event_lon, dtype=float))
    road_latitude = np.asarray(road_lat, dtype=float)
    road_longitude = canonical_longitude(np.asarray(road_lon, dtype=float))
    if event_latitude.ndim != 1 or event_longitude.shape != event_latitude.shape:
        raise ValueError(
            "event latitude and longitude must share one-dimensional shape"
        )
    _validate_regular_road_axis(road_latitude, "latitude")
    _validate_regular_road_axis(np.asarray(road_lon, dtype=float), "longitude")

    lat_index = np.rint(
        (event_latitude - road_latitude[0]) / ROAD_GRID_STEP_DEG
    ).astype(np.int64)
    lon_index = np.rint(
        (event_longitude - road_longitude[0]) / ROAD_GRID_STEP_DEG
    ).astype(np.int64)
    valid_index = (
        (lat_index >= 0)
        & (lat_index < road_latitude.size)
        & (lon_index >= 0)
        & (lon_index < road_longitude.size)
    )
    event_index = np.flatnonzero(valid_index)
    if event_index.size == 0:
        raise ValueError("no event centroids fall within the road grid")
    candidate_lat_index = lat_index[event_index]
    candidate_lon_index = lon_index[event_index]
    coincident = (
        np.abs(event_latitude[event_index] - road_latitude[candidate_lat_index])
        <= tolerance_deg
    ) & (
        _periodic_absolute_difference(
            event_longitude[event_index], road_longitude[candidate_lon_index]
        )
        <= tolerance_deg
    )
    event_index = event_index[coincident]
    lat_index = candidate_lat_index[coincident]
    lon_index = candidate_lon_index[coincident]
    if event_index.size == 0:
        raise ValueError("no event centroids coincide with road-cell centres")

    flat_road_index = lat_index * road_longitude.size + lon_index
    if np.unique(flat_road_index).size != flat_road_index.size:
        raise ValueError("multiple event centroids map to the same road-cell centre")
    order = np.argsort(flat_road_index, kind="stable")
    return event_index[order], lat_index[order], lon_index[order]


def build_road_overlap(
    roads_path: Path,
    hazard: dict[str, Any],
) -> tuple[xr.Dataset, dict[str, Any]]:
    """Sample paired event hazards at frozen 0.1-degree road-cell centres."""

    roads_path = roads_path.resolve()
    if not roads_path.is_file():
        raise FileNotFoundError(roads_path)
    with xr.open_dataset(roads_path) as roads:
        road_lat = _require_1d_coordinate(roads, "lat", "lat").astype(float)
        road_lon = _require_1d_coordinate(roads, "lon", "lon").astype(float)
        _validate_regular_road_axis(road_lat, "latitude")
        _validate_regular_road_axis(road_lon, "longitude")
        if not math.isclose(
            float(roads.attrs.get("grid_resolution_degrees", np.nan)),
            ROAD_GRID_STEP_DEG,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("road grid-resolution metadata is not 0.1 degrees")
        road_classes = _require_1d_coordinate(roads, "road_class", "road_class").astype(
            np.int8
        )
        if not np.array_equal(road_classes, ROAD_CLASS_IDS):
            raise ValueError("road_class coordinate must be exactly 0..4")
        class_meanings = str(roads["road_class"].attrs.get("flag_meanings", ""))
        if tuple(class_meanings.split()) != ROAD_CLASS_NAMES:
            raise ValueError("road class meanings differ from the frozen five classes")
        if "road_length_by_class" not in roads.variables:
            raise ValueError("road file lacks road_length_by_class")
        by_class = roads["road_length_by_class"]
        if by_class.dims != ("road_class", "lat", "lon"):
            raise ValueError("road_length_by_class has unexpected dimensions")
        if str(by_class.attrs.get("units", "")) != "km":
            raise ValueError("road_length_by_class must be in km")

        event_index, lat_index, lon_index = match_event_centroids_to_road_centres(
            hazard["lat"], hazard["lon"], road_lat, road_lon
        )
        sampled = by_class.isel(
            lat=xr.DataArray(lat_index, dims="road_cell"),
            lon=xr.DataArray(lon_index, dims="road_cell"),
        ).transpose("road_class", "road_cell")
        road_length_by_class = np.asarray(sampled.values, dtype=float)
        if not np.all(np.isfinite(road_length_by_class)) or np.any(
            road_length_by_class < 0.0
        ):
            raise ValueError("sampled road lengths must be finite and non-negative")

        if "road_length" in roads.variables:
            road_total = roads["road_length"]
            if (
                road_total.dims != ("lat", "lon")
                or str(road_total.attrs.get("units", "")) != "km"
            ):
                raise ValueError("road_length has unexpected dimensions or units")
            sampled_total = np.asarray(
                road_total.isel(
                    lat=xr.DataArray(lat_index, dims="road_cell"),
                    lon=xr.DataArray(lon_index, dims="road_cell"),
                ).values,
                dtype=float,
            )
            if not np.allclose(
                sampled_total,
                road_length_by_class.sum(axis=0),
                rtol=2e-5,
                atol=2e-4,
            ):
                raise ValueError("road_length differs from the sum of five classes")

    event_maximum_wind = np.asarray(hazard["event_maximum_wind_m_s"], dtype=float)[
        event_index
    ]
    event_total_rainfall = np.asarray(hazard["event_total_rainfall_mm"], dtype=float)[
        event_index
    ]
    sampled_lat = road_lat[lat_index]
    sampled_lon = canonical_longitude(road_lon[lon_index])
    road_cell = np.arange(event_index.size, dtype=np.int32)
    overlap = xr.Dataset(
        data_vars={
            "road_length_by_class": (
                ("road_class", "road_cell"),
                road_length_by_class.astype(np.float32),
            ),
            "total_road_length": (
                "road_cell",
                road_length_by_class.sum(axis=0).astype(np.float32),
            ),
            WIND_VARIABLE: ("road_cell", event_maximum_wind.astype(np.float32)),
            RAIN_VARIABLE: ("road_cell", event_total_rainfall.astype(np.float32)),
            "lat": ("road_cell", sampled_lat.astype(np.float32)),
            "lon": ("road_cell", sampled_lon.astype(np.float32)),
            "global_road_latitude_index": ("road_cell", lat_index.astype(np.int32)),
            "global_road_longitude_index": ("road_cell", lon_index.astype(np.int32)),
            "source_event_centroid": (
                "road_cell",
                np.asarray(hazard["centroid"])[event_index].astype(np.int32),
            ),
        },
        coords={"road_cell": road_cell, "road_class": ROAD_CLASS_IDS},
        attrs={
            "event_id": hazard["event_id"],
            "product_type": "descriptive road exposure overlay; not loss or damage",
            "evaluation_support": (
                "common moving <=300-km wind-rain evaluation support sampled at "
                "coincident 0.1-degree road-cell centres; not the full C15 r0 footprint"
            ),
            "hazard_sampling_method": (
                "nearest 0.05-degree event centroid at each coincident 0.1-degree "
                "road-cell centre"
            ),
            "spatial_interpolation_applied": 0,
            "domain_edge_extrapolation_applied": 0,
            "periodic_longitude_normalisation": "[-180,180)",
            "hazard_thresholds_applied": 0,
            "fragility_or_damage_function_applied": 0,
            "wind_averaging_period": hazard["wind_averaging_period"],
        },
    )
    overlap["road_class"].attrs.update(
        long_name="motor road hierarchy class",
        flag_values=ROAD_CLASS_IDS,
        flag_meanings=" ".join(ROAD_CLASS_NAMES),
    )
    overlap["road_length_by_class"].attrs.update(
        units="km", long_name="mapped motor-road length by hierarchy class"
    )
    overlap["total_road_length"].attrs.update(
        units="km", long_name="total mapped motor-road length in road cell"
    )
    overlap[WIND_VARIABLE].attrs.update(
        units="m s-1",
        averaging_period=hazard["wind_averaging_period"],
        long_name="event maximum model-native near-surface wind speed",
    )
    overlap[RAIN_VARIABLE].attrs.update(
        units="mm", long_name="event-total rainfall over all native one-hour nodes"
    )
    overlap["lat"].attrs.update(units="degrees_north", standard_name="latitude")
    overlap["lon"].attrs.update(units="degrees_east", standard_name="longitude")
    metadata = {
        "road_cell_count": int(event_index.size),
        "source_event_centroid_count": int(np.asarray(hazard["centroid"]).size),
        "fine_event_centroids_used": int(event_index.size),
        "fine_event_centroids_not_at_road_centres": int(
            np.asarray(hazard["centroid"]).size - event_index.size
        ),
        "latitude_range_degrees_north": [
            float(sampled_lat.min()),
            float(sampled_lat.max()),
        ],
        "longitude_range_degrees_east": [
            float(sampled_lon.min()),
            float(sampled_lon.max()),
        ],
    }
    return overlap, metadata


def weighted_empirical_quantile(
    values: np.ndarray, weights: np.ndarray, probability: float
) -> float | None:
    """Inverse weighted empirical CDF: first value reaching q times total weight."""

    value = np.asarray(values, dtype=float)
    weight = np.asarray(weights, dtype=float)
    if value.shape != weight.shape or value.ndim != 1:
        raise ValueError("weighted quantile inputs must share one-dimensional shape")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("quantile probability must be between zero and one")
    valid = np.isfinite(value) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return None
    value = value[valid]
    weight = weight[valid]
    order = np.argsort(value, kind="stable")
    value = value[order]
    cumulative = np.cumsum(weight[order])
    target = probability * cumulative[-1]
    index = int(np.searchsorted(cumulative, target, side="left"))
    return float(value[min(index, value.size - 1)])


def weighted_pearson_correlation(
    first: np.ndarray, second: np.ndarray, weights: np.ndarray
) -> float | None:
    """Road-length-weighted Pearson correlation for paired finite hazards."""

    x = np.asarray(first, dtype=float)
    y = np.asarray(second, dtype=float)
    weight = np.asarray(weights, dtype=float)
    if x.shape != y.shape or x.shape != weight.shape or x.ndim != 1:
        raise ValueError("weighted correlation inputs must share one-dimensional shape")
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weight) & (weight > 0.0)
    if np.count_nonzero(valid) < 2:
        return None
    x = x[valid]
    y = y[valid]
    weight = weight[valid]
    total = float(weight.sum())
    mean_x = float(np.sum(weight * x) / total)
    mean_y = float(np.sum(weight * y) / total)
    centered_x = x - mean_x
    centered_y = y - mean_y
    variance_x = float(np.sum(weight * centered_x**2) / total)
    variance_y = float(np.sum(weight * centered_y**2) / total)
    if variance_x <= 0.0 or variance_y <= 0.0:
        return None
    covariance = float(np.sum(weight * centered_x * centered_y) / total)
    return covariance / math.sqrt(variance_x * variance_y)


def summarize_by_road_class(overlap: xr.Dataset) -> dict[str, Any]:
    """Return threshold-free, road-length-weighted paired-hazard summaries."""

    wind = np.asarray(overlap[WIND_VARIABLE].values, dtype=float)
    rainfall = np.asarray(overlap[RAIN_VARIABLE].values, dtype=float)
    lengths = np.asarray(overlap["road_length_by_class"].values, dtype=float)
    rows: list[dict[str, Any]] = []
    quantiles = (("p50", 0.50), ("p90", 0.90), ("p95", 0.95), ("p99", 0.99))
    for class_id, class_name in enumerate(ROAD_CLASS_NAMES):
        weight = lengths[class_id]
        positive = weight > 0.0
        total = float(weight.sum())
        row: dict[str, Any] = {
            "road_class": class_id,
            "road_class_name": class_name,
            "road_cell_count_with_positive_length": int(np.count_nonzero(positive)),
            "road_length_in_sampled_event_footprint_km": total,
        }
        if total > 0.0:
            row["road_length_weighted_mean_max_wind_m_s"] = float(
                np.sum(weight * wind) / total
            )
            row["road_length_weighted_mean_event_total_rainfall_mm"] = float(
                np.sum(weight * rainfall) / total
            )
            for label, probability in quantiles:
                row[f"road_length_weighted_{label}_max_wind_m_s"] = (
                    weighted_empirical_quantile(wind, weight, probability)
                )
                row[f"road_length_weighted_{label}_event_total_rainfall_mm"] = (
                    weighted_empirical_quantile(rainfall, weight, probability)
                )
            row["road_length_weighted_pearson_wind_rain"] = (
                weighted_pearson_correlation(wind, rainfall, weight)
            )
        else:
            for name in (
                "road_length_weighted_mean_max_wind_m_s",
                "road_length_weighted_mean_event_total_rainfall_mm",
                *(
                    item
                    for label, _ in quantiles
                    for item in (
                        f"road_length_weighted_{label}_max_wind_m_s",
                        f"road_length_weighted_{label}_event_total_rainfall_mm",
                    )
                ),
                "road_length_weighted_pearson_wind_rain",
            ):
                row[name] = None
        rows.append(row)
    return {
        "schema_version": "1.0",
        "event_id": str(overlap.attrs["event_id"]),
        "definition": {
            "support": (
                "mapped road length in 0.1-degree cells whose centres coincide "
                "with a 0.05-degree event centroid inside the common moving "
                "<=300-km wind-rain evaluation support; not the full C15 r0 footprint"
            ),
            "weight": "road_length_by_class in km",
            "quantile": "inverse weighted empirical CDF",
            "joint_statistic": "road-length-weighted Pearson wind-rain correlation",
            "hazard_thresholds_applied": False,
            "damage_or_loss_model_applied": False,
        },
        "road_classes": rows,
    }


def write_netcdf(dataset: xr.Dataset, path: Path) -> None:
    encoding: dict[str, dict[str, Any]] = {}
    for name, variable in dataset.data_vars.items():
        if variable.dtype.kind not in {"U", "S", "O"}:
            encoding[name] = {"zlib": True, "complevel": 4, "shuffle": True}
    dataset.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding=encoding)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def write_summary_csv(path: Path, summary: dict[str, Any]) -> None:
    rows = summary["road_classes"]
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roads-nc", type=Path, required=True)
    parser.add_argument("--rainfall-nc", type=Path, required=True)
    parser.add_argument("--windfield-nc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-event-id", default=EXPECTED_EVENT_ID)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"atomic output target exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        hazard = load_hazard_contracts(
            args.rainfall_nc,
            args.windfield_nc,
            expected_event_id=args.expected_event_id,
        )
        overlap, overlap_metadata = build_road_overlap(args.roads_nc, hazard)
        summary = summarize_by_road_class(overlap)

        overlap_path = temporary / "lin_event0_road_grid_joint_exposure.nc"
        summary_json_path = (
            temporary / "lin_event0_road_class_joint_exposure_summary.json"
        )
        summary_csv_path = (
            temporary / "lin_event0_road_class_joint_exposure_summary.csv"
        )
        write_netcdf(overlap, overlap_path)
        write_json(summary_json_path, summary)
        write_summary_csv(summary_csv_path, summary)

        runner_path = Path(__file__).resolve()
        input_paths = {
            "roads": args.roads_nc.resolve(),
            "rainfall": args.rainfall_nc.resolve(),
            "windfield": args.windfield_nc.resolve(),
        }
        artifact_paths = {
            "road_grid_joint_exposure": overlap_path,
            "road_class_summary_json": summary_json_path,
            "road_class_summary_csv": summary_csv_path,
        }
        manifest = {
            "schema_version": "1.0",
            "script_version": SCRIPT_VERSION,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "completed",
            "event_id": args.expected_event_id,
            "runner": {
                "path": str(runner_path),
                "bytes": runner_path.stat().st_size,
                "sha256": sha256(runner_path),
            },
            "inputs": {
                name: {
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for name, path in input_paths.items()
            },
            "overlay_contract": {
                "road_grid_resolution_degrees": ROAD_GRID_STEP_DEG,
                "hazard_sampling": (
                    "coincident-centre nearest sampling from 0.05-degree event field"
                ),
                "evaluation_support": (
                    "common moving <=300-km wind-rain support; not the full C15 r0 "
                    "footprint"
                ),
                "coordinate_tolerance_degrees": COORDINATE_TOLERANCE_DEG,
                "longitude_normalisation": "periodic [-180,180)",
                "spatial_interpolation": False,
                "domain_edge_extrapolation": False,
                "hazard_thresholds": False,
                "fragility_damage_or_loss_model": False,
                "wind_averaging_period": hazard["wind_averaging_period"],
            },
            "result_summary": overlap_metadata,
            "artifacts": {
                name: {
                    "relative_path": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for name, path in artifact_paths.items()
            },
        }
        manifest_path = temporary / "lin_event0_road_overlap.manifest.json"
        write_json(manifest_path, manifest)
        os.replace(temporary, output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": "completed",
                    "event_id": args.expected_event_id,
                    **overlap_metadata,
                    "output_dir": str(output_dir),
                },
                indent=2,
                sort_keys=True,
            )
        )
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
