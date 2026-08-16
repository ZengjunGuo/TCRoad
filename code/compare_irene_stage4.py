#!/usr/bin/env python3
"""Build the minimal Stage IV comparison data for the Irene C15--TCR case.

The comparison window is fixed by hourly *ending time*: 2011-08-21 01 UTC
through 2011-08-30 00 UTC (216 one-hour accumulations).  Each TCR 0.05-degree
centroid in the approximate published visual map box is assigned one fixed
nearest Stage IV native-grid cell in spherical space.  Hourly Stage IV missing
values are then propagated through that fixed mapping.  No score threshold or
pass/fail decision is made here.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
from typing import Any

import numpy as np
from netCDF4 import Dataset, num2date
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree
import shapely
import shapefile
from shapely.geometry import shape
from shapely.ops import unary_union


START_END_TIME = datetime(2011, 8, 21, 1, tzinfo=timezone.utc)
FINAL_END_TIME = datetime(2011, 8, 30, 0, tzinfo=timezone.utc)
EXPECTED_HOURS = 216
COMMON_SUPPORT_RADIUS_KM = 300.0
EARTH_RADIUS_KM = 6371.0088

FIG1_LAT_MIN = 33.0
FIG1_LAT_MAX = 39.0
FIG1_LON_MIN = -81.0
FIG1_LON_MAX = -73.0

STAGE4_NX = 1121
STAGE4_NY = 881
STAGE4_DX_M = 4763.0
STAGE4_DY_M = 4763.0
STAGE4_EARTH_RADIUS_M = 6_367_470.0
STAGE4_LAT_TS_DEG = 60.0
STAGE4_LON_0_DEG = -105.0
STAGE4_FIRST_LAT_DEG = 23.117
STAGE4_FIRST_LON_DEG = -119.023
STAGE4_MISSING_VALUE = 9999.0

SHAPEFILE_PARTS = (".shp", ".shx", ".dbf", ".prj", ".cpg")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def expected_end_times() -> list[datetime]:
    values = [START_END_TIME + timedelta(hours=index) for index in range(EXPECTED_HOURS)]
    if values[-1] != FINAL_END_TIME:
        raise AssertionError("internal Stage IV comparison window error")
    return values


def extract_stage4_members(month_tar_path: Path) -> list[tuple[datetime, bytes]]:
    """Extract and decompress exactly the 216 required hourly GRIB messages."""

    requested = expected_end_times()
    by_day: dict[str, list[datetime]] = {}
    for end_time in requested:
        by_day.setdefault(end_time.strftime("%Y%m%d"), []).append(end_time)

    records: list[tuple[datetime, bytes]] = []
    seen: set[datetime] = set()
    with tarfile.open(month_tar_path, "r:") as outer:
        outer_names = outer.getnames()
        if len(outer_names) != 31 or set(outer_names) != {
            f"ST4.201108{day:02d}" for day in range(1, 32)
        }:
            raise ValueError("the frozen August 2011 Stage IV outer tar is not 31 daily members")
        for day, day_times in by_day.items():
            outer_name = f"ST4.{day}"
            extracted = outer.extractfile(outer_name)
            if extracted is None:
                raise FileNotFoundError(f"missing daily Stage IV member: {outer_name}")
            daily_bytes = extracted.read()
            with tarfile.open(fileobj=io.BytesIO(daily_bytes), mode="r:") as daily:
                inner_names = set(daily.getnames())
                for end_time in day_times:
                    inner_name = f"ST4.{end_time.strftime('%Y%m%d%H')}.01h.Z"
                    if inner_name not in inner_names:
                        raise FileNotFoundError(f"missing hourly Stage IV member: {inner_name}")
                    compressed = daily.extractfile(inner_name)
                    if compressed is None:
                        raise FileNotFoundError(f"cannot read hourly Stage IV member: {inner_name}")
                    result = subprocess.run(
                        ["uncompress", "-c"],
                        input=compressed.read(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )
                    if not result.stdout.startswith(b"GRIB"):
                        raise ValueError(f"decompressed member is not GRIB: {inner_name}")
                    records.append((end_time, result.stdout))
                    seen.add(end_time)
    if len(records) != EXPECTED_HOURS or seen != set(requested):
        raise ValueError(f"required Stage IV window is incomplete: {len(records)} records")
    records.sort(key=lambda item: item[0])
    if [item[0] for item in records] != requested:
        raise ValueError("Stage IV ending-time sequence is not strictly hourly")
    return records


def stage4_lon_lat() -> tuple[np.ndarray, np.ndarray, CRS]:
    crs = CRS.from_proj4(
        "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-105 "
        "+R=6367470 +units=m +no_defs"
    )
    forward = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    inverse = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    x0, y0 = forward.transform(STAGE4_FIRST_LON_DEG, STAGE4_FIRST_LAT_DEG)
    x = x0 + np.arange(STAGE4_NX, dtype=float) * STAGE4_DX_M
    y = y0 + np.arange(STAGE4_NY, dtype=float) * STAGE4_DY_M
    x_grid, y_grid = np.meshgrid(x, y)
    lon, lat = inverse.transform(x_grid, y_grid)
    return np.asarray(lon), np.asarray(lat), crs


def read_grib_values(grib_bytes: bytes, temp_path: Path) -> tuple[np.ndarray, dict[str, int]]:
    """Read one bitmap-aware Stage IV GRIB through the standard ecCodes binding."""

    from eccodes import (
        codes_get,
        codes_get_values,
        codes_grib_new_from_file,
        codes_release,
    )

    temp_path.write_bytes(grib_bytes)
    with temp_path.open("rb") as handle:
        gid = codes_grib_new_from_file(handle)
        if gid is None:
            raise ValueError("ecCodes found no GRIB message")
        try:
            keys = {
                key: int(codes_get(gid, key))
                for key in (
                    "edition",
                    "Ni",
                    "Nj",
                    "dataDate",
                    "dataTime",
                    "validityDate",
                    "validityTime",
                    "startStep",
                    "endStep",
                    "indicatorOfUnitOfTimeRange",
                    "bitmapPresent",
                )
            }
            if codes_get(gid, "shortName") != "tp":
                raise ValueError("Stage IV message is not Total Precipitation")
            if codes_get(gid, "gridType") != "polar_stereographic":
                raise ValueError("Stage IV grid is not polar_stereographic")
            if codes_get(gid, "stepType") != "accum":
                raise ValueError("Stage IV field is not an accumulation")
            if codes_get(gid, "units") != "kg m**-2":
                raise ValueError("unexpected Stage IV precipitation units")
            values = np.asarray(codes_get_values(gid), dtype=np.float32)
        finally:
            codes_release(gid)
    if keys["edition"] != 1 or keys["Ni"] != STAGE4_NX or keys["Nj"] != STAGE4_NY:
        raise ValueError(f"unexpected Stage IV GRIB geometry: {keys}")
    if keys["startStep"] != 0 or keys["endStep"] != 1:
        raise ValueError(f"Stage IV member is not a one-hour accumulation: {keys}")
    if keys["indicatorOfUnitOfTimeRange"] != 1:
        raise ValueError(f"Stage IV GRIB time-range unit is not hours: {keys}")
    values = values.reshape(STAGE4_NY, STAGE4_NX)
    values[values == STAGE4_MISSING_VALUE] = np.nan
    if np.any(values[np.isfinite(values)] < 0):
        raise ValueError("negative Stage IV precipitation encountered")
    return values, keys


def great_circle_km(
    center_lon: float, center_lat: float, lon: np.ndarray, lat: np.ndarray
) -> np.ndarray:
    center_lat_rad = np.radians(center_lat)
    lat_rad = np.radians(lat)
    dlat = lat_rad - center_lat_rad
    dlon = np.radians((lon - center_lon + 180.0) % 360.0 - 180.0)
    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(center_lat_rad) * np.cos(lat_rad) * np.sin(dlon / 2.0) ** 2
    )
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def unit_sphere_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon_rad = np.radians(np.asarray(lon, dtype=float))
    lat_rad = np.radians(np.asarray(lat, dtype=float))
    cos_lat = np.cos(lat_rad)
    return np.column_stack(
        (cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad))
    )


def fixed_nearest_stage4_mapping(
    stage4_lon: np.ndarray,
    stage4_lat: np.ndarray,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map Stage IV native cells to TCR centroids once in spherical space."""

    source_xyz = unit_sphere_xyz(stage4_lon.ravel(), stage4_lat.ravel())
    tree = cKDTree(source_xyz)
    chord, flat_index = tree.query(unit_sphere_xyz(target_lon, target_lat), k=1)
    flat_index = np.asarray(flat_index, dtype=np.int64)
    row, col = np.unravel_index(flat_index, stage4_lon.shape)
    angular = 2.0 * np.arcsin(np.clip(np.asarray(chord, dtype=float) / 2.0, 0.0, 1.0))
    distance_km = EARTH_RADIUS_KM * angular
    return flat_index, np.asarray(row), np.asarray(col), distance_km


def accumulate_mapped_stage4(
    records: list[tuple[datetime, bytes]],
    nearest_flat_index: np.ndarray,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
    track_lon: np.ndarray,
    track_lat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Accumulate raw context and paired 300-km support on TCR centroids."""

    count_target = target_lon.size
    raw_total = np.zeros(count_target, dtype=np.float64)
    raw_valid_count = np.zeros(count_target, dtype=np.uint16)
    common_total = np.zeros(count_target, dtype=np.float64)
    common_valid_count = np.zeros(count_target, dtype=np.uint16)
    common_expected_count = np.zeros(count_target, dtype=np.uint16)
    if track_lon.shape != (EXPECTED_HOURS,) or track_lat.shape != (EXPECTED_HOURS,):
        raise ValueError("common-support track slice must contain exactly 216 centers")
    with tempfile.TemporaryDirectory(prefix="irene_stage4_grib_") as temp_dir:
        temp_path = Path(temp_dir) / "hour.grb"
        for index, (expected, grib_bytes) in enumerate(records):
            values, keys = read_grib_values(grib_bytes, temp_path)
            actual = datetime.strptime(
                f"{keys['validityDate']:08d}{keys['validityTime']:04d}", "%Y%m%d%H%M"
            ).replace(tzinfo=timezone.utc)
            if actual != expected:
                raise ValueError(f"Stage IV valid time mismatch: {actual} != {expected}")
            mapped = values.ravel()[nearest_flat_index]
            valid = np.isfinite(mapped)
            raw_total[valid] += mapped[valid]
            raw_valid_count[valid] += 1

            within = great_circle_km(
                float(track_lon[index]), float(track_lat[index]), target_lon, target_lat
            ) <= COMMON_SUPPORT_RADIUS_KM
            common_expected_count[within] += 1
            common_valid = within & valid
            common_total[common_valid] += mapped[common_valid]
            common_valid_count[common_valid] += 1

    # A partial 216-hour sum is not a 216-hour total.  Preserve the count and
    # explicitly mask the context total unless all 216 hours are present.
    raw_total[raw_valid_count != EXPECTED_HOURS] = np.nan
    common_total[
        (common_expected_count == 0) | (common_valid_count != common_expected_count)
    ] = np.nan
    return (
        raw_total.astype(np.float32),
        raw_valid_count,
        common_total.astype(np.float32),
        common_valid_count,
        common_expected_count,
    )


def read_land_geometry(shapefile_path: Path) -> tuple[Any, list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    for suffix in SHAPEFILE_PARTS:
        part = shapefile_path.with_suffix(suffix)
        if not part.is_file():
            raise FileNotFoundError(f"missing Natural Earth shapefile component: {part}")
        records.append(artifact_record(part))
    projection_text = shapefile_path.with_suffix(".prj").read_text()
    if "WGS_1984" not in projection_text:
        raise ValueError("Natural Earth land shapefile is not WGS 84 longitude/latitude")
    with shapefile.Reader(str(shapefile_path)) as source:
        geometries = [shape(item.__geo_interface__) for item in source.shapes()]
    return unary_union(geometries), records


def read_tcr_event(
    tcr_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    with Dataset(tcr_path) as source:
        required = {"rainfall_rate", "lat", "lon"}
        if not required.issubset(source.variables):
            raise ValueError(f"TCR file missing variables: {required - set(source.variables)}")
        time_variable = source.variables["time"]
        times = time_variable[:]
        if times.size != EXPECTED_HOURS + 1:
            raise ValueError(f"TCR output must have 217 nodes, found {times.size}")
        decoded = num2date(
            times,
            units=time_variable.units,
            calendar=getattr(time_variable, "calendar", "standard"),
            only_use_cftime_datetimes=False,
        )
        decoded_utc = [
            datetime(
                value.year,
                value.month,
                value.day,
                value.hour,
                value.minute,
                value.second,
                value.microsecond,
                tzinfo=timezone.utc,
            )
            for value in decoded
        ]
        expected_nodes = [
            START_END_TIME - timedelta(hours=1) + timedelta(hours=index)
            for index in range(EXPECTED_HOURS + 1)
        ]
        if decoded_utc != expected_nodes:
            raise ValueError("TCR time coordinate is not 2011-08-21 00 through 2011-08-30 00 hourly")
        rainfall_variable = source.variables["rainfall_rate"]
        if getattr(rainfall_variable, "units", None) != "mm h-1":
            raise ValueError("TCR rainfall_rate units are not 'mm h-1'")
        first_rate = np.ma.asarray(rainfall_variable[0, :], dtype=float)
        first_rate = np.asarray(np.ma.filled(first_rate, np.nan), dtype=float)
        if not np.all(np.isfinite(first_rate)) or not np.all(first_rate == 0.0):
            raise ValueError("TCR first rainfall-rate node is not finite and identically zero")
        # The public CLIMADA angular-wind implementation zeros the first track
        # node.  The following 216 one-hour nodes correspond one-for-one to the
        # 216 Stage IV ending times frozen above.
        rates = np.ma.asarray(rainfall_variable[1:, :], dtype=float)
        rates = np.asarray(np.ma.filled(rates, np.nan), dtype=float)
        if rates.shape[0] != EXPECTED_HOURS:
            raise ValueError("TCR comparison slice is not exactly 216 hours")
        if not np.all(np.isfinite(rates)):
            raise ValueError("TCR 216-hour comparison slice contains non-finite values")
        lat = np.asarray(source.variables["lat"][:], dtype=float)
        lon = np.asarray(source.variables["lon"][:], dtype=float)
        date_strings = [value.isoformat() for value in decoded_utc]
    if np.any(rates[np.isfinite(rates)] < 0):
        raise ValueError("negative TCR rain rate encountered")
    return lon, lat, np.nansum(rates, axis=0), rates, date_strings


def read_track_centers(track_path: Path) -> tuple[np.ndarray, np.ndarray]:
    with Dataset(track_path) as source:
        time_variable = source.variables["time"]
        decoded = num2date(
            time_variable[:],
            units=time_variable.units,
            calendar=getattr(time_variable, "calendar", "standard"),
            only_use_cftime_datetimes=False,
        )
        decoded_utc = [
            datetime(
                value.year,
                value.month,
                value.day,
                value.hour,
                value.minute,
                value.second,
                value.microsecond,
                tzinfo=timezone.utc,
            )
            for value in decoded
        ]
        lon = np.asarray(source.variables["lon"][1:217], dtype=float)
        lat = np.asarray(source.variables["lat"][1:217], dtype=float)
    expected_nodes = [
        START_END_TIME - timedelta(hours=1) + timedelta(hours=index)
        for index in range(EXPECTED_HOURS + 1)
    ]
    if decoded_utc != expected_nodes:
        raise ValueError("track time coordinate is not 2011-08-21 00 through 2011-08-30 00 hourly")
    if lon.shape != (EXPECTED_HOURS,) or lat.shape != (EXPECTED_HOURS,):
        raise ValueError("hourly track must provide centers 1:217 for common support")
    return lon, lat


def comparison_statistics(observed: np.ndarray, modelled: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
    obs = np.asarray(observed[mask], dtype=float)
    mod = np.asarray(modelled[mask], dtype=float)
    if obs.size == 0 or not np.all(np.isfinite(obs)) or not np.all(np.isfinite(mod)):
        raise ValueError("comparison support is empty or non-finite")
    difference = mod - obs
    correlation = float(np.corrcoef(obs, mod)[0, 1]) if obs.size > 1 else None
    return {
        "comparison_cell_count": int(obs.size),
        "stage4_mean_mm": float(np.mean(obs)),
        "stage4_median_mm": float(np.median(obs)),
        "stage4_p95_mm": float(np.percentile(obs, 95)),
        "stage4_max_mm": float(np.max(obs)),
        "tcr_mean_mm": float(np.mean(mod)),
        "tcr_median_mm": float(np.median(mod)),
        "tcr_p95_mm": float(np.percentile(mod, 95)),
        "tcr_max_mm": float(np.max(mod)),
        "mean_bias_tcr_minus_stage4_mm": float(np.mean(difference)),
        "mean_absolute_error_mm": float(np.mean(np.abs(difference))),
        "root_mean_square_error_mm": float(np.sqrt(np.mean(difference**2))),
        "pearson_correlation": correlation,
    }


def write_output(
    path: Path,
    lon: np.ndarray,
    lat: np.ndarray,
    original_tcr_centroid: np.ndarray,
    stage4_row: np.ndarray,
    stage4_col: np.ndarray,
    stage4_distance_km: np.ndarray,
    stage4_context: np.ndarray,
    stage4_raw_valid_count: np.ndarray,
    stage4_common: np.ndarray,
    stage4_common_valid_count: np.ndarray,
    common_expected_count: np.ndarray,
    tcr_context: np.ndarray,
    tcr_common: np.ndarray,
    land: np.ndarray,
    support: np.ndarray,
) -> None:
    fill = np.float32(-9999.0)
    with Dataset(path, "w", format="NETCDF4") as target:
        target.createDimension("centroid", lon.size)
        target.setncatts(
            {
                "title": "Hurricane Irene Stage IV and C15-TCR comparison data in an approximate published visual bounding box",
                "reconstruction_scope": "full-lifecycle public reconstruction; not an exact reproduction of the published Figure 1",
                "time_window_definition": "216 one-hour accumulations by ending time",
                "time_window_start_end_time_utc": "2011-08-21T01:00:00Z",
                "time_window_final_end_time_utc": "2011-08-30T00:00:00Z",
                "approximate_visual_bbox": "33-39N, 81-73W",
                "mapping": "one fixed spherical-ECEF nearest Stage IV native-grid cell per original TCR 0.05-degree centroid; no interpolation or smoothing",
                "comparison_support": "approximate visual bbox AND Natural Earth 110m land AND Stage IV raw_valid_count=216; paired sums use only hours with great-circle distance to hourly track center <=300 km and require every such Stage IV hour to be finite",
                "decision_rule": "none; descriptive comparison only",
            }
        )
        variables = {
            "lon": ("f4", ("centroid",), lon, "degrees_east"),
            "lat": ("f4", ("centroid",), lat, "degrees_north"),
            "original_tcr_centroid": ("i4", ("centroid",), original_tcr_centroid, "1"),
            "nearest_stage4_row": ("i4", ("centroid",), stage4_row, "1"),
            "nearest_stage4_col": ("i4", ("centroid",), stage4_col, "1"),
            "nearest_stage4_distance": ("f4", ("centroid",), stage4_distance_km, "km"),
            "stage4_raw_216h_total_context_only": ("f4", ("centroid",), stage4_context, "mm"),
            "stage4_raw_valid_hour_count": ("u2", ("centroid",), stage4_raw_valid_count, "hours"),
            "stage4_common_support_total": ("f4", ("centroid",), stage4_common, "mm"),
            "stage4_common_support_valid_hour_count": ("u2", ("centroid",), stage4_common_valid_count, "hours"),
            "common_support_expected_hour_count": ("u2", ("centroid",), common_expected_count, "hours"),
            "tcr_raw_216h_total_context_only": ("f4", ("centroid",), tcr_context, "mm"),
            "tcr_common_support_total": ("f4", ("centroid",), tcr_common, "mm"),
            "natural_earth_land_mask": ("u1", ("centroid",), land.astype(np.uint8), "1"),
            "formal_comparison_support_mask": ("u1", ("centroid",), support.astype(np.uint8), "1"),
        }
        for name, (dtype, dims, values, units) in variables.items():
            kwargs: dict[str, Any] = {"zlib": True, "complevel": 4, "shuffle": True}
            if dtype.startswith("f"):
                kwargs["fill_value"] = fill
            variable = target.createVariable(name, dtype, dims, **kwargs)
            if dtype.startswith("f"):
                variable[:] = np.ma.masked_invalid(values)
            else:
                variable[:] = values
            variable.units = units


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage4-month-tar", type=Path, required=True)
    parser.add_argument("--tcr-rainfall", type=Path, required=True)
    parser.add_argument("--tcr-hourly-track", type=Path, required=True)
    parser.add_argument("--natural-earth-land-shp", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    for path in (
        args.stage4_month_tar,
        args.tcr_rainfall,
        args.tcr_hourly_track,
        args.natural_earth_land_shp,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError(f"atomic output target already exists: {args.output_dir}")
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
    try:
        hourly = extract_stage4_members(args.stage4_month_tar)
        stage4_lon, stage4_lat, stage4_crs = stage4_lon_lat()
        track_lon, track_lat = read_track_centers(args.tcr_hourly_track)
        all_tcr_lon, all_tcr_lat, all_tcr_total, all_tcr_rates, tcr_time_strings = read_tcr_event(
            args.tcr_rainfall
        )
        in_bbox = (
            np.isfinite(all_tcr_lon)
            & np.isfinite(all_tcr_lat)
            & np.isfinite(all_tcr_total)
            & (all_tcr_lat >= FIG1_LAT_MIN)
            & (all_tcr_lat <= FIG1_LAT_MAX)
            & (all_tcr_lon >= FIG1_LON_MIN)
            & (all_tcr_lon <= FIG1_LON_MAX)
        )
        original_tcr_centroid = np.flatnonzero(in_bbox).astype(np.int32)
        target_lon = all_tcr_lon[in_bbox]
        target_lat = all_tcr_lat[in_bbox]
        tcr_rates = all_tcr_rates[:, in_bbox]
        tcr_context = np.sum(tcr_rates, axis=0, dtype=np.float64).astype(np.float32)

        nearest_flat, nearest_row, nearest_col, nearest_distance = (
            fixed_nearest_stage4_mapping(
                stage4_lon, stage4_lat, target_lon, target_lat
            )
        )
        (
            stage4_context,
            stage4_raw_valid_count,
            stage4_common,
            stage4_common_valid_count,
            common_expected_count,
        ) = accumulate_mapped_stage4(
            hourly,
            nearest_flat,
            target_lon,
            target_lat,
            track_lon,
            track_lat,
        )

        tcr_common = np.zeros(target_lon.size, dtype=np.float64)
        tcr_common_valid_count = np.zeros(target_lon.size, dtype=np.uint16)
        tcr_common_expected_count = np.zeros(target_lon.size, dtype=np.uint16)
        for index in range(EXPECTED_HOURS):
            within = great_circle_km(
                float(track_lon[index]),
                float(track_lat[index]),
                target_lon,
                target_lat,
            ) <= COMMON_SUPPORT_RADIUS_KM
            valid = within & np.isfinite(tcr_rates[index])
            tcr_common_expected_count[within] += 1
            tcr_common_valid_count[valid] += 1
            tcr_common[valid] += tcr_rates[index, valid]
        if not np.array_equal(tcr_common_expected_count, common_expected_count):
            raise ValueError(
                "independently calculated TCR and Stage IV geometric common-support counts differ"
            )
        tcr_common[
            (tcr_common_expected_count == 0)
            | (tcr_common_valid_count != tcr_common_expected_count)
        ] = np.nan
        tcr_common = tcr_common.astype(np.float32)

        land_geometry, shapefile_records = read_land_geometry(args.natural_earth_land_shp)
        land = shapely.intersects_xy(land_geometry, target_lon, target_lat)
        support = (
            land
            & (stage4_raw_valid_count == EXPECTED_HOURS)
            & (common_expected_count > 0)
            & (stage4_common_valid_count == common_expected_count)
            & (tcr_common_valid_count == common_expected_count)
            & np.isfinite(stage4_common)
            & np.isfinite(tcr_common)
        )

        stats = comparison_statistics(stage4_common, tcr_common, support)
        output_nc = temporary / "irene_stage4_c15_tcr_approx_fig1_bbox_comparison.nc"
        write_output(
            output_nc,
            target_lon,
            target_lat,
            original_tcr_centroid,
            nearest_row,
            nearest_col,
            nearest_distance,
            stage4_context,
            stage4_raw_valid_count,
            stage4_common,
            stage4_common_valid_count,
            common_expected_count,
            tcr_context,
            tcr_common,
            land,
            support,
        )
        summary = {
            "status": "complete",
            "decision": "none_descriptive_comparison_only",
            "window": {
                "definition": "hourly accumulation ending time",
                "first": START_END_TIME.isoformat().replace("+00:00", "Z"),
                "last": FINAL_END_TIME.isoformat().replace("+00:00", "Z"),
                "hour_count": EXPECTED_HOURS,
                "tcr_node_slice": "rainfall_rate[1:217]",
            },
            "domain": {
                "approximate_published_visual_bbox": {
                    "lat_min": FIG1_LAT_MIN,
                    "lat_max": FIG1_LAT_MAX,
                    "lon_min": FIG1_LON_MIN,
                    "lon_max": FIG1_LON_MAX,
                },
                "land_mask": "Natural Earth 1:110m land polygons, boundary-inclusive",
                "complete_stage4_coverage_required": EXPECTED_HOURS,
                "raw_216h_complete_coverage_filter_role": (
                    "conservative data-completeness filter for the fixed nearest Stage IV "
                    "cell; paired physical support is independently defined by the hourly "
                    "300-km mask"
                ),
                "formal_pairing": "hourly great-circle distance to track center <=300 km before accumulation",
                "stage4_to_tcr_mapping": "one fixed spherical-ECEF nearest native Stage IV grid cell per original TCR centroid; no interpolation",
                "stage4_grid_crs_wkt": stage4_crs.to_wkt(),
                "stage4_grid_shape": [STAGE4_NY, STAGE4_NX],
                "bbox_tcr_centroid_count": int(target_lon.size),
                "bbox_land_tcr_centroid_count": int(np.count_nonzero(land)),
                "bbox_land_complete_stage4_tcr_centroid_count": int(
                    np.count_nonzero(land & (stage4_raw_valid_count == EXPECTED_HOURS))
                ),
                "nearest_stage4_distance_km": {
                    "minimum": float(np.min(nearest_distance)),
                    "median": float(np.median(nearest_distance)),
                    "maximum": float(np.max(nearest_distance)),
                },
                "common_support_hour_count_range": [
                    int(np.min(common_expected_count)),
                    int(np.max(common_expected_count)),
                ],
            },
            "statistics": stats,
            "inputs": {
                "stage4_month_tar": artifact_record(args.stage4_month_tar),
                "tcr_rainfall": artifact_record(args.tcr_rainfall),
                "tcr_hourly_track": artifact_record(args.tcr_hourly_track),
                "natural_earth_shapefile_components": shapefile_records,
            },
            "tcr_source_time_coordinate_count": len(tcr_time_strings),
            "software": {
                "script_sha256": sha256(Path(__file__).resolve()),
                "numpy": np.__version__,
                "shapely": shapely.__version__,
                "projection_constants": {
                    "earth_radius_m": STAGE4_EARTH_RADIUS_M,
                    "latitude_true_scale_deg": STAGE4_LAT_TS_DEG,
                    "central_longitude_deg": STAGE4_LON_0_DEG,
                    "dx_m": STAGE4_DX_M,
                    "dy_m": STAGE4_DY_M,
                },
            },
        }
        summary_path = temporary / "irene_stage4_c15_tcr_approx_fig1_bbox_summary.json"
        summary["outputs"] = {
            "netcdf": {
                "relative_path": output_nc.name,
                "final_path": str((args.output_dir / output_nc.name).resolve()),
                "bytes": output_nc.stat().st_size,
                "sha256": sha256(output_nc),
            },
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, args.output_dir)
        print(json.dumps({"status": "complete", "output_dir": str(args.output_dir), **stats}, sort_keys=True))
    except Exception:
        import shutil

        shutil.rmtree(temporary, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
