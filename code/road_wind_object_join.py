#!/usr/bin/env python3
"""Join historical compact C15 event-max wind onto valued OSM objects.

This is the object-level seam for the frozen wind-asset kernel.  It does not
price anything.  Dollars stay in ``road_wind_asset_impact``.

Wind source
-----------
Only ``lin_road_domain_300km_v1`` compact files, variable
``event_maximum_near_surface_wind_speed`` (model-native C15, 1995–2014).
The eight ``METHOD_DOMAIN_PENDING`` event positions have no compact file
and contribute no wind.  Future SSP windows are not read.

Sampling
--------
The way representative ``lon``/``lat`` (same point as replacement cost /
OSM extract mid-point) is snapped to the compact 0.05° cell after
periodic longitude wrapping to ``[-180, 180)``.  The cell's compact
event-max wind is used if that cell exists in the event footprint;
otherwise the way has no wind for that event (non-finite ``v_c15_ms``).
No interpolation, no 0.1° overlay.

Usage
-----
  python3 road_wind_object_join.py join-event ways.csv --compact 00000.nc \\
      --output joined.csv
  python3 road_wind_object_join.py score-event ways.csv --compact 00000.nc \\
      --output impact.csv [--trees crowther.tif]
  python3 road_wind_object_join.py score-historical \\
      --valued-dir DIR --extract-dir DIR --compact-dir DIR \\
      --trees crowther.tif --output-dir DIR
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence

import numpy as np
from netCDF4 import Dataset

from road_wind_asset_impact import (
    CONTRACT_NAME as KERNEL_CONTRACT,
    HISTORICAL_WINDOW_YEARS,
    TREE_DENSITY_CAP_KM2,
    apply_rows,
    attach_crowther_density,
    bridge_collapse_usd,
    cleanup_usd,
    escobedo_cleanup_usd_per_km,
    load_income_levels,
    tree_break_c15_ms,
    tree_fail_prob,
)


SCRIPT_VERSION = "1.0.0"
HAZARD_RUN = "lin_road_domain_300km_v1"
COMPACT_WIND_VARIABLE = "event_maximum_near_surface_wind_speed"
COMPACT_GRID_STEP_DEG = 0.05
OSM_SNAPSHOT = "planet-260803"

# Frozen leftovers: no compact footprint, no wind, no dollars, weights kept
# in the climate sample and not used here.
METHOD_DOMAIN_PENDING_EVENT_POSITIONS: frozenset[int] = frozenset(
    {11902, 11944, 12357, 50194, 62311, 68925, 72126, 86977}
)


def canonical_longitude(longitude_deg: np.ndarray) -> np.ndarray:
    """Map finite longitudes periodically to ``[-180, 180)``. NaN stays NaN."""

    longitude = np.asarray(longitude_deg, dtype=float)
    out = np.full(longitude.shape, np.nan, dtype=float)
    finite = np.isfinite(longitude)
    if np.any(finite):
        out[finite] = (longitude[finite] + 180.0) % 360.0 - 180.0
    return out


def compact_cell_keys(
    lat_deg: np.ndarray, lon_deg: np.ndarray, *, step: float = COMPACT_GRID_STEP_DEG
) -> tuple[np.ndarray, np.ndarray]:
    """Integer 0.05° cell indices after periodic longitude wrap."""

    lat = np.asarray(lat_deg, dtype=float)
    lon = canonical_longitude(lon_deg)
    ilat = np.full(lat.shape, np.iinfo(np.int32).min, dtype=np.int32)
    ilon = np.full(lat.shape, np.iinfo(np.int32).min, dtype=np.int32)
    finite = np.isfinite(lat) & np.isfinite(lon)
    if np.any(finite):
        ilat[finite] = np.rint(lat[finite] / step).astype(np.int32)
        ilon[finite] = np.rint(lon[finite] / step).astype(np.int32)
    return ilat, ilon


def pack_cell_key(ilat: int, ilon: int) -> int:
    return (int(ilat) << 32) ^ (int(ilon) & 0xFFFFFFFF)


def load_compact_wind_index(
    compact_path: Path,
) -> tuple[dict[int, float], dict[str, Any]]:
    """Map packed cell key → event-max C15 (m s-1) for one compact file."""

    compact_path = compact_path.resolve()
    if not compact_path.is_file():
        raise FileNotFoundError(compact_path)
    with Dataset(compact_path) as dataset:
        if COMPACT_WIND_VARIABLE not in dataset.variables:
            raise ValueError(
                f"{compact_path} lacks {COMPACT_WIND_VARIABLE}"
            )
        if "lat" not in dataset.variables or "lon" not in dataset.variables:
            raise ValueError(f"{compact_path} lacks lat/lon")
        lat = np.asarray(dataset.variables["lat"][:], dtype=float)
        lon = np.asarray(dataset.variables["lon"][:], dtype=float)
        wind = np.asarray(dataset.variables[COMPACT_WIND_VARIABLE][:], dtype=float)
        event_position = int(getattr(dataset, "event_position", -1))
        event_id = str(getattr(dataset, "event_id", ""))
        weight = float(
            getattr(
                dataset,
                "event_weight_climate_fixed_effect_ht_analysis_yr",
                math.nan,
            )
        )
    if lat.shape != wind.shape or lon.shape != wind.shape:
        raise ValueError("compact lat/lon/wind shapes differ")
    if event_position in METHOD_DOMAIN_PENDING_EVENT_POSITIONS:
        raise ValueError(
            f"event_position {event_position} is METHOD_DOMAIN_PENDING; "
            "it must not contribute wind"
        )
    ilat, ilon = compact_cell_keys(lat, lon)
    index: dict[int, float] = {}
    for i in range(wind.size):
        if not np.isfinite(wind[i]):
            continue
        if ilat[i] == np.iinfo(np.int32).min:
            continue
        index[pack_cell_key(int(ilat[i]), int(ilon[i]))] = float(wind[i])
    metadata = {
        "path": str(compact_path),
        "event_position": event_position,
        "event_id": event_id,
        "event_weight_climate_fixed_effect_ht_analysis_yr": weight,
        "centroid_count": int(wind.size),
        "indexed_cell_count": len(index),
        "hazard_run": HAZARD_RUN,
        "wind_variable": COMPACT_WIND_VARIABLE,
        "grid_step_deg": COMPACT_GRID_STEP_DEG,
        "script_version": SCRIPT_VERSION,
    }
    return index, metadata


def sample_compact_index(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    index: Mapping[int, float],
) -> np.ndarray:
    """Nearest compact cell wind, or NaN if the way is outside the footprint."""

    lon = np.asarray(lon_deg, dtype=float)
    lat = np.asarray(lat_deg, dtype=float)
    if lon.shape != lat.shape:
        raise ValueError("lon and lat must share shape")
    ilat, ilon = compact_cell_keys(lat, lon)
    out = np.full(lat.shape, np.nan, dtype=float)
    flat_lat = np.ravel(ilat)
    flat_lon = np.ravel(ilon)
    flat_out = np.ravel(out)
    sentinel = np.iinfo(np.int32).min
    for i in range(flat_lat.size):
        if flat_lat[i] == sentinel:
            continue
        value = index.get(pack_cell_key(int(flat_lat[i]), int(flat_lon[i])))
        if value is not None:
            flat_out[i] = value
    return out


def join_event_wind(
    rows: Iterable[dict[str, Any]],
    compact_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Write ``v_c15_ms`` from one compact footprint onto valued-way rows."""

    table = [dict(row) for row in rows]
    index, metadata = load_compact_wind_index(compact_path)
    lon = np.asarray([_finite_or_nan(row.get("lon")) for row in table], dtype=float)
    lat = np.asarray([_finite_or_nan(row.get("lat")) for row in table], dtype=float)
    wind = sample_compact_index(lon, lat, index)
    matched = 0
    for row, value in zip(table, wind):
        row["v_c15_ms"] = float(value)
        row["wind_event_position"] = metadata["event_position"]
        row["wind_event_id"] = metadata["event_id"]
        if math.isfinite(value):
            matched += 1
    metadata = dict(metadata)
    metadata["ways"] = len(table)
    metadata["ways_with_footprint_wind"] = matched
    metadata["ways_outside_footprint"] = len(table) - matched
    return table, metadata


def attach_extract_coordinates(
    valued_rows: Iterable[dict[str, Any]],
    extract_rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Fill lon/lat from the OSM extract when the valued table omitted them."""

    xy: dict[str, tuple[str, str]] = {}
    for row in extract_rows:
        way_id = str(row.get("way_id", ""))
        if way_id:
            xy[way_id] = (str(row.get("lon", "")), str(row.get("lat", "")))
    attached: list[dict[str, Any]] = []
    for row in valued_rows:
        record = dict(row)
        have_lon = math.isfinite(_finite_or_nan(record.get("lon")))
        have_lat = math.isfinite(_finite_or_nan(record.get("lat")))
        if not (have_lon and have_lat):
            pair = xy.get(str(record.get("way_id", "")))
            if pair is not None:
                record["lon"] = pair[0]
                record["lat"] = pair[1]
        attached.append(record)
    return attached


def _finite_or_nan(value: Any) -> float:
    if value is None or value == "":
        return float("nan")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("no rows to write")
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_compact_fixture(
    path: Path,
    *,
    lat: Sequence[float],
    lon: Sequence[float],
    wind: Sequence[float],
    event_position: int = 0,
    event_id: str = "fixture-event",
    weight: float = 0.001,
) -> Path:
    """Write a tiny compact-like NetCDF for tests. Not a production file."""

    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    lat_a = np.asarray(lat, dtype=np.float32)
    lon_a = np.asarray(lon, dtype=np.float32)
    wind_a = np.asarray(wind, dtype=np.float32)
    if lat_a.shape != lon_a.shape or lat_a.shape != wind_a.shape:
        raise ValueError("fixture lat/lon/wind must share shape")
    with Dataset(path, "w", format="NETCDF4") as dataset:
        dataset.createDimension("centroid", lat_a.size)
        for name, values, units in (
            ("lat", lat_a, "degrees_north"),
            ("lon", lon_a, "degrees_east"),
            (COMPACT_WIND_VARIABLE, wind_a, "m s-1"),
        ):
            variable = dataset.createVariable(name, "f4", ("centroid",))
            variable[:] = values
            variable.units = units
        centroid = dataset.createVariable("centroid", "i4", ("centroid",))
        centroid[:] = np.arange(lat_a.size, dtype=np.int32)
        dataset.event_id = event_id
        dataset.event_position = np.int64(event_position)
        dataset.event_weight_climate_fixed_effect_ht_analysis_yr = np.float64(
            weight
        )
    return path


def sample_crowther_many(
    geotiff: Path, lon_deg: np.ndarray, lat_deg: np.ndarray
) -> np.ndarray:
    """Nearest Crowther pixel at each lon/lat. NaN if outside or nodata."""

    from osgeo import gdal

    lon = np.asarray(lon_deg, dtype=float)
    lat = np.asarray(lat_deg, dtype=float)
    if lon.shape != lat.shape:
        raise ValueError("lon and lat must share shape")
    dataset = gdal.Open(str(geotiff), gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(geotiff)
    try:
        transform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        inverse = gdal.InvGeoTransform(transform)
        if inverse is None:
            raise ValueError("Crowther GeoTIFF geotransform is not invertible")
        width = dataset.RasterXSize
        height = dataset.RasterYSize
        out = np.full(lon.shape, np.nan, dtype=float)
        finite = np.isfinite(lon) & np.isfinite(lat)
        if not np.any(finite):
            return out
        px = inverse[0] + inverse[1] * lon + inverse[2] * lat
        py = inverse[3] + inverse[4] * lon + inverse[5] * lat
        cols = np.floor(px).astype(np.int64)
        rows = np.floor(py).astype(np.int64)
        inside = finite & (cols >= 0) & (rows >= 0) & (cols < width) & (rows < height)
        if not np.any(inside):
            return out
        # Read the whole band once when many points hit it; else point reads.
        n_inside = int(np.count_nonzero(inside))
        if n_inside >= 64:
            array = np.asarray(band.ReadAsArray(), dtype=np.float32)
            values = np.asarray(array[rows[inside], cols[inside]], dtype=float)
        else:
            values = []
            for col, row in zip(cols[inside], rows[inside]):
                values.append(float(band.ReadAsArray(int(col), int(row), 1, 1)[0, 0]))
            values = np.asarray(values, dtype=float)
        if nodata is not None:
            values = np.where(values == float(nodata), np.nan, values)
        out[inside] = values
        return out
    finally:
        dataset = None


def score_event_rows(
    rows: list[dict[str, Any]],
    compact_path: Path,
    repo: Path,
    *,
    trees: Path | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    """Join one compact event, optionally Crowther, then the frozen kernel."""

    if trees is not None:
        attached = attach_crowther_density(rows, trees)
    else:
        attached = [dict(row) for row in rows]
    joined, join_meta = join_event_wind(attached, compact_path)
    impacted, totals = apply_rows(joined, repo)
    totals = dict(totals)
    totals["join"] = join_meta
    totals["kernel_contract"] = KERNEL_CONTRACT
    totals["join_script_version"] = SCRIPT_VERSION
    return impacted, totals, join_meta


def _iter_paired_shards(
    valued_dir: Path, extract_dir: Path
) -> Iterator[tuple[Path, Path]]:
    valued_dir = valued_dir.resolve()
    extract_dir = extract_dir.resolve()
    valued = sorted(valued_dir.glob("ways-*.valued.csv"))
    if not valued:
        raise FileNotFoundError(f"no ways-*.valued.csv in {valued_dir}")
    for valued_path in valued:
        stem = valued_path.name[: -len(".valued.csv")]
        extract_path = extract_dir / f"{stem}.csv"
        if not extract_path.is_file():
            raise FileNotFoundError(extract_path)
        yield valued_path, extract_path


def prepare_object_arrays(
    valued_dir: Path,
    extract_dir: Path,
    trees: Path,
) -> dict[str, np.ndarray]:
    """Load valued shards + extract lon/lat + Crowther into object arrays."""

    way_ids: list[int] = []
    iso3: list[str] = []
    highway: list[str] = []
    lons: list[float] = []
    lats: list[float] = []
    length: list[float] = []
    replacement: list[float] = []
    accepted: list[int] = []
    bridge: list[int] = []
    tunnel: list[int] = []
    pairs = list(_iter_paired_shards(valued_dir, extract_dir))
    for shard_i, (valued_path, extract_path) in enumerate(pairs, 1):
        print(
            f"[prepare-objects] shard {shard_i}/{len(pairs)} {valued_path.name}",
            flush=True,
        )
        xy: dict[str, tuple[float, float]] = {}
        with extract_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                xy[str(row["way_id"])] = (
                    _finite_or_nan(row.get("lon")),
                    _finite_or_nan(row.get("lat")),
                )
        with valued_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                way_id = str(row.get("way_id", ""))
                lon = _finite_or_nan(row.get("lon"))
                lat = _finite_or_nan(row.get("lat"))
                if not (math.isfinite(lon) and math.isfinite(lat)):
                    pair = xy.get(way_id)
                    if pair is not None:
                        lon, lat = pair
                way_ids.append(int(float(way_id)) if way_id else -1)
                iso3.append(str(row.get("iso3", "")).upper())
                highway.append(str(row.get("highway", "")))
                lons.append(lon)
                lats.append(lat)
                length.append(_finite_or_nan(row.get("length_km")) or 0.0)
                replacement.append(_finite_or_nan(row.get("replacement_usd")) or 0.0)
                accepted.append(1 if str(row.get("accepted", "1")).strip() not in {"0", "false", "False"} else 0)
                bridge.append(1 if str(row.get("is_bridge", "0")).strip() in {"1", "true", "True"} else 0)
                tunnel.append(1 if str(row.get("is_tunnel", "0")).strip() in {"1", "true", "True"} else 0)
    lon_a = np.asarray(lons, dtype=np.float64)
    lat_a = np.asarray(lats, dtype=np.float64)
    tree = sample_crowther_many(trees, lon_a, lat_a)
    return {
        "way_id": np.asarray(way_ids, dtype=np.int64),
        "iso3": np.asarray(iso3),
        "highway": np.asarray(highway),
        "lon": lon_a,
        "lat": lat_a,
        "length_km": np.asarray(length, dtype=np.float64),
        "replacement_usd": np.asarray(replacement, dtype=np.float64),
        "accepted": np.asarray(accepted, dtype=np.uint8),
        "is_bridge": np.asarray(bridge, dtype=np.uint8),
        "is_tunnel": np.asarray(tunnel, dtype=np.uint8),
        "tree_dens_km2": tree,
    }


def build_cell_to_ways(lon: np.ndarray, lat: np.ndarray) -> dict[int, np.ndarray]:
    ilat, ilon = compact_cell_keys(lat, lon)
    buckets: dict[int, list[int]] = {}
    sentinel = np.iinfo(np.int32).min
    for i in range(ilat.size):
        if ilat[i] == sentinel:
            continue
        key = pack_cell_key(int(ilat[i]), int(ilon[i]))
        buckets.setdefault(key, []).append(i)
    return {key: np.asarray(idx, dtype=np.int32) for key, idx in buckets.items()}


def list_compact_event_positions(compact_dir: Path) -> list[int]:
    positions: list[int] = []
    for path in compact_dir.glob("[0-9][0-9][0-9][0-9][0-9].nc"):
        positions.append(int(path.stem))
    positions.sort()
    return positions


def _tree_fail_prob_array(density: np.ndarray) -> np.ndarray:
    """Same rule as ``tree_fail_prob``, on an array."""

    values = np.asarray(density, dtype=np.float64)
    probability = np.zeros(values.shape, dtype=np.float64)
    valid = np.isfinite(values) & (values > 0.0)
    probability[valid] = np.minimum(values[valid], TREE_DENSITY_CAP_KM2) / TREE_DENSITY_CAP_KM2
    return probability


def score_historical(
    *,
    valued_dir: Path,
    extract_dir: Path,
    compact_dir: Path,
    trees: Path,
    output_dir: Path,
    repo: Path,
    only_event_position: Sequence[int] | None = None,
    max_events: int | None = None,
) -> dict[str, Any]:
    """Score the historical compact set onto valued objects with the frozen kernel."""

    output_dir.mkdir(parents=True, exist_ok=True)
    objects = prepare_object_arrays(valued_dir, extract_dir, trees)
    n_ways = int(objects["way_id"].size)
    cell_to_ways = build_cell_to_ways(objects["lon"], objects["lat"])
    incomes = load_income_levels(repo)
    income = np.array(
        [incomes.get(str(iso), "UMC") for iso in objects["iso3"]], dtype=object
    )
    tree_prob = _tree_fail_prob_array(objects["tree_dens_km2"])
    cleanup_eligible = (
        (objects["accepted"] == 1)
        & (objects["is_tunnel"] == 0)
        & (objects["is_bridge"] == 0)
    )
    cut = tree_break_c15_ms()
    cleanup_per_km = escobedo_cleanup_usd_per_km("central")

    positions = list_compact_event_positions(compact_dir)
    if only_event_position:
        wanted = {int(pos) for pos in only_event_position}
        positions = [pos for pos in positions if pos in wanted]
    skipped_pending = sorted(METHOD_DOMAIN_PENDING_EVENT_POSITIONS)
    positions = [
        pos for pos in positions if pos not in METHOD_DOMAIN_PENDING_EVENT_POSITIONS
    ]
    if max_events is not None:
        positions = positions[: int(max_events)]

    cleanup_weighted = np.zeros(n_ways, dtype=np.float64)
    cleanup_sum = np.zeros(n_ways, dtype=np.float64)
    bridge_weighted = np.zeros(n_ways, dtype=np.float64)
    bridge_sum = np.zeros(n_ways, dtype=np.float64)
    max_wind = np.full(n_ways, np.nan, dtype=np.float64)
    events_with_wind = np.zeros(n_ways, dtype=np.int32)
    bridge_peaks: dict[int, list[float]] = {}
    bridge_weights: dict[int, list[float]] = {}

    event_rows: list[dict[str, Any]] = []
    for rank, position in enumerate(positions, 1):
        compact_path = compact_dir / f"{position:05d}.nc"
        index, meta = load_compact_wind_index(compact_path)
        weight = float(meta["event_weight_climate_fixed_effect_ht_analysis_yr"])
        if not math.isfinite(weight):
            weight = 0.0
        hit_blocks: list[np.ndarray] = []
        wind_blocks: list[np.ndarray] = []
        for key, wind in index.items():
            members = cell_to_ways.get(key)
            if members is None:
                continue
            hit_blocks.append(members)
            wind_blocks.append(np.full(members.shape, float(wind), dtype=np.float64))
        if hit_blocks:
            idx = np.concatenate(hit_blocks)
            v_hit = np.concatenate(wind_blocks)
        else:
            idx = np.empty(0, dtype=np.int32)
            v_hit = np.empty(0, dtype=np.float64)
        ways_hit = int(idx.size)
        if ways_hit:
            previous = max_wind[idx]
            max_wind[idx] = np.where(np.isfinite(previous), np.maximum(previous, v_hit), v_hit)
            events_with_wind[idx] += 1
            is_bridge = objects["is_bridge"][idx] == 1
            for i, wind in zip(idx[is_bridge], v_hit[is_bridge]):
                bridge_peaks.setdefault(int(i), []).append(float(wind))
                bridge_weights.setdefault(int(i), []).append(weight)
            pavement = cleanup_eligible[idx] & (v_hit >= cut)
            dollars = np.zeros(ways_hit, dtype=np.float64)
            dollars[pavement] = (
                objects["length_km"][idx][pavement]
                * cleanup_per_km
                * tree_prob[idx][pavement]
            )
            np.add.at(cleanup_sum, idx, dollars)
            np.add.at(cleanup_weighted, idx, dollars * weight)
            event_cleanup = float(dollars.sum())
        else:
            event_cleanup = 0.0
        event_rows.append(
            {
                "event_position": position,
                "event_id": meta["event_id"],
                "event_weight_climate_fixed_effect_ht_analysis_yr": weight,
                "centroid_count": meta["centroid_count"],
                "ways_in_footprint": ways_hit,
                "cleanup_usd": event_cleanup,
                "bridge_usd_pending_second_pass": True,
            }
        )
        if rank % 200 == 0 or rank == len(positions):
            print(
                f"[score-historical] event {rank}/{len(positions)} "
                f"pos={position} ways_hit={ways_hit} cleanup={event_cleanup:.3f}",
                flush=True,
            )

    collapsed_bridges = 0
    for i, peaks in bridge_peaks.items():
        weights = bridge_weights[i]
        highway = str(objects["highway"][i])
        income_level = str(income[i])
        replacement = float(objects["replacement_usd"][i])
        accepted = bool(objects["accepted"][i])
        is_tunnel = bool(objects["is_tunnel"][i])
        for wind, weight in zip(peaks, weights):
            dollars = bridge_collapse_usd(
                replacement_usd=replacement,
                is_bridge=True,
                is_tunnel=is_tunnel,
                accepted=accepted,
                highway=highway,
                income_level=income_level,
                v_c15_ms=float(wind),
                historical_peaks_ms=peaks,
                window_years=HISTORICAL_WINDOW_YEARS,
            )
            if dollars > 0.0:
                collapsed_bridges += 1
            bridge_sum[i] += dollars
            bridge_weighted[i] += weight * dollars

    nonzero = (
        (cleanup_sum > 0.0)
        | (bridge_sum > 0.0)
        | np.isfinite(max_wind)
    )
    hit_idx = np.flatnonzero(nonzero)
    way_path = output_dir / "way_wind_asset.csv"
    with way_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "way_id",
                "iso3",
                "highway",
                "lon",
                "lat",
                "length_km",
                "replacement_usd",
                "accepted",
                "is_bridge",
                "is_tunnel",
                "tree_dens_km2",
                "max_v_c15_ms",
                "events_with_footprint_wind",
                "wind_cleanup_usd_sum",
                "wind_cleanup_usd_weighted",
                "wind_bridge_usd_sum",
                "wind_bridge_usd_weighted",
                "wind_asset_usd_sum",
            ],
        )
        writer.writeheader()
        for i in hit_idx:
            writer.writerow(
                {
                    "way_id": int(objects["way_id"][i]),
                    "iso3": str(objects["iso3"][i]),
                    "highway": str(objects["highway"][i]),
                    "lon": float(objects["lon"][i]),
                    "lat": float(objects["lat"][i]),
                    "length_km": float(objects["length_km"][i]),
                    "replacement_usd": float(objects["replacement_usd"][i]),
                    "accepted": int(objects["accepted"][i]),
                    "is_bridge": int(objects["is_bridge"][i]),
                    "is_tunnel": int(objects["is_tunnel"][i]),
                    "tree_dens_km2": float(objects["tree_dens_km2"][i]),
                    "max_v_c15_ms": float(max_wind[i]),
                    "events_with_footprint_wind": int(events_with_wind[i]),
                    "wind_cleanup_usd_sum": float(cleanup_sum[i]),
                    "wind_cleanup_usd_weighted": float(cleanup_weighted[i]),
                    "wind_bridge_usd_sum": float(bridge_sum[i]),
                    "wind_bridge_usd_weighted": float(bridge_weighted[i]),
                    "wind_asset_usd_sum": float(cleanup_sum[i] + bridge_sum[i]),
                }
            )

    event_path = output_dir / "event_wind_asset.csv"
    _write_csv(event_path, event_rows if event_rows else [
        {
            "event_position": -1,
            "note": "no compact events selected",
        }
    ])

    summary = {
        "schema_version": "1.0",
        "script_version": SCRIPT_VERSION,
        "kernel_contract": KERNEL_CONTRACT,
        "hazard_run": HAZARD_RUN,
        "wind_variable": COMPACT_WIND_VARIABLE,
        "grid_step_deg": COMPACT_GRID_STEP_DEG,
        "osm_snapshot": OSM_SNAPSHOT,
        "historical_window_years": HISTORICAL_WINDOW_YEARS,
        "compact_events_scored": len(positions),
        "method_domain_pending_event_positions": sorted(
            METHOD_DOMAIN_PENDING_EVENT_POSITIONS
        ),
        "method_domain_pending_contribute_wind": False,
        "method_domain_pending_contribute_dollars": False,
        "pending_ids_skipped": skipped_pending,
        "ways_in_valued_extract": n_ways,
        "ways_with_any_footprint_or_loss": int(hit_idx.size),
        "cleanup_usd_sum": float(cleanup_sum.sum()),
        "cleanup_usd_weighted": float(cleanup_weighted.sum()),
        "bridge_usd_sum": float(bridge_sum.sum()),
        "bridge_usd_weighted": float(bridge_weighted.sum()),
        "wind_asset_usd_sum": float(cleanup_sum.sum() + bridge_sum.sum()),
        "cleanup_ways": int(np.count_nonzero(cleanup_sum > 0.0)),
        "bridges_with_any_footprint": len(bridge_peaks),
        "bridge_event_collapses": collapsed_bridges,
        "crowther": str(trees.resolve()),
        "valued_dir": str(valued_dir.resolve()),
        "extract_dir": str(extract_dir.resolve()),
        "compact_dir": str(compact_dir.resolve()),
        "artifacts": {
            "way_wind_asset": str(way_path),
            "event_wind_asset": str(event_path),
        },
        "notes": {
            "cleanup": (
                "per-event L_km * 4979 * 1{V_C15 >= 25.3 m/s} * P(N); "
                "sum is over scored compact events; weighted uses the Lin "
                "climate-fixed-effect year weight on the compact file"
            ),
            "bridge": (
                "collapse only if gust > gmtra class V* and empirical RP "
                "from the 20-year compact peaks at that bridge exceeds "
                "design RP; pending-eight peaks are absent"
            ),
            "not_a_0p1deg_grid": True,
            "pavement_replacement_times_mdr": False,
        },
    }
    summary_path = output_dir / "historical_wind_asset.summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    join_p = sub.add_parser("join-event", help="write v_c15_ms from one compact file")
    join_p.add_argument("input_csv", type=Path)
    join_p.add_argument("--compact", type=Path, required=True)
    join_p.add_argument("--output", type=Path, required=True)

    score_p = sub.add_parser(
        "score-event",
        help="join one compact file then run the frozen dollar kernel",
    )
    score_p.add_argument("input_csv", type=Path)
    score_p.add_argument("--compact", type=Path, required=True)
    score_p.add_argument("--output", type=Path, required=True)
    score_p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    score_p.add_argument("--trees", type=Path, default=None)

    hist_p = sub.add_parser(
        "score-historical",
        help="score lin_road_domain_300km_v1 compact files onto valued objects",
    )
    hist_p.add_argument("--valued-dir", type=Path, required=True)
    hist_p.add_argument("--extract-dir", type=Path, required=True)
    hist_p.add_argument("--compact-dir", type=Path, required=True)
    hist_p.add_argument("--trees", type=Path, required=True)
    hist_p.add_argument("--output-dir", type=Path, required=True)
    hist_p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    hist_p.add_argument("--only-event-position", action="append", type=int, default=[])
    hist_p.add_argument("--max-events", type=int, default=None)

    args = parser.parse_args(argv)

    if args.command == "join-event":
        rows, meta = join_event_wind(_read_csv(args.input_csv), args.compact)
        _write_csv(args.output, rows)
        summary_path = Path(str(args.output) + ".summary.json")
        summary_path.write_text(
            json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(meta, indent=2, sort_keys=True))
        return 0

    if args.command == "score-event":
        impacted, totals, _join = score_event_rows(
            _read_csv(args.input_csv),
            args.compact,
            args.repo,
            trees=args.trees,
        )
        _write_csv(args.output, impacted)
        summary_path = Path(str(args.output) + ".summary.json")
        summary_path.write_text(
            json.dumps(totals, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(totals, indent=2, sort_keys=True))
        return 0

    score_historical(
        valued_dir=args.valued_dir,
        extract_dir=args.extract_dir,
        compact_dir=args.compact_dir,
        trees=args.trees,
        output_dir=args.output_dir,
        repo=args.repo,
        only_event_position=args.only_event_position or None,
        max_events=args.max_events,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
