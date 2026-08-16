#!/usr/bin/env python3
"""Screen the 100,000-track catalogue for storms that enter the 300 km road domain.

Gori et al. (2022) compute physics for every storm in the analysis set.  For this
global-road study the analysis set is every accepted Lin track that comes within
the frozen TCR support of mapped motor roads.  The coarse 1-degree occupancy
tiles are dilated conservatively so the screen can only over-include events.

The output sample is ordered by source_track_index.  Event-level r0 is later
bound to that same order from the 100,000-draw sequence.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any
import uuid

from netCDF4 import Dataset, chartostring
import numpy as np


SCRIPT_VERSION = "1.0.0"
SCHEMA_VERSION = "tcr-road-domain-sample-v1"
JOINT_SUPPORT_KM = 300.0
EARTH_RADIUS_KM = 6371.0088
KM_PER_DEG_LAT = math.pi * EARTH_RADIUS_KM / 180.0
CHUNK_EVENTS = 2_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _nc_strings(variable: Any) -> np.ndarray:
    values = variable[:]
    if values.dtype.kind in {"U", "O"}:
        return np.asarray(values, dtype=str)
    if values.dtype.kind == "S" and values.ndim == 2:
        return np.asarray(chartostring(values), dtype=str)
    return np.asarray(values).astype(str)


def dilated_road_mask(occupancy_path: Path, support_km: float) -> dict[str, Any]:
    occupancy = np.load(occupancy_path)
    tile_deg = float(occupancy["tile_deg"])
    nlon = int(occupancy["nlon"])
    nlat = int(occupancy["nlat"])
    tile_id = np.asarray(occupancy["tile_id"], dtype=np.int32)
    occupied = np.zeros((nlat, nlon), dtype=bool)
    occupied[tile_id // nlon, tile_id % nlon] = True
    # Extra half-tile keeps the 300 km disk inside the dilated cells at the equator.
    radius_tiles = int(math.ceil((support_km + 0.5 * tile_deg * KM_PER_DEG_LAT) / (tile_deg * KM_PER_DEG_LAT)))
    padded = np.pad(occupied, ((radius_tiles, radius_tiles), (0, 0)), mode="constant")
    padded = np.concatenate(
        [padded[:, -radius_tiles:], padded, padded[:, :radius_tiles]], axis=1
    )
    window = 2 * radius_tiles + 1
    dilated = np.lib.stride_tricks.sliding_window_view(padded, (window, window))
    mask = np.any(dilated, axis=(-2, -1))
    return {
        "mask": mask,
        "tile_deg": tile_deg,
        "nlon": nlon,
        "nlat": nlat,
        "occupied_tile_count": int(occupied.sum()),
        "dilated_tile_count": int(mask.sum()),
        "radius_tiles": radius_tiles,
    }


def screen_tracks(
    tracks_path: Path,
    mask: np.ndarray,
    *,
    tile_deg: float,
    nlon: int,
    nlat: int,
) -> np.ndarray:
    with Dataset(tracks_path) as tracks:
        n_trk = len(tracks.dimensions["n_trk"])
        n_time = len(tracks.dimensions["time"])
        lat_var = tracks.variables["lat_trks"]
        lon_var = tracks.variables["lon_trks"]
        hits = np.zeros(n_trk, dtype=bool)
        for start in range(0, n_trk, CHUNK_EVENTS):
            stop = min(n_trk, start + CHUNK_EVENTS)
            lat = np.asarray(lat_var[start:stop, :], dtype=np.float64)
            lon = np.asarray(lon_var[start:stop, :], dtype=np.float64)
            finite = np.isfinite(lat) & np.isfinite(lon)
            if not np.any(finite):
                continue
            lon_wrapped = np.where(finite, np.mod(lon + 180.0, 360.0), 0.0)
            lat_clipped = np.where(finite, np.clip(lat, -90.0, 90.0 - 1e-12), 0.0)
            ix = np.floor(lon_wrapped / tile_deg).astype(np.int32) % nlon
            iy = np.floor((lat_clipped + 90.0) / tile_deg).astype(np.int32)
            iy = np.clip(iy, 0, nlat - 1)
            inside = mask[iy, ix] & finite
            hits[start:stop] = np.any(inside, axis=1)
        return hits, n_trk, n_time


def write_sample(
    weighted_path: Path,
    selected_track_index: np.ndarray,
    output_nc: Path,
    manifest_path: Path,
    *,
    tracks_path: Path,
    occupancy_path: Path,
    screen_record: dict[str, Any],
) -> dict[str, Any]:
    selected = np.asarray(selected_track_index, dtype=np.int64)
    selected.sort()
    with Dataset(weighted_path) as source:
        source_index = np.asarray(source.variables["source_track_index"][:], dtype=np.int64)
        if source_index.size != 100_000 or not np.array_equal(
            source_index, np.arange(100_000, dtype=np.int64)
        ):
            raise ValueError("weighted catalogue is not source_track_index 0..99999")
        keep = selected
        event_ids = _nc_strings(source.variables["event_id"])[keep]
        output_nc.parent.mkdir(parents=True, exist_ok=True)
        temporary = output_nc.with_name(f".{output_nc.name}.{uuid.uuid4().hex}.tmp")
        with Dataset(temporary, "w", format="NETCDF4") as output:
            output.createDimension("event", keep.size)
            def copy_var(name: str, dtype=None) -> None:
                values = source.variables[name][keep]
                if name == "event_id" or dtype is str:
                    variable = output.createVariable(name, str, ("event",))
                    variable[:] = np.asarray(_nc_strings(source.variables[name])[keep], dtype=object)
                    return
                if np.asarray(values).dtype.kind in {"U", "O", "S"}:
                    variable = output.createVariable(name, str, ("event",))
                    variable[:] = np.asarray(_nc_strings(source.variables[name])[keep], dtype=object)
                    return
                store = values if dtype is None else np.asarray(values, dtype=dtype)
                variable = output.createVariable(name, store.dtype, ("event",))
                variable[:] = store

            copy_var("event_id")
            copy_var("source_track_index", np.int32)
            positions = output.createVariable(
                "source_catalogue_event_position", "i4", ("event",)
            )
            positions[:] = selected.astype(np.int32)
            copy_var("task_year", np.int32)
            copy_var("lin_seed_month", np.int32)
            copy_var("lin_seed_basin")
            copy_var("threshold_genesis_region")
            copy_var("threshold_genesis_datetime")
            copy_var("threshold_lysis_datetime")
            copy_var("threshold_genesis_native_index", np.int32)
            copy_var("threshold_lysis_native_index", np.int32)
            weight = output.createVariable(
                "event_weight_climate_fixed_effect_ht_analysis_yr", "f8", ("event",)
            )
            weight[:] = np.asarray(
                source.variables["event_weight_climate_fixed_effect_primary_yr"][keep],
                dtype=np.float64,
            )
            weight.long_name = (
                "climate-window annual occurrence weight; copied from "
                "event_weight_climate_fixed_effect_primary_yr"
            )
            output.title = "Road-domain Lin events within dilated 300 km motor-road support"
            output.schema_version = SCHEMA_VERSION
            output.script_version = SCRIPT_VERSION
            output.status = "FROZEN_ROAD_DOMAIN"
            output.joint_support_km = JOINT_SUPPORT_KM
            output.selection_rule = (
                "keep accepted tracks whose hourly centre occupies a 1-degree "
                "tile after conservative dilation of mapped motor-road occupancy"
            )
            output.source_weighted_catalogue_sha256 = sha256(weighted_path)
            output.source_tracks_sha256 = sha256(tracks_path)
            output.source_occupancy_sha256 = sha256(occupancy_path)

        os.replace(temporary, output_nc)

    weight_sum = float(
        np.asarray(
            Dataset(output_nc).variables[
                "event_weight_climate_fixed_effect_ht_analysis_yr"
            ][:]
        ).sum()
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "script_version": SCRIPT_VERSION,
        "status": "FROZEN_ROAD_DOMAIN",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_set": (
            "every accepted Lin track that enters the frozen 300 km TCR "
            "support of mapped motor roads; Gori-style full physics set"
        ),
        "not_a_random_subsample": True,
        "joint_support_km": JOINT_SUPPORT_KM,
        "screen": screen_record,
        "event_count": int(selected.size),
        "weight_sum_climate_fixed_effect_primary_yr": weight_sum,
        "artifacts": {
            "sample_netcdf": {
                "path": str(output_nc),
                "bytes": output_nc.stat().st_size,
                "sha256": sha256(output_nc),
            }
        },
        "sources": {
            "weighted_catalogue": {
                "path": str(weighted_path),
                "sha256": sha256(weighted_path),
            },
            "hourly_tracks": {"path": str(tracks_path), "sha256": sha256(tracks_path)},
            "road_occupancy": {
                "path": str(occupancy_path),
                "sha256": sha256(occupancy_path),
            },
        },
    }
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=manifest_path.parent, delete=False
    ) as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temp = Path(handle.name)
    os.replace(temp, manifest_path)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--weighted-catalogue", type=Path, required=True)
    parser.add_argument("--road-occupancy", type=Path, required=True)
    parser.add_argument("--output-nc", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output_nc.exists() or args.manifest.exists():
        raise FileExistsError("road-domain sample already exists")
    occupancy = dilated_road_mask(args.road_occupancy, JOINT_SUPPORT_KM)
    hits, n_trk, n_time = screen_tracks(
        args.tracks,
        occupancy["mask"],
        tile_deg=occupancy["tile_deg"],
        nlon=occupancy["nlon"],
        nlat=occupancy["nlat"],
    )
    selected = np.flatnonzero(hits).astype(np.int64)
    screen_record = {
        "track_count": int(n_trk),
        "hourly_nodes_per_track": int(n_time),
        "chunk_events": CHUNK_EVENTS,
        "occupied_tile_count": occupancy["occupied_tile_count"],
        "dilated_tile_count": occupancy["dilated_tile_count"],
        "radius_tiles": occupancy["radius_tiles"],
        "tile_deg": occupancy["tile_deg"],
        "selected_count": int(selected.size),
        "excluded_count": int(n_trk - selected.size),
    }
    manifest = write_sample(
        args.weighted_catalogue,
        selected,
        args.output_nc,
        args.manifest,
        tracks_path=args.tracks,
        occupancy_path=args.road_occupancy,
        screen_record=screen_record,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
