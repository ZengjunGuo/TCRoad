#!/usr/bin/env python3
"""Run the public-data Hurricane Irene C15--TCR reconstruction.

The preparation path is deliberately independent of CLIMADA so that the
frozen IBTrACS/NCEP inputs, the Knaff--Zehr RMW completion, the environmental
averages, and the moving 300-km grid can be audited before a rainfall run.

This is a method-faithful public reconstruction of Xi et al. (2020), not a
claim that the authors' original byte-identical inputs or private code are
available.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

import numpy as np
import xarray as xr
from netCDF4 import Dataset, chartostring, num2date


SID = "2011233N15301"
ATCF_ID = "AL092011"
TIME_START = "2011-08-21 00:00:00"
TIME_END = "2011-08-30 00:00:00"
N_SIX_HOURLY = 37
N_ONE_HOURLY = 217

EARTH_RADIUS_KM = 6371.0088
KT_TO_MPS = 0.514444
NM_TO_KM = 1.852

GRID_RESOLUTION_DEG = 0.05
MAX_DISTANCE_EYE_KM = 300.0
Q_DISK_RADIUS_KM = 200.0
SHEAR_ANNULUS_INNER_KM = 600.0
SHEAR_ANNULUS_OUTER_KM = 800.0

E_PRECIP = 0.9
LOWER_TROPOSPHERE_HEIGHT_M = 4000.0
RHO_AIR_OVER_RHO_LIQUID = 0.0012
MAX_W_FOREGROUND_MPS = 7.0
RADIAL_STEP_M = 2000.0
MIN_DRAG_COEFFICIENT = 0.001

FROZEN_ELEVATION_BYTES = 3_054_640
FROZEN_ELEVATION_SHA256 = (
    "de8142fe9f50d0cfbd944884ee945bb355b09dfce2f214879e608f87ae0f0951"
)
FROZEN_C_DRAG_BYTES = 999_184
FROZEN_C_DRAG_SHA256 = (
    "1c3f3b525f0c2a9e73f6fe6d3ba3caf7436699f8caeb3edb69548e04fd3f4a42"
)

# Required by CLIMADA's generic track schema. The C15 provider and TCR branch
# must not use this placeholder as a physical input (checked by the adapter).
CLIMADA_SCHEMA_ENVIRONMENTAL_PRESSURE_HPA = 1010.0

KZ_DOI = "10.1175/WAF965.1"
KZ_RMW_FORMULA = (
    "Rmax_km = 66.785 - 0.09102 * Vmax_kt "
    "+ 1.0619 * (abs(latitude_deg) - 25)"
)


def sha256(path: Path) -> str:
    """Return the SHA-256 of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_dump(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _decoded_strings(variable: Any, index: Any = slice(None)) -> np.ndarray:
    return np.asarray(chartostring(variable[index])).astype(str)


def _filled_float(variable: Any, index: Any) -> np.ndarray:
    values = np.ma.asarray(variable[index], dtype=float)
    return np.asarray(np.ma.filled(values, np.nan), dtype=float)


def _datetime64(text: str) -> np.datetime64:
    return np.datetime64(text.replace(" ", "T"), "ns")


def _iso_utc(value: np.datetime64) -> str:
    text = np.datetime_as_string(value.astype("datetime64[s]"), unit="s")
    return text + "Z"


def validate_frozen_package(input_dir: Path) -> dict[str, Any]:
    """Validate every frozen raw artifact against its manifest."""

    input_dir = input_dir.resolve()
    manifest_path = input_dir / "irene2011_inputs.manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"missing input manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "pass":
        raise ValueError("frozen input manifest does not have status='pass'")
    case = manifest.get("case", {})
    if case.get("sid") != SID or case.get("atcf_id") != ATCF_ID:
        raise ValueError(f"unexpected case identity in manifest: {case}")

    records: dict[str, Any] = {}
    for key, record in manifest.get("artifacts", {}).items():
        relative = Path(str(record["relative_path"]))
        path = (input_dir / relative).resolve()
        try:
            path.relative_to(input_dir)
        except ValueError as exc:
            raise ValueError(f"artifact path escapes input directory: {relative}") from exc
        if not path.is_file():
            raise FileNotFoundError(f"missing frozen artifact: {path}")
        actual_bytes = path.stat().st_size
        actual_sha = sha256(path)
        if actual_bytes != int(record["bytes"]):
            raise ValueError(
                f"byte count mismatch for {key}: {actual_bytes} != {record['bytes']}"
            )
        if actual_sha != record["sha256"]:
            raise ValueError(
                f"SHA-256 mismatch for {key}: {actual_sha} != {record['sha256']}"
            )
        records[key] = {
            "path": path,
            "relative_path": str(relative),
            "bytes": actual_bytes,
            "sha256": actual_sha,
            "source_url": record.get("source_url"),
        }

    expected = {
        "ibtracs_na",
        "ncep_r1_shum925",
        "ncep_r1_uwnd200",
        "ncep_r1_uwnd850",
        "ncep_r1_vwnd200",
        "ncep_r1_vwnd850",
    }
    if set(records) != expected:
        raise ValueError(
            f"unexpected frozen artifact set: missing={expected - set(records)}, "
            f"extra={set(records) - expected}"
        )
    return {
        "manifest_path": manifest_path,
        "manifest_sha256": sha256(manifest_path),
        "manifest": manifest,
        "artifacts": records,
    }


def extract_irene_six_hourly(ibtracs_path: Path) -> dict[str, Any]:
    """Extract the original USA/NHC six-hourly Irene sequence."""

    with Dataset(ibtracs_path) as data:
        sid = _decoded_strings(data.variables["sid"])
        positions = np.flatnonzero(sid == SID)
        if positions.size != 1:
            raise ValueError(f"expected one {SID}, found {positions.size}")
        storm_index = int(positions[0])

        time_text = _decoded_strings(data.variables["iso_time"], storm_index)
        selected: list[int] = []
        for index, text in enumerate(time_text):
            if not text:
                continue
            stamp = datetime.strptime(text, "%Y-%m-%d %H:%M:%S")
            if (
                datetime(2011, 8, 21) <= stamp <= datetime(2011, 8, 30)
                and stamp.hour in (0, 6, 12, 18)
                and stamp.minute == 0
                and stamp.second == 0
            ):
                selected.append(index)
        selected_array = np.asarray(selected, dtype=int)
        selected_text = time_text[selected_array]
        if (
            selected_array.size != N_SIX_HOURLY
            or selected_text[0] != TIME_START
            or selected_text[-1] != TIME_END
        ):
            raise ValueError(f"unexpected Irene standard-time sequence: {selected_text}")
        times = np.asarray([_datetime64(text) for text in selected_text])
        if not np.all(np.diff(times).astype("timedelta64[h]").astype(int) == 6):
            raise ValueError("Irene standard-time sequence is not exactly six-hourly")

        atcf_id = _decoded_strings(
            data.variables["usa_atcf_id"], (storm_index, selected_array)
        )
        basin = _decoded_strings(data.variables["basin"], (storm_index, selected_array))
        iflag = _decoded_strings(data.variables["iflag"], (storm_index, selected_array))
        status = _decoded_strings(
            data.variables["usa_status"], (storm_index, selected_array)
        )
        if not np.all(atcf_id == ATCF_ID):
            raise ValueError("Irene rows do not all have the expected USA ATCF identity")
        if not np.all(basin == "NA"):
            raise ValueError("Irene rows do not all have basin='NA'")
        if not np.all(np.char.startswith(iflag, "O")):
            raise ValueError("a selected IBTrACS row is not flagged as original")

        lat = _filled_float(data.variables["usa_lat"], (storm_index, selected_array))
        lon = _filled_float(data.variables["usa_lon"], (storm_index, selected_array))
        vmax_kt = _filled_float(
            data.variables["usa_wind"], (storm_index, selected_array)
        )
        central_pressure_hpa = _filled_float(
            data.variables["usa_pres"], (storm_index, selected_array)
        )
        rmw_raw_nmi = _filled_float(
            data.variables["usa_rmw"], (storm_index, selected_array)
        )

    for name, values in {
        "usa_lat": lat,
        "usa_lon": lon,
        "usa_wind": vmax_kt,
        "usa_pres": central_pressure_hpa,
    }.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} contains missing values in the frozen sequence")

    rmw_missing = ~np.isfinite(rmw_raw_nmi) | (rmw_raw_nmi <= 0)
    expected_missing = np.asarray([33, 34, 35, 36])
    if not np.array_equal(np.flatnonzero(rmw_missing), expected_missing):
        raise ValueError(
            "the frozen Irene RMW missing positions changed: "
            f"{np.flatnonzero(rmw_missing).tolist()}"
        )

    rmw_filled_nmi = rmw_raw_nmi.copy()
    kz_records: list[dict[str, Any]] = []
    for index in np.flatnonzero(rmw_missing):
        if vmax_kt[index] <= 15:
            raise ValueError("Knaff--Zehr Eq. (6) is only valid for Vmax > 15 kt")
        rmw_km = (
            66.785
            - 0.09102 * vmax_kt[index]
            + 1.0619 * (abs(lat[index]) - 25.0)
        )
        if not np.isfinite(rmw_km) or rmw_km <= 0:
            raise ValueError(f"invalid Knaff--Zehr RMW at index {index}: {rmw_km}")
        rmw_filled_nmi[index] = rmw_km / NM_TO_KM
        kz_records.append(
            {
                "six_hourly_index": int(index),
                "time_utc": _iso_utc(times[index]),
                "vmax_kt": float(vmax_kt[index]),
                "latitude_deg_north": float(lat[index]),
                "rmw_estimated_km": float(rmw_km),
                "rmw_estimated_nmi": float(rmw_filled_nmi[index]),
            }
        )

    return {
        "storm_index": storm_index,
        "time": times,
        "time_text": selected_text,
        "lat": lat,
        "lon": lon,
        "vmax_kt": vmax_kt,
        "central_pressure_hpa": central_pressure_hpa,
        "rmw_raw_nmi": rmw_raw_nmi,
        "rmw_filled_nmi": rmw_filled_nmi,
        "rmw_observed": ~rmw_missing,
        "rmw_estimated": rmw_missing,
        "iflag": iflag,
        "status": status,
        "basin": basin,
        "atcf_id": atcf_id,
        "knaff_zehr_records": kz_records,
    }


def _ncep_time_axis(data: Dataset) -> np.ndarray:
    variable = data.variables["time"]
    decoded = num2date(
        variable[:],
        units=variable.units,
        calendar=getattr(variable, "calendar", "standard"),
    )
    return np.asarray(
        [_datetime64(value.strftime("%Y-%m-%d %H:%M:%S")) for value in decoded]
    )


def _read_ncep_field(
    path: Path, variable_name: str, pressure_hpa: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(path) as data:
        required = {"time", "level", "lat", "lon", variable_name}
        if not required.issubset(data.variables):
            raise ValueError(f"missing variables in {path.name}: {required - set(data.variables)}")
        level = np.asarray(data.variables["level"][:], dtype=float)
        if level.shape != (1,) or not np.isclose(level[0], pressure_hpa):
            raise ValueError(f"unexpected pressure level in {path.name}: {level}")
        lat = np.asarray(data.variables["lat"][:], dtype=float)
        lon = np.asarray(data.variables["lon"][:], dtype=float)
        time = _ncep_time_axis(data)
        field = np.ma.asarray(data.variables[variable_name][:, 0, :, :], dtype=float)
        field = np.asarray(np.ma.filled(field, np.nan), dtype=float)
        if not np.all(np.isfinite(field)):
            raise ValueError(f"non-finite values in {path.name}:{variable_name}")
        expected_units = "kg/kg" if variable_name == "shum" else "m/s"
        if getattr(data.variables[variable_name], "units", None) != expected_units:
            raise ValueError(
                f"unexpected units in {path.name}:{variable_name}: "
                f"{getattr(data.variables[variable_name], 'units', None)!r}"
            )
    return time, lat, lon, field


def _great_circle_grid_km(
    center_lat_deg: float,
    center_lon_deg: float,
    grid_lat_deg: np.ndarray,
    grid_lon_deg: np.ndarray,
) -> np.ndarray:
    """Great-circle distance from one center to a rectilinear lat/lon grid."""

    center_lat = np.radians(center_lat_deg)
    lat = np.radians(np.asarray(grid_lat_deg, dtype=float))[:, None]
    delta_lat = lat - center_lat
    delta_lon_deg = (
        np.asarray(grid_lon_deg, dtype=float)[None, :] - center_lon_deg + 180.0
    ) % 360.0 - 180.0
    delta_lon = np.radians(delta_lon_deg)
    haversine = (
        np.sin(0.5 * delta_lat) ** 2
        + np.cos(center_lat) * np.cos(lat) * np.sin(0.5 * delta_lon) ** 2
    )
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(
        np.sqrt(np.clip(haversine, 0.0, 1.0))
    )


def extract_environmental_series(
    frozen: dict[str, Any], storm: dict[str, Any]
) -> dict[str, Any]:
    """Average NCEP fields at six-hourly centers using spherical masks."""

    artifacts = frozen["artifacts"]
    q_time, q_lat, q_lon, q_field = _read_ncep_field(
        artifacts["ncep_r1_shum925"]["path"], "shum", 925.0
    )
    u200_time, shear_lat, shear_lon, u200 = _read_ncep_field(
        artifacts["ncep_r1_uwnd200"]["path"], "uwnd", 200.0
    )
    u850_time, u850_lat, u850_lon, u850 = _read_ncep_field(
        artifacts["ncep_r1_uwnd850"]["path"], "uwnd", 850.0
    )
    v200_time, v200_lat, v200_lon, v200 = _read_ncep_field(
        artifacts["ncep_r1_vwnd200"]["path"], "vwnd", 200.0
    )
    v850_time, v850_lat, v850_lon, v850 = _read_ncep_field(
        artifacts["ncep_r1_vwnd850"]["path"], "vwnd", 850.0
    )

    for name, time in {
        "q925": q_time,
        "u200": u200_time,
        "u850": u850_time,
        "v200": v200_time,
        "v850": v850_time,
    }.items():
        if not np.array_equal(time, storm["time"]):
            raise ValueError(f"{name} times do not match the Irene six-hourly track")
    for name, (lat, lon) in {
        "u850": (u850_lat, u850_lon),
        "v200": (v200_lat, v200_lon),
        "v850": (v850_lat, v850_lon),
    }.items():
        if not np.array_equal(lat, shear_lat) or not np.array_equal(lon, shear_lon):
            raise ValueError(f"{name} grid does not match the u200 shear grid")

    q925 = np.empty(N_SIX_HOURLY, dtype=float)
    ushear = np.empty(N_SIX_HOURLY, dtype=float)
    vshear = np.empty(N_SIX_HOURLY, dtype=float)
    q_count = np.empty(N_SIX_HOURLY, dtype=np.int32)
    shear_count = np.empty(N_SIX_HOURLY, dtype=np.int32)
    u_difference = u200 - u850
    v_difference = v200 - v850

    for index, (center_lat, center_lon) in enumerate(
        zip(storm["lat"], storm["lon"], strict=True)
    ):
        q_distance = _great_circle_grid_km(center_lat, center_lon, q_lat, q_lon)
        q_mask = q_distance <= Q_DISK_RADIUS_KM
        shear_distance = _great_circle_grid_km(
            center_lat, center_lon, shear_lat, shear_lon
        )
        shear_mask = (
            (shear_distance >= SHEAR_ANNULUS_INNER_KM)
            & (shear_distance <= SHEAR_ANNULUS_OUTER_KM)
        )
        q_count[index] = int(np.count_nonzero(q_mask))
        shear_count[index] = int(np.count_nonzero(shear_mask))
        if q_count[index] == 0 or shear_count[index] == 0:
            raise ValueError(
                f"empty environmental averaging mask at {_iso_utc(storm['time'][index])}: "
                f"q_count={q_count[index]}, shear_count={shear_count[index]}"
            )
        # Xi et al. state that values are averaged within the disk/annulus,
        # without specifying latitude-area weighting. We therefore use the
        # literal equal-gridpoint arithmetic mean and record that choice.
        q925[index] = float(np.mean(q_field[index][q_mask]))
        ushear[index] = float(np.mean(u_difference[index][shear_mask]))
        vshear[index] = float(np.mean(v_difference[index][shear_mask]))

    for name, values in {"q925": q925, "ushear": ushear, "vshear": vshear}.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"non-finite derived environmental series: {name}")

    return {
        "q925_kgkg": q925,
        "ushear_mps": ushear,
        "vshear_mps": vshear,
        "q_gridpoint_count": q_count,
        "shear_gridpoint_count": shear_count,
        "spatial_weighting": "equal_gridpoint_arithmetic_mean",
        "distance_metric": "great_circle_haversine",
    }


def _numeric_time_hours(times: np.ndarray) -> np.ndarray:
    origin = times[0].astype("datetime64[s]")
    return (times.astype("datetime64[s]") - origin).astype("timedelta64[s]").astype(
        float
    ) / 3600.0


def _linear_interpolate(
    source_times: np.ndarray, source_values: np.ndarray, target_times: np.ndarray
) -> np.ndarray:
    return np.interp(
        _numeric_time_hours(target_times),
        _numeric_time_hours(source_times),
        np.asarray(source_values, dtype=float),
    )


def build_one_hourly_track(
    storm: dict[str, Any], environment: dict[str, Any]
) -> xr.Dataset:
    """Linearly interpolate storm and six-hourly environmental parameters to 1 h."""

    start = storm["time"][0].astype("datetime64[h]")
    end = storm["time"][-1].astype("datetime64[h]")
    time = np.arange(start, end + np.timedelta64(1, "h"), np.timedelta64(1, "h"))
    time = time.astype("datetime64[ns]")
    if time.size != N_ONE_HOURLY:
        raise ValueError(f"expected {N_ONE_HOURLY} one-hourly nodes, got {time.size}")

    unwrapped_lon = np.degrees(np.unwrap(np.radians(storm["lon"])))
    lon = _linear_interpolate(storm["time"], unwrapped_lon, time)
    lon = (lon + 180.0) % 360.0 - 180.0

    estimated_source = np.asarray(storm["rmw_estimated"], dtype=bool)
    source_hours = _numeric_time_hours(storm["time"])
    target_hours = _numeric_time_hours(time)
    right = np.searchsorted(source_hours, target_hours, side="left")
    right = np.clip(right, 0, source_hours.size - 1)
    exact = source_hours[right] == target_hours
    left = np.where(exact, right, np.maximum(right - 1, 0))
    rmw_uses_estimate = estimated_source[left] | estimated_source[right]

    vmax_kt = _linear_interpolate(storm["time"], storm["vmax_kt"], time)
    dataset = xr.Dataset(
        data_vars={
            "lat": ("time", _linear_interpolate(storm["time"], storm["lat"], time)),
            "lon": ("time", lon),
            "time_step": ("time", np.ones(time.size, dtype=float)),
            "max_sustained_wind": ("time", vmax_kt),
            "central_pressure": (
                "time",
                _linear_interpolate(
                    storm["time"], storm["central_pressure_hpa"], time
                ),
            ),
            "environmental_pressure": (
                "time",
                np.full(time.size, CLIMADA_SCHEMA_ENVIRONMENTAL_PRESSURE_HPA),
            ),
            "radius_max_wind": (
                "time",
                _linear_interpolate(
                    storm["time"], storm["rmw_filled_nmi"], time
                ),
            ),
            # CLIMADA-Petals names this compatibility slot q950. Its values
            # here are explicitly the Xi et al. 925-hPa disk means.
            "q950": (
                "time",
                _linear_interpolate(
                    storm["time"], environment["q925_kgkg"], time
                ),
            ),
            "ushear": (
                "time",
                _linear_interpolate(
                    storm["time"], environment["ushear_mps"], time
                ),
            ),
            "vshear": (
                "time",
                _linear_interpolate(
                    storm["time"], environment["vshear_mps"], time
                ),
            ),
            "basin": ("time", np.full(time.size, "NA", dtype="U2")),
            "rmw_uses_knaff_zehr_estimate": (
                "time",
                rmw_uses_estimate.astype(np.uint8),
            ),
        },
        coords={"time": time},
        attrs={
            "sid": SID,
            "atcf_id": ATCF_ID,
            "name": "IRENE",
            "orig_event_flag": 1,
            "category": _saffir_simpson_category(float(np.max(vmax_kt))),
            "max_sustained_wind_unit": "knots",
            "central_pressure_unit": "mbar",
            "radius_max_wind_unit": "nmile",
            "temporal_resolution_hours": 1,
            "environmental_pressure_schema_only_unused_by_c15_tcr": 1,
            "environmental_pressure_schema_placeholder_hpa": (
                CLIMADA_SCHEMA_ENVIRONMENTAL_PRESSURE_HPA
            ),
            "reconstruction_scope": "method_faithful_public_reconstruction",
        },
    )
    dataset["lat"].attrs.update(units="degrees_north", source="IBTrACS usa_lat")
    dataset["lon"].attrs.update(units="degrees_east", source="IBTrACS usa_lon")
    dataset["time_step"].attrs.update(units="h")
    dataset["max_sustained_wind"].attrs.update(
        units="knots", source="IBTrACS usa_wind"
    )
    dataset["central_pressure"].attrs.update(
        units="hPa", source="IBTrACS usa_pres"
    )
    dataset["environmental_pressure"].attrs.update(
        units="hPa",
        source="CLIMADA DEF_ENV_PRESSURE schema placeholder",
        schema_only_unused_by_c15_tcr=1,
    )
    dataset["radius_max_wind"].attrs.update(
        units="nmile",
        source="IBTrACS usa_rmw; missing values only: Knaff-Zehr Eq. (6)",
    )
    dataset["q950"].attrs.update(
        units="kg kg-1",
        source="NCEP/NCAR Reanalysis 1 shum at 925 hPa",
        source_pressure_hpa=925,
        compatibility_slot_name="q950",
        spatial_average="r <= 200 km at six-hourly centers before 1-h interpolation",
    )
    for component in ("ushear", "vshear"):
        dataset[component].attrs.update(
            units="m s-1",
            source="NCEP/NCAR Reanalysis 1 (200 hPa minus 850 hPa)",
            spatial_average=(
                "600 <= r <= 800 km at six-hourly centers before 1-h interpolation"
            ),
        )
    return dataset


def _saffir_simpson_category(maximum_wind_kt: float) -> int:
    if maximum_wind_kt < 34:
        return -1
    if maximum_wind_kt < 64:
        return 0
    if maximum_wind_kt < 83:
        return 1
    if maximum_wind_kt < 96:
        return 2
    if maximum_wind_kt < 113:
        return 3
    if maximum_wind_kt < 137:
        return 4
    return 5


def wrap_lon_deg(lon: np.ndarray | float) -> np.ndarray:
    """Map longitude to [-180, 180)."""

    return np.mod(np.asarray(lon, dtype=float) + 180.0, 360.0) - 180.0


def shortest_lon_delta_deg(lon: np.ndarray | float, center_lon: float) -> np.ndarray:
    """East-positive shortest longitude difference onto ``center_lon``."""

    return np.mod(np.asarray(lon, dtype=float) - float(center_lon) + 180.0, 360.0) - 180.0


def unwrap_longitude_deg(lon: np.ndarray) -> np.ndarray:
    """Make a track longitude series continuous across the date line."""

    wrapped = wrap_lon_deg(lon)
    if wrapped.size == 0:
        return wrapped
    out = np.empty_like(wrapped)
    out[0] = wrapped[0]
    for i in range(1, wrapped.size):
        out[i] = out[i - 1] + shortest_lon_delta_deg(wrapped[i], out[i - 1])
    return out


def build_moving_union_grid(track: xr.Dataset) -> xr.Dataset:
    """Create global-anchored 0.05-degree points within 300 km of any track node.

    Longitude is periodic. Tracks that cross the antimeridian are kept; the
    300 km / 0.05° union and great-circle test are unchanged.
    """

    resolution = GRID_RESOLUTION_DEG
    track_lat = np.asarray(track["lat"].values, dtype=float)
    track_lon = wrap_lon_deg(np.asarray(track["lon"].values, dtype=float))
    if track_lat.size == 0:
        raise ValueError("moving 300-km grid is empty")
    unwrapped_lon = unwrap_longitude_deg(track_lon)

    latitude_buffer = MAX_DISTANCE_EYE_KM / 110.0
    max_abs_buffered_lat = min(89.0, np.max(np.abs(track_lat)) + latitude_buffer)
    longitude_buffer = MAX_DISTANCE_EYE_KM / (
        110.0 * math.cos(math.radians(max_abs_buffered_lat))
    )
    lat_index_min = math.floor((np.min(track_lat) - latitude_buffer) / resolution)
    lat_index_max = math.ceil((np.max(track_lat) + latitude_buffer) / resolution)
    lon_index_min = math.floor((np.min(unwrapped_lon) - longitude_buffer) / resolution)
    lon_index_max = math.ceil((np.max(unwrapped_lon) + longitude_buffer) / resolution)
    latitude = np.arange(lat_index_min, lat_index_max + 1, dtype=int) * resolution
    longitude_unwrapped = np.arange(lon_index_min, lon_index_max + 1, dtype=int) * resolution
    longitude = wrap_lon_deg(longitude_unwrapped)
    union = np.zeros((latitude.size, longitude.size), dtype=bool)

    for center_lat, center_lon in zip(track_lat, track_lon, strict=True):
        local_lat_buffer = MAX_DISTANCE_EYE_KM / 110.0
        local_lon_buffer = MAX_DISTANCE_EYE_KM / (
            110.0
            * max(math.cos(math.radians(abs(center_lat) + local_lat_buffer)), 0.05)
        )
        lat_slice = np.flatnonzero(np.abs(latitude - center_lat) <= local_lat_buffer)
        lon_slice = np.flatnonzero(
            np.abs(shortest_lon_delta_deg(longitude, center_lon)) <= local_lon_buffer
        )
        if lat_slice.size == 0 or lon_slice.size == 0:
            continue
        distance = _great_circle_grid_km(
            center_lat,
            center_lon,
            latitude[lat_slice],
            longitude[lon_slice],
        )
        local_mask = distance <= MAX_DISTANCE_EYE_KM
        union[np.ix_(lat_slice, lon_slice)] |= local_mask

    row, column = np.nonzero(union)
    if row.size == 0:
        raise ValueError("moving 300-km grid is empty")
    centroid_lat = latitude[row]
    centroid_lon = wrap_lon_deg(longitude_unwrapped[column])
    lat_idx = np.rint(centroid_lat / resolution).astype(np.int32)
    lon_idx = np.rint(centroid_lon / resolution).astype(np.int32)
    _, unique = np.unique(np.stack([lat_idx, lon_idx], axis=1), axis=0, return_index=True)
    unique = np.sort(unique)
    centroid_lat = centroid_lat[unique]
    centroid_lon = centroid_lon[unique]
    lat_idx = lat_idx[unique]
    lon_idx = lon_idx[unique]
    grid = xr.Dataset(
        data_vars={
            "lat": ("centroid", centroid_lat),
            "lon": ("centroid", centroid_lon),
            "global_latitude_index": ("centroid", lat_idx),
            "global_longitude_index": ("centroid", lon_idx),
        },
        coords={"centroid": np.arange(centroid_lat.size, dtype=np.int32)},
        attrs={
            "grid_resolution_degrees": resolution,
            "maximum_distance_to_hourly_eye_km": MAX_DISTANCE_EYE_KM,
            "construction": (
                "union of global-anchored 0.05-degree points within 300 km "
                "great-circle distance of any one-hourly track center; "
                "longitude is periodic across the antimeridian"
            ),
        },
    )
    grid["lat"].attrs["units"] = "degrees_north"
    grid["lon"].attrs["units"] = "degrees_east"
    return grid


def build_six_hourly_dataset(
    storm: dict[str, Any], environment: dict[str, Any]
) -> xr.Dataset:
    rmw_raw = np.asarray(storm["rmw_raw_nmi"], dtype=float)
    dataset = xr.Dataset(
        data_vars={
            "lat": ("time", storm["lat"]),
            "lon": ("time", storm["lon"]),
            "vmax_kt": ("time", storm["vmax_kt"]),
            "central_pressure_hpa": ("time", storm["central_pressure_hpa"]),
            "rmw_raw_nmi": ("time", rmw_raw),
            "rmw_filled_nmi": ("time", storm["rmw_filled_nmi"]),
            "rmw_observed": ("time", storm["rmw_observed"].astype(np.uint8)),
            "rmw_knaff_zehr_estimated": (
                "time",
                storm["rmw_estimated"].astype(np.uint8),
            ),
            "q925_kgkg": ("time", environment["q925_kgkg"]),
            "ushear_200_minus_850_mps": ("time", environment["ushear_mps"]),
            "vshear_200_minus_850_mps": ("time", environment["vshear_mps"]),
            "q_disk_gridpoint_count": (
                "time",
                environment["q_gridpoint_count"],
            ),
            "shear_annulus_gridpoint_count": (
                "time",
                environment["shear_gridpoint_count"],
            ),
            "ibtracs_iflag": ("time", storm["iflag"]),
            "usa_status": ("time", storm["status"]),
        },
        coords={"time": storm["time"]},
        attrs={
            "sid": SID,
            "atcf_id": ATCF_ID,
            "ibtracs_storm_index": storm["storm_index"],
            "rmw_missing_completion": "Knaff-Zehr Eq. (6), missing rows only",
            "rmw_formula": KZ_RMW_FORMULA,
            "rmw_reference_doi": KZ_DOI,
            "environment_spatial_weighting": environment["spatial_weighting"],
            "environment_distance_metric": environment["distance_metric"],
        },
    )
    return dataset


def _netcdf_encoding(dataset: xr.Dataset) -> dict[str, dict[str, Any]]:
    encoding: dict[str, dict[str, Any]] = {}
    for name, variable in dataset.data_vars.items():
        if variable.dtype.kind not in {"U", "S", "O"}:
            encoding[name] = {"zlib": True, "complevel": 4, "shuffle": True}
    return encoding


def _write_netcdf(dataset: xr.Dataset, path: Path) -> None:
    dataset.to_netcdf(
        path,
        engine="netcdf4",
        format="NETCDF4",
        encoding=_netcdf_encoding(dataset),
    )


def _run_public_reconstruction(
    track: xr.Dataset,
    grid: xr.Dataset,
    elevation_tif: Path,
    c_drag_tif: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    """The single integration point for the reviewed C15--CLIMADA adapter."""

    try:
        from c15_climada_tcr import (  # type: ignore[import-not-found]
            assert_environmental_pressure_schema_only,
            run_tcr_public_reconstruction,
        )
    except ImportError as exc:
        raise RuntimeError(
            "the reviewed code/c15_climada_tcr.py adapter is not importable"
        ) from exc

    # Static/code-path contract check: no second science run is performed.
    assert_environmental_pressure_schema_only()
    result = run_tcr_public_reconstruction(
        track=track,
        centroid_lat=np.asarray(grid["lat"].values, dtype=float),
        centroid_lon=np.asarray(grid["lon"].values, dtype=float),
        elevation_tif=elevation_tif,
        c_drag_tif=c_drag_tif,
        e_precip=E_PRECIP,
        lower_troposphere_height_m=LOWER_TROPOSPHERE_HEIGHT_M,
        rho_air_over_rho_liquid=RHO_AIR_OVER_RHO_LIQUID,
        max_w_foreground=MAX_W_FOREGROUND_MPS,
        res_radial_m=RADIAL_STEP_M,
        min_c_drag=MIN_DRAG_COEFFICIENT,
        max_dist_eye_km=MAX_DISTANCE_EYE_KM,
    )
    if not isinstance(result, dict):
        raise TypeError("C15--TCR adapter result must be a mapping")
    rainrate = np.asarray(result.get("rainfall_rate_mm_h"), dtype=float)
    expected_shape = (track.sizes["time"], grid.sizes["centroid"])
    if rainrate.shape != expected_shape:
        raise ValueError(
            f"unexpected rain-rate shape: {rainrate.shape} != {expected_shape}"
        )
    if not np.all(np.isfinite(rainrate)) or np.any(rainrate < 0):
        raise ValueError("adapter returned non-finite or negative rain rates")
    metadata = result.get("metadata", {})
    if not isinstance(metadata, dict):
        raise TypeError("C15--TCR adapter metadata must be a mapping")
    return rainrate, metadata


def build_rainfall_output(
    track: xr.Dataset,
    grid: xr.Dataset,
    rainrate_mm_h: np.ndarray,
    adapter_metadata: dict[str, Any],
) -> xr.Dataset:
    rates = np.asarray(rainrate_mm_h, dtype=np.float32)
    time_step_hours = np.asarray(track["time_step"].values, dtype=float)
    if not np.all(time_step_hours == 1.0):
        raise ValueError("Irene reconstruction requires uniform one-hour time steps")
    # This is exactly the public Petals accumulation convention in
    # _compute_rain_sparse: sum(rainrate * time_step) over every track node.
    event_total = np.sum(
        rates * time_step_hours[:, None], axis=0, dtype=np.float64
    ).astype(np.float32)
    window = 24
    cumulative = np.concatenate(
        [
            np.zeros((1, rates.shape[1]), dtype=np.float64),
            np.cumsum(rates, axis=0, dtype=np.float64),
        ],
        axis=0,
    )
    rolling_24h = cumulative[window:] - cumulative[:-window]
    maximum_24h = np.max(rolling_24h, axis=0).astype(np.float32)
    dataset = xr.Dataset(
        data_vars={
            "rainfall_rate": (("time", "centroid"), rates),
            "event_total_rainfall": ("centroid", event_total),
            "maximum_24h_rainfall": ("centroid", maximum_24h),
            "lat": ("centroid", np.asarray(grid["lat"].values, dtype=float)),
            "lon": ("centroid", np.asarray(grid["lon"].values, dtype=float)),
        },
        coords={
            "time": track["time"].values,
            "centroid": grid["centroid"].values,
        },
        attrs={
            "sid": SID,
            "atcf_id": ATCF_ID,
            "model": "official C15 wind profile with public CLIMADA-Petals TCR skeleton",
            "scope": "raw method-faithful public reconstruction; no bias correction",
            "adapter_metadata_json": json.dumps(adapter_metadata, sort_keys=True),
        },
    )
    dataset["rainfall_rate"].attrs.update(units="mm h-1")
    dataset["event_total_rainfall"].attrs.update(
        units="mm",
        accumulation="full 2011-08-21T00:00Z to 2011-08-30T00:00Z lifecycle",
        numerical_semantics=(
            "CLIMADA-Petals node sum of rainfall_rate * time_step over all 217 nodes; "
            "public angular-wind implementation zeros the first profile row"
        ),
    )
    dataset["maximum_24h_rainfall"].attrs.update(
        units="mm", accumulation="maximum of all 24 consecutive one-hour rain-rate nodes"
    )
    dataset["lat"].attrs["units"] = "degrees_north"
    dataset["lon"].attrs["units"] = "degrees_east"
    return dataset


def _artifact_record(path: Path, output_root: Path) -> dict[str, Any]:
    return {
        "relative_path": str(path.relative_to(output_root)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def _validate_frozen_static_field(
    path: Path, *, label: str, expected_bytes: int, expected_sha256: str
) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"missing frozen {label} GeoTIFF: {resolved}")
    actual_bytes = resolved.stat().st_size
    actual_sha256 = sha256(resolved)
    if actual_bytes != expected_bytes or actual_sha256 != expected_sha256:
        raise ValueError(
            f"wrong {label} GeoTIFF identity: bytes={actual_bytes}, "
            f"sha256={actual_sha256}; expected bytes={expected_bytes}, "
            f"sha256={expected_sha256}"
        )
    return {
        "path": str(resolved),
        "bytes": actual_bytes,
        "sha256": actual_sha256,
        "frozen_identity_pass": True,
    }


def _base_manifest(
    frozen: dict[str, Any],
    storm: dict[str, Any],
    environment: dict[str, Any],
    track: xr.Dataset,
    grid: xr.Dataset,
    prepare_only: bool,
) -> dict[str, Any]:
    raw_inputs = {
        key: {
            "relative_path_in_frozen_package": record["relative_path"],
            "bytes": record["bytes"],
            "sha256": record["sha256"],
            "source_url": record["source_url"],
        }
        for key, record in frozen["artifacts"].items()
    }
    return {
        "schema_version": "1.0",
        "status": "prepared_only" if prepare_only else "completed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": {
            "name": "IRENE",
            "sid": SID,
            "atcf_id": ATCF_ID,
            "six_hourly_time_start_utc": _iso_utc(storm["time"][0]),
            "six_hourly_time_end_utc": _iso_utc(storm["time"][-1]),
            "six_hourly_time_count": int(storm["time"].size),
            "one_hourly_time_count": int(track.sizes["time"]),
        },
        "scope": {
            "description": "method-faithful public reconstruction of Xi et al. (2020)",
            "paper_identical_claim": False,
            "bias_correction_applied": False,
            "published_table_window_status": "PUBLISHED_TABLE_CONFLICT",
            "published_table_window_utc": [
                "2011-08-21T00:00:00Z",
                "2011-08-24T00:00:00Z",
            ],
            "published_table_conflict_note": (
                "Xi et al. (2020) Table 1 prints this Irene window, but it is "
                "incompatible with the North Carolina/Virginia landfall domain shown "
                "in Fig. 1. This runner therefore preserves the full frozen IBTrACS "
                "lifecycle and does not claim an exact Fig. 1 accumulation-window match."
            ),
        },
        "runner": {
            "relative_path": "TCRoad/code/run_irene_c15_climada.py",
            "sha256": sha256(Path(__file__).resolve()),
        },
        "frozen_input_manifest": {
            "filename": frozen["manifest_path"].name,
            "sha256": frozen["manifest_sha256"],
        },
        "raw_inputs": raw_inputs,
        "rmw_completion": {
            "observed_count": int(np.count_nonzero(storm["rmw_observed"])),
            "estimated_count": int(np.count_nonzero(storm["rmw_estimated"])),
            "estimated_only_where_missing": True,
            "time_extrapolation_used": False,
            "climada_estimate_rmw_used": False,
            "reference": "Knaff and Zehr, Weather and Forecasting 22, Eq. (6)",
            "doi": KZ_DOI,
            "formula": KZ_RMW_FORMULA,
            "records": storm["knaff_zehr_records"],
        },
        "environmental_preprocessing": {
            "operation_order": (
                "great-circle spatial average at each 6-hourly center, then linear "
                "interpolation of the resulting scalar/vector series to 1 hour"
            ),
            "q925": {
                "source_pressure_hpa": 925,
                "disk_radius_km": Q_DISK_RADIUS_KM,
                "gridpoint_count_range": [
                    int(np.min(environment["q_gridpoint_count"])),
                    int(np.max(environment["q_gridpoint_count"])),
                ],
                "climada_compatibility_slot": "q950",
            },
            "deep_layer_shear": {
                "definition": "(u200-u850, v200-v850)",
                "annulus_km": [
                    SHEAR_ANNULUS_INNER_KM,
                    SHEAR_ANNULUS_OUTER_KM,
                ],
                "gridpoint_count_range": [
                    int(np.min(environment["shear_gridpoint_count"])),
                    int(np.max(environment["shear_gridpoint_count"])),
                ],
            },
            "distance_metric": environment["distance_metric"],
            "spatial_weighting": environment["spatial_weighting"],
        },
        "grid": {
            "resolution_degrees": GRID_RESOLUTION_DEG,
            "moving_union_radius_km": MAX_DISTANCE_EYE_KM,
            "centroid_count": int(grid.sizes["centroid"]),
            "latitude_range": [
                float(grid["lat"].min()),
                float(grid["lat"].max()),
            ],
            "longitude_range": [
                float(grid["lon"].min()),
                float(grid["lon"].max()),
            ],
        },
        "tcr_parameters": {
            "precipitation_efficiency": E_PRECIP,
            "lower_troposphere_height_m": LOWER_TROPOSPHERE_HEIGHT_M,
            "rho_air_over_rho_liquid": RHO_AIR_OVER_RHO_LIQUID,
            "maximum_foreground_vertical_velocity_m_s": MAX_W_FOREGROUND_MPS,
            "radial_derivative_step_m": RADIAL_STEP_M,
            "minimum_drag_coefficient": MIN_DRAG_COEFFICIENT,
            "maximum_distance_to_eye_km": MAX_DISTANCE_EYE_KM,
        },
        "temporal_accumulation": {
            "event_total": (
                "CLIMADA-Petals _compute_rain_sparse convention: sum of each "
                "one-hour-node rain rate multiplied by its 1-h time_step; the public "
                "angular-wind path zeros the first profile row"
            ),
            "maximum_24h": "maximum sum over 24 consecutive one-hour rain-rate nodes",
        },
        "climada_schema_placeholder": {
            "environmental_pressure_hpa": CLIMADA_SCHEMA_ENVIRONMENTAL_PRESSURE_HPA,
            "schema_only_unused_by_c15_tcr": True,
            "not_xi_environmental_input": True,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or run the public-data Hurricane Irene C15--TCR reconstruction."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--elevation-tif", type=Path)
    parser.add_argument("--c-drag-tif", type=Path)
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="validate and preprocess inputs without importing or running CLIMADA",
    )
    args = parser.parse_args()
    if not args.prepare_only:
        if args.elevation_tif is None or args.c_drag_tif is None:
            parser.error(
                "science runs require explicit --elevation-tif and --c-drag-tif; "
                "there are no static-field defaults"
            )
        for label, path in {
            "elevation": args.elevation_tif,
            "drag coefficient": args.c_drag_tif,
        }.items():
            if not path.is_file():
                parser.error(f"{label} GeoTIFF does not exist: {path}")
    return args


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(
            f"output directory already exists; choose a new atomic target: {output_dir}"
        )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        frozen = validate_frozen_package(args.input_dir)
        storm = extract_irene_six_hourly(frozen["artifacts"]["ibtracs_na"]["path"])
        environment = extract_environmental_series(frozen, storm)
        track = build_one_hourly_track(storm, environment)
        grid = build_moving_union_grid(track)
        six_hourly = build_six_hourly_dataset(storm, environment)

        six_hourly_path = temporary / "irene_six_hourly_public_inputs.nc"
        track_path = temporary / "irene_one_hourly_climada_track.nc"
        grid_path = temporary / "irene_moving_300km_grid.nc"
        _write_netcdf(six_hourly, six_hourly_path)
        _write_netcdf(track, track_path)
        _write_netcdf(grid, grid_path)

        manifest = _base_manifest(
            frozen, storm, environment, track, grid, args.prepare_only
        )
        artifacts = {
            "six_hourly_inputs": _artifact_record(six_hourly_path, temporary),
            "one_hourly_track": _artifact_record(track_path, temporary),
            "moving_union_grid": _artifact_record(grid_path, temporary),
        }

        if not args.prepare_only:
            elevation_tif = args.elevation_tif.resolve()
            c_drag_tif = args.c_drag_tif.resolve()
            manifest["static_fields"] = {
                "elevation_tif": _validate_frozen_static_field(
                    elevation_tif,
                    label="topography_land_360as",
                    expected_bytes=FROZEN_ELEVATION_BYTES,
                    expected_sha256=FROZEN_ELEVATION_SHA256,
                ),
                "c_drag_tif": _validate_frozen_static_field(
                    c_drag_tif,
                    label="c_drag_500",
                    expected_bytes=FROZEN_C_DRAG_BYTES,
                    expected_sha256=FROZEN_C_DRAG_SHA256,
                ),
            }
            rainrate, adapter_metadata = _run_public_reconstruction(
                track, grid, elevation_tif, c_drag_tif
            )
            rainfall = build_rainfall_output(
                track, grid, rainrate, adapter_metadata
            )
            rainfall_path = temporary / "irene_c15_tcr_raw_rainfall.nc"
            _write_netcdf(rainfall, rainfall_path)
            artifacts["raw_rainfall"] = _artifact_record(rainfall_path, temporary)
            manifest["adapter_metadata"] = adapter_metadata
            manifest["raw_result_summary"] = {
                "maximum_rainfall_rate_mm_h": float(rainfall["rainfall_rate"].max()),
                "maximum_event_total_mm": float(
                    rainfall["event_total_rainfall"].max()
                ),
                "maximum_24h_mm": float(rainfall["maximum_24h_rainfall"].max()),
            }

        manifest["artifacts"] = artifacts
        manifest_path = temporary / "irene_c15_tcr_run.manifest.json"
        _json_dump(manifest_path, manifest)
        os.replace(temporary, output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": manifest["status"],
                    "output_dir": str(output_dir),
                    "manifest": str(output_dir / manifest_path.name),
                    "six_hourly_count": N_SIX_HOURLY,
                    "one_hourly_count": N_ONE_HOURLY,
                    "grid_centroid_count": int(grid.sizes["centroid"]),
                    "rmw_knaff_zehr_estimated_count": int(
                        np.count_nonzero(storm["rmw_estimated"])
                    ),
                },
                indent=2,
                sort_keys=True,
            )
        )
    finally:
        if temporary is not None and temporary.exists():
            shutil.rmtree(temporary)


if __name__ == "__main__":
    main()
