#!/usr/bin/env python3
"""Prepare one frozen Lin event for the public C15--CLIMADA TCR chain.

This program is intentionally *prepare-only*.  It resolves every synthetic
environment quantity that is explicitly available from the frozen Lin and
Emanuel v6.4 public code, records unavailable quantities as blockers, and
never calls CLIMADA or computes a hazard field.

The event is recovered from the native one-hour 100,000-track catalogue by
``source_track_index``.  A hash-linked subset of that track may be supplied for
local audit runs; it is not a replacement catalogue.

Following Gori et al. (2022), one outer radius of vanishing wind is read from
the immutable event-level catalogue and held fixed for the storm lifetime.
The official C15 ``r0input`` solver combines that radius with each hourly Lin
``v_trks`` circular intensity and Coriolis magnitude to infer the hourly RMW.
No alternative size predictor, splice, clipping, or replacement is reachable
from this production prepare path.
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
import tempfile
from typing import Any

from netCDF4 import Dataset, chartostring, date2num, num2date
import numpy as np
from scipy.interpolate import RectBivariateSpline

from c15_climada_tcr import C15FixedR0WindProfileProvider, C15_SOURCE_DOI


SCRIPT_VERSION = "2.0.0"

RMW_PROVIDER_DECISION_STATUS = "FROZEN_GORI2022_CHAVAS2016_C15_R0INPUT"

EXPECTED_EVENT_POSITION = 0
EXPECTED_EVENT_ID = "stream0000-year1995-track000002"
EXPECTED_SOURCE_TRACK_INDEX = 2
EXPECTED_NATIVE_START = 56
EXPECTED_NATIVE_STOP = 149
EXPECTED_TARGET_COUNT = 94

FROZEN_SAMPLE_SHA256 = (
    "856ed368466cf4f8a1f0b8e351bcc8f44eae32d9d60c55623c6ad2217275d1af"
)
FROZEN_FULL_TRACK_SHA256 = (
    "5cccb10168df5a15144c1ae3bb97c1533c8722c03d2e5fecb599fd7496258fae"
)
FROZEN_CMIP6_TA_SHA256 = (
    "11d146367fc3f22588d265e2ccd79b1634693bf7c9840b0e44bce2afde33a10d"
)

MPS_TO_KNOTS = 3600.0 / 1852.0
KNOTS_TO_MPS = 1852.0 / 3600.0
TRANSLATION_SMOOTHING_FACTOR = 0.4
EMANUEL_VDRIFT_MPS = 1.5

GORI2022_DOI = "10.1038/s41558-021-01272-7"
CHAVAS2016_DOI = "10.1175/JCLI-D-15-0731.1"
EARTH_ANGULAR_VELOCITY_S = 7.2921159e-5
EMANUEL_V64_SOURCE = "vendor/Emanuel_TCR/scripts_ver6.4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def nc_strings(variable: Any) -> np.ndarray:
    values = variable[:]
    if values.dtype.kind in {"U", "O"}:
        return np.asarray(values, dtype=str)
    if values.dtype.kind == "S" and values.ndim == 2:
        return np.asarray(chartostring(values), dtype=str)
    return np.asarray(values).astype(str)


def nc_scalar(variable: Any, index: int) -> Any:
    value = variable[index]
    if np.ma.is_masked(value):
        raise ValueError(f"masked scalar in {variable.name} at event {index}")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8").strip()
    return value


def validate_identity(path: Path, expected_sha: str, label: str) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"missing {label}: {path}")
    actual = sha256(path)
    if actual != expected_sha:
        raise ValueError(
            f"{label} SHA-256 mismatch: expected {expected_sha}, got {actual}"
        )
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": actual}


def load_event_identity(
    sample_path: Path,
    event_position: int,
    *,
    enforce_event0: bool = True,
    expected_sample_sha256: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    sample_record = validate_identity(
        sample_path,
        FROZEN_SAMPLE_SHA256 if expected_sample_sha256 is None else expected_sample_sha256,
        "event sample",
    )
    with Dataset(sample_path) as sample:
        if enforce_event0 and event_position != EXPECTED_EVENT_POSITION:
            raise ValueError(
                f"this candidate is frozen to event position {EXPECTED_EVENT_POSITION}"
            )
        event_id = str(nc_scalar(sample.variables["event_id"], event_position))
        if enforce_event0 and event_id != EXPECTED_EVENT_ID:
            raise ValueError(f"unexpected event identity: sample={event_id!r}")
        source_track_index = int(
            nc_scalar(sample.variables["source_track_index"], event_position)
        )
        if enforce_event0 and source_track_index != EXPECTED_SOURCE_TRACK_INDEX:
            raise ValueError("sample source-track identity changed")

        native_start = int(
            nc_scalar(
                sample.variables["threshold_genesis_native_index"], event_position
            )
        )
        native_stop = int(
            nc_scalar(
                sample.variables["threshold_lysis_native_index"], event_position
            )
        )
        if enforce_event0 and (native_start, native_stop) != (
            EXPECTED_NATIVE_START,
            EXPECTED_NATIVE_STOP,
        ):
            raise ValueError("event-0 native threshold window changed")

        identity = {
            "event_position": event_position,
            "event_id": event_id,
            "source_track_index": source_track_index,
            "source_catalogue_event_position": int(
                nc_scalar(
                    sample.variables["source_catalogue_event_position"],
                    event_position,
                )
            ),
            "task_year": int(
                nc_scalar(sample.variables["task_year"], event_position)
            ),
            "seed_month": int(
                nc_scalar(sample.variables["lin_seed_month"], event_position)
            ),
            "seed_basin": str(
                nc_scalar(sample.variables["lin_seed_basin"], event_position)
            ),
            "threshold_genesis_region": str(
                nc_scalar(
                    sample.variables["threshold_genesis_region"], event_position
                )
            ),
            "threshold_genesis_datetime": str(
                nc_scalar(
                    sample.variables["threshold_genesis_datetime"], event_position
                )
            ),
            "threshold_lysis_datetime": str(
                nc_scalar(
                    sample.variables["threshold_lysis_datetime"], event_position
                )
            ),
            "native_start": native_start,
            "native_stop": native_stop,
            "event_weight_climate_fixed_effect_ht_analysis_yr": float(
                nc_scalar(
                    sample.variables[
                        "event_weight_climate_fixed_effect_ht_analysis_yr"
                    ],
                    event_position,
                )
            ),
        }
    return identity, {
        "sample": sample_record,
        "size_contract_note": (
            "One immutable event-level r0 is supplied by the frozen catalogue; "
            "official C15 r0input infers RMW at every event node."
        ),
    }


def load_fixed_r0_for_event(
    catalogue_path: Path,
    catalogue_manifest_path: Path,
    identity: dict[str, Any],
    *,
    expected_sample_sha256: str = FROZEN_SAMPLE_SHA256,
) -> tuple[float, dict[str, Any]]:
    """Load one exact event-position/id match from the immutable r0 catalogue."""

    catalogue_path = catalogue_path.resolve()
    catalogue_manifest_path = catalogue_manifest_path.resolve()
    if not catalogue_path.is_file() or not catalogue_manifest_path.is_file():
        raise FileNotFoundError("fixed-r0 catalogue NetCDF and manifest are required")
    manifest = json.loads(catalogue_manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "FROZEN_IMMUTABLE":
        raise ValueError("fixed-r0 catalogue manifest is not immutable and complete")
    artifact = manifest.get("artifacts", {}).get("fixed_r0_catalogue_netcdf", {})
    catalogue_sha = sha256(catalogue_path)
    if (
        artifact.get("sha256") != catalogue_sha
        or int(artifact.get("bytes", -1)) != catalogue_path.stat().st_size
    ):
        raise ValueError("fixed-r0 catalogue identity differs from its manifest")
    source_sample = manifest.get("source_sample", {})
    if source_sample.get("sha256") != expected_sample_sha256:
        raise ValueError("fixed-r0 catalogue is not linked to the frozen sample")

    position = int(identity["event_position"])
    event_id = str(identity["event_id"])
    with Dataset(catalogue_path) as catalogue:
        required = {"event_position", "event_id", "outer_radius_m"}
        missing = required - set(catalogue.variables)
        if missing:
            raise ValueError(f"fixed-r0 catalogue lacks variables: {sorted(missing)}")
        if str(getattr(catalogue, "status", "")) != "FROZEN_IMMUTABLE":
            raise ValueError("fixed-r0 catalogue NetCDF is not immutable")
        if str(getattr(catalogue, "source_sample_sha256", "")) != expected_sample_sha256:
            raise ValueError("fixed-r0 catalogue NetCDF has the wrong source sample")
        count = len(catalogue.dimensions["event"])
        if not 0 <= position < count:
            raise IndexError(f"event position {position} is outside fixed-r0 catalogue")
        stored_position = int(nc_scalar(catalogue.variables["event_position"], position))
        stored_id = str(nc_scalar(catalogue.variables["event_id"], position))
        outer_radius_m = float(
            nc_scalar(catalogue.variables["outer_radius_m"], position)
        )
        attrs: dict[str, Any] = {}
        for name in (
                "schema_version",
                "scientific_source_doi",
                "scientific_source_table",
                "distribution",
                "lognormal_mu_ln_km",
                "lognormal_sigma_ln_km",
                "distribution_contract_sha256",
                "rng_bit_generator",
                "rng_seed_decimal",
                "draw_order",
                "outer_radius_m_sequence_sha256",
                "event_outer_radius_binding_sha256",
                "truncation_applied",
                "rejection_or_resampling_applied",
                "clipping_applied",
            ):
            if hasattr(catalogue, name):
                value = getattr(catalogue, name)
                attrs[name] = value.item() if isinstance(value, np.generic) else value
    if stored_position != position or stored_id != event_id:
        raise ValueError(
            "fixed-r0 catalogue identity mismatch: "
            f"requested ({position}, {event_id!r}), stored "
            f"({stored_position}, {stored_id!r})"
        )
    if not np.isfinite(outer_radius_m) or outer_radius_m <= 0.0:
        raise ValueError("fixed-r0 catalogue returned a non-positive/non-finite radius")
    return outer_radius_m, {
        "method": "Gori2022 event-fixed outer radius with official C15 r0input",
        "gori2022_doi": GORI2022_DOI,
        "chavas2016_doi": CHAVAS2016_DOI,
        "nature_source_identical_sampler_claim": False,
        "catalogue": {
            "path": str(catalogue_path),
            "bytes": catalogue_path.stat().st_size,
            "sha256": catalogue_sha,
        },
        "catalogue_manifest": {
            "path": str(catalogue_manifest_path),
            "bytes": catalogue_manifest_path.stat().st_size,
            "sha256": sha256(catalogue_manifest_path),
        },
        "event_position": position,
        "event_id": event_id,
        "outer_radius_m": outer_radius_m,
        "outer_radius_km": outer_radius_m / 1000.0,
        "catalogue_attributes": attrs,
        "truncation_applied": False,
        "rejection_or_resampling_applied": False,
        "clipping_applied": False,
    }


def _filled_vector(variable: Any, selector: Any = slice(None)) -> np.ndarray:
    return np.asarray(
        np.ma.filled(variable[selector], np.nan), dtype=np.float64
    )


def load_track_window(track_path: Path, identity: dict[str, Any]) -> dict[str, Any]:
    """Load all finite native nodes and identify the target event window."""

    track_path = track_path.resolve()
    if not track_path.is_file():
        raise FileNotFoundError(track_path)
    file_sha = sha256(track_path)
    start = int(identity["native_start"])
    stop = int(identity["native_stop"])
    required_variables = (
        "lon_trks",
        "lat_trks",
        "v_trks",
        "vmax_trks",
        "u250_trks",
        "v250_trks",
        "u850_trks",
        "v850_trks",
        "time",
    )
    with Dataset(track_path) as tracks:
        missing = [name for name in required_variables if name not in tracks.variables]
        if missing:
            raise ValueError(f"track input is missing variables: {missing}")

        is_full_catalogue = "n_trk" in tracks.dimensions
        if is_full_catalogue:
            if file_sha != FROZEN_FULL_TRACK_SHA256:
                raise ValueError("full 100,000-track catalogue SHA-256 changed")
            source_index = int(identity["source_track_index"])
            native_available = np.arange(len(tracks.dimensions["time"]), dtype=np.int64)

            def read(name: str) -> np.ndarray:
                variable = tracks.variables[name]
                selector = (
                    (source_index, slice(None))
                    if variable.ndim == 2
                    else slice(None)
                )
                return _filled_vector(variable, selector)

            source_full_sha = file_sha
            snapshot_source = False
            first_finite_lat = float(read("lat_trks")[np.flatnonzero(np.isfinite(read("lat_trks")))[0]])
        else:
            source_full_sha = str(
                getattr(tracks, "source_full_track_sha256", "")
            )
            if source_full_sha != FROZEN_FULL_TRACK_SHA256:
                raise ValueError(
                    "local track snapshot is not linked to the frozen full catalogue"
                )
            if int(getattr(tracks, "source_track_index", -1)) != int(
                identity["source_track_index"]
            ):
                raise ValueError("local track snapshot has the wrong source track")
            if "native_index" not in tracks.variables:
                raise ValueError("local track snapshot lacks native_index")
            native_available = np.asarray(
                tracks.variables["native_index"][:], dtype=np.int64
            )

            def read(name: str) -> np.ndarray:
                return _filled_vector(tracks.variables[name])

            snapshot_source = True
            first_finite_lat = float(getattr(tracks, "first_finite_lat_deg"))

        if not np.all(np.diff(native_available) == 1):
            raise ValueError("track native indices are not a contiguous one-hour sequence")
        source_latitude = read("lat_trks")
        finite_positions = np.flatnonzero(np.isfinite(source_latitude))
        if finite_positions.size < 3 or not np.all(np.diff(finite_positions) == 1):
            raise ValueError("source track needs at least three contiguous finite nodes")
        required_native = native_available[
            finite_positions[0] : finite_positions[-1] + 1
        ]
        if start < required_native[0] or stop > required_native[-1]:
            raise ValueError("event target window lies outside finite source track")
        positions = np.searchsorted(native_available, required_native)
        if (
            np.any(positions >= native_available.size)
            or not np.array_equal(native_available[positions], required_native)
        ):
            raise ValueError("track input lacks required target translation context")

        arrays: dict[str, np.ndarray] = {}
        for name in required_variables[:-1]:
            values = read(name)
            arrays[name] = values[positions]
            if not np.all(np.isfinite(arrays[name])):
                raise ValueError(f"non-finite {name} in target/context window")
        time_all = read("time")
        arrays["time"] = time_all[positions]
        if not np.all(np.isfinite(arrays["time"])):
            raise ValueError("non-finite native time values")
        if not np.allclose(np.diff(arrays["time"]), 3600.0, rtol=0, atol=1e-6):
            raise ValueError("Lin event is not native one-hour data")

    target_selector = (required_native >= start) & (required_native <= stop)
    target_native = required_native[target_selector]
    if target_native.size != stop - start + 1:
        raise ValueError("unexpected target count")
    return {
        "path": str(track_path),
        "bytes": track_path.stat().st_size,
        "sha256": file_sha,
        "source_full_track_sha256": source_full_sha,
        "is_hash_linked_local_snapshot": snapshot_source,
        "first_finite_lat_deg": first_finite_lat,
        "native_index": required_native,
        "target_selector": target_selector,
        "target_native_index": target_native,
        **arrays,
    }


def official_utrans_adapted_to_native_hour(
    track: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Mechanically apply v6.4 ``utrans`` to a uniform one-hour Lin track.

    The public routine uses a centred two-node span and one second-difference
    smoothing pass.  Only the dimensional factor changes from its native 2 h
    event-set spacing to the documented 1 h Lin spacing.  No 2 h downsampling
    or new filter is introduced.
    """

    lon = np.asarray(track["lon_trks"], dtype=float)
    lat = np.asarray(track["lat_trks"], dtype=float)
    native = np.asarray(track["native_index"], dtype=np.int64)
    target = np.asarray(track["target_native_index"], dtype=np.int64)
    dt_hours = float(np.median(np.diff(np.asarray(track["time"], dtype=float))) / 3600)
    if not math.isclose(dt_hours, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("this candidate is frozen to native one-hour Lin tracks")

    delta_lon = ((lon[2:] - lon[:-2] + 180.0) % 360.0) - 180.0
    factor = 60.0 / (2.0 * dt_hours)
    if lon.size < 3:
        raise ValueError("official utrans requires at least three native nodes")
    raw_u = np.full(lon.shape, np.nan, dtype=float)
    raw_v = np.full(lon.shape, np.nan, dtype=float)
    raw_u[1:-1] = factor * np.cos(np.deg2rad(lat[1:-1])) * delta_lon
    raw_v[1:-1] = factor * (lat[2:] - lat[:-2])

    # Exact endpoint extrapolation in Emanuel v6.4 utrans.m.  It is dormant
    # for the reviewed event-0 interior window, but makes native index 0 and
    # the final available 15-day node computable without invented padding.
    raw_u[0] = 2.0 * raw_u[1] - raw_u[min(2, lon.size - 1)]
    raw_v[0] = 2.0 * raw_v[1] - raw_v[min(2, lon.size - 1)]
    raw_u[-1] = 2.0 * raw_u[-2] - raw_u[max(lon.size - 3, 0)]
    raw_v[-1] = 2.0 * raw_v[-2] - raw_v[max(lon.size - 3, 0)]

    smooth_u = raw_u.copy()
    smooth_v = raw_v.copy()
    smooth_u[1:-1] = raw_u[1:-1] + TRANSLATION_SMOOTHING_FACTOR * (
        raw_u[:-2] + raw_u[2:] - 2.0 * raw_u[1:-1]
    )
    smooth_v[1:-1] = raw_v[1:-1] + TRANSLATION_SMOOTHING_FACTOR * (
        raw_v[:-2] + raw_v[2:] - 2.0 * raw_v[1:-1]
    )
    target_positions = np.searchsorted(native, target)
    output_u = smooth_u[target_positions]
    output_v = smooth_v[target_positions]
    if not np.all(np.isfinite(output_u)) or not np.all(np.isfinite(output_v)):
        raise ValueError("official utrans stencil is incomplete for target window")
    return output_u, output_v


def emanuel_v64_shear(
    translation_u_knots: np.ndarray,
    translation_v_knots: np.ndarray,
    u850_ms: np.ndarray,
    v850_ms: np.ndarray,
    target_lat_deg: np.ndarray,
    first_finite_lat_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    hemisphere_sign = float(np.sign(first_finite_lat_deg))
    if hemisphere_sign == 0.0:
        raise ValueError("v6.4 drift sign is undefined for an equatorial first node")
    u850_knots = np.asarray(u850_ms, dtype=float) * MPS_TO_KNOTS
    v850_knots = np.asarray(v850_ms, dtype=float) * MPS_TO_KNOTS
    drift_knots = EMANUEL_VDRIFT_MPS * MPS_TO_KNOTS * hemisphere_sign
    ushear = 5.0 * KNOTS_TO_MPS * (
        np.asarray(translation_u_knots, dtype=float) - u850_knots
    )
    vshear = 5.0 * KNOTS_TO_MPS * (
        np.asarray(translation_v_knots, dtype=float)
        - drift_knots * np.cos(np.deg2rad(target_lat_deg))
        - v850_knots
    )
    return ushear, vshear


def sample_t600_lin_public_adapter(
    ta_path: Path,
    target_lon_deg: np.ndarray,
    target_lat_deg: np.ndarray,
    task_year: int,
    seed_month: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Use Lin's frozen month-15 temporal and kx=ky=1 spatial interpolation.

    Lin v1.1 does not persist ``T600`` in its track output.  This adapter uses
    the same public interpolation pattern Lin applies to its monthly
    thermodynamic fields: interpolate the monthly GCM field to
    ``datetime(year, month, 15)`` once, then bilinearly evaluate that fixed
    monthly field along the event.  It must not be described as unpublished
    Nature production code.
    """

    ta_record = validate_identity(ta_path, FROZEN_CMIP6_TA_SHA256, "CMIP6 ta input")
    with Dataset(ta_path) as source:
        required = {"ta", "time", "plev", "lat", "lon"}
        if not required.issubset(source.variables):
            raise ValueError(f"CMIP6 ta input lacks {required - set(source.variables)}")
        ta = source.variables["ta"]
        if getattr(ta, "units", None) != "K":
            raise ValueError(f"T600 source units must be K, got {getattr(ta, 'units', None)!r}")
        plev_var = source.variables["plev"]
        plev = np.asarray(plev_var[:], dtype=float)
        plev_units = str(getattr(plev_var, "units", ""))
        plev_hpa = plev / 100.0 if plev_units == "Pa" else plev
        if plev_units not in {"Pa", "hPa", "millibars"}:
            raise ValueError(f"unsupported pressure units: {plev_units!r}")
        pressure_matches = np.flatnonzero(np.isclose(plev_hpa, 600.0))
        if pressure_matches.size != 1:
            raise ValueError("CMIP6 ta input does not have one exact 600-hPa level")
        pressure_index = int(pressure_matches[0])

        time_var = source.variables["time"]
        calendar = str(getattr(time_var, "calendar", "standard"))
        units = str(time_var.units)
        time_numeric = np.asarray(time_var[:], dtype=float)
        target_datetime = datetime(task_year, seed_month, 15)
        target_numeric = float(date2num(target_datetime, units=units, calendar=calendar))
        right = int(np.searchsorted(time_numeric, target_numeric, side="left"))
        exact = right < time_numeric.size and math.isclose(
            time_numeric[right], target_numeric, rel_tol=0, abs_tol=1e-12
        )
        # An exact first or last time-axis value is valid.  The old ordering
        # rejected an exact January boundary before testing equality.
        boundary_month_field = False
        if exact:
            left = right
            weight_right = 0.0
        else:
            if right == 0:
                first_date = num2date(
                    time_numeric[0], units=units, calendar=calendar
                )
                if (int(first_date.year), int(first_date.month)) != (
                    task_year,
                    seed_month,
                ):
                    raise ValueError(
                        "month-15 T600 target precedes a different source month"
                    )
                # The CMIP monthly coordinate is the source interval midpoint
                # (1995-01-16 12:00 here), whereas Lin asks for the January
                # month-15 field.  Use that same January monthly field at the
                # truncated dataset boundary; do not extrapolate into 1994.
                left = right = 0
                weight_right = 0.0
                boundary_month_field = True
            elif right >= time_numeric.size:
                last_date = num2date(
                    time_numeric[-1], units=units, calendar=calendar
                )
                if (int(last_date.year), int(last_date.month)) != (
                    task_year,
                    seed_month,
                ):
                    raise ValueError(
                        "month-15 T600 target follows a different source month"
                    )
                left = right = time_numeric.size - 1
                weight_right = 0.0
                boundary_month_field = True
            else:
                left = right - 1
                weight_right = (target_numeric - time_numeric[left]) / (
                    time_numeric[right] - time_numeric[left]
                )
        field_left = np.asarray(
            np.ma.filled(ta[left, pressure_index, :, :], np.nan), dtype=float
        )
        if left == right:
            field = field_left
        else:
            field_right = np.asarray(
                np.ma.filled(ta[right, pressure_index, :, :], np.nan), dtype=float
            )
            field = (1.0 - weight_right) * field_left + weight_right * field_right
        lat = np.asarray(source.variables["lat"][:], dtype=float)
        lon = np.asarray(source.variables["lon"][:], dtype=float)

    if not np.all(np.isfinite(field)):
        raise ValueError("interpolated 600-hPa monthly field contains non-finite values")
    if lat[1] < lat[0]:
        lat = lat[::-1]
        field = field[::-1, :]
    if not np.all(np.diff(lat) > 0) or not np.all(np.diff(lon) > 0):
        raise ValueError("T600 source grid must be rectilinear and increasing")
    if lon.size < 2 or lon[-1] - lon[0] >= 360.0:
        raise ValueError("T600 longitude axis is not a non-duplicated periodic grid")
    # Rectilinear periodic extension of the same monthly field.  This changes
    # longitude indexing only, not the interpolated field or physical method.
    lon_extended = np.concatenate(([lon[-1] - 360.0], lon, [lon[0] + 360.0]))
    field_extended = np.concatenate(
        (field[:, -1:], field, field[:, :1]), axis=1
    )
    target_lon = (
        (np.asarray(target_lon_deg, dtype=float) - lon[0]) % 360.0 + lon[0]
    )
    interpolator = RectBivariateSpline(
        lon_extended, lat, field_extended.T, kx=1, ky=1
    )
    t600 = np.asarray(
        interpolator.ev(target_lon, np.asarray(target_lat_deg, dtype=float)),
        dtype=float,
    )
    if not np.all(np.isfinite(t600)) or np.any(t600 <= 100.0):
        raise ValueError("invalid event T600 series")
    metadata = {
        **ta_record,
        "source_variable": "ta",
        "source_pressure_hpa": 600.0,
        "source_units": "K",
        "fixed_environment_datetime": target_datetime.isoformat() + "Z",
        "time_interpolation": "linear to datetime(task_year, seed_month, 15)",
        "time_left_index": left,
        "time_right_index": right,
        "time_right_weight": float(weight_right),
        "boundary_month_field_used": boundary_month_field,
        "space_interpolation": (
            "periodic-longitude RectBivariateSpline kx=1 ky=1, as frozen Lin "
            "mat.interp2_fx"
        ),
        "held_fixed_monthly_field_through_event": True,
        "claim_boundary": "public-code-derived adapter; not Nature production source identity",
    }
    return t600, metadata


def emanuel_v64_qs950(
    t600_k: np.ndarray, circular_wind_ms: np.ndarray
) -> np.ndarray:
    """Exact vector form of official Emanuel v6.4 ``qs900b.m``."""

    t600 = np.asarray(t600_k, dtype=float)
    # Preserve the public interface: raingen supplies circular wind in knots;
    # qs900b converts it back to m/s before the warm-core entropy term.
    vmax_knots = np.asarray(circular_wind_ms, dtype=float) * MPS_TO_KNOTS
    vmax_ms = vmax_knots * KNOTS_TO_MPS
    cp = 1005.0
    rv = 491.0
    rd = 287.0
    lv = 2.5e6
    pref_hpa = 950.0
    tc600 = np.maximum(t600 - 273.15, -50.0)
    es600 = 6.112 * np.exp(17.67 * tc600 / (243.5 + tc600))
    q600 = 0.622 * es600 / (600.0 - es600)
    temperature = t600 + 20.0
    qs950 = np.zeros_like(t600)
    for _ in range(5):
        tc = temperature - 273.15
        es = 6.122 * np.exp(17.67 * tc / (243.5 + tc))
        qs950 = 0.622 * es / (pref_hpa - es)
        equation_error = (
            cp * np.log(temperature / t600)
            + lv * (qs950 / temperature - q600 / t600)
            - rd * np.log(pref_hpa / 600.0)
            - 0.016 * vmax_ms**2
        )
        derivative = (
            cp * temperature
            + lv * qs950 * ((lv / rv) / temperature - 1.0)
        ) / temperature**2
        temperature = temperature - equation_error / derivative
    # MATLAB assigns the qs value computed at the start of the fifth Newton
    # iteration, before updating T for a sixth evaluation.
    if not np.all(np.isfinite(qs950)) or np.any(qs950 <= 0.0):
        raise ValueError("official qs900b reconstruction returned invalid q950")
    return qs950


def derive_c15_r0input_rmw(
    circular_wind_ms: np.ndarray,
    latitude_deg: np.ndarray,
    outer_radius_m: float,
    *,
    provider: C15FixedR0WindProfileProvider | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Infer every hourly RMW from official C15 with one fixed event r0."""

    intensity = np.asarray(circular_wind_ms, dtype=float)
    latitude = np.asarray(latitude_deg, dtype=float)
    if intensity.shape != latitude.shape or intensity.ndim != 1:
        raise ValueError("C15 intensity and latitude must share one-dimensional shape")
    if not np.all(np.isfinite(intensity + latitude)) or np.any(intensity <= 0.0):
        raise ValueError("C15 hourly intensity/latitude must be finite and positive")
    coriolis = np.abs(
        2.0 * EARTH_ANGULAR_VELOCITY_S * np.sin(np.deg2rad(latitude))
    )
    if np.any(coriolis <= 0.0):
        raise ValueError("C15 r0input requires non-equatorial event nodes")
    active_provider = provider or C15FixedR0WindProfileProvider(outer_radius_m)
    if not math.isclose(
        float(active_provider.outer_radius_m),
        float(outer_radius_m),
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise ValueError("injected C15 provider has the wrong event outer radius")
    rmw_m = np.asarray(
        [
            active_provider.profile_for(float(vmax), float(fcor)).radius_max_wind_m
            for vmax, fcor in zip(intensity, coriolis, strict=True)
        ],
        dtype=float,
    )
    if rmw_m.shape != intensity.shape or not np.all(np.isfinite(rmw_m)):
        raise RuntimeError("official C15 r0input returned invalid hourly RMW")
    if np.any(rmw_m <= 0.0) or np.any(rmw_m >= float(outer_radius_m)):
        raise RuntimeError("official C15 r0input returned RMW outside (0, r0)")
    hits, misses, size = active_provider.cache_info
    return rmw_m / 1000.0, {
        "method": "official C15 r0input with one fixed event outer radius",
        "c15_source_doi": C15_SOURCE_DOI,
        "c15_input_mode": "r0input",
        "intensity": "Lin v_trks circular/azimuthal maximum wind",
        "coriolis": "absolute Coriolis magnitude at each track latitude",
        "outer_radius_m": float(outer_radius_m),
        "outer_radius_fixed_for_event_lifetime": True,
        "hourly_rmw_is_solver_output": True,
        "hourly_rmw_range_km": [float(np.min(rmw_m)) / 1000.0, float(np.max(rmw_m)) / 1000.0],
        "profile_cache_hits": hits,
        "profile_cache_misses": misses,
        "profile_cache_size": size,
        "clipping_applied": False,
        "replacement_provider_applied": False,
    }


def write_prepare_dataset(
    path: Path,
    identity: dict[str, Any],
    track: dict[str, Any],
    translation_u_knots: np.ndarray,
    translation_v_knots: np.ndarray,
    ushear_ms: np.ndarray,
    vshear_ms: np.ndarray,
    t600_k: np.ndarray,
    q950: np.ndarray,
    rmw_km: np.ndarray,
    outer_radius_m: float,
) -> None:
    target = np.asarray(track["target_selector"], dtype=bool)
    native_index = np.asarray(track["target_native_index"], dtype=np.int32)
    count = native_index.size
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    with Dataset(temporary, "w", format="NETCDF4") as output:
        output.title = "Prepare-only Lin event for public C15-CLIMADA TCR"
        output.script_version = SCRIPT_VERSION
        output.event_id = identity["event_id"]
        output.event_position_in_frozen_sample = np.int64(identity["event_position"])
        output.source_track_index = np.int64(identity["source_track_index"])
        output.status = "prepared_ready"
        output.climada_track_ready = np.int64(1)
        output.blocker_codes_json = "[]"
        output.rmw_provider_decision_status = RMW_PROVIDER_DECISION_STATUS
        output.outer_radius_m = np.float64(outer_radius_m)
        output.outer_radius_fixed_for_event_lifetime = np.int64(1)
        output.c15_size_input_mode = "r0input"
        output.no_hazard_model_called = np.int64(1)
        output.createDimension("time", count)

        def variable(name: str, values: np.ndarray, units: str, **attrs: Any) -> None:
            data = np.asarray(values)
            dtype = "u1" if data.dtype.kind == "b" else "f8"
            var = output.createVariable(name, dtype, ("time",), zlib=True, complevel=4)
            var[:] = data.astype(np.uint8) if dtype == "u1" else data
            var.units = units
            for key, value in attrs.items():
                setattr(var, key, value)

        native = output.createVariable("native_index", "i4", ("time",))
        native[:] = native_index
        native.long_name = "native one-hour index in frozen 100,000-track catalogue"
        seconds = output.createVariable("time_seconds_from_seed", "f8", ("time",))
        seconds[:] = np.asarray(track["time"], dtype=float)[target]
        seconds.units = "s"
        variable("lon", track["lon_trks"][target], "degrees_east")
        variable("lat", track["lat_trks"][target], "degrees_north")
        variable(
            "circular_wind",
            track["v_trks"][target],
            "m s-1",
            source="Lin v_trks; C15 and qs900b circular-wind input",
        )
        for name in ("u250_trks", "v250_trks", "u850_trks", "v850_trks"):
            variable(name, track[name][target], "m s-1", source="frozen Lin track")
        variable(
            "translation_u",
            translation_u_knots * KNOTS_TO_MPS,
            "m s-1",
            source="Emanuel v6.4 utrans mechanically adapted to native 1 h",
        )
        variable(
            "translation_v",
            translation_v_knots * KNOTS_TO_MPS,
            "m s-1",
            source="Emanuel v6.4 utrans mechanically adapted to native 1 h",
        )
        variable(
            "ushear",
            ushear_ms,
            "m s-1",
            source="Emanuel v6.4 raingen: 5*(translation-u850)",
        )
        variable(
            "vshear",
            vshear_ms,
            "m s-1",
            source="Emanuel v6.4 raingen including fixed 1.5 m s-1 beta drift",
        )
        variable(
            "t600",
            t600_k,
            "K",
            source="CMIP6 ta600 via frozen Lin month-15 and bilinear interpolation pattern",
        )
        variable(
            "q950",
            q950,
            "kg kg-1",
            source="official Emanuel v6.4 qs900b.m; actual output pressure 950 hPa",
        )
        variable(
            "radius_max_wind",
            rmw_km,
            "km",
            source="official C15 r0input using Lin v_trks, fixed event r0, and |f|",
        )
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--fixed-r0-catalogue", type=Path, required=True)
    parser.add_argument("--fixed-r0-manifest", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--cmip6-ta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--event-position", type=int, default=0)
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="required safety flag; this program has no hazard-run mode",
    )
    args = parser.parse_args()
    if not args.prepare_only:
        parser.error("--prepare-only is required; this candidate cannot run a hazard")
    return args


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"atomic output target already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        identity, provenance = load_event_identity(
            args.sample, args.event_position, enforce_event0=False
        )
        outer_radius_m, outer_radius_metadata = load_fixed_r0_for_event(
            args.fixed_r0_catalogue,
            args.fixed_r0_manifest,
            identity,
        )
        track = load_track_window(args.tracks, identity)
        target = np.asarray(track["target_selector"], dtype=bool)
        target_lat = np.asarray(track["lat_trks"], dtype=float)[target]
        target_lon = np.asarray(track["lon_trks"], dtype=float)[target]
        target_circular = np.asarray(track["v_trks"], dtype=float)[target]
        target_u850 = np.asarray(track["u850_trks"], dtype=float)[target]
        target_v850 = np.asarray(track["v850_trks"], dtype=float)[target]

        translation_u_knots, translation_v_knots = (
            official_utrans_adapted_to_native_hour(track)
        )
        ushear_ms, vshear_ms = emanuel_v64_shear(
            translation_u_knots,
            translation_v_knots,
            target_u850,
            target_v850,
            target_lat,
            float(track["first_finite_lat_deg"]),
        )
        t600_k, t600_metadata = sample_t600_lin_public_adapter(
            args.cmip6_ta,
            target_lon,
            target_lat,
            int(identity["task_year"]),
            int(identity["seed_month"]),
        )
        q950 = emanuel_v64_qs950(t600_k, target_circular)
        provider = C15FixedR0WindProfileProvider(outer_radius_m)
        rmw_km, rmw_metadata = derive_c15_r0input_rmw(
            target_circular,
            target_lat,
            outer_radius_m,
            provider=provider,
        )

        dataset_path = temporary / "lin_event0_public_inputs_prepare_only.nc"
        write_prepare_dataset(
            dataset_path,
            identity,
            track,
            translation_u_knots,
            translation_v_knots,
            ushear_ms,
            vshear_ms,
            t600_k,
            q950,
            rmw_km,
            outer_radius_m,
        )
        runner_path = Path(__file__).resolve()
        manifest = {
            "schema_version": "1.0",
            "script_version": SCRIPT_VERSION,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "prepared_ready",
            "prepare_only": True,
            "hazard_model_called": False,
            "science_run_authorized": True,
            "rmw_provider_decision_status": RMW_PROVIDER_DECISION_STATUS,
            "event": identity,
            "frozen_inputs": {
                **provenance,
                "track": {
                    key: track[key]
                    for key in (
                        "path",
                        "bytes",
                        "sha256",
                        "source_full_track_sha256",
                        "is_hash_linked_local_snapshot",
                    )
                },
                "outer_radius": outer_radius_metadata,
                "rmw_provider": rmw_metadata,
                "cmip6_ta": t600_metadata,
            },
            "availability_and_units": {
                "translation": {
                    "available": True,
                    "source": "computed from native one-hour Lin lon/lat",
                    "output_units": "knots in public interface; also persisted as m s-1",
                },
                "v_trks": {
                    "available": True,
                    "units": "m s-1 from frozen Lin README/code contract",
                    "role": "circular wind for C15 and Emanuel qs900b",
                },
                "u850_v850": {
                    "available": True,
                    "units": "m s-1 from frozen Lin README/code contract",
                    "role": "Emanuel v6.4 synthetic shear adapter",
                },
                "u250_v250": {
                    "available": True,
                    "units": "m s-1",
                    "role": "audited only; deliberately not substituted for v6.4 shear",
                },
                "t600": {
                    "available_after_public_adapter": True,
                    "units": "K",
                    "source": "frozen CMIP6 ta at 600 hPa",
                },
                "rmw": {
                    "available_hours": int(rmw_km.size),
                    "missing_hours": 0,
                    "units": "km",
                    "range": [float(np.min(rmw_km)), float(np.max(rmw_km))],
                },
                "outer_radius": {
                    "units": "m",
                    "value": float(outer_radius_m),
                    "fixed_for_event_lifetime": True,
                    "source": "immutable event-level fixed-r0 catalogue",
                },
            },
            "public_method_contract": {
                "t600_sampling": (
                    "Lin public pattern: linear interpolation of the monthly 600-hPa "
                    "GCM temperature to datetime(task_year,seed_month,15), then "
                    "RectBivariateSpline kx=1 ky=1 along the event; one fixed monthly "
                    "field for the whole event"
                ),
                "q950": (
                    "official Emanuel v6.4 qs900b.m: T600 plus circular Lin v_trks, "
                    "five Newton updates; despite the variable name q900, pref=950 hPa"
                ),
                "translation": (
                    "official Emanuel v6.4 utrans centred difference and smfac=0.4; "
                    "dimensional factor mechanically changed from public native 2 h "
                    "to frozen Lin native 1 h; no downsampling"
                ),
                "translation_claim_boundary": "PUBLIC_CODE_DERIVED_TIMESTEP_ADAPTATION",
                "shear": (
                    "official Emanuel v6.4 raingen: ush=5*(utrans-u850), "
                    "vsh=5*(vtrans-1.5 m s-1 hemisphere_sign*cos(lat)-v850)"
                ),
                "direct_250_minus_850_used_for_tcr": False,
                "irene_ncep_environment_reused": False,
                "rmw": (
                    "official C15 r0input at every event hour, using Lin v_trks "
                    "circular intensity, one event-fixed outer radius, and |f|"
                ),
                "wind_profile_for_future_science_run": "official C15 through reviewed c15_climada_tcr.py only",
            },
            "blockers": [],
            "runner": {
                "path": str(runner_path),
                "sha256": sha256(runner_path),
            },
            "artifacts": {
                "prepare_dataset": {
                    "relative_path": dataset_path.name,
                    "bytes": dataset_path.stat().st_size,
                    "sha256": sha256(dataset_path),
                }
            },
        }
        manifest_path = temporary / "lin_event0_public_inputs_prepare_only.manifest.json"
        json_dump(manifest_path, manifest)
        os.replace(temporary, output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": manifest["status"],
                    "event_id": identity["event_id"],
                    "target_hour_count": EXPECTED_TARGET_COUNT,
                    "rmw_defined_hour_count": int(rmw_km.size),
                    "rmw_missing_hour_count": 0,
                    "rmw_range_km": [float(np.min(rmw_km)), float(np.max(rmw_km))],
                    "outer_radius_m": float(outer_radius_m),
                    "blocker_codes": [],
                    "t600_range_k": [float(np.min(t600_k)), float(np.max(t600_k))],
                    "q950_range_kgkg": [float(np.min(q950)), float(np.max(q950))],
                    "output_dir": str(output_dir),
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
