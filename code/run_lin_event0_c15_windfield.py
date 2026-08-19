#!/usr/bin/env python3
"""Generate the first frozen Lin event's two-dimensional C15 wind field.

This runner consumes the reviewed prepare artifact and the already-materialized
moving 0.05-degree / 300-km union grid.  At each native one-hour node it:

1. evaluates the official C15 ``r0input`` near-surface azimuthal profile with
   the event-fixed outer radius by linear interpolation at great-circle radii;
2. maps that profile to the cyclonic tangential vector using one
   CLIMADA 6.1 majority-hemisphere sign for the whole track
   (counter-clockwise if northern nodes win or tie, clockwise if
   southern nodes strictly outnumber northern nodes); and
3. adds the spatially uniform Lin--Chavas (2012) surface-background vector,
   ``0.55 * translation``, rotated 20 degrees cyclonically.

No surface-wind reduction factor, inflow angle, radial decay of the background
wind, empirical taper, or wind-averaging-period conversion is applied.  In
particular, the frozen Lin ``v_trks`` averaging period is not documented, so
this output is labelled model-native near-surface wind and is not called a
10-minute sustained wind field.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Protocol

import numpy as np
import xarray as xr

from c15_climada_tcr import (
    C15FixedR0WindProfileProvider,
    C15ProfileDomainError,
    C15_PYTHON3_ADAPTER_SHA256,
    C15_R0INPUT_DEFAULTS,
    C15_SOURCE_DOI,
)
from run_irene_c15_climada import (
    EARTH_RADIUS_KM,
    GRID_RESOLUTION_DEG,
    MAX_DISTANCE_EYE_KM,
    build_moving_union_grid,
)
from run_lin_event0_c15_climada import (
    EXPECTED_EVENT_ID,
    build_climada_track,
    json_dump,
    sha256,
    validate_prepare_artifact,
    write_netcdf,
)


SCRIPT_VERSION = "2.1.0"
CLIMADA_MAJORITY_HEMISPHERE_RULE = (
    "climada_6.1.0_tctrack_to_si_majority_node_count"
)
LIN_CHAVAS_2012_DOI = "10.1029/2011JD017126"
EARTH_ANGULAR_VELOCITY_S = 7.2921159e-5
BACKGROUND_REDUCTION_FACTOR = 0.55
BACKGROUND_CCW_ROTATION_DEG = 20.0


class _Profile(Protocol):
    radius_m: np.ndarray
    gradient_wind_ms: np.ndarray
    outer_radius_m: float


class _ProfileProvider(Protocol):
    outer_radius_m: float

    def profile_for(
        self, vmax_ms: float, coriolis_s: float
    ) -> _Profile: ...


def spherical_distance_and_outward_bearing(
    center_lat_deg: float,
    center_lon_deg: float,
    point_lat_deg: np.ndarray,
    point_lon_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return distance (m) and point-local outward bearing from the center.

    The bearing is the forward direction of the center-to-point great circle
    *at the point*, clockwise from local north.  It is therefore the reverse
    point-to-center initial bearing plus pi, not the initial bearing measured
    in the center's tangent plane.  Longitudinal differences are wrapped to
    ``[-pi, pi]``.  Bearing at a point exactly coincident with the center is set
    to zero; the C15 vortex speed at that radius is zero, so its direction has
    no physical effect.
    """

    latitude = np.asarray(point_lat_deg, dtype=float)
    longitude = np.asarray(point_lon_deg, dtype=float)
    if latitude.ndim != 1 or longitude.shape != latitude.shape:
        raise ValueError(
            "point latitude and longitude must share one-dimensional shape"
        )
    inputs = np.concatenate(
        [
            np.asarray([center_lat_deg, center_lon_deg], dtype=float),
            latitude,
            longitude,
        ]
    )
    if not np.all(np.isfinite(inputs)):
        raise ValueError("spherical coordinates must be finite")

    center_lat = math.radians(float(center_lat_deg))
    point_lat = np.deg2rad(latitude)
    delta_lat = point_lat - center_lat
    delta_lon = (np.deg2rad(longitude - float(center_lon_deg)) + np.pi) % (
        2.0 * np.pi
    ) - np.pi
    haversine = (
        np.sin(0.5 * delta_lat) ** 2
        + math.cos(center_lat) * np.cos(point_lat) * np.sin(0.5 * delta_lon) ** 2
    )
    central_angle = 2.0 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
    distance_m = EARTH_RADIUS_KM * 1000.0 * central_angle

    # Calculate the initial point-to-center direction in the point's local
    # tangent plane, then reverse it to obtain the locally outward continuation
    # of the center-to-point great circle.
    reverse_delta_lon = -delta_lon
    toward_center_eastward = np.sin(reverse_delta_lon) * math.cos(center_lat)
    toward_center_northward = np.cos(point_lat) * math.sin(center_lat) - np.sin(
        point_lat
    ) * math.cos(center_lat) * np.cos(reverse_delta_lon)
    bearing_rad = np.arctan2(
        -toward_center_eastward,
        -toward_center_northward,
    )
    coincident = central_angle <= np.finfo(float).eps
    bearing_rad[coincident] = 0.0
    return distance_m, bearing_rad


def climada_majority_hemisphere_sign(
    lat_deg: np.ndarray,
) -> tuple[float, int, int]:
    """Return CLIMADA 6.1 ``tctrack_to_si`` majority-hemisphere ``latsign``.

    Frozen official rule (CLIMADA core ``v6.1.0``)::

        hemisphere = "N"
        if count(lat < 0) > count(lat > 0):
            hemisphere = "S"
        latsign = +1.0 if hemisphere == "N" else -1.0

    Equator nodes (``lat == 0``) count in neither side. A tie, including an
    all-equator track, stays Northern Hemisphere. The whole track uses one
    constant sign; Coriolis remains per-timestep ``sin(lat)``. Tracks are not
    split, dropped, or reweighted.
    """

    lat = np.asarray(lat_deg, dtype=float)
    if lat.ndim != 1 or lat.size == 0:
        raise ValueError("track latitudes must be a non-empty 1-D array")
    if not np.all(np.isfinite(lat)):
        raise ValueError("track latitudes must be finite")
    northern = int(np.count_nonzero(lat > 0.0))
    southern = int(np.count_nonzero(lat < 0.0))
    hemisphere_sign = -1.0 if southern > northern else 1.0
    return hemisphere_sign, northern, southern


def lin_chavas_background_wind(
    translation_u_ms: np.ndarray,
    translation_v_ms: np.ndarray,
    *,
    hemisphere_sign: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the published 0.55 / 20-degree cyclonic transform."""

    translation_u = np.asarray(translation_u_ms, dtype=float)
    translation_v = np.asarray(translation_v_ms, dtype=float)
    if translation_u.shape != translation_v.shape:
        raise ValueError("translation components must share shape")
    if not np.all(np.isfinite(translation_u + translation_v)):
        raise ValueError("translation components must be finite")
    if hemisphere_sign not in (-1.0, 1.0):
        raise ValueError("hemisphere_sign must be -1 or +1")
    angle = math.radians(BACKGROUND_CCW_ROTATION_DEG * hemisphere_sign)
    cosine = math.cos(angle)
    sine = math.sin(angle)
    background_u = BACKGROUND_REDUCTION_FACTOR * (
        cosine * translation_u - sine * translation_v
    )
    background_v = BACKGROUND_REDUCTION_FACTOR * (
        sine * translation_u + cosine * translation_v
    )
    return background_u, background_v


def validate_moving_grid(
    prepared: xr.Dataset,
    prepare_manifest: dict[str, Any],
    moving_grid_path: Path,
) -> tuple[xr.Dataset, dict[str, Any]]:
    """Verify that an input grid is exactly the frozen moving-grid product."""

    moving_grid_path = moving_grid_path.resolve()
    if not moving_grid_path.is_file():
        raise FileNotFoundError(moving_grid_path)
    observed = xr.open_dataset(moving_grid_path)
    try:
        required = {
            "centroid",
            "lat",
            "lon",
            "global_latitude_index",
            "global_longitude_index",
        }
        missing = required - set(observed.variables)
        if missing:
            raise ValueError(f"moving grid lacks variables: {sorted(missing)}")
        if not math.isclose(
            float(observed.attrs.get("grid_resolution_degrees", np.nan)),
            GRID_RESOLUTION_DEG,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("moving grid resolution is not the frozen 0.05 degrees")
        if not math.isclose(
            float(observed.attrs.get("maximum_distance_to_hourly_eye_km", np.nan)),
            MAX_DISTANCE_EYE_KM,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("moving grid radius is not the frozen 300 km")

        track = build_climada_track(prepared, prepare_manifest)
        expected = build_moving_union_grid(track)
        if observed.sizes.get("centroid") != expected.sizes["centroid"]:
            raise ValueError("moving grid centroid count differs from the frozen grid")
        for name in required:
            if not np.array_equal(observed[name].values, expected[name].values):
                raise ValueError(
                    f"moving grid differs from frozen construction: {name}"
                )
        observed.load()
    except Exception:
        observed.close()
        raise

    return observed, {
        "path": str(moving_grid_path),
        "bytes": moving_grid_path.stat().st_size,
        "sha256": sha256(moving_grid_path),
        "exact_frozen_construction_recomputed": True,
    }


def compute_wind_field(
    prepared: xr.Dataset,
    track: xr.Dataset,
    grid: xr.Dataset,
    *,
    provider: _ProfileProvider | None = None,
    expected_event_id: str = EXPECTED_EVENT_ID,
) -> tuple[xr.Dataset, dict[str, Any]]:
    """Compute hourly model-native near-surface vectors and event maximum."""

    required = {
        "lat",
        "lon",
        "circular_wind",
        "radius_max_wind",
        "translation_u",
        "translation_v",
    }
    missing = required - set(prepared.variables)
    if missing:
        raise ValueError(f"prepare artifact lacks wind-field inputs: {sorted(missing)}")
    for name in required:
        values = np.asarray(prepared[name].values, dtype=float)
        if values.shape != (prepared.sizes["time"],) or not np.all(np.isfinite(values)):
            raise ValueError(f"invalid prepared wind-field variable: {name}")

    center_lat = np.asarray(prepared["lat"].values, dtype=float)
    center_lon = np.asarray(prepared["lon"].values, dtype=float)
    hemisphere_sign, northern_nodes, southern_nodes = (
        climada_majority_hemisphere_sign(center_lat)
    )
    circular_wind = np.asarray(prepared["circular_wind"].values, dtype=float)
    prepared_rmw_m = (
        np.asarray(prepared["radius_max_wind"].values, dtype=float) * 1000.0
    )
    if np.any(circular_wind <= 0.0) or np.any(prepared_rmw_m <= 0.0):
        raise ValueError("C15 wind speed and RMW must be positive")
    outer_radius_m = float(prepared.attrs.get("outer_radius_m", np.nan))
    if not np.isfinite(outer_radius_m) or outer_radius_m <= 0.0:
        raise ValueError("prepare artifact lacks the fixed event outer radius")

    point_lat = np.asarray(grid["lat"].values, dtype=float)
    point_lon = np.asarray(grid["lon"].values, dtype=float)
    if point_lat.ndim != 1 or point_lon.shape != point_lat.shape:
        raise ValueError("moving-grid coordinates must share centroid shape")

    translation_u = np.asarray(prepared["translation_u"].values, dtype=float)
    translation_v = np.asarray(prepared["translation_v"].values, dtype=float)
    background_u, background_v = lin_chavas_background_wind(
        translation_u, translation_v, hemisphere_sign=hemisphere_sign
    )
    active_provider = provider or C15FixedR0WindProfileProvider(outer_radius_m)
    if not math.isclose(
        float(active_provider.outer_radius_m),
        outer_radius_m,
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise ValueError("wind provider outer radius differs from prepare artifact")

    time_count = prepared.sizes["time"]
    centroid_count = point_lat.size
    wind_u = np.full((time_count, centroid_count), np.nan, dtype=np.float32)
    wind_v = np.full_like(wind_u, np.nan)
    wind_speed = np.full_like(wind_u, np.nan)
    active_domain = np.zeros((time_count, centroid_count), dtype=np.uint8)
    event_maximum = np.full(centroid_count, np.nan, dtype=np.float32)
    outer_radius_km = np.empty(time_count, dtype=float)
    active_counts = np.empty(time_count, dtype=np.int64)

    for time_index in range(time_count):
        distance_m, bearing_rad = spherical_distance_and_outward_bearing(
            center_lat[time_index], center_lon[time_index], point_lat, point_lon
        )
        active = distance_m <= MAX_DISTANCE_EYE_KM * 1000.0
        if not np.any(active):
            raise ValueError(f"moving domain is empty at time index {time_index}")
        coriolis_s = (
            2.0
            * EARTH_ANGULAR_VELOCITY_S
            * math.sin(math.radians(center_lat[time_index]))
        )
        profile = active_provider.profile_for(circular_wind[time_index], coriolis_s)
        if not math.isclose(
            float(profile.radius_max_wind_m),
            float(prepared_rmw_m[time_index]),
            rel_tol=0.0,
            abs_tol=1e-6,
        ):
            raise RuntimeError(
                "wind provider RMW differs from the prepared official r0input output"
            )
        query_radius_m = distance_m[active]
        if np.any(query_radius_m > float(profile.outer_radius_m)):
            raise C15ProfileDomainError(
                "active 300-km wind-field query exceeds official C15 r0 at "
                f"time_index={time_index}; no clipping or tail is permitted"
            )
        vortex_speed = np.interp(
            query_radius_m,
            np.asarray(profile.radius_m, dtype=float),
            np.asarray(profile.gradient_wind_ms, dtype=float),
        )
        # Point-local outward bearing is clockwise from north.  Cyclonic
        # tangential flow is counter-clockwise in NH and clockwise in SH.
        vortex_u = (
            -hemisphere_sign * vortex_speed * np.cos(bearing_rad[active])
        )
        vortex_v = hemisphere_sign * vortex_speed * np.sin(bearing_rad[active])
        row_u = vortex_u + background_u[time_index]
        row_v = vortex_v + background_v[time_index]
        row_speed = np.hypot(row_u, row_v)
        if not np.all(np.isfinite(row_speed)):
            raise ValueError(f"non-finite active wind at time index {time_index}")

        wind_u[time_index, active] = row_u.astype(np.float32)
        wind_v[time_index, active] = row_v.astype(np.float32)
        wind_speed[time_index, active] = row_speed.astype(np.float32)
        active_domain[time_index, active] = 1
        if time_index == 0:
            event_maximum[active] = row_speed.astype(np.float32)
        else:
            event_maximum[active] = np.fmax(
                event_maximum[active], row_speed.astype(np.float32)
            )
        outer_radius_km[time_index] = float(profile.outer_radius_m) / 1000.0
        active_counts[time_index] = int(np.count_nonzero(active))

    if np.any(~np.isfinite(event_maximum)):
        raise ValueError(
            "moving grid contains a centroid outside every hourly 300-km domain"
        )
    output = xr.Dataset(
        data_vars={
            "near_surface_wind_u": (("time", "centroid"), wind_u),
            "near_surface_wind_v": (("time", "centroid"), wind_v),
            "near_surface_wind_speed": (("time", "centroid"), wind_speed),
            "active_300km_domain": (("time", "centroid"), active_domain),
            "event_maximum_near_surface_wind_speed": (
                "centroid",
                event_maximum,
            ),
            "lat": ("centroid", point_lat),
            "lon": ("centroid", point_lon),
        },
        coords={
            "time": track["time"].values,
            "centroid": np.asarray(grid["centroid"].values),
        },
        attrs={
            "event_id": expected_event_id,
            "model": "official C15 plus Lin-Chavas-2012 surface background wind",
            "wind_level": "model-native near-surface",
            "wind_averaging_period": "unspecified in frozen Lin v_trks source",
            "ten_minute_sustained_wind_claim": 0,
            "wind_averaging_period_conversion_applied": 0,
            "lin_chavas_2012_doi": LIN_CHAVAS_2012_DOI,
            "c15_code_archive_doi": C15_SOURCE_DOI,
            "background_wind_reduction_factor": BACKGROUND_REDUCTION_FACTOR,
            "background_wind_rotation_degrees": (
                BACKGROUND_CCW_ROTATION_DEG * hemisphere_sign
            ),
            "background_wind_rotation": (
                "cyclonic: counter-clockwise in NH, clockwise in SH"
            ),
            "background_wind_spatial_form": "uniform within each active hourly domain",
            "evaluation_domain": (
                "joint C15-TCR rainfall-analysis domain: <=300 km from each hourly center"
            ),
            "evaluation_domain_claim_boundary": (
                "analysis support shared with rainfall; not a physical C15 wind cutoff, "
                "not the C15 r0 footprint, and not the Nature 250-km POI event criterion"
            ),
            "surface_wind_reduction_factor_applied_to_c15": 0,
            "inflow_angle_applied_to_c15": 0,
            "radial_background_decay_or_taper_applied": 0,
            "c15_interpolation": "linear on official radial profile",
            "c15_size_input_mode": "r0input with one fixed event outer radius",
            "outer_radius_m": outer_radius_m,
            "inactive_hourly_cells": "NaN; not evaluated outside moving 300-km domain",
        },
    )
    for name in (
        "near_surface_wind_u",
        "near_surface_wind_v",
        "near_surface_wind_speed",
        "event_maximum_near_surface_wind_speed",
    ):
        output[name].attrs["units"] = "m s-1"
        output[name].attrs[
            "averaging_period"
        ] = "unspecified in frozen Lin v_trks source; no conversion applied"
    output["near_surface_wind_u"].attrs["positive"] = "eastward"
    output["near_surface_wind_v"].attrs["positive"] = "northward"
    output["active_300km_domain"].attrs.update(
        units="1",
        definition="1 where spherical distance to that hourly center is <=300 km",
    )
    output["lat"].attrs["units"] = "degrees_north"
    output["lon"].attrs["units"] = "degrees_east"

    cache_info = getattr(active_provider, "cache_info", None)
    provider_metadata: dict[str, Any] = {
        "c15_source_doi": C15_SOURCE_DOI,
        "c15_python3_adapter_sha256": C15_PYTHON3_ADAPTER_SHA256,
        "c15_input_mode": "r0input",
        "c15_outer_radius_m": outer_radius_m,
        "c15_r0input_parameters": dict(
            getattr(active_provider, "parameters", C15_R0INPUT_DEFAULTS)
        ),
        "linear_profile_interpolation": True,
        "c15_outer_radius_km_range": [
            float(np.min(outer_radius_km)),
            float(np.max(outer_radius_km)),
        ],
        "active_centroid_count_range": [
            int(np.min(active_counts)),
            int(np.max(active_counts)),
        ],
        "background_u_m_s_range": [
            float(np.min(background_u)),
            float(np.max(background_u)),
        ],
        "background_v_m_s_range": [
            float(np.min(background_v)),
            float(np.max(background_v)),
        ],
        "hemisphere_sign": hemisphere_sign,
        "hemisphere_rule": CLIMADA_MAJORITY_HEMISPHERE_RULE,
        "hemisphere_northern_node_count": northern_nodes,
        "hemisphere_southern_node_count": southern_nodes,
        "cyclonic_background_rotation_degrees": (
            BACKGROUND_CCW_ROTATION_DEG * hemisphere_sign
        ),
    }
    if cache_info is not None:
        hits, misses, size = cache_info
        provider_metadata["profile_cache"] = {
            "hits": int(hits),
            "misses": int(misses),
            "size": int(size),
        }
    return output, provider_metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-nc", type=Path, required=True)
    parser.add_argument("--prepare-manifest", type=Path, required=True)
    parser.add_argument("--moving-grid", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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
    prepared: xr.Dataset | None = None
    grid: xr.Dataset | None = None
    try:
        prepared, prepare_manifest, prepare_provenance = validate_prepare_artifact(
            args.prepare_nc, args.prepare_manifest
        )
        track = build_climada_track(prepared, prepare_manifest)
        grid, grid_provenance = validate_moving_grid(
            prepared, prepare_manifest, args.moving_grid
        )
        provider = C15FixedR0WindProfileProvider(
            float(prepared.attrs["outer_radius_m"])
        )
        wind, provider_metadata = compute_wind_field(
            prepared, track, grid, provider=provider
        )
        wind_path = temporary / "lin_event0_c15_lin_chavas_windfield.nc"
        write_netcdf(wind, wind_path)
        runner_path = Path(__file__).resolve()
        manifest = {
            "schema_version": "1.0",
            "script_version": SCRIPT_VERSION,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "completed",
            "event_id": EXPECTED_EVENT_ID,
            "runner": {
                "path": str(runner_path),
                "bytes": runner_path.stat().st_size,
                "sha256": sha256(runner_path),
            },
            "inputs": {
                **prepare_provenance,
                "moving_union_grid": grid_provenance,
            },
            "method_contract": {
                "axisymmetric_vortex": (
                    "official C15 r0input with one event-fixed outer radius, "
                    "evaluated by linear interpolation; cyclonic sign is the "
                    "CLIMADA 6.1 majority-hemisphere latsign"
                ),
                "c15_amplitude": "frozen Lin v_trks in m s-1",
                "c15_radius_max_wind": (
                    "official r0input hourly output persisted in prepare artifact"
                ),
                "c15_outer_radius_m": float(prepared.attrs["outer_radius_m"]),
                "surface_background": (
                    "0.55 times prepared translation vector, rotated 20 degrees "
                    "cyclonically by the CLIMADA 6.1 majority-hemisphere latsign, "
                    "spatially uniform"
                ),
                "lin_chavas_2012_doi": LIN_CHAVAS_2012_DOI,
                "moving_domain": "spherical distance <=300 km at each native hour",
                "moving_domain_claim_boundary": (
                    "joint rainfall-analysis support only; not a physical C15 wind "
                    "cutoff, not the full C15 r0 wind footprint, and not Nature's "
                    "250-km POI hazard-event selection distance"
                ),
                "excluded": [
                    "surface wind reduction factor",
                    "inflow angle",
                    "radial decay of translation background",
                    "radial taper or empirical wind tail",
                    "wind averaging-period conversion",
                ],
                "wind_averaging_period": (
                    "unspecified in frozen Lin v_trks source; output is not labelled "
                    "10-minute sustained wind"
                ),
            },
            "provider": provider_metadata,
            "result_summary": {
                "time_count": int(wind.sizes["time"]),
                "centroid_count": int(wind.sizes["centroid"]),
                "maximum_model_native_near_surface_wind_m_s": float(
                    wind["event_maximum_near_surface_wind_speed"].max()
                ),
            },
            "artifacts": {
                "windfield": {
                    "relative_path": wind_path.name,
                    "bytes": wind_path.stat().st_size,
                    "sha256": sha256(wind_path),
                }
            },
        }
        manifest_path = temporary / "lin_event0_c15_windfield.manifest.json"
        json_dump(manifest_path, manifest)
        os.replace(temporary, output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": "completed",
                    "event_id": EXPECTED_EVENT_ID,
                    **manifest["result_summary"],
                    "output_dir": str(output_dir),
                },
                indent=2,
                sort_keys=True,
            )
        )
    finally:
        if prepared is not None:
            prepared.close()
        if grid is not None:
            grid.close()
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
