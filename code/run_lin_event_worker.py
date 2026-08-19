#!/usr/bin/env python3
"""Run one frozen-sample event through the reviewed four-stage public chain.

The worker is deliberately orchestration glue.  It calls the existing prepare,
C15--TCR rainfall, C15--Lin--Chavas wind, and descriptive road-overlay
functions; no scientific equation is duplicated here.  One outer radius is
read from the immutable event-level catalogue and the same official C15
``r0input`` provider is reused by prepare, rainfall, and wind.  A finite-domain
or official-solver failure is published as ``METHOD_DOMAIN_PENDING`` for the
batch audit queue, without clipping, replacement, exclusion, resampling, or
weight modification.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

from netCDF4 import Dataset
import numpy as np
import xarray as xr

from c15_climada_tcr import (
    C15FixedR0WindProfileProvider,
    C15ProfileDomainError,
    assert_environmental_pressure_schema_only,
    run_tcr_public_reconstruction,
)
from overlay_lin_event0_hazards_on_roads import (
    build_road_overlap,
    load_hazard_contracts,
    summarize_by_road_class,
    write_summary_csv,
)
from prepare_lin_event0_c15_climada import (
    derive_c15_r0input_rmw,
    emanuel_v64_qs950,
    emanuel_v64_shear,
    load_event_identity,
    load_fixed_r0_for_event,
    load_track_window,
    official_utrans_adapted_to_native_hour,
    sample_t600_lin_public_adapter,
    write_prepare_dataset,
)
from run_irene_c15_climada import (
    E_PRECIP,
    FROZEN_C_DRAG_BYTES,
    FROZEN_C_DRAG_SHA256,
    FROZEN_ELEVATION_BYTES,
    FROZEN_ELEVATION_SHA256,
    LOWER_TROPOSPHERE_HEIGHT_M,
    MAX_DISTANCE_EYE_KM,
    MAX_W_FOREGROUND_MPS,
    MIN_DRAG_COEFFICIENT,
    RADIAL_STEP_M,
    RHO_AIR_OVER_RHO_LIQUID,
    build_moving_union_grid,
)
from run_lin_event0_c15_climada import (
    build_climada_track,
    build_rainfall_output,
    sha256,
    validate_static_field,
    write_netcdf,
)
from run_lin_event0_c15_windfield import compute_wind_field


SCRIPT_VERSION = "2.0.0"

# Official C15 r0-input solver failures. IndexError is the empty V_ER11
# profile from ER11_radprof_raw; RuntimeError is the wrapped solver path.
OFFICIAL_R0INPUT_SOLVER_FAILURE_TYPES = (RuntimeError, IndexError)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def publish_method_domain_audit(
    args: argparse.Namespace,
    identity: dict[str, Any],
    *,
    audit_code: str,
    detail: str,
    outer_radius_m: float | None,
) -> None:
    """Publish the only durable worker product for an unresolved method case."""

    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_json(
        args.output_dir / "scientific_audit.json",
        {
            "schema_version": "1.0",
            "status": "scientific_audit_required",
            "method_status": "METHOD_DOMAIN_PENDING",
            "audit_code": audit_code,
            "audit_detail": detail,
            "event_id": identity["event_id"],
            "event_position": int(identity["event_position"]),
            "outer_radius_m": outer_radius_m,
            "resolution_applied": False,
            "event_excluded": False,
            "event_weight_modified": False,
            "clipping_applied": False,
            "tail_or_taper_applied": False,
            "rejection_or_resampling_applied": False,
            "replacement_provider_applied": False,
            "event_weight_climate_fixed_effect_ht_analysis_yr": float(
                identity["event_weight_climate_fixed_effect_ht_analysis_yr"]
            ),
        },
    )


def artifact(path: Path) -> dict[str, Any]:
    return {
        "relative_path": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def build_compact_hazard_footprint(
    rainfall: xr.Dataset,
    wind: xr.Dataset,
    event: dict[str, Any],
    *,
    right_censored: bool,
) -> xr.Dataset:
    """Drop hourly arrays while preserving the reversible spatial footprint."""

    if not np.array_equal(rainfall["centroid"].values, wind["centroid"].values):
        raise ValueError("rainfall and wind centroid identifiers differ")
    for name in ("lat", "lon"):
        if not np.allclose(
            rainfall[name].values, wind[name].values, rtol=0.0, atol=1e-10
        ):
            raise ValueError(f"rainfall and wind {name} coordinates differ")
    footprint = xr.Dataset(
        data_vars={
            "lat": ("centroid", np.asarray(rainfall["lat"].values, dtype=np.float32)),
            "lon": ("centroid", np.asarray(rainfall["lon"].values, dtype=np.float32)),
            "event_total_rainfall": (
                "centroid",
                np.asarray(rainfall["event_total_rainfall"].values, dtype=np.float32),
            ),
            "maximum_24h_rainfall": (
                "centroid",
                np.asarray(rainfall["maximum_24h_rainfall"].values, dtype=np.float32),
            ),
            "event_maximum_near_surface_wind_speed": (
                "centroid",
                np.asarray(
                    wind["event_maximum_near_surface_wind_speed"].values,
                    dtype=np.float32,
                ),
            ),
        },
        coords={"centroid": np.asarray(rainfall["centroid"].values)},
        attrs={
            "event_id": event["event_id"],
            "event_position": int(event["event_position"]),
            "event_weight_climate_fixed_effect_ht_analysis_yr": float(
                event["event_weight_climate_fixed_effect_ht_analysis_yr"]
            ),
            "right_censored_at_15_day_limit": int(right_censored),
            "cumulative_rainfall_interpretation": (
                "lower bound over available 15-day-limited window"
                if right_censored
                else "complete threshold-window accumulation"
            ),
            "hourly_fields_included": 0,
            "hazard_thresholds_applied": 0,
            "damage_or_loss_model_applied": 0,
            "wind_averaging_period": wind.attrs["wind_averaging_period"],
            "outer_radius_m": float(event["outer_radius_m"]),
            "outer_radius_fixed_for_event_lifetime": 1,
            "fixed_r0_catalogue_sha256": str(
                event["fixed_r0_catalogue_sha256"]
            ),
            "fixed_r0_distribution_contract_sha256": str(
                event["fixed_r0_distribution_contract_sha256"]
            ),
        },
    )
    footprint["lat"].attrs["units"] = "degrees_north"
    footprint["lon"].attrs["units"] = "degrees_east"
    footprint["event_total_rainfall"].attrs["units"] = "mm"
    footprint["maximum_24h_rainfall"].attrs["units"] = "mm"
    footprint["event_maximum_near_surface_wind_speed"].attrs.update(
        units="m s-1",
        averaging_period=wind.attrs["wind_averaging_period"],
    )
    return footprint


def _prepare_event(
    args: argparse.Namespace, temporary: Path
) -> tuple[
    dict[str, Any],
    xr.Dataset,
    dict[str, Any],
    Path,
    Path,
    bool,
    C15FixedR0WindProfileProvider,
]:
    sample_sha = sha256(args.sample)
    identity, provenance = load_event_identity(
        args.sample,
        args.event_position,
        enforce_event0=False,
        expected_sample_sha256=sample_sha,
    )
    outer_radius_m, outer_radius_metadata = load_fixed_r0_for_event(
        args.fixed_r0_catalogue,
        args.fixed_r0_manifest,
        identity,
        expected_sample_sha256=sample_sha,
    )
    catalogue_attrs = outer_radius_metadata["catalogue_attributes"]
    identity.update(
        outer_radius_m=float(outer_radius_m),
        fixed_r0_catalogue_sha256=outer_radius_metadata["catalogue"]["sha256"],
        fixed_r0_distribution_contract_sha256=str(
            catalogue_attrs["distribution_contract_sha256"]
        ),
        fixed_r0_event_outer_radius_binding_sha256=str(
            catalogue_attrs["event_outer_radius_binding_sha256"]
        ),
    )
    required_outer_radius_m = MAX_DISTANCE_EYE_KM * 1000.0 + RADIAL_STEP_M
    if outer_radius_m < required_outer_radius_m:
        publish_method_domain_audit(
            args,
            identity,
            audit_code="C15_FIXED_R0_BELOW_302KM_NUMERICAL_DOMAIN",
            detail=(
                f"fixed r0={outer_radius_m:.9g} m is smaller than the frozen "
                f"{MAX_DISTANCE_EYE_KM:g}-km active analysis radius plus the "
                f"{RADIAL_STEP_M / 1000.0:g}-km radial derivative support"
            ),
            outer_radius_m=outer_radius_m,
        )
        raise SystemExit(20)
    track_data = load_track_window(args.tracks, identity)
    target = np.asarray(track_data["target_selector"], dtype=bool)
    target_lat = np.asarray(track_data["lat_trks"], dtype=float)[target]
    target_lon = np.asarray(track_data["lon_trks"], dtype=float)[target]
    target_circular = np.asarray(track_data["v_trks"], dtype=float)[target]
    provider = C15FixedR0WindProfileProvider(outer_radius_m)
    try:
        rmw_km, rmw_metadata = derive_c15_r0input_rmw(
            target_circular,
            target_lat,
            outer_radius_m,
            provider=provider,
        )
    except OFFICIAL_R0INPUT_SOLVER_FAILURE_TYPES as error:
        # Official ER11/C15 adapter can raise IndexError when V_ER11 is empty.
        # That is the same frozen case as RuntimeError: METHOD_DOMAIN_PENDING,
        # no resampling, no clip, no replacement solver.
        publish_method_domain_audit(
            args,
            identity,
            audit_code="C15_R0INPUT_SOLVER_FAILURE",
            detail=f"{type(error).__name__}: {error}",
            outer_radius_m=outer_radius_m,
        )
        raise SystemExit(20) from error

    translation_u_knots, translation_v_knots = (
        official_utrans_adapted_to_native_hour(track_data)
    )
    target_u850 = np.asarray(track_data["u850_trks"], dtype=float)[target]
    target_v850 = np.asarray(track_data["v850_trks"], dtype=float)[target]
    ushear_ms, vshear_ms = emanuel_v64_shear(
        translation_u_knots,
        translation_v_knots,
        target_u850,
        target_v850,
        target_lat,
        float(track_data["first_finite_lat_deg"]),
    )
    t600_k, t600_metadata = sample_t600_lin_public_adapter(
        args.cmip6_ta,
        target_lon,
        target_lat,
        int(identity["task_year"]),
        int(identity["seed_month"]),
    )
    q950 = emanuel_v64_qs950(t600_k, target_circular)

    prepare_dir = temporary / "prepare"
    prepare_dir.mkdir()
    prepare_nc = prepare_dir / "lin_event_public_inputs_prepare_only.nc"
    write_prepare_dataset(
        prepare_nc,
        identity,
        track_data,
        translation_u_knots,
        translation_v_knots,
        ushear_ms,
        vshear_ms,
        t600_k,
        q950,
        rmw_km,
        outer_radius_m,
    )
    with Dataset(args.tracks) as source_tracks:
        track_native_last_index = len(source_tracks.dimensions["time"]) - 1
    right_censored = int(identity["native_stop"]) == track_native_last_index
    prepare_manifest = {
        "schema_version": "1.0",
        "script_version": SCRIPT_VERSION,
        "generated_at_utc": utc_now(),
        "status": "prepared_ready",
        "prepare_only": True,
        "hazard_model_called": False,
        "blockers": [],
        "event": {
            **identity,
            "right_censored_at_15_day_limit": right_censored,
            "cumulative_hazard_interpretation": (
                "lower_bound_over_available_window"
                if right_censored
                else "complete_threshold_window_accumulation"
            ),
        },
        "frozen_inputs": {
            **provenance,
            "track": {
                key: track_data[key]
                for key in (
                    "path",
                    "bytes",
                    "sha256",
                    "source_full_track_sha256",
                    "is_hash_linked_local_snapshot",
                )
            },
            "cmip6_ta": t600_metadata,
            "outer_radius": outer_radius_metadata,
            "rmw_provider": rmw_metadata,
        },
        "artifacts": {"prepare_dataset": artifact(prepare_nc)},
    }
    prepare_manifest_path = prepare_dir / "lin_event_prepare_only.manifest.json"
    write_json(prepare_manifest_path, prepare_manifest)
    prepared = xr.open_dataset(prepare_nc)
    return (
        identity,
        prepared,
        prepare_manifest,
        prepare_nc,
        prepare_manifest_path,
        right_censored,
        provider,
    )


def run_event(args: argparse.Namespace, temporary: Path) -> None:
    (
        identity,
        prepared,
        prepare_manifest,
        prepare_nc,
        prepare_manifest_path,
        right_censored,
        provider,
    ) = _prepare_event(args, temporary)
    event_id = str(identity["event_id"])
    track: xr.Dataset | None = None
    grid: xr.Dataset | None = None
    rainfall: xr.Dataset | None = None
    wind: xr.Dataset | None = None
    overlap: xr.Dataset | None = None
    footprint: xr.Dataset | None = None
    try:
        track = build_climada_track(
            prepared, prepare_manifest, expected_event_id=event_id
        )
        assert_environmental_pressure_schema_only()
        grid = build_moving_union_grid(track)
        elevation = validate_static_field(
            args.elevation_tif,
            FROZEN_ELEVATION_BYTES,
            FROZEN_ELEVATION_SHA256,
            "topography",
        )
        c_drag = validate_static_field(
            args.c_drag_tif,
            FROZEN_C_DRAG_BYTES,
            FROZEN_C_DRAG_SHA256,
            "drag coefficient",
        )

        rain_dir = temporary / "rain"
        rain_dir.mkdir()
        track_path = rain_dir / "lin_event_one_hourly_climada_track.nc"
        grid_path = rain_dir / "lin_event_moving_300km_grid.nc"
        write_netcdf(track, track_path)
        write_netcdf(grid, grid_path)
        result = run_tcr_public_reconstruction(
            track=track,
            centroid_lat=np.asarray(grid["lat"].values, dtype=float),
            centroid_lon=np.asarray(grid["lon"].values, dtype=float),
            elevation_tif=Path(elevation["path"]),
            c_drag_tif=Path(c_drag["path"]),
            e_precip=E_PRECIP,
            lower_troposphere_height_m=LOWER_TROPOSPHERE_HEIGHT_M,
            rho_air_over_rho_liquid=RHO_AIR_OVER_RHO_LIQUID,
            max_w_foreground=MAX_W_FOREGROUND_MPS,
            res_radial_m=RADIAL_STEP_M,
            min_c_drag=MIN_DRAG_COEFFICIENT,
            max_dist_eye_km=MAX_DISTANCE_EYE_KM,
            provider=provider,
        )
        rainfall = build_rainfall_output(
            track,
            grid,
            np.asarray(result["rainfall_rate_mm_h"], dtype=float),
            result["metadata"],
            expected_event_id=event_id,
        )
        rainfall.attrs["right_censored_at_15_day_limit"] = int(right_censored)
        rainfall.attrs["cumulative_rainfall_interpretation"] = (
            "lower bound over available 15-day-limited window"
            if right_censored
            else "complete threshold-window accumulation"
        )
        rainfall_path = rain_dir / "lin_event_c15_tcr_raw_rainfall.nc"
        write_netcdf(rainfall, rainfall_path)
        write_json(
            rain_dir / "lin_event_c15_tcr_run.manifest.json",
            {
                "schema_version": "1.0",
                "script_version": SCRIPT_VERSION,
                "generated_at_utc": utc_now(),
                "status": "completed",
                "event_id": event_id,
                "prepare_provenance": {
                    "prepare_netcdf": artifact(prepare_nc),
                    "prepare_manifest": artifact(prepare_manifest_path),
                },
                "static_fields": {
                    "elevation_tif": elevation,
                    "c_drag_tif": c_drag,
                },
                "adapter_metadata": result["metadata"],
                "raw_result_summary": {
                    "time_count": int(rainfall.sizes["time"]),
                    "centroid_count": int(rainfall.sizes["centroid"]),
                    "maximum_rainfall_rate_mm_h": float(
                        rainfall["rainfall_rate"].max()
                    ),
                    "maximum_event_total_mm": float(
                        rainfall["event_total_rainfall"].max()
                    ),
                    "maximum_24h_mm": float(
                        rainfall["maximum_24h_rainfall"].max()
                    ),
                    "right_censored_at_15_day_limit": right_censored,
                },
                "artifacts": {
                    "raw_rainfall": artifact(rainfall_path),
                    "climada_track": artifact(track_path),
                    "moving_union_grid": artifact(grid_path),
                },
            },
        )

        wind_dir = temporary / "wind"
        wind_dir.mkdir()
        wind, wind_metadata = compute_wind_field(
            prepared,
            track,
            grid,
            provider=provider,
            expected_event_id=event_id,
        )
        wind_path = wind_dir / "lin_event_c15_lin_chavas_windfield.nc"
        write_netcdf(wind, wind_path)
        write_json(
            wind_dir / "lin_event_c15_windfield.manifest.json",
            {
                "schema_version": "1.0",
                "script_version": SCRIPT_VERSION,
                "generated_at_utc": utc_now(),
                "status": "completed",
                "event_id": event_id,
                "provider": wind_metadata,
                "result_summary": {
                    "time_count": int(wind.sizes["time"]),
                    "centroid_count": int(wind.sizes["centroid"]),
                    "maximum_model_native_near_surface_wind_m_s": float(
                        wind["event_maximum_near_surface_wind_speed"].max()
                    ),
                },
                "artifacts": {"windfield": artifact(wind_path)},
            },
        )

        road_dir = temporary / "road"
        road_dir.mkdir()
        hazard = load_hazard_contracts(
            rainfall_path, wind_path, expected_event_id=event_id
        )
        overlap, overlap_metadata = build_road_overlap(args.roads_nc, hazard)
        overlap.attrs.update(
            outer_radius_m=float(identity["outer_radius_m"]),
            outer_radius_fixed_for_event_lifetime=1,
            fixed_r0_catalogue_sha256=str(
                identity["fixed_r0_catalogue_sha256"]
            ),
            fixed_r0_distribution_contract_sha256=str(
                identity["fixed_r0_distribution_contract_sha256"]
            ),
        )
        summary = summarize_by_road_class(overlap)
        overlap_path = road_dir / "lin_event_road_grid_joint_exposure.nc"
        summary_path = road_dir / "lin_event_road_class_joint_exposure_summary.json"
        summary_csv_path = road_dir / "lin_event_road_class_joint_exposure_summary.csv"
        write_netcdf(overlap, overlap_path)
        write_json(summary_path, summary)
        write_summary_csv(summary_csv_path, summary)
        write_json(
            road_dir / "lin_event_road_overlap.manifest.json",
            {
                "schema_version": "1.0",
                "script_version": SCRIPT_VERSION,
                "generated_at_utc": utc_now(),
                "status": "completed",
                "event_id": event_id,
                "result_summary": overlap_metadata,
                "artifacts": {
                    "road_grid_joint_exposure": artifact(overlap_path),
                    "road_class_summary_json": artifact(summary_path),
                    "road_class_summary_csv": artifact(summary_csv_path),
                },
            },
        )

        compact_dir = temporary / "compact"
        compact_dir.mkdir()
        footprint = build_compact_hazard_footprint(
            rainfall, wind, identity, right_censored=right_censored
        )
        footprint_path = compact_dir / "lin_event_compact_hazard_footprint.nc"
        write_netcdf(footprint, footprint_path)
        write_json(
            compact_dir / "lin_event_compact_hazard_footprint.manifest.json",
            {
                "schema_version": "1.0",
                "script_version": SCRIPT_VERSION,
                "generated_at_utc": utc_now(),
                "status": "completed",
                "event_id": event_id,
                "right_censored_at_15_day_limit": right_censored,
                "artifact": artifact(footprint_path),
            },
        )
    finally:
        for dataset in (footprint, overlap, wind, rainfall, grid, track, prepared):
            if dataset is not None:
                dataset.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--fixed-r0-catalogue", type=Path, required=True)
    parser.add_argument("--fixed-r0-manifest", type=Path, required=True)
    parser.add_argument("--event-position", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--cmip6-ta", type=Path, required=True)
    parser.add_argument("--roads-nc", type=Path, required=True)
    parser.add_argument("--elevation-tif", type=Path, required=True)
    parser.add_argument("--c-drag-tif", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.sample = args.sample.resolve()
    args.fixed_r0_catalogue = args.fixed_r0_catalogue.resolve()
    args.fixed_r0_manifest = args.fixed_r0_manifest.resolve()
    args.tracks = args.tracks.resolve()
    args.cmip6_ta = args.cmip6_ta.resolve()
    args.roads_nc = args.roads_nc.resolve()
    args.elevation_tif = args.elevation_tif.resolve()
    args.c_drag_tif = args.c_drag_tif.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.output_dir.exists():
        raise FileExistsError(f"atomic output target exists: {args.output_dir}")
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent
        )
    )
    try:
        try:
            run_event(args, temporary)
        except C15ProfileDomainError as error:
            sample_sha = sha256(args.sample)
            identity, _ = load_event_identity(
                args.sample,
                args.event_position,
                enforce_event0=False,
                expected_sample_sha256=sample_sha,
            )
            outer_radius_m, _ = load_fixed_r0_for_event(
                args.fixed_r0_catalogue,
                args.fixed_r0_manifest,
                identity,
                expected_sample_sha256=sample_sha,
            )
            publish_method_domain_audit(
                args,
                identity,
                audit_code="C15_ACTIVE_QUERY_OUTSIDE_FIXED_R0",
                detail=str(error),
                outer_radius_m=outer_radius_m,
            )
            raise SystemExit(20) from error
        os.replace(temporary, args.output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": "completed",
                    "event_position": args.event_position,
                    "output_dir": str(args.output_dir),
                },
                sort_keys=True,
            )
        )
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
