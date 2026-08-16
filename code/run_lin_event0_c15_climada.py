#!/usr/bin/env python3
"""Run or dry-run the first frozen Lin event through public C15--TCR.

The runner consumes only the prepare artifact created by
``prepare_lin_event0_c15_climada.py``.  It reuses the reviewed C15--CLIMADA
adapter, the same 0.05-degree moving 300-km union-grid construction, and the
same hourly event-total / maximum-24-hour accumulation semantics as the Irene
public reconstruction.

``--dry-run`` validates and materializes the CLIMADA track and union grid but
does not import CLIMADA or compute a hazard.  A science run requires explicit
frozen CLIMADA static fields and is intended for controlled server execution,
not this local candidate-verification step.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

import numpy as np
import xarray as xr


CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from c15_climada_tcr import (  # noqa: E402
    C15FixedR0WindProfileProvider,
    assert_environmental_pressure_schema_only,
    run_tcr_public_reconstruction,
)
from run_irene_c15_climada import (  # noqa: E402
    E_PRECIP,
    FROZEN_C_DRAG_BYTES,
    FROZEN_C_DRAG_SHA256,
    FROZEN_ELEVATION_BYTES,
    FROZEN_ELEVATION_SHA256,
    GRID_RESOLUTION_DEG,
    LOWER_TROPOSPHERE_HEIGHT_M,
    MAX_DISTANCE_EYE_KM,
    MAX_W_FOREGROUND_MPS,
    MIN_DRAG_COEFFICIENT,
    RADIAL_STEP_M,
    RHO_AIR_OVER_RHO_LIQUID,
    build_moving_union_grid,
)


SCRIPT_VERSION = "2.0.0"
EXPECTED_EVENT_ID = "stream0000-year1995-track000002"
EXPECTED_PREPARE_STATUS = "prepared_ready"
EXPECTED_TIME_COUNT = 94
EXPECTED_NATIVE_START = 56
EXPECTED_NATIVE_STOP = 149
MPS_TO_KNOTS = 3600.0 / 1852.0
KM_TO_NMI = 1.0 / 1.852

# CLIMADA requires these generic schema fields before conversion to SI.  The
# reviewed C15--TCR path statically asserts that neither is consumed.
CENTRAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA = 1000.0
ENVIRONMENTAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA = 1010.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def validate_prepare_artifact(
    prepare_nc: Path, prepare_manifest: Path
) -> tuple[xr.Dataset, dict[str, Any], dict[str, Any]]:
    prepare_nc = prepare_nc.resolve()
    prepare_manifest = prepare_manifest.resolve()
    if not prepare_nc.is_file() or not prepare_manifest.is_file():
        raise FileNotFoundError("prepare NetCDF and manifest are both required")
    manifest = json.loads(prepare_manifest.read_text())
    if (
        manifest.get("status") != EXPECTED_PREPARE_STATUS
        or manifest.get("prepare_only") is not True
        or manifest.get("hazard_model_called") is not False
        or manifest.get("blockers") != []
    ):
        raise ValueError("prepare manifest is not a ready, hazard-free artifact")
    artifact = manifest.get("artifacts", {}).get("prepare_dataset", {})
    actual_sha = sha256(prepare_nc)
    if (
        artifact.get("sha256") != actual_sha
        or int(artifact.get("bytes", -1)) != prepare_nc.stat().st_size
    ):
        raise ValueError("prepare NetCDF identity does not match its manifest")
    dataset = xr.open_dataset(prepare_nc)
    required = {
        "native_index",
        "time_seconds_from_seed",
        "lon",
        "lat",
        "circular_wind",
        "radius_max_wind",
        "q950",
        "ushear",
        "vshear",
    }
    missing = required - set(dataset.variables)
    if missing:
        dataset.close()
        raise ValueError(f"prepare artifact lacks variables: {sorted(missing)}")
    if dataset.sizes.get("time") != EXPECTED_TIME_COUNT:
        dataset.close()
        raise ValueError("unexpected first-event time count")
    native = np.asarray(dataset["native_index"].values, dtype=int)
    if not np.array_equal(
        native, np.arange(EXPECTED_NATIVE_START, EXPECTED_NATIVE_STOP + 1)
    ):
        dataset.close()
        raise ValueError("unexpected first-event native index sequence")
    finite_names = required - {"native_index"}
    for name in finite_names:
        if not np.all(np.isfinite(np.asarray(dataset[name].values, dtype=float))):
            dataset.close()
            raise ValueError(f"non-finite prepare variable: {name}")
    if np.any(np.asarray(dataset["circular_wind"].values) <= 0.0):
        dataset.close()
        raise ValueError("C15 circular winds must be positive")
    if np.any(np.asarray(dataset["radius_max_wind"].values) <= 0.0):
        dataset.close()
        raise ValueError("C15 r0input-derived RMW must be positive")
    outer_radius_m = float(dataset.attrs.get("outer_radius_m", np.nan))
    if not np.isfinite(outer_radius_m) or outer_radius_m <= 0.0:
        dataset.close()
        raise ValueError("prepare artifact lacks a positive event outer radius")
    if int(dataset.attrs.get("outer_radius_fixed_for_event_lifetime", 0)) != 1:
        dataset.close()
        raise ValueError("prepare artifact does not freeze outer radius for the event")
    event_id = str(dataset.attrs.get("event_id", ""))
    if event_id != EXPECTED_EVENT_ID:
        dataset.close()
        raise ValueError(f"unexpected event id: {event_id!r}")
    provenance = {
        "prepare_netcdf": {
            "path": str(prepare_nc),
            "bytes": prepare_nc.stat().st_size,
            "sha256": actual_sha,
        },
        "prepare_manifest": {
            "path": str(prepare_manifest),
            "bytes": prepare_manifest.stat().st_size,
            "sha256": sha256(prepare_manifest),
        },
    }
    return dataset, manifest, provenance


def build_climada_track(
    prepared: xr.Dataset,
    prepare_manifest: dict[str, Any],
    *,
    expected_event_id: str = EXPECTED_EVENT_ID,
) -> xr.Dataset:
    seconds = np.asarray(prepared["time_seconds_from_seed"].values, dtype=float)
    if not np.allclose(np.diff(seconds), 3600.0, rtol=0.0, atol=1e-6):
        raise ValueError("first event must be native one-hour data")
    event = prepare_manifest["event"]
    genesis = np.datetime64(
        str(event["threshold_genesis_datetime"]).replace(" ", "T"), "ns"
    )
    time = genesis + np.rint(seconds - seconds[0]).astype("timedelta64[s]")
    count = seconds.size
    circular_knots = (
        np.asarray(prepared["circular_wind"].values, dtype=float) * MPS_TO_KNOTS
    )
    rmw_nmi = np.asarray(prepared["radius_max_wind"].values, dtype=float) * KM_TO_NMI
    track = xr.Dataset(
        data_vars={
            "lat": ("time", np.asarray(prepared["lat"].values, dtype=float)),
            "lon": ("time", np.asarray(prepared["lon"].values, dtype=float)),
            "time_step": ("time", np.ones(count, dtype=float)),
            # C15 maximum wind: Lin circular/azimuthal v_trks, never vmax_trks.
            "max_sustained_wind": ("time", circular_knots),
            # Official C15 r0input RMW derived from fixed event r0.
            "radius_max_wind": ("time", rmw_nmi),
            "q950": ("time", np.asarray(prepared["q950"].values, dtype=float)),
            "ushear": ("time", np.asarray(prepared["ushear"].values, dtype=float)),
            "vshear": ("time", np.asarray(prepared["vshear"].values, dtype=float)),
            "central_pressure": (
                "time",
                np.full(count, CENTRAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA),
            ),
            "environmental_pressure": (
                "time",
                np.full(count, ENVIRONMENTAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA),
            ),
            "basin": (
                "time",
                np.full(count, str(event["threshold_genesis_region"]), dtype="U2"),
            ),
            "native_index": (
                "time", np.asarray(prepared["native_index"].values, dtype=np.int32)
            ),
        },
        coords={"time": time},
        attrs={
            "sid": expected_event_id,
            "name": expected_event_id,
            "orig_event_flag": 0,
            "category": 0,
            "max_sustained_wind_unit": "knots",
            "central_pressure_unit": "hPa",
            "radius_max_wind_unit": "nmile",
            "schema_pressure_placeholders_unused_by_c15_tcr": 1,
            "rmw_provider": "official C15 r0input with fixed event outer radius",
            "outer_radius_m": float(prepared.attrs["outer_radius_m"]),
            "outer_radius_fixed_for_event_lifetime": 1,
            "reconstruction_scope": "method_faithful_public_reconstruction",
        },
    )
    track["max_sustained_wind"].attrs.update(
        units="knots",
        source="Lin v_trks circular wind converted from m s-1",
        explicitly_not="Lin vmax_trks",
    )
    track["radius_max_wind"].attrs.update(
        units="nmile",
        source="official C15 r0input using Lin v_trks, fixed event r0, and |f|",
    )
    for name in ("central_pressure", "environmental_pressure"):
        track[name].attrs.update(
            units="hPa",
            schema_placeholder_only=1,
            unused_by_c15_tcr=1,
        )
    return track


def build_rainfall_output(
    track: xr.Dataset,
    grid: xr.Dataset,
    rainrate_mm_h: np.ndarray,
    adapter_metadata: dict[str, Any],
    *,
    expected_event_id: str = EXPECTED_EVENT_ID,
) -> xr.Dataset:
    rates = np.asarray(rainrate_mm_h, dtype=np.float32)
    expected = (track.sizes["time"], grid.sizes["centroid"])
    if rates.shape != expected:
        raise ValueError(f"rain-rate shape {rates.shape} != {expected}")
    if not np.all(np.isfinite(rates)) or np.any(rates < 0.0):
        raise ValueError("rainfall rates must be finite and non-negative")
    event_total = np.sum(rates, axis=0, dtype=np.float64).astype(np.float32)
    # For an event shorter than 24 native hours, the maximum available
    # accumulation is its whole available window.  No padding or rate
    # extrapolation is introduced.
    window = min(24, rates.shape[0])
    cumulative = np.concatenate(
        [
            np.zeros((1, rates.shape[1]), dtype=np.float64),
            np.cumsum(rates, axis=0, dtype=np.float64),
        ],
        axis=0,
    )
    maximum_24h = np.max(
        cumulative[window:] - cumulative[:-window], axis=0
    ).astype(np.float32)
    output = xr.Dataset(
        data_vars={
            "rainfall_rate": (("time", "centroid"), rates),
            "event_total_rainfall": ("centroid", event_total),
            "maximum_24h_rainfall": ("centroid", maximum_24h),
            "lat": ("centroid", np.asarray(grid["lat"].values, dtype=float)),
            "lon": ("centroid", np.asarray(grid["lon"].values, dtype=float)),
        },
        coords={"time": track["time"].values, "centroid": grid["centroid"].values},
        attrs={
            "event_id": expected_event_id,
            "model": "official C15 with public CLIMADA-Petals TCR skeleton",
            "scope": "method-faithful public reconstruction",
            "adapter_metadata_json": json.dumps(adapter_metadata, sort_keys=True),
            "maximum_24h_window_hours_used": window,
            "event_shorter_than_24h": int(rates.shape[0] < 24),
            "short_event_accumulation_semantics": (
                "whole available event window; no temporal padding or extrapolation"
            ),
        },
    )
    output["rainfall_rate"].attrs["units"] = "mm h-1"
    output["event_total_rainfall"].attrs.update(
        units="mm",
        accumulation=f"sum of all {rates.shape[0]} native one-hour rain-rate nodes",
    )
    output["maximum_24h_rainfall"].attrs.update(
        units="mm",
        accumulation=(
            f"maximum of all {window} consecutive one-hour nodes; 24 unless "
            "the available event window is shorter"
        ),
    )
    return output


def write_netcdf(dataset: xr.Dataset, path: Path) -> None:
    encoding: dict[str, dict[str, Any]] = {}
    for name, variable in dataset.data_vars.items():
        if variable.dtype.kind not in {"U", "S", "O"}:
            encoding[name] = {"zlib": True, "complevel": 4, "shuffle": True}
    dataset.to_netcdf(
        path, engine="netcdf4", format="NETCDF4", encoding=encoding
    )


def validate_static_field(
    path: Path, expected_bytes: int, expected_sha: str, label: str
) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"missing {label}: {path}")
    actual_sha = sha256(path)
    if path.stat().st_size != expected_bytes or actual_sha != expected_sha:
        raise ValueError(f"wrong frozen {label} identity")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": actual_sha}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-nc", type=Path, required=True)
    parser.add_argument("--prepare-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--elevation-tif", type=Path)
    parser.add_argument("--c-drag-tif", type=Path)
    args = parser.parse_args()
    if not args.dry_run and (args.elevation_tif is None or args.c_drag_tif is None):
        parser.error("science mode requires explicit --elevation-tif and --c-drag-tif")
    return args


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"atomic output target exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    prepared: xr.Dataset | None = None
    try:
        prepared, prepare_manifest, provenance = validate_prepare_artifact(
            args.prepare_nc, args.prepare_manifest
        )
        track = build_climada_track(prepared, prepare_manifest)
        # Static, hash-pinned assertion: pressure placeholders are generic
        # CLIMADA schema inputs and are not read by the selected C15--TCR path.
        assert_environmental_pressure_schema_only()
        grid = build_moving_union_grid(track)
        track_path = temporary / "lin_event0_one_hourly_climada_track.nc"
        grid_path = temporary / "lin_event0_moving_300km_grid.nc"
        write_netcdf(track, track_path)
        write_netcdf(grid, grid_path)

        artifacts: dict[str, dict[str, Any]] = {
            "climada_track": {
                "relative_path": track_path.name,
                "bytes": track_path.stat().st_size,
                "sha256": sha256(track_path),
            },
            "moving_union_grid": {
                "relative_path": grid_path.name,
                "bytes": grid_path.stat().st_size,
                "sha256": sha256(grid_path),
            },
        }
        manifest: dict[str, Any] = {
            "schema_version": "1.0",
            "script_version": SCRIPT_VERSION,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "dry_run_completed" if args.dry_run else "completed",
            "dry_run": bool(args.dry_run),
            "hazard_model_called": not args.dry_run,
            "event_id": EXPECTED_EVENT_ID,
            "runner": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "prepare_provenance": provenance,
            "track_contract": {
                "time_count": int(track.sizes["time"]),
                "native_index_range": [EXPECTED_NATIVE_START, EXPECTED_NATIVE_STOP],
                "max_sustained_wind": "Lin v_trks circular wind converted to knots",
                "radius_max_wind": "official C15 r0input output, converted km to nmile",
                "outer_radius_m": float(prepared.attrs["outer_radius_m"]),
                "outer_radius_fixed_for_event_lifetime": True,
                "q950_and_shear": "byte-linked prepare artifact",
                "central_pressure_hpa": CENTRAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA,
                "environmental_pressure_hpa": ENVIRONMENTAL_PRESSURE_SCHEMA_PLACEHOLDER_HPA,
                "pressure_role": "schema placeholder only; statically asserted unused by C15-TCR",
            },
            "grid": {
                "resolution_degrees": GRID_RESOLUTION_DEG,
                "moving_radius_km": MAX_DISTANCE_EYE_KM,
                "centroid_count": int(grid.sizes["centroid"]),
                "latitude_range": [float(grid["lat"].min()), float(grid["lat"].max())],
                "longitude_range": [float(grid["lon"].min()), float(grid["lon"].max())],
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
                "event_total": "sum rainfall_rate over all native one-hour nodes",
                "maximum_24h": "maximum sum over 24 consecutive native one-hour nodes",
            },
            "artifacts": artifacts,
        }

        if not args.dry_run:
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
            manifest["static_fields"] = {"elevation_tif": elevation, "c_drag_tif": c_drag}
            provider = C15FixedR0WindProfileProvider(
                float(prepared.attrs["outer_radius_m"])
            )
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
            rainrate = np.asarray(result["rainfall_rate_mm_h"], dtype=float)
            rainfall = build_rainfall_output(track, grid, rainrate, result["metadata"])
            rainfall_path = temporary / "lin_event0_c15_tcr_raw_rainfall.nc"
            write_netcdf(rainfall, rainfall_path)
            artifacts["raw_rainfall"] = {
                "relative_path": rainfall_path.name,
                "bytes": rainfall_path.stat().st_size,
                "sha256": sha256(rainfall_path),
            }
            manifest["adapter_metadata"] = result["metadata"]
            manifest["raw_result_summary"] = {
                "maximum_rainfall_rate_mm_h": float(rainfall["rainfall_rate"].max()),
                "maximum_event_total_mm": float(rainfall["event_total_rainfall"].max()),
                "maximum_24h_mm": float(rainfall["maximum_24h_rainfall"].max()),
            }

        manifest_path = temporary / "lin_event0_c15_tcr_run.manifest.json"
        json_dump(manifest_path, manifest)
        os.replace(temporary, output_dir)
        temporary = None
        print(
            json.dumps(
                {
                    "status": manifest["status"],
                    "hazard_model_called": manifest["hazard_model_called"],
                    "event_id": EXPECTED_EVENT_ID,
                    "time_count": int(track.sizes["time"]),
                    "grid_centroid_count": int(grid.sizes["centroid"]),
                    "output_dir": str(output_dir),
                },
                indent=2,
                sort_keys=True,
            )
        )
    finally:
        if prepared is not None:
            prepared.close()
        if temporary is not None and temporary.exists():
            shutil.rmtree(temporary)


if __name__ == "__main__":
    main()
