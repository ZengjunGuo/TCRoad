#!/usr/bin/env python3
"""Generate one Lin v1.1 GLx5000/year catalogue for an arbitrary climate window.

Historical stream-0 used a frozen 1000/2500/5000 spatial-convergence ladder
and SHA-locked 1995-2014 environment.  Future SSP windows cannot use that
runner.  This script reuses the same formal runtime, namelist physics, RNG
scheme, land-mask stage, GL downscaling entry point, and track QA, but binds
the climate identity (experiment + years) from the accepted prepared-input
and environment manifests.

It does not regenerate the environment.  It does not rerun the historical
1000/2500 ladder.  Accepted tracks per year remain a fixed quota, not a
physical annual frequency.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

import run_lin_tracks_historical_pilot as core


SCRIPT_VERSION = "1.0.0"
MODEL = "MPI-ESM1-2-LR"
MEMBER = "r1i1p1f1"
BASIN = "GL"
MAIN_QUOTA = 5000
DEFAULT_YEAR_WORKERS = 8
DEFAULT_TRACK_TIMEOUT_SECONDS = 36 * 60 * 60
WIND_TIME_AXIS = {
    "statistic_period": "complete_calendar_month",
    "coordinate": "month_day_15_00_00",
    "downstream_use": "linear interpolation at track seed month day 15",
}


def bind_window(experiment: str, start_year: int, end_year: int) -> None:
    """Point the historical generation helpers at one climate window."""
    if end_year < start_year:
        raise ValueError(f"invalid year window {start_year}-{end_year}")
    core.MODEL = MODEL
    core.EXPERIMENT = experiment
    core.MEMBER = MEMBER
    core.PREFIX = f"{MODEL}_{experiment}_{MEMBER}"
    core.START_YEAR = start_year
    core.END_YEAR = end_year
    core.YEARS = tuple(range(start_year, end_year + 1))
    core.BASIN = BASIN


def git_state(project: Path) -> tuple[str, str]:
    commit = subprocess.check_output(
        ["git", "-C", str(project), "rev-parse", "HEAD"], text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(project), "status", "--short"], text=True
    ).strip()
    return commit, status


def derived_roots(project: Path, experiment: str, start_year: int, end_year: int) -> tuple[Path, Path]:
    window = f"{start_year}-{end_year}"
    label = f"{MODEL}_{experiment}_{MEMBER}_{window}_GLx{MAIN_QUOTA}peryear_stream{core.STREAM_ID}"
    work_root = project / "scratch/lin_track_window_work" / label
    publish_root = (
        project
        / "data/lin/tracks"
        / MODEL
        / experiment
        / MEMBER
        / window
        / f"GLx{MAIN_QUOTA}peryear_stream{core.STREAM_ID}"
    )
    return work_root, publish_root


def require_global_coverage(qa: dict, label: str) -> None:
    for key in (
        "all_seven_basins_in_every_year",
        "all_twelve_months_in_every_year",
        "all_seven_basins_in_main_catalogue",
        "all_twelve_months_in_main_catalogue",
    ):
        if qa.get(key) is not True:
            raise ValueError(f"{label} fails global coverage gate: {key}")


def validate_environment(manifest_path: Path, project: Path) -> tuple[dict, Path, Path, dict]:
    environment = core.read_pass_manifest(manifest_path, "climate-window environment")
    if environment.get("request") != core.expected_request():
        raise ValueError(f"environment request mismatch: {environment.get('request')}")
    if environment.get("configuration", {}).get("wind_time_axis") != WIND_TIME_AXIS:
        raise ValueError("environment lacks the accepted mid-month wind-time contract")

    outputs = environment.get("outputs", {})
    if set(outputs) != {"wind", "thermo"}:
        raise ValueError("environment manifest does not define exactly wind and thermo")
    wind_path = core.verify_record(outputs["wind"], project, "accepted wind environment")
    thermo_path = core.verify_record(outputs["thermo"], project, "accepted thermo environment")
    expected_time = core.expected_monthly_time()
    expected_months = len(core.YEARS) * 12

    with xr.open_dataset(wind_path) as wind:
        if dict(wind.sizes) != {"lon": 192, "lat": 96, "time": expected_months}:
            raise ValueError(f"wrong accepted-wind dimensions: {dict(wind.sizes)}")
        if not np.array_equal(wind.time.values, expected_time):
            raise ValueError("accepted wind is not the complete mid-month product")
        qa = outputs["wind"].get("qa", {})
        if qa.get("nonpositive_eigenvalue_count") != 0:
            raise ValueError("accepted wind manifest does not pass covariance PD")
        if qa.get("covariance_domain") != "entire 96x192 global grid":
            raise ValueError("accepted wind covariance gate was not global")
        if qa.get("covariance_matrix_count") != expected_months * 96 * 192:
            raise ValueError("accepted wind covariance count mismatch")
    with xr.open_dataset(thermo_path) as thermo:
        if dict(thermo.sizes) != {"time": expected_months, "lat": 96, "lon": 192}:
            raise ValueError(f"wrong accepted-thermo dimensions: {dict(thermo.sizes)}")
        if not np.array_equal(thermo.time.values, expected_time):
            raise ValueError("accepted thermo is not the complete mid-month product")
        for name in thermo.data_vars:
            if not np.isfinite(thermo[name].values).all():
                raise ValueError(f"accepted thermo contains non-finite values: {name}")

    upstream_runtime = environment.get("upstream", {}).get("runtime_manifest", {})
    core.verify_record(upstream_runtime, project, "environment-generation runtime manifest")
    return environment, wind_path, thermo_path, upstream_runtime


def preflight(
    project: Path,
    prepared_root: Path,
    environment_manifest: Path,
    formal_runtime: Path,
    python: Path,
) -> dict:
    if core.STREAM_ID != 0:
        raise ValueError("climate-window catalogues use RNG stream 0")
    if core.BASE_SEED != 20260809:
        raise ValueError("climate-window catalogues use the frozen base seed 20260809")
    environment, wind_path, thermo_path, environment_runtime = validate_environment(
        environment_manifest, project
    )
    prepared_record = core.validate_prepared(prepared_root, environment, project)
    runtime_record, runtime_manifest = core.validate_runtime(formal_runtime, project)
    runtime_binding = core.validate_environment_runtime_binding(
        environment_runtime, runtime_record, runtime_manifest
    )
    commit, status = git_state(project)
    return {
        "git_commit": commit,
        "git_status": status,
        "environment": environment,
        "wind_path": wind_path,
        "thermo_path": thermo_path,
        "environment_record": core.file_record(environment_manifest, project),
        "prepared_record": prepared_record,
        "runtime_record": runtime_record,
        "runtime_binding": runtime_binding,
        "python": python,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--start-year", required=True, type=int)
    parser.add_argument("--end-year", required=True, type=int)
    parser.add_argument("--prepared-root", required=True, type=Path)
    parser.add_argument("--environment-manifest", required=True, type=Path)
    parser.add_argument("--formal-runtime", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--year-workers", type=int, default=DEFAULT_YEAR_WORKERS)
    parser.add_argument(
        "--track-timeout-seconds",
        type=int,
        default=DEFAULT_TRACK_TIMEOUT_SECONDS,
    )
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    experiment = args.experiment.strip()
    if experiment not in {"ssp126", "ssp245", "ssp370", "ssp585", "historical"}:
        raise ValueError(f"unsupported experiment: {experiment}")
    bind_window(experiment, args.start_year, args.end_year)

    project = args.project.resolve()
    prepared_root = core.inside_project(args.prepared_root, project)
    environment_manifest = core.inside_project(args.environment_manifest, project)
    formal_runtime = core.inside_project(args.formal_runtime, project)
    python = args.python.absolute()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise FileNotFoundError(f"Python launcher is not executable: {python}")
    n_years = len(core.YEARS)
    if args.year_workers < 1 or args.year_workers > n_years:
        raise ValueError(f"year-workers must be in 1..{n_years}")
    if args.track_timeout_seconds <= 0:
        raise ValueError("track-timeout-seconds must be positive")

    work_root, publish_root = derived_roots(
        project, experiment, args.start_year, args.end_year
    )
    work_root = core.inside_project(work_root, project, allow_existing=False)
    publish_root = core.inside_project(publish_root, project, allow_existing=False)
    core.ensure_disjoint([work_root, publish_root])

    upstream = preflight(
        project,
        prepared_root,
        environment_manifest,
        formal_runtime,
        python,
    )
    if args.preflight_only:
        print(
            json.dumps(
                {
                    "status": "pass",
                    "mode": "read_only_preflight",
                    "request": {
                        **core.expected_request(),
                        "basin_mode": BASIN,
                        "accepted_tracks_per_year": MAIN_QUOTA,
                        "rng_base_seed": core.BASE_SEED,
                        "rng_stream_id": core.STREAM_ID,
                    },
                    "git_commit": upstream["git_commit"],
                    "year_workers": args.year_workers,
                    "environment_regeneration": False,
                    "work_root": str(work_root),
                    "publish_root": str(publish_root),
                },
                sort_keys=True,
            )
        )
        return

    accepted_environment_hashes = {
        "wind": core.sha256(upstream["wind_path"]),
        "thermo": core.sha256(upstream["thermo_path"]),
    }
    work_root.parent.mkdir(parents=True, exist_ok=True)
    publish_root.parent.mkdir(parents=True, exist_ok=True)
    if work_root.parent.stat().st_dev != publish_root.parent.stat().st_dev:
        raise ValueError("work and publication roots must share one filesystem")
    work_root.mkdir()

    try:
        main_run = core.run_once(
            work_root / f"main_{MAIN_QUOTA}",
            formal_runtime,
            upstream["wind_path"],
            upstream["thermo_path"],
            python,
            args.track_timeout_seconds,
            MAIN_QUOTA,
            args.year_workers,
        )
        require_global_coverage(
            main_run["qa"],
            f"{experiment} {args.start_year}-{args.end_year} GLx{MAIN_QUOTA}/year",
        )
        if accepted_environment_hashes != {
            "wind": core.sha256(upstream["wind_path"]),
            "thermo": core.sha256(upstream["thermo_path"]),
        }:
            raise ValueError("accepted environment changed during track generation")
        core.validate_runtime(formal_runtime, project)

        output_staging = work_root / "output"
        record_staging = work_root / "record"
        output_staging.mkdir()
        record_staging.mkdir()
        final_output = publish_root / "output"
        final_record = publish_root / "record"
        final_track_staging = output_staging / main_run["track"].name
        shutil.move(str(main_run["track"]), final_track_staging)
        shutil.copy2(main_run["archived_namelist"], output_staging / "namelist.py")
        shutil.copy2(main_run["land_log"], record_staging / "main_5000_land.log")
        shutil.copy2(main_run["track_log"], record_staging / "main_5000_tracks.log")

        manifest = {
            "schema_version": 1,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "status": "pass",
            "scope": (
                f"{args.start_year}-{args.end_year} GL {experiment} stream-0 "
                f"GLx{MAIN_QUOTA}/year catalogue; accepted quotas are not "
                "physical annual TC counts"
            ),
            "script": {
                "path": core.relative(Path(__file__).resolve(), project),
                "version": SCRIPT_VERSION,
                "sha256": core.sha256(Path(__file__).resolve()),
                "generation_helper": core.file_record(Path(core.__file__), project),
            },
            "git_commit": upstream["git_commit"],
            "git_status_at_start": upstream["git_status"],
            "request": {
                **core.expected_request(),
                "basin_mode": BASIN,
                "accepted_tracks_per_year": MAIN_QUOTA,
                "rng_base_seed": core.BASE_SEED,
                "rng_stream_id": core.STREAM_ID,
            },
            "upstream": {
                "prepared": upstream["prepared_record"],
                "environment_manifest": upstream["environment_record"],
                "environment_outputs": {
                    "wind": core.file_record(upstream["wind_path"], project),
                    "thermo": core.file_record(upstream["thermo_path"], project),
                },
                "formal_runtime": upstream["runtime_record"],
                "environment_runtime_binding": upstream["runtime_binding"],
            },
            "parallelism": {
                "logical_task": "one GL calendar year",
                "task_count": n_years,
                "year_workers": args.year_workers,
                "maximum_useful_year_workers": n_years,
                "dask_result_order": (
                    f"input year order {args.start_year}..{args.end_year}"
                ),
                "rng_isolation": (
                    "SeedSequence(base,stream,GL,calendar_year) resets RNG at each year task"
                ),
            },
            "configuration": {
                "effective_namelist": main_run["effective_namelist"],
                "thread_environment": {
                    "OMP_NUM_THREADS": 1,
                    "OPENBLAS_NUM_THREADS": 1,
                    "MKL_NUM_THREADS": 1,
                    "NUMEXPR_NUM_THREADS": 1,
                    "PYTHONHASHSEED": 0,
                },
                "environment_reuse": {
                    "wind_and_thermo_symlinked_into_each_isolated_run": True,
                    "environment_generation_called": False,
                    "source_hashes_verified_before_and_after": True,
                },
            },
            "sampling_semantics": {
                "accepted_tracks_per_year_is_fixed_quota": True,
                "accepted_track_count_is_not_physical_annual_frequency": True,
                "basin_or_month_quota_imposed": False,
                "same_rng_stream_as_historical_baseline": True,
                "climate_identity_enters_only_through_environment": True,
                "qualified_seed_definition": main_run["qa"]["lin_seed_count_definition"],
                "raw_proposal_count_available": False,
                "frequency_calibration_required_downstream": True,
                "authoritative_native_track_support": main_run["qa"]["canonical_support"],
                "tc_event_window": (
                    "inclusive first through last timestep with vmax >= 18 m/s"
                ),
            },
            "stages": {"main_5000": main_run["stages"]},
            "generated_land_masks": {
                **main_run["land_qa"],
                "retained_after_publication": False,
            },
            "outputs": {
                "track": {
                    **core.file_record(
                        final_track_staging,
                        project,
                        published_path=final_output / final_track_staging.name,
                    ),
                    "qa": main_run["qa"],
                },
                "namelist": core.file_record(
                    output_staging / "namelist.py",
                    project,
                    published_path=final_output / "namelist.py",
                ),
            },
            "publishing": {
                "atomic_unit": core.relative(publish_root, project),
                "method": "single same-filesystem rename of completed work root",
                "transient_run_directories_removed_before_publish": True,
                "failure_retains_single_work_root": True,
            },
            "logs": {
                "main_5000_land": core.file_record(
                    record_staging / "main_5000_land.log",
                    project,
                    published_path=final_record / "main_5000_land.log",
                ),
                "main_5000_tracks": core.file_record(
                    record_staging / "main_5000_tracks.log",
                    project,
                    published_path=final_record / "main_5000_tracks.log",
                ),
            },
            "software": {
                "python": platform.python_version(),
                "python_executable": str(python),
                "numpy": np.__version__,
            },
        }
        manifest_path = record_staging / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        shutil.rmtree(main_run["run_root"])
        if {path.name for path in work_root.iterdir()} != {"output", "record"}:
            raise ValueError("completed work root contains unexpected publication entries")
        success_payload = json.dumps(
            {
                "status": "pass",
                "publish_root": str(publish_root),
                "track_sha256": manifest["outputs"]["track"]["sha256"],
                "accepted_tracks": MAIN_QUOTA * n_years,
                "environment_regenerated": False,
            },
            sort_keys=True,
        )
        work_root.rename(publish_root)
    except Exception:
        if work_root.exists():
            print(f"failed isolated work retained at {work_root}", file=sys.stderr)
        raise
    print(success_payload)


if __name__ == "__main__":
    main()
