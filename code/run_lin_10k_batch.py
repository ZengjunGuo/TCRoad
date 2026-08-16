#!/usr/bin/env python3
"""Restartable orchestration for the frozen 10,000-event hazard sample.

This file is a scheduling and compaction layer only.  It does not implement,
modify, or post-process any wind, rainfall, or road-exposure equation.  A
separate event worker must call the reviewed prepare, C15--TCR rainfall,
C15--Lin--Chavas wind, and road-overlay modules and leave their existing
manifests below the supplied event output directory.

Every event in the selected shard is attempted, including events whose lysis
node is the final native node of the 15-day Lin track.  Such events are marked
``right_censored_at_15_day_limit`` and their cumulative quantities are labelled
as lower bounds over the available window; they are never screened out.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import traceback
from typing import Any, Iterator, Sequence
import uuid

from netCDF4 import Dataset, chartostring
import numpy as np


SCRIPT_VERSION = "2.0.0"
EXPECTED_SAMPLE_COUNT = 10_000
ROAD_CLASS_NAMES = ("highways", "primary", "secondary", "tertiary", "local")

PREPARE_MANIFEST_GLOB = "**/*prepare_only.manifest.json"
RAIN_MANIFEST_GLOB = "**/*tcr_run.manifest.json"
WIND_MANIFEST_GLOB = "**/*windfield.manifest.json"
ROAD_MANIFEST_GLOB = "**/*road_overlap.manifest.json"
ROAD_SUMMARY_GLOB = "**/*road_class_joint_exposure_summary.json"
ROAD_OVERLAP_NC_GLOB = "**/*road_grid_joint_exposure.nc"
COMPACT_FOOTPRINT_NC_GLOB = "**/*compact_hazard_footprint.nc"
SCIENTIFIC_AUDIT_REPORT_GLOB = "**/scientific_audit.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Publish one complete JSON document with a same-directory rename."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
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


def _nc_strings(variable: Any) -> np.ndarray:
    values = variable[:]
    if values.dtype.kind in {"U", "O"}:
        return np.asarray(values, dtype=str)
    if values.dtype.kind == "S" and values.ndim == 2:
        return np.asarray(chartostring(values), dtype=str)
    return np.asarray(values).astype(str)


def load_sample_events(
    sample_path: Path,
    *,
    track_native_last_index: int,
    require_10k: bool = True,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Read orchestration metadata without changing the frozen sample."""

    sample_path = sample_path.resolve()
    if not sample_path.is_file():
        raise FileNotFoundError(sample_path)
    required = {
        "event_id",
        "source_track_index",
        "source_catalogue_event_position",
        "task_year",
        "threshold_genesis_native_index",
        "threshold_lysis_native_index",
        "event_weight_climate_fixed_effect_ht_analysis_yr",
    }
    with Dataset(sample_path) as sample:
        missing = required - set(sample.variables)
        if missing:
            raise ValueError(f"frozen sample lacks variables: {sorted(missing)}")
        count = len(sample.dimensions["event"])
        if require_10k and count != EXPECTED_SAMPLE_COUNT:
            raise ValueError(f"expected {EXPECTED_SAMPLE_COUNT} events, found {count}")
        event_ids = _nc_strings(sample.variables["event_id"])
        arrays = {
            name: np.asarray(sample.variables[name][:])
            for name in required - {"event_id"}
        }

    events: list[dict[str, Any]] = []
    for position in range(count):
        start = int(arrays["threshold_genesis_native_index"][position])
        stop = int(arrays["threshold_lysis_native_index"][position])
        if start < 0 or stop < start or stop > track_native_last_index:
            raise ValueError(
                f"invalid native window at event position {position}: {start}..{stop}"
            )
        right_censored = stop == track_native_last_index
        events.append(
            {
                "event_position": position,
                "event_id": str(event_ids[position]),
                "source_track_index": int(arrays["source_track_index"][position]),
                "source_catalogue_event_position": int(
                    arrays["source_catalogue_event_position"][position]
                ),
                "task_year": int(arrays["task_year"][position]),
                "native_start": start,
                "native_stop": stop,
                "available_hour_count": stop - start + 1,
                "event_weight_climate_fixed_effect_ht_analysis_yr": float(
                    arrays["event_weight_climate_fixed_effect_ht_analysis_yr"][position]
                ),
                "right_censored_at_15_day_limit": right_censored,
                "cumulative_hazard_interpretation": (
                    "lower_bound_over_available_window"
                    if right_censored
                    else "complete_threshold_window_accumulation"
                ),
            }
        )
    return events, {
        "path": str(sample_path),
        "bytes": sample_path.stat().st_size,
        "sha256": sha256(sample_path),
        "event_count": count,
        "track_native_last_index": int(track_native_last_index),
        "right_censored_event_count": sum(
            int(event["right_censored_at_15_day_limit"]) for event in events
        ),
    }


def validate_fixed_r0_catalogue(
    catalogue_path: Path,
    manifest_path: Path,
    events: Sequence[dict[str, Any]],
    sample_record: dict[str, Any],
) -> dict[str, Any]:
    """Verify the immutable one-to-one event/r0 binding before any shard runs."""

    catalogue_path = catalogue_path.resolve()
    manifest_path = manifest_path.resolve()
    if not catalogue_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("fixed-r0 catalogue NetCDF and manifest are required")
    manifest = _load_json(manifest_path)
    if manifest.get("status") != "FROZEN_IMMUTABLE":
        raise ValueError("fixed-r0 catalogue manifest is not immutable and complete")
    artifact = manifest.get("artifacts", {}).get("fixed_r0_catalogue_netcdf", {})
    catalogue_sha = sha256(catalogue_path)
    if (
        artifact.get("sha256") != catalogue_sha
        or int(artifact.get("bytes", -1)) != catalogue_path.stat().st_size
    ):
        raise ValueError("fixed-r0 catalogue identity differs from its manifest")
    if manifest.get("source_sample", {}).get("sha256") != sample_record["sha256"]:
        raise ValueError("fixed-r0 catalogue manifest is linked to another sample")

    with Dataset(catalogue_path) as catalogue:
        required = {"event_position", "event_id", "outer_radius_m"}
        missing = required - set(catalogue.variables)
        if missing:
            raise ValueError(f"fixed-r0 catalogue lacks variables: {sorted(missing)}")
        if str(getattr(catalogue, "status", "")) != "FROZEN_IMMUTABLE":
            raise ValueError("fixed-r0 catalogue NetCDF is not immutable")
        if str(getattr(catalogue, "source_sample_sha256", "")) != sample_record["sha256"]:
            raise ValueError("fixed-r0 catalogue NetCDF is linked to another sample")
        positions = np.asarray(catalogue.variables["event_position"][:], dtype=np.int64)
        event_ids = _nc_strings(catalogue.variables["event_id"])
        outer_radius_m = np.asarray(
            catalogue.variables["outer_radius_m"][:], dtype=np.float64
        )
        attrs = {
            name: str(getattr(catalogue, name))
            for name in (
                "source_event_id_order_sha256",
                "distribution_contract_sha256",
                "outer_radius_m_sequence_sha256",
                "event_outer_radius_binding_sha256",
            )
            if hasattr(catalogue, name)
        }
    count = len(events)
    expected_ids = np.asarray([str(event["event_id"]) for event in events])
    if positions.shape != (count,) or not np.array_equal(
        positions, np.arange(count, dtype=np.int64)
    ):
        raise ValueError("fixed-r0 catalogue event positions are not exact sample order")
    if event_ids.shape != (count,) or not np.array_equal(event_ids, expected_ids):
        raise ValueError("fixed-r0 catalogue event IDs are not exact sample order")
    if outer_radius_m.shape != (count,) or not np.all(np.isfinite(outer_radius_m)):
        raise ValueError("fixed-r0 catalogue radius vector is invalid")
    if np.any(outer_radius_m <= 0.0):
        raise ValueError("fixed-r0 catalogue contains a non-positive radius")
    required_outer_radius_m = 302_000.0
    below_domain = np.flatnonzero(outer_radius_m < required_outer_radius_m)
    # Events below 302 km stay in the sample with unchanged weights and are
    # published as METHOD_DOMAIN_PENDING by the worker.  They are not
    # dropped, clipped, or redrawn here.

    exact = manifest.get("exact_catalogue_sequence", {})
    method = manifest.get("method", {})
    expected_attrs = {
        "source_event_id_order_sha256": manifest.get("source_sample", {}).get(
            "event_id_order_sha256"
        ),
        "distribution_contract_sha256": method.get(
            "distribution_contract_sha256"
        ),
        "outer_radius_m_sequence_sha256": exact.get(
            "outer_radius_m_sequence_sha256"
        ),
        "event_outer_radius_binding_sha256": exact.get(
            "event_outer_radius_binding_sha256"
        ),
    }
    if any(not value for value in expected_attrs.values()) or attrs != expected_attrs:
        raise ValueError("fixed-r0 catalogue contract hashes disagree with manifest")
    for event, radius in zip(events, outer_radius_m, strict=True):
        event["outer_radius_m"] = float(radius)
        event["fixed_r0_catalogue_sha256"] = catalogue_sha
        event["fixed_r0_distribution_contract_sha256"] = attrs[
            "distribution_contract_sha256"
        ]
        event["fixed_r0_event_outer_radius_binding_sha256"] = attrs[
            "event_outer_radius_binding_sha256"
        ]
    return {
        "path": str(catalogue_path),
        "bytes": catalogue_path.stat().st_size,
        "sha256": catalogue_sha,
        "manifest_path": str(manifest_path),
        "manifest_bytes": manifest_path.stat().st_size,
        "manifest_sha256": sha256(manifest_path),
        "event_count": count,
        "minimum_outer_radius_m": float(np.min(outer_radius_m)),
        "maximum_outer_radius_m": float(np.max(outer_radius_m)),
        "minimum_required_outer_radius_m": required_outer_radius_m,
        "count_below_required_outer_radius": int(below_domain.size),
        **attrs,
    }


def select_shard(
    events: Sequence[dict[str, Any]], shard_index: int, shard_count: int
) -> list[dict[str, Any]]:
    """Use sample position only, so selection is deterministic and order-free."""

    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("require 0 <= shard_index < shard_count")
    return [
        event
        for event in events
        if int(event["event_position"]) % shard_count == shard_index
    ]


def representative_qa_positions(
    events: Sequence[dict[str, Any]],
    *,
    require_frozen_10k_identities: bool = True,
) -> dict[str, int]:
    """Return coverage cases; these are QA selections, not science gates."""

    if not events:
        raise ValueError("representative QA selection requires at least one event")
    by_position = {int(event["event_position"]): event for event in events}
    required_positions = {
        "minimum_fixed_r0": min(
            events, key=lambda event: float(event["outer_radius_m"])
        )["event_position"],
        "maximum_fixed_r0": max(
            events, key=lambda event: float(event["outer_radius_m"])
        )["event_position"],
        # These positions were frozen from the exact 10k/full-track audit.
        "lowest_absolute_track_latitude": 1732,
        "southern_hemisphere": 1,
        "periodic_longitude_seam_crossing": 2649,
        "right_censored_at_15_day_limit": 6,
        "shorter_than_24h": 103,
        "former_nonpositive_size_predictor_case": 953,
    }
    result = {name: int(position) for name, position in required_positions.items()}
    missing = sorted(set(result.values()) - set(by_position))
    if missing:
        raise ValueError(
            "representative QA positions are absent from the frozen sample: "
            f"{missing}"
        )
    predicates = {
        "southern_hemisphere": (
            str(by_position[result["southern_hemisphere"]]["event_id"])
            == "stream0000-year1995-track000007"
        ),
        "periodic_longitude_seam_crossing": (
            str(by_position[result["periodic_longitude_seam_crossing"]]["event_id"])
            == "stream0000-year2000-track026613"
        ),
        "right_censored_at_15_day_limit": bool(
            by_position[result["right_censored_at_15_day_limit"]][
                "right_censored_at_15_day_limit"
            ]
        ),
        "shorter_than_24h": (
            int(by_position[result["shorter_than_24h"]]["available_hour_count"])
            < 24
        ),
        "former_nonpositive_size_predictor_case": (
            str(
                by_position[result["former_nonpositive_size_predictor_case"]][
                    "event_id"
                ]
            )
            == "stream0000-year1996-track009518"
        ),
    }
    if not require_frozen_10k_identities:
        return {
            "minimum_fixed_r0": int(required_positions["minimum_fixed_r0"]),
            "maximum_fixed_r0": int(required_positions["maximum_fixed_r0"]),
        }
    failed = sorted(name for name, valid in predicates.items() if not valid)
    if failed:
        raise ValueError(f"representative QA identity contract changed: {failed}")
    return result


def _exactly_one(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one {pattern!r} below {root}, found {len(matches)}"
        )
    return matches[0]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return payload


def atomic_copy_verified(source: Path, target: Path) -> dict[str, Any]:
    """Publish an immutable spatial artifact, safely resuming after a crash."""

    source = source.resolve()
    source_sha = sha256(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if target.stat().st_size != source.stat().st_size or sha256(target) != source_sha:
            raise FileExistsError(f"existing immutable artifact differs: {target}")
    else:
        temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
        try:
            shutil.copy2(source, temporary)
            if temporary.stat().st_size != source.stat().st_size:
                raise IOError("temporary artifact copy has wrong byte count")
            if sha256(temporary) != source_sha:
                raise IOError("temporary artifact copy has wrong SHA-256")
            os.replace(temporary, target)
        finally:
            if temporary.exists():
                temporary.unlink()
    return {"path": str(target), "bytes": target.stat().st_size, "sha256": source_sha}


def compact_completed_event(
    event_dir: Path,
    event: dict[str, Any],
    *,
    run_fingerprint: str,
) -> dict[str, Any]:
    """Compact existing stage products without recalculating any science."""

    event_dir = event_dir.resolve()
    paths = {
        "prepare": _exactly_one(event_dir, PREPARE_MANIFEST_GLOB),
        "rainfall": _exactly_one(event_dir, RAIN_MANIFEST_GLOB),
        "wind": _exactly_one(event_dir, WIND_MANIFEST_GLOB),
        "road_overlay": _exactly_one(event_dir, ROAD_MANIFEST_GLOB),
        "road_summary": _exactly_one(event_dir, ROAD_SUMMARY_GLOB),
        "road_overlap_nc": _exactly_one(event_dir, ROAD_OVERLAP_NC_GLOB),
        "compact_hazard_footprint_nc": _exactly_one(
            event_dir, COMPACT_FOOTPRINT_NC_GLOB
        ),
    }
    payloads = {
        name: _load_json(paths[name])
        for name in ("prepare", "rainfall", "wind", "road_overlay", "road_summary")
    }
    event_id = str(event["event_id"])
    expected_status = {
        "prepare": "prepared_ready",
        "rainfall": "completed",
        "wind": "completed",
        "road_overlay": "completed",
    }
    for stage, status in expected_status.items():
        manifest = payloads[stage]
        if manifest.get("status") != status:
            raise ValueError(f"{stage} status is not {status!r}")
        manifest_event_id = (
            manifest.get("event", {}).get("event_id")
            if stage == "prepare"
            else manifest.get("event_id")
        )
        if str(manifest_event_id) != event_id:
            raise ValueError(f"{stage} event_id differs from the sample")

    road_summary = payloads["road_summary"]
    if str(road_summary.get("event_id")) != event_id:
        raise ValueError("road summary event_id differs from the sample")
    definition = road_summary.get("definition", {})
    if definition.get("hazard_thresholds_applied") is not False:
        raise ValueError("road summary must remain threshold-free")
    if definition.get("damage_or_loss_model_applied") is not False:
        raise ValueError("road summary must remain loss-model-free")
    rows = road_summary.get("road_classes")
    if not isinstance(rows, list) or len(rows) != len(ROAD_CLASS_NAMES):
        raise ValueError("road summary must contain exactly five road classes")
    observed_names = tuple(str(row.get("road_class_name")) for row in rows)
    if observed_names != ROAD_CLASS_NAMES:
        raise ValueError(f"unexpected road-class order: {observed_names}")

    return {
        "schema_version": "1.0",
        "status": "completed",
        "completed_at_utc": utc_now(),
        "run_fingerprint": run_fingerprint,
        "event": dict(event),
        "interpretation": {
            "descriptive_exposure_only": True,
            "hazard_thresholds_applied": False,
            "damage_or_loss_model_applied": False,
            "right_censored_accumulation": (
                "available-window lower bound"
                if bool(event["right_censored_at_15_day_limit"])
                else "not right censored"
            ),
        },
        "road_exposure_definition": definition,
        "road_classes": rows,
        "stage_diagnostics": {
            stage: {
                "manifest_sha256": sha256(paths[stage]),
                "manifest_bytes": paths[stage].stat().st_size,
                "result_summary": payloads[stage].get("result_summary")
                or payloads[stage].get("raw_result_summary"),
            }
            for stage in ("prepare", "rainfall", "wind", "road_overlay")
        },
        "spatial_artifact_sources": {
            name: {
                "path": str(paths[name]),
                "bytes": paths[name].stat().st_size,
                "sha256": sha256(paths[name]),
            }
            for name in ("road_overlap_nc", "compact_hazard_footprint_nc")
        },
    }


@contextmanager
def event_lock(path: Path) -> Iterator[None]:
    """Prevent duplicate work while allowing automatic crash recovery."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _record_path(output_root: Path, event_position: int) -> Path:
    return output_root / "events" / f"{event_position:05d}.json"


def _completed_for_fingerprint(path: Path, fingerprint: str) -> bool:
    if not path.is_file():
        return False
    try:
        record = _load_json(path)
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return (
        record.get("status") == "completed"
        and record.get("run_fingerprint") == fingerprint
    )


def _scientific_audit_for_fingerprint(path: Path, fingerprint: str) -> bool:
    if not path.is_file():
        return False
    try:
        record = _load_json(path)
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return (
        record.get("status") == "METHOD_DOMAIN_PENDING"
        and record.get("run_fingerprint") == fingerprint
    )


def _load_worker_scientific_audit(
    event_output: Path, event: dict[str, Any]
) -> dict[str, Any] | None:
    """Read a structured worker audit request; never infer science from stderr."""

    matches = sorted(event_output.glob(SCIENTIFIC_AUDIT_REPORT_GLOB))
    if not matches:
        return None
    if len(matches) != 1:
        raise ValueError("event worker emitted multiple scientific-audit reports")
    payload = _load_json(matches[0])
    if payload.get("status") != "scientific_audit_required":
        raise ValueError("worker scientific-audit report has invalid status")
    if str(payload.get("event_id")) != str(event["event_id"]):
        raise ValueError("worker scientific-audit report has wrong event_id")
    if int(payload.get("event_position", -1)) != int(event["event_position"]):
        raise ValueError("worker scientific-audit report has wrong event_position")
    if not str(payload.get("audit_code", "")):
        raise ValueError("worker scientific-audit report lacks audit_code")
    if payload.get("method_status") != "METHOD_DOMAIN_PENDING":
        raise ValueError("worker scientific-audit report has wrong method_status")
    return {
        "report": payload,
        "report_sha256": sha256(matches[0]),
        "report_bytes": matches[0].stat().st_size,
    }


def _worker_command(
    worker: Path,
    sample: Path,
    fixed_r0_catalogue: Path,
    fixed_r0_manifest: Path,
    event: dict[str, Any],
    event_output_dir: Path,
    worker_args: Sequence[str],
) -> list[str]:
    command = [
        sys.executable,
        str(worker),
        "--sample",
        str(sample),
        "--fixed-r0-catalogue",
        str(fixed_r0_catalogue),
        "--fixed-r0-manifest",
        str(fixed_r0_manifest),
        "--event-position",
        str(event["event_position"]),
        "--output-dir",
        str(event_output_dir),
    ]
    command.extend(worker_args)
    return command


def run_one_event(
    *,
    event: dict[str, Any],
    sample: Path,
    fixed_r0_catalogue: Path,
    fixed_r0_manifest: Path,
    worker: Path,
    worker_args: Sequence[str],
    output_root: Path,
    scratch_root: Path,
    run_fingerprint: str,
    retain_full_fields: bool,
) -> str:
    """Run, validate, compact, and then release one event's temporary fields."""

    position = int(event["event_position"])
    record_path = _record_path(output_root, position)
    lock_path = output_root / "locks" / f"{position:05d}.lock"
    with event_lock(lock_path):
        if _completed_for_fingerprint(record_path, run_fingerprint):
            return "skipped_completed"
        if _scientific_audit_for_fingerprint(record_path, run_fingerprint):
            return "skipped_scientific_audit"

        attempt_id = uuid.uuid4().hex
        attempt_started = utc_now()
        scratch_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f"event-{position:05d}-", dir=scratch_root
        ) as temporary_text:
            temporary = Path(temporary_text)
            event_output = temporary / "pipeline"
            command = _worker_command(
                worker,
                sample,
                fixed_r0_catalogue,
                fixed_r0_manifest,
                event,
                event_output,
                worker_args,
            )
            process: subprocess.CompletedProcess[str] | None = None
            try:
                process = subprocess.run(
                    command,
                    check=False,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                if process.returncode != 0:
                    raise RuntimeError(
                        f"event worker exited with status {process.returncode}"
                    )
                completed = compact_completed_event(
                    event_output, event, run_fingerprint=run_fingerprint
                )
                completed.update(
                    attempt_id=attempt_id,
                    attempt_started_at_utc=attempt_started,
                )
                spatial_sources = completed.pop("spatial_artifact_sources")
                completed["permanent_spatial_artifacts"] = {
                    "road_overlap_nc": atomic_copy_verified(
                        Path(spatial_sources["road_overlap_nc"]["path"]),
                        output_root / "road_overlap" / f"{position:05d}.nc",
                    ),
                    "compact_hazard_footprint_nc": atomic_copy_verified(
                        Path(
                            spatial_sources["compact_hazard_footprint_nc"]["path"]
                        ),
                        output_root
                        / "compact_hazard_footprint"
                        / f"{position:05d}.nc",
                    ),
                }
                if retain_full_fields:
                    retained = (
                        output_root
                        / "qa_events"
                        / f"{position:05d}-{event['event_id']}"
                        / run_fingerprint
                    )
                    if retained.exists():
                        raise FileExistsError(
                            f"retained QA target already exists: {retained}"
                        )
                    retained.parent.mkdir(parents=True, exist_ok=True)
                    staging = retained.with_name(f".{retained.name}.{attempt_id}.tmp")
                    shutil.copytree(event_output, staging)
                    os.replace(staging, retained)
                    completed["retained_full_fields"] = str(retained)
                else:
                    completed["retained_full_fields"] = None
                atomic_json(record_path, completed)
                atomic_json(
                    output_root
                    / "attempts"
                    / f"{position:05d}"
                    / f"{attempt_id}.json",
                    completed,
                )
                return "completed"
            except Exception as error:
                scientific_audit = _load_worker_scientific_audit(
                    event_output, event
                )
                failure = {
                    "schema_version": "1.0",
                    "status": (
                        "METHOD_DOMAIN_PENDING"
                        if scientific_audit is not None
                        else "failed"
                    ),
                    "failed_at_utc": utc_now(),
                    "attempt_started_at_utc": attempt_started,
                    "attempt_id": attempt_id,
                    "run_fingerprint": run_fingerprint,
                    "event": dict(event),
                    "error_type": type(error).__name__,
                    "error_message": str(error),
                    "traceback": traceback.format_exc(),
                    "worker_returncode": (
                        None if process is None else int(process.returncode)
                    ),
                    "worker_stdout": None if process is None else process.stdout,
                    "worker_stderr": None if process is None else process.stderr,
                    "temporary_fields_retained": False,
                }
                if scientific_audit is not None:
                    failure["scientific_audit"] = scientific_audit
                atomic_json(record_path, failure)
                atomic_json(
                    output_root
                    / "attempts"
                    / f"{position:05d}"
                    / f"{attempt_id}.json",
                    failure,
                )
                if scientific_audit is not None:
                    atomic_json(
                        output_root / "audit_queue" / f"{position:05d}.json",
                        failure,
                    )
                    return "method_domain_pending"
                return "failed"


def build_run_fingerprint(
    sample_record: dict[str, Any],
    fixed_r0_record: dict[str, Any],
    worker: Path,
    worker_args: Sequence[str],
) -> str:
    payload = {
        "batch_script_sha256": sha256(Path(__file__).resolve()),
        "sample_sha256": sample_record["sha256"],
        "fixed_r0_catalogue_sha256": fixed_r0_record["sha256"],
        "fixed_r0_manifest_sha256": fixed_r0_record["manifest_sha256"],
        "fixed_r0_distribution_contract_sha256": fixed_r0_record[
            "distribution_contract_sha256"
        ],
        "fixed_r0_event_outer_radius_binding_sha256": fixed_r0_record[
            "event_outer_radius_binding_sha256"
        ],
        "track_native_last_index": sample_record["track_native_last_index"],
        "worker_sha256": sha256(worker),
        "worker_args": list(worker_args),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--fixed-r0-catalogue", type=Path, required=True)
    parser.add_argument("--fixed-r0-manifest", type=Path, required=True)
    parser.add_argument("--event-worker", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--track-native-last-index", type=int, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--worker-arg", action="append", default=[])
    parser.add_argument(
        "--retain-qa-event-position", action="append", type=int, default=[]
    )
    parser.add_argument(
        "--allow-non-10k-sample",
        action="store_true",
        help="accept a road-domain or other non-10k event sample",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample = args.sample.resolve()
    fixed_r0_catalogue = args.fixed_r0_catalogue.resolve()
    fixed_r0_manifest = args.fixed_r0_manifest.resolve()
    worker = args.event_worker.resolve()
    output_root = args.output_root.resolve()
    scratch_root = args.scratch_root.resolve()
    if not worker.is_file():
        raise FileNotFoundError(worker)
    events, sample_record = load_sample_events(
        sample,
        track_native_last_index=args.track_native_last_index,
        require_10k=not args.allow_non_10k_sample,
    )
    fixed_r0_record = validate_fixed_r0_catalogue(
        fixed_r0_catalogue,
        fixed_r0_manifest,
        events,
        sample_record,
    )
    qa_coverage_positions = representative_qa_positions(
        events, require_frozen_10k_identities=not args.allow_non_10k_sample
    )
    selected = select_shard(events, args.shard_index, args.shard_count)
    selected_positions = {int(event["event_position"]) for event in selected}
    retain = set(args.retain_qa_event_position)
    unknown_qa = retain - selected_positions
    if unknown_qa:
        raise ValueError(
            f"QA positions are outside this shard: {sorted(unknown_qa)}"
        )
    fingerprint = build_run_fingerprint(
        sample_record, fixed_r0_record, worker, args.worker_arg
    )
    counts = {
        "completed": 0,
        "failed": 0,
        "method_domain_pending": 0,
        "skipped_completed": 0,
        "skipped_scientific_audit": 0,
    }
    for event in selected:
        outcome = run_one_event(
            event=event,
            sample=sample,
            fixed_r0_catalogue=fixed_r0_catalogue,
            fixed_r0_manifest=fixed_r0_manifest,
            worker=worker,
            worker_args=args.worker_arg,
            output_root=output_root,
            scratch_root=scratch_root,
            run_fingerprint=fingerprint,
            retain_full_fields=int(event["event_position"]) in retain,
        )
        counts[outcome] += 1

    records = []
    for event in selected:
        path = _record_path(output_root, int(event["event_position"]))
        if path.is_file():
            record = _load_json(path)
            records.append(
                {
                    "event_position": int(event["event_position"]),
                    "event_id": event["event_id"],
                    "status": record.get("status"),
                    "run_fingerprint": record.get("run_fingerprint"),
                }
            )
    records.sort(key=lambda item: item["event_position"])
    shard_manifest = {
        "schema_version": "1.0",
        "script_version": SCRIPT_VERSION,
        "generated_at_utc": utc_now(),
        "run_fingerprint": fingerprint,
        "sample": sample_record,
        "fixed_r0_catalogue": fixed_r0_record,
        "representative_qa_coverage_positions": qa_coverage_positions,
        "representative_qa_note": (
            "coverage selections only; no event is screened or weighted by QA status"
        ),
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "selected_event_count": len(selected),
        "right_censored_event_count_in_shard": sum(
            int(event["right_censored_at_15_day_limit"]) for event in selected
        ),
        "attempt_outcomes_this_invocation": counts,
        "events": records,
    }
    shard_path = (
        output_root
        / "shards"
        / f"shard-{args.shard_index:05d}-of-{args.shard_count:05d}.json"
    )
    atomic_json(shard_path, shard_manifest)
    print(json.dumps(shard_manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
