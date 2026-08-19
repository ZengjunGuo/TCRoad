"""Contract tests for the 10,000-event scheduling and compaction layer."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

from netCDF4 import Dataset
import numpy as np


CODE_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from run_lin_10k_batch import (  # noqa: E402
    ROAD_CLASS_NAMES,
    atomic_json,
    compact_completed_event,
    load_sample_events,
    run_one_event,
    representative_qa_positions,
    select_event_positions,
    select_shard,
    sha256,
    validate_fixed_r0_catalogue,
)
from build_lin_fixed_r0_catalog import build_catalogue  # noqa: E402


def _sample_fixture(path: Path) -> None:
    with Dataset(path, "w", format="NETCDF4") as sample:
        sample.createDimension("event", 4)

        def numeric(name: str, dtype: str, values: list[float | int]) -> None:
            variable = sample.createVariable(name, dtype, ("event",))
            variable[:] = np.asarray(values)

        event_id = sample.createVariable("event_id", str, ("event",))
        event_id[:] = np.asarray(["e0", "e1", "e2", "e3"], dtype=object)
        numeric("source_track_index", "i8", [2, 7, 8, 9])
        numeric("source_catalogue_event_position", "i8", [2, 7, 8, 9])
        numeric("task_year", "i8", [1995, 1995, 1996, 1996])
        numeric("threshold_genesis_native_index", "i4", [56, 24, 30, 40])
        numeric("threshold_lysis_native_index", "i4", [149, 360, 110, 360])
        numeric(
            "event_weight_climate_fixed_effect_ht_analysis_yr",
            "f8",
            [0.1, 0.2, 0.3, 0.4],
        )


def _fixed_r0_fixture(sample: Path, root: Path) -> tuple[Path, Path]:
    catalogue = root / "fixed_r0.nc"
    manifest = root / "fixed_r0.manifest.json"
    build_catalogue(
        sample,
        catalogue,
        manifest,
        expected_sample_sha256=sha256(sample),
        expected_count=4,
    )
    return catalogue, manifest


class LinTenThousandBatchContractTest(unittest.TestCase):
    def test_right_censored_events_are_flagged_retained_and_lower_bounds(self):
        with tempfile.TemporaryDirectory() as directory:
            sample = Path(directory) / "sample.nc"
            _sample_fixture(sample)
            events, record = load_sample_events(
                sample, track_native_last_index=360, require_10k=False
            )

        self.assertEqual(len(events), 4)
        self.assertEqual(record["right_censored_event_count"], 2)
        self.assertFalse(events[0]["right_censored_at_15_day_limit"])
        self.assertTrue(events[1]["right_censored_at_15_day_limit"])
        self.assertEqual(
            events[1]["cumulative_hazard_interpretation"],
            "lower_bound_over_available_window",
        )
        # Right-censored events participate in the same position-only sharding;
        # there is no hazard or censoring screen.
        shard = select_shard(events, shard_index=1, shard_count=2)
        self.assertEqual([event["event_id"] for event in shard], ["e1", "e3"])
        self.assertTrue(
            all(event["right_censored_at_15_day_limit"] for event in shard)
        )

    def test_select_event_positions_requires_every_requested_id(self):
        events = [{"event_position": 10}, {"event_position": 11}, {"event_position": 12}]
        selected = select_event_positions(events, [12, 10])
        self.assertEqual([int(item["event_position"]) for item in selected], [10, 12])
        with self.assertRaises(ValueError):
            select_event_positions(events, [10, 99])

    def test_sharding_is_order_independent_and_exhaustive(self):
        events = [{"event_position": position} for position in range(37)]
        selected = [
            event["event_position"]
            for shard_index in range(6)
            for event in select_shard(events[::-1], shard_index, 6)
        ]
        self.assertEqual(sorted(selected), list(range(37)))
        self.assertEqual(len(selected), len(set(selected)))

    def test_fixed_r0_catalogue_exactly_binds_sample_order_and_enters_events(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.nc"
            _sample_fixture(sample)
            catalogue, manifest = _fixed_r0_fixture(sample, root)
            events, sample_record = load_sample_events(
                sample, track_native_last_index=360, require_10k=False
            )
            record = validate_fixed_r0_catalogue(
                catalogue, manifest, events, sample_record
            )

        self.assertEqual(record["event_count"], 4)
        self.assertTrue(record["distribution_contract_sha256"])
        self.assertGreater(events[0]["outer_radius_m"], 0.0)
        self.assertEqual(
            events[0]["fixed_r0_catalogue_sha256"], record["sha256"]
        )

    def test_representative_qa_positions_cover_frozen_edge_cases(self):
        events = []
        for position in range(9842):
            events.append(
                {
                    "event_position": position,
                    "event_id": f"e{position}",
                    "outer_radius_m": 800_000.0,
                    "right_censored_at_15_day_limit": False,
                    "available_hour_count": 48,
                }
            )
        events[4894]["outer_radius_m"] = 308_736.0
        events[9841]["outer_radius_m"] = 2_335_623.0
        events[1]["event_id"] = "stream0000-year1995-track000007"
        events[2649]["event_id"] = "stream0000-year2000-track026613"
        events[6]["right_censored_at_15_day_limit"] = True
        events[103]["available_hour_count"] = 9
        events[953]["event_id"] = "stream0000-year1996-track009518"

        selected = representative_qa_positions(events)
        self.assertEqual(selected["minimum_fixed_r0"], 4894)
        self.assertEqual(selected["maximum_fixed_r0"], 9841)
        self.assertEqual(selected["lowest_absolute_track_latitude"], 1732)
        self.assertEqual(selected["periodic_longitude_seam_crossing"], 2649)
        self.assertEqual(selected["former_nonpositive_size_predictor_case"], 953)

    def test_event0_compaction_is_exact_for_all_five_road_classes(self):
        source_event_dir = PROJECT_DIR / "validation" / "lin_event0"
        source_summary_path = (
            source_event_dir
            / "road_overlap"
            / "lin_event0_road_class_joint_exposure_summary.json"
        )
        source_summary = json.loads(source_summary_path.read_text())
        prepare_manifest = json.loads(
            (
                source_event_dir
                / "lin_event0_public_inputs_prepare_only.manifest.json"
            ).read_text()
        )
        event = {
            **prepare_manifest["event"],
            "right_censored_at_15_day_limit": False,
            "cumulative_hazard_interpretation": (
                "complete_threshold_window_accumulation"
            ),
        }
        with tempfile.TemporaryDirectory() as directory:
            event_dir = Path(directory) / "event0"
            shutil.copytree(source_event_dir, event_dir)
            (event_dir / "lin_event0_compact_hazard_footprint.nc").write_bytes(
                b"compact-fixture"
            )
            compact = compact_completed_event(
                event_dir, event, run_fingerprint="event0-equivalence-test"
            )

        self.assertEqual(compact["road_classes"], source_summary["road_classes"])
        self.assertEqual(
            compact["road_exposure_definition"], source_summary["definition"]
        )
        self.assertEqual(
            tuple(row["road_class_name"] for row in compact["road_classes"]),
            ROAD_CLASS_NAMES,
        )
        self.assertFalse(compact["interpretation"]["hazard_thresholds_applied"])
        self.assertFalse(
            compact["interpretation"]["damage_or_loss_model_applied"]
        )

    def test_atomic_json_never_leaves_a_partial_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events" / "00000.json"
            atomic_json(path, {"status": "failed", "attempt": 1})
            atomic_json(path, {"status": "completed", "attempt": 2})
            payload = json.loads(path.read_text())
            leftovers = list(path.parent.glob(".*.tmp"))

        self.assertEqual(payload, {"status": "completed", "attempt": 2})
        self.assertEqual(leftovers, [])

    def test_structured_fixed_r0_domain_failure_enters_audit_queue(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            worker = root / "audit_worker.py"
            worker.write_text(
                """import argparse, json
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--sample')
p.add_argument('--fixed-r0-catalogue')
p.add_argument('--fixed-r0-manifest')
p.add_argument('--event-position', type=int)
p.add_argument('--output-dir')
a = p.parse_args()
o = Path(a.output_dir)
o.mkdir(parents=True)
json.dump(
    {
        'status': 'scientific_audit_required',
        'method_status': 'METHOD_DOMAIN_PENDING',
        'event_id': 'e15',
        'event_position': a.event_position,
        'audit_code': 'C15_FIXED_R0_BELOW_302KM_NUMERICAL_DOMAIN',
        'resolution_applied': False,
    },
    open(o / 'scientific_audit.json', 'w'),
)
raise SystemExit(3)
""",
                encoding="utf-8",
            )
            sample = root / "sample.nc"
            sample.write_bytes(b"sample identity only")
            fixed_r0_catalogue = root / "fixed_r0.nc"
            fixed_r0_catalogue.write_bytes(b"catalogue identity only")
            fixed_r0_manifest = root / "fixed_r0.manifest.json"
            fixed_r0_manifest.write_text("{}")
            output = root / "output"
            scratch = root / "scratch"
            event = {
                "event_position": 15,
                "event_id": "e15",
                "right_censored_at_15_day_limit": False,
                "cumulative_hazard_interpretation": (
                    "complete_threshold_window_accumulation"
                ),
            }
            outcome = run_one_event(
                event=event,
                sample=sample,
                fixed_r0_catalogue=fixed_r0_catalogue,
                fixed_r0_manifest=fixed_r0_manifest,
                worker=worker,
                worker_args=[],
                output_root=output,
                scratch_root=scratch,
                run_fingerprint="same-science-contract",
                retain_full_fields=False,
            )
            record = json.loads((output / "events" / "00015.json").read_text())
            queued = json.loads(
                (output / "audit_queue" / "00015.json").read_text()
            )
            resumed = run_one_event(
                event=event,
                sample=sample,
                fixed_r0_catalogue=fixed_r0_catalogue,
                fixed_r0_manifest=fixed_r0_manifest,
                worker=worker,
                worker_args=[],
                output_root=output,
                scratch_root=scratch,
                run_fingerprint="same-science-contract",
                retain_full_fields=False,
            )

        self.assertEqual(outcome, "method_domain_pending")
        self.assertEqual(record["status"], "METHOD_DOMAIN_PENDING")
        self.assertEqual(
            record["scientific_audit"]["report"]["audit_code"],
            "C15_FIXED_R0_BELOW_302KM_NUMERICAL_DOMAIN",
        )
        self.assertFalse(
            record["scientific_audit"]["report"]["resolution_applied"]
        )
        self.assertEqual(queued, record)
        self.assertEqual(resumed, "skipped_scientific_audit")


if __name__ == "__main__":
    unittest.main()
