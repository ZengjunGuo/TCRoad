#!/usr/bin/env python3
"""Apply the 2025 unit-cost book to an extracted motor-road table.

Input CSV must contain: way_id, highway, length_km, and iso3.
Optional: lanes, surface, bridge, tunnel, slope_deg.

A first pass computes country × class lane medians from tagged rows, then
each row is valued.  Two totals are always printed: all motor roads, and
all motor roads except local.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import statistics
import sys
from typing import Any


CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from assign_road_countries import assign_rows as assign_iso3_from_points  # noqa: E402
from build_road_unit_cost_book import assemble_book  # noqa: E402
from road_replacement_value import classify_highway, parse_lanes, replacement_cost  # noqa: E402

UNCLASSIFIED_SHARE_FLAG = 0.60


def first_pass_lane_medians(rows: list[dict[str, str]]) -> dict[str, dict[str, float]]:
    buckets: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        iso = row.get("iso3", "")
        road_class, is_link, _ = classify_highway(row.get("highway", ""))
        if road_class is None or is_link or not iso:
            continue
        lanes = parse_lanes(row.get("lanes") or None)
        if lanes is None:
            continue
        buckets[iso][road_class].append(lanes)
    return {
        iso: {cls: statistics.median(vals) for cls, vals in classes.items() if vals}
        for iso, classes in buckets.items()
    }


def qa_replacement_totals(valued: list[dict[str, Any]]) -> dict[str, Any]:
    """Country and class QA, including the GIRI unclassified-share flag."""

    by_country: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "length_km": 0.0,
            "unclassified_km": 0.0,
            "local_km": 0.0,
            "replacement_usd": 0.0,
            "ways": 0.0,
        }
    )
    by_class_km: dict[str, float] = defaultdict(float)
    for row in valued:
        if not int(row.get("accepted") or 0):
            continue
        iso = row.get("iso3") or "UNK"
        length = float(row["length_km"])
        bucket = by_country[iso]
        bucket["length_km"] += length
        bucket["replacement_usd"] += float(row["replacement_usd"])
        bucket["ways"] += 1.0
        road_class = row.get("road_class") or ""
        by_class_km[road_class] += length
        if road_class == "local":
            bucket["local_km"] += length
        highway = str(row.get("highway") or "")
        if highway == "unclassified" or (
            road_class == "local" and highway.endswith("unclassified")
        ):
            bucket["unclassified_km"] += length

    flags = []
    country_table = {}
    for iso, bucket in sorted(by_country.items()):
        share = (
            bucket["unclassified_km"] / bucket["length_km"]
            if bucket["length_km"] > 0.0
            else 0.0
        )
        flagged = share > UNCLASSIFIED_SHARE_FLAG
        country_table[iso] = {
            "length_km": bucket["length_km"],
            "unclassified_km": bucket["unclassified_km"],
            "unclassified_share": share,
            "local_km": bucket["local_km"],
            "replacement_usd": bucket["replacement_usd"],
            "ways": int(bucket["ways"]),
            "unclassified_share_flag": flagged,
        }
        if flagged:
            flags.append(iso)
    return {
        "by_country": country_table,
        "by_class_km": dict(by_class_km),
        "unclassified_flag_threshold": UNCLASSIFIED_SHARE_FLAG,
        "countries_flagged_unclassified": flags,
    }


def apply_rows(rows: list[dict[str, str]], repo: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = assign_iso3_from_points(rows)
    book, meta = assemble_book(repo)
    medians = first_pass_lane_medians(rows)
    for iso, table in medians.items():
        if iso in book.countries:
            book.countries[iso].lane_median_by_class = table

    valued = []
    totals = {
        "accepted_ways": 0,
        "rejected_ways": 0,
        "length_km": 0.0,
        "length_km_no_local": 0.0,
        "replacement_usd": 0.0,
        "replacement_usd_low": 0.0,
        "replacement_usd_high": 0.0,
        "replacement_usd_no_local": 0.0,
        "by_class_usd": defaultdict(float),
        "by_book_usd": defaultdict(float),
    }
    for row in rows:
        result = replacement_cost(
            book,
            length_km=float(row["length_km"]),
            highway=row.get("highway", ""),
            iso3=row.get("iso3", ""),
            lanes=row.get("lanes") or None,
            surface=row.get("surface") or None,
            bridge=row.get("bridge") or None,
            tunnel=row.get("tunnel") or None,
            slope_deg=row.get("slope_deg") or None,
        )
        record = {
            "way_id": row.get("way_id", ""),
            "iso3": row.get("iso3", ""),
            "highway": row.get("highway", ""),
            "accepted": int(result.accepted),
            "reason": result.reason,
            "road_class": result.road_class,
            "is_link": int(result.is_link),
            "is_bridge": int(result.is_bridge),
            "is_tunnel": int(result.is_tunnel),
            "surface": result.surface,
            "lanes_used": result.lanes_used if result.lanes_used is not None else "",
            "lanes_source": result.lanes_source,
            "terrain_class": result.terrain_class,
            "work_type": result.work_type,
            "price_book": result.price_book,
            "length_km": result.length_km,
            "usd_per_km": result.usd_per_km,
            "replacement_usd": result.replacement_usd,
            "replacement_usd_low": result.replacement_usd_low,
            "replacement_usd_high": result.replacement_usd_high,
        }
        valued.append(record)
        if not result.accepted:
            totals["rejected_ways"] += 1
            continue
        totals["accepted_ways"] += 1
        totals["length_km"] += result.length_km
        totals["replacement_usd"] += result.replacement_usd
        totals["replacement_usd_low"] += result.replacement_usd_low
        totals["replacement_usd_high"] += result.replacement_usd_high
        totals["by_class_usd"][result.road_class] += result.replacement_usd
        totals["by_book_usd"][result.price_book] += result.replacement_usd
        if result.road_class != "local":
            totals["length_km_no_local"] += result.length_km
            totals["replacement_usd_no_local"] += result.replacement_usd

    totals["by_class_usd"] = dict(totals["by_class_usd"])
    totals["by_book_usd"] = dict(totals["by_book_usd"])
    totals["unit_cost_meta"] = meta
    totals["countries_with_lane_median"] = {
        iso: table for iso, table in medians.items()
    }
    totals.update(qa_replacement_totals(valued))
    return valued, totals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=ROOT)
    args = parser.parse_args()
    with args.input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    valued, totals = apply_rows(rows, args.repo)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(valued[0].keys()) if valued else ["way_id"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(valued)
    summary = args.output.with_suffix(args.output.suffix + ".summary.json")
    summary.write_text(json.dumps(totals, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: totals[k] for k in totals if k != "countries_with_lane_median"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
