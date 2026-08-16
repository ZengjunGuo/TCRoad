#!/usr/bin/env python3
"""Extract motor-road ways from the frozen planet-260803 PBF.

Run this on the server that holds
`data/osm/derived/planet-260803/roads.osm.pbf` (or the parent planet).
It writes row-group parquet/CSV shards with tags, length and a
representative point.  It does not assign countries or dollars.

Requires osmium/pyosmium.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, TextIO


ACCEPTED = {
    "motorway",
    "motorway_link",
    "trunk",
    "trunk_link",
    "primary",
    "primary_link",
    "secondary",
    "secondary_link",
    "tertiary",
    "tertiary_link",
    "residential",
    "unclassified",
}

CSV_FIELDS = [
    "way_id",
    "highway",
    "lanes",
    "surface",
    "bridge",
    "tunnel",
    "lit",
    "n_nodes",
    "length_km",
    "lon",
    "lat",
]


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    radius = 6371.0088
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2.0) ** 2
    return 2.0 * radius * math.asin(min(1.0, math.sqrt(a)))


def way_length_and_mid(coords: list[tuple[float, float]]) -> tuple[float, float, float]:
    if len(coords) < 2:
        lon, lat = coords[0]
        return 0.0, lon, lat
    total = 0.0
    parts = []
    for (lon1, lat1), (lon2, lat2) in zip(coords, coords[1:]):
        segment = haversine_km(lon1, lat1, lon2, lat2)
        parts.append((total, segment, lon1, lat1, lon2, lat2))
        total += segment
    if total <= 0.0:
        lon, lat = coords[len(coords) // 2]
        return 0.0, lon, lat
    half = 0.5 * total
    for start, segment, lon1, lat1, lon2, lat2 in parts:
        if start + segment >= half and segment > 0.0:
            t = (half - start) / segment
            return total, lon1 + t * (lon2 - lon1), lat1 + t * (lat2 - lat1)
    lon, lat = coords[-1]
    return total, lon, lat


class MotorRoadHandler:
    def __init__(self, writer: Any, stats: dict[str, int]) -> None:
        self.writer = writer
        self.stats = stats

    def way(self, way) -> None:  # pyosmium callback
        highway = way.tags.get("highway")
        self.stats["ways_seen"] += 1
        if highway not in ACCEPTED:
            return
        if not way.is_closed() and way.nodes is None:
            return
        try:
            coords = [(float(n.lon), float(n.lat)) for n in way.nodes if n.location.valid()]
        except Exception:
            self.stats["ways_missing_nodes"] += 1
            return
        if len(coords) < 2:
            self.stats["ways_too_short"] += 1
            return
        length_km, lon, lat = way_length_and_mid(coords)
        if length_km <= 0.0:
            self.stats["ways_zero_length"] += 1
            return
        self.writer.writerow(
            {
                "way_id": int(way.id),
                "highway": highway,
                "lanes": way.tags.get("lanes", ""),
                "surface": way.tags.get("surface", ""),
                "bridge": way.tags.get("bridge", ""),
                "tunnel": way.tags.get("tunnel", ""),
                "lit": way.tags.get("lit", ""),
                "n_nodes": len(coords),
                "length_km": f"{length_km:.6f}",
                "lon": f"{lon:.6f}",
                "lat": f"{lat:.6f}",
            }
        )
        self.stats["ways_written"] += 1
        self.stats["length_km_total"] += length_km


def extract_with_pyosmium(pbf: Path, dest: TextIO) -> dict[str, int]:
    import osmium  # type: ignore

    writer = csv.DictWriter(dest, fieldnames=CSV_FIELDS)
    writer.writeheader()
    stats = {
        "ways_seen": 0,
        "ways_written": 0,
        "ways_missing_nodes": 0,
        "ways_too_short": 0,
        "ways_zero_length": 0,
        "length_km_total": 0.0,
    }

    class Handler(osmium.SimpleHandler, MotorRoadHandler):
        def __init__(self) -> None:
            osmium.SimpleHandler.__init__(self)
            MotorRoadHandler.__init__(self, writer, stats)

    Handler().apply_file(str(pbf), locations=True)
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pbf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.pbf.is_file():
        raise FileNotFoundError(args.pbf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        stats = extract_with_pyosmium(args.pbf, handle)
    manifest = args.output.with_suffix(args.output.suffix + ".manifest.json")
    payload = {
        "pbf": str(args.pbf),
        "output": str(args.output),
        "osm_snapshot": "planet-260803",
        "osm_snapshot_date": "2026-08-03",
        "accepted_highway": sorted(ACCEPTED),
        "stats": stats,
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
