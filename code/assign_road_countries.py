#!/usr/bin/env python3
"""Assign ISO3 country codes to road points using Natural Earth 50m admin-0.

Stdlib shapefile reader + even-odd point-in-polygon.  No GDAL: this machine
has mixed Homebrew/Anaconda PROJ databases that make OGR raise on a simple
layer scan.

Natural Earth 50m omits some harbour pixels (e.g. lower Manhattan).  Coastal
misses are empty ISO3, not a wrong country.  A 10 m DEM-side pass can fill
those later; do not invent a nearest-country guess tonight.
"""

from __future__ import annotations

import argparse
import csv
import struct
from pathlib import Path
from typing import Optional


CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
DEFAULT_SHAPE = (
    ROOT / "vendor" / "naturalearth" / "50m" / "ne_50m_admin_0_countries.shp"
)

SHAPE_POLYGON = 5
SHAPE_POLYGON_Z = 15
SHAPE_POLYGON_M = 25


class CountryIndex:
    def __init__(self, shapefile: Path = DEFAULT_SHAPE) -> None:
        shapefile = shapefile.resolve()
        if not shapefile.is_file():
            raise FileNotFoundError(shapefile)
        records = _read_dbf_iso(shapefile.with_suffix(".dbf"))
        geometries = _read_polygons(shapefile)
        if len(records) != len(geometries):
            raise RuntimeError(
                f"DBF/SHP length mismatch: {len(records)} vs {len(geometries)}"
            )
        self._features: list[tuple[tuple[float, float, float, float], list[list[tuple[float, float]]], str]] = []
        for iso, geom in zip(records, geometries):
            if not iso or geom is None:
                continue
            box, rings = geom
            self._features.append((box, rings, iso))

    def iso3(self, lon: float, lat: float) -> str:
        if not _finite(lon) or not _finite(lat):
            return ""
        if lat < -90.0 or lat > 90.0:
            return ""
        lon = ((float(lon) + 180.0) % 360.0) - 180.0
        lat = float(lat)
        for (minx, maxx, miny, maxy), rings, iso in self._features:
            if lon < minx or lon > maxx or lat < miny or lat > maxy:
                continue
            if _point_in_rings(lon, lat, rings):
                return iso
        return ""


def _read_dbf_iso(path: Path) -> list[str]:
    data = path.read_bytes()
    nrecords = struct.unpack_from("<I", data, 4)[0]
    header_len = struct.unpack_from("<H", data, 8)[0]
    rec_len = struct.unpack_from("<H", data, 10)[0]
    fields = []
    pos = 32
    while pos < header_len - 1 and data[pos] != 0x0D:
        name = data[pos : pos + 11].split(b"\x00", 1)[0].decode("latin1")
        length = data[pos + 16]
        fields.append((name, length))
        pos += 32
    iso_off = adm_off = None
    cursor = 1  # skip deletion flag
    for name, length in fields:
        if name == "ISO_A3":
            iso_off = (cursor, length)
        elif name == "ADM0_A3":
            adm_off = (cursor, length)
        cursor += length
    if iso_off is None:
        raise RuntimeError("ISO_A3 missing from Natural Earth DBF")
    out = []
    body = data[header_len:]
    for i in range(nrecords):
        rec = body[i * rec_len : (i + 1) * rec_len]
        iso = _clean_iso(rec[iso_off[0] : iso_off[0] + iso_off[1]].decode("latin1"))
        if not iso and adm_off is not None:
            iso = _clean_iso(rec[adm_off[0] : adm_off[0] + adm_off[1]].decode("latin1"))
        out.append(iso)
    return out


def _read_polygons(path: Path):
    data = path.read_bytes()
    shx = path.with_suffix(".shx")
    if not shx.is_file():
        raise FileNotFoundError(shx)
    index = shx.read_bytes()
    count = (len(index) - 100) // 8
    offsets = []
    for i in range(count):
        off_words = struct.unpack_from(">i", index, 100 + 8 * i)[0]
        offsets.append(off_words * 2)
    geoms = []
    for off in offsets:
        rec_len_words = struct.unpack_from(">i", data, off + 4)[0]
        payload = data[off + 8 : off + 8 + rec_len_words * 2]
        geoms.append(_parse_polygon_payload(payload))
    return geoms


def _parse_polygon_payload(payload: bytes):
    if len(payload) < 4:
        return None
    shape_type = struct.unpack_from("<i", payload, 0)[0]
    if shape_type in (0,):
        return None
    if shape_type not in {SHAPE_POLYGON, SHAPE_POLYGON_Z, SHAPE_POLYGON_M}:
        return None
    minx, miny, maxx, maxy = struct.unpack_from("<4d", payload, 4)
    nparts, npoints = struct.unpack_from("<ii", payload, 36)
    parts = list(struct.unpack_from(f"<{nparts}i", payload, 44))
    pts_off = 44 + 4 * nparts
    points = [
        struct.unpack_from("<2d", payload, pts_off + 16 * i) for i in range(npoints)
    ]
    parts.append(npoints)
    rings = [points[parts[i] : parts[i + 1]] for i in range(nparts)]
    return (minx, maxx, miny, maxy), rings


def _point_in_rings(lon: float, lat: float, rings: list[list[tuple[float, float]]]) -> bool:
    # Even-odd over all rings: holes flip the bit, which is what we want.
    inside = False
    for ring in rings:
        if _even_odd(lon, lat, ring):
            inside = not inside
    return inside


def _even_odd(x: float, y: float, ring: list[tuple[float, float]]) -> bool:
    if len(ring) < 3:
        return False
    inside = False
    x1, y1 = ring[0]
    for x2, y2 in ring[1:]:
        if (y1 > y) != (y2 > y):
            at_x = (x2 - x1) * (y - y1) / (y2 - y1 + 0.0) + x1
            if x < at_x:
                inside = not inside
        x1, y1 = x2, y2
    return inside


def _clean_iso(value: object) -> str:
    text = str(value or "").strip().upper()
    if len(text) != 3 or text in {"-99", "NAN", "NONE", "NULL"}:
        return ""
    return text


def _finite(value: object) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number == number


def assign_rows(
    rows: list[dict[str, str]], index: Optional[CountryIndex] = None
) -> list[dict[str, str]]:
    """Fill empty iso3 from lon/lat. Existing iso3 is left unchanged."""

    need = any(
        not (row.get("iso3") or "").strip()
        and _finite(row.get("lon"))
        and _finite(row.get("lat"))
        for row in rows
    )
    if not need:
        return rows
    lookup = index or CountryIndex()
    for row in rows:
        if (row.get("iso3") or "").strip():
            continue
        if not _finite(row.get("lon")) or not _finite(row.get("lat")):
            continue
        row["iso3"] = lookup.iso3(float(row["lon"]), float(row["lat"]))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shapefile", type=Path, default=DEFAULT_SHAPE)
    args = parser.parse_args()
    with args.input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assign_rows(rows, CountryIndex(args.shapefile))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["iso3"]
    if "iso3" not in fields:
        fields.append("iso3")
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
