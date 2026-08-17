#!/usr/bin/env python3
"""One program for motor-road 2025 replacement cost.

This is the only road-valuation script.  It:
  1. builds the 2025 unit-cost book,
  2. tells which country a lon/lat point is in,
  3. prices a table of road segments,
  4. can extract motor roads from the frozen OSM PBF on the server.

Usage:
  python3 road_replacement_value.py write-book
  python3 road_replacement_value.py value roads.csv --output valued.csv
  python3 road_replacement_value.py extract planet.osm.pbf --output ways.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import struct
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, TextIO


CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
SCRIPT_VERSION = "1.1.0"

HIGHWAY_TO_CLASS: dict[str, str] = {
    "motorway": "motorway",
    "motorway_link": "motorway",
    "trunk": "trunk",
    "trunk_link": "trunk",
    "primary": "primary",
    "primary_link": "primary",
    "secondary": "secondary",
    "secondary_link": "secondary",
    "tertiary": "tertiary",
    "tertiary_link": "tertiary",
    "residential": "local",
    "unclassified": "local",
}

LINK_TAGS = {
    "motorway_link",
    "trunk_link",
    "primary_link",
    "secondary_link",
    "tertiary_link",
}

EXCLUDED_HIGHWAY = {
    "footway",
    "path",
    "cycleway",
    "steps",
    "pedestrian",
    "track",
    "service",
    "living_street",
    "construction",
    "proposed",
    "bridleway",
    "corridor",
    "busway",
    "raceway",
    "road",
    "escape",
    "rest_area",
    "services",
}

PAVED_SURFACES = {
    "paved",
    "asphalt",
    "concrete",
    "concrete:plates",
    "concrete:lanes",
    "chipseal",
    "paving_stones",
    "sett",
    "cobblestone",
    "unhewn_cobblestone",
    "metal",
    "wood",
    "compacted",
}

UNPAVED_SURFACES = {
    "unpaved",
    "gravel",
    "fine_gravel",
    "dirt",
    "earth",
    "ground",
    "grass",
    "sand",
    "mud",
    "pebblestone",
    "rock",
    "stone",
    "ground",
    "dirt/sand",
}

DEFAULT_LANES = {
    "motorway": 4,
    "trunk": 2,
    "primary": 2,
    "secondary": 2,
    "tertiary": 2,
    "local": 1,
}

WORK_TYPE_4L = "New 4L Expressway"
WORK_TYPE_2L = "New 2L Highway"
WORK_TYPE_1L = "New 1L Road"

PRICE_BOOK_ROCKS = "rocks_region"
PRICE_BOOK_NATIONAL = "national"
PRICE_BOOK_EUROPE = "europe_gdp_scaled"
PRICE_BOOK_HIC_SCALE = "hic_gdp_scaled"

LANE_STEP = 0.25
LANE_FLOOR = 0.5


@dataclass(frozen=True)
class CostBand:
    central: float
    low: float
    high: float
    source: str = ""
    n_support: int = 0
    fill_rule: str = "direct"


@dataclass
class CountryRecord:
    iso3: str
    name: str
    wb_region: str
    rocks_region: str
    income_level: str
    price_book: str
    paved_share: float
    paved_share_source: str
    gdp_pc_ppp_2015: Optional[float] = None
    gdp_pc_ppp_latest: Optional[float] = None
    eu28_ratio_2015: Optional[float] = None
    usa_ratio_latest: Optional[float] = None
    lane_median_by_class: dict[str, float] = field(default_factory=dict)
    paved_fraction_by_class: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class ReplacementResult:
    accepted: bool
    reason: str
    road_class: str
    is_link: bool
    is_bridge: bool
    is_tunnel: bool
    surface: str
    lanes_used: Optional[float]
    lanes_source: str
    terrain_class: str
    terrain_multiplier: float
    work_type: str
    price_book: str
    usd_per_km: float
    usd_per_km_low: float
    usd_per_km_high: float
    length_km: float
    replacement_usd: float
    replacement_usd_low: float
    replacement_usd_high: float


def classify_highway(highway: Optional[str]) -> tuple[Optional[str], bool, str]:
    """Return (road_class, is_link, reason). road_class is None if excluded."""

    tag = (highway or "").strip().lower()
    if not tag:
        return None, False, "missing_highway_tag"
    if tag in EXCLUDED_HIGHWAY:
        return None, False, f"excluded:{tag}"
    road_class = HIGHWAY_TO_CLASS.get(tag)
    if road_class is None:
        return None, False, f"unmapped_highway:{tag}"
    return road_class, tag in LINK_TAGS, "ok"


def parse_lanes(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        number = float(value)
        if number > 0.0 and number < 32.0:
            return number
        return None
    text = str(value).strip().lower()
    if not text or text in {"none", "nan", "null"}:
        return None
    if ";" in text:
        parts = []
        for piece in text.split(";"):
            parsed = parse_lanes(piece)
            if parsed is not None:
                parts.append(parsed)
        if not parts:
            return None
        return sum(parts) / len(parts)
    try:
        number = float(text.replace(",", "."))
    except ValueError:
        return None
    if number > 0.0 and number < 32.0:
        return number
    return None


def classify_surface(surface: Optional[str]) -> str:
    """Return paved, unpaved, or unknown."""

    tag = (surface or "").strip().lower()
    if not tag:
        return "unknown"
    head = tag.split(";")[0].strip()
    if head in PAVED_SURFACES or head.startswith("asphalt") or head.startswith("concrete"):
        return "paved"
    if head in UNPAVED_SURFACES:
        return "unpaved"
    return "unknown"


def is_truthy_osm(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"", "no", "false", "0", "none"}:
        return False
    return True


def terrain_class_from_slope(slope_deg: Optional[float]) -> tuple[str, float]:
    if slope_deg is None:
        return "unknown", 1.0
    try:
        slope = float(slope_deg)
    except (TypeError, ValueError):
        return "unknown", 1.0
    if slope != slope:  # NaN
        return "unknown", 1.0
    if slope < 0.0:
        return "unknown", 1.0
    if slope <= 10.0:
        return "plain", 1.0
    if slope <= 25.0:
        return "rolling", 1.573394495412844
    return "steep", 1.944954128440367


def default_lanes_for(road_class: str, is_link: bool) -> int:
    if is_link:
        return 1
    return DEFAULT_LANES[road_class]


def work_type_for(
    road_class: str, *, lanes: Optional[float], surface: str
) -> str:
    if road_class == "motorway":
        return WORK_TYPE_4L
    if road_class == "trunk":
        if lanes is not None and lanes >= 4.0:
            return WORK_TYPE_4L
        return WORK_TYPE_2L
    if road_class in {"primary", "secondary"}:
        return WORK_TYPE_2L
    if road_class == "tertiary":
        if surface == "unpaved":
            return WORK_TYPE_1L
        return WORK_TYPE_2L
    return WORK_TYPE_1L


def lane_factor(lanes_used: float, default_lanes: float, step: float = LANE_STEP) -> float:
    factor = 1.0 + step * (lanes_used - default_lanes)
    return max(LANE_FLOOR, factor)


def _scale_band(band: CostBand, scale: float) -> CostBand:
    return CostBand(
        central=band.central * scale,
        low=band.low * scale,
        high=band.high * scale,
        source=band.source,
        n_support=band.n_support,
        fill_rule=band.fill_rule,
    )


class UnitCostBook:
    """Country-aware 2025 USD/km lookup."""

    def __init__(
        self,
        countries: Mapping[str, CountryRecord],
        rocks: Mapping[str, Mapping[str, CostBand]],
        national: Mapping[str, Mapping[str, CostBand]],
        europe_baseline_usd: Mapping[str, CostBand],
        usa_baseline: Mapping[str, CostBand],
        bridge: CostBand,
        tunnel: CostBand,
        hic_clip: tuple[float, float] = (0.4, 1.3),
    ) -> None:
        self.countries = dict(countries)
        self.rocks = {region: dict(types) for region, types in rocks.items()}
        self.national = {iso: dict(classes) for iso, classes in national.items()}
        self.europe_baseline_usd = dict(europe_baseline_usd)
        self.usa_baseline = dict(usa_baseline)
        self.bridge = bridge
        self.tunnel = tunnel
        self.hic_clip = hic_clip

    def country(self, iso3: str) -> CountryRecord:
        try:
            return self.countries[iso3]
        except KeyError as exc:
            raise KeyError(f"unknown country {iso3!r}") from exc

    def lookup_pavement(
        self,
        iso3: str,
        road_class: str,
        work_type: str,
    ) -> tuple[CostBand, str]:
        record = self.country(iso3)
        if record.price_book == PRICE_BOOK_NATIONAL:
            table = self.national.get(iso3)
            if table is None or road_class not in table:
                raise KeyError(f"national book missing {iso3} {road_class}")
            return table[road_class], PRICE_BOOK_NATIONAL

        if record.price_book == PRICE_BOOK_EUROPE:
            baseline = self.europe_baseline_usd[road_class]
            ratio = record.eu28_ratio_2015
            if ratio is None:
                ratio = 1.0
            return _scale_band(baseline, ratio), PRICE_BOOK_EUROPE

        if record.price_book == PRICE_BOOK_HIC_SCALE:
            baseline = self.usa_baseline[road_class]
            ratio = record.usa_ratio_latest
            if ratio is None:
                ratio = 1.0
            lo, hi = self.hic_clip
            ratio = min(hi, max(lo, ratio))
            return _scale_band(baseline, ratio), PRICE_BOOK_HIC_SCALE

        region = record.rocks_region
        try:
            band = self.rocks[region][work_type]
        except KeyError as exc:
            raise KeyError(
                f"ROCKS book missing {region} {work_type} for {iso3}"
            ) from exc
        return band, PRICE_BOOK_ROCKS


def replacement_cost(
    book: UnitCostBook,
    *,
    length_km: float,
    highway: str,
    iso3: str,
    lanes: Any = None,
    surface: Any = None,
    bridge: Any = None,
    tunnel: Any = None,
    slope_deg: Any = None,
) -> ReplacementResult:
    """Value one segment. Rejected roads return accepted=False and zero dollars."""

    length = float(length_km)
    if not (length > 0.0) or length != length:
        return _rejected("non_positive_length", length_km=0.0)

    road_class, is_link, reason = classify_highway(highway)
    if road_class is None:
        return _rejected(reason, length_km=length)

    try:
        record = book.country(iso3)
    except KeyError:
        return _rejected(f"unknown_country:{iso3}", length_km=length, road_class=road_class)

    if is_truthy_osm(bridge):
        band = book.bridge
        return _accepted(
            road_class=road_class,
            is_link=is_link,
            is_bridge=True,
            is_tunnel=False,
            surface="structure",
            lanes_used=None,
            lanes_source="not_used_for_bridge",
            terrain_class="not_used_for_structure",
            terrain_multiplier=1.0,
            work_type="bridge",
            price_book="giri_structure",
            band=band,
            length_km=length,
        )
    if is_truthy_osm(tunnel):
        band = book.tunnel
        return _accepted(
            road_class=road_class,
            is_link=is_link,
            is_bridge=False,
            is_tunnel=True,
            surface="structure",
            lanes_used=None,
            lanes_source="not_used_for_tunnel",
            terrain_class="not_used_for_structure",
            terrain_multiplier=1.0,
            work_type="tunnel",
            price_book="giri_structure",
            band=band,
            length_km=length,
        )

    parsed_lanes = parse_lanes(lanes)
    parent_default_lanes = float(DEFAULT_LANES[road_class])
    link_or_class_default = float(default_lanes_for(road_class, is_link))
    if parsed_lanes is not None:
        lanes_used = parsed_lanes
        lanes_source = "osm_tag"
    elif (not is_link) and road_class in record.lane_median_by_class:
        lanes_used = float(record.lane_median_by_class[road_class])
        lanes_source = "country_class_median"
    else:
        lanes_used = link_or_class_default
        lanes_source = "link_default" if is_link else "class_default"

    surface_state = classify_surface(None if surface is None else str(surface))
    class_paved = record.paved_fraction_by_class.get(road_class)
    if class_paved is None:
        class_paved = 1.0 if record.paved_share >= 0.5 else record.paved_share

    terrain_name, terrain_mult = terrain_class_from_slope(
        None if slope_deg is None else float(slope_deg) if _is_number(slope_deg) else None
    )

    if surface_state == "unknown":
        paved_band, book_name = book.lookup_pavement(
            iso3, road_class, work_type_for(road_class, lanes=lanes_used, surface="paved")
        )
        unpaved_band, _ = book.lookup_pavement(
            iso3, road_class, work_type_for(road_class, lanes=lanes_used, surface="unpaved")
        )
        weight = min(1.0, max(0.0, float(class_paved)))
        band = CostBand(
            central=weight * paved_band.central + (1.0 - weight) * unpaved_band.central,
            low=weight * paved_band.low + (1.0 - weight) * unpaved_band.low,
            high=weight * paved_band.high + (1.0 - weight) * unpaved_band.high,
            source=f"expected_surface:{weight:.3f} paved; {paved_band.source}",
            n_support=paved_band.n_support,
            fill_rule="expected_surface",
        )
        work = (
            work_type_for(road_class, lanes=lanes_used, surface="paved")
            if weight >= 0.5
            else work_type_for(road_class, lanes=lanes_used, surface="unpaved")
        )
        surface_label = f"unknown_expected_{weight:.2f}_paved"
    else:
        work = work_type_for(road_class, lanes=lanes_used, surface=surface_state)
        band, book_name = book.lookup_pavement(iso3, road_class, work)
        surface_label = surface_state

    factor = lane_factor(lanes_used, parent_default_lanes)
    scaled = _scale_band(band, factor * terrain_mult)
    return _accepted(
        road_class=road_class,
        is_link=is_link,
        is_bridge=False,
        is_tunnel=False,
        surface=surface_label,
        lanes_used=lanes_used,
        lanes_source=lanes_source,
        terrain_class=terrain_name,
        terrain_multiplier=terrain_mult,
        work_type=work,
        price_book=book_name,
        band=scaled,
        length_km=length,
    )


def _is_number(value: Any) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number == number


def _rejected(
    reason: str,
    *,
    length_km: float,
    road_class: str = "",
) -> ReplacementResult:
    return ReplacementResult(
        accepted=False,
        reason=reason,
        road_class=road_class,
        is_link=False,
        is_bridge=False,
        is_tunnel=False,
        surface="",
        lanes_used=None,
        lanes_source="",
        terrain_class="",
        terrain_multiplier=1.0,
        work_type="",
        price_book="",
        usd_per_km=0.0,
        usd_per_km_low=0.0,
        usd_per_km_high=0.0,
        length_km=length_km,
        replacement_usd=0.0,
        replacement_usd_low=0.0,
        replacement_usd_high=0.0,
    )


def _accepted(
    *,
    road_class: str,
    is_link: bool,
    is_bridge: bool,
    is_tunnel: bool,
    surface: str,
    lanes_used: Optional[float],
    lanes_source: str,
    terrain_class: str,
    terrain_multiplier: float,
    work_type: str,
    price_book: str,
    band: CostBand,
    length_km: float,
) -> ReplacementResult:
    return ReplacementResult(
        accepted=True,
        reason="ok",
        road_class=road_class,
        is_link=is_link,
        is_bridge=is_bridge,
        is_tunnel=is_tunnel,
        surface=surface,
        lanes_used=lanes_used,
        lanes_source=lanes_source,
        terrain_class=terrain_class,
        terrain_multiplier=terrain_multiplier,
        work_type=work_type,
        price_book=price_book,
        usd_per_km=band.central,
        usd_per_km_low=band.low,
        usd_per_km_high=band.high,
        length_km=length_km,
        replacement_usd=band.central * length_km,
        replacement_usd_low=band.low * length_km,
        replacement_usd_high=band.high * length_km,
    )

# --- country from lon/lat (Natural Earth 50m) ---
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

# --- OSM extract helpers ---
EXTRACT_HIGHWAY = {
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
        if highway not in EXTRACT_HIGHWAY:
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

# --- 2025 unit-cost book ---
NEW_WORK_TYPES = (WORK_TYPE_4L, WORK_TYPE_2L, WORK_TYPE_1L)

WB_TO_ROCKS_REGION = {
    "East Asia & Pacific": "East Asia and Pacific",
    "Europe & Central Asia": "Europe and Central Asia",
    "Latin America & Caribbean": "Latin America and Caribbean",
    "Middle East & North Africa": "Middle East and North Africa",
    "Middle East, North Africa, Afghanistan & Pakistan": "Middle East and North Africa",
    "North America": "North America",
    "South Asia": "South Asia",
    "Sub-Saharan Africa": "Sub-Saharan Africa",
}

# ROCKS 2018 coded these two as South Asia. The 2026 WDI region list moved
# them into the MENA group; keep the ROCKS geography for unit-cost lookup.
ROCKS_REGION_OVERRIDE = {
    "AFG": "South Asia",
    "PAK": "South Asia",
}

EUROPE_HIC_ISO3 = {
    "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
    "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
    "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE", "GBR",
    "CHE", "NOR", "ISL", "LIE", "AND", "MCO", "SMR",
}

NAMED_NATIONAL = {"USA", "CAN", "AUS", "NZL", "JPN", "KOR", "SGP", "TWN", "HKG", "CHN"}

EU28_2015 = [
    "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
    "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
    "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE", "GBR",
]

# Manual extras not always present / complete in the WDI country list.
MANUAL_COUNTRIES = {
    "TWN": {
        "name": "Taiwan",
        "wb_region": "East Asia & Pacific",
        "income_level": "HIC",
    },
    "XKX": {
        "name": "Kosovo",
        "wb_region": "Europe & Central Asia",
        "income_level": "UMC",
    },
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def usa_deflator_index(raw_defl: list) -> dict[int, float]:
    rows = raw_defl[1]
    out: dict[int, float] = {}
    for row in rows:
        if row.get("value") is None:
            continue
        out[int(row["date"])] = float(row["value"])
    return out


def inflate_usd(amount: float, year: int, to_year: int, deflator: dict[int, float]) -> float:
    if year not in deflator or to_year not in deflator:
        raise KeyError(f"deflator missing {year} or {to_year}")
    return amount * deflator[to_year] / deflator[year]


def latest_by_country(indicator_rows: list, *, year: int | None = None) -> dict[str, tuple[int, float]]:
    best: dict[str, tuple[int, float]] = {}
    for row in indicator_rows:
        iso = row.get("countryiso3code") or ""
        if len(iso) != 3 or row.get("value") is None:
            continue
        y = int(row["date"])
        if year is not None and y != year:
            continue
        value = float(row["value"])
        if iso not in best or y > best[iso][0]:
            best[iso] = (y, value)
    return best


def read_rocks_new_build(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            work = row["work_type"]
            if work not in NEW_WORK_TYPES:
                continue
            try:
                year = int(float(row["year"]))
                unit = float(row["unit_cost_million_usd_per_km"])
            except (TypeError, ValueError):
                continue
            if unit <= 0.0 or not math.isfinite(unit):
                continue
            rows.append(
                {
                    "country": row["country"],
                    "region": row["region"],
                    "year": year,
                    "work_type": work,
                    "unit_million": unit,
                }
            )
    return rows


def drop_outliers(values: list[float]) -> list[float]:
    """Drop 1-row absurdities: < 0.05 or > 30 million USD/km, or 5× the median."""

    cleaned = [v for v in values if 0.05 <= v <= 30.0]
    if len(cleaned) < 3:
        return cleaned
    med = statistics.median(cleaned)
    if med <= 0:
        return cleaned
    return [v for v in cleaned if v <= 5.0 * med]


def rocks_medians(
    rows: list[dict[str, Any]], deflator: dict[int, float], to_year: int
) -> dict[str, dict[str, CostBand]]:
    buckets: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        try:
            inflated = inflate_usd(row["unit_million"], row["year"], to_year, deflator)
        except KeyError:
            continue
        buckets[(row["region"], row["work_type"])].append(inflated)

    regional: dict[str, dict[str, CostBand]] = defaultdict(dict)
    global_pool: dict[str, list[float]] = defaultdict(list)
    for (region, work), values in buckets.items():
        kept = drop_outliers(values)
        if not kept:
            continue
        global_pool[work].extend(kept)
        regional[region][work] = _band_from_values(
            kept, source=f"ROCKS Actual {work} {region}, inflated to {to_year} USD"
        )

    global_band = {
        work: _band_from_values(
            drop_outliers(vals),
            source=f"global LMIC ROCKS Actual {work}, inflated to {to_year} USD",
        )
        for work, vals in global_pool.items()
        if drop_outliers(vals)
    }

    all_regions = sorted(
        set(WB_TO_ROCKS_REGION.values()) | {row["region"] for row in rows}
    )
    filled: dict[str, dict[str, CostBand]] = {}
    for region in all_regions:
        if region == "North America":
            continue
        filled[region] = {}
        for work in NEW_WORK_TYPES:
            if region in regional and work in regional[region] and regional[region][work].n_support >= 3:
                filled[region][work] = regional[region][work]
            elif region in regional and work in regional[region] and regional[region][work].n_support >= 2:
                band = regional[region][work]
                filled[region][work] = CostBand(
                    central=band.central,
                    low=band.low,
                    high=band.high,
                    source=band.source + " (n<3, keep regional)",
                    n_support=band.n_support,
                    fill_rule="regional_n2",
                )
            elif work in global_band:
                g = global_band[work]
                filled[region][work] = CostBand(
                    central=g.central,
                    low=g.low,
                    high=g.high,
                    source=g.source + f"; fill for {region}",
                    n_support=g.n_support,
                    fill_rule="global_lmic_median",
                )
    return filled


def _band_from_values(values: list[float], source: str) -> CostBand:
    values = sorted(values)
    med = statistics.median(values)
    if len(values) >= 4:
        low = statistics.quantiles(values, n=4)[0]
        high = statistics.quantiles(values, n=4)[2]
    else:
        low = min(values)
        high = max(values)
    # million USD -> USD
    return CostBand(
        central=med * 1.0e6,
        low=low * 1.0e6,
        high=high * 1.0e6,
        source=source,
        n_support=len(values),
        fill_rule="direct",
    )


def expand_national_anchors(anchors: dict[str, Any], to_usd_scale: float = 1.0) -> dict[str, dict[str, CostBand]]:
    raw = anchors["national_2025_usd_per_km"]
    out: dict[str, dict[str, CostBand]] = {}
    pending_copies: list[tuple[str, str, float]] = []
    for iso, spec in raw.items():
        if "copy_from" in spec:
            pending_copies.append((iso, spec["copy_from"], float(spec.get("scale", 1.0))))
            continue
        out[iso] = {}
        for road_class, item in spec.items():
            if road_class == "note":
                continue
            out[iso][road_class] = CostBand(
                central=float(item["central"]) * to_usd_scale,
                low=float(item["low"]) * to_usd_scale,
                high=float(item["high"]) * to_usd_scale,
                source=item.get("source", ""),
                n_support=1,
                fill_rule="national_anchor",
            )
    for iso, src, scale in pending_copies:
        if src not in out:
            raise KeyError(f"{iso} copies missing {src}")
        out[iso] = {
            road_class: CostBand(
                central=band.central * scale,
                low=band.low * scale,
                high=band.high * scale,
                source=f"scaled {scale} from {src}; {band.source}",
                n_support=band.n_support,
                fill_rule="national_scaled",
            )
            for road_class, band in out[src].items()
        }
    return out


def europe_baseline_usd(anchors: dict[str, Any], deflator: dict[int, float]) -> dict[str, CostBand]:
    factor = anchors["eurusd_2015"] * (deflator[2025] / deflator[2015])
    out = {}
    for road_class, item in anchors["europe_2015_eur_per_km"].items():
        out[road_class] = CostBand(
            central=float(item["central"]) * factor,
            low=float(item["low"]) * factor,
            high=float(item["high"]) * factor,
            source=item.get("source", "") + f"; ×{factor:.4f} to 2025 USD",
            n_support=1,
            fill_rule="europe_2015_converted",
        )
    return out


def income_code(value: str) -> str:
    mapping = {
        "High income": "HIC",
        "Upper middle income": "UMC",
        "Lower middle income": "LMC",
        "Low income": "LIC",
    }
    return mapping.get(value, "UMC")


def choose_price_book(iso3: str, income: str, region: str) -> str:
    if iso3 in NAMED_NATIONAL:
        return PRICE_BOOK_NATIONAL
    if iso3 in EUROPE_HIC_ISO3 and income == "HIC":
        return PRICE_BOOK_EUROPE
    # EU members that WDI still lists as UMC (none expected) stay ROCKS unless named.
    if income == "HIC":
        return PRICE_BOOK_HIC_SCALE
    return PRICE_BOOK_ROCKS


def build_countries(
    wb_countries: list[dict[str, Any]],
    ppp_2015: dict[str, tuple[int, float]],
    ppp_latest: dict[str, tuple[int, float]],
    anchors: dict[str, Any],
) -> dict[str, CountryRecord]:
    paved_default = anchors["income_group_paved_default"]
    usa_ppp = ppp_latest.get("USA", (None, None))[1]
    eu_vals = [ppp_2015[c][1] for c in EU28_2015 if c in ppp_2015]
    eu_mean = sum(eu_vals) / len(eu_vals)

    records: dict[str, CountryRecord] = {}
    source_rows = []
    for row in wb_countries:
        region_name = (row.get("region") or {}).get("value")
        if region_name in (None, "Aggregates"):
            continue
        iso = row["id"]
        if len(iso) != 3:
            continue
        source_rows.append(
            {
                "id": iso,
                "name": row.get("name", iso),
                "region": region_name,
                "income": (row.get("incomeLevel") or {}).get("value", ""),
            }
        )
    for iso, extra in MANUAL_COUNTRIES.items():
        if iso not in {r["id"] for r in source_rows}:
            source_rows.append(
                {
                    "id": iso,
                    "name": extra["name"],
                    "region": extra["wb_region"],
                    "income": extra["income_level"]
                    if extra["income_level"] in {"High income", "Upper middle income", "Lower middle income", "Low income"}
                    else {
                        "HIC": "High income",
                        "UMC": "Upper middle income",
                        "LMC": "Lower middle income",
                        "LIC": "Low income",
                    }[extra["income_level"]],
                }
            )

    for row in source_rows:
        iso = row["id"]
        income = income_code(row["income"])
        wb_region = (row["region"] or "").strip()
        rocks_region = ROCKS_REGION_OVERRIDE.get(
            iso, WB_TO_ROCKS_REGION.get(wb_region, "")
        )
        book = choose_price_book(iso, income, wb_region)
        paved = float(paved_default[income])
        gdp_2015 = ppp_2015[iso][1] if iso in ppp_2015 else None
        gdp_latest = ppp_latest[iso][1] if iso in ppp_latest else None
        eu_ratio = (gdp_2015 / eu_mean) if gdp_2015 else None
        usa_ratio = (gdp_latest / usa_ppp) if (gdp_latest and usa_ppp) else None
        records[iso] = CountryRecord(
            iso3=iso,
            name=row["name"],
            wb_region=wb_region,
            rocks_region=rocks_region,
            income_level=income,
            price_book=book,
            paved_share=paved,
            paved_share_source=f"income_group_default:{income}",
            gdp_pc_ppp_2015=gdp_2015,
            gdp_pc_ppp_latest=gdp_latest,
            eu28_ratio_2015=eu_ratio,
            usa_ratio_latest=usa_ratio,
        )
    return records


def assemble_book(
    repo: Path,
    *,
    construction_premium_2018_2025: float = 1.0,
) -> tuple[UnitCostBook, dict[str, Any]]:
    raw = repo / "data" / "valuation" / "raw"
    anchors = load_json(repo / "data" / "valuation" / "national_anchors.json")
    deflator = usa_deflator_index(load_json(raw / "wb_defl_usa.json"))
    if construction_premium_2018_2025 != 1.0:
        # Stretch only the 2018-2025 increment. Years before 2018 keep GDP path
        # to 2018, then apply premium × remaining GDP path.
        base_2018 = deflator[2018]
        extra = construction_premium_2018_2025
        for year in list(deflator):
            if year > 2018:
                deflator[year] = base_2018 + extra * (deflator[year] - base_2018)

    wb_countries = load_json(raw / "wb_countries.json")[1]
    ppp_rows = load_json(raw / "wb_gdp_ppp.json")[1]
    ppp_2015 = latest_by_country(ppp_rows, year=2015)
    ppp_latest = latest_by_country(ppp_rows)

    rocks_rows = read_rocks_new_build(
        repo / "data" / "rocks" / "rocks_2018_actual_construction_like_rows.csv"
    )
    rocks = rocks_medians(rocks_rows, deflator, 2025)
    national = expand_national_anchors(anchors)
    europe = europe_baseline_usd(anchors, deflator)
    usa = national["USA"]
    countries = build_countries(wb_countries, ppp_2015, ppp_latest, anchors)

    d2023 = deflator[2023]
    d2025 = deflator[2025]
    structure_scale = d2025 / d2023
    bridge = CostBand(
        central=anchors["giri_bridge_usd_per_km_2023"] * structure_scale,
        low=anchors["giri_bridge_usd_per_km_2023"] * structure_scale * 0.5,
        high=anchors["giri_bridge_usd_per_km_2023"] * structure_scale * 2.0,
        source=anchors["giri_source"] + "; inflated 2023→2025 US GDP deflator",
        n_support=1,
        fill_rule="giri",
    )
    tunnel = CostBand(
        central=anchors["giri_tunnel_usd_per_km_2023"] * structure_scale,
        low=anchors["giri_tunnel_usd_per_km_2023"] * structure_scale * 0.5,
        high=anchors["giri_tunnel_usd_per_km_2023"] * structure_scale * 2.0,
        source=anchors["giri_source"] + "; inflated 2023→2025 US GDP deflator",
        n_support=1,
        fill_rule="giri",
    )
    clip = tuple(anchors["hic_gdp_scale_clip"])
    book = UnitCostBook(
        countries=countries,
        rocks=rocks,
        national=national,
        europe_baseline_usd=europe,
        usa_baseline=usa,
        bridge=bridge,
        tunnel=tunnel,
        hic_clip=clip,
    )
    meta = {
        "script_version": SCRIPT_VERSION,
        "price_year": 2025,
        "construction_premium_2018_2025": construction_premium_2018_2025,
        "usa_deflator_2018": deflator[2018],
        "usa_deflator_2025": deflator[2025],
        "usa_deflator_factor_2018_2025": deflator[2025] / deflator[2018],
        "n_countries": len(countries),
        "n_rocks_new_build_rows": len(rocks_rows),
        "rocks_regions": sorted(rocks),
        "named_national": sorted(national),
        "eu28_2015_ppp_mean": sum(ppp_2015[c][1] for c in EU28_2015 if c in ppp_2015)
        / sum(1 for c in EU28_2015 if c in ppp_2015),
    }
    return book, meta


def write_ledger(book: UnitCostBook, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "iso3",
        "name",
        "wb_region",
        "rocks_region",
        "income_level",
        "price_book",
        "paved_share",
        "paved_share_source",
        "gdp_pc_ppp_2015",
        "gdp_pc_ppp_latest",
        "eu28_ratio_2015",
        "usa_ratio_latest",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for iso in sorted(book.countries):
            rec = book.countries[iso]
            writer.writerow({name: getattr(rec, name) for name in fields})


def write_rocks_table(book: UnitCostBook, path: Path) -> None:
    fields = [
        "region",
        "work_type",
        "n_support",
        "fill_rule",
        "usd_per_km_2025",
        "usd_per_km_low",
        "usd_per_km_high",
        "source",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for region in sorted(book.rocks):
            for work in NEW_WORK_TYPES:
                band = book.rocks[region][work]
                writer.writerow(
                    {
                        "region": region,
                        "work_type": work,
                        "n_support": band.n_support,
                        "fill_rule": band.fill_rule,
                        "usd_per_km_2025": f"{band.central:.2f}",
                        "usd_per_km_low": f"{band.low:.2f}",
                        "usd_per_km_high": f"{band.high:.2f}",
                        "source": band.source,
                    }
                )

# --- price a table of roads ---
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
    rows = assign_rows(rows)
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
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_book = sub.add_parser("write-book", help="Write the 2025 country unit-cost tables")
    p_book.add_argument("--repo", type=Path, default=ROOT)
    p_book.add_argument("--construction-premium", type=float, default=1.0)

    p_val = sub.add_parser("value", help="Price a CSV of road segments")
    p_val.add_argument("input_csv", type=Path)
    p_val.add_argument("--output", type=Path, required=True)
    p_val.add_argument("--repo", type=Path, default=ROOT)

    p_ex = sub.add_parser("extract", help="Extract motor roads from an OSM PBF")
    p_ex.add_argument("pbf", type=Path)
    p_ex.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    if args.cmd == "write-book":
        book, meta = assemble_book(
            args.repo, construction_premium_2018_2025=args.construction_premium
        )
        out = args.repo / "data" / "valuation"
        write_ledger(book, out / "country_ledger.csv")
        write_rocks_table(book, out / "rocks_unit_cost_2025.csv")
        (out / "unit_cost_book.manifest.json").write_text(
            json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(meta, indent=2, sort_keys=True))
        return 0

    if args.cmd == "value":
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
        printable = {k: totals[k] for k in totals if k != "countries_with_lane_median"}
        print(json.dumps(printable, indent=2))
        return 0

    if args.cmd == "extract":
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
            "accepted_highway": sorted(EXTRACT_HIGHWAY),
            "stats": stats,
        }
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    raise SystemExit("unknown command")


if __name__ == "__main__":
    raise SystemExit(main())
