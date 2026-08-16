#!/usr/bin/env python3
"""Object-level 2025 replacement cost for one motor-road segment.

This module is the valuation kernel.  It does not read OSM files and does not
talk to the network.  `build_road_unit_cost_book.py` builds the country book;
this file applies the frozen contract in
`methods/ROAD_ASSET_VALUATION_CONTRACT.md`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional


SCRIPT_VERSION = "1.0.0"

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
