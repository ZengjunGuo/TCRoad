#!/usr/bin/env python3
"""Build the frozen 2025 country unit-cost book from ROCKS, WDI and anchors."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
import sys
from typing import Any


CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from road_replacement_value import (  # noqa: E402
    CostBand,
    CountryRecord,
    PRICE_BOOK_EUROPE,
    PRICE_BOOK_HIC_SCALE,
    PRICE_BOOK_NATIONAL,
    PRICE_BOOK_ROCKS,
    SCRIPT_VERSION,
    UnitCostBook,
    WORK_TYPE_1L,
    WORK_TYPE_2L,
    WORK_TYPE_4L,
)


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument(
        "--construction-premium",
        type=float,
        default=1.0,
        help="Extra 2018-2025 construction residual on the US GDP deflator (1.0 = main case).",
    )
    args = parser.parse_args()
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


if __name__ == "__main__":
    raise SystemExit(main())
