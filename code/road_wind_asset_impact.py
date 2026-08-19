#!/usr/bin/env python3
"""Object-level tropical-cyclone wind asset loss.

Follows Koks et al. (2019) using the public gmtra operational rules.
Cleanup dollars use Escobedo et al. (2009) priced 2025 USD, not the Koks
SI assumed cleanup bands.  Pavement replacement cost is never multiplied
by a wind MDR.

Usage:
  python3 road_wind_asset_impact.py constants
  python3 road_wind_asset_impact.py apply valued.csv --output impact.csv
  python3 road_wind_asset_impact.py hash-trees crowther.tif
  python3 road_wind_asset_impact.py sample-trees valued.csv --trees crowther.tif --output valued.trees.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

from road_replacement_value import assemble_book, classify_highway


SCRIPT_VERSION = "1.0.0"
CONTRACT_NAME = "WIND_ASSET_IMPACT_CONTRACT.md"

G_INLAND_3S_10MIN = 1.66
G_OFFLAND_3S_10MIN = 1.52
G_INLAND_3S_1MIN = 1.49

TREE_BREAK_GUST_KMH = 151.0
MS_PER_KMH = 1.0 / 3.6

ESCOBEDO_VOL_M3_PER_30P5M = {"low": 0.59, "central": 3.40, "high": 17.47}
ESCOBEDO_USD_PER_M3_2005 = 28.25
DEFLATOR_USA_2005 = 77.3945446986751
DEFLATOR_USA_2025 = 122.361624831936
STREET_SEGMENT_M = 30.5
TREE_DENSITY_CAP_KM2 = 10000.0
HISTORICAL_WINDOW_YEARS = 20.0

# gmtra parallel.py wind_threshs row mid-points (3-s gust, km/h).
BRIDGE_GUST_VSTAR_KMH = {
    "primary": 362.5,
    "secondary": 337.5,
    "other": 312.5,
}

# gmtra parallel.py design_tables, years.
DESIGN_RP_YEARS = {
    "HIC": {"primary": 200.0, "secondary": 100.0, "other": 50.0},
    "UMC": {"primary": 100.0, "secondary": 50.0, "other": 20.0},
    "LMC": {"primary": 50.0, "secondary": 20.0, "other": 10.0},
}

CROWTHER_SOURCE = (
    "Crowther et al. (2015) Nature 525:201-205; Yale EliScholar "
    "yale_fes_data/1 biome WGS84 GeoTIFF Revision_01"
)


def usa_deflator_2005_to_2025() -> float:
    return DEFLATOR_USA_2025 / DEFLATOR_USA_2005


def escobedo_cleanup_usd_per_km(band: str = "central") -> float:
    if band not in ESCOBEDO_VOL_M3_PER_30P5M:
        raise ValueError(f"unknown Escobedo band: {band!r}")
    volume_per_km = ESCOBEDO_VOL_M3_PER_30P5M[band] * (1000.0 / STREET_SEGMENT_M)
    usd_2005 = volume_per_km * ESCOBEDO_USD_PER_M3_2005
    return usd_2005 * usa_deflator_2005_to_2025()


def c15_threshold_ms(gust_kmh: float, gust_factor: float = G_INLAND_3S_10MIN) -> float:
    if gust_factor <= 0.0:
        raise ValueError("gust factor must be positive")
    return (float(gust_kmh) / float(gust_factor)) * MS_PER_KMH


def tree_break_c15_ms(gust_factor: float = G_INLAND_3S_10MIN) -> float:
    return c15_threshold_ms(TREE_BREAK_GUST_KMH, gust_factor)


def tree_fail_prob(tree_density_km2: float) -> float:
    """gmtra regional_cyclone: drop non-positive density, cap 10000, /10000."""

    density = float(tree_density_km2)
    if not math.isfinite(density) or density <= 0.0:
        return 0.0
    if density >= TREE_DENSITY_CAP_KM2:
        return 1.0
    return density / TREE_DENSITY_CAP_KM2


def koks_asset_class(highway: str) -> Optional[str]:
    road_class, _is_link, _reason = classify_highway(highway)
    if road_class is None:
        return None
    if road_class in {"motorway", "trunk", "primary"}:
        return "primary"
    if road_class == "secondary":
        return "secondary"
    return "other"


def gmtra_income_group(income_level: str) -> str:
    if income_level == "HIC":
        return "HIC"
    if income_level == "UMC":
        return "UMC"
    return "LMC"


def design_rp_years(income_level: str, asset_class: str) -> float:
    table = DESIGN_RP_YEARS[gmtra_income_group(income_level)]
    if asset_class not in table:
        raise ValueError(f"unknown gmtra class: {asset_class!r}")
    return table[asset_class]


def bridge_vstar_c15_ms(
    asset_class: str, gust_factor: float = G_INLAND_3S_10MIN
) -> float:
    if asset_class not in BRIDGE_GUST_VSTAR_KMH:
        raise ValueError(f"unknown gmtra class: {asset_class!r}")
    return c15_threshold_ms(BRIDGE_GUST_VSTAR_KMH[asset_class], gust_factor)


def empirical_rp_years(
    event_ms: float,
    historical_peaks_ms: Sequence[float],
    window_years: float = HISTORICAL_WINDOW_YEARS,
) -> float:
    """Return period of *event_ms* from catalogue peaks at the same site.

    RP = window / M, M = count of peaks >= event.  If M = 0 the event is
    rarer than the catalogue and RP is infinite.
    """

    if window_years <= 0.0:
        raise ValueError("window_years must be positive")
    if not math.isfinite(event_ms) or event_ms < 0.0:
        raise ValueError("event wind must be a finite non-negative speed")
    count = 0
    for peak in historical_peaks_ms:
        value = float(peak)
        if math.isfinite(value) and value >= event_ms:
            count += 1
    if count == 0:
        return math.inf
    return window_years / count


def cleanup_usd(
    *,
    length_km: float,
    tree_density_km2: float,
    v_c15_ms: float,
    is_tunnel: bool,
    is_bridge: bool,
    accepted: bool,
    band: str = "central",
    gust_factor: float = G_INLAND_3S_10MIN,
) -> float:
    if not accepted or is_tunnel or is_bridge:
        return 0.0
    if not math.isfinite(v_c15_ms) or v_c15_ms < tree_break_c15_ms(gust_factor):
        return 0.0
    probability = tree_fail_prob(tree_density_km2)
    if probability <= 0.0:
        return 0.0
    return float(length_km) * escobedo_cleanup_usd_per_km(band) * probability


def bridge_collapse_usd(
    *,
    replacement_usd: float,
    is_bridge: bool,
    is_tunnel: bool,
    accepted: bool,
    highway: str,
    income_level: str,
    v_c15_ms: float,
    historical_peaks_ms: Sequence[float],
    window_years: float = HISTORICAL_WINDOW_YEARS,
    gust_factor: float = G_INLAND_3S_10MIN,
) -> float:
    if not accepted or is_tunnel or not is_bridge:
        return 0.0
    asset_class = koks_asset_class(highway)
    if asset_class is None:
        return 0.0
    if not math.isfinite(v_c15_ms) or v_c15_ms <= bridge_vstar_c15_ms(
        asset_class, gust_factor
    ):
        return 0.0
    rp = empirical_rp_years(v_c15_ms, historical_peaks_ms, window_years)
    if rp <= design_rp_years(income_level, asset_class):
        return 0.0
    return float(replacement_usd)


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in {"", "0", "false", "no", "nan"}:
        return False
    return text in {"1", "true", "yes"} or bool(text)


def _as_float(value: Any, default: float = float("nan")) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_peak_list(value: Any) -> list[float]:
    if value is None or value == "":
        return []
    if isinstance(value, (list, tuple)):
        return [float(item) for item in value]
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        return [float(item) for item in json.loads(text)]
    return [float(piece) for piece in text.replace(",", ";").split(";") if piece.strip()]


def impact_row(
    row: dict[str, Any],
    *,
    income_level: str,
    band: str = "central",
    gust_factor: float = G_INLAND_3S_10MIN,
    window_years: float = HISTORICAL_WINDOW_YEARS,
) -> dict[str, Any]:
    accepted = _as_bool(row.get("accepted", True))
    is_bridge = _as_bool(row.get("is_bridge", row.get("bridge")))
    is_tunnel = _as_bool(row.get("is_tunnel", row.get("tunnel")))
    highway = str(row.get("highway", ""))
    length_km = _as_float(row.get("length_km"), 0.0)
    replacement = _as_float(row.get("replacement_usd"), 0.0)
    tree_density = _as_float(row.get("tree_dens_km2"))
    wind_ms = _as_float(row.get("v_c15_ms"))
    peaks = parse_peak_list(row.get("bridge_hist_peaks_ms"))

    cleanup = cleanup_usd(
        length_km=length_km,
        tree_density_km2=tree_density,
        v_c15_ms=wind_ms,
        is_tunnel=is_tunnel,
        is_bridge=is_bridge,
        accepted=accepted,
        band=band,
        gust_factor=gust_factor,
    )
    collapse = bridge_collapse_usd(
        replacement_usd=replacement,
        is_bridge=is_bridge,
        is_tunnel=is_tunnel,
        accepted=accepted,
        highway=highway,
        income_level=income_level,
        v_c15_ms=wind_ms,
        historical_peaks_ms=peaks,
        window_years=window_years,
        gust_factor=gust_factor,
    )
    out = dict(row)
    out["tree_fail_prob"] = tree_fail_prob(tree_density)
    out["wind_cleanup_usd"] = cleanup
    out["wind_bridge_usd"] = collapse
    out["wind_asset_usd"] = cleanup + collapse
    out["koks_asset_class"] = koks_asset_class(highway) or ""
    out["gmtra_income_group"] = gmtra_income_group(income_level)
    return out


def load_income_levels(repo: Path) -> dict[str, str]:
    book, _meta = assemble_book(repo)
    return {
        iso: record.income_level
        for iso, record in book.countries.items()
        if record.income_level
    }


def apply_rows(
    rows: Iterable[dict[str, Any]],
    repo: Path,
    *,
    band: str = "central",
    gust_factor: float = G_INLAND_3S_10MIN,
    window_years: float = HISTORICAL_WINDOW_YEARS,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    incomes = load_income_levels(repo)
    impacted: list[dict[str, Any]] = []
    totals = {
        "ways": 0,
        "cleanup_usd": 0.0,
        "bridge_usd": 0.0,
        "wind_asset_usd": 0.0,
        "cleanup_ways": 0,
        "collapsed_bridges": 0,
    }
    for row in rows:
        iso3 = str(row.get("iso3", "")).upper()
        income = incomes.get(iso3, "UMC")
        record = impact_row(
            row,
            income_level=income,
            band=band,
            gust_factor=gust_factor,
            window_years=window_years,
        )
        impacted.append(record)
        totals["ways"] += 1
        totals["cleanup_usd"] += record["wind_cleanup_usd"]
        totals["bridge_usd"] += record["wind_bridge_usd"]
        totals["wind_asset_usd"] += record["wind_asset_usd"]
        if record["wind_cleanup_usd"] > 0.0:
            totals["cleanup_ways"] += 1
        if record["wind_bridge_usd"] > 0.0:
            totals["collapsed_bridges"] += 1
    totals["gust_factor"] = gust_factor
    totals["escobedo_band"] = band
    totals["tree_break_c15_ms"] = tree_break_c15_ms(gust_factor)
    totals["cleanup_usd_per_km"] = escobedo_cleanup_usd_per_km(band)
    totals["script_version"] = SCRIPT_VERSION
    totals["contract"] = CONTRACT_NAME
    return impacted, totals


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sample_open_dataset(dataset: Any, lon: float, lat: float) -> float:
    """Nearest pixel on an open GDAL dataset. NaN if outside or nodata."""

    from osgeo import gdal

    if not math.isfinite(lon) or not math.isfinite(lat):
        return float("nan")
    transform = dataset.GetGeoTransform()
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    inverse = gdal.InvGeoTransform(transform)
    if inverse is None:
        raise ValueError("Crowther GeoTIFF geotransform is not invertible")
    px, py = gdal.ApplyGeoTransform(inverse, float(lon), float(lat))
    col = int(math.floor(px))
    row = int(math.floor(py))
    if col < 0 or row < 0 or col >= dataset.RasterXSize or row >= dataset.RasterYSize:
        return float("nan")
    value = float(band.ReadAsArray(col, row, 1, 1)[0, 0])
    if nodata is not None and value == float(nodata):
        return float("nan")
    return value


def sample_crowther_density(
    geotiff: Path, lon: float, lat: float
) -> float:
    """Nearest Crowther pixel at lon/lat. NaN if outside the raster."""

    from osgeo import gdal

    dataset = gdal.Open(str(geotiff), gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(geotiff)
    try:
        return _sample_open_dataset(dataset, lon, lat)
    finally:
        dataset = None


def attach_crowther_density(
    rows: Iterable[dict[str, Any]], geotiff: Path
) -> list[dict[str, Any]]:
    """Write ``tree_dens_km2`` from Crowther at each row's lon/lat."""

    from osgeo import gdal

    dataset = gdal.Open(str(geotiff), gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(geotiff)
    try:
        attached: list[dict[str, Any]] = []
        for row in rows:
            record = dict(row)
            record["tree_dens_km2"] = _sample_open_dataset(
                dataset,
                _as_float(row.get("lon")),
                _as_float(row.get("lat")),
            )
            attached.append(record)
        return attached
    finally:
        dataset = None


def frozen_constants() -> dict[str, Any]:
    return {
        "g_inland_3s_10min": G_INLAND_3S_10MIN,
        "tree_break_gust_kmh": TREE_BREAK_GUST_KMH,
        "tree_break_c15_kmh": TREE_BREAK_GUST_KMH / G_INLAND_3S_10MIN,
        "tree_break_c15_ms": tree_break_c15_ms(),
        "cleanup_usd_per_km": {
            band: escobedo_cleanup_usd_per_km(band)
            for band in ("low", "central", "high")
        },
        "cleanup_usd_per_km_rounded": {
            band: int(round(escobedo_cleanup_usd_per_km(band)))
            for band in ("low", "central", "high")
        },
        "deflator_2005_to_2025": usa_deflator_2005_to_2025(),
        "tree_density_cap_km2": TREE_DENSITY_CAP_KM2,
        "historical_window_years": HISTORICAL_WINDOW_YEARS,
        "bridge_gust_vstar_kmh": dict(BRIDGE_GUST_VSTAR_KMH),
        "design_rp_years": DESIGN_RP_YEARS,
        "crowther_source": CROWTHER_SOURCE,
        "si_gust_factors": {
            "off_land_3s_10min": G_OFFLAND_3S_10MIN,
            "inland_3s_1min": G_INLAND_3S_1MIN,
        },
        "script_version": SCRIPT_VERSION,
        "contract": CONTRACT_NAME,
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("no rows to write")
    fieldnames: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("constants", help="print frozen contract numbers as JSON")

    apply_p = sub.add_parser("apply", help="score a valued-way table that already has wind and trees")
    apply_p.add_argument("input_csv", type=Path)
    apply_p.add_argument("--output", type=Path, required=True)
    apply_p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    apply_p.add_argument("--band", choices=("low", "central", "high"), default="central")
    apply_p.add_argument("--gust-factor", type=float, default=G_INLAND_3S_10MIN)
    apply_p.add_argument("--window-years", type=float, default=HISTORICAL_WINDOW_YEARS)

    hash_p = sub.add_parser("hash-trees", help="SHA-256 a Crowther GeoTIFF into a manifest")
    hash_p.add_argument("geotiff", type=Path)
    hash_p.add_argument("--output", type=Path, required=True)

    sample_p = sub.add_parser(
        "sample-trees",
        help="write tree_dens_km2 from a Crowther GeoTIFF onto a valued-way table",
    )
    sample_p.add_argument("input_csv", type=Path)
    sample_p.add_argument("--trees", type=Path, required=True)
    sample_p.add_argument("--output", type=Path, required=True)

    args = parser.parse_args(argv)

    if args.command == "constants":
        print(json.dumps(frozen_constants(), indent=2, sort_keys=True))
        return 0

    if args.command == "hash-trees":
        geotiff = args.geotiff.resolve()
        if not geotiff.is_file():
            raise FileNotFoundError(geotiff)
        payload = {
            "path": str(geotiff),
            "sha256": sha256_file(geotiff),
            "bytes": geotiff.stat().st_size,
            "source": CROWTHER_SOURCE,
            "script_version": SCRIPT_VERSION,
            "contract": CONTRACT_NAME,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    if args.command == "sample-trees":
        geotiff = args.trees.resolve()
        if not geotiff.is_file():
            raise FileNotFoundError(geotiff)
        attached = attach_crowther_density(_read_csv(args.input_csv), geotiff)
        _write_csv(args.output, attached)
        finite = [
            float(row["tree_dens_km2"])
            for row in attached
            if math.isfinite(float(row["tree_dens_km2"]))
        ]
        summary = {
            "ways": len(attached),
            "sampled_finite": len(finite),
            "sampled_nodata_or_outside": len(attached) - len(finite),
            "trees": str(geotiff),
            "script_version": SCRIPT_VERSION,
            "contract": CONTRACT_NAME,
        }
        summary_path = Path(str(args.output) + ".summary.json")
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    rows = _read_csv(args.input_csv)
    impacted, totals = apply_rows(
        rows,
        args.repo,
        band=args.band,
        gust_factor=args.gust_factor,
        window_years=args.window_years,
    )
    _write_csv(args.output, impacted)
    summary_path = Path(str(args.output) + ".summary.json")
    summary_path.write_text(json.dumps(totals, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(totals, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
