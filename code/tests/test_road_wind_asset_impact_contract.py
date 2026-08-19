"""Contract tests for Koks / gmtra wind asset loss."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


CODE_DIR = Path(__file__).resolve().parents[1]
REPO = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from road_wind_asset_impact import (  # noqa: E402
    BRIDGE_GUST_VSTAR_KMH,
    G_INLAND_3S_10MIN,
    TREE_BREAK_GUST_KMH,
    apply_rows,
    attach_crowther_density,
    bridge_collapse_usd,
    bridge_vstar_c15_ms,
    c15_threshold_ms,
    cleanup_usd,
    design_rp_years,
    empirical_rp_years,
    escobedo_cleanup_usd_per_km,
    frozen_constants,
    gmtra_income_group,
    koks_asset_class,
    main,
    sample_crowther_density,
    tree_break_c15_ms,
    tree_fail_prob,
)


class WindAssetImpactTests(unittest.TestCase):
    def test_escobedo_central_rounds_to_contract_4979(self) -> None:
        self.assertEqual(int(round(escobedo_cleanup_usd_per_km("central"))), 4979)
        self.assertEqual(int(round(escobedo_cleanup_usd_per_km("low"))), 864)
        self.assertEqual(int(round(escobedo_cleanup_usd_per_km("high"))), 25583)

    def test_tree_break_threshold_uses_inland_1_66(self) -> None:
        expected_kmh = TREE_BREAK_GUST_KMH / G_INLAND_3S_10MIN
        self.assertAlmostEqual(expected_kmh, 91.0, places=1)
        self.assertAlmostEqual(tree_break_c15_ms(), expected_kmh / 3.6, places=6)
        self.assertGreater(tree_break_c15_ms(1.49), tree_break_c15_ms(1.66))

    def test_gmtra_density_factor(self) -> None:
        self.assertEqual(tree_fail_prob(0.0), 0.0)
        self.assertEqual(tree_fail_prob(-1.0), 0.0)
        self.assertEqual(tree_fail_prob(float("nan")), 0.0)
        self.assertAlmostEqual(tree_fail_prob(10.0), 0.001)
        self.assertAlmostEqual(tree_fail_prob(1000.0), 0.1)
        self.assertEqual(tree_fail_prob(10000.0), 1.0)
        self.assertEqual(tree_fail_prob(50000.0), 1.0)

    def test_koks_class_maps_motorway_trunk_to_primary(self) -> None:
        self.assertEqual(koks_asset_class("motorway"), "primary")
        self.assertEqual(koks_asset_class("trunk"), "primary")
        self.assertEqual(koks_asset_class("primary"), "primary")
        self.assertEqual(koks_asset_class("secondary"), "secondary")
        self.assertEqual(koks_asset_class("tertiary"), "other")
        self.assertEqual(koks_asset_class("residential"), "other")
        self.assertIsNone(koks_asset_class("footway"))

    def test_design_rp_follows_gmtra_tables(self) -> None:
        self.assertEqual(design_rp_years("HIC", "primary"), 200.0)
        self.assertEqual(design_rp_years("UMC", "secondary"), 50.0)
        self.assertEqual(design_rp_years("LIC", "other"), 10.0)
        self.assertEqual(gmtra_income_group("LIC"), "LMC")

    def test_cleanup_zero_below_threshold_or_on_structures(self) -> None:
        wind = tree_break_c15_ms() + 1.0
        kwargs = dict(
            length_km=2.0,
            tree_density_km2=10000.0,
            v_c15_ms=wind,
            is_tunnel=False,
            is_bridge=False,
            accepted=True,
        )
        full = cleanup_usd(**kwargs)
        self.assertAlmostEqual(full, 2.0 * escobedo_cleanup_usd_per_km("central"))
        self.assertEqual(cleanup_usd(**{**kwargs, "v_c15_ms": tree_break_c15_ms() - 0.01}), 0.0)
        self.assertEqual(cleanup_usd(**{**kwargs, "is_tunnel": True}), 0.0)
        self.assertEqual(cleanup_usd(**{**kwargs, "is_bridge": True}), 0.0)
        self.assertEqual(cleanup_usd(**{**kwargs, "tree_density_km2": 0.0}), 0.0)
        self.assertEqual(cleanup_usd(**{**kwargs, "accepted": False}), 0.0)
        half = cleanup_usd(**{**kwargs, "tree_density_km2": 5000.0})
        self.assertAlmostEqual(half, 0.5 * full)

    def test_cleanup_does_not_use_replacement_cost(self) -> None:
        wind = tree_break_c15_ms() + 10.0
        loss = cleanup_usd(
            length_km=1.0,
            tree_density_km2=10000.0,
            v_c15_ms=wind,
            is_tunnel=False,
            is_bridge=False,
            accepted=True,
        )
        self.assertLess(loss, 20000.0)
        self.assertNotAlmostEqual(loss, 12_000_000.0)

    def test_bridge_collapses_only_when_gust_and_rp_both_clear(self) -> None:
        vstar = bridge_vstar_c15_ms("primary")
        replacement = 4_146_699.67
        base = dict(
            replacement_usd=replacement,
            is_bridge=True,
            is_tunnel=False,
            accepted=True,
            highway="primary",
            income_level="HIC",
            window_years=20.0,
        )
        # Unique extreme gust: rarer than 1/200.
        self.assertEqual(
            bridge_collapse_usd(
                **base, v_c15_ms=vstar + 1.0, historical_peaks_ms=[]
            ),
            replacement,
        )
        # Frequent at this site: 20 events in 20 years → RP = 1 year.
        frequent = [vstar + 5.0] * 20
        self.assertEqual(
            bridge_collapse_usd(
                **base, v_c15_ms=vstar + 5.0, historical_peaks_ms=frequent
            ),
            0.0,
        )
        self.assertEqual(
            bridge_collapse_usd(
                **base, v_c15_ms=vstar - 1.0, historical_peaks_ms=[]
            ),
            0.0,
        )
        pavement = dict(base)
        pavement["is_bridge"] = False
        self.assertEqual(
            bridge_collapse_usd(
                **pavement, v_c15_ms=vstar + 1.0, historical_peaks_ms=[]
            ),
            0.0,
        )

    def test_empirical_rp_infinite_when_unseen(self) -> None:
        self.assertTrue(math.isinf(empirical_rp_years(50.0, [])))
        self.assertAlmostEqual(empirical_rp_years(50.0, [50.0, 40.0], 20.0), 20.0)
        self.assertAlmostEqual(empirical_rp_years(40.0, [50.0, 40.0], 20.0), 10.0)

    def test_apply_fixture_table(self) -> None:
        wind = tree_break_c15_ms() + 2.0
        rows = [
            {
                "way_id": "US-I10",
                "iso3": "USA",
                "highway": "motorway",
                "accepted": 1,
                "is_bridge": 0,
                "is_tunnel": 0,
                "length_km": 2.0,
                "replacement_usd": 24_000_000.0,
                "tree_dens_km2": 10000,
                "v_c15_ms": wind,
            },
            {
                "way_id": "JP-bridge",
                "iso3": "JPN",
                "highway": "primary",
                "accepted": 1,
                "is_bridge": 1,
                "is_tunnel": 0,
                "length_km": 0.4,
                "replacement_usd": 4_146_699.67,
                "tree_dens_km2": 8000,
                "v_c15_ms": bridge_vstar_c15_ms("primary") + 1.0,
                "bridge_hist_peaks_ms": "",
            },
            {
                "way_id": "FOOT",
                "iso3": "USA",
                "highway": "footway",
                "accepted": 0,
                "is_bridge": 0,
                "is_tunnel": 0,
                "length_km": 1.0,
                "replacement_usd": 0.0,
                "tree_dens_km2": 10000,
                "v_c15_ms": wind,
            },
        ]
        impacted, totals = apply_rows(rows, REPO)
        by_id = {row["way_id"]: row for row in impacted}
        self.assertAlmostEqual(
            by_id["US-I10"]["wind_cleanup_usd"],
            2.0 * escobedo_cleanup_usd_per_km("central"),
        )
        self.assertEqual(by_id["US-I10"]["wind_bridge_usd"], 0.0)
        self.assertEqual(by_id["JP-bridge"]["wind_cleanup_usd"], 0.0)
        self.assertEqual(by_id["JP-bridge"]["wind_bridge_usd"], 4_146_699.67)
        self.assertEqual(by_id["FOOT"]["wind_asset_usd"], 0.0)
        self.assertEqual(totals["cleanup_ways"], 1)
        self.assertEqual(totals["collapsed_bridges"], 1)
        self.assertEqual(totals["contract"], "WIND_ASSET_IMPACT_CONTRACT.md")

    def test_cli_constants_and_apply(self) -> None:
        constants = json.loads(self._run_cli(["constants"]))
        self.assertEqual(constants["g_inland_3s_10min"], 1.66)
        self.assertEqual(constants["cleanup_usd_per_km_rounded"]["central"], 4979)
        self.assertAlmostEqual(
            constants["bridge_gust_vstar_kmh"]["primary"],
            BRIDGE_GUST_VSTAR_KMH["primary"],
        )
        wind = tree_break_c15_ms() + 1.0
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src = tmp_path / "in.csv"
            out = tmp_path / "out.csv"
            src.write_text(
                "way_id,iso3,highway,accepted,is_bridge,is_tunnel,length_km,"
                "replacement_usd,tree_dens_km2,v_c15_ms\n"
                f"A,USA,residential,1,0,0,1.0,600000,10000,{wind}\n",
                encoding="utf-8",
            )
            main(["apply", str(src), "--output", str(out), "--repo", str(REPO)])
            text = out.read_text(encoding="utf-8")
            self.assertIn("wind_cleanup_usd", text)
            summary = json.loads((tmp_path / "out.csv.summary.json").read_text(encoding="utf-8"))
            self.assertGreater(summary["cleanup_usd"], 0.0)

    def _run_cli(self, args: list[str]) -> str:
        from io import StringIO
        from contextlib import redirect_stdout

        buffer = StringIO()
        with redirect_stdout(buffer):
            code = main(args)
        self.assertEqual(code, 0)
        return buffer.getvalue()

    def test_c15_threshold_rejects_nonpositive_factor(self) -> None:
        with self.assertRaises(ValueError):
            c15_threshold_ms(151.0, 0.0)

    def test_hash_trees_writes_sha256_and_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            blob = tmp_path / "tiny.tif"
            blob.write_bytes(b"II*\x00not-a-real-geotiff-but-hashable")
            manifest = tmp_path / "crowther.manifest.json"
            self.assertEqual(main(["hash-trees", str(blob), "--output", str(manifest)]), 0)
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertIn("sha256", payload)
            self.assertEqual(len(payload["sha256"]), 64)
            self.assertGreater(payload["bytes"], 0)
            self.assertIn("Crowther", payload["source"])

    def _write_tiny_crowther(self, path: Path) -> Path:
        import numpy as np
        from osgeo import gdal, osr

        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(str(path), 4, 2, 1, gdal.GDT_Float64)
        # lon [100,104), lat (8,10] — pixel (100.5, 9.5) is col 0 row 0
        dataset.SetGeoTransform((100.0, 1.0, 0.0, 10.0, 0.0, -1.0))
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        band = dataset.GetRasterBand(1)
        band.SetNoDataValue(-9999.0)
        band.WriteArray(
            np.array(
                [[1000.0, 2000.0, 3000.0, -9999.0], [4000.0, 5000.0, 6000.0, 7000.0]],
                dtype="float64",
            )
        )
        band.FlushCache()
        dataset = None
        return path

    def test_sample_crowther_reads_pixel_not_invented_value(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            geotiff = self._write_tiny_crowther(Path(tmp) / "tiny.tif")
            self.assertEqual(sample_crowther_density(geotiff, 100.5, 9.5), 1000.0)
            self.assertEqual(sample_crowther_density(geotiff, 101.1, 9.1), 2000.0)
            self.assertEqual(sample_crowther_density(geotiff, 100.2, 8.2), 4000.0)
            self.assertTrue(math.isnan(sample_crowther_density(geotiff, 103.2, 9.2)))
            self.assertTrue(math.isnan(sample_crowther_density(geotiff, 0.0, 0.0)))
            attached = attach_crowther_density(
                [
                    {"way_id": "A", "lon": "100.5", "lat": "9.5"},
                    {"way_id": "B", "lon": "", "lat": "9.5"},
                ],
                geotiff,
            )
            self.assertEqual(attached[0]["tree_dens_km2"], 1000.0)
            self.assertTrue(math.isnan(float(attached[1]["tree_dens_km2"])))

    def test_cli_sample_trees_writes_density_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            geotiff = self._write_tiny_crowther(tmp_path / "tiny.tif")
            src = tmp_path / "ways.csv"
            out = tmp_path / "ways.trees.csv"
            src.write_text(
                "way_id,lon,lat\nA,100.5,9.5\nB,103.2,9.2\n",
                encoding="utf-8",
            )
            main(
                [
                    "sample-trees",
                    str(src),
                    "--trees",
                    str(geotiff),
                    "--output",
                    str(out),
                ]
            )
            rows = list(__import__("csv").DictReader(out.open(encoding="utf-8")))
            self.assertEqual(float(rows[0]["tree_dens_km2"]), 1000.0)
            self.assertTrue(math.isnan(float(rows[1]["tree_dens_km2"])))
            summary = json.loads((tmp_path / "ways.trees.csv.summary.json").read_text())
            self.assertEqual(summary["sampled_finite"], 1)
            self.assertEqual(summary["sampled_nodata_or_outside"], 1)

    def test_apply_representative_fixture_zero_when_below_cut_or_no_trees(self) -> None:
        fixture = (
            REPO / "data" / "impact" / "fixtures" / "representative_wind_apply.csv"
        )
        self.assertTrue(fixture.is_file(), "representative apply fixture must exist")
        rows = []
        import csv

        with fixture.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        impacted, totals = apply_rows(rows, REPO)
        by_id = {row["way_id"]: row for row in impacted}
        self.assertIn("cleanup_usd", totals)
        self.assertIn("bridge_usd", totals)
        self.assertEqual(by_id["PH-sparse"]["wind_cleanup_usd"], 0.0)
        self.assertEqual(by_id["CN-below"]["wind_cleanup_usd"], 0.0)
        self.assertEqual(by_id["US-I10"]["wind_bridge_usd"], 0.0)
        self.assertLess(
            by_id["US-I10"]["wind_cleanup_usd"],
            0.01 * float(by_id["US-I10"]["replacement_usd"]),
        )
        self.assertGreater(by_id["US-I10"]["wind_cleanup_usd"], 0.0)
        self.assertGreater(by_id["JP-bridge"]["wind_bridge_usd"], 0.0)


if __name__ == "__main__":
    unittest.main()
