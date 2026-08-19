"""Contract tests for compact C15 → valued-object wind join."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


CODE_DIR = Path(__file__).resolve().parents[1]
REPO = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from road_wind_asset_impact import (  # noqa: E402
    apply_rows,
    escobedo_cleanup_usd_per_km,
    main as kernel_main,
    tree_break_c15_ms,
)
from road_wind_object_join import (  # noqa: E402
    COMPACT_WIND_VARIABLE,
    METHOD_DOMAIN_PENDING_EVENT_POSITIONS,
    _tree_fail_prob_array,
    attach_extract_coordinates,
    canonical_longitude,
    join_event_wind,
    main as join_main,
    sample_compact_index,
    sample_crowther_many,
    score_historical,
    write_compact_fixture,
)


class WindObjectJoinTests(unittest.TestCase):
    def test_pending_eight_are_the_frozen_ids(self) -> None:
        self.assertEqual(
            METHOD_DOMAIN_PENDING_EVENT_POSITIONS,
            frozenset({11902, 11944, 12357, 50194, 62311, 68925, 72126, 86977}),
        )

    def test_periodic_longitude_wraps_dateline(self) -> None:
        wrapped = canonical_longitude(np.asarray([185.45, -174.55, 0.0]))
        self.assertAlmostEqual(wrapped[0], -174.55, places=5)
        self.assertAlmostEqual(wrapped[1], -174.55, places=5)
        self.assertAlmostEqual(wrapped[2], 0.0, places=12)

    def _tiny_crowther(self, path: Path) -> Path:
        from osgeo import gdal, osr

        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(str(path), 4, 2, 1, gdal.GDT_Float64)
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

    def test_join_reads_compact_event_max_not_a_mock_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            compact = write_compact_fixture(
                tmp_path / "00000.nc",
                lat=[20.00, 20.00, 20.05],
                lon=[0.00, 0.05, 0.00],
                wind=[31.25, 18.0, 22.0],
                event_position=0,
                event_id="stream0000-year1995-track000000",
            )
            rows = [
                {"way_id": "ON", "lon": "0.01", "lat": "20.01", "highway": "motorway"},
                {"way_id": "OFF", "lon": "10.0", "lat": "0.0", "highway": "motorway"},
                {"way_id": "DATELINE", "lon": "-174.55", "lat": "20.0", "highway": "primary"},
            ]
            # dateline way is not in this compact; add a Pacific centroid
            compact_pacific = write_compact_fixture(
                tmp_path / "00001.nc",
                lat=[20.00],
                lon=[185.45],
                wind=[27.5],
                event_position=1,
                event_id="pacific",
            )
            joined, meta = join_event_wind(rows[:2], compact)
            by_id = {row["way_id"]: row for row in joined}
            self.assertEqual(meta["wind_variable"], COMPACT_WIND_VARIABLE)
            self.assertAlmostEqual(float(by_id["ON"]["v_c15_ms"]), 31.25, places=5)
            self.assertTrue(math.isnan(float(by_id["OFF"]["v_c15_ms"])))
            self.assertEqual(meta["ways_with_footprint_wind"], 1)
            self.assertEqual(meta["ways_outside_footprint"], 1)

            pacific, _ = join_event_wind([rows[2]], compact_pacific)
            self.assertAlmostEqual(float(pacific[0]["v_c15_ms"]), 27.5, places=5)

    def test_pending_compact_is_rejected_not_scored(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            compact = write_compact_fixture(
                Path(tmp) / "11902.nc",
                lat=[20.0],
                lon=[0.0],
                wind=[40.0],
                event_position=11902,
            )
            with self.assertRaises(ValueError) as error:
                join_event_wind(
                    [{"way_id": "A", "lon": "0.0", "lat": "20.0"}], compact
                )
            self.assertIn("METHOD_DOMAIN_PENDING", str(error.exception))

    def test_attach_extract_coordinates_uses_same_way_id(self) -> None:
        valued = [{"way_id": "37", "highway": "residential", "accepted": "1"}]
        extract = [{"way_id": "37", "lon": "-1.821422", "lat": "52.554016"}]
        attached = attach_extract_coordinates(valued, extract)
        self.assertEqual(attached[0]["lon"], "-1.821422")
        self.assertEqual(attached[0]["lat"], "52.554016")

    def test_crowther_batch_is_nearest_pixel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            geotiff = self._tiny_crowther(Path(tmp) / "tiny.tif")
            values = sample_crowther_many(
                geotiff, np.asarray([100.5, 103.2]), np.asarray([9.5, 9.2])
            )
            self.assertEqual(values[0], 1000.0)
            self.assertTrue(math.isnan(values[1]))

    def test_score_event_entry_writes_cleanup_and_bridge_usd(self) -> None:
        wind = tree_break_c15_ms() + 2.0
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            compact = write_compact_fixture(
                tmp_path / "00000.nc",
                lat=[20.0, -20.0],
                lon=[0.0, 140.0],
                wind=[wind, 80.0],
                event_position=7,
            )
            src = tmp_path / "ways.csv"
            out = tmp_path / "impact.csv"
            src.write_text(
                "way_id,iso3,highway,accepted,is_bridge,is_tunnel,length_km,"
                "replacement_usd,tree_dens_km2,lon,lat\n"
                "US-I10,USA,motorway,1,0,0,2.0,24000000,10000,0.0,20.0\n"
                "JP-bridge,JPN,primary,1,1,0,0.4,4146699.67,8000,140.0,-20.0\n"
                "MISS,USA,residential,1,0,0,1.0,600000,10000,50.0,50.0\n",
                encoding="utf-8",
            )
            self.assertEqual(
                join_main(
                    [
                        "score-event",
                        str(src),
                        "--compact",
                        str(compact),
                        "--output",
                        str(out),
                        "--repo",
                        str(REPO),
                    ]
                ),
                0,
            )
            rows = list(csv.DictReader(out.open(encoding="utf-8")))
            by_id = {row["way_id"]: row for row in rows}
            self.assertIn("wind_cleanup_usd", by_id["US-I10"])
            self.assertIn("wind_bridge_usd", by_id["JP-bridge"])
            self.assertAlmostEqual(
                float(by_id["US-I10"]["wind_cleanup_usd"]),
                2.0 * escobedo_cleanup_usd_per_km("central"),
            )
            self.assertEqual(float(by_id["US-I10"]["wind_bridge_usd"]), 0.0)
            self.assertEqual(float(by_id["JP-bridge"]["wind_cleanup_usd"]), 0.0)
            self.assertGreater(float(by_id["JP-bridge"]["wind_bridge_usd"]), 0.0)
            self.assertEqual(float(by_id["MISS"]["wind_cleanup_usd"]), 0.0)
            self.assertTrue(math.isnan(float(by_id["MISS"]["v_c15_ms"])))
            summary = json.loads((tmp_path / "impact.csv.summary.json").read_text())
            self.assertIn("cleanup_usd", summary)
            self.assertIn("bridge_usd", summary)
            # Applying the kernel twice on the joined table is consistent.
            joined, _ = join_event_wind(list(csv.DictReader(src.open())), compact)
            first, totals_a = apply_rows(joined, REPO)
            second, totals_b = apply_rows(joined, REPO)
            self.assertEqual(totals_a["cleanup_usd"], totals_b["cleanup_usd"])
            self.assertEqual(totals_a["bridge_usd"], totals_b["bridge_usd"])
            self.assertEqual(
                first[0]["wind_cleanup_usd"], second[0]["wind_cleanup_usd"]
            )

    def test_kernel_apply_still_the_dollar_entry_after_join(self) -> None:
        wind = tree_break_c15_ms() + 1.0
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            compact = write_compact_fixture(
                tmp_path / "e.nc",
                lat=[20.0],
                lon=[0.0],
                wind=[wind],
            )
            src = tmp_path / "ways.csv"
            joined = tmp_path / "joined.csv"
            out = tmp_path / "out.csv"
            src.write_text(
                "way_id,iso3,highway,accepted,is_bridge,is_tunnel,length_km,"
                "replacement_usd,tree_dens_km2,lon,lat\n"
                f"A,USA,residential,1,0,0,1.0,600000,10000,0.0,20.0\n",
                encoding="utf-8",
            )
            join_main(["join-event", str(src), "--compact", str(compact), "--output", str(joined)])
            kernel_main(["apply", str(joined), "--output", str(out), "--repo", str(REPO)])
            summary = json.loads((tmp_path / "out.csv.summary.json").read_text())
            self.assertGreater(summary["cleanup_usd"], 0.0)
            text = out.read_text(encoding="utf-8")
            self.assertIn("wind_cleanup_usd", text)
            self.assertIn("wind_bridge_usd", text)

    def test_sample_compact_index_misses_empty_index(self) -> None:
        wind = sample_compact_index(np.asarray([0.0]), np.asarray([0.0]), {})
        self.assertTrue(math.isnan(float(wind[0])))

    def test_vector_tree_prob_matches_kernel(self) -> None:
        from road_wind_asset_impact import tree_fail_prob

        values = np.asarray([0.0, -1.0, float("nan"), 10.0, 10000.0, 50000.0])
        array = _tree_fail_prob_array(values)
        for i, value in enumerate(values):
            self.assertEqual(array[i], tree_fail_prob(float(value)))

    def test_score_historical_skips_pending_and_uses_compact_wind(self) -> None:
        wind = tree_break_c15_ms() + 3.0
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            valued_dir = root / "valued"
            extract_dir = root / "extract"
            compact_dir = root / "compact"
            out_dir = root / "out"
            valued_dir.mkdir()
            extract_dir.mkdir()
            compact_dir.mkdir()
            (valued_dir / "ways-0000.valued.csv").write_text(
                "way_id,iso3,highway,accepted,reason,road_class,is_link,is_bridge,"
                "is_tunnel,surface,lanes_used,lanes_source,terrain_class,work_type,"
                "price_book,length_km,usd_per_km,replacement_usd,"
                "replacement_usd_low,replacement_usd_high\n"
                "1,USA,motorway,1,ok,motorway,0,0,0,paved,4,osm_tag,unknown,"
                "New 4L Expressway,national,2.0,12000000,24000000,1,1\n",
                encoding="utf-8",
            )
            (extract_dir / "ways-0000.csv").write_text(
                "way_id,highway,lanes,surface,bridge,tunnel,lit,n_nodes,length_km,lon,lat,iso3\n"
                "1,motorway,4,paved,0,0,,2,2.0,0.0,20.0,USA\n",
                encoding="utf-8",
            )
            write_compact_fixture(
                compact_dir / "00000.nc",
                lat=[20.0],
                lon=[0.0],
                wind=[wind],
                event_position=0,
                event_id="keep",
                weight=0.002,
            )
            write_compact_fixture(
                compact_dir / "11902.nc",
                lat=[20.0],
                lon=[0.0],
                wind=[80.0],
                event_position=11902,
                event_id="pending",
            )
            geotiff = self._tiny_crowther(root / "trees.tif")
            # Crowther fixture does not cover lon=0; tree density is NaN → cleanup 0
            # unless we put the way on the tiny raster. Use lon=100.5 lat=9.5.
            (extract_dir / "ways-0000.csv").write_text(
                "way_id,highway,lanes,surface,bridge,tunnel,lit,n_nodes,length_km,lon,lat,iso3\n"
                "1,motorway,4,paved,0,0,,2,2.0,100.5,9.5,USA\n",
                encoding="utf-8",
            )
            write_compact_fixture(
                compact_dir / "00000.nc",
                lat=[9.50],
                lon=[100.50],
                wind=[wind],
                event_position=0,
                event_id="keep",
                weight=0.002,
            )
            write_compact_fixture(
                compact_dir / "11902.nc",
                lat=[9.50],
                lon=[100.50],
                wind=[80.0],
                event_position=11902,
                event_id="pending",
            )
            summary = score_historical(
                valued_dir=valued_dir,
                extract_dir=extract_dir,
                compact_dir=compact_dir,
                trees=geotiff,
                output_dir=out_dir,
                repo=REPO,
            )
            self.assertEqual(summary["compact_events_scored"], 1)
            self.assertIn(11902, summary["method_domain_pending_event_positions"])
            self.assertFalse(summary["method_domain_pending_contribute_dollars"])
            self.assertGreater(summary["cleanup_usd_sum"], 0.0)
            self.assertEqual(summary["bridge_usd_sum"], 0.0)
            ways = list(csv.DictReader((out_dir / "way_wind_asset.csv").open(encoding="utf-8")))
            self.assertAlmostEqual(float(ways[0]["max_v_c15_ms"]), wind, places=5)
            self.assertAlmostEqual(
                float(ways[0]["wind_cleanup_usd_sum"]),
                2.0 * escobedo_cleanup_usd_per_km("central") * 0.1,
                places=5,
            )
            events = list(csv.DictReader((out_dir / "event_wind_asset.csv").open()))
            self.assertEqual(len(events), 1)
            self.assertEqual(int(events[0]["event_position"]), 0)


if __name__ == "__main__":
    unittest.main()
