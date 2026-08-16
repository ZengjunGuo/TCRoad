"""Contract tests for 2025 object-level road replacement cost."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


CODE_DIR = Path(__file__).resolve().parents[1]
REPO = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

from apply_road_replacement_value import apply_rows, qa_replacement_totals  # noqa: E402
from assign_road_countries import CountryIndex  # noqa: E402
from build_road_unit_cost_book import assemble_book  # noqa: E402
from extract_osm_motor_roads import haversine_km, way_length_and_mid  # noqa: E402
from road_replacement_value import (  # noqa: E402
    PRICE_BOOK_EUROPE,
    PRICE_BOOK_NATIONAL,
    PRICE_BOOK_ROCKS,
    WORK_TYPE_2L,
    classify_highway,
    parse_lanes,
    replacement_cost,
)


class RoadReplacementValueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.book, cls.meta = assemble_book(REPO)

    def test_planet_year_is_documented_as_20260803(self) -> None:
        contract = (REPO / "methods" / "ROAD_ASSET_VALUATION_CONTRACT.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("planet-260803", contract)
        self.assertIn("2026-08-03", contract)
        self.assertIn("Do not rewind OSM", contract)

    def test_excludes_footways_and_keeps_motor_roads(self) -> None:
        self.assertEqual(classify_highway("motorway")[0], "motorway")
        self.assertTrue(classify_highway("motorway_link")[1])
        self.assertIsNone(classify_highway("footway")[0])
        self.assertIsNone(classify_highway("service")[0])
        self.assertIsNone(classify_highway("track")[0])

    def test_parse_lanes_rejects_junk(self) -> None:
        self.assertEqual(parse_lanes("4"), 4.0)
        self.assertEqual(parse_lanes("2;3"), 2.5)
        self.assertIsNone(parse_lanes("none"))
        self.assertIsNone(parse_lanes("99"))

    def test_philippines_primary_uses_rocks_not_us_price(self) -> None:
        result = replacement_cost(
            self.book,
            length_km=3.2,
            highway="primary",
            iso3="PHL",
            lanes=None,
            surface="asphalt",
        )
        self.assertTrue(result.accepted)
        self.assertEqual(result.price_book, PRICE_BOOK_ROCKS)
        self.assertEqual(result.work_type, WORK_TYPE_2L)
        self.assertEqual(result.lanes_source, "class_default")
        usa = replacement_cost(
            self.book,
            length_km=3.2,
            highway="primary",
            iso3="USA",
            lanes=2,
            surface="asphalt",
        )
        self.assertEqual(usa.price_book, PRICE_BOOK_NATIONAL)
        self.assertGreater(usa.usd_per_km, result.usd_per_km)

    def test_us_motorway_is_national_book(self) -> None:
        result = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="USA",
            lanes=4,
            surface="asphalt",
        )
        self.assertTrue(result.accepted)
        self.assertEqual(result.price_book, PRICE_BOOK_NATIONAL)
        self.assertAlmostEqual(result.usd_per_km, 12_000_000.0, delta=1.0)

    def test_germany_scales_from_european_baseline_not_rocks_eca(self) -> None:
        de = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="DEU",
            lanes=4,
            surface="asphalt",
        )
        bg = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="BGR",
            lanes=4,
            surface="asphalt",
        )
        self.assertEqual(de.price_book, PRICE_BOOK_EUROPE)
        # Bulgaria is EU but upper-middle in some years; if HIC it is Europe,
        # otherwise ROCKS. Either way Germany must be far above ECA ROCKS 4L.
        eca_4l = self.book.rocks["Europe and Central Asia"]["New 4L Expressway"]
        self.assertGreater(de.usd_per_km, 2.0 * eca_4l.central)
        if bg.price_book == PRICE_BOOK_EUROPE:
            self.assertGreater(de.usd_per_km, bg.usd_per_km)

    def test_link_is_not_priced_as_full_four_lane(self) -> None:
        main = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="USA",
            lanes=4,
            surface="asphalt",
        )
        link = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway_link",
            iso3="USA",
            lanes=None,
            surface="asphalt",
        )
        self.assertTrue(link.is_link)
        # One-lane ramp vs 4-lane default: 25%/lane would be 0.25, floored at 0.5.
        self.assertAlmostEqual(link.usd_per_km, 0.5 * main.usd_per_km, delta=1.0)

    def test_bridge_overrides_pavement(self) -> None:
        road = replacement_cost(
            self.book,
            length_km=0.4,
            highway="primary",
            iso3="PHL",
            surface="asphalt",
        )
        bridge = replacement_cost(
            self.book,
            length_km=0.4,
            highway="primary",
            iso3="PHL",
            surface="asphalt",
            bridge="yes",
        )
        self.assertTrue(bridge.is_bridge)
        self.assertGreater(bridge.usd_per_km, road.usd_per_km)

    def test_steep_terrain_nearly_doubles_pavement(self) -> None:
        flat = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="CHN",
            lanes=4,
            surface="asphalt",
            slope_deg=2.0,
        )
        steep = replacement_cost(
            self.book,
            length_km=1.0,
            highway="motorway",
            iso3="CHN",
            lanes=4,
            surface="asphalt",
            slope_deg=30.0,
        )
        self.assertAlmostEqual(steep.usd_per_km / flat.usd_per_km, 1.945, places=2)

    def test_country_lane_median_beats_global_default(self) -> None:
        record = self.book.countries["KEN"]
        record.lane_median_by_class["primary"] = 1.0
        result = replacement_cost(
            self.book,
            length_km=1.0,
            highway="primary",
            iso3="KEN",
            lanes=None,
            surface="asphalt",
        )
        default = replacement_cost(
            self.book,
            length_km=1.0,
            highway="primary",
            iso3="UGA",
            lanes=None,
            surface="asphalt",
        )
        self.assertEqual(result.lanes_source, "country_class_median")
        self.assertLess(result.usd_per_km, default.usd_per_km)

    def test_south_asia_two_lane_does_not_use_42000_outlier(self) -> None:
        band = self.book.rocks["South Asia"]["New 2L Highway"]
        self.assertGreater(band.central, 200_000.0)
        self.assertNotEqual(band.fill_rule, "direct")

    def test_gravel_resurfacing_is_absent_from_the_book(self) -> None:
        for region, types in self.book.rocks.items():
            self.assertNotIn("Gravel Resurfacing", types)
            self.assertIn("New 4L Expressway", types)

    def test_china_uses_national_book(self) -> None:
        result = replacement_cost(
            self.book,
            length_km=10.0,
            highway="motorway",
            iso3="CHN",
            lanes=4,
            surface="asphalt",
        )
        self.assertEqual(result.price_book, PRICE_BOOK_NATIONAL)
        self.assertAlmostEqual(result.replacement_usd, 110_000_000.0, delta=1.0)

    def test_unknown_surface_is_an_expectation_not_a_world_average(self) -> None:
        record = self.book.countries["NGA"]
        record.paved_fraction_by_class["tertiary"] = 0.25
        result = replacement_cost(
            self.book,
            length_km=1.0,
            highway="tertiary",
            iso3="NGA",
            surface=None,
        )
        self.assertTrue(result.surface.startswith("unknown_expected_"))
        self.assertEqual(result.work_type, "New 1L Road")

    def test_way_midpoint_length_is_spherical(self) -> None:
        # One degree of latitude is about 111.2 km.
        length, lon, lat = way_length_and_mid([(120.0, 10.0), (120.0, 11.0)])
        self.assertAlmostEqual(length, haversine_km(120.0, 10.0, 120.0, 11.0), places=6)
        self.assertAlmostEqual(lat, 10.5, places=3)
        self.assertAlmostEqual(lon, 120.0, places=6)

    def test_apply_reports_main_and_no_local_totals(self) -> None:
        rows = [
            {
                "way_id": "1",
                "highway": "motorway",
                "length_km": "2.0",
                "iso3": "USA",
                "lanes": "4",
                "surface": "asphalt",
                "bridge": "",
                "tunnel": "",
            },
            {
                "way_id": "2",
                "highway": "residential",
                "length_km": "10.0",
                "iso3": "USA",
                "lanes": "1",
                "surface": "asphalt",
                "bridge": "",
                "tunnel": "",
            },
            {
                "way_id": "3",
                "highway": "footway",
                "length_km": "1.0",
                "iso3": "USA",
                "lanes": "",
                "surface": "",
                "bridge": "",
                "tunnel": "",
            },
        ]
        valued, totals = apply_rows(rows, REPO)
        self.assertEqual(totals["accepted_ways"], 2)
        self.assertEqual(totals["rejected_ways"], 1)
        self.assertAlmostEqual(totals["replacement_usd"], 2.0 * 12_000_000.0 + 10.0 * 600_000.0)
        self.assertAlmostEqual(totals["replacement_usd_no_local"], 2.0 * 12_000_000.0)
        self.assertEqual(sum(int(row["accepted"]) for row in valued), 2)

    def test_country_index_uses_natural_earth_not_a_stub(self) -> None:
        index = CountryIndex()
        # Inland Kansas, not a 50 m coastline sliver: NE 50m drops some
        # harbour pixels (Manhattan, Miami Beach) into the ocean.
        self.assertEqual(index.iso3(-97.0, 38.0), "USA")
        self.assertEqual(index.iso3(120.9842, 14.5995), "PHL")
        self.assertEqual(index.iso3(116.4074, 39.9042), "CHN")
        self.assertEqual(index.iso3(139.6917, 35.6895), "JPN")
        self.assertEqual(index.iso3(-160.0, 0.0), "")

    def test_unclassified_share_flag_uses_giri_threshold(self) -> None:
        valued = [
            {
                "accepted": 1,
                "iso3": "ZZZ",
                "highway": "unclassified",
                "road_class": "local",
                "length_km": 80.0,
                "replacement_usd": 80.0,
            },
            {
                "accepted": 1,
                "iso3": "ZZZ",
                "highway": "primary",
                "road_class": "primary",
                "length_km": 20.0,
                "replacement_usd": 20.0,
            },
            {
                "accepted": 1,
                "iso3": "USA",
                "highway": "motorway",
                "road_class": "motorway",
                "length_km": 10.0,
                "replacement_usd": 10.0,
            },
        ]
        qa = qa_replacement_totals(valued)
        self.assertIn("ZZZ", qa["countries_flagged_unclassified"])
        self.assertNotIn("USA", qa["countries_flagged_unclassified"])
        self.assertGreater(qa["by_country"]["ZZZ"]["unclassified_share"], 0.60)


if __name__ == "__main__":
    unittest.main()
