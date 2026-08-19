"""Tests for the shared paper-map helpers used by both world-map plates."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest


CODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_DIR))

from paper_map_style import (  # noqa: E402
    FIG_H_MM,
    FIG_W_MM,
    GRIDSPEC,
    INSET_SPECS,
    MM_TO_IN,
    UNIT_COST_BOUNDS,
    UNIT_COST_COLORS,
    UNIT_COST_LABELS,
    align_plate_to_globe,
    country_mean_unit_cost,
    discrete_ramp_origin,
    draw_discrete_ramp,
    inset_caption,
    inset_title_line_y,
    shared_panel_label_x,
    equal_inset_boxes_matching_edges,
    equal_panel_boxes,
    fitted_globe_position,
    inset_data_aspect,
    inset_spec_by_name,
    make_figure,
    make_gridspec,
    point_in_inset_bbox,
    robinson,
    robinson_global_aspect,
    stock_class,
    style_globe,
    unit_cost_class,
    usd_per_km,
)


class PaperMapStyleTests(unittest.TestCase):
    def test_shared_canvas_is_nature_double_column(self):
        self.assertEqual(FIG_W_MM, 183.0)
        self.assertEqual(FIG_H_MM, 150.0)

    def test_inset_boxes_are_the_three_baseline_basins(self):
        names = [title for title, _, _ in INSET_SPECS]
        self.assertEqual(names, ["Gulf Coast", "Pearl River Delta", "Bengal Delta"])
        gulf = inset_spec_by_name("Gulf Coast")[1]
        self.assertEqual(gulf, (-96.0, -89.5, 26.5, 31.9))

    def test_usd_per_km_is_cost_divided_by_length(self):
        self.assertEqual(usd_per_km(12_000_000.0, 8.0), 1_500_000.0)
        self.assertIsNone(usd_per_km(12_000_000.0, 0.0))
        self.assertIsNone(usd_per_km(12_000_000.0, -1.0))
        self.assertIsNone(usd_per_km("x", 2.0))
        self.assertIsNone(usd_per_km(1.0, float("nan")))

    def test_usd_per_km_uses_the_caller_supplied_totals(self):
        # Do not hard-code a country total. Drive the real helper.
        cost = 49230426730519.38
        length = 52249542.670469426
        mean = usd_per_km(cost, length)
        self.assertIsNotNone(mean)
        self.assertAlmostEqual(mean * length, cost, delta=1.0)

    def test_unit_cost_class_follows_the_shared_bounds(self):
        self.assertEqual(unit_cost_class(0.0), 0)
        self.assertEqual(unit_cost_class(UNIT_COST_BOUNDS[1] - 1.0), 0)
        self.assertEqual(unit_cost_class(UNIT_COST_BOUNDS[1]), 1)
        self.assertEqual(unit_cost_class(UNIT_COST_BOUNDS[-2]), 5)
        self.assertIsNone(unit_cost_class(None))
        self.assertIsNone(unit_cost_class(-1.0))

    def test_stock_class_is_in_trillion_usd(self):
        self.assertEqual(stock_class(0.04e12), 0)
        self.assertEqual(stock_class(0.20e12), 2)
        self.assertEqual(stock_class(4.0e12), 5)
        self.assertIsNone(stock_class("bad"))

    def test_point_in_inset_bbox_uses_the_shared_gulf_box(self):
        _, gulf, _ = inset_spec_by_name("Gulf Coast")
        self.assertTrue(point_in_inset_bbox(-92.0, 29.0, gulf))
        self.assertFalse(point_in_inset_bbox(114.0, 22.5, gulf))
        _, prd, _ = inset_spec_by_name("Pearl River Delta")
        self.assertTrue(point_in_inset_bbox(114.0, 22.5, prd))
        self.assertFalse(point_in_inset_bbox(-92.0, 29.0, prd))

    def test_country_mean_unit_cost_uses_usd_per_km(self):
        table = {
            "AAA": {"replacement_usd": 250.0, "length_km": 5.0},
            "BBB": {"replacement_usd": 10.0, "length_km": 0.0},
        }
        means = country_mean_unit_cost(table)
        self.assertEqual(means["AAA"], usd_per_km(250.0, 5.0))
        self.assertNotIn("BBB", means)

    def test_point_in_inset_bbox_rejects_reversed_edges(self):
        with self.assertRaises(ValueError):
            point_in_inset_bbox(0.0, 0.0, (10.0, 0.0, 0.0, 1.0))

    def test_discrete_ramp_is_centered_on_the_shared_axis(self):
        n = len(UNIT_COST_COLORS)
        chip_w = 0.048
        self.assertAlmostEqual(discrete_ramp_origin(n, chip_w, 0.5), 0.5 - n * chip_w / 2.0)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        x0, x1, y0, y1 = draw_discrete_ramp(
            ax, UNIT_COST_COLORS, UNIT_COST_LABELS, cx=0.5, y=0.30, chip_w=chip_w, chip_h=0.20,
        )
        self.assertAlmostEqual(x0, 0.5 - n * chip_w / 2.0)
        self.assertAlmostEqual(x1, 0.5 + n * chip_w / 2.0)
        self.assertEqual(len(ax.patches), n)
        self.assertAlmostEqual((x0 + x1) / 2.0, 0.5)
        plt.close(fig)

    def test_fitted_globe_position_matches_data_aspect_and_stays_centered(self):
        fig_size = (FIG_W_MM * MM_TO_IN, FIG_H_MM * MM_TO_IN)
        box = (0.05, 0.50, 0.93, 0.40)
        aspect = 2.0
        x0, y0, width, height = fitted_globe_position(box, fig_size, aspect)
        display_aspect = (width * fig_size[0]) / (height * fig_size[1])
        self.assertAlmostEqual(display_aspect, aspect, places=9)
        self.assertAlmostEqual(x0 + width / 2.0, box[0] + box[2] / 2.0)
        self.assertLess(width, box[2])
        self.assertAlmostEqual(height, box[3])
        self.assertAlmostEqual(y0, box[1])

    def test_equal_panel_boxes_share_the_caller_supplied_edges(self):
        left, right = 0.17, 0.83
        boxes = equal_panel_boxes(left, right, 0.04, 0.28, n=3, wspace=0.055)
        self.assertEqual(len(boxes), 3)
        self.assertAlmostEqual(boxes[0][0], left)
        last_x, _, last_w, _ = boxes[-1]
        self.assertAlmostEqual(last_x + last_w, right)
        widths = [box[2] for box in boxes]
        self.assertTrue(all(abs(w - widths[0]) < 1e-12 for w in widths))

    def test_equal_inset_boxes_keep_oval_edges_for_shared_inset_aspects(self):
        left, right = 0.16, 0.84
        aspects = [inset_data_aspect(spec[1]) for spec in INSET_SPECS]
        fig_size = (FIG_W_MM * MM_TO_IN, FIG_H_MM * MM_TO_IN)
        boxes = equal_inset_boxes_matching_edges(
            left, right, 0.046, 0.30, aspects, fig_size, GRIDSPEC["wspace"],
        )
        self.assertAlmostEqual(boxes[0][0], left)
        last_x, _, last_w, _ = boxes[-1]
        self.assertAlmostEqual(last_x + last_w, right)
        for box, aspect in zip(boxes, aspects):
            display_aspect = (box[2] * fig_size[0]) / (box[3] * fig_size[1])
            self.assertAlmostEqual(display_aspect, aspect, places=9)

    def test_inset_caption_is_letter_dot_space_title(self):
        self.assertEqual(inset_caption("b", "Gulf Coast"), "b. Gulf Coast")
        self.assertEqual(inset_caption("d", "Bengal Delta"), "d. Bengal Delta")

    def test_shared_panel_label_x_is_the_first_inset_left_edge(self):
        boxes = [(0.162, 0.05, 0.24, 0.22), (0.42, 0.05, 0.24, 0.21), (0.68, 0.05, 0.24, 0.20)]
        self.assertEqual(shared_panel_label_x(boxes), boxes[0][0])

    def test_inset_title_line_y_uses_the_tallest_frame(self):
        boxes = [(0.16, 0.05, 0.24, 0.20), (0.42, 0.048, 0.24, 0.22), (0.68, 0.06, 0.24, 0.19)]
        pad = 0.012
        self.assertAlmostEqual(inset_title_line_y(boxes, pad=pad), 0.048 + 0.22 + pad)

    def test_align_plate_to_globe_pins_insets_to_the_oval_edges(self):
        import matplotlib
        matplotlib.use("Agg")
        try:
            import cartopy.crs as ccrs
        except ImportError:
            self.skipTest("cartopy is required for the shipped oval-fit helper")
        fig = make_figure()
        gs = make_gridspec(fig)
        globe = fig.add_subplot(gs[0, :], projection=robinson())
        style_globe(globe)
        legend = fig.add_subplot(gs[1, :])
        insets = [
            fig.add_subplot(gs[2, 2 * j:2 * j + 2], projection=ccrs.PlateCarree())
            for j in range(3)
        ]
        left, right = align_plate_to_globe(
            globe, legend, insets, [spec[1] for spec in INSET_SPECS],
        )
        gpos = globe.get_position()
        self.assertAlmostEqual(gpos.x0, left)
        self.assertAlmostEqual(gpos.x1, right)
        self.assertAlmostEqual(insets[0].get_position().x0, left)
        self.assertAlmostEqual(insets[-1].get_position().x1, right)
        self.assertLess(legend.get_position().x0, left)
        self.assertGreater(legend.get_position().x1, right)
        fig_w, fig_h = fig.get_size_inches()
        display_aspect = (gpos.width * fig_w) / (gpos.height * fig_h)
        self.assertAlmostEqual(display_aspect, robinson_global_aspect(), places=6)
        import matplotlib.pyplot as plt
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
