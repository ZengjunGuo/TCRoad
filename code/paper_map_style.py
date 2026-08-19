"""Shared Nature-family world-map language for both paper plates.

Both ``plot_global_tc_roads_figure`` and ``plot_road_replacement_cost_figure``
import canvas, projection, inset boxes, locators, and scale bars from here
so the two plates cannot drift.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, fontManager, findfont
from matplotlib.patches import Rectangle

try:
    import cartopy.crs as ccrs
except ImportError:  # pragma: no cover - tests do not need Cartopy
    ccrs = None


# Academic Figure Skill Typography Baseline — COPY VERBATIM
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    "hatch.linewidth": 0.28,
    "hatch.color": "#222222",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
BLACK = "#222222"
GREY = "#999999"
LOCATOR_BLUE = "#004C73"
OCEAN = "#F4F8FB"
LAND = "#F8F7F2"
MAP_FACE = "#F4F4F2"
COAST = "#666666"

MM_TO_IN = 1.0 / 25.4
FIG_W_MM = 183.0
FIG_H_MM = 150.0
PANEL_FONT = 8.0
TITLE_FONT = 7.5
LEGEND_TITLE_FONT = 7.0
LEGEND_ITEM_FONT = 6.5
INSET_TITLE_FONT = 7.2
SCALE_FONT = 6.2

# Same three basins on both plates.
INSET_SPECS = [
    ("Gulf Coast", (-96.0, -89.5, 26.5, 31.9), 200),
    ("Pearl River Delta", (112.6, 115.0, 21.8, 23.8), 50),
    ("Bengal Delta", (88.5, 92.8, 20.8, 24.4), 100),
]

# Shared gridspec: hero, external legend strip, three equal insets.
# left/right are identical for the globe row and the inset row.
# Legend row is tall enough for a title-above-swatches pair plus a note row.
GRIDSPEC = dict(
    nrows=3,
    ncols=6,
    height_ratios=[1.78, 0.70, 1.00],
    width_ratios=[1, 1, 1, 1, 1, 1],
    left=0.050,
    right=0.980,
    bottom=0.046,
    top=0.938,
    hspace=0.13,
    wspace=0.055,
)

ROAD_NAMES = ["Highways", "Primary roads", "Secondary roads", "Tertiary roads", "Local roads"]
ROAD_COLORS = ["#A80000", "#FF0000", "#FFAA00", "#A8A800", "#FFFF00"]
ROAD_WIDTHS = [0.88, 0.68, 0.50, 0.34, 0.20]
ROAD_ALPHAS = [0.96, 0.90, 0.80, 0.68, 0.54]
# Object-level insets have 3–5×10^5 ways; lower local alpha keeps overplotting from
# collapsing the sequential unit-cost ramp into a single dark wash.
COST_ROAD_ALPHAS = [0.96, 0.88, 0.70, 0.40, 0.16]
NO_DATA_COLOR = "#E4E4E4"
INSET_LAND = "#EDE8DF"
INSET_WATER = "#E4EEF4"

# Country-mean / segment unit cost, million 2025 USD per km.
UNIT_COST_BOUNDS = np.array([0.0, 0.30e6, 0.60e6, 1.00e6, 1.50e6, 2.50e6, np.inf])
UNIT_COST_COLORS = ["#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5", "#08306B"]
UNIT_COST_LABELS = ["<0.30", "0.30–0.60", "0.60–1.00", "1.00–1.50", "1.50–2.50", "≥2.50"]

STOCK_BOUNDS_T = np.array([0.0, 0.05, 0.15, 0.40, 1.00, 3.00, np.inf])
STOCK_COLORS = ["#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5", "#08306B"]
STOCK_LABELS = ["<0.05", "0.05–0.15", "0.15–0.40", "0.40–1.0", "1.0–3.0", "≥3.0"]


def configure_runtime_assets(font_dir, cartopy_data=None):
    """Register the frozen Liberation Sans bundle (and optional Cartopy data)."""
    font_dir = Path(font_dir).resolve()
    font_files = sorted(font_dir.rglob("LiberationSans*.ttf"))
    if not font_files:
        raise FileNotFoundError(f"No Liberation Sans TTF files found under {font_dir}")
    for font_path in font_files:
        fontManager.addfont(str(font_path))
    resolved = Path(findfont(
        FontProperties(family="Liberation Sans"),
        fallback_to_default=False,
    )).resolve()
    if font_dir not in resolved.parents:
        raise RuntimeError(f"Liberation Sans resolved outside frozen font bundle: {resolved}")
    mpl.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "Helvetica"]
    if cartopy_data is not None:
        import cartopy
        cartopy_data = Path(cartopy_data).resolve()
        if not cartopy_data.is_dir():
            raise FileNotFoundError(f"Cartopy data directory not found: {cartopy_data}")
        cartopy.config["data_dir"] = str(cartopy_data)
    return resolved


def robinson():
    if ccrs is None:
        raise ImportError("cartopy is required to construct the shared Robinson(0) globe")
    return ccrs.Robinson(central_longitude=0)


def make_figure():
    """Exact Nature double-column canvas shared by both plates."""
    return plt.figure(
        figsize=(FIG_W_MM * MM_TO_IN, FIG_H_MM * MM_TO_IN),
        facecolor="white",
    )


def make_gridspec(fig):
    return fig.add_gridspec(**GRIDSPEC)


def robinson_global_aspect():
    """Width ÷ height of a set_global Robinson(0) globe in projection metres."""
    proj = robinson()
    x0, x1 = proj.x_limits
    y0, y1 = proj.y_limits
    width = float(x1) - float(x0)
    height = float(y1) - float(y0)
    if height <= 0.0:
        raise ValueError("Robinson y-limits do not span a positive height")
    return width / height


def inset_data_aspect(bbox):
    """PlateCarree data aspect (Δlon / Δlat) of an INSET_SPECS box."""
    if len(bbox) != 4:
        raise ValueError("bbox must be (xmin, xmax, ymin, ymax)")
    xmin, xmax, ymin, ymax = (float(v) for v in bbox)
    if xmax <= xmin or ymax <= ymin:
        raise ValueError("bbox edges are reversed or empty")
    return (xmax - xmin) / (ymax - ymin)


def fitted_globe_position(axes_bbox, fig_size_in, data_aspect):
    """Tight-fit a rectangle of ``data_aspect`` inside ``axes_bbox``.

    ``axes_bbox`` is ``(x0, y0, width, height)`` in figure fraction.
    ``fig_size_in`` is ``(fig_w, fig_h)`` in inches. The returned box has
    the same display aspect as ``data_aspect`` and stays centered.
    """
    if data_aspect <= 0.0:
        raise ValueError("data_aspect must be positive")
    x0, y0, width, height = (float(v) for v in axes_bbox)
    fig_w, fig_h = (float(v) for v in fig_size_in)
    if width <= 0.0 or height <= 0.0 or fig_w <= 0.0 or fig_h <= 0.0:
        raise ValueError("axes bbox and figure size must be positive")
    box_aspect = (width * fig_w) / (height * fig_h)
    if box_aspect > data_aspect:
        new_width = height * fig_h * data_aspect / fig_w
        new_x0 = x0 + (width - new_width) / 2.0
        return new_x0, y0, new_width, height
    new_height = width * fig_w / data_aspect / fig_h
    new_y0 = y0 + (height - new_height) / 2.0
    return x0, new_y0, width, new_height


def equal_panel_boxes(left, right, bottom, top, n=3, wspace=None):
    """n equal boxes whose outer left/right are exactly ``left`` / ``right``.

    ``wspace`` uses the matplotlib convention: gap ÷ average box width.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if right <= left or top <= bottom:
        raise ValueError("right/top must exceed left/bottom")
    if wspace is None:
        wspace = float(GRIDSPEC["wspace"])
    wspace = float(wspace)
    if wspace < 0.0:
        raise ValueError("wspace must be non-negative")
    span = right - left
    width = span / (n + (n - 1) * wspace)
    gap = wspace * width
    height = top - bottom
    boxes = []
    x = left
    for _ in range(n):
        boxes.append((x, bottom, width, height))
        x += width + gap
    return boxes


def equal_inset_boxes_matching_edges(left, right, row_bottom, row_top, aspects, fig_size_in, wspace=None):
    """Three (or n) map boxes that share ``left``/``right`` and keep each aspect.

    Widths stay equal so the row still meets the globe oval. Heights are
    derived from each PlateCarree aspect and top-aligned in the row.
    """
    aspects = [float(a) for a in aspects]
    if any(a <= 0.0 for a in aspects):
        raise ValueError("aspects must be positive")
    n = len(aspects)
    slots = equal_panel_boxes(left, right, row_bottom, row_top, n=n, wspace=wspace)
    fig_w, fig_h = (float(v) for v in fig_size_in)
    max_h = row_top - row_bottom
    boxes = []
    for (x, _y, width, _h), aspect in zip(slots, aspects):
        height = (width * fig_w) / (aspect * fig_h)
        if height > max_h + 1e-9:
            raise ValueError(
                f"inset row height {max_h:.4f} is shorter than {height:.4f} "
                "needed to keep the globe's left/right edges"
            )
        boxes.append((x, row_top - height, width, height))
    return boxes


def fit_axes_to_robinson_oval(ax):
    """Resize a set_global Robinson axes so the rectangle is the oval bbox."""
    ax.set_global()
    pos = ax.get_position()
    x0, y0, width, height = fitted_globe_position(
        (pos.x0, pos.y0, pos.width, pos.height),
        ax.figure.get_size_inches(),
        robinson_global_aspect(),
    )
    ax.set_position([x0, y0, width, height])
    ax.set_aspect("equal", adjustable="box")
    return ax.get_position()


def align_strip_to_edges(ax, left, right):
    """Keep an axes on its current y-band but pin it to ``[left, right]``."""
    if right <= left:
        raise ValueError("right must exceed left")
    pos = ax.get_position()
    ax.set_position([left, pos.y0, right - left, pos.height])
    return ax.get_position()


def align_plate_to_globe(globe_ax, legend_ax, inset_axes, inset_bboxes=None):
    """Fit the globe to its oval, then pin legend + inset row to those edges.

    Returns ``(left, right)`` in figure fraction — the shared visible x-edges
    for the oval, the inset frames, and panel letters a/b.
    """
    fit_axes_to_robinson_oval(globe_ax)
    pos = globe_ax.get_position()
    left, right = float(pos.x0), float(pos.x1)
    # Legend stays on the full gridspec band so ramp labels do not collide.
    # The strip is still centered under the oval because the oval is centered.
    inset_axes = list(inset_axes)
    if inset_axes:
        y0 = min(ax.get_position().y0 for ax in inset_axes)
        y1 = max(ax.get_position().y1 for ax in inset_axes)
        if inset_bboxes is None:
            inset_bboxes = [spec[1] for spec in INSET_SPECS[:len(inset_axes)]]
        if len(inset_bboxes) != len(inset_axes):
            raise ValueError("inset_bboxes must match inset_axes")
        aspects = [inset_data_aspect(bbox) for bbox in inset_bboxes]
        boxes = equal_inset_boxes_matching_edges(
            left, right, y0, y1, aspects, globe_ax.figure.get_size_inches(),
            GRIDSPEC["wspace"],
        )
        for ax, box, bbox in zip(inset_axes, boxes, inset_bboxes):
            ax.set_extent(bbox, crs=ccrs.PlateCarree())
            ax.set_position(box)
            # Keep the box we just computed; the extent already matches its aspect.
            ax.set_aspect("equal", adjustable="datalim")
    return left, right


def usd_per_km(replacement_usd, length_km):
    """Object- or country-level mean: replacement cost ÷ motor length."""
    try:
        cost = float(replacement_usd)
        length = float(length_km)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(cost) or not math.isfinite(length) or length <= 0.0:
        return None
    return cost / length


def _class_index(value, bounds):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or number < 0.0:
        return None
    for i in range(len(bounds) - 1):
        if bounds[i] <= number < bounds[i + 1]:
            return i
    return len(bounds) - 2


def unit_cost_class(value_usd_per_km):
    """Discrete class for the shared unit-cost colorbar (0..5, or None)."""
    return _class_index(value_usd_per_km, UNIT_COST_BOUNDS)


def stock_class(replacement_usd):
    """Discrete class for total replacement cost in trillion 2025 USD."""
    try:
        trillion = float(replacement_usd) / 1.0e12
    except (TypeError, ValueError):
        return None
    return _class_index(trillion, STOCK_BOUNDS_T)


def point_in_inset_bbox(lon, lat, bbox):
    """True if (lon, lat) lies in an ``INSET_SPECS`` box (xmin, xmax, ymin, ymax)."""
    if len(bbox) != 4:
        raise ValueError("bbox must be (xmin, xmax, ymin, ymax)")
    xmin, xmax, ymin, ymax = (float(v) for v in bbox)
    if xmin > xmax or ymin > ymax:
        raise ValueError("bbox edges are reversed")
    try:
        x = float(lon)
        y = float(lat)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(x) or not math.isfinite(y):
        return False
    return xmin <= x <= xmax and ymin <= y <= ymax


def country_mean_unit_cost(by_country):
    """Map ISO3 → replacement_usd / length_km for every country with length."""
    out = {}
    for iso, row in by_country.items():
        value = usd_per_km(row.get("replacement_usd"), row.get("length_km"))
        if value is not None:
            out[iso] = value
    return out


def inset_spec_by_name(name):
    for title, bbox, scale_km in INSET_SPECS:
        if title == name:
            return title, bbox, scale_km
    raise KeyError(name)


def style_globe(ax):
    """Coast, land, and oval spine used by both world maps."""
    ax.set_global()
    ax.set_facecolor(MAP_FACE)
    ax.spines["geo"].set_edgecolor("#444444")
    ax.spines["geo"].set_linewidth(0.6)
    ax.spines["geo"].set_visible(True)
    return ax


def style_inset(ax, land=None, water=None):
    """10 m land/coast/borders — the scale Nature insets actually print at."""
    if ccrs is None:
        raise ImportError("cartopy is required to style an inset")
    import cartopy.feature as cfeature

    land = INSET_LAND if land is None else land
    water = INSET_WATER if water is None else water
    ax.set_facecolor(water)
    land_artist = ax.add_feature(
        cfeature.LAND.with_scale("10m"), facecolor=land, edgecolor="none", zorder=1,
    )
    if hasattr(land_artist, "set_rasterized"):
        land_artist.set_rasterized(True)
    ax.coastlines(resolution="10m", color="#6F6F6F", linewidth=0.35, zorder=3)
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), edgecolor="#AAAAAA", linewidth=0.25, zorder=3)
    ax.spines["geo"].set_edgecolor(BLACK)
    ax.spines["geo"].set_linewidth(0.75)
    return ax


def draw_scale_bar(ax, bbox, length_km, x_fraction=0.08):
    """Inset scale bar — same geometry on both plates."""
    xmin, xmax, ymin, ymax = bbox
    lat0 = ymin + 0.12 * (ymax - ymin)
    lon0 = xmin + x_fraction * (xmax - xmin)
    deg = length_km / (111.32 * max(math.cos(math.radians(lat0)), 0.15))
    nseg = 2
    for j in range(nseg):
        x0 = lon0 + deg * j / nseg
        x1 = lon0 + deg * (j + 1) / nseg
        ax.plot(
            [x0, x1], [lat0, lat0], transform=ccrs.PlateCarree(),
            color=BLACK if j % 2 == 0 else "white", lw=2.2,
            solid_capstyle="butt", zorder=20,
        )
        ax.plot(
            [x0, x1], [lat0, lat0], transform=ccrs.PlateCarree(),
            color=BLACK, lw=0.35, solid_capstyle="butt", zorder=21,
        )
    tick_h = 0.025 * (ymax - ymin)
    for x in (lon0, lon0 + deg / 2, lon0 + deg):
        ax.plot(
            [x, x], [lat0 - tick_h, lat0 + tick_h],
            transform=ccrs.PlateCarree(), color=BLACK, lw=0.45, zorder=22,
        )
    halo = [path_effects.withStroke(linewidth=1.5, foreground="white")]
    ax.text(
        lon0, lat0 - 0.055 * (ymax - ymin), "0",
        transform=ccrs.PlateCarree(), fontsize=SCALE_FONT,
        ha="center", va="top", zorder=22, path_effects=halo,
    )
    ax.text(
        lon0 + deg, lat0 - 0.055 * (ymax - ymin), f"{length_km} km",
        transform=ccrs.PlateCarree(), fontsize=SCALE_FONT,
        ha="center", va="top", zorder=22, path_effects=halo,
    )


def add_numbered_locator(global_ax, bbox, locator_number):
    """Numbered locator box on the Robinson globe."""
    xmin, xmax, ymin, ymax = bbox
    halo = Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin, transform=ccrs.PlateCarree(),
        fill=False, edgecolor="white", linewidth=1.55, zorder=14,
    )
    global_ax.add_patch(halo)
    rect = Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin, transform=ccrs.PlateCarree(),
        fill=False, edgecolor=LOCATOR_BLUE, linewidth=0.72, zorder=15,
    )
    global_ax.add_patch(rect)
    global_ax.text(
        xmax, ymax, str(locator_number), transform=ccrs.PlateCarree(),
        ha="left", va="bottom", fontsize=6.0, fontweight="bold", color=LOCATOR_BLUE,
        bbox={
            "boxstyle": "circle,pad=0.10",
            "facecolor": "white",
            "edgecolor": LOCATOR_BLUE,
            "linewidth": 0.55,
        },
        zorder=16, clip_on=False,
    )


def label_hero(fig, letter, title, x=None):
    """Panel-a letter sits on the shared visible left edge of the oval."""
    x = GRIDSPEC["left"] if x is None else float(x)
    fig.text(x, 0.972, letter, fontsize=PANEL_FONT, fontweight="bold", ha="left", va="top")
    fig.text(
        x + 0.022, 0.972, title,
        fontsize=TITLE_FONT, fontweight="bold", ha="left", va="top",
    )


def shared_panel_label_x(inset_boxes):
    """Figure-x that letters a and b must share: the first inset's left edge."""
    if not inset_boxes:
        raise ValueError("inset_boxes must contain at least one (x0, y0, width, height)")
    x0 = float(inset_boxes[0][0])
    if not math.isfinite(x0):
        raise ValueError("inset left edge must be finite")
    return x0


def inset_title_line_y(inset_boxes, pad=0.012):
    """One figure-y for every inset letter and title, above the tallest frame."""
    if not inset_boxes:
        raise ValueError("inset_boxes must contain at least one (x0, y0, width, height)")
    tops = []
    for box in inset_boxes:
        if len(box) != 4:
            raise ValueError("each box must be (x0, y0, width, height)")
        _x0, y0, _width, height = (float(v) for v in box)
        tops.append(y0 + height)
    return max(tops) + float(pad)


def top_align_axes(axes):
    """Give every axes the same top edge so inset titles can share one row."""
    axes = list(axes)
    if not axes:
        return []
    axes[0].figure.canvas.draw()
    positions = [ax.get_position() for ax in axes]
    top = max(pos.y1 for pos in positions)
    aligned = []
    for ax, pos in zip(axes, positions):
        ax.set_position([pos.x0, top - pos.height, pos.width, pos.height])
        aligned.append(ax.get_position().bounds)
    return aligned


def inset_caption(letter, title):
    """Single-line inset label: ``b. Gulf Coast``."""
    return f"{letter}. {title}"


def label_inset_row(fig, inset_axes, letters, titles, pad=0.010):
    """Place b/c/d and their titles as one string on one figure-y.

    Returns the shared left x (letter a must use this) and the title-row y.
    """
    if len(inset_axes) != len(letters) or len(letters) != len(titles):
        raise ValueError("inset_axes, letters, and titles must have the same length")
    boxes = top_align_axes(inset_axes)
    x_shared = shared_panel_label_x(boxes)
    y_line = inset_title_line_y(boxes, pad=pad)
    for (x0, _y0, _w, _h), letter, title in zip(boxes, letters, titles):
        fig.text(
            x0, y_line, inset_caption(letter, title),
            fontsize=INSET_TITLE_FONT, fontweight="bold", ha="left", va="bottom",
        )
    return x_shared, y_line


def prepare_legend_ax(ax):
    """Empty, full-width strip under the globe. Both plates use this canvas."""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_navigate(False)
    return ax


def discrete_ramp_origin(n, chip_w, cx=0.5):
    """Left edge of a centered n-chip ramp. Exposed so tests can check centering."""
    if n <= 0:
        raise ValueError("n must be positive")
    if chip_w <= 0:
        raise ValueError("chip_w must be positive")
    return cx - n * chip_w / 2.0


def draw_discrete_ramp(ax, colors, labels, *, cx=0.5, y=0.30, chip_w=0.048, chip_h=0.20):
    """Nature-style adjacent chips with a label under every chip.

    Title is drawn by the caller *above* this ramp, never at the chip y.
    Returns (x0, x1, y, y + chip_h) in axes coordinates.
    """
    colors = list(colors)
    labels = list(labels)
    if len(colors) != len(labels):
        raise ValueError("colors and labels must have the same length")
    n = len(colors)
    x0 = discrete_ramp_origin(n, chip_w, cx)
    for i, (color, label) in enumerate(zip(colors, labels)):
        ax.add_patch(Rectangle(
            (x0 + i * chip_w, y), chip_w, chip_h, transform=ax.transAxes,
            facecolor=color, edgecolor="#666666", linewidth=0.32, clip_on=False,
        ))
        ax.text(
            x0 + (i + 0.5) * chip_w, y - 0.045, label,
            transform=ax.transAxes, fontsize=LEGEND_ITEM_FONT,
            ha="center", va="top", color=BLACK, clip_on=False,
        )
    return x0, x0 + n * chip_w, y, y + chip_h


def draw_line_keys(ax, colors, names, widths=None, *, cx=0.5, y=0.78, item_w=0.175):
    """Centered row of class line samples. Title is drawn by the caller above."""
    colors = list(colors)
    names = list(names)
    if len(colors) != len(names):
        raise ValueError("colors and names must have the same length")
    n = len(colors)
    x0 = discrete_ramp_origin(n, item_w, cx)
    for i, (color, name) in enumerate(zip(colors, names)):
        x = x0 + i * item_w
        lw = 2.2 if widths is None else max(float(widths[i]) * 1.6, 1.1)
        rgb = mpl.colors.to_rgb(color)
        lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        line, = ax.plot(
            [x, x + 0.028], [y, y], transform=ax.transAxes,
            color=color, lw=lw, solid_capstyle="butt", clip_on=False,
        )
        if lum > 0.80:
            line.set_path_effects([
                path_effects.Stroke(linewidth=lw + 0.85, foreground="#777777"),
                path_effects.Normal(),
            ])
        ax.text(
            x + 0.034, y, name, transform=ax.transAxes,
            fontsize=LEGEND_ITEM_FONT, ha="left", va="center", color=BLACK,
            clip_on=False,
        )
    return x0, x0 + n * item_w


def draw_legend_notes(ax, notes, *, cx=0.5, y=0.14, gap=0.038):
    """Centered second-row notes (hatch / swatch / line). Each note is a dict.

    Required keys: kind ('hatch'|'swatch'|'line'), label.
    Optional: color, linewidth.
    """
    if not notes:
        return 0.5, 0.5
    widths = []
    for note in notes:
        widths.append(0.024 + 0.006 + 0.0066 * len(note["label"]))
    total = sum(widths) + gap * (len(notes) - 1)
    x = cx - total / 2.0
    for note, width in zip(notes, widths):
        kind = note["kind"]
        color = note.get("color", "#888888")
        if kind == "hatch":
            ax.add_patch(Rectangle(
                (x, y - 0.07), 0.022, 0.14, transform=ax.transAxes,
                facecolor=note.get("facecolor", "#F4F4F4"), edgecolor="#222222",
                linewidth=0.3, hatch="//////", clip_on=False,
            ))
        elif kind == "swatch":
            ax.add_patch(Rectangle(
                (x, y - 0.07), 0.022, 0.14, transform=ax.transAxes,
                facecolor=color, edgecolor="#666666", linewidth=0.32, clip_on=False,
            ))
        elif kind == "line":
            ax.plot(
                [x, x + 0.028], [y, y], transform=ax.transAxes,
                color=color, lw=note.get("linewidth", 2.1),
                solid_capstyle="butt", clip_on=False,
            )
        else:
            raise ValueError(f"unknown note kind: {kind}")
        label_x = x + (0.034 if kind == "line" else 0.028)
        ax.text(
            label_x, y, note["label"], transform=ax.transAxes,
            fontsize=LEGEND_ITEM_FONT, ha="left", va="center", color=BLACK,
            clip_on=False,
        )
        x += width + gap
    return cx - total / 2.0, cx + total / 2.0


def save_plate(fig, output_stem):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{output_stem}.pdf", bbox_inches=None, dpi=450, facecolor="white")
        fig.savefig(f"{output_stem}.png", bbox_inches=None, dpi=450, facecolor="white")
    plt.close(fig)
    return Path(f"{output_stem}.pdf"), Path(f"{output_stem}.png")
