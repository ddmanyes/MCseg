"""
_panel_style.py
===============
Shared drawing helpers for the manuscript figure panels.

Font sizes here are tuned for panels that get placed into the merged Keynote
figures (fig_total.key) at roughly half-page width — a panel rendered at 183 mm
and then scaled down needs noticeably larger type than a standalone figure.
"""

from __future__ import annotations

from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

VISIUM_UM_PX = 0.2737   # Visium HD fullres (µm/px)


def add_scale_bar(ax: Axes, img_w_px: int, fontsize: int = 12,
                  scale_um: int = 50, um_per_px: float = VISIUM_UM_PX,
                  bar_color: str = "#111111", label: str | None = None,
                  corner: str = "left") -> None:
    """Draw a bottom-corner scale bar on an opaque white plate.

    The plate keeps the bar legible over any background (dark cell outlines,
    bright eosin) and visually matches the legend box in the opposite corner.
    Geometry is in axes fractions, so the plate stays a constant size
    regardless of the image aspect ratio.

    `label` overrides the default "{scale_um} µm" caption (e.g. "1 mm").
    `corner` is "left" or "right" (bottom-left / bottom-right placement).
    """
    bar_frac = (scale_um / um_per_px) / img_w_px   # bar length, axes fraction
    y0 = 0.055                                     # bar height, axes fraction
    pad_x, plate_h = 0.012, 0.145
    x0 = 0.022 if corner == "left" else 1.0 - 0.022 - bar_frac

    ax.add_patch(Rectangle((x0 - pad_x, y0 - 0.035), bar_frac + 2 * pad_x, plate_h,
                           transform=ax.transAxes, facecolor="white",
                           edgecolor="#bbbbbb", linewidth=0.8, alpha=0.92,
                           zorder=9, clip_on=False))
    ax.plot([x0, x0 + bar_frac], [y0, y0], transform=ax.transAxes,
            color=bar_color, linewidth=3.5, solid_capstyle="butt", zorder=10)
    ax.text(x0 + bar_frac / 2, y0 + 0.028, label or f"{scale_um} µm",
            transform=ax.transAxes, color=bar_color, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold", zorder=10)
