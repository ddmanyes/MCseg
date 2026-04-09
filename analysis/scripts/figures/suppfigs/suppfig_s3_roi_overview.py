"""
gen_suppfig_s3_roi_overview.py
==============================
Supp. Fig. S3 — ROI location overview
2-panel: Left = LUAD (6 ROIs), Right = CRC (14 ROIs + TLS)

Output: manuscript/supplementary/SuppFigS3.png
"""

from __future__ import annotations

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from pathlib import Path
from PIL import Image

# ── paths ─────────────────────────────────────────────────────────────────

MANUSCRIPT_DIR = Path("/Volumes/SSD/plan_a/manuscript")
OUT_DIR        = MANUSCRIPT_DIR / "supplementary"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LUAD_HIRES = Path("/Volumes/SSD/plan_a/tissue sample/LUAD/visium/binned_outputs/square_002um/spatial/tissue_hires_image.png")
CRC_HIRES  = Path("/Volumes/SSD/plan_a/tissue sample/CRC/visium/official_v4/binned_outputs/binned_outputs/square_002um/spatial/tissue_hires_image.png")

LUAD_STATE   = Path("/Volumes/SSD/plan_a/xenium_he_seg/config/state.json")
CRC_ROI_INFO = Path("/Volumes/SSD/plan_a/crc_he_seg/results/rois/roi_info.json")

# scale factors (fullres → hires)
LUAD_SCALEF = 0.1386642
CRC_SCALEF  = 0.07973422

# µm per fullres pixel
PIXEL_UM = 0.2737

# ── LUAD ROI metadata ─────────────────────────────────────────────────────

LUAD_ROI_META = {
    "roi1": dict(label="ROI 1", desc="Tumor boundary",       color="#e41a1c"),
    "2":    dict(label="ROI 2", desc="Tumor stroma",          color="#377eb8"),
    "3":    dict(label="ROI 3", desc="Mixed (tumor+stroma)",  color="#4daf4a"),
    "4":    dict(label="ROI 4", desc="Normal–tumor interface", color="#ff7f00"),
    "5":    dict(label="ROI 5", desc="Alveolar region",       color="#984ea3"),
    "6":    dict(label="ROI 6", desc="Tumor core",            color="#a65628"),
}

# ── CRC ROI metadata (roi1–roi14 for TAS, roi15 = TLS) ────────────────────

import matplotlib
CRC_COLORS = [matplotlib.colormaps["tab20"](i / 20) for i in range(0, 28, 2)]  # 14 distinct colours
CRC_TAS_ROIS = [f"roi{i}" for i in range(1, 15)]  # roi1 … roi14

# ── style ──────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":    ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":      7,
    "savefig.dpi":    300,
    "savefig.bbox":   "tight",
    "savefig.facecolor": "white",
})

MM_TO_IN = 1 / 25.4


def load_hires(path: Path) -> np.ndarray:
    img = np.array(Image.open(path).convert("RGB"))
    return img


def draw_roi_box(ax, x_full, y_full, w_full, h_full, scalef, color, label,
                 label_side="right", label_offset=(4, 0)):
    """Draw a ROI rectangle on a hires axes. Coordinates in fullres px."""
    x = x_full * scalef
    y = y_full * scalef
    w = w_full * scalef
    h = h_full * scalef

    rect = Rectangle((x, y), w, h,
                      linewidth=1.2, edgecolor=color, facecolor="none",
                      zorder=5)
    ax.add_patch(rect)

    # Determine label position
    if label_side == "right":
        lx = x + w + label_offset[0]
        ly = y + h / 2 + label_offset[1]
        ha = "left"
    elif label_side == "left":
        lx = x - label_offset[0]
        ly = y + h / 2 + label_offset[1]
        ha = "right"
    elif label_side == "top":
        lx = x + w / 2 + label_offset[0]
        ly = y - label_offset[1]
        ha = "center"
    else:  # bottom
        lx = x + w / 2 + label_offset[0]
        ly = y + h + label_offset[1]
        ha = "center"

    ax.text(lx, ly, label, color=color, fontsize=6.5, fontweight="bold",
            ha=ha, va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                      alpha=0.75, linewidth=0))


def add_scale_bar(ax, img_w_full, img_h_full, scalef, pixel_um,
                  bar_um=2000, margin=0.03):
    """Add a scale bar (bottom-right)."""
    bar_px_full = bar_um / pixel_um
    bar_px_hires = bar_px_full * scalef
    x0 = img_w_full * scalef * (1 - margin) - bar_px_hires
    y0 = img_h_full * scalef * (1 - margin)

    ax.plot([x0, x0 + bar_px_hires], [y0, y0],
            color="white", linewidth=2.5, solid_capstyle="butt", zorder=10)
    ax.text(x0 + bar_px_hires / 2, y0 - img_h_full * scalef * 0.012,
            f"{bar_um // 1000} mm",
            ha="center", va="bottom", fontsize=6.5, color="white",
            fontweight="bold", zorder=10)


def main():
    # ── load images ──────────────────────────────────────────────
    luad_img = load_hires(LUAD_HIRES)   # (H, W, 3) = (6000, 3705)
    crc_img  = load_hires(CRC_HIRES)    # (H, W, 3) = (3886, 6000)

    luad_H, luad_W = luad_img.shape[:2]   # 6000, 3705
    crc_H,  crc_W  = crc_img.shape[:2]    # 3886, 6000

    # fullres dimensions
    luad_W_full = luad_W / LUAD_SCALEF
    luad_H_full = luad_H / LUAD_SCALEF
    crc_W_full  = crc_W  / CRC_SCALEF
    crc_H_full  = crc_H  / CRC_SCALEF

    # ── load ROI coords ──────────────────────────────────────────
    with open(LUAD_STATE) as f:
        luad_state = json.load(f)
    luad_rois = {r["name"]: r for r in luad_state["rois"]}

    with open(CRC_ROI_INFO) as f:
        crc_roi_info = json.load(f)

    # ── figure layout ────────────────────────────────────────────
    # LUAD is portrait (3705×6000 hires); CRC is landscape (6000×3886 hires)
    # Use gridspec with width_ratio proportional to W
    luad_aspect = luad_W / luad_H   # ~0.617
    crc_aspect  = crc_W  / crc_H    # ~1.544

    total_w_in = 183 * MM_TO_IN   # double column
    luad_w_in  = total_w_in * luad_W / (luad_W + crc_W * (luad_H / crc_H))
    # Simpler: give each panel equal height, width proportional to aspect
    panel_h_in = 90 * MM_TO_IN
    luad_w_in  = panel_h_in * luad_aspect
    crc_w_in   = panel_h_in * crc_aspect
    total_w_in = luad_w_in + crc_w_in + 10 * MM_TO_IN  # +gap

    fig = plt.figure(figsize=(total_w_in, panel_h_in + 8 * MM_TO_IN))
    gs  = fig.add_gridspec(1, 2,
                           width_ratios=[luad_w_in, crc_w_in],
                           wspace=0.04,
                           left=0.01, right=0.99,
                           top=0.93, bottom=0.01)
    ax_luad = fig.add_subplot(gs[0])
    ax_crc  = fig.add_subplot(gs[1])

    # ── LUAD panel ───────────────────────────────────────────────
    ax_luad.imshow(luad_img, origin="upper")
    ax_luad.set_xlim(0, luad_W)
    ax_luad.set_ylim(luad_H, 0)
    ax_luad.axis("off")
    ax_luad.set_title("LUAD — Lung Adenocarcinoma", fontsize=8,
                      fontweight="bold", pad=4)
    ax_luad.text(-0.02, 1.04, "a", transform=ax_luad.transAxes,
                 fontsize=9, fontweight="bold", va="top", ha="right")

    # Label sides to avoid crowding (manual tuning)
    luad_label_side = {
        "roi1": ("right",  (4, 0)),
        "2":    ("right",  (4, 0)),
        "3":    ("left",   (4, 0)),
        "4":    ("left",   (4, 0)),
        "5":    ("right",  (4, 0)),
        "6":    ("right",  (4, 0)),
    }

    for name, meta in LUAD_ROI_META.items():
        if name not in luad_rois:
            continue
        r    = luad_rois[name]
        side, offset = luad_label_side.get(name, ("right", (4, 0)))
        draw_roi_box(ax_luad,
                     r["x"], r["y"], r["width_px"], r["height_px"],
                     LUAD_SCALEF, meta["color"], meta["label"],
                     label_side=side, label_offset=offset)

    add_scale_bar(ax_luad, luad_W_full, luad_H_full, LUAD_SCALEF, PIXEL_UM,
                  bar_um=2000)

    # LUAD legend
    luad_handles = [
        mpatches.Patch(color=m["color"], label=f'{m["label"]} – {m["desc"]}')
        for m in LUAD_ROI_META.values()
    ]
    ax_luad.legend(handles=luad_handles, loc="upper left",
                   fontsize=5.5, framealpha=0.9, frameon=True,
                   edgecolor="#ccc", handlelength=1.2, borderpad=0.5)

    # ── CRC panel ────────────────────────────────────────────────
    ax_crc.imshow(crc_img, origin="upper")
    ax_crc.set_xlim(0, crc_W)
    ax_crc.set_ylim(crc_H, 0)
    ax_crc.axis("off")
    ax_crc.set_title("CRC — Colorectal Cancer", fontsize=8,
                     fontweight="bold", pad=4)
    ax_crc.text(-0.02, 1.04, "b", transform=ax_crc.transAxes,
                fontsize=9, fontweight="bold", va="top", ha="right")

    crc_handles = []
    for i, roi_name in enumerate(CRC_TAS_ROIS):
        if roi_name not in crc_roi_info:
            continue
        r     = crc_roi_info[roi_name]
        color = CRC_COLORS[i]
        x0, y0 = r["x0"], r["y0"]
        w_full  = r["x1"] - r["x0"]
        h_full  = r["y1"] - r["y0"]

        # Decide label side based on position in image
        x_hires = x0 * CRC_SCALEF
        side = "right" if x_hires < crc_W * 0.6 else "left"

        draw_roi_box(ax_crc, x0, y0, w_full, h_full,
                     CRC_SCALEF, color, roi_name.upper(),
                     label_side=side, label_offset=(4, 0))
        crc_handles.append(mpatches.Patch(color=color, label=roi_name.upper()))

    # TLS roi15 — special marker (star outline)
    if "roi15" in crc_roi_info:
        r15    = crc_roi_info["roi15"]
        tls_color = "#d62728"
        x0, y0 = r15["x0"], r15["y0"]
        w_full  = r15["x1"] - r15["x0"]
        h_full  = r15["y1"] - r15["y0"]
        x_h = (x0 + w_full / 2) * CRC_SCALEF
        y_h = (y0 + h_full / 2) * CRC_SCALEF

        rect = Rectangle((x0 * CRC_SCALEF, y0 * CRC_SCALEF),
                          w_full * CRC_SCALEF, h_full * CRC_SCALEF,
                          linewidth=1.5, edgecolor=tls_color,
                          facecolor="none", linestyle="--", zorder=5)
        ax_crc.add_patch(rect)
        ax_crc.text(x_h, y0 * CRC_SCALEF - 6, "TLS",
                    color=tls_color, fontsize=6.5, fontweight="bold",
                    ha="center", va="bottom", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec="none", alpha=0.75))
        crc_handles.append(mpatches.Patch(color=tls_color, label="TLS (ROI 15)",
                                          linestyle="--",
                                          fill=False, edgecolor=tls_color))

    add_scale_bar(ax_crc, crc_W_full, crc_H_full, CRC_SCALEF, PIXEL_UM,
                  bar_um=5000)

    ax_crc.legend(handles=crc_handles, loc="lower left",
                  fontsize=5.0, framealpha=0.85, frameon=True,
                  edgecolor="#ccc", handlelength=1.0, borderpad=0.4,
                  ncol=3)

    # ── save ─────────────────────────────────────────────────────
    out_path = OUT_DIR / "SuppFigS3.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out_path}")

    # Remove old S3a/S3b if we are replacing them
    for old in ["SuppFigS3a.png", "SuppFigS3b.png"]:
        p = OUT_DIR / old
        if p.exists():
            p.unlink()
            print(f"  removed {old}")


if __name__ == "__main__":
    main()
