"""
ENACT CRC benchmark overlay figure
- H&E image background
- MCseg cell contours coloured by CellTypist broad_label (Epithelial / Stromal / Immune)
- GT centroids coloured by gt_label (matched = solid circle, unmatched = open circle)
Zoom into a representative region containing all three cell types.
"""
import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_dilation
from pathlib import Path

RESULT_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_crc_f1")
OUT_PATH    = RESULT_DIR / "fig_celltype_overlay.png"

CROP_X0 = 5154   # ENACT local coord origin of he_crop (col)
CROP_Y0 = 4635   # BTF row origin of he_crop (row)

TYPE_COLORS = {
    "epithelial cells": "#E74C3C",
    "stromal cells":    "#2ECC71",
    "immune cells":     "#3498DB",
}
DEFAULT_COLOR = "#AAAAAA"

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading H&E crop…")
he = tifffile.imread(str(RESULT_DIR / "he_crop.tif"))
H, W = he.shape[:2]

print("Loading MCseg mask…")
mask = np.load(str(RESULT_DIR / "mcseg_mask.npy")).astype(np.int32)

print("Loading CellTypist labels…")
ct = pd.read_csv(RESULT_DIR / "celltypist_labels.csv")
id2type = dict(zip(ct["cell_id"].values, ct["broad_label"].str.lower().values))

print("Loading GT matched…")
gt = pd.read_csv(RESULT_DIR / "gt_matched.csv")
gt = gt.dropna(subset=["cell_x", "cell_y"])
gt["col_local"] = (gt["cell_x"] - CROP_X0).astype(int)
gt["row_local"] = (gt["cell_y"] - CROP_Y0).astype(int)

# ── Find zoom region with best mix of all three cell types ─────────────────────
print("Finding representative zoom region…")
ZOOM_SIZE = 800   # smaller zoom = larger apparent cells at same DPI

best_score = -1
best_r0, best_c0 = 0, 0
step = 200
for r0 in range(0, H - ZOOM_SIZE, step):
    for c0 in range(0, W - ZOOM_SIZE, step):
        sub = mask[r0:r0+ZOOM_SIZE, c0:c0+ZOOM_SIZE]
        cell_ids = np.unique(sub[sub > 0])
        n_cells = len(cell_ids)
        if n_cells < 30:
            continue
        types = [id2type.get(int(cid), "") for cid in cell_ids]
        n_epi = sum(1 for t in types if t == "epithelial cells")
        n_str = sum(1 for t in types if t == "stromal cells")
        n_imm = sum(1 for t in types if t == "immune cells")
        if n_epi == 0 or n_str == 0 or n_imm == 0:
            continue
        # favour density + balance
        score = min(n_epi, n_str, n_imm) * 3 + n_cells
        if score > best_score:
            best_score = score
            best_r0, best_c0 = r0, c0

r0, c0 = best_r0, best_c0
r1, c1 = r0 + ZOOM_SIZE, c0 + ZOOM_SIZE
print(f"Zoom region: rows {r0}:{r1}, cols {c0}:{c1}  (score={best_score})")

he_z = he[r0:r1, c0:c1]

# Read mask with margin so boundary cells get complete contours
MARGIN = 150
mr0 = max(0, r0 - MARGIN)
mc0 = max(0, c0 - MARGIN)
mr1 = min(H, r1 + MARGIN)
mc1 = min(W, c1 + MARGIN)
mask_padded = mask[mr0:mr1, mc0:mc1]

# ── Build contour RGBA layer on padded region, then crop ───────────────────────
print("Computing cell contours…")
boundary_padded = find_boundaries(mask_padded, mode="thick")

# offset to get back to zoom coords
off_r = r0 - mr0
off_c = c0 - mc0

contour_rgba = np.zeros((ZOOM_SIZE, ZOOM_SIZE, 4), dtype=np.float32)
struct = np.ones((3, 3), dtype=bool)   # 3×3 dilation → ~2px thick contour
for cid in np.unique(mask_padded[mask_padded > 0]):
    ctype = id2type.get(int(cid), "")
    hex_c = TYPE_COLORS.get(ctype, DEFAULT_COLOR)
    r_c = int(hex_c[1:3], 16) / 255
    g_c = int(hex_c[3:5], 16) / 255
    b_c = int(hex_c[5:7], 16) / 255
    edge_pad = binary_dilation(boundary_padded & (mask_padded == cid), structure=struct)
    # crop back to display area
    edge_z = edge_pad[off_r:off_r+ZOOM_SIZE, off_c:off_c+ZOOM_SIZE]
    contour_rgba[edge_z] = [r_c, g_c, b_c, 1.0]

mask_z = mask_padded[off_r:off_r+ZOOM_SIZE, off_c:off_c+ZOOM_SIZE]

# ── GT centroids in zoom ───────────────────────────────────────────────────────
in_zoom = gt[
    (gt["row_local"] >= r0) & (gt["row_local"] < r1) &
    (gt["col_local"] >= c0) & (gt["col_local"] < c1)
].copy()
in_zoom["row_z"] = in_zoom["row_local"] - r0
in_zoom["col_z"] = in_zoom["col_local"] - c0

# ── Plot ───────────────────────────────────────────────────────────────────────
print("Rendering figure…")
fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=200)

for ax, show_gt, title in zip(
    axes,
    [False, True],
    ["MCseg cell-type contours", "MCseg contours + GT centroids"],
):
    ax.imshow(he_z, interpolation="nearest")
    ax.imshow(contour_rgba, interpolation="nearest")

    if show_gt:
        for label, grp in in_zoom.groupby("gt_label"):
            col = TYPE_COLORS.get(label.lower(), DEFAULT_COLOR)
            m  = grp[grp["mcseg_cell_id"] > 0]
            um = grp[grp["mcseg_cell_id"] == 0]
            if len(m):
                ax.scatter(m["col_z"], m["row_z"], s=10, c=col,
                           marker="o", linewidths=0, alpha=0.9, zorder=5)
            if len(um):
                ax.scatter(um["col_z"], um["row_z"], s=14, c="none",
                           edgecolors=col, marker="o", linewidths=0.8,
                           alpha=0.7, zorder=5)

    ax.set_title(title, fontsize=11, pad=6)
    ax.axis("off")

legend_patches = [
    mpatches.Patch(color=c, label=lbl.capitalize())
    for lbl, c in TYPE_COLORS.items()
]
gt_solid = plt.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="grey", markersize=5,
                      label="GT matched (solid)")
gt_open  = plt.Line2D([0], [0], marker="o", color="grey",
                      markerfacecolor="none", markersize=5,
                      linewidth=0.8, label="GT unmatched (open)")

axes[0].legend(handles=legend_patches, fontsize=8, loc="lower right",
               framealpha=0.75, title="Contour colour")
axes[1].legend(handles=legend_patches + [gt_solid, gt_open],
               fontsize=8, loc="lower right", framealpha=0.75,
               title="Contour / GT")

fig.suptitle(
    f"ENACT CRC benchmark — MCseg cell-type contour overlay\n"
    f"zoom {ZOOM_SIZE}×{ZOOM_SIZE} px  |  rows {r0}–{r1}, cols {c0}–{c1}",
    fontsize=11, y=1.01,
)
plt.tight_layout()
fig.savefig(str(OUT_PATH), dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_PATH}")
plt.close()
