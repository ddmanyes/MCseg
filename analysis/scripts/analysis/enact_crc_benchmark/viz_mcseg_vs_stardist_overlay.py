"""
MCseg vs StarDist contour comparison overlay
Left:  H&E + StarDist nuclear contours (from ENACT cells_df.csv WKT polygons)
Right: H&E + MCseg whole-cell contours (from mcseg_mask.npy)
Same zoom region, single colour per method.
"""
import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_dilation
from shapely import wkt
from pathlib import Path

RESULT_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_crc_f1")
ENACT_CRC  = Path("/Volumes/SSD/plan_a/tissue sample/ENACT_supporting_files/public_data/human_colorectal")
OUT_PATH   = RESULT_DIR / "fig_mcseg_vs_stardist_overlay.png"

CROP_X0 = 5154   # ENACT local x → he_crop col: col = cell_x - CROP_X0
CROP_Y0 = 4635   # BTF row       → he_crop row: row = cell_y - CROP_Y0

MCSEG_COLOR    = "#E74C3C"   # red
STARDIST_COLOR = "#00FFFF"   # cyan — high contrast on pink H&E

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading H&E crop…")
he = tifffile.imread(str(RESULT_DIR / "he_crop.tif"))
H, W = he.shape[:2]

print("Loading MCseg mask…")
mask = np.load(str(RESULT_DIR / "mcseg_mask.npy")).astype(np.int32)

print("Loading StarDist cells_df…")
cells_df = pd.read_csv(ENACT_CRC / "paper_results" / "cells_df.csv")
cells_df["row_crop"] = cells_df["cell_y"] - CROP_Y0
cells_df["col_crop"] = cells_df["cell_x"] - CROP_X0

# ── Find zoom region: dense cells, both methods present ───────────────────────
print("Finding zoom region…")
ZOOM_SIZE = 800
MARGIN    = 150

best_score = -1
best_r0, best_c0 = 0, 0
step = 200
for r0 in range(0, H - ZOOM_SIZE, step):
    for c0 in range(0, W - ZOOM_SIZE, step):
        sub = mask[r0:r0+ZOOM_SIZE, c0:c0+ZOOM_SIZE]
        n_mcseg = len(np.unique(sub[sub > 0]))
        in_z = cells_df[
            (cells_df["row_crop"] >= r0) & (cells_df["row_crop"] < r0+ZOOM_SIZE) &
            (cells_df["col_crop"] >= c0) & (cells_df["col_crop"] < c0+ZOOM_SIZE)
        ]
        n_sd = len(in_z)
        if n_mcseg < 30 or n_sd < 20:
            continue
        score = n_mcseg + n_sd
        if score > best_score:
            best_score = score
            best_r0, best_c0 = r0, c0

r0, c0 = best_r0, best_c0
r1, c1 = r0 + ZOOM_SIZE, c0 + ZOOM_SIZE
print(f"Zoom: rows {r0}:{r1}, cols {c0}:{c1}  (score={best_score})")

he_z = he[r0:r1, c0:c1]

# ── MCseg contours (with margin for complete boundaries) ───────────────────────
print("Computing MCseg contours…")
mr0 = max(0, r0 - MARGIN); mc0 = max(0, c0 - MARGIN)
mr1 = min(H, r1 + MARGIN); mc1 = min(W, c1 + MARGIN)
mask_padded   = mask[mr0:mr1, mc0:mc1]
boundary_pad  = find_boundaries(mask_padded, mode="thick")
off_r, off_c  = r0 - mr0, c0 - mc0

struct = np.ones((3, 3), dtype=bool)
mcseg_rgba = np.zeros((ZOOM_SIZE, ZOOM_SIZE, 4), dtype=np.float32)
r_c = int(MCSEG_COLOR[1:3], 16) / 255
g_c = int(MCSEG_COLOR[3:5], 16) / 255
b_c = int(MCSEG_COLOR[5:7], 16) / 255
for cid in np.unique(mask_padded[mask_padded > 0]):
    edge_pad = binary_dilation(boundary_pad & (mask_padded == cid), structure=struct)
    edge_z   = edge_pad[off_r:off_r+ZOOM_SIZE, off_c:off_c+ZOOM_SIZE]
    mcseg_rgba[edge_z] = [r_c, g_c, b_c, 1.0]

# ── StarDist polygon patches ───────────────────────────────────────────────────
print("Building StarDist polygon patches…")
sd_in = cells_df[
    (cells_df["row_crop"] >= r0 - MARGIN) & (cells_df["row_crop"] < r1 + MARGIN) &
    (cells_df["col_crop"] >= c0 - MARGIN) & (cells_df["col_crop"] < c1 + MARGIN)
]

sd_patches = []
for geom_wkt in sd_in["geometry"]:
    try:
        geom = wkt.loads(geom_wkt)
    except Exception:
        continue
    if geom.geom_type != "Polygon":
        continue
    coords   = np.array(geom.exterior.coords)
    col_arr  = coords[:, 0] - CROP_X0 - c0
    row_arr  = coords[:, 1] - CROP_Y0 - r0
    sd_patches.append(MplPolygon(np.column_stack([col_arr, row_arr]), closed=True))

print(f"  StarDist patches in view: {len(sd_patches)}")

# ── Plot ───────────────────────────────────────────────────────────────────────
print("Rendering…")
fig, axes = plt.subplots(1, 2, figsize=(16, 8), dpi=200)

# Left — StarDist
ax = axes[0]
ax.imshow(he_z, interpolation="nearest")
pc = PatchCollection(sd_patches, edgecolors=STARDIST_COLOR, facecolors="none",
                     linewidths=1.5, alpha=1.0)
ax.add_collection(pc)
ax.set_xlim(0, ZOOM_SIZE)
ax.set_ylim(ZOOM_SIZE, 0)
n_sd_z = len(cells_df[
    (cells_df["row_crop"] >= r0) & (cells_df["row_crop"] < r1) &
    (cells_df["col_crop"] >= c0) & (cells_df["col_crop"] < c1)
])
ax.set_title("StarDist (ENACT) — nuclear boundaries", fontsize=11, pad=6)
ax.axis("off")
ax.legend(handles=[mpatches.Patch(edgecolor=STARDIST_COLOR, facecolor="none",
                                   label=f"StarDist  n={n_sd_z:,}")],
          fontsize=9, loc="lower right", framealpha=0.75)

# Right — MCseg
ax = axes[1]
ax.imshow(he_z, interpolation="nearest")
ax.imshow(mcseg_rgba, interpolation="nearest")
n_mc_z = len(np.unique(mask[r0:r1, c0:c1])) - 1
ax.set_title("MCseg v2 — whole-cell boundaries (Voronoi)", fontsize=11, pad=6)
ax.axis("off")
ax.legend(handles=[mpatches.Patch(color=MCSEG_COLOR,
                                   label=f"MCseg  n={n_mc_z:,}")],
          fontsize=9, loc="lower right", framealpha=0.75)

fig.suptitle(
    "MCseg vs StarDist: segmentation boundary comparison — CRC H&E\n"
    f"zoom {ZOOM_SIZE}×{ZOOM_SIZE} px  |  he_crop rows {r0}–{r1}, cols {c0}–{c1}",
    fontsize=11, y=1.01,
)
plt.tight_layout()
fig.savefig(str(OUT_PATH), dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_PATH}")
plt.close()
