"""
ENACT Breast benchmark — qualitative generalisability overlay (SuppFigS12)
2-panel figure:
  Left:  H&E + MCseg cell contours coloured by CellTypist broad label
  Right: same + GT region boundaries (Tumour / Stroma dashed outlines)
"""
import json
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

RESULT_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_breast_f1")
GEOJSON    = Path("/Volumes/SSD/plan_a/tissue sample/ENACT_additional_samples"
                  "/human_breast/annotations"
                  "/Visium_HD_Human_Breast_Cancer_Fresh_Frozen-wsi-3588_1641_27772_23646.geojson")
OUT_S12    = RESULT_DIR / "fig_breast_overlay_s12.png"

CROP_X0 = 3588
CROP_Y0 = 1641

TYPE_COLORS = {
    "epithelial": "#E74C3C",
    "stromal":    "#2ECC71",
    "immune":     "#3498DB",
}
GT_COLORS = {
    "Tumor":  "#E74C3C",
    "Stroma": "#2ECC71",
}
DEFAULT_COLOR = "#AAAAAA"

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading H&E crop…")
he   = tifffile.imread(str(RESULT_DIR / "he_crop_breast.tif"))
H, W = he.shape[:2]

print("Loading MCseg mask…")
mask = np.load(str(RESULT_DIR / "mcseg_mask_breast.npy")).astype(np.int32)

print("Loading CellTypist labels…")
ct      = pd.read_csv(RESULT_DIR / "celltypist_labels_breast.csv")
id2type = dict(zip(ct["cell_id"].values, ct["broad_label"].str.lower().values))

print("Loading GT matched…")
gt         = pd.read_csv(RESULT_DIR / "gt_matched_breast.csv")
gt_matched = gt[gt["gt_label"].isin(["epithelial", "stromal"])].copy()

print("Loading GT GeoJSON polygons…")
with open(str(GEOJSON)) as f:
    gj = json.load(f)

def _geojson_to_crop_rings(feature) -> list[np.ndarray]:
    """Return list of (N,2) arrays [[col_crop, row_crop], ...] for each ring."""
    geom   = feature["geometry"]
    polys  = geom["coordinates"] if geom["type"] == "Polygon" else []
    if geom["type"] == "MultiPolygon":
        polys = [ring for mp in geom["coordinates"] for ring in mp]
    rings  = []
    for ring in polys:
        if len(ring) < 3:
            continue
        pts = np.array(ring, dtype=float)
        if pts.ndim == 3:       # some rings are [[x,y]] nested
            pts = pts[:, 0, :]
        # GeoJSON x=col, y=row; convert to crop-local coordinates
        col_crop = pts[:, 0] - CROP_X0
        row_crop = pts[:, 1] - CROP_Y0
        rings.append(np.column_stack([col_crop, row_crop]))
    return rings

gt_regions: dict[str, list[np.ndarray]] = {}
for feat in gj["features"]:
    name = feat["properties"]["classification"]["name"]   # "Tumor" or "Stroma"
    gt_regions[name] = _geojson_to_crop_rings(feat)

# ── Find zoom region ───────────────────────────────────────────────────────────
print("Finding representative zoom region…")
ZOOM_SIZE   = 800
best_score  = -1
best_r0 = best_c0 = 0
step = 200

for r0 in range(0, H - ZOOM_SIZE, step):
    for c0 in range(0, W - ZOOM_SIZE, step):
        sub   = gt_matched[
            (gt_matched["r"] >= r0) & (gt_matched["r"] < r0 + ZOOM_SIZE) &
            (gt_matched["c"] >= c0) & (gt_matched["c"] < c0 + ZOOM_SIZE)
        ]
        n_epi = (sub["gt_label"] == "epithelial").sum()
        n_str = (sub["gt_label"] == "stromal").sum()
        if n_epi < 10 or n_str < 10:
            continue
        score = min(n_epi, n_str) * 2 + len(sub)
        if score > best_score:
            best_score = score
            best_r0, best_c0 = r0, c0

r0, c0 = best_r0, best_c0
r1, c1 = r0 + ZOOM_SIZE, c0 + ZOOM_SIZE
print(f"Zoom: rows {r0}:{r1}, cols {c0}:{c1}  (score={best_score})")

he_z = he[r0:r1, c0:c1]

# ── Build MCseg contour RGBA ───────────────────────────────────────────────────
MARGIN = 150
mr0 = max(0, r0 - MARGIN);  mc0 = max(0, c0 - MARGIN)
mr1 = min(H, r1 + MARGIN);  mc1 = min(W, c1 + MARGIN)
mask_padded     = mask[mr0:mr1, mc0:mc1]
boundary_padded = find_boundaries(mask_padded, mode="thick")
off_r, off_c    = r0 - mr0, c0 - mc0

contour_rgba = np.zeros((ZOOM_SIZE, ZOOM_SIZE, 4), dtype=np.float32)
struct       = np.ones((3, 3), dtype=bool)
print("Computing contours…")
for cid in np.unique(mask_padded[mask_padded > 0]):
    ctype = id2type.get(int(cid), "")
    hex_c = TYPE_COLORS.get(ctype, DEFAULT_COLOR)
    rc, gc, bc = (int(hex_c[1:3], 16)/255,
                  int(hex_c[3:5], 16)/255,
                  int(hex_c[5:7], 16)/255)
    edge_pad = binary_dilation(boundary_padded & (mask_padded == cid), structure=struct)
    edge_z   = edge_pad[off_r:off_r+ZOOM_SIZE, off_c:off_c+ZOOM_SIZE]
    contour_rgba[edge_z] = [rc, gc, bc, 1.0]

# ── Plot ───────────────────────────────────────────────────────────────────────
print("Rendering figure…")
fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=200)

for ax, show_gt_boundary, title in zip(
    axes,
    [False, True],
    ["MCseg cell-type contours",
     "MCseg contours + GT region boundaries"],
):
    ax.imshow(he_z, interpolation="nearest")
    ax.imshow(contour_rgba, interpolation="nearest")

    if show_gt_boundary:
        for region_name, rings in gt_regions.items():
            col = GT_COLORS.get(region_name, "#FFFFFF")
            for ring in rings:
                # transform ring coords to zoom-local
                col_z = ring[:, 0] - c0
                row_z = ring[:, 1] - r0
                # only draw if ring intersects zoom window
                if (col_z.max() < 0 or col_z.min() > ZOOM_SIZE or
                        row_z.max() < 0 or row_z.min() > ZOOM_SIZE):
                    continue
                ax.plot(col_z, row_z, color=col, linewidth=1.2,
                        linestyle="--", alpha=0.85, zorder=6)

    ax.set_title(title, fontsize=11, pad=6)
    ax.axis("off")

# Legend
legend_handles = [
    mpatches.Patch(color=c, label=lbl.capitalize())
    for lbl, c in TYPE_COLORS.items()
] + [
    plt.Line2D([0], [0], color=GT_COLORS["Tumor"],  linewidth=1.5,
               linestyle="--", label="GT: Tumour boundary"),
    plt.Line2D([0], [0], color=GT_COLORS["Stroma"], linewidth=1.5,
               linestyle="--", label="GT: Stroma boundary"),
]
axes[0].legend(
    handles=[mpatches.Patch(color=c, label=lbl.capitalize())
             for lbl, c in TYPE_COLORS.items()],
    fontsize=8, loc="lower right", framealpha=0.75, title="Contour colour",
)
axes[1].legend(
    handles=legend_handles,
    fontsize=8, loc="lower right", framealpha=0.75,
)

fig.suptitle(
    "MCseg generalisation to breast cancer (Visium HD Human Breast Cancer Fresh Frozen)\n"
    f"Representative region {ZOOM_SIZE}×{ZOOM_SIZE} px — no parameter re-tuning",
    fontsize=11, y=1.01,
)
plt.tight_layout()
fig.savefig(str(OUT_S12), dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_S12}")
plt.close()
