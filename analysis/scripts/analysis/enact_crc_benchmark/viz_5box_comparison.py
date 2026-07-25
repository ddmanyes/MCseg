"""
viz_5box_comparison.py
----------------------
5-region × 5-method H&E overlay comparison.
Generates a 5×5 grid figure: rows = 5 ROI boxes, cols = 5 methods.

Usage:
    cd /Volumes/SSD/plan_a
    uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/viz_5box_comparison.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries

PLAN_A     = Path("/Volumes/SSD/plan_a")
F1_DIR     = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
RESULT_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
HE_CROP    = F1_DIR / "he_crop.tif"
OUT_FIG    = RESULT_DIR / "overlays" / "fig_5box_comparison.png"

METHODS = {
    "StarDist": {"mask": RESULT_DIR / "stardist_mask.npy",       "color": "#FFD700"},
    "MCseg":    {"mask": F1_DIR     / "mcseg_mask_7pass.npy",    "color": "#FFD700"},
    "SR":       {"mask": RESULT_DIR / "sr_mask.npy",             "color": "#FFD700"},
    "NUC":      {"mask": RESULT_DIR / "cellpose_nuc_mask.npy",   "color": "#FFD700"},
    "ProSeg":   {"mask": RESULT_DIR / "proseg_mask.npy",         "color": "#FFD700"},
}

PATCH_SIZE = 1024   # px per box

# 5 boxes: (col_center, row_center, label)
# Estimated from overview image visual inspection
BOXES = [
    (3236,  1523, "Box 1 — Mucosa / lumen edge"),
    (3570,  4937, "Box 2 — Stroma / nerve"),
    (2805,  7450, "Box 3 — Tumor (lower-left)"),
    (7259,  7450, "Box 4 — Glands (lower-right)"),
    (1640, 10227, "Box 5 — Vessels / necrosis"),
]


def hex_to_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def overlay_patch(he: np.ndarray, mask: np.ndarray, color: str,
                  boundary_width: int = 3) -> np.ndarray:
    from scipy.ndimage import binary_dilation
    bnd = find_boundaries(mask, mode="inner")
    if boundary_width > 1:
        bnd = binary_dilation(bnd, iterations=boundary_width - 1)
    img = he.astype(np.float32) / 255.0
    img[bnd] = hex_to_rgb(color)
    return np.clip(img, 0, 1)


def get_patch(arr: np.ndarray, cy: int, cx: int, H: int, W: int):
    half = PATCH_SIZE // 2
    r0 = max(0, cy - half)
    r1 = min(H, r0 + PATCH_SIZE)
    c0 = max(0, cx - half)
    c1 = min(W, c0 + PATCH_SIZE)
    return arr[r0:r1, c0:c1]


def main() -> None:
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)

    print("Loading H&E crop...")
    he = tifffile.imread(str(HE_CROP))
    if he.ndim == 3 and he.shape[-1] == 4:
        he = he[..., :3]
    H, W = he.shape[:2]

    print("Loading masks...")
    masks = {name: np.load(str(info["mask"])) for name, info in METHODS.items()}

    n_boxes   = len(BOXES)
    n_methods = len(METHODS)
    fig, axes = plt.subplots(n_boxes, n_methods,
                             figsize=(4 * n_methods, 4 * n_boxes))

    for row_i, (cx, cy, box_label) in enumerate(BOXES):
        print(f"  {box_label}")
        he_patch = get_patch(he, cy, cx, H, W)

        for col_j, (name, info) in enumerate(METHODS.items()):
            mask_patch = get_patch(masks[name], cy, cx, H, W)
            img = overlay_patch(he_patch, mask_patch, info["color"])
            ax = axes[row_i, col_j]
            ax.imshow(img)
            ax.axis("off")

            if row_i == 0:
                ax.set_title(name, fontsize=12, fontweight="bold", pad=6)
            if col_j == 0:
                ax.set_ylabel(box_label, fontsize=8, rotation=90,
                              labelpad=6, va="center")

    fig.suptitle("5-region × 5-method segmentation comparison",
                 fontsize=14, fontweight="bold", y=1.002)
    plt.tight_layout()
    fig.savefig(str(OUT_FIG), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {OUT_FIG}")


if __name__ == "__main__":
    main()
