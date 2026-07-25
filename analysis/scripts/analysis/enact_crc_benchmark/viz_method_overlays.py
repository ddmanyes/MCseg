"""
viz_method_overlays.py
----------------------
Generate H&E overlay images for each segmentation method.

Outputs per method:
  1. Full-image overview (downsampled to ~2000px wide)
  2. Zoomed patch (1024×1024px from a representative region)

Usage:
    cd /Volumes/SSD/plan_a
    uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/viz_method_overlays.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries
from skimage.transform import resize

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A      = Path("/Volumes/SSD/plan_a")
F1_DIR      = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
RESULT_DIR  = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
HE_CROP     = F1_DIR / "he_crop.tif"
OUT_DIR     = RESULT_DIR / "overlays"

METHODS = {
    "StarDist": {"mask": RESULT_DIR / "stardist_mask.npy",    "color": "#FFD700"},
    "MCseg":    {"mask": F1_DIR     / "mcseg_mask_7pass.npy", "color": "#FF4444"},
    "SR":       {"mask": RESULT_DIR / "sr_mask.npy",          "color": "#44AAFF"},
    "NUC":      {"mask": RESULT_DIR / "cellpose_nuc_mask.npy", "color": "#44DD88"},
    "ProSeg":   {"mask": RESULT_DIR / "proseg_mask.npy",      "color": "#CC88FF"},
}

# Patch center (approximate tumor gland region in crop local coords)
PATCH_CY, PATCH_CX = 7000, 5000
PATCH_SIZE = 1024

OVERVIEW_WIDTH = 2000   # px, for full-image downsampled view


def hex_to_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def draw_contours_on_ax(ax, he_patch: np.ndarray, mask_patch: np.ndarray,
                        color: str, title: str) -> None:
    # find_boundaries correctly separates adjacent labeled cells;
    # find_contours(level=0.5) only finds bg→fg edges, missing inter-cell boundaries
    bnd = find_boundaries(mask_patch, mode="inner")
    overlay = he_patch.copy().astype(np.float32) / 255.0
    rgb = hex_to_rgb(color)
    overlay[bnd] = rgb
    ax.imshow(np.clip(overlay, 0, 1))
    n_cells = int(mask_patch.max())
    ax.set_title(f"{title}\n({n_cells:,} cells in patch)", fontsize=9)
    ax.axis("off")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading H&E crop...")
    he = tifffile.imread(str(HE_CROP))
    if he.ndim == 3 and he.shape[-1] == 4:
        he = he[..., :3]
    H, W = he.shape[:2]
    print(f"  H&E shape: {he.shape}")

    # Patch bounds
    r0 = max(0, PATCH_CY - PATCH_SIZE // 2)
    r1 = min(H, r0 + PATCH_SIZE)
    c0 = max(0, PATCH_CX - PATCH_SIZE // 2)
    c1 = min(W, c0 + PATCH_SIZE)
    he_patch = he[r0:r1, c0:c1]

    # Overview downscale factor
    scale = OVERVIEW_WIDTH / W
    ov_h = int(H * scale)
    ov_w = OVERVIEW_WIDTH
    he_overview = (resize(he, (ov_h, ov_w), anti_aliasing=True) * 255).astype(np.uint8)

    # ── 1. Combined patch figure (all 5 methods in one figure) ──────────────
    print("Generating combined patch figure...")
    fig, axes = plt.subplots(1, 5, figsize=(22, 5))
    for ax, (name, info) in zip(axes, METHODS.items()):
        mask = np.load(str(info["mask"]))
        mask_patch = mask[r0:r1, c0:c1]
        draw_contours_on_ax(ax, he_patch, mask_patch, info["color"], name)
        del mask

    fig.suptitle(
        f"Cell segmentation comparison — H&E patch ({PATCH_SIZE}×{PATCH_SIZE}px)\n"
        f"crop coords: row {r0}–{r1}, col {c0}–{c1}",
        fontsize=10, fontweight="bold",
    )
    plt.tight_layout()
    out = OUT_DIR / "patch_all_methods.png"
    fig.savefig(str(out), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")

    # ── 2. Individual full-image overviews (one per method) ─────────────────
    print("Generating full-image overviews...")
    for name, info in METHODS.items():
        print(f"  [{name}]")
        mask = np.load(str(info["mask"]))
        mask_ov = (resize(mask.astype(np.float32), (ov_h, ov_w),
                          anti_aliasing=False, order=0)).astype(np.int32)

        fig, ax = plt.subplots(figsize=(10, 14))
        draw_contours_on_ax(ax, he_overview, mask_ov, info["color"], name)
        ax.set_title(f"{name} — full crop overview\n({int(mask.max()):,} cells)", fontsize=11)
        plt.tight_layout()
        out = OUT_DIR / f"overview_{name.lower()}.png"
        fig.savefig(str(out), dpi=150, bbox_inches="tight")
        plt.close(fig)
        del mask, mask_ov
        print(f"    Saved: {out.name}")

    # ── 3. Combined overview figure (all 5 side-by-side) ────────────────────
    print("Generating combined overview figure...")
    fig, axes = plt.subplots(1, 5, figsize=(40, 12))
    for ax, (name, info) in zip(axes, METHODS.items()):
        mask = np.load(str(info["mask"]))
        mask_ov = (resize(mask.astype(np.float32), (ov_h, ov_w),
                          anti_aliasing=False, order=0)).astype(np.int32)
        draw_contours_on_ax(ax, he_overview, mask_ov, info["color"], name)
        del mask, mask_ov

    fig.suptitle("Cell segmentation comparison — full crop overview", fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = OUT_DIR / "overview_all_methods.png"
    fig.savefig(str(out), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")

    print(f"\nAll outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
