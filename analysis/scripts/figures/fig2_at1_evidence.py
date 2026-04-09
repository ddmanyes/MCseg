"""
36_fig2f_at1_evidence.py
========================
Generate fig2f: AT1 detection evidence panel (replaces simple bar chart)

Layout:
  Left (60%): Spatial overlay — H&E + Xenium AT1 fill (red) + MCseg borders (blue)
  Right (40%): Waterfall chart — Xenium AT1 → MCseg geometric coverage → RNA detection

Key message:
  MCseg geometrically covers 71% of Xenium AT1 area,
  yet only 2.1% of MCseg cells express AGER/RTKN2.
  → Data sparsity, not segmentation miss.

Output: manuscript/figures/fig2/fig2f.png
Run:
  cd /Volumes/SSD/plan_a
  uv run python manuscript/scripts/36_fig2f_at1_evidence.py
"""

from __future__ import annotations
import sys, warnings
warnings.filterwarnings("ignore")
import os
os.environ["OMP_NUM_THREADS"] = "1"

import cv2
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import tifffile, zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from skimage.segmentation import find_boundaries
from pathlib import Path

PROJECT_ROOT = Path("/Volumes/SSD/plan_a/xenium_he_seg")
XENIUM_DIR   = Path("/Volumes/SSD/plan_a/tissue sample/LUAD/xenium")
BTF_PATH     = (Path("/Volumes/SSD/plan_a/tissue sample/LUAD/visium") /
                "Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2_tissue_image.btf")
MASK_DIR     = PROJECT_ROOT / "results" / "masks"
H5AD_DIR     = PROJECT_ROOT / "results" / "visiumhd" / "visiumhd_cells"
OUT_FIG2E    = Path("/Volumes/SSD/plan_a/manuscript/figures/fig2/fig2e.png")  # spatial overlay
OUT_FIG2F    = Path("/Volumes/SSD/plan_a/manuscript/figures/fig2/fig2f.png")  # rate chart

sys.path.insert(0, str(PROJECT_ROOT))
from backend.src.utils.alignment import load_alignment_matrix, he_pixel_to_xe_um, xe_um_to_he_pixel

ROI10_X0, ROI10_Y0 = 7562, 19440
ROI10_W,  ROI10_H  = 3194, 1587
MARGIN_UM = 15.0

XENIUM_GENE_SETS = {
    "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "LAMP3"],
    "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
    "AT1 Pneumocyte":      ["AGER", "RTKN2"],
}

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "pdf.fonttype": 42,
})
MM = 1 / 25.4


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_he_crop() -> np.ndarray:
    with tifffile.TiffFile(str(BTF_PATH)) as tif:
        store = tif.aszarr()
        z = zarr.open(store, mode="r")
        arr = z if not isinstance(z, zarr.Group) else z[0]
        crop = np.array(arr[ROI10_Y0:ROI10_Y0 + ROI10_H,
                            ROI10_X0:ROI10_X0 + ROI10_W])
    if crop.ndim == 3 and crop.shape[2] == 4:
        crop = crop[:, :, :3]
    return crop


def normalize_log1p(X_sparse) -> np.ndarray:
    X = sp.csr_matrix(X_sparse)
    totals = np.array(X.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    X = X.multiply(1e4 / totals[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X.toarray()


def score_hybrid(X_log, var_names, gene_sets):
    n = X_log.shape[0]
    gidx = {g: i for i, g in enumerate(var_names)}
    ct = np.full(n, "Unassigned", dtype=object)
    at1 = np.zeros(n, dtype=bool)
    for g in gene_sets.get("AT1 Pneumocyte", []):
        if g in gidx:
            at1 |= (X_log[:, gidx[g]] > 0)
    ct[at1] = "AT1 Pneumocyte"
    wta = [t for t in gene_sets if t != "AT1 Pneumocyte"]
    sc  = np.zeros((n, len(wta)), dtype=np.float32)
    for ci, t in enumerate(wta):
        for g in gene_sets[t]:
            if g in gidx:
                sc[:, ci] = np.maximum(sc[:, ci], X_log[:, gidx[g]])
    rem = ~at1
    ct[rem] = np.where(sc[rem].max(1) > 0,
                       np.array(wta)[sc[rem].argmax(1)], "Unassigned")
    return ct


def build_xenium_at1_mask() -> tuple[np.ndarray, int]:
    """Returns (AT1_binary_mask H×W, n_xenium_at1_cells)."""
    import scanpy as sc
    print("  [Xenium] alignment + bbox...")
    M = load_alignment_matrix(XENIUM_DIR)
    corners = np.array([[ROI10_X0, ROI10_Y0],
                        [ROI10_X0 + ROI10_W, ROI10_Y0],
                        [ROI10_X0 + ROI10_W, ROI10_Y0 + ROI10_H],
                        [ROI10_X0, ROI10_Y0 + ROI10_H]], dtype=float)
    um = he_pixel_to_xe_um(corners, M)
    xmin, xmax = um[:, 0].min() - MARGIN_UM, um[:, 0].max() + MARGIN_UM
    ymin, ymax = um[:, 1].min() - MARGIN_UM, um[:, 1].max() + MARGIN_UM

    cells_df = pd.read_parquet(XENIUM_DIR / "cells.parquet")
    in_roi   = ((cells_df.x_centroid >= xmin) & (cells_df.x_centroid <= xmax) &
                (cells_df.y_centroid >= ymin) & (cells_df.y_centroid <= ymax))
    sel_ids  = cells_df[in_roi]["cell_id"].values

    print("  [Xenium] rasterizing boundaries...")
    cb = pd.read_parquet(XENIUM_DIR / "cell_boundaries.parquet")
    cb = cb[cb["cell_id"].isin(sel_ids)]

    mask_xe = np.zeros((ROI10_H, ROI10_W), dtype=np.int32)
    cell_id_list = []
    for idx, (cid, grp) in enumerate(cb.groupby("cell_id"), start=1):
        verts = grp[["vertex_x", "vertex_y"]].values.astype(np.float64)
        pts   = (xe_um_to_he_pixel(verts, M) - [ROI10_X0, ROI10_Y0]
                 ).astype(np.int32).reshape(-1, 1, 2)
        cv2.fillPoly(mask_xe, [pts], color=idx)
        cell_id_list.append(str(cid))

    print("  [Xenium] typing cells...")
    adata_full = sc.read_10x_h5(str(XENIUM_DIR / "cell_feature_matrix.h5"))
    adata_full.obs_names = [s.decode() if isinstance(s, bytes) else s
                            for s in adata_full.obs_names]
    keep = [obs in set(str(x) for x in sel_ids) for obs in adata_full.obs_names]
    adata_xe = adata_full[keep].copy()
    X_log    = normalize_log1p(adata_xe.X)
    ct       = score_hybrid(X_log, np.array(adata_xe.var_names), XENIUM_GENE_SETS)
    obs2type = {obs: lbl for obs, lbl in zip(adata_xe.obs_names, ct)}

    id_max   = int(mask_xe.max())
    lut      = np.full(id_max + 1, "Unassigned", dtype=object)
    for i, cid_str in enumerate(cell_id_list, start=1):
        if i <= id_max:
            lut[i] = obs2type.get(cid_str, "Unassigned")

    at1_binary = (lut[mask_xe] == "AT1 Pneumocyte")
    n_at1 = int((ct == "AT1 Pneumocyte").sum())
    print(f"  Xenium AT1: n={n_at1}, pixels={at1_binary.sum():,}")
    return at1_binary, n_at1, mask_xe, lut


# ── Panel A: spatial overlay ───────────────────────────────────────────────────

def make_spatial_overlay(ax, at1_binary, mcg_mask):
    """Clean mask-only view: Xenium AT1 fill + MCseg boundaries on white background."""
    H, W = ROI10_H, ROI10_W

    # White background
    canvas = np.ones((H, W, 3), dtype=np.float32)

    # Xenium AT1 fill: cyan
    at1_rgb = np.array([0.0, 0.76, 0.80])   # #00C2CC
    canvas[at1_binary] = at1_rgb

    # MCseg boundaries: thickened (3px) dark outlines
    from skimage.morphology import dilation, disk
    bounds     = find_boundaries(mcg_mask, mode="outer").astype(np.uint8)
    bounds_fat = dilation(bounds, disk(2)).astype(bool)   # ~3px wide
    canvas[bounds_fat] = [0.15, 0.15, 0.15]   # near-black

    comp = (canvas * 255).clip(0, 255).astype(np.uint8)
    ax.imshow(comp, origin="upper", interpolation="nearest")
    ax.axis("off")

    # Legend
    handles = [
        mpatches.Patch(facecolor="#00C2CC", edgecolor="none",
                       label="Xenium GT  AT1 cells"),
        mpatches.Patch(facecolor="white", edgecolor="#262626", linewidth=1.5,
                       label="MCseg v2  cell boundaries"),
    ]
    ax.legend(handles=handles, fontsize=6.5, loc="lower right",
              framealpha=0.90, edgecolor="#bbb",
              handlelength=1.2, borderpad=0.5, labelspacing=0.3)


# ── Panel B: waterfall chart ───────────────────────────────────────────────────

def make_waterfall(ax, n_xe_total, n_xe_at1, n_mcg_in_at1,
                   n_mcg_total, n_mcg_at1):
    """Three-step attrition waterfall."""

    steps = [
        ("Xenium GT\nAT1 cells",
         n_xe_at1,
         f"n = {n_xe_at1}",
         "#B2182B"),
        ("MCseg cells\n≥50% overlap\nwith AT1 area",
         n_mcg_in_at1,
         f"n = {n_mcg_in_at1}",
         "#D6604D"),
        ("MCseg AT1\ndetected\n(AGER⁺ or RTKN2⁺)",
         n_mcg_at1,
         f"n = {n_mcg_at1}",
         "#F4A582"),
    ]

    bar_w = 0.42
    x_pos = [0, 1, 2]
    max_n = max(s[1] for s in steps)

    for x, (label, val, annot, color) in zip(x_pos, steps):
        ax.bar(x, val, width=bar_w, color=color,
               edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(x, val + max_n * 0.025, annot,
                ha="center", va="bottom", fontsize=7,
                color="#333", linespacing=1.4)

    # Connector bands + drop % annotation
    for i in range(len(steps) - 1):
        y_lo = min(steps[i][1], steps[i + 1][1])
        y_hi = max(steps[i][1], steps[i + 1][1])
        x_l  = x_pos[i] + bar_w / 2
        x_r  = x_pos[i + 1] - bar_w / 2
        ax.fill_betweenx([y_lo, y_hi], x_l, x_r,
                          color="#EEEEEE", alpha=0.7, zorder=2)
        ax.plot([x_l, x_r], [steps[i][1], steps[i][1]],
                color="#AAAAAA", lw=0.6, zorder=3)
        ax.plot([x_l, x_r], [steps[i + 1][1], steps[i + 1][1]],
                color="#AAAAAA", lw=0.6, zorder=3)
        drop_pct = (steps[i][1] - steps[i + 1][1]) / steps[i][1] * 100
        ax.text((x_l + x_r) / 2, (steps[i][1] + steps[i + 1][1]) / 2,
                f"▼ {drop_pct:.0f}%",
                ha="center", va="center", fontsize=7, color="#555",
                fontweight="bold", zorder=4)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([s[0] for s in steps], fontsize=7)
    ax.set_ylabel("Number of cells", fontsize=7)
    ax.set_ylim(0, max_n * 1.30)
    ax.set_xlim(-0.55, 2.55)
    ax.grid(axis="y", alpha=0.18, lw=0.5, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("[1] Loading MCseg mask & h5ad...")
    mcg_mask = np.load(str(MASK_DIR / "vhd_roi10_v12.npy"))
    adata    = ad.read_h5ad(str(H5AD_DIR / "roi10_v12.h5ad"))
    X_log    = normalize_log1p(adata.X)
    var_names = np.array(adata.var_names)

    mcg_gene_sets = {
        "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"],
        "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
        "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
        "AT1 Pneumocyte":      ["AGER", "RTKN2"],
    }
    mcg_ct      = score_hybrid(X_log, var_names, mcg_gene_sets)
    n_mcg_total = adata.n_obs
    n_mcg_at1   = int((mcg_ct == "AT1 Pneumocyte").sum())
    print(f"  MCseg: {n_mcg_total} cells, AT1={n_mcg_at1}")

    print("[2] Building Xenium AT1 mask...")
    at1_binary, n_xe_at1, xe_mask, xe_lut = build_xenium_at1_mask()
    n_xe_total = int((xe_lut[xe_mask] != "Unassigned").sum())
    # use actual cell count from Xenium obs
    # n_xe_total = number of Xenium cells (obs)
    # get from xe_lut unique non-background
    n_xe_total_cells = int(xe_mask.max())   # number of rasterized cells

    print("[3] Computing MCseg cells overlapping AT1 territory (>=50% area in AT1)...")
    mcg_labels_in_at1 = mcg_mask[at1_binary]
    nonzero_labels, overlap_px = np.unique(
        mcg_labels_in_at1[mcg_labels_in_at1 > 0], return_counts=True)
    # keep only cells where >=50% of their area falls inside AT1
    cell_total_px = np.array([(mcg_mask == lbl).sum() for lbl in nonzero_labels])
    overlap_frac  = overlap_px / cell_total_px
    mcg_in_at1_50 = nonzero_labels[overlap_frac >= 0.50]
    n_mcg_in_at1  = int(len(mcg_in_at1_50))
    print(f"  MCseg cells with >=50% area in AT1: {n_mcg_in_at1} "
          f"(any overlap: {len(nonzero_labels)})")

    # AT1 pixel recall
    at1_covered = (at1_binary & (mcg_mask > 0)).sum()
    at1_recall  = at1_covered / at1_binary.sum() if at1_binary.sum() > 0 else 0.0
    print(f"  AT1 pixel recall: {at1_recall:.1%}")

    print("[4] Rendering fig2e — spatial overlay (standalone)...")
    asp    = ROI10_H / ROI10_W
    fig_w  = 183 * MM
    fig_e, ax_e = plt.subplots(1, 1, figsize=(fig_w, fig_w * asp),
                                gridspec_kw=dict(left=0.01, right=0.99,
                                                 top=0.99, bottom=0.01))
    make_spatial_overlay(ax_e, at1_binary, mcg_mask)
    ax_e.text(-0.005, 1.01, "e", transform=ax_e.transAxes,
              fontsize=12, fontweight="bold", va="bottom")
    OUT_FIG2E.parent.mkdir(parents=True, exist_ok=True)
    fig_e.savefig(OUT_FIG2E, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_e)
    print(f"  ✅ {OUT_FIG2E}")

    print("[5] Rendering fig2f — AT1 detection rate chart (standalone)...")
    fig_f, ax_f = plt.subplots(1, 1, figsize=(90 * MM, 90 * MM),
                                gridspec_kw=dict(left=0.17, right=0.95,
                                                 top=0.92, bottom=0.20))
    # % rate bars: Xenium AT1, MCseg geometric, MCseg RNA
    rates  = [n_xe_at1 / n_xe_total_cells * 100,
              n_mcg_in_at1 / n_mcg_total * 100,
              n_mcg_at1 / n_mcg_total * 100]
    labels = ["Xenium GT\nAT1 cells",
              "MCseg cells\n≥50% overlap\nwith AT1 area",
              "MCseg AT1\ndetected\n(AGER⁺ or RTKN2⁺)"]
    ns     = [n_xe_at1, n_mcg_in_at1, n_mcg_at1]
    colors = ["#B2182B", "#D6604D", "#F4A582"]

    x_pos = [0, 1, 2]
    bar_w = 0.42
    max_r = max(rates)
    for x, r, lbl, n, c in zip(x_pos, rates, labels, ns, colors):
        ax_f.bar(x, r, width=bar_w, color=c, edgecolor="white", linewidth=0.8, zorder=3)
        ax_f.text(x, r + max_r * 0.025, f"{r:.1f}%\n(n={n})",
                  ha="center", va="bottom", fontsize=6.5, color="#333", linespacing=1.4)

    ax_f.set_xticks(x_pos)
    ax_f.set_xticklabels(labels, fontsize=6.5)
    ax_f.set_ylabel("% of cells", fontsize=7)
    ax_f.set_ylim(0, max_r * 1.45)
    ax_f.set_xlim(-0.55, 2.55)
    ax_f.grid(axis="y", alpha=0.18, lw=0.5, zorder=0)
    ax_f.spines["top"].set_visible(False)
    ax_f.spines["right"].set_visible(False)
    ax_f.text(-0.18, 1.04, "f", transform=ax_f.transAxes,
              fontsize=12, fontweight="bold", va="top")
    OUT_FIG2F.parent.mkdir(parents=True, exist_ok=True)
    fig_f.savefig(OUT_FIG2F, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_f)
    print(f"  ✅ {OUT_FIG2F}")
    print(f"\n── Key numbers ──────────────────────────────────")
    print(f"  Xenium AT1:               {n_xe_at1} cells")
    print(f"  MCseg cells in AT1 area:  {n_mcg_in_at1} cells ({n_mcg_in_at1/n_mcg_total:.1%})")
    print(f"  MCseg AT1 via RNA:        {n_mcg_at1} cells ({n_mcg_at1/n_mcg_total:.1%})")
    print(f"  AT1 pixel recall:         {at1_recall:.1%}")


if __name__ == "__main__":
    main()
