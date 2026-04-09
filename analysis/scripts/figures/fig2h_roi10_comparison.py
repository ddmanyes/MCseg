"""
34_fig2e_roi10_platform_comparison.py
======================================
Compare MCseg v2 vs Xenium GT in ROI10 (Normal Alveolar):
  1. Cell type composition (MCseg hybrid vs Xenium winner-take-all)
  2. Pixel-level area overlap (IoU, coverage)
  3. KEY METRIC: Among Xenium GT AT1 cells, what fraction overlaps MCseg regions?
     → Distinguishes segmentation miss vs data sparsity

Outputs:
  - manuscript/data/roi10_platform_comparison.csv
  - manuscript/figures/fig2/fig2e.png   (replaces old dotplot)

Run:
  cd /Volumes/SSD/plan_a
  uv run python manuscript/scripts/34_fig2e_roi10_platform_comparison.py
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from skimage.segmentation import find_boundaries
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT  = Path("/Volumes/SSD/plan_a/xenium_he_seg")
XENIUM_DIR    = Path("/Volumes/SSD/plan_a/tissue sample/LUAD/xenium")
MASK_DIR      = PROJECT_ROOT / "results" / "masks"
H5AD_DIR      = PROJECT_ROOT / "results" / "visiumhd" / "visiumhd_cells"
OUT_FIG       = Path("/Volumes/SSD/plan_a/manuscript/figures/fig2/fig2f.png")
OUT_CSV       = Path("/Volumes/SSD/plan_a/manuscript/data/roi10_platform_comparison.csv")

sys.path.insert(0, str(PROJECT_ROOT))
from backend.src.utils.alignment import load_alignment_matrix, he_pixel_to_xe_um, xe_um_to_he_pixel

# ROI10 definition (H&E fullres coords, consistent with 33_fig2_spatial_maps.py)
ROI10_X0, ROI10_Y0 = 7562, 19440
ROI10_W,  ROI10_H  = 3194, 1587
MARGIN_UM = 15.0

# ── Gene sets (must match 33_fig2_spatial_maps.py) ────────────────────────────
MCSEG_GENE_SETS = {
    "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"],
    "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
    "AT1 Pneumocyte":      ["AGER", "RTKN2"],          # threshold-based only
}

# Xenium uses different panel (no SFTPC/SFTPB/SFTPA1; uses LAMP3 for AT2)
XENIUM_GENE_SETS = {
    "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "LAMP3"],
    "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
    "AT1 Pneumocyte":      ["AGER", "RTKN2"],
}

PALETTE = {
    "AT2 Pneumocyte":      "#2166AC",
    "Alveolar Macrophage": "#D4841A",
    "Endothelial":         "#1A7340",
    "AT1 Pneumocyte":      "#B2182B",
    "Unassigned":          "#8C8C8C",
}

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "pdf.fonttype": 42,
    "axes.linewidth": 0.8,
})
MM = 1 / 25.4


# ── Helpers ───────────────────────────────────────────────────────────────────

def normalize_log1p(X_sparse) -> np.ndarray:
    X = sp.csr_matrix(X_sparse)
    totals = np.array(X.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    X = X.multiply(1e4 / totals[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X.toarray()


def score_hybrid(X_log: np.ndarray, var_names: np.ndarray,
                  gene_sets: dict[str, list[str]]) -> np.ndarray:
    """
    AT1: threshold (AGER > 0 OR RTKN2 > 0) applied first.
    Others: winner-take-all on remaining cells.
    """
    n_cells = X_log.shape[0]
    gene_idx = {g: i for i, g in enumerate(var_names)}
    cell_type = np.full(n_cells, "Unassigned", dtype=object)

    # Step 1: AT1 threshold
    at1_mask = np.zeros(n_cells, dtype=bool)
    for g in gene_sets.get("AT1 Pneumocyte", []):
        if g in gene_idx:
            at1_mask |= (X_log[:, gene_idx[g]] > 0)
    cell_type[at1_mask] = "AT1 Pneumocyte"

    # Step 2: winner-take-all for remaining
    wta_types = [t for t in gene_sets if t != "AT1 Pneumocyte"]
    remaining  = ~at1_mask
    scores     = np.zeros((n_cells, len(wta_types)), dtype=np.float32)
    for ci, tname in enumerate(wta_types):
        for g in gene_sets[tname]:
            if g in gene_idx:
                scores[:, ci] = np.maximum(scores[:, ci], X_log[:, gene_idx[g]])
    winner_idx = scores.argmax(axis=1)
    max_score  = scores.max(axis=1)
    wta_result = np.where(max_score > 0, np.array(wta_types)[winner_idx], "Unassigned")
    cell_type[remaining] = wta_result[remaining]
    return cell_type


# ── 1. Load MCseg data ─────────────────────────────────────────────────────────

def load_mcseg() -> tuple[np.ndarray, dict[str, int], np.ndarray]:
    """Returns (binary_mask, type_counts, cell_ids_per_pixel info)."""
    print("[MCseg] Loading mask & h5ad...")
    mask_full = np.load(str(MASK_DIR / "vhd_roi10_v12.npy"))
    adata     = ad.read_h5ad(str(H5AD_DIR / "roi10_v12.h5ad"))
    gene_names = np.array(adata.var_names)
    X_log      = normalize_log1p(adata.X)

    cell_type  = score_hybrid(X_log, gene_names, MCSEG_GENE_SETS)
    cell_ids   = adata.obs["cell_id"].values.astype(int)

    # build instance mask in local ROI10 coords
    # The stored mask is already in ROI10 local pixels
    mask_local = mask_full  # shape (ROI10_H, ROI10_W)

    # per-cell type binary mask
    id_max = int(mask_local.max())
    # LUT: cell_id → cell_type
    lut = np.full(id_max + 1, "Unassigned", dtype=object)
    for i, cid in enumerate(cell_ids):
        if 0 < cid <= id_max:
            lut[cid] = cell_type[i]

    binary = (mask_local > 0)
    all_types = list(MCSEG_GENE_SETS.keys()) + ["Unassigned"]
    type_counts = {t: int((cell_type == t).sum()) for t in all_types}
    total = sum(type_counts.values())
    print(f"  MCseg total cells: {total}")
    print(f"  Breakdown: { {k:v for k,v in type_counts.items() if v>0} }")
    return binary, type_counts, mask_local, lut


# ── 2. Generate Xenium GT mask & type labels ───────────────────────────────────

def load_xenium_gt() -> tuple[np.ndarray, dict[str, int], np.ndarray, np.ndarray]:
    """
    Rasterize Xenium cell_boundaries for ROI10 cells.
    Returns (binary_mask, type_counts, instance_mask, cell_type_lut).
    """
    print("[Xenium] Loading alignment matrix & filtering ROI10 cells...")
    M = load_alignment_matrix(XENIUM_DIR)

    corners_px = np.array([
        [ROI10_X0,           ROI10_Y0],
        [ROI10_X0 + ROI10_W, ROI10_Y0],
        [ROI10_X0 + ROI10_W, ROI10_Y0 + ROI10_H],
        [ROI10_X0,           ROI10_Y0 + ROI10_H],
    ], dtype=float)
    corners_um = he_pixel_to_xe_um(corners_px, M)
    xe_xmin = corners_um[:, 0].min() - MARGIN_UM
    xe_xmax = corners_um[:, 0].max() + MARGIN_UM
    xe_ymin = corners_um[:, 1].min() - MARGIN_UM
    xe_ymax = corners_um[:, 1].max() + MARGIN_UM

    cells_df = pd.read_parquet(XENIUM_DIR / "cells.parquet")
    in_roi   = ((cells_df.x_centroid >= xe_xmin) & (cells_df.x_centroid <= xe_xmax) &
                (cells_df.y_centroid >= xe_ymin) & (cells_df.y_centroid <= xe_ymax))
    roi_cells = cells_df[in_roi]
    sel_ids   = roi_cells["cell_id"].values
    print(f"  Xenium cells in ROI10: {len(sel_ids)}")

    # Rasterize boundaries
    print("[Xenium] Rasterizing cell_boundaries.parquet ...")
    import scanpy as sc
    cb     = pd.read_parquet(XENIUM_DIR / "cell_boundaries.parquet")
    cb     = cb[cb["cell_id"].isin(sel_ids)]

    mask_xe    = np.zeros((ROI10_H, ROI10_W), dtype=np.int32)
    cell_id_list = []
    for idx, (cell_id, grp) in enumerate(cb.groupby("cell_id"), start=1):
        verts_um    = grp[["vertex_x", "vertex_y"]].values.astype(np.float64)
        verts_he    = xe_um_to_he_pixel(verts_um, M)
        verts_local = verts_he - np.array([ROI10_X0, ROI10_Y0])
        pts = verts_local.astype(np.int32).reshape(-1, 1, 2)
        cv2.fillPoly(mask_xe, [pts], color=idx)
        cell_id_list.append(str(cell_id))
    print(f"  Rasterized {len(cell_id_list)} cells")

    # Xenium expression → cell type
    print("[Xenium] Loading cell_feature_matrix.h5 for typing...")
    adata_full = sc.read_10x_h5(str(XENIUM_DIR / "cell_feature_matrix.h5"))
    adata_full.obs_names = [s.decode() if isinstance(s, bytes) else s
                            for s in adata_full.obs_names]
    keep_set  = set(str(x) for x in sel_ids)
    keep_mask = [obs in keep_set for obs in adata_full.obs_names]
    adata_xe  = adata_full[keep_mask].copy()
    X_log_xe  = normalize_log1p(adata_xe.X)

    xe_cell_type = score_hybrid(X_log_xe, np.array(adata_xe.var_names), XENIUM_GENE_SETS)
    obs_to_type  = {obs: lbl for obs, lbl in zip(adata_xe.obs_names, xe_cell_type)}

    # Build LUT (rasterized label index → cell type)
    id_max   = int(mask_xe.max())
    lut_type = np.full(id_max + 1, "Unassigned", dtype=object)
    for lbl_idx, cell_id_str in enumerate(cell_id_list, start=1):
        if lbl_idx <= id_max:
            lut_type[lbl_idx] = obs_to_type.get(cell_id_str, "Unassigned")

    binary_xe = (mask_xe > 0)
    all_types = list(XENIUM_GENE_SETS.keys()) + ["Unassigned"]
    type_counts = {t: int((xe_cell_type == t).sum()) for t in all_types}
    total = sum(type_counts.values())
    print(f"  Xenium total cells: {total}")
    print(f"  Breakdown: { {k:v for k,v in type_counts.items() if v>0} }")
    return binary_xe, type_counts, mask_xe, lut_type


# ── 3. Compute overlap metrics ─────────────────────────────────────────────────

def compute_overlap_metrics(
        mcg_binary: np.ndarray,
        xe_binary:  np.ndarray,
        xe_mask:    np.ndarray,
        xe_lut:     np.ndarray,
        mcg_mask:   np.ndarray,
) -> dict:
    """
    Key metrics:
      - overall IoU between MCseg and Xenium GT binary masks
      - AT1-specific: % of Xenium AT1 pixels covered by MCseg
      - MCseg pixels that overlap Xenium AT1 (numerator for recall)
    """
    intersection = (mcg_binary & xe_binary).sum()
    union        = (mcg_binary | xe_binary).sum()
    iou          = intersection / union if union > 0 else 0.0

    xe_area      = xe_binary.sum()
    mcg_area     = mcg_binary.sum()
    coverage_xe  = intersection / xe_area  if xe_area  > 0 else 0.0  # how much of Xenium MCseg covers
    coverage_mcg = intersection / mcg_area if mcg_area > 0 else 0.0  # how much of MCseg overlaps Xenium

    # AT1-specific: pixels belonging to Xenium AT1 cells
    xe_at1_binary = np.isin(xe_lut[xe_mask], ["AT1 Pneumocyte"])
    at1_pixels    = xe_at1_binary[xe_mask > 0]  # only within Xenium cells
    # Rebuild per-pixel AT1 mask on full image
    at1_img = (xe_lut[xe_mask] == "AT1 Pneumocyte")  # bool H×W

    at1_total_pixels     = at1_img.sum()
    at1_covered_by_mcg   = (at1_img & mcg_binary).sum()
    at1_recall_pixel     = at1_covered_by_mcg / at1_total_pixels if at1_total_pixels > 0 else 0.0

    print(f"\n── Overlap Metrics ──────────────────────────────")
    print(f"  Overall IoU (MCseg ∩ Xenium / MCseg ∪ Xenium): {iou:.3f}")
    print(f"  MCseg coverage of Xenium area:                  {coverage_xe:.1%}")
    print(f"  Xenium area inside MCseg:                       {coverage_mcg:.1%}")
    print(f"  AT1 pixel recall (Xenium AT1 covered by MCseg): {at1_recall_pixel:.1%}")
    print(f"  Xenium AT1 pixels: {at1_total_pixels:,}, covered: {at1_covered_by_mcg:,}")

    return dict(
        iou=float(iou),
        coverage_xe=float(coverage_xe),
        coverage_mcg=float(coverage_mcg),
        at1_pixel_recall=float(at1_recall_pixel),
        at1_total_pixels=int(at1_total_pixels),
        at1_covered_pixels=int(at1_covered_by_mcg),
    )


# ── 4. Save CSV ────────────────────────────────────────────────────────────────

def save_csv(mcg_counts: dict, xe_counts: dict, overlap: dict) -> None:
    all_types = list(MCSEG_GENE_SETS.keys()) + ["Unassigned"]
    mcg_total = sum(mcg_counts.values())
    xe_total  = sum(xe_counts.values())

    rows = []
    for t in all_types:
        rows.append({
            "cell_type":       t,
            "mcseg_n":         mcg_counts.get(t, 0),
            "mcseg_pct":       round(mcg_counts.get(t, 0) / mcg_total * 100, 2),
            "xenium_n":        xe_counts.get(t, 0),
            "xenium_pct":      round(xe_counts.get(t, 0) / xe_total * 100, 2),
        })

    df_types = pd.DataFrame(rows)

    df_overlap = pd.DataFrame([{
        "metric":                   "overall_iou",
        "value":                    round(overlap["iou"], 4),
        "description":              "Pixel IoU between MCseg and Xenium GT binary masks",
    }, {
        "metric":                   "xenium_area_covered_by_mcseg",
        "value":                    round(overlap["coverage_xe"], 4),
        "description":              "Fraction of Xenium-covered pixels also covered by MCseg",
    }, {
        "metric":                   "mcseg_area_within_xenium",
        "value":                    round(overlap["coverage_mcg"], 4),
        "description":              "Fraction of MCseg pixels falling inside a Xenium cell",
    }, {
        "metric":                   "at1_pixel_recall",
        "value":                    round(overlap["at1_pixel_recall"], 4),
        "description":              "Fraction of Xenium GT AT1 pixels covered by MCseg",
    }])

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    # Write two sheets to same CSV with separator
    with open(OUT_CSV, "w") as f:
        f.write("# ROI10 Platform Comparison: MCseg v2 vs Xenium GT\n")
        f.write("# Section 1: Cell type composition\n")
    df_types.to_csv(OUT_CSV, mode="a", index=False)
    with open(OUT_CSV, "a") as f:
        f.write("\n# Section 2: Pixel overlap metrics\n")
    df_overlap.to_csv(OUT_CSV, mode="a", index=False)
    print(f"\n✅ CSV saved: {OUT_CSV}")


# ── 5. Generate figure ─────────────────────────────────────────────────────────

def make_figure(mcg_counts: dict, xe_counts: dict, overlap: dict) -> None:
    all_types = ["AT1 Pneumocyte", "AT2 Pneumocyte", "Alveolar Macrophage",
                 "Endothelial", "Unassigned"]
    mcg_total = sum(mcg_counts.values())
    xe_total  = sum(xe_counts.values())

    mcg_pct = [mcg_counts.get(t, 0) / mcg_total * 100 for t in all_types]
    xe_pct  = [xe_counts.get(t,  0) / xe_total  * 100 for t in all_types]
    colors  = [PALETTE[t] for t in all_types]

    fig_w = 90 * MM
    fig_h = 80 * MM
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h),
                              gridspec_kw=dict(width_ratios=[3, 1],
                                               left=0.13, right=0.97,
                                               top=0.88, bottom=0.10,
                                               wspace=0.55))

    # ── Panel A: grouped bar chart ──
    ax = axes[0]
    n = len(all_types)
    x = np.arange(n)
    w = 0.35

    bars1 = ax.bar(x - w/2, mcg_pct, w, color=colors, edgecolor="white",
                   linewidth=0.5, zorder=3, label="MCseg v2")
    bars2 = ax.bar(x + w/2, xe_pct,  w, color=colors, edgecolor="black",
                   linewidth=0.5, linestyle="--", zorder=3, label="Xenium GT",
                   alpha=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(
        ["AT1", "AT2", "Macro-\nphage", "Endo-\nthelial", "Un-\nassigned"],
        fontsize=6.5)
    ax.set_ylabel("% of cells", fontsize=7)
    ax.set_ylim(0, max(max(mcg_pct), max(xe_pct)) * 1.25)
    ax.grid(axis="y", alpha=0.2, lw=0.5, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # AT1 annotation: show absolute counts
    at1_idx = 0
    ax.annotate(f"n={mcg_counts.get('AT1 Pneumocyte', 0)}",
                xy=(at1_idx - w/2, mcg_pct[at1_idx]),
                xytext=(0, 3), textcoords="offset points",
                ha="center", fontsize=5.5, color="#B2182B")
    ax.annotate(f"n={xe_counts.get('AT1 Pneumocyte', 0)}",
                xy=(at1_idx + w/2, xe_pct[at1_idx]),
                xytext=(0, 3), textcoords="offset points",
                ha="center", fontsize=5.5, color="#B2182B")

    ax.set_title("Cell-type composition\nROI10 (Normal Alveolar)",
                 fontsize=7, fontweight="bold", pad=4)

    # legend
    from matplotlib.patches import Patch
    leg = [Patch(facecolor="#aaa", edgecolor="white", label=f"MCseg v2 (n={mcg_total:,})"),
           Patch(facecolor="#aaa", edgecolor="black", linewidth=0.8,
                 alpha=0.6, label=f"Xenium GT (n={xe_total:,})")]
    ax.legend(handles=leg, fontsize=5.5, loc="upper right",
              framealpha=0.85, edgecolor="#ccc", handlelength=1.2)

    # panel label
    ax.text(-0.18, 1.06, "e", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")

    # ── Panel B: overlap stats ──
    ax2 = axes[1]
    ax2.axis("off")

    metrics = [
        ("Overall IoU",          f"{overlap['iou']:.2f}"),
        ("MCseg covers\nXenium area",  f"{overlap['coverage_xe']:.0%}"),
        ("AT1 pixel\nrecall",    f"{overlap['at1_pixel_recall']:.0%}"),
    ]

    y_pos  = 0.88
    ax2.text(0.05, 0.97, "Area overlap\nmetrics",
             transform=ax2.transAxes,
             fontsize=7, fontweight="bold", va="top")
    for label, val in metrics:
        ax2.text(0.05, y_pos, label,
                 transform=ax2.transAxes,
                 fontsize=6, va="top", color="#555")
        ax2.text(0.95, y_pos - 0.04, val,
                 transform=ax2.transAxes,
                 fontsize=8, fontweight="bold", va="top", ha="right",
                 color="#222")
        y_pos -= 0.25

    # AT1 recall highlight box
    recall_pct = overlap["at1_pixel_recall"] * 100
    ax2.add_patch(plt.Rectangle((0.0, 0.03), 1.0, 0.18,
                                 transform=ax2.transAxes,
                                 facecolor="#FCE8E8", edgecolor="#B2182B",
                                 linewidth=0.8, zorder=0))
    ax2.text(0.5, 0.17, f"AT1 area\ncaptured by MCseg",
             transform=ax2.transAxes, fontsize=5.5, ha="center", va="top", color="#B2182B")
    ax2.text(0.5, 0.09, f"{recall_pct:.0f}%",
             transform=ax2.transAxes, fontsize=10, fontweight="bold",
             ha="center", va="top", color="#B2182B")

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✅ Figure saved: {OUT_FIG}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("ROI10 Platform Comparison: MCseg v2 vs Xenium GT")
    print("=" * 60)

    mcg_binary, mcg_counts, mcg_mask, mcg_lut = load_mcseg()
    xe_binary,  xe_counts,  xe_mask,  xe_lut  = load_xenium_gt()

    overlap = compute_overlap_metrics(
        mcg_binary, xe_binary, xe_mask, xe_lut, mcg_mask)

    save_csv(mcg_counts, xe_counts, overlap)
    make_figure(mcg_counts, xe_counts, overlap)

    print("\n── Summary ──────────────────────────────────────────")
    print(f"  MCseg v2  : {sum(mcg_counts.values()):,} cells  "
          f"| AT1 = {mcg_counts.get('AT1 Pneumocyte', 0)} ({mcg_counts.get('AT1 Pneumocyte', 0)/sum(mcg_counts.values())*100:.1f}%)")
    print(f"  Xenium GT : {sum(xe_counts.values()):,} cells  "
          f"| AT1 = {xe_counts.get('AT1 Pneumocyte', 0)} ({xe_counts.get('AT1 Pneumocyte', 0)/sum(xe_counts.values())*100:.1f}%)")
    print(f"  AT1 pixel recall (Xenium AT1 area covered by MCseg): "
          f"{overlap['at1_pixel_recall']:.1%}")
    print(f"  → If recall is HIGH but MCseg AT1 detection is low (~2%),")
    print(f"    this confirms DATA SPARSITY as the dominant mechanism.")


if __name__ == "__main__":
    main()
