"""
enact_C2_crc_f1_wba.py
=======================
Phase 2: MCseg + Weight-by-Area (WBA) bin attribution — ENACT CRC F1 comparison

ENACT WBA principle:
  Each 2µm bin has a 7×7 px footprint (bin_spacing ≈ 7.29 px → BIN_HALF = 3).
  Fractionally assign each bin's transcripts to overlapping cells by the
  proportion of footprint pixels that fall in each cell mask label.
  This mirrors what ENACT does with StarDist cells.

Comparison matrix:
  C1  MCseg + lookup (winner-take-all, 1 px at bin centre)  → already computed
  C2  MCseg + WBA    (fractional 7×7 footprint)             → this script

StarDist + WBA baseline from ENACT paper: micro F1 = 0.708 / weighted F1 = 0.758

Reuses cached outputs from C1:
  results/enact_crc_f1/mcseg_mask.npy
  results/enact_crc_f1/he_crop.tif  (not re-read)
  results/enact_crc_f1/f1_summary.csv  (C1 F1 for 3-way comparison figure)

Output: submission_bioinformatics/results/enact_crc_f1_wba/
"""

from __future__ import annotations

import gc
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ─── Paths ────────────────────────────────────────────────────────────────────

PLAN_A    = Path("/Volumes/SSD/plan_a")
ENACT_CRC = PLAN_A / "tissue sample" / "ENACT_supporting_files" / "public_data" / "human_colorectal"
H5_PATH   = ENACT_CRC / "input_files" / "filtered_feature_bc_matrix.h5"
TP_PATH   = ENACT_CRC / "input_files" / "tissue_positions.parquet"
EVAL_CSV  = ENACT_CRC / "paper_results" / "chunks" / "weighted_by_area" / "sargent_results" / "eval" / "cell_annotation_eval.csv"

MSSEG_ROOT  = PLAN_A / "MSseg"
C1_DIR      = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
RESULTS_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1_wba"

# Coordinate constants (same as C1)
CROP_X0, CROP_X1 = 5154, 15242
CROP_Y0, CROP_Y1 = 4635, 18599
COL_OFFSET = 40598

# WBA footprint: bin spacing ≈ 7.29 px → half-window = 3 → 7×7 = 49 px
BIN_HALF = 3

# StarDist + WBA baseline (ENACT paper)
STARDIST_MICRO_F1    = 0.708
STARDIST_WEIGHTED_F1 = 0.758

LABEL_MAP: dict[str, str] = {
    "CMS1": "epithelial cells", "CMS2": "epithelial cells",
    "CMS3": "epithelial cells", "CMS4": "epithelial cells",
    "Goblet cells": "epithelial cells",
    "Mature Enterocytes type 1": "epithelial cells",
    "Mature Enterocytes type 2": "epithelial cells",
    "Stem-like/TA": "epithelial cells",
    "Intermediate": "epithelial cells",
    "Proliferating": "epithelial cells",
    "Myofibroblasts": "stromal cells",
    "Pericytes": "stromal cells",
    "Smooth muscle cells": "stromal cells",
    "Stromal 1": "stromal cells", "Stromal 2": "stromal cells", "Stromal 3": "stromal cells",
    "Lymphatic ECs": "stromal cells", "Proliferative ECs": "stromal cells",
    "Stalk-like ECs": "stromal cells", "Enteric glial cells": "stromal cells",
    "CD19+CD20+ B": "immune cells", "CD4+ T cells": "immune cells",
    "CD8+ T cells": "immune cells", "Regulatory T cells": "immune cells",
    "NK cells": "immune cells", "IgA+ Plasma": "immune cells",
    "IgG+ Plasma": "immune cells", "Mast cells": "immune cells",
    "Pro-inflammatory": "immune cells", "SPP1+": "immune cells",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


# ─── Step 1: Load GT ──────────────────────────────────────────────────────────

def load_gt() -> pd.DataFrame:
    log.info("Step 1: 載入 ENACT expert GT")
    gt = pd.read_csv(EVAL_CSV)[["cell_x", "cell_y", "gt_label"]].copy()
    gt["cell_x"] = gt["cell_x"].astype(np.float32)
    gt["cell_y"] = gt["cell_y"].astype(np.float32)
    log.info(f"  GT cells: {len(gt):,}")
    return gt


# ─── Step 2: Load cached MCseg mask ───────────────────────────────────────────

def load_mask() -> np.ndarray:
    mask_path = C1_DIR / "mcseg_mask.npy"
    if not mask_path.exists():
        raise FileNotFoundError(
            f"MCseg mask not found: {mask_path}\n"
            "Run enact_C1_crc_f1_comparison.py first."
        )
    log.info(f"Step 2: 載入 MCseg mask from C1: {mask_path.name}")
    mask = np.load(str(mask_path))
    log.info(f"  mask shape: {mask.shape}  cells: {int(mask.max()):,}")
    return mask


# ─── Step 3: WBA attribution ──────────────────────────────────────────────────

def run_wba_attribution(mask: np.ndarray, wba_csv: Path) -> pd.DataFrame:
    """
    Weight-by-Area bin attribution.

    For each Visium HD 2µm bin in the crop:
      1. Extract a (2*BIN_HALF+1)² window centred at the bin's pixel.
      2. Count how many pixels in that window belong to each cell_id > 0.
      3. Emit one row per (barcode, cell_id) with fractional weight = count / total_foreground_pixels
         (background pixels with cell_id == 0 are excluded from the denominator).

    Returns DataFrame: barcode, cell_id, weight (float, sums to 1 per barcode).
    """
    if wba_csv.exists():
        log.info(f"Step 3: 載入已存在的 WBA attribution: {wba_csv.name}")
        return pd.read_csv(wba_csv)

    log.info("Step 3: WBA attribution (7×7 footprint per bin)")
    t0 = time.time()

    tp = pd.read_parquet(str(TP_PATH), columns=[
        "barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    tp = tp[tp["in_tissue"] == 1]

    in_crop = (
        (tp["pxl_col_in_fullres"] >= CROP_X0 + COL_OFFSET) &
        (tp["pxl_col_in_fullres"] <  CROP_X1 + COL_OFFSET) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0) &
        (tp["pxl_row_in_fullres"] <  CROP_Y1)
    )
    tp_crop = tp[in_crop].copy().reset_index(drop=True)
    log.info(f"  bins in crop: {len(tp_crop):,}")

    row_c = (tp_crop["pxl_row_in_fullres"].values - CROP_Y0).astype(np.int32)
    col_c = (tp_crop["pxl_col_in_fullres"].values - COL_OFFSET - CROP_X0).astype(np.int32)
    H, W  = mask.shape

    barcodes_out: list[str] = []
    cell_ids_out: list[int] = []
    weights_out:  list[float] = []

    n_boundary = 0
    CHUNK = 200_000

    for start in range(0, len(tp_crop), CHUNK):
        end = min(start + CHUNK, len(tp_crop))
        barcodes = tp_crop["barcode"].values[start:end]
        rows_b   = row_c[start:end]
        cols_b   = col_c[start:end]

        for i in range(end - start):
            r, c = rows_b[i], cols_b[i]
            r0 = max(0, r - BIN_HALF)
            r1 = min(H, r + BIN_HALF + 1)
            c0 = max(0, c - BIN_HALF)
            c1 = min(W, c + BIN_HALF + 1)

            if (r1 - r0) < (2 * BIN_HALF + 1) or (c1 - c0) < (2 * BIN_HALF + 1):
                n_boundary += 1

            window = mask[r0:r1, c0:c1]
            flat   = window.ravel()
            ids, counts = np.unique(flat[flat > 0], return_counts=True)

            if len(ids) == 0:
                continue

            total = float(counts.sum())
            bc    = barcodes[i]
            for cid, cnt in zip(ids, counts):
                barcodes_out.append(bc)
                cell_ids_out.append(int(cid))
                weights_out.append(cnt / total)

        if (start // CHUNK) % 10 == 0:
            log.info(
                f"  processed {end:,}/{len(tp_crop):,} bins  "
                f"({end/len(tp_crop):.0%})  rows so far: {len(barcodes_out):,}"
            )

    log.info(f"  boundary bins (partial window): {n_boundary:,}")
    log.info(f"  total (barcode, cell_id) pairs: {len(barcodes_out):,}")
    log.info(f"  elapsed: {time.time()-t0:.0f}s")

    wba = pd.DataFrame({
        "barcode": barcodes_out,
        "cell_id": cell_ids_out,
        "weight":  weights_out,
    })
    wba.to_csv(str(wba_csv), index=False)
    log.info(f"  儲存: {wba_csv.name}")
    return wba


# ─── Step 4: AnnData + CellTypist with WBA ────────────────────────────────────

def build_anndata_wba(wba: pd.DataFrame, celltypist_csv: Path) -> pd.DataFrame:
    """
    Weighted pseudo-bulk aggregation: for each cell_id, sum transcripts
    weighted by WBA weight.  Then CellTypist annotate.
    """
    if celltypist_csv.exists():
        log.info(f"Step 4: 載入已存在的 CellTypist WBA 結果: {celltypist_csv.name}")
        return pd.read_csv(celltypist_csv)

    log.info("Step 4: 建立 WBA AnnData + CellTypist 標注")
    import scipy.sparse as sp

    log.info(f"  讀取 h5 矩陣: {H5_PATH.name}")
    adata_full = sc.read_10x_h5(str(H5_PATH))
    adata_full.var_names_make_unique()
    log.info(f"  全片: {adata_full.n_obs:,} barcodes × {adata_full.n_vars:,} genes")

    barcodes_in_wba = wba["barcode"].unique()
    mask_obs = adata_full.obs_names.isin(barcodes_in_wba)
    adata_sub = adata_full[mask_obs].copy()
    del adata_full
    gc.collect()
    log.info(f"  barcodes in WBA matched in h5: {adata_sub.n_obs:,}")

    # barcode → index in adata_sub
    bc_to_idx = {bc: i for i, bc in enumerate(adata_sub.obs_names)}

    # unique cells
    unique_cells = np.sort(wba["cell_id"].unique()).astype(np.int32)
    n_cells      = len(unique_cells)
    cell_to_idx  = {int(c): i for i, c in enumerate(unique_cells)}
    log.info(f"  unique cells with WBA attribution: {n_cells:,}")

    # Build sparse weight matrix W (n_cells × n_barcodes_in_sub)
    # Only include rows where barcode appears in adata_sub
    wba_valid = wba[wba["barcode"].isin(bc_to_idx)].copy()
    row_idx = wba_valid["cell_id"].map(cell_to_idx).values.astype(np.int32)
    col_idx = wba_valid["barcode"].map(bc_to_idx).values.astype(np.int32)
    vals    = wba_valid["weight"].values.astype(np.float32)

    W = sp.csr_matrix(
        (vals, (row_idx, col_idx)),
        shape=(n_cells, adata_sub.n_obs),
    )
    log.info(f"  W matrix: {W.shape}  nnz: {W.nnz:,}")

    X = adata_sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X_agg = W @ X   # n_cells × n_genes

    adata_cells = sc.AnnData(
        X=X_agg.tocsr() if sp.issparse(X_agg) else sp.csr_matrix(X_agg),
        var=adata_sub.var.copy(),
    )
    adata_cells.obs_names = [str(c) for c in unique_cells]
    del adata_sub, W, X, X_agg
    gc.collect()

    sc.pp.normalize_total(adata_cells, target_sum=1e4)
    sc.pp.log1p(adata_cells)
    log.info(f"  AnnData: {adata_cells.n_obs:,} cells × {adata_cells.n_vars:,} genes")

    import celltypist
    log.info("  CellTypist 標注 (Human_Colorectal_Cancer.pkl, majority_voting=False)...")
    predictions = celltypist.annotate(
        adata_cells,
        model="Human_Colorectal_Cancer.pkl",
        majority_voting=False,
    )
    ct_labels = predictions.predicted_labels["predicted_labels"].values

    df = pd.DataFrame({
        "cell_id":          unique_cells,
        "celltypist_label": ct_labels,
        "broad_label":      [LABEL_MAP.get(lbl, "other") for lbl in ct_labels],
    })

    unmapped = df[df["broad_label"] == "other"]["celltypist_label"].value_counts()
    if len(unmapped) > 0:
        log.warning(
            f"  未映射標籤 ({len(unmapped)} types, "
            f"{df['broad_label'].eq('other').sum()} cells):\n{unmapped.to_string()}"
        )

    df.to_csv(str(celltypist_csv), index=False)
    log.info(f"  儲存: {celltypist_csv.name}")
    log.info(f"  broad_label 分佈:\n{df['broad_label'].value_counts().to_string()}")
    return df


# ─── Step 5: Spatial matching (same as C1) ────────────────────────────────────

def match_gt_to_mcseg(
    gt: pd.DataFrame,
    mask: np.ndarray,
    cell_labels: pd.DataFrame,
) -> pd.DataFrame:
    log.info("Step 5: 空間匹配 (GT centroid → MCseg mask lookup)")
    cell_id_to_label = cell_labels.set_index("cell_id")["broad_label"].to_dict()
    H, W = mask.shape

    cx, cy = gt["cell_x"].values, gt["cell_y"].values
    in_bounds = (
        (cx >= CROP_X0) & (cx < CROP_X1) &
        (cy >= CROP_Y0) & (cy < CROP_Y1)
    )
    row_local = np.where(in_bounds, (cy - CROP_Y0).astype(np.int32).clip(0, H - 1), 0)
    col_local = np.where(in_bounds, (cx - CROP_X0).astype(np.int32).clip(0, W - 1), 0)
    cell_ids  = np.where(in_bounds, mask[row_local, col_local], 0)
    pred_labels = np.array([
        cell_id_to_label.get(int(cid), "unmatched") if cid > 0 else "unmatched"
        for cid in cell_ids
    ])

    gt = gt.copy()
    gt["mcseg_cell_id"] = cell_ids
    gt["pred_label"]    = pred_labels

    n_matched = (pred_labels != "unmatched").sum()
    log.info(f"  GT cells matched: {n_matched:,} / {len(gt):,} ({n_matched/len(gt):.1%})")
    return gt


# ─── Step 6: F1 + 3-way comparison figure ─────────────────────────────────────

def compute_f1_with_comparison(gt_matched: pd.DataFrame) -> dict:
    from sklearn.metrics import f1_score, precision_recall_fscore_support, classification_report
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    log.info("Step 6: 計算 F1 + 3-way comparison figure")

    matched = gt_matched[gt_matched["pred_label"] != "unmatched"]
    valid   = matched[matched["pred_label"] != "other"]
    n_total, n_matched, n_valid = len(gt_matched), len(matched), len(valid)

    y_true  = valid["gt_label"].values
    y_pred  = valid["pred_label"].values
    classes = ["epithelial cells", "stromal cells", "immune cells"]

    micro_f1    = f1_score(y_true, y_pred, labels=classes, average="micro",    zero_division=0)
    weighted_f1 = f1_score(y_true, y_pred, labels=classes, average="weighted", zero_division=0)
    macro_f1    = f1_score(y_true, y_pred, labels=classes, average="macro",    zero_division=0)

    log.info(f"  WBA micro F1    = {micro_f1:.3f}  (StarDist: {STARDIST_MICRO_F1})")
    log.info(f"  WBA weighted F1 = {weighted_f1:.3f}  (StarDist: {STARDIST_WEIGHTED_F1})")
    log.info(f"  WBA macro F1    = {macro_f1:.3f}")
    log.info(f"\n{classification_report(y_true, y_pred, labels=classes, zero_division=0)}")

    prec, rec, f1_per, sup = precision_recall_fscore_support(
        y_true, y_pred, labels=classes, zero_division=0
    )
    pd.DataFrame({
        "class": classes, "precision": prec, "recall": rec, "f1": f1_per, "support": sup,
    }).to_csv(str(RESULTS_DIR / "f1_per_class_wba.csv"), index=False)

    summary = {
        "n_gt_total": n_total, "n_matched": n_matched, "n_valid": n_valid,
        "match_rate": n_matched / n_total, "valid_rate": n_valid / n_total,
        "micro_f1": micro_f1, "weighted_f1": weighted_f1, "macro_f1": macro_f1,
        "stardist_micro_f1": STARDIST_MICRO_F1, "stardist_weighted_f1": STARDIST_WEIGHTED_F1,
        "delta_micro": micro_f1 - STARDIST_MICRO_F1,
        "delta_weighted": weighted_f1 - STARDIST_WEIGHTED_F1,
    }
    pd.DataFrame([summary]).to_csv(str(RESULTS_DIR / "f1_summary_wba.csv"), index=False)

    # Load C1 results for comparison
    c1_summary_path = C1_DIR / "f1_summary.csv"
    c1_micro, c1_weighted = STARDIST_MICRO_F1, STARDIST_WEIGHTED_F1  # fallback
    if c1_summary_path.exists():
        c1 = pd.read_csv(c1_summary_path).iloc[0]
        c1_micro    = float(c1["micro_f1"])
        c1_weighted = float(c1["weighted_f1"])

    # 3-way comparison bar chart
    methods   = ["StarDist\n(ENACT WBA)", "MCseg\n+lookup", "MCseg\n+WBA"]
    micro_v   = [STARDIST_MICRO_F1,    c1_micro,    micro_f1]
    weighted_v = [STARDIST_WEIGHTED_F1, c1_weighted, weighted_f1]
    colours_m = ["#4e79a7", "#f28e2b", "#59a14f"]
    colours_w = ["#76b7b2", "#e15759", "#b07aa1"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4))

    ax = axes[0]
    x, w = np.arange(3), 0.35
    bars1 = ax.bar(x - w / 2, micro_v,    w, label="Micro F1",    color=colours_m)
    bars2 = ax.bar(x + w / 2, weighted_v, w, label="Weighted F1", color=colours_w)
    ax.set_ylim(0.60, 0.88)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=9)
    ax.set_ylabel("F1 Score")
    ax.set_title("Cell-type Annotation F1 — 3-way Comparison\n(ENACT CRC Expert GT, n=20,991)")
    ax.legend(fontsize=8)
    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.003, f"{h:.3f}",
                ha="center", va="bottom", fontsize=8)
    ax.axhline(STARDIST_MICRO_F1, ls="--", color="steelblue", lw=0.8, alpha=0.5)

    # Confusion matrix for WBA
    from sklearn.metrics import confusion_matrix
    cm      = confusion_matrix(y_true, y_pred, labels=classes)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    short   = ["Epithelial", "Stromal", "Immune"]
    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Greens",
                xticklabels=short, yticklabels=short, ax=axes[1],
                cbar_kws={"label": "Proportion"})
    axes[1].set_xlabel("Predicted")
    axes[1].set_ylabel("True (GT)")
    axes[1].set_title(f"MCseg+WBA Confusion Matrix (Micro F1={micro_f1:.3f})")

    plt.tight_layout()
    out_fig = RESULTS_DIR / "fig_f1_comparison_3way.png"
    fig.savefig(str(out_fig), dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  Figure saved: {out_fig.name}")

    return summary


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_total = time.time()

    gt         = load_gt()
    mask       = load_mask()
    wba        = run_wba_attribution(mask, RESULTS_DIR / "wba_attribution.csv")
    cell_labels = build_anndata_wba(wba, RESULTS_DIR / "celltypist_labels_wba.csv")
    gt_matched  = match_gt_to_mcseg(gt, mask, cell_labels)
    gt_matched.to_csv(str(RESULTS_DIR / "gt_matched_wba.csv"), index=False)
    summary     = compute_f1_with_comparison(gt_matched)

    elapsed = time.time() - t_total
    log.info(f"\n{'='*60}")
    log.info(f"完成 ({elapsed/60:.1f} min)  結果: {RESULTS_DIR}")
    log.info(f"MCseg+WBA  micro F1    = {summary['micro_f1']:.3f}")
    log.info(f"MCseg+WBA  weighted F1 = {summary['weighted_f1']:.3f}")
    log.info(f"Δ vs StarDist: micro {summary['delta_micro']:+.3f}  weighted {summary['delta_weighted']:+.3f}")
    log.info("="*60)


if __name__ == "__main__":
    main()
