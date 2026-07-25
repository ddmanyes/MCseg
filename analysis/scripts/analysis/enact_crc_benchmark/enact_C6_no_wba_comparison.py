"""
enact_C6_no_wba_comparison.py
Compare +WBA vs no-WBA (simple direct sum) for each segmentation method.

WBA  : weight = 1/n_bins  (average bin expression per cell)
noWBA: weight = 1         (sum all bins per cell, normalize_total handles scale)

Since CellTypist always runs normalize_total (10k) + log1p before annotation,
the difference between sum vs average should be minimal for cells with adequate
RNA. This script empirically verifies that claim.

Reads pre-existing masks from:
  results/enact_method_comparison/{stardist,mcseg,sr,cellpose_nuc,proseg}_mask.npy
"""

from __future__ import annotations
import logging
import time
from pathlib import Path

import anndata
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import f1_score, classification_report

# ── Paths (same as C3) ────────────────────────────────────────────────────────
PLAN_A       = Path("/Volumes/SSD/plan_a")
ENACT_CRC    = PLAN_A / "tissue sample" / "ENACT_supporting_files" / "public_data" / "human_colorectal"
H5_PATH      = ENACT_CRC / "input_files" / "filtered_feature_bc_matrix.h5"
TP_PATH      = ENACT_CRC / "input_files" / "tissue_positions.parquet"
CELLTYPIST_EVAL = (
    ENACT_CRC / "paper_results" / "chunks" / "weighted_by_area"
    / "celltypist_results" / "eval" / "cell_annotation_eval.csv"
)
RESULT_DIR   = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
ENACT_F1     = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
MCSEG_MASK   = ENACT_F1 / "mcseg_mask_7pass.npy"

CELLTYPIST_MODEL = "Human_Colorectal_Cancer.pkl"

COL_OFFSET        = 40598
CROP_X0, CROP_X1  = 5154, 15242
CROP_Y0, CROP_Y1  = 4635, 18599
CROP_H = CROP_Y1 - CROP_Y0
CROP_W = CROP_X1 - CROP_X0
BTF_COL_MIN = CROP_X0 + COL_OFFSET
BTF_COL_MAX = CROP_X1 + COL_OFFSET
BTF_ROW_MIN = CROP_Y0
BTF_ROW_MAX = CROP_Y1

CELLTYPIST_BROAD: dict[str, str] = {
    "CMS1": "epithelial cells", "CMS2": "epithelial cells",
    "CMS3": "epithelial cells", "Goblet cells": "epithelial cells",
    "Intermediate": "epithelial cells", "Proliferating": "epithelial cells",
    "Stem-like/TA": "epithelial cells",
    "Myofibroblasts": "stromal cells", "Lymphatic ECs": "stromal cells",
    "Pericytes": "stromal cells", "Smooth muscle cells": "stromal cells",
    "Stalk-like ECs": "stromal cells", "Enteric glial cells": "stromal cells",
    "Stromal 1": "stromal cells", "Stromal 2": "stromal cells",
    "Stromal 3": "stromal cells", "Proliferative ECs": "stromal cells",
    "CD4+ T cells": "immune cells", "CD8+ T cells": "immune cells",
    "CD19+CD20+ B": "immune cells", "IgA+ Plasma": "immune cells",
    "IgG+ Plasma": "immune cells", "NK cells": "immune cells",
    "Pro-inflammatory": "immune cells", "Regulatory T cells": "immune cells",
    "SPP1+": "immune cells", "Mast cells": "immune cells",
    "T follicular helper cells": "other", "T helper 17 cells": "other",
    "Tip-like ECs": "other", "Unknown": "other",
    "cDC": "other", "gamma delta T cells": "other",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_tissue_positions() -> pd.DataFrame:
    tp = pd.read_parquet(str(TP_PATH), columns=[
        "barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    in_crop = (
        (tp["pxl_col_in_fullres"] >= BTF_COL_MIN) &
        (tp["pxl_col_in_fullres"] <  BTF_COL_MAX) &
        (tp["pxl_row_in_fullres"] >= BTF_ROW_MIN) &
        (tp["pxl_row_in_fullres"] <  BTF_ROW_MAX)
    )
    tp_crop = tp[in_crop].copy()
    tp_crop["row_c"] = (tp_crop["pxl_row_in_fullres"] - CROP_Y0).astype(np.int32)
    tp_crop["col_c"] = (tp_crop["pxl_col_in_fullres"] - BTF_COL_MIN).astype(np.int32)
    log.info(f"  tissue_positions in crop: {len(tp_crop):,}")
    return tp_crop


def load_adata(tp_crop: pd.DataFrame) -> sc.AnnData:
    adata_full = sc.read_10x_h5(str(H5_PATH))
    adata_full.var_names_make_unique()
    crop_barcodes = set(tp_crop["barcode"].values)
    keep = np.array([b in crop_barcodes for b in adata_full.obs_names])
    adata = adata_full[keep].copy()
    log.info(f"  adata shape: {adata.shape}")
    return adata


def load_gt_cells() -> pd.DataFrame:
    df = pd.read_csv(str(CELLTYPIST_EVAL), usecols=["cell_x", "cell_y", "gt_label"])
    df["cell_x"] = df["cell_x"] - CROP_X0
    df["cell_y"] = df["cell_y"] - CROP_Y0
    log.info(f"  GT cells: {len(df):,}")
    return df


# ── Attribution ───────────────────────────────────────────────────────────────

def compute_attribution(mask: np.ndarray, tp_crop: pd.DataFrame, use_wba: bool) -> pd.DataFrame:
    """
    use_wba=True  → weight = 1/n_bins  (average expression per cell)
    use_wba=False → weight = 1.0       (sum expression; normalize_total handles scale)
    """
    row_c = tp_crop["row_c"].values.clip(0, CROP_H - 1)
    col_c = tp_crop["col_c"].values.clip(0, CROP_W - 1)
    cell_ids = mask[row_c, col_c]
    df = pd.DataFrame({"barcode": tp_crop["barcode"].values, "cell_id": cell_ids})
    df = df[df["cell_id"] > 0]
    if use_wba:
        n_bins = df.groupby("cell_id").size().rename("n_bins")
        df = df.join(n_bins, on="cell_id")
        df["weight"] = 1.0 / df["n_bins"]
    else:
        df["weight"] = 1.0
    return df[["barcode", "cell_id", "weight"]]


def aggregate_expr_to_adata(attr: pd.DataFrame, adata: sc.AnnData) -> anndata.AnnData:
    import scipy.sparse as sp
    barcode_to_idx = {b: i for i, b in enumerate(adata.obs_names)}
    src_idxs = np.array([barcode_to_idx.get(b, -1) for b in attr["barcode"].values])
    valid    = src_idxs >= 0
    src_v    = src_idxs[valid]
    cell_v   = attr["cell_id"].values[valid]
    weight_v = attr["weight"].values[valid].astype(np.float32)
    unique_cells = np.unique(cell_v)
    cell_to_row  = {c: i for i, c in enumerate(unique_cells)}
    row_v = np.array([cell_to_row[c] for c in cell_v], dtype=np.int32)
    W = sp.csr_matrix(
        (weight_v, (row_v, src_v)),
        shape=(len(unique_cells), adata.shape[0]),
        dtype=np.float32,
    )
    X_src = adata.X
    if not sp.issparse(X_src):
        X_src = sp.csr_matrix(X_src)
    X_cell = W @ X_src
    return anndata.AnnData(
        X=X_cell.astype(np.float32),
        obs=pd.DataFrame(index=[str(c) for c in unique_cells]),
        var=adata.var.copy(),
    )


def celltypist_annotate(cell_adata: anndata.AnnData) -> pd.DataFrame:
    sc.pp.normalize_total(cell_adata, target_sum=1e4)
    sc.pp.log1p(cell_adata)
    predictions = celltypist.annotate(cell_adata, model=CELLTYPIST_MODEL, majority_voting=False)
    pred_df = predictions.predicted_labels.copy()
    pred_df.index = cell_adata.obs_names
    pred_df["broad_label"] = pred_df["predicted_labels"].map(CELLTYPIST_BROAD).fillna("other")
    return pd.DataFrame({
        "cell_id":    pred_df.index.astype(int).values,
        "pred_label": pred_df["broad_label"].values,
    })


def match_gt(gt: pd.DataFrame, mask: np.ndarray, annotations: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    H, W = mask.shape
    ann_map = annotations.set_index("cell_id")["pred_label"].to_dict()
    gt = gt.dropna(subset=["cell_x", "cell_y"]).reset_index(drop=True)
    row = gt["cell_y"].round().astype(int).clip(0, H - 1).values
    col = gt["cell_x"].round().astype(int).clip(0, W - 1).values
    cell_ids = mask[row, col]
    pred = [ann_map.get(int(cid), "unmatched") if cid > 0 else "unmatched" for cid in cell_ids]
    result = gt.copy()
    result["pred_label"] = pred
    n_matched = int((np.array(cell_ids) > 0).sum())
    return result, n_matched


def compute_record(matched: pd.DataFrame, n_matched: int, n_gt: int, name: str) -> dict:
    classes = ["epithelial cells", "stromal cells", "immune cells"]
    wf1 = f1_score(matched["gt_label"], matched["pred_label"],
                   labels=classes, average="weighted", zero_division=0)
    matched_only = matched[matched["pred_label"] != "unmatched"]
    subset_wf1 = f1_score(
        matched_only["gt_label"], matched_only["pred_label"],
        labels=classes, average="weighted", zero_division=0,
    )
    report = classification_report(matched["gt_label"], matched["pred_label"],
                                   labels=classes, output_dict=True, zero_division=0)
    match_rate = n_matched / n_gt
    log.info(f"  match={match_rate*100:.1f}%  subset_wF1={subset_wf1:.3f}  wF1={wf1:.3f}")
    return {
        "method": name,
        "match_rate": match_rate,
        "subset_wf1": subset_wf1,
        "weighted_f1": wf1,
        **{
            f"{cls.split()[0]}_{k}": report.get(cls, {}).get(metric, 0)
            for cls in classes
            for k, metric in [("p", "precision"), ("r", "recall"), ("f1", "f1-score")]
        },
    }


def eval_mask(mask: np.ndarray, tp_crop: pd.DataFrame, adata: sc.AnnData,
              gt: pd.DataFrame, name: str, use_wba: bool) -> dict:
    log.info(f"=== {name} ===")
    t0 = time.time()
    attr        = compute_attribution(mask, tp_crop, use_wba=use_wba)
    cell_adata  = aggregate_expr_to_adata(attr, adata)
    annotations = celltypist_annotate(cell_adata)
    matched, n_matched = match_gt(gt, mask, annotations)
    record = compute_record(matched, n_matched, len(gt), name)
    log.info(f"  Done in {time.time()-t0:.0f}s")
    return record


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    log.info("Loading shared resources...")
    tp_crop = load_tissue_positions()
    adata   = load_adata(tp_crop)
    gt      = load_gt_cells()

    masks: dict[str, Path] = {
        "StarDist": RESULT_DIR / "stardist_mask.npy",
        "MCseg":    MCSEG_MASK,
        "SR":       RESULT_DIR / "sr_mask.npy",
        "NUC":      RESULT_DIR / "cellpose_nuc_mask.npy",
        "ProSeg":   RESULT_DIR / "proseg_mask.npy",
    }

    records = []
    for method, mask_path in masks.items():
        if not mask_path.exists():
            log.warning(f"  {method}: mask not found — skipping")
            continue
        mask = np.load(str(mask_path))
        for use_wba in [True, False]:
            tag = f"{method}+WBA" if use_wba else f"{method} (no WBA)"
            records.append(eval_mask(mask, tp_crop, adata, gt, tag, use_wba))
        del mask

    df = pd.DataFrame(records)
    out_csv = RESULT_DIR / "no_wba_comparison.csv"
    df.to_csv(out_csv, index=False)
    log.info(f"\nSaved: {out_csv}")

    # ── Results table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print(f"{'Method':<24}  {'match%':>7}  {'subset_wF1':>10}  {'weighted_F1':>11}")
    print("-" * 65)
    for method in masks:
        for use_wba in [True, False]:
            tag = f"{method}+WBA" if use_wba else f"{method} (no WBA)"
            row = df[df["method"] == tag]
            if row.empty:
                continue
            r = row.iloc[0]
            print(f"{tag:<24}  {r['match_rate']*100:>6.1f}%  {r['subset_wf1']:>10.3f}  {r['weighted_f1']:>11.3f}")
        print()

    # ── Delta table ───────────────────────────────────────────────────────────
    print("=" * 65)
    print("Delta: +WBA minus no-WBA  (positive = WBA better)")
    print("-" * 65)
    print(f"{'Method':<12}  {'Δmatch%':>8}  {'Δsubset_wF1':>12}  {'Δweighted_F1':>13}")
    print("-" * 65)
    for method in masks:
        wba_row   = df[df["method"] == f"{method}+WBA"]
        nowba_row = df[df["method"] == f"{method} (no WBA)"]
        if wba_row.empty or nowba_row.empty:
            continue
        w = wba_row.iloc[0]
        n = nowba_row.iloc[0]
        dm   = (w["match_rate"] - n["match_rate"]) * 100
        dwf1 = w["subset_wf1"] - n["subset_wf1"]
        df1  = w["weighted_f1"] - n["weighted_f1"]
        print(f"{method:<12}  {dm:>+8.1f}  {dwf1:>+12.3f}  {df1:>+13.3f}")


if __name__ == "__main__":
    main()
