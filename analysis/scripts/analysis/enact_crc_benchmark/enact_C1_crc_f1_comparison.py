"""
enact_C1_crc_f1_comparison.py
==============================
Phase 1: MCseg vs StarDist F1 比較 (ENACT CRC expert GT)

流程:
  Step 1  讀 cell_annotation_eval.csv (GT 細胞, gt_label, cell_x/y)
  Step 2  從 CRC BTF crop 出 H&E patch (tifffile region 讀取)
  Step 3  MCseg 分割 (deployment-mode, Voronoi d=8, tiled)
  Step 4  Bin attribution (tissue_positions + mask lookup)
  Step 5  AnnData 建立 + CellTypist 細胞型態標注
  Step 6  空間匹配: GT centroid → MCseg mask lookup
  Step 7  計算 F1 (micro + weighted) vs StarDist baseline 0.708/0.758

輸出: submission_bioinformatics/results/enact_crc_f1/
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
import tifffile
import zarr

# ─── Paths ────────────────────────────────────────────────────────────────────

PLAN_A    = Path("/Volumes/SSD/plan_a")
ENACT_CRC = PLAN_A / "tissue sample" / "ENACT_supporting_files" / "public_data" / "human_colorectal"
BTF_PATH  = ENACT_CRC / "input_files" / "Visium_HD_Human_Colon_Cancer_tissue_image.btf"
H5_PATH   = ENACT_CRC / "input_files" / "filtered_feature_bc_matrix.h5"
TP_PATH   = ENACT_CRC / "input_files" / "tissue_positions.parquet"
EVAL_CSV  = ENACT_CRC / "paper_results" / "chunks" / "weighted_by_area" / "sargent_results" / "eval" / "cell_annotation_eval.csv"

MSSEG_ROOT  = PLAN_A / "MSseg"
RESULTS_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"

# GT coordinate range in ENACT-local BTF pixels (= cell_x/cell_y) + 200px buffer
CROP_X0, CROP_X1 = 5154, 15242   # ENACT local col  (= BTF col)
CROP_Y0, CROP_Y1 = 4635, 18599   # ENACT local row  (= pxl_row_in_fullres directly)

# tissue_positions.pxl_col_in_fullres = ENACT_local_col + COL_OFFSET
# Empirically verified: offset = 40598.0 exactly (std=0), row offset = 0
COL_OFFSET = 40598

# ENACT StarDist baseline (confirmed from local eval files)
STARDIST_MICRO_F1    = 0.708
STARDIST_WEIGHTED_F1 = 0.758

# CellTypist Human_Colorectal_Cancer.pkl → 3 ENACT classes
LABEL_MAP: dict[str, str] = {
    # epithelial
    "CMS1":                      "epithelial cells",
    "CMS2":                      "epithelial cells",
    "CMS3":                      "epithelial cells",
    "CMS4":                      "epithelial cells",
    "Goblet cells":               "epithelial cells",
    "Mature Enterocytes type 1":  "epithelial cells",
    "Mature Enterocytes type 2":  "epithelial cells",
    "Stem-like/TA":               "epithelial cells",
    "Intermediate":               "epithelial cells",
    "Proliferating":              "epithelial cells",
    # stromal
    "Myofibroblasts":     "stromal cells",
    "Pericytes":          "stromal cells",
    "Smooth muscle cells":"stromal cells",
    "Stromal 1":          "stromal cells",
    "Stromal 2":          "stromal cells",
    "Stromal 3":          "stromal cells",
    "Lymphatic ECs":      "stromal cells",
    "Proliferative ECs":  "stromal cells",
    "Stalk-like ECs":     "stromal cells",
    "Enteric glial cells":"stromal cells",
    # immune
    "CD19+CD20+ B":       "immune cells",
    "CD4+ T cells":       "immune cells",
    "CD8+ T cells":       "immune cells",
    "Regulatory T cells": "immune cells",
    "NK cells":           "immune cells",
    "IgA+ Plasma":        "immune cells",
    "IgG+ Plasma":        "immune cells",
    "Mast cells":         "immune cells",
    "Pro-inflammatory":   "immune cells",
    "SPP1+":              "immune cells",
}

# ─── Logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── Step 1: Load GT ──────────────────────────────────────────────────────────

def load_gt(eval_csv: Path) -> pd.DataFrame:
    log.info("Step 1: 載入 ENACT expert GT")
    gt = pd.read_csv(eval_csv)
    log.info(f"  GT cells: {len(gt):,}")
    log.info(f"  gt_label distribution:\n{gt['gt_label'].value_counts().to_string()}")
    gt = gt[["cell_x", "cell_y", "gt_label"]].copy()
    gt["cell_x"] = gt["cell_x"].astype(np.float32)
    gt["cell_y"] = gt["cell_y"].astype(np.float32)
    return gt


# ─── Step 2: Crop H&E from BTF ────────────────────────────────────────────────

def crop_he_from_btf(btf_path: Path, crop_tif: Path) -> np.ndarray:
    """讀取 BTF 指定區域，儲存為 crop_tif，回傳 RGB array。"""
    if crop_tif.exists():
        log.info(f"Step 2: 載入已存在的 H&E crop: {crop_tif.name}")
        img = tifffile.imread(str(crop_tif))
        if img.ndim == 3 and img.shape[-1] == 4:
            img = img[..., :3]
        log.info(f"  crop shape: {img.shape}")
        return img

    btf_col0 = CROP_X0 + COL_OFFSET
    btf_col1 = CROP_X1 + COL_OFFSET
    log.info(f"Step 2: 從 BTF 讀取 H&E crop (row {CROP_Y0}:{CROP_Y1}, BTF col {btf_col0}:{btf_col1})")
    t0 = time.time()
    # Use zarr lazy region-read to avoid loading the full 10 GB BTF into RAM
    with tifffile.TiffFile(str(btf_path)) as tif:
        store = tif.aszarr()
    z = zarr.open(store, mode="r")
    # zarr store shape is (H, W, C) for single-level or (levels, H, W, C) for pyramid
    arr = z[0] if z.ndim == 4 else z
    img = np.asarray(arr[CROP_Y0:CROP_Y1, btf_col0:btf_col1])
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  crop shape: {img.shape}  ({time.time()-t0:.0f}s)")

    tifffile.imwrite(str(crop_tif), img, compression="zlib")
    log.info(f"  儲存: {crop_tif.name}")
    return img


# ─── Step 3: MCseg Segmentation ───────────────────────────────────────────────

def run_mcseg(img: np.ndarray, mask_npy: Path) -> np.ndarray:
    """Deployment-mode MCseg v2 (voronoi_distance=8, no GT guidance)."""
    if mask_npy.exists():
        log.info(f"Step 3: 載入已存在的 MCseg 遮罩: {mask_npy.name}")
        mask = np.load(str(mask_npy))
        log.info(f"  mask shape: {mask.shape}  cells: {int(mask.max()):,}")
        return mask

    sys.path.insert(0, str(MSSEG_ROOT / "backend"))
    from src.segmentation.cellpose_runner import run_tiled_mcseg_v2

    cfg = {
        "use_gpu":               True,
        "batch_size":            2,
        "dia_small":             13.0,
        "dia_mid":               17.0,
        "dia_large":             22.0,
        "use_hematoxylin":       True,
        "use_cpsam":             False,
        "voronoi_distance":      8,      # deployment mode
        "flow_threshold":        0.4,
        "cellprob_threshold":    -2.0,
        "min_size":              20,
        "max_size":              6000,
        "clahe_clip_limit":      3.0,
        "use_transcript_rescue": False,
    }

    log.info("Step 3: MCseg v2 tiled 分割 (deployment-mode, voronoi_d=8)")

    def _progress(p: float, msg: str) -> None:
        log.info(f"  [{p*100:.0f}%] {msg}")

    mask = run_tiled_mcseg_v2(
        img,
        cfg,
        tile_size=1024,
        overlap=128,
        progress_callback=_progress,
    )
    np.save(str(mask_npy), mask)
    log.info(f"  儲存遮罩: {mask_npy.name}  cells: {int(mask.max()):,}")
    return mask


# ─── Step 4: Bin Attribution ──────────────────────────────────────────────────

def run_bin_attribution(mask: np.ndarray, tp_path: Path) -> pd.DataFrame:
    """
    對 crop 範圍內的 Visium HD 2µm bins 做 mask lookup。
    回傳 DataFrame: barcode, cell_id
    """
    log.info("Step 4: Bin attribution (mask lookup)")
    tp = pd.read_parquet(str(tp_path), columns=[
        "barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    tp = tp[tp["in_tissue"] == 1]

    # Convert pxl_col → ENACT local col: enact_col = pxl_col - COL_OFFSET
    # pxl_row == ENACT local row directly (row offset = 0)
    in_crop = (
        (tp["pxl_col_in_fullres"] >= CROP_X0 + COL_OFFSET) &
        (tp["pxl_col_in_fullres"] <  CROP_X1 + COL_OFFSET) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0) &
        (tp["pxl_row_in_fullres"] <  CROP_Y1)
    )
    tp_crop = tp[in_crop].copy()
    log.info(f"  bins in crop: {len(tp_crop):,}")

    row_local = (tp_crop["pxl_row_in_fullres"].values - CROP_Y0).astype(np.int32)
    col_local = (tp_crop["pxl_col_in_fullres"].values - COL_OFFSET - CROP_X0).astype(np.int32)
    H, W = mask.shape
    row_local = row_local.clip(0, H - 1)
    col_local = col_local.clip(0, W - 1)

    tp_crop = tp_crop.copy()
    tp_crop["cell_id"] = mask[row_local, col_local]

    attributed = tp_crop[tp_crop["cell_id"] > 0]
    log.info(f"  attributed bins: {len(attributed):,} ({len(attributed)/len(tp_crop):.1%})")
    return attributed[["barcode", "cell_id"]].reset_index(drop=True)


# ─── Step 5: AnnData + CellTypist ─────────────────────────────────────────────

def build_anndata_and_annotate(
    attribution: pd.DataFrame,
    h5_path: Path,
    celltypist_csv: Path,
) -> pd.DataFrame:
    """
    建立 cell × gene AnnData，執行 CellTypist，
    回傳 DataFrame: cell_id, celltypist_label, broad_label
    """
    if celltypist_csv.exists():
        log.info(f"Step 5: 載入已存在的 CellTypist 結果: {celltypist_csv.name}")
        return pd.read_csv(celltypist_csv)

    log.info("Step 5: 建立 AnnData + CellTypist 標注")

    log.info(f"  讀取 h5 矩陣: {h5_path.name}")
    adata_full = sc.read_10x_h5(str(h5_path))
    adata_full.var_names_make_unique()
    log.info(f"  全片: {adata_full.n_obs:,} barcodes × {adata_full.n_vars:,} genes")

    barcodes_in_crop = attribution["barcode"].values
    mask_obs = adata_full.obs_names.isin(barcodes_in_crop)
    adata_crop = adata_full[mask_obs].copy()
    del adata_full
    gc.collect()
    log.info(f"  crop barcodes matched: {adata_crop.n_obs:,}")

    barcode_to_cell = attribution.set_index("barcode")["cell_id"]
    adata_crop.obs["cell_id"] = barcode_to_cell.reindex(adata_crop.obs_names).values

    import scipy.sparse as sp
    cell_ids = adata_crop.obs["cell_id"].values.astype(np.int32)
    valid = cell_ids > 0
    adata_valid = adata_crop[valid]
    cell_ids_v = cell_ids[valid]
    unique_cells = np.unique(cell_ids_v)
    n_cells = len(unique_cells)
    log.info(f"  unique cells with RNA: {n_cells:,}")

    cell_id_to_idx = {int(c): i for i, c in enumerate(unique_cells)}
    rows = np.array([cell_id_to_idx[int(c)] for c in cell_ids_v])
    cols = np.arange(len(cell_ids_v))
    A = sp.csr_matrix(
        (np.ones(len(cell_ids_v), dtype=np.float32), (rows, cols)),
        shape=(n_cells, adata_valid.n_obs),
    )
    X_agg = A @ adata_valid.X

    adata_cells = sc.AnnData(
        X=X_agg.tocsr() if sp.issparse(X_agg) else sp.csr_matrix(X_agg),
        var=adata_valid.var.copy(),
    )
    adata_cells.obs_names = [str(c) for c in unique_cells]
    del adata_crop, adata_valid, A
    gc.collect()

    sc.pp.normalize_total(adata_cells, target_sum=1e4)
    sc.pp.log1p(adata_cells)
    log.info(f"  AnnData: {adata_cells.n_obs:,} cells × {adata_cells.n_vars:,} genes")

    import celltypist
    log.info("  CellTypist 標注 (Human_Colorectal_Cancer.pkl)...")
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


# ─── Step 6: Spatial Matching ─────────────────────────────────────────────────

def match_gt_to_mcseg(
    gt: pd.DataFrame,
    mask: np.ndarray,
    cell_labels: pd.DataFrame,
) -> pd.DataFrame:
    """
    GT centroid → mask lookup → MCseg cell_id → predicted label.
    Returns gt DataFrame with added mcseg_cell_id and pred_label columns.
    """
    log.info("Step 6: 空間匹配 (GT centroid → MCseg mask lookup)")

    cell_id_to_label = cell_labels.set_index("cell_id")["broad_label"].to_dict()
    H, W = mask.shape

    cx = gt["cell_x"].values
    cy = gt["cell_y"].values
    in_bounds = (
        (cx >= CROP_X0) & (cx < CROP_X1) &
        (cy >= CROP_Y0) & (cy < CROP_Y1)
    )
    n_oor = (~in_bounds).sum()
    if n_oor > 0:
        log.warning(f"  {n_oor} GT cells outside crop bounds → marked unmatched")

    row_local = np.where(in_bounds, (cy - CROP_Y0).astype(np.int32).clip(0, H - 1), 0)
    col_local = np.where(in_bounds, (cx - CROP_X0).astype(np.int32).clip(0, W - 1), 0)

    cell_ids_for_gt = np.where(in_bounds, mask[row_local, col_local], 0)

    pred_labels = np.array([
        cell_id_to_label.get(int(cid), "unmatched") if cid > 0 else "unmatched"
        for cid in cell_ids_for_gt
    ])

    gt = gt.copy()
    gt["mcseg_cell_id"] = cell_ids_for_gt
    gt["pred_label"] = pred_labels

    n_matched = (pred_labels != "unmatched").sum()
    log.info(f"  GT cells matched: {n_matched:,} / {len(gt):,} ({n_matched/len(gt):.1%})")
    return gt


# ─── Step 7: F1 Calculation ───────────────────────────────────────────────────

def compute_f1(gt_matched: pd.DataFrame, results_dir: Path) -> dict:
    """計算 F1 並與 StarDist baseline 比較。"""
    from sklearn.metrics import (
        classification_report,
        confusion_matrix,
        f1_score,
        precision_recall_fscore_support,
    )
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    log.info("Step 7: 計算 F1 (vs StarDist baseline 0.708 / 0.758)")

    matched = gt_matched[gt_matched["pred_label"] != "unmatched"].copy()
    n_total   = len(gt_matched)
    n_matched = len(matched)
    log.info(f"  匹配率: {n_matched:,}/{n_total:,} ({n_matched/n_total:.1%})")

    valid = matched[matched["pred_label"] != "other"].copy()
    n_valid = len(valid)
    log.info(f"  有效標籤細胞: {n_valid:,} ({n_valid/n_total:.1%})")

    y_true = valid["gt_label"].values
    y_pred  = valid["pred_label"].values
    classes = ["epithelial cells", "stromal cells", "immune cells"]

    micro_f1    = f1_score(y_true, y_pred, labels=classes, average="micro",    zero_division=0)
    weighted_f1 = f1_score(y_true, y_pred, labels=classes, average="weighted", zero_division=0)
    macro_f1    = f1_score(y_true, y_pred, labels=classes, average="macro",    zero_division=0)

    log.info(f"  MCseg micro F1    = {micro_f1:.3f}  (StarDist: {STARDIST_MICRO_F1})")
    log.info(f"  MCseg weighted F1 = {weighted_f1:.3f}  (StarDist: {STARDIST_WEIGHTED_F1})")
    log.info(f"  MCseg macro F1    = {macro_f1:.3f}")
    log.info(f"\n{classification_report(y_true, y_pred, labels=classes, zero_division=0)}")

    prec, rec, f1_per, sup = precision_recall_fscore_support(
        y_true, y_pred, labels=classes, zero_division=0
    )
    pd.DataFrame({
        "class": classes,
        "precision": prec,
        "recall":    rec,
        "f1":        f1_per,
        "support":   sup,
    }).to_csv(str(results_dir / "f1_per_class.csv"), index=False)

    summary = {
        "n_gt_total":            n_total,
        "n_matched":             n_matched,
        "n_valid":               n_valid,
        "match_rate":            n_matched / n_total,
        "valid_rate":            n_valid / n_total,
        "micro_f1":              micro_f1,
        "weighted_f1":           weighted_f1,
        "macro_f1":              macro_f1,
        "stardist_micro_f1":     STARDIST_MICRO_F1,
        "stardist_weighted_f1":  STARDIST_WEIGHTED_F1,
        "delta_micro":           micro_f1    - STARDIST_MICRO_F1,
        "delta_weighted":        weighted_f1 - STARDIST_WEIGHTED_F1,
    }
    pd.DataFrame([summary]).to_csv(str(results_dir / "f1_summary.csv"), index=False)

    # Figure: bar chart + confusion matrix
    cm      = confusion_matrix(y_true, y_pred, labels=classes)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    short   = ["Epithelial", "Stromal", "Immune"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    x, w = np.arange(2), 0.35
    micro_v    = [STARDIST_MICRO_F1,    micro_f1]
    weighted_v = [STARDIST_WEIGHTED_F1, weighted_f1]
    bars1 = ax.bar(x - w/2, micro_v,    w, label="Micro F1",    color=["#4e79a7", "#f28e2b"])
    bars2 = ax.bar(x + w/2, weighted_v, w, label="Weighted F1", color=["#76b7b2", "#e15759"])
    ax.set_ylim(0, 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(["StarDist\n(ENACT baseline)", "MCseg v2\n(this work)"])
    ax.set_ylabel("F1 Score")
    ax.set_title("MCseg vs StarDist: Cell-type Annotation F1\n(ENACT CRC Expert GT, n=20,991)")
    ax.legend()
    for bar in bars1 + bars2:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.005, f"{h:.3f}",
                ha="center", va="bottom", fontsize=8)
    ax.axhline(0.708, ls="--", color="gray", lw=0.8, alpha=0.5)

    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Blues",
                xticklabels=short, yticklabels=short, ax=axes[1],
                cbar_kws={"label": "Proportion"})
    axes[1].set_xlabel("Predicted")
    axes[1].set_ylabel("True (GT)")
    axes[1].set_title(f"MCseg Confusion Matrix (Micro F1={micro_f1:.3f})")

    plt.tight_layout()
    fig.savefig(str(results_dir / "fig_f1_comparison.png"), dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  Figure saved: fig_f1_comparison.png")

    return summary


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_total = time.time()

    gt          = load_gt(EVAL_CSV)
    img         = crop_he_from_btf(BTF_PATH, RESULTS_DIR / "he_crop.tif")
    mask        = run_mcseg(img, RESULTS_DIR / "mcseg_mask.npy")
    del img
    gc.collect()
    attribution = run_bin_attribution(mask, TP_PATH)
    cell_labels = build_anndata_and_annotate(attribution, H5_PATH,
                                             RESULTS_DIR / "celltypist_labels.csv")
    gt_matched  = match_gt_to_mcseg(gt, mask, cell_labels)
    gt_matched.to_csv(str(RESULTS_DIR / "gt_matched.csv"), index=False)
    summary     = compute_f1(gt_matched, RESULTS_DIR)

    elapsed = time.time() - t_total
    log.info(f"\n{'='*60}")
    log.info(f"完成 ({elapsed/60:.1f} min)  結果: {RESULTS_DIR}")
    log.info(f"MCseg  micro F1    = {summary['micro_f1']:.3f}  (StarDist: {STARDIST_MICRO_F1})")
    log.info(f"MCseg  weighted F1 = {summary['weighted_f1']:.3f}  (StarDist: {STARDIST_WEIGHTED_F1})")
    log.info(f"Δ micro = {summary['delta_micro']:+.3f}  Δ weighted = {summary['delta_weighted']:+.3f}")
    if summary["micro_f1"] >= STARDIST_MICRO_F1:
        log.info("→ MCseg ≥ StarDist: 可宣稱相當或優越的細胞型態分類性能")
    else:
        log.info("→ MCseg < StarDist: 重點轉向 MCseg-wba 組合效果")
    log.info("="*60)


if __name__ == "__main__":
    main()
