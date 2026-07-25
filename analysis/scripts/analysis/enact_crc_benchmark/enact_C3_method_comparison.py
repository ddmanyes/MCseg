"""
enact_C3_method_comparison.py
Method C: 5-method cell-type annotation F1 comparison on ENACT CRC region

Methods compared (all evaluated against 20,991 GT cells):
  1. StarDist + WBA  — ENACT's own CellTypist result (directly reused)
  2. MCseg   + WBA  — mcseg_mask_7pass.npy
  3. SR      + WBA  — Space Ranger cell_segmentations.geojson
  4. NUC     + WBA  — cellpose_nuc_mask.npy
  5. ProSeg  + WBA  — proseg_mask.npy

Annotation: CellTypist (Human_Colorectal_Cancer.pkl)
  — consistent with ENACT's default pipeline
  — GT labels from ENACT celltypist eval (same source)

Usage:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/enact_C3_method_comparison.py
"""

from __future__ import annotations
import json
import logging
import time
from pathlib import Path

import anndata
import celltypist
import cv2
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import f1_score, classification_report

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A       = Path("/Volumes/SSD/plan_a")
ENACT_CRC    = PLAN_A / "tissue sample" / "ENACT_supporting_files" / "public_data" / "human_colorectal"
H5_PATH      = ENACT_CRC / "input_files" / "filtered_feature_bc_matrix.h5"
TP_PATH      = ENACT_CRC / "input_files" / "tissue_positions.parquet"
CELLTYPIST_EVAL = (
    ENACT_CRC / "paper_results" / "chunks" / "weighted_by_area"
    / "celltypist_results" / "eval" / "cell_annotation_eval.csv"
)
CELLS_DF     = ENACT_CRC / "paper_results" / "cells_df.csv"

SR_CELL_GJ   = PLAN_A / "tissue sample" / "CRC" / "visium" / "official_v4" / "segmented_outputs" / "segmented_outputs" / "cell_segmentations.geojson"

RESULT_DIR   = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
ENACT_F1     = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
MCSEG_MASK   = ENACT_F1 / "mcseg_mask_7pass.npy"
MCSEG_WBA    = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1_wba" / "wba_attribution.csv"

CELLTYPIST_MODEL = "Human_Colorectal_Cancer.pkl"

# ─── ENACT crop constants ──────────────────────────────────────────────────────
COL_OFFSET        = 40598
CROP_X0, CROP_X1  = 5154, 15242
CROP_Y0, CROP_Y1  = 4635, 18599
CROP_H = CROP_Y1 - CROP_Y0   # 13964
CROP_W = CROP_X1 - CROP_X0   # 10088

BTF_COL_MIN = CROP_X0 + COL_OFFSET
BTF_COL_MAX = CROP_X1 + COL_OFFSET
BTF_ROW_MIN = CROP_Y0
BTF_ROW_MAX = CROP_Y1

# CellTypist fine label → broad class (derived from celltypist_labels.csv)
CELLTYPIST_BROAD: dict[str, str] = {
    "CMS1":                      "epithelial cells",
    "CMS2":                      "epithelial cells",
    "CMS3":                      "epithelial cells",
    "Goblet cells":              "epithelial cells",
    "Intermediate":              "epithelial cells",
    "Proliferating":             "epithelial cells",
    "Stem-like/TA":              "epithelial cells",
    "Myofibroblasts":            "stromal cells",
    "Lymphatic ECs":             "stromal cells",
    "Pericytes":                 "stromal cells",
    "Smooth muscle cells":       "stromal cells",
    "Stalk-like ECs":            "stromal cells",
    "Enteric glial cells":       "stromal cells",
    "Stromal 1":                 "stromal cells",
    "Stromal 2":                 "stromal cells",
    "Stromal 3":                 "stromal cells",
    "Proliferative ECs":         "stromal cells",
    "CD4+ T cells":              "immune cells",
    "CD8+ T cells":              "immune cells",
    "CD19+CD20+ B":              "immune cells",
    "IgA+ Plasma":               "immune cells",
    "IgG+ Plasma":               "immune cells",
    "NK cells":                  "immune cells",
    "Pro-inflammatory":          "immune cells",
    "Regulatory T cells":        "immune cells",
    "SPP1+":                     "immune cells",
    "Mast cells":                "immune cells",
    "T follicular helper cells": "other",
    "T helper 17 cells":         "other",
    "Tip-like ECs":              "other",
    "Unknown":                   "other",
    "cDC":                       "other",
    "gamma delta T cells":       "other",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


# ─── Shared data loaders ──────────────────────────────────────────────────────

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
    """Load GT cells from ENACT CellTypist eval (20,991 cells, same as ENACT paper)."""
    df = pd.read_csv(str(CELLTYPIST_EVAL), usecols=["cell_x", "cell_y", "gt_label"])
    df["cell_x"] = df["cell_x"] - CROP_X0
    df["cell_y"] = df["cell_y"] - CROP_Y0
    log.info(f"  GT cells: {len(df):,}  gt_label: {df['gt_label'].value_counts().to_dict()}")
    return df


# ─── WBA: bin→cell attribution ────────────────────────────────────────────────

def compute_wba_from_mask(mask: np.ndarray, tp_crop: pd.DataFrame) -> pd.DataFrame:
    row_c = tp_crop["row_c"].values.clip(0, CROP_H - 1)
    col_c = tp_crop["col_c"].values.clip(0, CROP_W - 1)
    cell_ids = mask[row_c, col_c]
    df = pd.DataFrame({"barcode": tp_crop["barcode"].values, "cell_id": cell_ids})
    df = df[df["cell_id"] > 0]
    n_bins = df.groupby("cell_id").size().rename("n_bins")
    df = df.join(n_bins, on="cell_id")
    df["weight"] = 1.0 / df["n_bins"]
    return df[["barcode", "cell_id", "weight"]]


# ─── CellTypist annotation ────────────────────────────────────────────────────

def aggregate_expr_to_adata(wba: pd.DataFrame, adata: sc.AnnData) -> anndata.AnnData:
    """Build cell × gene AnnData from WBA-weighted bin expression (sparse W @ X)."""
    import scipy.sparse as sp

    barcode_to_idx = {b: i for i, b in enumerate(adata.obs_names)}

    src_idxs = np.array([barcode_to_idx.get(b, -1) for b in wba["barcode"].values])
    valid    = src_idxs >= 0
    src_v    = src_idxs[valid]
    cell_v   = wba["cell_id"].values[valid]
    weight_v = wba["weight"].values[valid].astype(np.float32)

    unique_cells = np.unique(cell_v)
    cell_to_row  = {c: i for i, c in enumerate(unique_cells)}
    row_v        = np.array([cell_to_row[c] for c in cell_v], dtype=np.int32)

    # W: (n_cells × n_bins) sparse weight matrix
    n_cells = len(unique_cells)
    n_bins  = adata.shape[0]
    W = sp.csr_matrix(
        (weight_v, (row_v, src_v)),
        shape=(n_cells, n_bins),
        dtype=np.float32,
    )

    X_src = adata.X
    if not sp.issparse(X_src):
        X_src = sp.csr_matrix(X_src)
    X_cell = W @ X_src   # (n_cells × n_genes)

    cell_adata = anndata.AnnData(
        X=X_cell.astype(np.float32),
        obs=pd.DataFrame(index=[str(c) for c in unique_cells]),
        var=adata.var.copy(),
    )
    log.info(f"  Cell AnnData: {cell_adata.shape}")
    return cell_adata


def celltypist_annotate(cell_adata: anndata.AnnData) -> pd.DataFrame:
    """Run CellTypist on cell AnnData. Returns DataFrame with cell_id, pred_label."""
    sc.pp.normalize_total(cell_adata, target_sum=1e4)
    sc.pp.log1p(cell_adata)

    predictions = celltypist.annotate(
        cell_adata,
        model=CELLTYPIST_MODEL,
        majority_voting=False,
    )
    pred_df = predictions.predicted_labels.copy()
    pred_df.index = cell_adata.obs_names
    pred_df["broad_label"] = pred_df["predicted_labels"].map(CELLTYPIST_BROAD).fillna("other")

    cell_id_arr = pred_df.index.astype(int).values
    result = pd.DataFrame({
        "cell_id":    cell_id_arr,
        "pred_label": pred_df["broad_label"].values,
    })
    n_mapped = (pred_df["broad_label"] != "other").sum()
    log.info(f"  CellTypist: {len(result):,} cells  ({n_mapped:,} mapped to epi/stromal/immune)")
    return result


# ─── GT matching: mask lookup ─────────────────────────────────────────────────

def match_gt_by_mask_lookup(gt: pd.DataFrame, mask: np.ndarray,
                             annotations: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    H, W = mask.shape
    ann_map = annotations.set_index("cell_id")["pred_label"].to_dict()

    gt = gt.dropna(subset=["cell_x", "cell_y"]).reset_index(drop=True)
    row = gt["cell_y"].round().astype(int).clip(0, H - 1).values
    col = gt["cell_x"].round().astype(int).clip(0, W - 1).values
    cell_ids = mask[row, col]

    pred = [
        ann_map.get(int(cid), "unmatched") if cid > 0 else "unmatched"
        for cid in cell_ids
    ]
    result = gt.copy()
    result["pred_label"] = pred
    n_matched = sum(1 for cid in cell_ids if cid > 0)
    log.info(f"  Match rate: {n_matched}/{len(gt)} = {100*n_matched/len(gt):.1f}%")
    return result, n_matched


# ─── F1 computation ───────────────────────────────────────────────────────────

def compute_f1(gt_labels: pd.Series, pred_labels: pd.Series, name: str) -> dict:
    classes = ["epithelial cells", "stromal cells", "immune cells"]
    report = classification_report(gt_labels, pred_labels,
                                    labels=classes, output_dict=True, zero_division=0)
    wf1 = f1_score(gt_labels, pred_labels, labels=classes,
                    average="weighted", zero_division=0)
    log.info(f"  Weighted F1 [{name}]: {wf1:.3f}")
    for cls in classes:
        m = report.get(cls, {})
        log.info(f"    {cls.split()[0]:12s}: P={m.get('precision',0):.3f}  R={m.get('recall',0):.3f}  F1={m.get('f1-score',0):.3f}")
    return {"method": name, "weighted_f1": wf1, **{
        f"{cls.split()[0]}_{k}": report.get(cls, {}).get(metric, 0)
        for cls in classes
        for k, metric in [("p", "precision"), ("r", "recall"), ("f1", "f1-score")]
    }}


# ─── Method evaluators ────────────────────────────────────────────────────────

def eval_stardist_celltypist() -> dict:
    """Load StarDist+WBA+CellTypist results directly from ENACT (reference)."""
    log.info("=== StarDist (ENACT CellTypist — reference) ===")
    df = pd.read_csv(str(CELLTYPIST_EVAL))
    df["pred"] = df["pred_label_clean"].replace({"no label": "unmatched"})
    record = compute_f1(df["gt_label"], df["pred"], "StarDist+WBA (ref)")
    record["match_rate"] = (df["pred"] != "unmatched").mean()
    matched_only = df[df["pred"] != "unmatched"]
    record["subset_wf1"] = f1_score(
        matched_only["gt_label"], matched_only["pred"],
        labels=["epithelial cells", "stromal cells", "immune cells"],
        average="weighted", zero_division=0,
    )
    return record


def cells_df_to_mask() -> np.ndarray:
    from shapely import wkt as shapely_wkt
    log.info("  Loading StarDist cells_df.csv...")
    df = pd.read_csv(str(CELLS_DF), usecols=["geometry", "cell_x", "cell_y"])
    in_crop = (
        (df["cell_x"] >= CROP_X0) & (df["cell_x"] < CROP_X1) &
        (df["cell_y"] >= CROP_Y0) & (df["cell_y"] < CROP_Y1)
    )
    df = df[in_crop].reset_index(drop=True)
    log.info(f"  StarDist cells in crop: {len(df):,}")

    mask = np.zeros((CROP_H, CROP_W), dtype=np.int32)
    for cell_id, row in enumerate(df.itertuples(), start=1):
        try:
            geom = shapely_wkt.loads(row.geometry)
        except Exception:
            continue
        if geom.geom_type == "Polygon":
            polys = [geom.exterior.coords]
        elif geom.geom_type == "MultiPolygon":
            polys = [g.exterior.coords for g in geom.geoms]
        else:
            continue
        for coords in polys:
            pts = np.array(coords, dtype=np.float64)
            local = np.stack([
                pts[:, 0] - CROP_X0,
                pts[:, 1] - CROP_Y0,
            ], axis=1).astype(np.int32).reshape((-1, 1, 2))
            cv2.fillPoly(mask, [local], color=cell_id)

    log.info(f"  StarDist mask: {mask.max():,} cells")
    return mask


def eval_method_with_mask(mask: np.ndarray, wba: pd.DataFrame,
                           adata: sc.AnnData, gt: pd.DataFrame,
                           method_name: str) -> dict:
    """mask → WBA → CellTypist → mask lookup GT match → F1."""
    log.info(f"=== {method_name} ===")
    t0 = time.time()
    cell_adata  = aggregate_expr_to_adata(wba, adata)
    annotations = celltypist_annotate(cell_adata)
    matched, n_matched = match_gt_by_mask_lookup(gt, mask, annotations)
    match_rate = n_matched / len(gt)

    record = compute_f1(matched["gt_label"], matched["pred_label"], method_name)
    record["match_rate"] = match_rate

    matched_only = matched[matched["pred_label"] != "unmatched"]
    classes = ["epithelial cells", "stromal cells", "immune cells"]
    subset_wf1 = f1_score(
        matched_only["gt_label"], matched_only["pred_label"],
        labels=classes, average="weighted", zero_division=0,
    )
    record["subset_wf1"] = subset_wf1
    log.info(f"  Subset wF1 (matched only, n={len(matched_only):,}): {subset_wf1:.3f}")
    log.info(f"  Done in {time.time()-t0:.0f}s")
    return record


# ─── Figure ───────────────────────────────────────────────────────────────────

def plot_results(records: list[dict], out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    df = pd.DataFrame(records)
    df.to_csv(str(out_dir / "method_comparison_f1.csv"), index=False)

    methods = df["method"].tolist()
    x = np.arange(len(methods))
    palette = ["#4A90D9", "#E07B54", "#7DC07D", "#9B59B6", "#F39C12", "#1ABC9C"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=150)

    ax = axes[0]
    ax.bar(x, df["match_rate"].values * 100,
           color=palette[:len(methods)], alpha=0.9, edgecolor="white", width=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Coverage rate (%)")
    ax.set_title("GT cell coverage rate\n(% of GT centroids matched to a predicted cell)")
    ax.set_ylim(0, 110)
    for xi, v in enumerate(df["match_rate"].values):
        ax.text(xi, v * 100 + 1.0, f"{v*100:.1f}%", ha="center", va="bottom",
                fontsize=9, fontweight="bold")

    ax2 = axes[1]
    ax2.bar(x, df["subset_wf1"].values, color=palette[:len(methods)],
            alpha=0.9, edgecolor="white", width=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(methods, rotation=20, ha="right", fontsize=9)
    ax2.set_ylabel("Weighted F1 (matched cells only)")
    ax2.set_title("Annotation accuracy — subset wF1\n(evaluated only on matched GT cells)")
    ax2.set_ylim(0, 1.0)
    for xi, v in enumerate(df["subset_wf1"].values):
        ax2.text(xi, v + 0.01, f"{v:.3f}", ha="center", va="bottom",
                 fontsize=9, fontweight="bold")

    plt.suptitle(
        "ENACT CRC Region — Cell-Type Annotation Benchmark\n"
        "(N=20,991 GT cells, CellTypist Human_Colorectal_Cancer.pkl, WBA attribution)",
        fontsize=11,
    )
    plt.tight_layout()
    out = out_dir / "fig_method_comparison_f1.png"
    plt.savefig(str(out), dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"  Figure saved: {out.name}")


def _geojson_to_mask(data: dict, label: str) -> np.ndarray:
    mask = np.zeros((CROP_H, CROP_W), dtype=np.int32)
    cell_id = 1
    n_skipped = 0
    for feat in data["features"]:
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            rings_list = [geom["coordinates"]]
        elif geom["type"] == "MultiPolygon":
            rings_list = geom["coordinates"]
        else:
            continue
        first_ring = np.array(rings_list[0][0])
        cx, cy = first_ring[:, 0].mean(), first_ring[:, 1].mean()
        if not (BTF_COL_MIN <= cx < BTF_COL_MAX and BTF_ROW_MIN <= cy < BTF_ROW_MAX):
            n_skipped += 1
            continue
        for rings in rings_list:
            for ring in rings:
                pts = np.array(ring, dtype=np.float64)
                local = np.stack([
                    pts[:, 0] - BTF_COL_MIN,
                    pts[:, 1] - CROP_Y0,
                ], axis=1).astype(np.int32).reshape((-1, 1, 2))
                cv2.fillPoly(mask, [local], color=cell_id)
        cell_id += 1
    log.info(f"  {label}: {cell_id-1:,} cells in crop  ({n_skipped:,} outside)")
    return mask


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading shared resources...")
    tp_crop = load_tissue_positions()
    adata   = load_adata(tp_crop)
    gt      = load_gt_cells()

    results = []

    # 1. StarDist — ENACT CellTypist (reference, directly reused)
    results.append(eval_stardist_celltypist())

    # 2. StarDist — same CellTypist scorer (fair comparison baseline)
    sd_mask_path = RESULT_DIR / "stardist_mask.npy"
    if sd_mask_path.exists():
        log.info("Loading cached StarDist mask...")
        sd_mask = np.load(str(sd_mask_path))
    else:
        sd_mask = cells_df_to_mask()
        np.save(str(sd_mask_path), sd_mask)
    sd_wba = compute_wba_from_mask(sd_mask, tp_crop)
    results.append(eval_method_with_mask(sd_mask, sd_wba, adata, gt, "StarDist+WBA"))
    del sd_mask

    # 3. MCseg (7-pass)
    log.info("Loading MCseg mask + WBA...")
    mcseg_mask = np.load(str(MCSEG_MASK))
    mcseg_wba  = pd.read_csv(str(MCSEG_WBA))
    results.append(eval_method_with_mask(mcseg_mask, mcseg_wba, adata, gt, "MCseg+WBA"))
    del mcseg_mask

    # 4. SR
    sr_mask_path = RESULT_DIR / "sr_mask.npy"
    if sr_mask_path.exists():
        sr_mask = np.load(str(sr_mask_path))
    else:
        with open(str(SR_CELL_GJ)) as f:
            data = json.load(f)
        sr_mask = _geojson_to_mask(data, "SR")
        np.save(str(sr_mask_path), sr_mask)
    sr_wba = compute_wba_from_mask(sr_mask, tp_crop)
    results.append(eval_method_with_mask(sr_mask, sr_wba, adata, gt, "SR+WBA"))
    del sr_mask

    # 5. NUC (Cellpose nuclei)
    nuc_mask = np.load(str(RESULT_DIR / "cellpose_nuc_mask.npy"))
    nuc_wba  = compute_wba_from_mask(nuc_mask, tp_crop)
    results.append(eval_method_with_mask(nuc_mask, nuc_wba, adata, gt, "NUC+WBA"))
    del nuc_mask

    # 6. ProSeg
    proseg_mask_path = RESULT_DIR / "proseg_mask.npy"
    if proseg_mask_path.exists():
        log.info("Loading ProSeg mask...")
        proseg_mask = np.load(str(proseg_mask_path))
        proseg_wba  = compute_wba_from_mask(proseg_mask, tp_crop)
        results.append(eval_method_with_mask(proseg_mask, proseg_wba, adata, gt, "ProSeg+WBA"))
        del proseg_mask
    else:
        log.info("ProSeg mask not found — skipping")

    log.info("\n" + "=" * 55)
    log.info("SUMMARY")
    log.info("=" * 55)
    for r in results:
        log.info(f"  {r['method']:22s}  match={r['match_rate']*100:.1f}%  subset_wF1={r['subset_wf1']:.3f}")

    plot_results(results, RESULT_DIR)


if __name__ == "__main__":
    main()
