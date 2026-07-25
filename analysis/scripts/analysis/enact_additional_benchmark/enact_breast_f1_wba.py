"""
enact_breast_f1_wba.py
======================
ENACT Breast Cancer F1 Benchmark — MCseg + WBA vs StarDist baseline

Dataset: Visium HD FF Human Breast Cancer (10x Genomics, Space Ranger 3.1.1)
GT:      GeoJSON region annotations (Tumor / Stroma) from ENACT paper
         → 2-class: epithelial vs stromal

Pipeline (mirrors enact_C2_crc_f1_wba.py):
  Step 1  Load GT: parse GeoJSON polygons + assign region label to each StarDist cell
  Step 2  Compute StarDist baseline F1 from ENACT merged_results.csv
  Step 3  Crop H&E from tissue_image.tif at GT bounding box
  Step 4  MCseg segmentation (deployment-mode, voronoi_d=8, tiled)
  Step 5  WBA attribution (7×7 footprint, tissue_positions.parquet)
  Step 6  AnnData + CellTypist annotation (Cells_Adult_Breast.pkl)
  Step 7  Spatial matching: GT region → MCseg centroid → compare predicted label
  Step 8  Compute F1 + comparison figure

Coordinate system: Space Ranger 3.1.1 full-res (25840 × 34009 px, 0.2738 µm/px)
GT region:  x0=3588, y0=1641, x1=27772, y1=23646  (from GeoJSON filename)
"""

from __future__ import annotations

import gc
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tifffile

# ─── Paths ────────────────────────────────────────────────────────────────────

PLAN_A   = Path("/Volumes/SSD/plan_a")
BSC_311  = PLAN_A / "tissue sample" / "bsc" / "binned_outputs_311" / "binned_outputs" / "square_002um"
TIF_PATH = PLAN_A / "tissue sample" / "Visium_HD_Human_Breast_Cancer_Fresh_Frozen_tissue_image.tif"
H5_PATH  = BSC_311 / "filtered_feature_bc_matrix.h5"
TP_PATH  = BSC_311 / "spatial" / "tissue_positions.parquet"
GJ_PATH  = (PLAN_A / "tissue sample" / "ENACT_additional_samples" / "human_breast" /
            "annotations" / "Visium_HD_Human_Breast_Cancer_Fresh_Frozen-wsi-3588_1641_27772_23646.geojson")
STARDIST_CSV = (PLAN_A / "tissue sample" / "ENACT_additional_samples" /
                "human_breast" / "merged_results.csv")

MSSEG_ROOT  = PLAN_A / "MSseg"
RESULTS_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_breast_f1"

# GT bounding box (WSI pixel coords, 3.1.1 space)
CROP_X0, CROP_X1 = 3588, 27772   # pxl_col
CROP_Y0, CROP_Y1 = 1641, 23646   # pxl_row

# WBA footprint: bin spacing ≈ 7.30 px → BIN_HALF = 3 → 7×7 window
BIN_HALF = 3

# CellTypist label → broad class (Cells_Adult_Breast.pkl actual labels)
LABEL_MAP: dict[str, str] = {
    # epithelial — luminal & basal
    "LummHR-SCGB":    "epithelial",
    "LummHR-active":  "epithelial",
    "LummHR-major":   "epithelial",
    "Lumsec-HLA":     "epithelial",
    "Lumsec-KIT":     "epithelial",
    "Lumsec-basal":   "epithelial",
    "Lumsec-lac":     "epithelial",
    "Lumsec-major":   "epithelial",
    "Lumsec-myo":     "epithelial",
    "Lumsec-prol":    "epithelial",
    "basal":          "epithelial",
    # stromal — fibroblast, endothelial, pericyte, SMC
    "Fibro-SFRP4":     "stromal",
    "Fibro-major":     "stromal",
    "Fibro-matrix":    "stromal",
    "Fibro-prematrix": "stromal",
    "Lymph-immune":    "stromal",
    "Lymph-major":     "stromal",
    "Lymph-valve1":    "stromal",
    "Lymph-valve2":    "stromal",
    "Vas-arterial":    "stromal",
    "Vas-capillary":   "stromal",
    "Vas-venous":      "stromal",
    "pericytes":       "stromal",
    "vsmc":            "stromal",
    # immune
    "CD4-Th":          "immune",
    "CD4-Th-like":     "immune",
    "CD4-Treg":        "immune",
    "CD4-activated":   "immune",
    "CD4-naive":       "immune",
    "CD8-Tem":         "immune",
    "CD8-Trm":         "immune",
    "CD8-activated":   "immune",
    "GD":              "immune",
    "Macro-IFN":       "immune",
    "Macro-lipo":      "immune",
    "Macro-m1":        "immune",
    "Macro-m1-CCL":    "immune",
    "Macro-m2":        "immune",
    "Macro-m2-CXCL":   "immune",
    "Mast":            "immune",
    "Mono-classical":      "immune",
    "Mono-non-classical":  "immune",
    "NK":              "immune",
    "NK-ILCs":         "immune",
    "NKT":             "immune",
    "Neutrophil":      "immune",
    "T_prol":          "immune",
    "b_naive":         "immune",
    "bmem_switched":   "immune",
    "bmem_unswitched": "immune",
    "cDC1":            "immune",
    "cDC2":            "immune",
    "mDC":             "immune",
    "mye-prol":        "immune",
    "pDC":             "immune",
    "plasma_IgA":      "immune",
    "plasma_IgG":      "immune",
}

# ENACT merged_results.csv cell_type → broad class (breast-specific ENACT labels)
STARDIST_LABEL_MAP: dict[str, str] = {
    "Cancer stem cell":           "epithelial",
    "Luminal progenitor":         "epithelial",
    "Basal epithelial cell":      "epithelial",
    "Basal cell":                 "epithelial",
    "Procr+ basal cell":          "epithelial",
    "Luminal epithelial cell":    "epithelial",
    "Hormone sensing progenitor": "epithelial",
    "Luminal progenitor cell":    "epithelial",
    "Myoepithelial cell":         "stromal",
    "Immune cell":                "immune",
}

GT_REGION_MAP: dict[str, str] = {
    "Tumor":  "epithelial",
    "Stroma": "stromal",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


# ─── GT helpers ───────────────────────────────────────────────────────────────

def load_gt_regions() -> list[tuple[object, str]]:
    try:
        from shapely.geometry import shape
    except ImportError:
        import subprocess; subprocess.run(["uv", "add", "shapely"], check=True)
        from shapely.geometry import shape

    log.info("Step 1: 載入 GeoJSON GT 區域")
    with open(str(GJ_PATH)) as f:
        gj = json.load(f)

    regions = []
    for feat in gj["features"]:
        name  = feat["properties"]["classification"]["name"]
        label = GT_REGION_MAP.get(name)
        if label is None:
            continue
        regions.append((shape(feat["geometry"]), label))
        log.info(f"  {name} → {label}")
    return regions


def assign_gt_labels(df: pd.DataFrame, regions: list,
                     x_col: str = "cell_x", y_col: str = "cell_y") -> pd.Series:
    from shapely.geometry import Point
    from shapely.strtree import STRtree

    polys  = [r[0] for r in regions]
    lmap   = {id(p): r[1] for r, p in zip(regions, polys)}
    tree   = STRtree(polys)

    xs, ys = df[x_col].values, df[y_col].values
    labels = ["unmatched"] * len(df)
    for i, (x, y) in enumerate(zip(xs, ys)):
        pt = Point(x, y)
        for idx in tree.query(pt):
            if polys[idx].contains(pt):
                labels[i] = lmap[id(polys[idx])]
                break
    return pd.Series(labels, index=df.index)


# ─── Step 2: StarDist baseline ────────────────────────────────────────────────

def compute_stardist_f1(regions: list) -> dict:
    log.info("Step 2: StarDist baseline F1 from merged_results.csv")
    sd = pd.read_csv(str(STARDIST_CSV))
    log.info(f"  StarDist cells: {len(sd):,}")

    sd["gt_label"] = assign_gt_labels(sd, regions)
    in_region = sd[sd["gt_label"] != "unmatched"].copy()
    log.info(f"  cells in GT regions: {len(in_region):,} ({len(in_region)/len(sd):.1%})")

    in_region["broad_label"] = in_region["cell_type"].map(STARDIST_LABEL_MAP).fillna("other")
    valid = in_region[in_region["broad_label"] != "other"]
    log.info(f"  cells with valid broad label: {len(valid):,}")

    from sklearn.metrics import f1_score, classification_report
    classes = sorted({c for c in ["epithelial", "stromal"] if c in valid["gt_label"].values})
    if not valid.empty and classes:
        micro    = f1_score(valid["gt_label"], valid["broad_label"], labels=classes,
                            average="micro", zero_division=0)
        weighted = f1_score(valid["gt_label"], valid["broad_label"], labels=classes,
                            average="weighted", zero_division=0)
        log.info(f"  StarDist micro F1={micro:.3f}  weighted={weighted:.3f}")
        log.info(f"\n{classification_report(valid['gt_label'], valid['broad_label'], labels=classes, zero_division=0)}")
    else:
        micro = weighted = 0.0
        log.warning("  insufficient valid cells")

    summary = {"micro_f1": micro, "weighted_f1": weighted, "n_valid": len(valid)}
    pd.DataFrame([summary]).to_csv(str(RESULTS_DIR / "stardist_f1_breast.csv"), index=False)
    return summary


# ─── Step 3: Crop H&E ─────────────────────────────────────────────────────────

def crop_he(crop_tif: Path) -> np.ndarray:
    if crop_tif.exists():
        log.info(f"Step 3: 載入 H&E crop: {crop_tif.name}")
        img = tifffile.imread(str(crop_tif))
        return img[..., :3] if img.ndim == 3 and img.shape[-1] == 4 else img

    log.info(f"Step 3: 從 TIF 讀取 crop ({CROP_Y1-CROP_Y0} × {CROP_X1-CROP_X0} px)")
    t0 = time.time()
    img = tifffile.imread(str(TIF_PATH))[CROP_Y0:CROP_Y1, CROP_X0:CROP_X1]
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  shape: {img.shape}  ({time.time()-t0:.0f}s)")
    tifffile.imwrite(str(crop_tif), img, compression="zlib")
    return img


# ─── Step 4: MCseg segmentation ───────────────────────────────────────────────

def run_mcseg(img: np.ndarray, mask_npy: Path) -> np.ndarray:
    if mask_npy.exists():
        log.info(f"Step 4: 載入 MCseg 遮罩: {mask_npy.name}")
        mask = np.load(str(mask_npy))
        log.info(f"  shape: {mask.shape}  cells: {int(mask.max()):,}")
        return mask

    sys.path.insert(0, str(MSSEG_ROOT / "backend"))
    from src.segmentation.cellpose_runner import run_tiled_mcseg_v2

    cfg = {
        "use_gpu": True, "batch_size": 2,
        "dia_small": 13.0, "dia_mid": 17.0, "dia_large": 22.0,
        "use_hematoxylin": True, "use_cpsam": False,
        "voronoi_distance": 8, "flow_threshold": 0.4,
        "cellprob_threshold": -2.0, "min_size": 20, "max_size": 6000,
        "clahe_clip_limit": 3.0, "use_transcript_rescue": False,
    }
    log.info("Step 4: MCseg v2 tiled (voronoi_d=8)")

    mask = run_tiled_mcseg_v2(
        img, cfg, tile_size=1024, overlap=128,
        progress_callback=lambda p, m: log.info(f"  [{p*100:.0f}%] {m}"),
    )
    np.save(str(mask_npy), mask)
    log.info(f"  儲存: {mask_npy.name}  cells: {int(mask.max()):,}")
    return mask


# ─── Step 5: WBA attribution ──────────────────────────────────────────────────

def run_wba(mask: np.ndarray, wba_csv: Path) -> pd.DataFrame:
    if wba_csv.exists():
        log.info(f"Step 5: 載入 WBA: {wba_csv.name}")
        return pd.read_csv(wba_csv)

    log.info("Step 5: WBA attribution (7×7 footprint)")
    t0 = time.time()
    tp = pd.read_parquet(str(TP_PATH),
                         columns=["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    tp = tp[tp["in_tissue"] == 1]

    in_crop = (
        (tp["pxl_col_in_fullres"] >= CROP_X0) & (tp["pxl_col_in_fullres"] < CROP_X1) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0) & (tp["pxl_row_in_fullres"] < CROP_Y1)
    )
    tp_crop = tp[in_crop].copy().reset_index(drop=True)
    log.info(f"  bins in crop: {len(tp_crop):,}")

    row_c = (tp_crop["pxl_row_in_fullres"].values - CROP_Y0).astype(np.int32)
    col_c = (tp_crop["pxl_col_in_fullres"].values - CROP_X0).astype(np.int32)
    H, W  = mask.shape

    barcodes_out: list[str] = []
    cell_ids_out: list[int] = []
    weights_out:  list[float] = []
    CHUNK = 200_000

    for start in range(0, len(tp_crop), CHUNK):
        end      = min(start + CHUNK, len(tp_crop))
        barcodes = tp_crop["barcode"].values[start:end]
        rows_b   = row_c[start:end]
        cols_b   = col_c[start:end]
        for i in range(end - start):
            r, c = rows_b[i], cols_b[i]
            window = mask[max(0, r-BIN_HALF):min(H, r+BIN_HALF+1),
                         max(0, c-BIN_HALF):min(W, c+BIN_HALF+1)]
            flat   = window.ravel()
            ids, counts = np.unique(flat[flat > 0], return_counts=True)
            if not len(ids):
                continue
            total = float(counts.sum())
            for cid, cnt in zip(ids, counts):
                barcodes_out.append(barcodes[i])
                cell_ids_out.append(int(cid))
                weights_out.append(cnt / total)
        if (start // CHUNK) % 5 == 0:
            log.info(f"  {end:,}/{len(tp_crop):,} bins ({end/len(tp_crop):.0%})")

    log.info(f"  WBA pairs: {len(barcodes_out):,}  {time.time()-t0:.0f}s")
    wba = pd.DataFrame({"barcode": barcodes_out, "cell_id": cell_ids_out, "weight": weights_out})
    wba.to_csv(str(wba_csv), index=False)
    return wba


# ─── Step 6: AnnData + CellTypist ─────────────────────────────────────────────

def build_anndata_wba(wba: pd.DataFrame, ct_csv: Path) -> pd.DataFrame:
    if ct_csv.exists():
        log.info(f"Step 6: 載入 CellTypist: {ct_csv.name}")
        return pd.read_csv(ct_csv)

    log.info("Step 6: AnnData + CellTypist (Cells_Adult_Breast.pkl)")
    import scipy.sparse as sp

    adata_full = sc.read_10x_h5(str(H5_PATH))
    adata_full.var_names_make_unique()
    mask_obs  = adata_full.obs_names.isin(wba["barcode"].unique())
    adata_sub = adata_full[mask_obs].copy()
    del adata_full; gc.collect()
    log.info(f"  barcodes matched: {adata_sub.n_obs:,}")

    bc_to_idx    = {bc: i for i, bc in enumerate(adata_sub.obs_names)}
    unique_cells = np.sort(wba["cell_id"].unique()).astype(np.int32)
    cell_to_idx  = {int(c): i for i, c in enumerate(unique_cells)}
    wba_v = wba[wba["barcode"].isin(bc_to_idx)].copy()

    W = sp.csr_matrix(
        (wba_v["weight"].values.astype(np.float32),
         (wba_v["cell_id"].map(cell_to_idx).values.astype(np.int32),
          wba_v["barcode"].map(bc_to_idx).values.astype(np.int32))),
        shape=(len(unique_cells), adata_sub.n_obs),
    )
    X = adata_sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X_agg = W @ X

    adata_cells = sc.AnnData(
        X=X_agg.tocsr() if sp.issparse(X_agg) else sp.csr_matrix(X_agg),
        var=adata_sub.var.copy(),
    )
    adata_cells.obs_names = [str(c) for c in unique_cells]
    del adata_sub, W, X, X_agg; gc.collect()
    sc.pp.normalize_total(adata_cells, target_sum=1e4)
    sc.pp.log1p(adata_cells)

    import celltypist
    preds = celltypist.annotate(adata_cells, model="Cells_Adult_Breast.pkl", majority_voting=False)
    ct_labels = preds.predicted_labels["predicted_labels"].values

    df = pd.DataFrame({
        "cell_id":          unique_cells,
        "celltypist_label": ct_labels,
        "broad_label":      [LABEL_MAP.get(lbl, "other") for lbl in ct_labels],
    })
    df.to_csv(str(ct_csv), index=False)
    log.info(f"  broad_label:\n{df['broad_label'].value_counts().to_string()}")
    return df


# ─── Step 7: spatial match ────────────────────────────────────────────────────

def match_gt_to_mcseg(mask: np.ndarray, cell_labels: pd.DataFrame, regions: list) -> pd.DataFrame:
    log.info("Step 7: MCseg centroids → GT region assignment")
    cell_id_to_label = cell_labels.set_index("cell_id")["broad_label"].to_dict()

    rows_arr, cols_arr = np.where(mask > 0)
    cell_id_arr = mask[rows_arr, cols_arr]
    cent = (pd.DataFrame({"cell_id": cell_id_arr, "r": rows_arr, "c": cols_arr})
              .groupby("cell_id")[["r", "c"]].mean().reset_index())
    cent["cell_x"] = cent["c"] + CROP_X0
    cent["cell_y"] = cent["r"] + CROP_Y0
    log.info(f"  MCseg cells: {len(cent):,}")

    cent["gt_label"]   = assign_gt_labels(cent, regions, x_col="cell_x", y_col="cell_y")
    cent["pred_label"] = cent["cell_id"].map(cell_id_to_label).fillna("other")
    log.info(f"  in GT regions: {(cent['gt_label'] != 'unmatched').sum():,}")
    return cent


# ─── Step 8: F1 + figure ──────────────────────────────────────────────────────

def compute_f1_and_figure(cent: pd.DataFrame, sd: dict) -> dict:
    from sklearn.metrics import f1_score, precision_recall_fscore_support, confusion_matrix, classification_report
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    log.info("Step 8: F1 + figure")
    matched = cent[cent["gt_label"] != "unmatched"]
    valid   = matched[matched["pred_label"] != "other"]
    classes = sorted({c for c in ["epithelial", "stromal"] if c in valid["gt_label"].unique()})

    if valid.empty or not classes:
        log.warning("  insufficient valid cells")
        return {}

    y_true, y_pred = valid["gt_label"].values, valid["pred_label"].values
    micro    = f1_score(y_true, y_pred, labels=classes, average="micro",    zero_division=0)
    weighted = f1_score(y_true, y_pred, labels=classes, average="weighted", zero_division=0)
    log.info(f"  MCseg micro={micro:.3f}  weighted={weighted:.3f}  (StarDist micro={sd['micro_f1']:.3f})")
    log.info(f"\n{classification_report(y_true, y_pred, labels=classes, zero_division=0)}")

    prec, rec, f1_per, sup = precision_recall_fscore_support(
        y_true, y_pred, labels=classes, zero_division=0)
    pd.DataFrame({"class": classes, "precision": prec, "recall": rec, "f1": f1_per, "support": sup}
    ).to_csv(str(RESULTS_DIR / "f1_per_class_breast.csv"), index=False)

    summary = {
        "n_gt_total": len(cent), "n_matched": len(matched), "n_valid": len(valid),
        "micro_f1": micro, "weighted_f1": weighted,
        "stardist_micro_f1": sd["micro_f1"], "stardist_weighted_f1": sd["weighted_f1"],
        "delta_micro": micro - sd["micro_f1"],
    }
    pd.DataFrame([summary]).to_csv(str(RESULTS_DIR / "f1_summary_breast.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    x, w = np.arange(2), 0.35
    micro_v  = [sd["micro_f1"],    micro]
    weight_v = [sd["weighted_f1"], weighted]
    bars1 = ax.bar(x - w/2, micro_v,  w, label="Micro F1",    color=["#4e79a7", "#59a14f"])
    bars2 = ax.bar(x + w/2, weight_v, w, label="Weighted F1", color=["#76b7b2", "#b07aa1"])
    ax.set_ylim(0, 1.0)
    ax.set_xticks(x); ax.set_xticklabels(["StarDist\n(ENACT WBA)", "MCseg\n+WBA"], fontsize=9)
    ax.set_ylabel("F1 Score")
    ax.set_title(f"Cell-type F1 — Breast Cancer ENACT GT\n(n={len(valid):,})")
    ax.legend(fontsize=8)
    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.01, f"{h:.3f}",
                ha="center", va="bottom", fontsize=8)

    cm      = confusion_matrix(y_true, y_pred, labels=classes)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Greens",
                xticklabels=[c.capitalize() for c in classes],
                yticklabels=[c.capitalize() for c in classes],
                ax=axes[1], cbar_kws={"label": "Proportion"})
    axes[1].set_xlabel("Predicted"); axes[1].set_ylabel("True (GT)")
    axes[1].set_title(f"MCseg+WBA Confusion (Micro F1={micro:.3f})")
    plt.tight_layout()
    fig.savefig(str(RESULTS_DIR / "fig_f1_breast.png"), dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  Figure saved: fig_f1_breast.png")
    return summary


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    regions = load_gt_regions()
    sd      = compute_stardist_f1(regions)
    img     = crop_he(RESULTS_DIR / "he_crop_breast.tif")
    mask    = run_mcseg(img, RESULTS_DIR / "mcseg_mask_breast.npy")
    wba     = run_wba(mask, RESULTS_DIR / "wba_attribution_breast.csv")
    labels  = build_anndata_wba(wba, RESULTS_DIR / "celltypist_labels_breast.csv")
    cent    = match_gt_to_mcseg(mask, labels, regions)
    cent.to_csv(str(RESULTS_DIR / "gt_matched_breast.csv"), index=False)
    summary = compute_f1_and_figure(cent, sd)

    log.info(f"\n{'='*60}")
    log.info(f"完成 ({(time.time()-t0)/60:.1f} min)  → {RESULTS_DIR}")
    if summary:
        log.info(f"StarDist  micro F1 = {summary['stardist_micro_f1']:.3f}")
        log.info(f"MCseg+WBA micro F1 = {summary['micro_f1']:.3f}  (Δ {summary['delta_micro']:+.3f})")
    log.info("="*60)


if __name__ == "__main__":
    main()
