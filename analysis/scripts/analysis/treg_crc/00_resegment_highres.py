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
BTF_PATH  = Path(r"K:\plan_a\tissue sample\CRC\visium\official_v4\Visium_HD_Human_Colon_Cancer_tissue_image.btf")
H5_PATH   = Path(r"K:\plan_a\tissue sample\ENACT_supporting_files\public_data\human_colorectal\input_files\filtered_feature_bc_matrix.h5")
TP_PATH   = Path(r"K:\plan_a\tissue sample\ENACT_supporting_files\public_data\human_colorectal\input_files\tissue_positions.parquet")

MSSEG_ROOT  = Path(r"K:\plan_a\MSseg")
RESULTS_DIR = Path(r"K:\plan_a\submission_bioinformatics\analysis\treg_crc\results\highres_seg")

# GT coordinate range in ENACT-local BTF pixels (= cell_x/cell_y) + 200px buffer
CROP_X0, CROP_X1 = 5154, 15242   # ENACT local col  (= BTF col)
CROP_Y0, CROP_Y1 = 4635, 18599   # ENACT local row  (= pxl_row_in_fullres directly)

COL_OFFSET = 40598

LABEL_MAP = {
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

# ─── Step 1: Crop H&E from BTF ────────────────────────────────────────────────
def crop_he_from_btf(btf_path: Path, crop_tif: Path) -> np.ndarray:
    if crop_tif.exists():
        log.info(f"載入已存在的 H&E crop: {crop_tif.name}")
        img = tifffile.imread(str(crop_tif))
        if img.ndim == 3 and img.shape[-1] == 4:
            img = img[..., :3]
        log.info(f"  crop shape: {img.shape}")
        return img

    btf_col0 = CROP_X0 + COL_OFFSET
    btf_col1 = CROP_X1 + COL_OFFSET
    log.info(f"從 BTF 讀取 H&E crop (row {CROP_Y0}:{CROP_Y1}, BTF col {btf_col0}:{btf_col1})")
    t0 = time.time()
    with tifffile.TiffFile(str(btf_path)) as tif:
        store = tif.aszarr()
    z = zarr.open(store, mode="r")
    arr = z[0] if z.ndim == 4 else z
    img = np.asarray(arr[CROP_Y0:CROP_Y1, btf_col0:btf_col1])
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  crop shape: {img.shape}  ({time.time()-t0:.0f}s)")

    tifffile.imwrite(str(crop_tif), img, compression="zlib")
    log.info(f"  儲存: {crop_tif.name}")
    return img

# ─── Step 2: MCseg Segmentation ───────────────────────────────────────────────
def run_mcseg(img: np.ndarray, mask_npy: Path) -> np.ndarray:
    if mask_npy.exists():
        log.info(f"載入已存在的 MCseg 遮罩: {mask_npy.name}")
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
        "use_cpsam":             True,
        "voronoi_distance":      8,
        "flow_threshold":        0.4,
        "cellprob_threshold":    -2.0,
        "min_size":              20,
        "max_size":              6000,
        "clahe_clip_limit":      3.0,
        "use_transcript_rescue": False,
    }

    log.info("MCseg v2 tiled 分割 (deployment-mode, voronoi_d=8)")
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

# ─── Step 3: Bin Attribution ──────────────────────────────────────────────────
def run_bin_attribution(mask: np.ndarray, tp_path: Path) -> pd.DataFrame:
    log.info("Bin attribution (mask lookup)")
    tp = pd.read_parquet(str(tp_path), columns=[
        "barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    tp = tp[tp["in_tissue"] == 1]

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

# ─── Step 4: AnnData + CellTypist ─────────────────────────────────────────────
def build_anndata_and_annotate(
    attribution: pd.DataFrame,
    h5_path: Path,
    celltypist_csv: Path,
) -> pd.DataFrame:
    if celltypist_csv.exists():
        log.info(f"載入已存在的 CellTypist 結果: {celltypist_csv.name}")
        return pd.read_csv(celltypist_csv)

    log.info("建立 AnnData + CellTypist 標注")
    log.info(f"  讀取 h5 矩陣: {h5_path.name}")
    adata_full = sc.read_10x_h5(str(h5_path))
    adata_full.var_names_make_unique()
    log.info(f"  全片: {adata_full.n_obs:,} barcodes × {adata_full.n_vars:,} genes")

    barcodes_in_crop = attribution["barcode"].values
    mask_obs = adata_full.obs_names.isin(barcodes_in_crop)
    adata_crop = adata_full[mask_obs].copy()
    del adata_full
    gc.collect()

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
    
    # Extract confidence score (max probability for the predicted class)
    prob_matrix = predictions.probability_matrix
    conf_scores = prob_matrix.max(axis=1).values


    df = pd.DataFrame({
        "cell_id":          unique_cells,
        "celltypist_label": ct_labels,
        "conf_score":       conf_scores,
        "broad_label":      [LABEL_MAP.get(lbl, "other") for lbl in ct_labels],
    })

    df.to_csv(str(celltypist_csv), index=False)
    log.info(f"  儲存: {celltypist_csv.name}")
    log.info(f"  broad_label 分佈:\n{df['broad_label'].value_counts().to_string()}")
    return df


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    img         = crop_he_from_btf(BTF_PATH, RESULTS_DIR / "he_crop.tif")
    mask        = run_mcseg(img, RESULTS_DIR / "mcseg_mask.npy")
    del img
    gc.collect()
    attribution = run_bin_attribution(mask, TP_PATH)
    _           = build_anndata_and_annotate(attribution, H5_PATH,
                                             RESULTS_DIR / "celltypist_labels.csv")
    log.info("全部高解析分割與標註流程完成！")

if __name__ == "__main__":
    main()
