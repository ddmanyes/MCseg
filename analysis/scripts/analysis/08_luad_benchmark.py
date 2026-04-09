"""
10_luad_validation_v2.py
========================
LUAD 效標效度驗證 v2：6 種本質上不同的分割方法 × TAS_v2 Spearman 相關。

v1 的問題：同一 Cellpose detection + 不同膨脹距離
  → A1 增幅 163%，NED 僅 4.2%，TAS 被 Capture 單調主導 → ρ = -0.607

v2 修正：測試 6 種本質上不同的演算法，製造真正的跨方法 Capture / D1 / NED 差異：

  M0: nuclei model, dia=8, 無膨脹        (AP 預期 ~0.20-0.30)
  M1: cyto3 dia=10, CLAHE, expand=3     (AP 預期 ~0.35-0.45)
  M2: cyto3 dia=17, CLAHE, expand=5     (AP 預期 ~0.50-0.55)
  M3: cyto3 dia=17, CLAHE, Voronoi=12   (AP 預期 ~0.58-0.62)
  M4: cyto3 多尺度[13,17,22], Voronoi=14 (AP 預期 ~0.62-0.65)
  M5: 色彩解卷積+多尺度, Adap-Voronoi=18 (AP 預期 ~0.65, V12-style)

TAS_v2：
  Capture = mean(A1_norm, A2_norm, A3_norm)
  Purity  = √(C1_norm × NED_norm)
  Core    = √(Capture × Purity)
  TAS_v2  = 0.65×Core + 0.20×D1_norm + 0.15×E1_norm

中間結果（masks）快取於 results/cache/luad_v2/，方便重跑。

輸出：
  results/metrics/luad_validation_v2.csv
  results/figures/fig_spearman_v2.png
"""

from __future__ import annotations

import sys
import time
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS   = cfg["paths"]
PLT_CFG = cfg["plotting"]
TAS_CFG = cfg["tas"]

FIGURES_DIR = ROOT / PATHS["figures_dir"]
METRICS_DIR = ROOT / PATHS["metrics_dir"]
CACHE_DIR   = ROOT / "results" / "cache" / "luad_v2"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
METRICS_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

DPI = PLT_CFG["dpi"]

# ── LUAD 資料路徑 ──────────────────────────────────────────────────────────

AUTORESEARCH_DIR = Path("/Volumes/SSD/plan_a/autoresearch_seg")
LUAD_DATA_DIR    = AUTORESEARCH_DIR / "data"
LUAD_BASE        = Path("/Volumes/SSD/plan_a/tissue sample/LUAD")
LUAD_POSITIONS   = LUAD_BASE / "visium/binned_outputs/square_002um/spatial/tissue_positions.parquet"
LUAD_MATRIX_H5   = LUAD_BASE / "visium/binned_outputs/square_002um/filtered_feature_bc_matrix.h5"

# LUAD patch ROI 座標（由 autoresearch_seg/prepare.py 計算）
COL_START = 8733
ROW_START = 4343
COL_END   = 10207
ROW_END   = 5817

# LUAD 特異性基因對（生物學不可能共現）
LUAD_IMPOSSIBLE_PAIRS = [
    ["EPCAM", "CD3E"],    # 上皮 vs T cell
    ["NAPSA", "CD68"],    # 肺腺癌特異 vs 巨噬細胞
    ["SFTPC", "NKG7"],   # 肺泡 II 型 vs NK/T cell
    ["KRT7",  "CD3E"],   # 腺癌角蛋白 vs T cell
]
LUAD_IMMUNE_MARKERS = ["CD3E", "NKG7", "MS4A1"]

# TAS 權重
UMI_REF  = TAS_CFG["umi_ref"]           # 1000
GENE_REF = TAS_CFG["genes_ref"]          # 300
W_CORE   = TAS_CFG["weight_core_tas"]    # 0.65
W_BIO    = TAS_CFG["weight_biology"]     # 0.20
W_IMM    = TAS_CFG["weight_immune"]      # 0.15
MIN_UMI  = cfg["qc"]["min_umi"]          # 50


# ── 6 種方法定義 ──────────────────────────────────────────────────────────

METHODS = [
    {
        "name":            "M0: nuclei d=8",
        "model":           "nuclei",
        "multi_passes":    None,          # None = single pass
        "dia":             8,
        "cellprob":        -1.0,
        "flow_thresh":     0.4,
        "use_clahe":       False,
        "use_color_deconv":False,
        "expand_mode":     "none",
        "expand_d":        0,
        "desc":            "Nuclei model tiny detection, no expansion",
    },
    {
        "name":            "M1: cyto3 d=10 exp=3",
        "model":           "cyto3",
        "multi_passes":    None,
        "dia":             10,
        "cellprob":        -2.0,
        "flow_thresh":     0.4,
        "use_clahe":       True,
        "use_color_deconv":False,
        "expand_mode":     "expand",
        "expand_d":        3,
        "desc":            "Small cyto3, minimal expand_labels",
    },
    {
        "name":            "M2: cyto3 d=17 exp=5",
        "model":           "cyto3",
        "multi_passes":    None,
        "dia":             17,
        "cellprob":        -2.0,
        "flow_thresh":     0.4,
        "use_clahe":       True,
        "use_color_deconv":False,
        "expand_mode":     "expand",
        "expand_d":        5,
        "desc":            "Standard cyto3, simple expand_labels",
    },
    {
        "name":            "M3: cyto3 d=17 Vor=12",
        "model":           "cyto3",
        "multi_passes":    None,
        "dia":             17,
        "cellprob":        -2.0,
        "flow_thresh":     0.4,
        "use_clahe":       True,
        "use_color_deconv":False,
        "expand_mode":     "voronoi",
        "expand_d":        12,
        "desc":            "cyto3 + Voronoi-constrained expansion",
    },
    {
        "name":            "M4: multi-dia Vor=14",
        "model":           "cyto3",
        "multi_passes":    [13, 17, 22],  # list of diameters
        "dia":             17,            # fallback single dia
        "cellprob":        -2.0,
        "flow_thresh":     0.4,
        "use_clahe":       True,
        "use_color_deconv":False,
        "expand_mode":     "voronoi",
        "expand_d":        14,
        "desc":            "Multi-scale cyto3 ensemble + Voronoi",
    },
    {
        "name":            "M5: colordeconv Vor=18",
        "model":           "cyto3",
        "multi_passes":    [13, 17, 22],
        "dia":             17,
        "cellprob":        -2.0,
        "flow_thresh":     0.4,
        "use_clahe":       True,
        "use_color_deconv":True,          # adds hematoxylin pass
        "expand_mode":     "adaptive_voronoi",
        "expand_d":        18,
        "desc":            "V12-style: color deconv + multi-scale + adaptive Voronoi",
    },
]


# ── 影像處理工具（複用自 segment_best.py）──────────────────────────────────

def apply_clahe(img: np.ndarray, clip_limit: float = 3.0, tile_size: int = 8) -> np.ndarray:
    lab = cv2.cvtColor(img, cv2.COLOR_RGB2LAB)
    l, a, b = cv2.split(lab)
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(tile_size, tile_size))
    cl = clahe.apply(l)
    return cv2.cvtColor(cv2.merge((cl, a, b)), cv2.COLOR_LAB2RGB)


def color_deconvolution_he(img: np.ndarray) -> np.ndarray:
    """萃取 Hematoxylin 通道（Ruifrok & Johnston）。"""
    img_float = img.astype(np.float64) + 1.0
    od = -np.log(img_float / 256.0)
    he_matrix = np.array([
        [0.6500286, 0.7041680, 0.2860126],
        [0.0728940, 0.9904310, 0.1155140],
        [0.2688350, 0.5706770, 0.7768750],
    ])
    for i in range(3):
        n = np.linalg.norm(he_matrix[i])
        if n > 0:
            he_matrix[i] /= n
    stains = od.reshape(-1, 3) @ np.linalg.inv(he_matrix).T
    hema = stains.reshape(img.shape)[:, :, 0]
    hema = np.clip(hema, 0, None)
    h_max = np.percentile(hema, 99.5)
    if h_max > 0:
        hema = np.clip(hema / h_max, 0, 1)
    return (hema * 255).astype(np.uint8)


def create_tissue_mask(img: np.ndarray) -> np.ndarray:
    gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    tissue = gray < 220
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
    tissue = cv2.morphologyEx(tissue.astype(np.uint8), cv2.MORPH_CLOSE, k)
    k2 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    tissue = cv2.dilate(tissue, k2, iterations=2)
    return tissue.astype(bool)


def voronoi_expand(mask: np.ndarray, max_distance: int,
                   tissue_mask: np.ndarray = None) -> np.ndarray:
    from scipy.ndimage import distance_transform_edt
    binary = mask > 0
    if not binary.any():
        return mask.copy()
    dist, nearest = distance_transform_edt(~binary, return_indices=True)
    expanded = mask[nearest[0], nearest[1]]
    expanded[dist > max_distance] = 0
    if tissue_mask is not None:
        expanded[~tissue_mask] = 0
    return expanded.astype(np.int32)


def adaptive_voronoi_expand(mask: np.ndarray, base_distance: int,
                             tissue_mask: np.ndarray = None) -> np.ndarray:
    from scipy.ndimage import distance_transform_edt
    binary = mask > 0
    if not binary.any():
        return mask.copy()
    dist, nearest = distance_transform_edt(~binary, return_indices=True)
    expanded = mask[nearest[0], nearest[1]]
    density = cv2.GaussianBlur(binary.astype(np.float32), (0, 0), sigmaX=25)
    d_max = density.max()
    if d_max > 0:
        density /= d_max
    adaptive_max = base_distance * (1.5 - 0.9 * density)
    expanded[dist > adaptive_max] = 0
    if tissue_mask is not None:
        expanded[~tissue_mask] = 0
    return expanded.astype(np.int32)


def merge_masks_fast(base: np.ndarray, new: np.ndarray,
                     max_overlap: float = 0.15, min_size: int = 15) -> np.ndarray:
    """非重疊合併：將 new 中未被 base 覆蓋的細胞加入 base。"""
    next_id = int(base.max()) + 1
    base_occ = base > 0
    for nid in np.unique(new):
        if nid == 0:
            continue
        pix = new == nid
        cnt = pix.sum()
        if cnt < min_size:
            continue
        if (base_occ & pix).sum() / cnt < max_overlap:
            empty = pix & (~base_occ)
            if empty.sum() >= min_size:
                base[empty] = next_id
                base_occ[empty] = True
                next_id += 1
    return base


def clean_mask(mask: np.ndarray, min_size: int = 10,
               max_size: int = 6000) -> np.ndarray:
    ids, cnts = np.unique(mask, return_counts=True)
    remove = ids[((cnts < min_size) | (cnts > max_size)) & (ids > 0)]
    if len(remove):
        mask[np.isin(mask, remove)] = 0
    return mask


# ── Cellpose 分割（依方法設定）─────────────────────────────────────────────

def run_method(he_patch: np.ndarray, mcfg: dict) -> np.ndarray:
    """依 mcfg 執行 Cellpose + 膨脹，回傳最終遮罩。"""
    from cellpose import models
    from skimage.segmentation import expand_labels

    print(f"\n  [{mcfg['name']}] {mcfg['desc']}")

    # 1. 前處理
    if mcfg["use_clahe"]:
        img = apply_clahe(he_patch)
    else:
        img = he_patch.copy()
    tissue_mask = create_tissue_mask(he_patch)

    # Hematoxylin（M5 用）
    hema = None
    if mcfg["use_color_deconv"]:
        h_raw = color_deconvolution_he(he_patch)
        hema = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(h_raw)

    # 2. 偵測
    t0 = time.time()
    model = models.CellposeModel(model_type=mcfg["model"], gpu=True)

    nucleus_mask = np.zeros(he_patch.shape[:2], dtype=np.int32)

    if mcfg["multi_passes"] is not None:
        # 多尺度偵測
        for dia in mcfg["multi_passes"]:
            m, _, _ = model.eval(
                img,
                diameter=float(dia),
                channels=[0, 0],
                flow_threshold=mcfg["flow_thresh"],
                cellprob_threshold=mcfg["cellprob"],
                min_size=8,
                augment=(dia == 17),
                resample=True,
            )
            nucleus_mask = merge_masks_fast(nucleus_mask, m.astype(np.int32))
            print(f"    cyto3 dia={dia}: {m.max()} cells → merged total: {nucleus_mask.max()}")

        # 色彩解卷積加一個 hematoxylin pass
        if hema is not None:
            m_h, _, _ = model.eval(
                hema,
                diameter=17.0,
                channels=[0, 0],
                flow_threshold=mcfg["flow_thresh"],
                cellprob_threshold=mcfg["cellprob"],
                min_size=8,
                augment=True,
                resample=True,
            )
            nucleus_mask = merge_masks_fast(nucleus_mask, m_h.astype(np.int32))
            print(f"    cyto3 hema: {m_h.max()} cells → merged total: {nucleus_mask.max()}")
    else:
        # 單尺度偵測
        # M0 (nuclei) → 使用灰度影像
        if mcfg["model"] == "nuclei":
            gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
            input_img = gray
        else:
            input_img = img

        m, _, _ = model.eval(
            input_img,
            diameter=float(mcfg["dia"]),
            channels=[0, 0],
            flow_threshold=mcfg["flow_thresh"],
            cellprob_threshold=mcfg["cellprob"],
            min_size=5 if mcfg["model"] == "nuclei" else 10,
            augment=(mcfg["model"] != "nuclei"),
            resample=True,
        )
        nucleus_mask = m.astype(np.int32)
        print(f"    {mcfg['model']} dia={mcfg['dia']}: {nucleus_mask.max()} cells")

    print(f"    偵測耗時 {time.time()-t0:.0f}s，共 {nucleus_mask.max()} 個細胞核")

    # 組織遮罩過濾
    nucleus_mask[~tissue_mask] = 0

    # 3. 膨脹
    mode = mcfg["expand_mode"]
    d    = mcfg["expand_d"]

    if mode == "expand" and d > 0:
        final_mask = expand_labels(nucleus_mask, distance=d)
        final_mask[~tissue_mask] = 0
        final_mask = final_mask.astype(np.int32)
    elif mode == "voronoi" and d > 0:
        final_mask = voronoi_expand(nucleus_mask, max_distance=d, tissue_mask=tissue_mask)
    elif mode == "adaptive_voronoi" and d > 0:
        final_mask = adaptive_voronoi_expand(nucleus_mask, base_distance=d,
                                              tissue_mask=tissue_mask)
    else:
        final_mask = nucleus_mask.copy()

    final_mask = clean_mask(final_mask.astype(np.int32))
    print(f"    膨脹後：{final_mask.max()} 個細胞")
    return final_mask


# ── AP@0.5 ────────────────────────────────────────────────────────────────

def compute_ap05(pred: np.ndarray, gt: np.ndarray) -> float:
    from cellpose import metrics
    if pred.shape != gt.shape:
        pred = cv2.resize(pred, (gt.shape[1], gt.shape[0]),
                          interpolation=cv2.INTER_NEAREST)
    ap, _, _, _ = metrics.average_precision(gt, pred, threshold=[0.5])
    return float(ap[0])


# ── NED ───────────────────────────────────────────────────────────────────

def compute_ned(mask: np.ndarray, cell_adata, n_hvgs: int = 1000) -> float:
    """
    Neighbor Expression Divergence：相鄰細胞對的 Hellinger 距離均值。
    高 NED → 邊界清晰（好）；低 NED → 轉錄本混入（差）。
    Returns np.nan 若資料不足。
    """
    from scipy.ndimage import grey_dilation

    if cell_adata is None or cell_adata.n_obs < 5 or int(mask.max()) < 2:
        return np.nan

    X = cell_adata.X
    if sp.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float32)
    else:
        X = np.array(X, dtype=np.float32)

    if X.shape[1] > n_hvgs:
        hvg_idx = np.argpartition(X.var(axis=0), -n_hvgs)[-n_hvgs:]
        X = X[:, hvg_idx]

    row_sums = np.maximum(X.sum(axis=1, keepdims=True), 1e-10)
    X_prob = (X / row_sums).astype(np.float32)

    obs_ids = cell_adata.obs["cell_id"].values
    cid2row = {int(c): i for i, c in enumerate(obs_ids)}

    struct   = np.ones((3, 3), dtype=np.int32)
    dilated  = grey_dilation(mask.astype(np.int32), footprint=struct)
    bnd      = (mask > 0) & (dilated != mask)
    ci       = mask[bnd].astype(np.int32)
    cj       = dilated[bnd].astype(np.int32)
    valid    = (cj > 0) & (cj != ci)
    ci, cj   = ci[valid], cj[valid]
    if len(ci) == 0:
        return np.nan

    pairs = np.unique(np.sort(np.stack([ci, cj], axis=1), axis=1), axis=0)
    known = set(cid2row.keys())
    ok    = np.array([int(a) in known and int(b) in known
                      for a, b in pairs])
    pairs = pairs[ok]
    if len(pairs) < 5:
        return np.nan

    if len(pairs) > 3000:
        rng   = np.random.default_rng(42)
        pairs = pairs[rng.choice(len(pairs), 3000, replace=False)]

    i_idx = np.array([cid2row[int(a)] for a in pairs[:, 0]])
    j_idx = np.array([cid2row[int(b)] for b in pairs[:, 1]])
    Xi    = X_prob[i_idx]
    Xj    = X_prob[j_idx]
    hell  = np.sqrt(np.sum((np.sqrt(np.maximum(Xi, 0)) -
                            np.sqrt(np.maximum(Xj, 0))) ** 2, axis=1) / 2)
    return float(np.clip(np.mean(hell), 0, 1))


# ── D1 per-method ─────────────────────────────────────────────────────────

def compute_d1(cell_adata, n_hvgs: int = 500, n_pcs: int = 10) -> float:
    """
    Per-method Silhouette score：HVG → normalize → PCA → Leiden → silhouette。
    Returns (sil+1)/2 ∈ [0,1]，失敗時返回 np.nan。
    """
    import scanpy as sc
    from sklearn.metrics import silhouette_score

    if cell_adata is None or cell_adata.n_obs < 50:
        return np.nan
    try:
        adata = cell_adata.copy()
        sc.pp.filter_cells(adata, min_counts=MIN_UMI)
        if adata.n_obs < 50:
            return np.nan

        n_h = min(n_hvgs, adata.n_vars - 1, adata.n_obs - 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(adata, n_top_genes=n_h,
                                        flavor="seurat_v3", span=0.5)
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        n_pc = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=n_pc)
        sc.pp.neighbors(adata, n_neighbors=min(10, adata.n_obs // 5), n_pcs=n_pc)
        sc.tl.leiden(adata, resolution=0.3, flavor="igraph",
                     directed=False, n_iterations=2)

        if adata.obs["leiden"].nunique() < 2:
            return np.nan

        sil = silhouette_score(adata.obsm["X_pca"][:, :n_pc],
                               adata.obs["leiden"])
        return float(np.clip((sil + 1) / 2, 0, 1))

    except Exception as exc:
        print(f"    [D1] 計算失敗（{exc}），返回 NaN")
        return np.nan


# ── VHD 資料載入 ──────────────────────────────────────────────────────────

def load_luad_vhd():
    import scanpy as sc

    print("  載入 LUAD tissue positions...")
    tp      = pd.read_parquet(LUAD_POSITIONS)
    roi_bins = tp[
        (tp["pxl_col_in_fullres"] >= COL_START) &
        (tp["pxl_col_in_fullres"] <  COL_END) &
        (tp["pxl_row_in_fullres"] >= ROW_START) &
        (tp["pxl_row_in_fullres"] <  ROW_END) &
        (tp["in_tissue"] == 1)
    ].copy()
    roi_bins["local_col"] = roi_bins["pxl_col_in_fullres"] - COL_START
    roi_bins["local_row"] = roi_bins["pxl_row_in_fullres"] - ROW_START
    print(f"  ROI 內 in-tissue bins: {len(roi_bins)}")

    print("  載入全片表達矩陣 (h5)...")
    adata_full = sc.read_10x_h5(str(LUAD_MATRIX_H5))
    adata_full.var_names_make_unique()

    bc2idx     = {bc: i for i, bc in enumerate(adata_full.obs_names)}
    valid_bc   = [bc for bc in roi_bins["barcode"] if bc in bc2idx]
    roi_indices = [bc2idx[bc] for bc in valid_bc]
    roi_bins_valid = roi_bins.set_index("barcode").loc[valid_bc]
    return roi_bins_valid, adata_full, roi_indices


# ── Bin → Cell AnnData 聚合 ───────────────────────────────────────────────

def attribute_and_aggregate(mask, roi_bins, adata_full, roi_indices):
    import anndata as ad
    import scanpy as sc

    lrows = roi_bins["local_row"].values.astype(int).clip(0, mask.shape[0] - 1)
    lcols = roi_bins["local_col"].values.astype(int).clip(0, mask.shape[1] - 1)
    cids  = mask[lrows, lcols]

    adata_roi      = adata_full[roi_indices].copy()
    adata_roi.obs["cell_id"] = cids.astype(np.int32)

    valid = cids > 0
    if valid.sum() == 0:
        return None

    adata_v    = adata_roi[valid]
    cids_v     = cids[valid]
    u_cells    = np.unique(cids_v)
    cid2i      = {c: i for i, c in enumerate(u_cells)}
    row_idx    = np.array([cid2i[c] for c in cids_v])
    A          = sp.csr_matrix(
        (np.ones(len(cids_v), dtype=np.float32),
         (row_idx, np.arange(len(cids_v)))),
        shape=(len(u_cells), adata_v.n_obs)
    )
    Xv = adata_v.X
    if not sp.issparse(Xv):
        Xv = sp.csr_matrix(Xv)
    X_cell = A @ Xv

    cell_adata = ad.AnnData(
        X   = X_cell.astype(np.float32),
        var = adata_roi.var.copy(),
        obs = pd.DataFrame({"cell_id": u_cells},
                           index=[f"cell_{c}" for c in u_cells])
    )
    sc.pp.calculate_qc_metrics(cell_adata, percent_top=None,
                                log1p=False, inplace=True)
    mt = cell_adata.var_names.str.upper().str.startswith("MT-")
    if mt.sum() > 0:
        mt_expr = np.asarray(cell_adata[:, mt].X.sum(axis=1)).ravel()
        tot     = np.asarray(cell_adata.X.sum(axis=1)).ravel()
        cell_adata.obs["pct_mt"] = np.where(tot > 0, mt_expr / tot * 100, 0.0)
    else:
        cell_adata.obs["pct_mt"] = 0.0

    cell_adata.obs.rename(columns={
        "total_counts":     "n_umis",
        "n_genes_by_counts":"n_genes"
    }, inplace=True)
    return cell_adata


# ── TAS 計算 ─────────────────────────────────────────────────────────────

def compute_metrics(mask, roi_bins, adata_full, roi_indices,
                    n_total_bins: int) -> dict:
    """計算單一 mask 的全部 TAS_v2 指標。"""

    # A1
    lrows = roi_bins["local_row"].values.astype(int).clip(0, mask.shape[0] - 1)
    lcols = roi_bins["local_col"].values.astype(int).clip(0, mask.shape[1] - 1)
    a1 = float((mask[lrows, lcols] > 0).sum() / n_total_bins)

    # Cell-level AnnData
    cell_adata = attribute_and_aggregate(mask, roi_bins, adata_full, roi_indices)
    if cell_adata is None or cell_adata.n_obs < 5:
        return {k: np.nan for k in
                ["a1","a2","a3","c1","ned","d1","e1",
                 "capture","purity","core_tas","d1_norm","e1_norm","tas_v2"]}

    umis  = cell_adata.obs["n_umis"].values
    genes = cell_adata.obs["n_genes"].values
    a2    = float(np.median(umis))
    a3    = float(np.median(genes))

    # C1: artificial co-expression
    rates = []
    for ga, gb in LUAD_IMPOSSIBLE_PAIRS:
        if ga not in cell_adata.var_names or gb not in cell_adata.var_names:
            continue
        X    = cell_adata.X
        ia   = cell_adata.var_names.get_loc(ga)
        ib   = cell_adata.var_names.get_loc(gb)
        ca   = np.asarray(X[:, ia].todense()).ravel() if sp.issparse(X) else X[:, ia]
        cb   = np.asarray(X[:, ib].todense()).ravel() if sp.issparse(X) else X[:, ib]
        rates.append(float(((ca > 0) & (cb > 0)).sum() / cell_adata.n_obs))
    c1 = float(np.mean(rates)) if rates else np.nan

    # NED
    try:
        ned = compute_ned(mask, cell_adata)
    except Exception as exc:
        print(f"    [NED] 失敗({exc})，設 NaN")
        ned = np.nan

    # D1 per-method clustering silhouette
    d1 = compute_d1(cell_adata)

    # E1: immune marker survival
    present_m = [m for m in LUAD_IMMUNE_MARKERS if m in cell_adata.var_names]
    if present_m:
        X = cell_adata.X
        immune = np.zeros(cell_adata.n_obs, dtype=bool)
        for m in present_m:
            col = X[:, cell_adata.var_names.get_loc(m)]
            if sp.issparse(col):
                col = np.asarray(col.todense()).ravel()
            immune |= (col > 0)
        tot_imm = immune.sum()
        e1 = float((immune & (umis >= MIN_UMI)).sum() / max(tot_imm, 1)) \
             if tot_imm > 0 else np.nan
    else:
        e1 = np.nan

    # 歸一化
    a1_n   = float(np.clip(a1, 0, 1))
    a2_n   = float(np.clip(np.log(max(a2, 1)) / np.log(UMI_REF), 0, 1))
    a3_n   = float(np.clip(a3 / GENE_REF, 0, 1))
    c1_n   = float(np.clip(1 - (c1 if not np.isnan(c1) else 0), 0, 1))
    ned_n  = float(np.clip(ned, 0, 1)) if not np.isnan(ned) else c1_n
    d1_n   = float(np.clip(d1, 0, 1)) if not np.isnan(d1) else 0.5
    e1_n   = float(np.clip(e1, 0, 1)) if not np.isnan(e1) else 0.5

    capture = float(np.mean([a1_n, a2_n, a3_n]))
    purity  = float(np.sqrt(c1_n * ned_n))        # NED + C1 幾何平均
    core    = float(np.sqrt(capture * purity))
    tas_v2  = float(W_CORE * core + W_BIO * d1_n + W_IMM * e1_n)

    return {
        "a1": a1, "a2": a2, "a3": a3, "c1": c1, "ned": ned, "d1": d1, "e1": e1,
        "a1_n": a1_n, "a2_n": a2_n, "a3_n": a3_n,
        "c1_n": c1_n, "ned_n": ned_n,
        "capture": capture, "purity": purity, "core_tas": core,
        "d1_norm": d1_n, "e1_norm": e1_n,
        "tas_v2": tas_v2,
        "n_cells": int(mask.max()),
    }


# ── 主流程 ────────────────────────────────────────────────────────────────

def run_validation():
    print("=== LUAD 效標效度驗證 v2（6 種真正不同的分割方法）===")

    print("\n1. 載入基礎資料...")
    he_patch = np.load(LUAD_DATA_DIR / "he_patch.npy")
    gt_mask  = np.load(LUAD_DATA_DIR / "gt_mask.npy")
    roi_bins, adata_full, roi_indices = load_luad_vhd()
    n_total = len(roi_bins)
    print(f"  HE patch: {he_patch.shape}, GT cells: {gt_mask.max()}")

    print("\n2. 逐方法執行分割 + 指標計算...")
    results = []

    for mcfg in METHODS:
        name_safe = mcfg["name"].replace(" ", "_").replace(":", "").replace("/", "_")
        mask_cache = CACHE_DIR / f"{name_safe}.npy"

        # 載入快取 or 重新運行
        if mask_cache.exists():
            mask = np.load(mask_cache)
            print(f"\n  [{mcfg['name']}] [CACHE] 載入 {mask_cache.name}，{mask.max()} 個細胞")
        else:
            t_seg = time.time()
            mask  = run_method(he_patch, mcfg)
            np.save(mask_cache, mask)
            print(f"    → 已快取 {mask_cache.name}（{time.time()-t_seg:.0f}s）")

        # AP@0.5
        t0   = time.time()
        ap05 = compute_ap05(mask, gt_mask)

        # TAS_v2 指標
        metrics = compute_metrics(mask, roi_bins, adata_full, roi_indices, n_total)

        elapsed = time.time() - t0
        row = {
            "method": mcfg["name"],
            "desc":   mcfg["desc"],
            "ap05":   ap05,
            **metrics,
            "elapsed_s": elapsed,
        }
        results.append(row)

        print(f"  {'':3s} AP@0.5={ap05:.4f}  TAS_v2={metrics['tas_v2']:.4f}  "
              f"Capture={metrics['capture']:.3f}  NED={metrics.get('ned',np.nan):.3f}  "
              f"D1={metrics.get('d1',np.nan):.3f}  N_cells={metrics.get('n_cells',0)}  "
              f"[{elapsed:.0f}s]")

    df = pd.DataFrame(results)
    out_csv = METRICS_DIR / "luad_validation_v2.csv"
    df.to_csv(out_csv, index=False)

    # Spearman
    valid = df.dropna(subset=["ap05", "tas_v2"])
    print(f"\n=== Spearman 相關 (n={len(valid)}) ===")
    if len(valid) >= 4:
        rho, pval = stats.spearmanr(valid["ap05"], valid["tas_v2"])
        verdict = (
            "✅ ρ ≥ 0.786，達 p<0.05，具效標效度" if abs(rho) >= 0.786
            else "⚠️  趨勢存在但未達顯著" if abs(rho) >= 0.60
            else "❌ 相關不足（需重新設計）"
        )
        print(f"  ρ = {rho:.3f},  p = {pval:.4f}  {verdict}")

        # 各指標的 Spearman
        for col in ["capture", "core_tas", "d1_norm", "tas_v2"]:
            if col in valid.columns:
                r, p = stats.spearmanr(valid["ap05"], valid[col])
                print(f"    ρ(AP@0.5, {col:12s}) = {r:+.3f}  p={p:.4f}")
    else:
        rho, pval = np.nan, np.nan
        print("  ⚠️  有效樣本不足")

    # 排名比較
    print("\n=== 方法排名 ===")
    display_cols = ["method", "ap05", "tas_v2", "capture", "purity", "d1_norm", "n_cells"]
    display_cols = [c for c in display_cols if c in df.columns]
    print(df[display_cols].sort_values("ap05").to_string(index=False))

    # 圖表
    _plot_results(df, rho, pval)

    print(f"\n✅ 10_luad_validation_v2.py 完成")
    print(f"   CSV  → {out_csv}")
    print(f"   圖表 → {FIGURES_DIR / 'fig_spearman_v2.png'}")
    return df, rho, pval


def _plot_results(df: pd.DataFrame, rho: float, pval: float):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    colors = plt.cm.RdYlGn(np.linspace(0.1, 0.9, len(df)))
    methods = df["method"].tolist()

    # Panel 1: AP@0.5 vs TAS_v2 scatter
    ax = axes[0]
    for i, (_, row) in enumerate(df.iterrows()):
        ax.scatter(row["ap05"], row["tas_v2"], color=colors[i], s=120, zorder=3)
        ax.annotate(row["method"], (row["ap05"], row["tas_v2"]),
                    textcoords="offset points", xytext=(5, 3), fontsize=8)
    lbl = f"ρ = {rho:.3f},  p = {pval:.4f}" if not np.isnan(rho) else ""
    ax.set_xlabel("AP@0.5（Xenium GT）", fontsize=11)
    ax.set_ylabel("TAS_v2", fontsize=11)
    ax.set_title(f"AP@0.5 vs TAS_v2\n{lbl}", fontsize=11)
    ax.grid(True, alpha=0.3)

    # Panel 2: 各分項指標
    ax = axes[1]
    x = np.arange(len(df))
    w = 0.2
    bars_ap   = ax.bar(x - 1.5*w, df["ap05"],    w, label="AP@0.5",   color="steelblue")
    bars_cap  = ax.bar(x - 0.5*w, df["capture"],  w, label="Capture",  color="forestgreen")
    bars_d1   = ax.bar(x + 0.5*w, df["d1_norm"],  w, label="D1_norm",  color="darkorange")
    bars_tas  = ax.bar(x + 1.5*w, df["tas_v2"],   w, label="TAS_v2",   color="crimson")
    ax.set_xticks(x)
    ax.set_xticklabels([m.split(":")[0] for m in methods], fontsize=9)
    ax.set_ylabel("Score", fontsize=11)
    ax.set_title("各方法分項指標比較", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # Panel 3: Ranking comparison
    ax = axes[2]
    rank_ap   = df["ap05"].rank().values
    rank_tas  = df["tas_v2"].rank().values
    for i, name in enumerate(methods):
        ax.plot([0, 1], [rank_ap[i], rank_tas[i]],
                "o-", color=colors[i], linewidth=2, markersize=8, label=name)
    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["AP@0.5 rank", "TAS_v2 rank"], fontsize=10)
    ax.set_ylabel("Rank (1=lowest)", fontsize=11)
    ax.set_title("排名一致性（完美相關 = 平行線）", fontsize=11)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        f"LUAD Criterion Validity v2 — 6 真正不同方法\n"
        f"Spearman ρ(AP@0.5, TAS_v2) = {rho:.3f}  (p={pval:.4f})",
        fontsize=13, y=1.02
    )
    plt.tight_layout()
    out = FIGURES_DIR / "fig_spearman_v2.png"
    plt.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  fig_spearman_v2.png → {out}")


if __name__ == "__main__":
    run_validation()
