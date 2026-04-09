"""
03_metrics_ac.py
================
計算 A 類（捕獲）與 C 類（純淨）指標：

  A1: Transcript capture rate（in-tissue bins 落在遮罩內的比例）
  A2: Median UMI/cell
  A3: Median genes/cell

  C1: Artificial co-expression rate（4 個生物不可能基因對）
  B1: Debris rate（UMI < 50 比例，輔助參考，不進 TAS）
  B2: QC yield（通過 QC 的細胞數，輔助參考，不進 TAS）

輸出：results/metrics/metrics_ac.csv
"""

from __future__ import annotations

import sys
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from pathlib import Path
from tqdm import tqdm

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS   = cfg["paths"]
DATA    = cfg["data"]
QC_CFG  = cfg["qc"]
PROFILE = cfg["tissue_profiles"]["CRC"]

ATTRIBUTION_DIR = ROOT / PATHS["attribution_dir"]
ANNDATA_DIR     = ROOT / PATHS["anndata_dir"]
METRICS_DIR     = ROOT / PATHS["metrics_dir"]
METRICS_DIR.mkdir(parents=True, exist_ok=True)

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

METHODS = DATA["methods"]
IMPOSSIBLE_PAIRS = PROFILE["impossible_pairs"]  # List of [gene_a, gene_b]
MASKS_DIR = ROOT / PATHS["masks_dir"]

PX_SIZE_UM = DATA["pixel_size_um"]   # 0.2738 µm/px
AREA_SCALE = PX_SIZE_UM ** 2         # µm² per pixel (= 0.07497 µm²)

# ── A1：需要 attribution parquet（含 in-tissue bins）─────────────────────

def compute_a1(method: str, roi_name: str) -> float:
    """Transcript capture rate (FTC) = (cell_id > 0 bins) / (all in-tissue bins in ROI)."""
    attr_path = ATTRIBUTION_DIR / f"{method}_{roi_name}.parquet"
    if not attr_path.exists():
        return np.nan
    attr = pd.read_parquet(attr_path)
    if len(attr) == 0:
        return np.nan
    return float((attr["cell_id"] > 0).sum() / len(attr))


def compute_umi_density(mask: np.ndarray, adata: ad.AnnData) -> float:
    """
    UMI density = median(total_UMI / cell_area_µm²) across all QC-passed cells.

    優於 FTC 之處：按面積歸一化，消除大細胞天然覆蓋更多 bins 的偏差。
    SR 的大多邊形在此指標中不再自動領先；NUC 的緊密核遮罩也能獲得公平評估。

    Returns: median UMI density [UMI/µm²]，若無法計算則 nan
    """
    if mask is None or adata is None or adata.n_obs == 0:
        return np.nan

    obs_cell_ids = adata.obs["cell_id"].values.astype(int)
    umis = adata.obs["n_umis"].values

    # 向量化計算 mask 中每個 cell_id 的像素數
    max_id = int(mask.max())
    if max_id < 1:
        return np.nan
    pixel_counts = np.bincount(mask.ravel(), minlength=max_id + 1)

    # 對應到 adata 中每顆細胞
    areas_px = np.where(
        obs_cell_ids <= max_id,
        pixel_counts[obs_cell_ids.clip(0, max_id)],
        0
    )
    areas_um2 = areas_px.astype(float) * AREA_SCALE

    valid = (areas_um2 > 0) & (umis > 0)
    if valid.sum() == 0:
        return np.nan

    density = umis[valid].astype(float) / areas_um2[valid]
    return float(np.median(density))


# ── A2, A3, C1, B1, B2：從 h5ad 計算 ────────────────────────────────────

def compute_metrics_from_adata(adata: ad.AnnData) -> dict:
    """
    計算 A2, A3, C1, B1, B2。
    輸入：cell-level AnnData（已含 obs.n_umis, obs.n_genes）
    """
    n_cells = adata.n_obs
    umis  = adata.obs["n_umis"].values
    genes = adata.obs["n_genes"].values

    # A2：中位 UMI/cell
    a2 = float(np.median(umis))

    # A3：中位 genes/cell
    a3 = float(np.median(genes))

    # B1：Debris rate（UMI < min_umi）
    b1 = float((umis < QC_CFG["min_umi"]).sum() / n_cells) if n_cells > 0 else np.nan

    # B2：通過 QC 的細胞數
    qc_pass = (
        (umis >= QC_CFG["min_umi"]) &
        (genes >= QC_CFG["min_genes"]) &
        (adata.obs["pct_mt"].values <= QC_CFG["max_pct_mt"])
    )
    b2 = int(qc_pass.sum())

    # C1：Artificial co-expression（4 對不可能基因對）
    c1 = compute_c1(adata)

    return {"a2": a2, "a3": a3, "b1": b1, "b2": b2, "c1": c1, "n_cells_raw": n_cells}


def compute_ned(mask: np.ndarray, adata: ad.AnnData, n_hvgs: int = 1000) -> float:
    """
    Neighbor Expression Divergence (NED)：相鄰細胞對之間的 Hellinger 距離均值。

    高 NED = 相鄰細胞表達差異大 = 邊界清晰（好）
    低 NED = 相鄰細胞表達相似 = 轉錄本混入（差，過度擴張或大細胞方法）

    - C1（artificial co-expression）對稀疏資料不敏感（僅偵測特定基因對共現）
    - NED 使用全部 HVG，對各種規模的轉錄本混入都靈敏
    - 對 CRC 異質性組織（上皮+免疫+基質）效果尤佳

    Returns: NED ∈ [0,1]，1 = 最好（邊界清晰）
    """
    from scipy.ndimage import grey_dilation

    if mask is None or adata is None or adata.n_obs < 5 or int(mask.max()) < 2:
        return np.nan

    # 表達矩陣 → float32 dense
    X = adata.X
    if sp.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float32)
    else:
        X = np.array(X, dtype=np.float32)

    # 選 HVGs（方差最高）
    if X.shape[1] > n_hvgs:
        gene_var = X.var(axis=0)
        hvg_idx = np.argpartition(gene_var, -n_hvgs)[-n_hvgs:]
        X = X[:, hvg_idx]

    # L1 歸一化（機率分布）
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-10)
    X_prob = (X / row_sums).astype(np.float32)

    # cell_id → adata row 索引映射（adata.obs["cell_id"] 對應 mask 中的 cell ID）
    obs_cell_ids = adata.obs["cell_id"].values
    cid_to_row = {int(c): i for i, c in enumerate(obs_cell_ids)}

    # 找相鄰細胞對（grey_dilation max → boundary where mask < dilated）
    struct = np.ones((3, 3), dtype=np.int32)
    dilated = grey_dilation(mask.astype(np.int32), footprint=struct)
    bnd = (mask > 0) & (dilated != mask)
    ci = mask[bnd].astype(np.int32)
    cj = dilated[bnd].astype(np.int32)
    valid = (cj > 0) & (cj != ci)
    ci, cj = ci[valid], cj[valid]

    if len(ci) == 0:
        return np.nan

    # 去重
    pairs_arr = np.unique(np.sort(np.stack([ci, cj], axis=1), axis=1), axis=0)

    # 過濾：兩者都在 adata 中
    known = set(cid_to_row.keys())
    mask_a = np.array([int(a) in known for a in pairs_arr[:, 0]])
    mask_b = np.array([int(b) in known for b in pairs_arr[:, 1]])
    pairs_arr = pairs_arr[mask_a & mask_b]

    if len(pairs_arr) < 5:
        return np.nan

    # 最多 3000 對（隨機取樣）
    if len(pairs_arr) > 3000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(pairs_arr), 3000, replace=False)
        pairs_arr = pairs_arr[idx]

    i_idx = np.array([cid_to_row[int(a)] for a in pairs_arr[:, 0]])
    j_idx = np.array([cid_to_row[int(b)] for b in pairs_arr[:, 1]])

    Xi = X_prob[i_idx]
    Xj = X_prob[j_idx]

    # Hellinger distance：√(Σ(√pi − √pj)² / 2)，範圍 [0,1]
    sqrt_i = np.sqrt(np.maximum(Xi, 0))
    sqrt_j = np.sqrt(np.maximum(Xj, 0))
    hell = np.sqrt(np.sum((sqrt_i - sqrt_j) ** 2, axis=1) / 2)

    return float(np.clip(np.mean(hell), 0, 1))


def compute_c1(adata: ad.AnnData) -> float:
    """
    計算 artificial co-expression rate：
    生物學不可能的基因對同時在同一細胞中表達（> 0）的細胞比例。
    對 4 個基因對取平均。

    返回 mean(co-expression rate over all pairs)；
    若某基因不在 panel 中則跳過該對。
    """
    import scipy.sparse as sp

    rates = []
    for gene_a, gene_b in IMPOSSIBLE_PAIRS:
        # 確認基因存在
        if gene_a not in adata.var_names or gene_b not in adata.var_names:
            continue

        idx_a = adata.var_names.get_loc(gene_a)
        idx_b = adata.var_names.get_loc(gene_b)

        X = adata.X
        if sp.issparse(X):
            col_a = np.asarray(X[:, idx_a].todense()).ravel()
            col_b = np.asarray(X[:, idx_b].todense()).ravel()
        else:
            col_a = X[:, idx_a]
            col_b = X[:, idx_b]

        # 兩者都 > 0 的細胞比例
        coexpr = ((col_a > 0) & (col_b > 0)).sum() / adata.n_obs
        rates.append(float(coexpr))

    return float(np.mean(rates)) if rates else np.nan


# ── 主流程 ────────────────────────────────────────────────────────────────

def run_metrics_ac():
    records = []

    for method in METHODS:
        print(f"\n[{method.upper()}] 計算 A/C 指標...")

        for roi_name in tqdm(ROI_INFO.keys(), desc=method):
            h5ad_path = ANNDATA_DIR / f"{method}_{roi_name}.h5ad"
            if not h5ad_path.exists():
                print(f"  ⚠️  {method}_{roi_name}.h5ad 不存在，跳過")
                continue

            # A1（從 attribution parquet 計算）
            a1 = compute_a1(method, roi_name)

            # A2, A3, C1, B1, B2（從 h5ad 計算）
            adata = sc.read_h5ad(h5ad_path)
            metrics = compute_metrics_from_adata(adata)

            # NED 和 UMI density（從 mask + h5ad 計算）
            mask_path = MASKS_DIR / f"{method}_{roi_name}.npy"
            if mask_path.exists():
                mask_arr = np.load(mask_path)
                ned = compute_ned(mask_arr, adata)
                umi_density = compute_umi_density(mask_arr, adata)
            else:
                ned = np.nan
                umi_density = np.nan

            records.append({
                "method":     method,
                "roi":        roi_name,
                "a1_capture": a1,
                "a1_umi_density":  umi_density,
                "a2_median_umi":   metrics["a2"],
                "a3_median_genes": metrics["a3"],
                "c1_coexpr":       metrics["c1"],
                "ned":             ned,
                "b1_debris_rate":  metrics["b1"],
                "b2_qc_yield":     metrics["b2"],
                "n_cells_raw":     metrics["n_cells_raw"],
            })

            print(
                f"  {roi_name}: FTC={a1:.3f} density={umi_density:.3f} A2={metrics['a2']:.0f} "
                f"A3={metrics['a3']:.0f} C1={metrics['c1']:.4f} NED={ned:.3f}"
            )

    df = pd.DataFrame(records)
    out = METRICS_DIR / "metrics_ac.csv"
    df.to_csv(out, index=False)
    print(f"\n✅ 03_metrics_ac.py 完成")
    print(f"   輸出：{out}")
    print(df.groupby("method")[["a1_capture","a2_median_umi","a3_median_genes","c1_coexpr","ned"]].mean().round(4))


if __name__ == "__main__":
    run_metrics_ac()
