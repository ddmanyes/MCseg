"""
04_metrics_de.py
================
計算 D 類（分群品質）與 E 類（免疫細胞存活）指標：

  D1: Silhouette score（per-method，合併所有 ROI 細胞 → Leiden + PCA）
      設計依據：同組織 14 ROI 細胞合併後數量充足（~14k），Silhouette 更穩定；
      per-ROI 版本因每 ROI 僅 50–500 顆細胞而雜訊過高（Friedman ns）。
      D1 廣播至各 ROI，在 TAS 中作為 method-level 固定偏移量。
  E1: Small immune cell survival rate（CD3E / MS4A1 / NKG7 陽性且通過 QC 的比例）

輸出：results/metrics/metrics_de.csv
"""

from __future__ import annotations

import sys
import json
import warnings
import numpy as np
import pandas as pd
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
TAS_CFG = cfg["tas"]
PROFILE = cfg["tissue_profiles"]["CRC"]

ANNDATA_DIR = ROOT / PATHS["anndata_dir"]
METRICS_DIR = ROOT / PATHS["metrics_dir"]
METRICS_DIR.mkdir(parents=True, exist_ok=True)

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

METHODS             = DATA["methods"]
IMMUNE_MARKERS      = PROFILE["small_immune_markers"]  # [CD3E, CD3D, MS4A1, NKG7, GNLY, PTPRC]
LEIDEN_RESOLUTION   = TAS_CFG["leiden_resolution"]     # 0.5
MIN_UMI             = QC_CFG["min_umi"]                # 50
E1_SPECIFICITY      = TAS_CFG.get("e1_marker_specificity", 0.01)  # marker_UMI/total_UMI 閾值


# ── D1：Silhouette（per-method，所有 ROI 合併）────────────────────────────

def compute_d1_silhouette_method(method: str) -> float:
    """
    合併同方法所有 ROI 的細胞後計算 Silhouette score。
    QC 過濾 → 合併 → seurat_v3 HVG → normalize → PCA → Leiden → Silhouette

    採用 global 計算的原因：
    - 14 ROI 來自同一組織，per-ROI 僅 50–500 顆細胞，Silhouette 雜訊極高
    - 合併後細胞數 ~14k，分群品質評估更穩定可靠
    """
    from sklearn.metrics import silhouette_score

    adatas = []
    for roi_name in ROI_INFO.keys():
        h5ad_path = ANNDATA_DIR / f"{method}_{roi_name}.h5ad"
        if not h5ad_path.exists():
            continue
        adata_roi = sc.read_h5ad(h5ad_path)
        # QC 過濾
        sc.pp.filter_cells(adata_roi, min_counts=MIN_UMI)
        sc.pp.filter_cells(adata_roi, min_genes=QC_CFG["min_genes"])
        if "pct_mt" in adata_roi.obs.columns:
            adata_roi = adata_roi[adata_roi.obs["pct_mt"] <= QC_CFG["max_pct_mt"]].copy()
        if adata_roi.n_obs > 0:
            adatas.append(adata_roi)

    if not adatas:
        warnings.warn(f"[{method}] 無有效 h5ad 檔案，D1 = nan")
        return np.nan

    # 合併所有 ROI
    adata = ad.concat(adatas, merge="same")

    if adata.n_obs < 50:
        warnings.warn(f"[{method}] 合併後細胞數不足（{adata.n_obs}），D1 = nan")
        return np.nan

    # HVG（seurat_v3 需要原始計數；細胞數少時 LOESS 可能奇異，fallback 至 seurat）
    n_hvg = min(2000, adata.n_vars)
    hvg_done = False
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3", span=0.3)
        hvg_done = True
    except ValueError:
        warnings.warn(f"[{method}] seurat_v3 HVG near-singular，fallback seurat")

    if hvg_done:
        adata = adata[:, adata.var["highly_variable"]].copy()

    # 標準化流程
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if not hvg_done:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    n_comps = min(30, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_comps)

    # Leiden clustering
    n_pcs = min(20, n_comps)
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, flavor="igraph",
                 directed=False, n_iterations=2)

    n_clusters = adata.obs["leiden"].nunique()
    if n_clusters < 2:
        warnings.warn(f"[{method}] Leiden 只產生 {n_clusters} 個 cluster，D1 = nan")
        return np.nan

    # Silhouette on PCA embedding
    pca_coords = adata.obsm["X_pca"][:, :n_pcs]
    labels = adata.obs["leiden"].astype(int).values
    d1 = float(silhouette_score(pca_coords, labels, sample_size=min(5000, adata.n_obs)))

    print(f"  [{method}] D1={d1:.4f} ({adata.n_obs} cells, {n_clusters} clusters)")
    return d1


# ── E1：小免疫細胞存活率（per-ROI，再取各 ROI 平均）────────────────────

def compute_e1_per_roi(adata: ad.AnnData) -> float:
    """
    E1 = (small_immune cells & UMI >= min_umi) / max(small_immune cells, 1)

    small_immune（v1.4）：
      - 擴展至 6 個 marker：CD3E, CD3D, MS4A1, NKG7, GNLY, PTPRC（方向 C）
      - 特異性過濾：marker_UMI / total_UMI >= E1_SPECIFICITY（方向 A）
        → 排除 ambient RNA 造成的假陽性（1 UMI 噪音）
    防零除：分母使用 max(..., 1)
    """
    import scipy.sparse as sp

    present_markers = [m for m in IMMUNE_MARKERS if m in adata.var_names]
    if not present_markers:
        return np.nan

    X = adata.X
    total_umi = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)

    # 累加所有 immune marker 的 UMI
    marker_umi = np.zeros(adata.n_obs, dtype=float)
    for marker in present_markers:
        idx = adata.var_names.get_loc(marker)
        col = X[:, idx]
        if sp.issparse(col):
            col = np.asarray(col.todense()).ravel()
        marker_umi += col.astype(float)

    # 特異性過濾：marker_UMI / total_UMI >= 閾值（過濾 ambient noise）
    specificity = marker_umi / np.maximum(total_umi, 1)
    small_immune = specificity >= E1_SPECIFICITY

    total_immune = small_immune.sum()
    if total_immune == 0:
        return np.nan

    # UMI 存活過濾
    umis = adata.obs["n_umis"].values if "n_umis" in adata.obs.columns else total_umi
    survived = (small_immune & (umis >= MIN_UMI)).sum()
    return float(survived / max(total_immune, 1))


# ── 主流程 ────────────────────────────────────────────────────────────────

def run_metrics_de():
    records = []

    print("\n計算 D1 Silhouette（per-method global）...")
    d1_per_method = {}
    for method in tqdm(METHODS, desc="D1"):
        d1_per_method[method] = compute_d1_silhouette_method(method)

    print("\n計算 E1 Small Immune Survival（per-ROI）...")
    for method in METHODS:
        for roi_name in ROI_INFO.keys():
            h5ad_path = ANNDATA_DIR / f"{method}_{roi_name}.h5ad"

            e1 = np.nan
            if h5ad_path.exists():
                adata = sc.read_h5ad(h5ad_path)
                e1 = compute_e1_per_roi(adata)

            records.append({
                "method":              method,
                "roi":                 roi_name,
                "d1_silhouette":       d1_per_method[method],
                "e1_immune_survival":  e1,
            })

    df_combined = pd.DataFrame(records)
    out = METRICS_DIR / "metrics_de.csv"
    df_combined.to_csv(out, index=False)

    print(f"\n✅ 04_metrics_de.py 完成")
    print(f"   輸出：{out}")
    print("\nD1 Silhouette (global per method):")
    print(df_combined.groupby("method")["d1_silhouette"].first())
    print("\nE1 Small Immune Survival (per-ROI mean) per method:")
    print(df_combined.groupby("method")["e1_immune_survival"].mean())


if __name__ == "__main__":
    run_metrics_de()
