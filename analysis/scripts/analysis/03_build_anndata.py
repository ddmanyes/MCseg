"""
02_build_anndata.py
===================
將 Visium HD 稀疏矩陣（cell×gene）按 cell_id 聚合：
  barcode × gene → cell_id × gene (group-by sum)

輸入：
  - filtered_feature_bc_matrix.h5（全片）
  - attribution/{method}_{roi}.parquet（bin → cell_id）

輸出：results/anndata/{method}_{roi}.h5ad
  - obs: cell_id, n_umis, n_genes, pct_mt
  - var: gene_ids, gene_names
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
PATHS = cfg["paths"]
DATA  = cfg["data"]
QC    = cfg["qc"]

ATTRIBUTION_DIR = ROOT / PATHS["attribution_dir"]
ANNDATA_DIR     = ROOT / PATHS["anndata_dir"]
ANNDATA_DIR.mkdir(parents=True, exist_ok=True)

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

METHODS = DATA["methods"]


def load_visiumhd_matrix() -> ad.AnnData:
    """讀取全片 filtered_feature_bc_matrix.h5（backed=False，全載入記憶體）。"""
    h5_path = Path(PATHS["visiumhd_matrix"])
    print(f"讀取表達矩陣: {h5_path.name}")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    print(f"  {adata.n_obs:,} barcodes × {adata.n_vars:,} genes")
    return adata


def groupby_sum(adata_roi: ad.AnnData) -> ad.AnnData:
    """
    按 obs["cell_id"] 分組加總表達量（向量化稀疏矩陣運算）。
    只保留 cell_id > 0 的 bins。
    """
    obs = adata_roi.obs.copy()
    cell_ids = obs["cell_id"].values.astype(np.int32)

    # 只保留有歸屬的 bins
    valid = cell_ids > 0
    if valid.sum() == 0:
        return None

    adata_valid = adata_roi[valid]
    cell_ids_valid = cell_ids[valid]

    unique_cells = np.unique(cell_ids_valid)
    n_cells = len(unique_cells)
    cell_id_to_idx = {cid: i for i, cid in enumerate(unique_cells)}

    # 建立聚合矩陣 A：shape (n_cells, n_bins_valid)，每行對應一顆細胞
    row_indices = np.array([cell_id_to_idx[c] for c in cell_ids_valid])
    col_indices = np.arange(len(cell_ids_valid))
    data_ones   = np.ones(len(cell_ids_valid), dtype=np.float32)
    A = sp.csr_matrix(
        (data_ones, (row_indices, col_indices)),
        shape=(n_cells, adata_valid.n_obs)
    )

    # 聚合表達量：(n_cells, n_genes) = A @ X
    X_valid = adata_valid.X
    if not sp.issparse(X_valid):
        X_valid = sp.csr_matrix(X_valid)
    X_cell = A @ X_valid   # → (n_cells, n_genes)

    # 建立 cell-level AnnData
    cell_adata = ad.AnnData(
        X   = X_cell.astype(np.float32),
        var = adata_roi.var.copy(),
        obs = pd.DataFrame(
            {"cell_id": unique_cells},
            index=[f"cell_{c}" for c in unique_cells]
        )
    )

    # 計算基本 QC
    sc.pp.calculate_qc_metrics(cell_adata, percent_top=None, log1p=False, inplace=True)
    mt_genes = cell_adata.var_names.str.upper().str.startswith("MT-")
    if mt_genes.sum() > 0:
        mt_expr = np.asarray(cell_adata[:, mt_genes].X.sum(axis=1)).ravel()
        total   = np.asarray(cell_adata.X.sum(axis=1)).ravel()
        cell_adata.obs["pct_mt"] = np.where(total > 0, mt_expr / total * 100, 0.0)
    else:
        cell_adata.obs["pct_mt"] = 0.0

    # 重命名 QC 欄位
    cell_adata.obs.rename(columns={
        "total_counts": "n_umis",
        "n_genes_by_counts": "n_genes"
    }, inplace=True)

    return cell_adata


def run_build():
    print("載入全片表達矩陣...")
    adata_full = load_visiumhd_matrix()

    # 建立 barcode → row index 的快速查找
    bc_to_idx = {bc: i for i, bc in enumerate(adata_full.obs_names)}
    print(f"  建立 barcode 索引完成")

    for method in METHODS:
        print(f"\n[{method.upper()}] 建立 AnnData...")

        for roi_name in tqdm(ROI_INFO.keys(), desc=method):
            dst = ANNDATA_DIR / f"{method}_{roi_name}.h5ad"
            if dst.exists():
                print(f"  {roi_name} 已存在，跳過")
                continue

            attr_path = ATTRIBUTION_DIR / f"{method}_{roi_name}.parquet"
            if not attr_path.exists():
                print(f"  ⚠️  attribution 不存在：{attr_path.name}，跳過")
                continue

            # 讀取歸屬結果
            attr = pd.read_parquet(attr_path)

            # 找出對應的 row indices（barcode 篩選）
            valid_bc = [bc for bc in attr["barcode"] if bc in bc_to_idx]
            if not valid_bc:
                print(f"  ⚠️  {roi_name} 無有效 barcode")
                continue

            roi_indices = [bc_to_idx[bc] for bc in valid_bc]
            adata_roi = adata_full[roi_indices].copy()

            # 將 cell_id 資訊加入 obs
            attr_indexed = attr.set_index("barcode")
            adata_roi.obs["cell_id"] = attr_indexed.loc[
                adata_roi.obs_names, "cell_id"
            ].values.astype(np.int32)

            # Group-by sum 聚合
            cell_adata = groupby_sum(adata_roi)
            if cell_adata is None:
                print(f"  ⚠️  {roi_name} 無有效細胞，跳過")
                continue

            cell_adata.uns["method"] = method
            cell_adata.uns["roi"]    = roi_name
            cell_adata.write_h5ad(dst)

            print(
                f"  {roi_name}: {cell_adata.n_obs} 顆細胞, "
                f"中位 UMI={np.median(cell_adata.obs['n_umis']):.0f}, "
                f"中位 genes={np.median(cell_adata.obs['n_genes']):.0f}"
            )

    print("\n✅ 02_build_anndata.py 完成")
    print(f"   輸出：{ANNDATA_DIR}")


if __name__ == "__main__":
    run_build()
