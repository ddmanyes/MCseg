"""
預處理模組
提供 QC 過濾、正規化、高變異基因選擇等功能
"""

import logging
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

logger = logging.getLogger("pipeline.preprocessing")


class Preprocessor:
    """
    預處理器，統一處理所有資料集的 QC 與正規化
    """

    def __init__(self, config: dict[str, Any]):
        """
        初始化預處理器

        Parameters
        ----------
        config : dict
            配置字典，需包含 preprocessing 區塊
        """
        self.config = config
        self.params = config.get("preprocessing", {})

    def calculate_qc_metrics(self, adata: ad.AnnData, qc_params: dict = None) -> ad.AnnData:
        """
        計算 QC 指標
        """
        logger.info("計算 QC 指標...")

        params = qc_params if qc_params else self.params.get("cellular", {})

        # 標記粒線體基因
        mito_prefix = params.get("mito_prefix", "mt-")
        adata.var["mt"] = adata.var_names.str.lower().str.startswith(mito_prefix.lower())

        # 計算 QC 指標
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True
        )

        logger.info(f"  - 粒線體基因數: {adata.var['mt'].sum()}")
        logger.info(f"  - 平均 UMI/細胞: {adata.obs['total_counts'].mean():.1f}")
        logger.info(f"  - 平均基因/細胞: {adata.obs['n_genes_by_counts'].mean():.1f}")

        return adata

    def filter_cells(self, adata: ad.AnnData, qc_params: dict = None) -> ad.AnnData:
        """
        過濾低品質細胞
        """
        n_before = adata.n_obs

        params = qc_params if qc_params else self.params.get("cellular", {})

        min_genes = params.get("min_genes", 20)
        max_genes = params.get("max_genes", 8000)
        max_pct_mito = params.get("max_pct_mito", 20)
        min_counts = params.get("min_counts")

        logger.info("過濾細胞...")
        logger.info(f"  - 參數: min_genes={min_genes}, max_genes={max_genes}, max_pct_mito={max_pct_mito}, min_counts={min_counts}")
        logger.info(f"  - 過濾前: {n_before:,} 細胞")

        # 基因數過濾
        sc.pp.filter_cells(adata, min_genes=min_genes)

        # UMI 總數過濾
        if min_counts is not None:
            sc.pp.filter_cells(adata, min_counts=min_counts)

        # 最大基因數過濾 (排除 doublets)
        if max_genes:
            adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :].copy()

        # 粒線體過濾
        if "pct_counts_mt" in adata.obs.columns:
            adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito, :].copy()

        n_after = adata.n_obs
        logger.info(f"  - 過濾後: {n_after:,} 細胞 (移除 {n_before - n_after:,})")

        return adata

    def filter_genes(self, adata: ad.AnnData, qc_params: dict = None) -> ad.AnnData:
        """
        過濾低表達基因
        """
        n_before = adata.n_vars

        params = qc_params if qc_params else self.params.get("cellular", {})
        min_cells = params.get("min_cells", 3)

        logger.info("過濾基因...")
        sc.pp.filter_genes(adata, min_cells=min_cells)

        n_after = adata.n_vars
        logger.info(f"  - 過濾前: {n_before:,} 基因")
        logger.info(f"  - 過濾後: {n_after:,} 基因 (移除 {n_before - n_after:,})")

        return adata

    def normalize(self, adata: ad.AnnData) -> ad.AnnData:
        """
        正規化資料
        """
        params = self.params.get("normalization", {})
        target_sum = params.get("target_sum", 10000)

        logger.info("正規化...")

        # 保存原始計數
        adata.layers["counts"] = adata.X.copy()

        # 正規化
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)

        logger.info(f"  - 目標總數: {target_sum:,}")

        return adata

    def select_hvg(self, adata: ad.AnnData) -> ad.AnnData:
        """
        選擇高變異基因
        """
        params = self.params.get("hvg", {})
        n_top_genes = params.get("n_top_genes", 2000)
        flavor = params.get("flavor", "seurat_v3")

        logger.info("選擇高變異基因...")

        # seurat_v3 需要原始計數
        if flavor == "seurat_v3" and "counts" in adata.layers:
            # 檢查是否有足夠的基因
            if adata.n_vars < n_top_genes:
                logger.warning(f"基因數 ({adata.n_vars}) 少於 HVG 目標 ({n_top_genes})，跳過 HVG 篩選")
                adata.var["highly_variable"] = True
                return adata

            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor=flavor,
                layer="counts"
            )
        else:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor="seurat"
            )

        n_hvg = adata.var["highly_variable"].sum()
        logger.info(f"  - 選擇 {n_hvg:,} 個高變異基因")

        return adata

    def run_pca(self, adata: ad.AnnData) -> ad.AnnData:
        """
        執行 PCA 降維

        Parameters
        ----------
        adata : AnnData
            輸入資料

        Returns
        -------
        AnnData
            加入 PCA 結果的資料
        """
        n_pcs = self.config.get("clustering", {}).get("n_pcs", 50)

        logger.info("執行 PCA...")

        # 確保 n_pcs 不超過 min(n_obs, n_vars) - 1（svd 數學限制）
        max_pcs = min(adata.n_obs, adata.n_vars) - 1
        if n_pcs > max_pcs:
            logger.warning(
                f"n_pcs={n_pcs} 超過資料維度上限 {max_pcs}（{adata.n_obs} 細胞 × {adata.n_vars} 基因），"
                f"自動降至 {max_pcs}"
            )
            n_pcs = max(1, max_pcs)

        # OOM 防護: 針對大型資料集跳過密集化縮放 (Dense Scaling)
        if adata.n_obs > 100000:
            logger.warning(f"資料集過大 ({adata.n_obs:,} > 100k)，跳過 sc.pp.scale 以避免記憶體不足")
            sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
        else:
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, n_comps=n_pcs)

        logger.info(f"  - 主成分數: {n_pcs}")

        return adata

    def preprocess(
        self,
        adata: ad.AnnData,
        run_pca: bool = True,
        qc_key: str = "cellular"
    ) -> ad.AnnData:
        """
        執行完整預處理流程
        """
        dataset_name = adata.uns.get("dataset_name", "Unknown")
        logger.info("=" * 60)
        logger.info(f"預處理: {dataset_name} (Mode: {qc_key})")
        logger.info("=" * 60)

        # 獲取特定 QC 參數
        qc_params = self.params.get(qc_key, {})
        if not qc_params:
            logger.warning(f"找不到 '{qc_key}' 預處理設定，將使用預設 cellular 設定")
            qc_params = self.params.get("cellular", {})

        adata = self.calculate_qc_metrics(adata, qc_params)
        adata = self.filter_cells(adata, qc_params)
        adata = self.filter_genes(adata, qc_params)
        adata = self.normalize(adata)
        adata = self.select_hvg(adata)

        if run_pca:
            adata = self.run_pca(adata)

        logger.info("=" * 60)
        logger.info(f"預處理完成: {adata.n_obs:,} 細胞, {adata.n_vars:,} 基因")
        logger.info("=" * 60)

        return adata


def compute_sparsity(adata: ad.AnnData) -> float:
    """
    計算資料稀疏度 (零值比例)

    Parameters
    ----------
    adata : AnnData
        輸入資料

    Returns
    -------
    float
        稀疏度 (0-1)
    """
    if hasattr(adata.X, "toarray"):
        # 稀疏矩陣
        n_zeros = adata.X.shape[0] * adata.X.shape[1] - adata.X.nnz
    else:
        n_zeros = (adata.X == 0).sum()

    total = adata.X.shape[0] * adata.X.shape[1]
    return n_zeros / total
