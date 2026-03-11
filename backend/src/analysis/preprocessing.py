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


def _remove_zero_var_hvg(adata: ad.AnnData) -> None:
    """scale 後將 HVG 中方差為零的基因標記為非 HVG，避免 ARPACK near-singularity 錯誤。
    僅計算 HVG 子集的方差，不對全矩陣密集化，避免 OOM。
    """
    if "highly_variable" not in adata.var:
        return
    import scipy.sparse as sp
    hvg_idx = np.where(adata.var["highly_variable"].values)[0]
    if len(hvg_idx) == 0:
        return
    X_hvg = adata.X[:, hvg_idx]
    # 分批計算每欄方差，避免一次密集化整個 HVG 矩陣造成 OOM
    if sp.issparse(X_hvg):
        # 稀疏矩陣：E[X^2] - E[X]^2，逐欄計算
        mean_sq = np.asarray(X_hvg.power(2).mean(axis=0)).ravel()
        mean_   = np.asarray(X_hvg.mean(axis=0)).ravel()
        var_vals = mean_sq - mean_ ** 2
    else:
        var_vals = X_hvg.var(axis=0)
    zero_var_mask = var_vals == 0
    n_zero = int(zero_var_mask.sum())
    if n_zero:
        col_loc = adata.var.columns.get_loc("highly_variable")
        adata.var.iloc[hvg_idx[zero_var_mask], col_loc] = False
        logger.warning(f"移除 {n_zero} 個零方差 HVG（scale 後方差為 0），避免 PCA 奇異矩陣")


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

            try:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes,
                    flavor=flavor,
                    layer="counts"
                )
            except (ValueError, BaseException) as e:
                # seurat_v3 內部使用 LOESS 迴歸，當資料稀疏或平均 UMI 極低時會產生
                # 「near singularities」錯誤（bytes 格式）。Fallback 至 seurat flavor
                logger.warning(
                    f"seurat_v3 HVG 失敗（{e!r}），"
                    "資料可能過於稀疏（平均 UMI 低），改用 seurat flavor"
                )
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=min(n_top_genes, adata.n_vars),
                    flavor="seurat"
                )
        else:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=min(n_top_genes, adata.n_vars),
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
        else:
            sc.pp.scale(adata, max_value=10)

        # scale（或跳過 scale）後：清除殘餘 NaN/Inf（避免奇異矩陣）
        if hasattr(adata.X, "toarray"):
            pass  # sparse 理論上不含 NaN
        else:
            import numpy as _np
            if _np.isnan(adata.X).any() or _np.isinf(adata.X).any():
                adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=10.0, neginf=-10.0)
                logger.warning("scale 後偵測到 NaN/Inf，已替換為 0/±10")

        # scale（或跳過 scale）後移除零方差 HVG，避免 ARPACK near-singularity 錯誤
        _remove_zero_var_hvg(adata)

        # 重新計算 n_pcs 上限（HVG 數量可能因零方差過濾而減少）
        n_hvg = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else adata.n_vars
        n_pcs = min(n_pcs, min(adata.n_obs, n_hvg) - 1)
        n_pcs = max(1, n_pcs)

        try:
            sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
        except BaseException as arpack_err:
            # ARPACK 對近奇異矩陣可能失敗。
            # 注意：對稀疏資料 svd_solver='randomized' 會被 scanpy 忽略並再次使用 ARPACK，
            # 因此改用 zero_center=False 強制走 sklearn TruncatedSVD（隨機化算法，絕不使用 ARPACK）
            logger.warning(
                f"ARPACK PCA 失敗（{arpack_err!r}），改用 TruncatedSVD (zero_center=False)"
            )
            try:
                sc.tl.pca(adata, n_comps=n_pcs, zero_center=False)
            except BaseException as fallback_err:
                # 最後防線：進一步減少 n_pcs 再試
                safe_pcs = max(1, n_pcs // 2)
                logger.warning(
                    f"TruncatedSVD (n_pcs={n_pcs}) 失敗（{fallback_err!r}），"
                    f"縮減至 n_pcs={safe_pcs} 再試"
                )
                sc.tl.pca(adata, n_comps=safe_pcs, zero_center=False)

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
