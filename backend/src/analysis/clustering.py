"""
分析模組
提供聚類、空間分析、比較分析等功能
"""

import logging
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics import silhouette_score

logger = logging.getLogger("pipeline.clustering")


class Analyzer:
    """
    分析器，執行聚類與空間分析
    """

    def __init__(self, config: dict[str, Any]):
        """
        初始化分析器

        Parameters
        ----------
        config : dict
            配置字典
        """
        self.config = config
        self.clustering_params = config.get("clustering", {})
        self.spatial_params = config.get("spatial", {})

    def run_clustering(self, adata: ad.AnnData, target_clusters: int = None) -> ad.AnnData:
        """
        執行聚類分析 (鄰域圖 + UMAP + Leiden/Louvain)

        Parameters
        ----------
        adata : AnnData
            預處理後的資料
        target_clusters : int, optional
            目標群落數量，若設定則自動調整 resolution

        Returns
        -------
        AnnData
            加入聚類結果的資料
        """
        n_neighbors = self.clustering_params.get("n_neighbors", 15)
        n_pcs = self.clustering_params.get("n_pcs", 50)
        resolution = self.clustering_params.get("resolution", 0.8)
        method = self.clustering_params.get("method", "leiden")

        dataset_name = adata.uns.get("dataset_name", "Unknown")
        logger.info(f"執行聚類分析: {dataset_name}")

        # 建立鄰域圖
        logger.info(f"  - 建立鄰域圖 (n_neighbors={n_neighbors})")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

        # UMAP
        logger.info("  - 計算 UMAP")
        min_dist = self.clustering_params.get("min_dist", 0.3)
        sc.tl.umap(adata, min_dist=min_dist)

        # 聚類與解析度優化
        if target_clusters:
            logger.info(f"  - 自動調整 resolution 以接近 {target_clusters} 個群落...")
            resolution = self._optimize_resolution(adata, target_clusters, method)

        logger.info(f"  - 執行 {method} 聚類 (resolution={resolution:.2f})")
        if method == "leiden":
            sc.tl.leiden(adata, resolution=resolution, key_added="cluster")
        else:
            sc.tl.louvain(adata, resolution=resolution, key_added="cluster")

        n_clusters = adata.obs["cluster"].nunique()
        logger.info(f"  - 識別 {n_clusters} 個群落")

        # === 修正 IORegistryError (uns['log1p']['base'] = None) ===
        # h5py 序列化 None 時可能導致 encoding_type='null' 錯誤
        if "log1p" in adata.uns and isinstance(adata.uns["log1p"], dict):
            if adata.uns["log1p"].get("base") is None:
                logger.info("  - 修正 uns['log1p']['base'] = None 以避免儲存錯誤")
                del adata.uns["log1p"]["base"]

        return adata

    def _optimize_resolution(self, adata: ad.AnnData, target_n: int, method: str) -> float:
        """
        二分搜尋法尋找最佳 resolution
        """
        best_res = 1.0
        low = 0.1
        high = 3.0
        best_diff = float("inf")

        for i in range(10): # 最多嘗試 10 次
            mid = (low + high) / 2

            # 使用暫存 key 避免覆蓋
            temp_key = f"cluster_temp_{i}"
            if method == "leiden":
                sc.tl.leiden(adata, resolution=mid, key_added=temp_key)
            else:
                sc.tl.louvain(adata, resolution=mid, key_added=temp_key)

            n_clusters = adata.obs[temp_key].nunique()
            diff = abs(n_clusters - target_n)

            # 更新最佳結果
            if diff < best_diff:
                best_diff = diff
                best_res = mid

            logger.info(f"    - Iter {i+1}: Res={mid:.3f}, Clusters={n_clusters}, Target={target_n}")

            if n_clusters == target_n:
                best_res = mid
                break

            if n_clusters < target_n:
                low = mid
            else:
                high = mid

        return best_res

    def compute_silhouette(self, adata: ad.AnnData) -> float:
        """
        計算 Silhouette Score

        Parameters
        ----------
        adata : AnnData
            含有聚類結果的資料

        Returns
        -------
        float
            Silhouette Score
        """
        if "X_pca" not in adata.obsm:
            logger.warning("缺少 PCA 結果，無法計算 Silhouette Score")
            return np.nan

        if "cluster" not in adata.obs:
            logger.warning("缺少聚類結果，無法計算 Silhouette Score")
            return np.nan

        # 使用 PCA 空間計算
        X = adata.obsm["X_pca"]
        labels = adata.obs["cluster"].astype(int).values

        score = silhouette_score(X, labels, sample_size=min(10000, len(labels)))
        logger.info(f"  - Silhouette Score: {score:.4f}")

        return score

    def find_markers(
        self,
        adata: ad.AnnData,
        groupby: str = "cluster",
        n_genes: int = 25,
    ) -> pd.DataFrame:
        """
        尋找差異表達基因

        Parameters
        ----------
        adata : AnnData
            含有聚類結果的資料
        groupby : str
            分組欄位
        n_genes : int
            每群取幾個基因

        Returns
        -------
        pd.DataFrame
            差異表達基因表
        """
        logger.info(f"尋找差異表達基因 (groupby={groupby})...")

        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            key_added="rank_genes"
        )

        # 轉換為 DataFrame
        result = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes")

        logger.info(f"  - 完成")

        return result


def compare_datasets(
    datasets: dict[str, ad.AnnData],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    比較多個資料集的 QC 指標 (含細胞大小)

    Parameters
    ----------
    datasets : dict
        資料集字典 {name: AnnData}
    config : dict
        配置字典

    Returns
    -------
    pd.DataFrame
        比較結果表格
    """
    logger.info("=" * 60)
    logger.info("資料集比較")
    logger.info("=" * 60)

    metrics = []

    for key, adata in datasets.items():
        # 計算總計數
        if hasattr(adata.X, "sum"):
            total_counts = float(adata.X.sum())
        else:
            total_counts = float(np.sum(adata.X))

        # 計算每細胞統計
        if hasattr(adata.X, "toarray"):
            row_sums = np.array(adata.X.sum(axis=1)).flatten()
        else:
            row_sums = np.array(adata.X.sum(axis=1)).flatten()

        # 計算基因數
        if hasattr(adata.X, "toarray"):
            genes_per_cell = np.array((adata.X > 0).sum(axis=1)).flatten()
        else:
            genes_per_cell = np.array((adata.X > 0).sum(axis=1)).flatten()

        # 計算稀疏度
        if hasattr(adata.X, "nnz"):
            sparsity = 1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])
        else:
            sparsity = (adata.X == 0).sum() / (adata.X.shape[0] * adata.X.shape[1])

        # 計算細胞/Bin 大小 (µm^2)
        sizes = np.array([])
        if "proseg" in key.lower():
            if "surface_area" in adata.obs:
                sizes = adata.obs["surface_area"].values
            elif "area" in adata.obs:
                sizes = adata.obs["area"].values
            else:
                logger.warning(f"{key}: 找不到 surface_area 或 area 欄位，無法計算大小")
        elif "bin_2um" in key.lower():
            sizes = np.full(adata.n_obs, 4.0)  # 2µm x 2µm
        elif "bin_8um" in key.lower():
            sizes = np.full(adata.n_obs, 64.0) # 8µm x 8µm

        # 計算大小統計量
        if len(sizes) > 0:
            size_mean = np.mean(sizes)
            size_median = np.median(sizes)
            size_q1 = np.percentile(sizes, 25)
            size_q3 = np.percentile(sizes, 75)
        else:
            size_mean = size_median = size_q1 = size_q3 = np.nan

        dataset_info = config.get("datasets", {}).get(key, {})
        name = dataset_info.get("name", key)

        metrics.append({
            "Dataset": name,
            "Key": key,
            "Cells/Bins": adata.n_obs,
            "Genes": adata.n_vars,
            "Total UMIs": total_counts,
            "Mean UMIs/Cell": row_sums.mean(),
            "Median UMIs/Cell": np.median(row_sums),
            "Mean Genes/Cell": genes_per_cell.mean(),
            "Median Genes/Cell": np.median(genes_per_cell),
            "Sparsity": sparsity,
            # 新增大小統計
            "Mean Size (µm²)": size_mean,
            "Median Size (µm²)": size_median,
            "Q1 Size (µm²)": size_q1,
            "Q3 Size (µm²)": size_q3,
        })

    df = pd.DataFrame(metrics)

    # 格式化顯示
    logger.info("\n" + df.to_string(index=False))

    return df


def statistical_test(
    adata1: ad.AnnData,
    adata2: ad.AnnData,
    metric: str = "total_counts",
) -> dict[str, float]:
    """
    對兩個資料集進行統計檢定

    Parameters
    ----------
    adata1 : AnnData
        資料集 1
    adata2 : AnnData
        資料集 2
    metric : str
        比較的指標欄位

    Returns
    -------
    dict
        統計結果
    """
    if metric not in adata1.obs or metric not in adata2.obs:
        return {"error": f"指標 {metric} 不存在"}

    values1 = adata1.obs[metric].values
    values2 = adata2.obs[metric].values

    # Mann-Whitney U test
    stat, pvalue = stats.mannwhitneyu(values1, values2, alternative="two-sided")

    # 效果量 (rank-biserial correlation)
    n1, n2 = len(values1), len(values2)
    effect_size = 1 - (2 * stat) / (n1 * n2)

    return {
        "statistic": stat,
        "p_value": pvalue,
        "effect_size": effect_size,
        "n1": n1,
        "n2": n2,
    }
