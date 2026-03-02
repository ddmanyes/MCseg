"""
Stage 4: 完整分析流水線
整合 QC → 正規化 → HVG → PCA → UMAP → Leiden → 標記基因
"""
import logging
from pathlib import Path
from typing import Any

import scanpy as sc
import anndata as ad

from backend.src.analysis.preprocessing import Preprocessor
from backend.src.analysis.clustering import Analyzer
from backend.src.utils.config import resolve_path

logger = logging.getLogger("pipeline.analysis")


def run_analysis_pipeline(config: dict[str, Any]) -> ad.AnnData:
    """
    執行完整分析流程。

    Parameters
    ----------
    config : pipeline.yaml 配置字典

    Returns
    -------
    AnnData 含有聚類結果
    """
    paths = config["paths"]
    analysis_cfg = config.get("analysis", {})

    # 確定輸入 h5ad（Proseg 輸出）
    proseg_dir = resolve_path(paths["proseg_dir"])
    input_h5ad = proseg_dir / "processed_proseg_cyto.h5ad"

    if not input_h5ad.exists():
        raise FileNotFoundError(f"找不到 Proseg 輸出：{input_h5ad}")

    logger.info(f"載入資料：{input_h5ad}")
    adata = sc.read_h5ad(str(input_h5ad))
    adata.uns["dataset_name"] = "proseg_cyto"
    logger.info(f"  {adata.n_obs:,} 細胞, {adata.n_vars:,} 基因")

    # 預處理
    preprocessor = Preprocessor(analysis_cfg)
    adata = preprocessor.preprocess(adata, run_pca=True, qc_key="cellular")

    # 聚類
    analyzer = Analyzer(analysis_cfg)
    adata = analyzer.run_clustering(adata)

    # 儲存結果
    output_dir = resolve_path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "clustered_final.h5ad"
    adata.write_h5ad(str(out_path))
    logger.info(f"已儲存：{out_path}")

    # 產生圖表
    fig_dir = resolve_path(paths["figure_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)
    _save_figures(adata, fig_dir, analysis_cfg)

    return adata


def _save_figures(adata: ad.AnnData, fig_dir: Path, analysis_cfg: dict) -> None:
    """產生並儲存標準圖表"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dpi = analysis_cfg.get("visualization", {}).get("figure_dpi", 300)

    # UMAP
    if "X_umap" in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.umap(adata, color="cluster", ax=axes[0], show=False, title="Leiden Clusters")
        if "total_counts" in adata.obs:
            sc.pl.umap(adata, color="total_counts", ax=axes[1], show=False, title="Total Counts")
        plt.tight_layout()
        fig.savefig(str(fig_dir / "umap.png"), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"已儲存：{fig_dir / 'umap.png'}")

    # Leiden 分布
    if "cluster" in adata.obs:
        fig, ax = plt.subplots(figsize=(8, 5))
        cluster_counts = adata.obs["cluster"].value_counts().sort_index()
        cluster_counts.plot(kind="bar", ax=ax)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Cell Count")
        ax.set_title("Cells per Cluster")
        plt.tight_layout()
        fig.savefig(str(fig_dir / "cluster_distribution.png"), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
