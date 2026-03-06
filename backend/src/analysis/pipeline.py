"""
Stage 4: 完整分析流水線
整合 QC → 正規化 → HVG → PCA → UMAP → Leiden → 標記基因
提供三步驟分段執行：run_qc_step / run_umap_step / run_heatmap_step
"""
import base64
import logging
from pathlib import Path
from typing import Any

import scanpy as sc
import anndata as ad

from backend.src.analysis.preprocessing import Preprocessor
from backend.src.analysis.clustering import Analyzer
from backend.src.utils.config import resolve_path

logger = logging.getLogger("pipeline.analysis")


# ─────────────────────────── helper ──────────────────────────────

def _encode_image(path: Path) -> str:
    """將圖片檔案 base64 編碼為字串。"""
    return base64.b64encode(path.read_bytes()).decode()


def _fix_log1p(adata: ad.AnnData) -> None:
    """修正 uns['log1p']['base']=None 導致 h5py 序列化錯誤。"""
    if "log1p" in adata.uns and isinstance(adata.uns["log1p"], dict):
        if adata.uns["log1p"].get("base") is None:
            adata.uns["log1p"].pop("base", None)


# ─────────────────────── helpers ──────────────────────────────────

def _generate_overlay_images(
    pre_spatial: "np.ndarray",
    pre_obs_names: list,
    post_obs_names_set: set,
    he_path: Path,
    fig_dir: Path,
    pixel_size_um: float,
) -> dict[str, str]:
    """在 H&E 底圖上疊加細胞重心，產生 QC 前後比較圖。

    Returns
    -------
    dict  {"pre_qc": base64_preview_png, "post_qc": base64_preview_png}
    HD 版本（300 DPI）同時存檔至 fig_dir，但不放入回傳值（供下載）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile
    import numpy as np

    figures: dict[str, str] = {}

    if not he_path.exists():
        logger.warning(f"找不到 H&E 影像：{he_path}，跳過疊圖生成")
        return figures

    he = tifffile.imread(str(he_path))
    H, W = he.shape[:2]

    # µm → HE 像素
    x_px = pre_spatial[:, 0] / pixel_size_um
    y_px = pre_spatial[:, 1] / pixel_size_um

    kept_mask    = np.array([n in post_obs_names_set for n in pre_obs_names])
    removed_mask = ~kept_mask
    n_total  = len(pre_obs_names)
    n_kept   = int(kept_mask.sum())
    n_removed = int(removed_mask.sum())

    # marker size：preview dpi=150 → s=12（≈4px radius），HD dpi=300 → s=3（≈2pt）
    for dpi, suffix, s in [(150, "", 12), (300, "_hd", 3)]:
        figsize = (W / dpi, H / dpi)

        # ── Pre-QC overlay ──
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(he)
        ax.scatter(x_px, y_px, s=s, c="cyan", alpha=0.6, linewidths=0,
                   label=f"All cells ({n_total:,})")
        ax.set_title(f"Pre-QC  ·  {n_total:,} cells", fontsize=9)
        ax.axis("off")
        ax.legend(loc="upper right", fontsize=7, markerscale=2, framealpha=0.5)
        path = fig_dir / f"overlay_pre_qc{suffix}.png"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        if not suffix:
            figures["pre_qc"] = _encode_image(path)
        logger.info(f"已儲存 {path.name}")

        # ── Post-QC overlay（kept=綠，removed=紅）──
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(he)
        if removed_mask.any():
            ax.scatter(x_px[removed_mask], y_px[removed_mask], s=s, c="red",
                       alpha=0.5, linewidths=0, label=f"Removed ({n_removed:,})")
        ax.scatter(x_px[kept_mask], y_px[kept_mask], s=s, c="lime",
                   alpha=0.7, linewidths=0, label=f"Kept ({n_kept:,})")
        ax.set_title(
            f"Post-QC  ·  kept {n_kept:,} / removed {n_removed:,}", fontsize=9)
        ax.axis("off")
        ax.legend(loc="upper right", fontsize=7, markerscale=2, framealpha=0.5)
        path = fig_dir / f"overlay_post_qc{suffix}.png"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        if not suffix:
            figures["post_qc"] = _encode_image(path)
        logger.info(f"已儲存 {path.name}")

    return figures


# ─────────────────────── Step 1: QC + PCA ─────────────────────────

def run_qc_step(config: dict[str, Any]) -> dict[str, str]:
    """
    Step 1：QC 前處理 + PCA。

    流程：載入 proseg_cells.h5ad → 計算 QC 指標 → 繪 violin / scatter →
    過濾細胞基因 → normalize → HVG → PCA → 繪 elbow → 儲存 qc_preprocessed.h5ad

    Returns
    -------
    dict  {chart_name: base64_png}  violin / scatter / elbow
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = config["paths"]
    analysis_cfg = config.get("analysis", {})

    out_base = resolve_path(paths["output_dir"]) / "roi"
    roi_name = config.get("rois", [{"name": "text"}])[0].get("name", "text")
    input_h5ad = out_base / roi_name / "proseg_cells.h5ad"
    if not input_h5ad.exists():
        raise FileNotFoundError(f"找不到 Proseg 輸出：{input_h5ad}")

    fig_dir = resolve_path(paths["figure_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)
    output_dir = resolve_path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"QC Step：載入 {input_h5ad}")
    adata = sc.read_h5ad(str(input_h5ad))
    adata.uns["dataset_name"] = "proseg_cells"
    logger.info(f"  {adata.n_obs:,} 細胞, {adata.n_vars:,} 基因")

    preprocessor = Preprocessor(analysis_cfg)
    qc_params = analysis_cfg.get("preprocessing", {}).get("cellular", {})

    # ── 1. QC 指標計算（過濾前，以顯示完整分布） ──
    adata = preprocessor.calculate_qc_metrics(adata, qc_params)

    figures: dict[str, str] = {}

    # ── 2. Violin（過濾前） ──
    qc_keys = [k for k in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if k in adata.obs.columns]
    if qc_keys:
        titles = {
            "total_counts": "Total UMI",
            "n_genes_by_counts": "Genes per Cell",
            "pct_counts_mt": "% Mitochondrial",
        }
        thresholds: dict[str, list[tuple]] = {
            "n_genes_by_counts": [
                (qc_params.get("min_genes", 20), "red", f"min={qc_params.get('min_genes', 20)}"),
                (qc_params.get("max_genes", 8000), "orange", f"max={qc_params.get('max_genes', 8000)}"),
            ],
            "pct_counts_mt": [
                (qc_params.get("max_pct_mito", 20), "red", f"max={qc_params.get('max_pct_mito', 20)}%"),
            ],
        }
        if qc_params.get("min_counts"):
            thresholds["total_counts"] = [
                (qc_params["min_counts"], "red", f"min={qc_params['min_counts']}"),
            ]

        fig, axes = plt.subplots(1, len(qc_keys), figsize=(5 * len(qc_keys), 4))
        if len(qc_keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, qc_keys):
            sc.pl.violin(adata, [key], ax=ax, show=False)
            ax.set_title(titles.get(key, key))
            for val, color, label in thresholds.get(key, []):
                ax.axhline(val, color=color, linestyle="--", linewidth=1.2, alpha=0.8, label=label)
            if thresholds.get(key):
                ax.legend(fontsize=8)
            # 鎖定 y 軸：下限=0，上限夾至 99th percentile × 1.5
            data_99 = float(adata.obs[key].quantile(0.99))
            _, y_top = ax.get_ylim()
            new_top = data_99 * 1.5 if y_top > data_99 * 3 else y_top
            ax.set_ylim(bottom=0, top=new_top)
        plt.tight_layout()
        violin_path = fig_dir / "qc_violin.png"
        fig.savefig(str(violin_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["violin"] = _encode_image(violin_path)
        logger.info("已儲存 qc_violin.png")

    # ── 3. Scatter（counts vs genes，以粒線體著色） ──
    if "total_counts" in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
        color_by = "pct_counts_mt" if "pct_counts_mt" in adata.obs.columns else None
        fig, ax = plt.subplots(figsize=(7, 5))
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color=color_by, ax=ax, show=False)
        ax.set_title("UMI vs Genes" + (" (colored by % Mito)" if color_by else ""))
        plt.tight_layout()
        scatter_path = fig_dir / "qc_scatter.png"
        fig.savefig(str(scatter_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["scatter"] = _encode_image(scatter_path)
        logger.info("已儲存 qc_scatter.png")

    # ── 4. 過濾 + normalize + HVG + PCA ──
    # 過濾前保存空間資訊，供疊圖比較使用
    _pre_obs_names = list(adata.obs_names)
    _pre_spatial   = adata.obsm["spatial"].copy() if "spatial" in adata.obsm else None

    adata = preprocessor.filter_cells(adata, qc_params)
    adata = preprocessor.filter_genes(adata, qc_params)

    # 細胞數過少（< 3）時無法繼續——提前拋出有意義的錯誤
    if adata.n_obs < 3:
        raise ValueError(
            f"QC 過濾後僅剩 {adata.n_obs} 顆細胞，無法繼續分析。\n"
            "建議降低 min_counts / min_genes 門檻，或檢查 Proseg 分析結果品質。\n"
            f"目前參數：min_genes={qc_params.get('min_genes')}, "
            f"min_counts={qc_params.get('min_counts')}, "
            f"max_pct_mito={qc_params.get('max_pct_mito')}"
        )
    adata = preprocessor.normalize(adata)
    adata = preprocessor.select_hvg(adata)
    adata = preprocessor.run_pca(adata)

    # ── 5. PCA Elbow（直接用 matplotlib，避免舊版 scanpy 不支援 ax 參數） ──
    if "X_pca" in adata.obsm and "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        variance_ratio = adata.uns["pca"]["variance_ratio"]
        n_shown = min(50, len(variance_ratio))
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(range(1, n_shown + 1), variance_ratio[:n_shown], "o-", markersize=4)
        ax.set_xlabel("PC")
        ax.set_ylabel("Variance Ratio")
        ax.set_title("PCA Variance Ratio (Elbow Plot)")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        elbow_path = fig_dir / "pca_elbow.png"
        fig.savefig(str(elbow_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["elbow"] = _encode_image(elbow_path)
        logger.info("已儲存 pca_elbow.png")

    # ── 6. 疊圖（H&E + 細胞重心，QC 前後比較） ──
    if _pre_spatial is not None:
        he_path = out_base / roi_name / "he_crop.tif"
        pixel_size_um = config.get("rois", [{}])[0].get("pixel_size_um", 0.2737)
        overlay_figs = _generate_overlay_images(
            _pre_spatial, _pre_obs_names, set(adata.obs_names),
            he_path, fig_dir, pixel_size_um,
        )
        figures.update(overlay_figs)

    # ── 7. 儲存前處理結果 ──
    _fix_log1p(adata)
    qc_h5ad = output_dir / "qc_preprocessed.h5ad"
    adata.write_h5ad(str(qc_h5ad))
    logger.info(f"QC Step 完成：{adata.n_obs:,} 細胞剩餘，已儲存 {qc_h5ad}")

    return figures


# ──────────────────── Step 2: UMAP 多解析度 ────────────────────────

def run_umap_step(
    config: dict[str, Any],
    resolutions: list[float],
    n_pcs: int = 30,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
) -> dict[str, str]:
    """
    Step 2：讀取前處理結果，計算 UMAP + Leiden（多解析度）。

    Returns
    -------
    dict  {str(resolution): base64_png, "grid": base64_png}
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    fig_dir = resolve_path(paths["figure_dir"])

    qc_h5ad = output_dir / "qc_preprocessed.h5ad"
    if not qc_h5ad.exists():
        raise FileNotFoundError(f"找不到前處理資料，請先執行 QC 步驟：{qc_h5ad}")

    logger.info(f"UMAP Step：載入 {qc_h5ad}")
    adata = sc.read_h5ad(str(qc_h5ad))

    # 確保 n_pcs 不超過已計算的 PC 數
    available_pcs = adata.obsm["X_pca"].shape[1] if "X_pca" in adata.obsm else n_pcs
    n_pcs_use = min(n_pcs, available_pcs)

    logger.info(f"建立 KNN 圖 (n_neighbors={n_neighbors}, n_pcs={n_pcs_use})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_use)

    logger.info(f"計算 UMAP (min_dist={min_dist})")
    sc.tl.umap(adata, min_dist=min_dist)

    resolutions = sorted(set(resolutions))
    figures: dict[str, str] = {}

    # ── Leiden 各解析度 + 個別圖 ──
    for res in resolutions:
        key = f"leiden_{res}"
        logger.info(f"  Leiden resolution={res}")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()

        fig_single, ax_single = plt.subplots(figsize=(7, 6))
        sc.pl.umap(adata, color=key, ax=ax_single, show=False,
                   title=f"Resolution = {res}  ({n_clusters} clusters)")
        plt.tight_layout()
        single_path = fig_dir / f"umap_res{res}.png"
        fig_single.savefig(str(single_path), dpi=150, bbox_inches="tight")
        plt.close(fig_single)
        figures[str(res)] = _encode_image(single_path)

    # ── 合併 Grid 圖 ──
    n = len(resolutions)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig_grid, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows), squeeze=False)
    for i, res in enumerate(resolutions):
        key = f"leiden_{res}"
        row, col = divmod(i, ncols)
        ax = axes[row][col]
        n_clusters = adata.obs[key].nunique()
        sc.pl.umap(adata, color=key, ax=ax, show=False,
                   title=f"Res={res} ({n_clusters} clusters)")
    for i in range(n, nrows * ncols):
        row, col = divmod(i, ncols)
        axes[row][col].set_visible(False)
    plt.tight_layout()
    grid_path = fig_dir / "umap_grid.png"
    fig_grid.savefig(str(grid_path), dpi=150, bbox_inches="tight")
    plt.close(fig_grid)
    figures["grid"] = _encode_image(grid_path)

    # ── 儲存含所有 Leiden 結果的 h5ad ──
    _fix_log1p(adata)
    umap_h5ad = output_dir / "umap_computed.h5ad"
    adata.write_h5ad(str(umap_h5ad))
    logger.info(f"UMAP Step 完成，已儲存 {umap_h5ad}")

    return figures


# ────────────────────── Step 3: Heatmap ──────────────────────────

def run_heatmap_step(
    config: dict[str, Any],
    resolution: float,
    n_top_genes: int = 20,
) -> dict[str, str]:
    """
    Step 3：針對指定解析度同時產生兩張 marker gene 圖表。

    - heatmap：seaborn.clustermap，顯示所有 HVGs，行（cluster）與列（基因）均有樹枝圖
    - dotplot：sc.pl.dotplot，每 cluster 取 n_top_genes 個 marker 基因，
              點大小 = 表達細胞比例，顏色 = 平均表達量

    Returns
    -------
    dict  {"heatmap": base64_png, "dotplot": base64_png}
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.sparse import issparse

    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    fig_dir = resolve_path(paths["figure_dir"])

    umap_h5ad = output_dir / "umap_computed.h5ad"
    if not umap_h5ad.exists():
        raise FileNotFoundError(f"找不到 UMAP 資料，請先執行 UMAP 步驟：{umap_h5ad}")

    logger.info(f"Heatmap Step：載入 {umap_h5ad}")
    adata = sc.read_h5ad(str(umap_h5ad))

    leiden_key = f"leiden_{resolution}"
    available = [c for c in adata.obs.columns if c.startswith("leiden_")]
    if leiden_key not in adata.obs.columns:
        raise ValueError(
            f"找不到 resolution={resolution} 的 Leiden 結果。可用：{available}"
        )

    clusters = adata.obs[leiden_key].cat.categories.tolist()

    # ── 1. Heatmap：所有 HVGs，seaborn clustermap，行列均有樹枝圖 ────
    logger.info("Heatmap：使用全部 HVGs 建立 seaborn clustermap")

    if "highly_variable" in adata.var.columns:
        gene_list = adata.var_names[adata.var["highly_variable"]].tolist()
    else:
        gene_list = adata.var_names.tolist()

    # 計算各 cluster 的平均表達量（n_clusters × n_genes）
    X_sub = adata[:, gene_list].X
    if issparse(X_sub):
        X_sub = X_sub.toarray()

    cluster_means: dict[str, Any] = {}
    for cluster in clusters:
        mask = (adata.obs[leiden_key] == cluster).values
        cluster_means[str(cluster)] = X_sub[mask].mean(axis=0)

    df_mean = pd.DataFrame(cluster_means, index=gene_list).T  # shape: (n_clusters, n_genes)

    # min-max scale per gene (column)
    col_min = df_mean.min(axis=0)
    col_max = df_mean.max(axis=0)
    df_scaled = (df_mean - col_min) / (col_max - col_min + 1e-8)

    n_genes_total = len(gene_list)
    n_clusters = len(clusters)
    hm_h = max(4, n_clusters * 0.5 + 2)
    hm_w = max(14, min(n_genes_total * 0.12 + 3, 60))  # 最寬 60 吋
    show_gene_labels = n_genes_total <= 150

    # 只有 ≥2 個 cluster / gene 才能做層次聚類；否則關閉對應樹枝圖避免報錯
    do_row_cluster = n_clusters >= 2
    do_col_cluster = n_genes_total >= 2

    logger.info(
        f"clustermap：{n_clusters} clusters × {n_genes_total} genes，"
        f"show_labels={show_gene_labels}，row_cluster={do_row_cluster}，col_cluster={do_col_cluster}"
    )

    plt.close("all")
    g = sns.clustermap(
        df_scaled,
        cmap="viridis",
        figsize=(hm_w, hm_h),
        yticklabels=True,              # cluster 標籤
        xticklabels=show_gene_labels,  # 基因標籤（數量少時才顯示）
        row_cluster=do_row_cluster,    # cluster 樹枝圖（左側 / 行）
        col_cluster=do_col_cluster,    # 基因樹枝圖（上方 / 列）
        linewidths=0,
        cbar_kws={"label": "Scaled mean expr.", "shrink": 0.5},
    )
    if not show_gene_labels:
        g.ax_heatmap.set_xlabel(f"Genes ({n_genes_total} HVGs)")
    g.ax_heatmap.set_ylabel("Cluster")
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)

    fig_path_heatmap = fig_dir / "heatmap.png"
    g.savefig(str(fig_path_heatmap), dpi=150, bbox_inches="tight")
    plt.close("all")
    logger.info(f"heatmap.png 已儲存 → {fig_path_heatmap}")

    # ── 2. Dotplot：每 cluster 取 n_top_genes 個 marker 基因 ─────────
    logger.info(f"Dotplot：計算 rank_genes_groups (groupby={leiden_key}, n_genes={n_top_genes})")
    sc.tl.rank_genes_groups(adata, groupby=leiden_key, method="wilcoxon", n_genes=n_top_genes)

    top_genes_per_cluster: dict[str, list[str]] = {}
    seen: set[str] = set()
    for cluster in clusters:
        cluster_genes = list(adata.uns["rank_genes_groups"]["names"][cluster][:n_top_genes])
        unique_genes = [g for g in cluster_genes if g not in seen]
        seen.update(unique_genes)
        if unique_genes:
            top_genes_per_cluster[f"C{cluster}"] = unique_genes

    total_dot_genes = sum(len(v) for v in top_genes_per_cluster.values())
    dp_w = max(12, total_dot_genes * 0.3)
    dp_h = max(5, n_clusters * 0.6 + 3)

    fig_path_dotplot = fig_dir / "dotplot.png"
    try:
        dp = sc.pl.dotplot(
            adata,
            var_names=top_genes_per_cluster,
            groupby=leiden_key,
            show=False,
            return_fig=True,
            standard_scale="var",
            dendrogram=True,
            figsize=(dp_w, dp_h),
        )
        dp.savefig(str(fig_path_dotplot), dpi=150, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        logger.warning(f"dotplot return_fig 失敗（{e}），改用 plt.savefig fallback")
        plt.close("all")
        sc.pl.dotplot(
            adata,
            var_names=top_genes_per_cluster,
            groupby=leiden_key,
            show=False,
            standard_scale="var",
            dendrogram=True,
        )
        plt.gcf().set_size_inches(dp_w, dp_h)
        plt.savefig(str(fig_path_dotplot), dpi=150, bbox_inches="tight")
        plt.close("all")

    logger.info(f"dotplot.png 已儲存 → {fig_path_dotplot}")

    return {
        "heatmap": _encode_image(fig_path_heatmap),
        "dotplot": _encode_image(fig_path_dotplot),
    }


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
    out_base = resolve_path(paths["output_dir"]) / "roi"
    roi_name = config.get("rois", [{"name": "text"}])[0].get("name", "text")
    input_h5ad = out_base / roi_name / "proseg_cells.h5ad"

    if not input_h5ad.exists():
        raise FileNotFoundError(f"找不到 Proseg 輸出：{input_h5ad}")

    logger.info(f"載入資料：{input_h5ad}")
    adata = sc.read_h5ad(str(input_h5ad))
    adata.uns["dataset_name"] = "proseg_cells"
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
