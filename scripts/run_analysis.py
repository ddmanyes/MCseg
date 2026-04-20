"""
Scanpy 分析一鍵執行器
QC 過濾 → 標準化 → PCA → UMAP → Leiden 分群
輸出 results/analysis/cellpose_cells_clustered.h5ad
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc


def _auto_resolution(n_cells: int) -> float:
    if n_cells < 5000:
        return 0.8
    if n_cells < 50000:
        return 0.5
    if n_cells < 500000:
        return 0.3
    return 0.2


def run_analysis(
    input_path:  Path | str,
    output_dir:  Path | str,
    min_genes:   int   = 150,
    max_pct_mt:  float = 12.0,
    min_counts:  int   = 80,
    resolution:  float | None = None,
    n_pcs:       int | None = None,
) -> sc.AnnData:
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_path)
    print(f"載入：{adata.n_obs} cells × {adata.n_vars} genes")

    # QC 指標
    mt_prefix = "MT-" if any(g.startswith("MT-") for g in adata.var_names) else "mt-"
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    before = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= min_genes) &
        (adata.obs["pct_counts_mt"]     <= max_pct_mt) &
        (adata.obs["total_counts"]      >= min_counts)
    ].copy()
    print(f"QC 過濾：{before} → {adata.n_obs} cells (移除 {before - adata.n_obs})")

    if adata.n_obs < 10:
        raise RuntimeError(f"QC 後細胞數不足 ({adata.n_obs})，請放寬閾值")

    # 標準化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_hvg = min(2000, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, svd_solver="arpack")
    var_ratio = adata.uns["pca"]["variance_ratio"]
    cumvar    = np.cumsum(var_ratio)
    auto_pcs  = int(np.argmax(cumvar >= 0.85)) + 1
    used_pcs  = n_pcs if n_pcs else max(10, min(auto_pcs, 50))
    print(f"PCA：使用 {used_pcs} PCs（累積變異 {cumvar[used_pcs-1]:.1%}）")

    # UMAP + Leiden
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=used_pcs)
    sc.tl.umap(adata, min_dist=0.5)

    used_res = resolution if resolution is not None else _auto_resolution(adata.n_obs)
    sc.tl.leiden(adata, resolution=used_res, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    print(f"Leiden resolution={used_res}：{n_clusters} clusters")

    out_path = output_dir / "cellpose_cells_clustered.h5ad"
    adata.write(out_path)
    print(f"[OK] 已儲存：{out_path}")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description="MSseg Scanpy 分析流程")
    parser.add_argument("--input",       required=True,        help="輸入 h5ad 路徑")
    parser.add_argument("--output-dir",  required=True,        help="輸出目錄")
    parser.add_argument("--min-genes",   type=int,   default=150)
    parser.add_argument("--max-pct-mt",  type=float, default=12.0)
    parser.add_argument("--min-counts",  type=int,   default=80)
    parser.add_argument("--resolution",  type=float, default=None)
    parser.add_argument("--n-pcs",       type=int,   default=None)
    args = parser.parse_args()

    run_analysis(
        input_path = args.input,
        output_dir = args.output_dir,
        min_genes  = args.min_genes,
        max_pct_mt = args.max_pct_mt,
        min_counts = args.min_counts,
        resolution = args.resolution,
        n_pcs      = args.n_pcs,
    )


if __name__ == "__main__":
    main()
