"""
Stage 2.5：Proseg RNA 重分配

使用固定參數（max_dist=20, compactness=0.06, dilation=5）執行 Proseg MCMC，
沿用 Cellpose 分割遮罩作為細胞邊界，僅重新分配 RNA 至細胞。

核心流程：
  1. 讀取 adata_002um.h5ad 與 segmentation_masks.npy
  2. 將 2µm bins 展開為偽轉錄本 CSV
     （每個 bin 的每個基因計數展開為一行，位置 = bin 中心 µm）
  3. 對遮罩套用 dilation，作為 Proseg cell_id 初始分配
  4. 執行 proseg binary（MCMC RNA 重分配）
  5. 解析輸出（counts.csv.gz + cells.csv + genes.csv）
  6. 儲存 proseg_cells.h5ad（cells × genes）

輸出欄位：
  obs: centroid_x, centroid_y（ROI 局部 µm）+ Proseg 細胞 metadata
  obsm['spatial']: [centroid_x, centroid_y]（ROI 局部 µm）
  var: 同 Proseg 輸出（基因名稱）
"""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

import numpy as np

from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.proseg_rna")

# ── 固定預設參數 ─────────────────────────────────────────────────────
DEFAULT_MAX_DIST    = 20.0
DEFAULT_COMPACTNESS = 0.06
DEFAULT_DILATION    = 5      # px，對遮罩的膨脹半徑
DEFAULT_SAMPLES     = 500
DEFAULT_BURNIN      = 150
DEFAULT_RECORDED    = 150


# ── 工具函式 ─────────────────────────────────────────────────────────

def _get_proseg_bin(config: dict) -> str:
    """取得 Proseg 執行檔路徑（優先 config，次選 ~/.cargo/bin/proseg）。"""
    bin_path = config.get("paths", {}).get("proseg_bin", "") or "~/.cargo/bin/proseg"
    return str(Path(bin_path).expanduser())


def _get_spatial_coords(adata) -> np.ndarray:
    """從 adata 取得全域 fullres px 座標，回傳 (N, 2): [x_col, y_row]。"""
    if "spatial" in adata.obsm:
        return adata.obsm["spatial"]
    if "pxl_col_in_fullres" in adata.obs.columns and "pxl_row_in_fullres" in adata.obs.columns:
        return np.stack([
            adata.obs["pxl_col_in_fullres"].values.astype(float),
            adata.obs["pxl_row_in_fullres"].values.astype(float),
        ], axis=1)
    raise ValueError("adata_002um.h5ad 缺少空間座標（obsm['spatial'] 或 obs 欄位）")


def _dilate_mask(seg_mask: np.ndarray, radius: int) -> np.ndarray:
    """對分割遮罩套用正方形膨脹核（保留 cell_id）。"""
    import cv2
    if radius <= 0:
        return seg_mask
    ksize = 2 * radius + 1
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (ksize, ksize))
    return cv2.dilate(seg_mask.astype(np.float32), kernel).astype(np.int32)


# ── 核心函式 1：bins → 偽轉錄本 CSV ─────────────────────────────────

def bins_to_transcript_csv(
    adata_path: Path,
    mask_path: Path,
    roi_x_px: int,
    roi_y_px: int,
    out_csv: Path,
    pixel_size_um: float = VISIUM_UM_PX,
    dilation_px: int = DEFAULT_DILATION,
) -> int:
    """
    將 Visium HD 2µm bins 展開為 Proseg 格式的偽轉錄本 CSV。

    每個 bin 的每個基因計數（非零）展開為 count 條 CSV 行：
        x（ROI local µm）、y、gene、qv（=40）、cell_id（dilated mask）、z（=0）

    Returns
    -------
    int
        寫入的總轉錄本行數
    """
    import scanpy as sc
    import scipy.sparse as sp
    import pandas as pd

    logger.info(f"載入 adata：{adata_path}")
    adata = sc.read_h5ad(str(adata_path))
    n_bins, n_genes = adata.n_obs, adata.n_vars
    logger.info(f"  {n_bins:,} bins × {n_genes:,} genes")

    logger.info(f"載入分割遮罩：{mask_path}")
    seg_mask = np.load(str(mask_path))
    H, W = seg_mask.shape
    n_cells = int(seg_mask.max())
    logger.info(f"  遮罩 {W}×{H}，細胞數：{n_cells}")

    # 膨脹遮罩以擴大 cell_id 查找範圍
    logger.info(f"膨脹遮罩（radius={dilation_px}px）...")
    lookup_mask = _dilate_mask(seg_mask, dilation_px)

    # bin 全域 px → ROI 局部 px
    coords = _get_spatial_coords(adata)
    col_px = (coords[:, 0] - roi_x_px).round().astype(int)  # x
    row_px = (coords[:, 1] - roi_y_px).round().astype(int)  # y
    valid = (col_px >= 0) & (col_px < W) & (row_px >= 0) & (row_px < H)
    n_valid = int(valid.sum())
    n_out = n_bins - n_valid
    logger.info(f"  有效 bins（ROI 內）：{n_valid:,}，ROI 外：{n_out:,}")
    if n_bins > 0 and n_out / n_bins > 0.3:
        logger.warning(f"⚠️ ROI 外 bins 超 30%，請檢查 ROI 偏移設定（x={roi_x_px}, y={roi_y_px}）")

    # 轉換為稀疏 COO 格式
    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    X_coo = X.tocoo()
    gene_names = np.array(list(adata.var_names))

    # 向量化：篩選有效 COO 項目
    valid_nnz = valid[X_coo.row]               # bool: 每個非零 entry 的 bin 是否在 ROI 內
    bin_idx  = X_coo.row[valid_nnz]
    gene_idx = X_coo.col[valid_nnz]
    counts   = X_coo.data[valid_nnz].astype(int)

    # 取得每個 entry 的位置（ROI 局部 µm）與 cell_id
    x_um     = col_px[bin_idx].astype(float) * pixel_size_um
    y_um     = row_px[bin_idx].astype(float) * pixel_size_um
    cell_ids = lookup_mask[row_px[bin_idx], col_px[bin_idx]].astype(int)
    genes    = gene_names[gene_idx]

    # 依計數展開（每個 count 對應一行偽轉錄本）
    total_tx = int(counts.sum())
    logger.info(f"展開計數（共 {total_tx:,} 偽轉錄本）...")
    x_exp     = np.repeat(x_um,     counts)
    y_exp     = np.repeat(y_um,     counts)

    # 隨機 jitter：±half_bin µm，與 pipeline_2 Zarr builder 邏輯一致
    # 避免同一 bin 內所有偽轉錄本座標完全重疊，確保 Proseg MCMC 空間分配正常運作
    half_bin = pixel_size_um / 2.0
    rng = np.random.default_rng(42)
    x_exp = x_exp + rng.uniform(-half_bin, half_bin, size=total_tx)
    y_exp = y_exp + rng.uniform(-half_bin, half_bin, size=total_tx)
    logger.info(f"  套用 jitter ±{half_bin:.3f} µm（seed=42）")
    cid_exp   = np.repeat(cell_ids, counts)
    gene_exp  = np.repeat(genes,    counts)

    # 寫出 CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"寫出偽轉錄本 CSV：{out_csv}")
    df = pd.DataFrame({
        "x":       x_exp,
        "y":       y_exp,
        "gene":    gene_exp,
        "qv":      np.full(total_tx, 40, dtype=np.int8),
        "cell_id": cid_exp,
        "z":       np.zeros(total_tx, dtype=np.float32),
    })
    df.to_csv(str(out_csv), index=False)

    n_in_cell = int((cid_exp > 0).sum())
    logger.info(
        f"  寫出 {total_tx:,} 行；落在細胞內：{n_in_cell:,}（{n_in_cell/total_tx*100:.1f}%）"
    )
    return total_tx


# ── 核心函式 2：執行 Proseg binary ───────────────────────────────────

def run_proseg_binary(
    csv_path: Path,
    out_dir: Path,
    proseg_bin: str,
    coordinate_scale: float = VISIUM_UM_PX,
    max_dist: float = DEFAULT_MAX_DIST,
    compactness: float = DEFAULT_COMPACTNESS,
    samples: int = DEFAULT_SAMPLES,
    burnin: int = DEFAULT_BURNIN,
    recorded: int = DEFAULT_RECORDED,
) -> dict:
    """
    呼叫 proseg binary，回傳輸出檔案路徑字典。
    使用 Smart Resume：若輸出已存在則跳過執行。
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "polygons": out_dir / "proseg_results.json",
        "counts":   out_dir / "counts.csv.gz",
        "cells":    out_dir / "cells.csv",
        "genes":    out_dir / "genes.csv",
    }

    if all(p.exists() for p in outputs.values()):
        logger.info("偵測到 Proseg 輸出已存在，跳過執行（Smart Resume）")
        return outputs

    cmd = [
        proseg_bin,
        "--overwrite",
        "--output-cell-polygons",       str(outputs["polygons"]),
        "--output-counts",              str(outputs["counts"]),
        "--output-counts-fmt",          "csv-gz",
        "--output-cell-metadata",       str(outputs["cells"]),
        "--output-cell-metadata-fmt",   "csv",
        "--output-gene-metadata",       str(outputs["genes"]),
        "--output-gene-metadata-fmt",   "csv",
        "--coordinate-scale",           str(coordinate_scale),
        "--gene-column",                "gene",
        "--x-column",                   "x",
        "--y-column",                   "y",
        "--z-column",                   "z",
        "--cell-id-column",             "cell_id",
        "--cell-id-unassigned",         "0",
        "--ignore-z-coord",
        "--min-qv",                     "0",
        "--max-transcript-nucleus-distance", str(max_dist),
        "--cell-compactness",           str(compactness),
        "--samples",                    str(samples),
        "--burnin-samples",             str(burnin),
        "--recorded-samples",           str(recorded),
        "--nuclear-reassignment-prob",  "0",
        "--prior-seg-reassignment-prob", "0",
        "--enforce-connectivity",
        str(csv_path),
    ]
    logger.info("執行 Proseg：" + " ".join(cmd))
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Proseg 執行失敗（exit {result.returncode}）：{result.stderr[:1000]}")
        raise RuntimeError(
            f"Proseg 執行失敗（exit {result.returncode}）：{result.stderr[:500]}"
        )
    logger.debug(result.stdout)
    logger.info("Proseg 執行完成")

    for name, path in outputs.items():
        if not path.exists():
            raise FileNotFoundError(f"Proseg 輸出缺少：{name} ({path})")

    return outputs


# ── 核心函式 3：組建 AnnData ──────────────────────────────────────────

def assemble_proseg_anndata(
    proseg_outputs: dict,
    original_var,   # adata.var DataFrame（供對齊基因 metadata）
) -> "anndata.AnnData":
    """解析 Proseg 輸出（counts + cells + genes），組建 cells × genes AnnData。"""
    import anndata as ad
    import pandas as pd
    from scipy.io import mmread
    from scipy.sparse import csr_matrix
    import gzip

    logger.info("解析 Proseg 輸出...")

    cells_df = pd.read_csv(str(proseg_outputs["cells"]))
    genes_df = pd.read_csv(str(proseg_outputs["genes"]))
    logger.info(f"  cells: {len(cells_df):,}，genes: {len(genes_df):,}")

    # 讀取計數矩陣（可能是 gzip MTX）
    counts_path = proseg_outputs["counts"]
    with open(str(counts_path), "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(str(counts_path), "rt") as fh:
            X = mmread(fh).tocsr()
    else:
        X = mmread(str(counts_path)).tocsr()

    # 轉置（部分 Proseg 版本輸出 genes×cells）
    if X.shape[0] == len(genes_df) and X.shape[1] == len(cells_df):
        logger.info("  偵測到 genes×cells 格式，轉置為 cells×genes")
        X = X.T

    if X.shape != (len(cells_df), len(genes_df)):
        raise ValueError(
            f"計數矩陣形狀 {X.shape} 與 metadata 不符"
            f"（cells={len(cells_df)}, genes={len(genes_df)}）"
        )

    # 設定索引
    if "cell" in cells_df.columns:
        cells_df = cells_df.set_index("cell")
    elif "cell_id" in cells_df.columns:
        cells_df = cells_df.set_index("cell_id")

    if "gene" in genes_df.columns:
        genes_df = genes_df.set_index("gene")
    elif "gene_id" in genes_df.columns:
        genes_df = genes_df.set_index("gene_id")

    # 對齊 var（儘量保留原始 adata.var metadata）
    gene_index = genes_df.index.tolist()
    try:
        common = [g for g in gene_index if g in original_var.index]
        if len(common) == len(gene_index):
            var = original_var.loc[gene_index].copy()
        else:
            logger.info(f"  部分基因不在原始 var 中（{len(gene_index)-len(common)} 個），改用 Proseg genes 元資料")
            var = genes_df
    except Exception:
        var = genes_df

    result = ad.AnnData(X=csr_matrix(X), obs=cells_df, var=var)

    # obsm['spatial'] — 使用 Proseg 輸出的重心（ROI 局部 µm）
    for cx_col, cy_col in [("centroid_x", "centroid_y"), ("cx", "cy")]:
        if cx_col in cells_df.columns and cy_col in cells_df.columns:
            result.obsm["spatial"] = cells_df[[cx_col, cy_col]].values.astype(float)
            break

    logger.info(
        f"  完成：{result.n_obs:,} cells × {result.n_vars:,} genes，"
        f"中位 counts/cell = {int(np.median(np.array(X.sum(axis=1)).ravel()))}"
    )
    return result


# ── 多 ROI Pipeline（供 API 層呼叫）────────────────────────────────

def run_proseg_rna_pipeline(config: dict, roi_name: Optional[str] = None):
    """
    對指定 ROI（或全部 ROI）執行 Proseg RNA 重分配，儲存 proseg_cells.h5ad。

    Parameters
    ----------
    config   : pipeline.yaml 設定字典
    roi_name : 若指定，只跑單一 ROI；否則跑全部
    """
    from backend.src.utils.config import resolve_path

    paths    = config.get("paths", {})
    rois     = config.get("rois", [])
    out_base = resolve_path(paths.get("output_dir", "results/analysis")) / "roi"

    stage25      = config.get("proseg", {}).get("stage25", {})
    proseg_bin   = _get_proseg_bin(config)
    max_dist     = float(stage25.get("max_dist",          DEFAULT_MAX_DIST))
    compactness  = float(stage25.get("compactness",       DEFAULT_COMPACTNESS))
    dilation_px  = int(stage25.get("dilation",            DEFAULT_DILATION))
    samples      = int(stage25.get("samples",             DEFAULT_SAMPLES))
    burnin       = int(stage25.get("burnin_samples",      DEFAULT_BURNIN))
    recorded     = int(stage25.get("recorded_samples",    DEFAULT_RECORDED))

    if roi_name:
        rois = [r for r in rois if r.get("name") == roi_name]
        if not rois:
            raise ValueError(f"找不到 ROI：{roi_name}")

    if not rois:
        raise ValueError("pipeline.yaml 未設定任何 ROI")

    for roi in rois:
        rn = roi.get("name", "")
        if not rn:
            logger.warning("跳過無名稱的 ROI 設定")
            continue

        roi_dir    = out_base / rn
        adata_path = roi_dir / "adata_002um.h5ad"
        mask_path  = roi_dir / "segmentation_masks.npy"
        out_path   = roi_dir / "proseg_cells.h5ad"
        work_dir   = roi_dir / "_proseg_work"

        if not adata_path.exists():
            logger.warning(f"[{rn}] 找不到 adata_002um.h5ad，跳過（請先完成 Stage 0）")
            continue
        if not mask_path.exists():
            logger.warning(f"[{rn}] 找不到 segmentation_masks.npy，跳過（請先完成 Stage 1）")
            continue

        roi_x_px      = int(roi.get("x", 0))
        roi_y_px      = int(roi.get("y", 0))
        pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

        logger.info(f"[{rn}] Step 1/3：bins → 偽轉錄本 CSV（dilation={dilation_px}px）...")
        work_dir.mkdir(parents=True, exist_ok=True)
        csv_path = work_dir / "transcripts_for_proseg.csv"

        n_tx = bins_to_transcript_csv(
            adata_path=adata_path,
            mask_path=mask_path,
            roi_x_px=roi_x_px,
            roi_y_px=roi_y_px,
            out_csv=csv_path,
            pixel_size_um=pixel_size_um,
            dilation_px=dilation_px,
        )
        logger.info(f"[{rn}]   {n_tx:,} 偽轉錄本")

        logger.info(
            f"[{rn}] Step 2/3：Proseg MCMC"
            f"（max_dist={max_dist}, compactness={compactness}）..."
        )
        proseg_outputs = run_proseg_binary(
            csv_path=csv_path,
            out_dir=work_dir,
            proseg_bin=proseg_bin,
            coordinate_scale=pixel_size_um,
            max_dist=max_dist,
            compactness=compactness,
            samples=samples,
            burnin=burnin,
            recorded=recorded,
        )

        logger.info(f"[{rn}] Step 3/3：組建 proseg_cells.h5ad...")
        import scanpy as sc
        adata_orig   = sc.read_h5ad(str(adata_path))
        adata_cells  = assemble_proseg_anndata(proseg_outputs, adata_orig.var)

        roi_dir.mkdir(parents=True, exist_ok=True)
        adata_cells.write_h5ad(str(out_path))
        logger.info(f"[{rn}] 完成：{out_path}（{adata_cells.n_obs:,} cells）")

    logger.info("所有 ROI Proseg RNA 重分配完成")
