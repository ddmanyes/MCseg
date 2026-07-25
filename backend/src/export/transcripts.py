"""從 Visium HD 2µm bins 生成 transcript-like CSV，供 Xenium Explorer 視覺化。

不依賴 FastAPI，可獨立測試。
"""
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger("pipeline.export.transcripts")


def generate_visiumhd_transcripts(
    adata_002um_path: Path,
    roi_cfg: dict,
    out_path: Path,
    pixel_size_um: float,
) -> "Optional[Path]":
    """
    從 Visium HD 2µm bin AnnData 生成 transcript-like CSV，供 Xenium Explorer 視覺化。

    資料來源：roi_out_dir/adata_002um.h5ad（bins 已裁切至 ROI 範圍）
    座標：ROI 局部 µm（原點 = ROI 左上角），由 obsm['spatial'] fullres px 換算。
    輸出：每個非零 (bin, gene) 對應一行 (x, y, gene)。

    Returns: out_path（寫出成功），或 None（失敗）。
    """
    import scanpy as sc
    import scipy.sparse as sp
    import numpy as np
    import pandas as pd

    logger.info(f"從 Visium HD 2µm bins 生成 transcripts：{adata_002um_path}")
    try:
        adata = sc.read_h5ad(str(adata_002um_path))

        # bin 空間位置：obsm['spatial'] = (n_bins, 2)，fullres px，col/x 在前
        spatial = adata.obsm["spatial"]
        roi_x0 = float(roi_cfg.get("x", 0))
        roi_y0 = float(roi_cfg.get("y", 0))

        # 轉換為 ROI 局部 µm
        x_local = (spatial[:, 0] - roi_x0) * pixel_size_um
        y_local = (spatial[:, 1] - roi_y0) * pixel_size_um

        # 取出稀疏矩陣的非零位置 (bin_idx, gene_idx)
        X = adata.X
        csr = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        rows, cols = csr.nonzero()

        if len(rows) == 0:
            logger.warning("2µm bin 矩陣無非零 entries，跳過 transcripts 層")
            return None

        gene_names = np.array(adata.var_names)
        df = pd.DataFrame({
            "x": x_local[rows],
            "y": y_local[rows],
            "gene": gene_names[cols],
        })

        # Xenium Explorer 載入 transcripts 層時無分頁，超過 ~500 萬行會卡很久。
        # 隨機取樣至上限，保留空間分佈的代表性。
        MAX_TX_ROWS = 5_000_000
        if len(df) > MAX_TX_ROWS:
            logger.warning(
                f"轉錄點共 {len(df):,} 行，超過上限 {MAX_TX_ROWS:,}，"
                f"隨機取樣（seed=42）以加快 Xenium Explorer 載入速度。"
            )
            df = df.sample(n=MAX_TX_ROWS, random_state=42).reset_index(drop=True)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(str(out_path), index=False)
        logger.info(f"已寫出 {len(df):,} 個 Visium HD 轉錄點至 {out_path}")

        # 診斷：確認轉錄點座標範圍與 ROI 一致
        roi_w_um = roi_cfg.get("width_px", 0) * pixel_size_um
        roi_h_um = roi_cfg.get("height_px", 0) * pixel_size_um
        x_out = df["x"].values; y_out = df["y"].values
        logger.info(
            f"轉錄點座標範圍 x=[{x_out.min():.1f}, {x_out.max():.1f}] µm，"
            f"y=[{y_out.min():.1f}, {y_out.max():.1f}] µm"
        )
        logger.info(f"ROI 物理尺寸 {roi_w_um:.1f} × {roi_h_um:.1f} µm")
        if x_out.max() > roi_w_um * 1.1 or y_out.max() > roi_h_um * 1.1:
            logger.warning(
                "⚠️ 轉錄點超出 ROI 範圍！可能是 roi_x0/roi_y0 未正確套用，"
                "請確認 adata_002um.h5ad 的 obsm['spatial'] 使用 global fullres px 座標"
            )
        return out_path

    except Exception as exc:
        logger.warning(f"Visium HD transcripts 生成失敗（繼續執行）：{exc}")
        return None
