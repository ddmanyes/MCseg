"""
Stage 0: ROI 裁切模組
提供 H&E BTF/TIFF、Visium HD AnnData、Xenium 資料的 ROI 裁切功能

來源參考：
- xenium_visiumhd_comparison/scripts/02_baseline/cluster_8um.py
- xenium_visiumhd_comparison/scripts/03_pipeline/segment_cellpose.py
- xenium_visiumhd_comparison/scripts/03_pipeline/segment_proseg_zarr.py
- xenium_visiumhd_comparison/scripts/04_xenium/cluster_xenium.py
"""
from __future__ import annotations

import base64
import io
import json
import logging
import pathlib
from pathlib import Path
from typing import Callable, Optional

import cv2
import numpy as np

from backend.src.utils.constants import XENIUM_UM_PX, VISIUM_UM_PX

logger = logging.getLogger("pipeline.roi")


# ── 座標工具 ─────────────────────────────────────────────────

def roi_fullres_to_um(roi: dict) -> tuple[float, float, float, float]:
    """
    將 fullres pixel ROI 轉為 µm 座標。
    支援格式 A（fullres px）和格式 B（µm 直接）。

    Returns: (x0_um, y0_um, w_um, h_um)
    """
    if "x_um" in roi:
        return roi["x_um"], roi["y_um"], roi["width_um"], roi["height_um"]
    pxum = roi.get("pixel_size_um", VISIUM_UM_PX)
    return (
        roi["x"] * pxum,
        roi["y"] * pxum,
        roi["width_px"] * pxum,
        roi["height_px"] * pxum,
    )


def roi_to_fullres_px(roi: dict) -> tuple[int, int, int, int]:
    """
    將 ROI 統一轉為 fullres pixel 座標。

    Returns: (x0, y0, w, h)
    """
    if "x" in roi:
        return roi["x"], roi["y"], roi["width_px"], roi["height_px"]
    pxum = roi.get("pixel_size_um", VISIUM_UM_PX)
    x0 = round(roi["x_um"] / pxum)
    y0 = round(roi["y_um"] / pxum)
    w  = round(roi["width_um"] / pxum)
    h  = round(roi["height_um"] / pxum)
    return x0, y0, w, h


# ── H&E 影像裁切（BTF/TIFF tile-based）───────────────────────

def read_btf_crop(
    btf_path: str | Path,
    x0: int,
    y0: int,
    w: int,
    h: int,
    margin: int = 0,
) -> tuple[np.ndarray, int, int]:
    """
    Tile-based 讀取 BTF/TIFF ROI crop（嚴禁全圖載入）。

    Parameters
    ----------
    btf_path : path to BTF/TIFF file
    x0, y0, w, h : ROI in fullres pixel coordinates
    margin : extra pixels around ROI

    Returns
    -------
    (crop_rgb, actual_x0, actual_y0)
    """
    import tifffile

    with tifffile.TiffFile(str(btf_path)) as tf:
        page = tf.pages[0]
        img_h, img_w = page.imagelength, page.imagewidth
        TW = getattr(page, "tilewidth",  512)
        TH = getattr(page, "tilelength", 512)
        n_tiles_x = (img_w + TW - 1) // TW

        # 計算帶 margin 的 ROI 邊界
        fx0 = max(0, x0 - margin)
        fy0 = max(0, y0 - margin)
        fx1 = min(img_w, x0 + w + margin)
        fy1 = min(img_h, y0 + h + margin)

        tx0 = fx0 // TW;  tx1 = (fx1 + TW - 1) // TW
        ty0 = fy0 // TH;  ty1 = (fy1 + TH - 1) // TH

        offsets    = page.tags.get("TileOffsets", None)
        bytecounts = page.tags.get("TileByteCounts", None)

        canvas = np.zeros(((ty1 - ty0) * TH, (tx1 - tx0) * TW, 3), dtype=np.uint8)

        if offsets and bytecounts:
            # Uncompressed tiled TIFF（BTF）
            offsets_v    = offsets.value
            bytecounts_v = bytecounts.value
            with open(str(btf_path), "rb") as fh:
                for ty in range(ty0, ty1):
                    for tx in range(tx0, tx1):
                        tidx = ty * n_tiles_x + tx
                        if tidx >= len(offsets_v):
                            continue
                        fh.seek(offsets_v[tidx])
                        raw = np.frombuffer(fh.read(bytecounts_v[tidx]), np.uint8)
                        if raw.size == TH * TW * 3:
                            canvas[
                                (ty - ty0) * TH:(ty - ty0 + 1) * TH,
                                (tx - tx0) * TW:(tx - tx0 + 1) * TW,
                            ] = raw.reshape(TH, TW, 3)
        else:
            # 普通 TIFF 無 TileOffsets tag，無法做 tile-based 讀取
            # 全圖載入（page.asarray()）可能導致 OOM（10-80 GB 檔案）
            raise NotImplementedError(
                f"TIFF 檔 '{btf_path}' 沒有 TileOffsets tag，無法進行 tile-based 讀取。\n"
                "請將影像轉換為 BigTIFF tiled 格式：\n"
                "  tifffile.imwrite(out, img, bigtiff=True, tile=(512, 512), compression='zlib')"
            )

        # 精確裁切
        cx0 = fx0 - tx0 * TW
        cy0 = fy0 - ty0 * TH
        crop = canvas[cy0: cy0 + (fy1 - fy0), cx0: cx0 + (fx1 - fx0)]

    return crop, fx0, fy0


# ── Visium HD AnnData 裁切 ──────────────────────────────────

def subset_anndata_roi(adata, roi: dict, binned_dir: str | Path = None):
    """
    根據 ROI fullres pixel 座標，裁切 Visium HD AnnData。

    期望 adata.obs 含有：
    - pxl_col_in_fullres (X)
    - pxl_row_in_fullres (Y)
    """
    x0, y0, w, h = roi_to_fullres_px(roi)
    x1, y1 = x0 + w, y0 + h

    if binned_dir is not None:
        import json
        import pathlib
        import numpy as np
        reg_path = pathlib.Path(binned_dir) / "spatial" / "H1-WGR3TC4-D1-fiducials-image-registration.json"
        if reg_path.exists():
            with open(reg_path) as f:
                data = json.load(f)
            if "cytAssistInfo" in data and "transformImages" in data["cytAssistInfo"]:
                try:
                    T1 = np.array(data["cytAssistInfo"]["transformImages"])
                    T2 = np.array(data["transform"])
                    T2_inv = np.linalg.inv(T2)
                    M = T2_inv @ T1
                    
                    pts = np.array([
                        [x0, y0, 1], [x1, y0, 1],
                        [x0, y1, 1], [x1, y1, 1]
                    ]).T
                    res = M @ pts
                    xs = res[0, :] / res[2, :]
                    ys = res[1, :] / res[2, :]
                    mapped_x0, mapped_x1 = min(xs), max(xs)
                    mapped_y0, mapped_y1 = min(ys), max(ys)
                    logger.info(f"CytAssist 座標自動修正: 映射至內部 X={mapped_x0:.1f}~{mapped_x1:.1f}, Y={mapped_y0:.1f}~{mapped_y1:.1f}")
                    
                    # Compute Raw TIFF coordinates for all bins
                    M_inv = np.linalg.inv(M)
                    all_cols = adata.obs["pxl_col_in_fullres"].values
                    all_rows = adata.obs["pxl_row_in_fullres"].values
                    all_pts = np.vstack([all_cols, all_rows, np.ones_like(all_cols)])
                    raw_res = M_inv @ all_pts
                    adata.obs["raw_tiff_col"] = raw_res[0, :] / raw_res[2, :]
                    adata.obs["raw_tiff_row"] = raw_res[1, :] / raw_res[2, :]
                    
                    # Update global obsm['spatial'] to align with Raw TIFF
                    adata.obsm["spatial"] = np.stack([
                        adata.obs["raw_tiff_col"].values,
                        adata.obs["raw_tiff_row"].values,
                    ], axis=1)
                    
                    # Update mask to use the mapped coordinates for filtering
                    mask = (
                        (adata.obs["pxl_col_in_fullres"] >= mapped_x0) &
                        (adata.obs["pxl_col_in_fullres"] <  mapped_x1) &
                        (adata.obs["pxl_row_in_fullres"] >= mapped_y0) &
                        (adata.obs["pxl_row_in_fullres"] <  mapped_y1)
                    )
                except Exception as e:
                    logger.error(f"CytAssist 矩陣轉換失敗: {e}")
                    mask = (
                        (adata.obs["pxl_col_in_fullres"] >= x0) &
                        (adata.obs["pxl_col_in_fullres"] <  x1) &
                        (adata.obs["pxl_row_in_fullres"] >= y0) &
                        (adata.obs["pxl_row_in_fullres"] <  y1)
                    )
        else:
            mask = (
                (adata.obs["pxl_col_in_fullres"] >= x0) &
                (adata.obs["pxl_col_in_fullres"] <  x1) &
                (adata.obs["pxl_row_in_fullres"] >= y0) &
                (adata.obs["pxl_row_in_fullres"] <  y1)
            )
    else:
        mask = (
            (adata.obs["pxl_col_in_fullres"] >= x0) &
            (adata.obs["pxl_col_in_fullres"] <  x1) &
            (adata.obs["pxl_row_in_fullres"] >= y0) &
            (adata.obs["pxl_row_in_fullres"] <  y1)
        )
    sub = adata[mask].copy()
    logger.info(f"AnnData ROI 裁切：{mask.sum()} / {len(adata)} bins")
    return sub


def load_visium_adata(binned_dir: str | Path, bin_size: str = "002"):
    """
    讀取 Visium HD binned AnnData 並加入空間座標。

    Parameters
    ----------
    binned_dir : square_002um 或 square_008um 目錄
    """
    import scanpy as sc
    import pandas as pd

    binned_dir = Path(binned_dir)
    h5_path = binned_dir / "filtered_feature_bc_matrix.h5"
    pos_path = binned_dir / "spatial" / "tissue_positions.parquet"

    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    if pos_path.exists():
        pos_df = pd.read_parquet(str(pos_path)).set_index("barcode")
        pos_df = pos_df[pos_df["in_tissue"] == 1]
        common = adata.obs_names.intersection(pos_df.index)
        adata = adata[common].copy()
        adata.obs["pxl_col_in_fullres"] = pos_df.loc[common, "pxl_col_in_fullres"].values
        adata.obs["pxl_row_in_fullres"] = pos_df.loc[common, "pxl_row_in_fullres"].values
        # 同步設定 obsm['spatial']：[x=col, y=row]（fullres px）
        adata.obsm["spatial"] = np.stack([
            adata.obs["pxl_col_in_fullres"].values.astype(float),
            adata.obs["pxl_row_in_fullres"].values.astype(float),
        ], axis=1)

    logger.info(f"載入 Visium {bin_size}µm: {adata.n_obs:,} bins")
    return adata


# ── Xenium 資料裁切 ─────────────────────────────────────────

def load_xenium_adata(xenium_outs: str | Path):
    """載入 Xenium cell_feature_matrix + 細胞座標"""
    import scanpy as sc
    import pandas as pd

    xenium_outs = Path(xenium_outs)
    h5_path = xenium_outs / "cell_feature_matrix.h5"
    cells_path = xenium_outs / "cells.parquet"

    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    if cells_path.exists():
        cells_df = pd.read_parquet(str(cells_path)).set_index("cell_id")
        cells_df.index = cells_df.index.astype(str)
        common = adata.obs_names.intersection(cells_df.index)
        adata = adata[common].copy()
        for col in ["x_centroid", "y_centroid", "nucleus_area", "cell_area"]:
            if col in cells_df.columns:
                adata.obs[col] = cells_df.loc[common, col].values

    logger.info(f"載入 Xenium: {adata.n_obs:,} cells")
    return adata


def subset_xenium_roi(adata, roi: dict):
    """
    根據 ROI 裁切 Xenium AnnData（µm 座標）。
    支援格式 A（fullres px → µm 轉換）和格式 B（直接 µm）。
    """
    if "x_um" in roi:
        x0_um, y0_um, w_um, h_um = roi_fullres_to_um(roi)
    else:
        pxum = roi.get("pixel_size_um", VISIUM_UM_PX)
        x0_um = roi["x"] * pxum
        y0_um = roi["y"] * pxum
        x1_um = (roi["x"] + roi["width_px"]) * pxum
        y1_um = (roi["y"] + roi["height_px"]) * pxum
        w_um = x1_um - x0_um
        h_um = y1_um - y0_um

    x1_um = x0_um + w_um
    y1_um = y0_um + h_um

    mask = (
        (adata.obs["x_centroid"] >= x0_um) &
        (adata.obs["x_centroid"] <  x1_um) &
        (adata.obs["y_centroid"] >= y0_um) &
        (adata.obs["y_centroid"] <  y1_um)
    )
    sub = adata[mask].copy()
    logger.info(f"Xenium ROI 裁切：{mask.sum()} cells")
    return sub


def load_xenium_transcripts_roi(
    xenium_outs: str | Path,
    roi: dict,
) -> "pd.DataFrame":
    """
    使用 pyarrow predicate pushdown 讀取 Xenium transcripts.parquet 並裁切至 ROI。
    """
    import pyarrow.parquet as pq

    xenium_outs = Path(xenium_outs)
    tx_path = xenium_outs / "transcripts.parquet"

    x0_um, y0_um, w_um, h_um = roi_fullres_to_um(roi)
    x1_um, y1_um = x0_um + w_um, y0_um + h_um

    filters = [
        ("x_location", ">=", x0_um), ("x_location", "<", x1_um),
        ("y_location", ">=", y0_um), ("y_location", "<", y1_um),
    ]
    cols = ["feature_name", "x_location", "y_location", "z_location", "qv"]
    df = pq.read_table(str(tx_path), filters=filters, columns=cols).to_pandas()

    logger.info(f"Xenium transcripts ROI: {len(df):,} 轉錄點")
    return df


def rasterize_nucleus_mask(
    xenium_outs: str | Path,
    roi: dict,
) -> np.ndarray:
    """
    從 nucleus_boundaries.parquet 光柵化核遮罩（int32）。
    使用 pyarrow predicate pushdown 加速讀取。

    Returns: (H, W) int32 label 矩陣
    """
    import pyarrow.parquet as pq

    xenium_outs = Path(xenium_outs)
    nuc_path = xenium_outs / "nucleus_boundaries.parquet"

    x0_um, y0_um, w_um, h_um = roi_fullres_to_um(roi)
    x1_um, y1_um = x0_um + w_um, y0_um + h_um

    # 步驟 1：用 bbox 篩選在 ROI 內的 cell_id
    filters_bbox = [
        ("vertex_x", ">=", x0_um), ("vertex_x", "<", x1_um),
        ("vertex_y", ">=", y0_um), ("vertex_y", "<", y1_um),
    ]
    df_in = pq.read_table(str(nuc_path), filters=filters_bbox, columns=["cell_id"]).to_pandas()
    cells_in_roi = df_in["cell_id"].unique().tolist()

    mask_h = round(h_um / XENIUM_UM_PX)
    mask_w = round(w_um / XENIUM_UM_PX)
    mask = np.zeros((mask_h, mask_w), dtype=np.int32)

    if not cells_in_roi:
        return mask

    # 步驟 2：取完整多邊形並光柵化
    df_roi = pq.read_table(
        str(nuc_path),
        filters=[("cell_id", "in", cells_in_roi)],
    ).to_pandas()

    id_to_lbl = {cid: i + 1 for i, cid in enumerate(df_roi["cell_id"].unique())}

    for cell_id, grp in df_roi.groupby("cell_id"):
        label = id_to_lbl[cell_id]
        pts_x = ((grp["vertex_x"].values - x0_um) / XENIUM_UM_PX).astype(np.int32)
        pts_y = ((grp["vertex_y"].values - y0_um) / XENIUM_UM_PX).astype(np.int32)
        pts = np.stack([pts_x, pts_y], axis=1).reshape(-1, 1, 2)
        cv2.fillPoly(mask, [pts], color=label)

    logger.info(f"核遮罩光柵化完成：{len(id_to_lbl)} 個核，shape={mask.shape}")
    return mask


# ── 預覽產圖 ────────────────────────────────────────────────

def get_overview(config: dict) -> dict:
    """
    取得 H&E hires 縮圖的 base64 與座標轉換參數 (scalef)，
    供前端互動式選取 ROI。
    """
    import json as _json
    from pathlib import Path as _Path
    import cv2
    import base64

    binned_002 = _Path(config["paths"].get("binned_002", ""))
    scalef_path = binned_002 / "spatial" / "scalefactors_json.json"
    scalef = 0.1
    mpp = VISIUM_UM_PX
    if scalef_path.exists():
        scalef_data = _json.loads(scalef_path.read_text())
        scalef = scalef_data.get("tissue_hires_scalef", 0.1)
        mpp = scalef_data.get("microns_per_pixel", VISIUM_UM_PX)

    hires_path = binned_002 / "spatial" / "tissue_hires_image.png"
    if not hires_path.exists():
        raise FileNotFoundError(f"找不到 hires 縮圖：{hires_path}")

    img = cv2.imread(str(hires_path))
    h_px, w_px = img.shape[:2]
    
    # 編碼為 base64 JPEG
    _, buf = cv2.imencode(".jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 85])
    img_b64 = base64.b64encode(buf.tobytes()).decode()

    return {
        "image_b64": img_b64,
        "width_hires": w_px,
        "height_hires": h_px,
        "scalef": scalef,
        "microns_per_pixel": mpp
    }


# ── 主裁切執行器 ────────────────────────────────────────────

class RoiExtractor:
    """統一執行所有 ROI 裁切任務"""

    def __init__(self, config: dict):
        self.config = config
        self.rois = config.get("rois", [])

    def run_all(self) -> None:
        """對所有定義的 ROI 執行裁切"""
        logger.info(f"開始裁切 {len(self.rois)} 個 ROI")

        for roi in self.rois:
            name = roi.get("name", "unnamed")
            logger.info(f"處理 ROI: {name}")
            try:
                self._extract_visium(roi)
                self._extract_he_crop(roi)
            except Exception as e:
                logger.error(f"ROI '{name}' 裁切失敗：{e}")

    def _extract_visium(self, roi: dict) -> None:
        """裁切 Visium HD AnnData（2µm 和 8µm）"""
        from pathlib import Path
        import os

        paths = self.config["paths"]
        out_dir = Path(paths["output_dir"]) / "roi" / roi["name"]
        out_dir.mkdir(parents=True, exist_ok=True)

        for bin_size, dir_key in [("002", "binned_002"), ("008", "binned_008")]:
            binned_dir = paths.get(dir_key, "")
            if not binned_dir or not Path(binned_dir).exists():
                logger.warning(f"  {dir_key} 路徑不存在，跳過")
                continue
            adata = load_visium_adata(binned_dir, bin_size)
            sub = subset_anndata_roi(adata, roi, binned_dir=binned_dir)
            out_path = out_dir / f"adata_{bin_size}um.h5ad"
            sub.write_h5ad(str(out_path))
            logger.info(f"  已儲存：{out_path} ({sub.n_obs:,} bins)")

    def _extract_he_crop(self, roi: dict) -> None:
        """裁切 H&E 影像 ROI"""
        from pathlib import Path
        import tifffile

        he_path = Path(self.config["paths"]["he_image"])
        if not he_path.exists():
            logger.warning(f"  H&E 影像不存在：{he_path}，跳過")
            return

        out_dir = Path(self.config["paths"]["output_dir"]) / "roi" / roi["name"]
        out_dir.mkdir(parents=True, exist_ok=True)

        x0, y0, w, h = roi_to_fullres_px(roi)
        crop, ax0, ay0 = read_btf_crop(he_path, x0, y0, w, h)

        out_path = out_dir / "he_crop.tif"
        tifffile.imwrite(str(out_path), crop)
        logger.info(f"  已儲存 H&E crop：{out_path} shape={crop.shape}")
