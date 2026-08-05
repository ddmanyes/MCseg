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
import logging
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
        # strip-based TIFF 的 tilewidth 屬性存在但為 0，getattr 的預設值不會生效 ——
        # 直接往下算會得到 ZeroDivisionError，錯誤訊息看不出真正原因
        TW = getattr(page, "tilewidth",  512) or 0
        TH = getattr(page, "tilelength", 512) or 0
        if not TW or not TH:
            raise NotImplementedError(
                f"TIFF 檔 '{btf_path}' 不是 tiled 格式（tilewidth={TW}, tilelength={TH}），"
                "無法進行 tile-based 讀取。請改用 read_image_crop() 自動分派 strip 讀取。"
            )
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


def read_image_crop(
    img_path: Path,
    x0: int,
    y0: int,
    w: int,
    h: int,
    margin: int = 0,
) -> tuple[np.ndarray, int, int]:
    """
    讀取 TIFF ROI crop，自動分派 tiled / strip 兩種版面（皆不全圖載入）。

    Space Ranger 輸出的 H&E 有 tiled（BTF）與 strip-based（一般 TIFF）兩種，
    對後者呼叫 `read_btf_crop` 會失敗。所有裁切點一律走這個函數。

    Returns
    -------
    (crop_rgb, actual_x0, actual_y0)
        tiled 版面會對齊 tile 邊界，實際原點可能與請求值不同 —— 下游座標
        必須以回傳的 origin 為準，用請求值會整批偏移。
    """
    import tifffile

    with tifffile.TiffFile(str(img_path)) as tf:
        is_tiled = bool(tf.pages[0].tags.get("TileOffsets"))

    if is_tiled:
        return read_btf_crop(img_path, x0, y0, w, h, margin=margin)

    from backend.src.roi.tile_server import read_strip_crop

    fx0 = max(0, x0 - margin)
    fy0 = max(0, y0 - margin)
    crop = read_strip_crop(img_path, fx0, fy0, w + (x0 - fx0) + margin, h + (y0 - fy0) + margin)
    return crop, fx0, fy0


# ── Visium HD AnnData 裁切 ──────────────────────────────────

def _resolve_roi_bin_transform(config: dict) -> tuple[np.ndarray, str]:
    """
    取得「Space Ranger fullres px → he_image（ROI 畫框所在的 raw TIFF）px」的變換。

    與 Stage 2 全片計數（`fullslide.pipeline.resolve_bin_to_image_transform`）走
    同一份決策邏輯（對位 JSON 優先、近似縮放 fallback，見該函式 docstring）——
    Stage 0 互動式 ROI 裁切過去自行維護一份只認單一寫死檔名
    （`H1-WGR3TC4-D1-fiducials-image-registration.json`，dpcp01 樣本殘留）的轉換，
    換一個 mpp≠0.2737 的樣本就靜默退化成「兩個座標系當同一個」，每次都篩到
    0 個 bin（2026-08-05，康育 SDS-D0D1D2 實測）。
    """
    import tifffile
    from backend.src.fullslide.pipeline import resolve_dense_bin_to_image_transform

    he_path = Path(config.get("paths", {}).get("he_image", ""))
    if not he_path.exists():
        raise FileNotFoundError(f"找不到 H&E 影像：{he_path}")

    with tifffile.TiffFile(str(he_path)) as tf:
        page = tf.pages[0]
        full_shape = (int(page.imagelength), int(page.imagewidth))

    return resolve_dense_bin_to_image_transform(config, full_shape)


def subset_anndata_roi(adata, roi: dict, config: dict | None = None):
    """
    根據 ROI fullres pixel 座標，裁切 Visium HD AnnData。

    ROI 的 `x`/`y`/`width_px`/`height_px` 是 UI 在 `he_image`（raw TIFF）上畫框
    得到的座標；`adata.obs` 的 `pxl_col_in_fullres`/`pxl_row_in_fullres` 是
    **Space Ranger 自己的 fullres 座標系**，兩者未必同尺度（CytAssist 樣本常見，
    可能差 1.5～2 倍甚至更多，取決於掃描機解析度）。篩選前一律先用
    `_resolve_roi_bin_transform` 把 bin 投影到 he_image 座標系再比對。

    投影後的座標會寫回 `obsm['spatial']`——下游（Stage 2 `counter.py`、
    `export/transcripts.py` 等）都假設這個值跟 ROI 的 `x`/`y` 同一座標系，這裡
    是唯一該做轉換的地方，下游不需要、也不應該再各自轉一次（CLAUDE.md §11 DRY）。

    Parameters
    ----------
    config : 完整 pipeline config（需含 `paths.he_image`／`paths.binned_002`／
        `alignment.*`）。為 `None` 時無法解析變換，退回「兩座標系同尺度」的
        假設並記警告——只有 mpp 剛好等於專案預設 `VISIUM_UM_PX`（如官方 CRC
        樣本）時這個假設才成立，其餘情況會靜默篩到 0 bin。
    """
    x0, y0, w, h = roi_to_fullres_px(roi)
    x1, y1 = x0 + w, y0 + h

    sr_col = adata.obs["pxl_col_in_fullres"].values.astype(float)
    sr_row = adata.obs["pxl_row_in_fullres"].values.astype(float)

    img_col, img_row = sr_col, sr_row
    if config is not None:
        try:
            transform, source = _resolve_roi_bin_transform(config)
            pts = np.vstack([sr_col, sr_row, np.ones_like(sr_col)])
            proj = transform @ pts
            img_col = proj[0, :] / proj[2, :]
            img_row = proj[1, :] / proj[2, :]
            logger.info(f"ROI 對位變換：{source}")
        except Exception as e:   # noqa: BLE001 — 診斷輔助不得中斷裁切（CLAUDE.md §11）
            logger.warning(
                f"ROI 對位變換解析失敗（{type(e).__name__}: {e}），"
                "退回兩座標系同尺度假設（僅官方 mpp=0.2737 樣本才正確）"
            )
    else:
        logger.warning("subset_anndata_roi 未收到 config，無法解析對位變換，退回兩座標系同尺度假設")

    mask = (
        (img_col >= x0) & (img_col < x1) &
        (img_row >= y0) & (img_row < y1)
    )
    sub = adata[mask].copy()
    sub.obsm["spatial"] = np.stack([img_col[mask], img_row[mask]], axis=1)
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

    def run_all(self) -> dict:
        """
        對所有定義的 ROI 執行裁切。

        回傳 ``{"total", "ok", "failed"}``；``failed`` 為失敗的 ROI 名稱清單。
        呼叫端必須依此判定成敗 —— 單一 ROI 失敗只記 log 不中斷，
        但整體狀態不得回報「完成」（否則畫面顯示成功、實際沒有產出）。
        """
        logger.info(f"開始裁切 {len(self.rois)} 個 ROI")
        if not self.rois:
            raise ValueError("ROI 清單為空，沒有任何可裁切的區域")

        failed: list[str] = []
        for roi in self.rois:
            name = roi.get("name", "unnamed")
            logger.info(f"處理 ROI: {name}")
            try:
                self._extract_visium(roi)
                self._extract_he_crop(roi)
            except Exception as e:
                logger.error(f"ROI '{name}' 裁切失敗：{e}")
                failed.append(name)

        return {"total": len(self.rois), "ok": len(self.rois) - len(failed), "failed": failed}

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
            sub = subset_anndata_roi(adata, roi, config=self.config)
            out_path = out_dir / f"adata_{bin_size}um.h5ad"
            sub.write_h5ad(str(out_path))
            logger.info(f"  已儲存：{out_path} ({sub.n_obs:,} bins)")
            if sub.n_obs == 0:
                # 框在組織外／座標系搞錯時會產生空 h5ad，要等到 Stage 2 才爆炸
                logger.warning(
                    f"  ⚠️ ROI '{roi['name']}' 的 {bin_size}um 沒有命中任何 bin —— "
                    f"請確認框選範圍落在組織上，且座標為 fullres pixel"
                )

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
        crop, ax0, ay0 = read_image_crop(he_path, x0, y0, w, h)

        out_path = out_dir / "he_crop.tif"
        tifffile.imwrite(str(out_path), crop)
        logger.info(f"  已儲存 H&E crop：{out_path} shape={crop.shape}")
