"""Stage 4: Browser 格式匯出 API（Pipeline 3 版本，使用 Cellpose mask 轉多邊形）"""
import asyncio
import logging
from pathlib import Path
from typing import Optional
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.constants import VISIUM_UM_PX
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.export")

_xenium_status = {"status": "idle", "progress": 0.0, "message": ""}
_loupe_status  = {"status": "idle", "progress": 0.0, "message": ""}
_result_status: dict = {"status": "idle", "progress": 0.0, "message": ""}
_result_images: dict = {}
_xenium_lock   = asyncio.Lock()
_loupe_lock    = asyncio.Lock()
_result_lock   = asyncio.Lock()


class ExportRequest(BaseModel):
    input_h5ad: str = ""   # 空字串 = 使用 config 預設輸出路徑
    output_dir: str = ""


@router.get("/status/xenium")
async def xenium_status():
    return _xenium_status


@router.get("/status/loupe")
async def loupe_status():
    return _loupe_status


# ──────────────────────────────────────────────────────────────────────────────
# 核心：Cellpose mask → GeoJSON 多邊形
# ──────────────────────────────────────────────────────────────────────────────

def _mask_to_geojson(
    mask_path: Path,
    pixel_size_um: float,
    min_area_px: int = 20,
) -> dict:
    """
    將 segmentation_masks.npy 轉換為 GeoJSON FeatureCollection。

    座標：ROI 局部 µm（原點 = ROI 左上角），與 cellpose_cells.h5ad obsm['spatial'] 一致。
    使用 regionprops 取 bounding box 後在小 patch 上做輪廓偵測，
    避免 O(n_cells × H×W) 的全圖掃描。
    """
    import numpy as np
    from skimage import measure

    seg_mask = np.load(str(mask_path))

    features = []
    # regionprops 一次性計算 bounding box + area，避免逐細胞全圖掃描
    for prop in measure.regionprops(seg_mask):
        if prop.area < min_area_px:
            continue
        cid = prop.label
        r0, c0, r1, c1 = prop.bbox

        # 在 bounding box patch 上找輪廓（比全圖快數個量級）
        cell_crop = (seg_mask[r0:r1, c0:c1] == cid).astype(np.uint8)
        padded = np.pad(cell_crop, 1, mode="constant")
        contours = measure.find_contours(padded, 0.5)
        if not contours:
            continue

        contour = max(contours, key=len)
        # 還原 padding(1) + bounding box offset，再轉成 (x, y) µm
        xy_um = np.column_stack([
            (contour[:, 1] - 1 + c0) * pixel_size_um,   # col → x
            (contour[:, 0] - 1 + r0) * pixel_size_um,   # row → y
        ])

        if not np.allclose(xy_um[0], xy_um[-1]):
            xy_um = np.vstack([xy_um, xy_um[0]])

        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [xy_um.tolist()],
            },
            "properties": {
                "full_id":   str(int(cid)),
                "cell_id":   int(cid),
            },
        })

    logger.info(f"  生成 {len(features)} 個 Cellpose 多邊形")
    return {"type": "FeatureCollection", "features": features}


# ──────────────────────────────────────────────────────────────────────────────
# Xenium 匯出
# ──────────────────────────────────────────────────────────────────────────────


async def _run_xenium(config: dict, req: ExportRequest):
    global _xenium_status
    set_current_stage("export")
    _xenium_status = {"status": "running", "progress": 0.0, "message": "匯出至 Xenium Explorer..."}
    try:
        from backend.src.export.xenium_exporter import XeniumExporter
        import json

        paths = config.get("paths", {})
        rois  = config.get("rois", [{}])

        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir      = resolve_path(paths.get("export_dir", "results/export"))
        data_root_dir   = resolve_path(paths.get("data_root", "."))

        # Find h5ad_path first
        h5ad_path = None
        if req.input_h5ad:
            p = Path(req.input_h5ad)
            if p.is_absolute():
                h5ad_path = p
            else:
                # 嘗試從 output_dir 與 data_root 兩處查找
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError(f"找不到指定的 h5ad：{req.input_h5ad}\n搜尋位置：{candidate_a} 或 {candidate_b}")
        else:
            for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad", "cellpose_cells.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
            # 新路徑：roi/{roi_name}/ 底下
            if h5ad_path is None:
                for roi in rois:
                    roi_name = roi.get("name", "")
                    if not roi_name:
                        continue
                    roi_dir = output_dir_base / "roi" / roi_name
                    for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                        p = roi_dir / candidate
                        if p.exists():
                            h5ad_path = p
                            break
                    if h5ad_path is not None:
                        break
        if h5ad_path is None:
            raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 3 分析。\n搜尋位置：{output_dir_base}")

        # Check if the h5ad is merged mode by inspecting obs_names
        import scanpy as sc
        adata_head = sc.read_h5ad(str(h5ad_path), backed="r")
        first_obs = adata_head.obs_names[0] if len(adata_head) > 0 else ""
        is_merged_mode = "__" in first_obs
        active_roi = adata_head.uns.get("active_roi", None) if "active_roi" in adata_head.uns else None
        del adata_head # Release backed file handle

        combined_poly_path: "Path | None" = None

        if is_merged_mode:
            logger.info(f"合併模式（{len(rois)} 個 ROI），根據 h5ad 的 obs_names 判定")
            all_features: list = []
            for roi in rois:
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                mask_path   = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

                if mask_path.exists():
                    logger.info(f"  [{rn}] 使用 MCseg v2 遮罩生成多邊形")
                    roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                else:
                    logger.warning(f"  [{rn}] 找不到 segmentation_masks.npy，跳過")
                    continue

                # 加入全域座標偏移
                roi_x_um = roi.get("x", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                roi_y_um = roi.get("y", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    _shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)
                logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

            combined_poly_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)

            roi_pixel_size_um = rois[0].get("pixel_size_um", VISIUM_UM_PX) if rois else VISIUM_UM_PX

            he_image_path = resolve_path(paths.get("he_image", "")) if paths.get("he_image") else None
            he_crop_bounds = None
            if he_image_path and rois:
                _x0 = min(r.get("x", 0) for r in rois)
                _y0 = min(r.get("y", 0) for r in rois)
                _x1 = max(r.get("x", 0) + r.get("width_px", 0) for r in rois)
                _y1 = max(r.get("y", 0) + r.get("height_px", 0) for r in rois)
                he_crop_bounds = (_x0, _y0, _x1, _y1)

        else:
            roi_name    = active_roi or (rois[-1].get("name", "") if rois else "")
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            mask_path   = roi_out_dir / "segmentation_masks.npy"
            roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
            pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))

            if not mask_path.exists():
                raise FileNotFoundError(f"找不到 {roi_name} 的 segmentation_masks.npy，請先完成 Stage 1")
            logger.info(f"單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成多邊形...")
            roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
            combined_poly_path = roi_out_dir / "cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump(roi_geo, f)

            roi_pixel_size_um = pixel_size_um
            he_image_path     = roi_out_dir / "he_crop.tif"
            if not he_image_path.exists():
                logger.warning(f"單 ROI 模式找不到 he_crop.tif，將不用底圖匯出")
                he_image_path = None
            he_crop_bounds    = None

            # ── 從 Visium HD 2µm bins 生成轉錄點 ──────────────────────────
            transcripts_csv_path = None
            adata_002um_path = roi_out_dir / "adata_002um.h5ad"
            if adata_002um_path.exists():
                transcripts_csv_path = _generate_visiumhd_transcripts(
                    adata_002um_path,
                    roi_cfg,
                    roi_out_dir / "transcripts_roi.csv",
                    pixel_size_um,
                )
            else:
                logger.info("未找到 adata_002um.h5ad，不匯出 transcripts 層。")

        if req.output_dir:
            out_dir = Path(req.output_dir)
        else:
            out_dir = export_dir / "xenium" if is_merged_mode else roi_out_dir / "export_xenium"

        exporter = XeniumExporter(
            zarr_path=None,
            poly_json_path=combined_poly_path if (combined_poly_path and combined_poly_path.exists()) else None,
            transcripts_csv_path=transcripts_csv_path if not is_merged_mode else None,
            pixel_size_um=roi_pixel_size_um,
            he_image_path=he_image_path,
            he_crop_bounds=he_crop_bounds,
        )
        await asyncio.get_running_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _xenium_status = {"status": "done", "progress": 1.0, "message": "Xenium 匯出完成"}
    except Exception as e:
        logger.error(f"Xenium 匯出失敗：{e}", exc_info=True)
        _xenium_status = {"status": "error", "progress": 0.0, "message": "Xenium 匯出失敗，請查閱 log"}


# ──────────────────────────────────────────────────────────────────────────────
# Loupe 匯出
# ──────────────────────────────────────────────────────────────────────────────

async def _run_loupe(config: dict, req: ExportRequest):
    global _loupe_status
    set_current_stage("export")
    _loupe_status = {"status": "running", "progress": 0.0, "message": "匯出至 Loupe Browser..."}
    try:
        from backend.src.export.loupe_exporter import LoupeExporter

        paths    = config.get("paths", {})
        rois     = config.get("rois", [{}])
        roi_name = rois[0].get("name", "") if rois else ""

        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir      = resolve_path(paths.get("export_dir", "results/export"))
        data_root_dir   = resolve_path(paths.get("data_root", "."))
        whitelist       = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

        # Find h5ad_path first
        h5ad_path = None
        if req.input_h5ad:
            p = Path(req.input_h5ad)
            if p.is_absolute():
                h5ad_path = p
            else:
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError(f"找不到指定的 h5ad：{req.input_h5ad}")
        else:
            for candidate in ["clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
            # 新路徑：roi/{roi_name}/ 底下
            if h5ad_path is None:
                for roi in rois:
                    roi_name = roi.get("name", "")
                    if not roi_name:
                        continue
                    roi_dir = output_dir_base / "roi" / roi_name
                    for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                        p = roi_dir / candidate
                        if p.exists():
                            h5ad_path = p
                            break
                    if h5ad_path is not None:
                        break
            if h5ad_path is None:
                raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 3 分析。搜尋位置：{output_dir_base}")

        # Check if the h5ad is merged mode by inspecting obs_names
        import scanpy as sc
        adata_head = sc.read_h5ad(str(h5ad_path), backed="r")
        first_obs = adata_head.obs_names[0] if len(adata_head) > 0 else ""
        is_merged_mode = "__" in first_obs
        active_roi = adata_head.uns.get("active_roi", None) if "active_roi" in adata_head.uns else None
        del adata_head

        poly_json_path: "Path | None" = None
        import json

        if is_merged_mode:
            logger.info("Loupe 匯出：合併模式，產生 combined_cellpose_polygons.json")
            all_features: list = []
            for roi in rois:
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                mask_path     = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

                if mask_path.exists():
                    roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                else:
                    logger.warning(f"  [{rn}] 找不到 segmentation_masks.npy，跳過")
                    continue

                roi_x_um = roi.get("x", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                roi_y_um = roi.get("y", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    _shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)

            poly_json_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(poly_json_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)
        else:
            roi_name    = active_roi or (rois[-1].get("name", "") if rois else "")
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            mask_path   = roi_out_dir / "segmentation_masks.npy"
            roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
            pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))

            if mask_path.exists():
                logger.info(f"Loupe 匯出：單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成")
                roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                poly_json_path = roi_out_dir / "cellpose_polygons.json"
                with open(poly_json_path, "w", encoding="utf-8") as f:
                    json.dump(roi_geo, f)
            else:
                logger.warning(f"找不到 {roi_name} 的 segmentation_masks.npy，將不匯出多邊形層")

        if req.output_dir:
            out_dir = Path(req.output_dir)
        else:
            out_dir = export_dir / "loupe" if is_merged_mode else roi_out_dir / "export_loupe"

        exporter = LoupeExporter(
            poly_json_path=poly_json_path if poly_json_path and poly_json_path.exists() else None,
            whitelist_path=resolve_path(whitelist) if whitelist else None,
        )
        await asyncio.get_running_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _loupe_status = {"status": "done", "progress": 1.0, "message": "Loupe 匯出完成"}
    except Exception as e:
        logger.error(f"Loupe 匯出失敗：{e}", exc_info=True)
        _loupe_status = {"status": "error", "progress": 0.0, "message": "Loupe 匯出失敗，請查閱 log"}


@router.post("/xenium")
async def export_xenium(req: ExportRequest, background_tasks: BackgroundTasks):
    global _xenium_status
    async with _xenium_lock:
        if _xenium_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        _xenium_status = {"status": "running", "progress": 0.0, "message": "準備匯出..."}
        background_tasks.add_task(_run_xenium, config, req)
    return {"status": "ok", "message": "Xenium 匯出已啟動"}


@router.post("/loupe")
async def export_loupe(req: ExportRequest, background_tasks: BackgroundTasks):
    global _loupe_status
    async with _loupe_lock:
        if _loupe_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        _loupe_status = {"status": "running", "progress": 0.0, "message": "準備匯出..."}
        background_tasks.add_task(_run_loupe, config, req)
    return {"status": "ok", "message": "Loupe 匯出已啟動"}


# ──────────────────────────────────────────────────────────────────────────────
# 工具函式
# ──────────────────────────────────────────────────────────────────────────────

def _generate_visiumhd_transcripts(
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

        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(str(out_path), index=False)
        logger.info(f"已寫出 {len(df):,} 個 Visium HD 轉錄點至 {out_path}")
        return out_path

    except Exception as exc:
        logger.warning(f"Visium HD transcripts 生成失敗（繼續執行）：{exc}")
        return None


def _shift_geojson_coords(feat: dict, dx: float, dy: float):
    """In-place 平移 GeoJSON feature 的座標。"""
    def _shift(coords):
        if not coords:
            return coords
        if isinstance(coords[0], (int, float)):
            return [coords[0] + dx, coords[1] + dy] + list(coords[2:])
        return [_shift(c) for c in coords]

    geom = feat.get("geometry", {})
    if geom and geom.get("coordinates") is not None:
        geom["coordinates"] = _shift(geom["coordinates"])


# ──────────────────────────────────────────────────────────────────────────────
# Result Visualizations（標註後視覺化）
# ──────────────────────────────────────────────────────────────────────────────

@router.get("/result_status")
async def get_result_status():
    return _result_status


@router.get("/result_images")
async def get_result_images():
    global _result_status, _result_images
    if not _result_images:
        try:
            fig_dir = resolve_path(load_config()["paths"]["figure_dir"])
            loaded: dict[str, str] = {}
            import base64
            # 固定圖
            for key, fname in [
                ("result_umap",     "result_umap.png"),
                ("result_dotplot",  "result_dotplot.png"),
                ("result_heatmap",  "result_heatmap.png"),
            ]:
                p = fig_dir / fname
                if p.exists():
                    loaded[key] = base64.b64encode(p.read_bytes()).decode()
            # 空間圖（支援單/多 ROI）
            import re as _re
            for p in sorted(fig_dir.glob("result_spatial*.png")):
                m = _re.match(r"(result_spatial(?:_filled)?(?:_.+)?)\.png$", p.name)
                if m:
                    loaded[m.group(1)] = base64.b64encode(p.read_bytes()).decode()
            if loaded:
                _result_images = loaded
                _result_status = {"status": "done", "progress": 1.0, "message": "已從磁碟載入結果圖"}
        except Exception as e:
            logger.warning(f"從磁碟載入結果圖失敗：{e}")
    if not _result_images:
        return {"status": "error", "message": "結果圖尚未產生，請先執行「生成結果圖」"}
    return {"status": "ok", "data": _result_images}


@router.post("/generate_result")
async def generate_result(background_tasks: BackgroundTasks):
    global _result_status
    async with _result_lock:
        if _result_status.get("status") == "running":
            return {"status": "running", "message": "已在執行中"}
        config = load_config()
        _result_status = {"status": "running", "progress": 0.0, "message": "生成結果視覺化中..."}
        background_tasks.add_task(_run_generate_result, config)
    return {"status": "started"}


async def _run_generate_result(config: dict):
    global _result_status, _result_images
    set_current_stage("export")
    try:
        from backend.src.analysis.pipeline import run_result_visualizations
        result = await asyncio.get_running_loop().run_in_executor(
            None, run_result_visualizations, config
        )
        _result_images = result
        _result_status = {
            "status": "done",
            "progress": 1.0,
            "message": f"結果圖生成完成（{len(result)} 張）",
        }
    except Exception as e:
        import traceback
        logger.error(f"結果視覺化失敗：{e!r}\n{traceback.format_exc()}")
        _result_status = {"status": "error", "progress": 0.0, "message": "結果視覺化失敗，請查閱 log"}
