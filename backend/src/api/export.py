"""Stage 4: Browser 格式匯出 API（Pipeline 3 版本，使用 Cellpose mask 轉多邊形）"""
import asyncio
import logging
from pathlib import Path
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.export.geometry import mask_to_geojson, shift_geojson_coords
from backend.src.export.inputs import resolve_export_inputs
from backend.src.export.transcripts import generate_visiumhd_transcripts
from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.constants import VISIUM_UM_PX
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.export")

# 自動尋找分析結果 h5ad 的優先序（兩種匯出格式各自不同，刻意不統一：
# 統一會改變 Xenium 實際拿到的 h5ad，屬產品層決策）
XENIUM_CANDIDATES = ("umap_computed.h5ad", "qc_preprocessed.h5ad", "cellpose_cells.h5ad")
LOUPE_CANDIDATES  = ("clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad")

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
        export_dir = resolve_path(paths.get("export_dir", "results/export"))

        inputs          = resolve_export_inputs(config, req.input_h5ad, XENIUM_CANDIDATES)
        h5ad_path       = inputs.h5ad_path
        output_dir_base = inputs.output_dir
        rois            = inputs.rois
        is_merged_mode  = inputs.is_merged
        active_roi      = inputs.active_roi

        if is_merged_mode:
            # ── 多 ROI 模式：每個 ROI 獨立輸出一個 Xenium bundle ──────────────
            # 使用各 ROI 的 he_crop.tif（已裁切，座標從 (0,0) 開始），
            # 避免全域座標偏移的複雜性，確保影像與多邊形對齊。
            logger.info(f"合併模式（{len(rois)} 個 ROI）：每個 ROI 獨立匯出 Xenium bundle")
            import scanpy as sc

            adata_full = sc.read_h5ad(str(h5ad_path))
            logger.info(f"載入完整 h5ad：{len(adata_full)} 個細胞")

            exported_dirs: list[str] = []
            n_rois = len(rois)

            for roi_idx, roi in enumerate(rois):
                rn            = roi.name
                roi_out_dir   = roi.out_dir
                mask_path     = roi.mask_path
                pixel_size_um = roi.pixel_size_um

                _xenium_status = {
                    "status": "running",
                    "progress": roi_idx / n_rois,
                    "message": f"ROI {rn}（{roi_idx + 1}/{n_rois}）匯出中…",
                }

                if not mask_path.exists():
                    logger.warning(f"  [{rn}] 找不到 segmentation_masks.npy，跳過")
                    continue

                # 1. 生成 ROI 局部 µm 的 GeoJSON（無全域偏移）
                logger.info(f"  [{rn}] 生成多邊形…")
                roi_geo = mask_to_geojson(mask_path, pixel_size_um)
                poly_path = roi_out_dir / "cellpose_polygons.json"
                with open(poly_path, "w", encoding="utf-8") as f:
                    json.dump(roi_geo, f)
                logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

                # 2. 生成 ROI 局部 µm 的轉錄點 CSV（無全域偏移）
                tx_path = None
                adata_002um_path = roi_out_dir / "adata_002um.h5ad"
                if adata_002um_path.exists():
                    tx_path = generate_visiumhd_transcripts(
                        adata_002um_path,
                        roi.cfg,
                        roi_out_dir / "transcripts_roi.csv",
                        pixel_size_um,
                    )

                # 3. 從完整 h5ad 取出此 ROI 的子集，重命名 obs_names 為 "cell_N" 格式
                #    以匹配 GeoJSON 的 full_id（mask 輸出為純數字 "N"）
                roi_col = adata_full.obs.get("roi", adata_full.obs.get("roi_name", None))
                if roi_col is not None:
                    roi_mask_bool = roi_col.astype(str) == str(rn)
                else:
                    # fallback：透過 obs_names 前綴篩選
                    roi_mask_bool = adata_full.obs_names.str.startswith(f"{rn}__")
                adata_roi = adata_full[roi_mask_bool].copy()

                if len(adata_roi) == 0:
                    logger.warning(f"  [{rn}] h5ad 中無此 ROI 的細胞，跳過")
                    continue

                # 將 "1__cell_7" → "cell_7"，讓 exporter 的 "cell_N" fallback 對應上
                renamed = []
                for nm in adata_roi.obs_names:
                    if "__cell_" in nm:
                        renamed.append(f"cell_{nm.split('__cell_')[1]}")
                    else:
                        logger.warning(
                            f"  [{rn}] obs_name '{nm}' 不含 '__cell_'，"
                            f"保留原名（可能與 GeoJSON full_id 不符）"
                        )
                        renamed.append(nm)
                adata_roi.obs_names = renamed
                roi_h5ad_path = roi_out_dir / "export_subset.h5ad"
                adata_roi.write_h5ad(str(roi_h5ad_path))
                logger.info(f"  [{rn}] 子集 h5ad：{len(adata_roi)} 個細胞")

                # 4. H&E 底圖：使用已裁切好的 he_crop.tif（座標從 (0,0) 開始）
                he_path = roi_out_dir / "he_crop.tif"

                # 5. 執行匯出，完成後清理臨時 h5ad
                roi_xenium_dir = export_dir / "xenium" / f"roi_{rn}"
                exporter = XeniumExporter(
                    zarr_path=None,
                    poly_json_path=poly_path if poly_path.exists() else None,
                    transcripts_csv_path=tx_path if (tx_path and tx_path.exists()) else None,
                    pixel_size_um=pixel_size_um,
                    he_image_path=he_path if he_path.exists() else None,
                    he_crop_bounds=None,  # he_crop.tif 已裁切，無需偏移
                )
                try:
                    await asyncio.get_running_loop().run_in_executor(
                        None, exporter.export, roi_h5ad_path, roi_xenium_dir,
                    )
                    exported_dirs.append(str(roi_xenium_dir))
                    logger.info(f"  [{rn}] Xenium bundle 完成：{roi_xenium_dir}")
                except Exception as roi_exc:
                    logger.error(f"  [{rn}] Xenium 匯出失敗，跳過此 ROI：{roi_exc}", exc_info=True)
                finally:
                    # 臨時子集 h5ad 無論成功或失敗均清除
                    if roi_h5ad_path.exists():
                        try:
                            roi_h5ad_path.unlink()
                        except OSError:
                            pass

            _xenium_status = {
                "status": "done",
                "progress": 1.0,
                "message": f"Xenium 匯出完成（{len(exported_dirs)} 個 ROI bundle）",
            }
            return

        else:
            roi_name    = active_roi or (rois[-1].name if rois else "")
            roi_in      = inputs.find_roi(roi_name)
            roi_out_dir = roi_in.out_dir if roi_in else output_dir_base
            mask_path   = roi_in.mask_path if roi_in else roi_out_dir / "segmentation_masks.npy"
            roi_cfg       = roi_in.cfg if roi_in else {}
            pixel_size_um = roi_in.pixel_size_um if roi_in else VISIUM_UM_PX

            if not mask_path.exists():
                raise FileNotFoundError(f"找不到 {roi_name} 的 segmentation_masks.npy，請先完成 Stage 1")
            logger.info(f"單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成多邊形...")
            roi_geo = mask_to_geojson(mask_path, pixel_size_um)
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
                transcripts_csv_path = generate_visiumhd_transcripts(
                    adata_002um_path,
                    roi_cfg,
                    roi_out_dir / "transcripts_roi.csv",
                    pixel_size_um,
                )
            else:
                logger.info("未找到 adata_002um.h5ad，不匯出 transcripts 層。")

        # ── 單 ROI 模式：直接匯出 ──────────────────────────────────────────────
        if req.output_dir:
            out_dir = Path(req.output_dir)
        else:
            out_dir = roi_out_dir / "export_xenium"

        exporter = XeniumExporter(
            zarr_path=None,
            poly_json_path=combined_poly_path if (combined_poly_path and combined_poly_path.exists()) else None,
            transcripts_csv_path=transcripts_csv_path if (transcripts_csv_path and transcripts_csv_path.exists()) else None,
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

        paths      = config.get("paths", {})
        export_dir = resolve_path(paths.get("export_dir", "results/export"))
        whitelist  = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

        inputs          = resolve_export_inputs(config, req.input_h5ad, LOUPE_CANDIDATES)
        h5ad_path       = inputs.h5ad_path
        output_dir_base = inputs.output_dir
        rois            = inputs.rois
        is_merged_mode  = inputs.is_merged
        active_roi      = inputs.active_roi

        poly_json_path: "Path | None" = None
        import json

        if is_merged_mode:
            logger.info("Loupe 匯出：合併模式，產生 combined_cellpose_polygons.json")
            all_features: list = []
            for roi in rois:
                rn            = roi.name
                mask_path     = roi.mask_path
                pixel_size_um = roi.pixel_size_um

                if mask_path.exists():
                    roi_geo = mask_to_geojson(mask_path, pixel_size_um)
                else:
                    logger.warning(f"  [{rn}] 找不到 segmentation_masks.npy，跳過")
                    continue

                roi_x_um = roi.cfg.get("x", 0) * pixel_size_um
                roi_y_um = roi.cfg.get("y", 0) * pixel_size_um
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)

            poly_json_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(poly_json_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)
        else:
            roi_name    = active_roi or (rois[-1].name if rois else "")
            roi_in      = inputs.find_roi(roi_name)
            roi_out_dir = roi_in.out_dir if roi_in else output_dir_base
            mask_path   = roi_in.mask_path if roi_in else roi_out_dir / "segmentation_masks.npy"
            pixel_size_um = roi_in.pixel_size_um if roi_in else VISIUM_UM_PX

            if mask_path.exists():
                logger.info(f"Loupe 匯出：單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成")
                roi_geo = mask_to_geojson(mask_path, pixel_size_um)
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
