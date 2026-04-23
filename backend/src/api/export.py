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
                resolved = p.resolve()
                allowed = [output_dir_base.resolve(), data_root_dir.resolve()]
                if not any(str(resolved).startswith(str(r)) for r in allowed):
                    raise ValueError("input_h5ad 路徑必須位於 output_dir 或 data_root 底下")
                h5ad_path = resolved
            else:
                # 嘗試從 output_dir 與 data_root 兩處查找
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError("找不到指定的 h5ad，請確認路徑正確")
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
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                mask_path     = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

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
                roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                poly_path = roi_out_dir / "cellpose_polygons.json"
                with open(poly_path, "w", encoding="utf-8") as f:
                    json.dump(roi_geo, f)
                logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

                # 2. 生成 ROI 局部 µm 的轉錄點 CSV（無全域偏移）
                tx_path = None
                adata_002um_path = roi_out_dir / "adata_002um.h5ad"
                if adata_002um_path.exists():
                    tx_path = _generate_visiumhd_transcripts(
                        adata_002um_path,
                        roi,
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
                resolved = p.resolve()
                allowed = [output_dir_base.resolve(), data_root_dir.resolve()]
                if not any(str(resolved).startswith(str(r)) for r in allowed):
                    raise ValueError("input_h5ad 路徑必須位於 output_dir 或 data_root 底下")
                h5ad_path = resolved
            else:
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError("找不到指定的 h5ad，請確認路徑正確")
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
