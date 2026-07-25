"""Stage 4 匯出的編排：masks + adata → Xenium / Loupe bundle。

「一次匯出要做哪些事、以什麼順序」的唯一擁有者。純同步、不依賴 FastAPI 與
event loop——呼叫端（API 背景任務、未來的 CLI）自行決定要不要丟進 thread。
進度靠注入的 progress callback 回報，這裡不知道 HTTP 狀態表示法長什麼樣。
"""
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from backend.src.export.geometry import mask_to_geojson, shift_geojson_coords
from backend.src.export.inputs import resolve_export_inputs
from backend.src.export.transcripts import generate_visiumhd_transcripts
from backend.src.utils.config import resolve_path
from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.export.jobs")

# (完成比例 0.0~1.0, 給人看的訊息)
ProgressFn = Callable[[float, str], None]

# 自動尋找分析結果 h5ad 的優先序（兩種格式各自不同，刻意不統一：
# 統一會改變 Xenium 實際拿到的 h5ad，屬產品層決策）
XENIUM_CANDIDATES = ("umap_computed.h5ad", "qc_preprocessed.h5ad", "cellpose_cells.h5ad")
LOUPE_CANDIDATES  = ("clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad")


@dataclass(frozen=True)
class ExportResult:
    """一次匯出的產出。"""
    output_dirs: list[Path]
    is_merged: bool


def _noop_progress(fraction: float, message: str) -> None:
    pass


def run_xenium_export(
    config: dict,
    input_h5ad: str = "",
    output_dir: str = "",
    progress: "Optional[ProgressFn]" = None,
) -> ExportResult:
    """匯出至 Xenium Explorer。

    多 ROI（合併）模式下每個 ROI 各出一個 bundle，單一 ROI 失敗只跳過該 ROI；
    單 ROI 模式則直接匯出一個 bundle。
    """
    from backend.src.export.xenium_exporter import XeniumExporter

    report = progress or _noop_progress

    export_dir = resolve_path(config.get("paths", {}).get("export_dir", "results/export"))
    inputs = resolve_export_inputs(config, input_h5ad, XENIUM_CANDIDATES)

    if inputs.is_merged:
        return _run_xenium_merged(inputs, export_dir, report, XeniumExporter)
    return _run_xenium_single(inputs, output_dir, report, XeniumExporter)


def _run_xenium_merged(inputs, export_dir: Path, report: ProgressFn, XeniumExporter) -> ExportResult:
    """多 ROI 模式：每個 ROI 獨立輸出一個 Xenium bundle。

    使用各 ROI 的 he_crop.tif（已裁切，座標從 (0,0) 開始），
    避免全域座標偏移的複雜性，確保影像與多邊形對齊。
    """
    import scanpy as sc

    logger.info(f"合併模式（{len(inputs.rois)} 個 ROI）：每個 ROI 獨立匯出 Xenium bundle")
    adata_full = sc.read_h5ad(str(inputs.h5ad_path))
    logger.info(f"載入完整 h5ad：{len(adata_full)} 個細胞")

    exported_dirs: list[Path] = []
    n_rois = len(inputs.rois)

    for roi_idx, roi in enumerate(inputs.rois):
        rn = roi.name
        report(roi_idx / n_rois, f"ROI {rn}（{roi_idx + 1}/{n_rois}）匯出中…")

        if not roi.mask_path.exists():
            logger.warning(f"  [{rn}] 找不到 segmentation_masks.npy，跳過")
            continue

        # 1. 生成 ROI 局部 µm 的 GeoJSON（無全域偏移）
        logger.info(f"  [{rn}] 生成多邊形…")
        roi_geo = mask_to_geojson(roi.mask_path, roi.pixel_size_um)
        poly_path = roi.out_dir / "cellpose_polygons.json"
        with open(poly_path, "w", encoding="utf-8") as f:
            json.dump(roi_geo, f)
        logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

        # 2. 生成 ROI 局部 µm 的轉錄點 CSV（無全域偏移）
        tx_path = None
        adata_002um_path = roi.out_dir / "adata_002um.h5ad"
        if adata_002um_path.exists():
            tx_path = generate_visiumhd_transcripts(
                adata_002um_path,
                roi.cfg,
                roi.out_dir / "transcripts_roi.csv",
                roi.pixel_size_um,
            )

        # 3. 從完整 h5ad 取出此 ROI 的子集，重命名 obs_names 為 "cell_N" 格式
        #    以匹配 GeoJSON 的 full_id（mask 輸出為純數字 "N"）
        adata_roi = _subset_roi(adata_full, rn)
        if len(adata_roi) == 0:
            logger.warning(f"  [{rn}] h5ad 中無此 ROI 的細胞，跳過")
            continue

        roi_h5ad_path = roi.out_dir / "export_subset.h5ad"
        adata_roi.write_h5ad(str(roi_h5ad_path))
        logger.info(f"  [{rn}] 子集 h5ad：{len(adata_roi)} 個細胞")

        # 4. H&E 底圖：使用已裁切好的 he_crop.tif（座標從 (0,0) 開始）
        he_path = roi.out_dir / "he_crop.tif"

        # 5. 執行匯出，完成後清理臨時 h5ad
        roi_xenium_dir = export_dir / "xenium" / f"roi_{rn}"
        exporter = XeniumExporter(
            zarr_path=None,
            poly_json_path=poly_path if poly_path.exists() else None,
            transcripts_csv_path=tx_path if (tx_path and tx_path.exists()) else None,
            pixel_size_um=roi.pixel_size_um,
            he_image_path=he_path if he_path.exists() else None,
            he_crop_bounds=None,  # he_crop.tif 已裁切，無需偏移
        )
        try:
            exporter.export(roi_h5ad_path, roi_xenium_dir)
            exported_dirs.append(roi_xenium_dir)
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

    return ExportResult(output_dirs=exported_dirs, is_merged=True)


def _subset_roi(adata_full, roi_name: str):
    """取出單一 ROI 的細胞，並把 "1__cell_7" 還原成 "cell_7"。

    exporter 以 GeoJSON 的 full_id（mask 輸出為純數字）對應細胞，
    合併後的 obs_names 帶 ROI 前綴，必須先剝掉才對得上。
    """
    roi_col = adata_full.obs.get("roi", adata_full.obs.get("roi_name", None))
    if roi_col is not None:
        roi_mask_bool = roi_col.astype(str) == str(roi_name)
    else:
        # fallback：透過 obs_names 前綴篩選
        roi_mask_bool = adata_full.obs_names.str.startswith(f"{roi_name}__")
    adata_roi = adata_full[roi_mask_bool].copy()
    if len(adata_roi) == 0:
        return adata_roi

    renamed = []
    for nm in adata_roi.obs_names:
        if "__cell_" in nm:
            renamed.append(f"cell_{nm.split('__cell_')[1]}")
        else:
            logger.warning(
                f"  [{roi_name}] obs_name '{nm}' 不含 '__cell_'，"
                f"保留原名（可能與 GeoJSON full_id 不符）"
            )
            renamed.append(nm)
    adata_roi.obs_names = renamed
    return adata_roi


def _run_xenium_single(inputs, output_dir: str, report: ProgressFn, XeniumExporter) -> ExportResult:
    """單 ROI 模式：直接匯出一個 bundle。"""
    roi_name = inputs.active_roi or (inputs.rois[-1].name if inputs.rois else "")
    roi = inputs.find_roi(roi_name)
    roi_out_dir   = roi.out_dir if roi else inputs.output_dir
    mask_path     = roi.mask_path if roi else roi_out_dir / "segmentation_masks.npy"
    roi_cfg       = roi.cfg if roi else {}
    pixel_size_um = roi.pixel_size_um if roi else VISIUM_UM_PX

    if not mask_path.exists():
        raise FileNotFoundError(f"找不到 {roi_name} 的 segmentation_masks.npy，請先完成 Stage 1")

    report(0.1, f"單 ROI 模式（{roi_name}）：生成多邊形…")
    logger.info(f"單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成多邊形...")
    roi_geo = mask_to_geojson(mask_path, pixel_size_um)
    poly_path = roi_out_dir / "cellpose_polygons.json"
    with open(poly_path, "w", encoding="utf-8") as f:
        json.dump(roi_geo, f)

    he_image_path = roi_out_dir / "he_crop.tif"
    if not he_image_path.exists():
        logger.warning("單 ROI 模式找不到 he_crop.tif，將不用底圖匯出")
        he_image_path = None

    # ── 從 Visium HD 2µm bins 生成轉錄點 ──────────────────────────
    transcripts_csv_path = None
    adata_002um_path = roi_out_dir / "adata_002um.h5ad"
    if adata_002um_path.exists():
        report(0.4, "生成轉錄點…")
        transcripts_csv_path = generate_visiumhd_transcripts(
            adata_002um_path,
            roi_cfg,
            roi_out_dir / "transcripts_roi.csv",
            pixel_size_um,
        )
    else:
        logger.info("未找到 adata_002um.h5ad，不匯出 transcripts 層。")

    out_dir = Path(output_dir) if output_dir else roi_out_dir / "export_xenium"

    report(0.6, "寫出 Xenium bundle…")
    exporter = XeniumExporter(
        zarr_path=None,
        poly_json_path=poly_path if poly_path.exists() else None,
        transcripts_csv_path=transcripts_csv_path if (transcripts_csv_path and transcripts_csv_path.exists()) else None,
        pixel_size_um=pixel_size_um,
        he_image_path=he_image_path,
        he_crop_bounds=None,
    )
    exporter.export(inputs.h5ad_path, out_dir)
    return ExportResult(output_dirs=[out_dir], is_merged=False)


def run_loupe_export(
    config: dict,
    input_h5ad: str = "",
    output_dir: str = "",
    progress: "Optional[ProgressFn]" = None,
) -> ExportResult:
    """匯出至 Loupe Browser。

    多 ROI（合併）模式下所有 ROI 的多邊形平移到全域座標後併成一份 GeoJSON；
    單 ROI 模式則直接用該 ROI 的局部座標。
    """
    from backend.src.export.loupe_exporter import LoupeExporter

    report = progress or _noop_progress

    export_dir = resolve_path(config.get("paths", {}).get("export_dir", "results/export"))
    whitelist  = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

    inputs = resolve_export_inputs(config, input_h5ad, LOUPE_CANDIDATES)

    poly_json_path: "Optional[Path]" = None
    roi_out_dir = inputs.output_dir

    if inputs.is_merged:
        report(0.2, "合併各 ROI 多邊形…")
        poly_json_path = _combine_roi_polygons(inputs)
    else:
        roi_name = inputs.active_roi or (inputs.rois[-1].name if inputs.rois else "")
        roi = inputs.find_roi(roi_name)
        roi_out_dir = roi.out_dir if roi else inputs.output_dir
        mask_path   = roi.mask_path if roi else roi_out_dir / "segmentation_masks.npy"
        pixel_size_um = roi.pixel_size_um if roi else VISIUM_UM_PX

        if mask_path.exists():
            report(0.2, f"單 ROI 模式（{roi_name}）：生成多邊形…")
            logger.info(f"Loupe 匯出：單 ROI 模式（{roi_name}），從 MCseg v2 遮罩生成")
            roi_geo = mask_to_geojson(mask_path, pixel_size_um)
            poly_json_path = roi_out_dir / "cellpose_polygons.json"
            with open(poly_json_path, "w", encoding="utf-8") as f:
                json.dump(roi_geo, f)
        else:
            logger.warning(f"找不到 {roi_name} 的 segmentation_masks.npy，將不匯出多邊形層")

    if output_dir:
        out_dir = Path(output_dir)
    else:
        out_dir = export_dir / "loupe" if inputs.is_merged else roi_out_dir / "export_loupe"

    report(0.6, "寫出 Loupe 檔案…")
    exporter = LoupeExporter(
        poly_json_path=poly_json_path if (poly_json_path and poly_json_path.exists()) else None,
        whitelist_path=resolve_path(whitelist) if whitelist else None,
    )
    exporter.export(inputs.h5ad_path, out_dir)
    return ExportResult(output_dirs=[out_dir], is_merged=inputs.is_merged)


def _combine_roi_polygons(inputs) -> Path:
    """把各 ROI 的多邊形平移到全域座標後併成一份 GeoJSON，full_id 加上 ROI 前綴。"""
    logger.info("Loupe 匯出：合併模式，產生 combined_cellpose_polygons.json")
    all_features: list = []

    for roi in inputs.rois:
        if not roi.mask_path.exists():
            logger.warning(f"  [{roi.name}] 找不到 segmentation_masks.npy，跳過")
            continue

        roi_geo = mask_to_geojson(roi.mask_path, roi.pixel_size_um)
        roi_x_um = roi.cfg.get("x", 0) * roi.pixel_size_um
        roi_y_um = roi.cfg.get("y", 0) * roi.pixel_size_um
        for feat in roi_geo.get("features", []):
            orig_id = feat["properties"].get("full_id", "")
            feat["properties"]["full_id"] = f"{roi.name}__{orig_id}"
            shift_geojson_coords(feat, roi_x_um, roi_y_um)
            all_features.append(feat)

    poly_json_path = inputs.output_dir / "combined_cellpose_polygons.json"
    with open(poly_json_path, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": all_features}, f)
    return poly_json_path
