"""Stage 4 匯出的輸入解析。

Xenium 與 Loupe 兩條匯出路徑共用同一套「輸入從哪來」的規則：
分析結果 h5ad 的尋找與路徑白名單驗證、單/多 ROI 模式判定、
每個 ROI 的輸出目錄與遮罩路徑。這裡是這套規則的唯一擁有者。

只收斂兩邊相同的部分——僅 Xenium 用得到的 he_crop.tif / adata_002um.h5ad
由呼叫端自己從 RoiInputs.out_dir 接，避免介面長出 for_xenium 之類的旗標。
"""
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from backend.src.utils.config import resolve_path
from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.export.inputs")

# 單/多 ROI 模式都會先掃 output_dir 根目錄的候選檔（由呼叫端指定優先序），
# 找不到才往 roi/{name}/ 底下找這兩個。
ROI_FALLBACK_CANDIDATES = ("umap_computed.h5ad", "qc_preprocessed.h5ad")


@dataclass(frozen=True)
class RoiInputs:
    """單一 ROI 的匯出輸入路徑。"""
    name: str
    out_dir: Path
    pixel_size_um: float
    mask_path: Path
    cfg: dict = field(repr=False, default_factory=dict)   # 原始 config 條目（x/y/width_px…）


@dataclass(frozen=True)
class ExportInputs:
    """一次匯出所需的全部輸入。"""
    h5ad_path: Path
    is_merged: bool
    active_roi: "Optional[str]"
    rois: list[RoiInputs]
    output_dir: Path

    def find_roi(self, name: str) -> "Optional[RoiInputs]":
        """依名稱取 ROI；找不到回 None（呼叫端自行決定退路）。"""
        return next((r for r in self.rois if r.name == name), None)


def resolve_export_inputs(
    config: dict,
    input_h5ad: str = "",
    candidates: Sequence[str] = ROI_FALLBACK_CANDIDATES,
) -> ExportInputs:
    """
    解析匯出輸入。

    Parameters
    ----------
    config : dict
        pipeline config（需要 paths.output_dir / paths.data_root 與 rois）
    input_h5ad : str
        使用者指定的 h5ad（空字串 = 自動尋找）。絕對路徑必須位於
        output_dir 或 data_root 底下，否則 ValueError。
    candidates : Sequence[str]
        自動尋找時 output_dir 根目錄的候選檔名優先序（Xenium 與 Loupe 各自不同）。

    Raises
    ------
    ValueError
        指定的絕對路徑不在允許目錄底下。
    FileNotFoundError
        指定的檔案不存在，或自動尋找不到任何分析結果。
    """
    paths = config.get("paths", {})
    rois_cfg = config.get("rois", [{}])

    output_dir = resolve_path(paths.get("output_dir", "results/analysis"))
    data_root = resolve_path(paths.get("data_root", "."))

    h5ad_path = _resolve_h5ad_path(input_h5ad, output_dir, data_root, rois_cfg, candidates)
    is_merged, active_roi = _inspect_export_mode(h5ad_path)

    rois: list[RoiInputs] = []
    for roi in rois_cfg:
        name = roi.get("name", "")
        if not name:
            continue
        roi_out_dir = output_dir / "roi" / name
        rois.append(RoiInputs(
            name=name,
            out_dir=roi_out_dir,
            pixel_size_um=float(roi.get("pixel_size_um", VISIUM_UM_PX)),
            mask_path=roi_out_dir / "segmentation_masks.npy",
            cfg=dict(roi),
        ))

    return ExportInputs(
        h5ad_path=h5ad_path,
        is_merged=is_merged,
        active_roi=active_roi,
        rois=rois,
        output_dir=output_dir,
    )


def _resolve_h5ad_path(
    input_h5ad: str,
    output_dir: Path,
    data_root: Path,
    rois_cfg: list,
    candidates: Sequence[str],
) -> Path:
    """找出要匯出的 h5ad；指定路徑優先，否則依候選清單自動尋找。"""
    if input_h5ad:
        p = Path(input_h5ad)
        if p.is_absolute():
            resolved = p.resolve()
            allowed = [output_dir.resolve(), data_root.resolve()]
            if not any(str(resolved).startswith(str(a)) for a in allowed):
                raise ValueError("input_h5ad 路徑必須位於 output_dir 或 data_root 底下")
            return resolved
        # 相對路徑：output_dir 與 data_root 兩處查找
        for base in (output_dir, data_root):
            candidate = base / p
            if candidate.exists():
                return candidate
        raise FileNotFoundError("找不到指定的 h5ad，請確認路徑正確")

    for name in candidates:
        p = output_dir / name
        if p.exists():
            return p

    # 新路徑：roi/{roi_name}/ 底下
    for roi in rois_cfg:
        roi_name = roi.get("name", "")
        if not roi_name:
            continue
        roi_dir = output_dir / "roi" / roi_name
        for name in ROI_FALLBACK_CANDIDATES:
            p = roi_dir / name
            if p.exists():
                return p

    raise FileNotFoundError(
        f"找不到分析結果 h5ad，請先執行 Stage 3 分析。\n搜尋位置：{output_dir}"
    )


def _inspect_export_mode(h5ad_path: Path) -> "tuple[bool, Optional[str]]":
    """以 backed 模式讀 h5ad，判斷是否為多 ROI 合併結果並取出 active_roi。"""
    import scanpy as sc

    adata_head = sc.read_h5ad(str(h5ad_path), backed="r")
    try:
        first_obs = adata_head.obs_names[0] if len(adata_head) > 0 else ""
        is_merged = "__" in first_obs
        active_roi = adata_head.uns.get("active_roi") if "active_roi" in adata_head.uns else None
    finally:
        del adata_head   # 釋放 backed file handle
    return is_merged, active_roi
