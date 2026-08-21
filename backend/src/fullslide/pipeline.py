"""
全片流程核心函式（CLI 與 API 共用）
=====================================

本模組收斂「MCseg v2 遮罩 + Visium HD bins → cells×genes」這條配方，
使 `cli/segment.py`（指令列全片流程）與 `api/cellpose_count.py`（GUI 全圖計數）
共用同一份實作，避免邏輯分歧（CLAUDE.md §11 DRY）。

大多數函式為純函式（不寫 log、不做快取），快取與進度回報留給呼叫端：

| 函式 | 職責 |
|------|------|
| `bin_attribution`   | 2µm bins → cell_id 對應表 |
| `aggregate_cells`   | 依對應表把 bins 聚合成 cells×genes 原始 counts |
| `add_centroids`     | 補細胞重心（裁切局部 px、全片 fullres px、µm） |
| `resolve_pixel_size`| 取樣本實際 µm/px（scalefactors 優先於預設常數） |

`run_full_slide_segmentation` 是例外——它是「開影像 → 裁切 → 分割 → 存檔 →
寫 metadata」這整套全片分割編排的唯一擁有者（架構深化 P8，下沉自
`api/segmentation.py`），純同步、不依賴 asyncio/FastAPI，進度透過注入的
callback 回報，呼叫端（API 背景任務、CLI）自行決定要不要丟進 thread/executor。
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Optional

import numpy as np

logger = logging.getLogger("pipeline.fullslide")

if TYPE_CHECKING:  # pragma: no cover - 僅供型別檢查，避免匯入期拉進重量級套件
    import anndata as ad
    import pandas as pd


def resolve_bin_to_mask_scale(
    config: dict, mask_shape: tuple[int, int]
) -> tuple[float, float]:
    """
    推導 Space Ranger fullres px → 遮罩(raw TIFF) px 的**分軸**縮放。

    為何需要：Space Ranger 的 "fullres" 座標系是校準到**餵給它的那張影像**，
    未必等於 MSseg 分割所用的 raw TIFF。例如 dpcp01 樣本的 SR fullres 是
    0.5464 µm/px 而 TIFF 約 0.2737 → 相差近 2 倍。若不縮放，每個 bin 都會
    落在約一半的位置，**不報錯、不警告**（EP 2026-07-21 曾因此重跑計數）。

    推導不變量：`tissue_hires_image.png` 尺寸 ÷ `tissue_hires_scalef`
    ＝ SR fullres 畫布尺寸。此法算出的 (11266, 23552) 與 EP 已驗證的
    `segmentation_masks_fullslide_vfr_CORRECTED.npy` shape 完全一致。

    **必須分軸**：實測 dpcp01 為 col 1.9088、row 2.0000（差 4.6%）。強套等向
    2.0 會把最右側 bin 映射到 22484 px，超出 TIFF 寬度 21504 近千像素。

    TIFF 的 `XResolution` 標籤不可用作依據 —— 實測為 96 dpi 的通用預設值。

    Parameters
    ----------
    config : dict
        需含 `paths.binned_002`；`alignment.bin_to_mask_scale` 若指定 `[sx, sy]`
        則優先採用（供自動推導失準時人工覆寫）。
    mask_shape : tuple[int, int]
        遮罩的 `(height, width)`，即 raw TIFF 的像素尺寸。

    Returns
    -------
    tuple[float, float]
        `(scale_x, scale_y)`；無法推導時回退 `(1.0, 1.0)`。
    """
    import json

    override = (config.get("alignment") or {}).get("bin_to_mask_scale")
    if override:
        return float(override[0]), float(override[1])

    binned = config.get("paths", {}).get("binned_002", "")
    if not binned:
        return 1.0, 1.0

    spatial = Path(binned) / "spatial"
    sf_path = spatial / "scalefactors_json.json"
    hires_path = spatial / "tissue_hires_image.png"
    if not (sf_path.exists() and hires_path.exists()):
        return 1.0, 1.0

    try:
        scalef = float(json.loads(sf_path.read_text(encoding="utf-8"))["tissue_hires_scalef"])
        from PIL import Image

        Image.MAX_IMAGE_PIXELS = None   # hires 圖仍可能觸發 decompression bomb 警戒
        with Image.open(hires_path) as im:
            hires_w, hires_h = im.size
    except (ValueError, OSError, KeyError, TypeError):
        return 1.0, 1.0

    if scalef <= 0 or hires_w <= 0 or hires_h <= 0:
        return 1.0, 1.0

    sr_w, sr_h = hires_w / scalef, hires_h / scalef
    mask_h, mask_w = mask_shape
    return mask_w / sr_w, mask_h / sr_h


def resolve_bin_to_image_transform(
    config: dict, mask_shape: tuple[int, int]
) -> tuple[np.ndarray | None, tuple[float, float], str]:
    """
    決定「SR fullres px → 分割影像 px」的變換，回傳 `(transform, scale, source)`。

    優先序：

    1. **對位 JSON 組合**（幾何正確）：`inv(H_new) @ H_old`。
       `H_old` 由 `pick_source_alignment` 依 `scalefactors.microns_per_pixel` 選出；
       `H_new` 取 `alignment.extra_alignment_json`，未指定時從 `spatial/*.json`
       依「近似縮放反推的影像 mpp」自動比對。
    2. **近似分軸縮放**（`resolve_bin_to_mask_scale`）：找不到 JSON 時的回退。
       它把 bin 橫向壓縮進影像寬度，**在 SR 畫布相對影像有 padding 時幾何上是錯的**
       —— 回退時 log 會標明其為近似值。

    Returns
    -------
    tuple[np.ndarray | None, tuple[float, float], str]
        `transform` 為 3×3 或 None（此時採用 `scale`）；`source` 為人類可讀的來源說明。
    """
    scale = resolve_bin_to_mask_scale(config, mask_shape)
    align_cfg = config.get("alignment") or {}
    residual = _resolve_residual_alignment(align_cfg)

    if not align_cfg.get("use_alignment_json", True):
        base = np.diag([scale[0], scale[1], 1.0])
        if residual is None:
            return None, scale, f"近似縮放（設定停用對位 JSON）scale={scale[0]:.4f}, {scale[1]:.4f}"
        return (
            residual.to_3x3() @ base,
            (1.0, 1.0),
            f"近似縮放（停用對位 JSON）＋人工修正 scale={scale[0]:.4f}, {scale[1]:.4f}",
        )

    try:
        transform, source = _compose_from_alignment_json(config, scale)
    except (ValueError, OSError) as e:
        logger.warning(f"對位 JSON 不可用（{e}）→ 回退近似縮放")
        transform, source = None, ""

    if transform is None:
        logger.warning(
            f"未採用對位 JSON，改用近似縮放 scale=({scale[0]:.4f}, {scale[1]:.4f})。"
            "此值由 hires 尺寸 ÷ scalef 推導，**為近似值** —— SR 畫布相對影像有 padding "
            "時會把 bin 橫向壓縮。建議提供 Loupe 對位 JSON（alignment.extra_alignment_json）。"
        )
        if residual is None:
            return None, scale, f"近似縮放 scale={scale[0]:.4f}, {scale[1]:.4f}"
        transform = np.diag([scale[0], scale[1], 1.0])
        source = f"近似縮放 scale={scale[0]:.4f}, {scale[1]:.4f}"

    if residual is not None:
        transform = residual.to_3x3() @ transform
        tx, ty = residual.array[0, 2], residual.array[1, 2]
        logger.info(f"疊加人工/估計修正：平移 ({tx:+.1f}, {ty:+.1f}) px")
        source = f"{source} ＋人工修正"

    return transform, (1.0, 1.0) if residual is not None else scale, source


def resolve_dense_bin_to_image_transform(
    config: dict, mask_shape: tuple[int, int]
) -> tuple[np.ndarray, str]:
    """
    `resolve_bin_to_image_transform` 的稠密版本：一律回傳 3×3 矩陣。

    呼叫端若不需要區分「有變換」vs「只有近似縮放」這兩種語意（多數只想要
    一個能直接拿去乘的矩陣），改叫這個函式即可、不用自己再正規化一次
    `None → np.diag(scale)`。目前 `api/registration.py` 與 `roi/extractor.py`
    共用此函式（CLAUDE.md §11 DRY —— 過去這段正規化在多處各自重寫過）。
    """
    transform, scale, source = resolve_bin_to_image_transform(config, mask_shape)
    if transform is None:
        transform = np.diag([scale[0], scale[1], 1.0])
    return transform, source


def _resolve_residual_alignment(align_cfg: dict):
    """
    取 `alignment.matrix` 的殘餘修正；**未啟用或為單位矩陣時回 None**。

    預設 `enabled: false` —— 估計錯誤而被靜默套用，比不修正更糟。使用者必須
    看過 residual 後手動啟用（`/api/registration/apply` 只寫入不啟用）。
    """
    if not align_cfg.get("enabled", False):
        return None
    try:
        from backend.src.registration.align import AffineAlignment

        al = AffineAlignment.from_dict(align_cfg)
    except (ValueError, TypeError, KeyError) as e:
        logger.warning(f"alignment.matrix 無法解析（{e}），忽略人工修正")
        return None
    return None if al.is_identity() else al


def _compose_from_alignment_json(
    config: dict, approx_scale: tuple[float, float]
) -> tuple[np.ndarray | None, str]:
    """由 spatial/ 內（與設定指定）的對位 JSON 組出 bin→影像 homography。"""
    from backend.src.registration.alignment import (
        compose_bin_to_image,
        load_alignment,
        matrix_scale,
        pick_source_alignment,
        validate_pair,
    )

    binned = config.get("paths", {}).get("binned_002", "")
    if not binned:
        return None, ""
    spatial = Path(binned) / "spatial"

    align_cfg = config.get("alignment") or {}
    extra = align_cfg.get("extra_alignment_json")

    candidates = sorted(spatial.glob("*.json")) if spatial.exists() else []
    extra_path = Path(extra) if extra else None
    if extra_path is not None and not extra_path.exists():
        # 設定指向一個不存在的檔案時**退回自動偵測**，而不是讓整條 JSON 路徑陣亡。
        # （曾因 YAML 少一個空格而讓 `null#...` 被解析成字串，靜默關掉整個功能。）
        logger.warning(
            f"alignment.extra_alignment_json 指向的檔案不存在（{extra_path.name}），"
            "改為自動從 spatial/ 偵測"
        )
        extra_path = None
    if extra_path is not None and extra_path not in candidates:
        candidates.append(extra_path)
    if not candidates:
        return None, ""

    target_mpp = resolve_pixel_size(config)
    h_old, others = pick_source_alignment(candidates, target_mpp)

    if extra_path is not None:
        h_new = load_alignment(extra_path)
    else:
        # 自動偵測 H_new：近似縮放反推的影像 mpp ≈ mpp_SR / scale。
        # 容差放寬到 5%，因為 approx_scale 本身就是近似值 —— 這裡只需**辨識出是哪一份檔**，
        # 精確幾何由 JSON 本身提供。
        implied = h_old.mpp / float(np.sqrt(abs(approx_scale[0] * approx_scale[1])))
        h_new = None
        # h_old 自己也是候選 —— 分割影像可能就是產生 pxl_*_in_fullres 的那一張
        for al in [h_old, *others]:
            if abs(al.mpp - implied) / implied < 0.05:
                h_new = al
                break
        if h_new is None:
            if others:
                detail = "、".join(f"{al.name}→{al.mpp:.4f}" for al in others)
                logger.info(
                    f"無對位檔的 mpp 接近分割影像推估值 {implied:.4f}（候選：{detail}）"
                )
            return None, ""

    if h_new.path == h_old.path:
        # 分割影像就是產生 pxl_*_in_fullres 的那張 → 無需變換
        return np.eye(3), f"對位 JSON {h_old.name}（同一影像，變換為單位矩陣）"

    validate_pair(h_old, h_new)
    m = compose_bin_to_image(h_old, h_new)
    logger.info(
        f"採用對位 JSON 組合變換：{h_old.name}（mpp {h_old.mpp:.4f}）→ "
        f"{h_new.name}（mpp {h_new.mpp:.4f}），等效縮放 {matrix_scale(m):.5f}"
    )
    return m, f"對位 JSON {h_old.name} → {h_new.name}"


def read_tissue_bins(tp_path: str | Path) -> "pd.DataFrame":  # noqa: F821
    """讀 `tissue_positions.parquet` 並只留 `in_tissue == 1` 的 bins。"""
    import pandas as pd

    tp = pd.read_parquet(
        str(tp_path),
        columns=["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    return tp[tp["in_tissue"] == 1].copy()


def map_bins_to_mask(
    tp: "pd.DataFrame",  # noqa: F821
    mask_shape: tuple[int, int],
    crop_y0: int,
    crop_x0: int,
    *,
    scale: tuple[float, float] = (1.0, 1.0),
    transform: np.ndarray | None = None,
    alignment=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """
    把 bins 的 SR fullres 座標映射到遮罩局部 `(row, col)`。

    ```text
    SR fullres px ──transform / ×scale──> 影像(raw TIFF) px ──−crop origin──> 遮罩局部 px
    ```

    **先變換、後扣原點** —— 原點是遮罩空間的量。此函式由 `bin_attribution` 與
    覆蓋率 QC 共用；兩者若各自實作這段換算，遲早會漂移成不同的座標系
    （CLAUDE.md §11）。

    Returns
    -------
    tuple
        `(row, col, in_bounds, desc)`；`desc` 為人類可讀的變換來源說明。
    """
    src_row = tp["pxl_row_in_fullres"].values.astype(float)
    src_col = tp["pxl_col_in_fullres"].values.astype(float)

    if alignment is not None and not alignment.is_identity():
        # 殘餘修正疊在主變換之後 → 先把主變換也表達成 3×3 再左乘
        base = (
            np.asarray(transform, dtype=float) if transform is not None
            else np.diag([scale[0], scale[1], 1.0])
        )
        transform = alignment.to_3x3() @ base
        scale = (1.0, 1.0)   # 已併入 transform，避免重複套用

    if transform is not None:
        if scale != (1.0, 1.0):
            logger.warning(
                f"同時指定 transform 與 scale={scale}；採用 transform（幾何正確），忽略 scale。"
            )
        x, y = _apply_homography(np.asarray(transform, dtype=float), src_col, src_row)
        desc = "homography"
    else:
        x, y = src_col * scale[0], src_row * scale[1]
        desc = f"scale={scale}"

    row = np.rint(y - crop_y0).astype(np.int64)
    col = np.rint(x - crop_x0).astype(np.int64)

    h, w = mask_shape
    in_bounds = (row >= 0) & (row < h) & (col >= 0) & (col < w)
    return row, col, in_bounds, desc


def bin_attribution(
    mask: np.ndarray,
    tp_path: str | Path,
    crop_y0: int,
    crop_x0: int,
    out_path: str | Path | None = None,
    scale: tuple[float, float] = (1.0, 1.0),
    transform: np.ndarray | None = None,
    alignment=None,
) -> "pd.DataFrame":
    """
    將 Visium HD 2µm bins 對應到 MCseg v2 細胞遮罩。

    座標換算順序：**先變換、後扣原點**（原點是遮罩空間的量）。

    ```text
    SR fullres px ──transform / ×scale──> 影像(raw TIFF) px ──−crop origin──> 遮罩局部 px
    ```

    Parameters
    ----------
    mask : np.ndarray
        (H, W) int 遮罩，像素值 = 細胞 ID，0 為背景。座標原點為裁切左上角。
    tp_path : str | Path
        `tissue_positions.parquet` 路徑。
    crop_y0, crop_x0 : int
        裁切左上角在 raw TIFF 座標系中的位置。
    out_path : str | Path | None
        給定時額外寫出 parquet 快取。
    scale : tuple[float, float]
        `(scale_x, scale_y)`，SR fullres → 遮罩 px 的分軸縮放。
        由 `resolve_bin_to_mask_scale` 取得；`(1.0, 1.0)` 表示兩者同座標系。
    transform : np.ndarray | None
        3×3 homography（作用於 `(x, y) = (col, row)` 齊次座標），由
        `registration.alignment.compose_bin_to_image` 取得。這是**幾何正確**
        的路徑；`scale` 僅為找不到對位 JSON 時的近似回退。兩者同時給定時
        `transform` 優先並記 warning。
    alignment : AffineAlignment | None
        殘餘位移／仿射修正（`registration.align.AffineAlignment`），疊加在
        `transform`／`scale` **之後**（即影像 px 空間內的微調）。預設 None →
        行為完全不變。

    Returns
    -------
    pd.DataFrame
        欄位 `barcode`、`cell_id`；只含 `in_tissue == 1`、落在遮罩範圍內、
        且落在細胞內（`cell_id > 0`）的 bins。
        `.attrs["coverage"]` 記錄涵蓋率（見 `_coverage_attrs`）。

    Notes
    -----
    越界的 bin 會被**排除**而非夾到邊緣 —— 夾邊會把界外 bin 的 RNA 誤記到
    邊界細胞上，且不留痕跡。落在影像外是真實限制（dpcp01 正確變換下 SR 右緣
    約 980 px 寬確實不在高解析 TIFF 內），必須如實回報而非壓縮硬塞。
    """
    import pandas as pd

    tp = read_tissue_bins(tp_path)
    row, col, in_bounds, desc = map_bins_to_mask(
        tp, mask.shape, crop_y0, crop_x0,
        scale=scale, transform=transform, alignment=alignment,
    )
    h, w = mask.shape

    cell_ids = np.zeros(len(tp), dtype=mask.dtype)
    cell_ids[in_bounds] = mask[row[in_bounds], col[in_bounds]]
    tp["cell_id"] = cell_ids

    n_total = len(tp)
    n_oob = int((~in_bounds).sum())
    if n_oob and n_total and n_oob / n_total > 0.3 and _mask_covers_most_bins(
        row, col, (h, w)
    ):
        logger.warning(
            f"⚠️ {n_oob:,}/{n_total:,}（{n_oob/n_total:.1%}）個 bin 落在遮罩範圍外。"
            f"{desc}、crop=({crop_x0}, {crop_y0})、遮罩 {w}×{h}px —— "
            f"請確認 binned_outputs 是否為對應此影像的 CytAssist 註冊版本。"
        )

    attr = tp[tp["cell_id"] > 0][["barcode", "cell_id"]].reset_index(drop=True)
    attr.attrs["coverage"] = {
        "n_total": n_total,
        "n_in_bounds": int(in_bounds.sum()),
        "n_assigned": int(len(attr)),
        "frac_out_of_image": (n_oob / n_total) if n_total else 0.0,
    }
    _warn_if_hit_rate_implausible(mask, attr.attrs["coverage"], desc)

    if out_path is not None:
        attr.to_parquet(str(out_path), index=False)
    return attr


def _mask_covers_most_bins(
    row: np.ndarray, col: np.ndarray, mask_shape: tuple[int, int], min_frac: float = 0.5
) -> bool:
    """
    遮罩是否涵蓋 bins 的大部分範圍（＝這是整片，不是一個小裁切窗格）。

    「越界 bin 很多」只有在**整片**模式下才代表異常。裁一個 768×768 的窗格出來，
    全片 400 多萬個 bin 當然有 100% 落在窗格外 —— 對這種情形發出「請確認
    binned_outputs 是否為對應此影像的註冊版本」只會把人導去查沒問題的東西。
    """
    h, w = mask_shape
    span_y = float(row.max() - row.min()) if len(row) else 0.0
    span_x = float(col.max() - col.min()) if len(col) else 0.0
    if span_y <= 0 or span_x <= 0:
        return True
    return (h / span_y) >= min_frac and (w / span_x) >= min_frac


# 命中率（落在細胞內的 bin ÷ 界內的 bin）的合理下限。
#
# 這條檢查回答一個問題：**你分割的那張圖，是不是就是餵給 Space Ranger 的那張圖？**
# 是 → `pxl_*_in_fullres` 直接可用；否（例如另外掃了更高解析的圖）→ 需要縮放。
# 判斷錯了**不會報錯** —— 錯的座標仍讓 bin 乖乖落在影像界內，只是落在錯的細胞上。
#
# 門檻取自實測（2026-07-25，五個樣本）：
#   正確慣例：dpcp01 20.9%、SDS_D02 25.2%、SDS_D35 27.6%、CRC ROI 45–53%
#   錯誤慣例：0.8% – 3.9%
# 兩者差 6 倍以上，中間空隙很大，取 8%。
MIN_PLAUSIBLE_HIT_RATE = 0.08

# 低於此 bin 數不做判斷（統計量不可信，且小型測試遮罩會誤觸）
_HIT_RATE_MIN_BINS = 1000


def _warn_if_hit_rate_implausible(mask: np.ndarray, coverage: dict, desc: str) -> None:
    """
    命中率過低時警告，並區分「座標錯」與「分割沒產出細胞」。

    區分這兩者很重要：兩者的症狀相同（幾乎沒有 bin 拿到 cell_id），但要查的
    地方完全不同 —— 一個是去確認影像來源，一個是去看分割參數。
    """
    n_in = coverage["n_in_bounds"]
    if n_in < _HIT_RATE_MIN_BINS:
        return

    hit_rate = coverage["n_assigned"] / n_in
    coverage["hit_rate"] = hit_rate
    if hit_rate >= MIN_PLAUSIBLE_HIT_RATE:
        return

    # 遮罩本身有沒有細胞？只看 max()，memmap 上也便宜
    if int(mask.max()) == 0:
        logger.warning(
            f"⚠️ 遮罩內沒有任何細胞（max label = 0）—— 這是**分割**沒有產出結果，"
            f"不是座標問題。請檢查分割參數與輸入影像。"
        )
        return

    logger.warning(
        f"⚠️ 只有 {coverage['n_assigned']:,}/{n_in:,}（{hit_rate:.1%}）個界內 bin 落在細胞上，"
        f"遠低於正常值（實測 20–50%）。遮罩本身有細胞，所以這通常是**座標對不上**：\n"
        f"    請確認「拿來分割的影像」就是「餵給 Space Ranger 的影像」。\n"
        f"    若是另外掃描的高解析影像，需要縮放 —— "
        f"scale = scalefactors 的 microns_per_pixel ÷ 分割影像的 µm/px。\n"
        f"    目前採用的變換：{desc}"
    )


def _apply_homography(
    m: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """對 `(x, y)` 套用 3×3 homography（含齊次除法）。"""
    denom = m[2, 0] * x + m[2, 1] * y + m[2, 2]
    denom = np.where(denom == 0, np.nan, denom)   # 退化點 → NaN，後續判界時自然排除
    xt = (m[0, 0] * x + m[0, 1] * y + m[0, 2]) / denom
    yt = (m[1, 0] * x + m[1, 1] * y + m[1, 2]) / denom
    return np.nan_to_num(xt, nan=-1.0), np.nan_to_num(yt, nan=-1.0)


def aggregate_cells(attribution: "pd.DataFrame", h5_path: str | Path) -> "ad.AnnData":
    """
    依 attribution（barcode → cell_id）把 2µm bins 聚合成 cells×genes 原始 counts。

    不做 normalize —— 保留原始 counts 供下游自由運用。

    Parameters
    ----------
    attribution : pd.DataFrame
        `bin_attribution` 的輸出（欄位 `barcode`、`cell_id`）。
    h5_path : str | Path
        `filtered_feature_bc_matrix.h5` 路徑。

    Returns
    -------
    ad.AnnData
        X = 原始 counts（稀疏 CSR），`obs_names` 為 cell_id 字串，
        另含 `obs['cell_id']`（int）與 `obs['n_bins']`。
    """
    import gc

    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp

    adata_full = sc.read_10x_h5(str(h5_path))
    adata_full.var_names_make_unique()

    keep = adata_full.obs_names.isin(attribution["barcode"].values)
    adata_crop = adata_full[keep].copy()
    del adata_full
    gc.collect()

    barcode_to_cell = attribution.set_index("barcode")["cell_id"]
    cell_ids = barcode_to_cell.reindex(adata_crop.obs_names).values.astype(np.int32)
    valid = cell_ids > 0
    adata_valid = adata_crop[valid]
    cell_ids_v = cell_ids[valid]

    unique_cells = np.unique(cell_ids_v)
    n_cells = len(unique_cells)

    if n_cells == 0:
        # 沒有任何 bin 落在細胞上 —— 幾乎都是座標系錯配（bin_attribution 會先警告
        # 越界比例）。這裡必須給出可讀訊息，而不是讓下游對空陣列取 max 拋
        # `zero-size array to reduction operation maximum`。
        raise ValueError(
            "沒有任何 bin 對應到細胞，無法聚合 cells×genes。"
            "常見原因：bins 與分割影像座標系錯配（請檢查對位 JSON／裁切原點），"
            "或裁切窗格落在組織之外。"
        )

    # 向量化 LUT：O(max_id) 建立、O(n) 查詢，比 dict 快 10-100x
    lut = np.zeros(int(unique_cells.max()) + 1, dtype=np.int64)
    lut[unique_cells] = np.arange(n_cells)
    rows = lut[cell_ids_v]
    cols = np.arange(len(cell_ids_v))

    A = sp.csr_matrix(
        (np.ones(len(cell_ids_v), dtype=np.float32), (rows, cols)),
        shape=(n_cells, adata_valid.n_obs),
    )
    x_agg = A @ adata_valid.X

    cells = ad.AnnData(
        X=x_agg.tocsr() if sp.issparse(x_agg) else sp.csr_matrix(x_agg),
        var=adata_valid.var.copy(),
    )
    cells.obs_names = [str(int(c)) for c in unique_cells]
    cells.obs["cell_id"] = unique_cells.astype(int)
    cells.obs["n_bins"] = np.asarray(A.sum(axis=1)).ravel().astype(int)

    del adata_crop, adata_valid, A
    gc.collect()
    return cells


def add_centroids(
    cells: "ad.AnnData",
    mask: np.ndarray,
    pixel_size_um: float,
    origin_xy: tuple[int, int] = (0, 0),
) -> None:
    """
    就地補上細胞重心（三種座標系並存）。

    | 欄位 | 座標系 | 用途 |
    |------|--------|------|
    | `obs['centroid_x_px']` / `_y_px`           | 裁切局部 px | 遮罩內定位（維持 CLI 既有語意） |
    | `obs['centroid_x_fullres']` / `_y_fullres` | 全片 fullres px | 跨 ROI／全片框選（Stage 3.5 區域選取） |
    | `obsm['spatial']`                          | µm（由局部換算） | 匯出與繪圖（維持既有語意） |

    fullres 欄位在此處落地，而非留給下游每次讀 metadata sidecar 自行補償 ——
    後者是重複且易錯的轉換（座標系錯配會讓框選偏移一整個裁切原點）。

    Parameters
    ----------
    cells : ad.AnnData
        `aggregate_cells` 的輸出，需含 `obs['cell_id']`。就地修改。
    mask : np.ndarray
        (H, W) 遮罩，座標原點為裁切左上角。
    pixel_size_um : float
        µm/px（建議由 `resolve_pixel_size` 取得樣本實際值）。
    origin_xy : tuple[int, int]
        裁切左上角於原始影像 fullres 座標系的位置 `(x0, y0)`。
    """
    from scipy.ndimage import center_of_mass

    unique_cells = cells.obs["cell_id"].values.astype(np.int64)
    cen = np.asarray(
        center_of_mass(mask > 0, labels=mask, index=unique_cells.tolist()),
        dtype=float,
    )
    cy_px, cx_px = cen[:, 0], cen[:, 1]
    x0, y0 = origin_xy

    cells.obs["centroid_x_px"] = cx_px
    cells.obs["centroid_y_px"] = cy_px
    cells.obs["centroid_x_fullres"] = cx_px + x0
    cells.obs["centroid_y_fullres"] = cy_px + y0
    cells.obsm["spatial"] = np.stack([cx_px * pixel_size_um, cy_px * pixel_size_um], axis=1)


# ── 全圖遮罩 metadata sidecar ────────────────────────────────────────────────

FULL_SEG_MASK_FILENAME = "full_image_segmentation_masks.npy"
FULL_SEG_META_FILENAME = "full_image_segmentation_meta.json"


def write_full_seg_meta(
    output_dir: str | Path,
    *,
    crop_x0: int,
    crop_y0: int,
    width: int,
    height: int,
    n_cells: int,
    pixel_size_um: float,
    passes: int,
    image_width: int | None = None,
    image_height: int | None = None,
) -> Path:
    """
    寫出全圖遮罩的 metadata sidecar。

    遮罩本身只有局部座標，下游（Stage 2 全圖計數、Stage 3.5 框選）必須靠
    `crop_x0/crop_y0` 才能還原回原始影像 fullres 座標系。

    `image_width/height` 為**來源影像**的完整尺寸（非本次裁切窗格）。座標變換
    的推導以它為準 —— 用窗格尺寸會算出荒謬的縮放。
    """
    import json
    from datetime import datetime, timezone

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "crop_x0": int(crop_x0),
        "crop_y0": int(crop_y0),
        "width": int(width),
        "height": int(height),
        "n_cells": int(n_cells),
        "pixel_size_um": float(pixel_size_um),
        "passes": int(passes),
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    if image_width and image_height:
        meta["image_width"] = int(image_width)
        meta["image_height"] = int(image_height)
    path = output_dir / FULL_SEG_META_FILENAME
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return path


def read_full_seg_meta(output_dir: str | Path) -> dict | None:
    """讀取 metadata sidecar；不存在或損壞時回傳 None（由呼叫端決定如何處理）。"""
    import json

    path = Path(output_dir) / FULL_SEG_META_FILENAME
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return None


def _resolve_image_shape(
    config: dict, meta: dict | None, mask_shape: tuple[int, int]
) -> tuple[int, int]:
    """
    取**來源影像**的 `(height, width)`，供座標變換推導使用。

    取用順序：sidecar 記錄的影像尺寸 → 直接開 `paths.he_image` 讀 → 退回遮罩尺寸
    ＋裁切原點（僅在前兩者都不可得時；此時只是下限，會 log 警告）。
    """
    if meta and meta.get("image_width") and meta.get("image_height"):
        return int(meta["image_height"]), int(meta["image_width"])

    he = config.get("paths", {}).get("he_image", "")
    if he and Path(he).exists():
        try:
            from backend.src.utils.slide_reader import open_slide

            w, h = open_slide(he).dimensions
            return int(h), int(w)
        except (ValueError, OSError) as e:
            logger.warning(f"無法讀取 he_image 尺寸（{e}），改以遮罩尺寸推估")

    oy = int(meta["crop_y0"]) if meta else 0
    ox = int(meta["crop_x0"]) if meta else 0
    if ox or oy:
        logger.warning(
            "找不到來源影像，改以「裁切原點＋遮罩尺寸」當作影像尺寸 —— "
            "這只是下限，若遮罩並非影像右下角，座標變換會失準。"
        )
    return mask_shape[0] + oy, mask_shape[1] + ox


def resolve_full_count_inputs(config: dict) -> tuple[dict | None, str | None]:
    """
    解析全圖 RNA 計數所需的輸入，回傳 `(inputs, error)`。

    成功時 `inputs` 含 `mask_path`、`tp_path`、`h5_path`、`origin_xy`、
    `pixel_size_um`、`meta_missing`；失敗時 `inputs` 為 None 且 `error` 為
    可直接回傳給前端的訊息（不含絕對路徑）。

    `meta_missing=True` 代表找不到 sidecar，原點退為 (0, 0)。這是**警示而非
    錯誤** —— 舊版遮罩沒有 sidecar，且全圖模式的原點本來就是 (0, 0)。
    """
    from backend.src.utils.config import resolve_path

    paths = config.get("paths", {})
    output_dir = resolve_path(paths["output_dir"])
    mask_path = output_dir / FULL_SEG_MASK_FILENAME
    if not mask_path.exists():
        return None, "找不到全圖分割遮罩，請先完成全圖分割（Stage 1）"

    binned = paths.get("binned_002", "")
    tp_path = Path(binned) / "spatial" / "tissue_positions.parquet"
    if not tp_path.exists():
        return None, "找不到 tissue_positions.parquet，請確認 paths.binned_002 設定"
    h5_path = Path(binned) / "filtered_feature_bc_matrix.h5"
    if not h5_path.exists():
        return None, "找不到 filtered_feature_bc_matrix.h5，請確認 paths.binned_002 設定"

    meta = read_full_seg_meta(output_dir)
    origin_xy = (0, 0) if meta is None else (int(meta["crop_x0"]), int(meta["crop_y0"]))

    # 遮罩尺寸只讀 .npy header（mmap 不載入陣列本體，全片可達數 GB）
    mask_shape = np.load(str(mask_path), mmap_mode="r").shape
    # ⚠️ 變換要以**來源影像**的尺寸推導，不是遮罩的尺寸。遮罩可能只是影像的一個
    # 裁切窗格（Stage 1 的 crop 座標），拿它當畫布會算出荒謬的縮放
    # （實測 800×800 窗格 → scale (0.071, 0.034)，且連帶讓 JSON 的 H_new
    # 自動偵測失準）。遮罩在影像中的位置由 origin_xy 負責，與縮放無關。
    image_shape = _resolve_image_shape(config, meta, mask_shape)
    transform, scale, transform_source = resolve_bin_to_image_transform(config, image_shape)

    return {
        "mask_path": mask_path,
        "tp_path": tp_path,
        "h5_path": h5_path,
        "origin_xy": origin_xy,
        # transform 為主（幾何正確）；scale 僅在 transform 為 None 時生效
        "transform": transform,
        "scale": (1.0, 1.0) if transform is not None else scale,
        "transform_source": transform_source,
        # 與 ROI 路徑（counter.py）採同一設定，否則全圖與 ROI 結果不可比
        "dilation_px": int((config.get("rna_counting") or {}).get("dilation_px", 0)),
        "pixel_size_um": resolve_pixel_size(config),
        "meta_missing": meta is None,
        "output_dir": output_dir,
    }, None


def resolve_pixel_size(config: dict) -> float:
    """
    取樣本實際的 µm/px：`scalefactors_json.json` 的 `microns_per_pixel` 優先，
    缺失／損壞時回退 `constants.VISIUM_UM_PX`。

    注意：此值只影響 µm 換算與匯出比例尺，**不影響 bin attribution**
    （後者為純像素運算）。
    """
    import json

    from backend.src.utils.constants import VISIUM_UM_PX

    binned = config.get("paths", {}).get("binned_002", "")
    if binned:
        sf = Path(binned) / "spatial" / "scalefactors_json.json"
        if sf.exists():
            try:
                mpp = json.loads(sf.read_text(encoding="utf-8")).get("microns_per_pixel")
                if mpp:
                    return float(mpp)
            except (ValueError, OSError, TypeError):
                pass   # 損壞不可中斷流程，回退預設常數
    return float(VISIUM_UM_PX)


# ─── 全片分割編排（架構深化 P8，下沉自 api/segmentation.py）──────────────────

# 上界（crop_x1/crop_y1）專用哨兵：代表「取到影像邊界」，與 CLI 的
# `--crop-y1 -1` / `--btf-col1 -1` 語意一致。下界不接受 -1（視為負值錯誤，
# 由 API 層的 `validate_crop` 檢查——那是請求驗證職責，留在 API 層）。
_SENTINEL_FULL = -1


def _upper_is_open(v: Optional[int]) -> bool:
    """上界是否為開放（None 或 -1 → 取影像邊界）。"""
    return v is None or v == _SENTINEL_FULL


def resolve_crop_window(
    crop_x0: Optional[int], crop_y0: Optional[int],
    crop_x1: Optional[int], crop_y1: Optional[int],
    w_img: int, h_img: int,
) -> tuple[int, int, int, int]:
    """把裁切參數展開為實際切片邊界 `(x0, y0, x1, y1)`，並夾制於影像範圍內。

    吃 4 個 plain scalar（不是 pydantic 物件）——CLI 與 API 都能直接呼叫，
    不需要讓這個領域函式反向依賴 API 層的請求型別。
    """
    x0 = min(int(crop_x0 or 0), w_img)
    y0 = min(int(crop_y0 or 0), h_img)
    x1 = w_img if _upper_is_open(crop_x1) else min(int(crop_x1), w_img)
    y1 = h_img if _upper_is_open(crop_y1) else min(int(crop_y1), h_img)
    return x0, y0, x1, y1


def check_full_seg_mask_memory(width: int, height: int, max_load_gb: float) -> None:
    """檢查全片分割輸出遮罩（int32，寫檔前需完整存在於記憶體一次）是否超限。

    共用給 API 全片路徑（`run_full_slide_segmentation`）與 CLI 全片路徑
    （`cli.segment.step_segment`）——兩邊窗格大小的來源不同（前者用串流讀取
    的裁切窗格，後者用已載入的 `img.shape`），但同一份 int32 遮罩的記憶體
    上限判斷邏輯只該有一份。超限時 raise `MemoryError`，不得默默 OOM。
    """
    mask_gb = height * width * 4 / 1024 ** 3
    if mask_gb > max_load_gb:
        raise MemoryError(
            f"窗格 {width}×{height}px 的分割遮罩約 {mask_gb:.1f} GB（int32），"
            f"超過安全上限（max_load_gb = {max_load_gb:g} GB）。"
            f"請縮小裁切範圍，或調高該設定值。"
        )


def apply_full_seg_safety_clamp(
    seg_cfg: dict,
    *,
    force_disable_cpsam: bool = True,
    use_cpsam: Optional[bool] = None,
) -> dict:
    """套用 MPS `batch_size ≤ 2` 安全鉗制與 cpsam 停用決策，回傳一份新 cfg。

    batch_size 鉗制對所有呼叫端**無條件**生效（見 `docs/adr/0005`）。
    `force_disable_cpsam=True`（Web UI 全片路徑的預設，因為全圖模式 cpsam 極耗
    記憶體）時無視 `use_cpsam` 一律關閉；CLI 呼叫端應傳 `force_disable_cpsam=
    False`，讓使用者透過 `--cpsam` 旗標的既有行為保留（CLI 裁切窗格由使用者
    自行控制大小，不像 Web UI 全片按鈕預設面對整張未知大小的切片）。
    """
    seg_cfg_safe = dict(seg_cfg)
    seg_cfg_safe["batch_size"] = min(int(seg_cfg_safe.get("batch_size", 2)), 2)
    if force_disable_cpsam:
        seg_cfg_safe["use_cpsam"] = False
    elif use_cpsam is not None:
        seg_cfg_safe["use_cpsam"] = use_cpsam
    return seg_cfg_safe


@dataclass
class FullSegResult:
    """一次全片分割的產出。"""
    mask_path: Path
    n_cells: int
    image_width: int
    image_height: int
    crop_x0: int
    crop_y0: int
    is_full_image: bool


def run_full_slide_segmentation(
    config: dict,
    crop_x0: Optional[int] = None,
    crop_y0: Optional[int] = None,
    crop_x1: Optional[int] = None,
    crop_y1: Optional[int] = None,
    use_cpsam: Optional[bool] = None,
    progress: Optional[Callable[[float, str], None]] = None,
) -> FullSegResult:
    """全片 MCseg v2 分割：開影像 → 裁切窗格 → OOM 防護 → tiled 串流分割 → 存檔 + metadata。

    純同步、不依賴 asyncio/FastAPI；呼叫端（API 背景任務、CLI）自行決定要不要
    丟進 thread/executor。遇錯直接 raise（`FileNotFoundError`/`MemoryError`），
    由呼叫端接住並轉成各自的錯誤表示法（API 轉 `_full_status`，CLI 轉 exit code）。

    MPS `batch_size ≤ 2` 安全鉗制對 API 與 CLI **無條件**生效，即使明確請求更高
    的 batch_size 也會被降到 2（見 `docs/adr/0005-mps-batch-size-clamp-applies-to-cli.md`）。
    """
    import gc

    from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2
    from backend.src.utils.config import resolve_path
    from backend.src.utils.slide_reader import open_slide

    report = progress or (lambda p, msg: None)

    paths = config.get("paths", {})
    output_dir = resolve_path(paths["output_dir"])
    btf_path = paths.get("he_image", "")
    seg_cfg = config.get("segmentation", {}).get("mcseg_v2", {})
    full_cfg = config.get("full_seg", {})
    max_load_gb = float(full_cfg.get("max_load_gb", 6.0))

    if not btf_path or not Path(btf_path).exists():
        raise FileNotFoundError(f"找不到 BTF/TIFF：{btf_path}  請在 config paths.he_image 指定")

    report(0.02, "開啟影像（串流讀取）...")
    # SlideReader 統一 BTF/TIFF 與 NDPI/SVS；read_region 只解壓被請求的 tile，
    # 整張影像不進 RAM —— 這是全片（>6 GB）跑得動的前提。
    reader = open_slide(btf_path)
    w_img, h_img = reader.dimensions

    rx0, ry0, rx1, ry1 = resolve_crop_window(crop_x0, crop_y0, crop_x1, crop_y1, w_img, h_img)
    crop_w, crop_h = rx1 - rx0, ry1 - ry0

    def tile_reader(x, y, w, h, _r=reader, _ox=rx0, _oy=ry0):
        """裁切窗格內的相對座標 → 影像絕對座標。"""
        tile = _r.read_region(_ox + x, _oy + y, w, h)
        return tile[..., :3] if tile.ndim == 3 and tile.shape[-1] > 3 else tile

    # 影像本身已改為串流讀取，不再受 RAM 限制；剩下的硬限制是**輸出遮罩**
    # （int32，4 bytes/px）—— 它在寫檔前必須完整存在於記憶體一次。
    check_full_seg_mask_memory(crop_w, crop_h, max_load_gb)

    is_full_image = (rx0, ry0, rx1, ry1) == (0, 0, w_img, h_img)
    scope_label = "全圖" if is_full_image else f"窗格 ({rx0},{ry0})"
    report(0.05, f"{scope_label} 尺寸 {crop_w}×{crop_h}px，開始 tiled 串流分割...")

    # MPS 安全設定：tile_size=1024, batch_size≤2（無條件套用，docs/adr/0005）
    tile_size = int(full_cfg.get("tile_size", 1024))
    overlap = int(full_cfg.get("overlap", 128))
    seg_cfg_safe = apply_full_seg_safety_clamp(
        seg_cfg,
        force_disable_cpsam=full_cfg.get("force_disable_cpsam", True),
        use_cpsam=use_cpsam,
    )

    final_mask = run_tiled_mcseg_v2(
        cfg=seg_cfg_safe,
        tile_size=tile_size,
        overlap=overlap,
        progress_callback=report,
        tile_reader=tile_reader,
        full_shape=(crop_h, crop_w),
        # 標籤圖落在 memmap 並記錄進度：中斷後可自上次完成的 tile 列續跑
        work_dir=output_dir,
    )
    gc.collect()

    out_path = output_dir / "full_image_segmentation_masks.npy"
    output_dir.mkdir(parents=True, exist_ok=True)
    np.save(str(out_path), final_mask)
    n_cells = int(len(np.unique(final_mask)) - 1)

    # metadata sidecar：下游（Stage 2 全圖計數）需靠它把遮罩局部座標
    # 還原回原始影像 fullres 座標系
    write_full_seg_meta(
        output_dir,
        image_width=w_img,
        image_height=h_img,
        crop_x0=rx0,
        crop_y0=ry0,
        width=int(final_mask.shape[1]),
        height=int(final_mask.shape[0]),
        n_cells=n_cells,
        pixel_size_um=resolve_pixel_size(config),
        passes=7 if seg_cfg_safe.get("use_cpsam") else 4,
    )

    report(1.0, f"{scope_label}分割完成：{n_cells:,} 個細胞  →  {out_path.name}")

    return FullSegResult(
        mask_path=out_path,
        n_cells=n_cells,
        image_width=w_img,
        image_height=h_img,
        crop_x0=rx0,
        crop_y0=ry0,
        is_full_image=is_full_image,
    )
