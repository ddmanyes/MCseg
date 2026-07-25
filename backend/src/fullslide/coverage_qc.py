"""
全片分割覆蓋率 QC
=================

回答一個「總數看不出來」的問題：**分割在哪些區域失敗了？**

動機（EP 紀錄 `99f6fad4`，真實 dpcp01）：全片分割在 Day0/Day3 區域只分出
**380 顆邊緣碎片**，而原始 8µm 資料該區有 **31,969 個 bins**。EP 已確認原始掃描
涵蓋完整 —— 是分割在該區域失敗。這件事**不會**反映在全片總命中率上
（1,545,019 bins / 35.3%）：一個區域整片漏掉，只讓總數低幾個百分點，看起來
像是正常的損耗。

做法：比較兩個彼此獨立的量

* **bin 密度**（有多少 `in_tissue` bins）→ 那裡「應該」有多少組織
* **細胞密度**（有多少細胞）→ 那裡「實際」分出多少細胞

`in_tissue == 1` 是 Space Ranger 的組織判定，與我們的分割完全獨立 —— 所以
「bins 很多、細胞很少」是**分割失敗**的證據，而非組織稀疏。

分析單位有兩層，**以切片層級為主要結論**：

| 層級 | 函式 | 用途 |
|------|------|------|
| 組織切片（連通元件） | `compute_section_coverage` | **主要** —— 抓整片缺口 |
| 固定網格 | `compute_coverage_qc` | 次要 —— 在有問題的切片內定位 |

為什麼不是只用網格：見下方 `SECTION_LOW_RATIO` 一節前的說明 —— dpcp01 實測
證明**固定網格抓不到目標失敗**，切片層級一次就分辨出來。

兩種失敗模式都要抓：

1. **整片漏掉**：`cells_per_1k_bins` 遠低於全片中位數
2. **只剩碎片**：細胞數不算少，但 `median_cell_area_px` 遠低於中位數
   —— 這正是 Day0/Day3 那 380 顆「邊緣碎片」的樣子

**遮罩可能是數 GB 的 memmap**，因此細胞統計一律分塊進行，整份不進 RAM。
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.fullslide.coverage")

# 網格邊長（px）。2048 ≈ 560µm，足以容納數千個 bins 與數百顆細胞 ——
# 太小則統計量雜訊過大，太大則失去定位失敗區域的能力。
DEFAULT_GRID_PX = 2048

# 低於「全片中位數 × 此比例」即標記。0.3 是刻意寬鬆的值：目標是抓出
# Day0/Day3 那種近乎整片漏掉的區域，不是挑出每一個略低於平均的格子。
DEFAULT_LOW_RATIO = 0.3

# bins 太少的格子（切片邊緣、組織孤島）統計量本身不可信，不納入中位數也不標記。
DEFAULT_MIN_BINS = 200

# ⚠️ 中位數相對門檻有一個結構性盲點：**若分割在全片都失敗，中位數本身就是失敗值，
# 於是一格都不會被標記，QC 反而回報「全部通過」。**（測試
# `test_globally_low_density_warns_even_when_nothing_flagged` 釘住此行為。）
#
# 因此另設一個絕對下限。**此值為實測校準，不是推導值** ——
# 初版用「2µm bin 面積 4µm²、~50µm² 細胞覆蓋 12 bins → 約 80 cells/1k bins」
# 推導出 15.0，實測後發現**錯了約 8 倍**：該推導預設每個 bin 都落在細胞內，
# 但實際上只有約 21% 的 bin 落在分割細胞內（正是本專案一路記載的 20–35% 涵蓋率）。
# 15.0 會對**完全正常**的切片發出「全片失敗」警告。
#
# 實測（dpcp01 全片，EP 遮罩 `segmentation_masks_fullslide.npy`）：
#   切片層級 6.0–13.8、網格層級中位數 8.8–10.5 cells/1k bins
# 取 3.0（約實測中位數的 1/3）作為下限 —— 用途是抓「整片都不對」，
# 不是評價分割品質。組織與流程差異很大，需要時由呼叫端覆寫。
DEFAULT_MIN_GLOBAL_DENSITY = 3.0

# `min_bins` 不可用固定值：一個 2048px 網格在 TIFF 空間可容納約 78,000 個 bin
# （bin 間距 = 2µm / 0.2732 µm/px ≈ 7.3 px），固定 200 相當於滿格的 0.26%，
# 等於幾乎不過濾。改為「佔理論滿格的比例」，隨網格大小自動縮放。
DEFAULT_MIN_BIN_FRACTION = 0.15

# 2µm bin 在影像空間的間距（px）。用於推算一個網格的理論滿格 bin 數。
# 預設以 Visium fullres 0.2737 µm/px 計；呼叫端可依實際影像 mpp 覆寫。
DEFAULT_BIN_PITCH_PX = 2.0 / VISIUM_UM_PX


def _grid_index(coord: np.ndarray, grid_px: int) -> np.ndarray:
    return coord // grid_px


# ─────────────────────────────────────────────────────────────
# 組織切片層級（**主要**分析單位）
# ─────────────────────────────────────────────────────────────
#
# ⚠️ 為什麼不是只用固定網格：dpcp01 實測（2026-07-25）證明**固定網格抓不到目標
# 失敗**。已知 Day0 區段整體只有 6.02 cells/1k bins（其他三段 11.2–13.8，
# 相對中位數 0.52×），但切成 2048/1024px 網格後，Day0 存活網格的密度中位數變成
# 10.0–10.6 —— 與全片中位數**完全相同**（1.00–1.01×），任何 `low_ratio` 都標記
# 不到半格。
#
# 原因有二：
#   1. 網格把一個切片切成數十塊，缺口被同段內的正常區域稀釋；
#   2. 網格會落在切片之間的空隙上，產生「幾乎沒有 bin 卻有細胞」的格子
#      （實測密度高達 5488/1k bins），把分布尾巴拉得極寬 —— p25/median 只有 0.34，
#      於是「中位數 × 0.3」這條線本來就穿過分布的正常部分，
#      標記率高達 16–28% 卻與真實缺口無關。
#
# 以**組織切片**為單位則一次就分辨出來。切片是生物上有意義的單位，而且天然
# 避開了上述兩個問題。固定網格保留下來，用途改為「在已知有問題的切片內定位」。

SECTION_RASTER_PX = 64      # 切片偵測用的佔用圖解析度（影像 px / 格）
MIN_SECTION_BINS = 5000     # 小於此 bin 數的連通塊視為雜點，不算一個切片

# 切片層級的標記門檻（相對中位數）。**實測校準值**，dpcp01 全片四個切片：
#
#   section 0 (Day3) 13.81/1k = 1.19×   section 1 (Day2) 12.11/1k = 1.04×
#   section 2 (Day1) 11.15/1k = 0.96×   section 3 (Day0)  6.01/1k = 0.52×  ← 已知缺口
#
# 正常切片與缺口切片之間有 0.52 ↔ 0.96 這個很寬的空隙，門檻取其中間值 0.7
# （而非剛好卡住 0.52 的 0.55）—— 對兩側都留足餘裕，避免過擬合單一樣本。
#
# ⚠️ 校準基礎僅**一張切片、四個切片段**。拿到更多樣本應重新檢視。
SECTION_LOW_RATIO = 0.7


def detect_tissue_sections(
    row: np.ndarray,
    col: np.ndarray,
    mask_shape: tuple[int, int],
    *,
    raster_px: int = SECTION_RASTER_PX,
    min_section_bins: int = MIN_SECTION_BINS,
    close_iter: int = 2,
) -> list[dict]:
    """
    由 bin 的落點找出彼此分離的組織切片（2D 連通元件）。

    做法：把 bin 落點打到 `raster_px` 解析度的佔用圖上 → 形態學閉運算把同一塊
    組織內的空隙補起來 → 連通元件標記 → 濾掉過小的雜點。

    比「x 方向直方圖找空隙」通用：切片若在 2D 上錯落排列（非單純左右並排），
    1D 投影會把它們黏成一塊。

    Returns
    -------
    list[dict]
        每個切片一項，含 `label`、`x0/x1/y0/y1`（影像座標）、`n_bins`、
        以及 `bin_index`（屬於該切片的 bin 在輸入陣列中的索引）。
    """
    from scipy import ndimage

    h, w = mask_shape
    nr, nc = int(np.ceil(h / raster_px)), int(np.ceil(w / raster_px))
    ry, rx = row // raster_px, col // raster_px
    valid = (ry >= 0) & (ry < nr) & (rx >= 0) & (rx < nc)

    occ = np.zeros((nr, nc), dtype=bool)
    occ[ry[valid], rx[valid]] = True
    occ = ndimage.binary_closing(occ, np.ones((3, 3), bool), iterations=close_iter)

    labels, n = ndimage.label(occ)
    if n == 0:
        return []

    # 每個 bin 屬於哪個連通元件
    bin_label = np.zeros(len(row), dtype=np.int32)
    bin_label[valid] = labels[ry[valid], rx[valid]]

    sections = []
    for lab in range(1, n + 1):
        idx = np.nonzero(bin_label == lab)[0]
        if len(idx) < min_section_bins:
            continue
        sections.append({
            "section": len(sections),
            "x0": int(col[idx].min()), "x1": int(col[idx].max()) + 1,
            "y0": int(row[idx].min()), "y1": int(row[idx].max()) + 1,
            "n_bins": int(len(idx)),
            "bin_index": idx,
        })
    # 由左至右、上而下排序，讓 section 編號穩定可引用
    sections.sort(key=lambda s: (s["x0"], s["y0"]))
    for i, s in enumerate(sections):
        s["section"] = i
    return sections


def compute_section_coverage(
    mask,
    row: np.ndarray,
    col: np.ndarray,
    *,
    low_ratio: float = SECTION_LOW_RATIO,
    block: int = 4096,
    **detect_kwargs,
) -> "pd.DataFrame":  # noqa: F821
    """
    以**組織切片**為單位比較 bin 密度與細胞密度。

    `low_ratio` 預設比網格版的 0.3 寬鬆，因為切片層級的離散度小得多 ——
    見 `SECTION_LOW_RATIO` 的校準說明。

    ⚠️ **切片數少時中位數本身脆弱**：dpcp01 只有 4 個切片，若其中 2 個失敗，
    中位數就會被拖到失敗值而使相對比較失效。這是 `min_global_density`
    絕對下限存在的理由 —— 兩個機制互補，不可只留一個。
    """
    import pandas as pd

    sections = detect_tissue_sections(row, col, mask.shape, **detect_kwargs)
    if not sections:
        return pd.DataFrame(columns=[
            "section", "x0", "y0", "x1", "y1", "n_bins", "n_cells",
            "cells_per_1k_bins", "median_cell_area_px", "flagged", "flag_reason",
        ])

    rows = []
    for s in sections:
        n_cells, areas = _cell_stats_in_box(
            mask, s["y0"], s["y1"], s["x0"], s["x1"], block=block
        )
        rows.append({
            "section": s["section"],
            "x0": s["x0"], "y0": s["y0"], "x1": s["x1"], "y1": s["y1"],
            "n_bins": s["n_bins"],
            "n_cells": n_cells,
            "cells_per_1k_bins": 1000.0 * n_cells / s["n_bins"],
            "median_cell_area_px": float(np.median(areas)) if len(areas) else 0.0,
        })

    df = pd.DataFrame(rows)
    med = float(df["cells_per_1k_bins"].median())
    df["flagged"] = df["cells_per_1k_bins"] < med * low_ratio
    df["flag_reason"] = np.where(
        df["flagged"],
        "切片細胞密度遠低於其他切片（分割可能在此切片失敗）", "",
    )
    df.attrs["summary"] = {
        "n_sections": int(len(df)),
        "n_flagged": int(df["flagged"].sum()),
        "median_cells_per_1k_bins": med,
        "threshold_cells_per_1k_bins": med * low_ratio,
        "low_ratio": low_ratio,
    }
    if df["flagged"].any():
        bad = ", ".join(
            f"section {r.section} @({r.x0},{r.y0}) {r.cells_per_1k_bins:.1f}/1k"
            for r in df[df["flagged"]].itertuples()
        )
        logger.warning(
            f"⚠️ 切片覆蓋率：{int(df['flagged'].sum())}/{len(df)} 個切片密度異常低"
            f"（中位數 {med:.1f}，門檻 {med*low_ratio:.1f}）—— {bad}"
        )
    else:
        logger.info(
            f"切片覆蓋率：{len(df)} 個切片皆正常，中位數 {med:.1f} cells/1k bins"
        )
    return df


def _cell_stats_in_box(mask, y0: int, y1: int, x0: int, x1: int, *, block: int = 4096):
    """分塊統計一個矩形範圍內的細胞數與面積（memmap 安全）。"""
    counts = np.zeros(1, dtype=np.int64)
    for by in range(y0, y1, block):
        blk = np.asarray(mask[by:min(by + block, y1), x0:x1])
        c = np.bincount(blk.ravel())
        if len(c) > len(counts):
            grown = np.zeros(len(c), dtype=np.int64)
            grown[: len(counts)] = counts
            counts = grown
        counts[: len(c)] += c
    labels = np.nonzero(counts)[0]
    labels = labels[labels > 0]
    return len(labels), counts[labels]


def compute_coverage_qc(
    mask,
    tp_path: str | Path,
    crop_y0: int = 0,
    crop_x0: int = 0,
    *,
    scale: tuple[float, float] = (1.0, 1.0),
    transform: np.ndarray | None = None,
    alignment=None,
    grid_px: int = DEFAULT_GRID_PX,
    min_bins: int | None = None,
    low_ratio: float = DEFAULT_LOW_RATIO,
    min_global_density: float = DEFAULT_MIN_GLOBAL_DENSITY,
    bin_pitch_px: float = DEFAULT_BIN_PITCH_PX,
) -> "pd.DataFrame":  # noqa: F821
    """
    逐網格比較 bin 密度與細胞密度，標記分割可能失敗的區域。

    Parameters
    ----------
    mask : np.ndarray or np.memmap
        (H, W) int 遮罩，座標原點為裁切左上角。可為 memmap（不會整份載入）。
    tp_path : str | Path
        `tissue_positions.parquet` 路徑。
    crop_y0, crop_x0 : int
        裁切左上角在來源影像座標系中的位置。
    scale, transform, alignment
        與 `bin_attribution` 相同語意 —— **必須傳入同一組值**，否則 QC 與實際
        計數會落在不同座標系，QC 反而變成雜訊來源。
    grid_px : int
        網格邊長（px）。
    min_bins : int
        低於此 bin 數的格子不納入中位數、也不標記（統計量不可信）。
    low_ratio : float
        `cells_per_1k_bins` 或 `median_cell_area_px` 低於中位數 × 此比例即標記。

    Returns
    -------
    pd.DataFrame
        一列一格，欄位見 `_EMPTY_COLUMNS`；已依 `cells_per_1k_bins` 升冪排序
        （最可疑的在最前面）。`.attrs["summary"]` 含全片摘要。
    """
    import pandas as pd

    from backend.src.fullslide.pipeline import map_bins_to_mask, read_tissue_bins

    h, w = mask.shape
    if min_bins is None:
        # 隨網格面積縮放：固定值在不同 grid_px 下意義天差地遠
        # （2048px 網格滿格約 78,000 bins，固定 200 等於幾乎不過濾）
        min_bins = max(int((grid_px / bin_pitch_px) ** 2 * DEFAULT_MIN_BIN_FRACTION), 1)

    tp = read_tissue_bins(tp_path)
    row, col, in_bounds, desc = map_bins_to_mask(
        tp, (h, w), crop_y0, crop_x0,
        scale=scale, transform=transform, alignment=alignment,
    )

    n_grid_y = int(np.ceil(h / grid_px))
    n_grid_x = int(np.ceil(w / grid_px))

    # ── 每格的 bin 數（只算落在遮罩範圍內的）
    gy = _grid_index(row[in_bounds], grid_px)
    gx = _grid_index(col[in_bounds], grid_px)
    flat = gy * n_grid_x + gx
    bins_per_grid = np.bincount(flat, minlength=n_grid_y * n_grid_x)

    # ── 每格的細胞統計（分塊，memmap 安全）
    cells_per_grid, area_sum_per_grid, areas_by_grid = _cell_stats_per_grid(
        mask, grid_px, n_grid_y, n_grid_x
    )

    rows = []
    for idx in range(n_grid_y * n_grid_x):
        n_bins = int(bins_per_grid[idx])
        n_cells = int(cells_per_grid[idx])
        gy_i, gx_i = divmod(idx, n_grid_x)
        areas = areas_by_grid.get(idx)
        rows.append({
            "grid_y": gy_i,
            "grid_x": gx_i,
            "x0": gx_i * grid_px + crop_x0,
            "y0": gy_i * grid_px + crop_y0,
            "n_bins": n_bins,
            "n_cells": n_cells,
            "cells_per_1k_bins": (1000.0 * n_cells / n_bins) if n_bins else 0.0,
            "median_cell_area_px": float(np.median(areas)) if areas is not None and len(areas) else 0.0,
            "mean_cell_area_px": (float(area_sum_per_grid[idx] / n_cells) if n_cells else 0.0),
        })

    df = pd.DataFrame(rows, columns=list(_EMPTY_COLUMNS))
    df = _flag_suspicious(df, min_bins=min_bins, low_ratio=low_ratio)

    df.attrs["summary"] = _summarise(df, desc, grid_px, min_bins, low_ratio)
    df.attrs["summary"]["min_global_density"] = min_global_density
    df.attrs["summary"]["global_density_ok"] = _check_global_density(
        df.attrs["summary"], min_global_density
    )
    _log_summary(df.attrs["summary"])
    return df.sort_values(
        ["flagged", "cells_per_1k_bins"], ascending=[False, True]
    ).reset_index(drop=True)


_EMPTY_COLUMNS = (
    "grid_y", "grid_x", "x0", "y0", "n_bins", "n_cells",
    "cells_per_1k_bins", "median_cell_area_px", "mean_cell_area_px",
)


def _cell_stats_per_grid(mask, grid_px: int, n_grid_y: int, n_grid_x: int):
    """
    逐網格統計細胞數與面積，分塊讀取（memmap 安全）。

    以**網格本身為分塊單位**：每格獨立 `np.unique`，因此橫跨網格邊界的細胞會在
    兩格各記一次。以 2048px 網格、約 17px 細胞徑估算，邊界細胞佔比約 3%，
    且對每格影響相當 —— 對「密度比」這個相對指標無實質影響，但絕對細胞數
    會略高於遮罩的 `max label`，引用時需留意。
    """
    n_cells = np.zeros(n_grid_y * n_grid_x, dtype=np.int64)
    area_sum = np.zeros(n_grid_y * n_grid_x, dtype=np.int64)
    areas_by_grid: dict[int, np.ndarray] = {}

    h, w = mask.shape
    for gy in range(n_grid_y):
        for gx in range(n_grid_x):
            y0, x0 = gy * grid_px, gx * grid_px
            block = np.asarray(mask[y0:min(y0 + grid_px, h), x0:min(x0 + grid_px, w)])
            counts = np.bincount(block.ravel())
            labels = np.nonzero(counts)[0]
            labels = labels[labels > 0]
            idx = gy * n_grid_x + gx
            n_cells[idx] = len(labels)
            if len(labels):
                areas = counts[labels]
                area_sum[idx] = int(areas.sum())
                areas_by_grid[idx] = areas
    return n_cells, area_sum, areas_by_grid


def _flag_suspicious(df, *, min_bins: int, low_ratio: float):
    """
    標記可疑格子，並記錄**為什麼**可疑。

    只用有足夠 bins 的格子算中位數 —— 讓切片邊緣的稀疏格子拉低基準，
    等於讓真正的失敗區域看起來正常。
    """
    considered = df["n_bins"] >= min_bins
    df["considered"] = considered

    if not considered.any():
        df["flagged"] = False
        df["flag_reason"] = ""
        df.attrs["thresholds"] = {}
        return df

    med_density = float(df.loc[considered, "cells_per_1k_bins"].median())
    med_area = float(df.loc[considered, "median_cell_area_px"].median())
    thr_density = med_density * low_ratio
    thr_area = med_area * low_ratio

    low_density = considered & (df["cells_per_1k_bins"] < thr_density)
    # 「只剩碎片」：細胞面積異常小。細胞數為 0 的格子已由 low_density 抓到，
    # 不重複標記（面積 0 並非碎片，而是完全沒有）。
    fragmented = considered & (df["n_cells"] > 0) & (df["median_cell_area_px"] < thr_area)

    df["flagged"] = low_density | fragmented
    df["flag_reason"] = np.select(
        [low_density & fragmented, low_density, fragmented],
        ["密度過低＋細胞碎片化", "細胞密度過低", "細胞碎片化（面積過小）"],
        default="",
    )
    df.attrs["thresholds"] = {
        "median_cells_per_1k_bins": med_density,
        "median_cell_area_px": med_area,
        "threshold_cells_per_1k_bins": thr_density,
        "threshold_cell_area_px": thr_area,
        "low_ratio": low_ratio,
        "min_bins": min_bins,
    }
    return df


def _summarise(df, desc: str, grid_px: int, min_bins: int, low_ratio: float) -> dict:
    considered = df["considered"]
    flagged = df["flagged"]
    return {
        "transform": desc,
        "grid_px": grid_px,
        "n_grid_total": int(len(df)),
        "n_grid_considered": int(considered.sum()),
        "n_grid_flagged": int(flagged.sum()),
        "flagged_bins": int(df.loc[flagged, "n_bins"].sum()),
        "total_bins_considered": int(df.loc[considered, "n_bins"].sum()),
        "min_bins": min_bins,
        "low_ratio": low_ratio,
        **df.attrs.get("thresholds", {}),
    }


def _check_global_density(s: dict, min_global_density: float) -> bool:
    """
    中位數本身是否合理。用來補中位數相對門檻的盲點：全片都失敗時沒有格子會
    「相對」偏低，逐格標記會全數通過。
    """
    if not s["n_grid_considered"]:
        return False
    return s.get("median_cells_per_1k_bins", 0.0) >= min_global_density


def _log_summary(s: dict) -> None:
    if not s["n_grid_considered"]:
        logger.warning(
            f"覆蓋率 QC：沒有任何網格達到 min_bins={s['min_bins']} —— "
            "可能是 bins 與遮罩座標系錯配（請檢查對位設定），或裁切窗格落在組織外。"
        )
        return

    if not s["global_density_ok"]:
        logger.warning(
            f"⚠️ 覆蓋率 QC：**全片**細胞密度中位數僅 "
            f"{s['median_cells_per_1k_bins']:.1f} cells/1k bins，低於下限 "
            f"{s['min_global_density']:.0f} —— 逐格標記在此情形下不可靠"
            "（中位數本身就是失敗值，沒有格子會「相對」偏低）。"
            "請先確認分割參數與對位設定，而非只看被標記的格子。"
        )

    if s["n_grid_flagged"]:
        frac = s["flagged_bins"] / max(s["total_bins_considered"], 1)
        logger.warning(
            f"⚠️ 覆蓋率 QC：{s['n_grid_flagged']}/{s['n_grid_considered']} 個網格被標記，"
            f"涵蓋 {s['flagged_bins']:,} bins（占已評估的 {frac:.1%}）。"
            f"中位數 {s['median_cells_per_1k_bins']:.1f} cells/1k bins，"
            f"標記門檻 {s['threshold_cells_per_1k_bins']:.1f}。"
        )
    else:
        logger.info(
            f"覆蓋率 QC：{s['n_grid_considered']} 個網格全部通過，"
            f"中位數 {s['median_cells_per_1k_bins']:.1f} cells/1k bins。"
        )


def run_coverage_qc_from_config(config: dict, **kwargs) -> tuple:
    """
    由 config 解析全片 QC 所需輸入並執行，回傳 `(sections_df, grid_df)`。

    **切片層級是主要結論**（實測證明網格層級抓不到目標失敗，見本檔上方說明）；
    網格層級保留為「在有問題的切片內定位」用。兩份分別寫出
    `coverage_qc_sections.parquet` 與 `coverage_qc.parquet`。

    重用 `resolve_full_count_inputs` —— 座標變換、裁切原點、對位來源一律與
    Stage 2 實際計數走同一條解析路徑，QC 才有意義。
    """
    from backend.src.fullslide.pipeline import (
        map_bins_to_mask,
        read_tissue_bins,
        resolve_full_count_inputs,
    )

    inputs, error = resolve_full_count_inputs(config)
    if inputs is None:
        raise ValueError(error)

    mask = np.load(str(inputs["mask_path"]), mmap_mode="r")
    crop_y0, crop_x0 = inputs["origin_xy"][1], inputs["origin_xy"][0]
    common = dict(scale=inputs["scale"], transform=inputs["transform"])

    grid_df = compute_coverage_qc(
        mask, inputs["tp_path"], crop_y0=crop_y0, crop_x0=crop_x0,
        **common, **kwargs,
    )
    grid_df.attrs["summary"]["transform_source"] = inputs["transform_source"]

    tp = read_tissue_bins(inputs["tp_path"])
    row, col, in_bounds, _ = map_bins_to_mask(
        tp, mask.shape, crop_y0, crop_x0, **common
    )
    sections_df = compute_section_coverage(mask, row[in_bounds], col[in_bounds])

    out_dir = Path(inputs["output_dir"])
    grid_df.to_parquet(str(out_dir / "coverage_qc.parquet"), index=False)
    sections_df.to_parquet(str(out_dir / "coverage_qc_sections.parquet"), index=False)
    logger.info("覆蓋率 QC 已寫出：coverage_qc_sections.parquet（主要）＋ coverage_qc.parquet")
    return sections_df, grid_df
