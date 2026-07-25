"""
全片分割覆蓋率 QC
=================

回答一個「總數看不出來」的問題：**分割在哪些區域失敗了？**

動機（EP 紀錄 `99f6fad4`，真實 dpcp01）：全片分割在 Day0/Day3 區域只分出
**380 顆邊緣碎片**，而原始 8µm 資料該區有 **31,969 個 bins**。EP 已確認原始掃描
涵蓋完整 —— 是分割在該區域失敗。這件事**不會**反映在全片總命中率上
（1,545,019 bins / 35.3%）：一個區域整片漏掉，只讓總數低幾個百分點，看起來
像是正常的損耗。

做法：把影像切成固定網格，逐格比較

* **bin 密度**（該格有多少 `in_tissue` bins）→ 那裡「應該」有多少組織
* **細胞密度**（該格有多少細胞）→ 那裡「實際」分出多少細胞

`in_tissue == 1` 是 Space Ranger 的組織判定，與我們的分割完全獨立 —— 所以
「bins 很多、細胞很少」是**分割失敗**的證據，而非組織稀疏。

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
# 因此另設一個絕對下限。推導：2µm bin 的面積為 4µm²，一顆 ~50µm² 的細胞約覆蓋
# 12 個 bin → 健康密度約 1000/12 ≈ 80 cells/1k bins。取 15 作為下限（約其 1/5），
# 刻意保守 —— 它的用途是抓「整片都不對」，不是評價分割品質。組織差異很大，
# 需要時由呼叫端覆寫。
DEFAULT_MIN_GLOBAL_DENSITY = 15.0


def _grid_index(coord: np.ndarray, grid_px: int) -> np.ndarray:
    return coord // grid_px


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
    min_bins: int = DEFAULT_MIN_BINS,
    low_ratio: float = DEFAULT_LOW_RATIO,
    min_global_density: float = DEFAULT_MIN_GLOBAL_DENSITY,
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


def run_coverage_qc_from_config(config: dict, **kwargs) -> "pd.DataFrame":  # noqa: F821
    """
    由 config 解析全片 QC 所需輸入並執行，結果寫出 `coverage_qc.parquet`。

    重用 `resolve_full_count_inputs` —— 座標變換、裁切原點、對位來源一律與
    Stage 2 實際計數走同一條解析路徑，QC 才有意義。
    """
    from backend.src.fullslide.pipeline import resolve_full_count_inputs

    inputs, error = resolve_full_count_inputs(config)
    if inputs is None:
        raise ValueError(error)

    mask = np.load(str(inputs["mask_path"]), mmap_mode="r")
    df = compute_coverage_qc(
        mask,
        inputs["tp_path"],
        crop_y0=inputs["origin_xy"][1],
        crop_x0=inputs["origin_xy"][0],
        scale=inputs["scale"],
        transform=inputs["transform"],
        **kwargs,
    )
    df.attrs["summary"]["transform_source"] = inputs["transform_source"]

    out_path = Path(inputs["output_dir"]) / "coverage_qc.parquet"
    df.to_parquet(str(out_path), index=False)
    logger.info(f"覆蓋率 QC 已寫出：{out_path.name}")
    return df
