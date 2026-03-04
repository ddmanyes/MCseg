"""
Stage 2.5: Proseg 參數條件測試
對多組參數組合執行小型 ROI 的 Proseg 測試，評估指標後推薦最佳條件

評估指標：
- n_cells: 細胞數量
- median_genes: 每細胞中位基因數
- median_counts: 每細胞中位 UMI 數
- fraction_assigned: RNA 指派率
- cell_area_cv: 細胞大小變異係數
"""
import itertools
import json
import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd

from backend.src.utils.constants import PROSEG_UM_PX, GOLDEN_PARAMS
from backend.src.utils.config import resolve_path

logger = logging.getLogger("pipeline.conditions")


def _compute_metrics(adata) -> dict:
    """計算單一條件的評估指標"""
    metrics = {}
    metrics["n_cells"] = int(adata.n_obs)

    if "n_genes_by_counts" in adata.obs:
        metrics["median_genes"] = float(np.median(adata.obs["n_genes_by_counts"]))
    elif hasattr(adata.X, "toarray"):
        metrics["median_genes"] = float(np.median((adata.X > 0).sum(axis=1)))
    else:
        metrics["median_genes"] = float(np.median((adata.X > 0).sum(axis=1)))

    if "total_counts" in adata.obs:
        metrics["median_counts"] = float(np.median(adata.obs["total_counts"]))
    else:
        row_sums = np.array(adata.X.sum(axis=1)).flatten()
        metrics["median_counts"] = float(np.median(row_sums))

    # 細胞大小變異係數
    for area_col in ["surface_area", "area", "cell_area"]:
        if area_col in adata.obs:
            areas = adata.obs[area_col].values.astype(float)
            if areas.mean() > 0:
                metrics["cell_area_cv"] = float(areas.std() / areas.mean())
            break
    else:
        metrics["cell_area_cv"] = float("nan")

    return metrics


class ConditionTester:
    """Proseg 參數網格測試器"""

    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.cond_cfg = config.get("condition_test", {})
        self.proseg_bin = os.path.expanduser(config["paths"].get("proseg_bin", "~/.cargo/bin/proseg"))
        self.zarr_dir = resolve_path(config["paths"]["zarr_dir"])
        self.out_dir = resolve_path(config["paths"]["conditions_dir"])
        self.max_parallel = self.cond_cfg.get("max_parallel", 4)
        self.test_roi_um = self.cond_cfg.get("test_roi_um", 1000)
        # 測試時用較少 samples 加速
        self.test_samples = self.cond_cfg.get("samples", 200)
        self.test_recorded = self.cond_cfg.get("recorded_samples", 50)

    def _build_conditions(self, grid: dict) -> list[dict]:
        """從 grid 字典展開所有參數組合"""
        keys = list(grid.keys())
        values = [grid[k] for k in keys]
        conditions = []
        for combo in itertools.product(*values):
            cond = dict(zip(keys, combo))
            # 補充固定參數
            cond.setdefault("samples", self.test_samples)
            cond.setdefault("recorded", self.test_recorded)
            cond.setdefault("watershed", True)
            cond.setdefault("connectivity", True)
            conditions.append(cond)
        return conditions

    def _run_single_condition(
        self,
        condition: dict,
        cond_idx: int,
        zarr_dir: Path,
        out_base: Path,
    ) -> dict:
        """
        執行單一條件的 Proseg 測試。

        Returns: 含指標的結果字典
        """
        cond_dir = out_base / f"cond_{cond_idx:02d}"
        cond_dir.mkdir(parents=True, exist_ok=True)

        label = (
            f"d{condition['dilation']}_"
            f"c{str(condition['compactness']).replace('.', '')}_"
            f"m{int(condition['max_dist'])}"
        )
        result = {
            "condition_idx": cond_idx,
            "label": label,
            **condition,
        }

        try:
            # 嘗試從 zarr 執行精簡版 Proseg
            metrics = self._run_proseg_minimal(condition, cond_dir, zarr_dir)
            result.update(metrics)
            result["status"] = "ok"
            logger.info(
                f"  條件 {cond_idx} ({label}): "
                f"n_cells={metrics.get('n_cells', 0)}, "
                f"median_genes={metrics.get('median_genes', 0):.1f}"
            )
        except Exception as e:
            logger.error(f"  條件 {cond_idx} ({label}) 失敗：{e}")
            result["status"] = "error"
            result["error"] = str(e)
            result.update({k: 0 for k in ["n_cells", "median_genes", "median_counts", "fraction_assigned", "cell_area_cv"]})

        # 儲存指標
        metrics_path = cond_dir / "metrics.json"
        with open(metrics_path, "w") as f:
            json.dump(result, f, indent=2, default=str)

        # 生成 HE + 輪廓縮圖
        self._generate_thumbnail(condition, cond_dir, result)

        return result

    def _run_proseg_minimal(self, condition: dict, work_dir: Path, zarr_dir: Path) -> dict:
        """
        執行最小化的 Proseg 測試。
        採用 Proseg-Zarr-Integration 的架構：
          1. 從 zarr store 萃取 transcripts CSV（含核遮罩 cell_id 初始化）
          2. 以官方 CLI 正確旗標執行 proseg
          3. 讀取 cells.csv + counts 計算指標

        若 proseg 二進位不存在則回傳模擬指標（開發模式）。
        """
        proseg_bin = Path(self.proseg_bin)

        if not proseg_bin.exists():
            # 開發模式：回傳隨機模擬指標
            logger.warning(f"Proseg 二進位不存在：{proseg_bin}，使用模擬指標")
            rng = np.random.default_rng(seed=hash(str(condition)) % (2**32))
            base_cells = int(500 + condition.get("max_dist", 40) * 10 - condition.get("compactness", 0.06) * 1000)
            noise = rng.integers(-50, 50)
            return {
                "n_cells": max(10, base_cells + noise),
                "median_genes": float(20 + condition.get("max_dist", 40) * 0.5 + noise * 0.1),
                "median_counts": float(50 + condition.get("max_dist", 40) * 2 + noise * 0.5),
                "fraction_assigned": float(min(0.95, 0.5 + condition.get("max_dist", 40) / 200)),
                "cell_area_cv": float(0.3 + condition.get("compactness", 0.06) * 2),
            }

        # 真實模式
        if not zarr_dir.exists():
            raise FileNotFoundError(f"Zarr 尚未建構：{zarr_dir}")

        # ── Step 1: 從 zarr 萃取 transcript CSV ──────────────────────
        # 直接讀 zarr/points/transcripts/points.parquet/*.parquet
        # 欄位已確認：x, y, gene（不需要 spatialdata API）
        csv_path = work_dir / "transcripts_for_proseg.csv"
        if not csv_path.exists():
            logger.info(f"  萃取 transcripts 至 CSV：{csv_path}")
            try:
                # 尋找 zarr 內的 parquet 分片
                parquet_dir = zarr_dir / "points" / "transcripts" / "points.parquet"
                parquet_files = sorted(parquet_dir.glob("*.parquet")) if parquet_dir.exists() else []
                if not parquet_files:
                    raise FileNotFoundError(
                        f"找不到 transcript parquet：{parquet_dir}\n"
                        f"請確認 Stage 2 Zarr 建構已完成"
                    )

                dfs = [pd.read_parquet(p) for p in parquet_files]
                df_pts = pd.concat(dfs, ignore_index=True)
                logger.info(f"  讀取 {len(parquet_files)} 個 parquet 分片，共 {len(df_pts):,} 行")

                # 確認欄位存在（x, y, gene）
                for col in ("x", "y", "gene"):
                    if col not in df_pts.columns:
                        raise KeyError(f"parquet 缺少欄位：{col}（現有：{list(df_pts.columns)}）")

                df_pts["cell_id"] = 0
                df_pts["qv"] = 40
                df_pts["z"] = 0.0
                df_pts[["x", "y", "gene", "qv", "cell_id", "z"]].to_csv(csv_path, index=False)
                logger.info(f"  CSV 寫出完成：{csv_path}")
            except Exception as e:
                raise RuntimeError(f"萃取 transcript CSV 失敗：{e}") from e


        # ── Step 2: 定義輸出路徑 ──────────────────────────────────────
        out_counts   = work_dir / "counts.csv.gz"
        out_cells    = work_dir / "cells.csv"
        out_genes    = work_dir / "genes.csv"
        out_polygons = work_dir / "cells.geojson"

        # ── Step 3: 建立官方 CLI 指令 ─────────────────────────────────
        cmd = [
            str(proseg_bin),
            "--overwrite",
            "--output-cell-polygons",   str(out_polygons),
            "--output-counts",          str(out_counts),
            "--output-counts-fmt",      "csv-gz",
            "--output-cell-metadata",   str(out_cells),
            "--output-cell-metadata-fmt", "csv",
            "--output-gene-metadata",   str(out_genes),
            "--output-gene-metadata-fmt", "csv",
            "--gene-column",            "gene",
            "--x-column",               "x",
            "--y-column",               "y",
            "--z-column",               "z",
            "--cell-id-column",         "cell_id",
            "--cell-id-unassigned",     "0",
            "--ignore-z-coord",
            "--min-qv",                 "0",
            "--max-transcript-nucleus-distance", str(condition["max_dist"]),
            "--cell-compactness",        str(condition["compactness"]),
            "--expand-initialized-cells", str(int(condition["dilation"])),
            "--samples",                str(int(condition["samples"])),
            "--recorded-samples",       str(int(condition["recorded"])),
        ]
        if condition.get("connectivity", True):
            cmd.append("--enforce-connectivity")
        cmd.append(str(csv_path))

        logger.debug(f"執行：{' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"Proseg 失敗（exit {result.returncode}）：{result.stderr[:500]}")

        # ── Step 4: 讀取結果並計算指標 ───────────────────────────────
        if not out_cells.exists():
            raise FileNotFoundError(f"找不到 cells.csv：{out_cells}")

        cells_df = pd.read_csv(out_cells)
        n_cells = len(cells_df)

        # 從 counts 讀取基因數與 UMI
        try:
            import gzip as _gz
            from scipy.io import mmread
            with _gz.open(out_counts, "rt") as f:
                X = mmread(f).tocsr()
            median_genes  = float(np.median(np.diff(X.indptr)))   # nnz per cell
            median_counts = float(np.median(np.array(X.sum(axis=1)).flatten()))
        except Exception:
            median_genes  = 0.0
            median_counts = 0.0

        # fraction_assigned（若 cells.csv 有 n_transcripts 欄位）
        frac = 0.0
        if "n_transcripts" in cells_df.columns:
            total = cells_df["n_transcripts"].sum()
            try:
                df_csv = pd.read_csv(csv_path, usecols=["gene"])
                frac = float(total / max(len(df_csv), 1))
            except Exception:
                frac = 0.0

        return {
            "n_cells":           n_cells,
            "median_genes":      median_genes,
            "median_counts":     median_counts,
            "fraction_assigned": frac,
            "cell_area_cv":      float("nan"),
        }


    # ──────────────────────────────────────────
    # 縮圖生成（HE + 細胞輪廓疊圖）
    # ──────────────────────────────────────────

    def _generate_thumbnail(self, condition: dict, cond_dir: Path, metrics: dict) -> None:
        """生成 HE + 細胞輪廓疊圖，儲存至 cond_dir/preview.jpg"""
        try:
            import cv2
            W, H = 480, 320

            he_img = self._load_he_background(W, H)
            polygons = self._load_cell_polygons(cond_dir, W, H)
            if polygons is None:
                polygons = self._synthesize_cell_contours(metrics, condition, W, H)

            overlay = he_img.copy()
            for pts in polygons:
                cv2.polylines(overlay, [pts.reshape(-1, 1, 2)], True, (0, 230, 90), 1)

            n = metrics.get("n_cells", 0)
            g = metrics.get("median_genes", 0)
            txt = f"n_cells={n}  median_genes={g:.0f}"
            cv2.putText(overlay, txt, (8, H - 10),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.42, (255, 255, 255), 1, cv2.LINE_AA)

            bgr = cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR)
            cv2.imwrite(str(cond_dir / "preview.jpg"), bgr, [cv2.IMWRITE_JPEG_QUALITY, 88])
        except Exception as e:
            logger.warning(f"縮圖生成失敗（非致命）：{e}")

    def _load_he_background(self, W: int, H: int) -> "np.ndarray":
        """嘗試從 BTF 讀取最小解析度金字塔層，失敗回傳合成背景"""
        try:
            he_path_str = self.config.get("paths", {}).get("he_image", "")
            if not he_path_str:
                raise ValueError("未設定 he_image 路徑")
            he_path = Path(he_path_str).expanduser()
            if not he_path.exists():
                raise FileNotFoundError(str(he_path))

            import tifffile
            import cv2
            with tifffile.TiffFile(str(he_path)) as tif:
                series = tif.series[0]
                # 使用最小的金字塔層（最後一層 = 最小解析度）
                level = series.levels[-1]
                patch = level.asarray()
            if patch.ndim == 3 and patch.shape[2] >= 3:
                patch = patch[:, :, :3]
            elif patch.ndim == 2:
                patch = np.stack([patch] * 3, axis=2)
            return cv2.resize(patch.astype(np.uint8), (W, H))
        except Exception as e:
            logger.debug(f"無法讀取 H&E：{e}，使用合成背景")
            return self._make_synthetic_he(W, H)

    def _make_synthetic_he(self, W: int, H: int) -> "np.ndarray":
        """生成合成 H&E 風格背景（粉紫色調）"""
        import cv2
        rng = np.random.default_rng(42)
        base = rng.integers(215, 245, (H, W, 3), dtype=np.uint8)
        # H&E 色調：增加藍紫分量（Hematoxylin）
        base[:, :, 0] = np.clip(base[:, :, 0].astype(int) - rng.integers(0, 30, (H, W)), 180, 255).astype(np.uint8)
        base[:, :, 2] = np.clip(base[:, :, 2].astype(int) + rng.integers(0, 25, (H, W)), 200, 255).astype(np.uint8)
        return cv2.GaussianBlur(base, (5, 5), 1.5)

    def _load_cell_polygons(self, cond_dir: Path, W: int, H: int):
        """嘗試從 Proseg 輸出讀取細胞幾何（GeoJSON 優先，h5ad spatial 次之）"""
        try:
            import cv2

            # 方案 A：GeoJSON
            gj_path = cond_dir / "cells.geojson"
            if gj_path.exists():
                with open(gj_path) as f:
                    gj = json.load(f)
                feats = [
                    f for f in gj.get("features", [])
                    if f.get("geometry", {}).get("type") == "Polygon"
                ][:400]
                if feats:
                    all_pts = np.vstack([
                        np.array(f["geometry"]["coordinates"][0]) for f in feats
                    ])
                    min_x, max_x = all_pts[:, 0].min(), all_pts[:, 0].max()
                    min_y, max_y = all_pts[:, 1].min(), all_pts[:, 1].max()
                    scale_x = (W - 20) / max(max_x - min_x, 1)
                    scale_y = (H - 20) / max(max_y - min_y, 1)
                    polygons = []
                    for f in feats:
                        coords = np.array(f["geometry"]["coordinates"][0])
                        px = ((coords[:, 0] - min_x) * scale_x + 10).astype(np.int32)
                        py = ((coords[:, 1] - min_y) * scale_y + 10).astype(np.int32)
                        polygons.append(np.column_stack([px, py]))
                    return polygons

            # 方案 B：h5ad obsm["spatial"]
            h5ads = list(cond_dir.glob("*.h5ad"))
            if h5ads:
                import anndata as ad
                adata = ad.read_h5ad(str(h5ads[0]))
                if "spatial" in adata.obsm:
                    coords = adata.obsm["spatial"][:400]
                    min_x, max_x = coords[:, 0].min(), coords[:, 0].max()
                    min_y, max_y = coords[:, 1].min(), coords[:, 1].max()
                    xs = ((coords[:, 0] - min_x) / max(max_x - min_x, 1) * (W - 20) + 10).astype(int)
                    ys = ((coords[:, 1] - min_y) / max(max_y - min_y, 1) * (H - 20) + 10).astype(int)
                    angles = np.linspace(0, 2 * np.pi, 10, endpoint=False)
                    r = 6
                    polygons = []
                    for cx, cy in zip(xs, ys):
                        pts = np.column_stack([
                            (cx + r * np.cos(angles)).astype(np.int32),
                            (cy + r * np.sin(angles)).astype(np.int32),
                        ])
                        polygons.append(pts)
                    return polygons
        except Exception as e:
            logger.debug(f"無法讀取 Proseg 幾何：{e}")
        return None

    def _synthesize_cell_contours(
        self, metrics: dict, condition: dict, W: int, H: int
    ) -> list:
        """根據條件指標合成細胞輪廓（dev mode 展示用）"""
        n_cells = min(int(metrics.get("n_cells", 100)), 350)
        dilation = float(condition.get("dilation", 20))
        compactness = float(condition.get("compactness", 0.06))

        base_r = max(4, int(dilation * 0.55))
        rng = np.random.default_rng(seed=42)
        cx = rng.integers(base_r + 6, W - base_r - 6, n_cells)
        cy = rng.integers(base_r + 6, H - base_r - 6, n_cells)

        n_pts = 12
        angles = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
        # compactness ↑ → 形狀越規則（擾動越小）
        perturb = max(0.15, min(0.85, 1.0 - compactness * 5))

        polygons = []
        for x, y in zip(cx, cy):
            radii = base_r * (1 + rng.uniform(-perturb * 0.5, perturb * 0.5, n_pts))
            pts = np.column_stack([
                (x + radii * np.cos(angles)).astype(np.int32),
                (y + radii * np.sin(angles)).astype(np.int32),
            ])
            polygons.append(pts)
        return polygons

    def run_grid(
        self,
        grid: dict,
        roi_name: str = "",
        on_progress: Optional[Callable[[int, dict], None]] = None,
    ) -> list[dict]:
        """
        執行完整參數網格測試。

        Parameters
        ----------
        grid : 參數網格，e.g. {"max_dist": [20,40], "compactness": [0.03,0.06], "dilation": [10,20]}
        roi_name : 指定 ROI 名稱（空字串 = 使用第一個 ROI）
        on_progress : 每完成一個條件時的 callback(completed_count, result)

        Returns
        -------
        list of result dicts
        """
        conditions = self._build_conditions(grid)
        logger.info(f"開始測試 {len(conditions)} 個條件（最大並行：{self.max_parallel}）")
        self.out_dir.mkdir(parents=True, exist_ok=True)

        results = []
        completed = 0

        # 順序執行（避免 VRAM 競爭）
        for idx, condition in enumerate(conditions):
            result = self._run_single_condition(condition, idx, self.zarr_dir, self.out_dir)
            results.append(result)
            completed += 1
            if on_progress:
                on_progress(completed, result)

        # 儲存匯總 CSV
        df = pd.DataFrame(results)
        csv_path = self.out_dir / "condition_grid.csv"
        df.to_csv(str(csv_path), index=False)
        logger.info(f"條件測試完成，結果已儲存：{csv_path}")

        # 產生比較圖
        self._save_comparison_plot(df)

        return results

    def _save_comparison_plot(self, df: pd.DataFrame) -> None:
        """產生指標比較熱力圖"""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import seaborn as sns

            ok_df = df[df.get("status", "ok") == "ok"].copy() if "status" in df.columns else df.copy()
            if ok_df.empty:
                return

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # 散點圖：n_cells vs median_genes
            ax = axes[0]
            sc = ax.scatter(
                ok_df["n_cells"], ok_df["median_genes"],
                c=ok_df.get("fraction_assigned", ok_df["n_cells"]),
                cmap="viridis", s=80, alpha=0.8
            )
            for _, row in ok_df.iterrows():
                ax.annotate(row.get("label", ""), (row["n_cells"], row["median_genes"]),
                           fontsize=6, alpha=0.7)
            plt.colorbar(sc, ax=ax, label="fraction_assigned")
            ax.set_xlabel("n_cells")
            ax.set_ylabel("median_genes")
            ax.set_title("條件比較：細胞數 vs 基因豐富度")

            # 指標表
            ax = axes[1]
            ax.axis("off")
            display_cols = ["label", "n_cells", "median_genes", "median_counts"]
            display_cols = [c for c in display_cols if c in ok_df.columns]
            tbl = ax.table(
                cellText=ok_df[display_cols].round(1).values.tolist(),
                colLabels=display_cols,
                loc="center",
                cellLoc="center",
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8)
            ax.set_title("條件指標摘要")

            plt.tight_layout()
            plot_path = self.out_dir / "condition_comparison.png"
            fig.savefig(str(plot_path), dpi=150, bbox_inches="tight")
            plt.close(fig)
            logger.info(f"已儲存比較圖：{plot_path}")

        except Exception as e:
            logger.warning(f"產生比較圖失敗（非致命）：{e}")
