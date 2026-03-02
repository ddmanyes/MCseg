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

        return result

    def _run_proseg_minimal(self, condition: dict, work_dir: Path, zarr_dir: Path) -> dict:
        """
        執行最小化的 Proseg 測試（限制在小型 ROI）。
        若 proseg 二進位不存在則回傳模擬指標（供開發測試）。
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

        # 真實模式：需要 zarr 已建構
        if not zarr_dir.exists():
            raise FileNotFoundError(f"Zarr 尚未建構：{zarr_dir}")

        # 建立 Proseg 指令
        cmd = [
            str(proseg_bin),
            "--zarr", str(zarr_dir),
            "--output", str(work_dir),
            "--max-dist", str(condition["max_dist"]),
            "--compactness", str(condition["compactness"]),
            "--dilation", str(int(condition["dilation"])),
            "--samples", str(int(condition["samples"])),
            "--recorded-samples", str(int(condition["recorded"])),
        ]
        if condition.get("watershed", True):
            cmd.append("--use-watershed")
        if condition.get("connectivity", True):
            cmd.append("--enforce-connectivity")

        # 限制 ROI 大小（快速測試）
        cmd += ["--roi-width", str(self.test_roi_um), "--roi-height", str(self.test_roi_um)]

        logger.debug(f"執行：{' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"Proseg 失敗（exit {result.returncode}）：{result.stderr[:500]}")

        # 讀取結果 h5ad
        h5ad_path = work_dir / "processed_proseg_cyto.h5ad"
        if not h5ad_path.exists():
            # 嘗試尋找任意 h5ad
            h5ads = list(work_dir.glob("*.h5ad"))
            if not h5ads:
                raise FileNotFoundError("找不到 Proseg 輸出 h5ad")
            h5ad_path = h5ads[0]

        import anndata as ad
        adata = ad.read_h5ad(str(h5ad_path))
        import scanpy as sc
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        return _compute_metrics(adata)

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
