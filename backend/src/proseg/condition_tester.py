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
        metrics["cell_area_cv"] = 0.0

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
        
        import gc
        gc.collect()

        return result

    def _run_proseg_minimal(self, condition: dict, work_dir: Path, zarr_dir: Path) -> dict:
        """
        執行最小化的 Proseg 測試。
        改以 ProsegPipeline 原生呼叫，確保引入正規 nucleus mask 避免 MCMC 無引導爆炸記憶體。
        """
        from backend.src.proseg.pipeline import ProsegPipeline
        
        # 開發模式：無可執行檔則模擬
        proseg_bin = Path(self.proseg_bin)
        if not proseg_bin.exists():
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

        scale_um_px = self.config.get("proseg", {}).get("constants", {}).get("scale_um_px", 0.2645833)
        rois = self.config.get("rois", [])
        if rois:
            scale_um_px = rois[0].get("pixel_size_um", scale_um_px)

        # 自動探測外部 cyto_mask.npy (如果 roi 下有)
        cyto_npy = None
        output_dir_base = Path(self.config["paths"].get("output_dir", "results/analysis"))
        roi_base = output_dir_base / "roi"
        if roi_base.exists():
            for d in roi_base.iterdir():
                cm = d / "cyto_mask.npy"
                if cm.exists():
                    cyto_npy = str(cm)
                    logger.info(f"自動探測到外部細胞質遮罩：{cm}")
                    break

        pipeline = ProsegPipeline(
            zarr_path=str(zarr_dir),
            output_dir=str(work_dir),
            max_dist=condition.get("max_dist", 40.0),
            compactness=condition.get("compactness", 0.06),
            dilation_radius=condition.get("dilation", 20),
            samples=condition.get("samples", 200),
            burnin_samples=min(150, int(condition.get("samples", 200) * 0.5)),
            recorded_samples=condition.get("recorded", 50),
            coordinate_scale=scale_um_px,
            padding=50,  # 稍微給 padding
            nucleus_label_name="cellpose_nuclei",
            use_cyto_mask_from_zarr=True,       # 優先從 Zarr，失敗則退回 cyto_mask_path
            cyto_mask_path=cyto_npy,            # [新增] 自動探測的外部遮罩
            cyto_label_name="eosin_cyto",
            use_watershed=condition.get("watershed", True),
            enforce_connectivity=condition.get("connectivity", True)
        )
        
        pipeline.run_full_pipeline()
        
        out_counts = work_dir / "counts.csv.gz"
        out_cells = work_dir / "cells.csv"
        csv_path = work_dir / "transcripts_for_proseg.csv"
        
        if not out_cells.exists():
            raise FileNotFoundError(f"找不到 cells.csv：{out_cells}")

        cells_df = pd.read_csv(out_cells)
        n_cells = len(cells_df)

        try:
            import gzip as _gz
            from scipy.io import mmread
            with _gz.open(out_counts, "rt") as f:
                X = mmread(f).tocsr()
            median_genes  = float(np.median(np.diff(X.indptr)))
            median_counts = float(np.median(np.array(X.sum(axis=1)).flatten()))
        except Exception:
            median_genes  = 0.0
            median_counts = 0.0

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
            "cell_area_cv":      0.0,
        }


    # ──────────────────────────────────────────
    # 縮圖生成（HE + 細胞輪廓疊圖）
    # ──────────────────────────────────────────

    def _generate_thumbnail(self, condition: dict, cond_dir: Path, metrics: dict) -> None:
        """生成 HE + 細胞輪廓疊圖（以 ROI crop TIF 為背景），儲存至 cond_dir/preview.jpg"""
        try:
            import cv2

            # 讀取 ROI crop TIF（正確的局部 H&E）
            he_img = self._load_roi_crop_background()
            orig_h, orig_w = he_img.shape[:2]

            # 讀取 Proseg 真實多邊形（um 座標），轉換到原始 ROI pixel 空間
            polygons = self._load_proseg_polygons_px(cond_dir, orig_w, orig_h)

            # 等比例縮放 H&E 至顯示尺寸
            DISPLAY_W = 686
            DISPLAY_H = 398
            scale_x = DISPLAY_W / max(orig_w, 1)
            scale_y = DISPLAY_H / max(orig_h, 1)
            he_disp = cv2.resize(he_img.astype(np.uint8), (DISPLAY_W, DISPLAY_H))

            overlay = he_disp.copy()
            if polygons:
                for pts in polygons:
                    # 縮放多邊形座標
                    pts_scaled = pts.copy().astype(np.float32)
                    pts_scaled[:, 0] = pts_scaled[:, 0] * scale_x
                    pts_scaled[:, 1] = pts_scaled[:, 1] * scale_y
                    pts_int = pts_scaled.astype(np.int32)
                    cv2.polylines(overlay, [pts_int.reshape(-1, 1, 2)], True, (0, 230, 90), 1)
            else:
                # 無真實多邊形時合成圓形示意
                synth = self._synthesize_cell_contours(metrics, condition, DISPLAY_W, DISPLAY_H)
                for pts in synth:
                    cv2.polylines(overlay, [pts.reshape(-1, 1, 2)], True, (0, 200, 255), 1)

            # 文字標註
            n = metrics.get("n_cells", 0)
            g = metrics.get("median_genes", 0)
            label_str = condition.get("label", "")
            txt = f"{label_str}  Cells:{n}  Genes:{g:.0f}"
            cv2.putText(overlay, txt, (8, DISPLAY_H - 10),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 180), 1, cv2.LINE_AA)

            bgr = cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR)
            cv2.imwrite(str(cond_dir / "preview.jpg"), bgr, [cv2.IMWRITE_JPEG_QUALITY, 90])
            logger.debug(f"縮圖已儲存：{cond_dir / 'preview.jpg'}")
        except Exception as e:
            logger.warning(f"縮圖生成失敗（非致命）：{e}")
            import traceback
            logger.debug(traceback.format_exc())

    def _load_roi_crop_background(self) -> "np.ndarray":
        """載入 ROI crop TIF（如 results/analysis/roi/text/he_crop.tif），失敗回傳合成背景"""
        try:
            import tifffile
            output_dir = self.config.get("paths", {}).get("output_dir", "results/analysis")
            rois = self.config.get("rois", [])
            roi_name = rois[0].get("name", "text") if rois else "text"
            from backend.src.utils.config import resolve_path
            he_crop_path = resolve_path(output_dir) / "roi" / roi_name / "he_crop.tif"
            if not he_crop_path.exists():
                raise FileNotFoundError(str(he_crop_path))
            with tifffile.TiffFile(str(he_crop_path)) as tif:
                patch = tif.asarray()
            if patch.ndim == 3 and patch.shape[0] in (3, 4):  # (C, H, W) → (H, W, C)
                patch = np.moveaxis(patch[:3], 0, -1)
            elif patch.ndim == 3:
                patch = patch[:, :, :3]
            elif patch.ndim == 2:
                patch = np.stack([patch] * 3, axis=2)
            logger.debug(f"ROI crop loaded: {patch.shape} from {he_crop_path}")
            return patch.astype(np.uint8)
        except Exception as e:
            logger.debug(f"無法讀取 ROI crop：{e}，使用合成背景")
            return self._make_synthetic_he(686, 398)

    def _load_proseg_polygons_px(self, cond_dir: Path, img_w: int, img_h: int) -> list:
        """
        從 proseg_results.json (gzipped GeoJSON) 讀取細胞多邊形，
        將 um 座標轉換為 ROI 像素座標。
        """
        import gzip as _gz
        scale_um_px = self.config.get("proseg", {}).get("constants", {}).get("scale_um_px", 0.2737)
        rois = self.config.get("rois", [])
        if rois:
            scale_um_px = rois[0].get("pixel_size_um", scale_um_px)
        try:
            json_path = cond_dir / "proseg_results.json"
            if not json_path.exists():
                return []
            # 自動偵測是否為 gzip
            with open(json_path, "rb") as f:
                magic = f.read(2)
            opener = _gz.open if magic == b"\x1f\x8b" else open
            mode = "rb" if magic == b"\x1f\x8b" else "r"
            with opener(json_path, mode) as f:
                data = json.loads(f.read().decode("utf-8") if mode == "rb" else f.read())

            polygons_px = []
            features = data.get("features", [])
            for feat in features:
                geom = feat.get("geometry", {})
                ctype = geom.get("type", "")
                coords = geom.get("coordinates", [])
                if ctype == "MultiPolygon":
                    ring = coords[0][0]  # 取第一個多邊形的外環
                elif ctype == "Polygon":
                    ring = coords[0]
                else:
                    continue
                pts = np.array(ring, dtype=np.float32)
                if pts.ndim == 3:
                    pts = pts[0]
                # um → pixel (直接除以 pixel_size_um)
                pts_px = pts / scale_um_px
                pts_px = np.clip(pts_px, 0, [img_w - 1, img_h - 1]).astype(np.int32)
                if len(pts_px) >= 3:
                    polygons_px.append(pts_px)
            logger.debug(f"載入 {len(polygons_px)} 個真實多邊形 (scale={scale_um_px} um/px)")
            return polygons_px
        except Exception as e:
            logger.debug(f"無法讀取 Proseg 幾何：{e}")
            return []

    def _make_synthetic_he(self, W: int, H: int) -> "np.ndarray":
        """生成合成 H&E 風格背景（粉紫色調）"""
        import cv2
        rng = np.random.default_rng(42)
        base = rng.integers(215, 245, (H, W, 3), dtype=np.uint8)
        base[:, :, 0] = np.clip(base[:, :, 0].astype(int) - rng.integers(0, 30, (H, W)), 180, 255).astype(np.uint8)
        base[:, :, 2] = np.clip(base[:, :, 2].astype(int) + rng.integers(0, 25, (H, W)), 200, 255).astype(np.uint8)
        return cv2.GaussianBlur(base, (5, 5), 1.5)

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

        # 決定真實的 Zarr 目錄 (針對特定 ROI)
        if not roi_name:
            rois = self.config.get("rois", [])
            if not rois:
                raise ValueError("未定義配置檔案中的 ROI")
            roi_name = rois[0]["name"]
        
        roi_zarr_dir = self.zarr_dir / roi_name / "proseg_integrated.zarr"
        logger.info(f"使用 Zarr 測試目標：{roi_zarr_dir}")

        results = []
        completed = 0

        # 順序執行（避免 VRAM 競爭）
        for idx, condition in enumerate(conditions):
            result = self._run_single_condition(condition, idx, roi_zarr_dir, self.out_dir)
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
        """產生指標比較矩陣拼圖 (取代原本的散點圖)，提供直觀的 H&E + 輪廓網格視野"""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import matplotlib.image as mpimg

            ok_df = df[df.get("status", "ok") == "ok"].copy() if "status" in df.columns else df.copy()
            if ok_df.empty:
                logger.warning("沒有成功的條件結果可供繪圖。")
                return

            n_conds = len(ok_df)
            cols = min(4, n_conds) # 最多 4 欄
            rows = max(1, (n_conds + cols - 1) // cols)

            fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4), dpi=150)
            
            # 展平 axes 處理單一 row 的情況
            if n_conds == 1:
                axes = np.array([axes])
            axes = axes.flatten()

            for i, (idx, row) in enumerate(ok_df.iterrows()):
                ax = axes[i]
                cond_idx = row["condition_idx"]
                cond_dir = self.out_dir / f"cond_{cond_idx:02d}"
                preview_path = cond_dir / "preview.jpg"
                
                if preview_path.exists():
                    img = mpimg.imread(str(preview_path))
                    ax.imshow(img)
                else:
                    ax.text(0.5, 0.5, "預覽圖生成失敗", ha="center", va="center")
                    
                title = f"{row.get('label', '')}\nCells: {row.get('n_cells', 0)} | Genes: {row.get('median_genes', 0):.0f}"
                color = "green" if row.get("median_genes", 0) > 30 else "black"
                ax.set_title(title, fontsize=10, fontweight="bold", color=color)
                ax.axis("off")

            # 隱藏多餘的空白子圖
            for j in range(len(ok_df), len(axes)):
                axes[j].axis("off")

            plt.tight_layout()
            plot_path = self.out_dir / "condition_comparison.png"
            fig.savefig(str(plot_path), bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"已儲存比較圖：{plot_path}")

        except Exception as e:
            logger.warning(f"產生比較圖失敗（非致命）：{e}")
