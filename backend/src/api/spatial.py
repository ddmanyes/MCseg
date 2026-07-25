"""Stage 3.5：空間基因表達探索 API"""
import base64
import io
import logging
from pathlib import Path
from typing import Literal, Optional

import matplotlib
matplotlib.use("Agg")  # must be set before any other matplotlib import

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from backend.src.utils.config import load_config
from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.spatial")
router = APIRouter()

CMAPS = ["viridis", "magma", "plasma", "inferno", "Reds", "Blues", "YlOrRd"]


MIN_REGION_CELLS = 10


class RegionRequest(BaseModel):
    """框選區域（座標一律為**全片 fullres px**，與 RegionSelector 一致）。"""

    roi_name: Optional[str] = None
    region: Optional[dict] = None                  # {x0, y0, x1, y1}
    polygon: Optional[list[list[float]]] = None    # [[x, y], ...]


class GenePlotRequest(BaseModel):
    roi_name: Optional[str] = None
    genes: list[str]
    mode: Literal["dot", "contour", "set"] = "dot"
    set_name: Optional[str] = None   # Gene Set 模式的標題
    point_size: int = 6
    cmap: str = "viridis"
    alpha: float = 0.8
    dpi: int = 300
    # 框選區域（全片 fullres px）；None = 全部細胞
    region: Optional[dict] = None
    polygon: Optional[list[list[float]]] = None


# ── 路徑輔助 ────────────────────────────────────────────────────────────────

def _get_roi_h5ad(config: dict, roi_name: Optional[str]) -> Path:
    from backend.src.utils.config import resolve_path
    from backend.src.analysis.pipeline import _get_analysis_h5ad_dir
    output_dir = resolve_path(config["paths"]["output_dir"])
    h5ad_dir = _get_analysis_h5ad_dir(config, output_dir, roi_name)
    for candidate in [h5ad_dir / "qc_preprocessed.h5ad",
                       output_dir / "roi" / (roi_name or "") / "cellpose_cells.h5ad"]:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"找不到 h5ad（{roi_name}），請先完成 QC 步驟。")


def _get_he_crop(config: dict, roi_name: str) -> Path:
    from backend.src.utils.config import resolve_path
    output_dir = resolve_path(config["paths"]["output_dir"])
    path = output_dir / "roi" / roi_name / "he_crop.tif"
    if not path.exists():
        raise FileNotFoundError(f"找不到 H&E 裁切影像：{path}")
    return path


def _get_mask_path(config: dict, roi_name: str) -> Path:
    from backend.src.utils.config import resolve_path
    output_dir = resolve_path(config["paths"]["output_dir"])
    return output_dir / "roi" / roi_name / "segmentation_masks.npy"


def _get_pixel_size(config: dict, roi_name: Optional[str]) -> float:
    rois = config.get("rois", [])
    if roi_name:
        roi = next((r for r in rois if r.get("name") == roi_name), None)
        if roi:
            return float(roi.get("pixel_size_um", VISIUM_UM_PX))
    return VISIUM_UM_PX


def _local_spatial(adata, config: dict, roi_name: Optional[str]):
    import numpy as np
    spatial = adata.obsm["spatial"].copy()
    merge_mode = config.get("analysis", {}).get("merge_rois", False)
    if merge_mode and roi_name:
        rois = config.get("rois", [])
        roi_cfg = next((r for r in rois if r.get("name") == roi_name), None)
        if roi_cfg and roi_cfg.get("x") is not None:
            px_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))
            spatial[:, 0] -= float(roi_cfg["x"]) * px_um
            spatial[:, 1] -= float(roi_cfg["y"]) * px_um
    return spatial


def _fullres_xy(adata, config: dict, roi_name: Optional[str]):
    """
    取細胞質心的**全片 fullres px** 座標，三層回退：

    | 順序 | 來源 | 補償 |
    |------|------|------|
    | ① | `obs['centroid_x_fullres']` / `_y_fullres`（P0-7 全圖流程） | 無，已是全片座標 |
    | ② | `obs['centroid_x_px']` / `_y_px`（ROI 流程，裁切局部） | ＋ ROI 原點 |
    | ③ | `obsm['spatial']`（µm，舊資料） | ÷ pixel_size ＋ ROI 原點 |

    **為何必須統一**：RegionSelector 在全片座標上框選，而 ②③ 是裁切局部座標。
    直接比對會整批偏移一整個裁切原點（真實 ROI 可達數萬 px）—— 選到的是完全
    錯誤的區域，而且不報錯。
    """
    import numpy as np

    obs = adata.obs
    if "centroid_x_fullres" in obs.columns and "centroid_y_fullres" in obs.columns:
        return np.stack([
            obs["centroid_x_fullres"].values.astype(float),
            obs["centroid_y_fullres"].values.astype(float),
        ], axis=1)

    roi_cfg = next(
        (r for r in config.get("rois", []) if r.get("name") == roi_name), None
    ) if roi_name else None
    ox = float(roi_cfg["x"]) if roi_cfg and roi_cfg.get("x") is not None else 0.0
    oy = float(roi_cfg["y"]) if roi_cfg and roi_cfg.get("y") is not None else 0.0

    if "centroid_x_px" in obs.columns and "centroid_y_px" in obs.columns:
        return np.stack([
            obs["centroid_x_px"].values.astype(float) + ox,
            obs["centroid_y_px"].values.astype(float) + oy,
        ], axis=1)

    if "spatial" not in adata.obsm:
        raise ValueError("細胞資料缺少質心座標（centroid_*_fullres / centroid_*_px / obsm['spatial']）")

    px_um = _get_pixel_size(config, roi_name)
    return np.asarray(adata.obsm["spatial"], dtype=float)[:, :2] / px_um + [ox, oy]


def _filter_by_region(
    adata,
    region: Optional[dict],
    polygon: Optional[list],
    config: dict,
    roi_name: Optional[str],
):
    """
    以全片 fullres px 的 bbox 或多邊形過濾細胞，回傳布林遮罩。

    `region` 與 `polygon` 同時給定時 **polygon 優先**（較精確者勝）。
    兩者皆 None 時回傳全 True（未框選 = 全選）。邊界上的點計入。
    """
    import numpy as np

    if region is None and polygon is None:
        return np.ones(adata.n_obs, dtype=bool)

    pts = _fullres_xy(adata, config, roi_name)

    if polygon is not None:
        if len(polygon) < 3:
            logger.warning(f"多邊形頂點數不足（{len(polygon)}），視為未選取任何細胞")
            return np.zeros(adata.n_obs, dtype=bool)
        from matplotlib.path import Path as MplPath

        # radius 極小正值：讓正好落在邊上的點被判為「在內」
        return MplPath(np.asarray(polygon, dtype=float)).contains_points(pts, radius=1e-9)

    # 使用者可能從右下往左上拉框 → 正規化角點，否則會回傳空集合
    x0, x1 = sorted((float(region["x0"]), float(region["x1"])))
    y0, y1 = sorted((float(region["y0"]), float(region["y1"])))
    return (
        (pts[:, 0] >= x0) & (pts[:, 0] <= x1) &
        (pts[:, 1] >= y0) & (pts[:, 1] <= y1)
    )


def _get_expr(adata, genes: list[str]):
    """取得基因表現矩陣，回傳 (n_cells,) 陣列（多基因取平均 log1p）。"""
    import numpy as np
    from scipy.sparse import issparse

    stacked = []
    for gene in genes:
        x = adata[:, gene].X
        if issparse(x):
            x = x.toarray().flatten()
        else:
            x = np.asarray(x).flatten()
        stacked.append(np.log1p(x))
    return np.mean(stacked, axis=0)  # (n_cells,)


import re as _re
_CELL_ID_RE = _re.compile(r'^cell_(\d+)$')

def _parse_cell_id(name: str) -> Optional[int]:
    """從 obs_name 取出整數 cell ID。
    支援 "cell_42"（single-ROI）及 "roi1__cell_42"（merge）格式。
    非 cell_N 格式（如 barcode）回傳 None。
    """
    base = name.split("__")[-1]
    m = _CELL_ID_RE.match(base)
    return int(m.group(1)) if m else None


def _style_ax(ax, title: str):
    """統一圖表樣式。"""
    ax.set_facecolor("#1a1a2e")
    ax.set_title(title, color="white", fontsize=11, fontweight="bold")
    ax.set_xlabel("x (px)", color="gray", fontsize=7)
    ax.set_ylabel("y (px)", color="gray", fontsize=7)
    ax.tick_params(colors="gray", labelsize=6)
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")


def _add_colorbar(fig, ax, sc_plot):
    import matplotlib.pyplot as plt
    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log1p(counts)", color="white", fontsize=8)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")


# ── Contour 渲染核心 ─────────────────────────────────────────────────────────

def _render_contour(ax, he_img, mask, adata, expr: "np.ndarray", cmap: str, alpha: float):
    """Solid fill (semi-transparent) + dark inner border, H&E background."""
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from skimage.segmentation import find_boundaries

    if he_img is not None:
        ax.imshow(he_img, aspect="equal")

    # cell_id → expression LUT（支援 "cell_N" 及 "roi__cell_N" 兩種格式）
    parsed = [_parse_cell_id(n) for n in adata.obs_names]
    if any(cid is None for cid in parsed):
        raise ValueError(
            f"Contour 模式需要 obs_names 為 'cell_N' 格式，"
            f"但偵測到：'{adata.obs_names[0]}'。請改用 barcode 資料的 Dot 模式。"
        )
    cell_ids = np.array(parsed, dtype=int)
    max_id = int(mask.max())
    expr_lut = np.full(max_id + 1, np.nan)
    for i, cid in enumerate(cell_ids):
        if cid <= max_id:
            expr_lut[cid] = expr[i]

    valid_expr = expr_lut[~np.isnan(expr_lut)]
    vmax = max(float(np.nanpercentile(valid_expr, 95)) if valid_expr.size > 0 else 1.0, 1e-6)
    cmap_obj = plt.get_cmap(cmap)
    norm = mcolors.Normalize(vmin=0, vmax=vmax)

    H, W = mask.shape

    # ── Layer 1: semi-transparent expression fill for all cell pixels ──────
    fill = np.zeros((H, W, 4), dtype=np.float32)
    cy, cx = np.where(mask > 0)
    ccids = mask[cy, cx]
    valid = ccids <= max_id
    cy, cx, ccids = cy[valid], cx[valid], ccids[valid]
    cexpr = expr_lut[ccids]
    has_expr = ~np.isnan(cexpr)
    if has_expr.any():
        colors = cmap_obj(norm(cexpr[has_expr])).astype(np.float32)
        colors[:, 3] = alpha * 0.55
        fill[cy[has_expr], cx[has_expr]] = colors
    ax.imshow(fill, aspect="equal", interpolation="nearest")

    # ── Layer 2: dark inner border for cell separation ──────────────
    borders = find_boundaries(mask, mode="inner")
    border_overlay = np.zeros((H, W, 4), dtype=np.float32)
    border_overlay[borders] = [0.0, 0.0, 0.0, 0.85]
    ax.imshow(border_overlay, aspect="equal", interpolation="nearest")

    dummy = ax.scatter([], [], c=[], cmap=cmap, vmin=0, vmax=vmax)
    return dummy, vmax


# ── 單一 ROI 渲染到 ax（供 merge grid 迴圈使用）────────────────────────────

def _render_roi_ax(ax, fig, adata_roi, config: dict, roi_name: str,
                   genes: list[str], mode: str, cmap: str, alpha: float, point_size: int):
    """將單一 ROI 的表現量渲染到指定 ax。"""
    import numpy as np
    import tifffile

    pixel_size = _get_pixel_size(config, roi_name)
    spatial = _local_spatial(adata_roi, config, roi_name)
    x_px = spatial[:, 0] / pixel_size
    y_px = spatial[:, 1] / pixel_size

    try:
        he_img = tifffile.imread(str(_get_he_crop(config, roi_name)))
    except FileNotFoundError:
        he_img = None
    except Exception as e:
        logger.warning(f"H&E 載入失敗 ({roi_name}): {e}")
        he_img = None

    expr = _get_expr(adata_roi, genes)

    if mode in ("contour", "set"):
        mask_path = _get_mask_path(config, roi_name)
        if mask_path.exists():
            mask = np.load(str(mask_path))
            dummy, _ = _render_contour(ax, he_img, mask, adata_roi, expr, cmap, alpha)
            _add_colorbar(fig, ax, dummy)
            return

    # dot / scatter fallback
    import matplotlib.pyplot as plt
    if he_img is not None:
        ax.imshow(he_img, aspect="equal")
    import matplotlib.colors as mcolors
    vmax = max(float(np.percentile(expr[expr > 0], 95)) if (expr > 0).any() else 1.0, 1e-6)
    sc = ax.scatter(x_px, y_px, c=expr, cmap=cmap, s=point_size, alpha=alpha,
                    vmin=0, vmax=vmax, linewidths=0, rasterized=True)
    _add_colorbar(fig, ax, sc)


# ── API Endpoints ────────────────────────────────────────────────────────────

@router.get("/gene_list")
async def get_gene_list(roi_name: Optional[str] = None):
    config = load_config()
    try:
        import anndata
        merge_mode = config.get("analysis", {}).get("merge_rois", False)
        h5ad_path = _get_roi_h5ad(config, roi_name)
        adata = anndata.read_h5ad(h5ad_path, backed="r")
        genes = sorted(adata.var_names.tolist())
        has_roi_col = "roi" in adata.obs.columns
        adata.file.close()
        return {
            "status": "ok",
            "data": {
                "genes": genes,
                "n_genes": len(genes),
                "merge_mode": merge_mode and has_roi_col,
            },
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"gene_list 失敗：{e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/gene_plot")
async def post_gene_plot(req: GenePlotRequest):
    if not req.genes:
        raise HTTPException(status_code=400, detail="請至少選擇一個基因")
    if req.mode in ("dot", "contour") and len(req.genes) > 4:
        raise HTTPException(status_code=400, detail="dot/contour 模式最多 4 個基因")
    if req.cmap not in CMAPS:
        raise HTTPException(status_code=400, detail=f"不支援的 colormap，可用：{CMAPS}")

    config = load_config()
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        import tifffile
        import scanpy as sc

        h5ad_path = _get_roi_h5ad(config, req.roi_name)
        adata = sc.read_h5ad(h5ad_path)

        invalid = [g for g in req.genes if g not in adata.var_names]
        if invalid:
            raise HTTPException(status_code=400, detail=f"基因不存在：{invalid}")

        if req.region is not None or req.polygon is not None:
            sel = _filter_by_region(adata, req.region, req.polygon, config, req.roi_name)
            n_sel = int(sel.sum())
            # 細胞太少時做出的圖只會誤導（單一細胞的表現不代表任何空間結構）
            if n_sel < MIN_REGION_CELLS:
                raise HTTPException(
                    status_code=400,
                    detail=f"選取區域細胞數不足（{n_sel} < {MIN_REGION_CELLS}），請放大選取範圍",
                )
            logger.info(f"框選區域：{n_sel}/{adata.n_obs} 個細胞")
            adata = adata[sel].copy()

        merge_mode = config.get("analysis", {}).get("merge_rois", False)
        is_merge = merge_mode and "roi" in adata.obs.columns

        # ════════════════════════════════════════════════════════════════
        # Merge ROI mode：依 obs["roi"] 分組，各畫一格
        # ════════════════════════════════════════════════════════════════
        if is_merge:
            roi_names = sorted(adata.obs["roi"].unique().tolist())
            n_rois = len(roi_names)
            n_panels = n_rois if req.mode == "set" else n_rois * len(req.genes)
            if n_panels > 12:
                raise HTTPException(status_code=400,
                    detail=f"Grid 太大（{n_panels} 格），請減少基因數或改用 Gene Set 模式")

            if req.mode == "set":
                # 一個基因集 → n_rois 個 subplot
                ncols = min(n_rois, 3)
                nrows = (n_rois + ncols - 1) // ncols
                fig, axes = plt.subplots(nrows, ncols,
                                         figsize=(7 * ncols, 5.5 * nrows),
                                         squeeze=False, facecolor="#1a1a2e")
                fig.patch.set_facecolor("#1a1a2e")
                set_title = req.set_name or " + ".join(req.genes)

                for i, rname in enumerate(roi_names):
                    ax = axes[i // ncols][i % ncols]
                    adata_roi = adata[adata.obs["roi"] == rname]
                    _render_roi_ax(ax, fig, adata_roi, config, rname,
                                   req.genes, req.mode, req.cmap, req.alpha, req.point_size)
                    _style_ax(ax, f"{set_title}\n{rname}")

                for j in range(n_rois, nrows * ncols):
                    axes[j // ncols][j % ncols].set_visible(False)

            else:
                # dot / contour → rows = genes, cols = ROIs
                n_genes = len(req.genes)
                fig, axes = plt.subplots(n_genes, n_rois,
                                         figsize=(7 * n_rois, 5.5 * n_genes),
                                         squeeze=False, facecolor="#1a1a2e")
                fig.patch.set_facecolor("#1a1a2e")

                for gi, gene in enumerate(req.genes):
                    for ri, rname in enumerate(roi_names):
                        ax = axes[gi][ri]
                        adata_roi = adata[adata.obs["roi"] == rname]
                        _render_roi_ax(ax, fig, adata_roi, config, rname,
                                       [gene], req.mode, req.cmap, req.alpha, req.point_size)
                        _style_ax(ax, f"{gene}\n{rname}")

        # ════════════════════════════════════════════════════════════════
        # Single ROI mode
        # ════════════════════════════════════════════════════════════════
        else:
            spatial = _local_spatial(adata, config, req.roi_name)
            pixel_size = _get_pixel_size(config, req.roi_name)
            x_px = spatial[:, 0] / pixel_size
            y_px = spatial[:, 1] / pixel_size

            he_img = None
            if req.roi_name:
                try:
                    he_img = tifffile.imread(str(_get_he_crop(config, req.roi_name)))
                except FileNotFoundError:
                    pass

            if req.mode == "set":
                expr = _get_expr(adata, req.genes)
                set_title = req.set_name or " + ".join(req.genes)
                gene_label = " | ".join(req.genes)

                fig, ax = plt.subplots(1, 1, figsize=(8, 6), facecolor="#1a1a2e")
                fig.patch.set_facecolor("#1a1a2e")

                mask_path = _get_mask_path(config, req.roi_name or "") if req.roi_name else None
                if mask_path and mask_path.exists():
                    mask = np.load(str(mask_path))
                    dummy, _ = _render_contour(ax, he_img, mask, adata, expr, req.cmap, req.alpha)
                    _add_colorbar(fig, ax, dummy)
                else:
                    vmax = max(float(np.percentile(expr[expr > 0], 95)) if (expr > 0).any() else 1.0, 1e-6)
                    if he_img is not None:
                        ax.imshow(he_img, aspect="equal")
                    sc_plot = ax.scatter(x_px, y_px, c=expr, cmap=req.cmap,
                                         s=req.point_size, alpha=req.alpha,
                                         vmin=0, vmax=vmax, linewidths=0, rasterized=True)
                    _add_colorbar(fig, ax, sc_plot)
                _style_ax(ax, set_title)
                ax.set_xlabel(f"Genes: {gene_label}", color="gray", fontsize=7)

            elif req.mode == "dot":
                n = len(req.genes)
                ncols = min(n, 2)
                nrows = (n + 1) // 2
                fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows),
                                         squeeze=False, facecolor="#1a1a2e")
                fig.patch.set_facecolor("#1a1a2e")

                for i, gene in enumerate(req.genes):
                    ax = axes[i // ncols][i % ncols]
                    if he_img is not None:
                        ax.imshow(he_img, aspect="equal")
                    expr = _get_expr(adata, [gene])
                    vmax = max(float(np.percentile(expr[expr > 0], 95)) if (expr > 0).any() else 1.0, 1e-6)
                    sc_plot = ax.scatter(x_px, y_px, c=expr, cmap=req.cmap,
                                         s=req.point_size, alpha=req.alpha,
                                         vmin=0, vmax=vmax, linewidths=0, rasterized=True)
                    _add_colorbar(fig, ax, sc_plot)
                    _style_ax(ax, gene)

                for j in range(n, nrows * ncols):
                    axes[j // ncols][j % ncols].set_visible(False)

            else:  # contour
                mask_path = _get_mask_path(config, req.roi_name or "")
                if not mask_path.exists():
                    raise HTTPException(status_code=404, detail="找不到 segmentation_masks.npy，請先執行分割步驟")
                mask = np.load(str(mask_path))

                n = len(req.genes)
                ncols = min(n, 2)
                nrows = (n + 1) // 2
                fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5.5 * nrows),
                                         squeeze=False, facecolor="#1a1a2e")
                fig.patch.set_facecolor("#1a1a2e")

                for i, gene in enumerate(req.genes):
                    ax = axes[i // ncols][i % ncols]
                    expr = _get_expr(adata, [gene])
                    dummy, _ = _render_contour(ax, he_img, mask, adata, expr, req.cmap, req.alpha)
                    _add_colorbar(fig, ax, dummy)
                    _style_ax(ax, gene)

                for j in range(n, nrows * ncols):
                    axes[j // ncols][j % ncols].set_visible(False)

        plt.tight_layout(pad=1.2)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=req.dpi, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        plt.close("all")
        buf.seek(0)

        return {
            "status": "ok",
            "data": {
                "image_b64": base64.b64encode(buf.read()).decode(),
                "genes": req.genes,
                "mode": req.mode,
                "n_cells": int(adata.n_obs),
                "roi_name": req.roi_name,
            },
        }

    except HTTPException:
        raise
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"gene_plot 失敗：{e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/region_stats")
async def post_region_stats(req: RegionRequest):
    """
    回報框選區域的細胞組成，讓使用者**先看組成再決定做圖**。

    回傳 `n_cells`、`median_counts`、`median_genes`、`cluster_counts`
    （取 `obs['leiden']`，無則 `obs['celltypist_label']`，皆無則空 dict）。
    """
    config = load_config()
    try:
        import numpy as np
        import scanpy as sc
        from scipy.sparse import issparse

        adata = sc.read_h5ad(_get_roi_h5ad(config, req.roi_name))
        sel = _filter_by_region(adata, req.region, req.polygon, config, req.roi_name)
        n_cells = int(sel.sum())
        if n_cells == 0:
            return {
                "status": "ok",
                "data": {
                    "n_cells": 0, "median_counts": 0.0, "median_genes": 0.0,
                    "cluster_counts": {}, "roi_name": req.roi_name,
                },
            }

        sub = adata[sel]
        x = sub.X
        counts = np.asarray(x.sum(axis=1)).ravel() if issparse(x) else np.asarray(x).sum(axis=1)
        genes = (
            np.asarray((x > 0).sum(axis=1)).ravel() if issparse(x)
            else (np.asarray(x) > 0).sum(axis=1)
        )

        cluster_counts: dict[str, int] = {}
        for col in ("leiden", "celltypist_label"):
            if col in sub.obs.columns:
                cluster_counts = {
                    str(k): int(v) for k, v in sub.obs[col].value_counts().items() if int(v) > 0
                }
                break

        return {
            "status": "ok",
            "data": {
                "n_cells": n_cells,
                "median_counts": float(np.median(counts)),
                "median_genes": float(np.median(genes)),
                "cluster_counts": cluster_counts,
                "roi_name": req.roi_name,
            },
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"region_stats 失敗：{e}")
        raise HTTPException(status_code=500, detail=str(e))
