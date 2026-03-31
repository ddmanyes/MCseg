"""
Stage 4: 完整分析流水線
整合 QC → 正規化 → HVG → PCA → UMAP → Leiden → 標記基因
提供三步驟分段執行：run_qc_step / run_umap_step / run_heatmap_step
"""
import base64
import logging
from pathlib import Path
from typing import Any

import numpy as np
import scanpy as sc
import anndata as ad

from backend.src.analysis.preprocessing import Preprocessor
from backend.src.utils.constants import VISIUM_UM_PX
from backend.src.analysis.clustering import Analyzer
from backend.src.utils.config import resolve_path

logger = logging.getLogger("pipeline.analysis")


# ─────────────────── TME 基因評分常數（CRC 預設，作為 fallback）──────────────
# 正常使用時由 tissue profile (config/profiles/*.yaml) 的 tme_panels 取代。
# 僅在 config=None 或 profile 未定義 tme_panels 時使用。

TME_PANELS: dict[str, list[str]] = {
    "T_exhausted":  ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "CTLA4"],
    "T_effector":   ["GZMB", "PRF1", "IFNG", "TNF", "NKG7", "GNLY"],
    "Treg":         ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF9"],
    "Macro_M2_TAM": ["CD163", "MRC1", "FOLR2", "TREM2", "C1QA", "C1QB"],
    "Macro_M1":     ["CD80", "CD86", "IL1B", "CXCL10", "CXCL9", "TNF"],
    "cDC_mature":   ["CCR7", "LAMP3", "FSCN1", "MARCKSL1"],
    "Plasma":       ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"],
    "NK_cytotox":   ["GNLY", "NKG7", "KLRB1", "KLRD1", "NCR1"],
}

IMMUNE_PANEL_MAP: dict[str, list[str]] = {
    "cd8":        ["T_exhausted", "T_effector"],
    "cd4":        ["Treg", "T_effector"],
    "regulatory": ["Treg"],
    "t cell":     ["T_exhausted", "T_effector", "Treg"],
    "macrophage": ["Macro_M2_TAM", "Macro_M1"],
    "monocyte":   ["Macro_M2_TAM", "Macro_M1"],
    "cdc":        ["cDC_mature"],
    "dendritic":  ["cDC_mature"],
    "plasma":     ["Plasma"],
    "nk":         ["NK_cytotox"],
}

_PANEL_STATE_LABELS: dict[str, str] = {
    "T_exhausted":  "exhausted",
    "T_effector":   "effector",
    "Treg":         "Treg",
    "Macro_M2_TAM": "M2/TAM",
    "Macro_M1":     "M1",
    "cDC_mature":   "mature DC",
    "Plasma":       "plasma",
    "NK_cytotox":   "cytotoxic",
}


def _build_tme_config(
    config: "dict[str, Any] | None",
) -> "tuple[dict[str, list[str]], dict[str, list[str]], dict[str, str]]":
    """
    從 tissue profile config 動態建立 TME 設定。

    若 config 包含 tme_panels（由 config/profiles/*.yaml 載入），使用 profile 定義。
    否則 fallback 到模組常數（CRC 預設值）。

    Returns
    -------
    (tme_panels, immune_panel_map, panel_state_labels)
    """
    tme_cfg = (config or {}).get("tme_panels") if config else None
    if tme_cfg:
        panels = {name: data.get("genes", []) for name, data in tme_cfg.items()}
        state_labels = {name: data.get("state_label", name) for name, data in tme_cfg.items()}
        ip_map = (config or {}).get("immune_panel_map", IMMUNE_PANEL_MAP)
        profile_name = ((config or {}).get("global") or {}).get("tissue_profile", "unknown")
        logger.info(f"  [TME] 使用 profile tme_panels（{profile_name}）：{list(panels.keys())}")
        return panels, ip_map, state_labels
    # Fallback
    logger.info("  [TME] 使用預設 CRC TME_PANELS（未設定 tissue profile）")
    return TME_PANELS, IMMUNE_PANEL_MAP, _PANEL_STATE_LABELS


def _panel_to_state(panel_name: str) -> str:
    return _PANEL_STATE_LABELS.get(panel_name, panel_name)


def _relevant_panels_for_label(immune_label: str) -> list[str]:
    """根據 CellTypist 大類標籤，回傳應跑的 panel 名稱列表（使用模組預設 map）。"""
    label_lower = immune_label.lower()
    panels: list[str] = []
    for keyword, panel_list in IMMUNE_PANEL_MAP.items():
        if keyword in label_lower:
            panels.extend(panel_list)
    seen: set[str] = set()
    result = []
    for p in panels:
        if p not in seen:
            seen.add(p)
            result.append(p)
    return result


def refine_immune_labels(
    adata_ct,
    leiden_key: str,
    cluster_info: dict[str, dict],
    score_threshold: float = 0.3,
    uncertain_threshold: float = 0.7,
    config: "dict[str, Any] | None" = None,
) -> dict[str, dict]:
    """
    對 source="immune" 的 cluster 執行兩件事：

    A) conf >= uncertain_threshold：
       - 跑 sc.tl.score_genes() 匹配功能狀態 panel
       - 若最高分 > score_threshold，在標籤後附加 [狀態]

    C) conf < uncertain_threshold：
       - 標記 uncertain=True（前端顯示警告色）

    cluster_info 在原地修改後回傳。

    Parameters
    ----------
    config : dict, optional
        pipeline config（含 tissue profile tme_panels）。
        若為 None 則 fallback 到模組層級的 CRC 預設常數。
    """
    # 從 profile 或 fallback 取得動態 TME 設定
    tme_panels, immune_panel_map, panel_state_labels = _build_tme_config(config)

    available_genes = set(adata_ct.var_names)
    leiden_vals = adata_ct.obs[leiden_key].astype(str).values

    # 預計算所有有效 panel 的基因評分
    panel_scores: dict[str, np.ndarray] = {}
    for panel_name, genes in tme_panels.items():
        valid = [g for g in genes if g in available_genes]
        if len(valid) < 2:
            logger.info(f"  [gene score] {panel_name}: 可用基因 {len(valid)}/{len(genes)}，跳過")
            continue
        sc.tl.score_genes(adata_ct, gene_list=valid, score_name=f"_tme_{panel_name}", use_raw=False)
        panel_scores[panel_name] = adata_ct.obs[f"_tme_{panel_name}"].values
        logger.info(f"  [gene score] {panel_name}: {len(valid)} genes 評分完成")

    if not panel_scores:
        logger.warning("  [gene score] 無任何 panel 可用（基因數量不足），跳過精細化")
        for info in cluster_info.values():
            if info.get("source") == "immune":
                info["uncertain"] = info.get("confidence", 0) < uncertain_threshold
                info["state"] = None
                info["state_score"] = 0.0
        return cluster_info

    for cluster, info in cluster_info.items():
        if info.get("source") != "immune":
            info["uncertain"] = False
            info["state"] = None
            info["state_score"] = 0.0
            continue

        conf = info.get("confidence", 0.0)

        # C：低信心 → 標記 uncertain，不做基因評分
        if conf < uncertain_threshold:
            info["uncertain"] = True
            info["state"] = None
            info["state_score"] = 0.0
            continue

        info["uncertain"] = False
        immune_label = info.get("label", "")
        mask = leiden_vals == cluster

        # 決定要跑哪些 panel（使用 profile immune_panel_map）
        label_lower = immune_label.lower()
        target_panels: list[str] = []
        for keyword, panel_list in immune_panel_map.items():
            if keyword in label_lower:
                target_panels.extend(panel_list)
        # 去重保序
        seen_panels: set[str] = set()
        target_panels = [p for p in target_panels if not (p in seen_panels or seen_panels.add(p))]
        if not target_panels:
            target_panels = list(panel_scores.keys())  # 無對應 → 跑全部

        # A：選出得分最高的 panel
        best_panel: str | None = None
        best_score = score_threshold  # 低於此值視為無顯著狀態
        for pname in target_panels:
            if pname not in panel_scores:
                continue
            avg = float(np.mean(panel_scores[pname][mask]))
            if avg > best_score:
                best_score = avg
                best_panel = pname

        if best_panel:
            state_label = panel_state_labels.get(best_panel, best_panel)
            info["label"] = f"{immune_label} [{state_label}]"
            info["state"] = best_panel
            info["state_score"] = round(best_score, 4)
            logger.info(
                f"  Cluster {cluster}: {immune_label} → [{state_label}]  score={best_score:.3f}"
            )
        else:
            info["state"] = None
            info["state_score"] = 0.0

    return cluster_info


# ─────────────────────── multi-ROI merge ──────────────────────────

def merge_all_rois(config: dict[str, Any]) -> Path:
    """
    合併所有 ROI 的 cellpose_cells.h5ad，並將 local µm 座標加回全局偏移。

    同一 H&E + Visium 來源不需要 batch correction。
    輸出：{output_dir}/merged_all_rois.h5ad

    Returns
    -------
    Path  輸出檔路徑
    """
    paths = config["paths"]
    rois = config.get("rois", [])
    out_base = resolve_path(paths["output_dir"]) / "roi"
    output_dir = resolve_path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    adatas = []
    for roi in rois:
        roi_name = roi.get("name", "")
        h5ad = out_base / roi_name / "cellpose_cells.h5ad"
        if not h5ad.exists():
            logger.warning(f"找不到 {roi_name}/cellpose_cells.h5ad，跳過")
            continue

        adata = sc.read_h5ad(str(h5ad))
        adata.obs["roi"] = roi_name

        # 加前綴避免多 ROI 合併後 obs_names 重複（防止後續 Xenium 匯出 Reindex 錯誤）
        adata.obs_names = [f"{roi_name}__{name}" for name in adata.obs_names]

        # 加回全局座標偏移（local µm → global µm）
        if "spatial" in adata.obsm and roi.get("x") is not None:
            pixel_size_um = roi.get("pixel_size_um", VISIUM_UM_PX)
            x_offset_um = float(roi["x"]) * pixel_size_um
            y_offset_um = float(roi["y"]) * pixel_size_um
            coords = adata.obsm["spatial"].copy()
            coords[:, 0] += x_offset_um
            coords[:, 1] += y_offset_um
            adata.obsm["spatial"] = coords
            logger.info(
                f"  {roi_name}: {adata.n_obs:,} 細胞，"
                f"偏移 x+{x_offset_um:.1f} µm, y+{y_offset_um:.1f} µm"
            )
        else:
            logger.info(f"  {roi_name}: {adata.n_obs:,} 細胞（無座標偏移資訊）")

        adatas.append(adata)

    if not adatas:
        raise ValueError("未找到任何 ROI 的 cellpose_cells.h5ad，請先完成 Stage 2 RNA 計數")

    if len(adatas) == 1:
        merged = adatas[0]
        logger.info("只有 1 個 ROI，直接使用")
    else:
        merged = ad.concat(adatas, join="outer", fill_value=0)
        logger.info(
            f"ROI 合併完成：{merged.n_obs:,} 細胞，{merged.n_vars:,} 基因，"
            f"{len(adatas)} 個 ROI"
        )

    out_path = output_dir / "merged_all_rois.h5ad"
    merged.write_h5ad(str(out_path))
    logger.info(f"已儲存合併結果：{out_path}")
    return out_path


# ─────────────────────────── helper ──────────────────────────────

def _encode_image(path: Path) -> str:
    """將圖片檔案 base64 編碼為字串。"""
    return base64.b64encode(path.read_bytes()).decode()


def _fix_log1p(adata: ad.AnnData) -> None:
    """修正 uns['log1p']['base']=None 導致 h5py 序列化錯誤。"""
    if "log1p" in adata.uns and isinstance(adata.uns["log1p"], dict):
        if adata.uns["log1p"].get("base") is None:
            adata.uns["log1p"].pop("base", None)


# ─────────────────────── helpers ──────────────────────────────────

def _generate_overlay_images(
    pre_spatial: "np.ndarray",
    pre_obs_names: list,
    post_obs_names_set: set,
    he_path: Path,
    fig_dir: Path,
    pixel_size_um: float,
    input_source: str = "cellpose",
) -> dict[str, str]:
    """在 H&E 底圖上疊加細胞重心或遮罩輪廓，產生 QC 前後比較圖。

    Returns
    -------
    dict  {"pre_qc": base64_preview_png, "post_qc": base64_preview_png}
    HD 版本（300 DPI）同時存檔至 fig_dir，但不放入回傳值（供下載）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile
    import numpy as np
    from skimage.segmentation import find_boundaries

    figures: dict[str, str] = {}

    if not he_path.exists():
        logger.warning(f"找不到 H&E 影像：{he_path}，跳過疊圖生成")
        return figures

    he = tifffile.imread(str(he_path))
    H, W = he.shape[:2]

    x_px = pre_spatial[:, 0] / pixel_size_um
    y_px = pre_spatial[:, 1] / pixel_size_um

    kept_mask    = np.array([n in post_obs_names_set for n in pre_obs_names])
    removed_mask = ~kept_mask
    n_total  = len(pre_obs_names)
    n_kept   = int(kept_mask.sum())
    n_removed = int(removed_mask.sum())

    mask_path = he_path.parent / "segmentation_masks.npy"
    has_mask = mask_path.exists()
    
    he_pre = he.copy()
    he_post = he.copy()
    
    if he_pre.ndim == 2:
        he_pre = np.stack([he_pre]*3, axis=-1)
        he_post = np.stack([he_post]*3, axis=-1)
    elif he_pre.ndim == 3 and he_pre.shape[-1] == 4:
        he_pre = he_pre[..., :3]
        he_post = he_post[..., :3]

    if has_mask:
        import cv2
        import json
        logger.info(f"為 {he_path.parent.name} 繪製細胞輪廓 ({input_source})")
        seg_mask = np.load(str(mask_path))
        def _get_id(n):
            import re
            m = re.search(r'\d+$', n.split("__")[-1])
            return int(m.group()) if m else -1
            
        proseg_polys = {}
        if input_source == "proseg":
            json_path = he_path.parent / "_proseg_work" / "proseg_results.json"
            if json_path.exists():
                import gzip
                with open(json_path, "rb") as f:
                    magic = f.read(2)
                if magic == b"\x1f\x8b":
                    with gzip.open(json_path, "rt", encoding="utf-8") as f:
                        js = json.load(f)
                else:
                    with open(json_path, "r", encoding="utf-8") as f:
                        js = json.load(f)
                for feat in js.get("features", []):
                    cid = feat["properties"]["cell"]
                    g_type = feat["geometry"].get("type", "Polygon")
                    coords = feat["geometry"].get("coordinates", [])
                    rings = []
                    if g_type == "Polygon" and coords:
                        pts = (np.array(coords[0]) / pixel_size_um).astype(np.int32)
                        rings.append(pts)
                    elif g_type == "MultiPolygon" and coords:
                        for poly in coords:
                            if poly:
                                pts = (np.array(poly[0]) / pixel_size_um).astype(np.int32)
                                rings.append(pts)
                    if rings:
                        proseg_polys[cid] = rings
        
        fallback_all_ids, fallback_kept_ids, fallback_removed_ids = [], [], []
        poly_all, poly_kept, poly_removed = [], [], []

        for name, k, r in zip(pre_obs_names, kept_mask, removed_mask):
            cid = _get_id(name)
            if input_source == "proseg" and cid in proseg_polys:
                poly_all.extend(proseg_polys[cid])
                if k: poly_kept.extend(proseg_polys[cid])
                if r: poly_removed.extend(proseg_polys[cid])
            else:
                fallback_all_ids.append(cid)
                if k: fallback_kept_ids.append(cid)
                if r: fallback_removed_ids.append(cid)

        if poly_all:
            cv2.polylines(he_pre, poly_all, True, (0, 255, 255), 1)
        if poly_kept:
            cv2.polylines(he_post, poly_kept, True, (0, 255, 0), 1)
        if poly_removed:
            cv2.polylines(he_post, poly_removed, True, (255, 0, 0), 1)

        if fallback_all_ids:
            b_all = find_boundaries(np.isin(seg_mask, fallback_all_ids), mode="thick")
            he_pre[b_all] = [255, 255, 0] if input_source == "proseg" else [0, 255, 255]
        if fallback_kept_ids:
            b_kept = find_boundaries(np.isin(seg_mask, fallback_kept_ids), mode="thick")
            he_post[b_kept] = [0, 128, 128] if input_source == "proseg" else [0, 255, 0]
        if fallback_removed_ids:
            b_removed = find_boundaries(np.isin(seg_mask, fallback_removed_ids), mode="thick")
            he_post[b_removed] = [255, 128, 0] if input_source == "proseg" else [255, 0, 0]

    for dpi, suffix, s in [(150, "", 12), (300, "_hd", 3)]:
        figsize = (W / dpi, H / dpi)

        # ── Pre-QC overlay ──
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(he_pre)
        if not has_mask:
            ax.scatter(x_px, y_px, s=s, c="cyan", alpha=0.6, linewidths=0, label=f"All cells ({n_total:,})")
        else:
            if input_source == "proseg":
                if poly_all: ax.plot([], [], "o", color="cyan", markersize=4, label=f"Proseg All ({len(poly_all):,})")
                if fallback_all_ids: ax.plot([], [], "o", color="gold", markersize=4, label=f"Fallback All ({len(fallback_all_ids):,})")
            else:
                ax.plot([], [], "o", color="cyan", markersize=4, label=f"All cells ({n_total:,})")
            
        ax.set_title(f"Pre-QC  ·  {n_total:,} cells", fontsize=9)
        ax.axis("off")
        ax.legend(loc="upper right", fontsize=7, markerscale=2 if not has_mask else 1, framealpha=0.5)
        path = fig_dir / f"overlay_pre_qc{suffix}.png"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        if not suffix:
            figures["pre_qc"] = _encode_image(path)
        logger.info(f"已儲存 {path.name}")

        # ── Post-QC overlay ──
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(he_post)
        if not has_mask:
            if removed_mask.any():
                ax.scatter(x_px[removed_mask], y_px[removed_mask], s=s, c="red", alpha=0.5, linewidths=0, label=f"Removed ({n_removed:,})")
            ax.scatter(x_px[kept_mask], y_px[kept_mask], s=s, c="lime", alpha=0.7, linewidths=0, label=f"Kept ({n_kept:,})")
        else:
            if input_source == "proseg":
                if poly_kept: ax.plot([], [], "o", color="lime", markersize=4, label=f"Proseg Kept ({len(poly_kept):,})")
                if poly_removed: ax.plot([], [], "o", color="red", markersize=4, label=f"Proseg Removed ({len(poly_removed):,})")
                if fallback_kept_ids: ax.plot([], [], "o", color="teal", markersize=4, label=f"Fallback Kept ({len(fallback_kept_ids):,})")
                if fallback_removed_ids: ax.plot([], [], "o", color="darkorange", markersize=4, label=f"Fallback Removed ({len(fallback_removed_ids):,})")
            else:
                if removed_mask.any():
                    ax.plot([], [], "o", color="red", markersize=4, label=f"Removed ({n_removed:,})")
                ax.plot([], [], "o", color="lime", markersize=4, label=f"Kept ({n_kept:,})")
            
        ax.set_title(f"Post-QC  ·  kept {n_kept:,} / removed {n_removed:,}", fontsize=9)
        ax.axis("off")
        ax.legend(loc="upper right", fontsize=7, markerscale=2 if not has_mask else 1, framealpha=0.5)
        path = fig_dir / f"overlay_post_qc{suffix}.png"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        if not suffix:
            figures["post_qc"] = _encode_image(path)
        logger.info(f"已儲存 {path.name}")

    return figures


def _save_roi_overlay_images(
    pre_spatial: "np.ndarray",
    pre_obs_names: list,
    post_obs_names_set: set,
    he_path: Path,
    fig_dir: Path,
    pixel_size_um: float,
    roi_name: str,
    input_source: str = "cellpose",
) -> None:
    """為單一 ROI 生成 QC 前後細胞重心或遮罩輪廓疊圖（供多 ROI 比較視圖），150 DPI。

    儲存至 fig_dir/overlay_{roi_name}_pre_qc.png 及 overlay_{roi_name}_post_qc.png。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile
    import numpy as np
    from skimage.segmentation import find_boundaries

    if not he_path.exists():
        logger.warning(f"找不到 H&E：{he_path}，跳過 {roi_name} per-ROI 疊圖")
        return

    he = tifffile.imread(str(he_path))
    H, W = he.shape[:2]

    x_px = pre_spatial[:, 0] / pixel_size_um
    y_px = pre_spatial[:, 1] / pixel_size_um

    kept_mask    = np.array([n in post_obs_names_set for n in pre_obs_names])
    removed_mask = ~kept_mask
    n_total   = len(pre_obs_names)
    n_kept    = int(kept_mask.sum())
    n_removed = int(removed_mask.sum())

    mask_path = he_path.parent / "segmentation_masks.npy"
    has_mask = mask_path.exists()
    
    he_pre = he.copy()
    he_post = he.copy()
    
    if he_pre.ndim == 2:
        he_pre = np.stack([he_pre]*3, axis=-1)
        he_post = np.stack([he_post]*3, axis=-1)
    elif he_pre.ndim == 3 and he_pre.shape[-1] == 4:
        he_pre = he_pre[..., :3]
        he_post = he_post[..., :3]

    if has_mask:
        import cv2
        import json
        logger.info(f"為 {he_path.parent.name} 繪製細胞輪廓 ({input_source})")
        seg_mask = np.load(str(mask_path))
        def _get_id(n):
            import re
            m = re.search(r'\d+$', n.split("__")[-1])
            return int(m.group()) if m else -1
            
        proseg_polys = {}
        if input_source == "proseg":
            json_path = he_path.parent / "_proseg_work" / "proseg_results.json"
            if json_path.exists():
                import gzip
                with open(json_path, "rb") as f:
                    magic = f.read(2)
                if magic == b"\x1f\x8b":
                    with gzip.open(json_path, "rt", encoding="utf-8") as f:
                        js = json.load(f)
                else:
                    with open(json_path, "r", encoding="utf-8") as f:
                        js = json.load(f)
                for feat in js.get("features", []):
                    cid = feat["properties"]["cell"]
                    g_type = feat["geometry"].get("type", "Polygon")
                    coords = feat["geometry"].get("coordinates", [])
                    rings = []
                    if g_type == "Polygon" and coords:
                        pts = (np.array(coords[0]) / pixel_size_um).astype(np.int32)
                        rings.append(pts)
                    elif g_type == "MultiPolygon" and coords:
                        for poly in coords:
                            if poly:
                                pts = (np.array(poly[0]) / pixel_size_um).astype(np.int32)
                                rings.append(pts)
                    if rings:
                        proseg_polys[cid] = rings
        
        fallback_all_ids, fallback_kept_ids, fallback_removed_ids = [], [], []
        poly_all, poly_kept, poly_removed = [], [], []

        for name, k, r in zip(pre_obs_names, kept_mask, removed_mask):
            cid = _get_id(name)
            if input_source == "proseg" and cid in proseg_polys:
                poly_all.extend(proseg_polys[cid])
                if k: poly_kept.extend(proseg_polys[cid])
                if r: poly_removed.extend(proseg_polys[cid])
            else:
                fallback_all_ids.append(cid)
                if k: fallback_kept_ids.append(cid)
                if r: fallback_removed_ids.append(cid)

        if poly_all:
            cv2.polylines(he_pre, poly_all, True, (0, 255, 255), 1)
        if poly_kept:
            cv2.polylines(he_post, poly_kept, True, (0, 255, 0), 1)
        if poly_removed:
            cv2.polylines(he_post, poly_removed, True, (255, 0, 0), 1)

        if fallback_all_ids:
            b_all = find_boundaries(np.isin(seg_mask, fallback_all_ids), mode="thick")
            he_pre[b_all] = [255, 255, 0] if input_source == "proseg" else [0, 255, 255]
        if fallback_kept_ids:
            b_kept = find_boundaries(np.isin(seg_mask, fallback_kept_ids), mode="thick")
            he_post[b_kept] = [0, 128, 128] if input_source == "proseg" else [0, 255, 0]
        if fallback_removed_ids:
            b_removed = find_boundaries(np.isin(seg_mask, fallback_removed_ids), mode="thick")
            he_post[b_removed] = [255, 128, 0] if input_source == "proseg" else [255, 0, 0]

    dpi = 150
    figsize = (W / dpi, H / dpi)
    s = 12

    # Pre-QC
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(he_pre)
    if not has_mask:
        ax.scatter(x_px, y_px, s=s, c="cyan", alpha=0.6, linewidths=0, label=f"All cells ({n_total:,})")
    else:
        if input_source == "proseg":
            if poly_all: ax.plot([], [], "o", color="cyan", markersize=4, label=f"Proseg All ({len(poly_all):,})")
            if fallback_all_ids: ax.plot([], [], "o", color="gold", markersize=4, label=f"Fallback All ({len(fallback_all_ids):,})")
        else:
            ax.plot([], [], "o", color="cyan", markersize=4, label=f"All cells ({n_total:,})")
    ax.set_title(f"Pre-QC  ·  {n_total:,} cells", fontsize=9)
    ax.axis("off")
    ax.legend(loc="upper right", fontsize=7, markerscale=2 if not has_mask else 1, framealpha=0.5)
    pre_path = fig_dir / f"overlay_{roi_name}_pre_qc.png"
    fig.savefig(str(pre_path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    logger.info(f"已儲存 {pre_path.name}")

    # Post-QC
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(he_post)
    if not has_mask:
        if removed_mask.any():
            ax.scatter(x_px[removed_mask], y_px[removed_mask], s=s, c="red", alpha=0.5, linewidths=0, label=f"Removed ({n_removed:,})")
        ax.scatter(x_px[kept_mask], y_px[kept_mask], s=s, c="lime", alpha=0.7, linewidths=0, label=f"Kept ({n_kept:,})")
    else:
        if input_source == "proseg":
            if poly_kept: ax.plot([], [], "o", color="lime", markersize=4, label=f"Proseg Kept ({len(poly_kept):,})")
            if poly_removed: ax.plot([], [], "o", color="red", markersize=4, label=f"Proseg Removed ({len(poly_removed):,})")
            if fallback_kept_ids: ax.plot([], [], "o", color="teal", markersize=4, label=f"Fallback Kept ({len(fallback_kept_ids):,})")
            if fallback_removed_ids: ax.plot([], [], "o", color="darkorange", markersize=4, label=f"Fallback Removed ({len(fallback_removed_ids):,})")
        else:
            if removed_mask.any():
                ax.plot([], [], "o", color="red", markersize=4, label=f"Removed ({n_removed:,})")
            ax.plot([], [], "o", color="lime", markersize=4, label=f"Kept ({n_kept:,})")
    ax.set_title(f"Post-QC  ·  kept {n_kept:,} / removed {n_removed:,}", fontsize=9)
    ax.axis("off")
    ax.legend(loc="upper right", fontsize=7, markerscale=2 if not has_mask else 1, framealpha=0.5)
    post_path = fig_dir / f"overlay_{roi_name}_post_qc.png"
    fig.savefig(str(post_path), dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    logger.info(f"已儲存 {post_path.name}")

def _generate_roi_comparison_grid(
    roi_names: list[str],
    fig_dir: Path,
) -> str | None:
    """將所有 ROI 的 pre/post QC 疊圖並排成一張比較圖。

    讀取已存在的 overlay_{roi_name}_pre_qc.png 及 overlay_{roi_name}_post_qc.png，
    排列為 (N_rois × 2) grid：左欄=Pre-QC，右欄=Post-QC。

    Returns
    -------
    str | None  base64 PNG 字串；若無任何圖片可讀取則回傳 None。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    valid_roi_names = []
    for rn in roi_names:
        pre  = fig_dir / f"overlay_{rn}_pre_qc.png"
        post = fig_dir / f"overlay_{rn}_post_qc.png"
        if pre.exists() and post.exists():
            valid_roi_names.append(rn)
        else:
            logger.warning(f"_generate_roi_comparison_grid：找不到 {rn} 的疊圖，跳過")

    if not valid_roi_names:
        logger.warning("_generate_roi_comparison_grid：無任何 ROI 疊圖可合成比較圖")
        return None

    n_rows = len(valid_roi_names)
    fig, axes = plt.subplots(
        n_rows, 2,
        figsize=(14, 6 * n_rows),
        squeeze=False,
    )

    for row, rn in enumerate(valid_roi_names):
        pre_img  = mpimg.imread(str(fig_dir / f"overlay_{rn}_pre_qc.png"))
        post_img = mpimg.imread(str(fig_dir / f"overlay_{rn}_post_qc.png"))

        ax_pre  = axes[row][0]
        ax_post = axes[row][1]

        ax_pre.imshow(pre_img)
        ax_pre.set_title(f"{rn}  ·  Pre-QC", fontsize=10, fontweight="bold")
        ax_pre.axis("off")

        ax_post.imshow(post_img)
        ax_post.set_title(f"{rn}  ·  Post-QC", fontsize=10, fontweight="bold")
        ax_post.axis("off")

    plt.suptitle("QC 前後細胞輪廓比較（所有 ROI）", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()

    out_path = fig_dir / "roi_comparison_grid.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"已儲存多 ROI 輪廓比較圖：{out_path.name}")
    return _encode_image(out_path)


# ─────────────────────── Step 1: QC + PCA ─────────────────────────

def _get_analysis_h5ad_dir(config: dict, output_dir: Path) -> Path:
    """
    決定 qc_preprocessed.h5ad / umap_computed.h5ad 的儲存目錄（純路徑計算，不建立目錄）。
    - 單 ROI 模式：output_dir/roi/{roi_name}/
    - 合併模式：output_dir/
    """
    analysis_cfg = config.get("analysis", {})
    rois = config.get("rois", [{"name": "test"}])
    merge_mode = analysis_cfg.get("merge_rois", False) and len(rois) > 1
    if merge_mode:
        return output_dir
    roi_name = rois[0].get("name", "test")
    return output_dir / "roi" / roi_name

def run_qc_step(config: dict[str, Any]) -> dict[str, str]:
    """
    Step 1：QC 前處理 + PCA。

    流程：載入 cellpose_cells.h5ad → 計算 QC 指標 → 繪 violin / scatter →
    過濾細胞基因 → normalize → HVG → PCA → 繪 elbow → 儲存 qc_preprocessed.h5ad

    Returns
    -------
    dict  {chart_name: base64_png}  violin / scatter / elbow
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = config["paths"]
    analysis_cfg = config.get("analysis", {})
    rois = config.get("rois", [{"name": "test"}])

    out_base = resolve_path(paths["output_dir"]) / "roi"
    output_dir = resolve_path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = resolve_path(paths["figure_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 決定輸入 h5ad 來源（cellpose / proseg）
    input_source = analysis_cfg.get("input_source", "cellpose")
    input_h5ad_name = f"{input_source}_cells.h5ad"

    # 多 ROI 合併模式
    merge_mode = analysis_cfg.get("merge_rois", False) and len(rois) > 1
    if merge_mode:
        merged_path = output_dir / "merged_all_rois.h5ad"
        # 永遠重新合併，確保 obs_names 採用最新的 {roi_name}__ 前綴格式
        logger.info("合併模式：重新執行 merge_all_rois() 以確保 obs_names 格式正確...")
        merge_all_rois(config)
        input_h5ad = merged_path
        roi_name = None   # 合併模式不使用單一 ROI H&E 疊圖
        logger.info(f"QC Step（合併模式）：載入 {input_h5ad}")
    else:
        roi_name = rois[0].get("name", "test")
        input_h5ad = out_base / roi_name / input_h5ad_name
        if not input_h5ad.exists():
            if input_source == "proseg":
                fallback = out_base / roi_name / "cellpose_cells.h5ad"
                if fallback.exists():
                    logger.warning(
                        f"找不到 {input_h5ad_name}，自動 fallback 至 cellpose_cells.h5ad"
                    )
                    input_h5ad = fallback
                    input_h5ad_name = "cellpose_cells.h5ad"
                else:
                    raise FileNotFoundError(
                        f"找不到 {input_h5ad_name} 也無 cellpose_cells.h5ad：{input_h5ad}"
                    )
            else:
                raise FileNotFoundError(f"找不到 Cellpose 計數輸出：{input_h5ad}")
        logger.info(f"QC Step：載入 {input_h5ad}（input_source={input_source}）")

    adata = sc.read_h5ad(str(input_h5ad))
    adata.uns["dataset_name"] = input_h5ad_name.replace(".h5ad", "")
    if roi_name is not None:
        adata.uns["active_roi"] = roi_name
    logger.info(f"  {adata.n_obs:,} 細胞, {adata.n_vars:,} 基因")

    preprocessor = Preprocessor(analysis_cfg)
    qc_params = analysis_cfg.get("preprocessing", {}).get("cellular", {})

    # ── 1. QC 指標計算（過濾前，以顯示完整分布） ──
    adata = preprocessor.calculate_qc_metrics(adata, qc_params)

    figures: dict[str, str] = {}

    # ── 2. Violin（過濾前） ──
    qc_keys = [k for k in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if k in adata.obs.columns]
    if qc_keys:
        titles = {
            "total_counts": "Total UMI",
            "n_genes_by_counts": "Genes per Cell",
            "pct_counts_mt": "% Mitochondrial",
        }
        thresholds: dict[str, list[tuple]] = {
            "n_genes_by_counts": [
                (qc_params.get("min_genes", 20), "red", f"min={qc_params.get('min_genes', 20)}"),
                (qc_params.get("max_genes", 8000), "orange", f"max={qc_params.get('max_genes', 8000)}"),
            ],
            "pct_counts_mt": [
                (qc_params.get("max_pct_mito", 20), "red", f"max={qc_params.get('max_pct_mito', 20)}%"),
            ],
        }
        if qc_params.get("min_counts"):
            thresholds["total_counts"] = [
                (qc_params["min_counts"], "red", f"min={qc_params['min_counts']}"),
            ]

        fig, axes = plt.subplots(1, len(qc_keys), figsize=(5 * len(qc_keys), 4))
        if len(qc_keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, qc_keys):
            sc.pl.violin(adata, [key], ax=ax, show=False)
            ax.set_title(titles.get(key, key))
            for val, color, label in thresholds.get(key, []):
                ax.axhline(val, color=color, linestyle="--", linewidth=1.2, alpha=0.8, label=label)
            if thresholds.get(key):
                ax.legend(fontsize=8)
            # 鎖定 y 軸：下限=0，上限夾至 99th percentile × 1.5
            data_99 = float(adata.obs[key].quantile(0.99))
            _, y_top = ax.get_ylim()
            new_top = data_99 * 1.5 if y_top > data_99 * 3 else y_top
            ax.set_ylim(bottom=0, top=new_top)
        plt.tight_layout()
        violin_path = fig_dir / "qc_violin.png"
        fig.savefig(str(violin_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["violin"] = _encode_image(violin_path)
        logger.info("已儲存 qc_violin.png")

    # ── 3. Scatter（counts vs genes，以粒線體著色） ──
    if "total_counts" in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
        color_by = "pct_counts_mt" if "pct_counts_mt" in adata.obs.columns else None
        fig, ax = plt.subplots(figsize=(7, 5))
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color=color_by, ax=ax, show=False)
        ax.set_title("UMI vs Genes" + (" (colored by % Mito)" if color_by else ""))
        plt.tight_layout()
        scatter_path = fig_dir / "qc_scatter.png"
        fig.savefig(str(scatter_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["scatter"] = _encode_image(scatter_path)
        logger.info("已儲存 qc_scatter.png")

    # ── 4. 過濾 + normalize + HVG + PCA ──
    # 過濾前保存空間資訊，供疊圖比較使用
    _pre_obs_names = list(adata.obs_names)
    _pre_spatial   = adata.obsm["spatial"].copy() if "spatial" in adata.obsm else None

    adata = preprocessor.filter_cells(adata, qc_params)
    adata = preprocessor.filter_genes(adata, qc_params)

    # 細胞數過少（< 3）時無法繼續——提前拋出有意義的錯誤
    if adata.n_obs < 3:
        raise ValueError(
            f"QC 過濾後僅剩 {adata.n_obs} 顆細胞，無法繼續分析。\n"
            "建議降低 min_counts / min_genes 門檻，或檢查 Cellpose 計數結果品質。\n"
            f"目前參數：min_genes={qc_params.get('min_genes')}, "
            f"min_counts={qc_params.get('min_counts')}, "
            f"max_pct_mito={qc_params.get('max_pct_mito')}"
        )
    adata = preprocessor.normalize(adata)
    adata = preprocessor.select_hvg(adata)
    adata = preprocessor.run_pca(adata)

    # ── 5. PCA Elbow（直接用 matplotlib，避免舊版 scanpy 不支援 ax 參數） ──
    if "X_pca" in adata.obsm and "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        variance_ratio = adata.uns["pca"]["variance_ratio"]
        n_shown = min(50, len(variance_ratio))
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(range(1, n_shown + 1), variance_ratio[:n_shown], "o-", markersize=4)
        ax.set_xlabel("PC")
        ax.set_ylabel("Variance Ratio")
        ax.set_title("PCA Variance Ratio (Elbow Plot)")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        elbow_path = fig_dir / "pca_elbow.png"
        fig.savefig(str(elbow_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        figures["elbow"] = _encode_image(elbow_path)
        logger.info("已儲存 pca_elbow.png")

    # ── 6. 疊圖（H&E + 細胞重心，QC 前後比較）──
    if _pre_spatial is not None and roi_name is not None:
        # 單一 ROI 模式：生成標準疊圖（回傳給前端）
        he_path = out_base / roi_name / "he_crop.tif"
        pixel_size_um = rois[0].get("pixel_size_um", VISIUM_UM_PX)
        overlay_figs = _generate_overlay_images(
            _pre_spatial, _pre_obs_names, set(adata.obs_names),
            he_path, fig_dir, pixel_size_um, input_source,
        )
        figures.update(overlay_figs)
        # 同時儲存為 per-ROI 命名（供 /roi_overlays 比較視圖）
        _save_roi_overlay_images(
            _pre_spatial, _pre_obs_names, set(adata.obs_names),
            he_path, fig_dir, pixel_size_um, roi_name, input_source,
        )
    elif _pre_spatial is not None and merge_mode:
        # 合併模式：各 ROI 有獨立 H&E 座標系，逐一生成 per-ROI 疊圖
        import numpy as np
        post_set_full = set(adata.obs_names)
        first_roi_done = False
        for roi_cfg in rois:
            rn = roi_cfg.get("name", "")
            prefix = f"{rn}__"
            roi_indices = [i for i, n in enumerate(_pre_obs_names) if n.startswith(prefix)]
            if not roi_indices:
                logger.warning(f"合併模式疊圖：找不到 '{prefix}' 前綴的細胞，跳過 {rn}")
                continue
            idx_arr = np.array(roi_indices)
            # merge_all_rois() 在合併時已加入全域座標偏移（µm），
            # 需減去各 ROI 的全域偏移，還原為 he_crop.tif 的本地座標（µm）
            roi_spatial_global = _pre_spatial[idx_arr]
            px_um = roi_cfg.get("pixel_size_um", VISIUM_UM_PX)
            if roi_cfg.get("x") is not None:
                x_offset_um = float(roi_cfg["x"]) * px_um
                y_offset_um = float(roi_cfg["y"]) * px_um
                roi_spatial_local = roi_spatial_global.copy()
                roi_spatial_local[:, 0] -= x_offset_um
                roi_spatial_local[:, 1] -= y_offset_um
                neg_x = (roi_spatial_local[:, 0] < 0).sum()
                neg_y = (roi_spatial_local[:, 1] < 0).sum()
                if neg_x > 0 or neg_y > 0:
                    logger.warning(f"⚠️ {rn} 座標偏移後出現負值：x={neg_x}, y={neg_y}")
            else:
                roi_spatial_local = roi_spatial_global
            roi_pre_names = [_pre_obs_names[i] for i in roi_indices]
            roi_post_set  = {n for n in post_set_full if n.startswith(prefix)}
            roi_he  = out_base / rn / "he_crop.tif"
            try:
                # 第一個成功的 ROI：同時加入 figures（供 ChartView 顯示）
                if not first_roi_done and roi_he.exists():
                    overlay_figs = _generate_overlay_images(
                        roi_spatial_local, roi_pre_names, roi_post_set,
                        roi_he, fig_dir, px_um, input_source,
                    )
                    figures.update(overlay_figs)
                    first_roi_done = True
                    logger.info(f"合併模式：已使用 {rn} 生成主疊圖（pre_qc / post_qc）")
                # 每個 ROI 都儲存獨立命名的疊圖（供 /roi_overlays 比較視圖）
                _save_roi_overlay_images(
                    roi_spatial_local, roi_pre_names, roi_post_set,
                    roi_he, fig_dir, px_um, rn, input_source,
                )
            except Exception as e:
                logger.warning(f"合併模式疊圖生成失敗（{rn}）：{e}")

        # 所有 ROI 疊圖完成後，生成並排輪廓比較圖
        roi_names_list = [r.get("name", "") for r in rois if r.get("name")]
        comparison_b64 = _generate_roi_comparison_grid(roi_names_list, fig_dir)
        if comparison_b64 is not None:
            figures["roi_comparison"] = comparison_b64
            logger.info("合併模式：已生成多 ROI 輪廓比較圖（roi_comparison）")

    # ── 7. 儲存前處理結果 ──
    _fix_log1p(adata)
    qc_h5ad_dir = _get_analysis_h5ad_dir(config, output_dir)
    qc_h5ad_dir.mkdir(parents=True, exist_ok=True)
    qc_h5ad = qc_h5ad_dir / "qc_preprocessed.h5ad"
    adata.write_h5ad(str(qc_h5ad))
    logger.info(f"QC Step 完成：{adata.n_obs:,} 細胞剩餘，已儲存 {qc_h5ad}")

    return figures


# ──────────────────── Step 2: UMAP 多解析度 ────────────────────────

def run_umap_step(
    config: dict[str, Any],
    resolutions: list[float],
    n_pcs: int = 30,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
) -> dict[str, str]:
    """
    Step 2：讀取前處理結果，計算 UMAP + Leiden（多解析度）。

    Returns
    -------
    dict  {str(resolution): base64_png, "grid": base64_png}
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    fig_dir = resolve_path(paths["figure_dir"])
    h5ad_dir = _get_analysis_h5ad_dir(config, output_dir)

    qc_h5ad = h5ad_dir / "qc_preprocessed.h5ad"
    if not qc_h5ad.exists():
        raise FileNotFoundError(f"找不到前處理資料，請先執行 QC 步驟：{qc_h5ad}")

    logger.info(f"UMAP Step：載入 {qc_h5ad}")
    adata = sc.read_h5ad(str(qc_h5ad))

    # 確保 n_pcs 不超過已計算的 PC 數
    available_pcs = adata.obsm["X_pca"].shape[1] if "X_pca" in adata.obsm else n_pcs
    n_pcs_use = min(n_pcs, available_pcs)

    logger.info(f"建立 KNN 圖 (n_neighbors={n_neighbors}, n_pcs={n_pcs_use})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_use)

    logger.info(f"計算 UMAP (min_dist={min_dist})")
    sc.tl.umap(adata, min_dist=min_dist)

    resolutions = sorted(set(resolutions))
    figures: dict[str, str] = {}

    # ── Leiden 各解析度 + 個別圖 ──
    for res in resolutions:
        key = f"leiden_{res}"
        logger.info(f"  Leiden resolution={res}")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()

        fig_single, ax_single = plt.subplots(figsize=(7, 6))
        sc.pl.umap(
            adata, color=key, ax=ax_single, show=False,
            title=f"Resolution = {res}  ({n_clusters} clusters)",
            frameon=False, legend_loc="on data",
            legend_fontsize=10, legend_fontoutline=2, size=15, alpha=0.8
        )
        plt.tight_layout()
        single_path = fig_dir / f"umap_res{res}.png"
        fig_single.savefig(str(single_path), dpi=150, bbox_inches="tight")
        plt.close(fig_single)
        figures[str(res)] = _encode_image(single_path)

    # ── 合併 Grid 圖 ──
    n = len(resolutions)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig_grid, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows), squeeze=False)
    for i, res in enumerate(resolutions):
        key = f"leiden_{res}"
        row, col = divmod(i, ncols)
        ax = axes[row][col]
        n_clusters = adata.obs[key].nunique()
        sc.pl.umap(
            adata, color=key, ax=ax, show=False,
            title=f"Res={res} ({n_clusters} clusters)",
            frameon=False, legend_loc="on data",
            legend_fontsize=8, legend_fontoutline=1.5, size=10, alpha=0.8
        )
    for i in range(n, nrows * ncols):
        row, col = divmod(i, ncols)
        axes[row][col].set_visible(False)
    plt.tight_layout()
    grid_path = fig_dir / "umap_grid.png"
    fig_grid.savefig(str(grid_path), dpi=150, bbox_inches="tight")
    plt.close(fig_grid)
    figures["grid"] = _encode_image(grid_path)

    # ── 儲存含所有 Leiden 結果的 h5ad ──
    _fix_log1p(adata)
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    umap_h5ad = h5ad_dir / "umap_computed.h5ad"
    adata.write_h5ad(str(umap_h5ad))
    logger.info(f"UMAP Step 完成，已儲存 {umap_h5ad}")

    return figures


# ────────────────────── Step 3: CellTypist 標註 ───────────────────

# 支援的 CellTypist 模型清單（顯示名稱 → 模型檔名）
CELLTYPIST_MODELS: dict[str, str] = {
    "Human CRC（大腸癌）":    "Human_Colorectal_Cancer.pkl",
    "Human Lung Atlas":       "Human_Lung_Atlas.pkl",
    "Immune（精細）":         "Immune_All_Low.pkl",
    "Immune（粗分類）":       "Immune_All_High.pkl",
    "Pan Cancer（泛癌）":     "Pan_Cancer.pkl",
}


def get_cluster_ids(config: dict[str, Any], resolution: float) -> tuple[list[str], dict[str, str]]:
    """取得指定 resolution 的 cluster ID 列表（字串排序）。"""
    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    umap_h5ad = _get_analysis_h5ad_dir(config, output_dir) / "umap_computed.h5ad"
    if not umap_h5ad.exists():
        raise FileNotFoundError(f"找不到 UMAP 資料：{umap_h5ad}")

    adata = sc.read_h5ad(str(umap_h5ad))
    leiden_key = f"leiden_{resolution}"
    if leiden_key not in adata.obs.columns:
        available = [c for c in adata.obs.columns if c.startswith("leiden_")]
        raise ValueError(f"找不到 resolution={resolution}，可用：{available}")

    ids = sorted(adata.obs[leiden_key].astype(str).unique().tolist(), key=lambda x: int(x))

    # 若已有套用的標籤，一起回傳
    existing_labels: dict[str, str] = {}
    if (
        adata.uns.get("cell_type_resolution") == resolution
        and "cell_type_labels" in adata.uns
    ):
        existing_labels = {str(k): str(v) for k, v in adata.uns["cell_type_labels"].items()}

    return ids, existing_labels


def _load_celltypist_model(model_name: str):
    """載入 CellTypist 模型，本地不存在時自動下載。"""
    from celltypist import models as ct_models
    try:
        return ct_models.Model.load(model_name)
    except Exception:
        logger.info(f"  本地找不到模型，嘗試下載 {model_name}...")
        ct_models.download_models(model=model_name)
        return ct_models.Model.load(model_name)


def _tier3_immune_subtype(
    celltypist,
    adata_ct,
    leiden_key: str,
    cluster_info: dict[str, dict],
    tier3_model: str = "Immune_All_Low.pkl",
    tier3_conf_threshold: float = 0.6,
) -> dict[str, dict]:
    """
    Tier 3：對 source=immune 且非 uncertain 的 cluster，
    改用 Immune_All_Low.pkl（98 種亞型）重新標註，
    獲得比 Immune_All_High（~21 種）更精細的亞型。

    例：
      Tier1「T cells」 → Tier3「MAIT cells」/ 「CD8-TEMRA」/ 「Tfh」
      Tier1「Macrophages」 → Tier3「Alveolar macrophages」/ 「Iron-recycling macrophages」
      Tier1「cDC2」 → Tier3「cDC2_CD1C」/ 「cDC2_CLEC10A」

    tier3_conf_threshold：低於此值只紀錄不替換標籤，保留 Tier1 結果。
    """
    immune_clusters = [
        c for c, info in cluster_info.items()
        if info.get("source") == "immune" and not info.get("uncertain", False)
    ]
    if not immune_clusters:
        logger.info("  [Tier3] 無符合條件的免疫 cluster，跳過")
        return cluster_info

    logger.info(f"  [Tier3] 載入 {tier3_model}（高解析度 ~98 亞型）...")
    try:
        model = _load_celltypist_model(tier3_model)
    except Exception as e:
        logger.warning(f"  [Tier3] 模型載入失敗，跳過：{e}")
        return cluster_info

    logger.info(f"  [Tier3] 對 {len(immune_clusters)} 個免疫 cluster 執行精細標註...")
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,
        over_clustering=leiden_key,
    )
    adata_pred = predictions.to_adata()
    leiden_vals = adata_ct.obs[leiden_key].astype(str).values

    for cluster in immune_clusters:
        mask = leiden_vals == cluster
        mv_labels = adata_pred.obs.loc[mask, "majority_voting"]
        conf_col = adata_pred.obs.loc[mask, "conf_score"] if "conf_score" in adata_pred.obs.columns else None

        if len(mv_labels) == 0:
            continue

        tier3_label = str(mv_labels.value_counts().index[0])
        tier3_conf = float(conf_col.mean()) if conf_col is not None else 0.0

        cluster_info[cluster]["tier3_label"] = tier3_label
        cluster_info[cluster]["tier3_conf"] = round(tier3_conf, 4)

        if tier3_conf >= tier3_conf_threshold:
            # 保留 Tier2 附加的功能狀態 suffix
            existing_state = cluster_info[cluster].get("state")
            state_suffix = f" [{_panel_to_state(existing_state)}]" if existing_state else ""
            cluster_info[cluster]["label"] = f"{tier3_label}{state_suffix}"
            logger.info(
                f"  Cluster {cluster}: Tier3 → {tier3_label}{state_suffix}"
                f"  (conf={tier3_conf:.2f} ✓)"
            )
        else:
            logger.info(
                f"  Cluster {cluster}: Tier3 → {tier3_label}"
                f"  (conf={tier3_conf:.2f} < {tier3_conf_threshold}，保留 Tier1 標籤)"
            )

    return cluster_info


def _annotate_one_model(
    celltypist,
    adata_ct,
    model_name: str,
    leiden_key: str,
    clusters: list[str],
) -> dict[str, tuple[str, float]]:
    """
    執行單一模型的 majority voting，回傳每個 cluster 的 (標籤, 平均信心分數)。
    """
    logger.info(f"  載入模型：{model_name}")
    model = _load_celltypist_model(model_name)
    logger.info(f"  執行 majority voting (over_clustering={leiden_key})")
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,
        over_clustering=leiden_key,
    )
    adata_pred = predictions.to_adata()

    results: dict[str, tuple[str, float]] = {}
    leiden_vals = adata_ct.obs[leiden_key].astype(str).values
    for cluster in clusters:
        mask = leiden_vals == cluster
        mv_labels = adata_pred.obs.loc[mask, "majority_voting"]
        conf_scores = adata_pred.obs.loc[mask, "conf_score"] if "conf_score" in adata_pred.obs.columns else None
        if len(mv_labels) == 0:
            results[cluster] = ("Unknown", 0.0)
            continue
        majority = str(mv_labels.value_counts().index[0])
        avg_conf = float(conf_scores.mean()) if conf_scores is not None else 0.0
        results[cluster] = (majority, avg_conf)
    return results


def run_celltypist_annotation(
    config: dict[str, Any],
    resolution: float,
    model_name: str = "Human_Colorectal_Cancer.pkl",
    mode: str = "dual",
    immune_conf_threshold: float = 0.5,
    score_threshold: float = 0.3,
    uncertain_threshold: float = 0.7,
    enable_tier3: bool = False,
    tier3_conf_threshold: float = 0.6,
) -> dict[str, dict]:
    """
    使用 CellTypist 對每個 Leiden cluster 進行多數投票標註。

    mode="dual"   : 先跑 Immune_All_High，再跑 model_name（CRC/組織），按信心選勝者
    mode="single" : 只跑 model_name

    Returns
    -------
    dict  {cluster_id: {"label", "confidence", "source",
                         "immune_label", "immune_conf",   # dual only
                         "crc_label",   "crc_conf"}}      # dual only
    """
    try:
        import celltypist
    except ImportError:
        raise ImportError("請先安裝 celltypist：uv add celltypist")

    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    umap_h5ad = _get_analysis_h5ad_dir(config, output_dir) / "umap_computed.h5ad"
    if not umap_h5ad.exists():
        raise FileNotFoundError(f"找不到 UMAP 資料：{umap_h5ad}")

    logger.info(f"CellTypist 標註：載入 {umap_h5ad}  [mode={mode}]")
    adata = sc.read_h5ad(str(umap_h5ad))

    leiden_key = f"leiden_{resolution}"
    if leiden_key not in adata.obs.columns:
        raise ValueError(f"找不到 resolution={resolution} 的 Leiden 結果")

    # CellTypist 需要 log1p(CPM)——從 raw counts 重建
    if "counts" in adata.layers:
        import anndata as _ad
        adata_ct = _ad.AnnData(
            X=adata.layers["counts"].copy(),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )
        adata_ct.var_names_make_unique()
        adata_ct.obs_names_make_unique()
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)
        logger.info("  從 layers['counts'] 重建 log1p(CPM) 資料")
    else:
        adata_ct = adata.copy()
        adata_ct.var_names_make_unique()
        adata_ct.obs_names_make_unique()
        logger.warning("  找不到 layers['counts']，直接使用 adata.X（可能為已縮放資料）")

    clusters = sorted(adata_ct.obs[leiden_key].astype(str).unique(), key=lambda x: int(x))

    if mode == "dual":
        immune_model = "Immune_All_High.pkl"
        logger.info(f"  [雙模型] 免疫模型：{immune_model}，組織模型：{model_name}")
        logger.info(f"  [雙模型] 免疫信心閾值：{immune_conf_threshold}")

        immune_results = _annotate_one_model(celltypist, adata_ct, immune_model, leiden_key, clusters)
        crc_results = _annotate_one_model(celltypist, adata_ct, model_name, leiden_key, clusters)

        cluster_info: dict[str, dict] = {}
        for cluster in clusters:
            immune_label, immune_conf = immune_results[cluster]
            crc_label, crc_conf = crc_results[cluster]

            # 免疫優先：非 Non-immune 且信心超過閾值
            if immune_label != "Non-immune" and immune_conf >= immune_conf_threshold:
                source, label, confidence = "immune", immune_label, immune_conf
            else:
                source, label, confidence = "crc", crc_label, crc_conf

            cluster_info[cluster] = {
                "label": label,
                "confidence": round(confidence, 4),
                "source": source,
                "immune_label": immune_label,
                "immune_conf": round(immune_conf, 4),
                "crc_label": crc_label,
                "crc_conf": round(crc_conf, 4),
            }
            logger.info(
                f"  Cluster {cluster:>2} → [{source:>6}] {label}"
                f"  (immune={immune_label} {immune_conf:.2f} | crc={crc_label} {crc_conf:.2f})"
            )
    else:
        single_results = _annotate_one_model(celltypist, adata_ct, model_name, leiden_key, clusters)
        cluster_info = {}
        for cluster in clusters:
            label, confidence = single_results[cluster]
            cluster_info[cluster] = {
                "label": label,
                "confidence": round(confidence, 4),
                "source": "single",
            }
            logger.info(f"  Cluster {cluster:>2} → {label}  (conf={confidence:.2f})")

    # Tier 2 (A+C)：基因評分精細化（僅 dual mode）
    if mode == "dual":
        logger.info("  [Tier2] 對免疫 cluster 執行基因評分精細化...")
        cluster_info = refine_immune_labels(
            adata_ct,
            leiden_key,
            cluster_info,
            score_threshold=score_threshold,
            uncertain_threshold=uncertain_threshold,
            config=config,
        )

    # Tier 3：Immune_All_Low.pkl 亞型精細化（選用）
    if mode == "dual" and enable_tier3:
        logger.info("  [Tier3] 啟用精細免疫亞型標註（Immune_All_Low.pkl）...")
        cluster_info = _tier3_immune_subtype(
            celltypist,
            adata_ct,
            leiden_key,
            cluster_info,
            tier3_model="Immune_All_Low.pkl",
            tier3_conf_threshold=tier3_conf_threshold,
        )

    logger.info(f"CellTypist 標註完成：{len(cluster_info)} 個 cluster")
    return cluster_info


def apply_cluster_labels(
    config: dict[str, Any],
    resolution: float,
    labels: dict[str, str],
) -> None:
    """
    將 cluster 標籤套用到 umap_computed.h5ad。
    結果存入 adata.obs["cell_type"]，並記錄於 adata.uns["cell_type_labels"]。
    """
    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    umap_h5ad = _get_analysis_h5ad_dir(config, output_dir) / "umap_computed.h5ad"
    if not umap_h5ad.exists():
        raise FileNotFoundError(f"找不到 UMAP 資料：{umap_h5ad}")

    logger.info(f"套用標籤：resolution={resolution}，{len(labels)} 個 cluster")
    adata = sc.read_h5ad(str(umap_h5ad))

    leiden_key = f"leiden_{resolution}"
    if leiden_key not in adata.obs.columns:
        raise ValueError(f"找不到 resolution={resolution} 的 Leiden 結果")

    adata.obs["cell_type"] = (
        adata.obs[leiden_key].astype(str).map(labels).fillna("Unknown")
    )
    adata.uns["cell_type_resolution"] = resolution
    adata.uns["cell_type_labels"] = labels

    _fix_log1p(adata)
    adata.write_h5ad(str(umap_h5ad))
    logger.info(f"標籤已套用並儲存：{umap_h5ad}")


# ────────────────────── Step 4: Heatmap ──────────────────────────

def run_heatmap_step(
    config: dict[str, Any],
    resolution: float,
    n_top_genes: int = 20,
    n_heatmap_genes: int = 50,
) -> dict[str, str]:
    """
    Step 3：針對指定解析度同時產生兩張 marker gene 圖表。

    - heatmap：seaborn.clustermap，顯示方差最高的 n_heatmap_genes 個 HVGs，
              行（cluster）與列（基因）均有樹枝圖
    - dotplot：sc.pl.dotplot，每 cluster 取 n_top_genes 個 marker 基因，
              點大小 = 表達細胞比例，顏色 = 平均表達量

    Returns
    -------
    dict  {"heatmap": base64_png, "dotplot": base64_png}
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.sparse import issparse

    paths = config["paths"]
    output_dir = resolve_path(paths["output_dir"])
    fig_dir = resolve_path(paths["figure_dir"])

    umap_h5ad = _get_analysis_h5ad_dir(config, output_dir) / "umap_computed.h5ad"
    if not umap_h5ad.exists():
        raise FileNotFoundError(f"找不到 UMAP 資料，請先執行 UMAP 步驟：{umap_h5ad}")

    logger.info(f"Heatmap Step：載入 {umap_h5ad}")
    adata = sc.read_h5ad(str(umap_h5ad))

    leiden_key = f"leiden_{resolution}"
    available = [c for c in adata.obs.columns if c.startswith("leiden_")]
    if leiden_key not in adata.obs.columns:
        raise ValueError(
            f"找不到 resolution={resolution} 的 Leiden 結果。可用：{available}"
        )

    clusters = adata.obs[leiden_key].cat.categories.tolist()

    # 若有套用 cell_type 標籤（且對應同一 resolution），用標籤取代數字
    cell_type_map: dict[str, str] | None = None
    if (
        adata.uns.get("cell_type_resolution") == resolution
        and "cell_type_labels" in adata.uns
    ):
        cell_type_map = {str(k): str(v) for k, v in adata.uns["cell_type_labels"].items()}
        logger.info(f"  使用 cell_type 標籤覆寫 y 軸：{cell_type_map}")

    # ── 1. Heatmap：取方差最高的 n_heatmap_genes 個 HVGs ────────────
    logger.info(f"Heatmap：從 HVGs 中選方差最高的 {n_heatmap_genes} 個基因")

    if "highly_variable" in adata.var.columns:
        hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()
    else:
        hvg_names = adata.var_names.tolist()

    # 計算每個 HVG 在全部細胞的方差，取 top N
    if len(hvg_names) > n_heatmap_genes:
        import scipy.sparse as sp
        X_hvg = adata[:, hvg_names].X
        if sp.issparse(X_hvg):
            # 稀疏矩陣：E[X^2] - E[X]^2
            mean_sq = np.array(X_hvg.power(2).mean(axis=0)).ravel()
            mean_   = np.array(X_hvg.mean(axis=0)).ravel()
            gene_var = mean_sq - mean_ ** 2
        else:
            gene_var = X_hvg.var(axis=0)
        top_idx = np.argsort(gene_var)[::-1][:n_heatmap_genes]
        gene_list = [hvg_names[i] for i in sorted(top_idx)]
        logger.info(f"  HVG 總數 {len(hvg_names)}，選用方差最高的 {len(gene_list)} 個")
    else:
        gene_list = hvg_names
        logger.info(f"  HVG 總數 {len(hvg_names)}（≤ {n_heatmap_genes}），全部使用")

    # 計算各 cluster 的平均表達量（n_clusters × n_genes）
    X_sub = adata[:, gene_list].X
    if issparse(X_sub):
        X_sub = X_sub.toarray()

    cluster_means: dict[str, Any] = {}
    for cluster in clusters:
        mask = (adata.obs[leiden_key] == cluster).values
        cluster_means[str(cluster)] = X_sub[mask].mean(axis=0)

    df_mean = pd.DataFrame(cluster_means, index=gene_list).T  # shape: (n_clusters, n_genes)

    # 若有 cell_type 標籤，重新命名 index，並將同名 cell type 合併（加權平均）
    if cell_type_map:
        df_mean.index = [cell_type_map.get(str(c), f"C{c}") for c in clusters]
        # 計算每個 cluster 的 cell 數（用於加權平均）
        cluster_sizes = {
            str(c): int((adata.obs[leiden_key] == c).sum()) for c in clusters
        }
        size_series = pd.Series(
            [cluster_sizes[str(c)] for c in clusters], index=df_mean.index
        )
        # 若有重複名稱，做加權平均後合併
        if df_mean.index.duplicated().any():
            df_mean["_n"] = size_series.values
            df_mean = (
                df_mean.groupby(df_mean.index)
                .apply(lambda g: pd.Series(
                    np.average(g.drop(columns="_n").values, axis=0, weights=g["_n"].values),
                    index=gene_list,
                ))
            )
            logger.info(f"  同名 cell type 合併後：{df_mean.shape[0]} 列")

    # Z-score per gene (column)：讓顏色反映「該 cluster 相對其他 cluster 的高低」
    # 截斷至 [-2.5, 2.5] 避免個別 outlier 撐爆色域導致其他欄全變同色
    col_mean = df_mean.mean(axis=0)
    col_std  = df_mean.std(axis=0).replace(0, 1)   # std=0 的基因設為 1 防除零
    df_scaled = ((df_mean - col_mean) / col_std).clip(-2.5, 2.5)

    n_genes_total = len(gene_list)
    n_clusters = len(clusters)
    # 寬度：每個基因 0.25 吋，最窄 10 吋，最寬 80 吋
    hm_w = max(10, min(n_genes_total * 0.25 + 3, 80))
    # 高度：每個 cluster 0.6 吋，同時確保最小縱橫比不低於 1:4（避免過扁）
    hm_h = max(4, n_clusters * 0.6 + 2, hm_w / 4)
    show_gene_labels = n_genes_total <= 80

    # 只有 ≥2 個 cluster / gene 才能做層次聚類；否則關閉對應樹枝圖避免報錯
    do_row_cluster = n_clusters >= 2
    do_col_cluster = n_genes_total >= 2

    logger.info(
        f"clustermap：{n_clusters} clusters × {n_genes_total} genes，"
        f"show_labels={show_gene_labels}，row_cluster={do_row_cluster}，col_cluster={do_col_cluster}"
    )

    plt.close("all")
    g = sns.clustermap(
        df_scaled,
        cmap="RdBu_r",
        vmin=-2.5, vmax=2.5,
        figsize=(hm_w, hm_h),
        yticklabels=True,              # cluster 標籤
        xticklabels=show_gene_labels,  # 基因標籤（數量少時才顯示）
        row_cluster=do_row_cluster,    # cluster 樹枝圖（左側 / 行）
        col_cluster=do_col_cluster,    # 基因樹枝圖（上方 / 列）
        linewidths=0,
        cbar_kws={"label": "Z-score (mean expr.)", "shrink": 0.5},
    )
    if not show_gene_labels:
        g.ax_heatmap.set_xlabel(f"Genes ({n_genes_total} HVGs)")
    g.ax_heatmap.set_ylabel("Cluster")
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)

    fig_path_heatmap = fig_dir / "heatmap.png"
    g.savefig(str(fig_path_heatmap), dpi=150, bbox_inches="tight")
    plt.close("all")
    logger.info(f"heatmap.png 已儲存 → {fig_path_heatmap}")

    # ── 2. Dotplot：每 cluster 取 n_top_genes 個 marker 基因 ─────────
    logger.info(f"Dotplot：計算 rank_genes_groups (groupby={leiden_key}, n_genes={n_top_genes})")
    sc.tl.rank_genes_groups(adata, groupby=leiden_key, method="wilcoxon", n_genes=n_top_genes)

    top_genes_per_cluster: dict[str, list[str]] = {}
    seen: set[str] = set()
    for cluster in clusters:
        cluster_genes = list(adata.uns["rank_genes_groups"]["names"][cluster][:n_top_genes])
        unique_genes = [g for g in cluster_genes if g not in seen]
        seen.update(unique_genes)
        if unique_genes:
            top_genes_per_cluster[f"C{cluster}"] = unique_genes

    total_dot_genes = sum(len(v) for v in top_genes_per_cluster.values())
    dp_w = max(12, total_dot_genes * 0.3)
    dp_h = max(5, n_clusters * 0.6 + 3)

    fig_path_dotplot = fig_dir / "dotplot.png"
    try:
        dp = sc.pl.dotplot(
            adata,
            var_names=top_genes_per_cluster,
            groupby=leiden_key,
            show=False,
            return_fig=True,
            standard_scale="var",
            dendrogram=True,
            figsize=(dp_w, dp_h),
        )
        dp.savefig(str(fig_path_dotplot), dpi=150, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        logger.warning(f"dotplot return_fig 失敗（{e}），改用 plt.savefig fallback")
        plt.close("all")
        sc.pl.dotplot(
            adata,
            var_names=top_genes_per_cluster,
            groupby=leiden_key,
            show=False,
            standard_scale="var",
            dendrogram=True,
        )
        plt.gcf().set_size_inches(dp_w, dp_h)
        plt.savefig(str(fig_path_dotplot), dpi=150, bbox_inches="tight")
        plt.close("all")

    logger.info(f"dotplot.png 已儲存 → {fig_path_dotplot}")

    return {
        "heatmap": _encode_image(fig_path_heatmap),
        "dotplot": _encode_image(fig_path_dotplot),
    }


def run_analysis_pipeline(config: dict[str, Any]) -> ad.AnnData:
    """
    執行完整分析流程。

    Parameters
    ----------
    config : pipeline.yaml 配置字典

    Returns
    -------
    AnnData 含有聚類結果
    """
    paths = config["paths"]
    analysis_cfg = config.get("analysis", {})

    # 確定輸入 h5ad（Cellpose 計數輸出）
    out_base = resolve_path(paths["output_dir"]) / "roi"
    roi_name = config.get("rois", [{"name": "text"}])[0].get("name", "text")
    input_h5ad = out_base / roi_name / "cellpose_cells.h5ad"

    if not input_h5ad.exists():
        raise FileNotFoundError(f"找不到 Cellpose 計數輸出：{input_h5ad}")

    logger.info(f"載入資料：{input_h5ad}")
    adata = sc.read_h5ad(str(input_h5ad))
    adata.uns["dataset_name"] = "cellpose_cells"
    logger.info(f"  {adata.n_obs:,} 細胞, {adata.n_vars:,} 基因")

    # 預處理
    preprocessor = Preprocessor(analysis_cfg)
    adata = preprocessor.preprocess(adata, run_pca=True, qc_key="cellular")

    # 聚類
    analyzer = Analyzer(analysis_cfg)
    adata = analyzer.run_clustering(adata)

    # 儲存結果
    output_dir = resolve_path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "clustered_final.h5ad"
    adata.write_h5ad(str(out_path))
    logger.info(f"已儲存：{out_path}")

    # 產生圖表
    fig_dir = resolve_path(paths["figure_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)
    _save_figures(adata, fig_dir, analysis_cfg)

    return adata


def _save_figures(adata: ad.AnnData, fig_dir: Path, analysis_cfg: dict) -> None:
    """產生並儲存標準圖表"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dpi = analysis_cfg.get("visualization", {}).get("figure_dpi", 300)

    # UMAP
    if "X_umap" in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.umap(
            adata, color="cluster", ax=axes[0], show=False, title="Leiden Clusters",
            frameon=False, legend_loc="on data", legend_fontsize=10, legend_fontoutline=2, size=15, alpha=0.8
        )
        if "total_counts" in adata.obs:
            sc.pl.umap(
                adata, color="total_counts", ax=axes[1], show=False, title="Total Counts",
                frameon=False, cmap="magma", size=15, alpha=0.8
            )
        plt.tight_layout()
        fig.savefig(str(fig_dir / "umap.png"), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"已儲存：{fig_dir / 'umap.png'}")

    # Leiden 分布
    if "cluster" in adata.obs:
        fig, ax = plt.subplots(figsize=(8, 5))
        cluster_counts = adata.obs["cluster"].value_counts().sort_index()
        cluster_counts.plot(kind="bar", ax=ax)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Cell Count")
        ax.set_title("Cells per Cluster")
        plt.tight_layout()
        fig.savefig(str(fig_dir / "cluster_distribution.png"), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
