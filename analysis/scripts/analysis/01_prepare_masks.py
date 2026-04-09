"""
00_prepare_masks.py
===================
收集 3 種分割方法的遮罩（7 個 CRC ROI）：
  1. V12          → results/masks/v12_roi{i}.npy   （直接複製）
  2. Space Ranger → results/masks/sr_roi{i}.npy    （GeoJSON 光柵化）
  3. Cellpose v3  → results/masks/cpv3_roi{i}.npy  （最佳參數重跑）

輸出：uint32 mask，shape=(1000,1000)，0 = 未歸屬，>0 = cell_id
"""

from __future__ import annotations

import sys
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm

# ── 允許 import 上游模組 ────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS = cfg["paths"]
DATA  = cfg["data"]
CP3   = cfg["cellpose_v3_best"]

MASKS_DIR = ROOT / PATHS["masks_dir"]
MASKS_DIR.mkdir(parents=True, exist_ok=True)

# ROI 定義
with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

PX_SIZE = DATA["pixel_size_um"]   # 0.2738 µm/px

# ═══════════════════════════════════════════════════════════════════════════
# 1. V12 遮罩（直接複製）
# ═══════════════════════════════════════════════════════════════════════════

def prepare_v12():
    print("\n[V12] 複製預測遮罩...")
    src_dir = Path(PATHS["v12_masks_dir"])
    for roi_name in ROI_INFO.keys():
        src = src_dir / f"{roi_name}_pred_mask.npy"
        dst = MASKS_DIR / f"v12_{roi_name}.npy"
        if dst.exists():
            continue
        if not src.exists():
            print(f"  ⚠️  {roi_name} 找不到預測遮罩 ({src})，跳過")
            continue
        mask = np.load(src).astype(np.uint32)
        np.save(dst, mask)
        n_cells = int(mask.max())
        print(f"  {roi_name}: {n_cells} 顆細胞 → {dst.name}")


# ═══════════════════════════════════════════════════════════════════════════
# 2. Space Ranger 遮罩（GeoJSON 光柵化）
# ═══════════════════════════════════════════════════════════════════════════

def prepare_spaceranger():
    """從 GeoJSON 多邊形光柵化到 1000×1000 ROI 遮罩。"""
    from shapely.geometry import shape
    from skimage.draw import polygon as sk_polygon

    geojson_path = Path(PATHS["cell_seg_geojson"])
    print(f"\n[SR] 讀取 GeoJSON: {geojson_path.name}")
    print("  （首次讀取約需 30 秒，~1M 細胞）")

    import json
    with open(geojson_path, encoding="utf-8") as f:
        geojson = json.load(f)

    features = geojson.get("features", [])
    print(f"  共 {len(features):,} 個 feature")

    for roi_name, roi in tqdm(ROI_INFO.items(), desc="光柵化 ROI"):
        dst = MASKS_DIR / f"sr_{roi_name}.npy"
        if dst.exists():
            print(f"  {roi_name} 已存在，跳過")
            continue

        x0, y0 = roi["x0"], roi["y0"]
        x1, y1 = roi["x1"], roi["y1"]
        H = y1 - y0
        W = x1 - x0
        mask = np.zeros((H, W), dtype=np.uint32)
        cell_id = 1

        for feat in features:
            geom = shape(feat["geometry"])
            # GeoJSON 座標是 (x, y) = (col, row) in fullres pixels
            gx0, gy0, gx1, gy1 = geom.bounds
            # 快速跳過 bbox 不交叉的細胞
            if gx1 < x0 or gx0 >= x1 or gy1 < y0 or gy0 >= y1:
                continue

            # 取多邊形頂點（只處理外環）
            if geom.geom_type == "Polygon":
                coords = np.array(geom.exterior.coords)
            elif geom.geom_type == "MultiPolygon":
                coords = np.array(max(geom.geoms, key=lambda g: g.area).exterior.coords)
            else:
                continue

            # 轉換到 ROI 局部座標（row = y - y0, col = x - x0）
            rows = coords[:, 1] - y0   # y → row
            cols = coords[:, 0] - x0   # x → col

            # 過濾到 ROI 範圍內
            rr, cc = sk_polygon(rows, cols, shape=(H, W))
            if len(rr) == 0:
                continue

            mask[rr, cc] = cell_id
            cell_id += 1

        np.save(dst, mask)
        print(f"  {roi_name}: {cell_id-1} 顆細胞 → {dst.name}")


# ═══════════════════════════════════════════════════════════════════════════
# 3. Cellpose v3 Best Params（重跑）
# ═══════════════════════════════════════════════════════════════════════════

def prepare_cellpose_v3():
    """使用 best_params_crc_v3.json 的最佳參數對每個 ROI 執行 Cellpose。"""
    from cellpose import models
    from scipy.ndimage import label as ndlabel
    from scipy.ndimage import zoom
    import cv2

    # 讀取最佳參數（從 json 確保最新值）
    with open(PATHS["cellpose_best_params"]) as f:
        bp = json.load(f)

    dia_small_px = bp["dia_small_um"] / PX_SIZE
    dia_large_px = bp["dia_large_um"] / PX_SIZE
    flow_thr     = bp["flow_threshold"]
    prob_thr     = bp["cellprob_threshold"]
    dilation_px  = int(bp["dilation_px"])
    frag_thr     = int(bp["fragment_threshold"])
    clahe_clip   = bp["clahe_clip_limit"]

    print(f"\n[CPv3] Cellpose best params:")
    print(f"  dia_small={dia_small_px:.1f}px  dia_large={dia_large_px:.1f}px")
    print(f"  flow={flow_thr:.3f}  prob={prob_thr:.3f}  dil={dilation_px}px")

    model = models.CellposeModel(gpu=True, pretrained_model="cyto3")

    rois_dir = Path(PATHS["crc_he_rois_dir"])

    for roi_name, roi in tqdm(ROI_INFO.items(), desc="Cellpose v3"):
        dst = MASKS_DIR / f"cpv3_{roi_name}.npy"
        if dst.exists():
            print(f"  {roi_name} 已存在，跳過")
            continue

        # 載入 HE 影像（fullres ROI patch）
        he_npy = rois_dir / f"{roi_name}_he.npy"
        he_png = rois_dir / f"{roi_name}_he.png"
        if he_npy.exists():
            img = np.load(he_npy)
        elif he_png.exists():
            img = cv2.cvtColor(cv2.imread(str(he_png)), cv2.COLOR_BGR2RGB)
        else:
            print(f"  ⚠️  {roi_name} 找不到 HE 影像，跳過")
            continue

        if img.dtype != np.uint8:
            img = (img * 255).clip(0, 255).astype(np.uint8)

        # CLAHE 增強
        img_clahe = _apply_clahe(img, clahe_clip)

        # Cellpose 小直徑（核）
        masks_small, _, _ = model.eval(
            img_clahe, diameter=dia_small_px,
            flow_threshold=flow_thr, cellprob_threshold=prob_thr,
            channels=[0, 0]
        )

        # Cellpose 大直徑（細胞體）
        masks_large, _, _ = model.eval(
            img_clahe, diameter=dia_large_px,
            flow_threshold=flow_thr, cellprob_threshold=prob_thr,
            channels=[0, 0]
        )

        # LOGIC_A 合併
        mask = _merge_logic_a(masks_small, masks_large, frag_thr)

        # Voronoi 擴張（expand_labels 等效）
        mask = _voronoi_expand(mask, max_dist=dilation_px)

        np.save(dst, mask.astype(np.uint32))
        print(f"  {roi_name}: {int(mask.max())} 顆細胞 → {dst.name}")


# ── 輔助函式 ─────────────────────────────────────────────────────────────

def _apply_clahe(img: np.ndarray, clip_limit: float = 2.0) -> np.ndarray:
    """對 RGB HE 影像的 L 通道（LAB）應用 CLAHE。"""
    import cv2
    lab = cv2.cvtColor(img, cv2.COLOR_RGB2LAB)
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(8, 8))
    lab[:, :, 0] = clahe.apply(lab[:, :, 0])
    return cv2.cvtColor(lab, cv2.COLOR_LAB2RGB)


def _merge_logic_a(masks_small: np.ndarray, masks_large: np.ndarray,
                   fragment_threshold: int = 50) -> np.ndarray:
    """
    LOGIC_A 合併：large cell 區域內的 small cells ≤1 個顯著 → 用 large cell 取代。
    直接複製自 crc_he_seg/backend/src/segmentation/cellpose_runner.py。
    """
    uids, counts = np.unique(masks_small, return_counts=True)
    small_areas = dict(zip(uids, counts))

    merged = masks_small.copy().astype(np.int32)
    large_ids = np.unique(masks_large)
    large_ids = large_ids[large_ids > 0]
    next_id = int(merged.max()) + 1

    for lid in large_ids:
        region = masks_large == lid
        overlapping = np.unique(masks_small[region])
        overlapping = overlapping[overlapping > 0]
        significant = [s for s in overlapping if small_areas.get(s, 0) >= fragment_threshold]
        if len(significant) <= 1:
            merged[region] = next_id
            next_id += 1

    return merged


def _voronoi_expand(mask: np.ndarray, max_dist: int = 6) -> np.ndarray:
    """
    Voronoi 式擴張：每個像素分配給最近的細胞，距離不超過 max_dist px。
    等效於 skimage.segmentation.expand_labels(mask, distance=max_dist)。
    """
    try:
        from skimage.segmentation import expand_labels
        return expand_labels(mask, distance=max_dist)
    except ImportError:
        from scipy.ndimage import distance_transform_edt, label as ndlabel
        # Fallback：binary dilation 近似
        from scipy.ndimage import binary_dilation, generate_binary_structure
        struct = generate_binary_structure(2, 2)
        expanded = mask.copy()
        for _ in range(max_dist):
            expanded_next = expanded.copy()
            bg = expanded == 0
            if not bg.any():
                break
            dilated = binary_dilation(expanded > 0, struct)
            # 只擴張到背景
            for cid in np.unique(expanded)[1:]:
                cell_dilated = binary_dilation(expanded == cid, struct)
                fill_region = cell_dilated & bg
                expanded_next[fill_region] = cid
            expanded = expanded_next
        return expanded


# ═══════════════════════════════════════════════════════════════════════════
# 4. visiumHD_pipeline_3 預設參數
# ═══════════════════════════════════════════════════════════════════════════

def prepare_p3():
    """
    visiumHD_pipeline_3 預設參數（pipeline.yaml + CRC profile 合併後）：
      - model: nuclei（pipeline.yaml 覆蓋 crc.yaml 的 cyto2）
      - LOGIC_A 策略：dia_small=5µm, dia_large=60µm
      - flow_threshold=0.8, cellprob_threshold=-2.0, fragment_threshold=20
      - clahe_clip_limit=1.0, normalize_stains=false
      - 6px Voronoi dilation（對應 rna_counting.dilation_px=6）
    """
    from cellpose import models
    import cv2

    # ── 從 visiumHD_pipeline_3/config/pipeline.yaml 讀取 ─────────────────
    import yaml as _yaml
    p3_cfg_path = Path("/Volumes/SSD/plan_a/visiumHD_pipeline_3/config/pipeline.yaml")
    p3_cfg = _yaml.safe_load(p3_cfg_path.read_text())
    seg = p3_cfg["segmentation"]
    strat = seg["strategy"]

    model_type    = seg["cellpose_model"]["model_type"]   # nuclei
    dia_small_px  = strat["dia_small"] / PX_SIZE          # 5.0 µm → ~18.3 px
    dia_large_px  = strat["dia_large"] / PX_SIZE          # 60.0 µm → ~219 px
    flow_thr      = strat["flow_threshold"]               # 0.8
    prob_thr      = strat["cellprob_threshold"]           # -2.0
    frag_thr      = strat["fragment_threshold"]           # 20
    clahe_clip    = seg["preprocessing"]["clahe_clip_limit"]  # 1.0
    dilation_px   = p3_cfg["rna_counting"]["dilation_px"] # 6

    print(f"\n[P3] visiumHD_pipeline_3 預設參數:")
    print(f"  model={model_type}  dia_small={dia_small_px:.1f}px  dia_large={dia_large_px:.1f}px")
    print(f"  flow={flow_thr}  prob={prob_thr}  frag={frag_thr}  dilation={dilation_px}px")

    model = models.CellposeModel(gpu=True, pretrained_model=model_type)
    rois_dir = Path(PATHS["crc_he_rois_dir"])

    for roi_name, roi in tqdm(ROI_INFO.items(), desc="P3 Cellpose"):
        dst = MASKS_DIR / f"p3_{roi_name}.npy"
        if dst.exists():
            print(f"  {roi_name} 已存在，跳過")
            continue

        # 載入 HE 影像
        he_npy = rois_dir / f"{roi_name}_he.npy"
        he_png = rois_dir / f"{roi_name}_he.png"
        if he_npy.exists():
            img = np.load(he_npy)
        elif he_png.exists():
            img = cv2.cvtColor(cv2.imread(str(he_png)), cv2.COLOR_BGR2RGB)
        else:
            print(f"  ⚠️  {roi_name} 找不到 HE 影像，跳過")
            continue

        if img.dtype != np.uint8:
            img = (img * 255).clip(0, 255).astype(np.uint8)

        # CLAHE 增強（clahe_clip=1.0）
        img_proc = _apply_clahe(img, clahe_clip)

        # Cellpose 小直徑
        masks_small, _, _ = model.eval(
            img_proc, diameter=dia_small_px,
            flow_threshold=flow_thr, cellprob_threshold=prob_thr,
            channels=[0, 0]
        )

        # Cellpose 大直徑
        masks_large, _, _ = model.eval(
            img_proc, diameter=dia_large_px,
            flow_threshold=flow_thr, cellprob_threshold=prob_thr,
            channels=[0, 0]
        )

        # LOGIC_A 合併
        mask = _merge_logic_a(masks_small, masks_large, frag_thr)

        # Voronoi 擴張（rna_counting.dilation_px=6，對應 pipeline 實際 RNA 歸屬範圍）
        mask = _voronoi_expand(mask, max_dist=dilation_px)

        np.save(dst, mask.astype(np.uint32))
        print(f"  {roi_name}: {int(mask.max())} 顆細胞 → {dst.name}")


# ═══════════════════════════════════════════════════════════════════════════
# 5. 裸核基礎版（nuclei basic）
# ═══════════════════════════════════════════════════════════════════════════

def prepare_nuc():
    """
    裸核基礎版：Cellpose nuclei 單次偵測，無 LOGIC_A，無 Voronoi dilation。
      - model: nuclei
      - diameter: 15µm（CRC 典型核徑，~54.8px）
      - flow_threshold: 0.4（Cellpose default）
      - cellprob_threshold: 0.0（Cellpose default，較嚴格）
      - 無 LOGIC_A，無 Voronoi 擴張
    用途：對照基準，顯示 LOGIC_A + 參數調整的實際增益
    """
    from cellpose import models
    import cv2

    dia_px   = 15.0 / PX_SIZE  # 15µm → ~54.8px
    flow_thr = 0.4
    prob_thr = 0.0

    print(f"\n[NUC] 裸核基礎版:")
    print(f"  model=nuclei  dia={dia_px:.1f}px  flow={flow_thr}  prob={prob_thr}  dilation=0px")

    model = models.CellposeModel(gpu=True, pretrained_model="nuclei")
    rois_dir = Path(PATHS["crc_he_rois_dir"])

    for roi_name, roi in tqdm(ROI_INFO.items(), desc="NUC Cellpose"):
        dst = MASKS_DIR / f"nuc_{roi_name}.npy"
        if dst.exists():
            print(f"  {roi_name} 已存在，跳過")
            continue

        he_npy = rois_dir / f"{roi_name}_he.npy"
        he_png = rois_dir / f"{roi_name}_he.png"
        if he_npy.exists():
            img = np.load(he_npy)
        elif he_png.exists():
            img = cv2.cvtColor(cv2.imread(str(he_png)), cv2.COLOR_BGR2RGB)
        else:
            print(f"  ⚠️  {roi_name} 找不到 HE 影像，跳過")
            continue

        if img.dtype != np.uint8:
            img = (img * 255).clip(0, 255).astype(np.uint8)

        # 無 CLAHE，無 normalize — 真正的「基礎版」
        mask, _, _ = model.eval(
            img, diameter=dia_px,
            flow_threshold=flow_thr, cellprob_threshold=prob_thr,
            channels=[0, 0]
        )

        np.save(dst, mask.astype(np.uint32))
        print(f"  {roi_name}: {int(mask.max())} 顆細胞 → {dst.name}")


# ═══════════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=["v12", "sr", "cpv3", "p3", "nuc", "all"], default="all")
    args = parser.parse_args()

    if args.method in ("v12", "all"):
        prepare_v12()

    if args.method in ("sr", "all"):
        prepare_spaceranger()

    if args.method in ("cpv3", "all"):
        prepare_cellpose_v3()

    if args.method in ("p3", "all"):
        prepare_p3()

    if args.method in ("nuc", "all"):
        prepare_nuc()

    print("\n✅ 00_prepare_masks.py 完成")
    print(f"   輸出：{MASKS_DIR}")
