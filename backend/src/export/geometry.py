"""匯出用的幾何轉換：Cellpose mask → GeoJSON 多邊形、座標平移。

不依賴 FastAPI 與 config，可獨立測試。
"""
import logging
from pathlib import Path

logger = logging.getLogger("pipeline.export.geometry")


def mask_to_geojson(
    mask_path: Path,
    pixel_size_um: float,
    min_area_px: int = 20,
) -> dict:
    """
    將 segmentation_masks.npy 轉換為 GeoJSON FeatureCollection。

    座標：ROI 局部 µm（原點 = ROI 左上角），與 cellpose_cells.h5ad obsm['spatial'] 一致。
    使用 regionprops 取 bounding box 後在小 patch 上做輪廓偵測，
    避免 O(n_cells × H×W) 的全圖掃描。
    """
    import numpy as np
    from skimage import measure

    seg_mask = np.load(str(mask_path))

    features = []
    # regionprops 一次性計算 bounding box + area，避免逐細胞全圖掃描
    for prop in measure.regionprops(seg_mask):
        if prop.area < min_area_px:
            continue
        cid = prop.label
        r0, c0, r1, c1 = prop.bbox

        # 在 bounding box patch 上找輪廓（比全圖快數個量級）
        cell_crop = (seg_mask[r0:r1, c0:c1] == cid).astype(np.uint8)
        padded = np.pad(cell_crop, 1, mode="constant")
        contours = measure.find_contours(padded, 0.5)
        if not contours:
            continue

        contour = max(contours, key=len)
        # 還原 padding(1) + bounding box offset，再轉成 (x, y) µm
        xy_um = np.column_stack([
            (contour[:, 1] - 1 + c0) * pixel_size_um,   # col → x
            (contour[:, 0] - 1 + r0) * pixel_size_um,   # row → y
        ])

        if not np.allclose(xy_um[0], xy_um[-1]):
            xy_um = np.vstack([xy_um, xy_um[0]])

        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [xy_um.tolist()],
            },
            "properties": {
                "full_id":   str(int(cid)),
                "cell_id":   int(cid),
            },
        })

    logger.info(f"  生成 {len(features)} 個 Cellpose 多邊形")
    return {"type": "FeatureCollection", "features": features}


def shift_geojson_coords(feat: dict, dx: float, dy: float) -> None:
    """In-place 平移 GeoJSON feature 的座標。"""
    def _shift(coords):
        if not coords:
            return coords
        if isinstance(coords[0], (int, float)):
            return [coords[0] + dx, coords[1] + dy] + list(coords[2:])
        return [_shift(c) for c in coords]

    geom = feat.get("geometry", {})
    if geom and geom.get("coordinates") is not None:
        geom["coordinates"] = _shift(geom["coordinates"])
