"""
空間座標系深模組 (SpatialCoordinateTransformer)
提供 Visium HD / 顯微鏡多模態（CytAssist, raw H&E TIFF, TissueFAX IF）之雙向座標投影與旋轉對齊。

設計準則：
1. 封裝所有 scale factor 計算（scalefactors_json.json）、TIFF 標頭、90° 旋轉與仿射矩陣。
2. 消除呼叫端手動乘/除 scale 的淺模組摩擦。
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Sequence, Union

import numpy as np

logger = logging.getLogger("mcseg.spatial_coords")

CoordinateSpace = Literal["virtual_fullres", "raw_tiff", "native_if", "tissue_hires", "tissue_lowres"]


@dataclass
class SpatialCoordinateTransformer:
    """
    多模態空間座標投影器。
    
    支援空間：
    - virtual_fullres: 10x CytAssist 捕獲窗像素空間（例如 3200x3000）
    - raw_tiff / native_if: 顯微鏡原生高解析像素空間（例如 27000x27000 或 13706x45045）
    - tissue_hires: 10x 縮圖空間
    """
    vfr_width: float = 0.0
    vfr_height: float = 0.0
    raw_width: float = 0.0
    raw_height: float = 0.0
    hires_scalef: float = 1.0
    lowres_scalef: float = 1.0
    rotation_deg: int = 0  # 支援 0, 90, 180, 270 順時針旋轉
    scale_raw_to_vfr: float = field(init=False)
    scale_vfr_to_raw: float = field(init=False)

    def __post_init__(self):
        # 規範化旋轉角度 (0, 90, 180, 270)
        self.rotation_deg = self.rotation_deg % 360
        if self.rotation_deg not in (0, 90, 180, 270):
            raise ValueError(f"rotation_deg 必須為 90 的倍數 (0, 90, 180, 270)，得到: {self.rotation_deg}")

        if self.vfr_width > 0 and self.raw_width > 0:
            scale_w = self.raw_width / self.vfr_width
            scale_h = self.raw_height / self.vfr_height if self.vfr_height > 0 else scale_w
            self.scale_vfr_to_raw = (scale_w + scale_h) / 2.0
            self.scale_raw_to_vfr = 1.0 / self.scale_vfr_to_raw
        else:
            self.scale_vfr_to_raw = 1.0
            self.scale_raw_to_vfr = 1.0

    @classmethod
    def from_visium_dir(
        cls,
        binned_dir: str | Path,
        image_path: str | Path | None = None,
        rotation_deg: int = 0,
    ) -> SpatialCoordinateTransformer:
        """從 Visium HD 輸出目錄與影像檔案自動解析所有尺度因數與維度。"""
        from PIL import Image

        binned_path = Path(binned_dir)
        spatial_dir = binned_path / "spatial" if (binned_path / "spatial").exists() else binned_path

        sf_path = spatial_dir / "scalefactors_json.json"
        hires_path = spatial_dir / "tissue_hires_image.png"

        hires_scalef = 1.0
        lowres_scalef = 1.0
        vfr_w, vfr_h = 0.0, 0.0

        if sf_path.exists():
            with open(sf_path) as f:
                sf = json.load(f)
            hires_scalef = float(sf.get("tissue_hires_scalef", 1.0))
            lowres_scalef = float(sf.get("tissue_lowres_scalef", 1.0))

            if hires_path.exists():
                with Image.open(hires_path) as img:
                    w_h, h_h = img.size
                vfr_w = w_h / hires_scalef
                vfr_h = h_h / hires_scalef

        raw_w, raw_h = 0.0, 0.0
        if image_path is not None:
            img_p = Path(image_path)
            if img_p.exists():
                suffix = img_p.suffix.lower()
                if suffix in (".tif", ".tiff", ".btf"):
                    try:
                        import tifffile
                        with tifffile.TiffFile(str(img_p)) as tf:
                            page = tf.pages[0]
                            raw_h, raw_w = float(page.imagelength), float(page.imagewidth)
                    except Exception:
                        with Image.open(img_p) as img:
                            raw_w, raw_h = float(img.width), float(img.height)
                else:
                    with Image.open(img_p) as img:
                        raw_w, raw_h = float(img.width), float(img.height)

        logger.info(
            f"Initialized SpatialCoordinateTransformer: VFR ({vfr_w:.0f}x{vfr_h:.0f}) -> "
            f"Raw ({raw_w:.0f}x{raw_h:.0f}) | Scale: {raw_w/vfr_w if vfr_w > 0 else 1.0:.4f} | Rot: {rotation_deg}°"
        )

        return cls(
            vfr_width=vfr_w,
            vfr_height=vfr_h,
            raw_width=raw_w,
            raw_height=raw_h,
            hires_scalef=hires_scalef,
            lowres_scalef=lowres_scalef,
            rotation_deg=rotation_deg,
        )

    def _get_scale_factor(self, from_space: str, to_space: str) -> float:
        """取得兩個座標系之間的縮放純量。"""
        if from_space == to_space:
            return 1.0

        # 將 from 轉為 vfr 的比例
        to_vfr = {
            "virtual_fullres": 1.0,
            "raw_tiff": self.scale_raw_to_vfr,
            "native_if": self.scale_raw_to_vfr,
            "tissue_hires": 1.0 / self.hires_scalef if self.hires_scalef > 0 else 1.0,
            "tissue_lowres": 1.0 / self.lowres_scalef if self.lowres_scalef > 0 else 1.0,
        }

        # 將 vfr 轉為 to 的比例
        from_vfr = {
            "virtual_fullres": 1.0,
            "raw_tiff": self.scale_vfr_to_raw,
            "native_if": self.scale_vfr_to_raw,
            "tissue_hires": self.hires_scalef,
            "tissue_lowres": self.lowres_scalef,
        }

        if from_space not in to_vfr or to_space not in from_vfr:
            raise ValueError(f"未知的座標空間轉換: {from_space} -> {to_space}")

        return to_vfr[from_space] * from_vfr[to_space]

    def transform_points(
        self,
        points: np.ndarray | Sequence[Sequence[float]],
        from_space: CoordinateSpace = "native_if",
        to_space: CoordinateSpace = "virtual_fullres",
        apply_rotation: bool = True,
    ) -> np.ndarray:
        """
        轉換點座標陣列 (N, 2) [X, Y]。
        """
        pts = np.asarray(points, dtype=np.float64)
        if pts.size == 0:
            return pts

        is_1d = (pts.ndim == 1)
        if is_1d:
            pts = pts.reshape(1, -1)

        scale = self._get_scale_factor(from_space, to_space)
        transformed = pts * scale

        if apply_rotation and self.rotation_deg != 0:
            # 在目標空間中處理旋轉
            # 目標畫布長寬
            if to_space in ("raw_tiff", "native_if"):
                W, H = self.raw_width, self.raw_height
            else:
                W, H = self.vfr_width, self.vfr_height

            x = transformed[:, 0]
            y = transformed[:, 1]

            if self.rotation_deg == 90:
                # 順時針 90度: (x, y) -> (H - y, x)
                new_x = H - y
                new_y = x
            elif self.rotation_deg == 180:
                # 180度: (x, y) -> (W - x, H - y)
                new_x = W - x
                new_y = H - y
            elif self.rotation_deg == 270:
                # 順時針 270度: (x, y) -> (y, W - x)
                new_x = y
                new_y = W - x
            else:
                new_x, new_y = x, y

            transformed = np.column_stack([new_x, new_y])

        return transformed[0] if is_1d else transformed

    def transform_polygons(
        self,
        polygons: list[list[Sequence[float]] | dict | Any],
        from_space: CoordinateSpace = "native_if",
        to_space: CoordinateSpace = "virtual_fullres",
    ) -> list[list[list[float]]]:
        """
        轉換多邊形頂點串列。
        支援 GeoJSON coordinates, Polygon 物件或原生頂點串列。
        """
        transformed_polys = []
        for poly in polygons:
            if hasattr(poly, "exterior"):  # Shapely Polygon
                coords = np.array(poly.exterior.coords)
            elif isinstance(poly, dict) and "coordinates" in poly:
                coords = np.array(poly["coordinates"][0])
            elif isinstance(poly, list):
                coords = np.array(poly[0] if isinstance(poly[0][0], (list, tuple)) else poly)
            else:
                coords = np.array(poly)

            t_coords = self.transform_points(coords, from_space, to_space)
            transformed_polys.append([t_coords.tolist()])

        return transformed_polys

    def transform_bounding_box(
        self,
        bbox: tuple[float, float, float, float] | Sequence[float],
        from_space: CoordinateSpace = "native_if",
        to_space: CoordinateSpace = "virtual_fullres",
    ) -> tuple[float, float, float, float]:
        """
        轉換 (min_x, min_y, max_x, max_y) 邊界框。
        """
        min_x, min_y, max_x, max_y = bbox
        corners = np.array([
            [min_x, min_y],
            [max_x, min_y],
            [max_x, max_y],
            [min_x, max_y]
        ])
        t_corners = self.transform_points(corners, from_space, to_space)
        return (
            float(t_corners[:, 0].min()),
            float(t_corners[:, 1].min()),
            float(t_corners[:, 0].max()),
            float(t_corners[:, 1].max()),
        )
