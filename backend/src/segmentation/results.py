"""
深度產物封裝模組 (SegmentationResult)
統一管理多邊形導出、形態學 QC、巨觀疊圖、1:1 ROI 特寫生成與轉錄組 Bin 歸屬。
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import cv2
import numpy as np
from PIL import Image

logger = logging.getLogger("mcseg.results")


@dataclass
class SegmentationResult:
    """
    全切片 / ROI 分割成果深度物件。
    """
    sample_id: str
    cells: list[dict[str, Any]]  # 每筆含 "cell_id", "polygon", "centroid", "area", "solidity"
    image_shape: tuple[int, int]  # (Height, Width)
    metadata: dict[str, Any] = field(default_factory=dict)
    runtime_seconds: float = 0.0

    def __len__(self) -> int:
        return len(self.cells)

    @property
    def total_cells(self) -> int:
        return len(self.cells)

    def to_geojson(
        self,
        output_path: str | Path,
        target_space: str = "native",
        transformer: Any = None,
    ) -> Path:
        """
        導出 Space Ranger v4.0+ / Xenium RFC 7946 標準 GeoJSON。
        """
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)

        features = []
        for cell in self.cells:
            poly = cell["polygon"]
            if hasattr(poly, "exterior"):
                coords = [list(poly.exterior.coords)]
            elif isinstance(poly, list):
                coords = [poly] if not isinstance(poly[0][0], (list, tuple)) else poly
            else:
                coords = poly

            if target_space != "native" and transformer is not None:
                coords = transformer.transform_polygons([coords], from_space="native_if", to_space=target_space)[0]

            feature = {
                "type": "Feature",
                "id": str(cell.get("cell_id", len(features))),
                "geometry": {
                    "type": "Polygon",
                    "coordinates": coords,
                },
                "properties": {
                    "cell_id": int(cell.get("cell_id", len(features))),
                    "area_px": float(cell.get("area", 0.0)),
                    "solidity": float(cell.get("solidity", 0.0)),
                },
            }
            features.append(feature)

        geojson_obj = {
            "type": "FeatureCollection",
            "features": features,
        }

        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(geojson_obj, f)

        logger.info(f"Saved GeoJSON: {out_p} ({len(features)} cells, {out_p.stat().st_size / (1024*1024):.2f} MB)")
        return out_p

    def compute_qc_metrics(self) -> dict[str, Any]:
        """
        計算單細胞形態學與分佈品質控制指標。
        """
        if not self.cells:
            return {"total_cells": 0}

        areas = np.array([c.get("area", 0.0) for c in self.cells])
        solidities = np.array([c.get("solidity", 0.0) for c in self.cells])

        qc = {
            "sample_id": self.sample_id,
            "total_cells": len(self.cells),
            "area_mean_px": float(np.mean(areas)),
            "area_median_px": float(np.median(areas)),
            "area_std_px": float(np.std(areas)),
            "area_iqr_px": float(np.percentile(areas, 75) - np.percentile(areas, 25)),
            "solidity_mean": float(np.mean(solidities)),
            "solidity_median": float(np.median(solidities)),
            "runtime_seconds": float(self.runtime_seconds),
            "raw_instances": self.metadata.get("raw_instances", len(self.cells)),
            "suppressed_duplicates": self.metadata.get("suppressed_duplicates", 0),
            "dedup_suppression_rate": float(
                self.metadata.get("suppressed_duplicates", 0) / max(1, self.metadata.get("raw_instances", 1))
            ),
        }
        return qc

    def to_summary_json(self, output_path: str | Path) -> Path:
        """導出摘要與 QC JSON。"""
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        summary = self.compute_qc_metrics()
        summary.update(self.metadata)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        return out_p

    def render_macro_overlay(
        self,
        image_path: str | Path,
        output_path: str | Path,
        target_max_dim: int = 4000,
    ) -> Path:
        """
        繪製全玻片巨觀疊圖（綠色細胞邊界）。
        """
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)

        H, W = self.image_shape
        scale = target_max_dim / max(W, H)
        target_w, target_h = int(W * scale), int(H * scale)

        with Image.open(image_path) as img:
            thumb = img.resize((target_w, target_h), Image.Resampling.BILINEAR)
            thumb_np = np.array(thumb)

        if thumb_np.ndim == 2:
            thumb_np = cv2.cvtColor(thumb_np, cv2.COLOR_GRAY2RGB)
        elif thumb_np.shape[2] > 3:
            thumb_np = thumb_np[:, :, :3]

        overlay = thumb_np.copy()
        for cell in self.cells:
            poly = cell["polygon"]
            coords = poly[0] if isinstance(poly, list) and isinstance(poly[0][0], (list, tuple)) else poly
            scaled_poly = (np.array(coords) * scale).astype(np.int32)
            cv2.polylines(overlay, [scaled_poly], isClosed=True, color=(0, 255, 0), thickness=1)

        result = cv2.addWeighted(thumb_np, 0.5, overlay, 0.5, 0)
        Image.fromarray(result).save(out_p, quality=85)
        logger.info(f"Saved Macro Overlay: {out_p}")
        return out_p
