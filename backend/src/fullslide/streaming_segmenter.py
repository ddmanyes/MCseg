"""
統一全切片串流分割引擎 (StreamingSlideSegmenter)
合流 H&E 與 IF 雙軌切片處理，統一 50x 縮圖遮罩過濾、PIL Lazy Windowed 切塊、MPS/CUDA GPU 調度與 KDTree 去重。
"""

from __future__ import annotations

import gc
import logging
import time
from pathlib import Path
from typing import Any, Callable, Literal, Sequence

import cv2
import numpy as np
import torch
from cellpose import models
from PIL import Image
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
from shapely.ops import unary_union

from backend.src.segmentation.results import SegmentationResult

# 解除 PIL 大圖像素限制
Image.MAX_IMAGE_PIXELS = None

logger = logging.getLogger("mcseg.streaming_segmenter")


class StreamingSlideSegmenter:
    """
    全切片串流分割引擎。
    """

    def __init__(
        self,
        tile_size: int = 1024,
        stride: int = 900,
        downsample_factor: int = 50,
        tissue_threshold: float = 0.01,
        device: Literal["auto", "mps", "cuda", "cpu"] = "auto",
        model_type: str = "cpsam",
    ):
        self.tile_size = tile_size
        self.stride = stride
        self.downsample_factor = downsample_factor
        self.tissue_threshold = tissue_threshold
        self.model_type = model_type

        # 硬體加速檢測
        self.device, self.use_gpu = self._setup_device(device)
        self.model = None

    def _setup_device(self, requested_device: str) -> tuple[torch.device, bool]:
        if requested_device == "auto":
            if torch.backends.mps.is_available():
                dev = torch.device("mps")
                use_gpu = True
                logger.info("Apple Silicon MPS detected and enabled.")
            elif torch.cuda.is_available():
                dev = torch.device("cuda")
                use_gpu = True
                logger.info("CUDA GPU detected and enabled.")
            else:
                dev = torch.device("cpu")
                use_gpu = False
                logger.info("Using CPU fallback.")
        elif requested_device in ("mps", "cuda"):
            dev = torch.device(requested_device)
            use_gpu = True
        else:
            dev = torch.device("cpu")
            use_gpu = False
        return dev, use_gpu

    def _load_model(self):
        if self.model is None:
            logger.info(f"Loading Cellpose {self.model_type} model on {self.device}...")
            self.model = models.CellposeModel(gpu=self.use_gpu, model_type=self.model_type)
            logger.info("Cellpose model loaded successfully.")

    def generate_thumbnail_tissue_mask(
        self, image_path: str | Path
    ) -> tuple[np.ndarray, int, int]:
        """
        在 0.5 秒內生成 50x 縮圖組織前景遮罩。
        """
        with Image.open(image_path) as img:
            W, H = img.size
            thumb_w = max(10, W // self.downsample_factor)
            thumb_h = max(10, H // self.downsample_factor)
            thumb = img.resize((thumb_w, thumb_h), Image.Resampling.BILINEAR)
            thumb_np = np.array(thumb)

        if thumb_np.ndim == 2:
            gray = thumb_np
        elif thumb_np.shape[2] >= 3:
            # 兼容 IF 螢光與 H&E
            gray = np.maximum.reduce([thumb_np[:, :, 0], thumb_np[:, :, 1], thumb_np[:, :, 2]])
        else:
            gray = thumb_np[:, :, 0]

        blurred = cv2.GaussianBlur(gray, (5, 5), 0)
        _, thresh = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        fg_mask = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)

        tissue_ratio = np.count_nonzero(fg_mask) / fg_mask.size
        logger.info(
            f"Thumbnail Mask: {W}x{H} px | Tissue Coverage: {tissue_ratio*100:.1f}% | "
            f"Background Skipped: {(1-tissue_ratio)*100:.1f}%"
        )
        return fg_mask, W, H

    def _is_tissue_tile(
        self, fg_mask: np.ndarray, x: int, y: int
    ) -> bool:
        ds = self.downsample_factor
        gx1, gx2 = x // ds, (x + self.tile_size) // ds
        gy1, gy2 = y // ds, (y + self.tile_size) // ds
        tile_fg = fg_mask[gy1:gy2, gx1:gx2]
        if tile_fg.size == 0:
            return False
        return (np.count_nonzero(tile_fg) / tile_fg.size) >= self.tissue_threshold

    def segment_slide(
        self,
        image_path: str | Path,
        sample_id: str = "sample",
        stain_type: Literal["if", "he"] = "if",
        skip_background: bool = True,
        max_tiles: int | None = None,
    ) -> SegmentationResult:
        """
        執行全切片串流分割。
        """
        start_t = time.time()
        self._load_model()
        image_path = Path(image_path)

        fg_mask, W, H = self.generate_thumbnail_tissue_mask(image_path)
        tiles_x = (W + self.stride - 1) // self.stride
        tiles_y = (H + self.stride - 1) // self.stride
        total_grid = tiles_x * tiles_y

        all_cells = []
        centroids = []
        active_tiles = 0
        skipped_tiles = 0

        logger.info(f"Starting Streaming Slide Segmentation: {sample_id} ({total_grid} total grid tiles)")

        with Image.open(image_path) as slide_img:
            for ty in range(tiles_y):
                for tx in range(tiles_x):
                    x = tx * self.stride
                    y = ty * self.stride

                    if skip_background and not self._is_tissue_tile(fg_mask, x, y):
                        skipped_tiles += 1
                        continue

                    active_tiles += 1
                    crop_w = min(self.tile_size, W - x)
                    crop_h = min(self.tile_size, H - y)
                    tile_crop = slide_img.crop((x, y, x + crop_w, y + crop_h))
                    tile_np = np.array(tile_crop)

                    if tile_np.shape[0] != self.tile_size or tile_np.shape[1] != self.tile_size:
                        pad_img = np.zeros((self.tile_size, self.tile_size, tile_np.shape[2]), dtype=tile_np.dtype)
                        pad_img[:crop_h, :crop_w] = tile_np
                        tile_np = pad_img

                    # Cellpose 推論
                    masks, _, _ = self.model.eval(
                        tile_np,
                        diameter=None,
                        channels=[0, 0],
                        augment=False,
                        min_size=15,
                    )

                    unique_ids = np.unique(masks)
                    for uid in unique_ids:
                        if uid == 0:
                            continue
                        cell_mask = (masks == uid).astype(np.uint8)
                        contours, _ = cv2.findContours(cell_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                        if not contours:
                            continue
                        cnt = contours[0]
                        if len(cnt) < 4:
                            continue

                        M = cv2.moments(cnt)
                        if M["m00"] == 0:
                            continue
                        cx = (M["m10"] / M["m00"]) + x
                        cy = (M["m01"] / M["m00"]) + y

                        if cx >= W or cy >= H:
                            continue

                        global_cnt = cnt.squeeze() + np.array([x, y])
                        if global_cnt.ndim != 2 or len(global_cnt) < 4:
                            continue

                        area = cv2.contourArea(cnt)
                        hull = cv2.convexHull(cnt)
                        hull_area = cv2.contourArea(hull)
                        solidity = float(area / hull_area) if hull_area > 0 else 1.0

                        all_cells.append({
                            "cell_id": len(all_cells),
                            "polygon": global_cnt.tolist(),
                            "centroid": [float(cx), float(cy)],
                            "area": float(area),
                            "solidity": solidity,
                        })
                        centroids.append([float(cx), float(cy)])

                    if max_tiles and active_tiles >= max_tiles:
                        break
                if max_tiles and active_tiles >= max_tiles:
                    break

        raw_instances = len(all_cells)
        logger.info(f"Tile processing complete: {active_tiles} active tiles, {skipped_tiles} skipped, {raw_instances} raw cells.")

        # KDTree 全域空間去重
        if len(centroids) > 0:
            logger.info("Running KDTree spatial deduplication...")
            kdtree = cKDTree(np.array(centroids))
            pairs = kdtree.query_pairs(r=8.0)
            suppressed = set()
            for idx1, idx2 in pairs:
                if idx1 in suppressed or idx2 in suppressed:
                    continue
                # 保留面積較大者
                if all_cells[idx1]["area"] < all_cells[idx2]["area"]:
                    suppressed.add(idx1)
                else:
                    suppressed.add(idx2)

            unique_cells = [c for i, c in enumerate(all_cells) if i not in suppressed]
            n_suppressed = len(suppressed)
        else:
            unique_cells = []
            n_suppressed = 0

        runtime_s = time.time() - start_t
        logger.info(f"KDTree complete: Retained {len(unique_cells)} unique cells (suppressed {n_suppressed}). Runtime: {runtime_s:.2f}s")

        return SegmentationResult(
            sample_id=sample_id,
            cells=unique_cells,
            image_shape=(H, W),
            metadata={
                "raw_instances": raw_instances,
                "suppressed_duplicates": n_suppressed,
                "active_tiles": active_tiles,
                "skipped_tiles": skipped_tiles,
                "stain_type": stain_type,
            },
            runtime_seconds=runtime_s,
        )
