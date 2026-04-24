"""
DZI (Deep Zoom Image) tile server for OpenSeadragon ROI viewer.

Serves tiles from a pyramidal BTF/TIFF file:
- High zoom (scale < 4): reads directly from BTF via read_btf_crop (tile-based, memory-safe)
- Low zoom  (scale >= 4): samples from a pre-built raw-TIFF thumbnail (correct coordinate system)

Previous design sampled from tissue_hires_image.png (Space Ranger's registered image), whose
scalef is calibrated to the Space Ranger internal fullres coordinate system — NOT the raw TIFF.
That mismatch caused overview tiles to show a misregistered view, so ROI coordinates drawn on
the overview didn't map to the same location in the raw TIFF crop.

Fix: build a thumbnail directly from the raw TIFF at init time so all zoom levels share the
same coordinate system. The thumbnail is cached to disk (alongside the TIFF) on first build.
"""
from __future__ import annotations

import io
import logging
import math
from pathlib import Path
from typing import Optional

import numpy as np
from PIL import Image

from .extractor import read_btf_crop

logger = logging.getLogger("pipeline.roi")

TILE_SIZE   = 256
OVERLAP     = 1
THUMB_SCALE = 32   # 1/32 downsample for overview: 21504×47104 → ~672×1472 px


class DZITileServer:
    """Serves Deep Zoom Image tiles from a pyramidal BTF/TIFF."""

    full_width:   int
    full_height:  int
    max_level:    int
    _thumb_arr:   Optional[np.ndarray]
    _thumb_scale: int

    def __init__(
        self,
        btf_path:   str,
        _hires_path: Optional[str] = None,   # kept for API compat, no longer used
        _scalef:     float = 0.1,           # kept for API compat, no longer used
    ):
        self.btf_path     = Path(btf_path)
        self._thumb_scale = THUMB_SCALE
        self._thumb_arr   = None

        import tifffile
        with tifffile.TiffFile(str(btf_path)) as tf:
            page = tf.pages[0]
            self.full_height = page.imagelength
            self.full_width  = page.imagewidth

        self.max_level = math.ceil(math.log2(max(self.full_width, self.full_height)))

        # Build/load thumbnail from raw TIFF (single consistent coordinate system)
        try:
            self._thumb_arr = _load_or_build_thumb(self.btf_path, THUMB_SCALE)
        except Exception as e:
            logger.warning(f"縮圖建立失敗，低倍概覽將停用：{e}")

    # ── DZI descriptor ────────────────────────────────────────────────────────

    def get_dzi(self) -> str:
        return (
            '<?xml version="1.0" encoding="utf-8"?>\n'
            f'<Image TileSize="{TILE_SIZE}" Overlap="{OVERLAP}" Format="jpeg"'
            ' xmlns="http://schemas.microsoft.com/deepzoom/2008">\n'
            f'  <Size Width="{self.full_width}" Height="{self.full_height}"/>\n'
            '</Image>'
        )

    # ── Tile generation ────────────────────────────────────────────────────────

    def get_tile(self, level: int, tx: int, ty: int) -> bytes:
        """Return JPEG bytes for DZI tile at (level, tx, ty)."""
        scale = 2 ** (self.max_level - level)

        x0 = max(0, tx * TILE_SIZE * scale - OVERLAP * scale)
        y0 = max(0, ty * TILE_SIZE * scale - OVERLAP * scale)
        x1 = min(self.full_width,  (tx * TILE_SIZE + TILE_SIZE + OVERLAP) * scale)
        y1 = min(self.full_height, (ty * TILE_SIZE + TILE_SIZE + OVERLAP) * scale)
        w, h = x1 - x0, y1 - y0

        if w <= 0 or h <= 0:
            return self._blank_tile()

        target_w = math.ceil(w / scale)
        target_h = math.ceil(h / scale)

        if scale >= 4 and self._thumb_arr is not None:
            crop = _crop_from_thumb(self._thumb_arr, self._thumb_scale,
                                    x0, y0, x1, y1)
        else:
            crop, _, _ = read_btf_crop(self.btf_path, x0, y0, w, h)

        img = Image.fromarray(crop).convert('RGB')
        if img.size != (target_w, target_h):
            img = img.resize((target_w, target_h), Image.LANCZOS)

        buf = io.BytesIO()
        img.save(buf, 'JPEG', quality=85)
        return buf.getvalue()

    # ── Helpers ────────────────────────────────────────────────────────────────

    @staticmethod
    def _blank_tile() -> bytes:
        buf = io.BytesIO()
        Image.new('RGB', (1, 1), (20, 20, 20)).save(buf, 'JPEG')
        return buf.getvalue()


# ── Thumbnail builder (module-level so it can be called independently) ────────

def _load_or_build_thumb(btf_path: Path, scale: int) -> np.ndarray:
    """
    Load a cached raw-TIFF thumbnail, or build and cache it.

    The thumbnail preserves the raw TIFF coordinate system: thumb pixel (thx, thy)
    corresponds to fullres pixel (thx*scale, thy*scale).  This is different from the
    Space Ranger hires image, which is registered to SpaceRanger's internal fullres
    coordinate space (and therefore has a different origin and scale factor).
    """
    cache = btf_path.with_suffix(f".thumb{scale}.npy")
    if cache.exists():
        logger.info(f"載入縮圖快取：{cache}")
        return np.load(str(cache))

    logger.info(f"建立 raw TIFF 縮圖（scale=1/{scale}）…  首次執行需數秒，之後載入快取")

    import tifffile

    with tifffile.TiffFile(str(btf_path)) as tf:
        page      = tf.pages[0]
        H         = page.imagelength
        W         = page.imagewidth
        TH        = getattr(page, "tilelength", 512)
        TW        = getattr(page, "tilewidth",  512)
        n_tiles_x = (W + TW - 1) // TW
        n_tiles_y = (H + TH - 1) // TH

        offsets_tag    = page.tags.get("TileOffsets")
        bytecounts_tag = page.tags.get("TileByteCounts")

        if not offsets_tag or not bytecounts_tag:
            # Non-tiled TIFF: full load (warn: may be slow / large)
            logger.warning("TIFF 無 TileOffsets，全圖載入以建立縮圖（可能耗時）")
            img   = page.asarray()
            thumb = np.ascontiguousarray(img[::scale, ::scale])
            del img
            np.save(str(cache), thumb)
            return thumb

        offsets    = offsets_tag.value
        bytecounts = bytecounts_tag.value

    thumb_h = (H + scale - 1) // scale
    thumb_w = (W + scale - 1) // scale
    thumb   = np.zeros((thumb_h, thumb_w, 3), dtype=np.uint8)

    with open(str(btf_path), "rb") as fh:
        for ty in range(n_tiles_y):
            for tx in range(n_tiles_x):
                tidx = ty * n_tiles_x + tx
                if tidx >= len(offsets) or offsets[tidx] == 0:
                    continue

                fh.seek(offsets[tidx])
                raw = np.frombuffer(fh.read(bytecounts[tidx]), np.uint8)
                if raw.size != TH * TW * 3:
                    logger.warning(
                        f"Tile ({tx},{ty}) 大小不符（{raw.size} != {TH*TW*3}），跳過。"
                        "如影像非 uint8 RGB，縮圖可能不完整。"
                    )
                    continue
                tile = raw.reshape(TH, TW, 3)

                # Fullres region this tile covers
                fy0, fx0 = ty * TH, tx * TW
                fy1 = min(H, fy0 + TH)
                fx1 = min(W, fx0 + TW)

                # Sample positions: multiples of `scale` that fall in this tile
                gy_start = ((fy0 + scale - 1) // scale) * scale
                gx_start = ((fx0 + scale - 1) // scale) * scale
                if gy_start >= fy1 or gx_start >= fx1:
                    continue

                gy = np.arange(gy_start, fy1, scale)
                gx = np.arange(gx_start, fx1, scale)

                ly = gy - fy0      # local rows within tile
                lx = gx - fx0      # local cols within tile
                thy = gy // scale  # thumbnail row indices
                thx = gx // scale  # thumbnail col indices

                # Clamp to tile bounds (edge tiles may be padded but data ends earlier)
                valid_y = ly < TH
                valid_x = lx < TW
                ly, thy = ly[valid_y], thy[valid_y]
                lx, thx = lx[valid_x], thx[valid_x]

                # Vectorised assignment
                ly2d, lx2d   = np.ix_(ly,  lx)
                thy2d, thx2d = np.ix_(thy, thx)
                thumb[thy2d, thx2d] = tile[ly2d, lx2d]

    np.save(str(cache), thumb)
    logger.info(f"縮圖快取已儲存：{cache}  shape={thumb.shape}")
    return thumb


def _crop_from_thumb(
    thumb: np.ndarray,
    scale: int,
    x0: int, y0: int,
    x1: int, y1: int,
) -> np.ndarray:
    """Crop a region from the raw-TIFF thumbnail using fullres coordinates."""
    hx0 = max(0, int(x0 / scale))
    hy0 = max(0, int(y0 / scale))
    hx1 = min(thumb.shape[1], math.ceil(x1 / scale))
    hy1 = min(thumb.shape[0], math.ceil(y1 / scale))
    crop = thumb[hy0:hy1, hx0:hx1]
    if crop.size == 0:
        return np.zeros((1, 1, 3), dtype=np.uint8)
    return crop
