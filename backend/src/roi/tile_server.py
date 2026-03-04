"""
DZI (Deep Zoom Image) tile server for OpenSeadragon ROI viewer.

Serves tiles from a pyramidal BTF/TIFF file:
- High zoom (scale < 4): reads directly from BTF via read_btf_crop (tile-based, memory-safe)
- Low zoom  (scale >= 4): resamples from preloaded hires PNG to avoid huge fullres reads
"""
from __future__ import annotations

import io
import math
from pathlib import Path
from typing import Optional

import numpy as np
from PIL import Image

from .extractor import read_btf_crop

TILE_SIZE = 256
OVERLAP   = 1


class DZITileServer:
    """Serves Deep Zoom Image tiles from a pyramidal BTF/TIFF."""

    full_width:  int
    full_height: int
    max_level:   int
    _hires_arr:  Optional[np.ndarray]
    _scalef:     float   # hires px / fullres px

    def __init__(
        self,
        btf_path:   str,
        hires_path: Optional[str] = None,
        scalef:     float = 0.1,
    ):
        self.btf_path = Path(btf_path)
        self._scalef  = scalef
        self._hires_arr = None

        import tifffile
        with tifffile.TiffFile(str(btf_path)) as tf:
            page = tf.pages[0]
            self.full_height = page.imagelength
            self.full_width  = page.imagewidth

        self.max_level = math.ceil(math.log2(max(self.full_width, self.full_height)))

        if hires_path:
            p = Path(hires_path)
            if p.exists():
                self._hires_arr = np.array(Image.open(p).convert('RGB'))

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
        # Number of fullres pixels per DZI pixel at this level
        scale = 2 ** (self.max_level - level)

        # Tile region in fullres pixels (with overlap on each side)
        x0 = max(0, tx * TILE_SIZE * scale - OVERLAP * scale)
        y0 = max(0, ty * TILE_SIZE * scale - OVERLAP * scale)
        x1 = min(self.full_width,  (tx * TILE_SIZE + TILE_SIZE + OVERLAP) * scale)
        y1 = min(self.full_height, (ty * TILE_SIZE + TILE_SIZE + OVERLAP) * scale)
        w, h = x1 - x0, y1 - y0

        if w <= 0 or h <= 0:
            return self._blank_tile()

        target_w = math.ceil(w / scale)
        target_h = math.ceil(h / scale)

        if scale >= 4 and self._hires_arr is not None:
            # Low zoom: sample from preloaded hires image
            crop = self._crop_from_hires(x0, y0, x1, y1)
        else:
            # High zoom: read tile-by-tile from BTF
            crop, _, _ = read_btf_crop(self.btf_path, x0, y0, w, h)

        img = Image.fromarray(crop).convert('RGB')
        if img.size != (target_w, target_h):
            img = img.resize((target_w, target_h), Image.LANCZOS)

        buf = io.BytesIO()
        img.save(buf, 'JPEG', quality=85)
        return buf.getvalue()

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _crop_from_hires(self, x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
        """Convert fullres coords to hires coords and crop."""
        sf  = self._scalef
        arr = self._hires_arr
        hx0 = max(0, int(x0 * sf))
        hy0 = max(0, int(y0 * sf))
        hx1 = min(arr.shape[1], math.ceil(x1 * sf))
        hy1 = min(arr.shape[0], math.ceil(y1 * sf))
        crop = arr[hy0:hy1, hx0:hx1]
        # Guard against zero-size crop
        if crop.size == 0:
            return np.zeros((1, 1, 3), dtype=np.uint8)
        return crop

    @staticmethod
    def _blank_tile() -> bytes:
        buf = io.BytesIO()
        Image.new('RGB', (1, 1), (20, 20, 20)).save(buf, 'JPEG')
        return buf.getvalue()
