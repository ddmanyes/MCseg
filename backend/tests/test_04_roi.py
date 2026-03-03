"""Test 4: ROI Extractor — BTF tile-based 裁切測試"""
import pytest
import numpy as np
from pathlib import Path

from .conftest import CRC_BTF, CRC_BINNED_002


class TestBtfTileRead:
    """BTF tile-based 讀取（CLAUDE.md §5 規範）"""

    def test_read_btf_metadata(self):
        """讀取 BTF metadata（不載入影像）"""
        import tifffile
        with tifffile.TiffFile(CRC_BTF) as tf:
            page = tf.pages[0]
            # Just verify the page metadata is accessible
            assert page.shape[0] > 10000
            assert page.shape[1] > 10000
            assert len(page.shape) >= 2  # at least H x W

    def test_btf_tile_properties(self):
        """BTF 有 tile 結構，可用於 tile-based 讀取"""
        import tifffile
        with tifffile.TiffFile(CRC_BTF) as tf:
            page = tf.pages[0]
            assert page.is_tiled, "BTF should be a tiled TIFF"
            tw, th = page.tilewidth, page.tilelength
            assert tw > 0 and th > 0
            print(f"  Tile: {tw}x{th}, Image: {page.shape}")


class TestRoiExtractorImport:
    """ROI Extractor 模組可匯入"""

    def test_import_extractor(self):
        from backend.src.roi.extractor import RoiExtractor
        assert RoiExtractor is not None

    def test_extractor_init(self, config_dict):
        """Extractor 可用 config 初始化"""
        from backend.src.roi.extractor import RoiExtractor
        ext = RoiExtractor(config_dict)
        assert ext is not None


class TestScalefactors:
    """Scalefactors 讀取與驗證"""

    def test_read_scalefactors(self):
        """讀取 CRC 2µm spatial 的 scalefactors"""
        import json
        sf_path = CRC_BINNED_002 / "spatial" / "scalefactors_json.json"
        with open(sf_path) as f:
            sf = json.load(f)
        # Visium HD 應有 microns_per_pixel
        mpp = sf.get("microns_per_pixel", None)
        if mpp:
            assert 0.1 < mpp < 1.0, f"Unexpected microns_per_pixel: {mpp}"
            print(f"  microns_per_pixel = {mpp}")

    def test_crc_pixel_size(self):
        """CRC pixel_size 應為 ~0.2737"""
        import json
        sf_path = CRC_BINNED_002 / "spatial" / "scalefactors_json.json"
        with open(sf_path) as f:
            sf = json.load(f)
        mpp = sf.get("microns_per_pixel", 0)
        assert abs(mpp - 0.2737) < 0.01, (
            f"Expected ~0.2737, got {mpp}"
        )
