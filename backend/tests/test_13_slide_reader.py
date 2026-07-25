"""Test 13: SlideReader 抽象層（backend/src/utils/slide_reader.py）

合成資料為主。NDPI/SVS 走 tifffile 原生的金字塔 series；BTF/TIFF 沿用既有的
tile/strip 讀取路徑。
"""
import numpy as np
import pytest


def _write_tiled_btf(path, size=512, tile=256, seed=0):
    """寫一張未壓縮的 tiled BigTIFF（read_btf_crop 只支援未壓縮 tile）。"""
    import tifffile

    rng = np.random.default_rng(seed)
    img = rng.integers(0, 255, (size, size, 3), dtype=np.uint8)
    tifffile.imwrite(str(path), img, bigtiff=True, tile=(tile, tile), photometric="rgb")
    return img


def _write_pyramid(path, size=512, levels=3, seed=1):
    """寫一張帶金字塔的 tiled TIFF（模擬 SVS/NDPI 的多層結構）。"""
    import tifffile

    rng = np.random.default_rng(seed)
    base = rng.integers(0, 255, (size, size, 3), dtype=np.uint8)
    with tifffile.TiffWriter(str(path), bigtiff=True) as tw:
        tw.write(base, tile=(256, 256), photometric="rgb", subifds=levels - 1)
        img = base
        for _ in range(levels - 1):
            img = img[::2, ::2]
            tw.write(img, tile=(256, 256), photometric="rgb", subfiletype=1)
    return base


# ── 分派 ────────────────────────────────────────────────────────────────────

class TestOpenSlideDispatch:
    """依副檔名選對 reader"""

    @pytest.mark.parametrize("suffix", [".btf", ".tif", ".tiff"])
    def test_tiff_suffixes_use_tiff_reader(self, tmp_path, suffix):
        from backend.src.utils.slide_reader import TiffSlideReader, open_slide

        p = tmp_path / f"img{suffix}"
        _write_tiled_btf(p)

        assert isinstance(open_slide(p), TiffSlideReader)

    @pytest.mark.parametrize("suffix", [".ndpi", ".svs", ".mrxs"])
    def test_slide_suffixes_use_pyramid_reader(self, tmp_path, suffix, monkeypatch):
        from backend.src.utils import slide_reader as sr

        p = tmp_path / f"slide{suffix}"
        _write_pyramid(p)

        reader = sr.open_slide(p)
        assert isinstance(reader, sr.PyramidSlideReader)

    def test_unknown_suffix_raises(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        p = tmp_path / "notes.txt"
        p.write_text("hello", encoding="utf-8")

        with pytest.raises(ValueError) as exc:
            open_slide(p)
        assert "notes.txt" in str(exc.value)

    def test_missing_file_raises(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        with pytest.raises(FileNotFoundError):
            open_slide(tmp_path / "nope.ndpi")


# ── TiffSlideReader（BTF/TIFF）──────────────────────────────────────────────

class TestTiffSlideReader:
    def test_dimensions_are_width_height(self, tmp_path):
        """`dimensions` 依 OpenSlide 慣例回 (W, H)，不是 numpy 的 (H, W)。"""
        from backend.src.utils.slide_reader import open_slide

        import tifffile
        rng = np.random.default_rng(0)
        img = rng.integers(0, 255, (256, 512, 3), dtype=np.uint8)   # H=256, W=512
        p = tmp_path / "img.btf"
        tifffile.imwrite(str(p), img, bigtiff=True, tile=(256, 256), photometric="rgb")

        assert open_slide(p).dimensions == (512, 256)

    def test_read_region_matches_source(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        img = _write_tiled_btf(tmp_path / "img.btf", size=512)
        crop = open_slide(tmp_path / "img.btf").read_region(100, 60, 128, 96)

        assert crop.shape == (96, 128, 3)
        np.testing.assert_array_equal(crop, img[60:156, 100:228])

    def test_level_count_is_one_for_flat_tiff(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        _write_tiled_btf(tmp_path / "img.btf")
        assert open_slide(tmp_path / "img.btf").level_count == 1


# ── PyramidSlideReader（NDPI/SVS/MRXS）──────────────────────────────────────

class TestPyramidSlideReader:
    def test_level_count_and_dimensions(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        _write_pyramid(tmp_path / "slide.svs", size=512, levels=3)
        reader = open_slide(tmp_path / "slide.svs")

        assert reader.dimensions == (512, 512)
        assert reader.level_count == 3
        assert reader.level_dimensions[1] == (256, 256)

    def test_read_region_level0_matches_source(self, tmp_path):
        from backend.src.utils.slide_reader import open_slide

        base = _write_pyramid(tmp_path / "slide.svs", size=512)
        crop = open_slide(tmp_path / "slide.svs").read_region(64, 32, 100, 80)

        assert crop.shape == (80, 100, 3)
        np.testing.assert_array_equal(crop, base[32:112, 64:164])

    def test_read_region_higher_level_uses_pyramid(self, tmp_path):
        """level>0 的座標為**該層自身**的座標系（與 OpenSlide 不同，故明確測）。"""
        from backend.src.utils.slide_reader import open_slide

        base = _write_pyramid(tmp_path / "slide.svs", size=512)
        crop = open_slide(tmp_path / "slide.svs").read_region(10, 20, 64, 48, level=1)

        assert crop.shape == (48, 64, 3)
        np.testing.assert_array_equal(crop, base[::2, ::2][20:68, 10:74])

    def test_read_region_clips_to_bounds(self, tmp_path):
        """越界請求須裁到邊界而非拋錯（tile server 邊緣 tile 會這樣要）。"""
        from backend.src.utils.slide_reader import open_slide

        _write_pyramid(tmp_path / "slide.svs", size=512)
        crop = open_slide(tmp_path / "slide.svs").read_region(480, 480, 128, 128)

        assert crop.shape == (32, 32, 3)

    def test_best_level_for_downsample(self, tmp_path):
        """挑「不小於目標倍率」的最近一層 —— 挑太小的層會放大而糊掉。"""
        from backend.src.utils.slide_reader import open_slide

        _write_pyramid(tmp_path / "slide.svs", size=512, levels=3)
        reader = open_slide(tmp_path / "slide.svs")

        assert reader.best_level_for_downsample(1) == 0
        assert reader.best_level_for_downsample(2) == 1
        assert reader.best_level_for_downsample(3) == 1
        assert reader.best_level_for_downsample(4) == 2
        assert reader.best_level_for_downsample(64) == 2      # 超過金字塔深度 → 最底層
