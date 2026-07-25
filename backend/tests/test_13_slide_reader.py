"""Test 13: SlideReader 抽象層（backend/src/utils/slide_reader.py）

合成資料為主。NDPI/SVS 走 tifffile 原生的金字塔 series；BTF/TIFF 沿用既有的
tile/strip 讀取路徑。
"""
from pathlib import Path

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

    def test_read_region_handles_jpeg_tiles(self, tmp_path):
        """SVS/NDPI 的 tile 是 JPEG 壓縮的 —— 解碼路徑必須走得通（含 JPEG tables）。"""
        import tifffile
        from scipy.ndimage import gaussian_filter

        from backend.src.utils.slide_reader import open_slide

        rng = np.random.default_rng(9)
        smooth = gaussian_filter(rng.random((512, 512)), 8)
        base = np.repeat((smooth * 255).astype(np.uint8)[:, :, None], 3, axis=2)

        p = tmp_path / "slide.ndpi"
        with tifffile.TiffWriter(str(p), bigtiff=True) as tw:
            tw.write(base, tile=(256, 256), photometric="rgb", compression="jpeg", subifds=1)
            tw.write(base[::2, ::2], tile=(256, 256), photometric="rgb",
                     compression="jpeg", subfiletype=1)

        crop = open_slide(p).read_region(100, 100, 128, 128)

        assert crop.shape == (128, 128, 3)
        # JPEG 有損 → 比對平均值而非逐位元
        assert abs(float(crop.mean()) - float(base[100:228, 100:228].mean())) < 5.0

    def test_mpp_from_resolution_tag(self, tmp_path):
        """NDPI/SVS 通常以 CENTIMETER 記錄解析度 → 換算成 µm/px。"""
        import tifffile

        from backend.src.utils.slide_reader import open_slide

        p = tmp_path / "slide.svs"
        img = np.zeros((256, 256, 3), dtype=np.uint8)
        # 44248 px/cm ≈ 0.226 µm/px（Hamamatsu 40x）
        tifffile.imwrite(str(p), img, tile=(256, 256), photometric="rgb",
                         resolution=(44248, 44248), resolutionunit="CENTIMETER")

        assert open_slide(p).mpp == pytest.approx(0.226, abs=0.001)


# ── tile server 整合 ────────────────────────────────────────────────────────

class TestTileServerUsesSlideReader:
    """DZITileServer 改由 SlideReader 取得尺寸與影像資料"""

    def test_dzi_dimensions_from_reader(self, tmp_path):
        from backend.src.roi.tile_server import DZITileServer

        _write_tiled_btf(tmp_path / "img.btf", size=512)
        srv = DZITileServer(str(tmp_path / "img.btf"))

        assert (srv.full_width, srv.full_height) == (512, 512)
        assert 'Width="512"' in srv.get_dzi()

    def test_get_tile_returns_jpeg(self, tmp_path):
        from backend.src.roi.tile_server import DZITileServer

        _write_tiled_btf(tmp_path / "img.btf", size=512)
        srv = DZITileServer(str(tmp_path / "img.btf"))
        data = srv.get_tile(srv.max_level, 0, 0)

        assert data[:2] == b"\xff\xd8"          # JPEG SOI

    def test_pyramid_slide_serves_tiles(self, tmp_path):
        """NDPI/SVS 也要能出 tile（過去只支援 BTF）。"""
        from backend.src.roi.tile_server import DZITileServer

        _write_pyramid(tmp_path / "slide.svs", size=512, levels=3)
        srv = DZITileServer(str(tmp_path / "slide.svs"))

        assert (srv.full_width, srv.full_height) == (512, 512)
        assert srv.get_tile(srv.max_level, 0, 0)[:2] == b"\xff\xd8"

    def test_thumb_uses_pyramid_when_available(self, tmp_path):
        """自帶金字塔時，縮圖直接取自金字塔層，不逐 tile 掃描全圖。"""
        from backend.src.roi.tile_server import _load_or_build_thumb, _thumb_from_pyramid

        _write_pyramid(tmp_path / "slide.svs", size=512, levels=3)

        fast = _thumb_from_pyramid(tmp_path / "slide.svs", 4)
        assert fast is not None and fast.shape[:2] == (128, 128)

        thumb = _load_or_build_thumb(tmp_path / "slide.svs", 4)
        assert thumb.shape[:2] == (128, 128)
        assert (tmp_path / "slide.thumb4.npy").exists()      # 已快取

    def test_thumb_falls_back_for_flat_tiff(self, tmp_path):
        """單層 BTF 沒有金字塔 → 回 None，改走既有的逐 tile 掃描。"""
        from backend.src.roi.tile_server import _load_or_build_thumb, _thumb_from_pyramid

        _write_tiled_btf(tmp_path / "img.btf", size=512)

        assert _thumb_from_pyramid(tmp_path / "img.btf", 4) is None
        assert _load_or_build_thumb(tmp_path / "img.btf", 4).shape[:2] == (128, 128)


# ── 資料掃描 ────────────────────────────────────────────────────────────────

class TestDiscoveryRecognisesSlideFormats:
    """/api/data/scan 要認得 NDPI/SVS/MRXS"""

    def test_ndpi_is_discovered_without_size_threshold(self, tmp_path):
        """NDPI 本身就是全片掃描檔 —— 不該像一般 TIFF 那樣要求 > 100MB。"""
        from backend.src.utils.discovery import scan_data_root

        (tmp_path / "scan.ndpi").write_bytes(b"\x00" * 1024)

        result = scan_data_root(str(tmp_path))

        assert result.he_image is not None
        assert "scan.ndpi" in result.he_image.label
        assert "NDPI" in result.he_image.label

    def test_small_plain_tiff_is_still_ignored(self, tmp_path):
        """小型 TIFF 多半是縮圖／輔助圖，維持既有的體積門檻。"""
        from backend.src.utils.discovery import scan_data_root

        (tmp_path / "thumb.tif").write_bytes(b"\x00" * 1024)

        assert scan_data_root(str(tmp_path)).he_image is None

    def test_svs_preferred_over_large_plain_tiff(self, tmp_path):
        """全片掃描格式優先於同目錄下的大型一般 TIFF。"""
        from backend.src.utils.discovery import scan_data_root

        (tmp_path / "big.tif").write_bytes(b"\x00" * (101 * 1024 * 1024))
        (tmp_path / "slide.svs").write_bytes(b"\x00" * 2048)

        result = scan_data_root(str(tmp_path))

        assert "slide.svs" in result.he_image.label
        assert any("big.tif" in f.label for f in result.extra_files)


# ── 真實 NDPI 的兩個陷阱 ────────────────────────────────────────────────────

NDPI_REAL = Path(
    "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/"
    "20251125_TGIA_VisiumHD/20250612_MQ250428/D1-2_40X_3_10.53.23.ndpi"
)


class TestNdpiQuirks:
    """真實 NDPI 與一般金字塔 TIFF 的兩處結構差異

    2026-07-25 以真實檔（119040×41472、6 層、2.17 GB）驗證時抓到：

    1. **series 有 focal plane 前導軸**（`ZYXS`，shape `(3, 41472, 119040, 3)`）
       —— 以 `series.shape[0:2]` 取尺寸會得到 `(41472, 3)`，而且
       `read_region` 只會安靜地回空陣列，不報錯。
    2. **segment 是無標頭的 JPEG scan data** —— `decode()` 不傳 `jpegheader`
       會得到 `Jpeg8Error: Not a JPEG file`。
    """

    def test_leading_axis_series_uses_page_dimensions(self, tmp_path):
        """series 帶前導軸時，尺寸仍須取自 page（合成 3-plane 檔重現）。"""
        import tifffile

        from backend.src.utils.slide_reader import open_slide

        rng = np.random.default_rng(11)
        plane = rng.integers(0, 255, (256, 512, 3), dtype=np.uint8)   # H=256, W=512
        p = tmp_path / "stack.svs"
        with tifffile.TiffWriter(str(p), bigtiff=True) as tw:
            for _ in range(3):        # 3 個 focal plane → series 取得前導軸
                tw.write(plane, tile=(256, 256), photometric="rgb")

        reader = open_slide(p)

        assert reader.dimensions == (512, 256), (
            f"取到 {reader.dimensions}，可能誤用了 series.shape 的前導軸"
        )
        assert reader.read_region(0, 0, 64, 32).shape == (32, 64, 3)

    def test_jpegheader_is_passed_to_decode(self):
        """page 帶 jpegheader 時必須轉交解碼器（NDPI 的硬需求）。"""
        from backend.src.utils import slide_reader as sr

        seen: dict = {}

        class _FakePage:
            is_tiled = True
            tilelength, tilewidth = 8, 16
            imagelength, imagewidth = 8, 16
            samplesperpixel = 3
            dtype = np.uint8
            dataoffsets = (0,)
            databytecounts = (4,)
            jpegheader = b"\xff\xd8FAKEHEADER"

            @property
            def decode(self):
                def _decode(data, index, **kw):
                    seen.update(kw)
                    return np.zeros((1, 8, 16, 3), np.uint8), (0,), (1, 8, 16, 3)
                return _decode

        class _FakeFH:
            def seek(self, *_): pass
            def read(self, n): return b"\x00" * n

        sr._read_page_region(_FakeFH(), _FakePage(), 0, 0, 16, 8)

        assert seen.get("jpegheader") == b"\xff\xd8FAKEHEADER"

    def test_no_jpegheader_kw_for_plain_tiff(self, tmp_path):
        """一般 TIFF 沒有 jpegheader 屬性時不得硬塞該參數。"""
        from backend.src.utils import slide_reader as sr

        seen: list = []

        class _FakePage:
            is_tiled = True
            tilelength, tilewidth = 8, 16
            imagelength, imagewidth = 8, 16
            samplesperpixel = 3
            dtype = np.uint8
            dataoffsets = (0,)
            databytecounts = (4,)
            jpegheader = None

            @property
            def decode(self):
                def _decode(data, index, **kw):
                    seen.append(kw)
                    return np.zeros((1, 8, 16, 3), np.uint8), (0,), (1, 8, 16, 3)
                return _decode

        class _FakeFH:
            def seek(self, *_): pass
            def read(self, n): return b"\x00" * n

        sr._read_page_region(_FakeFH(), _FakePage(), 0, 0, 16, 8)

        assert seen == [{}]

    @pytest.mark.skipif(not NDPI_REAL.exists(), reason="需要真實 NDPI 檔案")
    def test_real_ndpi_end_to_end(self):
        """真實 Hamamatsu NDPI：尺寸、金字塔、mpp、各層區域讀取。"""
        from backend.src.utils.slide_reader import open_slide

        r = open_slide(NDPI_REAL)

        assert r.dimensions == (119040, 41472)
        assert r.level_count == 6
        assert r.mpp == pytest.approx(0.2305, abs=0.005)   # Hamamatsu 40x

        # 各層都要讀得出實際影像（不是空陣列、不是全平）
        for level in (0, 2, 5):
            crop = r.read_region(25000 >> level, 9000 >> level, 256, 256, level=level)
            assert crop.shape == (256, 256, 3)
            assert crop.std() > 1
