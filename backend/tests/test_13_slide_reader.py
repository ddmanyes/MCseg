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

NDPI_DIR = Path(
    "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/"
    "20251125_TGIA_VisiumHD/20250612_MQ250428"
)
NDPI_REAL = NDPI_DIR / "D1-2_40X_3_10.53.23.ndpi"

# 同一批掃描的兩張切片。兩檔並存才能確認「6 層 / 3 focal plane / jpegheader」
# 是這台掃描機的一致特性，而非單一檔案的巧合。
NDPI_REAL_FILES = [
    (NDPI_REAL, (119040, 41472)),
    (NDPI_DIR / "D2-4_40X_3_10.55.49.ndpi", (126720, 43776)),
]


def _densest_tissue_xy(reader, block=32):
    """回傳「組織最密區塊」中心的 level 0 座標。

    不可用組織質心 —— 一張切片常有兩塊分離的組織，質心會落在兩者之間的空白，
    讀出來一片近白（std < 3），看起來像 reader 壞了。
    """
    top = reader.level_count - 1
    lw, lh = reader.level_dimensions[top]
    gray = reader.read_region(0, 0, lw, lh, level=top).mean(axis=2)
    tissue = (gray < 200).astype(np.float32)

    ny, nx = lh // block, lw // block
    dens = tissue[:ny * block, :nx * block].reshape(ny, block, nx, block).mean(axis=(1, 3))
    by, bx = np.unravel_index(dens.argmax(), dens.shape)
    ds = reader.dimensions[0] / lw
    return int((bx + 0.5) * block * ds), int((by + 0.5) * block * ds)


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

    @pytest.mark.parametrize(
        ("path", "dims"), NDPI_REAL_FILES,
        ids=[p.stem.split("_")[0] for p, _ in NDPI_REAL_FILES],
    )
    def test_real_ndpi_end_to_end(self, path, dims):
        """真實 Hamamatsu NDPI：尺寸、金字塔、mpp、各層在組織上讀得到紋理。"""
        if not path.exists():
            pytest.skip(f"需要真實 NDPI 檔案：{path.name}")

        from backend.src.utils.slide_reader import open_slide

        r = open_slide(path)

        assert r.dimensions == dims
        assert r.level_count == 6
        assert r.mpp == pytest.approx(0.2305, abs=0.005)   # Hamamatsu 40x
        assert r.level_downsamples == pytest.approx([1, 2, 4, 8, 16, 32], rel=0.01)

        # 各層都要讀得出實際影像（不是空陣列、不是全平）。位置必須挑在組織上 ——
        # 切片大半是空白，隨手挑的座標會 std < 3 而讓這個斷言失去意義。
        cx, cy = _densest_tissue_xy(r)
        for level in range(r.level_count):
            s = r.dimensions[0] / r.level_dimensions[level][0]
            crop = r.read_region(int(cx / s) - 128, int(cy / s) - 128, 256, 256, level=level)
            assert crop.shape == (256, 256, 3), f"L{level} 尺寸不符"
            assert crop.std() > 5, f"L{level} 在組織上仍無紋理（可能回了空/全平資料）"

    @pytest.mark.parametrize(
        ("path", "dims"), NDPI_REAL_FILES,
        ids=[p.stem.split("_")[0] for p, _ in NDPI_REAL_FILES],
    )
    def test_real_ndpi_levels_are_spatially_aligned(self, path, dims):
        """跨層座標一致性：L0 降採樣 2× 必須與同位置的 L1 高度相關。

        這是比「有回東西」強得多的檢查 —— focal plane 軸序或 tile 偏移算錯
        （bug 4 那一類）仍可能回傳形狀正確、有紋理的資料，只是**位置是錯的**。
        實測 D1-2 corr=+0.971、D2-4 corr=+0.991。
        """
        if not path.exists():
            pytest.skip(f"需要真實 NDPI 檔案：{path.name}")

        from backend.src.utils.slide_reader import open_slide

        r = open_slide(path)
        cx, cy = _densest_tissue_xy(r)

        l0 = r.read_region(cx - 512, cy - 512, 1024, 1024, level=0).mean(axis=2)
        l1 = r.read_region(cx // 2 - 256, cy // 2 - 256, 512, 512, level=1).mean(axis=2)
        l0_ds = l0.reshape(512, 2, 512, 2).mean(axis=(1, 3))

        corr = np.corrcoef(l0_ds.ravel(), l1.ravel())[0, 1]
        assert corr > 0.8, f"L0 與 L1 內容不相關（corr={corr:+.3f}）→ 座標映射有誤"


# ── roi/extractor.py::read_image_crop() 對 NDPI 的支援（2026-08-05）────────
#
# 舊版 read_image_crop() 只認 tiled/strip TIFF，NDPI 落進 strip 分支後把
# JPEG 壓縮資料當成未壓縮 bytes reshape 而壞掉。這裡釘住：①NDPI 走新分支不
# 報錯、在組織上讀得到紋理，②既有 tiled BTF 路徑完全不受影響（零回歸）。

class TestReadImageCropNdpi:
    def test_ndpi_crop_reads_real_texture(self):
        """真實 NDPI：read_image_crop() 在組織密集區塊裁得出有紋理的內容。"""
        if not NDPI_REAL.exists():
            pytest.skip(f"需要真實 NDPI 檔案：{NDPI_REAL.name}")

        from backend.src.roi.extractor import read_image_crop
        from backend.src.utils.slide_reader import open_slide

        r = open_slide(NDPI_REAL)
        cx, cy = _densest_tissue_xy(r)

        crop, ax0, ay0 = read_image_crop(NDPI_REAL, cx - 128, cy - 128, 256, 256)

        assert crop.shape == (256, 256, 3)
        assert crop.std() > 5, "NDPI crop 在組織密集區仍無紋理（可能 reshape 錯位/回空資料）"
        assert (ax0, ay0) == (cx - 128, cy - 128), "NDPI 不做 tile 邊界對齊，原點應等於請求值"

    def test_ndpi_crop_respects_margin(self):
        """margin 的換算式對 NDPI 分支與既有 strip 分支必須一致（同一套算式）。"""
        if not NDPI_REAL.exists():
            pytest.skip(f"需要真實 NDPI 檔案：{NDPI_REAL.name}")

        from backend.src.roi.extractor import read_image_crop
        from backend.src.utils.slide_reader import open_slide

        r = open_slide(NDPI_REAL)
        cx, cy = _densest_tissue_xy(r)

        crop, ax0, ay0 = read_image_crop(NDPI_REAL, cx, cy, 128, 128, margin=32)

        assert crop.shape == (128 + 2 * 32, 128 + 2 * 32, 3)
        assert (ax0, ay0) == (cx - 32, cy - 32)

    def test_tiled_btf_crop_unaffected_by_ndpi_branch(self, tmp_path):
        """既有 tiled BTF 路徑不受新增的 NDPI 分支影響（零回歸）。"""
        from backend.src.roi.extractor import read_btf_crop, read_image_crop

        img = _write_tiled_btf(tmp_path / "slide.btf", size=512, tile=256, seed=7)

        via_dispatch = read_image_crop(tmp_path / "slide.btf", 50, 60, 100, 80)
        via_direct = read_btf_crop(tmp_path / "slide.btf", 50, 60, 100, 80)

        assert via_dispatch[1:] == via_direct[1:]          # 對齊後的 origin 一致
        np.testing.assert_array_equal(via_dispatch[0], via_direct[0])


# ── _to_rgb() 高位元深度縮放（2026-08-05，診斷有勝 IF 樣本時發現）───────────
#
# 舊版對非 uint8 輸入直接 .astype(uint8)：16-bit 螢光影像常見數值（例如
# ~40000）會回捲成隨機小數，不是變暗——這裡釘住「線性縮放，不是截斷」，
# 以及「既有 uint8 H&E 路徑完全不受影響」兩件事。

class TestToRgbBitDepth:
    def test_uint16_is_linearly_rescaled_not_truncated(self):
        from backend.src.utils.slide_reader import _to_rgb

        # 40000 / 65535 * 255 ≈ 155；naive astype(uint8) 會回捲成 40000 % 256 = 160
        # ——兩個數字刻意選得接近，用來確認測到的真的是縮放邏輯而非巧合
        arr = np.full((4, 4), 40000, dtype=np.uint16)
        out = _to_rgb(arr)

        assert out.dtype == np.uint8
        assert out.shape == (4, 4, 3)
        expected = round(40000 / 65535 * 255)
        assert out[0, 0, 0] == expected

    def test_uint16_max_value_maps_to_255(self):
        from backend.src.utils.slide_reader import _to_rgb

        arr = np.array([[0, 65535]], dtype=np.uint16)
        out = _to_rgb(arr)

        assert out[0, 0, 0] == 0
        assert out[0, 1, 0] == 255

    def test_uint8_input_unchanged(self):
        """既有 H&E RGB（uint8）路徑必須零回歸——不進入縮放分支。"""
        from backend.src.utils.slide_reader import _to_rgb

        rng = np.random.default_rng(3)
        arr = rng.integers(0, 255, (8, 8, 3), dtype=np.uint8)
        out = _to_rgb(arr)

        np.testing.assert_array_equal(out, arr)

    def test_uint16_grayscale_still_replicated_to_three_channels(self):
        from backend.src.utils.slide_reader import _to_rgb

        arr = np.full((4, 4), 65535, dtype=np.uint16)
        out = _to_rgb(arr)

        assert out.shape == (4, 4, 3)
        assert (out == 255).all()
