"""
SlideReader 抽象層 —— 統一 BTF/TIFF 與 NDPI/SVS/MRXS 的區域讀取
================================================================

**為何不是 openslide-python**：openslide 需要 `brew install openslide`（macOS）
或手動放 DLL（Windows），與本專案「ExFAT 隨身碟、macOS ⇄ Windows 兩台機器」
的情境衝突。

**為何也不是 tiffslide**：實測 `tiffslide 2.5.0` 在本專案的
`tifffile 2026.2.24` 上直接 ImportError（`cannot import name 'ZarrTiffStore'`
—— 該 API 已被移除）。降級 tifffile 會影響整條既有的 BTF 讀取路徑，代價太大。

**採用**：`tifffile` 原生支援 —— 它以副檔名 `.ndpi` 辨識 NDPI 並提供
`_series_ndpi` / `_ndpi_load_pages`，SVS/MRXS 則走一般的 subIFD 金字塔。
`series[0].levels` 即各層影像，可直接切片讀取區域。

座標慣例（與 OpenSlide 相同）：

- `dimensions` 回 **(width, height)**（不是 numpy 的 (H, W)）
- `read_region(x, y, w, h, level)` 的 `(x, y)` 為**該 level 自身**的座標系
  （OpenSlide 是 level-0 座標；此處刻意不同，因為 MSseg 的呼叫端都是先算好
  該層的座標再讀，換算兩次反而容易出錯 —— 測試中明確釘死此語意）
- 回傳 `(h, w, 3)` uint8 RGB（丟棄 alpha）
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np

logger = logging.getLogger("pipeline.utils.slide_reader")

# 走既有 tile/strip 讀取路徑（無金字塔）
FLAT_TIFF_SUFFIXES = frozenset({".btf", ".tif", ".tiff"})
# 走 tifffile 金字塔 series（自帶多解析度層）
PYRAMID_SUFFIXES = frozenset({".ndpi", ".svs", ".mrxs"})


@runtime_checkable
class SlideReader(Protocol):
    """全片影像讀取器的最小介面。"""

    @property
    def dimensions(self) -> tuple[int, int]:
        """level 0 的 `(width, height)`。"""
        ...

    @property
    def level_count(self) -> int:
        ...

    @property
    def mpp(self) -> float | None:
        """µm/px；影像未記載時為 None（**不要**猜測預設值）。"""
        ...

    def read_region(
        self, x: int, y: int, w: int, h: int, level: int = 0
    ) -> np.ndarray:
        ...


def open_slide(path: str | Path) -> SlideReader:
    """依副檔名分派 reader。"""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"找不到影像檔：{path.name}")

    suffix = path.suffix.lower()
    if suffix in FLAT_TIFF_SUFFIXES:
        return TiffSlideReader(path)
    if suffix in PYRAMID_SUFFIXES:
        return PyramidSlideReader(path)
    raise ValueError(
        f"不支援的影像格式：{path.name}"
        f"（支援 {sorted(FLAT_TIFF_SUFFIXES | PYRAMID_SUFFIXES)}）"
    )


def _to_rgb(arr: np.ndarray) -> np.ndarray:
    """統一成 (h, w, 3) uint8：灰階補成三通道、RGBA 丟棄 alpha、高位元深度線性縮放。

    螢光顯微影像（IF：DAPI/marker channel）常見 16-bit——直接 `.astype(uint8)`
    是**截斷／回捲**，不是縮放，數值 > 255 會變成隨機雜訊而不是變暗
    （2026-08-05 診斷有勝 IF 樣本時發現：這裡原本就是這樣直接截斷）。
    非 uint8 的整數型別一律先依 dtype 的理論最大值線性縮放到 0–255 再轉型；
    既有 uint8 H&E 資料完全不受影響（`dtype != np.uint8` 才會進入這段）。
    """
    arr = np.asarray(arr)
    if arr.dtype != np.uint8:
        if np.issubdtype(arr.dtype, np.integer):
            max_val = float(np.iinfo(arr.dtype).max)
        else:
            max_val = float(arr.max()) or 1.0
        arr = np.round(arr.astype(np.float32) / max_val * 255.0).clip(0, 255).astype(np.uint8)
    if arr.ndim == 2:
        arr = np.repeat(arr[:, :, None], 3, axis=2)
    elif arr.shape[2] > 3:
        arr = arr[:, :, :3]
    return arr.astype(np.uint8, copy=False)


class TiffSlideReader:
    """BTF / 一般 TIFF：包裝既有的 `read_btf_crop`（tiled）與 `read_strip_crop`。"""

    def __init__(self, path: str | Path):
        import tifffile

        self.path = Path(path)
        with tifffile.TiffFile(str(self.path)) as tf:
            page = tf.pages[0]
            self._height = int(page.imagelength)
            self._width = int(page.imagewidth)
            self._is_tiled = bool(page.tags.get("TileOffsets"))

    @property
    def dimensions(self) -> tuple[int, int]:
        return self._width, self._height

    @property
    def level_dimensions(self) -> list[tuple[int, int]]:
        return [self.dimensions]

    @property
    def level_count(self) -> int:
        return 1

    @property
    def mpp(self) -> float | None:
        # 這類檔案的 XResolution 實測為 96 dpi 的通用預設值（≈264 µm/px，
        # 對組織切片顯然是假的）→ 一律回 None，由 scalefactors 決定。
        return None

    def best_level_for_downsample(self, downsample: float) -> int:
        return 0

    def read_region(self, x: int, y: int, w: int, h: int, level: int = 0) -> np.ndarray:
        if level != 0:
            raise ValueError(f"{self.path.name} 無金字塔，僅有 level 0")

        from backend.src.roi.extractor import read_btf_crop
        from backend.src.roi.tile_server import read_strip_crop

        if self._is_tiled:
            crop, _, _ = read_btf_crop(self.path, x, y, w, h)
            return _to_rgb(crop)
        return _to_rgb(read_strip_crop(self.path, x, y, w, h))


class PyramidSlideReader:
    """NDPI / SVS / MRXS：走 tifffile 的金字塔 series，逐層切片讀取。"""

    def __init__(self, path: str | Path):
        import tifffile

        self.path = Path(path)
        self._tf = tifffile.TiffFile(str(self.path))
        series = self._tf.series[0]
        # 沒有 subIFD 金字塔的檔案，levels 只有一層 —— 仍可正常運作
        self._levels = list(getattr(series, "levels", None) or [series])
        # 尺寸取自 page 而非 series.shape：真實 NDPI 的 series 是 'ZYXS'
        # （3 個 focal plane 的前導軸，實測 shape=(3, 41472, 119040, 3)），
        # 拿 shape[0]/shape[1] 會得到 (41472, 3) 這種荒謬尺寸，而且
        # read_region 只會安靜地回空陣列。page 的 imagewidth/imagelength 永遠是對的。
        self._dims = [
            (int(lv.pages[0].imagewidth), int(lv.pages[0].imagelength))   # (W, H)
            for lv in self._levels
        ]

    def close(self) -> None:
        self._tf.close()

    def __enter__(self) -> "PyramidSlideReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    @property
    def dimensions(self) -> tuple[int, int]:
        return self._dims[0]

    @property
    def level_dimensions(self) -> list[tuple[int, int]]:
        return list(self._dims)

    @property
    def level_count(self) -> int:
        return len(self._levels)

    @property
    def level_downsamples(self) -> list[float]:
        w0 = self._dims[0][0]
        return [w0 / w for w, _ in self._dims]

    @property
    def mpp(self) -> float | None:
        """由 TIFF resolution tag 推導 µm/px；NDPI/SVS 通常記載為 CENTIMETER。"""
        import tifffile

        page = self._tf.pages[0]
        try:
            res = page.tags["XResolution"].value
            unit = page.tags["ResolutionUnit"].value
        except (KeyError, AttributeError):
            return None

        value = res[0] / res[1] if isinstance(res, tuple) else float(res)
        if not value:
            return None
        if unit == tifffile.RESUNIT.CENTIMETER:
            return 10_000.0 / value          # px/cm → µm/px
        if unit == tifffile.RESUNIT.INCH:
            return 25_400.0 / value          # px/inch → µm/px
        return None

    def best_level_for_downsample(self, downsample: float) -> int:
        """
        挑出**不小於**目標倍率的最近一層。

        挑過小的層會被放大而糊掉；寧可挑解析度較高的層再自行降採樣。
        目標倍率超過金字塔深度時回最底層。
        """
        best = 0
        for i, ds in enumerate(self.level_downsamples):
            if ds <= downsample + 1e-9:
                best = i
        return best

    def read_region(self, x: int, y: int, w: int, h: int, level: int = 0) -> np.ndarray:
        if not 0 <= level < self.level_count:
            raise ValueError(f"level {level} 超出範圍（0..{self.level_count - 1}）")

        lw, lh = self._dims[level]
        x0, y0 = max(0, int(x)), max(0, int(y))
        x1, y1 = min(lw, x0 + int(w)), min(lh, y0 + int(h))
        if x1 <= x0 or y1 <= y0:
            return np.zeros((0, 0, 3), dtype=np.uint8)

        page = self._levels[level].pages[0]
        return _to_rgb(_read_page_region(self._tf.filehandle, page, x0, y0, x1, y1))


def _read_page_region(fh, page, x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
    """
    只解壓被請求區域覆蓋到的 segment（tile 或 strip），組出該區域。

    刻意**不走 tifffile 的 `aszarr()`** —— 它要求 `zarr >= 3`，而本專案的
    anndata/scanpy 依賴仍在 zarr 2.x（實測 `ValueError: zarr 2.18.7 < 3 is
    not supported`）。`page.decode()` 是同一份解碼器（含 JPEG tables 處理），
    直接用它即可，不必動整條依賴鏈。

    也**不可**用 `page.asarray()`：那會把整層載入 RAM，NDPI/SVS 的 level 0
    動輒數十 GB。
    """
    offsets = page.dataoffsets
    bytecounts = page.databytecounts
    if not offsets:
        raise ValueError(f"{getattr(page, 'name', 'page')} 沒有影像資料 offsets，無法讀取區域")

    samples = page.samplesperpixel
    height, width = int(page.imagelength), int(page.imagewidth)

    if page.is_tiled:
        th, tw = int(page.tilelength), int(page.tilewidth)
    else:
        th, tw = int(getattr(page, "rowsperstrip", height) or height), width
    n_across = (width + tw - 1) // tw

    # ⚠️ 必須在迴圈外先取出 `page.decode`。它是 cached_property，第一次取值會
    # 讀取 JPEG tables **而移動檔案指標**（實測 +10 bytes）。若寫成
    # `page.decode(fh.read(...), i)`，Python 會先解析屬性、再讀資料 ——
    # 於是第一塊 tile 從 offset+10 開始讀，得到「Not a JPEG file」。
    decode = page.decode

    # NDPI 的 segment 是**沒有標頭的 JPEG scan data**，必須把 page 的 jpegheader
    # 一起餵進解碼器；不傳會得到 `Jpeg8Error: Not a JPEG file`。
    # 一般 tiled TIFF / SVS 沒有這個屬性，傳 None 即可。
    decode_kw = {}
    jpegheader = getattr(page, "jpegheader", None)
    if jpegheader:
        decode_kw["jpegheader"] = jpegheader

    out = np.zeros((y1 - y0, x1 - x0, samples), dtype=page.dtype)
    for ty in range(y0 // th, (y1 + th - 1) // th):
        for tx in range(x0 // tw, (x1 + tw - 1) // tw):
            index = ty * n_across + tx
            if index >= len(offsets) or not bytecounts[index]:
                continue                       # 稀疏 TIFF 的空白 segment
            fh.seek(offsets[index])
            seg, _, shape = decode(fh.read(bytecounts[index]), index, **decode_kw)
            seg = np.asarray(seg).reshape(shape)[0]   # (depth, h, w, s) → (h, w, s)

            sy0, sx0 = ty * th, tx * tw
            # 取 segment 與請求區域的交集，兩邊各自換算成局部索引
            iy0, iy1 = max(y0, sy0), min(y1, sy0 + seg.shape[0])
            ix0, ix1 = max(x0, sx0), min(x1, sx0 + seg.shape[1])
            if iy1 <= iy0 or ix1 <= ix0:
                continue
            out[iy0 - y0:iy1 - y0, ix0 - x0:ix1 - x0] = seg[
                iy0 - sy0:iy1 - sy0, ix0 - sx0:ix1 - sx0
            ]
    return out
