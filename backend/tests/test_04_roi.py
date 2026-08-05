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


# ── subset_anndata_roi()：SR fullres px 與 he_image px 非 1:1 時的座標轉換 ─────
#
# 2026-08-05 康育 SDS-D0D1D2 樣本：scalefactors mpp=0.4201（不是專案預設的
# 0.2737），沒有任何 registration json，舊版 subset_anndata_roi 永遠假設兩個
# 座標系同尺度 → 每次都篩到 0 bin。這裡用合成資料重現「軸向縮放 2 倍」的情境，
# 釘住「不傳 config 時退回舊行為、傳 config 時套用 resolve_bin_to_image_transform」
# 兩條路徑，避免未來重構又漂移回naive比較。

def _write_synthetic_scale_sample(tmp_path, scale=2.0):
    """
    造一組「he_image（raw TIFF）／binned_002（scalefactors + hires 圖）」，
    兩者間已知有 `scale` 倍的縮放差（無 registration json，會走
    `resolve_bin_to_mask_scale` 的近似縮放 fallback）。

    he_image 尺寸故意用非正方形（100 W × 200 H），避免只測單一軸就誤放過
    x/y 軸序寫反的 bug。
    """
    import json
    import tifffile
    from PIL import Image

    img_w, img_h = 100, 200          # he_image（raw TIFF）px
    sr_w, sr_h = img_w / scale, img_h / scale   # SR fullres px（推導出的畫布）

    he_path = tmp_path / "he.tiff"
    tifffile.imwrite(str(he_path), np.zeros((img_h, img_w, 3), dtype=np.uint8))

    binned_dir = tmp_path / "binned_002"
    spatial_dir = binned_dir / "spatial"
    spatial_dir.mkdir(parents=True)

    hires_scalef = 0.5   # 任意值；hires 尺寸 ÷ scalef 必須還原出 sr_w/sr_h
    hires_w, hires_h = round(sr_w * hires_scalef), round(sr_h * hires_scalef)
    Image.new("RGB", (hires_w, hires_h)).save(spatial_dir / "tissue_hires_image.png")

    (spatial_dir / "scalefactors_json.json").write_text(json.dumps({
        "microns_per_pixel": 0.4201352316127554,
        "tissue_hires_scalef": hires_scalef,
        "tissue_lowres_scalef": hires_scalef / 10,
    }), encoding="utf-8")

    config = {
        "paths": {
            "he_image": str(he_path),
            "binned_002": str(binned_dir),
        },
        "alignment": {
            "use_alignment_json": True,
            "extra_alignment_json": None,
            "bin_to_mask_scale": None,
            "enabled": False,
            "matrix": [[1, 0, 0], [0, 1, 0]],
        },
    }
    return config


def _make_adata(coords_col_row):
    """造一個只含 obs 空間欄位＋obsm['spatial'] 的最小 AnnData（不需要真的表現量）。"""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    n = len(coords_col_row)
    cols = np.array([c for c, _ in coords_col_row], dtype=float)
    rows = np.array([r for _, r in coords_col_row], dtype=float)
    obs = pd.DataFrame({
        "pxl_col_in_fullres": cols,
        "pxl_row_in_fullres": rows,
    }, index=[f"bc{i}" for i in range(n)])
    adata = ad.AnnData(X=sp.csr_matrix((n, 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = np.stack([cols, rows], axis=1)
    return adata


class TestSubsetAnndataRoiTransform:
    """`subset_anndata_roi` 在 SR fullres px 與 he_image px 非 1:1 時的行為"""

    def test_applies_transform_when_config_given(self, tmp_path):
        """傳入 config 時，ROI 篩選應該用 resolve_bin_to_image_transform 換算後的座標，
        不是直接拿 SR fullres px 跟 ROI（he_image px）比。"""
        from backend.src.roi.extractor import subset_anndata_roi

        config = _write_synthetic_scale_sample(tmp_path, scale=2.0)

        # SR fullres px（scale=2 換算後的 he_image px 在註解標出）
        adata = _make_adata([
            (10, 10),   # → image (20, 20)  在 ROI 內
            (20, 20),   # → image (40, 40)  在 ROI 內（邊界 <45）
            (5, 5),     # → image (10, 10)  在 ROI 外（<15）
            (45, 95),   # → image (90, 190) 遠在 ROI 外
        ])
        roi = {"name": "t", "tissue": "test", "x": 15, "y": 15, "width_px": 30, "height_px": 30}

        sub = subset_anndata_roi(adata, roi, config=config)

        assert sub.n_obs == 2
        assert set(sub.obs_names) == {"bc0", "bc1"}
        # obsm['spatial'] 寫回後應該是 he_image px（下游 counter.py／export 都假設這個座標系）
        got = {tuple(row) for row in sub.obsm["spatial"]}
        assert got == {(20.0, 20.0), (40.0, 40.0)}

    def test_naive_1to1_would_pick_wrong_set(self, tmp_path):
        """反向驗證：不做轉換（config=None）在同一組資料上會篩到不同（錯誤）的集合，
        證明①測試真的在測轉換有沒有生效，②修復前的行為確實是錯的。"""
        from backend.src.roi.extractor import subset_anndata_roi

        adata = _make_adata([(10, 10), (20, 20), (5, 5), (45, 95)])
        roi = {"name": "t", "tissue": "test", "x": 15, "y": 15, "width_px": 30, "height_px": 30}

        sub = subset_anndata_roi(adata, roi, config=None)

        # 沒有轉換：直接拿 SR px 比對 ROI px → 只有 bc1 (20,20) 落在 [15,45)
        assert sub.n_obs == 1
        assert set(sub.obs_names) == {"bc1"}

    def test_identity_scale_matches_naive_result(self, tmp_path):
        """mpp 剛好等於專案預設值（scale≈1）時，傳 config 前後結果應一致 —— 這是
        CRC 官方樣本一路以來「沒出過事」的那個情境，修復不能改變它的行為。"""
        from backend.src.roi.extractor import subset_anndata_roi

        config = _write_synthetic_scale_sample(tmp_path, scale=1.0)
        coords = [(10, 10), (20, 20), (5, 5), (45, 95)]
        roi = {"name": "t", "tissue": "test", "x": 15, "y": 15, "width_px": 30, "height_px": 30}

        with_cfg = subset_anndata_roi(_make_adata(coords), roi, config=config)
        without_cfg = subset_anndata_roi(_make_adata(coords), roi, config=None)

        assert set(with_cfg.obs_names) == set(without_cfg.obs_names) == {"bc1"}
