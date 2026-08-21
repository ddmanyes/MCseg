"""Test 10: 全片流程共用模組（backend/src/fullslide/pipeline.py）

全部使用合成資料 —— 不需真實 CRC/BTF 資料，Windows 亦可執行。
"""
from pathlib import Path

import numpy as np
import pytest


# ── 合成資料輔助 ─────────────────────────────────────────────────────────────

def _make_mask() -> np.ndarray:
    """10×10 遮罩：label 1 佔 rows 0-4、label 2 佔 rows 5-9。"""
    mask = np.zeros((10, 10), dtype=np.int32)
    mask[0:5, :] = 1
    mask[5:10, :] = 2
    return mask


def _write_tissue_positions(path, rows_cols, in_tissue=None):
    """寫出 tissue_positions.parquet；rows_cols 為 [(row, col), ...]。"""
    import pandas as pd

    n = len(rows_cols)
    if in_tissue is None:
        in_tissue = [1] * n
    df = pd.DataFrame({
        "barcode": [f"BC{i}" for i in range(n)],
        "in_tissue": in_tissue,
        "pxl_row_in_fullres": [rc[0] for rc in rows_cols],
        "pxl_col_in_fullres": [rc[1] for rc in rows_cols],
    })
    df.to_parquet(str(path), index=False)
    return df


# ── bin_attribution ─────────────────────────────────────────────────────────

class TestBinAttribution:
    """2µm bins → cell_id 對應"""

    def test_bin_attribution_maps_bins_to_cells(self, tmp_path):
        """4 個 bin 分別落在 label 1 / label 2 區域，須對應正確 cell_id。"""
        from backend.src.fullslide.pipeline import bin_attribution

        mask = _make_mask()
        tp = tmp_path / "tissue_positions.parquet"
        # BC0,BC1 → rows 1,3（label 1）；BC2,BC3 → rows 6,8（label 2）
        _write_tissue_positions(tp, [(1, 2), (3, 4), (6, 2), (8, 4)])

        attr = bin_attribution(mask, tp, crop_y0=0, crop_x0=0)

        assert list(attr.columns) == ["barcode", "cell_id"]
        mapping = dict(zip(attr["barcode"], attr["cell_id"]))
        assert mapping == {"BC0": 1, "BC1": 1, "BC2": 2, "BC3": 2}

    def test_bin_attribution_excludes_out_of_tissue(self, tmp_path):
        """in_tissue=0 的 bin 必須被排除。"""
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(1, 1), (2, 2)], in_tissue=[1, 0])

        attr = bin_attribution(_make_mask(), tp, crop_y0=0, crop_x0=0)

        assert attr["barcode"].tolist() == ["BC0"]

    def test_bin_attribution_applies_crop_origin(self, tmp_path):
        """裁切原點須自 fullres 座標扣除。"""
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        # fullres (row=106, col=202) - origin (100, 200) = local (6, 2) → label 2
        _write_tissue_positions(tp, [(106, 202)])

        attr = bin_attribution(_make_mask(), tp, crop_y0=100, crop_x0=200)

        assert attr["cell_id"].tolist() == [2]

    def test_bin_attribution_drops_background_bins(self, tmp_path):
        """落在 cell_id=0（背景）的 bin 不應出現在結果中。"""
        from backend.src.fullslide.pipeline import bin_attribution

        mask = np.zeros((10, 10), dtype=np.int32)
        mask[0:2, 0:2] = 1          # 只有左上角是細胞
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(0, 0), (9, 9)])   # 後者落在背景

        attr = bin_attribution(mask, tp, crop_y0=0, crop_x0=0)

        assert attr["barcode"].tolist() == ["BC0"]

    def test_bin_attribution_applies_scale(self, tmp_path):
        """SR fullres px → 遮罩 px 的縮放須在扣除裁切原點**之前**套用。

        真實情境：dpcp01 樣本的 SR fullres 是 0.5464 µm/px，raw TIFF 約 0.2737，
        遮罩在 TIFF 空間 → bin 座標必須先乘 ~2 才對得上。
        """
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        # SR 座標 (row=3, col=1) × scale 2 → 遮罩 (6, 2) → label 2
        _write_tissue_positions(tp, [(3, 1)])

        attr = bin_attribution(_make_mask(), tp, 0, 0, scale=(2.0, 2.0))

        assert attr["cell_id"].tolist() == [2]

    def test_bin_attribution_scale_is_per_axis(self, tmp_path):
        """縮放須可分軸指定（真實樣本 col 1.9088 vs row 2.0000 並非等向）。"""
        from backend.src.fullslide.pipeline import bin_attribution

        mask = np.zeros((10, 10), dtype=np.int32)
        mask[8, 2] = 7
        tp = tmp_path / "tp.parquet"
        # SR (row=4, col=2) × (scale_x=1.0, scale_y=2.0) → 遮罩 (8, 2)
        _write_tissue_positions(tp, [(4, 2)])

        attr = bin_attribution(mask, tp, 0, 0, scale=(1.0, 2.0))

        assert attr["cell_id"].tolist() == [7]

    def test_bin_attribution_scale_applied_before_crop_origin(self, tmp_path):
        """順序必須是「先縮放、後減原點」——原點是遮罩空間的量。"""
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        # SR (row=53, col=51) × 2 = (106, 102)；減原點 (100, 100) → (6, 2) → label 2
        _write_tissue_positions(tp, [(53, 51)])

        attr = bin_attribution(_make_mask(), tp, crop_y0=100, crop_x0=100, scale=(2.0, 2.0))

        assert attr["cell_id"].tolist() == [2]

    def test_bin_attribution_excludes_out_of_bounds_instead_of_clamping(self, tmp_path):
        """越界的 bin 須被排除，不可夾到邊緣（否則會誤記到邊界細胞）。"""
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        # (5, 5) 在界內（label 2）；(500, 500) 遠在界外
        _write_tissue_positions(tp, [(5, 5), (500, 500)])

        attr = bin_attribution(_make_mask(), tp, 0, 0)

        assert attr["barcode"].tolist() == ["BC0"]

    def test_bin_attribution_with_homography_equals_scale_for_diagonal(self, tmp_path):
        """對角 homography 與同值 `scale` 路徑須產生完全相同的結果。

        這保證 P0.5 引入 homography 後**不破壞既有呼叫端** —— 現行的分軸縮放
        只是 3×3 的對角特例。
        """
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(3, 1), (1, 2), (4, 4), (0, 0)])
        mask = _make_mask()

        by_scale = bin_attribution(mask, tp, 0, 0, scale=(1.9, 2.0))
        by_matrix = bin_attribution(
            mask, tp, 0, 0, transform=np.diag([1.9, 2.0, 1.0])
        )

        assert by_scale["barcode"].tolist() == by_matrix["barcode"].tolist()
        assert by_scale["cell_id"].tolist() == by_matrix["cell_id"].tolist()

    def test_bin_attribution_homography_applies_translation(self, tmp_path):
        """homography 的平移項須作用於 (x, y) = (col, row)，不可軸序顛倒。"""
        from backend.src.fullslide.pipeline import bin_attribution

        mask = np.zeros((10, 10), dtype=np.int32)
        mask[7, 3] = 5
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(2, 1)])   # (row=2, col=1)

        # x = col + 2 = 3、y = row + 5 = 7 → 命中 label 5
        m = np.array([[1.0, 0.0, 2.0], [0.0, 1.0, 5.0], [0.0, 0.0, 1.0]])
        attr = bin_attribution(mask, tp, 0, 0, transform=m)

        assert attr["cell_id"].tolist() == [5]

    def test_bin_attribution_transform_overrides_scale(self, tmp_path):
        """兩者同時給定時以 transform 為準（幾何正確者優先）。"""
        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(3, 1)])   # ×2 → (6, 2) → label 2

        attr = bin_attribution(
            _make_mask(), tp, 0, 0,
            scale=(50.0, 50.0),                       # 這個會讓 bin 飛出界
            transform=np.diag([2.0, 2.0, 1.0]),
        )
        assert attr["cell_id"].tolist() == [2]

    def test_coverage_attrs_reports_out_of_image_fraction(self, tmp_path):
        """`.attrs['coverage']` 須如實回報落在影像外的比例。

        真實情境：dpcp01 在正確變換下，SR 右緣約 980 TIFF px 寬的 bin 確實不在
        高解析影像內。那是資料的真實限制，必須回報而非靠壓縮硬塞進去。
        """
        from backend.src.fullslide.pipeline import bin_attribution

        mask = _make_mask()
        tp = tmp_path / "tp.parquet"
        # 4 個 in_tissue bin：2 個在界內（其中 1 個落在細胞上）、2 個界外
        mask_bg = mask.copy()
        mask_bg[0:5, :] = 0                  # label 1 區改為背景
        _write_tissue_positions(tp, [(6, 1), (1, 1), (900, 900), (-500, 0)])

        attr = bin_attribution(mask_bg, tp, 0, 0)
        cov = attr.attrs["coverage"]

        assert cov["n_total"] == 4
        assert cov["n_in_bounds"] == 2
        assert cov["n_assigned"] == 1        # 界內但落在背景的那個不計入
        assert cov["frac_out_of_image"] == pytest.approx(0.5)

    def test_bin_attribution_writes_cache_when_out_path_given(self, tmp_path):
        """給 out_path 時須寫出 parquet，且內容與回傳一致。"""
        import pandas as pd

        from backend.src.fullslide.pipeline import bin_attribution

        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(1, 1)])
        out = tmp_path / "attr.parquet"

        attr = bin_attribution(_make_mask(), tp, 0, 0, out_path=out)

        assert out.exists()
        pd.testing.assert_frame_equal(pd.read_parquet(str(out)), attr)


# ── aggregate_cells ─────────────────────────────────────────────────────────

def _write_h5(path, barcodes, gene_names, counts):
    """寫出可被 sc.read_10x_h5 讀取的替身：直接寫 h5ad 並回傳路徑。

    aggregate_cells 內部用 scanpy 讀 10x h5；測試改以 monkeypatch 注入
    AnnData，避免建構 10x HDF5 二進位格式。
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    adata = ad.AnnData(
        X=sp.csr_matrix(np.asarray(counts, dtype=np.float32)),
    )
    adata.obs_names = list(barcodes)
    adata.var_names = list(gene_names)
    adata.write_h5ad(str(path))
    return path


class TestHitRateSanityCheck:
    """命中率健全性檢查

    回答「分割的圖是不是 Space Ranger 的輸入圖」。判斷錯了**不會報錯** ——
    錯的座標仍讓 bin 落在影像界內，只是落在錯的細胞上，越界檢查完全抓不到。

    2026-07-25 實測：正確慣例命中 20.9–52.9%，錯誤慣例 0.8–3.9%。
    """

    @staticmethod
    def _dense_mask(size=600, cell=8, pitch=16):
        """鋪滿方形細胞的遮罩，正常座標下命中率遠高於門檻。"""
        m = np.zeros((size, size), dtype=np.int32)
        lbl = 1
        for y in range(0, size - cell, pitch):
            for x in range(0, size - cell, pitch):
                m[y:y + cell, x:x + cell] = lbl
                lbl += 1
        return m

    @staticmethod
    def _bins_on_cells(n=5000, cell=8, pitch=16, size=600, offset=0):
        """落在細胞中心的 bin 座標；`offset` 用來模擬座標系錯配。"""
        out = []
        for y in range(0, size - cell, pitch):
            for x in range(0, size - cell, pitch):
                out.append((y + cell // 2 + offset, x + cell // 2 + offset))
                if len(out) >= n:
                    return out
        return out

    def test_correct_coords_do_not_warn(self, tmp_path, caplog):
        """座標正確時不得警告。"""
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        mask = self._dense_mask()
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, self._bins_on_cells())

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            attr = bin_attribution(mask, tp, 0, 0)

        assert attr.attrs["coverage"]["hit_rate"] > 0.9
        assert "座標對不上" not in caplog.text

    def test_shifted_coords_warn_about_coordinates(self, tmp_path, caplog):
        """bin 全部落在界內卻幾乎不命中 → 必須警告，且指向座標。

        這正是我 2026-07-25 親身踩到的情形：拿 scale=1 去讀一張需要 ×1.54 的
        遮罩，界內率 100%、命中率 3.3%，而當時**沒有任何警告**。
        """
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        mask = self._dense_mask()
        tp = tmp_path / "tp.parquet"
        # 位移到細胞之間的空隙（cell=8、pitch=16 → 偏移 8 剛好全部落空）
        _write_tissue_positions(tp, self._bins_on_cells(offset=8))

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            attr = bin_attribution(mask, tp, 0, 0)

        cov = attr.attrs["coverage"]
        assert cov["n_in_bounds"] == cov["n_total"], "bin 應全部落在界內（越界檢查抓不到）"
        assert cov["hit_rate"] < 0.08
        assert "座標對不上" in caplog.text
        assert "Space Ranger" in caplog.text

    def test_empty_mask_blames_segmentation_not_coordinates(self, tmp_path, caplog):
        """遮罩沒有任何細胞時，訊息須指向分割而非座標。

        兩者症狀相同但要查的地方完全不同 —— 講錯會把人帶去查錯的東西。
        """
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        mask = np.zeros((600, 600), dtype=np.int32)
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, self._bins_on_cells())

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            bin_attribution(mask, tp, 0, 0)

        assert "分割" in caplog.text
        assert "座標對不上" not in caplog.text

    def test_small_crop_does_not_warn_about_out_of_bounds(self, tmp_path, caplog):
        """裁切小窗格時不得對「越界 bin 很多」發警告。

        實跑重現：從全片裁 768×768，全片 4,378,840 個 bin 有 100% 落在窗格外，
        於是跳出「請確認 binned_outputs 是否為對應此影像的註冊版本」——
        資料完全正常，卻把人導去查沒問題的東西。
        """
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        # 遮罩只有 300×300，但 bins 散布在 6000×6000 的全片範圍
        mask = self._dense_mask(size=300)
        rows_cols = [(y, x) for y in range(0, 6000, 40) for x in range(0, 6000, 40)]
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, rows_cols)

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            bin_attribution(mask, tp, 0, 0)

        assert "落在遮罩範圍外" not in caplog.text

    def test_full_slide_still_warns_about_out_of_bounds(self, tmp_path, caplog):
        """整片模式下，越界比例過高仍須警告（原有行為不可退化）。"""
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        # 遮罩涵蓋 bins 的範圍（整片），但 bins 被推到一半在界外
        mask = self._dense_mask(size=600)
        rows_cols = [(y, x) for y in range(0, 1000, 12) for x in range(0, 1000, 12)]
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, rows_cols)

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            bin_attribution(mask, tp, 0, 0)

        assert "落在遮罩範圍外" in caplog.text

    def test_small_bin_count_is_not_judged(self, tmp_path, caplog):
        """bin 太少時不做判斷（統計不可信，小型 ROI 不該被誤觸）。"""
        import logging

        from backend.src.fullslide.pipeline import bin_attribution

        mask = self._dense_mask()
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, self._bins_on_cells(n=50, offset=8))

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide"):
            attr = bin_attribution(mask, tp, 0, 0)

        assert "hit_rate" not in attr.attrs["coverage"]
        assert "座標對不上" not in caplog.text


class TestAggregateCells:
    """bins → cells×genes 聚合"""

    def test_aggregate_cells_sums_bins_per_cell(self, tmp_path, monkeypatch):
        """3 bins → 2 cells：counts 須加總，n_bins 須為 [2, 1]。"""
        import pandas as pd
        import scanpy as sc

        from backend.src.fullslide import pipeline

        h5 = _write_h5(
            tmp_path / "m.h5ad",
            barcodes=["BC0", "BC1", "BC2"],
            gene_names=["GeneA", "GeneB"],
            counts=[[1, 0], [2, 3], [5, 5]],
        )
        # aggregate_cells 內部呼叫 sc.read_10x_h5；測試改讀 h5ad
        monkeypatch.setattr(pipeline_sc(), "read_10x_h5", lambda p: sc.read_h5ad(str(p)))

        attr = pd.DataFrame({"barcode": ["BC0", "BC1", "BC2"], "cell_id": [7, 7, 9]})
        cells = pipeline.aggregate_cells(attr, h5)

        assert cells.n_obs == 2
        assert cells.n_vars == 2
        assert cells.obs["cell_id"].tolist() == [7, 9]
        assert cells.obs["n_bins"].tolist() == [2, 1]
        # cell 7 = BC0 + BC1 = [3, 3]；cell 9 = BC2 = [5, 5]
        assert cells.X.toarray().tolist() == [[3.0, 3.0], [5.0, 5.0]]

    def test_aggregate_cells_ignores_unlisted_barcodes(self, tmp_path, monkeypatch):
        """attribution 未列出的 barcode 不得進入結果。"""
        import pandas as pd
        import scanpy as sc

        from backend.src.fullslide import pipeline

        h5 = _write_h5(
            tmp_path / "m.h5ad",
            barcodes=["BC0", "BC1"],
            gene_names=["GeneA"],
            counts=[[4], [9]],
        )
        monkeypatch.setattr(pipeline_sc(), "read_10x_h5", lambda p: sc.read_h5ad(str(p)))

        attr = pd.DataFrame({"barcode": ["BC0"], "cell_id": [3]})
        cells = pipeline.aggregate_cells(attr, h5)

        assert cells.n_obs == 1
        assert cells.X.toarray().tolist() == [[4.0]]


def pipeline_sc():
    """取得 aggregate_cells 實際使用的 scanpy 模組物件（供 monkeypatch）。"""
    import scanpy as sc
    return sc


# ── add_centroids ───────────────────────────────────────────────────────────

def _cells_with_ids(cell_ids):
    """建立僅含 obs['cell_id'] 的最小 AnnData。"""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    cells = ad.AnnData(X=sp.csr_matrix(np.ones((len(cell_ids), 1), dtype=np.float32)))
    cells.obs["cell_id"] = np.asarray(cell_ids, dtype=int)
    return cells


class TestAddCentroids:
    """細胞重心（裁切局部 px / 全片 fullres px / µm）"""

    def test_add_centroids_sets_local_and_fullres(self):
        """origin_xy 須只影響 fullres 欄位，局部座標維持 CLI 原語意。"""
        from backend.src.fullslide.pipeline import add_centroids

        mask = np.ones((10, 10), dtype=np.int32)   # 單一細胞，重心 = (4.5, 4.5)
        cells = _cells_with_ids([1])

        add_centroids(cells, mask, pixel_size_um=0.2737, origin_xy=(100, 200))

        assert cells.obs["centroid_x_px"].iloc[0] == pytest.approx(4.5)
        assert cells.obs["centroid_y_px"].iloc[0] == pytest.approx(4.5)
        assert cells.obs["centroid_x_fullres"].iloc[0] == pytest.approx(104.5)
        assert cells.obs["centroid_y_fullres"].iloc[0] == pytest.approx(204.5)

    def test_add_centroids_default_origin_is_zero(self):
        """未給 origin_xy 時 fullres 應等於局部座標。"""
        from backend.src.fullslide.pipeline import add_centroids

        mask = np.ones((4, 4), dtype=np.int32)
        cells = _cells_with_ids([1])

        add_centroids(cells, mask, pixel_size_um=0.2737)

        assert cells.obs["centroid_x_fullres"].iloc[0] == pytest.approx(
            cells.obs["centroid_x_px"].iloc[0]
        )

    def test_add_centroids_spatial_is_um_from_local(self):
        """obsm['spatial'] 須由**局部**座標換算（維持既有匯出語意）。"""
        from backend.src.fullslide.pipeline import add_centroids

        mask = np.ones((10, 10), dtype=np.int32)
        cells = _cells_with_ids([1])
        px_um = 0.2737

        add_centroids(cells, mask, pixel_size_um=px_um, origin_xy=(100, 200))

        assert cells.obsm["spatial"][0, 0] == pytest.approx(4.5 * px_um)
        assert cells.obsm["spatial"][0, 1] == pytest.approx(4.5 * px_um)

    def test_add_centroids_multiple_cells_ordered_by_cell_id(self):
        """多細胞時重心須與 obs['cell_id'] 逐列對應。"""
        from backend.src.fullslide.pipeline import add_centroids

        mask = np.zeros((10, 10), dtype=np.int32)
        mask[0:2, 0:2] = 5      # 重心 (0.5, 0.5)
        mask[8:10, 8:10] = 9    # 重心 (8.5, 8.5)
        cells = _cells_with_ids([5, 9])

        add_centroids(cells, mask, pixel_size_um=1.0)

        assert cells.obs["centroid_x_px"].tolist() == pytest.approx([0.5, 8.5])
        assert cells.obs["centroid_y_px"].tolist() == pytest.approx([0.5, 8.5])


# ── run_full 裁切座標驗證 ────────────────────────────────────────────────────

class TestFullSegCropValidation:
    """POST /api/segmentation/run_full 的裁切座標驗證（純驗證邏輯，不觸發分割）"""

    def test_crop_x1_not_greater_than_x0_is_rejected(self):
        from backend.src.api.segmentation import FullSegParams, validate_crop

        err = validate_crop(FullSegParams(crop_x0=100, crop_x1=50))

        assert err is not None
        assert "crop_x1 必須大於 crop_x0" in err

    def test_crop_y1_not_greater_than_y0_is_rejected(self):
        from backend.src.api.segmentation import FullSegParams, validate_crop

        err = validate_crop(FullSegParams(crop_y0=800, crop_y1=800))

        assert err is not None
        assert "crop_y1 必須大於 crop_y0" in err

    def test_negative_origin_is_rejected(self):
        from backend.src.api.segmentation import FullSegParams, validate_crop

        err = validate_crop(FullSegParams(crop_x0=-1, crop_x1=100))

        assert err is not None
        assert "不可為負" in err

    def test_full_image_defaults_are_valid(self):
        """全部為 None（＝全圖）須通過驗證。"""
        from backend.src.api.segmentation import FullSegParams, validate_crop

        assert validate_crop(FullSegParams()) is None

    def test_sentinel_minus_one_is_valid(self):
        """-1 代表影像邊界，須視為合法。"""
        from backend.src.api.segmentation import FullSegParams, validate_crop

        assert validate_crop(FullSegParams(crop_x0=0, crop_x1=-1, crop_y0=0, crop_y1=-1)) is None

    def test_valid_window_passes(self):
        from backend.src.api.segmentation import FullSegParams, validate_crop

        assert validate_crop(
            FullSegParams(crop_x0=100, crop_x1=1124, crop_y0=200, crop_y1=1224)
        ) is None


class TestResolveCropWindow:
    """裁切座標 → 實際切片邊界（-1/None 展開為影像邊界）

    架構深化 P8：`resolve_crop_window` 下沉到 `fullslide.pipeline`，簽章從
    `FullSegParams`（pydantic）改成 4 個 plain scalar，CLI 也能直接呼叫，
    不需要讓領域層反向依賴 API 層的請求型別。
    """

    def test_none_expands_to_full_image(self):
        from backend.src.fullslide.pipeline import resolve_crop_window

        x0, y0, x1, y1 = resolve_crop_window(None, None, None, None, w_img=500, h_img=400)

        assert (x0, y0, x1, y1) == (0, 0, 500, 400)

    def test_minus_one_expands_to_full_image(self):
        from backend.src.fullslide.pipeline import resolve_crop_window

        x0, y0, x1, y1 = resolve_crop_window(10, 20, -1, -1, w_img=500, h_img=400)

        assert (x0, y0, x1, y1) == (10, 20, 500, 400)

    def test_window_is_clamped_to_image_bounds(self):
        """超出影像邊界的請求須被夾住，而非產生越界切片。"""
        from backend.src.fullslide.pipeline import resolve_crop_window

        x0, y0, x1, y1 = resolve_crop_window(0, 0, 9999, 9999, w_img=500, h_img=400)

        assert (x1, y1) == (500, 400)


# ── metadata sidecar ────────────────────────────────────────────────────────

FULL_SEG_META_KEYS = {
    "crop_x0", "crop_y0", "width", "height",
    "n_cells", "pixel_size_um", "passes", "created_at",
}


class TestFullSegMeta:
    """全圖遮罩的 metadata sidecar"""

    def test_write_full_seg_meta_schema(self, tmp_path):
        """sidecar 須含全部 8 個欄位。"""
        import json

        from backend.src.fullslide.pipeline import write_full_seg_meta

        path = write_full_seg_meta(
            tmp_path, crop_x0=100, crop_y0=200, width=512, height=256,
            n_cells=42, pixel_size_um=0.2737, passes=4,
        )
        meta = json.loads(path.read_text())

        assert set(meta) == FULL_SEG_META_KEYS
        assert meta["crop_x0"] == 100
        assert meta["crop_y0"] == 200
        assert meta["n_cells"] == 42
        assert meta["passes"] == 4

    def test_read_full_seg_meta_roundtrip(self, tmp_path):
        from backend.src.fullslide.pipeline import read_full_seg_meta, write_full_seg_meta

        write_full_seg_meta(
            tmp_path, crop_x0=7, crop_y0=9, width=10, height=10,
            n_cells=1, pixel_size_um=0.2737, passes=7,
        )
        meta = read_full_seg_meta(tmp_path)

        assert meta["crop_x0"] == 7
        assert meta["crop_y0"] == 9

    def test_read_full_seg_meta_missing_returns_none(self, tmp_path):
        """sidecar 不存在時回傳 None，由呼叫端決定如何處理。"""
        from backend.src.fullslide.pipeline import read_full_seg_meta

        assert read_full_seg_meta(tmp_path) is None


# ── resolve_pixel_size ──────────────────────────────────────────────────────

class TestResolvePixelSize:
    """樣本實際 µm/px 優先於預設常數"""

    def test_prefers_scalefactors_value(self, tmp_path):
        import json

        from backend.src.fullslide.pipeline import resolve_pixel_size

        spatial = tmp_path / "spatial"
        spatial.mkdir()
        (spatial / "scalefactors_json.json").write_text(
            json.dumps({"microns_per_pixel": 0.4321, "tissue_hires_scalef": 0.1})
        )
        config = {"paths": {"binned_002": str(tmp_path)}}

        assert resolve_pixel_size(config) == pytest.approx(0.4321)

    def test_falls_back_to_constant_when_missing(self, tmp_path):
        from backend.src.fullslide.pipeline import resolve_pixel_size
        from backend.src.utils.constants import VISIUM_UM_PX

        config = {"paths": {"binned_002": str(tmp_path)}}

        assert resolve_pixel_size(config) == pytest.approx(VISIUM_UM_PX)

    def test_falls_back_when_key_absent_in_json(self, tmp_path):
        """json 存在但沒有 microns_per_pixel 時仍須回退，不可拋錯。"""
        import json

        from backend.src.fullslide.pipeline import resolve_pixel_size
        from backend.src.utils.constants import VISIUM_UM_PX

        spatial = tmp_path / "spatial"
        spatial.mkdir()
        (spatial / "scalefactors_json.json").write_text(json.dumps({"tissue_hires_scalef": 0.1}))
        config = {"paths": {"binned_002": str(tmp_path)}}

        assert resolve_pixel_size(config) == pytest.approx(VISIUM_UM_PX)

    def test_falls_back_on_corrupt_json(self, tmp_path):
        """壞掉的 json 不可讓流程中斷。"""
        from backend.src.fullslide.pipeline import resolve_pixel_size
        from backend.src.utils.constants import VISIUM_UM_PX

        spatial = tmp_path / "spatial"
        spatial.mkdir()
        (spatial / "scalefactors_json.json").write_text("{ not json")
        config = {"paths": {"binned_002": str(tmp_path)}}

        assert resolve_pixel_size(config) == pytest.approx(VISIUM_UM_PX)


# ── bin → 遮罩 縮放推導 ─────────────────────────────────────────────────────

def _write_spatial(binned_dir, hires_size, scalef, mpp=0.5464):
    """建立 spatial/ 目錄：scalefactors_json.json + tissue_hires_image.png。"""
    import json

    from PIL import Image

    sp = binned_dir / "spatial"
    sp.mkdir(parents=True, exist_ok=True)
    (sp / "scalefactors_json.json").write_text(json.dumps({
        "microns_per_pixel": mpp,
        "tissue_hires_scalef": scalef,
    }))
    Image.new("RGB", hires_size).save(sp / "tissue_hires_image.png")
    return sp


class TestResolveBinToMaskScale:
    """SR fullres → 遮罩(TIFF) px 的縮放推導

    推導不變量：`tissue_hires_image.png` 尺寸 ÷ `tissue_hires_scalef`
    ＝ Space Ranger fullres 畫布尺寸。此法已對照 EP 的 vfr_CORRECTED
    遮罩 shape (23552, 11266) 驗證一致。
    """

    def test_derives_per_axis_scale_from_hires_and_tiff_shape(self, tmp_path):
        """重現 dpcp01 真實數字：col 1.9088、row 2.0000（非等向）。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_mask_scale

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544)
        config = {"paths": {"binned_002": str(binned)}}

        sx, sy = resolve_bin_to_mask_scale(config, mask_shape=(47104, 21504))

        assert sx == pytest.approx(1.9088, abs=1e-3)
        assert sy == pytest.approx(2.0, abs=1e-3)

    def test_identity_when_shapes_match(self, tmp_path):
        """SR fullres 與遮罩同尺寸時應為 1.0（CRC 這類樣本）。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_mask_scale

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(1000, 2000), scalef=0.1)
        config = {"paths": {"binned_002": str(binned)}}

        sx, sy = resolve_bin_to_mask_scale(config, mask_shape=(20000, 10000))

        assert (sx, sy) == pytest.approx((1.0, 1.0))

    def test_config_override_wins(self, tmp_path):
        """alignment.bin_to_mask_scale 明確指定時優先於自動推導。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_mask_scale

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544)
        config = {
            "paths": {"binned_002": str(binned)},
            "alignment": {"bin_to_mask_scale": [1.5, 1.5]},
        }

        assert resolve_bin_to_mask_scale(config, mask_shape=(47104, 21504)) == pytest.approx((1.5, 1.5))

    def test_falls_back_to_identity_when_spatial_missing(self, tmp_path):
        """缺 spatial 檔案時回退 (1.0, 1.0)，不可拋錯中斷流程。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_mask_scale

        config = {"paths": {"binned_002": str(tmp_path / "nope")}}

        assert resolve_bin_to_mask_scale(config, mask_shape=(100, 100)) == pytest.approx((1.0, 1.0))


class TestResolveBinToImageTransform:
    """對位 JSON 優先、近似縮放回退"""

    @staticmethod
    def _write_align(path, *, scale_transform=0.2, scale_images=0.1, serial="H1-X", area="D1"):
        import json as _json

        path.write_text(_json.dumps({
            "serialNumber": serial,
            "area": area,
            "transform": [scale_transform, 0, 0, 0, scale_transform, 0, 0, 0, 1],
            "cytAssistInfo": {
                "transformImages": [scale_images, 0, 0, 0, scale_images, 0, 0, 0, 1],
            },
        }), encoding="utf-8")
        return path

    def test_auto_composes_from_spatial_jsons(self, tmp_path):
        """spatial/ 內同時有 old（mpp 符合 scalefactors）與 new（≈影像 mpp）時自動組合。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        # 近似縮放推得 ≈ (1.9088, 2.0) → 影像 mpp ≈ 0.5464 / 1.954 ≈ 0.2796
        sp = _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544, mpp=0.5464)
        self._write_align(sp / "old.json", scale_transform=0.2, scale_images=0.2 * 0.5464)
        self._write_align(sp / "new.json", scale_transform=0.2, scale_images=0.2 * 0.2732)
        config = {"paths": {"binned_002": str(binned)}}

        transform, scale, source = resolve_bin_to_image_transform(
            config, mask_shape=(47104, 21504)
        )

        assert transform is not None
        np.testing.assert_allclose(transform, np.diag([2.0, 2.0, 1.0]), rtol=1e-6)
        assert "old.json" in source and "new.json" in source

    def test_explicit_extra_json_is_used_as_h_new(self, tmp_path):
        """alignment.extra_alignment_json 明確指定時直接當作 H_new（不再自動偵測）。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        sp = _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544, mpp=0.5464)
        self._write_align(sp / "old.json", scale_transform=0.2, scale_images=0.2 * 0.5464)
        loupe = self._write_align(
            tmp_path / "loupe.json", scale_transform=0.2, scale_images=0.2 * 0.1366
        )
        config = {
            "paths": {"binned_002": str(binned)},
            "alignment": {"extra_alignment_json": str(loupe)},
        }

        transform, _, source = resolve_bin_to_image_transform(config, mask_shape=(47104, 21504))

        assert transform is not None
        np.testing.assert_allclose(transform, np.diag([4.0, 4.0, 1.0]), rtol=1e-6)
        assert "loupe.json" in source

    def test_falls_back_to_scale_without_alignment_json(self, tmp_path):
        """沒有對位 JSON 時回退近似縮放，且來源說明須標明為近似值。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544)
        config = {"paths": {"binned_002": str(binned)}}

        transform, scale, source = resolve_bin_to_image_transform(
            config, mask_shape=(47104, 21504)
        )

        assert transform is None
        assert scale[0] == pytest.approx(1.9088, abs=1e-3)
        assert "近似" in source

    def test_config_can_disable_alignment_json(self, tmp_path):
        """use_alignment_json=false 時即使有 JSON 也走近似縮放（供對照除錯）。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        sp = _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544, mpp=0.5464)
        self._write_align(sp / "old.json", scale_transform=0.2, scale_images=0.2 * 0.5464)
        self._write_align(sp / "new.json", scale_transform=0.2, scale_images=0.2 * 0.2732)
        config = {
            "paths": {"binned_002": str(binned)},
            "alignment": {"use_alignment_json": False},
        }

        transform, _, source = resolve_bin_to_image_transform(config, mask_shape=(47104, 21504))

        assert transform is None
        assert "停用" in source

    def test_identity_when_segmented_on_source_image(self, tmp_path):
        """只有 H_old（分割就用那張圖）→ 單位矩陣，不做任何縮放。"""
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        sp = _write_spatial(binned, hires_size=(1000, 2000), scalef=0.1, mpp=0.5464)
        self._write_align(sp / "old.json", scale_transform=0.2, scale_images=0.2 * 0.5464)
        config = {"paths": {"binned_002": str(binned)}}

        transform, _, source = resolve_bin_to_image_transform(config, mask_shape=(20000, 10000))

        assert transform is not None
        np.testing.assert_allclose(transform, np.eye(3), atol=1e-12)
        assert "單位矩陣" in source


# ── 全圖計數輸入解析 ─────────────────────────────────────────────────────────

class TestFullCountInputs:
    """/api/count/run_full 的輸入檢查（純函式，不啟動背景任務）"""

    def test_missing_mask_is_reported(self, tmp_path):
        """未跑全圖分割時須明確要求先完成分割。"""
        from backend.src.fullslide.pipeline import resolve_full_count_inputs

        config = {"paths": {"output_dir": str(tmp_path), "binned_002": str(tmp_path)}}

        inputs, err = resolve_full_count_inputs(config)

        assert inputs is None
        assert "請先完成全圖分割" in err

    def test_missing_tissue_positions_is_reported(self, tmp_path):
        """遮罩有了但缺 tissue_positions 時須指出缺哪個檔。"""
        from backend.src.fullslide.pipeline import (
            FULL_SEG_MASK_FILENAME,
            resolve_full_count_inputs,
        )

        out = tmp_path / "out"
        out.mkdir()
        np.save(str(out / FULL_SEG_MASK_FILENAME), np.ones((4, 4), dtype=np.int32))
        config = {"paths": {"output_dir": str(out), "binned_002": str(tmp_path / "nope")}}

        inputs, err = resolve_full_count_inputs(config)

        assert inputs is None
        assert "tissue_positions" in err

    def test_resolves_all_paths_and_origin(self, tmp_path):
        """齊備時回傳遮罩/tp/h5 路徑與 sidecar 的裁切原點。"""
        from backend.src.fullslide.pipeline import (
            FULL_SEG_MASK_FILENAME,
            resolve_full_count_inputs,
            write_full_seg_meta,
        )

        out = tmp_path / "out"
        out.mkdir()
        np.save(str(out / FULL_SEG_MASK_FILENAME), np.ones((4, 4), dtype=np.int32))
        write_full_seg_meta(
            out, crop_x0=11, crop_y0=22, width=4, height=4,
            n_cells=1, pixel_size_um=0.2737, passes=4,
        )
        binned = tmp_path / "b002"
        (binned / "spatial").mkdir(parents=True)
        (binned / "spatial" / "tissue_positions.parquet").write_bytes(b"stub")
        (binned / "filtered_feature_bc_matrix.h5").write_bytes(b"stub")
        config = {"paths": {"output_dir": str(out), "binned_002": str(binned)}}

        inputs, err = resolve_full_count_inputs(config)

        assert err is None
        assert inputs["origin_xy"] == (11, 22)
        assert inputs["scale"] == pytest.approx((1.0, 1.0))   # 無 spatial 檔 → 回退
        assert inputs["dilation_px"] == 0                     # config 未給 → 0
        assert inputs["mask_path"].name == FULL_SEG_MASK_FILENAME
        assert inputs["tp_path"].name == "tissue_positions.parquet"
        assert inputs["h5_path"].name == "filtered_feature_bc_matrix.h5"

    def test_dilation_px_comes_from_rna_counting_config(self, tmp_path):
        """擴張距離須與 ROI 路徑共用 rna_counting.dilation_px，否則兩者不可比。"""
        from backend.src.fullslide.pipeline import (
            FULL_SEG_MASK_FILENAME,
            resolve_full_count_inputs,
        )

        out = tmp_path / "out"
        out.mkdir()
        np.save(str(out / FULL_SEG_MASK_FILENAME), np.ones((4, 4), dtype=np.int32))
        binned = tmp_path / "b002"
        (binned / "spatial").mkdir(parents=True)
        (binned / "spatial" / "tissue_positions.parquet").write_bytes(b"stub")
        (binned / "filtered_feature_bc_matrix.h5").write_bytes(b"stub")
        config = {
            "paths": {"output_dir": str(out), "binned_002": str(binned)},
            "rna_counting": {"dilation_px": 6},
        }

        inputs, err = resolve_full_count_inputs(config)

        assert err is None
        assert inputs["dilation_px"] == 6

    def test_origin_defaults_to_zero_without_sidecar(self, tmp_path):
        """sidecar 缺失時原點視為 (0,0)，並附帶警示旗標而非直接失敗。"""
        from backend.src.fullslide.pipeline import (
            FULL_SEG_MASK_FILENAME,
            resolve_full_count_inputs,
        )

        out = tmp_path / "out"
        out.mkdir()
        np.save(str(out / FULL_SEG_MASK_FILENAME), np.ones((4, 4), dtype=np.int32))
        binned = tmp_path / "b002"
        (binned / "spatial").mkdir(parents=True)
        (binned / "spatial" / "tissue_positions.parquet").write_bytes(b"stub")
        (binned / "filtered_feature_bc_matrix.h5").write_bytes(b"stub")
        config = {"paths": {"output_dir": str(out), "binned_002": str(binned)}}

        inputs, err = resolve_full_count_inputs(config)

        assert err is None
        assert inputs["origin_xy"] == (0, 0)
        assert inputs["meta_missing"] is True


# ── 端點層級 ────────────────────────────────────────────────────────────────

@pytest.mark.asyncio
class TestFullSlideEndpoints:
    """/api/segmentation/run_full 與 /api/count/run_full 的請求層行為"""

    async def _client(self):
        from httpx import ASGITransport, AsyncClient

        from backend.main import app
        return AsyncClient(transport=ASGITransport(app=app), base_url="http://test")

    async def test_run_full_rejects_bad_crop(self):
        """裁切座標非法時須在啟動背景任務前就回錯。"""
        async with await self._client() as c:
            r = await c.post(
                "/api/segmentation/run_full",
                json={"crop_x0": 100, "crop_x1": 50},
            )
        assert r.status_code == 200
        body = r.json()
        assert body["status"] == "error"
        assert "crop_x1 必須大於 crop_x0" in body["message"]

    async def test_run_full_rejects_negative_origin(self):
        async with await self._client() as c:
            r = await c.post("/api/segmentation/run_full", json={"crop_y0": -5})
        assert r.json()["status"] == "error"

    async def test_count_full_status_shape(self):
        """/api/count/full_status 須回傳標準狀態欄位。"""
        async with await self._client() as c:
            r = await c.get("/api/count/full_status")
        body = r.json()
        assert set(body) >= {"status", "progress", "message"}

    async def test_count_run_full_requires_mask(self, monkeypatch, tmp_path):
        """未跑全圖分割時 POST /api/count/run_full 須要求先完成分割。"""
        from backend.src.api import cellpose_count

        monkeypatch.setattr(
            cellpose_count, "load_config",
            lambda: {"paths": {"output_dir": str(tmp_path), "binned_002": str(tmp_path)}},
        )
        async with await self._client() as c:
            r = await c.post("/api/count/run_full")
        body = r.json()
        assert body["status"] == "error"
        assert "請先完成全圖分割" in body["message"]


class TestBinAttributionAlignment:
    """`alignment`（殘餘修正）疊加在主變換之後"""

    def test_bin_attribution_with_translation_shifts_assignment(self, tmp_path):
        from backend.src.fullslide.pipeline import bin_attribution
        from backend.src.registration.align import AffineAlignment

        mask = np.zeros((10, 10), dtype=np.int32)
        mask[8, 3] = 9
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(3, 1)])   # (row=3, col=1)

        # 無修正 → 落在背景（被丟棄）
        assert bin_attribution(mask, tp, 0, 0).empty
        # 平移 (x+2, y+5) → (col=3, row=8) → 命中 label 9
        al = AffineAlignment(matrix=[[1, 0, 2], [0, 1, 5]])
        assert bin_attribution(mask, tp, 0, 0, alignment=al)["cell_id"].tolist() == [9]

    def test_alignment_stacks_on_top_of_scale(self, tmp_path):
        """修正作用在**影像 px 空間**：先縮放、再平移。"""
        from backend.src.fullslide.pipeline import bin_attribution
        from backend.src.registration.align import AffineAlignment

        mask = np.zeros((20, 20), dtype=np.int32)
        mask[11, 7] = 4
        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(3, 2)])   # ×2 → (6, 4)；+ (3, 5) → (11, 7)

        al = AffineAlignment(matrix=[[1, 0, 3], [0, 1, 5]])
        attr = bin_attribution(mask, tp, 0, 0, scale=(2.0, 2.0), alignment=al)

        assert attr["cell_id"].tolist() == [4]

    def test_identity_alignment_is_noop(self, tmp_path):
        """單位矩陣不得改變任何結果（enabled 但未估計時的常態）。"""
        from backend.src.fullslide.pipeline import bin_attribution
        from backend.src.registration.align import AffineAlignment

        tp = tmp_path / "tp.parquet"
        _write_tissue_positions(tp, [(1, 1), (6, 6)])
        mask = _make_mask()

        a = bin_attribution(mask, tp, 0, 0)
        b = bin_attribution(mask, tp, 0, 0, alignment=AffineAlignment.identity())

        assert a["cell_id"].tolist() == b["cell_id"].tolist()


class TestResidualAlignmentConfig:
    """`alignment.enabled` 預設 false → 零回歸"""

    def test_disabled_by_default(self, tmp_path):
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(1000, 2000), scalef=0.1)
        config = {
            "paths": {"binned_002": str(binned)},
            "alignment": {"matrix": [[1, 0, 500], [0, 1, 500]]},   # 未 enabled
        }

        transform, scale, _ = resolve_bin_to_image_transform(config, mask_shape=(20000, 10000))

        assert transform is None
        assert scale == pytest.approx((1.0, 1.0))

    def test_enabled_residual_is_composed(self, tmp_path):
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(1000, 2000), scalef=0.1)
        config = {
            "paths": {"binned_002": str(binned)},
            "alignment": {"enabled": True, "matrix": [[1, 0, 12], [0, 1, -7]]},
        }

        transform, scale, source = resolve_bin_to_image_transform(
            config, mask_shape=(20000, 10000)
        )

        assert transform is not None
        assert transform[0, 2] == pytest.approx(12.0)
        assert transform[1, 2] == pytest.approx(-7.0)
        assert scale == pytest.approx((1.0, 1.0))   # 已併入 transform，不可重複套用
        assert "人工修正" in source


# ── P5：串流 tiled 分割 ──────────────────────────────────────────────────────

def _tiled_cfg():
    return {
        "use_gpu": False, "batch_size": 1,
        "dia_small": 13.0, "dia_mid": 17.0, "dia_large": 22.0,
        "use_hematoxylin": False, "use_cpsam": False,
        "voronoi_distance": 4, "min_size": 20, "max_size": 6000,
        "clahe_clip_limit": 3.0,
    }


def _tissue_image(size=384, seed=0):
    """整片都算「組織」（灰階 < 220）的合成影像。"""
    rng = np.random.default_rng(seed)
    return rng.integers(30, 180, (size, size, 3), dtype=np.uint8)


class TestTiledStreaming:
    """`tile_reader` 串流讀取：結果須與傳整張 ndarray 完全相同"""

    def test_run_tiled_accepts_tile_reader(self, fake_cellpose):
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        img = _tissue_image(384)
        cfg = _tiled_cfg()

        by_array = run_tiled_mcseg_v2(img, cfg, tile_size=128, overlap=32)
        by_reader = run_tiled_mcseg_v2(
            cfg=cfg, tile_size=128, overlap=32,
            tile_reader=lambda x, y, w, h: img[y:y + h, x:x + w],
            full_shape=img.shape[:2],
        )

        assert np.array_equal(by_array, by_reader)
        assert by_array.max() > 0

    def test_requires_exactly_one_source(self, fake_cellpose):
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        img = _tissue_image(128)
        with pytest.raises(ValueError):
            run_tiled_mcseg_v2(cfg=_tiled_cfg())                    # 兩者皆無
        with pytest.raises(ValueError):
            run_tiled_mcseg_v2(img, _tiled_cfg(),
                               tile_reader=lambda *a: img)          # 兩者皆有

    def test_tile_reader_requires_full_shape(self, fake_cellpose):
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        with pytest.raises(ValueError) as exc:
            run_tiled_mcseg_v2(cfg=_tiled_cfg(), tile_reader=lambda *a: None)
        assert "full_shape" in str(exc.value)


class TestTiledMemmapAndResume:
    """memmap 落地與中斷續跑"""

    def test_memmap_result_matches_in_memory(self, fake_cellpose, tmp_path):
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        img = _tissue_image(384)
        cfg = _tiled_cfg()

        in_ram = run_tiled_mcseg_v2(img, cfg, tile_size=128, overlap=32)
        on_disk = run_tiled_mcseg_v2(img, cfg, tile_size=128, overlap=32, work_dir=tmp_path)

        assert np.array_equal(in_ram, on_disk)

    def test_temp_files_removed_only_after_success(self, fake_cellpose, tmp_path):
        """暫存檔在最終遮罩產出後才刪；提早刪會讓中斷的執行無從續跑。"""
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        run_tiled_mcseg_v2(_tissue_image(256), _tiled_cfg(),
                           tile_size=128, overlap=32, work_dir=tmp_path)

        assert list(tmp_path.glob("tmp_labels_*.npy")) == []
        assert list(tmp_path.glob("tmp_expanded_*.npy")) == []
        assert not (tmp_path / "full_seg_progress.json").exists()

    def test_resume_skips_done_tiles(self, fake_cellpose, tmp_path, monkeypatch):
        """續跑時已完成的 tile 不可再讀圖（reader 不應被呼叫於該區域）。"""
        import backend.src.segmentation.cellpose_runner as cr

        img = _tissue_image(256)
        cfg = _tiled_cfg()

        # 第一趟：跑完第一列後「中斷」（以例外模擬）
        boom = RuntimeError("simulated interruption")
        calls: list[tuple[int, int, int, int]] = []
        real_clahe = cr.apply_clahe

        n_tiles = {"count": 0}

        def flaky_clahe(image, **kw):
            n_tiles["count"] += 1
            if n_tiles["count"] > 2:      # 128px tile、256px 影像 → 一列 2 塊
                raise boom
            return real_clahe(image, **kw)

        monkeypatch.setattr(cr, "apply_clahe", flaky_clahe)
        with pytest.raises(RuntimeError):
            cr.run_tiled_mcseg_v2(img, cfg, tile_size=128, overlap=32, work_dir=tmp_path)

        progress = tmp_path / "full_seg_progress.json"
        assert progress.exists(), "中斷前應已寫出第一列的進度"

        # 第二趟：續跑 —— 已完成的 tile 不得再進 reader
        monkeypatch.setattr(cr, "apply_clahe", real_clahe)

        def counting_reader(x, y, w, h):
            calls.append((x, y, w, h))
            return img[y:y + h, x:x + w]

        cr.run_tiled_mcseg_v2(
            cfg=cfg, tile_size=128, overlap=32, work_dir=tmp_path,
            tile_reader=counting_reader, full_shape=img.shape[:2],
        )

        # Phase 2 也會用 reader 重算組織遮罩（整塊 256×256），故以 tile 尺寸區分
        phase1 = [(x, y) for x, y, w, h in calls if w <= 160 and h <= 160]
        assert all(y > 0 for _, y in phase1), f"第一列已完成，不應重新讀取：{phase1}"
        assert len(phase1) == 2, f"應只重跑第二列的 2 塊，實得 {phase1}"

    def test_changed_config_invalidates_resume(self, fake_cellpose, tmp_path):
        """參數變更後舊進度必須作廢 —— 否則會拼出混合兩組參數的錯誤遮罩。"""
        import json

        from backend.src.segmentation.cellpose_runner import (
            _load_seg_progress,
            _save_seg_progress,
            config_hash,
        )

        p = tmp_path / "full_seg_progress.json"
        h1 = config_hash({"dia_mid": 17.0}, 100, 100)
        h2 = config_hash({"dia_mid": 22.0}, 100, 100)
        _save_seg_progress(p, h1, {(0, 0), (0, 1)}, 42)

        assert h1 != h2
        assert _load_seg_progress(p, h1)["current_max"] == 42
        assert _load_seg_progress(p, h2) is None
        assert json.loads(p.read_text())["done_tiles"] == [[0, 0], [0, 1]]

    def test_corrupt_progress_falls_back_to_restart(self, tmp_path):
        from backend.src.segmentation.cellpose_runner import _load_seg_progress

        p = tmp_path / "full_seg_progress.json"
        p.write_text("{ not json", encoding="utf-8")

        assert _load_seg_progress(p, "abc123") is None


class TestBlockedPostprocessing:
    """分塊 Voronoi 與全域面積過濾"""

    def test_blocked_voronoi_matches_whole_image(self):
        """分塊結果須與整圖 Voronoi 完全相同（margin ≥ 擴張距離即為精確）。"""
        from backend.src.segmentation.cellpose_runner import (
            blocked_voronoi,
            voronoi_expand,
        )

        rng = np.random.default_rng(3)
        labels = np.zeros((256, 256), dtype=np.int32)
        for i in range(1, 21):
            cy, cx = rng.integers(10, 246, 2)
            labels[cy - 3:cy + 3, cx - 3:cx + 3] = i

        whole = voronoi_expand(labels, max_distance=6)
        out = np.zeros_like(labels)
        blocked_voronoi(labels, out, max_distance=6, block=64, margin=16)

        assert np.array_equal(whole, out)

    def test_blocked_voronoi_no_seam(self):
        """跨塊邊界的細胞必須維持同一個 label，不可被切成兩個。"""
        from backend.src.segmentation.cellpose_runner import blocked_voronoi

        labels = np.zeros((128, 128), dtype=np.int32)
        labels[60:68, 60:68] = 7           # 正好跨在 block=64 的邊界上

        out = np.zeros_like(labels)
        blocked_voronoi(labels, out, max_distance=5, block=64, margin=16)

        assert set(np.unique(out)) == {0, 7}
        assert out[64, 64] == 7            # 邊界兩側同一 label

    def test_blocked_voronoi_rejects_small_margin(self):
        """margin < 擴張距離時結果不再精確 —— 必須明確報錯而非靜默出錯。"""
        from backend.src.segmentation.cellpose_runner import blocked_voronoi

        labels = np.zeros((32, 32), dtype=np.int32)
        with pytest.raises(ValueError):
            blocked_voronoi(labels, np.zeros_like(labels),
                            max_distance=10, block=16, margin=4)

    def test_global_label_sizes_counts_across_blocks(self):
        """跨塊細胞的面積須加總，否則會被誤判為過小而刪除。"""
        from backend.src.segmentation.cellpose_runner import global_label_sizes

        labels = np.zeros((128, 128), dtype=np.int32)
        labels[60:68, 60:68] = 5           # 跨 block=64 邊界，共 64 px

        sizes = global_label_sizes(labels, block=64)

        assert sizes[5] == 64

    def test_clean_and_relabel_blocked_matches_whole_image(self):
        from backend.src.segmentation.cellpose_runner import (
            clean_and_relabel_blocked,
            clean_mask,
            relabel_sequential,
        )

        rng = np.random.default_rng(4)
        labels = np.zeros((128, 128), dtype=np.int32)
        for i in range(1, 15):
            cy, cx = rng.integers(5, 120, 2)
            r = int(rng.integers(1, 6))
            labels[cy - r:cy + r, cx - r:cx + r] = i

        expected = relabel_sequential(clean_mask(labels.copy(), min_size=20, max_size=60))
        got = labels.copy()
        n = clean_and_relabel_blocked(got, min_size=20, max_size=60, block=64)

        assert np.array_equal(expected, got)
        assert n == int(got.max())


class TestRunFullStreaming:
    """`/api/segmentation/run_full` 走串流路徑（不再整圖進 RAM）"""

    @pytest.mark.asyncio
    async def test_run_full_streams_from_slide_reader(self, fake_cellpose, tmp_path, monkeypatch):
        """端到端：合成 BTF → 串流 tiled 分割 → 產出遮罩與 sidecar。"""
        import json

        import tifffile

        import backend.src.api.segmentation as seg

        rng = np.random.default_rng(0)
        img = rng.integers(30, 180, (384, 384, 3), dtype=np.uint8)
        btf = tmp_path / "slide.btf"
        tifffile.imwrite(str(btf), img, bigtiff=True, tile=(128, 128), photometric="rgb")

        out = tmp_path / "out"
        config = {
            "paths": {"he_image": str(btf), "output_dir": str(out), "binned_002": ""},
            "segmentation": {"mcseg_v2": _tiled_cfg()},
            "full_seg": {"tile_size": 128, "overlap": 32, "max_load_gb": 6.0,
                         "force_disable_cpsam": True},
        }
        monkeypatch.setattr(seg, "resolve_path", lambda p: Path(p))

        await seg._run_full_segmentation(config, seg.FullSegParams())

        assert seg._full_status["status"] == "done", seg._full_status
        mask = np.load(str(out / "full_image_segmentation_masks.npy"))
        assert mask.shape == (384, 384)
        assert mask.max() > 0

        meta = json.loads((out / "full_image_segmentation_meta.json").read_text())
        assert (meta["width"], meta["height"]) == (384, 384)
        # 暫存檔已清掉（最終遮罩產出後才刪）
        assert list(out.glob("tmp_labels_*.npy")) == []

    @pytest.mark.asyncio
    async def test_run_full_rejects_oversized_mask(self, tmp_path, monkeypatch):
        """遮罩本身仍受 RAM 限制（int32 4 bytes/px）→ 超限須明確報錯。"""
        import tifffile

        import backend.src.api.segmentation as seg

        img = np.zeros((256, 256, 3), dtype=np.uint8)
        btf = tmp_path / "slide.btf"
        tifffile.imwrite(str(btf), img, bigtiff=True, tile=(128, 128), photometric="rgb")

        config = {
            "paths": {"he_image": str(btf), "output_dir": str(tmp_path / "out")},
            "segmentation": {"mcseg_v2": _tiled_cfg()},
            "full_seg": {"tile_size": 128, "overlap": 32,
                         "max_load_gb": 1e-6, "force_disable_cpsam": True},
        }
        monkeypatch.setattr(seg, "resolve_path", lambda p: Path(p))

        await seg._run_full_segmentation(config, seg.FullSegParams())

        assert seg._full_status["status"] == "error"


class TestSegmentationApiLayerIsThin:
    """回歸：全片分割編排邏輯只能有一份（`fullslide.pipeline.run_full_slide_segmentation`）"""

    def test_no_orchestration_left_in_run_full_segmentation(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "segmentation.py").read_text(
            encoding="utf-8"
        )
        start = src.index("async def _run_full_segmentation")
        end = src.index("@router", start)
        body = src[start:end]

        for gone in ("mask_gb", "tile_reader", "seg_cfg_safe"):
            assert gone not in body, f"{gone} 不應再出現在 _run_full_segmentation（改呼叫 run_full_slide_segmentation）"

    def test_api_delegates_to_run_full_slide_segmentation(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "segmentation.py").read_text(
            encoding="utf-8"
        )
        assert "run_full_slide_segmentation" in src


class TestMemmapLabelsOnDisk:
    """標籤圖確實落在磁碟（全片 21504×47104 int32 ≈ 4 GB，不可留在 RAM）"""

    def test_labels_are_memmapped_during_run(self, fake_cellpose, tmp_path, monkeypatch):
        import backend.src.segmentation.cellpose_runner as cr

        img = _tissue_image(256)
        seen: dict = {}
        real_clahe = cr.apply_clahe

        def peeking_clahe(image, **kw):
            # 執行中：暫存檔應已存在，且大小等於 H×W×4（int32）
            files = list(tmp_path.glob("tmp_labels_*.npy"))
            if files:
                seen["path"] = files[0]
                seen["size"] = files[0].stat().st_size
            return real_clahe(image, **kw)

        monkeypatch.setattr(cr, "apply_clahe", peeking_clahe)
        cr.run_tiled_mcseg_v2(img, _tiled_cfg(), tile_size=128, overlap=32, work_dir=tmp_path)

        assert "path" in seen, "執行期間應存在 memmap 暫存檔"
        assert seen["size"] >= 256 * 256 * 4      # npy header 之外即為 int32 陣列


class TestCroppedMaskTransform:
    """裁切窗格的遮罩不可被當成整片畫布去推導縮放

    真實案例（dpcp01，2026-07-25 實測）：800×800 的裁切窗格讓
    `resolve_bin_to_mask_scale` 算出 scale (0.071, 0.034)，並連帶讓對位 JSON 的
    `H_new` 自動偵測推估出 mpp 11.13 而全部落空 —— 任何用裁切座標跑的全圖分割，
    計數都會整批錯位。遮罩在影像中的位置由 `origin_xy` 負責，與縮放無關。
    """

    @staticmethod
    def _setup(tmp_path, mask_shape, crop_xy, *, with_image=True, meta_dims=False):
        import json as _json

        import tifffile

        from backend.src.fullslide.pipeline import write_full_seg_meta

        binned = tmp_path / "b002"
        _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544, mpp=0.5464)
        (binned / "filtered_feature_bc_matrix.h5").write_bytes(b"stub")
        _write_tissue_positions(binned / "spatial" / "tissue_positions.parquet", [(1, 1)])

        out = tmp_path / "out"
        out.mkdir(parents=True, exist_ok=True)
        np.save(str(out / "full_image_segmentation_masks.npy"),
                np.zeros(mask_shape, dtype=np.int32))

        he = tmp_path / "slide.btf"
        if with_image:
            tifffile.imwrite(str(he), np.zeros((47104, 21504, 3), dtype=np.uint8),
                             bigtiff=True, tile=(256, 256), photometric="rgb")

        kw = {"image_width": 21504, "image_height": 47104} if meta_dims else {}
        write_full_seg_meta(out, crop_x0=crop_xy[0], crop_y0=crop_xy[1],
                            width=mask_shape[1], height=mask_shape[0],
                            n_cells=0, pixel_size_um=0.2737, passes=4, **kw)

        return {"paths": {"output_dir": str(out), "binned_002": str(binned),
                          "he_image": str(he) if with_image else ""},
                "alignment": {"use_alignment_json": False}}

    def test_cropped_mask_derives_same_scale_as_full_mask(self, tmp_path):
        """800×800 窗格與整片遮罩，推導出的縮放必須相同。"""
        from backend.src.fullslide.pipeline import resolve_full_count_inputs

        # 影像真實尺寸 47104×21504；遮罩只是其中 800×800 的一塊
        cropped = self._setup(tmp_path / "a", (800, 800), (18384, 21712))
        full = self._setup(tmp_path / "b", (47104, 21504), (0, 0))

        inp_c, _ = resolve_full_count_inputs(cropped)
        inp_f, _ = resolve_full_count_inputs(full)

        assert inp_c["scale"] == pytest.approx(inp_f["scale"])
        assert inp_c["scale"][0] == pytest.approx(1.9088, abs=1e-3)
        assert inp_c["origin_xy"] == (18384, 21712)

    def test_sidecar_image_dims_take_priority(self, tmp_path):
        """sidecar 記了影像尺寸就直接用，不必再開影像檔。"""
        from backend.src.fullslide.pipeline import resolve_full_count_inputs

        cfg = self._setup(tmp_path, (800, 800), (18384, 21712),
                          with_image=False, meta_dims=True)
        inp, err = resolve_full_count_inputs(cfg)

        assert err is None
        assert inp["scale"][0] == pytest.approx(1.9088, abs=1e-3)

    def test_falls_back_to_origin_plus_mask_without_image(self, tmp_path):
        """影像與 sidecar 尺寸都沒有時退為「原點＋遮罩」，並且不得拋錯。"""
        from backend.src.fullslide.pipeline import _resolve_image_shape

        shape = _resolve_image_shape(
            {"paths": {"he_image": ""}},
            {"crop_x0": 18384, "crop_y0": 21712},
            (800, 800),
        )
        assert shape == (22512, 19184)
