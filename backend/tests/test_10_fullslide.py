"""Test 10: 全片流程共用模組（backend/src/fullslide/pipeline.py）

全部使用合成資料 —— 不需真實 CRC/BTF 資料，Windows 亦可執行。
"""
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
