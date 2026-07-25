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
