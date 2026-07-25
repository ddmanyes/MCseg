"""全片分割覆蓋率 QC（`backend/src/fullslide/coverage_qc.py`）

來自 EP 紀錄 `99f6fad4` 的真實失敗：全片分割在 Day0/Day3 區域只分出 380 顆
邊緣碎片，而原始 8µm 該區有 31,969 bins（原始掃描涵蓋完整）。這件事在全片
總命中率上看不出來 —— 一個區域整片漏掉只讓總數低幾個百分點。

本檔用合成資料重現兩種失敗模式，並釘住「不該誤報」的情形。
"""
import numpy as np
import pandas as pd
import pytest

GRID = 64          # 測試用小網格，讓合成遮罩維持在合理尺寸
CELL = 8           # 合成細胞邊長


def _write_tp(path, rows_cols):
    pd.DataFrame({
        "barcode": [f"BC{i}" for i in range(len(rows_cols))],
        "in_tissue": [1] * len(rows_cols),
        "pxl_row_in_fullres": [rc[0] for rc in rows_cols],
        "pxl_col_in_fullres": [rc[1] for rc in rows_cols],
    }).to_parquet(str(path), index=False)
    return path


def _bins_in_grid(gy, gx, n, grid=GRID):
    """在第 (gy, gx) 格內產生 n 個均勻散佈的 bin 座標。"""
    side = int(np.ceil(np.sqrt(n)))
    step = max(grid // (side + 1), 1)
    out = []
    for i in range(side):
        for j in range(side):
            if len(out) >= n:
                return out
            out.append((gy * grid + (i + 1) * step, gx * grid + (j + 1) * step))
    return out


def _fill_cells(mask, gy, gx, *, cell=CELL, next_label=1, grid=GRID):
    """在第 (gy, gx) 格內鋪滿 `cell`×`cell` 的方形細胞，回傳下一個可用 label。"""
    for y in range(gy * grid, (gy + 1) * grid - cell + 1, cell * 2):
        for x in range(gx * grid, (gx + 1) * grid - cell + 1, cell * 2):
            mask[y:y + cell, x:x + cell] = next_label
            next_label += 1
    return next_label


class TestCoverageQcDetectsFailures:
    """兩種失敗模式都要抓到"""

    def test_detects_region_with_bins_but_no_cells(self, tmp_path):
        """整片漏掉：格 (0,1) 有一樣多的 bins 卻完全沒有細胞 → 必須標記。

        這是 Day0/Day3 的簡化版。`in_tissue` 是 Space Ranger 的判定，與我們的
        分割獨立 —— 所以「bins 很多、細胞為 0」是分割失敗的證據。
        """
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID * 2), dtype=np.int32)
        _fill_cells(mask, 0, 0)          # 只有左格有細胞，右格全空

        tp = _write_tp(tmp_path / "tp.parquet",
                       _bins_in_grid(0, 0, 300) + _bins_in_grid(0, 1, 300))

        df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=100)

        flagged = df[df["flagged"]]
        assert len(flagged) == 1, f"預期只標記 1 格，實際 {len(flagged)}"
        assert (flagged.iloc[0]["grid_y"], flagged.iloc[0]["grid_x"]) == (0, 1)
        assert flagged.iloc[0]["n_cells"] == 0
        assert flagged.iloc[0]["n_bins"] >= 100
        assert "密度過低" in flagged.iloc[0]["flag_reason"]

    def test_detects_fragmented_region(self, tmp_path):
        """只剩碎片：細胞數不少，但面積遠小於其他格 → 必須標記。

        Day0/Day3 的 380 顆「邊緣碎片」正是這個樣子 —— 只看細胞數不會覺得異常。
        """
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID * 2), dtype=np.int32)
        nxt = _fill_cells(mask, 0, 0, cell=8)
        # 右格：一樣多的細胞，但每顆只有 1×1 px（碎片）
        for y in range(GRID * 0 + 2, GRID - 2, 16):
            for x in range(GRID + 2, GRID * 2 - 2, 16):
                mask[y, x] = nxt
                nxt += 1

        tp = _write_tp(tmp_path / "tp.parquet",
                       _bins_in_grid(0, 0, 300) + _bins_in_grid(0, 1, 300))

        df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=100)

        right = df[(df["grid_y"] == 0) & (df["grid_x"] == 1)].iloc[0]
        assert right["n_cells"] > 0, "此格確實有細胞，不是整片漏掉"
        assert right["flagged"], "細胞碎片化未被標記"
        assert "碎片" in right["flag_reason"]

    def test_uniform_slide_flags_nothing(self, tmp_path):
        """密度一致時不得誤報 —— 會亂叫的 QC 沒人會理它。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID * 2, GRID * 2), dtype=np.int32)
        nxt = 1
        rows_cols = []
        for gy in range(2):
            for gx in range(2):
                nxt = _fill_cells(mask, gy, gx, next_label=nxt)
                rows_cols += _bins_in_grid(gy, gx, 300)

        df = compute_coverage_qc(mask, _write_tp(tmp_path / "tp.parquet", rows_cols),
                                 grid_px=GRID, min_bins=100)

        assert not df["flagged"].any(), (
            f"密度一致卻誤報：\n{df[df['flagged']][['grid_y','grid_x','flag_reason']]}"
        )


class TestCoverageQcStatistics:
    """統計量本身的正確性"""

    def test_sparse_edge_grids_excluded_from_median(self, tmp_path):
        """bins 太少的格子不納入中位數，也不標記。

        否則切片邊緣的稀疏格會拉低基準，讓真正失敗的區域看起來正常。
        """
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID * 2), dtype=np.int32)
        _fill_cells(mask, 0, 0)

        # 右格只有 5 個 bins（邊緣孤島），且沒有細胞
        tp = _write_tp(tmp_path / "tp.parquet",
                       _bins_in_grid(0, 0, 300) + _bins_in_grid(0, 1, 5))

        df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=100)

        right = df[(df["grid_y"] == 0) & (df["grid_x"] == 1)].iloc[0]
        assert not right["considered"], "bins 不足的格子不該納入評估"
        assert not right["flagged"], "bins 不足的格子不該被標記"
        assert df.attrs["summary"]["n_grid_considered"] == 1

    def test_bin_counts_match_grid_assignment(self, tmp_path):
        """每格 bin 數必須等於實際落在該格的 bin 數。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID * 2), dtype=np.int32)
        mask[0:CELL, 0:CELL] = 1

        left, right = _bins_in_grid(0, 0, 120), _bins_in_grid(0, 1, 47)
        tp = _write_tp(tmp_path / "tp.parquet", left + right)

        df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=1)

        got = {(r.grid_y, r.grid_x): r.n_bins for r in df.itertuples()}
        assert got[(0, 0)] == len(left)
        assert got[(0, 1)] == len(right)

    def test_out_of_bounds_bins_are_not_counted(self, tmp_path):
        """落在遮罩外的 bins 不得計入任何格 —— 否則會虛報一個高 bin 密度的格子。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID), dtype=np.int32)
        mask[0:CELL, 0:CELL] = 1

        inside = _bins_in_grid(0, 0, 50)
        outside = [(GRID + 500, GRID + 500)] * 30      # 遠在遮罩之外
        tp = _write_tp(tmp_path / "tp.parquet", inside + outside)

        df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=1)

        assert int(df["n_bins"].sum()) == len(inside)

    def test_crop_origin_is_applied(self, tmp_path):
        """裁切原點須先扣除，且回報的 x0/y0 要是**來源影像**座標。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID), dtype=np.int32)
        _fill_cells(mask, 0, 0)

        # bins 位於來源影像的 (1000+, 2000+)，遮罩原點是 (2000, 1000) = (x, y)
        rows_cols = [(1000 + r, 2000 + c) for r, c in _bins_in_grid(0, 0, 200)]
        tp = _write_tp(tmp_path / "tp.parquet", rows_cols)

        df = compute_coverage_qc(mask, tp, crop_y0=1000, crop_x0=2000,
                                 grid_px=GRID, min_bins=100)

        assert int(df["n_bins"].sum()) == 200, "扣除原點後 bins 應全部落在格內"
        row = df.iloc[0]
        assert (row["x0"], row["y0"]) == (2000, 1000)

    def test_memmap_mask_is_not_fully_loaded(self, tmp_path):
        """遮罩為 memmap 時要能直接處理（全片可達數 GB）。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        path = tmp_path / "mask.npy"
        arr = np.zeros((GRID * 2, GRID * 2), dtype=np.int32)
        _fill_cells(arr, 0, 0)
        np.save(str(path), arr)

        mm = np.load(str(path), mmap_mode="r")
        assert isinstance(mm, np.memmap)

        tp = _write_tp(tmp_path / "tp.parquet", _bins_in_grid(0, 0, 300))
        df = compute_coverage_qc(mm, tp, grid_px=GRID, min_bins=100)

        assert df[(df["grid_y"] == 0) & (df["grid_x"] == 0)].iloc[0]["n_cells"] > 0


    def test_globally_low_density_warns_even_when_nothing_flagged(self, tmp_path, caplog):
        """全片都失敗時，逐格標記會全數通過 —— 必須另外警告。

        這是中位數相對門檻的結構性盲點：沒有格子會「相對」偏低，因為中位數
        本身就是失敗值。若不處理，QC 會在最該叫的時候回報「全部通過」。
        """
        import logging

        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        # 四格都只有 1 顆細胞、各 300 bins → 密度一致但極低（3.3 cells/1k bins）
        mask = np.zeros((GRID * 2, GRID * 2), dtype=np.int32)
        rows_cols = []
        label = 1
        for gy in range(2):
            for gx in range(2):
                mask[gy * GRID:gy * GRID + CELL, gx * GRID:gx * GRID + CELL] = label
                label += 1
                rows_cols += _bins_in_grid(gy, gx, 300)

        tp = _write_tp(tmp_path / "tp.parquet", rows_cols)

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide.coverage"):
            df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=100)

        assert not df["flagged"].any(), "密度一致，逐格標記本來就不會叫"
        assert df.attrs["summary"]["global_density_ok"] is False
        assert "全片" in caplog.text

    def test_healthy_density_does_not_trigger_global_warning(self, tmp_path):
        """密度正常時不得誤觸全片警告。"""
        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID * 2, GRID * 2), dtype=np.int32)
        nxt, rows_cols = 1, []
        for gy in range(2):
            for gx in range(2):
                nxt = _fill_cells(mask, gy, gx, next_label=nxt)
                rows_cols += _bins_in_grid(gy, gx, 300)

        df = compute_coverage_qc(mask, _write_tp(tmp_path / "tp.parquet", rows_cols),
                                 grid_px=GRID, min_bins=100)

        assert df.attrs["summary"]["global_density_ok"] is True


class TestCoverageQcCoordinateSystem:
    """QC 必須與實際計數走同一條座標解析路徑"""

    def test_shares_transform_helper_with_bin_attribution(self, tmp_path):
        """QC 與 `bin_attribution` 對同一組 bins 必須算出相同的遮罩座標。

        兩者若各自實作座標換算，遲早會漂移成不同座標系 —— 那時 QC 會開始
        對「其實正常」的區域亂叫，或漏掉真正的失敗。
        """
        from backend.src.fullslide.pipeline import map_bins_to_mask, read_tissue_bins

        rows_cols = _bins_in_grid(0, 0, 40)
        tp_path = _write_tp(tmp_path / "tp.parquet", rows_cols)
        transform = np.array([[2.0, 0, -5.0], [0, 2.0, 3.0], [0, 0, 1.0]])

        tp = read_tissue_bins(tp_path)
        row, col, in_bounds, desc = map_bins_to_mask(
            tp, (GRID * 4, GRID * 4), 7, 11, transform=transform
        )

        # 手算對照：先套 homography，再扣原點
        exp_x = 2.0 * np.array([c for _, c in rows_cols]) - 5.0
        exp_y = 2.0 * np.array([r for r, _ in rows_cols]) + 3.0
        np.testing.assert_array_equal(col, np.rint(exp_x - 11).astype(np.int64))
        np.testing.assert_array_equal(row, np.rint(exp_y - 7).astype(np.int64))
        assert desc == "homography"

    def test_zero_considered_grids_warns_about_alignment(self, tmp_path, caplog):
        """完全沒有格子達到 min_bins 時，訊息須指向座標系錯配。

        這是對位設定錯誤最常見的表現（右緣窗格整批落空），使用者需要被告知
        去查對位而不是以為分割壞了。
        """
        import logging

        from backend.src.fullslide.coverage_qc import compute_coverage_qc

        mask = np.zeros((GRID, GRID), dtype=np.int32)
        tp = _write_tp(tmp_path / "tp.parquet", [(GRID + 900, GRID + 900)] * 50)

        with caplog.at_level(logging.WARNING, logger="pipeline.fullslide.coverage"):
            df = compute_coverage_qc(mask, tp, grid_px=GRID, min_bins=100)

        assert df.attrs["summary"]["n_grid_considered"] == 0
        assert not df["flagged"].any()
        assert "座標系錯配" in caplog.text


class TestCoverageQcFromConfig:
    """config 入口：與 Stage 2 共用 resolve_full_count_inputs"""

    def test_writes_parquet_and_records_transform_source(self, tmp_path, monkeypatch):
        from backend.src.fullslide import coverage_qc

        out = tmp_path / "out"
        out.mkdir()
        mask = np.zeros((GRID, GRID), dtype=np.int32)
        _fill_cells(mask, 0, 0)
        np.save(str(out / "mask.npy"), mask)
        tp = _write_tp(tmp_path / "tp.parquet", _bins_in_grid(0, 0, 300))

        def fake_inputs(config):
            return {
                "mask_path": out / "mask.npy",
                "tp_path": tp,
                "h5_path": tmp_path / "x.h5",
                "origin_xy": (0, 0),
                "transform": None,
                "scale": (1.0, 1.0),
                "transform_source": "測試用近似縮放",
                "dilation_px": 0,
                "pixel_size_um": 0.2737,
                "meta_missing": False,
                "output_dir": out,
            }, None

        monkeypatch.setattr(
            "backend.src.fullslide.pipeline.resolve_full_count_inputs", fake_inputs
        )

        df = coverage_qc.run_coverage_qc_from_config(
            {}, grid_px=GRID, min_bins=100
        )

        assert (out / "coverage_qc.parquet").exists()
        assert df.attrs["summary"]["transform_source"] == "測試用近似縮放"
        assert pd.read_parquet(str(out / "coverage_qc.parquet")).shape[0] == len(df)

    @pytest.mark.asyncio
    async def test_api_returns_only_flagged_grids(self, tmp_path, monkeypatch):
        """API 只回傳被標記的格子與摘要 —— 全片可有數百格。"""
        from backend.src.api import cellpose_count

        out = tmp_path / "out"
        out.mkdir()
        mask = np.zeros((GRID, GRID * 2), dtype=np.int32)
        _fill_cells(mask, 0, 0)           # 右格全空 → 應被標記
        np.save(str(out / "mask.npy"), mask)
        tp = _write_tp(tmp_path / "tp.parquet",
                       _bins_in_grid(0, 0, 300) + _bins_in_grid(0, 1, 300))

        monkeypatch.setattr(
            "backend.src.fullslide.pipeline.resolve_full_count_inputs",
            lambda config: ({
                "mask_path": out / "mask.npy", "tp_path": tp,
                "h5_path": tmp_path / "x.h5", "origin_xy": (0, 0),
                "transform": None, "scale": (1.0, 1.0),
                "transform_source": "近似縮放", "dilation_px": 0,
                "pixel_size_um": 0.2737, "meta_missing": False, "output_dir": out,
            }, None),
        )
        monkeypatch.setattr(cellpose_count, "load_config", lambda: {})

        resp = await cellpose_count.run_coverage_qc(grid_px=GRID, min_bins=100)

        assert resp["status"] == "ok"
        assert resp["data"]["summary"]["n_grid_flagged"] == 1
        assert len(resp["data"]["flagged"]) == 1
        assert resp["data"]["flagged"][0]["grid_x"] == 1

    @pytest.mark.asyncio
    async def test_api_reports_missing_mask_as_error(self, monkeypatch):
        """前置條件不足時回 error 而非拋例外。"""
        from backend.src.api import cellpose_count

        monkeypatch.setattr(
            "backend.src.fullslide.pipeline.resolve_full_count_inputs",
            lambda config: (None, "找不到全圖分割遮罩，請先完成全圖分割（Stage 1）"),
        )
        monkeypatch.setattr(cellpose_count, "load_config", lambda: {})

        resp = await cellpose_count.run_coverage_qc()

        assert resp["status"] == "error"
        assert "找不到全圖分割遮罩" in resp["message"]

    def test_propagates_resolve_error(self, tmp_path, monkeypatch):
        """缺遮罩等前置條件不足時，須回報 resolve 的原始訊息。"""
        from backend.src.fullslide import coverage_qc

        monkeypatch.setattr(
            "backend.src.fullslide.pipeline.resolve_full_count_inputs",
            lambda config: (None, "找不到全圖分割遮罩，請先完成全圖分割（Stage 1）"),
        )

        with pytest.raises(ValueError, match="找不到全圖分割遮罩"):
            coverage_qc.run_coverage_qc_from_config({})
