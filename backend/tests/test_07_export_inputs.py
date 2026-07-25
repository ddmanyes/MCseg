"""Test 7: 匯出層的幾何轉換與輸入解析（不需真實資料）"""
import numpy as np
import pytest


# ─── geometry ───────────────────────────────────────────────

@pytest.fixture
def toy_mask(tmp_path):
    """40x40 label mask：一顆 10x10 細胞（label 1）、一顆 2x2 小雜訊（label 2）"""
    mask = np.zeros((40, 40), dtype=np.int32)
    mask[5:15, 5:15] = 1
    mask[30:32, 30:32] = 2
    path = tmp_path / "segmentation_masks.npy"
    np.save(str(path), mask)
    return path


class TestMaskToGeoJSON:
    """mask → GeoJSON 的座標與過濾規則"""

    def test_filters_small_cells(self, toy_mask):
        """面積小於 min_area_px 的細胞被剔除（2x2=4 px < 20）"""
        from backend.src.export.geometry import mask_to_geojson

        geo = mask_to_geojson(toy_mask, pixel_size_um=0.5)
        assert geo["type"] == "FeatureCollection"
        assert len(geo["features"]) == 1
        assert geo["features"][0]["properties"]["cell_id"] == 1

    def test_min_area_override_keeps_both(self, toy_mask):
        """放寬門檻後兩顆都留下"""
        from backend.src.export.geometry import mask_to_geojson

        geo = mask_to_geojson(toy_mask, pixel_size_um=0.5, min_area_px=1)
        assert len(geo["features"]) == 2

    def test_coords_in_um_and_closed(self, toy_mask):
        """座標換算成 µm（原點 = ROI 左上角），且多邊形首尾閉合"""
        from backend.src.export.geometry import mask_to_geojson

        px = 0.5
        ring = mask_to_geojson(toy_mask, pixel_size_um=px)["features"][0]["geometry"]["coordinates"][0]
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        # 細胞佔 col/row 5..14 → µm 範圍約 2.5 ~ 7.0
        assert min(xs) == pytest.approx(5 * px, abs=px)
        assert max(xs) == pytest.approx(14 * px, abs=px)
        assert min(ys) == pytest.approx(5 * px, abs=px)
        assert max(ys) == pytest.approx(14 * px, abs=px)
        assert ring[0] == pytest.approx(ring[-1])


class TestShiftGeoJSONCoords:
    """全域偏移（合併模式的 Loupe 匯出用）"""

    def test_shifts_nested_ring_in_place(self):
        from backend.src.export.geometry import shift_geojson_coords

        feat = {"geometry": {"type": "Polygon", "coordinates": [[[0.0, 0.0], [1.0, 2.0]]]}}
        shift_geojson_coords(feat, 10.0, 100.0)
        assert feat["geometry"]["coordinates"] == [[[10.0, 100.0], [11.0, 102.0]]]

    def test_missing_geometry_is_noop(self):
        from backend.src.export.geometry import shift_geojson_coords

        feat = {"properties": {"full_id": "1"}}
        shift_geojson_coords(feat, 1.0, 1.0)   # 不應拋錯
        assert "geometry" not in feat


# ─── inputs ─────────────────────────────────────────────────

def _write_h5ad(path, obs_names, active_roi=None):
    """寫一個最小 h5ad 供 resolver 的 backed 讀取使用"""
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(
        X=np.ones((len(obs_names), 2), dtype=np.float32),
        obs=pd.DataFrame(index=list(obs_names)),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    if active_roi is not None:
        adata.uns["active_roi"] = active_roi
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(path))
    return path


@pytest.fixture
def export_config(tmp_path):
    """兩個 ROI 的假 config，paths 指向 tmp_path"""
    return {
        "paths": {
            "output_dir": str(tmp_path / "analysis"),
            "data_root": str(tmp_path / "data"),
        },
        "rois": [
            {"name": "roi_a", "x": 100, "y": 200, "pixel_size_um": 0.5},
            {"name": "roi_b", "x": 300, "y": 400},   # 無 pixel_size_um → 用預設
        ],
    }


class TestResolveExportInputs:
    """h5ad 尋找、白名單驗證、模式判定、ROI 路徑組裝"""

    def test_candidate_priority(self, tmp_path, export_config):
        """依呼叫端給的候選順序取第一個存在的"""
        from backend.src.export.inputs import resolve_export_inputs

        out = tmp_path / "analysis"
        _write_h5ad(out / "qc_preprocessed.h5ad", ["cell_1"])
        _write_h5ad(out / "umap_computed.h5ad", ["cell_1"])

        inputs = resolve_export_inputs(
            export_config, "", ("umap_computed.h5ad", "qc_preprocessed.h5ad")
        )
        assert inputs.h5ad_path.name == "umap_computed.h5ad"

    def test_roi_subdir_fallback_keeps_active_roi(self, tmp_path, export_config):
        """根目錄找不到時往 roi/{name}/ 找；回歸測試：解析過程不得污染 ROI 清單

        （原 api/export.py:394 的 fallback 迴圈會覆蓋外層 roi_name，
        使單 ROI 模式後續可能拿到錯誤 ROI 的 mask）
        """
        from backend.src.export.inputs import resolve_export_inputs

        out = tmp_path / "analysis"
        _write_h5ad(out / "roi" / "roi_b" / "qc_preprocessed.h5ad", ["cell_1"], active_roi="roi_b")

        inputs = resolve_export_inputs(export_config, "", ("clustered_final.h5ad",))
        assert inputs.h5ad_path.parent.name == "roi_b"
        assert inputs.active_roi == "roi_b"
        assert [r.name for r in inputs.rois] == ["roi_a", "roi_b"]
        # 單 ROI 分支會依 active_roi 取 ROI，必須拿到 roi_b 自己的 mask
        assert inputs.find_roi(inputs.active_roi).mask_path == out / "roi" / "roi_b" / "segmentation_masks.npy"

    def test_merged_mode_detected_by_obs_names(self, tmp_path, export_config):
        """obs_names 含 '__' → 合併模式"""
        from backend.src.export.inputs import resolve_export_inputs

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["roi_a__cell_1", "roi_b__cell_2"])

        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))
        assert inputs.is_merged is True

    def test_single_mode_detected(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["cell_1", "cell_2"], active_roi="roi_a")

        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))
        assert inputs.is_merged is False
        assert inputs.active_roi == "roi_a"

    def test_roi_inputs_paths_and_pixel_size(self, tmp_path, export_config):
        """ROI 路徑組裝集中一處；未指定 pixel_size_um 時回退 VISIUM_UM_PX"""
        from backend.src.export.inputs import resolve_export_inputs
        from backend.src.utils.constants import VISIUM_UM_PX

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["cell_1"])

        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))
        roi_a, roi_b = inputs.rois
        assert roi_a.out_dir == out / "roi" / "roi_a"
        assert roi_a.mask_path == out / "roi" / "roi_a" / "segmentation_masks.npy"
        assert roi_a.pixel_size_um == 0.5
        assert roi_b.pixel_size_um == pytest.approx(VISIUM_UM_PX)
        assert roi_a.cfg["x"] == 100   # 原始 config 條目供 transcripts 用

    def test_absolute_path_outside_allowed_dirs_rejected(self, tmp_path, export_config):
        """絕對路徑必須位於 output_dir 或 data_root 底下"""
        from backend.src.export.inputs import resolve_export_inputs

        outsider = tmp_path / "elsewhere" / "evil.h5ad"
        _write_h5ad(outsider, ["cell_1"])

        with pytest.raises(ValueError, match="output_dir 或 data_root"):
            resolve_export_inputs(export_config, str(outsider), ("umap_computed.h5ad",))

    def test_absolute_path_inside_data_root_accepted(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs

        inside = tmp_path / "data" / "mine.h5ad"
        _write_h5ad(inside, ["cell_1"])

        inputs = resolve_export_inputs(export_config, str(inside), ("umap_computed.h5ad",))
        assert inputs.h5ad_path == inside.resolve()

    def test_missing_relative_path_raises(self, export_config):
        from backend.src.export.inputs import resolve_export_inputs

        with pytest.raises(FileNotFoundError, match="找不到指定的 h5ad"):
            resolve_export_inputs(export_config, "nope.h5ad", ("umap_computed.h5ad",))

    def test_nothing_found_raises(self, export_config):
        from backend.src.export.inputs import resolve_export_inputs

        with pytest.raises(FileNotFoundError, match="請先執行 Stage 3"):
            resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))


class TestApiExportIsThinnerNow:
    """回歸：匯出的核心轉換不得再住在 API 層"""

    def test_pure_helpers_moved_out(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "export.py").read_text(encoding="utf-8")
        for gone in ("def _mask_to_geojson", "def _generate_visiumhd_transcripts", "def _shift_geojson_coords"):
            assert gone not in src, f"{gone} 應已搬進 backend/src/export/"
        assert "resolve_export_inputs" in src

    def test_h5ad_lookup_not_duplicated(self):
        """白名單驗證只該出現在 resolver 一處"""
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "export.py").read_text(encoding="utf-8")
        assert "必須位於 output_dir 或 data_root 底下" not in src
