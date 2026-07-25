"""Test 8: 匯出編排 export/jobs.py（以假 exporter 取代重量級寫檔，不需真實資料）"""
import json

import numpy as np
import pytest


def _write_h5ad(path, obs_names, active_roi=None, roi_col=None):
    import anndata as ad
    import pandas as pd

    obs = pd.DataFrame(index=list(obs_names))
    if roi_col is not None:
        obs["roi"] = roi_col
    adata = ad.AnnData(
        X=np.ones((len(obs_names), 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    if active_roi is not None:
        adata.uns["active_roi"] = active_roi
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(path))
    return path


def _write_mask(path, label_boxes):
    """label_boxes: {label: (r0, r1, c0, c1)}"""
    mask = np.zeros((40, 40), dtype=np.int32)
    for label, (r0, r1, c0, c1) in label_boxes.items():
        mask[r0:r1, c0:c1] = label
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(str(path), mask)
    return path


class FakeExporter:
    """記錄建構參數與 export() 呼叫，取代 Xenium/Loupe 的真實寫檔。"""
    instances: list = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.calls: list = []
        FakeExporter.instances.append(self)

    def export(self, h5ad_path, out_dir):
        self.calls.append((h5ad_path, out_dir))
        return out_dir


@pytest.fixture(autouse=True)
def _reset_fake():
    FakeExporter.instances = []
    yield
    FakeExporter.instances = []


@pytest.fixture
def export_config(tmp_path):
    return {
        "paths": {
            "output_dir": str(tmp_path / "analysis"),
            "data_root": str(tmp_path / "data"),
            "export_dir": str(tmp_path / "export"),
        },
        "rois": [
            {"name": "roi_a", "x": 100, "y": 200, "pixel_size_um": 0.5},
            {"name": "roi_b", "x": 300, "y": 400, "pixel_size_um": 0.5},
        ],
    }


class TestXeniumSingleRoi:
    """單 ROI 模式：進度回報、產出路徑、傳給 exporter 的參數"""

    def _inputs(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["cell_1", "cell_2"], active_roi="roi_a")
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})
        return resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))

    def test_returns_result_and_reports_progress(self, tmp_path, export_config):
        from backend.src.export.jobs import _run_xenium_single

        seen: list = []
        result = _run_xenium_single(
            self._inputs(tmp_path, export_config), "", lambda f, m: seen.append((f, m)), FakeExporter
        )

        assert result.is_merged is False
        assert result.output_dirs == [tmp_path / "analysis" / "roi" / "roi_a" / "export_xenium"]
        assert [f for f, _ in seen] == sorted(f for f, _ in seen)      # 進度單調遞增
        assert all(0.0 <= f <= 1.0 for f, _ in seen)
        assert seen, "應至少回報一次進度"

    def test_polygons_written_and_passed_to_exporter(self, tmp_path, export_config):
        from backend.src.export.jobs import _run_xenium_single

        inputs = self._inputs(tmp_path, export_config)
        _run_xenium_single(inputs, "", lambda f, m: None, FakeExporter)

        poly = tmp_path / "analysis" / "roi" / "roi_a" / "cellpose_polygons.json"
        assert poly.exists()
        assert len(json.loads(poly.read_text(encoding="utf-8"))["features"]) == 1

        exporter = FakeExporter.instances[-1]
        assert exporter.kwargs["poly_json_path"] == poly
        assert exporter.kwargs["pixel_size_um"] == 0.5
        assert exporter.kwargs["he_image_path"] is None       # 無 he_crop.tif
        assert exporter.calls[0][0] == inputs.h5ad_path

    def test_output_dir_override(self, tmp_path, export_config):
        from backend.src.export.jobs import _run_xenium_single

        custom = tmp_path / "custom_out"
        result = _run_xenium_single(
            self._inputs(tmp_path, export_config), str(custom), lambda f, m: None, FakeExporter
        )
        assert result.output_dirs == [custom]

    def test_missing_mask_raises(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs
        from backend.src.export.jobs import _run_xenium_single

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["cell_1"], active_roi="roi_a")
        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))

        with pytest.raises(FileNotFoundError, match="segmentation_masks.npy"):
            _run_xenium_single(inputs, "", lambda f, m: None, FakeExporter)


class TestSubsetRoi:
    """合併 h5ad 取單一 ROI 子集並還原 obs_names"""

    def test_strips_roi_prefix(self):
        from backend.src.export.jobs import _subset_roi

        import anndata as ad
        import pandas as pd

        adata = ad.AnnData(
            X=np.ones((3, 2), dtype=np.float32),
            obs=pd.DataFrame({"roi": ["roi_a", "roi_a", "roi_b"]},
                             index=["roi_a__cell_1", "roi_a__cell_7", "roi_b__cell_2"]),
            var=pd.DataFrame(index=["G1", "G2"]),
        )
        subset = _subset_roi(adata, "roi_a")
        assert list(subset.obs_names) == ["cell_1", "cell_7"]

    def test_falls_back_to_obs_name_prefix(self):
        from backend.src.export.jobs import _subset_roi

        import anndata as ad
        import pandas as pd

        adata = ad.AnnData(
            X=np.ones((2, 2), dtype=np.float32),
            obs=pd.DataFrame(index=["roi_a__cell_1", "roi_b__cell_2"]),
            var=pd.DataFrame(index=["G1", "G2"]),
        )
        subset = _subset_roi(adata, "roi_b")
        assert list(subset.obs_names) == ["cell_2"]

    def test_unknown_roi_returns_empty(self):
        from backend.src.export.jobs import _subset_roi

        import anndata as ad
        import pandas as pd

        adata = ad.AnnData(
            X=np.ones((1, 2), dtype=np.float32),
            obs=pd.DataFrame(index=["roi_a__cell_1"]),
            var=pd.DataFrame(index=["G1", "G2"]),
        )
        assert len(_subset_roi(adata, "roi_zzz")) == 0


class TestCombineRoiPolygons:
    """Loupe 合併模式：全域偏移 + ROI 前綴"""

    def test_shifts_to_global_coords_and_prefixes_id(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs
        from backend.src.export.jobs import _combine_roi_polygons

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["roi_a__cell_1", "roi_b__cell_1"])
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})
        _write_mask(out / "roi" / "roi_b" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})

        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))
        assert inputs.is_merged is True

        combined = _combine_roi_polygons(inputs)
        features = json.loads(combined.read_text(encoding="utf-8"))["features"]

        assert combined.name == "combined_cellpose_polygons.json"
        assert {f["properties"]["full_id"] for f in features} == {"roi_a__1", "roi_b__1"}

        # roi_a 偏移 = x(100) * 0.5 = 50 µm；roi_b = 300 * 0.5 = 150 µm
        by_id = {f["properties"]["full_id"]: f for f in features}
        xs_a = [p[0] for p in by_id["roi_a__1"]["geometry"]["coordinates"][0]]
        xs_b = [p[0] for p in by_id["roi_b__1"]["geometry"]["coordinates"][0]]
        assert min(xs_a) == pytest.approx(50 + 5 * 0.5, abs=0.5)
        assert min(xs_b) == pytest.approx(150 + 5 * 0.5, abs=0.5)

    def test_missing_mask_is_skipped(self, tmp_path, export_config):
        from backend.src.export.inputs import resolve_export_inputs
        from backend.src.export.jobs import _combine_roi_polygons

        out = tmp_path / "analysis"
        _write_h5ad(out / "umap_computed.h5ad", ["roi_a__cell_1", "roi_b__cell_1"])
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})

        inputs = resolve_export_inputs(export_config, "", ("umap_computed.h5ad",))
        features = json.loads(_combine_roi_polygons(inputs).read_text(encoding="utf-8"))["features"]
        assert {f["properties"]["full_id"] for f in features} == {"roi_a__1"}


class TestLoupeJob:
    """run_loupe_export 端到端（假 exporter）"""

    def test_single_roi_writes_polygons_and_reports(self, tmp_path, export_config, monkeypatch):
        from backend.src.export import loupe_exporter
        from backend.src.export.jobs import run_loupe_export

        monkeypatch.setattr(loupe_exporter, "LoupeExporter", FakeExporter)

        out = tmp_path / "analysis"
        _write_h5ad(out / "clustered_final.h5ad", ["cell_1"], active_roi="roi_a")
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})

        seen: list = []
        result = run_loupe_export(export_config, progress=lambda f, m: seen.append((f, m)))

        assert result.is_merged is False
        assert result.output_dirs == [out / "roi" / "roi_a" / "export_loupe"]
        assert seen, "應回報進度"
        assert FakeExporter.instances[-1].kwargs["poly_json_path"] == out / "roi" / "roi_a" / "cellpose_polygons.json"

    def test_merged_uses_export_dir(self, tmp_path, export_config, monkeypatch):
        from backend.src.export import loupe_exporter
        from backend.src.export.jobs import run_loupe_export

        monkeypatch.setattr(loupe_exporter, "LoupeExporter", FakeExporter)

        out = tmp_path / "analysis"
        _write_h5ad(out / "clustered_final.h5ad", ["roi_a__cell_1", "roi_b__cell_1"])
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})
        _write_mask(out / "roi" / "roi_b" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})

        result = run_loupe_export(export_config)

        assert result.is_merged is True
        assert result.output_dirs == [tmp_path / "export" / "loupe"]
        assert (out / "combined_cellpose_polygons.json").exists()

    def test_progress_callback_is_optional(self, tmp_path, export_config, monkeypatch):
        """不給 callback 也要能跑（CLI 用途）"""
        from backend.src.export import loupe_exporter
        from backend.src.export.jobs import run_loupe_export

        monkeypatch.setattr(loupe_exporter, "LoupeExporter", FakeExporter)

        out = tmp_path / "analysis"
        _write_h5ad(out / "clustered_final.h5ad", ["cell_1"], active_roi="roi_a")
        _write_mask(out / "roi" / "roi_a" / "segmentation_masks.npy", {1: (5, 15, 5, 15)})

        assert run_loupe_export(export_config).output_dirs


class TestApiLayerIsThin:
    """回歸：API 只剩狀態與執行緒調度"""

    def test_no_orchestration_left_in_api(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "export.py").read_text(encoding="utf-8")
        for gone in ("mask_to_geojson", "resolve_export_inputs", "generate_visiumhd_transcripts",
                     "XeniumExporter", "LoupeExporter", "segmentation_masks.npy"):
            assert gone not in src, f"{gone} 不應再出現在 API 層"

    def test_api_delegates_to_jobs(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "export.py").read_text(encoding="utf-8")
        assert "run_xenium_export" in src
        assert "run_loupe_export" in src

    def test_jobs_module_has_no_event_loop_dependency(self):
        """領域層純同步：不得 import asyncio / fastapi"""
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "export" / "jobs.py").read_text(encoding="utf-8")
        assert "import asyncio" not in src
        assert "fastapi" not in src
