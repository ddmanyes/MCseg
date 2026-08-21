"""測試 msseg-segment / msseg-analyze skill 支援腳本"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# 確保專案根目錄在 sys.path
ROOT = Path(__file__).parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ─── 共用 Fixtures ────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_project(tmp_path):
    """建立暫時專案目錄結構（模擬 MSseg 執行環境）"""
    (tmp_path / "results" / "qc").mkdir(parents=True)
    (tmp_path / "results" / "analysis").mkdir(parents=True)
    (tmp_path / "results" / "export").mkdir(parents=True)
    (tmp_path / "config").mkdir()
    (tmp_path / "scripts").mkdir()
    return tmp_path


@pytest.fixture
def fake_tissue_positions(tmp_project):
    """建立假 tissue_positions.parquet：格狀密集 bins 確保覆蓋率計算可找到有效 ROI"""
    binned_dir = tmp_project / "binned_002"
    binned_dir.mkdir()

    # 在 2000–5000 範圍內每 16px 放一個 bin（格狀，~35,000 bins）
    step = 16
    rs, cs = np.meshgrid(
        np.arange(2000, 5000, step),
        np.arange(2000, 5000, step),
        indexing="ij",
    )
    rows = rs.ravel()
    cols = cs.ravel()
    n = len(rows)
    df = pd.DataFrame({
        "barcode": [f"AAAA{i:06d}-1" for i in range(n)],
        "in_tissue": [1] * n,
        "array_row": (rows // 16).astype(int),
        "array_col": (cols // 16).astype(int),
        "pxl_row_in_fullres": rows.astype(int),
        "pxl_col_in_fullres": cols.astype(int),
    })
    df.to_parquet(binned_dir / "tissue_positions.parquet", index=False)
    return binned_dir


@pytest.fixture
def fake_masks(tmp_project):
    """建立 2 個假 ROI 遮罩（100×100, 5 個細胞）"""
    qc_dir = tmp_project / "results" / "qc"
    masks = {}
    for name in ["qc_roi_1", "qc_roi_2"]:
        mask = np.zeros((100, 100), dtype=np.int32)
        # 放入 5 個 10×10 的方形細胞
        for cell_id, (r, c) in enumerate([(5,5),(5,50),(50,5),(50,50),(25,25)], start=1):
            mask[r:r+10, c:c+10] = cell_id
        np.save(qc_dir / f"{name}_nuc.npy", mask)
        np.save(qc_dir / f"{name}_mcseg.npy", mask)
        masks[name] = mask
    return masks


@pytest.fixture
def fake_adata(tmp_project):
    """建立假 AnnData（10 個 bins，20 個基因，含空間座標）"""
    n_bins, n_genes = 10, 20
    X = sp.csr_matrix(np.random.poisson(2, (n_bins, n_genes)).astype(np.float32))
    gene_names = [f"GENE{i:02d}" for i in range(n_genes)]
    # 加入測試用不可能基因對
    gene_names[0] = "EPCAM"
    gene_names[1] = "CD3E"
    gene_names[2] = "MUC2"
    gene_names[3] = "NKG7"

    obs = pd.DataFrame({
        "barcode": [f"AAAA{i:06d}-1" for i in range(n_bins)],
        "in_tissue": [1] * n_bins,
        "pxl_row_in_fullres": np.random.randint(0, 90, n_bins),
        "pxl_col_in_fullres": np.random.randint(0, 90, n_bins),
    }).set_index("barcode")

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))
    adata_path = tmp_project / "binned_002" / "adata_002um.h5ad"
    adata_path.parent.mkdir(exist_ok=True)
    adata.write(adata_path)
    return adata, adata_path


@pytest.fixture
def fake_qc_metrics_csv(tmp_project):
    """建立假 qc_metrics.csv"""
    df = pd.DataFrame([
        {"roi": "qc_roi_1", "method": "nuc",   "n_cells": 120, "ftc": 0.70, "ned": 0.681, "coexp_rate": 0.025},
        {"roi": "qc_roi_1", "method": "mcseg", "n_cells": 145, "ftc": 0.82, "ned": 0.724, "coexp_rate": 0.021},
        {"roi": "qc_roi_2", "method": "nuc",   "n_cells": 110, "ftc": 0.68, "ned": 0.675, "coexp_rate": 0.028},
        {"roi": "qc_roi_2", "method": "mcseg", "n_cells": 138, "ftc": 0.80, "ned": 0.718, "coexp_rate": 0.024},
    ])
    csv_path = tmp_project / "results" / "qc_metrics.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


# ─── A1: roi_sampler ─────────────────────────────────────────────────────────

class TestRoiSampler:
    def test_roi_sampler_produces_output(self, tmp_project, fake_tissue_positions):
        """roi_sampler 應產出 ≥1 個 ROI，且每個 ROI 有 coverage 欄位"""
        from scripts.roi_sampler import sample_rois

        out_path = tmp_project / "results" / "qc_rois.json"
        result = sample_rois(
            binned_dir=fake_tissue_positions,
            out_path=out_path,
        )

        assert "rois" in result
        assert len(result["rois"]) >= 1
        assert out_path.exists()
        for roi in result["rois"]:
            assert "coverage" in roi
            assert 0 < roi["coverage"] <= 1.0
            assert "x" in roi and "y" in roi
            assert "width_px" in roi and "height_px" in roi


# ─── A2: seg_quality ─────────────────────────────────────────────────────────

class TestSegQuality:
    def test_seg_quality_creates_masks(self, tmp_project):
        """seg_quality 對假影像應產出 _nuc.npy 與 _mcseg.npy"""
        import tifffile
        from scripts.seg_quality import run_seg_quality

        qc_dir = tmp_project / "results" / "qc"

        # 假 qc_rois.json
        rois_data = {
            "rois": [{"name": "qc_roi_1", "x": 0, "y": 0,
                      "width_px": 64, "height_px": 64,
                      "pixel_size_um": 0.2737, "coverage": 0.9}],
            "threshold_used": 0.6,
            "timestamp": "2026-01-01T00:00:00",
        }
        rois_json = tmp_project / "results" / "qc_rois.json"
        rois_json.write_text(json.dumps(rois_data), encoding="utf-8")

        # 假 TIFF（64×64 灰色影像）
        fake_img = np.ones((128, 128, 3), dtype=np.uint8) * 200
        btf_path = tmp_project / "fake.tif"
        tifffile.imwrite(str(btf_path), fake_img)

        mcseg_cfg = {
            "use_gpu": False, "batch_size": 1,
            "dia_small": 13.0, "dia_mid": 17.0, "dia_large": 22.0,
            "use_hematoxylin": False, "use_cpsam": False,
            "voronoi_distance": 0, "clahe_clip_limit": 1.0,
            "min_size": 5, "max_size": 6000,
            "flow_threshold": 0.4, "cellprob_threshold": -2.0,
            "use_transcript_rescue": False,
        }

        run_seg_quality(
            rois_json=rois_json,
            he_path=btf_path,
            mcseg_cfg=mcseg_cfg,
            out_dir=qc_dir,
        )

        assert (qc_dir / "qc_roi_1_nuc.npy").exists(), "_nuc.npy 不存在"
        assert (qc_dir / "qc_roi_1_mcseg.npy").exists(), "_mcseg.npy 不存在"


# ─── A3: qc_metrics ──────────────────────────────────────────────────────────
#
# `TestQcMetrics::test_qc_metrics_columns` 已刪除（架構深化 P6-6 收尾時清理）：
# `3649e58` 把 `scripts/qc_metrics.py` 的舊 `compute_metrics(rois_json=...)` 介面
# 整個改寫成新的 `compute_metrics_from_geojson()`（geojson 導向），但沒有同步
# 更新這條測試，留下一個測不存在函式的孤兒測試。新介面目前沒有對應測試覆蓋，
# 不在本次架構深化範圍內，未來若要動 `qc_metrics.py` 建議先補上。

# ─── A4: write_handoff ───────────────────────────────────────────────────────

class TestWriteHandoff:
    def test_handoff_report_keys(self, tmp_project):
        """handoff_report.json 應包含六個必要頂層鍵"""
        from scripts.write_handoff import write_handoff

        # 建立假 qc_metrics.csv
        csv_path = tmp_project / "results" / "qc_metrics.csv"
        pd.DataFrame([
            {"roi": "qc_roi_1", "method": "nuc",   "n_cells": 120, "ftc": 0.70, "ned": 0.681, "coexp_rate": 0.025},
            {"roi": "qc_roi_1", "method": "mcseg", "n_cells": 145, "ftc": 0.82, "ned": 0.724, "coexp_rate": 0.021},
        ]).to_csv(csv_path, index=False)

        out_path = tmp_project / "results" / "handoff_report.json"
        write_handoff(
            metrics_csv=csv_path,
            masks_dir=str(tmp_project / "results" / "masks"),
            binned_dir=str(tmp_project / "binned_002"),
            tissue_profile="crc",
            out_path=out_path,
        )

        assert out_path.exists(), "handoff_report.json 不存在"
        report = json.loads(out_path.read_text())

        required_keys = {
            "segmentation_complete", "roi_qc",
            "recommended_analysis_params",
            "masks_dir", "binned_dir", "tissue_profile",
        }
        assert required_keys.issubset(report.keys()), \
            f"缺少頂層鍵：{required_keys - set(report.keys())}"
        assert report["segmentation_complete"] is True


# ─── B1: build_full_adata ────────────────────────────────────────────────────

class TestBuildFullAdata:
    def test_build_adata_shape(self, tmp_project, fake_adata):
        """build_full_adata 應產出 n_obs>=1 的 h5ad"""
        import yaml

        adata, adata_path = fake_adata

        # 建立假遮罩（50×50, 3 cells），對應 fake_adata 的 bin 座標
        mask = np.zeros((100, 100), dtype=np.int32)
        mask[5:20, 5:20] = 1
        mask[5:20, 50:65] = 2
        mask[50:65, 5:20] = 3
        mask_path = tmp_project / "results" / "masks" / "qc_roi_1_mcseg.npy"
        mask_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(mask_path, mask)

        report = {
            "segmentation_complete": True,
            "masks_dir": str(mask_path.parent),
            "binned_dir": str(adata_path.parent),
            "tissue_profile": "crc",
            "roi_qc": {"ned_mcseg": 0.72, "ned_nuc": 0.68,
                       "ned_delta": 0.04, "ftc_mean": 0.82, "coexp_mean": 0.022},
            "recommended_analysis_params": {"min_genes": 1, "max_pct_mt": 20, "min_counts": 1},
            "n_rois_evaluated": 1,
        }
        (tmp_project / "results" / "handoff_report.json").write_text(json.dumps(report))

        rois_data = {"rois": [{"name": "qc_roi_1", "x": 0, "y": 0,
                               "width_px": 100, "height_px": 100, "pixel_size_um": 0.2737}]}
        (tmp_project / "results" / "qc_rois.json").write_text(json.dumps(rois_data))

        cfg = {
            "global": {"tissue_profile": "crc"},
            "paths": {
                "masks_dir": str(mask_path.parent),
                "binned_002": str(adata_path.parent),
            },
        }
        (tmp_project / "config" / "pipeline.yaml").write_text(yaml.dump(cfg))

        from scripts.build_full_adata import build_adata

        out_path = tmp_project / "results" / "analysis" / "cellpose_cells.h5ad"
        result = build_adata(
            handoff_json=tmp_project / "results" / "handoff_report.json",
            rois_json=tmp_project / "results" / "qc_rois.json",
            out_path=out_path,
        )

        assert out_path.exists(), "cellpose_cells.h5ad 不存在"
        assert result.n_obs >= 1, "輸出 AnnData 沒有細胞"
        assert result.n_vars == 20, f"基因數應為 20，得到 {result.n_vars}"


# ─── B2: run_analysis ────────────────────────────────────────────────────────

class TestRunAnalysis:
    def test_run_analysis_umap_exists(self, tmp_project, monkeypatch):
        """run_analysis 輸出應包含 X_umap 與 leiden 欄位"""
        import importlib.util, yaml, subprocess

        # 建立 100 細胞假 AnnData（需足夠多基因才能跑 HVG）
        np.random.seed(0)
        n_cells, n_genes = 100, 300
        X = sp.csr_matrix(np.random.poisson(3, (n_cells, n_genes)).astype(np.float32))
        adata = ad.AnnData(X=X,
                           obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_cells)]),
                           var=pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)]))
        input_path = tmp_project / "results" / "analysis" / "cellpose_cells.h5ad"
        adata.write(input_path)

        monkeypatch.chdir(tmp_project)

        result = __import__("subprocess").run(
            ["uv", "run", "python", str(ROOT / "scripts" / "run_analysis.py"),
             "--input", str(input_path),
             "--output-dir", str(tmp_project / "results" / "analysis"),
             "--min-genes", "1", "--max-pct-mt", "100", "--min-counts", "1",
             "--resolution", "0.5"],
            capture_output=True, text=True, cwd=str(ROOT)
        )
        assert result.returncode == 0, f"run_analysis.py 失敗：\n{result.stderr}"

        out_path = tmp_project / "results" / "analysis" / "cellpose_cells_clustered.h5ad"
        assert out_path.exists(), "clustered h5ad 不存在"
        clustered = ad.read_h5ad(out_path)
        assert "X_umap" in clustered.obsm, "缺少 X_umap"
        assert "leiden" in clustered.obs.columns, "缺少 leiden 欄位"


# ─── B4: export_mcseg ────────────────────────────────────────────────────────
#
# `TestExportMcseg::test_export_cli_both_format` 已刪除（架構深化 P6-6 收尾時清理）：
# `3649e58` 刪除了整個 `scripts/export_mcseg.py`（連同 `skills/scripts/` 那份鏡像），
# 但沒有同步移除這條測試，留下一個測不存在檔案的孤兒測試。匯出功能看起來已被
# 同一次 commit 新增的 `generate_segmentation_overlays.py` / `render_registered_
# cytassist_overlays.py` 取代，不在本次架構深化範圍內。
