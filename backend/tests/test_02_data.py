"""Test 2: Data Files — CRC 資料完整性驗證"""
import pytest
from pathlib import Path

from .conftest import CRC_BTF, CRC_BINNED_002, CRC_BINNED_008, CRC_XENIUM


class TestCrcDataIntegrity:
    """驗證 CRC 資料檔案存在且結構正確"""

    # --- H&E 影像 ---

    def test_btf_exists(self):
        """H&E BTF 影像存在"""
        assert CRC_BTF.exists(), f"BTF not found: {CRC_BTF}"

    def test_btf_size(self):
        """BTF 大小合理（>1 GB）"""
        size_gb = CRC_BTF.stat().st_size / (1024**3)
        assert size_gb > 1.0, f"BTF too small: {size_gb:.2f} GB"

    def test_btf_is_tiff(self):
        """BTF 為有效的 TIFF 格式"""
        import tifffile
        with tifffile.TiffFile(CRC_BTF) as tf:
            page = tf.pages[0]
            assert page.shape[0] > 1000, "Image height too small"
            assert page.shape[1] > 1000, "Image width too small"
            assert page.tilewidth > 0, "Not a tiled TIFF"
            assert page.tilelength > 0, "Not a tiled TIFF"

    # --- Visium HD Binned Outputs ---

    def test_binned_002_exists(self):
        """square_002um 目錄存在"""
        assert CRC_BINNED_002.exists(), f"Not found: {CRC_BINNED_002}"
        assert CRC_BINNED_002.is_dir()

    def test_binned_008_exists(self):
        """square_008um 目錄存在"""
        assert CRC_BINNED_008.exists(), f"Not found: {CRC_BINNED_008}"
        assert CRC_BINNED_008.is_dir()

    def test_binned_002_has_h5(self):
        """square_002um 包含 filtered_feature_bc_matrix.h5"""
        h5 = CRC_BINNED_002 / "filtered_feature_bc_matrix.h5"
        assert h5.exists(), f"Missing: {h5}"

    def test_binned_002_has_spatial(self):
        """square_002um 包含 spatial/ 目錄"""
        spatial = CRC_BINNED_002 / "spatial"
        assert spatial.exists() and spatial.is_dir()

    def test_scalefactors_json(self):
        """spatial/scalefactors_json.json 存在"""
        sf = CRC_BINNED_002 / "spatial" / "scalefactors_json.json"
        assert sf.exists(), f"Missing: {sf}"
        import json
        with open(sf) as f:
            data = json.load(f)
        assert "microns_per_pixel" in data or "spot_diameter_fullres" in data

    def test_tissue_positions(self):
        """spatial/tissue_positions.parquet 存在"""
        tp = CRC_BINNED_002 / "spatial" / "tissue_positions.parquet"
        assert tp.exists(), f"Missing: {tp}"

    # --- Xenium Outs ---

    def test_xenium_outs_exists(self):
        """Xenium outs 目錄存在"""
        assert CRC_XENIUM.exists(), f"Not found: {CRC_XENIUM}"

    def test_xenium_transcripts(self):
        """Xenium transcripts.parquet 存在"""
        tx = CRC_XENIUM / "transcripts.parquet"
        assert tx.exists(), f"Missing: {tx}"

    def test_xenium_cells(self):
        """Xenium cells.parquet 存在"""
        cells = CRC_XENIUM / "cells.parquet"
        assert cells.exists(), f"Missing: {cells}"

    def test_xenium_experiment(self):
        """experiment.xenium 存在"""
        exp = CRC_XENIUM / "experiment.xenium"
        assert exp.exists(), f"Missing: {exp}"

    def test_xenium_morphology(self):
        """morphology.ome.tif 存在"""
        morph = CRC_XENIUM / "morphology.ome.tif"
        assert morph.exists(), f"Missing: {morph}"

    def test_xenium_nucleus_boundaries(self):
        """nucleus_boundaries.parquet 存在"""
        nb = CRC_XENIUM / "nucleus_boundaries.parquet"
        assert nb.exists(), f"Missing: {nb}"
