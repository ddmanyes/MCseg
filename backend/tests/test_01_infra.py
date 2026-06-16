"""Test 1: Config / Constants / Discovery — 基礎設施測試"""
import pytest
from pathlib import Path


# ─── Config ─────────────────────────────────────────────────

class TestConfig:
    """pipeline.yaml 載入與存取"""

    def test_load_config(self, config_dict):
        """config 可成功載入"""
        assert isinstance(config_dict, dict)
        assert "global" in config_dict
        assert "paths" in config_dict

    def test_paths_section(self, config_dict):
        """paths 區塊包含必要欄位"""
        paths = config_dict["paths"]
        required = ["data_root", "he_image", "binned_002", "binned_008", "xenium_outs"]
        for key in required:
            assert key in paths, f"Missing key: {key}"

    def test_paths_not_empty(self, config_dict):
        """CRC 資料路徑已填入（不為空）"""
        paths = config_dict["paths"]
        assert paths["he_image"], "he_image path is empty"
        assert paths["binned_002"], "binned_002 path is empty"
        assert paths["binned_008"], "binned_008 path is empty"

    def test_rois_defined(self, config_dict):
        """ROI 列表已定義（數量取決於 state.json，至少 1 個）"""
        rois = config_dict.get("rois", [])
        assert len(rois) >= 1, f"Expected >= 1 ROI, got {len(rois)}"

    def test_roi_structure(self, config_dict):
        """每個 ROI 包含必要欄位"""
        rois = config_dict["rois"]
        for roi in rois:
            assert "name" in roi
            assert "x" in roi
            assert "y" in roi
            assert "width_px" in roi
            assert "height_px" in roi

    def test_segmentation_config(self, config_dict):
        """分割設定使用 MCseg v2（cyto3 多直徑集成）"""
        seg = config_dict["segmentation"]
        assert "mcseg_v2" in seg, "缺少 segmentation.mcseg_v2 區塊"
        mcseg = seg["mcseg_v2"]
        for key in ("dia_small", "dia_mid", "dia_large", "voronoi_distance"):
            assert key in mcseg, f"mcseg_v2 缺少必要欄位：{key}"

    def test_analysis_mito_prefix(self, config_dict):
        """人類資料必須用大寫 MT-"""
        prefix = config_dict["analysis"]["preprocessing"]["cellular"]["mito_prefix"]
        assert prefix == "MT-", f"Expected 'MT-', got '{prefix}'"


# ─── Constants ──────────────────────────────────────────────

class TestConstants:
    """物理常數集中管理"""

    def test_constants_importable(self):
        from backend.src.utils.constants import XENIUM_UM_PX, VISIUM_UM_PX, PROSEG_UM_PX
        assert XENIUM_UM_PX == 0.2125
        assert VISIUM_UM_PX == 0.2737
        assert PROSEG_UM_PX == 0.2645833

    def test_golden_params(self):
        from backend.src.utils.constants import GOLDEN_PARAMS
        assert "dilation" in GOLDEN_PARAMS
        assert "max_dist" in GOLDEN_PARAMS
        assert "compactness" in GOLDEN_PARAMS


# ─── Discovery ──────────────────────────────────────────────

class TestDiscovery:
    """自動掃描邏輯"""

    def test_scan_crc_root(self, crc_data_root):
        """掃描 CRC 資料根目錄"""
        from backend.src.utils.discovery import scan_data_root
        result = scan_data_root(str(crc_data_root))
        assert result is not None

    def test_scan_finds_btf(self, crc_data_root):
        """掃描能找到 BTF 影像"""
        from backend.src.utils.discovery import scan_data_root
        result = scan_data_root(str(crc_data_root))
        assert result.he_image is not None, "Should find H&E BTF image"
        assert "btf" in result.he_image.path.lower() or "tif" in result.he_image.path.lower()

    def test_scan_finds_binned(self, crc_data_root):
        """掃描能找到 binned 目錄"""
        from backend.src.utils.discovery import scan_data_root
        result = scan_data_root(str(crc_data_root))
        assert result.binned_002 is not None, "Should find square_002um"
        assert result.binned_008 is not None, "Should find square_008um"

    def test_scan_finds_xenium(self, crc_data_root):
        """掃描能找到 Xenium outs"""
        from backend.src.utils.discovery import scan_data_root
        result = scan_data_root(str(crc_data_root))
        assert result.xenium_outs is not None, "Should find Xenium outs"

    def test_scan_nonexistent_dir(self):
        """掃描不存在的目錄應回傳空結果 + warning"""
        from backend.src.utils.discovery import scan_data_root
        result = scan_data_root("/nonexistent/path/12345")
        assert result.he_image is None
        assert result.binned_002 is None
        assert len(result.warnings) > 0
