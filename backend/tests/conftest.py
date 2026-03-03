"""pytest fixtures — 共用測試夾具"""
import pytest
from pathlib import Path

# === CRC 資料路徑 ===
CRC_ROOT = Path("/Volumes/SSD/plan_a/tissue sample/CRC")
CRC_VISIUM = CRC_ROOT / "visium/official_v4"
CRC_XENIUM = CRC_ROOT / "xenium/official_v1_addon/outs"

CRC_BTF = CRC_VISIUM / "Visium_HD_Human_Colon_Cancer_tissue_image.btf"
CRC_BINNED_002 = CRC_VISIUM / "binned_outputs/binned_outputs/square_002um"
CRC_BINNED_008 = CRC_VISIUM / "binned_outputs/binned_outputs/square_008um"


@pytest.fixture
def crc_data_root():
    """CRC 資料根目錄"""
    return CRC_ROOT


@pytest.fixture
def crc_btf_path():
    """CRC H&E BTF 影像路徑"""
    return CRC_BTF


@pytest.fixture
def crc_binned_002():
    """CRC Visium HD 2µm binned 目錄"""
    return CRC_BINNED_002


@pytest.fixture
def crc_binned_008():
    """CRC Visium HD 8µm binned 目錄"""
    return CRC_BINNED_008


@pytest.fixture
def crc_xenium_outs():
    """CRC Xenium outs 目錄"""
    return CRC_XENIUM


@pytest.fixture
def config_dict():
    """載入 pipeline.yaml 並回傳 dict"""
    from backend.src.utils.config import load_config
    return load_config()


@pytest.fixture
def crc_tumor_roi():
    """CRC 腫瘤邊界 ROI 參數"""
    return {
        "name": "CRC_tumor_boundary",
        "tissue": "CRC",
        "x": 43490,
        "y": 12515,
        "width_px": 6569,
        "height_px": 4791,
        "pixel_size_um": 0.2737,
    }
