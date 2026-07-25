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


class FakeCellposeModel:
    """以固定規則產生 label 的假模型 —— 讓分割流程可在無 GPU/無模型下測試。

    規則：把輸入切成 32×32 網格，每格中央 16×16 給一個 label。結果只取決於
    tile 尺寸，因此「同一份影像走不同讀取路徑」必須得到相同輸出。
    """

    def __init__(self, *args, **kwargs):
        pass

    def eval(self, img, diameter=None, **kwargs):
        import numpy as np

        h, w = img.shape[:2]
        m = np.zeros((h, w), dtype=np.int32)
        lbl = 1
        for y in range(0, h - 31, 32):
            for x in range(0, w - 31, 32):
                m[y + 8:y + 24, x + 8:x + 24] = lbl
                lbl += 1
        return m, None, None


@pytest.fixture
def fake_cellpose(monkeypatch):
    """替換 cellpose 模型與 GPU 偵測，讓分割流程可離線執行。"""
    from cellpose import core, models

    monkeypatch.setattr(models, "CellposeModel", FakeCellposeModel)
    monkeypatch.setattr(core, "use_gpu", lambda *a, **k: False)
    return FakeCellposeModel


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
