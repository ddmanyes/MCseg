"""Test 18: 全片分割編排邏輯（`fullslide.pipeline.run_full_slide_segmentation`）

架構深化 P8：把 `api/segmentation.py::_run_full_segmentation` 裡的 OOM 防護、
MPS batch_size 鉗制、force_disable_cpsam 決策，下沉成純同步、不依賴 FastAPI
的領域函式，讓 API 與 CLI 共用同一套安全防護（見 docs/adr/0005）。
"""
from pathlib import Path

import numpy as np
import pytest


def _tiled_cfg(**overrides):
    cfg = {
        "use_gpu": False, "batch_size": 1,
        "dia_small": 13.0, "dia_mid": 17.0, "dia_large": 22.0,
        "use_hematoxylin": False, "use_cpsam": False,
        "voronoi_distance": 4, "min_size": 20, "max_size": 6000,
        "clahe_clip_limit": 3.0,
    }
    cfg.update(overrides)
    return cfg


def _write_btf(path, size=256):
    import tifffile

    rng = np.random.default_rng(0)
    img = rng.integers(30, 180, (size, size, 3), dtype=np.uint8)
    tifffile.imwrite(str(path), img, bigtiff=True, tile=(128, 128), photometric="rgb")
    return path


def _base_config(tmp_path, btf_path, **full_seg_overrides):
    full_seg = {"tile_size": 128, "overlap": 32, "max_load_gb": 6.0,
                "force_disable_cpsam": True}
    full_seg.update(full_seg_overrides)
    return {
        "paths": {"he_image": str(btf_path), "output_dir": str(tmp_path / "out"),
                   "binned_002": ""},
        "segmentation": {"mcseg_v2": _tiled_cfg()},
        "full_seg": full_seg,
    }


class TestOomGuard:
    def test_oversized_mask_raises_memory_error(self, tmp_path, monkeypatch):
        """遮罩超過 max_load_gb 必須明確報錯，而非默默 OOM"""
        from backend.src.fullslide import pipeline
        from backend.src.utils import config as cfg_mod

        monkeypatch.setattr(cfg_mod, "resolve_path", lambda p: Path(p))

        btf = _write_btf(tmp_path / "slide.btf")
        config = _base_config(tmp_path, btf, max_load_gb=1e-6)

        with pytest.raises(MemoryError, match="請縮小"):
            pipeline.run_full_slide_segmentation(config)


class TestMpsSafetyClamp:
    def test_batch_size_clamped_to_two_even_when_requested_higher(self, tmp_path, monkeypatch):
        """batch_size 無條件鉗制到 2（docs/adr/0005：API 與 CLI 一致套用）"""
        from backend.src.fullslide import pipeline
        from backend.src.segmentation import cellpose_runner
        from backend.src.utils import config as cfg_mod

        monkeypatch.setattr(cfg_mod, "resolve_path", lambda p: Path(p))

        seen = {}

        def spy(cfg=None, **kwargs):
            seen["batch_size"] = cfg.get("batch_size")
            return np.zeros((128, 128), dtype=np.int32)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        btf = _write_btf(tmp_path / "slide.btf", size=128)
        config = _base_config(tmp_path, btf)
        config["segmentation"]["mcseg_v2"]["batch_size"] = 8

        pipeline.run_full_slide_segmentation(config)

        assert seen["batch_size"] == 2


class TestForceDisableCpsam:
    """`full_seg.force_disable_cpsam` 與請求的 `use_cpsam` 之間的三種分支"""

    def _run_and_capture_use_cpsam(self, tmp_path, monkeypatch, *, force_disable, use_cpsam_arg):
        from backend.src.fullslide import pipeline
        from backend.src.segmentation import cellpose_runner
        from backend.src.utils import config as cfg_mod

        monkeypatch.setattr(cfg_mod, "resolve_path", lambda p: Path(p))

        seen = {}

        def spy(cfg=None, **kwargs):
            seen["use_cpsam"] = cfg.get("use_cpsam")
            return np.zeros((128, 128), dtype=np.int32)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        btf = _write_btf(tmp_path / "slide.btf", size=128)
        config = _base_config(tmp_path, btf, force_disable_cpsam=force_disable)

        pipeline.run_full_slide_segmentation(config, use_cpsam=use_cpsam_arg)
        return seen["use_cpsam"]

    def test_force_disable_true_ignores_requested_use_cpsam(self, tmp_path, monkeypatch):
        """force_disable_cpsam=True 時無視請求的 use_cpsam，強制關閉"""
        result = self._run_and_capture_use_cpsam(
            tmp_path, monkeypatch, force_disable=True, use_cpsam_arg=True
        )
        assert result is False

    def test_force_disable_false_honors_requested_true(self, tmp_path, monkeypatch):
        """force_disable_cpsam=False 時，明確請求的 use_cpsam=True 生效"""
        result = self._run_and_capture_use_cpsam(
            tmp_path, monkeypatch, force_disable=False, use_cpsam_arg=True
        )
        assert result is True

    def test_force_disable_false_and_unset_falls_back_to_profile(self, tmp_path, monkeypatch):
        """force_disable_cpsam=False 且未指定 use_cpsam 時，沿用 mcseg_v2.use_cpsam"""
        result = self._run_and_capture_use_cpsam(
            tmp_path, monkeypatch, force_disable=False, use_cpsam_arg=None
        )
        # _tiled_cfg() 的 mcseg_v2.use_cpsam 預設 False
        assert result is False
