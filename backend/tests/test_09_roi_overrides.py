"""Test 9: ROI 分割參數覆寫的驗證入口（不需真實資料）"""
import json
import pytest


class TestValidateRoiOverrides:
    """validate_roi_overrides 是「哪些欄位可覆寫」的對外唯一入口"""

    def test_splits_clean_and_invalid(self):
        from backend.src.segmentation.cellpose_runner import validate_roi_overrides

        clean, invalid = validate_roi_overrides({
            "dia_mid": 30,
            "voronoi_distance": 5,
            "not_a_field": 1,
            "another_bogus": "x",
        })
        assert clean == {"dia_mid": 30, "voronoi_distance": 5}
        assert invalid == ["not_a_field", "another_bogus"]   # 保持輸入順序

    def test_none_values_dropped_from_clean_but_not_invalid(self):
        """值為 None = 前端沒填，不算未知欄位，但也不套用"""
        from backend.src.segmentation.cellpose_runner import validate_roi_overrides

        clean, invalid = validate_roi_overrides({"dia_small": None, "min_size": 40})
        assert clean == {"min_size": 40}
        assert invalid == []

    def test_empty_input(self):
        from backend.src.segmentation.cellpose_runner import validate_roi_overrides

        assert validate_roi_overrides({}) == ({}, [])


class TestMergeRoiParams:
    """_merge_roi_params 透過同一個驗證入口套用覆寫"""

    def test_overrides_applied_to_mcseg_section(self):
        from backend.src.segmentation.cellpose_runner import _merge_roi_params

        base = {"mcseg_v2": {"dia_mid": 20, "min_size": 10}, "other": {"keep": True}}
        merged = _merge_roi_params(base, {"dia_mid": 35})

        assert merged["mcseg_v2"]["dia_mid"] == 35
        assert merged["mcseg_v2"]["min_size"] == 10
        assert merged["other"] == {"keep": True}

    def test_does_not_mutate_input(self):
        from backend.src.segmentation.cellpose_runner import _merge_roi_params

        base = {"mcseg_v2": {"dia_mid": 20}}
        _merge_roi_params(base, {"dia_mid": 99})
        assert base["mcseg_v2"]["dia_mid"] == 20

    def test_unknown_fields_ignored_and_logged(self, caplog):
        from backend.src.segmentation.cellpose_runner import _merge_roi_params

        with caplog.at_level("WARNING", logger="pipeline.segmentation"):
            merged = _merge_roi_params({"mcseg_v2": {}}, {"bogus": 1, "dia_mid": 25})

        assert "bogus" not in merged["mcseg_v2"]
        assert merged["mcseg_v2"]["dia_mid"] == 25
        assert "bogus" in caplog.text

    def test_none_override_does_not_clobber_global(self):
        from backend.src.segmentation.cellpose_runner import _merge_roi_params

        merged = _merge_roi_params({"mcseg_v2": {"dia_mid": 20}}, {"dia_mid": None})
        assert merged["mcseg_v2"]["dia_mid"] == 20


class TestApiUsesPublicEntry:
    """回歸：API 不得再 import 領域私有符號"""

    def test_api_does_not_import_private_field_set(self):
        from pathlib import Path

        src = (Path(__file__).resolve().parents[1] / "src" / "api" / "segmentation.py").read_text(encoding="utf-8")
        assert "_ROI_OVERRIDE_FIELDS" not in src
        assert "validate_roi_overrides" in src


class TestSaveStateKeyReplacesNotMerges:
    """回歸：save_state_key 必須整份取代，不能沿用 save_state 的深合併

    背景：put_roi_overrides 曾經呼叫 save_state({"roi_seg_overrides": body})，
    經 _deep_merge 遞迴合併——送 {}（全部重設）或少一個 key（刪除單一 ROI）
    都不會真的清空/刪除，回應仍是成功但後端參數其實沒變（見
    docs/brainstorming/tauri_desktop_packaging.md Phase 0）。
    """

    @pytest.fixture
    def state_path(self, tmp_path, monkeypatch):
        from backend.src.utils import config as config_module

        path = tmp_path / "state.json"
        monkeypatch.setattr(config_module, "_STATE_PATH", path)
        return path

    def test_empty_dict_clears_existing_overrides(self, state_path):
        """對照組：save_state（深合併）送空 dict 不會清空——證明舊行為確實是 bug"""
        from backend.src.utils.config import save_state, save_state_key

        save_state({"roi_seg_overrides": {"CRC_1": {"dia_mid": 20}}})
        save_state({"roi_seg_overrides": {}})
        assert json.loads(state_path.read_text())["roi_seg_overrides"] == {
            "CRC_1": {"dia_mid": 20}
        }  # 深合併對空 dict 沒有 key 可覆寫，舊值原封不動殘留（舊 bug 的行為）

        save_state_key("roi_seg_overrides", {"CRC_1": {"dia_mid": 20}})
        save_state_key("roi_seg_overrides", {})
        assert json.loads(state_path.read_text())["roi_seg_overrides"] == {}

    def test_removed_key_actually_deleted(self, state_path):
        """單一 ROI 的覆寫被移除後，重新寫入時不應該殘留"""
        from backend.src.utils.config import save_state_key

        save_state_key(
            "roi_seg_overrides",
            {"CRC_1": {"dia_mid": 20}, "CRC_2": {"dia_mid": 25}},
        )
        save_state_key("roi_seg_overrides", {"CRC_2": {"dia_mid": 25}})

        data = json.loads(state_path.read_text())["roi_seg_overrides"]
        assert "CRC_1" not in data
        assert data == {"CRC_2": {"dia_mid": 25}}

    def test_does_not_disturb_other_state_keys(self, state_path):
        """整份取代只動指定的 key，state.json 其他欄位不受影響"""
        from backend.src.utils.config import save_state, save_state_key

        save_state({"paths": {"data_root": "/some/path"}})
        save_state_key("roi_seg_overrides", {"CRC_1": {"dia_mid": 20}})

        data = json.loads(state_path.read_text())
        assert data["paths"] == {"data_root": "/some/path"}
        assert data["roi_seg_overrides"] == {"CRC_1": {"dia_mid": 20}}
