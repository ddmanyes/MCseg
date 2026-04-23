"""
配置管理：讀取 pipeline.yaml（靜態設定）+ state.json（執行期狀態）

- pipeline.yaml：使用者設定，程式永不覆寫
- config/state.json：程式執行時寫入的動態狀態（paths、rois、params 等）
- load_config()：回傳兩者深度合併的結果（state 優先）
- save_state(updates)：只更新 state.json 的對應欄位
"""
import json
import logging
import os
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger("pipeline.config")

_DEFAULT_CONFIG_PATH = Path(__file__).parents[3] / "config" / "pipeline.yaml"
_STATE_PATH = Path(__file__).parents[3] / "config" / "state.json"
_PROFILES_DIR = Path(__file__).parents[3] / "config" / "profiles"


def _deep_merge(base: dict, override: dict) -> dict:
    """Deep merge：override 的值覆蓋 base，巢狀 dict 遞迴合併。"""
    result = base.copy()
    for k, v in override.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = v
    return result


def load_state() -> dict[str, Any]:
    """讀取 state.json（不存在或 JSON 損壞時回傳空字典）。"""
    if not _STATE_PATH.exists():
        return {}
    try:
        with open(_STATE_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        logger.warning(f"state.json 損壞，忽略並重置：{e}")
        return {}


def save_state(updates: dict[str, Any]) -> None:
    """
    將動態狀態寫入 state.json（原子寫入，防止寫到一半損壞）。
    只傳入需要更新的部分，其餘已有狀態不受影響。
    """
    state = load_state()
    state = _deep_merge(state, updates)
    _STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = _STATE_PATH.with_suffix(".tmp")
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    tmp_path.replace(_STATE_PATH)  # 原子替換，防止寫到一半損壞
    logger.info(f"已更新 state.json：{list(updates.keys())}")


def _load_profile(profile_name: str) -> dict[str, Any]:
    """
    載入 config/profiles/{profile_name}.yaml。

    若檔案不存在，回傳空字典並記錄警告。
    """
    import re as _re
    if not _re.match(r'^[\w\-]+$', str(profile_name)):
        logger.warning(f"tissue_profile 名稱含非法字元，略過：{profile_name!r}")
        return {}
    profile_path = _PROFILES_DIR / f"{profile_name}.yaml"
    if not profile_path.exists():
        logger.warning(f"找不到 tissue profile：{profile_path}，使用空 profile")
        return {}
    with open(profile_path, "r", encoding="utf-8") as f:
        profile = yaml.safe_load(f) or {}
    logger.info(f"已載入 tissue profile：{profile_name} ({profile_path.name})")
    return profile


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """
    讀取 pipeline.yaml 並與 tissue profile + state.json 合併後回傳。

    合併順序（後者覆寫前者）：
      profile.yaml（組織基底）← pipeline.yaml（專案覆寫）← state.json（執行期覆寫）

    Parameters
    ----------
    path : str or Path, optional
        YAML 路徑。若為 None 則使用 config/pipeline.yaml。
    """
    config_path = Path(path) if path else _DEFAULT_CONFIG_PATH

    if not config_path.exists():
        raise FileNotFoundError(f"找不到設定檔：{config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        pipeline_cfg = yaml.safe_load(f)

    # 載入 tissue profile（若 global.tissue_profile 有指定）
    profile_name = (pipeline_cfg.get("global") or {}).get("tissue_profile") or \
                   pipeline_cfg.get("tissue_profile")
    if profile_name:
        profile = _load_profile(str(profile_name))
        # pipeline.yaml 覆寫 profile（profile 為基底）
        config = _deep_merge(profile, pipeline_cfg)
    else:
        config = pipeline_cfg

    # state.json 覆寫最終 config
    state = load_state()
    if state:
        config = _deep_merge(config, state)

    logger.info(f"已載入設定檔：{config_path}")
    return config


def resolve_path(p: str, base: Path | None = None) -> Path:
    """
    將設定檔中的路徑轉為絕對路徑。
    - 若已是絕對路徑則直接回傳
    - 否則相對於 base（預設為專案根目錄）解析
    """
    path = Path(os.path.expanduser(p))
    if path.is_absolute():
        return path
    base = base or _DEFAULT_CONFIG_PATH.parents[1]
    return (base / path).resolve()


def get_roi_list(config: dict[str, Any]) -> list[dict]:
    """從 config 中取得 ROI 清單，自動過濾空項目。"""
    return [r for r in config.get("rois", []) if r.get("name")]
