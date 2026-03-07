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
    """讀取 state.json（不存在時回傳空字典）。"""
    if not _STATE_PATH.exists():
        return {}
    with open(_STATE_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def save_state(updates: dict[str, Any]) -> None:
    """
    將動態狀態寫入 state.json。
    只傳入需要更新的部分，其餘已有狀態不受影響。
    """
    state = load_state()
    state = _deep_merge(state, updates)
    _STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(_STATE_PATH, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    logger.info(f"已更新 state.json：{list(updates.keys())}")


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """
    讀取 pipeline.yaml 並與 state.json 合併後回傳。
    state.json 的值優先於 pipeline.yaml。

    Parameters
    ----------
    path : str or Path, optional
        YAML 路徑。若為 None 則使用 config/pipeline.yaml。
    """
    config_path = Path(path) if path else _DEFAULT_CONFIG_PATH

    if not config_path.exists():
        raise FileNotFoundError(f"找不到設定檔：{config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

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
