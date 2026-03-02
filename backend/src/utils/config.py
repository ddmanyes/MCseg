"""
配置管理：讀取、驗證 pipeline.yaml
"""
import logging
import os
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger("pipeline.config")

_DEFAULT_CONFIG_PATH = Path(__file__).parents[3] / "config" / "pipeline.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """
    讀取 pipeline.yaml 並回傳配置字典。

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

    logger.info(f"已載入設定檔：{config_path}")
    return config


def save_config(config: dict[str, Any], path: str | Path | None = None) -> None:
    """將配置寫回 YAML 檔案。"""
    config_path = Path(path) if path else _DEFAULT_CONFIG_PATH
    config_path.parent.mkdir(parents=True, exist_ok=True)

    with open(config_path, "w", encoding="utf-8") as f:
        yaml.dump(config, f, allow_unicode=True, sort_keys=False, default_flow_style=False)

    logger.info(f"已儲存設定檔：{config_path}")


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
