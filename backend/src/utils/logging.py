"""
統一日誌設定 — 同時輸出到 console、WebSocket 廣播佇列，以及磁碟 log 檔
"""
import asyncio
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Callable


# ── 全域 WebSocket 廣播佇列（stage → list of queues）────────
_stage_queues: dict[str, list[asyncio.Queue]] = defaultdict(list)

_current_stage: str = "global"


def set_current_stage(stage: str) -> None:
    global _current_stage
    _current_stage = stage


def register_ws_queue(stage: str, queue: asyncio.Queue) -> None:
    _stage_queues[stage].append(queue)


def unregister_ws_queue(stage: str, queue: asyncio.Queue) -> None:
    if queue in _stage_queues[stage]:
        _stage_queues[stage].remove(queue)


class WebSocketHandler(logging.Handler):
    """將 log 記錄推入對應 stage 的 WebSocket 佇列。"""

    def __init__(self, stage_getter: Callable[[], str]):
        super().__init__()
        self._stage_getter = stage_getter

    def emit(self, record: logging.LogRecord) -> None:
        msg = self.format(record)
        stage = self._stage_getter()
        queues = _stage_queues.get(stage, []) + _stage_queues.get("global", [])
        for q in queues:
            try:
                q.put_nowait({"type": "log", "stage": stage, "message": msg, "level": record.levelname})
            except asyncio.QueueFull:
                pass


def setup_logging(log_level: str = "INFO", log_dir: str | Path | None = None, keep_logs: int = 30) -> None:
    """
    初始化 pipeline 日誌系統。
    呼叫一次即可，FastAPI 啟動時呼叫。

    Parameters
    ----------
    log_level : str
        日誌等級（DEBUG / INFO / WARNING / ERROR）
    log_dir : str or Path, optional
        log 檔存放目錄。None 時使用專案根目錄下的 logs/
    keep_logs : int
        保留最新幾份 log 檔，超過時自動刪除最舊的
    """
    from datetime import datetime

    level = getattr(logging, log_level.upper(), logging.INFO)

    fmt = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fmt_console = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    root = logging.getLogger("pipeline")
    root.setLevel(level)
    root.handlers.clear()

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt_console)
    root.addHandler(ch)

    # WebSocket handler
    ws_handler = WebSocketHandler(stage_getter=lambda: _current_stage)
    ws_handler.setLevel(level)
    ws_handler.setFormatter(fmt_console)
    root.addHandler(ws_handler)

    # File handler — 每次啟動產生一個新 log 檔
    _log_dir = Path(log_dir) if log_dir else Path(__file__).parents[3] / "logs"
    _log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = _log_dir / f"msseg_{timestamp}.log"

    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.DEBUG)   # 檔案永遠記錄 DEBUG，便於追蹤
    fh.setFormatter(fmt)
    root.addHandler(fh)

    root.info(f"Log file: {log_path}")

    # 超過 keep_logs 份時刪除最舊的
    existing = sorted(_log_dir.glob("msseg_*.log"))
    for old in existing[:-keep_logs]:
        try:
            old.unlink()
        except OSError:
            pass
