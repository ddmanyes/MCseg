"""
統一日誌設定 — 同時輸出到 console 和 WebSocket 廣播佇列
"""
import asyncio
import logging
import sys
from collections import defaultdict
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


def setup_logging(log_level: str = "INFO") -> None:
    """
    初始化 pipeline 日誌系統。
    呼叫一次即可，FastAPI 啟動時呼叫。
    """
    level = getattr(logging, log_level.upper(), logging.INFO)

    fmt = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    root = logging.getLogger("pipeline")
    root.setLevel(level)
    root.handlers.clear()

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    root.addHandler(ch)

    # WebSocket handler
    ws_handler = WebSocketHandler(stage_getter=lambda: _current_stage)
    ws_handler.setLevel(level)
    ws_handler.setFormatter(fmt)
    root.addHandler(ws_handler)
