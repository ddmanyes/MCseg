"""
對位 QC 疊圖
============

把 H&E 影像與 Visium bin 質心疊在同一張 patch 上輸出 PNG —— 對位有沒有問題，
肉眼一看便知，不必先相信任何估計數值。
"""
from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger("pipeline.registration.qc")
