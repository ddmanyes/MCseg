"""
對位誤差的量化（bin 密度光柵化 + 次像素位移估計）
================================================

**定位（P0.5 之後）**：主要的座標對位由 `alignment.py` 的組合 homography 負責
（幾何正確、來源可追溯）。本模組是**驗證手段** —— 回答「JSON 套下去之後，
RNA 與影像還剩多少殘餘偏移？」20-50 px 的殘餘位移不會觸發
`bin_attribution` 的 30% 越界警告，卻足以讓 RNA 落到隔壁細胞。
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("pipeline.registration.align")
