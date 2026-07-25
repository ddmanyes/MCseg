"""
全片流程核心函式（CLI 與 API 共用）
=====================================

本模組收斂「MCseg v2 遮罩 + Visium HD bins → cells×genes」這條配方，
使 `cli/segment.py`（指令列全片流程）與 `api/cellpose_count.py`（GUI 全圖計數）
共用同一份實作，避免邏輯分歧（CLAUDE.md §11 DRY）。

函式一律為純函式（不寫 log、不做快取），快取與進度回報留給呼叫端：

| 函式 | 職責 |
|------|------|
| `bin_attribution`   | 2µm bins → cell_id 對應表 |
| `aggregate_cells`   | 依對應表把 bins 聚合成 cells×genes 原始 counts |
| `add_centroids`     | 補細胞重心（裁切局部 px、全片 fullres px、µm） |
| `resolve_pixel_size`| 取樣本實際 µm/px（scalefactors 優先於預設常數） |
"""
from __future__ import annotations
