"""
所有物理常數集中管理
禁止在其他模組硬編碼這些數值
"""

# ── 像素尺寸（µm/px）──────────────────────────────────────────
XENIUM_UM_PX: float = 0.2125
"""Xenium morphology 影像像素尺寸（µm/px）"""

XENIUM_NM_PX: float = 212.5
"""Xenium morphology 影像像素尺寸（nm/px）"""

PROSEG_UM_PX: float = 0.2645833
"""Proseg 輸出座標系像素尺寸（µm/px）"""

PROSEG_NM_PX: float = 264.5833
"""Proseg 輸出座標系像素尺寸（nm/px）"""

VISIUM_UM_PX: float = 0.2737
"""Visium HD fullres 像素尺寸（µm/px），LUAD 和 CRC 均適用"""

# ── Proseg 黃金參數（可被 config.yaml 覆寫）──────────────────
GOLDEN_PARAMS: dict = {
    "dilation":     20,      # 核遮罩擴張半徑（px）
    "max_dist":     40.0,    # RNA 到核心的最大距離（µm）
    "compactness":  0.06,    # 細胞緊密度正則化
    "samples":      500,     # MCMC 迭代次數
    "recorded":     150,     # 記錄樣本數
    "watershed":    True,    # 使用分水嶺種子分配
    "connectivity": True,    # 強制細胞連通性
}

# ── Visium HD 縮放因子（tissue_hires）────────────────────────
# 讀自 scalefactors_json.json，此處為已知典型值（供參考）
HIRES_SCALEF_CRC:  float = 0.07973422
HIRES_SCALEF_LUAD: float = 0.1386642

# ── Xenium Browser pixel_size（修補用）───────────────────────
XENIUM_EXPLORER_HARDCODED_UM_PX: float = 0.2125
"""spatialdata_xenium_explorer 硬編碼的像素尺寸，需要修補"""
