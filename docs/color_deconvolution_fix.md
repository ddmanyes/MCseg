# H&E 色彩解卷積轉置修正（2026-09-20）

- 根因：染色向量按 row 排列，`OD = C @ M`；逆運算應為
  `C = OD @ inv(M)`，原本額外的 `.T` 會混淆染色濃度。
- 修正：僅移除 `color_deconvolution_he` 的逆矩陣轉置並補上方向註解。
- 重現：已知 H/E/residual 濃度經正向 Beer-Lambert 模型合成 RGB，
  分別測試浮點與 uint8 量化輸入。舊版兩項還原測試均失敗，部分非零 H
  被壓成零；修正版均通過。另驗證實際 CLAHE 的輸出尺寸、uint8 與 RGBA 相容性。
- 驗證：下列命令共 112 passed；第一次執行因 Numba 快取目錄不可寫而有
  2 failed / 4 errors，指定可寫暫存快取後全部通過。

```bash
NUMBA_CACHE_DIR=/private/tmp/msseg-numba-cache .venv/bin/python -m pytest \
  backend/tests/test_17_color_deconvolution.py \
  backend/tests/test_09_roi_overrides.py \
  backend/tests/test_10_fullslide.py \
  backend/tests/test_14_cli.py -q
```

本次未重跑 b00、模型推論或歷史分割結果，亦未重新評估 AP。
此修正會改變使用 Hematoxylin 通道的後續輸入；歷史結果不會自動更新。
