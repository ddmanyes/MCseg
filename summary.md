# Visium HD Pipeline 2 - Stage 1 Segmentation 功能開發與修復進度總結

## 1. Eosin 細胞質遮罩 (Cyto Mask) 邏輯重構與修復

- **問題**：原本使用的 `R - B` (Red 減 Blue) 通道相減法來判斷 Eosin 組織分佈，對於判定組織與空腔背景的準確率極低。且在正式分割 `run_segmentation_rois` 中，Watershed 結果錯誤覆寫了原本 Cellpose 辛苦算出的細胞核/膜成果，導致偵測出的細胞數量從數百顆掉到 77 顆以內。
- **解法**：
    1. **改用「亮度」判斷 (Brightness Method)**：對於影像中所有像素取 RGB 最大值 $max(R, G, B)$，大於 `(255 - Eosin_BG_Threshold)` 的視為純空腔背景，其餘保留為組織細胞質。
    2. **分離最終遮罩**：取消對 `final_masks` (即 `segmentation_masks.npy`) 的覆寫與破壞。Cellpose / Logic A 產出的 `final_masks` 僅負責傳遞高品質細胞邊界，而 Eosin 亮度判斷另外生成二元矩陣並獨立儲存為 `cyto_mask.npy`，專供下游 Proseg 定位使用。

## 2. 預覽介面 (Preview UI) 大幅擴充

- **快速 Patch 預覽整合**：在 512x512 小塊預覽中，統一以 Tab 切換顯示：
  - H&E + 綠色邊界
  - Macenko 前處理影像 (Cellpose 真實輸入)
  - 光流方向圖 (Flows dP)
  - Cyto 遮罩
- **完整分割結果預覽支援多圖**：在正式跑完全圖/全 ROI 的任務後，除了能在「完整分割結果預覽」查看 H&E 疊圖外，現在新增：
  - Eosin 細胞質背景圖 (cyto_mask)
  - 小尺寸光流方向圖 (Flows dP)
- **前後端快取機制**：更新 Backend GET `/preview` 端點，支援自動偵測與讀取 `cyto_mask.npy` 與 `flows_preview.jpg` 轉為 base64 傳回給前端切換。

## 3. 面板參數 UI 體驗優化

- **新增滑鼠懸停提示 (Tooltips)**：對所有的 `NumberInput` 與 `Toggle` 元件加上說明支援。
- **翻譯與定義**：為每一個 segmentation 參數補齊繁體中文的實戰意義與調整建議（例如 Batch Size, Flow Threshold, Eosin BG Threshold, Logic A 雙尺寸設計等）。這在實驗與參數調整的過程中，幫助釐清各數值的增減邏輯。

## 下一步/注意事項

由於移除了破壞性 Watershed 邏輯，目前 `segmentation_masks.npy` 細胞數量恢復正常。可以利用新增的 `Flows` 預覽來判斷 Cellpose 找核的品質，並觀察 `Eosin Cytoplasm Mask` 是否能精準框出深粉紅色組織並濾掉背景。確認滿意參數後即可推進到 Proseg 分子指派階段。
