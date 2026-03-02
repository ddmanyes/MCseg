# visiumHD_pipeline_2 — 前端建置與驗證計畫

> 建立日期：2026-03-03
> 狀態：🔄 進行中

---

## 專案現況

後端模組（FastAPI + 7 個 API 路由 + WebSocket log 串流）已**全部完成**。前端原始碼也已完整撰寫（18 個檔案），結構如下：

| 類別 | 檔案 |
|------|------|
| 入口 | `main.tsx`, `App.tsx`, `index.html`, `index.css` |
| 佈局 | `Sidebar.tsx`, `Header.tsx` |
| 共用元件 | `StageCard.tsx`, `Terminal.tsx` |
| 頁面（7 個） | `Stage0_ROI.tsx` ~ `Stage5_Export.tsx` |
| 狀態管理 | `pipelineStore.ts` (Zustand) |
| API 層 | `client.ts` (Axios) |
| WebSocket | `useStageLog.ts` (custom hook) |
| 型別 | `pipeline.ts` |
| 配置 | `vite.config.ts`, `tailwind.config.js`, `postcss.config.js`, `tsconfig.json` |

---

## 執行清單

### Phase 1：安裝與編譯

- [ ] 在 `frontend/` 目錄執行 `npm install`
- [ ] 執行 `npx tsc --noEmit` 檢查 TypeScript 編譯
- [ ] 修復任何型別錯誤或缺少的依賴

### Phase 2：啟動驗證

- [ ] 執行 `npm run dev` 啟動 Vite dev server（port 3000）
- [ ] 使用瀏覽器訪問 `http://localhost:3000`
- [ ] 驗證各頁面可正常載入，Sidebar 導航運作
- [ ] 深色主題正確渲染

### Phase 3：Git 初始化

- [ ] 在專案根目錄執行 `git init` + 初始提交

---

## 驗證方式

### 自動測試

```bash
# TypeScript 編譯檢查
cd frontend && npx tsc --noEmit

# Vite build 檢查（確認無建構錯誤）
cd frontend && npm run build
```

### 瀏覽器驗證

1. 啟動 `npm run dev`（port 3000）
2. 訪問 `http://localhost:3000`
3. 確認重點：
   - Sidebar 左側欄可見，含 Stage 0 ~ Stage 5
   - 點擊每個 Stage 連結可切換頁面
   - Header 顯示當前 Stage 標題
   - 深色主題正確渲染

---

## 變更紀錄

| 日期 | 變更內容 |
|------|---------|
| 2026-03-03 | 初始建立，開始前端建置 |
