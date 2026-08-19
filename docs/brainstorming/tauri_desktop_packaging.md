# MCseg 桌面應用程式 (Tauri Desktop App) 完整設計規格與實作計畫

- **狀態**：計畫已就緒，等待下次啟動實作 (Plan Ready, Saved for Next Session)
- **記錄日期**：2026-08-19
- **核心目標**：將現有 MCseg 空間轉錄體分析流水線封裝為跨平台原生桌面應用程式（`.dmg` / `.msi`），為實驗室研究人員與生信分析師提供開箱即用、免終端機指令、具備原生檔案對話框與系統通知的一體化體驗。

---

## 🏛️ 系統架構圖 (System Architecture)

```mermaid
graph TD
    subgraph "MCseg 桌面應用程式 (Tauri Desktop App)"
        A["原生視窗 (macOS / Windows)"] <-->|Tauri IPC Event| B["Tauri Rust 核心"]
        B -->|展示| C["MCseg React 前端 (UI)"]
        B -->|原生檔案對話框| F["系統 Finder / 檔案總管"]
        B -->|系統原生通知| G["macOS / Windows 橫幅通知"]
        B -->|Sidecar 生命週期管理| D["微型 uv 隔離環境 (AppData)"]
        D -->|啟動 & 監聽| E["Python FastAPI 後端<br>(PyTorch MPS/CUDA + Cellpose + Scanpy)"]
    end
    C <-->|HTTP / WebSocket (Port 8001)| E
```

---

## 📋 5 大定案規格 (Finalized Specifications)

1. **執行環境分發模式 (Runtime Distribution)**
   * **輕量外殼發布**：安裝包（`.dmg` / `.msi`）維持在 **~20 MB**，不預先打包數 GB 的 PyTorch 與 CUDA 二進位檔。
   * **首次自動環境配置**：
     * 首次啟動時在使用者本機應用程式支援目錄（`~/Library/Application Support/MCseg/` 或 `%APPDATA%/MCseg/`）建立專屬虛擬環境。
     * 利用極速套件管理器 `uv` 自動下載與安裝 PyTorch（依硬體自動適配 MPS/CUDA）、Cellpose 與 Scanpy。
     * 支援軟體快速版本升級，不需每次重新下載數 GB 檔案。

2. **目標作業系統與硬體加速策略 (OS & Hardware Acceleration)**
   * **支援平台**：**macOS (Apple Silicon M 系列)** 與 **Windows 10/11 (NVIDIA CUDA / CPU)** 雙平台。
   * **硬體自動適配**：
     * **macOS**：自動偵測 Apple Silicon 並啟用 `mps` Metal 加速。
     * **Windows**：自動檢測 NVIDIA GPU 與驅動，自動安裝對應 CUDA 版本的 PyTorch；無專用 GPU 則自動回退為多執行緒 CPU 模式。

3. **首次啟動與安裝導覽 (First-Time Setup Wizard)**
   * **視覺化安裝精靈 (Setup Wizard)**：
     * 首次開啟時呈現歡迎與初始化引導卡片。
     * 步驟視覺化：`[✓] 檢查系統硬體與 GPU` $\rightarrow$ `[🔄] 初始化 Python 與生信環境` $\rightarrow$ `[ ] 啟動 MCseg 運算引擎`。
     * 具備動態進度百分比條與可展開的即時日誌抽屜（Log Drawer），遇到網路超時可一鍵重試。
     * 完成後自動平滑過渡進入 Stage 0 數據載入頁面。

4. **本機檔案瀏覽與路徑存取 (Native File Dialog Integration)**
   * **雙模混合輸入**：
     * 點擊各輸入框旁的「📁 瀏覽 (Browse)」按鈕直接呼叫 macOS Finder / Windows 檔案總管原生對話框。
     * 檔案選取後自動將本機絕對路徑（如 `/Volumes/SSD/sample_crc/tissue_hires_image.png`）回填入輸入框。
     * 保留輸入框文字可編輯性，支援手動貼上路徑。

5. **背景行程生命週期與防呆管理 (Lifecycle & Resource Cleanup)**
   * **運算中防呆退出**：當後台有耗時分析正在運行（如 Cellpose 分割、UMAP 計算），關閉視窗時跳出警告對話框確認是否強制終止。
   * **徹底釋放顯存**：確認退出時，Tauri 自動發送 SIGTERM/SIGKILL 終止 Python 後端行程，釋放 Port 8001、系統 RAM 與 GPU 顯存，不留殭屍行程。
   * **重啟防衝突**：啟動前自動檢查並清理殘留的 8001 端口。

---

## 🛠️ 下次啟動之實作執行計畫 (Actionable Implementation Plan)

### 階段一：Tauri 專案骨架與環境初始化 (Phase 1: Tauri Infrastructure)
- [ ] 安裝 `@tauri-apps/cli` 與 Rust Tauri 工具鏈。
- [ ] 執行 `tauri init` 初始化 `src-tauri/` 配置目錄。
- [ ] 配置 `tauri.conf.json`：
  - App Identifier: `com.mcseg.app`
  - 視窗預設尺寸：`1440x920`（最小限制 `1024x720`）
  - 前端建置命令：`npm run build`，靜態產物目錄：`../dist`
  - 啟用安全權限與必要外掛（`dialog`, `notification`, `shell`, `process`）。
- [ ] 設計與配置 MCseg 原生應用程式圖示（App Icons: `.icns` for Mac, `.ico` for Windows）。

### 階段二：微型自動環境配置器與 Setup Wizard (Phase 2: Bootstrapper & Setup UI)
- [ ] 編寫跨平台微型環境安裝器 `scripts/bootstrap_env.py`：
  - 自動下載 Standalone `uv` 執行檔。
  - 執行 `uv venv` 建立專屬虛擬環境。
  - 偵測硬體平台並執行 `uv pip install` 安裝 FastAPI, PyTorch (MPS/CUDA), Cellpose, Scanpy, Anndata, Libvips, OpenCV。
- [ ] 於前端建立 `frontend/src/components/setup/SetupWizard.tsx`：
  - 接收後台環境建立進度（百分比、當前步驟、即時日誌串流）。
  - 當環境檢驗合格時，自動切換至 Stage 0。

### 階段三：原生對話框與系統通知整合 (Phase 3: Native Dialogs & Notifications)
- [ ] 安裝與配置 `@tauri-apps/plugin-dialog`：
  - 封裝通用 Hook `useFileDialog()`。
  - 升級 `Stage0_DataSetup.tsx`、`Stage0_ROI.tsx` 與自訂 Marker CSV 選取器，支援原生 Finder / 檔案總管選取。
- [ ] 安裝與配置 `@tauri-apps/plugin-notification`：
  - 當 Cellpose 分割、QC 評估、UMAP 計算、CellTypist 標註完成時，調用系統橫幅通知。

### 階段四：Rust 生命週期調度與進程守護 (Phase 4: Sidecar Process Lifecycle)
- [ ] 在 `src-tauri/src/main.rs` 實作：
  - 啟動時自動呼叫後端 Uvicorn 服務。
  - 監聽 `tauri::WindowEvent::CloseRequested` 事件。
  - 檢查當前 API 分析狀態；若處於 `running` 狀態則跳出原生確認對話框。
  - 確認關閉時，發送信號安全終止 Python 行程並釋放 Port 8001 與顯存。

### 階段五：編譯打包與驗證 (Phase 5: Build & Packaging)
- [ ] 執行 `npm run tauri build` 編譯 macOS Universal / Apple Silicon `.dmg` 安裝包。
- [ ] 驗證安裝流程（安裝精靈 $\rightarrow$ 數據載入 $\rightarrow$ 空間分割 $\rightarrow$ 下游分析 $\rightarrow$ 關閉退出）。
