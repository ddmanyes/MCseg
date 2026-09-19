# MCseg

[English](README.md) | **繁體中文**

### AI 代理引導的工作流程搜尋，用於免寫程式的細胞分割與空間轉錄體轉錄本歸屬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](pyproject.toml)

**MCseg（Multiple Cellpose Segmentation）** 是在本機執行的免寫程式分析平台，將 **Visium HD 的 H&E 影像與 2 µm 空間表現量網格（bins）轉成細胞層級的空間轉錄體資料**。平台以同一個網頁介面串接感興趣區域（ROI）選取、細胞分割、轉錄本歸屬、品質控制、分群、細胞型別註解、空間視覺化與結果匯出；也提供 CLI 進行全切片批次處理。

AI 代理協助方法開發階段的候選流程搜尋。**日常分析在本機執行開發後保留的流程，不需要重新執行 AI 代理搜尋、使用外部語言模型 API，或提供 Xenium 參考資料。** 首次安裝與模型下載需要網路。

[桌面安裝](#桌面安裝windows-與-macos) · [快速開始](#快速開始) · [工作流程](#工作流程) · [研究結果](#研究結果) · [命令列](#命令列操作) · [操作指南（英文）](docs/usage.md) · [重現性](#重現性) · [引用](#引用)

## 工作流程

<p align="center">
  <img src="docs/fig1_development_deployment.png" width="1000" alt="MCseg 開發與部署：研究者定義操作元件，由 AI 代理透過 Xenium 參考評分搜尋並保留分割流程；使用者在本機介面匯入影像、分割、RNA 計數、分析與匯出。右側顯示 H&E、MCseg 分割遮罩與 Xenium 參考邊界。">
</p>

**方法開發與日常部署分為兩個階段。** (a) 研究者定義操作元件、參考資料、目標與執行限制，由代理提出並評估候選流程，再由研究者審查保留的流程。(b) 使用者透過本機平台分析自己的影像與空間表現量資料。(c) 代表性的 H&E 影像、MCseg 分割遮罩與 Xenium 參考邊界。Xenium 邊界是計算得到的參考，不是人工描繪的完整細胞真值。

<details>
<summary>圖片來源與模型名稱</summary>

此圖由文章原始圖稿直接裁切，未重新取樣或更改面板內容。原圖中的「spsam」為 `cpsam` 的拼字誤植，模型說明見[模型實作](#模型實作)。裁切來源與範圍記錄於 [figure-source.md](docs/figure-source.md)。

</details>

| 步驟 | 操作 | 主要產出 |
| --- | --- | --- |
| 設定資料 | 選取 H&E 與對應的 Space Ranger 輸出 | 已驗證的輸入路徑 |
| 選取區域 | 繪製 ROI，或使用全切片 CLI | 影像裁切與空間座標 |
| 分割細胞 | 執行多輪整合分割與受約束擴張 | 細胞分割遮罩（masks） |
| 歸屬表現量 | 將空間 bins 分配至遮罩並加總計數 | 細胞 × 基因 AnnData 矩陣 |
| 分析 | 篩選細胞、計算 PCA／UMAP／Leiden，並以 CellTypist 註解 | 細胞特徵與型別標籤 |
| 探索與匯出 | 檢視空間表現量並匯出結果 | AnnData、Xenium Explorer 或 Loupe Browser 輸出 |

詳細操作請見[介面導覽](docs/usage.md#interface-tour)、[逐步操作](docs/usage.md#usage-guide)與[輸出結構](docs/usage.md#output-structure)（英文）。

## 研究結果

以下數值取自 2026 年 9 月 18 日文章稿件。各項評估的資料與條件不同，不應合併解讀為單一整體效能排名。

| 評估 | 文章結果 | 解讀範圍 |
| --- | --- | --- |
| **LUAD：固定參數，6 個開發 ROI** | PQ **0.472 ± 0.072**；Optuna 調參的 2Cseg 為 **0.432 ± 0.037** | PQ 絕對增加 0.040，約 9% 相對提升；這些 ROI 參與了方法開發 |
| **LUAD：參考引導校準** | PQ **0.554 ± 0.063** | 使用 Xenium 遮罩逐 ROI 選擇擴張策略與距離的上限分析，不代表日常部署效能 |
| **CRC：專家審閱的 ENACT 參考** | MCseg 與 ENACT 的 micro-F1 分別為 **0.805、0.723** | 僅比較兩方法共同涵蓋的 **10,275** 個參考細胞；MCseg 涵蓋全部 20,991 個參考中心點的 65.1% |
| **CRC：轉錄本品質，15 個 ROI** | MCseg 與 Space Ranger 的 NED 為 **0.727 vs 0.712**；互斥細胞譜系共表現率為 **0.49% vs 0.67%** | sign-flip p 分別為 0.008、0.010；UMI 密度接近（**11.6 vs 11.7 UMIs/µm²**），但轉錄本捕獲比例較低（**0.737 vs 0.934**） |
| **新鮮冷凍乳癌：固定流程** | **96,876** 個細胞、FTC **0.514**、每個細胞 UMI 中位數 **1,502**、NED **0.519** | 未重新搜尋架構的跨組織應用，不代表已證實在乳癌資料中優於其他方法 |

**指標解讀。** PQ 同時考量邊界吻合程度與偵測完整性。FTC 表示組織 UMI 中被分配至細胞遮罩的比例。NED 衡量相鄰遮罩間的表現量分離程度，不是絕對幾何準確度。NED 較高或細胞譜系混合較少，均不足以單獨證明分割較好；UMI 密度接近也不等於轉錄本捕獲量相同。

**驗證範圍。** LUAD 幾何結果屬於開發集估計。CRC 專家審閱區域與 15 個 CRC ROI 不重疊，但來自同一組織切片；其參考中心點最初由 StarDist 偵測，再經人工審閱。跨組織驗證仍有限，擴張設定會受組織形態與影像尺度影響。相關資料見[結果範例](docs/usage.md#example-results)與 [analysis](analysis/) 目錄。

## 快速開始

### 輸入資料與系統需求

請準備 H&E 影像及**與影像相符**的 Space Ranger Visium HD 輸出：

- H&E 影像，通常為分塊 BigTIFF（`.btf`、`.tif` 或 `.tiff`）。
- **2 µm bins** 的 `tissue_positions.parquet` 與 `filtered_feature_bc_matrix.h5`。
- 所選流程需要的空間中繼資料。另行掃描的影像必須先完成對位，再進行轉錄本歸屬；見[影像格式與對位](docs/usage.md#supported-image-formats)。

從原始碼安裝需要 **Python ≥3.10**、**uv**，以及網頁介面使用的 **Node.js/npm**；只使用 CLI 則不需要 Node.js。建議準備至少 16 GB RAM，大型影像需要更多記憶體；實際記憶體與磁碟需求取決於資料及安裝的依賴。文章分析使用 Apple Silicon 的 MPS 加速。系統支援 CPU，但未系統性測試 CPU 執行時間；文章中使用 MPS 的幾何評估約需每個 ROI 20–40 分鐘。

### 桌面安裝：Windows 與 macOS

桌面安裝包內含 MCseg 介面、後端程式與 `uv` 環境管理工具，**不需要自行安裝 Python、Node.js、Rust 或 uv**。首次啟動會下載並建立 Python／PyTorch／Cellpose 環境。這不是離線安裝包：請保持網路連線，預留至少 15 GB 可用磁碟空間供初始化使用，並額外保留資料與分析結果的空間。

| 平台 | 桌面安裝包 | 處理器架構 |
| --- | --- | --- |
| Windows 10/11 | [0.2.0 舊版安裝檔](https://github.com/ddmanyes/MCseg/releases/download/v0.8.0/mcseg_0.2.0_x64-setup.exe)；清理後的 0.2.1 尚待建置 | Intel／AMD x64 |
| macOS 12+ | [mcseg_0.2.1_aarch64.dmg（預發行版）](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.1/mcseg_0.2.1_aarch64.dmg) | Apple Silicon（M 系列） |

**發行狀態：**[桌面版 0.2.1](https://github.com/ddmanyes/MCseg/releases/tag/desktop-v0.2.1) 提供清理後的 macOS 安裝包、[SHA-256 校驗檔](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.1/mcseg_0.2.1_aarch64.dmg.sha256)與[建置來源紀錄](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.1/mcseg_0.2.1_aarch64.dmg.build.json)。此安裝包已通過建置、內容檢查、ad-hoc 簽章完整性及靜態憑證掃描；尚未驗證乾淨機器上的首次安裝與完整分析流程。Windows 0.2.1 尚未建置。舊版 Windows 0.2.0 仍含開發用輸入路徑預設值，使用前請改選自己的資料。舊版 macOS 0.2.0 曾包含開發分析狀態，已由此候選版本取代。

桌面版與原始碼／Python 套件（`0.8.0`）使用不同版本編號。GitHub 的「Source code」壓縮檔不是桌面安裝包。此處的 macOS 安裝包僅適用於 Apple Silicon，未提供 Intel Mac 安裝檔。

#### Windows

1. 雙擊 **`mcseg_0.2.0_x64-setup.exe`**，依安裝精靈完成安裝。
2. 若 Microsoft Defender SmartScreen 顯示無法辨識的應用程式，先確認安裝檔來自 MCseg 維護者，再於系統政策允許時選擇 **其他資訊 → 仍要執行**。
3. 從開始功能表啟動 **MCseg**。環境建立與分析引擎啟動期間，請保持初始化視窗開啟。
4. 初始化完成後會開啟主介面。選取自己的資料，再依[操作指南](docs/usage.md#usage-guide)開始分析。

#### macOS（Apple Silicon）

1. 開啟 **`mcseg_0.2.1_aarch64.dmg`**，將 **MCseg** 拖入 **Applications（應用程式）**。
2. 從 Applications 啟動。若 macOS 阻擋未經公證的版本，先確認來源，再於系統提供此選項時使用 **系統設定 → 隱私權與安全性 → 強制打開**，並確認提示。
3. 建立 Python 環境與下載依賴期間，請保持網路連線與初始化視窗開啟。
4. 初始化完成後會開啟主介面；之後啟動會重用已建立的環境。

#### 首次啟動疑難排解

| 狀況 | 檢查方式 |
| --- | --- |
| 初始化看似停住 | 展開初始化紀錄，檢查網路、磁碟空間與套件是否仍在下載。先排除紀錄中的錯誤，再選擇 **Retry**。 |
| 連接埠 8001 已被占用 | 關閉自己啟動的其他 MCseg／後端工作階段；若為其他程式占用，先確認程式身分再處理。 |
| 安裝包不適用於電腦 | 確認 Windows x64 與 Apple Silicon macOS 的差別；其他環境可從原始碼安裝。 |

### 從原始碼安裝

適合開發者、CLI 使用者，或沒有對應桌面安裝包的系統。

安裝 [uv](https://docs.astral.sh/uv/getting-started/installation/) 與 [Node.js](https://nodejs.org/)，將專案複製到原生檔案系統（macOS 使用 APFS、Windows 使用 NTFS，或 Linux 原生檔案系統）：

```bash
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync
npm --prefix frontend install
```

在專案根目錄啟動後端：

```bash
uv run uvicorn backend.main:app --host 127.0.0.1 --port 8001
```

於第二個終端機切換到同一專案根目錄，啟動前端：

```bash
npm --prefix frontend run dev
```

開啟 **[http://localhost:3000](http://localhost:3000)**，選取資料並依[操作指南](docs/usage.md#usage-guide)開始分析。上述指令也適用於 PowerShell。首次下載依賴與模型可能需要較長時間。

<details>
<summary>啟動腳本與外接硬碟</summary>

專案另提供 `start.sh` 與 `start.ps1`。使用前請先檢查內容：目前的 macOS shell 啟動腳本會將專案環境重建為指向 `~/.venvs/msseg` 的符號連結，並終止占用 8001／3000 連接埠的程序。在原生檔案系統上使用上述兩個終端機的啟動方式即可，不需要這些額外處理。

若資料位於 ExFAT 外接硬碟，建議將程式與 Python 環境放在系統磁碟，再讀取外接硬碟的資料。Windows 的 ExFAT 安裝可能需要設定 `UV_LINK_MODE=copy`；詳見[疑難排解](docs/usage.md#troubleshooting)。

</details>

## 命令列操作

執行指令保留既有名稱 **`msseg-segment`**。

```bash
uv run msseg-segment \
  --btf /path/to/image.btf \
  --tp /path/to/tissue_positions.parquet \
  --h5 /path/to/filtered_feature_bc_matrix.h5 \
  --out /path/to/output \
  --tissue crc \
  --cpsam
```

此指令會執行分割、bin 歸屬、細胞矩陣彙整與 CellTypist 註解。加入 `--export-xenium` 可匯出 Explorer 套件，`--skip-celltypist` 可略過註解，`--no-gpu` 則使用 CPU。若只提供 `--btf` 與 `--out`，會僅執行分割。`--cpsam` 開啟額外三輪推論，組成七輪配置；主要與額外輪次均使用 `cpsam`。

```bash
uv run msseg-segment --help
```

詳見[完整選項與 PowerShell 範例](docs/usage.md#cli-no-ui-whole-slide-pipeline)。

## 方法與設定

保留的流程結合 CLAHE 前處理、不同影像表示與直徑設定下的多輪 Cellpose `cpsam` 推論、依優先順序整合遮罩、可選的轉錄本密度補救，以及 Voronoi 約束邊界擴張。應用程式提供較精簡的配置，也能開啟額外輪次，使用七輪流程。

轉錄本歸屬會將 bin 中心點轉換到影像／遮罩座標，再將計數加總為稀疏的細胞 × 基因矩陣。正確的影像對位、像素尺度與 ROI 位移十分重要。RNA 計數階段還可透過 `rna_counting.dilation_px` 額外擴張遮罩；這會改變歸屬範圍，比較結果或重現基準時必須記錄。

設定檔位於 [`config/pipeline.yaml`](config/pipeline.yaml) 與 [`config/profiles/`](config/profiles/)。組織設定檔提供起始值，pipeline 與執行時設定可能覆寫這些數值。每份資料都應檢查最終遮罩與實際生效參數。

### 模型實作

MCseg 的多輪整合分割使用 **Cellpose `cpsam`**。各輪差異在影像表示、直徑與細胞機率閾值（cell-probability threshold），並非混合不同模型家族。七輪配置包含三輪不同直徑、一輪 hematoxylin 輸入，以及三輪額外 `cpsam` 推論。[`_load_primary_model`](backend/src/segmentation/cellpose_runner.py) 載入器會在執行紀錄中寫入實際模型權重路徑。

為便於重現，請保留分析版本、環境、實際模型權重與生效設定。上述研究結果取自文章，本次 README 更新未重新執行科學分析。

## 重現性

| 資源 | 內容 |
| --- | --- |
| [`analysis/scripts/`](analysis/scripts/) | 基準分析與製圖腳本 |
| [`analysis/data/`](analysis/data/) | 已納入版本控制的指標與摘要表格 |
| [`analysis/supplementary/`](analysis/supplementary/) | 補充說明與表格；使用前請核對與稿件的版本一致性 |
| [`docs/autoResearch/`](docs/autoResearch/) | 開發提示詞、執行器與起始範本 |
| [`backend/src/`](backend/src/) | 部署使用的分析實作 |

代理引導的開發迴圈參考 AutoResearch 方法：研究者選定候選操作、參考資料、評分方式、提示詞與執行限制，由代理提出並評估可執行流程，再由研究者審查保留的配置。開發範本與日常分析分開，執行範本需要另行設定 API；提供範本不代表已完整封存每次歷史搜尋紀錄。

重現分析時，請保留 Git 版本、解析後的依賴、模型權重識別資訊、輸入資料與座標、生效的分割／計數參數，以及執行紀錄。專案套件沿用歷史名稱 `msseg`。以目前原始碼或清理後的 macOS 0.2.1 安裝包全新安裝時，預設輸入路徑與 ROI 選取為空，請先設定自己的資料；既有安裝可能保留先前儲存的設定。

### 資料取得

- **LUAD：**10x Genomics dataset portal 的配對 Visium HD 與 Xenium Prime 資料；使用六個開發 ROI。
- **CRC：**[GEO GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318)；包含 15 個轉錄本評估 ROI，以及同一切片上另一個經專家審閱的 ENACT 區域。
- **乳癌：**10x Genomics dataset portal 的公開新鮮冷凍 Visium HD 資料，用於固定流程的跨組織應用。

文章規劃將處理後的 AnnData 與分割遮罩存入 Zenodo。目前此處尚未提供公開典藏識別碼或文章 DOI。

## 引用

若使用 MCseg，請引用以下稿件：

> Chan, C.-R., Chang, N.-W., Wang, C.-Y., Tan, H.-Y., and Lin, S.-J. (2026). **MCseg: AI agent-guided workflow search for no-code cell segmentation and transcript attribution in spatial transcriptomics.** Manuscript.

Chan 與 Chang 為共同第一作者。文章規劃先發表於 bioRxiv，正式上線並取得 DOI 後再更新為預印本引用。

## 支援與授權

操作細節請見[操作指南（英文）](docs/usage.md)。如遇問題，請至 [GitHub Issues](https://github.com/ddmanyes/MCseg/issues) 提供作業系統、Git 版本、套件版本、指令／設定及相關紀錄。

MCseg 採用 [MIT License](LICENSE)。第三方元件資訊請見 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)，桌面建置與驗證方式見[桌面建置指南（英文）](docs/desktop-build.md)。
