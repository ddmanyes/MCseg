//! Phase 2：環境配置器（bootstrap）。
//!
//! 設計核心（見 docs/brainstorming/tauri_desktop_packaging.md Phase 2）：
//! - 不用 Python 腳本做 bootstrap ——避免「用 Python 腳本去裝 Python」的雞生蛋問題。
//!   `uv` 本身是獨立 Rust binary，`uv sync` 會自動裝好符合 pyproject.toml
//!   `requires-python` 的直譯器，不需要目標機器預先有任何 Python。
//! - `uv` 執行檔直接當 Tauri sidecar 包進安裝檔（見 tauri.conf.json
//!   bundle.externalBin），不在執行期下載，拿掉最脆弱的「第一步就要有網路」風險。
//! - 不覆寫 UV_CACHE_DIR：讓 uv 用它自己平台原生的預設快取路徑
//!   （macOS/Linux `~/.cache/uv`、Windows `%LOCALAPPDATA%\uv\cache`），
//!   才是「沿用既有快取、不另開一份」這個決定在跨平台下唯一正確的做法——
//!   寫死 `~/.cache/uv` 在 Windows 上其實是錯的路徑，不是 uv 的原生慣例。
//! - port 8001 佔用防護：啟動前一律先打 /api/health 確認，回應
//!   `service == "mcseg"` 才視為「我方既有行程」直接附加；回應了但不是
//!   我方、或格式對不上，一律回報錯誤、绝不主動砍掉佔用該 port 的行程
//!   （對照：`start.sh` 舊有的「先 kill 佔用 port 的行程」邏輯絕對不能原樣
//!   搬進共用工作站會用的桌面殼，會有砍掉別人正在跑的分析工作的風險）。

use serde::Serialize;
use std::path::{Path, PathBuf};
use std::time::Duration;
use tauri::{AppHandle, Emitter, Manager};
use tauri_plugin_shell::ShellExt;
use tokio::process::Command as TokioCommand;

const BACKEND_PORT: u16 = 8001;
const HEALTH_URL: &str = "http://127.0.0.1:8001/api/health";

#[derive(Serialize, Clone)]
struct StepEvent<'a> {
    step: &'a str,
    status: &'a str, // "active" | "done" | "error"
}

#[derive(Serialize, Clone)]
struct MessageEvent<'a> {
    message: &'a str,
}

#[derive(Serialize, Clone)]
struct LogEvent {
    line: String,
}

fn emit_step(app: &AppHandle, step: &str, status: &str) {
    let _ = app.emit("bootstrap://step", StepEvent { step, status });
}

fn emit_message(app: &AppHandle, message: &str) {
    let _ = app.emit("bootstrap://message", MessageEvent { message });
    log::info!("[bootstrap] {message}");
}

fn emit_log(app: &AppHandle, line: impl Into<String>) {
    let line = line.into();
    log::info!("[bootstrap:uv] {line}");
    let _ = app.emit("bootstrap://log", LogEvent { line });
}

fn emit_error(app: &AppHandle, message: &str) {
    log::error!("[bootstrap] {message}");
    let _ = app.emit("bootstrap://error", MessageEvent { message });
}

/// 健康檢查結果：區分「是我方後端」「port 被別的東西佔用」「port 是空的」三種情況。
/// 三選一都要各自處理，不能把後兩者混為一談——這正是先前審查抓到的
/// 「共用工作站砍掉別人行程」風險的根本防線。
enum HealthState {
    Ours,
    OccupiedByOther,
    Free,
}

async fn probe_health() -> HealthState {
    let client = match reqwest::Client::builder().timeout(Duration::from_millis(800)).build() {
        Ok(c) => c,
        Err(_) => return HealthState::Free,
    };
    match client.get(HEALTH_URL).send().await {
        Ok(resp) if resp.status().is_success() => {
            match resp.json::<serde_json::Value>().await {
                Ok(v) if v.get("service").and_then(|s| s.as_str()) == Some("mcseg") => {
                    HealthState::Ours
                }
                _ => HealthState::OccupiedByOther,
            }
        }
        Ok(_) => HealthState::OccupiedByOther,
        // 連線被拒絕／逾時：port 沒人聽，視為空的。
        Err(_) => HealthState::Free,
    }
}

fn app_support_dir(app: &AppHandle) -> Result<PathBuf, String> {
    app.path()
        .app_data_dir()
        .map(|p| p.join("app"))
        .map_err(|e| format!("無法解析應用程式資料目錄：{e}"))
}

fn resource_app_dir(app: &AppHandle) -> Result<PathBuf, String> {
    app.path()
        .resource_dir()
        .map(|p| p.join("app"))
        .map_err(|e| format!("無法解析內嵌資源目錄：{e}"))
}

/// 把唯讀的 Resources/app（隨安裝檔一起打包）複製到可寫入的
/// Application Support/MCseg/app。只在目標不存在，或版本標記檔對不上
/// 目前 App 版本時才複製——重複啟動不必每次都重新複製整包原始碼。
fn sync_resources(app: &AppHandle) -> Result<PathBuf, String> {
    let dst = app_support_dir(app)?;
    let src = resource_app_dir(app)?;
    log::info!("[bootstrap] app_support_dir = {}", dst.display());
    log::info!("[bootstrap] resource_app_dir = {} (exists={})", src.display(), src.exists());
    let version = app.package_info().version.to_string();
    let marker = dst.join(".resource_version");

    let up_to_date = marker
        .exists()
        .then(|| std::fs::read_to_string(&marker).ok())
        .flatten()
        .map(|v| v.trim() == version)
        .unwrap_or(false);

    if up_to_date {
        return Ok(dst);
    }

    std::fs::create_dir_all(&dst).map_err(|e| format!("建立應用程式資料目錄失敗：{e}"))?;

    let mut opts = fs_extra::dir::CopyOptions::new();
    opts.overwrite = true;
    opts.content_only = true; // 複製 src 底下的內容到 dst，不是複製 src 這個資料夾本身

    for entry in ["backend", "config", "frontend"] {
        let s = src.join(entry);
        if s.exists() {
            log::info!("[bootstrap] 複製 {entry} ...");
            let d = dst.join(entry);
            std::fs::create_dir_all(&d).map_err(|e| format!("建立 {entry} 目錄失敗：{e}"))?;
            fs_extra::dir::copy(&s, &d, &opts).map_err(|e| format!("複製 {entry} 失敗：{e}"))?;
            log::info!("[bootstrap] 複製 {entry} 完成");
        } else {
            log::warn!("[bootstrap] 資源來源不存在，跳過：{}", s.display());
        }
    }
    for file in ["pyproject.toml", "uv.lock"] {
        let s = src.join(file);
        if s.exists() {
            std::fs::copy(&s, dst.join(file)).map_err(|e| format!("複製 {file} 失敗：{e}"))?;
            log::info!("[bootstrap] 複製 {file} 完成");
        } else {
            log::warn!("[bootstrap] 資源來源不存在，跳過：{}", s.display());
        }
    }

    std::fs::write(&marker, &version).map_err(|e| format!("寫入版本標記失敗：{e}"))?;
    Ok(dst)
}

fn venv_python_dir(app_dir: &Path) -> PathBuf {
    if cfg!(target_os = "windows") {
        app_dir.join(".venv").join("Scripts")
    } else {
        app_dir.join(".venv").join("bin")
    }
}

fn uvicorn_bin(app_dir: &Path) -> PathBuf {
    let bin = venv_python_dir(app_dir).join("uvicorn");
    if cfg!(target_os = "windows") {
        bin.with_extension("exe")
    } else {
        bin
    }
}

/// 跑 `uv sync --project <app_dir> --no-install-project`，逐行把
/// stdout/stderr 轉成 bootstrap://log 事件。用 sidecar（隨裝隨附的 uv
/// 執行檔），不倚賴系統 PATH 上有沒有裝 uv——這正是拿掉「首次啟動就要
/// 下載 uv」風險的落地方式。
///
/// ⚠️ `--no-install-project` 是實機測試才踩到、必須加的旗標：`uv sync`
/// 預設會把 `msseg` 這個專案本身也建成 editable 安裝，而 hatchling 的
/// metadata 驗證要求 `pyproject.toml` 的 `readme = "README.md"` 指到的
/// 檔案存在——但打包資源沒有複製 README.md，會直接讓 `uv sync` 失敗
/// （`OSError: Readme file does not exist: README.md`）。修法不是補
/// README.md 進資源清單，而是根本不需要把 `msseg` 裝成套件：
/// `backend.main:app` 本來就是用複製過去的原始碼＋`current_dir` 直接
/// 跑，只需要它的**依賴**裝好，不需要 `msseg` 本身可 import。這樣也
/// 避免以後 pyproject.toml 的打包 metadata 需求（README/LICENSE/…）
/// 又跟資源複製清單各自漂移出新的落差。
async fn run_uv_sync(app: &AppHandle, app_dir: &Path) -> Result<(), String> {
    let sidecar = app
        .shell()
        .sidecar("uv")
        .map_err(|e| format!("找不到內嵌的 uv 執行檔：{e}"))?
        .args(["sync", "--project", &app_dir.to_string_lossy(), "--no-install-project"])
        // 刻意不設 UV_CACHE_DIR，讓 uv 用該平台原生預設快取路徑。
        .current_dir(app_dir);

    let (mut rx, _child) = sidecar.spawn().map_err(|e| format!("啟動 uv sync 失敗：{e}"))?;

    use tauri_plugin_shell::process::CommandEvent;
    while let Some(event) = rx.recv().await {
        match event {
            CommandEvent::Stdout(line) | CommandEvent::Stderr(line) => {
                emit_log(app, String::from_utf8_lossy(&line).trim_end().to_string());
            }
            CommandEvent::Error(e) => {
                return Err(format!("uv sync 執行期錯誤：{e}"));
            }
            CommandEvent::Terminated(payload) => {
                return match payload.code {
                    Some(0) => Ok(()),
                    code => Err(format!("uv sync 失敗（exit code {code:?}），詳見上方日誌")),
                };
            }
            _ => {}
        }
    }
    Ok(())
}

/// 啟動後端（不含 --reload，正式版不需要檔案監聽）。以子行程方式常駐——
/// 只要 Tauri App 行程還活著（含縮到系統列的狀態），後端就跟著活著。
///
/// `kill_on_drop(false)` 讓這個 handle 被 drop 時不會誤殺子行程（縮到
/// 系統列時我們就是要它繼續活著）；但這也代表「真的要結束 App」時需要
/// 另外主動終止它——handle 存進 `lifecycle::AppState`，由 Phase 4 的
/// 結束流程（`lifecycle::shutdown_and_exit`）取出來 `.kill()`。單純
/// `kill_on_drop(false)` 而不存 handle，會變成永遠殺不掉、只能手動
/// `pkill` 的孤兒行程（Phase 2 實機測試時就踩過一次）。
async fn spawn_backend(app: &AppHandle, app_dir: &Path) -> Result<(), String> {
    let uvicorn = uvicorn_bin(app_dir);
    if !uvicorn.exists() {
        return Err(format!("找不到 uvicorn 執行檔：{}", uvicorn.display()));
    }
    let child = TokioCommand::new(uvicorn)
        .args(["backend.main:app", "--host", "127.0.0.1", "--port", &BACKEND_PORT.to_string()])
        .current_dir(app_dir)
        .kill_on_drop(false)
        .spawn()
        .map_err(|e| format!("啟動後端行程失敗：{e}"))?;
    crate::lifecycle::store_backend_child(app, child).await;
    Ok(())
}

async fn wait_for_health(timeout: Duration) -> bool {
    let start = std::time::Instant::now();
    while start.elapsed() < timeout {
        if let HealthState::Ours = probe_health().await {
            return true;
        }
        tokio::time::sleep(Duration::from_millis(300)).await;
    }
    false
}

/// Setup Wizard 的唯一入口。boot-ui 載入後立刻呼叫一次；使用者按「重試」
/// 也是呼叫這個。每個階段都會 emit 對應的 step/message/log 事件，
/// 最終成功會 emit `bootstrap://ready`，boot-ui 收到後把視窗導向後端。
#[tauri::command]
pub async fn run_bootstrap(app: AppHandle) -> Result<(), String> {
    emit_step(&app, "check_backend", "active");
    match probe_health().await {
        HealthState::Ours => {
            emit_step(&app, "check_backend", "done");
            emit_step(&app, "python_env", "done");
            emit_step(&app, "start_engine", "done");
            emit_message(&app, "偵測到既有的分析引擎正在執行，直接接上。");
            let _ = app.emit("bootstrap://ready", ());
            return Ok(());
        }
        HealthState::OccupiedByOther => {
            let msg = format!(
                "Port {BACKEND_PORT} 已被其他程式佔用（回應內容不是 MCseg）。\
                 請先關閉佔用該 port 的程式後再重試——為了不誤終止別人可能正在\
                 執行的工作，這裡不會自動嘗試關閉它。"
            );
            emit_step(&app, "check_backend", "error");
            emit_error(&app, &msg);
            return Err(msg);
        }
        HealthState::Free => {
            emit_step(&app, "check_backend", "done");
        }
    }

    emit_step(&app, "python_env", "active");
    emit_message(&app, "正在準備 Python 與分析環境（首次啟動較久）...");

    let app_dir = match sync_resources(&app) {
        Ok(d) => d,
        Err(e) => {
            // ⚠️ 這裡曾經是裸 `?`，會讓錯誤原路傳回 invoke() 的 Promise
            // reject（boot-ui 的 .catch 仍抓得到），但完全不會走
            // emit_error()——log 檔案看不到任何線索,只有真的盯著 GUI
            // 視窗才看得到「重試」按鈕跳出來。這個 bug 是這輪實機測試
            // 才踩到的：本機沒有可視化畫面可看，只能靠 log 檔案判斷,
            // 裸 `?` 直接讓這條路徑變成對我來說的「靜默卡住」。
            emit_step(&app, "python_env", "error");
            emit_error(&app, &e);
            return Err(e);
        }
    };

    if let Err(e) = run_uv_sync(&app, &app_dir).await {
        emit_step(&app, "python_env", "error");
        emit_error(&app, &e);
        return Err(e);
    }
    emit_step(&app, "python_env", "done");

    emit_step(&app, "start_engine", "active");
    emit_message(&app, "正在啟動 MCseg 運算引擎...");
    if let Err(e) = spawn_backend(&app, &app_dir).await {
        emit_step(&app, "start_engine", "error");
        emit_error(&app, &e);
        return Err(e);
    }

    if wait_for_health(Duration::from_secs(30)).await {
        emit_step(&app, "start_engine", "done");
        let _ = app.emit("bootstrap://ready", ());
        Ok(())
    } else {
        let msg = "MCseg 運算引擎啟動逾時（30 秒內未回應健康檢查）".to_string();
        emit_step(&app, "start_engine", "error");
        emit_error(&app, &msg);
        Err(msg)
    }
}
