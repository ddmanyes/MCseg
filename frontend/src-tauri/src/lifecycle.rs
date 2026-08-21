//! Phase 4：視窗生命週期、系統列、結束前的運算中防呆確認。
//!
//! 設計依據（拍板決定，見 docs/brainstorming/tauri_desktop_packaging.md）：
//! - 長任務背景執行 UX = 縮到系統列繼續跑（不是與 GUI 完全解耦的
//!   job-queue 模式）。所以視窗的「關閉」按鈕一律只是隱藏視窗，後端
//!   行程不受影響；只有系統列選單的「結束 MCseg」才是真的退出。
//! - 「結束」前必須確認目前有沒有工作在跑（呼叫 Phase 4 新增的
//!   `GET /api/system/busy` 彙總端點），運算中才跳確認對話框，避免
//!   使用者手滑結束掉正在跑的長任務（可能長達數十小時）。

use std::sync::Arc;
use tauri::menu::{Menu, MenuItem, PredefinedMenuItem};
use tauri::tray::TrayIconBuilder;
use tauri::{AppHandle, Manager, WindowEvent};
use tauri_plugin_dialog::{DialogExt, MessageDialogButtons, MessageDialogKind};
use tokio::process::Child;
use tokio::sync::Mutex;

/// App 層級的共享狀態。目前只放後端子行程的 handle——`bootstrap.rs` 用
/// `kill_on_drop(false)` 讓後端在「縮到系統列」時不會被誤殺，但這也代表
/// 「真的要結束」時沒有一個現成的 handle 可以主動終止它，必須自己存著。
#[derive(Default)]
pub struct AppState {
    pub backend_child: Mutex<Option<Child>>,
}

/// 供 bootstrap.rs 在成功啟動後端後呼叫，把 Child handle 存進 managed state。
pub async fn store_backend_child(app: &AppHandle, child: Child) {
    let state = app.state::<Arc<AppState>>();
    *state.backend_child.lock().await = Some(child);
}

async fn is_backend_busy() -> bool {
    let client = match reqwest::Client::builder()
        .timeout(std::time::Duration::from_secs(3))
        .build()
    {
        Ok(c) => c,
        Err(_) => return false,
    };
    match client.get("http://127.0.0.1:8001/api/system/busy").send().await {
        Ok(resp) => resp
            .json::<serde_json::Value>()
            .await
            .ok()
            .and_then(|v| v.get("busy").and_then(|b| b.as_bool()))
            .unwrap_or(false),
        // 連健康檢查都打不到，代表後端根本沒在跑，談不上「busy」。
        Err(_) => false,
    }
}

/// 真正結束整個 App：終止後端子行程，然後退出。
async fn shutdown_and_exit(app: AppHandle) {
    let state = app.state::<Arc<AppState>>();
    let mut guard = state.backend_child.lock().await;
    if let Some(mut child) = guard.take() {
        log::info!("[lifecycle] 結束前終止後端行程...");
        if let Err(e) = child.kill().await {
            log::warn!("[lifecycle] 終止後端行程失敗（可能已經結束）：{e}");
        }
    }
    drop(guard);
    app.exit(0);
}

/// 系統列選單「結束 MCseg」的處理：先問後端忙不忙，忙的話跳原生確認
/// 對話框；不忙或使用者確認結束，才真的終止後端＋退出。
fn handle_quit_request(app: AppHandle) {
    tauri::async_runtime::spawn(async move {
        if is_backend_busy().await {
            let app_for_dialog = app.clone();
            app.dialog()
                .message("MCseg 目前有分析工作正在執行中（可能耗時數小時），確定要結束嗎？結束後這個工作會被中止且無法恢復。")
                .title("確定要結束 MCseg？")
                .kind(MessageDialogKind::Warning)
                .buttons(MessageDialogButtons::OkCancelCustom(
                    "結束並中止工作".into(),
                    "取消".into(),
                ))
                .show(move |confirmed| {
                    if confirmed {
                        tauri::async_runtime::spawn(shutdown_and_exit(app_for_dialog));
                    }
                });
        } else {
            shutdown_and_exit(app).await;
        }
    });
}

/// 建立系統列圖示與選單（顯示視窗 / 結束）。視窗關閉事件也在這裡一併
/// 註冊：關閉按鈕一律「縮到系統列」（隱藏視窗，行程不受影響），不會
/// 觸發上面的結束確認流程——那是系統列選單「結束」專屬的路徑。
pub fn setup(app: &AppHandle) -> tauri::Result<()> {
    app.manage(Arc::new(AppState::default()));

    let show_item = MenuItem::with_id(app, "show", "顯示視窗", true, None::<&str>)?;
    let quit_item = MenuItem::with_id(app, "quit", "結束 MCseg", true, None::<&str>)?;
    let separator = PredefinedMenuItem::separator(app)?;
    let menu = Menu::with_items(app, &[&show_item, &separator, &quit_item])?;

    TrayIconBuilder::new()
        .icon(app.default_window_icon().cloned().ok_or(tauri::Error::InvalidIcon(
            std::io::Error::new(std::io::ErrorKind::NotFound, "no default window icon"),
        ))?)
        .menu(&menu)
        .tooltip("MCseg")
        .on_menu_event(|app, event| match event.id().as_ref() {
            "show" => {
                if let Some(w) = app.get_webview_window("main") {
                    let _ = w.unminimize();
                    let _ = w.show();
                    let _ = w.set_focus();
                }
            }
            "quit" => handle_quit_request(app.clone()),
            _ => {}
        })
        .build(app)?;

    if let Some(window) = app.get_webview_window("main") {
        let app_for_close = app.clone();
        window.on_window_event(move |event| {
            if let WindowEvent::CloseRequested { api, .. } = event {
                // 一律縮到系統列，不直接關閉——長任務可能還在背景跑。
                // 真正的結束只走系統列選單的「結束 MCseg」。
                api.prevent_close();
                if let Some(w) = app_for_close.get_webview_window("main") {
                    let _ = w.hide();
                }
            }
        });
    }

    Ok(())
}
