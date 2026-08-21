mod bootstrap;
mod lifecycle;

use tauri::Manager;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        // 單一實例鎖必須是第一個註冊的 plugin（Tauri 官方要求）。第二次啟動
        // 嘗試時，這個 callback 會在「已經在跑」的那個實例裡被呼叫，而不是
        // 真的另外開一個新行程／新的後端——這是防止共用工作站上兩人各自
        // 開一次 App、各自啟動一份後端互撞 port 8001 的第一道防線。
        .plugin(tauri_plugin_single_instance::init(|app, _argv, _cwd| {
            if let Some(w) = app.get_webview_window("main") {
                let _ = w.unminimize();
                let _ = w.show();
                let _ = w.set_focus();
            }
        }))
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_notification::init())
        .setup(|app| {
            // ⚠️ 原本只在 debug build 註冊 log plugin，release build 完全沒有
            // 任何 log 輸出——這台機器沒有畫面可以看 GUI，release 模式下
            // bootstrap 卡住時完全無從診斷。log plugin 預設就會同時寫
            // stdout 與 app 的 log 目錄檔案，兩種模式都保留，對正式使用者
            // 之後回報問題（附上 log 檔）也有用，不是只為了這次除錯。
            app.handle().plugin(
                tauri_plugin_log::Builder::default()
                    .level(log::LevelFilter::Info)
                    .build(),
            )?;
            // 系統列圖示 + 選單 + 視窗關閉事件（縮到系統列而非真的關閉）。
            // 見 lifecycle.rs 開頭註解說明設計依據。
            lifecycle::setup(app.handle())?;
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![bootstrap::run_bootstrap])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
