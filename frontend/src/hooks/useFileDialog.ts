import { isTauri } from '@tauri-apps/api/core'

/**
 * 原生檔案/資料夾選取，僅在 Tauri 桌面殼內可用。
 *
 * 這個前端同時服務兩種情境：純瀏覽器（`npm run dev` / `start.sh`，沒有
 * 任何 Tauri API 可用）與桌面殼（Tauri 視窗載入同一份前端）。用
 * `isTauri()` 偵測執行環境——瀏覽器情境呼叫端會拿到 `null`，維持既有的
 * 自建資料夾瀏覽器（`FolderBrowser`，見 DataSetup.tsx）行為不變；只有
 * 真的在 Tauri 裡才會跳出 Finder / 檔案總管原生對話框。
 *
 * 使用者取消原生對話框時也回傳 `null`——呼叫端不應該把「取消」誤判成
 * 「原生對話框不可用」而又跳回自建瀏覽器，這兩者要分開處理（見
 * DataSetup.tsx 的呼叫端邏輯）。
 */
export function useFileDialog() {
  const pickDirectory = async (defaultPath?: string): Promise<string | null> => {
    if (!isTauri()) return null
    const { open } = await import('@tauri-apps/plugin-dialog')
    const result = await open({ directory: true, multiple: false, defaultPath })
    return typeof result === 'string' ? result : null
  }

  const pickFile = async (
    extensions: string[],
    defaultPath?: string,
  ): Promise<string | null> => {
    if (!isTauri()) return null
    const { open } = await import('@tauri-apps/plugin-dialog')
    const result = await open({
      directory: false,
      multiple: false,
      defaultPath,
      filters: extensions.length ? [{ name: '檔案', extensions }] : undefined,
    })
    return typeof result === 'string' ? result : null
  }

  return { isNative: isTauri(), pickDirectory, pickFile }
}
