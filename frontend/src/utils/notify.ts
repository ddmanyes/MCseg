import { isTauri } from '@tauri-apps/api/core'

let permissionChecked = false
let permissionGranted = false

async function ensurePermission(): Promise<boolean> {
  if (permissionChecked) return permissionGranted
  const { isPermissionGranted, requestPermission } = await import(
    '@tauri-apps/plugin-notification'
  )
  permissionGranted = await isPermissionGranted()
  if (!permissionGranted) {
    const result = await requestPermission()
    permissionGranted = result === 'granted'
  }
  permissionChecked = true
  return permissionGranted
}

/**
 * 系統通知，僅在 Tauri 桌面殼內生效。純瀏覽器情境（`npm run dev` /
 * `start.sh`）直接 no-op——分析階段耗時可能長達數十小時，桌面殼縮到
 * 系統列後靠這個通知使用者「跑完了」，但瀏覽器分頁本來就有自己的分頁
 * 標題/favicon 提示機制，不需要在這裡另外處理。
 *
 * 失敗（例如使用者拒絕通知權限）刻意吞掉，不拋錯——通知本來就是錦上
 * 添花，不該讓分析流程因為通知發送失敗而中斷或彈錯誤訊息。
 */
export async function notify(title: string, body: string): Promise<void> {
  if (!isTauri()) return
  try {
    const granted = await ensurePermission()
    if (!granted) return
    const { sendNotification } = await import('@tauri-apps/plugin-notification')
    sendNotification({ title, body })
  } catch {
    // 通知失敗不影響分析流程本身，刻意靜默。
  }
}
