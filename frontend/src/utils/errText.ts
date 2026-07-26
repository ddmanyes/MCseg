import axios from 'axios'

/**
 * 把任意 catch 到的東西轉成一行可讀訊息。
 *
 * 優先序：後端 `{"status":"error","message":...}` → HTTP 狀態 → 網路層 → 通用。
 * 後端統一回應格式見 CLAUDE.md §9。
 */
export function errText(e: unknown): string {
  if (axios.isAxiosError(e)) {
    const body = e.response?.data as { message?: string; detail?: string } | undefined
    if (typeof body?.message === 'string' && body.message) return body.message
    if (typeof body?.detail === 'string' && body.detail) return body.detail
    if (e.response) return `HTTP ${e.response.status}`
    // 沒有 response = 請求根本沒送達（後端未啟動／網路中斷）
    return e.code === 'ECONNABORTED' ? '請求逾時' : '無法連線到後端'
  }
  if (e instanceof Error) return e.message
  return String(e)
}

/**
 * 後端「業務層失敗」偵測。
 *
 * ⚠️ 本專案的 API 失敗時回的是 **HTTP 200 + `{"status":"error"}`**（CLAUDE.md §9），
 * 所以 axios 不會 reject，單靠 `.catch()` 抓不到 —— 例如
 * `PUT /segmentation/roi_overrides` 遇到未知 ROI 名稱時就是這種回應。
 * 任何「寫入型」呼叫都必須額外過這一關，否則畫面會顯示成功而後端根本沒存。
 *
 * @returns 失敗訊息；成功時回 null
 */
export function apiErrorMessage(res: { data?: unknown }): string | null {
  const body = res.data as { status?: string; message?: string; detail?: string } | undefined
  if (body?.status !== 'error') return null
  return body.message ?? body.detail ?? '後端回報失敗但未提供訊息'
}
