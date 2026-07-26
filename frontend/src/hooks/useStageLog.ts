import { useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'

/** 指數退避：1s → 2s → 4s … 上限 30s */
const RECONNECT_BASE_MS = 1000
const RECONNECT_MAX_MS = 30000
/** 連續解析失敗幾次才判定協定對不上（單一畸形訊框吞掉是合理的） */
const PARSE_FAIL_REPORT_AT = 5

/**
 * 訂閱某個 stage 的 log WebSocket。
 *
 * 後端 `--reload` 每次改碼都會重啟，舊版沒有 onclose 處理 —— 連線一斷
 * Terminal 就從此不再有輸出，也不說自己斷了。這裡加上自動重連，
 * 並把連線狀態寫進 Terminal（`[WARNING]` 前綴會顯示為黃色）。
 */
export default function useStageLog(stage: string) {
  const appendLog = usePipelineStore((s) => s.appendLog)

  useEffect(() => {
    let ws: WebSocket | null = null
    let retryTimer: ReturnType<typeof setTimeout> | undefined
    let attempt = 0
    let parseFailures = 0
    let parseFailureReported = false
    let disposed = false          // unmount：本地主動關閉，不該再重連

    const connect = () => {
      if (disposed) return
      const protocol = location.protocol === 'https:' ? 'wss' : 'ws'
      ws = new WebSocket(`${protocol}://${location.host}/ws/log/${stage}`)

      ws.onopen = () => {
        if (attempt > 0) appendLog(stage, '[INFO] 已重新連上後端 log 串流')
        attempt = 0
      }

      ws.onmessage = (e) => {
        try {
          const msg = JSON.parse(e.data as string)
          parseFailures = 0
          if (msg.type === 'log') appendLog(stage, `[${msg.level}] ${msg.message}`)
        } catch {
          // 單一畸形訊框忽略；持續失敗代表協定對不上，必須留痕
          parseFailures += 1
          if (parseFailures >= PARSE_FAIL_REPORT_AT && !parseFailureReported) {
            parseFailureReported = true
            appendLog(stage, `[WARNING] log 串流連續 ${parseFailures} 則訊息無法解析為 JSON，前後端協定可能不一致`)
          }
        }
      }

      // onerror 之後瀏覽器一定會再觸發 onclose，重連只掛在 onclose 避免重複排程
      ws.onclose = () => {
        if (disposed) return
        const delay = Math.min(RECONNECT_BASE_MS * 2 ** attempt, RECONNECT_MAX_MS)
        attempt += 1
        appendLog(stage, `[WARNING] log 串流已中斷，${Math.round(delay / 1000)} 秒後重新連線（第 ${attempt} 次）`)
        retryTimer = setTimeout(connect, delay)
      }
    }

    connect()

    return () => {
      disposed = true
      clearTimeout(retryTimer)   // 未取消的待重連 timer 會造成洩漏／連線風暴
      ws?.close()
    }
  }, [stage, appendLog])
}
