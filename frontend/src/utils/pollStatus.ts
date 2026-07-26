import { errText } from './errText'

export interface PollStatus {
  status: string
  progress?: number
  message?: string
}

/** 連續失敗幾次才判定失聯（× 2s 間隔 ≈ 10 秒） */
export const POLL_MAX_FAILURES = 5

export interface StatusPollOptions {
  /** 取回目前狀態；回傳 null/undefined 視為「本輪沒有新狀態」，不算失敗 */
  fetchStatus: () => Promise<PollStatus | null | undefined>
  /** 拿到狀態時呼叫。status !== 'running' 會自動停止輪詢 */
  onStatus: (d: PollStatus) => void
  /** 連續失敗達上限、輪詢已停止時呼叫 */
  onLost: (message: string) => void
  /** 單次失敗（尚未達上限）時呼叫，用來顯示「重連中 n/N」 */
  onTransientFailure?: (consecutive: number, max: number) => void
  intervalMs?: number
  maxFailures?: number
}

/**
 * 帶重試上限的狀態輪詢。
 *
 * 取代 `catch { clearInterval(...) }` —— 一次網路抖動就永久停止輪詢，
 * 進度條會凍在最後一個數字，使用者分不清「卡住」還是「跑完」。
 *
 * 回傳 interval handle，呼叫端仍需在 unmount 時 `clearInterval`。
 */
export function startStatusPoll(opts: StatusPollOptions): ReturnType<typeof setInterval> {
  const {
    fetchStatus, onStatus, onLost, onTransientFailure,
    intervalMs = 2000, maxFailures = POLL_MAX_FAILURES,
  } = opts

  let failures = 0
  // 用容器持有 handle：tick 在 setInterval 回傳前就被定義，需延後取值
  const timer: { id?: ReturnType<typeof setInterval> } = {}

  const tick = async () => {
    try {
      const d = await fetchStatus()
      failures = 0            // 成功即重置，只有「連續」失敗才算失聯
      if (!d) return
      onStatus(d)
      if (d.status !== 'running') clearInterval(timer.id)
    } catch (e: unknown) {
      failures += 1
      if (failures >= maxFailures) {
        clearInterval(timer.id)
        onLost(
          `已失去與後端的連線（連續 ${failures} 次查詢失敗：${errText(e)}）。` +
          `後端任務可能仍在執行 —— 請重新整理頁面確認狀態。`
        )
      } else {
        onTransientFailure?.(failures, maxFailures)
      }
    }
  }

  timer.id = setInterval(() => void tick(), intervalMs)
  return timer.id
}
