import { useQuery } from '@tanstack/react-query'
import { useEffect, useRef } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import { notify } from '../utils/notify'

// 桌面殼縮到系統列後，這是使用者知道「跑完了」的唯一方式（分析階段
// 可能長達數十小時，見 msseg-tauri 打包計畫）。純瀏覽器情境 notify()
// 自己會 no-op，這裡不用額外判斷環境。
const STAGE_LABELS: Record<string, string> = {
  roi: 'ROI 裁切',
  segmentation: 'MCseg 分割',
  count: 'RNA 計數',
  analysis: 'Scanpy 分析（QC/UMAP/CellTypist）',
  spatial: '空間圖表探索',
  xenium: 'Xenium 匯出',
  loupe: 'Loupe 匯出',
}

/**
 * TanStack Query 封裝：替換 setInterval 輪詢模式。
 * - status === 'running' 時每隔 interval ms 重新拉取
 * - 其他狀態停止輪詢（不浪費請求）
 * - 元件 unmount 時自動清除（無洩漏）
 * - 結果同步至 Zustand，讓 Sidebar 狀態點可讀到
 * - 回傳 refetch 以便 handleRun 後立即觸發第一次拉取
 */
export function useStageStatus(
  stage: string,
  queryFn: () => Promise<any>,
  interval = 3000,
) {
  const updateStage = usePipelineStore((s) => s.updateStage)
  const prevStatus = useRef<string | undefined>(undefined)

  const query = useQuery({
    queryKey: [stage, 'status'],
    queryFn: async () => (await queryFn()).data,
    refetchInterval: (q) =>
      q.state.data?.status === 'running' ? interval : false,
    staleTime: interval - 100,
  })

  useEffect(() => {
    if (!query.data) return
    updateStage(stage, query.data)

    // 只在「這次 session 裡真的從 running 轉換過來」才通知——避免頁面
    // 一載入、query 第一次拿到本來就已經是 done/error 的舊狀態，就誤發
    // 一次通知（那並不是「剛跑完」，是「本來就跑完了」）。
    const wasRunning = prevStatus.current === 'running'
    const newStatus = query.data.status
    if (wasRunning && newStatus === 'done') {
      const label = STAGE_LABELS[stage] ?? stage
      void notify('MCseg 完成', `${label} 已完成`)
    } else if (wasRunning && newStatus === 'error') {
      const label = STAGE_LABELS[stage] ?? stage
      void notify('MCseg 發生錯誤', `${label} 執行失敗，請查看記錄`)
    }
    prevStatus.current = newStatus
  }, [query.data, stage, updateStage])

  return query
}
