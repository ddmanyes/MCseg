import { useQuery } from '@tanstack/react-query'
import { useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'

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

  const query = useQuery({
    queryKey: [stage, 'status'],
    queryFn: async () => (await queryFn()).data,
    refetchInterval: (q) =>
      q.state.data?.status === 'running' ? interval : false,
    staleTime: interval - 100,
  })

  useEffect(() => {
    if (query.data) updateStage(stage, query.data)
  }, [query.data, stage, updateStage])

  return query
}
