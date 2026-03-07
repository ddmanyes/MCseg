import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runProseg, getProsegStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

export default function Stage3_Proseg() {
  useStageLog('proseg')
  const { stages, updateStage, recommendedCondition } = usePipelineStore()
  const stage = stages['proseg']
  const { refetch: refetchStatus } = useStageStatus('proseg', getProsegStatus, 5000)

  const params = recommendedCondition ?? { max_dist: 40, compactness: 0.06, dilation: 20 }

  const handleRun = async () => {
    updateStage('proseg', { status: 'running', progress: 0, message: '啟動 Proseg...' })
    await runProseg({ max_dist: params.max_dist, compactness: params.compactness, dilation: params.dilation })
    refetchStatus()
  }

  return (
    <div className="space-y-4">
      <StageCard title="Proseg RNA 重新分配" status={stage.status}
                 progress={stage.progress} message={stage.message} onRun={handleRun} runLabel="執行 Proseg">
        <div className="text-sm text-gray-400 space-y-1">
          {recommendedCondition && <p className="text-yellow-400 text-xs">使用條件測試推薦參數</p>}
          <p>max_dist: {params.max_dist} µm</p>
          <p>compactness: {params.compactness}</p>
          <p>dilation: {params.dilation} px</p>
        </div>
      </StageCard>
      <Terminal stage="proseg" />
    </div>
  )
}
