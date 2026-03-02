import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runSegmentation, getSegmentationStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'

export default function Stage1_Segmentation() {
  useStageLog('segmentation')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['segmentation']

  const handleRun = async () => {
    updateStage('segmentation', { status: 'running', progress: 0, message: '啟動 Cellpose...' })
    await runSegmentation()
    const poll = setInterval(async () => {
      const s = await getSegmentationStatus()
      updateStage('segmentation', s.data)
      if (s.data.status !== 'running') clearInterval(poll)
    }, 3000)
  }

  return (
    <div className="space-y-4">
      <StageCard title="細胞分割（Cellpose + Logic A）" status={stage.status}
                 progress={stage.progress} message={stage.message} onRun={handleRun} runLabel="執行分割">
        <div className="text-sm text-gray-400 space-y-1">
          <p>模型：nuclei（雙尺寸合併 Logic A）</p>
          <p>Macenko 色彩標準化：啟用</p>
          <p>Eosin Watershed 擴張：啟用</p>
        </div>
      </StageCard>
      <Terminal stage="segmentation" />
    </div>
  )
}
