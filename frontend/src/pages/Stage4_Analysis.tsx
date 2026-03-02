import { useState } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runAnalysis, getAnalysisStatus, getUmap } from '../api/client'
import useStageLog from '../hooks/useStageLog'

export default function Stage4_Analysis() {
  useStageLog('analysis')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['analysis']
  const [umapImg, setUmapImg] = useState<string | null>(null)

  const handleRun = async () => {
    updateStage('analysis', { status: 'running', progress: 0, message: '執行聚類...' })
    await runAnalysis()
    const poll = setInterval(async () => {
      const s = await getAnalysisStatus()
      updateStage('analysis', s.data)
      if (s.data.status !== 'running') {
        clearInterval(poll)
        getUmap().then(r => { if (r.data.data) setUmapImg(r.data.data.image_b64) })
      }
    }, 5000)
  }

  return (
    <div className="space-y-4">
      <StageCard title="下游聚類分析（Scanpy + Leiden）" status={stage.status}
                 progress={stage.progress} message={stage.message} onRun={handleRun} runLabel="執行分析">
        <div className="text-sm text-gray-400 space-y-1">
          <p>流程：QC → normalize → HVG → PCA → UMAP → Leiden</p>
        </div>
      </StageCard>

      {umapImg && (
        <div className="bg-surface-card rounded-xl border border-surface-border p-4">
          <p className="text-sm font-medium text-gray-300 mb-2">UMAP</p>
          <img src={`data:image/png;base64,${umapImg}`} className="max-w-full rounded" alt="UMAP" />
        </div>
      )}

      <Terminal stage="analysis" />
    </div>
  )
}
