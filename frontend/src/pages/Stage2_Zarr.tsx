import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { buildZarr, getZarrStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'

export default function Stage2_Zarr() {
  useStageLog('zarr')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['zarr']

  const handleRun = async () => {
    updateStage('zarr', { status: 'running', progress: 0, message: '建構 Zarr...' })
    await buildZarr()
    const poll = setInterval(async () => {
      const s = await getZarrStatus()
      updateStage('zarr', s.data)
      if (s.data.status !== 'running') clearInterval(poll)
    }, 3000)
  }

  return (
    <div className="space-y-4">
      <StageCard title="Zarr 建構（SpatialData OME-Zarr）" status={stage.status}
                 progress={stage.progress} message={stage.message} onRun={handleRun} runLabel="建構 Zarr">
        <div className="text-sm text-gray-400 space-y-1">
          <p>整合：H&E 影像 + 核遮罩 + 2µm/8µm matrix</p>
          <p>格式：SpatialData OME-Zarr（多尺度）</p>
        </div>
      </StageCard>
      <Terminal stage="zarr" />
    </div>
  )
}
