import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { buildZarr, getZarrStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

export default function Stage2_Zarr() {
  useStageLog('zarr')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['zarr']
  const { refetch: refetchStatus } = useStageStatus('zarr', getZarrStatus, 3000)

  const handleRun = async () => {
    updateStage('zarr', { status: 'running', progress: 0, message: '建構 Zarr...' })
    await buildZarr()
    refetchStatus()
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
