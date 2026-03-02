import { useState } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { exportXenium, exportLoupe, getXeniumStatus, getLoupeStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'

export default function Stage5_Export() {
  useStageLog('export')
  const { stages, updateStage } = usePipelineStore()
  const xenium = stages['xenium']
  const loupe  = stages['loupe']
  const [outDir, setOutDir] = useState('')

  const handleXenium = async () => {
    updateStage('xenium', { status: 'running', progress: 0, message: '匯出至 Xenium Explorer...' })
    await exportXenium({ output_dir: outDir })
    const poll = setInterval(async () => {
      const s = await getXeniumStatus()
      updateStage('xenium', s.data)
      if (s.data.status !== 'running') clearInterval(poll)
    }, 3000)
  }

  const handleLoupe = async () => {
    updateStage('loupe', { status: 'running', progress: 0, message: '匯出至 Loupe Browser...' })
    await exportLoupe({ output_dir: outDir })
    const poll = setInterval(async () => {
      const s = await getLoupeStatus()
      updateStage('loupe', s.data)
      if (s.data.status !== 'running') clearInterval(poll)
    }, 3000)
  }

  return (
    <div className="space-y-4">
      <div className="bg-surface-card rounded-xl border border-surface-border p-4">
        <label className="text-xs text-gray-400">輸出目錄（空白 = config 預設）</label>
        <input value={outDir} onChange={e => setOutDir(e.target.value)}
               placeholder="results/export/..."
               className="w-full mt-1 px-3 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none" />
      </div>

      <StageCard title="Xenium Explorer 格式匯出" status={xenium.status}
                 progress={xenium.progress} message={xenium.message}
                 onRun={handleXenium} runLabel="匯出 Xenium">
        <div className="text-sm text-gray-400 space-y-1">
          <p>輸出：morphology.ome.tif + transcripts.zarr + cell_boundaries</p>
          <p className="text-yellow-400 text-xs">自動修補 experiment.xenium pixel_size Bug</p>
        </div>
      </StageCard>

      <StageCard title="Loupe Browser 格式匯出" status={loupe.status}
                 progress={loupe.progress} message={loupe.message}
                 onRun={handleLoupe} runLabel="匯出 Loupe">
        <div className="text-sm text-gray-400 space-y-1">
          <p>輸出：.cloupe + cell_boundaries.geojson</p>
          <p className="text-yellow-400 text-xs">自動指派 10X 白名單條碼</p>
        </div>
      </StageCard>

      <Terminal stage="export" />
    </div>
  )
}
