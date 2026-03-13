import { useState } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { exportXenium, exportLoupe, getXeniumStatus, getLoupeStatus } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

export default function Stage5_Export() {
  useStageLog('export')
  const { stages, updateStage } = usePipelineStore()
  const xenium = stages['xenium']
  const loupe  = stages['loupe']
  const [outDir, setOutDir] = useState('')
  const [inputH5ad, setInputH5ad] = useState('')
  const [maskSource, setMaskSource] = useState('auto')
  const { refetch: refetchXenium } = useStageStatus('xenium', getXeniumStatus, 3000)
  const { refetch: refetchLoupe }  = useStageStatus('loupe',  getLoupeStatus,  3000)

  const handleXenium = async () => {
    updateStage('xenium', { status: 'running', progress: 0, message: '匯出至 Xenium Explorer...' })
    await exportXenium({ output_dir: outDir || undefined, input_h5ad: inputH5ad || undefined, mask_source: maskSource })
    refetchXenium()
  }

  const handleLoupe = async () => {
    updateStage('loupe', { status: 'running', progress: 0, message: '匯出至 Loupe Browser...' })
    await exportLoupe({ output_dir: outDir || undefined, input_h5ad: inputH5ad || undefined, mask_source: maskSource })
    refetchLoupe()
  }

  return (
    <div className="space-y-4">
      <div className="bg-surface-card rounded-xl border border-surface-border p-4">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <div>
            <label className="text-xs text-gray-400">來源 h5ad 檔案名稱（選填）</label>
            <input value={inputH5ad} onChange={e => setInputH5ad(e.target.value)}
                   placeholder="roi/2/proseg_cells.h5ad"
                   className="w-full mt-1 px-3 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none" />
          </div>
          <div>
            <label className="text-xs text-gray-400">輸出目錄（選填）</label>
            <input value={outDir} onChange={e => setOutDir(e.target.value)}
                   placeholder="預設：export_xenium"
                   className="w-full mt-1 px-3 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none" />
          </div>
          <div>
            <label className="text-xs text-gray-400">強制細胞遮罩來源</label>
            <div className="flex gap-2 mt-1">
              {['auto', 'cellpose', 'proseg'].map(src => (
                <button key={src} onClick={() => setMaskSource(src)}
                        className={`flex-1 py-1.5 rounded text-sm transition-colors border ${
                          maskSource === src ? 'bg-primary/20 border-primary text-primary-light' 
                          : 'bg-surface border-surface-border text-gray-400 hover:text-gray-200'
                        }`}>
                  {src.toUpperCase()}
                </button>
              ))}
            </div>
          </div>
        </div>
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
