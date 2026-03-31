import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import {
  exportXenium, exportLoupe,
  getXeniumStatus, getLoupeStatus,
  generateResult, getResultStatus, getResultImages,
} from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'
import { useQuery } from '@tanstack/react-query'

type TabId = 'spatial' | 'umap' | 'dotplot' | 'heatmap'

const TABS: { id: TabId; label: string }[] = [
  { id: 'spatial',  label: '空間分型圖' },
  { id: 'umap',     label: 'UMAP' },
  { id: 'dotplot',  label: 'Dotplot' },
  { id: 'heatmap',  label: 'Heatmap' },
]

export default function Stage4_Export() {
  useStageLog('export')
  const { stages, updateStage } = usePipelineStore()
  const xenium = stages['xenium']
  const loupe  = stages['loupe']
  const [outDir, setOutDir] = useState('')
  const [inputH5ad, setInputH5ad] = useState('')
  const { refetch: refetchXenium } = useStageStatus('xenium', getXeniumStatus, 3000)
  const { refetch: refetchLoupe }  = useStageStatus('loupe',  getLoupeStatus,  3000)

  const [activeTab, setActiveTab] = useState<TabId>('spatial')
  const [resultRunning, setResultRunning] = useState(false)
  const [resultMessage, setResultMessage] = useState('')
  const [activeResultRoi, setActiveResultRoi] = useState<string>('')

  const { data: resultImagesData, refetch: refetchResultImages } = useQuery({
    queryKey: ['result_images'],
    queryFn: () => getResultImages().then(r => r.data),
    staleTime: 0,
  })

  const { data: resultStatusData, refetch: refetchResultStatus } = useQuery({
    queryKey: ['result_status'],
    queryFn: () => getResultStatus().then(r => r.data),
    refetchInterval: resultRunning ? 2000 : false,
  })

  useEffect(() => {
    if (!resultStatusData) return
    const s = resultStatusData as { status: string; message?: string }
    if (s.status === 'done') {
      setResultRunning(false)
      setResultMessage('結果圖已生成')
      refetchResultImages()
    } else if (s.status === 'error') {
      setResultRunning(false)
      setResultMessage(s.message ?? '生成失敗')
    } else if (s.status === 'running') {
      setResultMessage(s.message ?? '生成中...')
    }
  }, [resultStatusData])

  const handleGenerateResult = async () => {
    setResultRunning(true)
    setResultMessage('生成結果圖...')
    await generateResult()
    refetchResultStatus()
  }

  const handleXenium = async () => {
    updateStage('xenium', { status: 'running', progress: 0, message: '匯出至 Xenium Explorer...' })
    await exportXenium({ output_dir: outDir || undefined, input_h5ad: inputH5ad || undefined })
    refetchXenium()
  }

  const handleLoupe = async () => {
    updateStage('loupe', { status: 'running', progress: 0, message: '匯出至 Loupe Browser...' })
    await exportLoupe({ output_dir: outDir || undefined, input_h5ad: inputH5ad || undefined })
    refetchLoupe()
  }

  const images: Record<string, string> = (resultImagesData as { status?: string; data?: Record<string, string> })?.data ?? {}

  // 偵測是否存在多 ROI 的 UMAP/dotplot/heatmap 結果
  const roiUmapKeys = Object.keys(images).filter(k => k.startsWith('result_umap'))
  const hasMultiRoi = roiUmapKeys.some(k => k !== 'result_umap')
  const resultRois = hasMultiRoi
    ? roiUmapKeys.filter(k => k !== 'result_umap').map(k => k.replace('result_umap_', '')).sort()
    : ['']

  // 當 images 更新時重置 activeResultRoi 到第一個 ROI
  useEffect(() => {
    if (resultRois.length > 0) setActiveResultRoi(resultRois[0])
  }, [images])

  // 收集所有空間圖，依 ROI 分組：{ roiSuffix: { outline?: string, filled?: string } }
  const spatialGroups: Record<string, { outline?: string; filled?: string }> = {}
  for (const key of Object.keys(images)) {
    const m = key.match(/^result_spatial(_filled)?(_(.+))?$/)
    if (!m) continue
    const isFilled = !!m[1]
    const roi = m[3] ?? ''
    if (!spatialGroups[roi]) spatialGroups[roi] = {}
    if (isFilled) spatialGroups[roi].filled = images[key]
    else spatialGroups[roi].outline = images[key]
  }
  const roiKeys = Object.keys(spatialGroups).sort()
  const hasImages = Object.keys(images).length > 0

  return (
    <div className="space-y-4">
      {/* Result Visualizations Section */}
      <div className="bg-surface-card rounded-xl border border-surface-border p-4">
        <div className="flex items-center justify-between mb-3">
          <h3 className="text-sm font-semibold text-gray-200">結果圖表</h3>
          <button
            onClick={handleGenerateResult}
            disabled={resultRunning}
            className="px-3 py-1.5 text-xs rounded bg-primary text-white hover:bg-primary/80 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {resultRunning ? '生成中...' : '生成結果圖'}
          </button>
        </div>

        {resultMessage && (
          <p className="text-xs text-gray-400 mb-3">{resultMessage}</p>
        )}

        {hasImages ? (
          <>
            {/* Tab bar */}
            <div className="flex gap-1 mb-4 border-b border-surface-border">
              {TABS.map(tab => (
                <button
                  key={tab.id}
                  onClick={() => setActiveTab(tab.id)}
                  className={`px-3 py-1.5 text-xs rounded-t transition-colors ${
                    activeTab === tab.id
                      ? 'bg-primary text-white border-b-2 border-primary'
                      : 'text-gray-400 hover:text-gray-200'
                  }`}
                >
                  {tab.label}
                </button>
              ))}
            </div>

            {/* Spatial tab */}
            {activeTab === 'spatial' && (
              <div className="space-y-6">
                {roiKeys.map(roi => (
                  <div key={roi}>
                    {roiKeys.length > 1 && (
                      <p className="text-xs font-semibold text-gray-300 mb-2">ROI: {roi || '(default)'}</p>
                    )}
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      {spatialGroups[roi].outline && (
                        <div>
                          <p className="text-xs text-gray-400 mb-1">彩色輪廓（H&amp;E 透明）</p>
                          <img src={`data:image/png;base64,${spatialGroups[roi].outline}`}
                               className="w-full rounded border border-surface-border" alt="spatial outline" />
                        </div>
                      )}
                      {spatialGroups[roi].filled && (
                        <div>
                          <p className="text-xs text-gray-400 mb-1">實心填色 + 輪廓</p>
                          <img src={`data:image/png;base64,${spatialGroups[roi].filled}`}
                               className="w-full rounded border border-surface-border" alt="spatial filled" />
                        </div>
                      )}
                    </div>
                  </div>
                ))}
                {roiKeys.length === 0 && (
                  <p className="text-xs text-gray-500">空間圖尚未生成</p>
                )}
              </div>
            )}

            {/* UMAP tab */}
            {activeTab === 'umap' && (
              <div>
                {hasMultiRoi && (
                  <div className="flex gap-1 mb-3">
                    {resultRois.map(roi => (
                      <button key={roi} onClick={() => setActiveResultRoi(roi)}
                        className={`px-2 py-1 text-xs rounded ${activeResultRoi === roi ? 'bg-primary/30 text-primary' : 'text-gray-400 hover:text-gray-200'}`}>
                        {roi}
                      </button>
                    ))}
                  </div>
                )}
                {(() => {
                  const uKey = hasMultiRoi ? `result_umap_${activeResultRoi}` : 'result_umap'
                  return images[uKey]
                    ? <img src={`data:image/png;base64,${images[uKey]}`}
                           className="w-full rounded border border-surface-border" alt="umap annotated" />
                    : <p className="text-xs text-gray-500">UMAP 尚未生成</p>
                })()}
              </div>
            )}

            {/* Dotplot tab */}
            {activeTab === 'dotplot' && (
              <div>
                {hasMultiRoi && (
                  <div className="flex gap-1 mb-3">
                    {resultRois.map(roi => (
                      <button key={roi} onClick={() => setActiveResultRoi(roi)}
                        className={`px-2 py-1 text-xs rounded ${activeResultRoi === roi ? 'bg-primary/30 text-primary' : 'text-gray-400 hover:text-gray-200'}`}>
                        {roi}
                      </button>
                    ))}
                  </div>
                )}
                {(() => {
                  const dKey = hasMultiRoi ? `result_dotplot_${activeResultRoi}` : 'result_dotplot'
                  return images[dKey]
                    ? <img src={`data:image/png;base64,${images[dKey]}`}
                           className="w-full rounded border border-surface-border" alt="dotplot" />
                    : <p className="text-xs text-gray-500">Dotplot 尚未生成</p>
                })()}
              </div>
            )}

            {/* Heatmap tab */}
            {activeTab === 'heatmap' && (
              <div>
                {hasMultiRoi && (
                  <div className="flex gap-1 mb-3">
                    {resultRois.map(roi => (
                      <button key={roi} onClick={() => setActiveResultRoi(roi)}
                        className={`px-2 py-1 text-xs rounded ${activeResultRoi === roi ? 'bg-primary/30 text-primary' : 'text-gray-400 hover:text-gray-200'}`}>
                        {roi}
                      </button>
                    ))}
                  </div>
                )}
                {(() => {
                  const hKey = hasMultiRoi ? `result_heatmap_${activeResultRoi}` : 'result_heatmap'
                  return images[hKey]
                    ? <img src={`data:image/png;base64,${images[hKey]}`}
                           className="w-full rounded border border-surface-border" alt="heatmap" />
                    : <p className="text-xs text-gray-500">Heatmap 尚未生成</p>
                })()}
              </div>
            )}
          </>
        ) : (
          <p className="text-xs text-gray-500">尚未生成結果圖，請點擊「生成結果圖」</p>
        )}
      </div>

      {/* Export Options */}
      <div className="bg-surface-card rounded-xl border border-surface-border p-4">
        <h3 className="text-sm font-semibold text-gray-200 mb-3">匯出設定</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div>
            <label className="text-xs text-gray-400">來源 h5ad 檔案名稱（選填）</label>
            <input value={inputH5ad} onChange={e => setInputH5ad(e.target.value)}
                   placeholder="預設：umap_computed.h5ad"
                   className="w-full mt-1 px-3 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none" />
          </div>
          <div>
            <label className="text-xs text-gray-400">輸出目錄（選填）</label>
            <input value={outDir} onChange={e => setOutDir(e.target.value)}
                   placeholder="預設：export_xenium"
                   className="w-full mt-1 px-3 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none" />
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
