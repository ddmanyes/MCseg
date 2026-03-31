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
import { useT } from '../i18n'

type TabId = 'spatial' | 'umap' | 'dotplot' | 'heatmap'

export default function Stage4_Export() {
  useStageLog('export')
  const { stages, updateStage } = usePipelineStore()
  const xenium = stages['xenium']
  const loupe  = stages['loupe']
  const { refetch: refetchXenium } = useStageStatus('xenium', getXeniumStatus, 3000)
  const { refetch: refetchLoupe }  = useStageStatus('loupe',  getLoupeStatus,  3000)
  const t = useT()

  const TABS: { id: TabId; label: string }[] = [
    { id: 'spatial',  label: t('stage4.tab.spatial') },
    { id: 'umap',     label: 'UMAP' },
    { id: 'dotplot',  label: 'Dotplot' },
    { id: 'heatmap',  label: 'Heatmap' },
  ]

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
      setResultMessage(t('stage4.result.done'))
      refetchResultImages()
    } else if (s.status === 'error') {
      setResultRunning(false)
      setResultMessage(s.message ?? t('stage4.result.failed'))
    } else if (s.status === 'running') {
      setResultMessage(s.message ?? t('stage4.result.generating'))
    }
  }, [resultStatusData])

  const handleGenerateResult = async () => {
    setResultRunning(true)
    setResultMessage(t('stage4.result.starting'))
    await generateResult()
    refetchResultStatus()
  }

  const handleXenium = async () => {
    updateStage('xenium', { status: 'running', progress: 0, message: t('stage4.xenium.starting') })
    await exportXenium({})
    refetchXenium()
  }

  const handleLoupe = async () => {
    updateStage('loupe', { status: 'running', progress: 0, message: t('stage4.loupe.starting') })
    await exportLoupe({})
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
          <h3 className="text-sm font-semibold text-gray-200">{t('stage4.result.title')}</h3>
          <button
            onClick={handleGenerateResult}
            disabled={resultRunning}
            className="px-3 py-1.5 text-xs rounded bg-primary text-white hover:bg-primary/80 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {resultRunning ? t('stage4.result.generating') : t('stage4.result.generate')}
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
                          <p className="text-xs text-gray-400 mb-1">{t('stage4.spatial.outline')}</p>
                          <img src={`data:image/png;base64,${spatialGroups[roi].outline}`}
                               className="w-full rounded border border-surface-border" alt="spatial outline" />
                        </div>
                      )}
                      {spatialGroups[roi].filled && (
                        <div>
                          <p className="text-xs text-gray-400 mb-1">{t('stage4.spatial.filled')}</p>
                          <img src={`data:image/png;base64,${spatialGroups[roi].filled}`}
                               className="w-full rounded border border-surface-border" alt="spatial filled" />
                        </div>
                      )}
                    </div>
                  </div>
                ))}
                {roiKeys.length === 0 && (
                  <p className="text-xs text-gray-500">{t('common.not_generated')}</p>
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
                    : <p className="text-xs text-gray-500">{t('common.not_generated')}</p>
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
                    : <p className="text-xs text-gray-500">{t('common.not_generated')}</p>
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
                    : <p className="text-xs text-gray-500">{t('common.not_generated')}</p>
                })()}
              </div>
            )}
          </>
        ) : (
          <p className="text-xs text-gray-500">{t('stage4.result.no_images')}</p>
        )}
      </div>


      <StageCard title={t('stage4.xenium.title')} status={xenium.status}
                 progress={xenium.progress} message={xenium.message}
                 onRun={handleXenium} runLabel={t('stage4.xenium.run')}>
        <div className="text-sm text-gray-400 space-y-1">
          <p>{t('stage4.xenium.output_desc')}</p>
          <p className="text-yellow-400 text-xs">{t('stage4.xenium.note_text')}</p>
        </div>
      </StageCard>

      <StageCard title={t('stage4.loupe.title')} status={loupe.status}
                 progress={loupe.progress} message={loupe.message}
                 onRun={handleLoupe} runLabel={t('stage4.loupe.run')}>
        <div className="text-sm text-gray-400 space-y-1">
          <p>{t('stage4.loupe.output_desc')}</p>
          <p className="text-yellow-400 text-xs">{t('stage4.loupe.note_text')}</p>
        </div>
      </StageCard>

      <Terminal stage="export" />
    </div>
  )
}
