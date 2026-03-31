import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { getCellposeCountStatus, runCellposeCount, listCountRois } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'
import { useT } from '../i18n'

interface RoiCountInfo {
  name: string
  has_mask: boolean
  has_count: boolean
}

export default function Stage2_Count() {
  useStageLog('count')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['count']
  const { refetch: refetchStatus } = useStageStatus('count', getCellposeCountStatus, 3000)
  const [roiInfos, setRoiInfos] = useState<RoiCountInfo[]>([])
  const t = useT()

  useEffect(() => {
    listCountRois().then(res => {
      if (res.data?.data) setRoiInfos(res.data.data)
    }).catch(() => {})
  }, [stage.status])

  const handleRunAll = async () => {
    updateStage('count', { status: 'running', progress: 0, message: t('stage2.run_all') + '...' })
    try {
      await runCellposeCount(null)
    } catch (e) {
      updateStage('count', { status: 'error', message: `API error: ${e}` })
      return
    }
    refetchStatus()
  }

  const handleRunSingle = async (roiName: string) => {
    updateStage('count', { status: 'running', progress: 0, message: `RNA Count (${roiName})...` })
    try {
      await runCellposeCount(roiName)
    } catch (e) {
      updateStage('count', { status: 'error', message: `API error: ${e}` })
      return
    }
    refetchStatus()
  }

  const readyRois   = roiInfos.filter(r => r.has_mask)
  const doneRois    = roiInfos.filter(r => r.has_count)
  const missingRois = roiInfos.filter(r => !r.has_mask)

  return (
    <div className="space-y-4">
      <StageCard
        title={t('stage2.title')}
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
      >
        {/* Info box */}
        <div className="mt-3 p-3 rounded-lg bg-blue-900/20 border border-blue-700/40 text-xs text-blue-300 space-y-1">
          <p className="font-semibold text-blue-200">MCseg v2 RNA Counting</p>
          <ul className="list-disc pl-4 space-y-0.5 text-blue-400">
            <li>Input: <code>adata_002um.h5ad</code> (Stage 0), <code>segmentation_masks.npy</code> (Stage 1)</li>
            <li>Output: <code>cellpose_cells.h5ad</code> (cells × genes, for Stage 3)</li>
          </ul>
        </div>

        {/* Run button */}
        <div className="mt-4 flex items-center justify-between">
          <div className="text-xs text-gray-500">
            {doneRois.length > 0 && (
              <span className="text-green-400">✓ {doneRois.length} ROI(s) {t('stage2.done_rois')}</span>
            )}
            {readyRois.length > 0 && doneRois.length === 0 && (
              <span>{readyRois.length} ROI(s) {t('stage2.ready_rois')}</span>
            )}
          </div>
          <button
            onClick={handleRunAll}
            disabled={stage.status === 'running' || readyRois.length === 0}
            className="px-4 py-1.5 text-sm rounded-lg font-medium transition-colors
                       bg-brand-primary text-white hover:bg-brand-primary/90
                       disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {stage.status === 'running' ? t('common.running') : t('stage2.run_all')}
          </button>
        </div>
      </StageCard>

      {/* ROI status table */}
      {roiInfos.length > 0 && (
        <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-3">
          <h3 className="text-sm font-semibold text-gray-200">ROI Count Status</h3>

          {missingRois.length > 0 && (
            <div className="text-xs text-yellow-400 bg-yellow-900/20 rounded px-3 py-2">
              ⚠️ {t('stage2.missing_rois')}: {missingRois.map(r => r.name).join(', ')}
            </div>
          )}

          <div className="overflow-x-auto">
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="text-gray-500 border-b border-gray-700">
                  <th className="text-left py-2 pr-4 font-medium">ROI</th>
                  <th className="text-center py-2 px-3 font-medium">Mask</th>
                  <th className="text-center py-2 px-3 font-medium">Count</th>
                  <th className="text-center py-2 px-3 font-medium">{t('stage2.run_single')}</th>
                </tr>
              </thead>
              <tbody>
                {roiInfos.map(roi => (
                  <tr key={roi.name} className="border-b border-gray-800/60">
                    <td className="py-2 pr-4 text-gray-200 font-medium">{roi.name}</td>
                    <td className="py-2 px-3 text-center">
                      {roi.has_mask
                        ? <span className="text-green-400">✓</span>
                        : <span className="text-gray-600">—</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      {roi.has_count
                        ? <span className="text-green-400">✓ cellpose_cells.h5ad</span>
                        : <span className="text-gray-600">—</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      <button
                        onClick={() => handleRunSingle(roi.name)}
                        disabled={!roi.has_mask || stage.status === 'running'}
                        className="px-2 py-0.5 rounded text-xs font-medium transition-colors
                                   bg-blue-700/40 border border-blue-600 text-blue-300
                                   hover:bg-blue-600/60
                                   disabled:opacity-40 disabled:cursor-not-allowed"
                      >
                        {t('common.run')}
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      <Terminal stage="count" />
    </div>
  )
}
