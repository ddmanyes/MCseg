import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { getProsegRNAStatus, runProsegRNA, listProsegRNARois } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

interface RoiProsegInfo {
  name: string
  has_adata: boolean
  has_mask: boolean
  has_proseg: boolean
}

export default function Stage25_ProsegRNA() {
  useStageLog('proseg_rna')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['proseg_rna']
  const { refetch: refetchStatus } = useStageStatus('proseg_rna', getProsegRNAStatus, 3000)
  const [roiInfos, setRoiInfos] = useState<RoiProsegInfo[]>([])

  useEffect(() => {
    listProsegRNARois().then(res => {
      if (res.data?.data) setRoiInfos(res.data.data)
    }).catch(() => {})
  }, [stage.status])

  const handleRunAll = async () => {
    updateStage('proseg_rna', { status: 'running', progress: 0, message: 'Proseg RNA 重分配（全部 ROI）...' })
    try {
      await runProsegRNA(null)
    } catch (e) {
      updateStage('proseg_rna', { status: 'error', message: `API 呼叫失敗：${e}` })
      return
    }
    refetchStatus()
  }

  const handleRunSingle = async (roiName: string) => {
    updateStage('proseg_rna', { status: 'running', progress: 0, message: `Proseg RNA 重分配（${roiName}）...` })
    try {
      await runProsegRNA(roiName)
    } catch (e) {
      updateStage('proseg_rna', { status: 'error', message: `API 呼叫失敗：${e}` })
      return
    }
    refetchStatus()
  }

  const readyRois   = roiInfos.filter(r => r.has_adata && r.has_mask)
  const doneRois    = roiInfos.filter(r => r.has_proseg)
  const missingRois = roiInfos.filter(r => !r.has_adata || !r.has_mask)

  return (
    <div className="space-y-4">
      <StageCard
        title="Proseg RNA 重分配（Stage 2.5）"
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
      >
        {/* 說明 */}
        <div className="mt-3 p-3 rounded-lg bg-purple-900/20 border border-purple-700/40 text-xs text-purple-300 space-y-1">
          <p className="font-semibold text-purple-200">可選階段：Proseg MCMC RNA 重分配</p>
          <p>
            沿用 Cellpose 分割遮罩（細胞邊界不變），以 Proseg MCMC 對
            Visium HD bins 進行 RNA 重分配，提升細胞邊界處的分配精確度。
          </p>
          <ul className="list-disc pl-4 space-y-0.5 text-purple-400">
            <li>輸入：<code>adata_002um.h5ad</code>、<code>segmentation_masks.npy</code>（沿用 Cellpose）</li>
            <li>固定參數：max_dist=20 µm、compactness=0.06、dilation=5 px</li>
            <li>輸出：<code>proseg_cells.h5ad</code>（cells × genes，供 Stage 3 選用）</li>
            <li>Stage 3 可選擇使用 Cellpose 或 Proseg 版本的 RNA 計數</li>
          </ul>
          <p className="text-purple-500 text-[11px]">
            ⚠️ 需要 proseg binary（<code>~/.cargo/bin/proseg</code> 或 config 中設定的路徑）
          </p>
        </div>

        {/* 執行按鈕 */}
        <div className="mt-4 flex items-center justify-between">
          <div className="text-xs text-gray-500">
            {doneRois.length > 0 && (
              <span className="text-green-400">✓ {doneRois.length} 個 ROI 已完成 Proseg RNA 重分配</span>
            )}
            {readyRois.length > 0 && doneRois.length === 0 && (
              <span>{readyRois.length} 個 ROI 就緒</span>
            )}
          </div>
          <button
            onClick={handleRunAll}
            disabled={stage.status === 'running' || readyRois.length === 0}
            className="px-4 py-1.5 text-sm rounded-lg font-medium transition-colors
                       bg-purple-700 text-white hover:bg-purple-600
                       disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {stage.status === 'running' ? 'MCMC 執行中...' : '全部 ROI 執行 Proseg RNA'}
          </button>
        </div>
      </StageCard>

      {/* ROI 狀態表格 */}
      {roiInfos.length > 0 && (
        <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-3">
          <h3 className="text-sm font-semibold text-gray-200">ROI Proseg RNA 狀態</h3>

          {missingRois.length > 0 && (
            <div className="text-xs text-yellow-400 bg-yellow-900/20 rounded px-3 py-2">
              ⚠️ 以下 ROI 缺少必要輸入，請先完成 Stage 0/1：
              {missingRois.map(r => r.name).join('、')}
            </div>
          )}

          <div className="overflow-x-auto">
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="text-gray-500 border-b border-gray-700">
                  <th className="text-left py-2 pr-4 font-medium">ROI</th>
                  <th className="text-center py-2 px-3 font-medium">Bins (S0)</th>
                  <th className="text-center py-2 px-3 font-medium">Mask (S1)</th>
                  <th className="text-center py-2 px-3 font-medium">Proseg 結果</th>
                  <th className="text-center py-2 px-3 font-medium">單獨執行</th>
                </tr>
              </thead>
              <tbody>
                {roiInfos.map(roi => (
                  <tr key={roi.name} className="border-b border-gray-800/60">
                    <td className="py-2 pr-4 text-gray-200 font-medium">{roi.name}</td>
                    <td className="py-2 px-3 text-center">
                      {roi.has_adata
                        ? <span className="text-green-400">✓</span>
                        : <span className="text-gray-600">—</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      {roi.has_mask
                        ? <span className="text-green-400">✓</span>
                        : <span className="text-gray-600">—</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      {roi.has_proseg
                        ? <span className="text-green-400">✓ proseg_cells.h5ad</span>
                        : <span className="text-gray-600">尚未執行</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      <button
                        onClick={() => handleRunSingle(roi.name)}
                        disabled={!roi.has_adata || !roi.has_mask || stage.status === 'running'}
                        className="px-2 py-0.5 rounded text-xs font-medium transition-colors
                                   bg-purple-700/40 border border-purple-600 text-purple-300
                                   hover:bg-purple-600/60
                                   disabled:opacity-40 disabled:cursor-not-allowed"
                      >
                        執行
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <p className="text-xs text-gray-600">
            ⓘ 此步驟為選用，可直接跳過進入 Stage 3（將使用 Cellpose 直接計數結果）
          </p>
        </div>
      )}

      <Terminal stage="proseg_rna" />
    </div>
  )
}
