import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'
import { getProsegRNAStatus, runProsegRNA, listProsegRNARois, getProsegComparison } from '../api/client'

interface RoiProsegInfo {
  name: string
  has_adata: boolean
  has_mask: boolean
  has_proseg: boolean
}

function ComparisonModal({ roiName, onClose }: { roiName: string; onClose: () => void }) {
  const [data, setData] = useState<{ he: string; cellpose: string; proseg: string; width: number; height: number } | null>(null)
  const [loading, setLoading] = useState(true)
  const [showHe, setShowHe] = useState(true)
  const [showCellpose, setShowCellpose] = useState(true)
  const [showProseg, setShowProseg] = useState(true)

  useEffect(() => {
    getProsegComparison(roiName).then(res => {
      if (res.data?.data) setData(res.data.data)
      setLoading(false)
    }).catch(() => setLoading(false))
  }, [roiName])

  if (loading) {
    return (
      <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-sm">
        <div className="text-purple-400 animate-pulse">載入比較資料中...</div>
      </div>
    )
  }

  return (
    <div className="fixed inset-0 z-50 flex flex-col items-center justify-center bg-black/90 p-4 sm:p-8">
      <div className="relative w-full max-w-5xl bg-surface-card rounded-2xl border border-surface-border overflow-hidden flex flex-col max-h-full shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between p-4 border-b border-surface-border bg-surface">
          <div>
            <h3 className="text-lg font-bold text-white flex items-center gap-2">
              <span className="w-2 h-2 rounded-full bg-purple-500"></span>
                分群邊界比較：{roiName}
            </h3>
            <p className="text-xs text-gray-400">交叉確認 Cellpose 與 Proseg 的分群結果</p>
          </div>
          <button onClick={onClose} className="p-2 hover:bg-white/10 rounded-full transition-colors">
            <svg className="w-6 h-6 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Controls */}
        <div className="p-3 bg-surface-card/50 border-b border-surface-border flex flex-wrap gap-4 items-center">
            <div className="text-[10px] uppercase tracking-wider text-gray-500 font-bold mr-2">圖層切換</div>
            
            <button 
                onClick={() => setShowHe(!showHe)}
                className={`flex items-center gap-2 px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                    showHe ? 'bg-gray-700 border-gray-500 text-white' : 'bg-transparent border-gray-700 text-gray-500'
                }`}
            >
                <div className={`w-2 h-2 rounded-full ${showHe ? 'bg-gray-300' : 'bg-gray-700'}`}></div>
                HE 影像
            </button>

            <button 
                onClick={() => setShowCellpose(!showCellpose)}
                className={`flex items-center gap-2 px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                    showCellpose ? 'bg-cyan-900/40 border-cyan-700 text-cyan-300' : 'bg-transparent border-gray-700 text-gray-500'
                }`}
            >
                <div className={`w-2 h-2 rounded-full ${showCellpose ? 'bg-cyan-400' : 'bg-gray-700'}`}></div>
                Cellpose 遮罩 (青色)
            </button>

            <button 
                onClick={() => setShowProseg(!showProseg)}
                className={`flex items-center gap-2 px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                    showProseg ? 'bg-red-900/40 border-red-700 text-red-300' : 'bg-transparent border-gray-700 text-gray-500'
                }`}
            >
                <div className={`w-2 h-2 rounded-full ${showProseg ? 'bg-red-500' : 'bg-gray-700'}`}></div>
                Proseg 輪廓 (紅色)
            </button>
        </div>

        {/* Viewer Area */}
        <div className="flex-1 overflow-auto bg-[#0a0a0a] relative flex items-center justify-center p-8">
            <div className="relative shadow-2xl" style={{ 
                width: data ? 'auto' : 0, 
                height: data ? 'auto' : 0,
                maxWidth: '100%',
                maxHeight: '70vh'
            }}>
                {data && (
                    <>
                        <img 
                            src={`data:image/jpeg;base64,${data.he}`} 
                            className={`max-w-full h-auto rounded transition-opacity duration-300 block ${showHe ? 'opacity-100' : 'opacity-0'}`}
                            style={{ imageRendering: 'auto' }}
                            alt="HE" 
                        />
                        {data.cellpose && (
                            <img 
                                src={`data:image/png;base64,${data.cellpose}`} 
                                className={`absolute inset-0 w-full h-full pointer-events-none transition-opacity duration-300 block ${showCellpose ? 'opacity-100' : 'opacity-0'}`}
                                style={{ imageRendering: 'pixelated' }}
                                alt="Cellpose" 
                            />
                        )}
                        {data.proseg && (
                            <img 
                                src={`data:image/png;base64,${data.proseg}`} 
                                className={`absolute inset-0 w-full h-full pointer-events-none transition-opacity duration-300 block ${showProseg ? 'opacity-100' : 'opacity-0'}`}
                                style={{ imageRendering: 'pixelated' }}
                                alt="Proseg" 
                            />
                        )}
                    </>
                )}
            </div>
            
            {!data && <div className="text-gray-600 italic">無可用比較影像</div>}
        </div>

        {/* Footer info */}
        <div className="p-3 bg-surface border-t border-surface-border text-[10px] text-gray-500 flex justify-between">
            <span>建議事項：Proseg 通常會依據 RNA 分布對邊界進行微調 (Dilation)</span>
            <span>{data?.width} x {data?.height} px</span>
        </div>
      </div>
    </div>
  )
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
  const [comparingRoi, setComparingRoi] = useState<string | null>(null)

  return (
    <div className="space-y-4">
      {comparingRoi && (
        <ComparisonModal
          roiName={comparingRoi}
          onClose={() => setComparingRoi(null)}
        />
      )}
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
                  <th className="text-center py-2 px-3 font-medium">比較與驗證</th>
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
                        ? <span className="text-green-400">✓</span>
                        : <span className="text-gray-600">—</span>}
                    </td>
                    <td className="py-2 px-3 text-center">
                      <button
                        onClick={() => setComparingRoi(roi.name)}
                        disabled={!roi.has_mask}
                        className="px-2 py-0.5 rounded text-xs font-medium transition-colors
                                   bg-cyan-900/40 border border-cyan-800 text-cyan-300
                                   hover:bg-cyan-800/60
                                   disabled:opacity-40 disabled:cursor-not-allowed flex items-center gap-1 mx-auto"
                      >
                        <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
                        </svg>
                        檢視輪廓
                      </button>
                    </td>
                    <td className="py-2 px-3 text-center">
                      <button
                        onClick={() => handleRunSingle(roi.name)}
                        disabled={!roi.has_adata || !roi.has_mask || stage.status === 'running'}
                        className="px-2 py-0.5 rounded text-xs font-medium transition-colors
                                   bg-purple-700/40 border border-purple-600 text-purple-300
                                   hover:bg-purple-600/60
                                   disabled:opacity-40 disabled:cursor-not-allowed mx-auto"
                      >
                        再次執行
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
