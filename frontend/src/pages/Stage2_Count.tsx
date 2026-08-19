import { useState, useEffect, useRef } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { getCellposeCountStatus, runCellposeCount, listCountRois, runFullCount, getFullCountStatus, runCoverageQc } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'
import { useT } from '../i18n'
import { errText } from '../utils/errText'
import { startStatusPoll, type PollStatus } from '../utils/pollStatus'

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
  const [roiListWarn, setRoiListWarn] = useState<string | null>(null)
  const t = useT()

  useEffect(() => {
    // B 類讀取：失敗只是清單空白，做可見降級即可
    listCountRois().then(res => {
      if (res.data?.data) setRoiInfos(res.data.data)
      setRoiListWarn(null)
    }).catch((e: unknown) => {
      setRoiListWarn(t('stage2.warn.roi_load_failed', { err: errText(e) }))
    })
  }, [stage.status])

  const handleRunAll = async () => {
    updateStage('count', { status: 'running', progress: 0, message: t('stage2.run_all') + '...' })
    try {
      await runCellposeCount(null)
    } catch (e: unknown) {
      updateStage('count', { status: 'error', message: errText(e) })
      return
    }
    void refetchStatus()
  }

  const handleRunSingle = async (roiName: string) => {
    updateStage('count', { status: 'running', progress: 0, message: `RNA Count (${roiName})...` })
    try {
      await runCellposeCount(roiName)
    } catch (e: unknown) {
      updateStage('count', { status: 'error', message: errText(e) })
      return
    }
    void refetchStatus()
  }

  // ── 全圖計數（吃 Stage 1 的 full_image_segmentation_masks.npy）─────────────
  const [fullStatus, setFullStatus] = useState<PollStatus | null>(null)
  const [fullRetry, setFullRetry] = useState<string | null>(null)
  const fullPollRef = useRef<ReturnType<typeof setInterval>>()

  const startFullPoll = () => {
    clearInterval(fullPollRef.current)
    setFullRetry(null)
    fullPollRef.current = startStatusPoll({
      fetchStatus: async () => {
        const res = await getFullCountStatus()
        return (res.data?.data ?? res.data) as PollStatus | null
      },
      onStatus: d => { setFullRetry(null); setFullStatus(d) },
      onTransientFailure: (n, max) => setFullRetry(t('stage2.msg.retrying', { n, max })),
      onLost: msg => { setFullRetry(null); setFullStatus({ status: 'error', message: msg }) },
    })
  }

  useEffect(() => {
    getFullCountStatus().then(res => {
      const d = res.data?.data ?? res.data
      if (d) {
        setFullStatus(d)
        if (d.status === 'running') startFullPoll()
      }
    }).catch((e: unknown) => {
      setFullRetry(t('stage2.warn.full_status_failed', { err: errText(e) }))
    })
    return () => clearInterval(fullPollRef.current)
  }, [])

  const handleRunFullCount = async () => {
    setFullStatus({ status: 'running', progress: 0, message: t('stage2.full_count.starting') })
    try {
      const res = await runFullCount()
      if (res.data?.status === 'error') {
        setFullStatus({ status: 'error', message: res.data.message })
        return
      }
      startFullPoll()
    } catch (e: any) {
      setFullStatus({ status: 'error', message: e?.response?.data?.message ?? 'API error' })
    }
  }

  // ── 分割覆蓋率 QC（切片層級為主要結論）──────────────────────────────────
  interface QcSection {
    section: number
    x0: number; y0: number; x1: number; y1: number
    n_bins: number
    n_cells: number
    cells_per_1k_bins: number
    median_cell_area_px: number
    flagged: boolean
    flag_reason: string
  }
  interface QcResult {
    sections: QcSection[]
    sections_summary: { median_cells_per_1k_bins?: number; n_flagged?: number }
    grid_summary: { global_density_ok?: boolean; n_grid_flagged?: number }
    grid_flagged: unknown[]
  }

  const [qc, setQc] = useState<QcResult | null>(null)
  const [qcBusy, setQcBusy] = useState(false)
  const [qcError, setQcError] = useState<string | null>(null)

  /**
   * QC 回應結構防呆。
   *
   * 後端若回非預期結構（例如 SPA catch-all 把 index.html 當成 JSON 回來），
   * 直接 setQc 會讓下方 `qc.sections.map` 在 render 期間拋錯 → 整頁空白。
   * 這裡先驗形狀，形狀不對就走 qcError 路徑。
   */
  const isQcResult = (d: unknown): d is QcResult =>
    !!d && typeof d === 'object' && Array.isArray((d as QcResult).sections)

  const handleRunQc = async () => {
    setQcBusy(true); setQcError(null)
    try {
      const res = await runCoverageQc()
      if (res.data?.status === 'error') setQcError(res.data.message)
      else if (isQcResult(res.data?.data)) { setQc(res.data.data); setQcError(null) }
      else setQcError(t('stage2.err.qc_structure'))
    } catch (e: unknown) {
      setQcError(errText(e))
    } finally {
      setQcBusy(false)
    }
  }


  const qcMedian = qc?.sections_summary?.median_cells_per_1k_bins ?? 0

  const readyRois   = roiInfos.filter(r => r.has_mask)
  const doneRois    = roiInfos.filter(r => r.has_count)
  const missingRois = roiInfos.filter(r => !r.has_mask)

  return (
    <div className="space-y-6">
      {/* ── 主流程：ROI 空間細胞歸屬（標準工作流） ── */}
      <StageCard
        title={t('stage2.title')}
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
        onRun={handleRunAll}
        runLabel={t('stage2.run_all')}
        disabled={readyRois.length === 0}
      >
        {/* Info box */}
        <div className="mt-2 p-3 rounded-lg bg-blue-950/30 border border-blue-800/40 text-xs text-blue-300 space-y-1.5">
          <p className="font-semibold text-blue-200">MCseg v2 — Spatial Cell Attribution (Bin-to-Cell Mapping)</p>
          <p className="text-blue-300/90">{t('stage2.subtitle')}</p>
          <ul className="list-disc pl-4 space-y-0.5 text-blue-400 font-mono text-[11px]">
            <li>{t('stage2.info.input')}</li>
            <li>{t('stage2.info.output')}</li>
          </ul>
        </div>

        {roiListWarn && (
          <p className="mt-3 text-xs text-amber-400/80">ⓘ {roiListWarn}</p>
        )}

        {/* ROI Status Table */}
        <div className="mt-4 border-t border-surface-border pt-4 space-y-3">
          <div className="flex items-center justify-between">
            <h4 className="text-xs font-semibold text-gray-300 uppercase tracking-wider">
              {t('stage2.table.title')}
            </h4>
            <div className="text-xs">
              {doneRois.length > 0 && (
                <span className="text-emerald-400 font-medium">✓ {doneRois.length} {t('stage2.done_rois')}</span>
              )}
              {readyRois.length > 0 && doneRois.length === 0 && (
                <span className="text-blue-400">{readyRois.length} {t('stage2.ready_rois')}</span>
              )}
            </div>
          </div>

          {missingRois.length > 0 && (
            <div className="text-xs text-yellow-400 bg-yellow-900/20 border border-yellow-800/40 rounded px-3 py-2">
              ⚠️ {missingRois.length} {t('stage2.missing_rois')}: {missingRois.map(r => r.name).join(', ')}
            </div>
          )}

          {roiInfos.length === 0 ? (
            <p className="text-xs text-gray-500 py-3">{t('stage0.no_rois')}</p>
          ) : (
            <div className="overflow-x-auto rounded-lg border border-surface-border bg-surface/50">
              <table className="w-full text-xs border-collapse">
                <thead>
                  <tr className="text-gray-400 border-b border-surface-border bg-surface/80">
                    <th className="text-left py-2.5 px-3 font-medium">{t('stage2.table.roi')}</th>
                    <th className="text-center py-2.5 px-3 font-medium">{t('stage2.table.mask')}</th>
                    <th className="text-center py-2.5 px-3 font-medium">{t('stage2.table.count')}</th>
                    <th className="text-center py-2.5 px-3 font-medium">{t('stage2.table.action')}</th>
                  </tr>
                </thead>
                <tbody>
                  {roiInfos.map(roi => (
                    <tr key={roi.name} className="border-b border-surface-border/50 hover:bg-white/[0.02] transition-colors">
                      <td className="py-2.5 px-3 text-gray-200 font-semibold">{roi.name}</td>
                      <td className="py-2.5 px-3 text-center">
                        {roi.has_mask
                          ? <span className="inline-flex items-center gap-1 text-emerald-400 font-medium">✓ segmentation_masks.npy</span>
                          : <span className="text-gray-600">—</span>}
                      </td>
                      <td className="py-2.5 px-3 text-center">
                        {roi.has_count
                          ? <span className="inline-flex items-center gap-1 text-emerald-400 font-medium">✓ cellpose_cells.h5ad</span>
                          : <span className="text-gray-600">—</span>}
                      </td>
                      <td className="py-2.5 px-3 text-center">
                        <button
                          onClick={() => handleRunSingle(roi.name)}
                          disabled={!roi.has_mask || stage.status === 'running'}
                          className="px-3 py-1 rounded text-xs font-medium transition-all
                                     bg-primary/20 border border-primary/40 text-primary hover:bg-primary hover:text-white
                                     disabled:opacity-40 disabled:cursor-not-allowed"
                        >
                          {t('stage2.run_single')}
                        </button>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {/* Prominent Action Bar */}
          <div className="pt-2 flex items-center justify-end">
            <button
              onClick={handleRunAll}
              disabled={stage.status === 'running' || readyRois.length === 0}
              className="px-6 py-2.5 text-sm rounded-lg font-semibold transition-all
                         bg-primary hover:bg-primary-dark text-white shadow-lg shadow-blue-500/25
                         disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
            >
              <span>{stage.status === 'running' ? t('stage2.running') : `⚡ ${t('stage2.run_all')}`}</span>
            </button>
          </div>
        </div>
      </StageCard>

      {/* ── 進階全片處理區塊 (Fullslide & QC) ── */}
      <div className="rounded-xl border border-surface-border bg-surface-card/60 p-5 space-y-4">
        <div className="flex items-center gap-2.5 border-b border-surface-border pb-3">
          <span className="px-2 py-0.5 text-[11px] font-bold bg-indigo-950 text-indigo-300 border border-indigo-700/50 rounded uppercase tracking-wider">
            {t('common.advanced')}
          </span>
          <div>
            <h3 className="text-sm font-semibold text-gray-200">{t('stage2.advanced.title')}</h3>
            <p className="text-xs text-gray-400 mt-0.5">{t('stage2.advanced.subtitle')}</p>
          </div>
        </div>

        {/* 全圖計數卡片 */}
        <div className="rounded-lg border border-surface-border/80 bg-surface/60 p-4 space-y-3">
          <div className="flex items-center justify-between gap-4">
            <div>
              <h4 className="text-sm font-semibold text-gray-300">{t('stage2.full_count.title')}</h4>
              <p className="text-xs text-gray-500 mt-0.5">
                {t('stage2.full_count.description')}
                <code className="mx-1 text-yellow-400 font-mono">fullslide/cells.h5ad</code>
              </p>
            </div>
            <button
              onClick={handleRunFullCount}
              disabled={fullStatus?.status === 'running'}
              className="shrink-0 px-4 py-1.5 text-xs rounded-lg font-medium transition-colors
                         bg-indigo-700/70 hover:bg-indigo-600 border border-indigo-600/60 disabled:opacity-50 disabled:cursor-not-allowed text-indigo-100"
            >
              {fullStatus?.status === 'running' ? t('common.running') : t('stage2.full_count.run')}
            </button>
          </div>

          {fullRetry && (
            <p className="text-xs text-amber-400/80">ⓘ {fullRetry}</p>
          )}

          {fullStatus && fullStatus.status !== 'idle' && (
            <div className={`rounded-lg px-3 py-2 text-xs font-mono space-y-1
              ${fullStatus.status === 'error' ? 'bg-red-900/30 text-red-300 border border-red-800'
                : fullStatus.status === 'done' ? 'bg-emerald-900/30 text-emerald-300 border border-emerald-800'
                : 'bg-gray-800 text-gray-300 border border-gray-700'}`}
            >
              <div className="flex items-center justify-between">
                <span>
                  {fullStatus.status === 'running' && '⏳ '}
                  {fullStatus.status === 'done' && '✓ '}
                  {fullStatus.status === 'error' && '✗ '}
                  {fullStatus.message ?? fullStatus.status}
                </span>
                {fullStatus.progress != null && fullStatus.status === 'running' && (
                  <span className="text-gray-400">{Math.round(fullStatus.progress * 100)}%</span>
                )}
              </div>
              {fullStatus.status === 'running' && fullStatus.progress != null && (
                <div className="w-full bg-gray-700 rounded-full h-1.5 overflow-hidden">
                  <div
                    className="bg-indigo-500 h-full transition-all duration-500"
                    style={{ width: `${Math.round(fullStatus.progress * 100)}%` }}
                  />
                </div>
              )}
            </div>
          )}
        </div>

        {/* 分割覆蓋率 QC 卡片 */}
        <div className="rounded-lg border border-surface-border/80 bg-surface/60 p-4 space-y-3">
          <div className="flex items-start justify-between gap-4">
            <div>
              <h4 className="text-sm font-semibold text-gray-300">{t('stage2.qc.title')}</h4>
              <p className="text-xs text-gray-500 mt-0.5 max-w-2xl">{t('stage2.qc.description')}</p>
            </div>
            <button
              onClick={handleRunQc}
              disabled={qcBusy}
              className="shrink-0 px-4 py-1.5 text-xs rounded-lg font-medium transition-colors
                         bg-teal-800/70 hover:bg-teal-700 border border-teal-600/60 disabled:opacity-50 disabled:cursor-not-allowed text-teal-100"
            >
              {qcBusy ? t('stage2.qc.running') : t('stage2.qc.run')}
            </button>
          </div>

          {qcError && (
            <div className="rounded-lg px-3 py-2 text-xs bg-red-900/30 text-red-300 border border-red-800">
              ✗ {qcError}
            </div>
          )}

          {qc && (
            <div className="space-y-2">
              {qc.grid_summary?.global_density_ok === false && (
                <div className="rounded-lg px-3 py-2 text-xs bg-yellow-900/30 text-yellow-300 border border-yellow-700">
                  ⚠️ {t('stage2.qc.global_warn')}
                </div>
              )}

              <div className={`rounded-lg px-3 py-2 text-xs border
                ${(qc.sections_summary?.n_flagged ?? 0) > 0
                  ? 'bg-yellow-900/30 text-yellow-300 border-yellow-700'
                  : 'bg-emerald-900/30 text-emerald-300 border-emerald-800'}`}
              >
                {(qc.sections_summary?.n_flagged ?? 0) > 0
                  ? `⚠️ ${qc.sections_summary.n_flagged} ${t('stage2.qc.flagged')}`
                  : `✓ ${t('stage2.qc.ok')}（${qc.sections.length} ${t('stage2.qc.sections')}）`}
                {(qc.grid_summary?.n_grid_flagged ?? 0) > 0 && (
                  <span className="text-gray-400 ml-2">
                    · {qc.grid_summary.n_grid_flagged} {t('stage2.qc.grid_hint')}
                  </span>
                )}
              </div>

              <div className="overflow-x-auto">
                <table className="w-full text-xs border-collapse">
                  <thead>
                    <tr className="text-gray-500 border-b border-gray-700">
                      <th className="text-left py-2 pr-4 font-medium">{t('stage2.qc.section')}</th>
                      <th className="text-left py-2 px-3 font-medium">{t('stage2.qc.position')}</th>
                      <th className="text-right py-2 px-3 font-medium">{t('stage2.qc.bins')}</th>
                      <th className="text-right py-2 px-3 font-medium">{t('stage2.qc.cells')}</th>
                      <th className="text-right py-2 px-3 font-medium">{t('stage2.qc.density')}</th>
                      <th className="text-right py-2 px-3 font-medium">{t('stage2.qc.vs_median')}</th>
                    </tr>
                  </thead>
                  <tbody>
                    {qc.sections.map(s => (
                      <tr
                        key={s.section}
                        className={`border-b border-gray-800/60 ${s.flagged ? 'bg-yellow-900/15' : ''}`}
                        title={s.flag_reason || undefined}
                      >
                        <td className="py-2 pr-4 font-medium text-gray-200">
                          {s.flagged && <span className="text-yellow-400 mr-1">⚠️</span>}
                          {s.section}
                        </td>
                        <td className="py-2 px-3 text-gray-400 font-mono">
                          ({s.x0.toLocaleString()}, {s.y0.toLocaleString()})
                        </td>
                        <td className="py-2 px-3 text-right text-gray-300">{s.n_bins.toLocaleString()}</td>
                        <td className="py-2 px-3 text-right text-gray-300">{s.n_cells.toLocaleString()}</td>
                        <td className={`py-2 px-3 text-right font-mono ${s.flagged ? 'text-yellow-300' : 'text-gray-300'}`}>
                          {s.cells_per_1k_bins.toFixed(2)}
                        </td>
                        <td className={`py-2 px-3 text-right font-mono ${s.flagged ? 'text-yellow-300' : 'text-gray-500'}`}>
                          {qcMedian > 0 ? `${(s.cells_per_1k_bins / qcMedian).toFixed(2)}×` : '—'}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              {qc.sections.some(s => s.flagged) && (
                <p className="text-xs text-gray-500">
                  {qc.sections.filter(s => s.flagged).map(s => `#${s.section}: ${s.flag_reason}`).join(' · ')}
                </p>
              )}
            </div>
          )}
        </div>
      </div>

      <Terminal stage="count" />
    </div>
  )
}

