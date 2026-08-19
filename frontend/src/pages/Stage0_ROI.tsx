import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { listRois, addRoi, deleteRoi, runRoiExtract, getRoiStatus, getConfig, getRoiBounds, type RnaBounds } from '../api/client'
import type { RoiDefinition } from '../types/pipeline'
import useStageLog from '../hooks/useStageLog'
import RoiSelector from '../components/roi/RoiSelector'
import AlignmentPanel from '../components/shared/AlignmentPanel'
import { useStageStatus } from '../hooks/useStageStatus'
import { useT } from '../i18n'
import { errText, apiErrorMessage } from '../utils/errText'

export default function Stage0_ROI() {
  useStageLog('roi')
  const { stages, updateStage, rois, setRois } = usePipelineStore()
  const stage = stages['roi']
  const { refetch: refetchStatus } = useStageStatus('roi', getRoiStatus, 2000)
  const [form, setForm] = useState<Partial<RoiDefinition>>({ pixel_size_um: 0.2737 })
  const [formError, setFormError] = useState<string | null>(null)
  const [configPixelSize, setConfigPixelSize] = useState<number>(0.2737)
  const [captureBounds, setCaptureBounds] = useState<RnaBounds | null>(null)
  // B 類讀取降級：載入失敗只是清單／預填值取不到，不擋人
  const [loadWarn, setLoadWarn] = useState<string | null>(null)
  const t = useT()

  useEffect(() => {
    listRois()
      .then(r => setRois(r.data.data ?? []))
      .catch((e: unknown) => setLoadWarn(t('stage0.warn.load_failed', { err: errText(e) })))
    getRoiBounds()
      .then(r => { if (r.data?.data) setCaptureBounds(r.data.data) })
      .catch(() => {})
    // 從 config 讀取 pixel_size_um（Data Setup 掃描時寫入）
    getConfig().then((r: any) => {
      const ps = r.data?.data?.global?.pixel_size_um
      if (ps && typeof ps === 'number' && ps > 0) {
        setConfigPixelSize(ps)
        setForm(f => ({ ...f, pixel_size_um: ps }))
      }
    }).catch((e: unknown) => {
      setLoadWarn(t('stage0.warn.pixel_size_failed', { err: errText(e), size: configPixelSize }))
    })
  }, [])


  const handleRun = async () => {
    updateStage('roi', { status: 'running', progress: 0, message: t('stage0.running') })
    try {
      // A 類寫入：HTTP 200 也可能是失敗（例如 ROI 清單為空），必須檢查 body
      const res = await runRoiExtract()
      const msg = apiErrorMessage(res)
      if (msg) {
        updateStage('roi', { status: 'error', progress: 0, message: msg })
        setFormError(msg)
        return
      }
      setFormError(null)
      void refetchStatus()
    } catch (e: unknown) {
      updateStage('roi', { status: 'error', message: errText(e) })
    }
  }

  // A 類寫入：刪除失敗必須報錯，否則 ROI 只是「畫面上消失」
  const handleDelete = async (name: string) => {
    try {
      const res = await deleteRoi(name)
      const msg = apiErrorMessage(res)
      if (msg) { setFormError(t('stage0.err.delete_failed', { name, msg })); return }
      const updated = await listRois()
      setRois(updated.data.data ?? [])
      setFormError(null)
    } catch (e: unknown) {
      setFormError(t('stage0.err.delete_failed', { name, msg: errText(e) }))
    }
  }

  const handleAdd = async () => {
    if (!form.name) { setFormError(t('stage0.form.name') + ' required'); return }
    if (!form.tissue) { setFormError(t('stage0.form.tissue') + ' required'); return }
    setFormError(null)
    // A 類寫入：失敗必須報錯，否則使用者以為 ROI 已新增
    try {
      const res = await addRoi(form as RoiDefinition)
      // HTTP 200 也可能是失敗（`{"status":"error"}`）
      const msg = apiErrorMessage(res)
      if (msg) { setFormError(t('stage0.err.add_failed', { msg })); return }
      const updated = await listRois()
      setRois(updated.data.data ?? [])
      setForm({ pixel_size_um: configPixelSize })
    } catch (e: unknown) {
      setFormError(t('stage0.err.add_failed_retry', { err: errText(e) }))
    }
  }


  const isOutOfBounds = Boolean(
    captureBounds &&
    form.x != null && form.width_px != null &&
    form.y != null && form.height_px != null &&
    (
      form.x + form.width_px < captureBounds.min_x ||
      form.x > captureBounds.max_x ||
      form.y + form.height_px < captureBounds.min_y ||
      form.y > captureBounds.max_y
    )
  )

  return (
    <div className="space-y-4">
      <StageCard title={t('stage0.title')} status={stage.status} progress={stage.progress}
        message={stage.message} onRun={handleRun} runLabel={t('stage0.run')}>

        {/* ROI list */}
        <div className="space-y-2">
          <p className="text-xs text-gray-400 font-medium uppercase tracking-wide">{t('stage0.defined_rois')}</p>
          {rois.length === 0 && <p className="text-sm text-gray-500">{t('stage0.no_rois')}</p>}
          {rois.map(roi => (
            <div key={roi.name} className="flex items-center justify-between bg-surface/50 rounded-lg px-3 py-2">
              <div>
                <span className="text-sm font-medium text-gray-200">{roi.name}</span>
                <span className="text-xs text-gray-400 ml-2">({roi.tissue})</span>
                {'x' in roi && (
                  <span className="text-xs text-gray-500 ml-2">
                    x={roi.x}, y={roi.y}, w={roi.width_px}, h={roi.height_px}
                  </span>
                )}
              </div>
              <button
                onClick={() => void handleDelete(roi.name)}
                className="text-red-400 hover:text-red-300 text-xs"
              >
                {t('common.delete')}
              </button>
            </div>
          ))}
        </div>

        {/* Interactive ROI selector */}
        <div className="border-t border-surface-border pt-4">
          <div className="flex items-center justify-between mb-3">
            <p className="text-xs text-gray-400 font-medium uppercase tracking-wide">{t('stage0.interactive')}</p>
            {captureBounds && (
              <span className="text-xs font-mono text-emerald-400 bg-emerald-950/40 border border-emerald-800/40 px-2 py-0.5 rounded">
                {t('stage0.capture_bounds.badge', {
                  min_x: Math.round(captureBounds.min_x),
                  max_x: Math.round(captureBounds.max_x),
                  min_y: Math.round(captureBounds.min_y),
                  max_y: Math.round(captureBounds.max_y),
                  bins: captureBounds.total_bins.toLocaleString(),
                })}
              </span>
            )}
          </div>
          <RoiSelector
            existingRois={rois as any}
            captureBounds={captureBounds}
            onSelect={(roi) => setForm(f => ({ ...f, ...roi }))}
          />
        </div>

        {/* 對位檢查 */}
        <div className="border-t border-surface-border pt-4">
          <p className="text-xs text-gray-400 font-medium uppercase tracking-wide">{t('align.title')}</p>
          <p className="text-xs text-gray-500 mb-3">{t('align.subtitle')}</p>
          <AlignmentPanel />
        </div>

        {/* Add ROI form */}
        <div className="border-t border-surface-border pt-4">
          <p className="text-xs text-gray-400 font-medium uppercase tracking-wide mb-3">{t('stage0.add_roi')}</p>
          <div className="grid grid-cols-2 gap-3">
            {([
              ['name',          t('stage0.form.name'),      'text'],
              ['tissue',        t('stage0.form.tissue'),     'text'],
              ['x',             t('stage0.form.x'),          'number'],
              ['y',             t('stage0.form.y'),          'number'],
              ['width_px',      t('stage0.form.width'),      'number'],
              ['height_px',     t('stage0.form.height'),     'number'],
              ['pixel_size_um', 'Pixel Size (µm/px)',        'number'],
            ] as [string, string, string][]).map(([key, label, type]) => (
              <div key={key} className={key === 'pixel_size_um' ? 'col-span-2' : ''}>
                <label className="text-xs text-gray-400">{label}</label>
                <input
                  type={type}
                  value={(form as any)[key] ?? ''}
                  onChange={e => {
                    setFormError(null)
                    setForm(f => ({
                      ...f,
                      [key]: type === 'number' ? Number(e.target.value) : e.target.value,
                    }))
                  }}
                  className="w-full mt-1 px-2 py-1.5 bg-surface border border-surface-border rounded text-sm text-gray-200 focus:border-primary focus:outline-none"
                />
              </div>
            ))}
          </div>

          {isOutOfBounds && (
            <div className="mt-3 p-2.5 bg-amber-950/40 border border-amber-500/40 rounded text-xs text-amber-300 flex items-start gap-2">
              <span>{t('stage0.warn.out_of_bounds')}</span>
            </div>
          )}

          <button
            onClick={handleAdd}
            className="mt-3 px-4 py-1.5 bg-surface-border hover:bg-surface-border/80 rounded text-sm text-gray-200 transition-colors"
          >
            + {t('stage0.add_roi')}
          </button>
          {formError && <p className="mt-2 text-xs text-red-400">{formError}</p>}
          {loadWarn && <p className="mt-2 text-xs text-amber-400/80">ⓘ {loadWarn}</p>}
        </div>
      </StageCard>

      <Terminal stage="roi" />
    </div>
  )
}
