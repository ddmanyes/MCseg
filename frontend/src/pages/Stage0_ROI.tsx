import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { listRois, addRoi, deleteRoi, runRoiExtract, getRoiStatus } from '../api/client'
import type { RoiDefinition } from '../types/pipeline'
import useStageLog from '../hooks/useStageLog'
import RoiSelector from '../components/roi/RoiSelector'
import { useStageStatus } from '../hooks/useStageStatus'
import { useT } from '../i18n'

export default function Stage0_ROI() {
  useStageLog('roi')
  const { stages, updateStage, rois, setRois } = usePipelineStore()
  const stage = stages['roi']
  const { refetch: refetchStatus } = useStageStatus('roi', getRoiStatus, 2000)
  const [form, setForm] = useState<Partial<RoiDefinition>>({ pixel_size_um: 0.2737 })
  const [formError, setFormError] = useState<string | null>(null)
  const t = useT()

  useEffect(() => {
    listRois().then(r => setRois(r.data.data ?? []))
  }, [])

  const handleRun = async () => {
    updateStage('roi', { status: 'running', progress: 0, message: t('stage0.running') })
    try {
      await runRoiExtract()
      refetchStatus()
    } catch (e: any) {
      updateStage('roi', { status: 'error', message: e.message })
    }
  }

  const handleAdd = async () => {
    if (!form.name) { setFormError(t('stage0.form.name') + ' required'); return }
    if (!form.tissue) { setFormError(t('stage0.form.tissue') + ' required'); return }
    setFormError(null)
    await addRoi(form as RoiDefinition)
    const updated = await listRois()
    setRois(updated.data.data ?? [])
    setForm({ pixel_size_um: 0.2737 })
  }

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
                onClick={() => deleteRoi(roi.name).then(() => listRois().then(r => setRois(r.data.data ?? [])))}
                className="text-red-400 hover:text-red-300 text-xs"
              >
                {t('common.delete')}
              </button>
            </div>
          ))}
        </div>

        {/* Interactive ROI selector */}
        <div className="border-t border-surface-border pt-4">
          <p className="text-xs text-gray-400 font-medium uppercase tracking-wide mb-3">{t('stage0.interactive')}</p>
          <RoiSelector
            existingRois={rois as any}
            onSelect={(roi) => setForm(f => ({ ...f, ...roi }))}
          />
        </div>

        {/* Add ROI form */}
        <div className="border-t border-surface-border pt-4">
          <p className="text-xs text-gray-400 font-medium uppercase tracking-wide mb-3">{t('stage0.add_roi')}</p>
          <div className="grid grid-cols-2 gap-3">
            {([
              ['name',      t('stage0.form.name'),   'text'],
              ['tissue',    t('stage0.form.tissue'),  'text'],
              ['x',         t('stage0.form.x'),       'number'],
              ['y',         t('stage0.form.y'),       'number'],
              ['width_px',  t('stage0.form.width'),   'number'],
              ['height_px', t('stage0.form.height'),  'number'],
            ] as [string, string, string][]).map(([key, label, type]) => (
              <div key={key}>
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
          <button
            onClick={handleAdd}
            className="mt-3 px-4 py-1.5 bg-surface-border hover:bg-surface-border/80 rounded text-sm text-gray-200 transition-colors"
          >
            + {t('stage0.add_roi')}
          </button>
          {formError && <p className="mt-2 text-xs text-red-400">{formError}</p>}
        </div>
      </StageCard>

      <Terminal stage="roi" />
    </div>
  )
}
