import { useEffect, useState } from 'react'
import {
  getRegistrationEstimate,
  makeRegistrationQcPatches,
  getRegistrationQcImages,
  applyRegistration,
  getAlignmentJson,
  setAlignmentJson,
} from '../../api/client'
import { useT } from '../../i18n'

/**
 * 對位檢查面板：QC 疊圖 + 殘餘位移估計 + 套用變換。
 *
 * 「套用」刻意分成兩步（先估計、看過 residual 再確認）—— 估計錯誤而被靜默
 * 套用，比完全不修正更糟。
 */

interface Estimate {
  dy: number
  dx: number
  spread: number | null
  transform_source: string
}

interface QcImage {
  name: string
  data: string
}

export default function AlignmentPanel() {
  const t = useT()
  const [estimate, setEstimate] = useState<Estimate | null>(null)
  const [images, setImages] = useState<QcImage[]>([])
  const [busy, setBusy] = useState<'' | 'estimate' | 'qc' | 'apply'>('')
  const [error, setError] = useState('')
  const [applied, setApplied] = useState('')
  const [pending, setPending] = useState<{ matrix: number[][]; estimated_error: number | null } | null>(null)

  // 對位 JSON（Loupe 重新對位）路徑
  const [jsonPath, setJsonPath] = useState('')
  const [currentJsonPath, setCurrentJsonPath] = useState<string | null>(null)
  const [jsonBusy, setJsonBusy] = useState(false)
  const [jsonError, setJsonError] = useState('')
  const [jsonResult, setJsonResult] = useState('')

  useEffect(() => {
    getAlignmentJson()
      .then(r => {
        const p = r.data?.data?.path ?? null
        setCurrentJsonPath(p)
        setJsonPath(p ?? '')
      })
      .catch(() => {})   // B 類讀取降級：讀不到只是欄位空白，不擋人
  }, [])

  const handleSetJson = async () => {
    if (!jsonPath.trim()) return
    setJsonBusy(true)
    setJsonError('')
    setJsonResult('')
    try {
      const r = await setAlignmentJson(jsonPath.trim())
      if (r.data.status !== 'ok') { setJsonError(r.data.message ?? 'Error'); return }
      setCurrentJsonPath(jsonPath.trim())
      setJsonResult(r.data.message ?? '')
    } catch (e: any) {
      setJsonError(e.response?.data?.detail ?? e.message ?? 'Unknown error')
    } finally {
      setJsonBusy(false)
    }
  }

  const handleClearJson = async () => {
    setJsonBusy(true)
    setJsonError('')
    setJsonResult('')
    try {
      const r = await setAlignmentJson('')
      if (r.data.status !== 'ok') { setJsonError(r.data.message ?? 'Error'); return }
      setCurrentJsonPath(null)
      setJsonPath('')
      setJsonResult(t('align.json_clear'))
    } catch (e: any) {
      setJsonError(e.response?.data?.detail ?? e.message ?? 'Unknown error')
    } finally {
      setJsonBusy(false)
    }
  }


  const run = async (kind: 'estimate' | 'qc' | 'apply', enable = false) => {
    setBusy(kind)
    setError('')
    try {
      if (kind === 'estimate') {
        const r = await getRegistrationEstimate()
        if (r.data.status === 'ok') setEstimate(r.data.data)
        else setError(r.data.message ?? 'Error')
      } else if (kind === 'qc') {
        const made = await makeRegistrationQcPatches({ n: 3, size: 512 })
        if (made.data.status !== 'ok') { setError(made.data.message ?? 'Error'); return }
        const got = await getRegistrationQcImages()
        setImages(got.data.data ?? [])
      } else {
        const r = await applyRegistration({ enable })
        if (r.data.status !== 'ok') { setError(r.data.message ?? 'Error'); return }
        setPending({ matrix: r.data.data.matrix, estimated_error: r.data.data.estimated_error })
        setApplied(r.data.message ?? '')
      }
    } catch (e: any) {
      setError(e.response?.data?.detail ?? e.message ?? 'Unknown error')
    } finally {
      setBusy('')
    }
  }

  const btn = 'px-3 py-1.5 rounded text-xs font-medium border border-surface-border text-gray-300 hover:border-primary hover:text-primary transition-colors disabled:opacity-40'

  return (
    <div className="space-y-3">
      <div className="bg-surface/50 rounded px-3 py-2 space-y-2">
        <p className="text-xs text-gray-400 font-medium">{t('align.json_title')}</p>
        <p className="text-xs text-gray-500">{t('align.json_subtitle')}</p>
        <p className="text-xs text-gray-500">
          {t('align.json_current')}：{' '}
          <span className="text-gray-300 font-mono">{currentJsonPath ?? t('align.json_none')}</span>
        </p>
        <div className="flex items-center gap-2">
          <input
            type="text"
            value={jsonPath}
            onChange={e => { setJsonPath(e.target.value); setJsonError(''); setJsonResult('') }}
            placeholder={t('align.json_placeholder')}
            className="flex-1 px-2 py-1.5 bg-surface border border-surface-border rounded text-xs font-mono text-gray-200 focus:border-primary focus:outline-none"
          />
          <button
            className={btn}
            disabled={jsonBusy || !jsonPath.trim()}
            onClick={() => void handleSetJson()}
          >
            {jsonBusy ? t('align.json_applying') : t('align.json_apply')}
          </button>
          {(currentJsonPath || jsonPath) && (
            <button
              className={btn}
              disabled={jsonBusy}
              onClick={() => void handleClearJson()}
            >
              {t('align.json_clear')}
            </button>
          )}
        </div>

        {jsonError && <p className="text-xs text-red-400">{jsonError}</p>}
        {jsonResult && <p className="text-xs text-green-400">{jsonResult}</p>}
      </div>

      <div className="flex items-center gap-2 flex-wrap">
        <button className={btn} disabled={busy !== ''} onClick={() => run('estimate')}>
          {busy === 'estimate' ? t('align.estimating') : t('align.estimate')}
        </button>
        <button className={btn} disabled={busy !== ''} onClick={() => run('qc')}>
          {busy === 'qc' ? t('align.rendering') : t('align.qc_patches')}
        </button>
        <button className={btn} disabled={busy !== ''} onClick={() => run('apply', false)}>
          {t('align.compute_affine')}
        </button>
      </div>

      {error && <p className="text-xs text-red-400">{error}</p>}

      {estimate && (
        <div className="text-xs text-gray-400 bg-surface/50 rounded px-3 py-2 space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>dy: <b className="text-gray-200">{estimate.dy.toFixed(1)} px</b></span>
            <span>dx: <b className="text-gray-200">{estimate.dx.toFixed(1)} px</b></span>
            <span>
              spread:{' '}
              <b className={estimate.spread != null && estimate.spread < 3 ? 'text-green-400' : 'text-yellow-400'}>
                {estimate.spread == null ? 'n/a' : `${estimate.spread.toFixed(1)} px`}
              </b>
            </span>
          </div>
          <p className="text-gray-500">{estimate.transform_source}</p>
          <p className="text-gray-500">{t('align.spread_hint')}</p>
        </div>
      )}

      {pending && (
        <div className="text-xs text-gray-400 bg-surface/50 rounded px-3 py-2 space-y-2">
          <p className="font-mono text-[11px] text-gray-300">
            [[{pending.matrix[0].map(v => v.toFixed(4)).join(', ')}],<br />
            &nbsp;[{pending.matrix[1].map(v => v.toFixed(4)).join(', ')}]]
          </p>
          <p>
            residual:{' '}
            <b className="text-gray-200">
              {pending.estimated_error == null ? 'n/a' : pending.estimated_error.toFixed(2)}
            </b>
          </p>
          <p className="text-yellow-400">{applied}</p>
          <button
            className={btn}
            disabled={busy !== ''}
            onClick={() => {
              // run 內部已有 try/catch → setError
              if (window.confirm(t('align.confirm_apply'))) void run('apply', true)
            }}
          >
            {t('align.apply')}
          </button>
        </div>
      )}

      {images.length > 0 && (
        <div className="grid grid-cols-3 gap-2">
          {images.map(im => (
            <img key={im.name} src={im.data} alt={im.name} className="w-full rounded border border-surface-border" />
          ))}
        </div>
      )}
    </div>
  )
}
