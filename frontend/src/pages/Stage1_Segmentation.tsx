import { useState, useEffect, useRef, useCallback } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runSegmentation, getSegmentationStatus, getSegmentationPreview, runSegmentationPreview, previewPreproc, getRoiSegOverrides, saveRoiSegOverrides, listRois, runFullSegmentation, getFullSegStatus } from '../api/client'
import { useT } from '../i18n'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

// MCseg v2 ROI 個別參數覆寫（欄位與 SegmentationParams 對應）
interface RoiOverride {
  dia_small?: number | null
  dia_mid?: number | null
  dia_large?: number | null
  use_hematoxylin?: boolean | null
  use_cpsam?: boolean | null
  voronoi_distance?: number | null
  flow_threshold?: number | null
  cellprob_threshold?: number | null
  clahe_clip_limit?: number | null
}

// MCseg v2 全域分割參數
interface SegParams {
  use_gpu: boolean
  batch_size: number
  // cyto3 直徑
  dia_small: number
  dia_mid: number
  dia_large: number
  // 可選 pass
  use_hematoxylin: boolean
  use_cpsam: boolean
  // cpsam 7-pass 獨立規格（僅 use_cpsam=true 生效）
  dia_cpsam_auto: number
  dia_cpsam_small: number
  cellprob_cpsam_auto: number
  cellprob_cpsam_small: number
  cellprob_cpsam_hema: number
  // 後處理
  voronoi_distance: number
  min_size: number
  max_size: number
  // Cellpose
  flow_threshold: number
  cellprob_threshold: number
  clahe_clip_limit: number
  // 轉錄本補救
  use_transcript_rescue: boolean
}

const DEFAULT_PARAMS: SegParams = {
  use_gpu: true,
  batch_size: 4,
  dia_small: 13.0,
  dia_mid: 17.0,
  dia_large: 22.0,
  use_hematoxylin: true,
  use_cpsam: false,
  dia_cpsam_auto: 0.0,
  dia_cpsam_small: 16.0,
  cellprob_cpsam_auto: -1.0,
  cellprob_cpsam_small: -3.0,
  cellprob_cpsam_hema: -1.0,
  voronoi_distance: 9,
  min_size: 20,
  max_size: 6000,
  flow_threshold: 0.4,
  cellprob_threshold: -2.0,
  clahe_clip_limit: 3.0,
  use_transcript_rescue: true,
}

function Tooltip({ text }: { text: string }) {
  return (
    <span className="group relative inline-flex items-center ml-1 cursor-help">
      <span className="text-xs text-gray-600 hover:text-gray-400 transition-colors">❓</span>
      <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 w-64 px-3 py-2
                      bg-gray-900 border border-gray-600 rounded-lg text-xs text-gray-300
                      leading-relaxed shadow-2xl pointer-events-none z-50
                      opacity-0 group-hover:opacity-100 transition-opacity whitespace-normal">
        {text}
      </span>
    </span>
  )
}

function NumberInput({
  label, value, onChange, step = 1, min, max, hint, tooltip,
}: {
  label: string; value: number; onChange: (v: number) => void
  step?: number; min?: number; max?: number; hint?: string; tooltip?: string
}) {
  // 本地字串 state：允許鍵盤輸入中間狀態（如 "-1."、"0."）不被 React 重新渲染覆寫
  const [localStr, setLocalStr] = useState(String(value))
  const isEditing = useRef(false)

  useEffect(() => {
    // 只在非編輯狀態（例如快速預設切換）才同步父元件值
    if (!isEditing.current) setLocalStr(String(value))
  }, [value])

  return (
    <div className="flex items-center justify-between gap-3">
      <div className="flex-1 flex items-center">
        <span className="text-sm text-gray-300">{label}</span>
        {hint && <span className="text-xs text-gray-500 ml-1">({hint})</span>}
        {tooltip && <Tooltip text={tooltip} />}
      </div>
      <input
        type="number"
        value={localStr}
        step={step}
        min={min}
        max={max}
        onFocus={() => { isEditing.current = true }}
        onChange={e => {
          setLocalStr(e.target.value)
          const n = parseFloat(e.target.value)
          if (!isNaN(n)) onChange(n)
        }}
        onBlur={() => {
          isEditing.current = false
          // 離開輸入框時若為無效值，還原為父元件值
          if (isNaN(parseFloat(localStr))) setLocalStr(String(value))
        }}
        className="w-24 px-2 py-1 text-sm text-right bg-gray-800 border border-gray-600
                   rounded text-gray-100 focus:outline-none focus:border-blue-500"
      />
    </div>
  )
}

function Toggle({ label, value, onChange, hint, tooltip }: {
  label: string; value: boolean; onChange: (v: boolean) => void; hint?: string; tooltip?: string
}) {
  return (
    <div className="flex items-center justify-between">
      <div className="flex items-center">
        <span className="text-sm text-gray-300">{label}</span>
        {hint && <span className="text-xs text-gray-500 ml-1">({hint})</span>}
        {tooltip && <Tooltip text={tooltip} />}
      </div>
      <button
        onClick={() => onChange(!value)}
        className={`relative w-10 h-5 rounded-full transition-colors ${value ? 'bg-blue-600' : 'bg-gray-600'
          }`}
      >
        <span className={`absolute top-0.5 w-4 h-4 bg-white rounded-full shadow transition-transform ${value ? 'translate-x-5' : 'translate-x-0.5'
          }`} />
      </button>
    </div>
  )
}

function Section({ title, children }: { title: string; children: React.ReactNode }) {
  return (
    <div>
      <h4 className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">{title}</h4>
      <div className="space-y-2.5">{children}</div>
    </div>
  )
}

/** ROI 覆寫表格的數字輸入格——需要本地 string state 才能流暢鍵盤輸入 */
function RoiNumCell({
  val, defaultVal, onUpdate, onClear, step = 0.5, min, max,
}: {
  val: number | null; defaultVal: number
  onUpdate: (v: number) => void; onClear: () => void
  step?: number; min?: number; max?: number
}) {
  const [localStr, setLocalStr] = useState(val != null ? String(val) : '')
  const isEditing = useRef(false)

  useEffect(() => {
    if (!isEditing.current && val != null) setLocalStr(String(val))
    if (val == null) setLocalStr('')
  }, [val])

  if (val == null) {
    return (
      <button
        onClick={() => onUpdate(defaultVal)}
        className="mx-auto block px-2 py-0.5 rounded bg-gray-700/60 text-gray-500 hover:bg-gray-600 hover:text-gray-200 transition-colors"
      >
        全域
      </button>
    )
  }
  return (
    <div className="flex items-center justify-center gap-0.5">
      <input
        type="number"
        value={localStr}
        step={step}
        min={min}
        max={max}
        onFocus={() => { isEditing.current = true }}
        onChange={e => {
          setLocalStr(e.target.value)
          const n = parseFloat(e.target.value)
          if (!isNaN(n)) onUpdate(n)
        }}
        onBlur={() => {
          isEditing.current = false
          if (isNaN(parseFloat(localStr))) setLocalStr(String(val))
        }}
        className="w-16 px-1 py-0.5 bg-gray-800 border border-blue-500 rounded text-gray-100 text-xs text-right focus:outline-none"
      />
      <button
        onClick={onClear}
        className="text-gray-500 hover:text-red-400 leading-none px-0.5"
        title="還原為全域"
      >×</button>
    </div>
  )
}

export default function Stage1_Segmentation() {
  useStageLog('segmentation')
  const { stages, updateStage, rois, setRois } = usePipelineStore()
  const stage = stages['segmentation']
  const { refetch: refetchStatus } = useStageStatus('segmentation', getSegmentationStatus, 3000)
  const t = useT()
  const [params, setParams] = useState<SegParams>(DEFAULT_PARAMS)
  const [showParams, setShowParams] = useState(false)
  const [previewSrc, setPreviewSrc] = useState<string | null>(null)
  const [previewFlows, setPreviewFlows] = useState<string | null>(null)
  const [previewTab, setPreviewTab] = useState<'overlay' | 'flows'>('overlay')
  const [previewRoi, setPreviewRoi] = useState('')
  const [previewAvail, setPreviewAvail] = useState<string[]>([])
  const [previewNCells, setPreviewNCells] = useState<number | null>(null)
  const previewImgRef = useRef<HTMLImageElement>(null)
  const [previewOrigSize, setPreviewOrigSize] = useState<{ w: number; h: number } | null>(null)
  const [previewHover, setPreviewHover] = useState<{
    dx: number; dy: number; ix: number; iy: number; nearRight: boolean; nearBottom: boolean
  } | null>(null)
  const [previewClickMsg, setPreviewClickMsg] = useState('')

  // ── Quick preview state ─────────────────────────────────────────────────
  const [prevRoi, setPrevRoi] = useState('')
  const [prevX, setPrevX] = useState(0)
  const [prevY, setPrevY] = useState(0)
  const [prevPatchSize, setPrevPatchSize] = useState(512)
  const [prevLoading, setPrevLoading] = useState(false)
  const [quickSrc, setQuickSrc] = useState<string | null>(null)
  const [quickClahe, setQuickClahe] = useState<string | null>(null)
  const [quickFlows, setQuickFlows] = useState<string | null>(null)
  const [quickTab, setQuickTab] = useState<'overlay' | 'clahe' | 'flows'>('overlay')
  const [quickInfo, setQuickInfo] = useState<{ n_cells: number; roi_name: string; patch_info: string } | null>(null)
  const [quickError, setQuickError] = useState<string | null>(null)
  const [preprocSrc, setPreprocSrc] = useState<string | null>(null)
  const [preprocInfo, setPreprocInfo] = useState<string | null>(null)
  const [preprocLoading, setPreprocLoading] = useState(false)

  // ── 全圖分割狀態 ──────────────────────────────────────────────────────────
  const [fullSegStatus, setFullSegStatus] = useState<{ status: string; progress?: number; message?: string } | null>(null)
  const fullSegPollRef = useRef<ReturnType<typeof setInterval>>()

  const startFullSegPoll = () => {
    clearInterval(fullSegPollRef.current)
    fullSegPollRef.current = setInterval(async () => {
      try {
        const res = await getFullSegStatus()
        const d = res.data?.data ?? res.data
        setFullSegStatus(d)
        if (d?.status !== 'running') clearInterval(fullSegPollRef.current)
      } catch { clearInterval(fullSegPollRef.current) }
    }, 2000)
  }

  const handleRunFullSeg = async () => {
    setFullSegStatus({ status: 'running', progress: 0, message: '啟動全圖分割...' })
    try {
      await runFullSegmentation()
      startFullSegPoll()
    } catch (e: any) {
      setFullSegStatus({ status: 'error', message: e?.response?.data?.message ?? '啟動失敗' })
    }
  }

  useEffect(() => {
    // 載入時查一次目前狀態
    getFullSegStatus().then(res => {
      const d = res.data?.data ?? res.data
      if (d) {
        setFullSegStatus(d)
        if (d.status === 'running') startFullSegPoll()
      }
    }).catch(() => { })
    return () => clearInterval(fullSegPollRef.current)
  }, [])

  // ── ROI 個別參數覆寫 ────────────────────────────────────────────────────
  const [roiOverrides, setRoiOverrides] = useState<Record<string, RoiOverride>>({})
  const saveTimerRef = useRef<ReturnType<typeof setTimeout>>()

  useEffect(() => {
    getRoiSegOverrides().then(res => {
      if (res.data?.data) setRoiOverrides(res.data.data)
    }).catch(() => { })
    // 重整頁面後 store 清空，自動補載 ROI 清單
    if (rois.length === 0) {
      listRois().then(res => {
        if (res.data?.data) setRois(res.data.data)
      }).catch(() => { })
    }
  }, [])

  const _persistOverrides = (next: Record<string, RoiOverride>) => {
    clearTimeout(saveTimerRef.current)
    saveTimerRef.current = setTimeout(() => {
      saveRoiSegOverrides(next as Record<string, Record<string, unknown>>).catch(() => { })
    }, 600)
  }

  const updateRoiField = (roiName: string, field: keyof RoiOverride, value: unknown) => {
    const next = { ...roiOverrides, [roiName]: { ...roiOverrides[roiName], [field]: value } }
    setRoiOverrides(next)
    _persistOverrides(next)
  }

  const clearRoiField = (roiName: string, field: keyof RoiOverride) => {
    const next = { ...roiOverrides, [roiName]: { ...roiOverrides[roiName], [field]: null } }
    setRoiOverrides(next)
    _persistOverrides(next)
  }

  const resetRoiRow = (roiName: string) => {
    const next = { ...roiOverrides }
    delete next[roiName]
    setRoiOverrides(next)
    saveRoiSegOverrides(next as Record<string, Record<string, unknown>>).catch(() => { })
  }

  const resetAllOverrides = () => {
    setRoiOverrides({})
    saveRoiSegOverrides({}).catch(() => { })
  }

  const set = <K extends keyof SegParams>(key: K, value: SegParams[K]) =>
    setParams(prev => ({ ...prev, [key]: value }))

  const [runningRoi, setRunningRoi] = useState<string | null>(null)  // null = 全部, name = 單一 ROI

  const _buildCleanOverrides = () => Object.fromEntries(
    Object.entries(roiOverrides).filter(([, ov]) => Object.values(ov).some(v => v != null))
  )

  const handleRunAll = async () => {
    updateStage('segmentation', { status: 'running', progress: 0, message: '啟動 MCseg v2（全部 ROI）...' })
    setRunningRoi(null)
    await runSegmentation({ ...params, mode: 'roi', roi_overrides: _buildCleanOverrides() })
    refetchStatus()
  }

  const handleRunSingleRoi = async (roiName: string) => {
    updateStage('segmentation', { status: 'running', progress: 0, message: `啟動 MCseg v2（${roiName}）...` })
    setRunningRoi(roiName)
    await runSegmentation({ ...params, mode: 'roi', roi_overrides: _buildCleanOverrides(), target_roi: roiName })
    refetchStatus()
  }

  const handlePreview = async (roi?: string) => {
    const res = await getSegmentationPreview(roi)
    const d = res.data?.data
    if (d?.image_b64) {
      setPreviewSrc(`data:image/jpeg;base64,${d.image_b64}`)
      setPreviewFlows(d.flows_b64 ? `data:image/jpeg;base64,${d.flows_b64}` : null)
      setPreviewTab('overlay')
      if (d.available_rois) setPreviewAvail(d.available_rois)
      if (d.roi) setPreviewRoi(d.roi)
      if (d.n_cells != null) setPreviewNCells(d.n_cells)
      if (d.orig_w && d.orig_h) setPreviewOrigSize({ w: d.orig_w, h: d.orig_h })
      setPreviewClickMsg('')
    }
  }

  const handlePreviewImgMove = useCallback((e: React.MouseEvent<HTMLDivElement>) => {
    const img = previewImgRef.current
    if (!img || !previewOrigSize) return
    const rect = img.getBoundingClientRect()
    const dx = e.clientX - rect.left
    const dy = e.clientY - rect.top
    setPreviewHover({
      dx, dy,
      ix: Math.round(dx / rect.width * previewOrigSize.w),
      iy: Math.round(dy / rect.height * previewOrigSize.h),
      nearRight: dx > rect.width * 0.65,
      nearBottom: dy > rect.height * 0.75,
    })
  }, [previewOrigSize])

  const handlePreviewImgClick = useCallback(() => {
    if (!previewHover) return
    setPrevRoi(previewRoi)
    setPrevX(previewHover.ix)
    setPrevY(previewHover.iy)
    setPreviewClickMsg(`已選取 (x=${previewHover.ix}, y=${previewHover.iy}) → 快速預覽座標已更新`)
  }, [previewHover, previewRoi])

  // 參數改變時清除過時的快速預覽圖，避免誤以為舊圖是新參數的結果
  useEffect(() => {
    setQuickSrc(null)
    setQuickClahe(null)
    setQuickFlows(null)
    setQuickInfo(null)
    setQuickError(null)
    setPreprocSrc(null)
  }, [params])

  const handleQuickPreview = async () => {
    setPrevLoading(true)
    setQuickError(null)
    setQuickSrc(null)
    setPreprocSrc(null)
    setQuickClahe(null)
    setQuickFlows(null)
    setQuickInfo(null)
    try {
      const res = await runSegmentationPreview({
        roi_name: prevRoi || undefined,
        x: prevX,
        y: prevY,
        patch_size: prevPatchSize,
        use_gpu: params.use_gpu,
        dia_small: params.dia_small,
        dia_mid: params.dia_mid,
        dia_large: params.dia_large,
        use_hematoxylin: params.use_hematoxylin,
        use_cpsam: params.use_cpsam,
        voronoi_distance: params.voronoi_distance,
        flow_threshold: params.flow_threshold,
        cellprob_threshold: params.cellprob_threshold,
        clahe_clip_limit: params.clahe_clip_limit,
      })
      const d = res.data
      if (d?.status === 'ok' && d.data?.image_b64) {
        setQuickSrc(`data:image/jpeg;base64,${d.data.image_b64}`)
        if (d.data.clahe_b64) setQuickClahe(`data:image/jpeg;base64,${d.data.clahe_b64}`)
        if (d.data.flows_b64) setQuickFlows(`data:image/jpeg;base64,${d.data.flows_b64}`)
        setQuickTab('overlay')
        setQuickInfo({ n_cells: d.data.n_cells, roi_name: d.data.roi_name, patch_info: d.data.patch_info })
      } else {
        setQuickError(d?.message ?? '預覽失敗')
      }
    } catch (e: any) {
      setQuickError(e?.response?.data?.message ?? e.message ?? '請求錯誤')
    } finally {
      setPrevLoading(false)
    }
  }

  const handlePreprocPreview = async () => {
    setPreprocLoading(true)
    setQuickError(null)
    setQuickSrc(null)
    setQuickClahe(null)
    setQuickFlows(null)
    setPreprocSrc(null)
    setQuickInfo(null)
    try {
      const res = await previewPreproc({
        roi_name: prevRoi || undefined,
        x: prevX,
        y: prevY,
        patch_size: prevPatchSize,
        clahe_clip_limit: params.clahe_clip_limit,
      })
      const d = res.data
      if (d?.status === 'ok' && d.data?.image_b64) {
        setQuickSrc(`data:image/jpeg;base64,${d.data.image_b64}`)
        if (d.data.clahe_b64) setQuickClahe(`data:image/jpeg;base64,${d.data.clahe_b64}`)
        setQuickTab('clahe') // 前處理預設看 CLAHE 分頁
        setQuickInfo({
          n_cells: 0, // 前處理不會算細胞數
          roi_name: d.data.roi_name || prevRoi || '',
          patch_info: `${d.data.method} · ${d.data.patch_info}`
        })
      } else {
        setQuickError(d?.message ?? '前處理預覽失敗')
      }
    } catch (e: any) {
      setQuickError(e?.response?.data?.message ?? e.message ?? '請求錯誤')
    } finally {
      setPreprocLoading(false)
    }
  }

  return (
    <div className="space-y-4">
      <StageCard
        title={t('stage1.title')}
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
      >
        {/* 參數面板標題列 */}
        <div className="mt-3 flex items-center justify-between">
          <button
            onClick={() => setShowParams(v => !v)}
            className="flex items-center gap-1.5 text-xs text-gray-400 hover:text-gray-200 transition-colors select-none"
          >
            <span className={`transition-transform duration-200 ${showParams ? 'rotate-90' : ''}`}>▶</span>
            <span>{showParams ? t('stage1.hide_params') : t('stage1.show_params')}</span>
          </button>
          {/* 折疊時顯示關鍵參數摘要 */}
          {!showParams && (
            <div className="flex gap-3 text-xs text-gray-500">
              <span>dia {params.dia_small}/{params.dia_mid}/{params.dia_large}px</span>
              <span>Voronoi {params.voronoi_distance}px</span>
              <span>Flow {params.flow_threshold}</span>
              <span className={params.use_cpsam ? 'text-blue-400' : ''}>cpsam {params.use_cpsam ? 'ON' : 'off'}</span>
            </div>
          )}
        </div>

        {/* 參數面板（可折疊）*/}
        {showParams && <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-6">

          {/* 左欄 */}
          <div className="space-y-5">
            <Section title={t('stage1.sec.model')}>
              <Toggle label={t('stage1.param.gpu')} value={params.use_gpu} onChange={v => set('use_gpu', v)}
                tooltip="啟用 CUDA GPU 加速推論。建議開啟；無 GPU 時可關閉改用 CPU（速度約慢 10-20 倍）。" />
              <NumberInput label="Batch Size" value={params.batch_size}
                onChange={v => set('batch_size', v)} min={1} max={16}
                tooltip="GPU 批次大小。GPU 記憶體越大可設越高（建議 4-8）。記憶體不足請降低。" />
            </Section>

            <Section title={t('stage1.sec.diameters')}>
              <NumberInput label={t('stage1.param.dia_small')} value={params.dia_small}
                onChange={v => set('dia_small', v)} step={0.5} min={4} max={40} hint="px"
                tooltip="cyto3 小細胞 pass 的預期直徑。用來補救主 pass 漏掉的小細胞（如淋巴細胞）。預設 13px。" />
              <NumberInput label={t('stage1.param.dia_mid')} value={params.dia_mid}
                onChange={v => set('dia_mid', v)} step={0.5} min={8} max={50} hint="px"
                tooltip="cyto3 主要 pass 的預期細胞直徑（此 pass 結果作為集成基底）。H&E 細胞核通常 15-20px。預設 17px。" />
              <NumberInput label={t('stage1.param.dia_large')} value={params.dia_large}
                onChange={v => set('dia_large', v)} step={0.5} min={12} max={80} hint="px"
                tooltip="cyto3 大細胞 pass，補救大型細胞（如上皮細胞、巨噬細胞）。預設 22px。" />
              <Toggle label={t('stage1.param.hematoxylin')} value={params.use_hematoxylin}
                onChange={v => set('use_hematoxylin', v)}
                tooltip="額外對 Ruifrok H&E 分離的 Hematoxylin 通道跑一次 cyto3（dia=主要直徑）。可補充 H 通道清晰但 RGB 不佳的細胞核。建議開啟。" />
              <Toggle label={t('stage1.param.cpsam')} value={params.use_cpsam}
                onChange={v => set('use_cpsam', v)}
                tooltip="額外加入 Cellpose SAM 模型（cpsam）的 3 個 pass。可提升小型/不規則細胞的召回率，但會顯著增加運算時間（約 2-3 倍）。預設關閉。" />
              {params.use_cpsam && (
                <div className="mt-2 pl-3 border-l-2 border-blue-500/40 space-y-3">
                  <p className="text-xs text-blue-400/80">cpsam 7-pass 進階規格（論文 Pass 5/6/7）</p>
                  <NumberInput label={t('stage1.param.dia_cpsam_auto')} value={params.dia_cpsam_auto}
                    onChange={v => set('dia_cpsam_auto', v)} step={1} min={0} max={60} hint="px"
                    tooltip="Pass 5/7 的 cpsam 直徑。0 = Cellpose 自動偵測（約 30px）。預設 0（auto）。" />
                  <NumberInput label={t('stage1.param.dia_cpsam_small')} value={params.dia_cpsam_small}
                    onChange={v => set('dia_cpsam_small', v)} step={0.5} min={4} max={40} hint="px"
                    tooltip="Pass 6 的 cpsam 固定直徑。預設 16px。" />
                  <NumberInput label={t('stage1.param.cellprob_cpsam_auto')} value={params.cellprob_cpsam_auto}
                    onChange={v => set('cellprob_cpsam_auto', v)} step={0.5} min={-6} max={4}
                    tooltip="Pass 5（CLAHE-RGB, auto dia）的 cellprob 閾值。預設 -1.0。" />
                  <NumberInput label={t('stage1.param.cellprob_cpsam_small')} value={params.cellprob_cpsam_small}
                    onChange={v => set('cellprob_cpsam_small', v)} step={0.5} min={-6} max={4}
                    tooltip="Pass 6（CLAHE-RGB, dia=16）的 cellprob 閾值。預設 -3.0。" />
                  <NumberInput label={t('stage1.param.cellprob_cpsam_hema')} value={params.cellprob_cpsam_hema}
                    onChange={v => set('cellprob_cpsam_hema', v)} step={0.5} min={-6} max={4}
                    tooltip="Pass 7（Hematoxylin, auto dia）的 cellprob 閾值。預設 -1.0。" />
                </div>
              )}
            </Section>
          </div>

          {/* 右欄 */}
          <div className="space-y-5">
            <Section title={t('stage1.sec.voronoi')}>
              <NumberInput label={t('stage1.param.voronoi_dist')} value={params.voronoi_distance}
                onChange={v => set('voronoi_distance', v)} min={3} max={25} hint="px"
                tooltip="Voronoi 擴張的最大距離（像素）。每個細胞向外擴張至多此距離，填補細胞間隙中的 RNA bins。擴張不重疊（Voronoi 性質）。預設 9px ≈ 2.5µm。" />
              <NumberInput label={t('stage1.param.min_size')} value={params.min_size}
                onChange={v => set('min_size', v)} min={5} max={100} hint="px²"
                tooltip="小於此面積的細胞視為雜訊並移除。" />
              <NumberInput label={t('stage1.param.max_size')} value={params.max_size}
                onChange={v => set('max_size', v)} step={500} min={500} max={20000} hint="px²"
                tooltip="大於此面積的細胞（如組織碎片）視為雜訊並移除。" />
              <Toggle label={t('stage1.param.transcript_rescue')} value={params.use_transcript_rescue}
                onChange={v => set('use_transcript_rescue', v)}
                tooltip="從 vhd_pseudo_transcripts.csv 尋找 Cellpose 遺漏的細胞位置（高轉錄本密度但無遮罩的區域）。需要對應的 CSV 檔案存在，若不存在則自動跳過。" />
            </Section>

            <Section title={t('stage1.sec.cellpose_qc')}>
              <NumberInput label={t('stage1.param.flow')} value={params.flow_threshold}
                onChange={v => set('flow_threshold', v)} step={0.05} min={0} max={2}
                tooltip="dP 光流誤差容忍閾值。值越大 → 召回率高但邊界品質低；值越小 → 精確度高。預設 0.4。" />
              <NumberInput label={t('stage1.param.cellprob')} value={params.cellprob_threshold}
                onChange={v => set('cellprob_threshold', v)} step={0.5} min={-6} max={4}
                tooltip="細胞存在機率閾值。值越低（如 -3）→ 更容易偵測細胞（高召回率）；值越高 → 只接受確定的細胞。預設 -2。" />
              <NumberInput label={t('stage1.param.clahe')} value={params.clahe_clip_limit}
                onChange={v => set('clahe_clip_limit', v)} step={0.5} min={0.5} max={8}
                hint="一般組織建議 3.0"
                tooltip="CLAHE 對比度增強強度。值越高 → 邊界更清晰，但雜訊也放大。預設 3.0 適合多數 H&E 場景。" />
            </Section>

            {/* Quick presets */}
            <Section title={t('stage1.sec.presets')}>
              <div className="flex gap-2 flex-wrap">
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS })}
                  className="px-3 py-1 text-xs bg-blue-800 hover:bg-blue-700 rounded text-gray-200"
                >
                  {t('stage1.preset.crc')}
                </button>
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS, dia_small: 10, dia_mid: 14, dia_large: 18, voronoi_distance: 7 })}
                  className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-200"
                >
                  {t('stage1.preset.sparse')}
                </button>
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS, dia_small: 16, dia_mid: 22, dia_large: 30, voronoi_distance: 11, cellprob_threshold: -1.0 })}
                  className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-200"
                >
                  {t('stage1.preset.dense')}
                </button>
              </div>
            </Section>
          </div>
        </div>}
      </StageCard>

      {/* ── ROI 個別參數覆寫 ──────────────────────────────────────────────── */}
      {rois.length > 0 && (
        <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-3">
          <div className="flex items-center justify-between">
            <div>
              <h3 className="text-sm font-semibold text-gray-200">{t('stage1.roi_override')}</h3>
            </div>
            <div className="flex items-center gap-2">
              <button
                onClick={resetAllOverrides}
                className="text-xs text-gray-500 hover:text-gray-300 transition-colors"
              >
                {t('stage1.roi_override.reset_all')}
              </button>
              <button
                onClick={handleRunAll}
                disabled={stage.status === 'running'}
                className="px-4 py-1.5 text-sm rounded-lg font-medium transition-colors
                           bg-brand-primary text-white hover:bg-brand-primary/90
                           disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {stage.status === 'running' && runningRoi === null ? t('common.running') : t('stage1.run_all')}
              </button>
            </div>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="text-gray-500 border-b border-gray-700">
                  <th className="text-left py-2 pr-4 font-medium w-36">ROI</th>
                  <th className="text-center py-2 px-2 font-medium">小徑</th>
                  <th className="text-center py-2 px-2 font-medium">主徑</th>
                  <th className="text-center py-2 px-2 font-medium">大徑</th>
                  <th className="text-center py-2 px-2 font-medium">Voronoi</th>
                  <th className="text-center py-2 px-2 font-medium">Flow</th>
                  <th className="text-center py-2 px-2 font-medium">Cell Prob</th>
                  <th className="text-center py-2 px-2 font-medium">重新分割</th>
                  <th className="py-2 px-2 w-8"></th>
                </tr>
              </thead>
              <tbody>
                {rois.map(roi => {
                  const ov: RoiOverride = roiOverrides[roi.name] ?? {}
                  const hasAny = Object.values(ov).some(v => v != null)

                  return (
                    <tr
                      key={roi.name}
                      className={`border-b border-gray-800/60 ${hasAny ? 'border-l-2 border-l-blue-500' : 'border-l-2 border-l-transparent'}`}
                    >
                      <td className="py-2 pr-4 pl-1">
                        <span className="text-gray-200 font-medium">{roi.name}</span>
                        {roi.tissue && <span className="text-gray-600 ml-1.5 text-xs">{roi.tissue}</span>}
                      </td>

                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.dia_small ?? null} defaultVal={params.dia_small} step={0.5} min={4} max={40}
                          onUpdate={v => updateRoiField(roi.name, 'dia_small', v)}
                          onClear={() => clearRoiField(roi.name, 'dia_small')} />
                      </td>
                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.dia_mid ?? null} defaultVal={params.dia_mid} step={0.5} min={8} max={50}
                          onUpdate={v => updateRoiField(roi.name, 'dia_mid', v)}
                          onClear={() => clearRoiField(roi.name, 'dia_mid')} />
                      </td>
                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.dia_large ?? null} defaultVal={params.dia_large} step={0.5} min={12} max={80}
                          onUpdate={v => updateRoiField(roi.name, 'dia_large', v)}
                          onClear={() => clearRoiField(roi.name, 'dia_large')} />
                      </td>
                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.voronoi_distance ?? null} defaultVal={params.voronoi_distance} step={1} min={3} max={25}
                          onUpdate={v => updateRoiField(roi.name, 'voronoi_distance', v)}
                          onClear={() => clearRoiField(roi.name, 'voronoi_distance')} />
                      </td>
                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.flow_threshold ?? null} defaultVal={params.flow_threshold} step={0.05} min={0} max={2}
                          onUpdate={v => updateRoiField(roi.name, 'flow_threshold', v)}
                          onClear={() => clearRoiField(roi.name, 'flow_threshold')} />
                      </td>
                      <td className="py-2 px-2">
                        <RoiNumCell val={ov.cellprob_threshold ?? null} defaultVal={params.cellprob_threshold} step={0.5} min={-6} max={6}
                          onUpdate={v => updateRoiField(roi.name, 'cellprob_threshold', v)}
                          onClear={() => clearRoiField(roi.name, 'cellprob_threshold')} />
                      </td>
                      {/* 單 ROI 重新分割 */}
                      <td className="py-2 px-2 text-center">
                        <button
                          onClick={() => handleRunSingleRoi(roi.name)}
                          disabled={stage.status === 'running'}
                          title={`只重跑 ${roi.name}`}
                          className={`px-2 py-0.5 rounded text-xs font-medium transition-colors
                            ${stage.status === 'running' && runningRoi === roi.name
                              ? 'bg-blue-800/60 text-blue-300 cursor-not-allowed'
                              : stage.status === 'running'
                                ? 'opacity-40 cursor-not-allowed bg-gray-700 text-gray-400'
                                : 'bg-blue-700/40 border border-blue-600 text-blue-300 hover:bg-blue-600/60'
                            }`}
                        >
                          {stage.status === 'running' && runningRoi === roi.name ? t('common.running') : t('stage1.run_single')}
                        </button>
                      </td>

                      {/* 重置整列 */}
                      <td className="py-2 px-2">
                        {hasAny && (
                          <button
                            onClick={() => resetRoiRow(roi.name)}
                            className="text-gray-600 hover:text-red-400 transition-colors text-xs"
                            title="重置此 ROI 所有覆寫"
                          >
                            重置
                          </button>
                        )}
                      </td>
                    </tr>
                  )
                })}
              </tbody>
            </table>
          </div>

          <p className="text-xs text-gray-600">
            ⓘ 點擊「全域」啟用個別設定（顯示藍框輸入）；點擊「×」還原；有覆寫的 ROI 列顯示藍色左邊線
          </p>
        </div>
      )}

      {/* ── 快速 Patch 預覽 ─────────────────────────────────────────────────── */}
      <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-4">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-sm font-semibold text-gray-200">{t('stage1.quick_preview')}</h3>
          </div>
          <div className="flex gap-2">
            <button
              onClick={handlePreprocPreview}
              disabled={preprocLoading}
              className="px-4 py-1.5 text-sm rounded bg-gray-600 text-gray-200 font-medium
                         hover:bg-gray-500 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              {preprocLoading ? t('common.loading') : `⚡ ${t('stage1.quick_preview.preproc')}`}
            </button>
            <button
              onClick={handleQuickPreview}
              disabled={prevLoading}
              className="px-4 py-1.5 text-sm rounded bg-primary text-black font-medium
                         hover:bg-primary/90 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              {prevLoading ? t('common.running') : `🔬 ${t('stage1.quick_preview.run')}`}
            </button>
          </div>
        </div>

        {/* 參數列 */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
          <div>
            <label className="text-xs text-gray-400 mb-1 block">ROI</label>
            <select
              value={prevRoi}
              onChange={e => setPrevRoi(e.target.value)}
              className="w-full px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                         text-gray-100 focus:outline-none focus:border-blue-500"
            >
              <option value="">自動選取</option>
              {rois.filter(r => r.name).map(r => (
                <option key={r.name} value={r.name}>{r.name}</option>
              ))}
            </select>
          </div>

          <div>
            <label className="text-xs text-gray-400 mb-1 block">{t('stage1.quick_preview.x')}</label>
            <input
              type="number" value={prevX} min={0}
              onChange={e => setPrevX(Math.max(0, parseInt(e.target.value) || 0))}
              className="w-full px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                         text-gray-100 text-right focus:outline-none focus:border-blue-500"
            />
          </div>

          <div>
            <label className="text-xs text-gray-400 mb-1 block">{t('stage1.quick_preview.y')}</label>
            <input
              type="number" value={prevY} min={0}
              onChange={e => setPrevY(Math.max(0, parseInt(e.target.value) || 0))}
              className="w-full px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                         text-gray-100 text-right focus:outline-none focus:border-blue-500"
            />
          </div>

          <div>
            <label className="text-xs text-gray-400 mb-1 block">Patch Size</label>
            <select
              value={prevPatchSize}
              onChange={e => setPrevPatchSize(parseInt(e.target.value))}
              className="w-full px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                         text-gray-100 focus:outline-none focus:border-blue-500"
            >
              <option value={256}>256 × 256</option>
              <option value={512}>512 × 512</option>
              <option value={1024}>1024 × 1024</option>
            </select>
          </div>
        </div>

        {/* 錯誤訊息 */}
        {quickError && (
          <p className="text-xs text-red-400 bg-red-900/20 rounded px-3 py-2">{quickError}</p>
        )}

        {/* 預覽結果 */}
        {quickSrc && (
          <div className="space-y-2">
            {quickInfo && (
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-4 text-xs text-gray-400">
                  {quickInfo.n_cells > 0 && (
                    <span className="text-green-400 font-medium">✓ {quickInfo.n_cells} {t('stage1.preview.cells')}</span>
                  )}
                  <span>ROI: {quickInfo.roi_name}</span>
                  <span>{quickInfo.patch_info}</span>
                </div>
                {/* 圖層切換 Tab */}
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  {(['overlay', 'clahe', 'flows'] as const).map(tab => {
                    const labels: Record<typeof tab, string> = {
                      overlay: 'H&E + 邊界',
                      clahe: 'CLAHE 前處理',
                      flows: 'Flow 方向圖'
                    }
                    const available = tab === 'overlay' || (tab === 'clahe' && !!quickClahe) || (tab === 'flows' && !!quickFlows)
                    return (
                      <button
                        key={tab}
                        disabled={!available}
                        onClick={() => setQuickTab(tab)}
                        className={`px-3 py-1 transition-colors disabled:opacity-30 disabled:cursor-not-allowed ${quickTab === tab
                          ? 'bg-blue-600 text-white'
                          : 'bg-gray-800 text-gray-400 hover:bg-gray-700'
                          }`}
                      >
                        {labels[tab]}
                      </button>
                    )
                  })}
                </div>
              </div>
            )}
            {/* 說明文字 */}
            <p className="text-xs text-gray-500">
              {quickTab === 'overlay' && 'H&E 原圖 + 綠色細胞邊界（MCseg v2 Voronoi 集成）'}
              {quickTab === 'clahe' && 'CLAHE 局部對比增強（Cellpose 實際輸入）'}
              {quickTab === 'flows' && 'Cellpose 小尺寸 dP 光流方向圖（色相 = 方向，飽和度 = 強度）；白線 = 細胞邊界'}
            </p>
            <div className="rounded-lg overflow-hidden border border-surface-border"
              style={{ imageRendering: 'pixelated' }}>
              <img
                src={
                  quickTab === 'clahe' ? (quickClahe ?? quickSrc) :
                    quickTab === 'flows' ? (quickFlows ?? quickSrc) :
                      quickSrc
                }
                alt={quickTab}
                className="w-full"
              />
            </div>
          </div>
        )}
      </div>

      {/* 完整分割後的遮罩預覽 */}
      {stage.status === 'done' && (
        <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-3">
          <div className="flex items-center justify-between">
            <div>
              <h3 className="text-sm font-semibold text-gray-200">{t('stage1.preview.title')}</h3>
            </div>
            <div className="flex items-center gap-2">
              {previewAvail.length > 0 && (
                <select
                  value={previewRoi}
                  onChange={e => handlePreview(e.target.value)}
                  className="px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                             text-gray-100 focus:outline-none focus:border-blue-500"
                >
                  {previewAvail.map(r => <option key={r} value={r}>{r}</option>)}
                </select>
              )}
              <button
                onClick={() => handlePreview(previewRoi || undefined)}
                className="px-4 py-1.5 text-sm rounded bg-gray-700 hover:bg-gray-600 text-gray-200"
              >
                {previewSrc ? t('stage1.preview.reload') : t('stage1.preview.load')}
              </button>
            </div>
          </div>
          {previewSrc && (
            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <span className="text-xs text-green-400 font-medium">
                  {previewNCells != null && `✓ ${previewNCells.toLocaleString()} ${t('stage1.preview.cells')} · ROI: ${previewRoi}`}
                </span>
                {/* Tab 切換 */}
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  <button
                    onClick={() => setPreviewTab('overlay')}
                    className={`px-3 py-1 transition-colors ${previewTab === 'overlay' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
                  >{t('stage1.preview.overlay')}</button>
                  <button
                    disabled={!previewFlows}
                    onClick={() => setPreviewTab('flows')}
                    className={`px-3 py-1 transition-colors disabled:opacity-30 disabled:cursor-not-allowed ${previewTab === 'flows' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
                  >{t('stage1.preview.flows')}</button>
                </div>
              </div>
              <p className="text-xs text-gray-500">
                {previewTab === 'overlay' && 'H&E 原圖 + 綠色細胞邊界（來自已存遮罩）'}
                {previewTab === 'flows' && 'Cellpose 小尺寸 dP 光流方向圖（色相 = 方向，飽和度 = 強度）'}
              </p>
              {/* 互動預覽圖：懸停顯示座標，點擊填入快速測試 */}
              <div
                className="relative rounded-lg overflow-hidden border border-surface-border cursor-crosshair select-none"
                onMouseMove={handlePreviewImgMove}
                onMouseLeave={() => setPreviewHover(null)}
                onClick={handlePreviewImgClick}
              >
                <img
                  ref={previewImgRef}
                  src={previewTab === 'flows' ? (previewFlows ?? previewSrc!) : previewSrc!}
                  alt="segmentation preview"
                  className="w-full block"
                />
                {/* 十字線 */}
                {previewHover && (
                  <>
                    <div className="absolute top-0 bottom-0 w-px bg-yellow-400/70 pointer-events-none"
                      style={{ left: previewHover.dx }} />
                    <div className="absolute left-0 right-0 h-px bg-yellow-400/70 pointer-events-none"
                      style={{ top: previewHover.dy }} />
                    {/* 座標 badge */}
                    <div
                      className="absolute bg-black/85 text-yellow-300 text-xs font-mono px-2 py-1 rounded pointer-events-none whitespace-nowrap z-10"
                      style={{
                        left: previewHover.nearRight ? previewHover.dx - 8 : previewHover.dx + 8,
                        top: previewHover.nearBottom ? previewHover.dy - 30 : previewHover.dy + 8,
                        transform: previewHover.nearRight ? 'translateX(-100%)' : undefined,
                      }}
                    >
                      x={previewHover.ix}, y={previewHover.iy}
                    </div>
                  </>
                )}
                {/* 提示：未懸停時 */}
                {!previewHover && previewOrigSize && (
                  <div className="absolute bottom-2 left-1/2 -translate-x-1/2 bg-black/65 text-gray-400 text-xs px-3 py-1 rounded-full pointer-events-none whitespace-nowrap">
                    懸停顯示座標 · 點擊選取 → 小片段測試
                  </div>
                )}
              </div>
              {/* 點擊後回饋訊息 */}
              {previewClickMsg && (
                <p className="text-xs text-yellow-400 bg-yellow-900/20 rounded px-3 py-1.5">
                  ✓ {previewClickMsg}，可至上方「快速 Patch 預覽」調整參數後執行測試
                </p>
              )}
            </div>
          )}
        </div>
      )}

      {/* ── 全圖分割區塊 ─────────────────────────────────────────────────── */}
      <div className="rounded-xl border border-surface-border bg-surface-card p-4 space-y-3">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-sm font-semibold text-gray-200">{t('stage1.full_seg.title')}</h3>
            <p className="text-xs text-gray-500 mt-0.5">
              {t('stage1.full_seg.description')}
              <code className="mx-1 text-yellow-400">full_image_segmentation_masks.npy</code>
            </p>
          </div>
          <button
            onClick={handleRunFullSeg}
            disabled={fullSegStatus?.status === 'running'}
            className="shrink-0 px-4 py-2 text-sm rounded-lg font-medium transition-colors
                       bg-indigo-600 hover:bg-indigo-500 disabled:opacity-50 disabled:cursor-not-allowed text-white"
          >
            {fullSegStatus?.status === 'running' ? t('stage1.full_seg.running') : t('stage1.full_seg.run')}
          </button>
        </div>

        {fullSegStatus && (
          <div className={`rounded-lg px-3 py-2 text-xs font-mono space-y-1
            ${fullSegStatus.status === 'error' ? 'bg-red-900/30 text-red-300 border border-red-800'
              : fullSegStatus.status === 'done' ? 'bg-green-900/30 text-green-300 border border-green-800'
              : 'bg-gray-800 text-gray-300 border border-gray-700'}`}
          >
            <div className="flex items-center justify-between">
              <span>
                {fullSegStatus.status === 'running' && '⏳ '}
                {fullSegStatus.status === 'done' && '✓ '}
                {fullSegStatus.status === 'error' && '✗ '}
                {fullSegStatus.message ?? fullSegStatus.status}
              </span>
              {fullSegStatus.progress != null && fullSegStatus.status === 'running' && (
                <span className="text-gray-400">{Math.round(fullSegStatus.progress * 100)}%</span>
              )}
            </div>
            {fullSegStatus.status === 'running' && fullSegStatus.progress != null && (
              <div className="w-full bg-gray-700 rounded-full h-1.5 overflow-hidden">
                <div
                  className="bg-indigo-500 h-full rounded-full transition-all duration-500"
                  style={{ width: `${Math.round(fullSegStatus.progress * 100)}%` }}
                />
              </div>
            )}
          </div>
        )}
      </div>

      <Terminal stage="segmentation" />
    </div>
  )
}
