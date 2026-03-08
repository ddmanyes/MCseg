import { useState, useEffect, useRef } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runSegmentation, getSegmentationStatus, getSegmentationPreview, runSegmentationPreview, previewPreproc, getRoiSegOverrides, saveRoiSegOverrides } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

interface RoiOverride {
  model_type?: string | null
  dia_small?: number | null
  dia_large?: number | null
  flow_threshold?: number | null
  cellprob_threshold?: number | null
  fragment_threshold?: number | null
  eosin_bg_threshold?: number | null
}

interface SegParams {
  mode: 'roi' | 'full'
  model_type: 'cyto2' | 'cyto3' | 'nuclei'
  use_gpu: boolean
  batch_size: number
  dia_small: number
  dia_large: number
  flow_threshold: number
  cellprob_threshold: number
  fragment_threshold: number
  normalize_stains: boolean
  clahe_clip_limit: number
  enable_eosin_watershed: boolean
  eosin_bg_threshold: number
  block_size: number
  overlap: number
}

const DEFAULT_PARAMS: SegParams = {
  mode: 'roi',
  model_type: 'cyto2',
  use_gpu: true,
  batch_size: 4,
  dia_small: 30.0,
  dia_large: 60.0,
  flow_threshold: 0.4,
  cellprob_threshold: -1.0,
  fragment_threshold: 200,
  normalize_stains: true,
  clahe_clip_limit: 1.0,
  enable_eosin_watershed: true,
  eosin_bg_threshold: 80,
  block_size: 2048,
  overlap: 256,
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
  return (
    <div className="flex items-center justify-between gap-3">
      <div className="flex-1 flex items-center">
        <span className="text-sm text-gray-300">{label}</span>
        {hint && <span className="text-xs text-gray-500 ml-1">({hint})</span>}
        {tooltip && <Tooltip text={tooltip} />}
      </div>
      <input
        type="number"
        value={value}
        step={step}
        min={min}
        max={max}
        onChange={e => onChange(parseFloat(e.target.value))}
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

export default function Stage1_Segmentation() {
  useStageLog('segmentation')
  const { stages, updateStage, rois } = usePipelineStore()
  const stage = stages['segmentation']
  const { refetch: refetchStatus } = useStageStatus('segmentation', getSegmentationStatus, 3000)
  const [params, setParams] = useState<SegParams>(DEFAULT_PARAMS)
  const [previewSrc, setPreviewSrc] = useState<string | null>(null)
  const [previewCyto, setPreviewCyto] = useState<string | null>(null)
  const [previewFlows, setPreviewFlows] = useState<string | null>(null)
  const [previewTab, setPreviewTab] = useState<'overlay' | 'cyto' | 'flows'>('overlay')
  const [previewRoi, setPreviewRoi] = useState('')
  const [previewAvail, setPreviewAvail] = useState<string[]>([])
  const [previewNCells, setPreviewNCells] = useState<number | null>(null)

  // ── Quick preview state ─────────────────────────────────────────────────
  const [prevRoi, setPrevRoi] = useState('')
  const [prevX, setPrevX] = useState(0)
  const [prevY, setPrevY] = useState(0)
  const [prevPatchSize, setPrevPatchSize] = useState(512)
  const [prevLoading, setPrevLoading] = useState(false)
  const [quickSrc, setQuickSrc] = useState<string | null>(null)
  const [quickMacenko, setQuickMacenko] = useState<string | null>(null)
  const [quickFlows, setQuickFlows] = useState<string | null>(null)
  const [quickCyto, setQuickCyto] = useState<string | null>(null) // 新增 cyto state
  const [quickTab, setQuickTab] = useState<'overlay' | 'macenko' | 'flows' | 'cyto'>('overlay') // 加入 cyto type
  const [quickInfo, setQuickInfo] = useState<{ n_cells: number; roi_name: string; patch_info: string } | null>(null)
  const [quickError, setQuickError] = useState<string | null>(null)
  const [preprocSrc, setPreprocSrc] = useState<string | null>(null)
  const [preprocInfo, setPreprocInfo] = useState<string | null>(null)
  const [preprocLoading, setPreprocLoading] = useState(false)

  // ── ROI 個別參數覆寫 ────────────────────────────────────────────────────
  const [roiOverrides, setRoiOverrides] = useState<Record<string, RoiOverride>>({})
  const saveTimerRef = useRef<ReturnType<typeof setTimeout>>()

  useEffect(() => {
    getRoiSegOverrides().then(res => {
      if (res.data?.data) setRoiOverrides(res.data.data)
    }).catch(() => { })
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
    updateStage('segmentation', { status: 'running', progress: 0, message: '啟動 Cellpose（全部 ROI）...' })
    setRunningRoi(null)
    await runSegmentation({ ...params, roi_overrides: _buildCleanOverrides() })
    refetchStatus()
  }

  const handleRunSingleRoi = async (roiName: string) => {
    updateStage('segmentation', { status: 'running', progress: 0, message: `啟動 Cellpose（${roiName}）...` })
    setRunningRoi(roiName)
    await runSegmentation({ ...params, roi_overrides: _buildCleanOverrides(), target_roi: roiName })
    refetchStatus()
  }

  const handlePreview = async (roi?: string) => {
    const res = await getSegmentationPreview(roi)
    const d = res.data?.data
    if (d?.image_b64) {
      setPreviewSrc(`data:image/jpeg;base64,${d.image_b64}`)
      setPreviewCyto(d.cyto_b64 ? `data:image/jpeg;base64,${d.cyto_b64}` : null)
      setPreviewFlows(d.flows_b64 ? `data:image/jpeg;base64,${d.flows_b64}` : null)
      setPreviewTab('overlay')
      if (d.available_rois) setPreviewAvail(d.available_rois)
      if (d.roi) setPreviewRoi(d.roi)
      if (d.n_cells != null) setPreviewNCells(d.n_cells)
    }
  }

  const handleQuickPreview = async () => {
    setPrevLoading(true)
    setQuickError(null)
    setQuickSrc(null)
    setPreprocSrc(null)
    setQuickMacenko(null)
    setQuickFlows(null)
    setQuickCyto(null)
    setQuickInfo(null)
    try {
      const res = await runSegmentationPreview({
        roi_name: prevRoi || undefined,
        x: prevX,
        y: prevY,
        patch_size: prevPatchSize,
        model_type: params.model_type,
        use_gpu: params.use_gpu,
        dia_small: params.dia_small,
        dia_large: params.dia_large,
        flow_threshold: params.flow_threshold,
        cellprob_threshold: params.cellprob_threshold,
        fragment_threshold: params.fragment_threshold,
        normalize_stains: params.normalize_stains,
        clahe_clip_limit: params.clahe_clip_limit,
        enable_eosin_watershed: params.enable_eosin_watershed,
        eosin_bg_threshold: params.eosin_bg_threshold
      })
      const d = res.data
      if (d?.status === 'ok' && d.data?.image_b64) {
        setQuickSrc(`data:image/jpeg;base64,${d.data.image_b64}`)
        if (d.data.macenko_b64) setQuickMacenko(`data:image/jpeg;base64,${d.data.macenko_b64}`)
        if (d.data.flows_b64) setQuickFlows(`data:image/jpeg;base64,${d.data.flows_b64}`)
        if (d.data.cyto_b64) setQuickCyto(`data:image/jpeg;base64,${d.data.cyto_b64}`)
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
    setQuickMacenko(null)
    setQuickFlows(null)
    setQuickCyto(null)
    setPreprocSrc(null)
    setQuickInfo(null)
    try {
      const res = await previewPreproc({
        roi_name: prevRoi || undefined,
        x: prevX,
        y: prevY,
        patch_size: prevPatchSize,
        normalize_stains: params.normalize_stains,
        clahe_clip_limit: params.clahe_clip_limit,
        enable_eosin_watershed: params.enable_eosin_watershed,
        eosin_bg_threshold: params.eosin_bg_threshold
      })
      const d = res.data
      if (d?.status === 'ok' && d.data?.image_b64) {
        setQuickSrc(`data:image/jpeg;base64,${d.data.image_b64}`)
        if (d.data.macenko_b64) setQuickMacenko(`data:image/jpeg;base64,${d.data.macenko_b64}`)
        if (d.data.cyto_b64) setQuickCyto(`data:image/jpeg;base64,${d.data.cyto_b64}`)
        setQuickTab('macenko') // 前處理預設看 macenko 分頁
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
        title="細胞分割（Cellpose + Logic A）"
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
      >
        {/* 參數面板 */}
        <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-6">

          {/* 左欄 */}
          <div className="space-y-5">
            <Section title="模型設定">
              <div className="flex items-center justify-between">
                <div className="flex items-center">
                  <span className="text-sm text-gray-300">執行範圍</span>
                  <span className="text-xs text-gray-500 ml-1">
                    {params.mode === 'roi' ? '（各 ROI 裁切圖）' : '（完整 BTF，耗時）'}
                  </span>
                  <Tooltip text="ROI 模式：對每張已裁切的 he_crop.tif 分別執行，速度快、記憶體低。全圖模式：對整張 BTF 大圖執行，需大量記憶體，耗時較長。" />
                </div>
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  <button
                    onClick={() => set('mode', 'roi')}
                    className={`px-3 py-1 transition-colors ${params.mode === 'roi' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'
                      }`}
                  >ROI</button>
                  <button
                    onClick={() => set('mode', 'full')}
                    className={`px-3 py-1 transition-colors ${params.mode === 'full' ? 'bg-orange-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'
                      }`}
                  >全圖</button>
                </div>
              </div>
              <div className="flex items-center justify-between">
                <div className="flex items-center">
                  <span className="text-sm text-gray-300">模型類型</span>
                  <Tooltip text="cyto2：適合上皮細胞、狹長形細胞（如大腸癌）。cyto3：Cellpose 3 通用模型，廣泛場景適用。nuclei：只偵測細胞核，適合核仁清晰的場景。" />
                </div>
                <select
                  value={params.model_type}
                  onChange={e => set('model_type', e.target.value as SegParams['model_type'])}
                  className="px-2 py-1 text-sm bg-gray-800 border border-gray-600 rounded
                             text-gray-100 focus:outline-none focus:border-blue-500"
                >
                  <option value="cyto2">cyto2（上皮/狹長細胞）</option>
                  <option value="cyto3">cyto3（通用）</option>
                  <option value="nuclei">nuclei（細胞核）</option>
                </select>
              </div>
              <Toggle label="GPU 加速" value={params.use_gpu} onChange={v => set('use_gpu', v)}
                tooltip="啟用 CUDA GPU 加速推論。GPU 約比 CPU 快 10-20 倍。若無 GPU 可關閉改用 CPU。" />
              <NumberInput label="Batch Size" value={params.batch_size}
                onChange={v => set('batch_size', v)} min={1} max={16}
                tooltip="每次同時送入 Cellpose 的影像 patch 數量。GPU 記憶體越大可設越高（建議 4-8）。記憶體不足時請降低。" />
            </Section>

            <Section title="Logic A 雙尺寸策略">
              <NumberInput label="小細胞直徑" value={params.dia_small}
                onChange={v => set('dia_small', v)} step={0.5} min={4} max={50} hint="px"
                tooltip="Cellpose 小尺寸推論的預期細胞直徑（像素）。設定接近真實細胞核直徑。偏小值 → 傾向過度分割；偏大值 → 傾向欠分割。" />
              <NumberInput label="大細胞直徑" value={params.dia_large}
                onChange={v => set('dia_large', v)} step={0.5} min={10} max={100} hint="px"
                tooltip="Cellpose 大尺寸推論的預期細胞直徑。邏輯 A 策略：若大尺寸 mask 內只包含 ≤1 個小尺寸 mask，則保留大尺寸（避免過度分割）；否則保留小尺寸細胞（避免欠分割）。" />
              <NumberInput label="Flow Threshold" value={params.flow_threshold}
                onChange={v => set('flow_threshold', v)} step={0.05} min={0} max={2}
                tooltip="Cellpose dP 光流向量誤差容忍閾值。值越大 → 接受品質更差的細胞邊界（召回率高但精確度低）。值越小 → 只接受高品質邊界。預設 0.4 適合多數場景。" />
              <NumberInput label="Cell Prob Threshold" value={params.cellprob_threshold}
                onChange={v => set('cellprob_threshold', v)} step={0.5} min={-6} max={6}
                tooltip="Cellpose 細胞存在機率的閾值。值越低（如 -3）→ 更容易偵測到細胞（高召回率）；值越高（如 2）→ 只有非常確定的細胞才被接受（高精確度）。預設 -1 適合細胞較密集的場景。" />
              <NumberInput label="Fragment Threshold" value={params.fragment_threshold}
                onChange={v => set('fragment_threshold', v)} min={0} max={500} hint="px²"
                tooltip="LOGIC A 合併時過濾碎片的面積閾值。面積小於此值的 mask 視為雜訊碎片並忽略。單位為像素平方（px²）。值越大 → 過濾掉更多小碎片；設 0 則不過濾。" />
            </Section>
          </div>

          {/* 右欄 */}
          <div className="space-y-5">
            <Section title="前處理">
              <Toggle label="Macenko 色彩標準化" value={params.normalize_stains}
                onChange={v => set('normalize_stains', v)}
                tooltip="Macenko 染色標準化：將 H&E 影像分解為 Hematoxylin（藍紫）與 Eosin（粉紅）兩個成分，並提取 Hematoxylin 通道作為 Cellpose 的灰階輸入。可消除不同批次染色差異。若影像已是灰階或染色品質差可關閉。" />
              <NumberInput label="CLAHE Clip Limit" value={params.clahe_clip_limit}
                onChange={v => set('clahe_clip_limit', v)} step={0.5} min={0.5} max={8}
                hint="細長/破碎核用1.0，一般細胞用2.0"
                tooltip="CLAHE（對比度限幅直方圖均等化）的裁剪限制值。值越高 → 對比度增強越強烈（細胞邊界更清晰，但雜訊也放大）。細長或稀疏的細胞核建議 1.0；一般圓形細胞核建議 2.0。" />
            </Section>

            <Section title="後處理">
              <Toggle label="Eosin Watershed 擴張" value={params.enable_eosin_watershed}
                onChange={v => set('enable_eosin_watershed', v)}
                tooltip="Nuclear Shield 演算法：利用 Eosin/亮度遮罩將 Cellpose 偵測到的細胞核輪廓，透過分水嶺演算法往外擴張至細胞質邊界。可生成包含細胞質的完整細胞 mask，供 Proseg 空間轉錄體定位使用。" />
              <NumberInput label="Eosin BG Threshold" value={params.eosin_bg_threshold}
                onChange={v => set('eosin_bg_threshold', v)} min={0} max={255}
                hint="灰階"
                tooltip="判斷背景（空腔）的亮度閾值。影像中 max(R,G,B) > (255 - 此值) 的像素視為空腔背景，Proseg 擴張被禁止進入此區域。值越大 → 判定更嚴苛（只排除極白的空腔）；值越小 → 排除更多偏亮區域。建議從 20-40 開始調試，搭配 Eosin 遮罩預覽確認。" />
            </Section>

            <Section title="分塊設定">
              <NumberInput label="Block Size" value={params.block_size}
                onChange={v => set('block_size', v)} step={256} min={512} max={8192} hint="px"
                tooltip="將大圖切分成多個小塊分別執行 Cellpose 的每塊邊長（像素）。值越大 → 需要更多 GPU 記憶體；值越小 → 切塊更多、邊界處理開銷更大。建議依 GPU 記憶體調整，常見值為 2048-4096。" />
              <NumberInput label="Overlap" value={params.overlap}
                onChange={v => set('overlap', v)} step={32} min={64} max={512} hint="px"
                tooltip="相鄰分塊的重疊區域大小（像素）。重疊可避免分塊邊界處細胞被截斷。建議設定約等於最大細胞直徑的 2-3 倍（一般 128-256 px 已足夠）。" />
            </Section>

            {/* 快速預設 */}
            <Section title="快速預設">
              <div className="flex gap-2 flex-wrap">
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS, model_type: 'cyto2', dia_small: 30, dia_large: 60, flow_threshold: 0.4, fragment_threshold: 200 })}
                  className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-200"
                >
                  CRC 上皮
                </button>
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS, model_type: 'nuclei', dia_small: 25, dia_large: 50, flow_threshold: 0.4, fragment_threshold: 200 })}
                  className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-200"
                >
                  LUAD 細胞核
                </button>
                <button
                  onClick={() => setParams({ ...DEFAULT_PARAMS, model_type: 'nuclei', dia_small: 15, dia_large: 40, flow_threshold: 0.3, cellprob_threshold: -3, fragment_threshold: 100 })}
                  className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-200"
                >
                  LUAD 正常肺
                </button>
              </div>
            </Section>
          </div>
        </div>
      </StageCard>

      {/* ── ROI 個別參數覆寫 ──────────────────────────────────────────────── */}
      {rois.length > 0 && (
        <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-3">
          <div className="flex items-center justify-between">
            <div>
              <h3 className="text-sm font-semibold text-gray-200">ROI 個別參數覆寫 &amp; 執行分割</h3>
              <p className="text-xs text-gray-500 mt-0.5">
                針對特定 ROI 覆寫參數後單獨重跑，或全部一起執行
              </p>
            </div>
            <div className="flex items-center gap-2">
              <button
                onClick={resetAllOverrides}
                className="text-xs text-gray-500 hover:text-gray-300 transition-colors"
              >
                全部重置
              </button>
              <button
                onClick={handleRunAll}
                disabled={stage.status === 'running'}
                className="px-4 py-1.5 text-sm rounded-lg font-medium transition-colors
                           bg-brand-primary text-white hover:bg-brand-primary/90
                           disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {stage.status === 'running' && runningRoi === null ? '執行中...' : '全部執行分割'}
              </button>
            </div>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="text-gray-500 border-b border-gray-700">
                  <th className="text-left py-2 pr-4 font-medium w-36">ROI</th>
                  <th className="text-center py-2 px-2 font-medium">模型</th>
                  <th className="text-center py-2 px-2 font-medium">小徑 (px)</th>
                  <th className="text-center py-2 px-2 font-medium">大徑 (px)</th>
                  <th className="text-center py-2 px-2 font-medium">Flow</th>
                  <th className="text-center py-2 px-2 font-medium">Cell Prob</th>
                  <th className="text-center py-2 px-2 font-medium">Fragment</th>
                  <th className="text-center py-2 px-2 font-medium">重新分割</th>
                  <th className="py-2 px-2 w-8"></th>
                </tr>
              </thead>
              <tbody>
                {rois.map(roi => {
                  const ov: RoiOverride = roiOverrides[roi.name] ?? {}
                  const hasAny = Object.values(ov).some(v => v != null)

                  // 小工具：顯示「全域」按鈕或數字輸入
                  const NumCell = ({
                    field, defaultVal, step = 0.5, min, max,
                  }: {
                    field: keyof RoiOverride; defaultVal: number; step?: number; min?: number; max?: number
                  }) => {
                    const val = ov[field]
                    return val != null ? (
                      <div className="flex items-center justify-center gap-0.5">
                        <input
                          type="number"
                          value={val as number}
                          step={step}
                          min={min}
                          max={max}
                          onChange={e => updateRoiField(roi.name, field, parseFloat(e.target.value))}
                          className="w-16 px-1 py-0.5 bg-gray-800 border border-blue-500 rounded text-gray-100 text-xs text-right focus:outline-none"
                        />
                        <button
                          onClick={() => clearRoiField(roi.name, field)}
                          className="text-gray-500 hover:text-red-400 leading-none px-0.5"
                          title="還原為全域"
                        >×</button>
                      </div>
                    ) : (
                      <button
                        onClick={() => updateRoiField(roi.name, field, defaultVal)}
                        className="mx-auto block px-2 py-0.5 rounded bg-gray-700/60 text-gray-500 hover:bg-gray-600 hover:text-gray-200 transition-colors"
                      >
                        全域
                      </button>
                    )
                  }

                  return (
                    <tr
                      key={roi.name}
                      className={`border-b border-gray-800/60 ${hasAny ? 'border-l-2 border-l-blue-500' : 'border-l-2 border-l-transparent'}`}
                    >
                      <td className="py-2 pr-4 pl-1">
                        <span className="text-gray-200 font-medium">{roi.name}</span>
                        {roi.tissue && <span className="text-gray-600 ml-1.5 text-xs">{roi.tissue}</span>}
                      </td>

                      {/* 模型類型 */}
                      <td className="py-2 px-2">
                        {ov.model_type != null ? (
                          <div className="flex items-center justify-center gap-0.5">
                            <select
                              value={ov.model_type}
                              onChange={e => updateRoiField(roi.name, 'model_type', e.target.value)}
                              className="px-1 py-0.5 bg-gray-800 border border-blue-500 rounded text-gray-100 text-xs focus:outline-none"
                            >
                              <option value="cyto2">cyto2</option>
                              <option value="cyto3">cyto3</option>
                              <option value="nuclei">nuclei</option>
                            </select>
                            <button
                              onClick={() => clearRoiField(roi.name, 'model_type')}
                              className="text-gray-500 hover:text-red-400 leading-none px-0.5"
                              title="還原為全域"
                            >×</button>
                          </div>
                        ) : (
                          <button
                            onClick={() => updateRoiField(roi.name, 'model_type', params.model_type)}
                            className="mx-auto block px-2 py-0.5 rounded bg-gray-700/60 text-gray-500 hover:bg-gray-600 hover:text-gray-200 transition-colors"
                          >
                            全域
                          </button>
                        )}
                      </td>

                      <td className="py-2 px-2">
                        <NumCell field="dia_small" defaultVal={params.dia_small} step={0.5} min={4} max={50} />
                      </td>
                      <td className="py-2 px-2">
                        <NumCell field="dia_large" defaultVal={params.dia_large} step={0.5} min={10} max={100} />
                      </td>
                      <td className="py-2 px-2">
                        <NumCell field="flow_threshold" defaultVal={params.flow_threshold} step={0.05} min={0} max={2} />
                      </td>
                      <td className="py-2 px-2">
                        <NumCell field="cellprob_threshold" defaultVal={params.cellprob_threshold} step={0.5} min={-6} max={6} />
                      </td>
                      <td className="py-2 px-2">
                        <NumCell field="fragment_threshold" defaultVal={params.fragment_threshold} step={10} min={0} max={500} />
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
                          {stage.status === 'running' && runningRoi === roi.name ? '執行中...' : '執行分割'}
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
            <h3 className="text-sm font-semibold text-gray-200">快速 Patch 預覽</h3>
            <p className="text-xs text-gray-500 mt-0.5">
              從 he_crop.tif 取一小塊國 Cellpose，不需先執行完整分割 · 使用上方當前參數
            </p>
          </div>
          <div className="flex gap-2">
            <button
              onClick={handlePreprocPreview}
              disabled={preprocLoading}
              className="px-4 py-1.5 text-sm rounded bg-gray-600 text-gray-200 font-medium
                         hover:bg-gray-500 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              {preprocLoading ? '處理中...' : '⚡ 前處理預覽'}
            </button>
            <button
              onClick={handleQuickPreview}
              disabled={prevLoading}
              className="px-4 py-1.5 text-sm rounded bg-primary text-black font-medium
                         hover:bg-primary/90 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              {prevLoading ? '執行中...' : '🔬 分割預覽'}
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
            <label className="text-xs text-gray-400 mb-1 block">X 起點（px）</label>
            <input
              type="number" value={prevX} min={0}
              onChange={e => setPrevX(Math.max(0, parseInt(e.target.value) || 0))}
              className="w-full px-2 py-1.5 text-sm bg-gray-800 border border-gray-600 rounded
                         text-gray-100 text-right focus:outline-none focus:border-blue-500"
            />
          </div>

          <div>
            <label className="text-xs text-gray-400 mb-1 block">Y 起點（px）</label>
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
                    <span className="text-green-400 font-medium">✓ {quickInfo.n_cells} 個細胞</span>
                  )}
                  <span>ROI: {quickInfo.roi_name}</span>
                  <span>{quickInfo.patch_info}</span>
                </div>
                {/* 圖層切換 Tab */}
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  {(['overlay', 'macenko', 'flows', 'cyto'] as const).map(tab => {
                    const labels: Record<typeof tab, string> = {
                      overlay: 'H&E + 邊界',
                      macenko: 'Macenko 前處理',
                      flows: 'Flow 方向圖',
                      cyto: 'Eosin 細胞質遮罩'
                    }
                    const available = tab === 'overlay' || (tab === 'macenko' && !!quickMacenko) || (tab === 'flows' && !!quickFlows) || (tab === 'cyto' && !!quickCyto)
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
              {quickTab === 'overlay' && 'H&E 原圖 + 綠色細胞邊界（LOGIC_A 合併）'}
              {quickTab === 'macenko' && 'Macenko Hematoxylin 萃取 → CLAHE 增強（Cellpose 實際輸入）'}
              {quickTab === 'flows' && 'Cellpose 小尺寸 dP 光流方向圖（色相 = 方向，飽和度 = 強度）；白線 = 細胞邊界'}
              {quickTab === 'cyto' && 'Eosin 背景預先過濾遮罩 (紅底 = 細胞質)'}
            </p>
            <div className="rounded-lg overflow-hidden border border-surface-border"
              style={{ imageRendering: 'pixelated' }}>
              <img
                src={
                  quickTab === 'macenko' ? (quickMacenko ?? quickSrc) :
                    quickTab === 'flows' ? (quickFlows ?? quickSrc) :
                      quickTab === 'cyto' ? (quickCyto ?? quickSrc) :
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
              <h3 className="text-sm font-semibold text-gray-200">完整分割結果預覽</h3>
              <p className="text-xs text-gray-500 mt-0.5">H&E + 綠色細胞邊界疊圖（來自已存遮罩）</p>
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
                {previewSrc ? '重新載入' : '載入預覽'}
              </button>
            </div>
          </div>
          {previewSrc && (
            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <span className="text-xs text-green-400 font-medium">
                  {previewNCells != null && `✓ ${previewNCells.toLocaleString()} 個細胞 · ROI: ${previewRoi}`}
                </span>
                {/* Tab 切換 */}
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  <button
                    onClick={() => setPreviewTab('overlay')}
                    className={`px-3 py-1 transition-colors ${previewTab === 'overlay' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
                  >H&E + 邊界</button>
                  <button
                    disabled={!previewFlows}
                    onClick={() => setPreviewTab('flows')}
                    className={`px-3 py-1 transition-colors disabled:opacity-30 disabled:cursor-not-allowed ${previewTab === 'flows' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
                  >光流方向 (Flows)</button>
                  <button
                    disabled={!previewCyto}
                    onClick={() => setPreviewTab('cyto')}
                    className={`px-3 py-1 transition-colors disabled:opacity-30 disabled:cursor-not-allowed ${previewTab === 'cyto' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
                  >Eosin 細胞質遮罩</button>
                </div>
              </div>
              <p className="text-xs text-gray-500">
                {previewTab === 'overlay' && 'H&E 原圖 + 綠色細胞邊界（來自已存遮罩）'}
                {previewTab === 'flows' && 'Cellpose 小尺寸 dP 光流方向圖（色相 = 方向，飽和度 = 強度）'}
                {previewTab === 'cyto' && 'Eosin 背景過濾遮罩（亮色 = 組織，暗色 = 空腔）'}
              </p>
              <div className="rounded-lg overflow-hidden border border-surface-border">
                <img src={previewTab === 'cyto' ? (previewCyto ?? previewSrc!) : previewTab === 'flows' ? (previewFlows ?? previewSrc!) : previewSrc!} alt="segmentation preview" className="w-full" />
              </div>
            </div>
          )}
        </div>
      )}

      <Terminal stage="segmentation" />
    </div>
  )
}
