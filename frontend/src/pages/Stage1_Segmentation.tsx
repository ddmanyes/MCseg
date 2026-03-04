import { useState } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runSegmentation, getSegmentationStatus, getSegmentationPreview, runSegmentationPreview } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

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
  enable_eosin_watershed: true,
  eosin_bg_threshold: 80,
  block_size: 2048,
  overlap: 256,
}

function NumberInput({
  label, value, onChange, step = 1, min, max, hint,
}: {
  label: string; value: number; onChange: (v: number) => void
  step?: number; min?: number; max?: number; hint?: string
}) {
  return (
    <div className="flex items-center justify-between gap-3">
      <div className="flex-1">
        <span className="text-sm text-gray-300">{label}</span>
        {hint && <span className="text-xs text-gray-500 ml-1">({hint})</span>}
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

function Toggle({ label, value, onChange, hint }: {
  label: string; value: boolean; onChange: (v: boolean) => void; hint?: string
}) {
  return (
    <div className="flex items-center justify-between">
      <div>
        <span className="text-sm text-gray-300">{label}</span>
        {hint && <span className="text-xs text-gray-500 ml-1">({hint})</span>}
      </div>
      <button
        onClick={() => onChange(!value)}
        className={`relative w-10 h-5 rounded-full transition-colors ${
          value ? 'bg-blue-600' : 'bg-gray-600'
        }`}
      >
        <span className={`absolute top-0.5 w-4 h-4 bg-white rounded-full shadow transition-transform ${
          value ? 'translate-x-5' : 'translate-x-0.5'
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
  const [previewSrc, setPreviewSrc]         = useState<string | null>(null)
  const [previewRoi, setPreviewRoi]         = useState('')
  const [previewAvail, setPreviewAvail]     = useState<string[]>([])
  const [previewNCells, setPreviewNCells]   = useState<number | null>(null)

  // ── Quick preview state ─────────────────────────────────────────────────
  const [prevRoi, setPrevRoi]            = useState('')
  const [prevX, setPrevX]                = useState(0)
  const [prevY, setPrevY]                = useState(0)
  const [prevPatchSize, setPrevPatchSize] = useState(512)
  const [prevLoading, setPrevLoading]    = useState(false)
  const [quickSrc, setQuickSrc]          = useState<string | null>(null)
  const [quickMacenko, setQuickMacenko]  = useState<string | null>(null)
  const [quickFlows, setQuickFlows]      = useState<string | null>(null)
  const [quickTab, setQuickTab]          = useState<'overlay' | 'macenko' | 'flows'>('overlay')
  const [quickInfo, setQuickInfo]        = useState<{ n_cells: number; roi_name: string; patch_info: string } | null>(null)
  const [quickError, setQuickError]      = useState<string | null>(null)

  const set = <K extends keyof SegParams>(key: K, value: SegParams[K]) =>
    setParams(prev => ({ ...prev, [key]: value }))

  const handleRun = async () => {
    updateStage('segmentation', { status: 'running', progress: 0, message: '啟動 Cellpose...' })
    await runSegmentation(params)
    refetchStatus()
  }

  const handlePreview = async (roi?: string) => {
    const res = await getSegmentationPreview(roi)
    const d = res.data?.data
    if (d?.image_b64) {
      setPreviewSrc(`data:image/jpeg;base64,${d.image_b64}`)
      if (d.available_rois) setPreviewAvail(d.available_rois)
      if (d.roi)            setPreviewRoi(d.roi)
      if (d.n_cells != null) setPreviewNCells(d.n_cells)
    }
  }

  const handleQuickPreview = async () => {
    setPrevLoading(true)
    setQuickError(null)
    setQuickSrc(null)
    setQuickMacenko(null)
    setQuickFlows(null)
    setQuickInfo(null)
    try {
      const res = await runSegmentationPreview({
        roi_name:           prevRoi || undefined,
        x:                  prevX,
        y:                  prevY,
        patch_size:         prevPatchSize,
        model_type:         params.model_type,
        use_gpu:            params.use_gpu,
        dia_small:          params.dia_small,
        dia_large:          params.dia_large,
        flow_threshold:     params.flow_threshold,
        cellprob_threshold: params.cellprob_threshold,
        fragment_threshold: params.fragment_threshold,
        normalize_stains:   params.normalize_stains,
      })
      const d = res.data
      if (d?.status === 'ok' && d.data?.image_b64) {
        setQuickSrc(`data:image/jpeg;base64,${d.data.image_b64}`)
        if (d.data.macenko_b64) setQuickMacenko(`data:image/jpeg;base64,${d.data.macenko_b64}`)
        if (d.data.flows_b64)   setQuickFlows(`data:image/jpeg;base64,${d.data.flows_b64}`)
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

  return (
    <div className="space-y-4">
      <StageCard
        title="細胞分割（Cellpose + Logic A）"
        status={stage.status}
        progress={stage.progress}
        message={stage.message}
        onRun={handleRun}
        runLabel="執行分割"
      >
        {/* 參數面板 */}
        <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-6">

          {/* 左欄 */}
          <div className="space-y-5">
            <Section title="模型設定">
              <div className="flex items-center justify-between">
                <div>
                  <span className="text-sm text-gray-300">執行範圍</span>
                  <span className="text-xs text-gray-500 ml-1">
                    {params.mode === 'roi' ? '（各 ROI 裁切圖）' : '（完整 BTF，耗時）'}
                  </span>
                </div>
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  <button
                    onClick={() => set('mode', 'roi')}
                    className={`px-3 py-1 transition-colors ${
                      params.mode === 'roi' ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'
                    }`}
                  >ROI</button>
                  <button
                    onClick={() => set('mode', 'full')}
                    className={`px-3 py-1 transition-colors ${
                      params.mode === 'full' ? 'bg-orange-600 text-white' : 'bg-gray-800 text-gray-400 hover:bg-gray-700'
                    }`}
                  >全圖</button>
                </div>
              </div>
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-300">模型類型</span>
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
              <Toggle label="GPU 加速" value={params.use_gpu} onChange={v => set('use_gpu', v)} />
              <NumberInput label="Batch Size" value={params.batch_size}
                onChange={v => set('batch_size', v)} min={1} max={16} />
            </Section>

            <Section title="Logic A 雙尺寸策略">
              <NumberInput label="小細胞直徑" value={params.dia_small}
                onChange={v => set('dia_small', v)} step={0.5} min={4} max={50} hint="px" />
              <NumberInput label="大細胞直徑" value={params.dia_large}
                onChange={v => set('dia_large', v)} step={0.5} min={10} max={100} hint="px" />
              <NumberInput label="Flow Threshold" value={params.flow_threshold}
                onChange={v => set('flow_threshold', v)} step={0.05} min={0} max={2} />
              <NumberInput label="Cell Prob Threshold" value={params.cellprob_threshold}
                onChange={v => set('cellprob_threshold', v)} step={0.5} min={-6} max={6} />
              <NumberInput label="Fragment Threshold" value={params.fragment_threshold}
                onChange={v => set('fragment_threshold', v)} min={0} max={500} hint="px²" />
            </Section>
          </div>

          {/* 右欄 */}
          <div className="space-y-5">
            <Section title="前處理">
              <Toggle label="Macenko 色彩標準化" value={params.normalize_stains}
                onChange={v => set('normalize_stains', v)} />
            </Section>

            <Section title="後處理">
              <Toggle label="Eosin Watershed 擴張" value={params.enable_eosin_watershed}
                onChange={v => set('enable_eosin_watershed', v)} />
              <NumberInput label="Eosin BG Threshold" value={params.eosin_bg_threshold}
                onChange={v => set('eosin_bg_threshold', v)} min={0} max={255}
                hint="灰階" />
            </Section>

            <Section title="分塊設定">
              <NumberInput label="Block Size" value={params.block_size}
                onChange={v => set('block_size', v)} step={256} min={512} max={8192} hint="px" />
              <NumberInput label="Overlap" value={params.overlap}
                onChange={v => set('overlap', v)} step={32} min={64} max={512} hint="px" />
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

      {/* ── 快速 Patch 預覽 ─────────────────────────────────────────────────── */}
      <div className="rounded-xl bg-surface border border-surface-border p-4 space-y-4">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-sm font-semibold text-gray-200">快速 Patch 預覽</h3>
            <p className="text-xs text-gray-500 mt-0.5">
              從 he_crop.tif 取一小塊跑 Cellpose，不需先執行完整分割 · 使用上方當前參數
            </p>
          </div>
          <button
            onClick={handleQuickPreview}
            disabled={prevLoading}
            className="px-4 py-1.5 text-sm rounded bg-primary text-black font-medium
                       hover:bg-primary/90 disabled:opacity-50 disabled:cursor-not-allowed
                       transition-colors"
          >
            {prevLoading ? '執行中...' : '執行預覽'}
          </button>
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
                  <span className="text-green-400 font-medium">✓ {quickInfo.n_cells} 個細胞</span>
                  <span>ROI: {quickInfo.roi_name}</span>
                  <span>{quickInfo.patch_info}</span>
                </div>
                {/* 圖層切換 Tab */}
                <div className="flex rounded overflow-hidden border border-gray-600 text-xs">
                  {(['overlay', 'macenko', 'flows'] as const).map(tab => {
                    const labels: Record<typeof tab, string> = {
                      overlay: 'H&E + 邊界',
                      macenko: 'Macenko 前處理',
                      flows:   'Flow 方向圖',
                    }
                    const available = tab === 'overlay' || (tab === 'macenko' && !!quickMacenko) || (tab === 'flows' && !!quickFlows)
                    return (
                      <button
                        key={tab}
                        disabled={!available}
                        onClick={() => setQuickTab(tab)}
                        className={`px-3 py-1 transition-colors disabled:opacity-30 disabled:cursor-not-allowed ${
                          quickTab === tab
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
              {quickTab === 'overlay'  && 'H&E 原圖 + 綠色細胞邊界（LOGIC_A 合併）'}
              {quickTab === 'macenko'  && 'Macenko Hematoxylin 萃取 → CLAHE 增強（Cellpose 實際輸入）'}
              {quickTab === 'flows'    && 'Cellpose 小尺寸 dP 光流方向圖（色相 = 方向，飽和度 = 強度）；白線 = 細胞邊界'}
            </p>
            <div className="rounded-lg overflow-hidden border border-surface-border"
                 style={{ imageRendering: 'pixelated' }}>
              <img
                src={
                  quickTab === 'macenko' ? (quickMacenko ?? quickSrc) :
                  quickTab === 'flows'   ? (quickFlows   ?? quickSrc) :
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
              {previewNCells != null && (
                <span className="text-xs text-green-400 font-medium">
                  ✓ {previewNCells.toLocaleString()} 個細胞 · ROI: {previewRoi}
                </span>
              )}
              <div className="rounded-lg overflow-hidden border border-surface-border">
                <img src={previewSrc} alt="segmentation preview" className="w-full" />
              </div>
            </div>
          )}
        </div>
      )}

      <Terminal stage="segmentation" />
    </div>
  )
}
