import { useState, useEffect, useCallback, useRef } from 'react'
import { useQuery, useQueryClient } from '@tanstack/react-query'
import { usePipelineStore } from '../stores/pipelineStore'
import Terminal from '../components/shared/Terminal'
import useStageLog from '../hooks/useStageLog'
import {
  runQC, getQCStatus, getQCImages,
  runUMAPExplore, getUMAPExploreStatus, getUMAPImages,
  runHeatmap, getHeatmapStatus, getHeatmapImage,
  getConfig, getOverlayHdUrl, getAvailableRois, getRoiOverlays,
  getClusterInfo, getCelltypistModels,
  runAnnotate, getAnnotateStatus, applyLabels,
  getRawHistogram,
} from '../api/client'
import QcHistogram, { type HistogramMetric } from '../components/shared/QcHistogram'

// ── 小工具 ────────────────────────────────────────────────────────

function SectionHeader({ step, title, subtitle }: { step: number; title: string; subtitle?: string }) {
  return (
    <div className="flex items-center gap-3 mb-4">
      <div className="w-7 h-7 rounded-full bg-brand-primary flex items-center justify-center text-sm font-bold text-white shrink-0">
        {step}
      </div>
      <div>
        <h3 className="font-semibold text-gray-100 text-base">{title}</h3>
        {subtitle && <p className="text-xs text-gray-400 mt-0.5">{subtitle}</p>}
      </div>
    </div>
  )
}

function StatusBadge({ status, message }: { status: string; message: string }) {
  const colors: Record<string, string> = {
    idle: 'text-gray-400',
    running: 'text-yellow-400',
    done: 'text-green-400',
    error: 'text-red-400',
  }
  const icons: Record<string, string> = {
    idle: '○', running: '⟳', done: '✓', error: '✗',
  }
  return (
    <span className={`text-xs font-medium ${colors[status] ?? 'text-gray-400'}`}>
      {icons[status] ?? '○'} {message}
    </span>
  )
}

function RunButton({
  label, onClick, status, disabled,
}: { label: string; onClick: () => void; status: string; disabled?: boolean }) {
  const isRunning = status === 'running'
  return (
    <button
      onClick={onClick}
      disabled={isRunning || disabled}
      className={`px-4 py-1.5 rounded-lg text-sm font-medium transition-colors ${
        isRunning || disabled
          ? 'bg-gray-700 text-gray-500 cursor-not-allowed'
          : 'bg-primary text-white hover:bg-primary-dark'
      }`}
    >
      {isRunning ? '執行中...' : label}
    </button>
  )
}

function NumberField({
  label, value, onChange, step = 1, min, hint,
}: { label: string; value: number; onChange: (v: number) => void; step?: number; min?: number; hint?: string }) {
  return (
    <div>
      <label className="block text-xs text-gray-400 mb-1">{label}</label>
      <input
        type="number"
        step={step}
        min={min}
        className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
        value={value}
        onChange={e => onChange(parseFloat(e.target.value) || 0)}
      />
      {hint && <p className="text-xs text-gray-500 mt-0.5">{hint}</p>}
    </div>
  )
}

// fullWidthKeys：指定哪些 tab key 要全寬顯示，其餘縮小為 50%
// 未傳入時預設全部全寬
function ChartView({
  images, tabs, fullWidthKeys,
}: {
  images: Record<string, string>
  tabs: { key: string; label: string }[]
  fullWidthKeys?: string[]
}) {
  const available = tabs.filter(t => images[t.key])
  const [active, setActive] = useState(available[0]?.key ?? '')
  useEffect(() => {
    if (available.length && !images[active]) setActive(available[0].key)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [images, active])
  if (!available.length) return null

  const isFullWidth = fullWidthKeys ? fullWidthKeys.includes(active) : true

  return (
    <div className="mt-4">
      <div className="flex gap-1 border-b border-surface-border mb-3">
        {available.map(t => (
          <button
            key={t.key}
            onClick={() => setActive(t.key)}
            className={`px-3 py-1.5 text-xs font-medium rounded-t transition-colors ${
              active === t.key
                ? 'bg-surface-card text-brand-primary border-b-2 border-brand-primary'
                : 'text-gray-400 hover:text-gray-200'
            }`}
          >
            {t.label}
          </button>
        ))}
      </div>
      {images[active] && (
        <img
          src={`data:image/png;base64,${images[active]}`}
          className={`rounded-lg ${isFullWidth ? 'max-w-full' : 'max-w-[50%]'}`}
          alt={active}
        />
      )}
    </div>
  )
}

// ── ROI 輪廓比較面板 ─────────────────────────────────────────────

type RoiOverlayData = Record<string, { pre_qc?: string; post_qc?: string }>

function RoiContourPanel({ data }: { data: RoiOverlayData }) {
  const roiNames = Object.keys(data)
  const [activeRoi, setActiveRoi] = useState(roiNames[0] ?? '')
  const [open, setOpen] = useState(false)

  if (!roiNames.length) return null

  const current = data[activeRoi] ?? {}
  const stages: { key: 'pre_qc' | 'post_qc'; label: string }[] = [
    { key: 'pre_qc',  label: 'Pre-QC（全部細胞）' },
    { key: 'post_qc', label: 'Post-QC（保留 / 移除）' },
  ]

  return (
    <div className="mt-4 border border-surface-border rounded-lg overflow-hidden">
      {/* 標題列（可收折）*/}
      <button
        onClick={() => setOpen(o => !o)}
        className="w-full flex items-center justify-between px-4 py-2.5 bg-surface-darker hover:bg-surface-highlight transition-colors text-left"
      >
        <span className="text-xs font-medium text-gray-300">
          細胞輪廓比較（QC 前後）
        </span>
        <span className="text-gray-500 text-xs">{open ? '▲ 收起' : '▼ 展開'}</span>
      </button>

      {open && (
        <div className="bg-surface-card p-4">
          {/* ROI tabs（多 ROI 時才顯示）*/}
          {roiNames.length > 1 && (
            <div className="flex gap-1 border-b border-surface-border mb-3">
              {roiNames.map(name => (
                <button
                  key={name}
                  onClick={() => setActiveRoi(name)}
                  className={`px-3 py-1.5 text-xs font-medium rounded-t transition-colors ${
                    activeRoi === name
                      ? 'bg-surface-darker text-brand-primary border-b-2 border-brand-primary'
                      : 'text-gray-400 hover:text-gray-200'
                  }`}
                >
                  {name}
                </button>
              ))}
            </div>
          )}

          {/* 2 欄並排：Pre / Post */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {stages.map(({ key, label }) =>
              current[key] ? (
                <div key={key}>
                  <p className="text-xs text-gray-400 mb-1">{label}</p>
                  <img
                    src={`data:image/png;base64,${current[key]}`}
                    className="rounded w-full"
                    alt={`${activeRoi} ${key}`}
                  />
                </div>
              ) : null
            )}
          </div>

          {/* 跨 ROI 摘要列 */}
          {roiNames.length > 1 && (
            <div className="mt-3 text-xs text-gray-500">
              共 {roiNames.length} 個 ROI：{roiNames.join('、')}。點擊上方 Tab 切換。
            </div>
          )}
        </div>
      )}
    </div>
  )
}

// ── 主元件 ────────────────────────────────────────────────────────

export default function Stage4_Analysis() {
  useStageLog('analysis')
  const { updateStage } = usePipelineStore()
  const queryClient = useQueryClient()

  // ── 分析來源選擇 ──
  const [analysisMode, setAnalysisMode] = useState<'single' | 'merge'>('single')
  const [selectedRoi, setSelectedRoi] = useState<string>('')
  const [inputSource, setInputSource] = useState<'cellpose' | 'proseg'>('cellpose')

  const { data: availableRoisData } = useQuery({
    queryKey: ['available_rois'],
    queryFn: async () => (await getAvailableRois()).data,
  })
  const availableRois: { name: string; available: boolean; has_cellpose?: boolean; has_proseg?: boolean }[] =
    availableRoisData?.data ?? []
  const hasMultipleRois = availableRois.length > 1
  const hasProsegRois   = availableRois.some(r => r.has_proseg)

  // ── QC 參數 ──
  const [qcParams, setQcParams] = useState({
    min_genes: 10, max_genes: 8000,
    min_counts: 5, min_cells: 3,
    max_pct_mito: 80.0, n_top_genes: 2000, n_pcs: 30,
    min_complexity: 0.8,
  })

  // ── UMAP 參數 ──
  const [umapParams, setUmapParams] = useState({
    n_pcs: 30, n_neighbors: 15, min_dist: 0.3,
  })
  const [resolutionInput, setResolutionInput] = useState('0.3, 0.5, 0.8')

  // ── 標註參數 ──
  const [annotateRes, setAnnotateRes] = useState<string>('')
  const [annotateModel, setAnnotateModel] = useState('Human_Colorectal_Cancer.pkl')
  const [annotateMode, setAnnotateMode] = useState<'dual' | 'single'>('dual')
  const [immuneConfThreshold, setImmuneConfThreshold] = useState(0.5)
  const [celltypistModels, setCelltypistModels] = useState<Record<string, string>>({})
  const [clusterLabels, setClusterLabels] = useState<Record<string, string>>({})   // cluster_id → label
  const [clusterMeta, setClusterMeta] = useState<Record<string, {
    confidence: number
    source: string
    uncertain?: boolean
    state?: string | null
    state_score?: number
    tier3_label?: string
    tier3_conf?: number
    immune_label?: string
    immune_conf?: number
    crc_label?: string
    crc_conf?: number
  }>>({})
  const [scoreThreshold, setScoreThreshold] = useState(0.3)
  const [uncertainThreshold, setUncertainThreshold] = useState(0.7)
  const [enableTier3, setEnableTier3] = useState(false)
  const [tier3ConfThreshold, setTier3ConfThreshold] = useState(0.6)
  const [labelApplied, setLabelApplied] = useState(false)

  // ── Heatmap 參數 ──
  const [selectedRes, setSelectedRes] = useState<string>('')
  const [nTopGenes, setNTopGenes] = useState(20)
  const [nHeatmapGenes, setNHeatmapGenes] = useState(50)

  // ── 原始分布直方圖 ──
  const [histData, setHistData] = useState<{
    n_cells: number
    metrics: Record<string, HistogramMetric>
  } | null>(null)
  const [histLoading, setHistLoading] = useState(false)
  const [histError, setHistError] = useState('')
  const [applyLabelError, setApplyLabelError] = useState('')
  const [logScales, setLogScales] = useState<Record<string, boolean>>({})

  // ── 儲存圖表 ──
  const [qcImages, setQcImages] = useState<Record<string, string>>({})
  const [umapImages, setUmapImages] = useState<Record<string, string>>({})
  const [heatmapImages, setHeatmapImages] = useState<Record<string, string>>({})
  const [roiOverlays, setRoiOverlays] = useState<RoiOverlayData>({})

  // 初始化：自動選第一個有資料的 ROI
  useEffect(() => {
    if (!selectedRoi && availableRois.length > 0) {
      const first = availableRois.find(r => r.available)
      if (first) setSelectedRoi(first.name)
    }
  }, [availableRois])

  // ── 從設定載入初始值 ──
  useEffect(() => {
    getConfig().then(res => {
      const ana = res.data?.data?.analysis
      if (!ana) return
      const c = ana.preprocessing?.cellular ?? {}
      const h = ana.preprocessing?.hvg ?? {}
      const cl = ana.clustering ?? {}
      setQcParams(p => ({
        min_genes: c.min_genes ?? p.min_genes,
        max_genes: c.max_genes ?? p.max_genes,
        min_counts: c.min_counts ?? p.min_counts,
        min_cells: c.min_cells ?? p.min_cells,
        max_pct_mito: c.max_pct_mito ?? p.max_pct_mito,
        min_complexity: c.min_complexity ?? p.min_complexity,
        n_top_genes: h.n_top_genes ?? p.n_top_genes,
        n_pcs: cl.n_pcs ?? p.n_pcs,
      }))
      setUmapParams(p => ({
        n_pcs: cl.n_pcs ?? p.n_pcs,
        n_neighbors: cl.n_neighbors ?? p.n_neighbors,
        min_dist: cl.min_dist ?? p.min_dist,
      }))
    })
  }, [])

  // ── 載入 CellTypist 模型清單 ──
  useEffect(() => {
    getCelltypistModels().then(r => {
      if (r.data?.data) setCelltypistModels(r.data.data)
    })
  }, [])

  // ── TanStack Query: 四步驟 status 輪詢 ──
  const { data: qcSt } = useQuery({
    queryKey: ['qc_status'],
    queryFn: async () => (await getQCStatus()).data,
    refetchInterval: (q) => q.state.data?.status === 'running' ? 2000 : false,
  })

  const { data: umapSt } = useQuery({
    queryKey: ['umap_explore_status'],
    queryFn: async () => (await getUMAPExploreStatus()).data,
    refetchInterval: (q) => q.state.data?.status === 'running' ? 2000 : false,
  })

  const { data: heatSt } = useQuery({
    queryKey: ['heatmap_status'],
    queryFn: async () => (await getHeatmapStatus()).data,
    refetchInterval: (q) => q.state.data?.status === 'running' ? 2000 : false,
  })

  const { data: annotSt } = useQuery({
    queryKey: ['annotate_status'],
    queryFn: async () => (await getAnnotateStatus()).data,
    refetchInterval: (q) => q.state.data?.status === 'running' ? 2000 : false,
  })

  // ── 各步驟完成後自動拉取圖表（包含頁面初次載入 / 後端重啟後恢復） ──
  useEffect(() => {
    if (qcSt?.status === 'done') {
      getQCImages().then(r => { if (r.data.data) setQcImages(r.data.data) })
      getRoiOverlays().then(r => { if (r.data?.data) setRoiOverlays(r.data.data) })
    }
  }, [qcSt?.status])

  useEffect(() => {
    if (umapSt?.status === 'done') {
      getUMAPImages().then(r => {
        if (r.data.data) {
          setUmapImages(r.data.data)
          const keys = Object.keys(r.data.data).filter(k => k !== 'grid')
          if (keys.length) setSelectedRes(keys[0])
        }
      })
    }
  }, [umapSt?.status])

  useEffect(() => {
    if (heatSt?.status === 'done') {
      getHeatmapImage().then(r => { if (r.data.data) setHeatmapImages(r.data.data) })
      updateStage('analysis', { status: 'done', progress: 1, message: '分析完成' })
    }
  }, [heatSt?.status])

  // CellTypist 完成後自動填入建議標籤 + 信心分數
  useEffect(() => {
    if (annotSt?.status === 'done' && annotSt?.suggestions) {
      const suggestions = annotSt.suggestions as Record<string, {
        label: string; confidence: number; source: string
        immune_label?: string; immune_conf?: number
        crc_label?: string; crc_conf?: number
      } | string>  // 相容舊格式
      if (Object.keys(suggestions).length === 0) return
      const labels: Record<string, string> = {}
      const meta: Record<string, { confidence: number; source: string; immune_label?: string; immune_conf?: number; crc_label?: string; crc_conf?: number }> = {}
      Object.entries(suggestions).forEach(([cluster, info]) => {
        if (typeof info === 'string') {
          labels[cluster] = info
          meta[cluster] = { confidence: 0, source: 'single' }
        } else {
          labels[cluster] = info.label
          meta[cluster] = {
            confidence: info.confidence,
            source: info.source,
            immune_label: info.immune_label,
            immune_conf: info.immune_conf,
            crc_label: info.crc_label,
            crc_conf: info.crc_conf,
          }
        }
      })
      setClusterLabels(prev => ({ ...prev, ...labels }))
      setClusterMeta(meta)
    }
  }, [annotSt?.status])

  // ── 可用的 resolution 列表（供 Heatmap / 標註下拉）——必須在用到它的 useEffect 之前宣告 ──
  const availableResolutions = Object.keys(umapImages).filter(k => k !== 'grid').sort()

  // UMAP 圖表載入後（availableResolutions 更新），自動設定 annotateRes 並載入 cluster 資訊
  useEffect(() => {
    if (availableResolutions.length === 0 || annotateRes) return
    const res = availableResolutions[0]
    setAnnotateRes(res)
    getClusterInfo(parseFloat(res)).then(r => {
      if (r.data?.data) {
        const { cluster_ids, existing_labels } = r.data.data
        const init: Record<string, string> = {}
        cluster_ids.forEach((id: string) => { init[id] = existing_labels[id] ?? '' })
        setClusterLabels(init)
        if (Object.values(existing_labels).some(v => v)) setLabelApplied(true)
      }
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [availableResolutions.join(',')])

  // ── 解析 resolution 文字輸入 ──
  const parseResolutions = useCallback((): number[] => {
    return resolutionInput
      .split(/[,\s]+/)
      .map(s => parseFloat(s.trim()))
      .filter(v => !isNaN(v) && v > 0)
  }, [resolutionInput])

  // ── 直方圖 handlers ──
  const handleLoadHist = async (src?: 'cellpose' | 'proseg') => {
    setHistLoading(true)
    setHistError('')
    try {
      const mergeFlag = analysisMode === 'merge' && hasMultipleRois
      const res = await getRawHistogram(
        mergeFlag ? undefined : (selectedRoi || undefined),
        mergeFlag,
        src ?? inputSource,
      )
      if (res.data?.status === 'ok') {
        setHistData(res.data.data)
      } else {
        setHistError(res.data?.message ?? '載入失敗')
      }
    } catch (e: unknown) {
      setHistError(e instanceof Error ? e.message : '載入失敗')
    } finally {
      setHistLoading(false)
    }
  }

  const handleApplyMad = () => {
    if (!histData) return
    const m = histData.metrics
    setQcParams(p => ({
      ...p,
      ...(m.total_counts    ? { min_counts: Math.ceil(m.total_counts.mad_min) }     : {}),
      ...(m.n_genes_by_counts ? {
        min_genes: Math.ceil(m.n_genes_by_counts.mad_min),
        max_genes: Math.ceil(m.n_genes_by_counts.mad_max),
      } : {}),
    }))
  }

  // ── 處理函式 ──
  const handleRunQC = async () => {
    // 若下游已有 UMAP/標註結果，提示使用者重跑會使其失效
    if (umapSt?.status === 'done' || heatSt?.status === 'done' || labelApplied) {
      const ok = window.confirm(
        '⚠️ 已有 UMAP / 熱圖 / 標註結果。\n\n重跑 QC 將使下游結果失效，需重新執行 UMAP → 熱圖 → 標註。\n\n確定繼續？'
      )
      if (!ok) return
    }
    setQcImages({})
    setUmapImages({})
    setHeatmapImages({})
    setRoiOverlays({})
    setClusterLabels({})
    setClusterMeta({})
    setLabelApplied(false)
    updateStage('analysis', { status: 'running', progress: 0, message: 'QC 前處理中...' })
    const mergeFlag = analysisMode === 'merge' && hasMultipleRois
    await runQC({
      ...qcParams,
      merge_rois: mergeFlag,
      roi_name: mergeFlag ? undefined : (selectedRoi || undefined),
      input_source: inputSource,
    })
  }

  const handleRunUMAP = async () => {
    const resolutions = parseResolutions()
    if (!resolutions.length) return
    setUmapImages({})        // 清空舊圖
    setHeatmapImages({})
    updateStage('analysis', { status: 'running', progress: 0, message: 'UMAP 計算中...' })
    await runUMAPExplore({ ...umapParams, resolutions })
  }

  const handleRunHeatmap = async () => {
    if (!selectedRes) return
    setHeatmapImages({})     // 清空舊圖
    updateStage('analysis', { status: 'running', progress: 0, message: '熱圖產生中...' })
    await runHeatmap({ resolution: parseFloat(selectedRes), n_top_genes: nTopGenes, n_heatmap_genes: nHeatmapGenes })
  }

  const handleRunAnnotate = async () => {
    if (!annotateRes) return
    // 先清空舊的 meta，使 UI 進入等待狀態
    setClusterMeta({})
    await runAnnotate({
      resolution: parseFloat(annotateRes),
      model_name: annotateModel,
      mode: annotateMode,
      immune_conf_threshold: immuneConfThreshold,
      score_threshold: scoreThreshold,
      uncertain_threshold: uncertainThreshold,
      enable_tier3: enableTier3,
      tier3_conf_threshold: tier3ConfThreshold,
    })
    // 強制 refetch：解決第二次點擊時 refetchInterval 因 status='done' 停止 poll 的問題
    queryClient.invalidateQueries({ queryKey: ['annotate_status'] })
  }

  const handleApplyLabels = async () => {
    if (!annotateRes || !Object.keys(clusterLabels).length) return
    setApplyLabelError('')
    try {
      const r = await applyLabels({ resolution: parseFloat(annotateRes), labels: clusterLabels })
      if (r.data?.status === 'ok') {
        setLabelApplied(true)
      } else {
        setApplyLabelError(r.data?.message ?? '套用失敗')
      }
    } catch (e: unknown) {
      setApplyLabelError(e instanceof Error ? e.message : '套用失敗')
    }
  }

  const handleAnnotateResChange = (res: string) => {
    setAnnotateRes(res)
    setClusterLabels({})
    setClusterMeta({})
    setLabelApplied(false)
    setApplyLabelError('')
    getClusterInfo(parseFloat(res)).then(r => {
      if (r.data?.data) {
        const { cluster_ids, existing_labels } = r.data.data
        const init: Record<string, string> = {}
        cluster_ids.forEach((id: string) => { init[id] = existing_labels[id] ?? '' })
        setClusterLabels(init)
        if (Object.values(existing_labels).some((v: unknown) => v)) setLabelApplied(true)
      }
    }).catch((e: unknown) => {
      setApplyLabelError(`載入 cluster 資訊失敗：${e instanceof Error ? e.message : String(e)}`)
    })
  }

  const qcDone   = qcSt?.status   === 'done'
  const umapDone = umapSt?.status === 'done'

  return (
    <div className="space-y-4">

      {/* ═══════════════════════════════════════════
          區塊 1：QC 前處理
      ═══════════════════════════════════════════ */}
      <div className="bg-surface-card rounded-xl border border-surface-border p-5">
        <div className="flex items-center justify-between mb-1">
          <SectionHeader
            step={1}
            title="QC 前處理"
            subtitle="QC → normalize → HVG → PCA"
          />
          <div className="flex items-center gap-4">
            <StatusBadge
              status={qcSt?.status ?? 'idle'}
              message={qcSt?.message ?? '尚未執行'}
            />
            <RunButton label="執行前處理" onClick={handleRunQC} status={qcSt?.status ?? 'idle'} />
          </div>
        </div>

        {/* 分析來源選擇 */}
        <div className="bg-surface-darker rounded-lg border border-surface-border p-4 mt-4">
          <p className="text-xs text-gray-400 mb-3 font-medium">分析來源</p>
          <div className="flex flex-col gap-2">
            <label className="flex items-center gap-3 cursor-pointer">
              <input
                type="radio"
                name="analysisMode"
                value="single"
                checked={analysisMode === 'single'}
                onChange={() => setAnalysisMode('single')}
                className="accent-brand-primary"
              />
              <span className="text-sm text-gray-200">單一 ROI</span>
              {analysisMode === 'single' && availableRois.length > 0 && (
                <select
                  value={selectedRoi}
                  onChange={e => setSelectedRoi(e.target.value)}
                  className="ml-2 bg-surface-card border border-gray-600 rounded px-2 py-0.5 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
                >
                  {availableRois.map(r => (
                    <option key={r.name} value={r.name} disabled={!r.available}>
                      {r.name}{r.available ? '' : ' (尚未執行)'}
                    </option>
                  ))}
                </select>
              )}
            </label>
            <label className={`flex items-center gap-3 ${hasMultipleRois ? 'cursor-pointer' : 'opacity-40 cursor-not-allowed'}`}>
              <input
                type="radio"
                name="analysisMode"
                value="merge"
                checked={analysisMode === 'merge'}
                onChange={() => setAnalysisMode('merge')}
                disabled={!hasMultipleRois}
                className="accent-brand-primary"
              />
              <span className="text-sm text-gray-200">
                合併所有 ROI
                {hasMultipleRois
                  ? <span className="text-xs text-gray-400 ml-2">（{availableRois.filter(r => r.available).length} 個可用，同一 H&E + Visium，無需 batch correction）</span>
                  : <span className="text-xs text-gray-500 ml-2">（需 ≥ 2 個 ROI）</span>
                }
              </span>
            </label>
          </div>
        </div>

        {/* RNA 計數來源選擇 */}
        <div className="bg-surface-darker rounded-lg border border-surface-border p-4 mt-3">
          <p className="text-xs text-gray-400 mb-2 font-medium">RNA 計數來源</p>
          <div className="flex gap-4">
            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="radio"
                name="inputSource"
                value="cellpose"
                checked={inputSource === 'cellpose'}
                onChange={() => {
                  setInputSource('cellpose')
                  if (histData) handleLoadHist('cellpose')
                }}
                className="accent-brand-primary"
              />
              <span className="text-sm text-gray-200">Cellpose 直接計數</span>
              <span className="text-[11px] text-gray-500 ml-1">（cellpose_cells.h5ad）</span>
            </label>
            <label className={`flex items-center gap-2 ${hasProsegRois ? 'cursor-pointer' : 'opacity-40 cursor-not-allowed'}`}>
              <input
                type="radio"
                name="inputSource"
                value="proseg"
                checked={inputSource === 'proseg'}
                onChange={() => {
                  setInputSource('proseg')
                  if (histData) handleLoadHist('proseg')
                }}
                disabled={!hasProsegRois}
                className="accent-purple-500"
              />
              <span className="text-sm text-gray-200">Proseg RNA 重分配</span>
              <span className="text-[11px] text-gray-500 ml-1">（proseg_cells.h5ad）</span>
              {!hasProsegRois && (
                <span className="text-[11px] text-yellow-600 ml-1">← 請先完成 Stage 2.5</span>
              )}
            </label>
          </div>
        </div>

        {/* ── 原始分布預覽 ── */}
        <div className="bg-surface-darker rounded-lg border border-surface-border p-4 mt-4">
          <div className="flex items-center justify-between mb-3">
            <div>
              <p className="text-xs font-medium text-gray-300">原始分布預覽</p>
              <p className="text-xs text-gray-500 mt-0.5">
                觀察細胞數量分布後再設定 QC 閾值｜
                <span className="text-red-400">紅線</span> = min，
                <span className="text-orange-400">橙線</span> = max，
                <span className="text-green-500">綠虛線</span> = MAD建議
              </p>
            </div>
            <div className="flex items-center gap-2 shrink-0">
              {histData && (
                <button
                  onClick={handleApplyMad}
                  className="px-3 py-1 rounded text-xs font-medium bg-green-800/40 border border-green-700 text-green-300 hover:bg-green-700/50 transition-colors"
                >
                  套用 MAD 建議值
                </button>
              )}
              <button
                onClick={() => handleLoadHist()}
                disabled={histLoading}
                className="px-3 py-1.5 rounded-lg text-xs font-medium bg-surface-card border border-gray-600 text-gray-300 hover:border-brand-primary hover:text-brand-primary disabled:opacity-50 transition-colors"
              >
                {histLoading ? '載入中...' : histData ? '重新載入' : '載入原始分布'}
              </button>
            </div>
          </div>

          {histError && (
            <p className="text-xs text-red-400 mb-2">{histError}</p>
          )}

          {histData && (
            <div className={`grid gap-6 ${Object.keys(histData.metrics).length >= 4 ? 'grid-cols-1 md:grid-cols-2 lg:grid-cols-3' : 'grid-cols-1 md:grid-cols-3'}`}>
              {/* Transcripts Per Cell */}
              {histData.metrics.total_counts && (
                <QcHistogram
                  metric={histData.metrics.total_counts}
                  minVal={qcParams.min_counts}
                  maxVal={null}
                  showMin={true}
                  showMax={false}
                  onMinChange={v => setQcParams(p => ({ ...p, min_counts: v ?? 0 }))}
                  onMaxChange={() => {}}
                  logScale={!!logScales.total_counts}
                  onLogScaleToggle={() => setLogScales(s => ({ ...s, total_counts: !s.total_counts }))}
                  totalCells={histData.n_cells}
                />
              )}

              {/* Genes Per Cell */}
              {histData.metrics.n_genes_by_counts && (
                <QcHistogram
                  metric={histData.metrics.n_genes_by_counts}
                  minVal={qcParams.min_genes}
                  maxVal={qcParams.max_genes}
                  showMin={true}
                  showMax={true}
                  onMinChange={v => setQcParams(p => ({ ...p, min_genes: v ?? 0 }))}
                  onMaxChange={v => setQcParams(p => ({ ...p, max_genes: v ?? 99999 }))}
                  logScale={!!logScales.n_genes_by_counts}
                  onLogScaleToggle={() => setLogScales(s => ({ ...s, n_genes_by_counts: !s.n_genes_by_counts }))}
                  totalCells={histData.n_cells}
                />
              )}

              {/* Cell Size (Proseg only) */}
              {histData.metrics.cell_area && (
                <QcHistogram
                  metric={histData.metrics.cell_area}
                  minVal={null}
                  maxVal={null}
                  showMin={false}
                  showMax={false}
                  showMad={true}
                  onMinChange={() => {}}
                  onMaxChange={() => {}}
                  logScale={!!logScales.cell_area}
                  onLogScaleToggle={() => setLogScales(s => ({ ...s, cell_area: !s.cell_area }))}
                  totalCells={histData.n_cells}
                />
              )}

              {/* Complexity Score */}
              {histData.metrics.complexity && (
                <QcHistogram
                  metric={histData.metrics.complexity}
                  minVal={qcParams.min_complexity}
                  maxVal={null}
                  showMin={true}
                  showMax={false}
                  showMad={true}
                  onMinChange={v => setQcParams(p => ({ ...p, min_complexity: v ?? 0 }))}
                  onMaxChange={() => {}}
                  logScale={!!logScales.complexity}
                  onLogScaleToggle={() => setLogScales(s => ({ ...s, complexity: !s.complexity }))}
                  totalCells={histData.n_cells}
                />
              )}

              {/* Mitochondrial % */}
              {histData.metrics.pct_counts_mt && (
                <QcHistogram
                  metric={histData.metrics.pct_counts_mt}
                  minVal={null}
                  maxVal={qcParams.max_pct_mito}
                  showMin={false}
                  showMax={true}
                  showMad={true}
                  onMinChange={() => {}}
                  onMaxChange={v => setQcParams(p => ({ ...p, max_pct_mito: v ?? 100 }))}
                  logScale={!!logScales.pct_counts_mt}
                  onLogScaleToggle={() => setLogScales(s => ({ ...s, pct_counts_mt: !s.pct_counts_mt }))}
                  totalCells={histData.n_cells}
                />
              )}
            </div>
          )}
        </div>

        {/* QC 參數 */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-3 bg-surface-darker p-4 rounded-lg border border-surface-border mt-4">
          <NumberField label="最低基因數 (min_genes)" value={qcParams.min_genes}
            onChange={v => setQcParams(p => ({ ...p, min_genes: v }))} min={0} />
          <NumberField label="最高基因數 (max_genes)" value={qcParams.max_genes}
            onChange={v => setQcParams(p => ({ ...p, max_genes: v }))} min={0} />
          <NumberField label="最低 UMI 數 (min_counts)" value={qcParams.min_counts}
            onChange={v => setQcParams(p => ({ ...p, min_counts: v }))} min={0} />
          <NumberField label="最低細胞數/基因 (min_cells)" value={qcParams.min_cells}
            onChange={v => setQcParams(p => ({ ...p, min_cells: v }))} min={0} />
          <NumberField label="粒線體上限 % (max_pct_mito)" value={qcParams.max_pct_mito}
            onChange={v => setQcParams(p => ({ ...p, max_pct_mito: v }))} step={0.5} min={0} />
          <NumberField label="最低複雜度 (min_complexity)" value={qcParams.min_complexity}
            onChange={v => setQcParams(p => ({ ...p, min_complexity: v }))} step={0.01} min={0} hint="建議 0.8 以上" />
          <NumberField label="HVG 數量 (n_top_genes)" value={qcParams.n_top_genes}
            onChange={v => setQcParams(p => ({ ...p, n_top_genes: v }))} min={100} />
          <NumberField label="PCA 維度數 (n_pcs)" value={qcParams.n_pcs}
            onChange={v => setQcParams(p => ({ ...p, n_pcs: v }))} min={10} />
        </div>

        {/* QC 圖表 */}
        {Object.keys(qcImages).length > 0 && (
          <>
            <ChartView
              images={qcImages}
              tabs={[
                { key: 'violin',         label: '小提琴圖 (QC 分布)' },
                { key: 'scatter',        label: '散佈圖 (UMI vs Genes)' },
                { key: 'elbow',          label: 'PCA Elbow' },
                { key: 'pre_qc',         label: 'H&E 疊圖 (QC 前)' },
                { key: 'post_qc',        label: 'H&E 疊圖 (QC 後)' },
                { key: 'roi_comparison', label: 'ROI 輪廓比較（全部）' },
              ]}
            />
            {/* HD 存檔下載 */}
            {(qcImages['pre_qc'] || qcImages['post_qc']) && (
              <div className="mt-2 flex gap-3">
                <span className="text-xs text-gray-400 self-center">HD 存檔 (300 DPI)：</span>
                {qcImages['pre_qc'] && (
                  <a
                    href={getOverlayHdUrl('pre_qc')}
                    download="overlay_pre_qc_hd.png"
                    className="text-xs text-brand-primary hover:underline"
                  >
                    下載 Pre-QC HD
                  </a>
                )}
                {qcImages['post_qc'] && (
                  <a
                    href={getOverlayHdUrl('post_qc')}
                    download="overlay_post_qc_hd.png"
                    className="text-xs text-brand-primary hover:underline"
                  >
                    下載 Post-QC HD
                  </a>
                )}
              </div>
            )}
          </>
        )}

        {/* ROI 輪廓比較（QC 完成後出現，多 ROI 或單 ROI 均支援）*/}
        {Object.keys(roiOverlays).length > 0 && (
          <RoiContourPanel data={roiOverlays} />
        )}
      </div>

      {/* ═══════════════════════════════════════════
          區塊 2：UMAP 解析
      ═══════════════════════════════════════════ */}
      <div className={`bg-surface-card rounded-xl border border-surface-border p-5 transition-opacity ${!qcDone ? 'opacity-50' : ''}`}>
        <div className="flex items-center justify-between mb-1">
          <SectionHeader
            step={2}
            title="UMAP 解析"
            subtitle="建立 KNN 圖 → UMAP → Leiden（支援多組解析度同時比較）"
          />
          <div className="flex items-center gap-4">
            <StatusBadge
              status={umapSt?.status ?? 'idle'}
              message={umapSt?.message ?? '尚未執行'}
            />
            <RunButton
              label="執行 UMAP"
              onClick={handleRunUMAP}
              status={umapSt?.status ?? 'idle'}
              disabled={!qcDone}
            />
          </div>
        </div>

        {/* UMAP 參數 */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-3 bg-surface-darker p-4 rounded-lg border border-surface-border mt-4">
          <NumberField label="KNN 鄰居數 (n_neighbors)" value={umapParams.n_neighbors}
            onChange={v => setUmapParams(p => ({ ...p, n_neighbors: v }))} min={2}
            hint="建議：15–30；細胞多可調至 30–50" />
          <NumberField label="使用 PC 數 (n_pcs)" value={umapParams.n_pcs}
            onChange={v => setUmapParams(p => ({ ...p, n_pcs: v }))} min={5}
            hint="建議：20–30；參考 Elbow 圖選擇" />
          <NumberField label="最小距離 (min_dist)" value={umapParams.min_dist}
            onChange={v => setUmapParams(p => ({ ...p, min_dist: v }))} step={0.05} min={0.01}
            hint="建議：0.1–0.3；越小叢集越緊密" />
          <div>
            <label className="block text-xs text-gray-400 mb-1">解析度列表（逗號分隔）</label>
            <input
              type="text"
              className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
              value={resolutionInput}
              placeholder="例：0.3, 0.5, 0.8, 1.2"
              onChange={e => setResolutionInput(e.target.value)}
            />
          </div>
        </div>

        {/* UMAP 圖表:個別 + Grid */}
        {Object.keys(umapImages).length > 0 && (() => {
          const resTabs = availableResolutions.map(r => ({ key: r, label: `Res = ${r}` }))
          const allTabs = [...resTabs, { key: 'grid', label: '全覽 Grid' }]
          return (
            <ChartView images={umapImages} tabs={allTabs} fullWidthKeys={['grid']} />
          )
        })()}
      </div>

      {/* ═══════════════════════════════════════════
          區塊 3：熱圖輸出
      ═══════════════════════════════════════════ */}
      <div className={`bg-surface-card rounded-xl border border-surface-border p-5 transition-opacity ${!umapDone ? 'opacity-50' : ''}`}>
        <div className="flex items-center justify-between mb-1">
          <SectionHeader
            step={3}
            title="熱圖輸出"
            subtitle="觀察 Heatmap + Dotplot 了解各 cluster 的 marker gene，再進行細胞類型標註"
          />
          <div className="flex items-center gap-4">
            <StatusBadge
              status={heatSt?.status ?? 'idle'}
              message={heatSt?.message ?? '尚未執行'}
            />
            <RunButton
              label="產生熱圖"
              onClick={handleRunHeatmap}
              status={heatSt?.status ?? 'idle'}
              disabled={!umapDone || !selectedRes}
            />
          </div>
        </div>

        {/* 熱圖參數 */}
        <div className="grid grid-cols-2 md:grid-cols-3 gap-3 bg-surface-darker p-4 rounded-lg border border-surface-border mt-4">
          <div>
            <label className="block text-xs text-gray-400 mb-1">解析度選擇</label>
            <select
              className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
              value={selectedRes}
              onChange={e => setSelectedRes(e.target.value)}
              disabled={!availableResolutions.length}
            >
              {availableResolutions.length === 0 && (
                <option value="">（請先執行 UMAP）</option>
              )}
              {availableResolutions.map(r => (
                <option key={r} value={r}>Resolution = {r}</option>
              ))}
            </select>
          </div>
          <NumberField
            label="熱圖基因數 (n_heatmap_genes)"
            value={nHeatmapGenes}
            onChange={setNHeatmapGenes}
            min={10}
          />
          <NumberField
            label="Dotplot 每群基因數 (n_top_genes)"
            value={nTopGenes}
            onChange={setNTopGenes}
            min={5}
          />
        </div>

        {/* 熱圖 + 點圖圖表 */}
        {Object.keys(heatmapImages).length > 0 && (
          <ChartView
            images={heatmapImages}
            tabs={[
              { key: 'heatmap', label: 'Heatmap（細胞等級分布）' },
              { key: 'dotplot', label: 'Dotplot（表達比例 + 強度）' },
            ]}
          />
        )}
      </div>

      {/* ═══════════════════════════════════════════
          區塊 4（原3）：細胞類型標註
      ═══════════════════════════════════════════ */}
      <div className={`bg-surface-card rounded-xl border border-surface-border p-5 transition-opacity ${!umapDone ? 'opacity-50' : ''}`}>
        <div className="flex items-center justify-between mb-1">
          <SectionHeader
            step={4}
            title="細胞類型標註"
            subtitle="觀察 Heatmap/Dotplot 後，以 CellTypist 自動建議或手動填寫標籤；套用後可重跑熱圖令 y 軸顯示細胞名稱"
          />
          <div className="flex items-center gap-4">
            {labelApplied && (
              <span className="text-xs text-green-400 font-medium">✓ 標籤已套用</span>
            )}
          </div>
        </div>

        {/* 解析度 + 模型 + 模式設定 */}
        <div className="grid grid-cols-2 md:grid-cols-3 gap-3 bg-surface-darker p-4 rounded-lg border border-surface-border mt-4">
          <div>
            <label className="block text-xs text-gray-400 mb-1">標註解析度</label>
            <select
              className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
              value={annotateRes}
              onChange={e => handleAnnotateResChange(e.target.value)}
              disabled={!availableResolutions.length}
            >
              {availableResolutions.length === 0 && <option value="">（請先執行 UMAP）</option>}
              {availableResolutions.map(r => (
                <option key={r} value={r}>Resolution = {r}</option>
              ))}
            </select>
          </div>

          <div>
            <label className="block text-xs text-gray-400 mb-1">組織模型（CRC / 器官）</label>
            <select
              className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
              value={annotateModel}
              onChange={e => setAnnotateModel(e.target.value)}
            >
              {Object.entries(celltypistModels).map(([label, filename]) => (
                <option key={filename} value={filename}>{label}</option>
              ))}
              {Object.keys(celltypistModels).length === 0 && (
                <option value="Human_Colorectal_Cancer.pkl">Human CRC（大腸癌）</option>
              )}
            </select>
          </div>

          <div>
            <label className="block text-xs text-gray-400 mb-1">標註模式</label>
            <select
              className="w-full bg-surface-highlight border border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary"
              value={annotateMode}
              onChange={e => setAnnotateMode(e.target.value as 'dual' | 'single')}
            >
              <option value="dual">雙模型（免疫 + 組織）— 推薦</option>
              <option value="single">單模型（僅組織模型）</option>
            </select>
          </div>

          {annotateMode === 'dual' && (<>
            <div>
              <label className="block text-xs text-gray-400 mb-1">
                免疫識別閾值：<span className="text-gray-200">{immuneConfThreshold.toFixed(2)}</span>
                <span className="ml-1 text-gray-600">（低於此值→組織模型）</span>
              </label>
              <input
                type="range" min="0.1" max="0.9" step="0.05"
                value={immuneConfThreshold}
                onChange={e => setImmuneConfThreshold(parseFloat(e.target.value))}
                className="w-full accent-brand-primary"
              />
              <div className="flex justify-between text-xs text-gray-600 mt-0.5">
                <span>寬鬆</span><span>嚴格</span>
              </div>
            </div>
            <div>
              <label className="block text-xs text-gray-400 mb-1">
                信心警告閾值：<span className="text-gray-200">{uncertainThreshold.toFixed(2)}</span>
                <span className="ml-1 text-gray-600">（低於此值→橘色警告）</span>
              </label>
              <input
                type="range" min="0.3" max="0.95" step="0.05"
                value={uncertainThreshold}
                onChange={e => setUncertainThreshold(parseFloat(e.target.value))}
                className="w-full accent-orange-400"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-400 mb-1">
                基因評分閾值：<span className="text-gray-200">{scoreThreshold.toFixed(2)}</span>
                <span className="ml-1 text-gray-600">（低於此值→不附加功能狀態）</span>
              </label>
              <input
                type="range" min="0.1" max="0.8" step="0.05"
                value={scoreThreshold}
                onChange={e => setScoreThreshold(parseFloat(e.target.value))}
                className="w-full accent-purple-400"
              />
            </div>
            {/* Tier 3 開關 */}
            <div className="col-span-full border-t border-gray-700 pt-3 mt-1">
              <div className="flex items-center gap-3 mb-2">
                <input
                  type="checkbox"
                  id="tier3-toggle"
                  checked={enableTier3}
                  onChange={e => setEnableTier3(e.target.checked)}
                  className="w-4 h-4 accent-indigo-400 cursor-pointer"
                />
                <label htmlFor="tier3-toggle" className="text-xs text-gray-300 cursor-pointer">
                  <span className="font-medium text-indigo-300">Tier 3：啟用精細免疫亞型</span>
                  <span className="ml-2 text-gray-500">（Immune_All_Low.pkl — 98 種亞型，首次使用需下載）</span>
                </label>
              </div>
              {enableTier3 && (
                <div className="ml-7">
                  <label className="block text-xs text-gray-400 mb-1">
                    Tier3 替換閾值：<span className="text-gray-200">{tier3ConfThreshold.toFixed(2)}</span>
                    <span className="ml-1 text-gray-600">（低於此值只顯示不替換標籤）</span>
                  </label>
                  <input
                    type="range" min="0.3" max="0.9" step="0.05"
                    value={tier3ConfThreshold}
                    onChange={e => setTier3ConfThreshold(parseFloat(e.target.value))}
                    className="w-full accent-indigo-400"
                  />
                  <p className="text-xs text-yellow-600 mt-1">
                    ⚠ Tier 3 會對每個免疫 cluster 額外跑一次模型，標註時間約增加 50%
                  </p>
                </div>
              )}
            </div>
          </>)}

          <div className={`flex flex-col justify-end ${annotateMode === 'dual' ? '' : 'col-start-3'}`}>
            <RunButton
              label="CellTypist 自動標註"
              onClick={handleRunAnnotate}
              status={annotSt?.status ?? 'idle'}
              disabled={!umapDone || !annotateRes}
            />
            {annotSt?.message && (
              <p className="text-xs text-gray-500 mt-1">{annotSt.message}</p>
            )}
          </div>
        </div>

        {/* Cluster 標籤表格 */}
        {Object.keys(clusterLabels).length > 0 && (
          <div className="mt-4 bg-surface-darker rounded-lg border border-surface-border p-4">
            <div className="flex items-center justify-between mb-3">
              <p className="text-xs text-gray-400 font-medium">
                Cluster 標籤（共 {Object.keys(clusterLabels).length} 個）
                <span className="ml-2 text-gray-500">— 可直接編輯，標註來源以色塊標示</span>
              </p>
              <div className="flex flex-col items-end gap-1">
                <button
                  onClick={handleApplyLabels}
                  disabled={!umapDone || !annotateRes}
                  className="px-4 py-1.5 rounded-lg text-sm font-medium bg-primary text-white hover:bg-primary-dark disabled:bg-gray-700 disabled:text-gray-500 disabled:cursor-not-allowed transition-colors"
                >
                  套用標籤
                </button>
                {applyLabelError && (
                  <p className="text-xs text-red-400">{applyLabelError}</p>
                )}
              </div>
            </div>

            {/* 圖例 */}
            {Object.keys(clusterMeta).length > 0 && (
              <div className="flex flex-wrap gap-x-4 gap-y-1 mb-3 text-xs text-gray-500">
                <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-purple-400 inline-block"/>免疫（Immune_All_High）</span>
                <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-sky-400 inline-block"/>組織模型</span>
                <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-orange-400 inline-block"/>需確認（信心不足）</span>
                <span className="ml-2">信心：<span className="text-green-400">●</span>≥{uncertainThreshold.toFixed(1)} / <span className="text-yellow-400">●</span>中 / <span className="text-red-400">●</span>低</span>
              </div>
            )}

            <div className="grid grid-cols-1 md:grid-cols-2 gap-2 max-h-80 overflow-y-auto pr-1">
              {Object.entries(clusterLabels)
                .sort(([a], [b]) => parseInt(a) - parseInt(b))
                .map(([clusterId, label]) => {
                  const meta = clusterMeta[clusterId]
                  const conf = meta?.confidence ?? 0
                  const source = meta?.source ?? ''
                  const uncertain = meta?.uncertain ?? false
                  const state = meta?.state
                  const stateScore = meta?.state_score ?? 0

                  const sourceDot = uncertain
                    ? 'bg-orange-400'
                    : source === 'immune' ? 'bg-purple-400'
                    : source === 'crc' ? 'bg-sky-400'
                    : 'bg-gray-500'

                  const confColor = conf >= uncertainThreshold
                    ? 'text-green-400'
                    : conf >= 0.5 ? 'text-yellow-400'
                    : 'text-red-400'

                  const rowBorder = uncertain
                    ? 'border border-orange-900/50 bg-orange-950/20 rounded'
                    : ''

                  const tooltipLines = [
                    meta?.tier3_label ? `Tier3 精細亞型: ${meta.tier3_label} (conf=${(meta.tier3_conf ?? 0).toFixed(2)})` : '',
                    meta?.immune_label ? `Tier1 免疫: ${meta.immune_label} (${(meta.immune_conf ?? 0).toFixed(2)})` : '',
                    meta?.crc_label ? `Tier1 組織: ${meta.crc_label} (${(meta.crc_conf ?? 0).toFixed(2)})` : '',
                    state ? `Tier2 功能狀態: ${state} (score=${stateScore.toFixed(3)})` : '',
                    uncertain ? '⚠ 信心不足，建議手動確認' : '',
                  ].filter(Boolean).join('\n')

                  return (
                    <div key={clusterId} className={`flex items-center gap-2 px-1 py-0.5 ${rowBorder}`}>
                      <span className="text-xs text-gray-400 shrink-0 w-8 text-right font-mono">
                        {clusterId}
                      </span>
                      {meta && (
                        <span
                          className={`w-2 h-2 rounded-full shrink-0 ${sourceDot} cursor-help`}
                          title={tooltipLines}
                        />
                      )}
                      <input
                        type="text"
                        value={label}
                        placeholder="輸入細胞類型..."
                        className={`flex-1 bg-surface-highlight border rounded px-2 py-0.5 text-xs text-gray-200 focus:outline-none focus:border-brand-primary ${
                          uncertain ? 'border-orange-700' : 'border-gray-600'
                        }`}
                        onChange={e => setClusterLabels(prev => ({ ...prev, [clusterId]: e.target.value }))}
                      />
                      {/* 功能狀態 badge */}
                      {state && (
                        <span className="text-xs px-1 py-0.5 rounded bg-purple-900/60 text-purple-300 shrink-0 font-mono" title={`基因評分: ${stateScore.toFixed(3)}`}>
                          {state.replace('_', ' ')}
                        </span>
                      )}
                      {/* 信心分數 */}
                      {meta && !state && (
                        <span className={`text-xs font-mono shrink-0 w-10 text-right ${confColor}`} title={tooltipLines}>
                          {conf.toFixed(2)}
                        </span>
                      )}
                      {/* uncertain 警告圖示 */}
                      {uncertain && (
                        <span className="text-orange-400 shrink-0 text-xs" title="信心不足，建議手動確認">⚠</span>
                      )}
                    </div>
                  )
                })}
            </div>
          </div>
        )}
      </div>

      <Terminal stage="analysis" />
    </div>
  )
}

