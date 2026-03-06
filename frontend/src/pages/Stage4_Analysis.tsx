import { useState, useEffect, useCallback } from 'react'
import { useQuery } from '@tanstack/react-query'
import { usePipelineStore } from '../stores/pipelineStore'
import Terminal from '../components/shared/Terminal'
import useStageLog from '../hooks/useStageLog'
import {
  runQC, getQCStatus, getQCImages,
  runUMAPExplore, getUMAPExploreStatus, getUMAPImages,
  runHeatmap, getHeatmapStatus, getHeatmapImage,
  getConfig, getOverlayHdUrl,
} from '../api/client'

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
  label, value, onChange, step = 1, min,
}: { label: string; value: number; onChange: (v: number) => void; step?: number; min?: number }) {
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
    </div>
  )
}

function ChartView({ images, tabs, compact = false }: { images: Record<string, string>; tabs: { key: string; label: string }[]; compact?: boolean }) {
  const available = tabs.filter(t => images[t.key])
  const [active, setActive] = useState(available[0]?.key ?? '')
  useEffect(() => { if (available.length && !images[active]) setActive(available[0].key) }, [images])
  if (!available.length) return null

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
          className={`rounded-lg ${compact ? 'max-w-[50%]' : 'max-w-full'}`}
          alt={active}
        />
      )}
    </div>
  )
}

// UMAP 專用：個別解析度圖縮小為 50%，Grid 全寬
function UMAPChartView({ images, tabs }: { images: Record<string, string>; tabs: { key: string; label: string }[] }) {
  const available = tabs.filter(t => images[t.key])
  const [active, setActive] = useState(available[0]?.key ?? '')
  useEffect(() => { if (available.length && !images[active]) setActive(available[0].key) }, [images])
  if (!available.length) return null

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
          className={`rounded-lg ${active === 'grid' ? 'max-w-full' : 'max-w-[50%]'}`}
          alt={active}
        />
      )}
    </div>
  )
}

// ── 主元件 ────────────────────────────────────────────────────────

export default function Stage4_Analysis() {
  useStageLog('analysis')
  const { updateStage } = usePipelineStore()

  // ── QC 參數 ──
  const [qcParams, setQcParams] = useState({
    min_genes: 10, max_genes: 8000,
    min_counts: 5, min_cells: 3,
    max_pct_mito: 80.0, n_top_genes: 2000, n_pcs: 30,
  })

  // ── UMAP 參數 ──
  const [umapParams, setUmapParams] = useState({
    n_pcs: 30, n_neighbors: 15, min_dist: 0.3,
  })
  const [resolutionInput, setResolutionInput] = useState('0.3, 0.5, 0.8')

  // ── Heatmap 參數 ──
  const [selectedRes, setSelectedRes] = useState<string>('')
  const [nTopGenes, setNTopGenes] = useState(20)

  // ── 儲存圖表 ──
  const [qcImages, setQcImages] = useState<Record<string, string>>({})
  const [umapImages, setUmapImages] = useState<Record<string, string>>({})
  const [heatmapImages, setHeatmapImages] = useState<Record<string, string>>({})

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

  // ── TanStack Query: 三步驟 status 輪詢 ──
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

  // ── 各步驟完成後自動拉取圖表（包含頁面初次載入 / 後端重啟後恢復） ──
  useEffect(() => {
    if (qcSt?.status === 'done') {
      getQCImages().then(r => { if (r.data.data) setQcImages(r.data.data) })
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
    }
    if (heatSt?.status === 'done') {
      updateStage('analysis', { status: 'done', progress: 1, message: '分析完成' })
    }
  }, [heatSt?.status])

  // ── 解析 resolution 文字輸入 ──
  const parseResolutions = useCallback((): number[] => {
    return resolutionInput
      .split(/[,\s]+/)
      .map(s => parseFloat(s.trim()))
      .filter(v => !isNaN(v) && v > 0)
  }, [resolutionInput])

  // ── 可用的 resolution 列表（供 Heatmap 下拉） ──
  const availableResolutions = Object.keys(umapImages).filter(k => k !== 'grid').sort()

  // ── 處理函式 ──
  const handleRunQC = async () => {
    setQcImages({})          // 清空舊圖，等新結果回來再重新拉取
    setUmapImages({})        // QC 重跑時下游也要重置
    setHeatmapImages({})
    updateStage('analysis', { status: 'running', progress: 0, message: 'QC 前處理中...' })
    await runQC(qcParams)
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
    await runHeatmap({ resolution: parseFloat(selectedRes), n_top_genes: nTopGenes })
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
                { key: 'violin',  label: '小提琴圖 (QC 分布)' },
                { key: 'scatter', label: '散佈圖 (UMI vs Genes)' },
                { key: 'elbow',   label: 'PCA Elbow' },
                { key: 'pre_qc',  label: 'H&E 疊圖 (QC 前)' },
                { key: 'post_qc', label: 'H&E 疊圖 (QC 後)' },
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
            onChange={v => setUmapParams(p => ({ ...p, n_neighbors: v }))} min={2} />
          <NumberField label="使用 PC 數 (n_pcs)" value={umapParams.n_pcs}
            onChange={v => setUmapParams(p => ({ ...p, n_pcs: v }))} min={5} />
          <NumberField label="最小距離 (min_dist)" value={umapParams.min_dist}
            onChange={v => setUmapParams(p => ({ ...p, min_dist: v }))} step={0.05} min={0.01} />
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
            <UMAPChartView images={umapImages} tabs={allTabs} />
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
            subtitle="選擇解析度 → rank_genes_groups → heatmap + dotplot（y 軸：cluster，x 軸：強表現基因）"
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
            label="每群顯示基因數 (n_top_genes)"
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

      <Terminal stage="analysis" />
    </div>
  )
}

