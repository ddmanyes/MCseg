import { useState, useEffect, useRef, useMemo, Component, type ReactNode } from 'react'
import { useQuery } from '@tanstack/react-query'
import { getSpatialGeneList, postSpatialGenePlot, getAvailableRois } from '../api/client'
import { useT } from '../i18n'

class ErrorBoundary extends Component<{ children: ReactNode }, { error: string }> {
  state = { error: '' }
  static getDerivedStateFromError(e: Error) { return { error: e.message } }
  render() {
    if (this.state.error) return (
      <div className="p-6 text-red-400 text-sm font-mono whitespace-pre-wrap">
        <p className="font-bold mb-2">Render error:</p>
        {this.state.error}
      </div>
    )
    return this.props.children
  }
}

const CMAPS = ['viridis', 'magma', 'plasma', 'inferno', 'Reds', 'Blues', 'YlOrRd']
const MAX_GENES = 4

type PlotMode = 'contour' | 'set'

interface PresetGeneSet {
  label: string
  genes: string[]
}

const PRESET_GENE_SETS: PresetGeneSet[] = [
  { label: 'CD4 T',       genes: ['CD4', 'IL7R', 'CCR7', 'TCF7'] },
  { label: 'CD8 T',       genes: ['CD8A', 'CD8B', 'GZMB', 'PRF1'] },
  { label: 'NK',          genes: ['NKG7', 'GNLY', 'KLRD1', 'NCR1'] },
  { label: 'Macrophage',  genes: ['CD68', 'CSF1R', 'MRC1', 'ITGAM'] },
  { label: 'B cell',      genes: ['CD19', 'MS4A1', 'CD79A', 'CD79B'] },
  { label: 'Tumor',       genes: ['EPCAM', 'KRT8', 'KRT18', 'MUC1'] },
  { label: 'Fibroblast',  genes: ['ACTA2', 'FAP', 'PDGFRA', 'COL1A1'] },
  { label: 'Endothelial', genes: ['PECAM1', 'VWF', 'CDH5', 'CLDN5'] },
]

function SpatialExplorerInner() {
  const t = useT()

  const [selectedRoi, setSelectedRoi] = useState<string>('')
  const [selectedGenes, setSelectedGenes] = useState<string[]>([])
  const [searchText, setSearchText]     = useState('')
  const [showDropdown, setShowDropdown] = useState(false)
  const [mode, setMode]                 = useState<PlotMode>('contour')
  const [setName, setSetName]           = useState('')
  const [pointSize, setPointSize]       = useState(6)
  const [cmap, setCmap]                 = useState('viridis')
  const [alpha, setAlpha]               = useState(0.8)
  const [plotting, setPlotting]         = useState(false)
  const [plotImage, setPlotImage]       = useState<string | null>(null)
  const [plotError, setPlotError]       = useState('')
  const [plotInfo, setPlotInfo]         = useState<{ n_cells?: number } | null>(null)
  const searchRef = useRef<HTMLDivElement>(null)

  // ── Load available ROIs ───────────────────────────────────────
  const { data: roisData } = useQuery({
    queryKey: ['available_rois'],
    queryFn: () => getAvailableRois().then(r => r.data),
    staleTime: 30000,
  })
  const availableRois: string[] = useMemo(() =>
    ((roisData as any)?.data ?? [])
      .filter((r: any) => r.available)
      .map((r: any) => r.name as string),
    [roisData]
  )

  useEffect(() => {
    if (!selectedRoi && availableRois.length > 0) setSelectedRoi(availableRois[0])
  }, [availableRois, selectedRoi])

  // ── Load gene list ────────────────────────────────────────────
  const { data: geneListData, isLoading: geneLoading } = useQuery({
    queryKey: ['spatial_gene_list', selectedRoi],
    queryFn: () => getSpatialGeneList(selectedRoi || undefined).then(r => r.data),
    enabled: true,
    staleTime: 60000,
  })
  const allGenes: string[] = (geneListData as any)?.data?.genes ?? []
  const isMergeMode: boolean = (geneListData as any)?.data?.merge_mode ?? false

  // Search: starts-with priority, then contains, cap at 30
  const filteredGenes = (() => {
    if (searchText.length === 0) return allGenes.slice(0, 30)
    const q = searchText.toLowerCase()
    const starts   = allGenes.filter(g =>  g.toLowerCase().startsWith(q))
    const contains = allGenes.filter(g => !g.toLowerCase().startsWith(q) && g.toLowerCase().includes(q))
    return [...starts, ...contains].slice(0, 30)
  })()

  // Close dropdown on outside click
  useEffect(() => {
    const handler = (e: MouseEvent) => {
      if (searchRef.current && !searchRef.current.contains(e.target as Node))
        setShowDropdown(false)
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [])

  const addGene = (gene: string) => {
    if (selectedGenes.includes(gene) || selectedGenes.length >= MAX_GENES) return
    setSelectedGenes(prev => [...prev, gene])
    setSearchText('')
    setShowDropdown(false)
  }

  const removeGene = (gene: string) => setSelectedGenes(prev => prev.filter(g => g !== gene))

  const applyPreset = (preset: PresetGeneSet) => {
    // Case-insensitive match against actual gene list (handles mouse Title Case)
    const resolved = preset.genes
      .map(g => allGenes.find(ag => ag.toLowerCase() === g.toLowerCase()))
      .filter((g): g is string => g !== undefined)
      .slice(0, MAX_GENES)
    if (resolved.length === 0) return
    setSelectedGenes(resolved)
    setSetName(preset.label)
    setMode('set')
    setPlotImage(null)
  }

  const handleSavePng = () => {
    if (!plotImage) return
    const roiPart = isMergeMode ? 'all_rois' : (selectedRoi || 'roi')
    const filename = [
      roiPart,
      mode === 'set' ? (setName || 'geneset') : selectedGenes.join('_'),
      mode,
    ].join('_') + '.png'
    const a = document.createElement('a')
    a.href = `data:image/png;base64,${plotImage}`
    a.download = filename
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
  }

  const handlePlot = async () => {
    if (selectedGenes.length === 0) return
    setPlotting(true)
    setPlotError('')
    setPlotImage(null)
    try {
      const res = await postSpatialGenePlot({
        roi_name: isMergeMode ? undefined : (selectedRoi || undefined),
        genes: selectedGenes,
        mode,
        set_name: mode === 'set' && setName ? setName : undefined,
        point_size: pointSize,
        cmap,
        alpha,
      })
      const d = res.data
      if (d.status === 'ok') {
        setPlotImage(d.data.image_b64)
        setPlotInfo({ n_cells: d.data.n_cells })
      } else {
        setPlotError(d.message ?? 'Error')
      }
    } catch (e: any) {
      setPlotError(e.response?.data?.detail ?? e.message ?? 'Unknown error')
    } finally {
      setPlotting(false)
    }
  }

  const modeOptions: { value: PlotMode; label: string }[] = [
    { value: 'contour', label: t('spatial.mode.contour') },
    { value: 'set',     label: t('spatial.mode.set') },
  ]

  return (
    <div className="space-y-4 max-w-5xl mx-auto">
      {/* Header */}
      <div>
        <h2 className="text-base font-semibold text-gray-100">{t('spatial.title')}</h2>
        <p className="text-xs text-gray-500 mt-0.5">{t('spatial.subtitle')}</p>
      </div>

      {/* Controls */}
      <div className="bg-surface-card rounded-xl border border-surface-border p-4 space-y-4">

        {/* Row 1: ROI + Gene search */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {/* ROI selector — hidden in merge mode */}
          {!isMergeMode ? (
            <div>
              <label className="block text-xs text-gray-400 mb-1">{t('spatial.roi')}</label>
              <select
                value={selectedRoi}
                onChange={e => { setSelectedRoi(e.target.value); setSelectedGenes([]); setPlotImage(null); setPlotInfo(null) }}
                className="w-full bg-surface border border-surface-border rounded px-3 py-1.5 text-sm text-gray-200 focus:border-primary focus:outline-none"
              >
                {availableRois.length === 0
                  ? <option value="">—</option>
                  : availableRois.map(r => <option key={r} value={r}>{r}</option>)
                }
              </select>
            </div>
          ) : (
            <div className="flex items-center gap-2 px-3 py-1.5 rounded border border-surface-border bg-surface">
              <span className="text-xs text-primary">⊞</span>
              <span className="text-xs text-gray-400">All ROIs (merged) — each ROI shown as a separate panel</span>
            </div>
          )}

          {/* Gene search */}
          <div ref={searchRef} className="relative">
            <label className="block text-xs text-gray-400 mb-1">{t('spatial.gene_search')}</label>
            <input
              type="text"
              value={searchText}
              onChange={e => { setSearchText(e.target.value); setShowDropdown(true) }}
              onFocus={() => setShowDropdown(true)}
              placeholder={geneLoading ? t('spatial.loading_genes') : t('spatial.gene_placeholder')}
              disabled={geneLoading || selectedGenes.length >= MAX_GENES}
              className="w-full bg-surface border border-surface-border rounded px-3 py-1.5 text-sm text-gray-200 focus:border-primary focus:outline-none disabled:opacity-50"
            />
            {showDropdown && allGenes.length > 0 && filteredGenes.length > 0 && (
              <div className="absolute z-20 top-full mt-1 w-full bg-surface-card border border-surface-border rounded shadow-xl max-h-52 overflow-y-auto">
                {filteredGenes.map(gene => (
                  <button
                    key={gene}
                    onClick={() => addGene(gene)}
                    disabled={selectedGenes.includes(gene)}
                    className="w-full text-left px-3 py-1.5 text-sm text-gray-200 hover:bg-surface-border disabled:text-gray-600 disabled:cursor-not-allowed"
                  >
                    {gene}
                    {selectedGenes.includes(gene) && <span className="ml-2 text-green-400 text-xs">✓</span>}
                  </button>
                ))}
              </div>
            )}
          </div>
        </div>

        {/* Preset gene sets */}
        <div>
          <label className="block text-xs text-gray-400 mb-1.5">{t('spatial.presets')}</label>
          <div className="flex flex-wrap gap-1.5">
            {PRESET_GENE_SETS.map(ps => (
              <button
                key={ps.label}
                onClick={() => applyPreset(ps)}
                className="px-2.5 py-0.5 rounded-full text-xs border border-surface-border text-gray-400 hover:border-primary hover:text-primary transition-colors"
              >
                {ps.label}
              </button>
            ))}
          </div>
        </div>

        {/* Selected gene chips */}
        {selectedGenes.length > 0 && (
          <div className="flex flex-wrap gap-2">
            <span className="text-xs text-gray-500 self-center">{t('spatial.selected')}:</span>
            {selectedGenes.map(gene => (
              <span key={gene}
                className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-primary/20 border border-primary/40 text-primary"
              >
                {gene}
                <button onClick={() => removeGene(gene)} className="hover:text-red-400 leading-none">×</button>
              </span>
            ))}
          </div>
        )}

        {/* Row 2: Mode toggle + visual params */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {/* Mode toggle */}
          <div>
            <label className="block text-xs text-gray-400 mb-1">{t('spatial.mode')}</label>
            <div className="flex rounded overflow-hidden border border-surface-border">
              {modeOptions.map(opt => (
                <button
                  key={opt.value}
                  onClick={() => setMode(opt.value)}
                  className={`flex-1 py-1.5 text-xs font-medium transition-colors ${
                    mode === opt.value
                      ? 'bg-primary text-black'
                      : 'bg-surface text-gray-400 hover:text-gray-200'
                  }`}
                >
                  {opt.label}
                </button>
              ))}
            </div>
          </div>

          {/* Set name input (visible only in set mode) */}
          {mode === 'set' && (
            <div>
              <label className="block text-xs text-gray-400 mb-1">{t('spatial.set_name')}</label>
              <input
                type="text"
                value={setName}
                onChange={e => setSetName(e.target.value)}
                placeholder={t('spatial.set_name_ph')}
                className="w-full bg-surface border border-surface-border rounded px-3 py-1.5 text-sm text-gray-200 focus:border-primary focus:outline-none"
              />
            </div>
          )}
        </div>

        {/* Row 3: visual params */}
        <div className="grid grid-cols-3 gap-4">
          <div>
            <label className="block text-xs text-gray-400 mb-1">
              {t('spatial.point_size')}: {pointSize}
            </label>
            <input type="range" min={2} max={20} value={pointSize}
              onChange={e => setPointSize(Number(e.target.value))}
              className="w-full accent-primary" />
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">{t('spatial.colormap')}</label>
            <select value={cmap} onChange={e => setCmap(e.target.value)}
              className="w-full bg-surface border border-surface-border rounded px-2 py-1.5 text-sm text-gray-200 focus:border-primary focus:outline-none">
              {CMAPS.map(c => <option key={c} value={c}>{c}</option>)}
            </select>
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">
              {t('spatial.alpha')}: {alpha.toFixed(1)}
            </label>
            <input type="range" min={0.1} max={1.0} step={0.1} value={alpha}
              onChange={e => setAlpha(Number(e.target.value))}
              className="w-full accent-primary" />
          </div>
        </div>

        {/* Plot button */}
        <div className="flex items-center gap-3">
          <button
            onClick={handlePlot}
            disabled={plotting || selectedGenes.length === 0}
            className="px-5 py-2 rounded text-sm font-medium bg-primary text-black hover:bg-primary/80 disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
          >
            {plotting ? t('spatial.plotting') : t('spatial.plot')}
          </button>
          {plotInfo?.n_cells != null && !plotting && (
            <span className="text-xs text-gray-500">{plotInfo.n_cells.toLocaleString()} {t('spatial.cells')}</span>
          )}
          {plotError && (
            <span className="text-xs text-red-400">{plotError}</span>
          )}
        </div>
      </div>

      {/* Result image */}
      <div className="bg-surface-card rounded-xl border border-surface-border p-4">
        {plotImage ? (
          <>
            <div className="flex justify-end mb-2">
              <button
                onClick={handleSavePng}
                className="px-3 py-1 rounded text-xs font-medium border border-surface-border text-gray-300 hover:border-primary hover:text-primary transition-colors"
              >
                ↓ {t('spatial.save_png')}
              </button>
            </div>
            <img
              src={`data:image/png;base64,${plotImage}`}
              className="w-full rounded"
              alt="spatial gene expression"
            />
          </>
        ) : (
          <p className="text-xs text-gray-500 text-center py-12">
            {plotting ? t('spatial.plotting') : t('spatial.no_result')}
          </p>
        )}
      </div>
    </div>
  )
}

export default function Stage35_SpatialExplorer() {
  return (
    <ErrorBoundary>
      <SpatialExplorerInner />
    </ErrorBoundary>
  )
}
