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

interface PresetGroup {
  group: string
  sets: PresetGeneSet[]
}

const PRESET_GROUPS: PresetGroup[] = [
  {
    group: 'Immune / Tumor',
    sets: [
      { label: 'CD4 T',       genes: ['CD4', 'IL7R', 'CCR7', 'TCF7'] },
      { label: 'CD8 T',       genes: ['CD8A', 'CD8B', 'GZMB', 'PRF1'] },
      { label: 'NK',          genes: ['NKG7', 'GNLY', 'KLRD1', 'NCR1'] },
      { label: 'Macrophage',  genes: ['CD68', 'CSF1R', 'MRC1', 'ITGAM'] },
      { label: 'B cell',      genes: ['CD19', 'MS4A1', 'CD79A', 'CD79B'] },
      { label: 'Tumor',       genes: ['EPCAM', 'KRT8', 'KRT18', 'MUC1'] },
      { label: 'Fibroblast',  genes: ['ACTA2', 'FAP', 'PDGFRA', 'COL1A1'] },
      { label: 'Endothelial', genes: ['PECAM1', 'VWF', 'CDH5', 'CLDN5'] },
    ],
  },
  {
    group: 'Hair Follicle',
    sets: [
      { label: 'Bulge (HFSC)',  genes: ['Cd34', 'Lgr5', 'Krt15', 'Sox9'] },
      { label: 'ORS',           genes: ['Krt17', 'Krt14', 'Sox9'] },
      { label: 'Matrix',        genes: ['Lef1', 'Mki67', 'Shh'] },
      { label: 'IRS',           genes: ['Krt71', 'Krt25', 'Tgm3'] },
      { label: 'Dermal Papilla',genes: ['P2ry1', 'Sox2', 'Alpl', 'Corin'] },
    ],
  },
]

// Flat list kept for backward-compat with applyPreset logic
const PRESET_GENE_SETS: PresetGeneSet[] = PRESET_GROUPS.flatMap(g => g.sets)

const CUSTOM_SETS_KEY   = 'msseg_custom_gene_sets'
const OVERRIDE_SETS_KEY = 'msseg_preset_overrides'   // { [label]: string[] }

function loadCustomSets(): PresetGeneSet[] {
  try { return JSON.parse(localStorage.getItem(CUSTOM_SETS_KEY) ?? '[]') } catch { return [] }
}
function loadOverrides(): Record<string, string[]> {
  try { return JSON.parse(localStorage.getItem(OVERRIDE_SETS_KEY) ?? '{}') } catch { return {} }
}

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

  // ── Custom & override presets ─────────────────────────────────
  const [customSets, setCustomSets]       = useState<PresetGeneSet[]>(loadCustomSets)
  const [overrides, setOverrides]         = useState<Record<string, string[]>>(loadOverrides)
  const [showNewSet, setShowNewSet]       = useState(false)
  const [newSetName, setNewSetName]       = useState('')
  const [editingLabel, setEditingLabel]   = useState<string | null>(null)  // non-null = edit mode

  // All built-in labels (for override detection)
  const builtinLabels = useMemo(() => new Set(PRESET_GENE_SETS.map(s => s.label)), [])

  // Effective genes for a built-in set (override if present)
  const effectiveGenes = (label: string, defaultGenes: string[]) =>
    overrides[label] ?? defaultGenes

  const openEditPreset = (ps: PresetGeneSet) => {
    setSelectedGenes(effectiveGenes(ps.label, ps.genes))
    setNewSetName(ps.label)
    setEditingLabel(ps.label)
    setShowNewSet(true)
  }

  const openEditCustom = (ps: PresetGeneSet) => {
    setSelectedGenes([...ps.genes])
    setNewSetName(ps.label)
    setEditingLabel(ps.label)
    setShowNewSet(true)
  }

  const saveSet = () => {
    if (!newSetName.trim() || selectedGenes.length === 0) return
    const name = newSetName.trim()

    if (editingLabel !== null && builtinLabels.has(editingLabel)) {
      // Editing a built-in set → save as override (keep same label even if renamed)
      const updatedOverrides = { ...overrides, [editingLabel]: [...selectedGenes] }
      // If user renamed it, also store under new label
      if (name !== editingLabel) updatedOverrides[name] = [...selectedGenes]
      setOverrides(updatedOverrides)
      localStorage.setItem(OVERRIDE_SETS_KEY, JSON.stringify(updatedOverrides))
    } else if (editingLabel !== null) {
      // Editing a custom set → replace in-place
      const updated = customSets.map(s =>
        s.label === editingLabel ? { label: name, genes: [...selectedGenes] } : s
      )
      setCustomSets(updated)
      localStorage.setItem(CUSTOM_SETS_KEY, JSON.stringify(updated))
    } else {
      // New custom set
      const updated = [...customSets, { label: name, genes: [...selectedGenes] }]
      setCustomSets(updated)
      localStorage.setItem(CUSTOM_SETS_KEY, JSON.stringify(updated))
    }

    setNewSetName('')
    setEditingLabel(null)
    setShowNewSet(false)
  }

  const cancelEdit = () => { setShowNewSet(false); setNewSetName(''); setEditingLabel(null) }

  const resetOverride = (label: string) => {
    const updated = { ...overrides }
    delete updated[label]
    setOverrides(updated)
    localStorage.setItem(OVERRIDE_SETS_KEY, JSON.stringify(updated))
  }

  const deleteCustomSet = (label: string) => {
    const updated = customSets.filter(s => s.label !== label)
    setCustomSets(updated)
    localStorage.setItem(CUSTOM_SETS_KEY, JSON.stringify(updated))
  }

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
              placeholder={geneLoading ? t('spatial.loading_genes') : selectedGenes.length >= MAX_GENES ? `已達上限 (${MAX_GENES})，請先移除基因` : t('spatial.gene_placeholder')}
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
        <div className="space-y-2">
          <label className="block text-xs text-gray-400">{t('spatial.presets')}</label>
          {PRESET_GROUPS.map(group => (
            <div key={group.group} className="flex flex-wrap items-center gap-1.5">
              <span className="text-xs text-gray-600 w-28 shrink-0">{group.group}</span>
              {group.sets.map(ps => {
                const isOverridden = !!overrides[ps.label]
                return (
                  <span key={ps.label} className="inline-flex items-center gap-0 group/ps">
                    <button
                      onClick={() => applyPreset({ ...ps, genes: effectiveGenes(ps.label, ps.genes) })}
                      className={`px-2.5 py-0.5 rounded-l-full text-xs border transition-colors ${
                        isOverridden
                          ? 'border-amber-600/50 text-amber-400 hover:border-primary hover:text-primary'
                          : 'border-surface-border text-gray-400 hover:border-primary hover:text-primary'
                      }`}
                      title={isOverridden ? `已修改：${overrides[ps.label]?.join(', ')}` : ps.genes.join(', ')}
                    >
                      {ps.label}{isOverridden && <span className="ml-1 text-amber-500">·</span>}
                    </button>
                    <button
                      onClick={() => openEditPreset(ps)}
                      className="px-1 py-0.5 text-xs border border-l-0 border-surface-border text-gray-700 hover:border-primary hover:text-primary opacity-0 group-hover/ps:opacity-100 transition-all"
                      title="編輯基因"
                    >✎</button>
                    {isOverridden && (
                      <button
                        onClick={() => resetOverride(ps.label)}
                        className="px-1 py-0.5 rounded-r-full text-xs border border-l-0 border-amber-600/40 text-amber-700 hover:text-red-400 hover:border-red-500 transition-colors"
                        title="恢復預設"
                      >↺</button>
                    )}
                    {!isOverridden && <span className="rounded-r-full border border-l-0 border-surface-border opacity-0 group-hover/ps:opacity-0" />}
                  </span>
                )
              })}
            </div>
          ))}

          {/* Custom presets */}
          {(customSets.length > 0 || showNewSet) && (
            <div className="flex flex-wrap items-center gap-1.5">
              <span className="text-xs text-gray-600 w-28 shrink-0">Custom</span>
              {customSets.map(ps => (
                <span key={ps.label} className="inline-flex items-center gap-0.5">
                  <button
                    onClick={() => applyPreset(ps)}
                    className="px-2.5 py-0.5 rounded-l-full text-xs border border-surface-border text-gray-400 hover:border-primary hover:text-primary transition-colors"
                  >
                    {ps.label}
                  </button>
                  <button
                    onClick={() => openEditCustom(ps)}
                    className="px-1 py-0.5 text-xs border border-l-0 border-surface-border text-gray-600 hover:border-primary hover:text-primary transition-colors"
                    title="編輯"
                  >✎</button>
                  <button
                    onClick={() => deleteCustomSet(ps.label)}
                    className="px-1.5 py-0.5 rounded-r-full text-xs border border-l-0 border-surface-border text-gray-600 hover:border-red-500 hover:text-red-400 transition-colors"
                    title="刪除"
                  >×</button>
                </span>
              ))}
            </div>
          )}

          {/* New / edit set form */}
          {showNewSet ? (
            <div className="flex items-center gap-2 mt-1">
              <span className="text-xs text-gray-600 w-28 shrink-0">
                {editingLabel ? '編輯 Set' : '新增 Set'}
              </span>
              <input
                type="text"
                value={newSetName}
                onChange={e => setNewSetName(e.target.value)}
                onKeyDown={e => { if (e.key === 'Enter') saveSet(); if (e.key === 'Escape') cancelEdit() }}
                placeholder="Set 名稱"
                autoFocus
                className="w-32 px-2 py-0.5 bg-surface border border-surface-border rounded text-xs text-gray-200 focus:border-primary focus:outline-none"
              />
              <span className="text-xs text-gray-500">← {selectedGenes.length > 0 ? selectedGenes.join(', ') : '請先選基因'}</span>
              <button
                onClick={saveSet}
                disabled={!newSetName.trim() || selectedGenes.length === 0}
                className="px-2.5 py-0.5 rounded-full text-xs border border-primary text-primary disabled:opacity-40 disabled:cursor-not-allowed"
              >儲存</button>
              <button onClick={cancelEdit} className="text-xs text-gray-600 hover:text-gray-400">取消</button>
            </div>
          ) : (
            <div className="flex items-center gap-1.5">
              <span className="w-28 shrink-0" />
              <button
                onClick={() => { setEditingLabel(null); setNewSetName(''); setShowNewSet(true) }}
                className="px-2.5 py-0.5 rounded-full text-xs border border-dashed border-surface-border text-gray-600 hover:border-primary hover:text-primary transition-colors"
              >+ 新增 Set</button>
            </div>
          )}
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
