import { useState, useEffect, useCallback } from 'react'
import { scanData, applyData, getDataStatus, browseDir, getOutputDir } from '../api/client'
import { FolderSearch, Check, AlertTriangle, HardDrive, FileSearch, FolderOpen, ChevronRight, ArrowUp, File, X } from 'lucide-react'
import { useT } from '../i18n'

// ── Types ──────────────────────────────────────────────────

interface DiscoveredFile {
    path: string
    label: string
    size_bytes?: number
    size_human?: string
}

interface ScanResult {
    data_root: string
    he_image: DiscoveredFile | null
    binned_002: DiscoveredFile | null
    binned_008: DiscoveredFile | null
    extra_files: { path: string; label: string; size_human: string }[]
    warnings: string[]
}

interface PathStatus {
    path: string
    configured: boolean
}

interface BrowseItem {
    name: string
    path: string
    type: 'dir' | 'file'
    children?: number
    size?: number
    size_human?: string
}

interface BrowseData {
    current: string
    parent: string | null
    items: BrowseItem[]
}

// ── Folder Browser Modal ───────────────────────────────────

function FolderBrowser({
    onSelect,
    onClose,
}: {
    onSelect: (path: string) => void
    onClose: () => void
}) {
    const t = useT()
    const [browseData, setBrowseData] = useState<BrowseData | null>(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState('')
    const [pathInput, setPathInput] = useState('')

    const navigate = useCallback(async (path: string) => {
        setLoading(true)
        setError('')
        try {
            const r = await browseDir(path)
            if (r.data.status === 'ok') {
                setBrowseData(r.data.data)
                setPathInput(r.data.data.current)
            } else {
                setError(r.data.message)
            }
        } catch (e: unknown) {
            setError(e instanceof Error ? e.message : '連線失敗')
        } finally {
            setLoading(false)
        }
    }, [])

    useEffect(() => {
        navigate('~')
    }, [navigate])

    const handlePathSubmit = () => {
        if (pathInput.trim()) navigate(pathInput.trim())
    }

    const QUICK_LOCATIONS = [
        { label: '🏠 Home', path: '~' },
        { label: '💾 /Volumes', path: '/Volumes' },
        { label: '/', path: '/' },
    ]

    return (
        <div className="fixed inset-0 bg-black/60 backdrop-blur-sm z-50 flex items-center justify-center p-6">
            <div className="bg-surface-card border border-surface-border rounded-2xl w-full max-w-2xl max-h-[80vh] flex flex-col shadow-2xl">
                {/* Header */}
                <div className="flex items-center justify-between px-5 py-3.5 border-b border-surface-border">
                    <div className="flex items-center gap-2">
                        <FolderOpen className="w-4 h-4 text-primary" />
                        <h3 className="font-semibold text-gray-200 text-sm">{t('data.browser.title')}</h3>
                    </div>
                    <button onClick={onClose} className="text-gray-500 hover:text-gray-300 transition-colors">
                        <X className="w-4 h-4" />
                    </button>
                </div>

                {/* 快捷位置 + 路徑輸入 */}
                <div className="px-5 py-2.5 border-b border-surface-border bg-surface/50 space-y-2">
                    {/* 快捷按鈕 */}
                    <div className="flex items-center gap-1.5">
                        {QUICK_LOCATIONS.map(loc => (
                            <button
                                key={loc.path}
                                onClick={() => navigate(loc.path)}
                                className="px-2.5 py-0.5 bg-surface border border-surface-border rounded text-xs text-gray-300 hover:border-primary/60 hover:text-primary transition-colors"
                            >
                                {loc.label}
                            </button>
                        ))}
                    </div>
                    {/* 路徑輸入（始終顯示） */}
                    <div className="flex items-center gap-2">
                        <input
                            value={pathInput}
                            onChange={e => setPathInput(e.target.value)}
                            onKeyDown={e => {
                                if (e.key === 'Enter') handlePathSubmit()
                            }}
                            placeholder="/Volumes/SSD/plan_a/tissue sample/CRC"
                            className="flex-1 px-2 py-1 bg-surface border border-surface-border rounded text-xs text-gray-200 font-mono focus:outline-none focus:border-primary"
                        />
                        <button
                            onClick={handlePathSubmit}
                            className="px-3 py-1 bg-primary text-white rounded text-xs hover:bg-primary-dark transition-colors flex-shrink-0"
                        >{t('data.browser.jump')}</button>
                    </div>
                </div>

                {/* Content */}

                <div className="flex-1 overflow-y-auto p-2 min-h-[300px]">
                    {loading && <div className="text-sm text-gray-500 p-4 text-center">{t('common.loading')}</div>}
                    {error && <div className="text-sm text-red-400 p-4 text-center">{error}</div>}

                    {browseData && !loading && (
                        <div className="space-y-0.5">
                            {/* Go up */}
                            {browseData.parent && (
                                <button
                                    onClick={() => navigate(browseData.parent!)}
                                    className="w-full flex items-center gap-2.5 px-3 py-2 rounded-lg hover:bg-surface-border/50 text-left transition-colors group"
                                >
                                    <ArrowUp className="w-4 h-4 text-gray-500 group-hover:text-primary" />
                                    <span className="text-sm text-gray-400 group-hover:text-gray-200">..</span>
                                </button>
                            )}

                            {browseData.items.map((item) => (
                                <button
                                    key={item.path}
                                    onClick={() => item.type === 'dir' ? navigate(item.path) : undefined}
                                    className={`w-full flex items-center gap-2.5 px-3 py-2 rounded-lg text-left transition-colors ${item.type === 'dir'
                                        ? 'hover:bg-surface-border/50 cursor-pointer group'
                                        : 'opacity-50 cursor-default'
                                        }`}
                                >
                                    {item.type === 'dir' ? (
                                        <FolderOpen className="w-4 h-4 text-yellow-500/80 flex-shrink-0" />
                                    ) : (
                                        <File className="w-4 h-4 text-gray-600 flex-shrink-0" />
                                    )}
                                    <span className={`text-sm flex-1 truncate ${item.type === 'dir' ? 'text-gray-300 group-hover:text-gray-100' : 'text-gray-500'
                                        }`}>
                                        {item.name}
                                    </span>
                                    {item.type === 'dir' && item.children != null && (
                                        <span className="text-xs text-gray-600">{item.children} {t('data.browser.items')}</span>
                                    )}
                                    {item.type === 'file' && item.size_human && (
                                        <span className="text-xs text-gray-600">{item.size_human}</span>
                                    )}
                                    {item.type === 'dir' && (
                                        <ChevronRight className="w-3.5 h-3.5 text-gray-600 group-hover:text-gray-400" />
                                    )}
                                </button>
                            ))}

                            {browseData.items.length === 0 && (
                                <p className="text-sm text-gray-600 text-center py-6">{t('data.browser.empty')}</p>
                            )}
                        </div>
                    )}
                </div>

                {/* Footer */}
                <div className="flex items-center justify-between px-5 py-3 border-t border-surface-border">
                    <p className="text-xs text-gray-500 truncate max-w-[60%]">
                        {browseData?.current ?? ''}
                    </p>
                    <button
                        onClick={() => browseData && onSelect(browseData.current)}
                        disabled={!browseData}
                        className="px-5 py-1.5 bg-primary text-white rounded-lg text-sm font-medium hover:bg-primary-dark transition-colors disabled:opacity-40"
                    >
                        {t('data.browser.select')}
                    </button>
                </div>
            </div>
        </div>
    )
}

// ── Main Page ──────────────────────────────────────────────

export default function DataSetup() {
    const t = useT()

    const FILE_LABELS: Record<string, string> = {
        he_image: t('data.file.he_image'),
        binned_002: t('data.file.binned_002'),
        binned_008: t('data.file.binned_008'),
    }

    const [dataRoot, setDataRoot] = useState('')
    const [scanning, setScanning] = useState(false)
    const [scanResult, setScanResult] = useState<ScanResult | null>(null)
    const [scanError, setScanError] = useState('')
    const [applying, setApplying] = useState(false)
    const [applied, setApplied] = useState(false)
    const [pathStatus, setPathStatus] = useState<Record<string, PathStatus>>({})
    const [showBrowser, setShowBrowser] = useState(false)

    // output dir state
    const [outputDir, setOutputDir] = useState('')
    const [showOutputBrowser, setShowOutputBrowser] = useState(false)
    const [savingOutput, setSavingOutput] = useState(false)
    const [savedOutput, setSavedOutput] = useState(false)
    const [outputError, setOutputError] = useState('')

    // 載入目前配置狀態
    useEffect(() => {
        getDataStatus().then((r: { data: { status: string; data: Record<string, PathStatus> } }) => {
            if (r.data.status === 'ok') setPathStatus(r.data.data)
        }).catch(() => { })
        getOutputDir().then((r: { data: { status: string; data: { output_dir: string; resolved: string } } }) => {
            if (r.data.status === 'ok') setOutputDir(r.data.data.resolved)
        }).catch(() => { })
    }, [applied])

    const handleScan = async () => {
        if (!dataRoot.trim()) return
        setScanning(true)
        setScanResult(null)
        setScanError('')
        setApplied(false)
        try {
            const r = await scanData({ data_root: dataRoot })
            if (r.data.status === 'ok') {
                setScanResult(r.data.data)
            } else {
                setScanError(r.data.message ?? '掃描失敗')
            }
        } catch (e: unknown) {
            setScanError(
                e instanceof Error ? e.message : t('data.error.connection')
            )
        } finally {
            setScanning(false)
        }
    }

    const handleApply = async () => {
        if (!scanResult) return
        setApplying(true)
        try {
            const paths: Record<string, string> = {}
            if (scanResult.he_image) paths.he_image = scanResult.he_image.path
            if (scanResult.binned_002) paths.binned_002 = scanResult.binned_002.path
            if (scanResult.binned_008) paths.binned_008 = scanResult.binned_008.path
            await applyData(paths)
            setApplied(true)
        } catch {
            // noop
        } finally {
            setApplying(false)
        }
    }

    const handleBrowseSelect = (path: string) => {
        setDataRoot(path)
        setShowBrowser(false)
    }

    const handleOutputBrowseSelect = (path: string) => {
        setOutputDir(path)
        setShowOutputBrowser(false)
    }

    const handleSaveOutput = async () => {
        if (!outputDir.trim()) return
        setSavingOutput(true)
        setSavedOutput(false)
        setOutputError('')
        try {
            const r = await applyData({ output_dir: outputDir.trim() })
            if (r.data.status === 'ok') {
                setSavedOutput(true)
                setTimeout(() => setSavedOutput(false), 3000)
            } else {
                setOutputError(r.data.message ?? '儲存失敗')
            }
        } catch (e: unknown) {
            setOutputError(e instanceof Error ? e.message : '無法連線至後端')
        } finally {
            setSavingOutput(false)
        }
    }

    const foundCount = scanResult
        ? [scanResult.he_image, scanResult.binned_002, scanResult.binned_008].filter(Boolean).length
        : 0

    const configuredCount = Object.values(pathStatus).filter(s => s.configured).length

    return (
        <div className="space-y-4">
            {/* Folder Browser Modal */}
            {showBrowser && (
                <FolderBrowser
                    onSelect={handleBrowseSelect}
                    onClose={() => setShowBrowser(false)}
                />
            )}
            {/* Folder Browser Modal - output dir */}
            {showOutputBrowser && (
                <FolderBrowser
                    onSelect={handleOutputBrowseSelect}
                    onClose={() => setShowOutputBrowser(false)}
                />
            )}

            {/* 目前配置狀態 */}
            <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
                <div className="flex items-center gap-2">
                    <HardDrive className="w-4 h-4 text-primary" />
                    <h3 className="font-semibold text-gray-200">{t('data.current.config')}</h3>
                    <span className="text-xs text-gray-500 ml-auto">{configuredCount}/4 {t('data.configured_count_suffix')}</span>
                </div>
                <div className="grid grid-cols-1 gap-2">
                    {Object.entries(FILE_LABELS).map(([key, label]) => {
                        const s = pathStatus[key]
                        return (
                            <div key={key} className="flex items-center gap-3 px-3 py-2 bg-surface/50 rounded-lg">
                                {s?.configured
                                    ? <Check className="w-4 h-4 text-green-400 flex-shrink-0" />
                                    : <AlertTriangle className="w-4 h-4 text-yellow-500 flex-shrink-0" />}
                                <div className="min-w-0 flex-1">
                                    <p className="text-sm text-gray-300">{label}</p>
                                    {s?.configured
                                        ? <p className="text-xs text-gray-500 truncate">{s.path}</p>
                                        : <p className="text-xs text-yellow-500/70">{t('data.current.not_configured')}</p>}
                                </div>
                            </div>
                        )
                    })}
                </div>
            </div>

            {/* 掃描區域 */}
            <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
                <div className="flex items-center gap-2">
                    <FolderSearch className="w-4 h-4 text-primary" />
                    <h3 className="font-semibold text-gray-200">{t('data.scan.title')}</h3>
                </div>
                <p className="text-xs text-gray-400">
                    {t('data.scan.description')}
                </p>
                {scanError && (
                    <div className="flex items-center gap-2 px-3 py-2 bg-red-900/20 border border-red-700/40 rounded-lg">
                        <AlertTriangle className="w-4 h-4 text-red-400 flex-shrink-0" />
                        <p className="text-xs text-red-400">{scanError}</p>
                    </div>
                )}
                <div className="flex gap-2">
                    <button
                        onClick={() => setShowBrowser(true)}
                        className="px-3 py-2 bg-surface border border-surface-border rounded-lg text-sm text-gray-300 hover:bg-surface-border hover:text-gray-100 transition-colors flex items-center gap-1.5 flex-shrink-0"
                        title="瀏覽資料夾"
                    >
                        <FolderOpen className="w-4 h-4" />
                        {t('data.scan.browse')}
                    </button>
                    <input
                        value={dataRoot}
                        onChange={e => setDataRoot(e.target.value)}
                        placeholder="/path/to/data/root"
                        onKeyDown={e => e.key === 'Enter' && handleScan()}
                        className="flex-1 px-3 py-2 bg-surface border border-surface-border rounded-lg text-sm text-gray-200 placeholder-gray-600 focus:border-primary focus:outline-none font-mono"
                    />
                    <button
                        onClick={handleScan}
                        disabled={scanning || !dataRoot.trim()}
                        className={`px-5 py-2 rounded-lg text-sm font-medium transition-colors flex items-center gap-2 flex-shrink-0 ${scanning || !dataRoot.trim()
                            ? 'bg-gray-700 text-gray-500 cursor-not-allowed'
                            : 'bg-primary text-white hover:bg-primary-dark'
                            }`}
                    >
                        <FileSearch className="w-4 h-4" />
                        {scanning ? t('data.scan.scanning') : t('data.scan.scan')}
                    </button>
                </div>
            </div>

            {/* 輸出目錄設定 */}
            <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
                <div className="flex items-center gap-2">
                    <HardDrive className="w-4 h-4 text-primary" />
                    <h3 className="font-semibold text-gray-200">{t('data.output.title')}</h3>
                </div>
                <p className="text-xs text-gray-400">
                    {t('data.output.description')} <span className="font-mono text-gray-300">roi/</span>, <span className="font-mono text-gray-300">analysis/</span>.
                </p>
                {outputError && (
                    <div className="flex items-center gap-2 px-3 py-2 bg-red-900/20 border border-red-700/40 rounded-lg">
                        <AlertTriangle className="w-4 h-4 text-red-400 flex-shrink-0" />
                        <p className="text-xs text-red-400">{outputError}</p>
                    </div>
                )}
                <div className="flex gap-2">
                    <button
                        onClick={() => setShowOutputBrowser(true)}
                        className="px-3 py-2 bg-surface border border-surface-border rounded-lg text-sm text-gray-300 hover:bg-surface-border hover:text-gray-100 transition-colors flex items-center gap-1.5 flex-shrink-0"
                    >
                        <FolderOpen className="w-4 h-4" />
                        {t('data.scan.browse')}
                    </button>
                    <input
                        value={outputDir}
                        onChange={e => { setOutputDir(e.target.value); setSavedOutput(false) }}
                        placeholder="/Volumes/SSD/plan_a/visiumHD_pipeline_3/results/analysis"
                        onKeyDown={e => e.key === 'Enter' && handleSaveOutput()}
                        className="flex-1 px-3 py-2 bg-surface border border-surface-border rounded-lg text-sm text-gray-200 placeholder-gray-600 focus:border-primary focus:outline-none font-mono"
                    />
                    <button
                        onClick={handleSaveOutput}
                        disabled={savingOutput || !outputDir.trim()}
                        className={`px-5 py-2 rounded-lg text-sm font-medium transition-colors flex items-center gap-2 flex-shrink-0 ${
                            savedOutput
                                ? 'bg-green-900/40 text-green-400'
                                : savingOutput || !outputDir.trim()
                                    ? 'bg-gray-700 text-gray-500 cursor-not-allowed'
                                    : 'bg-primary text-white hover:bg-primary-dark'
                        }`}
                    >
                        <Check className="w-4 h-4" />
                        {savedOutput ? t('common.saved') : savingOutput ? t('common.saving') : t('common.save')}
                    </button>
                </div>
            </div>

            {/* 掃描結果 */}
            {scanResult && (
                <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
                    <div className="flex items-center justify-between">
                        <h3 className="font-semibold text-gray-200">
                            掃描結果 — 找到 {foundCount}/4 項
                        </h3>
                        {foundCount > 0 && (
                            <button
                                onClick={handleApply}
                                disabled={applying || applied}
                                className={`px-4 py-1.5 rounded-lg text-sm font-medium transition-colors ${applied
                                    ? 'bg-green-900/40 text-green-400 cursor-default'
                                    : applying
                                        ? 'bg-gray-700 text-gray-500 cursor-not-allowed'
                                        : 'bg-primary text-white hover:bg-primary-dark'
                                    }`}
                            >
                                {applied ? t('data.results.applied') : applying ? t('data.results.applying') : t('data.results.apply')}
                            </button>
                        )}
                    </div>

                    <div className="space-y-2">
                        {(['he_image', 'binned_002', 'binned_008'] as const).map(key => {
                            const item = scanResult[key]
                            return (
                                <div key={key} className={`flex items-center gap-3 px-3 py-2.5 rounded-lg ${item ? 'bg-green-900/10 border border-green-700/30' : 'bg-surface/50'
                                    }`}>
                                    {item
                                        ? <Check className="w-4 h-4 text-green-400 flex-shrink-0" />
                                        : <AlertTriangle className="w-4 h-4 text-gray-600 flex-shrink-0" />}
                                    <div className="min-w-0 flex-1">
                                        <p className="text-sm text-gray-300">{FILE_LABELS[key]}</p>
                                        {item ? (
                                            <div className="flex items-center gap-2">
                                                <p className="text-xs text-gray-500 truncate font-mono">{item.path}</p>
                                                {item.size_human && (
                                                    <span className="text-xs text-gray-600 flex-shrink-0">({item.size_human})</span>
                                                )}
                                            </div>
                                        ) : (
                                            <p className="text-xs text-gray-600">{t('data.results.not_found')}</p>
                                        )}
                                    </div>
                                </div>
                            )
                        })}
                    </div>

                    {/* 警告 */}
                    {scanResult.warnings.length > 0 && (
                        <div className="text-xs text-yellow-500/70 space-y-1 pt-1">
                            {scanResult.warnings.map((w, i) => (
                                <p key={i}>⚠ {w}</p>
                            ))}
                        </div>
                    )}

                    {/* 備選檔案 */}
                    {scanResult.extra_files.length > 0 && (
                        <div className="pt-2">
                            <p className="text-xs text-gray-500 mb-1">其他發現的檔案：</p>
                            {scanResult.extra_files.map((f, i) => (
                                <p key={i} className="text-xs text-gray-600 font-mono truncate">{f.label}: {f.path}</p>
                            ))}
                        </div>
                    )}
                </div>
            )}
        </div>
    )
}
