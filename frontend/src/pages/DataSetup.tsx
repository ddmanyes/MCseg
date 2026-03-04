import { useState, useEffect, useCallback } from 'react'
import { scanData, applyData, getDataStatus, browseDir } from '../api/client'
import { FolderSearch, Check, AlertTriangle, HardDrive, FileSearch, FolderOpen, ChevronRight, ArrowUp, File, X } from 'lucide-react'

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
    xenium_outs: DiscoveredFile | null
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

const FILE_LABELS: Record<string, string> = {
    he_image: 'H&E 影像（BTF/TIFF）',
    binned_002: 'Visium HD 2µm（square_002um）',
    binned_008: 'Visium HD 8µm（square_008um）',
    xenium_outs: 'Xenium Outs（可選）',
}

// ── Folder Browser Modal ───────────────────────────────────

function FolderBrowser({
    onSelect,
    onClose,
}: {
    onSelect: (path: string) => void
    onClose: () => void
}) {
    const [browseData, setBrowseData] = useState<BrowseData | null>(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState('')
    const [editingPath, setEditingPath] = useState(false)
    const [pathInput, setPathInput] = useState('')

    const navigate = useCallback(async (path: string) => {
        setLoading(true)
        setError('')
        try {
            const r = await browseDir(path)
            if (r.data.status === 'ok') {
                setBrowseData(r.data.data)
                setEditingPath(false)
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


    return (
        <div className="fixed inset-0 bg-black/60 backdrop-blur-sm z-50 flex items-center justify-center p-6">
            <div className="bg-surface-card border border-surface-border rounded-2xl w-full max-w-2xl max-h-[80vh] flex flex-col shadow-2xl">
                {/* Header */}
                <div className="flex items-center justify-between px-5 py-3.5 border-b border-surface-border">
                    <div className="flex items-center gap-2">
                        <FolderOpen className="w-4 h-4 text-primary" />
                        <h3 className="font-semibold text-gray-200 text-sm">選擇資料目錄</h3>
                    </div>
                    <button onClick={onClose} className="text-gray-500 hover:text-gray-300 transition-colors">
                        <X className="w-4 h-4" />
                    </button>
                </div>

                {/* Current path bar — 點擊可直接輸入路徑 */}
                <div className="px-5 py-2 border-b border-surface-border bg-surface/50">
                    {editingPath ? (
                        <div className="flex items-center gap-2">
                            <input
                                autoFocus
                                value={pathInput}
                                onChange={e => setPathInput(e.target.value)}
                                onKeyDown={e => {
                                    if (e.key === 'Enter') handlePathSubmit()
                                    if (e.key === 'Escape') setEditingPath(false)
                                }}
                                placeholder="/Volumes/SSD/plan_a/tissue sample/CRC"
                                className="flex-1 px-2 py-0.5 bg-surface border border-primary/50 rounded text-xs text-gray-200 font-mono focus:outline-none focus:border-primary"
                            />
                            <button
                                onClick={handlePathSubmit}
                                className="px-2 py-0.5 bg-primary text-white rounded text-xs hover:bg-primary-dark transition-colors"
                            >跳轉</button>
                            <button
                                onClick={() => setEditingPath(false)}
                                className="text-gray-500 hover:text-gray-300 text-xs"
                            >✕</button>
                        </div>
                    ) : (
                        <button
                            onClick={() => { setPathInput(browseData?.current ?? ''); setEditingPath(true) }}
                            className="w-full text-left group flex items-center gap-1"
                            title="點擊直接輸入路徑"
                        >
                            <p className="text-xs text-gray-400 font-mono truncate group-hover:text-gray-200 transition-colors">
                                {browseData?.current ?? '載入中...'}
                            </p>
                            <span className="text-xs text-gray-600 group-hover:text-primary transition-colors ml-1 flex-shrink-0">✎</span>
                        </button>
                    )}
                </div>

                {/* Content */}

                <div className="flex-1 overflow-y-auto p-2 min-h-[300px]">
                    {loading && <div className="text-sm text-gray-500 p-4 text-center">載入中...</div>}
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
                                        <span className="text-xs text-gray-600">{item.children} 項</span>
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
                                <p className="text-sm text-gray-600 text-center py-6">此目錄為空</p>
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
                        選擇此目錄
                    </button>
                </div>
            </div>
        </div>
    )
}

// ── Main Page ──────────────────────────────────────────────

export default function DataSetup() {
    const [dataRoot, setDataRoot] = useState('')
    const [scanning, setScanning] = useState(false)
    const [scanResult, setScanResult] = useState<ScanResult | null>(null)
    const [applying, setApplying] = useState(false)
    const [applied, setApplied] = useState(false)
    const [pathStatus, setPathStatus] = useState<Record<string, PathStatus>>({})
    const [showBrowser, setShowBrowser] = useState(false)

    // 載入目前配置狀態
    useEffect(() => {
        getDataStatus().then((r: { data: { status: string; data: Record<string, PathStatus> } }) => {
            if (r.data.status === 'ok') setPathStatus(r.data.data)
        }).catch(() => { })
    }, [applied])

    const handleScan = async () => {
        if (!dataRoot.trim()) return
        setScanning(true)
        setScanResult(null)
        setApplied(false)
        try {
            const r = await scanData({ data_root: dataRoot })
            if (r.data.status === 'ok') setScanResult(r.data.data)
        } catch {
            // noop
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
            if (scanResult.xenium_outs) paths.xenium_outs = scanResult.xenium_outs.path
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

    const foundCount = scanResult
        ? [scanResult.he_image, scanResult.binned_002, scanResult.binned_008, scanResult.xenium_outs].filter(Boolean).length
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

            {/* 目前配置狀態 */}
            <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
                <div className="flex items-center gap-2">
                    <HardDrive className="w-4 h-4 text-primary" />
                    <h3 className="font-semibold text-gray-200">目前資料配置</h3>
                    <span className="text-xs text-gray-500 ml-auto">{configuredCount}/4 已設定</span>
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
                                        : <p className="text-xs text-yellow-500/70">尚未設定</p>}
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
                    <h3 className="font-semibold text-gray-200">自動掃描資料目錄</h3>
                </div>
                <p className="text-xs text-gray-400">
                    選擇或輸入資料根目錄路徑，系統將自動尋找 SpaceRanger 和 Xenium 輸出檔案。
                </p>
                <div className="flex gap-2">
                    <button
                        onClick={() => setShowBrowser(true)}
                        className="px-3 py-2 bg-surface border border-surface-border rounded-lg text-sm text-gray-300 hover:bg-surface-border hover:text-gray-100 transition-colors flex items-center gap-1.5 flex-shrink-0"
                        title="瀏覽資料夾"
                    >
                        <FolderOpen className="w-4 h-4" />
                        瀏覽
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
                        {scanning ? '掃描中...' : '掃描'}
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
                                {applied ? '✓ 已套用' : applying ? '套用中...' : '套用至 Config'}
                            </button>
                        )}
                    </div>

                    <div className="space-y-2">
                        {(['he_image', 'binned_002', 'binned_008', 'xenium_outs'] as const).map(key => {
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
                                            <p className="text-xs text-gray-600">未找到</p>
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
