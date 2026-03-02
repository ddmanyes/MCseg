import { useState, useEffect } from 'react'
import { scanData, applyData, getDataStatus } from '../api/client'
import { FolderSearch, Check, AlertTriangle, HardDrive, FileSearch } from 'lucide-react'

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

const FILE_LABELS: Record<string, string> = {
    he_image: 'H&E 影像（BTF/TIFF）',
    binned_002: 'Visium HD 2µm（square_002um）',
    binned_008: 'Visium HD 8µm（square_008um）',
    xenium_outs: 'Xenium Outs（可選）',
}

export default function DataSetup() {
    const [dataRoot, setDataRoot] = useState('')
    const [scanning, setScanning] = useState(false)
    const [scanResult, setScanResult] = useState<ScanResult | null>(null)
    const [applying, setApplying] = useState(false)
    const [applied, setApplied] = useState(false)
    const [pathStatus, setPathStatus] = useState<Record<string, PathStatus>>({})

    // 載入目前配置狀態
    useEffect(() => {
        getDataStatus().then(r => {
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
        } catch (e: any) {
            console.error(e)
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
        } catch (e: any) {
            console.error(e)
        } finally {
            setApplying(false)
        }
    }

    const foundCount = scanResult
        ? [scanResult.he_image, scanResult.binned_002, scanResult.binned_008, scanResult.xenium_outs].filter(Boolean).length
        : 0

    const configuredCount = Object.values(pathStatus).filter(s => s.configured).length

    return (
        <div className="space-y-4">
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
                    輸入資料根目錄路徑，系統將自動尋找 SpaceRanger 和 Xenium 輸出檔案。
                </p>
                <div className="flex gap-2">
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
                        className={`px-5 py-2 rounded-lg text-sm font-medium transition-colors flex items-center gap-2 ${scanning || !dataRoot.trim()
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
