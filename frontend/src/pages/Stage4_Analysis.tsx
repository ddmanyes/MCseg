import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runAnalysis, getAnalysisStatus, getUmap, getConfig } from '../api/client'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'

export default function Stage4_Analysis() {
  useStageLog('analysis')
  const { stages, updateStage } = usePipelineStore()
  const stage = stages['analysis']
  const { refetch: refetchStatus } = useStageStatus('analysis', getAnalysisStatus, 5000)
  const [umapImg, setUmapImg] = useState<string | null>(null)

  // Parameters state
  const [params, setParams] = useState({
    min_genes: 20,
    min_counts: 100,
    max_genes: 8000,
    min_cells: 3,
    max_pct_mito: 20,
    resolution: 0.5,
    n_pcs: 30,
    n_neighbors: 15,
    min_dist: 0.3
  })

  // Load backend config
  useEffect(() => {
    getConfig().then(res => {
      if (res.data?.data?.analysis) {
        const ana = res.data.data.analysis
        setParams({
          min_genes: ana.preprocessing?.cellular?.min_genes ?? 20,
          min_counts: ana.preprocessing?.cellular?.min_counts ?? 100,
          max_genes: ana.preprocessing?.cellular?.max_genes ?? 8000,
          min_cells: ana.preprocessing?.cellular?.min_cells ?? 3,
          max_pct_mito: ana.preprocessing?.cellular?.max_pct_mito ?? 20,
          resolution: ana.clustering?.resolution ?? 0.5,
          n_pcs: ana.clustering?.n_pcs ?? 30,
          n_neighbors: ana.clustering?.n_neighbors ?? 15,
          min_dist: ana.clustering?.min_dist ?? 0.3
        })
      }
    })
  }, [])

  // 分析完成後自動載入 UMAP
  useEffect(() => {
    if (stage.status === 'done') {
      getUmap().then(r => { if (r.data.data) setUmapImg(r.data.data.image_b64) })
    }
  }, [stage.status])

  const handleRun = async () => {
    updateStage('analysis', { status: 'running', progress: 0, message: '執行聚類...' })
    await runAnalysis(params)
    refetchStatus()
  }

  return (
    <div className="space-y-4">
      <StageCard title="下游聚類分析（Scanpy + Leiden）" status={stage.status}
        progress={stage.progress} message={stage.message} onRun={handleRun} runLabel="執行分析">
        <div className="text-sm text-gray-400 space-y-1 mb-4">
          <p>流程：QC → normalize → HVG → PCA → UMAP → Leiden</p>
        </div>

        {/* 參數設定區塊 */}
        <div className="grid grid-cols-2 md:grid-cols-5 gap-4 mt-4 bg-surface-darker p-4 rounded-lg border border-surface-border">
          <div>
            <label className="block text-xs text-gray-400 mb-1">聚類解析度 (Resolution)</label>
            <input
              type="number" step="0.1"
              className="w-full bg-surface-highlight border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary border"
              value={params.resolution}
              onChange={e => setParams({ ...params, resolution: parseFloat(e.target.value) || 0 })}
            />
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">PCA 維度數 (n_pcs)</label>
            <input
              type="number"
              className="w-full bg-surface-highlight border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary border"
              value={params.n_pcs}
              onChange={e => setParams({ ...params, n_pcs: parseInt(e.target.value) || 0 })}
            />
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">最低 UMI 數 (min_counts)</label>
            <input
              type="number"
              className="w-full bg-surface-highlight border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary border"
              value={params.min_counts}
              onChange={e => setParams({ ...params, min_counts: parseInt(e.target.value) || 0 })}
            />
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">最低基因數 (min_genes)</label>
            <input
              type="number"
              className="w-full bg-surface-highlight border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary border"
              value={params.min_genes}
              onChange={e => setParams({ ...params, min_genes: parseInt(e.target.value) || 0 })}
            />
          </div>
          <div>
            <label className="block text-xs text-gray-400 mb-1">粒線體上限 (%)</label>
            <input
              type="number"
              className="w-full bg-surface-highlight border-gray-600 rounded px-2 py-1 text-sm text-gray-200 focus:outline-none focus:border-brand-primary border"
              value={params.max_pct_mito}
              onChange={e => setParams({ ...params, max_pct_mito: parseFloat(e.target.value) || 0 })}
            />
          </div>
        </div>
      </StageCard>

      {umapImg && (
        <div className="bg-surface-card rounded-xl border border-surface-border p-4">
          <p className="text-sm font-medium text-gray-300 mb-2">UMAP</p>
          <img src={`data:image/png;base64,${umapImg}`} className="max-w-full rounded" alt="UMAP" />
        </div>
      )}

      <Terminal stage="analysis" />
    </div>
  )
}
