import { useState, useEffect, useCallback } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import {
  runConditions, getConditionsStatus, getConditionsResults,
  getConditionsRecommend, getConditionThumbnail,
} from '../api/client'
import {
  ScatterChart, Scatter, XAxis, YAxis, CartesianGrid,
  Tooltip, ResponsiveContainer,
} from 'recharts'
import useStageLog from '../hooks/useStageLog'
import { useStageStatus } from '../hooks/useStageStatus'
import type { ConditionResult } from '../types/pipeline'

// ── 排名標章 ────────────────────────────────────────────────────
const MEDALS = ['🥇', '🥈', '🥉']

function score(r: ConditionResult) {
  return r.n_cells * r.median_genes
}

function MetricBadge({ label, value, unit = '' }: { label: string; value: React.ReactNode; unit?: string }) {
  return (
    <div className="flex flex-col items-center">
      <span className="text-xs text-gray-500">{label}</span>
      <span className="text-sm font-semibold text-gray-200">{value}{unit && <span className="text-xs text-gray-400 ml-0.5">{unit}</span>}</span>
    </div>
  )
}

// ── Top-3 縮圖卡片 ───────────────────────────────────────────────
function ThumbnailCard({
  rank, result, isRecommended, thumbnail, onApply,
}: {
  rank: number
  result: ConditionResult
  isRecommended: boolean
  thumbnail: string | null
  onApply: () => void
}) {
  return (
    <div className={`rounded-xl border p-3 flex flex-col gap-2 transition-colors
      ${isRecommended
        ? 'border-green-600/50 bg-green-900/10'
        : 'border-gray-700 bg-gray-800/40'}`}
    >
      {/* 標題列 */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-1.5">
          <span className="text-lg leading-none">{MEDALS[rank]}</span>
          <span className="text-xs font-mono text-gray-300 truncate">{result.label}</span>
        </div>
        {isRecommended && (
          <span className="text-xs bg-green-800/60 text-green-300 px-1.5 py-0.5 rounded">推薦</span>
        )}
      </div>

      {/* 縮圖 */}
      <div className="relative rounded overflow-hidden bg-gray-900 aspect-[3/2]">
        {thumbnail
          ? <img src={`data:image/jpeg;base64,${thumbnail}`} alt={result.label}
            className="w-full h-full object-cover" />
          : <div className="absolute inset-0 flex items-center justify-center">
            <span className="text-xs text-gray-600 animate-pulse">載入縮圖...</span>
          </div>
        }
      </div>

      {/* 指標格 */}
      <div className="grid grid-cols-3 gap-1 py-1 border-y border-gray-700/50">
        <MetricBadge label="n_cells" value={result.n_cells.toLocaleString()} />
        <MetricBadge label="median_genes" value={result.median_genes.toFixed(1)} />
        <MetricBadge label="score" value={(score(result) / 1000).toFixed(1)} unit="k" />
      </div>

      {/* 參數列 */}
      <div className="flex flex-wrap gap-x-3 gap-y-0.5 text-xs text-gray-400">
        <span>max_dist <b className="text-gray-200">{result.max_dist}</b> µm</span>
        <span>compact <b className="text-gray-200">{result.compactness}</b></span>
        <span>dilation <b className="text-gray-200">{result.dilation}</b> px</span>
      </div>
      {typeof result.fraction_assigned === 'number' && result.fraction_assigned > 0 && (
        <div className="text-xs text-gray-500">
          RNA 指派率 <span className="text-gray-300">{(result.fraction_assigned * 100).toFixed(1)}%</span>
          {typeof result.cell_area_cv === 'number' && !isNaN(result.cell_area_cv) &&
            <> · CV <span className="text-gray-300">{result.cell_area_cv.toFixed(2)}</span></>}
        </div>
      )}

      <button
        onClick={onApply}
        className="w-full py-1.5 text-xs rounded bg-blue-700 hover:bg-blue-600 text-white transition-colors"
      >
        套用此條件 →
      </button>
    </div>
  )
}

// ── 結果表格 ────────────────────────────────────────────────────
type SortKey = 'rank' | 'n_cells' | 'median_genes' | 'median_counts' | 'fraction_assigned' | 'cell_area_cv' | 'score'

function ResultsTable({
  results, recommendedIdx,
}: { results: ConditionResult[]; recommendedIdx: number | null }) {
  const [sortKey, setSortKey] = useState<SortKey>('score')
  const [asc, setAsc] = useState(false)

  const sorted = [...results].sort((a, b) => {
    const va = sortKey === 'score' ? score(a) : (a as any)[sortKey] ?? 0
    const vb = sortKey === 'score' ? score(b) : (b as any)[sortKey] ?? 0
    return asc ? va - vb : vb - va
  })

  const Th = ({ label, k }: { label: string; k: SortKey }) => (
    <th
      className="px-3 py-2 text-left text-xs font-medium text-gray-400 cursor-pointer hover:text-gray-200 select-none whitespace-nowrap"
      onClick={() => { if (sortKey === k) setAsc(!asc); else { setSortKey(k); setAsc(false) } }}
    >
      {label}
      {sortKey === k && <span className="ml-1 text-blue-400">{asc ? '↑' : '↓'}</span>}
    </th>
  )

  const rankMap = new Map(
    [...results].sort((a, b) => score(b) - score(a)).map((r, i) => [r.condition_idx, i + 1])
  )

  return (
    <div className="overflow-x-auto rounded-xl border border-gray-700">
      <table className="w-full text-sm">
        <thead className="bg-gray-800/60">
          <tr>
            <Th label="#" k="rank" />
            <th className="px-3 py-2 text-left text-xs font-medium text-gray-400">Label</th>
            <Th label="n_cells" k="n_cells" />
            <Th label="median_genes" k="median_genes" />
            <Th label="median_counts" k="median_counts" />
            <Th label="RNA指派%" k="fraction_assigned" />
            <Th label="CV" k="cell_area_cv" />
            <Th label="Score" k="score" />
            <th className="px-3 py-2 text-xs text-gray-400">params</th>
          </tr>
        </thead>
        <tbody>
          {sorted.map(r => {
            const rank = rankMap.get(r.condition_idx) ?? 99
            const isRec = r.condition_idx === recommendedIdx
            return (
              <tr key={r.condition_idx}
                className={`border-t border-gray-700/50 transition-colors
                  ${isRec ? 'bg-green-900/10' : 'hover:bg-gray-800/40'}`}
              >
                <td className="px-3 py-2 text-gray-400 font-mono text-xs">
                  {rank <= 3 ? MEDALS[rank - 1] : rank}
                </td>
                <td className="px-3 py-2 font-mono text-xs text-gray-300 whitespace-nowrap">
                  {r.label}
                  {isRec && <span className="ml-1.5 text-xs text-green-400">★</span>}
                </td>
                <td className="px-3 py-2 text-right text-gray-200">{r.n_cells.toLocaleString()}</td>
                <td className="px-3 py-2 text-right text-gray-200">{r.median_genes.toFixed(1)}</td>
                <td className="px-3 py-2 text-right text-gray-300">{r.median_counts?.toFixed(1) ?? '—'}</td>
                <td className="px-3 py-2 text-right text-gray-300">
                  {r.fraction_assigned > 0 ? `${(r.fraction_assigned * 100).toFixed(1)}%` : '—'}
                </td>
                <td className="px-3 py-2 text-right text-gray-300">
                  {!isNaN(r.cell_area_cv) ? r.cell_area_cv.toFixed(2) : '—'}
                </td>
                <td className="px-3 py-2 text-right font-semibold text-blue-300">
                  {(score(r) / 1000).toFixed(1)}k
                </td>
                <td className="px-3 py-2 text-xs text-gray-500 whitespace-nowrap">
                  d={r.max_dist} c={r.compactness} dil={r.dilation}
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}

// ── 主頁面 ───────────────────────────────────────────────────────
export default function Stage2b_ConditionTest() {
  useStageLog('conditions')
  const { stages, updateStage, conditionResults, setConditionResults,
    recommendedCondition, setRecommendedCondition } = usePipelineStore()
  const stage = stages['conditions']
  const { refetch: refetchStatus } = useStageStatus('conditions', getConditionsStatus, 3000)

  const [maxDist, setMaxDist] = useState([20, 40])
  const [compactness, setCompactness] = useState([0.03, 0.06])
  const [dilation, setDilation] = useState([10, 20])
  const [thumbnails, setThumbnails] = useState<Record<number, string>>({})

  const okResults = conditionResults.filter(r => r.status === 'ok')
  const top3 = [...okResults].sort((a, b) => score(b) - score(a)).slice(0, 3)
  const recIdx = recommendedCondition?.condition_idx ?? null

  // 載入縮圖
  const fetchThumbnail = useCallback(async (idx: number) => {
    if (thumbnails[idx]) return
    try {
      const res = await getConditionThumbnail(idx)
      if (res.data?.status === 'ok') {
        setThumbnails(prev => ({ ...prev, [idx]: res.data.data.image_b64 }))
      }
    } catch (_) { /* 靜默失敗 */ }
  }, [thumbnails])

  useEffect(() => {
    top3.forEach(r => fetchThumbnail(r.condition_idx))
  }, [top3.map(r => r.condition_idx).join(',')])  // eslint-disable-line

  // 初始載入既有結果
  useEffect(() => {
    getConditionsResults().then(r => {
      const data = r.data.data ?? []
      setConditionResults(data)
    })
    getConditionsRecommend().then(r => {
      if (r.data.status === 'ok') setRecommendedCondition(r.data.data)
    })
  }, [])

  // 測試完成後自動刷新結果
  useEffect(() => {
    if (stage.status === 'done') {
      getConditionsResults().then(r => setConditionResults(r.data.data ?? []))
      getConditionsRecommend().then(r => {
        if (r.data.status === 'ok') setRecommendedCondition(r.data.data)
      })
    }
  }, [stage.status])

  const handleRun = async () => {
    const total = maxDist.length * compactness.length * dilation.length
    updateStage('conditions', { status: 'running', progress: 0, message: '開始條件測試...', completed: 0, total })
    await runConditions({ max_dist: maxDist, compactness, dilation, quick_mode: true })
    refetchStatus()
  }

  const applyCondition = (r: ConditionResult) => {
    setRecommendedCondition(r)
  }

  return (
    <div className="space-y-5">
      {/* ── 參數選擇卡 ── */}
      <StageCard
        title="Proseg 參數條件測試"
        status={stage.status}
        progress={stage.progress}
        message={`${stage.completed ?? 0} / ${stage.total ?? 0} 條件完成`}
        onRun={handleRun}
        runLabel="開始測試"
      >
        <div className="grid grid-cols-3 gap-6 mt-2 text-sm">
          {([
            { label: 'max_dist (µm)', val: maxDist, set: setMaxDist, options: [10, 15, 20, 30, 40, 50] },
            { label: 'compactness', val: compactness, set: setCompactness, options: [0.03, 0.06, 0.1, 0.2, 0.3] },
            { label: 'dilation (px)', val: dilation, set: setDilation, options: [5, 10, 20, 30] },
          ] as const).map(({ label, val, set, options }) => (
            <div key={label}>
              <p className="text-xs text-gray-400 mb-2 font-medium">{label}</p>
              {(options as readonly number[]).map(opt => (
                <label key={opt} className="flex items-center gap-2 text-gray-300 mb-1.5 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={(val as number[]).includes(opt)}
                    onChange={e =>
                      (set as (fn: (v: number[]) => number[]) => void)(
                        v => e.target.checked ? [...v, opt] : v.filter(x => x !== opt)
                      )
                    }
                    className="rounded accent-blue-500"
                  />
                  {opt}
                </label>
              ))}
            </div>
          ))}
        </div>
        <p className="mt-3 text-xs text-gray-500">
          共 {maxDist.length} × {compactness.length} × {dilation.length} = <b className="text-gray-400">{maxDist.length * compactness.length * dilation.length}</b> 個條件
        </p>
      </StageCard>

      {/* ── 推薦條件橫幅 ── */}
      {recommendedCondition && (
        <div className="bg-green-900/15 border border-green-700/40 rounded-xl px-4 py-3 flex items-center gap-4">
          <span className="text-2xl">🏆</span>
          <div className="flex-1">
            <p className="text-xs text-green-400 font-semibold mb-0.5">目前推薦條件（已套用至 Stage 3）</p>
            <p className="text-sm text-gray-200 font-mono">
              max_dist={recommendedCondition.max_dist} µm · compactness={recommendedCondition.compactness} · dilation={recommendedCondition.dilation} px
            </p>
            <p className="text-xs text-gray-400 mt-0.5">
              n_cells={recommendedCondition.n_cells?.toLocaleString()} · median_genes={recommendedCondition.median_genes?.toFixed(1)} · score={(score(recommendedCondition) / 1000).toFixed(1)}k
            </p>
          </div>
        </div>
      )}

      {/* ── Top 3 縮圖 ── */}
      {top3.length > 0 && (
        <div>
          <h3 className="text-sm font-semibold text-gray-300 mb-3">前三名條件 — H&amp;E 細胞輪廓疊圖</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {top3.map((r, i) => (
              <ThumbnailCard
                key={r.condition_idx}
                rank={i}
                result={r}
                isRecommended={r.condition_idx === recIdx}
                thumbnail={thumbnails[r.condition_idx] ?? null}
                onApply={() => applyCondition(r)}
              />
            ))}
          </div>
        </div>
      )}

      {/* ── 完整結果表格 ── */}
      {okResults.length > 0 && (
        <div>
          <h3 className="text-sm font-semibold text-gray-300 mb-3">完整結果（Score = n_cells × median_genes）</h3>
          <ResultsTable results={okResults} recommendedIdx={recIdx} />
        </div>
      )}

      {/* ── 散點圖 ── */}
      {okResults.length > 0 && (
        <div className="bg-gray-800/30 rounded-xl border border-gray-700 p-4">
          <p className="text-sm font-medium text-gray-300 mb-3">細胞數 vs 基因豐富度</p>
          <ResponsiveContainer width="100%" height={260}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" stroke="#374151" />
              <XAxis dataKey="n_cells" name="n_cells" stroke="#9ca3af"
                tick={{ fontSize: 11 }}
                label={{ value: 'n_cells', position: 'insideBottom', offset: -4, fill: '#9ca3af', fontSize: 11 }} />
              <YAxis dataKey="median_genes" name="median_genes" stroke="#9ca3af"
                tick={{ fontSize: 11 }}
                label={{ value: 'median_genes', angle: -90, position: 'insideLeft', fill: '#9ca3af', fontSize: 11 }} />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} content={({ payload }) => {
                if (!payload?.length) return null
                const d = payload[0].payload as ConditionResult
                return (
                  <div className="bg-gray-900 border border-gray-600 rounded p-2 text-xs text-gray-300 shadow-lg">
                    <p className="font-semibold text-gray-200 mb-1">{d.label}</p>
                    <p>n_cells: <b>{d.n_cells.toLocaleString()}</b></p>
                    <p>median_genes: <b>{d.median_genes.toFixed(1)}</b></p>
                    <p>score: <b>{(score(d) / 1000).toFixed(1)}k</b></p>
                    <p className="text-gray-500 mt-1">d={d.max_dist} c={d.compactness} dil={d.dilation}</p>
                  </div>
                )
              }} />
              <Scatter
                data={okResults}
                fill="#3b82f6"
                fillOpacity={0.75}
                stroke="#60a5fa"
                strokeWidth={0.5}
              />
              {/* 標示推薦點 */}
              {recommendedCondition && (
                <Scatter
                  data={[recommendedCondition]}
                  fill="#22c55e"
                  fillOpacity={1}
                  stroke="#4ade80"
                  strokeWidth={1}
                  r={7}
                />
              )}
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      )}

      <Terminal stage="conditions" />
    </div>
  )
}
