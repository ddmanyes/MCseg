import { useState, useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'
import StageCard from '../components/shared/StageCard'
import Terminal from '../components/shared/Terminal'
import { runConditions, getConditionsStatus, getConditionsResults, getConditionsRecommend } from '../api/client'
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts'
import useStageLog from '../hooks/useStageLog'

export default function Stage2b_ConditionTest() {
  useStageLog('conditions')
  const { stages, updateStage, conditionResults, setConditionResults, recommendedCondition, setRecommendedCondition } = usePipelineStore()
  const stage = stages['conditions']
  const [maxDist, setMaxDist] = useState([20, 40])
  const [compactness, setCompactness] = useState([0.03, 0.06])
  const [dilation, setDilation] = useState([10, 20])

  useEffect(() => {
    getConditionsResults().then(r => setConditionResults(r.data.data ?? []))
    getConditionsRecommend().then(r => { if (r.data.status === 'ok') setRecommendedCondition(r.data.data) })
  }, [])

  const handleRun = async () => {
    updateStage('conditions', { status: 'running', progress: 0, message: '開始條件測試...', completed: 0, total: maxDist.length * compactness.length * dilation.length })
    await runConditions({ max_dist: maxDist, compactness, dilation, quick_mode: true })
    const poll = setInterval(async () => {
      const s = await getConditionsStatus()
      updateStage('conditions', s.data)
      if (s.data.status !== 'running') {
        clearInterval(poll)
        const r = await getConditionsResults()
        setConditionResults(r.data.data ?? [])
        const rec = await getConditionsRecommend()
        if (rec.data.status === 'ok') setRecommendedCondition(rec.data.data)
      }
    }, 3000)
  }

  const okResults = conditionResults.filter(r => r.status === 'ok')

  return (
    <div className="space-y-4">
      <StageCard title="Proseg 參數條件測試" status={stage.status} progress={stage.progress}
                 message={`${stage.completed ?? 0} / ${stage.total ?? 0} 條件完成`}
                 onRun={handleRun} runLabel="開始測試">
        <div className="grid grid-cols-3 gap-4 text-sm">
          {[
            { label: 'max_dist (µm)', val: maxDist, set: setMaxDist, options: [20, 30, 40, 50] },
            { label: 'compactness', val: compactness, set: setCompactness, options: [0.03, 0.06, 0.1] },
            { label: 'dilation (px)', val: dilation, set: setDilation, options: [10, 20, 30] },
          ].map(({ label, val, set, options }) => (
            <div key={label}>
              <p className="text-xs text-gray-400 mb-2">{label}</p>
              {options.map(opt => (
                <label key={opt} className="flex items-center gap-2 text-gray-300 mb-1">
                  <input type="checkbox" checked={val.includes(opt as never)}
                    onChange={e => set(v => e.target.checked ? [...v, opt as never] : v.filter(x => x !== opt))}
                    className="rounded" />
                  {opt}
                </label>
              ))}
            </div>
          ))}
        </div>
      </StageCard>

      {/* 推薦條件 */}
      {recommendedCondition && (
        <div className="bg-green-900/20 border border-green-700/40 rounded-xl p-4">
          <p className="text-xs text-green-400 font-medium mb-1">推薦最佳條件</p>
          <p className="text-sm text-gray-200">
            max_dist={recommendedCondition.max_dist} µm ·
            compactness={recommendedCondition.compactness} ·
            dilation={recommendedCondition.dilation} px
          </p>
          <p className="text-xs text-gray-400 mt-1">
            n_cells={recommendedCondition.n_cells} · median_genes={recommendedCondition.median_genes?.toFixed(1)}
          </p>
        </div>
      )}

      {/* 散點圖 */}
      {okResults.length > 0 && (
        <div className="bg-surface-card rounded-xl border border-surface-border p-4">
          <p className="text-sm font-medium text-gray-300 mb-3">細胞數 vs 基因豐富度</p>
          <ResponsiveContainer width="100%" height={260}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" stroke="#3a3a5c" />
              <XAxis dataKey="n_cells" name="n_cells" stroke="#9ca3af" tick={{ fontSize: 11 }} label={{ value: 'n_cells', position: 'insideBottom', offset: -5, fill: '#9ca3af', fontSize: 11 }} />
              <YAxis dataKey="median_genes" name="median_genes" stroke="#9ca3af" tick={{ fontSize: 11 }} />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} content={({ payload }) => {
                if (!payload?.length) return null
                const d = payload[0].payload
                return (
                  <div className="bg-surface-card border border-surface-border rounded p-2 text-xs text-gray-300">
                    <p className="font-medium">{d.label}</p>
                    <p>n_cells: {d.n_cells}</p>
                    <p>median_genes: {d.median_genes?.toFixed(1)}</p>
                  </div>
                )
              }} />
              <Scatter data={okResults} fill="#3b82f6" fillOpacity={0.8} />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      )}

      <Terminal stage="conditions" />
    </div>
  )
}
