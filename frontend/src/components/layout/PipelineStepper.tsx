import { useNavigate, useLocation } from 'react-router-dom'
import { clsx } from 'clsx'
import { usePipelineStore } from '../../stores/pipelineStore'

// Stage 定義：path, label, stage key (對應 Zustand), 依賴的前置 stage key
const STEPS = [
  { path: '/data',        label: 'Setup',   stage: 'data',         dep: null },
  { path: '/roi',         label: 'ROI',     stage: 'roi',          dep: null },
  { path: '/segmentation',label: 'Seg',     stage: 'segmentation', dep: 'roi' },
  { path: '/zarr',        label: 'Zarr',    stage: 'zarr',         dep: 'segmentation' },
  { path: '/conditions',  label: '2.5',     stage: 'conditions',   dep: 'zarr' },
  { path: '/proseg',      label: 'Proseg',  stage: 'proseg',       dep: 'zarr' },
  { path: '/analysis',    label: 'Analysis',stage: 'analysis',     dep: 'proseg' },
  { path: '/export',      label: 'Export',  stage: 'export',       dep: 'analysis' },
]

export default function PipelineStepper() {
  const navigate = useNavigate()
  const { pathname } = useLocation()
  const stages = usePipelineStore((s) => s.stages)

  // 判斷是否鎖定：若有依賴且依賴未完成則鎖定
  const isLocked = (dep: string | null) => {
    if (!dep) return false
    return stages[dep]?.status !== 'done'
  }

  return (
    <div className="flex items-center gap-0 px-4 py-2 bg-surface-card border-b border-surface-border overflow-x-auto">
      {STEPS.map((step, idx) => {
        const status = stages[step.stage]?.status ?? 'idle'
        const active  = pathname === step.path
        const locked  = isLocked(step.dep)
        const done    = status === 'done'
        const error   = status === 'error'
        const running = status === 'running'

        return (
          <div key={step.path} className="flex items-center flex-shrink-0">
            {/* 連接線（第一個不顯示） */}
            {idx > 0 && (
              <div className={clsx('h-px w-6 mx-1 flex-shrink-0', {
                'bg-green-500': done,
                'bg-gray-600': !done,
              })} />
            )}

            {/* 步驟按鈕 */}
            <button
              onClick={() => !locked && navigate(step.path)}
              disabled={locked}
              title={locked ? `請先完成 ${STEPS.find(s => s.stage === step.dep)?.label ?? '前一步驟'}` : step.label}
              className={clsx(
                'flex flex-col items-center gap-0.5 px-2 py-1 rounded transition-colors select-none',
                {
                  'cursor-not-allowed opacity-35': locked,
                  'cursor-pointer hover:bg-surface-border': !locked,
                }
              )}
            >
              {/* 圓圈圖示 */}
              <div className={clsx(
                'w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold border-2 transition-colors',
                {
                  'bg-green-500 border-green-400 text-black':  done && !active,
                  'bg-red-500   border-red-400   text-white':  error && !active,
                  'bg-yellow-500 border-yellow-400 text-black animate-pulse': running && !active,
                  'bg-primary   border-primary   text-black ring-2 ring-primary/30': active,
                  'bg-gray-700  border-gray-600  text-gray-400': !done && !error && !running && !active && !locked,
                  'bg-gray-800  border-gray-700  text-gray-600': locked,
                }
              )}>
                {done && !active ? '✓' : error && !active ? '✕' : idx}
              </div>

              {/* 標籤 */}
              <span className={clsx('text-[10px] leading-tight whitespace-nowrap', {
                'text-primary font-semibold': active,
                'text-green-400':  done && !active,
                'text-red-400':    error && !active,
                'text-yellow-400': running && !active,
                'text-gray-500':   !done && !error && !running && !active && !locked,
                'text-gray-600':   locked,
              })}>
                {step.label}
              </span>
            </button>
          </div>
        )
      })}
    </div>
  )
}
