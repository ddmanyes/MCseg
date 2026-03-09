import { NavLink } from 'react-router-dom'
import { clsx } from 'clsx'
import { usePipelineStore } from '../../stores/pipelineStore'
import type { StageStatus } from '../../types/pipeline'

const STAGES = [
  { path: '/data',         label: '📂',      sub: '資料設定',    stage: 'data',         dep: null },
  { path: '/roi',          label: 'Stage 0', sub: 'ROI 裁切',    stage: 'roi',          dep: null },
  { path: '/segmentation', label: 'Stage 1', sub: '細胞分割',    stage: 'segmentation', dep: 'roi' },
  { path: '/count',        label: 'Stage 2', sub: 'RNA 計數',    stage: 'count',        dep: 'segmentation' },
  { path: '/analysis',     label: 'Stage 3', sub: '下游分析',    stage: 'analysis',     dep: 'count' },
  { path: '/export',       label: 'Stage 4', sub: 'Browser 匯出',stage: 'export',       dep: 'analysis' },
]

function StatusDot({ status }: { status: StageStatus }) {
  return (
    <span className={clsx('w-2 h-2 rounded-full flex-shrink-0', {
      'bg-gray-500': status === 'idle',
      'bg-yellow-400 animate-pulse': status === 'running',
      'bg-green-400': status === 'done',
      'bg-red-400': status === 'error',
    })} />
  )
}

export default function Sidebar() {
  const stages = usePipelineStore((s) => s.stages)

  const isLocked = (dep: string | null) => {
    if (!dep) return false
    return stages[dep]?.status !== 'done'
  }

  return (
    <aside className="w-52 bg-surface-card border-r border-surface-border flex flex-col py-4">
      <div className="px-4 mb-6">
        <h1 className="text-sm font-bold text-primary leading-tight">VisiumHD</h1>
        <p className="text-xs text-gray-400">Pipeline v3</p>
      </div>
      <nav className="flex-1 space-y-1 px-2">
        {STAGES.map(({ path, label, sub, stage, dep }) => {
          const locked = isLocked(dep)
          const depLabel = dep ? STAGES.find(s => s.stage === dep)?.sub : null

          if (locked) {
            return (
              <div
                key={path}
                className="flex items-center gap-3 px-3 py-2.5 rounded-lg text-sm
                           opacity-40 cursor-not-allowed select-none"
                title={`請先完成「${depLabel}」`}
              >
                <StatusDot status={stages[stage]?.status ?? 'idle'} />
                <div className="flex-1">
                  <div className="font-mono text-xs text-gray-500">{label}</div>
                  <div className="leading-tight text-gray-400">{sub}</div>
                </div>
                <span className="text-[10px] text-gray-600">🔒</span>
              </div>
            )
          }

          return (
            <NavLink
              key={path}
              to={path}
              className={({ isActive }) =>
                clsx('flex items-center gap-3 px-3 py-2.5 rounded-lg text-sm transition-colors', {
                  'bg-primary/20 text-primary font-medium': isActive,
                  'text-gray-400 hover:bg-surface-border hover:text-gray-200': !isActive,
                })
              }
            >
              <StatusDot status={stages[stage]?.status ?? 'idle'} />
              <div>
                <div className="font-mono text-xs text-gray-500">{label}</div>
                <div className="leading-tight">{sub}</div>
              </div>
            </NavLink>
          )
        })}
      </nav>
    </aside>
  )
}
