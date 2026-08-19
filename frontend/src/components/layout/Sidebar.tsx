import { NavLink } from 'react-router-dom'
import { clsx } from 'clsx'
import { usePipelineStore } from '../../stores/pipelineStore'
import { useT } from '../../i18n'
import type { StageStatus } from '../../types/pipeline'

const STAGES = [
  { path: '/data',         label: '📂',      tKey: 'sidebar.data',     stage: 'data',         dep: null },
  { path: '/roi',          label: 'Stage 0', tKey: 'sidebar.roi',      stage: 'roi',          dep: null },
  { path: '/segmentation', label: 'Stage 1', tKey: 'sidebar.seg',      stage: 'segmentation', dep: 'roi' },
  { path: '/count',        label: 'Stage 2', tKey: 'sidebar.count',    stage: 'count',        dep: 'segmentation' },
  { path: '/analysis',     label: 'Stage 3', tKey: 'sidebar.analysis', stage: 'analysis',     dep: 'count' },
  { path: '/export',       label: 'Stage 4', tKey: 'sidebar.export',   stage: 'export',       dep: 'analysis' },
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
  const t = useT()

  const isLocked = (dep: string | null) => {
    if (!dep) return false
    return stages[dep]?.status !== 'done'
  }

  return (
    <aside className="w-52 bg-surface-card border-r border-surface-border flex flex-col py-4">
      <div className="px-4 mb-6">
        <h1 className="text-sm font-bold text-primary leading-tight">MCseg</h1>
        <p className="text-xs text-gray-400">MCseg v2</p>
      </div>
      <nav className="flex-1 space-y-1 px-2">
        {STAGES.map(({ path, label, tKey, stage, dep }) => {
          const locked = isLocked(dep)
          const depItem = dep ? STAGES.find(s => s.stage === dep) : null
          const depLabel = depItem ? t(depItem.tKey) : ''

          if (locked) {
            return (
              <div
                key={path}
                className="flex items-center gap-3 px-3 py-2.5 rounded-lg text-sm
                           opacity-40 cursor-not-allowed select-none"
                title={t('sidebar.complete_dep', { dep: depLabel })}
              >
                <StatusDot status={stages[stage]?.status ?? 'idle'} />
                <div className="flex-1">
                  <div className="font-mono text-xs text-gray-500">{label}</div>
                  <div className="leading-tight text-gray-400">{t(tKey)}</div>
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
                clsx(
                  'flex items-center gap-3 px-3 py-2.5 rounded-lg text-sm transition-colors',
                  {
                    'bg-primary/10 text-primary font-medium border border-primary/20': isActive,
                    'text-gray-300 hover:bg-surface-border hover:text-white': !isActive,
                  }
                )
              }
            >
              <StatusDot status={stages[stage]?.status ?? 'idle'} />
              <div className="flex-1">
                <div className="font-mono text-xs text-gray-500">{label}</div>
                <div className="leading-tight">{t(tKey)}</div>
              </div>
            </NavLink>
          )
        })}
      </nav>
    </aside>
  )
}
