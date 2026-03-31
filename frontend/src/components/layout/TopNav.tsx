import { NavLink } from 'react-router-dom'
import { clsx } from 'clsx'
import { usePipelineStore } from '../../stores/pipelineStore'
import type { StageStatus } from '../../types/pipeline'
import { Microscope, Settings, Loader2 } from 'lucide-react'

const STAGES = [
  { path: '/data',         idx: '⬡',  title: 'Setup',      stage: 'data',         dep: null },
  { path: '/roi',          idx: '0',   title: 'ROI',        stage: 'roi',          dep: null },
  { path: '/segmentation', idx: '1',   title: 'Seg',        stage: 'segmentation', dep: 'roi' },
  { path: '/count',        idx: '2',   title: 'Count',      stage: 'count',        dep: 'segmentation' },
  { path: '/analysis',     idx: '3',   title: 'Analysis',   stage: 'analysis',     dep: 'count' },
  { path: '/export',       idx: '4',   title: 'Export',     stage: 'export',       dep: 'analysis' },
] as const

function StageCircle({ status, isActive, idx }: { status: StageStatus; isActive: boolean; idx: string }) {
  return (
    <div className={clsx(
      'w-[18px] h-[18px] rounded-full flex items-center justify-center text-[9px] font-bold border flex-shrink-0 transition-colors',
      isActive
        ? 'bg-blue-500 border-blue-300 text-white'
        : status === 'done'
          ? 'bg-green-500 border-green-400 text-black'
          : status === 'running'
            ? 'bg-amber-500 border-amber-400 text-black animate-pulse'
            : status === 'error'
              ? 'bg-red-500 border-red-400 text-white'
              : 'bg-transparent border-gray-600 text-gray-500',
    )}>
      {!isActive && status === 'done' ? '✓' : !isActive && status === 'error' ? '✕' : idx}
    </div>
  )
}

export default function TopNav() {
  const stages = usePipelineStore(s => s.stages)

  const isLocked = (dep: string | null) =>
    dep !== null && stages[dep as keyof typeof stages]?.status !== 'done'

  const runningStage = STAGES.find(s => stages[s.stage]?.status === 'running')

  return (
    <header className="h-11 flex items-center flex-shrink-0 border-b border-white/[0.06]"
            style={{ background: 'linear-gradient(180deg, #12121a 0%, #0f0f17 100%)' }}>

      {/* ── Logo ───────────────────────────────────────────── */}
      <div className="flex items-center gap-2 px-4 h-full border-r border-white/[0.06] flex-shrink-0">
        <Microscope className="w-4 h-4 text-blue-400" strokeWidth={1.5} />
        <span className="text-sm font-semibold tracking-wide text-gray-100">MSseg</span>
        <span className="text-[10px] text-gray-600 font-mono leading-none px-1 py-0.5
                         bg-white/5 rounded border border-white/10">v1</span>
      </div>

      {/* ── Stage tabs ─────────────────────────────────────── */}
      <nav className="flex items-center h-full flex-1 px-2 overflow-x-auto
                      [&::-webkit-scrollbar]:hidden [-ms-overflow-style:none] [scrollbar-width:none]">
        {STAGES.map((s, i) => {
          const status  = stages[s.stage]?.status ?? 'idle'
          const locked  = isLocked(s.dep)
          const depLabel = s.dep ? STAGES.find(x => x.stage === s.dep)?.title : null

          const inner = (isActive: boolean) => (
            <div className={clsx(
              'flex items-center gap-1.5 px-3 h-7 rounded-md text-xs font-medium',
              'transition-all duration-150 select-none flex-shrink-0',
              isActive
                ? 'bg-blue-600/20 border border-blue-500/50 text-blue-200'
                : locked
                  ? 'opacity-30 cursor-not-allowed text-gray-600'
                  : status === 'done'
                    ? 'text-green-300 hover:bg-green-900/20 border border-transparent hover:border-green-700/30 cursor-pointer'
                    : status === 'running'
                      ? 'text-amber-300 border border-amber-600/40 bg-amber-900/15 cursor-pointer'
                      : status === 'error'
                        ? 'text-red-300 border border-red-600/40 bg-red-900/15 cursor-pointer'
                        : 'text-gray-400 hover:text-gray-200 hover:bg-white/[0.04] border border-transparent cursor-pointer',
            )}>
              <StageCircle status={status} isActive={isActive} idx={s.idx} />
              <span className="tracking-wide">{s.title}</span>
            </div>
          )

          const divider = i > 0 && (
            <div className="w-3 flex-shrink-0 flex items-center justify-center">
              <div className={clsx('h-px w-full transition-colors',
                stages[STAGES[i - 1].stage]?.status === 'done' ? 'bg-green-700/50' : 'bg-white/[0.06]'
              )} />
            </div>
          )

          return (
            <div key={s.path} className="flex items-center">
              {divider}
              {locked ? (
                <div title={`請先完成「${depLabel}」`}>
                  {inner(false)}
                </div>
              ) : (
                <NavLink to={s.path}>
                  {({ isActive }) => inner(isActive)}
                </NavLink>
              )}
            </div>
          )
        })}
      </nav>

      {/* ── Right: global running status ───────────────────── */}
      {runningStage && (
        <div className="flex items-center gap-1.5 px-2.5 h-6 mx-2 rounded-full
                        bg-amber-900/25 border border-amber-600/30 flex-shrink-0">
          <Loader2 className="w-3 h-3 text-amber-400 animate-spin flex-shrink-0" />
          <span className="text-[10px] text-amber-300 max-w-[160px] truncate">
            {stages[runningStage.stage]?.message || runningStage.title}
          </span>
        </div>
      )}

      {/* ── Settings icon ──────────────────────────────────── */}
      <button
        className="p-2 mx-1 mr-2 rounded-lg text-gray-600 hover:text-gray-300
                   hover:bg-white/[0.05] transition-colors flex-shrink-0"
        title="設定"
        onClick={() => {/* future: open settings panel */}}
      >
        <Settings className="w-4 h-4" strokeWidth={1.5} />
      </button>
    </header>
  )
}
