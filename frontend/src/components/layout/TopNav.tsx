import { useState, useRef, useEffect } from 'react'
import { NavLink } from 'react-router-dom'
import { clsx } from 'clsx'
import { usePipelineStore } from '../../stores/pipelineStore'
import { useLanguageStore } from '../../stores/languageStore'
import { useT } from '../../i18n'
import type { StageStatus } from '../../types/pipeline'
import { Microscope, Settings, Loader2, Languages } from 'lucide-react'

const STAGE_KEYS = [
  { path: '/data',         idx: '⬡',  tKey: 'nav.stage.setup',    stage: 'data',         dep: null },
  { path: '/roi',          idx: '0',   tKey: 'nav.stage.roi',      stage: 'roi',          dep: null },
  { path: '/segmentation', idx: '1',   tKey: 'nav.stage.seg',      stage: 'segmentation', dep: 'roi' },
  { path: '/count',        idx: '2',   tKey: 'nav.stage.count',    stage: 'count',        dep: 'segmentation' },
  { path: '/analysis',     idx: '3',   tKey: 'nav.stage.analysis', stage: 'analysis',     dep: 'count' },
  { path: '/spatial',      idx: '✦',   tKey: 'nav.stage.spatial',  stage: 'spatial',      dep: 'analysis' },
  { path: '/export',       idx: '4',   tKey: 'nav.stage.export',   stage: 'export',       dep: 'analysis' },
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
  const { lang, setLang } = useLanguageStore()
  const t = useT()
  const [showLangMenu, setShowLangMenu] = useState(false)
  const menuRef = useRef<HTMLDivElement>(null)

  // Close dropdown when clicking outside
  useEffect(() => {
    const handler = (e: MouseEvent) => {
      if (menuRef.current && !menuRef.current.contains(e.target as Node)) {
        setShowLangMenu(false)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [])

  const isLocked = (dep: string | null) =>
    dep !== null && stages[dep as keyof typeof stages]?.status !== 'done'

  const runningStage = STAGE_KEYS.find(s => stages[s.stage]?.status === 'running')

  return (
    <header className="h-11 flex items-center flex-shrink-0 border-b border-white/[0.06]"
            style={{ background: 'linear-gradient(180deg, #12121a 0%, #0f0f17 100%)' }}>

      {/* ── Logo ───────────────────────────────────────────── */}
      <div className="flex items-center gap-2 px-4 h-full border-r border-white/[0.06] flex-shrink-0">
        <Microscope className="w-4 h-4 text-blue-400" strokeWidth={1.5} />
        <span className="text-sm font-semibold tracking-wide text-gray-100">MCseg</span>
        <span className="text-[10px] text-gray-600 font-mono leading-none px-1 py-0.5
                         bg-white/5 rounded border border-white/10">dev</span>
      </div>

      {/* ── Stage tabs ─────────────────────────────────────── */}
      <nav className="flex items-center h-full flex-1 px-2 overflow-x-auto
                      [&::-webkit-scrollbar]:hidden [-ms-overflow-style:none] [scrollbar-width:none]">
        {STAGE_KEYS.map((s, i) => {
          const status  = stages[s.stage]?.status ?? 'idle'
          const locked  = isLocked(s.dep)
          const depLabel = s.dep ? STAGE_KEYS.find(x => x.stage === s.dep)?.tKey : null
          const title = t(s.tKey)

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
              <span className="tracking-wide">{title}</span>
            </div>
          )

          const divider = i > 0 && (
            <div className="w-3 flex-shrink-0 flex items-center justify-center">
              <div className={clsx('h-px w-full transition-colors',
                stages[STAGE_KEYS[i - 1].stage]?.status === 'done' ? 'bg-green-700/50' : 'bg-white/[0.06]'
              )} />
            </div>
          )

          return (
            <div key={s.path} className="flex items-center">
              {divider}
              {locked ? (
                <div title={depLabel ? `${lang === 'zh' ? '請先完成「' : 'Complete '}${t(depLabel)}${lang === 'zh' ? '」' : ' first'}` : ''}>
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
            {stages[runningStage.stage]?.message || t(runningStage.tKey)}
          </span>
        </div>
      )}

      {/* ── Language + Settings ─────────────────────────────── */}
      <div className="relative flex-shrink-0 mr-1" ref={menuRef}>
        <button
          onClick={() => setShowLangMenu(v => !v)}
          className={clsx(
            'flex items-center gap-1.5 px-2 py-1.5 mx-1 rounded-lg text-xs transition-colors',
            showLangMenu
              ? 'bg-white/[0.08] text-gray-200'
              : 'text-gray-500 hover:text-gray-300 hover:bg-white/[0.05]'
          )}
          title={t('nav.settings')}
        >
          <Languages className="w-4 h-4" strokeWidth={1.5} />
          <span className="font-medium">{lang === 'zh' ? '中' : 'EN'}</span>
          <Settings className="w-3.5 h-3.5 opacity-60" strokeWidth={1.5} />
        </button>

        {showLangMenu && (
          <div className="absolute right-0 top-full mt-1 w-36 rounded-xl border border-white/[0.08]
                          bg-gray-900 shadow-2xl shadow-black/60 overflow-hidden z-50">
            <div className="px-3 pt-2.5 pb-1.5">
              <p className="text-[10px] font-semibold text-gray-500 uppercase tracking-wider">
                {lang === 'zh' ? '介面語言' : 'Interface Language'}
              </p>
            </div>
            {(['zh', 'en'] as const).map(l => (
              <button
                key={l}
                onClick={() => { setLang(l); setShowLangMenu(false) }}
                className={clsx(
                  'w-full flex items-center justify-between px-3 py-2 text-sm transition-colors',
                  lang === l
                    ? 'bg-blue-600/20 text-blue-300'
                    : 'text-gray-300 hover:bg-white/[0.05]'
                )}
              >
                <span>{l === 'zh' ? '中文' : 'English'}</span>
                {lang === l && <span className="text-blue-400 text-xs">✓</span>}
              </button>
            ))}
          </div>
        )}
      </div>
    </header>
  )
}
