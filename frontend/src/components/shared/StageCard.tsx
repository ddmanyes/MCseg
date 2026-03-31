import { clsx } from 'clsx'
import type { StageStatus } from '../../types/pipeline'
import { useT } from '../../i18n'

interface StageCardProps {
  title: string
  status: StageStatus
  progress?: number
  message?: string
  onRun?: () => void
  children?: React.ReactNode
  runLabel?: string
  disabled?: boolean
}

export default function StageCard({
  title, status, progress = 0, message = '', onRun, children,
  runLabel, disabled = false,
}: StageCardProps) {
  const t = useT()
  const isRunning = status === 'running'
  const label = runLabel ?? t('common.run')

  return (
    <div className="bg-surface-card rounded-xl border border-surface-border p-5 space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="font-semibold text-gray-200">{title}</h3>
        {onRun && (
          <button
            onClick={onRun}
            disabled={isRunning || disabled}
            className={clsx(
              'px-4 py-1.5 rounded-lg text-sm font-medium transition-colors',
              isRunning || disabled
                ? 'bg-gray-700 text-gray-500 cursor-not-allowed'
                : 'bg-primary text-white hover:bg-primary-dark'
            )}
          >
            {isRunning ? t('common.running') : label}
          </button>
        )}
      </div>

      {/* Progress bar */}
      {isRunning && (
        <div>
          <div className="flex justify-between text-xs text-gray-400 mb-1">
            <span>{message}</span>
            <span>{Math.round(progress * 100)}%</span>
          </div>
          <div className="h-1.5 bg-surface-border rounded-full overflow-hidden">
            <div
              className="h-full bg-primary rounded-full transition-all duration-300"
              style={{ width: `${progress * 100}%` }}
            />
          </div>
        </div>
      )}

      {/* Status badge */}
      {status !== 'idle' && status !== 'running' && (
        <div className={clsx('text-xs px-2 py-1 rounded inline-block', {
          'bg-green-900/40 text-green-400': status === 'done',
          'bg-red-900/40 text-red-400': status === 'error',
        })}>
          {status === 'done'
            ? `${t('common.done')} — ${message}`
            : `${t('common.error')}：${message}`}
        </div>
      )}

      {children}
    </div>
  )
}
