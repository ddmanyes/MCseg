import { useEffect, useRef } from 'react'
import { Terminal as XTerm } from '@xterm/xterm'
import { FitAddon } from '@xterm/addon-fit'
import '@xterm/xterm/css/xterm.css'
import { usePipelineStore } from '../../stores/pipelineStore'

interface TerminalProps {
  stage: string
  height?: string
}

export default function Terminal({ stage, height = '12rem' }: TerminalProps) {
  const logs = usePipelineStore((s) => s.stages[stage]?.logs ?? [])
  const containerRef = useRef<HTMLDivElement>(null)
  const termRef = useRef<XTerm | null>(null)
  const writtenCount = useRef(0)

  // 初始化 xterm（mount 時執行一次）
  useEffect(() => {
    if (!containerRef.current) return

    const term = new XTerm({
      theme: {
        background: '#0d1117',
        foreground: '#98c379',
        cursor: '#98c379',
        cursorAccent: '#0d1117',
        selectionBackground: '#3e4451',
      },
      fontSize: 12,
      fontFamily: 'ui-monospace, "Cascadia Code", "SF Mono", monospace',
      scrollback: 1000,
      disableStdin: true,
      convertEol: true,
    })

    const fit = new FitAddon()
    term.loadAddon(fit)
    term.open(containerRef.current)
    fit.fit()

    termRef.current = term
    writtenCount.current = 0

    // 視窗 resize 時自動調整
    const onResize = () => fit.fit()
    window.addEventListener('resize', onResize)

    return () => {
      window.removeEventListener('resize', onResize)
      term.dispose()
      termRef.current = null
    }
  }, [])

  // 增量寫入新增的 log 行（避免重繪全部）
  useEffect(() => {
    const term = termRef.current
    if (!term) return

    const newLines = logs.slice(writtenCount.current)
    for (const line of newLines) {
      if (line.includes('[ERROR]') || line.includes('ERROR'))
        term.writeln('\x1b[31m' + line + '\x1b[0m')       // 紅色
      else if (line.includes('[WARNING]') || line.includes('WARNING'))
        term.writeln('\x1b[33m' + line + '\x1b[0m')       // 黃色
      else if (line.includes('[DEBUG]') || line.includes('DEBUG'))
        term.writeln('\x1b[90m'  + line + '\x1b[0m')       // 暗灰
      else
        term.writeln('\x1b[32m'  + line + '\x1b[0m')       // 綠色（INFO）
    }
    writtenCount.current = logs.length
  }, [logs])

  return (
    <div
      ref={containerRef}
      className="rounded-lg overflow-hidden border border-surface-border"
      style={{ height }}
    />
  )
}
