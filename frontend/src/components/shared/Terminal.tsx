import { useEffect, useRef } from 'react'
import { usePipelineStore } from '../../stores/pipelineStore'

interface TerminalProps {
  stage: string
  maxLines?: number
}

export default function Terminal({ stage, maxLines = 200 }: TerminalProps) {
  const logs = usePipelineStore((s) => s.stages[stage]?.logs ?? [])
  const bottomRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' })
  }, [logs])

  return (
    <div className="bg-black/60 rounded-lg border border-surface-border h-48 overflow-y-auto p-3 font-mono text-xs text-green-400">
      {logs.slice(-maxLines).map((line, i) => (
        <div key={i} className="whitespace-pre-wrap leading-relaxed">{line}</div>
      ))}
      <div ref={bottomRef} />
    </div>
  )
}
