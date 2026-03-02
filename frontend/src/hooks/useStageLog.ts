import { useEffect } from 'react'
import { usePipelineStore } from '../stores/pipelineStore'

export default function useStageLog(stage: string) {
  const appendLog = usePipelineStore((s) => s.appendLog)

  useEffect(() => {
    const protocol = location.protocol === 'https:' ? 'wss' : 'ws'
    const ws = new WebSocket(`${protocol}://${location.host}/ws/log/${stage}`)

    ws.onmessage = (e) => {
      try {
        const msg = JSON.parse(e.data)
        if (msg.type === 'log') appendLog(stage, `[${msg.level}] ${msg.message}`)
      } catch {}
    }

    return () => ws.close()
  }, [stage, appendLog])
}
