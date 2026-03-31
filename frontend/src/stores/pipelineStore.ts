import { create } from 'zustand'
import type { TaskStatus, RoiDefinition, LogMessage } from '../types/pipeline'

interface StageState extends TaskStatus {
  logs: string[]
}

const defaultStage = (): StageState => ({
  status: 'idle', progress: 0, message: '', logs: []
})

interface PipelineStore {
  stages: Record<string, StageState>
  rois: RoiDefinition[]
  updateStage: (stage: string, update: Partial<StageState>) => void
  appendLog: (stage: string, msg: string) => void
  setRois: (rois: RoiDefinition[]) => void
}

export const usePipelineStore = create<PipelineStore>((set) => ({
  stages: {
    roi: defaultStage(), segmentation: defaultStage(), count: defaultStage(),
    analysis: defaultStage(), spatial: defaultStage(), xenium: defaultStage(), loupe: defaultStage(),
  },
  rois: [],
  updateStage: (stage, update) => set((s) => ({
    stages: { ...s.stages, [stage]: { ...s.stages[stage], ...update } }
  })),

  appendLog: (stage, msg) => set((s) => ({
    stages: {
      ...s.stages,
      [stage]: {
        ...s.stages[stage],
        logs: [...(s.stages[stage]?.logs ?? []).slice(-500), msg],
      }
    }
  })),

  setRois: (rois) => set({ rois }),

}))
