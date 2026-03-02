import { create } from 'zustand'
import type { TaskStatus, ConditionResult, RoiDefinition, LogMessage } from '../types/pipeline'

interface StageState extends TaskStatus {
  logs: string[]
}

const defaultStage = (): StageState => ({
  status: 'idle', progress: 0, message: '', logs: []
})

interface PipelineStore {
  stages: Record<string, StageState>
  rois: RoiDefinition[]
  conditionResults: ConditionResult[]
  recommendedCondition: ConditionResult | null

  updateStage: (stage: string, update: Partial<StageState>) => void
  appendLog: (stage: string, msg: string) => void
  setRois: (rois: RoiDefinition[]) => void
  setConditionResults: (results: ConditionResult[]) => void
  setRecommendedCondition: (cond: ConditionResult | null) => void
}

export const usePipelineStore = create<PipelineStore>((set) => ({
  stages: {
    roi: defaultStage(), segmentation: defaultStage(), zarr: defaultStage(),
    conditions: defaultStage(), proseg: defaultStage(), analysis: defaultStage(),
    xenium: defaultStage(), loupe: defaultStage(),
  },
  rois: [],
  conditionResults: [],
  recommendedCondition: null,

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
  setConditionResults: (results) => set({ conditionResults: results }),
  setRecommendedCondition: (cond) => set({ recommendedCondition: cond }),
}))
