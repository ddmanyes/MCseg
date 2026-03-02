export type StageStatus = 'idle' | 'running' | 'done' | 'error'

export interface TaskStatus {
  status: StageStatus
  progress: number
  message: string
  completed?: number
  total?: number
}

export interface RoiDefinition {
  name: string
  tissue: string
  x?: number
  y?: number
  width_px?: number
  height_px?: number
  pixel_size_um?: number
  x_um?: number
  y_um?: number
  width_um?: number
  height_um?: number
}

export interface ConditionResult {
  condition_idx: number
  label: string
  max_dist: number
  compactness: number
  dilation: number
  n_cells: number
  median_genes: number
  median_counts: number
  fraction_assigned: number
  cell_area_cv: number
  status: string
}

export interface LogMessage {
  type: 'log' | 'ping'
  stage: string
  message: string
  level: string
}
