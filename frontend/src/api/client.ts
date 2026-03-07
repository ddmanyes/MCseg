import axios from 'axios'
import type { RoiDefinition } from '../types/pipeline'

const api = axios.create({ baseURL: '/api' })

// Health
export const getHealth = () => api.get('/health')

// Config
export const getConfig = () => api.get('/config')

// Data Setup
export const scanData = (body: { data_root: string }) => api.post('/data/scan', body)
export const applyData = (paths: object) => api.post('/data/apply', paths)
export const getDataStatus = () => api.get('/data/status')
export const browseDir = (path: string) => api.get('/data/browse', { params: { path } })
export const getDiskStatus = () => api.get('/data/disk-status')

// Stage 0: ROI
export const listRois = () => api.get('/roi/list')
export const addRoi = (roi: RoiDefinition) => api.post('/roi/add', roi)
export const deleteRoi = (name: string) => api.delete(`/roi/${name}`)
export const getRoiOverview = () => api.get('/roi/overview')
export const runRoiExtract = () => api.post('/roi/extract')
export const getRoiStatus = () => api.get('/roi/status')

// Stage 1: Segmentation
export const runSegmentation = (params?: object) => api.post('/segmentation/run', params ?? {})
export const getSegmentationStatus = () => api.get('/segmentation/status')
export const getSegmentationPreview = (roi?: string) => api.get('/segmentation/preview', { params: roi ? { roi_name: roi } : {} })
export const runSegmentationPreview = (body: object) => api.post('/segmentation/run_preview', body)
export const previewPreproc = (body: object) => api.post('/segmentation/preview_preproc', body)

// Stage 2: Zarr
export const buildZarr = () => api.post('/zarr/build')
export const getZarrStatus = () => api.get('/zarr/status')

// Stage 2.5: Conditions
export const runConditions = (body: object) => api.post('/conditions/run', body)
export const getConditionsStatus = () => api.get('/conditions/status')
export const getConditionsResults = () => api.get('/conditions/results')
export const getConditionsRecommend = () => api.get('/conditions/recommend')
export const getConditionThumbnail = (idx: number) => api.get(`/conditions/thumbnail/${idx}`)
export const getConditionThumbnailHd = (idx: number) => api.get(`/conditions/thumbnail_hd/${idx}`)

// Stage 3: Proseg
export const runProseg = (params?: { max_dist?: number; compactness?: number; dilation?: number }) =>
  api.post('/proseg/run', params ?? {})
export const getProsegStatus = () => api.get('/proseg/status')

// Stage 4: Analysis (舊版整合執行)
export const runAnalysis = (params?: object) => api.post('/analysis/run', params ?? {})
export const getAnalysisStatus = () => api.get('/analysis/status')
export const getUmap = () => api.get('/analysis/umap')

// Stage 4: Step 1 — QC 前處理
export const runQC = (params?: object) => api.post('/analysis/run_qc', params ?? {})
export const getQCStatus = () => api.get('/analysis/qc_status')
export const getQCImages = () => api.get('/analysis/qc_images')
export const getOverlayHdUrl = (name: 'pre_qc' | 'post_qc') => `/api/analysis/overlay_hd/${name}`

// Stage 4: Step 2 — UMAP 多解析度
export const runUMAPExplore = (params?: object) => api.post('/analysis/run_umap', params ?? {})
export const getUMAPExploreStatus = () => api.get('/analysis/umap_status')
export const getUMAPImages = () => api.get('/analysis/umap_images')

// Stage 4: Step 3 — Heatmap
export const runHeatmap = (params: object) => api.post('/analysis/run_heatmap', params)
export const getHeatmapStatus = () => api.get('/analysis/heatmap_status')
export const getHeatmapImage = () => api.get('/analysis/heatmap')

// Stage 5: Export
export const exportXenium = (body: object) => api.post('/export/xenium', body)
export const exportLoupe = (body: object) => api.post('/export/loupe', body)
export const getXeniumStatus = () => api.get('/export/status/xenium')
export const getLoupeStatus = () => api.get('/export/status/loupe')
