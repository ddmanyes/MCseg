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
export const getOutputDir = () => api.get('/data/output-dir')
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
export const getRoiSegOverrides = () => api.get('/segmentation/roi_overrides')
export const saveRoiSegOverrides = (overrides: Record<string, Record<string, unknown>>) =>
  api.put('/segmentation/roi_overrides', overrides)

// Stage 2: Cellpose RNA 計數
export const runCellposeCount = (roiName: string | null) =>
  api.post('/count/run', roiName ? { roi_name: roiName } : {})
export const getCellposeCountStatus = () => api.get('/count/status')
export const listCountRois = () => api.get('/count/available_rois')

// Stage 2.5: Proseg RNA 重分配
export const runProsegRNA = (roiName: string | null) =>
  api.post('/proseg_rna/run', roiName ? { roi_name: roiName } : {})
export const getProsegRNAStatus = () => api.get('/proseg_rna/status')
export const listProsegRNARois = () => api.get('/proseg_rna/available_rois')
export const getProsegComparison = (roiName: string) => api.get(`/proseg_rna/comparison/${roiName}`)

// Stage 3: Analysis (舊版整合執行)
export const runAnalysis = (params?: object) => api.post('/analysis/run', params ?? {})
export const getAnalysisStatus = () => api.get('/analysis/status')
export const getUmap = () => api.get('/analysis/umap')

// Stage 3: 原始分布直方圖
export const getRawHistogram = (roiName?: string, mergeRois?: boolean, source?: string) =>
  api.get('/analysis/raw_histogram', { params: { roi_name: roiName, merge_rois: mergeRois, source } })

// Stage 3: Step 1 — QC 前處理
export const runQC = (params?: object) => api.post('/analysis/run_qc', params ?? {})
export const getQCStatus = () => api.get('/analysis/qc_status')
export const getQCImages = () => api.get('/analysis/qc_images')
export const getOverlayHdUrl = (name: 'pre_qc' | 'post_qc') => `/api/analysis/overlay_hd/${name}`
export const getAvailableRois = () => api.get('/analysis/available_rois')
export const getRoiOverlays = () => api.get('/analysis/roi_overlays')

// Stage 3: Step 2 — UMAP 多解析度
export const runUMAPExplore = (params?: object) => api.post('/analysis/run_umap', params ?? {})
export const getUMAPExploreStatus = () => api.get('/analysis/umap_status')
export const getUMAPImages = () => api.get('/analysis/umap_images')

// Stage 3: Step 3 — 細胞類型標註
export const getClusterInfo = (resolution: number) =>
  api.get('/analysis/cluster_info', { params: { resolution } })
export const getCelltypistModels = () => api.get('/analysis/celltypist_models')
export const runAnnotate = (params: object) => api.post('/analysis/annotate', params)
export const getAnnotateStatus = () => api.get('/analysis/annotate_status')
export const getAnnotateSuggestions = () => api.get('/analysis/annotate_suggestions')
export const applyLabels = (params: object) => api.post('/analysis/apply_labels', params)

// Stage 3: Step 4 — Heatmap
export const runHeatmap = (params: object) => api.post('/analysis/run_heatmap', params)
export const getHeatmapStatus = () => api.get('/analysis/heatmap_status')
export const getHeatmapImage = () => api.get('/analysis/heatmap')

// Stage 4: Export
export const exportXenium = (body: object) => api.post('/export/xenium', body)
export const exportLoupe = (body: object) => api.post('/export/loupe', body)
export const getXeniumStatus = () => api.get('/export/status/xenium')
export const getLoupeStatus = () => api.get('/export/status/loupe')
