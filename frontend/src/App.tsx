import { useEffect } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import TopNav from './components/layout/TopNav'
import DataSetup from './pages/DataSetup'
import Stage0_ROI from './pages/Stage0_ROI'
import Stage1_Segmentation from './pages/Stage1_Segmentation'
import Stage2_Count from './pages/Stage2_Count'
import Stage3_Analysis from './pages/Stage3_Analysis'
import Stage35_SpatialExplorer from './pages/Stage35_SpatialExplorer'
import Stage4_Export from './pages/Stage4_Export'
import { getDiskStatus } from './api/client'
import { usePipelineStore } from './stores/pipelineStore'

export default function App() {
  const { updateStage } = usePipelineStore()

  // 啟動時掃描磁碟，恢復各 Stage 完成狀態
  useEffect(() => {
    getDiskStatus().then(res => {
      const d = res.data?.data
      if (!d) return
      if (d.roi?.done)          updateStage('roi',          { status: 'done', message: `Done (${d.roi.roi_names?.length ?? 0} ROIs)` })
      if (d.segmentation?.done) updateStage('segmentation', { status: 'done', message: 'Segmentation masks found' })
      if (d.count?.done)        updateStage('count',        { status: 'done', message: 'RNA counting complete' })
      if (d.analysis?.done)     updateStage('analysis',     { status: 'done', message: 'Analysis complete' })
    }).catch(() => {/* 靜默失敗 */})
  }, [])

  return (
    <div className="flex flex-col h-screen overflow-hidden bg-surface">
      <TopNav />
      <main className="flex-1 overflow-y-auto p-6">
        <Routes>
          <Route path="/" element={<Navigate to="/data" replace />} />
          <Route path="/data" element={<DataSetup />} />
          <Route path="/roi" element={<Stage0_ROI />} />
          <Route path="/segmentation" element={<Stage1_Segmentation />} />
          <Route path="/count" element={<Stage2_Count />} />
          <Route path="/analysis" element={<Stage3_Analysis />} />
          <Route path="/spatial" element={<Stage35_SpatialExplorer />} />
          <Route path="/export" element={<Stage4_Export />} />
        </Routes>
      </main>
    </div>
  )
}

