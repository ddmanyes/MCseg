import { useEffect } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import Sidebar from './components/layout/Sidebar'
import Header from './components/layout/Header'
import PipelineStepper from './components/layout/PipelineStepper'
import DataSetup from './pages/DataSetup'
import Stage0_ROI from './pages/Stage0_ROI'
import Stage1_Segmentation from './pages/Stage1_Segmentation'
import Stage2_Zarr from './pages/Stage2_Zarr'
import Stage2b_ConditionTest from './pages/Stage2b_ConditionTest'
import Stage3_Proseg from './pages/Stage3_Proseg'
import Stage4_Analysis from './pages/Stage4_Analysis'
import Stage5_Export from './pages/Stage5_Export'
import { getDiskStatus } from './api/client'
import { usePipelineStore } from './stores/pipelineStore'

export default function App() {
  const { updateStage } = usePipelineStore()

  // 啟動時掃描磁碟，恢復各 Stage 完成狀態
  useEffect(() => {
    getDiskStatus().then(res => {
      const d = res.data?.data
      if (!d) return
      if (d.roi?.done)          updateStage('roi',          { status: 'done', message: `已完成（${d.roi.roi_names?.length ?? 0} 個 ROI）` })
      if (d.segmentation?.done) updateStage('segmentation', { status: 'done', message: '分割遮罩已存在' })
      if (d.zarr?.done)         updateStage('zarr',         { status: 'done', message: 'Zarr 已建構' })
      if (d.proseg?.done)       updateStage('proseg',       { status: 'done', message: 'Proseg 已執行' })
      if (d.analysis?.done)     updateStage('analysis',     { status: 'done', message: '分析已完成' })
    }).catch(() => {/* 靜默失敗，不影響正常使用 */})
  }, [])

  return (
    <div className="flex h-screen overflow-hidden bg-surface">
      <Sidebar />
      <div className="flex-1 flex flex-col overflow-hidden">
        <Header />
        <PipelineStepper />
        <main className="flex-1 overflow-y-auto p-6">
          <Routes>
            <Route path="/" element={<Navigate to="/data" replace />} />
            <Route path="/data" element={<DataSetup />} />
            <Route path="/roi" element={<Stage0_ROI />} />
            <Route path="/segmentation" element={<Stage1_Segmentation />} />
            <Route path="/zarr" element={<Stage2_Zarr />} />
            <Route path="/conditions" element={<Stage2b_ConditionTest />} />
            <Route path="/proseg" element={<Stage3_Proseg />} />
            <Route path="/analysis" element={<Stage4_Analysis />} />
            <Route path="/export" element={<Stage5_Export />} />
          </Routes>
        </main>
      </div>
    </div>
  )
}
