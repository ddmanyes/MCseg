import { Routes, Route, Navigate } from 'react-router-dom'
import Sidebar from './components/layout/Sidebar'
import Header from './components/layout/Header'
import DataSetup from './pages/DataSetup'
import Stage0_ROI from './pages/Stage0_ROI'
import Stage1_Segmentation from './pages/Stage1_Segmentation'
import Stage2_Zarr from './pages/Stage2_Zarr'
import Stage2b_ConditionTest from './pages/Stage2b_ConditionTest'
import Stage3_Proseg from './pages/Stage3_Proseg'
import Stage4_Analysis from './pages/Stage4_Analysis'
import Stage5_Export from './pages/Stage5_Export'

export default function App() {
  return (
    <div className="flex h-screen overflow-hidden bg-surface">
      <Sidebar />
      <div className="flex-1 flex flex-col overflow-hidden">
        <Header />
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
