import { useLocation } from 'react-router-dom'

const TITLES: Record<string, string> = {
  '/data': '📂 資料設定',
  '/roi': 'Stage 0 — ROI 定義與裁切',
  '/segmentation': 'Stage 1 — 細胞分割（Cellpose）',
  '/zarr': 'Stage 2 — Zarr 建構',
  '/conditions': 'Stage 2.5 — Proseg 參數條件測試',
  '/proseg': 'Stage 3 — Proseg RNA 重新分配',
  '/analysis': 'Stage 4 — 下游聚類分析',
  '/export': 'Stage 5 — Browser 格式匯出',
}

export default function Header() {
  const { pathname } = useLocation()
  return (
    <header className="h-12 border-b border-surface-border flex items-center px-6 bg-surface-card flex-shrink-0">
      <h2 className="text-sm font-semibold text-gray-200">
        {TITLES[pathname] ?? 'VisiumHD Pipeline 2'}
      </h2>
    </header>
  )
}
