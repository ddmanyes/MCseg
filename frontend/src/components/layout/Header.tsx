import { useLocation } from 'react-router-dom'

const TITLES: Record<string, string> = {
  '/data': '📂 資料設定',
  '/roi': 'Stage 0 — ROI 定義與裁切',
  '/segmentation': 'Stage 1 — 細胞分割（MCseg v2）',
  '/count': 'Stage 2 — RNA 計數',
  '/analysis': 'Stage 3 — 下游聚類分析',
  '/export': 'Stage 4 — Browser 格式匯出',
}

export default function Header() {
  const { pathname } = useLocation()
  return (
    <header className="h-12 border-b border-surface-border flex items-center px-6 bg-surface-card flex-shrink-0">
      <h2 className="text-sm font-semibold text-gray-200">
        {TITLES[pathname] ?? 'MCseg'}
      </h2>
    </header>
  )
}
