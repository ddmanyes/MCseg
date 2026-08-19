import { useLocation } from 'react-router-dom'
import { useT } from '../../i18n'

export default function Header() {
  const { pathname } = useLocation()
  const t = useT()

  const titles: Record<string, string> = {
    '/data': t('header.data'),
    '/roi': t('header.roi'),
    '/segmentation': t('header.seg'),
    '/count': t('header.count'),
    '/analysis': t('header.analysis'),
    '/spatial': t('header.spatial'),
    '/export': t('header.export'),
  }

  return (
    <header className="h-12 border-b border-surface-border flex items-center px-6 bg-surface-card flex-shrink-0">
      <h2 className="text-sm font-semibold text-gray-200">
        {titles[pathname] ?? 'MCseg'}
      </h2>
    </header>
  )
}

