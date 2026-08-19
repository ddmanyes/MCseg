import { Component, type ErrorInfo, type ReactNode } from 'react'
import { useLanguageStore } from '../../stores/languageStore'
import { translations } from '../../i18n/translations'

interface Props {
  children: ReactNode
  /** 顯示在錯誤卡片上的位置標籤，例如頁面路徑 */
  label?: string
}

interface State {
  error: string
  stack: string
}

/**
 * 全域 render 例外攔截。
 *
 * 沒有這層的話，任何一頁 render 期間拋錯（例如後端回傳結構不符、
 * `qc.sections.map` 打在 undefined 上）會炸掉整個 React tree → 全白畫面，
 * 使用者連錯在哪一頁都看不到。
 *
 * 在 `App.tsx` 以 `key={location.pathname}` 掛載，換頁時會重新 mount，
 * 因此壞掉的頁面不會把整個 app 永久卡在錯誤畫面。
 */
export default class ErrorBoundary extends Component<Props, State> {
  state: State = { error: '', stack: '' }

  static getDerivedStateFromError(e: Error): State {
    return { error: e.message || String(e), stack: e.stack ?? '' }
  }

  componentDidCatch(e: Error, info: ErrorInfo) {
    // 保留完整堆疊到 console，畫面上只顯示摘要
    console.error('[ErrorBoundary]', this.props.label ?? '', e, info.componentStack)
  }

  render() {
    if (this.state.error) {
      const lang = useLanguageStore.getState().lang
      const title = translations['error.boundary.title']?.[lang] ?? '✗ An error occurred while rendering this page'
      const subtitle = translations['error.boundary.subtitle']?.[lang] ?? 'Other pages are still working normally. Switch pages to leave this screen. See browser console for details.'

      return (
        <div className="m-4 rounded-xl border border-red-800 bg-red-900/20 p-4 text-sm text-red-300 space-y-2">
          <p className="font-semibold text-red-200">
            {title}{this.props.label ? ` (${this.props.label})` : ''}
          </p>
          <p className="font-mono text-xs whitespace-pre-wrap break-all">{this.state.error}</p>
          <p className="text-xs text-red-400/80">
            {subtitle}
          </p>
        </div>
      )
    }
    return this.props.children
  }
}

