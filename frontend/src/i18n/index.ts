import { useLanguageStore } from '../stores/languageStore'
import { translations } from './translations'

/** Returns a translation function `t(key, fallbackOrVars?)` bound to the current language. */
export function useT() {
  const lang = useLanguageStore(s => s.lang)
  return (key: string, fallbackOrVars?: string | Record<string, string | number>): string => {
    const entry = translations[key]
    let text = entry ? (entry[lang] ?? entry.zh ?? key) : key
    if (typeof fallbackOrVars === 'string') {
      if (!entry) text = fallbackOrVars
    } else if (typeof fallbackOrVars === 'object' && fallbackOrVars !== null) {
      for (const [k, v] of Object.entries(fallbackOrVars)) {
        text = text.replace(new RegExp(`\\{${k}\\}`, 'g'), String(v))
      }
    }
    return text
  }
}

export function useLang() {
  return useLanguageStore(s => s.lang)
}

