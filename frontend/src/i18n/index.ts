import { useLanguageStore } from '../stores/languageStore'
import { translations } from './translations'

/** Returns a translation function `t(key, fallback?)` bound to the current language. */
export function useT() {
  const lang = useLanguageStore(s => s.lang)
  return (key: string, fallback?: string): string => {
    const entry = translations[key]
    if (!entry) return fallback ?? key
    return entry[lang] ?? entry.zh ?? fallback ?? key
  }
}

export function useLang() {
  return useLanguageStore(s => s.lang)
}
