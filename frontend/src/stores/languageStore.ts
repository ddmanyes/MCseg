import { create } from 'zustand'
import { persist } from 'zustand/middleware'

export type Lang = 'zh' | 'en'

interface LanguageStore {
  lang: Lang
  setLang: (lang: Lang) => void
}

export const useLanguageStore = create<LanguageStore>()(
  persist(
    (set) => ({
      lang: 'zh',
      setLang: (lang) => set({ lang }),
    }),
    { name: 'msseg-lang' }
  )
)
