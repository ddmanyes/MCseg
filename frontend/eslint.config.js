import js from '@eslint/js'
import globals from 'globals'
import reactHooks from 'eslint-plugin-react-hooks'
import tseslint from 'typescript-eslint'

export default tseslint.config(
  // `._*` 是 macOS 在 ExFAT 外接碟上產生的 AppleDouble 雜訊，不是原始碼
  { ignores: ['dist', 'node_modules', '**/._*'] },
  {
    files: ['src/**/*.{ts,tsx}'],
    extends: [
      js.configs.recommended,
      ...tseslint.configs.recommendedTypeChecked,
    ],
    languageOptions: {
      ecmaVersion: 2020,
      globals: globals.browser,
      parserOptions: {
        project: ['./tsconfig.json'],
        tsconfigRootDir: import.meta.dirname,
      },
    },
    plugins: {
      'react-hooks': reactHooks,
    },
    rules: {
      ...reactHooks.configs.recommended.rules,

      // === 本計畫的三條核心規則（靜默出錯防線）===
      'no-empty': ['error', { allowEmptyCatch: false }],
      'no-empty-function': 'off',
      '@typescript-eslint/no-empty-function': 'error',
      '@typescript-eslint/no-floating-promises': 'error',

      // 既有程式碼大量使用 any／未型別化的 API 回應，先降為警告避免淹沒上面三條
      '@typescript-eslint/no-explicit-any': 'off',
      '@typescript-eslint/no-unsafe-assignment': 'off',
      '@typescript-eslint/no-unsafe-member-access': 'off',
      '@typescript-eslint/no-unsafe-call': 'off',
      '@typescript-eslint/no-unsafe-argument': 'off',
      '@typescript-eslint/no-unsafe-return': 'off',
      '@typescript-eslint/no-unused-vars': ['warn', { argsIgnorePattern: '^_' }],
      '@typescript-eslint/restrict-template-expressions': 'off',
      '@typescript-eslint/no-misused-promises': [
        'error',
        { checksVoidReturn: { attributes: false } },
      ],
    },
  },
)
