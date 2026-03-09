/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: {
    extend: {
      colors: {
        primary:  { DEFAULT: '#3b82f6', dark: '#1d4ed8' },
        surface:  {
          DEFAULT: '#0f0f17',   // Xenium Explorer 風格深黑背景
          card:    '#161621',   // 卡片背景
          border:  '#ffffff0f', // 白色 6% 邊框
        },
      },
    },
  },
  plugins: [],
}
