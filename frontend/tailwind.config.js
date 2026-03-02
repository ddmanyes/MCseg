/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: {
    extend: {
      colors: {
        primary: { DEFAULT: '#3b82f6', dark: '#1d4ed8' },
        surface: { DEFAULT: '#1e1e2e', card: '#2a2a3e', border: '#3a3a5c' },
      },
    },
  },
  plugins: [],
}
