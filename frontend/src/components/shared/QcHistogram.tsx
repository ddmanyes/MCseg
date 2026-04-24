/**
 * QcHistogram — 純 SVG 互動直方圖
 *
 * 功能：
 *  - 可拖曳的 min / max 閾值線（紅 / 橙）
 *  - MAD 建議線（綠虛線）
 *  - Log / Linear Y 軸切換
 *  - 即時通過率（閾值內的 bar 高亮）
 *  - 百分位統計列
 */

import { useState, useRef, useCallback, useEffect } from 'react'

// ── SVG layout constants ──────────────────────────────────────────
const PAD_L = 44
const PAD_R = 12
const PAD_T = 10
const PAD_B = 30
const CHART_H = 130
const TOTAL_H = CHART_H + PAD_T + PAD_B

// ── Types ─────────────────────────────────────────────────────────

export interface HistogramMetric {
  label: string
  unit: string
  bin_edges: number[]
  counts: number[]
  mad_min: number
  mad_max: number
  p5: number
  p50: number
  p95: number
  p99: number
  mean: number
}

interface Props {
  metric: HistogramMetric
  /** 目前的 min 閾值（null = 無限制） */
  minVal: number | null
  /** 目前的 max 閾值（null = 無限制） */
  maxVal: number | null
  onMinChange: (v: number | null) => void
  onMaxChange: (v: number | null) => void
  showMin?: boolean
  showMax?: boolean
  showMad?: boolean
  logScale?: boolean
  onLogScaleToggle?: () => void
  totalCells: number
}

// ── Helpers ───────────────────────────────────────────────────────

function fmtNum(v: number): string {
  if (!isFinite(v)) return '–'
  if (v >= 1e6) return (v / 1e6).toFixed(1) + 'M'
  if (v >= 1e3) return (v / 1e3).toFixed(1) + 'K'
  return Number.isInteger(v) ? String(v) : v.toFixed(1)
}

// ── Component ─────────────────────────────────────────────────────

export default function QcHistogram({
  metric,
  minVal,
  maxVal,
  onMinChange,
  onMaxChange,
  showMin = true,
  showMax = true,
  showMad = true,
  logScale = false,
  onLogScaleToggle,
  totalCells,
}: Props) {
  const svgRef = useRef<SVGSVGElement>(null)
  const [svgWidth, setSvgWidth] = useState(360)
  const dragging = useRef<'min' | 'max' | null>(null)

  // Responsive width via ResizeObserver
  useEffect(() => {
    const el = svgRef.current
    if (!el) return
    const ro = new ResizeObserver(entries => {
      setSvgWidth(Math.floor(entries[0].contentRect.width))
    })
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  const { bin_edges, counts } = metric
  const n = counts?.length ?? 0
  const chartW = svgWidth - PAD_L - PAD_R

  // 提供安全預設值讓 hooks 可以在資料無效時仍正常執行（Rules of Hooks 要求 hooks 不能在條件分支後）
  const isInvalid = !bin_edges || !counts || n === 0
    || bin_edges[n] === undefined
    || bin_edges[0] < 0
    || bin_edges[n] <= bin_edges[0]
  const xMin = isInvalid ? 0   : bin_edges[0]
  const xMax = isInvalid ? 100 : bin_edges[n]

  // X scale helpers — linear
  const xScaleLin  = useCallback((v: number) => PAD_L + ((v - xMin) / (xMax - xMin)) * chartW, [xMin, xMax, chartW])
  const xInvertLin = useCallback((px: number) => xMin + ((px - PAD_L) / chartW) * (xMax - xMin), [xMin, xMax, chartW])

  // X scale helpers — log₁₀(x+1)
  const logXMin = Math.log10(xMin + 1)
  const logXMax = Math.log10(xMax + 1)
  const xScaleLog  = useCallback((v: number) =>
    PAD_L + ((Math.log10(v + 1) - logXMin) / (logXMax - logXMin)) * chartW,
  [logXMin, logXMax, chartW])
  const xInvertLog = useCallback((px: number) => {
    const logV = logXMin + ((px - PAD_L) / chartW) * (logXMax - logXMin)
    return Math.max(0, Math.pow(10, logV) - 1)
  }, [logXMin, logXMax, chartW])

  const xScale  = logScale ? xScaleLog  : xScaleLin
  const xInvert = logScale ? xInvertLog : xInvertLin

  // Y scale
  const safeCounts = isInvalid ? [] : counts
  const yData = logScale ? safeCounts.map(c => Math.log10(c + 1)) : safeCounts
  const yMax  = Math.max(...yData, 1)
  const yScale = (v: number) => PAD_T + CHART_H - (v / yMax) * CHART_H

  // Pass rate (approximate via bins)
  const totalInBins = safeCounts.reduce((a, b) => a + b, 0)
  const passInBins = safeCounts.reduce((sum, c, i) => {
    const center = (bin_edges[i] + bin_edges[i + 1]) / 2
    const okMin = minVal === null || center >= minVal
    const okMax = maxVal === null || center <= maxVal
    return sum + (okMin && okMax ? c : 0)
  }, 0)
  const passRate  = totalInBins > 0 ? (passInBins / totalInBins) * 100 : 100
  const passCount = Math.round((passRate / 100) * totalCells)

  // Mouse helpers for dragging
  const getSvgPx = (e: React.MouseEvent<SVGSVGElement>) => {
    const rect = svgRef.current!.getBoundingClientRect()
    return e.clientX - rect.left
  }

  const handleMouseMove = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    if (!dragging.current) return
    const raw = xInvert(getSvgPx(e))
    const clamped = Math.max(xMin, Math.min(xMax, raw))
    const range = xMax - xMin
    const rounded = range <= 10 ? Number(clamped.toFixed(3)) : Math.round(clamped)
    if (dragging.current === 'min') onMinChange(rounded)
    else onMaxChange(rounded)
  }, [xMin, xMax, xInvert, onMinChange, onMaxChange])

  const handleMouseUp = useCallback(() => { dragging.current = null }, [])

  // ── 資料無效時顯示佔位符（所有 hooks 已在上方完成，此處可安全 return）
  if (isInvalid) {
    return (
      <div className="flex flex-col gap-1.5">
        <span className="text-xs font-semibold text-gray-200">{metric.label}</span>
        <div className="text-xs text-gray-500 italic py-4 text-center">無資料（細胞計數全為 0）</div>
      </div>
    )
  }

  // X-axis ticks
  const xTicks = logScale
    ? (() => {
        const ticks: { x: number; label: string }[] = []
        for (let d = 0; d <= 5; d++) {
          const v = d === 0 ? 0 : Math.pow(10, d) - 1
          if (v > xMax * 1.05) break
          const px = xScaleLog(v)
          if (px >= PAD_L - 2 && px <= PAD_L + chartW + 2)
            ticks.push({ x: px, label: d === 0 ? '0' : fmtNum(Math.pow(10, d)) })
        }
        return ticks
      })()
    : Array.from({ length: 5 }, (_, i) => {
        const v = xMin + (i / 4) * (xMax - xMin)
        return { x: xScaleLin(v), label: fmtNum(v) }
      })

  // Y-axis ticks
  const yTicks = logScale
    ? (() => {
        const ticks: { y: number; label: string }[] = []
        for (let d = 0; d <= 5; d++) {
          const count = d === 0 ? 1 : Math.pow(10, d)  // 1, 10, 100, 1K, 10K...
          const logVal = Math.log10(count + 1)
          if (logVal > yMax * 1.05) break
          const yPos = yScale(logVal)
          if (yPos >= PAD_T - 2 && yPos <= PAD_T + CHART_H + 2)
            ticks.push({ y: yPos, label: fmtNum(count) })
        }
        return ticks
      })()
    : [0, 0.33, 0.66, 1].map(f => ({
        y: PAD_T + CHART_H * (1 - f),
        label: fmtNum(f * yMax),
      }))

  // Threshold line renderer
  const ThresholdLine = ({
    val, color, type,
  }: { val: number | null; color: string; type: 'min' | 'max' }) => {
    if (val === null) return null
    const lx = xScale(val)
    if (lx < PAD_L - 2 || lx > PAD_L + chartW + 2) return null
    const labelX = type === 'min' ? lx + 5 : lx - 5
    const labelAnchor = type === 'min' ? 'start' : 'end'
    return (
      <g>
        <line x1={lx} y1={PAD_T} x2={lx} y2={PAD_T + CHART_H}
          stroke={color} strokeWidth={1.5} />
        <text x={labelX} y={PAD_T + 11} fill={color} fontSize={9} textAnchor={labelAnchor}>
          {type === 'min' ? 'min' : 'max'}={fmtNum(val)}
        </text>
        {/* Drag handle */}
        <circle
          cx={lx} cy={PAD_T + CHART_H / 2} r={7}
          fill={color} stroke="white" strokeWidth={1.5}
          style={{ cursor: 'ew-resize' }}
          onMouseDown={e => { e.stopPropagation(); e.preventDefault(); dragging.current = type }}
        />
      </g>
    )
  }

  const rateColor = passRate >= 80 ? 'text-green-400' : passRate >= 50 ? 'text-yellow-400' : 'text-red-400'

  return (
    <div className="flex flex-col gap-1.5">
      {/* Header row */}
      <div className="flex items-center justify-between gap-2 flex-wrap">
        <span className="text-xs font-semibold text-gray-200">
          {metric.label}{metric.unit ? <span className="text-gray-500 font-normal"> ({metric.unit})</span> : ''}
        </span>
        <div className="flex items-center gap-3">
          <span className="text-xs text-gray-400">
            <span className={rateColor}>{passRate.toFixed(1)}%</span>
            {' '}通過 ({passCount.toLocaleString()} / {totalCells.toLocaleString()})
          </span>
          {onLogScaleToggle && (
            <button
              onClick={onLogScaleToggle}
              className={`px-2 py-0.5 rounded text-xs border transition-colors ${
                logScale
                  ? 'bg-brand-primary/20 border-brand-primary text-brand-primary'
                  : 'border-gray-600 text-gray-500 hover:text-gray-300'
              }`}
            >
              Log
            </button>
          )}
        </div>
      </div>

      {/* SVG chart */}
      <svg
        ref={svgRef}
        className="w-full select-none cursor-default"
        style={{ height: TOTAL_H }}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
      >
        {/* Chart background */}
        <rect x={PAD_L} y={PAD_T} width={chartW} height={CHART_H}
          fill="rgba(15,20,30,0.6)" rx={3} />

        {/* Y-axis grid + ticks */}
        {yTicks.map((t, i) => (
          <g key={i}>
            <line x1={PAD_L} y1={t.y} x2={PAD_L + chartW} y2={t.y}
              stroke="#374151" strokeWidth={0.5} />
            <text x={PAD_L - 4} y={t.y + 3.5} textAnchor="end" fill="#6b7280" fontSize={8}>
              {t.label}
            </text>
          </g>
        ))}

        {/* Bars */}
        {counts.map((count, i) => {
          const x  = xScale(bin_edges[i])
          const x2 = xScale(bin_edges[i + 1])
          const bw = Math.max(1, x2 - x - 0.5)
          const h  = yScale(yData[i])
          const barH = CHART_H + PAD_T - h
          const center = (bin_edges[i] + bin_edges[i + 1]) / 2
          const inPass = (minVal === null || center >= minVal) && (maxVal === null || center <= maxVal)
          return (
            <rect
              key={i} x={x} y={h} width={bw} height={Math.max(0, barH)}
              fill={inPass ? '#38bdf8' : '#475569'}
              opacity={inPass ? 0.82 : 0.28}
            />
          )
        })}

        {/* MAD suggestion lines */}
        {showMad && metric.mad_min > xMin && metric.mad_min < xMax && (
          <g>
            <line
              x1={xScale(metric.mad_min)} y1={PAD_T}
              x2={xScale(metric.mad_min)} y2={PAD_T + CHART_H}
              stroke="#4ade80" strokeWidth={1} strokeDasharray="4,3" opacity={0.65}
            />
          </g>
        )}
        {showMad && metric.mad_max > xMin && metric.mad_max < xMax && (
          <line
            x1={xScale(metric.mad_max)} y1={PAD_T}
            x2={xScale(metric.mad_max)} y2={PAD_T + CHART_H}
            stroke="#4ade80" strokeWidth={1} strokeDasharray="4,3" opacity={0.65}
          />
        )}

        {/* Threshold lines (rendered on top) */}
        {showMin && <ThresholdLine val={minVal} color="#f87171" type="min" />}
        {showMax && <ThresholdLine val={maxVal} color="#fb923c" type="max" />}

        {/* X-axis line */}
        <line
          x1={PAD_L} y1={PAD_T + CHART_H}
          x2={PAD_L + chartW} y2={PAD_T + CHART_H}
          stroke="#4b5563" strokeWidth={0.5}
        />

        {/* X-axis ticks + labels */}
        {xTicks.map((t, i) => (
          <g key={i}>
            <line x1={t.x} y1={PAD_T + CHART_H} x2={t.x} y2={PAD_T + CHART_H + 4}
              stroke="#6b7280" strokeWidth={0.5} />
            <text x={t.x} y={PAD_T + CHART_H + 14} textAnchor="middle" fill="#9ca3af" fontSize={8}>
              {t.label}
            </text>
          </g>
        ))}

        {/* Y-axis label */}
        <text
          transform={`rotate(-90) translate(${-(PAD_T + CHART_H / 2)},${PAD_L - 33})`}
          textAnchor="middle" fill="#6b7280" fontSize={8}
        >
          {logScale ? 'cells (log)' : 'cells'}
        </text>
      </svg>

      {/* Stats + MAD row */}
      <div className="flex gap-3 text-xs text-gray-500 flex-wrap">
        <span>P50: <span className="text-gray-300">{fmtNum(metric.p50)}</span></span>
        <span>P5: <span className="text-gray-300">{fmtNum(metric.p5)}</span></span>
        <span>P95: <span className="text-gray-300">{fmtNum(metric.p95)}</span></span>
        <span>P99: <span className="text-gray-300">{fmtNum(metric.p99)}</span></span>
        {showMad && (
          <span className="text-green-600/80">
            MAD建議: {fmtNum(metric.mad_min)}–{fmtNum(metric.mad_max)}
          </span>
        )}
      </div>
    </div>
  )
}
