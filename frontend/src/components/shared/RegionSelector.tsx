import { useEffect, useRef, useState } from 'react'
import OpenSeadragon from 'openseadragon'
import { clsx } from 'clsx'
import { useT } from '../../i18n'

/**
 * 在全片影像上框選區域（矩形或多邊形）。
 *
 * 輸出座標一律為**全片 fullres px**，與後端 `_filter_by_region` 的座標系一致。
 * 由 Stage 0 的 ROI 畫框邏輯提取而來，供 Stage 0（定義 ROI）與
 * Stage 3.5（框選區域做圖）共用。
 */

export interface BoxOverlay {
  name?: string
  x: number
  y: number
  width_px: number
  height_px: number
}

export type RegionSelection =
  | { type: 'bbox'; x: number; y: number; width_px: number; height_px: number }
  | { type: 'polygon'; points: [number, number][] }

interface Props {
  onChange: (sel: RegionSelection | null) => void
  mode?: 'bbox' | 'polygon' | 'both'
  overlays?: BoxOverlay[]
  dziUrl?: string
  tilesUrl?: string
  height?: string
}

export default function RegionSelector({
  onChange,
  mode = 'bbox',
  overlays = [],
  dziUrl = '/api/roi/dzi',
  tilesUrl = '/api/roi/dzi_files/',
  height = '26rem',
}: Props) {
  const containerRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<OpenSeadragon.Viewer | null>(null)
  const [tool, setTool] = useState<'pan' | 'bbox' | 'polygon'>('pan')
  const [drawing, setDrawing] = useState(false)
  const [drawBox, setDrawBox] = useState<{ x: number; y: number; w: number; h: number } | null>(null)
  const [poly, setPoly] = useState<{ x: number; y: number }[]>([])
  const [ready, setReady] = useState(false)
  const startPt = useRef<{ x: number; y: number } | null>(null)
  const t = useT()

  const allowBbox = mode === 'bbox' || mode === 'both'
  const allowPoly = mode === 'polygon' || mode === 'both'

  // ── 初始化 OSD（mount 一次）────────────────────────────────────
  useEffect(() => {
    if (!containerRef.current) return
    let active = true

    const stamp = Date.now()
    fetch(`${dziUrl}?t=${stamp}`)
      .then(res => res.text())
      .then(xmlStr => {
        if (!active || !containerRef.current) return
        const doc = new DOMParser().parseFromString(xmlStr, 'application/xml')
        const image = doc.getElementsByTagName('Image')[0]
        const size = doc.getElementsByTagName('Size')[0]
        if (!image || !size) {
          console.error('Failed to parse DZI XML', xmlStr)
          return
        }

        const viewer = new OpenSeadragon.Viewer({
          element: containerRef.current,
          showNavigationControl: false,
          showNavigator: true,
          navigatorPosition: 'BOTTOM_RIGHT',
          gestureSettingsMouse: { clickToZoom: false, dblClickToZoom: true, scrollToZoom: true },
          tileSources: {
            Image: {
              xmlns: 'http://schemas.microsoft.com/deepzoom/2008',
              Url: tilesUrl,
              Format: image.getAttribute('Format') || 'jpeg',
              Overlap: (Number(image.getAttribute('Overlap')) || 1).toString(),
              TileSize: (Number(image.getAttribute('TileSize')) || 256).toString(),
              Size: {
                Width: String(Number(size.getAttribute('Width'))),
                Height: String(Number(size.getAttribute('Height'))),
              },
            },
          },
        })

        // 注入 cache-buster：換樣本後 tile URL 不變，不加會拿到上一片的快取影像
        const patchTiledImage = () => {
          const tiledImage = viewer.world.getItemAt(0)
          if (tiledImage && tiledImage.source && !(tiledImage.source as any)._patched) {
            const original = tiledImage.source.getTileUrl.bind(tiledImage.source)
            tiledImage.source.getTileUrl = (level: number, x: number, y: number) =>
              `${original(level, x, y)}?t=${stamp}`
            ;(tiledImage.source as any)._patched = true
          }
        }
        if (viewer.world.getItemCount() > 0) patchTiledImage()
        else viewer.world.addHandler('add-item', patchTiledImage)

        viewer.addHandler('open', () => { if (active) setReady(true) })
        viewerRef.current = viewer
      })
      .catch(e => console.error('Failed to fetch DZI:', e))

    return () => {
      active = false
      if (viewerRef.current) {
        viewerRef.current.destroy()
        viewerRef.current = null
      }
      setReady(false)
    }
  }, [dziUrl, tilesUrl])

  // ── 重繪既有覆蓋框 ─────────────────────────────────────────────
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return

    viewer.clearOverlays()
    for (const box of overlays) {
      if (box.x == null || box.width_px == null) continue

      const el = document.createElement('div')
      el.style.border = '2px solid #61dafb'
      el.style.background = 'rgba(97, 218, 251, 0.08)'
      el.style.pointerEvents = 'none'

      const label = document.createElement('span')
      label.textContent = box.name ?? ''
      label.style.cssText = [
        'position:absolute', 'top:2px', 'left:4px',
        'font-size:10px', 'color:#61dafb', 'font-weight:600',
        'text-shadow:0 0 4px #000', 'white-space:nowrap',
      ].join(';')
      el.appendChild(label)

      viewer.addOverlay({
        element: el,
        location: viewer.viewport.imageToViewportRectangle(
          new OpenSeadragon.Rect(box.x, box.y, box.width_px, box.height_px),
        ),
      })
    }
  }, [overlays, ready])

  // 螢幕座標 → 全片影像座標
  const toImage = (cx: number, cy: number) => {
    const viewer = viewerRef.current!
    return viewer.viewport.viewportToImageCoordinates(
      viewer.viewport.pointFromPixel(new OpenSeadragon.Point(cx, cy)),
    )
  }

  const localPt = (e: React.PointerEvent<HTMLDivElement>) => {
    const rect = containerRef.current!.getBoundingClientRect()
    return { x: e.clientX - rect.left, y: e.clientY - rect.top }
  }

  // ── bbox ──────────────────────────────────────────────────────
  const handlePointerDown = (e: React.PointerEvent<HTMLDivElement>) => {
    if (tool !== 'bbox') return
    const p = localPt(e)
    startPt.current = p
    setDrawBox({ x: p.x, y: p.y, w: 0, h: 0 })
    setDrawing(true)
    ;(e.target as HTMLElement).setPointerCapture(e.pointerId)
  }

  const handlePointerMove = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!drawing || !startPt.current) return
    const p = localPt(e)
    const dx = p.x - startPt.current.x
    const dy = p.y - startPt.current.y
    setDrawBox({
      x: dx >= 0 ? startPt.current.x : p.x,
      y: dy >= 0 ? startPt.current.y : p.y,
      w: Math.abs(dx),
      h: Math.abs(dy),
    })
  }

  const handlePointerUp = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!drawing || !drawBox || !startPt.current) return
    ;(e.target as HTMLElement).releasePointerCapture(e.pointerId)
    setDrawing(false)

    if (drawBox.w > 10 && drawBox.h > 10 && viewerRef.current) {
      const tl = toImage(drawBox.x, drawBox.y)
      const br = toImage(drawBox.x + drawBox.w, drawBox.y + drawBox.h)
      onChange({
        type: 'bbox',
        x: Math.round(Math.min(tl.x, br.x)),
        y: Math.round(Math.min(tl.y, br.y)),
        width_px: Math.round(Math.abs(br.x - tl.x)),
        height_px: Math.round(Math.abs(br.y - tl.y)),
      })
    }
    setDrawBox(null)
  }

  // ── polygon ───────────────────────────────────────────────────
  const handlePolyClick = (e: React.PointerEvent<HTMLDivElement>) => {
    if (tool !== 'polygon') return
    setPoly(pts => [...pts, localPt(e)])
  }

  const finishPolygon = () => {
    if (poly.length < 3 || !viewerRef.current) return
    const points = poly.map(p => {
      const ip = toImage(p.x, p.y)
      return [Math.round(ip.x), Math.round(ip.y)] as [number, number]
    })
    onChange({ type: 'polygon', points })
    setPoly([])
  }

  const clearSelection = () => {
    setPoly([])
    setDrawBox(null)
    onChange(null)
  }

  const toolBtn = (value: 'pan' | 'bbox' | 'polygon', label: string) => (
    <button
      key={value}
      onClick={() => { setTool(value); setPoly([]) }}
      className={clsx(
        'px-3 py-1 rounded text-xs font-medium transition-colors',
        tool === value ? 'bg-primary text-black' : 'bg-surface-border text-gray-300 hover:bg-surface-border/80',
      )}
    >
      {label}
    </button>
  )

  return (
    <div className="space-y-2">
      <div className="flex items-center gap-2 flex-wrap">
        {toolBtn('pan', t('roi.pan_mode'))}
        {allowBbox && toolBtn('bbox', t('roi.draw_mode'))}
        {allowPoly && toolBtn('polygon', t('region.polygon_mode'))}

        {tool === 'polygon' && (
          <>
            <button
              onClick={finishPolygon}
              disabled={poly.length < 3}
              className="px-3 py-1 rounded text-xs font-medium bg-surface-border text-gray-300 hover:bg-surface-border/80 disabled:opacity-40"
            >
              {t('region.finish_polygon')} ({poly.length})
            </button>
            <button
              onClick={clearSelection}
              className="px-3 py-1 rounded text-xs font-medium bg-surface-border text-gray-300 hover:bg-surface-border/80"
            >
              {t('region.clear')}
            </button>
          </>
        )}

        <span className="text-xs text-gray-500">
          {tool === 'pan'
            ? t('roi.hint.pan')
            : tool === 'bbox'
              ? t('roi.hint.draw')
              : t('region.hint.polygon')}
        </span>
        {!ready && <span className="text-xs text-yellow-400 animate-pulse">{t('roi.loading')}</span>}
      </div>

      <div className="relative rounded-lg overflow-hidden border border-surface-border">
        <div ref={containerRef} className="w-full bg-black" style={{ height }} />

        {/* 繪製模式：覆蓋層攔截 pointer events，阻止 OSD pan */}
        {tool !== 'pan' && (
          <div
            className="absolute inset-0 cursor-crosshair"
            style={{ touchAction: 'none' }}
            onPointerDown={tool === 'bbox' ? handlePointerDown : handlePolyClick}
            onPointerMove={handlePointerMove}
            onPointerUp={handlePointerUp}
            onPointerCancel={handlePointerUp}
          >
            {drawing && drawBox && (
              <div
                className="absolute border-2 border-primary bg-primary/10 pointer-events-none"
                style={{ left: drawBox.x, top: drawBox.y, width: drawBox.w, height: drawBox.h }}
              >
                <span className="absolute -top-5 left-0 bg-primary text-black text-[10px] px-1 rounded-t font-mono whitespace-nowrap">
                  {drawBox.w} × {drawBox.h} px
                </span>
              </div>
            )}

            {tool === 'polygon' && poly.length > 0 && (
              <svg className="absolute inset-0 w-full h-full pointer-events-none">
                <polygon
                  points={poly.map(p => `${p.x},${p.y}`).join(' ')}
                  fill="rgba(97, 218, 251, 0.12)"
                  stroke="#61dafb"
                  strokeWidth={2}
                />
                {poly.map((p, i) => (
                  <circle key={i} cx={p.x} cy={p.y} r={3} fill="#61dafb" />
                ))}
              </svg>
            )}
          </div>
        )}
      </div>
    </div>
  )
}
