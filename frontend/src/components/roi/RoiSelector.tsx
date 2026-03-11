import { useEffect, useRef, useState } from 'react'
import OpenSeadragon from 'openseadragon'
import { clsx } from 'clsx'

interface RoiBox {
  name?: string
  x: number
  y: number
  width_px: number
  height_px: number
}

interface Props {
  onSelect: (roi: Omit<RoiBox, 'name'>) => void
  existingRois?: RoiBox[]
}

export default function RoiSelector({ onSelect, existingRois = [] }: Props) {
  const containerRef = useRef<HTMLDivElement>(null)
  const viewerRef    = useRef<OpenSeadragon.Viewer | null>(null)
  const [mode, setMode]       = useState<'pan' | 'draw'>('pan')
  const [drawing, setDrawing] = useState(false)
  const [drawBox, setDrawBox] = useState<{ x: number; y: number; w: number; h: number } | null>(null)
  const [ready, setReady]     = useState(false)
  const startPt = useRef<{ x: number; y: number } | null>(null)

  // ── 初始化 OSD（mount 一次）────────────────────────────────────
  useEffect(() => {
    if (!containerRef.current) return
    let active = true

    const t = Date.now();
    fetch(`/api/roi/dzi?t=${t}`)
      .then(res => res.text())
      .then(xmlStr => {
        if (!active || !containerRef.current) return;
        const parser = new DOMParser();
        const doc = parser.parseFromString(xmlStr, "application/xml");
        const image = doc.getElementsByTagName("Image")[0];
        const size = doc.getElementsByTagName("Size")[0];

        if (!image || !size) {
            console.error("Failed to parse DZI XML", xmlStr);
            return;
        }

        const width = Number(size.getAttribute("Width"));
        const height = Number(size.getAttribute("Height"));
        const tileSize = Number(image.getAttribute("TileSize")) || 256;
        const overlap = Number(image.getAttribute("Overlap")) || 1;
        const format = image.getAttribute("Format") || "jpeg";

        const viewer = new OpenSeadragon.Viewer({
          element: containerRef.current,
          showNavigationControl: false,
          showNavigator: true,
          navigatorPosition: 'BOTTOM_RIGHT',
          gestureSettingsMouse: {
            clickToZoom:    false,
            dblClickToZoom: true,
            scrollToZoom:   true,
          },
          tileSources: {
            Image: {
              xmlns: "http://schemas.microsoft.com/deepzoom/2008",
              Url: `/api/roi/dzi_files/`,
              Format: format,
              Overlap: overlap.toString(),
              TileSize: tileSize.toString(),
              Size: {
                Width: width.toString(),
                Height: height.toString()
              }
            }
          }
        })

        // Monkey-patch getTileUrl to inject the cache-buster query parameter
        const patchTiledImage = () => {
          const tiledImage = viewer.world.getItemAt(0);
          if (tiledImage && tiledImage.source && !(tiledImage.source as any)._patched) {
            const originalGetTileUrl = tiledImage.source.getTileUrl.bind(tiledImage.source);
            tiledImage.source.getTileUrl = function(level: number, x: number, y: number) {
              const url = originalGetTileUrl(level, x, y);
              return `${url}?t=${t}`;
            };
            (tiledImage.source as any)._patched = true;
          }
        };

        if (viewer.world.getItemCount() > 0) {
          patchTiledImage();
        } else {
          viewer.world.addHandler('add-item', patchTiledImage);
        }

        viewer.addHandler('open', () => {
           if (active) setReady(true);
        });

        viewerRef.current = viewer;
      })
      .catch(e => console.error("Failed to fetch DZI:", e));

    return () => {
      active = false;
      if (viewerRef.current) {
        viewerRef.current.destroy()
        viewerRef.current = null
      }
      setReady(false)
    }
  }, [])

  // ── 重繪既有 ROI 覆蓋框 ────────────────────────────────────────
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return

    viewer.clearOverlays()
    for (const roi of existingRois) {
      if (roi.x == null || roi.width_px == null) continue

      const el = document.createElement('div')
      el.style.border        = '2px solid #61dafb'
      el.style.background    = 'rgba(97, 218, 251, 0.08)'
      el.style.pointerEvents = 'none'

      const label = document.createElement('span')
      label.textContent = roi.name ?? ''
      label.style.cssText = [
        'position:absolute', 'top:2px', 'left:4px',
        'font-size:10px', 'color:#61dafb', 'font-weight:600',
        'text-shadow:0 0 4px #000', 'white-space:nowrap',
      ].join(';')
      el.appendChild(label)

      const vpRect = viewer.viewport.imageToViewportRectangle(
        new OpenSeadragon.Rect(roi.x, roi.y, roi.width_px, roi.height_px)
      )
      viewer.addOverlay({ element: el, location: vpRect })
    }
  }, [existingRois, ready])

  // ── Pointer handlers（Draw mode 專用）──────────────────────────
  const handlePointerDown = (e: React.PointerEvent<HTMLDivElement>) => {
    const rect = containerRef.current!.getBoundingClientRect()
    const cx = e.clientX - rect.left
    const cy = e.clientY - rect.top
    startPt.current = { x: cx, y: cy }
    setDrawBox({ x: cx, y: cy, w: 0, h: 0 })
    setDrawing(true)
    ;(e.target as HTMLElement).setPointerCapture(e.pointerId)
  }

  const handlePointerMove = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!drawing || !startPt.current) return
    const rect = containerRef.current!.getBoundingClientRect()
    const cx   = e.clientX - rect.left
    const cy   = e.clientY - rect.top
    const dx   = cx - startPt.current.x
    const dy   = cy - startPt.current.y
    setDrawBox({
      x: dx >= 0 ? startPt.current.x : cx,
      y: dy >= 0 ? startPt.current.y : cy,
      w: Math.abs(dx),
      h: Math.abs(dy),
    })
  }

  const handlePointerUp = (e: React.PointerEvent<HTMLDivElement>) => {
    if (!drawing || !drawBox || !startPt.current) return
    ;(e.target as HTMLElement).releasePointerCapture(e.pointerId)
    setDrawing(false)

    if (drawBox.w > 10 && drawBox.h > 10 && viewerRef.current) {
      const viewer = viewerRef.current
      const toImg  = (cx: number, cy: number) =>
        viewer.viewport.viewportToImageCoordinates(
          viewer.viewport.pointFromPixel(new OpenSeadragon.Point(cx, cy))
        )

      const tl = toImg(drawBox.x, drawBox.y)
      const br = toImg(drawBox.x + drawBox.w, drawBox.y + drawBox.h)

      onSelect({
        x:         Math.round(Math.min(tl.x, br.x)),
        y:         Math.round(Math.min(tl.y, br.y)),
        width_px:  Math.round(Math.abs(br.x - tl.x)),
        height_px: Math.round(Math.abs(br.y - tl.y)),
      })
    }
    setDrawBox(null)
  }

  // ── Render ─────────────────────────────────────────────────────
  return (
    <div className="space-y-2">
      {/* 工具列 */}
      <div className="flex items-center gap-3">
        <button
          onClick={() => setMode(m => (m === 'pan' ? 'draw' : 'pan'))}
          className={clsx(
            'px-3 py-1 rounded text-xs font-medium transition-colors',
            mode === 'draw'
              ? 'bg-primary text-black'
              : 'bg-surface-border text-gray-300 hover:bg-surface-border/80',
          )}
        >
          {mode === 'pan' ? '✋ Pan 模式' : '✏️ 畫 ROI 模式'}
        </button>
        <span className="text-xs text-gray-500">
          {mode === 'pan'
            ? '拖曳平移 · 滾輪縮放 · 雙擊放大'
            : '按住拖曳框選 ROI 範圍，放開後座標自動填入表單'}
        </span>
        {!ready && (
          <span className="text-xs text-yellow-400 animate-pulse">載入影像中...</span>
        )}
      </div>

      {/* Viewer + Draw overlay */}
      <div className="relative rounded-lg overflow-hidden border border-surface-border">
        <div
          ref={containerRef}
          className="w-full bg-black"
          style={{ height: '26rem' }}
        />

        {/* Draw 模式：覆蓋層攔截 pointer events，阻止 OSD pan */}
        {mode === 'draw' && (
          <div
            className="absolute inset-0 cursor-crosshair"
            style={{ touchAction: 'none' }}
            onPointerDown={handlePointerDown}
            onPointerMove={handlePointerMove}
            onPointerUp={handlePointerUp}
            onPointerCancel={handlePointerUp}
          >
            {drawing && drawBox && (
              <div
                className="absolute border-2 border-primary bg-primary/10 pointer-events-none"
                style={{
                  left:   drawBox.x,
                  top:    drawBox.y,
                  width:  drawBox.w,
                  height: drawBox.h,
                }}
              >
                <span className="absolute -top-5 left-0 bg-primary text-black text-[10px] px-1 rounded-t font-mono whitespace-nowrap">
                  {drawBox.w} × {drawBox.h} px
                </span>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  )
}
