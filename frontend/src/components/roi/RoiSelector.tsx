import React, { useState, useRef, useEffect } from 'react'

interface RoiSelectorProps {
    imageB64: string
    widthHires: number
    heightHires: number
    scalef: number
    mpp: number
    onSelect: (roi: { x: number; y: number; width_px: number; height_px: number }) => void
}

export default function RoiSelector({ imageB64, widthHires, heightHires, scalef, mpp, onSelect }: RoiSelectorProps) {
    const containerRef = useRef<HTMLDivElement>(null)
    const imageRef = useRef<HTMLImageElement>(null)

    const [isDrawing, setIsDrawing] = useState(false)
    const [startPos, setStartPos] = useState({ x: 0, y: 0 })
    const [currentBox, setCurrentBox] = useState<{ x: number, y: number, w: number, h: number } | null>(null)

    // 顯示比例 (圖片實際顯示寬度 vs hires 原圖寬度)
    const [displayScale, setDisplayScale] = useState(1)

    useEffect(() => {
        const updateScale = () => {
            if (imageRef.current) {
                setDisplayScale(imageRef.current.clientWidth / widthHires)
            }
        }
        window.addEventListener('resize', updateScale)
        // 稍微延遲以確保圖片載入完畢
        setTimeout(updateScale, 100)
        return () => window.removeEventListener('resize', updateScale)
    }, [widthHires, imageB64])

    const handlePointerDown = (e: React.PointerEvent) => {
        if (!imageRef.current) return
        const rect = imageRef.current.getBoundingClientRect()
        const x = e.clientX - rect.left
        const y = e.clientY - rect.top
        setStartPos({ x, y })
        setCurrentBox({ x, y, w: 0, h: 0 })
        setIsDrawing(true)
            ; (e.target as HTMLElement).setPointerCapture(e.pointerId)
    }

    const handlePointerMove = (e: React.PointerEvent) => {
        if (!isDrawing || !imageRef.current) return
        const rect = imageRef.current.getBoundingClientRect()
        const currentX = Math.max(0, Math.min(e.clientX - rect.left, rect.width))
        const currentY = Math.max(0, Math.min(e.clientY - rect.top, rect.height))

        setCurrentBox({
            x: Math.min(startPos.x, currentX),
            y: Math.min(startPos.y, currentY),
            w: Math.abs(currentX - startPos.x),
            h: Math.abs(currentY - startPos.y)
        })
    }

    const handlePointerUp = (e: React.PointerEvent) => {
        if (!isDrawing) return
        setIsDrawing(false)
            ; (e.target as HTMLElement).releasePointerCapture(e.pointerId)

        if (currentBox && currentBox.w > 10 && currentBox.h > 10) {
            // 轉換座標：畫面 px -> hires px -> fullres px
            const hiresX = currentBox.x / displayScale
            const hiresY = currentBox.y / displayScale
            const hiresW = currentBox.w / displayScale
            const hiresH = currentBox.h / displayScale

            onSelect({
                x: Math.round(hiresX / scalef),
                y: Math.round(hiresY / scalef),
                width_px: Math.round(hiresW / scalef),
                height_px: Math.round(hiresH / scalef)
            })
        }
    }

    // 計算全圖大約 mm
    const fullresW = widthHires / scalef
    const fullresH = heightHires / scalef
    const mmW = (fullresW * mpp) / 1000
    const mmH = (fullresH * mpp) / 1000

    return (
        <div className="flex flex-col gap-2">
            <div className="text-xs text-gray-400 bg-surface/50 p-2 rounded">
                <p>H&E Hires: {widthHires} × {heightHires} px | 轉換比例: {scalef.toFixed(5)}</p>
                <p>組織尺寸估計: ~{mmW.toFixed(1)} × {mmH.toFixed(1)} mm</p>
                <p className="text-primary mt-1">💡 請在下方影像按住拖曳，畫出要裁切的範圍 (ROI)</p>
            </div>

            <div
                ref={containerRef}
                className="relative border border-primary/20 rounded-lg overflow-hidden bg-black/40 w-full max-w-2xl select-none"
                style={{ touchAction: 'none' }}
            >
                <img
                    ref={imageRef}
                    src={`data:image/jpeg;base64,${imageB64}`}
                    alt="H&E Overview"
                    className="w-full h-auto block select-none pointer-events-none"
                    onLoad={() => {
                        if (imageRef.current) setDisplayScale(imageRef.current.clientWidth / widthHires)
                    }}
                />

                {/* Drawing Overlay */}
                <div
                    className="absolute inset-0 cursor-crosshair"
                    onPointerDown={handlePointerDown}
                    onPointerMove={handlePointerMove}
                    onPointerUp={handlePointerUp}
                    onPointerCancel={handlePointerUp}
                >
                    {currentBox && (
                        <div
                            className="absolute border-2 border-primary bg-primary/10"
                            style={{
                                left: currentBox.x,
                                top: currentBox.y,
                                width: currentBox.w,
                                height: currentBox.h
                            }}
                        >
                            <div className="absolute -top-6 left-0 bg-primary text-black text-[10px] px-1 font-bold whitespace-nowrap rounded-t opacity-80">
                                {Math.round(currentBox.w / displayScale / scalef)}x{Math.round(currentBox.h / displayScale / scalef)} px
                            </div>
                        </div>
                    )}
                </div>
            </div>
        </div>
    )
}
