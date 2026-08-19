import RegionSelector, { type BoxOverlay } from '../shared/RegionSelector'
import type { RnaBounds } from '../../api/client'

/**
 * Stage 0 的 ROI 畫框器。
 *
 * 畫框邏輯已提取至 `shared/RegionSelector`（Stage 0 與 Stage 3.5 共用）；
 * 本元件僅保留 Stage 0 既有的 `onSelect` 介面，行為不變。
 */

interface Props {
  onSelect: (roi: Omit<BoxOverlay, 'name'>) => void
  existingRois?: BoxOverlay[]
  captureBounds?: RnaBounds | null
}

export default function RoiSelector({ onSelect, existingRois = [], captureBounds }: Props) {
  return (
    <RegionSelector
      mode="bbox"
      overlays={existingRois}
      captureBounds={captureBounds}
      onChange={sel => {
        if (sel?.type !== 'bbox') return
        onSelect({ x: sel.x, y: sel.y, width_px: sel.width_px, height_px: sel.height_px })
      }}
    />
  )
}

