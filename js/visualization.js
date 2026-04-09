'use strict';

/**
 * MCseg – Canvas Visualization Layer
 *
 * Manages rendering of the source image, segmentation overlays, and borders.
 * Supports interactive zoom (mouse-wheel + pinch) and pan (drag).
 */
class Visualizer {
  /**
   * @param {HTMLCanvasElement} canvas
   */
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx    = canvas.getContext('2d');

    // Viewport state
    this.zoom = 1;
    this.panX = 0;
    this.panY = 0;

    // View mode: 'original' | 'overlay' | 'borders' | 'mask'
    this.mode = 'original';

    /** @type {ImageData|null} */
    this._srcImgData  = null;
    /** @type {ImageData|null} */
    this._overlayData = null;
    /** @type {ImageData|null} */
    this._borderData  = null;

    // Cached offscreen canvases (rebuilt when data changes)
    this._srcOC     = null;
    this._overlayOC = null;
    this._borderOC  = null;

    this._dragging = false;
    this._lastX    = 0;
    this._lastY    = 0;

    this._bindEvents();
  }

  // ─── Public ────────────────────────────────────────────────────────────────

  /** Load or replace the source H&E image. */
  setImage(imageData) {
    this._srcImgData = imageData;
    this._srcOC      = this._makeOffscreen(imageData);
    this.fitToCanvas();
    this.render();
  }

  /** Load segmentation result from an MCsegEngine instance. */
  setSegmentation(engine) {
    const w = engine.width, h = engine.height;
    this._overlayData = new ImageData(engine.renderOverlay(0.65), w, h);
    this._borderData  = new ImageData(engine.renderBorders(),     w, h);
    this._overlayOC   = this._makeOffscreen(this._overlayData);
    this._borderOC    = this._makeOffscreen(this._borderData);
    this.render();
  }

  /** Switch view mode and re-render. */
  setMode(mode) {
    this.mode = mode;
    this.render();
  }

  /** Fit the image to fill the canvas with 5 % margin. */
  fitToCanvas() {
    if (!this._srcImgData) return;
    const cw = this.canvas.width, ch = this.canvas.height;
    const iw = this._srcImgData.width, ih = this._srcImgData.height;
    this.zoom = Math.min(cw / iw, ch / ih) * 0.95;
    this.panX = (cw - iw * this.zoom) / 2;
    this.panY = (ch - ih * this.zoom) / 2;
  }

  /**
   * Zoom in/out around a focal point (canvas-space coordinates).
   * @param {number} factor  Zoom multiplier (> 1 = zoom in)
   * @param {number} cx      Canvas-space focal X
   * @param {number} cy      Canvas-space focal Y
   */
  zoomBy(factor, cx, cy) {
    const newZoom = Math.max(0.05, Math.min(40, this.zoom * factor));
    const ratio   = newZoom / this.zoom;
    this.panX     = cx - ratio * (cx - this.panX);
    this.panY     = cy - ratio * (cy - this.panY);
    this.zoom     = newZoom;
    this.render();
    return newZoom;
  }

  /** Full render pass. */
  render() {
    const { ctx, canvas, zoom, panX, panY } = this;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (!this._srcOC) return;

    const iw = this._srcImgData.width;
    const ih = this._srcImgData.height;
    const dw = iw * zoom;
    const dh = ih * zoom;

    // Smooth when zoomed out, pixelated when zoomed in (>= 4×)
    ctx.imageSmoothingEnabled = zoom < 4;
    ctx.imageSmoothingQuality = 'high';

    switch (this.mode) {
      case 'original':
        ctx.drawImage(this._srcOC, panX, panY, dw, dh);
        break;

      case 'overlay':
        ctx.drawImage(this._srcOC, panX, panY, dw, dh);
        if (this._overlayOC) ctx.drawImage(this._overlayOC, panX, panY, dw, dh);
        break;

      case 'borders':
        ctx.drawImage(this._srcOC, panX, panY, dw, dh);
        if (this._borderOC) ctx.drawImage(this._borderOC, panX, panY, dw, dh);
        break;

      case 'mask':
        // Black background + full-opacity mask
        ctx.fillStyle = '#111';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        if (this._overlayOC) {
          // Draw mask at full opacity by drawing a temporary copy with alpha = 1
          ctx.drawImage(this._overlayOC, panX, panY, dw, dh);
        }
        break;
    }
  }

  /**
   * Convert canvas-space (x, y) to image-space coordinates.
   * @returns {{ ix: number, iy: number }}
   */
  canvasToImage(cx, cy) {
    return {
      ix: Math.round((cx - this.panX) / this.zoom),
      iy: Math.round((cy - this.panY) / this.zoom),
    };
  }

  // ─── Private ───────────────────────────────────────────────────────────────

  /** Build an HTMLCanvasElement pre-painted with the given ImageData. */
  _makeOffscreen(imgData) {
    const oc  = document.createElement('canvas');
    oc.width  = imgData.width;
    oc.height = imgData.height;
    oc.getContext('2d').putImageData(imgData, 0, 0);
    return oc;
  }

  _bindEvents() {
    const c = this.canvas;

    // ── Wheel zoom ─────────────────────────────────────────────────────────
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const rect = c.getBoundingClientRect();
      const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
      const newZoom = this.zoomBy(factor, e.clientX - rect.left, e.clientY - rect.top);
      this._emitZoom(newZoom);
    }, { passive: false });

    // ── Mouse drag pan ─────────────────────────────────────────────────────
    c.addEventListener('mousedown', (e) => {
      this._dragging = true;
      this._lastX    = e.clientX;
      this._lastY    = e.clientY;
      c.style.cursor = 'grabbing';
    });
    window.addEventListener('mousemove', (e) => {
      // Report cursor position
      const rect = c.getBoundingClientRect();
      const { ix, iy } = this.canvasToImage(e.clientX - rect.left, e.clientY - rect.top);
      this._emitCursor(ix, iy);

      if (!this._dragging) return;
      this.panX += e.clientX - this._lastX;
      this.panY += e.clientY - this._lastY;
      this._lastX = e.clientX;
      this._lastY = e.clientY;
      this.render();
    });
    window.addEventListener('mouseup', () => {
      this._dragging = false;
      c.style.cursor = 'crosshair';
    });

    // ── Pinch zoom ─────────────────────────────────────────────────────────
    let lastPinchDist = 0;
    c.addEventListener('touchstart', (e) => {
      if (e.touches.length === 2) {
        lastPinchDist = Math.hypot(
          e.touches[0].clientX - e.touches[1].clientX,
          e.touches[0].clientY - e.touches[1].clientY
        );
      }
    }, { passive: true });

    c.addEventListener('touchmove', (e) => {
      if (e.touches.length !== 2) return;
      e.preventDefault();
      const dist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY
      );
      const rect = c.getBoundingClientRect();
      const cx   = (e.touches[0].clientX + e.touches[1].clientX) / 2 - rect.left;
      const cy   = (e.touches[0].clientY + e.touches[1].clientY) / 2 - rect.top;
      if (lastPinchDist) {
        const newZoom = this.zoomBy(dist / lastPinchDist, cx, cy);
        this._emitZoom(newZoom);
      }
      lastPinchDist = dist;
    }, { passive: false });
  }

  _emitZoom(zoom) {
    const pct = Math.round(zoom * 100);
    document.dispatchEvent(new CustomEvent('mcseg:zoom', { detail: { pct } }));
  }

  _emitCursor(ix, iy) {
    document.dispatchEvent(new CustomEvent('mcseg:cursor', { detail: { ix, iy } }));
  }
}
