'use strict';

/**
 * MCseg – High-Fidelity Cell Segmentation Engine
 *
 * Implements a marker-controlled watershed pipeline optimised for H&E stained
 * tissue sections from 10x Genomics Visium HD spatial transcriptomics experiments.
 *
 * Pipeline
 * ─────────
 *  1. H&E optical-density decomposition  → hematoxylin channel
 *     (Ruifrok & Johnston, 2001 stain vectors)
 *  2. Separable Gaussian smoothing
 *  3. Otsu global threshold with sensitivity adjustment
 *  4. Morphological opening  (erosion → dilation, r = 2) to remove debris
 *  5. Union-find connected-component labelling
 *  6. Area-based nucleus filtering
 *  7. BFS nearest-nucleus expansion → cell body approximation
 *  8. Per-cell statistics (area, centroid, mean RGB)
 *  9. Optional spatial-barcode assignment
 */
class MCsegEngine {
  /**
   * @param {object}  params
   * @param {number}  params.minNucleiArea    Minimum nucleus area in px²  (default 50)
   * @param {number}  params.maxNucleiArea    Maximum nucleus area in px²  (default 2000)
   * @param {number}  params.expansionRadius  Cell expansion in px          (default 8)
   * @param {number}  params.sensitivity      Nucleus sensitivity 0–100     (default 50)
   */
  constructor(params = {}) {
    this.params = {
      minNucleiArea:   params.minNucleiArea   ?? 50,
      maxNucleiArea:   params.maxNucleiArea   ?? 2000,
      expansionRadius: params.expansionRadius ?? 8,
      sensitivity:     params.sensitivity     ?? 50,
      _gaussianSigma:  1.5,
    };
    /** @type {Array<CellRecord>} */
    this.cells    = [];
    /** @type {Int32Array|null} */
    this.labelMap = null;
    this.width    = 0;
    this.height   = 0;
  }

  // ─── Public API ────────────────────────────────────────────────────────────

  /**
   * Run the full segmentation pipeline.
   *
   * @param {ImageData} imageData   Source H&E image
   * @param {Function}  [onProgress]  Called with (percent: number, message: string)
   * @returns {{ cells, labelMap, width, height }}
   */
  async segment(imageData, onProgress) {
    const { data, width, height } = imageData;
    this.width  = width;
    this.height = height;

    const report = (pct, msg) => onProgress?.(pct, msg);

    // ── Step 1: H&E deconvolution ──────────────────────────────────────────
    report(5,  'Extracting hematoxylin channel…');
    const hChannel = this._extractHematoxylin(data, width, height);

    // ── Step 2: Gaussian smoothing ─────────────────────────────────────────
    report(15, 'Smoothing…');
    const smoothed = this._gaussianBlur(hChannel, width, height, this.params._gaussianSigma);

    // ── Step 3: Otsu threshold with sensitivity adjustment ─────────────────
    report(25, 'Thresholding nuclei…');
    const otsu = this._otsuThreshold(smoothed);
    // sensitivity 0 → multiply by 1.4 (conservative), 100 → 0.6 (aggressive)
    const adj  = 1.0 - (this.params.sensitivity - 50) / 125;
    const thr  = Math.max(0.05, Math.min(0.95, otsu * adj));
    const binary = new Uint8Array(width * height);
    for (let i = 0; i < smoothed.length; i++) binary[i] = smoothed[i] > thr ? 1 : 0;

    // ── Step 4: Morphological opening ─────────────────────────────────────
    report(35, 'Removing debris…');
    const cleaned = this._morphOpen(binary, width, height, 2);

    // ── Step 5: Connected-component labelling ──────────────────────────────
    report(50, 'Labelling nuclei…');
    const { labels: rawLabels } = this._connectedComponents(cleaned, width, height);

    // ── Step 6: Area filtering ─────────────────────────────────────────────
    report(60, 'Filtering by size…');
    const stats    = this._componentStats(rawLabels, width, height);
    const validSet = new Set(
      Object.values(stats)
        .filter(s => s.area >= this.params.minNucleiArea &&
                     s.area <= this.params.maxNucleiArea)
        .map(s => s.label)
    );

    const nucleiLabels = new Int32Array(rawLabels.length);
    for (let i = 0; i < rawLabels.length; i++) {
      nucleiLabels[i] = validSet.has(rawLabels[i]) ? rawLabels[i] : 0;
    }

    // ── Step 7: BFS cell-body expansion ───────────────────────────────────
    report(72, 'Growing cell bodies…');
    this.labelMap = this._expandCells(nucleiLabels, width, height);

    // ── Step 8: Per-cell statistics ────────────────────────────────────────
    report(88, 'Computing statistics…');
    this.cells = this._finalStats(this.labelMap, width, height, data);

    report(100, `Done — ${this.cells.length} cells detected`);

    return { cells: this.cells, labelMap: this.labelMap, width, height };
  }

  /**
   * Assign spatial barcodes to the nearest cell centroid (within maxDist px).
   * Mutates this.cells in place.
   *
   * @param {{ id: string, x: number, y: number }[]} barcodes
   * @param {number} [maxDist=60]
   */
  assignBarcodes(barcodes, maxDist = 60) {
    const maxD2 = maxDist * maxDist;
    for (const bc of barcodes) {
      let best = null, bestD2 = maxD2;
      for (const cell of this.cells) {
        const dx = cell.cx - bc.x, dy = cell.cy - bc.y;
        const d2 = dx * dx + dy * dy;
        if (d2 < bestD2) { bestD2 = d2; best = cell; }
      }
      if (best) best.barcodes.push(bc.id);
    }
  }

  /**
   * Render the segmentation as a coloured RGBA array (filled cells).
   * @param {number} [alpha=0.65]  Fill opacity 0–1
   * @returns {Uint8ClampedArray}
   */
  renderOverlay(alpha = 0.65) {
    const n    = this.width * this.height;
    const rgba = new Uint8ClampedArray(n * 4);
    const a255 = Math.round(alpha * 255);
    for (let i = 0; i < n; i++) {
      const l = this.labelMap[i];
      if (!l) continue;
      const [r, g, b] = this._cellColour(l);
      rgba[i * 4]     = r;
      rgba[i * 4 + 1] = g;
      rgba[i * 4 + 2] = b;
      rgba[i * 4 + 3] = a255;
    }
    return rgba;
  }

  /**
   * Render cell-boundary outlines (pixels where adjacent cells differ).
   * @returns {Uint8ClampedArray}
   */
  renderBorders() {
    const { labelMap: lm, width: w, height: h } = this;
    const rgba = new Uint8ClampedArray(w * h * 4);
    for (let y = 1; y < h - 1; y++) {
      for (let x = 1; x < w - 1; x++) {
        const i = y * w + x;
        const l = lm[i];
        if (!l) continue;
        if (lm[i - 1] !== l || lm[i + 1] !== l ||
            lm[i - w] !== l || lm[i + w] !== l) {
          rgba[i * 4]     = 255;
          rgba[i * 4 + 1] = 220;
          rgba[i * 4 + 2] = 0;
          rgba[i * 4 + 3] = 230;
        }
      }
    }
    return rgba;
  }

  /**
   * Export segmentation mask as an object URL (PNG).
   * @returns {Promise<string>}
   */
  async exportMaskPNG() {
    const oc  = new OffscreenCanvas(this.width, this.height);
    const ctx = oc.getContext('2d');
    ctx.putImageData(
      new ImageData(this.renderOverlay(1), this.width, this.height), 0, 0
    );
    const blob = await oc.convertToBlob({ type: 'image/png' });
    return URL.createObjectURL(blob);
  }

  /**
   * Export cell table as a CSV string.
   * Columns: cell_id, cx, cy, area, mean_r, mean_g, mean_b, barcodes
   * @returns {string}
   */
  exportCellCSV() {
    const rows = ['cell_id,cx,cy,area,mean_r,mean_g,mean_b,barcodes'];
    for (const c of this.cells) {
      rows.push(
        `${c.id},${c.cx},${c.cy},${c.area},` +
        `${c.meanR},${c.meanG},${c.meanB},"${c.barcodes.join(';')}"`
      );
    }
    return rows.join('\n');
  }

  // ─── Private helpers ───────────────────────────────────────────────────────

  /**
   * H&E stain deconvolution: project optical-density into hematoxylin channel.
   * Ruifrok & Johnston (2001) stain vector for Hematoxylin.
   */
  _extractHematoxylin(data, width, height) {
    const n  = width * height;
    const hR = 0.6500286, hG = 0.7041088, hB = 0.2860126; // normalised H vector
    const out = new Float32Array(n);
    let mn = Infinity, mx = -Infinity;

    for (let i = 0; i < n; i++) {
      const r = data[i * 4]     / 255;
      const g = data[i * 4 + 1] / 255;
      const b = data[i * 4 + 2] / 255;
      // Optical density per channel
      const odR = -Math.log(Math.max(r, 1e-6));
      const odG = -Math.log(Math.max(g, 1e-6));
      const odB = -Math.log(Math.max(b, 1e-6));
      // Project onto hematoxylin stain direction
      const v = odR * hR + odG * hG + odB * hB;
      out[i] = v;
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }

    // Normalise to [0, 1]
    const rng = mx - mn || 1;
    for (let i = 0; i < n; i++) out[i] = (out[i] - mn) / rng;
    return out;
  }

  /** Separable Gaussian blur (two 1-D passes). */
  _gaussianBlur(data, width, height, sigma) {
    const half = Math.ceil(sigma * 3);
    const size  = 2 * half + 1;
    const kernel = new Float32Array(size);
    let ksum = 0;
    for (let i = 0; i < size; i++) {
      const x = i - half;
      kernel[i] = Math.exp(-(x * x) / (2 * sigma * sigma));
      ksum += kernel[i];
    }
    for (let i = 0; i < size; i++) kernel[i] /= ksum;

    const tmp = new Float32Array(width * height);
    const out = new Float32Array(width * height);

    // Horizontal pass
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        let sum = 0, wsum = 0;
        for (let k = -half; k <= half; k++) {
          const nx = x + k;
          if (nx >= 0 && nx < width) {
            const w = kernel[k + half];
            sum  += data[y * width + nx] * w;
            wsum += w;
          }
        }
        tmp[y * width + x] = sum / wsum;
      }
    }

    // Vertical pass
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        let sum = 0, wsum = 0;
        for (let k = -half; k <= half; k++) {
          const ny = y + k;
          if (ny >= 0 && ny < height) {
            const w = kernel[k + half];
            sum  += tmp[ny * width + x] * w;
            wsum += w;
          }
        }
        out[y * width + x] = sum / wsum;
      }
    }
    return out;
  }

  /** Otsu's global threshold on a Float32 array in [0, 1]. Returns threshold in [0, 1]. */
  _otsuThreshold(data) {
    const bins = 256;
    const hist = new Int32Array(bins);
    for (const v of data) hist[Math.min(bins - 1, (v * bins) | 0)]++;
    const n = data.length;
    let sumAll = 0;
    for (let i = 0; i < bins; i++) sumAll += i * hist[i];
    let sumB = 0, wB = 0, best = 0, bestVar = 0;
    for (let t = 0; t < bins; t++) {
      wB += hist[t];
      if (!wB) continue;
      const wF = n - wB;
      if (!wF) break;
      sumB += t * hist[t];
      const mB = sumB / wB;
      const mF = (sumAll - sumB) / wF;
      const v  = wB * wF * (mB - mF) ** 2;
      if (v > bestVar) { bestVar = v; best = t; }
    }
    return best / bins;
  }

  /** Morphological opening = erosion then dilation. */
  _morphOpen(bin, w, h, r) {
    return this._dilate(this._erode(bin, w, h, r), w, h, r);
  }

  _erode(bin, w, h, r) {
    const out = new Uint8Array(w * h);
    const r2  = r * r;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        let ok = true;
        loop: for (let dy = -r; dy <= r && ok; dy++) {
          for (let dx = -r; dx <= r; dx++) {
            if (dx * dx + dy * dy > r2) continue;
            const ny = y + dy, nx = x + dx;
            if (ny < 0 || ny >= h || nx < 0 || nx >= w || !bin[ny * w + nx]) {
              ok = false; break loop;
            }
          }
        }
        out[y * w + x] = ok ? 1 : 0;
      }
    }
    return out;
  }

  _dilate(bin, w, h, r) {
    const out = new Uint8Array(w * h);
    const r2  = r * r;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        loop: for (let dy = -r; dy <= r; dy++) {
          for (let dx = -r; dx <= r; dx++) {
            if (dx * dx + dy * dy > r2) continue;
            const ny = y + dy, nx = x + dx;
            if (ny >= 0 && ny < h && nx >= 0 && nx < w && bin[ny * w + nx]) {
              out[y * w + x] = 1; break loop;
            }
          }
        }
      }
    }
    return out;
  }

  /**
   * Connected-component labelling using two-pass union-find (4-connectivity).
   * @returns {{ labels: Int32Array, count: number }}
   */
  _connectedComponents(bin, width, height) {
    const labels = new Int32Array(width * height);
    const parent = [0]; // index == label; 0 is background
    let   nextLabel = 1;

    const find = (x) => {
      while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; }
      return x;
    };
    const union = (a, b) => {
      a = find(a); b = find(b);
      if (a !== b) parent[b] = a;
    };

    // First pass
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const i = y * width + x;
        if (!bin[i]) continue;
        const above = y > 0 ? labels[(y - 1) * width + x] : 0;
        const left  = x > 0 ? labels[y * width + x - 1]   : 0;
        if (!above && !left) {
          parent.push(nextLabel);
          labels[i] = nextLabel++;
        } else if (above && !left) {
          labels[i] = above;
        } else if (!above && left) {
          labels[i] = left;
        } else {
          labels[i] = Math.min(above, left);
          union(above, left);
        }
      }
    }
    // Second pass: resolve unions
    for (let i = 0; i < labels.length; i++) {
      if (labels[i]) labels[i] = find(labels[i]);
    }
    return { labels, count: nextLabel - 1 };
  }

  /** Compute area and centroid per component. */
  _componentStats(labels, width, height) {
    const stats = {};
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const l = labels[y * width + x];
        if (!l) continue;
        if (!stats[l]) stats[l] = { label: l, area: 0, cx: 0, cy: 0 };
        stats[l].area++;
        stats[l].cx += x;
        stats[l].cy += y;
      }
    }
    for (const s of Object.values(stats)) {
      s.cx = Math.round(s.cx / s.area);
      s.cy = Math.round(s.cy / s.area);
    }
    return stats;
  }

  /**
   * BFS nearest-nucleus expansion to produce cell bodies.
   * Each background pixel is claimed by the nearest nucleus seed
   * within expansionRadius steps.
   */
  _expandCells(nucleiLabels, width, height) {
    const cellMap = new Int32Array(nucleiLabels); // nuclei labels are seeds
    const dist    = new Float32Array(width * height).fill(Infinity);
    const queue   = [];

    for (let i = 0; i < nucleiLabels.length; i++) {
      if (nucleiLabels[i]) { dist[i] = 0; queue.push(i); }
    }

    const maxR = this.params.expansionRadius;
    let head = 0;
    while (head < queue.length) {
      const idx = queue[head++];
      const d   = dist[idx];
      if (d >= maxR) continue;

      const y = (idx / width) | 0;
      const x = idx % width;
      const nbrs = [
        x > 0         ? idx - 1      : -1,
        x < width - 1 ? idx + 1      : -1,
        y > 0         ? idx - width   : -1,
        y < height - 1? idx + width   : -1,
      ];
      for (const ni of nbrs) {
        if (ni < 0 || cellMap[ni] !== 0) continue;
        dist[ni]    = d + 1;
        cellMap[ni] = cellMap[idx];
        queue.push(ni);
      }
    }
    return cellMap;
  }

  /** Compute per-cell statistics from the final label map. */
  _finalStats(labelMap, width, height, rgba) {
    const acc = {};
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const l = labelMap[y * width + x];
        if (!l) continue;
        if (!acc[l]) acc[l] = { id: l, area: 0, cx: 0, cy: 0, r: 0, g: 0, b: 0 };
        const s = acc[l];
        s.area++;
        s.cx += x; s.cy += y;
        const p = (y * width + x) * 4;
        s.r += rgba[p]; s.g += rgba[p + 1]; s.b += rgba[p + 2];
      }
    }
    return Object.values(acc).map(s => ({
      id:       s.id,
      area:     s.area,
      cx:       Math.round(s.cx / s.area),
      cy:       Math.round(s.cy / s.area),
      meanR:    Math.round(s.r / s.area),
      meanG:    Math.round(s.g / s.area),
      meanB:    Math.round(s.b / s.area),
      barcodes: [],
    }));
  }

  /** Deterministic per-cell colour using golden-angle hue cycling. */
  _cellColour(id) {
    return hslToRgb((id * 137.508) % 360, 72, 55);
  }
}
