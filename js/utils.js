'use strict';

/**
 * MCseg – Shared Utilities
 * Small helpers used by both app.js and segmentation.js.
 */

/**
 * Parse a Visium HD spatial-barcode CSV / TSV.
 * Recognised column names (case-insensitive):
 *   barcode / barcodes / barcode_id
 *   x / col / pxl_col_in_fullres
 *   y / row / pxl_row_in_fullres
 *
 * @param {string} text  Raw file text
 * @returns {{ id: string, x: number, y: number }[]}
 */
function parseSpatialCSV(text) {
  const sep   = text.includes('\t') ? '\t' : ',';
  const lines = text.trim().split(/\r?\n/);
  const header = lines[0].split(sep).map(h => h.trim().toLowerCase().replace(/['"]/g, ''));

  const colIdx = (candidates) => {
    for (const c of candidates) {
      const i = header.indexOf(c);
      if (i !== -1) return i;
    }
    return -1;
  };

  const xi = colIdx(['pxl_col_in_fullres', 'x', 'col', 'x_coord', 'xcoord']);
  const yi = colIdx(['pxl_row_in_fullres', 'y', 'row', 'y_coord', 'ycoord']);
  const bi = colIdx(['barcode', 'barcodes', 'barcode_id', 'id']);

  if (xi < 0 || yi < 0) {
    throw new Error(
      `Cannot find x/y columns in CSV. ` +
      `Header: [${header.join(', ')}]. ` +
      `Expected: x / pxl_col_in_fullres, y / pxl_row_in_fullres.`
    );
  }

  const barcodes = [];
  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;
    const cols = line.split(sep).map(c => c.replace(/['"]/g, '').trim());
    const x = parseFloat(cols[xi]);
    const y = parseFloat(cols[yi]);
    if (!isFinite(x) || !isFinite(y)) continue;
    barcodes.push({
      id: bi >= 0 ? cols[bi] : `bc_${i}`,
      x,
      y,
    });
  }
  return barcodes;
}

/**
 * Trigger a file download in the browser.
 * @param {string} url       Object URL or data-URL
 * @param {string} filename  Suggested filename
 */
function downloadURL(url, filename) {
  const a = document.createElement('a');
  a.href     = url;
  a.download = filename;
  a.style.display = 'none';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}

/**
 * Convert an HSL colour to an [r, g, b] triple (each 0–255).
 * @param {number} h  Hue, 0–360
 * @param {number} s  Saturation, 0–100
 * @param {number} l  Lightness, 0–100
 * @returns {[number, number, number]}
 */
function hslToRgb(h, s, l) {
  h /= 360; s /= 100; l /= 100;
  const a = s * Math.min(l, 1 - l);
  const f = (n) => {
    const k = (n + h * 12) % 12;
    return Math.round((l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1))) * 255);
  };
  return [f(0), f(8), f(4)];
}

/**
 * Generate N visually distinct colours using the golden-angle hue step.
 * @param {number} n
 * @returns {[number, number, number][]}
 */
function generatePalette(n) {
  const palette = [];
  for (let i = 1; i <= n; i++) {
    palette.push(hslToRgb((i * 137.508) % 360, 72, 55));
  }
  return palette;
}
