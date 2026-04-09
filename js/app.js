'use strict';

/**
 * MCseg – Application Controller
 *
 * Orchestrates file loading, parameter controls, segmentation runs,
 * visualization updates, and data exports.
 */

// ─── State ───────────────────────────────────────────────────────────────────
const state = {
  /** @type {ImageData|null} */  imageData: null,
  /** @type {MCsegEngine|null} */ engine:    null,
  /** @type {Visualizer|null} */  viz:       null,
  /** @type {Array}            */ barcodes:  [],
};

// ─── DOM helpers ─────────────────────────────────────────────────────────────
const $ = (id) => document.getElementById(id);

// ─── Init ─────────────────────────────────────────────────────────────────────
window.addEventListener('DOMContentLoaded', () => {
  const canvas = $('mainCanvas');
  resizeCanvas(canvas);
  state.viz = new Visualizer(canvas);

  window.addEventListener('resize', () => {
    resizeCanvas(canvas);
    state.viz.fitToCanvas();
    state.viz.render();
  });

  document.addEventListener('mcseg:cursor', (e) => {
    const { ix, iy } = e.detail;
    $('cursorPos').textContent = `x: ${ix}  y: ${iy}`;
  });
  document.addEventListener('mcseg:zoom', (e) => {
    $('zoomLevelDisp').textContent = `${e.detail.pct}%`;
  });

  bindUploads();
  bindParams();
  bindRunButton();
  bindViewButtons();
  bindZoomButtons();
  bindExports();
  bindTheme();
});

function resizeCanvas(canvas) {
  const container = $('canvasContainer');
  canvas.width  = container.clientWidth;
  canvas.height = container.clientHeight;
}

// ─── File uploads ─────────────────────────────────────────────────────────────
function bindUploads() {
  // H&E image
  const imageZone = $('imageUpload');
  const imageFile = $('imageFile');
  imageZone.addEventListener('click', (e) => {
    // Prevent double-firing when the hidden <input> is clicked directly
    if (e.target !== imageFile) imageFile.click();
  });
  imageZone.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); imageFile.click(); }
  });
  imageFile.addEventListener('change', (e) => {
    if (e.target.files[0]) loadImageFile(e.target.files[0]);
  });
  makeDropZone(imageZone, (f) => {
    if (f.type.startsWith('image/') || /\.(tif|tiff)$/i.test(f.name)) loadImageFile(f);
  });

  // Spatial barcodes
  const spatialZone = $('spatialUpload');
  const spatialFile = $('spatialFile');
  spatialZone.addEventListener('click', (e) => {
    if (e.target !== spatialFile) spatialFile.click();
  });
  spatialZone.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); spatialFile.click(); }
  });
  spatialFile.addEventListener('change', (e) => {
    if (e.target.files[0]) loadSpatialFile(e.target.files[0]);
  });
  makeDropZone(spatialZone, (f) => {
    if (/\.(csv|tsv)$/i.test(f.name)) loadSpatialFile(f);
  });
}

function makeDropZone(el, onFile) {
  el.addEventListener('dragenter', (e) => { e.preventDefault(); el.classList.add('drag-over'); });
  el.addEventListener('dragover',  (e) => { e.preventDefault(); el.classList.add('drag-over'); });
  el.addEventListener('dragleave', (e) => {
    // Only remove if actually leaving the zone (not entering a child)
    if (!el.contains(e.relatedTarget)) el.classList.remove('drag-over');
  });
  el.addEventListener('drop', (e) => {
    e.preventDefault();
    el.classList.remove('drag-over');
    const f = e.dataTransfer.files[0];
    if (f) onFile(f);
  });
}

async function loadImageFile(file) {
  setStatus(`Loading ${file.name}…`);
  try {
    // createImageBitmap handles PNG, JPEG, WebP; TIFF requires browser support
    const bitmap  = await createImageBitmap(file);
    const { width, height } = bitmap;

    // Warn / auto-downsample images larger than 3000×3000
    const MAX_SIDE = 3000;
    let targetW = width, targetH = height;
    if (width > MAX_SIDE || height > MAX_SIDE) {
      const scale = MAX_SIDE / Math.max(width, height);
      targetW = Math.round(width  * scale);
      targetH = Math.round(height * scale);
      setStatus(`Image too large (${width}×${height}); downsampling to ${targetW}×${targetH}…`);
    }

    const oc  = new OffscreenCanvas(targetW, targetH);
    const ctx = oc.getContext('2d');
    ctx.drawImage(bitmap, 0, 0, targetW, targetH);
    bitmap.close();

    state.imageData = ctx.getImageData(0, 0, targetW, targetH);

    // Show image
    $('emptyState').classList.add('hidden');
    state.viz.setImage(state.imageData);

    // Update upload zone label
    $('imageLabel').textContent = file.name;
    $('imageUpload').classList.add('loaded');

    $('runBtn').disabled = false;
    setStatus(`Image loaded: ${targetW}×${targetH} px`);
    $('zoomLevelDisp').textContent = `${Math.round(state.viz.zoom * 100)}%`;
  } catch (err) {
    setStatus(`Error loading image: ${err.message}`);
    console.error('[MCseg] Image load error:', err);
  }
}

function loadSpatialFile(file) {
  setStatus(`Loading ${file.name}…`);
  const reader = new FileReader();
  reader.onload = (e) => {
    try {
      state.barcodes = parseSpatialCSV(e.target.result);
      $('spatialLabel').textContent = file.name;
      $('spatialUpload').classList.add('loaded');
      setStatus(`Loaded ${state.barcodes.length.toLocaleString()} barcodes from ${file.name}`);
    } catch (err) {
      setStatus(`CSV parse error: ${err.message}`);
      console.error('[MCseg] CSV parse error:', err);
    }
  };
  reader.onerror = () => setStatus('Error reading spatial file.');
  reader.readAsText(file);
}

// ─── Parameter sliders ────────────────────────────────────────────────────────
function bindParams() {
  const pairs = [
    ['minArea',    'minAreaVal'],
    ['maxArea',    'maxAreaVal'],
    ['expansion',  'expansionVal'],
    ['sensitivity','sensitivityVal'],
  ];
  for (const [slId, outId] of pairs) {
    $(slId).addEventListener('input', () => {
      $(outId).value = $(slId).value;
    });
  }
}

// ─── Run segmentation ─────────────────────────────────────────────────────────
function bindRunButton() {
  $('runBtn').addEventListener('click', runSegmentation);
}

async function runSegmentation() {
  if (!state.imageData) return;

  $('runBtn').disabled = true;
  showProgress(0, 'Starting…');

  const params = {
    minNucleiArea:   +$('minArea').value,
    maxNucleiArea:   +$('maxArea').value,
    expansionRadius: +$('expansion').value,
    sensitivity:     +$('sensitivity').value,
  };

  state.engine = new MCsegEngine(params);

  try {
    const result = await runAsync(() =>
      state.engine.segment(state.imageData, (pct, msg) => {
        showProgress(pct, msg);
        updateProgressAttr(pct);
      })
    );

    // Barcode assignment
    if (state.barcodes.length > 0) {
      state.engine.assignBarcodes(state.barcodes);
    }

    // Update visualizer
    state.viz.setSegmentation(state.engine);
    state.viz.setMode('overlay');
    setActiveViewBtn('viewOverlay');

    // Update stats panel
    const totalBarcodes = result.cells.reduce((s, c) => s + c.barcodes.length, 0);
    const avgArea = result.cells.length
      ? Math.round(result.cells.reduce((s, c) => s + c.area, 0) / result.cells.length)
      : 0;
    $('cellCount').textContent   = result.cells.length.toLocaleString();
    $('avgArea').textContent     = avgArea.toLocaleString();
    $('barcodeCount').textContent = totalBarcodes.toLocaleString();

    $('exportMaskBtn').disabled = false;
    $('exportCsvBtn').disabled  = false;

    setStatus(`Segmentation complete — ${result.cells.length.toLocaleString()} cells`);
  } catch (err) {
    setStatus(`Segmentation error: ${err.message}`);
    console.error('[MCseg] Segmentation error:', err);
  } finally {
    hideProgress();
    $('runBtn').disabled = false;
  }
}

/**
 * Wrap an async function with a minimal yield so the browser can repaint
 * before the synchronous segmentation loop begins.
 */
function runAsync(fn) {
  return new Promise((resolve, reject) => {
    setTimeout(() => fn().then(resolve).catch(reject), 20);
  });
}

// ─── View buttons ─────────────────────────────────────────────────────────────
function bindViewButtons() {
  const modes = {
    viewOriginal: 'original',
    viewOverlay:  'overlay',
    viewBorders:  'borders',
    viewMask:     'mask',
  };
  for (const [id, mode] of Object.entries(modes)) {
    $(id)?.addEventListener('click', () => {
      state.viz?.setMode(mode);
      setActiveViewBtn(id);
    });
  }
}

function setActiveViewBtn(activeId) {
  for (const id of ['viewOriginal', 'viewOverlay', 'viewBorders', 'viewMask']) {
    $(id)?.classList.toggle('active', id === activeId);
  }
}

// ─── Zoom buttons ─────────────────────────────────────────────────────────────
function bindZoomButtons() {
  const mid = () => [
    $('canvasContainer').clientWidth  / 2,
    $('canvasContainer').clientHeight / 2,
  ];
  $('zoomIn')?.addEventListener('click', () => {
    const [cx, cy] = mid();
    const z = state.viz?.zoomBy(1.25, cx, cy);
    if (z) $('zoomLevelDisp').textContent = `${Math.round(z * 100)}%`;
  });
  $('zoomOut')?.addEventListener('click', () => {
    const [cx, cy] = mid();
    const z = state.viz?.zoomBy(0.8, cx, cy);
    if (z) $('zoomLevelDisp').textContent = `${Math.round(z * 100)}%`;
  });
  $('zoomFit')?.addEventListener('click', () => {
    state.viz?.fitToCanvas();
    state.viz?.render();
    if (state.viz) {
      $('zoomLevelDisp').textContent = `${Math.round(state.viz.zoom * 100)}%`;
    }
  });
}

// ─── Exports ──────────────────────────────────────────────────────────────────
function bindExports() {
  $('exportMaskBtn')?.addEventListener('click', async () => {
    if (!state.engine) return;
    try {
      setStatus('Generating mask PNG…');
      const url = await state.engine.exportMaskPNG();
      downloadURL(url, 'MCseg_mask.png');
      setTimeout(() => URL.revokeObjectURL(url), 10_000);
      setStatus('Mask exported.');
    } catch (err) {
      setStatus(`Export error: ${err.message}`);
    }
  });

  $('exportCsvBtn')?.addEventListener('click', () => {
    if (!state.engine) return;
    const csv  = state.engine.exportCellCSV();
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
    const url  = URL.createObjectURL(blob);
    downloadURL(url, 'MCseg_cells.csv');
    setTimeout(() => URL.revokeObjectURL(url), 10_000);
    setStatus(`Exported ${state.engine.cells.length.toLocaleString()} cells.`);
  });
}

// ─── Theme ────────────────────────────────────────────────────────────────────
function bindTheme() {
  const btn = $('themeToggle');
  // Respect system preference on first load
  if (window.matchMedia?.('(prefers-color-scheme: dark)').matches) {
    document.body.classList.add('dark');
    btn.textContent = '☀️';
  }
  btn?.addEventListener('click', () => {
    const isDark = document.body.classList.toggle('dark');
    btn.textContent = isDark ? '☀️' : '🌙';
  });
}

// ─── Progress helpers ─────────────────────────────────────────────────────────
function showProgress(pct, label) {
  const container = $('progressContainer');
  container.classList.remove('hidden');
  $('progressFill').style.width = `${pct}%`;
  $('progressLabel').textContent = label ?? '';
  setStatus(label ?? '');
}

function updateProgressAttr(pct) {
  $('progressFill').setAttribute('aria-valuenow', pct);
}

function hideProgress() {
  setTimeout(() => $('progressContainer').classList.add('hidden'), 1200);
}

function setStatus(msg) {
  $('statusMsg').textContent = msg ?? '';
}
