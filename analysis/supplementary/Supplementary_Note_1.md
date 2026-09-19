# Supplementary Note 1: MCseg workflow and implementation reference

This note describes the current repository implementation. Historical benchmark runs must be reproduced with their original revision, model weights, and effective settings; the defaults below are not a substitute for those run records.

## Development and deployment

Researchers defined candidate operations, registered Xenium reference data, AP@0.5 scoring, prompts, and execution limits. An AI agent proposed and evaluated candidate workflows in a constrained sandbox, followed by researcher review. The retained workflow runs locally without a language-model API or a new architecture search. The files in `docs/autoResearch/` are adaptation templates, not proof of a complete historical execution archive.

## Preprocessing and model

The engine applies CLAHE (`clip_limit=3.0`, tile size 8 by default), a tissue mask, and HED-based hematoxylin extraction. All passes use **Cellpose cpsam**; variations are in image representation, diameter, probability threshold, and inference settings. The primary loader records the resolved model-weight path.

## Pass configuration

The following table reflects the CRC profile and current implementation. `c` denotes the effective `cellprob_threshold`, normally -2.0 for CRC. The actual execution order is mid, small, large, hematoxylin, then the additional passes.

| Pass | Model | Input | Diameter | cellprob threshold |
| --- | --- | --- | --- | --- |
| Primary mid | cpsam | CLAHE RGB | `dia_mid` (CRC: 17 px) | c (-2.0) |
| Primary small | cpsam | CLAHE RGB | `dia_small` (CRC: 13 px) | c - 1 (-3.0) |
| Primary large | cpsam | CLAHE RGB | `dia_large` (CRC: 22 px) | c + 1 (-1.0) |
| Primary hematoxylin | cpsam | Hematoxylin replicated across channels | `dia_mid` (CRC: 17 px) | c (-2.0) |
| Additional RGB auto | cpsam | CLAHE RGB | `dia_cpsam_auto` (0 maps to `None`) | `cellprob_cpsam_auto` (-1.0) |
| Additional RGB small | cpsam | CLAHE RGB | `dia_cpsam_small` (16 px) | `cellprob_cpsam_small` (-3.0) |
| Additional hematoxylin | cpsam | Hematoxylin replicated across channels | `dia_cpsam_auto` (0 maps to `None`) | `cellprob_cpsam_hema` (-1.0) |

`flow_threshold` defaults to 0.4. Hematoxylin passes require `use_hematoxylin=true`; additional passes require `use_cpsam=true`. With hematoxylin enabled, the shorter configuration has four passes and the expanded configuration has seven. The additional flag does not switch model families. Primary mid/hematoxylin passes use augmentation, while the remaining passes do not; primary passes enable resampling and additional passes disable it.

## Integration, rescue, and expansion

The engine merges masks in a defined priority order using overlap filtering; it does not select a higher-confidence detection from a generic confidence score. Optional transcript-density rescue depends on available transcript-density input. Voronoi-constrained expansion and area filtering follow. CRC and LUAD profiles specify different expansion distances (9 and 8 px, respectively); pixel size and effective overrides must accompany reported distances. Area filtering defaults to 20–6000 px².

Review the source for the exact operation ordering and overlap rule:

- [Segmentation engine](../../backend/src/segmentation/cellpose_runner.py)
- [CRC profile](../../config/profiles/crc.yaml) and [LUAD profile](../../config/profiles/luad.yaml)
- [Pipeline configuration](../../config/pipeline.yaml)

The configuration merge order is tissue profile, pipeline settings, then runtime state. Per-ROI overrides may further change effective parameters. Do not infer the settings of a historical benchmark solely from a present-day default.

## Transcript attribution

Bin centroids are transformed into image/ROI coordinates and assigned to corresponding mask labels. Counts are aggregated into a sparse cell-by-gene matrix. The application's counting stage supports additional expansion through `rna_counting.dilation_px`; record it because it changes the attribution geometry. Unassigned bins and the denominator used for FTC must be handled consistently across comparisons.

## Evaluation scope

The September 18 manuscript reports fixed-parameter LUAD development-set PQ of **0.472 ± 0.072**, versus **0.432 ± 0.037** for tuned 2Cseg. Separately selecting expansion strategy and distance per ROI with Xenium masks gives **0.554 ± 0.063**. This latter value is a reference-guided upper bound, not ordinary deployment performance or an independent test-set estimate.

For the CRC results, UMI density and fractional transcript capture are distinct. NED is a reference-free expression-separation measure and should be interpreted with geometry, capture, mask area, and lineage mixing. See the [README results table](../../README.md#research-results).

## Reproduction checklist

Retain the source revision, dependency lock, resolved model-weight identity, dataset identifiers and ROI coordinates, effective segmentation/attribution parameters, and logs. The desktop `BUILD_MANIFEST.json` records packaged file hashes and source revision; it does not archive model weights or independently validate the manuscript's numerical results.
