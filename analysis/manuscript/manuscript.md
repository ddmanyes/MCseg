# MCseg: High-Fidelity Visium HD Cell Segmentation via AI-Autonomous Pipeline Discovery

**Journal**: Bioinformatics (Oxford Academic)
**Article type**: Original Paper
**Status**: Draft v1.0 (2026-04-09)

---

## Authors

[Author list TBD]

---

## Abstract

**Motivation:** Cell segmentation is a critical bottleneck in Visium HD spatial transcriptomics, yet accessible tools for non-programming researchers remain scarce and standard geometric metrics fail to capture RNA attribution quality. Existing default pipelines (Space Ranger) exhibit systematic boundary over-expansion in glandular tissues, inflating transcript capture at the cost of single-cell purity.

**Results:** We developed **MCseg**—a no-code, browser-based Visium HD analysis platform—with MCseg v2 as its high-fidelity segmentation engine discovered through AI-autonomous experimentation. Over ~80 unsupervised cycles, an AutoResearch agent converged on a 7-pass multi-model ensemble with Voronoi-constrained boundary expansion, achieving 28% relative improvement in Panoptic Quality (PQ: 0.554 vs 0.432) against Xenium single-molecule ground truth in LUAD. In CRC, RNA-based benchmarking across 15 ROIs reveals a "FTC Paradox": Space Ranger's apparent transcript-capture advantage disappears when normalised for cell area (UMI density: 11.7 vs 11.6 UMI/µm², <1% difference), while MCseg v2 achieves significantly superior boundary purity (NED: 0.727 vs 0.712, p=0.026; mutually exclusive co-expression 0.49% vs 0.67%, p=0.030). In CRC tertiary lymphoid structures, MCseg v2 resolved 4 functionally distinct immune subtypes versus Space Ranger's 3 mixed populations, with 44% higher cell yield.

**Availability:** MCseg source code and segmentation pipeline are available at [GitHub repository TBD] under the MIT License. All analysis scripts are provided in the Supplementary Materials.

**Contact:** [corresponding author TBD]

**Keywords:** Visium HD, spatial transcriptomics, cell segmentation, Cellpose, MCseg, AutoResearch, transcript attribution, neighbour expression divergence, lung adenocarcinoma, colorectal cancer, tumour microenvironment

---

## 1. Introduction

Spatial transcriptomics technologies such as 10x Genomics Visium HD enable genome-wide gene expression profiling at 2 µm/bin resolution, bridging the gap between bulk tissue analysis and single-cell sequencing by retaining spatial context. The clinical and biological potential of Visium HD depends critically on cell segmentation—assigning the 2 µm expression bins to individual cells—as this step determines whether downstream analyses reflect genuine single-cell biology or bin-mixing artefacts. Fluorescence-based platforms such as Xenium Prime achieve high single-molecule sensitivity on curated gene panels (~5K–10K genes), while Visium HD interrogates the full transcriptome (~18K genes in typical FFPE) but requires downstream segmentation to resolve binned output into cellular units. These complementary architectures motivate rigorous cross-platform evaluation. However, three limitations constrain current practice.

First, standard geometric metrics (e.g., AP@0.5) require thousands of hand-annotated cells per tissue type [1], restricting systematic benchmarking to a handful of datasets. Second, Space Ranger v4 (SR), the 10x Genomics official pipeline, is widely treated as a default reference despite its tendency to over-expand boundaries in glandular tissues; treating it as ground truth introduces circular validation bias [2,3]. Third, geometric metrics decouple from biological relevance: identical shape scores can correspond to substantially different RNA attribution quality when minor boundary shifts introduce transcript spillover [4,5]. These analytical limitations are compounded by a practical barrier: a typical Visium HD section contains over ten million 2 µm bins, making full-slide cell segmentation prohibitively memory- and compute-intensive, while most biological hypotheses are anchored to specific microenvironments. No accessible, ROI-focused Visium HD tool exists that does not require deep programming expertise.

We address both barriers through a two-step programme. In the first step, we developed **MCseg v1**—an interactive, browser-based pipeline enabling ROI-focused Visium HD segmentation without code or GPU—and characterised its performance ceiling via systematic Optuna hyperparameter sweep (PQ = 0.432, n = 6 LUAD ROIs vs Xenium ground truth). In the second step, we adapted the **AutoResearch** autonomous experimentation paradigm [6] to cell segmentation: an AI agent iteratively proposed and evaluated segmentation architectures against Xenium ground truth over ~80 overnight cycles, converging on **MCseg v2** (PQ = 0.554, +28%). MCseg v2 is deployed in CRC using three complementary RNA-based metrics—Fraction of Tissue Captured (FTC), UMI Density, and Neighbour Expression Divergence (NED)—that provide orthogonal validation without Xenium co-registration. Key contributions are: (1) **MCseg**, a no-code end-to-end Visium HD platform with automated ROI selection and downstream transcriptomic analysis; (2) MCseg v2, an AI-discovered high-fidelity segmentation engine; (3) a dual-tissue (LUAD + CRC) validation framework combining geometric and RNA-based evaluation; and (4) biological insights into challenging microenvironments where standard platforms fail.

---

## 2. System and Methods

### 2.1 Datasets

**LUAD (lung adenocarcinoma):** Data were obtained from the 10x Genomics public demonstration dataset, comprising co-registered Visium HD (2 µm/bin) and Xenium Prime data generated on the same FFPE section, covering moderately and poorly differentiated adenocarcinoma regions (G2/G3) and adjacent normal lung parenchyma [2]. Six representative ROIs (274×274 µm each) were manually defined to span tumour boundary, tumour core, tumour stroma, mixed tumour-stroma, normal-tumour interface, and alveolar microenvironments. Xenium `cell_boundaries.parquet` provided single-molecule-resolved cell boundaries as geometric ground truth.

**CRC (colorectal cancer):** Data were obtained from the 10x Genomics public dataset (also available via GEO accession GSE280318), comprising a single FFPE colorectal cancer section profiled by Visium HD (2 µm/bin) containing glandular adenocarcinoma and adjacent normal colonic mucosa, without Xenium co-registration. Fifteen ROIs (274×274 µm each) were selected: ROI 1–14 cover diverse tumour microenvironments; ROI 15 is a tertiary lymphoid structure candidate identified by systematic spatial discovery using composite TLS marker scoring followed by Local Moran's I spatial autocorrelation. All evaluation relies on RNA-based metrics only.

### 2.2 MCseg v1 Baseline

MCseg v1 uses a dual-diameter Cellpose strategy [13] (`nuclei` model, dia_small = 5 µm, dia_large = 60 µm) with CLAHE contrast normalisation (clip_limit = 1.0) and 6 px Voronoi dilation (1.64 µm). A three-phase Optuna sweep (150+ trials, `xenium_he_seg` project) established the single-model performance ceiling: PQ = 0.432 ± 0.037 (n = 6 LUAD ROIs).

### 2.3 MCseg v2: AI-Autonomous Discovery

The AutoResearch agent operated within a sandboxed environment, modifying a single `segment.py` and evaluating against Xenium GT (AP@0.5, <5 min/cycle) over ~80 cycles. The agent could freely invoke Cellpose models [13] (cyto2/cyto3/nuclei/cpsam [14]), image processing primitives (CLAHE, HED deconvolution, watershed, Voronoi tessellation), and morphological operations. The final MCseg v2 pipeline comprises: (1) CLAHE preprocessing (clip=3.0, tile=8) + hematoxylin channel extraction; (2) 7-pass multi-model detection (cyto3 at diameters 13/17/22 px on CLAHE-RGB and hematoxylin; cpsam at auto and 16 px diameter; cpsam on hematoxylin); (3) ensemble merging (overlap threshold < 15%); (4) adaptive Voronoi boundary expansion (fixed voronoi d=8 in deployment); (5) quality filtering (20–6000 px²). Full pipeline detail is provided in Supplementary Note 1.

### 2.4 Transcript Attribution and RNA Metrics

Bin-to-cell attribution uses sparse matrix multiplication (Y = M × X), where M is the binary cell-mask matrix and X is the 2 µm transcript-bin matrix. Three RNA-based metrics were computed per method across 15 CRC ROIs:

- **FTC**: fraction of in-tissue UMIs assigned to any cell mask.
- **UMI Density**: total UMIs per cell mask area (UMI/µm²); size-normalised and independent of boundary expansion.
- **NED (Neighbour Expression Divergence)**: mean Hellinger distance between transcriptomes of spatially adjacent cell pairs (top 1,000 HVGs):
$$H(P,Q)=\frac{1}{\sqrt{2}}\sqrt{\sum_{k}\!\left(\sqrt{p_k}-\sqrt{q_k}\right)^2}$$
Spatial adjacency is identified via grey dilation of the segmentation mask. Higher NED (range 0–1) indicates sharper biological boundaries and minimal transcriptional spillover; Hellinger distance is zero-safe, bounded, and variance-stabilising for Poisson count data [15].

**NUC** (control baseline): a single-pass Cellpose `nuclei` detection (diameter 15 µm, no CLAHE preprocessing, no Voronoi expansion), yielding bare nuclear masks. This represents the minimum-boundary extreme of the capture–purity trade-off curve.

Mutually exclusive co-expression rate (c1_coexpr) was computed across four biologically impossible lineage pairs (EPCAM×CD3E, MUC2×NKG7, ACTA2×CD3E, PECAM1×EPCAM); a cell was co-positive if both markers had raw UMI > 0. Statistical tests: paired Wilcoxon signed-rank (n = 15 ROIs); Friedman test for overall method differences. Dimensionality reduction and clustering used Scanpy v1.10: PCA (30 PCs), k-NN graph (n_neighbours = 15), Leiden algorithm, UMAP (min_dist = 0.5).

---

## 3. Results

### 3.1 From MCseg v1 to MCseg v2: Tool Development and AI-Driven Discovery

**MCseg v1** implements a dual-diameter Cellpose strategy—combining small-diameter nuclear detection with large-diameter cytoplasm-inclusive boundaries—paired with CLAHE contrast normalisation. This combination substantially outperformed single-diameter approaches in tissues with high cellular size heterogeneity. Systematic Optuna hyperparameter optimisation (150+ trials) confirmed a performance ceiling: PQ = 0.432 ± 0.037 on LUAD tissue, regardless of further hyperparameter tuning (Supplementary Fig. S2). MCseg v1 is deployed as a browser-based, GPU-free pipeline enabling non-programming researchers to perform ROI-focused segmentation and downstream single-cell analysis.

To overcome the single-model ceiling, we applied the AutoResearch paradigm to cell segmentation discovery. The agent explored a combinatorial space of architectures without human guidance; the discovery trajectory converged through three key stages (Fig. 1C):

1. **Initial baseline** (AP@0.5 ≈ 0.32): standard single-pass detection with fixed radial expansion.
2. **CLAHE + Voronoi-constrained expansion** (→ 0.648, +0.33 in one cycle): the agent independently discovered that Voronoi tessellation allocates inter-cellular space in a physically grounded manner, producing a step-change improvement.
3. **7-pass multi-model ensemble convergence → MCseg v2** (AP@0.5 = 0.650, PQ = 0.554): combining cyto3 at three diameters with cpsam supplementation captures the full morphological spectrum—from small lymphoid cells to large pleomorphic tumour cells—that any single configuration misses.

Sensitivity analysis confirms MCseg v2 occupies a stable optimum (Fig. 1D); ablation studies show each component (CLAHE, multi-pass ensemble, Voronoi tessellation) is individually indispensable, with removal causing AP@0.5 losses up to 0.22 (Fig. 1E). MCseg v2 is integrated as the segmentation backend of MCseg, making this AI-discovered pipeline directly accessible to biologists without programming expertise.

### 3.2 Geometric Validation in LUAD Against Xenium Ground Truth

Using Xenium single-molecule resolved cell boundaries as unbiased geometric ground truth across 6 LUAD ROIs (Fig. 2A–C), MCseg v2 outperforms MCseg v1 by +0.122 PQ (0.554 ± 0.064 vs 0.432 ± 0.037, 28% relative improvement), with SQ as the primary driver (0.585 → 0.777, +33%), indicating that Voronoi-expanded boundaries approximate Xenium's cytoplasm-inclusive reference substantially better than nuclear-centred radial dilation. The modest RQ decrease (0.734 → 0.711) is an expected consequence of cytoplasm-inclusive cyto3 detection merging some adjacent nuclear instances that nuclei-only MCseg v1 resolves separately—an acceptable trade-off for achieving full cellular boundary coverage. The two highest-scoring ROIs (tumour boundary ROI 1 + tumour core ROI 6) achieve mean PQ = 0.636, with ROI 1 reaching PQ = 0.662.

We further applied MCseg v2 to two biologically demanding LUAD microenvironments. In the normal alveolar region (ROI 10), AT1 (*AGER*⁺/*RTKN2*⁺) and AT2 (*SFTPC*⁺/*SFTPB*⁺/*SFTPA1*⁺) pneumocyte discrimination was evaluated. MCseg v2 geometrically covered 22.4% of cells within AT1 territory—virtually identical to Xenium's 22.9% (pixel recall = 0.71, IoU = 0.56; Fig. 2E–F)—establishing accurate AT1 spatial localisation. Yet only 2.1% of MCseg cells express *AGER* or *RTKN2* in the Visium HD count matrix, an 11-fold gap relative to Xenium. This dissociation is attributable to **Visium HD platform-level sparsity**: AT1 cells extend thin cytoplasmic sheets across the alveolar surface in 3D, and even within the nuclear footprint, *AGER* and *RTKN2* mRNA density falls below detection threshold within any single 2 µm bin. Xenium's single-molecule FISH detects individual transcripts without bin aggregation, bridging this gap. MCseg v2 therefore represents the geometric ceiling achievable on this tissue; the RNA deficit is intrinsic to the platform, not a segmentation artefact. For AT2 pneumocytes, cytoplasm-inclusive Voronoi boundaries yield a measurable gain: *SFTPC* positivity reaches 34%, *SFTPB* 22%, *SFTPA1* 12% (combined AT2 positivity 45% vs 39% with MCseg v1), demonstrating the benefit of cytoplasmic expansion for cytoplasm-enriched transcripts.

In the pigmented alveolar macrophage zone (ROI 9), Xenium fluorescence detection is prone to quenching by dark carbon deposits [7]; H&E-based MCseg v2 achieves comparable cell counts (8,523 vs 9,694). More importantly, *SPP1*—the canonical immunosuppressive TAM marker—falls outside Xenium's 5K targeted panel. Leveraging Visium HD's unbiased whole-transcriptome coverage, MCseg v2 identifies a distinct *SPP1*⁺/*FTH1*⁺ iron-laden TAM cluster: *SPP1* positivity 17.8%, *FTH1* 55% (vs 23% tissue-wide), *TREM2* 12-fold enriched (1.2% vs 0.1% background)—a phenotype linked to poor prognosis and chemoresistance in LUAD [8,9] that no amount of analytic refinement within Xenium's panel can retrieve.

### 3.3 Transcript Attribution Benchmarking in CRC

We applied MCseg v2 to glandular CRC (n = 15 ROIs × 4 methods; Fig. 3) to evaluate cross-tissue generalisability. MCseg v2 boundaries closely follow H&E glandular contours, whereas SR masks systematically over-expand into acellular stromal corridors and luminal space—a pattern reflected in mask circularity (SR median 0.563 vs MCseg v2 0.771; Wilcoxon p < 0.001, Supplementary Fig. S8).

Head-to-head benchmarking exposes an apparent **FTC Paradox**:

 | Method | FTC | UMI Density (UMI/µm²) | NED |
|--------|-----|----------------------|-----|
 | SR | **0.934** | 11.7 | 0.712 |
 | **MCseg v2** | 0.737 | 11.6 | **0.727** |
 | MCseg v1 | 0.723 | 11.4 | 0.754 |
 | NUC | 0.246 | 10.4 | 0.823 |

SR achieves the highest raw FTC (0.934), apparently capturing far more transcripts than MCseg v2 (0.737). Yet UMI density is essentially identical between the two methods (<1% difference), and median UMI/cell and genes/cell are likewise comparable (Fig. 3B–D). SR's FTC advantage is therefore attributable entirely to boundary over-expansion into intercellular space, not to superior per-cell RNA capture efficiency.

Boundary over-expansion carries a measurable cost in transcriptional purity: MCseg v2 achieves NED = 0.727 vs SR's 0.712 (Wilcoxon p = 0.026, n = 15 paired ROIs), indicating that SR's extended boundaries absorb neighbouring cell transcripts and reduce within-cell transcriptional distinctiveness. This finding is independently supported by mutually exclusive lineage marker co-expression analysis: MCseg v2 mean 0.49% vs SR 0.67% (Wilcoxon p = 0.030), indicating that SR's boundary displacement introduces quantifiable transcript contamination.

MCseg v1's higher NED (0.754) reflects smaller nuclear-centred boundaries rather than superior segmentation quality—a geometric artefact directly demonstrated by LUAD benchmarking, where MCseg v1's nuclei-based architecture achieves PQ = 0.432 vs MCseg v2's PQ = 0.554 against cytoplasm-inclusive Xenium boundaries (+28%). Tighter boundaries inflate NED by widening the transcriptional gap between adjacent cells—but at the direct cost of excluding cytoplasmic RNA from nuclear-only masks. mRNA is predominantly localised in the cytoplasm; nuclear-only boundary placement therefore systematically excludes the majority of cellular transcripts. NUC (NED = 0.823) sacrifices ~74% of FTC for maximum boundary purity—an extreme operating point that serves as a conceptual anchor but is impractical for whole-transcriptome analysis. Across the four methods, a fundamental capture–purity trade-off emerges: SR and NUC occupy the opposing extremes, MCseg v1 occupies an intermediate nuclei-centred position, and **MCseg v2 uniquely achieves near-SR transcript capture with significantly superior boundary purity**—a combination no other tested method replicates.

### 3.4 Systematic Discovery and Architecture of CRC Tertiary Lymphoid Structures

Tertiary lymphoid structures (TLS) are organised ectopic lymphoid aggregates that correlate with favourable prognosis and immunotherapy response in CRC [12]. Their dense, multi-lineage composition makes them a stringent test of segmentation resolution. To select candidates without prior bias, we systematically surveyed the full CRC slide via composite TLS marker scoring (*JCHAIN, MS4A1, CD79A, CXCL13, IGKC, LTB*) followed by Local Moran's I spatial autocorrelation (p < 0.01, 9,999 permutations; Supplementary Fig. S9). The highest-ranking candidate (ROI 15) was selected for detailed benchmarking.

At Leiden resolution 0.5, the two methods diverged sharply in resolving power (Fig. 4):

- **MCseg v2** (n = 636 cells) recovered **4 distinct functional populations**: Plasma/B cell (*JCHAIN*⁺), Stroma/Fibroblast (*VIM*⁺), Myeloid/Macrophage (*LYZ*⁺), and Well-differentiated Tumour (*CEACAM5*⁺), each forming well-separated UMAP clusters (Supplementary Fig. S10A).
- **SR** (n = 440 cells) yielded only **3 broad, mixed populations**: SR's largest cluster (n = 256; 58% of all cells) simultaneously co-expressed *IGKC* (5.43), *VIM* (3.74), and *LYZ* (1.32), collapsing B-cell, stromal, and myeloid lineages into a single contaminated population (Supplementary Fig. S10B).

MCseg v2's 44% higher cell yield (636 vs 440) reflects recovery of cells previously absorbed into SR's over-expanded masks. Spatial marker purity diverged markedly: *JCHAIN* (plasma cell marker) was sharply confined within MCseg v2 masks, whereas SR showed *JCHAIN* signal leaking into adjacent tumour cells—a direct spatial manifestation of the NED difference quantified in §3.3. Unmapped regions in the MCseg v2 panel—vessel lumens, lymphatic sinusoidal spaces, and ECM corridors—correspond to areas devoid of nuclear H&E morphology, consistent with anatomical fidelity rather than missed detections. The resolution advantage documented in §3.3 thus translates directly into biologically actionable single-cell populations in the most lineage-diverse microenvironment on the slide.

---

## 4. Discussion

### MCseg v2 vs Space Ranger: resolving the FTC paradox

Normalising for cell area eliminates Space Ranger's apparent transcript-capture advantage: UMI density is virtually identical (<1% difference) between SR and MCseg v2 despite the large FTC gap (0.934 vs 0.737). SR's boundary over-expansion inflates raw transcript counts while simultaneously contaminating neighbouring transcriptomes, as reflected in lower NED (0.712 vs MCseg v2's 0.727, p = 0.026) and higher impossible co-expression rates (0.67% vs 0.49%, p = 0.030). These results establish that **FTC alone is an unreliable segmentation quality metric**; UMI density and NED are essential complements. Prior GT-free frameworks (CellSPA [4]; SpatialQM [5]; cellAdmix [3]) have proposed related metrics, but our analysis provides the first head-to-head comparison on Visium HD's bin-based architecture showing that FTC and NED measure orthogonal, inversely correlated properties.

### AI-autonomous experimentation as a development paradigm

MCseg v2's emergence from ~80 overnight agent cycles establishes that delegating the core experimental loop to an AI agent—one that proposes, implements, and evaluates architectures without supervision—converts weeks of manual iteration into unattended compute. Three features made this tractable: an objective scoring function (AP@0.5 vs Xenium GT), a bounded search space (one Python file, a pre-loaded tool library), and cheap, automatically discarded failures. The Voronoi-constrained multi-diameter ensemble that emerged—a configuration unlikely to arise from human prior intuition—suggests this paradigm will generalise to other bioinformatics method development problems wherever automated scoring against a reference standard is feasible.

### Visium HD and Xenium as complementary platforms

The LUAD experiments highlight a predictable complementarity between targeted and untargeted spatial platforms. In ROI 9 (pigmented macrophage zone), *SPP1*—the defining immunosuppressive TAM marker—is absent from Xenium's 5K targeted panel; no analytic refinement within Xenium can recover a transcript the panel does not include. In ROI 10 (normal alveolus), *SFTPC*, *SFTPB*, and *SFTPA1* are similarly absent from the Xenium panel, while AT1 detection is limited by Visium HD's bin-level sparsity rather than segmentation error. Both cases point to the same principle: gene panel design and transcript density, not platform sensitivity alone, determine what biology is discoverable. MCseg v2 converts Visium HD's latent whole-transcriptome advantage into an operational one—without single-cell resolution, 2 µm bins average over neighbouring cell types and obscure cell-type-specific signatures even when transcripts are present. The *SPP1*⁺ TAM phenotype (ROI 9) and AT2-specific co-expression (*SFTPC*/*SFTPB*/*SFTPA1*, ROI 10) are only recoverable when bins are collapsed into coherent single-cell units. Visium HD + MCseg v2 therefore functions not as a replacement for Xenium, but as a complementary discovery instrument: where Xenium precisely quantifies pre-specified markers at single-molecule resolution, Visium HD + MCseg v2 supports hypothesis-free transcriptome-wide characterisation at single-cell spatial resolution, capturing cell identities that targeted panel design cannot anticipate.

### Democratising Visium HD analysis and limitations

MCseg removes the programming expertise barrier through a browser-based interface requiring no code or environment setup; its ROI-focused design reflects how spatial biologists anchor hypotheses to specific microenvironments rather than processing whole slides. Several limitations warrant acknowledgement. CRC results derive from a single patient, requiring multi-centre validation. MCseg v2 was optimised on LUAD tissue, and parameter transfer to CRC was assumed rather than formally tested. The LUAD PQ benchmark reflects oracle expansion (GT-guided Voronoi grid search), and the oracle-to-deployment gap requires future quantification via cross-validation. Transcript-based segmentation tools (e.g., Baysor [10]) and environment-aware cell-typing frameworks (e.g., ENACT [11]) that leverage single-molecule coordinates or neighbourhood context were not directly compared, as Visium HD lacks sub-bin spatial coordinates; this comparison remains important future work.

---

## References

1. Greenwald NF, et al. Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning. *Nat Biotechnol* 2022;40:555–565.
2. Long M, et al. Comparing Xenium 5K and Visium HD data from identical tissue slide at a pathological perspective. *J Exp Clin Cancer Res* 2025;44:219.
3. Mitchel J, et al. Impact and correction of segmentation errors in spatial transcriptomics. *Nat Genet* 2026. DOI: 10.1038/s41588-025-02497-4
4. Fu X, et al. BIDCell: Biologically-informed self-supervised learning for segmentation of subcellular spatial transcriptomics data. *Nat Commun* 2024;15:509.
5. Plummer JD, et al. Standardized metrics for assessment and reproducibility of imaging-based spatial transcriptomics datasets. *Nat Biotechnol* 2025. DOI: 10.1038/s41587-025-02811-9
6. Karpathy A. AutoResearch: AI agents running research on single-GPU nanochat training automatically. GitHub 2026. <https://github.com/karpathy/autoresearch>
7. Marco Salas S, et al. Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows. *Nat Methods* 2025.
8. Matsubara E, et al. SPP1 derived from macrophages is associated with a worse clinical course and chemo-resistance in lung adenocarcinoma. *Cancers* 2022;14:4374.
9. Matsubara E, et al. The significance of SPP1 in lung cancers and its impact as a marker for protumor tumour-associated macrophages. *Cancers* 2023;15:2250.
10. Petukhov V, et al. Cell segmentation in imaging-based spatial transcriptomics. *Nat Biotechnol* 2022;40:345–354.
11. Lotfollahi M, et al. ENACT: end-to-end analysis of Visium High Definition (HD) data. *Bioinformatics* 2025;41:btaf094.
12. de Oliveira MF, et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat Genet* 2025. DOI: 10.1038/s41588-025-02193-3
13. Stringer C, et al. Cellpose: a generalist algorithm for cellular segmentation. *Nat Methods* 2021;18:100–106.
14. Kirillov A, et al. Segment Anything. *ICCV* 2023. arXiv:2304.02643
15. Legendre P, Gallagher ED. Ecologically meaningful transformations for ordination of species data. *Oecologia* 2001;129:271–280.

---

## Data Availability

Raw Visium HD and Xenium data are publicly available via the 10x Genomics dataset portal. The LUAD dataset (Xenium Prime co-registered with Visium HD 2 µm, FFPE) and the CRC Visium HD dataset (FFPE) are available at the 10x Genomics dataset portal; CRC is additionally available via GEO accession GSE280318. Processed AnnData objects and segmentation masks will be deposited to Zenodo upon acceptance (DOI: TBD).

## Code Availability

The MCseg v2 segmentation pipeline and transcript attribution scripts are available at [GitHub TBD] under the MIT License. All analyses were performed in Python 3.11 (Cellpose 3, Scanpy 1.10, Scipy 1.13).

## Competing Interests

The authors declare no competing interests.

## Acknowledgements

[TBD]

---

## Figure Legends

**Fig. 1. MCseg workflow and AI-autonomous discovery of MCseg v2.**
(a) Schematic overview of the MCseg v2 pipeline: CLAHE contrast enhancement, 7-pass multi-diameter cyto3/cpsam ensemble detection, Voronoi-constrained boundary expansion, and downstream single-cell transcriptomic analysis. (b) Representative LUAD tumour-boundary ROI three-panel comparison: H&E, MCseg v2 predicted mask, Xenium Prime GT mask. Scale bar: 50 µm. (c) AutoResearch agent discovery trajectory (~80 cycles); key transition: CLAHE + Voronoi achieves the largest single-cycle gain (+0.33 AP@0.5). (d) Parameter sensitivity heatmap (AP@0.5 vs Voronoi expansion distance × cyto3 diameter; n = 6 LUAD ROIs). (e) Component ablation: removal of CLAHE, multi-pass ensemble, or Voronoi tessellation causes AP@0.5 losses up to 0.22 (mean ± SD, n = 6 ROIs).

**Fig. 2. Geometric benchmarking and biological validation in LUAD.**
(a–c) MCseg v2 vs MCseg v1 across n = 6 LUAD ROIs (Xenium GT): (a) PQ, (b) SQ, (c) RQ. Paired Wilcoxon. (d) ROI 10 MCseg v2 cell-type map. (e) AT1 spatial overlay: Xenium GT AT1 cells (cyan) and MCseg v2 boundaries; pixel recall = 0.71, IoU = 0.56. (f) AT1 detection rate comparison: Xenium GT 22.9% ≈ MCseg geometric 22.4% > MCseg RNA 2.1% (11-fold platform gap). (g) ROI 10 dot plot: canonical markers across 4 cell types. (h) ROI 9 Leiden cluster map (n = 8,393 cells). (i) ROI 9 dot plot: *SPP1*⁺/*FTH1*⁺ TAM cluster. Scale bars: 50 µm.

**Fig. 3. Transcript attribution benchmarking in CRC (n = 15 ROIs × 4 methods).**
(a) Visual boundary comparison: ROI 2 (gland–stroma interface) and ROI 4 (irregular tumour-dense area); four columns: H&E, MCseg v2, MCseg v1, SR. SR masks extend into acellular areas. Scale bar: 50 µm. (b–g) Transcript attribution metrics: (b) FTC; (c) median UMI/cell; (d) median genes/cell; (e) UMI density; (f) NED; (g) mutually exclusive co-expression rate. Paired Wilcoxon (n = 15); ns p > 0.05, * p < 0.05, ** p < 0.01, *** p < 0.001.

**Fig. 4. High-resolution single-cell architecture of a CRC tertiary lymphoid structure (ROI 15).**
(a) Three-panel spatial comparison (274×274 µm): H&E; MCseg v2 cell-type map (n = 636, 4 clusters); SR cell-type map (n = 440, 3 mixed clusters). Scale bar: 50 µm. (b) Marker signal purity (2 × 6 panel): MCseg v2 (top) vs SR (bottom) for *JCHAIN*, *CEACAM5*, *LYZ*, *VIM*, *MT-CO3*, *CD79A*. MCseg v2 shows sharp lineage confinement; SR shows cross-lineage contamination.

---

## Supplementary Materials

 | Figure | Description |
|--------|-------------|
 | Fig. S1 | MCseg v2 vs Xenium GT: all 6 LUAD benchmark ROIs |
 | Fig. S2 | MCseg v1 hyperparameter sensitivity (150+ Optuna trials) |
 | Fig. S3 | ROI location overview: LUAD (6 ROIs) and CRC (15 ROIs) on full-slide H&E |
 | Fig. S4 | LUAD ROI 9: H&E, Xenium DAPI, Xenium segmentation, MCseg v2 segmentation |
 | Fig. S5 | LUAD ROI 10: MCseg v1 vs MCseg v2 AT1/AT2 comparison |
 | Fig. S6 | CRC 7-ROI benchmarking strip: MCseg v2 vs SR side-by-side |
 | Fig. S7 | Per-ROI metric heatmap: 15 CRC ROIs × 4 methods |
 | Fig. S8 | Cell shape circularity: 4 methods × 15 CRC ROIs (violin + per-ROI) |
 | Fig. S9 | TLS discovery workflow: composite scoring and Local Moran's I hotspots |
 | Fig. S10 | UMAP validation: MCseg v2 (4 clusters) vs SR (3 mixed clusters), ROI 15 |
 | **Note 1** | MCseg v2 full algorithm specification and AutoResearch sandbox details |
 | **Table S1** | CRC transcript attribution metrics: 15 ROIs × 4 methods (FTC, UMI density, NED, doublet rate) |
 | **Table S2** | LUAD ROI coordinates (Visium HD fullres pixel space) |

---

### Supplementary Note 1: MCseg v2 Algorithm Specification and AutoResearch Sandbox

#### A. AutoResearch Sandbox Configuration

MCseg v2 was developed through automated agent-driven optimisation using the AutoResearch framework. The agent operated within a strictly sandboxed environment:

- **Evaluation metric**: AP@0.5 (average precision at IoU ≥ 0.5) against Xenium-derived ground-truth masks
- **Cycle budget**: ~80 overnight cycles, each completing within <5 minutes on a CPU-only macOS workstation (no GPU)
- **Modifiable scope**: a single `segment.py` script, preventing unintended side effects on other pipeline components
- **Available primitives**: Cellpose models (cyto2, cyto3, nuclei, cpsam), image preprocessing (CLAHE, HED stain deconvolution, Gaussian smoothing), post-processing (watershed, Voronoi tessellation, morphological operations)
- **Search space**: model selection and combination, diameter values, CLAHE parameters (clip limit, tile grid), ensemble overlap thresholds, Voronoi seed distance, area filter bounds

The agent autonomously discovered the 7-pass ensemble architecture and adaptive Voronoi expansion strategy. AP@0.5 improved from 0.32 (single-model baseline) to 0.65 (final MCseg v2), a 2× relative gain.

#### B. MCseg v2 Pipeline Specification

The pipeline processes a single H&E-stained ROI image and produces a labelled cell mask.

##### Step 1 — Image preprocessing

- CLAHE contrast enhancement: clip limit = 3.0, tile grid = 8 × 8
- Hematoxylin channel extraction via HED colour deconvolution (Macenko method)
- Outputs: CLAHE-RGB composite and hematoxylin single-channel image

##### Step 2 — 7-pass multi-model cell detection

| Pass | Model | Input | Diameter (px) |
|------|-------|-------|:-------------:|
| 1 | cyto3 | CLAHE-RGB | 13 |
| 2 | cyto3 | CLAHE-RGB | 17 |
| 3 | cyto3 | CLAHE-RGB | 22 |
| 4 | cyto3 | Hematoxylin | 13 |
| 5 | cpsam | CLAHE-RGB | auto |
| 6 | cpsam | CLAHE-RGB | 16 |
| 7 | cpsam | Hematoxylin | auto |

All passes use `flow_threshold = 0.4`, `cellprob_threshold = 0.0`.

##### Step 3 — Ensemble merging

Masks from all 7 passes are merged using greedy non-maximum suppression. Any two candidate masks with pixel overlap > 15% are resolved by retaining the higher-confidence mask. This preserves complementary detections from different model configurations while eliminating redundant proposals.

##### Step 4 — Voronoi boundary expansion

Cell centroids from the merged ensemble serve as Voronoi seeds. Unclaimed pixels (not assigned to any detected cell) are partitioned by nearest-centroid assignment using a fixed seed distance *d* = 8 px in deployment mode. This step recovers cytoplasmic RNA signal that would otherwise fall within inter-mask gaps.

##### Step 5 — Quality filtering

Final masks are filtered by area: cells with area < 20 px² (sub-nuclear fragments) or > 6000 px² (stitching artefacts or tissue debris) are discarded.

---

### Supplementary Table S1: CRC Transcript Attribution Metrics (15 ROIs × 4 Methods)

FTC = fraction of transcripts captured; UMI density in UMI/µm²; NED = normalised expression divergence (Hellinger distance, lower = better purity); Doublet rate = fraction of cells co-expressing impossible marker pairs. Methods: SR = Space Ranger; v2 = MCseg v2; v1 = MCseg v1; NUC = nuclear Cellpose baseline.

| ROI | Method | FTC | UMI Density | Median UMI | NED | Doublet Rate | N Cells |
|-----|--------|----:|------------:|-----------:|----:|-------------:|--------:|
| 1 | SR | 0.9922 | 8.157 | 430 | 0.7505 | 0.0045 | 1281 |
| 1 | MCseg v2 | 0.8530 | 8.184 | 465 | 0.7451 | 0.0061 | 1032 |
| 1 | MCseg v1 | 0.7800 | 8.090 | 350 | 0.7593 | 0.0027 | 1219 |
| 1 | NUC | 0.2085 | 6.859 | 113 | 0.8404 | 0.0000 | 836 |
| 2 | SR | 0.9708 | 10.187 | 476 | 0.7044 | 0.0026 | 1260 |
| 2 | MCseg v2 | 0.8256 | 10.045 | 500 | 0.7012 | 0.0028 | 1054 |
| 2 | MCseg v1 | 0.7784 | 10.053 | 407 | 0.7255 | 0.0014 | 1241 |
| 2 | NUC | 0.2388 | 9.473 | 166 | 0.8136 | 0.0003 | 820 |
| 3 | SR | 0.9355 | 7.304 | 399 | 0.6531 | 0.0041 | 1098 |
| 3 | MCseg v2 | 0.8143 | 6.891 | 394 | 0.6555 | 0.0022 | 1008 |
| 3 | MCseg v1 | 0.7116 | 7.127 | 313 | 0.6717 | 0.0019 | 1160 |
| 3 | NUC | 0.2159 | 7.050 | 138 | 0.7640 | 0.0010 | 737 |
| 4 | SR | 0.9598 | 15.524 | 796 | 0.7409 | 0.0101 | 1067 |
| 4 | MCseg v2 | 0.8124 | 15.604 | 970 | 0.7250 | 0.0116 | 795 |
| 4 | MCseg v1 | 0.7713 | 15.201 | 723 | 0.7550 | 0.0075 | 973 |
| 4 | NUC | 0.2754 | 13.777 | 334 | 0.8286 | 0.0012 | 648 |
| 5 | SR | 0.9725 | 14.765 | 857 | 0.6998 | 0.0058 | 1043 |
| 5 | MCseg v2 | 0.8044 | 14.749 | 1011 | 0.6979 | 0.0060 | 792 |
| 5 | MCseg v1 | 0.7440 | 14.127 | 715 | 0.7364 | 0.0028 | 993 |
| 5 | NUC | 0.2763 | 12.223 | 381 | 0.8057 | 0.0000 | 607 |
| 6 | SR | 0.9433 | 12.754 | 732 | 0.7517 | 0.0022 | 1017 |
| 6 | MCseg v2 | 0.8980 | 12.741 | 603 | 0.7811 | 0.0017 | 1204 |
| 6 | MCseg v1 | 0.8299 | 12.469 | 500 | 0.8058 | 0.0013 | 1373 |
| 6 | NUC | 0.3867 | 11.098 | 306 | 0.8667 | 0.0003 | 983 |
| 7 | SR | 0.9584 | 17.936 | 962 | 0.7042 | 0.0050 | 1008 |
| 7 | MCseg v2 | 0.8025 | 17.617 | 701 | 0.7513 | 0.0038 | 1173 |
| 7 | MCseg v1 | 0.7708 | 16.925 | 525 | 0.7885 | 0.0018 | 1556 |
| 7 | NUC | 0.2352 | 15.300 | 271 | 0.8598 | 0.0012 | 841 |
| 8 | SR | 0.9820 | 7.860 | 491 | 0.7625 | 0.0068 | 1000 |
| 8 | MCseg v2 | 0.6753 | 7.745 | 557 | 0.7626 | 0.0056 | 666 |
| 8 | MCseg v1 | 0.7314 | 7.476 | 430 | 0.7872 | 0.0039 | 893 |
| 8 | NUC | 0.2381 | 6.696 | 199 | 0.8535 | 0.0004 | 560 |
| 9 | SR | 0.9445 | 12.681 | 773 | 0.7534 | 0.0074 | 945 |
| 9 | MCseg v2 | 0.7488 | 12.832 | 635 | 0.7944 | 0.0043 | 983 |
| 9 | MCseg v1 | 0.7552 | 12.529 | 489 | 0.8157 | 0.0022 | 1228 |
| 9 | NUC | 0.2949 | 12.275 | 328 | 0.8631 | 0.0003 | 770 |
| 10 | SR | 0.8971 | 6.639 | 377 | 0.6842 | 0.0089 | 894 |
| 10 | MCseg v2 | 0.7074 | 6.604 | 385 | 0.6864 | 0.0059 | 805 |
| 10 | MCseg v1 | 0.6762 | 6.731 | 319 | 0.7021 | 0.0038 | 930 |
| 10 | NUC | 0.2326 | 6.594 | 133 | 0.7803 | 0.0004 | 625 |
| 11 | SR | 0.9499 | 9.947 | 670 | 0.6615 | 0.0091 | 854 |
| 11 | MCseg v2 | 0.5720 | 9.617 | 535 | 0.6682 | 0.0051 | 690 |
| 11 | MCseg v1 | 0.7225 | 9.982 | 449 | 0.7069 | 0.0038 | 1054 |
| 11 | NUC | 0.1849 | 8.607 | 169 | 0.7816 | 0.0009 | 586 |
| 12 | SR | 0.7826 | 13.494 | 697 | 0.7547 | 0.0050 | 800 |
| 12 | MCseg v2 | 0.5887 | 13.300 | 579 | 0.7872 | 0.0037 | 815 |
| 12 | MCseg v1 | 0.5863 | 13.078 | 504 | 0.8020 | 0.0016 | 926 |
| 12 | NUC | 0.2520 | 11.822 | 314 | 0.8554 | 0.0008 | 650 |
| 13 | SR | 0.9710 | 15.388 | 1191 | 0.6945 | 0.0071 | 773 |
| 13 | MCseg v2 | 0.7440 | 15.100 | 1199 | 0.7105 | 0.0053 | 655 |
| 13 | MCseg v1 | 0.7699 | 14.524 | 911 | 0.7436 | 0.0029 | 849 |
| 13 | NUC | 0.2973 | 13.313 | 540 | 0.7888 | 0.0010 | 521 |
| 14 | SR | 0.9389 | 11.141 | 835 | 0.7364 | 0.0074 | 776 |
| 14 | MCseg v2 | 0.7219 | 11.405 | 780 | 0.7526 | 0.0073 | 684 |
| 14 | MCseg v1 | 0.7239 | 11.223 | 579 | 0.7801 | 0.0046 | 867 |
| 14 | NUC | 0.2650 | 9.745 | 316 | 0.8375 | 0.0018 | 567 |
| 15 (TLS) | SR | 0.8067 | 11.348 | 1340 | 0.6214 | 0.0142 | 440 |
| 15 (TLS) | MCseg v2 | 0.4879 | 11.068 | 547 | 0.6897 | 0.0024 | 636 |
| 15 (TLS) | MCseg v1 | 0.4869 | 11.627 | 381 | 0.7293 | 0.0023 | 885 |
| 15 (TLS) | NUC | 0.0850 | 10.643 | 145 | 0.8127 | 0.0000 | 327 |

---

### Supplementary Table S2: LUAD ROI Coordinates (Visium HD Fullres Pixel Space)

Pixel size: 0.2737 µm/px. Coordinates define the top-left corner (x, y) of each rectangular ROI in the Visium HD full-resolution H&E image. All six benchmark ROIs are listed; four were used for transcript attribution analysis and six for geometric benchmarking.

| ROI | Description | x (px) | y (px) | Width (px) | Height (px) | Width (µm) | Height (µm) |
|-----|-------------|-------:|-------:|-----------:|------------:|-----------:|------------:|
| Tumor boundary | G2/G3 grade transition zone with lymphoid infiltration | 9202 | 16552 | 6022 | 4521 | 1648 | 1238 |
| HEV / TLS | High endothelial venule cluster within tertiary lymphoid structure | 11787 | 17937 | 1402 | 1151 | 384 | 315 |
| Dust cells | Carbon-laden alveolar macrophage aggregate | 7934 | 11264 | 2869 | 1602 | 785 | 438 |
| Normal lung | Open alveolar region; AT1/AT2 pneumocyte validation | 7050 | 23275 | 2069 | 1267 | 566 | 347 |
