---
name: sci-plot
description: >-
  Generate publication-grade scientific figures (Nature, Cell, Science, Nature Medicine,
  Cancer Cell standard) in Python (matplotlib/seaborn/scanpy) or R (ggplot2/ggrastr/ComplexHeatmap).
  Use this skill whenever the user asks to plot biological, biomedical, clinical, single-cell,
  spatial, omics, or statistical data; choose an appropriate chart type; design a color palette
  for a figure or entire manuscript; export Adobe Illustrator-ready 1:1 millimeter vector PDF/SVG
  with editable Arial typography; or produce hybrid vector+raster figures for high-density data
  (UMAP, volcano plots, large heatmaps).
---

# Scientific Publication Plotting (`sci-plot`)

Distilled empirically from **5,000+ high-IF original research articles** (deep vector/color/typography audit of 167 Tier-1 CNS papers + large-scale chart-selection audit of **4,968 post-2020 tumor single-cell/spatial/multi-omics papers** covering 38,306 captions and 141,092 subpanels + vector label-tilt measurement of 504 axes across 294 Tier-1 PDFs), 172 chart recipes from `scidraw.msttnote.cn`, and official *Nature*/*Cell Press*/*Science* artwork specifications.

## Core Non-Negotiable Rules

1. **1:1 Physical Print Scale (mm, ready for Adobe Illustrator)**:
   - Never plot on default `7×7 in` / `6×4 in` canvases expecting to shrink in AI.
   - Always size each panel directly in **millimeters** to fit target journal grids (`Nature`: 89 / 120 / **180 mm**; `Cell`: 85 / 114 / **174 mm**; `Science`: 57 / 121 / **184 mm**; gutter 2.5–3 mm).
   - Consult [sizing-and-aspect-ratios.md](./references/sizing-and-aspect-ratios.md) for exact slot width, height, and aspect ratio (`W:H`) per chart type.

2. **Arial Typography Ladder, Statistical Symbols & Single-Line AI Text (Locked Default: `cns_compact`)**:
   - **Always default to `cns_compact`** unless the user explicitly asks for `cell_standard` or `arial_8pt`:
     - Font family: **Arial** (fallback Helvetica), including italics, Greek letters, and math.
     - Panel letter (`a, b, c` Nature / `A, B, C` Cell & Science): **8 pt bold** (Cell/subjournals: 9 pt bold)
     - Axis titles / legend titles / short panel titles: **7 pt**
     - Tick labels / category labels / legend items / stat annotations (`p`, `r`, `n`): **6 pt**
     - Dense matrix labels (heatmap/dotplot genes, KM risk tables): **5–5.5 pt** (floor **5 pt**)
     - Axis & frame stroke width: **0.4 pt**; tick length: **2.0 pt**
   - **Statistical Symbols, Lowercase Italic Greek & Exact Spacing** (mandatory — use `fmt_p`, `fmt_q`, `fmt_corr`, `fmt_n`, `fmt_greek`, `fmt_log2fc`, `fmt_neglog10p` in `sci_plot_py.py` / `sci_plot_r.R`):
     - Statistical variables & Greek letters are **lowercase italic**: *p*, *q*, *r*, *n*, *t*, *z*, *ρ*, *α*, *β*, *γ*, *δ*, *μ*, *σ*, *χ*² (only *R*² uses uppercase italic *R*; human genes italic uppercase e.g. *SPP1*; mouse genes italic title-case e.g. *Spp1*).
     - **Single space around operators**: `*p* = 0.004`, `*p* < 0.0001`, `*p* = 1.4 × 10⁻⁹`, `*ρ* = 0.68, *p* = 0.003 (*n* = 18)`, `2.3 ± 0.4`.
     - **Zero space** before `%` (`75%`) or inside parentheses (`(*n* = 18)`); **one space** before physical units (`50 *μ*m`, `24 h`).
   - **Prevent Adobe Illustrator Text Fragmentation (Single-Line Live Text Guarantee)**:
     - **Python**: ALWAYS export via `save_ai_figure(..., save_svg=True, coalesce_text=True)`:
       - Automatically runs `coalesce_pdf_text_runs()` on `.pdf` (fixes per-character `Td` splits + optical `TJ` kerning splits) AND `coalesce_svg_text_runs()` on `.svg` (rewrites matplotlib's scrambled mathtext into ONE clean semantic `<text>` with inline `<tspan>` and zero inner `x`/`y` splits!).
     - **R**: Export via `svglite(..., fix_text_size = FALSE)` or `cairo_pdf` + `--fix-pdf`; NEVER enable `showtext` or `ggtext` for plain labels.
     - **Illustrator-side 1-Click Guarantee for PDF subscripts & multi-style runs**: Because Adobe Illustrator's PDF parser natively splits PDF `Ts` subscripts (`log₂`, `−log₁₀(p_adj)`) and 3+ inline font switches (`Validation (n = 18)`, `ρ = 0.68, p = 0.002`) upon PDF import, ship and recommend [scripts/MergeSciPlotText.jsx](./scripts/MergeSciPlotText.jsx) (`File -> Scripts -> MergeSciPlotText`), which merges all horizontal & 90°-rotated fragments inside Illustrator into native single-line `PointText` objects with true `Arial Italic` and native `baselineShift` subscripts/superscripts in < 0.3 s!

3. **Hybrid Vector + Raster Layering (AI-Safe High-Density Export)**:
   - **Auto-rasterize dense data layers at 600 dpi** while preserving **axes, spines, ticks, legends, colorbar ticks, leader lines, annotations, and Arial text as 100% editable vector objects**:
     - Scatter / UMAP / t-SNE / spatial / volcano / single-cell violin jitter with **> 2,000 points** → rasterize point layer (`rasterized=True` in Python; `ggrastr::rasterise(..., dpi = 600)` in R).
     - Heatmaps / correlation matrices with **> 2,500 cells** → rasterize matrix body (`rasterized=True` in Python `imshow`/`pcolormesh`; `use_raster = TRUE, raster_quality = 5` in `ComplexHeatmap`).
     - Sparse plots (≤ 2,000 points, bar/box/KM/forest/dotplot) → pure vector.
   - **Font embedding guardrails** (detailed in [hybrid-vector-raster-ai.md](./references/hybrid-vector-raster-ai.md)):
     - Python: `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = 'none'`, custom Arial `mathtext`.
     - R: `grDevices::cairo_pdf(..., family = "Arial")` or `svglite::svglite(...)` (NEVER enable `showtext` during editable PDF/SVG export, as it outlines fonts into uneditable paths!).

4. **Compactness First (Fit 3–5 Crisp Panels per 180 mm Row)**:
   - **Eliminate default right-side legends** on multi-panel rows: use (a) direct on-plot labels with short leader lines, (b) legend tucked inside empty plot quadrant with no box, (c) horizontal multi-column legend below/above, or (d) shared color mapping across panels (omit redundant legends).
   - **UMAP / t-SNE (Default: Zemin Zhang Lab "Cloudy / Watercolor" Hybrid Style)**:
     - **Corner L-axes**: ALWAYS draw horizontal (`UMAP1`) and vertical (`UMAP2`) mini-arrows with **EXACTLY EQUAL PHYSICAL LENGTH in mm** (e.g. `8–9 mm` via `add_umap_corner_axes(..., arrow_mm=9.0)` or R `coord_fixed(ratio = 1)`).
     - **Density-calibrated dot size (avoid sparse grainy look)**: Scale point size inversely with cell count (`≥150k`: R `0.06` / Py `0.14`; `50k–150k`: R `0.16` / Py `0.45`; `15k–50k`: R `0.32` / Py `0.95`; `<15k`: R `0.55` / Py `1.8`). For tight sculptural Zemin Zhang islands, compute UMAP with `min.dist = 0.2–0.3, spread = 0.8–1.0`.
     - **Few clusters (`K ≤ 12`, e.g., major lineages or one compartment)**: Use `plot_cloudy_umap(..., label_style="colored")` (`sci_plot_py.py` / `sci_plot_r.R`) — shuffled order + 600-dpi rasterized dots + equal-length vector L-axes + **100% vector same-color direct cluster labels** (subtle white halo, zero legend needed!).
     - **Many subclusters (`K ≥ 15`, e.g., 20–91 fine subclusters across major lineages)**: Use `plot_atlas_umap_composite()` (`sci_plot_py.py` / `sci_plot_r.R`) — combines (1) clutter-free 600-dpi cohesive subcluster UMAP with equal-length vector L-axes, (2) top-left circular clipped 600-dpi major-lineage inset with crisp `0.45 pt` vector ring, and (3) tight bottom 100% vector multi-column subtype legend (`~2.55 mm/row` pitch!) with colored major-lineage header banners — ALL inside **ONE single 1:1 mm AI-editable PDF** (zero manual PNG/text stitching!).
   - **Minimal text inside figure canvas**: concise axis labels, ticks, short group names, sample sizes (`n`), gene names, and numeric statistics (`P`, `r`, `HR`, `AUC`) only — never explanatory sentences, biological conclusions, or multi-line narrative boxes inside plots.

5. **Default CNS Color System & Two-Tier Color Architecture**:
   - **Locked Default Color System** (empirically #1 across 45/45 Tier-1 CNS papers — Blue & Red 100%, Amber 93%, Sage/Teal 91%, Plum 89%; auto-activated in `sci_plot_py.py` and `sci_plot_r.R`):
     - **Default Categorical Cycle (`CNS_DEFAULT`)**: `#3B6EA5` (Steel Blue) → `#C84B4B` (Brick Coral-Red) → `#2A9D8F` (Teal) → `#E08E45` (Warm Amber) → `#7A5195` (Muted Plum) → `#6BAED6` (Sky Blue) → `#D6729B` (Dusty Rose) → `#7DB977` (Sage) → `#E2C062` (Gold) → `#8C6239` (Taupe).
     - **Default Diverging Heatmap (`rdbu_soft`)**: Soft Blue–OffWhite–Crimson (`#2166AC → #F7F7F7 → #B2182B`).
     - **Default Sequential Expression (`gray_crimson`)**: LightGray–Salmon–Crimson (`#EFEFEF → #FB6A4A → #67000D`) — use `gray_indigo` or `white_charcoal` if an adjacent heatmap already uses red.
     - **Default Neutral / Non-significant**: `#D9D9D9` (NS/background), `#595959` (Vehicle/Control charcoal).
   - Keep chromatic ink restrained (~5–10% of figure area; desaturated HSL `S ≈ 0.44–0.63`, `L ≈ 0.57–0.71`), and lock global entities in a shared project `palette.json` across Fig. 1–6.
   - Consult [color-architecture.md](./references/color-architecture.md) for lineage hue-families and specialized alternatives.

6. **Concise CNS Axis Labels, Smart Label Tilting & Anti-Flatness Visual Craft**:
   - Consult [axis-labels-and-annotations.md](./references/axis-labels-and-annotations.md), [chart-selector.md](./references/chart-selector.md), and [micro-aesthetics.md](./references/micro-aesthetics.md) on every figure:
   - **Axis labels**: **1–3 words** (max 4), **≤ 22 chars**, **Sentence case**, units in `()` with one leading space, **ZERO `snake_case`**; **delete self-explanatory x-axis titles** on group box/bar/dot/heatmaps; never repeat panel title on y-axis.
   - **Smart X-Axis Label Tilting (Space Efficiency — measured on 504 CNS axes)**: ALWAYS call `smart_rotate_xticklabels(ax)` (Python) / `smart_x_labels()` (R) on categorical axes (`R = max_label_w_mm / pitch_mm`): `R < 1.15` → `0°`; `1.15 ≤ R < 3.80` → **`+45°` upward-right anchored** (`ha="right", va="top", rotation_mode="anchor"`, used in 96.9% of CNS tilted axes; NEVER `-45°`!); `R ≥ 3.80` (dense DotPlot/Heatmap) → **`+90°` bottom-to-top** (94.7% of CNS vertical axes); `max_len > 22` (long pathways) → **flip horizontally (`barh`) instead of tilting**!
   - **Clean Gridless White Background (Mandatory Default)**: ALL Cartesian plots (Bar, Box, Violin, Raincloud, Lollipop, Diverging Bar, Paired Slope, DotPlot, Spatial Boundary Curve, Ripley's *K*, Milo Beeswarm, Swimmer, SHAP, Scatter, Volcano, KM, Forest) MUST default to a **pure-white, zero-gridline L-axis canvas** (`axes.grid=False` / `theme_classic()`, top & right spines removed; `show_grid=False` by default in every helper). Use only single crisp dashed threshold/baseline reference lines (`x = 0`, `HR = 1`, `−30%/+20%` RECIST, tumor boundary `0 μm`) where statistically meaningful—NEVER background coordinate meshes (only polar Radar charts keep concentric rings).
   - **Dual-Track Chart Craft (Conventional Polished vs. Eye-Catching Upgrades)**: Never settle for flat default bar/box charts! Consult [chart-selector.md](./references/chart-selector.md), [single-cell-and-repertoire.md](./references/single-cell-and-repertoire.md), and [clinical-and-oncology-charts.md](./references/clinical-and-oncology-charts.md) to choose between **Track A (Polished Conventional)** (`plot_enhanced_bar` with `outlined_pastel` fill + `0.85 pt` border + white-rimmed dots on a clean gridless L-axis) and **Track B (Eye-Catching Tumor-Omics Upgrades)**:
     - *Enrichment & Stats*: `plot_lollipop()` (`×3.63`), `plot_diverging_enrichment()` (opposite-embedded labels), `plot_raincloud()`, `plot_split_corr_matrix()`, `plot_upset()`, `plot_radar_profile()`.
     - *Single-Cell & Repertoire*: `plot_cloudy_umap()`, `plot_atlas_umap_composite()`, `plot_marker_dotplot()` (grouped module brackets + dark dot rims), `plot_split_violin()`, `plot_milo_beeswarm()`, `plot_multi_cluster_volcano()`, `plot_roe_bubble_heatmap()`, `plot_pseudotime_cascade()`, `plot_chord_diagram()`, `plot_alluvial_flow()`.
     - *Clinical, Survival & ML*: `plot_km_survival()` (with exact aligned `Number at risk` table), `plot_cox_forest()` (zebra rows + log HR table), `plot_clinical_oncoprint()` (top TMB + alterations + right `%` + bottom clinical tracks), `plot_swimmer_plot()`, `plot_waterfall_response()`, `plot_spider_response()`, `plot_paired_slope()`, `plot_roc_curves()`, `plot_shap_beeswarm()`.
   - **Spatial Omics Visualization Standard** (empirically calibrated from **1,035 tumor spatial-omics papers + 588 frontier 2025–2026 spatial papers** — see [spatial-omics-visualization.md](./references/spatial-omics-visualization.md)): ALWAYS lock `1:1` physical aspect ratio, remove pixel coordinate axes, stamp a true vector scale bar (`add_spatial_scalebar`: `500 μm–1 mm` macro / `50–100 μm` ROI zoom), avoid opaque spots hiding H&E (pair clean H&E side-by-side or contour overlay), use **Focal-Lineage Spotlight** (muted `#E4E4E4` tissue footprint + 2–5 focal niches/lineages instead of rainbow confetti!), and leverage built-in core & 2025–2026 frontier spatial helpers (`plot_spatial_tissue_zoom`, `plot_spatial_spots`, `plot_spatial_boundary_gradient`, `plot_spatial_neighborhood_matrix`, `plot_spatial_lr_dual`, **`plot_spatial_streamline`** [COMMOT/Spateo signaling flux], **`plot_spatial_bivariate_map`** [3×3 2D color square], **`plot_spatial_delaunay_contacts`** [physical cell-contact mesh], **`plot_spatial_ripley_curve`** [Besag's `L(r) − r` + CSR envelope], **`plot_spatial_niche_network`** [CN border-adjacency graph], **`plot_spatial_3d_stack`** [2.5D isometric serial slices], and **`plot_spatial_kde_contours`**).
   - **Bar/Box/Point/Stroke micro-specs** (`CNS_SPEC` in `sci_plot_py.py` / `sci_plot_r.R`): bar width `0.60` with **zero bottom gap at `y = 0`**; hollow white-fill box width `0.52` (`showfliers=False`/`outlier.shape=NA` when jitter overlaid, `median_lw = 0.9 pt` > `box_lw = 0.55 pt` > `axis_lw = 0.40 pt` > `ref_lw = 0.35 pt`, data curves `1.05 pt`); tight jitter `±0.14`; ribbons `alpha = 0.18`; paired lines `0.35 pt, alpha = 0.30`; DotPlot dots with `0.3 pt` dark outline; Z-score heatmaps clipped symmetrically at `[-2, +2]`; 3–5 ticks per axis with true minus `−`; monotonic sorting for unordered categories; max 2–3 stacked brackets per panel.
   - **Biological italics**: human genes ALL-CAPS italic (*SPP1*, *TP53*), mouse genes Title-case italic (*Spp1*, *Trp53*), proteins upright (PD-1, CD8, IFN-*γ*), Latin phrases italic (*in vivo*, *in vitro*, *ex vivo*, *de novo*).

7. **Rigorous Statistical Significance Testing & Anti-Pseudoreplication**:
   - Consult [statistical-tests.md](./references/statistical-tests.md) and use `auto_stat_test()` / `auto_corr_test()` (`sci_plot_py.py` / `sci_plot_r.R`) whenever adding *p*-values or brackets:
   - **Always two-sided** by default; **match paired vs unpaired** structure (paired Wilcoxon signed-rank / paired *t*-test for matched Pre/Post or Tumor/Normal); **never test raw single-cell counts (`n = 50,000` cells) to claim patient-level differences** (aggregate to patient/sample-level pseudobulk or use LMM!).
   - **Test selection**: 2 groups → Mann–Whitney U (skewed/omics default) or Welch's *t*-test (normal/controlled wet-lab replicates); ≥3 groups → Kruskal–Wallis + BH-adjusted Dunn/pairwise Wilcoxon (or Welch ANOVA + Games–Howell/Tukey); categorical counts → two-sided Fisher's exact (expected < 5) or $\chi^2$; correlation → Spearman *ρ* (omics default) or Pearson *r* (linear normal); survival → Log-rank + Cox HR (95% CI).
   - **Always report** exact test name, two-sided, paired/unpaired, multiple-testing adjustment (BH FDR / Holm), and exact biological `n` in the accompanying caption/summary.

---

## Standard Workflow (Follow Step by Step)

### Step 1: Profile Data & Propose Dual-Track Chart Options (Interactive Default)
1. **Inspect the data structure first**: sample size `n` and biological unit (patients/mice vs raw single cells), number of groups/categories `K`, paired vs unpaired design, unidirectional vs bidirectional (`±`) values, secondary dimensions (gene count, FDR, clinical stage/response), and max category label length (`> 22 chars` → horizontal flip).
2. **Consult [chart-selector.md](./references/chart-selector.md)** to map the data + biological goal to the matching tumor-omics task.
3. **Interactive Dual-Track Proposal Rule (Default when the user provides data + goal WITHOUT naming a chart type)**:
   - Before writing plotting code, concisely propose **2–3 concrete options** tailored to the data profile:
     - **Option A — Conventional Polished (`Track A`)**: Top CNS workhorse with anti-flatness visual polish (e.g., `outlined_pastel` bars + white-rimmed dots, hollow boxplots on a crisp gridless pure-white L-axis).
     - **Option B (and optional C) — Eye-Catching High-Craft (`Track B`)**: 1–2 empirically rising Tier-1 CNS alternatives supported by the data's dimensions (e.g., Lollipop, Opposite-Embedded Diverging Enrichment, Raincloud, Multi-Cluster Track Volcano, $R_{o/e}$ Matrix, Paired Slope, Split-Triangle Correlation, RECIST Waterfall), marking the top recommendation and (for multi-panel rows) balancing **1–2 Track B focal panels + 2–3 clean Track A panels per 180 mm row**.
   - Let the user choose (`A` / `B` / `C` / `"compare all"`), then immediately execute Steps 2–5.
   - **Fast-path bypass**: If the user already explicitly requests a specific chart type (e.g., *"plot a lollipop"* / *"use the most eye-catching chart directly"*), skip the proposal and proceed straight to plotting.

### Step 2: Plan Physical Dimensions, Aspect Ratio & Full-Page Overflow Guard (`sci-paper-writing` Linkage)
Read [sizing-and-aspect-ratios.md](./references/sizing-and-aspect-ratios.md). Decide:
- Target journal full-page canvas (`180 × 220 mm` Nature / `174 × 215 mm` Cell / `184 × 218 mm` Science) and column fraction (`1/4`, `1/3`, `1/2`, `2/3`, `full`).
- Exact `(width_mm, height_mm)` for each panel (or consume exact dimensions from a shared `storyboard.json` generated by `sci-paper-writing`'s `plan-figures` CLI) so the plot box hits the empirical ideal aspect ratio.
- **Full-Page Overflow Check (Cross-Skill Handoff)**: When assembling a multi-panel Main Figure, NEVER shrink typography below `5 pt` just to cram >10 panels onto one page! If stacked rows exceed usable page height (`> 215–220 mm`), immediately (a) consolidate redundant panels via Track-B high-density helpers (`plot_multi_cluster_volcano`, `plot_split_violin`, `plot_diverging_enrichment`, `plot_roe_bubble_heatmap`, `plot_atlas_umap_composite`), OR (b) hand back to `sci-paper-writing` (`audit_manuscript.py plan-figures --rebalance`) to demote QC/gating/control/breakdown panels into paired Extended Data / Supplementary Figures or split across consecutive Main Figures before rendering!
- Typography preset (`cns_compact` by default unless user specifies `cell_standard` or `arial_8pt`).

### Step 3: Assign Colors (Global Consistency + Orthogonal Scales)
Read [color-architecture.md](./references/color-architecture.md):
- Check if the workspace already contains a `palette.json` / color mapping file; reuse exact hex codes for matching entities.
- Separate categorical groups, lineage hue-families, diverging Z-score/logFC scales (`RdBu`), and sequential expression scales (`LightGray -> Crimson/Indigo` or `Viridis`).

### Step 4: Generate Plot Using Reusable Helper Kit
- **Python**: import functions from [scripts/sci_plot_py.py](./scripts/sci_plot_py.py) (`apply_cns_style`, `mm_to_in`, `add_umap_corner_axes`, `rasterize_dense_artists`, `save_ai_figure`, `PALETTES`, `verify_ai_pdf`).
- **R**: source [scripts/sci_plot_r.R](./scripts/sci_plot_r.R) (`theme_cns`, `pt_to_lwd`, `pt_to_ggtext`, `umap_corner_axes`, `save_ai_pdf`, `CNS_PALETTES`).
- Refer to [hybrid-vector-raster-ai.md](./references/hybrid-vector-raster-ai.md) for edge cases (Scanpy, Seaborn clustermap, ComplexHeatmap, Seurat, `ggrepel`, rotated labels, broken axes).

### Step 5: Mandatory Programmatic & Visual Verification
Before finishing:
1. Export both editable **`.pdf`** (primary AI deliverable; optionally `.svg`) and a **200–300 dpi `.png`** preview in the project output directory.
2. Run `verify_ai_pdf("path/to/fig.pdf")` (or inspect with PyMuPDF) to confirm:
   - Exact physical `width_mm × height_mm`
   - Live embedded `Arial` / `Helvetica` vector text spans (count > 0, within target pt range)
   - `ai_single_line_text_pass == True` (`fragmented_td_blocks == 0` AND `kerning_split_tj_blocks == 0`, guaranteeing zero broken multi-element labels when opened in Illustrator)
   - Embedded high-DPI raster image present ONLY when dense layer (`> 2,000` pts or `> 2,500` cells) was rasterized
   - Reasonable PDF file size (typically `< 2 MB` even for 100k-cell UMAPs)
3. Visually inspect the rendered PNG preview using `read_file` to verify zero overlapping text, clean whitespace, legible ticks, and balanced aspect ratio. Fix any crowding before delivering.
