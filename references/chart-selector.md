# Data-to-Chart Decision Matrix: Conventional vs. Eye-Catching Upgrades (Tumor-Omics Edition)

Empirically calibrated from **4,968 parsed post-2020 high-IF tumor single-cell / spatial / multi-omics research articles** (`38,306` figure captions, `141,092` subpanel chunks across 28 flagship journals including *Nature*, *Cell*, *Science*, *Cancer Cell*, *Nat Med*, *Nat Cancer*, *Nat Genet*, *Nat Immunol*, *Immunity*, *Cancer Discov*, *Cell Rep Med*, *Nat Commun*) plus 172 recipes from `scidraw.msttnote.cn`.

> **Design Philosophy (告别“平淡” — Dual-Track Selection)**:
> Every scientific data type supports **two valid CNS publication tracks**:
> - **Track A — 常规稳健主力 (`Conventional Workhorse`)**: Highest reviewer familiarity — MUST apply our **Visual Craft Polish** (outlined-pastel fill + white-rimmed dots, subtle `#EFEFEF` grid, opposite-embedded labels, hollow boxplots) so it never looks flat or generic.
> - **Track B — 小众吸睛高阶 (`Eye-Catching / High-Craft Upgrade`)**: Empirically validated, rapidly rising in 2024–2026 Tier-1 CNS papers (`×1.6` to `×3.6` growth!). Use 1–2 eye-catching focal panels per main figure to elevate visual storytelling without sacrificing clarity.

---

## 1. Empirical Tumor-Omics Data Type × Chart Matrix (`n = 4,968` papers)

| Tumor-Omics Data Type | Track A: 常规稳健主力 (`Conventional` + Craft Polish) | Track B: 小众吸睛高阶 (`Eye-Catching` — Empirical 2020→2026 Trend) | Built-In Helper (`sci_plot_py` / `sci_plot_r`) | Rasterize? |
| :--- | :--- | :--- | :--- | :---: |
| **1. Pathway / GSEA Enrichment** (GO, KEGG, Hallmark, Reactome) | **Sorted Horizontal Bar** (`8.8%`) or **Enrichment Bubble/DotPlot** (`6.8%`, `×1.66` growth) — *upgrade*: trim `HALLMARK_` prefixes, add subtle vertical grid | **Lollipop / Cleveland Dot** (`×3.63` explosion! size = gene count, color = `−log₁₀(FDR)`), **Bidirectional Diverging Bar/Lollipop with Opposite-Side Embedded Labels** (*Nat Med* / *Cancer Cell* signature!), **Multi-Pathway GSEA Running Curve** (`4.8%` Tier-1), **EnrichmentMap / cnetplot Network** (`×1.62`) | `plot_lollipop()` / `plot_diverging_enrichment()` | No (vector) |
| **2. Single-Cell Atlas & Cell States** (`>2k` cells) | **Cloudy UMAP / t-SNE** (`59.5%` of papers; `37.6%` of atlas panels in Tier-1) with equal-length mini L-axes + direct labels (`K ≤ 12`) or circular major inset (`K ≥ 15`) | **Density Contour / 2D KDE UMAP** (`4.0%` Tier-1; gray full-atlas background + dark focal condition contours), **Milo Differential-Abundance Neighborhood Graph**, **Multi-Stage Small Multiples** | `plot_cloudy_umap()` / `plot_atlas_umap_composite()` | **Yes (600 dpi)** |
| **3. Cell Composition & Tissue Preference** (proportions across stages/tissues — see [single-cell-and-repertoire.md](./single-cell-and-repertoire.md)) | **100% Stacked Bar** (`×1.97` growth, `9.8%` Tier-1) or **Grouped Boxplot** (`13.6%` Tier-1) | **Milo Differential-Abundance Beeswarm** (anti-pseudoreplication gold standard!), **Zemin Zhang $R_{o/e}$ Tissue-Preference Matrix** (`×2.18` growth!), **Alluvial / Sankey Flow** (`×2.26` growth, `7.9%` Tier-1), **Ternary Plot** | `plot_milo_beeswarm()` / `plot_roe_bubble_heatmap()` / `plot_alluvial_flow()` | No (vector) |
| **4. Differential Expression & Cluster Markers** (`1k–30k` genes) | **Volcano Plot** (`28.9%`), **Clustered Heatmap** (`23.8%` Tier-1), or **Grouped Marker DotPlot** (`26.8%` Tier-1!) | **Grouped Marker DotPlot with Top Module Brackets & `0.3 pt` Dark Rims**, **Single-Cell Multi-Cluster Track Volcano** (`scRNAtoolVis` style), **Two-Condition Split-Violin**, **Fold-Change vs. Fold-Change 4-Quadrant Scatter** | `plot_marker_dotplot()` / `plot_multi_cluster_volcano()` / `plot_split_violin()` | **Yes (>2k pts, 600 dpi)** |
| **5. Clinical / Experimental Group Comparisons** (CR vs PR vs PD, Tumor vs Normal, scores) | **Hollow-Fill Boxplot + White-Rimmed Jitter** (`46.6%` Tier-1) or **Enhanced Outlined-Pastel Bar + Dots** (`42.6%` Tier-1) | **CNS Raincloud Plot** (half-violin density + slim hollow box + raindrop dots), **Paired Pre/Post Slope + Flanking Half-Boxes** (`15.9%` Tier-1, `×2.02` growth!), **Beeswarm / Sina + Median Bar**, **Dumbbell Plot** | `plot_enhanced_bar()` / `plot_raincloud()` / `plot_paired_slope()` | No (unless `>2k` cells) |
| **6. Genomics & Clinical Cohort Landscape** (WES/WGS/scDNA/Clinical Metadata — see [clinical-and-oncology-charts.md](./clinical-and-oncology-charts.md)) | **Multi-Track Clinical OncoPrint / Sample Inventory** (`11.1%` / `8.2%` Tier-1: Top TMB + Alterations + Right `%` + Bottom Stage/Response/MSI!) & **CNV Heatmap** (`8.4%`) | **Fishplot / Clonal Phylogeny Tree** (`5.4%` Tier-1!), **Mutation Lollipop Domain Track** (`3.7%`), **Circos SV/CNV Ring** (`7.2%` Tier-1) | `plot_clinical_oncoprint()` / `plot_lollipop()` | CNV >2.5k cells: **Yes** |
| **7. Cell–Cell Communication & TCR/BCR Repertoire** | **Ligand–Receptor DotPlot** (`12.6%`) or **Stacked Clonal Expansion Tier Bar** (`9.4%` Tier-1) | **Circos / Chord Diagram** (`12.6%` Tier-1, `×1.79` growth!), **Directed Signaling Network Graph**, **Alluvial / Sankey Clonal Sharing Flow** (`7.9%` Tier-1, `×2.26`) | `plot_chord_diagram()` / `plot_alluvial_flow()` | No (vector) |
| **8. Spatial Transcriptomics & mIF Neighborhoods** (Visium, Xenium, CosMx, CODEX, mIF — see [spatial-omics-visualization.md](./spatial-omics-visualization.md)) | **Paired H&E + Clean Hexagonal Spot / Single-Cell Map** (`54.5%` Tier-1!) + **Neighborhood Enrichment Heatmap** (`14.1%` Tier-1) + vector scale bar (`50–1000 μm`) | **Macro Tissue → Micro Single-Cell Voronoi/Polygon ROI Zoom** (`43.6%` Tier-1!), **Focal-Lineage Spotlight** (`×2.53` growth!), **Visium Pie-Chart Deconvolution** (`14.8%`), **Signed Boundary Distance Curve** (`10.4%`), **2025–26 Frontier**: Signaling Streamlines, 3×3 Bivariate Map, Delaunay Contact Graph, Ripley's `L(r) − r`, CN Adjacency Network, 2.5D Isometric 3D Stack | `plot_spatial_tissue_zoom()` / `plot_spatial_spots()` / `plot_spatial_boundary_gradient()` / `plot_spatial_streamline()` / `plot_spatial_bivariate_map()` / `plot_spatial_3d_stack()` | **Yes (600 dpi)** |
| **9. Survival, Clinical Outcomes & Treatment Response** | **Kaplan–Meier Curve + Exact Aligned `Number at Risk` Table** (`43.0%` Tier-1!) + **Tumor Growth Curve** | **Multivariable Cox Forest Plot** (`14.2%` Tier-1! zebra rows + left table + right log HR), **RECIST Waterfall** (`8.4%` Tier-1), **Oncology Swimmer Plot** (`2.9%` Tier-1), **Longitudinal Spider / ctDNA Plot** (`1.2%`) | `plot_km_survival()` / `plot_cox_forest()` / `plot_swimmer_plot()` / `plot_waterfall_response()` / `plot_spider_response()` | No (vector) |
| **10. Trajectory & Pseudotime Dynamics** | **UMAP + Vector Trajectory / RNA Velocity Arrows** (`33.7%` Tier-1) | **Smoothed Pseudotime Gene-Cascade Heatmap + Top State Ribbon** (`8.4%` Tier-1!), **Ridge / Joyplot along Pseudotime**, **Streamgraph Composition along Pseudotime** | `plot_pseudotime_cascade()` | Cell BG / large heatmap: **Yes** |
| **11. Correlation & Multi-Omics Association** | **Scatter + Regression + 95% CI Ribbon** (`15.8%` Tier-1) or **Correlation Heatmap** (`16.9%` Tier-1) | **Split-Triangle Correlation Matrix** (upper circles / lower exact *ρ* + stars / bold diagonal!), **Hallmark Multi-Axis Radar / Spider Profile** (`1.7%` Tier-1, `×1.77`), **Mantel Link Network** | `plot_split_corr_matrix()` / `plot_radar_profile()` | Scatter >2k pts: **Yes** |
| **12. ML Classifiers, Biomarkers, Overlaps & Perturbation** | **Multi-Model ROC / PR Curve + AUROC (95% CI)** (`7.5%` Tier-1, strict 1:1 square) + **CRISPR / Drug Sensitivity Ranked Plot** (`10.3%` Tier-1) | **SHAP Summary Beeswarm** (`×1.97` growth!), **UpSet Plot** for `≥ 4` sets (`3.7%` Tier-1, `×1.79` — always beat 4-circle Venn!) | `plot_roc_curves()` / `plot_shap_beeswarm()` / `plot_upset()` | No (vector) |

---

## 2. How to Upgrade "Flat / Plain" Conventional Charts (Anti-Flatness Playbook)

When sticking with a conventional chart (Track A), apply these 4 high-craft upgrades immediately so it looks *Nature*/*Cancer Cell* publication-grade rather than default matplotlib/ggplot2:

1. **Upgrade Plain Bar Plots (`plot_enhanced_bar` / `plot_diverging_enrichment`)**:
   - **Never** draw thick opaque rectangles with zero point overlay and blank flat background.
   - Use **`outlined_pastel`** (`alpha = 0.28` tinted fill + `0.85 pt` crisp same-hue outline) + **white-rimmed jittered replicate dots** (`edgecolor="white", lw=0.3 pt`) + subtle `#EFEFEF` (`0.30 pt`) horizontal guide gridlines + zero bottom gap at `y = 0`.
   - For bidirectional pathway/GSEA bars, **embed pathway labels on the empty opposite side of `x = 0`** (`plot_diverging_enrichment`) — cuts 40% wasted margin and immediately looks editorial!
   - Or swap bars outright for a **Lollipop plot (`plot_lollipop`)** encoding a 2nd metric (gene count or `−log₁₀(FDR)`) in terminal circle size/color!
2. **Upgrade Plain Boxplots (`plot_raincloud` / `plot_paired_slope`)**:
   - For `n ≥ 15` per group (clinical cohorts, signature scores), use a **Raincloud plot** (half-violin density cloud + ultra-slim hollow boxplot + left-side white-rimmed raindrop dots).
   - For matched Pre/Post or Normal/Tumor pairs, use **Paired Slope + Flanking Half-Boxes** with connecting lines tinted subtly by trajectory direction (coral-red increase vs steel-blue decrease).
3. **Upgrade Dense Multi-Cluster Marker Figures (`plot_multi_cluster_volcano` / `plot_roe_bubble_heatmap`)**:
   - Instead of 8 separate tiny volcano plots, draw **ONE Multi-Cluster Track Volcano** with colored cluster ribbon tiles along `y = 0` and italic top marker callouts.
   - Instead of cluttered stacked bars for tissue enrichment, draw a **Zemin Zhang $R_{o/e}$ tile/bubble matrix** with `+++/++/+` overlays.
4. **Apply Smart Space-Saving Label Tilting (`smart_rotate_xticklabels` / `smart_x_labels`)**:
   - Whenever horizontal x-labels crowd (`overflow ratio R ≥ 1.15`), tilt cleanly to **`+45°` right-anchored** (or `+90°` for ultra-dense DotPlot/Heatmap gene columns) instead of shrinking fonts below 5 pt or letting labels collide!

---

## 3. Anti-Pitfall Quick Rules

- **Never use pie charts for > 4 categories** — use a sorted horizontal lollipop/bar, 100% stacked bar, or $R_{o/e}$ matrix.
- **Never use Venn diagrams for ≥ 4 sets** — switch immediately to an **UpSet plot** (`×1.79` growth in Tier-1 CNS).
- **Never use bar charts without overlaid individual points** for biological replicates (`n ≤ 30`).
- **Don't over-decorate every panel**: Pair **1–2 eye-catching Track B focal panels** (e.g., Cloudy UMAP, Multi-Cluster Volcano, Raincloud, Diverging Lollipop, Split-Triangle Corr) with **2–3 clean Track A supporting panels** per 180 mm figure row for balanced CNS visual hierarchy.
