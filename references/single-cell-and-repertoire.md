# Single-Cell, Immune Repertoire & Multi-Omics Visualization Playbook

Empirically calibrated from **747 Tier-1 CNS 2025–2026 tumor papers** (*Nature*, *Cell*, *Science*, *Cancer Cell*, *Nat Med*, *Nat Cancer*, *Nat Genet*, *Nat Immunol*, *Immunity*, *Cancer Discov*, *Cell Rep Med*).

---

## 1. Complete Single-Cell & Repertoire Data-to-Chart Matrix (`2025–2026 Tier-1 CNS`)

| Single-Cell / Repertoire Readout | Tier-1 % (`25–26`) | Track A: Polished Conventional | Track B: Eye-Catching Niche Upgrade | Built-In Helper |
| :--- | :---: | :--- | :--- | :--- |
| **1. Cell-State Atlas & Subclusters** (`10k–500k` cells) | `63.8%` | Zemin Zhang Cloudy UMAP (`K ≤ 12`: same-color direct labels; `K ≥ 15`: circular major inset + tight grouped vector legend) | **Multi-Condition Density Contour UMAP** (gray full atlas + dark focal contours) or **Pseudobulk PCA Trajectory** | `plot_cloudy_umap()` / `plot_atlas_umap_composite()` |
| **2. Multi-Cluster Marker Expression** (`5–30` clusters × `15–50` genes) | **`26.8%`** | Clustered Heatmap (`[-2, +2]` `rdbu_soft` clipped) | **Grouped Marker DotPlot with Top Module Brackets** (`0.3 pt` dark dot rims so pale dots stay sharp + top color-coded lineage/module bars + italic `+45°/+90°` gene ticks!) | `plot_marker_dotplot()` |
| **3. Two-Condition Expression Across Multiple Clusters** | `15.9%` | Dodged boxplots across facets | **Compact Split-Violin** (left half = Group A, right half = Group B + internal median/IQR bar — cuts panel width by 50%!) | `plot_split_violin()` |
| **4. Differential Abundance / Subtype Enrichment** | `13.6%` | Patient-level pseudobulk boxplot or 100% stacked bar | **Milo Differential-Abundance Beeswarm** (neighborhood log₂FC colored red/blue/NS — gold-standard anti-pseudoreplication!) OR **Zemin Zhang $R_{o/e}$ Matrix** (`+++/++/+`) | `plot_milo_beeswarm()` / `plot_roe_bubble_heatmap()` |
| **5. Multi-Cluster Differential Expression** (`4–10` clusters) | `14.0%` | Separate volcano plots per cluster | **Single-Cell Multi-Cluster Track Volcano** (colored cluster ribbon along `y = 0` + mirrored up/down gene clouds + italic top markers!) | `plot_multi_cluster_volcano()` |
| **6. Cell–Cell Communication (CellChat / CellPhoneDB / NicheNet)** | **`12.6%`** | Ligand–Receptor DotPlot (`size = −log₁₀P`, `color = prob`) | **Circular Chord / Circos Diagram** (`≤ 8` lineages, proportional outer arcs + smooth Bézier sender→receiver ribbons!) | `plot_chord_diagram()` |
| **7. TCR / BCR Clonal Expansion & Cross-Compartment Sharing** | **`9.4%`** | Stacked Clonal Expansion Tier Bar (`Hyperexpanded >100 \| Large \| Medium \| Small \| Single`) | **Alluvial / Sankey Clonal Transition Flow** (`PBMC → Adjacent → Primary → Metastasis` or `Pre-ICB → On-ICB`) + **STARTRAC Migration/Transition Scatter** | `plot_alluvial_flow()` |
| **8. Pseudotime & Differentiation Dynamics** (`Monocle 3 / Slingshot / scVelo`) | **`8.4%`** | UMAP colored by pseudotime + arrows | **Smoothed Pseudotime Gene-Cascade Heatmap** (genes sorted monotonically by peak timing + top cell-state ribbon + right italic TF callouts!) | `plot_pseudotime_cascade()` |
| **9. Single-Cell CNV & Malignant Subclones** (`InferCNV / CopyKAT / Numbat`) | **`8.4%`** | Chromosome-ordered CNV Heatmap (600-dpi raster body + chromosome arm alternating top strips) | **Clonal Evolution Fishplot / Phylogenetic Tree** synced to spatial or longitudinal sample colors | Raster `imshow` + `rdbu_soft` |
| **10. scATAC-seq / Multiome Chromatin Accessibility** | **`4.3%`** | Motif enrichment DotPlot / `chromVAR` heatmap | **Normalized Pseudo-Bulk Coverage Tracks + Peak-to-Gene Co-Accessibility Arcs** (aligned above gene body exon/intron track!) | Custom vector track stack |
| **11. CRISPR / Perturb-seq Functional Screens** | **`10.3%`** | Volcano / Ranked `MAGeCK` β-score plot | **Hockey-Stick Ranked Lollipop** (top sensitizers vs resistance hits labeled in italic) + **Perturbation-by-Module Heatmap** | `plot_lollipop()` |

---

## 2. Single-Cell Visual Craft Rules (Zero Pseudoreplication & Zero Clutter)

1. **Always Separate Cell-Level Visualizations from Patient-Level Inference**:
   - UMAP, DotPlot, Split-Violin, and Pseudotime Cascade display single-cell structure; **never** slap raw `n = 50,000 cells` Mann–Whitney *p*-values on patient-group comparisons — use **Milo neighborhoods (`plot_milo_beeswarm`)** or **patient-level pseudobulk (`plot_raincloud` / `plot_paired_slope`)**.
2. **DotPlot Readability Fix**:
   - Pale dots (`Z-score ≈ 0` or low expression) disappear against white backgrounds unless EVERY dot has a crisp **`0.28–0.35 pt` dark rim (`edgecolor="#333333"`)**.
   - Always group marker genes into biological modules with crisp colored top header strips (`Naive`, `Cytotoxic`, `Exhausted`, `Treg`) so readers decode 25 genes in 2 seconds.
3. **Chord & Alluvial Restraint**:
   - Limit Chord diagrams to **≤ 8 lineages** and keep only the top 15–25 strongest edges (`alpha = 0.48–0.65`) so ribbons stay sculptural rather than becoming an unreadable ball of yarn.
   - Limit Alluvial/Sankey flows to **2–4 columns** and **≤ 7 states per column** with clean right/left direct labels.
4. **Exact Physical-mm UMAP Corner L-Axes (`X == Y` Across Every Panel)**:
   - **Never** draw UMAP corner arrows in `axes fraction` (e.g., `frac = 0.17`) or raw data units (`max_span * 0.17`): axes-fraction arrows stretch into **unequal X vs Y lengths (`X ≠ Y`)** whenever a panel is non-square or reshaped post-hoc by `constrained_layout`, `tight_layout()`, `colorbar()`, or `set_aspect("equal")`, and scale inconsistently (5 mm vs 22 mm) between small subplots and full-width atlases.
   - **Always use physical-mm draw-time transforms (`arrow_mm = 8.0` default)**:
     - **Python**: `add_umap_corner_axes(ax, arrow_mm=8.0)` (built into `plot_cloudy_umap` & `plot_atlas_umap_composite`) anchors the origin in `ax.transAxes` and computes both arrowheads at draw time via `matplotlib.transforms.offset_copy(ax.transAxes, fig=fig, x=arrow_in, y=0, units="inches")` and `(0, arrow_in)` — guaranteeing exact `8.00 × 8.00 mm` X and Y lengths regardless of panel geometry or post-hoc colorbars.
     - **R**: `add_umap_corner_axes_r(arrow_mm = 8.0)` (built into `plot_cloudy_umap` & `plot_atlas_umap_composite`) draws `grid::segmentsGrob` + `grid::textGrob` with `grid::unit(arrow_mm, "mm")` inside `annotation_custom(..., xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)` — guaranteeing exact `8.00 × 8.00 mm` X and Y lengths even inside uneven `patchwork` layouts!
