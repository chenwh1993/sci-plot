# Spatial Omics Visualization Playbook (CNS Tumor-Omics Standard)

Empirically distilled from **1,035 parsed 2020–2026 high-IF tumor spatial-omics original research articles** (`365` Tier-1 CNS flagship papers across *Nature*, *Cell*, *Cancer Cell*, *Nat Med*, *Nat Cancer*, *Nat Genet*, *Nat Methods*, *Nat Biotechnol*, *Cancer Discov*, *Cell Rep Med*, *Nat Commun*), covering Visium, Visium HD, Xenium, CosMx, MERFISH, Stereo-seq, Slide-seq, CODEX/PhenoCycler, IMC, MIBI, GeoMx DSP, and mIF.

---

## 1. Seven Non-Negotiable Spatial Figure Rules

1. **Strict 1:1 Physical Aspect Ratio — Zero Pixel Coordinate Axes**:
   - ALWAYS lock `ax.set_aspect("equal")` (Python) / `coord_fixed(ratio = 1)` (R) so tissue anatomy is never stretched or squashed.
   - ALWAYS delete raw pixel ticks, axis numbers (`0, 5000, 10000 px`), and bounding boxes (`theme_spatial_clean()`).
2. **Mandatory Vector Physical Scale Bar (Empirically Calibrated by Field of View)**:
   - Replace coordinate axes with a crisp vector scale bar (`add_spatial_scalebar()`) calibrated to true μm/mm:
     - **Whole-slide macro slice** (Visium, Xenium, CODEX whole tissue): **`500 μm`**, **`1 mm`**, or **`2 mm`**
     - **Regional FOV / IMC / CosMx / medium zoom**: **`100 μm`** (`#1` most frequent across 1,035 papers!) or **`200 μm`**
     - **High-magnification cellular ROI inset**: **`50 μm`** (`#2` overall!) or **`20 μm`**
   - Place in lower-right or lower-left corner (`#1A1A1A` bar on white canvas; `#FFFFFF` bar on dark-field fluorescence or over a subtle translucent dark pill on busy H&E).
3. **Never Bury Histology Under Opaque Spots (Paired H&E Rule — `54.5%` of Tier-1 CNS)**:
   - Default `Seurat::SpatialDimPlot(alpha = 1)` hides H&E completely, while `alpha = 0.4` makes colors muddy.
   - **CNS Gold Standard**: Show **clean H&E (left)** side-by-side at 1:1 alignment with **clean spatial domain/expression map (right)**, OR overlay crisp **vector niche boundary contours** over softly desaturated H&E (`alpha = 0.85` H&E + border-only contours).
4. **Macro Tissue → Micro Cellular ROI Zoom-In (`43.6%` of Tier-1 CNS!)**:
   - Almost half of flagship spatial figures pair a **whole-tissue macro map** (with a crisp `0.65 pt` solid/dashed ROI box) with **1–2 magnified cellular ROI insets** (`50–100 μm` scale bar) showing true cell polygons/Voronoi tessellation or ligand–receptor contact.
   - Use `plot_spatial_tissue_zoom()` to generate both panels cleanly in one figure.
5. **Focal-Lineage Spotlight Over Muted Tissue Footprint (Beat "Rainbow Confetti"!)**:
   - Plotting 20+ colorful cell types simultaneously on one tissue slice creates unreadable visual noise.
   - **CNS Trick**: Render non-focal background cells/stroma in soft uniform gray (`#E4E4E4`) to preserve exact tissue silhouette, and color ONLY the **2–5 focal cell types or spatial niches** in desaturated CNS hues (`#C84B4B`, `#3B6EA5`, `#2A9D8F`, `#E08E45`, `#7A5195`)!
6. **Strict Light-Canvas vs. Dark-Field Palette Separation**:
   - **Light Canvas (`#FFFFFF` background)** — Default for Visium, Visium HD, Stereo-seq, Xenium/CosMx cell-type & niche maps, Voronoi maps, and continuous gene expression (`gray_crimson` `#EFEFEF → #67000D` or `gray_indigo` — NEVER `jet`/`rainbow`!).
   - **Dark-Field Canvas (`#0A0A0A` background)** — Reserved for **multiplex fluorescence micrographs** (mIF, CODEX/PhenoCycler, CyCIF, IMC, MIBI, RNA-FISH). Use **colorblind-safe orthogonal fluorophores** (NEVER red + green alone!):
     - Nuclei / DAPI: `#3B82F6` (dim blue, low brightness so markers pop)
     - Tumor / PanCK: `#E5E7EB` (crisp silver-white) or `#38BDF8` (cyan)
     - Focal Marker 1 (e.g., CD8): `#22C55E` (emerald) or `#F59E0B` (amber)
     - Focal Marker 2 (e.g., CXCL13 / PD-L1): `#F43F5E` (magenta-rose)
     - Focal Marker 3 (e.g., SPP1 / Macrophage): `#A855F7` (violet) or `#FACC15` (yellow)
7. **Hybrid 600-dpi Raster Tissue + 100% Vector Overlays (AI-Safe)**:
   - Always rasterize high-density cell/spot clouds (`> 2,000` cells/spots) and H&E/fluorescence bitmaps at **600 dpi**, while keeping **ROI zoom boxes, leader lines, tissue boundary contours, scale bars, colorbars, and Arial labels 100% editable vectors**!

---

## 2. Platform-to-Visualization Decision Matrix (`n = 1,035` Spatial Tumor Papers)

| Spatial Platform | Resolution & Geometry | Track A: 常规高稳展示 (`Polished Conventional`) | Track B: 小众吸睛高阶 (`Eye-Catching Tier-1 CNS`) | Built-In Helper |
| :--- | :--- | :--- | :--- | :--- |
| **10x Visium** (`19.3%`) | `55 μm` hexagonal spots (multi-cell) | Side-by-side **H&E + clean hexagonal spot map** (`gray_crimson` continuous or desaturated niche domains) | **Spot Pie-Chart Deconvolution Map** (SPOTlight/RCTD/Cell2location miniature pies per spot, `14.8%` Tier-1!), **Signed Boundary Distance Gradient Curve** (`10.4%` Tier-1), **Spatial Niche Domain Map** (`×2.53` growth!) | `plot_spatial_spots(mode="pie_deconv")` / `plot_spatial_boundary_gradient()` |
| **Visium HD / Stereo-seq / Slide-seq** (`9.0%`) | `2–16 μm` subcellular/near-single-cell bins | High-DPI rasterized bin map (`600 dpi`) + matched H&E + `500 μm` scale bar | **Macro Tissue → Micro ROI Zoom** (`43.6%` Tier-1) + **Focal-Lineage Spotlight** (muted `#E4E4E4` tissue + 3–4 focal niches) + **Spatial Gradient / Flow Vectors** (`7.4%` Tier-1!) | `plot_spatial_tissue_zoom()` |
| **Xenium / CosMx / MERFISH** (`24.4%`) | True single-cell `(x, y)` + cell segmentation polygons | Rasterized single-cell scatter (`600 dpi`) colored by lineage (`K ≤ 8`) + `500 μm` scale bar | **Clipped Voronoi / Single-Cell Polygon Tissue Map** (`13.7%` Tier-1, crisp `0.18 pt` cell borders inside magnified ROI!), **Dual-Channel Ligand–Receptor Colocalization Spotlight** (`13.4%` Tier-1), **Spatial Neighborhood Z-Score Matrix** (`×1.69` growth) | `plot_spatial_tissue_zoom()` / `plot_spatial_lr_dual()` / `plot_spatial_neighborhood_matrix()` |
| **CODEX / IMC / MIBI / mIF** (`33.7%`) | Multiplex protein imaging + segmented cells | **Dark-field multiplex composite (`#0A0A0A`)** paired with **digital cell-phenotype map** + channel strip + `100 μm` scale bar | **Cellular Neighborhood (CN) Voronoi Map** (`×2.53` growth!) + **Neighborhood Co-Enrichment / Avoidance Matrix** (`squidpy`/IMC style) + **Signed Invasive-Margin Distance Curve** (`−200 μm` Core ← `0` Margin → `+200 μm` Stroma) | `plot_spatial_neighborhood_matrix()` / `plot_spatial_boundary_gradient()` |
| **GeoMx DSP** (`7.1%`) | Compartmentalized AOIs (PanCK⁺ Tumor vs CD45⁺ Stroma) | Paired fluorescence ROI thumbnail + **Paired Pre/Post or Tumor–Stroma Slope Plot** (`plot_paired_slope`) | **Diverging Opposite-Embedded Pathway Plot** (Tumor vs Stroma AOIs) + **Multi-Compartment Raincloud** | `plot_paired_slope()` / `plot_diverging_enrichment()` |

---

## 3. Recommended 180 mm Spatial Main Figure Blueprint (Nature / Cancer Cell)

To tell a complete, non-flat spatial story in **one full-width 180 mm figure (`180 × 106 mm`)**:
- **Top row (Spatial Anatomy)**:
  - **a** (`~56 mm`): **Whole-Tissue Macro Map** (H&E or Focal-Lineage Spotlight at 600 dpi + dashed ROI box + `500 μm` scale bar)
  - **b** (`~56 mm`): **Magnified Single-Cell Voronoi / Polygon ROI Zoom** (exact cell-cell niche interface + `50 μm` scale bar)
  - **c** (`~58 mm`): **Visium Spot Deconvolution Pie Map** OR **Dual-Channel Ligand–Receptor Colocalization**
- **Bottom row (Quantitative Spatial Proof)**:
  - **d** (`~62 mm`): **Signed Distance to Tumor–Stroma Boundary Curve** (`−200` Core ← `0` Margin → `+200 μm` Stroma with `95% CI` ribbons!)
  - **e** (`~54 mm`): **Spatial Neighborhood Co-Enrichment / Avoidance Z-Score Matrix** (`[-3, +3]` diverging + `*`/`**`/`***` + `+45°` labels)
  - **f** (`~54 mm`): **Spatial Niche Clinical Comparison** (Raincloud or Paired Pre/Post Slope comparing Responder vs Non-Responder niche abundance!)

---

## 4. 2025–2026 Frontier Spatial Data Readouts & Eye-Catching Niche Visualizations

Empirically distilled from **`588` brand-new 2025–2026 high-IF tumor spatial-omics papers** (`194` Tier-1 CNS flagship papers). Beyond classic spot/cell maps, 2025–2026 papers present **12 specialized spatial data readouts** using these high-craft, eye-catching visualizations:

| 2025–2026 Spatial Data Readout (展示什么数据) | Tier-1 CNS % (`2025–26`) | Conventional Display (Track A) | Eye-Catching Niche Upgrade (Track B — 小众吸睛) | Built-In Helper |
| :--- | :---: | :--- | :--- | :--- |
| **1. Multi-Slice / 3D Tumor & TLS Depth Continuity** (serial Visium / Xenium / mIF sections) | **`11.9%`** | 3–4 separate flat 2D panels side by side | **2.5D Isometric Oblique Serial-Section Stack** (translucent frosted slice planes `z = 0..300 μm` + vertical alignment rails showing 3D tumor/TLS continuity — zero ugly `mplot3d` boxes!) | `plot_spatial_3d_stack()` |
| **2. Perivascular / Perineural / Hypoxic-Necrotic Micro-Axes** (distance to vessels, HEVs, nerves, necrosis) | **`9.3%`** | Bar plot comparing "Near" vs "Far" bins | **Signed Anatomical-Landmark Distance Curve** (`plot_spatial_boundary_gradient`) OR **Concentric Multi-Ring Annuli** (`0–25`, `25–50`, `50–100 μm` drawn on tissue + radial composition strip) | `plot_spatial_boundary_gradient()` |
| **3. Spatial Multi-Omics Co-Mapping** (Spatial RNA + MALDI/DESI Metabolomics, Spatial ATAC, Proteomics) | **`9.3%`** | Unaligned separate tissue plots | **Synchronized Twin-Slice Co-Registration** (identical tissue contour + orthogonal palettes: `gray_crimson` for RNA vs `gray_indigo` for metabolite/ATAC) OR **3×3 Bivariate Choropleth** | `plot_spatial_bivariate_map()` |
| **4. Multi-Scale Spatial Clustering vs. Exclusion Physics** (Ripley's *K* / Besag's *L*(*r*) / Cross-*K* / Moran's *I*) | **`8.2%`** | Single-radius nearest-neighbor boxplot | **Multi-Scale `L(r) − r` Point-Process Curve** (`r = 0..200 μm`) with **shaded CSR null-permutation envelope** (`> 0` Attraction vs `< 0` Exclusion across cohorts!) | `plot_spatial_ripley_curve()` |
| **5. Cellular Neighborhood (`CN`) Macro-Architecture & Borders** (which niches physically touch each other) | **`6.7%`** | Co-enrichment heatmap only | **CN Border-Adjacency Network Graph** (node size ∝ CN abundance, crimson edges = enriched physical border, steel-blue = avoided interface) paired with CN cell-type stacked composition | `plot_spatial_niche_network()` |
| **6. AI Pathology Foundation-Model Attention / Virtual ST** (H&E WSI tile attention weights, predicted expression) | **`6.7%`** | Coarse heatmap overlay | **Side-by-Side H&E + Smooth Tile Attention Overlay** (`gray_crimson` clipped `2nd–98th` pct) + **Per-Slide Moran's *I* vs. Prediction Scatter** | `plot_spatial_spots()` |
| **7. Subcellular Single-Molecule mRNA Polarity** (Xenium / CosMx / MERFISH / Visium HD `10–20 μm` zoom) | **`5.7%`** | Cell-level mean expression dot | **Segmented Cell-Membrane + Nucleus Polygon Map with Individual Single-Molecule Transcript Dots** (distinct high-contrast markers for nuclear vs cytoplasmic/synaptic transcripts) | `plot_spatial_tissue_zoom()` |
| **8. Dual Continuous Spatial Programs** (Hypoxia × IFN-γ, EMT × Glycolysis, Bivariate Moran's *R*) | **`5.2%`** | Two separate feature maps side by side | **3×3 Bivariate Choropleth Spatial Map** (blends 2 continuous programs into ONE tissue map with a crisp 3×3 2D square legend!) | `plot_spatial_bivariate_map()` |
| **9. Tumor–Immune Spatial Mixing / Entropy Archetypes** (Inflamed vs Excluded vs Desert, Shannon diversity) | **`5.2%`** | Categorical bar chart of subtypes | **Continuous Local Spatial Entropy Map** + **Cohort Spatial Archetype Barcode Strip** aligned with RECIST/OS survival tracks | `plot_spatial_kde_contours()` |
| **10. Spatial Signaling Flux & Differentiation Velocity** (`COMMOT`, `Spateo`, `Moscot`, `FlowSig`) | **`3.6%`** | Non-directional L–R dotplot | **Curved Vector Streamline Field (`streamplot`) Overlaid on Tissue** (directional signaling arrows flowing from sender niche to receiver niche across the invasive margin!) | `plot_spatial_streamline()` |
| **11. Spatial Genetic Subclones & Spatial TCR Clonotypes** (`CalicoST`, `Numbat`, Spatial V(D)J) | **`2.1%`** | Clonal barplot without spatial coordinates | **Spatial Clonal Territory Map + Paired Phylogenetic Tree** (exact color sync between phylogenetic branches and spatial subclone territories / infiltrating TCR clonotypes!) | `plot_spatial_tissue_zoom()` |
| **12. Physical Cell–Cell Contact Mesh & Density Hotspots** (`k`-NN / Delaunay physical contacts, KDE isocontours) | **`2.1%`** | Plain scatter plot | **Delaunay Physical Contact Graph** (`< 25–30 μm` edges: bold crimson heterotypic tumor–immune contacts over light-gray homotypic mesh!) & **Topographic KDE Isocontour Hotspots** over muted tissue | `plot_spatial_delaunay_contacts()` / `plot_spatial_kde_contours()` |
