# Physical Sizing, Aspect Ratios & Compact Multi-Panel Layout

Empirically calibrated from 846 main figures and 3,859 subpanel slots across 167 high-IF research papers + official *Nature*, *Cell Press*, and *Science* artwork guides.

## 1. Journal Canvas Grids (1:1 Print Scale)

| Journal Family | 1-Column | 1.5-Column | Full Width (2-Col / 3-Col) | Safe Max Height | Panel Labels |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Nature Portfolio** (*Nature*, *Nat Med*, *Nat Genet*, *Nat Cancer*, *Nat Immunol*, *Nat Commun*) | 89 mm | 120–136 mm | **180 mm** (max 183 mm) | 225 mm (max 247 mm) | **8 pt bold** `a, b, c` (subjournals: 9 pt) |
| **Cell Press** (*Cell*, *Cancer Cell*, *Immunity*, *Cell Metab*) | 85 mm | 114 mm | **174 mm** | 220 mm (max 230 mm) | **9 pt bold** `A, B, C` |
| **Science (AAAS)** (*Science*, *Sci Immunol*, *Sci Transl Med*) | 57 mm | 121 mm | **184 mm** | 220 mm | **10 pt bold** `A, B, C` |

> **Golden Conversion**: `inches = mm / 25.4`. Always pass exact `mm` converted to inches at export time so panels assemble in Adobe Illustrator at **100% scale** with zero distortion of font size or stroke width.

### Row Slot Budgets (180 mm Full Width, 2.5–3 mm Gutter)
- **2 panels / row**: `88 mm` each (or asymmetric `55 + 122 mm`, `70 + 107 mm`)
- **3 panels / row**: `58 mm` each (most common CNS row layout!)
- **4 panels / row**: `43 mm` each (standard for boxplots, bars, compact KM, volcano, scatters)
- **5 panels / row**: `34 mm` each (ultra-compact 2-group boxplots/bars)
- *(For Cell 174 mm, subtract ~1.5 mm per column: 85 / 56 / 41.5 / 32.5 mm.)*

---

## 2. Empirical Subpanel Dimensions & Aspect Ratios (`W / H`)

Use this lookup table to set exact `(width_mm, height_mm)` when generating a single standalone panel for Illustrator assembly:

| Chart Type | Recommended Single-Panel `W × H` (mm) | Plot-Box Aspect (`W:H`) | Empirical CNS Median (IQR) | Sizing Formula / Notes |
| :--- | :---: | :---: | :---: | :--- |
| **Boxplot / Bar (2 groups)** | **`30 × 40`** to **`36 × 44`** | `0.65–0.85 : 1` (tall) | W 37–45, H 41–57 | `W ≈ 14 + 7 × n_groups` mm; fits 4–5/row! |
| **Boxplot / Violin / Bar (3–5 groups)** | **`42 × 42`** to **`54 × 45`** | `0.9–1.15 : 1` | W 37–45, H 41–57 | Rotate x-labels 45° or wrap to 2 lines |
| **Stacked Bar (composition, 3–8 bars)** | **`42 × 46`** to **`58 × 46`** | `0.85–1.25 : 1` | W 45, H 57 | Horizontal 100% bar: `120–175 × 24–32 mm` |
| **Volcano / MA Plot** | **`46 × 44`** to **`58 × 50`** | `1.0–1.15 : 1` | W 54 (43–67), H 43 | Square plot box + compact top legend; 3–4/row |
| **Kaplan-Meier Survival (curve only)** | **`50 × 40`** to **`58 × 42`** | `1.25–1.4 : 1` | W 56 (45–85), H 39 | Legend inside upper-right; `P` lower-left |
| **Kaplan-Meier + Risk Table** | **`44 × 52`** to **`56 × 54`** | `0.85–1.05 : 1` | W 44–56, H 50–55 | Fits **4/row** across 180 mm (42 mm each!) |
| **Line / Tumor Growth / Time Course** | **`54 × 40`** to **`65 × 44`** | `1.35–1.5 : 1` | W 53 (50–72), H 42 | Direct endpoint labels or shared top legend |
| **Scatter / Correlation / ROC** | **`45 × 44`** to **`56 × 52`** | `1.0 : 1` (strict square) | W 71 (46–100), H 43 | Tuck `r`, `P`, `n` inside empty top-left/bottom-right |
| **UMAP / t-SNE (direct-labeled, mini L-axis)** | **`46 × 46`** to **`58 × 54`** | `1.0 : 1` (point cloud) | W 63 (52–86), H 46 | **3/row** easily! Small multiples: `24 × 24 mm` each |
| **UMAP / t-SNE (with bottom horizontal legend)** | **`56 × 56`** to **`60 × 62`** | `1.0 : 1` + 10 mm legend | W 56–60, H 58 | Put 2–3 cols of legend underneath, never on right |
| **GSEA / Pathway Enrichment** | **`48 × 40`** to **`72 × 46`** | `1.2–1.6 : 1` | W 42 (39–55), H 37 | Horizontal bar width depends on pathway label length (keep names ≤35 chars!) |
| **Marker DotPlot (`G` genes × `C` clusters)** | **Matrix-driven** (typ. `85–140 × 42–60`) | `1.5–2.8 : 1` | W 74 (63–105), H 39 | `W ≈ 24 + 2.8 × G` mm; `H ≈ 18 + 2.6 × C` mm |
| **Heatmap (`M` cols × `N` rows)** | **Matrix-driven** (typ. `65–180 × 45–80`) | `1.2–3.5 : 1` (or `0.7:1` tall) | W 81 (54–105), H 49 | Tall 15–25 gene heatmap: `38–48 × 55 mm`; 60–80 gene panel: `175 × 50 mm` |
| **Forest Plot (Cox HR)** | **`85 × 65`** to **`110 × 95`** | `0.9–1.3 : 1` | W 88, H 70 | `H ≈ 16 + 3.8 × n_rows` mm; alternating zebra shading |
| **Hero Real-Image Strip (`hero_pathology_strip` / `hero_mif_codex_panel` / `hero_radiology_mri_ct`)** | **`174–180 × 65–92`** | `2.0–2.8 : 1` (**Full-Width Hero Band!**) | **35–48% of plate height!** | Multi-condition H&E/IHC/mIF overview + magnified ROI + channels, or Pre/Post clinical MRI/CT showcase |
| **2/3-Width Real Visual Hero (`pathology_he_ihc` / `mif_hero` / `radiology_mri_ct` / `gross_tumor_ivis`)** | **`115–122 × 58–82`** | `1.5–2.0 : 1` (**2/3 Row Hero!**) | Paired w/ 1/3 quant column | High-impact representative H&E/IHC/mIF/MRI/excised tumors with ruler; never shrink `< 80 mm` wide! |
| **1/2-Width Single-FOV Histology / Organoid / Blot (`histology` / `ihc` / `organoid_microscopy` / `western_blot`)** | **`86–90 × 52–62`** | `1.3–1.7 : 1` | **1/2 Row** | Minimum crisp readable scale bar & cellular detail (`600 dpi` raster + vector scale bar/kDa) |

---

## 3. Typography & Stroke Presets (Exact Real Print Units)

| Element | `cns_compact` **(Default)** | `cell_standard` | `arial_8pt` |
| :--- | :---: | :---: | :---: |
| Font family | `Arial` | `Arial` | `Arial` |
| Panel label (`a/A`) | **8 pt bold** (9 pt Cell) | **9 pt bold** | **10 pt bold** |
| Axis title / legend title / panel title | **7 pt** | **7.5 pt** | **8 pt** |
| Axis tick / category / legend item / stats | **6 pt** | **6.5–7 pt** | **7.5–8 pt** |
| Dense heatmap/dotplot gene names / risk table | **5–5.5 pt** | **5.5–6 pt** | **6 pt** |
| Axis spine & boxplot stroke width | **0.4 pt** | **0.45 pt** | **0.5 pt** |
| Tick length (`major.size`) | **2.0 pt** | **2.2 pt** | **2.5 pt** |
| Outer figure margin | **1.5–2.0 mm** | **2.0 mm** | **2.5 mm** |
| Legend swatch size | **2.2–2.5 mm** | **2.5 mm** | **3.0 mm** |

> **ggplot2 Unit Conversion Warning**:
> - `element_text(size = 6)` → **6 pt** directly.
> - `geom_text(size = 6 / ggplot2::.pt)` → **6 pt** (`ggplot2::.pt = 2.845`).
> - `element_line(linewidth = 0.4 / 2.13)` → **0.4 pt** (`linewidth = 0.19`! Default `0.5` is `1.07 pt` = 3× thicker than Nature!).

---

## 4. Six Anti-Crowding Mechanics (How CNS Fits 4 Panels in One Row)

1. **Ban right-side legends on 3–5 panel rows**: Move legend (a) directly onto clusters/lines, (b) inside an empty plot quadrant (`frameon=False`), (c) horizontally above/below in 2–4 columns, or (d) omit completely when a previous panel in the same figure defines the exact same color mapping.
2. **Use mini L-shaped corner axes on UMAP / t-SNE**: Replaces 4-sided box + numeric coordinates, recovering ~12 mm of width and height per panel.
3. **Strictly clean gridless white background (`show_grid = False` default)**: Remove top & right spines and **ALL background coordinate gridlines** (`axes.grid = False` / `theme_classic()` / `panel.grid = element_blank()`). Only keep (a) single crisp dashed reference/threshold lines (`y = 0`, `x = 0`, `HR = 1`, volcano cutoffs, tumor border `0 μm`) at `0.35–0.45 pt` gray, and (b) polar concentric rings on Radar charts. Never draw gray coordinate grids on Bar, Box, Violin, Raincloud, Lollipop, DotPlot, Swimmer, or Trajectory panels by default.
4. **Compact tick labels & strategic rotation**: Keep ≤5 numeric ticks per axis; rotate 3–5 word x-labels **45°** (`ha='right'`), dense gene names / >6 groups **90°**, or wrap 2-word labels onto 2 lines (`"Brain\nmetastases"`).
5. **Tight padding everywhere**: `axes.labelpad = 2.5 pt`, `xtick/ytick.major.pad = 2.0 pt`, `legend.handletextpad = 0.35`, `legend.columnspacing = 0.8`.
6. **Zero prose inside plots**: Put only concise labels, `n`, and exact `P`/`r`/`HR`. Move all interpretive descriptions to the caption.

---

## 5. Strict `1-Panel = 1-Unique-Letter` Indexing Law (`90.4%` CNS Standard)

Empirically audited across **`8,852` panel descriptions in `2,091` CNS Main Figures**: **`90.4%`** strictly assign **ONE unique panel letter (`a, b, c...` / `A, B, C...`) to ONE independent analytical plot**, while multi-chart bundling under a single letter occurs in only **`3.0%`**.

- **Hard Rule (Default)**: Every distinct analytical plot (`UMAP`, `boxplot`, `barplot`, `volcano`, `heatmap`, `KM`, `forest`, `dotplot`, `scatter`, `spatial map`) MUST receive its **own unique sequential top-left bold Arial letter** via `label_all_panels()` (Python) or `tag_panels_1to1()` (R).
- **Strictly Forbidden**:
  - NEVER let `Fig. 1a` cover a UMAP + proportion barplot (`a` = UMAP, `b` = barplot!).
  - NEVER let `Fig. 2b` cover 3 unrelated boxplots or 2 different statistical assays (`b, c, d` separately!).
  - NEVER bundle a representative microscopy photo + its statistical boxplot under one letter when space permits — **strongly prefer `a` (representative image) + `b` (cohort quantification)** so text references (`(Fig. 3a, b)`) and Illustrator assembly stay modular.
- **The ONLY 3 Narrow Allowed Exceptions**:
  1. **Overview + Magnified ROI Zoom / Split Fluorescence Channels** of the **exact same tissue FOV** (connected by a crisp ROI box or labeled `DAPI / CD8 / SPP1 / Merge`).
  2. **Identical-Assay Small-Multiple Strip** sharing identical axes, chart type, and a single shared colorbar/scale bar (e.g., 3 FeaturePlots on the same UMAP, or `Vehicle vs Treated` H&E/MRI images of identical scale).
  3. **Structurally Inseparable Composite Widget** (`KM + aligned Number-at-risk table`, `OncoPrint + top TMB + clinical tracks`, `plot_atlas_umap_composite` with major-lineage inset).

---

## 6. Visual Hero Zone Layouts for Real Imaging / Pathology / Radiology (`35–55%` Area Reservation)

Across `2,091` CNS Main Figures (`65.1%` of manuscripts; `42.4%` of *Cancer Cell* and `39.4%` of *Nature* Main Figures), real visual evidence (**MRI/CT/PET, H&E, IHC, mIF/CODEX/IMC, gross excised tumors/IVIS, organoids**) serves as the primary visual anchor ("seeing is believing") and routinely occupies **`35–55%` of the Main Figure canvas**.

Use one of these **3 CNS Hero Grid Blueprints** (`create_visual_hero_layout()` in Python / `compose_visual_hero_r()` in R):

1. **Blueprint A — Full-Width Top Hero Band (`40–48%` height) + Quantitative Rows Below** *(Best for multi-group mIF/CODEX, H&E overview + zooms, or longitudinal Pre/Post MRI/CT)*:
   - Row 1 (`180 × 75–90 mm`): **`a`** Full-width real visual showcase (e.g., Pre vs Post MRI + H&E + magnified boundary zoom, OR 3-group mIF + zoom).
   - Rows 2–3 (`180 × 115–130 mm`): **`b, c, d, e, f, g`** Crisp 1-letter-each statistical validations, spatial boundary gradients, and KM survival.
2. **Blueprint B — Asymmetric Left `2/3` Visual Hero + Right `1/3` Stacked Quant** *(Best for high-impact single-cohort tissue showcase paired immediately with its exact statistical proof)*:
   - Top block (`180 × 85–95 mm`): Left `118 × 90 mm` **`a`** (Large H&E/mIF/gross-tumor hero image) | Right `58 × 43 mm` **`b`** (Quantification 1) stacked over `58 × 43 mm` **`c`** (Quantification 2).
   - Bottom rows (`180 × 110 mm`): **`d, e, f, g`** Downstream mechanistic & clinical validation panels.
3. **Blueprint C — Twin Visual Anchor Plate** *(Best for *in vivo* therapy / clinical translation figures)*:
   - Row 1 (`180 × 64 mm`): **`a`** Gross excised tumors with ruler / IVIS bioluminescence (`115 × 64 mm`) + **`b`** Tumor growth curve (`60 × 64 mm`).
   - Row 2 (`180 × 68 mm`): **`c`** Representative IHC/mIF of Control vs Treated (`115 × 68 mm`) + **`d`** Quantification boxplot (`60 × 68 mm`).
   - Row 3 (`180 × 60 mm`): **`e, f, g`** Flow/immune infiltration + survival curves.

---

## 7. Anti-Sparse Plate Calibration, High-Substance Density & Anti-Half-Empty ED Rules

Crowded plates look unreadable, **but overly sparse plates (`> 76%` pure white or half-empty `4`-panel Extended Data plates) look thin and less substantial than real flagship CNS figures.** Balance BOTH bounds (full details in [journal-figure-standards.md](./journal-figure-standards.md)):

| Density / Layout Metric | Multi-Modal Main Fig Target | Pure Stat/Clinical Main Fig | Extended Data / Supp Fig Target | Automated Red Flag (`--audit-plate`) |
| :--- | :---: | :---: | :---: | :--- |
| **Pure-White Pixel Ratio (`RGB > 246`)** | **`52–72%`** | **`≤ 75%`** | **`≤ 78%`** | **`> 76%`** (Main) / **`> 80%`** (ED/Supp) |
| **Page Height Utilization (`H / H_max`)** | **`82–98%`** (`180–218 mm`) | **`82–98%`** | **`75–98%`** (`165–215 mm`) | **`< 65%`** (half-empty plate!) |
| **Panels per Extended Data / Supp Plate** | `8–10` (Main) | `8–10` (Main) | **`6–8` panels** | **`≤ 4` panels** + `< 65%` height → **merge twin ED plates!** |
| **Internal Horizontal White Void** | **`< 8 mm`** | **`< 8 mm`** | **`< 8 mm`** | **`≥ 10 mm`** empty band between legend/title & image! |
| **Visual / High-Density Anchor per Main Fig** | **≥ 1 real image Hero** | **≥ 1 rich matrix/OncoPrint/KM+Risk** | Paired visual/QC proof | `0` real visuals AND `0` complex matrices (`all scatter/bar`) |

### Five Quick Fixes When a Plate Looks Too Empty or "PPT-Like"
1. **Upgrade Bare Charts**: Always add **95% CI bands + aligned `Number at risk` table + `HR (95% CI)`** to every KM survival curve (`plot_km_survival`); never stretch a 3–4 bar univariate chart across `≥ 58 mm`; consolidate `≥ 3` tiny scatter plots into faceted 2D-contour strips or a correlation/forest summary.
2. **Promote Buried Real Visuals & Audit Unused Assets**: Check `100%` of extracted microscopy / MRI / H&E / spatial ROI files against assembled plates so zero real images are accidentally dropped or buried in ED while a Main Fig sits `80%` white.
3. **Merge Twin Half-Empty `4`-Panel ED Plates**: Combine two thematically related `4`-panel (`115–125 mm`, `83–90%` white) ED plates into ONE full-height `7–8` panel plate (`185–215 mm`), freeing slots for single-channel mIF splits, antibody controls, mouse toxicity/organ H&E, or uncropped blots.
4. **Ban Plain Code Text-Box Schematics**: Never render Study Designs (`Fig. 1a`) or Final Working Models as matplotlib rounded rectangles + straight arrows — draw biological cell/organ/receptor illustrations in Illustrator / BioRender (or reserve clean vector frames with `add_visual_placeholder_frame()`).
5. **Purge In-Figure AI / Notebook Labels**: Never leave `audit`, `crosswalk`, `parsimonious`, `habitat`, `OOF`, `Quarantined`, `phase cliff`, `Adversarial`, or pipeline cluster IDs (`c02_Fib_*`) inside plot titles or axis labels.
