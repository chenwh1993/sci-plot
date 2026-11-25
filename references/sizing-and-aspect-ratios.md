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
