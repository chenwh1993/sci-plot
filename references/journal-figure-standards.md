# Cross-Journal Artwork Format & Anti-Sparse Visual-Density Standards (`sci-plot`)

> **Purpose**: Two essential plotting guarantees in one reference:
> 1. **Exact Journal Artwork & Format Matrix** — never submit a `174 mm` Cell Press uppercase-`A` plate to *Nature* (`180 mm`, lowercase `a`), or vice versa!
> 2. **Anti-Sparse Plate & Visual-Density Calibration Laws** — prevent plots from looking "thin, empty, or PPT-like" compared with real flagship *Nature / Cell / Cancer Cell / Nature Cell Biology* plates (`n = 1,215` published figure architectures).

---

## Part 1 — Cross-Journal Artwork, Canvas, Typography & Figure-Architecture Matrix

| Artwork / Layout Standard | **Nature** *(Flagship)* | **Nature Cell Biology / Nat Cancer / Nat Immunol** | **Nature Medicine** | **Cell / Cancer Cell / Immunity** *(Cell Press)* | **Science / Sci Transl Med / Sci Adv** | **Cancer Discov / JCI / Clin Cancer Res** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Full-Page Width (`2-col`)** | **`180 mm`** *(max 183)* | **`180 mm`** | **`180 mm`** | **`174 mm`** | **`184 mm`** | **`178 mm`** |
| **1-Column / 1.5-Column Width** | `89 mm` / `120 mm` | `89 mm` / `120 mm` | `89 mm` / `120 mm` | `85 mm` / `114 mm` | `57 mm` / `121 mm` | `85 mm` / `118 mm` |
| **Usable Plate Height (`H_max`)** | **`220 mm`** *(full 247)* | **`220 mm`** | **`215 mm`** | **`215 mm`** *(full 230)* | **`218 mm`** *(full 235)* | **`210 mm`** |
| **Main Figure Count** | **`5`** (`4–5`, hard `6`) | **`7`** (`6–7`) | **`5`** (`4–6`) | **`7`** (**strict cap `7`**) | **`6`** (`5–6`; *STM* `6–7`) | **`6–7`** (*CCR* `5`) |
| **Ideal Panels / Main Fig** | **`8`** (`6–10`) | **`8–10`** (`6–11`) | **`6`** (`5–8`) | **`9–10`** (`7–11`) | **`9–10`** (`7–11`) | **`8–9`** (*CCR* `5–7`) |
| **Extended Data / Supp Architecture** | **`8–10` Extended Data Figs** (`6–8` panels, full-page!) + `~3` Supp | **`8–10` Extended Data Figs** (`6–8` panels, **`170–215 mm`**!) | **`7–10` ED Figs** (`5–6` panels) + Supp | **Strict `1:1` `Figure S1–S7`** (`7–8` panels, `174×215 mm`!) | **Consolidated Supp PDF** (`9–14` figs, `5–6` panels) | **`8–11` Supp Figs** (`4–6` panels) |
| **Panel Letter Style & Size** | **`a, b, c`** (`8 pt` bold Arial) | **`a, b, c`** (`8–9 pt` bold Arial) | **`a, b, c`** (`8 pt` bold Arial) | **`A, B, C`** (`9 pt` bold Arial) | **`A, B, C`** (`10 pt` bold Arial) | **`A, B, C`** (`9 pt` bold Arial) |
| **Typography Preset** | `cns_compact` (`6–7 pt`, floor `5 pt`) | `cns_compact` (`6–7 pt`, floor `5 pt`) | `cns_compact` (`6–7 pt`) | `cell_standard` (`6.5–7.5 pt`, floor `5.5 pt`) | `cns_compact` (`6–7 pt`) | `cell_standard` (`6.5–7.5 pt`) |
| **Spelling Inside Plot Labels** | **UK** (`Tumour`, `Haematoxylin`, `Oestrogen`) | **UK** (*NCB* strictly UK) | **US** (`Tumor`) | **US** (`Tumor`, `Hematoxylin`) | **US** (`Tumor`) | **US** (`Tumor`) |
| **Color Mode & Raster DPI** | RGB (web/AI), **≥300 dpi** (600 dpi dense raster) | RGB, **≥300–600 dpi** | RGB, **≥300–600 dpi** | RGB, **≥300–600 dpi** | RGB, **≥300–600 dpi** | RGB / CMYK at proof, **≥300–600 dpi** |
| **Mandatory Scale Bar / Source Data** | Exact numeric scale bar (`100 µm` / `1 cm`) + **Source Data** | Exact numeric scale bar + **Source Data** | Exact numeric scale bar + **Source Data** | Exact numeric scale bar + **Graphical Abstract** | Exact numeric scale bar + **Data File S1** (*STM* for `n < 20`) | Exact numeric scale bar |

---

## Part 2 — The Six Anti-Sparse Plate & Visual-Density Calibration Laws

Even when a manuscript has `7` Main Figures and `60+` panels, plates still look "thinner" than real *Nature / Cell / Cancer Cell* figures whenever they suffer from **high white-pixel emptiness, missing visual anchors, repetitive bare univariate charts, half-empty 4-panel Extended Data plates, or PPT-style code text-box schematics**. Enforce all six laws before final assembly:

### Law 1 — Plate Whitespace & Height-Utilization Targets (Anti-All-White Stat-Dump Law)
- **Target Pure-White Pixel Ratio (`RGB > 246`)**:
  - **Multi-modal Main Figure (with real microscopy/MRI/H&E/spatial map/heatmap anchor)**: **`52–72%` white** (flag **`> 76%` white** as sparse!).
  - **Statistical / Clinical Main Figure**: **`≤ 75%` white** (achieve substance via rich multi-dimensional encodings, NOT clutter or gridlines!).
  - **Extended Data / Supplementary Figure**: **`≤ 78%` white** AND **`≥ 75%` full-page height (`165–218 mm`)** (flag **`> 80%` white** or **`< 65%` height** as half-empty!).
- **Every Main Figure Needs ≥ 1 High-Density Visual Anchor**: Never let an entire Main Figure become `8–10` bare univariate scatterplots, simple barplots, and 2-line KM curves (`78–82%` pure white). Anchor every Main Figure with at least ONE of:
  1. **Real visual imaging Hero zone (`35–55%` area)**: clinical MRI/CT, H&E, mIF/CODEX/IHC, Visium/Xenium tissue ROI, gross excised tumors, or Western blot; OR
  2. **High-density data visualization anchor**: multi-track clinical OncoPrint, clustered multi-block heatmap, 2D spatial boundary/streamline/bivariate map, grouped marker dotplot with module brackets, multi-cluster volcano, or biological vector illustration.

### Law 2 — Zero Dropped / Buried Real-Visual Assets Audit Gate
- Before finalizing any multi-figure set, cross-check **100% of extracted / available real imaging files** (`MRI`, `H&E`, `IHC`, `mIF/CODEX`, `mouse IF`, `spatial tissue ROIs`, `WB`, `gross tumor photos`) against the assembled plates:
  1. **Zero orphaned assets**: Never leave extracted microscopy/IF images unused on disk while an *in vivo* or spatial figure lacks tissue-level proof!
  2. **Promote flagship visual proof out of Extended Data**: If a Main Figure sits `> 76%` white with pure scatter/bar plots while its paired Extended Data holds a striking representative mIF or high-power spatial ROI proving that exact figure's core claim, **promote that visual anchor into the Main Figure** and demote a redundant scatter/bar panel to ED!

### Law 3 — Minimum Information-Density Standards per Chart (Anti-Bare-Plot Law)
- **Kaplan–Meier Survival (`plot_km_survival`)**: ALWAYS include **95% CI shaded ribbons + aligned `Number at risk` table + Cox `HR (95% CI)` + two-sided log-rank `p`**. **NEVER** render bare two-line KM curves without risk tables or CI bands in Main or flagship ED figures!
- **Scatter / Correlation Plots**: NEVER chain `≥ 3` consecutive tiny (`n < 40`) single-variable scatter plots across a row. Either add 2D density contours + marginal/subgroup cues, merge into a multi-cohort forest/lollipop summary, or pair with a correlation matrix.
- **Univariate Bar Charts**: NEVER let a 3–4 bar chart occupy `≥ 1/3` row width (`≥ 58 mm`) of mostly white space! Compress to a crisp `38–44 mm` compact dot-CI / lollipop plot or combine across cohorts/strata.
- **Forest Plots (`plot_cox_forest`)**: Target `≥ 6–8` informative rows or paired unadjusted/adjusted dual-markers with aligned numeric `HR/OR (95% CI)` + `p` table columns; avoid sparse 3–4 row forest plots stretched across half a page.

### Law 4 — Anti-Half-Empty Extended Data Consolidation (`6–8` Panels / Full-Page `170–215 mm`)
- In *Nature / NCB* (cap `10` ED Figs) and *Cell Press* (strict `7` S-Figs), **NEVER scatter `4`-panel half-height (`110–130 mm`, `82–90%` white) plates across multiple slots** when two thematically related ED figures (e.g., two 4-panel external validation plates both supporting `Fig. 6`) can be **merged into ONE substantial `7–8` panel full-page plate (`185–215 mm`)**!
- Consolidating twin 4-panel ED plates immediately (a) raises visual substance to the *Nature / NCB* median (`6–8` panels/ED), and (b) frees dedicated ED slots for high-value reviewer armor (single-channel mIF splits, antibody controls, mouse body weight & normal-organ H&E toxicity, uncropped Western blots).

### Law 5 — Zero Dead Internal Voids & Ban Plain Code Text-Box Schematics
- **Eliminate Internal White Voids (`≥ 10 mm` empty vertical bands)**: Never leave `15–25 mm` of empty white gap between a top legend/title and a microscopy strip, nor shrink an atlas UMAP to a tiny floating island surrounded by `35 mm` blank margins!
- **Ban Matplotlib/Python Text-Box Flowcharts for Study Designs (`Fig. 1a`) & Final Working Models (`Last Panel`)**: Plain rounded rectangles filled with bullet text + straight arrows look like software architecture slides, not *Nature / Cell* biology! Always craft Study Design schematics and Final Mechanistic Working Models as **biological cell/tissue/organ/receptor illustrations** (BioRender, Adobe Illustrator, or crisp biological SVG vector art — using `add_visual_placeholder_frame()` if assembling vectors later in Illustrator).

### Law 6 — In-Figure Anti-AI & 1:1 Text↔Figure Terminology Hygiene
- Every panel title, axis label, and annotation inside a figure MUST pass the **exact same anti-AI and anti-notebook-jargon filter** as the manuscript text (`sci-paper-writing` `14-anti-ai-prose-stat-restraint-and-rule-harmony.md`):
  - **NEVER print internal pipeline/cohort/cluster codes** inside publication figures (`c02_Fib_POSTN`, `m01_Macro_SPP1`, `PipeV2_DropAll`, `MalignantGate`, `Primary40`, `ValSet15`, `OOF`, `Quarantined Val`).
  - **NEVER print non-CNS AI/tech-blog terms** inside figure titles (`audit`, `crosswalk`, `parsimonious`, `habitat`, `out-of-fold`, `phase cliff`, `Adversarial`, `triage cascade`).
  - Run `python scripts/sci_plot_py.py --audit-plate <assembled_dir_or_pdf> --journal nature_sister` before final delivery to automatically check **pixel whitespace %**, **page height utilization %**, **dead vertical white bands**, **in-figure AI/jargon text leaks**, and **orphaned visual assets**!
