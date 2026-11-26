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

### Law 1 — Dual-Tier Plate Whitespace & Height Targets (Strict on Main, Relaxed & Clean on Supp/ED!)
- **Main Figures (High-Substance Flagship Standard)**:
  - **Multi-modal Main Figure (with real microscopy/MRI/H&E/spatial map/heatmap anchor)**: **`52–72%` white** (flag **`> 76%` white** as sparse!), height **`≥ 72–82%`**.
  - **Statistical / Clinical Main Figure**: **`≤ 75%` white** (achieve substance via genuine multi-dimensional encodings, NOT decorative clutter or gridlines!).
  - **Every Main Figure Needs ≥ 1 High-Density Visual Anchor**: Never let an entire Main Figure become `8–10` bare univariate scatterplots, simple barplots, and 2-line KM curves (`78–82%` pure white). Anchor every Main Figure with at least ONE of:
    1. **Real visual imaging Hero zone (`35–55%` area)**: clinical MRI/CT, H&E, mIF/CODEX/IHC, Visium/Xenium tissue ROI, gross excised tumors, or Western blot; OR
    2. **High-density data visualization anchor**: multi-track clinical OncoPrint, clustered multi-block heatmap, 2D spatial boundary/streamline/bivariate map, grouped marker dotplot with module brackets, multi-cluster volcano, or biological vector illustration.
- **Extended Data / Supplementary Figures (Relaxed, Clean & Honest — ZERO Artificial Padding!)**:
  - Supplementary / Extended Data figures exist for QC, subgroup splits, sensitivity checks, and validation scatter/box plots — **clarity and visual simplicity beat pixel density!**
  - **Relaxed whitespace & height**: **`≤ 86%` pure white** is completely natural and acceptable on Supp/ED plates (only flag **`> 88%` white** or **`< 50%` page height** `< 110 mm`).
  - **NEVER artificially pad Supplementary/ED figures** with decorative 2D KDE contour lines, unnecessary shading, or ornamental frames just to lower white-pixel `%`!

### Law 2 — Zero Dropped / Buried Real-Visual Assets Audit Gate
- Before finalizing any multi-figure set, cross-check **100% of extracted / available real imaging files** (`MRI`, `H&E`, `IHC`, `mIF/CODEX`, `mouse IF`, `spatial tissue ROIs`, `WB`, `gross tumor photos`) against the assembled plates:
  1. **Zero orphaned assets**: Never leave extracted microscopy/IF images unused on disk while an *in vivo* or spatial figure lacks tissue-level proof!
  2. **Promote flagship visual proof out of Extended Data**: If a Main Figure sits `> 76%` white with pure scatter/bar plots while its paired Extended Data holds a striking representative mIF or high-power spatial ROI proving that exact figure's core claim, **promote that visual anchor into the Main Figure** and demote a redundant scatter/bar panel to ED!

### Law 3 — Honest Information Density per Chart (Ban Decorative KDE Contours on Small-`n` Scatters!)
- **Kaplan–Meier Survival (`plot_km_survival`)**: ALWAYS include **95% CI shaded ribbons + aligned `Number at risk` table + Cox `HR (95% CI)` + two-sided log-rank `p`** on Main and flagship clinical KM curves (compact clean styling is fine on secondary subgroup Supp KMs).
- **Scatter / Correlation Plots (Zero Padding Contours!)**:
  - **Default clean scatter (MANDATORY for `n < 300`, especially in Supplementary / Extended Data figures)**: pure-white L-axis + crisp scatter points + clean linear/LOESS fit (`± 95% CI` ribbon optional) + concise `ρ, p` annotation. **NEVER slap decorative 2D KDE contour lines onto small- or medium-`n` (`n < 300`) patient/sample scatter plots merely to "fill whitespace"!**
  - **When 2D KDE contours ARE genuinely justified**: ONLY on massive overplotted datasets (`n ≥ 300–500+` overlapping single cells/spots where raw points merge into an unreadable solid blob) or true 2D spatial-density maps.
  - **Main-Figure scatter consolidation**: In Main Figures, avoid chaining `≥ 3` repetitive small-`n` single-variable scatter plots across a row — instead summarize multiple correlations into a **clean correlation matrix, lollipop, or forest plot** (never by padding tiny scatter plots with fake contours!).
- **Univariate Bar Charts**: In Main Figures, avoid stretching a 3–4 bar chart across `≥ 58 mm` of mostly white space; compress to a crisp `38–44 mm` compact dot-CI / lollipop plot or group across cohorts/strata.
- **Forest Plots (`plot_cox_forest`)**: Target `≥ 6–8` informative rows or paired unadjusted/adjusted dual-markers with aligned numeric `HR/OR (95% CI)` + `p` table columns; avoid sparse 3–4 row forest plots stretched across half a page.

### Law 4 — Pragmatic Extended Data / Supplementary Consolidation (Keep Supp Clean & Spacious!)
- Supplementary / Extended Data figures do **NOT** need to be packed as tightly as Main Figures — a clean, airy `55–85%`-height plate with `4–7` well-spaced panels is completely fine!
- Only consolidate when **`≥ 2` ultra-sparse (`≤ 3–4` panel, `< 50%` height `< 110 mm`)** plates split the exact same validation theme (e.g., two tiny 3-panel external-cohort plates both supporting `Fig. 6`) — never force-stretch panels or add ornamental contours just to inflate a Supp plate's pixel density.

### Law 5 — Zero Dead Internal Voids & Ban Plain Code Text-Box Schematics
- **Eliminate Internal White Voids (`≥ 10 mm` empty vertical bands)**: Never leave `15–25 mm` of empty white gap between a top legend/title and a microscopy strip, nor shrink an atlas UMAP to a tiny floating island surrounded by `35 mm` blank margins!
- **Ban Matplotlib/Python Text-Box Flowcharts for Study Designs (`Fig. 1a`) & Final Working Models (`Last Panel`)**: Plain rounded rectangles filled with bullet text + straight arrows look like software architecture slides, not *Nature / Cell* biology! Always craft Study Design schematics and Final Mechanistic Working Models as **biological cell/tissue/organ/receptor illustrations** (BioRender, Adobe Illustrator, or crisp biological SVG vector art — using `add_visual_placeholder_frame()` if assembling vectors later in Illustrator).

### Law 6 — In-Figure Anti-AI, UK/US Spelling Purity & 1:1 Text↔Figure Terminology Hygiene
- Every panel title, axis label, and annotation inside a figure MUST pass the **exact same anti-AI, spelling, and anti-notebook-jargon filter** as the manuscript text (`sci-paper-writing` `14-anti-ai-prose-stat-restraint-and-rule-harmony.md`):
  - **NEVER print internal pipeline/cohort/cluster codes** inside publication figures (`c02_Fib_POSTN`, `m01_Macro_SPP1`, `PipeV2_DropAll`, `MalignantGate`, `Primary40`, `ValSet15`, `OOF`, `Quarantined Val`).
  - **NEVER print non-CNS AI/tech-blog terms** inside figure titles (`audit`, `crosswalk`, `parsimonious`, `habitat`, `out-of-fold`, `phase cliff`, `Adversarial`, `triage cascade`).
  - **Match journal UK vs US spelling INSIDE figure axes & legends**: in UK-spelling *Nature / NCB / Nat Commun* figures, write `Standardised`, `Tumour`, `Normalised`, `Haematoxylin`, `Signalling` — never let a stray `Standardized Z` slip onto an axis! Proofread subtype legends for upstream typos (`Neutro` not `Nutri`, consistent `Mast` casing).

### Law 7 — Multi-Figure Experimental-Design, Assay-Timeline & Spatial-Workflow Schematics (`52–59%` of CNS Papers!)
Empirically across **`360` published CNS manuscripts (`2,990` Main captions)**, **`52–59%` include schematics in `≥ 2` Main Figures** (`25–41%` in `≥ 3`; mean `1.8–2.3` schematics/paper). Do NOT stop at `Fig. 1a` + the final `Working Model` — give readers instant visual orientation in every complex experimental or spatial figure using **compact `22–38 mm` biological vector strips**:
- **`Fig. 1` (`58–62%`)**: Overall clinical cohort & multimodal workflow (`88 × 46 mm`).
- **Spatial / Multi-Modal Figures (`21–37%`)**: Compact **`spatial_workflow_schematic`** (`76 × 38 mm` — e.g., vessel distance rings `≤50 µm` vs `100–200 µm`, co-registered H&E→mIHC workflow).
- ***In Vitro* Cellular Assay Figures (`20–36%`)**: Compact **`cell_assay_schematic`** (`64 × 34 mm` — e.g., EMT induction, co-culture, CRISPR KO + drug IC50 setup).
- ***In Vivo* Mouse / Perturbation Figures (`26–45%`)**: Compact **`mouse_timeline_schematic`** (`88 × 28 mm` or `176 × 24 mm` horizontal strip — e.g., `Day 0 orthotopic implant → Day 14 3T MRI + start ADC/inhibitor → Day 17 MRI → Day 21 MRI & RNA-seq harvest`).
- **Final Working Model (`24–44%`)**: Biological cell/vessel/receptor synthesis diagram (`176 × 42–48 mm`).

### Law 8 — Zero Phantom Scale Bars & 100% Schematic ↔ Data-Panel Numerical Cross-Check
- **Scale-Bar Truthfulness Gate**: NEVER write `(scale bar, 2 cm)` or `(scale bar, 5 mm)` in a Figure Legend unless a calibrated physical scale bar is **physically drawn and visibly verified on the rendered image**! For clinical/preclinical MRI or microscopy panels without an embedded bar: stamp a verified vector bar via `add_spatial_scalebar()` using exact DICOM/OME-TIFF pixel spacing, OR truthfully report exact FOV dimensions (`FOV, 18 × 18 cm`) in the legend.
- **100% Schematic ↔ Data-Panel Cross-Check**: Before final delivery, audit EVERY numerical distance (`µm`), metric (`AUC`, `C-index`, `HR`, `%`, `P`, `n`), gene symbol, and arrow/icon inside `Fig. 1a`, intermediate timelines, and the final `Working Model` against the exact data panels — guaranteeing **zero numerical mismatch** (e.g., never misattributing one gene's distance to another) and **zero icon-on-text overlap**!
- Run `python scripts/sci_plot_py.py --audit-plate <assembled_dir_or_pdf> --journal nature_sister` before final delivery to automatically check **pixel whitespace %**, **page height utilization %**, **dead vertical white bands**, **in-figure AI/jargon/spelling/typo leaks**, and **orphaned visual assets**!
