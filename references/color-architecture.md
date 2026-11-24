# Color Architecture: Single-Figure Harmony & Whole-Manuscript Consistency

Empirically calibrated from 134 main figures across 32 CNS/Tier-1 papers + 149 discrete and 39 gradient palettes from `scidraw.msttnote.cn` + Bang Wong (*Nature Methods* Points of View) colorblind-safe standards.

## 1. Quantitative CNS Color Laws

- **Chromatic Ink Ratio**: Only **5–10% (median 7.5%)** of a CNS figure page is colored; **~90% is crisp white background + black/charcoal text and neutral grays**.
- **Controlled Saturation & Lightness**: Real CNS palettes sit at **HSL Saturation `0.44–0.63`** and **Lightness `0.57–0.71`** — noticeably softer and cleaner than raw neon RGB or high-saturation web defaults (`S > 0.85`). Reserve high saturation strictly for focal highlights.
- **Colorblind-Safe Baseline**: Never pair pure red `#FF0000` with pure green `#00FF00` on heatmaps or scatter plots. Use **Blue ↔ Red/Coral**, **Teal ↔ Terracotta**, or **Purple ↔ Gold/Amber**.

---

## 2. Single-Figure 4-Layer Color Separation Rule

Within ANY multi-panel main figure, assign distinct visual domains so colors never collide in meaning:

1. **Layer 1 — Primary Categorical Comparison (2–6 groups)**:
   - Shared verbatim across every panel in the figure that compares those groups (so downstream panels need zero legends).
2. **Layer 2 — Background / Non-Significant Neutral Anchor**:
   - All non-significant genes/clusters, reference atlas backgrounds, and baseline controls stay in **neutral gray**: `#E6E6E6` (light background cells), `#BDBDBD` (non-significant dots), `#595959` (reference/vehicle bar).
3. **Layer 3 — Lineage Hue-Family Hierarchy (for 10–40 single-cell subclusters)**:
   - Assign each major lineage a hue family, then generate 3–8 tonal steps (light→dark / warm→cool) for its subclusters:
     - `CD4 T` → Purple/Lavender family (`#7A5195`, `#9B7BB5`, `#BC9FD4`, `#5E3C75`)
     - `CD8 T` → Amber/Ochre/Orange family (`#D97724`, `#EDA04C`, `#F4C27F`, `#A65214`)
     - `NK` → Leaf/Emerald Green family (`#2E8B57`, `#63B37F`, `#9CD4AD`)
     - `B / Plasma` → Sky/Cyan/Rose family (`#2B9EB3`, `#70C4D4`, `#D86C96`, `#EFA2BE`)
     - `Myeloid` → Teal/Mint family (`#1F8A70`, `#57B59C`, `#94D5C3`, `#145E4C`)
     - `Stromal / CAF` → Warm Brown/Taupe (`#8C6239`, `#B58D67`, `#D6B89A`)
     - `Malignant / Tumor` → Crimson/Brick Red (`#B22222` / `#C84B4B`)
4. **Layer 4 — Orthogonal Continuous & Diverging Scales**:
   - **Diverging** (Z-score, log2FC, correlation `-1..+1`, spatial interaction): Soft **Blue–White–Red** (`#2166AC → #F7F7F7 → #B2182B`) or **Teal–White–Burgundy** (`#2B7A78 → #FFFFFF → #9E2A2B`).
   - **Sequential 1** (expression, gene activity, density): **LightGray → Crimson** (`#EFEFEF → #F4A582 → #B2182B`), **LightGray → Indigo** (`#EFEFEF → #8C96C6 → #2D3184`), or **Viridis / Magma**.
   - **Sequential 2 (when the SAME figure already has a colorful heatmap!)**: Use **Monochrome White → Charcoal (`#F7F7F7 → #252525`)** or **White → Warm Gold (`#FFFFFF → #FDB813`)** so adjacent heatmaps never visually fight (classic *Nature Medicine* trick).

---

## 3. Whole-Manuscript Color Architecture (`palette.json`)

For any multi-figure paper (Fig. 1–6 + Extended Data), create and commit ONE `palette.json` in the project directory BEFORE generating Fig. 1:

```json
{
  "condition": {
    "Normal": "#6B8E9B",
    "Control": "#595959",
    "Primary_Tumor": "#3B6EA5",
    "Metastasis": "#C84B4B"
  },
  "response": {
    "Responder": "#2A9D8F",
    "Non_Responder": "#D1495B"
  },
  "lineage": {
    "Malignant": "#C84B4B",
    "CD4_T": "#7A5195",
    "CD8_T": "#E08E45",
    "NK": "#4C9F70",
    "B_cell": "#3D85C6",
    "Plasma": "#D6729B",
    "Myeloid": "#2A9D8F",
    "Fibroblast": "#8C6239",
    "Endothelial": "#6C7A89"
  },
  "continuous": {
    "diverging": ["#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"],
    "expression": ["#EFEFEF", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"],
    "score": "viridis"
  }
}
```

**Global Consistency Rules**:
- **Never change** the hex code of a core biological arm (`Tumor` vs `Normal`, `Responder` vs `Non-responder`, major cell lineages) between Fig. 1 and Fig. 6.
- **Unify continuous scales** paper-wide: pick ONE diverging heatmap scale and ONE expression scale for the entire manuscript.
- **Isolate figure-local variables**: if Fig. 3 introduces a brand-new 3-level spatial zone (`TC / IM / S`), pick colors that do NOT duplicate `Responder/Non-responder` (e.g., Plum `#54278F`, Steel `#3B5B92`, Teal `#2CA089`).

---

## 4. Locked Default CNS Color System & Curated Alternatives

### Locked System Default (Auto-Applied by `sci_plot_py.py` & `sci_plot_r.R`)
Empirically #1 across 45/45 Tier-1 CNS papers (**Blue 100%, Red 100%, Amber 93%, Sage/Teal 91%, Plum 89%, Rose 82%, Gold 89%**):

- **Default Categorical (`CNS_DEFAULT`)**:
  - 2 groups: `["#3B6EA5", "#C84B4B"]` (Steel Blue vs Brick Coral-Red)
  - 3 groups: `["#3B6EA5", "#C84B4B", "#2A9D8F"]` (+ CNS Teal)
  - 4 groups: `["#3B6EA5", "#C84B4B", "#2A9D8F", "#E08E45"]` (+ Warm Ochre-Amber)
  - Full 10-color core: `["#3B6EA5", "#C84B4B", "#2A9D8F", "#E08E45", "#7A5195", "#6BAED6", "#D6729B", "#7DB977", "#E2C062", "#8C6239"]`
  - Extended 11–20 pastel companions: `["#9ECAE1", "#E89A9A", "#8CD1C5", "#F2C28B", "#B39BC8", "#A8D8EA", "#EBB2C8", "#B5D9A8", "#EFE09E", "#BFA084"]`
- **Default Diverging (`rdbu_soft`)**: `["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]`
- **Default Sequential Expression (`gray_crimson`)**: `["#EFEFEF", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"]` (secondary conflict-free options: `gray_indigo`, `white_charcoal`, `viridis`)
- **Default Volcano (`volcano`)**: `{"down": "#3B6EA5", "ns": "#D3D3D3", "up": "#C84B4B"}`
- **Default Neutral Anchors**: `#D9D9D9` (background / NS), `#595959` (Vehicle / Control charcoal)

### Curated Specialized Alternatives
- **Binary Contrast**:
  - `cns_binary_nature`: `["#737DBB", "#E28187"]` (Soft Periwinkle vs Muted Coral — *Nature* 2025 style)
  - `cns_binary_clinical`: `["#2A9D8F", "#D1495B"]` (Teal Responder vs Crimson Non-responder)
  - `cns_binary_warm`: `["#8C6239", "#F2B441"]` (Warm Brown vs Amber Gold — *Science* 2026 MSS/MSI style)
- **3-Group Spatial Zones**:
  - `cns_triad_spatial`: `["#54278F", "#3B5B92", "#2CA089"]` (Plum / Steel / Teal)
- **4-Group Ordered Progression (Ctrl → Low → Mid → High / Timecourse)**:
  - `cns_seq_red`: `["#4D4D4D", "#F4C8B4", "#E58C7E", "#B93B3B"]` (*Science* 2024 treatment ramp)
  - `cns_seq_blue`: `["#BDBDBD", "#9ECAE1", "#4292C6", "#08519C"]`
- **Colorblind-Safe Strict 8-Color**:
  - `okabe_ito_soft`: `["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8C6239", "#595959"]`
- **Alternate Diverging Heatmap**:
  - `teal_plum`: `["#1B7837", "#A6DBA0", "#F7F7F7", "#C2A5CF", "#762A83"]`
