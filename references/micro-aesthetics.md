# CNS Micro-Aesthetics Specification (Alpha, Point/Bar/Box Sizes, Strokes, Ticks & Sorting)

The subtle numerical parameters that separate unmistakable *Nature* / *Cell* / *Science* figures from default `matplotlib` / `ggplot2` plots.

## 1. Stroke Width Hierarchy (Visual Figure-Ground Separation)

| Visual Layer | Exact Width (`pt`) | ggplot2 `linewidth` (`pt / 2.13`) | Color & Style |
| :--- | :---: | :---: | :--- |
| **Paired sample connecting lines** | **0.30–0.35 pt** | `0.14–0.16` | `#8C8C8C` or group color, `alpha = 0.25–0.35` |
| **Reference / cutoff dashed lines** (`y=0`, `HR=1`, volcano) | **0.35 pt** | `0.16` | `#777777`, fine dash `(2.5, 1.5)` |
| **Significance brackets & leader lines** (`ggrepel`) | **0.35–0.40 pt** | `0.16–0.19` | `#1A1A1A` (brackets) / `#666666` (leader lines) |
| **Axis spines, ticks, colorbar border** | **0.40 pt** | `0.19` | `#1A1A1A`, outward ticks `2.0 pt` |
| **Error bars (SEM / SD / 95% CI)** | **0.45–0.50 pt** | `0.21–0.23` | `#1A1A1A` or group color; cap width `0.18–0.22` |
| **Boxplot outline & whiskers** | **0.55–0.60 pt** | `0.26–0.28` | Group color (`fill = "white"` hollow box!) |
| **Boxplot median bar** | **0.90–1.00 pt** | `0.42–0.47` | Group color (1.6× thicker than box outline!) |
| **Primary data curves (KM survival, ROC, growth, trend)** | **0.95–1.20 pt** | `0.45–0.56` | Group color (2.5× thicker than axes!) |

---

## 2. Bar Width, Bar Spacing (柱间距), Box & Violin Geometry

| Parameter | CNS Standard Value | Resulting Gap / Spacing Rule | Anti-Pattern to Avoid |
| :--- | :--- | :--- | :--- |
| **Single-series vertical bar** | **`width = 0.58–0.62`** | Inter-bar gap = **`0.38–0.42`** (~65% of bar width) | Default `0.8–0.9` looks like thick brick wall |
| **Grouped / dodged bars (2–4 arms)** | **`width = 0.66`, `dodge = 0.75`**<br>*(Py: cluster `0.72`, inner slit `0.04`)* | **Intra-group slit `~0.04` << Inter-group gap `~0.25–0.28`** (Gestalt proximity!) | Equal spacing between & within groups (hard to group visually) |
| **100% stacked bar (≤ 15 bars)** | **`width = 0.75`** | Tight gap **`0.25`** (~1/3 bar width) aids cross-bar trend tracking | Wide `0.50` gaps breaking composition continuity |
| **100% stacked cohort bar (> 20 patients)** | **`width = 0.88–0.95`** (or `1.0`) | Hairline slit `0.05–0.12` (or seamless ribbon) | Massive white stripes dominating patient strips |
| **Horizontal bar (pathways/enrichment)** | **`width = 0.65–0.70`** | Row gap **`0.30–0.35`** (~50% of bar height) | Cramped touching bars or huge empty gaps |
| **Outer left/right x-padding** | **`add = 0.55–0.65`** (`xlim(-0.6, k-0.4)`) | Balanced breathing room beside y-axis and right edge | Bars glued to y-axis or huge wasted side margins |
| **Boxplot width & spacing** | **`width = 0.48–0.56`** (`dodge = 0.74`) | Clean gap `~0.46`; **ALWAYS `outlier.shape = NA` / `showfliers = False`** when jitter overlaid | Double-plotting outliers on top of jitter |
| **Bar / Box fill style** | **Box**: hollow white fill + colored `0.55 pt` outline (or soft `alpha = 0.25`)<br>**Bar**: soft solid fill `alpha = 0.85` (borderless or `0.35 pt` outline) | Clean figure-ground contrast | Harsh 100% neon fill + thick `1.5 pt` black border |
| **Bar Y-axis baseline** | **Zero bottom gap** (`expand = c(0, 0.16)` / `ylim(bottom=0)`) | Bars sit flush on x-axis | Floating bars above `y = 0` (ggplot2 default!) |
| **Top headroom for brackets/labels** | **`+14%` to `+22%`** above max data value | Zero bracket–title collision | Brackets colliding with panel title |

---

## 3. Point Size, Outline & Alpha Matrix

| Plot Type & Point Count | Python `s` (pt²) | R `size` (mm) | `alpha` | Edge / Outline | Jitter / Layout |
| :--- | :---: | :---: | :---: | :--- | :--- |
| **Bar / Box biological replicates** (`n = 3–30`) | `7–12` | `0.9–1.3` | `0.65–0.85` | `lw = 0` (colored) or `0.3 pt` white/dark rim (`shape = 21`) | Tight jitter **`±0.12–0.16`** or `quasirandom(width = 0.20)` |
| **Cohort scatter / correlation** (`n = 30–500`) | `9–16` | `1.1–1.5` | `0.72–0.85` | `0.25 pt` white rim separates overlapping dots | Exact data coords |
| **Volcano: NS background** (`>2k` genes) | `0.6–1.0` | `0.20–0.35` | `0.45–0.55` | `lw = 0` (`rasterized=True`) | `zorder = 1`, `#D3D3D3` |
| **Volcano: Significant Up/Down** | `1.4–2.4` | `0.40–0.65` | `0.75–0.85` | `lw = 0` (`rasterized=True`) | `zorder = 2`, `#C84B4B` / `#3B6EA5` |
| **Volcano: Labeled focal genes** | `4.5–7.0` | `0.90–1.20` | `1.00` | `0.35 pt` `#1A1A1A` rim (vector!) | `zorder = 4` + thin leader line |
| **UMAP (`≥150k` cells)** | `0.10–0.18` | `0.05–0.08` | `0.75–0.88` | `lw = 0` (`rasterized=True`, 600 dpi) | Shuffled row order! |
| **UMAP (`50k–150k` cells)** | `0.30–0.55` | `0.12–0.20` | `0.80–0.90` | `lw = 0` (`rasterized=True`, 600 dpi) | Shuffled row order! |
| **UMAP (`15k–50k` cells)** | `0.75–1.40` | `0.25–0.40` | `0.85–0.92` | `lw = 0` (`rasterized=True`, 600 dpi) | Shuffled row order! |
| **UMAP (`<15k` cells)** | `1.50–3.00` | `0.45–0.75` | `0.85–0.95` | `lw = 0` | Shuffled row order! |
| **Marker DotPlot** | `max_diam ≤ 80%` cell pitch | `range = c(0.4, 3.2)` | `1.00` | **`0.25–0.35 pt` `#333333` rim (`shape = 21`)** so pale dots stay visible! | Matrix grid |
| **CI / SEM shaded ribbons** | — | — | **`0.15–0.22`** | Borderless (`edgecolor="none"`) | Behind main curve (`zorder = 1`) |

---

## 4. Ticks, Scaling, Heatmaps, Brackets & Sorting Rules

1. **Tick Density & Numbers**:
   - Keep **3–5 major ticks per axis** (e.g., `0, 25, 50, 75, 100` instead of `0, 10, 20...100`).
   - Use **true typographic minus (`−` U+2212)** on negative ticks (`axes.unicode_minus = True` with Arial!).
   - Uniform decimal places on an axis (`0.0, 0.5, 1.0`); format large integers with commas (`10,000`) or `10⁴`, never `1e+05`.
2. **Heatmap Z-score Clipping & Colorbars**:
   - Always clamp Z-score heatmaps symmetrically at **`[-2, +2]`** (or `[-2.5, +2.5]` / 2nd–98th percentile) so single outliers never wash out the palette.
   - Keep colorbars slim (`2.5–3.2 mm` wide × `16–24 mm` tall), `0.35 pt` frame, 3–5 ticks (`−2, 0, 2`).
3. **Significance Bracket Hygiene**:
   - Max **2–3 stacked brackets** per panel (tip drop = `2–3%` of y-range). If >4 pairs exist, report one omnibus test (`Kruskal-Wallis, p = ...`) or compact letters (`a, b, c`) instead of a tall bracket tower.
   - Prefer **exact lowercase italic *p*-values** (`p = 0.004`, `p < 0.0001`) over bare stars in 2–3 group comparisons.
4. **Monotonic Ordering**:
   - Never leave unordered bars, pathways, cell types, or cancer types in alphabetical order — **sort by value/median/NES** (keep biological order only for ordered stages: `Normal → Primary → Metastasis`, `CR → PR → SD → PD`, timepoints).
5. **Biological & Statistical Italicization**:
   - **Human genes / mRNA**: ALL-CAPS Italic (*TP53*, *SPP1*, *CXCL13*, *PDCD1*)
   - **Mouse genes / mRNA**: Title-case Italic (*Spp1*, *Trp53*, *Krt19*)
   - **Proteins / markers**: Upright (PD-1, PD-L1, CD8, Ki67, EpCAM, IFN-*γ*, TNF-*α*)
   - **Latin phrases & species**: Italic (*in vivo*, *in vitro*, *ex vivo*, *de novo*)
   - **Stats & Greek**: Lowercase italic (*p*, *q*, *r*, *n*, *t*, *z*, *ρ*, *α*, *β*, *γ*, *μ*, *χ*²; *R*² uppercase italic *R*).
6. **Multi-Panel Row Alignment**:
   - Panels sharing one row should align cleanly along exact `top` and `bottom` axis baselines so x-axes form a continuous horizontal line across the 180 mm page.
