# Hybrid Vector + Raster Export Cookbook for Adobe Illustrator

How to keep dense scRNA-seq / spatial / volcano / heatmap data layers light and fast in Adobe Illustrator while keeping **every coordinate axis, tick, spine, legend, leader line, and Arial character 100% crisp and editable**.

## 1. Rasterization Decision Thresholds

| Plot Component | Threshold | Export Mode | Resolution |
| :--- | :---: | :--- | :---: |
| UMAP / t-SNE / Spatial cell cloud / Volcano / High-n Scatter | `n_points > 2,000` | **Rasterize point layer only** | **600 dpi** (300 dpi min for huge panels) |
| Single-cell Violin + Jitter / Beeswarm | `n_points > 2,000` | **Rasterize jitter points**; keep violin/box outlines vector | **600 dpi** |
| Expression / Correlation Heatmap | `n_rows × n_cols > 2,500` | **Rasterize matrix body**; keep dendrograms, labels, tracks, colorbar vector | **600 dpi** |
| Bar, Box (`n ≤ 2,000`), KM, Forest, DotPlot, Line, ROC | `≤ 2,000` elements | **Pure 100% Vector** | Vector |

---

## 2. Python (`matplotlib` / `seaborn` / `scanpy`) Recipes

### Global rcParams (Must Set Before Creating Figures)
```python
import matplotlib as mpl
mpl.rcParams.update({
    "pdf.fonttype": 42,              # TrueType live text in Illustrator (NEVER default Type 3!)
    "ps.fonttype": 42,
    "svg.fonttype": "none",          # Live <text> in SVG, never converted to <path>
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom",    # Keep italic gene names & P-values in Arial, not serif!
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "axes.unicode_minus": False,
})
```

### Selective Rasterization Patterns
- **UMAP / Scatter / Volcano**:
  ```python
  ax.scatter(x, y, c=colors, s=1.2, alpha=0.85, linewidths=0, rasterized=True)
  # Highlighted top genes / threshold lines / text stay rasterized=False (vector!)
  ax.axhline(y_cut, ls="--", lw=0.4, color="#666666")
  ax.text(gx, gy, "SPP1", fontsize=5.5, style="italic")
  fig.savefig("volcano.pdf", dpi=600)  # dpi ONLY affects rasterized=True artists!
  ```
- **Scanpy UMAP (`sc.pl.umap`)**:
  ```python
  sc.pl.umap(adata, color="cell_type", ax=ax, show=False, frameon=False)
  for coll in ax.collections:
      coll.set_rasterized(True)
  fig.savefig("umap.pdf", dpi=600)
  ```
- **Large Heatmap (`imshow` / `pcolormesh` / `seaborn.heatmap`)**:
  ```python
  im = ax.imshow(mat, aspect="auto", cmap=cmap, interpolation="nearest", rasterized=True)
  # OR seaborn:
  sns.heatmap(df, ax=ax, cmap=cmap, rasterized=True)
  # Row/col gene labels, ticks, colorbar frame & ticks stay 100% vector!
  fig.savefig("heatmap.pdf", dpi=600)
  ```
- **Z-order trick for mixed raster + vector axes**:
  Set `ax.set_rasterization_zorder(1)` and give dense scatter/mesh `zorder=0`, while lines/annotations/boxes get `zorder>=2`. Everything with `zorder < 1` becomes one single crisp 600-dpi bitmap layer underneath editable vector axes!

---

## 3. R (`ggplot2` / `ggrastr` / `ComplexHeatmap` / `Seurat`) Recipes

### Critical Font-Editability Trap in R
- **NEVER call `showtext::showtext_auto()` when saving PDF/SVG for Adobe Illustrator** — `showtext` converts every character into uneditable polygon outlines!
- **Preferred AI-Editable Devices**:
  1. `grDevices::cairo_pdf("fig.pdf", width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")` (Windows/Linux/macOS; embeds live TrueType Arial!).
  2. `svglite::svglite("fig.svg", width = w_mm / 25.4, height = h_mm / 25.4, system_fonts = list(sans = "Arial"))`.
  3. Fallback `grDevices::pdf("fig.pdf", width = w_mm / 25.4, height = h_mm / 25.4, family = "Helvetica", useDingbats = FALSE)` (always set `useDingbats = FALSE` so scatter dots never turn into garbled `AdobePiStd` glyphs in AI!).

### Selective Rasterization Patterns in R
- **`ggplot2` + `ggrastr`** (UMAP, Volcano, large Scatter, Cell-level Violin Jitter):
  ```r
  library(ggplot2)
  library(ggrastr)

  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = cell_type)) +
    rasterise(geom_point(size = 0.25, stroke = 0, alpha = 0.85), dpi = 600) +
    theme_cns()
  ggsave("umap.pdf", p, device = cairo_pdf, width = 52, height = 50, units = "mm")
  ```
- **`Seurat::DimPlot` / `FeaturePlot` (built-in raster support)**:
  ```r
  p <- DimPlot(obj, group.by = "cell_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 1.2) +
    theme_cns()
  ```
- **`ComplexHeatmap` (large gene × cell heatmaps)**:
  ```r
  library(ComplexHeatmap)
  ht <- Heatmap(
    mat,
    use_raster = TRUE,
    raster_by_magick = TRUE,   # Crisp high-quality rasterization (requires magick)
    raster_quality = 5,        # Or raster_device = "CairoPNG", raster_device_param = list(res = 600)
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 5.5, fontface = "italic"),
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 5.5),
    heatmap_legend_param = list(
      title_gp = gpar(fontfamily = "Arial", fontsize = 6.5, fontface = "bold"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
    )
  )
  cairo_pdf("heatmap.pdf", width = 120 / 25.4, height = 65 / 25.4, family = "Arial")
  draw(ht)
  dev.off()
  ```

---

## 4. Eliminating Broken / Multi-Element Text in Adobe Illustrator

Why do PDFs from matplotlib or R often open in Adobe Illustrator with one label split into 2–15 separate elements (e.g., `V` + `alidation`, `VEGF` + `A`, or `p`, `=`, `0`, `.`, `0`, `0`, `4` as 9 isolated objects)?

### Two Root Causes in PDF Streams
1. **Mathtext per-character `Td` shattering (Python matplotlib)**: Any label containing `$...$` (such as `$p$ = 0.004`, `$\mathit{\rho}$ = 0.68`, `log$_2$(Fold change)`) routes through matplotlib's `draw_mathtext()`, which emits a separate `dx dy Td (c) Tj` cursor jump for **every single character** — splitting even a 9-character `p = 0.004` into 9 separate text objects in AI.
2. **Optical kerning `TJ` splits (matplotlib & R `cairo_pdf`)**: For plain text containing kerned pairs (`Va`, `To`, `FA`, `Tu`, `We`), PDF writers emit `[ (V) 73.6 (alidation) ] TJ` or per-glyph positions. Illustrator splits these at kerning offsets into disconnected Point-Type objects.

### Built-In Three-Layer Fix
1. **PDF Stream Coalescer (`coalesce_pdf_text_runs()` in `sci_plot_py.py`)**:
   - Collapses optical kerning `TJ` splits (`[(V) 73.6 (alidation)] TJ` → `(Validation) Tj`) and merges per-character mathtext `Td` runs into one `BT...ET` (fixes 100% of plain text and 2-run labels like `*p* = 0.004`, `*p* < 0.0001`).
2. **Semantic SVG Coalescer (`coalesce_svg_text_runs()` in `sci_plot_py.py`, auto-run when `save_svg=True`)**:
   - Rewrites matplotlib's scrambled font-grouped SVG `<tspan x="...">` into ONE reading-order `<text>` element with inline `<tspan font-style="italic">` and `<tspan baseline-shift="...">` (ZERO `x`/`y` attributes on inner `<tspan>`s!) while keeping 600-dpi embedded `<image>` bitmaps for UMAP/Volcano/Heatmap.
3. **One-Click Adobe Illustrator Script (`scripts/MergeSciPlotText.jsx`) — 100% Guarantee for ANY PDF**:
   - Adobe Illustrator's own PDF-to-ATE importer inherently breaks PDF `Ts` subscripts (`log₂(Fold change)`, `−log₁₀(p_adj)`) and 3+ inline font switches (`Validation (n = 18)`, `ρ = 0.68, p = 0.002`) during PDF import — even when the PDF stream is 100% standard ISO 32000.
   - Running `File -> Scripts -> MergeSciPlotText.jsx` inside Illustrator merges ALL remaining horizontal and 90°-rotated fragments across the entire figure into **true single-line native Illustrator `PointText` objects** with preserved `Arial` / `Arial Italic` and native `baselineShift` subscripts/superscripts in < 0.3 s!
   - Works equally on PDFs from Python (`matplotlib`/`scanpy`), R (`ggplot2`/`cairo_pdf`/`ComplexHeatmap`), Prism, or BioRender.
- **R (`sci_plot_r.R`)**:
  - **Cleanest 1-line vector format**: `svglite::svglite(..., fix_text_size = FALSE, system_fonts = list(sans = "Arial"))` (NEVER use `showtext` or `ggtext` for plain labels).
  - **For PDF (`cairo_pdf`)**: Pass `fix_pdf_py = ".../sci_plot_py.py"` to `save_ai_pdf()`, or run `MergeSciPlotText.jsx` once opened in Illustrator.

---

## 5. Statistical Symbol, Italic Greek & Spacing Rules

| Category | Correct Format | Avoid | Helper (`sci_plot_py.py` / `sci_plot_r.R`) |
| :--- | :--- | :--- | :--- |
| *p*-value / *q*-value | *p* = 0.004; *p* < 0.0001; *q* = 0.012 | `P=0.004`, `p<0.001` (no space / uppercase) | `fmt_p(0.004)` / `fmt_p_expr(0.004)` |
| Scientific *p*-value | *p* = 1.4 × 10⁻⁹ | `p=1.4e-9`, `p = 1.4×10^-9` | `fmt_p(1.4e-9, sci=True)` |
| Correlation & sample size | *ρ* = 0.68, *p* = 0.003 (*n* = 18); *r* = 0.72 | `ρ=0.68,P=0.003 ( n=18 )` | `fmt_corr(0.68, p=0.003, n=18)` |
| Greek letters | Lowercase italic: *α*, *β*, *γ*, *δ*, *ρ*, *μ*, *σ*, *χ*² | Upright or mismatched serif Greek | `fmt_greek("gamma")` |
| Axes with log / *p* | log₂(Fold change); −log₁₀(*p*_adj) | `-log10(Padj)` (hyphen instead of `−`) | `fmt_log2fc()` / `fmt_neglog10p()` |
| Spacing rules | Single space around `=`, `<`, `>`, `×`, `±`; space before unit (`50 μm`, `24 h`); NO space before `%` (`75%`) or inside `()` | `50μm`, `75 %`, `( n = 18 )`, `2.3±0.4` | Enforced in all helpers |

---

## 6. Illustrator Assembly Checklist

When opening exported `.pdf` files in Adobe Illustrator:
1. **Place (`File -> Place`) or Open at 100% scale**: Artboard matches exact Nature (`180 mm`) or Cell (`174 mm`) width — zero manual rescaling needed.
2. **Release Clipping Masks sparingly**: Keep the rasterized point cloud grouped under its clipping frame while editing text, moving `ggrepel` labels, or tweaking legend alignment.
3. **Verify Live Fonts**: `Type -> Find/Replace Font` should list only `Arial` (Regular / Bold / Italic) with zero missing fonts and zero outlined text blocks.
4. **Final Journal Export**: Keep master editable `.ai`/`.pdf` with linked/embedded 600-dpi layers; submit PDF/EPS/TIFF per journal final production portal.
