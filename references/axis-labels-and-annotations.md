# Concise Axis Labels & Annotation Copywriting Guide

Empirically distilled from **1,016 rotated Y-axis labels** and **511 horizontal X-axis labels** mined verbatim from 60 Tier-1 CNS research articles (*Nature*, *Nature Medicine*, *Nature Genetics*, *Nature Immunology*, *Nature Cancer*, *Nature Methods*, *Cell*, *Cancer Cell*, *Immunity*).

## 1. Empirical CNS Benchmarks (Hard Limits)

| Metric | Real CNS Measured Distribution | Hard Skill Guideline |
| :--- | :--- | :--- |
| **Word count** | Median **2 words** (IQR `2–3`, P90 = **4 words**) | Target **1–3 words**; max **4 words** (wrap with `\n` if unavoidable) |
| **Character length** | Median **15 chars** (IQR `10–20`, P90 = **25 chars**) | Target **≤ 22 chars**; hard ceiling **≤ 28 chars** on one line |
| **Capitalization** | **Sentence case 92.6%** (`327` vs `26` Title Case) | **Always Sentence case** (only capitalize 1st word, acronyms, genes) |
| **Units format** | Parentheses **21.2%** vs brackets 3.2% (`6.6:1`) | Single space + parentheses: `Time (months)`, `Volume (mm³)`, `Proportion (%)` |
| **Code column names (`_`)** | Only 3.9% contain `_` | **Never** leak raw `snake_case` columns (`cell_prop`, `log2_fc`) onto axes |

---

## 2. Five Golden Rules to Eliminate Redundant Axis Labels

1. **Delete Self-Explanatory X-Axis Titles (`xlab = NULL` / `set_xlabel("")`)**:
   - When x-axis ticks are already clear categories (`Normal`, `Tumor`; `CR`, `PR`, `SD`; `MSS`, `MSI`; `Vehicle`, `AG120`; `Myeloid`, `CD8 T`, `Malignant`; italic gene symbols), **leave the x-axis title completely blank**.
   - Never print redundant placeholders like `Group`, `Condition`, `Sample type`, `Cell type`, `Cluster`, or `Gene`.
2. **Keep Methods & Normalization Details Out of Axis Labels**:
   - Write `Expression`, `Signature score`, or `Z-score` on the axis — move `log2(TPM + 1) batch-corrected normalized via Seurat AddModuleScore` to the Figure Legend / Methods!
3. **Never Duplicate the Panel Title on the Y-Axis**:
   - If the panel title already says `CD8+ T cells` or `SPP1`, the y-axis should only say `Frequency (% CD45⁺)` or `Expression` — NOT `SPP1 expression in CD8+ T cells`.
   - Alternatively, omit the panel title completely and put concise `SPP1 expression` or `CD8⁺ T cell (%)` on the y-axis!
4. **Strip Empty Filler Words**:
   - Drop `level of`, `value of`, `values`, `total number of`, `distribution of`, `comparison of`, `between X and Y`.
   - Prefer `Cell count` over `Total number of cells`; `Expression` over `Gene expression level`; `log₂(fold change)` over `Log2 fold change (Tumor vs Normal)`.
5. **Avoid Double Units & Redundant Scales**:
   - If Y ticks are `0, 25, 50, 75, 100` → label `Proportion (%)` or `Overall survival (%)` (NO `%` on every tick).
   - If Y ticks are `0.0, 0.25, 0.50, 0.75, 1.00` → label `Proportion`, `Cell fraction`, or `Survival probability` (NO `(%)`!).

---

## 3. Before vs After Rewrite Examples (Anti-Redundancy Clinic)

| Redundant / Poor Label (Avoid ❌) | Crisp CNS Label (Use ✅) | Why Better |
| :--- | :--- | :--- |
| `Log2 Fold Change of Gene Expression (Tumor vs. Normal)` | `log₂(fold change)` or `log₂FC` | Cuts 54 chars → 17 chars; comparison belongs in title/caption |
| `Negative Log10 Adjusted P-Value (-log10 Padj)` | `−log₁₀(p_adj)` or `−log₁₀(FDR)` | Removes duplicate verbal + math text; lowercase italic *p* |
| `Normalized Log1p Expression Level of SPP1 Gene` | Title: *SPP1* / Y: `Expression` (or Y: `*SPP1* expression`) | Drops transform prose & filler word `Level`/`Gene` |
| `Percentage of CD8+ T Cells Among All Immune Cells (%)` | `CD8⁺ T cells (% immune)` or `Frequency (% CD45⁺)` | Concise biological denominator in parentheses |
| `Probability of Overall Survival (%)` (with 0–1 ticks) | `Overall survival` or `Survival probability` | Fixes 0–1 vs `%` unit mismatch; drops filler |
| `Follow-up Time After Treatment Initiation (Months)` | `Time (months)` or `Months` | Matches top CNS Kaplan-Meier convention |
| `Uniform Manifold Approximation and Projection 1` | `UMAP1` (mini L-arrow) | Standard 5-char compact label |
| `cell_type_cluster` / `treatment_group` | *(omit x-axis title completely!)* | Self-explanatory ticks never need generic x-titles |
| `Gene Set Enrichment Normalized Score (NES)` | `NES` or `Normalized enrichment score` | Drops redundant acronym expansion |
| `Spearman Correlation Coefficient Value` | `Spearman ρ` or `Correlation (ρ)` | Direct lowercase italic *ρ* |
| `Subcutaneous Tumor Volume Measurement (mm3)` | `Tumor volume (mm³)` | Proper superscript `³`, drops filler `Measurement` |

---

## 4. Verbatim CNS Axis-Label Dictionary by Chart Type

| Chart Type | Standard Concise X-Axis | Standard Concise Y-Axis | Notes |
| :--- | :--- | :--- | :--- |
| **UMAP / t-SNE / PCA** | `UMAP1` / `t-SNE1` / `PC1 (24%)` | `UMAP2` / `t-SNE2` / `PC2 (12%)` | Use mini equal-length L-arrows on UMAP/t-SNE |
| **Volcano Plot** | `log₂(fold change)` or `log₂FC` | `−log₁₀(p_adj)` or `−log₁₀(p)` | Add tiny arrow cues below x-axis (`Normal ← → Tumor`) if needed |
| **Kaplan-Meier Survival** | `Time (months)` / `Time (years)` / `Days` | `Overall survival` / `Progression-free survival` / `Survival probability` | Use `(%)` ONLY if y-ticks are `0–100` |
| **Boxplot / Violin / Bar (Groups)** | *(Omit! Leave blank)* | `Signature score` / `Expression` / `Proportion (%)` / `Cell density (mm⁻²)` | Put group comparison `p` directly above brackets |
| **100% Stacked Bar** | *(Omit if samples/stages obvious)* | `Proportion` (0–1) or `Proportion (%)` (0–100) | Horizontal bar: X = `Proportion`, Y = *(omit)* |
| **DotPlot / Heatmap** | *(Omit — italic gene ticks speak for themselves)* | *(Omit — cluster ticks speak for themselves)* | Colorbar: `Z-score` or `Scaled expression`; Size: `Cells (%)` |
| **Pathway Enrichment (Bar/Dot)** | `NES` or `−log₁₀(p_adj)` or `Gene ratio` | *(Omit — pathway names are self-explanatory)* | Trim verbose prefixes (`HALLMARK_`, `GOBP_`) and convert to Sentence case |
| **GSEA Running Curve** | `Rank in gene list` | `Running enrichment score` | Compact 2-line wrap if panel width < 45 mm |
| **Forest Plot (Cox / Logistic)** | `Hazard ratio (95% CI)` or `Odds ratio (95% CI)` | *(Omit — variable table on left)* | Log-scale x-axis with dashed vertical line at `1.0` |
| **Tumor Growth / Time Course** | `Days post implantation` / `Time (days)` | `Tumor volume (mm³)` / `Tumor weight (mg)` | Use `Tumour` for UK journals (*Nature*), `Tumor` for US (*Cell*/*Science*) |
| **Correlation Scatter** | Short variable A (e.g. `CD8⁺ T score`) | Short variable B (e.g. `SPP1 expression`) | Tuck `ρ = 0.68, p = 0.002` inside empty plot corner |
| **ROC / PR Curve** | `1 − Specificity` (or `False positive rate`) | `Sensitivity` (or `True positive rate`) | Strict 1:1 square; tuck `AUC = 0.89` inside bottom-right |
| **Pseudotime** | `Pseudotime` | `Expression` (or gene name `*TCF7*`) | Smooth curve + ribbon |
