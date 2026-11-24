# Statistical Significance Testing Guide for CNS Figures

Choosing the wrong statistical test—or running cell-level tests on patient-level questions (pseudoreplication)—is one of the fastest ways a high-impact submission gets rejected. Follow this decision matrix whenever computing *p*-values, *q*-values, or significance brackets for a figure.

---

## 1. Core Statistical Guardrails (Non-Negotiable)

1. **Always Two-Sided by Default**: *Nature*, *Cell*, and *Science* require **two-sided** tests unless a strictly pre-registered directional hypothesis justifies a one-sided test.
2. **Prevent Single-Cell Pseudoreplication (CRITICAL)**:
   - **Never** run cell-level Wilcoxon / *t*-tests across thousands of single cells to claim **patient-level** or **group-level** differences (e.g., Responder vs Non-responder, Tumor vs Normal across patients) — this inflates *n* from 15 patients to 50,000 non-independent cells and yields spurious `p < 10⁻³⁰⁰`.
   - **Correct solution**: Aggregate each patient/biological replicate into **one summary value** (pseudobulk expression, sample-level cell proportion, per-patient median/mean module score) and test across biological replicates (`n = patients/mice`), or fit a **Linear Mixed-Effects Model (LMM)** with `(1 | patient_id)`.
3. **Never Omit Multiple-Testing Correction**:
   - **3–6 pairwise comparisons**: Holm–Bonferroni or Benjamini–Hochberg (BH) FDR (`p_adj`).
   - **Omics / pathway / multi-cluster screens (> 10 tests)**: **Benjamini–Hochberg (BH) FDR** (*p*_adj or *q*).
4. **Match Paired vs Unpaired Structure**:
   - Matched samples from the same patient/mouse (Pre vs Post, Tumor vs Adjacent Normal, Core vs Margin) **MUST** use a **paired** test (paired Wilcoxon signed-rank / paired *t*-test / LMM) and visually connect pairs with thin lines (`0.35 pt, alpha = 0.30`).
5. **Minimal On-Plot Text, Full Method in Caption**:
   - **On plot**: concise formatted statistic (`p = 0.004`, `p < 0.0001`, `ρ = 0.68, p = 0.002`, `Kruskal–Wallis, p = 1.4 × 10⁻⁹`).
   - **In caption / report**: always document **(a) exact test name**, **(b) two-sided**, **(c) paired/unpaired**, **(d) multiple-testing adjustment**, and **(e) exact biological unit & `n` per group**.

---

## 2. Master Test-Selection Matrix

| Experimental / Data Design | Default CNS Choice (Robust / Omics) | Parametric Alternative (When Justified) | Post-Hoc / Adjustment | Key Caveats |
| :--- | :--- | :--- | :--- | :--- |
| **2 independent groups** (continuous: expression, score, proportion, tumor size) | **Two-sided Mann–Whitney U** (Wilcoxon rank-sum) | **Two-sided Welch's *t*-test** (normal data; safer than Student's *t* as it does not assume equal variance; use Student's *t* only if normality + equal variance verified, e.g., `n = 3–6` controlled lab replicates) | None (1 test) or BH/Holm if repeated across genes/markers | Always show individual replicate points (`n ≤ 30`) |
| **2 paired groups** (Pre vs Post, Tumor vs Normal same patient) | **Two-sided Wilcoxon signed-rank test** | **Two-sided paired *t*-test** (when paired differences $\Delta$ are normal) | BH/Holm if multi-feature | Require complete pairs; plot subtle connecting lines |
| **≥ 3 independent groups** (stages, subtypes, multi-arm drugs) | **Kruskal–Wallis omnibus test** | **Welch one-way ANOVA** (or classical ANOVA if equal variance) | **Dunn's test (BH)** or pairwise Wilcoxon (BH); parametric: **Tukey HSD** / **Games–Howell**; vs single control: **Dunnett's test** | Only run pairwise brackets if omnibus test is significant or pre-specified; max 2–3 brackets on plot |
| **≥ 3 paired / longitudinal timepoints** (same patients/mice over time) | **Friedman test** (complete data) or **Linear Mixed Model (LMM)** (handles missing timepoints) | **Repeated-measures ANOVA** (with Greenhouse–Geisser sphericity correction) or **LMM** | Pairwise Wilcoxon signed-rank (BH) or `emmeans` contrasts (BH/Tukey) | For tumor growth curves, test `group × time` interaction via LMM or two-way RM-ANOVA (not separate daily *t*-tests!) |
| **Categorical counts / proportions** (2×2 or R×C: CR/PR/SD, metastasis +/-) | **Two-sided Fisher's exact test** (mandatory if any expected count < 5) | **Pearson's $\chi^2$ test** (when all expected counts ≥ 5) | Pairwise Fisher + BH for R×C | For ordered categories (Stage I→IV), consider **Cochran–Armitage trend test** |
| **Bivariate correlation** (2 continuous variables) | **Two-sided Spearman rank correlation (*ρ*)** (default for gene expression, signatures, skewed clinical metrics) | **Two-sided Pearson correlation (*r*)** (strictly linear, bivariate normal, no extreme outliers) | BH if correlating a matrix of variables | Report both *ρ* (or *r*) AND *p* (+ *n*) |
| **Survival / Time-to-event** (OS, PFS, DFS, RFS) | **Two-sided Log-rank test** (KM curves) + **Univariable / Multivariable Cox PH regression** (HR & 95% CI) | Parametric/RMST if curves cross (violating proportional hazards!) | BH if scanning many biomarkers | Check PH assumption (Schoenfeld residuals); show `Number at risk` table |
| **Differential expression (Bulk or scRNA Pseudobulk)** | **DESeq2 / edgeR / limma-voom** on raw summed pseudobulk counts per sample | — | **BH FDR (*p*_adj / *q*)** | Adjust for batch/patient covariates in design formula (`~ batch + condition`) |
| **Differential cell-type abundance (scRNA / spatial)** | **Compositional / sample-level tests**: `propeller` (logit/arcsin), `miloR` (neighborhoods), or sample-level Mann–Whitney / LMM | Beta regression / Dirichlet regression | **BH FDR** | Never run naive $\chi^2$ on raw pooled cell counts across patients! |
| **Pathway enrichment (GSEA / ORA)** | **fgsea** (permutation-based NES + BH *p*_adj) or **Hypergeometric / Fisher ORA** (with proper expressed-gene background!) | — | **BH FDR (*p*_adj / *q*)** | Always specify the actual expressed gene universe for ORA (never entire genome blindly) |

---

## 3. Normality & Variance Quick Decision Rule (for 2–4 Group Continuous Plots)

- **`n < 5` per group** (e.g., `n = 3–4` mice/organoids/qPCR replicates): Normality tests have almost zero power. Two-sided **Welch's *t*-test** (or Student's *t*-test for tightly controlled biological replicates) is standard in *Nature*/*Cell* wet-lab panels, paired with mean ± s.e.m. (or s.d.) and **all individual dots plotted**.
- **`5 ≤ n ≤ 50` per group** (most clinical cohorts, patient-level scRNA summaries): Run Shapiro–Wilk on each group (and Levene's test). If **any group fails normality (`p_shapiro < 0.05`)** OR the metric is a bounded proportion/score/gene expression → **use non-parametric Wilcoxon / Kruskal–Wallis**. If all groups are normal → **Welch's *t*-test / Welch ANOVA**.
- **`n > 50` per group**: Visual skewness/outlier check; prefer **Mann–Whitney U / Kruskal–Wallis** for skewed omics/clinical variables and **Welch's *t*** for symmetric continuous measurements.
- **Note**: In Python (`sci_plot_py.py`) and R (`sci_plot_r.R`), `auto_stat_test()` automates this check and explicitly returns the chosen test name, exact two-sided *p*-value, and formatted `fmt_p()` string.
