# Clinical Cohort, Oncology Trial, Survival & Translational ML Playbook

Empirically calibrated from **747 Tier-1 CNS 2025–2026 tumor papers** (*Nature Medicine*, *Cancer Cell*, *Nature*, *Cell*, *Nature Cancer*, *Cancer Discovery*, *Cell Reports Medicine*, *Science Translational Medicine*).

---

## 1. Clinical & Translational Oncology Data-to-Chart Matrix (`2025–2026 Tier-1 CNS`)

| Clinical / Translational Readout | Tier-1 % (`25–26`) | Track A: Polished Conventional | Track B: Eye-Catching CNS Upgrade | Built-In Helper |
| :--- | :---: | :--- | :--- | :--- |
| **1. Time-to-Event Survival (OS / PFS / RFS / DFS)** | **`43.0%`** | Kaplan–Meier step curves (`1.05 pt`) + censored `+` ticks + `log-rank p` | **KM Curve + 95% CI Ribbons + Exact X-Aligned `Number at Risk` Table** + tucked-in legend & `HR (95% CI)` inset | `plot_km_survival()` |
| **2. Paired Pre-Treatment vs. On/Post-Treatment Biopsies** | **`15.9%`** | Boxplot + thin gray connecting lines + paired Wilcoxon *p* | **Dual-Color Up/Down Paired Slope + Flanking Hollow Half-Boxes** (coral-red increase / steel-blue decrease trajectories!) | `plot_paired_slope()` |
| **3. Multivariable Prognostic / Subgroup Regression** | **`14.2%`** | Plain HR dot-and-whisker | **Publication Cox Forest Plot** (alternating `#F6F6F6` zebra stripes + left aligned `Variable \| HR (95% CI) \| p` table + right log-scale HR squares & dashed `HR = 1`) | `plot_cox_forest()` |
| **4. Patient Cohort Multi-Omic Inventory & Genomic Landscape** | **`11.1%` / `8.2%`** | Summary bar chart of clinical features | **Multi-Track Clinical Cohort OncoPrint** (Top TMB/assay count bar + Main mutation/assay tiles with split double-hits + Right `%` bar + Bottom Stage/Response/Subtype/MSI/Sex ribbons!) | `plot_clinical_oncoprint()` |
| **5. Best Radiographic Tumor Change (RECIST v1.1)** | **`8.4%`** | Bar plot sorted by response | **RECIST Waterfall Plot** (descending `% change` bars colored by `CR/PR/SD/PD` + dashed `+20%` PD and `−30%` PR cutoffs + optional bottom biomarker strip) | `plot_waterfall_response()` |
| **6. Predictive Biomarker / ML Classifier & Interpretability** | **`7.5%`** | Single ROC curve | **Multi-Model ROC + AUROC (95% CI) Inset** paired with **SHAP Summary Beeswarm** (colored high→low feature value!) & Calibration Curve | `plot_roc_curves()` / `plot_shap_beeswarm()` |
| **7. Multi-Cohort / Multi-Signature Overlaps (`≥ 4` sets)** | **`3.7%`** | Venn diagram (NEVER use for `≥ 4` sets!) | **UpSet Intersection Plot** (sorted vertical intersection bars + connected dot-matrix grid + left horizontal set-size bars!) | `plot_upset()` |
| **8. Longitudinal Trial Duration & Clinical Events** | **`2.9%`** | Summary median PFS bar | **Oncology Swimmer Plot** (horizontal patient duration bars + overlaid `CR`/`PR` onset triangles, `PD` diamonds, Ongoing `→` arrows + vertical landmark dashed line!) | `plot_swimmer_plot()` |
| **9. Multi-Signature Functional Profile Across Subtypes** | **`1.7%`** | Small grouped barplot | **Multi-Axis Radar / Spider Profile** (`5–8` hallmark axes, `2–3` overlaid translucent polygon series with `alpha = 0.16`!) | `plot_radar_profile()` |
| **10. Individual Longitudinal Tumor Burden / ctDNA Kinetics** | **`1.2%`** | Group mean line plot | **Clinical Trial Spider Plot** (patient `% change` trajectories over months colored by Responder vs Non-Responder + `+20%`/`−30%` RECIST lines!) | `plot_spider_response()` |

---

## 2. Clinical Oncology Visual Craft & Reporting Standards

1. **Exact Standardized Clinical Color Semantics (Lock Across All Figures in a Manuscript!)**:
   - **RECIST Response**: `CR` = `#2A9D8F` (Teal), `PR` = `#3B6EA5` (Steel Blue), `SD` = `#E08E45` (Warm Amber), `PD` = `#C84B4B` (Brick Red).
   - **Binary Clinical Outcome**: `Responder / Benefit` = `#2A9D8F` (Teal), `Non-responder / Resistance` = `#C84B4B` (Brick Red).
   - **OncoPrint Alterations**: `Missense` = `#3B6EA5`, `Truncating / Nonsense` = `#1A1A1A`, `Amplification` = `#C84B4B`, `Deep Deletion` = `#6BAED6`, `Splice / Fusion` = `#7A5195`, `Wild-type / NA` = `#EFEFEF`.
   - **Clinical Stage (Sequential)**: `Stage I` `#D9E6F2` → `Stage II` `#92B8DB` → `Stage III` `#4A80B4` → `Stage IV` `#1D4673`.
2. **Kaplan–Meier Non-Negotiables**:
   - Always plot `0.0–1.0` survival probability (or `0–100%` — NEVER mix `%` with `0–1` ticks!).
   - Always include **censored tick marks** (`|`), **`two-sided log-rank p`** in lower-left, and an **exact X-aligned `Number at risk` table** underneath!
3. **OncoPrint Sorting & Hierarchy**:
   - Sort genes top-to-bottom by descending alteration frequency; sort patients left-to-right by **mutual-exclusivity waterfall order** (or grouped by Responder vs Non-Responder / Molecular Subtype with clean white vertical gap separating clinical arms!).
   - Keep bottom clinical annotation strips slim (`~2.2–2.8 mm` tall per track) with crisp `0.4 pt` white cell gaps so 5–6 clinical tracks stay compact and razor-sharp.
4. **Forest Plot Readability**:
   - Always plot HR/OR on a **log scale** (so `0.5` and `2.0` are equidistant from `HR = 1.0`!).
   - Alternate subtle `#F7F7F7` / `#FFFFFF` zebra row backgrounds so the eye tracks effortlessly from the left text table (`Variable`, `HR (95% CI)`, `*p*`) across to the right whisker plot.
