# Do spatial tumor-microenvironment features add prognostic value in breast cancer?

A small, pre-specified comparison on the Jackson et al. 2020 (*Nature*) imaging mass cytometry (IMC) cohort. The comparison is between overall-survival models built from clinical variables alone, from clinical variables plus per-patient cell-type composition, and from both of those plus a small panel of spatial-interaction features.

**Short answer: no detectable added value on the primary endpoint.** On overall survival (N = 280, 79 deaths), adding composition and spatial features changed the out-of-fold C-index by less than 0.01, and every confidence interval includes zero. The clinical-only model reaches 0.741 and the full model 0.750. The spatial-over-composition gain is +0.002 (95% CI −0.005 to +0.008), and it averages −0.002 across 10 CV repeats. Breast-cancer-specific survival (a secondary sensitivity analysis) shows a small, unstable gain from spatial features, and one post-hoc variant of it passes a permutation test. Treat it as a hypothesis for an external cohort, not a finding (see [Robustness](#robustness-cv-repeats-and-a-permutation-null)).

Every number below is copied from files in `results/`. The tables are produced verbatim by `scripts/report_tables.py` → `results/tables.md`, and `tests/test_readme_numbers.py` checks that each table row in this README appears in that file.

---

## 1. Question

On one cohort, does each added data layer improve out-of-fold discrimination (Harrell's C) for overall survival? The three models are:

1. `clinical`: clinical variables only
2. `clinical+composition`: plus per-patient cell-type composition
3. `clinical+composition+spatial`: plus per-patient spatial-interaction features

## 2. Data

**Source.** Zenodo record [10.5281/zenodo.4607374](https://doi.org/10.5281/zenodo.4607374), "The Single-Cell Pathology Landscape of Breast Cancer". The earlier version, 3518284, holds a subset of the same files. `scripts/download_data.py` lists the record and verifies MD5 checksums. It downloads only:

| File | Size | Contents used |
|---|---|---|
| `singlecell_cluster_labels.zip` | 7.9 MB | `Basel_metaclusters.csv` (cell id → metacluster 1–27), `Metacluster_annotations.csv` |
| `singlecell_locations.zip` | 23.2 MB | `Basel_SC_locations.csv` (cell id, core, `Location_Center_X/Y`) |
| `SingleCell_and_Metadata.zip` (36.8 GB) | **not downloaded** | only `Data_publication/BaselTMA/Basel_PatientMetadata.csv` (0.16 MB), extracted by HTTP range requests (a few MB transferred, CRC-checked by `zipfile`) |

**Cohort: Basel only.** Basel has a standalone patient-metadata table with overall survival. Zurich clinical data exists only inside a 3.7 GB histoCAT MATLAB session, so Zurich is not used.

**Filtering.** 376 cores → 289 tumor cores (`diseasestatus == "tumor"`; adjacent normal-tissue cores excluded) → 281 patients → **N = 280** after dropping one patient with `OSmonth = 0`. The final cohort has 288 images and 762,776 labelled cells. 8 patients have two tumor cores. Clinical fields were verified to be identical across a patient's cores.

**Target.** Time = `OSmonth`. The event is `Patientstatus` ∈ {"death", "death by primary disease"}, giving **79 events**. "alive" and "alive w metastases" count as censored. Median follow-up time across all patients is 74 months.

**Secondary endpoint (sensitivity, agreed before modelling).** Breast-cancer-specific survival (DSS) counts only "death by primary disease" as an event (**57 events**); other deaths are censored.

**Patient-level aggregation.** For composition, cells from all of a patient's tumor cores are pooled before counting. Spatial features are computed per image, then averaged per patient with weights by cell count (NaN-aware).

## 3. Features

**Clinical (7).** Age, tumor size (mm), histological grade (1–3, ordinal), pN (0–3; "x" → missing), and ER, PR, HER2 (positive = 1). Missing values (pN 13, ER 1, PR 6) are median-imputed inside each training fold. Two choices deviated from the original plan, both forced by the data inside CV folds:
- *Grade is ordinal rather than one-hot.* Grade 1 has only 2 disease-specific deaths in 38 patients, so in some training folds the grade dummies are quasi-separated and their unpenalized coefficients diverge.
- *pM is not modelled.* Only 7/280 patients are M1, and some inner training folds contain a single M1 patient, so the unpenalized coefficient cannot be estimated.

pT was left out as redundant with tumor size, and treatment because it is decided after baseline.

**Composition (27).** Counts of the 27 published metaclusters per patient. The CLR transform adds a pseudocount of 0.5 cells, closes each row, takes logs, and centres each row on its own mean. CLR is row-wise and uses no information from other patients. The standardisation that follows is fit on training folds only.

**Spatial (11).** Per image:
- **Graph.** k-nearest neighbours, k = 8, on cell centroids, symmetrised to an undirected graph. IMC pixels are 1 µm, so coordinates are in µm. I chose kNN over a 20 µm radius because it adapts to local density and leaves no cell isolated.
- **Coarse classes** (from `Metacluster_annotations.csv`):
  - Tumor (all 14 "Tumor"-class metaclusters, including myoepithelial)
  - T (3, 5)
  - B (1 B cell, 2 T & B cells)
  - Macrophage (4, 6)
  - Endothelial (7)
  - Stroma (8–13)
- **Neighbour enrichment, 9 pairs.** Tumor–Tumor, Tumor–T, Tumor–Macrophage, Tumor–B, Tumor–Stroma, Tumor–Endothelial, T–Macrophage, T–B, Stroma–T. The statistic is the number of undirected edges between the two classes, converted to a z-score against 200 within-image label permutations (graph fixed, labels shuffled). It is set to missing if either class has fewer than 5 cells in the image; missing values are median-imputed in-fold.
- **`tumor_immune_nbr_frac`.** The mean fraction of immune cells (T, B, macrophage) among each tumor cell's 8 nearest neighbours.
- **`tumor_immune_mixing`.** A bounded variant of the Keren et al. 2018 mixing score: tumor–immune edges / (tumor–immune + immune–immune edges).

The enrichment statistic is unit-tested against a hand-computed exact null on a 4-node graph and against the closed-form permutation mean.

**Sensitivity variant (added after the first results; see §7).** The enrichment z-scores turned out to track image size. A permutation z grows roughly with the square root of the edge count, and the correlation with log(cell count) reaches |ρ| = 0.55 (table in §5). I therefore also ran the full model with **log₂((observed + 1)/(expected + 1))** for the 9 pairs, which is much less size-dependent.

## 4. Models and validation

- **Clinical model:** ridge Cox (`CoxPHSurvivalAnalysis`). α is chosen from 20 log-spaced values from 10⁻² to 10³.
- **Models with omics blocks:** elastic-net Cox (`CoxnetSurvivalAnalysis`, l1_ratio = 0.5). Clinical columns are unpenalized (`penalty_factor = 0`), so these models nest the clinical model. α is chosen from a 20-point path (α_max down to 0.01·α_max) derived from the training data of each outer fold.
- **Tuning effort is identical:** one penalty, 20 candidate values, inner 5-fold CV stratified by event, selected by mean Harrell's C.
- **Baseline:** random Gaussian risk scores, as a sanity check that C ≈ 0.5.
- **Outer CV:** 5-fold stratified by event, repeated with seeds 0, 1, 2, using the same folds for every model. Rows are patients, so all splits are patient splits.
- **Leakage control:** imputation, scaling and the α path are fit inside each training fold; CLR is row-wise. `tests/test_leakage_and_splits.py` spies on every fit call inside the nested CV and asserts that none ever sees an outer-test patient.
- **Metric:** Harrell's C on the **pooled** out-of-fold predictions of each repeat, then averaged over the 3 repeats.
- **Confidence intervals:** 1,000 patient bootstrap resamples, with the same resamples for every model, so model differences are paired. Percentile 95% CIs.
- **Robustness checks:** (a) 10 CV repeats (seeds 0–9) for the paired spatial-vs-composition delta; (b) a permutation null in which the spatial block is shuffled across patients, then the full nested CV is rerun (20 permutations × 3 seeds). The bootstrap above holds the out-of-fold predictions fixed and so ignores CV-partition and refitting variability; these two checks address that.

## 5. Results

### Overall survival (primary)

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.741 | 0.687 to 0.792 | 0.735, 0.755, 0.733 |
| clinical+composition | 0.748 | 0.690 to 0.801 | 0.740, 0.760, 0.743 |
| clinical+composition+spatial | 0.750 | 0.693 to 0.803 | 0.742, 0.760, 0.747 |
| clinical+composition+spatial (log O/E, sensitivity) | 0.744 | 0.688 to 0.796 | 0.742, 0.740, 0.750 |
| random baseline | 0.503 | 0.465 to 0.540 | 0.564, 0.484, 0.461 |

Paired bootstrap differences:

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.007 | -0.017 to +0.029 | 0.67 |
| clinical+composition+spatial minus clinical+composition | +0.002 | -0.005 to +0.008 | 0.70 |
| clinical+composition+spatial minus clinical | +0.009 | -0.014 to +0.031 | 0.74 |
| clinical+composition+spatial (log O/E, sensitivity) minus clinical+composition | -0.004 | -0.014 to +0.006 | 0.24 |

![OS C-index](figures/os_cindex.png)

### Breast-cancer-specific survival (sensitivity)

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.721 | 0.660 to 0.777 | 0.734, 0.737, 0.693 |
| clinical+composition | 0.723 | 0.660 to 0.780 | 0.753, 0.707, 0.710 |
| clinical+composition+spatial | 0.738 | 0.677 to 0.792 | 0.759, 0.739, 0.715 |
| clinical+composition+spatial (log O/E, sensitivity) | 0.751 | 0.691 to 0.803 | 0.768, 0.748, 0.736 |
| random baseline | 0.505 | 0.459 to 0.553 | 0.561, 0.464, 0.489 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.002 | -0.038 to +0.040 | 0.54 |
| clinical+composition+spatial minus clinical+composition | +0.014 | +0.000 to +0.031 | 0.97 |
| clinical+composition+spatial minus clinical | +0.017 | -0.021 to +0.052 | 0.81 |
| clinical+composition+spatial (log O/E, sensitivity) minus clinical+composition | +0.028 | +0.008 to +0.050 | 1.00 |

![DSS C-index](figures/dss_cindex.png)

### Robustness: CV repeats and a permutation null

OS:

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| clinical+composition+spatial | +0.002 | -0.002 ± 0.006 [-0.012, +0.007] | 4/10 | +0.002 ± 0.003 (+0.005) | 0.667 |
| clinical+composition+spatial (log O/E, sensitivity) | -0.004 | +0.001 ± 0.010 [-0.020, +0.016] | 7/10 | +0.001 ± 0.003 (+0.006) | 0.952 |

DSS:

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| clinical+composition+spatial | +0.014 | +0.002 ± 0.014 [-0.021, +0.032] | 7/10 | +0.007 ± 0.008 (+0.023) | 0.143 |
| clinical+composition+spatial (log O/E, sensitivity) | +0.028 | +0.011 ± 0.017 [-0.018, +0.041] | 8/10 | +0.005 ± 0.008 (+0.016) | 0.048 |

The permutation p is (1 + number of null deltas ≥ observed)/(1 + 20). With 20 permutations its smallest possible value is 1/21 ≈ 0.048.

### Kaplan–Meier (median split of out-of-fold risk)

The best OS model by C-index (`clinical+composition+spatial`) separates the groups strongly (log-rank p = 6.5e-09). However, the clinical-only model separates them at least as strongly (p = 8.3e-10), and the two models put 94% of patients in the same risk half. The separation in this figure is therefore clinical signal, not spatial signal. Risk is the mean out-of-fold rank across the 3 CV repeats.

![KM best model](figures/os_km_best_model.png)
![KM clinical model](figures/os_km_clinical_model.png)

### Spatial coefficients (exploratory)

Standardised log-hazard coefficients of the spatial features in `clinical+composition+spatial`, averaged over the 15 outer-fold fits, with the share of folds in which the elastic net kept each feature. These are fits with correlated inputs and shrinkage, so they **must not be read as effect estimates**.

OS:

| Spatial feature | mean coef | SD | selected in folds |
|---|---|---|---|
| `tumor_immune_mixing` | +0.045 | 0.077 | 33% |
| `enrich_T__Macrophage` | -0.037 | 0.062 | 33% |
| `enrich_Tumor__Macrophage` | +0.011 | 0.034 | 13% |
| `enrich_Tumor__T` | +0.007 | 0.025 | 13% |
| `enrich_Stroma__T` | +0.007 | 0.025 | 7% |
| `enrich_Tumor__Endothelial` | +0.003 | 0.011 | 7% |
| `enrich_Tumor__Tumor` | +0.000 | 0.000 | 0% |
| `enrich_Tumor__B` | +0.000 | 0.000 | 0% |
| `enrich_Tumor__Stroma` | +0.000 | 0.000 | 0% |
| `enrich_T__B` | +0.000 | 0.000 | 0% |
| `tumor_immune_nbr_frac` | +0.000 | 0.000 | 0% |

DSS:

| Spatial feature | mean coef | SD | selected in folds |
|---|---|---|---|
| `enrich_Tumor__T` | +0.091 | 0.066 | 87% |
| `enrich_T__Macrophage` | -0.069 | 0.074 | 80% |
| `tumor_immune_mixing` | +0.064 | 0.070 | 60% |
| `enrich_Tumor__Tumor` | -0.017 | 0.034 | 33% |
| `enrich_Tumor__Macrophage` | +0.011 | 0.019 | 33% |
| `enrich_Tumor__B` | +0.008 | 0.030 | 7% |
| `tumor_immune_nbr_frac` | +0.007 | 0.019 | 13% |
| `enrich_Tumor__Endothelial` | +0.000 | 0.000 | 0% |
| `enrich_Tumor__Stroma` | +0.000 | 0.000 | 0% |
| `enrich_Stroma__T` | +0.000 | 0.000 | 0% |
| `enrich_T__B` | +0.000 | 0.000 | 0% |

![OS spatial coefficients](figures/os_spatial_coefficients.png)
![DSS spatial coefficients](figures/dss_spatial_coefficients.png)

### Diagnostics

Penalty choice and sparsity across the 15 outer folds:

| Endpoint | Model | selected penalty α: median [min, max] | non-zero omics coefs: median [min, max] |
|---|---|---|---|
| OS | clinical | 88.6 [0.379, 1e+03] | 0 [0, 0] |
| OS | clinical+composition | 0.104 [0.0247, 0.146] | 1 [0, 15] |
| OS | clinical+composition+spatial | 0.128 [0.034, 0.183] | 1 [0, 19] |
| OS | clinical+composition+spatial (log O/E, sensitivity) | 0.102 [0.0262, 0.196] | 3 [0, 22] |
| DSS | clinical | 88.6 [0.0616, 1e+03] | 0 [0, 0] |
| DSS | clinical+composition | 0.0532 [0.00914, 0.12] | 8 [0, 21] |
| DSS | clinical+composition+spatial | 0.0566 [0.0256, 0.155] | 10 [0, 20] |
| DSS | clinical+composition+spatial (log O/E, sensitivity) | 0.0511 [0.0256, 0.146] | 11 [0, 21] |

Image-level dependence of spatial features on image size:

| Feature | Spearman ρ with log(cells), z-score version | log O/E version |
|---|---|---|
| `enrich_Tumor__Tumor` | +0.55 | -0.29 |
| `enrich_Tumor__T` | -0.44 | -0.03 |
| `enrich_Tumor__Macrophage` | -0.30 | -0.02 |
| `enrich_Tumor__B` | -0.21 | +0.06 |
| `enrich_Tumor__Stroma` | -0.41 | -0.21 |
| `enrich_Tumor__Endothelial` | -0.51 | -0.10 |
| `enrich_T__Macrophage` | +0.30 | +0.30 |
| `enrich_T__B` | -0.03 | -0.06 |
| `enrich_Stroma__T` | +0.39 | +0.38 |
| `tumor_immune_nbr_frac` | -0.17 | (same feature) |
| `tumor_immune_mixing` | +0.19 | (same feature) |

## 6. Interpretation

- **Primary endpoint (OS): a null result.** Clinical variables alone give C ≈ 0.74. Composition adds +0.007 and spatial features a further +0.002, both with CIs spanning zero. Over 10 CV repeats the spatial delta averages −0.002, and it sits inside the permutation null (p = 0.67). In the OS composition and spatial models the elastic net usually keeps only about one omics feature (median), which is consistent with little extra signal.
- **Secondary endpoint (DSS): weak, unstable hints.** With z-score features, the +0.014 gain seen on the 3 main CV repeats shrinks to +0.002 ± 0.014 over 10 repeats and is not separable from the permuted-feature null (p = 0.14). Notably, even permuted spatial blocks raise the DSS C-index by +0.007 on average, which shows how noisy the composition baseline is at 57 events. The log O/E variant does better: +0.011 ± 0.017 over 10 repeats, 8/10 repeats positive, permutation p = 0.048. That variant was added after seeing the first results, the endpoint is secondary, and about 8 paired comparisons were examined, so this does not survive any reasonable multiplicity correction. It is still the most interesting lead: spatial organisation of immune cells relative to tumor might carry disease-specific information that is diluted in OS by the 22 non-cancer deaths.
- **Kaplan–Meier separation is not evidence for spatial features.** The clinical-only model gives an equally strong split.
- **Exploratory coefficients.** For DSS, the features kept most often were Tumor–T enrichment (positive coefficient, i.e. *higher* hazard with more tumor–T contact than chance) and T–Macrophage enrichment (negative coefficient). The Tumor–T z-score correlates with image size (ρ = −0.44), so this sign may partly reflect core size rather than biology. I would not interpret it without a size-robust replication.

## 7. Deviations from the original plan, stated plainly

1. Grade was encoded as ordinal instead of one-hot, and pM was dropped. Both caused numerical divergence of unpenalized clinical coefficients within CV folds (§3). The decision was made from fit failures, not from C-index values. For transparency: a first OS run with one-hot grade and pM did complete before the DSS fits failed. Its output was printed to the console and not saved, so no numbers from it are reported here. It showed the same pattern: small, CI-spanning gains.
2. The log O/E spatial variant, the CV-repeat analysis, and the permutation null were **added after** seeing the first DSS result, because the z-score size dependence and the bootstrap's fixed-prediction assumption made me distrust it. They are reported as sensitivity and robustness analyses. They are not a re-definition of the primary model.
3. The DSS endpoint was agreed before modelling, as a sensitivity analysis.

## 8. Reproduction

```bash
uv sync                                   # Python 3.12, pinned deps (pyproject.toml + uv.lock)
uv run python scripts/download_data.py    # ~31 MB download + 0.16 MB range-extracted; checksums verified
uv run python scripts/run.py              # all results/ and figures/; ~5 min on a 20-thread laptop CPU (exact: runtime_seconds in results/summary.json)
uv run python scripts/report_tables.py    # results/tables.md (the tables above)
uv run pytest                             # 12 tests
```

**Determinism.** A clean rerun (`scripts/run.py --recompute-spatial`, which also recomputes the cached image-level spatial features) reproduced `results/tables.md` exactly, apart from the runtime line. Raw out-of-fold risk scores differed by at most ~4e-6 between runs, which is floating-point noise from multithreaded numerics.

`uv run python scripts/download_data.py --list` lists the Zenodo record and the 36.8 GB archive's members (reading only its central directory) without downloading anything. Configuration (seeds, k, permutations, grids) is in `CONFIG` in `scripts/run.py` and is saved to `results/summary.json` together with package versions and the runtime. Raw per-fold records (seed, fold, selected α, inner C, test-fold C, number of non-zero coefficients, test patient IDs) are in `results/{os,dss}_per_fold.csv`, and out-of-fold risks are in `results/{os,dss}_oof_risk.csv`.

**Tests** (`tests/`):
- **Leakage.** Imputer, scaler and CLR outputs depend only on training rows. A spy on every fit inside the nested CV confirms no fit sees an outer-test patient.
- **Patient-level splits.** Outer folds partition patients, and multi-image patients collapse to one row before CV.
- **Neighbour enrichment.** Exact hand-worked null on a 4-node path graph, the closed-form permutation mean, and segregated vs. checkerboard layouts.
- **Toy end-to-end.** Simulated images in which T-cell infiltration drives survival at fixed composition; the spatial model must beat the clinical and composition models by more than 0.1 C.
- **README integrity.** Every table row in this file appears in `results/tables.md`.

Repository layout: `scripts/download_data.py`, `scripts/run.py`, `scripts/report_tables.py`, `src/spatialsurv/{data,features,models,cv,eval,robustness,plots}.py`, `tests/`, `results/`, `figures/`.

## 9. Caveats and limitations

- **Single cohort, no external validation.** Nested CV estimates internal performance only. Zurich, the companion cohort, lacks accessible survival data here.
- **Overall survival includes non-cancer deaths** (22 of 79). These dilute any tumor-biology signal, which is why DSS was run. DSS treats competing deaths as censored, so it estimates a cause-specific hazard, not a cumulative incidence.
- **Small n relative to features.** 280 patients and 79 (OS) or 57 (DSS) events, against up to 7 + 27 + 11 = 45 features. Penalisation helps, but C-index differences of about 0.01 are below what this cohort can resolve. The 95% CIs of individual models are about ±0.05 wide.
- **One small TMA core per patient (rarely two)**, so the spatial features describe a tiny sample of each tumor and are noisy. Pairs involving rare classes are often undefined; Tumor–B is missing for 108/280 patients and median-imputed.
- **Spatial features are exploratory and partly confounded with image size** (z-score version). The panel of 11 was chosen a priori for interpretability; other summaries (e.g. cellular neighbourhoods, distance-based statistics, the paper's own community features) might behave differently.
- **Harrell's C** depends on the censoring distribution and summarises discrimination only. Calibration was not assessed.
- **Ridge α for the clinical model sometimes sits at the grid edge** (10³). Large ridge penalties mostly rescale the linear predictor, so the rank-based C-index is insensitive to them. I did not re-tune the grid after seeing results.
- **The cell-type labels are the published metaclusters**, fit once on all Basel cells, including future test patients. This is unsupervised and outcome-blind, so it is not outcome leakage, but it is a shared upstream step that a strictly prospective pipeline would not have.

## 10. Next steps

- **External validation** of the one lead: DSS with size-robust (log O/E) immune–tumor enrichment, pre-registered in a second IMC/CODEX cohort, e.g. the METABRIC IMC cohort of Danenberg et al. 2022, which has breast-cancer-specific survival.
- **Size-robust, multi-scale spatial statistics.** Radius-based graphs at several scales, Ripley's cross-K/L, and cellular-neighbourhood fractions, with image-size adjustment built in.
- **Multimodal integration** with clinical data, bulk expression or H&E-derived features, comparing late fusion (stacked risk scores) with joint penalised models, and judged with the same nested-CV-plus-permutation-null protocol used here.
- **Competing-risks modelling** (Fine–Gray or cause-specific Cox for both causes) instead of censoring non-cancer deaths.
- **Calibration and decision-curve analysis** if any layer clears the discrimination bar.

---

*Data: Jackson, H.W. et al. "The single-cell pathology landscape of breast cancer." Nature 578, 615–620 (2020).*
