# Do spatial tumor-microenvironment features add prognostic value in breast cancer?

A methods note on two imaging mass cytometry (IMC) breast cancer cohorts: Basel (Jackson et al. 2020, *Nature*) and METABRIC (Danenberg et al. 2022, *Nature Genetics*). The analysis asks whether per-patient cell-type composition, and then spatial-interaction features, improve out-of-fold survival prediction **beyond standard clinical variables**.
- **Part I:** discovery in Basel.
- **Part II:** a pre-registered external replication in METABRIC.
- **Part III:** re-typing the METABRIC cells with **CELESTA-Lite**, a Python port of the CELESTA spatial cell-typing engine, to see whether better-informed cell labels change the answer.
- **Part IV:** a separate, descriptive question. Using the official CELESTA R package and the Plevritis lab's co-location quotient (CLQ) code, do tumour subtype and genotype shape spatial organisation beyond cell composition? No survival outcomes are involved. See [Part IV](#part-iv-tumour-genotype-and-spatial-organisation).

CPU only. Each part reruns from one command.

| Basel (discovery) | METABRIC (pre-registered replication) |
|:---:|:---:|
| ![Basel headline](figures/headline.png) | ![METABRIC headline](figures/part2_metabric_headline.png) |

**Findings**

1. **Spatial features add no detectable prognostic value, in either cohort.**
   - **Basel**, overall survival (79 deaths): +0.002 C over composition (95% CI −0.005 to +0.008).
   - **METABRIC** (231 deaths and 156 breast-cancer deaths within 15 years): every spatial variant lands within ±0.006 of the composition model. Every 95% CI upper bound is at most +0.005, so gains larger than that are ruled out (§P2.2).
   - Composition doesn't help either: −0.007 in METABRIC.
2. **The one lead from Basel failed its pre-registered test.** Basel suggested a gain for size-robust spatial enrichment on breast-cancer-specific survival.
   - Hypothesis H1 and its decision rule were committed to git before any METABRIC outcome was linked.
   - Result in METABRIC: ΔC = −0.004 (95% CI −0.009 to +0.001), permutation p = 0.810. **Verdict: refuted.**
3. **Two popular "spatial" summaries are mostly composition in disguise.**
   - Across Basel images, the tumor–immune mixing score has ρ = −0.89 with immune-cell fraction, and the immune fraction among tumor cells' neighbours has ρ = +0.92.
   - Permutation z-scores track image size (|ρ| up to 0.55).
   - Removing composition and size from the spatial features, in-fold, does not reveal hidden signal in METABRIC (+0.002 for OS, −0.004 for DSS).
4. **The signal is clinical, and it transfers.**
   - A clinical Cox model trained on Basel reaches C = 0.708 on METABRIC, matching a clinical model cross-validated within METABRIC itself (0.704).
   - Adding Basel-trained composition *hurts* transfer (−0.025, 95% CI −0.044 to −0.006).
   - Calibration slopes show the penalised omics models are overfit (0.76–0.85 in METABRIC), while the heavily ridge-shrunk clinical model is under-confident (1.37–1.64).
5. **Re-typing the cells with CELESTA-Lite changes the cell labels, not the conclusion.**
   - The port's spatial step improves agreement with the published labels for the ambiguous cells it assigns (0.589 vs. 0.398 for markers alone).
   - Overall agreement is moderate (κ = 0.476).
   - Spatial features computed from CELESTA-Lite labels add nothing either (OS +0.000, DSS −0.006).

**What guards the result:**
- Nested CV, split by patient.
- Tests that spy on every model fit and fail if any fit sees a held-out patient.
- A shared patient bootstrap, plus 10 CV repeats and a permutation null.
- A pre-registered external replication, with amendments time-stamped in git.
- Dry runs on outcome-shuffled data, so bugs were fixed without seeing results.
- A test that every table row in this README is reproduced verbatim from saved results.

Deviations are listed in §7 (part I), in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md) amendments (part II), and in §IV.7 and [AMBIGUITIES_GENOTYPE.md](AMBIGUITIES_GENOTYPE.md) (part IV).

---

## 1. Question and relation to prior work

On a single cohort, does each added layer improve out-of-fold discrimination (Harrell's C) of overall survival? The layers are clinical → + cell-type composition → + spatial interactions.

This question is narrower than the claims of the source paper. **Jackson et al. 2020** profiled 35 markers in 720 IMC images from 352 patients, 281 of them with long-term survival data (the Basel cohort used here). They defined 18 single-cell pathology (SCP) subgroups from cell phenotypes and their community organisation, and these subgroups differed in outcome. That is a group-level association discovered on the whole cohort. It does not measure how much out-of-fold discrimination these data add once standard clinical variables are already in the model, which is the question here. A null result here is therefore not in tension with theirs.

The mixing score adapts **Keren et al. 2018** (MIBI-TOF, 41 triple-negative breast cancers), who described "mixed" vs. "compartmentalized" tumor–immune architecture and linked compartmentalization to survival. In this cohort my variant of the score is largely explained by immune content (finding 3), so it should be read alongside composition, not instead of it. **Danenberg et al. 2022** (693 METABRIC tumors) found recurrent multicellular TME structures, including a "suppressed expansion" structure that predicted poor outcome in ER-positive disease. That suggests prognostic spatial information may require richer descriptors and more patients than this pre-specified panel of 11 features on 280 patients.

## 2. Data

**Source.** Zenodo [10.5281/zenodo.4607374](https://doi.org/10.5281/zenodo.4607374); the earlier version, 3518284, holds a subset of the same files. `scripts/download_data.py` lists the record, verifies MD5 checksums, and fetches only:

| File | Size | Contents used |
|---|---|---|
| `singlecell_cluster_labels.zip` | 7.9 MB | `Basel_metaclusters.csv` (cell id → metacluster 1–27), `Metacluster_annotations.csv` |
| `singlecell_locations.zip` | 23.2 MB | `Basel_SC_locations.csv` (cell id, core, `Location_Center_X/Y`) |
| `SingleCell_and_Metadata.zip` (36.8 GB) | **not downloaded** | only `Data_publication/BaselTMA/Basel_PatientMetadata.csv` (0.16 MB), extracted with HTTP range requests (a few MB transferred, CRC-checked) |

**Cohort: Basel only.** Basel has a standalone patient-metadata table with survival. Zurich clinical data exists only inside a 3.7 GB histoCAT MATLAB session, so Zurich is not used.

**Filtering.** 376 cores → 289 tumor cores (`diseasestatus == "tumor"`; normal-tissue cores excluded) → 281 patients → **N = 280** after dropping one patient with `OSmonth = 0`. The final cohort has 288 images and 762,776 labelled cells. 8 patients have two tumor cores, and their clinical fields are identical across cores.

**Endpoints.**
- **Overall survival (OS, primary).** Time = `OSmonth`; event = any death, **79 events**. "alive" and "alive w metastases" are censored. Median follow-up time across all patients is 74 months.
- **Breast-cancer-specific survival (DSS, secondary; agreed before modelling).** Only "death by primary disease" counts, **57 events**; other deaths are censored.

**Patient-level aggregation.** For composition, cells from all of a patient's tumor cores are pooled. Spatial features are computed per image and averaged per patient, weighted by cell count.

## 3. Features

![Example cores at similar immune content but different tumor-immune mixing](figures/example_cores.png)

*Example cores: one dot per cell. Tumor (dark grey), stroma (light grey) and endothelium are context; immune classes are coloured. All three cores come from images between the 40th and 60th percentile of immune fraction, so the figure contrasts **arrangement** at similar **amount**.*

**Clinical (7).** Age, tumor size (mm), grade (1–3, ordinal), pN (0–3; "x" → missing), and ER, PR, HER2 (positive = 1). Missing values (pN 13, ER 1, PR 6) are median-imputed inside each training fold. pM and one-hot grade were dropped for numerical reasons (§7). pT was left out as redundant with tumor size, and treatment because it is decided after baseline.

**Composition (27).** Counts of the 27 published metaclusters per patient, CLR-transformed: add a pseudocount of 0.5 cells, close each row, log, and centre each row on its own mean. CLR is row-wise; the standardisation that follows is fit on training folds only.

**Spatial (11), pre-specified.** Per image:
- **Graph.** Undirected kNN graph, k = 8, on cell centroids. IMC pixels are 1 µm. I chose kNN over a 20 µm radius because it adapts to local density and leaves no cell isolated.
- **Coarse classes** (from the metacluster annotations):
  - Tumor (all 14 "Tumor"-class clusters)
  - T (3, 5)
  - B (1 B cell, 2 T & B cells)
  - Macrophage (4, 6)
  - Endothelial (7)
  - Stroma (8–13)
- **Neighbour enrichment, 9 pairs.** Tumor–Tumor, Tumor–T, Tumor–Macrophage, Tumor–B, Tumor–Stroma, Tumor–Endothelial, T–Macrophage, T–B, Stroma–T. The statistic is the undirected edge count between two classes, as a z-score against 200 within-image label permutations (graph fixed, labels shuffled, which preserves composition). It is missing if either class has fewer than 5 cells; missing values are imputed in-fold.
- **`tumor_immune_nbr_frac`.** The mean fraction of immune cells among each tumor cell's 8 neighbours.
- **`tumor_immune_mixing`.** Tumor–immune edges / (tumor–immune + immune–immune edges), a bounded variant of Keren et al.

The enrichment code is unit-tested against a hand-computed exact null and the closed-form permutation mean.

**Sensitivity variant (post hoc, §7).** z-scores grow with the number of edges in an image, so the full model was also fit with **log₂((observed + 1)/(expected + 1))** for the 9 pairs.

## 4. Models and validation

- **Clinical model:** ridge Cox (`CoxPHSurvivalAnalysis`), with α chosen from 20 log-spaced values between 10⁻² and 10³.
- **Models with omics blocks:** elastic-net Cox (`CoxnetSurvivalAnalysis`, l1_ratio 0.5) with clinical columns **unpenalized**, so they nest the clinical model. α is chosen from a 20-point path derived from each outer-training set.
- **Tuning effort is identical:** one penalty, 20 candidates, inner 5-fold CV stratified by event, selected by mean Harrell's C.
- **Baseline:** random Gaussian risk scores, as a sanity check that C ≈ 0.5.
- **Outer CV:** 5-fold stratified by event, repeated with 3 seeds, using the same folds for every model. Rows are patients, so all splits are by patient.
- **Leakage control:** imputation, scaling and the α path are fit inside training folds. A test spies on every fit inside the nested CV and asserts that none sees an outer-test patient.
- **Metric:** Harrell's C on the **pooled** out-of-fold predictions of each CV repeat, averaged over repeats.
- **Uncertainty, three ways:**
  1. A 1,000-resample patient bootstrap shared by all models (paired differences, percentile CIs). It holds the out-of-fold predictions fixed.
  2. 10 CV repeats, which capture partition variability.
  3. A permutation null: shuffle the spatial block across patients and rerun the full nested CV (20 permutations × 3 seeds). This shows how much C a penalised model gains from 11 features that carry no outcome information.

## 5. Results

### 5.1 Overall survival (primary)

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.741 | 0.687 to 0.792 | 0.735, 0.755, 0.733 |
| clinical+composition | 0.748 | 0.690 to 0.801 | 0.740, 0.760, 0.743 |
| clinical+composition+spatial | 0.750 | 0.693 to 0.803 | 0.742, 0.760, 0.747 |
| clinical+composition+spatial (log O/E, sensitivity) | 0.744 | 0.688 to 0.796 | 0.742, 0.740, 0.750 |
| random baseline | 0.503 | 0.465 to 0.540 | 0.564, 0.484, 0.461 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.007 | -0.017 to +0.029 | 0.67 |
| clinical+composition+spatial minus clinical+composition | +0.002 | -0.005 to +0.008 | 0.70 |
| clinical+composition+spatial minus clinical | +0.009 | -0.014 to +0.031 | 0.74 |
| clinical+composition+spatial (log O/E, sensitivity) minus clinical+composition | -0.004 | -0.014 to +0.006 | 0.24 |

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| clinical+composition+spatial | +0.002 | -0.002 ± 0.006 [-0.012, +0.007] | 4/10 | +0.002 ± 0.003 (+0.005) | 0.667 |
| clinical+composition+spatial (log O/E, sensitivity) | -0.004 | +0.001 ± 0.010 [-0.020, +0.016] | 7/10 | +0.001 ± 0.003 (+0.006) | 0.952 |

### 5.2 Breast-cancer-specific survival (secondary)

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

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| clinical+composition+spatial | +0.014 | +0.002 ± 0.014 [-0.021, +0.032] | 7/10 | +0.007 ± 0.008 (+0.023) | 0.143 |
| clinical+composition+spatial (log O/E, sensitivity) | +0.028 | +0.011 ± 0.017 [-0.018, +0.041] | 8/10 | +0.005 ± 0.008 (+0.016) | 0.048 |

The permutation p is (1 + number of null deltas ≥ observed)/(1 + 20), so its smallest possible value is 0.048. The bootstrap CIs alone would have suggested a DSS gain (lower bounds +0.000 and +0.008). The CV-repeat and permutation analyses show how much of that came from the particular 3 CV partitions: even shuffled spatial features "gain" +0.007 on average for DSS.

### 5.3 Kaplan–Meier (median split of out-of-fold risk)

The best OS model by C-index (`clinical+composition+spatial`) gives log-rank p = 6.5e-09. The clinical-only model gives p = 8.3e-10, and the two put 94% of patients in the same risk half. Risk is the mean out-of-fold rank over the 3 CV repeats.

| Best model | Clinical only |
|:---:|:---:|
| ![KM best model](figures/os_km_best_model.png) | ![KM clinical model](figures/os_km_clinical_model.png) |

### 5.4 How large a gain could this cohort detect?

Approximate detectable ΔC at 80% power and two-sided α = 0.05, computed as 2.80 × the SD of the paired bootstrap difference. The last column also adds CV-repeat variability.

| Endpoint | Comparison | SD of paired bootstrap ΔC | detectable ΔC (80% power) | with CV-repeat variability |
|---|---|---|---|---|
| OS | clinical+composition minus clinical | 0.0120 | 0.034 | n/a |
| OS | clinical+composition+spatial minus clinical+composition | 0.0035 | 0.010 | 0.014 |
| DSS | clinical+composition minus clinical | 0.0200 | 0.056 | n/a |
| DSS | clinical+composition+spatial minus clinical+composition | 0.0079 | 0.022 | 0.032 |

For OS this design would usually have detected a spatial gain of about 0.01–0.014 C. The estimate is +0.002 with an upper CI bound of +0.008, so the data argue against a gain of that size, not merely fail to find one. For composition over clinical, only gains of about 0.03 or more were detectable, so smaller composition effects cannot be excluded.

### 5.5 Are the spatial features measuring something other than composition and image size?

Image-level Spearman correlations (288 images):

| Feature | ρ with log(cells): z-score | ρ with log(cells): log O/E | ρ with immune fraction: z-score | ρ with immune fraction: log O/E |
|---|---|---|---|---|
| `enrich_Tumor__Tumor` | +0.55 | -0.29 | +0.24 | +0.67 |
| `enrich_Tumor__T` | -0.44 | -0.03 | -0.43 | -0.24 |
| `enrich_Tumor__Macrophage` | -0.30 | -0.02 | -0.15 | -0.21 |
| `enrich_Tumor__B` | -0.21 | +0.06 | -0.47 | -0.42 |
| `enrich_Tumor__Stroma` | -0.41 | -0.21 | +0.38 | +0.12 |
| `enrich_Tumor__Endothelial` | -0.51 | -0.10 | +0.22 | -0.06 |
| `enrich_T__Macrophage` | +0.30 | +0.30 | +0.11 | -0.34 |
| `enrich_T__B` | -0.03 | -0.06 | +0.31 | -0.08 |
| `enrich_Stroma__T` | +0.39 | +0.38 | -0.47 | -0.65 |
| `tumor_immune_nbr_frac` | -0.17 | (same feature) | +0.92 | (same feature) |
| `tumor_immune_mixing` | +0.19 | (same feature) | -0.89 | (same feature) |

- **Mixing and neighbour fraction are near-composition proxies.** They are ratios of edge counts, and edge counts scale with how many immune cells there are. When immune cells are sparse, each one mostly touches tumor, so "mixing" is high. When immune cells are dense, they touch each other.
- **The z-scores are size-dependent.** The O/E version removes most of that dependence for tumor-centred pairs, but several O/E features correlate *more* strongly with immune fraction (e.g. Tumor–Tumor +0.67, Stroma–T −0.65).
- **Implication.** No feature in this panel is a clean "arrangement independent of amount" measure. A future panel should be checked against both axes before modelling, e.g. by residualising on composition and size, or by using statistics conditioned on composition.

<details>
<summary><b>5.6 Exploratory: spatial coefficients and model sparsity</b> (click to expand)</summary>

These are standardised log-hazard coefficients of the spatial features in `clinical+composition+spatial`, averaged over the 15 outer-fold fits, with the share of folds in which the elastic net kept each feature. The inputs are correlated and shrunk, so **these are not effect estimates**. For DSS, the most consistently kept feature is Tumor–T enrichment, with a *positive* coefficient (more tumor–T contact than chance goes with higher hazard). Its z-score correlates with both image size (−0.44) and immune fraction (−0.43), so I would not interpret the sign.

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

![OS spatial coefficients](figures/os_spatial_coefficients.png)
![DSS spatial coefficients](figures/dss_spatial_coefficients.png)
![OS C-index, all models](figures/os_cindex.png)
![DSS C-index, all models](figures/dss_cindex.png)

</details>

## 6. Interpretation

- **The clinical model is hard to beat.** Age, size, grade, nodes and receptor status give C ≈ 0.74. Adding 27 composition features and 11 spatial features does not measurably improve on it for OS, and the elastic net usually keeps only about one omics feature (median) when predicting OS.
- **This does not show that spatial organisation is irrelevant to breast cancer outcome.** It shows that *this* panel, on single small cores from 280 patients, carries no detectable *incremental, out-of-fold* information about OS beyond clinical variables. It also shows that two of its most intuitive features largely re-measure composition.
- **The DSS hint is worth one pre-registered test elsewhere, not a claim.** Spatial organisation could plausibly matter more for cancer-specific death than for all-cause death, which includes 22 non-cancer deaths here. But the signal is unstable across CV partitions and comes from a post-hoc variant.

## 7. Deviations from the original plan

1. **Grade is ordinal, not one-hot, and pM is not modelled.** Grade 1 has only 2 disease-specific deaths in 38 patients, and only 7/280 patients are M1. In some CV training folds this produced quasi-separation, and the unpenalized clinical coefficients diverged (Coxnet raised numerical errors). The change was made because of fit failures, not C-index values. A first OS run with the original encoding did finish before the DSS failure. It was not saved, so no numbers from it are reported. It showed the same pattern of small, CI-spanning gains.
2. **Added after seeing first results:** the log O/E spatial variant, the 10 CV repeats, the permutation null, the detectable-effect calculation, and the composition/size confound check. They are robustness and sensitivity analyses and do not redefine the primary model. They were added because the first DSS result looked too good to trust.
3. **DSS** was agreed as a secondary endpoint before modelling.

## Part II: pre-registered external replication in METABRIC

### P2.1 Design

The analysis plan, [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md), was committed to git before any METABRIC outcome was linked to features. Three amendments followed, each committed before outcomes were linked and each giving its reason:
1. the cell-class mapping;
2. a correction to how the calibration slope is defined;
3. an ER-negative column rule and the CELESTA-Lite labels.

Before the real run, every code path was dry-run on data with **outcomes shuffled across patients**. In that run all models sat near C = 0.5, so the bugs it surfaced were fixed without any real result being visible.

- **Data.** METABRIC IMC single cells (Zenodo 10.5281/zenodo.5850952; one 383 MB member range-extracted from a 6.2 GB archive) and public cBioPortal clinical data (`brca_metabric`).
- **Cells.** Invasive-tumour images only; imaging-artefact cells removed.
- **Cohort.** Patients with at least 500 cells and some epithelium: **N = 492** (543 images, 904,202 cells).
- **Endpoints.** OS and breast-cancer-specific survival (DSS), administratively censored at 15 years, the horizon fixed in the plan. That leaves 231 OS events and 156 DSS events.
- **Clinical variables.** The same 7 as Basel. pN is derived from the positive-node count using the AJCC bands.
- **Cell classes.** The published phenotypes are mapped to Basel's 6 coarse classes; the table is in Amendment 1.
- **Models and protocol.** Identical to part I.

### P2.2 Within-METABRIC replication (A1) and the confirmatory test (H1)

Overall survival:

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.704 | 0.670 to 0.735 | 0.711, 0.695, 0.705 |
| clinical+composition | 0.697 | 0.662 to 0.730 | 0.700, 0.688, 0.702 |
| +spatial (z) | 0.697 | 0.662 to 0.729 | 0.701, 0.688, 0.701 |
| +spatial (log O/E) | 0.698 | 0.663 to 0.730 | 0.703, 0.692, 0.700 |
| +spatial (residualised) | 0.699 | 0.663 to 0.731 | 0.703, 0.691, 0.703 |
| +spatial (log O/E, CELESTA labels) | 0.697 | 0.662 to 0.729 | 0.698, 0.690, 0.702 |
| random baseline | 0.497 | 0.475 to 0.520 | 0.495, 0.487, 0.510 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | -0.007 | -0.020 to +0.007 | 0.16 |
| +spatial (z) minus clinical+composition | -0.000 | -0.003 to +0.003 | 0.49 |
| +spatial (log O/E) minus clinical+composition | +0.001 | -0.002 to +0.005 | 0.79 |
| +spatial (residualised) minus clinical+composition | +0.002 | -0.001 to +0.005 | 0.90 |
| +spatial (log O/E, CELESTA labels) minus clinical+composition | +0.000 | -0.003 to +0.003 | 0.49 |

Breast-cancer-specific survival:

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.708 | 0.669 to 0.746 | 0.700, 0.716, 0.709 |
| clinical+composition | 0.713 | 0.672 to 0.750 | 0.702, 0.718, 0.718 |
| +spatial (z) | 0.710 | 0.669 to 0.747 | 0.701, 0.715, 0.713 |
| +spatial (log O/E) | 0.709 | 0.667 to 0.747 | 0.702, 0.709, 0.716 |
| +spatial (residualised) | 0.708 | 0.667 to 0.746 | 0.702, 0.707, 0.716 |
| +spatial (log O/E, CELESTA labels) | 0.707 | 0.665 to 0.746 | 0.701, 0.711, 0.709 |
| random baseline | 0.495 | 0.469 to 0.522 | 0.491, 0.497, 0.497 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.004 | -0.010 to +0.020 | 0.74 |
| +spatial (z) minus clinical+composition | -0.003 | -0.008 to +0.002 | 0.14 |
| +spatial (log O/E) minus clinical+composition | -0.004 | -0.009 to +0.001 | 0.06 |
| +spatial (residualised) minus clinical+composition | -0.004 | -0.011 to +0.002 | 0.09 |
| +spatial (log O/E, CELESTA labels) minus clinical+composition | -0.006 | -0.011 to -0.001 | 0.01 |

**H1, the only confirmatory test in part II.** Under the pre-registered rule, H1 is supported only if (a) the CI lies above 0, (b) the permutation p ≤ 0.05, and (c) the mean over 10 CV repeats is > 0. It is refuted if ΔC ≤ 0.

| ΔC (log O/E minus composition, DSS) | 95% CI | permutation p | mean ΔC over 10 CV repeats | (a) CI > 0 | (b) p ≤ 0.05 | (c) repeat mean > 0 | verdict |
|---|---|---|---|---|---|---|---|
| -0.004 | -0.009 to +0.001 | 0.810 | -0.001 | no | no | no | **refuted** |

The Basel DSS lead (part I, §5.2) did not replicate. Every spatial variant, including the residualised and CELESTA-label versions, has its 10-repeat deltas inside the shuffled-feature null (headline figure, right). The CELESTA-label model's DSS bootstrap CI excludes 0 *below* zero; adding those features slightly hurt.

<details>
<summary>CV-repeat and permutation-null details for METABRIC (click to expand)</summary>

OS:

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | -0.000 | +0.001 ± 0.003 [-0.004, +0.006] | 4/10 | +0.001 ± 0.001 (+0.004) | 0.857 |
| +spatial (log O/E) | +0.001 | +0.001 ± 0.003 [-0.004, +0.004] | 5/10 | +0.001 ± 0.001 (+0.003) | 0.286 |
| +spatial (residualised) | +0.002 | +0.001 ± 0.003 [-0.007, +0.005] | 6/10 | +0.001 ± 0.001 (+0.002) | 0.190 |
| +spatial (log O/E, CELESTA labels) | +0.000 | +0.000 ± 0.002 [-0.002, +0.004] | 5/10 | +0.001 ± 0.001 (+0.004) | 0.810 |

DSS:

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | -0.003 | -0.002 ± 0.005 [-0.008, +0.011] | 2/10 | -0.002 ± 0.002 (+0.003) | 0.762 |
| +spatial (log O/E) | -0.004 | -0.001 ± 0.007 [-0.009, +0.016] | 3/10 | -0.002 ± 0.003 (+0.003) | 0.810 |
| +spatial (residualised) | -0.004 | +0.000 ± 0.006 [-0.011, +0.009] | 5/10 | -0.002 ± 0.002 (+0.001) | 0.857 |
| +spatial (log O/E, CELESTA labels) | -0.006 | -0.001 ± 0.006 [-0.010, +0.008] | 5/10 | -0.002 ± 0.002 (+0.001) | 0.857 |

</details>

### P2.3 Transfer: trained on Basel, tested once on METABRIC (A2)

These models use only features that exist in both cohorts, so composition is the 6 coarse classes. Scalers and imputers were fit on Basel only.

OS:

| Model (trained on Basel) | external C-index on METABRIC | 95% CI |
|---|---|---|
| clinical | 0.708 | 0.674 to 0.739 |
| clinical+composition | 0.683 | 0.649 to 0.717 |
| +spatial (z) | 0.684 | 0.651 to 0.717 |
| +spatial (log O/E) | 0.684 | 0.650 to 0.717 |

| Transfer comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | -0.025 | -0.044 to -0.006 | 0.00 |
| +spatial (z) minus clinical+composition | +0.001 | -0.002 to +0.005 | 0.73 |
| +spatial (log O/E) minus clinical+composition | +0.000 | -0.006 to +0.007 | 0.55 |

DSS:

| Model (trained on Basel) | external C-index on METABRIC | 95% CI |
|---|---|---|
| clinical | 0.702 | 0.661 to 0.738 |
| clinical+composition | 0.689 | 0.646 to 0.728 |
| +spatial (z) | 0.686 | 0.644 to 0.724 |
| +spatial (log O/E) | 0.688 | 0.646 to 0.727 |

| Transfer comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | -0.013 | -0.023 to -0.003 | 0.01 |
| +spatial (z) minus clinical+composition | -0.003 | -0.008 to +0.002 | 0.14 |
| +spatial (log O/E) minus clinical+composition | -0.001 | -0.014 to +0.012 | 0.44 |

- **The clinical model transfers intact.** Its external C equals what METABRIC achieves internally.
- **Basel composition does not transfer.** Adding it costs 0.025 (OS) and 0.013 (DSS), with CIs excluding 0. The likely reasons are that the coarse classes come from different panels and phenotyping pipelines in the two studies, and that composition coefficients learned on 280 patients are partly noise.
- **Spatial features neither help nor hurt transfer.**

### P2.4 Calibration and time-dependent AUC (C)

The calibration slope regresses the outcome on the out-of-fold linear predictor (1 is ideal; < 1 means overfit, too-extreme predictions). Time-dependent AUC is cumulative/dynamic, with IPCW.

METABRIC, OS:

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.37 (1.08 to 2.17) | 0.76 (0.71 to 0.81) | 0.76 (0.71 to 0.80) |
| clinical+composition | 0.76 (0.62 to 0.98) | 0.75 (0.70 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (z) | 0.77 (0.62 to 0.99) | 0.75 (0.69 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (log O/E) | 0.77 (0.63 to 0.98) | 0.75 (0.70 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (residualised) | 0.77 (0.63 to 1.00) | 0.75 (0.70 to 0.80) | 0.75 (0.71 to 0.80) |
| +spatial (log O/E, CELESTA labels) | 0.76 (0.62 to 0.99) | 0.75 (0.70 to 0.79) | 0.75 (0.70 to 0.80) |
| random baseline | -0.04 (-0.11 to 0.04) | 0.51 (0.48 to 0.54) | 0.49 (0.45 to 0.52) |

METABRIC, DSS:

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.64 (1.38 to 2.02) | 0.78 (0.72 to 0.83) | 0.73 (0.67 to 0.78) |
| clinical+composition | 0.85 (0.65 to 1.08) | 0.78 (0.73 to 0.83) | 0.74 (0.69 to 0.79) |
| +spatial (z) | 0.81 (0.62 to 1.02) | 0.78 (0.73 to 0.82) | 0.74 (0.68 to 0.79) |
| +spatial (log O/E) | 0.82 (0.63 to 1.04) | 0.78 (0.73 to 0.83) | 0.74 (0.68 to 0.78) |
| +spatial (residualised) | 0.81 (0.63 to 1.03) | 0.78 (0.72 to 0.82) | 0.74 (0.68 to 0.78) |
| +spatial (log O/E, CELESTA labels) | 0.81 (0.62 to 1.05) | 0.77 (0.72 to 0.82) | 0.73 (0.68 to 0.78) |
| random baseline | -0.03 (-0.13 to 0.06) | 0.50 (0.46 to 0.53) | 0.48 (0.44 to 0.52) |

- **Time-dependent AUC** tells the same story as Harrell's C: no layer improves on clinical.
- **Calibration** separates the models.
  - The omics models (slope ≈ 0.8) are mildly overfit: the elastic net keeps clinical columns unpenalised and adds noisy omics terms.
  - The clinical ridge model is *under*-confident (slope 1.37–1.64). Its inner CV often picks a very large ridge penalty, sometimes at the edge of the grid (part I, §9). That shrinks the linear predictor without changing its ranking, which is why discrimination is unaffected.
  - Basel shows the same pattern (generated tables: `results/part2/tables.md`, section C).
- **Practical note:** a discrimination-only comparison would have missed both problems.

### P2.5 Composition-conditioned ("residualised") spatial features (B) and the ER-negative subgroup (D)

- **Residualised spatial features.** Each spatial feature is replaced by its residual after regression on CLR composition and log cell count, fit in-fold (`ResidualizeSpatial`, unit-tested to remove a planted composition effect and keep a planted arrangement effect).
  - In METABRIC (pre-registered secondary): no gain (tables above).
  - In Basel (exploratory, since those data had already been seen): the DSS gain is +0.022 on the 3 main repeats but +0.001 ± 0.022 across 10 repeats. This repeats the part I pattern of a favourable partition, not a signal.
- **ER-negative subgroup** (exploratory, gated on ≥ 40 DSS events; 106 patients, 48 events). ER and PR were dropped as near-constant (Amendment 3a).
  - For DSS, z-score spatial features gain +0.036 on the 3 main repeats and +0.027 ± 0.026 across 10 repeats (9/10 positive), with permutation p = 0.095.
  - This is the only setting in the study where spatial features lean consistently positive, and immune organisation is biologically most plausible in ER-negative disease.
  - It is exploratory, one of many subgroup comparisons, and does not reach the permutation threshold, so it is reported as a lead and not a finding. Full tables are in `results/part2/tables.md`.

## Part III: CELESTA-Lite, spatially informed cell typing

[CELESTA](https://doi.org/10.1038/s41592-022-01498-z) (Zhang et al., *Nature Methods* 2022) assigns cell types from marker expression and, for ambiguous cells, from the types of their spatial neighbours (a Potts-model Markov random field). `src/celesta_lite/` is a vectorised Python port that follows the authors' reference R code (plevritis-lab/CELESTA):
- **Marker probabilities:** an equal-weight two-component mixture per marker.
- **Anchor cells:** cells assigned from markers alone (high markers EP ≥ 0.7, low markers EP ≤ 0.9).
- **Index cells:** the rest, assigned by a mean-field update u_ik ∝ exp(F_ik) · exp(β_ik Σ_{j∈kNN(i), s_j=k} p_jk), with β_ik = 5(1 − d_ik/h).
- **Signature updating** toward the cells already assigned.
- **Convergence** when fewer than 1% of cells change per iteration.

It is independent of, and not endorsed by, the original authors.

**Troubleshooting that mattered:**
1. **Anchor threshold.** A first version, built from the paper text, found zero anchor cells: it required a cell-type probability ≥ 0.5. The reference code's default threshold is 0.
2. **Marker-probability scale.** The reference maps expression to EP with a slope-1 sigmoid of (x − crossing point). That suits CODEX intensities, which span several arcsinh units. On METABRIC's dim IMC channels the sigmoid is nearly flat: CD3's two mixture components sit only 0.15 arcsinh units apart (`results/celesta/v1/marker_models.csv`). The port adds a scale-free option, the mixture posterior P(high | x) made monotone around the crossing point, and uses it here. The reference method remains the library default. Both were run on all cells with the v1 signature:

| EP method | artifact-filtered | cells assigned | agreement (assigned) | κ (assigned) | CELESTA T cells | T recall | mean CD3 EP, published T | mean CD3 EP, published Tumor |
|---|---|---|---|---|---|---|---|---|
| reference (sigmoid of x − crossing point) | 24.3% | 62.8% | 0.876 | 0.628 | 575 | 0.00 | 0.16 | 0.09 |
| posterior (this port) | 4.9% | 92.2% | 0.673 | 0.476 | 79,788 | 0.33 | 0.71 | 0.47 |

   With the reference EP, a quarter of all cells fall to the artifact filter (every marker EP below 0.4, or every one above 0.9), and the method finds almost no T cells. Its higher agreement comes from abstaining: it labels only the cells whose panCK or CD20 signal is unambiguous. The posterior EP restores a usable CD3 contrast (0.71 vs. 0.47), at the price of lower agreement on a much larger set of cells. Reproduce with the `v1_reference_ep` run.

   ![Reference vs posterior EP on one core](figures/celesta_ep_methods.png)

   *The same core as the example below. With the reference EP, the published T cells (orange) disappear and most stromal and immune cells end up Unknown. With the posterior EP they are typed again.*

![CELESTA-Lite example core](figures/celesta_example_core_v1.png)

**Applied to all 1,066,966 METABRIC tumour cells** (cohort-level marker mixtures, then the MRF per image). It runs in about 2 minutes on a laptop. Agreement is measured against the published phenotypes, which are themselves clustering output, so this is concordance and not accuracy:

| Signature | Labelling | cells assigned | agreement (assigned cells) | Cohen's κ (assigned cells) |
|---|---|---|---|---|
| v1 | markers only (argmax score) | 95.1% | 0.602 | 0.406 |
| v1 | CELESTA-Lite (markers + spatial MRF) | 92.2% | 0.673 | 0.476 |
| v2 | markers only (argmax score) | 96.2% | 0.600 | 0.406 |
| v2 | CELESTA-Lite (markers + spatial MRF) | 93.6% | 0.669 | 0.475 |

| Signature | non-anchor cells assigned by the MRF | agreement: CELESTA-Lite | agreement: markers only, same cells |
|---|---|---|---|
| v1 | 334,187 | 0.589 | 0.398 |
| v2 | 333,914 | 0.587 | 0.395 |

| Class (v1) | published cells | CELESTA-Lite cells | recall | precision |
|---|---|---|---|---|
| Tumor | 601,652 | 554,518 | 0.83 | 0.90 |
| T | 78,977 | 79,788 | 0.33 | 0.32 |
| B | 39,584 | 62,605 | 0.55 | 0.35 |
| Macrophage | 50,837 | 127,942 | 0.57 | 0.23 |
| Endothelial | 48,508 | 59,322 | 0.38 | 0.31 |
| Stroma | 247,408 | 99,652 | 0.28 | 0.71 |

- **The spatial step does what CELESTA claims.** On the cells markers alone cannot settle, it raises agreement from 0.398 to 0.589.
- **Overall concordance is moderate.** The signature uses only 7 markers, and two of them separate cell types weakly here: CD3 is dim, and CD68 is non-specific. T cells and macrophages are where the methods disagree most. v2, which adds the two fibroblast markers (FSP1, Podoplanin) that v1 omitted, changed almost nothing, so the signature was not tuned further against the published labels.
- **Downstream effect.** Spatial features computed from CELESTA-Lite labels correlate with the published-label features at ρ = 0.82–0.84 for tumor-centred summaries but only 0.18–0.44 for T-cell pairs (`results/part2/tables.md`). The METABRIC model built on them (pre-registered in Amendment 3) adds nothing (§P2.2). Better-informed cell labels did not uncover prognostic spatial signal.

Reproduce: `uv run python scripts/run_celesta_metabric.py`, which writes `results/celesta/{v1,v2,v1_reference_ep}/`.

## Part IV: tumour genotype and spatial organisation

Parts I–III asked whether spatial features predict survival. Part IV asks a different, descriptive question and makes no prognostic claim:
- **Q1:** Do breast cancer subtypes (PAM50, plus claudin-low) differ in how their cells are arranged, beyond how many cells of each type they contain?
- **Q2:** Is TP53 (and, secondarily, PIK3CA) mutation associated with a distinct arrangement within subtype?
- **Q3:** Do METABRIC's subtype differences replicate in Basel with receptor-based surrogate subtypes?

**Short answer.**
- **Subtype and arrangement (Q1):** after adjusting for composition and image size, no subtype–arrangement association survived FDR control (0 of 243 tests, against 7 unadjusted). The subtype differences seen in co-location here are differences in composition.
- **Hypotheses:** H1, H2 and Zhang et al.'s node-status question (d) all failed.
- **Replication:** one arrangement replicated across cohorts. Macrophage–endothelial co-localisation is higher in basal-like (METABRIC) and triple-negative (Basel) tumours.
- **Positive controls:** only partly reproduced. The nulls mean "not detected with this pipeline and sample", not "absent".

### IV.1 Approach

**Cell typing.** This part applies CELESTA (Zhang et al., *Nature Methods* 2022) as the cell-typing engine. It uses the authors' R package unmodified (plevritis-lab/CELESTA @ `e5e92ff`), not the Python port of Part III.
- **Markers.** The signature uses the lineage markers present in both panels: panCK, CD3, CD20, CD68, CD31-vWF and SMA.
- **Classes.** Round 1 gives Tumour, T, B, macrophage, endothelial and stroma. In METABRIC only, round 2 splits T cells into CD8 T, CD4 T and Treg; Basel has no CD4, CD8 or FOXP3.
- **Why the input is rescaled.** The package's expression probability is a slope-1 sigmoid in the input's units, and its `cofactor` argument is ignored in the code (it always uses 10). IMC input is therefore scaled before the package is called: arcsinh(x/1), each channel divided by its cohort 99th percentile × 10.
- **How the scaling was chosen.** On 20 seeded images per cohort, by a rule fixed in advance that uses only assignment rate, filtered fraction and class presence.
- **Thresholds.** Package defaults.

**Co-location.** The Plevritis lab's own CLQ code was used (plevritis-lab/Spatial_Permutation_and_Normalization @ `02b6da5`, the "colocatome" code), with its functions evaluated unmodified.
- **Formula.** CLQ_{b→a} = (C_{b→a}/N_a)/(N_b/(N−1)), where C counts type-b cells among the k nearest neighbours of type-a cells.
- **Scale.** The lab's implementation counts neighbours rather than fractions, so its CLQ is k × the textbook value. We report CLQ/k, for which 1 means random.
- **Neighbourhoods.** The lab code uses k nearest neighbours only, so k = 20 (the lab default), 10 and 40 replace the planned 20/40/80 µm radii.
- **Permutations.** 500 within-image label permutations at k = 20 give per-image significance and confirm the null is centred on 1.

**Models.** Patient-level linear models:

    CLQ ~ group + CLR composition + log typed cells + age + grade

- Benjamini–Hochberg is applied within each question, across all ordered pairs × k.
- Hypotheses, numeric predictions and pass rules were committed before any label was linked (`HYPOTHESES_GENOTYPE.md`, commit `e6e9389`), together with what prior work predicts (`PRIOR_EXPECTATIONS.md`).
- Every choice that the task or the packages left open is in `AMBIGUITIES_GENOTYPE.md`, including a table of which code (original or new) ran each step.

**Data.**
- **METABRIC:** 749 invasive-tumour images and 693 patients.
  - 588 patients link to cBioPortal subtype or sequencing (METABRIC_173 panel).
  - PAM50: LumA 208, LumB 152, claudin-low 88, HER2-enriched 55, Basal 42.
- **Basel:** 288 images and 280 patients, with per-cell intensities streamed from the Zenodo archive.
  - Surrogate subtypes: ER⁺/HER2⁻ 179, HER2⁺ 52, TN 47.

### IV.2 Cell typing

| Cohort | Class | Defining marker | Images with marker | Images assigned | Share |
|---|---|---|---|---|---|
| metabric | Tumour | panCK | 709 | 707 | 0.997 |
| metabric | T cell | CD3 | 725 | 696 | 0.960 |
| metabric | B cell | CD20 | 710 | 695 | 0.979 |
| metabric | Macrophage | CD68 | 737 | 737 | 1.000 |
| metabric | Endothelial | CD31-vWF | 746 | 746 | 1.000 |
| metabric | Stroma | SMA | 745 | 744 | 0.999 |
| metabric | CD8 T cell | CD8 | 728 | 718 | 0.986 |
| metabric | CD4 T cell | CD4 | 729 | 712 | 0.977 |
| metabric | Treg | FOXP3 | 742 | 720 | 0.970 |
| basel | Tumour | panCK | 277 | 277 | 1.000 |
| basel | T cell | CD3 | 287 | 287 | 1.000 |
| basel | B cell | CD20 | 286 | 286 | 1.000 |
| basel | Macrophage | CD68 | 285 | 285 | 1.000 |
| basel | Endothelial | CD31-vWF | 287 | 287 | 1.000 |
| basel | Stroma | SMA | 286 | 286 | 1.000 |
- **Acceptance:** every class is assigned in ≥ 96% of the images where its defining marker is present (at least 5 cells above the package's own GMM crossing point). The acceptance bar was 80%.
- **Concordance with the published labels:** reported as concordance only, without judging either labelling. At the round-1 level it is 0.73 in both cohorts (κ 0.55 and 0.57).
- **The METABRIC T-cell subclasses agree much less:**

| CELESTA class | Published class | CELESTA cells | Published cells | Share of CELESTA in published | Share of published recovered | CELESTA cells published as epithelial |
|---|---|---|---|---|---|---|
| CD8 T cell | CD8+ T cells | 33,005 | 34,588 | 0.492 | 0.469 | 8,393 |
| CD4 T cell | CD4+ T cells | 27,581 | 15,647 | 0.149 | 0.263 | 8,412 |
| Treg | TReg & TEx | 25,432 | 6,318 | 0.020 | 0.082 | 12,352 |
The CELESTA Treg class is the weak point: half of its cells are cells the published labels call epithelial. It was not retuned against the published labels. Treg-based results are always shown with the published-label refit beside them.

![Marker expression by CELESTA class](figures/genotype_marker_heatmap.png)

![Example cores per subtype](figures/genotype_example_cores.png)

*One core per PAM50 subtype, chosen as the core closest to the subtype's median composition, coloured by official CELESTA class.*

### IV.3 Hypotheses and positive controls

| ID | Test | Estimate | p | n | Outcome |
|---|---|---|---|---|---|
| H1 | METABRIC Tumour -> CD8 T cell, k = 20, Basal - LumA (adjusted) | -0.047 [-0.157, 0.064] | 0.406 | 368 | fail |
| H1-rev | METABRIC CD8 T cell -> Tumour, k = 20, Basal - LumA (adjusted) | -0.040 [-0.146, 0.066] | 0.458 | 368 | fail |
| H2 | METABRIC ER+, Tumour -> Macrophage, k = 20, TP53 mutant - wild-type (adjusted) | -0.005 [-0.065, 0.054] | 0.857 | 341 | fail |
| H3 | METABRIC round-1 contrasts with q < 0.05 (9) with the same sign in Basel | 9.000 | – | 9 | pass |
| H3-strict | share of replicated signs > 0.5 (one-sided binomial; needs >= 5 discoveries) | 1.000 | 0.002 | 9 | pass |
| P1 | METABRIC immune fraction: Basal - LumA and Her2 - LumA | 0.105 [0.033, 0.178] | 0.005 | 436 | fail |
| P2a | Tumour -> Tumour and Stroma -> Stroma, Basal - LumA > 0 | – | – | – | fail |
| P2b | B -> T and T -> B (round-1), Basal - LumA > 0 | – | – | – | pass |
| P2c | Treg -> CD8 T and Treg -> CD4 T, Basal - LumA > 0 (CELESTA) | – | – | – | fail |
| d | METABRIC Tumour -> Treg, k = 20, node-positive - node-negative (adjusted, CELESTA) | -0.006 [-0.069, 0.056] | 0.846 | 386 | fail |
| P2 | at least 2 of P2a-c pass | 1.000 | – | – | fail |
- **H1 fails.** H1 predicted that Basal tumours show more CD8 T–tumour co-localisation than LumA. The estimate is slightly negative at every k, adjusted or not (forest plot below).
- **H2 fails.** H2 predicted that TP53 mutation raises macrophage–tumour co-localisation in ER⁺ tumours.
- **H3 passes as pre-registered, but read it narrowly.** The 9 METABRIC discoveries are 3 ordered pairs at 3 values of k, so they are not independent. Macrophage↔Endothelial replicates in Basel. Macrophage→Tumour has the same sign but is about 0 in Basel (IV.5).
- **(d) fails with CELESTA labels.** Question (d) asks whether tumour–Treg co-localisation is higher in node-positive tumours. With the published Treg class the estimate points the way Zhang et al. found in head-and-neck cancer, but it is not significant after adjustment, in only 88 patients.

**P1, immune fraction (a composition measure):**

| Cohort | Labels | Contrast | Immune-fraction difference | p | n |
|---|---|---|---|---|---|
| METABRIC | CELESTA | LumB | 0.045 [-0.004, 0.094] | 0.069 | 436 |
| METABRIC | CELESTA | Her2 | 0.040 [-0.027, 0.108] | 0.242 | 436 |
| METABRIC | CELESTA | Basal | 0.105 [0.033, 0.178] | 0.005 | 436 |
| METABRIC | CELESTA | claudin-low | 0.120 [0.064, 0.177] | 3.8e-05 | 436 |
| METABRIC | published | LumB | 0.011 [-0.030, 0.052] | 0.599 | 463 |
| METABRIC | published | Her2 | 0.075 [0.019, 0.131] | 0.009 | 463 |
| METABRIC | published | Basal | 0.088 [0.026, 0.149] | 0.005 | 463 |
| METABRIC | published | claudin-low | 0.110 [0.063, 0.158] | 5.7e-06 | 463 |
| Basel | CELESTA | HER2+ | 0.029 [-0.027, 0.085] | 0.308 | 260 |
| Basel | CELESTA | TN | 0.114 [0.054, 0.175] | 2.3e-04 | 260 |
- **Reproduced:** Basal > LumA in METABRIC and TN > ER⁺/HER2⁻ in Basel, as Ali 2020, Denkert 2018 and Jackson 2020 report.
- **Not reproduced with CELESTA labels:** HER2-enriched > LumA (+0.040, CI includes 0). With the published labels it is reproduced (+0.075).
- **Verdict:** P1 fails by the pre-registered rule, and the failure comes from the CELESTA typing.

**P2, subtype-specific interactions reported by Ali 2020 and Danenberg 2022:**

| Control | Pair | Basal - LumA | p | n |
|---|---|---|---|---|
| P2a | Tumour -> Tumour | -0.260 [-0.879, 0.359] | 0.410 | 430 |
| P2a | Stroma -> Stroma | 0.329 [-0.268, 0.925] | 0.279 | 433 |
| P2b | B cell -> T cell | 0.122 [-0.313, 0.558] | 0.581 | 346 |
| P2b | T cell -> B cell | 0.059 [-0.395, 0.513] | 0.798 | 346 |
| P2c | Treg -> CD8 T cell | -0.036 [-0.482, 0.409] | 0.873 | 353 |
| P2c | Treg -> CD4 T cell | -0.051 [-0.603, 0.500] | 0.854 | 346 |
| P2c (published labels) | Treg -> CD8 T cell | 0.438 [-1.606, 2.483] | 0.670 |  85 |
| P2c (published labels) | Treg -> CD4 T cell | -0.352 [-0.934, 0.230] | 0.231 |  78 |
- **P2a, Ali's "separation phenotype":** not reproduced.
- **P2c, Danenberg's suppressed expansion:** not reproduced.
- **P2b, TLS-like B–T co-localisation:** in the predicted direction only, with wide intervals.
- **P2 fails overall.**
- **The translation is approximate.** These prior results come from touching-cell permutation tests (Ali) and TME-structure enrichment (Danenberg) on 32–57 clustered phenotypes. Mapping them onto six- to nine-class kNN CLQs loses detail.
- **P3 (BRCA1)** is descriptive only. Only 4–6 BRCA1-mutant patients have a defined Treg pair (`results/genotype/tables.md`).

![Effect sizes for H1–H3](figures/genotype_forest.png)

### IV.4 Arrangement beyond amount

| Question | Adjusted | Tests | q < 0.05 | p < 0.05 |
|---|---|---|---|---|
| Q1 subtype (F-test) | composition-adjusted | 243 | 0 | 24 |
| Q1 subtype (F-test) | unadjusted | 243 | 7 | 42 |
| Q2 TP53 | composition-adjusted | 243 | 0 | 33 |
| Q2 TP53 | unadjusted | 243 | 0 | 39 |
| Q2 PIK3CA | composition-adjusted | 243 | 0 | 10 |
| Q2 PIK3CA | unadjusted | 243 | 0 | 27 |
- **Subtype:** the 7 unadjusted subtype discoveries are B, macrophage and endothelial cells around tumour cells. None survives adjustment for CLR composition. The closest is Macrophage→Tumour, with adjusted p = 0.003 and q = 0.11.
- **Mutations:** no TP53 or PIK3CA association survives FDR control.
- **Strongest TP53 lead (exploratory, not FDR-significant):** fewer tumour cells around CD8 T cells in TP53-mutant tumours, consistent at all three k:

| Pair | k | TP53 mut - wt | p | q | n | n mutant |
|---|---|---|---|---|---|---|
| CD8 T cell -> Tumour | 40 | -0.091 [-0.150, -0.032] | 0.003 | 0.267 | 380 | 142 |
| CD8 T cell -> Tumour | 20 | -0.093 [-0.155, -0.031] | 0.003 | 0.267 | 380 | 142 |
| CD8 T cell -> Tumour | 10 | -0.091 [-0.155, -0.028] | 0.005 | 0.267 | 380 | 142 |
| Macrophage -> Tumour | 20 | -0.066 [-0.113, -0.020] | 0.005 | 0.267 | 439 | 155 |
| Macrophage -> Tumour | 40 | -0.061 [-0.104, -0.018] | 0.005 | 0.267 | 439 | 155 |
| Endothelial -> Tumour | 40 | -0.062 [-0.108, -0.016] | 0.009 | 0.309 | 441 | 156 |
| Macrophage -> Endothelial | 40 | 0.125 [0.031, 0.220] | 0.010 | 0.309 | 443 | 159 |
| Tumour -> CD8 T cell | 20 | -0.085 [-0.150, -0.020] | 0.010 | 0.309 | 380 | 142 |
![Subtype-by-pair CLQ heatmap](figures/genotype_heatmap.png)

*Mean patient CLQ (log2) by PAM50 subtype for every ordered pair at k = 20. No pair carries a significance mark: none reached q < 0.05 after adjustment.*

**Robustness:**

| Check | Value |
|---|---|
| Subtype-label permutation: share of 1,000 permutations with any BH discovery | 0.039 |
| Subtype-label permutation: mean null p-value | 0.503 |
| Permutation check calibrated (pre-registered rule) | yes |
| Q1 discoveries significant (p < 0.05) at both other k | no discoveries |
| Published labels, k = 20: Spearman correlation of subtype F across pairs | 0.184 |
| Published labels: Basal - LumA sign agreement (pairs defined in both) | 0.765 |
| Published labels: CELESTA k = 20 discoveries with published q < 0.05 | 0 of 0 |
| Lab convention (Unknown cells kept in graph), k = 20: Spearman correlation of F | 0.741 |
| Lab convention: CELESTA k = 20 discoveries also q < 0.05 | 0 of 0 |
- **The permutation check is calibrated.** Under 1,000 subtype-label shuffles, 3.9% of runs produce any discovery, and the null p-values average 0.50.
- **The published-label refit agrees only loosely** (F correlation 0.18; Basal − LumA sign agreement 0.77). The co-location results depend on the cell typing.
- **The lab's own convention changes little.** Keeping Unknown cells in the neighbour graph gives an F correlation of 0.74 with the primary run.

### IV.5 Replication in Basel

| METABRIC round-1 contrasts tested | q < 0.05 | Refitted in Basel | Same sign | Same sign and Basel p < 0.05 |
|---|---|---|---|---|
| 216 | 9 | 9 | 9 | 4 |
| Pair | k | METABRIC | METABRIC q | Basel | Basel p | Same sign |
|---|---|---|---|---|---|---|
| Endothelial -> Macrophage | 20 | Basal - LumA 0.428 [0.223, 0.632] | 0.010 | TN 0.178 [0.027, 0.328] | 0.021 | yes |
| Macrophage -> Endothelial | 10 | Basal - LumA 0.417 [0.201, 0.633] | 0.014 | TN 0.106 [-0.048, 0.260] | 0.175 | yes |
| Endothelial -> Macrophage | 10 | Basal - LumA 0.452 [0.207, 0.696] | 0.014 | TN 0.158 [0.008, 0.307] | 0.039 | yes |
| Macrophage -> Endothelial | 20 | Basal - LumA 0.362 [0.169, 0.555] | 0.014 | TN 0.125 [-0.028, 0.278] | 0.110 | yes |
| Endothelial -> Macrophage | 40 | Basal - LumA 0.325 [0.155, 0.494] | 0.014 | TN 0.168 [0.031, 0.306] | 0.017 | yes |
| Macrophage -> Endothelial | 40 | Basal - LumA 0.286 [0.129, 0.442] | 0.014 | TN 0.149 [0.009, 0.288] | 0.037 | yes |
| Macrophage -> Tumour | 10 | Basal - LumA -0.141 [-0.220, -0.062] | 0.015 | TN -0.008 [-0.077, 0.061] | 0.812 | yes |
| Macrophage -> Tumour | 20 | Basal - LumA -0.136 [-0.212, -0.059] | 0.015 | TN -0.008 [-0.081, 0.066] | 0.838 | yes |
| Macrophage -> Tumour | 40 | Basal - LumA -0.116 [-0.187, -0.045] | 0.033 | TN -0.010 [-0.082, 0.063] | 0.795 | yes |
![METABRIC vs Basel replication](figures/genotype_replication.png)

*Every METABRIC Basal − LumA and HER2-enriched − LumA contrast against its Basel TN and HER2⁺ counterpart, using the six classes both cohorts share. Red points are the METABRIC discoveries.*

- **Replicated:** macrophage–endothelial co-localisation is higher in basal-like and triple-negative tumours in both cohorts.
- **Not replicated:** the METABRIC finding of fewer macrophages around tumour cells in basal-like tumours is flat in Basel.

### IV.6 Agreement with prior work

| Prior result | This analysis | Verdict |
|---|---|---|
| More immune cells in basal-like/TN (Ali 2020, Denkert 2018, Jackson 2020) | immune fraction Basal +0.105 (METABRIC), TN +0.114 (Basel) | agrees |
| More immune cells in HER2-enriched (Denkert 2018) | +0.040 with CELESTA labels (p 0.242); +0.075 with published labels | partly agrees (depends on the cell typing) |
| Immune cells excluded from the tumour mass (Jackson 2020) | median Tumour -> T and Tumour -> Macrophage CLQ 0.50-0.88 across subtypes and cohorts | agrees |
| T cells sit near endothelium (Jackson 2020) | median T -> Endothelial CLQ 1.05-1.41 | agrees |
| Homotypic 'separation phenotype' in basal-like (Ali 2020) | Basal - LumA: Tumour -> Tumour -0.260, Stroma -> Stroma +0.329 (both CIs include 0) | not reproduced |
| TLS-like B-T structures in basal-like (Danenberg 2022) | Basal - LumA: B -> T +0.122, T -> B +0.059 (both CIs include 0) | not reproduced (same direction only) |
| Treg-T 'suppressed expansion' in ER-negative/basal and BRCA1 (Danenberg 2022) | Basal - LumA: Treg -> CD8 T -0.036, Treg -> CD4 T -0.051; BRCA1-mutant patients with a defined pair: 4 | not reproduced / BRCA1 not testable |
| Tumour-Treg co-localisation higher in node-positive disease (Zhang 2022, HNSCC) | CELESTA labels -0.006 (p 0.846); published labels +0.126 (p 0.118, n = 88) | not reproduced; no earlier breast-cancer test found |

New beyond these papers:
- Composition explains the subtype differences in co-location seen here (IV.4).
- Macrophage–endothelial co-localisation is higher in basal-like/TN tumours, and this replicates across cohorts (IV.5).
- TP53 within subtype has no FDR-significant arrangement effect.
- To our knowledge, this is the first test of the node-status question (d) in breast cancer, and it is null.

### IV.7 Deviations from the plan

1. **kNN instead of radii.** The lab's CLQ code supports k nearest neighbours only, so k = 10/20/40 replace 20/40/80 µm. "Leave-one-radius-out" became leave-one-k-out.
2. **Permutations at k = 20 only.** The first run, with permutations at all k, exhausted memory and was lost (AMBIGUITIES S2.1). The association models use observed CLQs, which are computed at every k.
3. **Statistics and figures in R, not Python.** This was the user's choice; Python is kept where existing repo code did the job.
4. **No T subsets in Basel.** Its panel has no CD4, CD8 or FOXP3, so replication uses the six shared classes.
5. **Danenberg's public R code was not run.** It works on their own TME-structure assignments, not on our classes. Ali 2020's analysis code is not public.
6. **H3 counts dependent contrasts,** so it is reported as one replicated arrangement (IV.5).

### IV.8 Limitations

- **Small sample of each tumour.** There is one small TMA core per patient (rarely two), so co-location describes only a few hundred µm of each tumour.
- **Targeted mutation panel.** Mutations come from the 173-gene METABRIC panel. BRCA1 (6 mutant patients) and TP53 within basal tumours (2 wild-type) cannot be tested.
- **Cross-sectional design.** The results are associations: nothing here is causal or prognostic.
- **Scaling choice.** The IMC scaling determines the CELESTA expression probabilities. It was fixed by a pre-declared rule, but other reasonable choices would give somewhat different labels.
- **Cell classes.** The shared marker panel is small, and the CELESTA Treg class is poorly concordant with the published labels.

Memo: [docs/GENOTYPE_MEMO.md](docs/GENOTYPE_MEMO.md).

Reproduce:

```bash
Rscript scripts/setup_r.R                                  # pin CELESTA + lab CLQ code; results/genotype/environment.json
uv run python scripts/download_data.py --genomics --basel-expression   # cBioPortal mutations; Basel intensities (~0.6 GB)
uv run python scripts/run_genotype_spatial.py              # cohort tables (Python loaders)
Rscript scripts/celesta_example.R                          # bundled CELESTA example vs. its reference output
Rscript scripts/celesta_scaling.R                          # scaling choice on 20 seeded images per cohort
for c in metabric basel; do Rscript scripts/celesta_official.R --cells data/processed/genotype/${c}_cells.csv.gz \
  --signature results/genotype/signature_${c}.csv --thresholds results/genotype/thresholds_${c}.csv \
  --out results/genotype/cell_labels_${c}.csv.gz --cofactor 1 --norm p99 --workers 18 \
  --status-out results/genotype/celesta_status_${c}.csv; done
Rscript scripts/celesta_report.R
# CLQ: the six clq_official.R runs listed in AMBIGUITIES_GENOTYPE.md S2.2 (about 1.7 h on 18 cores)
Rscript scripts/clq_features.R
Rscript scripts/genotype_models.R                          # (--dry-run: all labels shuffled)
Rscript scripts/genotype_report.R                          # results/genotype/tables.md + figures
```

## 8. Reproduction

```bash
uv sync                                              # Python 3.12, pinned dependencies (pyproject.toml + uv.lock)
uv run python scripts/download_data.py               # Basel: ~31 MB download + 0.16 MB range-extracted; checksums verified
uv run python scripts/run.py                         # part I: results/ and figures/ (a few minutes on a laptop CPU)
uv run python scripts/report_tables.py               # results/tables.md
uv run python scripts/download_data.py --metabric    # METABRIC: 383 MB range-extracted (850 MB on disk) + cBioPortal clinical
uv run python scripts/run_celesta_metabric.py        # part III: results/celesta/ (~2 min per run, 3 runs)
uv run python scripts/run_part2.py                   # part II: results/part2/ (~31 min; needs the CELESTA v1 labels)
uv run python scripts/report_part2.py                # results/part2/tables.md
uv run pytest                                        # 32 tests (the R suite runs inside one of them when R is installed)
```

- **Order matters for part II.** Its CELESTA-label model (Amendment 3) reads `results/celesta/v1/cell_labels.csv.gz`, so run part III first. The labels are committed, so `run_part2.py` also works on its own.
- **METABRIC provenance.** The Zenodo member is CRC-checked against the archive's central directory. The cBioPortal clinical data come from a live API, so the MD5 of the snapshot used here is recorded in `results/part2/data_manifest.json`. A later download that differs from it may shift part II numbers slightly.

- **Listing only.** `scripts/download_data.py --list` lists the record and the 36.8 GB archive's members by reading only its central directory.
- **Configuration and raw records.** Seeds, k, permutations and grids live in `CONFIG` in `scripts/run.py` and are saved to `results/summary.json` with package versions and runtime. Per-fold records (seed, fold, selected α, inner C, test C, non-zero coefficients, test patient IDs) are in `results/{os,dss}_per_fold.csv`. Out-of-fold risks are in `results/{os,dss}_oof_risk.csv`.
- **Determinism.** A run from empty `results/`, `figures/` and spatial-feature cache reproduced `results/tables.md` exactly, apart from the runtime line. Raw risk scores differ between runs only at about 1e-6, which is floating-point noise from multithreaded numerics.

**Tests** (`tests/`):
- **Leakage.** Imputer, scaler and CLR outputs depend only on training rows, and a spy on every fit inside the nested CV confirms no fit sees an outer-test patient.
- **Patient-level splits.** Outer folds partition patients, and multi-image patients collapse to one row before CV.
- **Neighbour enrichment.** An exact hand-worked null on a 4-node path graph, the closed-form permutation mean, and segregated vs. checkerboard layouts.
- **Toy end-to-end.** Simulated images in which T-cell infiltration drives survival at fixed composition; the spatial model must beat the clinical and composition models by more than 0.1 C.
- **Residualisation.** In-fold residualisation removes a planted composition effect, keeps a planted arrangement effect, and is fit on training rows only.
- **METABRIC harmonisation.** The pN banding from node counts, the 180-month censoring, the phenotype-to-class mapping (which reaches every Basel class and fails loudly on an unmapped phenotype), and the cell filters.
- **CELESTA-Lite.** The mixture fit and its crossing point, both EP methods (the posterior is invariant to rescaling the channel), cell-type scores, lineage rounds, the artifact filter, and a toy core where the spatial step must rescue cells whose markers are ambiguous.
- **README integrity.** Every result-table row here appears in `results/tables.md`, `results/part2/tables.md` or `results/genotype/tables.md`.
- **Part IV (R, `tests/R/`, run by `tests/test_r_scripts.py`).**
  - The Plevritis lab's `get_CLQ` on a hand-computed toy layout.
  - Random labelling gives CLQ/k of about 1, and within-image permutations are centred on 1.
  - Rare and absent types are handled.
  - Observed CLQs are identical with and without permutations.
  - The output format of the CELESTA bridge `scripts/celesta_official.R`.
  - Skipped when R is not installed (CI).

CI (`.github/workflows/tests.yml`) runs the tests on every push.

Layout:
- `scripts/`: `download_data`, `run` and `report_tables` (part I); `run_part2` and `report_part2` (part II); `run_celesta_metabric` (part III).
- Part IV scripts:
  - Python: `run_genotype_spatial.py` (cohort build).
  - R: `setup_r`, `celesta_example`, `celesta_official`, `celesta_scaling`, `celesta_report`, `clq_official`, `clq_features`, `genotype_models`, `genotype_report`.
  - Upstream code is pinned into git-ignored `third_party/`.
- `src/spatialsurv/`: `data`, `features`, `models`, `cv`, `eval`, `robustness`, `plots`, `metabric`.
- `src/celesta_lite/`: `preprocessing`, `mrf`, `graph`, `benchmark`.
- Also `tests/`, `results/` (`part2/`, `celesta/`) and `figures/`.

## 9. Caveats and limitations

**Scope of the null result**
- **It concerns incremental value over clinical variables, for this feature panel.** It does not say spatial structure is unrelated to outcome. Danenberg et al. report outcome associations for METABRIC TME structures. Those analyses ask a different question (association, not out-of-fold gain over a clinical model). The panel here is also small: pairwise enrichment, neighbour fractions and mixing on 6 coarse classes. Richer descriptors, such as cellular neighbourhoods, distance-based statistics or the papers' own community features, might behave differently.
- **One small TMA core per patient (rarely two), in both cohorts.** Spatial features describe a tiny sample of each tumor. Rare-class pairs are often undefined (Tumor–B is missing for 108/280 Basel patients and imputed).

**Statistics**
- **Basel is small relative to the feature count.** 280 patients and 79 or 57 events against up to 45 features. Per-model 95% CIs are about ±0.05 wide. METABRIC is larger (§P2.1), which is why its CIs exclude gains above +0.005.
- **Overall survival includes non-cancer deaths.** In Basel that is 22 of 79. DSS censors them, so it estimates a cause-specific hazard, not a cumulative incidence.
- **Harrell's C depends on the censoring distribution.** Calibration and time-dependent AUC are reported in §P2.4. Decision-curve analysis was not done, because no layer cleared the discrimination bar.
- **The clinical ridge α sometimes hits the grid edge** (10³). Large ridge penalties mostly rescale the linear predictor, to which the rank-based C is insensitive. They do affect calibration, which is the under-confidence seen in §P2.4. The grid was not re-tuned after seeing results.
- **The ER-negative lead is one of many exploratory comparisons.** Part I and part II together ran several endpoints, feature variants and a subgroup, and no multiplicity correction is applied to the exploratory ones. Only H1 was confirmatory.

**Cell labels**
- **The cell labels are the published clusterings,** made once on all cells of each cohort. They are outcome-blind, but they are a shared upstream step that a strictly prospective pipeline would not have. Cross-cohort class mapping (Amendment 1) is a judgement call, such as placing CD57⁺ cells with T cells.
- **The published labels are not ground truth.** CELESTA-Lite "agreement" is concordance between two automated methods.
- **CELESTA-Lite departs from the reference in three ways.** It fits marker mixtures once per cohort rather than per sample. It uses a posterior EP rather than the reference sigmoid (Part III). It uses an arcsinh cofactor of 1 for IMC counts. The signature was written from marker biology, and only v1 and v2 were tried. It was not tuned against the published labels.

## 10. Next steps

- **Test the ER-negative lead prospectively.** Run a pre-registered DSS analysis of z-score spatial features in an independent ER-negative or triple-negative cohort with cancer-specific follow-up, with the hypothesis and decision rule fixed in advance, as was done for H1.
- **Richer spatial descriptors under the same protocol.** Cellular neighbourhoods, multi-scale radius graphs and cross-type K/L functions, residualised on composition and image size in-fold as in §P2.5.
- **Competing-risks modelling** (cause-specific Cox for both causes, or Fine–Gray) instead of censoring non-cancer deaths.
- **Calibration of the transferred model.** Transfer (§P2.3) was scored on discrimination only, and the within-cohort clinical models are under-confident (§P2.4). The next step is to measure calibration of the Basel-trained model on METABRIC, and to test whether recalibrating it (baseline hazard and slope) on a small METABRIC subset makes absolute risk transfer too.
- **Validate CELESTA-Lite against an annotated reference.** For example, use manually gated cells, or CODEX data where the reference sigmoid EP is in its intended regime, before using its labels for biology.

## References

- Jackson, H.W. et al. The single-cell pathology landscape of breast cancer. *Nature* 578, 615–620 (2020).
- Keren, L. et al. A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging. *Cell* 174, 1373–1387 (2018).
- Danenberg, E. et al. Breast tumor microenvironment structures are associated with genomic features and clinical outcome. *Nature Genetics* 54, 660–669 (2022).
- Pölsterl, S. scikit-survival: a library for time-to-event analysis built on top of scikit-learn. *JMLR* 21(212), 1–6 (2020).
- Zhang, W. et al. Identification of cell types in multiplexed in situ images by combining protein expression and spatial information using CELESTA. *Nature Methods* 19, 759–769 (2022).
- Ali, H.R. et al. Imaging mass cytometry and multiplatform genomics define the phenogenomic landscape of breast cancer. *Nature Cancer* 1, 163–175 (2020).
- Denkert, C. et al. Tumour-infiltrating lymphocytes and prognosis in different subtypes of breast cancer: a pooled analysis of 3771 patients treated with neoadjuvant therapy. *Lancet Oncology* 19, 40–50 (2018).
- Leslie, T.F. & Kronenfeld, B.J. The colocation quotient: a new measure of spatial association between categorical subsets of points. *Geographical Analysis* 43, 306–326 (2011).
- Plevritis lab code: github.com/plevritis-lab/CELESTA (`e5e92ff`) and github.com/plevritis-lab/Spatial_Permutation_and_Normalization (`02b6da5`).
