# Analysis plan, part II (pre-registered)

Written and committed **before** any METABRIC outcome data were linked to features. At the time of writing, I had inspected only the METABRIC IMC file listing, the header of `SingleCells.csv`, and the list of cBioPortal clinical attribute names. The commit timestamp of this file is the pre-registration record. Any change after that commit goes in the *Amendments* section at the end, with its reason and its date.

Part I (the Basel analysis, README §1–10) is finished and is not modified by this plan. Its main open lead was a weak, post-hoc DSS signal for size-robust (log O/E) spatial enrichment, so part II tests that lead in an independent cohort.

## A. External cohort: METABRIC IMC (Danenberg et al. 2022)

**Data.**
- Imaging: Zenodo 10.5281/zenodo.5850952 (CC-BY-4.0). Only `SingleCells.csv` is used (one row per cell: `metabric_id`, `ImageNumber`, `cellPhenotype`, `is_epithelial`, `Location_Center_X/Y`). It is extracted by HTTP range request from the 6.2 GB archive; the images are not downloaded.
- Clinical: the public cBioPortal study `brca_metabric`, patient attributes.

**Cohort.** All patients with ≥ 1 IMC image, OS time > 0 and a known vital status. Patients are excluded if they have fewer than 500 cells in total or no epithelial cells, so that composition and spatial features are defined.

**Endpoints.**
- **OS:** `OS_MONTHS`; event = `VITAL_STATUS` ∈ {"Died of Disease", "Died of Other Causes"}.
- **DSS:** event = "Died of Disease"; other deaths are censored. Patients with unknown cause of death are excluded from DSS only.
- **Administrative censoring at 180 months for both endpoints.** METABRIC follow-up exceeds 25 years, whereas Basel's maximum is 233 months and its median 74. Very late deaths in this older cohort are dominated by non-cancer causes, and a common horizon makes the cohorts comparable. This choice is fixed here, before outcomes are seen.

**Clinical variables (harmonised to Basel's 7).**

| Basel | METABRIC (cBioPortal) | Mapping |
|---|---|---|
| age | `AGE_AT_DIAGNOSIS` | as is |
| tumor_size (mm) | `TUMOR_SIZE` | as is |
| grade (1–3) | `GRADE` | as is |
| pN (0–3) | `LYMPH_NODES_EXAMINED_POSITIVE` | 0 → 0, 1–3 → 1, 4–9 → 2, ≥ 10 → 3 (AJCC pN count bands) |
| ER, PR, HER2 | `ER_STATUS`, `PR_STATUS`, `HER2_STATUS` | Positive → 1, Negative → 0 |

Missing values are median-imputed inside training folds, exactly as in part I.

**Cell classes.** METABRIC `cellPhenotype` labels are mapped to the 6 coarse classes used in Basel (Tumor, T, B, Macrophage, Endothelial, Stroma):
- Cells with `is_epithelial == 1` → Tumor.
- Non-epithelial phenotypes are assigned by name: T cells (CD4⁺/CD8⁺/Treg/T-cell names), B cells, macrophages/myeloid (CD68⁺/CD163⁺/macrophage names), endothelial, and everything else stromal (fibroblasts, myofibroblasts, pericytes, unassigned non-epithelial).

The exact phenotype → class table is written in Amendment 1 **after viewing the phenotype names only** (no outcomes) and before any model is run. Ambiguous mixed-lineage phenotypes, such as "T cells & APCs", go to the lineage named first, and the table lists every such case.

**Features.**
- Clinical (7).
- Composition, two versions:
  - **coarse:** the 6 coarse classes, CLR-transformed. This version is identical in both cohorts and is used for transfer.
  - **native:** METABRIC's own phenotypes, CLR-transformed. Used for the within-METABRIC analysis only, mirroring Basel's 27 metaclusters.
- Spatial: the same 11 features with the same code and parameters (kNN k = 8, 200 permutations, minimum 5 cells per class), plus the log O/E variant. Image → patient aggregation is weighted by cell count.
- Coordinates are in pixels (1 µm for IMC).

### A1. Within-METABRIC replication of the 3-layer comparison (primary for part II)

The same protocol as part I:
- Ridge Cox for clinical; elastic-net Cox (l1_ratio 0.5) with clinical columns unpenalised for the others.
- 20-point α grids and inner 5-fold CV.
- Outer 5-fold CV × 3 seeds, pooled out-of-fold Harrell's C.
- 1,000 shared patient bootstrap resamples.
- 10 CV repeats and a 20-permutation spatial-block null.

Models (composition = native):
1. clinical
2. clinical + composition
3. clinical + composition + spatial (z)
4. clinical + composition + spatial (log O/E)

**Confirmatory hypothesis H1 (the Basel lead).** For **DSS**, model 4 has higher C than model 2. H1 is **supported** only if all three of the following hold:
- (a) the paired-bootstrap 95% CI for ΔC lies above 0;
- (b) the permutation p ≤ 0.05 (spatial block shuffled across patients, 20 permutations, 3 seeds);
- (c) the mean ΔC over 10 CV repeats is > 0.

It is **refuted** if the point estimate of ΔC is ≤ 0. Anything else is **inconclusive**. This is the only confirmatory test in part II; everything else is descriptive or exploratory.

**Primary descriptive question.** The same three-layer comparison for **OS** (models 1–3), reported with the same tables as part I.

### A2. Transfer: train on Basel, test on METABRIC

Models are fit on all 280 Basel patients, using the same inner-CV α selection, and evaluated once on METABRIC. Features must exist in both cohorts, so composition = coarse (6 classes):
1. clinical
2. clinical + coarse composition
3. clinical + coarse composition + spatial (z)
4. clinical + coarse composition + spatial (log O/E)

Scalers and imputers are fit on Basel only (strict transfer). Reported: external Harrell's C with patient-bootstrap 95% CIs and paired differences, for OS and DSS. The question is whether the layers' ranking from Basel transfers. **Expectation, stated in advance:** spatial z-scores depend on image size, and METABRIC image sizes may differ, so z-score transfer may degrade; the log O/E variant should transfer better. The protocol has a single run with no tuning on METABRIC outcomes.

## B. Composition-conditioned spatial features

Part I found that `tumor_immune_mixing` and `tumor_immune_nbr_frac` are ~90% explained by immune fraction, and that z-scores track image size. To isolate *arrangement beyond amount*, each spatial feature is replaced by its residual from a linear regression on [CLR composition, log(cell count)]. The regression is **fit on the training fold only** and applied to the test fold; it is part of the pipeline and covered by the leakage test.

- Model B: clinical + composition + residualised spatial (log O/E pairs + the two ratio features).
- Basel: **exploratory**, since part I's data have been seen.
- METABRIC: **pre-registered secondary**, for both OS and DSS; ΔC vs. model 2 with the same reporting as A1. No confirmatory claim.

## C. Calibration and time-dependent discrimination

For every model in part I and A1, on out-of-fold predictions:
- **Calibration slope.** Univariable Cox regression of the outcome on the out-of-fold linear predictor (standardised within each outer fold, so fold scales match); 1 is ideal, < 1 means overfitting. Reported with a bootstrap 95% CI.
- **Time-dependent AUC** at 60 and 120 months (cumulative/dynamic, IPCW; scikit-survival `cumulative_dynamic_auc`).

These are descriptive. They answer whether the clinical model's lead holds under other metrics and whether the omics models are miscalibrated by overfitting.

## D. Exploratory subgroup: ER-negative tumours (METABRIC)

Immune organisation is thought to matter most in ER-negative disease. Within METABRIC ER-negative patients, I run the A1 OS and DSS comparisons (models 2 vs. 3 and 4) with the same protocol, **if the number of DSS events is ≥ 40**; otherwise this section is skipped and that is reported. This is labelled exploratory, with no multiplicity correction and no claim.

## Multiplicity and reporting rules

- One confirmatory test (H1). All other p-values and CIs are descriptive.
- Every analysis listed here is reported whatever its result, including skipped analyses and the reason.
- README tables are generated from saved results by `scripts/report_tables.py` and checked by `tests/test_readme_numbers.py`, as in part I.
- Runtime: METABRIC is larger, so part II runs in its own script (`scripts/run_metabric.py`) and its runtime is reported.

## Amendments

### Amendment 1: phenotype → class table and cell filters (made before any model was run)

This amendment is based on the phenotype labels, the `is_epithelial` flag, and the column descriptions in `SingleCellsAnnotation.xlsx`. No outcome data were linked.

- **Images.** Only cells with `is_tumour == 1` (the 693 main invasive tumours) are used. Cells with `is_hotAggregate == 1` (imaging artefacts) are removed. Pixel size is 1 µm, per the annotation file.
- **All epithelial phenotypes → Tumor:** Basal, CD15⁺, CK8-18⁺ ER^hi, CK8-18^hi CXCL12^hi, CK8-18^hi ER^lo, CK⁺ CXCL12⁺, CK^lo ER^lo, CK^lo ER^med, CK^med ER^lo, ER^hi CXCL12⁺, Ep CD57⁺, Ep Ki67⁺, HER2⁺, MHC I & II^hi, MHC I^hi CD57⁺, MHC^hi CD15⁺.
- **Non-epithelial phenotypes:**

| METABRIC phenotype | Class | Note |
|---|---|---|
| CD4⁺ T cells, CD8⁺ T cells, T_Reg & T_Ex | T | |
| CD4⁺ T cells & APCs | T | mixed lineage; first-named lineage rule |
| CD57⁺ (non-epithelial) | T | CD57 marks NK cells and terminally differentiated T cells. Basel has no NK class and its CD57⁺ lymphocytes fall in T-cell clusters, so the closest Basel class is T. This is a deviation from the "everything else → Stroma" default, declared here. |
| B cells | B | |
| CD38⁺ lymphocytes | B | plasma cells/plasmablasts, B lineage |
| Macrophages, Macrophages & granulocytes, Granulocytes | Macrophage | myeloid |
| Endothelial | Endothelial | |
| Fibroblasts, Fibroblasts FSP1⁺, Myofibroblasts, Myofibroblasts PDPN⁺ | Stroma | |
| Ki67⁺ (non-epithelial, lineage unassigned) | Stroma | "unassigned non-epithelial" default |

The table was drafted from the first 500,000 rows. The loader raises an error if any non-epithelial phenotype in the full file is unmapped; in that case a further amendment will be added.

### Amendment 2: correction to section C (made before any model was run)

Section C said the linear predictor would be "standardised within each outer fold" before computing the calibration slope. That was a mistake. Standardising removes the scale of the linear predictor, which is exactly what the calibration slope measures, so the slope would reflect only discrimination. The calibration slope is computed on the **raw** out-of-fold linear predictor. Each fold's model sees features standardised on its own training data, so training-fold linear predictors are centred near 0 and pooling across folds is reasonable.

### Amendment 3: ER-negative clinical columns and a CELESTA-label analysis (made before any real METABRIC outcome was linked to features)

Status when this was written: the part II pipeline had only been dry-run on a truncated copy of the METABRIC cells with **outcomes shuffled across patients**, to find bugs. No real METABRIC outcome had been linked to any feature. The CELESTA-Lite cell typing (commit "Add CELESTA-Lite …") uses no outcome data.

**(a) Section D, clinical columns in the ER-negative subgroup.** ER is constant (0) by definition in this subgroup, and other binary clinical variables can be nearly constant. The outcome-shuffled dry run showed that an unpenalised clinical block then cannot be fitted, the same problem that removed pM in part I. Rule: within the subgroup, drop ER and any binary clinical column whose minority level has fewer than 10 patients. The rule uses covariates only, never outcomes. The dropped columns are reported.

**(b) New section E: spatial features from CELESTA-Lite labels (METABRIC; secondary, descriptive).** Part I and section A compute spatial features from published cell labels. Section E asks whether labels from an independent, spatially informed cell typer change the answer.

- **Labels.** CELESTA-Lite as committed: signature **v1** (the version written before any comparison with the published labels), posterior expression probabilities, arcsinh cofactor 1, cohort-level marker mixtures, and reference defaults for everything else (k = 5, bandwidth 100 µm, γ = 5). Signature v2 is not used here. Cells labelled Unknown are left out of the spatial graph.
- **Features.** The same 11 spatial features, log O/E version (as in H1), computed from CELESTA coarse classes with the same code and parameters as section A.
- **Model E.** clinical + native composition + spatial (log O/E, CELESTA labels), compared with model 2 (clinical + native composition), for OS and DSS. Reporting is the same as A1: bootstrap ΔC, 10 CV repeats, and a 20-permutation null. There is no confirmatory claim.
