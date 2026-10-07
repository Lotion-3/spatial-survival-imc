## Cohort

- cohort: Basel TMA (Jackson et al. 2020), tumor cores only
- n_tumor_patients: 281
- n_dropped_nonpositive_or_missing_time: 1
- n_patients: 280
- n_images: 288
- patients_with_2_images: 8
- n_cells: 762776
- events_os: 79
- events_dss: 57
- median_followup_months_all: 74.0
- missing_clinical: {'pN': 13, 'ER': 1, 'PR': 6}
- spatial_missing: {'enrich_Tumor__Tumor': 1, 'enrich_Tumor__T': 9, 'enrich_Tumor__Macrophage': 11, 'enrich_Tumor__B': 108, 'enrich_Tumor__Stroma': 1, 'enrich_Tumor__Endothelial': 24, 'enrich_T__Macrophage': 17, 'enrich_T__B': 107, 'enrich_Stroma__T': 8, 'logoe_Tumor__Tumor': 1, 'logoe_Tumor__T': 9, 'logoe_Tumor__Macrophage': 11, 'logoe_Tumor__B': 108, 'logoe_Tumor__Stroma': 1, 'logoe_Tumor__Endothelial': 24, 'logoe_T__Macrophage': 17, 'logoe_T__B': 107, 'logoe_Stroma__T': 8}
- n_composition_features: 27
- n_spatial_features: 11
- grade1_patients: 38
- grade1_dss_events: 2
- pM1_patients: 7

## Overall survival (primary)

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

## Breast-cancer-specific survival (sensitivity)

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

## Detectable effect size

| Endpoint | Comparison | SD of paired bootstrap ΔC | detectable ΔC (80% power) | with CV-repeat variability |
|---|---|---|---|---|
| OS | clinical+composition minus clinical | 0.0120 | 0.034 | n/a |
| OS | clinical+composition+spatial minus clinical+composition | 0.0035 | 0.010 | 0.014 |
| DSS | clinical+composition minus clinical | 0.0200 | 0.056 | n/a |
| DSS | clinical+composition+spatial minus clinical+composition | 0.0079 | 0.022 | 0.032 |

## Penalty and sparsity

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

## Spatial features vs. image size and immune content

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

## Kaplan-Meier

- best model: clinical+composition+spatial, log-rank p = 6.5e-09
- clinical-only model: log-rank p = 8.3e-10
- patients assigned to the same risk half by both: 94%

runtime_seconds: 383.1
