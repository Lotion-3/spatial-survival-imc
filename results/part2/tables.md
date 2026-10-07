## METABRIC cohort

- cohort: METABRIC IMC (Danenberg et al. 2022), invasive tumour images, artefact cells removed
- patients_with_images_after_cell_filters: 577
- patients_dropped_cells_lt_min_or_no_epithelium: 116
- patients_with_clinical: 492
- n_patients: 492
- n_images: 543
- n_cells: 904202
- events_os_180m: 231
- events_dss_180m: 156
- median_time_months: 102.85
- er_negative_patients: 106
- er_negative_dss_events: 48
- missing_clinical: {'tumor_size': 1, 'grade': 7, 'pN': 14}
- n_native_phenotypes: 32
- median_cells_per_patient: 1512.0

## A1 METABRIC OS

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

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | -0.000 | +0.001 ± 0.003 [-0.004, +0.006] | 4/10 | +0.001 ± 0.001 (+0.004) | 0.857 |
| +spatial (log O/E) | +0.001 | +0.001 ± 0.003 [-0.004, +0.004] | 5/10 | +0.001 ± 0.001 (+0.003) | 0.286 |
| +spatial (residualised) | +0.002 | +0.001 ± 0.003 [-0.007, +0.005] | 6/10 | +0.001 ± 0.001 (+0.002) | 0.190 |
| +spatial (log O/E, CELESTA labels) | +0.000 | +0.000 ± 0.002 [-0.002, +0.004] | 5/10 | +0.001 ± 0.001 (+0.004) | 0.810 |

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.37 (1.08 to 2.17) | 0.76 (0.71 to 0.81) | 0.76 (0.71 to 0.80) |
| clinical+composition | 0.76 (0.62 to 0.98) | 0.75 (0.70 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (z) | 0.77 (0.62 to 0.99) | 0.75 (0.69 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (log O/E) | 0.77 (0.63 to 0.98) | 0.75 (0.70 to 0.80) | 0.75 (0.70 to 0.80) |
| +spatial (residualised) | 0.77 (0.63 to 1.00) | 0.75 (0.70 to 0.80) | 0.75 (0.71 to 0.80) |
| +spatial (log O/E, CELESTA labels) | 0.76 (0.62 to 0.99) | 0.75 (0.70 to 0.79) | 0.75 (0.70 to 0.80) |
| random baseline | -0.04 (-0.11 to 0.04) | 0.51 (0.48 to 0.54) | 0.49 (0.45 to 0.52) |

## A1 METABRIC DSS

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

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | -0.003 | -0.002 ± 0.005 [-0.008, +0.011] | 2/10 | -0.002 ± 0.002 (+0.003) | 0.762 |
| +spatial (log O/E) | -0.004 | -0.001 ± 0.007 [-0.009, +0.016] | 3/10 | -0.002 ± 0.003 (+0.003) | 0.810 |
| +spatial (residualised) | -0.004 | +0.000 ± 0.006 [-0.011, +0.009] | 5/10 | -0.002 ± 0.002 (+0.001) | 0.857 |
| +spatial (log O/E, CELESTA labels) | -0.006 | -0.001 ± 0.006 [-0.010, +0.008] | 5/10 | -0.002 ± 0.002 (+0.001) | 0.857 |

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.64 (1.38 to 2.02) | 0.78 (0.72 to 0.83) | 0.73 (0.67 to 0.78) |
| clinical+composition | 0.85 (0.65 to 1.08) | 0.78 (0.73 to 0.83) | 0.74 (0.69 to 0.79) |
| +spatial (z) | 0.81 (0.62 to 1.02) | 0.78 (0.73 to 0.82) | 0.74 (0.68 to 0.79) |
| +spatial (log O/E) | 0.82 (0.63 to 1.04) | 0.78 (0.73 to 0.83) | 0.74 (0.68 to 0.78) |
| +spatial (residualised) | 0.81 (0.63 to 1.03) | 0.78 (0.72 to 0.82) | 0.74 (0.68 to 0.78) |
| +spatial (log O/E, CELESTA labels) | 0.81 (0.62 to 1.05) | 0.77 (0.72 to 0.82) | 0.73 (0.68 to 0.78) |
| random baseline | -0.03 (-0.13 to 0.06) | 0.50 (0.46 to 0.53) | 0.48 (0.44 to 0.52) |

## H1 (confirmatory)

| ΔC (log O/E minus composition, DSS) | 95% CI | permutation p | mean ΔC over 10 CV repeats | (a) CI > 0 | (b) p ≤ 0.05 | (c) repeat mean > 0 | verdict |
|---|---|---|---|---|---|---|---|
| -0.004 | -0.009 to +0.001 | 0.810 | -0.001 | no | no | no | **refuted** |

## A2 transfer OS

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

## A2 transfer DSS

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

## B Basel residualised OS (exploratory)

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.741 | 0.687 to 0.792 | 0.735, 0.755, 0.733 |
| clinical+composition | 0.748 | 0.690 to 0.801 | 0.740, 0.760, 0.743 |
| +spatial (residualised) | 0.746 | 0.689 to 0.798 | 0.740, 0.748, 0.749 |
| random baseline | 0.503 | 0.465 to 0.540 | 0.564, 0.484, 0.461 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.007 | -0.017 to +0.029 | 0.67 |
| +spatial (residualised) minus clinical+composition | -0.002 | -0.011 to +0.007 | 0.31 |

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (residualised) | -0.002 | +0.000 ± 0.010 [-0.017, +0.014] | 6/10 | +0.002 ± 0.002 (+0.005) | 1.000 |

## B Basel residualised DSS (exploratory)

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.721 | 0.660 to 0.777 | 0.734, 0.737, 0.693 |
| clinical+composition | 0.723 | 0.660 to 0.780 | 0.753, 0.707, 0.710 |
| +spatial (residualised) | 0.746 | 0.687 to 0.795 | 0.762, 0.755, 0.720 |
| random baseline | 0.505 | 0.459 to 0.553 | 0.561, 0.464, 0.489 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.002 | -0.038 to +0.040 | 0.54 |
| +spatial (residualised) minus clinical+composition | +0.022 | -0.000 to +0.048 | 0.97 |

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (residualised) | +0.022 | +0.001 ± 0.022 [-0.037, +0.048] | 5/10 | +0.003 ± 0.008 (+0.019) | 0.048 |

## C Basel calibration OS

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.30 (0.95 to 1.71) | 0.79 (0.72 to 0.85) | 0.74 (0.65 to 0.82) |
| clinical+composition | 0.79 (0.57 to 1.04) | 0.80 (0.73 to 0.86) | 0.74 (0.64 to 0.82) |
| clinical+composition+spatial | 0.80 (0.57 to 1.05) | 0.80 (0.73 to 0.87) | 0.75 (0.65 to 0.83) |
| clinical+composition+spatial (log O/E, sensitivity) | 0.75 (0.54 to 0.99) | 0.79 (0.72 to 0.86) | 0.74 (0.65 to 0.83) |
| random baseline | 0.01 (-0.10 to 0.13) | 0.50 (0.45 to 0.55) | 0.54 (0.48 to 0.61) |

## C Basel calibration DSS

| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |
|---|---|---|---|
| clinical | 1.27 (0.81 to 1.81) | 0.77 (0.70 to 0.83) | 0.74 (0.65 to 0.83) |
| clinical+composition | 0.65 (0.44 to 0.89) | 0.76 (0.68 to 0.84) | 0.70 (0.60 to 0.79) |
| clinical+composition+spatial | 0.71 (0.49 to 0.95) | 0.78 (0.71 to 0.85) | 0.73 (0.63 to 0.82) |
| clinical+composition+spatial (log O/E, sensitivity) | 0.73 (0.53 to 0.95) | 0.79 (0.72 to 0.85) | 0.76 (0.67 to 0.85) |
| random baseline | 0.01 (-0.14 to 0.15) | 0.51 (0.45 to 0.56) | 0.54 (0.47 to 0.60) |

## D ER-negative subgroup

- n: 106
- dss_events: 48
- run: True
- clinical_columns_dropped: ['ER', 'PR']

### ER-negative OS

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.621 | 0.561 to 0.681 | 0.617, 0.629, 0.617 |
| clinical+composition | 0.637 | 0.575 to 0.697 | 0.655, 0.632, 0.623 |
| +spatial (z) | 0.620 | 0.561 to 0.674 | 0.647, 0.581, 0.632 |
| +spatial (log O/E) | 0.622 | 0.559 to 0.683 | 0.623, 0.633, 0.610 |
| random baseline | 0.483 | 0.441 to 0.527 | 0.484, 0.457, 0.508 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | +0.015 | -0.029 to +0.061 | 0.74 |
| +spatial (z) minus clinical+composition | -0.017 | -0.055 to +0.025 | 0.23 |
| +spatial (log O/E) minus clinical+composition | -0.015 | -0.043 to +0.012 | 0.17 |

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | -0.017 | -0.011 ± 0.030 [-0.067, +0.036] | 4/10 | -0.008 ± 0.012 (+0.010) | 0.810 |
| +spatial (log O/E) | -0.015 | -0.011 ± 0.024 [-0.053, +0.020] | 4/10 | -0.010 ± 0.010 (+0.008) | 0.714 |

### ER-negative DSS

| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |
|---|---|---|---|
| clinical | 0.631 | 0.550 to 0.707 | 0.624, 0.645, 0.623 |
| clinical+composition | 0.626 | 0.552 to 0.695 | 0.640, 0.607, 0.630 |
| +spatial (z) | 0.662 | 0.587 to 0.731 | 0.647, 0.659, 0.681 |
| +spatial (log O/E) | 0.638 | 0.570 to 0.700 | 0.601, 0.638, 0.673 |
| random baseline | 0.494 | 0.452 to 0.537 | 0.492, 0.471, 0.518 |

| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |
|---|---|---|---|
| clinical+composition minus clinical | -0.005 | -0.073 to +0.065 | 0.45 |
| +spatial (z) minus clinical+composition | +0.036 | -0.013 to +0.086 | 0.93 |
| +spatial (log O/E) minus clinical+composition | +0.012 | -0.034 to +0.054 | 0.67 |

| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |
|---|---|---|---|---|---|
| +spatial (z) | +0.036 | +0.027 ± 0.026 [-0.022, +0.061] | 9/10 | -0.015 ± 0.026 (+0.064) | 0.095 |
| +spatial (log O/E) | +0.012 | +0.026 ± 0.028 [-0.039, +0.059] | 9/10 | -0.013 ± 0.025 (+0.034) | 0.190 |

## CELESTA-Lite vs published labels

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

| Spatial feature (v1 labels) | images | Spearman ρ, CELESTA-Lite vs. published labels |
|---|---|---|
| `enrich_Tumor__Tumor` | 684 | +0.84 |
| `enrich_Tumor__T` | 528 | +0.62 |
| `enrich_Tumor__Macrophage` | 578 | +0.66 |
| `enrich_Tumor__B` | 286 | +0.67 |
| `enrich_Tumor__Stroma` | 645 | +0.68 |
| `enrich_Tumor__Endothelial` | 639 | +0.46 |
| `enrich_T__Macrophage` | 524 | +0.18 |
| `enrich_T__B` | 301 | +0.35 |
| `enrich_Stroma__T` | 539 | +0.35 |
| `tumor_immune_nbr_frac` | 690 | +0.60 |
| `tumor_immune_mixing` | 725 | +0.82 |
| `logoe_Tumor__Tumor` | 684 | +0.82 |
| `logoe_Tumor__T` | 528 | +0.45 |
| `logoe_Tumor__Macrophage` | 578 | +0.62 |
| `logoe_Tumor__B` | 286 | +0.50 |
| `logoe_Tumor__Stroma` | 645 | +0.61 |
| `logoe_Tumor__Endothelial` | 639 | +0.52 |
| `logoe_T__Macrophage` | 524 | +0.23 |
| `logoe_T__B` | 301 | +0.44 |
| `logoe_Stroma__T` | 539 | +0.26 |

runtime_seconds: 1852.2
