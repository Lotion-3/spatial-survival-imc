# Part IV: ambiguities, decisions and code provenance

This file records every choice in Part IV ("tumour genotype and spatial organisation") that the task description, the packages or the data left open: parameters, output slots, data mismatches and label mappings. Entries are added stage by stage, and each is dated by the stage at which it was decided. Nothing here was chosen after seeing subtype or genotype associations with spatial features.

## Pipeline: which code does each step

| Step | Code (original / new) | Source and version | Language | Why |
|---|---|---|---|---|
| Basel download (locations, labels, metadata) | original repo code, unchanged | `scripts/download_data.py` (this repo) | Python | existing downloader; Zenodo MD5s verified |
| Basel per-cell marker intensities | new function in the existing downloader (`--basel-expression`) | Zenodo 10.5281/zenodo.4607374, `SC_dat.csv` member, CRC-checked | Python | CELESTA needs marker intensities, which the repo never downloaded; extends the existing range-request code |
| METABRIC IMC cells | original repo code, unchanged | `download_data.py --metabric`; files on disk match their recorded MD5s, so they were not re-fetched | Python | existing |
| METABRIC mutations, sequenced samples, gene panel | new function in the existing downloader (`--genomics`) | cBioPortal REST API, study `brca_metabric` (snapshot date in `results/genotype/data_manifest.json`) | Python | extends the existing cBioPortal client |
| Cohort build (cells, patients, receptor and PAM50 labels, mutation calls) | new script that calls the repo's loaders unchanged (`spatialsurv.metabric.load_metabric_cells`, `load_metabric_clinical`; `spatialsurv.data.load_cells`, `load_patients`) | `scripts/run_genotype_spatial.py` | Python | the loaders already implement the cohort filters used in Parts I–III |
| R environment and upstream code | new setup script | `scripts/setup_r.R`; pins in `results/genotype/environment.json` | R | installs and pins the original packages and repositories |
| CELESTA bundled example | original package code; new wrapper | CELESTA `e5e92ff` README "Usage" block, demo data `reg009` | R | stage-0 check; compared against the reference output shipped in the package's `tests/testthat` |
| Cell typing | **original package, unmodified** (`CreateCelestaObject`, `FilterCells`, `AssignCells`, `PlotExpProb`, `PlotCellsAnyCombination`); new CLI wrapper `scripts/celesta_official.R` | github.com/plevritis-lab/CELESTA @ `e5e92ffcef5ccada8053c2e042f34c446b6b3565` (the paper's github.com/plevritis/CELESTA resolves to the same repository and HEAD) | R | task requirement; the repo's Python port (`src/celesta_lite/`) is not used for any Part IV result |
| Input scaling (arcsinh, per-channel normalisation) | new (in `celesta_official.R`) | — | R | the package ignores its own `cofactor` argument (see S0.3), so scaling has to happen before the package is called |
| CLQ and per-image permutation test | **original lab functions, unmodified** (`get_CLQ`, `find_cell_type_neighbors`, `KNN_neighbors`, `significance_matrix_gen`), evaluated from the original file; new driver `scripts/clq_official.R` | github.com/plevritis-lab/Spatial_Permutation_and_Normalization @ `02b6da5e7bb8ebabb68a05daa78ae46c36ec0378` ("permutation analysis and normalization strategy for the colocatome paper") | R | the lab's own CLQ code exists, so no re-implementation was written. Its driver loops hard-code 18 cell types and k = 20 and glob files in the working directory, so only the function definitions are used (S0.6) |
| Two-sided permutation p-value per image | new (in `clq_official.R`) | — | R | the lab's rule gives a percentile, not a p-value |
| Tests | new | `tests/R/` (testthat), run by `tests/test_r_scripts.py` | R / Python | CLQ hand example, random-label null, rare-type handling, bridge output format |

## Stage 0: setup

**S0.1 Repositories and versions.**
- CELESTA is installed with `remotes::install_github("plevritis-lab/CELESTA@e5e92ff…")`. Package version 0.0.0.9000; last commit 2023-01-20.
- `git ls-remote https://github.com/plevritis/CELESTA.git HEAD` returns the same SHA, so the URL cited in the paper and the lab URL are the same repository.
- The CLQ repository is pinned at `02b6da5` (2024-10-09).
- Both are cloned into `third_party/`, which is git-ignored, by `scripts/setup_r.R`. That script checks that each working tree is clean (unmodified) and at the pinned SHA.

**S0.2 CELESTA bundled example.**
- There is no vignette. The README "Usage" block (`CreateCelestaObject` → `FilterCells(0.9, 0.4)` → `AssignCells` with the demo threshold vectors) was run on the package's demo data: 2,434 cells, 1.2 s.
- The final labels agree with the package's own reference output (`tests/testthat/project_title_final_cell_type_assignment.csv`) for **100% of cells**.
- The reference output contains 2 cells with final type "15" instead of a name. That is a quirk of `GetFinalInferredCellTypes` which we reproduce rather than fix.
- `PlotExpProb` writes one file per marker named after the marker. The demo marker names contain ":", which Windows does not allow in file names (ggsave wrote NTFS alternate data streams), so the example renames ":" to "_" in a copy used only for plotting. Our own marker names contain no ":".

**S0.3 Parameters verified against the CELESTA source (`R/CELESTA_functions.R` @ e5e92ff).**
- **`CreateCelestaObject(project_title, prior_marker_info, imaging_data_file, cofactor = 10, transform_type = 1, number_of_neighbors = 5, bandwidth = 100)`.**
  - **The `cofactor` argument is ignored.** The function calls `GetMarkerExpMatrix(..., cofactor = 10, ...)` with a literal 10.
  - We therefore arcsinh-transform and normalise in `celesta_official.R` and call the package with `transform_type = 0`. The marker columns of `prior_marker_info` (column 3 onward) must match imaging-data column names exactly.
  - Images with ≤ 10 cells return `NULL`; `celesta_official.R` reports these as status `too_few_cells`.
- **Expression probabilities (`CalcMarkerActivationProbability`).**
  - A 2-component, equal-proportion Gaussian mixture (`Rmixmod`, model `Gaussian_p_Lk_Ck`) is fitted per marker **on the cells of the object**. We build one object per image, so the fit is per image.
  - Cells with zero expression are partly removed before fitting when 10–90% of values are 0.
  - EP = logistic(1 · (x − x_root)), min–max rescaled to [0, 1] within the object, where x_root is the crossing point of the two weighted Gaussians.
  - The slope is fixed at 1 in the input's units, so **the input scale sets how sharp the EPs are**. This is why the scaling choice (Stage 1) matters.
- **`FitGmmModel` calls `set.seed(1)` on every marker.** This resets the global RNG. The CLQ permutations therefore run in a separate process (`clq_official.R`) with their own logged seeds.
- **`FilterCells(high_marker_threshold = 0.9, low_marker_threshold = 0.4)`.**
  - A cell is removed (all rounds set to NA) if **every** marker EP is < 0.4 or **every** marker EP is > 0.9.
  - The README calls the step optional; we run it with the defaults. Removed cells are reported as `Filtered`.
- **`AssignCells(max_iteration = 10, cell_change_threshold = 0.01, min_diff = 0, min_probability = 0, high_expression_threshold_anchor = rep(0.7, 50), low_expression_threshold_anchor = rep(0.9, 50), high_expression_threshold_index = rep(0.5, 50), low_expression_threshold_index = rep(1, 50), scale_factor = 5, bandwidth = 100, save_result = TRUE)`.**
  - Threshold vectors are indexed by **cell-type number** (the third field of `Lineage_level`).
  - `save_result = TRUE` writes a CSV into the working directory, so we pass `FALSE` and read the slot instead.
- **Lineage format `round_parent_number`.**
  - `FindCellsToCheck` takes `unique(Previous_cell_type)` of a round, so **each round can subdivide only one parent type**. The hierarchy has to be a chain: one parent per round.
- **Output slots.**
  - `final_cell_type_assignment` is a character matrix with columns `Round 1..R`, `Cell type number` and `Final cell type`.
  - Unassigned cells are "Unknown" with number 0. Filtered cells are NA.
  - A cell typed at round r but left unassigned at round r + 1 keeps its round-r label. For example, a T cell not resolved into CD8/CD4/Treg stays "T cell".
  - `anchor_cell_type_assignment` holds the anchor (marker-only) assignments per round.

**S0.4 Data mismatches and label mappings.**
- **METABRIC cells.** We use `load_metabric_cells` unchanged: invasive-tumour images (`is_tumour == 1`), hot-aggregate artefacts removed. That gives 749 images, 693 patients and 1,066,966 cells. The marker columns are read with the same filter, and the coordinates are asserted to be identical.
- **METABRIC IMC → cBioPortal.** `metabric_id` equals the cBioPortal `patientId` (e.g. MB-0000). In `brca_metabric` the `sampleId` equals the `patientId` for every sequenced sample, which is checked.
  - 588 of the 693 IMC patients link to a subtype or to sequencing (acceptance ≥ 400).
- **PAM50.**
  - `CLAUDIN_SUBTYPE` values LumA, LumB, Her2, Basal and claudin-low are the five study subtypes.
  - "Normal" (42 IMC patients) and "NC" (not classified), plus patients with no value (106), are excluded from the subtype analyses but kept for the mutation analyses.
  - The cBioPortal label "Her2" is the PAM50 HER2-enriched subtype.
- **Mutations.**
  - Panel: METABRIC_173 (173 genes) for all 2,433 sequenced samples. It contains TP53, PIK3CA, BRCA1, BRCA2 and CASP8, so **P3 (BRCA1) is run**.
  - A gene counts as mutated if the sample has a mutation of type Missense, Nonsense, Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site, Splice_Region, Translation_Start_Site or Nonstop. The API returned no Silent or UTR calls.
  - Wild-type = sequenced and no such mutation. Not sequenced = missing (125 IMC patients).
  - Splice_Region is included because cBioPortal counts it as a mutation. It is 123 of 13,183 calls cohort-wide.
- **Counts among IMC patients.**
  - TP53: 192 mutant / 376 wild-type.
  - PIK3CA: 240 / 328.
  - BRCA1: **6 / 562**. P3 has almost no power; this is reported, not hidden.
  - CASP8: 10 / 558.
  - TP53 by subtype (mutant / wild-type): Basal 38/2, HER2-enriched 38/17, LumA 28/177, LumB 42/103, claudin-low 35/46. **A TP53 effect within Basal tumours is not estimable (2 wild-type)**; within-subtype TP53 contrasts are possible in LumA, LumB, claudin-low and HER2-enriched.
- **Node status (METABRIC).** `LYMPH_NODES_EXAMINED_POSITIVE` > 0, through the existing `pn_from_count` banding (pN > 0).
- **Basel surrogate subtypes** (task definition: ER+/HER2−, HER2+, TN), from `ERStatus`, `PRStatus` and `HER2Status`:
  - HER2+ = HER2 positive, whatever the ER status.
  - ER+/HER2− = ER positive and HER2 negative.
  - TN = ER, PR and HER2 all negative.
  - The single ER−/PR+/HER2− patient fits none of these and is left missing.
  - The published `clinical_type` column (HR+HER2−, HR−HER2+, …) is kept for reference only.
- **Basel patients.** `spatialsurv.data.load_patients` is used unchanged. It drops patients without a positive OS time (a survival-analysis filter from Part I). The number dropped is recorded in `results/genotype/cohort_counts.json`. It is kept for consistency with Parts I–III; this study does not use outcomes.

**S0.5 Marker panels and the signature (decided before any cell typing beyond the one-image smoke test).**
- **Lineage markers present in both panels:** panCK, CK8-18, CK5, CD3, CD20, CD68, CD31-vWF and SMA. HER2 (3B5), Ki-67 and Histone H3 are also shared but are not lineage markers.
- **No pan-immune marker is shared.** Basel has CD45; METABRIC has only CD45RA/CD45RO. So there can be no "Immune" parent in a shared round 1.
- **No fibroblast marker beyond SMA is shared.** Basel has Vimentin and Fibronectin; METABRIC has PDGFRB, FSP1 and Podoplanin. SMA-negative fibroblasts can therefore only be "Unknown" in the shared signature. This is a limitation.
- **Shared round 1 (both cohorts):** Tumour, T cell, B cell, Macrophage, Endothelial, Stroma. The reasoning per row:
  - **Tumour: panCK 1; CD3, CD20, CD68, CD31-vWF 0; SMA NA.** Pan-cytokeratin marks epithelial cells; within invasive-tumour cores we call epithelial cells Tumour, as Parts II–III do. CK8-18 and CK5 are not used: luminal (CK8-18) and basal (CK5) tumours differ in them, and they would split tumour cells by subtype — the variable under study. SMA is NA because myoepithelial cells can be SMA-positive.
  - **T cell: CD3 1; panCK, CD20, CD68 0.** CD3 is the pan-T-cell marker. CD31-vWF and SMA are NA, being irrelevant to T cells.
  - **B cell: CD20 1; panCK, CD3 0.** CD68 is NA because of spill-over in dense lymphoid aggregates.
  - **Macrophage: CD68 1; panCK, CD3, CD20 0.**
  - **Endothelial: CD31-vWF 1; panCK, CD3, CD68 0.** CD20 is NA. SMA is NA because pericytes and vessel walls can be SMA-positive.
  - **Stroma: SMA 1; panCK, CD3, CD20, CD68, CD31-vWF 0.** This captures SMA-positive (myo)fibroblasts.
- **METABRIC round 2 (parent T cell, METABRIC only):**
  - **CD8 T: CD8 1, CD4 0, FOXP3 0.**
  - **CD4 T: CD4 1, CD8 0, FOXP3 0.**
  - **Treg: FOXP3 1, CD8 0, CD4 NA.** IMC CD4 is dim, and FOXP3 defines Tregs.
- Basel lacks CD8, CD4 and FOXP3, so cross-cohort comparisons use only the round-1 classes.
- **Thresholds:** the package defaults for every type (anchor high 0.7 / low 0.9; index high 0.5 / low 1). They are not tuned.

**S0.6 How the lab's CLQ code is run.**
- **What is used unmodified.** `Spatial_Perm_and_Normalization.R` is a script, not a package. Its top-level loops glob `*cell_type_assignment.csv` in the working directory, and its drivers hard-code 18 cell types (`rep(0, 18)`, `cell_types[1:18]`, `1:324`) and k = 20. Sourcing it as is would run those loops. `clq_official.R` therefore parses the original file and evaluates only the definitions of `get_CLQ`, `find_cell_type_neighbors`, `KNN_neighbors` and `significance_matrix_gen`, byte for byte. The driver logic is repeated in our script: label permutation by `sample()` of the whole label vector with coordinates fixed, 500 iterations, and zeros for absent types.
- **Formula as implemented.** `get_CLQ(cell_A = a, cell_B = b)` returns (C_ab / N_a) / (N_b / (N − 1)), with (N_b − 1) when a = b.
  - C_ab = the number of type-b cells among the k nearest neighbours (spdep `knearneigh`, excluding self) of all type-a cells.
  - N counts cells with a type (label not 0 or NA).
  - This is the CELESTA paper's CLQ_{b→a}.
- **C_ab / N_a is a count out of k, not a fraction.** So the lab's CLQ equals **k × the Leslie–Kronenfeld CLQ**, with expected value ≈ k under random labelling. The lab's own normalisation step confirms this: it warns if the mean CLQ differs from k by more than 3.
  - We report both `clq` (lab scale) and `clq_k = clq / k` (classic scale, null ≈ 1).
  - Models use `clq_k`, so effects are on the familiar scale; this is a linear rescaling, so tests are unchanged.
- **Rare types.** `get_CLQ` returns 0 when either type has ≤ 5 cells. Those pairs are set to NA, as the lab's `CLQ_normalization_by_sample` does. A literal 0 would otherwise be read as complete avoidance.
- **Neighbourhood.** The lab code supports only k nearest neighbours (no radius), so the task's radii (20/40/80 µm) do not apply. We use **k = 20 (lab default) plus k = 10 and k = 40**, and "leave-one-radius-out" becomes leave-one-k-out.
- **Unknown and filtered cells.**
  - In the lab's clinical example input, 14% of cells are Unknown (type 0) and stay in the kNN graph. They are neighbours of typed cells but are not counted in N, which deflates every CLQ by roughly the image's typed fraction. That would tie CLQ to each image's assignment rate.
  - **Primary analysis:** Unknown and Filtered cells are removed before the kNN graph is built (`--keep-unknown 0`), so the classic-scale null is 1.
  - **Sensitivity analysis:** the lab's convention (`--keep-unknown 1`) at k = 20.
- **Significance per image.**
  - The lab rule (`significance_matrix_gen`, run unmodified in a scratch directory) gives the observed CLQ's percentile among its 500 permutations, flagged if ≤ 0.05 or ≥ 0.95.
  - We add a two-sided permutation p-value. Per-image significance is descriptive; the association tests are patient-level.
  - **Quirk:** the lab rule computes the percentile as the first rank of the observed value among ties. A pair whose observed and permuted CLQs are all 0 would get percentile 1/500 and be flagged. Rare pairs are NA, so this does not arise in the outputs.
- **Seeds.** Image i (sorted by name) uses seed 20261019 + i, logged in `results/genotype/seeds_*.csv`. The same 500 permutations are reused for each k of an image.

**S0.7 Speed.**
- The lab's `get_CLQ` takes 0.055 s for a full 9 × 9 pass on a 1,500-cell image, so 500 permutations take about 0.5 min per image per k.
- CELESTA takes about 1 s per METABRIC image.
- All images of both cohorts fit the budget, so **no subset of patients is used**.

## Stage 1: cell typing

**S1.1 Signature.** As in S0.5: `results/genotype/signature_{metabric,basel}.csv` and `thresholds_{metabric,basel}.csv` (package defaults).

**S1.2 Input scaling: decision rule fixed before the grid was run (committed with stage 0).**
- **Grid:** arcsinh cofactor {1, 5} × per-channel normalisation {none, z-score, p99}. Normalisation statistics are computed over all cells of a cohort, never per image.
  - z-score: centre and divide by the cohort SD of the arcsinh value.
  - p99: divide by the cohort 99th percentile and multiply by 10, so every channel spans roughly the 0–10 range of the package's CODEX demo data after its arcsinh(x/10).
  - One setting is chosen for both cohorts.
- **Sample:** 20 images per cohort, drawn with `set.seed(20261019)` from the sorted image list.
- **Rule:**
  - Eligible settings have a filtered fraction < 10% in both cohorts, and every round-1 class assigned (≥ 1 cell) in ≥ 80% of the 20 images in both cohorts.
  - Among eligible settings, pick the highest mean assigned fraction (not Unknown or Filtered) across the two cohorts.
  - Ties go to fewer transforms (none, then p99, then z-score), then cofactor 5 before 1.
  - If no setting is eligible, the class-presence bar is lowered to the best value reached, and that is reported.
- **Plots:** the package's `PlotExpProb` is drawn for the chosen setting and for (5, none) on one image per cohort.
- **What the rule does not use:** published labels, subtype, genotype and outcome.
- **Why these criteria:** the task asks for "assignment rate and marker-probability plots look sensible". The filtered fraction guards against settings that push every EP to one extreme. Class presence guards against settings that assign many cells but lose a lineage.

**S1.3 Scaling outcome** (`results/genotype/scaling/`). Only one setting met the pre-declared rule: **arcsinh cofactor 1, channels divided by the cohort 99th percentile × 10**.

| Setting | Mean assigned | Worst-cohort filtered | Worst class presence | Eligible |
|---|---|---|---|---|
| (1, none) | 0.666 | 0.200 | 0.85 | no |
| (1, p99) | 0.822 | 0.086 | 0.95 | **yes** |
| (1, z-score) | 0.709 | 0.129 | 0.90 | no |
| (5, none) | 0.470 | 0.408 | 0.95 | no |
| (5, p99) | 0.753 | 0.129 | 0.95 | no |
| (5, z-score) | 0.606 | 0.225 | 0.95 | no |

- **Per cohort:** METABRIC assigned 77.4%, filtered 8.6%; Basel assigned 87.0%, filtered 5.9%.
- **Why the raw scale filters so much:** with no normalisation (the package's intended input), many cells have every marker EP < 0.4 and are filtered (20–41%). The package's slope-1 sigmoid assumes CODEX-like dynamic ranges, which IMC mean ion counts do not have; rescaling to that range is what p99 × 10 does.
- **Plot check:** the `PlotExpProb` maps for the chosen setting (`scaling/ep_*_1_p99/`) show panCK-high cells in cohesive nests and sparse CD3-high cells between them.
- **Fixed for the full run:** the choice is used unchanged for all images of both cohorts.

**S1.4 Full run.**
- **METABRIC:** 749 images, 1,066,966 cells, 90 s on 18 workers. 748 images typed; 1 image had ≤ 10 cells, so the package returns NULL and its 5 cells are "Skipped".
- **Basel:** 288 images, 762,776 cells, 42 s.
- **No subset was needed.**
- **Class shares** (`celesta/class_share.csv`):

| Cohort | Tumour | T (unresolved) | CD8 T | CD4 T | Treg | B | Macrophage | Endothelial | Stroma | Unknown | Filtered |
|---|---|---|---|---|---|---|---|---|---|---|---|
| METABRIC | 43.8% | 1.5% | 3.1% | 2.6% | 2.4% | 3.8% | 8.3% | 5.5% | 7.2% | 14.3% | 7.6% |
| Basel | 47.5% | 7.8% | – | – | – | 2.5% | 11.7% | 4.0% | 10.7% | 8.3% | 7.4% |

**S1.5 Acceptance** (`celesta/acceptance.csv`).
- **Rule:** each class is assigned (≥ 1 cell) in ≥ 80% of the images where its defining marker is present.
- **"Present"** means at least 5 cells above the cohort-level crossing point from the package's own `BuildSigmoidFunction` (its 2-component GMM), fitted once per cohort on 100,000 cells drawn with seed 20261019.
- **Result: all classes pass in both cohorts.** The lowest is Treg in METABRIC at 0.970. No class was merged.

**S1.6 Concordance with the published labels** (concordance only; neither labelling is ground truth).
- **Round-1 level** (T subsets merged into T), on cells that CELESTA assigned to a shared class:
  - METABRIC: agreement 0.725, Cohen's κ 0.545.
  - Basel: agreement 0.726, κ 0.571.
- **The T-cell subclasses agree much less** (`celesta/concordance_tsubsets_metabric.csv`):
  - CD8 T: 49% of CELESTA CD8 T cells are published "CD8⁺ T cells"; 47% of published CD8⁺ T cells are recovered.
  - CD4 T: 15% and 26%.
  - **Treg: 2% and 8%.** 48.6% of CELESTA "Tregs" are cells that Danenberg et al. call epithelial. The CELESTA Treg class is 2.4% of cells, against 0.6% for the published "T_Reg & T_Ex" class.
- **Likely cause:** the package min-max rescales EPs within each image, so every image has some cells near EP = 1 for FOXP3 and CD3. The default "low" anchor threshold of 0.9 also lets panCK-intermediate cells into the T branch.
- **Decision (made before any genotype or subtype was linked):**
  - The signature and thresholds are **not** retuned against the published labels. They were fixed in advance, and the published labels are the comparison, not a target.
  - Results that depend on the Treg class are marked as such: P3 (Treg–T CLQ by BRCA1) and question (d) (Tumour–Treg CLQ by node status). The pre-planned refit on the published labels, which have their own Treg class, is the check for them.
  - The CD8 T class, used for H1, is moderately concordant.
