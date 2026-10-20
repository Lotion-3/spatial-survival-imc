# Part IV: what prior work predicts

Written before any subtype or genotype label was linked to a spatial feature, and committed with HYPOTHESES_GENOTYPE.md. For each paper this file records what was read, what it reports that this analysis can represent, and the direction it predicts. The positive controls P1–P3 are fixed here, with their sources.

**Our cell classes.**
- **Both cohorts:** Tumour, T cell, B cell, Macrophage, Endothelial, Stroma (SMA⁺).
- **METABRIC only:** the T cells are further split into CD8 T, CD4 T and Treg (FOXP3⁺).

The METABRIC CELESTA Treg class has low concordance with the published Treg class (AMBIGUITIES S1.6), so Treg-based controls are also run on the published labels.

**Our spatial measure.** The CLQ of a ordered pair a→b. It is the Plevritis lab's CLQ divided by k: the number of type-b cells among the k nearest neighbours of type-a cells, relative to what random labelling would give (null = 1).

## Ali et al., *Nature Cancer* 2020 (METABRIC, 483 tumours, IMC + genomics)

Read: the accepted manuscript (Cambridge Apollo repository), main text and figure legends.

- **Composition (Results, "Genomic subtypes … diverse tumour ecosystems"; Fig. 4):** "We noted both T cell and macrophage enrichment among Basal-like and IntClust 10 tumours." Luminal B / IntClust 9 tumours were the only ER⁺ group with both macrophage and T-cell enrichment. Myofibroblasts were enriched in indolent ER⁺ tumours (Luminal A, IntClust 3 and 4⁺). → **Predicts:** a higher T-cell and macrophage fraction in Basal-like than in Luminal A tumours (P1), and a higher Stroma (SMA⁺) fraction in Luminal A.
- **Interactions by subtype (Results; Extended Data Figs. 5b and 6):** "Subtypes significantly enriched for interactions included HER2, Basal-like, and IntClust 10. In Basal-like and IntClust 10 tumours, we observed abundant homotypic relationships among both epithelial and stromal cells … a starker separation between compartments … Basal-like and IntClust 10 tumours were associated with more homotypic interactions, also suggestive of a 'separation phenotype'." → **Predicts (P2a):** Tumour→Tumour and Stroma→Stroma CLQ are higher in Basal-like than in Luminal A tumours.
- **Mutations (Fig. 5):** TP53 mutation was associated with ER⁻ basal, hypoxia-associated and HR⁻ Ki-67⁺ epithelial phenotypes and with fibroblast phenotypes 30 and 37 and myofibroblast phenotype 32. PIK3CA mutation was associated with HR⁺ CK7⁻ epithelial cells. These are composition (phenotype-proportion) results. **No interaction-by-mutation result is reported**, so TP53/PIK3CA spatial effects (H2) are new here.
- **Method:** neighbours were touching cells (CellProfiler object relationships), with 1,000 label permutations per image (the histoCAT approach). This differs from our kNN CLQ.
- **Code:** "In house image preprocessing scripts are available at github.com/BodenmillerGroup/imctools. Other analysis code is available from the authors upon request." No public analysis code, so the positive controls are reproduced from the paper.

## Danenberg et al., *Nature Genetics* 2022 (METABRIC, 693 tumours; the single-cell data used here)

Read: PMC7612730 full text; the analysis code in Zenodo 10.5281/zenodo.6036188 (R scripts, listed and read by range request).

- **Ten TME structures (Fig. 4; Extended Data Fig. 6):** FSP1⁺ enriched, active stroma, PDPN⁺ active stroma, vascular stroma, active IR, CD8⁺ & macrophages, granulocyte enriched, APC enriched, TLS-like, and suppressed expansion.
- **Subtype enrichment (Extended Data Fig. 8; code `ExtFig8/glmsTMEStructuresVTumourSubtype.R`):** a binomial GLM per structure × subtype (one versus rest), with BH correction within each taxonomy (PAM50, IntClust). The text reports:
  - ER status has the greatest influence on TME structure.
  - **Suppressed expansion** (regulatory and dysfunctional T cells co-occurring) is enriched in ER⁻ tumours (and IntClust 2) and "significantly depleted" in Luminal A.
  - **TLS-like** structures (B and T cells together) are characteristic of basal-like tumours.
  - **Active IR** is characteristic of Luminal A.
  - → **Predicts (P2b):** B→T and T→B co-localisation is higher in Basal-like than in Luminal A tumours (TLS-like).
  - → **Predicts (P2c):** Treg↔T-cell co-localisation is higher in Basal-like than in Luminal A tumours (suppressed expansion). CELESTA labels give Treg→CD8 T and Treg→CD4 T; the check is the same contrast on the published labels.
- **Mutations (Extended Data Fig. 9; code `ExtFig9/TMEStructuresVDriversGLMs.R`):** suppressed expansion is enriched for BRCA1 and CASP8 mutations. → **Predicts (P3):** Treg↔T co-localisation is higher in BRCA1-mutant tumours. Only **6 IMC patients are BRCA1-mutant** on the METABRIC_173 panel, so P3 can only be reported descriptively (AMBIGUITIES S0.4).
- **TP53:** TP53 mutation is better predicted by TME network properties (AUC 0.72) than by tumour-cell phenotypes (0.65) (Fig. 5). Some spatial signal for TP53 is therefore expected, but its direction and pairs are not reported.
- **Composition:** TME cell diversity is higher in ER⁻ subtypes (basal; IntClust 4, 5, 10) than in Luminal A (Fig. 2d).
- **Neighbour results by subtype:** "no neighbour results stratified by subtype, ER status or IntClust" appear in the text. The interface analyses (Fig. 3) are pooled: lymphocytes are depleted at the tumour–stroma interface, and endothelial and CD38⁺ cells are enriched perivascularly.
- **Lymph-node status:** not analysed.
- **Code use:** their code works on their own TME-structure assignments (32 phenotypes, graph communities via igraph). It does not compute our controls on our classes. We therefore use their reported directions, not their code (AMBIGUITIES pipeline table).

## Denkert et al., *Lancet Oncology* 2018 (3,771 patients, six neoadjuvant GBG trials)

Read: the abstract and published summaries (PubMed 29233559).

- **TIL levels by subtype:** stromal TILs were assessed on H&E (low 0–10%, intermediate 11–59%, high ≥ 60%). High TILs were most frequent in TNBC (273 of 906, 30%) and less frequent in HER2-positive and luminal-HER2-negative disease.
- **TILs and outcome:** a 10% TIL increase was associated with longer DFS in TNBC (HR 0.93, 95% CI 0.87–0.98) and HER2-positive disease (0.94, 0.89–0.99), but not in luminal-HER2-negative disease (1.02, 0.96–1.09), where higher TILs went with shorter OS.
- → **Predicts (P1):** a higher immune fraction in Basal-like, HER2-enriched, TN and HER2⁺ tumours than in Luminal A / ER⁺HER2⁻ tumours.
- TILs are a stromal composition measure (amount, not arrangement), so Denkert predicts nothing about CLQ. This study makes no prognostic claims.

## Jackson et al., *Nature* 2020 (Basel/Zurich IMC, 281 Basel tumours with follow-up; the Basel data used here)

Read: the accepted manuscript (Cambridge Apollo repository).

- **"Across all patients, immune cells were excluded from the cohesive tumor mass"** (Fig. 1c; Extended Data Fig. 4). → **Predicts:** Tumour→immune CLQ below 1 in every subtype.
- **"T cells were surrounded by endothelial cells"** (Fig. 2b). → **Predicts:** T→Endothelial CLQ above 1 (descriptive; not by subtype).
- **Stromal environments by subtype (Extended Data Figs. 7 and 8b):** "TNBCs … were associated with T cell-enriched (SEs 5 and 8) or macrophage-enriched (SE 3) stromal environments. HR+ tumors are likely to be immune-cold, but some contained rare and localized immune-enriched communities … the HR+ tumors never contained an exclusively immune-hot stromal [environment]." → **Predicts (P1, Basel replication):** a higher T and macrophage fraction in TN than in ER⁺/HER2⁻.
- Single-cell pathology subgroups cut across clinical subtypes. They are defined from tumour-cell phenotypes and communities, not from pairwise co-localisation.

## Keren et al., *Cell* 2018 (MIBI-TOF, 41 triple-negative tumours)

Read: from prior knowledge of the paper (not re-read here). The claims below are limited to its central, widely cited findings.

- **Three architectures:** "cold" (few immune cells), "mixed" (immune cells intermixed with tumour cells) and "compartmentalised" (tumour and immune cells in separate regions). Compartmentalised tumours had better overall survival.
- → **Predicts (descriptive only):** within Basal-like and TN tumours, tumour–immune CLQ should vary widely between patients. Mixed tumours sit near 1; compartmentalised tumours sit well below 1.
- Keren compares no subtypes, so it predicts no subtype difference and no H1 direction. We report the spread of Tumour→T CLQ within Basal/TN tumours; it gets no pass/fail.

## Zhang et al., *Nature Methods* 2022 (CELESTA), and question (d)

Read: the PMC9728133 full text.

- **The original finding:** with the CLQ CLQ_{b→a} = (C_{b→a}/N_a)/(N_b/(N−1)), "malignant cells and T-regulatory cells (Tregs)" were more co-localised in lymph-node-positive than node-negative HNSCC.
  - CODEX: 4 vs 4 patients, nominal P = 0.038 (two-sided t-test) (Fig. 5b–d).
  - IHC validation: cytokeratin–FOXP3 density correlation, 44 vs 32 patients, P = 0.011 (Fig. 5f–g).
  - Also reported as higher in N⁺ tumours: CD4 T–endothelial, CD8 T–CD4 T, and CD4 T–CD4 T co-localisation.
- **Breast cancer:**
  - Neither Ali 2020 nor Danenberg 2022 analyses node status spatially.
  - A literature search (web, October 2026; not exhaustive) found mIF studies of Treg density and of CD8–tumour proximity, but **no breast-cancer study testing Tumour–Treg co-localisation against node status**.
  - Question (d) is therefore treated as new. Its direction is taken from Zhang: Tumour↔Treg CLQ higher in node-positive tumours.

## Positive controls (fixed here)

| Control | Source | Our measure (METABRIC unless stated) | Predicted direction | Pass rule |
|---|---|---|---|---|
| P1 | Ali 2020 (Fig. 4 text), Denkert 2018, Jackson 2020 (ED 7–8) | immune fraction = (T + CD8 T + CD4 T + Treg + B + Macrophage) / typed cells, per patient | Basal > LumA and HER2-enriched > LumA | both contrasts > 0 with 95% CI excluding 0 (linear model with age, grade, log cell count) |
| P2a | Ali 2020, ED Figs. 5b, 6 ("separation phenotype") | CLQ Tumour→Tumour; CLQ Stroma→Stroma | Basal > LumA | composition-adjusted Basal − LumA coefficient > 0 for both (k = 20) |
| P2b | Danenberg 2022, ED Fig. 8 (TLS-like in basal) | CLQ B→T cell and T cell→B (T cell = all T classes pooled at round 1) | Basal > LumA | coefficient > 0 for both (k = 20) |
| P2c | Danenberg 2022, ED Fig. 8 (suppressed expansion enriched in ER⁻, depleted in LumA) | CLQ Treg→CD8 T and Treg→CD4 T | Basal > LumA | coefficient > 0 for both (k = 20); also on the published labels |
| P3 | Danenberg 2022, ED Fig. 9 (BRCA1) | CLQ Treg→CD8 T and Treg→CD4 T | BRCA1-mutant > wild-type | descriptive only (6 mutant patients); no pass/fail |

- **P2 overall:** P2 passes if at least two of P2a, P2b and P2c pass. This is the "direction is reproduced" requirement of the task, applied to the three subtype-specific interactions our classes can represent.
- **P2b pooling:** CLQ for P2b, and for the METABRIC side of the Basel replication, is computed on CELESTA's round-1 labels. There the T subsets are still "T cell", because CELESTA assigns them in round 2. The run uses the same lab function and settings. This gives the six classes shared with Basel.
