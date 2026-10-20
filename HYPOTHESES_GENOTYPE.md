# Part IV: pre-registered hypotheses and analysis plan

Committed before any subtype, receptor, mutation or node label was linked to a spatial feature or to composition. **This file is not edited after that commit.** Any later change to the analysis is reported as a deviation in the README and in AMBIGUITIES_GENOTYPE.md. Sources for the expectations are in PRIOR_EXPECTATIONS.md.

## Data and features (as fixed in stages 0–2)

**Cell types.** Official CELESTA, with the signatures `results/genotype/signature_{metabric,basel}.csv`.
- METABRIC has 9 classes: Tumour, T cell (unresolved), B cell, Macrophage, Endothelial, Stroma, CD8 T cell, CD4 T cell, Treg.
- "Round-1" labels (T subsets merged into T cell) give the 6 classes shared with Basel.

**CLQ.** The Plevritis lab's `get_CLQ`, divided by k, so the random-labelling null is 1.
- Unknown and Filtered cells are removed before the kNN graph is built.
- "CLQ A→B" means: B cells among the k nearest neighbours of A cells, relative to random labelling. In the output files the centre cell is `cell_A` and the neighbour is `cell_B`. This is the CELESTA paper's CLQ_{B→A}.
- k ∈ {10, 20, 40}; k = 20 is the lab default and the primary value.

**Patient features.**
- **CLQ:** the mean of image CLQs across the patient's images, weighted by the image's typed-cell count. Images where the pair is undefined (≤ 5 cells of either type) are left out.
- **Composition:** the CLR of the patient's class counts (pseudocount 0.5, as in `spatialsurv.features.clr`). All but the last component enter the models, because CLR components sum to zero.
- **Size:** log of the patient's typed-cell count.

**Inclusion.**
- Patients with ≥ 500 typed cells (the repo's METABRIC `MIN_CELLS`), with age and grade known.
- A pair's model uses the patients in whom that pair is defined, and is fitted only if ≥ 50 patients qualify.

**Labels.**
- **PAM50:** LumA (reference), LumB, Her2 (HER2-enriched), Basal, claudin-low. Normal-like, NC and missing are excluded from Q1.
- **Mutations:** non-silent on METABRIC_173; missing if the patient was not sequenced.
- **ER:** cBioPortal ER_STATUS.
- **Node status:** ≥ 1 positive lymph node.
- **Basel surrogate subtypes:** ER⁺/HER2⁻ (reference), HER2⁺, TN.

## Models

All models are ordinary least squares (R `lm`), one per ordered pair and k:

    CLQ ~ group + CLR composition + log(typed cells) + age + grade

Effects are coefficients (CLQ units) with 95% CIs.

**Unadjusted versions** drop the CLR terms and keep log cells, age and grade. They show how much of each association is composition.

**Q1, subtypes (METABRIC).**
- Model: group = PAM50 factor.
- Test per pair × k: the F-test of the subtype block against the model without it.
- BH across all 81 pairs × 3 k = 243 tests. Contrasts against LumA are reported with CIs.

**Q2, mutations (METABRIC, sequenced patients).**
- Model: group = TP53 (0/1), with PAM50 entered as a covariate. Normal-like, NC and missing form one "other" level, so those patients are kept.
- Test: the t-test of the TP53 coefficient. BH across 243 tests.
- **PIK3CA** is a separate family, with the same model and its own BH correction.
- Within-subtype estimates (TP53 in LumA, LumB, Her2 and claudin-low separately) are reported descriptively. Basal has 2 TP53-wild-type patients, so no within-Basal estimate is made.

**Replication family (H3).**
- METABRIC on the round-1 labels: 36 pairs × 3 k. The contrasts are Basal − LumA and Her2 − LumA, which have receptor-defined counterparts: TN − ER⁺/HER2⁻ and HER2⁺ − ER⁺/HER2⁻.
- BH across 36 × 3 × 2 = 216 contrast t-tests.
- Every contrast with q < 0.05 is refitted in Basel with the surrogate subtypes, at the same pair and k.

**Robustness.**
1. **Label permutation (METABRIC, Q1 family):** PAM50 labels are shuffled across patients 1,000 times (seed 20261020) and all 243 F-tests refitted.
   - **Calibrated if** ≤ 10% of permutations give any BH discovery (nominal ≈ 5%) **and** the mean null p-value is within 0.5 ± 0.02.
   - A permutation p-value is also reported for each observed F.
2. **Leave-one-k-out:** for each Q1 discovery, report whether the subtype F-test has p < 0.05 at both other k.
3. **Alternative cell typing:** Q1 is refitted at k = 20 on CLQs from the published labels (Danenberg phenotypes mapped to the same 9 classes; same lab code and settings). Report the agreement of F-statistics across pairs and of the Basal − LumA signs.

**Positive controls.** Defined in PRIOR_EXPECTATIONS.md:
- P1: immune fraction by subtype.
- P2a–c: Basal versus LumA for Tumour→Tumour and Stroma→Stroma; B↔T; Treg→CD8 T and Treg→CD4 T.
- P3: BRCA1, descriptive.

All are at k = 20 with the adjusted model, except P1, which is a composition model:

    immune fraction ~ subtype + log cells + age + grade

## Hypotheses and my numeric predictions

**Scale.** CLQ = 1 is random. A between-patient SD of tumour-centred heterotypic CLQs of about 0.2–0.3 is expected. These are my guesses before seeing the data.

| ID | Hypothesis (direction) | Exact test | Pass rule | My prediction | P(pass), my guess |
|---|---|---|---|---|---|
| H1 | Basal-like tumours show higher CD8 T-cell–tumour co-localisation than LumA, after adjusting for composition and image size | METABRIC, CLQ Tumour→CD8 T cell (CD8 T cells among tumour cells' neighbours), k = 20, Q1 model, Basal − LumA coefficient | coefficient > 0 with two-sided p < 0.05 (single pre-specified test; its Q1-family q is also reported) | +0.10 (plausible range −0.05 to +0.25) | 0.40 |
| H1-rev (secondary) | same, reverse direction | CLQ CD8 T cell→Tumour | as H1 | +0.08 | – |
| H2 | Within ER⁺ tumours, TP53-mutant tumours show higher macrophage–tumour co-localisation | METABRIC ER⁺ sequenced patients, CLQ Tumour→Macrophage, k = 20, Q2 model (TP53 + PAM50/other + CLR + log cells + age + grade) | TP53 coefficient > 0, p < 0.05 | +0.05 (range −0.05 to +0.15) | 0.30 |
| H3 | ≥ 1 METABRIC subtype–CLQ association with q < 0.05 has the same sign in Basel (surrogate subtypes) | replication family above | ≥ 1 same-sign | 6–15 METABRIC discoveries; sign agreement 70–85%; 2–5 with Basel p < 0.05 | 0.85 |
| H3-strict (secondary, mine) | sign agreement is better than chance | the same discoveries | one-sided binomial test of the same-sign share > 0.5, p < 0.05 (needs ≥ 5 discoveries) | as above | 0.55 |
| P1 | Basal > LumA and Her2 > LumA in immune fraction | composition model | both > 0, CI excluding 0 | Basal +0.08, Her2 +0.05 (absolute fraction) | 0.85 |
| P2a | Tumour→Tumour, Stroma→Stroma higher in Basal | Q1 model, k = 20 | both coefficients > 0 | +0.05 and +0.10 | 0.55 |
| P2b | B↔T higher in Basal | round-1 labels, Q1 model, k = 20 | B→T and T→B coefficients > 0 | +0.15 | 0.45 |
| P2c | Treg→CD8 T and Treg→CD4 T higher in Basal | Q1 model, k = 20; repeated on published labels | both > 0 | +0.10 | 0.40 (CELESTA Treg class has low concordance) |
| P2 overall | ≥ 2 of P2a–c pass | | | | 0.45 |
| (d) | Tumour↔Treg co-localisation higher in node-positive tumours (Zhang 2022 question) | METABRIC, CLQ Tumour→Treg and Treg→Tumour, k = 20, model: node-positive + PAM50/other + CLR + log cells + age + grade; CELESTA and published labels | Tumour→Treg coefficient > 0, p < 0.05 (CELESTA labels); published labels reported alongside | +0.03 (range −0.05 to +0.10) | 0.20 |
| Keren (descriptive) | within Basal / TN, Tumour→T CLQ varies widely | IQR of patient Tumour→T CLQ (round-1) in Basal vs LumA | none | Basal IQR wider | – |
| Jackson (descriptive) | immune cells are excluded from tumour mass; T cells sit near endothelium | median patient CLQ Tumour→T, Tumour→Macrophage < 1 in every subtype; T→Endothelial > 1 | none | yes, both | – |

## What would change my reading

**A pass for H1, H2 or (d) is supporting evidence only when:**
- the unadjusted and adjusted estimates do not reverse sign;
- the effect has the same sign at k = 10 and k = 40;
- for (d) and P2c, the published-label refit has the same sign.

**Interpretation of nulls.** With 42 Basal and 208 LumA patients, a single-core CLQ difference of about 0.1 is near the detectable limit. A null is reported as "not detected at this sample size", not as absence of an effect.

**Nothing here is prognostic.** The study is cross-sectional and associative.
