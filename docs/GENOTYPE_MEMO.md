# Memo: do breast-cancer subtype and genotype shape spatial organisation?

**Bottom line.** After accounting for how many cells of each type a tumour contains, we found no subtype or TP53/PIK3CA difference in how those cells are arranged. One arrangement replicated across both cohorts: macrophages and endothelial cells sit closer together in basal-like and triple-negative tumours. The published positive controls were only partly reproduced, so the nulls are best read as "not detectable with this pipeline in single cores", not as evidence of no effect. Nothing here is prognostic.

## Question

Using the official CELESTA cell typing and the Plevritis lab's co-location quotient (CLQ), we asked three things:
1. Do PAM50 subtypes (with claudin-low) differ in cell–cell co-location beyond composition?
2. Is TP53 or PIK3CA mutation associated with co-location within subtype?
3. Do the METABRIC subtype differences replicate in Basel with receptor-based subtypes?

Hypotheses and pass rules were committed before any label was linked.

## Data and methods

- **METABRIC IMC (Danenberg 2022):** 749 tumour images and 693 patients, with PAM50, the 173-gene mutation panel and clinical data from cBioPortal.
- **Basel IMC (Jackson 2020):** 288 images and 280 patients.
- **Cell typing:** the official CELESTA R package, unmodified, on the six lineage classes both panels share. METABRIC also has CD8, CD4 and Treg subsets.
- **Co-location:** the lab's own CLQ functions, unmodified (k = 10/20/40 nearest neighbours).
- **Models:** patient-level linear models adjusted for CLR composition, cell count, age and grade, with Benjamini–Hochberg within each question.

## Main associations

| Question | Result |
|---|---|
| Subtype → arrangement (243 tests) | 0 discoveries after composition adjustment; 7 before it (immune and endothelial cells around tumour cells). Composition explains them. |
| H1: more CD8 T–tumour co-location in Basal than LumA | No: −0.05 [−0.16, 0.06] |
| H2: TP53 and macrophage–tumour co-location in ER⁺ tumours | No: −0.01 [−0.06, 0.05] |
| TP53 and PIK3CA, all pairs | No FDR-significant effects. Strongest TP53 lead: fewer tumour cells around CD8 T cells in mutant tumours (−0.09, p = 0.003, q = 0.27; exploratory). |
| Node status and tumour–Treg co-location (Zhang 2022's question) | No with CELESTA labels (−0.01). Published labels: +0.13 [−0.03, 0.28], n = 88. |

## What replicated

- **Replicated in Basel:** macrophage–endothelial co-localisation higher in basal-like (METABRIC, +0.29 to +0.45) and triple-negative (Basel, +0.11 to +0.18) tumours, at all three k.
- **Not replicated:** METABRIC's lower macrophage–tumour co-location in basal-like tumours (Basel ≈ 0).
- **Known composition differences reproduced:** a higher immune fraction in basal-like (+0.11) and triple-negative (+0.11) tumours.
- **Not reproduced:**
  - The HER2-enriched immune excess. It is present with the published labels but not with CELESTA labels.
  - Ali 2020's homotypic "separation phenotype".
  - Danenberg 2022's Treg–T co-occurrence in ER⁻ tumours.

## Limitations

- **Single TMA cores.** Usually one small core per patient, so co-location samples a few hundred µm of each tumour.
- **Targeted mutation panel.** 173 genes. BRCA1 has 6 mutant patients, and only 2 basal tumours are TP53 wild-type, so those effects cannot be tested.
- **Cross-sectional association only.** No causal or prognostic interpretation.
- **IMC scaling choice.** The package's sigmoid assumes CODEX-like intensity ranges, so IMC counts were rescaled by a rule fixed in advance. Other reasonable scalings would shift the labels.
- **Cell typing.** The shared six-marker panel is small. The CELESTA Treg class agrees poorly with the published labels, and the published-label refit agrees only loosely overall (F correlation 0.18).
- **Translating prior results.** Earlier studies used 32–57 clustered phenotypes and different neighbourhood definitions, so mapping their results onto our classes is approximate.
- **Compute-driven deviation.** Per-image permutations (500) were run at k = 20 only. The association models use observed CLQs at every k.

Full details: README Part IV, `HYPOTHESES_GENOTYPE.md`, `PRIOR_EXPECTATIONS.md`, `AMBIGUITIES_GENOTYPE.md`.
