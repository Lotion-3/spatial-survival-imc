# Memo: do breast-cancer subtype and genotype shape spatial organisation?

**Summary.** This memo asks whether breast-cancer subtype and genotype shape how tumour, immune and stromal cells are arranged, beyond how many of each there are. Cells were typed with the official CELESTA R package (plevritis-lab/CELESTA @ `e5e92ff`), and co-location was measured with the Plevritis lab's CLQ code (Spatial_Permutation_and_Normalization @ `02b6da5`). Both were used unmodified.
- **Omnibus subtype tests:** none survives composition adjustment (0 of 243).
- **Pre-specified Basal − LumA contrasts:** three ordered pairs reached q < 0.05 at all k.
- **Replication:** of those three, macrophage–endothelial co-localisation replicated in Basel. It is higher in basal-like (METABRIC) and triple-negative (Basel) tumours.
- **Pre-registration:** hypotheses, numeric predictions and pass rules were committed before any label was linked.

The published positive controls were only partly reproduced with this setup, so the nulls are best read as "not detectable with this signature, scaling and single cores", not as evidence of no effect. Nothing here is prognostic.

## Question

Using CELESTA cell types and CLQ co-location, the analysis asks three questions:
1. Do PAM50 subtypes (with claudin-low) differ in cell–cell co-location beyond composition?
2. Is TP53 or PIK3CA mutation associated with co-location within subtype?
3. Do the METABRIC subtype differences replicate in Basel with receptor-based subtypes?

## Data and methods

- **METABRIC IMC (Danenberg 2022):** 749 tumour images, 693 patients, with PAM50, the 173-gene mutation panel and clinical data from cBioPortal.
- **Basel IMC (Jackson 2020):** 288 images, 280 patients.
- **Cell typing:** CELESTA on the six lineage classes both panels share; METABRIC also has CD8, CD4 and Treg subsets.
- **Input scaling:** IMC counts were rescaled before input so that marker intensities span the range the package's expression-probability step expects. The rule was fixed in advance (AMBIGUITIES_GENOTYPE.md).
- **Co-location:** the lab's CLQ functions at k = 10/20/40 nearest neighbours.
- **Models:** patient-level linear models adjusted for CLR composition, cell count, age and grade, with Benjamini–Hochberg within each pre-specified family.

## Main associations

| Question | Result |
|---|---|
| Subtype → arrangement, omnibus F-test (9 classes, 243 tests) | 0 discoveries after composition adjustment; 7 before (immune and endothelial cells around tumour cells), so composition explains these |
| Subtype → arrangement, pre-specified Basal/HER2-enriched − LumA contrasts (6 shared classes, 216 tests) | 9 with q < 0.05, all Basal − LumA. These are 3 ordered pairs, each significant at all three k: macrophage → endothelial, endothelial → macrophage (higher) and macrophage → tumour (lower). The omnibus test is broader and has a larger family, which is why it gives macrophage → tumour q = 0.11 against q = 0.015 here. |
| H1: more CD8 T–tumour co-location in Basal than LumA | not supported: −0.05 [−0.16, 0.06] |
| H2: TP53 and macrophage–tumour co-location in ER⁺ tumours | not supported: −0.01 [−0.06, 0.05] |
| TP53 and PIK3CA, all pairs | no FDR-significant effects. Strongest TP53 lead: fewer tumour cells around CD8 T cells in mutant tumours (−0.09, p = 0.003, q = 0.27; exploratory) |
| Node status and tumour–Treg co-location (the question from Zhang 2022) | not supported with this pipeline's labels (−0.01). Published labels: +0.13 [−0.03, 0.28], n = 88. We found no earlier breast-cancer test of this question; here it is null. |

## What replicated

- **Replicated in Basel:** macrophage–endothelial co-localisation is higher in basal-like (METABRIC, +0.29 to +0.45) and triple-negative (Basel, +0.11 to +0.18) tumours, at all three k.
- **Not replicated:** METABRIC's lower macrophage–tumour co-location in basal-like tumours (Basel ≈ 0).
- **Known composition differences reproduced:** a higher immune fraction in basal-like (+0.11) and triple-negative (+0.11) tumours.
- **Not reproduced:**
  - The HER2-enriched immune excess. It is reproduced with the published labels (+0.075) but not with this pipeline's labels (+0.040), so it depends on the cell typing.
  - Ali 2020's homotypic "separation phenotype".
  - Danenberg 2022's Treg–T co-occurrence in ER⁻ tumours.

## Limitations

- **Single TMA cores.** Usually one small core per patient, so co-location samples only a few hundred µm of each tumour.
- **Targeted mutation panel.** 173 genes. BRCA1 has 6 mutant patients and only 2 basal tumours are TP53 wild-type, so those effects cannot be tested.
- **Cross-sectional association only.** No causal or prognostic interpretation.
- **IMC scaling choice.** IMC mean ion counts span a narrower range than the CODEX-like intensities the expression-probability step matches, so they were rescaled by a pre-declared rule. Other reasonable scalings would shift the labels.
- **Cell typing with this setup.** The shared six-marker signature is small. With this signature and IMC scaling, the Treg class agrees poorly with the published labels (about half of its cells are published as epithelial). The published-label refit agrees only loosely overall (F correlation 0.18).
- **Mapping prior results.** Earlier studies used 32–57 clustered phenotypes and different neighbourhood definitions, so the comparison with our classes is approximate.
- **Compute-driven deviation.** Per-image permutations (500) were run at k = 20 only. The association models use observed CLQs at every k.

Full details: README Part IV, `HYPOTHESES_GENOTYPE.md`, `PRIOR_EXPECTATIONS.md`, `AMBIGUITIES_GENOTYPE.md`.
