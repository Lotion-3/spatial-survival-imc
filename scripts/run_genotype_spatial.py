"""Part IV cohort build: per-cell marker tables and patient tables for the R steps.

Data loading reuses src/spatialsurv (unchanged); everything downstream (CELESTA, CLQ, models,
figures) is R, reading the files written here (CSV with explicit column names).

    uv run python scripts/run_genotype_spatial.py              # both cohorts
    uv run python scripts/run_genotype_spatial.py metabric     # one cohort

Writes
    data/processed/genotype/{metabric,basel}_cells.csv.gz   cell_id, image, X, Y, markers..., published
    results/genotype/patients_{metabric,basel}.csv          one row per patient, labels and covariates
    results/genotype/cohort_counts.json                     stage-0 acceptance counts
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from spatialsurv.data import load_cells, load_patients
from spatialsurv.metabric import load_metabric_cells, load_metabric_clinical

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW_MB = ROOT / "data" / "raw_metabric"
PROC = ROOT / "data" / "processed" / "genotype"
OUT = ROOT / "results" / "genotype"

# Marker column names used in both cohorts (METABRIC SingleCells.csv name -> common name).
SHARED = {"panCK": "panCK", "CK8-18": "CK8-18", "CK5": "CK5", "CD3": "CD3", "CD20": "CD20", "CD68": "CD68",
          "CD31-vWF": "CD31-vWF", "SMA": "SMA", "HER2 (3B5)": "HER2", "Ki-67": "Ki-67", "Histone H3": "Histone H3"}
METABRIC_ONLY = {"CD8": "CD8", "CD4": "CD4", "FOXP3": "FOXP3", "CD45RA": "CD45RA", "CD45RO": "CD45RO",
                 "PDGFRB": "PDGFRB", "FSP1": "FSP1", "Podoplanin": "Podoplanin", "ER": "ER"}
BASEL_ONLY = ["CD45", "Vimentin", "Fibronectin"]

GENES = ["TP53", "PIK3CA", "BRCA1", "CASP8"]
# cBioPortal mutationType values counted as a mutation (non-silent, protein-altering or splice).
NON_SILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
              "In_Frame_Ins", "Splice_Site", "Splice_Region", "Translation_Start_Site", "Nonstop_Mutation"}
# Published labels at the CELESTA class resolution (robustness: alternative cell typing). Non-T classes
# follow the repo's PHENOTYPE_CLASS (part II, Amendment 1); T phenotypes are split where CELESTA splits.
COARSE_NAME = {"Tumor": "Tumour", "T": "T cell", "B": "B cell", "Macrophage": "Macrophage",
               "Endothelial": "Endothelial", "Stroma": "Stroma"}
T_FINE = {"CD8^{+} T cells": "CD8 T cell", "CD4^{+} T cells": "CD4 T cell", "T_{Reg} & T_{Ex}": "Treg"}
PAM50 = {"LumA": "LumA", "LumB": "LumB", "Her2": "Her2", "Basal": "Basal", "claudin-low": "claudin-low"}


def metabric_cells() -> pd.DataFrame:
    base = load_metabric_cells(RAW_MB)  # tumour images, hot aggregates removed, published coarse class
    cols = ["is_tumour", "is_hotAggregate", "Location_Center_X", "Location_Center_Y", *SHARED, *METABRIC_ONLY]
    m = pd.read_csv(RAW_MB / "SingleCells.csv", usecols=cols)
    m = m[(m["is_tumour"] == 1) & (m["is_hotAggregate"] == 0)].reset_index(drop=True)
    assert np.allclose(m["Location_Center_X"], base["x"]) and np.allclose(m["Location_Center_Y"], base["y"])
    m = m.rename(columns=SHARED | METABRIC_ONLY)
    out = pd.DataFrame({"cell_id": np.arange(len(base)), "image": base["core"], "PID": base["PID"],
                        "X": base["x"], "Y": base["y"]})
    out[list((SHARED | METABRIC_ONLY).values())] = m[list((SHARED | METABRIC_ONLY).values())]
    out["published"] = base["coarse"]
    out["published_phenotype"] = base["cellPhenotype"]
    fine = base["cellPhenotype"].map(T_FINE).where(base["is_epithelial"] == 0)
    out["published_fine"] = fine.fillna(base["coarse"].map(COARSE_NAME))
    return out


def basel_cells(cores: pd.DataFrame) -> pd.DataFrame:
    cells = load_cells(RAW, cores=list(cores["core"]))  # published metacluster labels + coordinates
    expr = pd.read_csv(RAW / "basel_expression.csv.gz")
    expr["id"] = expr["core"] + "_" + expr["CellId"].astype(str)
    c = cells.merge(expr.drop(columns=["core", "CellId"]), on="id", how="inner")
    c = c.merge(cores, on="core", how="left")
    out = pd.DataFrame({"cell_id": np.arange(len(c)), "image": c["core"], "PID": c["PID"].astype(str),
                        "X": c["x"], "Y": c["y"]})
    for m in list(SHARED.values()) + BASEL_ONLY:
        out[m] = c[m]
    out["published"] = c["coarse"]
    out["published_phenotype"] = c["cluster"].astype(str)
    out["published_fine"] = c["coarse"].map(COARSE_NAME)
    out["published_id"] = c["id"]
    out.attrs["n_location_cells"] = len(cells)
    return out


def mutation_status(sequenced: set[str], muts: list[dict]) -> pd.DataFrame:
    """1 = non-silent mutation, 0 = sequenced without one, NaN = not sequenced."""
    df = pd.DataFrame(muts)
    df = df[df["hugoGeneSymbol"].isin(GENES) & df["mutationType"].isin(NON_SILENT)]
    status = pd.DataFrame(0.0, index=sorted(sequenced), columns=GENES)
    for g in GENES:
        status.loc[status.index.isin(df.loc[df["hugoGeneSymbol"] == g, "sampleId"]), g] = 1.0
    return status


def metabric_patients(pids: pd.Index) -> pd.DataFrame:
    clin = load_metabric_clinical(RAW_MB)  # age, grade, pN, ER, PR, HER2 (existing harmonisation)
    raw = pd.DataFrame.from_dict(json.loads((RAW_MB / "metabric_clinical.json").read_text()), orient="index")
    seq = json.loads((RAW_MB / "metabric_sequenced_samples.json").read_text())
    muts = json.loads((RAW_MB / "metabric_mutations.json").read_text())
    panels = json.loads((RAW_MB / "metabric_gene_panels.json").read_text())
    samples = pd.Series(seq["sampleIds"])
    # brca_metabric has one sample per patient and sampleId == patientId (checked here)
    assert samples.isin(raw.index).mean() > 0.99, "sampleId != patientId in brca_metabric"
    status = mutation_status(set(samples), muts)
    for g in GENES:  # a gene missing from a sample's panel cannot be called wild-type
        on_panel = {s for s, p in seq["genePanelId"].items() if p is None or g in panels.get(p, [g])}
        status.loc[~status.index.isin(on_panel), g] = np.nan
    p = pd.DataFrame(index=pd.Index(sorted(pids), name="PID"))
    p["pam50"] = raw["CLAUDIN_SUBTYPE"].reindex(p.index).map(PAM50)  # Normal / NC / missing -> NaN
    p["pam50_raw"] = raw["CLAUDIN_SUBTYPE"].reindex(p.index)
    p["intclust"] = raw["INTCLUST"].reindex(p.index)
    p = p.join(clin[["age", "grade", "pN", "ER", "PR", "HER2"]])
    p["node_positive"] = np.where(p["pN"].isna(), np.nan, (p["pN"] > 0).astype(float))
    p["sequenced"] = p.index.isin(status.index)
    p = p.join(status)
    return p


def basel_patients() -> tuple[pd.DataFrame, pd.DataFrame]:
    pat, cores = load_patients(RAW)  # existing loader (drops patients without OS time; count recorded)
    meta = pd.read_csv(RAW / "Basel_PatientMetadata.csv")
    tumour = meta[meta["diseasestatus"] == "tumor"].groupby("PID").first()
    p = pd.DataFrame(index=pd.Index(pat.index.astype(str), name="PID"))
    p["age"], p["grade"], p["pN"] = pat["age"].to_numpy(), pat["grade"].to_numpy(), pat["pN"].to_numpy()
    p["ER"], p["PR"], p["HER2"] = pat["ER"].to_numpy(), pat["PR"].to_numpy(), pat["HER2"].to_numpy()
    p["node_positive"] = np.where(p["pN"].isna(), np.nan, (p["pN"] > 0).astype(float))
    # Receptor-based surrogate subtypes (task definition): HER2+ ; ER+/HER2- ; triple negative.
    # ER-/PR+/HER2- fits none of them and is left missing.
    sub = np.select([p["HER2"] == 1, (p["ER"] == 1) & (p["HER2"] == 0),
                     (p["ER"] == 0) & (p["PR"] == 0) & (p["HER2"] == 0)], ["HER2+", "ER+/HER2-", "TN"], default="")
    p["surrogate"] = pd.Series(sub, index=p.index).replace("", np.nan)
    p["clinical_type_published"] = tumour["clinical_type"].reindex(pat.index).to_numpy()
    cores = cores.assign(PID=cores["PID"].astype(str))
    p.attrs["n_dropped_time"] = pat.attrs["n_dropped_time"]
    return p, cores


def main(cohorts: list[str]) -> None:
    PROC.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "cohort_counts.json"
    counts: dict = json.loads(path.read_text()) if path.exists() else {}
    if "metabric" in cohorts:
        counts["metabric"] = build_metabric()
    if "basel" in cohorts:
        counts["basel"] = build_basel()
    path.write_text(json.dumps(counts, indent=2, default=str))
    print(json.dumps({k: counts[k] for k in cohorts}, indent=2, default=str))


def build_metabric() -> dict:
    mc = metabric_cells()
    mp = metabric_patients(pd.Index(mc["PID"].unique()))
    mc.to_csv(PROC / "metabric_cells.csv.gz", index=False, compression="gzip", float_format="%.6g")
    mp.to_csv(OUT / "patients_metabric.csv")
    linked = mp["pam50_raw"].notna() | mp["sequenced"]
    return {
        "imc_patients": int(len(mp)), "imc_images": int(mc["image"].nunique()), "cells": int(len(mc)),
        "linked_to_cbioportal_genomics_or_subtype": int(linked.sum()),
        "sequenced": int(mp["sequenced"].sum()),
        "pam50_raw": mp["pam50_raw"].value_counts(dropna=False).rename(index=str).to_dict(),
        **{f"{g}_status": mp[g].map({1.0: "mutant", 0.0: "wild-type"}).fillna("not called").value_counts().to_dict()
           for g in GENES},
        "TP53_by_pam50": pd.crosstab(mp["pam50"].fillna("other/NA"), mp["TP53"].map({1.0: "mut", 0.0: "wt"})
                                     .fillna("NA")).to_dict(orient="index"),
        "PIK3CA_by_pam50": pd.crosstab(mp["pam50"].fillna("other/NA"), mp["PIK3CA"].map({1.0: "mut", 0.0: "wt"})
                                       .fillna("NA")).to_dict(orient="index"),
    }


def build_basel() -> dict:
    bp, cores = basel_patients()
    bc = basel_cells(cores)
    bc.to_csv(PROC / "basel_cells.csv.gz", index=False, compression="gzip", float_format="%.6g")
    bp = bp.loc[bp.index.isin(bc["PID"].unique())]
    bp.to_csv(OUT / "patients_basel.csv")
    return {
        "patients_with_cells": int(len(bp)), "patients_dropped_no_os_time_by_loader": int(bp.attrs.get("n_dropped_time", 0)),
        "images": int(bc["image"].nunique()), "cells_with_expression": int(len(bc)),
        "cells_with_labels_and_location": int(bc.attrs["n_location_cells"]),
        "surrogate": bp["surrogate"].value_counts(dropna=False).rename(index=str).to_dict(),
    }


if __name__ == "__main__":
    import sys

    main(sys.argv[1:] or ["metabric", "basel"])
