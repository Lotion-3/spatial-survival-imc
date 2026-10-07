"""Load the Basel IMC cohort and build patient-level clinical/survival tables.

Column names below were taken from inspecting the raw files (see README, "Data").
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

META = "Basel_PatientMetadata.csv"
LABELS = "singlecell_cluster_labels/Cluster_labels/Basel_metaclusters.csv"
ANNOT = "singlecell_cluster_labels/Cluster_labels/Metacluster_annotations.csv"
LOCS = "singlecell_locations/Basel_SC_locations.csv"

# Coarse cell classes used by the spatial features (metacluster ids from Metacluster_annotations.csv).
# 2 ("T & B cells") is grouped with B cells as a lymphoid-aggregate class.
COARSE = {
    "B": [1, 2],
    "T": [3, 5],
    "Macrophage": [4, 6],
    "Endothelial": [7],
    "Stroma": [8, 9, 10, 11, 12, 13],
    "Tumor": list(range(14, 28)),
}
IMMUNE = ("B", "T", "Macrophage")

CLINICAL_COLS = ["age", "tumor_size", "grade_2", "grade_3", "pN", "pM", "ER", "PR", "HER2"]

_PN_MAP = {"0": 0, "0sl": 0, "0sn": 0, "1": 1, "1a": 1, "1mi": 1, "2": 2, "2a": 2, "3": 3, "3a": 3, "3b": 3}
_STATUS_MAP = {"positive": 1.0, "negative": 0.0}


def load_annotations(raw: Path) -> pd.DataFrame:
    ann = pd.read_csv(raw / ANNOT, sep=None, engine="python")
    ann.columns = [c.strip() for c in ann.columns]
    ann = ann.rename(columns={"Metacluster": "cluster", "Cell type": "cell_type", "Class": "cell_class"})
    ann["cell_type"] = ann["cell_type"].str.strip()
    # "Macrohage" is a typo in the source file
    ann["cell_type"] = ann["cell_type"].replace({"Macrohage": "Macrophage"})
    ann["label"] = ann["cluster"].astype(str).str.zfill(2) + "_" + ann["cell_type"]
    coarse = {c: name for name, ids in COARSE.items() for c in ids}
    ann["coarse"] = ann["cluster"].map(coarse)
    return ann


def load_cells(raw: Path, cores: list[str] | None = None) -> pd.DataFrame:
    """One row per labelled cell: core, x, y (pixels = um for IMC), metacluster, coarse class."""
    loc = pd.read_csv(raw / LOCS, usecols=["core", "id", "Location_Center_X", "Location_Center_Y"])
    lab = pd.read_csv(raw / LABELS)
    cells = loc.merge(lab, on="id", how="inner")  # drops ~1.3% unlabelled cells (mostly liver controls)
    if cores is not None:
        cells = cells[cells["core"].isin(cores)]
    ann = load_annotations(raw)
    cells = cells.merge(ann[["cluster", "coarse"]], on="cluster", how="left")
    cells = cells.rename(columns={"Location_Center_X": "x", "Location_Center_Y": "y"})
    return cells[["core", "id", "x", "y", "cluster", "coarse"]].reset_index(drop=True)


def load_patients(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (patients, cores).

    patients: one row per PID with survival targets and raw clinical features (NaN = missing).
    cores: tumor cores (core -> PID) used to build image features.
    Filtering: tumor cores only; drop patients with missing or non-positive OS time.
    """
    meta = pd.read_csv(raw / META)
    tumor = meta[meta["diseasestatus"] == "tumor"].copy()
    pat = tumor.groupby("PID").first()  # clinical fields are constant within PID (checked in tests)

    out = pd.DataFrame(index=pat.index)
    out["time"] = pd.to_numeric(pat["OSmonth"], errors="coerce")
    status = pat["Patientstatus"].astype(str)
    out["event_os"] = status.str.startswith("death")
    out["event_dss"] = status.eq("death by primary disease")
    out["age"] = pd.to_numeric(pat["age"], errors="coerce")
    out["tumor_size"] = pd.to_numeric(pat["tumor_size"], errors="coerce")
    grade = pd.to_numeric(pat["grade"], errors="coerce")
    out["grade_2"] = (grade == 2).astype(float).where(grade.notna())
    out["grade_3"] = (grade == 3).astype(float).where(grade.notna())
    out["pN"] = pat["PTNM_N"].astype(str).str.strip().map(_PN_MAP)  # 'x'/'X' (not assessed) -> NaN
    out["pM"] = pd.to_numeric(pat["PTNM_M"], errors="coerce")
    for src, dst in [("ERStatus", "ER"), ("PRStatus", "PR"), ("HER2Status", "HER2")]:
        out[dst] = pat[src].map(_STATUS_MAP)

    n_before = len(out)
    keep = out["time"].notna() & (out["time"] > 0)
    out = out[keep]
    cores = tumor.loc[tumor["PID"].isin(out.index), ["core", "PID"]].reset_index(drop=True)
    out.attrs["n_tumor_patients"] = n_before
    out.attrs["n_dropped_time"] = int((~keep).sum())
    return out, cores


def survival_array(time: np.ndarray, event: np.ndarray) -> np.ndarray:
    return np.array(list(zip(event.astype(bool), time.astype(float))), dtype=[("event", "?"), ("time", "<f8")])
