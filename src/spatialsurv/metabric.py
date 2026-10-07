"""METABRIC IMC cohort (Danenberg et al. 2022) harmonised to the Basel analysis (part II).

All choices here follow ANALYSIS_PLAN.md (pre-registered) and its amendments.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .data import CLINICAL_COLS

HORIZON_MONTHS = 180.0  # administrative censoring, plan section A
MIN_CELLS = 500

# cellPhenotype -> Basel coarse class for NON-epithelial cells (epithelial cells -> "Tumor").
# ANALYSIS_PLAN.md Amendment 1 (made after viewing phenotype names only; no outcomes).
PHENOTYPE_CLASS: dict[str, str] = {
    "CD4^{+} T cells": "T",
    "CD8^{+} T cells": "T",
    "T_{Reg} & T_{Ex}": "T",
    "CD4^{+} T cells & APCs": "T",
    "CD57^{+}": "T",
    "B cells": "B",
    "CD38^{+} lymphocytes": "B",
    "Macrophages": "Macrophage",
    "Macrophages & granulocytes": "Macrophage",
    "Granulocytes": "Macrophage",
    "Endothelial": "Endothelial",
    "Fibroblasts": "Stroma",
    "Fibroblasts FSP1^{+}": "Stroma",
    "Myofibroblasts": "Stroma",
    "Myofibroblasts PDPN^{+}": "Stroma",
    "Ki67^{+}": "Stroma",
}

_PN_BANDS = [(0, 0), (1, 1), (4, 2), (10, 3)]  # positive-node count lower bounds -> pN


def pn_from_count(n: float) -> float:
    """AJCC pN bands from the number of positive lymph nodes: 0, 1-3, 4-9, >=10."""
    if pd.isna(n):
        return np.nan
    out = 0
    for lo, band in _PN_BANDS:
        if n >= lo:
            out = band
    return float(out)


def load_metabric_clinical(raw: Path) -> pd.DataFrame:
    clin = json.loads((raw / "metabric_clinical.json").read_text())
    df = pd.DataFrame.from_dict(clin, orient="index")
    num = lambda c: pd.to_numeric(df.get(c), errors="coerce")  # noqa: E731
    out = pd.DataFrame(index=df.index)
    vital = df.get("VITAL_STATUS")
    os_time = num("OS_MONTHS")
    died = vital.isin(["Died of Disease", "Died of Other Causes"])
    known = vital.isin(["Living", "Died of Disease", "Died of Other Causes"])
    out["time"] = np.minimum(os_time, HORIZON_MONTHS)
    after = os_time > HORIZON_MONTHS  # events after the horizon become censored at it
    out["event_os"] = died & ~after
    out["event_dss"] = vital.eq("Died of Disease") & ~after
    out["vital_known"] = known
    out["age"] = num("AGE_AT_DIAGNOSIS")
    out["tumor_size"] = num("TUMOR_SIZE")
    out["grade"] = num("GRADE")
    out["pN"] = num("LYMPH_NODES_EXAMINED_POSITIVE").map(pn_from_count)
    for src, dst in [("ER_STATUS", "ER"), ("PR_STATUS", "PR"), ("HER2_STATUS", "HER2")]:
        out[dst] = df.get(src).map({"Positive": 1.0, "Negative": 0.0})
    assert all(c in out for c in CLINICAL_COLS)
    return out


def load_metabric_cells(raw: Path) -> pd.DataFrame:
    """Invasive-tumour images only (is_tumour == 1); artefact cells (is_hotAggregate) removed."""
    cols = ["ImageNumber", "metabric_id", "cellPhenotype", "is_epithelial", "is_tumour", "is_hotAggregate",
            "Location_Center_X", "Location_Center_Y"]
    c = pd.read_csv(raw / "SingleCells.csv", usecols=cols)
    c = c[(c["is_tumour"] == 1) & (c["is_hotAggregate"] == 0)]
    c = c.rename(columns={"metabric_id": "PID", "Location_Center_X": "x", "Location_Center_Y": "y"})
    c["core"] = "MB_img" + c["ImageNumber"].astype(str)
    c["coarse"] = np.where(c["is_epithelial"] == 1, "Tumor", c["cellPhenotype"].map(PHENOTYPE_CLASS))
    unmapped = c.loc[c["coarse"].isna(), "cellPhenotype"].unique()
    if len(unmapped):
        raise ValueError(f"unmapped non-epithelial phenotypes: {sorted(unmapped)}")
    return c[["core", "PID", "x", "y", "cellPhenotype", "is_epithelial", "coarse"]].reset_index(drop=True)
