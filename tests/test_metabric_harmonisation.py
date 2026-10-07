import json

import numpy as np
import pandas as pd
import pytest

from spatialsurv.data import COARSE
from spatialsurv.metabric import HORIZON_MONTHS, PHENOTYPE_CLASS, load_metabric_cells, load_metabric_clinical, pn_from_count


def test_pn_bands_match_ajcc_counts():
    # AJCC: pN0 = 0 nodes, pN1 = 1-3, pN2 = 4-9, pN3 = >=10
    got = [pn_from_count(n) for n in [0, 1, 3, 4, 9, 10, 25]]
    assert got == [0, 1, 1, 2, 2, 3, 3]
    assert np.isnan(pn_from_count(np.nan))


def test_clinical_censoring_and_endpoints(tmp_path):
    base = {"AGE_AT_DIAGNOSIS": "60", "TUMOR_SIZE": "20", "GRADE": "2", "LYMPH_NODES_EXAMINED_POSITIVE": "5",
            "ER_STATUS": "Positive", "PR_STATUS": "Negative", "HER2_STATUS": "Negative"}
    clin = {
        "early_cancer_death": base | {"OS_MONTHS": "50", "VITAL_STATUS": "Died of Disease"},
        "early_other_death": base | {"OS_MONTHS": "60", "VITAL_STATUS": "Died of Other Causes"},
        "late_cancer_death": base | {"OS_MONTHS": "200", "VITAL_STATUS": "Died of Disease"},
        "alive": base | {"OS_MONTHS": "100", "VITAL_STATUS": "Living"},
        "unknown": base | {"OS_MONTHS": "70"},
    }
    (tmp_path / "metabric_clinical.json").write_text(json.dumps(clin))
    d = load_metabric_clinical(tmp_path)
    assert d.loc["late_cancer_death", "time"] == HORIZON_MONTHS  # deaths after 15 years are censored at 15 years
    assert not d.loc["late_cancer_death", "event_os"] and not d.loc["late_cancer_death", "event_dss"]
    assert d.loc["early_cancer_death", "event_os"] and d.loc["early_cancer_death", "event_dss"]
    assert d.loc["early_other_death", "event_os"] and not d.loc["early_other_death", "event_dss"]  # DSS censors other causes
    assert not d.loc["alive", "event_os"]
    assert not d.loc["unknown", "vital_known"] and d.loc["alive", "vital_known"]
    assert d.loc["alive", "pN"] == 2 and d.loc["alive", "ER"] == 1 and d.loc["alive", "PR"] == 0


def test_phenotype_mapping_targets_basel_classes():
    assert set(PHENOTYPE_CLASS.values()) == set(COARSE) - {"Tumor"}  # every non-tumour Basel class is reachable


def _cells(phenotypes, epithelial, tumour=1, hot=0):
    n = len(phenotypes)
    return pd.DataFrame({"ImageNumber": [7] * n, "metabric_id": ["MB-1"] * n, "cellPhenotype": phenotypes,
                         "is_epithelial": epithelial, "is_tumour": tumour, "is_hotAggregate": hot,
                         "Location_Center_X": np.arange(n, dtype=float), "Location_Center_Y": 0.0})


def test_cell_filters_and_classes(tmp_path):
    c = _cells(["CK^{+} CXCL12^{+}", "CD8^{+} T cells", "Fibroblasts", "B cells", "Macrophages"],
               epithelial=[1, 0, 0, 0, 0], hot=[0, 0, 0, 0, 1])
    c = pd.concat([c, _cells(["CD8^{+} T cells"], [0], tumour=0)])  # a cell from a non-tumour image
    c.to_csv(tmp_path / "SingleCells.csv", index=False)
    out = load_metabric_cells(tmp_path)
    assert out["coarse"].tolist() == ["Tumor", "T", "Stroma", "B"]  # hot-aggregate and non-tumour-image cells dropped
    assert (out["core"] == "MB_img7").all()


def test_unmapped_phenotype_fails_loudly(tmp_path):
    _cells(["Something new"], [0]).to_csv(tmp_path / "SingleCells.csv", index=False)
    with pytest.raises(ValueError, match="unmapped"):
        load_metabric_cells(tmp_path)
