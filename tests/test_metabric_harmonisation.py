import numpy as np

from spatialsurv.metabric import pn_from_count


def test_pn_bands_match_ajcc_counts():
    # AJCC: pN0 = 0 nodes, pN1 = 1-3, pN2 = 4-9, pN3 = >=10
    got = [pn_from_count(n) for n in [0, 1, 3, 4, 9, 10, 25]]
    assert got == [0, 1, 1, 2, 2, 3, 3]
    assert np.isnan(pn_from_count(np.nan))
