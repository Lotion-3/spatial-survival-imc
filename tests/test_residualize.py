import numpy as np
import pandas as pd

from spatialsurv.data import survival_array
from spatialsurv.features import clr
from spatialsurv.models import FeatureSet, ResidualizeSpatial, make_model


def _toy(n=200, seed=0):
    rng = np.random.default_rng(seed)
    C = rng.poisson(lam=[300, 80, 20], size=(n, 3)).astype(float)
    Z = np.c_[clr(C), np.log(C.sum(axis=1) + 1)]
    arrangement = rng.normal(size=n)  # independent of composition by construction
    s_conf = Z @ np.array([2.0, -1.0, -1.0, 0.5])  # purely composition/size driven
    s_mixed = s_conf + arrangement
    return C, s_conf, s_mixed, arrangement


def test_removes_composition_and_keeps_arrangement():
    C, s_conf, s_mixed, arrangement = _toy()
    X = np.c_[s_conf, s_mixed, C]
    R = ResidualizeSpatial(n_spatial=2, pseudocount=0.5).fit(X).transform(X)
    assert np.abs(R[:, 0]).max() < 1e-6  # pure composition effect -> residual ~ 0
    # arrangement survives (not exactly 1: in a finite sample the regression also removes
    # arrangement's chance correlation with composition)
    assert np.corrcoef(R[:, 1], arrangement)[0, 1] > 0.97


def test_residualizer_fit_on_training_rows_only():
    C, _, s_mixed, _ = _toy()
    X = np.c_[s_mixed, C]
    tr = np.arange(150)
    a = ResidualizeSpatial(1).fit(X[tr])
    X2 = X.copy()
    X2[150:] = X2[150:] * 50 + 3  # perturb held-out rows only
    b = ResidualizeSpatial(1).fit(X2[tr])
    np.testing.assert_allclose(a.coef_, b.coef_)
    np.testing.assert_allclose(a.medians_, b.medians_)


def test_residualized_model_pipeline_runs_with_nans():
    C, _, s_mixed, _ = _toy(n=120)
    rng = np.random.default_rng(1)
    X = pd.DataFrame({"age": rng.normal(60, 10, 120), "c1": C[:, 0], "c2": C[:, 1], "c3": C[:, 2], "s1": s_mixed})
    X.loc[X.sample(frac=0.2, random_state=0).index, "s1"] = np.nan
    y = survival_array(rng.exponential(50, 120), rng.random(120) < 0.5)
    fs = FeatureSet("resid", ["age"], ["c1", "c2", "c3"], ["s1"], residualize_spatial=True)
    cfg = dict(clr_pseudocount=0.5, l1_ratio=0.5, n_alphas=5, alpha_min_ratio=0.05, coxnet_max_iter=100000)
    m = make_model(fs, cfg).fit(X, y)
    assert m.named_steps["pre"].transform(X).shape[1] == 1 + 3 + 1
