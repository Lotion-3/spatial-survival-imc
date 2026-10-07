import numpy as np
import pandas as pd
import pytest

import spatialsurv.cv as cvmod
from spatialsurv.data import survival_array
from spatialsurv.features import aggregate_to_patient, clr, composition_counts, SPATIAL_COLS
from spatialsurv.models import CLRTransformer, FeatureSet, make_model

CFG = dict(outer_folds=5, inner_folds=3, clr_pseudocount=0.5, l1_ratio=0.5, n_alphas=8,
           alpha_min_ratio=0.05, coxnet_max_iter=100000, ridge_log10_min=-1, ridge_log10_max=2)
CLIN, COMP, SPAT = ["age", "er"], ["c1", "c2", "c3"], ["s1", "s2"]


def toy_data(n=120, seed=0):
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(
        {
            "age": rng.normal(60, 10, n),
            "er": rng.integers(0, 2, n).astype(float),
            "c1": rng.poisson(100, n), "c2": rng.poisson(30, n), "c3": rng.poisson(5, n),
            "s1": rng.normal(size=n), "s2": rng.normal(size=n),
        },
        index=[f"P{i}" for i in range(n)],
    )
    X.loc[X.sample(frac=0.15, random_state=1).index, "age"] = np.nan
    X.loc[X.sample(frac=0.15, random_state=2).index, "s1"] = np.nan
    y = survival_array(rng.exponential(50, n), rng.random(n) < 0.4)
    return X, y


def test_preprocessing_params_depend_only_on_training_rows():
    X, y = toy_data()
    fs = FeatureSet("full", CLIN, COMP, SPAT)
    tr, te = np.arange(90), np.arange(90, 120)
    m1 = make_model(fs, CFG).fit(X.iloc[tr], y[tr])

    X2 = X.copy()  # wildly perturb the held-out rows only
    X2.iloc[te, :] = X2.iloc[te, :] * 1000 + 7
    m2 = make_model(fs, CFG).fit(X2.iloc[tr], y[tr])

    pre1, pre2 = m1.named_steps["pre"], m2.named_steps["pre"]
    for blk in ["clin", "comp", "spat"]:
        s1, s2 = pre1.named_transformers_[blk], pre2.named_transformers_[blk]
        np.testing.assert_allclose(s1[-1].mean_, s2[-1].mean_)
        np.testing.assert_allclose(s1[-1].scale_, s2[-1].scale_)
    # imputer medians equal the training-row medians, not the full-data medians
    imp = pre1.named_transformers_["clin"][0]
    np.testing.assert_allclose(imp.statistics_, X.iloc[tr][CLIN].median().to_numpy())
    assert not np.allclose(imp.statistics_[0], X2[CLIN].median().iloc[0])
    # scaler of the composition block = mean of training-row CLR values
    np.testing.assert_allclose(pre1.named_transformers_["comp"][-1].mean_,
                               clr(X.iloc[tr][COMP].to_numpy()).mean(axis=0))


def test_clr_is_rowwise():
    X, _ = toy_data()
    a = CLRTransformer().fit(X[COMP]).transform(X[COMP].iloc[:5])
    other = X[COMP].copy()
    other.iloc[5:] = 1
    b = CLRTransformer().fit(other).transform(other.iloc[:5])
    np.testing.assert_allclose(a, b)
    np.testing.assert_allclose(a.sum(axis=1), 0, atol=1e-12)


@pytest.mark.parametrize("fs", [FeatureSet("clin", CLIN), FeatureSet("full", CLIN, COMP, SPAT)])
def test_nested_cv_never_fits_on_test_patients(monkeypatch, fs):
    """Spy on every fit inside run_nested_cv: no fit (grid derivation, inner CV, final
    model) may see a patient from the current outer test fold."""
    X, y = toy_data()
    seen = []
    for name in ["alpha_grid", "fit_path", "fit_final"]:
        orig = getattr(cvmod, name)

        def spy(fs_, Xf, *a, _orig=orig, **k):
            seen.append(set(Xf.index))
            return _orig(fs_, Xf, *a, **k)

        monkeypatch.setattr(cvmod, name, spy)

    oof, records, _ = cvmod.run_nested_cv(fs, X, y, CFG, seed=0)
    assert np.isfinite(oof).all()
    test_sets = [set(r["test_ids"].split(";")) for r in records]
    # outer test folds partition the patients
    assert sum(len(t) for t in test_sets) == len(X)
    assert set().union(*test_sets) == set(X.index)
    # fits happen in fold order; each fold's fits (1 grid + inner folds + 1 final) avoid its test set
    per_fold = 1 + CFG["inner_folds"] + 1
    assert len(seen) == per_fold * CFG["outer_folds"]
    for k, test in enumerate(test_sets):
        for fitted_on in seen[k * per_fold:(k + 1) * per_fold]:
            assert fitted_on.isdisjoint(test)


def test_patient_level_folds_disjoint_and_stratified():
    rng = np.random.default_rng(0)
    event = rng.random(280) < 0.28
    folds = cvmod.patient_folds(event, 5, seed=3)
    all_test = np.concatenate([te for _, te in folds])
    assert sorted(all_test) == list(range(280))
    for tr, te in folds:
        assert set(tr).isdisjoint(te)
        assert abs(event[te].mean() - event.mean()) < 0.03


def test_multi_image_patients_collapse_to_one_row():
    """Images are aggregated to patients before CV, so a patient's two cores can never be
    split across train and test."""
    cores = pd.DataFrame({"core": ["a1", "a2", "b1"], "PID": [1, 1, 2]})
    img = pd.DataFrame({c: [1.0, 3.0, 5.0] for c in SPATIAL_COLS}, index=["a1", "a2", "b1"])
    img["n_cells"] = [100, 300, 50]
    pat = aggregate_to_patient(img, cores)
    assert list(pat.index) == [1, 2]
    assert pat.loc[1, SPATIAL_COLS[0]] == pytest.approx((1 * 100 + 3 * 300) / 400)  # cell-weighted
    cells = pd.DataFrame({"core": ["a1", "a2", "a2", "b1"], "cluster": [1, 1, 2, 2]})
    comp = composition_counts(cells, cores, [1, 2, 3])
    assert comp.loc[1].tolist() == [2, 1, 0] and comp.loc[2].tolist() == [0, 1, 0]
