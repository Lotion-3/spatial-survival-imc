"""Nested cross-validation at the patient level.

Rows of X are patients (images were aggregated to patients beforehand), so every split
is a patient split. Outer: stratified k-fold on the event indicator. Inner: stratified
k-fold within the outer training set, used only to choose the penalty.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sksurv.metrics import concordance_index_censored

from .models import FeatureSet, alpha_grid, coefficients, fit_final, fit_path, predict_risk


def cindex(y: np.ndarray, risk: np.ndarray) -> float:
    return float(concordance_index_censored(y["event"], y["time"], risk)[0])


def patient_folds(event: np.ndarray, n_splits: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(skf.split(np.zeros(len(event)), event.astype(int)))


def select_alpha(fs: FeatureSet, X: pd.DataFrame, y: np.ndarray, cfg: dict, seed: int) -> tuple[float, float, int]:
    """Inner CV over a grid derived from (outer-)training data.

    Returns (alpha, mean inner C, number of (inner fold, alpha) fits that failed)."""
    grid = alpha_grid(fs, X, y, cfg)
    scores = np.full((cfg["inner_folds"], len(grid)), np.nan)
    for f, (tr, va) in enumerate(patient_folds(y["event"], cfg["inner_folds"], seed)):
        preds = fit_path(fs, X.iloc[tr], y[tr], grid, cfg)
        for j, p in enumerate(preds):
            if p is not None:
                r = p(X.iloc[va])
                if np.all(np.isfinite(r)):
                    scores[f, j] = cindex(y[va], r)
    mean = np.nanmean(scores, axis=0)
    j = int(np.nanargmax(mean))
    return float(grid[j]), float(mean[j]), int(np.isnan(scores).sum())


def run_nested_cv(
    fs: FeatureSet, X: pd.DataFrame, y: np.ndarray, cfg: dict, seed: int
) -> tuple[np.ndarray, list[dict], list[np.ndarray]]:
    """Returns (out-of-fold risk, per-fold records, per-fold coefficient vectors)."""
    X = X[fs.columns]
    oof = np.full(len(X), np.nan)
    records, coefs = [], []
    for k, (tr, te) in enumerate(patient_folds(y["event"], cfg["outer_folds"], seed)):
        alpha, inner_c, n_failed = select_alpha(fs, X.iloc[tr], y[tr], cfg, seed=seed * 1000 + k)
        model = fit_final(fs, X.iloc[tr], y[tr], alpha, cfg)
        oof[te] = predict_risk(model, X.iloc[te])
        coef = coefficients(model)
        coefs.append(coef)
        records.append(
            dict(
                model=fs.name, seed=seed, fold=k, n_train=len(tr), n_test=len(te),
                events_test=int(y["event"][te].sum()), alpha_selected=alpha, alpha_fitted=model.chosen_alpha_, inner_cindex=inner_c, inner_fits_failed=n_failed,
                test_cindex=cindex(y[te], oof[te]),
                n_nonzero_omics=int(np.sum(np.abs(coef[len(fs.clinical):]) > 1e-10)),
                test_ids=";".join(map(str, X.index[te])),
            )
        )
    return oof, records, coefs


def random_baseline(y: np.ndarray, seed: int) -> np.ndarray:
    """Sanity check: risk scores carrying no information -> C-index ~ 0.5."""
    return np.random.default_rng(seed).normal(size=len(y))


def run_models(
    feature_sets: list[FeatureSet], X: pd.DataFrame, y: np.ndarray, cfg: dict, seeds: list[int]
) -> tuple[dict[str, np.ndarray], pd.DataFrame, dict[str, list[np.ndarray]]]:
    """Nested CV for every feature set and seed, plus the random baseline.

    Returns oof[model] (n_seeds, n), per-fold records, coefficient vectors per model.
    The same outer folds (per seed) are used for every feature set.
    """
    oof, records, coefs = {}, [], {}
    for fs in feature_sets:
        runs = [run_nested_cv(fs, X, y, cfg, s) for s in seeds]
        oof[fs.name] = np.stack([r[0] for r in runs])
        records += [rec for r in runs for rec in r[1]]
        coefs[fs.name] = [c for r in runs for c in r[2]]
    oof["random baseline"] = np.stack([random_baseline(y, s) for s in seeds])
    return oof, pd.DataFrame(records), coefs
