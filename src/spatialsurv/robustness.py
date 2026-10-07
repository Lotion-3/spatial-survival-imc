"""Robustness checks for the incremental value of the spatial block.

The patient bootstrap in eval.py holds the out-of-fold predictions fixed, so it does not
reflect variability from the CV partition or from refitting the models. Two checks:

1. CV-repeat stability: paired C-index deltas over many CV repeats (seeds).
2. Permutation null: shuffle the spatial block rows across patients (breaking any link to
   outcome while keeping clinical/composition intact), rerun the full nested CV, and
   compare the observed delta with the null distribution of deltas.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .cv import cindex, run_nested_cv
from .models import FeatureSet


def _oof(fs: FeatureSet, X: pd.DataFrame, y: np.ndarray, cfg: dict, seed: int) -> np.ndarray:
    return run_nested_cv(fs, X, y, cfg, seed)[0]


def _permuted(X: pd.DataFrame, cols: list[str], perm_seed: int) -> pd.DataFrame:
    Xp = X.copy()
    idx = np.random.default_rng(perm_seed).permutation(len(X))
    Xp[cols] = X[cols].to_numpy()[idx]
    return Xp


def cv_repeat_deltas(
    base: FeatureSet, others: list[FeatureSet], X, y, cfg, seeds: list[int], n_jobs: int = -1
) -> pd.DataFrame:
    """Per-seed C-index of each model and its delta vs. `base` (same folds per seed)."""
    fsets = [base] + others
    tasks = [(fs, s) for s in seeds for fs in fsets]
    oofs = Parallel(n_jobs=n_jobs)(delayed(_oof)(fs, X, y, cfg, s) for fs, s in tasks)
    c = {(fs.name, s): cindex(y, r) for (fs, s), r in zip(tasks, oofs)}
    rows = []
    for s in seeds:
        for fs in others:
            rows.append(dict(seed=s, model=fs.name, base=base.name, cindex=c[(fs.name, s)],
                             cindex_base=c[(base.name, s)], delta=c[(fs.name, s)] - c[(base.name, s)]))
    return pd.DataFrame(rows)


def permutation_null(
    base_oof: np.ndarray, fs: FeatureSet, X, y, cfg, seeds: list[int], n_perm: int, n_jobs: int = -1
) -> np.ndarray:
    """Null deltas C(fs with permuted spatial block) - C(base), shape (n_perm, n_seeds).

    Average over axis 1 for the null of the seed-averaged delta; use single entries to compare
    with single CV repeats. base_oof: (n_seeds, n) OOF risks of the base model for `seeds`
    (base is unaffected by the permutation, so it is not refit).
    """
    tasks = [(p, s) for p in range(n_perm) for s in seeds]
    oofs = Parallel(n_jobs=n_jobs)(
        delayed(_oof)(fs, _permuted(X, fs.spatial, 10_000 + p), y, cfg, s) for p, s in tasks
    )
    c_base = {s: cindex(y, r) for s, r in zip(seeds, base_oof)}
    d = np.zeros((n_perm, len(seeds)))
    for (p, s), r in zip(tasks, oofs):
        d[p, seeds.index(s)] = cindex(y, r) - c_base[s]
    return d
