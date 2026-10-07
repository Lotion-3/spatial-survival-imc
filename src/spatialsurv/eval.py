"""Pooled out-of-fold C-index with patient-bootstrap CIs and paired differences."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cv import cindex


def seed_mean_cindex(y: np.ndarray, oof: dict[str, np.ndarray], idx: np.ndarray | None = None) -> dict[str, float]:
    """For each model, C-index of pooled OOF predictions, averaged over CV repeats (seeds).

    oof[model] has shape (n_seeds, n_patients).
    """
    out = {}
    for m, r in oof.items():
        if idx is None:
            out[m] = float(np.mean([cindex(y, rs) for rs in r]))
        else:
            out[m] = float(np.mean([cindex(y[idx], rs[idx]) for rs in r]))
    return out


def bootstrap(y: np.ndarray, oof: dict[str, np.ndarray], n_boot: int, seed: int) -> pd.DataFrame:
    """Resample patients with replacement; the SAME resamples are used for every model,
    so differences between models are paired. Returns (n_boot, n_models)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    rows = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if y["event"][idx].sum() == 0:
            continue
        rows.append(seed_mean_cindex(y, oof, idx))
    return pd.DataFrame(rows)


def summarize(point: dict[str, float], boot: pd.DataFrame, comparisons: list[tuple[str, str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Percentile 95% CIs for each model and for paired differences (a minus b)."""
    models = pd.DataFrame(
        [
            dict(model=m, cindex=point[m], ci_low=boot[m].quantile(0.025), ci_high=boot[m].quantile(0.975))
            for m in point
        ]
    )
    diffs = []
    for a, b in comparisons:
        d = boot[a] - boot[b]
        diffs.append(
            dict(
                comparison=f"{a} minus {b}",
                delta=point[a] - point[b],
                ci_low=d.quantile(0.025),
                ci_high=d.quantile(0.975),
                frac_boot_gt0=float((d > 0).mean()),
            )
        )
    return models, pd.DataFrame(diffs)
