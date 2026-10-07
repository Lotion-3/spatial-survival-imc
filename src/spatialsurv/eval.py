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


# --------------------------------------------------------------------------- part II, section C


def calibration_slope(y: np.ndarray, lp: np.ndarray) -> float:
    """Coefficient of a univariable Cox model of the outcome on the out-of-fold linear
    predictor. 1 = well calibrated in relative risk; < 1 = predictions too extreme (overfit)."""
    import warnings

    from sksurv.linear_model import CoxPHSurvivalAnalysis

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow").fit(np.asarray(lp, float).reshape(-1, 1), y)
    return float(m.coef_[0])


def td_auc(y: np.ndarray, risk: np.ndarray, times: list[float]) -> np.ndarray:
    """Cumulative/dynamic AUC at the given times (IPCW, censoring estimated from y itself)."""
    from sksurv.metrics import cumulative_dynamic_auc

    return cumulative_dynamic_auc(y, y, np.asarray(risk, float), np.asarray(times, float))[0]


def extra_metrics(y: np.ndarray, oof: dict[str, np.ndarray], times: list[float], n_boot: int, seed: int) -> pd.DataFrame:
    """Calibration slope and time-dependent AUC per model (mean over CV repeats), with
    percentile 95% CIs from patient bootstrap resamples (same resamples for every model)."""

    def stats(idx):
        out = {}
        for m, r in oof.items():
            yi = y[idx]
            out[(m, "calibration_slope")] = float(np.mean([calibration_slope(yi, rs[idx]) for rs in r]))
            aucs = np.mean([td_auc(yi, rs[idx], times) for rs in r], axis=0)
            for t, a in zip(times, aucs):
                out[(m, f"auc_{int(t)}m")] = float(a)
        return out

    def safe(idx):
        try:
            return stats(idx)
        except ValueError:  # resample without events before a time point
            return None

    from joblib import Parallel, delayed

    point = stats(np.arange(len(y)))
    rng = np.random.default_rng(seed)
    idxs = [rng.integers(0, len(y), len(y)) for _ in range(n_boot)]
    boots = Parallel(n_jobs=-1)(delayed(safe)(i) for i in idxs)
    B = pd.DataFrame([b for b in boots if b is not None])
    rows = []
    for (m, metric), v in point.items():
        rows.append(dict(model=m, metric=metric, value=v, ci_low=B[(m, metric)].quantile(0.025),
                         ci_high=B[(m, metric)].quantile(0.975), n_boot=len(B)))
    return pd.DataFrame(rows)
