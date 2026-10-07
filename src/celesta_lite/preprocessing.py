"""Marker expression probabilities and cell-type scores (CELESTA stage 1).

Follows the reference CELESTA R implementation (plevritis-lab/CELESTA,
``CalcMarkerActivationProbability``, ``CalculateScores``):

1. Expression is arcsinh-transformed: :math:`x' = \\operatorname{arcsinh}(x / c)`.
2. Per marker, a two-component 1-D Gaussian mixture with **equal mixing proportions**
   and free variances is fit (a share of exact zeros is dropped first, as in the
   reference). :math:`x_c` is the point where the two weighted component densities
   are equal.
3. :math:`EP = \\sigma(x' - x_c)`, then min-max rescaled to [0, 1] across cells.
4. Cell-type scores :math:`Z_{ik} = 1 - \\frac{1}{|M_k|}\\sum_{m \\in M_k}(EP_{im} - SP_{km})^2`
   over the markers informative for type k, normalised over the candidate types:
   :math:`F_{ik} = Z_{ik} / \\sum_t Z_{it}`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = [
    "MarkerModel",
    "Signature",
    "arcsinh_transform",
    "artifact_mask",
    "cell_type_scores",
    "expression_probabilities",
    "fit_marker_models",
]


def arcsinh_transform(expression: pd.DataFrame, cofactor: float = 5.0) -> pd.DataFrame:
    """:math:`\\operatorname{arcsinh}(x / c)`, the standard variance stabiliser for cytometry counts."""
    if cofactor <= 0:
        raise ValueError("cofactor must be positive")
    if (expression < 0).any().any():
        raise ValueError("expression contains negative values")
    return np.arcsinh(expression / cofactor)


# --------------------------------------------------------------------------- signature


@dataclass(frozen=True)
class Signature:
    """Prior-knowledge matrix with lineage rounds.

    ``table`` rows are cell types. Required columns: ``round`` (int >= 1) and ``parent``
    (the cell type whose cells are refined in this round; empty/None for round 1).
    The remaining columns are markers with entries 1 (high), 0 (low) or NaN (ignored).
    All types in one round share one parent, as in the reference implementation.
    """

    table: pd.DataFrame

    def __post_init__(self) -> None:
        t = self.table
        for col in ("round", "parent"):
            if col not in t.columns:
                raise ValueError(f"signature needs a {col!r} column")
        if t.index.duplicated().any():
            raise ValueError("duplicated cell types")
        vals = t[self.markers].to_numpy(dtype=float)
        if not (np.isnan(vals) | np.isin(vals, [0.0, 1.0])).all():
            raise ValueError("marker entries must be 0, 1 or NaN")
        if np.isnan(vals).all(axis=1).any():
            raise ValueError("every cell type needs at least one informative marker")
        for r, g in t.groupby("round"):
            parents = set(g["parent"].fillna(""))
            if len(parents) != 1:
                raise ValueError(f"round {r} mixes parents {parents}")
            p = parents.pop()
            if r == 1 and p:
                raise ValueError("round-1 types must have no parent")
            if r > 1 and (p not in t.index or t.loc[p, "round"] >= r):
                raise ValueError(f"round {r} parent {p!r} must be a type from an earlier round")

    @property
    def markers(self) -> list[str]:
        return [c for c in self.table.columns if c not in ("round", "parent")]

    @property
    def rounds(self) -> list[int]:
        return sorted(int(r) for r in self.table["round"].unique())

    def types_in_round(self, r: int) -> list[str]:
        return list(self.table.index[self.table["round"] == r])

    def parent_of_round(self, r: int) -> str | None:
        p = self.table.loc[self.table["round"] == r, "parent"].iloc[0]
        return None if (p is None or p == "" or (isinstance(p, float) and np.isnan(p))) else str(p)

    def matrix(self, types: list[str]) -> np.ndarray:
        return self.table.loc[types, self.markers].to_numpy(dtype=float)


# --------------------------------------------------------------------------- per-marker GMM


@dataclass(frozen=True)
class MarkerModel:
    marker: str
    means: tuple[float, float]  # (low, high)
    variances: tuple[float, float]
    critical_point: float


def _drop_zeros(x: np.ndarray) -> np.ndarray:
    """Reference behaviour: remove part of the exact zeros before fitting (FitGmmModel)."""
    zero = np.flatnonzero(x == 0)
    frac = zero.size / x.size
    if frac < 0.1:
        return x
    if frac < 0.2:
        n_rm = int(np.floor(x.size * frac))
    elif frac < 0.5:
        n_rm = int(np.floor(x.size * (frac - 0.05)))
    elif frac <= 0.9:
        n_rm = int(np.ceil(x.size * (frac - 0.02)))
    else:
        n_rm = zero.size
    return np.delete(x, zero[:n_rm])


def _equal_weight_gmm(x: np.ndarray, n_iter: int = 300, tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """EM for a 2-component 1-D Gaussian mixture with mixing proportions fixed at 1/2."""
    mu = np.quantile(x, [0.25, 0.75]).astype(float)
    if mu[0] == mu[1]:
        mu = np.array([x.min(), x.max()], float)
    var = np.full(2, max(x.var() / 4.0, 1e-6))
    prev = -np.inf
    for _ in range(n_iter):
        logp = -0.5 * (np.log(2 * np.pi * var)[None, :] + (x[:, None] - mu[None, :]) ** 2 / var[None, :])
        m = logp.max(axis=1, keepdims=True)
        ll = float(np.sum(m[:, 0] + np.log(np.exp(logp - m).sum(axis=1))))
        r = np.exp(logp - m)
        r /= r.sum(axis=1, keepdims=True)
        nk = r.sum(axis=0) + 1e-12
        mu = (r * x[:, None]).sum(axis=0) / nk
        var = np.maximum((r * (x[:, None] - mu[None, :]) ** 2).sum(axis=0) / nk, 1e-6)
        if ll - prev < tol * abs(ll):
            break
        prev = ll
    order = np.argsort(mu)
    return mu[order], var[order]


def _crossing_point(mu: np.ndarray, var: np.ndarray) -> float:
    """Where the two (equal-weight) component densities are equal, between the means.

    Same quadratic as the reference ``BuildSigmoidFunction`` (weights cancel)."""
    (m1, m2), (v1, v2) = mu, var  # m1 = low component
    a = 0.5 / v2 - 0.5 / v1
    b = m1 / v1 - m2 / v2
    c = 0.5 * (-(m1**2) / v1 + m2**2 / v2) + 0.5 * np.log(v2 / v1)
    if abs(a) < 1e-12:  # equal variances: linear equation
        return float(-c / b) if b != 0 else float(0.5 * (m1 + m2))
    disc = b * b - 4 * a * c
    if disc < 0:
        return float(0.5 * (m1 + m2))
    roots = np.array([(-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)])
    inside = roots[(roots >= m1) & (roots <= m2)]
    return float(inside[0]) if inside.size else float(roots[np.argmin(np.abs(roots - 0.5 * (m1 + m2)))])


def fit_marker_models(expression: pd.DataFrame) -> dict[str, MarkerModel]:
    """Fit the equal-weight GMM and decision point for every marker column (transformed data)."""
    if expression.isna().any().any():
        raise ValueError("expression contains NaNs")
    out: dict[str, MarkerModel] = {}
    for m in expression.columns:
        x = _drop_zeros(expression[m].to_numpy(dtype=float))
        if x.size < 10 or np.unique(x).size < 2:
            raise ValueError(f"marker {m!r}: too few distinct values to fit a mixture")
        mu, var = _equal_weight_gmm(x)
        out[str(m)] = MarkerModel(str(m), (float(mu[0]), float(mu[1])), (float(var[0]), float(var[1])),
                                  _crossing_point(mu, var))
    return out


def expression_probabilities(expression: pd.DataFrame, models: dict[str, MarkerModel],
                             method: str = "reference") -> pd.DataFrame:
    """Per-marker expression probabilities.

    ``method="reference"``: :math:`\\sigma(x - x_c)` with slope 1, min-max rescaled to [0, 1]
    across cells, exactly as the reference implementation. The fixed slope assumes
    transformed intensities spanning several units (e.g. CODEX with arcsinh(x/10)).
    On dim channels whose whole range is a fraction of a unit, the sigmoid is nearly flat
    and EP is driven by the min-max of outliers.

    ``method="posterior"``: the mixture posterior :math:`P(\\text{high} \\mid x)` with equal
    priors, i.e. a sigmoid of the components' log-likelihood ratio. It equals 0.5 at
    :math:`x_c` and is invariant to the intensity scale. With unequal variances the broader
    component also wins in the opposite tail, so the raw posterior is made monotone around
    :math:`x_c`: a running maximum above it and a running minimum below it.
    """
    missing = [m for m in expression.columns if m not in models]
    if missing:
        raise KeyError(f"no fitted model for markers {missing}")
    x = expression.to_numpy(dtype=float)
    if method == "reference":
        xc = np.array([models[m].critical_point for m in expression.columns])
        y = 1.0 / (1.0 + np.exp(-(x - xc)))
        lo, hi = y.min(axis=0), y.max(axis=0)
        ep = np.divide(y - lo, hi - lo, out=np.zeros_like(y), where=hi > lo)
    elif method == "posterior":
        ep = np.empty_like(x)
        for j, m in enumerate(expression.columns):
            (m1, m2), (v1, v2) = models[m].means, models[m].variances
            xj = x[:, j]
            llr = (-0.5 * np.log(v2) - (xj - m2) ** 2 / (2 * v2)) - (-0.5 * np.log(v1) - (xj - m1) ** 2 / (2 * v1))
            p = 1.0 / (1.0 + np.exp(-np.clip(llr, -50, 50)))
            # Enforce monotonicity, anchored at the crossing point: a broader component also
            # wins in the *opposite* tail (e.g. a wide 'high' component dominates far below
            # the low mean). Running max to the right of x_c, running min to the left.
            order = np.argsort(xj, kind="stable")
            ps, xs = p[order], xj[order]
            right = xs >= models[m].critical_point
            ps[right] = np.maximum.accumulate(np.maximum(ps[right], 0.5))
            ps[~right] = np.minimum.accumulate(np.minimum(ps[~right], 0.5)[::-1])[::-1]
            p[order] = ps
            ep[:, j] = p
    else:
        raise ValueError(f"unknown method {method!r}")
    return pd.DataFrame(ep, index=expression.index, columns=expression.columns)


def artifact_mask(ep: pd.DataFrame, high: float = 0.9, low: float = 0.4) -> np.ndarray:
    """Reference ``FilterArtifactCells``: cells with *all* markers > high or *all* < low."""
    v = ep.to_numpy(dtype=float)
    return (v > high).all(axis=1) | (v < low).all(axis=1)


def cell_type_scores(ep: np.ndarray, prior: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Raw scores Z and normalised F for cells (rows of ``ep``) against ``prior`` (types x markers).

    NaN entries of ``prior`` are ignored for that type.
    """
    informative = ~np.isnan(prior)
    sq = (ep[:, None, :] - np.nan_to_num(prior)[None, :, :]) ** 2
    sq = np.where(informative[None], sq, 0.0)
    Z = 1.0 - sq.sum(axis=2) / informative.sum(axis=1)[None, :]
    s = Z.sum(axis=1, keepdims=True)
    F = np.divide(Z, s, out=np.full_like(Z, 1.0 / Z.shape[1]), where=s > 0)
    return Z, F
