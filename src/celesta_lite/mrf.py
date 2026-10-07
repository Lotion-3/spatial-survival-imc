"""CELESTA cell-type assignment: anchor cells, then iterative spatial (Potts/MRF) refinement.

Mirrors the reference R implementation (plevritis-lab/CELESTA, ``AssignCells``), vectorised:

For each lineage round (types sharing one parent):

1. **Scores.** :math:`F_{ik}` for candidate cells (cells of the parent type, or all
   non-artifact cells in round 1).
2. **Anchors.** Cell i gets :math:`k^* = \\arg\\max_k F_{ik}` if every marker that :math:`k^*`
   expects high has :math:`EP \\ge` ``anchor_high`` and every marker expected low has
   :math:`EP \\le` ``anchor_low`` (reference defaults 0.7 / 0.9).
3. **Index cells (mean field).** For still-unassigned cells

   .. math::
      u_{ik} \\propto \\exp(F_{ik}) \\, \\exp\\Big(\\beta_{ik} \\sum_{j \\in N(i),\\, s_j = k} p_{jk}\\Big),
      \\qquad \\beta_{ik} = \\gamma \\, (1 - d_{ik}/h) \\text{ if } d_{ik} < h \\text{ else } 0,

   where :math:`N(i)` are the cell's ``n_neighbors`` nearest cells, :math:`d_{ik}` is the
   distance to the nearest cell already assigned type k, γ = ``scale_factor`` (5) and
   h = ``bandwidth``. Cells whose argmax type passes the *index* marker thresholds
   (defaults 0.5 / 1.0) are assigned. The signature is then moved halfway towards the
   mean EP of the cells assigned to each type, and the scores are recomputed.
4. Iterate until the share of round candidates that change assignment in an iteration
   falls below ``cell_change_threshold`` (1%), or ``max_iteration`` is reached.

Cells never assigned are labelled ``"Unknown"``. Cells flagged by the artifact filter
(all markers > 0.9 or all < 0.4) are left out.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from .preprocessing import Signature, artifact_mask, cell_type_scores

__all__ = ["CelestaConfig", "CelestaResult", "run_celesta", "UNKNOWN"]

UNKNOWN = "Unknown"


@dataclass(frozen=True)
class CelestaConfig:
    """Defaults are the reference implementation's (``CreateCelestaObject``/``AssignCells``)."""

    n_neighbors: int = 5
    bandwidth: float = 100.0  # in coordinate units
    scale_factor: float = 5.0
    max_iteration: int = 10
    cell_change_threshold: float = 0.01
    anchor_high: float | dict[str, float] = 0.7
    anchor_low: float | dict[str, float] = 0.9
    index_high: float | dict[str, float] = 0.5
    index_low: float | dict[str, float] = 1.0
    filter_artifacts: bool = True
    artifact_high: float = 0.9
    artifact_low: float = 0.4

    def __post_init__(self) -> None:
        if self.n_neighbors < 1:
            raise ValueError("n_neighbors must be >= 1")
        if self.bandwidth <= 0 or self.scale_factor < 0:
            raise ValueError("bandwidth must be > 0 and scale_factor >= 0")
        if self.max_iteration < 1 or not 0 <= self.cell_change_threshold <= 1:
            raise ValueError("invalid iteration settings")

    @staticmethod
    def per_type(value: float | dict[str, float], types: list[str]) -> np.ndarray:
        if isinstance(value, dict):
            return np.array([value.get(t, np.nan) for t in types], float)
        return np.full(len(types), float(value))


@dataclass
class CelestaResult:
    labels: pd.Series  # final (deepest) type per cell, or UNKNOWN
    is_anchor: pd.Series  # anchor in the deepest round that assigned the cell
    is_artifact: pd.Series
    marker_only: pd.Series  # argmax of initial scores per round, no thresholds, no spatial term
    rounds: list[dict] = field(default_factory=list)  # per-round diagnostics


def _assign(cand: np.ndarray, P: np.ndarray, ep: np.ndarray, prior0: np.ndarray,
            hi_thr: np.ndarray, lo_thr: np.ndarray, current: np.ndarray) -> np.ndarray:
    """Reference ``AssignCellTypes``: argmax type, accepted only if the marker thresholds pass."""
    out = current.copy()
    if cand.size == 0:
        return out
    k = P[cand].argmax(axis=1)
    high = prior0[k] == 1.0
    low = prior0[k] == 0.0
    E = ep[cand]
    ok = np.all(~high | (E >= hi_thr[k, None]), axis=1) & np.all(~low | (E <= lo_thr[k, None]), axis=1)
    out[cand[ok]] = k[ok]
    return out


def _beta(coords: np.ndarray, cand_unassigned: np.ndarray, assign: np.ndarray, K: int,
          cfg: CelestaConfig) -> np.ndarray:
    beta = np.zeros((cand_unassigned.size, K))
    for k in range(K):
        pts = np.flatnonzero(assign == k)
        if pts.size == 0 or cand_unassigned.size == 0:
            continue
        d, _ = KDTree(coords[pts]).query(coords[cand_unassigned], k=1, distance_upper_bound=cfg.bandwidth)
        inside = np.isfinite(d)
        beta[inside, k] = cfg.scale_factor * (1.0 - d[inside] / cfg.bandwidth)
    return beta


def run_celesta(ep: pd.DataFrame, coords: np.ndarray, signature: Signature,
                config: CelestaConfig = CelestaConfig()) -> CelestaResult:
    """Assign cell types to the cells of one image.

    Args:
        ep: cells x markers expression probabilities (from :func:`expression_probabilities`),
            containing every signature marker.
        coords: (n, 2) cell centroids, aligned with ``ep``.
        signature: prior knowledge with lineage rounds.
        config: thresholds and MRF parameters.
    """
    missing = [m for m in signature.markers if m not in ep.columns]
    if missing:
        raise KeyError(f"signature markers missing from EP: {missing}")
    xy = np.asarray(coords, float)
    if xy.shape != (len(ep), 2):
        raise ValueError("coords must be (n_cells, 2) and aligned with ep")
    E = ep[signature.markers].to_numpy(dtype=float)
    n = len(E)
    artifact = artifact_mask(ep[signature.markers], config.artifact_high, config.artifact_low) \
        if config.filter_artifacts else np.zeros(n, bool)
    k_eff = min(config.n_neighbors, n - 1)
    nbrs = KDTree(xy).query(xy, k=k_eff + 1)[1][:, 1:]

    label = np.full(n, UNKNOWN, dtype=object)
    marker_only = np.full(n, UNKNOWN, dtype=object)
    anchor = np.zeros(n, bool)
    diags = []
    for r in signature.rounds:
        types = signature.types_in_round(r)
        K = len(types)
        parent = signature.parent_of_round(r)
        cand = np.flatnonzero(~artifact & ((label == parent) if parent else np.ones(n, bool)))
        info = dict(round=r, parent=parent, types=types, candidates=int(cand.size))
        if cand.size == 0:
            diags.append(info | dict(anchors=0, iterations=0, change=[]))
            continue
        prior0 = signature.matrix(types)
        prior = prior0.copy()
        a_hi, a_lo = config.per_type(config.anchor_high, types), config.per_type(config.anchor_low, types)
        i_hi, i_lo = config.per_type(config.index_high, types), config.per_type(config.index_low, types)

        S = np.zeros((n, K))
        S[cand] = cell_type_scores(E[cand], prior)[1]
        marker_only[cand] = np.array(types, dtype=object)[S[cand].argmax(axis=1)]
        P = S.copy()
        assign = np.full(n, -1)
        assign = _assign(cand, P, E, prior0, a_hi, a_lo, assign)
        is_anchor_r = assign >= 0
        info["anchors"] = int(is_anchor_r.sum())
        changes: list[float] = []
        it = 1
        while info["anchors"] > 0 and it < config.max_iteration:
            it += 1
            un = cand[assign[cand] < 0]
            if un.size == 0:
                break
            beta = _beta(xy, un, assign, K, config)
            nb = nbrs[un]  # (u, N)
            nb_type = assign[nb]
            sums = np.zeros((un.size, K))
            for k in range(K):
                sums[:, k] = np.where(nb_type == k, P[nb, k], 0.0).sum(axis=1)
            logu = S[un] + beta * sums
            u = np.exp(logu - logu.max(axis=1, keepdims=True))
            P[un] = u / u.sum(axis=1, keepdims=True)

            old = assign.copy()
            assign = _assign(un, P, E, prior0, i_hi, i_lo, assign)
            frac = float(np.sum(old != assign) / cand.size)
            changes.append(frac)
            if frac < config.cell_change_threshold:
                break
            # move the signature halfway towards the assigned cells' mean EP (NaN stays NaN)
            for k in range(K):
                members = np.flatnonzero(assign == k)
                if members.size:
                    prior[k] = (E[members].mean(axis=0) + prior0[k]) / 2.0
            un = cand[assign[cand] < 0]
            if un.size:
                S[un] = cell_type_scores(E[un], prior)[1]
        info.update(iterations=it - 1, change=changes, assigned=int((assign >= 0).sum()),
                    counts={t: int((assign == k).sum()) for k, t in enumerate(types)})
        diags.append(info)
        done = assign >= 0
        label[done] = np.array(types, dtype=object)[assign[done]]
        anchor[done] = is_anchor_r[done]

    idx = ep.index
    return CelestaResult(labels=pd.Series(label, index=idx), is_anchor=pd.Series(anchor, index=idx),
                         is_artifact=pd.Series(artifact, index=idx), marker_only=pd.Series(marker_only, index=idx),
                         rounds=diags)
