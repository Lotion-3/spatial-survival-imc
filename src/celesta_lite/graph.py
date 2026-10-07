"""Spatial neighbourhood graphs over cell centroids.

Two builders are provided:

* :func:`knn_graph`: each cell linked to its k nearest neighbours (CELESTA recommends
  5-10), optionally symmetrised. Robust to density differences and never leaves a cell
  isolated.
* :func:`delaunay_graph`: Delaunay triangulation, optionally with long edges pruned.
  Parameter-free, and adapts to local geometry.

Both return a :class:`SpatialGraph` that stores a sparse adjacency matrix with edge
distances, which the MRF engine uses for neighbour lookups and distance weighting.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.spatial import Delaunay, KDTree

__all__ = ["SpatialGraph", "delaunay_graph", "knn_graph", "neighbourhood_composition"]


@dataclass(frozen=True)
class SpatialGraph:
    """Sparse spatial graph. ``adjacency[i, j]`` holds the Euclidean distance of edge i-j."""

    coords: np.ndarray  # (n, 2)
    adjacency: sp.csr_matrix  # (n, n), distances on edges

    @property
    def n_cells(self) -> int:
        return int(self.coords.shape[0])

    def neighbours(self, i: int) -> np.ndarray:
        """Indices of the neighbours of cell ``i``."""
        a = self.adjacency
        return a.indices[a.indptr[i] : a.indptr[i + 1]]

    def degree(self) -> np.ndarray:
        return np.diff(self.adjacency.indptr)

    def weights(self, kind: str = "binary", bandwidth: float | None = None) -> sp.csr_matrix:
        """Edge weights derived from distances.

        Args:
            kind: ``"binary"`` (all edges 1), ``"inverse"`` (1/d) or ``"gaussian"``
                (exp(-d^2 / (2 h^2))).
            bandwidth: h for the Gaussian kernel (required for ``"gaussian"``).
        """
        w = self.adjacency.copy().astype(float)
        if kind == "binary":
            w.data[:] = 1.0
        elif kind == "inverse":
            w.data = 1.0 / np.maximum(w.data, 1e-12)
        elif kind == "gaussian":
            if bandwidth is None or bandwidth <= 0:
                raise ValueError("gaussian weights need a positive bandwidth")
            w.data = np.exp(-(w.data**2) / (2.0 * bandwidth**2))
        else:
            raise ValueError(f"unknown weight kind {kind!r}")
        return w


def _check_coords(coords: np.ndarray) -> np.ndarray:
    c = np.asarray(coords, dtype=float)
    if c.ndim != 2 or c.shape[1] != 2:
        raise ValueError(f"coords must have shape (n, 2), got {c.shape}")
    if not np.isfinite(c).all():
        raise ValueError("coords contain NaN or infinite values")
    if len(c) < 2:
        raise ValueError("need at least two cells to build a graph")
    return c


def _from_edges(c: np.ndarray, i: np.ndarray, j: np.ndarray, symmetric: bool) -> SpatialGraph:
    d = np.linalg.norm(c[i] - c[j], axis=1)
    d = np.maximum(d, 1e-12)  # keep coincident cells as explicit (non-zero) entries
    a = sp.coo_matrix((d, (i, j)), shape=(len(c), len(c))).tocsr()
    if symmetric:
        a = a.maximum(a.T).tocsr()
    a.sort_indices()
    return SpatialGraph(coords=c, adjacency=a)


def knn_graph(coords: np.ndarray, k: int = 8, symmetric: bool = True) -> SpatialGraph:
    """k-nearest-neighbour graph.

    Args:
        coords: (n, 2) cell centroids (same units as any distance thresholds used later).
        k: neighbours per cell, excluding the cell itself (capped at n - 1).
        symmetric: if True, keep edge i-j when either cell lists the other (union).
    """
    c = _check_coords(coords)
    if k < 1:
        raise ValueError("k must be >= 1")
    k_eff = min(k, len(c) - 1)
    _, idx = KDTree(c).query(c, k=k_eff + 1)
    i = np.repeat(np.arange(len(c)), k_eff)
    j = idx[:, 1:].ravel()
    return _from_edges(c, i, j, symmetric)


def delaunay_graph(coords: np.ndarray, max_edge_length: float | None = None) -> SpatialGraph:
    """Delaunay triangulation graph, optionally dropping edges longer than ``max_edge_length``."""
    c = _check_coords(coords)
    if len(c) < 3:
        raise ValueError("Delaunay triangulation needs at least three cells")
    tri = Delaunay(c)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    if max_edge_length is not None:
        keep = np.linalg.norm(c[e[:, 0]] - c[e[:, 1]], axis=1) <= max_edge_length
        e = e[keep]
    return _from_edges(c, e[:, 0], e[:, 1], symmetric=True)


def neighbourhood_composition(graph: SpatialGraph, labels: np.ndarray, n_types: int,
                              weights: sp.csr_matrix | None = None) -> np.ndarray:
    """Per-cell (weighted) fraction of neighbours of each type.

    Args:
        graph: spatial graph.
        labels: integer type per cell in ``[0, n_types)``; ``-1`` marks unassigned cells,
            which are ignored as neighbours.
        n_types: number of cell types.
        weights: optional edge weights (defaults to binary).

    Returns:
        (n, n_types) array whose rows sum to 1 (or 0 when no neighbour is assigned).
    """
    lab = np.asarray(labels)
    if lab.shape != (graph.n_cells,):
        raise ValueError("labels must have one entry per cell")
    if lab.max(initial=-1) >= n_types or lab.min(initial=0) < -1:
        raise ValueError("labels must lie in [-1, n_types)")
    w = graph.weights("binary") if weights is None else weights
    assigned = lab >= 0
    onehot = sp.csr_matrix(
        (np.ones(assigned.sum()), (np.flatnonzero(assigned), lab[assigned])), shape=(graph.n_cells, n_types)
    )
    counts = np.asarray((w @ onehot).todense())
    tot = counts.sum(axis=1, keepdims=True)
    return np.divide(counts, tot, out=np.zeros_like(counts), where=tot > 0)
