"""Per-image spatial features and per-patient composition / spatial aggregation.

Everything here is computed within a single image (or a single patient) and uses no
survival information and no other patients, so it can be computed once before CV
without leakage. Transforms with parameters fitted across patients (imputation,
scaling) live in models.py and are fit inside CV folds.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .data import COARSE, IMMUNE

TYPES = list(COARSE)  # coarse classes, fixed order
# Interpretable neighbour-enrichment pairs (unordered).
PAIRS = [
    ("Tumor", "Tumor"),
    ("Tumor", "T"),
    ("Tumor", "Macrophage"),
    ("Tumor", "B"),
    ("Tumor", "Stroma"),
    ("Tumor", "Endothelial"),
    ("T", "Macrophage"),
    ("T", "B"),
    ("Stroma", "T"),
]
MIN_CELLS_PER_TYPE = 5  # enrichment for a pair is NaN if either type has fewer cells in the image


def pair_name(a: str, b: str) -> str:
    return f"enrich_{a}__{b}"


SPATIAL_COLS = [pair_name(a, b) for a, b in PAIRS] + ["tumor_immune_nbr_frac", "tumor_immune_mixing"]


# --------------------------------------------------------------------------- graphs


def knn_neighbors(xy: np.ndarray, k: int) -> np.ndarray:
    """Directed kNN: (n, k') neighbour indices excluding self, k' = min(k, n - 1)."""
    k_eff = min(k, len(xy) - 1)
    _, idx = cKDTree(xy).query(xy, k=k_eff + 1)
    return idx[:, 1:]


def undirected_edges(nbrs: np.ndarray) -> np.ndarray:
    """Symmetrised kNN graph as an (m, 2) array of unique undirected edges with i < j."""
    n, k = nbrs.shape
    i = np.repeat(np.arange(n), k)
    j = nbrs.ravel()
    e = np.stack([np.minimum(i, j), np.maximum(i, j)], axis=1)
    return np.unique(e, axis=0)


# --------------------------------------------------------------------------- enrichment


def _pair_counts(edges: np.ndarray, labels: np.ndarray, n_types: int) -> np.ndarray:
    """Symmetric (T, T) matrix of undirected edge counts between types.

    Diagonal entries count homotypic edges; off-diagonals count each heterotypic edge once.
    labels may be 2-D (n_perm, n) to count many label permutations at once.
    """
    la = labels[..., edges[:, 0]]
    lb = labels[..., edges[:, 1]]
    lo, hi = np.minimum(la, lb), np.maximum(la, lb)
    flat = lo * n_types + hi
    if flat.ndim == 1:
        c = np.bincount(flat, minlength=n_types * n_types).reshape(n_types, n_types)
    else:
        offs = (np.arange(flat.shape[0]) * n_types * n_types)[:, None]
        c = np.bincount((flat + offs).ravel(), minlength=flat.shape[0] * n_types * n_types)
        c = c.reshape(flat.shape[0], n_types, n_types)
    return c + np.swapaxes(c, -1, -2) - c * np.eye(n_types, dtype=c.dtype)


def neighbor_enrichment(
    edges: np.ndarray,
    labels: np.ndarray,
    n_types: int,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Edge-count enrichment of each type pair vs. a within-image label-permutation null.

    Returns (observed, null_mean, null_sd), each (T, T). The graph is held fixed and cell
    labels are shuffled, which preserves composition and tissue architecture.
    z = (observed - null_mean) / null_sd.
    """
    obs = _pair_counts(edges, labels, n_types).astype(float)
    perms = np.stack([rng.permutation(labels) for _ in range(n_perm)])
    null = _pair_counts(edges, perms, n_types).astype(float)
    return obs, null.mean(axis=0), null.std(axis=0, ddof=1)


def image_spatial_features(
    xy: np.ndarray, coarse: np.ndarray, k: int, n_perm: int, rng: np.random.Generator
) -> dict[str, float]:
    """Spatial features for one image. coarse: array of coarse class names per cell."""
    type_idx = {t: i for i, t in enumerate(TYPES)}
    labels = np.array([type_idx[c] for c in coarse])
    counts = np.bincount(labels, minlength=len(TYPES))
    nbrs = knn_neighbors(xy, k)
    edges = undirected_edges(nbrs)
    obs, mu, sd = neighbor_enrichment(edges, labels, len(TYPES), n_perm, rng)

    feats: dict[str, float] = {}
    for a, b in PAIRS:
        ia, ib = type_idx[a], type_idx[b]
        ok = counts[ia] >= MIN_CELLS_PER_TYPE and counts[ib] >= MIN_CELLS_PER_TYPE and sd[ia, ib] > 0
        feats[pair_name(a, b)] = float((obs[ia, ib] - mu[ia, ib]) / sd[ia, ib]) if ok else np.nan

    is_immune = np.isin(labels, [type_idx[t] for t in IMMUNE])
    is_tumor = labels == type_idx["Tumor"]
    # Mean fraction of immune cells among each tumor cell's k nearest neighbours.
    feats["tumor_immune_nbr_frac"] = float(is_immune[nbrs[is_tumor]].mean()) if is_tumor.any() else np.nan
    # Tumor-immune mixing (after Keren et al. 2018), bounded variant:
    # tumor-immune edges / (tumor-immune + immune-immune edges). High = immune cells mixed
    # into tumor; low = immune cells segregated into their own compartments.
    ei, ej = edges[:, 0], edges[:, 1]
    n_ti = np.sum((is_tumor[ei] & is_immune[ej]) | (is_immune[ei] & is_tumor[ej]))
    n_ii = np.sum(is_immune[ei] & is_immune[ej])
    feats["tumor_immune_mixing"] = float(n_ti / (n_ti + n_ii)) if (n_ti + n_ii) > 0 else np.nan
    return feats


def spatial_features_by_image(cells: pd.DataFrame, k: int, n_perm: int, seed: int) -> pd.DataFrame:
    """One row per image (core) with SPATIAL_COLS and n_cells."""
    rng = np.random.default_rng(seed)
    rows = []
    for core, df in cells.groupby("core", sort=True):
        f = image_spatial_features(df[["x", "y"]].to_numpy(), df["coarse"].to_numpy(), k, n_perm, rng)
        f["core"], f["n_cells"] = core, len(df)
        rows.append(f)
    return pd.DataFrame(rows).set_index("core")


def aggregate_to_patient(img: pd.DataFrame, cores: pd.DataFrame) -> pd.DataFrame:
    """Cell-count-weighted mean of image features per patient (NaN-aware)."""
    df = img.join(cores.set_index("core")["PID"], how="inner")
    out = {}
    for pid, g in df.groupby("PID"):
        w = g["n_cells"].to_numpy(float)
        row = {}
        for c in SPATIAL_COLS:
            v = g[c].to_numpy(float)
            m = ~np.isnan(v)
            row[c] = float(np.average(v[m], weights=w[m])) if m.any() else np.nan
        out[pid] = row
    return pd.DataFrame.from_dict(out, orient="index")[SPATIAL_COLS]


# --------------------------------------------------------------------------- composition


def composition_counts(cells: pd.DataFrame, cores: pd.DataFrame, labels: list[int]) -> pd.DataFrame:
    """Per-patient cell counts per metacluster (cells pooled across a patient's cores).

    Raw counts are returned; the CLR transform (row-wise, with pseudocount) is applied in
    the model pipeline.
    """
    df = cells.merge(cores, on="core", how="inner")
    tab = pd.crosstab(df["PID"], df["cluster"]).reindex(columns=labels, fill_value=0)
    return tab


def clr(counts: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Centred log-ratio of closed compositions. Row-wise: each row uses only itself."""
    x = np.asarray(counts, dtype=float) + pseudocount
    x = x / x.sum(axis=1, keepdims=True)
    lx = np.log(x)
    return lx - lx.mean(axis=1, keepdims=True)
