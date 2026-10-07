import numpy as np
import pytest

from spatialsurv.features import image_spatial_features, knn_neighbors, neighbor_enrichment, undirected_edges


def test_path_graph_exact_null():
    """Path A-A-B-B (edges 0-1, 1-2, 2-3). Worked by hand:
    observed A-B edges = 1. Over the 6 equally likely placements of the two A labels,
    A-B edge counts are {1,3,2,2,3,1}: mean 2, variance 4/6, so z = (1-2)/sqrt(2/3) = -1.2247.
    Homotypic A-A: observed 1, null values {1,0,0,1,0,1} -> mean 1/2."""
    edges = np.array([[0, 1], [1, 2], [2, 3]])
    labels = np.array([0, 0, 1, 1])
    obs, mu, sd = neighbor_enrichment(edges, labels, 2, n_perm=40000, rng=np.random.default_rng(0))
    assert obs[0, 1] == obs[1, 0] == 1
    assert obs[0, 0] == 1 and obs[1, 1] == 1
    assert mu[0, 1] == pytest.approx(2.0, abs=0.02)
    assert mu[0, 0] == pytest.approx(0.5, abs=0.02)
    z = (obs[0, 1] - mu[0, 1]) / sd[0, 1]
    assert z == pytest.approx(-1 / np.sqrt(2 / 3), abs=0.03)


def test_null_mean_matches_closed_form():
    """Under label permutation, E[#A-B edges] = E * 2 nA nB / (n (n-1)),
    E[#A-A edges] = E * nA (nA-1) / (n (n-1))."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 100, size=(300, 2))
    labels = rng.permutation(np.r_[np.zeros(100, int), np.ones(80, int), np.full(120, 2)])
    edges = undirected_edges(knn_neighbors(xy, 6))
    _, mu, _ = neighbor_enrichment(edges, labels, 3, n_perm=4000, rng=rng)
    E, n = len(edges), 300
    assert mu[0, 1] == pytest.approx(E * 2 * 100 * 80 / (n * (n - 1)), rel=0.02)
    assert mu[0, 0] == pytest.approx(E * 100 * 99 / (n * (n - 1)), rel=0.02)


def test_undirected_edges_dedup():
    nbrs = np.array([[1], [0], [1]])  # 0->1, 1->0 (same edge), 2->1
    assert undirected_edges(nbrs).tolist() == [[0, 1], [1, 2]]


def test_segregated_vs_mixed_image_features():
    """Tumor and T cells in separate halves -> strongly negative Tumor-T enrichment and
    almost no immune neighbours of tumor cells; interleaved grid -> positive."""
    rng = np.random.default_rng(0)
    g = np.array([(x, y) for x in range(20) for y in range(20)], float)
    seg = np.where(g[:, 0] < 10, "Tumor", "T")
    mixed = np.where((g[:, 0] + g[:, 1]) % 2 == 0, "Tumor", "T")  # checkerboard
    fs = image_spatial_features(g, seg, k=4, n_perm=200, rng=rng)
    fm = image_spatial_features(g, mixed, k=4, n_perm=200, rng=rng)
    assert fs["enrich_Tumor__T"] < -5 < 5 < fm["enrich_Tumor__T"]
    assert fs["tumor_immune_nbr_frac"] < 0.1 and fm["tumor_immune_nbr_frac"] > 0.9
    # checkerboard: interior cells have no same-type 4-neighbours; only border cells pick up a
    # same-colour diagonal, so mixing is close to (not exactly) 1
    assert fm["tumor_immune_mixing"] > 0.85 and fs["tumor_immune_mixing"] < 0.15
    assert np.isnan(fs["enrich_Tumor__Macrophage"])  # absent type -> NaN, not 0
