import numpy as np
import pandas as pd
import pytest

from celesta_lite import (
    UNKNOWN,
    CelestaConfig,
    Signature,
    cell_type_scores,
    expression_probabilities,
    fit_marker_models,
    knn_graph,
    neighbourhood_composition,
    run_celesta,
)
from celesta_lite.preprocessing import _crossing_point, _equal_weight_gmm


def test_equal_weight_gmm_symmetric_crossing():
    rng = np.random.default_rng(0)
    x = np.r_[rng.normal(0, 0.5, 3000), rng.normal(4, 0.5, 3000)]
    mu, var = _equal_weight_gmm(x)
    assert mu == pytest.approx([0, 4], abs=0.05)
    assert _crossing_point(mu, var) == pytest.approx(2.0, abs=0.05)


def test_crossing_point_closed_form_unequal_variances():
    # equal weights: N(x;0,1) = N(x;3,4) solved by hand -> root between the means
    mu, var = np.array([0.0, 3.0]), np.array([1.0, 4.0])
    xc = _crossing_point(mu, var)
    pdf = lambda x, m, v: np.exp(-((x - m) ** 2) / (2 * v)) / np.sqrt(2 * np.pi * v)  # noqa: E731
    assert 0 < xc < 3 and pdf(xc, 0, 1) == pytest.approx(pdf(xc, 3, 4), rel=1e-9)


def test_expression_probabilities_rescaled_to_unit_interval():
    rng = np.random.default_rng(1)
    x = pd.DataFrame({"m": np.r_[rng.normal(0, 0.3, 500), rng.normal(3, 0.3, 500)]})
    ep = expression_probabilities(x, fit_marker_models(x))["m"]
    assert ep.min() == 0.0 and ep.max() == 1.0


def test_scores_formula():
    Z, F = cell_type_scores(np.array([[0.9, 0.2, 0.8]]), np.array([[1, 0, 1], [0, 1, np.nan]], float))
    assert Z[0, 0] == pytest.approx(1 - (0.01 + 0.04 + 0.04) / 3)
    assert Z[0, 1] == pytest.approx(1 - (0.81 + 0.64) / 2)
    assert F.sum() == pytest.approx(1.0)


def _sig(rows):
    return Signature(pd.DataFrame(rows).T)


def test_signature_round_validation():
    with pytest.raises(ValueError):  # round-2 parent must exist in an earlier round
        _sig({"A": {"round": 1, "parent": "", "m": 1}, "B": {"round": 2, "parent": "Z", "m": 0}})


def _toy_tissue(seed=0):
    """Tumour disc (PanCK+) inside stroma (Vim+). 25% of tumour cells are ambiguous
    (PanCK EP 0.52, Vim EP 0.55): not anchors, and their markers lean slightly towards
    Stroma, so marker-only scoring mislabels them."""
    rng = np.random.default_rng(seed)
    g = np.array([(x, y) for x in range(40) for y in range(40)], float) * 10.0
    tumour = np.linalg.norm(g - 200, axis=1) < 120
    ep = np.where(tumour[:, None], [0.95, 0.05], [0.05, 0.95]) + rng.normal(0, 0.02, (len(g), 2))
    amb = tumour & (rng.random(len(g)) < 0.25)
    ep[amb] = [0.52, 0.55]
    ep = np.clip(ep, 0, 1)
    sig = _sig({"Tumor": {"round": 1, "parent": "", "PanCK": 1, "Vim": 0},
                "Stroma": {"round": 1, "parent": "", "PanCK": 0, "Vim": 1}})
    return pd.DataFrame(ep, columns=["PanCK", "Vim"]), g, tumour, amb, sig


def test_spatial_refinement_rescues_ambiguous_cells():
    ep, xy, tumour, amb, sig = _toy_tissue()
    res = run_celesta(ep, xy, sig, CelestaConfig(bandwidth=50.0, filter_artifacts=False))
    lab = res.labels.to_numpy()
    assert not res.is_anchor[amb].any()  # ambiguous cells are not anchors
    assert (lab[amb] == "Tumor").mean() > 0.95  # but the MRF assigns them via their neighbours
    assert (lab[~tumour] == "Stroma").mean() > 0.99
    assert res.rounds[0]["anchors"] == int((~amb).sum())


def test_without_spatial_term_ambiguous_cells_are_mislabelled():
    ep, xy, tumour, amb, sig = _toy_tissue()
    res = run_celesta(ep, xy, sig, CelestaConfig(bandwidth=50.0, scale_factor=0.0, filter_artifacts=False))
    # gamma = 0 removes the spatial term: the ambiguous cells follow their markers -> Stroma
    assert (res.labels[amb] == "Stroma").mean() > 0.95
    assert (res.marker_only[amb] == "Stroma").all()


def test_lineage_rounds_refine_only_parent_cells():
    rng = np.random.default_rng(0)
    n = 400
    xy = rng.uniform(0, 1000, (n, 2))
    is_imm = np.arange(n) < 200
    is_t = is_imm & (np.arange(n) < 100)
    ep = pd.DataFrame({"CD45": np.where(is_imm, 0.95, 0.05), "CD3": np.where(is_t, 0.95, 0.05),
                       "PanCK": np.where(is_imm, 0.05, 0.95)})
    sig = _sig({"Immune": {"round": 1, "parent": "", "CD45": 1, "PanCK": 0, "CD3": np.nan},
                "Tumor": {"round": 1, "parent": "", "CD45": 0, "PanCK": 1, "CD3": np.nan},
                "T cell": {"round": 2, "parent": "Immune", "CD45": np.nan, "PanCK": np.nan, "CD3": 1},
                "Other immune": {"round": 2, "parent": "Immune", "CD45": np.nan, "PanCK": np.nan, "CD3": 0}})
    res = run_celesta(ep, xy, sig, CelestaConfig(filter_artifacts=False))
    lab = res.labels.to_numpy()
    assert (lab[is_t] == "T cell").all() and (lab[is_imm & ~is_t] == "Other immune").all()
    assert (lab[~is_imm] == "Tumor").all()
    assert res.rounds[1]["candidates"] == 200


def test_artifact_cells_are_not_assigned():
    ep = pd.DataFrame({"a": [0.1, 0.95, 0.2], "b": [0.1, 0.95, 0.9]})
    sig = _sig({"A": {"round": 1, "parent": "", "a": 1, "b": 0}, "B": {"round": 1, "parent": "", "a": 0, "b": 1}})
    res = run_celesta(ep, np.array([[0, 0], [1, 0], [2, 0]], float), sig)
    assert res.is_artifact.tolist() == [True, True, False]
    assert res.labels.iloc[0] == UNKNOWN and res.labels.iloc[2] == "B"


def test_graph_helpers_still_work():
    xy = np.c_[np.arange(4.0), np.zeros(4)]
    comp = neighbourhood_composition(knn_graph(xy, k=1), np.array([0, 1, -1, 1]), n_types=2)
    np.testing.assert_allclose(comp.sum(axis=1), [1, 1, 1, 0])


def test_posterior_ep_is_scale_free_and_monotone():
    rng = np.random.default_rng(2)
    dim = pd.DataFrame({"m": np.r_[rng.normal(0.10, 0.02, 3000), rng.normal(0.30, 0.08, 1000)]})
    bright = dim * 30.0  # same distribution on a 30x larger scale
    ep_dim = expression_probabilities(dim, fit_marker_models(dim), method="posterior")["m"]
    ep_bright = expression_probabilities(bright, fit_marker_models(bright), method="posterior")["m"]
    np.testing.assert_allclose(ep_dim, ep_bright, atol=1e-3)  # invariant to intensity scale
    s = ep_dim.to_numpy()[np.argsort(dim["m"].to_numpy())]
    assert np.all(np.diff(s) >= 0)  # monotone in intensity
    assert ep_dim[3000:].mean() > 0.8 and ep_dim[:3000].mean() < 0.2
    # the reference sigmoid (slope 1) is nearly flat on the dim scale before min-max
    ref = expression_probabilities(dim, fit_marker_models(dim), method="reference")["m"]
    assert ref[3000:].mean() - ref[:3000].mean() < ep_dim[3000:].mean() - ep_dim[:3000].mean()
