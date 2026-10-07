"""End-to-end toy: a known spatial pattern (T-cell infiltration into the tumor region)
drives survival while composition is held fixed and clinical variables are noise.
The spatial model must beat the clinical and composition models."""

import numpy as np
import pandas as pd

from spatialsurv.cv import run_models
from spatialsurv.data import survival_array
from spatialsurv.eval import seed_mean_cindex
from spatialsurv.features import SPATIAL_COLS, aggregate_to_patient, composition_counts, spatial_features_by_image
from spatialsurv.models import FeatureSet


def simulate(n_pat=150, n_cells=300, seed=0):
    rng = np.random.default_rng(seed)
    rows, cores = [], []
    infil = rng.uniform(0, 1, n_pat)  # fraction of T cells inside the tumor half
    for p in range(n_pat):
        n_t, n_tc, n_st = int(0.6 * n_cells), int(0.15 * n_cells), int(0.25 * n_cells)
        tumor = np.c_[rng.uniform(0, 100, n_t), rng.uniform(0, 200, n_t)]
        stroma = np.c_[rng.uniform(100, 200, n_st), rng.uniform(0, 200, n_st)]
        inside = rng.random(n_tc) < infil[p]
        tcell = np.c_[np.where(inside, rng.uniform(0, 100, n_tc), rng.uniform(100, 200, n_tc)),
                      rng.uniform(0, 200, n_tc)]
        xy = np.r_[tumor, tcell, stroma]
        lab = ["Tumor"] * n_t + ["T"] * n_tc + ["Stroma"] * n_st
        clus = [20] * n_t + [3] * n_tc + [10] * n_st
        core = f"img{p}"
        rows.append(pd.DataFrame({"core": core, "x": xy[:, 0], "y": xy[:, 1], "coarse": lab, "cluster": clus}))
        cores.append((core, f"P{p}"))
    cells = pd.concat(rows, ignore_index=True)
    cores = pd.DataFrame(cores, columns=["core", "PID"])
    # Higher infiltration -> lower hazard (strong effect), ~30% censoring.
    t_event = rng.exponential(1 / np.exp(-3.0 * (infil - 0.5)))
    t_cens = rng.exponential(2.5, n_pat)
    time, event = np.minimum(t_event, t_cens), t_event <= t_cens
    clin = pd.DataFrame({"c_age": rng.normal(size=n_pat), "c_grade": rng.integers(1, 4, n_pat).astype(float)},
                        index=cores["PID"])
    return cells, cores, clin, survival_array(time, event)


def test_spatial_model_wins_when_spatial_pattern_drives_survival():
    cells, cores, clin, y = simulate()
    img = spatial_features_by_image(cells, k=8, n_perm=50, seed=0)
    spat = aggregate_to_patient(img, cores)
    comp = composition_counts(cells, cores, [3, 10, 20])
    comp.columns = ["comp_T", "comp_Stroma", "comp_Tumor"]
    X = clin.join(comp).join(spat).loc[cores["PID"]]
    cfg = dict(outer_folds=5, inner_folds=3, clr_pseudocount=0.5, l1_ratio=0.5, n_alphas=10,
               alpha_min_ratio=0.01, coxnet_max_iter=100000, ridge_log10_min=-1, ridge_log10_max=2)
    clin_cols, comp_cols = list(clin.columns), list(comp.columns)
    fsets = [
        FeatureSet("clinical", clin_cols),
        FeatureSet("composition", clin_cols, comp_cols),
        FeatureSet("spatial", clin_cols, comp_cols, SPATIAL_COLS),
    ]
    oof, records, _ = run_models(fsets, X, y, cfg, seeds=[0])
    c = seed_mean_cindex(y, oof)
    assert c["spatial"] > 0.65
    assert c["spatial"] > c["composition"] + 0.1
    assert c["spatial"] > c["clinical"] + 0.1
    assert abs(c["random baseline"] - 0.5) < 0.08
