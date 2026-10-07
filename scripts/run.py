"""Reproduce every result and figure: `uv run python scripts/run.py`.

Expects data/raw populated by scripts/download_data.py.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
import sksurv

from spatialsurv.cv import cindex, run_models
from spatialsurv.data import CLINICAL_COLS, IMMUNE, load_annotations, load_cells, load_patients, survival_array
from spatialsurv.eval import bootstrap, seed_mean_cindex, summarize
from spatialsurv.features import (
    ALL_SPATIAL_COLS,
    SPATIAL_COLS,
    SPATIAL_COLS_LOGOE,
    aggregate_to_patient,
    composition_counts,
    spatial_features_by_image,
)
from spatialsurv.models import FeatureSet
from spatialsurv.plots import cindex_dotplot, example_cores_plot, headline_plot, km_plot, spatial_coef_plot
from spatialsurv.robustness import cv_repeat_deltas, permutation_null

ROOT = Path(__file__).resolve().parents[1]
RAW, PROC = ROOT / "data" / "raw", ROOT / "data" / "processed"
RES, FIG = ROOT / "results", ROOT / "figures"

CONFIG = dict(
    seeds=[0, 1, 2],
    outer_folds=5,
    inner_folds=5,
    knn_k=8,
    n_perm=200,
    spatial_seed=12345,
    clr_pseudocount=0.5,
    l1_ratio=0.5,
    n_alphas=20,
    alpha_min_ratio=0.01,
    coxnet_max_iter=100000,
    ridge_log10_min=-2.0,
    ridge_log10_max=3.0,
    n_boot=1000,
    boot_seed=2024,
    endpoints={"os": "event_os", "dss": "event_dss"},
    primary_endpoint="os",
    robustness_seeds=list(range(10)),  # CV repeats for the delta-stability check
    n_perm_null=20,  # spatial-block permutations (each run with the 3 main CV seeds)
)
M_CLIN, M_COMP, M_SPAT = "clinical", "clinical+composition", "clinical+composition+spatial"
M_SPAT_OE = "clinical+composition+spatial (log O/E, sensitivity)"


def slug(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", s).strip("_")


def build_features(cfg: dict, recompute_spatial: bool = False) -> tuple[pd.DataFrame, list[str], dict]:
    patients, cores = load_patients(RAW)
    cells = load_cells(RAW, cores["core"].tolist())
    ann = load_annotations(RAW)

    counts = composition_counts(cells, cores, ann["cluster"].tolist())
    counts.columns = [f"comp_{slug(l)}" for l in ann["label"]]

    PROC.mkdir(parents=True, exist_ok=True)
    cache = PROC / f"spatial_img_v2_k{cfg['knn_k']}_p{cfg['n_perm']}_s{cfg['spatial_seed']}.csv"
    if cache.exists() and not recompute_spatial:
        img = pd.read_csv(cache, index_col="core")
    else:
        img = spatial_features_by_image(cells, cfg["knn_k"], cfg["n_perm"], cfg["spatial_seed"])
        img.to_csv(cache)
    img.to_csv(RES / "spatial_features_by_image.csv")
    spatial = aggregate_to_patient(img, cores)

    X = patients.join(counts, how="inner").join(spatial, how="inner")
    assert len(X) == len(patients), "every patient must have cells"
    cohort = dict(
        cohort="Basel TMA (Jackson et al. 2020), tumor cores only",
        n_tumor_patients=patients.attrs["n_tumor_patients"],
        n_dropped_nonpositive_or_missing_time=patients.attrs["n_dropped_time"],
        n_patients=len(X),
        n_images=int(cores["core"].nunique()),
        patients_with_2_images=int((cores.groupby("PID").size() == 2).sum()),
        n_cells=int(len(cells)),
        events_os=int(X["event_os"].sum()),
        events_dss=int(X["event_dss"].sum()),
        median_followup_months_all=float(X["time"].median()),
        missing_clinical={c: int(X[c].isna().sum()) for c in CLINICAL_COLS if X[c].isna().any()},
        spatial_missing={c: int(X[c].isna().sum()) for c in ALL_SPATIAL_COLS if X[c].isna().any()},
        n_composition_features=counts.shape[1],
        n_spatial_features=len(SPATIAL_COLS),
        # sparse clinical categories behind the encoding decisions in README section 7
        grade1_patients=int((X["grade"] == 1).sum()),
        grade1_dss_events=int(X.loc[X["grade"] == 1, "event_dss"].sum()),
        pM1_patients=int((X["pM"] == 1).sum()),
    )
    return X, list(counts.columns), cohort


def analyse_endpoint(name: str, X: pd.DataFrame, comp_cols: list[str], cfg: dict) -> dict:
    y = survival_array(X["time"].to_numpy(), X[cfg["endpoints"][name]].to_numpy())
    fsets = [
        FeatureSet(M_CLIN, CLINICAL_COLS),
        FeatureSet(M_COMP, CLINICAL_COLS, comp_cols),
        FeatureSet(M_SPAT, CLINICAL_COLS, comp_cols, SPATIAL_COLS),
        FeatureSet(M_SPAT_OE, CLINICAL_COLS, comp_cols, SPATIAL_COLS_LOGOE),
    ]
    oof, records, coefs = run_models(fsets, X, y, cfg, cfg["seeds"])
    records.to_csv(RES / f"{name}_per_fold.csv", index=False)

    oof_df = pd.DataFrame({"PID": X.index, "time": y["time"], "event": y["event"]})
    for m, r in oof.items():
        for s, rs in zip(cfg["seeds"], r):
            oof_df[f"{m}|seed{s}"] = rs
    oof_df.to_csv(RES / f"{name}_oof_risk.csv", index=False)

    point = seed_mean_cindex(y, oof)
    per_seed = {m: [cindex(y, rs) for rs in r] for m, r in oof.items()}
    boot = bootstrap(y, oof, cfg["n_boot"], cfg["boot_seed"])
    boot.to_csv(RES / f"{name}_bootstrap.csv", index=False)
    comparisons = [(M_COMP, M_CLIN), (M_SPAT, M_COMP), (M_SPAT, M_CLIN), (M_SPAT_OE, M_COMP)]
    models, diffs = summarize(point, boot, comparisons)
    models["cindex_per_seed"] = [";".join(f"{v:.4f}" for v in per_seed[m]) for m in models["model"]]
    models.to_csv(RES / f"{name}_cindex.csv", index=False)
    diffs.to_csv(RES / f"{name}_paired_differences.csv", index=False)

    label = {"os": "Overall survival", "dss": "Breast-cancer-specific survival (sensitivity)"}[name]
    cindex_dotplot(models, FIG / f"{name}_cindex.png", f"{label}: N={len(X)}, events={int(y['event'].sum())}")

    out = dict(models=models, diffs=diffs, oof=oof, y=y, coefs=coefs, fsets=fsets)
    out["robustness"] = robustness(name, X, y, fsets, oof, models, cfg)
    return out


def robustness(name: str, X: pd.DataFrame, y: np.ndarray, fsets: list[FeatureSet], oof: dict,
               models: pd.DataFrame, cfg: dict) -> dict:
    by = {f.name: f for f in fsets}
    spatial_models = [by[M_SPAT], by[M_SPAT_OE]]
    rep = cv_repeat_deltas(by[M_COMP], spatial_models, X, y, cfg, cfg["robustness_seeds"])
    rep.to_csv(RES / f"{name}_cv_repeat_deltas.csv", index=False)
    point = models.set_index("model")["cindex"]
    out = {}
    for fs in spatial_models:
        d = rep.loc[rep["model"] == fs.name, "delta"]
        null_by_seed = permutation_null(oof[M_COMP], fs, X, y, cfg, cfg["seeds"], cfg["n_perm_null"])
        null = null_by_seed.mean(axis=1)
        obs = float(point[fs.name] - point[M_COMP])
        nd = pd.DataFrame(null_by_seed, columns=[f"null_delta_seed{s}" for s in cfg["seeds"]])
        nd.insert(0, "null_delta", null)
        nd.to_csv(RES / f"{name}_perm_null_{slug(fs.name)}.csv", index=False)
        out[fs.name] = dict(
            cv_repeats=len(d), delta_mean=float(d.mean()), delta_sd=float(d.std(ddof=1)),
            delta_min=float(d.min()), delta_max=float(d.max()), frac_repeats_gt0=float((d > 0).mean()),
            observed_delta_3seeds=obs, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
            null_max=float(null.max()),
            perm_p=float((1 + np.sum(null >= obs)) / (1 + len(null))),
        )
    pd.DataFrame(out).T.to_csv(RES / f"{name}_robustness.csv")
    return out


def spatial_coef_table(res: dict, comp_cols: list[str]) -> pd.DataFrame:
    fs = [f for f in res["fsets"] if f.name == M_SPAT][0]
    C = np.stack(res["coefs"][M_SPAT])  # (folds, features)
    cols = fs.columns
    rows = []
    for c in SPATIAL_COLS:
        v = C[:, cols.index(c)]
        rows.append(dict(feature=c, mean_coef=v.mean(), sd_coef=v.std(ddof=1),
                         selection_freq=float(np.mean(np.abs(v) > 1e-10))))
    return pd.DataFrame(rows).sort_values("mean_coef", key=np.abs, ascending=False).reset_index(drop=True)


def make_overview_figures(img: pd.DataFrame, cfg: dict) -> None:
    # Headline figure: discrimination per endpoint + spatial delta vs. permutation null.
    eps = {"OS": "os", "DSS": "dss"}
    spatial = {"z-score features": M_SPAT, "log O/E (sensitivity)": M_SPAT_OE}
    cidx = {k: pd.read_csv(RES / f"{v}_cindex.csv") for k, v in eps.items()}
    reps = {k: pd.read_csv(RES / f"{v}_cv_repeat_deltas.csv") for k, v in eps.items()}
    nulls = {
        k: {m: pd.read_csv(RES / f"{v}_perm_null_{slug(m)}.csv").filter(like="null_delta_seed").to_numpy().ravel()
            for m in spatial.values()}
        for k, v in eps.items()
    }
    headline_plot(cidx, reps, nulls, [M_CLIN, M_COMP, M_SPAT], spatial, FIG / "headline.png")

    # Example tissue. The mixing ratio also depends on how many immune cells there are
    # (sparse immune cells mostly touch tumor), so examples are drawn from images between the 40th
    # and 60th percentile of immune fraction, at the 10th/50th/90th percentile of mixing within that band.
    patients, cores = load_patients(RAW)
    cells = load_cells(RAW, cores["core"].tolist())
    n_cls = cells.groupby(["core", "coarse"]).size().unstack(fill_value=0)
    immune_frac = n_cls[list(IMMUNE)].sum(axis=1) / n_cls.sum(axis=1)
    feats = img.join(immune_frac.rename("immune_frac"))
    dep = feats[ALL_SPATIAL_COLS].corrwith(feats["immune_frac"], method="spearman").rename("spearman_vs_immune_frac")
    dep.to_csv(RES / "spatial_vs_immune_fraction.csv")
    ok = (n_cls["Tumor"] >= 200) & (n_cls[list(IMMUNE)].sum(axis=1) >= 50)
    cand = feats.loc[ok[ok].index].dropna(subset=["tumor_immune_mixing"])
    lo, hi = cand["immune_frac"].quantile([0.4, 0.6])
    cand = cand[cand["immune_frac"].between(lo, hi)]
    picks = []
    for q, lab in [(0.1, "low"), (0.5, "median"), (0.9, "high")]:
        target = cand["tumor_immune_mixing"].quantile(q)
        core = (cand["tumor_immune_mixing"] - target).abs().idxmin()
        r = cand.loc[core]
        picks.append((core, f"{lab} mixing: {r.tumor_immune_mixing:.2f}  (immune {r.immune_frac:.0%} of cells)\n"
                            f"Tumor-T enrichment z = {r.enrich_Tumor__T:.1f}; {int(r.n_cells)} cells"))
    example_cores_plot(cells, picks, FIG / "example_cores.png",
                       f"Example tumor cores at similar immune content ({lo:.0%}-{hi:.0%} immune cells), "
                       "low / median / high tumor-immune mixing")
    pd.DataFrame(picks, columns=["core", "caption"]).to_csv(RES / "example_cores.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recompute-spatial", action="store_true", help="ignore the cached image-level spatial features")
    args = ap.parse_args()
    t0 = time.time()
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    cfg = CONFIG
    X, comp_cols, cohort = build_features(cfg, args.recompute_spatial)
    X.to_csv(RES / "patient_features.csv")
    print(json.dumps(cohort, indent=2))

    summary = {"cohort": cohort}
    results = {}
    for ep in cfg["endpoints"]:
        t = time.time()
        results[ep] = analyse_endpoint(ep, X, comp_cols, cfg)
        print(f"\n[{ep}] ({time.time() - t:.0f}s)\n", results[ep]["models"].round(4).to_string(index=False))
        print(results[ep]["diffs"].round(4).to_string(index=False))
        summary[ep] = dict(
            cindex=results[ep]["models"].to_dict("records"),
            paired_differences=results[ep]["diffs"].to_dict("records"),
            robustness=results[ep]["robustness"],
        )
        print(pd.DataFrame(results[ep]["robustness"]).T.round(4).to_string())

    # Kaplan-Meier on the best model (primary endpoint, highest pooled OOF C-index).
    prim = results[cfg["primary_endpoint"]]
    mods = prim["models"][prim["models"]["model"].isin([M_CLIN, M_COMP, M_SPAT])]
    best = mods.sort_values("cindex", ascending=False)["model"].iloc[0]
    # Average of per-seed OOF risk ranks, so seeds are on a common scale.
    ranks = np.mean([pd.Series(r).rank().to_numpy() for r in prim["oof"][best]], axis=0)
    p = km_plot(prim["y"]["time"], prim["y"]["event"], ranks, FIG / "os_km_best_model.png",
                f"Best model ({best}): out-of-fold risk, median split")
    summary["km"] = dict(model=best, logrank_p=p, risk="mean OOF rank across 3 CV repeats")
    # Same plot for the clinical-only model, to show how much of the separation is clinical.
    ranks_c = np.mean([pd.Series(r).rank().to_numpy() for r in prim["oof"][M_CLIN]], axis=0)
    p_c = km_plot(prim["y"]["time"], prim["y"]["event"], ranks_c, FIG / "os_km_clinical_model.png",
                  "Clinical-only model: out-of-fold risk, median split")
    summary["km_clinical"] = dict(model=M_CLIN, logrank_p=p_c)
    hi_best, hi_clin = ranks > np.median(ranks), ranks_c > np.median(ranks_c)
    summary["km_group_agreement"] = float(np.mean(hi_best == hi_clin))

    for ep in cfg["endpoints"]:
        tab = spatial_coef_table(results[ep], comp_cols)
        tab.to_csv(RES / f"{ep}_spatial_coefficients.csv", index=False)
        summary[ep]["spatial_coefficients"] = tab.to_dict("records")
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        spatial_coef_plot(pd.read_csv(RES / f"{ep}_spatial_coefficients.csv"), FIG / f"{ep}_spatial_coefficients.png",
                          f"Spatial features in clinical+composition+spatial model ({lab}) - exploratory")

    # Diagnostic: permutation z-scores grow with the number of cells in an image, so check
    # how strongly each image-level spatial feature tracks image size.
    img = pd.read_csv(RES / "spatial_features_by_image.csv", index_col="core")
    diag = img[ALL_SPATIAL_COLS].corrwith(np.log(img["n_cells"]), method="spearman").rename("spearman_vs_log_ncells")
    diag.to_csv(RES / "spatial_vs_ncells.csv")
    summary["spatial_vs_ncells_spearman"] = diag.round(3).to_dict()

    make_overview_figures(img, cfg)

    summary["config"] = cfg
    summary["environment"] = dict(python=platform.python_version(), sklearn=sklearn.__version__,
                                  sksurv=sksurv.__version__, numpy=np.__version__, pandas=pd.__version__)
    summary["runtime_seconds"] = round(time.time() - t0, 1)
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"\nKM best={best} p={p:.3g}; total {summary['runtime_seconds']}s")


if __name__ == "__main__":
    main()
