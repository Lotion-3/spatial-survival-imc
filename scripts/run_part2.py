"""Part II (pre-registered in ANALYSIS_PLAN.md): METABRIC external replication (A1),
Basel -> METABRIC transfer (A2), composition-conditioned spatial features (B),
calibration / time-dependent AUC (C), ER-negative subgroup (D).

Run after scripts/run.py (it reuses Basel results/patient_features.csv and OOF risks):
    uv run python scripts/run_part2.py
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from spatialsurv.cv import cindex, run_models, select_alpha
from spatialsurv.data import CLINICAL_COLS, COARSE, survival_array
from spatialsurv.eval import bootstrap, extra_metrics, seed_mean_cindex, summarize
from spatialsurv.features import (
    ALL_SPATIAL_COLS,
    SPATIAL_COLS,
    SPATIAL_COLS_LOGOE,
    aggregate_to_patient,
    spatial_features_by_image,
)
from spatialsurv.metabric import MIN_CELLS, load_metabric_cells, load_metabric_clinical
from spatialsurv.models import FeatureSet, fit_final, predict_risk
from spatialsurv.plots import cindex_dotplot, headline_plot
from spatialsurv.robustness import cv_repeat_deltas, permutation_null

ROOT = Path(__file__).resolve().parents[1]
RAW_MB, PROC = ROOT / "data" / "raw_metabric", ROOT / "data" / "processed"
RES1, OUT, FIG = ROOT / "results", ROOT / "results" / "part2", ROOT / "figures"

# Identical modelling settings to part I (scripts/run.py CONFIG).
CONFIG = dict(
    seeds=[0, 1, 2], outer_folds=5, inner_folds=5, knn_k=8, n_perm=200, spatial_seed=12345,
    clr_pseudocount=0.5, l1_ratio=0.5, n_alphas=20, alpha_min_ratio=0.01, coxnet_max_iter=100000,
    ridge_log10_min=-2.0, ridge_log10_max=3.0, n_boot=1000, boot_seed=2024,
    robustness_seeds=list(range(10)), n_perm_null=20, auc_times=[60.0, 120.0],
    er_neg_min_dss_events=40,
)
ENDPOINTS = {"os": "event_os", "dss": "event_dss"}
M_CLIN, M_COMP = "clinical", "clinical+composition"
M_SPAT, M_OE, M_RES = "+spatial (z)", "+spatial (log O/E)", "+spatial (residualised)"
M_CEL = "+spatial (log O/E, CELESTA labels)"  # Amendment 3, section E
CELESTA_LABELS = ROOT / "results" / "celesta" / "v1" / "cell_labels.csv.gz"
CEL_COLS = [f"cel_{c}" for c in SPATIAL_COLS_LOGOE]


def slug(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", s).strip("_")


def coarse_counts(counts: pd.DataFrame, cluster_to_class) -> pd.DataFrame:
    out = pd.DataFrame(index=counts.index)
    for cls in COARSE:
        cols = [c for c in counts.columns if cluster_to_class(c) == cls]
        out[f"coarse_{cls}"] = counts[cols].sum(axis=1)
    return out


# --------------------------------------------------------------------------- data


def basel_table() -> tuple[pd.DataFrame, list[str], list[str]]:
    X = pd.read_csv(RES1 / "patient_features.csv", index_col=0)
    comp = [c for c in X.columns if c.startswith("comp_")]
    cls_of = {cid: name for name, ids in COARSE.items() for cid in ids}
    cc = coarse_counts(X[comp], lambda c: cls_of[int(c.split("_")[1])])
    return X.join(cc), comp, list(cc.columns)


def metabric_table(cfg: dict) -> tuple[pd.DataFrame, list[str], list[str], dict]:
    clin = load_metabric_clinical(RAW_MB)
    cells = load_metabric_cells(RAW_MB)
    if CELESTA_LABELS.exists():
        lab = pd.read_csv(CELESTA_LABELS, usecols=["core", "metabric_id", "celesta_coarse"])
        if len(lab) != len(cells) or not (lab["core"].to_numpy() == cells["core"].to_numpy()).all()                 or not (lab["metabric_id"].to_numpy() == cells["PID"].to_numpy()).all():
            raise ValueError("CELESTA labels are not aligned with the METABRIC cells")
        cells["celesta_coarse"] = lab["celesta_coarse"].to_numpy()
    n_cells = cells.groupby("PID").size()
    has_epi = cells.groupby("PID")["is_epithelial"].max() == 1
    keep_cells = n_cells.index[(n_cells >= MIN_CELLS) & has_epi]
    cells = cells[cells["PID"].isin(keep_cells)]

    cache = PROC / f"metabric_spatial_img_k{cfg['knn_k']}_p{cfg['n_perm']}_s{cfg['spatial_seed']}.csv"
    if cache.exists():
        img = pd.read_csv(cache, index_col="core")
    else:
        PROC.mkdir(parents=True, exist_ok=True)
        img = spatial_features_by_image(cells, cfg["knn_k"], cfg["n_perm"], cfg["spatial_seed"])
        img.to_csv(cache)
    cores = cells[["core", "PID"]].drop_duplicates()
    spatial = aggregate_to_patient(img, cores)
    if "celesta_coarse" in cells:  # section E: same features from CELESTA labels, Unknown cells dropped
        cache_e = PROC / f"metabric_celesta_spatial_img_k{cfg['knn_k']}_p{cfg['n_perm']}_s{cfg['spatial_seed']}.csv"
        if cache_e.exists():
            img_e = pd.read_csv(cache_e, index_col="core")
        else:
            known = cells[cells["celesta_coarse"] != "Unknown"]
            img_e = spatial_features_by_image(known.assign(coarse=known["celesta_coarse"]),
                                              cfg["knn_k"], cfg["n_perm"], cfg["spatial_seed"])
            img_e.to_csv(cache_e)
        sp_e = aggregate_to_patient(img_e, cores)[SPATIAL_COLS_LOGOE]
        sp_e.columns = CEL_COLS
        spatial = spatial.join(sp_e, how="left")

    native = pd.crosstab(cells["PID"], cells["cellPhenotype"])
    native.columns = [f"pheno_{slug(c)}" for c in native.columns]
    cc = pd.crosstab(cells["PID"], cells["coarse"]).reindex(columns=list(COARSE), fill_value=0)
    cc.columns = [f"coarse_{c}" for c in cc.columns]

    n_img_pat = cores["PID"].nunique()
    X = clin.join(native, how="inner").join(cc, how="inner").join(spatial, how="inner")
    n_joined = len(X)
    X = X[X["vital_known"] & (X["time"] > 0)]
    cohort = dict(
        cohort="METABRIC IMC (Danenberg et al. 2022), invasive tumour images, artefact cells removed",
        patients_with_images_after_cell_filters=int(n_img_pat),
        patients_dropped_cells_lt_min_or_no_epithelium=int(n_cells.size - len(keep_cells)),
        patients_with_clinical=int(n_joined),
        n_patients=int(len(X)),
        n_images=int(cores[cores["PID"].isin(X.index)]["core"].nunique()),
        n_cells=int(cells["PID"].isin(X.index).sum()),
        events_os_180m=int(X["event_os"].sum()),
        events_dss_180m=int(X["event_dss"].sum()),
        median_time_months=float(X["time"].median()),
        er_negative_patients=int((X["ER"] == 0).sum()),
        er_negative_dss_events=int(X.loc[X["ER"] == 0, "event_dss"].sum()),
        missing_clinical={c: int(X[c].isna().sum()) for c in CLINICAL_COLS if X[c].isna().any()},
        n_native_phenotypes=int(native.shape[1]),
        median_cells_per_patient=float(n_cells[X.index].median()),
    )
    img.loc[img.index.isin(cores.loc[cores["PID"].isin(X.index), "core"])].to_csv(OUT / "metabric_spatial_by_image.csv")
    return X, list(native.columns), list(cc.columns), cohort


# --------------------------------------------------------------------------- analyses


def compare(name: str, X: pd.DataFrame, y: np.ndarray, fsets: list[FeatureSet], base: str,
            spatial_models: list[str], cfg: dict, robust: bool = True) -> dict:
    """Part I protocol: nested CV x3 seeds, shared bootstrap, paired diffs vs `base`,
    CV-repeat deltas and permutation null for each spatial model."""
    oof, records, _ = run_models(fsets, X, y, cfg, cfg["seeds"])
    records.to_csv(OUT / f"{name}_per_fold.csv", index=False)
    point = seed_mean_cindex(y, oof)
    boot = bootstrap(y, oof, cfg["n_boot"], cfg["boot_seed"])
    comps = [(M_COMP, M_CLIN)] + [(m, base) for m in spatial_models]
    models, diffs = summarize(point, boot, comps)
    models["cindex_per_seed"] = [";".join(f"{cindex(y, r):.4f}" for r in oof[m]) for m in models["model"]]
    models.to_csv(OUT / f"{name}_cindex.csv", index=False)
    diffs.to_csv(OUT / f"{name}_paired_differences.csv", index=False)
    oof_df = pd.DataFrame({"PID": X.index, "time": y["time"], "event": y["event"]})
    for m, r in oof.items():
        for s, rs in zip(cfg["seeds"], r):
            oof_df[f"{m}|seed{s}"] = rs
    oof_df.to_csv(OUT / f"{name}_oof_risk.csv", index=False)
    rob = {}
    if robust:
        by = {f.name: f for f in fsets}
        rep = cv_repeat_deltas(by[base], [by[m] for m in spatial_models], X, y, cfg, cfg["robustness_seeds"])
        rep.to_csv(OUT / f"{name}_cv_repeat_deltas.csv", index=False)
        for m in spatial_models:
            nb = permutation_null(oof[base], by[m], X, y, cfg, cfg["seeds"], cfg["n_perm_null"])
            null = nb.mean(axis=1)
            nd = pd.DataFrame(nb, columns=[f"null_delta_seed{s}" for s in cfg["seeds"]])
            nd.insert(0, "null_delta", null)
            nd.to_csv(OUT / f"{name}_perm_null_{slug(m)}.csv", index=False)
            d = rep.loc[rep["model"] == m, "delta"]
            obs = point[m] - point[base]
            rob[m] = dict(cv_repeats=len(d), delta_mean=float(d.mean()), delta_sd=float(d.std(ddof=1)),
                          delta_min=float(d.min()), delta_max=float(d.max()),
                          frac_repeats_gt0=float((d > 0).mean()), observed_delta_3seeds=float(obs),
                          null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)), null_max=float(null.max()),
                          perm_p=float((1 + np.sum(null >= obs)) / (1 + len(null))))
        pd.DataFrame(rob).T.to_csv(OUT / f"{name}_robustness.csv")
    return dict(oof=oof, models=models, diffs=diffs, robustness=rob)


def h1_verdict(res: dict) -> dict:
    """Pre-registered decision rule for H1 (METABRIC, DSS, log O/E vs composition)."""
    d = res["diffs"].set_index("comparison").loc[f"{M_OE} minus {M_COMP}"]
    r = res["robustness"][M_OE]
    a, b, c = bool(d.ci_low > 0), bool(r["perm_p"] <= 0.05), bool(r["delta_mean"] > 0)
    verdict = "supported" if (a and b and c) else ("refuted" if d.delta <= 0 else "inconclusive")
    return dict(delta=float(d.delta), ci_low=float(d.ci_low), ci_high=float(d.ci_high), perm_p=r["perm_p"],
                cv_repeat_mean=r["delta_mean"], criterion_a_ci_above_0=a, criterion_b_perm_p_le_005=b,
                criterion_c_repeat_mean_gt0=c, verdict=verdict)


def transfer(Xb: pd.DataFrame, Xm: pd.DataFrame, coarse_cols: list[str], cfg: dict) -> dict:
    """A2: fit on all Basel patients (alpha by inner CV), evaluate once on METABRIC."""
    fsets = [
        FeatureSet(M_CLIN, CLINICAL_COLS),
        FeatureSet(M_COMP, CLINICAL_COLS, coarse_cols),
        FeatureSet(M_SPAT, CLINICAL_COLS, coarse_cols, SPATIAL_COLS),
        FeatureSet(M_OE, CLINICAL_COLS, coarse_cols, SPATIAL_COLS_LOGOE),
    ]
    out = {}
    for ep, col in ENDPOINTS.items():
        yb = survival_array(Xb["time"].to_numpy(), Xb[col].to_numpy())
        ym = survival_array(Xm["time"].to_numpy(), Xm[col].to_numpy())
        risk, alphas = {}, {}
        for fs in fsets:
            a, _, _ = select_alpha(fs, Xb[fs.columns], yb, cfg, seed=0)
            m = fit_final(fs, Xb[fs.columns], yb, a, cfg)
            risk[fs.name] = predict_risk(m, Xm[fs.columns])[None, :]
            alphas[fs.name] = m.chosen_alpha_
        point = seed_mean_cindex(ym, risk)
        boot = bootstrap(ym, risk, cfg["n_boot"], cfg["boot_seed"])
        models, diffs = summarize(point, boot, [(M_COMP, M_CLIN), (M_SPAT, M_COMP), (M_OE, M_COMP)])
        models["alpha"] = models["model"].map(alphas)
        models.to_csv(OUT / f"transfer_{ep}_cindex.csv", index=False)
        diffs.to_csv(OUT / f"transfer_{ep}_paired_differences.csv", index=False)
        out[ep] = dict(models=models.to_dict("records"), diffs=diffs.to_dict("records"))
    return out


def calibration(name: str, y: np.ndarray, oof: dict, cfg: dict) -> pd.DataFrame:
    t = extra_metrics(y, oof, cfg["auc_times"], cfg["n_boot"], cfg["boot_seed"])
    t.to_csv(OUT / f"{name}_calibration_auc.csv", index=False)
    return t


# --------------------------------------------------------------------------- main


def main() -> None:
    t0 = time.time()
    cfg = CONFIG
    OUT.mkdir(parents=True, exist_ok=True)
    summary: dict = {"config": cfg}

    # ---- METABRIC data
    Xm, native_cols, coarse_cols, cohort = metabric_table(cfg)
    Xm.to_csv(OUT / "metabric_patient_features.csv")
    summary["metabric_cohort"] = cohort
    print(json.dumps(cohort, indent=2))

    # ---- A1 (+ B in METABRIC, pre-registered secondary)
    fsets_m = [
        FeatureSet(M_CLIN, CLINICAL_COLS),
        FeatureSet(M_COMP, CLINICAL_COLS, native_cols),
        FeatureSet(M_SPAT, CLINICAL_COLS, native_cols, SPATIAL_COLS),
        FeatureSet(M_OE, CLINICAL_COLS, native_cols, SPATIAL_COLS_LOGOE),
        FeatureSet(M_RES, CLINICAL_COLS, native_cols, SPATIAL_COLS_LOGOE,
                   residualize_spatial=True),
        FeatureSet(M_CEL, CLINICAL_COLS, native_cols, CEL_COLS),
    ]
    a1 = {}
    for ep, col in ENDPOINTS.items():
        y = survival_array(Xm["time"].to_numpy(), Xm[col].to_numpy())
        t = time.time()
        a1[ep] = compare(f"metabric_{ep}", Xm, y, fsets_m, M_COMP, [M_SPAT, M_OE, M_RES, M_CEL], cfg)
        a1[ep]["calibration"] = calibration(f"metabric_{ep}", y, a1[ep]["oof"], cfg)
        print(f"\n[METABRIC {ep}] {time.time() - t:.0f}s\n", a1[ep]["models"].round(4).to_string(index=False))
        print(a1[ep]["diffs"].round(4).to_string(index=False))
        print(pd.DataFrame(a1[ep]["robustness"]).T.round(4).to_string())
        summary[f"metabric_{ep}"] = dict(cindex=a1[ep]["models"].to_dict("records"),
                                         paired_differences=a1[ep]["diffs"].to_dict("records"),
                                         robustness=a1[ep]["robustness"])
    summary["H1"] = h1_verdict(a1["dss"])
    print("\nH1:", summary["H1"])

    # ---- A2 transfer
    Xb, basel_comp, basel_coarse = basel_table()
    assert basel_coarse == coarse_cols
    summary["transfer"] = transfer(Xb, Xm, coarse_cols, cfg)
    for ep in ENDPOINTS:
        print(f"\n[transfer {ep}]\n", pd.DataFrame(summary["transfer"][ep]["models"]).round(4).to_string(index=False))

    # ---- B in Basel (exploratory) and C for Basel part I models
    fs_b = [
        FeatureSet(M_CLIN, CLINICAL_COLS),
        FeatureSet(M_COMP, CLINICAL_COLS, basel_comp),
        FeatureSet(M_RES, CLINICAL_COLS, basel_comp, SPATIAL_COLS_LOGOE,
                   residualize_spatial=True),
    ]
    for ep, col in ENDPOINTS.items():
        y = survival_array(Xb["time"].to_numpy(), Xb[col].to_numpy())
        rb = compare(f"basel_{ep}_residualised", Xb, y, fs_b, M_COMP, [M_RES], cfg)
        summary[f"basel_{ep}_residualised"] = dict(cindex=rb["models"].to_dict("records"),
                                                   paired_differences=rb["diffs"].to_dict("records"),
                                                   robustness=rb["robustness"])
        o = pd.read_csv(RES1 / f"{ep}_oof_risk.csv")
        names = sorted({c.split("|")[0] for c in o.columns if "|seed" in c})
        oof1 = {m: np.stack([o[f"{m}|seed{s}"].to_numpy() for s in cfg["seeds"]]) for m in names}
        calibration(f"basel_{ep}", survival_array(o["time"].to_numpy(), o["event"].to_numpy()), oof1, cfg)

    # ---- D: ER-negative subgroup (exploratory, gated on DSS events)
    er = Xm[Xm["ER"] == 0]
    # Amendment 3(a): drop ER and binary clinical columns whose minority level has < 10 patients
    clin_er = []
    for col in CLINICAL_COLS:
        v = er[col].dropna()
        if col == "ER" or (set(v.unique()) <= {0.0, 1.0} and min((v == 0).sum(), (v == 1).sum()) < 10):
            continue
        clin_er.append(col)
    n_ev = int(er["event_dss"].sum())
    summary["er_negative"] = dict(n=len(er), dss_events=n_ev, run=n_ev >= cfg["er_neg_min_dss_events"],
                                  clinical_columns_used=clin_er,
                                  clinical_columns_dropped=[c for c in CLINICAL_COLS if c not in clin_er])
    if summary["er_negative"]["run"]:
        for ep, col in ENDPOINTS.items():
            y = survival_array(er["time"].to_numpy(), er[col].to_numpy())
            fs_er = [FeatureSet(f.name, clin_er, f.composition, f.spatial) for f in fsets_m[:4]]
            rd = compare(f"metabric_erneg_{ep}", er, y, fs_er, M_COMP, [M_SPAT, M_OE], cfg)
            summary["er_negative"][ep] = dict(cindex=rd["models"].to_dict("records"),
                                              paired_differences=rd["diffs"].to_dict("records"),
                                              robustness=rd["robustness"])

    # ---- figures
    eps = {"OS": "os", "DSS": "dss"}
    spatial = {"z-score": M_SPAT, "log O/E": M_OE, "residualised": M_RES, "CELESTA labels": M_CEL}
    cidx = {k: pd.read_csv(OUT / f"metabric_{v}_cindex.csv") for k, v in eps.items()}
    reps = {k: pd.read_csv(OUT / f"metabric_{v}_cv_repeat_deltas.csv") for k, v in eps.items()}
    nulls = {k: {m: pd.read_csv(OUT / f"metabric_{v}_perm_null_{slug(m)}.csv").filter(like="null_delta_seed").to_numpy().ravel()
                 for m in spatial.values()} for k, v in eps.items()}
    headline_plot(cidx, reps, nulls, [M_CLIN, M_COMP, M_SPAT], spatial, FIG / "part2_metabric_headline.png")
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        t = pd.read_csv(OUT / f"transfer_{ep}_cindex.csv")
        cindex_dotplot(t, FIG / f"part2_transfer_{ep}.png",
                       f"Trained on Basel, tested on METABRIC ({lab}, 180-month horizon)")

    summary["runtime_seconds"] = round(time.time() - t0, 1)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"\npart II done in {summary['runtime_seconds']}s")


if __name__ == "__main__":
    main()
