"""Re-type the METABRIC IMC cells with CELESTA-Lite and compare against the published
phenotypes (Danenberg et al. 2022), then see how the spatial features change.

No survival outcomes are used here, so the pre-registered part II analysis
(ANALYSIS_PLAN.md) is unaffected.

    uv run python scripts/run_celesta_metabric.py          # both signature versions
    uv run python scripts/run_celesta_metabric.py v2       # one version
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.metrics import cohen_kappa_score, confusion_matrix  # noqa: E402

from celesta_lite import (  # noqa: E402
    UNKNOWN,
    CelestaConfig,
    Signature,
    arcsinh_transform,
    expression_probabilities,
    fit_marker_models,
    run_celesta,
)
from spatialsurv.features import SPATIAL_COLS, SPATIAL_COLS_LOGOE, spatial_features_by_image  # noqa: E402
from spatialsurv.metabric import PHENOTYPE_CLASS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw_metabric"
OUT = ROOT / "results" / "celesta"
FIG = ROOT / "figures"

CONFIG = dict(
    cofactor=1.0,  # arcsinh(x / 1): IMC mean ion counts are mostly < 5
    ep_method="posterior",  # see celesta_lite.expression_probabilities: the reference slope-1
    # sigmoid is nearly flat on these dim channels (CD3 spans ~0.2 arcsinh units)
    gmm_subsample=200_000,  # cells used to fit the cohort-level marker mixtures
    seed=0,
    celesta=dict(n_neighbors=5, bandwidth=100.0, scale_factor=5.0, max_iteration=10, cell_change_threshold=0.01),
    spatial=dict(knn_k=8, n_perm=200, seed=12345),
)

# Prior knowledge (cell type x marker): 1 = high, 0 = low, NaN = uninformative. Written from
# marker biology, not tuned on the published labels it is compared against.
NA = np.nan
_COLS = ["panCK", "CD3", "CD20", "CD68", "CD31-vWF", "SMA", "PDGFRB", "FSP1", "Podoplanin"]


def _signature(rows: dict) -> pd.DataFrame:
    t = pd.DataFrame.from_dict(rows, orient="index", columns=_COLS)
    t = t.loc[:, t.notna().any()]  # drop markers no type uses
    t.insert(0, "parent", "")
    t.insert(0, "round", 1)
    return t


#              panCK CD3 CD20 CD68 CD31 SMA PDGFRB FSP1 PDPN
_V1 = {
    "Tumor":          [1, 0, 0, 0, 0, NA, NA, NA, NA],
    "T cell":         [0, 1, 0, 0, NA, NA, NA, NA, NA],
    "B cell":         [0, 0, 1, NA, NA, NA, NA, NA, NA],
    "Macrophage":     [0, 0, 0, 1, NA, NA, NA, NA, NA],
    "Endothelial":    [0, 0, NA, 0, 1, NA, NA, NA, NA],
    "Myofibroblast":  [0, 0, NA, 0, 0, 1, NA, NA, NA],
    "Fibroblast":     [0, 0, NA, 0, 0, NA, 1, NA, NA],
}
# v2 (revised after v1 results, kept alongside v1): adds the two canonical cancer-associated
# fibroblast markers in the panel that v1 omitted. v1 called 35% of published fibroblasts
# Unknown and sent FSP1+ fibroblasts to Macrophage/B. FSP1 (S100A4) is also expressed by
# macrophages, hence CD68 low for that type.
_V2 = _V1 | {
    "FSP1+ fibroblast": [0, 0, NA, 0, 0, NA, NA, 1, NA],
    "PDPN+ fibroblast": [0, 0, NA, 0, 0, NA, NA, NA, 1],
}
SIGNATURES = {"v1": _signature(_V1), "v2": _signature(_V2)}
SIGNATURE = SIGNATURES["v2"]  # superset of markers, used for loading
TO_COARSE = {"Tumor": "Tumor", "T cell": "T", "B cell": "B", "Macrophage": "Macrophage",
             "Endothelial": "Endothelial", "Myofibroblast": "Stroma", "Fibroblast": "Stroma",
             "FSP1+ fibroblast": "Stroma", "PDPN+ fibroblast": "Stroma", UNKNOWN: UNKNOWN}
COARSE_ORDER = ["Tumor", "T", "B", "Macrophage", "Endothelial", "Stroma"]


def load_cells() -> pd.DataFrame:
    markers = [c for c in SIGNATURE.columns if c not in ("round", "parent")]
    cols = ["ImageNumber", "metabric_id", "cellPhenotype", "is_epithelial", "is_tumour", "is_hotAggregate",
            "Location_Center_X", "Location_Center_Y"] + markers
    c = pd.read_csv(RAW / "SingleCells.csv", usecols=cols)
    c = c[(c["is_tumour"] == 1) & (c["is_hotAggregate"] == 0)].reset_index(drop=True)
    c["published"] = np.where(c["is_epithelial"] == 1, "Tumor", c["cellPhenotype"].map(PHENOTYPE_CLASS))
    c["core"] = "MB_img" + c["ImageNumber"].astype(str)
    return c


def run_image(core: str, ep: pd.DataFrame, xy: np.ndarray, sig: Signature, cfg: CelestaConfig) -> pd.DataFrame:
    r = run_celesta(ep, xy, sig, cfg)
    return pd.DataFrame({"celesta": r.labels, "marker_only": r.marker_only, "anchor": r.is_anchor,
                         "artifact": r.is_artifact}, index=ep.index)


def agreement(truth: pd.Series, pred: pd.Series) -> dict:
    assigned = pred != UNKNOWN
    t, p = truth[assigned], pred[assigned]
    return dict(coverage=float(assigned.mean()), n_assigned=int(assigned.sum()),
                agreement_assigned=float((t == p).mean()),
                kappa_assigned=float(cohen_kappa_score(t, p, labels=COARSE_ORDER)),
                agreement_all_cells=float(((truth == pred) & assigned).mean()))


def per_class(truth: pd.Series, pred: pd.Series) -> pd.DataFrame:
    rows = []
    for k in COARSE_ORDER:
        tp = int(((truth == k) & (pred == k)).sum())
        rows.append(dict(cls=k, published_n=int((truth == k).sum()), celesta_n=int((pred == k).sum()),
                         recall=tp / max(int((truth == k).sum()), 1), precision=tp / max(int((pred == k).sum()), 1)))
    return pd.DataFrame(rows)


def example_figure(c: pd.DataFrame, core: str, out: Path) -> None:
    from spatialsurv.plots import CELL_COLOURS, INK, SURFACE  # same colours as the part I tissue figure

    d = c[c["core"] == core]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8))
    for ax, col, title in zip(axes, ["published", "marker_only_coarse", "celesta_coarse"],
                              ["Published phenotypes", "CELESTA-Lite, markers only", "CELESTA-Lite, markers + spatial MRF"]):
        for cls, colour in list(CELL_COLOURS.items()) + [(UNKNOWN, "#ffffff")]:
            m = d[col] == cls
            ax.scatter(d.loc[m, "Location_Center_X"], d.loc[m, "Location_Center_Y"], s=4, c=colour,
                       linewidths=0.3 if cls == UNKNOWN else 0, edgecolors="#8a8984" if cls == UNKNOWN else None,
                       rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
    handles = [plt.Line2D([], [], ls="", marker="o", ms=7, color=v, label=k) for k, v in CELL_COLOURS.items()]
    handles.append(plt.Line2D([], [], ls="", marker="o", ms=7, mfc="#ffffff", mec="#8a8984", label="Unknown"))
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False)
    fig.suptitle(f"{core}: cell typing compared", x=0.01, ha="left", fontsize=11, color=INK)
    fig.set_facecolor(SURFACE)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main(version: str) -> None:
    t0 = time.time()
    out = OUT / version
    out.mkdir(parents=True, exist_ok=True)
    cfg = CONFIG
    sig = Signature(SIGNATURES[version])
    c = load_cells()
    print(f"{len(c)} cells, {c['core'].nunique()} images ({time.time() - t0:.0f}s)")

    # Cohort-level marker mixtures (deviation from the reference's per-sample fit: ~1,400-cell
    # TMA cores often lack a cell type entirely, and a 2-component fit to a unimodal marker
    # distribution invents a spurious 'positive' population).
    X = arcsinh_transform(c[sig.markers], cfg["cofactor"])
    sub = X.sample(n=min(cfg["gmm_subsample"], len(X)), random_state=cfg["seed"])
    models = fit_marker_models(sub)
    ep = expression_probabilities(X, models, method=cfg["ep_method"])
    pd.DataFrame([dict(marker=m.marker, mean_low=m.means[0], mean_high=m.means[1], var_low=m.variances[0],
                       var_high=m.variances[1], critical_point=m.critical_point,
                       frac_cells_above=float((X[m.marker] > m.critical_point).mean())) for m in models.values()]
                 ).to_csv(out / "marker_models.csv", index=False)

    ccfg = CelestaConfig(**cfg["celesta"])
    groups = c.groupby("core").indices
    parts = Parallel(n_jobs=-1)(
        delayed(run_image)(core, ep.iloc[idx], c[["Location_Center_X", "Location_Center_Y"]].to_numpy()[idx], sig, ccfg)
        for core, idx in groups.items()
    )
    res = pd.concat(parts).sort_index()
    c = c.join(res)
    c["celesta_coarse"] = c["celesta"].map(TO_COARSE)
    c["marker_only_coarse"] = c["marker_only"].map(TO_COARSE)
    print(f"CELESTA done ({time.time() - t0:.0f}s)")

    # ---- agreement with the published phenotypes
    summary = {"config": cfg, "signature_version": version, "signature": SIGNATURES[version].replace({np.nan: None}).to_dict(orient="index"),
               "n_cells": int(len(c)), "n_images": int(c["core"].nunique()),
               "artifact_frac": float(c["artifact"].mean()), "anchor_frac": float(c["anchor"].mean())}
    summary["celesta_vs_published"] = agreement(c["published"], c["celesta_coarse"])
    summary["marker_only_vs_published"] = agreement(c["published"], c["marker_only_coarse"])
    # Cells assigned by the spatial MRF step (not anchors): does the spatial term help on them?
    idx_cells = (c["celesta"] != UNKNOWN) & ~c["anchor"]
    summary["index_cells"] = dict(
        n=int(idx_cells.sum()),
        celesta_agreement=float((c.loc[idx_cells, "celesta_coarse"] == c.loc[idx_cells, "published"]).mean()),
        marker_only_agreement_same_cells=float((c.loc[idx_cells, "marker_only_coarse"] == c.loc[idx_cells, "published"]).mean()),
        changed_from_marker_only=float((c.loc[idx_cells, "celesta"] != c.loc[idx_cells, "marker_only"]).mean()),
    )
    anchors = c["anchor"]
    summary["anchor_cells_agreement"] = float((c.loc[anchors, "celesta_coarse"] == c.loc[anchors, "published"]).mean())
    per_class(c["published"], c["celesta_coarse"]).to_csv(out / "per_class_celesta.csv", index=False)
    per_class(c["published"], c["marker_only_coarse"]).to_csv(out / "per_class_marker_only.csv", index=False)
    cm = confusion_matrix(c["published"], c["celesta_coarse"], labels=COARSE_ORDER + [UNKNOWN])
    pd.DataFrame(cm, index=[f"published {k}" for k in COARSE_ORDER + [UNKNOWN]],
                 columns=[f"celesta {k}" for k in COARSE_ORDER + [UNKNOWN]]).to_csv(out / "confusion_celesta.csv")
    c[["core", "metabric_id", "published", "celesta", "celesta_coarse", "marker_only_coarse", "anchor", "artifact"]] \
        .to_csv(out / "cell_labels.csv.gz", index=False, compression="gzip")
    print(json.dumps({k: summary[k] for k in ["celesta_vs_published", "marker_only_vs_published", "index_cells",
                                                "anchor_cells_agreement"]}, indent=2))

    # ---- downstream: spatial features from CELESTA labels vs. published labels (per image)
    sp = cfg["spatial"]
    base = c.rename(columns={"Location_Center_X": "x", "Location_Center_Y": "y"})
    f_pub = spatial_features_by_image(base.assign(coarse=base["published"]), sp["knn_k"], sp["n_perm"], sp["seed"])
    known = base[base["celesta_coarse"] != UNKNOWN]  # Unknown cells are dropped from the graph
    f_cel = spatial_features_by_image(known.assign(coarse=known["celesta_coarse"]), sp["knn_k"], sp["n_perm"], sp["seed"])
    f_pub.to_csv(out / "spatial_features_published.csv")
    f_cel.to_csv(out / "spatial_features_celesta.csv")
    rows = []
    for col in SPATIAL_COLS + [x for x in SPATIAL_COLS_LOGOE if x not in SPATIAL_COLS]:
        j = pd.concat([f_pub[col], f_cel[col]], axis=1, keys=["pub", "cel"]).dropna()
        rows.append(dict(feature=col, n_images=len(j),
                         spearman=float(j["pub"].corr(j["cel"], method="spearman")) if len(j) > 2 else np.nan))
    feat = pd.DataFrame(rows)
    feat.to_csv(out / "spatial_feature_concordance.csv", index=False)
    summary["spatial_feature_concordance"] = feat.round(3).to_dict(orient="records")
    print(feat.round(3).to_string(index=False))

    # ---- example image: the one closest to median CELESTA-vs-published agreement
    per_img = c.groupby("core").apply(lambda d: (d["celesta_coarse"] == d["published"]).mean(), include_groups=False)
    core = (per_img - per_img.median()).abs().idxmin()
    example_figure(c, core, FIG / f"celesta_example_core_{version}.png")
    summary["example_core"] = dict(core=core, agreement=float(per_img[core]))
    summary["per_image_agreement"] = dict(median=float(per_img.median()), q10=float(per_img.quantile(0.1)),
                                          q90=float(per_img.quantile(0.9)))
    summary["runtime_seconds"] = round(time.time() - t0, 1)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"done in {summary['runtime_seconds']}s")


if __name__ == "__main__":
    import sys

    for v in (sys.argv[1:] or list(SIGNATURES)):
        main(v)
