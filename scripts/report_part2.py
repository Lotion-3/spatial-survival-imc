"""Write results/part2/tables.md from saved part II results. No numbers are typed by hand.

Usage: uv run python scripts/report_part2.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from report_tables import cindex_table, diff_table, robustness_table, signed

OUT = Path(__file__).resolve().parents[1] / "results" / "part2"
CEL = Path(__file__).resolve().parents[1] / "results" / "celesta"


def celesta_tables() -> list[str]:
    rows = ["| Signature | Labelling | cells assigned | agreement (assigned cells) | Cohen's κ (assigned cells) |",
            "|---|---|---|---|---|"]
    idx = ["| Signature | non-anchor cells assigned by the MRF | agreement: CELESTA-Lite | agreement: markers only, same cells |",
           "|---|---|---|---|"]
    for v in ["v1", "v2"]:
        s = json.loads((CEL / v / "summary.json").read_text())
        for key, lab in [("marker_only_vs_published", "markers only (argmax score)"), ("celesta_vs_published", "CELESTA-Lite (markers + spatial MRF)")]:
            a = s[key]
            rows.append(f"| {v} | {lab} | {a['coverage']:.1%} | {a['agreement_assigned']:.3f} | {a['kappa_assigned']:.3f} |")
        i = s["index_cells"]
        idx.append(f"| {v} | {i['n']:,} | {i['celesta_agreement']:.3f} | {i['marker_only_agreement_same_cells']:.3f} |")
    pc = pd.read_csv(CEL / "v1" / "per_class_celesta.csv")
    cls = ["| Class (v1) | published cells | CELESTA-Lite cells | recall | precision |", "|---|---|---|---|---|"]
    for _, r in pc.iterrows():
        cls.append(f"| {r.cls} | {r.published_n:,} | {r.celesta_n:,} | {r.recall:.2f} | {r.precision:.2f} |")
    fc = pd.read_csv(CEL / "v1" / "spatial_feature_concordance.csv")
    feat = ["| Spatial feature (v1 labels) | images | Spearman ρ, CELESTA-Lite vs. published labels |", "|---|---|---|"]
    for _, r in fc.iterrows():
        feat.append(f"| `{r.feature}` | {r.n_images} | {r.spearman:+.2f} |")
    nl = "\n"
    return ["## CELESTA-Lite vs published labels", "", nl.join(rows), "", nl.join(idx), "", nl.join(cls), "",
            nl.join(feat), "", "## CELESTA-Lite: reference vs posterior expression probabilities (v1 signature)", "",
            ep_method_table()]


def ep_method_table() -> str:
    lines = ["| EP method | artifact-filtered | cells assigned | agreement (assigned) | κ (assigned) | CELESTA T cells | T recall "
             "| mean CD3 EP, published T | mean CD3 EP, published Tumor |",
             "|---|---|---|---|---|---|---|---|---|"]
    for run, lab in [("v1_reference_ep", "reference (sigmoid of x − crossing point)"), ("v1", "posterior (this port)")]:
        s = json.loads((CEL / run / "summary.json").read_text())
        a = s["celesta_vs_published"]
        pc = pd.read_csv(CEL / run / "per_class_celesta.csv").set_index("cls")
        ep = pd.read_csv(CEL / run / "mean_ep_by_published_class.csv").set_index("published")
        lines.append(f"| {lab} | {s['artifact_frac']:.1%} | {a['coverage']:.1%} | {a['agreement_assigned']:.3f} | "
                     f"{a['kappa_assigned']:.3f} | {pc.loc['T', 'celesta_n']:,} | {pc.loc['T', 'recall']:.2f} | "
                     f"{ep.loc['T', 'CD3']:.2f} | {ep.loc['Tumor', 'CD3']:.2f} |")
    return "\n".join(lines)


def transfer_table(ep: str) -> str:
    t = pd.read_csv(OUT / f"transfer_{ep}_cindex.csv")
    lines = ["| Model (trained on Basel) | external C-index on METABRIC | 95% CI |", "|---|---|---|"]
    for _, r in t.iterrows():
        lines.append(f"| {r.model} | {r.cindex:.3f} | {r.ci_low:.3f} to {r.ci_high:.3f} |")
    return "\n".join(lines)


def transfer_diff_table(ep: str) -> str:
    t = pd.read_csv(OUT / f"transfer_{ep}_paired_differences.csv")
    lines = ["| Transfer comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |", "|---|---|---|---|"]
    for _, r in t.iterrows():
        lines.append(f"| {r.comparison} | {signed(r.delta)} | {signed(r.ci_low)} to {signed(r.ci_high)} | {r.frac_boot_gt0:.2f} |")
    return "\n".join(lines)


def calib_table(name: str) -> str:
    t = pd.read_csv(OUT / f"{name}_calibration_auc.csv")
    wide = {}
    for _, r in t.iterrows():
        wide.setdefault(r.model, {})[r.metric] = f"{r.value:.2f} ({r.ci_low:.2f} to {r.ci_high:.2f})"
    metrics = ["calibration_slope", "auc_60m", "auc_120m"]
    lines = ["| Model | calibration slope (95% CI) | AUC at 60 months | AUC at 120 months |", "|---|---|---|---|"]
    for m, v in wide.items():
        lines.append(f"| {m} | " + " | ".join(v[k] for k in metrics) + " |")
    return "\n".join(lines)


def main() -> None:
    s = json.loads((OUT / "summary.json").read_text())
    out = ["## METABRIC cohort", ""] + [f"- {k}: {v}" for k, v in s["metabric_cohort"].items()]
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        n = f"metabric_{ep}"
        out += ["", f"## A1 METABRIC {lab}", "", cindex_table(n, OUT), "", diff_table(n, OUT), "",
                robustness_table(n, OUT), "", calib_table(n)]
    h = s["H1"]
    out += ["", "## H1 (confirmatory)", "",
            "| ΔC (log O/E minus composition, DSS) | 95% CI | permutation p | mean ΔC over 10 CV repeats | (a) CI > 0 | (b) p ≤ 0.05 | (c) repeat mean > 0 | verdict |",
            "|---|---|---|---|---|---|---|---|",
            f"| {signed(h['delta'])} | {signed(h['ci_low'])} to {signed(h['ci_high'])} | {h['perm_p']:.3f} | "
            f"{signed(h['cv_repeat_mean'])} | {'yes' if h['criterion_a_ci_above_0'] else 'no'} | "
            f"{'yes' if h['criterion_b_perm_p_le_005'] else 'no'} | {'yes' if h['criterion_c_repeat_mean_gt0'] else 'no'} | "
            f"**{h['verdict']}** |"]
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        out += ["", f"## A2 transfer {lab}", "", transfer_table(ep), "", transfer_diff_table(ep)]
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        n = f"basel_{ep}_residualised"
        out += ["", f"## B Basel residualised {lab} (exploratory)", "", cindex_table(n, OUT), "",
                diff_table(n, OUT), "", robustness_table(n, OUT)]
    for ep, lab in [("os", "OS"), ("dss", "DSS")]:
        out += ["", f"## C Basel calibration {lab}", "", calib_table(f"basel_{ep}")]
    e = s["er_negative"]
    out += ["", "## D ER-negative subgroup", "", f"- n: {e['n']}", f"- dss_events: {e['dss_events']}", f"- run: {e['run']}",
            f"- clinical_columns_dropped: {e['clinical_columns_dropped']}"]
    if e["run"]:
        for ep, lab in [("os", "OS"), ("dss", "DSS")]:
            n = f"metabric_erneg_{ep}"
            out += ["", f"### ER-negative {lab}", "", cindex_table(n, OUT), "", diff_table(n, OUT), "",
                    robustness_table(n, OUT)]
    out += [""] + celesta_tables()
    out += ["", f"runtime_seconds: {s['runtime_seconds']}"]
    (OUT / "tables.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'tables.md'}")


if __name__ == "__main__":
    main()
