"""Write results/part2/tables.md from saved part II results. No numbers are typed by hand.

Usage: uv run python scripts/report_part2.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from report_tables import cindex_table, diff_table, robustness_table, signed

OUT = Path(__file__).resolve().parents[1] / "results" / "part2"


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
    out += ["", "## D ER-negative subgroup", "", f"- n: {e['n']}", f"- dss_events: {e['dss_events']}", f"- run: {e['run']}"]
    if e["run"]:
        for ep, lab in [("os", "OS"), ("dss", "DSS")]:
            n = f"metabric_erneg_{ep}"
            out += ["", f"### ER-negative {lab}", "", cindex_table(n, OUT), "", diff_table(n, OUT), "",
                    robustness_table(n, OUT)]
    out += ["", f"runtime_seconds: {s['runtime_seconds']}"]
    (OUT / "tables.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'tables.md'}")


if __name__ == "__main__":
    main()
