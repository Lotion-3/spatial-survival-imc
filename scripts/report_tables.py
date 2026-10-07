"""Print the README's results tables (markdown) from results/. No numbers are typed by hand.

Usage: uv run python scripts/report_tables.py   (writes results/tables.md)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

RES = Path(__file__).resolve().parents[1] / "results"


def f3(x: float) -> str:
    return f"{x:.3f}"


def signed(x: float) -> str:
    return f"{x:+.3f}"


def cindex_table(ep: str, d: Path = RES) -> str:
    t = pd.read_csv(d / f"{ep}_cindex.csv")
    lines = ["| Model | C-index | 95% CI | per CV repeat (seeds 0, 1, 2) |", "|---|---|---|---|"]
    for _, r in t.iterrows():
        seeds = ", ".join(f"{float(v):.3f}" for v in str(r.cindex_per_seed).split(";"))
        lines.append(f"| {r.model} | {f3(r.cindex)} | {f3(r.ci_low)} to {f3(r.ci_high)} | {seeds} |")
    return "\n".join(lines)


def diff_table(ep: str, d: Path = RES) -> str:
    t = pd.read_csv(d / f"{ep}_paired_differences.csv")
    lines = ["| Comparison | ΔC | 95% CI | share of bootstrap Δ > 0 |", "|---|---|---|---|"]
    for _, r in t.iterrows():
        lines.append(f"| {r.comparison} | {signed(r.delta)} | {signed(r.ci_low)} to {signed(r.ci_high)} | {r.frac_boot_gt0:.2f} |")
    return "\n".join(lines)


def robustness_table(ep: str, d: Path = RES) -> str:
    t = pd.read_csv(d / f"{ep}_robustness.csv", index_col=0)
    lines = [
        "| Spatial model vs. clinical+composition | ΔC (3 main repeats) | ΔC over 10 CV repeats: mean ± SD [min, max] | repeats with Δ > 0 | permuted-spatial null ΔC: mean ± SD (max) | permutation p |",
        "|---|---|---|---|---|---|",
    ]
    for m, r in t.iterrows():
        lines.append(
            f"| {m} | {signed(r.observed_delta_3seeds)} | {signed(r.delta_mean)} ± {r.delta_sd:.3f} "
            f"[{signed(r.delta_min)}, {signed(r.delta_max)}] | {int(round(r.frac_repeats_gt0 * r.cv_repeats))}/{int(r.cv_repeats)} | "
            f"{signed(r.null_mean)} ± {r.null_sd:.3f} ({signed(r.null_max)}) | {r.perm_p:.3f} |"
        )
    return "\n".join(lines)


def coef_table(ep: str) -> str:
    t = pd.read_csv(RES / f"{ep}_spatial_coefficients.csv")
    lines = ["| Spatial feature | mean coef | SD | selected in folds |", "|---|---|---|---|"]
    for _, r in t.iterrows():
        lines.append(f"| `{r.feature}` | {signed(r.mean_coef)} | {r.sd_coef:.3f} | {r.selection_freq:.0%} |")
    return "\n".join(lines)


def sparsity_table() -> str:
    lines = ["| Endpoint | Model | selected penalty α: median [min, max] | non-zero omics coefs: median [min, max] |", "|---|---|---|---|"]
    for ep in ["os", "dss"]:
        t = pd.read_csv(RES / f"{ep}_per_fold.csv")
        for m, g in t.groupby("model", sort=False):
            a, nz = g.alpha_selected, g.n_nonzero_omics
            lines.append(f"| {ep.upper()} | {m} | {a.median():.3g} [{a.min():.3g}, {a.max():.3g}] | "
                         f"{nz.median():.0f} [{nz.min()}, {nz.max()}] |")
    return "\n".join(lines)


def confound_table() -> str:
    """Image-level Spearman correlation of each spatial feature with image size (log cell
    count) and with immune-cell fraction (a composition quantity), for both enrichment versions."""
    t = pd.read_csv(RES / "spatial_vs_ncells.csv", index_col=0).iloc[:, 0]
    f = pd.read_csv(RES / "spatial_vs_immune_fraction.csv", index_col=0).iloc[:, 0]
    z = [c for c in t.index if not c.startswith("logoe_")]
    lines = ["| Feature | ρ with log(cells): z-score | ρ with log(cells): log O/E | ρ with immune fraction: z-score | ρ with immune fraction: log O/E |",
             "|---|---|---|---|---|"]
    for c in z:
        oe = c.replace("enrich_", "logoe_")
        if oe != c:
            lines.append(f"| `{c}` | {t[c]:+.2f} | {t[oe]:+.2f} | {f[c]:+.2f} | {f[oe]:+.2f} |")
        else:
            lines.append(f"| `{c}` | {t[c]:+.2f} | (same feature) | {f[c]:+.2f} | (same feature) |")
    return "\n".join(lines)


def mde_table() -> str:
    """Approximate minimal detectable difference: (z_0.975 + z_0.80) * SD of the paired
    bootstrap ΔC, i.e. the true gain this cohort would detect with ~80% power at two-sided
    alpha = 0.05. Uses patient-sampling variability only (OOF predictions held fixed), so it
    is optimistic; the 'with CV-repeat variability' column adds the across-repeat SD of the
    spatial delta in quadrature, scaled to a 3-repeat average."""
    k = 1.959964 + 0.841621
    lines = ["| Endpoint | Comparison | SD of paired bootstrap ΔC | detectable ΔC (80% power) | with CV-repeat variability |",
             "|---|---|---|---|---|"]
    pairs = [("clinical+composition", "clinical"), ("clinical+composition+spatial", "clinical+composition")]
    for ep in ["os", "dss"]:
        b = pd.read_csv(RES / f"{ep}_bootstrap.csv")
        rep = pd.read_csv(RES / f"{ep}_cv_repeat_deltas.csv")
        for a, c in pairs:
            sd = (b[a] - b[c]).std(ddof=1)
            extra = "n/a"
            if a == "clinical+composition+spatial":
                sd_rep = rep.loc[rep["model"] == a, "delta"].std(ddof=1) / np.sqrt(3)
                extra = f"{k * np.sqrt(sd**2 + sd_rep**2):.3f}"
            lines.append(f"| {ep.upper()} | {a} minus {c} | {sd:.4f} | {k * sd:.3f} | {extra} |")
    return "\n".join(lines)


def main() -> None:
    s = json.loads((RES / "summary.json").read_text())
    nl = "\n"
    out = ["## Cohort", ""]
    out += [f"- {k}: {v}" for k, v in s["cohort"].items()]
    for ep, name in [("os", "Overall survival (primary)"), ("dss", "Breast-cancer-specific survival (sensitivity)")]:
        out += ["", f"## {name}", "", cindex_table(ep), "", diff_table(ep), "", robustness_table(ep), "", coef_table(ep)]
    out += ["", "## Detectable effect size", "", mde_table()]
    out += ["", "## Penalty and sparsity", "", sparsity_table()]
    out += ["", "## Spatial features vs. image size and immune content", "", confound_table()]
    out += ["", "## Kaplan-Meier", "",
            f"- best model: {s['km']['model']}, log-rank p = {s['km']['logrank_p']:.2g}",
            f"- clinical-only model: log-rank p = {s['km_clinical']['logrank_p']:.2g}",
            f"- patients assigned to the same risk half by both: {s['km_group_agreement']:.0%}",
            "", f"runtime_seconds: {s['runtime_seconds']}"]
    (RES / "tables.md").write_text(nl.join(out) + nl, encoding="utf-8")
    print(f"wrote {RES / 'tables.md'}")


if __name__ == "__main__":
    main()
