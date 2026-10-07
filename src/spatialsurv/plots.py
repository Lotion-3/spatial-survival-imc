"""Figures (static PNG). Colours: reference categorical slots 1-2, neutral ink for text."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import KaplanMeierFitter  # noqa: E402
from lifelines.statistics import logrank_test  # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
    }
)


def cindex_dotplot(table: pd.DataFrame, out: Path, title: str) -> None:
    """table: model, cindex, ci_low, ci_high (top-to-bottom in given order)."""
    t = table.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.4, 0.55 * len(t) + 1.2))
    ax.axvline(0.5, color=MUTED, lw=1, ls="--", zorder=0)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    for i, r in t.iterrows():
        c = MUTED if r["model"].startswith("random") else BLUE
        ax.plot([r.ci_low, r.ci_high], [i, i], color=c, lw=2, solid_capstyle="round", zorder=2)
        ax.plot(r.cindex, i, "o", ms=8, color=c, mec=SURFACE, mew=2, zorder=3)
        ax.text(r.ci_high + 0.008, i, f"{r.cindex:.3f}", va="center", color=INK2, fontsize=9)
    ax.set_yticks(range(len(t)), t["model"])
    ax.set_xlabel("Harrell's C-index (pooled out-of-fold, mean of 3 CV repeats; 95% bootstrap CI)")
    lo = min(0.4, t.ci_low.min() - 0.02)
    ax.set_xlim(lo, max(0.85, t.ci_high.max() + 0.06))
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    fig.savefig(out)
    plt.close(fig)


def km_plot(time: np.ndarray, event: np.ndarray, risk: np.ndarray, out: Path, title: str) -> float:
    """Median split of out-of-fold risk; returns log-rank p."""
    high = risk > np.median(risk)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.grid(color=GRID, lw=0.8, zorder=0)
    for mask, lab, c in [(~high, "Low predicted risk", BLUE), (high, "High predicted risk", ORANGE)]:
        kmf = KaplanMeierFitter(label=f"{lab} (n={mask.sum()}, events={int(event[mask].sum())})")
        kmf.fit(time[mask], event[mask])
        kmf.plot_survival_function(ax=ax, color=c, lw=2, ci_alpha=0.12, show_censors=True,
                                   censor_styles={"ms": 5, "marker": "|"})
    p = logrank_test(time[high], time[~high], event[high], event[~high]).p_value
    ax.set_xlabel("Months")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, loc="lower left")
    ax.set_title(f"{title}\nlog-rank p = {p:.2g}", loc="left", color=INK, fontsize=11)
    fig.savefig(out)
    plt.close(fig)
    return float(p)


def spatial_coef_plot(df: pd.DataFrame, out: Path, title: str) -> None:
    """df: feature, mean_coef, sd_coef, selection_freq (sorted)."""
    t = df.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(t) + 1.4))
    ax.axvline(0, color=MUTED, lw=1)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    for i, r in t.iterrows():
        c = MUTED if r.selection_freq == 0 else (ORANGE if r.mean_coef > 0 else BLUE)
        ax.plot([r.mean_coef - r.sd_coef, r.mean_coef + r.sd_coef], [i, i], color=c, lw=2,
                solid_capstyle="round", alpha=0.5)
        ax.plot(r.mean_coef, i, "o", ms=8, color=c, mec=SURFACE, mew=2)
    ax.set_yticks(range(len(t)), [f"{f}  ({s:.0%})" for f, s in zip(t.feature, t.selection_freq)])
    ax.set_xlabel("Standardised log-hazard coefficient (mean ± SD over 15 outer folds)\n"
                  "orange = higher value, higher hazard; grey = never selected; % = folds with non-zero coefficient")
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    fig.savefig(out)
    plt.close(fig)
