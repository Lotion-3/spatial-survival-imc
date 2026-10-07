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


CV_XLABEL = "Harrell's C-index (pooled out-of-fold, mean of 3 CV repeats; 95% bootstrap CI)"


def cindex_dotplot(table: pd.DataFrame, out: Path, title: str, xlabel: str = CV_XLABEL) -> None:
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
    ax.set_xlabel(xlabel)
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


AQUA = "#1baf7a"
# Tumor/stroma/endothelium are context (neutral greys); immune classes carry colour.
CELL_COLOURS = {
    "Stroma": "#e6e4de",
    "Tumor": "#a8a59c",
    "Endothelial": "#6f6e69",
    "Macrophage": BLUE,
    "T": ORANGE,
    "B": AQUA,
}


def example_cores_plot(cells: pd.DataFrame, picks: list[tuple[str, str]], out: Path, title: str) -> None:
    """picks: [(core, subtitle)]; cells: core, x, y, coarse."""
    fig, axes = plt.subplots(1, len(picks), figsize=(4.2 * len(picks), 4.9))
    for ax, (core, sub) in zip(np.atleast_1d(axes), picks):
        d = cells[cells["core"] == core]
        for cls, col in CELL_COLOURS.items():  # draw context first, immune cells on top
            m = d["coarse"] == cls
            ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=3.5 if cls in ("Stroma", "Tumor") else 5,
                       c=col, linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()  # image coordinates
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        x0, y1 = d["x"].min(), d["y"].max()
        ax.plot([x0, x0 + 100], [y1 + 25, y1 + 25], color=INK2, lw=2)
        ax.text(x0 + 50, y1 + 35, "100 µm", ha="center", va="top", color=INK2, fontsize=8)
        ax.set_title(sub, loc="left", color=INK, fontsize=9.5)
    handles = [plt.Line2D([], [], ls="", marker="o", ms=7, color=c, label=k) for k, c in CELL_COLOURS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=11)
    fig.subplots_adjust(top=0.80, wspace=0.08)
    fig.savefig(out)
    plt.close(fig)


def headline_plot(cidx: dict[str, pd.DataFrame], reps: dict[str, pd.DataFrame], nulls: dict[str, dict[str, np.ndarray]],
                  models: list[str], spatial_models: dict[str, str], out: Path) -> None:
    """Left: C-index with 95% CI per endpoint for the 3 main models.
    Right: spatial minus composition delta per CV repeat (dots) vs. permuted-spatial null (grey)."""
    eps = list(cidx)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 3.9), gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.55})
    short = {models[0]: "clinical", models[1]: "+ composition", models[2]: "+ composition + spatial"}
    # left panel
    yt, yl = [], []
    for e, ep in enumerate(eps):
        t = cidx[ep].set_index("model")
        for j, m in enumerate(models):
            y = e * (len(models) + 1.2) + j
            r = t.loc[m]
            a1.plot([r.ci_low, r.ci_high], [y, y], color=BLUE, lw=2, solid_capstyle="round")
            a1.plot(r.cindex, y, "o", ms=7, color=BLUE, mec=SURFACE, mew=1.5)
            a1.text(r.ci_high + 0.006, y, f"{r.cindex:.3f}", va="center", fontsize=8.5, color=INK2)
            yt.append(y)
            yl.append(f"{ep}: {short[m]}")
    a1.set_yticks(yt, yl)
    a1.invert_yaxis()
    a1.set_xlim(0.58, 0.85)
    a1.grid(axis="x", color=GRID, lw=0.8)
    a1.set_xlabel("Out-of-fold Harrell's C (95% patient-bootstrap CI)")
    a1.set_title("A  Discrimination by feature set", loc="left", color=INK, fontsize=10.5)
    # right panel
    rng = np.random.default_rng(0)
    yt, yl = [], []
    y = 0
    for ep in eps:
        for lab, m in spatial_models.items():
            nd = nulls[ep][m]
            a2.fill_betweenx([y - 0.32, y + 0.32], np.quantile(nd, 0.025), np.quantile(nd, 0.975), color=GRID, lw=0)
            a2.plot([nd.mean()] * 2, [y - 0.32, y + 0.32], color=MUTED, lw=1.5)
            d = reps[ep].loc[reps[ep]["model"] == m, "delta"].to_numpy()
            a2.plot(d, y + rng.uniform(-0.18, 0.18, len(d)), "o", ms=5, color=ORANGE, mec=SURFACE, mew=1, alpha=0.9)
            a2.plot([d.mean()] * 2, [y - 0.3, y + 0.3], color=ORANGE, lw=2.5)
            yt.append(y)
            yl.append(f"{ep}: {lab}")
            y += 1
        y += 0.5
    a2.axvline(0, color=INK2, lw=1)
    a2.set_yticks(yt, yl)
    a2.invert_yaxis()
    a2.grid(axis="x", color=GRID, lw=0.8)
    a2.set_xlabel("ΔC, spatial model minus clinical+composition")
    a2.set_title("B  Added value of spatial features", loc="left", color=INK, fontsize=10.5)
    handles = [plt.Line2D([], [], ls="", marker="o", ms=6, color=ORANGE, label="each of 10 CV repeats (bar = mean)"),
               plt.Rectangle((0, 0), 1, 1, color=GRID, label="null: spatial block shuffled across patients (95% range, grey bar = mean)")]
    a2.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.2), frameon=False, fontsize=8.5)
    fig.savefig(out)
    plt.close(fig)
