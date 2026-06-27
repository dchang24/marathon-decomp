"""nu-selection curves for every slice, laid out as a cohort x gender grid.

Reads the cross-slice rollup written by ``e01_nu_cv.py``
(``results/model_selection/baseline/nu_selection_all.csv``) and draws ONE figure
with a panel per slice:

    columns : gender   (M | W | B)
    rows    : cohort   (ALL | Po10 | WA)

Each panel plots the held-out K-fold CV log predictive density per cell vs nu
(log-x over the finite grid) -- the SELECTION criterion -- with one line per CV
solver, grid points as markers and Brent refinements as small dots. The nu=inf
(Gaussian) value is a dashed horizontal reference, and the selected nu* a dotted
vertical line. Only the cohort/gender combinations actually present in the
rollup get a panel; empty grid cells are blanked.

One figure ->  results/model_selection/baseline/loglik_vs_nu.png

This script is SELF-CONTAINED and needs no arguments.

Run::

    python scripts/02_model_selection/baseline/p01_opt_nu.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"
ROLLUP = OUT_ROOT / "nu_selection_all.csv"

COHORT_ORDER = ["ALL", "Po10", "WA"]   # rows
GENDER_ORDER = ["M", "W", "B"]         # columns
_GENDER_LABEL = {"M": "men", "W": "women", "B": "both"}
_SOLVER_COLOR = {"als": "#4363d8", "anderson": "#e6194b", "cv": "#4363d8"}
_FALLBACK = ["#3cb44b", "#911eb4", "#f58231"]


def _split_slice(name: str) -> tuple[str, str]:
    """'Po10_M' -> ('Po10', 'M'); cohort is everything before the last token."""
    cohort, _, gender = name.rpartition("_")
    return cohort, gender


def _finite(df: pd.DataFrame) -> pd.DataFrame:
    return df[np.isfinite(df["nu"])]


def _solver_color(solver: str, i: int) -> str:
    return _SOLVER_COLOR.get(solver, _FALLBACK[i % len(_FALLBACK)])


def _draw_panel(ax, sel: pd.DataFrame) -> None:
    # nu* dash lines often overlap, so the values go in a small per-panel legend
    # (the dashed vertical lines stay) rather than as on-axis text.
    nu_star: list[Line2D] = []
    for i, (solver, g) in enumerate(sel.groupby("solver")):
        col = _solver_color(solver, i)
        grid = _finite(g[g.source == "grid"]).sort_values("nu")
        brent = _finite(g[g.source == "brent"]).sort_values("nu")
        inf = g[(~np.isfinite(g["nu"])) & (g.source == "grid")]

        ax.plot(grid["nu"], grid["cv_per_cell"], "-o", color=col, ms=4,
                label=f"{solver} (grid)")
        if not brent.empty:
            ax.plot(brent["nu"], brent["cv_per_cell"], ".", color=col, alpha=0.6,
                    label=f"{solver} (brent)")
        if not inf.empty:
            ax.axhline(float(inf["cv_per_cell"].iloc[0]), ls="--", color=col,
                       lw=0.9, alpha=0.7, label=f"{solver} nu=inf")
        # selected nu = max CV logdens across grid + brent
        sel_nu = float(g.loc[g["cv_logdens"].idxmax(), "nu"])
        if np.isfinite(sel_nu):
            ax.axvline(sel_nu, ls=":", color=col, lw=1.2)
            lab = f"{solver}: nu*={sel_nu:.2f}"
        else:
            lab = f"{solver}: nu*=inf"
        nu_star.append(Line2D([], [], ls=":", color=col, lw=1.2, label=lab))
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.25, lw=0.4)
    # per-panel legend with the (possibly overlapping) nu* values, made from
    # proxy handles so they don't leak into the figure-level style legend.
    ax.legend(handles=nu_star, fontsize=8, loc="lower center", framealpha=0.9)


def main() -> None:
    if not ROLLUP.is_file():
        print(f"No rollup at {ROLLUP} — run e01_nu_cv.py first.")
        return
    df = pd.read_csv(ROLLUP)
    df[["cohort", "gender"]] = df["slice"].apply(lambda s: pd.Series(_split_slice(s)))

    cohorts = [c for c in COHORT_ORDER if c in set(df["cohort"])]
    genders = [g for g in GENDER_ORDER if g in set(df["gender"])]
    nrow, ncol = len(cohorts), len(genders)

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             squeeze=False, constrained_layout=True)
    legend_ax = None
    for r, cohort in enumerate(cohorts):
        for c, gender in enumerate(genders):
            ax = axes[r, c]
            sel = df[(df.cohort == cohort) & (df.gender == gender)]
            if sel.empty:
                ax.set_axis_off()
                continue
            _draw_panel(ax, sel)
            ax.set_title(f"{cohort} — {_GENDER_LABEL.get(gender, gender)}", fontsize=11)
            if c == 0:
                ax.set_ylabel("held-out CV logdens / cell")
            if r == nrow - 1:
                ax.set_xlabel("nu  (log; dashed = nu=inf Gaussian)")
            legend_ax = ax

    if legend_ax is not None:
        handles, labels = legend_ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside upper right", ncol=2, fontsize=8,
                   frameon=True)
    fig.suptitle("Baseline nu selection — held-out CV log predictive density vs nu",
                 fontsize=13)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "loglik_vs_nu.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
