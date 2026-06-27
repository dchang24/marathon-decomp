"""Across-nu v *magnitude* drift, per slice — companion to p02 (correlations).

Same cohort x gender layout and source (``v_xnu.csv``) as ``p02_vcorr_vs_nu.py``,
but instead of correlations the cells show the absolute change in the race factor
v between two nu fits (in log-time units, i.e. ~fractional time):

    upper triangle = mean_j |v_j(nu_i) - v_j(nu_j)|
    lower triangle = max_j  |v_j(nu_i) - v_j(nu_j)|
    diagonal       = 0

Each nu column is mean-centered first (the model's mean(v)=0 gauge), so the
differences are gauge-invariant. The colour scale is shared and logarithmic (mean
and max live on different magnitudes); cells are annotated with the value in
1e-3 log-units (so "2.0" ≈ a 0.2% time difference). One CV solver (anderson).

One figure ->  results/model_selection/baseline/vdiff_vs_nu.png

Self-contained; no arguments needed (VS Code "Run" works).

Run::

    python scripts/02_model_selection/baseline/p03_vdiff_vs_nu.py
    python scripts/02_model_selection/baseline/p03_vdiff_vs_nu.py --solver cv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"
COHORT_ORDER = ["ALL", "Po10", "WA"]
GENDER_ORDER = ["M", "W", "B"]
_GENDER_LABEL = {"M": "men", "W": "women", "B": "both"}


def _nu_label(nu: float) -> str:
    return "inf" if not np.isfinite(nu) else f"{nu:g}"


def _diff_matrices(V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """mean and max |v_i - v_j| over races, from a centered (races x nus) matrix."""
    n = V.shape[1]
    mean_d = np.zeros((n, n))
    max_d = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = np.abs(V[:, i] - V[:, j])
            mean_d[i, j] = d.mean()
            max_d[i, j] = d.max()
    return mean_d, max_d


def _panel(df_slice: pd.DataFrame, solver: str):
    g = df_slice[df_slice.solver == solver]
    if g.empty:
        g = df_slice[df_slice.solver == sorted(df_slice.solver.unique())[0]]
    nus = sorted(g["nu"].unique(), key=lambda x: (np.isfinite(x), -x))  # inf first
    piv = g.pivot(index="race_id", columns="nu", values="v").reindex(columns=nus).dropna()
    V = piv.to_numpy()
    V = V - V.mean(axis=0, keepdims=True)            # remove the mean(v)=0 gauge
    mean_d, max_d = _diff_matrices(V)
    comb = np.triu(mean_d, 1) + np.tril(max_d, -1)   # upper=mean, lower=max
    return comb, [_nu_label(n) for n in nus]


def _draw(ax, C: np.ndarray, labels: list[str], *, norm, cmap):
    n = len(labels)
    M = np.ma.masked_where(np.eye(n, dtype=bool), C)  # hide the zero diagonal
    im = ax.imshow(M, cmap=cmap, norm=norm)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=7); ax.set_yticklabels(labels, fontsize=7)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            v = C[i, j]
            t = norm(v) if v > 0 else 0.0
            ax.text(j, i, f"{v * 1e3:.1f}", ha="center", va="center", fontsize=6,
                    color="white" if t < 0.5 else "black")
    return im


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver", default="anderson")
    ap.add_argument("--cmap", default="magma_r")
    args = ap.parse_args()

    paths = sorted(OUT_ROOT.glob("*/v_xnu.csv"))
    if not paths:
        print(f"No v_xnu.csv under {OUT_ROOT} — run e01_nu_cv.py first.")
        return

    panels: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for p in paths:
        slug = p.parent.name
        cohort, gender = slug.split("_")[0], slug.split("_")[1]
        panels[(cohort, gender)] = _panel(pd.read_csv(p), args.solver)

    cohorts = [c for c in COHORT_ORDER if any(k[0] == c for k in panels)]
    genders = [g for g in GENDER_ORDER if any(k[1] == g for k in panels)]
    nrow, ncol = len(cohorts), len(genders)

    off = np.concatenate([C[~np.eye(C.shape[0], dtype=bool)] for C, _ in panels.values()])
    off = off[off > 0]
    norm = LogNorm(vmin=float(off.min()), vmax=float(off.max()))
    cmap = plt.get_cmap(args.cmap)

    fig, axes = plt.subplots(nrow, ncol, figsize=(3.8 * ncol, 3.6 * nrow),
                             squeeze=False, constrained_layout=True)
    im = None
    for r, cohort in enumerate(cohorts):
        for c, gender in enumerate(genders):
            ax = axes[r, c]
            key = (cohort, gender)
            if key not in panels:
                ax.set_axis_off()
                continue
            C, labels = panels[key]
            im = _draw(ax, C, labels, norm=norm, cmap=cmap)
            ax.set_title(f"{cohort} — {_GENDER_LABEL.get(gender, gender)}", fontsize=11)
            if r == nrow - 1:
                ax.set_xlabel("nu")
            if c == 0:
                ax.set_ylabel("nu")

    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.6, location="right")
        cb.set_label("|v(nu_i) - v(nu_j)|  (log-time units)")
    fig.suptitle("Baseline v drift across nu  —  upper △ mean |Δv|, lower △ max |Δv| "
                 "(cells ×1e-3 log-units)", fontsize=12)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "vdiff_vs_nu.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
