"""How stable is the fitted race factor v across the nu grid? — per slice.

Reads the full-data fitted v per nu (``v_xnu.csv``) for each slice and draws,
per slice, a heatmap of the across-nu agreement of v: cell (i, j) is the
correlation between v fitted at nu_i and at nu_j (same race set, mean-centered).
The matrix is split on the diagonal — **upper triangle = Pearson, lower triangle
= Spearman** — so both rank- and value-agreement are visible in one panel. The
diagonal is 1 by construction.

Panels are laid out as a cohort x gender grid (rows ALL/Po10/WA, cols M/W/B),
matching ``p01_opt_nu.py``. nu runs from inf (Gaussian) down to the grid floor.
One CV solver is used (anderson by default — solvers share the fixed point).
A shared colorbar is auto-ranged to the off-diagonal minimum so the (typically
tiny) differences near 1 are actually visible; cells are annotated with 1 - corr
in units of 1e-3 to read the gap directly.

One figure ->  results/model_selection/baseline/vcorr_vs_nu.png

Self-contained; no arguments needed.

Run::

    python scripts/02_model_selection/baseline/p02_vcorr_vs_nu.py
    python scripts/02_model_selection/baseline/p02_vcorr_vs_nu.py --solver cv
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"

COHORT_ORDER = ["ALL", "Po10", "WA"]   # rows
GENDER_ORDER = ["M", "W", "B"]         # columns
_GENDER_LABEL = {"M": "men", "W": "women", "B": "both"}


def _split_slice(name: str) -> tuple[str, str]:
    cohort, _, gender = name.rpartition("_")
    return cohort, gender


def _nu_label(nu: float) -> str:
    return "inf" if not np.isfinite(nu) else (f"{nu:g}")


def _corr_matrix(V: np.ndarray, method: str) -> np.ndarray:
    """Column-wise correlation matrix of V (races x nus). method in {pearson, spearman}."""
    M = V
    if method == "spearman":
        M = np.apply_along_axis(lambda c: pd.Series(c).rank().to_numpy(), 0, V)
    return np.corrcoef(M, rowvar=False)


def _combined(V: np.ndarray) -> np.ndarray:
    """Upper triangle = Pearson, lower = Spearman, diagonal = 1."""
    pear = _corr_matrix(V, "pearson")
    spear = _corr_matrix(V, "spearman")
    out = np.tril(spear, -1) + np.triu(pear, 1)
    np.fill_diagonal(out, 1.0)
    return out


def _panel_matrix(df_slice: pd.DataFrame, solver: str):
    """Return (combined-corr matrix, nu labels) for one slice + solver."""
    g = df_slice[df_slice.solver == solver]
    if g.empty:                      # fall back to whatever solver is present
        g = df_slice[df_slice.solver == sorted(df_slice.solver.unique())[0]]
    nus = sorted(g["nu"].unique(), key=lambda x: (np.isfinite(x), -x))  # inf first, then desc
    piv = g.pivot(index="race_id", columns="nu", values="v").reindex(columns=nus).dropna()
    V = piv.to_numpy()
    return _combined(V), [_nu_label(n) for n in nus]


def _draw_panel(ax, C: np.ndarray, labels: list[str], *, vmin, cmap):
    n = len(labels)
    im = ax.imshow(C, cmap=cmap, vmin=vmin, vmax=1.0)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=7); ax.set_yticklabels(labels, fontsize=7)
    # annotate off-diagonal cells with 1 - corr in units of 1e-3.
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            gap = (1.0 - C[i, j]) * 1e3
            txt = "0" if gap < 0.05 else f"{gap:.1f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                    color="white" if C[i, j] < (vmin + 1.0) / 2 else "black")
    return im


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver", default="anderson",
                    help="CV solver whose full-data v to use (default anderson).")
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    paths = sorted(OUT_ROOT.glob("*/v_xnu.csv"))
    if not paths:
        print(f"No v_xnu.csv under {OUT_ROOT} — run e01_nu_cv.py first.")
        return

    # load every slice; collect matrices first so the colorbar can be shared.
    panels: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for p in paths:
        name = p.parent.name
        # slug is e.g. ALL_M_14-25_mrc2; cohort/gender are the first two tokens.
        cohort, gender = name.split("_")[0], name.split("_")[1]
        df = pd.read_csv(p)
        panels[(cohort, gender)] = _panel_matrix(df, args.solver)

    cohorts = [c for c in COHORT_ORDER if any(k[0] == c for k in panels)]
    genders = [g for g in GENDER_ORDER if any(k[1] == g for k in panels)]
    nrow, ncol = len(cohorts), len(genders)

    # shared color floor = smallest off-diagonal corr across all panels.
    off_min = min(
        float(np.min(C[~np.eye(C.shape[0], dtype=bool)]))
        for C, _ in panels.values()
    )
    vmin = max(0.0, off_min - 0.001)
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
            im = _draw_panel(ax, C, labels, vmin=vmin, cmap=cmap)
            ax.set_title(f"{cohort} — {_GENDER_LABEL.get(gender, gender)}", fontsize=11)
            if r == nrow - 1:
                ax.set_xlabel("nu")
            if c == 0:
                ax.set_ylabel("nu")

    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.6, location="right")
        cb.set_label("corr(v(nu_i), v(nu_j))")
    fig.suptitle("Baseline v stability across nu  —  upper △ Pearson, lower △ Spearman "
                 "(cells: 1−corr ×1e-3)", fontsize=12)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "vcorr_vs_nu.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
