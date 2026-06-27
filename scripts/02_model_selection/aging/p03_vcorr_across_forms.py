"""How stable is the fitted race factor v across aging forms? -- per slice.

Mirror of ``scripts/02_model_selection/baseline/p02_vcorr_vs_nu.py`` but with the
candidate aging *forms* on the axes instead of nu. For each slice, at a fixed
(nu, solver), cell (i, j) is the correlation between v fitted under form_i and
form_j (same race set, mean-v=0 gauge). Split on the diagonal -- upper triangle
Pearson, lower triangle Spearman. Panels laid out as a cohort x gender grid.
Reads each ``{slug}/grid/v_xform.parquet`` -- no refit.

One figure -> results/model_selection/aging/vcorr_across_forms.png

Run::

    python scripts/02_model_selection/aging/p03_vcorr_across_forms.py
    python scripts/02_model_selection/aging/p03_vcorr_across_forms.py --nu inf
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
COHORT_ORDER = ["ALL", "Po10", "WA"]
GENDER_ORDER = ["M", "W", "B"]
_GLABEL = {"M": "men", "W": "women", "B": "both"}


def _short(cand: str) -> str:
    """spline5-gvarying -> S5v ; poly2-gscalar -> P2s ; ...-goff -> ...o."""
    basis, _, gtok = cand.partition("-g")
    head = ("P" if basis.startswith("poly") else "S") + "".join(ch for ch in basis if ch.isdigit())
    return head + (gtok[:1] if gtok else "")


def _form_key(cand: str) -> tuple[int, int, str]:
    basis = cand.split("-g")[0]
    return (0 if basis.startswith("poly") else 1,
            int("".join(ch for ch in basis if ch.isdigit()) or 0), cand)


def _combined(V: np.ndarray) -> np.ndarray:
    pear = np.corrcoef(V, rowvar=False)
    R = np.apply_along_axis(lambda c: pd.Series(c).rank().to_numpy(), 0, V)
    spear = np.corrcoef(R, rowvar=False)
    out = np.tril(spear, -1) + np.triu(pear, 1)
    np.fill_diagonal(out, 1.0)
    return out


def _panel(path: Path, nu: str, solver: str):
    df = pd.read_parquet(path)
    g = df[(df.nu == nu) & (df.solver == solver)]
    if g.empty:
        return None
    cands = sorted(g.cand.unique(), key=_form_key)
    piv = g.pivot(index="race_id", columns="cand", values="v").reindex(columns=cands).dropna()
    if piv.shape[1] < 2:
        return None
    return _combined(piv.to_numpy()), [_short(c) for c in cands]


def _draw(ax, C, labels, *, vmin, cmap):
    n = len(labels)
    im = ax.imshow(C, cmap=cmap, vmin=vmin, vmax=1.0)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=6, rotation=90); ax.set_yticklabels(labels, fontsize=6)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            gap = (1.0 - C[i, j]) * 1e3
            ax.text(j, i, "0" if gap < 0.05 else f"{gap:.1f}", ha="center", va="center",
                    fontsize=5, color="white" if C[i, j] < (vmin + 1.0) / 2 else "black")
    return im


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nu", default="8")
    ap.add_argument("--solver", default="anderson")
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    paths = sorted(OUT_ROOT.glob("*/grid/v_xform.parquet"))
    if not paths:
        print("No v_xform.parquet -- run e01_aging_grid first.")
        return

    panels: dict[tuple[str, str], tuple] = {}
    for p in paths:
        slug = p.parent.parent.name
        cohort, gender = slug.split("_")[0], slug.split("_")[1]
        res = _panel(p, args.nu, args.solver)
        if res is not None:
            panels[(cohort, gender)] = res
    if not panels:
        print(f"No panels for nu={args.nu} solver={args.solver}.")
        return

    cohorts = [c for c in COHORT_ORDER if any(k[0] == c for k in panels)]
    genders = [g for g in GENDER_ORDER if any(k[1] == g for k in panels)]
    off_min = min(float(np.min(C[~np.eye(C.shape[0], dtype=bool)])) for C, _ in panels.values())
    vmin = max(0.0, off_min - 0.001)
    cmap = plt.get_cmap(args.cmap)

    fig, axes = plt.subplots(len(cohorts), len(genders),
                             figsize=(4.2 * len(genders), 4.0 * len(cohorts)),
                             squeeze=False, constrained_layout=True)
    im = None
    for r, cohort in enumerate(cohorts):
        for cc, gender in enumerate(genders):
            ax = axes[r, cc]
            key = (cohort, gender)
            if key not in panels:
                ax.set_axis_off(); continue
            C, labels = panels[key]
            im = _draw(ax, C, labels, vmin=vmin, cmap=cmap)
            ax.set_title(f"{cohort} - {_GLABEL.get(gender, gender)}", fontsize=11)
    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.6, location="right")
        cb.set_label("corr(v(form_i), v(form_j))")
    fig.suptitle(f"v stability across aging forms  (nu={args.nu}, {args.solver})  "
                 "-- upper triangle Pearson, lower Spearman (cells: 1-corr x1e-3)", fontsize=12)

    out = OUT_ROOT / f"vcorr_across_forms_{args.nu}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
