"""Held-out CV log-density per cell across the candidate aging forms.

For one slice: CV log-density / held-out cell (the selection criterion) vs the
candidate form, one line per gamma_form, with BIC on a twin axis; nu as rows
(inf on top, 8 below). The CV-best form per row is starred. Reads e02's
``cv/form_selection.csv`` -- no refit.

One figure -> results/model_selection/aging/{slug}/cv/fig_cv_form.png

Run::

    python scripts/02_model_selection/aging/p02_cv_form.py --slice ALL_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"


def _basis_key(name: str) -> tuple[int, int]:
    # order poly2..6 then spline3..6 on the x-axis
    kind = 0 if name.startswith("poly") else 1
    return (kind, int("".join(ch for ch in name if ch.isdigit()) or 0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default=None, help="slice name or slug (default: first available).")
    args = ap.parse_args()

    glob = (f"{args.slice}*/cv/form_selection.csv" if args.slice else "*/cv/form_selection.csv")
    hits = sorted(OUT_ROOT.glob(glob))
    if not hits:
        print("No cv/form_selection.csv found -- run e02_aging_cv first.")
        return
    path = hits[0]
    slug = path.parent.parent.name
    df = pd.read_csv(path)

    nus = sorted(df.nu.unique(), key=lambda x: (x != "inf", x))   # inf first
    bases = sorted(df.basis.unique(), key=_basis_key)
    fig, axes = plt.subplots(len(nus), 1, figsize=(10, 4.2 * len(nus)),
                             squeeze=False, constrained_layout=True)
    for r, nu in enumerate(nus):
        ax = axes[r, 0]
        sub = df[df.nu == nu]
        ax2 = ax.twinx()
        for gf, gb in sub.groupby("gamma_form"):
            gb = gb.set_index("basis").reindex(bases).reset_index()
            ax.plot(gb.basis, gb.cv_per_cell, marker="o", label=f"gamma={gf}")
            ax2.plot(gb.basis, gb.full_bic, marker="x", ls=":", alpha=0.5)
        best = sub.loc[sub.cv_per_cell.idxmax()]
        ax.scatter([best.basis], [best.cv_per_cell], marker="*", s=320,
                   color="red", zorder=6, label=f"best: {best['cand']}")
        ax.set_title(f"nu = {nu}")
        ax.set_ylabel("CV log-density / cell")
        ax2.set_ylabel("BIC (dotted)", color="0.4")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.3)
    axes[-1, 0].set_xlabel("basis")
    fig.suptitle(f"{slug}  aging-form CV selection  (line: CV/cell  |  dotted: BIC)", fontsize=12)

    out = path.parent / "fig_cv_form.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
