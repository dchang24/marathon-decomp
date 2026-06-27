"""Overlay the fitted aging curves of all bases (poly | spline) for one slice.

Two panels (poly bases | spline bases) on shared axes, one curve per basis at a
fixed (nu, gamma_form, solver); poly2 drawn bold as the canonical reference. The
A_n density rug along the bottom shows where the data is (so tail wiggle in the
sparse high-A_n zone is obvious). Reads e01's ``grid/{curves,an_density}.parquet``
-- no refit.

One figure -> results/model_selection/aging/{slug}/grid/fig_aging_curves_{nu}_{gamma}.png

Run::

    python scripts/02_model_selection/aging/p01_aging_curves.py --slice ALL_M
    python scripts/02_model_selection/aging/p01_aging_curves.py --slice Po10_M --nu inf --gamma varying
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


def _resolve_slug(s: str) -> str | None:
    hits = sorted(p.parent.parent.name for p in OUT_ROOT.glob(f"{s}*/grid/curves.parquet"))
    return hits[0] if hits else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default=None, help="slice name or slug (default: first available).")
    ap.add_argument("--nu", default="8")
    ap.add_argument("--gamma", default="scalar", choices=["off", "scalar", "varying"])
    ap.add_argument("--solver", default="anderson")
    args = ap.parse_args()

    if args.slice is None:
        cands = sorted(p.parent.parent.name for p in OUT_ROOT.glob("*/grid/curves.parquet"))
        slug = cands[0] if cands else None
    else:
        slug = _resolve_slug(args.slice)
    if slug is None:
        print("No e01 curves.parquet found -- run e01_aging_grid first.")
        return

    grid_dir = OUT_ROOT / slug / "grid"
    c = pd.read_parquet(grid_dir / "curves.parquet")
    c = c[(c.nu == args.nu) & (c.gamma_form == args.gamma) & (c.solver == args.solver)]
    if c.empty:
        print(f"No curves for nu={args.nu} gamma={args.gamma} solver={args.solver} in {slug}")
        return
    dens = pd.read_parquet(grid_dir / "an_density.parquet")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, constrained_layout=True)
    for ax, kind, title in ((axes[0], "poly", "polynomial"), (axes[1], "spline", "natural cubic spline")):
        sub = c[c.basis.str.startswith(kind)]
        for basis, gb in sorted(sub.groupby("basis"), key=lambda kv: kv[0]):
            gb = gb.sort_values("A_n")
            ref = basis == "poly2"
            ax.plot(gb.A_n, gb.aging_curve, label=basis, lw=2.6 if ref else 1.4,
                    color="black" if ref else None, zorder=5 if ref else 2)
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
        # A_n density rug along the bottom
        ymin = c.aging_curve.min()
        h = 0.06 * (c.aging_curve.max() - ymin + 1e-9)
        w = dens.bin_hi - dens.bin_lo
        ax.bar(dens.bin_lo, h * dens["count"] / dens["count"].max(), width=w, bottom=ymin - h * 1.2,
               align="edge", color="0.8", zorder=1)
        ax.set_title(f"{title} bases"); ax.set_xlabel("A_n  (years since debut)")
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel("aging curve  theta . B(A_n)  (log-time)")
    fig.suptitle(f"{slug}  aging curves  (nu={args.nu}, gamma={args.gamma}, {args.solver})  "
                 f"-- poly2 bold; grey = A_n density", fontsize=12)

    out = grid_dir / f"fig_aging_curves_{args.nu}_{args.gamma}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
