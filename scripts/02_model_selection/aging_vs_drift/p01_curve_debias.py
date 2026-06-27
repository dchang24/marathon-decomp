"""Overlay the population aging curve at {no-d, +d} x {mrc2, mrc5}, one panel/slice.

The visual companion to section 6 of ``q01_grid_compare.py``: the deferred
aging-curve de-biasing check. For each slice it draws the four population aging
curves -- no-d (agingS4gv) and +d (full), at the everyone (mrc2) and
dedicated-runner (mrc5) field cuts -- centered to mean 0 on the shared A_n range
(one APC beta=0 gauge, so only curvature/peak is compared).

Hypothesis (read off the plot): the **no-d / everyone (mrc2)** curve is the
contaminated outlier -- casual runners' improvement trajectories steepen the
everyone-curve. Adding d_i (no-d -> +d at mrc2) OR restricting to dedicated
runners (mrc2 -> mrc5) both pull it onto the other three, which cluster.

Colour = field cut (mrc2 / mrc5); style = drift (no-d dashed / +d solid); the
contaminated no-d/mrc2 curve is drawn bold. Argument-free (VS Code play): every
slice with curves at >=2 mrc.

Output -> results/model_selection/aging_vs_drift/curve_debias.png

Run::

    python scripts/02_model_selection/aging_vs_drift/p01_curve_debias.py
    python scripts/02_model_selection/aging_vs_drift/p01_curve_debias.py --slices Po10_M ALL_W
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from baseline_common import slices as S                          # noqa: E402
# sibling module in this dir (importable despite the digit-prefixed parent).
from q01_grid_compare import OUT_ROOT, _centered_curves_on_common_grid  # noqa: E402

MRC_COLOR = {2: "#e6194b", 5: "#4363d8"}      # everyone / dedicated
VARIANT_STYLE = {"noD": (2.0, "--"), "withD": (2.0, "-")}
LABEL = {("noD", 2): "no-d, mrc2 (everyone)", ("withD", 2): "+d, mrc2",
         ("noD", 5): "no-d, mrc5 (dedicated)", ("withD", 5): "+d, mrc5"}


def _panel(ax, name: str, mrcs: list[int], nu: float) -> bool:
    grid, curves = _centered_curves_on_common_grid(name, mrcs, nu)
    if curves is None:
        ax.set_visible(False)
        return False
    for (m, variant), y in sorted(curves.items()):
        lw, ls = VARIANT_STYLE[variant]
        bold = (variant == "noD" and m == min(c for c, _ in curves))  # contaminated curve
        ax.plot(grid, y, color=MRC_COLOR.get(m, "#808000"), ls=ls,
                lw=lw + (1.4 if bold else 0.0), label=LABEL.get((variant, m), f"{variant} mrc{m}"),
                zorder=3 if bold else 2)
        ax.plot(grid[int(np.argmin(y))], y.min(), "o", ms=4,
                color=MRC_COLOR.get(m, "#808000"), zorder=4)
    ax.axhline(0.0, color="grey", lw=0.5, alpha=0.5)
    ax.set_title(name, fontsize=11, fontweight="bold")
    ax.set_xlabel("A_n  (yr since first race)")
    ax.set_ylabel("centered log-time")
    ax.grid(True, alpha=0.25, lw=0.4)
    ax.legend(fontsize=7, frameon=True)
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5])
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    if len(mrcs) < 2:
        ap.error("need >=2 --mrc values to overlay (e.g. --mrc 2 5).")

    ncol = 2 if len(names) > 1 else 1
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.0 * nrow),
                             squeeze=False, constrained_layout=True)
    flat = axes.ravel()
    n_ok = sum(_panel(flat[i], name, mrcs, args.nu) for i, name in enumerate(names))
    for j in range(len(names), len(flat)):
        flat[j].set_visible(False)

    if n_ok == 0:
        print("No slice had curves at >=2 mrc; nothing plotted.")
        plt.close(fig)
        return
    fig.suptitle(f"Aging-curve de-biasing: no-d vs +d, mrc {mrcs}  (nu={args.nu:g})",
                 fontsize=13)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "curve_debias.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}   ({n_ok}/{len(names)} panels)")


if __name__ == "__main__":
    main()
