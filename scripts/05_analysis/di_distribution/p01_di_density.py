"""Career-drift d_i distribution figures (fast; reads precomputed band).

Two figures off the production mrc2 AxD fit (eligible athletes, n_i >= floor):

  di_density.{png,pdf}   grid: rows = persistence floor n_i >= {3,5,10},
                         cols = cohort {ALL, Po10}; men vs women density of d_i,
                         x in %/yr (= d_i * 100). Dashed line at 0 splits
                         improvers (left) from decliners (right). Shaded band =
                         95% athlete-bootstrap envelope, **precomputed by q01**
                         (read from density_band.csv -- this script does NO
                         resampling). Run q01 with --n-band to build/grow the band.
  du_scatter.{png,pdf}   rows = sex, cols = cohort; hexbin of ability u (x) vs
                         drift d_i (y, %/yr) with the OLS line + Spearman corr.

Style is provisional pending a shared figure-style module.

Run::

    python scripts/05_analysis/di_distribution/p01_di_density.py
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

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

import di_common as DI  # noqa: E402

FLOORS = (3, 5, 10)
SEX_COLOR = {"M": "#2c7fb8", "W": "#d7301f"}
SEX_LABEL = {"M": "men", "W": "women"}
BAND_CSV = DI.OUT_ROOT / "density_band.csv"


def fig_density(fits: dict, band: pd.DataFrame, nu: float, mrc: int) -> None:
    nrow, ncol = len(FLOORS), len(DI.POPS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.6 * nrow),
                             sharex=True, sharey="row", constrained_layout=True)
    for i, floor in enumerate(FLOORS):
        for j, cohort in enumerate(DI.POPS):
            ax = axes[i, j]
            ax.axvline(0.0, color="0.5", lw=0.8, ls="--", zorder=1)
            ax.grid(True, color="0.92", lw=0.6, zorder=0)
            for sx in ("M", "W"):
                b = band[(band.cohort == cohort) & (band.sex == sx)
                         & (band.n_floor == floor)].sort_values("x")
                if b.empty:
                    continue
                x = b.x.to_numpy()
                if b.lo.notna().any():
                    ax.fill_between(x, b.lo, b.hi, color=SEX_COLOR[sx],
                                    alpha=0.20, lw=0, zorder=2)
                f = fits.get((cohort, sx))
                lbl = SEX_LABEL[sx]
                if f is not None:
                    t = f.table(min_n=floor)
                    lbl += f" (N={len(t):,}, improvers {np.mean(t.d < 0):.0%})"
                ax.plot(x, b.dens, color=SEX_COLOR[sx], lw=2.0, label=lbl, zorder=3)
            if i == 0:
                ax.set_title(cohort, fontsize=12)
            if j == 0:
                ax.set_ylabel(f"n_i >= {floor}\ndensity", fontsize=10)
            ax.legend(fontsize=7, frameon=False, loc="upper right")
    for ax in axes[-1, :]:
        ax.set_xlabel("career drift  d_i   (%/yr;  <0 = improver)")
    nb = int(band["n_draws"].max()) if not band.empty else 0
    fig.suptitle(f"Career-drift d_i distribution by sex and persistence floor "
                 f"(full nu{nu:g}, mrc{mrc} fit; band = 95% bootstrap, "
                 f"{nb} draws)", fontsize=13)
    for ext in ("png", "pdf"):
        out = DI.OUT_ROOT / f"di_density.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close(fig)


def fig_du_scatter(fits: dict, nu: float, mrc: int, min_n: int = 3) -> None:
    sexes = ["M", "W"]
    fig, axes = plt.subplots(len(sexes), len(DI.POPS),
                             figsize=(5.2 * len(DI.POPS), 4.4),
                             sharex="col", sharey=True, constrained_layout=True)
    for i, sx in enumerate(sexes):
        for j, cohort in enumerate(DI.POPS):
            ax = axes[i, j]
            f = fits.get((cohort, sx))
            if f is None:
                ax.set_axis_off()
                continue
            t = f.table(min_n=min_n)
            u = t.u.to_numpy()
            d = t.d.to_numpy() * 100.0
            ax.hexbin(u, d, gridsize=45, cmap="Blues", mincnt=1, linewidths=0)
            ax.axhline(0.0, color="0.5", lw=0.8, ls="--")
            b, a = np.polyfit(u, d, 1)
            xs = np.array([u.min(), u.max()])
            ax.plot(xs, a + b * xs, color="#d7301f", lw=1.8)
            rho = DI._spearman(t.d.to_numpy(), u)
            ax.text(0.03, 0.97, f"{cohort} {SEX_LABEL[sx]}\nrho(d,u)={rho:+.2f}",
                    transform=ax.transAxes, va="top", fontsize=9)
            if i == len(sexes) - 1:
                ax.set_xlabel("ability  u_i   (log;  <0 = faster)")
            if j == 0:
                ax.set_ylabel("drift d_i (%/yr)")
    fig.suptitle(f"Career drift vs ability (full nu{nu:g}, mrc{mrc} fit, n_i>={min_n})",
                 fontsize=13)
    for ext in ("png", "pdf"):
        out = DI.OUT_ROOT / f"du_scatter.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--mrc", type=int, default=2)
    args = ap.parse_args()

    fits = {}
    for cohort in DI.POPS:
        for sx in ("M", "W"):
            fits[(cohort, sx)] = DI.load_cached(cohort, sx, args.nu, args.mrc)
    if all(v is None for v in fits.values()):
        print("no fits found; nothing plotted.")
        return
    DI.OUT_ROOT.mkdir(parents=True, exist_ok=True)

    if BAND_CSV.is_file():
        fig_density(fits, pd.read_csv(BAND_CSV), args.nu, args.mrc)
    else:
        print(f"[skip density] {BAND_CSV} not found -- run q01 first (it writes the "
              "point density + band there, even with --n-band 0).")
    fig_du_scatter(fits, args.nu, args.mrc)


if __name__ == "__main__":
    main()
