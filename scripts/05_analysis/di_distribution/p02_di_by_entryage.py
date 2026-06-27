"""Career-drift d_i distribution BY ENTRY AGE (gauge-safe individual-progression check).

The varying-gamma block already lets the *global* progression curve depend on
entry age (the entry-age fan, paper 6.2). d_i is a different object: the
*individual* residual career slope, on top of that fan. This script asks whether
the d_i DISTRIBUTION itself shifts with entry age -- i.e. does any entry-age
structure survive into the individual deviations, or did gamma absorb it?

Gauge note (decisive for what is plotted). The level/mean of d_i is NOT identified
(a constant c added to every d_i is absorbed by the aging linear coef + u_i). So
absolute d_i (and "improver vs decliner") is gauge-pinned and NOT plotted as such.
But the gauge shifts EVERY athlete's d_i by the same c, so DIFFERENCES BETWEEN
entry-age bins -- and the slope of d_i on entry age -- are gauge-INVARIANT.
Everything here is therefore shown RELATIVE TO THE SLICE MEAN d_i: only bin-to-bin
shifts are interpretable, and 0 = the slice-average trajectory, not "no drift".

Expectation (supporting the paper claim): the bins overlap and slope(d ~ A_e) ~ 0
-> gamma captured the entry-age structure; d_i is idiosyncratic heterogeneity, not
leftover entry-age confound. A non-flat gradient would instead flag gamma underfit.

Two figures off the production mrc2 AxD fit (eligible athletes, n_i >= floor):

  di_by_entryage.{png,pdf}    grid rows = cohort {ALL, Po10}, cols = sex {M, W};
                              each panel overlays the d_i density for entry-age
                              quantile bins (centered on the slice mean), %/yr.
  di_entryage_gradient.{png,pdf}  same grid; bin-mean d_i (relative to slice mean)
                              vs bin-median entry age, 95% athlete-bootstrap CI,
                              with the continuous OLS slope(d ~ A_e) line + CI.

Also writes di_entryage_bins.csv (per bin) and di_entryage_slope.csv (per cell:
slope, bootstrap CI, Spearman) so the numbers are quotable.

Run::

    python scripts/05_analysis/di_distribution/p02_di_by_entryage.py
    python scripts/05_analysis/di_distribution/p02_di_by_entryage.py --floor 3 --nbins 5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402
from scipy import stats           # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

import di_common as DI            # noqa: E402

SEXES = ("M", "W")
SEX_LABEL = {"M": "men", "W": "women"}
CI = (2.5, 97.5)
AE_LO, AE_HI = 10.0, 90.0    # plausible entry-age band; drops missing-DOB sentinels


def _entryage_bins(A_e: np.ndarray, nbins: int) -> tuple[np.ndarray, np.ndarray]:
    """Quantile-bin entry age. Returns (bin_idx[N], edges[nbins+1]) with unique edges."""
    qs = np.linspace(0, 1, nbins + 1)
    edges = np.unique(np.quantile(A_e, qs))
    # clip so the rightmost point lands in the last bin
    idx = np.clip(np.digitize(A_e, edges[1:-1], right=False), 0, len(edges) - 2)
    return idx, edges


def _boot(d: np.ndarray, A_e: np.ndarray, bin_idx: np.ndarray, nb: int, *,
          n_boot: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Athlete bootstrap. Returns (slope_draws[n_boot], binmean_rel_draws[nb,n_boot]).

    Each draw: resample athletes, OLS slope of d on A_e, and each bin's mean d
    RELATIVE to that draw's overall mean (gauge-invariant centering).
    """
    N = len(d)
    rng = np.random.default_rng(seed)
    slopes = np.empty(n_boot)
    bmean = np.full((nb, n_boot), np.nan)
    for b in range(n_boot):
        s = rng.integers(0, N, N)
        ds, as_, bs = d[s], A_e[s], bin_idx[s]
        slopes[b] = np.polyfit(as_, ds, 1)[0]
        mu = ds.mean()
        for k in range(nb):
            m = bs == k
            if m.any():
                bmean[k, b] = ds[m].mean() - mu
    return slopes, bmean


def _cell_stats(f, floor: int, nbins: int, n_boot: int, seed: int) -> dict | None:
    """All per-cell quantities for one (cohort, sex) DriftFit."""
    t = f.table(min_n=floor)
    if len(t) < max(50, 5 * nbins):
        return None
    d = t.d.to_numpy()
    A_e = t.A_e.to_numpy()
    ok = np.isfinite(d) & np.isfinite(A_e) & (A_e >= AE_LO) & (A_e <= AE_HI)
    n_drop = int((~ok).sum())
    d, A_e = d[ok], A_e[ok]
    bin_idx, edges = _entryage_bins(A_e, nbins)
    nb = len(edges) - 1

    mu = d.mean()
    d_rel = d - mu                       # gauge-safe centering (slice mean -> 0)
    slope = float(np.polyfit(A_e, d, 1)[0])
    rho = DI._spearman(d, A_e)
    sl_draws, bm_draws = _boot(d, A_e, bin_idx, nb, n_boot=n_boot, seed=seed)
    sl_lo, sl_hi = np.percentile(sl_draws, CI)

    bins = []
    for k in range(nb):
        m = bin_idx == k
        lo, hi = np.percentile(bm_draws[k], CI)
        bins.append(dict(
            bin=k, age_lo=float(edges[k]), age_hi=float(edges[k + 1]),
            age_med=float(np.median(A_e[m])), n=int(m.sum()),
            mean_rel=float(d_rel[m].mean()), ci_lo=float(lo), ci_hi=float(hi),
        ))
    return dict(d_rel=d_rel, A_e=A_e, bin_idx=bin_idx, edges=edges, mu=mu,
                slope=slope, slope_lo=float(sl_lo), slope_hi=float(sl_hi),
                rho=rho, bins=bins, n=len(d), n_drop=n_drop)


# --------------------------------------------------------------------------- #
# figures                                                                     #
# --------------------------------------------------------------------------- #
def fig_density(cells: dict, nu: float, mrc: int, floor: int) -> None:
    fig, axes = plt.subplots(len(DI.POPS), len(SEXES),
                             figsize=(5.4 * len(SEXES), 3.0 * len(DI.POPS)),
                             sharex=True, sharey="row", constrained_layout=True)
    for i, cohort in enumerate(DI.POPS):
        for j, sx in enumerate(SEXES):
            ax = axes[i, j]
            ax.grid(True, color="0.92", lw=0.6, zorder=0)
            ax.axvline(0.0, color="0.5", lw=0.8, ls="--", zorder=1)
            c = cells.get((cohort, sx))
            if c is None:
                ax.set_axis_off()
                continue
            nb = len(c["edges"]) - 1
            cmap = plt.get_cmap("viridis")
            xr = np.percentile(c["d_rel"], [0.5, 99.5]) * 100.0
            grid = np.linspace(xr[0], xr[1], 256)
            for k in range(nb):
                m = c["bin_idx"] == k
                xk = c["d_rel"][m] * 100.0
                if xk.size < 20 or xk.std() < 1e-9:
                    continue
                dens = stats.gaussian_kde(xk)(grid)
                col = cmap(k / max(nb - 1, 1))
                amed = np.median(c["A_e"][m])
                ax.plot(grid, dens, color=col, lw=2.0,
                        label=f"entry ~{amed:.0f} (n={m.sum():,})", zorder=3)
            ax.text(0.03, 0.97, f"{cohort} {SEX_LABEL[sx]}", transform=ax.transAxes,
                    va="top", fontsize=10, fontweight="bold")
            ax.legend(fontsize=7, frameon=False, loc="upper right", title="entry age")
            if j == 0:
                ax.set_ylabel("density")
    for ax in axes[-1, :]:
        ax.set_xlabel("d_i - slice mean   (%/yr;  bin-to-bin shift is gauge-safe, "
                      "0 = slice-avg trajectory)")
    fig.suptitle(f"Individual career-drift d_i by entry age "
                 f"(full nu{nu:g}, mrc{mrc}, n_i>={floor}; centered on slice mean)",
                 fontsize=12)
    for ext in ("png", "pdf"):
        out = DI.OUT_ROOT / f"di_by_entryage.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close(fig)


def fig_gradient(cells: dict, nu: float, mrc: int, floor: int) -> None:
    fig, axes = plt.subplots(len(DI.POPS), len(SEXES),
                             figsize=(5.4 * len(SEXES), 3.0 * len(DI.POPS)),
                             sharex=True, sharey=True, constrained_layout=True)
    for i, cohort in enumerate(DI.POPS):
        for j, sx in enumerate(SEXES):
            ax = axes[i, j]
            ax.grid(True, color="0.92", lw=0.6, zorder=0)
            ax.axhline(0.0, color="0.5", lw=0.8, ls="--", zorder=1)
            c = cells.get((cohort, sx))
            if c is None:
                ax.set_axis_off()
                continue
            xs = np.array([b["age_med"] for b in c["bins"]])
            ys = np.array([b["mean_rel"] for b in c["bins"]]) * 100.0
            lo = np.array([b["ci_lo"] for b in c["bins"]]) * 100.0
            hi = np.array([b["ci_hi"] for b in c["bins"]]) * 100.0
            ax.errorbar(xs, ys, yerr=[ys - lo, hi - ys], fmt="o-", color="#2c7fb8",
                        lw=1.6, ms=5, capsize=3, zorder=3)
            # continuous OLS line through (mean A_e, 0), slope in %/yr per year
            ag = np.linspace(c["A_e"].min(), c["A_e"].max(), 50)
            ax.plot(ag, c["slope"] * 100.0 * (ag - c["A_e"].mean()),
                    color="#d7301f", lw=1.4, ls="--", zorder=2)
            sl10 = c["slope"] * 100.0 * 10.0          # %/yr per decade of entry age
            slo10, shi10 = c["slope_lo"] * 1000.0, c["slope_hi"] * 1000.0
            flat = "flat (gamma absorbed entry age)" if slo10 <= 0 <= shi10 else "gradient"
            ax.text(0.03, 0.97,
                    f"{cohort} {SEX_LABEL[sx]}\n"
                    f"slope {sl10:+.2f} [{slo10:+.2f},{shi10:+.2f}] %/yr per 10yr\n"
                    f"rho(d,A_e)={c['rho']:+.2f}  -> {flat}",
                    transform=ax.transAxes, va="top", fontsize=8)
            if j == 0:
                ax.set_ylabel("mean d_i - slice mean (%/yr)")
    for ax in axes[-1, :]:
        ax.set_xlabel("entry age (yr at debut)")
    fig.suptitle(f"Entry-age gradient of individual drift d_i "
                 f"(full nu{nu:g}, mrc{mrc}, n_i>={floor}; bin shifts & slope are "
                 f"gauge-safe; compare slope to the ~1.3%/yr individual SD)",
                 fontsize=12)
    for ext in ("png", "pdf"):
        out = DI.OUT_ROOT / f"di_entryage_gradient.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close(fig)


def write_csv(cells: dict, nu: float, mrc: int, floor: int) -> None:
    brows, srows = [], []
    for (cohort, sx), c in cells.items():
        if c is None:
            continue
        for b in c["bins"]:
            brows.append(dict(cohort=cohort, sex=sx, nu=nu, mrc=mrc, floor=floor, **b))
        srows.append(dict(cohort=cohort, sex=sx, nu=nu, mrc=mrc, floor=floor, n=c["n"],
                          slope_per_yr=c["slope"], slope_lo=c["slope_lo"],
                          slope_hi=c["slope_hi"], spearman_d_Ae=c["rho"]))
    if brows:
        pd.DataFrame(brows).to_csv(DI.OUT_ROOT / "di_entryage_bins.csv", index=False)
        pd.DataFrame(srows).to_csv(DI.OUT_ROOT / "di_entryage_slope.csv", index=False)
        print(f"wrote {DI.OUT_ROOT / 'di_entryage_bins.csv'}")
        print(f"wrote {DI.OUT_ROOT / 'di_entryage_slope.csv'}")
        print("\nslope(d_i ~ entry age)  [%/yr per 10yr entry age]:")
        for r in srows:
            print(f"  {r['cohort']:>4} {r['sex']}: "
                  f"{r['slope_per_yr'] * 1000:+.2f} "
                  f"[{r['slope_lo'] * 1000:+.2f}, {r['slope_hi'] * 1000:+.2f}]  "
                  f"rho={r['spearman_d_Ae']:+.3f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--mrc", type=int, default=2)
    ap.add_argument("--floor", type=int, default=5, help="n_i floor (default 5)")
    ap.add_argument("--nbins", type=int, default=4, help="entry-age quantile bins")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    cells: dict = {}
    for cohort in DI.POPS:
        for sx in SEXES:
            f = DI.load_cached(cohort, sx, args.nu, args.mrc)
            cells[(cohort, sx)] = None if f is None else _cell_stats(
                f, args.floor, args.nbins, args.n_boot, args.seed)
    if all(v is None for v in cells.values()):
        print("no fits found; nothing plotted.")
        return
    DI.OUT_ROOT.mkdir(parents=True, exist_ok=True)
    fig_density(cells, args.nu, args.mrc, args.floor)
    fig_gradient(cells, args.nu, args.mrc, args.floor)
    write_csv(cells, args.nu, args.mrc, args.floor)


if __name__ == "__main__":
    main()
