"""Plot the v_j-bias collapse: year-residualized delta vs v, no-d | +d.

One row per (slice, mrc), two columns -- v from the no-d (agingS4gv) fit and the
+d (full) fit -- scattering the year-partialled race tilt ``delta_j`` against the
year-partialled race factor ``v_j``, with the OLS line and partial-pearson
annotated. The bias signature is a clear positive tilt on the left (no-d) that
collapses toward flat on the right (+d).

Reads one delta estimator's table (``--estimator loo`` default, the credible
one) and recomputes the residuals on the fly (cheap; independent of q03's CSVs).

Output -> results/validation/d_i/vj_bias_<estimator>.png

Run::

    python scripts/04_validation/d_i/p01_vj_bias.py --slices Po10_W
    python scripts/04_validation/d_i/p01_vj_bias.py --estimator eb
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import load_slice, registry                # noqa: E402

from baseline_common import slices as S                          # noqa: E402
import delta_common as DC                                        # noqa: E402  (sibling)

FIT_COL = {"no-d": "#e6194b", "+d": "#4363d8"}


def _scatter(ax, dv, vv, color, title):
    ax.scatter(dv, vv, s=10, alpha=0.5, color=color, edgecolors="none")
    if dv.size >= 3 and np.var(dv) > 1e-18:
        b = DC.ols_slope(vv, dv)
        a = vv.mean() - b * dv.mean()
        xs = np.array([dv.min(), dv.max()])
        ax.plot(xs, a + b * xs, color="k", lw=1.2)
        ax.set_title(f"{title}   r={DC.pearson(dv, vv):+.3f}", fontsize=10)
    else:
        ax.set_title(f"{title}   r=-", fontsize=10)
    ax.axhline(0, color="grey", lw=0.4, alpha=0.5)
    ax.axvline(0, color="grey", lw=0.4, alpha=0.5)
    ax.grid(True, alpha=0.2, lw=0.4)


def _rows(name: str, mrc: int, nu: float, table):
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    adir, fdir = DC.aging_dir(spec, nu), DC.full_dir(spec, nu)
    sub = table[table["slug"] == slug].sort_values("race_idx") if table is not None else None
    if sub is None or not len(sub) or not (DC.present(adir) and DC.present(fdir)):
        return None
    fd = load_slice(spec)
    idx = sub["race_idx"].to_numpy()
    year = DC.race_year(fd)[idx]
    delta = DC.partial_on_year(sub["delta"].to_numpy(), year)
    v_nod = DC.partial_on_year(np.asarray(registry.load_fit(adir, fd).params["v"], float)[idx], year)
    v_d = DC.partial_on_year(np.asarray(registry.load_fit(fdir, fd).params["v"], float)[idx], year)
    return slug, delta, v_nod, v_d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5])
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--estimator", choices=["eb", "loo"], default="loo")
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    table = DC.read_delta(args.estimator)
    if table is None:
        ap.error(f"no delta_{args.estimator}.csv; run the producer first.")

    panels = []
    for name in names:
        for mrc in mrcs:
            r = _rows(name, mrc, args.nu, table)
            if r is not None:
                panels.append(r)
    if not panels:
        print("Nothing to plot (no slice had a delta row + both fits).")
        return

    nrow = len(panels)
    fig, axes = plt.subplots(nrow, 2, figsize=(9.0, 3.4 * nrow),
                             squeeze=False, constrained_layout=True)
    for i, (slug, delta, v_nod, v_d) in enumerate(panels):
        _scatter(axes[i][0], delta, v_nod, FIT_COL["no-d"], f"{slug}  no-d")
        _scatter(axes[i][1], delta, v_d, FIT_COL["+d"], f"{slug}  +d")
        axes[i][0].set_ylabel("v_j  (year-resid.)")
    for ax in axes[-1]:
        ax.set_xlabel(f"delta_{args.estimator}  (year-resid.)")

    fig.suptitle(f"v_j-bias collapse: delta_{args.estimator.upper()} vs v_j  "
                 f"(no-d | +d), nu={args.nu:g}", fontsize=12)
    DC.OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = DC.OUT_ROOT / f"vj_bias_{args.estimator}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}   ({len(panels)} panel row(s))")


if __name__ == "__main__":
    main()
