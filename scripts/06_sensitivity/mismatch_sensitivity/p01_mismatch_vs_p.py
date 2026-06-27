"""p01 -- the headline mismatch-sensitivity figures.

Main figure ``fig_mismatch_vs_p.png`` -- the same three rows plotted against TWO
x-axes, one per column, so the intended knob and the realised (gated) rate sit
side by side:
  (left col)  x = requested per-athlete rate ``p``        -- the intended knob;
  (right col) x = realised per-athlete rate ``real_ath_mean`` -- AFTER the join
              compatibility gating (same sex / |dlog-time|<=tau / same-day /
              yob), so it is below ``p`` wherever a merge found no partner.
Rows:
  (top)    Spearman of v-ranking vs the cold baseline (median + 2.5-97.5% band);
  (middle) median per-race |delta v| in minutes-at-3:00 (a GLOBAL summary);
  (bottom) per-race sensitivity ratio = mismatch SD / bootstrap SD on a log
           y-axis -- median (solid), 75th pct (dashed), and max (dotted) across
           races -- with a y=1 line (above 1 -> identity error exceeds ordinary
           sampling noise) and the p* crossing annotated (left column only).
           Omitted if q01 had no bootstrap reference.

Companion ``fig_susceptible_races.png`` -- the INDIVIDUAL view: at the largest p
per op, the top-N races by mismatch SD of v (the most susceptible races), with
each race's bootstrap SD overlaid as a diamond. A bar extending past its diamond
is a race the identity error moves more than ordinary sampling noise.

Reads q01's summary.parquet (+ meta.json for p*) and per_race.parquet. One set
of figures per group.

Run::
    python scripts/06_sensitivity/mismatch_sensitivity/p01_mismatch_vs_p.py --group results/sensitivity/mismatch_sensitivity/Po10_B_14-25_mrc2/full_nu8p00_best__7cde3824

    python scripts/06_sensitivity/mismatch_sensitivity/p01_mismatch_vs_p.py --group results/sensitivity/mismatch_sensitivity/ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

_COLORS = {"break": "#1f77b4", "join": "#d62728", "both": "#2ca02c"}
MIN_PER_3H = 180.0       # delta minutes ~= MIN_PER_3H * delta(log time)
TOPN = 12                # races shown in the susceptible-races figure


def _main_figure(summary: pd.DataFrame, pstar: dict, sdir: Path) -> Path:
    """Three rows x two columns: each metric vs requested p (left) and vs the
    realised per-athlete rate (right)."""
    has_sens = summary["sens_ratio_median"].notna().any()
    nrows = 3 if has_sens else 2
    # the two x-axes: (column key, axis label, whether p* vlines apply)
    xcols = [("p", "requested per-athlete rate p  (intended)", True),
             ("real_ath_mean", "realised per-athlete rate  (after gating)", False)]

    fig, axes = plt.subplots(nrows, 2, figsize=(12, 3 * nrows),
                             sharey="row", squeeze=False)

    for j, (xcol, xlab, is_p) in enumerate(xcols):
        for op, g in summary.groupby("op"):
            g = g.sort_values("p")
            x = g[xcol].to_numpy()
            c = _COLORS.get(op)
            axes[0][j].plot(x, g.spearman_median, "-o", color=c, label=op)
            axes[0][j].fill_between(x, g.spearman_lo, g.spearman_hi, color=c, alpha=0.15)
            axes[1][j].plot(x, g.median_abs_dv_min, "-o", color=c, label=op)
            if has_sens:
                axes[2][j].plot(x, g.sens_ratio_median, "-o", color=c, label=op)
                axes[2][j].plot(x, g.sens_ratio_p75, "--", color=c, lw=1)
                axes[2][j].plot(x, g.sens_ratio_max, ":", color=c, lw=1)
        axes[-1][j].set_xlabel(xlab)
        if has_sens:
            axes[2][j].axhline(1.0, ls="-", color="k", lw=0.8, alpha=0.6)
            axes[2][j].set_yscale("log")
            if is_p:
                for op, pv in pstar.items():
                    if pv is not None and np.isfinite(pv):
                        axes[2][j].axvline(pv, ls=":", color=_COLORS.get(op), lw=1)

    # row labels + legends on the left column only
    axes[0][0].set_ylabel("Spearman vs point fit")
    axes[0][0].legend(title="operation", fontsize=8)
    axes[1][0].set_ylabel("median |dv|  (min @ 3:00)")
    if has_sens:
        axes[2][0].set_ylabel("mismatch SD / bootstrap SD\n(per-race; med / p75 / max)")
        handles = [Line2D([0], [0], color="0.3", ls=ls, label=lab)
                   for ls, lab in [("-", "median"), ("--", "p75"), (":", "max")]]
        axes[2][0].legend(handles=handles, fontsize=7, loc="best")
    axes[0][0].set_title("vs intended rate", fontsize=9)
    axes[0][1].set_title("vs realised rate", fontsize=9)

    fig.suptitle(f"Identity-mismatch sensitivity -- {sdir.parent.parent.name}/"
                 f"{sdir.parent.name}", fontsize=9)
    fig.tight_layout()
    out = sdir / "fig_mismatch_vs_p.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def _susceptible_figure(per_race: pd.DataFrame, sdir: Path) -> Path | None:
    """Individual-race view: at the largest p per op, the top-N races by mismatch
    SD of v, with bootstrap SD overlaid. Naming the susceptible races."""
    if per_race.empty:
        return None
    ops = list(per_race["op"].unique())
    fig, axes = plt.subplots(1, len(ops), squeeze=False,
                             figsize=(4.8 * len(ops), 0.42 * TOPN + 1.8))
    has_boot = "boot_sd" in per_race.columns and per_race["boot_sd"].notna().any()
    for ax, op in zip(axes[0], ops):
        sub = per_race[per_race["op"] == op]
        pmax = sub["p"].max()
        sub = sub[sub["p"] == pmax].copy()
        sub["dv_sd_min"] = sub["dv_sd"] * MIN_PER_3H
        # most susceptible at the TOP of the axis -> sort asc, take tail
        sub = sub.sort_values("dv_sd").tail(TOPN)
        y = np.arange(len(sub))
        ax.barh(y, sub["dv_sd_min"], color=_COLORS.get(op), alpha=0.8,
                label="mismatch SD")
        if has_boot and sub["boot_sd"].notna().any():
            ax.plot(sub["boot_sd"] * MIN_PER_3H, y, "D", color="k", ms=5,
                    label="bootstrap SD")
        ax.set_yticks(y)
        ax.set_yticklabels([f"race {int(r)}" for r in sub["race_id"]], fontsize=7)
        ax.set_title(f"{op}  (p={pmax:g})", fontsize=9)
        ax.set_xlabel("per-race SD of v  (min @ 3:00)")
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Most susceptible races -- {sdir.parent.parent.name}/"
                 f"{sdir.parent.name}", fontsize=9)
    fig.tight_layout()
    out = sdir / "fig_susceptible_races.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()

    sdir = Path(args.group).resolve() / "summary"
    summary = pd.read_parquet(sdir / "summary.parquet")
    meta = json.loads((sdir / "meta.json").read_text()) if (sdir / "meta.json").is_file() else {}
    pstar = meta.get("p_star", {})

    out1 = _main_figure(summary, pstar, sdir)
    print(f"wrote {out1}")

    pr_path = sdir / "per_race.parquet"
    if pr_path.is_file():
        out2 = _susceptible_figure(pd.read_parquet(pr_path), sdir)
        if out2 is not None:
            print(f"wrote {out2}")


if __name__ == "__main__":
    main()
