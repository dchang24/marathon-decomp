"""
plot_finish_time_age_hist.py
============================
Descriptive figure for the paper's Dataset section (§2.1): the distribution of
marathon finish times, split by sex and stacked by age group.

Two side-by-side panels (left = men, right = women). Each is a histogram of
finish time in 5-minute bins, with every bar subdivided (stacked) into age-on-
race-day brackets: <=34, 35-44, 45-54, 55-64, >=65.

Data is the headline ALL_B analysis slice (the distilled dataset reported in the
paper), read via `load_slice` so the figure is consistent with the model:

  * finish time   = exp(y)               (y is log finish-time seconds)
  * race-day age  = A_e[row] + A_n        (entry age + years since first race),
                    where entry age uses the athlete's year-of-birth MIDPOINT
                    (0.5*(yob_min+yob_max)); finishes with no yob are dropped.

Output (PNG + PDF) -> data/misc/.

Usage
-----
    python scripts/01_data_prep/plot_finish_time_age_hist.py
    python scripts/01_data_prep/plot_finish_time_age_hist.py --tmin 120 --tmax 420
"""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from marathon_decomp import SliceSpec, load_slice
from marathon_decomp.config import MISC_DIR, PAPER_FIG_DIR

# Age brackets: edges fed to np.digitize, plus display labels (one more label
# than internal edges). Older = later index = darker colour.
AGE_EDGES = [35, 45, 55, 65]
AGE_LABELS = ["<=34", "35-44", "45-54", "55-64", ">=65"]
AGE_COLORS = ["#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tmin", type=float, default=120.0,
                    help="lower finish-time axis bound, minutes (default 120)")
    ap.add_argument("--tmax", type=float, default=420.0,
                    help="upper finish-time axis bound, minutes (default 420)")
    ap.add_argument("--bin", type=float, default=5.0,
                    help="finish-time bin width, minutes (default 5)")
    args = ap.parse_args()

    print("Loading headline slice ALL_B (mrc2) ...")
    fd = load_slice(SliceSpec(cohort="ALL", sex="ALL", min_race_count=2))

    t_min = np.exp(fd.y) / 60.0                    # finish time, minutes
    age = fd.A_e[fd.row_idx] + fd.A_n              # race-day age (yob midpoint)
    sex = np.asarray(fd.athlete_sex)[fd.row_idx]

    known = ~np.isnan(age)
    print(f"  finishes: {fd.N:,}; with known age: {known.sum():,} "
          f"({known.mean() * 100:.1f}%); dropping {(~known).sum():,} unknown-age.")

    edges = np.arange(args.tmin, args.tmax + args.bin, args.bin)
    age_grp = np.digitize(age, AGE_EDGES)          # 0..4

    import pathlib
    style_path = pathlib.Path(__file__).parent.parent / "paper.mplstyle"
    if style_path.exists():
        plt.style.use(str(style_path))
    else:
        # Fallback if run from project root
        plt.style.use("scripts/paper.mplstyle")

    fig, axes = plt.subplots(1, 2, sharey=True)
    for ax, (label, code) in zip(axes, [("(a) Men", "M"), ("(b) Women", "W")]):
        m = known & (sex == code)
        # one count array per age group, stacked
        stacks = [t_min[m & (age_grp == g)] for g in range(len(AGE_LABELS))]
        ax.hist(stacks, bins=edges, stacked=True, color=AGE_COLORS,
                label=AGE_LABELS, edgecolor="white", linewidth=0.2)
        ax.set_title(f"{label}  (N = {m.sum():,})", fontsize=11)
        ax.set_xlabel("Finish time (H:MM)")
        ax.set_xlim(args.tmin, args.tmax)
        # major ticks every 60 min, minor every 30 min
        major_ticks = np.arange(args.tmin, args.tmax + 1, 60)
        minor_ticks = np.arange(args.tmin, args.tmax + 1, 30)
        ax.set_xticks(major_ticks)
        ax.set_xticklabels([f"{int(x) // 60}:00" for x in major_ticks], fontsize=8)
        ax.set_xticks(minor_ticks, minor=True)
        ax.grid(axis="y", alpha=0.25)
        # vertical grid line at every 30-min tick
        for x in minor_ticks:
            ax.axvline(x, color="0.8", lw=0.5, zorder=0)

    axes[0].set_ylabel("Finishes")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[1].legend(handles[::-1], labels[::-1], title="Age on race day",
                   loc="upper right", frameon=True,
                   fontsize=8, title_fontsize=9)
    fig.tight_layout()

    out = MISC_DIR / "fig_finish_time_age_hist.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"[write] {out} (+ .pdf)")

    out_paper = PAPER_FIG_DIR / "fig_finish_time_age_hist.pdf"
    out_paper.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_paper, bbox_inches="tight")
    print(f"[write] {out_paper}")


if __name__ == "__main__":
    main()
