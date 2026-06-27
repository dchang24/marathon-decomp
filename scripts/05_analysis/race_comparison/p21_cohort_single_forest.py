"""p21 - single-panel (ALL only) series forest from q20's processed data.

A trimmed sibling of `p20_cohort_connected_forest.py`: draws ONLY the full
dataset (ALL_B) forest, with no Po10 panel and no cross-cohort connectors. Race
labels are moved back to the LEFT of the panel (conventional forest-plot layout),
so the figure reads left-to-right as label -> boxplot. Rows ordered by the
BOOTSTRAP MEDIAN (q20 rank); colour = within-slice M-vs-W contrast (red harder
for men, blue harder for women, grey = none). Reads + draws only -- compute is in
q20, no fitting/gauge math here.

Inputs (run q20 first): results/analysis/race_comparison/cohort_forest/
    cohort_forest_ALL_B.csv, cohort_forest_ALL_B_boot.parquet
Output: results/analysis/race_comparison/fig_cohort_single_forest.{png,pdf}
    (paper copy save is left COMMENTED OUT -- review the figure first).

Run::

    python scripts/05_analysis/race_comparison/p21_cohort_single_forest.py
    python scripts/05_analysis/race_comparison/p21_cohort_single_forest.py --style strip
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from race_common import OUT_ROOT, REF_MARATHON_MIN  # noqa: E402
from marathon_decomp.config import PAPER_FIG_DIR  # noqa: E402

CF_DIR = OUT_ROOT / "cohort_forest"
NEUTRAL = "0.75"

# ===========================================================================
# ADJUSTABLE PARAMETERS
#
# The base look is the shared scripts/paper.mplstyle, loaded in
# apply_paper_style() exactly as the other p* scripts do. The values below are
# this figure's DEPARTURES from that style: a dense 46-row forest with long race
# labels needs a tall canvas and type that runs smaller than the paper defaults.
# Each departure cites the paper.mplstyle default it overrides; "custom" marks
# figure-specific text that has no rcParam / style default at all.
# ===========================================================================

# --- Canvas (departs from paper.mplstyle figure.figsize: 6.1 x 4.0) ----------
FIGSIZE = (6.1, 8)   # tall single panel to fit 46 rows; default is 6.1 x 4.0
LBL_LEFT = 0.0          # custom: race labels sit left of this axes-x-fraction,
                        #   right-aligned hugging the left spine (0.0 = spine).
                        #   bbox_inches="tight" reclaims the left margin.

# --- Font sizes --------------------------------------------------------------
# (1) rcParam-backed text: FS_AXIS is pushed into rcParams in apply_paper_style()
#     (see the DEPARTURES block there), so axis + tick + colorbar text follow it.
FS_AXIS = 8.0           # axis/tick/colorbar text | paper.mplstyle: font.size 9, labelsize 9, ticks 8
# (2) figure-specific text: NO paper.mplstyle equivalent, set per-artist, and
#     deliberately well below the default body size (there are 46 of each).
FS_RACE = 6.0           # custom: left-side race-name labels
FS_COUNT = 5.0          # custom: (men/women) counts hugging each box

# --- Markers / spacing (geometry, not style) ---------------------------------
BOX_WIDTH = 0.62        # boxplot box height (data units = 1 row)
LBL_PAD = 7             # gap (points) between a race label and the spine
COUNT_PAD = 3           # gap (points) between a count and the whisker cap

# --- Colorbar (vertical inset, top-right empty corner) -----------------------
# Fast races sit at the top with their boxes on the LEFT, so the upper-right of
# the panel is blank -- drop a tall vertical colorbar in there instead of
# stealing a strip of figure height along the bottom.
CBAR_RECT = (0.80, 0.58, 0.045, 0.36)  # inset [x0, y0, w, h] in axes fraction
CBAR_NTICKS = 5                         # ticks across [-vmax, +vmax]
# ===========================================================================


def apply_paper_style() -> None:
    """Load the shared paper.mplstyle (as the other p* scripts do), then apply
    THIS figure's style departures.

    Departure: the figure is a single tall panel with 46 rows and long left
    labels, so the default 9 pt body type is too large and collides. We shrink
    the rcParam-backed text to FS_AXIS; the figure-specific labels (FS_RACE,
    FS_COUNT) are set per-artist. Everything else inherits paper.mplstyle.
    """
    style_path = Path(__file__).resolve().parents[2] / "paper.mplstyle"
    plt.style.use(str(style_path) if style_path.exists() else "scripts/paper.mplstyle")
    plt.rcParams.update({                 # DEPARTURES (paper.mplstyle default in comment)
        "font.size":       FS_AXIS,       # 9
        "axes.labelsize":  FS_AXIS,       # 9
        "xtick.labelsize": FS_AXIS,       # 8 (unchanged; coupled to one knob)
        "ytick.labelsize": FS_AXIS,       # 8 (unchanged; coupled to one knob)
    })


def _short(s: str) -> str:
    return s.replace("_marathon", "").replace("_", " ")


def load():
    a = pd.read_csv(CF_DIR / "cohort_forest_ALL_B.csv")
    bootA = pd.read_parquet(CF_DIR / "cohort_forest_ALL_B_boot.parquet")
    return a, bootA


def draw_marker(ax, y, data, color, style, rng):
    if style == "boxplot":
        bp = ax.boxplot([data], positions=[y], vert=False, widths=BOX_WIDTH,
                        patch_artist=True, showfliers=False, manage_ticks=False)
        for box in bp["boxes"]:
            box.set(facecolor=color, edgecolor="black", linewidth=0.5, alpha=0.85)
        for med in bp["medians"]:
            med.set(color="black", linewidth=1.0)
        for wk in bp["whiskers"] + bp["caps"]:
            wk.set(color="0.45", linewidth=0.8)
    elif style == "strip":
        ax.scatter(data, y + (rng.random(len(data)) - 0.5) * 0.55,
                   c=[color], s=5, alpha=0.4, edgecolors="none")
        ax.scatter([np.median(data)], [y], c=[color], s=32,
                   edgecolors="black", linewidths=0.6, zorder=4)
    else:  # ci_box
        p = np.percentile(data, [2.5, 25, 50, 75, 97.5])
        ax.hlines(y, p[0], p[4], color="0.5", lw=1.0)
        ax.add_patch(Rectangle((p[1], y - 0.3), p[3] - p[1], 0.6, facecolor=color,
                               edgecolor="black", lw=0.5, alpha=0.65))
        ax.scatter([p[2]], [y], c=[color], s=38, edgecolors="black",
                   linewidths=0.6, zorder=4)


def minutes_axis(ax):
    # exact log-time -> minutes for a 3:00 runner: v -> 180*(exp(v)-1); the bars
    # live in v (log) units, so this only relabels the top axis (a nonlinear
    # minute scale), not the bar positions.
    sec = ax.secondary_xaxis(
        "top",
        functions=(lambda v: REF_MARATHON_MIN * np.expm1(v),
                   lambda m: np.log1p(m / REF_MARATHON_MIN)))
    sec.set_xlabel("minutes change for 3h marathon", fontsize=FS_AXIS)
    # minor ticks every 0.5 min (= 30 s) on the (minute-valued) secondary axis
    sec.xaxis.set_minor_locator(MultipleLocator(0.5))
    sec.tick_params(which="minor", length=2.5, color="0.5")


def whisker_caps(data):
    """Matplotlib's drawn whisker ends: most extreme points within 1.5*IQR of
    the quartiles. Pin count labels to these (not data.min/max) so they hug the
    box instead of floating out at the bootstrap extremes."""
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lo = data[data >= q1 - 1.5 * iqr].min()
    hi = data[data <= q3 + 1.5 * iqr].max()
    return lo, hi


def panel(ax, df, boot, ymap, colmap, style, rng):
    for t, grp in df.groupby("tier"):                 # alternating tier bands
        if t % 2:
            ys = [ymap[s] for s in grp["series_key"]]
            ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color="0.88", alpha=0.55,
                       zorder=0)
    ax.axvline(0, color="0.8", lw=1.0, ls="--", zorder=1)
    for r in df.itertuples():
        s = r.series_key
        data = boot[s].to_numpy(np.float64)
        draw_marker(ax, ymap[s], data, colmap.get(s, NEUTRAL), style, rng)
        cnt = f"({int(r.n_men)}/{int(r.n_women)})"
        wlo, whi = whisker_caps(data)
        if r.boot_p50 < 0:            # fast: count right of the right whisker cap
            ax.annotate(cnt, (whi, ymap[s]), xytext=(COUNT_PAD, 0),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=FS_COUNT, color="0.35")
        else:                         # slow: count left of the left whisker cap
            ax.annotate(cnt, (wlo, ymap[s]), xytext=(-COUNT_PAD, 0),
                        textcoords="offset points", ha="right", va="center",
                        fontsize=FS_COUNT, color="0.35")
    ax.set_yticks([])
    ax.set_xlabel("median $v_j$  (log-time)", fontsize=FS_AXIS)
    ax.grid(axis="x", alpha=0.25)
    minutes_axis(ax)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--style", choices=["boxplot", "strip", "ci_box"],
                    default="boxplot")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    apply_paper_style()
    dfA, bootA = load()

    n = len(dfA)
    yA = {s: i for i, s in enumerate(dfA["series_key"])}
    lblA = {r.series_key: f"({r.country}) {_short(r.series_key)} [{int(r.k_editions)}]"
            for r in dfA.itertuples()}

    # diverging colour scale = within-slice M-vs-W contrast
    allr = dfA["sex_contrast_min3h"].to_numpy()
    allr = allr[np.isfinite(allr)]
    vmax = max(np.nanpercentile(np.abs(allr), 95), 1e-6) if allr.size else 1.0
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    colA = {r.series_key: cmap(norm(r.sex_contrast_min3h))
            for r in dfA.itertuples() if np.isfinite(r.sex_contrast_min3h)}

    fig, axA = plt.subplots(1, 1, figsize=FIGSIZE)
    panel(axA, dfA, bootA, yA, colA, args.style, rng)
    axA.set_ylim(n - 0.5, -0.5)

    # race labels on the LEFT of the panel, right-aligned, hugging the spine
    trA = blended_transform_factory(axA.transAxes, axA.transData)
    for r in dfA.itertuples():
        axA.annotate(lblA[r.series_key], (LBL_LEFT, yA[r.series_key]), xycoords=trA,
                     xytext=(-LBL_PAD, 0), textcoords="offset points", ha="right",
                     va="center", fontsize=FS_RACE, annotation_clip=False)

    axA.set_xlim([-0.023, 0.031])

    # vertical colorbar dropped into the blank top-right corner (see CBAR_RECT)
    cax = axA.inset_axes(CBAR_RECT)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=cax, orientation="vertical")
    cb.set_label("M vs W contrast (min at 3:00 marathon)", fontsize=FS_AXIS)
    cb.set_ticks(np.linspace(-vmax, vmax, CBAR_NTICKS))
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    cb.ax.tick_params(labelsize=FS_RACE, length=2.5)
    cb.outline.set_linewidth(0.6)

    out = OUT_ROOT / "fig_cohort_single_forest.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"[write] {out} (+ .pdf)  (ALL {n} series, style={args.style})")

    out_paper = PAPER_FIG_DIR / "fig_cohort_single_forest.pdf"
    out_paper.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_paper, bbox_inches="tight")
    print(f"[write] {out_paper}")


if __name__ == "__main__":
    main()
