"""p20 - draw the cohort connected forest from q20's processed data (compute is
in q20; p20 only reads + draws).

ALL (left) vs UK Po10 (right) series forests, native per-slice gauge, connectors
linking shared series. Rows ordered by the BOOTSTRAP MEDIAN (q20 rank), so the
boxplot centre line is monotone top->bottom. Colour = within-slice M-vs-W
contrast (red harder for men, blue harder for women, grey = none). Connectors
grey, dot-to-dot; thick = q04 ALL-vs-Po10 gauge-free mover, thin = raw shift.
In-panel (men/women) counts; gap labels: ALL "(country) name [editions]" right-
aligned at the dot, Po10 "name [editions]" left-aligned at the dot. See the
race_comparison README for the full visual spec and gauge handling.

Inputs (run q20 first): results/analysis/race_comparison/cohort_forest/
    cohort_forest_{ALL,Po10}_B.csv, _boot.parquet, cohort_forest_connectors.csv
Output: results/analysis/race_comparison/fig_cohort_connected_forest.{png,pdf}
    and paper/v1_final/figures/fig_cohort_connected_forest.pdf (paper copy)

Run::

    python scripts/05_analysis/race_comparison/p20_cohort_connected_forest.py
    python scripts/05_analysis/race_comparison/p20_cohort_connected_forest.py --style strip
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import MultipleLocator
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
# ADJUSTABLE PARAMETERS - tune everything about the figure here.
# ===========================================================================

# --- Figure size + horizontal spacing between the two panels -----------------
FIGSIZE = (16.5, 11.5)  # whole-figure (width, height) in inches
WSPACE = 0.45           # gap between the two panels (the middle connector zone);
                        #   SMALLER = panels closer together / narrower middle
LBL_A = 0.205           # ALL dot sits this far RIGHT of the left panel (axA axes
                        #   frac). This IS the gap between the left plot and its
                        #   race labels -> SMALLER = labels hug the left panel.
                        #   NB: axes-fraction units, so it scales with panel
                        #   width; shrinking WSPACE widens panels and inflates
                        #   this gap, so retune LBL_A after changing WSPACE.
LBL_B = 0.165           # Po10 label-zone width; Po10 dot this far LEFT of right
                        #   panel. Increase if long Po10 labels (e.g. "milton
                        #   keynes [11]") poke into the right panel.

# --- Font sizes (three text groups, independently tunable) -------------------
FS_AXIS = 8.0           # (1) axis labels, tick labels, colorbar label
FS_RACE = 5.5           # (2) middle race-name labels between the panels
FS_COUNT = 5.0          # (3) the (men/women) counts next to each boxplot
FS_TITLE = 10.0         # in-panel subplot title

# --- Markers, dots, connectors -----------------------------------------------
BOX_WIDTH = 0.62        # boxplot box height (data units = 1 row)
DOT_SIZE = 10           # connector dot marker size (points^2)
LBL_PAD = 7             # gap (points) between a race label and its dot
COUNT_PAD = 3           # gap (points) between a count and the whisker cap
CONN_LW = 1.0           # connector line width: thin = raw shift
CONN_LW_FLAG = 2.3      #                       thick = cohort-flagged mover
CONN_ALPHA = 0.5        # connector opacity: thin
CONN_ALPHA_FLAG = 0.95  #                    thick
CONN_COLOR = "0.55"     # connector colour: thin
CONN_COLOR_FLAG = "0.2" #                   thick (also the dot colour)
TITLE_Y = 0.995         # in-panel title vertical position (axes frac, 1 = top)

# --- Colorbar (horizontal, bottom) -------------------------------------------
CBAR_SHRINK = 0.5       # fraction of the available width the bar spans
CBAR_ASPECT = 55        # length / thickness (bigger = thinner bar)
CBAR_FRACTION = 0.04    # fraction of axes stolen for the bar
CBAR_PAD = 0.06         # gap between panels and the bar
# ===========================================================================


def apply_paper_style() -> None:
    """Load the shared paper.mplstyle, then shrink fonts for this dense figure."""
    style_path = Path(__file__).resolve().parents[2] / "paper.mplstyle"
    if style_path.exists():
        plt.style.use(str(style_path))
    else:                                              # fallback: run from root
        plt.style.use("scripts/paper.mplstyle")
    plt.rcParams.update({
        "font.size": FS_AXIS,
        "axes.titlesize": FS_TITLE,
        "axes.labelsize": FS_AXIS,
        "xtick.labelsize": FS_AXIS,
        "ytick.labelsize": FS_AXIS,
    })


def _short(s: str) -> str:
    return s.replace("_marathon", "").replace("_", " ")


def load():
    a = pd.read_csv(CF_DIR / "cohort_forest_ALL_B.csv")
    b = pd.read_csv(CF_DIR / "cohort_forest_Po10_B.csv")
    bootA = pd.read_parquet(CF_DIR / "cohort_forest_ALL_B_boot.parquet")
    bootB = pd.read_parquet(CF_DIR / "cohort_forest_Po10_B_boot.parquet")
    conn = pd.read_csv(CF_DIR / "cohort_forest_connectors.csv")
    return a, b, bootA, bootB, conn


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


def panel(ax, df, boot, ymap, colmap, style, title, rng):
    for t, grp in df.groupby("tier"):                 # alternating tier bands
        if t % 2:
            ys = [ymap[s] for s in grp["series_key"]]
            ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color="0.88", alpha=0.55,
                       zorder=0)
    ax.axvline(0, color="0.3", lw=1.0, ls="--", zorder=1)
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
    # subplot title sits INSIDE the panel (top-centre) to save vertical space
    ax.text(0.5, TITLE_Y, title, transform=ax.transAxes, ha="center", va="top",
            fontsize=FS_TITLE, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="0.7", alpha=0.85), zorder=7)
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
    dfA, dfB, bootA, bootB, conn = load()
    flagset = set(conn.loc[conn["cohort_flag"], "series_key"])

    nA, nB = len(dfA), len(dfB)
    nmax = max(nA, nB)
    offB = (nmax - nB) / 2.0
    yA = {s: i for i, s in enumerate(dfA["series_key"])}
    yB = {s: i + offB for i, s in enumerate(dfB["series_key"])}

    lblA = {r.series_key: f"({r.country}) {_short(r.series_key)} [{int(r.k_editions)}]"
            for r in dfA.itertuples()}
    lblB = {r.series_key: f"{_short(r.series_key)} [{int(r.k_editions)}]"
            for r in dfB.itertuples()}

    # shared diverging colour scale = within-slice M-vs-W contrast
    allr = np.concatenate([dfA["sex_contrast_min3h"].to_numpy(),
                           dfB["sex_contrast_min3h"].to_numpy()])
    allr = allr[np.isfinite(allr)]
    vmax = max(np.nanpercentile(np.abs(allr), 95), 1e-6) if allr.size else 1.0
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")

    def cmap_of(df):
        return {r.series_key: cmap(norm(r.sex_contrast_min3h))
                for r in df.itertuples() if np.isfinite(r.sex_contrast_min3h)}
    colA, colB = cmap_of(dfA), cmap_of(dfB)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=FIGSIZE,
                                   gridspec_kw=dict(wspace=WSPACE))
    panel(axA, dfA, bootA, yA, colA, args.style, "full dataset", rng)
    panel(axB, dfB, bootB, yB, colB, args.style, "UK Power of 10 subset", rng)
    for ax in (axA, axB):
        ax.set_ylim(nmax - 0.5, -0.5)

    trA = blended_transform_factory(axA.transAxes, axA.transData)
    trB = blended_transform_factory(axB.transAxes, axB.transData)

    # labels hug the dots: ALL right-aligned ENDING at its dot (label . ----),
    # Po10 left-aligned STARTING at its dot (---- . label). Dot = connector start.
    xa, xb = 1.0 + LBL_A, -LBL_B
    for r in dfA.itertuples():
        bold = "bold" if r.series_key in flagset else "normal"
        axA.annotate(lblA[r.series_key], (xa, yA[r.series_key]), xycoords=trA,
                     xytext=(-LBL_PAD, 0), textcoords="offset points", ha="right",
                     va="center", fontsize=FS_RACE, fontweight=bold, annotation_clip=False)
    for r in dfB.itertuples():
        bold = "bold" if r.series_key in flagset else "normal"
        axB.annotate(lblB[r.series_key], (xb, yB[r.series_key]), xycoords=trB,
                     xytext=(LBL_PAD, 0), textcoords="offset points", ha="left",
                     va="center", fontsize=FS_RACE, fontweight=bold, annotation_clip=False)

    shared = [s for s in dfA["series_key"] if s in yB]
    axA.scatter([xa] * len(shared), [yA[s] for s in shared], transform=trA,
                s=DOT_SIZE, color=CONN_COLOR_FLAG, clip_on=False, zorder=6)
    axB.scatter([xb] * len(shared), [yB[s] for s in shared], transform=trB,
                s=DOT_SIZE, color=CONN_COLOR_FLAG, clip_on=False, zorder=6)
    for s in shared:
        flagged = s in flagset
        con = ConnectionPatch(xyA=(xa, yA[s]), coordsA=trA, xyB=(xb, yB[s]),
                              coordsB=trB,
                              color=CONN_COLOR_FLAG if flagged else CONN_COLOR,
                              lw=CONN_LW_FLAG if flagged else CONN_LW,
                              alpha=CONN_ALPHA_FLAG if flagged else CONN_ALPHA,
                              zorder=5 if flagged else 2)
        fig.add_artist(con)

    axA.set_xlim([-0.023, 0.031])
    axB.set_xlim([-0.026, 0.038])

    # horizontal colorbar along the bottom so the two panels can span full width
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=[axA, axB], orientation="horizontal",
                      fraction=CBAR_FRACTION, pad=CBAR_PAD, aspect=CBAR_ASPECT,
                      shrink=CBAR_SHRINK)
    cb.set_label("M vs W contrast (min at 3:00 marathon)", fontsize=FS_AXIS)
    cb.ax.tick_params(labelsize=FS_AXIS)

    out = OUT_ROOT / "fig_cohort_connected_forest.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"[write] {out} (+ .pdf)  (ALL {nA} / Po10 {nB}, {len(shared)} shared, "
          f"{len(flagset)} cohort-flagged, style={args.style})")

    #out_paper = PAPER_FIG_DIR / "fig_cohort_connected_forest.pdf"
    #out_paper.parent.mkdir(parents=True, exist_ok=True)
    #fig.savefig(out_paper, bbox_inches="tight")
    #print(f"[write] {out_paper}")


if __name__ == "__main__":
    main()
