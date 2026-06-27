"""Paper figure for the 'spread of individual drift' paragraph (results 4.3).

Turns the dry d_i numbers into one two-panel figure off the production AxD point
fit (NO refit; d_i read straight from params["d"]):

  (A) Distribution of career-drift d_i -- the gauge-invariant headline. KDE per
      sex (cohort ALL, n_i >= floor), with +/-1 disattenuated-SD markers (the
      noise-corrected spread) and the right-skew annotated.

  (B) Real-runner trajectory fan. Holding ENTRY AGE FIXED to a single narrow band
      (default 35-40 yr, just above the masters threshold), we pick the REAL
      runners sitting at the p10/25/50/75/90 of d_i among well-sampled athletes
      (n_i >= floor) in that band, and plot each runner's ACTUAL finishes (markers)
      with the race effect v_j removed (an 'average race') plus the model's
      predicted trajectory (line). Everything is shown relative to each runner's
      own career-average performance, so the five are comparable and the overall
      ability level u_i drops out; markers scatter around their line by exactly
      the model residual.

Gauge note. Only the SPREAD of d_i is identified (the mean is pinned by the EB
prior). Panel A therefore reports shape only (SD, skew). In panel B the
per-runner de-meaning cancels the un-identified d_i-level shift, so the relative
trajectories + scatter are gauge-invariant.

Run::

    python scripts/05_analysis/di_distribution/p03_di_spread_fig.py
    python scripts/05_analysis/di_distribution/p03_di_spread_fig.py --fan-sex W --ae-lo 35 --ae-hi 40
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
from marathon_decomp.aging import (                               # noqa: E402
    aging_curve_on_grid,
    entry_age_curve_on_grid,
)
from marathon_decomp.config import PAPER_FIG_DIR                  # noqa: E402


# ===========================================================================
# ADJUSTABLE PARAMETERS
#
# Base look = the shared scripts/paper.mplstyle, loaded in apply_paper_style()
# exactly as the other p* scripts do. The knobs below retune THIS figure without
# touching the shared style. Font sizes split into (1) rcParam-backed text whose
# FS_* defaults MATCH paper.mplstyle (so the figure is unchanged out of the box
# -- edit to rescale) and (2) figure-specific text (panel tags, mini-table) that
# has no style default and is set per-artist. Each note cites the paper.mplstyle
# default it sets/overrides.
# ===========================================================================

# --- Canvas (paper.mplstyle figure.figsize default: 6.1 x 4.0) ---------------
FIGSIZE = (6.1, 3.5)         # two side-by-side panels

# --- Font sizes: (1) rcParam-backed (applied in apply_paper_style) -----------
FS_AXIS         = 9.0    # axis labels  | paper.mplstyle axes.labelsize: 9
FS_TICK         = 8.0    # tick labels  | paper.mplstyle xtick/ytick: 8
FS_LEGEND       = 6.0    # legend body  | paper.mplstyle legend.fontsize: 8 (departure)
FS_LEGEND_TITLE = 7.0    # legend title | paper.mplstyle legend.title_fontsize: 9 (departure)
# --- Font sizes: (2) figure-specific (custom; set per-artist) ----------------
FS_PANEL_TAG  = 12.0     # custom: the "(a)" / "(b)" panel tags
FS_TABLE_HEAD = 8.0      # custom: panel-A mini-table "men"/"women" headers
FS_TABLE      = 7.0      # custom: panel-A mini-table value rows

# --- Markers / lines (geometry, not style) -----------------------------------
KDE_LW     = 2.0         # panel-A KDE curve width
FAN_LW     = 1.6         # panel-B trajectory line (off-median percentiles)
FAN_LW_MED = 2.4         # panel-B trajectory line (median, highlighted)
FAN_MS     = 28          # panel-B finish-marker size (points^2)

# --- Legend padding (panel B; smaller = tighter box, masks less of the plot) -
LEG_HANDLETEXTPAD = 0.3  # gap marker <-> text   | matplotlib default 0.8
LEG_HANDLELENGTH  = 1.0  # legend handle width   | matplotlib default 2.0
LEG_BORDERPAD     = 0.3  # padding inside border | matplotlib default 0.4
LEG_LABELSPACING  = 0.25 # vertical gap, entries | matplotlib default 0.5

# --- Axis spacing / padding (rcParam-backed; applied in apply_paper_style) ----
AXES_LABELPAD = 2.0      # gap axis label <-> tick numbers | matplotlib default 4.0
TICK_PAD      = 1.5      # gap tick numbers <-> axis line  | matplotlib default 3.5

# --- Panel-A mini-table location (axes fraction; tune freely) ----------------
TBL_LAB_X = 0.23         # row-label column x (left-aligned)
TBL_MEN_X = 0.48         # men value column x (centred)
TBL_WOMEN_X = 0.64       # women value column x (centred)
TBL_TOP_Y = 0.20         # header ("men"/"women") row y
TBL_ROW_DY = 0.058       # vertical spacing between rows
# ===========================================================================

# --- Palette / markers (figure identity, not size tuning) --------------------
SEX_LABEL = {"M": "men", "W": "women"}
SEX_COLOR = {"M": "#2c7fb8", "W": "#d7301f"}
PCTLS = (10, 25, 50, 75, 90)
# diverging palette over the five percentiles (improver end -> decliner end)
PCTL_COLOR = {10: "#2166ac", 25: "#67a9cf", 50: "#444444",
              75: "#ef8a62", 90: "#b2182b"}
# non-directional markers only (no triangles): circle / square / diamond /
# pentagon / hexagon
PCTL_MARK = {10: "o", 25: "s", 50: "D", 75: "p", 90: "h"}


def apply_paper_style() -> None:
    """Load the shared paper.mplstyle (as the other p* scripts do), then push the
    tunable font knobs into rcParams.

    The FS_* rcParam defaults match paper.mplstyle, so the figure is unchanged
    out of the box; edit the constants above to rescale. Figure-specific text
    (panel tags, mini-table) is sized per-artist where it is drawn.
    """
    style_path = Path(__file__).resolve().parents[2] / "paper.mplstyle"
    plt.style.use(str(style_path) if style_path.exists() else "scripts/paper.mplstyle")
    plt.rcParams.update({                      # (paper.mplstyle default in comment)
        "axes.labelsize":        FS_AXIS,          # 9
        "xtick.labelsize":       FS_TICK,          # 8
        "ytick.labelsize":       FS_TICK,          # 8
        "legend.fontsize":       FS_LEGEND,        # 8
        "legend.title_fontsize": FS_LEGEND_TITLE,  # 9
        "axes.labelpad":         AXES_LABELPAD,    # matplotlib default 4.0
        "xtick.major.pad":       TICK_PAD,         # matplotlib default 3.5
        "ytick.major.pad":       TICK_PAD,         # matplotlib default 3.5
    })


# --------------------------------------------------------------------------- #
# panel A: distribution of d_i (spread + skew)                                #
# --------------------------------------------------------------------------- #
def panel_distribution(ax, cohort: str, floor: int) -> None:
    ax.grid(True, color="0.92", lw=0.6, zorder=0)
    ax.axvline(0.0, color="0.75", lw=0.8, ls="--", zorder=1)

    lo_hi, ymax = [], 0.0
    stat = {}     # sex -> (N, SD, skew)
    for sx in ("M", "W"):
        f = DI.load_cached(cohort, sx, 8.0, 2)
        if f is None:
            continue
        t = f.table(min_n=floor)
        d = t.d.to_numpy() * 100.0                       # %/yr
        pv = t.post_var.to_numpy()
        sd = np.sqrt(DI.disatt_var(t.d.to_numpy(), pv)) * 100.0   # disatt SD, %/yr
        sk = float(stats.skew(t.d.to_numpy()))
        stat[sx] = (len(d), sd, sk)
        lo_hi.append(np.percentile(d, [0.5, 99.5]))
        grid = np.linspace(d.min(), d.max(), 400)
        dens = stats.gaussian_kde(d)(grid)
        ymax = max(ymax, dens.max())
        ax.plot(grid, dens, color=SEX_COLOR[sx], lw=KDE_LW, zorder=3)

    if lo_hi:
        ax.set_xlim(min(p[0] for p in lo_hi), max(p[1] for p in lo_hi))
    ax.set_ylim(0.0, ymax * 1.10)
    ax.set_xlabel("career drift $d_i$  (%/yr;  $<0$ = improver)")
    ax.set_ylabel("density")
    ax.text(0.02, 0.98, "(a)", transform=ax.transAxes, ha="left", va="top",
            fontweight="bold", fontsize=FS_PANEL_TAG)

    # mini-table legend (item | men | women); position from the top-level
    # TBL_* constants so it is easy to tune
    ax.text(TBL_MEN_X, TBL_TOP_Y, "men", color=SEX_COLOR["M"], fontweight="bold",
            ha="center", transform=ax.transAxes, fontsize=FS_TABLE_HEAD)
    ax.text(TBL_WOMEN_X, TBL_TOP_Y, "women", color=SEX_COLOR["W"], fontweight="bold",
            ha="center", transform=ax.transAxes, fontsize=FS_TABLE_HEAD)
    rows = [
        ("N", lambda s: f"{stat[s][0]:,}"),
        ("SD (%/yr)", lambda s: f"{stat[s][1]:.2f}"),
        ("skew", lambda s: f"{stat[s][2]:+.2f}"),
    ]
    for r, (lab, fmt) in enumerate(rows):
        y = TBL_TOP_Y - (r + 1) * TBL_ROW_DY
        ax.text(TBL_LAB_X, y, lab, ha="left", transform=ax.transAxes, fontsize=FS_TABLE)
        ax.text(TBL_MEN_X, y, fmt("M"), ha="center", transform=ax.transAxes, fontsize=FS_TABLE)
        ax.text(TBL_WOMEN_X, y, fmt("W"), ha="center", transform=ax.transAxes, fontsize=FS_TABLE)


# --------------------------------------------------------------------------- #
# panel B: real-athlete trajectory fan at a fixed entry age                   #
# --------------------------------------------------------------------------- #
def _pick_pctl_athletes(d_sub: np.ndarray, sub_idx: np.ndarray) -> dict:
    """Real athlete (global index) nearest each d_i percentile, no repeats."""
    picks, used = {}, set()
    for p in PCTLS:
        target = np.percentile(d_sub, p)
        order = np.argsort(np.abs(d_sub - target))
        for k in order:
            if int(sub_idx[k]) not in used:
                picks[p] = int(sub_idx[k])
                used.add(int(sub_idx[k]))
                break
    return picks


def select_runners(f, floor: int, ae_lo: float, ae_hi: float) -> tuple[int, dict]:
    """(band size, {percentile -> global athlete index}) for the fixed-debut band."""
    sub = f.elig & (f.n_i >= floor) & np.isfinite(f.A_e) \
        & (f.A_e >= ae_lo) & (f.A_e <= ae_hi)
    sub_idx = np.flatnonzero(sub)
    return sub_idx.size, _pick_pctl_athletes(f.d[sub_idx], sub_idx)


def runner_traj(m, ae_i: float, d_i: float, an_bar: float, a) -> np.ndarray:
    """Predicted aging + g + drift contribution at career age(s) `a` (log units).

    The athlete level u_i is omitted; it cancels once the curve is anchored at
    debut (or de-meaned), which is how the panel and the deltas are displayed.
    """
    a = np.asarray(a, dtype=float)
    return (aging_curve_on_grid(m, a)
            + entry_age_curve_on_grid(m, a, ae_i - m.Ae_bar)
            + d_i * (a - an_bar))


def panel_fan(ax, f, picks: dict, n_sub: int, sex: str,
              ae_lo: float, ae_hi: float) -> None:
    if f is None:
        ax.set_axis_off()
        return
    m = f.model
    fd = f.fd
    resid = np.asarray(m.residuals(), dtype=float)       # y - fitted (log units)

    ax.grid(True, color="0.92", lw=0.6, zorder=0)
    ax.axhline(0.0, color="0.5", lw=0.8, ls="--", zorder=1)
    for p in PCTLS:
        ai = picks[p]
        rows = np.flatnonzero(fd.row_idx == ai)
        a_obs = fd.A_n[rows]
        ae_i = float(f.A_e[ai])
        d_i = float(f.d[ai])
        an_bar = float(a_obs.mean())                     # within-athlete centring

        def traj(a):
            return runner_traj(m, ae_i, d_i, an_bar, a)

        t0 = float(traj([0.0])[0])                       # debut level -> normalise to 0
        g = np.linspace(0.0, a_obs.max(), 100)
        y_line_obs = (traj(a_obs) - t0) * 100.0          # predicted at each finish
        y_obs = y_line_obs + resid[rows] * 100.0         # actual finish (line + resid)
        col, mk = PCTL_COLOR[p], PCTL_MARK[p]
        # faint stems tie each real finish to its predicted trajectory
        for xa, yl, yo in zip(a_obs, y_line_obs, y_obs):
            ax.plot([xa, xa], [yl, yo], color=col, lw=0.6, alpha=0.35, zorder=3)
        ax.plot(g, (traj(g) - t0) * 100.0, color=col,
                lw=FAN_LW_MED if p == 50 else FAN_LW,
                zorder=5 if p == 50 else 4)
        ax.scatter(a_obs, y_obs, s=FAN_MS, marker=mk, color=col, edgecolor="white",
                   linewidth=0.4, zorder=6,
                   label=f"p{p}: $d_i$={d_i * 100:+.2f}%/yr  (n={rows.size})")

    ax.set_xlabel("years since debut")
    ax.set_ylabel("finish time vs debut (%)   ($<0$ = faster)")
    ax.text(0.02, 0.98, "(b)", transform=ax.transAxes, ha="left", va="top",
            fontweight="bold", fontsize=FS_PANEL_TAG)
    leg = ax.legend(loc="lower left", framealpha=0.85, frameon=True,
                    edgecolor="none", title="$d_i$ percentile",
                    handletextpad=LEG_HANDLETEXTPAD, handlelength=LEG_HANDLELENGTH,
                    borderpad=LEG_BORDERPAD, labelspacing=LEG_LABELSPACING)
    leg.set_zorder(10)


# --------------------------------------------------------------------------- #
# console: the actual race records behind the five selected runners          #
# --------------------------------------------------------------------------- #
def _hms(sec: float) -> str:
    s = int(round(sec))
    return f"{s // 3600}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def print_records(f, picks: dict, sex: str, ae_lo: float, ae_hi: float) -> None:
    """Dump each selected runner's real finishes (console) and write an auditable
    per-runner career-delta summary (.md) backing the Figure-(b) numbers."""
    fd = f.fd
    m = f.model
    resid = np.asarray(m.residuals(), dtype=float)
    is_log = (fd.response_kind == "log_time")

    summ = []
    for p in PCTLS:
        ai = picks[p]
        rows = np.flatnonzero(fd.row_idx == ai)
        rows = rows[np.argsort(fd.race_date[fd.col_idx[rows]])]
        ae_i, d_i = float(f.A_e[ai]), float(f.d[ai])
        an_bar = float(fd.A_n[rows].mean())
        a_max = float(fd.A_n[rows].max())
        # predicted change debut->last career age in an average race (= Fig b line end)
        pct_end = float((runner_traj(m, ae_i, d_i, an_bar, [a_max])[0]
                         - runner_traj(m, ae_i, d_i, an_bar, [0.0])[0]) * 100.0)
        t_first = np.exp(fd.y[rows[0]]) if is_log else fd.y[rows[0]]
        t_last = np.exp(fd.y[rows[-1]]) if is_log else fd.y[rows[-1]]
        summ.append(dict(pctl=f"p{p}", aid=int(fd.athlete_ids[ai]), entry_age=ae_i,
                         d_i=d_i * 100.0, n=rows.size, career_yrs=a_max,
                         pred_change=pct_end, t_first=_hms(t_first), t_last=_hms(t_last)))

        print(f"\n=== p{p}  athlete_id={int(fd.athlete_ids[ai])}  "
              f"sex={fd.athlete_sex[ai]}  entry_age={ae_i:.1f}  "
              f"d_i={d_i * 100:+.2f}%/yr  n={rows.size} ===")
        print(f"  {'date':<11}{'race_id':>8}  {'series':<16}{'cty':<4}"
              f"{'yrs':>5}  {'actual':>9}{'pred':>9}  {'resid%':>7}")
        for r in rows:
            j = int(fd.col_idx[r])
            y = float(fd.y[r])
            yhat = y - float(resid[r])
            t_act = np.exp(y) if is_log else y
            t_pred = np.exp(yhat) if is_log else yhat
            d = pd.Timestamp(fd.race_date[j]).date()
            print(f"  {str(d):<11}{int(fd.race_ids[j]):>8}  "
                  f"{str(fd.race_series[j])[:15]:<16}{str(fd.race_country[j])[:3]:<4}"
                  f"{fd.A_n[r]:>5.1f}  {_hms(t_act):>9}{_hms(t_pred):>9}  "
                  f"{resid[r] * 100:>+7.2f}")

    # auditable markdown backing Figure (b)
    out = DI.OUT_ROOT / "fig_di_spread_records.md"
    lines = [
        f"# Figure (b) runners: {SEX_LABEL[sex]} debuting {ae_lo:.0f}-{ae_hi:.0f}",
        "",
        "Real athletes at each d_i percentile of the fixed-debut-age band (production "
        "AxD fit, no refit). `pred change debut->last` is the model's predicted "
        "finish-time change from debut to the runner's last career age in an average "
        "race (the endpoint of each line in Figure b); `<0` = faster. d_i x100 ~ %/yr.",
        "",
        "| pctl | athlete_id | entry age | d_i (%/yr) | n | career (yr) | "
        "pred change debut->last (%) | first finish | last finish |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for s in summ:
        lines.append(f"| {s['pctl']} | {s['aid']} | {s['entry_age']:.1f} | "
                     f"{s['d_i']:+.2f} | {s['n']} | {s['career_yrs']:.1f} | "
                     f"{s['pred_change']:+.1f} | {s['t_first']} | {s['t_last']} |")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="ALL")
    ap.add_argument("--fan-sex", default="M", choices=("M", "W"))
    ap.add_argument("--floor", type=int, default=5, help="n_i floor (default 5)")
    ap.add_argument("--ae-lo", type=float, default=35.0)
    ap.add_argument("--ae-hi", type=float, default=40.0)
    args = ap.parse_args()

    apply_paper_style()
    ff = DI.load_cached(args.cohort, args.fan_sex, 8.0, 2)
    n_sub, picks = (0, {}) if ff is None else select_runners(
        ff, args.floor, args.ae_lo, args.ae_hi)

    DI.OUT_ROOT.mkdir(parents=True, exist_ok=True)
    # records md first, so it is written even if a figure file is locked by a viewer
    if picks:
        print_records(ff, picks, args.fan_sex, args.ae_lo, args.ae_hi)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=FIGSIZE, constrained_layout=True)
    panel_distribution(axA, args.cohort, args.floor)
    panel_fan(axB, ff, picks, n_sub, args.fan_sex, args.ae_lo, args.ae_hi)

    stem = "fig_di_spread"
    for ext in ("png", "pdf"):
        out = DI.OUT_ROOT / f"{stem}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    out_paper = PAPER_FIG_DIR / f"{stem}.pdf"
    out_paper.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_paper, bbox_inches="tight")
    print(f"wrote {out_paper}")
    plt.close(fig)


if __name__ == "__main__":
    main()
