"""p01 - v_j (ALL_B full model) vs air temp & total elevation gain, 2-panel scatter.

Left panel : v_j vs air temp (temp_field, the q02 weather winner).
Right panel: v_j vs total elevation gain (total_gain_m).

Both panels overlay the **single 2D regression** fit
    v_j ~ b0 + b_temp * temp + b_gain * total_gain
fitted once on the complete cases. The line drawn in each panel is that fit with
the OTHER predictor held at its sample mean (the partial/marginal effect), so the
two lines come from one joint model, not two separate univariate fits.

The N races with the largest |residual| from the joint fit are annotated with a
small text label (series + 2-digit year) in both panels.

Eleven series are individually highlighted: the seven Abbott World Marathon
Majors (Tokyo, Boston, London, Berlin, Chicago, New York, Sydney) plus Valencia,
Paris, Stockholm, Amsterdam. All others share one light-grey style. Markers are
restricted to {o, s, p, h, H, D} and recycled (colour disambiguates repeats).

Uses the production ALL_B full-model beta=0-gauged v_j.
Output: results/analysis/covariate/p01_vj_temp_gain/fig_vj_temp_gain__ALL_B.png
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import covariate_common as C
from marathon_decomp.config import PAPER_FIG_DIR

SUBDIR = "p01_vj_temp_gain"
YCOL = f"v_{C.VAR_SELECT_SLICE}"     # ALL_B full v_j
TEMP, GAIN = "temp_field", "total_gain_m"

LAYOUT = {
    "figsize": (6.1, 4.2),           # figure size in inches (default paper width, slightly taller for legend)
    "marker_size_other": 15,         # scatter size for non-highlighted races
    "marker_size_highlight": 15,     # scatter size for highlighted races
    "outlier_fontsize": 5.0,         # fontsize for the outlier annotations
    "y_labelpad": 2,                 # padding for the left y-axis label
    
    # Legend settings
    "legend_bbox_y": 0.0,            # anchor y-position for legend (0 = bottom)
    "legend_ncol": 5,                # number of columns in legend
    "legend_fontsize": 8,            # fontsize for legend text
    
    # Subplot margins
    "margin_left": 0.09,
    "margin_right": 0.98,
    "margin_bottom": 0.28,           # bottom margin (increase if legend overlaps x-axis)
    "margin_top": 0.92,
    "wspace": 0.03,                  # horizontal space between the two subplots
}

# Markers are restricted to this set and recycled (colour disambiguates the
# repeats); cycled in HIGHLIGHT order below.
MARKERS = ["o", "s", "p", "h", "H", "D"]

# (series_key, label, colour) for the highlighted races; marker assigned by
# cycling MARKERS in order.
_HIGHLIGHT = [
    # --- 7 Abbott World Marathon Majors ---
    ("tokyo_marathon",   "Tokyo",     "#d62728"),
    ("boston_marathon",  "Boston",    "#1f77b4"),
    ("london_marathon",  "London",    "#2ca02c"),
    ("berlin_marathon",  "Berlin",    "#9467bd"),
    ("chicago_marathon", "Chicago",   "#ff7f0e"),
    ("nyc_marathon",     "New York",  "#8c564b"),
    ("sydney_marathon",  "Sydney",    "#e377c2"),
    # --- four additional flagged races ---
    ("valencia_marathon",  "Valencia",  "#17becf"),
    ("paris_marathon",     "Paris",     "#bcbd22"),
    ("stockholm_marathon", "Stockholm", "#393b79"),
    ("amsterdam_marathon", "Amsterdam", "#7f7f7f"),
]
HIGHLIGHT = [(k, lab, col, MARKERS[i % len(MARKERS)])
             for i, (k, lab, col) in enumerate(_HIGHLIGHT)]
REST_STYLE = dict(color="0.55", marker="o", s=LAYOUT["marker_size_other"], alpha=0.85,
                  edgecolors="none", zorder=1)


def fit_2d(df: pd.DataFrame):
    """Raw-units OLS v_j ~ temp + gain on complete cases. Returns (beta, r2, mask)."""
    y = df[YCOL].to_numpy(float)
    Xt = df[TEMP].to_numpy(float)
    Xg = df[GAIN].to_numpy(float)
    m = np.isfinite(y) & np.isfinite(Xt) & np.isfinite(Xg)
    X = np.column_stack([np.ones(m.sum()), Xt[m], Xg[m]])
    beta, _, r2 = C.ols(X, y[m])
    return beta, r2, m


N_OUTLIERS = 5                      # races labelled (largest +residual ABOVE fit)

# Per-label annotation overrides for outliers whose default top-right offset runs
# off the panel. Keyed by series_key -> kwargs passed to ax.annotate (xytext is in
# offset points; ha/va control the anchor). Anything not listed uses DEFAULT_ANNO.
DEFAULT_ANNO = dict(xytext=(3, 3), ha="left", va="baseline")
LABEL_OVERRIDES = {
    "milton_keynes_marathon": dict(xytext=(-3, 3), ha="right", va="baseline"),
    "stockholm_marathon":     dict(xytext=(8, -4), ha="center", va="top"),
    "london_marathon":     dict(xytext=(-1, -4), ha="right", va="top"),
}


def short_label(row: pd.Series) -> str:
    """Compact race tag for an outlier annotation, e.g. 'milton keynes 18'."""
    name = str(row["series_key"]).replace("_marathon", "").replace("_", " ")
    return f"{name} {int(row['year']) % 100:02d}"


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    beta, r2, m = fit_2d(df)
    b0, b_t, b_g = beta
    d = df[m].copy()
    t_mean, g_mean = d[TEMP].mean(), d[GAIN].mean()
    print(f"[fit] {YCOL} ~ temp + gain   R2={r2:.3f}   "
          f"b_temp={b_t:+.5f} log/C   b_gain={b_g:+.6f} log/m   n={int(m.sum())}")

    # residuals from the joint fit -> the N largest |resid| races get labelled
    X_all = np.column_stack([np.ones(len(d)), d[TEMP].to_numpy(float),
                             d[GAIN].to_numpy(float)])
    d["_resid"] = d[YCOL].to_numpy(float) - X_all @ beta
    # only the series ABOVE the trendline (positive residual) get labelled
    above = d[d["_resid"] > 0]
    outliers = above.loc[above["_resid"].sort_values(ascending=False).index[:N_OUTLIERS]]
    print(f"[outliers] labelled {len(outliers)}: "
          + ", ".join(f"{short_label(r)} ({r._resid:+.3f})"
                      for _, r in outliers.iterrows()))

    highlight_keys = {k for k, *_ in HIGHLIGHT}
    is_rest = ~d["series_key"].isin(highlight_keys)

    import pathlib
    style_path = pathlib.Path(__file__).resolve().parents[3] / "scripts" / "paper.mplstyle"
    if style_path.exists():
        plt.style.use(str(style_path))
    else:
        plt.style.use("scripts/paper.mplstyle")

    fig, axes = plt.subplots(1, 2, figsize=LAYOUT["figsize"], sharey=True)
    # each panel holds the OTHER predictor at its mean: X0 row = [1, temp, gain]
    panels = [
        (axes[0], TEMP, "air temperature (°C)",
         lambda xs: np.column_stack([np.ones_like(xs), xs, np.full_like(xs, g_mean)])),
        (axes[1], GAIN, "total elevation gain (m)",
         lambda xs: np.column_stack([np.ones_like(xs), np.full_like(xs, t_mean), xs])),
    ]

    for ax, xcol, xlabel, x0_fn in panels:
        ax.scatter(d.loc[is_rest, xcol], d.loc[is_rest, YCOL],
                   label="other races", **REST_STYLE)
        for key, lab, col, mk in HIGHLIGHT:
            sub = d[d["series_key"] == key]
            if sub.empty:
                continue
            ax.scatter(sub[xcol], sub[YCOL], color=col, marker=mk, s=LAYOUT["marker_size_highlight"],
                       edgecolors="black", linewidths=0.4, zorder=3, label=lab)

        # fitted (mean) line, other predictor held at its sample mean
        xs = np.linspace(d[xcol].min(), d[xcol].max(), 100)
        yhat = x0_fn(xs) @ beta
        ax.plot(xs, yhat, color="black", lw=2.0, zorder=4,
                label="2D fit (other var at mean)")

        # small text labels next to the outlier races (only on the weather subplot)
        if ax == axes[0]:
            for _, r in outliers.iterrows():
                anno = LABEL_OVERRIDES.get(r["series_key"], DEFAULT_ANNO)
                ax.annotate(short_label(r), (r[xcol], r[YCOL]),
                            textcoords="offset points",
                            fontsize=LAYOUT["outlier_fontsize"], color="0.15",
                            zorder=5, **anno)

        ax.set_xlabel(xlabel)
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel(r"race factor $v_j$  (higher = harder)", labelpad=LAYOUT["y_labelpad"])
    axes[0].set_title("(a) Weather")
    axes[1].set_title("(b) Course")

    axes[0].set_xlim([2, 27])
    axes[1].set_xlim([0, 450]) 

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, LAYOUT["legend_bbox_y"]),
               ncol=LAYOUT["legend_ncol"], frameon=True, fontsize=LAYOUT["legend_fontsize"], title="marathon",
               borderpad=0.3, labelspacing=0.3, columnspacing=0.8, handletextpad=0.4)
    
    fig.subplots_adjust(
        left=LAYOUT["margin_left"], right=LAYOUT["margin_right"], 
        bottom=LAYOUT["margin_bottom"], top=LAYOUT["margin_top"], 
        wspace=LAYOUT["wspace"]
    )

    out = C.out_path(SUBDIR, "fig_vj_temp_gain", C.VAR_SELECT_SLICE, "png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"[write] {out}")

    out_pdf = C.out_path(SUBDIR, "fig_vj_temp_gain", C.VAR_SELECT_SLICE, "pdf")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"[write] {out_pdf}")

    out_paper = PAPER_FIG_DIR / "fig_vj_temp_gain.pdf"
    out_paper.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_paper, bbox_inches="tight")
    print(f"[write] {out_paper}")


if __name__ == "__main__":
    main()
