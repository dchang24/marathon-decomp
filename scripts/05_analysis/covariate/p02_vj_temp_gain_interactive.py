"""p02 - interactive (HTML) twin of p01: v_j vs air temp & total elevation gain.

Same two-panel scatter + joint 2D fit line as ``p01_vj_temp_gain.py``, but
rendered as a self-contained HTML page (Plotly loaded from the CDN, like
``scripts/visualizations/race_dashboard.py`` -- no Python ``plotly`` dependency;
the figure is embedded as JSON). Hovering a marker shows the exact race identity
and its raw covariates:

    series + year, v_j, air temperature, WBGT (field / max),
    total + net elevation gain, course id

The 11 highlighted series reuse p01's colours; markers are the
{circle, square, pentagon, hexagon, hexagon2, diamond} set, recycled (colour
disambiguates repeats). All others are one light-grey style.

Output: results/analysis/covariate/p02_vj_temp_gain_interactive/
        fig_vj_temp_gain_interactive__ALL_B.html
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import covariate_common as C
# reuse p01's fit + highlight styling (single source of truth)
from p01_vj_temp_gain import _HIGHLIGHT, fit_2d, YCOL, TEMP, GAIN

SUBDIR = "p02_vj_temp_gain_interactive"
GAIN_NET = "net_gain_m"
WBGT_F, WBGT_M, CID = "wbgt_field", "wbgt_max", "course_id"

# matplotlib {o,s,p,h,H,D} -> Plotly symbol names, recycled across HIGHLIGHT order
MARKERS_PLOTLY = ["circle", "square", "pentagon", "hexagon", "hexagon2", "diamond"]
REST_COLOR = "rgba(140,140,140,0.55)"

HOVER = ("<b>%{customdata[0]} %{customdata[1]}</b><br>"
         "v_j = %{customdata[2]}<br>"
         "air temp = %{customdata[3]} C<br>"
         "WBGT field/max = %{customdata[4]} C<br>"
         "total gain = %{customdata[5]} m (net %{customdata[6]})<br>"
         "course id = %{customdata[7]}<extra></extra>")

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>v_j vs weather + course (interactive)</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  html, body { margin: 0; height: 100%; font-family: 'Inter', 'Segoe UI', system-ui, sans-serif; }
  #chart { width: 100%; height: 100%; }
</style>
</head>
<body>
<div id="chart"></div>
<script>
  const FIG = __FIG_JSON__;
  Plotly.newPlot('chart', FIG.data, FIG.layout, {responsive: true});
</script>
</body>
</html>
"""


def _num(x, fmt: str) -> str:
    return fmt.format(x) if pd.notna(x) else "n/a"


def build_customdata(d: pd.DataFrame, highlight: dict) -> list[list[str]]:
    """Per-row hover payload (pre-formatted strings; index order matches HOVER)."""
    cd = []
    for _, r in d.iterrows():
        lab = highlight.get(r["series_key"], (None,))[0]
        disp = lab or str(r["series_key"]).replace("_marathon", "").replace("_", " ").title()
        cid = str(int(r[CID])) if pd.notna(r[CID]) else "n/a"
        cd.append([
            disp, str(int(r["year"])),
            _num(r[YCOL], "{:+.4f}"),
            _num(r[TEMP], "{:.1f}"),
            f"{_num(r[WBGT_F], '{:.1f}')} / {_num(r[WBGT_M], '{:.1f}')}",
            _num(r[GAIN], "{:.0f}"),
            _num(r.get(GAIN_NET, np.nan), "{:+.0f}"),
            cid,
        ])
    return cd


def marker_trace(d, xcol, *, color, symbol, name, xaxis, yaxis,
                 showlegend, legendgroup, highlight, size=9, line_w=0.5):
    return {
        "type": "scatter", "mode": "markers",
        "x": d[xcol].astype(float).tolist(),
        "y": d[YCOL].astype(float).tolist(),
        "xaxis": xaxis, "yaxis": yaxis,
        "marker": {"color": color, "symbol": symbol, "size": size,
                   "line": {"color": "black", "width": line_w}},
        "customdata": build_customdata(d, highlight),
        "hovertemplate": HOVER,
        "name": name, "legendgroup": legendgroup, "showlegend": showlegend,
    }


def line_trace(xs, ys, *, xaxis, yaxis, showlegend):
    return {
        "type": "scatter", "mode": "lines",
        "x": [float(v) for v in xs], "y": [float(v) for v in ys],
        "xaxis": xaxis, "yaxis": yaxis,
        "line": {"color": "black", "width": 2.0},
        "name": "2D fit (other var at mean)", "legendgroup": "fit",
        "showlegend": showlegend, "hoverinfo": "skip",
    }


def build_figure(df: pd.DataFrame) -> dict:
    beta, r2, m = fit_2d(df)
    b0, b_t, b_g = (float(x) for x in beta)
    d = df[m].copy()
    n = int(m.sum())
    t_mean, g_mean = float(d[TEMP].mean()), float(d[GAIN].mean())

    highlight = {k: (lab, col, MARKERS_PLOTLY[i % len(MARKERS_PLOTLY)])
                 for i, (k, lab, col) in enumerate(_HIGHLIGHT)}
    is_rest = ~d["series_key"].isin(highlight)
    rest = d[is_rest]

    # (axis pair, x column, line y(xs)) for the two panels
    panels = [
        ("x", "y", TEMP, lambda xs: b0 + b_t * xs + b_g * g_mean),
        ("x2", "y2", GAIN, lambda xs: b0 + b_t * t_mean + b_g * xs),
    ]

    data = []
    for p, (xaxis, yaxis, xcol, yfun) in enumerate(panels):
        first = p == 0                       # legend entries only on panel 1
        data.append(marker_trace(rest, xcol, color=REST_COLOR, symbol="circle",
                                  name="other races", xaxis=xaxis, yaxis=yaxis,
                                  showlegend=first, legendgroup="other",
                                  highlight=highlight, size=7, line_w=0.0))
        for key, (lab, col, sym) in highlight.items():
            sub = d[d["series_key"] == key]
            if sub.empty:
                continue
            data.append(marker_trace(sub, xcol, color=col, symbol=sym, name=lab,
                                      xaxis=xaxis, yaxis=yaxis, showlegend=first,
                                      legendgroup=key, highlight=highlight))
        xs = np.linspace(d[xcol].min(), d[xcol].max(), 100)
        data.append(line_trace(xs, yfun(xs), xaxis=xaxis, yaxis=yaxis,
                               showlegend=first))

    ylab = "v_j  (ALL_B full model, beta=0 gauge)  -  higher = harder"
    layout = {
        "title": {"text": f"Race factor v_j vs air temperature and total elevation "
                          f"gain  (joint 2D fit, R^2={r2:.3f}, n={n})",
                  "x": 0.5, "font": {"size": 15}},
        "hovermode": "closest",
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "font": {"family": "Inter, system-ui, sans-serif", "size": 12},
        "xaxis": {"domain": [0.0, 0.46], "anchor": "y",
                  "title": {"text": "air temperature (C)"}, "gridcolor": "#eee"},
        "xaxis2": {"domain": [0.54, 1.0], "anchor": "y2",
                   "title": {"text": "total elevation gain (m)"}, "gridcolor": "#eee"},
        "yaxis": {"domain": [0.0, 1.0], "anchor": "x",
                  "title": {"text": ylab}, "zeroline": False, "gridcolor": "#eee"},
        "yaxis2": {"domain": [0.0, 1.0], "anchor": "x2", "matches": "y",
                   "showticklabels": False, "zeroline": False, "gridcolor": "#eee"},
        "legend": {"title": {"text": "series"}, "x": 1.02, "y": 0.5,
                   "font": {"size": 10}},
        "margin": {"l": 80, "r": 140, "t": 60, "b": 60},
    }
    print(f"[fit] {YCOL} ~ temp + gain   R2={r2:.3f}   n={n}   "
          f"(traces={len(data)})")
    return {"data": data, "layout": layout}


def main() -> None:
    df = pd.read_parquet(C.MERGED_PATH)
    fig = build_figure(df)
    html = HTML.replace("__FIG_JSON__",
                        json.dumps(fig, allow_nan=False, default=lambda x: None))
    out = C.out_path(SUBDIR, "fig_vj_temp_gain_interactive", C.VAR_SELECT_SLICE, "html")
    out.write_text(html, encoding="utf-8")
    print(f"[write] {out}")
    print(f"  open: file://{out.resolve()}")


if __name__ == "__main__":
    main()
