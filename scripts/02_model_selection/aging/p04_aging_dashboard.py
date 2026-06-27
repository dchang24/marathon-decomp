"""Interactive aging-form dashboard -- one self-contained HTML for every slice.

Recomputes nothing: reads e01's ``grid/curves.parquet`` + ``grid/an_density.parquet``
+ ``grid/slice_info.csv`` + ``grid/metrics.csv`` and e02's ``cv/form_selection.csv``
for every slice with e01 output, and bundles them into a single Plotly page.

What you get
------------
* **Left panel** -- the reconstructed curve (aging block + entry-age gamma term)
  in log-time for every parametric form of the chosen slice + nu. Colour = basis
  family, dash = gamma form (off=dotted, scalar=dashed, varying=solid);
  ``spline5-gvarying`` is drawn bold as the reference. The grey rug is the A_n
  density and the dotted vertical line marks the observed-A_n p95 (extrapolation
  zone to its right). Spline gvarying curves carry dot markers at their knot
  positions. Toggle any form on/off via its legend entry, the gamma/basis quick
  filters, or select-all / unselect-all -- and the on/off set PERSISTS when you
  change slice / nu / entry age. Axis limits are fixed per slice so curves never
  shift as you toggle.
* **Entry-age toggle** -- the gamma block is linear in entry age, so the curve is
  shown at four presets (35 / 45 / 55 / 65 yr). gamma=off curves are entry-age
  invariant.
* **Right panel** -- the model-comparison table for the same slice + nu:
  in-sample loglik / AIC / BIC / n_params (best-init, from the grid) and held-out
  CV log-density per cell + RMSE (from the CV stage). The CV-best form is flagged,
  and the best value in each column is bold.
* **Two-way linking** -- click a curve (or its legend, or a table row) and the
  matching table row + curve are highlighted together.

Controls: slice dropdown, nu radio (inf / 8), gamma + basis quick filters.

Output -> results/model_selection/aging/aging_dashboard.html  (open in a browser;
needs network access for the Plotly CDN, matching scripts/visualizations).

Run::

    python scripts/02_model_selection/aging/p04_aging_dashboard.py
    python scripts/02_model_selection/aging/p04_aging_dashboard.py --slices ALL_M Po10_M
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
OUT_HTML = OUT_ROOT / "aging_dashboard.html"
REFERENCE_CAND = "spline5-gvarying"   # form drawn bold on the curve panel
ENTRY_AGES = (35.0, 45.0, 55.0, 65.0)  # entry-age presets shown by the toggle
DEFAULT_ENTRY_AGE = 45.0


# -- helpers ---------------------------------------------------------------

def _nukey(x) -> str:
    """Canonical nu label shared across files: 'inf' or e.g. '8'."""
    s = str(x)
    if "inf" in s.lower():
        return "inf"
    return f"{float(x):g}"


def _basis_sort_key(name: str) -> tuple[int, int]:
    kind = 0 if name.startswith("poly") else 1
    n = int("".join(ch for ch in name if ch.isdigit()) or 0)
    return (kind, n)


def _cand_sort_key(cand: str) -> tuple:
    basis, _, g = cand.partition("-g")
    gamma_order = {"off": 0, "scalar": 1, "varying": 2}.get(g, 9)
    return (*_basis_sort_key(basis), gamma_order)


def _r(x, nd):
    """Round to nd places, mapping non-finite to None (JSON null)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return round(v, nd) if math.isfinite(v) else None


def _best_init_metrics(grid_dir: Path) -> pd.DataFrame:
    """Best-init (max loglik) in-sample row per (nu, cand) from grid/metrics.csv."""
    m = pd.read_csv(grid_dir / "metrics.csv")
    m["nu"] = m["nu"].map(_nukey)
    idx = m.groupby(["nu", "cand"])["loglik"].idxmax()
    keep = ["nu", "cand", "basis", "gamma_form", "loglik", "aic", "bic", "n_params"]
    return m.loc[idx, keep].reset_index(drop=True)


def _spline_knots(grid_dir: Path, cand: str) -> list | None:
    """Spline knot positions (A_n) for a saved fit, from its pickle's model_extra."""
    hits = sorted((grid_dir / "fits").glob(f"{cand}_*.pkl"))
    if not hits:
        return None
    with hits[0].open("rb") as f:
        payload = pickle.load(f)
    k = payload.get("model_extra", {}).get("spline_knots")
    return [_r(x, 3) for x in np.asarray(k)] if k is not None else None


def _cv_table(cv_dir: Path) -> pd.DataFrame | None:
    p = cv_dir / "form_selection.csv"
    if not p.is_file():
        return None
    f = pd.read_csv(p)
    f["nu"] = f["nu"].map(_nukey)
    return f[["nu", "cand", "cv_per_cell", "cv_logdens", "heldout_rmse",
              "n_test_total", "n_orphan_total"]]


# -- per-slice payload ------------------------------------------------------

def build_slice(slug: str) -> dict | None:
    grid_dir = OUT_ROOT / slug / "grid"
    if not (grid_dir / "curves.parquet").is_file():
        return None
    curves = pd.read_parquet(grid_dir / "curves.parquet")
    dens = pd.read_parquet(grid_dir / "an_density.parquet")
    info = pd.read_csv(grid_dir / "slice_info.csv").iloc[0]
    insample = _best_init_metrics(grid_dir)
    cv = _cv_table(OUT_ROOT / slug / "cv")

    curves["nu"] = curves["nu"].map(_nukey)
    nus = sorted(curves["nu"].unique(), key=lambda x: (x != "inf", x))

    # shared A_n grid (identical across nu/cand within a slice)
    a_grid = (curves[curves.nu == nus[0]].sort_values("A_n")
              .drop_duplicates("A_n")["A_n"].to_numpy())
    A_n = [_r(a, 3) for a in a_grid]

    cands = sorted(curves["cand"].unique(), key=_cand_sort_key)

    # merge in-sample + cv into a (nu, cand) -> metrics lookup
    tbl = insample.copy()
    if cv is not None:
        tbl = tbl.merge(cv, on=["nu", "cand"], how="left")

    # entry-age axis: the gamma block is exactly linear in (A_e - ae_mean), so the
    # stored gamma curve (computed at +ae_off yr) scales to any entry age. We carry
    # the aging curve + the unit gamma curve and combine in the browser.
    ae_mean = float(info.get("ae_mean")) if math.isfinite(float(info.get("ae_mean"))) else float("nan")
    ae_off = float(curves["ae_offset_yr"].dropna().iloc[0]) if curves["ae_offset_yr"].notna().any() else 10.0
    DISPLAY_AGES = list(ENTRY_AGES)

    curves_out: dict[str, dict[str, list]] = {}
    gamma_out: dict[str, dict[str, list | None]] = {}
    table_out: dict[str, dict[str, dict]] = {}
    ymin = ymax = None      # fixed y-range over ALL forms, nu, and display entry ages
    for nu in nus:
        cn = curves[curves.nu == nu]
        cn_aging = {c: g.sort_values("A_n")["aging_curve"].to_numpy()
                    for c, g in cn.groupby("cand")}
        cn_gamma = {c: g.sort_values("A_n")["gamma_curve"].to_numpy()
                    for c, g in cn.groupby("cand")}
        curves_out[nu] = {c: [_r(y, 6) for y in cn_aging[c]]
                          for c in cands if c in cn_aging}
        gamma_out[nu] = {}
        for c in cands:
            if c not in cn_aging:
                continue
            gv = cn_gamma[c]
            has_g = np.isfinite(gv).any()
            gamma_out[nu][c] = [_r(y, 6) for y in gv] if has_g else None
            # y-range over the displayed entry ages (gamma scales linearly)
            for age in DISPLAY_AGES:
                y = cn_aging[c] + (gv * ((age - ae_mean) / ae_off) if has_g else 0.0)
                _gmin = float(np.nanmin(y)); _gmax = float(np.nanmax(y))
                ymin = _gmin if ymin is None else min(ymin, _gmin)
                ymax = _gmax if ymax is None else max(ymax, _gmax)
        tn = tbl[tbl.nu == nu].set_index("cand")
        rows: dict[str, dict] = {}
        for c in cands:
            if c not in tn.index:
                continue
            row = tn.loc[c]
            rows[c] = dict(
                basis=str(row["basis"]), gamma=str(row["gamma_form"]),
                loglik=_r(row.get("loglik"), 1), aic=_r(row.get("aic"), 1),
                bic=_r(row.get("bic"), 1), n_params=int(row["n_params"]),
                cv_per_cell=_r(row.get("cv_per_cell"), 5),
                cv_logdens=_r(row.get("cv_logdens"), 1),
                heldout_rmse=_r(row.get("heldout_rmse"), 5),
            )
        table_out[nu] = rows

    # fixed axis ranges (5% y-pad) so curves don't shift on toggle / nu switch
    xmax = float(info.get("an_max")) if math.isfinite(float(info.get("an_max"))) else max(a_grid)
    if ymin is None:
        ymin, ymax = -1.0, 1.0
    pad = 0.05 * (ymax - ymin + 1e-9)
    yrange = [_r(ymin - pad, 5), _r(ymax + pad, 5)]
    xrange = [_r(-0.02 * xmax, 3), _r(xmax, 3)]

    # spline knot positions (gvarying forms only) -> dot markers on the curve
    knots = {}
    for c in cands:
        if c.startswith("spline") and c.endswith("gvarying"):
            kp = _spline_knots(grid_dir, c)
            if kp is not None:
                knots[c] = kp

    return dict(
        slug=slug, nus=nus, cands=cands, A_n=A_n, yrange=yrange, xrange=xrange,
        an_p95=_r(info.get("an_p95"), 3), an_max=_r(info.get("an_max"), 3),
        ae_mean=_r(ae_mean, 2), ae_off=_r(ae_off, 4),
        I=int(info["I"]), J=int(info["J"]), N=int(info["N"]),
        rug=dict(lo=[_r(x, 3) for x in dens["bin_lo"]],
                 hi=[_r(x, 3) for x in dens["bin_hi"]],
                 count=[int(c) for c in dens["count"]]),
        curves=curves_out, gamma=gamma_out, table=table_out, knots=knots,
    )


# -- HTML emission ----------------------------------------------------------

def render(payload: dict) -> str:
    data_json = json.dumps(payload, separators=(",", ":"))
    return _TEMPLATE.replace("/*__DATA__*/", data_json)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"],
                    help="slugs / slice names, or 'all' (every slice with e01 output).")
    args = ap.parse_args()

    all_slugs = sorted(p.parent.parent.name
                       for p in OUT_ROOT.glob("*/grid/curves.parquet"))
    if args.slices == ["all"]:
        slugs = all_slugs
    else:
        slugs = []
        for s in args.slices:
            slugs.extend([g for g in all_slugs if g == s or g.startswith(s)] or [])
        slugs = sorted(dict.fromkeys(slugs))
    if not slugs:
        print("No e01 grid/curves.parquet found -- run e01_aging_grid first.")
        return

    slices = {}
    for slug in slugs:
        print(f"  packing {slug}")
        d = build_slice(slug)
        if d is not None:
            slices[slug] = d
    if not slices:
        print("No usable slices.")
        return

    payload = dict(slices=slices, slug_order=list(slices.keys()),
                   reference=REFERENCE_CAND,
                   entry_ages=list(ENTRY_AGES), default_entry_age=DEFAULT_ENTRY_AGE)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(render(payload), encoding="utf-8")
    print(f"\nWrote {OUT_HTML}  ({len(slices)} slices, "
          f"{OUT_HTML.stat().st_size / 1024:.0f} KB)")


_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Aging-form dashboard</title>
<meta name="description" content="Interactive aging-curve dashboard with per-form metrics (loglik / AIC / BIC / CV) by slice and nu">
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  :root {
    --bg:#f0f2f5; --surface:#fff; --surface-alt:#f8f9fb;
    --border:#e2e5ea; --border-light:#eef0f3;
    --text:#1a1d23; --text-muted:#6b7280; --text-dim:#9ca3af;
    --accent:#4f46e5; --accent-light:#eef2ff; --accent-border:#c7d2fe;
    --green:#059669; --green-bg:#ecfdf5; --green-border:#a7f3d0;
    --radius:6px;
    --font:'Inter','Segoe UI',system-ui,-apple-system,sans-serif;
    --mono:'JetBrains Mono','Cascadia Code','Consolas',monospace;
  }
  * { box-sizing:border-box; margin:0; padding:0; }
  html,body { height:100%; overflow:hidden; }
  body { font-family:var(--font); font-size:13px; color:var(--text);
         background:var(--bg); display:flex; flex-direction:column; }

  .header { background:var(--surface); border-bottom:1px solid var(--border);
            padding:10px 16px; display:flex; flex-wrap:wrap; align-items:flex-end;
            gap:10px 20px; flex-shrink:0; box-shadow:0 1px 2px rgba(0,0,0,.05); }
  .header-title h1 { font-size:15px; font-weight:700; color:var(--accent); letter-spacing:-0.3px; }
  .header-title .subtitle { font-size:10px; color:var(--text-dim); text-transform:uppercase;
            font-weight:600; letter-spacing:0.5px; margin-top:1px; }
  .controls { display:flex; flex-wrap:wrap; gap:8px 16px; align-items:flex-end; }
  .ctrl-group { display:flex; flex-direction:column; gap:3px; }
  .ctrl-label { font-size:9px; font-weight:700; text-transform:uppercase;
            letter-spacing:0.6px; color:var(--text-dim); }
  .toggle-row { display:flex; gap:3px; }
  .toggle-btn { font-family:var(--font); font-size:11px; font-weight:600;
            padding:3px 10px; border-radius:4px; cursor:pointer;
            border:1.5px solid #d1d5db; background:#f3f4f6; color:#6b7280;
            transition:all .15s ease; user-select:none; }
  .toggle-btn:hover { border-color:var(--accent); color:var(--accent); }
  .toggle-btn.active { background:var(--accent-light); color:var(--accent);
            border-color:var(--accent-border); }
  select { font-family:var(--font); font-size:11px; font-weight:500;
            border:1.5px solid var(--border); border-radius:4px; padding:3px 6px;
            background:var(--surface-alt); color:var(--text); cursor:pointer; }
  select:hover { border-color:var(--accent); }
  .slice-info { margin-left:auto; font-size:11px; color:var(--text-muted);
            background:var(--accent-light); border:1px solid var(--accent-border);
            border-radius:var(--radius); padding:6px 12px; }
  .slice-info b { color:var(--accent); font-family:var(--mono); }

  .main { flex:1; display:flex; min-height:0; }
  .left { flex:1.55; min-width:0; padding:8px; }
  #chart { width:100%; height:100%; }
  .right { flex:1; min-width:380px; border-left:1px solid var(--border);
            background:var(--surface); display:flex; flex-direction:column; min-height:0; }
  .right-head { padding:10px 14px 6px; border-bottom:1px solid var(--border-light); }
  .right-head h2 { font-size:13px; font-weight:700; }
  .right-head p { font-size:10.5px; color:var(--text-muted); margin-top:2px; }
  .table-wrap { overflow:auto; flex:1; }
  table { border-collapse:collapse; width:100%; font-size:11.5px; }
  thead th { position:sticky; top:0; background:var(--surface-alt); z-index:2;
            text-align:right; padding:6px 8px; font-weight:700; color:var(--text-muted);
            border-bottom:1.5px solid var(--border); white-space:nowrap; cursor:pointer; }
  thead th.lft { text-align:left; }
  thead th:hover { color:var(--accent); }
  tbody td { padding:4px 8px; text-align:right; border-bottom:1px solid var(--border-light);
            font-family:var(--mono); white-space:nowrap; }
  tbody td.lft { text-align:left; font-family:var(--font); }
  tbody tr { cursor:pointer; }
  tbody tr:hover { background:var(--accent-light); }
  tbody tr.sel { background:var(--accent-light); box-shadow:inset 3px 0 0 var(--accent); }
  tbody tr.cvbest td.lft::after { content:" * CV-best"; color:var(--green);
            font-size:9px; font-weight:700; }
  td.best { font-weight:700; color:var(--accent); }
  .swatch { display:inline-block; width:10px; height:10px; border-radius:2px;
            margin-right:6px; vertical-align:middle; }
  .muted { color:var(--text-dim); }
  .foot { font-size:9.5px; color:var(--text-dim); padding:6px 14px;
            border-top:1px solid var(--border-light); }
</style>
</head>
<body>
<div class="header">
  <div class="header-title">
    <h1>Aging-form dashboard</h1>
    <div class="subtitle">curve shape &amp; model comparison by slice</div>
  </div>
  <div class="controls">
    <div class="ctrl-group">
      <span class="ctrl-label">Slice</span>
      <select id="sel-slice"></select>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">nu</span>
      <div class="toggle-row" id="nu-row"></div>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">entry age</span>
      <div class="toggle-row" id="age-row"></div>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">gamma</span>
      <div class="toggle-row" id="gamma-row">
        <span class="toggle-btn" data-g="off">off</span>
        <span class="toggle-btn active" data-g="scalar">scalar</span>
        <span class="toggle-btn" data-g="varying">varying</span>
      </div>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">basis</span>
      <div class="toggle-row" id="basis-row">
        <span class="toggle-btn active" data-b="poly">poly</span>
        <span class="toggle-btn active" data-b="spline">spline</span>
      </div>
    </div>
    <div class="ctrl-group">
      <span class="ctrl-label">all forms</span>
      <div class="toggle-row" id="all-row">
        <span class="toggle-btn" id="btn-all">select all</span>
        <span class="toggle-btn" id="btn-none">unselect all</span>
      </div>
    </div>
  </div>
  <div class="slice-info" id="slice-info"></div>
</div>

<div class="main">
  <div class="left"><div id="chart"></div></div>
  <div class="right">
    <div class="right-head">
      <h2>Model comparison <span id="tbl-nu" class="muted"></span></h2>
      <p>In-sample loglik / AIC / BIC (best init) &amp; held-out CV. Click a row or
         curve to link them. Bold = best in column; AIC/BIC lower is better,
         loglik &amp; CV/cell higher is better.</p>
    </div>
    <div class="table-wrap"><table id="tbl">
      <thead></thead><tbody></tbody>
    </table></div>
    <div class="foot">Toggle forms via the legend, gamma/basis quick filters, or
      select/unselect all -- the on/off set persists across slice / nu / entry-age
      changes. CV columns blank for gamma=off (not CV-scored).</div>
  </div>
</div>

<script>
const DATA = /*__DATA__*/;

// ---- palette: colour by basis family, dash by gamma ---------------------
const BASIS_COLOR = {
  poly2:"#1f77b4", poly3:"#ff7f0e", poly4:"#2ca02c", poly5:"#d62728", poly6:"#9467bd",
  spline3:"#17becf", spline4:"#bcbd22", spline5:"#e377c2", spline6:"#8c564b"
};
const GAMMA_DASH = { off:"dot", scalar:"dash", varying:"solid" };
const REF = DATA.reference;
// canonical form order: poly2..6, spline3..6; within each block off, scalar, varying
const BASIS_RANK = { poly2:0,poly3:1,poly4:2,poly5:3,poly6:4,
                     spline3:5,spline4:6,spline5:7,spline6:8 };
const GAMMA_RANK = { off:0, scalar:1, varying:2 };
function candRank(c){ return (BASIS_RANK[basisOf(c)]??99)*10 + (GAMMA_RANK[gammaOf(c)]??9); }

let state = {
  slug: DATA.slug_order[0],
  nu: null,
  entryAge: DATA.default_entry_age,   // entry age (yr) the curves are drawn at
  visible: {},        // cand -> bool, PERSISTS across slice / nu / age switches
  selForm: null,
  forms: [],          // ordered cand list currently plotted (== trace order)
};
// per-form visibility persists; init = scalar + varying shown, off hidden
(function initVisible(){
  const cands = DATA.slices[DATA.slug_order[0]].cands;
  cands.forEach(c => { state.visible[c] = (c.split("-g")[1] !== "off"); });
})();

const $ = id => document.getElementById(id);

function basisFamily(cand){ return cand.startsWith("poly") ? "poly" : "spline"; }
function gammaOf(cand){ return cand.split("-g")[1]; }
function basisOf(cand){ return cand.split("-g")[0]; }
function fmt(x, nd){ return (x===null||x===undefined) ? "" : Number(x).toFixed(nd); }

// ---- the per-slice/nu set of plotted forms (canonical order) ------------
function visibleForms(){
  const sl = DATA.slices[state.slug];
  return sl.cands.filter(c => (c in (sl.curves[state.nu]||{})));
}

// curve at the chosen entry age: aging + gamma * (age - ae_mean)/ae_off
function curveAt(sl, c){
  const aging = sl.curves[state.nu][c];
  const g = sl.gamma[state.nu][c];
  if (!g) return aging;                       // gamma=off -> entry age irrelevant
  const k = (state.entryAge - sl.ae_mean) / sl.ae_off;
  return aging.map((y,i) => y + g[i]*k);
}
function interpY(A, y, x){                     // linear interp on the (sorted) grid
  if (x<=A[0]) return y[0];
  if (x>=A[A.length-1]) return y[y.length-1];
  for (let i=1;i<A.length;i++){
    if (A[i]>=x){ const t=(x-A[i-1])/(A[i]-A[i-1]); return y[i-1]+t*(y[i]-y[i-1]); }
  }
  return y[y.length-1];
}
// merge spline knots into a polyline, marking only the knot points (size>0)
function withKnots(A, y, knots){
  const pts = A.map((x,i)=>({x, y:y[i], k:false}));
  knots.forEach(kx => pts.push({x:kx, y:interpY(A,y,kx), k:true}));
  pts.sort((a,b)=>a.x-b.x);
  return { x:pts.map(p=>p.x), y:pts.map(p=>p.y), size:pts.map(p=>p.k?7:0) };
}

// ---- chart --------------------------------------------------------------
function buildTraces(){
  const sl = DATA.slices[state.slug];
  const A = sl.A_n;
  const forms = visibleForms();
  state.forms = forms;
  const traces = [];

  // A_n density rug (faint, secondary y-axis so it stays at the bottom)
  const cmax = Math.max(...sl.rug.count, 1);
  const rx = [], ry = [];
  for (let i=0;i<sl.rug.lo.length;i++){
    const mid = (sl.rug.lo[i]+sl.rug.hi[i])/2;
    rx.push(mid); ry.push(sl.rug.count[i]/cmax);
  }
  traces.push({ x:rx, y:ry, type:"bar", name:"A_n density", yaxis:"y2",
    marker:{color:"rgba(150,150,150,0.25)"}, hoverinfo:"skip", showlegend:false });

  forms.forEach(c => {
    const grp = gammaOf(c);
    const col = BASIS_COLOR[basisOf(c)]||"#555";
    const isRef = (c === REF);
    const isSel = (c === state.selForm);
    const y = curveAt(sl, c);
    const knots = (sl.knots && sl.knots[c]) ? sl.knots[c] : null;  // spline gvarying only
    const tr = {
      type:"scatter", name:c,
      line:{ color:col, dash:GAMMA_DASH[grp]||"solid",
             width: isSel ? 5 : (isRef ? 3.4 : 1.7) },
      opacity: isSel ? 1 : 0.9,
      visible: state.visible[c] ? true : "legendonly",
      hovertemplate:"<b>"+c+"</b><br>A_n=%{x:.1f} yr<br>curve=%{y:.4f}<extra></extra>",
    };
    if (knots){
      const w = withKnots(A, y, knots);
      tr.x=w.x; tr.y=w.y; tr.mode="lines+markers";
      tr.marker={ size:w.size, color:col, symbol:"circle", line:{width:0} };
    } else {
      tr.x=A; tr.y=y; tr.mode="lines";
    }
    traces.push(tr);
  });
  return traces;
}

function layout(){
  const sl = DATA.slices[state.slug];
  const shapes = [{ type:"line", x0:0, x1:Math.max(...sl.A_n), y0:0, y1:0,
                    line:{color:"rgba(120,120,120,0.5)", width:1, dash:"dot"} }];
  if (sl.an_p95!==null){
    shapes.push({ type:"line", x0:sl.an_p95, x1:sl.an_p95, yref:"paper", y0:0, y1:1,
      line:{color:"rgba(200,80,80,0.5)", width:1, dash:"dot"} });
  }
  return {
    margin:{l:60,r:14,t:34,b:46},
    title:{ text:state.slug+"  -  aging + entry-age curve   (nu="+state.nu
            +", entry age="+(+state.entryAge).toFixed(0)+" yr)", font:{size:13} },
    xaxis:{ title:"A_n  (years since debut)", zeroline:false,
            gridcolor:"#eef0f3", range:sl.xrange.slice(), autorange:false },
    yaxis:{ title:"aging curve  (log-time)", gridcolor:"#eef0f3",
            range:sl.yrange.slice(), autorange:false },
    yaxis2:{ overlaying:"y", side:"right", range:[0,8], showgrid:false,
             showticklabels:false, fixedrange:true },
    legend:{ font:{size:10}, orientation:"v", x:1.01, y:1, bgcolor:"rgba(255,255,255,0.6)" },
    hovermode:"closest", plot_bgcolor:"#fff", paper_bgcolor:"#fff",
    annotations: sl.an_p95!==null ? [{ x:sl.an_p95, y:1, yref:"paper",
        text:"A_n p95", showarrow:false, font:{size:9,color:"#c05050"},
        xanchor:"left", yanchor:"top" }] : [],
    shapes,
  };
}

let _listenersWired = false;
function drawChart(){
  Plotly.react("chart", buildTraces(), layout(),
    {responsive:true, displayModeBar:true,
     modeBarButtonsToRemove:["lasso2d","select2d"]});
  if (_listenersWired) return;        // Plotly.react reuses the div -> attach once
  _listenersWired = true;
  const el = $("chart");
  el.on("plotly_click", ev => {
    const name = ev.points[0].data.name;
    if (name && name !== "A_n density"){ selectForm(name); }
  });
  el.on("plotly_legendclick", ev => {
    const tr = ev.data[ev.curveNumber];
    const name = tr.name;
    if (!name || name === "A_n density") return false;
    // record the toggle so it persists across slice / nu / age redraws
    state.visible[name] = !(tr.visible === true || tr.visible === undefined);
    syncPills();
    return true;   // let Plotly perform the visibility toggle
  });
}

// ---- table --------------------------------------------------------------
const COLS = [
  {k:"cand",       t:"form",      lft:true},
  {k:"n_params",   t:"k",         nd:0},
  {k:"loglik",     t:"loglik",    nd:1, hi:true},
  {k:"aic",        t:"AIC",       nd:1, lo:true},
  {k:"bic",        t:"BIC",       nd:1, lo:true},
  {k:"cv_per_cell",t:"CV/cell",   nd:5, hi:true},
  {k:"heldout_rmse",t:"CV rmse",  nd:5, lo:true},
];
let sortKey = "cand", sortDesc = false;   // default = canonical legend order

function tableRows(){
  const sl = DATA.slices[state.slug];
  const tbl = sl.table[state.nu] || {};
  return state.forms.filter(c => c in tbl).map(c => ({cand:c, ...tbl[c]}));
}

function bestPerCol(rows){
  const best = {};
  COLS.forEach(col => {
    if (!col.hi && !col.lo) return;
    let bi=-1, bv=null;
    rows.forEach((r,i) => {
      const v=r[col.k]; if (v===null||v===undefined) return;
      if (bv===null || (col.hi ? v>bv : v<bv)){ bv=v; bi=i; }
    });
    best[col.k]=bi;
  });
  return best;
}

function renderTable(){
  const rows = tableRows();
  // sort (cand -> canonical legend order; metrics -> numeric, nulls last)
  rows.sort((a,b)=>{
    if (sortKey==="cand")
      return sortDesc ? candRank(b.cand)-candRank(a.cand) : candRank(a.cand)-candRank(b.cand);
    const x=a[sortKey], y=b[sortKey];
    if (x===null||x===undefined) return 1;
    if (y===null||y===undefined) return -1;
    return sortDesc ? y-x : x-y;
  });
  const best = bestPerCol(rows);
  // CV-best (max cv_per_cell) for the flag
  let cvBest=null, cvBestV=null;
  rows.forEach(r=>{ const v=r.cv_per_cell; if (v!==null&&v!==undefined&&(cvBestV===null||v>cvBestV)){cvBestV=v; cvBest=r.cand;} });

  const thead = "<tr>"+COLS.map(c=>{
    const arrow = (c.k===sortKey)?(sortDesc?" ▼":" ▲"):"";
    return `<th class="${c.lft?'lft':''}" data-k="${c.k}">${c.t}${arrow}</th>`;
  }).join("")+"</tr>";
  $("tbl").querySelector("thead").innerHTML = thead;

  const body = rows.map((r)=>{
    const selCls = (r.cand===state.selForm)?"sel":"";
    const cvCls = (r.cand===cvBest)?"cvbest":"";
    const tds = COLS.map((c,ci)=>{
      if (c.k==="cand"){
        const sw=`<span class="swatch" style="background:${BASIS_COLOR[basisOf(r.cand)]||'#555'}"></span>`;
        return `<td class="lft">${sw}${r.cand}</td>`;
      }
      const v=r[c.k];
      const isBest = (best[c.k]===idxByCandSorted(rows,r.cand));
      const cls = isBest?"best":"";
      return `<td class="${cls}">${v===null||v===undefined?'<span class="muted">-</span>':fmt(v,c.nd)}</td>`;
    }).join("");
    return `<tr class="${selCls} ${cvCls}" data-cand="${r.cand}">${tds}</tr>`;
  }).join("");
  $("tbl").querySelector("tbody").innerHTML = body;

  // wire row clicks
  $("tbl").querySelectorAll("tbody tr").forEach(tr=>{
    tr.onclick = ()=> selectForm(tr.dataset.cand);
  });
  // header sort
  $("tbl").querySelectorAll("thead th").forEach(th=>{
    th.onclick = ()=>{
      const k=th.dataset.k;
      if (k===sortKey) sortDesc=!sortDesc;
      else { sortKey=k; sortDesc = (k!=="cand" && k!=="aic" && k!=="bic" && k!=="heldout_rmse"); }
      renderTable();
    };
  });
  $("tbl-nu").textContent = "(nu="+state.nu+")";
}

// helper: index of a cand within the *sorted* rows array, for best() match
function idxByCandSorted(rows, cand){ return rows.findIndex(r=>r.cand===cand); }

// ---- selection linking --------------------------------------------------
function selectForm(cand, fromLegend){
  state.selForm = (state.selForm===cand && !fromLegend) ? null : cand;
  // restyle line widths/opacity without full redraw
  const upd = {line_width:[], opacity:[]};
  // trace 0 is the rug
  upd.line_width.push(undefined); upd.opacity.push(undefined);
  state.forms.forEach(c=>{
    const isSel=(c===state.selForm), isRef=(c===REF);
    upd.line_width.push(isSel?5:(isRef?3.4:1.7));
    upd.opacity.push(isSel?1:0.9);
  });
  Plotly.restyle("chart", {"line.width":upd.line_width, "opacity":upd.opacity});
  renderTable();
}

// ---- visibility (persists in state.visible) -----------------------------
function applyVisibility(){
  const vis = state.forms.map(c => state.visible[c] ? true : "legendonly");
  const idx = vis.map((_,i)=>i+1);          // +1 for the rug trace
  if (idx.length) Plotly.restyle("chart", {visible: vis}, idx);
  syncPills();
}
// reflect "is every form in this group shown?" on the quick-filter pills
function syncPills(){
  $("gamma-row").querySelectorAll(".toggle-btn").forEach(b=>{
    const f = state.forms.filter(c=>gammaOf(c)===b.dataset.g);
    b.classList.toggle("active", f.length>0 && f.every(c=>state.visible[c]));
  });
  $("basis-row").querySelectorAll(".toggle-btn").forEach(b=>{
    const f = state.forms.filter(c=>basisFamily(c)===b.dataset.b);
    b.classList.toggle("active", f.length>0 && f.every(c=>state.visible[c]));
  });
}

// ---- header wiring ------------------------------------------------------
function buildSliceInfo(){
  const sl = DATA.slices[state.slug];
  $("slice-info").innerHTML =
    `I=<b>${sl.I.toLocaleString()}</b> &nbsp; J=<b>${sl.J.toLocaleString()}</b> `+
    `&nbsp; N=<b>${sl.N.toLocaleString()}</b> &nbsp; A_n p95=<b>${sl.an_p95}</b> yr`;
}

function buildNuRow(){
  const sl = DATA.slices[state.slug];
  if (!sl.nus.includes(state.nu)) state.nu = sl.nus.includes("8") ? "8" : sl.nus[0];
  $("nu-row").innerHTML = sl.nus.map(n =>
    `<span class="toggle-btn ${n===state.nu?'active':''}" data-nu="${n}">${n}</span>`).join("");
  $("nu-row").querySelectorAll(".toggle-btn").forEach(b=>{
    b.onclick=()=>{ state.nu=b.dataset.nu; state.selForm=null; refresh(); };
  });
}

function buildAgeRow(){
  $("age-row").innerHTML = DATA.entry_ages.map(a =>
    `<span class="toggle-btn ${a===state.entryAge?'active':''}" data-age="${a}">${(+a).toFixed(0)}</span>`).join("");
  $("age-row").querySelectorAll(".toggle-btn").forEach(b=>{
    b.onclick=()=>{ state.entryAge=parseFloat(b.dataset.age);
      $("age-row").querySelectorAll(".toggle-btn").forEach(x=>
        x.classList.toggle("active", parseFloat(x.dataset.age)===state.entryAge));
      drawChart(); applyVisibility(); };   // curves move with entry age; visibility persists
  });
}

function refresh(){
  buildNuRow();
  buildAgeRow();
  buildSliceInfo();
  drawChart();
  applyVisibility();
  renderTable();
}

function init(){
  // slice dropdown
  $("sel-slice").innerHTML = DATA.slug_order.map(s=>`<option value="${s}">${s}</option>`).join("");
  $("sel-slice").value = state.slug;
  $("sel-slice").onchange = e=>{ state.slug=e.target.value; state.selForm=null; refresh(); };

  // gamma / basis quick filters -> toggle whole group's persistent visibility
  $("gamma-row").querySelectorAll(".toggle-btn").forEach(b=>{
    b.onclick=()=>{ const pred=c=>gammaOf(c)===b.dataset.g;
      const target=!state.forms.filter(pred).every(c=>state.visible[c]);
      state.forms.filter(pred).forEach(c=>state.visible[c]=target);
      applyVisibility(); };
  });
  $("basis-row").querySelectorAll(".toggle-btn").forEach(b=>{
    b.onclick=()=>{ const pred=c=>basisFamily(c)===b.dataset.b;
      const target=!state.forms.filter(pred).every(c=>state.visible[c]);
      state.forms.filter(pred).forEach(c=>state.visible[c]=target);
      applyVisibility(); };
  });

  // select-all / unselect-all
  const setAll = on => {
    DATA.slices[state.slug].cands.forEach(c=>state.visible[c]=on);
    applyVisibility();
  };
  $("btn-all").onclick  = ()=>setAll(true);
  $("btn-none").onclick = ()=>setAll(false);

  refresh();
}
init();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
