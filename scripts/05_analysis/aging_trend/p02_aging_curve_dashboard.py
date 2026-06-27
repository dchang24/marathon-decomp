"""Interactive aging-curve dashboard (self-contained HTML).

Same content as the static `p01_aging_curve.py` figure -- left men | right women,
the aging curve f(A_n) fanned by entry age (35/45/55/65) with bootstrap bands --
but as one HTML page with three independent dropdowns:

  * **slice**  ALL | Po10
  * **drift (d_i)**  on (`full` = AxD) | off (`agingS4gv` = aging-only)
  * **min race count**  mrc2 | mrc5

Every (slice x drift x mrc) combination is precomputed in Python (curve math
reused from `p01_aging_curve.reconstruct`), serialized to JSON, and embedded in a
standalone page that redraws via Plotly.js (CDN) on dropdown change. Where a
bootstrap exists the 95% band is drawn; otherwise the curve is point-only (line,
no band) and the panel is annotated as such -- bootstrap currently covers all
mrc5 fits plus the ALL mrc2 aging-only fits.

No refit; reads `results/models/{slug}/{model}_nu8p00_best__*/` (+ its
`bootstrap/global_coeffs.parquet`).

Output:
    results/analysis/aging_trend/aging_curve_dashboard.html

Run::

    python scripts/05_analysis/aging_trend/p02_aging_curve_dashboard.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir (p01)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from p01_aging_curve import ENTRY_AGE_COLORS, OUT_ROOT, reconstruct  # noqa: E402

COHORTS = ["ALL", "Po10"]
# label -> registry model tag.  D on = AxD ("full"); D off = aging-only ("agingS4gv").
MODELS = {"on": "full", "off": "agingS4gv"}
MRCS = {"2": "mrc2", "5": "mrc5"}
DEFAULT_ENTRY_AGES = (35, 45, 55, 65)


def _hex_to_rgba(h: str, a: float) -> str:
    h = h.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"


def build_payload(*, nutag: str, data_version: str, entry_ages, scale: float) -> dict:
    """Nested dict: data[combo_key][sex] = {A_grid, has_boot, curves, rug}."""
    combos: dict = {}
    for cohort in COHORTS:
        for dlabel, model in MODELS.items():
            for mlabel, mrc in MRCS.items():
                combo_key = f"{cohort}|{dlabel}|{mlabel}"
                per_sex: dict = {}
                for sx in ("M", "W"):
                    slug = f"{cohort}_{sx}_14-25_{mrc}"
                    print(f"reconstructing {slug}  (d={dlabel}, {mrc}) ...")
                    fan = reconstruct(slug, model=model, nutag=nutag,
                                      data_version=data_version, allow_point_only=True)
                    if fan is None:
                        continue
                    curves = {}
                    for ea in entry_ages:
                        p, lo, hi = fan.at_entry(ea)
                        curves[str(int(ea))] = {
                            "point": (p * scale).tolist(),
                            "lo": None if lo is None else (lo * scale).tolist(),
                            "hi": None if hi is None else (hi * scale).tolist(),
                        }
                    centers, counts = fan.an_hist
                    per_sex[sx] = {
                        "A_grid": fan.A_grid.tolist(),
                        "has_boot": bool(fan.has_boot),
                        "mean_Ae": fan.mean_Ae,
                        "curves": curves,
                        "rug": {"centers": centers.tolist(),
                                "counts": counts.astype(float).tolist()},
                    }
                if per_sex:
                    combos[combo_key] = per_sex

    return {
        "combos": combos,
        "entry_ages": [int(e) for e in entry_ages],
        "colors": {str(int(e)): ENTRY_AGE_COLORS.get(int(e), "#444") for e in entry_ages},
        "colors_band": {str(int(e)): _hex_to_rgba(ENTRY_AGE_COLORS.get(int(e), "#444"), 0.18)
                        for e in entry_ages},
        "cohorts": COHORTS,
        "models": list(MODELS.keys()),
        "mrcs": list(MRCS.keys()),
    }


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Marathon aging curve by entry age</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {{ font-family: system-ui, Arial, sans-serif; margin: 18px; color: #222;
          background: #fff; }}
  #controls {{ margin-bottom: 10px; }}
  #controls label {{ margin-right: 18px; font-size: 14px; }}
  select {{ font-size: 14px; padding: 2px 4px; }}
  #plots {{ display: flex; gap: 8px; }}
  #plotM, #plotW {{ flex: 1 1 0; min-width: 0; height: 540px; }}
  h2 {{ font-size: 18px; margin: 0 0 4px 0; }}
  .note {{ color: #666; font-size: 12px; margin-top: 4px; }}
</style>
</head>
<body>
<h2>Marathon aging curve by entry age</h2>
<div id="controls">
  <label>slice
    <select id="sel-cohort"></select></label>
  <label>drift d<sub>i</sub>
    <select id="sel-model">
      <option value="on">on (AxD)</option>
      <option value="off">off (aging only)</option>
    </select></label>
  <label>min race count
    <select id="sel-mrc">
      <option value="5">mrc5</option>
      <option value="2">mrc2</option>
    </select></label>
</div>
<div id="plots">
  <div id="plotM"></div>
  <div id="plotW"></div>
</div>
<div class="note">x: career age A<sub>n</sub> (years since debut). y: aging block in
log-time x100 (~ % slowdown vs the athlete's own debut). beta=0 APC gauge; band = 95%
bootstrap where available. nu=8.</div>

<script>
const DATA = {payload_json};
const ENTRY = DATA.entry_ages;
const SEX_LABEL = {{M: "men", W: "women"}};

function buildTraces(panel) {{
  // panel: {{A_grid, has_boot, curves, rug}}  or undefined
  const traces = [];
  if (!panel) return traces;
  const x = panel.A_grid;
  for (const ea of ENTRY) {{
    const key = String(ea);
    const c = panel.curves[key];
    if (!c) continue;
    if (c.lo && c.hi) {{
      traces.push({{
        x: x.concat(x.slice().reverse()),
        y: c.hi.concat(c.lo.slice().reverse()),
        fill: "toself", fillcolor: DATA.colors_band[key],
        line: {{width: 0}}, hoverinfo: "skip", showlegend: false,
        type: "scatter", mode: "lines",
      }});
    }}
    traces.push({{
      x: x, y: c.point, type: "scatter", mode: "lines",
      name: "entry age " + ea,
      line: {{color: DATA.colors[key], width: 2.4}},
      hovertemplate: "A_n=%{{x:.1f}}<br>%{{y:.2f}}<extra>entry " + ea + "</extra>",
    }});
  }}
  return traces;
}}

function yRange(panels) {{
  let lo = Infinity, hi = -Infinity;
  for (const p of panels) {{
    if (!p) continue;
    for (const ea of ENTRY) {{
      const c = p.curves[String(ea)];
      if (!c) continue;
      for (const arr of [c.point, c.lo, c.hi]) {{
        if (!arr) continue;
        for (const v of arr) {{ if (v < lo) lo = v; if (v > hi) hi = v; }}
      }}
    }}
  }}
  if (!isFinite(lo)) {{ lo = -1; hi = 1; }}
  const pad = 0.08 * (hi - lo + 1e-9);
  return [lo - pad, hi + pad];
}}

// Fixed axis ranges across every slice/model/mrc combo so the axes do not jump
// when toggling. Includes both sexes (and all bands, for y).
const ALL_PANELS = (() => {{
  const all = [];
  for (const k of Object.keys(DATA.combos))
    for (const sx of ["M", "W"])
      if (DATA.combos[k][sx]) all.push(DATA.combos[k][sx]);
  return all;
}})();
const Y_RANGE = yRange(ALL_PANELS);
const X_RANGE = (() => {{
  let hi = 0;
  for (const p of ALL_PANELS) {{
    const g = p.A_grid;
    if (g.length) hi = Math.max(hi, g[g.length - 1]);
  }}
  return [0, hi];
}})();

function rugTrace(panel, yr) {{
  // light grey A_n density rug along the bottom
  if (!panel) return null;
  const ctr = panel.rug.centers, cnt = panel.rug.counts;
  const cmax = Math.max(...cnt) || 1;
  const h = 0.06 * (yr[1] - yr[0]);
  return {{
    x: ctr, y: cnt.map(v => h * v / cmax), base: yr[0],
    type: "bar", marker: {{color: "rgba(0,0,0,0.10)"}},
    width: (ctr.length > 1 ? (ctr[1] - ctr[0]) : 1),
    hoverinfo: "skip", showlegend: false,
  }};
}}

function render() {{
  const cohort = document.getElementById("sel-cohort").value;
  const model = document.getElementById("sel-model").value;
  const mrc = document.getElementById("sel-mrc").value;
  const key = [cohort, model, mrc].join("|");
  const combo = DATA.combos[key] || {{}};
  const pM = combo.M, pW = combo.W;
  const yr = Y_RANGE;

  for (const [div, sx, panel] of [["plotM", "M", pM], ["plotW", "W", pW]]) {{
    const traces = buildTraces(panel);
    const rug = rugTrace(panel, yr);
    if (rug) traces.unshift(rug);
    let title = cohort + " — " + SEX_LABEL[sx];
    if (panel && !panel.has_boot) title += "  (point only, no bootstrap)";
    if (!panel) title += "  (no fit)";
    const layout = {{
      title: {{text: title, font: {{size: 15}}}},
      margin: {{l: 60, r: 16, t: 40, b: 48}},
      xaxis: {{title: "career age A_n (years since debut)", zeroline: false,
               range: X_RANGE, showgrid: true, gridcolor: "#e6e6e6"}},
      yaxis: {{title: (sx === "M" ? "aging curve f(A_n)  (x100 ~ % slowdown)" : ""),
               range: yr, zeroline: true, zerolinecolor: "#bbb",
               showgrid: true, gridcolor: "#e6e6e6"}},
      showlegend: (sx === "M"),
      legend: {{x: 0.02, y: 0.98, bgcolor: "rgba(255,255,255,0.6)"}},
      shapes: [{{type: "line", x0: 0, x1: 1, xref: "paper",
                 y0: 0, y1: 0, line: {{color: "#999", width: 1, dash: "dash"}}}}],
    }};
    Plotly.react(div, traces, layout, {{responsive: true, displaylogo: false}});
  }}
}}

// populate slice dropdown then wire events
const selC = document.getElementById("sel-cohort");
for (const c of DATA.cohorts) {{
  const o = document.createElement("option"); o.value = c; o.text = c; selC.add(o);
}}
for (const id of ["sel-cohort", "sel-model", "sel-mrc"])
  document.getElementById(id).addEventListener("change", render);
render();
</script>
</body>
</html>
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--entry-ages", type=float, nargs="+", default=list(DEFAULT_ENTRY_AGES))
    ap.add_argument("--scale", type=float, default=100.0)
    args = ap.parse_args()

    payload = build_payload(nutag=args.nutag, data_version=args.data_version,
                            entry_ages=args.entry_ages, scale=args.scale)
    if not payload["combos"]:
        print("no fits found; nothing written.")
        return

    html = HTML_TEMPLATE.format(payload_json=json.dumps(payload))
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "aging_curve_dashboard.html"
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}  ({len(payload['combos'])} slice/model/mrc combos)")


if __name__ == "__main__":
    main()
