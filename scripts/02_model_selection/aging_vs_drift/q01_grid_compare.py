"""2x2 model-term grid: does aging absorb d_i, does d_i absorb aging, need both?

Compares the four registered nu=8 fits per (slice, mrc) -- the {no-aging, +aging}
x {no-d, +d} grid -- WITHOUT any refit and WITHOUT touching `delta` (the field
career-stage composition; that is the later q0x bias test). Strictly a read of
the already-registered fits:

    no-d, no-aging  = baseline_nu8p00_best
    no-d, +aging    = agingS4gv_nu8p00_best
    +d,  no-aging   = drift_nu8p00_best
    +d,  +aging     = full_nu8p00_best

Runs at every requested ``--mrc`` (default both 2 and 5), so the same grid is
read at the everyone (mrc2) and dedicated-runner (mrc5) field cuts.

Structural prior (why overlap is expected in only ONE direction): aging and d_i
share the regressor A_n (within-career age). aging is a *population, non-linear*
function of A_n; d_i is a *per-athlete, linear, within-athlete-centered* slope
(EB-shrunk to mean zero). They can overlap only in the linear-population-slope
direction -- which is APC-unidentified and gauged away. So we expect: aging
*curvature* untouched by adding d_i, and d_i *heterogeneity* untouched by adding
aging. Each section below tests one face of that.

Per-(slice, mrc) sections (printed AND mirrored to markdown):
  1. FIT 2x2      -- loglik / RSS / AIC / BIC per cell; marginal gains; the
                     non-additivity (interaction) term. Descriptive only.
  2. CONTRIBUTIONS -- full fit, per-obs aging contribution vs d_i contribution:
                     correlation + variance split. ~0 corr => distinct signal.
  3. DRIFT SHIFT   -- drift vs full: does aging ON shrink the individual drift?
  4. AGING SHIFT   -- aging vs full: does d ON change the population aging curve?
  5. STABILITY     -- 4x4 v_j / u_i agreement (Pearson, Spearman, max|delta|).

Per-slice CROSS-MRC section (when >=2 mrc available):
  6. CURVE DE-BIASING -- the deferred check: overlay the population aging curve at
                     {no-d, +d} x {mrc2, mrc5} on one gauge (centered to mean 0 on
                     the shared A_n range). Hypothesis: the no-d / everyone (mrc2)
                     curve is the contaminated outlier; adding d_i OR restricting
                     to dedicated runners (mrc5) both move it onto the others.
                     Reports peak ages + the gap-closing contrasts:
                       gap_noD   = ||agingMrcLo - agingMrcHi||   (everyone vs dedicated, no d)
                       gap_withD = ||fullMrcLo  - fullMrcHi||    (everyone vs dedicated, +d)
                       d_on_mrcLo / d_on_mrcHi = ||aging - full|| within each mrc
                     Expect gap_withD << gap_noD and d_on_mrcLo >> d_on_mrcHi.

Outputs (under results/model_selection/aging_vs_drift/):
  * grid_compare.md   -- the full console run, formatted
  * fit_grid.csv      -- per (slug, mrc, cell): loglik, rss, aic, bic, n_params
  * gains.csv         -- per (slug, mrc): marginal gains + interaction
  * block_shift.csv   -- per (slug, mrc): contribution + drift + aging shift scalars
  * stability.csv     -- per (slug, mrc, pair): pearson, spearman, max_abs
  * curve_debias.csv  -- per slice: cross-mrc peak ages + gap-closing contrasts

Argument-free (VS Code play): every slice, both mrc 2 and 5, that has the grid.

Run::

    python scripts/02_model_selection/aging_vs_drift/q01_grid_compare.py
    python scripts/02_model_selection/aging_vs_drift/q01_grid_compare.py --slices Po10_M
    python scripts/02_model_selection/aging_vs_drift/q01_grid_compare.py --mrc 5 --slices Po10_M ALL_W
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import ModelConfig, load_slice, registry   # noqa: E402
from marathon_decomp.config import RESULTS_DIR, display_path     # noqa: E402

from baseline_common import slices as S                          # noqa: E402
from aging_common.fitting import aging_cfg, aging_stem           # noqa: E402
from drift_common.fitting import drift_cfg, drift_stem           # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging_vs_drift"
MODELS_ROOT = RESULTS_DIR / "models"
CONV_ROOT = RESULTS_DIR / "convergence"

# grid cell -> (aging on?, d on?) and the human label of each axis position
CELLS = ["baseline", "aging", "drift", "full"]
CELL_LABEL = {"baseline": "no-d / no-aging", "aging": "no-d / +aging",
              "drift": "+d / no-aging", "full": "+d / +aging"}
AE_OFFSET_YR = 10.0
N_CURVE_GRID = 100   # shared A_n grid for the cross-mrc curve comparison


# ── locating the four registered fits ────────────────────────────────
def grid_cfg_stem(nu: float) -> dict[str, tuple[ModelConfig, str]]:
    base = ModelConfig(use_phi12=False, use_gamma=False, use_d=False, nu=float(nu))
    ag = aging_cfg(nu)
    dr = drift_cfg(nu, variant="drift")
    fu = drift_cfg(nu, variant="full")
    return {
        "baseline": (base, registry.model_stem(base, "best")),
        "aging": (ag, aging_stem(ag)),
        "drift": (dr, drift_stem("drift", nu)),
        "full": (fu, drift_stem("full", nu)),
    }


def fit_dirs(spec, nu: float) -> dict[str, Path]:
    parent = MODELS_ROOT / registry.slice_slug(spec)
    return {c: registry.fit_path(parent, stem, spec, cfg, resample_tag="base")
            for c, (cfg, stem) in grid_cfg_stem(nu).items()}


def all_present(dirs: dict[str, Path]) -> bool:
    return all((d / "fit.pkl").is_file() and (d / "manifest.json").is_file()
               for d in dirs.values())


def _manifest(d: Path) -> dict:
    return json.loads((d / "manifest.json").read_text())


# ── small stats helpers (no scipy) ───────────────────────────────────
def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 3 or a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _rank(a: np.ndarray) -> np.ndarray:
    return a.argsort().argsort().astype(float)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    return _pearson(_rank(np.asarray(a, float)), _rank(np.asarray(b, float)))


def _rms(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a, float) - np.asarray(b, float)) ** 2)))


# ── markdown / console reporting ─────────────────────────────────────
class Report:
    """Echoes everything to the console and accumulates a markdown mirror."""

    def __init__(self) -> None:
        self.md: list[str] = []

    def line(self, text: str = "") -> None:
        print(text)
        self.md.append(text)

    def head(self, text: str, level: int = 2) -> None:
        print(f"\n{'=' * 4} {text} {'=' * 4}" if level == 2 else f"\n-- {text} --")
        self.md.append("")
        self.md.append(f"{'#' * level} {text}")
        self.md.append("")

    def table(self, df: pd.DataFrame, *, index: bool = False, floatfmt: str = "{:.6g}") -> None:
        def fmt(x):
            if isinstance(x, float):
                return "-" if not np.isfinite(x) else floatfmt.format(x)
            return str(x)
        disp = df.copy()
        for c in disp.columns:
            disp[c] = disp[c].map(fmt)
        with pd.option_context("display.width", 200, "display.max_columns", 40):
            print(disp.to_string(index=index))
        self.md.append(_df_to_md(disp, index=index))
        self.md.append("")

    def save(self, path: Path, title: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# " + title + "\n\n" + "\n".join(self.md) + "\n", encoding="utf-8")
        print(f"\nmarkdown -> {path}")


def _md_cell(x) -> str:
    return str(x).replace("|", "\\|")   # escape pipes so cells don't split columns


def _df_to_md(df: pd.DataFrame, *, index: bool) -> str:
    cols = ([df.index.name or ""] if index else []) + [str(c) for c in df.columns]
    rows = []
    for idx, row in df.iterrows():
        cells = ([_md_cell(idx)] if index else []) + [_md_cell(v) for v in row]
        rows.append("| " + " | ".join(cells) + " |")
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    return "\n".join([header, sep] + rows)


def _grid2x2(values: dict[str, float], name: str) -> pd.DataFrame:
    """A 2x2 (rows no-d/+d, cols no-aging/+aging) frame for a per-cell scalar."""
    df = pd.DataFrame(
        [[values["baseline"], values["aging"]],
         [values["drift"], values["full"]]],
        index=["no-d", "+d"], columns=["no-aging", "+aging"])
    df.index.name = name
    return df


# ── per-(slice, mrc) analysis ────────────────────────────────────────
def run_cell_grid(name: str, mrc: int, nu: float, rep: Report) -> dict | None:
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    dirs = fit_dirs(spec, nu)
    if not all_present(dirs):
        missing = [c for c in CELLS if not (dirs[c] / "fit.pkl").is_file()]
        rep.line(f"[skip] {slug}: missing fits {missing}")
        return None

    rep.head(f"{slug}   nu={nu:g}", level=2)
    fd = load_slice(spec)
    rep.line(f"I={fd.I:,}  J={fd.J:,}  N={fd.N:,}   "
             f"cells: baseline / aging / drift / full")

    man = {c: _manifest(dirs[c]) for c in CELLS}
    models = {c: registry.load_fit(dirs[c], fd) for c in CELLS}
    tag = dict(slug=slug, mrc=mrc, nu=nu)
    csv: dict = {"fit_grid": [], "gains": {}, "block_shift": {}, "stability": []}

    # ---------------- 1. FIT 2x2 -------------------------------------
    rep.head("1. Fit 2x2 (descriptive: loglik / RSS / AIC / BIC)", level=3)
    metric_keys = dict(loglik="loglik_final", rss="rss_final", aic="aic",
                       bic="bic", n_params="n_params_naive")
    vals = {m: {c: float(man[c].get(k, np.nan)) for c in CELLS}
            for m, k in metric_keys.items()}
    fit_tbl = pd.DataFrame({"cell": CELLS,
                            "label": [CELL_LABEL[c] for c in CELLS],
                            **{m: [vals[m][c] for c in CELLS] for m in metric_keys}})
    rep.table(fit_tbl)
    for c in CELLS:
        csv["fit_grid"].append(dict(**tag, cell=c, **{m: vals[m][c] for m in metric_keys}))

    rep.line("\nloglik on the 2x2 grid:")
    rep.table(_grid2x2(vals["loglik"], "loglik"), index=True, floatfmt="{:.4f}")

    ll = vals["loglik"]
    gains = {
        "aging_gain | no-d": ll["aging"] - ll["baseline"],
        "aging_gain | +d": ll["full"] - ll["drift"],
        "d_gain | no-aging": ll["drift"] - ll["baseline"],
        "d_gain | +aging": ll["full"] - ll["aging"],
        "interaction (non-additivity)": ll["full"] + ll["baseline"] - ll["aging"] - ll["drift"],
    }
    rep.line("\nmarginal loglik gains (does a term survive adding the other?):")
    rep.table(pd.DataFrame({"effect": list(gains), "d_loglik": list(gains.values())}),
              floatfmt="{:.4f}")
    csv["gains"] = dict(**tag, **gains)

    # ---------------- 2. CONTRIBUTIONS (full fit) --------------------
    rep.head("2. Contribution overlap in the full fit (aging vs d_i)", level=3)
    terms = models["full"]._yhat_terms()
    aging_c = np.zeros(fd.N)
    if "aging" in terms:
        aging_c = aging_c + terms["aging"]
    if "gamma" in terms:
        aging_c = aging_c + terms["gamma"]
    drift_c = terms["d"]
    va, vd = float(np.var(aging_c)), float(np.var(drift_c))
    cov = float(np.cov(aging_c, drift_c)[0, 1])
    corr_c = _pearson(aging_c, drift_c)
    contrib = dict(corr_aging_drift=corr_c, var_aging=va, var_drift=vd,
                   cov=cov, sd_ratio_d_over_aging=float(np.sqrt(vd / va)) if va > 0 else np.nan)
    rep.table(pd.DataFrame({"quantity": list(contrib), "value": list(contrib.values())}),
              floatfmt="{:.6g}")
    rep.line("\n(corr near 0 => terms carry distinct signal => both needed; "
             "gauge-sensitive in the linear direction.)")

    # ---------------- 3. DRIFT SHIFT (drift vs full) -----------------
    rep.head("3. Drift block: does turning aging ON shrink d_i?", level=3)
    elig = np.asarray(models["full"].eligible_d, bool) & np.asarray(models["drift"].eligible_d, bool)
    d_dr = np.asarray(models["drift"].params["d"], float)[elig]
    d_fu = np.asarray(models["full"].params["d"], float)[elig]
    drift_shift = dict(
        omega_d2_drift=float(man["drift"].get("omega_d2", np.nan)),
        omega_d2_full=float(man["full"].get("omega_d2", np.nan)),
        n_eligible=int(elig.sum()),
        d_sd_drift=float(d_dr.std(ddof=1)) if elig.sum() > 1 else np.nan,
        d_sd_full=float(d_fu.std(ddof=1)) if elig.sum() > 1 else np.nan,
        frac_improver_drift=float((d_dr < 0).mean()) if elig.sum() else np.nan,
        frac_improver_full=float((d_fu < 0).mean()) if elig.sum() else np.nan,
        corr_d_drift_full=_pearson(d_dr, d_fu),
        max_abs_dd=float(np.max(np.abs(d_dr - d_fu))) if elig.sum() else np.nan,
    )
    rep.table(pd.DataFrame({"quantity": list(drift_shift), "value": list(drift_shift.values())}),
              floatfmt="{:.6g}")

    # ---------------- 4. AGING SHIFT (aging vs full) -----------------
    rep.head("4. Aging block: does turning d ON change the aging curve / gamma?", level=3)
    ca = _load_curve(slug, "aging_S4gv")
    cf = _load_curve(slug, "drift_full")
    aging_shift: dict = {}
    if ca is not None and cf is not None and np.allclose(ca["A_n"], cf["A_n"]):
        ya, yf = ca["aging_curve"].to_numpy(), cf["aging_curve"].to_numpy()
        An = ca["A_n"].to_numpy()
        aging_shift.update(
            curve_corr=_pearson(ya, yf),
            curve_max_abs_delta=float(np.max(np.abs(ya - yf))),
            peak_age_aging=float(An[int(np.argmin(ya))]),
            peak_age_full=float(An[int(np.argmin(yf))]),
        )
    else:
        rep.line("(curves.parquet missing or grid mismatch -- curve compare skipped)")
    g_ag = np.atleast_1d(np.asarray(models["aging"].params["gamma"], float))
    g_fu = np.atleast_1d(np.asarray(models["full"].params["gamma"], float))
    if g_ag.size == g_fu.size and g_ag.size:
        aging_shift["gamma_max_abs_delta"] = float(np.max(np.abs(g_ag - g_fu)))
        aging_shift["gamma_corr"] = _pearson(g_ag, g_fu) if g_ag.size >= 3 else np.nan
    if aging_shift:
        rep.table(pd.DataFrame({"quantity": list(aging_shift), "value": list(aging_shift.values())}),
                  floatfmt="{:.6g}")

    csv["block_shift"] = dict(**tag, **contrib, **drift_shift, **aging_shift)

    # ---------------- 5. STABILITY (v, u) ----------------------------
    rep.head("5. Factor stability across the grid (v_j and u_i)", level=3)
    for param in ("v", "u"):
        vecs = {c: np.asarray(models[c].params[param], float) for c in CELLS}
        vecs = {c: x - x.mean() for c, x in vecs.items()}   # mean-center (gauge)
        rep.line(f"\n{param}_j Pearson corr (4x4):")
        rep.table(_pair_matrix(vecs, _pearson), index=True, floatfmt="{:.4f}")
        rep.line(f"{param}_j Spearman corr (4x4):")
        rep.table(_pair_matrix(vecs, _spearman), index=True, floatfmt="{:.4f}")
        rep.line(f"{param}_j max|delta| (4x4):")
        rep.table(_pair_matrix(vecs, lambda a, b: float(np.max(np.abs(a - b)))),
                  index=True, floatfmt="{:.4g}")
        for a in CELLS:
            for b in CELLS:
                if a < b:
                    csv["stability"].append(dict(
                        **tag, param=param, cell_a=a, cell_b=b,
                        pearson=_pearson(vecs[a], vecs[b]),
                        spearman=_spearman(vecs[a], vecs[b]),
                        max_abs=float(np.max(np.abs(vecs[a] - vecs[b])))))

    return csv


# ── per-slice cross-mrc aging-curve de-biasing ───────────────────────
def _centered_curves_on_common_grid(name: str, mrcs: list[int], nu: float):
    """Load the population aging curve at (mrc, variant) for variant in {noD,withD}.

    Returns (grid, curves) where curves maps (mrc, variant) -> centered y on a
    shared A_n grid (mean-removed, so only shape/curvature is compared under the
    common APC beta=0 gauge). variant 'noD' = agingS4gv curve, 'withD' = full curve.
    """
    raw: dict[tuple[int, str], pd.DataFrame] = {}
    for m in mrcs:
        slug = registry.slice_slug(S.build_spec(name, min_race_count=m))
        for variant, stage in (("noD", "aging_S4gv"), ("withD", "drift_full")):
            c = _load_curve(slug, stage)
            if c is not None:
                raw[(m, variant)] = c
    if len({m for (m, _) in raw}) < 2:
        return None, None
    a_hi = min(float(c["A_n"].max()) for c in raw.values())
    grid = np.linspace(0.0, a_hi, N_CURVE_GRID)
    curves = {}
    for k, c in raw.items():
        y = np.interp(grid, c["A_n"].to_numpy(), c["aging_curve"].to_numpy())
        curves[k] = y - y.mean()
    return grid, curves


def cross_mrc_curve_debias(name: str, mrcs: list[int], nu: float, rep: Report) -> dict | None:
    grid, curves = _centered_curves_on_common_grid(name, mrcs, nu)
    if curves is None:
        return None
    have_mrc = sorted({m for (m, _) in curves})
    lo, hi = have_mrc[0], have_mrc[-1]

    rep.head(f"6. Aging-curve de-biasing across mrc -- {name} (mrc {lo} vs {hi})", level=3)
    rep.line(f"population aging curve, centered to mean 0 on shared A_n in "
             f"[0, {grid[-1]:.2f}] yr ({N_CURVE_GRID} pts), one APC beta=0 gauge.")

    # peak ages
    peak_rows = []
    for (m, variant), y in sorted(curves.items()):
        peak_rows.append(dict(mrc=m, variant=variant,
                              peak_age=float(grid[int(np.argmin(y))])))
    rep.line("\npeak age (A_n at fastest, yr) per curve:")
    rep.table(pd.DataFrame(peak_rows), floatfmt="{:.3f}")

    row = dict(slug=registry.slice_slug(S.build_spec(name, min_race_count=lo)),
               slice=name, nu=nu, mrc_lo=lo, mrc_hi=hi)
    for (m, variant), y in curves.items():
        row[f"peak_{variant}_mrc{m}"] = float(grid[int(np.argmin(y))])

    # gap-closing contrasts (only the pairs that exist)
    def has(*keys):
        return all(k in curves for k in keys)

    contrasts: dict[str, float] = {}
    if has((lo, "noD"), (hi, "noD")):
        contrasts["gap_noD (everyone vs dedicated, no d)"] = _rms(curves[(lo, "noD")], curves[(hi, "noD")])
    if has((lo, "withD"), (hi, "withD")):
        contrasts["gap_withD (everyone vs dedicated, +d)"] = _rms(curves[(lo, "withD")], curves[(hi, "withD")])
    if has((lo, "noD"), (lo, "withD")):
        contrasts[f"d_on_mrc{lo} (||aging-full|| @ everyone)"] = _rms(curves[(lo, "noD")], curves[(lo, "withD")])
    if has((hi, "noD"), (hi, "withD")):
        contrasts[f"d_on_mrc{hi} (||aging-full|| @ dedicated)"] = _rms(curves[(hi, "noD")], curves[(hi, "withD")])
    if contrasts:
        rep.line("\ncurve-distance contrasts (RMS on the shared grid):")
        rep.table(pd.DataFrame({"contrast": list(contrasts), "rms": list(contrasts.values())}),
                  floatfmt="{:.5g}")
        rep.line("\n(expect gap_withD << gap_noD and d_on_mrc"
                 f"{lo} >> d_on_mrc{hi} if no-d/everyone is the contaminated curve.)")
    for k, v in contrasts.items():
        row[k.split(" ")[0]] = v
    return row


def _load_curve(slug: str, stage: str) -> pd.DataFrame | None:
    p = CONV_ROOT / slug / stage / "curves.parquet"
    return pd.read_parquet(p) if p.is_file() else None


def _pair_matrix(vecs: dict[str, np.ndarray], fn) -> pd.DataFrame:
    M = pd.DataFrame(index=CELLS, columns=CELLS, dtype=float)
    for a in CELLS:
        for b in CELLS:
            M.loc[a, b] = fn(vecs[a], vecs[b])
    M.index.name = "cell"
    return M


# ── driver ───────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=["all"],
                    help="slice name(s) (e.g. Po10_M) or 'all' (default).")
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5],
                    help="min-race-count field cut(s) to read (default: 2 5).")
    ap.add_argument("--nu", type=float, default=8.0, help="Student-t nu (default 8).")
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    rep = Report()
    rep.line(f"2x2 aging-vs-drift term grid   nu={args.nu:g}   mrc={mrcs}   "
             f"slices requested: {', '.join(names)}")

    fit_grid, gains, block_shift, stability, curve_rows = [], [], [], [], []
    n_cells = 0
    for name in names:
        for mrc in mrcs:
            csv = run_cell_grid(name, mrc, args.nu, rep)
            if csv is None:
                continue
            n_cells += 1
            fit_grid += csv["fit_grid"]
            gains.append(csv["gains"])
            block_shift.append(csv["block_shift"])
            stability += csv["stability"]
        if len(mrcs) >= 2:
            cr = cross_mrc_curve_debias(name, mrcs, args.nu, rep)
            if cr is not None:
                curve_rows.append(cr)

    if n_cells == 0:
        rep.line("\nNo (slice, mrc) had all four fits. Run e01_fit_drift.py "
                 "(--variant both) for the target slices/mrc first.")
        rep.save(OUT_ROOT / "grid_compare.md", "Aging vs Drift -- 2x2 term grid")
        return

    # ---- tidy CSV rollups ------------------------------------------
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(fit_grid).to_csv(OUT_ROOT / "fit_grid.csv", index=False)
    pd.DataFrame(gains).to_csv(OUT_ROOT / "gains.csv", index=False)
    pd.DataFrame(block_shift).to_csv(OUT_ROOT / "block_shift.csv", index=False)
    pd.DataFrame(stability).to_csv(OUT_ROOT / "stability.csv", index=False)
    if curve_rows:
        pd.DataFrame(curve_rows).to_csv(OUT_ROOT / "curve_debias.csv", index=False)

    rep.line(f"\nCSV rollups -> {display_path(OUT_ROOT)}")
    rep.line("  fit_grid.csv  gains.csv  block_shift.csv  stability.csv"
             + ("  curve_debias.csv" if curve_rows else ""))
    rep.save(OUT_ROOT / "grid_compare.md", "Aging vs Drift -- 2x2 term grid")
    print(f"\nDone: {n_cells} (slice, mrc) cell(s) compared; "
          f"{len(curve_rows)} cross-mrc curve check(s).")


if __name__ == "__main__":
    main()
