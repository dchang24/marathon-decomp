"""Interactive race-factor dashboard with fit results + optional bootstrap CIs.

Displays the race factor v_j from ``results/models/`` for every slice×model
combination, with toggleable bootstrap error bars from
``results/bootstrap/`` when available.

Toggle axes
-----------
- **Cohort**: ALL / Po10
- **Sex**: M / W / B
- **Nu**: L2 & nu=6
- **Model terms**: A (aging) and D (athlete drift) — independent on/off buttons
  - A=off, D=off → ``baseline``
  - A=on,  D=off → ``agingS4gv`` (fits) / ``Axx`` (bootstrap)
  - A=on,  D=on  → ``AxD``
  - A=off, D=on  → not produced (greyed out)
- **MRC**: min-race-count 2 & 5
- **Error bars**: on/off (only available when bootstrap results exist)

Data sources
~~~~~~~~~~~~
- **Fits**: ``results/models/{slug}/{model_stem}__{hash}/fit.pkl``
  → race factors are extracted from the fit payload.
- **Bootstrap**: ``results/models/{slug}/{model_stem}__{hash}/bootstrap/``
  → percentile intervals from ``race_factors.parquet``.
- **Competitions**: ``data/race_results/v1/competitions.parquet``
  → race metadata (series_key, country, year).
- **External covariates**: ``results/analysis/covariate/01_merge/
  series_year_covariates.csv`` (from ``covariate/q01_merge.py``)
  → per-race air temp / WBGT (field, max) / total elevation gain, shown on hover
  and joined into the per-slice CSVs. Optional: if absent, hover/CSV omit them.

Outputs:
- ``results/visualizations/race_dashboard.html`` — the interactive page.
- ``results/visualizations/race_dashboard_csv/{slice_model}.csv`` — one CSV per
  slice/model with the numbers behind the bars (v_j, n_j, bootstrap CI, covariates).

Run::
    python -m scripts.visualizations.race_dashboard
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))

MODELS_DIR = _PROJECT_ROOT / "results" / "models"
COMPS_PATH = _PROJECT_ROOT / "data" / "race_results" / "competitions.parquet"
OUT_DIR = _PROJECT_ROOT / "results" / "visualizations"
# per series+year external covariates (written by covariate/q01_merge.py)
COV_PATH = (_PROJECT_ROOT / "results" / "analysis" / "covariate" / "01_merge"
            / "series_year_covariates.csv")
# one CSV per slice/model with the numbers shown in the dashboard
CSV_DIR = OUT_DIR / "race_dashboard_csv"
# external-covariate columns merged into hover + the per-slice CSVs
COV_COLS = ["temp_field", "wbgt_field", "wbgt_max", "total_gain_m", "net_gain_m"]

_slice_cache: dict = {}

# Minutes-at-a-3h-marathon scaling: v is in log-time, so the exact effect on a
# 3:00:00 (180 min) finish of a race factor v is 180 * (exp(v) - 1) minutes (used
# everywhere a difficulty *level* is shown, matching the paper); the first-order
# 180 * v is only used for log *contrasts* (differences of two v's).
REF_MARATHON_MIN = 180.0


# ── beta=0 APC gauge ────────────────────────────────────────────────
# Fits saved by fitter <= 1.2.0 sit OFF the beta=0 manifold: their per-iteration
# apply_apc_gauge_beta0() silently no-opped (a _t_j units bug, fixed in fitter
# 1.3.0; see src/marathon_decomp/models/model.py), so v_j carries a prior-pinned
# secular date tilt (slope(v ~ t) != 0). We re-impose beta=0 on the point fit and
# on every bootstrap replicate at read time, so the dashboard shows the
# era-relative production estimand. Idempotent (c ~= 0) once fits are >= 1.3.0.
def _frac_year(dates) -> np.ndarray:
    rd = pd.DatetimeIndex(pd.to_datetime(dates))
    return (rd.year + (rd.dayofyear - 1) / 365.25).to_numpy(np.float64)


def _beta0_gauge(v: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Subtract OLS date slope from v, re-centre to mean 0 (prediction-invariant)."""
    v = np.asarray(v, np.float64); t = np.asarray(t, np.float64)
    tc = t - t.mean(); denom = float(tc @ tc)
    if denom <= 0.0:
        return v - v.mean()
    c = float(((v - v.mean()) @ tc) / denom)
    out = v - c * t
    return out - out.mean()


def _beta0_gauge_rows(B: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Row-wise :func:`_beta0_gauge` for a (R, n) replicate matrix."""
    B = np.asarray(B, np.float64); t = np.asarray(t, np.float64)
    tc = t - t.mean(); denom = float(tc @ tc)
    if denom <= 0.0:
        return B - B.mean(axis=1, keepdims=True)
    c = ((B - B.mean(axis=1, keepdims=True)) @ tc) / denom
    out = B - np.outer(c, t)
    return out - out.mean(axis=1, keepdims=True)


_RACE_FRAC_YEAR: dict[int, float] | None = None


def _frac_year_for_races(race_ids: np.ndarray) -> np.ndarray:
    global _RACE_FRAC_YEAR
    if _RACE_FRAC_YEAR is None:
        c = pd.read_parquet(COMPS_PATH, columns=["race_id", "date"])
        _RACE_FRAC_YEAR = dict(zip(c["race_id"].to_numpy(),
                                   _frac_year(c["date"].to_numpy())))
    return np.array([_RACE_FRAC_YEAR[int(r)] for r in race_ids], np.float64)

PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
]

# percentile interval levels exposed in the page: label -> (lo_q, hi_q)
CI_LEVELS = {"95": (0.025, 0.975), "90": (0.05, 0.95), "68": (0.16, 0.84)}

# ── Model naming conventions ────────────────────────────────────────
# The user's toggle axes are (A, D).  The fits directory naming is:
#   baseline  → A=off, D=off   (use_s=false, use_d=false, no aging)
#   agingS4gv   → A=on,  D=off   (use_phi12=true, use_gamma=true, use_d=false)
#   AxD       → A=on,  D=on    (use_phi12=true, use_gamma=true, use_d=true)
# The bootstrap directory naming is:
#   baseline  → A=off, D=off
#   Axx       → A=on,  D=off  (source_fit_dir points to agingS4gv)
# The nu suffix is e.g. ``_L2`` or ``_nu8p00``.

# Model key = (A: bool, D: bool)
MODEL_PREFIXES_FITS = {
    (False, False): "baseline",
    (True, False):  "agingS4gv",
    (False, True):  "drift",
    (True, True):   "full",
}

MODEL_PREFIXES_BOOT = {
    (False, False): "baseline",
    (True, False):  "Axx",
    # (True, True): AxD bootstrap — not yet produced, but supported if it appears
}


def _nutag_label(tag: str) -> tuple[str, float]:
    """Parse 'L2' or 'nu8p00' into (display_label, nu_float)."""
    if tag == "L2":
        return "L2", float("inf")
    if tag.startswith("nu"):
        num = tag[2:].replace("p", ".")
        nu = float(num)
        num_str = num.rstrip("0").rstrip(".")
        return f"nu={num_str}", nu
    return tag, float("inf")


def _slice_from_slug(slug: str) -> tuple[str, str, str, int]:
    """Parse ``ALL_M_14-25_mrc2`` → (cohort='ALL', sex='M', yr='14-25', mrc=2)."""
    parts = slug.split("_")
    cohort = parts[0]
    sex = parts[1] if len(parts) >= 2 else "M"
    yr = parts[2] if len(parts) >= 3 else "14-25"
    mrc = 2  # default
    for p in parts:
        if p.startswith("mrc"):
            mrc = int(p[3:])
    return cohort, sex, yr, mrc


def _yr_to_dates(yr: str) -> tuple[int, int]:
    """Parse '14-25' → (2014, 2025); '21-25' → (2021, 2025)."""
    lo_str, hi_str = yr.split("-")
    lo = int(lo_str) + (2000 if int(lo_str) < 100 else 0)
    hi = int(hi_str) + (2000 if int(hi_str) < 100 else 0)
    return lo, hi


def _model_key_from_prefix(prefix: str) -> tuple[bool, bool] | None:
    """Map a fit directory prefix to (aging, drift) booleans."""
    if prefix == "baseline":
        return (False, False)
    if prefix.startswith("agingS4gv") or prefix.startswith("aging"):
        return (True, False)
    if prefix == "drift":
        return (False, True)
    if prefix == "full" or prefix == "AxD":
        return (True, True)
    if prefix == "Axx":
        return (True, False)
    return None


def _load_slice_cached(spec):
    key = repr(spec)
    if key not in _slice_cache:
        from marathon_decomp import load_slice
        _slice_cache[key] = load_slice(spec)
    return _slice_cache[key]



# ── Discover fits ────────────────────────────────────────────────────

def discover_fits() -> dict[str, dict[str, dict[str, Path]]]:
    """Walk FITS_DIR -> {slug: {model_tag: {nu_label: fit_dir}}}.

    model_tag is one of 'baseline', 'agingS4gv', 'AxD'.
    nu_label is 'L2' or 'nu=8'.
    """
    if not MODELS_DIR.is_dir():
        return {}
    found: dict[str, dict[str, dict[str, Path]]] = {}
    target_nus = {"L2", "nu8p00"}
    target_prefixes = {"baseline", "agingS4gv", "drift", "full", "AxD"}

    for slug_dir in sorted(p for p in MODELS_DIR.iterdir() if p.is_dir()):
        slug = slug_dir.name
        for fit_dir in sorted(p for p in slug_dir.iterdir() if p.is_dir()):
            if not (fit_dir / "fit.pkl").is_file():
                continue
            # stem before hash: e.g. "baseline_L2_anderson" or "agingS4gv_nu6p00_anderson"
            stem = fit_dir.name.split("__")[0]
            parts = stem.split("_")
            if len(parts) < 3:
                continue
            prefix = parts[0]
            nutag = parts[1]
            # solver = parts[2]  (anderson / als)

            if prefix not in target_prefixes:
                continue
            if nutag not in target_nus:
                continue

            # Prefer anderson solver; skip als if anderson already found
            nu_label, _ = _nutag_label(nutag)
            existing = found.get(slug, {}).get(prefix, {}).get(nu_label)
            if existing is not None:
                # keep anderson over als
                if "anderson" in fit_dir.name or "best" in fit_dir.name:
                    pass  # overwrite
                else:
                    continue

            found.setdefault(slug, {}).setdefault(prefix, {})[nu_label] = fit_dir
    return found


def load_fit_v(fit_dir: Path) -> pd.DataFrame | None:
    """Load fit.pkl and return a DataFrame with race_id, v_point, n_j."""
    pkl_path = fit_dir / "fit.pkl"
    if not pkl_path.is_file():
        return None
    with pkl_path.open("rb") as f:
        payload = pickle.load(f)
    v = payload.get("params", {}).get("v")
    spec = payload.get("spec")
    if v is None or spec is None:
        return None
    fd = _load_slice_cached(spec)
    counts = np.bincount(fd.col_idx, minlength=fd.J)
    v_g = _beta0_gauge(np.asarray(v, float), _frac_year(fd.race_date))
    return pd.DataFrame({
        "race_id": fd.race_ids.astype(int),
        "v_point": v_g,
        "n_j": counts.astype(int)
    })


def build_bootstrap_group(group_dir: Path) -> pd.DataFrame | None:
    """Per-race point estimate + percentile intervals from one bootstrap group."""
    rf_path = group_dir / "race_factors.parquet"
    if not rf_path.is_file():
        return None
    rf = pd.read_parquet(rf_path)
    # Pivot to (run_id, race_id) so each replicate can be re-gauged to beta=0 by
    # its OWN date slope before percentiles are taken (removes the wandering APC
    # tilt of <= 1.2.0 fits; otherwise it inflates early/late-window CIs).
    wide = rf.pivot(index="run_id", columns="race_id", values="v")
    race_ids = wide.columns.to_numpy()
    t = _frac_year_for_races(race_ids)
    boot = wide.loc[wide.index > 0]
    if boot.empty:
        return None
    B = boot.to_numpy(np.float64)
    if not np.isfinite(B).all():            # rare missing cell -> column mean
        B = np.where(np.isfinite(B), B, np.nanmean(B, axis=0))
    Bg = _beta0_gauge_rows(B, t)            # (R, n) gauged replicates

    cols = {
        "race_id": race_ids,
        "v_boot_med": np.median(Bg, axis=0),
        "v_sd": np.std(Bg, axis=0, ddof=1),
        "n_boot": np.full(len(race_ids), Bg.shape[0]),
    }
    for lbl, (qlo, qhi) in CI_LEVELS.items():
        cols[f"lo{lbl}"] = np.quantile(Bg, qlo, axis=0)
        cols[f"hi{lbl}"] = np.quantile(Bg, qhi, axis=0)
    return pd.DataFrame(cols)


# ── Master axis (shared across all payloads) ────────────────────────

def _fmt_num(x, fmt: str) -> str:
    """Format a possibly-NaN number for hover; em-dash when missing."""
    return fmt.format(x) if pd.notna(x) else "—"


def build_master_axis(all_race_ids: set[int], comps: pd.DataFrame,
                      cov: pd.DataFrame | None = None) -> dict:
    """Build a single canonical x-axis from the union of all race_ids.

    Returns a dict with races, years, countries, colors, ticktext,
    country_breaks — these are shared across every payload so that each
    race keeps a fixed position and colour regardless of which slice is
    shown. When ``cov`` (per-race external covariates) is supplied, also
    returns formatted ``cov_temp`` / ``cov_wbgt`` / ``cov_gain`` string arrays
    aligned to the axis, for the hover tooltip.
    """
    # race_id -> (air-temp, "WBGT field / max", total-gain) display strings
    cov_map: dict[int, tuple[str, str, str]] = {}
    if cov is not None:
        for _, c in cov.iterrows():
            cov_map[int(c["race_id"])] = (
                _fmt_num(c.get("temp_field"), "{:.1f}"),
                f"{_fmt_num(c.get('wbgt_field'), '{:.1f}')} / "
                f"{_fmt_num(c.get('wbgt_max'), '{:.1f}')}",
                _fmt_num(c.get("total_gain_m"), "{:.0f}"),
            )
    # Build metadata for every race in the union
    ids_df = pd.DataFrame({"race_id": sorted(all_race_ids)})
    m = ids_df.merge(comps, on="race_id", how="left")
    m["country"] = m["country"].fillna("ZZZ")
    m["series_key"] = m["series_key"].astype(str)
    m["year"] = m["year"].fillna(0).astype(int)
    m = m.sort_values(["country", "series_key", "year"]).reset_index(drop=True)

    series_unique = sorted(m["series_key"].unique())
    cmap = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(series_unique)}
    m["color"] = m["series_key"].map(cmap)

    keys, years, countries, colors, series_list = [], [], [], [], []
    cov_temp, cov_wbgt, cov_gain = [], [], []
    total_runners = []
    race_id_order: list[int] = []  # parallel array of race_ids (gaps get -1)
    breaks: list[int] = []
    prev = None
    for _, r in m.iterrows():
        if prev is not None and r["country"] != prev:
            keys.append(f"__gap_{len(breaks)}__")
            years.append(0); countries.append(""); colors.append("rgba(0,0,0,0)")
            series_list.append("")
            race_id_order.append(-1)
            cov_temp.append("—"); cov_wbgt.append("—"); cov_gain.append("—")
            total_runners.append(None)
            breaks.append(len(keys) - 1)
        keys.append(f"{r['series_key']}_{int(r['year'])}")
        years.append(int(r["year"])); countries.append(r["country"])
        colors.append(r["color"]); series_list.append(r["series_key"])
        race_id_order.append(int(r["race_id"]))
        ct, cw, cg = cov_map.get(int(r["race_id"]), ("—", "—", "—"))
        cov_temp.append(ct); cov_wbgt.append(cw); cov_gain.append(cg)
        total_runners.append(float(r["total_runners"]) if "total_runners" in r and pd.notna(r["total_runners"]) else None)
        prev = r["country"]

    ticktext = [""] * len(keys)
    i = 0
    while i < len(series_list):
        if series_list[i] == "":
            i += 1
            continue
        j = i
        while j + 1 < len(series_list) and series_list[j + 1] == series_list[i]:
            j += 1
        ticktext[(i + j) // 2] = series_list[i].replace("_", " ")
        i = j + 1

    # Build race_id → master index lookup (for align_to_master)
    rid_to_idx = {rid: idx for idx, rid in enumerate(race_id_order) if rid >= 0}

    return {
        "races": keys, "years": years, "countries": countries,
        "colors": colors, "ticktext": ticktext, "country_breaks": breaks,
        "ref_min": REF_MARATHON_MIN,
        "cov_temp": cov_temp, "cov_wbgt": cov_wbgt, "cov_gain": cov_gain,
        "total_runners": total_runners,
        # internal, not serialised to JS:
        "_race_id_order": race_id_order,
        "_rid_to_idx": rid_to_idx,
        "_n": len(keys),
    }


def align_to_master(df: pd.DataFrame, master: dict, *,
                    has_boot: bool = False) -> dict:
    """Align a per-fit DataFrame to the master axis, filling null for absent races."""
    rid_to_idx = master["_rid_to_idx"]
    n = master["_n"]

    # Index the dataframe by race_id for O(1) lookup
    df_indexed = df.set_index("race_id")

    val_cols = ["v_point", "n_j"]
    ci_cols = [f"lo{l}" for l in CI_LEVELS] + [f"hi{l}" for l in CI_LEVELS]
    boot_cols = ["v_boot_med", "v_sd", "n_boot"]
    if has_boot:
        val_cols += boot_cols + ci_cols

    arrs: dict[str, list] = {c: [None] * n for c in val_cols}

    for rid, idx in rid_to_idx.items():
        if rid not in df_indexed.index:
            continue
        row = df_indexed.loc[rid]
        # handle duplicate race_ids (shouldn't happen, but guard)
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        for c in val_cols:
            if c not in row.index:
                continue
            val = row[c]
            if pd.isna(val):
                continue
            if c in ("n_j", "n_boot"):
                arrs[c][idx] = int(val)
            else:
                arrs[c][idx] = float(val)

    result = {
        "v_point": arrs["v_point"],
        "n_j": arrs["n_j"],
        "has_boot": has_boot,
    }
    if has_boot:
        result["v_boot_med"] = arrs.get("v_boot_med")
        result["v_sd"] = arrs.get("v_sd")
        result["n_boot"] = arrs.get("n_boot")
        for l in CI_LEVELS:
            result[f"lo{l}"] = arrs.get(f"lo{l}")
            result[f"hi{l}"] = arrs.get(f"hi{l}")
    return result


# ── HTML template ────────────────────────────────────────────────────

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Race-factor dashboard</title>
<meta name="description" content="Interactive race-factor v_j dashboard with model term toggles and optional bootstrap confidence intervals">
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  :root {
    --bg: #f0f2f5; --surface: #ffffff; --surface-alt: #f8f9fb;
    --border: #e2e5ea; --border-light: #eef0f3;
    --text: #1a1d23; --text-muted: #6b7280; --text-dim: #9ca3af;
    --accent: #4f46e5; --accent-light: #eef2ff; --accent-border: #c7d2fe;
    --accent-hover: #4338ca;
    --on-btn: #4f46e5; --on-bg: #eef2ff; --on-border: #a5b4fc;
    --off-btn: #6b7280; --off-bg: #f3f4f6; --off-border: #d1d5db;
    --green: #059669; --green-bg: #ecfdf5; --green-border: #a7f3d0;
    --red: #dc2626; --red-bg: #fef2f2; --red-border: #fecaca;
    --orange: #d97706; --orange-bg: #fffbeb;
    --shadow-sm: 0 1px 2px rgba(0,0,0,.05);
    --shadow-md: 0 4px 6px -1px rgba(0,0,0,.07), 0 2px 4px -2px rgba(0,0,0,.05);
    --radius: 6px; --radius-lg: 10px;
    --font: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
    --mono: 'JetBrains Mono', 'Cascadia Code', 'Fira Code', 'Consolas', monospace;
  }
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

  * { box-sizing: border-box; margin: 0; padding: 0; }
  html, body { height: 100%; overflow: hidden; }
  body {
    font-family: var(--font); font-size: 13px; color: var(--text);
    background: var(--bg); display: flex; flex-direction: column;
  }

  /* ─── Header bar ─── */
  .header {
    background: var(--surface); border-bottom: 1px solid var(--border);
    box-shadow: var(--shadow-sm);
    padding: 10px 16px; display: flex; flex-wrap: wrap;
    align-items: flex-end; gap: 12px 20px; flex-shrink: 0;
  }
  .header-title h1 {
    font-size: 15px; font-weight: 700; color: var(--accent);
    letter-spacing: -0.3px;
  }
  .header-title .subtitle {
    font-size: 10px; color: var(--text-dim); text-transform: uppercase;
    font-weight: 600; letter-spacing: 0.5px; margin-top: 1px;
  }

  /* ─── Control groups ─── */
  .controls { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: flex-end; }
  .ctrl-group {
    display: flex; flex-direction: column; gap: 2px;
  }
  .ctrl-label {
    font-size: 9px; font-weight: 700; text-transform: uppercase;
    letter-spacing: 0.6px; color: var(--text-dim);
  }
  .ctrl-sep { border-left: 1px solid var(--border-light); margin: 0 4px; align-self: stretch; }

  /* ─── Toggle button pills ─── */
  .toggle-row { display: flex; gap: 3px; }
  .toggle-btn {
    font-family: var(--font); font-size: 11px; font-weight: 600;
    padding: 3px 10px; border-radius: 4px; cursor: pointer;
    border: 1.5px solid var(--off-border); background: var(--off-bg);
    color: var(--off-btn); transition: all 0.15s ease;
    user-select: none;
  }
  .toggle-btn:hover { border-color: var(--accent); color: var(--accent); }
  .toggle-btn.active {
    background: var(--on-bg); color: var(--on-btn);
    border-color: var(--on-border); box-shadow: 0 0 0 1px var(--accent-light);
  }
  .toggle-btn.unavail {
    opacity: 0.35; cursor: not-allowed;
    pointer-events: none;
  }

  /* ─── Selects ─── */
  select {
    font-family: var(--font); font-size: 11px; font-weight: 500;
    border: 1.5px solid var(--border); border-radius: 4px;
    padding: 3px 6px; background: var(--surface-alt); color: var(--text);
    cursor: pointer; transition: border-color 0.15s;
  }
  select:hover { border-color: var(--accent); }
  select:focus { outline: none; border-color: var(--accent); box-shadow: 0 0 0 2px var(--accent-light); }

  /* ─── Checkbox ─── */
  .chk-label {
    display: flex; align-items: center; gap: 5px;
    font-size: 11px; color: var(--text-muted); font-weight: 500;
    cursor: pointer;
  }
  .chk-label input[type="checkbox"] { accent-color: var(--accent); cursor: pointer; }

  /* ─── Info panel ─── */
  .info-panel {
    margin-left: auto; background: var(--accent-light);
    border: 1px solid var(--accent-border); border-radius: var(--radius);
    padding: 4px 12px; font-family: var(--mono); font-size: 10px;
    color: var(--accent-hover); min-width: 320px; max-height: 52px;
    overflow: hidden; display: flex; flex-direction: column; justify-content: center;
    line-height: 1.5;
  }

  /* ─── Chart area ─── */
  .chart-wrap {
    flex: 1; padding: 10px; min-height: 0;
  }
  #chart {
    background: var(--surface); border-radius: var(--radius-lg);
    box-shadow: var(--shadow-md); border: 1px solid var(--border-light);
    width: 100%; height: 100%;
  }
</style>
</head>
<body>

<div class="header">
  <div class="header-title">
    <h1>Race-factor dashboard</h1>
    <span class="subtitle">
      v<sub>j</sub> from fits &middot; bar = point fit &middot; optional bootstrap percentile CI
    </span>
  </div>

  <div class="ctrl-sep"></div>

  <div class="controls">
    <!-- Cohort -->
    <div class="ctrl-group">
      <span class="ctrl-label">Cohort</span>
      <div class="toggle-row" id="cohort-btns">
        <button class="toggle-btn active" data-val="ALL">ALL</button>
        <button class="toggle-btn" data-val="Po10">Po10</button>
        <button class="toggle-btn" data-val="WA">WA</button>
      </div>
    </div>

    <!-- Sex -->
    <div class="ctrl-group">
      <span class="ctrl-label">Sex</span>
      <div class="toggle-row" id="sex-btns">
        <button class="toggle-btn" data-val="M">M</button>
        <button class="toggle-btn" data-val="W">W</button>
        <button class="toggle-btn active" data-val="B">B</button>
      </div>
    </div>

    <!-- Nu -->
    <div class="ctrl-group">
      <span class="ctrl-label">Nu</span>
      <div class="toggle-row" id="nu-btns">
        <button class="toggle-btn" data-val="L2">L2</button>
        <button class="toggle-btn active" data-val="nu=8">ν=8</button>
      </div>
    </div>

    <!-- Year range -->
    <div class="ctrl-group">
      <span class="ctrl-label">Years</span>
      <div class="toggle-row" id="yr-btns">
        <button class="toggle-btn" data-val="14-20">2014–20</button>
        <button class="toggle-btn active" data-val="14-25">2014–25</button>
        <button class="toggle-btn" data-val="21-25">2021–25</button>
      </div>
    </div>

    <div class="ctrl-sep"></div>

    <!-- Model terms: A and D -->
    <div class="ctrl-group">
      <span class="ctrl-label">Model terms</span>
      <div class="toggle-row">
        <button class="toggle-btn" id="btn-A" onclick="toggleTerm('A')">A</button>
        <button class="toggle-btn" id="btn-D" onclick="toggleTerm('D')">D</button>
      </div>
    </div>

    <!-- MRC -->
    <div class="ctrl-group">
      <span class="ctrl-label">MRC</span>
      <div class="toggle-row" id="mrc-btns">
        <button class="toggle-btn active" data-val="2">2</button>
        <button class="toggle-btn" data-val="5">5</button>
      </div>
    </div>

    <div class="ctrl-sep"></div>

    <!-- Display options -->
    <div class="ctrl-group">
      <span class="ctrl-label">Display</span>
      <select id="disp-sel">
        <option value="log">v_j (log-time)</option>
        <option value="mult" selected>multiplier exp(v−mean)</option>
        <option value="min">minutes @ 3h marathon</option>
      </select>
    </div>

    <div class="ctrl-group">
      <span class="ctrl-label">CI level</span>
      <select id="ci-sel">
        <option value="95" selected>95%</option>
        <option value="90">90%</option>
        <option value="68">68%</option>
      </select>
    </div>

    <div class="ctrl-group">
      <span class="ctrl-label">y range</span>
      <select id="range-sel">
        <option value="auto">auto</option>
        <option value="0.05">±0.05</option>
        <option value="0.10">±0.10</option>
        <option value="0.15" selected>±0.15</option>
        <option value="0.20">±0.20</option>
      </select>
    </div>

    <label class="chk-label">
      <input id="err-chk" type="checkbox" checked> error bars
    </label>
  </div>

  <div class="info-panel" id="info-panel"></div>
</div>

<div class="chart-wrap">
  <div id="chart"></div>
</div>

<script>
// ─── Embedded data ──────────────────────────────────────────────────
// MASTER holds the shared x-axis: races, years, countries, colors, ticktext, country_breaks, ref_min
const MASTER = __MASTER_JSON__;
const FITS = __FITS_JSON__;
// FITS[dataKey] = {v_point, n_j, has_boot}  — value arrays aligned to MASTER.races
const BOOT = __BOOT_JSON__;
// BOOT[dataKey] = {v_boot_med, v_sd, n_boot, lo95, hi95, lo90, hi90, lo68, hi68}
// dataKey format: "{cohort}_{sex}_{yr}_{mrc}_{modelPrefix}_{nuLabel}"
const AVAIL_KEYS = __AVAIL_KEYS_JSON__;  // Set of available fit data keys
const BOOT_KEYS  = __BOOT_KEYS_JSON__;   // Set of keys that have bootstrap data

// ─── State ──────────────────────────────────────────────────────────
let state = {
  cohort: 'ALL', sex: 'B', nu: 'nu=8',
  A: true, D: true,
  yr: '14-25', mrc: '2',
};

function modelPrefix() {
  if (!state.A && !state.D) return 'baseline';
  if (state.A && !state.D) return 'agingS4gv';
  if (!state.A && state.D) return 'drift';
  if (state.A && state.D)  return 'full';
  return null;
}

function dataKey() {
  const pfx = modelPrefix();
  if (!pfx) return null;
  return `${state.cohort}_${state.sex}_${state.yr}_${state.mrc}_${pfx}_${state.nu}`;
}

// ─── Toggle button wiring ───────────────────────────────────────────
function wireToggleRow(containerId, stateKey) {
  const btns = document.querySelectorAll(`#${containerId} .toggle-btn`);
  btns.forEach(btn => {
    btn.addEventListener('click', () => {
      btns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      state[stateKey] = btn.dataset.val;
      updateUI();
    });
  });
}
wireToggleRow('cohort-btns', 'cohort');
wireToggleRow('sex-btns', 'sex');
wireToggleRow('nu-btns', 'nu');
wireToggleRow('yr-btns', 'yr');
wireToggleRow('mrc-btns', 'mrc');

function toggleTerm(term) {
  state[term] = !state[term];
  updateUI();
}

function updateUI() {
  // Update A/D button styles
  const btnA = document.getElementById('btn-A');
  const btnD = document.getElementById('btn-D');
  btnA.classList.toggle('active', state.A);
  btnD.classList.toggle('active', state.D);

  render();
}

// ─── Rendering ──────────────────────────────────────────────────────
function autoRange(values, padFrac, around) {
  const finite = values.filter(v => Number.isFinite(v));
  if (!finite.length) return null;
  let lo = Math.min(...finite), hi = Math.max(...finite);
  const pad = (hi - lo) * padFrac || Math.abs(hi) * padFrac || 0.01;
  if (around !== undefined) {
    const mm = Math.max(Math.abs(hi - around), Math.abs(around - lo)) + pad;
    return [around - mm, around + mm];
  }
  return [lo - pad, hi + pad];
}

function layout(yRange, center, ylab) {
  const breaks = MASTER.country_breaks || [];
  const N = MASTER.races.length;
  const shapes = breaks.map(idx => ({
    type:'line', xref:'x', yref:'y domain', x0:idx, x1:idx, y0:0, y1:1,
    line:{color:'#bbb', width:0.8, dash:'dash'}, layer:'above',
  }));
  shapes.push({ type:'line', xref:'x domain', yref:'y',
    x0:0, x1:1, y0:center, y1:center,
    line:{color:'#888',width:1,dash:'dash'}, layer:'above' });
  return {
    margin:{l:70,r:20,t:10,b:160},
    plot_bgcolor:'#fafafa', paper_bgcolor:'white', showlegend:false, bargap:0.15,
    font:{family:'Inter, ui-sans-serif, system-ui, sans-serif', size:10},
    shapes,
    xaxis:{ type:'category', categoryorder:'array', categoryarray:MASTER.races,
            tickmode:'array', tickvals:MASTER.races, ticktext:MASTER.ticktext,
            tickangle:40, range:[-0.5, N-0.5], gridcolor:'#eee' },
    yaxis:{ title:{text:ylab, standoff:8}, range:yRange, zeroline:false, gridcolor:'#eee' },
  };
}

function render() {
  const key = dataKey();
  const info = document.getElementById('info-panel');

  if (!key || !AVAIL_KEYS.includes(key)) {
    Plotly.purge('chart');
    const pfx = modelPrefix();
    info.innerHTML = pfx === null
      ? '<b>D without A</b> is not a valid model — enable A first'
      : `no fit data for <b>${key || '?'}</b>`;
    return;
  }

  const p = FITS[key];
  const hasBoot = BOOT_KEYS.includes(key);
  const boot = hasBoot ? BOOT[key] : null;

  const disp    = document.getElementById('disp-sel').value;
  const ciLvl   = document.getElementById('ci-sel').value;
  const showErr = document.getElementById('err-chk').checked && hasBoot;

  const center0 = p.v_point;
  const lo = (boot && boot['lo'+ciLvl]) || [];
  const hi = (boot && boot['hi'+ciLvl]) || [];

  // Transform
  let meanLog = 0;
  if (disp === 'mult') {
    const pres = center0.filter(Number.isFinite);
    meanLog = pres.length ? pres.reduce((a,b)=>a+b,0)/pres.length : 0;
  }
  const tf = x => (x===null||!Number.isFinite(x)) ? null
              : disp==='mult' ? Math.exp(x-meanLog)
              : disp==='min'  ? MASTER.ref_min*(Math.exp(x)-1)
              : x;
  const center = disp==='mult' ? 1 : 0;
  const ylab = disp==='mult' ? 'v as multiplier exp(v_j − mean)'
             : disp==='min'  ? 'race effect on a 3:00:00 finish (minutes)'
             : 'v_j (race difficulty, log-time units)';

  const y    = center0.map(tf);
  const yLo  = lo.map(tf);
  const yHi  = hi.map(tf);
  const ePlus  = y.map((v,i) => (showErr && Number.isFinite(v) && Number.isFinite(yHi[i])) ? Math.max(0,yHi[i]-v) : 0);
  const eMinus = y.map((v,i) => (showErr && Number.isFinite(v) && Number.isFinite(yLo[i])) ? Math.max(0,v-yLo[i]) : 0);

  // Custom data for hover
  const cd = MASTER.races.map((r,i)=>[
    r, MASTER.years[i]??'', MASTER.countries[i]??'',
    (p.n_j && Number.isFinite(p.n_j[i]))? p.n_j[i].toLocaleString():'—',
    hasBoot && boot.n_boot ? (boot.n_boot[i]??'—') : '—',
    Number.isFinite(p.v_point[i])? p.v_point[i].toFixed(4):'—',
    hasBoot && boot.v_boot_med ? (Number.isFinite(boot.v_boot_med[i])? boot.v_boot_med[i].toFixed(4):'—') : '—',
    Number.isFinite(lo[i])? lo[i].toFixed(4):'—',
    Number.isFinite(hi[i])? hi[i].toFixed(4):'—',
    hasBoot && boot.v_sd ? (Number.isFinite(boot.v_sd[i])? boot.v_sd[i].toExponential(2):'—') : '—',
    (Number.isFinite(hi[i])&&Number.isFinite(lo[i]))? (MASTER.ref_min*(Math.exp(hi[i])-Math.exp(lo[i]))).toFixed(2):'—',
    MASTER.cov_temp ? (MASTER.cov_temp[i] ?? '—') : '—',
    MASTER.cov_wbgt ? (MASTER.cov_wbgt[i] ?? '—') : '—',
    MASTER.cov_gain ? (MASTER.cov_gain[i] ?? '—') : '—',
    (MASTER.total_runners && Number.isFinite(MASTER.total_runners[i]))? MASTER.total_runners[i].toLocaleString():'—',
  ]);

  const BARGAP = 0.15;
  const plotW = Math.max(0, (document.getElementById('chart').clientWidth || 800) - 90);
  const barPx = plotW / Math.max(1, MASTER.races.length) * (1 - BARGAP);
  const capPx = Math.max(1, 0.5 * barPx);

  const hoverBase = '<b>%{x}</b> (%{customdata[1]}, %{customdata[2]})<br>'
    + 'runners n_j: %{customdata[3]} (total: %{customdata[14]})<br>'
    + 'v_j: %{customdata[5]}';
  // external covariates (air temp / WBGT field-max / total elevation gain)
  const covLine = '<br>air temp: %{customdata[11]} C'
    + ' &nbsp; WBGT f/m: %{customdata[12]} C'
    + '<br>elevation gain: %{customdata[13]} m';
  const hoverBoot = hoverBase
    + ' &nbsp; median: %{customdata[6]}<br>'
    + ciLvl + '% CI: [%{customdata[7]}, %{customdata[8]}] '
    + '(sd=%{customdata[9]})<br>'
    + 'CI width ≈ %{customdata[10]} min @ 3:00'
    + covLine + '<extra></extra>';
  const hoverNoBoot = hoverBase + covLine + '<extra></extra>';

  const trace = {
    x:MASTER.races, y, type:'bar',
    marker:{ color: MASTER.colors },
    customdata:cd,
    error_y:{ type:'data', symmetric:false, array:ePlus, arrayminus:eMinus,
              color:'rgba(40,40,40,0.55)', thickness:1, width:capPx, visible:showErr },
    hovertemplate: hasBoot ? hoverBoot : hoverNoBoot,
  };

  const rangeSel = document.getElementById('range-sel').value;
  let yRange;
  if (rangeSel === 'auto') {
    const env = showErr ? yLo.concat(yHi).concat(y) : y;
    yRange = autoRange(env, 0.10, center) ?? (disp==='mult'?[0.9,1.1]:[-0.2,0.2]);
  } else {
    const h = parseFloat(rangeSel);
    yRange = disp==='mult' ? [1-h, 1+h] : disp==='min' ? [MASTER.ref_min*(Math.exp(-h)-1), MASTER.ref_min*(Math.exp(h)-1)] : [-h, h];
  }

  Plotly.react('chart', [trace], layout(yRange, center, ylab), {responsive:true});

  // Info panel
  const nPresent = p.v_point.filter(v => v !== null && Number.isFinite(v)).length;
  const nTotal = MASTER.races.filter(r=>!r.startsWith('__gap_')).length;
  const pfx = modelPrefix();
  let infoHtml = `<b>${state.cohort}_${state.sex}</b> yr=${state.yr} mrc=${state.mrc} &nbsp; `
    + `model: <b>${pfx}</b> &nbsp; ν: <b>${state.nu}</b> &nbsp; `
    + `races: <b>${nPresent}</b>/${nTotal}`;
  if (hasBoot) {
    const widths = lo.map((l,i) => (Number.isFinite(l)&&Number.isFinite(hi[i]))? hi[i]-l : null)
                     .filter(Number.isFinite).sort((a,b)=>a-b);
    const medW = widths.length ? widths[Math.floor(widths.length/2)] : null;
    // exact minute width per race = 180*(exp(hi)-exp(lo)); median over races
    const minWidths = lo.map((l,i) => (Number.isFinite(l)&&Number.isFinite(hi[i]))? p.ref_min*(Math.exp(hi[i])-Math.exp(l)) : null)
                        .filter(Number.isFinite).sort((a,b)=>a-b);
    const medMinW = minWidths.length ? minWidths[Math.floor(minWidths.length/2)] : null;
    infoHtml += `<br>bootstrap CI (${ciLvl}%) median width: <b>${medW!==null? medW.toFixed(4):'—'}</b> log`
      + ` &nbsp;≈ <b>${medMinW!==null? medMinW.toFixed(1):'—'}</b> min @ 3:00`;
  } else {
    infoHtml += '<br><span style="color:#9ca3af">no bootstrap data for this combination</span>';
  }
  info.innerHTML = infoHtml;
}

// Wire up remaining controls
['disp-sel','ci-sel','range-sel'].forEach(id =>
  document.getElementById(id).addEventListener('change', render));
document.getElementById('err-chk').addEventListener('change', render);

// Initial render
updateUI();
</script>
</body>
</html>
"""


def write_slice_csvs(loaded: list[tuple[str, pd.DataFrame, bool, Path | None]],
                     comps: pd.DataFrame, cov: pd.DataFrame | None) -> None:
    """One CSV per slice/model with the numbers shown in the dashboard.

    Each row is a race: v_point + n_j (+ bootstrap median / sd / percentile
    intervals when present) + the external covariates. Written to CSV_DIR, one
    file per dataKey (``{cohort}_{sex}_{yr}_{mrc}_{model}_{nu}.csv``).
    """
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    cov_small = cov[["race_id"] + COV_COLS] if cov is not None else None
    # preferred column order (only those actually present are written)
    order = (["race_id", "series_key", "country", "year", "v_point", "n_j",
              "v_boot_med", "v_sd", "n_boot"]
             + [f"lo{l}" for l in CI_LEVELS] + [f"hi{l}" for l in CI_LEVELS]
             + COV_COLS)
    for dkey, df, _has_boot, _ in loaded:
        out = df.merge(comps, on="race_id", how="left")
        if cov_small is not None:
            out = out.merge(cov_small, on="race_id", how="left")
        cols = [c for c in order if c in out.columns]
        out = out[cols].sort_values(["country", "series_key", "year"])
        out.to_csv(CSV_DIR / f"{dkey}.csv", index=False)
    print(f"Wrote {len(loaded)} per-slice/model CSVs to {CSV_DIR}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    comps = pd.read_parquet(COMPS_PATH)[["race_id", "competition_id", "series_key", "country", "year"]]

    race_info_path = _PROJECT_ROOT / "data" / "misc" / "race_stats_summary.csv"
    if race_info_path.is_file():
        rinfo = pd.read_csv(race_info_path)[["comp_id", "total_runners"]]
        comps = comps.merge(rinfo, left_on="competition_id", right_on="comp_id", how="left")
    else:
        comps["total_runners"] = np.nan

    # external covariates (air temp / WBGT / elevation), optional
    if COV_PATH.is_file():
        cov = pd.read_csv(COV_PATH)
        print(f"Loaded covariates for {cov['race_id'].nunique()} races from {COV_PATH.name}")
    else:
        cov = None
        print(f"[warn] no covariate lookup at {COV_PATH}; hover/CSV will omit "
              f"covariates. Run scripts/05_analysis/covariate/q01_merge.py first.")

    print("Discovering fits...")
    fits_found = discover_fits()
    if not fits_found:
        raise SystemExit(f"No fits found under {MODELS_DIR}")


    # ── Pass 1: collect all race_ids + load DataFrames ───────────────
    _counts_cache: dict[str, dict[int, int]] = {}
    all_race_ids: set[int] = set()
    # Store (dkey, df, has_boot, boot_group_dir) for second pass
    loaded: list[tuple[str, pd.DataFrame, bool, Path | None]] = []

    for slug, models in sorted(fits_found.items()):
        cohort, sex, yr, mrc = _slice_from_slug(slug)
        short = f"{cohort}_{sex}"

        for model_prefix, nu_map in sorted(models.items()):
            mk = _model_key_from_prefix(model_prefix)
            if mk is None:
                continue

            for nu_label, fit_dir in sorted(nu_map.items()):
                dkey = f"{cohort}_{sex}_{yr}_{mrc}_{model_prefix}_{nu_label}"
                print(f"  Loading fit: {dkey} <- {fit_dir.name}")

                df = load_fit_v(fit_dir)
                if df is None:
                    print("    [skip] could not load fit.pkl")
                    continue


                # Collect race_ids for the master axis
                all_race_ids.update(df["race_id"].astype(int).tolist())

                # Check bootstrap availability
                boot_group_dir = fit_dir / "bootstrap"
                has_boot = boot_group_dir.is_dir()
                if has_boot:
                    boot_df = build_bootstrap_group(boot_group_dir)
                    if boot_df is not None:
                        df = df.merge(boot_df, on="race_id", how="left")
                    else:
                        has_boot = False
                        boot_group_dir = None

                loaded.append((dkey, df, has_boot, boot_group_dir))

    if not loaded:
        raise SystemExit("No usable fits found.")

    # ── Write one CSV per slice/model (the numbers shown in the dashboard) ─
    write_slice_csvs(loaded, comps, cov)

    # ── Build master axis from union of all race_ids ─────────────────
    print(f"\nBuilding master axis from {len(all_race_ids)} unique races...")
    master = build_master_axis(all_race_ids, comps, cov)
    print(f"  Master axis: {master['_n']} slots "
          f"({len(all_race_ids)} races + {len(master['country_breaks'])} gaps)")

    # ── Pass 2: align each payload to the master axis ────────────────
    fits_data: dict[str, dict] = {}
    boot_data: dict[str, dict] = {}
    avail_keys: list[str] = []
    boot_keys: list[str] = []

    for dkey, df, has_boot, boot_group_dir in loaded:
        payload = align_to_master(df, master, has_boot=has_boot)
        fits_data[dkey] = payload
        avail_keys.append(dkey)

        if has_boot:
            n = master["_n"]
            boot_data[dkey] = {
                "v_boot_med": payload.get("v_boot_med"),
                "v_sd": payload.get("v_sd"),
                "n_boot": payload.get("n_boot"),
            }
            for l in CI_LEVELS:
                boot_data[dkey][f"lo{l}"] = payload.get(f"lo{l}")
                boot_data[dkey][f"hi{l}"] = payload.get(f"hi{l}")
            boot_keys.append(dkey)
            print(f"    + bootstrap CI from {boot_group_dir.name}")

        n_present = sum(1 for v in payload["v_point"] if v is not None)
        print(f"    {dkey}: {n_present}/{len(all_race_ids)} races"
              f"{' (+ bootstrap)' if has_boot else ''}")

    # ── Strip bootstrap fields from fits_data (they're in boot_data) ─
    for dkey in fits_data:
        for field in ["v_boot_med", "v_sd", "n_boot"] + [f"lo{l}" for l in CI_LEVELS] + [f"hi{l}" for l in CI_LEVELS]:
            fits_data[dkey].pop(field, None)

    # ── Serialise master axis (drop internal fields) ─────────────────
    master_js = {
        k: v for k, v in master.items() if not k.startswith("_")
    }

    # ── Emit HTML ────────────────────────────────────────────────────
    def _json(obj):
        return json.dumps(obj, allow_nan=False, default=lambda x: None)

    html = (HTML
            .replace("__MASTER_JSON__", _json(master_js))
            .replace("__FITS_JSON__", _json(fits_data))
            .replace("__BOOT_JSON__", _json(boot_data))
            .replace("__AVAIL_KEYS_JSON__", _json(avail_keys))
            .replace("__BOOT_KEYS_JSON__", _json(boot_keys)))

    out_path = OUT_DIR / "race_dashboard.html"
    out_path.write_text(html, encoding="utf-8")
    print(f"\nWrote {out_path}")
    print(f"  Open in a browser: file://{out_path.resolve()}")
    print(f"\nAvailable keys ({len(avail_keys)}):")
    for k in sorted(avail_keys):
        boot_tag = " [+bootstrap]" if k in boot_keys else ""
        print(f"  {k}{boot_tag}")


if __name__ == "__main__":
    main()
