"""Shared loaders for the covariate analysis (weather + course vs v_j).

Production port of ``scratch/covariate_analysis_v2``. Two questions, two axes:

  * **variable selection** (q02): which of the available physical covariates
    (weather + course) most strongly tracks race difficulty? Tested against the
    single ALL_B **full**-model ``v_j`` (the richest cohort, the headline fit).
  * **regression** (q03): how well do the chosen weather + course conditions
    explain ``v_j`` across the six production slices (ALL/Po10 x M/W/B) and how
    that compares against the naive 'simple' race metrics in
    ``race_stats_summary.csv``.

All ``v_j`` are reloaded (**no refit**) from ``results/models/{slug}/`` at the
production operating point (full model, nu=8, mrc=2). Each fit's ``v_j`` is
re-expressed under the **G1 beta=0 APC gauge** (subtract the secular linear year
trend) so the comparison is trend-free and consistent across slices. Higher
``v_j`` = slower = harder race.

Data sources are the production exports:
  * weather   : ``data/weather/marathon_weather_features.csv``   (per competition_id)
  * elevation : ``data/course_profile/marathon_year_elevation.csv`` (per series_key, year)
  * simple    : ``data/misc/race_stats_summary.csv``               (per competition_id)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))   # scripts/ (for *_common packages)

from marathon_decomp import load_slice, registry                       # noqa: E402
from marathon_decomp.config import (                                   # noqa: E402
    DATA_DIR, MISC_DIR, WEATHER_DIR, COURSE_PROFILE_DIR, RESULTS_DIR)

from baseline_common import slices as S                                # noqa: E402
from baseline_common.fitting import baseline_cfg                       # noqa: E402
from aging_common.fitting import aging_cfg, aging_stem                 # noqa: E402
from drift_common.fitting import drift_cfg, drift_stem                 # noqa: E402

MODELS_ROOT = RESULTS_DIR / "models"
OUT_ROOT = RESULTS_DIR / "analysis" / "covariate"

NU = 8.0
MRC = 2
MODEL = "full"                       # production AxD point fit
# the four nested models, in increasing complexity (registered per slice at nu=8).
# Used by q08 to compare how each model's v_j tracks external covariates.
MODELS = ["baseline", "aging", "drift", "full"]

# the six regression slices (ALL/Po10 x M/W/B); full model each
SLICE_ORDER = ["ALL_B", "ALL_M", "ALL_W", "Po10_B", "Po10_M", "Po10_W"]
# the headline slice used for variable selection (q02)
VAR_SELECT_SLICE = "ALL_B"
# slice tag for the multi-slice outputs (q01/q03/q04/q07 span all six slices)
ALL_SLICES_TAG = "6slices"


# --------------------------------------------------------------------------- #
# output layout: one subdirectory per script, every file slice-labeled         #
# --------------------------------------------------------------------------- #
def out_path(subdir: str, stem: str, slice_tag: str, ext: str) -> Path:
    """Resolve ``OUT_ROOT/<subdir>/<stem>__<slice_tag>.<ext>``.

    Each runnable script writes into its own ``subdir`` (e.g. ``02_variable_selection``)
    and every output name carries the ``slice_tag`` it was computed on (a single
    slice like ``ALL_B``, or ``ALL_SLICES_TAG`` for the multi-slice tables). The
    subdirectory is created on demand.
    """
    d = OUT_ROOT / subdir
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{stem}__{slice_tag}.{ext}"


# canonical merged-table location (written by q01, read by q02-q07 + p01)
MERGED_PATH = OUT_ROOT / "01_merge" / f"covariate_merged__{ALL_SLICES_TAG}.parquet"

# focused per-race covariate lookup written by q01: one row per series+year with
# air temp + WBGT (field/max) + course elevation gain. Human-readable reference,
# also consumed by scripts/visualizations/race_dashboard.py for hover info.
SERIES_YEAR_COV_PATH = OUT_ROOT / "01_merge" / "series_year_covariates.csv"

# the covariate columns surfaced in that lookup (kept in this order)
SERIES_YEAR_COV_COLS = [
    "race_id", "competition_id", "series_key", "country", "year", "course_id",
    "temp_field", "wbgt_field", "wbgt_max", "total_gain_m", "net_gain_m",
]

# --- candidate covariates for variable selection (q02) -------------------- #
# weather: physical thermal / wind / precip load (the QC-only wbgt_spatial_spread
# and wbgt_anchor_sens diagnostics are excluded; wbgt_rising_slope kept as a
# within-race trajectory descriptor).
WEATHER_CANDIDATES = [
    ("temp_field", "air temp (C)"),
    ("wbgt_field", "WBGT field (C)"),
    ("wbgt_max", "WBGT max (C)"),
    ("wind_field", "wind (m/s)"),
    ("precip_total_mm", "precip (mm)"),
    ("wbgt_rising_slope_c_per_h", "WBGT rising slope (C/h)"),
]
COURSE_CANDIDATES = [
    ("total_gain_m", "total elevation gain (m)"),
    ("net_gain_m", "net elevation gain (m)"),
    ("abs_net_gain_m", "|net elevation gain| (m)"),
]

# default predictor sets for the joint regression (q03): the best weather var
# (per q02) paired with course gain, plus the principled-thermal WBGT variant.
COURSE_PREDICTOR = "total_gain_m"
WEATHER_SETS = [("temp_field", "temp"), ("wbgt_max", "WBGT")]


# --------------------------------------------------------------------------- #
# gauge helpers                                                                #
# --------------------------------------------------------------------------- #
def frac_year(race_date: np.ndarray) -> np.ndarray:
    """Fractional year from a datetime64 array (day resolution)."""
    d = race_date.astype("datetime64[D]")
    y0 = d.astype("datetime64[Y]")
    yr = y0.astype(int) + 1970
    day_in = (d - y0).astype("timedelta64[D]").astype(float)
    span = ((y0 + 1).astype("datetime64[D]") - y0.astype("datetime64[D]")).astype(float)
    return yr + day_in / span


def beta0_gauge(v: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, float]:
    """G1 beta=0: remove the secular linear year trend. Returns (v_gauged, c_beta).

    t is centered before the shift so the gauged v keeps its raw level (~0); the
    additive offset is gauge-arbitrary and irrelevant to any correlation/R2.
    """
    tm = t - t.mean()
    c_beta = float((tm @ (v - v.mean())) / (tm @ tm))   # slope(v ~ t)
    return v - c_beta * tm, c_beta


# --------------------------------------------------------------------------- #
# locating + loading a registered fit for a slice                              #
# --------------------------------------------------------------------------- #
def _model_cfg_stem(model: str, nu: float):
    """(ModelConfig, registry stem) for one of the four nested models.

    Each model is registered by a different fitting stage under its own stem:
    ``baseline_{nutag}_best`` / ``agingS4gv_{nutag}_best`` / ``drift_{nutag}_best``
    / ``full_{nutag}_best``. ``aging`` is the no-d (rank-1 + aging block) fit.
    """
    if model == "baseline":
        cfg = baseline_cfg(nu)
        return cfg, registry.model_stem(cfg, "best")
    if model == "aging":
        cfg = aging_cfg(nu, use_d=False)
        return cfg, aging_stem(cfg)
    if model in ("drift", "full"):
        cfg = drift_cfg(nu, variant=model)
        return cfg, drift_stem(model, nu)
    raise ValueError(f"model must be one of {MODELS}, got {model!r}")


def fit_dir(spec, model: str = MODEL, nu: float = NU) -> Path:
    """Deterministic registry path of `model`'s production fit for `spec`."""
    cfg, stem = _model_cfg_stem(model, nu)
    parent = MODELS_ROOT / registry.slice_slug(spec)
    return registry.fit_path(parent, stem, spec, cfg, resample_tag="base")


def full_fit_dir(spec, nu: float = NU) -> Path:
    """Back-compat shim: the production full-model fit dir (``fit_dir(.., 'full')``)."""
    return fit_dir(spec, MODEL, nu)


def load_slice_vj(slice_name: str, *, model: str = MODEL,
                  nu: float = NU, mrc: int = MRC) -> pd.DataFrame:
    """One row per race in `slice_name` with raw + beta=0-gauged `model` v_j.

    Columns: race_id, year_frac, n_{slice}, v_{slice}_raw, v_{slice} (gauged).
    `model` defaults to the production full fit (used to build the merged table);
    q08 passes the other models to compare them on a fixed slice.
    """
    spec = S.build_spec(slice_name, min_race_count=mrc)
    fdir = fit_dir(spec, model, nu)
    if not (fdir / "fit.pkl").is_file():
        raise SystemExit(f"missing {slice_name} {model} fit: {fdir}")
    fd = load_slice(spec)
    m = registry.load_fit(fdir, fd)
    assert np.array_equal(m.data.race_ids, fd.race_ids)
    t = frac_year(fd.race_date)
    v_raw = np.asarray(m.params["v"], dtype=float)
    v_g, c_beta = beta0_gauge(v_raw, t)
    print(f"[fit] {slice_name:8s} {model:8s} {fdir.name}  J={fd.J:4d}  "
          f"c_beta={c_beta:+.5f} log/yr")
    return pd.DataFrame({
        "race_id": fd.race_ids,
        "year_frac": t,
        f"n_{slice_name}": np.bincount(fd.col_idx, minlength=fd.J),
        f"v_{slice_name}_raw": v_raw,
        f"v_{slice_name}": v_g,
    })


def load_vj_boot_wide(slice_name: str, *, model: str = MODEL, nu: float = NU,
                      mrc: int = MRC
                      ) -> tuple[Path, np.ndarray, np.ndarray, np.ndarray] | None:
    """Re-gauged bootstrap replicate matrix of the `model` v_j.

    Reads the fit's bootstrap replicates (``bootstrap/race_factors.parquet``)
    and re-applies the same G1 beta=0 gauge to **each** replicate, so the result
    is the gauged quantity that q03/q04 (and the weather-correction in q07)
    operate on. `model` defaults to the production full fit; q08 passes the other
    models. Returns None if no bootstrap is present.

    Returns ``(fit_dir, race_ids, year_frac, Bg)`` where ``Bg`` is the (R, J)
    replicate-by-race matrix, aligned column-wise to ``race_ids`` (the slice's
    ``load_slice`` order). Replicate columns absent from the bootstrap table
    come back as NaN.
    """
    spec = S.build_spec(slice_name, min_race_count=mrc)
    fdir = fit_dir(spec, model, nu)
    bp = fdir / "bootstrap" / "race_factors.parquet"
    if not bp.is_file():
        return None
    fd = load_slice(spec)
    t = frac_year(fd.race_date)
    rf = pd.read_parquet(bp, columns=["run_id", "race_id", "v"])
    wide = (rf[rf["run_id"] > 0]
            .pivot(index="run_id", columns="race_id", values="v")
            .reindex(columns=fd.race_ids))
    B = wide.to_numpy(np.float64)                 # (R, J) replicate x race
    tm = t - t.mean()
    denom = float(tm @ tm)
    c = (B @ tm) / denom                          # per-replicate beta=0 slope (R,)
    Bg = B - np.outer(c, tm)                      # re-gauged replicates
    return fdir, np.asarray(fd.race_ids), t, Bg


def load_vj_boot_sd(slice_name: str, *, model: str = MODEL, nu: float = NU,
                    mrc: int = MRC) -> pd.DataFrame | None:
    """Per-race athlete-sampling SD of the beta=0-gauged `model` v_j.

    Thin wrapper over :func:`load_vj_boot_wide` (re-gauges each replicate, then
    takes the across-replicate SD per race). Returns None if no bootstrap is
    present. Columns: race_id, vsd_{slice} (log-time units, like v_j).
    """
    res = load_vj_boot_wide(slice_name, model=model, nu=nu, mrc=mrc)
    if res is None:
        return None
    fdir, race_ids, _t, Bg = res
    sd = np.nanstd(Bg, axis=0, ddof=1)
    R = int(np.isfinite(Bg).all(axis=1).sum())
    print(f"[boot] {slice_name:8s} {fdir.name}  R={R}  "
          f"vsd median={np.nanmedian(sd):.5f}")
    return pd.DataFrame({"race_id": race_ids, f"vsd_{slice_name}": sd})


def build_vj_wide(*, nu: float = NU, mrc: int = MRC) -> pd.DataFrame:
    """Wide per-race table: full-model v_j for each of the six slices.

    Keyed on the global race_id; NaN where a race is absent from a slice. The
    ALL_B slice is the race-set superset so its `year_frac` covers every race.
    """
    base: pd.DataFrame | None = None
    for s in SLICE_ORDER:
        d = load_slice_vj(s, nu=nu, mrc=mrc)
        if base is None:
            base = d
        else:
            base = base.merge(d.drop(columns="year_frac"), on="race_id", how="outer")
    assert base is not None
    return base


# --------------------------------------------------------------------------- #
# covariate joins                                                              #
# --------------------------------------------------------------------------- #
def load_year_elevation() -> pd.DataFrame:
    """Per (series_key, year) elevation; drop rows with no usable GPX."""
    e = pd.read_csv(COURSE_PROFILE_DIR / "marathon_year_elevation.csv")
    e = e[e["net_gain_m"].notna()].copy()       # drops status=="no_data"
    out = e[["series_key", "year", "course_id", "net_gain_m", "total_gain_m",
             "uncertainty_total_m", "uncertainty_net_m", "status"]].rename(
        columns={"uncertainty_total_m": "total_gain_sd_m",
                 "uncertainty_net_m": "net_gain_sd_m", "status": "elev_status"})
    out["abs_net_gain_m"] = out["net_gain_m"].abs()
    dup = out.duplicated(["series_key", "year"]).sum()
    if dup:
        print(f"[WARN] {dup} duplicate (series_key, year) elevation rows; keeping first.")
        out = out.drop_duplicates(["series_key", "year"])
    return out


def simple_metric_cols(ss: pd.DataFrame) -> list[str]:
    """All naive 'simple' race metrics in race_stats_summary: *_sec and *_pct.

    Field-contaminated proxies for race difficulty (finish-time percentiles,
    top-N elites, sub-N / BQ shares). Counts and ids are excluded.
    """
    return [c for c in ss.columns if c.endswith("_sec") or c.endswith("_pct")]


def build_merged(*, nu: float = NU, mrc: int = MRC) -> tuple[pd.DataFrame, list[str]]:
    """Full per-race table: six slices' v_j + weather + elevation + simple metrics."""
    df = build_vj_wide(nu=nu, mrc=mrc)

    comp = pd.read_parquet(DATA_DIR / "competitions.parquet")
    df = df.merge(comp[["race_id", "competition_id", "series_key", "country", "year"]],
                  on="race_id", how="left")

    # weather (per competition_id)
    wx = pd.read_csv(WEATHER_DIR / "marathon_weather_features.csv")
    wx_cols = ["competition_id", "wbgt_field", "wbgt_max", "temp_field", "wind_field",
               "precip_total_mm", "wbgt_rising_slope_c_per_h", "wbgt_spatial_spread",
               "wbgt_anchor_sens"]
    df = df.merge(wx[wx_cols], on="competition_id", how="left")

    # elevation (per series_key, year)
    df = df.merge(load_year_elevation(), on=["series_key", "year"], how="left")

    # simple metrics (per competition_id == comp_id)
    ss = pd.read_csv(MISC_DIR / "race_stats_summary.csv")
    metric_cols = simple_metric_cols(ss)
    df = df.merge(ss[["comp_id"] + metric_cols].rename(columns={"comp_id": "competition_id"}),
                  on="competition_id", how="left")
    return df, metric_cols


# --------------------------------------------------------------------------- #
# small stats utilities shared by q02 / q03                                    #
# --------------------------------------------------------------------------- #
def standardize(a: np.ndarray) -> np.ndarray:
    return (a - np.nanmean(a)) / np.nanstd(a)


def physical_course_id(df: pd.DataFrame) -> pd.Series:
    """Course-change-aware cluster id: series_key + course_id.

    The elevation join carries `course_id` (a series can re-route between years);
    clustering on (series_key, course_id) groups editions run on the *same
    physical course*, which is the level at which `total_gain_m` actually varies.
    """
    cid = pd.to_numeric(df["course_id"], errors="coerce").astype("Int64").astype(str)
    return df["series_key"].astype(str) + "#" + cid


def allb_design(df: pd.DataFrame, weather: str = "temp_field",
                course: str = "total_gain_m", target: str = "v_ALL_B"
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """Standardized (y, X[1,temp,gain], cluster_id, masked rows) for ALL_B.

    Complete cases on target + both predictors. Predictors and response are
    z-standardized so betas match q03/q04. Returns the physical-course cluster
    id aligned to the masked rows.
    """
    m = (df[target].notna() & df[weather].notna() & df[course].notna()).to_numpy()
    d = df.loc[m].copy()
    y = standardize(d[target].to_numpy(float))
    Xz = np.column_stack([standardize(d[weather].to_numpy(float)),
                          standardize(d[course].to_numpy(float))])
    X = np.column_stack([np.ones(len(d)), Xz])
    cl = physical_course_id(d).to_numpy()
    return y, X, cl, d


def ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """(beta, se, R2) for y = X beta; X already includes the intercept column."""
    XtX = X.T @ X
    beta = np.linalg.solve(XtX, X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    sigma2 = (resid @ resid) / (n - k)
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(XtX)))
    r2 = 1 - (resid @ resid) / ((y - y.mean()) ** 2).sum()
    return beta, se, r2


def robust_wls(X: np.ndarray, y: np.ndarray, prior_w: np.ndarray | None = None,
               *, huber: bool = False, c: float = 1.345,
               max_iter: int = 100, tol: float = 1e-9
               ) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """One IRLS routine covering all four regression variants.

    total weight = (prior weight) x (Huber weight on the current residual):
      * OLS         : prior_w=None, huber=False  (exactly `ols` above)
      * WLS         : prior_w=1/var, huber=False  (down-weight uncertain v_j)
      * Robust      : prior_w=None, huber=True    (Huber M-estimator, c=1.345)
      * Robust+WLS  : prior_w=1/var, huber=True

    The Huber scale is a per-iteration MAD of the (prior-weighted) residual.
    Returns (beta, se, weighted_R2, final_weights). `se`/`R2` are weighted with
    the final total weights; `final_weights` exposes the Huber down-weighting so
    callers can flag outliers. X includes the intercept column.
    """
    n, k = X.shape
    pw = np.ones(n) if prior_w is None else np.asarray(prior_w, float).copy()
    W = pw.copy()
    beta = np.linalg.solve((X * W[:, None]).T @ X, (X * W[:, None]).T @ y)
    for _ in range(max_iter):
        r = y - X @ beta
        if huber:
            s = 1.4826 * np.median(np.abs(r - np.median(r)))
            s = max(s, 1e-12)
            u = np.abs(r) / s
            wh = np.where(u <= c, 1.0, c / np.maximum(u, 1e-12))
        else:
            wh = np.ones(n)
        W = pw * wh
        XtW = X * W[:, None]
        beta_new = np.linalg.solve(XtW.T @ X, XtW.T @ y)
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    r = y - X @ beta
    XtWX = (X * W[:, None]).T @ X
    sigma2 = float((W * r * r).sum() / (n - k))
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(XtWX)))
    wmean = float((W * y).sum() / W.sum())
    ss_res = float((W * r * r).sum())
    ss_tot = float((W * (y - wmean) ** 2).sum())
    r2 = 1 - ss_res / ss_tot
    return beta, se, r2, W
