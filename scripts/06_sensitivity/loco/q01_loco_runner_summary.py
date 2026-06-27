"""q01 -- LOCO (runner side) gauge-aware summary: signal vs gauge vs noise.

Reads the fits saved by e01 (``baseline/`` + ``runner_{c}/``) and, for each
excluded country c, decomposes the race-factor movement on the common race set

    dv_j = v_j^(-c) - v_j^(baseline)        (+ = race looks harder without c)

into the part a cross-fit GAUGE move can explain and the part it cannot, then
compares both against the bootstrap sampling SD of v_j -- answering "is the
README's |dv| statistically meaningful, gauge convention, or noise?".

WHY dv CARRIES A GAUGE COMPONENT
--------------------------------
Both fits impose mean(v)=0 and the beta=0 APC gauge (slope(v ~ t) = 0) IN-FIT
(``apc_gauge_beta0``), each over its OWN race set. So when no races are lost
the level and linear-in-date parts of dv vanish by construction (kept here as
sanity columns); when races drop they pick up a small, exactly-conventional
set-difference re-gauge. The one genuinely uncontrolled cross-fit direction
that touches v is the date-QUADRATIC, Gq (docs/model_derivation.md S7.4-S7.7):

    v_j -> v_j - eps*t_j^2     paid by    f(A) -> f(A) + eps*A^2,
    d_i -> d_i + 2*eps*debut_i (debut-correlated drift),  u_i -> (constants)

Gq is pinned per fit by the min-norm d_i prior (d orthogonal to
{1, debut, entry_age} over that fit's eligible set) and the pin RELOCATES when
the athlete population changes -- which is exactly what a LOCO refit does. The
quadratic-in-date component of dv is therefore convention, not country signal.

Decomposition (OLS over common races, mirroring race_comparison/q01):

    dv_j = affine_j + c2*q_j + r_j,    q = t~^2 orthogonalized against {1, t~}

r_j (``dv_resid``) is the gauge-free movement, the only part attributable to
the country itself.

THE THREE NUMBERS COMPARED (all also in sec at a 3:00 marathon, x10800)
-----------------------------------------------------------------------
1. movement       med|dv| raw vs med|dv_resid|: how much of the headline |dv|
                  was smooth-in-date (removable) vs real.
2. gauge band     two estimates of the Gq uncertainty:
                    realized  |c2*q_j|              the quadratic the refit
                                                    actually picked up;
                    priced    delta_plaus*|q_j|     the R=1 plausible re-pin,
                              delta_plaus = omega_d / (2*SD(debut)) on the
                              baseline eligible set -- a quad gauge motion
                              whose drift cost equals the WHOLE fitted drift
                              spread (same pricing as race_comparison/q01 and
                              aging_trend/q02).
                  R_kill = 2*|c2|*SD(debut)/omega_d prices the FITTED c2:
                  << 1 means the observed quadratic is cheap to be pure gauge
                  (do not read it as signal); >> 1 means it cannot be gauge.
                  Aging cross-check: a pure Gq relocation by eps must surface
                  as +eps*A^2 in the aging curve, so eps_f (quad coefficient
                  of f_loco - f_base on a grid) should match -c2 if and only
                  if the quadratic in dv is the gauge re-pin
                  (``gauge_consistency`` = eps_f / -c2).
3. bootstrap SD   per-race SD of v over the SOURCE fit's athlete-weight
                  bootstrap, with every replicate residualized on
                  {1, t~, t~^2} exactly like dv (apples-to-apples; raw SD also
                  reported). NOTE the bootstrap is a smaller perturbation than
                  LOCO (reweights, never removes), so ratio >> 1 means "real
                  athlete-set sensitivity beyond sampling noise", not a bug.

HOME / AWAY (the README's headline finding, recomputed gauge-free)
------------------------------------------------------------------
home = common races with race_country == c, on raw dv and on dv_resid, with

    sd_boot_home = SD over replicates of the home-race mean of the
                   residualized bootstrap v (sampling noise OF the home-mean
                   statistic; the away mean is its mechanical complement under
                   mean(dv) ~= 0 and is not separately tested),
    gauge_home   = delta_plaus * |mean_home(q)|   (how far the priced re-pin
                   can move the home mean).

Verdict per country: NOISE (|z_home| <= 2) / SIGNAL (clears noise AND the
priced band) / SIGNAL* (clears noise but sits inside the priced worst-case
band) / NO-HOME / POINT-ONLY. For LOCO the starred caveat is conservative:
the REALIZED re-pin between this specific pair of fits is directly measured
by c2 (R_kill prices it), and when R_kill << 1 the actual gauge motion was
negligible -- delta_plaus is the cross-fit worst case, not an estimate.

PRACTICAL WARNING (why the affine removal is load-bearing, not pedantic):
saved fits can sit OFF the beta=0 manifold. The in-fit gauge is applied per
sweep, but Anderson mixing blends history states including the (ungauged)
warm start, and nothing re-applies the gauge post-fit -- so a warm-started
fit that converges in a few iterations keeps an arbitrary fraction of its
init's date tilt (the production ALL_B mrc2 source fit and the e01 baseline
sit at slope(v~t) = +0.00165/yr; the long LOCO refits at exactly 0). The
study header prints each baseline's measured slope as a diagnostic; any
common tilt it injects into dv lands in the affine bucket and is removed.

NO REFIT -- reads saved fits + the source fit's bootstrap only.

Outputs under each study dir, in ``summary/``:
    q01_summary.parquet/.csv   one row per excluded country
    q01_per_race.parquet       per (country, common race): dv decomposition,
                               gauge band, boot SDs, ratio, is_home
    q01_summary.md             readable table
    meta.json                  provenance + study-level gauge price
plus a stdout table.

Run::

    python scripts/06_sensitivity/loco/q01_loco_runner_summary.py        # all studies
    python scripts/06_sensitivity/loco/q01_loco_runner_summary.py \
        --study results/sensitivity/loco/ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from marathon_decomp import load_slice  # noqa: E402
from marathon_decomp.aging import aging_curve_from_payload  # noqa: E402
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "sensitivity" / "loco"
SEC3H = 180.0 * 60.0          # log-time -> seconds at a 3:00:00 marathon


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def _frac_years(dates) -> np.ndarray:
    rd = pd.DatetimeIndex(pd.to_datetime(dates))
    return (rd.year + (rd.dayofyear - 1) / 365.25).to_numpy(np.float64)


def load_bundle(fit_dir: Path) -> dict:
    """Saved fit dir -> payload + slice-derived per-race arrays (no model build)."""
    with open(fit_dir / "fit.pkl", "rb") as f:
        payload = pickle.load(f)
    fd = load_slice(payload["spec"], payload["data_version"])
    man = json.loads((fit_dir / "manifest.json").read_text())
    return dict(
        fit_dir=fit_dir, payload=payload, fd=fd, manifest=man,
        v=np.asarray(payload["params"]["v"], np.float64),
        race_ids=np.asarray(fd.race_ids, np.int64),
        t=_frac_years(fd.race_date),
        race_country=np.asarray(fd.race_country, dtype=object),
    )


def drift_price(bundle: dict) -> dict:
    """omega_d + SD(centered debut) over the d-eligible set of one fit.

    Replicates the Model constructor's eligibility (n_i >= d_min_n and career
    span >= d_min_span_years) so no model instantiation is needed. The Gq
    re-pin moves d_i by 2*eps*debut_i, so the drift cost of a quad gauge
    motion eps is R = 2*|eps|*SD(debut)/omega_d (race_comparison/q01 logic).
    """
    fd, cfg = bundle["fd"], bundle["payload"]["config"]
    row = np.asarray(fd.row_idx, np.int64)
    An = np.asarray(fd.A_n, np.float64)
    n_per = np.bincount(row, minlength=fd.I)
    An_min = np.full(fd.I, np.inf)
    An_max = np.full(fd.I, -np.inf)
    np.minimum.at(An_min, row, An)
    np.maximum.at(An_max, row, An)
    span = np.where(np.isfinite(An_max) & np.isfinite(An_min), An_max - An_min, 0.0)
    elig = ((n_per >= cfg.d_min_n) & (span >= cfg.d_min_span_years)
            if cfg.use_d else np.zeros(fd.I, bool))

    t_cell = bundle["t"][np.asarray(fd.col_idx, np.int64)]
    debut = np.full(fd.I, np.nan)
    debut[row] = t_cell - An                       # constant per athlete
    f_e = debut[elig & np.isfinite(debut)]
    sd_debut = float(np.std(f_e - f_e.mean(), ddof=0)) if f_e.size else 0.0
    omega_d = float(np.sqrt(max(float(bundle["payload"]["params"]["omega_d2"]), 0.0)))
    return dict(omega_d=omega_d, sd_debut=sd_debut, n_elig=int(elig.sum()))


def load_boot_wide(source_fit_dir: Path, race_ids: np.ndarray) -> np.ndarray | None:
    """Bootstrap replicates of v aligned to `race_ids` -> (R, n) or None.
    Races absent from the bootstrap table become NaN columns (NaN boot SD)."""
    p = source_fit_dir / "bootstrap" / "race_factors.parquet"
    if not p.is_file():
        return None
    rf = pd.read_parquet(p, columns=["run_id", "race_id", "v"])
    wide = (rf[rf["run_id"] > 0]
            .pivot(index="run_id", columns="race_id", values="v")
            .reindex(columns=race_ids))
    if wide.empty:
        return None
    return wide.to_numpy(np.float64)


# ---------------------------------------------------------------------------
# gauge decomposition
# ---------------------------------------------------------------------------

def date_design(t: np.ndarray) -> dict:
    """{1, t~, t~^2} design on the common races + the orthogonalized quad shape."""
    t_c = t - t.mean()
    X = np.column_stack([np.ones_like(t_c), t_c, t_c ** 2])
    Xp = np.linalg.pinv(X)                          # (3, n)
    X_aff = X[:, :2]
    coef_q, *_ = np.linalg.lstsq(X_aff, t_c ** 2, rcond=None)
    q = t_c ** 2 - X_aff @ coef_q                   # quad shape, orthogonal to {1, t~}
    return dict(t_c=t_c, X=X, Xp=Xp, q=q)


def decompose(dv: np.ndarray, des: dict) -> dict:
    """dv = affine + c2*q + resid (exact); c2 is the raw t~^2 coefficient."""
    beta = des["Xp"] @ dv
    c2 = float(beta[2])
    smooth = des["X"] @ beta
    dv_quad = c2 * des["q"]
    return dict(c2=c2, level=float(beta[0]), slope=float(beta[1]),
                dv_affine=smooth - dv_quad, dv_quad=dv_quad, dv_resid=dv - smooth)


def resid_rows(B: np.ndarray, des: dict) -> np.ndarray:
    """Residualize each replicate row of (R, n) on {1, t~, t~^2} (NaN-safe:
    rows are projected using finite columns only via masking the coefficients
    is unnecessary here -- bootstrap NaNs come as whole-race columns, which
    drop out of the lstsq fit when zero-filled with a mask)."""
    if not np.isnan(B).any():
        coef = B @ des["Xp"].T                      # (R, 3)
        return B - coef @ des["X"].T
    out = np.full_like(B, np.nan)
    fin_cols = ~np.isnan(B).any(axis=0)
    sub = date_design(des["t_c"][fin_cols] + 0.0)   # re-center on finite cols
    Bf = B[:, fin_cols]
    coef = Bf @ sub["Xp"].T
    out[:, fin_cols] = Bf - coef @ sub["X"].T
    return out


def aging_quad_eps(p_base: dict, p_loco: dict, a_max: float,
                   n_grid: int = 201) -> dict:
    """Quad coefficient of (f_loco - f_base) on [0, a_max] + identified remainder.

    Both curves anchor f(0)=0, so the gauge orbit of the difference is
    span{A, A^2}; eps_f is the A^2 coefficient (the Gq companion), and the
    remainder after removing span{A, A^2} is the IDENTIFIED aging change.
    """
    A = np.linspace(0.0, a_max, n_grid)
    df = aging_curve_from_payload(p_loco, A) - aging_curve_from_payload(p_base, A)
    Xa = np.column_stack([A, A * A])
    coef, *_ = np.linalg.lstsq(Xa, df, rcond=None)
    ident = df - Xa @ coef
    return dict(eps_f=float(coef[1]), d_aging_max=float(np.abs(df).max()),
                d_aging_ident_max=float(np.abs(ident).max()))


# ---------------------------------------------------------------------------
# per-country analysis
# ---------------------------------------------------------------------------

def analyse_country(c: str, base: dict, loco: dict, price: dict,
                    boot_full: np.ndarray | None) -> tuple[dict, pd.DataFrame]:
    common, iB, iL = np.intersect1d(base["race_ids"], loco["race_ids"],
                                    return_indices=True)
    dv = loco["v"][iL] - base["v"][iB]
    t = base["t"][iB]
    country = base["race_country"][iB]
    is_home = (country == c)

    des = date_design(t)
    dec = decompose(dv, des)
    var_dv = float(np.var(dv))
    shares = {k: (float(np.var(dec[f"dv_{k}"]) / var_dv) if var_dv > 1e-18 else np.nan)
              for k in ("affine", "quad", "resid")}

    # ---- gauge band: realized + priced --------------------------------
    gauge_real = np.abs(dec["c2"] * des["q"])
    delta_plaus = (price["omega_d"] / (2.0 * price["sd_debut"])
                   if price["sd_debut"] > 0 else np.inf)
    gauge_plaus = delta_plaus * np.abs(des["q"])
    R_kill = (2.0 * abs(dec["c2"]) * price["sd_debut"] / price["omega_d"]
              if price["omega_d"] > 0 else np.inf)

    # ---- aging cross-check --------------------------------------------
    ag = aging_quad_eps(base["payload"], loco["payload"], a_max=ARGS.a_max)
    g_b = np.atleast_1d(np.asarray(base["payload"]["params"]["gamma"], float))
    g_l = np.atleast_1d(np.asarray(loco["payload"]["params"]["gamma"], float))
    max_dgamma = float(np.abs(g_l - g_b).max()) if g_b.shape == g_l.shape else np.nan
    consistency = (ag["eps_f"] / -dec["c2"]) if abs(dec["c2"]) > 1e-12 else np.nan

    # ---- bootstrap yardstick (source fit, same races, same projection) -
    sd_boot_raw = np.full(common.size, np.nan)
    sd_boot_res = np.full(common.size, np.nan)
    sd_boot_home = np.nan
    if boot_full is not None:
        # boot_full is aligned to base["race_ids"]; subset to the common races
        B = boot_full[:, iB]
        sd_boot_raw = np.nanstd(B, axis=0, ddof=1)
        Br = resid_rows(B, des)
        sd_boot_res = np.nanstd(Br, axis=0, ddof=1)
        if is_home.any():
            hm = np.nanmean(Br[:, is_home], axis=1)
            sd_boot_home = float(np.nanstd(hm, ddof=1))

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.abs(dec["dv_resid"]) / sd_boot_res
    fin = np.isfinite(ratio)

    # ---- home / away ----------------------------------------------------
    n_home = int(is_home.sum())
    home_raw = float(dv[is_home].mean()) if n_home else np.nan
    home_res = float(dec["dv_resid"][is_home].mean()) if n_home else np.nan
    away_res = float(dec["dv_resid"][~is_home].mean()) if (~is_home).any() else np.nan
    gauge_home_plaus = (float(delta_plaus * abs(des["q"][is_home].mean()))
                        if n_home else np.nan)
    gauge_home_real = (float((dec["dv_affine"] + dec["dv_quad"])[is_home].mean())
                       if n_home else np.nan)
    z_home = (home_res / sd_boot_home
              if n_home and np.isfinite(sd_boot_home) and sd_boot_home > 0 else np.nan)

    if n_home == 0:
        verdict = "NO-HOME"
    elif not np.isfinite(z_home):
        verdict = "POINT-ONLY"
    elif abs(z_home) <= 2.0:
        verdict = "NOISE"
    elif abs(home_res) <= gauge_home_plaus:
        verdict = "SIGNAL*"          # clears noise; inside priced worst-case band
    else:
        verdict = "SIGNAL"

    # ---- rank stability (raw + era-relative) ----------------------------
    sp_raw = float(spearmanr(base["v"][iB], loco["v"][iL]).statistic)
    vb_era = base["v"][iB] - des["X"] @ (des["Xp"] @ base["v"][iB])
    vl_era = loco["v"][iL] - des["X"] @ (des["Xp"] @ loco["v"][iL])
    sp_era = float(spearmanr(vb_era, vl_era).statistic)

    man = loco["manifest"]
    rec = dict(
        excl_country=c,
        n_excluded_runners=int(man.get("n_excluded_runners", 0)),
        n_dropped_races=int(man.get("n_dropped_races", 0)),
        n_common=int(common.size), n_home=n_home,
        # sanity (should be ~0 when no races dropped: both gauges in-fit)
        dv_mean=float(dv.mean()), dv_slope=dec["slope"],
        # movement, raw vs gauge-free
        med_abs_dv_s=float(np.median(np.abs(dv))) * SEC3H,
        med_abs_resid_s=float(np.median(np.abs(dec["dv_resid"]))) * SEC3H,
        share_affine=shares["affine"], share_quad=shares["quad"],
        share_resid=shares["resid"],
        # gauge band
        c2=dec["c2"], R_kill=R_kill,
        gauge_real_med_s=float(np.median(gauge_real)) * SEC3H,
        gauge_real_max_s=float(gauge_real.max()) * SEC3H,
        delta_plaus=delta_plaus,
        gauge_plaus_med_s=float(np.median(gauge_plaus)) * SEC3H,
        gauge_plaus_max_s=float(gauge_plaus.max()) * SEC3H,
        # aging cross-check
        eps_f=ag["eps_f"], gauge_consistency=consistency,
        d_aging_max_s=ag["d_aging_max"] * SEC3H,
        d_aging_ident_max_s=ag["d_aging_ident_max"] * SEC3H,
        max_abs_dgamma=max_dgamma,
        # bootstrap comparison
        boot_raw_med_s=float(np.nanmedian(sd_boot_raw)) * SEC3H,
        boot_resid_med_s=float(np.nanmedian(sd_boot_res)) * SEC3H,
        ratio_med=float(np.median(ratio[fin])) if fin.any() else np.nan,
        ratio_p90=float(np.percentile(ratio[fin], 90)) if fin.any() else np.nan,
        frac_gt2=float((ratio[fin] > 2.0).mean()) if fin.any() else np.nan,
        # rank stability
        spearman_raw=sp_raw, spearman_era=sp_era,
        # home / away
        home_raw_s=home_raw * SEC3H if n_home else np.nan,
        home_resid_s=home_res * SEC3H if n_home else np.nan,
        away_resid_s=away_res * SEC3H if np.isfinite(away_res) else np.nan,
        sd_boot_home_s=sd_boot_home * SEC3H if np.isfinite(sd_boot_home) else np.nan,
        z_home=z_home,
        gauge_home_plaus_s=gauge_home_plaus * SEC3H if n_home else np.nan,
        gauge_home_real_s=gauge_home_real * SEC3H if n_home else np.nan,
        verdict=verdict,
    )

    per_race = pd.DataFrame(dict(
        excl_country=c, race_id=common,
        series_key=np.asarray(base["fd"].race_series, dtype=object)[iB],
        country=country, t=t, is_home=is_home,
        dv=dv, dv_affine=dec["dv_affine"], dv_quad=dec["dv_quad"],
        dv_resid=dec["dv_resid"], q_shape=des["q"],
        gauge_plaus=gauge_plaus, sd_boot_raw=sd_boot_raw,
        sd_boot_resid=sd_boot_res, ratio_resid_boot=ratio,
    ))
    return rec, per_race


# ---------------------------------------------------------------------------
# one study
# ---------------------------------------------------------------------------

def run_study(study: Path, countries_filter: list[str] | None) -> pd.DataFrame | None:
    meta = json.loads((study / "loco_meta.json").read_text())
    cdirs = sorted(d for d in study.glob("runner_*") if (d / "manifest.json").is_file())
    if countries_filter:
        cdirs = [d for d in cdirs if d.name[len("runner_"):] in countries_filter]
    if not (study / "baseline" / "manifest.json").is_file() or not cdirs:
        print(f"  [skip] {study}: missing baseline/ or runner_* fits")
        return None

    print(f"\n=== q01 LOCO summary: {study.parent.name}/{study.name} "
          f"({len(cdirs)} countries) ===")
    base = load_bundle(study / "baseline")
    price = drift_price(base)
    delta_plaus = (price["omega_d"] / (2.0 * price["sd_debut"])
                   if price["sd_debut"] > 0 else np.inf)
    print(f"    baseline J={base['race_ids'].size}  omega_d={price['omega_d']:.4e} "
          f"SD(debut)={price['sd_debut']:.2f} yr (n_elig={price['n_elig']:,}) "
          f"-> delta_plaus={delta_plaus:.3e} /yr^2")

    # Diagnostic: is the baseline actually ON the beta=0 manifold? (Warm-started
    # fits can keep their init's tilt -- see the module docstring.)
    t_c0 = base["t"] - base["t"].mean()
    base_slope = float(np.dot(base["v"], t_c0) / np.dot(t_c0, t_c0))
    if abs(base_slope) > 1e-6:
        print(f"    [warn] baseline is OFF the beta=0 manifold: slope(v~t) = "
              f"{base_slope:+.6f}/yr (~{abs(base_slope) * 6 * SEC3H:.0f} s at the "
              "window ends). This tilt enters every dv as a common affine "
              "artifact; the decomposition removes it.")

    src = Path(meta.get("source_fit_dir", ""))
    if not src.is_dir():
        src = (RESULTS_DIR / "models" / study.parent.name / study.name)
    boot_full = load_boot_wide(src, base["race_ids"])
    if boot_full is None:
        print(f"    [warn] no bootstrap under {src} -- boot columns will be NaN")
    else:
        print(f"    bootstrap: {boot_full.shape[0]} replicates from {src.name}")

    recs, frames = [], []
    for cdir in cdirs:
        c = cdir.name[len("runner_"):]
        loco = load_bundle(cdir)
        rec, pr = analyse_country(c, base, loco, price, boot_full)
        recs.append(rec)
        frames.append(pr)

    df = pd.DataFrame(recs).sort_values("n_excluded_runners",
                                        ascending=False).reset_index(drop=True)
    per_race = pd.concat(frames, ignore_index=True)

    out = study / "summary"
    out.mkdir(exist_ok=True)
    df.to_parquet(out / "q01_summary.parquet", index=False)
    num = df.select_dtypes(include="number").columns
    df_csv = df.copy()
    df_csv[num] = df_csv[num].round(6)
    df_csv.to_csv(out / "q01_summary.csv", index=False)
    per_race.to_parquet(out / "q01_per_race.parquet", index=False)
    (out / "meta.json").write_text(json.dumps(dict(
        study=str(study), source_fit_dir=str(src),
        n_boot_replicates=int(boot_full.shape[0]) if boot_full is not None else 0,
        omega_d=price["omega_d"], sd_debut=price["sd_debut"],
        n_elig=price["n_elig"], delta_plaus=delta_plaus, a_max=ARGS.a_max,
        baseline_slope_v_t=base_slope,
        updated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    ), indent=2))

    # ---- markdown -------------------------------------------------------
    cols = ["country", "n_drop", "n_home", "med|dv| s", "resid med s",
            "resid share", "gauge plaus med s", "R_kill", "boot med s",
            "ratio med", "home raw s", "home resid s", "z_home",
            "gauge home s", "verdict"]
    lines = [f"# q01 LOCO gauge-aware summary -- {study.parent.name}/{study.name}",
             "",
             f"delta_plaus = {delta_plaus:.3e} /yr^2 (R=1 Gq re-pin; omega_d = "
             f"{price['omega_d']:.3e}, SD(debut) = {price['sd_debut']:.2f} yr). "
             f"Baseline slope(v~t) = {base_slope:+.6f}/yr"
             + (" -- OFF the beta=0 manifold; its tilt is removed as part of the "
                "affine component." if abs(base_slope) > 1e-6 else ".")
             + " All seconds at a 3:00:00 marathon. dv = v(-country) - "
             "v(baseline); resid = dv after projecting out {1, t, t^2} (the "
             "gauge-movable directions); ratio = per-race |resid| / "
             "residualized bootstrap SD. Verdicts: SIGNAL clears 2x noise and "
             "the priced worst-case gauge band; SIGNAL* clears noise but sits "
             "inside that (conservative) band -- check R_kill: << 1 means the "
             "realized re-pin between these two fits was negligible.",
             "",
             "| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        def _f(x, fmt="{:.0f}"):
            return fmt.format(x) if pd.notna(x) else "-"
        lines.append("| " + " | ".join([
            r["excl_country"], f"{r['n_excluded_runners']:,}", f"{r['n_home']}",
            _f(r["med_abs_dv_s"]), _f(r["med_abs_resid_s"]),
            _f(r["share_resid"], "{:.2f}"),
            _f(r["gauge_plaus_med_s"], "{:.1f}"), _f(r["R_kill"], "{:.2f}"),
            _f(r["boot_resid_med_s"], "{:.0f}"), _f(r["ratio_med"], "{:.2f}"),
            _f(r["home_raw_s"], "{:+.0f}"), _f(r["home_resid_s"], "{:+.0f}"),
            _f(r["z_home"], "{:+.1f}"), _f(r["gauge_home_plaus_s"], "{:.1f}"),
            r["verdict"],
        ]) + " |")
    (out / "q01_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"    wrote {out / 'q01_summary.parquet'} (+csv/md, per_race)")

    # ---- stdout ----------------------------------------------------------
    hdr = (f"{'country':>8} {'n_drop':>7} {'n_hm':>4} | {'med|dv|':>7} {'resid':>6} "
           f"{'r_sh':>5} | {'g_plaus':>7} {'R_kill':>6} | {'boot':>5} {'r_med':>5} | "
           f"{'home_raw':>8} {'home_res':>8} {'z':>5} {'g_home':>6}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for _, r in df.iterrows():
        def _n(x, fmt):
            return fmt.format(x) if pd.notna(x) else " " * (len(fmt.format(0)) - 1) + "-"
        print(f"{r['excl_country']:>8} {r['n_excluded_runners']:>7,} {r['n_home']:>4} | "
              f"{r['med_abs_dv_s']:>7.0f} {r['med_abs_resid_s']:>6.0f} "
              f"{r['share_resid']:>5.2f} | "
              f"{r['gauge_plaus_med_s']:>7.1f} {r['R_kill']:>6.2f} | "
              f"{_n(r['boot_resid_med_s'], '{:>5.0f}')} {_n(r['ratio_med'], '{:>5.2f}')} | "
              f"{_n(r['home_raw_s'], '{:>8.0f}')} {_n(r['home_resid_s'], '{:>8.0f}')} "
              f"{_n(r['z_home'], '{:>5.1f}')} {_n(r['gauge_home_plaus_s'], '{:>6.1f}')}  "
              f"{r['verdict']}")
    print()
    return df


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def discover_studies(root: Path) -> list[Path]:
    return sorted(p.parent for p in root.glob("*/*/loco_meta.json"))


ARGS: argparse.Namespace


def main() -> None:
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", nargs="+", default=None,
                    help="study dir(s) (default: every dir under "
                         "results/sensitivity/loco with a loco_meta.json).")
    ap.add_argument("--countries", nargs="+", default=None,
                    help="restrict to these excluded countries.")
    ap.add_argument("--a-max", type=float, default=10.0,
                    help="career-age horizon for the aging-curve cross-check.")
    ARGS = ap.parse_args()

    studies = ([Path(s) for s in ARGS.study] if ARGS.study
               else discover_studies(OUT_ROOT))
    if not studies:
        raise SystemExit(f"no LOCO studies found under {OUT_ROOT}")
    for s in studies:
        if not (s / "loco_meta.json").is_file():
            raise SystemExit(f"not a LOCO study dir (no loco_meta.json): {s}")
        run_study(s, ARGS.countries)


if __name__ == "__main__":
    main()
