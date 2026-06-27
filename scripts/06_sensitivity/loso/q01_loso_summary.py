"""q01 -- LOSO gauge-aware summary: spillover signal vs gauge vs noise.

Race-side twin of ``../loco/q01_loco_runner_summary.py``. Reads the fits saved
by e01 (``baseline/`` + ``series_{label}/``) and, for each removed series
(venue), decomposes the SURVIVING races' factor movement

    dv_j = v_j^(-venue) - v_j^(baseline)      (+ = race looks harder without it)

into the part a cross-fit GAUGE move can explain and the part it cannot, then
compares both against the bootstrap sampling SD of v_j. The probe inverts vs
LOCO: the removed races are gone, so the measured quantity is the SPILLOVER
onto the survivors -- and the "home" stratum becomes the venue's
SAME-COUNTRY surviving races (the regional re-leveling of the README).

GAUGE BACKGROUND (same as the LOCO q01; see its docstring for the full story)
------------------------------------------------------------------------------
Both fits impose mean(v)=0 + the beta=0 APC gauge in-fit, but each over its
OWN race set -- and removing a venue CHANGES the set, so unlike runner-side
LOCO the affine part of dv is not zero even in exact arithmetic (a small,
exactly-conventional set-difference re-gauge). On top of that, saved fits can
sit OFF the beta=0 manifold entirely (warm-start + Anderson mixing + no
post-fit re-gauge -- the production ALL_B baselines sit at slope(v~t) ~
+0.0016/yr while long refits land at 0), which injects a large common tilt
into every dv. Both effects live in span{1, t, t~^2} and are removed exactly
by the decomposition

    dv_j = affine_j + c2*q_j + r_j,    q = t~^2 orthogonalized against {1, t~}

r_j (``dv_resid``) is the gauge-free spillover. The quadratic freedom (Gq) is
priced against omega_d (delta_plaus, R_kill) and cross-checked against its
aging-curve companion (eps_f ~ -c2 iff the quadratic is a pure re-pin); the
gauge-free residual is compared per-race with the source fit's bootstrap SD
(replicates residualized identically).

SAME-COUNTRY SPILLOVER (the README's headline, recomputed gauge-free)
---------------------------------------------------------------------
samec = surviving common races with race_country == the venue's country.
Reported on raw dv and dv_resid, with

    sd_boot_samec = SD over replicates of the same-country mean of the
                    residualized bootstrap v (sampling noise OF that mean),
    gauge_samec   = delta_plaus * |mean_samec(q)|  (priced worst-case motion).

Verdict per venue: NOISE (|z| <= 2) / SIGNAL (clears noise AND the priced
band) / SIGNAL* (clears noise, inside the conservative priced band -- check
R_kill: << 1 means the realized re-pin between the two fits was negligible) /
NO-SAMEC (venue is its country's only race in the design) / POINT-ONLY.

NO REFIT -- reads saved fits + the source fit's bootstrap only.

Outputs under each study dir, in ``summary/``:
    q01_summary.parquet/.csv   one row per removed series
    q01_per_race.parquet       per (series, surviving race): dv decomposition,
                               gauge band, boot SDs, ratio, is_samec
    q01_summary.md             readable table
    meta.json                  provenance + study-level gauge price
plus a stdout table.

Run::

    python scripts/06_sensitivity/loso/q01_loso_summary.py          # all studies
    python scripts/06_sensitivity/loso/q01_loso_summary.py \
        --study results/sensitivity/loso/ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b \
        --series boston london valencia
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

OUT_ROOT = RESULTS_DIR / "sensitivity" / "loso"
SEC3H = 180.0 * 60.0          # log-time -> seconds at a 3:00:00 marathon


# ---------------------------------------------------------------------------
# loading (mirrors loco/q01)
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
        race_series=np.asarray(fd.race_series, dtype=object),
    )


def drift_price(bundle: dict) -> dict:
    """omega_d + SD(centered debut) over the d-eligible set (see loco/q01)."""
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
    debut[row] = t_cell - An
    f_e = debut[elig & np.isfinite(debut)]
    sd_debut = float(np.std(f_e - f_e.mean(), ddof=0)) if f_e.size else 0.0
    omega_d = float(np.sqrt(max(float(bundle["payload"]["params"]["omega_d2"]), 0.0)))
    return dict(omega_d=omega_d, sd_debut=sd_debut, n_elig=int(elig.sum()))


def load_boot_wide(source_fit_dir: Path, race_ids: np.ndarray) -> np.ndarray | None:
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
# gauge decomposition (identical machinery to loco/q01)
# ---------------------------------------------------------------------------

def date_design(t: np.ndarray) -> dict:
    t_c = t - t.mean()
    X = np.column_stack([np.ones_like(t_c), t_c, t_c ** 2])
    Xp = np.linalg.pinv(X)
    X_aff = X[:, :2]
    coef_q, *_ = np.linalg.lstsq(X_aff, t_c ** 2, rcond=None)
    q = t_c ** 2 - X_aff @ coef_q
    return dict(t_c=t_c, X=X, Xp=Xp, q=q)


def decompose(dv: np.ndarray, des: dict) -> dict:
    beta = des["Xp"] @ dv
    c2 = float(beta[2])
    smooth = des["X"] @ beta
    dv_quad = c2 * des["q"]
    return dict(c2=c2, level=float(beta[0]), slope=float(beta[1]),
                dv_affine=smooth - dv_quad, dv_quad=dv_quad, dv_resid=dv - smooth)


def resid_rows(B: np.ndarray, des: dict) -> np.ndarray:
    if not np.isnan(B).any():
        coef = B @ des["Xp"].T
        return B - coef @ des["X"].T
    out = np.full_like(B, np.nan)
    fin_cols = ~np.isnan(B).any(axis=0)
    sub = date_design(des["t_c"][fin_cols] + 0.0)
    Bf = B[:, fin_cols]
    coef = Bf @ sub["Xp"].T
    out[:, fin_cols] = Bf - coef @ sub["X"].T
    return out


def aging_quad_eps(p_base: dict, p_loco: dict, a_max: float,
                   n_grid: int = 201) -> dict:
    A = np.linspace(0.0, a_max, n_grid)
    df = aging_curve_from_payload(p_loco, A) - aging_curve_from_payload(p_base, A)
    Xa = np.column_stack([A, A * A])
    coef, *_ = np.linalg.lstsq(Xa, df, rcond=None)
    ident = df - Xa @ coef
    return dict(eps_f=float(coef[1]), d_aging_max=float(np.abs(df).max()),
                d_aging_ident_max=float(np.abs(ident).max()))


# ---------------------------------------------------------------------------
# per-series analysis
# ---------------------------------------------------------------------------

def venue_country(base: dict, series_key: str) -> str | None:
    """Country of the removed venue, read off the BASELINE race table."""
    mask = base["race_series"] == series_key
    if not mask.any():
        return None
    vals, counts = np.unique(base["race_country"][mask].astype(str),
                             return_counts=True)
    return str(vals[np.argmax(counts)])


def analyse_series(label: str, base: dict, loso: dict, price: dict,
                   boot_full: np.ndarray | None) -> tuple[dict, pd.DataFrame]:
    man = loso["manifest"]
    series_key = man.get("series_key", label)
    vcountry = venue_country(base, series_key)

    common, iB, iL = np.intersect1d(base["race_ids"], loso["race_ids"],
                                    return_indices=True)
    dv = loso["v"][iL] - base["v"][iB]
    t = base["t"][iB]
    country = base["race_country"][iB]
    is_samec = (country == vcountry) if vcountry is not None \
        else np.zeros(common.size, bool)

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

    # ---- aging cross-check ---------------------------------------------
    ag = aging_quad_eps(base["payload"], loso["payload"], a_max=ARGS.a_max)
    g_b = np.atleast_1d(np.asarray(base["payload"]["params"]["gamma"], float))
    g_l = np.atleast_1d(np.asarray(loso["payload"]["params"]["gamma"], float))
    max_dgamma = float(np.abs(g_l - g_b).max()) if g_b.shape == g_l.shape else np.nan
    consistency = (ag["eps_f"] / -dec["c2"]) if abs(dec["c2"]) > 1e-12 else np.nan

    # ---- bootstrap yardstick --------------------------------------------
    sd_boot_raw = np.full(common.size, np.nan)
    sd_boot_res = np.full(common.size, np.nan)
    sd_boot_samec = np.nan
    if boot_full is not None:
        B = boot_full[:, iB]
        sd_boot_raw = np.nanstd(B, axis=0, ddof=1)
        Br = resid_rows(B, des)
        sd_boot_res = np.nanstd(Br, axis=0, ddof=1)
        if is_samec.any():
            sm = np.nanmean(Br[:, is_samec], axis=1)
            sd_boot_samec = float(np.nanstd(sm, ddof=1))

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.abs(dec["dv_resid"]) / sd_boot_res
    fin = np.isfinite(ratio)

    # ---- same-country spillover -----------------------------------------
    n_samec = int(is_samec.sum())
    samec_raw = float(dv[is_samec].mean()) if n_samec else np.nan
    samec_res = float(dec["dv_resid"][is_samec].mean()) if n_samec else np.nan
    rest_res = float(dec["dv_resid"][~is_samec].mean()) if (~is_samec).any() else np.nan
    gauge_samec_plaus = (float(delta_plaus * abs(des["q"][is_samec].mean()))
                         if n_samec else np.nan)
    gauge_samec_real = (float((dec["dv_affine"] + dec["dv_quad"])[is_samec].mean())
                        if n_samec else np.nan)
    z_samec = (samec_res / sd_boot_samec
               if n_samec and np.isfinite(sd_boot_samec) and sd_boot_samec > 0
               else np.nan)

    if n_samec == 0:
        verdict = "NO-SAMEC"
    elif not np.isfinite(z_samec):
        verdict = "POINT-ONLY"
    elif abs(z_samec) <= 2.0:
        verdict = "NOISE"
    elif abs(samec_res) <= gauge_samec_plaus:
        verdict = "SIGNAL*"          # clears noise; inside priced worst-case band
    else:
        verdict = "SIGNAL"

    # ---- rank stability (raw + era-relative) ----------------------------
    sp_raw = float(spearmanr(base["v"][iB], loso["v"][iL]).statistic)
    vb_era = base["v"][iB] - des["X"] @ (des["Xp"] @ base["v"][iB])
    vl_era = loso["v"][iL] - des["X"] @ (des["Xp"] @ loso["v"][iL])
    sp_era = float(spearmanr(vb_era, vl_era).statistic)

    rec = dict(
        excl_series=label, series_key=series_key, venue_country=vcountry,
        n_editions=int(man.get("n_editions", 0)),
        n_dropped_races=int(man.get("n_dropped_races", 0)),
        n_dropped_runners=int(man.get("n_dropped_runners", 0)),
        n_common=int(common.size), n_samec=n_samec,
        # sanity (set-difference re-gauge + any baseline off-manifold tilt)
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
        # same-country spillover
        samec_raw_s=samec_raw * SEC3H if n_samec else np.nan,
        samec_resid_s=samec_res * SEC3H if n_samec else np.nan,
        rest_resid_s=rest_res * SEC3H if np.isfinite(rest_res) else np.nan,
        sd_boot_samec_s=sd_boot_samec * SEC3H if np.isfinite(sd_boot_samec) else np.nan,
        z_samec=z_samec,
        gauge_samec_plaus_s=gauge_samec_plaus * SEC3H if n_samec else np.nan,
        gauge_samec_real_s=gauge_samec_real * SEC3H if n_samec else np.nan,
        verdict=verdict,
    )

    per_race = pd.DataFrame(dict(
        excl_series=label, race_id=common,
        series_key=np.asarray(base["fd"].race_series, dtype=object)[iB],
        country=country, t=t, is_samec=is_samec,
        dv=dv, dv_affine=dec["dv_affine"], dv_quad=dec["dv_quad"],
        dv_resid=dec["dv_resid"], q_shape=des["q"],
        gauge_plaus=gauge_plaus, sd_boot_raw=sd_boot_raw,
        sd_boot_resid=sd_boot_res, ratio_resid_boot=ratio,
    ))
    return rec, per_race


# ---------------------------------------------------------------------------
# one study
# ---------------------------------------------------------------------------

def run_study(study: Path, series_filter: list[str] | None) -> pd.DataFrame | None:
    meta = json.loads((study / "loso_meta.json").read_text())
    sdirs = sorted(d for d in study.glob("series_*") if (d / "manifest.json").is_file())
    if series_filter:
        keep = {s.removeprefix("series_").removesuffix("_marathon")
                for s in series_filter}
        sdirs = [d for d in sdirs if d.name[len("series_"):] in keep]
    if not (study / "baseline" / "manifest.json").is_file() or not sdirs:
        print(f"  [skip] {study}: missing baseline/ or series_* fits")
        return None

    print(f"\n=== q01 LOSO summary: {study.parent.name}/{study.name} "
          f"({len(sdirs)} series) ===")
    base = load_bundle(study / "baseline")
    price = drift_price(base)
    delta_plaus = (price["omega_d"] / (2.0 * price["sd_debut"])
                   if price["sd_debut"] > 0 else np.inf)
    print(f"    baseline J={base['race_ids'].size}  omega_d={price['omega_d']:.4e} "
          f"SD(debut)={price['sd_debut']:.2f} yr (n_elig={price['n_elig']:,}) "
          f"-> delta_plaus={delta_plaus:.3e} /yr^2")

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
    for sdir in sdirs:
        label = sdir.name[len("series_"):]
        loso = load_bundle(sdir)
        rec, pr = analyse_series(label, base, loso, price, boot_full)
        recs.append(rec)
        frames.append(pr)

    df = pd.DataFrame(recs).sort_values("n_dropped_runners",
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
    cols = ["series", "ctry", "n_drop", "n_sc", "med|dv| s", "resid med s",
            "resid share", "gauge plaus med s", "R_kill", "boot med s",
            "ratio med", "samec raw s", "samec resid s", "z", "gauge sc s",
            "verdict"]
    lines = [f"# q01 LOSO gauge-aware summary -- {study.parent.name}/{study.name}",
             "",
             f"delta_plaus = {delta_plaus:.3e} /yr^2 (R=1 Gq re-pin; omega_d = "
             f"{price['omega_d']:.3e}, SD(debut) = {price['sd_debut']:.2f} yr). "
             f"Baseline slope(v~t) = {base_slope:+.6f}/yr"
             + (" -- OFF the beta=0 manifold; its tilt is removed as part of the "
                "affine component." if abs(base_slope) > 1e-6 else ".")
             + " All seconds at a 3:00:00 marathon. dv = v(-series) - "
             "v(baseline) on SURVIVING common races; resid = dv after "
             "projecting out {1, t, t^2}; samec = surviving races in the "
             "venue's country; ratio = per-race |resid| / residualized "
             "bootstrap SD. Verdicts: SIGNAL clears 2x noise and the priced "
             "worst-case gauge band; SIGNAL* clears noise but sits inside that "
             "(conservative) band -- check R_kill: << 1 means the realized "
             "re-pin between these two fits was negligible.",
             "",
             "| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        def _f(x, fmt="{:.0f}"):
            return fmt.format(x) if pd.notna(x) else "-"
        lines.append("| " + " | ".join([
            r["excl_series"], str(r["venue_country"]),
            f"{r['n_dropped_runners']:,}", f"{r['n_samec']}",
            _f(r["med_abs_dv_s"]), _f(r["med_abs_resid_s"]),
            _f(r["share_resid"], "{:.2f}"),
            _f(r["gauge_plaus_med_s"], "{:.1f}"), _f(r["R_kill"], "{:.2f}"),
            _f(r["boot_resid_med_s"], "{:.0f}"), _f(r["ratio_med"], "{:.2f}"),
            _f(r["samec_raw_s"], "{:+.0f}"), _f(r["samec_resid_s"], "{:+.0f}"),
            _f(r["z_samec"], "{:+.1f}"), _f(r["gauge_samec_plaus_s"], "{:.1f}"),
            r["verdict"],
        ]) + " |")
    (out / "q01_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"    wrote {out / 'q01_summary.parquet'} (+csv/md, per_race)")

    # ---- stdout ----------------------------------------------------------
    hdr = (f"{'series':>13} {'ctry':>4} {'n_drop':>7} {'n_sc':>4} | {'med|dv|':>7} "
           f"{'resid':>6} {'r_sh':>5} | {'g_plaus':>7} {'R_kill':>6} | {'boot':>5} "
           f"{'r_med':>5} | {'sc_raw':>7} {'sc_res':>7} {'z':>5} {'g_sc':>6}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for _, r in df.iterrows():
        def _n(x, fmt):
            return fmt.format(x) if pd.notna(x) else " " * (len(fmt.format(0)) - 1) + "-"
        print(f"{r['excl_series']:>13} {str(r['venue_country']):>4} "
              f"{r['n_dropped_runners']:>7,} {r['n_samec']:>4} | "
              f"{r['med_abs_dv_s']:>7.0f} {r['med_abs_resid_s']:>6.0f} "
              f"{r['share_resid']:>5.2f} | "
              f"{r['gauge_plaus_med_s']:>7.1f} {r['R_kill']:>6.2f} | "
              f"{_n(r['boot_resid_med_s'], '{:>5.0f}')} {_n(r['ratio_med'], '{:>5.2f}')} | "
              f"{_n(r['samec_raw_s'], '{:>7.0f}')} {_n(r['samec_resid_s'], '{:>7.0f}')} "
              f"{_n(r['z_samec'], '{:>5.1f}')} {_n(r['gauge_samec_plaus_s'], '{:>6.1f}')}  "
              f"{r['verdict']}")
    print()
    return df


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def discover_studies(root: Path) -> list[Path]:
    return sorted(p.parent for p in root.glob("*/*/loso_meta.json"))


ARGS: argparse.Namespace


def main() -> None:
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", nargs="+", default=None,
                    help="study dir(s) (default: every dir under "
                         "results/sensitivity/loso with a loso_meta.json).")
    ap.add_argument("--series", nargs="+", default=None,
                    help="restrict to these series (short label or full key).")
    ap.add_argument("--a-max", type=float, default=10.0,
                    help="career-age horizon for the aging-curve cross-check.")
    ARGS = ap.parse_args()

    studies = ([Path(s) for s in ARGS.study] if ARGS.study
               else discover_studies(OUT_ROOT))
    if not studies:
        raise SystemExit(f"no LOSO studies found under {OUT_ROOT}")
    for s in studies:
        if not (s / "loso_meta.json").is_file():
            raise SystemExit(f"not a LOSO study dir (no loso_meta.json): {s}")
        run_study(s, ARGS.series)


if __name__ == "__main__":
    main()
