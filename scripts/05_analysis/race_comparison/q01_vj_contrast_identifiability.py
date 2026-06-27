"""Gauge bound on the cross-sex race-factor contrast  Dv_j = v_M,j - v_W,j.

Everything below is computed PER SHARED RACE j first; the summary table then
reports medians / percentiles OVER the shared races. Nothing is a slice-level
aggregate compared against another slice-level aggregate: `ratio_med` is the
median of the per-race ratios, not `shift_med / boot_med`.

SETUP AND NOTATION
------------------
Two independent production fits (same cohort, mrc, window; sex M and W).
  j = 1..n     races present in BOTH fits ("shared races")
  Dv_j         = v_M,j - v_W,j      the cross-sex contrast (log-time units)
  t_j          race date in fractional years;  t~_j = t_j - mean_j(t)
  f_i          athlete i's debut date (first race, years); for sex s,
               f~ = f centered on s's drift-eligible set, SD_s = SD(f~)
  omega_d_s    fitted drift-spread sqrt(params["omega_d2"]) of sex s's fit

WHY Dv CARRIES A CONVENTION PART (the gauge family)
---------------------------------------------------
Each fit is pinned to one representative of a flat-direction family
(docs/model_derivation.md S7): mean(v)=0 + beta=0 in-fit, and the min-norm
d_i prior (d orthogonal to {1, debut, entry_age} on that fit's eligible set).
Because the two sexes are pinned independently on different populations, the
contrast can be moved, with NO change to either fit's likelihood, by

    Dv_j  ->  Dv_j + a + b*t~_j + delta*t~_j^2                          (G)

  a (G0 level)    each mean(v)=0 is over its OWN race set; the contrast
                  level is absorbed in the sex-specific u levels.
                  Pure convention -- UNPRICED, any value achievable.
  b (G1 linear)   each beta=0 zeroes slope(v ~ t) over its OWN races; a real
                  M-W difference in secular trend IS this direction.
                  Pure convention -- UNPRICED.
  delta (Gq quad) v_j += delta*t_j^2 is a flat direction because, with
                  career age A_n = t_j - f_i,
                      delta*t^2 = delta*A_n^2 + 2*delta*A_n*f_i + delta*f_i^2,
                  so it is repaid by the aging-curve quadratic (A_n^2 term),
                  a per-athlete drift shift d_i -= 2*delta*f_i (the A_n-linear
                  term), and u_i (the constants). The d-prior pins it via
                  cov(d, debut)=0, so unlike a/b it has a PRICE: re-pinning
                  sex s by delta leaves a debut-correlated drift component of
                  SD = 2*|delta|*SD_s in d, measured against omega_d_s
                  (same pricing logic as aging_trend/q02 for the fan slope).

Only the residual of Dv after projecting out {1, t~, t~^2} is gauge-free.

STEP 1 -- DECOMPOSITION (per race)
----------------------------------
OLS over the shared races:

    Dv_j = beta0 + beta1*t~_j + beta2*t~_j^2 + r_j
    q_j  = t~_j^2 - [OLS projection of t~^2 on {1, t~}]     (q _|_ {1, t~})

    Dv_j = affine_j + beta2*q_j + r_j        (exact identity; affine_j is the
                                              part in span{1, t~})
    share_resid = Var_j(r) / Var_j(Dv)       gauge-free fraction = the ceiling
                                              on any reportable cross-sex signal

STEP 2 -- PRICE THE QUADRATIC FREEDOM (the only priced direction)
-----------------------------------------------------------------
    R_s(delta)  = 2*|delta| * SD_s / omega_d_s      drift cost of a re-pin,
                                                    as a fraction of the
                                                    fitted drift spread
    delta_kill  = -beta2                            re-pin that zeroes the
                                                    FITTED quad part
    R_min       = min(R_M(delta_kill), R_W(delta_kill))
                  q02 convention: < 0.5 FRAGILE (an undetectable drift-debut
                  correlation reproduces it), > 1.5 robust
    delta_plaus = max_s omega_d_s / (2*SD_s)        the R=1 re-pin via the
                                                    cheapest sex: a quad gauge
                                                    motion whose drift cost
                                                    equals the WHOLE fitted
                                                    drift spread
    shift_j     = delta_plaus * |q_j|               per-race plausible gauge
                                                    motion of Dv_j (largest at
                                                    the date extremes)

STEP 3 -- COMPARE WITH BOOTSTRAP (per race, then median over races)
-------------------------------------------------------------------
The two fits' athlete-weight bootstraps are independent, so replicate pairs
are formed by random pairing (n_boot draws of (b_M, b_W)):

    sd_boot_j = SD_b( v_M,j^(b_M) - v_W,j^(b_W) )   per-race sampling SD of Dv_j
    ratio_j   = shift_j / sd_boot_j                 gauge freedom vs sampling
                                                    noise AT RACE j

Table columns: shift_med = median_j(shift_j); boot_med = median_j(sd_boot_j);
ratio_med / ratio_p90 / ratio_max = median / 90th pct / max of ratio_j over j.
Also gauge_vs_boot_c2 = delta_plaus / SD_b(beta2^(b)) compares the freedom in
the quadratic COEFFICIENT to its sampling SD (a coefficient-level analogue).

READING THE VERDICT
-------------------
ratio_med >> 1 (GAUGE > BOOT): even the PRICED slice of the gauge freedom
moves per-race Dv_j by more than its bootstrap SD -- the "M and W v_j are not
directly comparable" warning is load-bearing, and bootstrap CIs on Dv_j
understate the true cross-sex uncertainty. ratio_med << 1: the priced part is
subdominant to sampling noise -- but a and b remain UNBOUNDED regardless, so
cross-sex v levels and date trends are unreportable either way; only the
residual pattern r_j could ever be quoted.

CAVEATS. (i) q is orthogonalized on the SHARED race set while each fit's own
gauges are imposed on its FULL race set -- adequate for a bound, not an exact
reprojection. (ii) Gq is exactly flat only for d-eligible athletes; for
ineligible (short-span) athletes it is approximate. (iii) Eligible sets (and
hence the price) differ sharply between mrc2 and mrc5 -- compare rows, do not
pool. (iv) Units are log-time; x180 ~= minutes at a 3:00 marathon.

NO REFIT -- reads saved point fits + their bootstraps.

Outputs (results/analysis/race_comparison/):
    vj_contrast_identifiability.parquet/.csv/.md   one summary row per (cohort, mrc)
    vj_contrast_per_race.parquet                   per shared race: decomposition,
                                                   plausible shift, boot SD, ratios
    + a readable table to stdout.

Run::

    python scripts/05_analysis/race_comparison/q01_vj_contrast_identifiability.py
    python scripts/05_analysis/race_comparison/q01_vj_contrast_identifiability.py --n-boot 4000
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import load_slice, registry  # noqa: E402
from marathon_decomp.config import DATA_DIR, RESULTS_DIR  # noqa: E402

MODELS_ROOT = RESULTS_DIR / "models"
OUT_ROOT = RESULTS_DIR / "analysis" / "race_comparison"

COHORTS = ["ALL", "Po10"]
MRCS = ["mrc2", "mrc5"]
REF_MARATHON_MIN = 180.0   # minutes at a 3:00 marathon, for readable units
NO_GAP_EPS = 1e-12


# --------------------------------------------------------------------------- #
# fit loading (find_fit_dir mirrors aging_trend/p01_aging_curve.py)           #
# --------------------------------------------------------------------------- #
def find_fit_dir(slug: str, model: str, nutag: str) -> Path | None:
    """First `{model}_{nutag}_best__{hash}` dir under results/models/{slug} with a fit.pkl."""
    sd = MODELS_ROOT / slug
    if not sd.is_dir():
        return None
    cands = sorted(sd.glob(f"{model}_{nutag}_best__*"))
    cands = [c for c in cands if (c / "fit.pkl").is_file()]
    return cands[0] if cands else None


def _frac_years(dates) -> np.ndarray:
    """datetime64 array -> fractional Julian-ish years (matches p01 convention)."""
    rd = pd.DatetimeIndex(pd.to_datetime(dates))
    return (rd.year + (rd.dayofyear - 1) / 365.25).to_numpy(np.float64)


def load_vfit(slug: str, model: str, nutag: str, data_version: str) -> dict | None:
    """Point fit -> {v, race_ids, t (years), fd, fit_dir, payload}."""
    fit_dir = find_fit_dir(slug, model, nutag)
    if fit_dir is None:
        print(f"  [skip] no point fit for {slug} ({model}_{nutag})")
        return None
    payload = pickle.load(open(fit_dir / "fit.pkl", "rb"))
    fd = load_slice(payload["spec"], payload.get("data_version", data_version))
    return dict(
        fit_dir=fit_dir, payload=payload, fd=fd,
        v=np.asarray(payload["params"]["v"], np.float64),
        race_ids=np.asarray(fd.race_ids, np.int64),
        t=_frac_years(fd.race_date),
    )


def load_boot_v(fit_dir: Path, shared_ids: np.ndarray) -> np.ndarray | None:
    """Bootstrap replicates of v on the shared races -> (R, n_shared) or None."""
    p = fit_dir / "bootstrap" / "race_factors.parquet"
    if not p.is_file():
        return None
    rf = pd.read_parquet(p, columns=["run_id", "race_id", "v"])
    wide = (rf[rf["run_id"] > 0]
            .pivot(index="run_id", columns="race_id", values="v")
            .reindex(columns=shared_ids))
    if wide.empty:
        return None
    arr = wide.to_numpy(np.float64)
    if np.isnan(arr).any():
        n_bad = int(np.isnan(arr).any(axis=0).sum())
        print(f"    [warn] {n_bad} shared races missing from bootstrap "
              f"{p.parent.parent.name}; dropping replicate NaNs by race is "
              "not possible -- those races get NaN boot SD")
    return arr


def drift_inputs(vfit: dict) -> dict:
    """omega_d and SD of centered debut date over the drift-eligible set.

    The Gq re-pin moves d_i by -2*delta*f_i; the min-norm prior keeps d
    orthogonal to debut over the eligible set, so the injected component that
    must live in d has SD = 2*|delta|*SD(debut - mean) there. Debut is data
    (first race date), so no DOB-known mask is needed (unlike q02's entry age).
    """
    fd = vfit["fd"]
    model = registry.load_fit(vfit["fit_dir"], fd)
    omega_d2 = float(model.params["omega_d2"])

    t_cell = vfit["t"][np.asarray(fd.col_idx, np.int64)]
    f_cell = t_cell - np.asarray(fd.A_n, np.float64)   # debut date, const per athlete
    debut = np.full(len(fd.athlete_ids), np.nan)
    debut[np.asarray(fd.row_idx, np.int64)] = f_cell    # any cell of the athlete

    elig = np.asarray(model.eligible_d, bool)
    f_e = debut[elig & np.isfinite(debut)]
    sd_f = float(np.std(f_e - f_e.mean(), ddof=0)) if f_e.size else 0.0
    return dict(omega_d=float(np.sqrt(max(omega_d2, 0.0))),
                sd_debut=sd_f, n_elig=int(elig.sum()))


# --------------------------------------------------------------------------- #
# decomposition of the contrast on the shared races                           #
# --------------------------------------------------------------------------- #
def poly_decompose(t: np.ndarray, dv: np.ndarray) -> dict:
    """Split dv into affine + orthogonalized-quadratic + residual over races.

    Returns the raw quadratic coefficient c2 (per yr^2, the Gq-movable
    quantity), the orthogonalized quadratic shape q (t~^2 residualized
    against {1, t~}), and the three components (affine includes the affine
    shadow of the t^2 term, so dv = affine + c2*q + resid exactly).
    """
    t_c = t - t.mean()
    X_aff = np.column_stack([np.ones_like(t_c), t_c])
    X3 = np.column_stack([X_aff, t_c ** 2])

    beta3, *_ = np.linalg.lstsq(X3, dv, rcond=None)
    c2 = float(beta3[2])
    smooth = X3 @ beta3
    resid = dv - smooth

    coef_q, *_ = np.linalg.lstsq(X_aff, t_c ** 2, rcond=None)
    q = t_c ** 2 - X_aff @ coef_q              # quad shape orthogonal to affine
    dv_quad = c2 * q
    dv_affine = smooth - dv_quad               # in span{1, t~}

    return dict(t_c=t_c, c2=c2, q=q,
                dv_affine=dv_affine, dv_quad=dv_quad, dv_resid=resid,
                slope=float(beta3[1]), level=float(dv.mean()))


def _verdict(ratio_med: float, has_boot: bool) -> str:
    if not has_boot:
        return "POINT-ONLY"
    if not np.isfinite(ratio_med):
        return "UNDETERMINED"
    if ratio_med > 1.0:
        return "GAUGE > BOOT"
    if ratio_med > 0.3:
        return "COMPARABLE"
    return "BOOT > GAUGE"


# --------------------------------------------------------------------------- #
# per-(cohort, mrc) analysis                                                  #
# --------------------------------------------------------------------------- #
def analyse(cohort: str, mrc: str, *, model: str, nutag: str, data_version: str,
            n_boot: int, rng) -> tuple[dict, pd.DataFrame] | None:
    vM = load_vfit(f"{cohort}_M_14-25_{mrc}", model, nutag, data_version)
    vW = load_vfit(f"{cohort}_W_14-25_{mrc}", model, nutag, data_version)
    if vM is None or vW is None:
        return None

    shared_ids, iM, iW = np.intersect1d(vM["race_ids"], vW["race_ids"],
                                        return_indices=True)
    if shared_ids.size < 10:
        print(f"  [skip] {cohort} {mrc}: only {shared_ids.size} shared races")
        return None
    dv = vM["v"][iM] - vW["v"][iW]
    t = vM["t"][iM]                            # same race_id -> same date

    dec = poly_decompose(t, dv)
    var_dv = float(np.var(dv))
    shares = {k: float(np.var(dec[f"dv_{k}"]) / var_dv) if var_dv > NO_GAP_EPS
              else np.nan for k in ("affine", "quad", "resid")}

    # ---- price the quadratic direction (the only priced one) --------------
    di_M, di_W = drift_inputs(vM), drift_inputs(vW)
    delta_kill = -dec["c2"]

    def _R(di, delta):
        return 2.0 * abs(delta) * di["sd_debut"] / di["omega_d"] \
            if di["omega_d"] > 0 else np.inf

    def _plaus(di):
        return di["omega_d"] / (2.0 * di["sd_debut"]) if di["sd_debut"] > 0 else np.inf

    R_M, R_W = _R(di_M, delta_kill), _R(di_W, delta_kill)
    delta_plaus = max(_plaus(di_M), _plaus(di_W))   # cheapest-sex R=1 re-pin
    shift_plaus = delta_plaus * np.abs(dec["q"])    # per-race quad gauge motion

    # ---- bootstrap comparison (random pairing of independent replicates) --
    bM = load_boot_v(vM["fit_dir"], shared_ids)
    bW = load_boot_v(vW["fit_dir"], shared_ids)
    has_boot = bM is not None and bW is not None
    sd_boot = np.full(shared_ids.size, np.nan)
    c2_boot_sd = np.nan
    if has_boot:
        idx_M = rng.integers(0, bM.shape[0], n_boot)
        idx_W = rng.integers(0, bW.shape[0], n_boot)
        dv_b = bM[idx_M] - bW[idx_W]               # (n_boot, n_shared)
        sd_boot = dv_b.std(axis=0, ddof=1)
        # quadratic coefficient sampling SD via one shared pseudo-inverse
        X3 = np.column_stack([np.ones_like(dec["t_c"]), dec["t_c"], dec["t_c"] ** 2])
        beta_b = dv_b @ np.linalg.pinv(X3).T       # (n_boot, 3)
        c2_boot_sd = float(beta_b[:, 2].std(ddof=1))

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = shift_plaus / sd_boot
        ratio_fit = np.abs(dec["dv_affine"] + dec["dv_quad"]) / sd_boot
    fin = np.isfinite(ratio)

    rec = dict(
        cohort=cohort, mrc=mrc, n_shared=int(shared_ids.size),
        t_span=float(t.max() - t.min()), has_boot=has_boot,
        # realized smooth components (convention residue as currently pinned)
        level=dec["level"], slope=dec["slope"], c2=dec["c2"],
        sd_dv=float(np.std(dv)),
        sd_affine=float(np.std(dec["dv_affine"])),
        sd_quad=float(np.std(dec["dv_quad"])),
        sd_resid=float(np.std(dec["dv_resid"])),
        share_affine=shares["affine"], share_quad=shares["quad"],
        share_resid=shares["resid"],
        # quadratic pricing
        delta_kill=delta_kill, R_M=R_M, R_W=R_W, R_min=min(R_M, R_W),
        omega_d_M=di_M["omega_d"], omega_d_W=di_W["omega_d"],
        sd_debut_M=di_M["sd_debut"], sd_debut_W=di_W["sd_debut"],
        n_elig_M=di_M["n_elig"], n_elig_W=di_W["n_elig"],
        delta_plaus=delta_plaus,
        shift_plaus_med=float(np.median(shift_plaus)),
        shift_plaus_max=float(shift_plaus.max()),
        shift_plaus_med_min=float(np.median(shift_plaus)) * REF_MARATHON_MIN,
        # bootstrap comparison
        sd_boot_med=float(np.nanmedian(sd_boot)) if has_boot else np.nan,
        c2_boot_sd=c2_boot_sd,
        gauge_vs_boot_c2=(delta_plaus / c2_boot_sd
                          if has_boot and c2_boot_sd > 0 else np.nan),
        ratio_med=float(np.median(ratio[fin])) if fin.any() else np.nan,
        ratio_p90=float(np.percentile(ratio[fin], 90)) if fin.any() else np.nan,
        ratio_max=float(ratio[fin].max()) if fin.any() else np.nan,
    )
    rec["verdict"] = _verdict(rec["ratio_med"], has_boot)

    per_race = pd.DataFrame(dict(
        cohort=cohort, mrc=mrc, race_id=shared_ids, t=t, dv=dv,
        dv_affine=dec["dv_affine"], dv_quad=dec["dv_quad"],
        dv_resid=dec["dv_resid"], q_shape=dec["q"],
        shift_plaus=shift_plaus, sd_boot=sd_boot,
        ratio_plaus_boot=ratio, ratio_fitted_smooth_boot=ratio_fit,
    ))
    return rec, per_race


# --------------------------------------------------------------------------- #
# main                                                                        #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="full", help="registry model tag (AxD = 'full').")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--cohorts", nargs="+", default=COHORTS)
    ap.add_argument("--mrcs", nargs="+", default=MRCS)
    ap.add_argument("--n-boot", type=int, default=4000,
                    help="random-pairing bootstrap draws for the contrast SDs.")
    ap.add_argument("--seed", type=int, default=20260611)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    recs, per_race_frames = [], []
    for cohort in args.cohorts:
        for mrc in args.mrcs:
            print(f"analysing {cohort} {mrc} ...")
            out = analyse(cohort, mrc, model=args.model, nutag=args.nutag,
                          data_version=args.data_version,
                          n_boot=args.n_boot, rng=rng)
            if out is not None:
                recs.append(out[0])
                per_race_frames.append(out[1])

    if not recs:
        raise SystemExit("no (cohort, mrc) pairs produced a result")

    df = pd.DataFrame(recs).sort_values(["cohort", "mrc"]).reset_index(drop=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_ROOT / "vj_contrast_identifiability.parquet", index=False)
    print(f"wrote {OUT_ROOT / 'vj_contrast_identifiability.parquet'}  ({len(df)} rows)")

    num = df.select_dtypes(include="number").columns
    df_csv = df.copy()
    df_csv[num] = df_csv[num].round(8)
    df_csv.to_csv(OUT_ROOT / "vj_contrast_identifiability.csv", index=False)
    print(f"wrote {OUT_ROOT / 'vj_contrast_identifiability.csv'}")

    pr = pd.concat(per_race_frames, ignore_index=True)
    # attach race metadata for readability
    comps_path = DATA_DIR / "competitions.parquet"
    if comps_path.is_file():
        comps = pd.read_parquet(comps_path)[["race_id", "series_key", "country", "year"]]
        pr = pr.merge(comps, on="race_id", how="left")
    pr.to_parquet(OUT_ROOT / "vj_contrast_per_race.parquet", index=False)
    print(f"wrote {OUT_ROOT / 'vj_contrast_per_race.parquet'}  ({len(pr)} rows)")

    # ---- readable Markdown -------------------------------------------------
    cols = ["Cohort", "Mrc", "n_shared", "sd(Dv)", "resid share",
            "delta_kill", "R_min", "delta_plaus", "shift_plaus med (min@3h)",
            "boot sd med", "ratio med", "ratio p90", "Verdict"]
    lines = ["# Cross-sex v_j contrast gauge bound (Dv = v_M - v_W on shared races)",
             "",
             "Level and linear-in-date parts of Dv are PURE CONVENTION (G0/G1,",
             "unpriced -- any value achievable); the table prices only the",
             "quadratic (Gq) freedom against omega_d. `ratio` compares the R=1",
             "plausible per-race quadratic gauge shift to the per-race bootstrap",
             "SD of the contrast.",
             "",
             "| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join([
            r["cohort"], r["mrc"][3:], f"{r['n_shared']}",
            f"{r['sd_dv']:.4f}",
            f"{r['share_resid']:.2f}" if pd.notna(r["share_resid"]) else "-",
            f"{r['delta_kill']:+.2e}", f"{r['R_min']:.2f}",
            f"{r['delta_plaus']:.2e}",
            f"{r['shift_plaus_med']:.4f} ({r['shift_plaus_med_min']:.2f})",
            f"{r['sd_boot_med']:.4f}" if pd.notna(r["sd_boot_med"]) else "-",
            f"{r['ratio_med']:.2f}" if pd.notna(r["ratio_med"]) else "-",
            f"{r['ratio_p90']:.2f}" if pd.notna(r["ratio_p90"]) else "-",
            r["verdict"],
        ]) + " |")
    note = (
        "\n_`sd(Dv)` = SD of the raw contrast over shared races (log units; x180 "
        "= minutes at 3:00). `resid share` = gauge-free fraction of the contrast "
        "variance (residual after projecting out {1, t, t^2}) -- the ceiling on "
        "any reportable cross-sex v_j pattern. `delta_kill` = quadratic-in-date "
        "coefficient that zeroes the fitted quadratic part; `R_min` = its drift "
        "cost as a fraction of the fitted drift spread, min over sexes (q02 "
        "convention: <0.5 fragile, >1.5 robust). `delta_plaus` = the R=1 "
        "quadratic re-pin via the cheapest sex; `shift_plaus med` = median "
        "per-race |delta_plaus * q(t)| it induces. `ratio` = that shift / "
        "per-race bootstrap SD of Dv: >1 means the priced gauge freedom alone "
        "exceeds sampling noise. The affine freedom is infinite regardless -- "
        "cross-sex v levels and date trends are never reportable._\n")
    (OUT_ROOT / "vj_contrast_identifiability.md").write_text(
        "\n".join(lines) + "\n" + note, encoding="utf-8")
    print(f"wrote {OUT_ROOT / 'vj_contrast_identifiability.md'}\n")

    # ---- stdout ------------------------------------------------------------
    hdr = (f"{'cohort':>6} {'mrc':>5} | {'n_sh':>5} {'sd(Dv)':>7} {'res_sh':>6} | "
           f"{'d_kill':>9} {'R_min':>6} {'d_plaus':>9} {'shift_med':>9} | "
           f"{'boot_med':>8} {'r_med':>6} {'r_p90':>6}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for _, r in df.iterrows():
        print(f"{r['cohort']:>6} {r['mrc']:>5} | {r['n_shared']:>5} "
              f"{r['sd_dv']:>7.4f} "
              f"{(r['share_resid'] if pd.notna(r['share_resid']) else float('nan')):>6.2f} | "
              f"{r['delta_kill']:>+9.2e} {r['R_min']:>6.2f} "
              f"{r['delta_plaus']:>9.2e} {r['shift_plaus_med']:>9.4f} | "
              f"{(r['sd_boot_med'] if pd.notna(r['sd_boot_med']) else float('nan')):>8.4f} "
              f"{(r['ratio_med'] if pd.notna(r['ratio_med']) else float('nan')):>6.2f} "
              f"{(r['ratio_p90'] if pd.notna(r['ratio_p90']) else float('nan')):>6.2f}  "
              f"{r['verdict']}")
    print()


if __name__ == "__main__":
    main()
