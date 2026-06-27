"""Gauge bound on the cross-sex entry-age FAN difference  Dg = g_M - g_W.

THE QUESTION. The figure shows entry-35 (blue) curves nearly coinciding across
sex while entry-65 (red) curves split (men appear to decline faster). Because
*all* entry-age dependence lives in the fan g(A_n) = gamma . B(A_n), that
observation is fundamentally a statement about the cross-sex fan difference

    Dg(A_n) = g_M(A_n) - g_W(A_n).

CAN WE CLAIM IT? The fan has a gauge freedom (the `Ggamma` direction, see
`docs/model_derivation.md` S7.4):

    g_s(A_n) -> g_s(A_n) + c_s * A_n          (a LINEAR-in-A_n shift)
    paid for by   d_i -> d_i - c_s * (A_e_i - mean_Ae_s)   on sex s's drifts.

Each fit's `d_i` prior pins c_s by the min-norm choice cov(d, A_e)=0 over the
eligible set -- a modelling assumption (drift uncorrelated with entry age), not
data. So the cross-sex linear offset delta = c_M - c_W is convention, and it
shifts Dg by delta * A_n. Consequences:

  * the CURVATURE (bend) of Dg is gauge-FREE -> identified;
  * the SLOPE (linear part) of Dg is gauge-pinned -> convention.

This script does two things per (cohort, mrc):

1. SPLIT  Dg(A_n) = delta_fit * A_n + Dg_curv(A_n)   (LS slope through origin;
   g(0)=0 so no intercept). Reports the value at A_max carried by each part and
   `curv_share` = |curv_end| / |Dg_end|. If the conclusion lives in `Dg_curv`,
   it is claimable with no bracket. Bootstrap CI (mrc5) says whether the
   curvature part even excludes 0.

2. BOUND the slope with omega_d. The delta that would ERASE the gap at A_max is
   delta_kill = -Dg(A_max) / A_max. Realising it injects a drift<->entry-age
   gradient of slope delta_kill into one sex, i.e. a drift component of SD
   |delta_kill| * SD(A_e_tilde) on the eligible set. Compare to the fitted drift
   spread omega_d:

       R_s = |delta_kill| * SD(A_e_tilde_s) / omega_d_s.

   R >> 1: erasing the gap needs a drift/entry-age correlation far larger than
           the entire fitted drift spread -> implausible -> claim ROBUST.
   R << 1: a correlation smaller than the drift noise floor (undetectable) would
           suffice -> claim FRAGILE (rests on the cov(d,A_e)=0 assumption).

   R is the same whether read per-unit-entry-age or at entry 65 (the entry-age
   factor cancels between the gap and the gauge knob), so it is ONE number per
   (cohort, mrc). The cheapest reattribution picks the sex with the smaller R,
   so the verdict keys on R_min = min(R_M, R_W).

NO REFIT -- everything is read from the saved point fits + their bootstraps.

Outputs (results/analysis/aging_trend/):
    fan_identifiability.parquet   all columns
    fan_identifiability.csv       rounded, spreadsheet-friendly
    fan_identifiability.md        readable table + verdicts
    + a readable table to stdout.

Run::

    python scripts/05_analysis/aging_trend/q02_fan_identifiability.py
    python scripts/05_analysis/aging_trend/q02_fan_identifiability.py --a-max 10 --n-boot 4000
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir (p01)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from p01_aging_curve import OUT_ROOT, find_fit_dir, reconstruct  # noqa: E402

from marathon_decomp import load_slice, registry  # noqa: E402

COHORTS = ["ALL", "Po10"]
MRCS = ["mrc2", "mrc5"]
ENTRY_REF = 65.0   # illustrative entry age for the "fan contribution" magnitude
NO_GAP_EPS = 1e-9  # |Dg_end| below this -> nothing to claim


# --------------------------------------------------------------------------- #
# drift-side inputs for the bound: omega_d and SD(A_e_tilde) over eligible set #
# --------------------------------------------------------------------------- #
def drift_inputs(slug: str, model_tag: str, nutag: str, data_version: str) -> dict | None:
    """Load the saved fit as a live model and return omega_d, SD(centered entry
    age over the drift-eligible & DOB-known athletes), and the eligible count."""
    fit_dir = find_fit_dir(slug, model_tag, nutag)
    if fit_dir is None:
        return None
    payload = pickle.load(open(fit_dir / "fit.pkl", "rb"))
    fd = load_slice(payload["spec"], payload.get("data_version", data_version))
    model = registry.load_fit(fit_dir, fd)

    omega_d2 = float(model.params["omega_d2"])
    A_e = np.asarray(fd.A_e, np.float64)
    mean_Ae = float(np.nanmean(A_e))
    a_tilde = A_e - mean_Ae                       # matches the model's centering
    elig = np.asarray(model.eligible_d, bool)
    known = np.isfinite(A_e)                       # unknown-DOB athletes have a_tilde=0
    mask = elig & known                            # the only athletes Ggamma can move
    sd_atilde = float(np.std(a_tilde[mask], ddof=0)) if mask.any() else 0.0
    return dict(omega_d=float(np.sqrt(max(omega_d2, 0.0))),
                sd_atilde=sd_atilde, mean_Ae=mean_Ae,
                n_elig=int(elig.sum()), n_elig_known=int(mask.sum()))


def _verdict(r_min: float, no_gap: bool, curv_excl0: bool) -> str:
    if no_gap:
        return "NO GAP"
    if curv_excl0:
        return "ROBUST (curvature)"
    if not np.isfinite(r_min):
        return "UNDETERMINED"
    if r_min < 0.5:
        return "FRAGILE"
    if r_min < 1.5:
        return "BORDERLINE"
    return "ROBUST"


def analyse(cohort: str, mrc: str, *, model_tag: str, nutag: str,
            data_version: str, a_max: float, n_boot: int, rng) -> dict | None:
    fan_M = reconstruct(f"{cohort}_M_14-25_{mrc}", model=model_tag, nutag=nutag,
                        data_version=data_version, a_max=a_max, allow_point_only=True)
    fan_W = reconstruct(f"{cohort}_W_14-25_{mrc}", model=model_tag, nutag=nutag,
                        data_version=data_version, a_max=a_max, allow_point_only=True)
    if fan_M is None or fan_W is None:
        print(f"  [skip] {cohort} {mrc}: missing M or W fit")
        return None

    A = fan_M.A_grid                               # shared grid (same a_max, n_grid)
    A_max = float(A[-1])
    AA = float(A @ A)

    # ---- point split of Dg = g_M - g_W into linear (gauge) + curvature (free)
    Dg = fan_M.fan_point - fan_W.fan_point         # per unit centered entry age
    delta_fit = float((A @ Dg) / AA)               # LS slope through origin
    Dg_end = float(Dg[-1])
    lin_end = delta_fit * A_max
    curv = Dg - delta_fit * A
    curv_end = float(curv[-1])
    curv_share = abs(curv_end) / abs(Dg_end) if abs(Dg_end) > NO_GAP_EPS else np.nan

    # ---- omega_d bound: delta that zeroes the gap at A_max, drift cost vs omega_d
    delta_kill = -Dg_end / A_max
    di_M = drift_inputs(f"{cohort}_M_14-25_{mrc}", model_tag, nutag, data_version)
    di_W = drift_inputs(f"{cohort}_W_14-25_{mrc}", model_tag, nutag, data_version)
    if di_M is None or di_W is None:
        return None
    R_M = abs(delta_kill) * di_M["sd_atilde"] / di_M["omega_d"] if di_M["omega_d"] > 0 else np.inf
    R_W = abs(delta_kill) * di_W["sd_atilde"] / di_W["omega_d"] if di_W["omega_d"] > 0 else np.inf
    R_min = min(R_M, R_W)

    # ---- bootstrap CI on Dg_end / slope / curvature (independent fits -> random pair)
    rec: dict = {}
    has_boot = fan_M.has_boot and fan_W.has_boot
    curv_excl0 = no_gap_boot = False
    if has_boot:
        bM, bW = fan_M.fan_boot, fan_W.fan_boot     # (R_s, G)
        iM = rng.integers(0, bM.shape[0], n_boot)
        iW = rng.integers(0, bW.shape[0], n_boot)
        dg_b = bM[iM] - bW[iW]                       # (n_boot, G)
        slope_b = (dg_b @ A) / AA
        dg_end_b = dg_b[:, -1]
        curv_end_b = dg_end_b - slope_b * A_max
        lo_dg, hi_dg = np.percentile(dg_end_b, [2.5, 97.5])
        lo_sl, hi_sl = np.percentile(slope_b, [2.5, 97.5])
        lo_cv, hi_cv = np.percentile(curv_end_b, [2.5, 97.5])
        curv_excl0 = (lo_cv > 0) or (hi_cv < 0)
        no_gap_boot = (lo_dg <= 0 <= hi_dg)
        rec.update(dg_end_lo=lo_dg, dg_end_hi=hi_dg,
                   slope_lo=lo_sl, slope_hi=hi_sl,
                   curv_end_lo=lo_cv, curv_end_hi=hi_cv)

    no_gap = (abs(Dg_end) <= NO_GAP_EPS) or no_gap_boot
    verdict = _verdict(R_min, no_gap, curv_excl0)

    # illustrative magnitude: fan-only contribution to the entry-ENTRY_REF curve,
    # using the across-sex mean entry age for the (entry - mean) lever (log + %).
    mean_Ae = 0.5 * (di_M["mean_Ae"] + di_W["mean_Ae"])
    lever = ENTRY_REF - mean_Ae
    fan_end_log = lever * Dg_end
    rec.update(
        cohort=cohort, mrc=mrc, a_max=A_max, has_boot=has_boot,
        dg_end=Dg_end, slope_delta=delta_fit, lin_end=lin_end,
        curv_end=curv_end, curv_share=curv_share, curv_excl0=curv_excl0,
        delta_kill=delta_kill,
        omega_d_M=di_M["omega_d"], omega_d_W=di_W["omega_d"],
        sd_atilde_M=di_M["sd_atilde"], sd_atilde_W=di_W["sd_atilde"],
        n_elig_M=di_M["n_elig"], n_elig_W=di_W["n_elig"],
        R_M=R_M, R_W=R_W, R_min=R_min,
        mean_Ae=mean_Ae, entry_ref=ENTRY_REF,
        fan_end_log_at_ref=fan_end_log,
        fan_end_pct_at_ref=(np.exp(fan_end_log) - 1.0) * 100.0,
        verdict=verdict,
    )
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="full", help="registry model tag (AxD = 'full').")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--a-max", type=float, default=10.0,
                    help="career-age (A_n) horizon at which the gap is read.")
    ap.add_argument("--n-boot", type=int, default=4000,
                    help="random-pairing bootstrap draws for the difference CIs.")
    ap.add_argument("--seed", type=int, default=20260610)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    recs = []
    for cohort in COHORTS:
        for mrc in MRCS:
            r = analyse(cohort, mrc, model_tag=args.model, nutag=args.nutag,
                        data_version=args.data_version, a_max=args.a_max,
                        n_boot=args.n_boot, rng=rng)
            if r is not None:
                recs.append(r)

    if not recs:
        raise SystemExit("no (cohort, mrc) pairs produced a result")

    df = pd.DataFrame(recs).sort_values(["cohort", "mrc"]).reset_index(drop=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_ROOT / "fan_identifiability.parquet", index=False)
    print(f"wrote {OUT_ROOT / 'fan_identifiability.parquet'}  ({len(df)} rows)")

    num = df.select_dtypes(include="number").columns
    df_csv = df.copy()
    df_csv[num] = df_csv[num].round(6)
    df_csv.to_csv(OUT_ROOT / "fan_identifiability.csv", index=False)
    print(f"wrote {OUT_ROOT / 'fan_identifiability.csv'}")

    # ---- readable Markdown ------------------------------------------------- #
    def ci(r, key):
        lo, hi = r.get(f"{key}_lo"), r.get(f"{key}_hi")
        return f" [{lo:+.4f}, {hi:+.4f}]" if pd.notna(lo) else ""

    cols = ["Cohort", "Mrc", "Dg(A_max)", "curv share", "curv excl 0",
            "delta_kill", "R_M", "R_W", "R_min", "fan@65 %", "Verdict"]
    lines = ["# Cross-sex fan-difference gauge bound (Dg = g_M - g_W)", "",
             "| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join([
            r["cohort"], r["mrc"][3:],
            f"{r['dg_end']:+.4f}" + ci(r, "dg_end"),
            f"{r['curv_share']:.2f}" if pd.notna(r["curv_share"]) else "-",
            ("yes" if r["curv_excl0"] else "no") if r["has_boot"] else "n/a",
            f"{r['delta_kill']:+.5f}",
            f"{r['R_M']:.2f}", f"{r['R_W']:.2f}", f"{r['R_min']:.2f}",
            f"{r['fan_end_pct_at_ref']:+.2f}",
            r["verdict"],
        ]) + " |")
    note = (
        "\n_`Dg(A_max)` = cross-sex fan difference per unit centered entry age at "
        "the A_max horizon (+ => men's fan higher), 95% random-pairing bootstrap "
        "CI (mrc5). `curv share` = |curvature part| / |Dg(A_max)| (the gauge-free "
        "fraction); `curv excl 0` = does the curvature-part CI exclude 0. "
        "`delta_kill` = the Ggamma offset that would zero the gap. "
        "`R_s = |delta_kill| * SD(entry_age_tilde_s) / omega_d_s` = size of the "
        "required drift<->entry-age correlation as a fraction of the fitted drift "
        "spread; verdict keys on `R_min` (cheapest reattribution). "
        "`fan@65 %` = illustrative fan-only contribution to the entry-65 curve "
        "gap at A_max. Verdict: FRAGILE R<0.5, BORDERLINE 0.5-1.5, ROBUST >=1.5, "
        "or ROBUST(curvature) if the gauge-free bend alone excludes 0._\n")
    (OUT_ROOT / "fan_identifiability.md").write_text(
        "\n".join(lines) + "\n" + note, encoding="utf-8")
    print(f"wrote {OUT_ROOT / 'fan_identifiability.md'}\n")

    # ---- stdout ------------------------------------------------------------ #
    hdr = (f"{'cohort':>6} {'mrc':>5} | {'Dg(Amax)':>9} {'curv_sh':>7} {'cv!=0':>5} | "
           f"{'d_kill':>8} {'R_M':>6} {'R_W':>6} {'R_min':>6} | {'fan@65%':>8}  verdict")
    print(hdr)
    print("-" * len(hdr))
    for _, r in df.iterrows():
        cv = ("yes" if r["curv_excl0"] else "no") if r["has_boot"] else "n/a"
        print(f"{r['cohort']:>6} {r['mrc']:>5} | "
              f"{r['dg_end']:>+9.4f} "
              f"{(r['curv_share'] if pd.notna(r['curv_share']) else float('nan')):>7.2f} "
              f"{cv:>5} | {r['delta_kill']:>+8.5f} "
              f"{r['R_M']:>6.2f} {r['R_W']:>6.2f} {r['R_min']:>6.2f} | "
              f"{r['fan_end_pct_at_ref']:>+8.2f}  {r['verdict']}")
    print()


if __name__ == "__main__":
    main()
