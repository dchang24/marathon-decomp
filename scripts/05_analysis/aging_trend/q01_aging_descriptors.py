"""Scalar descriptors of the entry-age aging curves (no refit).

For every (slice x sex x min-race-count x entry age) the gauged curve
f(A_n) (see `p01_aging_curve.reconstruct`, beta=0 APC gauge, A_n in [0, a_max])
is reduced to a few interpretable numbers:

  * peak_age        career age A_n of best performance (argmin f); chronological
                    peak = entry age + peak_age.
  * peak_improve    largest time reduction vs debut: -min(f) in log units, and
                    (1 - exp(min f)) * 100 as a %.
  * improve_span    career-age width of the improvement window: first A_n > peak
                    where f returns to 0 (right-censored at a_max if still < 0).
  * end_change      f at A_n = a_max (net change vs debut), log + %.
  * decline         end_change - min(f): rise from the trough to the horizon (%).

Bootstrap (any fit carrying an athlete-weight bootstrap -- now both mrc2 and mrc5)
gives a 2.5/97.5 CI per descriptor by recomputing it on each replicate curve
(theta/gamma vary; the gauge tilt c_beta is held at its point value, matching the
p01 band).

Outputs (results/analysis/aging_trend/):
    aging_descriptors.parquet   tidy, all columns + CIs
    aging_descriptors.csv       same, flat + rounded (spreadsheet-friendly)
    aging_descriptors.md        readable table, CIs inline (paper-friendly)
    + a readable table to stdout.

Run::

    python scripts/05_analysis/aging_trend/q01_aging_descriptors.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir (p01)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from p01_aging_curve import OUT_ROOT, reconstruct  # noqa: E402

COHORTS = ["ALL", "Po10"]
MRCS = ["mrc2", "mrc5"]
DEFAULT_ENTRY_AGES = (35, 45, 55, 65)
NO_IMPROVE = -1e-6  # min f above this -> treat as "no improvement"


def _curve_descriptors(A: np.ndarray, y: np.ndarray) -> dict:
    """Descriptors of one gauged curve y = f(A), with f(0)=0."""
    j = int(np.argmin(y))
    fmin = float(y[j])
    improved = fmin < NO_IMPROVE

    # improvement span: first crossing back up through 0 after the trough
    span = 0.0
    censored = False
    if improved:
        span = np.nan
        for i in range(j, len(y) - 1):
            if y[i] < 0 <= y[i + 1]:
                t = -y[i] / (y[i + 1] - y[i])
                span = float(A[i] + t * (A[i + 1] - A[i]))
                break
        if np.isnan(span):           # still below 0 at the horizon
            span, censored = float(A[-1]), True

    end = float(y[-1])
    return {
        "peak_age": float(A[j]) if improved else 0.0,
        "peak_improve_log": -fmin if improved else 0.0,
        "peak_improve_pct": (1.0 - np.exp(fmin)) * 100.0 if improved else 0.0,
        "improve_span": span,
        "span_censored": censored,
        "end_log": end,
        "end_pct": (np.exp(end) - 1.0) * 100.0,
        "decline_pct": (np.exp(end) - np.exp(fmin)) * 100.0,  # trough -> horizon
    }


def _ci(vals: np.ndarray, lo=2.5, hi=97.5) -> tuple[float, float]:
    return float(np.percentile(vals, lo)), float(np.percentile(vals, hi))


def descriptors_for(fan, entry_age: float, a_max: float, n_grid: int = 400) -> dict:
    """Point descriptors + bootstrap CIs for one fan at one entry age."""
    A = np.linspace(0.0, a_max, n_grid)
    # gauged curve(s) on a fine grid: rebuild from stored basis pieces via at_entry's
    # grid (reconstruct already fixed A_grid); resample by interpolation to A.
    gp, _, _ = fan.at_entry(entry_age)
    yp = np.interp(A, fan.A_grid, gp)
    d = _curve_descriptors(A, yp)

    out = {k: v for k, v in d.items()}
    if fan.has_boot:
        a = entry_age - fan.mean_Ae
        boot = (fan.raw_boot + fan.c_beta * fan.A_grid[None, :] + a * fan.fan_boot)
        rows = [_curve_descriptors(A, np.interp(A, fan.A_grid, b)) for b in boot]
        bdf = pd.DataFrame(rows)
        for key in ("peak_age", "peak_improve_pct", "improve_span", "end_pct", "decline_pct"):
            lo, hi = _ci(bdf[key].to_numpy())
            out[f"{key}_lo"], out[f"{key}_hi"] = lo, hi
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="full")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--entry-ages", type=float, nargs="+", default=list(DEFAULT_ENTRY_AGES))
    ap.add_argument("--a-max", type=float, default=10.0)
    args = ap.parse_args()

    recs = []
    for cohort in COHORTS:
        for sx in ("M", "W"):
            for mrc in MRCS:
                slug = f"{cohort}_{sx}_14-25_{mrc}"
                fan = reconstruct(slug, model=args.model, nutag=args.nutag,
                                  data_version=args.data_version, a_max=args.a_max,
                                  allow_point_only=True)
                if fan is None:
                    continue
                for ea in args.entry_ages:
                    d = descriptors_for(fan, ea, args.a_max)
                    d.update(cohort=cohort, sex=sx, mrc=mrc, entry_age=int(ea),
                             has_boot=fan.has_boot,
                             chrono_peak_age=ea + d["peak_age"])
                    recs.append(d)

    df = pd.DataFrame(recs)
    df = df.sort_values(["cohort", "sex", "entry_age", "mrc"]).reset_index(drop=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out = OUT_ROOT / "aging_descriptors.parquet"
    df.to_parquet(out, index=False)
    print(f"wrote {out}  ({len(df)} rows)")

    # ---- flat CSV (every column, rounded) -------------------------------- #
    order = ["cohort", "sex", "mrc", "entry_age", "has_boot",
             "peak_age", "chrono_peak_age", "peak_improve_log", "peak_improve_pct",
             "improve_span", "span_censored", "end_log", "end_pct", "decline_pct",
             "peak_age_lo", "peak_age_hi", "peak_improve_pct_lo", "peak_improve_pct_hi",
             "improve_span_lo", "improve_span_hi", "end_pct_lo", "end_pct_hi",
             "decline_pct_lo", "decline_pct_hi"]
    order = [c for c in order if c in df.columns]
    df_csv = df[order].copy()
    num = df_csv.select_dtypes(include="number").columns
    df_csv[num] = df_csv[num].round(3)
    csv_out = OUT_ROOT / "aging_descriptors.csv"
    df_csv.to_csv(csv_out, index=False)
    print(f"wrote {csv_out}")

    # ---- readable Markdown table (CIs inline for mrc5) ------------------- #
    def cell(r, key, dp=2):
        s = f"{r[key]:.{dp}f}"
        lo, hi = r.get(f"{key}_lo"), r.get(f"{key}_hi")
        if pd.notna(lo):
            s += f" [{lo:.1f}, {hi:.1f}]"
        return s

    cols = ["Cohort", "Sex", "Min races", "Entry age", "Peak age (A_n)",
            "Chrono peak", "Peak improve %", "Improve span (yr)",
            "End % (A_n=10)", "Decline % (trough->10)"]
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    any_censored = False
    for _, r in df.iterrows():
        span = cell(r, "improve_span")
        if bool(r["span_censored"]):
            span += " *"
            any_censored = True
        lines.append("| " + " | ".join([
            r["cohort"], r["sex"], r["mrc"][3:], str(int(r["entry_age"])),
            cell(r, "peak_age"), f"{r['chrono_peak_age']:.1f}",
            cell(r, "peak_improve_pct"), span,
            cell(r, "end_pct"), cell(r, "decline_pct"),
        ]) + " |")
    note = ("\n_Values are point estimates; `[lo, hi]` = 95% athlete-bootstrap CI "
            "(present wherever the fit carries a bootstrap -- now both mrc2 and "
            "mrc5). beta=0 APC gauge. "
            "`%` = (exp(f) - 1) x 100, time vs the athlete's own debut.")
    if any_censored:
        note += " `*` improve span right-censored: curve still below 0 at A_n=10."
    note += "_\n"
    md_out = OUT_ROOT / "aging_descriptors.md"
    md_out.write_text("# Aging-curve descriptors\n\n" + "\n".join(lines) + "\n" + note,
                      encoding="utf-8")
    print(f"wrote {md_out}\n")

    # readable table
    def fmt_ci(row, key):
        lo, hi = row.get(f"{key}_lo"), row.get(f"{key}_hi")
        if pd.isna(lo):
            return "       -      "
        return f"[{lo:5.1f},{hi:6.1f}]"

    for cohort in COHORTS:
        print(f"===== {cohort} =====")
        sub = df[df.cohort == cohort]
        hdr = (f"{'sex':>3} {'mrc':>5} {'entry':>5} | {'peak_age':>8} {'chr_peak':>8} "
               f"{'improve%':>8} {'span':>6} {'cens':>4} | {'end%':>7} {'decline%':>8}  "
               f"{'improve%_CI':>18}")
        print(hdr)
        for _, r in sub.sort_values(["sex", "entry_age", "mrc"]).iterrows():
            print(f"{r.sex:>3} {r.mrc:>5} {r.entry_age:>5d} | "
                  f"{r.peak_age:8.2f} {r.chrono_peak_age:8.1f} "
                  f"{r.peak_improve_pct:8.2f} {r.improve_span:6.2f} "
                  f"{str(bool(r.span_censored)):>4} | "
                  f"{r.end_pct:7.2f} {r.decline_pct:8.2f}  "
                  f"{fmt_ci(r, 'peak_improve_pct'):>18}")
        print()


if __name__ == "__main__":
    main()
