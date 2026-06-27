"""One-stop number lookup for the d_i-distribution paragraph (QC, no fit).

The d_i distribution paragraph reports only **gauge-invariant** quantities (the
level/mean of d_i is not identified -- pinned by the EB prior -- so anything that
reads the absolute zero is excluded). This pulls the few numbers worth quoting,
for one cohort (default ALL), both sexes, printing on screen + to a text file.

Reportable (location-free / slope-only):
  * heterogeneity magnitude:  sqrt(omega_d2) (EB prior SD) + disatt_sd
  * shape:                    skew (right tail of decliners)
  * sex comparison:           disatt_sd / skew M-W difference (both gauge-safe)
  * entry-age gradient:       slope(d_i ~ entry age), rho  (slope is shift-invariant)

NOT reportable (flagged): frac_improver, corr(d,u) -- level/u-side gauge-anchored.

Sources:
  results/analysis/di_distribution/di_descriptors.csv      (disatt_sd, skew, +CIs)
  results/analysis/di_distribution/di_entryage_slope.csv   (slope, rho; from p02)
  results/model_selection/athlete_drift/omega_profile/{slug}/profile_{nutag}.parquet

Output -> results/analysis/di_distribution/di_dist_summary_{cohort}.md
(human-readable Markdown).

Run::

    python scripts/05_analysis/di_distribution/q99_di_dist_summary.py
    python scripts/05_analysis/di_distribution/q99_di_dist_summary.py --cohort Po10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from report_md import render_markdown, write_markdown  # noqa: E402

from marathon_decomp import registry                # noqa: E402
from marathon_decomp.config import RESULTS_DIR      # noqa: E402
from baseline_common import slices as S             # noqa: E402
import di_common as DI                              # noqa: E402

DIST = DI.OUT_ROOT                                   # analysis/di_distribution
ADROOT = RESULTS_DIR / "model_selection" / "athlete_drift"
NAN = float("nan")
SEXES = ("M", "W")


def _nutag(nu: float) -> str:
    return f"nu{nu:.2f}".replace(".", "p")


def _desc(dd: pd.DataFrame, cohort: str, group: str, floor: int, col: str) -> float:
    m = (dd.cohort == cohort) & (dd.group == group) & (dd.n_floor == floor)
    sub = dd[m]
    return float(sub[col].iloc[0]) if len(sub) else NAN


def _omega(cohort: str, sex: str, nutag: str, mrc: int) -> float:
    slug = registry.slice_slug(S.build_spec(f"{cohort}_{sex}", min_race_count=mrc))
    p = ADROOT / "omega_profile" / slug / f"profile_{nutag}.parquet"
    if not p.is_file():
        return NAN
    pr = pd.read_parquet(p)
    free = pr[pr.is_free]
    return float(free.omega_d2.iloc[0]) if len(free) else NAN


def collect(cohort: str, nu: float, floor: int, mrc: int) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    nutag = _nutag(nu)
    dd = pd.read_csv(DIST / "di_descriptors.csv")
    slope_p = DIST / "di_entryage_slope.csv"
    sl = pd.read_csv(slope_p) if slope_p.is_file() else None

    # ============ A. HETEROGENEITY (gauge-invariant) ================
    A = "A. HETEROGENEITY MAGNITUDE  (gauge-invariant)"
    for sx in SEXES:
        om = _omega(cohort, sx, nutag, mrc)
        add(A, f"sqrt(omega_d2)=EB prior SD  [{sx}]",
            f"{np.sqrt(om) * 100:.2f} %/yr" if np.isfinite(om) else "n/a",
            f"athlete_drift/omega_profile/{cohort}_{sx}.../profile_{nutag}.parquet, is_free")
    for sx in SEXES:
        sd = _desc(dd, cohort, sx, floor, "disatt_sd")
        lo = _desc(dd, cohort, sx, floor, "disatt_sd_lo")
        hi = _desc(dd, cohort, sx, floor, "disatt_sd_hi")
        add(A, f"disatt SD (n>={floor})  [{sx}]",
            f"{sd * 100:.2f} [{lo * 100:.2f}, {hi * 100:.2f}] %/yr",
            "di_descriptors.csv, disatt_sd(+CI)")
    add(A, "  -> reading", "a 1-SD athlete drifts ~1.3-1.7 %/yr off the shared curve "
        "(noise-corrected; gauge-free)", "derived")

    # ============ B. SHAPE (gauge-invariant) ========================
    B = "B. SHAPE  (gauge-invariant)"
    for sx in SEXES:
        sk = _desc(dd, cohort, sx, floor, "skew")
        lo = _desc(dd, cohort, sx, floor, "skew_lo")
        hi = _desc(dd, cohort, sx, floor, "skew_hi")
        sig = "" if (lo <= 0 <= hi) else "  (CI excludes 0)"
        add(B, f"skew (n>={floor})  [{sx}]", f"{sk:+.3f} [{lo:+.3f}, {hi:+.3f}]{sig}",
            "di_descriptors.csv, skew(+CI)")
    add(B, "  -> reading", "modestly right-skewed: a tail of rapid decliners "
        "(d>0); small in magnitude", "derived")

    # ============ C. SEX COMPARISON (gauge-safe diff) ===============
    C = "C. SEX COMPARISON  (only gauge-safe SD & skew)"
    sd_lo = _desc(dd, cohort, "M-W", floor, "disatt_sd_lo")
    sd_hi = _desc(dd, cohort, "M-W", floor, "disatt_sd_hi")
    sd_d = _desc(dd, cohort, "M-W", floor, "disatt_sd")
    add(C, f"disatt SD  M-W (n>={floor})",
        f"{sd_d * 100:+.2f} [{sd_lo * 100:+.2f}, {sd_hi * 100:+.2f}] %/yr",
        "di_descriptors.csv, group=M-W")
    sk_lo = _desc(dd, cohort, "M-W", floor, "skew_lo")
    sk_hi = _desc(dd, cohort, "M-W", floor, "skew_hi")
    sk_d = _desc(dd, cohort, "M-W", floor, "skew")
    incl0 = "CI includes 0" if (sk_lo <= 0 <= sk_hi) else "CI excludes 0"
    add(C, f"skew  M-W (n>={floor})", f"{sk_d:+.3f} [{sk_lo:+.3f}, {sk_hi:+.3f}]  ({incl0})",
        "di_descriptors.csv, group=M-W")
    add(C, "  -> reading", "trajectory spread & shape essentially equal across sexes "
        "(no robust difference)", "derived")

    # ============ D. ENTRY-AGE GRADIENT (gauge-safe slope) ==========
    D = "D. ENTRY-AGE GRADIENT  (gauge-safe; slope is shift-invariant)"
    if sl is not None:
        for sx in SEXES:
            r = sl[(sl.cohort == cohort) & (sl.sex == sx)]
            if not len(r):
                continue
            r = r.iloc[0]
            sp = r.slope_per_yr * 1000  # %/yr per 10yr entry age
            lo, hi = r.slope_lo * 1000, r.slope_hi * 1000
            rho = r.spearman_d_Ae
            add(D, f"slope(d ~ entry age)  [{sx}]",
                f"{sp:+.2f} [{lo:+.2f}, {hi:+.2f}] %/yr per 10yr   rho={rho:+.3f} "
                f"(var expl ~{rho ** 2 * 100:.1f}%)",
                f"di_entryage_slope.csv ({cohort}_{sx})")
        add(D, "  -> reading", "gamma fan absorbs most of the entry-age effect; a faint "
            "residual remains (later debut -> decliner, Stones direction), an order of "
            "magnitude below the ~1.3 %/yr individual SD -> d_i is overwhelmingly "
            "idiosyncratic", "derived")
    else:
        add(D, "(di_entryage_slope.csv missing)", "run p02_di_by_entryage.py first", "-")

    # ============ E. NOT REPORTABLE (flag) ==========================
    E = "E. NOT REPORTABLE  (gauge/prior-anchored -- do NOT quote)"
    add(E, "frac_improver", "EXCLUDED: counts d<0, but the zero is prior-pinned "
        "(moves with the gauge); not comparable across slices", "caveat")
    add(E, "corr(d, u)", "EXCLUDED: gauge-dependent via the u-side transform "
        "(u += -c*mean_i A_n); additionally era-confounded", "caveat")

    return rows


def render(cohort: str, floor: int, rows: list[tuple[str, str, str, str]]) -> str:
    return render_markdown(
        f"d_i distribution summary (cohort {cohort}, n_i>={floor})",
        rows,
        subtitle=[
            "Only GAUGE-INVARIANT quantities are reportable: the level/mean of d_i is "
            "pinned by the EB prior, so shape (spread, skew) and slopes only.",
            "d_i < 0 = improver; frac_improver and corr(d,u) are level-anchored and "
            "explicitly flagged NOT reportable.",
        ],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="ALL", choices=list(DI.POPS))
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--floor", type=int, default=5)
    ap.add_argument("--mrc", type=int, default=2)
    args = ap.parse_args()

    rows = collect(args.cohort, args.nu, args.floor, args.mrc)
    report = render(args.cohort, args.floor, rows)
    print(report)

    out = write_markdown(DIST / f"di_dist_summary_{args.cohort}.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
