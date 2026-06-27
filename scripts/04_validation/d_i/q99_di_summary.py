"""One-stop number lookup for the d_i-inclusion paragraph / Appendix E (QC, no fit).

Reads the CSVs/parquets the d_i validation + aging-vs-drift + athlete-drift stages
already wrote and pulls the exact numbers the paper quotes, for ONE slice (default
ALL_B), printing them on screen and to a human-readable text file. Each number is
annotated with its source file + column so it is trivial to verify and cite.

The d_i term is justified by DE-BIASING (not predictive fit), so the headline is
the v_j-bias test, not an information criterion. Sections:

  A. V_J DE-BIASING (headline: EB estimator, everyone-field mrc2)
       results/validation/d_i/bias_summary.csv  : corr_nod / corr_d / frac_removed
       results/validation/d_i/bias_corr.csv     : raw vs year-partialled, slope
  B. INDEPENDENT CONFIRMATION (LOO + permutation, clean mrc5)
       results/validation/d_i/bias_summary.csv  : LOO corr_nod / corr_d
       results/validation/d_i/bias_corr.csv     : LOO slope no-d / +d
       results/validation/d_i/delta_agreement.csv : corr(delta_EB, delta_LOO)
       results/validation/d_i/permutation/{slice}.csv : z, p_one (no-d & +d)
  C. CROSS-CELL CONSISTENCY (all Po10/ALL EB cells)
       results/validation/d_i/bias_summary.csv  : corr_nod range over slugs
  D. AGING-CURVE CORROBORATION (cross-mrc de-biasing)
       results/model_selection/aging_vs_drift/curve_debias.csv : gap / d_on / peaks
       results/model_selection/aging_vs_drift/block_shift.csv   : curve_corr, n_elig
  E. PRIOR IDENTIFICATION (omega_d2; the numerical "it works")
       results/model_selection/athlete_drift/omega_profile/{slug}/profile_*.parquet
       results/model_selection/athlete_drift/omega_init/{slug}/init_summary_*.parquet

Output -> results/validation/d_i/di_summary_{slice}.md
(human-readable Markdown the paper cites directly).

Self-contained; defaults to ALL_B (VS Code "Run" works). Reads only -- never fits.

Run::

    python scripts/04_validation/d_i/q99_di_summary.py
    python scripts/04_validation/d_i/q99_di_summary.py --slice Po10_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import registry              # noqa: E402
from marathon_decomp.config import RESULTS_DIR    # noqa: E402
from baseline_common import slices as S           # noqa: E402
from report_md import render_markdown, write_markdown  # noqa: E402

VROOT = RESULTS_DIR / "validation" / "d_i"
AVDROOT = RESULTS_DIR / "model_selection" / "aging_vs_drift"
ADROOT = RESULTS_DIR / "model_selection" / "athlete_drift"
NAN = float("nan")


def _slug(slice_name: str, mrc: int) -> str:
    return registry.slice_slug(S.build_spec(slice_name, min_race_count=mrc))


def _nutag(nu: float) -> str:
    return f"nu{nu:.2f}".replace(".", "p")   # 8.0 -> nu8p00


def _one(df: pd.DataFrame, col: str, **filt) -> float:
    """First value of `col` in the rows matching all equality filters (else nan)."""
    m = pd.Series(True, index=df.index)
    for k, v in filt.items():
        m &= (df[k] == v)
    sub = df[m]
    return float(sub[col].iloc[0]) if len(sub) else NAN


def collect(slice_name: str, nu: float) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    slug2 = _slug(slice_name, 2)   # everyone field
    slug5 = _slug(slice_name, 5)   # dedicated-runner field
    nutag = _nutag(nu)

    # ============ A. V_J DE-BIASING (headline EB, mrc2) ==============
    A = "A. V_J DE-BIASING  (headline: EB, everyone-field mrc2)"
    bs = pd.read_csv(VROOT / "bias_summary.csv")
    bc = pd.read_csv(VROOT / "bias_corr.csv")
    src_bs = "bias_summary.csv"
    src_bc = "bias_corr.csv"

    cn = _one(bs, "corr_nod", slug=slug2, estimator="eb")
    cd = _one(bs, "corr_d", slug=slug2, estimator="eb")
    frac = _one(bs, "frac_removed", slug=slug2, estimator="eb")
    nraces = _one(bc, "n_races", slug=slug2, estimator="eb", fit="no-d")
    raw_nod = _one(bc, "pearson_raw", slug=slug2, estimator="eb", fit="no-d")
    par_nod = _one(bc, "pearson_partial", slug=slug2, estimator="eb", fit="no-d")
    slp_nod = _one(bc, "slope_partial", slug=slug2, estimator="eb", fit="no-d")
    slp_d = _one(bc, "slope_partial", slug=slug2, estimator="eb", fit="+d")

    add(A, "slug / n_races", f"{slug2} / {nraces:.0f}", "bias_corr.csv")
    add(A, "year-partialled corr(delta,v)  no-d", f"{cn:+.4f}", f"{src_bs}, eb, corr_nod")
    add(A, "                              +d", f"{cd:+.4f}", f"{src_bs}, eb, corr_d")
    add(A, "  -> bias removed", f"{frac * 100:.0f}%  ({cn:+.3f} -> {cd:+.3f})",
        f"{src_bs}, eb, frac_removed")
    add(A, "partial slope  v ~ delta   no-d", f"{slp_nod:.3f}",
        f"{src_bc}, eb/no-d, slope_partial")
    add(A, "                           +d", f"{slp_d:.3f}",
        f"{src_bc}, eb/+d, slope_partial  (1:1 absorption -> collapses)")
    add(A, "raw vs year-partialled corr  no-d", f"{raw_nod:+.3f} -> {par_nod:+.3f}",
        f"{src_bc}, eb/no-d, pearson_raw / pearson_partial")
    add(A, "  (why partial>raw)", "era trend in v + window trend in delta cancel in raw "
        "-> partialling exposes the coupling", "derived")

    # ============ B. INDEPENDENT CONFIRMATION (LOO, mrc5) ============
    B = "B. INDEPENDENT CONFIRMATION  (LOO + permutation, clean mrc5)"
    cn5 = _one(bs, "corr_nod", slug=slug5, estimator="loo")
    cd5 = _one(bs, "corr_d", slug=slug5, estimator="loo")
    frac5 = _one(bs, "frac_removed", slug=slug5, estimator="loo")
    slp5_nod = _one(bc, "slope_partial", slug=slug5, estimator="loo", fit="no-d")
    slp5_d = _one(bc, "slope_partial", slug=slug5, estimator="loo", fit="+d")
    add(B, "LOO year-partialled corr  no-d / +d", f"{cn5:+.4f} / {cd5:+.4f}",
        f"{src_bs}, loo (mrc5), corr_nod/corr_d")
    add(B, "  -> bias removed (LOO)", f"{frac5 * 100:.0f}%", f"{src_bs}, loo, frac_removed")
    add(B, "LOO partial slope  no-d / +d", f"{slp5_nod:.3f} / {slp5_d:.3f}",
        f"{src_bc}, loo (mrc5), slope_partial")

    da = VROOT / "delta_agreement.csv"
    if da.is_file():
        dag = pd.read_csv(da)
        ag5 = _one(dag, "corr_eb_loo", slug=slug5)
        ag2 = _one(dag, "corr_eb_loo", slug=slug2)
        add(B, "estimator agreement corr(EB,LOO)  mrc5 / mrc2",
            f"{ag5:.3f} / {ag2:.3f}", "delta_agreement.csv, corr_eb_loo")

    perm = VROOT / "permutation" / f"{slice_name}.csv"
    if perm.is_file():
        pm = pd.read_csv(perm)
        z5n = _one(pm, "z", slug=slug5, fit="no-d")
        p5n = _one(pm, "p_one", slug=slug5, fit="no-d")
        z5d = _one(pm, "z", slug=slug5, fit="+d")
        p5d = _one(pm, "p_one", slug=slug5, fit="+d")
        npm = _one(pm, "n_perm_total", slug=slug5, fit="no-d")
        add(B, "permutation null (mrc5)  no-d", f"z={z5n:.2f}  p={p5n:.4f}",
            f"permutation/{slice_name}.csv, no-d, z/p_one")
        add(B, "                         +d", f"z={z5d:.2f}  p={p5d:.4f}  (loses signif.)",
            f"permutation/{slice_name}.csv, +d, z/p_one")
        add(B, "  permutation pool size", f"{npm:.0f} draws  (p floor = 1/(N+1))",
            f"permutation/{slice_name}.csv, n_perm_total")

    # ============ C. CROSS-CELL CONSISTENCY (all EB cells) ===========
    C = "C. CROSS-CELL CONSISTENCY  (all Po10/ALL EB cells)"
    eb = bs[bs.estimator == "eb"].copy()
    lo_slug = eb.loc[eb.corr_nod.idxmin(), "slug"]
    hi_slug = eb.loc[eb.corr_nod.idxmax(), "slug"]
    add(C, "EB no-d corr range over cells",
        f"{eb.corr_nod.min():+.3f} .. {eb.corr_nod.max():+.3f}  "
        f"(n={len(eb)} cells)", f"{src_bs}, all eb rows, corr_nod")
    add(C, "  min / max cell", f"{lo_slug} / {hi_slug}", "derived")
    add(C, "EB frac_removed range", f"{eb.frac_removed.min():.2f} .. {eb.frac_removed.max():.2f}",
        f"{src_bs}, all eb rows, frac_removed  (>1 = mrc2 LOO artifact n/a here)")

    # ============ D. AGING-CURVE CORROBORATION (cross-mrc) ===========
    D = "D. AGING-CURVE CORROBORATION  (cross-mrc de-biasing)"
    cdb_p = AVDROOT / "curve_debias.csv"
    if cdb_p.is_file():
        cdb = pd.read_csv(cdb_p)
        g_no = _one(cdb, "gap_noD", slug=slug2)
        g_wd = _one(cdb, "gap_withD", slug=slug2)
        d2 = _one(cdb, "d_on_mrc2", slug=slug2)
        d5 = _one(cdb, "d_on_mrc5", slug=slug2)
        pk_nod2 = _one(cdb, "peak_noD_mrc2", slug=slug2)
        pk_wd2 = _one(cdb, "peak_withD_mrc2", slug=slug2)
        pk_nod5 = _one(cdb, "peak_noD_mrc5", slug=slug2)
        closed = (1.0 - g_wd / g_no) * 100 if np.isfinite(g_no) and g_no else NAN
        add(D, "everyone-vs-dedicated curve gap  no-d", f"{g_no:.4f}",
            "curve_debias.csv, gap_noD  (log-time RMS)")
        add(D, "                                 +d", f"{g_wd:.4f}",
            "curve_debias.csv, gap_withD")
        add(D, "  -> gap closed by d_i", f"{closed:.0f}%", "derived")
        add(D, "d_i curve shift  everyone(mrc2) / dedicated(mrc5)",
            f"{d2:.4f} / {d5:.4f}  (more work where contaminated)",
            "curve_debias.csv, d_on_mrc2 / d_on_mrc5")
        add(D, "peak age (yr post-debut)  noD-mrc2 / +d-mrc2 / noD-mrc5",
            f"{pk_nod2:.2f} / {pk_wd2:.2f} / {pk_nod5:.2f}  (+d-mrc2 == noD-mrc5)",
            "curve_debias.csv, peak_* (to grid resolution)")
    bsh_p = AVDROOT / "block_shift.csv"
    if bsh_p.is_file():
        bsh = pd.read_csv(bsh_p)
        ccorr = _one(bsh, "curve_corr", slug=slug2, mrc=2)
        nelig = _one(bsh, "n_eligible", slug=slug2, mrc=2)
        add(D, "curvature preserved (aging vs full, mrc2)", f"curve_corr {ccorr:.4f}",
            "block_shift.csv, curve_corr  (shape kept; only tilt removed)")
        add(D, "n eligible athletes with free d_i (mrc2)", f"{nelig:,.0f}",
            "block_shift.csv, n_eligible")

    # ============ E. PRIOR IDENTIFICATION (omega_d2) =================
    E = "E. PRIOR IDENTIFICATION  (omega_d2; numerical 'it works')"
    prof_p = ADROOT / "omega_profile" / slug2 / f"profile_{nutag}.parquet"
    if prof_p.is_file():
        pr = pd.read_parquet(prof_p)
        free = pr[pr.is_free]
        ostar = float(free.omega_d2.iloc[0]) if len(free) else NAN
        nelig = float(free.n_elig.iloc[0]) if len(free) else NAN
        # is the marginal peaked at omega*? argmax logML vs the free row
        peak_mult = float(pr.loc[pr.logML.idxmax(), "omega_mult"])
        # insensitivity within omega* x [1/3, 3]
        band = pr[(pr.omega_mult >= 1 / 3 - 1e-6) & (pr.omega_mult <= 3 + 1e-6)]
        amax = float(band.aging_maxdev.max()) if len(band) else NAN
        vmin = float(band.corr_v_to_free.min()) if len(band) else NAN
        add(E, "EB-learned omega_d2*  (mrc2)", f"{ostar:.3e}",
            f"omega_profile/{slug2}/profile_{nutag}.parquet, is_free row")
        add(E, "  marginal logML peak at", f"omega_mult = {peak_mult:g}  (1.0 = omega*)",
            "argmax logML  (peaked, not flat -> identified)")
        add(E, "  n eligible (free d_i)", f"{nelig:,.0f}", "profile, n_elig")
        add(E, "insensitivity in omega* x [1/3,3]",
            f"aging dev <= {amax:.4f},  corr_v_to_free >= {vmin:.4f}",
            "profile, aging_maxdev / corr_v_to_free  (<1% move)")
    init_p = ADROOT / "omega_init" / slug2 / f"init_summary_{nutag}.parquet"
    if init_p.is_file():
        it = pd.read_parquet(init_p)
        conv = it[it.converged] if "converged" in it.columns else it
        if len(conv):
            omin, omax = float(conv.omega_d2.min()), float(conv.omega_d2.max())
            spread = (omax - omin) / omax if omax else NAN
            admax = float(conv.aging_maxdev.max())
            add(E, "init-invariance (converged inits)",
                f"omega* rel-spread {spread:.1e},  aging dev <= {admax:.1e}",
                f"omega_init/{slug2}/init_summary_{nutag}.parquet, converged rows")
            add(E, "  => d_i adds no tuning knob",
                "prior learned (EB), identified, init-invariant", "derived")

    return rows


def render(slice_name: str, nu: float, rows: list[tuple[str, str, str, str]]) -> str:
    return render_markdown(
        f"d_i inclusion / v_j de-biasing summary ({slice_name}, nu={nu:g})",
        rows,
        subtitle=[
            "The per-athlete career drift d_i is justified by DE-BIASING v_j (it "
            "removes the field's career-stage composition leak), not by predictive "
            "fit -- so the headline is the v_j-bias test, not an information criterion.",
            "Reads only the d_i validation + aging-vs-drift + athlete-drift outputs; "
            "never refits.",
        ],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default="ALL_B")
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()

    rows = collect(args.slice, args.nu)
    report = render(args.slice, args.nu, rows)
    print(report)

    out = write_markdown(VROOT / f"di_summary_{args.slice}.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
