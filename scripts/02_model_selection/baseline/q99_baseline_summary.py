"""One-stop number lookup for the L2-vs-Student-t / nu=8 paragraph (QC, no fit).

Reads the CSVs the baseline selection + residual diagnostics already wrote and
pulls the exact numbers the paper quotes, for ONE slice (default ALL_M), printing
them on screen and to a human-readable text file. Each number is annotated with
its source file + column so it is trivial to verify and cite.

Sources (under results/model_selection/baseline/):
  * {slug}/param_sensitivity.csv : sigma2 per nu; v-vs-L2 spearman/pearson/|dv|.
  * {slug}/nu_selection.csv       : held-out CV log-density per cell per nu.
  * {slug}/selected_nu.csv        : grid-best / Brent-refined / selected nu.
  * nu_decision.csv               : per-slice 1-SE nu interval + recommended nu.
  * qq_plot/residual_diagnostics.csv : sigma2, skew, excess kurtosis,
                                       kurtosis-implied nu, QQ R^2 (Normal / t).

Output -> results/model_selection/baseline/baseline_summary_{slug}.md
(human-readable Markdown the paper cites directly).

Run::

    python scripts/02_model_selection/baseline/q02_baseline_summary.py
    python scripts/02_model_selection/baseline/q02_baseline_summary.py --slice ALL_B
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from baseline_common import slices as S  # noqa: E402
from report_md import render_markdown, write_markdown  # noqa: E402

ROOT = RESULTS_DIR / "model_selection" / "baseline"
INF = float("inf")


def _pick_solver(df: pd.DataFrame) -> str:
    """Prefer the 'anderson' solver (ALS/Anderson share the fixed point)."""
    sv = set(df["solver"].unique())
    return "anderson" if "anderson" in sv else sorted(sv)[0]


def _at(df: pd.DataFrame, nu: float, col: str) -> float:
    """Value of `col` at degrees-of-freedom `nu` (inf matched specially)."""
    m = ~np.isfinite(df["nu"]) if not np.isfinite(nu) else np.isclose(df["nu"], nu)
    sub = df[m]
    return float(sub[col].iloc[0]) if len(sub) else float("nan")


def collect(slug: str) -> list[tuple[str, str, str, str]]:
    """Return (section, label, value-string, source) rows for the report."""
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    # ── param_sensitivity: sigma2 + v-vs-L2 ──────────────────────────
    ps = pd.read_csv(ROOT / slug / "param_sensitivity.csv")
    sv = _pick_solver(ps)
    ps = ps[ps.solver == sv]
    s2_l2 = _at(ps, INF, "sigma2")
    s2_t8 = _at(ps, 8.0, "sigma2")
    spear8 = _at(ps, 8.0, "spearman_vs_L2")
    pear8 = _at(ps, 8.0, "pearson_vs_L2")
    maxdv8 = _at(ps, 8.0, "max_abs_dv")
    src_ps = f"{slug}/param_sensitivity.csv (solver={sv})"

    W = "WHY STUDENT-T"
    add(W, "sigma^2  (L2, nu=inf)", f"{s2_l2:.4e}", f"{src_ps}, nu=inf, col sigma2")
    add(W, "sigma^2  (Student-t, nu=8)", f"{s2_t8:.4e}", f"{src_ps}, nu=8, col sigma2")
    add(W, "  -> L2 variance inflation",
        f"{s2_l2 / s2_t8:.3f}x  ({(s2_l2 / s2_t8 - 1) * 100:.0f}% larger)", "derived")

    # ── residual_diagnostics: QQ R^2, skew, kurtosis ─────────────────
    rd = pd.read_csv(ROOT / "qq_plot" / "residual_diagnostics.csv")
    rd = rd[rd.slice == slug]
    src_rd = "qq_plot/residual_diagnostics.csv"

    def rd_val(model: str, col: str) -> float:
        sub = rd[rd.model == model]
        return float(sub[col].iloc[0]) if len(sub) else float("nan")

    L2M, T8M = "baseline (u+v), L2", "baseline (u+v), t8"
    add(W, "QQ R^2  baseline L2 vs Normal", f"{rd_val(L2M, 'qq_r2_normal'):.4f}",
        f"{src_rd}, model='{L2M}', qq_r2_normal")
    add(W, "QQ R^2  baseline t8 vs t(8)", f"{rd_val(T8M, 'qq_r2_ref'):.4f}",
        f"{src_rd}, model='{T8M}', qq_r2_ref")
    add(W, "skewness (baseline t8 resid)", f"{rd_val(T8M, 'skewness'):.3f}",
        f"{src_rd}, model='{T8M}', skewness")
    add(W, "excess kurtosis (baseline t8)", f"{rd_val(T8M, 'excess_kurtosis'):.3f}",
        f"{src_rd}, model='{T8M}', excess_kurtosis")
    add(W, "nu implied by kurtosis", f"{rd_val(T8M, 'nu_implied_by_kurtosis'):.3f}",
        f"{src_rd}, nu_implied_by_kurtosis  (=6/exk+4)")

    # ── in-sample fit + complexity: AIC / BIC / effective d.o.f. ──────
    # AIC = 2k - 2*loglik  =>  effective d.o.f. k = AIC/2 + loglik (no N needed).
    # nu/sigma2 are fixed hyperparameters (nu chosen on held-out CV), so k is
    # IDENTICAL for L2 and t8: the whole AIC/BIC drop is the likelihood gain, not
    # a complexity trade. Both log-liks are FULLY normalized densities of the same
    # response y (Gaussian vs Student-t constants), so they compare across nu.
    ll_l2, ll_t8 = _at(ps, INF, "full_loglik"), _at(ps, 8.0, "full_loglik")
    aic_l2, aic_t8 = _at(ps, INF, "full_aic"), _at(ps, 8.0, "full_aic")
    bic_l2, bic_t8 = _at(ps, INF, "full_bic"), _at(ps, 8.0, "full_bic")
    dof_l2, dof_t8 = aic_l2 / 2 + ll_l2, aic_t8 / 2 + ll_t8

    F = "IN-SAMPLE FIT & COMPLEXITY (AIC / BIC / eff. d.o.f.)"
    add(F, "log-lik   (L2, nu=inf)", f"{ll_l2:,.1f}", f"{src_ps}, nu=inf, full_loglik")
    add(F, "log-lik   (Student-t, nu=8)", f"{ll_t8:,.1f}", f"{src_ps}, nu=8, full_loglik")
    add(F, "effective d.o.f.  (L2 / t8)",
        f"{dof_l2:,.0f} / {dof_t8:,.0f}  (=AIC/2+loglik; equal, nu not counted)", "derived")
    add(F, "AIC       (L2 / t8)", f"{aic_l2:,.0f} / {aic_t8:,.0f}", f"{src_ps}, full_aic")
    add(F, "  -> dAIC (t8 - L2)", f"{aic_t8 - aic_l2:,.0f}  (t8 lower = better)", "derived")
    add(F, "BIC       (L2 / t8)", f"{bic_l2:,.0f} / {bic_t8:,.0f}", f"{src_ps}, full_bic")
    add(F, "  -> dBIC (t8 - L2)", f"{bic_t8 - bic_l2:,.0f}  (t8 lower; dof unchanged)", "derived")

    # ── nu_selection: held-out CV density per cell ───────────────────
    ns = pd.read_csv(ROOT / slug / "nu_selection.csv")
    ns = ns[(ns.solver == _pick_solver(ns)) & (ns.source == "grid")]
    cv_l2 = _at(ns, INF, "cv_per_cell")
    cv_8 = _at(ns, 8.0, "cv_per_cell")
    fin = ns[np.isfinite(ns["nu"])]
    cv_max_nu = float(fin.loc[fin.cv_per_cell.idxmax(), "nu"])
    src_ns = f"{slug}/nu_selection.csv (source=grid)"

    C = "CHOOSING NU (out-of-sample CV)"
    add(C, "held-out logdens/cell  (L2, nu=inf)", f"{cv_l2:.4f}", f"{src_ns}, nu=inf, cv_per_cell")
    add(C, "held-out logdens/cell  (nu=8)", f"{cv_8:.4f}", f"{src_ns}, nu=8, cv_per_cell")
    add(C, "  -> CV gain (nats/cell)", f"+{cv_8 - cv_l2:.4f}", "derived")
    add(C, "  => term earns its keep?",
        f"in-sample BIC {bic_t8 - bic_l2:,.0f} AND held-out CV +{cv_8 - cv_l2:.4f}  "
        "(both improve -> not in-sample overfit)", "derived")
    add(C, "  grid argmax nu", f"{cv_max_nu:g}", f"{src_ns}, argmax cv_per_cell")

    # ── selected_nu: grid / Brent / selected ─────────────────────────
    sn = pd.read_csv(ROOT / slug / "selected_nu.csv").iloc[0]
    add(C, "grid-best nu", f"{sn['grid_best_nu']:g}", f"{slug}/selected_nu.csv, grid_best_nu")
    add(C, "Brent-refined nu*", f"{sn['brent_nu']:.3f}", f"{slug}/selected_nu.csv, brent_nu")
    add(C, "selected (shared) nu", f"{sn['selected_nu']:g}", f"{slug}/selected_nu.csv, selected_nu")

    # ── nu_decision: per-slice 1-SE band + recommended shared nu ─────
    nd = pd.read_csv(ROOT / "nu_decision.csv")
    this = nd[nd.slice == slug]
    if len(this):
        lo, hi = float(this.nu_1se_lo.iloc[0]), float(this.nu_1se_hi.iloc[0])
        add(C, "1-SE acceptable nu interval", f"[{lo:g}, {hi:g}]",
            "nu_decision.csv, this slice, nu_1se_lo/hi")
    rec = nd[nd.slice == "__RECOMMENDED__"]
    if len(rec):
        add(C, "recommended shared nu (across slices)", f"{float(rec.nu_argmax.iloc[0]):g}",
            "nu_decision.csv, slice=__RECOMMENDED__, nu_argmax")

    add(C, "v ranking stability vs L2 @nu=8", f"Spearman {spear8:.4f}, Pearson {pear8:.4f}",
        f"{src_ps}, nu=8, spearman_vs_L2 / pearson_vs_L2")
    add(C, "  max |dv| vs L2 (log-time)", f"{maxdv8:.4f}", f"{src_ps}, nu=8, max_abs_dv")

    return rows


def render(slug: str, rows: list[tuple[str, str, str, str]]) -> str:
    return render_markdown(
        f"Baseline noise model: Student-t selection ({slug})",
        rows,
        subtitle=[
            "Rank-1 baseline log t = u_i + v_j; ALS/Anderson share the fixed point.",
            "Degrees of freedom nu fixed by held-out predictive density (the "
            "in-sample fit always prefers heavier tails); a single shared nu=8 is used.",
        ],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default="ALL_B")
    ap.add_argument("--mrc", type=int, default=None)
    args = ap.parse_args()

    slug = S.slug(S.build_spec(args.slice, min_race_count=args.mrc))
    rows = collect(slug)
    report = render(slug, rows)
    print(report)

    out = write_markdown(ROOT / f"baseline_summary_{slug}.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
