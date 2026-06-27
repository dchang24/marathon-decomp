"""One-stop number lookup for the aging-form-selection paragraph (QC, no fit).

Reads the rollup CSVs the aging form-selection grid + CV + curve-metric scripts
already wrote and pulls the exact numbers the paper quotes, for ONE slice
(default ALL_B, the headline cohort), printing them on screen and to a
human-readable text file. Each number is annotated with its source file +
column so it is trivial to verify and cite.

The decision hierarchy (see this dir's README): held-out CV log-density/cell is
primary, curve plausibility (sane peak age, no tail wiggle) is a guardrail, BIC
is a tie-break. Drift d_i is OFF during form selection so the global curve
isolates the aggregate aging signal.

Sources (under results/model_selection/aging/):
  * form_compare.csv        : per-form loglik/AIC/BIC + CV/cell, dCV_vs_best,
                              paired_z across folds; incl. the [baseline no-aging]
                              row scored on the SAME folds (directly paired).
  * best_form_all.csv       : per-slice CV-best form + its CV/cell.
  * curve_metrics_all.csv   : curve-shape guardrail (peak age, tail sign-changes /
                              non-monotone flag) at the mean entry age.
  * v_vs_reference_all.csv  : v agreement (Pearson/Spearman/|dv|) of each form
                              vs the poly2-gscalar reference -> v robust to form.

Output -> results/model_selection/aging/aging_summary_{slug}.md
(human-readable Markdown the paper cites directly).

Self-contained; defaults to ALL_B (VS Code "Run" works).

Run::

    python scripts/02_model_selection/aging/q99_aging_summary.py
    python scripts/02_model_selection/aging/q99_aging_summary.py --slice ALL_M
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

ROOT = RESULTS_DIR / "model_selection" / "aging"
NU = 8.0

# Production aging form, and the reference points it is judged against.
CHOSEN = "spline4-gvarying"   # the production form (natural cubic spline, 4 knots, varying gamma)
GSCALAR = "spline4-gscalar"   # same basis, single (scalar) entry-age coef
GOFF = "spline4-goff"         # same basis, no entry-age term at all
HIPOLY = "poly5-gvarying"     # lowest poly degree that matches spline flexibility -> tail wiggle
BASELINE = "[baseline no-aging]"
V_REF = "poly2-gscalar"       # reference form for the v-agreement table


def _row(df: pd.DataFrame, **eq) -> pd.Series:
    """Single matching row (NaN-safe), else an all-NaN-ish empty Series."""
    m = pd.Series(True, index=df.index)
    for k, v in eq.items():
        m &= np.isclose(df[k], v) if isinstance(v, float) else (df[k] == v)
    sub = df[m]
    return sub.iloc[0] if len(sub) else pd.Series(dtype=float)


def collect(slug: str) -> list[tuple[str, str, str, str]]:
    """Return (section, label, value-string, source) rows for the report."""
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    # ---- form_compare: CV/cell, dCV, paired_z, BIC per form ------------
    fc = pd.read_csv(ROOT / "form_compare.csv")
    fc = fc[(fc.slug == slug) & np.isclose(fc.nu, NU)]
    src_fc = "form_compare.csv (slug, nu=8)"

    def fc_val(form: str, col: str) -> float:
        return float(_row(fc, form=form).get(col, np.nan))

    cv_chosen = fc_val(CHOSEN, "cv_per_cell")
    cv_base = fc_val(BASELINE, "cv_per_cell")
    z_base = fc_val(BASELINE, "paired_z")
    dcv_chosen = fc_val(CHOSEN, "dCV_vs_best")
    z_chosen = fc_val(CHOSEN, "paired_z")
    bic_chosen = fc_val(CHOSEN, "bic")
    bic_goff = fc_val(GOFF, "bic")
    cv_gscalar = fc_val(GSCALAR, "cv_per_cell")
    z_gscalar = fc_val(GSCALAR, "paired_z")

    # CV-best form for this slice (from the rollup).
    bf = pd.read_csv(ROOT / "best_form_all.csv")
    bf_row = _row(bf, slug=slug, nu=NU)
    cv_best_form = str(bf_row.get("best_cand", "?"))
    cv_best = float(bf_row.get("cv_per_cell", np.nan))
    bic_cvbest = fc_val(cv_best_form, "bic")

    A = "A. AGING BLOCK vs NONE  (does an aging term help at all?)"
    add(A, "held-out logdens/cell  (no aging)", f"{cv_base:.4f}",
        f"{src_fc}, form='{BASELINE}', cv_per_cell")
    add(A, f"held-out logdens/cell  ({CHOSEN})", f"{cv_chosen:.4f}",
        f"{src_fc}, form='{CHOSEN}', cv_per_cell")
    add(A, "  -> CV gain from adding aging", f"+{cv_chosen - cv_base:.4f} nats/cell",
        "derived")
    add(A, "  paired z (no-aging vs CV-best, same folds)", f"{z_base:.1f}",
        f"{src_fc}, form='{BASELINE}', paired_z")
    bic_base = fc_val(BASELINE, "bic")
    add(A, "  in-sample BIC  (no aging / chosen)", f"{bic_base:,.0f} / {bic_chosen:,.0f}",
        f"{src_fc}, bic")
    add(A, "  => term earns its keep?",
        f"in-sample BIC {bic_chosen - bic_base:,.0f} AND held-out CV +{cv_chosen - cv_base:.4f}  "
        "(both improve)", "derived")

    B = "B. WHICH FORM  (differences among flexible forms are tiny)"
    add(B, "CV-best form on this slice", f"{cv_best_form}  (CV {cv_best:.4f})",
        "best_form_all.csv (slug, nu=8), best_cand / cv_per_cell")
    add(B, f"chosen {CHOSEN} CV/cell", f"{cv_chosen:.4f}",
        f"{src_fc}, form='{CHOSEN}', cv_per_cell")
    add(B, "  -> CV shortfall vs CV-best", f"{dcv_chosen:.2e} nats/cell  (z {z_chosen:.1f})",
        f"{src_fc}, dCV_vs_best / paired_z")
    # effective d.o.f. k = AIC/2 + loglik (no N needed); naive n_params = 'k' col.
    aic_chosen = fc_val(CHOSEN, "aic")
    ll_chosen = fc_val(CHOSEN, "loglik")
    dof_chosen = aic_chosen / 2 + ll_chosen
    add(B, f"AIC / BIC  {CHOSEN}", f"{aic_chosen:,.0f} / {bic_chosen:,.0f}",
        f"{src_fc}, form='{CHOSEN}', aic / bic")
    add(B, "  effective d.o.f. / naive n_params", f"{dof_chosen:,.0f} / {fc_val(CHOSEN, 'k'):,.0f}  "
        "(equal: no EB shrinkage, drift off)", f"derived / {src_fc}, k")
    add(B, f"BIC  CV-best ({cv_best_form})", f"{bic_cvbest:,.0f}",
        f"{src_fc}, form='{cv_best_form}', bic")
    add(B, "  -> BIC gap chosen vs CV-best", f"{bic_chosen - bic_cvbest:,.0f}  (scale ~1.8e6)",
        "derived")

    # ---- curve_metrics: tail-wiggle guardrail + peak age --------------
    cm = pd.read_csv(ROOT / "curve_metrics_all.csv")
    cm = cm[(cm.slug == slug) & np.isclose(cm.nu, NU) & (cm.entry_age_label == "mean")]
    src_cm = "curve_metrics_all.csv (slug, nu=8, entry_age=mean)"

    def cm_val(cand: str, col: str) -> float:
        return float(_row(cm, cand=cand).get(col, np.nan))

    C = "C. CURVE PLAUSIBILITY GUARDRAIL  (spline vs high-order polynomial)"
    add(C, f"tail non-monotone?  {HIPOLY}", f"{cm_val(HIPOLY, 'tail_nonmonotone'):.0f}  "
        f"(sign changes {cm_val(HIPOLY, 'tail_sign_changes'):.0f})",
        f"{src_cm}, cand='{HIPOLY}', tail_nonmonotone / tail_sign_changes")
    add(C, f"tail non-monotone?  {CHOSEN}", f"{cm_val(CHOSEN, 'tail_nonmonotone'):.0f}  "
        f"(sign changes {cm_val(CHOSEN, 'tail_sign_changes'):.0f})",
        f"{src_cm}, cand='{CHOSEN}', tail_nonmonotone / tail_sign_changes")
    add(C, f"peak (improvement) age  {CHOSEN}", f"{cm_val(CHOSEN, 'peak_age'):.1f} yr",
        f"{src_cm}, cand='{CHOSEN}', peak_age")

    D = "D. ENTRY-AGE INTERACTION  (varying gamma, after Stones 2019)"
    add(D, f"BIC  no entry-age term ({GOFF})", f"{bic_goff:,.0f}",
        f"{src_fc}, form='{GOFF}', bic")
    add(D, f"BIC  varying gamma ({CHOSEN})", f"{bic_chosen:,.0f}",
        f"{src_fc}, form='{CHOSEN}', bic")
    add(D, "  -> BIC improvement from entry-age term", f"{bic_goff - bic_chosen:,.0f}",
        "derived")
    add(D, f"CV/cell  scalar gamma ({GSCALAR})", f"{cv_gscalar:.4f}  (z {z_gscalar:.1f})",
        f"{src_fc}, form='{GSCALAR}', cv_per_cell / paired_z")
    add(D, f"CV/cell  varying gamma ({CHOSEN})", f"{cv_chosen:.4f}",
        f"{src_fc}, form='{CHOSEN}', cv_per_cell")

    # ---- v_vs_reference: v robust to the aging form -------------------
    vr = pd.read_csv(ROOT / "v_vs_reference_all.csv")
    vr_row = _row(vr[(vr.slug == slug) & np.isclose(vr.nu, NU)],
                  cand=CHOSEN, reference=V_REF)
    src_vr = f"v_vs_reference_all.csv (slug, nu=8, cand='{CHOSEN}', ref='{V_REF}')"

    E = "E. v ROBUST TO FORM  (the target output barely moves)"
    add(E, f"v Pearson   {CHOSEN} vs {V_REF}", f"{float(vr_row.get('pearson', np.nan)):.4f}",
        f"{src_vr}, pearson")
    add(E, f"v Spearman  {CHOSEN} vs {V_REF}", f"{float(vr_row.get('spearman', np.nan)):.4f}",
        f"{src_vr}, spearman")
    add(E, "  mean |dv| vs reference (log-time)", f"{float(vr_row.get('mean_abs_dv', np.nan)):.4f}",
        f"{src_vr}, mean_abs_dv")

    # ---- spline4 vs neighbours: simpler vs spline5/6, richer vs spline3/poly --
    # Sources: in-sample (loglik/aic/bic) from grid/metrics.csv (ALL forms),
    # CV/cell from cv/form_selection.csv (CV-scored forms only -- poly4/6 were NOT
    # CV-scored, so their CV shows n/a and only their in-sample BIC is comparable).
    # "frac captured" = share of the WHOLE aging block's held-out value-add a form
    # buys = (cv_form - cv_noaging) / (cv_CVbest - cv_noaging).
    block_val = cv_best - cv_base  # CV-best minus no-aging
    mt = pd.read_csv(ROOT / slug / "grid" / "metrics.csv")
    mt = mt[np.isclose(mt.nu, NU)]
    mt_best = mt.loc[mt.groupby("cand")["loglik"].idxmax()].set_index("cand")
    fsel = pd.read_csv(ROOT / slug / "cv" / "form_selection.csv")
    fsel = fsel[np.isclose(fsel.nu, NU)].set_index("cand")

    def _ladder(form: str):
        ll = float(mt_best.loc[form, "loglik"]); aic = float(mt_best.loc[form, "aic"])
        bic = float(mt_best.loc[form, "bic"]); dof = aic / 2 + ll
        cv = float(fsel.loc[form, "cv_per_cell"]) if form in fsel.index else np.nan
        frac = (cv - cv_base) / block_val if (np.isfinite(cv) and block_val) else np.nan
        return dof, bic, cv, frac

    # spline3 (simpler) .. spline6 (CV-best), then high-order polys (not CV-scored)
    ladder = ["spline3-gvarying", "spline4-gvarying", "spline5-gvarying",
              "spline6-gvarying", "poly4-gvarying", "poly6-gvarying"]
    G = "F. WHY spline4 (vs spline3 / spline5,6 / poly4,6)"
    for form in ladder:
        dof_f, bic_f, cv_f, frac = _ladder(form)
        star = "  <- CHOSEN" if form == CHOSEN else (
            "  <- CV-best" if form == cv_best_form else "")
        cv_s = f"CV {cv_f:.6f}  capt {frac * 100:5.1f}%" if np.isfinite(cv_f) else "CV   n/a  (not CV-scored)"
        add(G, f"{form}", f"dof {dof_f:,.0f}  BIC {bic_f:,.0f}  {cv_s}{star}",
            f"grid/metrics.csv + cv/form_selection.csv")
    # direct "how much better is spline4" contrasts
    dof3, bic3, cv3, fr3 = _ladder("spline3-gvarying")
    add(G, "  spline4 vs spline3 (too simple)",
        f"CV +{cv_chosen - cv3:.4f}/cell ({fr3 * 100:.0f}%->{(cv_chosen - cv_base) / block_val * 100:.0f}% capt), "
        f"BIC {bic_chosen - bic3:,.0f}", "derived")
    add(G, "  spline4 vs spline6 (CV-best)",
        f"CV {dcv_chosen:.2e}/cell (last {(1 - (cv_chosen - cv_base) / block_val) * 100:.1f}%), "
        f"BIC +{bic_chosen - bic_cvbest:,.0f}", "derived")
    _, bic_p4, _, _ = _ladder("poly4-gvarying")
    _, bic_p6, _, _ = _ladder("poly6-gvarying")
    add(G, "  spline4 vs poly4 / poly6 (in-sample BIC only)",
        f"dBIC {bic_chosen - bic_p4:+,.0f} / {bic_chosen - bic_p6:+,.0f}  "
        "(<0 = spline4 better; poly route rejected on tail shape, see C)", "derived")
    # v: simplifying spline6 -> spline4 barely moves the difficulty index.
    vx = pd.read_csv(ROOT / slug / "v_xform_agreement.csv")
    vx = vx[np.isclose(vx.nu, NU)]
    vmask = (((vx.form_a == CHOSEN) & (vx.form_b == "spline6-gvarying"))
             | ((vx.form_a == "spline6-gvarying") & (vx.form_b == CHOSEN)))
    vrow = vx[vmask]
    if len(vrow):
        add(G, "  v agreement spline4-gv vs spline6-gv",
            f"Pearson {float(vrow.pearson.iloc[0]):.4f}, Spearman {float(vrow.spearman.iloc[0]):.4f}",
            f"{slug}/v_xform_agreement.csv")
    # CV-best knot count is unstable across the ALL slices (none is spline4).
    bf_all = pd.read_csv(ROOT / "best_form_all.csv")
    bf8 = bf_all[np.isclose(bf_all.nu, NU) & bf_all.slug.str.startswith("ALL")]
    winners = ", ".join(f"{r.slice}:{r.best_cand}" for _, r in bf8.iterrows())
    add(G, "  CV-best form across ALL slices", winners,
        "best_form_all.csv (nu=8, ALL_*) -- unstable knot count, all within ~5e-4 of spline4")

    return rows


def render(slug: str, rows: list[tuple[str, str, str, str]]) -> str:
    return render_markdown(
        f"Aging progression form selection ({slug})",
        rows,
        subtitle=[
            f"Production form = {CHOSEN} (natural cubic spline, 4 knots, varying "
            "gamma); nu=8; drift d_i OFF during selection.",
            "Decision: held-out CV/cell primary -> curve plausibility guardrail "
            "-> in-sample BIC tie-break.",
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

    out = write_markdown(ROOT / f"aging_summary_{slug}.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
