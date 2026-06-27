"""One-stop number lookup for the section-4 covariate paragraph (QC, no fit).

Reads the CSVs the covariate scripts (q02-q08) already wrote under their per-script
subdirectories and pulls the exact numbers the paper's section 4 quotes, printing
them on screen and to a human-readable text file. Each number is annotated with
its source file + column so it is trivial to verify and cite. Mirrors
``02_model_selection/baseline/q99_baseline_summary.py``.

Sources (under results/analysis/covariate/):
  * 02_variable_selection/variable_selection__ALL_B.csv : best weather/course var.
  * 03_regression/regression__6slices.csv               : v_j ~ temp+gain R2 per
                                                          slice + simple-metric R2.
  * 08_model_comparison/model_comparison__ALL_B.csv      : R2 per model (baseline/
                                                          aging/drift/full).
  * 04_regression_robustness/regression_robustness__6slices.csv : 2x2 OLS/WLS/
                                                          Robust fits.
  * 05_cluster_robust/cluster_robust__ALL_B.csv          : cluster-robust z.
  * 06_mixed_model/mixed_model__ALL_B.csv                : random-course mixed model.

Output -> results/analysis/covariate/covariate_summary_{slice}.md
(human-readable Markdown the paper cites directly).

The single-slice analyses (q02/q05/q06/q08) are pinned to ALL_B; ``--slice`` only
selects which slice the q03/q04 blocks highlight (default ALL_B).

Run::

    python scripts/05_analysis/covariate/q99_covariate_summary.py
    python scripts/05_analysis/covariate/q99_covariate_summary.py --slice Po10_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/ (for report_md)

import covariate_common as C
from report_md import render_markdown, write_markdown  # noqa: E402

COV = C.OUT_ROOT
VAR_SLICE = C.VAR_SELECT_SLICE          # ALL_B; the pinned single-slice analyses

FILES = {
    "varsel":  COV / "02_variable_selection" / f"variable_selection__{VAR_SLICE}.csv",
    "regr":    COV / "03_regression" / f"regression__{C.ALL_SLICES_TAG}.csv",
    "robust":  COV / "04_regression_robustness"
               / f"regression_robustness__{C.ALL_SLICES_TAG}.csv",
    "cluster": COV / "05_cluster_robust" / f"cluster_robust__{VAR_SLICE}.csv",
    "mixed":   COV / "06_mixed_model" / f"mixed_model__{VAR_SLICE}.csv",
    "models":  COV / "08_model_comparison" / f"model_comparison__{VAR_SLICE}.csv",
}


def _src(p: Path) -> str:
    """File path relative to the covariate output root, for source annotations."""
    return str(p.relative_to(COV)).replace("\\", "/")


def collect(slice_name: str) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    def note_missing(section, key):
        add(section, "(missing)", f"{FILES[key].name} not found -- run the script",
            _src(FILES[key]))

    # == A. variable selection: EVERY candidate covariate (q02, ALL_B) =
    A = f"VARIABLE SELECTION  v_j vs each covariate  (q02, {VAR_SLICE})"
    if FILES["varsel"].is_file():
        vs = pd.read_csv(FILES["varsel"])
        src = _src(FILES["varsel"])
        for family in ("weather", "course"):
            sub = vs[vs.family == family].copy()
            sub = sub.reindex(sub.spearman.abs().sort_values(ascending=False).index)
            for i, (_, r) in enumerate(sub.iterrows()):
                mark = "  <= best" if i == 0 else ""
                add(A, f"{family:7s} {r.label}",
                    f"Pearson {r.pearson:+.3f}  Spearman {r.spearman:+.3f} "
                    f"[{r.spearman_lo:+.2f},{r.spearman_hi:+.2f}]  (n={int(r.n)}){mark}",
                    f"{src}, row {r.covariate}")
    else:
        note_missing(A, "varsel")

    # == B. joint regression v_j ~ weather+gain (q03) =================
    if FILES["regr"].is_file():
        rg = pd.read_csv(FILES["regr"])
        srcr = _src(FILES["regr"])
        sl = rg[rg.kind == "slice"]
        for wset, wlab in (("temp", "temp"), ("WBGT", "WBGT")):
            B = f"JOINT REGRESSION  v_j ~ {wlab}+gain  (q03)"
            sub = sl[sl.weather_set == wset].set_index("target")
            for s in C.SLICE_ORDER:
                key = f"v_{s}"
                if key not in sub.index:
                    continue
                r = sub.loc[key]
                mark = "  <= highlight" if s == slice_name else ""
                add(B, f"{s} R2 [95% CI]",
                    f"{r.r2:.3f} [{r.r2_lo:.3f},{r.r2_hi:.3f}]  "
                    f"b_{wlab} {r.beta_weather:+.3f}(z{r.z_weather:+.1f})  "
                    f"b_gain {r.beta_course:+.3f}(z{r.z_course:+.1f}){mark}",
                    f"{srcr}, weather_set={wset}, target={key}")

        # v_j vs simple metrics: same joint ~ temp+gain, directly comparable
        # (R2 with CI + BOTH partial betas, so temp and course line up with v_j).
        M = "v_j vs SIMPLE METRICS  (joint ~ temp+gain; directly comparable)"

        def reg_line(label, r, source):
            add(M, label,
                f"R2 {r.r2:.3f} [{r.r2_lo:.3f},{r.r2_hi:.3f}]  "
                f"b_temp {r.beta_weather:+.3f}(z{r.z_weather:+.1f})  "
                f"b_gain {r.beta_course:+.3f}(z{r.z_course:+.1f})", source)

        ref = sl[sl.weather_set == "temp"].set_index("target").loc[f"v_{VAR_SLICE}"]
        reg_line(f"v_{VAR_SLICE} (full model)  [REF]", ref,
                 f"{srcr}, target=v_{VAR_SLICE}")
        met = (rg[(rg.kind == "metric") & (rg.weather_set == "temp")]
               .sort_values("r2", ascending=False))
        for _, r in met.head(10).iterrows():
            reg_line(r.target, r, f"{srcr}, kind=metric, target={r.target}")

        vr2, best, ms = float(ref.r2), met.iloc[0], met.set_index("target")
        add(M, f"-> R2 ratio  v_{VAR_SLICE} / best metric",
            f"{vr2:.3f} / {best.r2:.3f} = {vr2 / best.r2:.1f}x  ({best.target})", "derived")
        for mname in ("mean_time_sec", "median_time_sec"):
            if mname in ms.index:
                add(M, f"-> R2 ratio  v_{VAR_SLICE} / {mname}",
                    f"{vr2:.3f} / {ms.loc[mname, 'r2']:.3f} = "
                    f"{vr2 / ms.loc[mname, 'r2']:.1f}x", "derived")
    else:
        note_missing("JOINT REGRESSION  (q03)", "regr")

    # == C. model comparison (q08, ALL_B) =============================
    MC = f"MODEL COMPARISON  v_j ~ temp+gain  (q08, {VAR_SLICE})"
    if FILES["models"].is_file():
        mc = pd.read_csv(FILES["models"]).set_index("model")
        srcm = _src(FILES["models"])
        for m in C.MODELS:
            if m not in mc.index:
                continue
            r = mc.loc[m]
            add(MC, f"{m:9s} R2 [95% CI]",
                f"{r.r2:.3f} [{r.r2_lo:.3f},{r.r2_hi:.3f}]  "
                f"rho(temp) {r.spearman_temp:.3f}", f"{srcm}, model={m}")
        if {"baseline", "aging", "drift", "full"}.issubset(mc.index):
            r2 = mc["r2"]
            add(MC, "-> aging gain (baseline->aging)",
                f"{r2['aging'] - r2['baseline']:+.3f}", "derived")
            add(MC, "-> d_i gain  (baseline->drift)",
                f"{r2['drift'] - r2['baseline']:+.3f}", "derived")
            add(MC, "-> d_i gain  (aging->full)",
                f"{r2['full'] - r2['aging']:+.3f}", "derived")
            add(MC, "-> best model", f"{r2.idxmax()} (R2 {r2.max():.3f})", "derived")
    else:
        note_missing(MC, "models")

    # == D. robustness: 2x2 fits (q04) ================================
    RB = f"ROBUSTNESS  v_j ~ temp+gain, 2x2 fits  (q04, {slice_name})"
    if FILES["robust"].is_file():
        rb = pd.read_csv(FILES["robust"])
        srcb = _src(FILES["robust"])
        a = rb[(rb.kind == "slice") & (rb.target == f"v_{slice_name}")].set_index("fit")
        for fk in ("OLS", "WLS", "Robust", "Robust+WLS"):
            if fk in a.index:
                r = a.loc[fk]
                add(RB, fk, f"b_temp {r.b_temp:+.3f}(z{r.z_temp:+.1f})  "
                    f"b_gain {r.b_gain:+.3f}(z{r.z_gain:+.1f})  R2 {r.r2:.3f}",
                    f"{srcb}, target=v_{slice_name}, fit={fk}")
    else:
        note_missing(RB, "robust")

    # == E. over-counting: cluster-robust SE (q05) ====================
    CL = f"OVER-COUNTING #1: cluster-robust SE  (q05, {VAR_SLICE})"
    if FILES["cluster"].is_file():
        cl = pd.read_csv(FILES["cluster"])
        srcc = _src(FILES["cluster"])

        def clt(fit, term):
            s = cl[(cl.fit == fit) & (cl.term == term)]
            return float(s.t.iloc[0]) if len(s) else float("nan")

        # the clustered fits (CR1/CR2/pairs_boot) carry dof = G-1; iid fits carry
        # n-k, so read G off a clustered row, not the max.
        cr = cl[cl.fit.isin(["CR1", "CR2"])]
        G = int(cr["dof"].iloc[0]) + 1 if len(cr) else -1
        n_ed = int(cl["dof"].max()) + 3             # iid dof = n - k (k=3)
        add(CL, "physical courses G (clusters)", f"{G}  (vs n={n_ed} editions)", srcc)
        add(CL, "temp z: iid -> CR2",
            f"{clt('iid', 'temp'):.1f} -> {clt('CR2', 'temp'):.1f}  (robust)", srcc)
        add(CL, "gain z: iid -> CR1 -> CR2",
            f"{clt('iid', 'gain'):.1f} -> {clt('CR1', 'gain'):.1f} -> "
            f"{clt('CR2', 'gain'):.1f}", srcc)
        add(CL, "gain z (WLS): iid -> CR2",
            f"{clt('WLS_iid', 'gain'):.1f} -> {clt('WLS_CR2', 'gain'):.1f}", srcc)
    else:
        note_missing(CL, "cluster")

    # == F. over-counting: random-course mixed model (q06) ============
    MX = f"OVER-COUNTING #2: random-course mixed model  (q06, {VAR_SLICE})"
    if FILES["mixed"].is_file():
        mx = pd.read_csv(FILES["mixed"]).set_index("term")
        srcx = _src(FILES["mixed"])
        for term in ("temp", "gain"):
            r = mx.loc[term]
            z_ols = r.b_ols / r.se_ols
            add(MX, f"{term}: OLS z -> mixed z",
                f"{z_ols:+.1f} -> {r.z:+.1f}  (SE x{r.se_ratio:.2f})",
                f"{srcx}, term={term}")
        if "_var" in mx.index:
            add(MX, "ICC (between-course share of resid var)",
                f"{mx.loc['_var', 'z']:.3f}", f"{srcx}, _var row (z=ICC)")
    else:
        note_missing(MX, "mixed")

    return rows


def render(slice_name: str, rows: list[tuple[str, str, str, str]]) -> str:
    return render_markdown(
        f"Covariate validation: v_j vs weather + course ({slice_name})",
        rows,
        subtitle=[
            "Does the race factor v_j track objective race-day conditions? Post-hoc "
            "on the fitted v_j (full production model, beta=0 gauge); no refit.",
            "Single-slice analyses (q02/q05/q06/q08) are pinned to ALL_B; the q03/q04 "
            "blocks highlight the chosen slice.",
        ],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default="ALL_B",
                    help="slice the q03/q04 blocks highlight (default ALL_B)")
    args = ap.parse_args()

    rows = collect(args.slice)
    report = render(args.slice, rows)
    print(report)

    out = write_markdown(COV / f"covariate_summary_{args.slice}.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
