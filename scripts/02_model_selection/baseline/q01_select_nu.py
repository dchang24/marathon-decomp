"""Quantitative single-nu decision for the rank-1 baseline (QC / diagnostic).

Turns the per-slice CV sweep into a *defensible* single shared nu, using two
criteria that need no new fitting (everything is read from the CSVs that
``e01_nu_cv.py`` already wrote):

  1. 1-SE rule on CV.  At each nu, CV(nu) = mean over the K folds of the
     held-out log predictive density PER CELL (``heldout_per_cell`` in
     ``cv_folds.csv``); SE(nu) = sd_fold / sqrt(K).  The acceptable set for a
     slice is {nu : CV(nu) >= max CV - SE_at_max}.  This widens each flat CV top
     into an interval [nu_lo, nu_hi] and guards against chasing fold noise.
     (Folds are stratified within-athlete, so the across-fold SE is a principled
     heuristic, not a strict CI — it slightly *over*states the spread.)

  2. v-stability plateau.  For each slice, 1 - corr(v(nu), v(nu_argmax)) on the
     full-data fits (``v_xnu.csv``) measures how much the *deliverable* race
     factor moves with nu.  If the deliverable is invariant across the candidate
     range, the exact nu is immaterial and a single shared value is justified on
     comparability grounds even when the per-slice 1-SE intervals don't all
     overlap.

Decision: a single shared nu = the grid value covered by the most slices' 1-SE
intervals (ties -> closest to the center of the per-slice argmaxes).  The script
reports, at that nu, the per-slice v-cost (max 1-corr to each slice's own
argmax-v) so the comparability trade is explicit.

Outputs -> results/model_selection/baseline/ :
  * nu_decision.csv   : per-slice argmax/Brent nu, 1-SE interval, v-plateau width,
                        v-cost at the recommended nu; plus the cross-slice summary.
  * nu_decision.png   : per-slice CV +- SE vs nu (cohort x gender grid), 1-SE
                        threshold line + acceptable band shaded, recommended nu marked.

Self-contained; no arguments needed (VS Code "Run" works).

Run::

    python scripts/02_model_selection/baseline/q01_select_nu.py
    python scripts/02_model_selection/baseline/q01_select_nu.py --solver cv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from report_md import render_markdown, write_markdown  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"
COHORT_ORDER = ["ALL", "Po10", "WA"]
GENDER_ORDER = ["M", "W", "B"]
_GENDER_LABEL = {"M": "men", "W": "women", "B": "both"}


METRICS = ("pearson_1mcorr", "spearman_1mcorr", "max_abs_dv", "mean_abs_dv")


def _one_minus_corr(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) < 1e-15 or np.std(b) < 1e-15:
        return float("nan")
    return float(1.0 - np.corrcoef(a, b)[0, 1])


def _v_metrics(a: np.ndarray, b: np.ndarray) -> dict:
    """v-cost between two v vectors, four ways. Correlations are shift-invariant;
    the |Δv| metrics center each vector first (the model's mean(v)=0 gauge)."""
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    d = (a - a.mean()) - (b - b.mean())
    return dict(pearson_1mcorr=_one_minus_corr(a, b),
                spearman_1mcorr=_one_minus_corr(ra, rb),
                max_abs_dv=float(np.max(np.abs(d))),
                mean_abs_dv=float(np.mean(np.abs(d))))


def analyse_slice(slug: str, solver: str) -> dict | None:
    """Per-slice 1-SE interval + v-plateau; returns a dict or None if missing."""
    sdir = OUT_ROOT / slug
    folds_p, v_p = sdir / "cv_folds.csv", sdir / "v_xnu.csv"
    if not folds_p.is_file() or not v_p.is_file():
        return None

    folds = pd.read_csv(folds_p)
    folds = folds[folds.solver == solver]
    if folds.empty:                       # fall back to any present solver
        folds = pd.read_csv(folds_p)
        solver = sorted(folds.solver.unique())[0]
        folds = folds[folds.solver == solver]

    agg = (folds.groupby("nu")["heldout_per_cell"]
           .agg(mean="mean", sd="std", k="count").reset_index())
    agg["se"] = agg["sd"] / np.sqrt(agg["k"])
    fin = agg[np.isfinite(agg["nu"])].sort_values("nu").reset_index(drop=True)

    i_best = int(fin["mean"].idxmax())
    nu_best = float(fin.loc[i_best, "nu"])
    thr = float(fin.loc[i_best, "mean"] - fin.loc[i_best, "se"])
    ok = fin[fin["mean"] >= thr]
    nu_lo, nu_hi = float(ok["nu"].min()), float(ok["nu"].max())

    # full-data v per nu (same solver); ref = v at the CV-argmax nu.
    v = pd.read_csv(v_p)
    v = v[v.solver == solver] if (v.solver == solver).any() else v
    piv = v.pivot(index="race_id", columns="nu", values="v").dropna()
    ref = piv[nu_best].to_numpy()
    vstab = {float(nu): _one_minus_corr(piv[nu].to_numpy(), ref) for nu in piv.columns}
    # plateau width: max 1-corr over the 1-SE-acceptable nus (how much v moves
    # across the whole CV-acceptable band).
    plateau = max((vstab[nu] for nu in vstab if nu_lo <= nu <= nu_hi and np.isfinite(nu)),
                  default=float("nan"))

    return dict(slug=slug, solver=solver, nu_best=nu_best, cv_best=float(fin["mean"].max()),
                se_best=float(fin.loc[i_best, "se"]), nu_lo=nu_lo, nu_hi=nu_hi,
                plateau_1mcorr=plateau, fin=fin, vstab=vstab, piv=piv)


def recommend(rows: list[dict]) -> tuple[float, dict]:
    """Single shared nu = grid value in the most 1-SE intervals.

    Ties -> closest to the mean of per-slice argmaxes. Also reports the
    intersection of all intervals (possibly empty) and the per-slice v-cost at
    the chosen nu (1-corr to that slice's own argmax-v).
    """
    grid = sorted({float(nu) for r in rows for nu in r["fin"]["nu"]})
    coverage = {g: sum(r["nu_lo"] <= g <= r["nu_hi"] for r in rows) for g in grid}
    center = float(np.mean([r["nu_best"] for r in rows]))
    best_cov = max(coverage.values())
    cands = [g for g, c in coverage.items() if c == best_cov]
    nu_star = min(cands, key=lambda g: abs(g - center))

    inter_lo = max(r["nu_lo"] for r in rows)
    inter_hi = min(r["nu_hi"] for r in rows)
    intersection = (inter_lo, inter_hi) if inter_lo <= inter_hi else None

    nan4 = {m: float("nan") for m in METRICS}
    vcost = {}
    for r in rows:
        piv, nb = r["piv"], r["nu_best"]
        if nu_star in piv.columns and nb in piv.columns:
            vcost[r["slug"]] = _v_metrics(piv[nu_star].to_numpy(), piv[nb].to_numpy())
        else:
            vcost[r["slug"]] = dict(nan4)
    max_vcost = {m: float(np.nanmax([vcost[s][m] for s in vcost])) for m in METRICS}
    return nu_star, dict(coverage=coverage, n_slices=len(rows), center=center,
                         intersection=intersection, vcost=vcost, max_vcost=max_vcost)


def _split(slug: str) -> tuple[str, str]:
    parts = slug.split("_")
    return parts[0], parts[1]


def plot(rows: list[dict], nu_star: float, info: dict, out: Path) -> None:
    keyed = {_split(r["slug"]): r for r in rows}
    cohorts = [c for c in COHORT_ORDER if any(k[0] == c for k in keyed)]
    genders = [g for g in GENDER_ORDER if any(k[1] == g for k in keyed)]
    nrow, ncol = len(cohorts), len(genders)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             squeeze=False, constrained_layout=True)
    for r_i, cohort in enumerate(cohorts):
        for c_i, gender in enumerate(genders):
            ax = axes[r_i, c_i]
            r = keyed.get((cohort, gender))
            if r is None:
                ax.set_axis_off()
                continue
            fin = r["fin"]
            thr = r["cv_best"] - r["se_best"]
            ax.errorbar(fin["nu"], fin["mean"], yerr=fin["se"], fmt="-o", ms=4,
                        color="#4363d8", capsize=2, lw=1.2)
            ax.axhline(thr, ls="--", color="#888", lw=0.9)
            ax.axvspan(r["nu_lo"], r["nu_hi"], color="#3cb44b", alpha=0.12)
            ax.axvline(nu_star, ls="-", color="#e6194b", lw=1.4)
            ax.set_xscale("log")
            ax.grid(True, which="both", alpha=0.25, lw=0.4)
            ax.set_title(f"{cohort} — {_GENDER_LABEL.get(gender, gender)}   "
                         f"[1-SE: {r['nu_lo']:g}–{r['nu_hi']:g}]", fontsize=9)
            if c_i == 0:
                ax.set_ylabel("CV logdens / cell  (±SE)")
            if r_i == nrow - 1:
                ax.set_xlabel("nu  (log)")
    inter = info["intersection"]
    inter_s = f"{inter[0]:g}–{inter[1]:g}" if inter else "empty"
    mv = info["max_vcost"]
    fig.suptitle(f"Single-nu decision  —  recommended nu*={nu_star:g}  "
                 f"(1-SE intersection: {inter_s}; max v-cost: 1−Pear={mv['pearson_1mcorr']:.1e}, "
                 f"max|Δv|={mv['max_abs_dv']:.1e})\n"
                 f"green = per-slice 1-SE band, dashed = max−SE, red = nu*", fontsize=11)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver", default="anderson",
                    help="CV solver to read (default anderson).")
    args = ap.parse_args()

    slugs = sorted(p.parent.name for p in OUT_ROOT.glob("*/cv_folds.csv"))
    rows = [r for s in slugs if (r := analyse_slice(s, args.solver)) is not None]
    if not rows:
        print(f"No cv_folds.csv under {OUT_ROOT} — run e01_nu_cv.py first.")
        return

    nu_star, info = recommend(rows)

    # ── report table (v-cost at nu* = v(nu*) vs v(this slice's argmax), 4 ways) ─
    hdr = (f"{'slice':<22}{'argmax':>7}{'1SElo':>6}{'1SEhi':>6}{'CV*':>8}{'SE':>8}"
           f"{'1-Pear':>9}{'1-Spear':>9}{'max|dv|':>9}{'mean|dv|':>9}")
    print(hdr); print("-" * len(hdr))
    out_rows = []
    for r in sorted(rows, key=lambda x: x["slug"]):
        vc = info["vcost"][r["slug"]]
        print(f"{r['slug']:<22}{r['nu_best']:>7g}{r['nu_lo']:>6g}{r['nu_hi']:>6g}"
              f"{r['cv_best']:>8.4f}{r['se_best']:>8.4f}"
              f"{vc['pearson_1mcorr']:>9.1e}{vc['spearman_1mcorr']:>9.1e}"
              f"{vc['max_abs_dv']:>9.1e}{vc['mean_abs_dv']:>9.1e}")
        out_rows.append(dict(slice=r["slug"], solver=r["solver"], nu_argmax=r["nu_best"],
                             nu_1se_lo=r["nu_lo"], nu_1se_hi=r["nu_hi"],
                             cv_best_per_cell=r["cv_best"], se_at_best=r["se_best"],
                             plateau_1mcorr_over_band=r["plateau_1mcorr"],
                             vcost_pearson_1mcorr=vc["pearson_1mcorr"],
                             vcost_spearman_1mcorr=vc["spearman_1mcorr"],
                             vcost_max_abs_dv=vc["max_abs_dv"],
                             vcost_mean_abs_dv=vc["mean_abs_dv"]))

    inter = info["intersection"]
    print("\n" + "=" * 60)
    print(f"per-slice argmax nu : {[r['nu_best'] for r in rows]}")
    print(f"1-SE interval intersection : "
          f"{f'[{inter[0]:g}, {inter[1]:g}]' if inter else 'EMPTY'}")
    print(f"coverage by grid nu (n slices whose 1-SE band contains it):")
    for g, c in sorted(info["coverage"].items()):
        print(f"    nu={g:<5g} : {c}/{info['n_slices']}")
    print(f"\n>>> recommended single shared nu* = {nu_star:g}")
    print("    max deliverable v-cost at nu* across slices (vs each slice's own argmax-v):")
    mv = info["max_vcost"]
    print(f"        1-corr pearson  : {mv['pearson_1mcorr']:.2e}")
    print(f"        1-corr spearman : {mv['spearman_1mcorr']:.2e}")
    print(f"        max |dv|        : {mv['max_abs_dv']:.2e}  (log-time units)")
    print(f"        mean |dv|       : {mv['mean_abs_dv']:.2e}  (log-time units)")
    if inter is None:
        print("    NOTE: 1-SE intervals do not all overlap - the single-nu choice "
              "rests on the v-plateau:\n          CV can resolve nearby nu, but the "
              "v deliverable barely moves (see max v-cost above).")
    print("=" * 60)

    # ── write artifacts ─────────────────────────────────────────────
    summary = dict(slice="__RECOMMENDED__", solver=args.solver, nu_argmax=nu_star,
                   nu_1se_lo=(inter[0] if inter else np.nan),
                   nu_1se_hi=(inter[1] if inter else np.nan),
                   cv_best_per_cell=np.nan, se_at_best=np.nan,
                   plateau_1mcorr_over_band=np.nan,
                   vcost_pearson_1mcorr=mv["pearson_1mcorr"],
                   vcost_spearman_1mcorr=mv["spearman_1mcorr"],
                   vcost_max_abs_dv=mv["max_abs_dv"],
                   vcost_mean_abs_dv=mv["mean_abs_dv"])
    df = pd.DataFrame(out_rows + [summary])
    csv_p = OUT_ROOT / "nu_decision.csv"
    df.to_csv(csv_p, index=False)
    png_p = OUT_ROOT / "nu_decision.png"
    plot(rows, nu_star, info, png_p)

    md_p = write_markdown(OUT_ROOT / "nu_decision.md",
                          _markdown(rows, nu_star, info))
    print(f"\nwrote {csv_p}\nwrote {png_p}\nwrote {md_p}")


def _markdown(rows: list[dict], nu_star: float, info: dict) -> str:
    """Render the single-nu decision (per-slice 1-SE bands + recommendation)."""
    md_rows: list[tuple[str, str, str, str]] = []

    A = "Per-slice nu (1-SE rule on held-out CV/cell)"
    for r in sorted(rows, key=lambda x: x["slug"]):
        md_rows.append((
            A, r["slug"],
            f"argmax {r['nu_best']:g}; 1-SE band [{r['nu_lo']:g}, {r['nu_hi']:g}]; "
            f"CV* {r['cv_best']:.4f} +- {r['se_best']:.4f}",
            f"{r['slug']}/cv_folds.csv (heldout_per_cell, K folds)"))

    inter = info["intersection"]
    mv = info["max_vcost"]
    best_cov = max(info["coverage"].values())
    B = "Recommended single shared nu"
    md_rows.append((B, "recommended nu*", f"{nu_star:g}",
                    "derived: grid value in the most 1-SE bands"))
    md_rows.append((B, "1-SE interval intersection (all slices)",
                    f"[{inter[0]:g}, {inter[1]:g}]" if inter else "empty",
                    "derived from per-slice 1-SE bands"))
    md_rows.append((B, "coverage at nu*",
                    f"{best_cov}/{info['n_slices']} slices' 1-SE bands contain nu*",
                    "derived"))
    md_rows.append((B, "max deliverable v-cost at nu* (vs each slice's argmax-v)",
                    f"1-Pearson {mv['pearson_1mcorr']:.2e}, 1-Spearman "
                    f"{mv['spearman_1mcorr']:.2e}, max|dv| {mv['max_abs_dv']:.2e} "
                    f"(log-time)", "derived from v_xnu.csv"))

    C = "Coverage by grid nu (# slices whose 1-SE band contains it)"
    for g, c in sorted(info["coverage"].items()):
        md_rows.append((C, f"nu = {g:g}", f"{c}/{info['n_slices']}", "derived"))

    notes = [
        "The single shared nu rests on the v-plateau when the per-slice 1-SE "
        "intervals do not all overlap: CV can resolve nearby nu, but the "
        "deliverable race factor v barely moves (see max v-cost above).",
        "Folds are stratified within-athlete, so the across-fold SE slightly "
        "*over*states the spread (a principled heuristic, not a strict CI).",
    ]
    return render_markdown(
        "Single-nu decision for the rank-1 baseline",
        md_rows,
        subtitle="held-out CV/cell (1-SE rule) + v-stability plateau across data "
                 "subsets; no new fitting (reads the e01_nu_cv CV sweep).",
        notes=notes,
    )


if __name__ == "__main__":
    main()
