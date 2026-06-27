"""nu-sweep + K-fold CV selection for the rank-1 baseline (Student-t d.o.f.).

Task 2 (selection). In-sample profile loglik is monotone in heavy tails (it pins
nu to the grid floor), so nu is chosen by held-out predictive density via
stratified-within-athlete K-fold CV. For each nu we fit on each fold's training
set (warm continuation down the nu grid) and score the held-out finishes;
CV(nu) = total held-out log predictive density. A Brent refinement then maximizes
CV over a continuous bracket around the finite-grid argmax.

This script is SELF-CONTAINED: it fits its own L2 anchor (the nu=inf grid point
and warm-start source) so it has no dependency on the essential L2 deliverable.
Its full-data fits are SELECTION-grade (single warm continuation, tol ~1e-9); the
absolute-best tight-tol deliverable models are produced separately by
``scripts/03_model_fit/baseline``. The only file this hands downstream is
``selected_nu.csv``, which ``03/e02_fit_baseline_t.py`` reads.

Outputs (CSV) -> results/model_selection/baseline/{slug}/ :
  * nu_selection.csv    : per nu (grid + Brent) in-sample loglik/AIC/BIC + CV
                          logdens (total + per-cell) + held-out RMSE + is_selected.
  * cv_folds.csv        : per (nu, fold) held-out breakdown.
  * param_sensitivity.csv: per nu, how the full-data fit moves vs L2 (sigma2,
                          v-vs-L2 pearson/spearman, hardest-10% Jaccard, rank
                          shift, mean/max |dv|) — the nu-sensitivity table.
  * v_xnu.csv           : full-data fitted v per (nu, race_id), long form.
  * selected_nu.csv     : the one-row selection summary (read by 03/e02).
Cross-slice rollups (nu_selection_all.csv, selected_nu_all.csv) at the dir root.

Run::

    python scripts/02_model_selection/baseline/e01_nu_cv.py                 # the six slices
    python scripts/02_model_selection/baseline/e01_nu_cv.py --slices Po10_M
    python scripts/02_model_selection/baseline/e01_nu_cv.py --slices Po10_M --cv-solvers als anderson
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import (  # noqa: E402
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    assign_folds,
    heldout_logdensity,
    load_slice,
    subset_fitdata,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from baseline_common import slices as S  # noqa: E402
from baseline_common.fitting import baseline_cfg  # noqa: E402

# ── tunables (CLI flags override) ────────────────────────────────────
NU_GRID = (15.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0)   # plus nu=inf (L2), descended warm
K = 5
FOLD_SEED = 0
TOL = 1e-9            # selection-grade (the tight-tol deliverable lives in 03)
MAX_ITER = 2000
MIN_TEST_N = 1        # >1 restricts held-out cells to athletes with >= this many races
BRENT_XATOL = 0.01    # tolerance on log(nu) for Brent refinement
# ─────────────────────────────────────────────────────────────────────

INF = float("inf")
OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _fitter(solver: str, *, init="mean", warmstart=None, tol=TOL, max_iter=MAX_ITER):
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init=init, warmstart=warmstart, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def _uv(m) -> dict[str, np.ndarray]:
    return {"u": m.params["u"].copy(), "v": m.params["v"].copy()}


def _fit_at_nu(fd, nu, warm_uv, solver, *, tol, max_iter):
    if warm_uv is None:
        m = _model_cls(solver)(fd, baseline_cfg(nu), _fitter(solver, init="mean",
                               tol=tol, max_iter=max_iter))
    else:
        m = _model_cls(solver)(fd, baseline_cfg(nu),
                               _fitter(solver, init="warmstart", warmstart=warm_uv,
                                       tol=tol, max_iter=max_iter))
    m.fit()
    return m


def _continuation_chain(fd, solver, nu_grid, *, tol, max_iter) -> dict[float, object]:
    """Fit nu=inf (mean init), then descend the grid warm-starting each nu from
    the previous nu's solution. Returns {nu: fitted model}."""
    out: dict[float, object] = {}
    m = _fit_at_nu(fd, INF, None, solver, tol=tol, max_iter=max_iter)
    out[INF] = m
    prev = _uv(m)
    for nu in nu_grid:
        m = _fit_at_nu(fd, nu, prev, solver, tol=tol, max_iter=max_iter)
        out[nu] = m
        prev = _uv(m)
    return out


def _nearest_grid_uv(nu: float, chain_uv: dict[float, dict]) -> dict:
    finite = [g for g in chain_uv if np.isfinite(g)]
    g = min(finite, key=lambda x: abs(x - nu))
    return chain_uv[g]


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    if np.std(ra) < 1e-12 or np.std(rb) < 1e-12:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _v_sensitivity(v: np.ndarray, v_ref: np.ndarray) -> dict:
    """How a fitted v differs from the L2 reference v (same race ordering)."""
    J = len(v)
    k = max(1, int(round(0.10 * J)))
    top = set(np.argsort(v)[-k:])            # hardest 10% (largest v)
    top_ref = set(np.argsort(v_ref)[-k:])
    jacc = len(top & top_ref) / len(top | top_ref) if (top | top_ref) else float("nan")
    rank = pd.Series(v).rank().to_numpy()
    rank_ref = pd.Series(v_ref).rank().to_numpy()
    return dict(
        pearson_vs_L2=_pearson(v, v_ref), spearman_vs_L2=_spearman(v, v_ref),
        top10_jaccard_vs_L2=jacc,
        max_abs_rank_shift_frac=float(np.max(np.abs(rank - rank_ref)) / J),
        mean_abs_dv=float(np.mean(np.abs(v - v_ref))),
        max_abs_dv=float(np.max(np.abs(v - v_ref))),
    )


def run_solver(solver, fd, fold, K, nu_grid, *, tol, max_iter):
    nus = [INF, *nu_grid]
    print(f"\n#### solver = {solver} ####", flush=True)

    # ── full-data continuation chain (selection-grade) ───────────────
    print("-- full-data continuation --", flush=True)
    full = _continuation_chain(fd, solver, nu_grid, tol=tol, max_iter=max_iter)
    full_uv = {nu: _uv(full[nu]) for nu in nus}
    v_ref = full[INF].params["v"]            # L2 reference for sensitivity

    insample: dict[float, dict] = {}
    sens_rows: list[dict] = []
    v_rows: list[dict] = []
    for nu in nus:
        m = full[nu]
        insample[nu] = dict(
            full_loglik=m.log_lik(), full_aic=m.aic(), full_bic=m.bic(),
            n_iter=m.fit_result.n_iter, converged=bool(m.fit_result.converged))
        sens_rows.append(dict(
            solver=solver, nu=nu, sigma2=float(m.params["sigma2"]),
            full_loglik=m.log_lik(), full_aic=m.aic(), full_bic=m.bic(),
            **_v_sensitivity(m.params["v"], v_ref)))
        for rid, vj in zip(m.data.race_ids, m.params["v"]):
            v_rows.append(dict(solver=solver, nu=nu, race_id=int(rid), v=float(vj)))
        print(f"   nu={nu:<5} loglik={m.log_lik():.3f}  aic={m.aic():.1f}  "
              f"bic={m.bic():.1f}  iters={m.fit_result.n_iter}", flush=True)

    # ── grid K-fold CV ───────────────────────────────────────────────
    train_fds = {f: subset_fitdata(fd, fold != f) for f in range(K)}
    fold_rows: list[dict] = []
    train_chains: dict[int, dict[float, dict]] = {}
    for f in range(K):
        print(f"-- fold {f} --", flush=True)
        chain = _continuation_chain(train_fds[f], solver, nu_grid, tol=tol, max_iter=max_iter)
        train_chains[f] = {nu: _uv(chain[nu]) for nu in nus}
        for nu in nus:
            tm = chain[nu]
            sc = heldout_logdensity(fd, fold == f, tm)
            fold_rows.append(dict(
                solver=solver, nu=nu, fold=f, train_loglik=tm.log_lik(),
                heldout_logdens=sc["sum_logdens"], heldout_per_cell=sc["mean_logdens"],
                heldout_rmse=sc["rmse"], n_test=sc["n_test"], n_orphan=sc["n_orphan"],
                n_iter=tm.fit_result.n_iter))
        best = max((r for r in fold_rows if r["fold"] == f), key=lambda r: r["heldout_logdens"])
        print(f"   best held-out nu={best['nu']}  logdens={best['heldout_logdens']:.1f}",
              flush=True)

    fold_df = pd.DataFrame(fold_rows)
    grid_agg = (fold_df.groupby("nu")
                .agg(cv_logdens=("heldout_logdens", "sum"),
                     n_test_total=("n_test", "sum"),
                     n_orphan_total=("n_orphan", "sum"),
                     heldout_rmse=("heldout_rmse", "mean"))
                .reset_index())
    grid_agg["cv_per_cell"] = grid_agg["cv_logdens"] / grid_agg["n_test_total"]

    sel_rows: list[dict] = []
    for nu in nus:
        a = grid_agg[grid_agg.nu == nu].iloc[0]
        sel_rows.append(dict(
            solver=solver, nu=nu, source="grid", **insample[nu],
            cv_logdens=float(a.cv_logdens), cv_per_cell=float(a.cv_per_cell),
            heldout_rmse=float(a.heldout_rmse), n_test_total=int(a.n_test_total)))
    grid_best_nu = float(grid_agg.loc[grid_agg.cv_logdens.idxmax(), "nu"])

    # ── Brent refinement on log(nu) around the finite-grid argmax ────
    finite_agg = grid_agg[np.isfinite(grid_agg.nu)]
    grid_best_finite = float(finite_agg.loc[finite_agg.cv_logdens.idxmax(), "nu"])
    brent_nu = None
    if not np.isfinite(grid_best_nu):
        print("-- Brent skipped: grid argmax is nu=inf (heavy tails not favored) --",
              flush=True)
    else:
        g = sorted(float(x) for x in nu_grid)
        i = g.index(grid_best_finite)
        lo, hi = g[max(0, i - 1)], g[min(len(g) - 1, i + 1)]
        print(f"-- Brent refine in [{lo:g}, {hi:g}] (log nu) --", flush=True)

        def neg_cv(log_nu: float) -> float:
            nu = float(np.exp(log_nu))
            tot, ntot = 0.0, 0
            for f in range(K):
                tm = _fit_at_nu(train_fds[f], nu, _nearest_grid_uv(nu, train_chains[f]),
                                solver, tol=tol, max_iter=max_iter)
                sc = heldout_logdensity(fd, fold == f, tm)
                tot += sc["sum_logdens"]; ntot += sc["n_test"]
            sel_rows.append(dict(solver=solver, nu=nu, source="brent",
                                 full_loglik=np.nan, full_aic=np.nan, full_bic=np.nan,
                                 n_iter=np.nan, converged=np.nan,
                                 cv_logdens=tot, cv_per_cell=tot / max(ntot, 1),
                                 heldout_rmse=np.nan, n_test_total=ntot))
            print(f"     brent nu={nu:.3f}  cv_logdens={tot:.2f}", flush=True)
            return -tot

        res = minimize_scalar(neg_cv, bounds=(np.log(lo), np.log(hi)),
                              method="bounded", options={"xatol": BRENT_XATOL})
        brent_nu = float(np.exp(res.x))
        print(f"   Brent nu* = {brent_nu:.3f}  cv_logdens={-res.fun:.2f}", flush=True)

    sel_df = pd.DataFrame(sel_rows)
    selected_nu = float(sel_df.loc[sel_df.cv_logdens.idxmax(), "nu"])
    sel_df["is_selected"] = np.isclose(sel_df["nu"], selected_nu)

    return dict(sel_df=sel_df, fold_df=fold_df, sens_df=pd.DataFrame(sens_rows),
                v_df=pd.DataFrame(v_rows), grid_best_nu=grid_best_nu,
                brent_nu=brent_nu, selected_nu=selected_nu)


def run_slice(name, *, mrc, date_lo, date_hi, min_runner, cv_solvers,
              nu_grid, K, fold_seed, min_test_n, tol, max_iter):
    spec = S.build_spec(name, min_race_count=mrc, date_lo=date_lo, date_hi=date_hi,
                        min_runner=min_runner)
    slug = S.slug(spec)
    out_dir = OUT_ROOT / slug
    out_dir.mkdir(parents=True, exist_ok=True)

    fd = load_slice(spec)
    fold = assign_folds(fd, K, seed=fold_seed, min_test_n=min_test_n)
    n_testable = int((fold >= 0).sum())
    print(f"\n=== {name}  ({slug})  K={K} solvers={cv_solvers} ===", flush=True)
    print(f"    I={fd.I:,} J={fd.J:,} N={fd.N:,}  testable cells={n_testable:,}", flush=True)

    sel_all, fold_all, sens_all, v_all = [], [], [], []
    headline: list[dict] = []
    for solver in cv_solvers:
        r = run_solver(solver, fd, fold, K, nu_grid, tol=tol, max_iter=max_iter)
        sel_all.append(r["sel_df"]); fold_all.append(r["fold_df"])
        sens_all.append(r["sens_df"]); v_all.append(r["v_df"])
        headline.append(dict(solver=solver, grid_best_nu=r["grid_best_nu"],
                             brent_nu=r["brent_nu"], selected_nu=r["selected_nu"]))

    sel_df = pd.concat(sel_all, ignore_index=True)
    sel_df.insert(0, "slug", slug)
    sel_df.insert(0, "slice", name)
    sel_df.to_csv(out_dir / "nu_selection.csv", index=False)
    pd.concat(fold_all, ignore_index=True).to_csv(out_dir / "cv_folds.csv", index=False)
    pd.concat(sens_all, ignore_index=True).to_csv(out_dir / "param_sensitivity.csv", index=False)
    pd.concat(v_all, ignore_index=True).to_csv(out_dir / "v_xnu.csv", index=False)

    # selected_nu.csv keyed on the PRIMARY solver (first in cv_solvers) — the one
    # 03/e02_fit_baseline_t reads (ALS and Anderson reach the same fixed point).
    h0 = headline[0]
    sel_one = pd.DataFrame([dict(slice=name, slug=slug, solver=h0["solver"],
                                 selected_nu=h0["selected_nu"],
                                 grid_best_nu=h0["grid_best_nu"],
                                 brent_nu=h0["brent_nu"], K=K, fold_seed=fold_seed)])
    sel_one.to_csv(out_dir / "selected_nu.csv", index=False)

    print(f"\nWrote {out_dir}", flush=True)
    for h in headline:
        print(f"  [{h['solver']}] grid-best nu={h['grid_best_nu']}  "
              f"Brent nu={h['brent_nu']}  => selected nu={h['selected_nu']}", flush=True)


def rebuild_rollups() -> None:
    """Rebuild the cross-slice rollups by globbing the per-slice CSVs on disk.

    Idempotent and order-independent: the rollups always reflect EVERY slice
    present under OUT_ROOT, so running the sweep for a new slice augments the
    rollups instead of overwriting them with only the current run.
    """
    for name, glob in (("nu_selection_all.csv", "*/nu_selection.csv"),
                       ("selected_nu_all.csv", "*/selected_nu.csv")):
        parts = [pd.read_csv(p) for p in sorted(OUT_ROOT.glob(glob))]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(OUT_ROOT / name, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", "--min-race-count", dest="mrc", type=int, default=None)
    S.add_spec_args(ap, with_mrc=False)
    ap.add_argument("--nu-grid", type=float, nargs="+", default=list(NU_GRID))
    ap.add_argument("--K", type=int, default=K)
    ap.add_argument("--seed", type=int, default=FOLD_SEED)
    ap.add_argument("--cv-solvers", nargs="+", default=["anderson", "cv"],
                    choices=["als", "anderson"],
                    help="solver(s) for CV (default anderson; pass both to cross-check).")
    ap.add_argument("--min-test-n", type=int, default=MIN_TEST_N)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    for name in names:
        run_slice(
            name, mrc=args.mrc, date_lo=args.date_lo, date_hi=args.date_hi,
            min_runner=args.min_runner, cv_solvers=args.cv_solvers,
            nu_grid=list(args.nu_grid), K=args.K, fold_seed=args.seed,
            min_test_n=args.min_test_n, tol=args.tol, max_iter=args.max_iter)

    # rebuilt from every per-slice CSV on disk, so a new slice augments (not
    # overwrites) the rollups.
    rebuild_rollups()
    print(f"\nRebuilt cross-slice rollups -> {OUT_ROOT}")
    print(f"Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
