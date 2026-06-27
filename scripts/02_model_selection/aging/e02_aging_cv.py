"""K-fold CV selection of the aging form (model selection, stage A2).

Picks the aging parametric form (basis x gamma) by held-out predictive density,
at fixed nu in {inf, 8.0}. In-sample loglik is monotone in flexibility (more
columns never hurt the fit), so the form is chosen out-of-sample; AIC/BIC are
reported alongside as a cheap cross-check.

For each candidate (basis x gamma_form) and each nu:
  * full-data fit  -- in-sample loglik / AIC / BIC (continuation inf -> 8 so the
    nu=8 fit starts from the converged Gaussian fixed point),
  * K-fold CV      -- fit on each fold's training cells (same continuation), then
    score the held-out cells with ``heldout_logdensity`` (which includes the
    train-fitted aging block, so the score sees the aging term).

Selection criterion = CV log predictive density per held-out cell, maximized.

CV AGE BASIS -- full-sample, fixed. Career age A_n and entry age A_e are computed
once on the whole slice and carried unchanged into every fold (the default
behaviour of ``subset_fitdata`` / ``heldout_logdensity``: ``A_n[test]`` is read
from the full ``fd``). A_n/A_e are covariates, not responses; theta/gamma/u/v are
still fit on training cells only, and the convention is identical across every
candidate form, so it cannot bias the relative form ranking.

PERSISTENCE: the full-data continuation fits are saved to ``cv/fits/*.pkl`` so the
in-sample side (loglik/AIC/BIC, curves, params) is recomputable without re-fit.
The per-fold training models are CV-internal and not persisted.

Outputs -> results/model_selection/aging/{slug}/cv/:
  * form_selection.csv -- one row per (nu, candidate): in-sample loglik/AIC/BIC +
    CV logdens (total + per-cell) + held-out rmse + n_test + n_orphan.
  * cv_folds.csv       -- per (nu, candidate, fold) held-out breakdown.
  * best_form.csv      -- argmax cv/cell per nu (convenience rollup).
  * fits/{cand}_{nutag}_{solver}.pkl -- full-data continuation fit per cell.

Run::

    python scripts/02_model_selection/aging/e02_aging_cv.py --slice Po10_W --bases poly2 spline4 --gamma scalar
    python scripts/02_model_selection/aging/e02_aging_cv.py --slice ALL_M
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/
sys.path.insert(0, str(Path(__file__).resolve().parent))       # this dir (grid.py)

from marathon_decomp import (  # noqa: E402
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    SaveSpec,
    assign_folds,
    heldout_logdensity,
    load_slice,
    subset_fitdata,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from baseline_common import slices as S  # noqa: E402
import grid as G  # noqa: E402

# -- tunables ---------------------------------------------------------------
DEFAULT_BASES = ["poly3", "poly4", "poly5", "poly6", "spline3", "spline4", "spline5"]
DEFAULT_GAMMA = ["scalar", "varying"]          # add "off" via --gamma to floor it
K = 5
FOLD_SEED = 0
MIN_TEST_N = 1
TOL = 1e-10
MAX_ITER = 2000
SOLVER = "anderson"
# ---------------------------------------------------------------------------

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _fitter(solver: str, *, init="mean", warmstart=None, tol=TOL, max_iter=MAX_ITER):
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init=init, warmstart=warmstart, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def _warm(m) -> dict[str, np.ndarray]:
    """Warm-start payload carrying the aging params too (s/d are off here)."""
    return {k: np.asarray(m.params[k]).copy()
            for k in ("u", "v", "theta_aging", "gamma")}


def _order_nus(nus) -> list[float]:
    """inf first, then finite descending -- the continuation order."""
    finite = sorted((x for x in nus if np.isfinite(x)), reverse=True)
    head = [float("inf")] if any(not np.isfinite(x) for x in nus) else []
    return head + finite


def _continuation(fd_, basis, gamma_form, nus, solver, *, tol, max_iter) -> dict[float, object]:
    """Fit the candidate down the nu order, warm-starting each nu from the prev."""
    out: dict[float, object] = {}
    prev = None
    for nu in nus:
        cfg = G.build_config(basis, gamma_form, nu)
        fcfg = (_fitter(solver, init="mean", tol=tol, max_iter=max_iter) if prev is None
                else _fitter(solver, init="warmstart", warmstart=prev, tol=tol, max_iter=max_iter))
        m = _model_cls(solver)(fd_, cfg, fcfg)
        m.fit()
        out[nu] = m
        prev = _warm(m)
    return out


def run(slice_name, bases, gammas, nus, *, K, seed, min_test_n, solver,
        tol, max_iter, spec_over):
    spec = S.build_spec(slice_name, **spec_over)
    slug = S.slug(spec)
    fd = load_slice(spec)
    fold = assign_folds(fd, K, seed=seed, min_test_n=min_test_n)
    n_testable = int((fold >= 0).sum())
    print(f"=== {slice_name}  ({slug})  aging-CV (K={K}, solver={solver}) ===", flush=True)
    print(f"    I={fd.I:,} J={fd.J:,} N={fd.N:,}  testable={n_testable:,}", flush=True)
    nus = _order_nus(nus)
    train_fds = {f: subset_fitdata(fd, fold != f) for f in range(K)}

    out_dir = OUT_ROOT / slug / "cv"
    fits_dir = out_dir / "fits"
    fits_dir.mkdir(parents=True, exist_ok=True)

    cands = [(G.BASES_BY_NAME[b], g) for b in bases for g in gammas]
    sel_rows: list[dict] = []
    fold_rows: list[dict] = []
    t0 = time.perf_counter()

    for basis, gamma_form in cands:
        cand = G.cand_label(basis, gamma_form)
        # full-data continuation -> in-sample metrics + persisted payloads.
        full = _continuation(fd, basis, gamma_form, nus, solver, tol=tol, max_iter=max_iter)
        insample = {nu: dict(
            full_loglik=full[nu].log_lik(), full_aic=full[nu].aic(),
            full_bic=full[nu].bic(), n_params=full[nu].n_params_naive(),
            n_iter=full[nu].fit_result.n_iter,
            converged=bool(full[nu].fit_result.converged)) for nu in nus}
        for nu in nus:
            full[nu].save(fits_dir / f"{G.fit_stem(basis, gamma_form, nu, solver)}.pkl",
                          what=SaveSpec(params=True))

        # K-fold CV
        cv = {nu: dict(logdens=0.0, n_test=0, n_orphan=0, rmse_sq=0.0) for nu in nus}
        for f in range(K):
            chain = _continuation(train_fds[f], basis, gamma_form, nus, solver,
                                  tol=tol, max_iter=max_iter)
            for nu in nus:
                tm = chain[nu]
                sc = heldout_logdensity(fd, fold == f, tm)
                fold_rows.append(dict(
                    slug=slug, solver=solver, nu=G.nu_label(nu), cand=cand, basis=basis.name,
                    gamma_form=gamma_form, fold=f, train_loglik=tm.log_lik(),
                    heldout_logdens=sc["sum_logdens"], heldout_per_cell=sc["mean_logdens"],
                    heldout_rmse=sc["rmse"], n_test=sc["n_test"], n_orphan=sc["n_orphan"],
                    n_iter=tm.fit_result.n_iter))
                cv[nu]["logdens"] += sc["sum_logdens"]
                cv[nu]["n_test"] += sc["n_test"]
                cv[nu]["n_orphan"] += sc["n_orphan"]
                cv[nu]["rmse_sq"] += sc["rmse"] ** 2 * sc["n_test"]

        for nu in nus:
            c = cv[nu]
            nt = max(c["n_test"], 1)
            sel_rows.append(dict(
                slug=slug, slice=slice_name, solver=solver, nu=G.nu_label(nu), cand=cand,
                basis=basis.name, gamma_form=gamma_form, **insample[nu],
                cv_logdens=c["logdens"], cv_per_cell=c["logdens"] / nt,
                heldout_rmse=float(np.sqrt(c["rmse_sq"] / nt)),
                n_test_total=c["n_test"], n_orphan_total=c["n_orphan"]))
            print(f"    {cand:18s} nu={G.nu_label(nu):>3s}  ll={insample[nu]['full_loglik']:.1f}"
                  f"  bic={insample[nu]['full_bic']:.0f}  cv/cell={c['logdens'] / nt:+.5f}",
                  flush=True)

    sel_df = pd.DataFrame(sel_rows)
    fold_df = pd.DataFrame(fold_rows)
    sel_df.to_csv(out_dir / "form_selection.csv", index=False)
    fold_df.to_csv(out_dir / "cv_folds.csv", index=False)

    best_rows = []
    for nu in nus:
        sub = sel_df[sel_df.nu == G.nu_label(nu)]
        b = sub.loc[sub.cv_per_cell.idxmax()]
        best_rows.append(dict(slug=slug, slice=slice_name, solver=solver, nu=G.nu_label(nu),
                              best_cand=b["cand"], cv_per_cell=float(b.cv_per_cell),
                              full_bic=float(b.full_bic), K=K, fold_seed=seed))
    pd.DataFrame(best_rows).to_csv(out_dir / "best_form.csv", index=False)

    print(f"\nWrote {out_dir}  ({time.perf_counter() - t0:.1f}s)")
    for r in best_rows:
        print(f"  nu={r['nu']:>3s}: best CV form = {r['best_cand']}  "
              f"(cv/cell={r['cv_per_cell']:+.5f}, bic={r['full_bic']:.0f})")
    return sel_df, fold_df


def rebuild_rollups() -> None:
    for name, glob in (("form_selection_all.csv", "*/cv/form_selection.csv"),
                       ("best_form_all.csv", "*/cv/best_form.csv")):
        parts = [pd.read_csv(p) for p in sorted(OUT_ROOT.glob(glob))]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(OUT_ROOT / name, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", "--slices", dest="slices", nargs="+", default=["Po10_M"])
    ap.add_argument("--bases", nargs="+", default=DEFAULT_BASES)
    ap.add_argument("--gamma", nargs="+", default=DEFAULT_GAMMA,
                    choices=["off", "scalar", "varying"])
    ap.add_argument("--nu", nargs="+", default=None,
                    help="nu values ('inf', '8', ...); default inf + 8.0.")
    ap.add_argument("--K", type=int, default=K)
    ap.add_argument("--seed", type=int, default=FOLD_SEED)
    ap.add_argument("--min-test-n", type=int, default=MIN_TEST_N)
    ap.add_argument("--solver", default=SOLVER, choices=["als", "anderson"])
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    S.add_spec_args(ap)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    nus = G.NU_GRID if args.nu is None else tuple(float(x) for x in args.nu)
    spec_over = dict(min_race_count=args.min_race_count, date_lo=args.date_lo,
                     date_hi=args.date_hi, min_runner=args.min_runner)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for name in names:
        run(name, args.bases, args.gamma, nus, K=args.K, seed=args.seed,
            min_test_n=args.min_test_n, solver=args.solver, tol=args.tol,
            max_iter=args.max_iter, spec_over=spec_over)
    rebuild_rollups()
    print(f"\nRebuilt cross-slice rollups -> {OUT_ROOT}")


if __name__ == "__main__":
    main()
