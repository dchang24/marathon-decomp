"""Aging-form grid fit (model selection, stage A1).

Fit the aging-only model (rank-1 ``u+v`` + phi-block + optional gamma; ``d`` off)
for every cell of the candidate grid

    basis (poly2..6, spline3..6)  x  gamma_form (off, scalar, varying)
                                  x  nu (inf, 8.0)
                                  x  solver (anderson [+ als])
                                  x  init (mean + N random restarts)

and record, per cell:
  * in-sample fit quality -- loglik, AIC, BIC (effective d.o.f.), naive n_params,
  * convergence -- n_iter, converged flag, wall time, and the init spread
    (loglik gap + max |dv| across inits: does the form converge to an
    init-invariant fixed point? splines with the `direct` inner solver are the
    thing to validate),
  * the reconstructed aging curve theta_aging @ B(A_n) on a common grid plus the
    gamma entry-age curve, so curve inspection needs no refit,
  * the best-init full-data v per race, for cross-form v comparison (q02).

This stage does NOT select a winner -- that is the held-out CV job of
``e02_aging_cv``. It only fits, scores in-sample, and reconstructs curves.

PERSISTENCE: this is the expensive step. In addition to the tidy tables below it
saves the full param payload of every best-init cell to ``grid/fits/*.pkl`` (the
source of truth); q01/q02 and any future metric recompute from these, never
re-fit.

Outputs (per slice) under ``results/model_selection/aging/{slug}/grid/``:
  * metrics.csv        -- one row per (nu, cand, gamma_form, solver, init).
  * convergence.csv    -- per (nu, cand, gamma_form, solver) init-spread rollup.
  * curves.parquet     -- best-init aging + gamma curve on a common A_n grid.
  * an_density.parquet -- A_n histogram (the rug, for spotting sparse edges).
  * v_xform.parquet    -- best-init full-data v per (nu, cand, solver, race_id).
  * slice_info.csv     -- ae_mean / entry-age + A_n percentiles / I,J,N.
  * fits/{cand}_{nutag}_{solver}.pkl -- best-init full param payload per cell.
Cross-slice rollups (metrics.csv, convergence.csv) at the dir root.

Run::

    python scripts/02_model_selection/aging/e01_aging_grid.py --smoke
    python scripts/02_model_selection/aging/e01_aging_grid.py --slices Po10_M --gamma scalar
    python scripts/02_model_selection/aging/e01_aging_grid.py          # full grid, all slices
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
    aging_curve_on_grid,
    entry_age_curve_on_grid,
    load_slice,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from baseline_common import slices as S  # noqa: E402
from baseline_common.inits import build_inits  # noqa: E402
import grid as G  # noqa: E402

# -- tunables (CLI flags override) ------------------------------------------
TOL = 1e-10            # outer-loop stopping tolerance (relative loglik)
MAX_ITER = 2000        # outer-iteration cap
N_RANDOM = 2           # random-perturbed-around-mean inits (convergence check)
RAND_JIT_U = 0.10      # std of the Gaussian perturbation on u (log-seconds)
RAND_JIT_V = 0.10      # std of the Gaussian perturbation on v (log-seconds)
RAND_SEED0 = 0
N_GRID = 200           # points on the A_n reconstruction grid
N_AN_BINS = 50         # bins for the A_n density rug
AE_OFFSET_YR = 10.0    # entry-age offset (yr above mean) for the stored gamma curve
# ---------------------------------------------------------------------------

STAGE = "grid"
OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
SOLVERS_DEFAULT = ("anderson",)


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _make_fitter(solver: str, warmstart: dict, *, max_iter: int, tol: float):
    # inner_solver defaults to "direct" -- exact dense solve, required for the
    # spline bases and fine for poly. We do not override it.
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init="warmstart", warmstart=warmstart, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def _curve_rows(model, *, A_grid, ae_off, slug, nu, cand, basis, gamma_form,
                solver, loglik, converged) -> list[dict]:
    """Reconstruct the aging curve (and gamma curve, if on) on `A_grid`."""
    aging = aging_curve_on_grid(model, A_grid)
    if model.config.use_gamma:
        gamma = entry_age_curve_on_grid(model, A_grid, ae_off)
    else:
        gamma = np.full_like(A_grid, np.nan)
    return [dict(
        slug=slug, nu=G.nu_label(nu), cand=cand, basis=basis, gamma_form=gamma_form,
        solver=solver, A_n=float(a), aging_curve=float(yc), gamma_curve=float(gc),
        ae_offset_yr=ae_off, loglik=loglik, converged=converged,
    ) for a, yc, gc in zip(A_grid, aging, gamma)]


def run_slice(name, *, bases, gamma_forms, nus, solvers, n_random, jit_u, jit_v,
              seed0, max_iter, tol, spec_over) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = S.build_spec(name, **spec_over)
    slug = S.slug(spec)
    print(f"\n=== {name}  ({slug}) ===", flush=True)
    fd = load_slice(spec)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    init_list = build_inits(fd, n_random=n_random, jit_u=jit_u, jit_v=jit_v, seed0=seed0)

    # common reconstruction grid + A_n density rug (slice-level, write once).
    A_n = fd.A_n[np.isfinite(fd.A_n)]
    A_max = float(A_n.max()) if A_n.size else 1.0
    an_p95 = float(np.percentile(A_n, 95)) if A_n.size else A_max
    A_grid = np.linspace(0.0, A_max, N_GRID)
    counts, edges = np.histogram(A_n, bins=N_AN_BINS)
    A_e = fd.A_e[np.isfinite(fd.A_e)]
    ae_mean = float(A_e.mean()) if A_e.size else float("nan")

    out_dir = OUT_ROOT / slug / STAGE
    fits_dir = out_dir / "fits"
    fits_dir.mkdir(parents=True, exist_ok=True)

    metric_rows: list[dict] = []
    conv_rows: list[dict] = []
    curve_rows: list[dict] = []
    v_rows: list[dict] = []

    for nu in nus:
        for basis in bases:
            for gamma_form in gamma_forms:
                cand = G.cand_label(basis, gamma_form)
                cfg = G.build_config(basis, gamma_form, nu)
                for solver in solvers:
                    ModelCls = _model_cls(solver)
                    final_v: dict[str, np.ndarray] = {}
                    lls: list[float] = []
                    all_conv = True
                    best: tuple[float, str, object] | None = None  # (ll, init, model)
                    for init_name, seed, ws in init_list:
                        model = ModelCls(fd, cfg, _make_fitter(solver, ws, max_iter=max_iter, tol=tol))
                        t0 = time.perf_counter()
                        res = model.fit()
                        dt = time.perf_counter() - t0
                        final_v[init_name] = model.params["v"].copy()
                        lls.append(model.log_lik())
                        all_conv = all_conv and bool(res.converged)
                        metric_rows.append(dict(
                            slug=slug, slice=name, nu=G.nu_label(nu), cand=cand,
                            basis=basis.name, gamma_form=gamma_form, solver=solver,
                            init=init_name, seed=seed,
                            loglik=model.log_lik(), aic=model.aic(), bic=model.bic(),
                            n_params=model.n_params_naive(),
                            n_iter=res.n_iter, converged=bool(res.converged), wall_s=dt))
                        if best is None or res.loglik_final > best[0]:
                            best = (res.loglik_final, init_name, model)

                    ll_best, init_best, model_best = best
                    vmat = np.array(list(final_v.values()))
                    max_dv = float(np.max(np.abs(vmat - vmat[0]))) if len(vmat) > 1 else 0.0
                    conv_rows.append(dict(
                        slug=slug, slice=name, nu=G.nu_label(nu), cand=cand,
                        basis=basis.name, gamma_form=gamma_form, solver=solver,
                        loglik_best=ll_best, best_init=init_best,
                        loglik_gap_across_inits=float(max(lls) - min(lls)),
                        max_abs_dv_across_inits=max_dv, all_converged=all_conv))

                    curve_rows.extend(_curve_rows(
                        model_best, A_grid=A_grid, ae_off=AE_OFFSET_YR, slug=slug,
                        nu=nu, cand=cand, basis=basis.name, gamma_form=gamma_form,
                        solver=solver, loglik=ll_best,
                        converged=bool(model_best.fit_result.converged)))
                    for rid, vj in zip(model_best.data.race_ids, model_best.params["v"]):
                        v_rows.append(dict(
                            slug=slug, nu=G.nu_label(nu), cand=cand, basis=basis.name,
                            gamma_form=gamma_form, solver=solver, race_id=int(rid),
                            v=float(vj)))

                    # source-of-truth payload for post-hoc recompute (no re-fit).
                    model_best.save(fits_dir / f"{G.fit_stem(basis, gamma_form, nu, solver)}.pkl",
                                    what=SaveSpec(params=True))

                    print(f"    nu={G.nu_label(nu):>3s}  {cand:18s} {solver:8s}"
                          f"  ll={ll_best:.3f}  bic={model_best.bic():.0f}"
                          f"  gap={float(max(lls) - min(lls)):.1e}  max|dv|={max_dv:.1e}",
                          flush=True)

    metrics_df = pd.DataFrame(metric_rows)
    conv_df = pd.DataFrame(conv_rows)
    pd.DataFrame(curve_rows).to_parquet(out_dir / "curves.parquet", index=False)
    pd.DataFrame(dict(slug=slug, bin_lo=edges[:-1], bin_hi=edges[1:], count=counts)
                 ).to_parquet(out_dir / "an_density.parquet", index=False)
    pd.DataFrame(v_rows).to_parquet(out_dir / "v_xform.parquet", index=False)
    metrics_df.to_csv(out_dir / "metrics.csv", index=False)
    conv_df.to_csv(out_dir / "convergence.csv", index=False)
    pd.DataFrame([dict(
        slug=slug, slice=name, I=fd.I, J=fd.J, N=fd.N, ae_mean=ae_mean,
        ae_p10=float(np.percentile(A_e, 10)) if A_e.size else float("nan"),
        ae_p50=float(np.percentile(A_e, 50)) if A_e.size else float("nan"),
        ae_p90=float(np.percentile(A_e, 90)) if A_e.size else float("nan"),
        an_p95=an_p95, an_max=A_max)]).to_csv(out_dir / "slice_info.csv", index=False)

    print(f"    wrote {out_dir}  (+ {len(list(fits_dir.glob('*.pkl')))} payloads)", flush=True)
    return metrics_df, conv_df


def rebuild_rollups() -> None:
    """Rebuild cross-slice rollups by globbing the per-slice CSVs on disk."""
    for name, glob in (("metrics.csv", "*/grid/metrics.csv"),
                       ("convergence.csv", "*/grid/convergence.csv")):
        parts = [pd.read_csv(p) for p in sorted(OUT_ROOT.glob(glob))]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(OUT_ROOT / name, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--bases", nargs="+", default=["all"],
                    help="basis names (poly2..6, spline3..6) or 'all'.")
    ap.add_argument("--gamma", nargs="+", default=["all"],
                    choices=[*G.GAMMA_FORMS, "all"], help="gamma forms or 'all'.")
    ap.add_argument("--nu", nargs="+", default=["all"],
                    help="nu values ('inf', '8', ...) or 'all'.")
    ap.add_argument("--solvers", nargs="+", default=list(SOLVERS_DEFAULT),
                    choices=["als", "anderson"])
    ap.add_argument("--n-random", type=int, default=N_RANDOM)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--tol", type=float, default=TOL)
    S.add_spec_args(ap)
    ap.add_argument("--smoke", action="store_true",
                    help="Po10_W, poly2/spline4, gamma=scalar, nu=8, 1 random.")
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    bases = G.BASES if args.bases == ["all"] else tuple(G.BASES_BY_NAME[b] for b in args.bases)
    gforms = G.GAMMA_FORMS if args.gamma == ["all"] else tuple(args.gamma)
    nus = G.NU_GRID if args.nu == ["all"] else tuple(float(x) for x in args.nu)
    n_random = args.n_random

    if args.smoke:
        names = ["Po10_W"]
        bases = (G.BASES_BY_NAME["poly2"], G.BASES_BY_NAME["spline4"])
        gforms = ("scalar",)
        nus = (8.0,)
        n_random = 1

    spec_over = dict(min_race_count=args.min_race_count, date_lo=args.date_lo,
                     date_hi=args.date_hi, min_runner=args.min_runner)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    for name in names:
        run_slice(name, bases=bases, gamma_forms=gforms, nus=nus, solvers=tuple(args.solvers),
                  n_random=n_random, jit_u=RAND_JIT_U, jit_v=RAND_JIT_V, seed0=RAND_SEED0,
                  max_iter=args.max_iter, tol=args.tol, spec_over=spec_over)

    rebuild_rollups()
    print(f"\nRebuilt cross-slice rollups -> {OUT_ROOT}")
    print(f"Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
