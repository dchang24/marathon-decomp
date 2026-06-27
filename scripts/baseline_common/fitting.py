"""Shared fitting primitives for the rank-1 baseline essential fits.

The two ``scripts/03_model_fit/baseline`` scripts (L2 / Student-t) are identical
apart from the noise level ``nu``, so the whole "fit the absolute-best model and
record the ALS-vs-Anderson convergence comparison" routine lives here as
:func:`run_essential_fit`.

Outputs (per slice ``{slug}`` = ``ALL_B_14-25_mrc2``):
  * essential model  -> results/models/{slug}/{stem}__{hash}/   (registry)
      the single overall-best fit across inits AND solvers (ALS/Anderson share
      the fixed point), registered with resample_tag="base" so each slice keeps
      exactly one L2 + one Student-t model.
  * convergence data -> results/convergence/{slug}/{stage}/      (CSV + parquet)
      init_spread.csv, loglik_traces.csv, v_agreement.csv,
      convergence_summary.csv, and the large v_traces.parquet.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from marathon_decomp import (
    AndersonFitterConfig,
    FitterConfig,
    ModelAnderson,
    Model,
    ModelConfig,
    SaveSpec,
    load_slice,
    registry,
)
from marathon_decomp.config import RESULTS_DIR

from baseline_common import inits as init_schemes

SOLVERS = ("als", "anderson")
MODELS_ROOT = RESULTS_DIR / "models"
CONV_ROOT = RESULTS_DIR / "convergence"


def baseline_cfg(nu: float) -> ModelConfig:
    """Rank-1 baseline (u_i + v_j); every optional term off, noise level nu."""
    return ModelConfig(use_phi12=False, use_gamma=False, use_d=False, nu=float(nu))


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _fitter(solver: str, ws: dict, *, max_iter: int, tol: float):
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init="warmstart", warmstart=ws, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def load_l2_warmstart(spec, fd) -> dict | None:
    """The slice's registered Gaussian (L2) best fit as a {u, v} warmstart.

    Reconstructs the deterministic registry path
    ``results/models/{slug}/baseline_L2_best__<hash>/`` and reloads its params.
    Returns None if that fit hasn't been produced yet (so callers can fall back
    to a cold start).
    """
    cfg2 = baseline_cfg(float("inf"))
    parent = MODELS_ROOT / registry.slice_slug(spec)
    fit_dir = registry.fit_path(parent, registry.model_stem(cfg2, "best"),
                                spec, cfg2, resample_tag="base")
    if not (fit_dir / "fit.pkl").is_file():
        return None
    m = registry.load_fit(fit_dir, fd)
    return {"u": np.asarray(m.params["u"]).copy(), "v": np.asarray(m.params["v"]).copy()}


def _stepping_warmstart(fd, base_ws: dict, stepping_nu: float, *,
                        tol: float, max_iter: int, solver: str = "anderson") -> dict:
    """Fit the baseline at ``stepping_nu`` warm-started from ``base_ws`` (the L2
    solution) and return its (u, v) — the 'stepping stone' init for the target nu."""
    model = _model_cls(solver)(fd, baseline_cfg(stepping_nu),
                               _fitter(solver, base_ws, max_iter=max_iter, tol=tol))
    model.fit()
    return {"u": model.params["u"].copy(), "v": model.params["v"].copy()}


def _augment_with_warm(init_list, spec, fd, *, n_random, jit_u, jit_v, seed0,
                       stepping_nu, tol, max_iter):
    """Append the L2-based warm-start strategies to the cold ``init_list``:
      * ``l2_warm``        — anchor at the L2 (u, v).
      * ``l2_rand0..k``    — perturbations around the L2 anchor.
      * ``step{nu}``       — fit ``stepping_nu`` warm from L2, then use that as
                             the warmstart for the target nu (the stepping stone).
    No-op (cold only) if the slice has no registered L2 fit yet.
    """
    l2 = load_l2_warmstart(spec, fd)
    if l2 is None:
        print("    [warn] warm_from_l2 set but no L2 fit on disk -> cold start only",
              flush=True)
        return init_list
    print(f"    + L2 warm strategies: l2_warm, l2_rand0..{n_random - 1}, "
          f"step{stepping_nu:g}", flush=True)
    out = list(init_list)
    out.append(("l2_warm", None, {"u": l2["u"].copy(), "v": l2["v"].copy()}))
    for k in range(n_random):
        seed = seed0 + 100 + k
        out.append((f"l2_rand{k}", seed,
                    init_schemes.perturb(l2, seed, jit_u=jit_u, jit_v=jit_v)))
    step = _stepping_warmstart(fd, l2, stepping_nu, tol=tol, max_iter=max_iter)
    out.append((f"step{stepping_nu:g}", None, step))
    return out


def _fit_capture(model, *, v_thin: int):
    """Fit `model`, capturing per-iter loglik/rss (all) and v (thinned).

    Returns (fit_result, v_arr[(T,J)], v_its, ll[], rss[], it_idx[], wall_s).
    """
    v_iters: list[np.ndarray] = []
    v_its: list[int] = []
    ll: list[float] = []
    rss: list[float] = []
    it_idx: list[int] = []

    def hook(m, it, loglik, rss_):
        ll.append(float(loglik))
        rss.append(float(rss_))
        it_idx.append(int(it))
        if it % v_thin == 0:
            v_iters.append(m.params["v"].copy())
            v_its.append(int(it))

    model.iter_hook = hook
    t0 = time.perf_counter()
    res = model.fit()
    dt = time.perf_counter() - t0
    v_arr = np.array(v_iters) if v_iters else np.empty((0, model.data.J))
    return res, v_arr, v_its, np.array(ll), np.array(rss), np.array(it_idx), dt


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def run_essential_fit(
    spec,
    nu: float,
    *,
    stage: str,
    n_random: int = 5,
    jit_u: float = 0.10,
    jit_v: float = 0.10,
    seed0: int = 0,
    max_iter: int = 2000,
    tol: float = 1e-12,
    v_thin: int = 1,
    keep_v_traces: bool = True,
    register: bool = True,
    study: str = "baseline",
    warm_from_l2: bool = False,
    stepping_nu: float = 15.0,
) -> dict:
    """Fit the absolute-best baseline at `nu` and record convergence diagnostics.

    `stage` is the convergence sub-dir label ("L2" or "nu_selected"). Fits
    mean + `n_random` perturbed restarts under BOTH solvers at tolerance `tol`,
    keeps the single overall-best (max loglik) as the registered essential model,
    and writes the full ALS-vs-Anderson convergence comparison.

    When `warm_from_l2` is set, the cold init set is augmented with warm-start
    strategies built off the slice's registered L2 fit (``l2_warm``, ``l2_rand*``,
    and a ``step{stepping_nu}`` stepping stone); see :func:`_augment_with_warm`.
    Falls back to cold start (with a warning) if no L2 fit is on disk.

    Returns a dict with the best model, its (solver, init), and the path written.
    """
    cfg = baseline_cfg(nu)
    sl = registry.slice_slug(spec)
    print(f"\n=== {sl}  nu={nu}  (stage={stage}) ===", flush=True)
    fd = load_slice(spec)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    # cold start (mean + random restarts) — always present, the default that
    # needs nothing precomputed.
    init_list = init_schemes.build_inits(
        fd, n_random=n_random, jit_u=jit_u, jit_v=jit_v, seed0=seed0)
    # warm strategies built off the already-fitted L2 model (skip for L2 itself).
    if warm_from_l2:
        init_list = _augment_with_warm(
            init_list, spec, fd, n_random=n_random, jit_u=jit_u, jit_v=jit_v,
            seed0=seed0, stepping_nu=stepping_nu, tol=tol, max_iter=max_iter)

    spread_rows: list[dict] = []
    ll_rows: list[dict] = []
    v_rows: list[dict] = []
    final_v: dict[tuple[str, str], np.ndarray] = {}
    # overall best across BOTH solvers + inits.
    best: tuple[float, str, str, int | None, object] | None = None  # (ll, solver, init, seed, model)

    for solver in SOLVERS:
        ModelCls = _model_cls(solver)
        for init_name, seed, ws in init_list:
            model = ModelCls(fd, cfg, _fitter(solver, ws, max_iter=max_iter, tol=tol))
            res, v_arr, v_its, ll, rss, it_idx, dt = _fit_capture(model, v_thin=v_thin)
            mean_v = float(model.params["v"].mean())
            final_v[(solver, init_name)] = model.params["v"].copy()

            spread_rows.append(dict(
                slice=sl, solver=solver, init=init_name, seed=seed,
                n_iter=res.n_iter, converged=res.converged,
                loglik_final=res.loglik_final, rss_final=res.rss_final,
                wall_s=dt, mean_v=mean_v,
            ))
            for t, l, r in zip(it_idx, ll, rss):
                ll_rows.append(dict(slice=sl, solver=solver, init=init_name,
                                    iter=int(t), loglik=float(l), rss=float(r)))
            if keep_v_traces:
                for t, vvec in zip(v_its, v_arr):
                    for j, vj in enumerate(vvec):
                        v_rows.append(dict(
                            slice=sl, solver=solver, init=init_name, iter=int(t),
                            race_idx=j, race_id=int(fd.race_ids[j]), v=float(vj)))

            if best is None or res.loglik_final > best[0]:
                best = (res.loglik_final, solver, init_name, seed, model)

            conv = "OK " if res.converged else "MAX"
            print(f"    {solver:8s} {init_name:8s} iters={res.n_iter:4d} {conv}"
                  f"  loglik={res.loglik_final:.6f}  {dt:.2f}s", flush=True)

    spread_df = pd.DataFrame(spread_rows)
    # loglik gap to the global best (across solvers + inits).
    spread_df["loglik_gap_to_best"] = spread_df["loglik_final"].max() - spread_df["loglik_final"]

    # converged-v agreement across inits, within each solver (pairwise).
    agree_rows: list[dict] = []
    for solver in SOLVERS:
        keys = [k for k in final_v if k[0] == solver]
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                va, vb = final_v[keys[a]], final_v[keys[b]]
                agree_rows.append(dict(
                    slice=sl, solver=solver, init_a=keys[a][1], init_b=keys[b][1],
                    corr_v=_safe_corr(va, vb), max_abs_dv=float(np.max(np.abs(va - vb)))))
    # cross-solver agreement of the two solvers' best-init v.
    best_v_by_solver: dict[str, np.ndarray] = {}
    for solver in SOLVERS:
        sub = spread_df[spread_df.solver == solver]
        bi = sub.loc[sub.loglik_final.idxmax(), "init"]
        best_v_by_solver[solver] = final_v[(solver, bi)]
    cross_corr = _safe_corr(best_v_by_solver["als"], best_v_by_solver["anderson"])
    cross_dv = float(np.max(np.abs(best_v_by_solver["als"] - best_v_by_solver["anderson"])))

    # per-solver convergence rollup (the ALS-vs-Anderson comparison).
    conv_rows: list[dict] = []
    for solver in SOLVERS:
        sub = spread_df[spread_df.solver == solver]
        conv_rows.append(dict(
            slice=sl, solver=solver, nu=float(nu),
            loglik_best=float(sub.loglik_final.max()),
            n_iter_best=int(sub.loc[sub.loglik_final.idxmax(), "n_iter"]),
            n_iter_median=float(sub["n_iter"].median()),
            n_iter_max=int(sub["n_iter"].max()),
            wall_s_best=float(sub.loc[sub.loglik_final.idxmax(), "wall_s"]),
            wall_s_median=float(sub["wall_s"].median()),
            all_converged=bool(sub["converged"].all()),
            loglik_gap_across_inits=float(sub["loglik_final"].max() - sub["loglik_final"].min()),
            max_abs_dv_across_inits=float(
                max((r["max_abs_dv"] for r in agree_rows if r["solver"] == solver), default=0.0)),
            cross_solver_corr_v=cross_corr, cross_solver_max_dv=cross_dv,
        ))
    conv_df = pd.DataFrame(conv_rows)

    # ── write convergence diagnostics ────────────────────────────────
    conv_dir = CONV_ROOT / sl / stage
    conv_dir.mkdir(parents=True, exist_ok=True)
    spread_df.to_csv(conv_dir / "init_spread.csv", index=False)
    pd.DataFrame(ll_rows).to_csv(conv_dir / "loglik_traces.csv", index=False)
    conv_df.to_csv(conv_dir / "convergence_summary.csv", index=False)
    if agree_rows:
        pd.DataFrame(agree_rows).to_csv(conv_dir / "v_agreement.csv", index=False)
    if keep_v_traces and v_rows:
        # the one genuinely large table -> parquet (raw data for later plotting).
        pd.DataFrame(v_rows).to_parquet(conv_dir / "v_traces.parquet", index=False)

    ll_best, s_best, i_best, seed_best, model = best
    for _, r in conv_df.iterrows():
        print(f"    [{r.solver}] n_iter(best/med/max)={r.n_iter_best}/{r.n_iter_median:g}/"
              f"{r.n_iter_max}  loglik gap across inits={r.loglik_gap_across_inits:.2e}",
              flush=True)
    print(f"    overall best: solver={s_best} init={i_best}  loglik={ll_best:.6f}", flush=True)
    print(f"    cross-solver v: corr={cross_corr:.6f}  max|dv|={cross_dv:.2e}", flush=True)
    print(f"    convergence -> {conv_dir}", flush=True)

    # ── register the single overall-best model ───────────────────────
    fit_dir = None
    if register:
        parent = MODELS_ROOT / sl
        stem = registry.model_stem(cfg, "best")   # baseline_L2_best / baseline_nu5p00_best
        fit_dir = registry.fit_path(parent, stem, spec, cfg, resample_tag="base")
        registry.register_fit(
            model, fit_dir, resample_tag="base", save=SaveSpec(params=True),
            study=study, best_solver=s_best, best_init=i_best,
            best_init_seed=(int(seed_best) if seed_best is not None else None),
            n_inits=len(init_list), tol=tol, stage=stage,
            cross_solver_corr_v=cross_corr, cross_solver_max_dv=cross_dv,
        )
        print(f"    model -> {fit_dir}", flush=True)

    return dict(model=model, best_solver=s_best, best_init=i_best,
                loglik_best=ll_best, fit_dir=fit_dir, conv_dir=conv_dir,
                spread_df=spread_df, conv_df=conv_df)
