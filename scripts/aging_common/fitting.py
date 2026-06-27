"""Shared fitting routine for the production aging-block essential fit.

Fits the **production aging form** (rank-1 ``u_i + v_j`` + parametric aging
block; per-athlete drift ``d_i`` off) at the settled operating point and keeps
the single overall-best fit (max loglik across inits) as the registered model.

Production aging form (the default everywhere downstream):
  natural cubic spline, **4 knots**, **varying** gamma, **nu=8**, Anderson.

The fit is **warm-started from the slice's registered baseline nu=8 model**
(``baseline_nu8p00_best``): its ``(u, v)`` seed the anchor init, with the aging
coefficients zero-filled (the fitter does that). ``--n-random`` perturbed
restarts around that anchor probe init-invariance of the fixed point.

Outputs (per slice ``{slug}`` = ``ALL_B_14-25_mrc2``):
  * essential model  -> results/models/{slug}/{stem}__{hash}/        (registry)
      e.g. ``agingS4gv_nu8p00_best__<hash>`` (resample_tag="base").
  * convergence data -> results/convergence/{slug}/{stage}/          (CSV+parquet)
      init_spread.csv, loglik_traces.csv, v_agreement.csv,
      convergence_summary.csv, v_traces.parquet (per-iter v_j),
      aging_coef_traces.parquet (per-iter theta_aging + gamma), and
      curves.parquet (best-init reconstructed aging + gamma curve on a grid).
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
    aging_curve_on_grid,
    entry_age_curve_on_grid,
    load_slice,
    registry,
)
from marathon_decomp.config import RESULTS_DIR

from baseline_common import inits as init_schemes

MODELS_ROOT = RESULTS_DIR / "models"
CONV_ROOT = RESULTS_DIR / "convergence"

# entry-age offset (yr above mean) for the stored gamma curve, matching the
# model-selection grid stage so the curves are directly comparable.
AE_OFFSET_YR = 10.0
N_GRID = 200            # points on the A_n reconstruction grid


def aging_cfg(nu: float, *, n_knots: int = 4, gamma_form: str = "varying",
              use_d: bool = False) -> ModelConfig:
    """Production aging config: rank-1 + spline aging block + varying gamma.

    Defaults reproduce the settled production form (spline-4, varying gamma);
    ``d`` is off (the per-athlete drift is a later production stage).
    """
    return ModelConfig(
        use_phi12=True,
        use_gamma=True,
        use_d=use_d,
        basis_kind="spline",
        n_knots=n_knots,
        gamma_form=gamma_form,
        nu=float(nu),
    )


def aging_stem(cfg: ModelConfig) -> str:
    """Readable model-folder stem ``aging{form}{g}_{nutag}_best``.

    ``form`` = ``S{n_knots}`` for spline / ``P{degree}`` for poly; ``g`` =
    ``gv`` (varying) / ``gs`` (scalar) / ``g0`` (off). E.g. the production
    spline-4 varying-gamma nu=8 fit -> ``agingS4gv_nu8p00_best``.
    """
    if cfg.basis_kind == "spline":
        form = f"S{cfg.n_knots}"
    else:
        form = f"P{cfg.degree}"
    g = "g0" if not cfg.use_gamma else ("gv" if cfg.gamma_form == "varying" else "gs")
    return f"aging{form}{g}_{registry.nu_tag(float(cfg.nu))}_best"


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _fitter(solver: str, ws: dict, *, max_iter: int, tol: float):
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init="warmstart", warmstart=ws, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def load_baseline_warmstart(spec, fd, nu: float) -> dict | None:
    """The slice's registered baseline nu model as a {u, v} warmstart.

    Reconstructs the deterministic registry path
    ``results/models/{slug}/baseline_nu<p>_best__<hash>/`` (the rank-1 u+v
    model at the same nu) and reloads its params. Returns None if absent.
    """
    cfg_base = ModelConfig(use_phi12=False, use_gamma=False, use_d=False, nu=float(nu))
    parent = MODELS_ROOT / registry.slice_slug(spec)
    fit_dir = registry.fit_path(parent, registry.model_stem(cfg_base, "best"),
                                spec, cfg_base, resample_tag="base")
    if not (fit_dir / "fit.pkl").is_file():
        return None
    m = registry.load_fit(fit_dir, fd)
    return {"u": np.asarray(m.params["u"]).copy(), "v": np.asarray(m.params["v"]).copy()}


def _fit_capture(model, *, v_thin: int):
    """Fit `model`, capturing per-iter loglik/rss (all) and v + aging coefs
    (thinned).

    Returns (fit_result, v_arr[(T,J)], theta_arr[(T,K)], gamma_arr[(T,G)],
    cap_its, ll[], rss[], it_idx[], wall_s).
    """
    v_iters: list[np.ndarray] = []
    theta_iters: list[np.ndarray] = []
    gamma_iters: list[np.ndarray] = []
    cap_its: list[int] = []
    ll: list[float] = []
    rss: list[float] = []
    it_idx: list[int] = []

    def hook(m, it, loglik, rss_):
        ll.append(float(loglik))
        rss.append(float(rss_))
        it_idx.append(int(it))
        if it % v_thin == 0:
            v_iters.append(m.params["v"].copy())
            theta_iters.append(np.atleast_1d(m.params["theta_aging"]).copy())
            gamma_iters.append(np.atleast_1d(m.params["gamma"]).copy())
            cap_its.append(int(it))

    model.iter_hook = hook
    t0 = time.perf_counter()
    res = model.fit()
    dt = time.perf_counter() - t0
    v_arr = np.array(v_iters) if v_iters else np.empty((0, model.data.J))
    theta_arr = np.array(theta_iters) if theta_iters else np.empty((0, model.K_basis))
    gamma_arr = np.array(gamma_iters) if gamma_iters else np.empty((0, 0))
    return (res, v_arr, theta_arr, gamma_arr, cap_its,
            np.array(ll), np.array(rss), np.array(it_idx), dt)


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _curve_rows(model, *, A_grid, ae_off, slug, nu) -> list[dict]:
    """Best-init reconstructed aging curve (+ gamma curve) on `A_grid`."""
    aging = aging_curve_on_grid(model, A_grid)
    gamma = entry_age_curve_on_grid(model, A_grid, ae_off)
    return [dict(slug=slug, nu=float(nu), A_n=float(a),
                 aging_curve=float(yc), gamma_curve=float(gc), ae_offset_yr=ae_off)
            for a, yc, gc in zip(A_grid, aging, gamma)]


def run_aging_fit(
    spec,
    nu: float = 8.0,
    *,
    stage: str = "aging_S4gv",
    n_knots: int = 4,
    gamma_form: str = "varying",
    use_d: bool = False,
    solvers: tuple[str, ...] = ("anderson",),
    n_random: int = 10,
    jit_u: float = 0.10,
    jit_v: float = 0.10,
    seed0: int = 0,
    max_iter: int = 2000,
    tol: float = 1e-12,
    v_thin: int = 1,
    keep_v_traces: bool = True,
    register: bool = True,
    study: str = "aging",
) -> dict:
    """Fit the absolute-best production aging model and record convergence.

    Warm-starts from the slice's registered baseline nu fit (anchor + `n_random`
    perturbed restarts around it), fits under each solver in `solvers` at
    tolerance `tol`, keeps the single overall-best (max loglik) as the
    registered model, and writes the per-iter v_j + aging-coef traces.

    Raises SystemExit if the baseline warm source is missing (warm start from
    the baseline is a hard requirement of this stage).
    """
    cfg = aging_cfg(nu, n_knots=n_knots, gamma_form=gamma_form, use_d=use_d)
    sl = registry.slice_slug(spec)
    print(f"\n=== {sl}  nu={nu}  spline{n_knots}-g{gamma_form}  (stage={stage}) ===",
          flush=True)
    fd = load_slice(spec)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    base_uv = load_baseline_warmstart(spec, fd, nu)
    if base_uv is None:
        raise SystemExit(
            f"ERROR: no baseline nu={nu:g} fit on disk for '{sl}'.\n"
            f"Run the baseline first:\n"
            f"    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py "
            f"--slices {sl.split('_14')[0]}")
    print(f"    warm source: baseline_{registry.nu_tag(nu)}_best  "
          f"(anchor 'warm' + rand0..{n_random - 1})", flush=True)

    # anchor = baseline (u, v); restarts perturb around it (aging coefs zero-fill).
    init_list = init_schemes.build_inits(
        fd, n_random=n_random, jit_u=jit_u, jit_v=jit_v, seed0=seed0, base=base_uv)

    # common A_n reconstruction grid for the final best-init curve.
    A_n = fd.A_n[np.isfinite(fd.A_n)]
    A_max = float(A_n.max()) if A_n.size else 1.0
    A_grid = np.linspace(0.0, A_max, N_GRID)

    spread_rows: list[dict] = []
    ll_rows: list[dict] = []
    v_rows: list[dict] = []
    coef_rows: list[dict] = []
    final_v: dict[tuple[str, str], np.ndarray] = {}
    best: tuple[float, str, str, int | None, object] | None = None  # (ll,solver,init,seed,model)

    for solver in solvers:
        ModelCls = _model_cls(solver)
        for init_name, seed, ws in init_list:
            model = ModelCls(fd, cfg, _fitter(solver, ws, max_iter=max_iter, tol=tol))
            (res, v_arr, theta_arr, gamma_arr, cap_its,
             ll, rss, it_idx, dt) = _fit_capture(model, v_thin=v_thin)
            final_v[(solver, init_name)] = model.params["v"].copy()

            spread_rows.append(dict(
                slice=sl, solver=solver, init=init_name, seed=seed,
                n_iter=res.n_iter, converged=res.converged,
                loglik_final=res.loglik_final, rss_final=res.rss_final,
                wall_s=dt, mean_v=float(model.params["v"].mean()),
            ))
            for t, l, r in zip(it_idx, ll, rss):
                ll_rows.append(dict(slice=sl, solver=solver, init=init_name,
                                    iter=int(t), loglik=float(l), rss=float(r)))
            # per-iter aging coefficients (theta_aging + gamma), long format.
            for t, th, gm in zip(cap_its, theta_arr, gamma_arr):
                for k, val in enumerate(th):
                    coef_rows.append(dict(slice=sl, solver=solver, init=init_name,
                                          iter=int(t), kind="theta", k=k, value=float(val)))
                for k, val in enumerate(gm):
                    coef_rows.append(dict(slice=sl, solver=solver, init=init_name,
                                          iter=int(t), kind="gamma", k=k, value=float(val)))
            if keep_v_traces:
                for t, vvec in zip(cap_its, v_arr):
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
    spread_df["loglik_gap_to_best"] = spread_df["loglik_final"].max() - spread_df["loglik_final"]

    # converged-v agreement across inits, within each solver (pairwise).
    agree_rows: list[dict] = []
    for solver in solvers:
        keys = [k for k in final_v if k[0] == solver]
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                va, vb = final_v[keys[a]], final_v[keys[b]]
                agree_rows.append(dict(
                    slice=sl, solver=solver, init_a=keys[a][1], init_b=keys[b][1],
                    corr_v=_safe_corr(va, vb), max_abs_dv=float(np.max(np.abs(va - vb)))))

    # cross-solver agreement (only when >1 solver was run).
    best_v_by_solver: dict[str, np.ndarray] = {}
    for solver in solvers:
        sub = spread_df[spread_df.solver == solver]
        bi = sub.loc[sub.loglik_final.idxmax(), "init"]
        best_v_by_solver[solver] = final_v[(solver, bi)]
    if len(solvers) > 1:
        s0, s1 = solvers[0], solvers[1]
        cross_corr = _safe_corr(best_v_by_solver[s0], best_v_by_solver[s1])
        cross_dv = float(np.max(np.abs(best_v_by_solver[s0] - best_v_by_solver[s1])))
    else:
        cross_corr, cross_dv = float("nan"), float("nan")

    # per-solver convergence rollup.
    conv_rows: list[dict] = []
    for solver in solvers:
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

    ll_best, s_best, i_best, seed_best, model = best

    # ── write convergence diagnostics ────────────────────────────────
    conv_dir = CONV_ROOT / sl / stage
    conv_dir.mkdir(parents=True, exist_ok=True)
    spread_df.to_csv(conv_dir / "init_spread.csv", index=False)
    pd.DataFrame(ll_rows).to_csv(conv_dir / "loglik_traces.csv", index=False)
    conv_df.to_csv(conv_dir / "convergence_summary.csv", index=False)
    if agree_rows:
        pd.DataFrame(agree_rows).to_csv(conv_dir / "v_agreement.csv", index=False)
    if coef_rows:
        pd.DataFrame(coef_rows).to_parquet(conv_dir / "aging_coef_traces.parquet", index=False)
    if keep_v_traces and v_rows:
        pd.DataFrame(v_rows).to_parquet(conv_dir / "v_traces.parquet", index=False)
    pd.DataFrame(_curve_rows(model, A_grid=A_grid, ae_off=AE_OFFSET_YR,
                             slug=sl, nu=nu)).to_parquet(conv_dir / "curves.parquet",
                                                         index=False)

    for _, r in conv_df.iterrows():
        print(f"    [{r.solver}] n_iter(best/med/max)={r.n_iter_best}/{r.n_iter_median:g}/"
              f"{r.n_iter_max}  loglik gap across inits={r.loglik_gap_across_inits:.2e}",
              flush=True)
    print(f"    overall best: solver={s_best} init={i_best}  loglik={ll_best:.6f}", flush=True)
    if len(solvers) > 1:
        print(f"    cross-solver v: corr={cross_corr:.6f}  max|dv|={cross_dv:.2e}", flush=True)
    print(f"    convergence -> {conv_dir}", flush=True)

    # ── register the single overall-best model ───────────────────────
    fit_dir = None
    if register:
        parent = MODELS_ROOT / sl
        stem = aging_stem(cfg)                       # agingS4gv_nu8p00_best
        fit_dir = registry.fit_path(parent, stem, spec, cfg, resample_tag="base")
        registry.register_fit(
            model, fit_dir, resample_tag="base", save=SaveSpec(params=True),
            study=study, best_solver=s_best, best_init=i_best,
            best_init_seed=(int(seed_best) if seed_best is not None else None),
            warm_source=f"baseline_{registry.nu_tag(nu)}_best",
            n_inits=len(init_list), tol=tol, stage=stage,
            cross_solver_corr_v=cross_corr, cross_solver_max_dv=cross_dv,
        )
        print(f"    model -> {fit_dir}", flush=True)

    return dict(model=model, best_solver=s_best, best_init=i_best,
                loglik_best=ll_best, fit_dir=fit_dir, conv_dir=conv_dir,
                spread_df=spread_df, conv_df=conv_df)
