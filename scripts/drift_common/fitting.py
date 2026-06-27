"""Production fit for the per-athlete drift models (`full` and `drift`).

Two variants at the settled operating point (Student-t nu=8, Anderson):

  * **full**  = baseline + aging (spline-4, varying gamma) + d_i
                -> warm-started from the slice's registered ``agingS4gv_nu8p00_best``
                   (u, v, theta_aging, gamma seeded; d zero-filled).
                -> stem ``full_nu8p00_best``.
  * **drift** = baseline + d_i  (no aging block)
                -> warm-started from the slice's registered ``baseline_nu8p00_best``
                   (u, v seeded; d zero-filled).
                -> stem ``drift_nu8p00_best``.

Mirrors ``aging_common.fitting.run_aging_fit``: a deterministic warm anchor +
``n_random`` perturbed restarts (probing init-invariance of the fixed point),
fit under each solver, keep the single overall-best (max loglik) as the
registered model, and write per-iter convergence traces. The EB drift prior
variance ``omega_d2`` is learned by default (type-II MLE) from
``omega_d2_init`` (model default 1e-4); ``freeze_eb`` locks it at the init.

Slice/eligibility/prior knobs (``omega_d2_init``, ``freeze_eb``, ``d_min_span``)
flow into the model config and therefore the **identity hash**, so exploratory
variants never collide with the production point.

Outputs (per slice ``{slug}``):
  * model       -> results/models/{slug}/{full|drift}_nu8p00_best__<hash>/
  * convergence -> results/convergence/{slug}/{stage}/  (init_spread, loglik+omega
                   traces, v_traces, convergence_summary, v_agreement; aging
                   coef_traces + curves for the `full` variant only)
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from marathon_decomp import (
    ModelConfig,
    SaveSpec,
    aging_curve_on_grid,
    entry_age_curve_on_grid,
    load_slice,
    registry,
)
from marathon_decomp.config import RESULTS_DIR

from baseline_common import inits as init_schemes
# low-level helpers shared with the aging stage
from aging_common.fitting import (
    AE_OFFSET_YR,
    N_GRID,
    _fitter,
    _model_cls,
    _safe_corr,
    aging_cfg,
    aging_stem,
    load_baseline_warmstart,
)

MODELS_ROOT = RESULTS_DIR / "models"
CONV_ROOT = RESULTS_DIR / "convergence"

VARIANTS = ("full", "drift")


def drift_cfg(nu: float, *, variant: str, omega_d2_init: float | None = None,
              freeze_eb: bool = False, d_min_span: float = 1e-3,
              n_knots: int = 4, gamma_form: str = "varying") -> ModelConfig:
    """Production drift config. ``full`` keeps the aging block; ``drift`` drops it.

    ``omega_d2_init=None`` -> the model default (1e-4); ``freeze_eb`` locks
    ``omega_d2`` at the init instead of learning it. ``d_min_span`` is the
    career-span eligibility floor in years (production: 1e-3, an effective-zero
    safeguard -- see ``02_model_selection/athlete_drift``).
    """
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")
    aging_on = variant == "full"
    return ModelConfig(
        use_phi12=aging_on, use_gamma=aging_on, use_d=True,
        basis_kind="spline", n_knots=n_knots, gamma_form=gamma_form,
        nu=float(nu),
        omega_d2_init=(None if omega_d2_init is None else float(omega_d2_init)),
        freeze_eb_prior=bool(freeze_eb),
        d_min_span_years=float(d_min_span),
    )


def drift_stem(variant: str, nu: float) -> str:
    """``full_nu8p00_best`` / ``drift_nu8p00_best``."""
    return f"{variant}_{registry.nu_tag(float(nu))}_best"


def load_aging_warmstart(spec, fd, nu: float, *, n_knots: int = 4,
                         gamma_form: str = "varying") -> dict | None:
    """The slice's registered ``agingS4gv_{nutag}_best`` no-d fit (u,v,theta,gamma).

    Reconstructed from the deterministic registry path (same cfg the aging stage
    registered under). Returns None if absent.
    """
    cfg_ag = aging_cfg(nu, n_knots=n_knots, gamma_form=gamma_form, use_d=False)
    parent = MODELS_ROOT / registry.slice_slug(spec)
    fit_dir = registry.fit_path(parent, aging_stem(cfg_ag), spec, cfg_ag,
                                resample_tag="base")
    if not (fit_dir / "fit.pkl").is_file():
        return None
    m = registry.load_fit(fit_dir, fd)
    return {k: np.asarray(m.params[k]).copy()
            for k in ("u", "v", "theta_aging", "gamma")}


def _fit_capture(model, *, v_thin: int):
    """Fit `model`, capturing per-iter loglik/rss/omega_d2 (all) + v (thinned).

    Returns (res, v_arr[(T,J)], cap_its, ll[], rss[], om[], it_idx[], wall_s).
    """
    v_iters: list[np.ndarray] = []
    cap_its: list[int] = []
    ll: list[float] = []
    rss: list[float] = []
    om: list[float] = []
    it_idx: list[int] = []

    def hook(m, it, loglik, rss_):
        ll.append(float(loglik)); rss.append(float(rss_))
        om.append(float(m.params["omega_d2"])); it_idx.append(int(it))
        if it % v_thin == 0:
            v_iters.append(m.params["v"].copy())
            cap_its.append(int(it))

    model.iter_hook = hook
    t0 = time.perf_counter()
    res = model.fit()
    dt = time.perf_counter() - t0
    v_arr = np.array(v_iters) if v_iters else np.empty((0, model.data.J))
    return (res, v_arr, cap_its, np.array(ll), np.array(rss), np.array(om),
            np.array(it_idx), dt)


def _d_summary(model) -> dict:
    """Best-fit drift summary for the manifest (eligible athletes only)."""
    elig = np.asarray(model.eligible_d, dtype=bool)
    d = np.asarray(model.params["d"], dtype=float)[elig]
    n = int(elig.sum())
    return dict(
        n_eligible_d=n,
        omega_d2=float(model.params["omega_d2"]),
        d_sd=float(d.std(ddof=1)) if n > 1 else float("nan"),
        frac_improver=float((d < 0).mean()) if n else float("nan"),
    )


def run_drift_fit(
    spec,
    variant: str = "full",
    nu: float = 8.0,
    *,
    omega_d2_init: float | None = None,
    freeze_eb: bool = False,
    d_min_span: float = 1e-3,
    n_knots: int = 4,
    gamma_form: str = "varying",
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
    study: str = "drift",
) -> dict:
    """Fit the absolute-best production drift model and record convergence.

    Warm-starts from the variant's source fit (``agingS4gv`` for ``full``,
    ``baseline`` for ``drift``), runs an anchor + `n_random` perturbed restarts
    under each solver, keeps the overall-best (max loglik) as the registered
    model, and writes per-iter v_j / loglik / omega_d2 traces.

    Raises SystemExit if the warm source is missing.
    """
    cfg = drift_cfg(nu, variant=variant, omega_d2_init=omega_d2_init,
                    freeze_eb=freeze_eb, d_min_span=d_min_span,
                    n_knots=n_knots, gamma_form=gamma_form)
    sl = registry.slice_slug(spec)
    stage = f"drift_{variant}"
    nutag = registry.nu_tag(float(nu))
    print(f"\n=== {sl}  nu={nu:g}  variant={variant}  "
          f"omega_init={'default' if omega_d2_init is None else omega_d2_init}  "
          f"freeze_eb={freeze_eb}  min_span={d_min_span:g}  (stage={stage}) ===",
          flush=True)
    fd = load_slice(spec)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    # --- warm source ------------------------------------------------
    if variant == "full":
        ws_src = load_aging_warmstart(spec, fd, nu, n_knots=n_knots, gamma_form=gamma_form)
        warm_name = f"agingS4gv_{nutag}_best"
        if ws_src is None:
            raise SystemExit(
                f"ERROR: no agingS4gv_{nutag}_best fit for '{sl}'. Run first:\n"
                f"    python scripts/03_model_fit/aging/e01_fit_aging.py --slices "
                f"{sl.split('_14')[0]}")
        base_uv = {"u": ws_src["u"], "v": ws_src["v"]}
        aging_seed = {"theta_aging": ws_src["theta_aging"], "gamma": ws_src["gamma"]}
    else:  # drift
        base_uv = load_baseline_warmstart(spec, fd, nu)
        warm_name = f"baseline_{nutag}_best"
        if base_uv is None:
            raise SystemExit(
                f"ERROR: no baseline_{nutag}_best fit for '{sl}'. Run first:\n"
                f"    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py "
                f"--slices {sl.split('_14')[0]}")
        aging_seed = {}
    print(f"    warm source: {warm_name}  (anchor 'warm' + rand0..{n_random - 1})",
          flush=True)

    # anchor = warm (u,v); restarts perturb (u,v). For `full`, every init keeps
    # the warm aging coefficients (un-perturbed); d zero-fills in the fitter.
    init_list = init_schemes.build_inits(
        fd, n_random=n_random, jit_u=jit_u, jit_v=jit_v, seed0=seed0, base=base_uv)
    for _, _, ws in init_list:
        ws.update(aging_seed)

    A_n = fd.A_n[np.isfinite(fd.A_n)]
    A_grid = np.linspace(0.0, float(A_n.max()) if A_n.size else 1.0, N_GRID)

    spread_rows: list[dict] = []
    ll_rows: list[dict] = []
    v_rows: list[dict] = []
    final_v: dict[tuple[str, str], np.ndarray] = {}
    best: tuple[float, str, str, int | None, object] | None = None

    for solver in solvers:
        ModelCls = _model_cls(solver)
        for init_name, seed, ws in init_list:
            model = ModelCls(fd, cfg, _fitter(solver, ws, max_iter=max_iter, tol=tol))
            res, v_arr, cap_its, ll, rss, om, it_idx, dt = _fit_capture(model, v_thin=v_thin)
            final_v[(solver, init_name)] = model.params["v"].copy()
            ds = _d_summary(model)

            spread_rows.append(dict(
                slice=sl, variant=variant, solver=solver, init=init_name, seed=seed,
                n_iter=res.n_iter, converged=res.converged,
                loglik_final=res.loglik_final, rss_final=res.rss_final,
                omega_d2=ds["omega_d2"], n_eligible_d=ds["n_eligible_d"],
                d_sd=ds["d_sd"], frac_improver=ds["frac_improver"], wall_s=dt,
            ))
            for t, l, r, o in zip(it_idx, ll, rss, om):
                ll_rows.append(dict(slice=sl, solver=solver, init=init_name,
                                    iter=int(t), loglik=float(l), rss=float(r),
                                    omega_d2=float(o)))
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
                  f"  loglik={res.loglik_final:.6f}  omega_d2={ds['omega_d2']:.4e}"
                  f"  elig_d={ds['n_eligible_d']}  {dt:.1f}s", flush=True)

    spread_df = pd.DataFrame(spread_rows)
    spread_df["loglik_gap_to_best"] = spread_df["loglik_final"].max() - spread_df["loglik_final"]

    # converged-v agreement across inits, per solver
    agree_rows: list[dict] = []
    for solver in solvers:
        keys = [k for k in final_v if k[0] == solver]
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                va, vb = final_v[keys[a]], final_v[keys[b]]
                agree_rows.append(dict(
                    slice=sl, solver=solver, init_a=keys[a][1], init_b=keys[b][1],
                    corr_v=_safe_corr(va, vb), max_abs_dv=float(np.max(np.abs(va - vb)))))

    best_v_by_solver: dict[str, np.ndarray] = {}
    for solver in solvers:
        sub = spread_df[spread_df.solver == solver]
        best_v_by_solver[solver] = final_v[(solver, sub.loc[sub.loglik_final.idxmax(), "init"])]
    if len(solvers) > 1:
        s0, s1 = solvers[0], solvers[1]
        cross_corr = _safe_corr(best_v_by_solver[s0], best_v_by_solver[s1])
        cross_dv = float(np.max(np.abs(best_v_by_solver[s0] - best_v_by_solver[s1])))
    else:
        cross_corr, cross_dv = float("nan"), float("nan")

    conv_rows: list[dict] = []
    for solver in solvers:
        sub = spread_df[spread_df.solver == solver]
        conv_rows.append(dict(
            slice=sl, variant=variant, solver=solver, nu=float(nu),
            loglik_best=float(sub.loglik_final.max()),
            n_iter_best=int(sub.loc[sub.loglik_final.idxmax(), "n_iter"]),
            n_iter_median=float(sub["n_iter"].median()),
            n_iter_max=int(sub["n_iter"].max()),
            wall_s_best=float(sub.loc[sub.loglik_final.idxmax(), "wall_s"]),
            all_converged=bool(sub["converged"].all()),
            loglik_gap_across_inits=float(sub["loglik_final"].max() - sub["loglik_final"].min()),
            omega_d2_spread=float(sub["omega_d2"].max() - sub["omega_d2"].min()),
            max_abs_dv_across_inits=float(
                max((r["max_abs_dv"] for r in agree_rows if r["solver"] == solver), default=0.0)),
            cross_solver_corr_v=cross_corr, cross_solver_max_dv=cross_dv,
        ))
    conv_df = pd.DataFrame(conv_rows)

    ll_best, s_best, i_best, seed_best, model = best

    # --- write convergence diagnostics ------------------------------
    conv_dir = CONV_ROOT / sl / stage
    conv_dir.mkdir(parents=True, exist_ok=True)
    spread_df.to_csv(conv_dir / "init_spread.csv", index=False)
    pd.DataFrame(ll_rows).to_csv(conv_dir / "loglik_traces.csv", index=False)
    conv_df.to_csv(conv_dir / "convergence_summary.csv", index=False)
    if agree_rows:
        pd.DataFrame(agree_rows).to_csv(conv_dir / "v_agreement.csv", index=False)
    if keep_v_traces and v_rows:
        pd.DataFrame(v_rows).to_parquet(conv_dir / "v_traces.parquet", index=False)
    if variant == "full":
        aging = aging_curve_on_grid(model, A_grid)
        gamma = entry_age_curve_on_grid(model, A_grid, AE_OFFSET_YR)
        pd.DataFrame([dict(slug=sl, nu=float(nu), A_n=float(a),
                           aging_curve=float(yc), gamma_curve=float(gc),
                           ae_offset_yr=AE_OFFSET_YR)
                      for a, yc, gc in zip(A_grid, aging, gamma)]
                     ).to_parquet(conv_dir / "curves.parquet", index=False)

    for _, r in conv_df.iterrows():
        print(f"    [{r.solver}] n_iter(best/med/max)={r.n_iter_best}/{r.n_iter_median:g}/"
              f"{r.n_iter_max}  loglik gap across inits={r.loglik_gap_across_inits:.2e}  "
              f"omega_d2 spread={r.omega_d2_spread:.2e}", flush=True)
    print(f"    overall best: solver={s_best} init={i_best}  loglik={ll_best:.6f}", flush=True)
    print(f"    convergence -> {conv_dir}", flush=True)

    # --- register the overall-best model ----------------------------
    fit_dir = None
    if register:
        parent = MODELS_ROOT / sl
        stem = drift_stem(variant, nu)
        fit_dir = registry.fit_path(parent, stem, spec, cfg, resample_tag="base")
        ds = _d_summary(model)
        registry.register_fit(
            model, fit_dir, resample_tag="base", save=SaveSpec(params=True),
            study=study, variant=variant, best_solver=s_best, best_init=i_best,
            best_init_seed=(int(seed_best) if seed_best is not None else None),
            warm_source=warm_name, n_inits=len(init_list), tol=tol, stage=stage,
            omega_d2=ds["omega_d2"], omega_d2_init=omega_d2_init,
            freeze_eb_prior=freeze_eb, d_min_span_years=d_min_span,
            n_eligible_d=ds["n_eligible_d"],
            cross_solver_corr_v=cross_corr, cross_solver_max_dv=cross_dv,
        )
        print(f"    model -> {fit_dir}", flush=True)

    return dict(model=model, variant=variant, best_solver=s_best, best_init=i_best,
                loglik_best=ll_best, fit_dir=fit_dir, conv_dir=conv_dir,
                spread_df=spread_df, conv_df=conv_df)
