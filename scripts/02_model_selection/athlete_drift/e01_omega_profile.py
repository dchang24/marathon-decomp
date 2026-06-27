"""Frozen-omega_d2 sweep: marginal-likelihood profile + result insensitivity.

The per-athlete drift prior is d_i ~ N(0, omega_d2). In production omega_d2 is
learned by the EB / type-II MLE fixed point ``omega <- mean(d_hat^2 + post_var)``.
This script asks two questions by SWEEPING omega_d2 at frozen values around the
learned optimum omega*:

  (A) IDENTIFICATION -- is omega_d2 identified (a peaked marginal), and does the
      EB-learned omega* sit at the peak? The *penalized* loglik the fitter tracks
      is NOT comparable across omega (it drops the -n_elig/2*log(omega) prior
      normalizer), so it rises monotonically in omega. The right object is the
      type-II MARGINAL likelihood, which integrates d_i out and self-regularizes.
      We compute its Laplace approximation (exact under Gaussian noise; the same
      IRLS linearization the fitter already uses under Student-t):

          logML(omega) = data_ll(d_hat)                 # data fit at the mode
                       - (1/(2 omega)) * sum d_hat^2     # ridge  (favors large omega)
                       - (n_elig/2) * log(omega)         # prior normalizer (favors small)
                       - (1/2) * sum log(s_den_i + 1/omega)

      The last two terms oppose the first -> an interior peak. Using the live
      model, ``post_var_i = 1/(s_den_i + 1/omega)``, so the det term is just
      ``+(1/2) sum log(post_var_i)`` -- no s_den reconstruction needed.

  (B) INSENSITIVITY -- do the *reported* estimands move with omega? We track v_j
      (corr + max|dv| to the free fit), the reconstructed aging & entry-age
      curves (max/RMS deviation to free), u_i, and the d-distribution summaries.
      Flat near omega* => freeze-vs-free and the exact omega value are immaterial.

Also recorded per omega for context: data/penalized loglik, naive & effective
AIC/BIC (naive uses the fixed param count -> monotone; effective uses the
shrinkage edf), edf_d, runtime, convergence. The naive ICs improving while logML
peaks is the headline contrast.

Each grid fit is a full AxD refit (aging spline-4 / varying-gamma / nu=8, the
production form), warm-started from the slice's registered ``agingS4gv`` no-d fit
(free fit) then from the free fit (frozen fits). These are THROWAWAY fits -- not
registered. Outputs go to
``results/model_selection/athlete_drift/omega_profile/{slug}/``.

Run::

    python scripts/02_model_selection/athlete_drift/e01_omega_profile.py            # Po10_M
    python scripts/02_model_selection/athlete_drift/e01_omega_profile.py --slices Po10_M Po10_W
    python scripts/02_model_selection/athlete_drift/e01_omega_profile.py --slices Po10_M --tol 1e-12
"""
from __future__ import annotations

import argparse
import dataclasses
import glob
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import (  # noqa: E402
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    ModelConfig,
    aging_curve_on_grid,
    entry_age_curve_on_grid,
    load_slice,
    registry,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

from baseline_common.slices import build_spec, resolve_names  # noqa: E402

MODELS_ROOT = RESULTS_DIR / "models"
OUT_ROOT = RESULTS_DIR / "model_selection" / "athlete_drift" / "omega_profile"

AE_OFFSET_YR = 10.0          # entry-age offset for the stored gamma curve (matches aging stage)
N_GRID = 200                 # A_n reconstruction grid
# Frozen omega_d2 = omega* * MULT. The 1.0 point is the freeze-machinery sanity
# gate (must reproduce the free fit); the +-100x extremes show the falloff.
MULTS = (1 / 100, 1 / 30, 1 / 10, 1 / 3, 1 / 2, 1, 2, 3, 10, 30, 100)


def _model_cls(solver: str):
    return ModelAnderson if solver == "anderson" else Model


def _fitter(solver: str, ws: dict, *, max_iter: int, tol: float):
    common = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                  init="warmstart", warmstart=ws, record_trace=False, verbose=0)
    return AndersonFitterConfig(**common) if solver == "anderson" else FitterConfig(**common)


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def load_aging_warmstart(slug: str, fd, nu: float) -> dict | None:
    """The slice's registered ``agingS4gv_{nutag}_best`` no-d fit as a warmstart."""
    stem = f"agingS4gv_{registry.nu_tag(float(nu))}_best"
    hits = sorted(glob.glob(str(MODELS_ROOT / slug / f"{stem}__*")))
    if not hits:
        return None
    m = registry.load_fit(Path(hits[0]), fd)
    return {k: np.asarray(m.params[k]).copy()
            for k in ("u", "v", "theta_aging", "gamma")}


def _metrics(model, res, fd, nu, *, A_grid, free_ref, dt) -> tuple[dict, tuple]:
    """One profile row + the (u, v, aging, gamma) curves for the free-ref / fan."""
    p = model.params
    r = model.residuals()
    N = fd.N
    sigma2 = float(p["sigma2"])
    omega = float(p["omega_d2"])
    elig = np.asarray(model.eligible_d, dtype=bool)
    n_elig = int(elig.sum())
    d = np.asarray(p["d"], float)
    u = np.asarray(p["u"], float)
    v = np.asarray(p["v"], float)

    # --- likelihoods -------------------------------------------------
    data_ll = float(Model._loglik_from_resid(r, sigma2, float(nu), N))
    eb_pen = float(model._eb_penalty(d, omega))     # = -1/(2 omega) sum d^2
    pen_ll = data_ll + eb_pen

    # --- Laplace marginal: det term = +1/2 sum log(post_var) ---------
    pv = np.maximum(np.asarray(model._post_var["d"], float)[elig], 1e-300)
    logML = (data_ll + eb_pen
             - 0.5 * n_elig * np.log(max(omega, 1e-300))
             + 0.5 * float(np.sum(np.log(pv))))

    # --- complexity / ICs --------------------------------------------
    edf = model.effective_dof()
    edf_d = float(edf.get("d", 0.0))
    free_omega = not (model.config.omega_d2_fixed is not None
                      or model.config.freeze_eb_prior)
    edf_total = float(sum(edf.values())) + 1.0 + (1.0 if free_omega else 0.0)  # +sigma2 (+omega)
    k_naive = int(model.n_params_naive())
    logN = float(np.log(N))
    aic_naive = 2 * k_naive - 2 * data_ll
    bic_naive = k_naive * logN - 2 * data_ll
    aic_eff = 2 * edf_total - 2 * data_ll
    bic_eff = edf_total * logN - 2 * data_ll

    # --- d distribution ----------------------------------------------
    de = d[elig]
    d_sd = float(de.std(ddof=1)) if n_elig > 1 else float("nan")
    frac_improver = float((de < 0).mean()) if n_elig else float("nan")
    z = de / np.sqrt(pv)
    frac_credible = float((np.abs(z) > 2).mean()) if n_elig else float("nan")
    corr_du = _safe_corr(de, u[elig])

    # --- curves ------------------------------------------------------
    aging = np.asarray(aging_curve_on_grid(model, A_grid), float)
    gamma_c = np.asarray(entry_age_curve_on_grid(model, A_grid, AE_OFFSET_YR), float)

    row = dict(
        omega_d2=omega, is_free=free_omega, n_elig=n_elig,
        data_loglik=data_ll, penalized_loglik=pen_ll, logML=logML,
        edf_d=edf_d, edf_total=edf_total, k_naive=k_naive,
        aic_naive=aic_naive, bic_naive=bic_naive,
        aic_eff=aic_eff, bic_eff=bic_eff,
        d_sd=d_sd, frac_improver=frac_improver, frac_credible=frac_credible,
        corr_du=corr_du, sigma2=sigma2,
        n_iter=int(res.n_iter), converged=bool(res.converged), wall_s=float(dt),
    )

    if free_ref is None:                         # this IS the free row
        row.update(corr_v_to_free=1.0, max_abs_dv=0.0, rms_dv=0.0,
                   corr_u_to_free=1.0, aging_maxdev=0.0, aging_rmsdev=0.0,
                   gamma_maxdev=0.0)
    else:
        uf, vf, agf, gcf = free_ref
        row.update(
            corr_v_to_free=_safe_corr(v, vf),
            max_abs_dv=float(np.max(np.abs(v - vf))),
            rms_dv=float(np.sqrt(np.mean((v - vf) ** 2))),
            corr_u_to_free=_safe_corr(u, uf),
            aging_maxdev=float(np.max(np.abs(aging - agf))),
            aging_rmsdev=float(np.sqrt(np.mean((aging - agf) ** 2))),
            gamma_maxdev=float(np.max(np.abs(gamma_c - gcf))),
        )
    return row, (u, v, aging, gamma_c)


def run_slice(name: str, *, nu: float, solver: str, tol: float,
              max_iter: int) -> None:
    spec = build_spec(name)
    slug = registry.slice_slug(spec)
    nutag = registry.nu_tag(float(nu))
    fd = load_slice(spec)
    print(f"\n=== {slug}  nu={nu:g}  solver={solver}  tol={tol:g} ===", flush=True)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    ws_aging = load_aging_warmstart(slug, fd, nu)
    if ws_aging is None:
        print(f"    [skip] no agingS4gv_{nutag}_best fit under {slug} "
              f"(run the aging production fit first)")
        return

    A_n = fd.A_n[np.isfinite(fd.A_n)]
    A_grid = np.linspace(0.0, float(A_n.max()) if A_n.size else 1.0, N_GRID)

    cfg_free = ModelConfig(use_phi12=True, use_gamma=True, use_d=True,
                           basis_kind="spline", n_knots=4, gamma_form="varying",
                           nu=float(nu))
    ModelCls = _model_cls(solver)

    # ---- free EB fit -> omega* + reference curves -------------------
    t0 = time.perf_counter()
    m_free = ModelCls(fd, cfg_free, _fitter(solver, ws_aging, max_iter=max_iter, tol=tol))
    res = m_free.fit()
    dt = time.perf_counter() - t0
    row, free_ref = _metrics(m_free, res, fd, nu, A_grid=A_grid, free_ref=None, dt=dt)
    omega_star = row["omega_d2"]
    row["omega_mult"] = 1.0
    rows = [row]
    curve_rows = [dict(slug=slug, nu=float(nu), omega_d2=omega_star, is_free=True,
                       omega_mult=1.0, A_n=float(a), aging=float(ag), gamma=float(gc))
                  for a, ag, gc in zip(A_grid, free_ref[2], free_ref[3])]
    print(f"    [free ] omega*={omega_star:.4e}  logML={row['logML']:.3f}  "
          f"edf_d={row['edf_d']:.1f}  iters={res.n_iter}  {dt:.1f}s", flush=True)

    # warm-start the frozen fits from the free fit (same basin, faster).
    ws_free = {k: np.asarray(m_free.params[k]).copy()
               for k in ("u", "v", "d", "theta_aging", "gamma")}

    # ---- frozen grid ------------------------------------------------
    for mult in MULTS:
        om = omega_star * mult
        cfg = dataclasses.replace(cfg_free, omega_d2_fixed=float(om))
        t0 = time.perf_counter()
        m = ModelCls(fd, cfg, _fitter(solver, ws_free, max_iter=max_iter, tol=tol))
        res = m.fit()
        dt = time.perf_counter() - t0
        row, _ = _metrics(m, res, fd, nu, A_grid=A_grid, free_ref=free_ref, dt=dt)
        row["omega_mult"] = float(mult)
        rows.append(row)
        aging_c = aging_curve_on_grid(m, A_grid)
        gamma_c = entry_age_curve_on_grid(m, A_grid, AE_OFFSET_YR)
        curve_rows.extend(
            dict(slug=slug, nu=float(nu), omega_d2=float(om), is_free=False,
                 omega_mult=float(mult), A_n=float(a), aging=float(ag), gamma=float(gc))
            for a, ag, gc in zip(A_grid, aging_c, gamma_c))
        print(f"    [x{mult:<5g}] omega={om:.4e}  logML={row['logML']:.3f}  "
              f"dataLL={row['data_loglik']:.1f}  edf_d={row['edf_d']:.1f}  "
              f"aging_maxdev={row['aging_maxdev']:.2e}  max|dv|={row['max_abs_dv']:.2e}  "
              f"iters={res.n_iter}  {dt:.1f}s", flush=True)

    # ---- write ------------------------------------------------------
    out = OUT_ROOT / slug
    out.mkdir(parents=True, exist_ok=True)
    prof = pd.DataFrame(rows).sort_values("omega_d2").reset_index(drop=True)
    prof.insert(0, "slug", slug)
    prof.insert(1, "nu", float(nu))
    prof.to_parquet(out / f"profile_{nutag}.parquet", index=False)
    pd.DataFrame(curve_rows).to_parquet(out / f"curves_{nutag}.parquet", index=False)

    # peak vs free
    peak = prof.loc[prof["logML"].idxmax()]
    print(f"    logML peak at omega={peak['omega_d2']:.4e} "
          f"(x{peak['omega_mult']:g}); free omega*={omega_star:.4e}  "
          f"dlogML(peak-free)={peak['logML'] - rows[0]['logML']:.3f}", flush=True)
    print(f"    wrote {out / f'profile_{nutag}.parquet'}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["Po10_M"],
                    help="named slices or 'all' (default: Po10_M).")
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--solver", choices=["anderson", "als"], default="anderson")
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--max-iter", type=int, default=2000)
    args = ap.parse_args()

    names = resolve_names(args.slices, ap)
    for name in names:
        run_slice(name, nu=args.nu, solver=args.solver, tol=args.tol,
                  max_iter=args.max_iter)


if __name__ == "__main__":
    main()
