"""omega_d2 INIT sensitivity: does the starting prior variance change anything?

Companion to ``e01_omega_profile.py``. There the prior variance was *frozen*;
here it is always **EB-learned** (free), and we vary only its **initial guess**
``omega_d2_init`` (the iter-0 value; production default is ``sigma2/1e4``). The
EB update ``omega <- mean(d_hat^2 + post_var)`` is an EM contraction, so the
converged ``omega*`` and every estimand should be **init-invariant** — the init
should affect only the iteration count (a bad start costs convergence time). This
script confirms that and quantifies the cost.

To isolate the omega-init effect, the ``(u, v, aging, gamma)`` blocks are
warm-started **identically** from the registered ``agingS4gv`` no-d fit for every
run; only ``omega_d2_init`` changes. A huge init starts the drift block
under-shrunk (large iter-1 d_i that can perturb u/v/aging); a tiny init starts
over-shrunk (d ~ 0). Either should recover to the same fixed point.

For each init we record the converged ``omega_d2``, convergence (n_iter, wall_s),
the likelihoods, and the deviation of every reported estimand (v_j, u_i, aging &
gamma curves, d) from the **reference** run (the production default init), plus
the full per-iter ``omega_d2`` and loglik **trajectory** (the headline: every
trajectory lands on the same omega*).

Outputs -> ``results/model_selection/athlete_drift/omega_init/{slug}/``
  * ``init_summary_{nutag}.parquet`` — one row per init.
  * ``init_traces_{nutag}.parquet``  — per (init, iter) omega_d2 + loglik.
Throwaway fits — not registered.

Run::

    python scripts/02_model_selection/athlete_drift/e02_omega_init.py            # Po10_M
    python scripts/02_model_selection/athlete_drift/e02_omega_init.py --slices Po10_M Po10_W
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/
sys.path.insert(0, str(Path(__file__).resolve().parent))      # this dir (for e01 reuse)

from marathon_decomp import ModelConfig, load_slice, registry  # noqa: E402
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

from baseline_common.slices import build_spec, resolve_names  # noqa: E402

# reuse the profile machinery (fitter builder, metric extractor, warmstart loader)
from e01_omega_profile import (  # noqa: E402
    N_GRID, _fitter, _metrics, _model_cls, load_aging_warmstart,
)

OUT_ROOT = RESULTS_DIR / "model_selection" / "athlete_drift" / "omega_init"

# Starting omega_d2 values spanning well below and above the typical omega*~3e-4.
# None = the production default (sigma2/1e4) and is used as the REFERENCE run.
INITS: tuple[float | None, ...] = (
    None, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
)


def _init_label(v: float | None) -> str:
    return "default" if v is None else f"{v:.0e}"


def _fit_with_trace(model):
    """Fit, capturing per-iter (iter, omega_d2, loglik)."""
    tr: list[tuple[int, float, float]] = []

    def hook(m, it, loglik, rss):
        tr.append((int(it), float(m.params["omega_d2"]), float(loglik)))

    model.iter_hook = hook
    t0 = time.perf_counter()
    res = model.fit()
    dt = time.perf_counter() - t0
    return res, dt, tr


def run_slice(name: str, *, nu: float, solver: str, tol: float, max_iter: int,
              inits=INITS) -> None:
    spec = build_spec(name)
    slug = registry.slice_slug(spec)
    nutag = registry.nu_tag(float(nu))
    fd = load_slice(spec)
    print(f"\n=== {slug}  nu={nu:g}  solver={solver}  tol={tol:g} ===", flush=True)
    print(f"    I={fd.I:,}  J={fd.J:,}  N={fd.N:,}", flush=True)

    ws_aging = load_aging_warmstart(slug, fd, nu)
    if ws_aging is None:
        print(f"    [skip] no agingS4gv_{nutag}_best fit under {slug}")
        return

    A_n = fd.A_n[np.isfinite(fd.A_n)]
    A_grid = np.linspace(0.0, float(A_n.max()) if A_n.size else 1.0, N_GRID)

    cfg0 = ModelConfig(use_phi12=True, use_gamma=True, use_d=True,
                       basis_kind="spline", n_knots=4, gamma_form="varying",
                       nu=float(nu))
    ModelCls = _model_cls(solver)

    rows: list[dict] = []
    trace_rows: list[dict] = []
    ref_curves = None
    ref_omega = None

    for init_val in inits:
        cfg = cfg0 if init_val is None else dataclasses.replace(
            cfg0, omega_d2_init=float(init_val))
        model = ModelCls(fd, cfg, _fitter(solver, ws_aging, max_iter=max_iter, tol=tol))
        res, dt, tr = _fit_with_trace(model)
        is_ref = ref_curves is None
        row, curves = _metrics(model, res, fd, nu, A_grid=A_grid,
                               free_ref=(None if is_ref else ref_curves), dt=dt)
        if is_ref:
            ref_curves = curves
            ref_omega = row["omega_d2"]
        row["omega_d2_init"] = (np.nan if init_val is None else float(init_val))
        row["init_label"] = _init_label(init_val)
        row["is_reference"] = is_ref
        row["omega_final_rel_to_ref"] = (
            float(row["omega_d2"] / ref_omega - 1.0) if ref_omega else np.nan)
        rows.append(row)
        for it, om, ll in tr:
            trace_rows.append(dict(slug=slug, nu=float(nu),
                                   init_label=row["init_label"], iter=it,
                                   omega_d2=om, loglik=ll))
        print(f"    init={row['init_label']:>14s}  ->  omega*={row['omega_d2']:.6e}  "
              f"(rel_ref={row['omega_final_rel_to_ref']:+.2e})  "
              f"logML={row['logML']:.3f}  aging_maxdev={row['aging_maxdev']:.2e}  "
              f"iters={res.n_iter}{'' if res.converged else ' MAX'}  {dt:.1f}s", flush=True)

    out = OUT_ROOT / slug
    out.mkdir(parents=True, exist_ok=True)
    summ = pd.DataFrame(rows)
    summ.insert(0, "slug", slug); summ.insert(1, "nu", float(nu))
    summ.to_parquet(out / f"init_summary_{nutag}.parquet", index=False)
    pd.DataFrame(trace_rows).to_parquet(out / f"init_traces_{nutag}.parquet", index=False)

    # headline spread across inits
    om = summ["omega_d2"].to_numpy()
    print(f"    converged omega* spread: min={om.min():.6e} max={om.max():.6e}  "
          f"(max rel-range {om.max()/om.min() - 1:.2e})", flush=True)
    print(f"    max aging dev across inits: {summ['aging_maxdev'].max():.2e}  "
          f"max |dv|: {summ['max_abs_dv'].max():.2e}  "
          f"iters range [{int(summ['n_iter'].min())}, {int(summ['n_iter'].max())}]",
          flush=True)
    print(f"    wrote {out / f'init_summary_{nutag}.parquet'}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["Po10_M"],
                    help="named slices or 'all' (default: Po10_M).")
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--solver", choices=["anderson", "als"], default="anderson")
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--inits", nargs="*", type=float, default=None,
                    help="omega_d2_init grid (floats); default uses the built-in "
                         "grid. The production default (sigma2/1e4) is always run "
                         "first as the reference. Use e.g. '--inits 1e-6 1e-5 1e-4 "
                         "1e-3' to skip the known 1e-8 stall on expensive slices.")
    args = ap.parse_args()

    inits = INITS if args.inits is None else (None, *(float(x) for x in args.inits))
    names = resolve_names(args.slices, ap)
    for name in names:
        run_slice(name, nu=args.nu, solver=args.solver, tol=args.tol,
                  max_iter=args.max_iter, inits=inits)


if __name__ == "__main__":
    main()
