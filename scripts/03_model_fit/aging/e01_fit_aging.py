"""Fit the production aging model (spline-4 + varying gamma) at nu=8.

The settled production aging form for **all runs from now on**: rank-1
``u_i + v_j`` plus a natural cubic spline aging block (4 knots) with a varying
entry-age gamma, Student-t nu=8, Anderson solver at a very tight tolerance
(1e-12). Per-athlete drift ``d_i`` stays off (a later production stage).

Each slice is **warm-started from its registered baseline nu=8 model**
(``baseline_nu8p00_best``; run ``03_model_fit/baseline/e02_fit_baseline_t.py``
first): the baseline ``(u, v)`` seed the anchor init and ``--n-random`` (default
10) perturbed restarts probe init-invariance. The aging coefficients zero-fill.

Writes (per slice):
  * essential model -> results/models/{slug}/agingS4gv_nu8p00_best__<hash>/
  * convergence     -> results/convergence/{slug}/aging_S4gv/
      per-iter v_j (v_traces.parquet) AND aging factors
      (aging_coef_traces.parquet: theta_aging + gamma per iter), plus the
      usual init-spread / loglik-trace / agreement tables and the best-init
      reconstructed aging+gamma curve (curves.parquet).

ZERO-ARG RUN (VS Code play): fits all eight slices, warm from baseline nu=8.

Run::

    python scripts/03_model_fit/aging/e01_fit_aging.py             # all 8
    python scripts/03_model_fit/aging/e01_fit_aging.py --slices Po10_M
    python scripts/03_model_fit/aging/e01_fit_aging.py --slices WA_M --min-runner 10
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from baseline_common import slices as S            # noqa: E402
from aging_common.fitting import run_aging_fit     # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--nu", type=float, default=8.0,
                    help="Student-t d.o.f. (production decision: nu=8).")
    ap.add_argument("--n-knots", type=int, default=4,
                    help="spline knots (production aging form: 4).")
    ap.add_argument("--gamma-form", default="varying", choices=["scalar", "varying"],
                    help="entry-age x elapsed-age block (production: varying).")
    ap.add_argument("--use-d", action="store_true",
                    help="also fit per-athlete drift d_i (off by default; "
                         "the aging stage isolates the global aging block).")
    ap.add_argument("--solvers", nargs="+", default=["anderson"],
                    choices=["als", "anderson"])
    ap.add_argument("--n-random", type=int, default=10,
                    help="random restarts around the baseline warm anchor.")
    ap.add_argument("--jit-u", type=float, default=0.10)
    ap.add_argument("--jit-v", type=float, default=0.10)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--no-v-traces", action="store_true")
    ap.add_argument("--no-register", action="store_true")
    S.add_spec_args(ap)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    t0 = time.perf_counter()
    for name in names:
        spec = S.build_spec(name, min_race_count=args.min_race_count,
                            date_lo=args.date_lo, date_hi=args.date_hi,
                            min_runner=args.min_runner)
        run_aging_fit(
            spec, args.nu, n_knots=args.n_knots, gamma_form=args.gamma_form,
            use_d=args.use_d, solvers=tuple(args.solvers), n_random=args.n_random,
            jit_u=args.jit_u, jit_v=args.jit_v, seed0=args.seed0,
            max_iter=args.max_iter, tol=args.tol,
            keep_v_traces=not args.no_v_traces, register=not args.no_register,
        )
    print(f"\nDone. Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
