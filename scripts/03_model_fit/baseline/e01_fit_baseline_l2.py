"""Fit the absolute-best Gaussian (L2, nu=inf) rank-1 baseline per slice.

Task 1. For each slice, fits the baseline (u_i + v_j, all optional terms off)
under both solvers (ALS, Anderson) from a mean init + N random restarts, at a
very tight tolerance (default 1e-12), and:
  * registers the single overall-best fit as the essential L2 model under
        results/models/{slug}/baseline_L2_best__<hash>/
  * writes the ALS-vs-Anderson convergence comparison (per-init spread,
    per-iter loglik/rss traces, v-traces, cross-init/-solver agreement) under
        results/convergence/{slug}/L2/

No upstream dependency — run anytime.

Run::

    python scripts/03_model_fit/baseline/e01_fit_baseline_l2.py             # the six slices
    python scripts/03_model_fit/baseline/e01_fit_baseline_l2.py --slices Po10_W
    python scripts/03_model_fit/baseline/e01_fit_baseline_l2.py --slices WA_M --min-runner 10
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from baseline_common import slices as S          # noqa: E402
from baseline_common.fitting import run_essential_fit  # noqa: E402

INF = float("inf")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"],
                    help="slice names or 'all' (default = the six in SLICE_ORDER).")
    ap.add_argument("--n-random", type=int, default=5)
    ap.add_argument("--jit-u", type=float, default=0.10)
    ap.add_argument("--jit-v", type=float, default=0.10)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--no-v-traces", action="store_true",
                    help="skip the large v_traces.parquet.")
    ap.add_argument("--no-register", action="store_true")
    S.add_spec_args(ap)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    t0 = time.perf_counter()
    for name in names:
        spec = S.build_spec(name, min_race_count=args.min_race_count,
                            date_lo=args.date_lo, date_hi=args.date_hi,
                            min_runner=args.min_runner)
        run_essential_fit(
            spec, INF, stage="L2", study="baseline_l2",
            n_random=args.n_random, jit_u=args.jit_u, jit_v=args.jit_v,
            seed0=args.seed0, max_iter=args.max_iter, tol=args.tol,
            keep_v_traces=not args.no_v_traces, register=not args.no_register,
        )
    print(f"\nDone. Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
