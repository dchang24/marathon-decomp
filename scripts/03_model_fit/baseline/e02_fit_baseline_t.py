"""Fit the absolute-best Student-t rank-1 baseline at the decided nu (=8).

Task 2 (final step). Refits the baseline at the chosen nu under BOTH solvers at a
very tight tolerance, comparing init strategies and keeping the single overall-best
fit. Besides the always-present **cold start** (mean + N random restarts — the
default when nothing else has been run), it leverages the already-fitted Gaussian
(L2) model with three warm strategies (``--cold-only`` disables them):

  1. ``l2_warm``   — anchor directly at the L2 (u, v).
  2. ``l2_rand*``  — N random restarts around the L2 anchor.
  3. ``step{S}``   — fit a 'stepping-stone' nu=S (``--stepping-nu``, default 15)
                     warm from L2, then warm-start the target nu from that.

Writes (per slice):
  * essential model -> results/models/{slug}/baseline_nu<p>_best__<hash>/
  * convergence     -> results/convergence/{slug}/nu_selected/  (the per-init
                       spread/traces let you compare the strategies; ``best_init``
                       in the manifest names the winner.)

ZERO-ARG RUN (VS Code ▶): fits ALL eight slices at nu=8 with warm starts on,
stepping-stone nu=15. No CLI needed. The L2 warm source is the registered
``baseline_L2_best`` per slice (run ``e01_fit_baseline_l2.py`` first); slices
without one silently fall back to cold start.

Run::

    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py             # all 8, nu=8, warm
    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py --slices Po10_M
    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py --cold-only  # no L2 warm
    python scripts/03_model_fit/baseline/e02_fit_baseline_t.py --nu 5 --stepping-nu 12
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR              # noqa: E402
from baseline_common import slices as S                     # noqa: E402
from baseline_common.fitting import run_essential_fit       # noqa: E402

SELECT_ROOT = RESULTS_DIR / "model_selection" / "baseline"


def _read_selected_nu(slug: str) -> float:
    path = SELECT_ROOT / slug / "selected_nu.csv"
    if not path.is_file():
        raise SystemExit(
            f"ERROR: no selected_nu.csv for '{slug}' at {path}.\n"
            f"Run the selection sweep first (same slice overrides):\n"
            f"    python scripts/02_model_selection/baseline/e01_nu_cv.py "
            f"--slices <name>")
    row = pd.read_csv(path).iloc[0]
    return float(row["selected_nu"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--nu", type=float, default=8.0,
                    help="fit at this explicit nu (decision: nu=8). "
                         "Pass None-like via --use-selected to read selected_nu.csv.")
    ap.add_argument("--use-selected", action="store_true",
                    help="read selected_nu.csv per slice instead of the fixed --nu.")
    ap.add_argument("--cold-only", action="store_true",
                    help="disable the L2 warm-start strategies (cold start only).")
    ap.add_argument("--stepping-nu", type=float, default=15.0,
                    help="stepping-stone nu for the L2->step->target warm chain.")
    ap.add_argument("--n-random", type=int, default=5)
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
        nu = _read_selected_nu(S.slug(spec)) if args.use_selected else args.nu
        run_essential_fit(
            spec, nu, stage="nu_selected", study="baseline_t",
            n_random=args.n_random, jit_u=args.jit_u, jit_v=args.jit_v,
            seed0=args.seed0, max_iter=args.max_iter, tol=args.tol,
            keep_v_traces=not args.no_v_traces, register=not args.no_register,
            warm_from_l2=not args.cold_only, stepping_nu=args.stepping_nu,
        )
    print(f"\nDone. Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
