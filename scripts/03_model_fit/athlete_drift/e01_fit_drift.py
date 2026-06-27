"""Fit the production per-athlete drift models at nu=8 (Anderson).

Two variants (choose with ``--variant``):

  * **full**  = baseline + aging (spline-4, varying gamma) + d_i
                -> warm from the slice's ``agingS4gv_nu8p00_best``
                -> results/models/{slug}/full_nu8p00_best__<hash>/
  * **drift** = baseline + d_i (no aging)
                -> warm from the slice's ``baseline_nu8p00_best``
                -> results/models/{slug}/drift_nu8p00_best__<hash>/

The EB drift prior variance ``omega_d2`` is learned (type-II MLE) from
``--omega-init`` (default = model default 1e-4); ``--freeze-eb`` locks it at the
init. ``--min-span`` is the career-span eligibility floor in years (production
default 1e-3, an effective-zero safeguard). These knobs enter the identity hash,
so non-default variants never overwrite the production point.

Each slice warm-starts from its variant's registered source fit, with
``--n-random`` (default 10) perturbed restarts probing init-invariance. Writes
per-iter convergence traces (loglik / omega_d2 / v_j) and registers the single
overall-best (max loglik) model.

ZERO-ARG RUN (VS Code play): fits the `full` variant on all eight slices.

Run::

    python scripts/03_model_fit/athlete_drift/e01_fit_drift.py                       # full, all 8
    python scripts/03_model_fit/athlete_drift/e01_fit_drift.py --variant both --slices Po10_M
    python scripts/03_model_fit/athlete_drift/e01_fit_drift.py --variant drift --slices ALL_B
    python scripts/03_model_fit/athlete_drift/e01_fit_drift.py --slices Po10_M --omega-init 1e-3 --freeze-eb
    python scripts/03_model_fit/athlete_drift/e01_fit_drift.py --slices WA_M --min-runner 10
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from baseline_common import slices as S          # noqa: E402
from drift_common.fitting import run_drift_fit    # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--variant", default="full", choices=["full", "drift", "both"],
                    help="full = baseline+aging+d_i; drift = baseline+d_i; "
                         "both fits each in turn (default: full).")
    ap.add_argument("--nu", type=float, default=8.0,
                    help="Student-t d.o.f. (production decision: nu=8).")
    # --- drift-prior knobs ---
    ap.add_argument("--omega-init", type=float, default=None,
                    help="iter-0 omega_d2 (EB prior variance); default = model "
                         "default 1e-4. Result is init-invariant; affects speed.")
    ap.add_argument("--freeze-eb", action="store_true",
                    help="freeze omega_d2 at --omega-init instead of EB-learning it.")
    ap.add_argument("--min-span", type=float, default=1e-3,
                    help="d_i career-span eligibility floor in years (default 1e-3).")
    # --- aging form (full variant only) ---
    ap.add_argument("--n-knots", type=int, default=4,
                    help="spline knots for the `full` aging block (production: 4).")
    ap.add_argument("--gamma-form", default="varying", choices=["scalar", "varying"])
    # --- solve / inits ---
    ap.add_argument("--solvers", nargs="+", default=["anderson"],
                    choices=["als", "anderson"])
    ap.add_argument("--n-random", type=int, default=10,
                    help="random restarts around the warm anchor (probe invariance).")
    ap.add_argument("--jit-u", type=float, default=0.10)
    ap.add_argument("--jit-v", type=float, default=0.10)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--no-v-traces", action="store_true")
    ap.add_argument("--no-register", action="store_true")
    S.add_spec_args(ap)
    args = ap.parse_args()

    variants = ["full", "drift"] if args.variant == "both" else [args.variant]
    names = S.resolve_names(args.slices, ap)
    t0 = time.perf_counter()
    for name in names:
        spec = S.build_spec(name, min_race_count=args.min_race_count,
                            date_lo=args.date_lo, date_hi=args.date_hi,
                            min_runner=args.min_runner)
        for variant in variants:
            run_drift_fit(
                spec, variant=variant, nu=args.nu,
                omega_d2_init=args.omega_init, freeze_eb=args.freeze_eb,
                d_min_span=args.min_span, n_knots=args.n_knots,
                gamma_form=args.gamma_form, solvers=tuple(args.solvers),
                n_random=args.n_random, jit_u=args.jit_u, jit_v=args.jit_v,
                seed0=args.seed0, max_iter=args.max_iter, tol=args.tol,
                keep_v_traces=not args.no_v_traces, register=not args.no_register,
            )
    print(f"\nDone. Total wall: {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
