"""e01 -- faithful full-model recovery on high-fill, small synthetic data.

Test A of the synthetic validation (see this dir's README). Generates data from
the model's OWN forward pass with ALL_B-scale ground truth (full DGP: aging +
entry-age gamma + per-athlete drift + Student-t noise), fits it back, and over
many seeds reports BOTH:

  * statistical recovery  -- gauge-aligned r/RMSE on u,v,d; sigma^2 ratio; the
    APC-invariant aging curvature error; omega_d^2 ratio; Brinker MRE.
  * numerical convergence -- converged flag, n_iter, the per-race stationarity
    residual (S5.1), oracle-dominance margin, monotone log-lik trace, and
    ALS == Anderson agreement at the fixed point.

The default designs keep FILL high (>=0.5) so the SNR(v) and d_i/v_j coupling
problems that plague tiny-but-sparse data are absent, letting the dataset be
genuinely small while still well covered:

  base  : I=400, J=40, fill=0.50  ->  N~8000,  n/race~200, n/runner~20
  tiny  : I=150, J=25, fill=0.60  ->  N~2250,  n/race~90,  n/runner~15

Run:
  python scripts/04_validation/synthetic/e01_recovery.py                  # all designs
  python scripts/04_validation/synthetic/e01_recovery.py --seeds 30 --nu 6
  python scripts/04_validation/synthetic/e01_recovery.py --designs base   # one design

Argument-free = both designs, 20 seeds, Gaussian + nu=6, Anderson; one ALS run
per (design, first seed) for the agreement check. Writes to
results/validation/synthetic/.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from synth_common import (
    OUT_ROOT,
    Report,
    aging_curvature,
    als_fitter,
    anderson_fitter,
    convergence_diagnostics,
    fit_with_trace,
    make_synthetic,
    model_config,
    recovery_metrics,
    upsert_csv,
)
from marathon_decomp.models.anderson import ModelAnderson
from marathon_decomp.models.model import Model

DESIGNS = {
    "base": dict(I=400, J=40, fill=0.50),
    "tiny": dict(I=150, J=25, fill=0.60),
}

KEY = ["design", "sampler", "nu", "seed"]


def _run_one(design: str, dims: dict, sampler: str, nu: float, seed: int) -> dict:
    fd, truth = make_synthetic(
        **dims, sampler=sampler, terms="full", noise_nu=nu, seed=seed,
    )
    # cap raised to 2000: the d_i/v_j BCD coupling is slowest at moderate fill
    # (q01: ~1120 iters at fill=0.5 vs 9 at fill=0.9), so 500 can truncate.
    model = ModelAnderson(fd, model_config("full", nu=nu),
                          anderson_fitter(max_outer_iter=2000))
    trace = fit_with_trace(model)

    row = dict(design=design, sampler=sampler, nu=nu, seed=seed,
               I=fd.I, J=fd.J, N=fd.N, fill=round(fd.N / (fd.I * fd.J), 4))
    row.update(recovery_metrics(model, fd, truth))      # curvature_err = full-model (confounded)
    row.update(convergence_diagnostics(model, fd, truth, trace))

    # Aging curvature is confounded with d_i in the full model (d_i*t_j ==
    # -2*phi2*b_i*t_j); it is identified only in the no-drift model, which is how
    # production characterizes aging (agingS4gv). Refit without drift for the
    # honest curvature-recovery number.
    nod = ModelAnderson(fd, model_config("full", nu=nu, use_d=False),
                        anderson_fitter(max_outer_iter=2000))
    nod.fit()
    row["curvature_err_nodrift"] = aging_curvature(nod) - truth.gt.curvature
    return row


def _als_anderson_agreement(design: str, dims: dict, nu: float, seed: int) -> dict:
    """Fit the same data with plain ALS and Anderson; compare at the fixed point."""
    fd, truth = make_synthetic(**dims, sampler="balanced", terms="full",
                               noise_nu=nu, seed=seed)
    cfg = model_config("full", nu=nu)
    als = Model(fd, cfg, als_fitter())
    als.fit()
    anderson = ModelAnderson(fd, cfg, anderson_fitter(max_outer_iter=1000))
    anderson.fit()
    r_v = float(np.corrcoef(als.params["v"], anderson.params["v"])[0, 1])
    return dict(design=design, nu=nu, seed=seed,
                als_n_iter=als.fit_result.n_iter,
                and_n_iter=anderson.fit_result.n_iter,
                agree_r_v=r_v,
                sigma2_ratio=float(als.params["sigma2"] / anderson.params["sigma2"]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", nargs="+", default=list(DESIGNS),
                    choices=list(DESIGNS))
    ap.add_argument("--samplers", nargs="+", default=["staggered", "balanced"],
                    choices=["staggered", "balanced", "random"])
    ap.add_argument("--nus", nargs="+", type=float, default=[float("inf"), 6.0])
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--seed0", type=int, default=2025)
    args = ap.parse_args()

    rep = Report("e01 -- full-model recovery (Test A)")
    rows: list[dict] = []
    for design in args.designs:
        dims = DESIGNS[design]
        for sampler in args.samplers:
            for nu in args.nus:
                for s in range(args.seeds):
                    rows.append(_run_one(design, dims, sampler, nu,
                                         args.seed0 + s))
    df = pd.DataFrame(rows)

    # --- aggregate: mean +- sd across seeds per (design, sampler, nu) ---
    metric_cols = [c for c in df.columns if c not in KEY + ["I", "J", "N", "fill"]]
    agg = (df.groupby(["design", "sampler", "nu"])[metric_cols]
             .agg(["mean", "std"]).reset_index())
    agg.columns = ["__".join(c).rstrip("_") for c in agg.columns]

    rep.h2("Per-seed rows")
    rep.line(f"{len(df)} fits over designs={args.designs}, "
             f"samplers={args.samplers}, nus={args.nus}, seeds={args.seeds}.")

    rep.h2("Recovery (mean across seeds; r_v/r_u gauge-fixed, raw shown too)")
    show = ["design", "sampler", "nu", "r_v__mean", "r_v_raw__mean",
            "r_u__mean", "r_d__mean", "sigma2_ratio__mean", "omega_d2_ratio__mean",
            "mre__mean"]
    rep.table(agg[[c for c in show if c in agg.columns]])

    rep.h2("Aging curvature error (full model confounds with drift; "
           "no-drift refit is the identified value; true=0.0044)")
    showk = ["design", "sampler", "nu", "curvature_err__mean",
             "curvature_err_nodrift__mean"]
    rep.table(agg[[c for c in showk if c in agg.columns]])

    rep.h2("Convergence (mean across seeds)")
    showc = ["design", "sampler", "nu", "converged__mean", "n_iter__mean",
             "stationarity_v__mean", "oracle_margin__mean", "mono_min_step__mean"]
    rep.table(agg[[c for c in showc if c in agg.columns]])

    # --- ALS vs Anderson agreement (one per design, Gaussian, first seed) ---
    rep.h2("ALS == Anderson at the fixed point")
    agree = pd.DataFrame([
        _als_anderson_agreement(d, DESIGNS[d], float("inf"), args.seed0)
        for d in args.designs
    ])
    rep.table(agree)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    upsert_csv(df, OUT_ROOT / "e01_recovery.csv", KEY)
    agg.to_csv(OUT_ROOT / "e01_recovery_agg.csv", index=False)
    agree.to_csv(OUT_ROOT / "e01_als_anderson.csv", index=False)
    rep.save(OUT_ROOT / "e01_recovery.md")
    rep.line(f"\nWrote -> {OUT_ROOT}")


if __name__ == "__main__":
    main()
