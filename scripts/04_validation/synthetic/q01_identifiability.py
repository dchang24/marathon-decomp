"""q01 -- identifiability / correctness probe.

Answers, with data, the question raised during e01 calibration: is the low raw
r_v solver collapse, or the documented gauge non-identification (S7.1-S7.3)?

Two separate claims are measured per condition:

  CORRECTNESS (noiseless): with sigma -> 0 the model represents the data
  exactly, so a correct solver must drive RSS -> 0 (machine precision) and the
  fitted PREDICTIONS must equal the truth pointwise. Parameter coordinates then
  differ from truth only by the gauge null-space. rss_fit and pred_maxerr are
  the gold-standard correctness numbers; they do NOT depend on any gauge.

  RECOVERY (noisy): with noise we can only recover statistically. We compare
  u,v in the SAME gauge applied to truth and fit:
    * v: residualize on [1, t_j]            (kills the APC + drift date-tilt)
    * u: residualize on [1, b_i, tbar_i]    (kills the cohort/debut tilts)
  and report raw vs gauge-fixed Pearson so the gap attributable to the flat
  directions is explicit.

Grid: terms {rank1, full} x fill {0.5, 0.9, 0.99} x noise {none, gauss, t6}.
Gaussian/large-nu isolates any Student-t nonconvexity. Sampler = balanced.

Run:
  python scripts/04_validation/synthetic/q01_identifiability.py
  python scripts/04_validation/synthetic/q01_identifiability.py --seeds 5
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from synth_common import (
    GT,
    OUT_ROOT,
    Report,
    _athlete_debut_mean_centered,
    _poly_resid,
    _race_year_centered,
    anderson_fitter,
    make_synthetic,
    model_config,
    replace,
)
from marathon_decomp.models.anderson import ModelAnderson

NOISE = {"none": ("inf", 1e-14), "gauss": ("inf", None), "t6": (6.0, None)}
FILLS = [0.5, 0.9, 0.99]
TERMS = ["rank1", "full"]


def _corr(a, b) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _run(terms: str, fill: float, noise: str, seed: int) -> dict:
    nu_s, s2 = NOISE[noise]
    nu = float(nu_s)
    gt = GT if s2 is None else replace(GT, sigma2=s2)
    fd, truth = make_synthetic(I=400, J=40, fill=fill, sampler="balanced",
                               terms=terms, noise_nu=nu, seed=seed, gt=gt)
    m = ModelAnderson(fd, model_config(terms, nu=nu),
                      anderson_fitter(tol=1e-10, max_outer_iter=2000))
    m.fit()

    yhat = m.predict()
    rss = float(np.sum((fd.y - yhat) ** 2))
    pred_maxerr = float(np.max(np.abs(yhat - truth.signal)))   # vs noiseless truth surface

    full = terms == "full"
    K = m.config.degree if full else 0      # aging basis degree -> date-poly gauge
    yr = _race_year_centered(fd)
    debut_c, meant_c = _athlete_debut_mean_centered(fd)
    vf, vt = m.params["v"], truth.v
    uf, ut = m.params["u"], truth.u

    return dict(
        terms=terms, fill=fill, noise=noise, seed=seed,
        n_iter=m.fit_result.n_iter, converged=m.fit_result.converged,
        rss_fit=rss, pred_maxerr=pred_maxerr,
        sigma2_ratio=float(m.params["sigma2"] / gt.sigma2),
        r_v_raw=_corr(vf - vf.mean(), vt - vt.mean()),
        r_v_gauge=_corr(_poly_resid(vf, (yr, K)), _poly_resid(vt, (yr, K))),
        r_u_raw=_corr(uf - uf.mean(), ut - ut.mean()),
        r_u_gauge=_corr(_poly_resid(uf, (debut_c, K), (meant_c, 1 if full else 0)),
                        _poly_resid(ut, (debut_c, K), (meant_c, 1 if full else 0))),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--seed0", type=int, default=2025)
    args = ap.parse_args()

    rows = [
        _run(terms, fill, noise, args.seed0 + s)
        for terms in TERMS
        for fill in FILLS
        for noise in NOISE
        for s in range(args.seeds)
    ]
    df = pd.DataFrame(rows)
    agg = (df.groupby(["terms", "fill", "noise"], sort=False)
             .agg(n_iter=("n_iter", "mean"),
                  rss_fit=("rss_fit", "mean"),
                  pred_maxerr=("pred_maxerr", "max"),
                  sigma2_ratio=("sigma2_ratio", "mean"),
                  r_v_raw=("r_v_raw", "mean"),
                  r_v_gauge=("r_v_gauge", "mean"),
                  r_u_raw=("r_u_raw", "mean"),
                  r_u_gauge=("r_u_gauge", "mean"))
             .reset_index())

    rep = Report("q01 -- identifiability / correctness probe")
    rep.line(f"I=400 J=40, balanced sampler, {args.seeds} seeds/cell. "
             "noise=none is sigma^2=1e-14 (correctness); rss_fit & pred_maxerr "
             "are gauge-free. r_*_gauge residualizes v on [1,year], u on "
             "[1,debut,mean-date].")
    rep.h2("Correctness (noiseless) + recovery (noisy)")
    rep.table(agg, floatfmt="{:.4g}")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_ROOT / "q01_identifiability.csv", index=False)
    rep.save(OUT_ROOT / "q01_identifiability.md")
    rep.line(f"\nWrote -> {OUT_ROOT}")


if __name__ == "__main__":
    main()
