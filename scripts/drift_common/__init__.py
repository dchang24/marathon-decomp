"""Shared production-fit machinery for the per-athlete drift models.

`run_drift_fit` fits one of two production variants at the settled operating
point and registers the overall-best (max-loglik) model:

    full   = baseline + aging (spline-4, varying gamma) + d_i   (warm: agingS4gv)
    drift  = baseline + d_i                                      (warm: baseline)

See `fitting.py`.
"""
from .fitting import drift_cfg, drift_stem, load_aging_warmstart, run_drift_fit  # noqa: F401
