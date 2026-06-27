"""Shared fitting primitives for the production aging-block fits.

Parallels ``baseline_common`` but for the aging stage (rank-1 ``u+v`` + the
settled production aging block: natural cubic spline, 4 knots, varying gamma).
Reuses ``baseline_common.slices`` / ``baseline_common.inits`` for the slice and
init plumbing; only the model config, the warm source, and the per-iter capture
(adds ``theta_aging`` / ``gamma`` to the recorded ``v`` traces) differ.
"""
