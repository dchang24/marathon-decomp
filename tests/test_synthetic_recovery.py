"""Parameter-recovery validation for Model and ModelAnderson.

Generates synthetic marathon data from the model's own forward pass with
calibrated parameters (matching the scale of the ALL_B_14-20 AxD nu6 reference
fit), then checks that the fitter recovers the true latent structure.

Three levels of validation:

1. Objective dominance — the fitted model achieves a log-likelihood on the
   training data at least as good as the true generating parameters.  This is
   the cleanest check: if the fitter genuinely minimises -loglik it must beat
   the oracle, because the oracle was optimal for the population, not this
   noisy sample.

2. Identifiable-quantity recovery — Pearson correlations of recovered vs true
   u_i, v_j, and d_i above thresholds calibrated to the design.

3. Scalar quantity accuracy — sigma² and the aging-curve curvature (the only
   gauge-invariant part for a poly-2 basis) within a loose absolute/relative band.

Design rationale
----------------
Recovery tests use I = 2000 athletes, J = 50 races, 10 % fill (N ≈ 10 000).
This gives n_per_race ≈ 200, SNR(v_j) ≈ v_std / (σ/√200) ≈ 10, and
n_per_athlete ≈ 5 (same as real production slices).  At this scale the full
model — including drift d_i — converges with Anderson in ≈ 170 iterations and
recovers all parameters to Pearson r > 0.95.

Why not I = 200?  With n_per_race ≈ 20, SNR(v_j) ≈ 3 and SE(v_j) exceeds the
signal std, so v_j recovery is inherently noisy regardless of convergence.  The
BCD slow-convergence from the d_i / v_j coupling also becomes the dominant
bottleneck: d_i · Ãn = d_i · (t_j − mean_i(t_j)) creates athlete-specific t_j
contributions; their per-race average is O(d_std / √n_j), which is comparable
to v_std when n_j = 20 but negligible when n_j = 200.  The β=0 gauge does NOT
fix this — it removes the slope of v_j on t_j from the aging block, but leaves
d_i unchanged.

Notes on identifiability:
  - mean(v) = 0 is imposed by gauge; we compare v after centring.
  - The linear slope of the aging curve is APC-unidentified (§4.2); we compare
    the curvature (2·phi2 for a raw poly-2 basis), which is gauge-invariant.
  - All recovery tests use apc_gauge_beta0=False during fitting to avoid
    per-iteration gauge shifts perturbing the convergence trajectory (the gauge
    is tested separately).
  - Debut date per athlete = min race date across that athlete's cells, matching
    real FitData semantics so A_n ≥ 0 without any filtering.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from marathon_decomp.data import FitData, SliceSpec
from marathon_decomp.models.anderson import AndersonFitterConfig, ModelAnderson
from marathon_decomp.models.model import (
    FitterConfig,
    Model,
    ModelConfig,
    _SEC_PER_YEAR,
)

# ---------------------------------------------------------------------------
# Reference parameter magnitudes (ALL_B_14-20 AxD nu6 anderson)
# ---------------------------------------------------------------------------
# u: mean=9.6046, std=0.1957  (log-seconds; exp(9.6) ≈ 14760 s ≈ 4.1 h)
# v: std=0.0192               (race difficulty spread)
# d: std=0.0414               (per-year drift slope, eligible athletes only)
# sigma2: 0.00084             (sigma ≈ 0.029 in log-time)
# omega_d2: 0.00192

TRUE_PARAMS = dict(
    u_mean=9.60, u_std=0.20,
    v_std=0.020,
    phi1=-0.010, phi2=0.0022,
    gamma=0.0015,
    d_std=np.sqrt(0.00192),      # ≈ 0.044 per year
    sigma2=0.00090,              # sigma ≈ 0.030
    omega_d2=0.00192,
)

SEED = 2025


# ---------------------------------------------------------------------------
# Cell sampler — connected bipartite spanning structure
# ---------------------------------------------------------------------------

def _sample_connected_cells(
    I: int, J: int, target_N: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (row, col) arrays for a connected bipartite athlete-race graph.

    Guarantees every athlete has ≥ 1 cell and every race has ≥ 1 athlete.
    """
    cells: set[tuple[int, int]] = set()
    pa = rng.permutation(I)
    pj = rng.permutation(J)
    for i in range(I):
        cells.add((int(pa[i]), int(pj[i % J])))
    for j in range(J):
        cells.add((int(pa[j % I]), int(pj[j])))
    while len(cells) < target_N:
        cells.add((int(rng.integers(0, I)), int(rng.integers(0, J))))
    arr = np.array(sorted(cells), dtype=np.int64)
    return arr[:, 0], arr[:, 1]


# ---------------------------------------------------------------------------
# Synthetic data factory (full DGP, including drift)
# ---------------------------------------------------------------------------

def _make_synthetic(
    I: int = 2000,
    J: int = 50,
    fill: float = 0.10,
    sigma2: float = TRUE_PARAMS["sigma2"],
    noise_nu: float = float("inf"),
    seed: int = SEED,
) -> tuple[FitData, dict]:
    """Build a synthetic FitData with all model terms, including drift d_i.

    Debut date per athlete = min race date in their cell set (FitData convention),
    so A_n ≥ 0 everywhere without filtering.  Returns (FitData, true_params).
    """
    rng = np.random.default_rng(seed)

    race_dates = pd.date_range("2015-01-01", "2024-12-31", periods=J).to_numpy()
    t_j = race_dates.astype("int64") * (1e-9 / _SEC_PER_YEAR)

    target_N = max(int(round(fill * I * J)), I + J + 10)
    row, col = _sample_connected_cells(I, J, target_N, rng)
    N = len(row)

    f_i = np.full(I, np.inf)
    np.minimum.at(f_i, row, t_j[col])
    A_n = t_j[col] - f_i[row]

    A_e = rng.uniform(18.0, 55.0, I)
    Ae_c_cell = (A_e - float(A_e.mean()))[row]

    n_per = np.bincount(row, minlength=I).astype(np.float64)
    An_bar_i = np.where(
        n_per > 0,
        np.bincount(row, weights=A_n, minlength=I) / np.maximum(n_per, 1.0),
        0.0,
    )
    An_tilde = A_n - An_bar_i[row]

    u_true = rng.normal(TRUE_PARAMS["u_mean"], TRUE_PARAMS["u_std"], I)
    v_raw = rng.normal(0.0, TRUE_PARAMS["v_std"], J)
    v_true = v_raw - v_raw.mean()
    d_true = rng.normal(0.0, TRUE_PARAMS["d_std"], I)
    phi1, phi2 = TRUE_PARAMS["phi1"], TRUE_PARAMS["phi2"]
    gamma_true = TRUE_PARAMS["gamma"]

    signal = (
        u_true[row] + v_true[col]
        + phi1 * A_n + phi2 * A_n ** 2
        + gamma_true * Ae_c_cell * A_n
        + d_true[row] * An_tilde
    )

    sigma = float(np.sqrt(sigma2))
    if np.isfinite(noise_nu):
        tau = rng.gamma(noise_nu / 2.0, 2.0 / noise_nu, N)
        noise = rng.normal(0.0, sigma / np.sqrt(tau))
    else:
        noise = rng.normal(0.0, sigma, N)

    y = signal + noise

    spec = SliceSpec(sex="M")
    fd = FitData(
        row_idx=row.astype(np.int64), col_idx=col.astype(np.int64),
        y=y.astype(np.float64), A_n=A_n.astype(np.float64),
        A_e=A_e.astype(np.float64),
        athlete_ids=np.arange(I, dtype=np.int64),
        athlete_sex=np.array(["M"] * I, dtype=object),
        athlete_country=np.array(["GB"] * I, dtype=object),
        race_ids=np.arange(J, dtype=np.int64),
        race_series=np.array(["synth"] * J, dtype=object),
        race_country=np.array(["GB"] * J, dtype=object),
        race_date=race_dates,
        I=I, J=J, N=N,
        spec=spec,
        data_version="synthetic", cache_key="synthetic",
        response_kind="log_time",
    )
    return fd, dict(
        u=u_true, v=v_true, d=d_true,
        phi1=phi1, phi2=phi2, gamma=gamma_true,
        sigma2=sigma2, omega_d2=TRUE_PARAMS["omega_d2"],
        signal=signal,
    )


# ---------------------------------------------------------------------------
# Model / fitter config helpers
# ---------------------------------------------------------------------------

def _model_config(**overrides) -> ModelConfig:
    """Full model config (all terms enabled) used throughout."""
    cfg = ModelConfig(
        use_phi12=True, use_gamma=True, use_d=True,
        basis_kind="poly", degree=2, orthogonalize=True,
        gamma_form="scalar",
        omega_d2_max=0.05,
        apc_gauge_beta0=False,
        runner_weight_n0=1.0,
    )
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


def _anderson_fitter(**overrides) -> AndersonFitterConfig:
    """Anderson fitter config: tol=1e-7 converges in ~170 iters at I=2000."""
    cfg = AndersonFitterConfig(max_outer_iter=500, tol=1e-7, verbose=0)
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _log_lik_gaussian(y: np.ndarray, yhat: np.ndarray, sigma2: float) -> float:
    N = len(y)
    r = y - yhat
    return float(-0.5 * (N * np.log(2 * np.pi * sigma2) + np.sum(r * r) / sigma2))


def _aging_curvature(model: Model) -> float:
    """Gauge-invariant curvature d²f/dA_n² at A_n = 2 yr.

    True value = 2·phi2 = 0.0044.  Central-difference cancels the linear term.
    """
    from marathon_decomp.aging import aging_curve_on_grid
    grid = np.array([1.90, 2.00, 2.10], dtype=np.float64)
    f = aging_curve_on_grid(model, grid)
    return float((f[2] - 2.0 * f[1] + f[0]) / 0.01)


# ---------------------------------------------------------------------------
# Test 1: Objective dominance (Gaussian)
# ---------------------------------------------------------------------------

def test_objective_dominance_gaussian():
    """Fitted model achieves higher log-likelihood than the true DGP parameters.

    The MLE on the noisy training sample must beat the oracle (population-optimal)
    parameters.  Fails only if the optimiser has a sign error or cannot descend.
    """
    fd, tp = _make_synthetic(noise_nu=float("inf"))

    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    ll_fitted = model.log_lik()
    ll_true = _log_lik_gaussian(fd.y, tp["signal"], tp["sigma2"])

    assert ll_fitted >= ll_true, (
        f"Fitter loglik={ll_fitted:.4f} < oracle loglik={ll_true:.4f}. "
        "The optimiser failed to beat the true DGP — convergence or sign error."
    )


# ---------------------------------------------------------------------------
# Test 2: u_i Pearson correlation (Gaussian)
# ---------------------------------------------------------------------------

def test_ui_recovery_gaussian():
    """u_i correlates very strongly with the true athlete abilities.

    With n_per_athlete ≈ 5, u_std=0.20, sigma≈0.030: theoretical r > 0.999.
    Threshold 0.99 is conservative, accounting for aging/drift block interactions.
    """
    fd, tp = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    r = float(np.corrcoef(tp["u"], model.params["u"])[0, 1])
    assert r > 0.99, f"u_i Pearson r = {r:.4f}, expected > 0.99."


# ---------------------------------------------------------------------------
# Test 3: v_j Pearson correlation (Gaussian)
# ---------------------------------------------------------------------------

def test_vj_recovery_gaussian():
    """v_j correlates strongly with the true race difficulties.

    At I=2000, n_per_race ≈ 200: SE(v_j) ≈ σ/√200 ≈ 0.002, SNR ≈ 10,
    theoretical Pearson r ≈ 0.995.  Threshold 0.90 is conservative.
    """
    fd, tp = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    r = float(np.corrcoef(tp["v"], model.params["v"])[0, 1])
    assert r > 0.90, f"v_j Pearson r = {r:.4f}, expected > 0.90."


# ---------------------------------------------------------------------------
# Test 4: d_i Pearson correlation (Gaussian)
# ---------------------------------------------------------------------------

def test_di_recovery_gaussian():
    """Eligible d_i correlates with the true drift slopes.

    With n_per_athlete ≈ 5 and An_tilde spanning ~4 yr per athlete:
    SE(d_i) ≈ σ / (std(An_tilde_i) · √5) ≈ 0.007, d_std ≈ 0.044, r_theory ≈ 0.99.
    Threshold 0.85 accounts for the EB shrinkage flattening the estimates.
    """
    fd, tp = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    elig = model.eligible_d
    r = float(np.corrcoef(tp["d"][elig], model.params["d"][elig])[0, 1])
    assert r > 0.85, f"d_i Pearson r (eligible) = {r:.4f}, expected > 0.85."


# ---------------------------------------------------------------------------
# Test 5: sigma² accuracy (Gaussian)
# ---------------------------------------------------------------------------

def test_sigma2_recovery_gaussian():
    """Estimated sigma² is within 40 % of the true value.

    The uncorrected MLE is biased down by ~(I+J)/N ≈ 20 %; the band (0.60, 1.30)
    accommodates this without requiring a degrees-of-freedom correction.
    """
    fd, tp = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    ratio = model.params["sigma2"] / tp["sigma2"]
    assert 0.60 <= ratio <= 1.30, (
        f"sigma² ratio fitted/true = {ratio:.3f}, expected in [0.60, 1.30]. "
        f"(fitted={model.params['sigma2']:.5f}, true={tp['sigma2']:.5f})"
    )


# ---------------------------------------------------------------------------
# Test 6: Aging curve curvature (gauge-invariant)
# ---------------------------------------------------------------------------

def test_aging_curvature_gaussian():
    """The fitted aging curve's curvature matches the true value to < 0.002.

    Curvature = 2·phi2 = 0.0044 is APC-gauge invariant.  Error < 0.002.
    """
    fd, _ = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    curv_fit = _aging_curvature(model)
    curv_true = 2.0 * TRUE_PARAMS["phi2"]   # 0.0044

    assert abs(curv_fit - curv_true) < 0.002, (
        f"Aging curvature: fitted={curv_fit:.5f}, true={curv_true:.5f}, "
        f"error={abs(curv_fit - curv_true):.5f} > 0.002"
    )


# ---------------------------------------------------------------------------
# Test 7: ALS and Anderson converge to the same fixed point
# ---------------------------------------------------------------------------

def test_als_anderson_agree():
    """Model (plain ALS) and ModelAnderson reach the same fixed point.

    Both share the same BCD step; Anderson only extrapolates faster.
    At the fixed point v_j and sigma² must be identical to high precision.
    Uses a smaller dataset (I=500) so plain ALS converges within the budget.
    """
    fd, _ = _make_synthetic(I=500, J=50, noise_nu=float("inf"), seed=SEED + 7)
    cfg = _model_config(nu=float("inf"))

    als = Model(fd, cfg, FitterConfig(max_outer_iter=1000, tol=1e-7, verbose=0))
    als.fit()

    anderson = ModelAnderson(fd, cfg, _anderson_fitter(max_outer_iter=1000, tol=1e-7))
    anderson.fit()

    r_v = float(np.corrcoef(als.params["v"], anderson.params["v"])[0, 1])
    assert r_v > 0.9999, (
        f"ALS vs Anderson v_j Pearson r = {r_v:.6f}, expected > 0.9999"
    )

    sigma_ratio = als.params["sigma2"] / anderson.params["sigma2"]
    assert 0.999 <= sigma_ratio <= 1.001, (
        f"sigma² ALS/Anderson ratio = {sigma_ratio:.5f}, expected ≈ 1.000"
    )


# ---------------------------------------------------------------------------
# Test 8: Student-t objective dominance
# ---------------------------------------------------------------------------

def test_objective_dominance_t():
    """Under Student-t noise (nu=6) the fitter still beats the true DGP params.

    Confirms IRLS correctly implements the t log-likelihood.
    """
    from scipy.special import gammaln

    fd, tp = _make_synthetic(noise_nu=6.0)
    model = ModelAnderson(fd, _model_config(nu=6.0), _anderson_fitter())
    model.fit()

    ll_fitted = model.log_lik()

    nu = 6.0
    s2 = tp["sigma2"]
    r = fd.y - tp["signal"]
    N = fd.N
    c = N * (
        gammaln((nu + 1) / 2) - gammaln(nu / 2)
        - 0.5 * np.log(nu * np.pi * s2)
    )
    ll_true = float(c - 0.5 * (nu + 1) * np.sum(np.log1p(r * r / (nu * s2))))

    assert ll_fitted >= ll_true, (
        f"Student-t fitter loglik={ll_fitted:.4f} < oracle loglik={ll_true:.4f}. "
        "IRLS may not be correctly implementing the t loss."
    )


# ---------------------------------------------------------------------------
# Test 9: Stationarity / gradient check at convergence
# ---------------------------------------------------------------------------

def test_stationarity_at_convergence():
    """At convergence the RSS gradient w.r.t. v_j is numerically small.

    Forward-difference on five randomly chosen v_j: perturbing by ε must not
    decrease RSS (optimality condition at an interior minimum).
    |ΔRSS|/ε < 0.20 is a loose bound; a converged fitter should be much tighter.
    """
    fd, _ = _make_synthetic(I=500, J=30, noise_nu=float("inf"), seed=SEED + 9)
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    model.fit()

    rng = np.random.default_rng(SEED + 99)
    eps = 1e-5
    rss0 = float(np.sum((fd.y - model.predict()) ** 2))
    probe_races = rng.choice(fd.J, size=5, replace=False)

    for j in probe_races:
        model.params["v"][j] += eps
        rss_plus = float(np.sum((fd.y - model.predict()) ** 2))
        model.params["v"][j] -= eps

        grad = (rss_plus - rss0) / eps
        assert abs(grad) < 0.20, (
            f"RSS gradient w.r.t. v[{j}] ≈ {grad:.4f}, expected |·| < 0.20. "
            "BCD may not have reached a stationary point."
        )


# ---------------------------------------------------------------------------
# Test 10: Convergence flag
# ---------------------------------------------------------------------------

def test_convergence_flag():
    """FitResult.converged is True within the iteration budget.

    Anderson at tol=1e-7 converges in ≈ 170 iterations for this design.
    """
    fd, _ = _make_synthetic(noise_nu=float("inf"))
    model = ModelAnderson(fd, _model_config(nu=float("inf")), _anderson_fitter())
    result = model.fit()

    assert result.converged, (
        f"ModelAnderson did not converge in {result.n_iter} iterations "
        "(expected ~170 at tol=1e-7 for I=2000, J=50, 10% fill)."
    )
