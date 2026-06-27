"""The aging-form candidate grid (model-selection stage).

This stage layers the **global aging block** on the bare rank-1 baseline
(``u_i + v_j``). Per-athlete drift ``d_i`` is held **off** so the fitted curve
isolates the aggregate aging signal with no competing per-athlete slope (caveat:
with ``d`` off the global curve absorbs the within-athlete aging it would
otherwise share with ``d`` -- noted in the analysis log when first run). The
stretch term ``s_j`` no longer exists in the model.

Three axes are swept:

  * **basis** -- the phi-block parametric form ``theta_aging @ B(A_n)``:
      - ``polyD``   raw monomials of degree ``D`` (QR-orthonormalized at fit
        time for conditioning; orthogonalization fixes *conditioning*, not
        boundary *wiggle* -- high-degree polynomials still oscillate in the
        sparse high-A_n tail; that is exactly what we want to see).
      - ``splineK`` natural cubic regression spline with ``K`` knots (linear
        beyond the boundary knots, so it extrapolates sanely).
  * **gamma_form** -- the entry-age x elapsed-age interaction block:
      - ``off``      no entry-age term (``use_gamma=False``)
      - ``scalar``   ``gamma * (A_e - mean) * A_n``            (canonical)
      - ``varying``  ``sum_k gamma_k * (A_e - mean) * B_k(A_n)``
  * **nu** -- Student-t d.o.f.: ``inf`` (Gaussian cross-check + warm-start
    anchor) and ``8.0`` (the settled production operating point).

``poly2 + gamma=scalar`` is the canonical baseline-equation aging block and
sits inside the grid as the reference point.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from marathon_decomp import ModelConfig

# -- axes -------------------------------------------------------------------


@dataclass(frozen=True)
class Basis:
    """One phi-block parametric form."""

    name: str           # readable tag, e.g. "poly3", "spline5"
    kind: str           # "poly" | "spline"
    degree: int = 2     # poly only
    n_knots: int = 5    # spline only


# poly1 (pure linear) is omitted -- it cannot represent the improvement-then-
# decline shape and is a known-wrong term.
POLY_DEGREES = (2, 3, 4, 5, 6)
SPLINE_KNOTS = (3, 4, 5, 6)

BASES: tuple[Basis, ...] = (
    *(Basis(f"poly{d}", "poly", degree=d) for d in POLY_DEGREES),
    *(Basis(f"spline{k}", "spline", n_knots=k) for k in SPLINE_KNOTS),
)
BASES_BY_NAME: dict[str, Basis] = {b.name: b for b in BASES}

GAMMA_FORMS = ("off", "scalar", "varying")
NU_GRID = (float("inf"), 8.0)

# The canonical baseline-equation aging block, for reference highlighting.
REFERENCE = ("poly2", "scalar")
REFERENCE_CAND = "poly2-gscalar"


# -- config construction ----------------------------------------------------


def build_config(basis: Basis, gamma_form: str, nu: float) -> ModelConfig:
    """Aging-only model config: rank-1 + phi-block + (optional) gamma.

    ``use_d`` is off by design for this stage; ``s_j`` no longer exists.
    """
    if gamma_form not in GAMMA_FORMS:
        raise ValueError(f"unknown gamma_form {gamma_form!r}")
    use_gamma = gamma_form != "off"
    # gamma_form Literal only accepts scalar|varying; default to scalar when off.
    gf = gamma_form if use_gamma else "scalar"
    return ModelConfig(
        use_phi12=True,
        use_gamma=use_gamma,
        use_d=False,
        basis_kind=basis.kind,
        degree=basis.degree,
        n_knots=basis.n_knots,
        gamma_form=gf,
        nu=float(nu),
    )


# -- labels -----------------------------------------------------------------


def nu_label(nu: float) -> str:
    """``inf`` for Gaussian else the numeric value -- for tidy table columns."""
    return "inf" if not math.isfinite(float(nu)) else f"{float(nu):g}"


def nu_tag(nu: float) -> str:
    """Filesystem-safe nu tag for fit-payload filenames: ``L2`` / ``nu8p00``."""
    return "L2" if not math.isfinite(float(nu)) else f"nu{float(nu):.2f}".replace(".", "p")


def cand_label(basis: Basis, gamma_form: str) -> str:
    """Short candidate label used in diagnostic tables/plots, e.g. spline5-gvarying."""
    return f"{basis.name}-g{gamma_form}"


def form_code(basis: Basis) -> str:
    """Compact parametric-form code: P{d} / S{k} (poly4 -> P4, spline5 -> S5)."""
    head = "P" if basis.kind == "poly" else "S"
    n = basis.degree if basis.kind == "poly" else basis.n_knots
    return f"{head}{n}"


def fit_stem(basis: Basis, gamma_form: str, nu: float, solver: str) -> str:
    """Readable per-cell payload filename stem ``{cand}_{nutag}_{solver}``."""
    return f"{cand_label(basis, gamma_form)}_{nu_tag(nu)}_{solver}"
