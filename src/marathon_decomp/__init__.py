"""marathon-decomp: latent factor decomposition of marathon finish times.

Top-level imports expose the user-facing surface so a script needs a single
line, e.g.:

    from marathon_decomp import SliceSpec, load_slice, Model, ModelConfig
"""
from . import registry
from .aging import (
    aging_curve_from_payload,
    aging_curve_on_grid,
    build_poly_basis,
    default_knots_from_An,
    entry_age_curve_from_payload,
    entry_age_curve_on_grid,
    ncs_basis,
)
from .base import BaseModel, FitResult, SaveSpec
from .crossval import (
    assign_folds,
    heldout_logdensity,
    subset_fitdata,
    t_logpdf,
)
from .data import FitData, SliceSpec, clear_cache, load_slice
from .resample import bayesian_athlete_weights, boot_cell_weights
from .models import (
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    ModelConfig,
)

__all__ = [
    # data
    "SliceSpec",
    "FitData",
    "load_slice",
    "clear_cache",
    # cross-validation
    "assign_folds",
    "subset_fitdata",
    "heldout_logdensity",
    "t_logpdf",
    # resampling / bootstrap
    "bayesian_athlete_weights",
    "boot_cell_weights",
    # base
    "BaseModel",
    "FitResult",
    "SaveSpec",
    # models + configs
    "Model",
    "ModelConfig",
    "FitterConfig",
    "ModelAnderson",
    "AndersonFitterConfig",
    # aging basis + curve reconstruction
    "build_poly_basis",
    "ncs_basis",
    "default_knots_from_An",
    "aging_curve_on_grid",
    "entry_age_curve_on_grid",
    "aging_curve_from_payload",
    "entry_age_curve_from_payload",
    # run registry for sweeps
    "registry",
]
