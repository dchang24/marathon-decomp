"""Model variants. `Model` is the canonical baseline+aging+drift model fit by
plain block-coordinate descent (ALS); `ModelAnderson` is the same model with an
Anderson-accelerated outer loop."""
from .model import FitterConfig, Model, ModelConfig
from .anderson import AndersonFitterConfig, ModelAnderson

__all__ = [
    "Model",
    "ModelConfig",
    "FitterConfig",
    "ModelAnderson",
    "AndersonFitterConfig",
]
