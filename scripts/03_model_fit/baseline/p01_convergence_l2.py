"""Per-slice L2 (Gaussian) convergence diagnostics: metric-vs-iter, ALS | Anderson.

Thin wrapper over ``baseline_common.convergence_plot`` for the ``L2`` stage.
Plots every slice that has L2 traces under ``results/convergence/{slug}/L2/``
(written by ``e01_fit_baseline_l2.py``) to ``.../L2/convergence.png``.

Argument-free (VS Code ▶ works). See the shared module for the figure spec.

Run::

    python scripts/03_model_fit/baseline/p01_convergence_l2.py            # all slices
    python scripts/03_model_fit/baseline/p01_convergence_l2.py --slices ALL_M_14-25_mrc2
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from baseline_common.convergence_plot import main  # noqa: E402

if __name__ == "__main__":
    main(stage="L2")
