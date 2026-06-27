"""Per-slice convergence diagnostics for the `full` drift fit (baseline+aging+d_i).

Thin wrapper over ``drift_common.convergence_plot`` for the ``drift_full`` stage.
Same metric-vs-iter, ALS | Anderson figure as the baseline plotter, plus the
EB-prior ``|omega_d2 - omega*|`` convergence row. The inits are the warm anchor
``warm`` (warm-started from ``agingS4gv_nu8p00_best``) plus the ``rand*``
perturbed restarts -- warm-started curves should sit on the log floor from the
first iteration, and every restart should land on the same fixed point.

Plots every slice with ``drift_full`` traces under
``results/convergence/{slug}/drift_full/`` to ``.../drift_full/convergence.png``.
Safe to run mid-sweep: a slice is plotted only once BOTH trace files exist.

Argument-free (VS Code play works).

Run::

    python scripts/03_model_fit/athlete_drift/p01_convergence_full.py            # all done slices
    python scripts/03_model_fit/athlete_drift/p01_convergence_full.py --slices Po10_M_14-25_mrc2
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from drift_common.convergence_plot import main  # noqa: E402

if __name__ == "__main__":
    main(stage="drift_full")
