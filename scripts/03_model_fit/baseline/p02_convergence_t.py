"""Per-slice Student-t convergence diagnostics: metric-vs-iter, ALS | Anderson.

Thin wrapper over ``baseline_common.convergence_plot`` for the ``nu_selected``
stage (the finite-nu essential fit). Same 6x2 figure as the L2 plotter, but the
lines now include the warm-start strategies tried by ``e02_fit_baseline_t.py``
(``mean``/``rand*`` cold, plus ``l2_warm``/``l2_rand*``/``step{nu}`` warm) — so
the figure doubles as the init-strategy comparison (warm starts should sit on the
log floor from iteration 1).

Plots every slice that has ``nu_selected`` traces under
``results/convergence/{slug}/nu_selected/`` to ``.../nu_selected/convergence.png``.
Safe to run mid-sweep: a slice is plotted only once BOTH trace files exist, so
in-progress slices are simply skipped.

Argument-free (VS Code ▶ works).

Run::

    python scripts/03_model_fit/baseline/p02_convergence_t.py            # all done slices
    python scripts/03_model_fit/baseline/p02_convergence_t.py --slices Po10_M_14-25_mrc2
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from baseline_common.convergence_plot import main  # noqa: E402

if __name__ == "__main__":
    main(stage="nu_selected")
