"""Shared convergence-plotter for the per-athlete drift fits (`full` / `drift`).

The drift stage (``drift_common.fitting.run_drift_fit``) writes the same per-init
diagnostic traces as the baseline stage, so this reuses the ``baseline_common``
metric machinery (reference solution, loglik/v error frames, log-y panel plotter)
and only adds the one diagnostic unique to this stage: the per-iteration
trajectory of the EB drift-prior variance ``omega_d2`` (type-II MLE-learned by
default; the quantity ``02_model_selection/athlete_drift`` characterized). The
figure is otherwise identical in spirit to the baseline one.

    stage="drift_full"  -> results/convergence/{slug}/drift_full/   (baseline+aging+d_i)
    stage="drift_drift" -> results/convergence/{slug}/drift_drift/  (baseline+d_i)

The thin runnable wrappers ``03_model_fit/athlete_drift/p01_convergence_full.py``
and ``p02_convergence_drift.py`` just call :func:`main` with their stage.

Figure: a 7-row x 2-column grid (columns = ALS | Anderson) of log-y metrics vs
outer iteration, one coloured line per init (the warm ``warm`` anchor + the
``rand*`` perturbed restarts). Rows:
  (1) loglik gap to best       (2) mean |v-v*|         (3) max |v-v*|
  (4) 1-corr(v,v*) pearson     (5) 1-corr spearman     (6) |omega_d2 - omega*|
  (7) rel. dloglik (+tol line)
The reference (v*, omega*) is the single max-final-loglik run across all
solvers x inits; its own curves rest on the log floor.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Reuse the baseline metric machinery verbatim — the trace schema is identical.
from baseline_common.convergence_plot import (  # noqa: E402
    OUT_ROOT,
    SOLVERS,
    _DEFAULT_TOL,
    _FLOOR,
    _PALETTE,
    _plot_metric,
    _stage_dir,
    _stage_nu,
    _trail_int,
    has_traces,
    loglik_metric_frame,
    reference,
    v_error_frame,
)


# ── discovery (drift stages live under their own sub-dirs) ───────────
def discover_slices(stage: str) -> list[str]:
    if not OUT_ROOT.is_dir():
        return []
    return sorted(d.name for d in OUT_ROOT.iterdir()
                  if d.is_dir() and has_traces(d.name, stage))


# ── omega_d2 convergence frame ───────────────────────────────────────
def omega_metric_frame(ll: pd.DataFrame, ref) -> pd.DataFrame:
    """Per (solver, init, iter): |omega_d2 - omega*|, omega* = ref run's final."""
    run = ll[(ll.solver == ref.solver) & (ll.init == ref.init)].sort_values("iter")
    omega_star = float(run["omega_d2"].iloc[-1])
    out = ll.sort_values("iter")[["solver", "init", "iter", "omega_d2"]].copy()
    out["omega_gap"] = (out["omega_d2"] - omega_star).abs()
    return out


# ── init ordering: warm anchor first, then rand0..randN ──────────────
def _init_order(inits: list[str]) -> list[str]:
    def key(s: str) -> tuple[int, int]:
        if s == "warm":
            return (0, -1)
        if s.startswith("rand"):
            return (1, _trail_int(s))
        return (2, _trail_int(s))
    return sorted(inits, key=key)


# ── plotting ─────────────────────────────────────────────────────────
def build_figure(slug, stage, ll_m, v_err, om_m, ref, *, tol, floor, nu=None):
    v_err = v_err.copy()
    v_err["one_minus_pearson"] = 1.0 - v_err["corr_pearson"]
    v_err["one_minus_spearman"] = 1.0 - v_err["corr_spearman"]

    inits = _init_order(sorted(set(ll_m["init"]) | set(v_err["init"])))
    colors = {init: _PALETTE[i % len(_PALETTE)] for i, init in enumerate(inits)}

    rows = [
        ("loglik gap to best", ll_m, "loglik_gap", floor),
        ("mean |v - v*|", v_err, "mae", floor),
        ("max |v - v*|", v_err, "maxae", floor),
        ("1 - corr (pearson)", v_err, "one_minus_pearson", floor),
        ("1 - corr (spearman)", v_err, "one_minus_spearman", floor),
        ("|omega_d2 - omega*|", om_m, "omega_gap", floor),
        ("rel. dloglik", ll_m, "rel_change", None),
    ]

    fig, axes = plt.subplots(len(rows), 2, figsize=(11, 17),
                             sharex="col", sharey="row", constrained_layout=True)
    last = len(rows) - 1
    for c, solver in enumerate(SOLVERS):
        axes[0, c].set_title(solver.upper(), fontsize=12, fontweight="bold")
        for r, (label, frame, col, fl) in enumerate(rows):
            ax = axes[r, c]
            _plot_metric(ax, frame[frame.solver == solver],
                         value=col, inits=inits, colors=colors, floor=fl)
            if c == 0:
                ax.set_ylabel(label, fontsize=10)
            ax.grid(True, which="both", alpha=0.25, lw=0.4)
        axes[last, c].axhline(tol, ls="--", color="crimson", lw=1.0)
        axes[last, c].set_xlabel("outer iteration")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper right", ncol=2, fontsize=8,
               frameon=True)
    nu_s = "" if nu is None else (f"   nu={'inf' if not np.isfinite(nu) else f'{nu:g}'}")
    fig.suptitle(f"{slug} - {stage}{nu_s}   |   ref = {ref.solver}/{ref.init} "
                 f"(loglik={ref.best_loglik:.4f})", fontsize=12)
    return fig


def run_slice(slug, stage, *, tol, floor, out_name="convergence.png") -> bool:
    if not has_traces(slug, stage):
        print(f"  [skip] {slug}: no {stage} traces", flush=True)
        return False
    d = _stage_dir(slug, stage)
    ll = pd.read_csv(d / "loglik_traces.csv")
    v = pd.read_parquet(d / "v_traces.parquet")
    ref = reference(ll, v)
    fig = build_figure(slug, stage, loglik_metric_frame(ll, ref.best_loglik),
                       v_error_frame(v, ref), omega_metric_frame(ll, ref), ref,
                       tol=tol, floor=floor, nu=_stage_nu(slug, stage))
    out = d / out_name
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}", flush=True)
    return True


def main(stage: str, *, out_name: str = "convergence.png") -> None:
    """Argument-parsing entry point used by the thin per-stage wrappers."""
    ap = argparse.ArgumentParser(description=f"convergence plots for stage={stage}")
    ap.add_argument("--slices", nargs="+", default=None,
                    help=f"slice slug(s) under results/convergence/. "
                         f"Default: every slice with {stage} traces.")
    ap.add_argument("--tol", type=float, default=_DEFAULT_TOL,
                    help="stopping tolerance marked on the rel-change row.")
    ap.add_argument("--floor", type=float, default=_FLOOR,
                    help="log-axis floor for gap/error/1-corr rows.")
    args = ap.parse_args()

    slugs = args.slices if args.slices else discover_slices(stage)
    if not slugs:
        print(f"No slices with {stage} traces under {OUT_ROOT}")
        return
    n_ok = sum(run_slice(s, stage, tol=args.tol, floor=args.floor, out_name=out_name)
               for s in slugs)
    print(f"\nDone: {n_ok}/{len(slugs)} slice figure(s) plotted.")
