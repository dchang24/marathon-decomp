"""Shared convergence-plotter for the rank-1 baseline essential fits.

Both the L2 (Gaussian) and the finite-nu (Student-t) essential fits write the
same per-iteration diagnostic traces (via ``run_essential_fit``), so the
"metric-vs-iteration, ALS | Anderson, one line per init" figure is built here
once and parametrized by the convergence sub-stage:

    stage="L2"          -> results/convergence/{slug}/L2/
    stage="nu_selected" -> results/convergence/{slug}/nu_selected/

The thin runnable wrappers ``03_model_fit/baseline/p01_convergence_l2.py`` and
``p02_convergence_t.py`` just call :func:`main` with their stage.

Figure: a 6-row x 2-column grid (columns = ALS | Anderson) of log-y metrics vs
outer iteration, one coloured line per init (cold ``mean``/``rand*`` and, for the
warm-started t-fit, ``l2_warm``/``l2_rand*``/``step{nu}``). Rows:
  (1) loglik gap to best   (2) mean |v-v*|   (3) max |v-v*|
  (4) 1-corr(v,v*) pearson (5) 1-corr spearman (6) rel. dloglik (+tol line)
The reference v* is the single max-final-loglik run across all solvers x inits;
its own curves rest on the log floor.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "convergence"
SOLVERS = ("als", "anderson")

_EPS = 1e-12
_DEFAULT_TOL = 1e-12
_FLOOR = 1e-10   # log-axis floor for curves that reach 0 (reference / converged)

# Enough distinct, high-contrast colours for the 13-init warm-started t-fit.
_PALETTE = [
    "#e6194b", "#f58231", "#3cb44b", "#4363d8", "#911eb4", "#42d4f4", "#f032e6",
    "#9a6324", "#469990", "#000075", "#808000", "#e6beff", "#aaffc3", "#fabed4",
    "#a9a9a9", "#bfef45",
]


# ── loaders ──────────────────────────────────────────────────────────
def _stage_dir(slug: str, stage: str) -> Path:
    return OUT_ROOT / slug / stage


def has_traces(slug: str, stage: str) -> bool:
    d = _stage_dir(slug, stage)
    return (d / "loglik_traces.csv").exists() and (d / "v_traces.parquet").exists()


def discover_slices(stage: str) -> list[str]:
    if not OUT_ROOT.is_dir():
        return []
    return sorted(d.name for d in OUT_ROOT.iterdir()
                  if d.is_dir() and has_traces(d.name, stage))


def _stage_nu(slug: str, stage: str) -> float | None:
    """The fitted nu, read from convergence_summary.csv (for the figure title)."""
    p = _stage_dir(slug, stage) / "convergence_summary.csv"
    if not p.is_file():
        return None
    df = pd.read_csv(p)
    return float(df["nu"].iloc[0]) if "nu" in df.columns and len(df) else None


# ── reference solution + metric frames ───────────────────────────────
@dataclass
class Reference:
    v_star: pd.Series        # index = race_idx, mean-centered
    best_loglik: float
    solver: str
    init: str


def reference(ll: pd.DataFrame, v: pd.DataFrame) -> Reference:
    finals = ll.sort_values("iter").groupby(["solver", "init"], as_index=False).last()
    best = finals.loc[finals["loglik"].idxmax()]
    solver, init = str(best["solver"]), str(best["init"])
    run = v[(v.solver == solver) & (v.init == init)]
    last_it = int(run["iter"].max())
    v_star = run[run["iter"] == last_it].set_index("race_idx")["v"].sort_index()
    v_star = v_star - v_star.mean()
    return Reference(v_star=v_star, best_loglik=float(best["loglik"]),
                     solver=solver, init=init)


def _ordinal_ranks(M: np.ndarray, axis: int = 1) -> np.ndarray:
    return M.argsort(axis=axis).argsort(axis=axis)


def v_error_frame(v: pd.DataFrame, ref: Reference) -> pd.DataFrame:
    """Per (solver, init, iter): mae, maxae, corr_pearson, corr_spearman vs v*."""
    cols = ref.v_star.index.to_numpy()
    vs_c = ref.v_star.to_numpy() - ref.v_star.to_numpy().mean()
    vs_rank = _ordinal_ranks(vs_c[None, :])[0]
    vs_rank_c = vs_rank - vs_rank.mean()
    vs_rank_den = np.sqrt((vs_rank_c ** 2).sum())

    rows: list[dict] = []
    for (solver, init), g in v.groupby(["solver", "init"], sort=False):
        piv = g.pivot(index="iter", columns="race_idx", values="v").reindex(columns=cols)
        iters = piv.index.to_numpy()
        M = piv.to_numpy()
        M = M - np.nanmean(M, axis=1, keepdims=True)
        err = M - vs_c[None, :]
        mae = np.abs(err).mean(axis=1)
        maxae = np.abs(err).max(axis=1)

        num = (M * vs_c[None, :]).sum(axis=1)
        den = np.sqrt((M ** 2).sum(axis=1) * (vs_c ** 2).sum())
        pearson = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)

        R = _ordinal_ranks(M).astype(float)
        R = R - R.mean(axis=1, keepdims=True)
        snum = (R * vs_rank_c[None, :]).sum(axis=1)
        sden = np.sqrt((R ** 2).sum(axis=1)) * vs_rank_den
        spearman = np.where(sden > 0, snum / np.where(sden > 0, sden, 1.0), np.nan)

        for k, it in enumerate(iters):
            rows.append(dict(solver=solver, init=init, iter=int(it),
                             mae=float(mae[k]), maxae=float(maxae[k]),
                             corr_pearson=float(pearson[k]),
                             corr_spearman=float(spearman[k])))
    return pd.DataFrame(rows)


def loglik_metric_frame(ll: pd.DataFrame, best_loglik: float) -> pd.DataFrame:
    out = []
    for (solver, init), g in ll.sort_values("iter").groupby(["solver", "init"], sort=False):
        gg = g.copy()
        gg["loglik_gap"] = best_loglik - gg["loglik"]
        gg["rel_change"] = gg["loglik"].diff().abs() / (gg["loglik"].abs() + _EPS)
        out.append(gg[["solver", "init", "iter", "loglik", "loglik_gap", "rel_change"]])
    return pd.concat(out, ignore_index=True)


# ── plotting ─────────────────────────────────────────────────────────
def _trail_int(s: str) -> int:
    digits = "".join(ch for ch in s if ch.isdigit())
    return int(digits) if digits else 0


def _init_order(inits: list[str]) -> list[str]:
    """Group cold (mean, rand*) then warm (l2_warm, l2_rand*, step*) inits."""
    def key(s: str) -> tuple[int, int]:
        if s == "mean":
            return (0, -1)
        if s.startswith("l2_rand"):
            return (3, _trail_int(s))
        if s.startswith("l2"):           # l2_warm
            return (2, -1)
        if s.startswith("step"):
            return (4, _trail_int(s))
        if s.startswith("rand"):
            return (1, _trail_int(s))
        return (5, _trail_int(s))
    return sorted(inits, key=key)


def _plot_metric(ax, df, *, value, inits, colors, floor):
    for init in inits:
        g = df[df.init == init].sort_values("iter")
        if g.empty:
            continue
        y = g[value].to_numpy(dtype=float)
        y = np.where(y > 0, y, np.nan) if floor is None else np.maximum(y, floor)
        ax.plot(g["iter"].to_numpy(), y, lw=1.2, color=colors[init], label=init)
    ax.set_yscale("log")


def build_figure(slug, stage, ll_m, v_err, ref, *, tol, floor, nu=None):
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
        ("rel. dloglik", ll_m, "rel_change", None),
    ]

    fig, axes = plt.subplots(len(rows), 2, figsize=(11, 15),
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
    nu_s = "" if nu is None else (f"   ν={'∞' if not np.isfinite(nu) else f'{nu:g}'}")
    fig.suptitle(f"{slug} — {stage}{nu_s}   |   ref = {ref.solver}/{ref.init} "
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
                       v_error_frame(v, ref), ref, tol=tol, floor=floor,
                       nu=_stage_nu(slug, stage))
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
