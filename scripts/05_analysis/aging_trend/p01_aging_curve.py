"""Entry-age aging-curve fan, one figure per population slice (ALL | Po10).

Two panels -- men (left) | women (right) -- each overlaying the fitted aging
curve f(A_n) for a handful of fixed *entry ages* (default 35/45/55/65), with a
bootstrap uncertainty band per curve. Reads the production AxD fit
(`full_nu8p00_best`) + its athlete-weight bootstrap (`bootstrap/global_coeffs.parquet`)
-- no refit. The production mrc2 fits now carry this bootstrap, so the default
(`--mrcs mrc2`) is the full-field production curve with bands; pass `--mrcs mrc5`
(or both) to overlay the dedicated-runner (>=5 finishes) curve.

Curve math (mirrors the archived `cohort_comparison/aging/aging_common.py`):

    f(A_n; A_e) = theta @ B(A_n)  +  c * A_n  +  (A_e - mean_Ae) * (gamma @ B(A_n))

  * `theta @ B`     spline aging block, anchored f(0)=0 at debut.
  * `c * A_n`       APC gauge tilt. The aging block's *linear* slope is gauge-
                    dependent (career age A_n = period - cohort); production uses
                    the **beta=0** gauge `c = slope(v ~ t_race)`. The tilt is the
                    same for every entry age, so it shifts all four curves
                    together -- the *fan* (spread by entry age) is gauge-invariant.
  * gamma fan       per-unit-centered-entry-age modulation; older debutants
                    decline faster (the entry-age x aging interaction).

x-axis is **career age** A_n (years since debut), the model-native variable, so
all curves share an axis and the fan is read directly. y is the aging block in
log-time (x100 ~ % slowdown vs the athlete's own debut). A grey A_n-density rug
marks where the data is (tail wiggle in the sparse high-A_n zone is extrapolation).

Style (linewidth / fontsize / palette) is deliberately provisional -- a shared
figure style will be applied later.

Outputs (png + vector pdf):
    results/analysis/aging_trend/{slice}_aging_curve.{png,pdf}

Run::

    python scripts/05_analysis/aging_trend/p01_aging_curve.py            # ALL + Po10
    python scripts/05_analysis/aging_trend/p01_aging_curve.py --cohort ALL
    python scripts/05_analysis/aging_trend/p01_aging_curve.py --entry-ages 35 50 65
"""
from __future__ import annotations

import argparse
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import load_slice  # noqa: E402
from marathon_decomp.aging import ncs_basis  # noqa: E402
from marathon_decomp.config import PAPER_FIG_DIR, RESULTS_DIR  # noqa: E402

MODELS_ROOT = RESULTS_DIR / "models"
OUT_ROOT = RESULTS_DIR / "analysis" / "aging_trend"

ENTRY_AGE_COLORS = {35: "#2c7fb8", 45: "#41ab5d", 55: "#fe9929", 65: "#d7301f"}
DEFAULT_ENTRY_AGES = (35, 45, 55, 65)
RUG_HEIGHT_FRAC = 0.25       # A_n-density rug height as a fraction of the y-range


def apply_paper_style() -> None:
    """Load the shared scripts/paper.mplstyle (falls back to the root-relative path)."""
    style_path = Path(__file__).resolve().parents[2] / "paper.mplstyle"
    if style_path.exists():
        plt.style.use(str(style_path))
    else:
        plt.style.use("scripts/paper.mplstyle")


# --------------------------------------------------------------------------- #
# locating the fit + bootstrap (current results/models layout)                #
# --------------------------------------------------------------------------- #
def find_fit_dir(slug: str, model: str, nutag: str) -> Path | None:
    """First `{model}_{nutag}_best__{hash}` dir under results/models/{slug} with a fit.pkl."""
    sd = MODELS_ROOT / slug
    if not sd.is_dir():
        return None
    cands = sorted(sd.glob(f"{model}_{nutag}_best__*"))
    cands = [c for c in cands if (c / "fit.pkl").is_file()]
    return cands[0] if cands else None


# --------------------------------------------------------------------------- #
# curve reconstruction                                                        #
# --------------------------------------------------------------------------- #
def _slope(x: np.ndarray, y: np.ndarray) -> float:
    """OLS slope of y on x, ignoring NaNs in x."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    xm = x - x.mean()
    return float((xm @ (y - y.mean())) / (xm @ xm))


@dataclass
class AgingFan:
    A_grid: np.ndarray          # career age (years since debut)
    raw_point: np.ndarray       # theta-curve at mean entry, point  (G,)
    raw_boot: np.ndarray        # (R, G) replicate theta-curves; (0, G) if point-only
    fan_point: np.ndarray       # gamma-curve per unit centered A_e, point  (G,)
    fan_boot: np.ndarray        # (R, G) replicate gamma-curves; (0, G) if point-only
    mean_Ae: float
    c_beta: float               # production beta=0 tilt
    an_hist: tuple              # (centers, counts) A_n density rug
    knots: np.ndarray           # spline knot locations (career age A_n)
    ci: tuple

    @property
    def has_boot(self) -> bool:
        return self.raw_boot.shape[0] > 0

    def at_entry(self, entry_age: float):
        """Gauged curve for a fixed chronological entry age.

        Returns ``(point, lo, hi)``; ``lo``/``hi`` are ``None`` when no
        bootstrap replicates are available.
        """
        a = entry_age - self.mean_Ae
        gp = self.raw_point + self.c_beta * self.A_grid + a * self.fan_point
        if not self.has_boot:
            return gp, None, None
        gb = self.raw_boot + self.c_beta * self.A_grid[None, :] + a * self.fan_boot
        return gp, np.percentile(gb, self.ci[0], axis=0), np.percentile(gb, self.ci[1], axis=0)


def reconstruct(slug: str, *, model: str, nutag: str, data_version: str,
                a_max: float | None = None, n_grid: int = 200,
                allow_point_only: bool = False,
                ci: tuple[float, float] = (2.5, 97.5)) -> AgingFan | None:
    fit_dir = find_fit_dir(slug, model, nutag)
    if fit_dir is None:
        print(f"  [skip] no point fit for {slug} ({model}_{nutag})")
        return None
    boot = fit_dir / "bootstrap" / "global_coeffs.parquet"
    has_boot = boot.is_file()
    if not has_boot and not allow_point_only:
        print(f"  [skip] no bootstrap for {slug} ({fit_dir.name})")
        return None

    payload = pickle.load(open(fit_dir / "fit.pkl", "rb"))
    knots = np.asarray(payload["model_extra"]["spline_knots"], np.float64)
    v = np.asarray(payload["params"]["v"], np.float64)

    if has_boot:
        gc = pd.read_parquet(boot, columns=["run_id", "block", "k", "value"])

        def _pivot(block: str) -> np.ndarray:  # (R+1, K); row 0 = point fit
            return (gc[gc["block"] == block]
                    .pivot(index="run_id", columns="k", values="value")
                    .sort_index().to_numpy(np.float64))

        theta = _pivot("theta_aging")
        gamma = _pivot("gamma")
    else:
        # point-only: row 0 = point, no replicates -> curve with no band
        theta = np.asarray(payload["params"]["theta_aging"], np.float64)[None, :]
        gamma = np.asarray(payload["params"]["gamma"], np.float64)[None, :]

    fd = load_slice(payload["spec"], payload.get("data_version", data_version))
    A_n = np.asarray(fd.A_n, np.float64)
    A_e = np.asarray(fd.A_e, np.float64)
    mean_Ae = float(np.nanmean(A_e))

    # beta=0 APC gauge: c = slope(v ~ race_year)
    rd = pd.to_datetime(fd.race_date)
    t_year = (rd.year + (rd.dayofyear - 1) / 365.25).to_numpy()
    c_beta = _slope(t_year, v)

    if a_max is None:
        a_max = float(np.quantile(A_n[A_n > 0], 0.975))
    A_grid = np.linspace(0.0, a_max, n_grid)
    B = ncs_basis(A_grid, knots)            # (G, K)
    raw = B @ theta.T                       # (G, R+1), f(0)=0
    fan = B @ gamma.T                       # (G, R+1), per unit centered A_e

    counts, edges = np.histogram(A_n[(A_n >= 0) & (A_n <= a_max)], bins=40)
    centers = 0.5 * (edges[:-1] + edges[1:])

    return AgingFan(
        A_grid=A_grid, raw_point=raw[:, 0], raw_boot=raw[:, 1:].T,
        fan_point=fan[:, 0], fan_boot=fan[:, 1:].T,
        mean_Ae=mean_Ae, c_beta=c_beta, an_hist=(centers, counts),
        knots=knots, ci=ci,
    )


# --------------------------------------------------------------------------- #
# plotting                                                                    #
# --------------------------------------------------------------------------- #
MRC_STYLE = {"mrc2": "-", "mrc5": "--"}  # solid = mrc2, dashed = mrc5


def _mrc_label(tag: str) -> str:
    """'mrc2' -> 'min race count 2'."""
    return f"min race count {tag[3:]}" if tag.startswith("mrc") else tag


def plot_cohort(cohort: str, *, model: str, nutag: str, mrcs,
                data_version: str, entry_ages, scale: float, a_max: float) -> None:
    sexes = [("M", "men"), ("W", "women")]
    # reconstruct every (sex, mrc); point-only allowed (falls back to point curves
    # if a requested fit happens to lack a bootstrap)
    fans: dict[tuple, AgingFan] = {}
    for sx, _ in sexes:
        for mrc in mrcs:
            slug = f"{cohort}_{sx}_14-25_{mrc}"
            print(f"reconstructing {slug} ...")
            fan = reconstruct(slug, model=model, nutag=nutag,
                              data_version=data_version, a_max=a_max,
                              allow_point_only=True)
            if fan is not None:
                fans[(sx, mrc)] = fan
    if not fans:
        print(f"  no fits available for {cohort}; nothing plotted.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(6.1, 3.5), sharex=True, sharey=True,
                             constrained_layout=True)

    # --- curves (color = entry age; solid mrc2 / dashed mrc5) ---------------
    for ax, (sx, _label) in zip(axes, sexes):
        ax.set_xlabel(r"career age  $A_n$  (years since debut)")
        ax.grid(True, color="0.9", lw=0.6, zorder=0)
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--", zorder=1)
        for mrc in mrcs:
            fan = fans.get((sx, mrc))
            if fan is None:
                continue
            ls = MRC_STYLE.get(mrc, "-")
            for ea in entry_ages:
                gp, lo, hi = fan.at_entry(ea)
                color = ENTRY_AGE_COLORS.get(ea)
                if lo is not None:
                    ax.fill_between(fan.A_grid, lo * scale, hi * scale, color=color,
                                    alpha=0.15, lw=0, zorder=2)
                ax.plot(fan.A_grid, gp * scale, color=color, lw=2.0, ls=ls, zorder=3)

    # --- A_n density rug: each panel scaled to ITS OWN peak so the (smaller)
    #     women sample is just as visible as the men's -----------------------
    rug_mrc = next((m for m in mrcs if any((sx, m) in fans for sx, _ in sexes)), mrcs[0])

    def _rug(sx):
        fan = fans.get((sx, rug_mrc)) or next(
            (fans[(sx, m)] for m in mrcs if (sx, m) in fans), None)
        return fan.an_hist if fan is not None else None

    ylo, yhi = axes[0].get_ylim()
    h = RUG_HEIGHT_FRAC * (yhi - ylo)
    for ax, (sx, _label) in zip(axes, sexes):
        r = _rug(sx)
        if r is None:
            continue
        centers, counts = r
        w = (centers[1] - centers[0]) if len(centers) > 1 else 1.0
        peak = counts.max() if counts.max() > 0 else 1.0   # per-panel normalize
        ax.bar(centers, h * counts / peak, width=w, bottom=ylo,
               align="center", color="0.8", zorder=0)
    axes[0].set_ylim(ylo, yhi)

    # --- per-panel subtitles in place of a figure title --------------------
    for ax, subtitle in zip(axes, ("(a) Men", "(b) Women")):
        ax.set_title(subtitle)

    # --- legend: entry ages top->bottom 65..35; mrc line styles only if >1 --
    color_h = [Line2D([0], [0], color=ENTRY_AGE_COLORS.get(ea), lw=2.4,
                      label=f"entry age {int(ea)}")
               for ea in sorted(entry_ages, reverse=True)]
    style_h = ([Line2D([0], [0], color="0.3", lw=2.0, ls=MRC_STYLE.get(m, "-"),
                       label=_mrc_label(m)) for m in mrcs] if len(mrcs) > 1 else [])
    knot_h = [Line2D([0], [0], color="0.55", lw=0.8, ls=(0, (4, 3)),
                     label="spline knot")]
    axes[0].legend(handles=color_h + style_h, fontsize=8, frameon=True,
                   loc="upper left")

    # y = predicted finish-time change vs debut = spline block + APC tilt +
    # entry-age (gamma) interaction. "Aging" would prejudge the sign (runners
    # improve early, decline later), so the label is direction-neutral; the
    # per-curve entry age is what the colour/legend encodes.
    axes[0].set_ylabel("finish-time change relative to debut (log scale)")
    # secondary axis on the right panel: exact % change = (exp(f) - 1) * 100,
    # mapping the plotted left-axis value (scale * f) back through f.
    sec = axes[1].secondary_yaxis(
        "right",
        functions=(lambda y: (np.exp(y / scale) - 1.0) * 100.0,
                   lambda p: scale * np.log1p(np.clip(p, -99.999, None) / 100.0)))
    sec.set_ylabel("% change in finish time, relative to debut")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out = OUT_ROOT / f"{cohort}_aging_curve.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")

    if cohort == 'ALL':
        out_paper = PAPER_FIG_DIR / f"fig_aging_curve_{cohort}.pdf"
        out_paper.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_paper, bbox_inches="tight")
        print(f"wrote {out_paper}")
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default=None, choices=["ALL", "Po10"],
                    help="population slice; default both ALL and Po10.")
    ap.add_argument("--model", default="full", help="registry model tag (AxD = 'full').")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--mrcs", nargs="+", default=["mrc2"],
                    help="min-race-count tags to overlay (solid=mrc2, dashed=mrc5); "
                         "default mrc2 only.")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--entry-ages", type=float, nargs="+", default=list(DEFAULT_ENTRY_AGES))
    ap.add_argument("--a-max", type=float, default=10.0,
                    help="career-age (A_n) grid upper bound; all curves end here.")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="left-axis y multiplier; 1 -> raw log-time units.")
    args = ap.parse_args()

    apply_paper_style()
    cohorts = [args.cohort] if args.cohort else ["ALL", "Po10"]
    for cohort in cohorts:
        plot_cohort(cohort, model=args.model, nutag=args.nutag, mrcs=args.mrcs,
                    data_version=args.data_version, entry_ages=args.entry_ages,
                    scale=args.scale, a_max=args.a_max)


if __name__ == "__main__":
    main()
