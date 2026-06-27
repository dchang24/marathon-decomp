"""Residual QQ + tail/skew diagnostics: Gaussian (L2) vs Student-t noise (QC).

The nu sweep (``e01_nu_cv.py``) shows held-out density prefers heavy tails with an
interior optimum at nu~8; this script shows *why* directly, on the fitted
residuals -- the motivation half of the L2-vs-Student-t story for the paper.

For the headline slice it loads the already-registered fits (NO refitting) and,
for each, forms the standardized residual ``z = r / sqrt(sigma2)`` -- the pivot
that is standard-Normal under L2 and standard-t_nu under Student-t -- then:

  * QQ figure (one panel per fit): empirical z-quantiles vs the model's own
    reference law. Baseline L2 is read against the Normal (fat tails peel off);
    baseline t8 and the production full t8 are read against t_nu (straight through
    the bulk + left tail, with the RIGHT tail bending -- the residual skew that
    Student-t, being symmetric, does not capture; cf. Griffin skew-t).
  * Diagnostics table: sigma2 (the L2 variance-inflation number), skewness,
    excess kurtosis, the kurtosis-implied nu (6/excess_kurt + 4) as an
    independent corroboration of the CV nu, and observed-vs-expected tail
    fractions |z|>{2,3,4} under both the Normal and the t_nu reference, split by
    sign to expose the right-skew.

Models read (per slice, by stem; hash globbed): ``baseline_L2_best`` (L2),
``baseline_nu8p00_best`` (t8), ``full_nu8p00_best`` (production AxD, t8). A
``full_L2_best`` / ``full_Linf_best`` is picked up automatically if ever minted,
adding the production-model L2 overlay; absent, the baseline pair carries the
explicit L2-vs-t contrast (consistent with nu being selected on the baseline).

Outputs -> results/model_selection/baseline/ :
  * residual_qq.png          : QQ panels (one per loaded fit).
  * residual_diagnostics.csv : one row per loaded fit, the shape numbers above.

Self-contained; defaults to ALL_B (VS Code "Run" works).

Run::

    python scripts/02_model_selection/baseline/q02_residual_qq.py
    python scripts/02_model_selection/baseline/q02_residual_qq.py --slice Po10_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import load_slice, registry  # noqa: E402
from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from baseline_common import slices as S  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "baseline"
MODELS_ROOT = RESULTS_DIR / "models"

# (stem prefix, short label, panel reference law) in display order. The hash
# suffix is globbed. "normal"/"t" picks the QQ reference; we always also draw the
# *other* law as a faint line so the contrast is visible in every panel.
SPECS = [
    ("baseline_L2_best", "baseline (u+v), L2", "normal"),
    ("baseline_nu8p00_best", "baseline (u+v), t8", "t"),
    ("full_L2_best", "production (AxD), L2", "normal"),
    ("full_Linf_best", "production (AxD), L2", "normal"),
    ("full_nu8p00_best", "production (AxD), t8", "t"),
]

TAIL_KS = (2.0, 3.0, 4.0)

# ── display tunables (edit these) ─────────────────────────────────────
# Cells with |residual| (log-time) above TRIM_ABS_R are non-physical: |r|>1.5
# is a >4.5x time error, beyond any real marathon (genuine hit-the-wall tops out
# ~2.5x, r~0.93). The production aging/d_i block extrapolates a few athlete-races
# to absurd predictions (r up to ~8 => >1000x); those would dominate the moments
# and the axis without telling us anything about the noise law. We exclude them
# from the PLOT and the headline ("physical") moments, but keep the raw all-cell
# moments + n_excluded in the CSV. (Student-t already neutralizes these cells in
# the fit via near-zero IRLS weight, so v_j/sigma2 are unaffected -- the trim is
# only for an honest residual-shape read.)
TRIM_ABS_R = 1.5
X_TICK_STEP = 5.0     # x major gridline spacing (both laws)
X_MINOR_STEP = 1.0    # x minor tick + gridline spacing
Y_TICK_STEP = 10.0    # y gridline spacing
OUT_SUBDIR = "qq_plot"

# ── per-slice axis control (edit per slice) ───────────────────────────
# Within one figure all panels share `ylim`; the L2 (Normal-ref) panel uses
# `xlim_normal`, the t-ref panels use `xlim_t` (the laws' quantiles live on
# different scales). `inset` = (x0, x1, y0, y1) zoom window for the top-left
# inset that magnifies the central bulk (x in the panel's own law-quantile
# units, y in standardized-residual z). Anything omitted falls back to DEFAULTS.
DEFAULTS = dict(
    ylim=(-20.0, 30.0),
    xlim_normal=(-5.0, 5.0),
    xlim_t=(-13.0, 13.0),
    inset=(-2.5, 2.5, -4.0, 4.0),
)
PER_SLICE: dict[str, dict] = {
    "ALL_B":  dict(ylim=(-20.0, 30.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-13.0, 13.0), inset=(-2.5, 2.5, -4.0, 4.0)),
    "ALL_M":  dict(ylim=(-20.0, 30.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-13.0, 13.0), inset=(-2.5, 2.5, -4.0, 4.0)),
    "ALL_W":  dict(ylim=(-20.0, 30.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-13.0, 13.0), inset=(-2.5, 2.5, -4.0, 4.0)),

    "Po10_B": dict(ylim=(-20.0, 22.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-10.0, 10.0), inset=(-2.5, 2.5, -3.0, 4.0)),
    "Po10_M": dict(ylim=(-20.0, 22.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-10.0, 10.0), inset=(-2.5, 2.5, -3.0, 4.0)),
    "Po10_W": dict(ylim=(-20.0, 20.0), xlim_normal=(-5.0, 5.0),
                   xlim_t=(-10.0, 10.0), inset=(-2.5, 2.5, -3.0, 4.0)),
}
SLICES = tuple(PER_SLICE)   # one figure each -> OUT_SUBDIR/qq_{slug}.png
# ──────────────────────────────────────────────────────────────────────


def _settings(name: str) -> dict:
    s = dict(DEFAULTS)
    s.update(PER_SLICE.get(name, {}))
    return s


def _find_fit_dir(slug: str, stem: str) -> Path | None:
    hits = sorted((MODELS_ROOT / slug).glob(f"{stem}__*"))
    return hits[0] if hits else None


def _nu_of(m) -> float:
    nu = m.params.get("nu", None)
    if nu is None:
        nu = getattr(m.config, "nu", float("inf"))
    return float(nu)


def _qq_grid(n: int) -> np.ndarray:
    """Probability grid for QQ points, denser in the tails (where the action is)."""
    body = np.linspace(0.01, 0.99, 120)
    tail = np.concatenate([
        np.geomspace(1.0 / n, 0.01, 40),          # lower tail down to ~1/n
        1.0 - np.geomspace(1.0 / n, 0.01, 40),    # upper tail
    ])
    return np.unique(np.clip(np.concatenate([body, tail]), 1.0 / n, 1.0 - 1.0 / n))


def _diagnostics(z: np.ndarray, nu: float, sigma2: float) -> dict:
    n = z.size
    skew = float(stats.skew(z))
    exk = float(stats.kurtosis(z, fisher=True))     # excess kurtosis
    # kurtosis-implied nu: for t_nu, excess kurtosis = 6/(nu-4) (nu>4) -> nu = 6/exk + 4
    nu_kurt = float(6.0 / exk + 4.0) if exk > 1e-9 else float("inf")
    out = dict(nu=nu, n=n, sigma2=sigma2, sigma=float(np.sqrt(sigma2)),
               skewness=skew, excess_kurtosis=exk, nu_implied_by_kurtosis=nu_kurt)
    ref_nu = stats.t(df=nu) if np.isfinite(nu) else stats.norm
    for k in TAIL_KS:
        obs_pos = float(np.mean(z > k))
        obs_neg = float(np.mean(z < -k))
        out[f"obs_frac_gt{k:g}"] = obs_pos + obs_neg
        out[f"obs_frac_pos_gt{k:g}"] = obs_pos      # right tail (hit-the-wall)
        out[f"obs_frac_neg_lt{k:g}"] = obs_neg      # left tail
        out[f"exp_frac_normal_gt{k:g}"] = float(2.0 * stats.norm.sf(k))
        out[f"exp_frac_ref_gt{k:g}"] = float(2.0 * ref_nu.sf(k))
    return out


def _qq_r2(z: np.ndarray, p: np.ndarray, dist) -> float:
    """Straightness of the QQ against `dist`: R^2 of empirical vs theoretical q."""
    emp = np.quantile(z, p)
    theo = dist.ppf(p)
    ok = np.isfinite(theo) & np.isfinite(emp)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(emp[ok], theo[ok])[0, 1] ** 2)


def _prep(z: np.ndarray, nu: float) -> dict:
    """Empirical + reference quantiles for one fit, on a shared probability grid."""
    p = _qq_grid(z.size)
    t_dist = stats.t(df=nu) if np.isfinite(nu) else stats.norm
    return dict(p=p, emp=np.quantile(z, p), norm_q=stats.norm.ppf(p),
                t_q=t_dist.ppf(p), t_dist=t_dist,
                r2_normal=_qq_r2(z, p, stats.norm), r2_t=_qq_r2(z, p, t_dist))


def _content(ax, d: dict, nu: float, ref: str, *, ms: float, labels: bool) -> None:
    """Draw the QQ content (points + model-law line + rival-law line) on `ax`.
    Limits must already be set; lines are drawn across the current x-range."""
    emp, norm_q, t_q = d["emp"], d["norm_q"], d["t_q"]
    primary_q = norm_q if ref == "normal" else t_q
    ax.plot(primary_q, emp, ".", ms=ms, color="#1f3b73", zorder=3,
            label="empirical residual z" if labels else None)
    # red solid = this panel's model law (points land on y=x if z follows it).
    lo, hi = ax.get_xlim()
    ax.plot([lo, hi], [lo, hi], "-", color="#e6194b", lw=1.4, zorder=2,
            label=(("Normal" if ref == "normal" else f"t({nu:g})") + " (model law, y=x)")
            if labels else None)
    # gray dashed = the RIVAL law drawn on these same axes: for each p the point
    # (this_law_q(p), rival_law_q(p)) -- where points WOULD lie if z followed the
    # rival law. Peels away from y=x in the tails.
    if ref == "normal" and np.isfinite(nu):
        ax.plot(norm_q, t_q, "--", color="#888", lw=1.1, zorder=1,
                label=f"if z were t({nu:g}) (heavy-tailed)" if labels else None)
    elif ref == "t":
        ax.plot(t_q, norm_q, "--", color="#888", lw=1.1, zorder=1,
                label="if z were Normal (thin-tailed)" if labels else None)


def _panel(ax, d: dict, nu: float, label: str, ref: str, *,
           xlim, ylim, inset, ylabel) -> None:
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_autoscale_on(False)
    _content(ax, d, nu, ref, ms=3, labels=True)
    ax.xaxis.set_major_locator(MultipleLocator(X_TICK_STEP))
    ax.xaxis.set_minor_locator(MultipleLocator(X_MINOR_STEP))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_STEP))
    r2 = (f"R$^2$ N {d['r2_normal']:.3f}"
          + (f" / t {d['r2_t']:.3f}" if np.isfinite(nu) else ""))
    ax.set_title(f"{label}   ({r2})", fontsize=9)
    ax.set_xlabel(("standard Normal" if ref == "normal" else f"t({nu:g})")
                  + " quantiles")
    if ylabel:
        ax.set_ylabel("standardized residual z")
    ax.grid(True, which="major", alpha=0.25, lw=0.4)
    ax.grid(True, which="minor", axis="x", alpha=0.12, lw=0.3)
    ax.legend(fontsize=7, loc="lower right")

    # top-left zoom inset over the central bulk.
    if inset is not None:
        x0, x1, y0, y1 = inset
        axins = ax.inset_axes([0.045, 0.55, 0.40, 0.40])
        axins.set_xlim(x0, x1); axins.set_ylim(y0, y1); axins.set_autoscale_on(False)
        _content(axins, d, nu, ref, ms=2, labels=False)
        axins.grid(True, alpha=0.25, lw=0.3)
        axins.tick_params(labelsize=6)
        ax.indicate_inset_zoom(axins, edgecolor="0.5", lw=0.8, alpha=0.6)


def render_slice(slice_name: str, mrc: int | None, out_dir: Path) -> list[dict]:
    """One QQ figure + the per-fit diagnostic rows for `slice_name`."""
    spec = S.build_spec(slice_name, min_race_count=mrc)
    slug = S.slug(spec)
    st = _settings(slice_name)
    print(f"=== {slug} ===", flush=True)

    # de-dup stems (full_L2 / full_Linf are aliases; keep the first that exists)
    panels: list[dict] = []   # {label, ref, nu, d}
    rows: list[dict] = []
    seen_dirs: set[Path] = set()
    fd = None
    for stem, label, ref in SPECS:
        fdir = _find_fit_dir(slug, stem)
        if fdir is None or fdir in seen_dirs:
            continue
        seen_dirs.add(fdir)
        if fd is None:
            fd = load_slice(spec)
            print(f"    I={fd.I:,} J={fd.J:,} N={fd.N:,}", flush=True)
        m = registry.load_fit(fdir, fd)
        # standardize by the model's own fitted sigma2 (the assumed scale), then
        # split off non-physical cells (|r|>TRIM_ABS_R) for an honest shape read.
        sigma2 = float(m.params["sigma2"])
        nu = _nu_of(m)
        r = m.residuals()
        z_all = r / np.sqrt(sigma2)
        keep = np.abs(r) <= TRIM_ABS_R
        z = z_all[keep]
        n_excl = int((~keep).sum())
        d = _prep(z, nu)
        panels.append(dict(label=label, ref=ref, nu=nu, d=d))
        rows.append(dict(
            slice=slug, model=label,
            **_diagnostics(z, nu, sigma2),               # headline = physical
            n_total=int(z_all.size), n_excluded=n_excl,
            frac_excluded=n_excl / z_all.size,
            trim_abs_r=TRIM_ABS_R, max_abs_z_kept=float(np.abs(z).max()),
            skewness_raw=float(stats.skew(z_all)),
            excess_kurtosis_raw=float(stats.kurtosis(z_all, fisher=True)),
            qq_r2_normal=d["r2_normal"], qq_r2_ref=d["r2_t"]))
        print(f"    {stem:24s} nu={nu:<4g} trimmed={n_excl:<4d} "
              f"QQ-R2(t)={d['r2_t']:.4f}", flush=True)

    if not panels:
        print(f"    no registered fits under {MODELS_ROOT / slug} -- skipped.")
        return rows

    ncol = len(panels)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 4.0),
                             squeeze=False, sharey=True)
    for c, pan in enumerate(panels):
        xlim = st["xlim_normal"] if pan["ref"] == "normal" else st["xlim_t"]
        _panel(axes[0, c], pan["d"], pan["nu"], pan["label"], pan["ref"],
               xlim=xlim, ylim=st["ylim"], inset=st["inset"], ylabel=(c == 0))
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"qq_{slug}.png"
    pdf = out_dir / f"qq_{slug}.pdf"
    fig.savefig(png, dpi=130, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {png}", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=list(SLICES))
    ap.add_argument("--mrc", type=int, default=None)
    args = ap.parse_args()

    out_dir = OUT_ROOT / OUT_SUBDIR
    all_rows: list[dict] = []
    for name in args.slices:
        all_rows.extend(render_slice(name, args.mrc, out_dir))

    if not all_rows:
        print("No fits found for any requested slice.")
        return

    df = pd.DataFrame(all_rows)
    csv = out_dir / "residual_diagnostics.csv"
    df.to_csv(csv, index=False)

    # ── console summary (ASCII only) ─────────────────────────────────
    print("\n" + "=" * 78)
    show = ["slice", "model", "nu", "sigma2", "n_excluded", "skewness",
            "excess_kurtosis", "nu_implied_by_kurtosis", "qq_r2_normal", "qq_r2_ref"]
    with pd.option_context("display.width", 220, "display.max_columns", 25):
        print(df[show].to_string(index=False))
    print("=" * 78)
    print(f"wrote {csv}")


if __name__ == "__main__":
    main()
