"""Visualize the frozen-omega_d2 sweep from ``e01_omega_profile.py``.

One 2x3 figure per slice:

    A  logML vs omega (peak) with data_loglik twin (monotone)   <- the headline
    B  aging-curve fan colored by omega (visual insensitivity)
    C  estimand drift from the free fit (aging / v_j / gamma)  vs omega
    D  effective AIC/BIC + edf_d vs omega
    E  d-distribution summaries (sd / frac_credible / corr(d,u)) vs omega
    F  runtime + iteration count vs omega

The vertical line / shaded band marks the EB-learned omega* and the +-3x
neighborhood. Reads ``profile_{nutag}.parquet`` + ``curves_{nutag}.parquet``
under ``results/model_selection/athlete_drift/omega_profile/{slug}/``. With no
``--slice`` it renders every slug that has a profile for the given nu.

Run::

    python scripts/02_model_selection/athlete_drift/p01_omega_profile.py
    python scripts/02_model_selection/athlete_drift/p01_omega_profile.py --nu 8 --slice Po10_M_14-25_mrc2
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
from matplotlib import cm  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import registry  # noqa: E402
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "athlete_drift" / "omega_profile"


def _omega_star(prof: pd.DataFrame) -> float:
    free = prof[prof["is_free"]]
    return float(free["omega_d2"].iloc[0]) if len(free) else float(
        prof.loc[prof["logML"].idxmax(), "omega_d2"])


def _band(ax, om_star: float, lo: float = 1 / 3, hi: float = 3.0) -> None:
    ax.axvspan(om_star * lo, om_star * hi, color="0.85", alpha=0.5, zorder=0)
    ax.axvline(om_star, color="k", lw=1.0, ls="--", zorder=1)


def render(slug: str, nutag: str) -> None:
    d = OUT_ROOT / slug
    pf = d / f"profile_{nutag}.parquet"
    if not pf.is_file():
        print(f"  [skip] no profile for {slug} ({nutag})")
        return
    prof = pd.read_parquet(pf).sort_values("omega_d2").reset_index(drop=True)
    curves = (pd.read_parquet(d / f"curves_{nutag}.parquet")
              if (d / f"curves_{nutag}.parquet").is_file() else None)
    om = prof["omega_d2"].to_numpy()
    om_star = _omega_star(prof)

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # --- A: logML peak + data_loglik monotone ------------------------
    a = ax[0, 0]
    a.plot(om, prof["logML"], "-o", color="#4363d8", lw=1.8, label="logML (marginal)")
    a.set_ylabel("logML (marginal)", color="#4363d8")
    a.tick_params(axis="y", labelcolor="#4363d8")
    a2 = a.twinx()
    a2.plot(om, prof["data_loglik"], "-s", color="0.45", lw=1.3, ms=4,
            label="data loglik")
    a2.set_ylabel("data loglik (monotone)", color="0.45")
    a2.tick_params(axis="y", labelcolor="0.45")
    _band(a, om_star)
    a.set_xscale("log")
    a.set_title("A. marginal peaks where EB lands; data-fit is monotone")
    a.set_xlabel("omega_d2 (frozen)")

    # --- B: aging-curve fan ------------------------------------------
    b = ax[0, 1]
    if curves is not None:
        fz = curves[~curves["is_free"]]
        mults = np.sort(fz["omega_mult"].unique())
        norm = Normalize(vmin=np.log10(mults.min()), vmax=np.log10(mults.max()))
        sm = cm.ScalarMappable(norm=norm, cmap="coolwarm")
        for mlt in mults:
            sub = fz[fz["omega_mult"] == mlt].sort_values("A_n")
            b.plot(sub["A_n"], sub["aging"], color=sm.to_rgba(np.log10(mlt)), lw=1.0)
        fr = curves[curves["is_free"]].sort_values("A_n")
        b.plot(fr["A_n"], fr["aging"], color="k", lw=2.2, ls="--", label="free (omega*)")
        cb = fig.colorbar(sm, ax=b, pad=0.01)
        cb.set_label("log10(omega / omega*)")
        b.legend(fontsize=8, loc="best")
    b.set_title("B. aging curve vs omega")
    b.set_xlabel("A_n (years since debut)")
    b.set_ylabel("aging curve (log-time)")

    # --- C: estimand drift from free ---------------------------------
    c = ax[0, 2]
    for col, lab, mk in [("aging_maxdev", "aging max|dev|", "o"),
                         ("max_abs_dv", "v_j max|dev|", "s"),
                         ("gamma_maxdev", "gamma max|dev|", "^")]:
        c.plot(om, prof[col].clip(lower=1e-9), f"-{mk}", lw=1.4, ms=4, label=lab)
    _band(c, om_star)
    c.set_xscale("log"); c.set_yscale("log")
    c.set_title("C. estimand drift from free fit (flat near omega*)")
    c.set_xlabel("omega_d2 (frozen)"); c.set_ylabel("max |deviation| (log-time)")
    c.legend(fontsize=8, loc="best")

    # --- D: ICs + edf_d ----------------------------------------------
    dd = ax[1, 0]
    dd.plot(om, prof["aic_eff"], "-o", color="#e6194b", lw=1.4, ms=4, label="AIC (eff)")
    dd.plot(om, prof["bic_eff"], "-s", color="#f58231", lw=1.4, ms=4, label="BIC (eff)")
    dd.set_ylabel("AIC / BIC (effective dof)")
    dd.legend(fontsize=8, loc="upper left")
    d2 = dd.twinx()
    d2.plot(om, prof["edf_d"], "-^", color="#3cb44b", lw=1.3, ms=4, label="edf_d")
    d2.set_ylabel("edf_d", color="#3cb44b")
    d2.tick_params(axis="y", labelcolor="#3cb44b")
    _band(dd, om_star)
    dd.set_xscale("log")
    dd.set_title("D. effective AIC/BIC + drift edf")
    dd.set_xlabel("omega_d2 (frozen)")

    # --- E: d distribution -------------------------------------------
    e = ax[1, 1]
    e.plot(om, prof["d_sd"], "-o", color="#911eb4", lw=1.4, ms=4, label="sd(d)")
    e.set_ylabel("sd(d_i)", color="#911eb4")
    e.tick_params(axis="y", labelcolor="#911eb4")
    e2 = e.twinx()
    e2.plot(om, prof["frac_credible"], "-s", color="#46f0f0", lw=1.4, ms=4,
            label="frac credible |z|>2")
    e2.plot(om, prof["corr_du"], "-^", color="#000075", lw=1.4, ms=4, label="corr(d,u)")
    e2.set_ylabel("frac credible / corr(d,u)")
    _band(e, om_star)
    e.set_xscale("log")
    e.set_title("E. drift distribution vs omega")
    e.set_xlabel("omega_d2 (frozen)")
    lines = [ln for ln in (e.get_lines() + e2.get_lines())
             if not ln.get_label().startswith("_")]
    e.legend(lines, [ln.get_label() for ln in lines], fontsize=8, loc="best")

    # --- F: runtime + iters ------------------------------------------
    f = ax[1, 2]
    f.plot(om, prof["wall_s"], "-o", color="#808000", lw=1.4, ms=4, label="wall_s")
    f.set_ylabel("wall_s", color="#808000")
    f.tick_params(axis="y", labelcolor="#808000")
    f2 = f.twinx()
    f2.plot(om, prof["n_iter"], "-s", color="0.4", lw=1.3, ms=4, label="n_iter")
    f2.set_ylabel("n_iter", color="0.4")
    _band(f, om_star)
    f.set_xscale("log")
    f.set_title("F. cost vs omega")
    f.set_xlabel("omega_d2 (frozen)")

    nrow = prof[prof["is_free"]].iloc[0] if prof["is_free"].any() else prof.iloc[0]
    fig.suptitle(
        f"{slug}  nu={nutag}  |  omega*={om_star:.3e}  edf_d*={nrow['edf_d']:.0f}  "
        f"n_elig={int(nrow['n_elig'])}  |  shaded = omega* x[1/3, 3]",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = d / f"omega_profile_{nutag}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default=None, help="slug (e.g. Po10_M_14-25_mrc2); "
                    "omit to render all with a profile.")
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()
    nutag = registry.nu_tag(float(args.nu))

    if args.slice:
        render(args.slice, nutag)
        return
    slugs = sorted(p.parent.name for p in OUT_ROOT.glob(f"*/profile_{nutag}.parquet"))
    if not slugs:
        print(f"no profiles found under {OUT_ROOT} for {nutag}")
        return
    print(f"rendering {len(slugs)} slice(s): {slugs}")
    for s in slugs:
        render(s, nutag)


if __name__ == "__main__":
    main()
