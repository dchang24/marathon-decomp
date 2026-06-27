"""Visualize the omega_d2 INIT sensitivity from ``e02_omega_init.py``.

One 2x2 figure per slice:

    A  omega_d2 trajectory vs iter (every init -> same omega*)   <- headline
    B  converged omega* and n_iter vs init (invariance + cost)
    C  estimand deviation from reference vs init (init-invariant results)
    D  loglik trajectory vs iter (all reach the same plateau)

Vertical/horizontal guides mark omega* and the production default init. Reads
``init_summary_{nutag}.parquet`` + ``init_traces_{nutag}.parquet`` under
``results/model_selection/athlete_drift/omega_init/{slug}/``. With no ``--slice``
renders every slug with a summary for the given nu.

Run::

    python scripts/02_model_selection/athlete_drift/p02_omega_init.py
    python scripts/02_model_selection/athlete_drift/p02_omega_init.py --slice Po10_M_14-25_mrc2 --nu 8
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

OUT_ROOT = RESULTS_DIR / "model_selection" / "athlete_drift" / "omega_init"


def _color_map(summ: pd.DataFrame):
    """label -> color; numeric inits on a viridis ramp, default in black."""
    num = summ[~summ["is_reference"]].copy()
    vals = num["omega_d2_init"].to_numpy()
    norm = Normalize(vmin=np.log10(np.nanmin(vals)), vmax=np.log10(np.nanmax(vals)))
    sm = cm.ScalarMappable(norm=norm, cmap="viridis")
    colors = {r["init_label"]: sm.to_rgba(np.log10(r["omega_d2_init"]))
              for _, r in num.iterrows()}
    for _, r in summ[summ["is_reference"]].iterrows():
        colors[r["init_label"]] = "k"
    return colors, sm


def render(slug: str, nutag: str) -> None:
    d = OUT_ROOT / slug
    sf = d / f"init_summary_{nutag}.parquet"
    if not sf.is_file():
        print(f"  [skip] no init summary for {slug} ({nutag})")
        return
    summ = pd.read_parquet(sf)
    traces = pd.read_parquet(d / f"init_traces_{nutag}.parquet")
    ref = summ[summ["is_reference"]].iloc[0]
    om_star = float(ref["omega_d2"])
    # resolved default init: the reference run's earliest recorded omega_d2 (just
    # after iter 0) -- robust to whatever the model's default_init resolves to.
    ref_tr = traces[traces["init_label"] == ref["init_label"]].sort_values("iter")
    default_init = float(ref_tr["omega_d2"].iloc[0]) if len(ref_tr) else om_star
    colors, sm = _color_map(summ)

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))

    # --- A: omega_d2 trajectory --------------------------------------
    a = ax[0, 0]
    for lab, g in traces.groupby("init_label"):
        g = g.sort_values("iter")
        a.plot(g["iter"], g["omega_d2"].clip(lower=1e-12),
               color=colors.get(lab, "0.5"),
               lw=2.2 if colors.get(lab) == "k" else 1.2,
               ls="--" if colors.get(lab) == "k" else "-",
               label=lab if colors.get(lab) == "k" else None)
    a.axhline(om_star, color="r", lw=0.9, ls=":", label=f"omega*={om_star:.2e}")
    a.set_yscale("log")
    a.set_xlabel("outer iteration"); a.set_ylabel("omega_d2")
    a.set_title("A. omega_d2 path by init (all -> omega*; tiny init stalls)")
    a.legend(fontsize=8, loc="best")
    cb = fig.colorbar(sm, ax=a, pad=0.01)
    cb.set_label("log10(omega_d2_init)  [default=black]")

    # --- B: converged omega* + n_iter vs init ------------------------
    b = ax[0, 1]
    num = summ[~summ["is_reference"]].sort_values("omega_d2_init")
    b.plot(num["omega_d2_init"], num["omega_d2"], "-o", color="#4363d8",
           label="converged omega*")
    b.axhline(om_star, color="r", lw=0.9, ls=":")
    b.axvline(default_init, color="0.5", lw=0.9, ls="--", label="default init")
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlabel("omega_d2_init"); b.set_ylabel("converged omega*", color="#4363d8")
    b.tick_params(axis="y", labelcolor="#4363d8")
    b2 = b.twinx()
    b2.plot(num["omega_d2_init"], num["n_iter"], "-s", color="#e6194b", label="n_iter")
    b2.set_ylabel("n_iter (cost)", color="#e6194b")
    b2.tick_params(axis="y", labelcolor="#e6194b")
    b.set_title("B. converged omega* flat; iteration cost is U-shaped")
    b.legend(fontsize=8, loc="lower right")

    # --- C: estimand deviation from reference vs init ----------------
    c = ax[1, 0]
    for col, lab, mk in [("aging_maxdev", "aging max|dev|", "o"),
                         ("max_abs_dv", "v_j max|dev|", "s"),
                         ("gamma_maxdev", "gamma max|dev|", "^")]:
        c.plot(num["omega_d2_init"], num[col].clip(lower=1e-12), f"-{mk}",
               lw=1.4, ms=4, label=lab)
    c.axvline(default_init, color="0.5", lw=0.9, ls="--")
    c.set_xscale("log"); c.set_yscale("log")
    c.set_xlabel("omega_d2_init"); c.set_ylabel("max |deviation| from reference")
    c.set_title("C. results init-invariant (dev ~ 0 where converged)")
    c.legend(fontsize=8, loc="best")

    # --- D: loglik trajectory ----------------------------------------
    dd = ax[1, 1]
    for lab, g in traces.groupby("init_label"):
        g = g.sort_values("iter")
        dd.plot(g["iter"], g["loglik"], color=colors.get(lab, "0.5"),
                lw=2.2 if colors.get(lab) == "k" else 1.2,
                ls="--" if colors.get(lab) == "k" else "-")
    dd.set_xlabel("outer iteration"); dd.set_ylabel("penalized loglik")
    dd.set_title("D. loglik path by init (same plateau)")

    fig.suptitle(
        f"{slug}  nu={nutag}  |  omega_d2 INIT sensitivity  |  omega*={om_star:.3e}  "
        f"default init={default_init:.2e}  |  converged spread "
        f"{summ['omega_d2'].max()/summ['omega_d2'].min() - 1:.1e}",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = d / f"omega_init_{nutag}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default=None, help="slug; omit to render all.")
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()
    nutag = registry.nu_tag(float(args.nu))

    if args.slice:
        render(args.slice, nutag)
        return
    slugs = sorted(p.parent.name for p in OUT_ROOT.glob(f"*/init_summary_{nutag}.parquet"))
    if not slugs:
        print(f"no init summaries under {OUT_ROOT} for {nutag}")
        return
    print(f"rendering {len(slugs)} slice(s): {slugs}")
    for s in slugs:
        render(s, nutag)


if __name__ == "__main__":
    main()
