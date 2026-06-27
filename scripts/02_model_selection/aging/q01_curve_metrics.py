"""Quantify aging-curve shape -- recomputed from persisted fits, no re-fit.

Reads the e01 ``grid/fits/*.pkl`` payloads (the source of truth) and, for each
fitted cell, reconstructs the curve at the mean entry age plus a few entry-age
percentiles (via ``aging_curve_from_payload`` / ``entry_age_curve_from_payload``,
which need only the pickle dict). Applies ``curve_metrics.curve_shape_metrics`` to
each curve: peak age/depth, improvement span, break-even age, plateau widths, and
the tail-wiggle flag.

This is a cheap QC pass -- edit it freely to add metrics without touching e01.

Output -> results/model_selection/aging/{slug}/curve_metrics.csv
       + cross-slice rollup curve_metrics_all.csv at the dir root.

Run::

    python scripts/02_model_selection/aging/q01_curve_metrics.py            # all slices with e01 output
    python scripts/02_model_selection/aging/q01_curve_metrics.py --slices ALL_M Po10_M
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/
sys.path.insert(0, str(Path(__file__).resolve().parent))       # this dir

from marathon_decomp import (  # noqa: E402
    aging_curve_from_payload,
    entry_age_curve_from_payload,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402
import curve_metrics as CM  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
N_GRID = 200            # MUST match e01_aging_grid.N_GRID


def _cand_from_config(cfg) -> tuple[str, str, str]:
    """(basis_name, gamma_token, cand) from a saved ModelConfig."""
    basis = (f"poly{cfg.degree}" if cfg.basis_kind == "poly"
             else f"spline{cfg.n_knots}")
    token = cfg.gamma_form if cfg.use_gamma else "off"
    return basis, token, f"{basis}-g{token}"


def _nu_label(nu: float) -> str:
    return "inf" if not np.isfinite(float(nu)) else f"{float(nu):g}"


def run_slice(slug: str) -> pd.DataFrame | None:
    grid_dir = OUT_ROOT / slug / "grid"
    fits = sorted((grid_dir / "fits").glob("*.pkl"))
    if not fits:
        print(f"  [skip] no fits under {grid_dir / 'fits'}")
        return None
    info = pd.read_csv(grid_dir / "slice_info.csv").iloc[0]
    ae_mean, an_max, an_p95 = float(info.ae_mean), float(info.an_max), float(info.an_p95)
    A_grid = np.linspace(0.0, an_max, N_GRID)

    # entry ages to evaluate: mean (gamma contributes 0) + p10/p50/p90.
    ages = {"mean": ae_mean, "p10": float(info.ae_p10),
            "p50": float(info.ae_p50), "p90": float(info.ae_p90)}

    rows: list[dict] = []
    for pkl in fits:
        solver = pkl.stem.split("_")[-1]
        with pkl.open("rb") as f:
            payload = pickle.load(f)
        cfg = payload["config"]
        basis, token, cand = _cand_from_config(cfg)
        aging = aging_curve_from_payload(payload, A_grid)
        for age_label, age in ages.items():
            if not np.isfinite(age):
                continue
            g = aging.copy()
            if cfg.use_gamma and age_label != "mean":
                g = g + entry_age_curve_from_payload(payload, A_grid, age - ae_mean)
            m = CM.curve_shape_metrics(g, A_grid, entry_age=age, an_p95=an_p95)
            rows.append(dict(
                slug=slug, nu=_nu_label(cfg.nu), cand=cand, basis=basis,
                gamma_form=token, solver=solver, entry_age_label=age_label,
                entry_age=age, **m))
    df = pd.DataFrame(rows)
    df.to_csv(OUT_ROOT / slug / "curve_metrics.csv", index=False)
    print(f"  wrote {OUT_ROOT / slug / 'curve_metrics.csv'}  ({len(df)} rows)")
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"],
                    help="slug names or 'all' (every slice with e01 output).")
    args = ap.parse_args()

    if args.slices == ["all"]:
        slugs = sorted(p.parent.parent.name for p in OUT_ROOT.glob("*/grid/slice_info.csv"))
    else:
        # accept either slice names or full slugs by matching the dir prefix.
        slugs = []
        for s in args.slices:
            hit = [p.parent.parent.name for p in OUT_ROOT.glob(f"{s}*/grid/slice_info.csv")]
            slugs.extend(hit or [s])

    parts = []
    for slug in slugs:
        print(f"=== {slug} ===")
        df = run_slice(slug)
        if df is not None:
            parts.append(df)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(OUT_ROOT / "curve_metrics_all.csv", index=False)
        print(f"\nRollup -> {OUT_ROOT / 'curve_metrics_all.csv'}")


if __name__ == "__main__":
    main()
