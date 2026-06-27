"""Descriptors + density band for the career-drift d_i distribution (no refit).

Per cohort (ALL, Po10) x sex (M, W) x career-persistence floor (n_i >= {3,5,10},
post-hoc subsets of the production mrc2 AxD fit), the fitted d_i sample is reduced
to scalar descriptors:

  * disatt_sd       disattenuated SD = sqrt(Var(d_hat) - mean post_var); the true
                    spread of career trajectories, EB posterior noise removed.
                    GAUGE-INVARIANT -- the headline heterogeneity number.
  * skew            distribution skew (right tail = rapid decliners). GAUGE-INVARIANT.
  * frac_improver   P(d < 0). EB-prior/gauge-anchored (see di_common caveat).
  * corr(d,u)       Spearman (+ high-n_i cut). Prior/gauge-anchored AND era-
                    confounded -- read with care.

This script owns ALL the resampling. Point estimates come straight from your
fitted d_i (NOTHING is refit). Two athlete-bootstrap products, both accumulating:

  1. scalar descriptor CIs   -> bootstrap/draws.csv         (--n-boot)
  2. density-curve band      -> bootstrap/density_draws.parquet + density_band.csv
                                (--n-band; KDE of each resample on a fixed grid)

The plotting script (p01) just READS density_band.csv -- it does no resampling.

ACCUMULATING DRAWS (run as many times as you like, never lose any). ``--n-boot``
and ``--n-band`` are the draws to **ADD** this run, not totals. Each draw k has a
reproducible seed ``SeedSequence([seed, crc32(cohort_sex_floor), k])`` and a run
continues at ``k = max_existing + 1``, so re-runs only grow the pool and CIs are
re-derived over the full pool each time. Point estimates are deterministic.

    q01 --n-boot 1000 --n-band 200      # scalars k=0..999, band k=0..199
    q01 --n-boot 1000 --n-band 200      # -> scalar pool 2000, band pool 400
    q01 --n-boot 0    --n-band 0        # just rebuild summaries from existing pools

Storage: draws.csv ~60 B/row (12 rows/draw). density_draws.parquet holds one row
per cell/draw with ~140 grid columns (compresses well; ~tens of MB at thousands of
band draws -- the band converges in a few hundred). Delete a file to reset it.

Outputs (results/analysis/di_distribution/):
    bootstrap/draws.csv             scalar-descriptor draw accumulator
    bootstrap/density_draws.parquet density-curve draw accumulator
    bootstrap/density_grid.csv      the fixed x-grid (keeps draws aligned across runs)
    di_descriptors.parquet / .csv   scalar summary, point + CI over the full pool
    di_descriptors.md               readable table, CIs inline
    density_band.csv                per (cohort,sex,floor,x): point density + 95% band

Run::

    python scripts/05_analysis/di_distribution/q01_di_descriptors.py --n-boot 1000 --n-band 200
"""
from __future__ import annotations

import argparse
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir (di_common)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

import di_common as DI  # noqa: E402

FLOORS = (3, 5, 10)
KEYS = ["disatt_sd", "skew", "frac_improver", "corr_du", "corr_du_highn"]
BOOT_DIR = DI.OUT_ROOT / "bootstrap"
DRAWS = BOOT_DIR / "draws.csv"
DRAW_COLS = ["cohort", "sex", "n_floor", "k"] + KEYS
DENS_DRAWS = BOOT_DIR / "density_draws.parquet"
GRID_FILE = BOOT_DIR / "density_grid.csv"
BAND_OUT = DI.OUT_ROOT / "density_band.csv"
BAND_GRID_N = 140


# --------------------------------------------------------------------------- #
# small helpers                                                                #
# --------------------------------------------------------------------------- #
def _crc(cohort: str, sex: str, floor: int) -> int:
    return zlib.crc32(f"{cohort}_{sex}_{floor}".encode())


def _ci(a: np.ndarray) -> tuple[float, float]:
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan, np.nan
    return float(np.percentile(a, DI.CI[0])), float(np.percentile(a, DI.CI[1]))


def _kde(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    x = x[np.isfinite(x)]
    if x.size < 5 or x.std() < 1e-9:
        return np.zeros_like(grid)
    return gaussian_kde(x)(grid)


# --------------------------------------------------------------------------- #
# scalar-descriptor accumulator                                                #
# --------------------------------------------------------------------------- #
def _read_draws() -> pd.DataFrame:
    return pd.read_csv(DRAWS) if DRAWS.is_file() else pd.DataFrame(columns=DRAW_COLS)


def _n_existing(df: pd.DataFrame, cohort: str, sex: str, floor: int) -> int:
    if df.empty:
        return 0
    m = (df.cohort == cohort) & (df.sex == sex) & (df.n_floor == floor)
    return int(df.loc[m, "k"].max()) + 1 if m.any() else 0


def _compute_draws(d, u, pv, ni, cohort, sex, floor, k0, n_add, seed) -> list[dict]:
    crc, Ne, rows = _crc(cohort, sex, floor), len(d), []
    for k in range(k0, k0 + n_add):
        idx = np.random.default_rng(np.random.SeedSequence([seed, crc, k])).integers(0, Ne, Ne)
        dd = DI.descriptors(d[idx], u[idx], pv[idx], ni[idx])
        rows.append({"cohort": cohort, "sex": sex, "n_floor": floor, "k": k,
                     **{kk: dd[kk] for kk in KEYS}})
    return rows


# --------------------------------------------------------------------------- #
# density-curve accumulator                                                    #
# --------------------------------------------------------------------------- #
def _dcols(G: int) -> list[str]:
    return [f"d{g:03d}" for g in range(G)]


def _get_grid(fits: dict) -> np.ndarray:
    """Fixed x-grid (%/yr) persisted so accumulated density draws stay aligned."""
    if GRID_FILE.is_file():
        return pd.read_csv(GRID_FILE)["x"].to_numpy(float)
    alld = np.concatenate([f.d[f.elig] * 100.0 for f in fits.values() if f is not None])
    lo, hi = np.percentile(alld, [1, 99])
    pad = 0.1 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, BAND_GRID_N)
    BOOT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"x": grid}).to_csv(GRID_FILE, index=False)
    return grid


def _read_dens() -> pd.DataFrame | None:
    return pd.read_parquet(DENS_DRAWS) if DENS_DRAWS.is_file() else None


def _n_existing_dens(df, cohort, sex, floor) -> int:
    if df is None or df.empty:
        return 0
    m = (df.cohort == cohort) & (df.sex == sex) & (df.n_floor == floor)
    return int(df.loc[m, "k"].max()) + 1 if m.any() else 0


def _compute_dens(d100, grid, cohort, sex, floor, k0, n_add, seed) -> list[dict]:
    crc, n, cols, rows = _crc(cohort, sex, floor), len(d100), _dcols(len(grid)), []
    for k in range(k0, k0 + n_add):
        idx = np.random.default_rng(np.random.SeedSequence([seed, crc, k])).integers(0, n, n)
        dens = _kde(d100[idx], grid)
        rows.append({"cohort": cohort, "sex": sex, "n_floor": floor, "k": k,
                     **dict(zip(cols, dens))})
    return rows


# --------------------------------------------------------------------------- #
# main                                                                         #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--mrc", type=int, default=2, help="base fit; floors are post-hoc n_i subsets.")
    ap.add_argument("--n-boot", type=int, default=1000,
                    help="scalar-descriptor draws to ADD this run (accumulates).")
    ap.add_argument("--n-band", type=int, default=200,
                    help="density-curve draws to ADD this run (accumulates; KDE per draw).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    # load every cell's fit up front (grid needs them all)
    fits = {}
    for cohort in DI.POPS:
        for sx in ("M", "W"):
            fits[(cohort, sx)] = DI.load_cached(cohort, sx, args.nu, args.mrc)
    avail = [c for c in DI.POPS if fits[(c, "M")] is not None and fits[(c, "W")] is not None]
    if not avail:
        print("no fits found; nothing written.")
        return
    grid = _get_grid(fits)

    draws_df = _read_draws()
    dens_df = _read_dens()
    new_scalar, new_dens = [], []
    arrays = {}  # (cohort,sx,floor) -> (d,u,pv,ni,d100)

    for cohort in avail:
        for sx in ("M", "W"):
            for floor in FLOORS:
                t = fits[(cohort, sx)].table(min_n=floor)
                d, u = t.d.to_numpy(), t.u.to_numpy()
                pv, ni = t.post_var.to_numpy(), t.n_i.to_numpy()
                arrays[(cohort, sx, floor)] = (d, u, pv, ni, d * 100.0)
                if args.n_boot > 0:
                    k0 = _n_existing(draws_df, cohort, sx, floor)
                    new_scalar += _compute_draws(d, u, pv, ni, cohort, sx, floor,
                                                 k0, args.n_boot, args.seed)
                if args.n_band > 0:
                    k0b = _n_existing_dens(dens_df, cohort, sx, floor)
                    new_dens += _compute_dens(d * 100.0, grid, cohort, sx, floor,
                                              k0b, args.n_band, args.seed)
                    print(f"  {cohort} {sx} n>={floor}: +{args.n_boot} scalar, "
                          f"+{args.n_band} band draws")

    # persist accumulators
    BOOT_DIR.mkdir(parents=True, exist_ok=True)
    if new_scalar:
        draws_df = pd.concat([draws_df, pd.DataFrame(new_scalar)], ignore_index=True)
        draws_df[["n_floor", "k"]] = draws_df[["n_floor", "k"]].astype(int)
        draws_df[KEYS] = draws_df[KEYS].astype(float)
        draws_df = draws_df.sort_values(["cohort", "sex", "n_floor", "k"]).reset_index(drop=True)
        draws_df.to_csv(DRAWS, index=False)
    if new_dens:
        add = pd.DataFrame(new_dens)
        dens_df = add if dens_df is None else pd.concat([dens_df, add], ignore_index=True)
        dens_df = dens_df.sort_values(["cohort", "sex", "n_floor", "k"]).reset_index(drop=True)
        dens_df.to_parquet(DENS_DRAWS, index=False)

    spool = ({k: g for k, g in draws_df.groupby(["cohort", "sex", "n_floor"])}
             if not draws_df.empty else {})
    dpool = ({k: g for k, g in dens_df.groupby(["cohort", "sex", "n_floor"])}
             if dens_df is not None and not dens_df.empty else {})

    # ---- scalar summary: point (from fit) + CI (over full pool) ----------- #
    def sarr(c, s, f, key):
        g = spool.get((c, s, f))
        return g[key].to_numpy(dtype=float) if g is not None else np.array([], float)

    recs = []
    for cohort in avail:
        for floor in FLOORS:
            pt = {}
            for sx in ("M", "W"):
                d, u, pv, ni, _ = arrays[(cohort, sx, floor)]
                p = DI.descriptors(d, u, pv, ni)
                pt[sx] = p
                row = dict(cohort=cohort, group=sx, n_floor=floor, n=len(d),
                           n_draws=len(sarr(cohort, sx, floor, "disatt_sd")))
                for key in KEYS:
                    row[key] = p[key]
                    row[f"{key}_lo"], row[f"{key}_hi"] = _ci(sarr(cohort, sx, floor, key))
                recs.append(row)
            mrow = dict(cohort=cohort, group="M-W", n_floor=floor, n=np.nan)
            nmin = min(len(sarr(cohort, "M", floor, "disatt_sd")),
                       len(sarr(cohort, "W", floor, "disatt_sd")))
            mrow["n_draws"] = nmin
            for key in KEYS:
                mrow[key] = pt["M"][key] - pt["W"][key]
                diff = (sarr(cohort, "M", floor, key)[:nmin] -
                        sarr(cohort, "W", floor, key)[:nmin]) if nmin else np.array([], float)
                mrow[f"{key}_lo"], mrow[f"{key}_hi"] = _ci(diff)
            recs.append(mrow)

    df = pd.DataFrame(recs)
    df.to_parquet(DI.OUT_ROOT / "di_descriptors.parquet", index=False)
    num = df.select_dtypes(include="number").columns
    dfc = df.copy(); dfc[num] = dfc[num].round(4)
    dfc.to_csv(DI.OUT_ROOT / "di_descriptors.csv", index=False)
    print(f"\nwrote {DI.OUT_ROOT / 'di_descriptors.csv'}  "
          f"(scalar pool up to {int(df['n_draws'].max()) if not df.empty else 0})")

    # ---- density band: point KDE + 95% band over the full density pool ---- #
    dcols = _dcols(len(grid))
    band_rows = []
    for cohort in avail:
        for sx in ("M", "W"):
            for floor in FLOORS:
                _, _, _, _, d100 = arrays[(cohort, sx, floor)]
                pt = _kde(d100, grid)
                g = dpool.get((cohort, sx, floor))
                if g is not None and len(g):
                    arr = g[dcols].to_numpy(dtype=float)
                    lo = np.percentile(arr, DI.CI[0], axis=0)
                    hi = np.percentile(arr, DI.CI[1], axis=0)
                    nd = len(g)
                else:
                    lo = hi = np.full(len(grid), np.nan)
                    nd = 0
                for gi, x in enumerate(grid):
                    band_rows.append(dict(cohort=cohort, sex=sx, n_floor=floor,
                                          x=x, dens=pt[gi], lo=lo[gi], hi=hi[gi],
                                          n_draws=nd))
    pd.DataFrame(band_rows).to_csv(BAND_OUT, index=False)
    nb = max((r["n_draws"] for r in band_rows), default=0)
    print(f"wrote {BAND_OUT}  (band pool up to {nb})")

    # ---- readable markdown (scalar table) --------------------------------- #
    def cell(r, k, dp=3):
        if pd.isna(r[k]):
            return "-"
        s = f"{r[k]:.{dp}f}"
        lo, hi = r.get(f"{k}_lo"), r.get(f"{k}_hi")
        if pd.notna(lo):
            s += f" [{lo:.{dp}f}, {hi:.{dp}f}]"
        return s

    cols = ["Cohort", "Group", "n_i>=", "N", "disatt SD", "skew",
            "frac improver", "corr(d,u)", "corr(d,u) hi-n"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        n = "-" if pd.isna(r["n"]) else f"{int(r['n']):,}"
        lines.append("| " + " | ".join([
            r["cohort"], r["group"], str(int(r["n_floor"])), n,
            cell(r, "disatt_sd", 4), cell(r, "skew", 2),
            cell(r, "frac_improver", 3), cell(r, "corr_du", 3),
            cell(r, "corr_du_highn", 3),
        ]) + " |")
    note = (
        f"\n_Pool: up to {int(df['n_draws'].max()) if not df.empty else 0} athlete-"
        "bootstrap draws (accumulated). `d_i < 0` = improver. Units: log finish-"
        "time per year of career (x100 ~ %/yr). `[lo, hi]` = 95% CI; `M-W` = "
        "difference. **disatt SD** + **skew** are gauge-invariant; **frac improver**"
        " + **corr(d,u)** are EB-prior/gauge-anchored (d_i level not separately "
        "identified from the aging slope) and corr(d,u) is era-confounded -- read "
        "those two with care. Floors = post-hoc n_i subsets of the mrc2 fit._\n")
    (DI.OUT_ROOT / "di_descriptors.md").write_text(
        "# Career-drift d_i descriptors\n\n" + "\n".join(lines) + "\n" + note,
        encoding="utf-8")
    print(f"wrote {DI.OUT_ROOT / 'di_descriptors.md'}")


if __name__ == "__main__":
    main()
