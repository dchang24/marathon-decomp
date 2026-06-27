"""Within-athlete permutation null for the v_j-bias correlation (refit-free).

q03 shows ``corr(delta_LOO, v_no-d) > 0`` after year-partialling. Is that real
within-athlete trajectory signal, or could the LOO construction manufacture it?
This puts a p-value on it by a permutation that destroys exactly the trajectory:

  * Observed: T_obs = year-partialled pearson( delta_LOO , v ).
  * Null draw: shuffle the no-d residual ``r`` **within each athlete** (preserving
    each athlete's residual multiset and the race assignment), which breaks the
    within-athlete an<->r pairing; recompute delta_LOO on the shuffled r; recorrelate
    with the (fixed) v. Under H0 (career-stage unrelated to performance) the slopes
    b_i are null, so delta is null and T ~ 0.

The field composition (the an_tilde values present at each race) and v_j are held
fixed across draws, so only the trajectory signal is nulled. Reported for both v
fits (no-d=agingS4gv, +d=full) off the same draws: we expect the no-d statistic
far in the upper tail and the +d statistic buried in the null.

EB is NOT permuted here -- d_i^full cannot be recomputed without a refit, so the
permutation test is LOO-only (the independent, refit-free estimator). This is the
inference step for the descriptive q03 result.

Outputs -- ONE SET OF FILES PER SLICE under results/validation/d_i/permutation/:
  * <slice>.draws.csv  -- the ACCUMULATOR: one row per (slug, fit, draw k) -> t_perm.
                          ``--n-perm`` ADDS draws (k continues past what's stored);
                          runs accumulate (200 now + 50 + 500 = 750 distinct draws).
  * <slice>.csv        -- summary per (slug, fit): T_obs + p recomputed over the FULL
                          pool (null_mean/sd, z, p_one, p_two, n_perm_total).
  * <slice>.md         -- the formatted run for that slice.

Each draw k is reproducible by index (SeedSequence([seed, crc32(slug), k])), so
p-values only ever get *more* accurate as draws accumulate -- never overwritten.

Run::

    python scripts/04_validation/d_i/q04_perm_null.py --slices Po10_W            # +1000 draws
    python scripts/04_validation/d_i/q04_perm_null.py --slices Po10_W --n-perm 200   # add 200 more
    python scripts/04_validation/d_i/q04_perm_null.py                            # all 8, mrc 2 & 5
"""
from __future__ import annotations

import argparse
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import load_slice, registry                # noqa: E402

from baseline_common import slices as S                          # noqa: E402
import delta_common as DC                                        # noqa: E402  (sibling)

FITS = ("no-d", "+d")
PERM_DIR = DC.OUT_ROOT / "permutation"   # one file per slice lives here
DRAW_COLS = ["slug", "fit", "k", "t_perm"]


def slice_summary_csv(name: str) -> Path:
    return PERM_DIR / f"{name}.csv"            # recomputed p-values per slug/fit


def slice_draws_csv(name: str) -> Path:
    return PERM_DIR / f"{name}.draws.csv"      # the accumulator: one row per draw


def slice_md(name: str) -> Path:
    return PERM_DIR / f"{name}.md"


def _slug_seed(slug: str) -> int:
    """Stable (process-independent) per-slug seed; Python's hash() is randomized."""
    return zlib.crc32(slug.encode()) & 0xFFFFFFFF


def _draw_rng(seed: int, slug: str, k: int) -> np.random.Generator:
    """RNG for draw k of `slug`: reproducible by index, independent across (slug, k).

    Keying on k (not a sequential stream) makes draw k identical regardless of the
    batch it was computed in, so accumulating 200 + 50 + 500 draws gives 750
    distinct, reproducible permutations with no overlap.
    """
    return np.random.default_rng(np.random.SeedSequence([seed, _slug_seed(slug), k]))


def _partial_corr(delta: np.ndarray, v_part: np.ndarray, year: np.ndarray) -> float:
    """corr(partial(delta), v_part); v already year-residualized (fixed across draws)."""
    return DC.pearson(DC.partial_on_year(delta, year), v_part)


def _existing_draws(name: str, slug: str) -> pd.DataFrame:
    """All draws already on disk for this slug (empty frame if none)."""
    p = slice_draws_csv(name)
    if not p.is_file():
        return pd.DataFrame(columns=DRAW_COLS)
    df = pd.read_csv(p)
    return df[df["slug"] == slug]


def run_cell(name: str, mrc: int, nu: float, n_new: int, seed: int,
             rep: DC.Report) -> bool:
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    adir, fdir = DC.aging_dir(spec, nu), DC.full_dir(spec, nu)
    if not (DC.present(adir) and DC.present(fdir)):
        return False

    fd = load_slice(spec)
    amodel = registry.load_fit(adir, fd)
    r = fd.y - DC.yhat(amodel)
    an = DC.an_tilde(fd)
    year = DC.race_year(fd)
    row = fd.row_idx

    # year-residualize each v once (fixed across permutation draws)
    v_part = {
        "no-d": DC.partial_on_year(np.asarray(amodel.params["v"], float), year),
        "+d": DC.partial_on_year(
            np.asarray(registry.load_fit(fdir, fd).params["v"], float), year),
    }
    delta_obs, _ = DC.loo_delta(fd, r, an)
    t_obs = {f: _partial_corr(delta_obs, v_part[f], year) for f in FITS}

    # continue the draw index past whatever is already stored for this slug
    ex_slug = _existing_draws(name, slug)
    n0 = int(ex_slug["k"].max()) + 1 if len(ex_slug) else 0

    rep.head(f"{slug}   nu={nu:g}   N={fd.N:,}  J={fd.J:,}   "
             f"(+{n_new} draws; had {n0})", level=2)
    print(f"    T_obs: no-d={t_obs['no-d']:+.4f}  +d={t_obs['+d']:+.4f} ... "
          f"drawing k={n0}..{n0 + n_new - 1}", flush=True)

    new_rows: list[dict] = []
    tick = max(1, n_new // 10)
    t0 = time.perf_counter()
    for j, k in enumerate(range(n0, n0 + n_new)):
        rng = _draw_rng(seed, slug, k)
        r_perm = DC.permute_within_athlete(r, row, rng)
        delta_p, _ = DC.loo_delta(fd, r_perm, an)
        for f in FITS:
            new_rows.append(dict(slug=slug, fit=f, k=int(k),
                                 t_perm=_partial_corr(delta_p, v_part[f], year)))
        if (j + 1) % tick == 0:
            dt = time.perf_counter() - t0
            print(f"      {j + 1}/{n_new} new draws  ({dt:.1f}s, "
                  f"{1000 * dt / (j + 1):.1f} ms/draw)", flush=True)

    # --- append draws to the slice's accumulator (idempotent on k) -------
    PERM_DIR.mkdir(parents=True, exist_ok=True)
    DC.upsert_csv(slice_draws_csv(name), pd.DataFrame(new_rows, columns=DRAW_COLS),
                  keys=["slug", "fit", "k"])
    draws_all = pd.read_csv(slice_draws_csv(name))
    draws_slug = draws_all[draws_all["slug"] == slug]

    # --- recompute p-values over the FULL accumulated pool ---------------
    out = []
    for f in FITS:
        tn = draws_slug[draws_slug["fit"] == f]["t_perm"].to_numpy()
        n_tot = int(tn.size)
        mu = float(tn.mean()) if n_tot else float("nan")
        sd = float(tn.std(ddof=1)) if n_tot > 1 else float("nan")
        z = (t_obs[f] - mu) / sd if np.isfinite(sd) and sd > 1e-12 else float("nan")
        p_one = (1 + int((tn >= t_obs[f]).sum())) / (n_tot + 1)
        p_two = (1 + int((np.abs(tn) >= abs(t_obs[f])).sum())) / (n_tot + 1)
        out.append(dict(slug=slug, slice=name, mrc=mrc, nu=float(nu), fit=f,
                        n_perm_total=n_tot, T_obs=t_obs[f], null_mean=mu,
                        null_sd=sd, z=z, p_one=p_one, p_two=p_two))
    rep.table(pd.DataFrame(out)[["fit", "n_perm_total", "T_obs", "null_mean",
                                 "null_sd", "z", "p_one", "p_two"]], floatfmt="{:.4f}")
    verdict = ("no-d significant (p_one={:.4g}, n={}), +d {}".format(
        out[0]["p_one"], out[0]["n_perm_total"],
        "null as expected" if out[1]["p_one"] > 0.05 else f"also low (p_one={out[1]['p_one']:.4g})"))
    rep.line(f"\n{verdict}")

    DC.upsert_csv(slice_summary_csv(name), pd.DataFrame(out), keys=["slug"])
    rep.save(slice_md(name), f"d_i v_j-bias permutation null -- {name}")
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5])
    ap.add_argument("--nu", type=float, default=8.0)
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="permutation draws to ADD this run (accumulates across "
                         "runs; p recomputed over the full pool).")
    ap.add_argument("--seed", type=int, default=0,
                    help="master seed; per-draw RNG = SeedSequence([seed, crc32(slug), k]).")
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    PERM_DIR.mkdir(parents=True, exist_ok=True)
    print(f"d_i v_j-bias PERMUTATION NULL (LOO, within-athlete residual shuffle)\n"
          f"nu={args.nu:g}   mrc={mrcs}   +{args.n_perm} draws/cell this run   "
          f"slices: {', '.join(names)}\n"
          f"per-slice draws + summary -> {PERM_DIR}", flush=True)

    n_slices = 0
    for name in names:
        # a fresh per-slice report -> permutation/<slice>.md
        rep = DC.Report()
        rep.line(f"d_i v_j-bias permutation null -- {name}   nu={args.nu:g}   "
                 f"n_perm={args.n_perm}")
        rep.line("H0: within-athlete career-stage unrelated to performance. Expect "
                 "the no-d statistic in the upper tail, the +d statistic in the null.")
        did = False
        for mrc in mrcs:
            # each cell upserts permutation/<slice>.csv + re-saves the .md on finish
            if run_cell(name, mrc, args.nu, args.n_perm, args.seed, rep):
                did = True
        if not did:
            rep.line(f"\n[skip] {name}: no (mrc) had both fits.")
            continue
        # per-slice headline from this slice's own accumulated summary
        full = pd.read_csv(slice_summary_csv(name)).sort_values(["slug", "fit"])
        rep.head(f"{name} summary (p over full accumulated pool)", level=2)
        rep.table(full[["slug", "fit", "n_perm_total", "T_obs", "z", "p_one", "p_two"]],
                  floatfmt="{:.4f}")
        rep.save(slice_md(name), f"d_i v_j-bias permutation null -- {name}")
        n_slices += 1

    if n_slices == 0:
        print("\nNo slice had both fits. Run the production fits first.", flush=True)
        return
    print(f"\nDone: {n_slices} slice file(s) under {PERM_DIR}; "
          f"+{args.n_perm} draws each cell this run.", flush=True)


if __name__ == "__main__":
    main()
