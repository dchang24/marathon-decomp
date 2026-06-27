"""Pairwise dominance + credible tiers for the single-slice series ranking.

The N-free replacement for q03's `P(top-N)`. `P(top-N)` only resolves series
sitting near the rank-N boundary (everything robustly inside saturates at 1,
outside at 0), so it throws away the interior structure and depends on an
arbitrary N. Here we use the full bootstrap instead.

PER-REPLICATE SERIES MEDIANS. As in q03, the series index is m_s = median_j(v_j)
over its editions (n_j >= --min-n), recomputed INSIDE each bootstrap replicate b
so shared-athlete correlation between a series' editions is preserved:

    M[b, s] = median_j( v_j^(b) )            (R replicates x S series)

DOMINANCE. For every ordered pair of series (s, s'),

    D[s, s'] = mean_b( M[b, s] < M[b, s'] ) = P( s faster than s' )

a full who-beats-whom structure with no N and no boundary. D[s,s'] near 1 means
s is almost surely faster; near 0.5 means the two are bootstrap-indistinguishable.

TIERS (--thr, default 0.9). Series are sorted fastest-first by point m_s and
swept into ordered tiers: a series opens a NEW tier iff the current tier's
LEADER (its fastest member) is confidently faster than it, D[leader, s] >= thr.
So within a tier no member is confidently slower than the tier leader, and each
new tier's leader is confidently slower than the previous leader. Tiers are the
small set of genuinely distinguishable groups; within-tier order is not
bootstrap-resolvable and is not claimed.

Also per series: n_faster_conf / n_slower_conf = how many other series it is
confidently (>= thr) faster / slower than -- a robust, N-free centrality score.

Outputs (results/analysis/race_comparison/{slice}/):
    series_dominance.csv          one row per series: m_s + 95% CI + tier +
                                  n_faster_conf/n_slower_conf, sorted fast->slow
    series_dominance_matrix.parquet   the full S x S D matrix (long form)
    series_boot_medians.parquet       per-replicate series medians (run_id x
                                      series_key), the raw material for box /
                                      violin marker styles
    + tier groupings to stdout.

NO REFIT -- reads the saved point fit + its bootstrap.

Run::

    python scripts/05_analysis/race_comparison/q05_series_dominance.py
    python scripts/05_analysis/race_comparison/q05_series_dominance.py \
        --cohorts ALL --sexes B --min-editions 3 --thr 0.9
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from marathon_decomp import load_slice  # noqa: E402
from marathon_decomp.config import DATA_DIR  # noqa: E402
from race_common import (  # noqa: E402
    OUT_ROOT, v_to_min, find_fit_dir, load_boot_wide, load_editions, pct,
)


def sex_counts(slug, model, nutag, data_version):
    """Per-series M / W finish-record counts inside the slice (fit component)."""
    fit_dir = find_fit_dir(slug, model, nutag)
    if fit_dir is None:
        return None
    payload = pickle.load(open(fit_dir / "fit.pkl", "rb"))
    fd = load_slice(payload["spec"], payload.get("data_version", data_version))
    ath = pd.read_parquet(DATA_DIR / "athletes.parquet",
                          columns=["athlete_id", "sex"])
    sex_by = dict(zip(ath["athlete_id"].to_numpy(), ath["sex"].to_numpy()))
    aid = np.asarray(fd.athlete_ids)[np.asarray(fd.row_idx, int)]
    obs_sex = np.array([sex_by.get(int(a), "?") for a in aid])
    rid = np.asarray(fd.race_ids)[np.asarray(fd.col_idx, int)]
    comps = pd.read_parquet(DATA_DIR / "competitions.parquet",
                            columns=["race_id", "series_key"])
    r2s = dict(zip(comps["race_id"].to_numpy(), comps["series_key"].to_numpy()))
    obs_series = np.array([r2s.get(int(r), "unknown") for r in rid])
    g = (pd.DataFrame({"series_key": obs_series, "sex": obs_sex})
         .groupby(["series_key", "sex"]).size().unstack(fill_value=0))
    for c in ("M", "W"):
        if c not in g:
            g[c] = 0
    return g[["M", "W"]].reset_index().rename(columns={"M": "n_men", "W": "n_women"})


# --------------------------------------------------------------------------- #
def series_median_matrix(slug, model, nutag, data_version, min_n, min_editions):
    """-> (sr point DataFrame, M (R,S) per-replicate series-median matrix or None)."""
    loaded = load_editions(slug, model, nutag, data_version)
    if loaded is None:
        return None
    ed, fit_dir = loaded
    ed = ed[ed["n_j"] >= min_n].reset_index(drop=True)
    sizes = ed.groupby("series_key").size()
    keep = sizes[sizes >= min_editions].index
    ed = ed[ed["series_key"].isin(keep)].reset_index(drop=True)
    if ed["series_key"].nunique() < 3:
        print(f"  [skip] {slug}: only {ed['series_key'].nunique()} eligible series")
        return None

    g = ed.groupby("series_key")
    sr = pd.DataFrame({
        "series_key": g.size().index,
        "country": g["country"].first().to_numpy(),
        "k_editions": g.size().to_numpy(int),
        "yr_min": g["year"].min().to_numpy(int),
        "yr_max": g["year"].max().to_numpy(int),
        "median_v": g["v"].median().to_numpy(np.float64),
    }).reset_index(drop=True)

    B = load_boot_wide(fit_dir, ed["race_id"].to_numpy(np.int64),
                       ed["v"].to_numpy(np.float64))
    if B is None:
        print(f"  [skip] {slug}: no bootstrap (dominance needs replicates)")
        return None
    cols = [np.flatnonzero((ed["series_key"] == s).to_numpy())
            for s in sr["series_key"]]
    M = np.column_stack([np.median(B[:, c], axis=1) for c in cols])   # (R, S)
    return sr, M


def assign_tiers(order, D, thr):
    """Sweep fastest-first; new tier when the tier leader dominates the series."""
    tiers = np.empty(len(order), int)
    t = 0
    leader = order[0]
    for k, s in enumerate(order):
        if k > 0 and D[leader, s] >= thr:   # leader confidently faster than s
            t += 1
            leader = s
        tiers[s] = t
    return tiers


def analyse(slug, *, model, nutag, data_version, min_n, min_editions, thr):
    got = series_median_matrix(slug, model, nutag, data_version,
                               min_n, min_editions)
    if got is None:
        return None
    sr, M = got
    S = M.shape[1]
    keys_orig = sr["series_key"].to_numpy().copy()   # M columns align to this

    # D[s, s'] = P(s faster than s') = P(m_s < m_s')
    # vectorized: for each pair compare columns across replicates
    D = np.empty((S, S))
    for s in range(S):
        D[s] = (M[:, [s]] < M).mean(axis=0)
    np.fill_diagonal(D, np.nan)

    m = sr["median_v"].to_numpy(np.float64)
    sr["min_at_3h"] = v_to_min(m)
    sr["m_lo95"], sr["m_hi95"] = pct(M, 2.5), pct(M, 97.5)
    # robust centrality: how many series each is confidently faster / slower than
    sr["n_faster_conf"] = np.nansum(D >= thr, axis=1).astype(int)        # rows: s faster than .
    sr["n_slower_conf"] = np.nansum(D <= (1 - thr), axis=1).astype(int)
    sr["n_boot"] = M.shape[0]

    order = np.argsort(m)                       # fastest first
    tiers = assign_tiers(order, D, thr)
    sr["tier"] = tiers + 1                       # 1 = fastest tier
    sr = sr.iloc[order].reset_index(drop=True)   # sort fast -> slow

    sr.insert(0, "slice", slug)
    sr.insert(1, "model", f"{model}_{nutag}")
    sr.insert(2, "thr", thr)

    sc = sex_counts(slug, model, nutag, data_version)
    if sc is not None:
        sr = sr.merge(sc, on="series_key", how="left")
        sr[["n_men", "n_women"]] = sr[["n_men", "n_women"]].fillna(0).astype(int)

    keys = sr["series_key"].to_numpy()           # in fast->slow order
    order_pos = {s_old: pos for pos, s_old in enumerate(order)}
    # long-form D in the sorted order for a tidy matrix file
    idx = order
    D_sorted = D[np.ix_(idx, idx)]
    long = pd.DataFrame({
        "slice": slug,
        "faster": np.repeat(keys, S),
        "slower": np.tile(keys, S),
        "p_faster": D_sorted.ravel(),
    })
    long = long[long["faster"] != long["slower"]].reset_index(drop=True)

    # per-replicate series-median matrix (for box/violin downstream)
    mdf = pd.DataFrame(M, columns=keys_orig)
    mdf.insert(0, "run_id", np.arange(M.shape[0], dtype=int))
    return sr, long, mdf


def print_tiers(sr):
    print(f"\n  tiers (thr applied) for {sr['slice'].iloc[0]}:")
    for t, grp in sr.groupby("tier"):
        names = ", ".join(f"{r.series_key}({r.min_at_3h:+.1f})"
                          for r in grp.itertuples())
        print(f"    Tier {t:>2} [{len(grp):>2}]: {names}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="+", default=["ALL"])
    ap.add_argument("--sexes", nargs="+", default=["B"])
    ap.add_argument("--mrcs", nargs="+", default=["mrc2"])
    ap.add_argument("--window", default="14-25")
    ap.add_argument("--model", default="full")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--min-n", type=int, default=0)
    ap.add_argument("--min-editions", type=int, default=3)
    ap.add_argument("--thr", type=float, default=0.9,
                    help="dominance probability for a confident faster/slower call.")
    args = ap.parse_args()

    n_done = 0
    for cohort in args.cohorts:
        for sex in args.sexes:
            for mrc in args.mrcs:
                slug = f"{cohort}_{sex}_{args.window}_{mrc}"
                print(f"analysing {slug} ...")
                out = analyse(slug, model=args.model, nutag=args.nutag,
                              data_version=args.data_version, min_n=args.min_n,
                              min_editions=args.min_editions, thr=args.thr)
                if out is None:
                    continue
                sr, long, mdf = out
                out_dir = OUT_ROOT / slug
                out_dir.mkdir(parents=True, exist_ok=True)
                num = sr.select_dtypes(include="number").columns
                sr[num] = sr[num].round(6)
                sr.to_csv(out_dir / "series_dominance.csv", index=False)
                long.to_parquet(out_dir / "series_dominance_matrix.parquet",
                                index=False)
                mdf.to_parquet(out_dir / "series_boot_medians.parquet",
                               index=False)
                print(f"  wrote {out_dir / 'series_dominance.csv'}  "
                      f"({len(sr)} series, {sr['tier'].max()} tiers)")
                print(f"  wrote {out_dir / 'series_dominance_matrix.parquet'}")
                print(f"  wrote {out_dir / 'series_boot_medians.parquet'}")
                print_tiers(sr)
                n_done += 1

    if not n_done:
        raise SystemExit("no slice produced a result")
    print()


if __name__ == "__main__":
    main()
