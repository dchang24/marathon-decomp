"""Identity-mismatch perturbation engine for the §5.7 sensitivity test.

The decomposition's only input is a **partition of finish records into
athletes** (the identity resolution baked by ``load_slice``). Real record
linkage makes two kinds of error against that partition:

* **recall loss / under-merging** -- one true runner's finishes are split across
  several athlete ids  ->  here the ``break`` operation;
* **precision loss / over-merging** -- two distinct runners are fused into one
  athlete id           ->  here the ``join`` operation.

Production sits at a deliberately conservative operating point (splink p=0.999,
rule e=1e-6), so this module treats the **production slice as ground truth** and
injects *additional* error at a controlled per-athlete rate ``p`` to ask how far
the fit moves. ``break`` and ``join`` are deliberately symmetric inverses (one
cluster -> two; two clusters -> one), both parametrised by the fraction ``p`` of
athletes perturbed; ``both`` composes them.

The perturbation acts on the already-loaded baseline ``FitData`` (no parquet
re-read, no union-find): it relabels each cell's athlete, then re-runs the
**post-identity** steps of ``load_slice`` -- iterative ``min_race_count`` /
``min_runners_per_race`` pruning, single-connected-component reduction,
re-indexing, and ``A_n`` / ``A_e`` recomputation. That logic is duplicated here
(rather than imported) on purpose: this engine must stay a faithful mirror of
``data.load_slice`` steps 5-10, and the duplication is small and self-contained.

Entry point::

    pert = perturb(fd, op="join", p=0.05, rng=np.random.default_rng(1))
    new_fd = pert.fd            # a fresh FitData ready for any fitter
    pert.record_frac           # realised fraction of finishes reassigned
    pert.common_race_ids       # races surviving in BOTH fd and new_fd

Only the production operating point is supported (no entry-age band): an active
``age_lo``/``age_hi`` raises, since that filter is not reproduced here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from marathon_decomp.data import FitData

# Seconds per Julian year -- must match data._SEC_PER_YEAR exactly.
_SEC_PER_YEAR = 365.25 * 86400.0

Op = Literal["break", "join", "both"]
SplitMode = Literal["datecut", "iid", "singleton"]

# Default finish-time compatibility gate for a join: two athletes may be fused
# only if their median log-finish-times differ by <= TAU_LOGTIME. log(1.20) is
# ~20%, matching the finish-time mismatch threshold used upstream in linkage, so
# we never fuse e.g. a 2:10 and a 4:30 runner.
TAU_LOGTIME = float(np.log(1.20))
# yob gate: a false merge of two runners born far apart is implausible. The
# intended rule is "block if their yob *ranges* are disjoint AND the gap exceeds
# YOB_GAP_MAX years". FitData carries only the yob *midpoint* (recovered from
# A_e), not [yob_min, yob_max], so the range collapses to a point and the rule
# reduces to: block when both midpoints are known and |delta midpoint| > YOB_GAP_MAX.
# That is strictly more conservative than the range rule (it cannot stay lenient
# toward two wide-uncertain-range athletes) -- the safe direction here.
YOB_GAP_MAX = 6.0


# ---------------------------------------------------------------------------
# Cell view of a FitData (everything a re-partition needs)
# ---------------------------------------------------------------------------

def _cells_frame(fd: FitData) -> pd.DataFrame:
    """One row per finish cell, carrying the fields needed to re-partition and
    rebuild. Per-athlete attributes (sex/country/yob) are broadcast onto cells;
    per-race attributes (series/country/date) likewise, so the rebuild never
    needs to touch the original parquets.

    yob midpoint is recovered from the baseline ``A_e`` (entry age = first race
    year - yob midpoint), inverting data.load_slice step 8.
    """
    row = fd.row_idx
    col = fd.col_idx

    # per-athlete first-race year, to invert A_e -> yob midpoint
    race_year = pd.to_datetime(fd.race_date).year.to_numpy().astype(np.float64)
    cell_year = race_year[col]
    # invert A_e: yob midpoint = athlete first-race year - entry age
    cell_first_year_ath = pd.Series(cell_year).groupby(row).transform("min").to_numpy()
    yob_mid = cell_first_year_ath - fd.A_e[row]      # NaN where A_e unknown

    return pd.DataFrame({
        "cell": np.arange(fd.N, dtype=np.int64),
        "ath0": row.astype(np.int64),                 # baseline (truth) athlete idx
        "race_idx0": col.astype(np.int64),
        "race_id": fd.race_ids[col].astype(np.int64),
        "t_sec": np.exp(fd.y) if fd.response_kind == "log_time" else fd.y.copy(),
        "race_date": fd.race_date[col],               # datetime64[ns]
        "sex": fd.athlete_sex[row],
        "country": fd.athlete_country[row],
        "yob_mid": yob_mid,
        "race_series": fd.race_series[col],
        "race_country": fd.race_country[col],
    })


# ---------------------------------------------------------------------------
# Per-athlete summaries used by the join compatibility test
# ---------------------------------------------------------------------------

def _athlete_summary(cells: pd.DataFrame, label_col: str) -> pd.DataFrame:
    """Per-athlete median-log-time / sex / yob summary that defines the join
    window. The same-day gate is now built separately and vectorised inside
    ``_join_labels`` (as an integer day-membership index), so this no longer
    materialises a per-athlete day ``frozenset`` -- that Python ``apply`` over
    every athlete was itself a large slice of the old join cost."""
    g = cells.groupby(label_col)
    return pd.DataFrame({
        "med_logt": g["t_sec"].apply(lambda s: float(np.median(np.log(s.to_numpy())))),
        "sex": g["sex"].first(),
        "yob_mid": g["yob_mid"].mean(),    # NaN if all unknown
        "n": g.size(),
    })


# ---------------------------------------------------------------------------
# break: split a fraction p of athletes into two sub-runners (recall loss)
# ---------------------------------------------------------------------------

def _break_labels(
    cells: pd.DataFrame, labels: np.ndarray, p: float, rng: np.random.Generator,
    mode: SplitMode, next_label: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Split each of a random fraction ``p`` of eligible athletes (n>=2) into two
    sub-runners. Returns (new_labels, moved_mask, next_label). ``moved_mask`` is
    per-cell True where the cell was reassigned to a fresh sub-label.

    Split mechanics (``mode``):
      datecut   -- sort the athlete's finishes by date, cut at a uniform interior
                  index -> an early-career and a late-career sub-runner (default;
                  the realistic recall failure of a career drifting in name/age).
      iid       -- assign each finish to side A/B by a fair coin.
      singleton -- detach exactly one random finish into its own sub-runner
                  (the "one lost race" failure mode).
    """
    labels = labels.copy()
    moved = np.zeros(len(labels), dtype=bool)
    cell_arr = cells["cell"].to_numpy()
    order = np.argsort(cell_arr)         # cells are 0..N-1 already, but be safe
    assert np.array_equal(cell_arr[order], np.arange(len(labels)))

    # eligible = current labels with >= 2 cells
    by_label: dict[int, np.ndarray] = {}
    for lab, idx in pd.Series(np.arange(len(labels))).groupby(labels).groups.items():
        by_label[int(lab)] = np.asarray(idx, dtype=np.int64)
    eligible = [lab for lab, idx in by_label.items() if len(idx) >= 2]
    if not eligible or p <= 0:
        return labels, moved, next_label

    k = int(round(p * len(eligible)))
    if k <= 0:
        return labels, moved, next_label
    chosen = rng.choice(np.asarray(eligible), size=min(k, len(eligible)), replace=False)

    dates = pd.to_datetime(cells["race_date"]).to_numpy()
    for lab in chosen:
        idx = by_label[int(lab)]
        n = len(idx)
        if mode == "datecut":
            o = idx[np.argsort(dates[idx], kind="stable")]
            cut = int(rng.integers(1, n))           # 1..n-1
            side_b = o[cut:]
        elif mode == "iid":
            coin = rng.integers(0, 2, size=n).astype(bool)
            if coin.all() or (~coin).all():         # force a real split
                coin[rng.integers(0, n)] = ~coin[rng.integers(0, n)]
            side_b = idx[coin]
        elif mode == "singleton":
            side_b = idx[[int(rng.integers(0, n))]]
        else:  # pragma: no cover
            raise ValueError(f"unknown split mode {mode!r}")
        labels[side_b] = next_label
        moved[side_b] = True
        next_label += 1

    return labels, moved, next_label


# ---------------------------------------------------------------------------
# join: fuse a fraction p of athletes with a compatible partner (precision loss)
# ---------------------------------------------------------------------------

def _join_labels(
    cells: pd.DataFrame, labels: np.ndarray, p: float, rng: np.random.Generator,
    tau: float, yob_gap_max: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Fuse ~p*(#athletes) athletes into compatible partners. Returns
    (new_labels, moved_mask, n_merges) where moved cells are those of the
    *absorbed* cluster (relabelled to its partner), and ``n_merges`` is the
    number of merges actually formed. ``n_merges`` can be BELOW the target
    ``round(p*I/2)`` when the compatibility gates (sex / time / same-day / yob)
    leave an athlete with no eligible partner -- that shortfall is exactly why
    the realised join error rate is lower than the requested ``p``.

    Compatibility (mirrors the linkage gates + load_slice gatekeeping):
      * same sex (hard);
      * |delta median log finish-time| <= tau (~20%, so no implausible fusions);
      * no shared calendar day (a real runner can't finish two races a day);
      * if both yob midpoints known, within ``yob_gap_max`` years (a >2y age
        gap blocks the merge; see YOB_GAP_MAX on the range-vs-midpoint proxy).

    Candidate pairs are generated by a sliding window over athletes sorted by
    median log-time (near-linear). The static gates (sex / time / yob) and the
    same-day gate are applied with vectorised numpy over the window rather than a
    Python comprehension with a per-candidate ``frozenset`` intersection; the
    greedy matching order and every RNG draw are byte-for-byte unchanged, so this
    returns a partition IDENTICAL to the original scalar version (asserted in
    scratch/verify_join_vec.py) at a fraction of the cost on the large ALL
    slices, where the old inner loop was O(visits x window) interpreted ops.
    """
    labels = labels.copy()
    moved = np.zeros(len(labels), dtype=bool)
    if p <= 0:
        return labels, moved, 0

    summ = _athlete_summary(cells.assign(_lab=labels), "_lab")
    summ = summ.sort_values("med_logt")
    labs = summ.index.to_numpy()
    med = summ["med_logt"].to_numpy()
    sex = summ["sex"].to_numpy()
    yob = summ["yob_mid"].to_numpy()
    I = len(labs)
    pos = {int(l): k for k, l in enumerate(labs)}

    n_target = int(round(p * I / 2.0))        # each merge consumes 2 athletes
    if n_target <= 0:
        return labels, moved, 0

    cells_by_label: dict[int, np.ndarray] = {
        int(l): np.asarray(idx, dtype=np.int64)
        for l, idx in pd.Series(np.arange(len(labels))).groupby(labels).groups.items()
    }

    # --- vectorised gate precompute (replaces the per-candidate `compatible`) --
    # sex equality via integer codes: factorize maps equal values (incl. two
    # NaN/None -> -1) to equal codes, exactly mirroring the old `sex[a] != sex[b]`.
    sex_codes = pd.factorize(pd.Series(sex))[0]
    # same-day gate as an integer day-membership index in med-sorted position
    # space: for athlete a, the set of positions sharing >=1 calendar day with a
    # is the union of `day_pos[d]` over a's day codes -- equal to the old
    # `bool(days[a] & days[b])` test, but built once and applied vectorised.
    day_code = (pd.to_datetime(cells["race_date"]).dt.normalize()
                .factorize(sort=False)[0].astype(np.int64))
    if (day_code < 0).any():                  # NaT would break the CSR indexing
        raise ValueError("perturb join: encountered an unparseable race_date")
    posn = pd.Series(labels).map(pos).to_numpy().astype(np.int64)   # cell -> position
    pairs = pd.DataFrame({"posn": posn, "day": day_code}).drop_duplicates()
    pn = pairs["posn"].to_numpy()
    dy = pairs["day"].to_numpy()
    # CSR 1: athlete position -> its distinct day codes
    o = np.argsort(pn, kind="stable")
    ath_ptr = np.zeros(I + 1, dtype=np.int64)
    np.add.at(ath_ptr, pn[o] + 1, 1)
    ath_ptr = np.cumsum(ath_ptr)
    ath_day_flat = dy[o]
    # CSR 2: day code -> athlete positions racing that day
    n_days = int(dy.max()) + 1
    o2 = np.argsort(dy, kind="stable")
    day_ptr = np.zeros(n_days + 1, dtype=np.int64)
    np.add.at(day_ptr, dy[o2] + 1, 1)
    day_ptr = np.cumsum(day_ptr)
    day_pos_flat = pn[o2]
    day_scratch = np.zeros(I, dtype=bool)     # persistent; reset only where touched
    _EMPTY = np.empty(0, dtype=np.int64)

    used = np.zeros(I, dtype=bool)
    visit = rng.permutation(I)
    n_merged = 0
    for a in visit:
        if used[a] or n_merged >= n_target:
            continue
        # search outward from a within the median-log-time window [med-tau, med+tau]
        lo = np.searchsorted(med, med[a] - tau, side="left")
        hi = np.searchsorted(med, med[a] + tau, side="right")
        win = np.arange(lo, hi)
        # static gates (same conjunction as the old `compatible`, vectorised)
        m = (win != a) & (~used[win]) & (sex_codes[win] == sex_codes[a])
        m &= np.abs(med[win] - med[a]) <= tau          # re-check (FP parity w/ old)
        ya = yob[a]
        if np.isfinite(ya):                            # both-known + gap > max -> block
            yw = yob[win]
            m &= ~(np.isfinite(yw) & (np.abs(yw - ya) > yob_gap_max))
        # same-day gate
        d0, d1 = ath_ptr[a], ath_ptr[a + 1]
        if d1 > d0:
            segs = [day_pos_flat[day_ptr[d]:day_ptr[d + 1]] for d in ath_day_flat[d0:d1]]
            blocked = np.concatenate(segs) if segs else _EMPTY
            day_scratch[blocked] = True
            m &= ~day_scratch[win]
            day_scratch[blocked] = False
        cand = win[m]                                  # ascending b, == old list order
        if cand.size == 0:
            continue
        b = int(rng.choice(cand))
        used[a] = used[b] = True
        # absorb b into a (a keeps its label; b's cells take a's label)
        la, lb = int(labs[a]), int(labs[b])
        idx_b = cells_by_label[lb]
        labels[idx_b] = la
        moved[idx_b] = True
        n_merged += 1

    return labels, moved, n_merged


# ---------------------------------------------------------------------------
# Rebuild a FitData from a (cells, labels) partition -- mirrors load_slice 5-10
# ---------------------------------------------------------------------------

def _keep_giant_component(label: np.ndarray, race_id: np.ndarray) -> np.ndarray:
    """Boolean cell mask keeping the largest connected component of the
    athlete-race bipartite graph (duplicates data._keep_giant_component)."""
    _, ci = np.unique(label, return_inverse=True)
    ur, ri = np.unique(race_id, return_inverse=True)
    A, R = len(np.unique(ci)), len(ur)
    g = coo_matrix((np.ones(len(ci)), (ci, A + ri)), shape=(A + R, A + R))
    n_comp, comp = connected_components(g, directed=False)
    if n_comp <= 1:
        return np.ones(len(label), dtype=bool)
    cell_comp = comp[ci]
    giant = np.bincount(cell_comp).argmax()
    return cell_comp == giant


def _rebuild_fitdata(fd: FitData, cells: pd.DataFrame, labels: np.ndarray) -> FitData:
    """Reconstruct a FitData from a relabelled cell set, applying the same
    post-identity pruning + indexing + A_n/A_e recomputation as load_slice."""
    spec = fd.spec
    if spec.age_lo is not None or spec.age_hi is not None:
        raise NotImplementedError(
            "perturb only supports the production point (no entry-age band); "
            "load_slice's age-band loop is not reproduced here.")

    df = cells[["race_id", "t_sec", "race_date", "sex", "country", "yob_mid",
                "race_series", "race_country"]].copy()
    df["label"] = labels

    mrc = spec.min_race_count
    mr = spec.min_runners_per_race
    while True:
        prev = len(df)
        if mrc > 1:
            c = df.groupby("label")["label"].transform("size")
            df = df[c >= mrc]
        if mr > 1:
            c = df.groupby("race_id")["race_id"].transform("size")
            df = df[c >= mr]
        if len(df) == prev:
            break

    if len(df) == 0:
        raise ValueError("perturbation pruned the slice to 0 cells")

    keep = _keep_giant_component(df["label"].to_numpy(), df["race_id"].to_numpy())
    df = df[keep].reset_index(drop=True)

    # re-index athletes (by first appearance) and races (sorted by id, as load_slice)
    ath_codes, ath_uniques = pd.factorize(df["label"], sort=False)
    race_codes, race_uniques = pd.factorize(df["race_id"], sort=True)
    df["row_idx"] = ath_codes
    df["col_idx"] = race_codes
    I = len(ath_uniques)
    J = len(race_uniques)
    N = len(df)

    # per-race side info aligned to race_uniques (first cell per race)
    race_first = df.drop_duplicates("col_idx").set_index("col_idx").sort_index()
    race_date = pd.to_datetime(race_first["race_date"]).to_numpy()   # length J
    race_series = race_first["race_series"].to_numpy(dtype=object)
    race_country = race_first["race_country"].to_numpy(dtype=object)

    # A_n: cell race date minus the athlete's first race date, in years
    cell_race_date = race_date[df["col_idx"].to_numpy()]
    first_race_date = (pd.Series(cell_race_date)
                       .groupby(df["row_idx"].to_numpy()).min()
                       .reindex(range(I)).to_numpy())
    A_n = ((cell_race_date - first_race_date[df["row_idx"].to_numpy()])
           .astype("timedelta64[s]").astype(np.float64) / _SEC_PER_YEAR)

    # per-athlete attributes aligned to row index
    ath_first = df.drop_duplicates("row_idx").set_index("row_idx").sort_index()
    athlete_sex = ath_first["sex"].to_numpy(dtype=object)
    athlete_country = ath_first["country"].to_numpy(dtype=object)
    yob_mid = df.groupby("row_idx")["yob_mid"].mean().reindex(range(I)).to_numpy()
    first_race_year = pd.to_datetime(first_race_date).year.astype("float64").to_numpy()
    A_e = first_race_year - yob_mid

    # athlete id = the new partition label (unique & stable). Original surrogate
    # ids are not preserved: after a split one parent maps to two new athletes,
    # so ids could no longer be unique -- and nothing downstream joins on them.
    athlete_ids = ath_uniques.astype(np.int64)

    t_sec = df["t_sec"].to_numpy(dtype=np.float64)
    y = np.log(t_sec) if spec.response == "log_time" else t_sec

    return FitData(
        row_idx=df["row_idx"].to_numpy(dtype=np.int64),
        col_idx=df["col_idx"].to_numpy(dtype=np.int64),
        y=y.astype(np.float64),
        A_n=A_n.astype(np.float64),
        A_e=A_e.astype(np.float64),
        athlete_ids=athlete_ids,
        athlete_sex=athlete_sex,
        athlete_country=athlete_country,
        race_ids=race_uniques.astype(np.int64),
        race_series=race_series,
        race_country=race_country,
        race_date=race_date,
        I=I, J=J, N=N,
        spec=spec,
        data_version=fd.data_version,
        cache_key=fd.cache_key,
        response_kind=fd.response_kind,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

@dataclass
class Perturbation:
    fd: FitData                      # rebuilt, perturbed slice
    op: Op
    p: float
    record_frac: float               # fraction of finishes reassigned vs truth
    n_break_moves: int               # cells moved by break
    n_join_moves: int                # cells moved by join
    n_break_splits: int              # athletes actually split (realised, pre-prune)
    n_join_merges: int               # merges actually formed (<= round(p*I/2))
    n_athletes_base: int             # baseline athlete count I the rates are over
    common_race_ids: np.ndarray      # race ids present in BOTH baseline & perturbed


def perturb(
    fd: FitData, op: Op, p: float, rng: np.random.Generator, *,
    split_mode: SplitMode = "datecut", tau: float = TAU_LOGTIME,
    yob_gap_max: float = YOB_GAP_MAX,
) -> Perturbation:
    """Apply one mismatch perturbation to ``fd`` and return a rebuilt slice.

    ``op``: 'break' (recall loss / split), 'join' (precision loss / merge), or
    'both'. ``p`` is the fraction of athletes perturbed *per operation*; 'both'
    applies join at rate ``p`` and THEN break at rate ``p`` independently, so it
    is ~p+p ~= 2p total edits (NOT a single budget of p split between the two) --
    hence 'both' at p ~= break at p + join at p (the ``record_frac`` is additive).
    The returned ``record_frac`` is the realised fraction of finishes whose
    cluster membership changed -- the record-level error rate, for plotting
    against the linkage system's own precision/recall.

    The *requested* rate ``p`` is an upper bound on what actually happens: the
    join compatibility gates (sex / time / same-day / yob) can leave an athlete
    with no eligible partner, so fewer merges form than the ``round(p*I/2)``
    target. The realised counts ``n_break_splits`` / ``n_join_merges`` (over the
    baseline athlete count ``n_athletes_base``) make the true per-athlete rate
    reconstructable per run; ``record_frac`` gives the realised finish-level
    rate. (break is essentially ungated, so its realised rate tracks ``p``.)
    """
    if op not in ("break", "join", "both"):
        raise ValueError(f"unknown op {op!r}")

    cells = _cells_frame(fd)
    labels = cells["ath0"].to_numpy().copy()
    # baseline athlete count -- the denominator for the realised per-athlete rate.
    n_ath_base = int(len(np.unique(labels)))
    next_label0 = int(labels.max()) + 1
    next_label = next_label0
    moved_join = np.zeros(fd.N, dtype=bool)
    moved_break = np.zeros(fd.N, dtype=bool)
    n_merges = 0

    if op in ("join", "both"):
        labels, moved_join, n_merges = _join_labels(
            cells, labels, p, rng, tau, yob_gap_max)
    if op in ("break", "both"):
        labels, moved_break, next_label = _break_labels(
            cells, labels, p, rng, split_mode, next_label)
    # one fresh sub-label is minted per split, so the label counter advanced by
    # exactly the number of athletes actually split.
    n_splits = next_label - next_label0

    new_fd = _rebuild_fitdata(fd, cells, labels)
    moved = moved_join | moved_break
    common = np.intersect1d(fd.race_ids, new_fd.race_ids)
    return Perturbation(
        fd=new_fd, op=op, p=float(p),
        record_frac=float(moved.mean()),
        n_break_moves=int(moved_break.sum()),
        n_join_moves=int(moved_join.sum()),
        n_break_splits=int(n_splits),
        n_join_merges=int(n_merges),
        n_athletes_base=n_ath_base,
        common_race_ids=common,
    )
