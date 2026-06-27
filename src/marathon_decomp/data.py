"""Data loading and slicing for marathon decomposition models.

Single entry point: `load_slice(spec, data_version)` returns a `FitData`.
Every downstream consumer (fitters, diagnostics, bootstrap) takes a `FitData`
— no other module should touch the parquets directly.

Identity is pre-resolved at export time (the extractor runs the splink+rule
union-find and gatekeeping), so `athlete_id` here is already the physical-runner
id — this module does no clustering.

Layout of one slice operation (in order):
    1. read 3 parquets
    2. filter athletes (cohort, sex, country, yob requirement)
    3. filter races (series, country, date)
    4. filter results, pick time column, apply finish-time bounds
    5. enforce min_race_count after all upstream filters
    6. keep the single largest connected component of the design graph
    7. re-index athletes 0..I-1 and races 0..J-1
    8. compute per-cell A_n and per-athlete A_e
    9. apply response transform (log_time | time)
   10. pack into FitData

NaN propagation: A_e is NaN for athletes without resolvable DOB. The fitter
masks those out of the gamma term; nothing in this module decides eligibility
beyond what the user asked for via the spec.
"""
from __future__ import annotations

import hashlib
import json
import pickle
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .config import DATA_DIR, CACHE_DIR

# Bump when load_slice's output semantics change (filtering, A_n/A_e, response)
# so stale on-disk FitData caches are invalidated automatically.
#   1 -> initial
#   2 -> add sex/yob/same-day cluster conflict gatekeeping
#   3 -> keep only the largest connected component of the design graph
#   4 -> identity now pre-resolved at export; loader no longer clusters/gatekeeps
LOADER_VERSION = "4"

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

Cohort = Literal["ALL", "WMM", "Po10", "WA", "NYRR"]
Sex = Literal["M", "W", "ALL"]
TimeKind = Literal["chip", "gun", "chip_or_gun"]
ResponseKind = Literal["log_time", "time"]

_COHORT_TO_FLAG: dict[str, str | None] = {
    "ALL": None,
    "WMM": "has_wmm",
    "Po10": "has_po10",
    "WA": "has_wa",
    "NYRR": "has_nyrr",
}

# Seconds per Julian year — used for A_n / A_e in *years*.
_SEC_PER_YEAR = 365.25 * 86400.0


# ---------------------------------------------------------------------------
# Spec & container dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SliceSpec:
    """User-facing knobs for selecting a subset of the raw export.

    Every filter is optional and has a permissive default. Defaults are
    chosen to be inclusive ("ALL", no range bound, no merging) so that the
    user is opting *in* to every restriction explicitly.
    """

    # --- athlete cohort & demographics ---
    cohort: Cohort = "ALL"
    sex: Sex = "ALL"
    countries_include: tuple[str, ...] | None = None    # athlete countries
    countries_exclude: tuple[str, ...] = ()
    require_yob: bool = False                            # drop athletes w/o DOB info
    # Entry-age band (years), evaluated as first-race-year minus the athlete's
    # yob midpoint. Drops whole athletes outside the band. Athletes without
    # any yob info are kept (use require_yob to drop those). None => no bound.
    age_lo: float | None = None
    age_hi: float | None = None

    # --- race subset ---
    races_include: tuple[str, ...] | None = None         # by series_key
    races_exclude: tuple[str, ...] = ()
    race_countries_include: tuple[str, ...] | None = None
    race_countries_exclude: tuple[str, ...] = ()
    date_lo: date | None = date(2014, 1, 1)
    date_hi: date | None = date(2025, 12, 31)

    # --- result-level filters ---
    time_kind: TimeKind = "chip"
    finish_time_lo_sec: float = 7100       # sligtly below 2 hours (new sub2 WR in 2026)
    finish_time_hi_sec: float = 10*3600.0  # 10 hours

    # --- athlete reliability ---
    min_race_count: int = 2                              # applied AFTER all above
    min_runners_per_race: int = 20                       # applied AFTER min_race_count

    # --- response ---
    response: ResponseKind = "log_time"

    def cache_key(self) -> str:
        """Stable short hash of this spec — for on-disk caching of FitData."""
        payload = json.dumps(asdict(self), default=str, sort_keys=True)
        return hashlib.sha256(payload.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class FitData:
    """Frozen cell-indexed dataset ready for a fitter.

    Conventions:
        - i in [0, I) indexes athletes; j in [0, J) indexes races.
        - row_idx, col_idx, y, A_n have length N (one per observed cell).
        - A_e has length I (per-athlete entry age in years; NaN where DOB unknown).
        - athlete_ids / race_ids hold the *upstream surrogate ids* so results
          can be joined back to the parquets. athlete_id is the export's
          pre-resolved cluster id (one physical runner).
    """

    # cell-level (length N)
    row_idx: np.ndarray            # int64
    col_idx: np.ndarray            # int64
    y: np.ndarray                  # float64
    A_n: np.ndarray                # float64, years since athlete's first race

    # per-athlete (length I)
    A_e: np.ndarray                # float64, entry age in years; NaN allowed
    athlete_ids: np.ndarray        # int64; pre-resolved cluster id
    athlete_sex: np.ndarray        # 'M' / 'W'
    athlete_country: np.ndarray    # str

    # per-race (length J)
    race_ids: np.ndarray           # int64, surrogate race_id
    race_series: np.ndarray        # str, series_key
    race_country: np.ndarray       # str
    race_date: np.ndarray          # datetime64[ns]

    # sizes
    I: int
    J: int
    N: int

    # provenance
    spec: SliceSpec
    data_version: str
    cache_key: str
    response_kind: ResponseKind    # tags what y is (log_time | time)


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def load_slice(
    spec: SliceSpec,
    data_version: str | None = None,
    data_root: Path | None = None,
    cache: bool = True,
) -> FitData:
    """
    Apply `spec` to the raw export at `data_root` (default config.DATA_DIR).

    The export is read directly from that directory (no version subdir). Identity
    is already resolved at export time, so this does no clustering. `data_version`
    is a provenance label only; when None it's taken from the export's manifest.

    With `cache=True` (default) the resulting FitData is memoized on disk under
    config.CACHE_DIR, keyed by the spec, the data version, a fingerprint of the
    export's manifest, and LOADER_VERSION. The entry is reused only when all of
    those match, so changing any param, re-exporting the data, or bumping the
    loader invalidates it automatically. Pass `cache=False` to bypass.
    """

    root = data_root or DATA_DIR
    if data_version is None:
        data_version = _read_data_version(root)

    cache_file = _cache_path(spec, data_version, root) if cache else None
    if cache_file is not None and cache_file.exists():
        with open(cache_file, "rb") as f:
            return pickle.load(f)

    athletes = pd.read_parquet(root / "athletes.parquet")
    competitions = pd.read_parquet(root / "competitions.parquet")
    results = pd.read_parquet(root / "results.parquet")

    # --- 2. athlete filter -------------------------------------------------
    keep_ath_mask = _athlete_mask(athletes, spec)
    athletes = athletes.loc[keep_ath_mask].reset_index(drop=True)

    # --- 3. race filter ----------------------------------------------------
    keep_race_mask = _race_mask(competitions, spec)
    competitions = competitions.loc[keep_race_mask].reset_index(drop=True)

    # Identity is pre-resolved at export: athlete_id is the cluster id (one
    # physical runner, already gatekept). Alias it as canonical_id so the
    # connectivity / indexing code below reads uniformly.
    athletes["canonical_id"] = athletes["athlete_id"]
    keep_canonical = set(athletes["canonical_id"].tolist())

    # --- 4. result-level filters ------------------------------------------
    results = results[results["athlete_id"].isin(keep_canonical)].copy()
    results["canonical_id"] = results["athlete_id"]
    results = results[results["race_id"].isin(set(competitions["race_id"].tolist()))]

    t_sec = _pick_time_column(results, spec.time_kind)
    results = results.assign(t_sec=t_sec).dropna(subset=["t_sec"])
    results = results[
        (results["t_sec"] >= spec.finish_time_lo_sec)
        & (results["t_sec"] <= spec.finish_time_hi_sec)
    ]

    # --- 5. age band, min_race_count and min_runners_per_race -------------
    # We filter iteratively until convergence because dropping an athlete can reduce a race's size
    # below the threshold, and dropping a race can reduce an athlete's race count below the threshold.
    # The entry-age band lives in the same loop: dropping an athlete's debut race can shift their
    # first race within the slice, so entry age is re-evaluated each pass.
    age_active = spec.age_lo is not None or spec.age_hi is not None
    if age_active:
        cid_to_yob = dict(zip(
            athletes["canonical_id"].to_numpy(),
            _yob_midpoint(athletes),
        ))
        rid_to_year = dict(zip(
            competitions["race_id"].to_numpy(),
            pd.to_datetime(competitions["date"]).dt.year.to_numpy(),
        ))

    while True:
        prev_len = len(results)

        # Entry-age band (first race within slice minus yob midpoint)
        if age_active:
            first_year = (
                pd.Series(results["race_id"].map(rid_to_year).to_numpy(dtype="float64"))
                  .groupby(results["canonical_id"].to_numpy())
                  .min()
            )
            yob = first_year.index.to_series().map(cid_to_yob).to_numpy(dtype="float64")
            entry_age = first_year.to_numpy() - yob
            ok = np.ones(len(entry_age), dtype=bool)
            if spec.age_lo is not None:
                ok &= (entry_age >= spec.age_lo) | np.isnan(entry_age)
            if spec.age_hi is not None:
                ok &= (entry_age <= spec.age_hi) | np.isnan(entry_age)
            keep_age = set(first_year.index.to_numpy()[ok].tolist())
            results = results[results["canonical_id"].isin(keep_age)]

        # Filter athletes
        if spec.min_race_count > 1:
            counts = results.groupby("canonical_id").size()
            keep_ath = counts[counts >= spec.min_race_count].index
            results = results[results["canonical_id"].isin(keep_ath)]
            
        # Filter races
        if spec.min_runners_per_race > 1:
            counts = results.groupby("race_id").size()
            keep_race = counts[counts >= spec.min_runners_per_race].index
            results = results[results["race_id"].isin(keep_race)]
            
        if len(results) == prev_len:
            break

    # --- 6. connectivity: keep a single connected design ------------------
    # Drop disconnected components (whole runners + races). Because components
    # are edge-disjoint, this never reduces a surviving runner's race count or a
    # surviving race's field, so the min-race/min-runner constraints above stay
    # satisfied and no re-iteration is needed.
    results, n_comp, drop_run, drop_race, drop_cell = _keep_giant_component(results)
    if n_comp > 1:
        print(f"[load_slice] WARNING: design graph split into {n_comp} connected "
              f"components; kept the largest and dropped {drop_run} runners / "
              f"{drop_race} races / {drop_cell} finishes across {n_comp - 1} "
              f"small component(s).", flush=True)

    # Align side tables with the converged results
    athletes = athletes[athletes["canonical_id"].isin(results["canonical_id"])].reset_index(drop=True)
    competitions = competitions[competitions["race_id"].isin(results["race_id"])].reset_index(drop=True)

    if len(results) == 0:
        raise ValueError(f"Slice {spec.cache_key()} produced 0 cells; loosen the filters.")

    # --- 7. re-index ------------------------------------------------------
    ath_codes, ath_uniques = pd.factorize(results["canonical_id"], sort=True)
    race_codes, race_uniques = pd.factorize(results["race_id"], sort=True)
    results = results.assign(row_idx=ath_codes, col_idx=race_codes)

    # Align per-athlete / per-race side tables to the new ordering.
    # (merge on a one-column frame so the column order from `ath_uniques`
    # is preserved as the new row order.)
    athletes_idx = (
        pd.DataFrame({"canonical_id": ath_uniques})
          .merge(athletes, on="canonical_id", how="left")
    )
    competitions_idx = (
        pd.DataFrame({"race_id": race_uniques})
          .merge(competitions, on="race_id", how="left")
    )

    I = len(athletes_idx)
    J = len(competitions_idx)
    N = len(results)

    # --- 8. A_n and A_e ---------------------------------------------------
    race_date = pd.to_datetime(competitions_idx["date"]).to_numpy()  # length J
    cell_race_date = race_date[results["col_idx"].to_numpy()]
    first_race_date = (
        pd.Series(cell_race_date)
          .groupby(results["row_idx"].to_numpy())
          .min()
          .reindex(range(I))
          .to_numpy()
    )
    A_n = (cell_race_date - first_race_date[results["row_idx"].to_numpy()]
           ).astype("timedelta64[s]").astype(np.float64) / _SEC_PER_YEAR

    # A_e: entry age = years between athlete's birth and their first race.
    # Use midpoint of yob bounds when available, NaN otherwise. Athletes
    # without yob info get NaN; downstream the fitter masks them out of the
    # gamma block.
    yob_mid = _yob_midpoint(athletes_idx)        # length I, float64, may be NaN
    first_race_year = (pd.to_datetime(first_race_date).year
                       .astype("float64").to_numpy())
    A_e = first_race_year - yob_mid

    # --- 9. response transform -------------------------------------------
    t_sec = results["t_sec"].to_numpy(dtype=np.float64)
    if spec.response == "log_time":
        y = np.log(t_sec)
    elif spec.response == "time":
        y = t_sec
    else:  # pragma: no cover - exhaustive Literal
        raise ValueError(f"unknown response transform {spec.response!r}")

    # --- 10. pack --------------------------------------------------------
    fit = FitData(
        row_idx=results["row_idx"].to_numpy(dtype=np.int64),
        col_idx=results["col_idx"].to_numpy(dtype=np.int64),
        y=y.astype(np.float64),
        A_n=A_n.astype(np.float64),
        A_e=A_e.astype(np.float64),
        athlete_ids=athletes_idx["canonical_id"].to_numpy(dtype=np.int64),
        athlete_sex=athletes_idx["sex"].to_numpy(dtype=object),
        athlete_country=athletes_idx["country"].to_numpy(dtype=object),
        race_ids=competitions_idx["race_id"].to_numpy(dtype=np.int64),
        race_series=competitions_idx["series_key"].to_numpy(dtype=object),
        race_country=competitions_idx["country"].to_numpy(dtype=object),
        race_date=race_date,
        I=I, J=J, N=N,
        spec=spec,
        data_version=data_version,
        cache_key=spec.cache_key(),
        response_kind=spec.response,
    )

    if cache_file is not None:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_file.with_suffix(".tmp")
        with open(tmp, "wb") as f:
            pickle.dump(fit, f, protocol=pickle.HIGHEST_PROTOCOL)
        tmp.replace(cache_file)  # atomic publish

    return fit


# ---------------------------------------------------------------------------
# On-disk FitData cache
# ---------------------------------------------------------------------------

def _data_fingerprint(root: Path) -> str:
    """Short hash identifying the export at `root`.

    Uses manifest.json (which records built_at / git_sha / row counts) when
    present, so re-exporting the data changes the fingerprint and invalidates
    cached slices. Falls back to parquet sizes + mtimes if there's no manifest.
    """
    manifest = root / "manifest.json"
    if manifest.exists():
        return hashlib.sha256(manifest.read_bytes()).hexdigest()[:16]
    h = hashlib.sha256()
    for name in ("athletes.parquet", "competitions.parquet", "results.parquet"):
        p = root / name
        if p.exists():
            st = p.stat()
            h.update(f"{name}:{st.st_size}:{int(st.st_mtime)}".encode())
    return h.hexdigest()[:16]


def _read_data_version(root: Path) -> str:
    """Provenance label for the export at `root` (manifest 'version' or dir name)."""
    manifest = root / "manifest.json"
    if manifest.exists():
        try:
            return str(json.loads(manifest.read_text()).get("version", root.name))
        except Exception:
            pass
    return root.name


def _cache_path(spec: SliceSpec, data_version: str, root: Path) -> Path:
    """Cache file for this (spec, data export, loader) combination."""
    payload = json.dumps({
        "spec": spec.cache_key(),
        "data_version": data_version,
        "data_fp": _data_fingerprint(root),
        "loader": LOADER_VERSION,
    }, sort_keys=True)
    key = hashlib.sha256(payload.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{data_version}_{key}.pkl"


def _keep_giant_component(results: pd.DataFrame):
    """Restrict `results` to the largest connected component of the athlete-race
    bipartite graph.

    The decomposition's athlete/race factors are only jointly identified within
    a connected component — separate components float on arbitrary relative
    scales — so the design must be a single component. Returns
    (filtered_results, n_components, n_drop_runners, n_drop_races, n_drop_cells).
    """
    cids = results["canonical_id"].to_numpy()
    rids = results["race_id"].to_numpy()
    _, ci = np.unique(cids, return_inverse=True)   # athlete-node ids 0..A-1
    ur, ri = np.unique(rids, return_inverse=True)  # race indices 0..R-1
    A, R = len(np.unique(ci)), len(ur)
    n = A + R
    g = coo_matrix((np.ones(len(ci)), (ci, A + ri)), shape=(n, n))
    n_comp, labels = connected_components(g, directed=False)
    if n_comp <= 1:
        return results, n_comp, 0, 0, 0
    cell_label = labels[ci]
    giant = np.bincount(cell_label).argmax()
    keep = cell_label == giant
    n_drop_runners = int((labels[:A] != giant).sum())
    n_drop_races = int((labels[A:] != giant).sum())
    n_drop_cells = int((~keep).sum())
    return results[keep], n_comp, n_drop_runners, n_drop_races, n_drop_cells


def clear_cache() -> int:
    """Delete all on-disk FitData cache entries. Returns the count removed."""
    if not CACHE_DIR.exists():
        return 0
    n = 0
    for f in list(CACHE_DIR.glob("*.pkl")) + list(CACHE_DIR.glob("*.tmp")):
        f.unlink()
        n += 1
    return n


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _athlete_mask(athletes: pd.DataFrame, spec: SliceSpec) -> np.ndarray:
    m = np.ones(len(athletes), dtype=bool)

    flag = _COHORT_TO_FLAG[spec.cohort]
    if flag is not None:
        m &= athletes[flag].to_numpy(dtype=bool)

    if spec.sex != "ALL":
        m &= (athletes["sex"].to_numpy() == spec.sex)

    if spec.countries_include is not None:
        m &= athletes["country"].isin(spec.countries_include).to_numpy()
    if spec.countries_exclude:
        m &= ~athletes["country"].isin(spec.countries_exclude).to_numpy()

    if spec.require_yob:
        # require *some* yob info; midpoint is computable if either bound exists
        has_any_yob = athletes[["yob_min", "yob_max"]].notna().any(axis=1).to_numpy()
        m &= has_any_yob

    return m


def _race_mask(competitions: pd.DataFrame, spec: SliceSpec) -> np.ndarray:
    m = np.ones(len(competitions), dtype=bool)

    if spec.races_include is not None:
        m &= competitions["series_key"].isin(spec.races_include).to_numpy()
    if spec.races_exclude:
        m &= ~competitions["series_key"].isin(spec.races_exclude).to_numpy()

    if spec.race_countries_include is not None:
        m &= competitions["country"].isin(spec.race_countries_include).to_numpy()
    if spec.race_countries_exclude:
        m &= ~competitions["country"].isin(spec.race_countries_exclude).to_numpy()

    if spec.date_lo is not None or spec.date_hi is not None:
        d = pd.to_datetime(competitions["date"]).dt.date.to_numpy()
        if spec.date_lo is not None:
            m &= (d >= spec.date_lo)
        if spec.date_hi is not None:
            m &= (d <= spec.date_hi)

    return m


def _pick_time_column(results: pd.DataFrame, kind: TimeKind) -> np.ndarray:
    chip = results["chip_time_sec"].to_numpy(dtype=np.float64)
    gun = results["gun_time_sec"].to_numpy(dtype=np.float64)
    if kind == "chip":
        return chip
    if kind == "gun":
        return gun
    if kind == "chip_or_gun":
        return np.where(np.isnan(chip), gun, chip)
    raise ValueError(f"unknown time_kind {kind!r}")


def _yob_midpoint(athletes: pd.DataFrame) -> np.ndarray:
    """Per-athlete year-of-birth midpoint. NaN where both bounds are missing.

    If only one bound is present we use it as-is (treated as a point estimate
    rather than throwing away the row).
    """
    lo = athletes["yob_min"].astype("Float64").to_numpy(dtype=np.float64, na_value=np.nan)
    hi = athletes["yob_max"].astype("Float64").to_numpy(dtype=np.float64, na_value=np.nan)
    both = ~np.isnan(lo) & ~np.isnan(hi)
    only_lo = ~np.isnan(lo) & np.isnan(hi)
    only_hi = np.isnan(lo) & ~np.isnan(hi)
    out = np.full(len(athletes), np.nan, dtype=np.float64)
    out[both] = 0.5 * (lo[both] + hi[both])
    out[only_lo] = lo[only_lo]
    out[only_hi] = hi[only_hi]
    return out
