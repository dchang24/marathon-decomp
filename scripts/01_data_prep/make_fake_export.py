"""
Build a **synthetic, shareable** copy of a real marathon export.

The real export under `data/race_results/` can't be published (it contains
scraped individual finish records). This script reads it and writes a perturbed
twin with an *identical schema* to `data/race_results_fake/`, so anyone can run
the analysis scripts unchanged by pointing `config.DATA_DIR` at the fake dir and
get qualitatively similar (but not real) results.

Privacy model (decision 2026-06-24, "simple"): perturb the **athlete** data
heavily, keep **race-level** data real. A dedicated party could still partially
reverse-engineer identity -- but they could do so from the public race pages
regardless -- so the trade for a trivial pipeline is accepted.

What's perturbed (seeded, reproducible):
    * finish times          multiplicative log-normal jitter (~6%)
    * yob_min / yob_max      integer jitter in {-1, 0, +1} (kept min <= max)
    * country                shared anonymization map over runner + race
                             countries -> opaque codes (AAA, AAB, ...). Keeps
                             home-race-advantage matching (one map for both).
    * master_athlete_id      overwritten with athlete_id (drops source-DB linkage)

What's preserved (race-level, real): race series_key + date + year + country
*code*, competition_id (a meaningless integer outside the private source DB, so
the covariate joins keep working untouched), ids/layout, sex, has_* flags, and
the overall athlete/race structure the decomposition depends on. No
`race_lookup.csv` is written (it maps to the private source DB).

Usage
-----
    python scripts/01_data_prep/make_fake_export.py
    python scripts/01_data_prep/make_fake_export.py --seed 7 --time-jitter 0.015
"""
from __future__ import annotations

import argparse
import json
import string
import time
from pathlib import Path

import numpy as np
import pandas as pd

from marathon_decomp.config import DATA_DIR, FAKE_DATA_DIR

PARQUETS = ("competitions.parquet", "athletes.parquet", "results.parquet")
_LETTERS = string.ascii_uppercase


def _jitter_times(s: pd.Series, rng: np.random.Generator, sigma: float) -> pd.Series:
    """Multiplicative log-normal jitter, preserving NaNs; rounded to 0.1s."""
    t = s.to_numpy(dtype="float64")
    ok = ~np.isnan(t)
    factor = np.exp(rng.normal(0.0, sigma, size=ok.sum()))
    out = t.copy()
    out[ok] = np.round(t[ok] * factor, 1)
    return pd.Series(out, index=s.index)


def _jitter_yob(athletes: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Jitter each yob bound by +/-1yr, preserving the original null pattern.

    Only enforce min <= max where both bounds are present (a one-sided bound is
    left one-sided, matching the real export's layout).
    """
    lo = athletes["yob_min"].astype("Float64").to_numpy(dtype="float64", na_value=np.nan)
    hi = athletes["yob_max"].astype("Float64").to_numpy(dtype="float64", na_value=np.nan)
    lo = lo + rng.integers(-1, 2, size=len(lo))   # NaN + int -> NaN
    hi = hi + rng.integers(-1, 2, size=len(hi))
    both = ~np.isnan(lo) & ~np.isnan(hi)
    new_lo = np.where(both, np.minimum(lo, hi), lo)
    new_hi = np.where(both, np.maximum(lo, hi), hi)
    out = athletes.copy()
    out["yob_min"] = pd.array(np.round(new_lo), dtype="Int64")
    out["yob_max"] = pd.array(np.round(new_hi), dtype="Int64")
    return out


def _country_codes(n: int) -> list[str]:
    """First `n` three-letter codes AAA, AAB, ... (base-26, distinct)."""
    codes = []
    i = 0
    while len(codes) < n:
        a, rem = divmod(i, 26 * 26)
        b, c = divmod(rem, 26)
        codes.append(_LETTERS[a] + _LETTERS[b] + _LETTERS[c])
        i += 1
    return codes


def _country_map(series_list: list[pd.Series], rng: np.random.Generator) -> dict:
    """Shared, shuffled country -> opaque-code map over the union of `series_list`.

    Codes are assigned in a seeded random order so they don't leak the original
    alphabetical ordering. NaNs are left unmapped (preserved as null downstream).
    """
    vals = pd.unique(pd.concat([s.dropna().astype(str) for s in series_list]))
    vals = sorted(vals)
    codes = _country_codes(len(vals))
    order = rng.permutation(len(vals))
    return {v: codes[order[i]] for i, v in enumerate(vals)}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--src", type=Path, default=DATA_DIR, help="Real export dir.")
    p.add_argument("--out", type=Path, default=FAKE_DATA_DIR, help="Fake export dir.")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--time-jitter", type=float, default=0.06,
                   help="Std-dev of the log-normal multiplicative time jitter "
                        "(~5-7%% perturbs finish times enough to block matching).")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()

    missing = [q for q in PARQUETS if not (args.src / q).exists()]
    if missing:
        raise SystemExit(f"source export {args.src} missing {missing}")
    existing = [q for q in PARQUETS if (args.out / q).exists()]
    if existing and not args.overwrite:
        raise SystemExit(f"output {args.out} already has {existing}. Use --overwrite.")
    args.out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    print(f"Faking export  {args.src} -> {args.out}  (seed={args.seed}, "
          f"time_jitter={args.time_jitter})")

    competitions = pd.read_parquet(args.src / "competitions.parquet")
    athletes = pd.read_parquet(args.src / "athletes.parquet")
    results = pd.read_parquet(args.src / "results.parquet")

    # competitions: keep race-level metadata real (series_key/date/year/
    # competition_id); only the country column is anonymized below.
    competitions = competitions.copy()

    # athletes: jitter yob, drop the source-DB athlete linkage.
    athletes = _jitter_yob(athletes, rng)
    athletes["master_athlete_id"] = athletes["athlete_id"]

    # country: one shared map over runner + race countries -> opaque codes.
    cmap = _country_map([athletes["country"], competitions["country"]], rng)
    athletes["country"] = athletes["country"].astype("object").map(cmap)
    competitions["country"] = competitions["country"].astype("object").map(cmap)
    print(f"  anonymized {len(cmap)} countries -> AAA codes")

    # results: jitter finish times.
    results = results.copy()
    results["chip_time_sec"] = _jitter_times(results["chip_time_sec"], rng, args.time_jitter)
    results["gun_time_sec"] = _jitter_times(results["gun_time_sec"], rng, args.time_jitter)

    competitions.to_parquet(args.out / "competitions.parquet", index=False)
    athletes.to_parquet(args.out / "athletes.parquet", index=False)
    results.to_parquet(args.out / "results.parquet", index=False)

    # Record *which* export this twin was built from, but never a machine path:
    # store it relative to the repo root (falls back to the bare dir name if the
    # source lives outside the repo). Keeps the shipped manifest free of usernames.
    _repo_root = Path(__file__).resolve().parents[2]
    try:
        _src_label = str(Path(args.src).resolve().relative_to(_repo_root)).replace("\\", "/")
    except ValueError:
        _src_label = Path(args.src).name

    manifest = {
        "synthetic": True,
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "source_export": _src_label,
        "seed": args.seed,
        "perturbation": {
            "time_jitter_lognormal_sigma": args.time_jitter,
            "yob_jitter_years": [-1, 0, 1],
            "country_anonymized": True,
            "n_countries": len(cmap),
            "dropped_source_ids": ["master_athlete_id"],
            "kept_real": ["competition_id", "series_key", "date", "year"],
        },
        "row_counts": {
            "competitions": int(len(competitions)),
            "athletes": int(len(athletes)),
            "results": int(len(results)),
        },
    }
    with open(args.out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    (args.out / "README.md").write_text(
        "# Marathon export — SYNTHETIC\n\n"
        "**This is fake data.** Generated by `scripts/01_data_prep/make_fake_export.py` "
        f"from a real export. Finish times carry ~{args.time_jitter:.1%} "
        "multiplicative jitter, yob bounds are jittered +/-1yr, country labels are "
        "anonymized to opaque codes, and the source-DB athlete linkage is dropped. "
        "Race-level metadata (series, date, year, competition_id) is real. Schema "
        "is identical to the real export so analysis scripts run unchanged by "
        "pointing `config.DATA_DIR` here. Do not interpret individual records as "
        "real.\n",
        encoding="utf-8",
    )
    print(f"Done. {len(results):,} results written to {args.out}")


if __name__ == "__main__":
    main()
