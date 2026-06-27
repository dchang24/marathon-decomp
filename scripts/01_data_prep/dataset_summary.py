"""
dataset_summary.py
==================
Descriptive summary of the marathon modeling export for the paper's Dataset
section (paper_draft §2). Reads the read-only export at `config.DATA_DIR`
(parquets + manifest) for corpus-level stats, and calls `load_slice` for the
**analysis** counts -- i.e. what actually enters the model after min-race-count,
min-field-size and largest-connected-component pruning.

Pure descriptive: no model is fitted. Output goes to the console AND to a
markdown file (tables ready to paste into the paper), plus a per-race
supplementary CSV.

Tables produced
---------------
1. Corpus overview        - raw export + manifest (scale + identity provenance)
2. Analysis-slice grid    - load_slice over {ALL,Po10,WA,WMM} x {B,M,W} x mrc
3. Temporal coverage      - races / finishes / distinct athletes per year
4. Geography              - top race countries + athlete home countries
5. Demographics           - sex, entry age, yob/dob availability, ID coverage
6. Finish-time dist.      - mean + percentiles by sex (hh:mm:ss)
7. Sparsity & connectivity- fill %, finishes/athlete + field-size quantiles
   (+ per_race_table.csv supplementary)

Usage
-----
    python scripts/01_data_prep/dataset_summary.py            # mrc 2 and 5
    python scripts/01_data_prep/dataset_summary.py --mrc 2    # mrc2 only (fast)
"""
from __future__ import annotations

import argparse
import json
from datetime import date

import numpy as np
import pandas as pd

from marathon_decomp import SliceSpec, load_slice
from marathon_decomp.config import DATA_DIR, REFERENCE_DIR, MISC_DIR, display_path

OUT_DIR = MISC_DIR

COHORTS = ["ALL", "Po10", "WA", "WMM"]
SEXES = [("B", "ALL"), ("M", "M"), ("W", "W")]  # (label, spec value)


# ---------------------------------------------------------------------------
# Small rendering helpers (accumulate markdown while printing to console)
# ---------------------------------------------------------------------------

class Report:
    """Collects markdown sections and echoes a plain-text version to stdout."""

    def __init__(self) -> None:
        self.md: list[str] = []

    def h(self, level: int, text: str) -> None:
        self.md.append(f"\n{'#' * level} {text}\n")
        print(f"\n{'=' * 70}\n{text}\n{'=' * 70}")

    def p(self, text: str) -> None:
        self.md.append(text + "\n")
        print(text)

    def kv_table(self, rows: list[tuple[str, object]], headers=("Field", "Value")) -> None:
        self.md.append(f"| {headers[0]} | {headers[1]} |")
        self.md.append("|---|---|")
        for k, v in rows:
            self.md.append(f"| {k} | {v} |")
        self.md.append("")
        width = max(len(str(k)) for k, _ in rows)
        for k, v in rows:
            print(f"  {str(k):<{width}}  {v}")

    def table(self, df: pd.DataFrame, align_right_from: int = 1) -> None:
        cols = list(df.columns)
        self.md.append("| " + " | ".join(cols) + " |")
        sep = []
        for i in range(len(cols)):
            sep.append("---:" if i >= align_right_from else "---")
        self.md.append("| " + " | ".join(sep) + " |")
        for _, r in df.iterrows():
            self.md.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        self.md.append("")
        print(df.to_string(index=False))

    def write(self, path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(self.md), encoding="utf-8")
        print(f"\n[written] {path}")


def fmt_time(seconds) -> str:
    """Seconds -> h:mm:ss."""
    if seconds is None or (isinstance(seconds, float) and np.isnan(seconds)):
        return "-"
    s = int(round(float(seconds)))
    return f"{s // 3600:d}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def fmt_int(x) -> str:
    return f"{int(x):,}"


def fmt_pct(x, d: int = 1) -> str:
    return f"{x:.{d}f}%"


# ---------------------------------------------------------------------------
# 1. Corpus overview (raw export + manifest)
# ---------------------------------------------------------------------------

def corpus_overview(rep: Report, athletes, competitions, results, manifest) -> None:
    rep.h(2, "1. Corpus overview (raw export)")
    rep.p("One row per cluster-merged athlete (>=2 finishes), per competition, "
          "per finish. Counts are the full export *before* any modeling slice.")

    years = pd.to_datetime(competitions["date"]).dt.year
    n_ath, n_race, n_fin = len(athletes), len(competitions), len(results)
    fill = n_fin / (n_ath * n_race) * 100.0

    fpa = results.groupby("athlete_id").size()
    field = results.groupby("race_id").size()

    rows = [
        ("Window", f"{int(years.min())}-{int(years.max())} ({int(years.max()) - int(years.min()) + 1} years)"),
        ("Competitions (race editions)", fmt_int(n_race)),
        ("Distinct race series", fmt_int(competitions["series_key"].nunique())),
        ("Race countries", fmt_int(competitions["country"].nunique())),
        ("Athletes (resolved clusters, >=2 finishes)", fmt_int(n_ath)),
        ("Finishes", fmt_int(n_fin)),
        ("Matrix fill (non-missing cells)", f"{fill:.3f}%  (missing {100 - fill:.2f}%)"),
        ("Finishes / athlete (mean, median, max)",
         f"{fpa.mean():.1f}, {int(fpa.median())}, {int(fpa.max())}"),
        ("Field size / race (mean, median, min-max)",
         f"{int(field.mean())}, {int(field.median())}, {int(field.min())}-{int(field.max())}"),
    ]
    if manifest:
        idy = manifest.get("identity", {})
        cc = idy.get("conflict_clusters", {})
        # The per-rule tallies (sex/yob/same-day) are NOT disjoint: a cluster can
        # trip more than one rule but is removed once, so total_removed is the
        # union and is < sum of the per-rule counts. Report the overlap so the
        # breakdown is self-consistent (sex + yob + same_day - total = "both").
        cc_total = cc.get("total_removed", 0)
        cc_overlap = max(0, cc.get("sex", 0) + cc.get("yob", 0)
                         + cc.get("same_day", 0) - cc_total)
        rows += [
            ("Source git SHA", f"`{manifest.get('git_sha', '?')[:12]}`"),
            ("Built at", manifest.get("built_at", "?")),
            ("Surrogate athlete-records pre-merge", fmt_int(idy.get("surrogate_athletes", 0))),
            ("Clusters after union-find", fmt_int(idy.get("clusters", 0))),
            ("Conflict clusters removed (sex/yob/same-day)",
             f"{cc_total} removed "
             f"({cc.get('sex', 0)} sex, {cc.get('yob', 0)} yob, "
             f"{cc.get('same_day', 0)} same-day; {cc_overlap} tripped >1 rule, "
             f"so the per-rule counts overlap and exceed the total)"),
        ]
    rep.kv_table(rows)


def competition_catalogue(rep: Report, headline, catalogue) -> None:
    """One row per race *series* (not edition) in the headline ALL_B analysis
    slice -- i.e. the distilled dataset reported in the paper, after the >=20
    repeat-finisher cut and the largest-connected-component reduction. Columns:
    country, display name, year span, #editions, and finishers (total / men /
    women) that actually enter the model. Dumped to a CSV for the LaTeX table.
    A `reference` column is left blank as a TODO for the literature citation
    analyzing each race.
    """
    rep.h(2, "1c. Competition catalogue (ALL_B analysis series)")
    rep.p("One row per race series (pooled over editions) in the headline "
          "ALL_B slice -- the distilled dataset reported in the paper. "
          "`editions` = number of yearly editions retained; `finishers` = "
          "repeat-runner finishes entering the model (total / men / women). "
          "`reference` is a TODO for the paper citation analyzing that race.")

    fd = headline
    per_fin = pd.DataFrame({
        "series_key": np.asarray(fd.race_series)[fd.col_idx],
        "country": np.asarray(fd.race_country)[fd.col_idx],
        "year": pd.to_datetime(np.asarray(fd.race_date)[fd.col_idx]).year,
        "sex": np.asarray(fd.athlete_sex)[fd.row_idx],
    })
    g = per_fin.groupby("series_key").agg(
        country=("country", "first"),
        yr_min=("year", "min"),
        yr_max=("year", "max"),
        editions=("year", "nunique"),
        finishers=("year", "size"),
        n_men=("sex", lambda s: int((s == "M").sum())),
        n_women=("sex", lambda s: int((s == "W").sum())),
    )

    # display name from the catalogue (fall back to a de-slugged series_key)
    names = (catalogue[catalogue["distance_m"] == 42195]
             [["series_key", "name"]].drop_duplicates("series_key")
             .set_index("series_key")["name"])
    deslug = (pd.Series(g.index, index=g.index).str.replace("_marathon", "")
              .str.replace("_", " ").str.title())
    g["marathon"] = pd.Series(g.index.map(names), index=g.index).fillna(deslug)
    g["reference"] = ""  # TODO: literature citation per race
    g = g.reset_index().sort_values(["country", "marathon"])

    df = pd.DataFrame({
        "country": g["country"],
        "marathon": g["marathon"],
        "years": [f"{int(a)}-{int(b)}" if a != b else f"{int(a)}"
                  for a, b in zip(g["yr_min"], g["yr_max"])],
        "editions": g["editions"].astype(int),
        "finishers": g["finishers"].astype(int).map(fmt_int),
        "men": g["n_men"].astype(int).map(fmt_int),
        "women": g["n_women"].astype(int).map(fmt_int),
        "reference": g["reference"],
    })
    rep.table(df)
    rep.p(f"\nTotal: {len(df)} series across {g['country'].nunique()} countries, "
          f"{int(g['editions'].sum())} editions, "
          f"{fmt_int(int(g['finishers'].sum()))} finishers "
          f"({fmt_int(int(g['n_men'].sum()))} men / "
          f"{fmt_int(int(g['n_women'].sum()))} women).")

    out = OUT_DIR / "competition_catalogue.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    g[["series_key", "country", "marathon", "yr_min", "yr_max", "editions",
       "finishers", "n_men", "n_women", "reference"]].to_csv(
        out, index=False, encoding="utf-8")
    rep.p(f"\nPer-series catalogue ({len(g)} series) -> `{display_path(out)}`")


NI_BUCKETS = ["2", "3", "4", "5", "6-10", ">10"]


def _ni_breakdown(counts: np.ndarray) -> dict[str, int]:
    """Number of athletes with exactly n_i = 2,3,4,5 finishes, 6-10, and >10."""
    counts = np.asarray(counts)
    out = {str(k): int((counts == k).sum()) for k in (2, 3, 4, 5)}
    out["6-10"] = int(((counts >= 6) & (counts <= 10)).sum())
    out[">10"] = int((counts > 10).sum())
    return out


def finishes_per_athlete(rep: Report, results, headline) -> None:
    rep.h(2, "1b. Athletes by finish count (n_i)")
    rep.p("Number of athletes with exactly n_i finishes (6-10 and >10 grouped). "
          "Full export vs the headline ALL_B analysis slice.")

    full = _ni_breakdown(results.groupby("athlete_id").size().to_numpy())
    slc = _ni_breakdown(np.bincount(headline.row_idx))
    n_full, n_slc = sum(full.values()), sum(slc.values())

    records = []
    for k in NI_BUCKETS:
        records.append({
            "n_i": k,
            "athletes (export)": fmt_int(full[k]),
            "% (export)": f"{full[k] / n_full * 100:.1f}",
            "athletes (ALL_B)": fmt_int(slc[k]),
            "% (ALL_B)": f"{slc[k] / n_slc * 100:.1f}",
        })
    records.append({
        "n_i": "total",
        "athletes (export)": fmt_int(n_full), "% (export)": "100.0",
        "athletes (ALL_B)": fmt_int(n_slc), "% (ALL_B)": "100.0",
    })
    rep.table(pd.DataFrame(records))


# ---------------------------------------------------------------------------
# 2. Analysis-slice grid (via load_slice)
# ---------------------------------------------------------------------------

def slice_grid(rep: Report, mrc_list: list[int]) -> None:
    rep.h(2, "2. Analysis slices (what enters the model)")
    rep.p("Per slice, after min-race-count, min-field-size (>=20) and "
          "largest-connected-component pruning. I = athletes, J = races, "
          "N = finishes; fill = N / (I*J). B = both sexes.")

    records = []
    for mrc in mrc_list:
        for cohort in COHORTS:
            for sex_label, sex_val in SEXES:
                spec = SliceSpec(cohort=cohort, sex=sex_val, min_race_count=mrc)
                try:
                    fd = load_slice(spec)
                except Exception as e:  # sparse slices may collapse
                    records.append({
                        "slice": f"{cohort}_{sex_label}", "mrc": mrc,
                        "I": "-", "J": "-", "N": "-", "fill%": "-",
                        "fin/ath (med)": "-", "field (med)": "-",
                    })
                    print(f"  [skip] {cohort}_{sex_label} mrc{mrc}: {e}")
                    continue
                fpa = np.bincount(fd.row_idx)
                field = np.bincount(fd.col_idx)
                records.append({
                    "slice": f"{cohort}_{sex_label}", "mrc": mrc,
                    "I": fmt_int(fd.I), "J": fmt_int(fd.J), "N": fmt_int(fd.N),
                    "fill%": f"{fd.N / (fd.I * fd.J) * 100:.2f}",
                    "fin/ath (med)": int(np.median(fpa)),
                    "field (med)": int(np.median(field)),
                })
    rep.table(pd.DataFrame(records))


# ---------------------------------------------------------------------------
# 3. Temporal coverage
# ---------------------------------------------------------------------------

def temporal_coverage(rep: Report, competitions, results) -> None:
    rep.h(2, "3. Temporal coverage (per year)")
    comp = competitions[["race_id", "year"]]
    res = results.merge(comp, on="race_id", how="left")
    g = res.groupby("year")
    df = pd.DataFrame({
        "year": sorted(comp["year"].unique()),
    })
    races = comp.groupby("year")["race_id"].nunique()
    fin = g.size()
    ath = g["athlete_id"].nunique()
    df["races"] = df["year"].map(races).map(fmt_int)
    df["finishes"] = df["year"].map(fin).map(fmt_int)
    df["distinct athletes"] = df["year"].map(ath).map(fmt_int)
    df["year"] = df["year"].astype(str)
    rep.table(df)


# ---------------------------------------------------------------------------
# 4. Geography
# ---------------------------------------------------------------------------

def geography(rep: Report, athletes, competitions, results, catalogue, top=15) -> None:
    rep.h(2, f"4. Geography (top {top} race countries)")
    comp = competitions.merge(
        results.groupby("race_id").size().rename("finishes"), on="race_id", how="left"
    )
    g = comp.groupby("country").agg(
        races=("race_id", "nunique"),
        series=("series_key", "nunique"),
        finishes=("finishes", "sum"),
    ).sort_values("finishes", ascending=False)
    df = g.head(top).reset_index()
    df["races"] = df["races"].map(fmt_int)
    df["series"] = df["series"].map(fmt_int)
    df["finishes"] = df["finishes"].map(fmt_int)
    rep.table(df)
    rep.p(f"\nTotal race countries: {comp['country'].nunique()}; "
          f"series: {comp['series_key'].nunique()}.")

    rep.h(3, f"4b. Athlete home countries (top {top})")
    ah = athletes["country"].value_counts().head(top).rename("athletes").reset_index()
    ah.columns = ["country", "athletes"]
    ah["athletes"] = ah["athletes"].map(fmt_int)
    rep.table(ah)
    rep.p(f"\nTotal distinct athlete home countries: {athletes['country'].nunique()}.")


# ---------------------------------------------------------------------------
# 5. Demographics (raw export + headline ALL_B slice for entry age)
# ---------------------------------------------------------------------------

def demographics(rep: Report, athletes, headline) -> None:
    rep.h(2, "5. Demographics")

    n = len(athletes)
    sex = athletes["sex"].value_counts()
    n_m, n_w = int(sex.get("M", 0)), int(sex.get("W", 0))
    yob_known = int(athletes[["yob_min", "yob_max"]].notna().any(axis=1).sum())
    span = (athletes["yob_max"] - athletes["yob_min"]).dropna()

    rows = [
        ("Athletes", fmt_int(n)),
        ("Men", f"{fmt_int(n_m)} ({n_m / n * 100:.1f}%)"),
        ("Women", f"{fmt_int(n_w)} ({n_w / n * 100:.1f}%)"),
        ("With any year-of-birth info", f"{fmt_int(yob_known)} ({yob_known / n * 100:.1f}%)"),
        ("With known DOB (has_dob)",
         f"{fmt_int(int(athletes['has_dob'].sum()))} ({athletes['has_dob'].mean() * 100:.1f}%)"),
        ("YOB-band width (median, p90 years)",
         f"{int(span.median())}, {int(span.quantile(0.9))}" if len(span) else "-"),
    ]
    rep.p("**Sex & year-of-birth (full export)**\n")
    rep.kv_table(rows)

    rep.p("\n**External-ID coverage (full export)**\n")
    idrows = []
    for col, label in [("has_wa", "World Athletics"), ("has_po10", "Power of 10"),
                       ("has_wmm", "World Marathon Majors"), ("has_nyrr", "NYRR")]:
        c = int(athletes[col].sum())
        idrows.append((label, f"{fmt_int(c)} ({c / n * 100:.1f}%)"))
    rep.kv_table(idrows)

    # Entry age from the headline analysis slice (A_e = first-race-year - yob mid)
    ae = headline.A_e
    ae = ae[~np.isnan(ae)]
    rep.p(f"\n**Entry age (years) -- headline {headline.spec.cohort}_B, "
          f"{len(ae):,}/{headline.I:,} with known yob**\n")
    qs = [5, 25, 50, 75, 95]
    rep.kv_table(
        [("mean", f"{ae.mean():.1f}")]
        + [(f"p{q}", f"{np.percentile(ae, q):.1f}") for q in qs],
        headers=("Stat", "Age"),
    )


# ---------------------------------------------------------------------------
# 6. Finish-time distribution by sex
# ---------------------------------------------------------------------------

def finish_times(rep: Report, athletes, results) -> None:
    rep.h(2, "6. Finish-time distribution (chip time, by sex)")
    t = results["chip_time_sec"].fillna(results["gun_time_sec"])
    res = results.assign(t=t).dropna(subset=["t"])
    res = res.merge(athletes[["athlete_id", "sex"]], on="athlete_id", how="left")

    qs = [1, 10, 25, 50, 75, 90]
    records = []
    for label, mask in [("All", slice(None)), ("Men", res["sex"] == "M"),
                        ("Women", res["sex"] == "W")]:
        arr = res.loc[mask, "t"].to_numpy() if label != "All" else res["t"].to_numpy()
        row = {"group": label, "n": fmt_int(len(arr)), "mean": fmt_time(arr.mean())}
        for q in qs:
            row[f"p{q}"] = fmt_time(np.percentile(arr, q))
        records.append(row)
    rep.table(pd.DataFrame(records))


# ---------------------------------------------------------------------------
# 7. Sparsity & connectivity (headline slice) + per-race supplementary
# ---------------------------------------------------------------------------

def sparsity_and_perrace(rep: Report, headline, competitions, race_info,
                         catalogue, raw_n_races, raw_n_ath) -> None:
    rep.h(2, "7. Sparsity & connectivity (headline ALL_B)")
    fd = headline
    fpa = np.bincount(fd.row_idx)
    field = np.bincount(fd.col_idx)
    rows = [
        ("Athletes I", fmt_int(fd.I)),
        ("Races J", fmt_int(fd.J)),
        ("Finishes N", fmt_int(fd.N)),
        ("Matrix fill", f"{fd.N / (fd.I * fd.J) * 100:.3f}%  (missing {100 - fd.N / (fd.I * fd.J) * 100:.2f}%)"),
        ("Finishes/athlete (p50, p90, max)",
         f"{int(np.median(fpa))}, {int(np.percentile(fpa, 90))}, {int(fpa.max())}"),
        ("Field size (min, p10, p50, p90, max)",
         f"{int(field.min())}, {int(np.percentile(field, 10))}, {int(np.median(field))}, "
         f"{int(np.percentile(field, 90))}, {int(field.max())}"),
        ("Giant-component drop (races)", f"{raw_n_races - fd.J} of {raw_n_races}"),
        ("Giant-component / mrc drop (athletes)", f"{raw_n_ath - fd.I:,} of {raw_n_ath:,}"),
    ]
    rep.kv_table(rows)

    # Per-race supplementary CSV: raw finishers vs entering ALL_B analysis.
    cat = (catalogue[catalogue["distance_m"] == 42195]
           [["series_key", "name", "city"]].drop_duplicates("series_key"))
    analysis_field = pd.Series(field, name="n_analysis")
    per_race = pd.DataFrame({
        "race_id": fd.race_ids,
        "series_key": fd.race_series,
        "country": fd.race_country,
        "year": pd.to_datetime(fd.race_date).year,
        "n_analysis_B": field,
    }).merge(cat, on="series_key", how="left")
    # raw exported finisher counts (pre-pruning) from race_info
    per_race = per_race.merge(
        race_info[["race_id", "n_ALL_both", "n_ALL_M", "n_ALL_W"]],
        on="race_id", how="left",
    ).rename(columns={"n_ALL_both": "n_raw_B", "n_ALL_M": "n_raw_M", "n_ALL_W": "n_raw_W"})
    per_race = per_race[[
        "race_id", "name", "city", "country", "year",
        "n_raw_B", "n_raw_M", "n_raw_W", "n_analysis_B",
    ]].sort_values(["country", "name", "year"])
    out = OUT_DIR / "per_race_table.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    per_race.to_csv(out, index=False, encoding="utf-8")
    rep.p(f"\nPer-race supplementary (raw vs analysis finishers, "
          f"{len(per_race)} editions) -> `{display_path(out)}`")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5],
                    help="min-race-count values for the slice grid (default 2 5)")
    ap.add_argument("--top", type=int, default=15, help="rows in geography tables")
    args = ap.parse_args()

    print(f"Reading export from {DATA_DIR}")
    athletes = pd.read_parquet(DATA_DIR / "athletes.parquet")
    competitions = pd.read_parquet(DATA_DIR / "competitions.parquet")
    results = pd.read_parquet(DATA_DIR / "results.parquet")
    race_info = pd.read_csv(DATA_DIR / "race_info.csv")
    catalogue = pd.read_csv(REFERENCE_DIR / "race_catalogue.csv")
    manifest = {}
    mpath = DATA_DIR / "manifest.json"
    if mpath.exists():
        manifest = json.loads(mpath.read_text())

    raw_n_races, raw_n_ath = len(competitions), len(athletes)

    print("Loading headline slice ALL_B (mrc2) ...")
    headline = load_slice(SliceSpec(cohort="ALL", sex="ALL", min_race_count=2))

    rep = Report()
    rep.md.append("# Marathon dataset summary")
    rep.md.append(f"\n_Generated by `scripts/01_data_prep/dataset_summary.py` "
                  f"from export `{manifest.get('git_sha', '?')[:12]}` "
                  f"(built {manifest.get('built_at', '?')})._\n")

    corpus_overview(rep, athletes, competitions, results, manifest)
    competition_catalogue(rep, headline, catalogue)
    finishes_per_athlete(rep, results, headline)
    slice_grid(rep, args.mrc)
    temporal_coverage(rep, competitions, results)
    geography(rep, athletes, competitions, results, catalogue, top=args.top)
    demographics(rep, athletes, headline)
    finish_times(rep, athletes, results)
    sparsity_and_perrace(rep, headline, competitions, race_info, catalogue,
                         raw_n_races, raw_n_ath)

    rep.write(OUT_DIR / "dataset_summary.md")


if __name__ == "__main__":
    main()
