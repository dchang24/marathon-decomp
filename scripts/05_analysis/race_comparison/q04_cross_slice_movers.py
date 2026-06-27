"""Cross-slice "big movers": where does a race rank DIFFERENTLY between two fits?

Compares the race factor v_j between two production fits over their SHARED
races, on the four axes requested:

    M vs W   within each cohort   (A = M, B = W)
    ALL vs Po10  within each sex  (A = ALL, B = Po10)

WHY A NAIVE RANK-DIFF IS WRONG (the gauge).  Each fit is pinned independently
(mean(v)=0 + beta=0 over its OWN races), so the two slices' v_j differ by a
flat family with no likelihood cost (docs/model_derivation.md S7):

    Dv_j = v_A,j - v_B,j  ->  Dv_j + a + b*t_j + delta*t_j^2        (G0,G1,Gq)

  - G0 (level): a constant. RANKING IGNORES IT (adds the same to every race).
  - G1 + Gq (linear + quadratic in race date): a date-dependent TILT that DOES
    reorder races. A pure-convention difference here would masquerade as a pile
    of spurious rank churn.

So we strip the whole smooth-in-date subspace {1, t, t^2} from Dv over the
shared races and keep only the gauge-free residual

    r_j = Dv_j - proj_{1,t,t^2}(Dv)_j

= "race j is IDIOSYNCRATICALLY harder in A than B, net of any smooth-in-date
story". (Deliberate cost: a genuine secular A-vs-B trend is inseparable from
convention here and is discarded; what survives is race-specific, which is the
interesting signal.)  Sign: r_j > 0 => race relatively SLOWER in A than B.

BIG, NOT NOISE.  A race is flagged only if it clears BOTH bars:
    material:  |180 * r_j|  >= --min-min   (minutes at a 3:00 marathon)
    real:      z_j = |r_j| / sd_boot_j  >= --z-min
where sd_boot_j is the per-race SD of the (residualized) contrast from random-
paired draws of the two fits' INDEPENDENT bootstraps. This structurally drops
the fluctuations you don't care about (a one-rank wobble, a P(topN) jitter):
they have small r_j and z_j < 2.  We do NOT compare on P(topN) at all (it is
N-dependent and jittery).  With no bootstrap for a pair, falls back to a
magnitude-only flag (point-only, winner's-curse-prone) and warns.

RANK VIEW (for reading only).  rank_A = within-overlap fastest-rank on v_A;
rank_B = fastest-rank on the gauge-aligned v_B* = v_A - r_j (i.e. v_B re-pinned
into A's gauge), so rank_A and rank_B differ ONLY by the gauge-free r_j. d_rank
= rank_B - rank_A is reported, but interpret it only where the row is flagged.

GRANULARITY.  Series headline + edition drill-down. The per-series contrast is
the MEDIAN of its editions' r_j (robust to one bad day), with sd_boot recomputed
inside each paired replicate (preserves shared-athlete correlation). Series need
>= --min-editions overlap editions.

Outputs, one dir per comparison (results/analysis/race_comparison/cross_slice/):
    {label}/cross_movers_series.csv      ALL overlap series, sorted big->small,
                                         `flag` marks the movers
    {label}/cross_movers_editions.csv    ALL overlap editions, same
    big_movers_summary.md                flagged series across all comparisons
    + stdout. label = e.g. MvsW_ALL_14-25_mrc2 / ALLvsPo10_M_14-25_mrc2.

NO REFIT -- reads saved point fits + their bootstraps.

Run::

    python scripts/05_analysis/race_comparison/q04_cross_slice_movers.py
    python scripts/05_analysis/race_comparison/q04_cross_slice_movers.py \
        --mrcs mrc2 --z-min 2 --min-min 0.5 --min-editions 3
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from race_common import (  # noqa: E402
    OUT_ROOT, REF_MARATHON_MIN, frac_years, load_boot_wide, load_editions,
    pct, poly_resid_matrix, rank_fast,
)

CROSS_ROOT = OUT_ROOT / "cross_slice"


# --------------------------------------------------------------------------- #
def _overlap(slugA: str, slugB: str, model: str, nutag: str,
             data_version: str):
    """Inner-join two slices on race_id -> (overlap df, fit_dir_A, fit_dir_B).

    Overlap df columns: race_id, series_key, country, date, year, v_A, v_B,
    n_j_A, n_j_B. Rows with a missing date are dropped (cannot be detrended).
    """
    a = load_editions(slugA, model, nutag, data_version)
    b = load_editions(slugB, model, nutag, data_version)
    if a is None or b is None:
        return None
    dfA, fitA = a
    dfB, fitB = b
    mA = dfA.rename(columns={"v": "v_A", "n_j": "n_j_A"})
    mB = dfB[["race_id", "v", "n_j"]].rename(columns={"v": "v_B", "n_j": "n_j_B"})
    ov = mA.merge(mB, on="race_id", how="inner")
    ov = ov[ov["date"].notna()].reset_index(drop=True)
    return ov, fitA, fitB


def analyse_pair(slugA: str, slugB: str, label: str, *, model: str, nutag: str,
                 data_version: str, n_boot: int, rng, z_min: float,
                 min_min: float, min_editions: int):
    loaded = _overlap(slugA, slugB, model, nutag, data_version)
    if loaded is None:
        return None
    ov, fitA, fitB = loaded
    if len(ov) < 10:
        print(f"  [skip] {label}: only {len(ov)} shared races")
        return None

    vA = ov["v_A"].to_numpy(np.float64)
    vB = ov["v_B"].to_numpy(np.float64)
    t = frac_years(ov["date"])
    dv = vA - vB

    # strip the convention-pinned smooth-in-date gauge family from the contrast
    Mp = poly_resid_matrix(t, deg=2)
    r = dv @ Mp                                   # gauge-free residual contrast
    share_resid = float(np.var(r) / np.var(dv)) if np.var(dv) > 0 else np.nan

    # bootstrap SD of the residual contrast (random pairing of the two fits)
    bA = load_boot_wide(fitA, ov["race_id"].to_numpy(np.int64), vA)
    bB = load_boot_wide(fitB, ov["race_id"].to_numpy(np.int64), vB)
    has_boot = bA is not None and bB is not None
    if has_boot:
        iA = rng.integers(0, bA.shape[0], n_boot)
        iB = rng.integers(0, bB.shape[0], n_boot)
        r_b = (bA[iA] - bB[iB]) @ Mp              # (n_boot, n) residualized
        sd_boot = r_b.std(axis=0, ddof=1)
        n_b = int(n_boot)
    else:
        print(f"    [note] no bootstrap for {label}; magnitude-only flag")
        r_b, sd_boot, n_b = None, np.full(len(ov), np.nan), 0

    minutes = REF_MARATHON_MIN * r
    with np.errstate(invalid="ignore"):
        z = np.abs(r) / sd_boot
    flag = (np.abs(minutes) >= min_min) & (
        (z >= z_min) if has_boot else np.ones(len(ov), bool))

    # rank view: A on v_A, B on the gauge-aligned v_B* = v_A - r
    rank_A = rank_fast(vA)
    rank_B = rank_fast(vA - r)

    ed = pd.DataFrame({
        "label": label, "slice_A": slugA, "slice_B": slugB,
        "model": f"{model}_{nutag}",
        "series_key": ov["series_key"], "country": ov["country"],
        "year": ov["year"], "race_id": ov["race_id"].to_numpy(np.int64),
        "n_j_A": ov["n_j_A"], "n_j_B": ov["n_j_B"],
        "v_A": vA, "v_B": vB, "dv": dv, "r": r, "r_min_at_3h": minutes,
        "sd_boot": sd_boot, "z": z, "n_boot": n_b,
        "rank_A": rank_A, "rank_B": rank_B, "d_rank": rank_B - rank_A,
        "flag": flag,
    })
    sort_key = ed["z"] if has_boot else ed["r_min_at_3h"].abs()
    ed = ed.assign(_k=sort_key.to_numpy()).sort_values(
        "_k", ascending=False, na_position="last").drop(columns="_k"
        ).reset_index(drop=True)

    sr = _series_table(ov, r, r_b, vA, label, slugA, slugB, model, nutag,
                       has_boot, n_b, z_min, min_min, min_editions)
    return ed, sr, dict(label=label, n_overlap=len(ov), sd_dv=float(np.std(dv)),
                        share_resid=share_resid, has_boot=has_boot,
                        n_flag_ed=int(flag.sum()),
                        n_flag_sr=int(sr["flag"].sum()) if sr is not None else 0,
                        n_series=0 if sr is None else len(sr))


def _series_table(ov, r, r_b, vA, label, slugA, slugB, model, nutag,
                  has_boot, n_b, z_min, min_min, min_editions):
    """Per-series median of r_j (medians recomputed inside each replicate)."""
    series = ov["series_key"].to_numpy()
    keep = [s for s, n in zip(*np.unique(series, return_counts=True))
            if n >= min_editions]
    if len(keep) < 3:
        print(f"    [note] {label}: only {len(keep)} series with "
              f">= {min_editions} overlap editions; no series table")
        return None
    cols = {s: np.flatnonzero(series == s) for s in keep}

    m = np.array([np.median(r[cols[s]]) for s in keep])
    med_vA = np.array([np.median(vA[cols[s]]) for s in keep])
    med_vBstar = np.array([np.median((vA - r)[cols[s]]) for s in keep])

    if has_boot:
        M = np.column_stack([np.median(r_b[:, cols[s]], axis=1) for s in keep])
        sd_boot = M.std(axis=0, ddof=1)
        m_lo, m_hi = pct(M, 2.5), pct(M, 97.5)
    else:
        sd_boot = np.full(len(keep), np.nan)
        m_lo = m_hi = np.full(len(keep), np.nan)

    minutes = REF_MARATHON_MIN * m
    with np.errstate(invalid="ignore"):
        z = np.abs(m) / sd_boot
    flag = (np.abs(minutes) >= min_min) & (
        (z >= z_min) if has_boot else np.ones(len(keep), bool))

    country = {s: ov.loc[cols[s][0], "country"] for s in keep}
    sr = pd.DataFrame({
        "label": label, "slice_A": slugA, "slice_B": slugB,
        "model": f"{model}_{nutag}",
        "series_key": keep, "country": [country[s] for s in keep],
        "k_editions": [len(cols[s]) for s in keep],
        "median_vA": med_vA, "median_vBstar": med_vBstar,
        "contrast_r": m, "r_min_at_3h": minutes,
        "m_lo95": m_lo, "m_hi95": m_hi, "sd_boot": sd_boot, "z": z,
        "n_boot": n_b,
        "rank_A": rank_fast(med_vA), "rank_B": rank_fast(med_vBstar),
        "flag": flag,
    })
    sr["d_rank"] = sr["rank_B"] - sr["rank_A"]
    sort_key = sr["z"] if has_boot else sr["r_min_at_3h"].abs()
    sr = sr.assign(_k=sort_key.to_numpy()).sort_values(
        "_k", ascending=False, na_position="last").drop(columns="_k"
        ).reset_index(drop=True)
    return sr


# --------------------------------------------------------------------------- #
def _round(df: pd.DataFrame) -> pd.DataFrame:
    num = df.select_dtypes(include="number").columns
    df = df.copy()
    df[num] = df[num].round(6)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mvw-cohorts", nargs="*", default=["ALL", "Po10"],
                    help="cohorts in which to compare M vs W (A=M, B=W); '' to skip.")
    ap.add_argument("--avp-sexes", nargs="*", default=["M", "W"],
                    help="sexes in which to compare ALL vs Po10 (A=ALL, B=Po10); '' to skip.")
    ap.add_argument("--mrcs", nargs="+", default=["mrc2"])
    ap.add_argument("--window", default="14-25")
    ap.add_argument("--model", default="full", help="registry model tag (AxD = 'full').")
    ap.add_argument("--nutag", default="nu8p00")
    ap.add_argument("--data-version", default="race_results")
    ap.add_argument("--n-boot", type=int, default=4000,
                    help="random-pairing bootstrap draws for the contrast SDs.")
    ap.add_argument("--seed", type=int, default=20260614)
    ap.add_argument("--z-min", type=float, default=2.0,
                    help="flag if |r_j| / sd_boot_j >= this (beyond sampling noise).")
    ap.add_argument("--min-min", type=float, default=0.5,
                    help="flag if |180*r_j| >= this many minutes at 3:00 (material).")
    ap.add_argument("--min-editions", type=int, default=3,
                    help="series needs >= this many shared editions to rank.")
    args = ap.parse_args()

    win = args.window
    pairs: list[tuple[str, str, str]] = []
    for mrc in args.mrcs:
        for c in args.mvw_cohorts:
            pairs.append((f"{c}_M_{win}_{mrc}", f"{c}_W_{win}_{mrc}",
                          f"MvsW_{c}_{win}_{mrc}"))
        for s in args.avp_sexes:
            pairs.append((f"ALL_{s}_{win}_{mrc}", f"Po10_{s}_{win}_{mrc}",
                          f"ALLvsPo10_{s}_{win}_{mrc}"))

    rng = np.random.default_rng(args.seed)
    summaries, flagged_series = [], []
    for slugA, slugB, label in pairs:
        print(f"analysing {label}  ({slugA}  vs  {slugB}) ...")
        out = analyse_pair(slugA, slugB, label, model=args.model,
                           nutag=args.nutag, data_version=args.data_version,
                           n_boot=args.n_boot, rng=rng, z_min=args.z_min,
                           min_min=args.min_min, min_editions=args.min_editions)
        if out is None:
            continue
        ed, sr, summ = out
        summaries.append(summ)

        out_dir = CROSS_ROOT / label
        out_dir.mkdir(parents=True, exist_ok=True)
        _round(ed).to_csv(out_dir / "cross_movers_editions.csv", index=False)
        print(f"  wrote {out_dir / 'cross_movers_editions.csv'}  "
              f"({len(ed)} rows, {summ['n_flag_ed']} flagged)")
        if sr is not None:
            _round(sr).to_csv(out_dir / "cross_movers_series.csv", index=False)
            print(f"  wrote {out_dir / 'cross_movers_series.csv'}  "
                  f"({len(sr)} rows, {summ['n_flag_sr']} flagged)")
            flagged_series.append(sr[sr["flag"]])

    if not summaries:
        raise SystemExit("no slice pair produced a result")

    # ---- combined summary of flagged series movers ------------------------- #
    CROSS_ROOT.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Cross-slice big movers (gauge-free, beyond-noise series rank shifts)",
        "",
        f"r_j = (v_A - v_B) with {{1, t, t^2}} removed over the shared races =",
        "the part of the difference no gauge convention can produce; r > 0 =",
        "race relatively SLOWER in A than B. Series contrast = median of edition",
        "r_j. Flagged: |180*r| >= "
        f"{args.min_min} min at 3:00 AND z = |r|/sd_boot >= {args.z_min}.",
        "P(top-N) is deliberately NOT used (N-dependent, jittery).",
        "",
    ]
    if flagged_series:
        allflag = pd.concat(flagged_series, ignore_index=True)
        allflag = allflag.reindex(
            allflag["z"].abs().sort_values(ascending=False, na_position="last").index)
        hdr = ["comparison", "series", "ctry", "k", "contrast r",
               "min@3:00", "z", "rank A->B (d)"]
        lines += ["| " + " | ".join(hdr) + " |",
                  "|" + "|".join(["---"] * len(hdr)) + "|"]
        for _, r in allflag.iterrows():
            zr = f"{r['z']:.1f}" if pd.notna(r["z"]) else "-"
            lines.append("| " + " | ".join([
                r["label"], r["series_key"], str(r["country"]),
                f"{int(r['k_editions'])}", f"{r['contrast_r']:+.4f}",
                f"{r['r_min_at_3h']:+.2f}", zr,
                f"{int(r['rank_A'])}->{int(r['rank_B'])} ({r['d_rank']:+d})",
            ]) + " |")
    else:
        lines.append("_No series cleared both the material and noise bars._")
    (CROSS_ROOT / "big_movers_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {CROSS_ROOT / 'big_movers_summary.md'}")

    # ---- stdout recap ------------------------------------------------------ #
    print(f"\n{'comparison':>26} {'n_ov':>5} {'sd(Dv)':>7} {'res_sh':>6} "
          f"{'flag_sr':>7} {'flag_ed':>7}")
    print("-" * 62)
    for s in summaries:
        print(f"{s['label']:>26} {s['n_overlap']:>5} {s['sd_dv']:>7.4f} "
              f"{s['share_resid']:>6.2f} {s['n_flag_sr']:>7} {s['n_flag_ed']:>7}")
    print()


if __name__ == "__main__":
    main()
