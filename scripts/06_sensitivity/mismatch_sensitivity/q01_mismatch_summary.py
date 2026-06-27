"""q01 -- summarise how the fit moves with the mismatch rate p (S5.7).

Reads e01's per-replicate factors and, for each (operation, p), compares every
perturbed replicate's race factor v against the unperturbed cold baseline (run 0)
on the **common race set** (races surviving in both -- break/join can prune small
races). Reports ranking agreement, |delta v| magnitude, and a bias-vs-variance
split, then -- if a bootstrap group (results/models/.../bootstrap) for the same
fit is supplied -- the headline **sensitivity ratio** = mismatch per-race SD /
bootstrap per-race SD, plus the interpolated p* where that ratio crosses 1 (the
error rate at which identity error overtakes ordinary sampling noise).

Sign convention: rank 1 = hardest (largest v), matching the bootstrap.
delta v is reported in log-time and as an approximate minutes-at-3:00 figure
(delta minutes ~= 180 * delta v).

Outputs under the e01 group dir, in a ``summary/`` subdir:
  summary.parquet   one row per (op, p): replicate-distribution of rank/dv
                    metrics, record_frac, and (optional) sensitivity ratio.
  summary.md        the same headline table as readable Markdown (the paper cites
                    this) -- sensitivity ratio + rank stability per (op, p) + p*.
  per_race.parquet  one row per (op, p, race_id): v_point, dv mean (bias) & SD.
  meta.json         p* crossings + provenance.

Run::
    python scripts/06_sensitivity/mismatch_sensitivity/q01_mismatch_summary.py --group results/sensitivity/mismatch_sensitivity/Po10_B_14-25_mrc2/full_nu8p00_best__7cde3824 --bootstrap results/models/Po10_B_14-25_mrc2/full_nu8p00_best__7cde3824/bootstrap

    python scripts/06_sensitivity/mismatch_sensitivity/q01_mismatch_summary.py --group results/sensitivity/mismatch_sensitivity/ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b --bootstrap results/models/ALL_B_14-25_mrc2/full_nu8p00_best__c6a5e58b/bootstrap
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, rankdata, spearmanr

LO, HI = 0.025, 0.975
MIN_PER_3H = 180.0       # ref pace: delta minutes ~= MIN_PER_3H * delta(log time)


def _summ(a: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return dict(median=np.nan, lo=np.nan, hi=np.nan)
    return dict(median=float(np.median(a)), lo=float(np.quantile(a, LO)),
                hi=float(np.quantile(a, HI)))


def _bootstrap_race_sd(boot_dir: Path) -> pd.Series | None:
    """Per-race SD of v across the bootstrap replicates (run_id>0)."""
    rf = boot_dir / "race_factors.parquet"
    if not rf.is_file():
        return None
    df = pd.read_parquet(rf, columns=["run_id", "race_id", "v"])
    df = df[df.run_id > 0]
    return df.groupby("race_id")["v"].std(ddof=1)


def compute(group_dir: Path, boot_dir: Path | None):
    rf = pd.read_parquet(group_dir / "race_factors.parquet")
    runs = pd.read_parquet(group_dir / "runs.parquet")

    base = rf[rf.op == "none"][["race_id", "v"]].set_index("race_id")["v"]
    pert = rf[rf.op != "none"]
    # Realised rates per (op, p): the record-level reassignment fraction and the
    # realised *per-athlete* rate (the honest counterpart of the requested knob p
    # after the join gating; a merge touches 2 athletes -> 2*n_join_merges).
    runs_p = runs[runs.op != "none"].copy()
    runs_p["real_ath"] = ((runs_p["n_break_splits"] + 2 * runs_p["n_join_merges"])
                          / runs_p["n_athletes_base"])
    rate = (runs_p.groupby(["op", "p"])
            .agg(record_frac_mean=("record_frac", "mean"),
                 real_ath_mean=("real_ath", "mean")))

    boot_sd = _bootstrap_race_sd(boot_dir) if boot_dir is not None else None

    summary_rows: list[dict] = []
    per_race_rows: list[dict] = []
    for (op, p), grp in pert.groupby(["op", "p"]):
        reps = grp.run_id.unique()
        spear, kend, pear, jac, shift, medabs, rms = ([] for _ in range(7))
        # accumulate per-race delta v across replicates for bias/variance
        dv_acc: dict[int, list[float]] = {}
        for r in reps:
            vr = grp[grp.run_id == r][["race_id", "v"]].set_index("race_id")["v"]
            common = base.index.intersection(vr.index)
            if len(common) < 5:
                continue
            v0 = base.loc[common].to_numpy()
            vp = vr.loc[common].to_numpy()
            spear.append(spearmanr(v0, vp).statistic)
            kend.append(kendalltau(v0, vp).statistic)
            pear.append(np.corrcoef(v0, vp)[0, 1])
            r0 = rankdata(-v0); rp = rankdata(-vp)
            shift.append(float(np.abs(r0 - rp).max() / len(common)))
            k = max(1, round(0.10 * len(common)))
            top0 = set(common[np.argsort(r0)[:k]]); topp = set(common[np.argsort(rp)[:k]])
            jac.append(len(top0 & topp) / len(top0 | topp))
            dv = vp - v0
            medabs.append(float(np.median(np.abs(dv))))
            rms.append(float(np.sqrt(np.mean(dv ** 2))))
            for rid, d in zip(common.to_numpy(), dv):
                dv_acc.setdefault(int(rid), []).append(float(d))

        # per-race bias (mean delta v) and variance (SD delta v) across replicates
        ratios, biases = [], []
        for rid, ds in dv_acc.items():
            ds = np.asarray(ds)
            sd = float(ds.std(ddof=1)) if len(ds) > 1 else np.nan
            biases.append(abs(float(ds.mean())))
            bsd = (float(boot_sd[rid]) if boot_sd is not None and rid in boot_sd.index
                   else np.nan)
            ratio = (sd / bsd if np.isfinite(sd) and np.isfinite(bsd) and bsd > 0
                     else np.nan)
            per_race_rows.append(dict(op=op, p=p, race_id=rid, v_point=float(base.loc[rid]),
                                      dv_mean=float(ds.mean()), dv_sd=sd, boot_sd=bsd,
                                      sens_ratio=ratio, n_rep=len(ds)))
            if np.isfinite(ratio):
                ratios.append(ratio)

        row = dict(
            op=op, p=float(p), n_rep=len(reps),
            **{f"spearman_{k}": v for k, v in _summ(spear).items()},
            **{f"kendall_{k}": v for k, v in _summ(kend).items()},
            pearson_median=_summ(pear)["median"],
            top10_jaccard_median=_summ(jac)["median"],
            max_rank_shift_median=_summ(shift)["median"],
            median_abs_dv_log=_summ(medabs)["median"],
            median_abs_dv_min=_summ(medabs)["median"] * MIN_PER_3H,
            rms_dv_log=_summ(rms)["median"],
            mean_abs_bias_log=(float(np.mean(biases)) if biases else np.nan),
            sens_ratio_median=(float(np.median(ratios)) if ratios else np.nan),
            sens_ratio_p75=(float(np.quantile(ratios, 0.75)) if ratios else np.nan),
            sens_ratio_max=(float(np.max(ratios)) if ratios else np.nan),
            n_ratio_races=len(ratios),
        )
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows).merge(
        rate.reset_index(), on=["op", "p"], how="left").sort_values(["op", "p"])
    per_race = pd.DataFrame(per_race_rows)

    # p* where the median sensitivity ratio first crosses 1, per op (linear interp)
    pstar = {}
    if boot_dir is not None:
        for op, g in summary.groupby("op"):
            g = g.sort_values("p")
            x = g.p.to_numpy(); y = g.sens_ratio_median.to_numpy()
            cross = np.where((y[:-1] < 1) & (y[1:] >= 1))[0]
            if cross.size:
                i = cross[0]
                t = (1 - y[i]) / (y[i + 1] - y[i]) if y[i + 1] != y[i] else 0.0
                pstar[op] = float(x[i] + t * (x[i + 1] - x[i]))
            else:
                pstar[op] = float("nan") if (y >= 1).any() else float("inf")
    return summary, per_race, pstar


def _pstar_str(pstar: dict) -> str:
    return ", ".join(
        f"{op} = {('inf' if np.isinf(v) else ('n/a' if not np.isfinite(v) else f'{v:.3f}'))}"
        for op, v in pstar.items()) or "n/a (no bootstrap supplied)"


def render_markdown(group_dir: Path, summary: pd.DataFrame, pstar: dict) -> str:
    """Render the headline (op, p) table as a readable GFM document.

    Built by hand (no tabulate dependency). Mirrors the loco/loso q01 summaries.
    """
    # (column, header, format) for the headline table.
    spec = [
        ("op", "op", lambda x: str(x)),
        ("p", "p (requested)", lambda x: f"{x:.3f}"),
        ("real_ath_mean", "realised per-athlete rate", lambda x: f"{x:.4f}"),
        ("spearman_median", "spearman(v, v0)", lambda x: f"{x:.4f}"),
        ("top10_jaccard_median", "top-10% Jaccard", lambda x: f"{x:.3f}"),
        ("median_abs_dv_min", "med abs(dv) min@3:00", lambda x: f"{x:.3f}"),
        ("sens_ratio_median", "sens_ratio med", lambda x: f"{x:.3f}"),
        ("sens_ratio_p75", "sens_ratio p75", lambda x: f"{x:.3f}"),
        ("sens_ratio_max", "sens_ratio max", lambda x: f"{x:.3f}"),
    ]
    head = "| " + " | ".join(h for _, h, _ in spec) + " |"
    sep = "|" + "|".join("---" for _ in spec) + "|"
    body = []
    for _, r in summary.sort_values(["op", "p"]).iterrows():
        body.append("| " + " | ".join(f(r[c]) for c, _, f in spec) + " |")

    lines = [
        f"# Identity-mismatch sensitivity -- {group_dir.parent.name}/{group_dir.name}",
        "",
        "Inject ADDITIONAL linkage error on top of the production resolution -- "
        "`op` = `break` (split one runner's records into two, mimicking a failed "
        "join), `join` (merge two different runners, mimicking a wrong match), or "
        "`both` -- at requested rate `p`, refit under Monte-Carlo replication, and "
        "compare the perturbed race factor `v_j` to the unperturbed baseline on the "
        "common race set (rank 1 = hardest).",
        "",
        "**sens_ratio** = mismatch per-race SD / bootstrap per-race SD: the "
        "mismatch-induced movement relative to the athlete-sampling noise the "
        "analysis already carries. `sens_ratio < 1` means identity error perturbs "
        "`v_j` by LESS than ordinary sampling noise. **p\\*** = interpolated rate "
        "where the median ratio crosses 1 (`inf` = never, within the tested range). "
        "`|dv|` is in minutes at a 3:00:00 marathon (`180 * delta log-time`).",
        "",
        f"**p\\* (rate where identity error = sampling noise):** {_pstar_str(pstar)}",
        "",
        head, sep, *body, "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", required=True, help="an e01 mismatch group dir.")
    ap.add_argument("--bootstrap", default=None,
                    help="matching bootstrap dir (results/models/.../bootstrap) for "
                         "the sensitivity ratio.")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    group_dir = Path(args.group).resolve()
    boot_dir = Path(args.bootstrap).resolve() if args.bootstrap else None
    out_dir = group_dir / "summary"
    if (out_dir / "summary.parquet").is_file() and not args.overwrite:
        print(f"skip (done): {out_dir} -- use --overwrite")
        return

    summary, per_race, pstar = compute(group_dir, boot_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_parquet(out_dir / "summary.parquet", index=False)
    (out_dir / "summary.md").write_text(
        render_markdown(group_dir, summary, pstar), encoding="utf-8")
    per_race.to_parquet(out_dir / "per_race.parquet", index=False)
    (out_dir / "meta.json").write_text(json.dumps(dict(
        source_group=str(group_dir), bootstrap_dir=(str(boot_dir) if boot_dir else None),
        p_star=pstar, min_per_3h=MIN_PER_3H,
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds")), indent=2))

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        cols = ["op", "p", "real_ath_mean", "record_frac_mean", "spearman_median",
                "median_abs_dv_min", "sens_ratio_median", "sens_ratio_p75",
                "sens_ratio_max"]
        print(summary[cols].to_string(index=False))
    if pstar:
        print("\np* (rate where identity error = sampling noise):",
              {k: (f"{v:.3f}" if np.isfinite(v) else str(v)) for k, v in pstar.items()})


if __name__ == "__main__":
    main()
