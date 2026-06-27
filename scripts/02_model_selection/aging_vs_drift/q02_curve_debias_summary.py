"""One-stop number lookup for the aging-curve de-biasing paragraph (QC, no fit).

Reads the section-6 cross-mrc peak/gap table that `q01_grid_compare.py` already
wrote (`curve_debias.csv`) and renders the exact numbers the paper quotes into a
readable Markdown file.

The check: on the full (everyone) field the population aging curve is flattened by
newcomers caught mid-improvement. Two independent corrections agree -- restricting
to dedicated runners (mrc5), OR adding the per-athlete drift d_i on the full field
(mrc2) -- both restore the steeper common shape. So the `+d` full-field curve
reproduces the dedicated-runner peak and closes most of the no-d gap between the
two field cuts.

Source: results/model_selection/aging_vs_drift/curve_debias.csv
  columns: peak_{noD,withD}_mrc{2,5} (peak/improvement career-age, yr post-debut),
           gap_{noD,withD} (RMS log-time distance between the everyone and
           dedicated-runner curves), d_on_mrc{2,5} (curve shift from adding d_i).

Output -> results/model_selection/aging_vs_drift/curve_debias.md

Self-contained; no arguments needed (VS Code "Run" works).

Run::

    python scripts/02_model_selection/aging_vs_drift/q02_curve_debias_summary.py
    python scripts/02_model_selection/aging_vs_drift/q02_curve_debias_summary.py --slice ALL_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from report_md import render_markdown, write_markdown  # noqa: E402

ROOT = RESULTS_DIR / "model_selection" / "aging_vs_drift"
SRC = "curve_debias.csv"


def _pct_closed(gap_no: float, gap_with: float) -> float:
    return (gap_no - gap_with) / gap_no * 100.0 if gap_no else float("nan")


def collect(df: pd.DataFrame, head: str) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    hr = df[df.slug == head]
    if len(hr):
        r = hr.iloc[0]
        A = f"A. Peak (improvement) career-age, yr post-debut  ({head})"
        add(A, "no-d, everyone (mrc2)  [the contaminated outlier]",
            f"{r.peak_noD_mrc2:.2f} yr", f"{SRC}, peak_noD_mrc2")
        add(A, "+d, everyone (mrc2)", f"{r.peak_withD_mrc2:.2f} yr",
            f"{SRC}, peak_withD_mrc2")
        add(A, "no-d, dedicated (mrc5)", f"{r.peak_noD_mrc5:.2f} yr",
            f"{SRC}, peak_noD_mrc5")
        add(A, "+d, dedicated (mrc5)", f"{r.peak_withD_mrc5:.2f} yr",
            f"{SRC}, peak_withD_mrc5")
        add(A, "-> +d on the full field reproduces the dedicated-runner peak",
            f"{r.peak_withD_mrc2:.2f} yr  (vs dedicated {r.peak_noD_mrc5:.2f} yr)",
            "derived")

        B = f"B. Gap between everyone and dedicated-runner curves (RMS log-time)  ({head})"
        add(B, "gap, no d_i", f"{r.gap_noD:.4f}", f"{SRC}, gap_noD")
        add(B, "gap, with d_i", f"{r.gap_withD:.4f}", f"{SRC}, gap_withD")
        add(B, "-> fraction of the no-d gap closed by d_i",
            f"{_pct_closed(r.gap_noD, r.gap_withD):.0f}%", "derived")

    C = "C. Cross-subset summary (gap no-d -> with-d, % closed)"
    for _, r in df.sort_values("slug").iterrows():
        add(C, r.slug,
            f"{r.gap_noD:.4f} -> {r.gap_withD:.4f}  ({_pct_closed(r.gap_noD, r.gap_withD):.0f}% closed); "
            f"+d peak {r.peak_withD_mrc2:.2f} yr",
            f"{SRC}, gap_noD/gap_withD/peak_withD_mrc2")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default="ALL_B", help="headline slice for the detail sections")
    args = ap.parse_args()

    path = ROOT / SRC
    if not path.is_file():
        print(f"{path} not found -- run q01_grid_compare.py first.")
        return
    df = pd.read_csv(path)
    head = next((s for s in df.slug if s.startswith(args.slice + "_")), df.slug.iloc[0])

    rows = collect(df, head)
    report = render_markdown(
        "Aging-curve de-biasing: dedicated-runner restriction vs adding d_i",
        rows,
        subtitle=[
            "On the full field the population aging curve is flattened by newcomers "
            "caught mid-improvement; d_i (full field) and the dedicated-runner cut "
            "are two independent corrections that restore the steeper common shape.",
            "nu=8; mrc2 = everyone, mrc5 = dedicated runners (>=5 finishes).",
        ],
    )
    print(report)
    out = write_markdown(ROOT / "curve_debias.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
