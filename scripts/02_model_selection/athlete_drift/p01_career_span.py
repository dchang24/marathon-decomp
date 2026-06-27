"""Career-span distribution by race-count cohort -- the d_i span-floor sizer.

The per-athlete drift term d_i is only eligible for athletes with n_i >= 3 races
AND a career span >= ``d_min_span_years`` (span = max_j A_n - min_j A_n; since
A_n is years-since-debut, min A_n = 0, so span = the athlete's last A_n). Before
sweeping that floor we need to SEE the data: how is career span distributed, and
how much of the d-eligible population (n_i >= 3) would any candidate floor drop?

This script loads each production slice at the production operating point
(min_race_count=2, all other SliceSpec defaults) ONCE, computes per-athlete n_i
and career span, and draws an ECDF of span stratified by n_i cohort {2,3,4,5,6+}.
ECDF (not a density) because the decision is a *fraction-below-threshold* read:
at any candidate floor you read off directly what share of each cohort it cuts.
Log-x because the action is all in the sub-year tail where the floor bites.

    columns : gender  (M | W | B)        rows : cohort (ALL | Po10 | WA)

Vertical guides mark candidate floors (12 / 24 / 36 weeks); n_i=2 (greyed) is
the never-eligible context cohort. A console table prints, per slice, the share
of n_i>=3 athletes below each floor (-> sweep values) and the n_i>=6 count
(-> split-sample-reliability feasibility).

NOTE: this is the mrc=2 production design. An mrc=5 *re-load* would give a
slightly different span distribution (different giant component / race set), so
read this as "what the production d_i eligibility sees", not as the mrc=5 subset.

One figure -> results/model_selection/athlete_drift/career_span_mrc.png

Self-contained; VS Code "Run" works.

Run::

    python scripts/02_model_selection/athlete_drift/p01_career_span.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp import load_slice  # noqa: E402
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

from baseline_common.slices import SLICE_ORDER, build_spec  # noqa: E402

OUT_DIR = RESULTS_DIR / "model_selection" / "athlete_drift"

COHORT_ORDER = ["ALL", "Po10", "WA"]          # rows
GENDER_ORDER = ["M", "W", "B"]                # columns
_GENDER_LABEL = {"M": "men", "W": "women", "B": "both"}

# n_i strata. "2" is never d-eligible (n_i>=3 required) -> greyed context.
NI_KEYS = ["2", "3", "4", "5", "6+"]
_NI_COLOR = {
    "2": "#9e9e9e",
    "3": "#4575b4",
    "4": "#74add1",
    "5": "#f46d43",
    "6+": "#a50026",
}

WEEK = 1.0 / 52.0
FLOORS_W = [12, 24, 36]                        # candidate span floors, weeks
SAFEGUARD_YR = 1e-3                            # effective "0" floor (~half a day)


def _split_slice(name: str) -> tuple[str, str]:
    cohort, _, gender = name.rpartition("_")
    return cohort, gender


def _ni_label(n_i: np.ndarray) -> np.ndarray:
    """Vectorized n_i -> stratum key in NI_KEYS."""
    out = np.full(n_i.shape, "6+", dtype=object)
    for k in ("2", "3", "4", "5"):
        out[n_i == int(k)] = k
    return out


def _span_and_count(fd) -> tuple[np.ndarray, np.ndarray]:
    """Per-athlete (career span in years, n_i)."""
    I = fd.I
    n_i = np.bincount(fd.row_idx, minlength=I).astype(int)
    a_max = np.full(I, -np.inf)
    a_min = np.full(I, np.inf)
    np.maximum.at(a_max, fd.row_idx, fd.A_n)
    np.minimum.at(a_min, fd.row_idx, fd.A_n)
    span = np.where(np.isfinite(a_max) & np.isfinite(a_min), a_max - a_min, 0.0)
    return span, n_i


def _ecdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xs = np.sort(x)
    ys = np.arange(1, xs.size + 1) / xs.size
    return xs, ys


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- load every slice once, collect per-athlete span + n_i ----------
    data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name in SLICE_ORDER:
        try:
            fd = load_slice(build_spec(name))
        except Exception as exc:  # noqa: BLE001 -- skip unbuildable slices
            print(f"[skip] {name}: {exc}")
            continue
        data[name] = _span_and_count(fd)
        print(f"[load] {name:8s}  I={fd.I:6d}  J={fd.J:4d}  N={fd.N:7d}")

    # ---- console sweep table (the floor-sizing read) --------------------
    print("\nShare of d-eligible (n_i>=3) athletes below each candidate floor")
    print("(these are the span-sweep values); n>=6 = split-sample feasibility\n")
    hdr = (f"{'slice':10s} {'n(>=3)':>8s} " + " ".join(f'<{w}w' .rjust(7)
           for w in FLOORS_W) + f" {'n(>=6)':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for name in SLICE_ORDER:
        if name not in data:
            continue
        span, n_i = data[name]
        elig = n_i >= 3
        se = span[elig]
        n3 = se.size
        if n3 == 0:
            continue
        fracs = [f"{(se < w * WEEK).mean():7.1%}" for w in FLOORS_W]
        n6 = int((n_i >= 6).sum())
        print(f"{name:10s} {n3:8d} " + " ".join(fracs) + f" {n6:8d}")

    # ---- figure: cohort(row) x gender(col) ECDF grid --------------------
    present = {_split_slice(n) for n in data}
    cohorts = [c for c in COHORT_ORDER if any(c == cc for cc, _ in present)]
    genders = [g for g in GENDER_ORDER if any(g == gg for _, gg in present)]

    fig, axes = plt.subplots(
        len(cohorts), len(genders), figsize=(4.6 * len(genders), 3.6 * len(cohorts)),
        squeeze=False, sharex=True, sharey=True,
    )

    for r, cohort in enumerate(cohorts):
        for c, gender in enumerate(genders):
            ax = axes[r][c]
            name = f"{cohort}_{gender}"
            if name not in data:
                ax.set_visible(False)
                continue
            span, n_i = data[name]
            keys = _ni_label(n_i)
            for k in NI_KEYS:
                m = keys == k
                if m.sum() == 0:
                    continue
                xs, ys = _ecdf(span[m])
                ax.step(np.clip(xs, SAFEGUARD_YR, None), ys, where="post",
                        color=_NI_COLOR[k], lw=1.8 if k != "2" else 1.2,
                        ls="--" if k == "2" else "-",
                        label=f"n_i={k} (n={int(m.sum())})")

            for w in FLOORS_W:
                ax.axvline(w * WEEK, color="0.6", lw=0.8, ls=":")
                ax.text(w * WEEK, 0.02, f"{w}w", rotation=90, va="bottom",
                        ha="right", fontsize=7, color="0.4")
            ax.axvline(SAFEGUARD_YR, color="0.8", lw=0.8, ls=":")

            ax.set_xscale("log")
            ax.set_xlim(1e-2, 12.0)
            ax.set_ylim(0, 1)
            ax.set_title(f"{cohort} - {_GENDER_LABEL.get(gender, gender)}",
                         fontsize=11)
            ax.grid(True, which="both", alpha=0.25)
            ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
            if r == len(cohorts) - 1:
                ax.set_xlabel("career span (years, log)")
            if c == 0:
                ax.set_ylabel("cumulative share of athletes")

    fig.suptitle("Career-span ECDF by race-count cohort (mrc=2 production design)\n"
                 "dotted verticals = candidate d_i span floors; n_i=2 never "
                 "d-eligible", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = OUT_DIR / "career_span_mrc.png"
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
