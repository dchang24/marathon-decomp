"""The named production slices + spec/slug/CLI plumbing.

Each slice is a ``SliceSpec`` at the production operating point (all SliceSpec
defaults — date window 2014..2025, min_race_count=2, min_runners_per_race=20,
chip time, log response). Only cohort and sex vary:

    {ALL, Po10} x {M (men), W (women), B (both sexes)}    -> 6 slices
    {WA}        x {M, W}                                  -> 2 slices
                                                         => SLICE_ORDER (eight)

``B`` maps to ``sex="ALL"`` (both sexes in one fit), not a separate cohort. The
WA cohort has no ``B`` (the both-sexes finish-time distribution is bimodal), so
``--slices all`` covers the eight = the six {ALL,Po10}×{W,M,B} plus WA_M, WA_W.

The readable per-slice output dir name is ``registry.slice_slug(spec)`` —
e.g. ``ALL_B_14-25_mrc2`` (cohort_gender_YY-YY_mrcN) — so the overrides below
flow straight into the results path.
"""
from __future__ import annotations

import argparse
import dataclasses
from datetime import date

from marathon_decomp import SliceSpec
from marathon_decomp import registry

# sex token -> SliceSpec.sex value
_SEX = {"M": "M", "W": "W", "B": "ALL"}

SLICES: dict[str, SliceSpec] = {
    f"{cohort}_{tok}": SliceSpec(cohort=cohort, sex=_SEX[tok])
    for cohort in ("ALL", "Po10")
    for tok in ("M", "W", "B")
}
SLICES.update({
    f"WA_{tok}": SliceSpec(cohort="WA", sex=_SEX[tok])
    for tok in ("M", "W")
})

# Stable display order — the default ``--slices all`` target. The six
# {ALL,Po10}×{W,M,B} slices plus WA_M / WA_W (WA has no B: its both-sex
# finish-time distribution is bimodal).
SLICE_ORDER = ["Po10_W", "Po10_M", "Po10_B", "ALL_W", "ALL_M", "ALL_B",
               "WA_M", "WA_W"]


def add_spec_args(ap: argparse.ArgumentParser, *, with_mrc: bool = True) -> None:
    """Register the standard slice-override flags on `ap`.

    All default to None -> fall through to the SliceSpec defaults. Pass
    ``with_mrc=False`` for scripts that define their own ``--mrc`` flag.
    """
    if with_mrc:
        ap.add_argument("--min-race-count", type=int, default=None,
                        help="athlete min race count (SliceSpec default 2).")
    ap.add_argument("--date-lo", type=int, default=None,
                    help="start year (default 2014).")
    ap.add_argument("--date-hi", type=int, default=None,
                    help="end year (default 2025).")
    ap.add_argument("--min-runner", type=int, default=None,
                    help="min finishers per race (SliceSpec default 20). Lower "
                         "for sparse cohorts, e.g. --min-runner 10 for WA.")


def build_spec(slice_name: str, *, min_race_count=None, date_lo=None,
               date_hi=None, min_runner=None) -> SliceSpec:
    """A named slice's SliceSpec with optional overrides (years -> date bounds).

    Returns the unmodified named spec when no override is given.
    """
    spec = SLICES[slice_name]
    repl: dict = {}
    if min_race_count is not None:
        repl["min_race_count"] = int(min_race_count)
    if date_lo is not None:
        repl["date_lo"] = date(int(date_lo), 1, 1)
    if date_hi is not None:
        repl["date_hi"] = date(int(date_hi), 12, 31)
    if min_runner is not None:
        repl["min_runners_per_race"] = int(min_runner)
    return dataclasses.replace(spec, **repl) if repl else spec


def slug(spec: SliceSpec) -> str:
    """Readable per-slice output dir name, e.g. ``ALL_B_14-25_mrc2``."""
    return registry.slice_slug(spec)


def resolve_names(names: list[str], ap: argparse.ArgumentParser) -> list[str]:
    """Expand ``["all"]`` to SLICE_ORDER; validate against SLICES."""
    out = SLICE_ORDER if names == ["all"] else names
    bad = [n for n in out if n not in SLICES]
    if bad:
        ap.error(f"unknown slice(s) {bad}; choose from {sorted(SLICES)} or 'all'")
    return out
