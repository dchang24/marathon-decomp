"""Producer 2/2: delta_j^LOO -- field career-stage tilt, independent of any d fit.

For each (slice, mrc) it reads the registered no-d ``agingS4gv_nu8p00_best`` fit,
takes its per-finish residual ``r_ij = y_ij - yhat^{no-d}_ij``, and estimates each
athlete's raw within-career drift slope by OLS of the residual on the
within-athlete-centered career age::

    b_i = sum_j ( an_tilde_ij * r_ij ) / sum_j ( an_tilde_ij^2 )

Then the per-race composite uses the **leave-one-out** slope (race j's own
residual removed via a rank-1 downdate) so delta_j shares no information with
v_j^{no-d}::

    b_i^{(-j)} = (S_ar_i - an_ij r_ij) / (S_aa_i - an_ij^2)
    delta_j^LOO = mean_{i in j} ( b_i^{(-j)} * an_tilde_ij )

Eligible athletes have n_i >= 2 and sum(an_tilde^2) > 0 (a non-degenerate career
span); ineligible athletes contribute 0, as in the model. An observation whose
downdated denominator is non-positive (the athlete's only career spread came from
race j) also contributes 0.

NOTE on the centering: ``an_tilde`` is centered on the athlete's *full* race set
(the model's regressor); the leave-one-out touches only the slope sums, not the
centering. This removes the direct r_ij -> slope tie that would otherwise inflate
the no-d correlation; the residual O(1/n_i) effect of race j on the centering is
left as a deliberate, documented approximation.

Writes/updates ``results/validation/d_i/delta_loo.csv`` (upsert keyed on slug).
Argument-free runs all 8 slices x mrc {2,5}; a cell without the `agingS4gv` fit
is skipped with a warning.

Run::

    python scripts/04_validation/d_i/q02_delta_loo.py --slices Po10_W
    python scripts/04_validation/d_i/q02_delta_loo.py                # all 8, mrc 2 & 5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp import load_slice, registry                # noqa: E402

from baseline_common import slices as S                          # noqa: E402
import delta_common as DC                                        # noqa: E402  (sibling)


def delta_loo_rows(name: str, mrc: int, nu: float) -> list[dict] | None:
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    adir = DC.aging_dir(spec, nu)
    if not DC.present(adir):
        print(f"[skip] {slug}: no agingS4gv (no-d) fit at {adir}")
        return None

    fd = load_slice(spec)
    model = registry.load_fit(adir, fd)
    r = fd.y - DC.yhat(model)                       # no-d residual, length N
    an = DC.an_tilde(fd)                            # drift regressor, length N

    delta, elig_ath = DC.loo_delta(fd, r, an)       # shared with q04
    col = fd.col_idx
    n_fin = np.bincount(col, minlength=fd.J)
    n_elig = np.bincount(col, weights=elig_ath[fd.row_idx].astype(float), minlength=fd.J)
    year = DC.race_year(fd)

    print(f"    {slug}: J={fd.J:,}  eligible athletes={int(elig_ath.sum()):,}/{fd.I:,}  "
          f"delta sd={delta.std():.4e}  "
          f"frac races |delta|>0={(np.abs(delta) > 1e-12).mean():.3f}")

    return [dict(slug=slug, slice=name, mrc=mrc, nu=float(nu),
                 race_idx=int(j), race_id=int(fd.race_ids[j]), year=float(year[j]),
                 delta=float(delta[j]), n_finishers=int(n_fin[j]),
                 n_eligible_in_race=int(round(n_elig[j])))
            for j in range(fd.J)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--mrc", type=int, nargs="+", default=[2, 5])
    ap.add_argument("--nu", type=float, default=8.0)
    args = ap.parse_args()

    names = S.resolve_names(args.slices, ap)
    mrcs = sorted(dict.fromkeys(args.mrc))
    print(f"delta_LOO producer   nu={args.nu:g}   mrc={mrcs}   slices: {', '.join(names)}")

    rows: list[dict] = []
    n_ok = 0
    for name in names:
        for mrc in mrcs:
            rr = delta_loo_rows(name, mrc, args.nu)
            if rr is not None:
                rows += rr
                n_ok += 1
    if not rows:
        print("\nNo (slice, mrc) had an agingS4gv fit; nothing written.")
        return
    p = DC.merge_delta(pd.DataFrame(rows), "loo")
    print(f"\nWrote/updated {p}   ({n_ok} slice-mrc cell(s) this run)")


if __name__ == "__main__":
    main()
