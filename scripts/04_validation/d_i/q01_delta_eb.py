"""Producer 1/2: delta_j^EB -- the field career-stage composition from `full`.

For each (slice, mrc) it reads the registered ``full_nu8p00_best`` fit and forms,
per race j, the mean of the model's own per-finish drift contribution::

    delta_j^EB = mean_{i in j} ( d_i^full * an_tilde_ij )

i.e. the per-race average of ``Model._yhat_terms()["d"]``. This is the field's
net career-stage tilt as the +d model sees it. Ineligible athletes have d=0 and
so contribute 0 (consistent with the model).

CAVEAT (documented again in q03's output): because v_j^full and d_i^full come
from the *same* fit, corr(delta^EB, v^full) ~ 0 is partly mechanical -- the fit
orthogonalizes them. That near-zero is necessary but not sufficient; the
independent leave-out estimator (q02) is what makes the no-d correlation a real
finding. Use both.

Writes/updates ``results/validation/d_i/delta_eb.csv`` (upsert keyed on slug, so
running one slice then another accumulates). Argument-free runs all 8 slices x
mrc {2,5}; a (slice, mrc) without the `full` fit is skipped with a warning.

Run::

    python scripts/04_validation/d_i/q01_delta_eb.py --slices Po10_W
    python scripts/04_validation/d_i/q01_delta_eb.py                 # all 8, mrc 2 & 5
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


def delta_eb_rows(name: str, mrc: int, nu: float) -> list[dict] | None:
    spec = S.build_spec(name, min_race_count=mrc)
    slug = registry.slice_slug(spec)
    fdir = DC.full_dir(spec, nu)
    if not DC.present(fdir):
        print(f"[skip] {slug}: no full fit at {fdir}")
        return None

    fd = load_slice(spec)
    model = registry.load_fit(fdir, fd)
    terms = model._yhat_terms()
    drift_c = np.asarray(terms.get("d", np.zeros(fd.N)), float)   # d_i * an_tilde

    col = fd.col_idx
    delta = DC.per_race_mean(drift_c, col, fd.J)
    n_fin = np.bincount(col, minlength=fd.J)
    elig_obs = np.asarray(model.eligible_d, bool)[fd.row_idx].astype(float)
    n_elig = np.bincount(col, weights=elig_obs, minlength=fd.J)
    year = DC.race_year(fd)

    print(f"    {slug}: J={fd.J:,}  delta sd={delta.std():.4e}  "
          f"frac races |delta|>0={(np.abs(delta) > 1e-12).mean():.3f}  "
          f"elig finishers={int(elig_obs.sum()):,}/{fd.N:,}")

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
    print(f"delta_EB producer   nu={args.nu:g}   mrc={mrcs}   slices: {', '.join(names)}")

    rows: list[dict] = []
    n_ok = 0
    for name in names:
        for mrc in mrcs:
            r = delta_eb_rows(name, mrc, args.nu)
            if r is not None:
                rows += r
                n_ok += 1
    if not rows:
        print("\nNo (slice, mrc) had a `full` fit; nothing written.")
        return
    p = DC.merge_delta(pd.DataFrame(rows), "eb")
    print(f"\nWrote/updated {p}   ({n_ok} slice-mrc cell(s) this run)")


if __name__ == "__main__":
    main()
