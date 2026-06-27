"""Does the aging form move the race factor v? -- recomputed from e01, no re-fit.

Reads ``grid/v_xform.parquet`` (best-init full-data v per form). At a fixed nu,
over the common race set (all forms share the slice's races, all in the
mean-v=0 gauge so directly comparable), computes:
  * pairwise across-form agreement -- Pearson + Spearman of v(form_a) vs v(form_b),
  * each form vs the canonical reference (poly2-gscalar) -- Pearson/Spearman,
    hardest-decile Jaccard, max rank shift, mean/max |dv| (log-time units).

Output -> results/model_selection/aging/{slug}/:
  * v_xform_agreement.csv  -- pairwise (form_a, form_b) corr, long form.
  * v_vs_reference.csv     -- each form vs poly2-gscalar.
Cross-slice rollups (v_vs_reference_all.csv) at the dir root.

Run::

    python scripts/02_model_selection/aging/q02_v_across_forms.py            # all slices, nu=8
    python scripts/02_model_selection/aging/q02_v_across_forms.py --nu inf --solver anderson
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "model_selection" / "aging"
REFERENCE_CAND = "poly2-gscalar"


def _spearman_cols(V: np.ndarray) -> np.ndarray:
    R = np.apply_along_axis(lambda c: pd.Series(c).rank().to_numpy(), 0, V)
    return np.corrcoef(R, rowvar=False)


def _vs_reference(v: np.ndarray, v_ref: np.ndarray) -> dict:
    J = len(v)
    k = max(1, int(round(0.10 * J)))
    top, top_ref = set(np.argsort(v)[-k:]), set(np.argsort(v_ref)[-k:])
    jacc = len(top & top_ref) / len(top | top_ref) if (top | top_ref) else float("nan")
    rank = pd.Series(v).rank().to_numpy()
    rank_ref = pd.Series(v_ref).rank().to_numpy()
    pear = float(np.corrcoef(v, v_ref)[0, 1]) if np.std(v) > 1e-12 else float("nan")
    spear = float(np.corrcoef(rank, rank_ref)[0, 1]) if np.std(rank) > 1e-12 else float("nan")
    return dict(pearson=pear, spearman=spear, top10_jaccard=jacc,
                max_abs_rank_shift_frac=float(np.max(np.abs(rank - rank_ref)) / J),
                mean_abs_dv=float(np.mean(np.abs(v - v_ref))),
                max_abs_dv=float(np.max(np.abs(v - v_ref))))


def run_slice(slug: str, nu_label: str, solver: str) -> pd.DataFrame | None:
    path = OUT_ROOT / slug / "grid" / "v_xform.parquet"
    if not path.is_file():
        print(f"  [skip] no {path}")
        return None
    df = pd.read_parquet(path)
    g = df[(df.nu == nu_label) & (df.solver == solver)]
    if g.empty:
        print(f"  [skip] {slug}: no rows for nu={nu_label} solver={solver}")
        return None
    piv = g.pivot(index="race_id", columns="cand", values="v").dropna()
    cands = list(piv.columns)
    V = piv.to_numpy()

    # pairwise agreement (long form)
    pear = np.corrcoef(V, rowvar=False)
    spear = _spearman_cols(V)
    pair_rows = [dict(slug=slug, nu=nu_label, solver=solver, form_a=cands[i],
                      form_b=cands[j], pearson=float(pear[i, j]), spearman=float(spear[i, j]))
                 for i in range(len(cands)) for j in range(len(cands))]
    pd.DataFrame(pair_rows).to_csv(OUT_ROOT / slug / "v_xform_agreement.csv", index=False)

    # vs the canonical reference
    ref_rows: list[dict] = []
    if REFERENCE_CAND in piv.columns:
        v_ref = piv[REFERENCE_CAND].to_numpy()
        for cand in cands:
            ref_rows.append(dict(slug=slug, nu=nu_label, solver=solver, cand=cand,
                                 reference=REFERENCE_CAND, **_vs_reference(piv[cand].to_numpy(), v_ref)))
    else:
        print(f"  [warn] {slug}: reference {REFERENCE_CAND} absent "
              f"(fit it in e01 for the vs-reference table)")
    ref_df = pd.DataFrame(ref_rows)
    ref_df.to_csv(OUT_ROOT / slug / "v_vs_reference.csv", index=False)
    print(f"  {slug}: {len(cands)} forms, {len(piv)} common races -> agreement + vs-ref written")
    return ref_df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", default=["all"])
    ap.add_argument("--nu", default="8", help="nu label to compare at (default 8; 'inf' for Gaussian).")
    ap.add_argument("--solver", default="anderson", choices=["als", "anderson"])
    args = ap.parse_args()

    if args.slices == ["all"]:
        slugs = sorted(p.parent.parent.name for p in OUT_ROOT.glob("*/grid/v_xform.parquet"))
    else:
        slugs = []
        for s in args.slices:
            hit = [p.parent.parent.name for p in OUT_ROOT.glob(f"{s}*/grid/v_xform.parquet")]
            slugs.extend(hit or [s])

    parts = []
    for slug in slugs:
        print(f"=== {slug} ===")
        df = run_slice(slug, args.nu, args.solver)
        if df is not None and not df.empty:
            parts.append(df)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(OUT_ROOT / "v_vs_reference_all.csv", index=False)
        print(f"\nRollup -> {OUT_ROOT / 'v_vs_reference_all.csv'}")


if __name__ == "__main__":
    main()
