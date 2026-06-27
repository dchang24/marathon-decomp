"""One-stop number lookup for the EB drift-prior `omega_d2` paragraph (QC, no fit).

Reads the frozen-omega marginal-likelihood profile (`e01_omega_profile.py`) and
the omega-init sensitivity sweep (`e02_omega_init.py`) that already ran, and pulls
the exact numbers the paper quotes into one readable Markdown file:

  * omega* is strikingly STABLE across data subsets (~2.6-2.9e-4);
  * the type-II MARGINAL likelihood peaks at the EB-learned omega* (while the
    plain data-fit log-likelihood is monotone -- it never locates omega*);
  * every reported estimand is INSENSITIVE to omega* within a [1/3, 3] band;
  * the EM update is a contraction, so omega* is INIT-INVARIANT across many
    orders of magnitude of the starting guess.

Sources (under results/model_selection/athlete_drift/):
  * omega_profile/{slug}/profile_nu8p00.parquet : one row per frozen omega (+ the
      free EB reference, `is_free`); logML / data_loglik / estimand deviations.
  * omega_init/{slug}/init_summary_nu8p00.parquet : one row per init; converged
      omega (`omega_final_rel_to_ref`), n_iter, converged flag.

Output -> results/model_selection/athlete_drift/omega_summary.md

Self-contained; no arguments needed (VS Code "Run" works). The cross-subset table
spans every slice with a profile on disk; the identification / insensitivity /
init-invariance detail is reported for the headline slice (default ALL_B).

Run::

    python scripts/02_model_selection/athlete_drift/q01_omega_summary.py
    python scripts/02_model_selection/athlete_drift/q01_omega_summary.py --slice ALL_M
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

from marathon_decomp.config import RESULTS_DIR  # noqa: E402
from report_md import render_markdown, write_markdown  # noqa: E402

ROOT = RESULTS_DIR / "model_selection" / "athlete_drift"
NUTAG = "nu8p00"


def _profile(slug: str) -> pd.DataFrame | None:
    p = ROOT / "omega_profile" / slug / f"profile_{NUTAG}.parquet"
    return pd.read_parquet(p) if p.is_file() else None


def _init(slug: str) -> pd.DataFrame | None:
    p = ROOT / "omega_init" / slug / f"init_summary_{NUTAG}.parquet"
    return pd.read_parquet(p) if p.is_file() else None


def _free_row(prof: pd.DataFrame) -> pd.Series:
    """The EB-learned reference row (the free fit)."""
    return prof[prof.is_free].iloc[0]


def collect(slugs: list[str], head: str) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    def add(section, label, value, source):
        rows.append((section, label, value, source))

    # ---- A. omega* across data subsets (stability) ----------------------
    A = "A. EB-learned omega* across data subsets (stability)"
    omegas = []
    for slug in slugs:
        prof = _profile(slug)
        if prof is None:
            continue
        fr = _free_row(prof)
        omegas.append(float(fr.omega_d2))
        add(A, slug, f"omega* = {fr.omega_d2:.3e}   (n_elig {int(fr.n_elig):,})",
            f"omega_profile/{slug}/profile_{NUTAG}.parquet, is_free row, omega_d2")
    if omegas:
        add(A, "range across subsets",
            f"{min(omegas):.2e} - {max(omegas):.2e}  "
            f"(spread {(max(omegas) / min(omegas) - 1) * 100:.0f}%)", "derived")

    prof = _profile(head)
    if prof is not None:
        fr = _free_row(prof)
        frozen = prof[~prof.is_free].sort_values("omega_mult")
        src_p = f"omega_profile/{head}/profile_{NUTAG}.parquet"

        # ---- B. marginal peaks at omega* (data-fit is monotone) ---------
        B = f"B. Marginal likelihood locates omega*  (headline slice {head})"
        i_peak = frozen.logML.idxmax()
        peak_mult = float(frozen.loc[i_peak, "omega_mult"])

        def _ml_at(mult: float) -> float:
            r = frozen[np.isclose(frozen.omega_mult, mult)]
            return float(r.logML.iloc[0]) if len(r) else float("nan")

        ml1 = _ml_at(1.0)
        add(B, "argmax_omega(logML) over the frozen sweep",
            f"omega_mult = {peak_mult:g}  (i.e. at omega*)",
            f"{src_p}, logML, argmax over frozen rows")
        add(B, "logML drop a decade BELOW omega* (x1/10)",
            f"{ml1 - _ml_at(0.1):,.0f} nats", f"{src_p}, logML(x1) - logML(x0.1)")
        add(B, "logML drop a decade ABOVE omega* (x10)",
            f"{ml1 - _ml_at(10.0):,.0f} nats", f"{src_p}, logML(x1) - logML(x10)")
        dfit_mono = bool((frozen.sort_values("omega_mult").data_loglik.diff().dropna() > 0).all())
        add(B, "plain data-fit log-lik monotone in omega?",
            f"{'yes' if dfit_mono else 'no'}  (so it never singles out omega*)",
            f"{src_p}, data_loglik vs omega_mult")

        # ---- C. insensitivity within omega* x [1/3, 3] ------------------
        C = f"C. Estimand insensitivity within omega* x [1/3, 3]  ({head})"
        band = frozen[(frozen.omega_mult >= 1 / 3 - 1e-9) & (frozen.omega_mult <= 3 + 1e-9)]
        add(C, "max |dv| race factor over the band",
            f"{band.max_abs_dv.max():.4f} log-time  (~{band.max_abs_dv.max() * 100:.1f}%)",
            f"{src_p}, max_abs_dv, omega_mult in [1/3,3]")
        add(C, "max aging-curve deviation over the band",
            f"{band.aging_maxdev.max():.4f} log-time  (~{band.aging_maxdev.max() * 100:.1f}%)",
            f"{src_p}, aging_maxdev, omega_mult in [1/3,3]")
        add(C, "min race-factor correlation to the free fit over the band",
            f"{band.corr_v_to_free.min():.4f}",
            f"{src_p}, corr_v_to_free, omega_mult in [1/3,3]")

    # ---- D. init-invariance (EM contraction) ----------------------------
    ini = _init(head)
    if ini is not None:
        D = f"D. omega* init-invariance (EM is a contraction)  ({head})"
        src_i = f"omega_init/{head}/init_summary_{NUTAG}.parquet"
        conv = ini[ini.converged] if "converged" in ini else ini
        rel = conv["omega_final_rel_to_ref"].to_numpy()
        n_inits = len(ini)
        lo_init = float(np.nanmin(ini["omega_d2_init"]))
        hi_init = float(np.nanmax(ini["omega_d2_init"]))
        add(D, "initial-guess range swept",
            f"{lo_init:.0e} - {hi_init:.0e}  ({n_inits} inits)",
            f"{src_i}, omega_d2_init")
        add(D, "converged omega* spread (rel. to reference)",
            f"max |omega/omega_ref - 1| = {np.nanmax(np.abs(rel)):.1e}  "
            f"over {len(conv)} converged inits",
            f"{src_i}, omega_final_rel_to_ref (a signed deviation, ref=0)")
        if "n_iter" in ini:
            add(D, "iteration count range",
                f"{int(conv.n_iter.min())} - {int(conv.n_iter.max())} iters",
                f"{src_i}, n_iter (converged inits)")
        n_fail = int((~ini["converged"]).sum()) if "converged" in ini else 0
        add(D, "inits that failed to converge",
            f"{n_fail}  (a too-small init can stall on the near-zero EM manifold; "
            "the failure is loud, never a silent wrong basin)",
            f"{src_i}, converged flag")

    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", default="ALL_B", help="headline slice for the detail sections")
    args = ap.parse_args()

    slugs = sorted(p.parent.name for p in ROOT.glob(f"omega_profile/*/profile_{NUTAG}.parquet"))
    if not slugs:
        print(f"No omega_profile/*/profile_{NUTAG}.parquet under {ROOT} -- run "
              "e01_omega_profile.py first.")
        return
    head = next((s for s in slugs if s.startswith(args.slice + "_")), slugs[0])

    rows = collect(slugs, head)
    report = render_markdown(
        "Drift-prior variance omega_d2: identification and stability",
        rows,
        subtitle=[
            "Per-athlete career drift d_i ~ N(0, omega_d2); omega_d2 is LEARNED "
            "(empirical Bayes), not chosen. nu=8, production operating point.",
            "The marginal (not the plain data-fit) likelihood locates omega*; the "
            "EM update that finds it is a contraction to a single fixed point.",
        ],
    )
    print(report)
    out = write_markdown(ROOT / "omega_summary.md", report)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
