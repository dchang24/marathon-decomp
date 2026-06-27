"""e01 -- leave-one-country-out (RUNNER side) sensitivity: drop, refit, save.

Does the race-difficulty ranking ``v_j`` (and the global aging block) depend on
any single nationality's *runners*? For each target country ``c`` we drop every
athlete whose ``country == c``, refit the SAME model config as a registered
production fit, and save the full result so any v_j / u_i / d_i / aging
comparison can be done later (post-fit qc/plot scripts are deferred -- this
script only produces the data they consume).

RUNNER side = exclude athletes by country (``SliceSpec.countries_exclude``).
``load_slice`` then re-runs the full downstream cascade (min_race_count,
min_runners_per_race, single connected component, re-index, A_n/A_e), so a refit
includes every race that becomes unidentifiable once ``c``'s runners leave. The
race-side sibling (drop races held in ``c``) is a later script.

DETERMINISTIC. Country exclusion has no randomness -- one refit per country, no
replicates, no seeds. Re-running a country that is already on disk is a no-op
(skipped) unless ``--overwrite`` (whole study) or ``--redo c`` (one country).

APPENDABLE BY CONSTRUCTION. Each country is its OWN registry-style fit dir:

    results/sensitivity/loco/{slice_slug}/{stem}__{key}/
        baseline/          cold refit of the FULL slice (matched-estimator anchor)
        runner_USA/        fit.pkl + manifest.json  (excl_country=USA, side=runner)
        runner_GBR/        ...
        loco_meta.json     source fit + config signature + countries done

Adding a new ``--countries`` value just creates a new ``runner_{c}/`` subdir;
nothing already computed is touched. The excluded country is recorded both
natively (``manifest.slice.countries_exclude``) and as explicit
``excl_country`` / ``side`` manifest columns.

The ``baseline/`` cold refit is the comparison anchor: the registered production
fit was warm-started, so comparing cold LOCO refits against it would
mis-attribute the cold-vs-warm gap to the country exclusion. ``loco_meta.json``
logs that gap (registered-vs-cold v Spearman + max|dv|) as a sanity check.

DISCOVERY (mirrors e10_bootstrap / mismatch e01). ``--fit`` may be a single fit
dir, a slice dir, or the ``results/models`` root; a non-existent trailing
component is matched as a prefix against siblings, so ``results/models/ALL``
selects every ``ALL_*`` slice. ``--pattern`` filters by leaf dir name.

Run::

    # smoke test: drop USA from ALL_W mrc5 production full fit
    python scripts/06_sensitivity/loco/e01_loco_runner_refit.py \
        --fit results/models/ALL_W_14-25_mrc5/full_nu8p00_best__a9bece20 \
        --countries USA
    # later: append GBR (USA untouched)
    python scripts/06_sensitivity/loco/e01_loco_runner_refit.py \
        --fit results/models/ALL_W_14-25_mrc5/full_nu8p00_best__a9bece20 \
        --countries GBR
    # or pick the N highest-volume countries in the slice automatically
    python scripts/06_sensitivity/loco/e01_loco_runner_refit.py \
        --fit results/models/ALL_W_14-25_mrc5/full_nu8p00_best__a9bece20 --top-k 10
"""
from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
import time
from dataclasses import replace
from datetime import datetime, timezone
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))      # scripts/

from marathon_decomp import (  # noqa: E402
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    load_slice,
    registry,
)
from marathon_decomp.config import RESULTS_DIR  # noqa: E402

OUT_ROOT = RESULTS_DIR / "sensitivity" / "loco"

# short solver name -> (model class, fitter-config class)
_SOLVERS: dict[str, tuple[type, type]] = {
    "als": (Model, FitterConfig),
    "anderson": (ModelAnderson, AndersonFitterConfig),
}
_NAME_TO_SHORT = {Model.NAME: "als", ModelAnderson.NAME: "anderson"}

TOL = 1e-11        # match the production bootstrap tolerance
MAX_ITER = 2000


# ---------------------------------------------------------------------------
# Discovery (mirrors scripts/06_sensitivity/mismatch_sensitivity/e01)
# ---------------------------------------------------------------------------

def _name_matches(name: str, patterns: list[str] | None) -> bool:
    if not patterns:
        return True
    for p in patterns:
        if any(ch in p for ch in "*?["):
            if fnmatch(name, p):
                return True
        elif p in name:
            return True
    return False


def resolve_roots(fit_arg: str) -> list[Path]:
    """Resolve `--fit` to existing search roots. A non-existent trailing
    component is treated as a prefix (or glob) against sibling dirs, so
    `results/models/ALL` matches every `ALL_*` slice."""
    p = Path(fit_arg)
    if p.exists():
        return [p.resolve()]
    parent, stem = p.parent, p.name
    if not parent.is_dir():
        raise SystemExit(f"--fit parent path does not exist: {parent}")
    is_glob = any(ch in stem for ch in "*?[")
    roots = sorted(
        d.resolve() for d in parent.iterdir()
        if d.is_dir() and (fnmatch(d.name, stem) if is_glob else d.name.startswith(stem))
    )
    if not roots:
        raise SystemExit(f"--fit matched no directories: {parent}/{stem}*")
    return roots


def discover_fits(root: Path, patterns: list[str] | None) -> list[Path]:
    """Every fit dir (holding a manifest.json) under `root`, matching `patterns`
    on the leaf dir name. An explicit single fit dir is never pattern-filtered."""
    if (root / "manifest.json").is_file():
        return [root]
    out = []
    for man in sorted(root.rglob("manifest.json")):
        d = man.parent
        # skip anything already living under a loco output tree
        if "loco" in d.parts:
            continue
        if _name_matches(d.name, patterns):
            out.append(d)
    return out


# ---------------------------------------------------------------------------
# Loading the source point fit
# ---------------------------------------------------------------------------

def _load_point_fit(fit_dir: Path):
    """Reload a registered point fit + its slice. Guards that saved factors still
    match the slice shape (a mismatch means the export changed -- re-fit first)."""
    man = json.loads((fit_dir / "manifest.json").read_text())
    with open(fit_dir / "fit.pkl", "rb") as f:
        payload = pickle.load(f)
    fd = load_slice(payload["spec"], payload["data_version"])
    model = registry.load_fit(fit_dir, fd)
    v = model.params.get("v")
    if v is not None and len(v) != fd.J:
        raise SystemExit(f"{fit_dir}: saved v has length {len(v)} but slice now has "
                         f"J={fd.J}; the data export changed -- re-fit first.")
    short = _NAME_TO_SHORT.get(man.get("solver"), "anderson")
    return model, fd, short, man, payload["data_version"]


# ---------------------------------------------------------------------------
# Country selection
# ---------------------------------------------------------------------------

def _country_counts(fd) -> pd.Series:
    """Per-country athlete counts in the slice, descending."""
    return (pd.Series(fd.athlete_country, dtype=object)
            .value_counts(dropna=False))


def _resolve_countries(fd, args) -> list[str]:
    counts = _country_counts(fd)
    if args.top_k:
        return [str(c) for c in counts.index[:args.top_k] if pd.notna(c)]
    return [str(c) for c in args.countries]


# ---------------------------------------------------------------------------
# Warm start
# ---------------------------------------------------------------------------

def _build_warmstart(src_params: dict, src_ath_ids, src_race_ids, tgt_fd) -> dict:
    """Map a source fit's params onto `tgt_fd`'s athlete/race ordering.

    The target athlete/race sets are subsets of the source (LOCO only drops
    runners + the races that lose their field), so every target id is present in
    the source. u_i / d_i are gathered by athlete_id; v_j by race_id; the global
    aging coeffs (theta_aging, gamma) carry over unchanged (basis dims depend
    only on the knot config). init='warmstart' consumes u, v, d, theta_aging,
    gamma; sigma2/omega_d2/nu re-estimate from there."""
    pos_ath = {int(a): i for i, a in enumerate(src_ath_ids)}
    pos_race = {int(r): j for j, r in enumerate(src_race_ids)}
    ai = np.fromiter((pos_ath[int(a)] for a in tgt_fd.athlete_ids), dtype=np.int64,
                     count=tgt_fd.I)
    rj = np.fromiter((pos_race[int(r)] for r in tgt_fd.race_ids), dtype=np.int64,
                     count=tgt_fd.J)
    warm: dict[str, np.ndarray] = {
        "u": np.asarray(src_params["u"], float)[ai],
        "v": np.asarray(src_params["v"], float)[rj],
    }
    if "d" in src_params:
        warm["d"] = np.asarray(src_params["d"], float)[ai]
    for k in ("theta_aging", "gamma"):
        if k in src_params:
            warm[k] = np.asarray(src_params[k], float).copy()
    return warm


def _atomic_rename_dir(tmp_dir: Path, out_dir: Path, *, retries: int = 8) -> None:
    """Publish `tmp_dir` as `out_dir`, retrying: OneDrive / a file indexer / AV
    can briefly hold a handle on a freshly written dir and raise PermissionError
    (WinError 5) on the rename. Back off and retry a few times before giving up."""
    for i in range(retries):
        try:
            if out_dir.exists():
                shutil.rmtree(out_dir)
            tmp_dir.rename(out_dir)
            return
        except PermissionError:
            if i == retries - 1:
                raise
            time.sleep(0.3 * (i + 1))


# ---------------------------------------------------------------------------
# Refit one slice + register
# ---------------------------------------------------------------------------

def _refit_and_register(spec, data_version, cfg, ModelCls, fcfg_cls, max_iter, tol,
                        out_dir: Path, *, excl_country: str, side: str,
                        n_excluded: int, base_J: int, warm_source=None) -> dict:
    """Load the (LOCO or full) slice, refit, register a full fit dir.

    `warm_source`, if given, is `(src_params, src_ath_ids, src_race_ids)`: the
    fit is warm-started from it (init-invariant, so it reaches the same optimum
    but converges fast). Else cold `mean` init. Writes into a temp dir and
    atomically renames on success, so a killed refit never leaves a half-written
    `runner_{c}/`. Returns a small status dict."""
    fd = load_slice(spec, data_version)
    if warm_source is not None:
        warm = _build_warmstart(warm_source[0], warm_source[1], warm_source[2], fd)
        fcfg = fcfg_cls(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                        init="warmstart", warmstart=warm, record_trace=False, verbose=0)
    else:
        fcfg = fcfg_cls(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
                        init="mean", record_trace=False, verbose=0)
    model = ModelCls(fd, cfg, fcfg)
    t0 = time.perf_counter()
    fr = model.fit()
    dt = time.perf_counter() - t0

    tmp_dir = out_dir.with_name(out_dir.name + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    registry.register_fit(
        model, tmp_dir,
        resample_tag=f"loco_{side}_{excl_country}",
        warn_on_overwrite=False,
        excl_country=excl_country, loco_side=side,
        n_excluded_runners=int(n_excluded),
        n_dropped_races=int(base_J - fd.J),
        wall_s=float(dt),
    )
    _atomic_rename_dir(tmp_dir, out_dir)
    return dict(model=model, fd=fd, fr=fr, wall_s=dt)


# ---------------------------------------------------------------------------
# Run one source fit
# ---------------------------------------------------------------------------

def run_fit(fit_dir: Path, args) -> str:
    base_model, fd, src_solver, man, data_version = _load_point_fit(fit_dir)
    solver = args.solver or src_solver
    ModelCls, fcfg_cls = _SOLVERS[solver]
    cfg = base_model.config

    key = fit_dir.name.rpartition("__")[2] or "nokey"
    stem = fit_dir.name.rpartition("__")[0] or fit_dir.name
    slug = fit_dir.parent.name
    out_dir = OUT_ROOT / slug / f"{stem}__{key}"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_path = out_dir / "loco_meta.json"

    if args.overwrite:
        for d in out_dir.iterdir():
            if d.is_dir():
                shutil.rmtree(d)

    countries = _resolve_countries(fd, args)
    counts = _country_counts(fd)

    print(f"\n=== LOCO runner refit: {slug}/{stem} (solver={solver}) ===", flush=True)
    print(f"    source I={fd.I:,} J={fd.J:,} N={fd.N:,}  "
          f"loglik={base_model.log_lik():.2f}", flush=True)
    print(f"    countries={countries}", flush=True)

    # ---- baseline: refit of the FULL slice (warm from the registered fit) -
    # Matched-estimator anchor: warm-started like every LOCO refit, so the
    # comparison isn't confounded by a cold-vs-warm gap. Same slice as the
    # source, so it lands at the registered optimum (gap logged as a check).
    base_dir = out_dir / "baseline"
    gap_spearman = gap_maxdv = float("nan")
    base_warm = (base_model.params, fd.athlete_ids, fd.race_ids)
    if not (base_dir / "manifest.json").is_file() or args.overwrite:
        st = _refit_and_register(base_model.data.spec, data_version, cfg, ModelCls,
                                 fcfg_cls, args.max_iter, args.tol, base_dir,
                                 excl_country="(none)", side="runner", n_excluded=0,
                                 base_J=fd.J, warm_source=base_warm)
        vw = np.asarray(base_model.params["v"]); vc = np.asarray(st["model"].params["v"])
        gap_spearman = float(spearmanr(vw, vc).statistic)
        gap_maxdv = float(np.abs(vw - vc).max())
        conv = "OK " if st["fr"].converged else "MAX"
        print(f"    baseline {conv} loglik={st['fr'].loglik_final:.2f} "
              f"(vs registered d={base_model.log_lik() - st['fr'].loglik_final:.2f}; "
              f"v spearman={gap_spearman:.5f} max|dv|={gap_maxdv:.4f}) "
              f"{st['wall_s']:.1f}s", flush=True)
    else:
        print("    baseline already present (skip)", flush=True)

    # Warm source for the LOCO refits = the on-disk baseline fit (its u/d/v are
    # the full-data solution; LOCO athletes/races are a subset of it).
    base_cold_model, base_fd, *_ = _load_point_fit(base_dir)
    loco_warm = (base_cold_model.params, base_fd.athlete_ids, base_fd.race_ids)

    # ---- per-country LOCO refits -----------------------------------------
    n_new = 0
    for c in countries:
        cdir = out_dir / f"runner_{c}"
        redo = args.overwrite or (args.redo and c in args.redo)
        if (cdir / "manifest.json").is_file() and not redo:
            print(f"    runner_{c}: already present (skip)", flush=True)
            continue
        n_excl = int(counts.get(c, 0))
        if n_excl == 0:
            print(f"    runner_{c}: 0 athletes with country=={c} in slice (skip)",
                  flush=True)
            continue
        loco_spec = replace(
            base_model.data.spec,
            countries_exclude=tuple(sorted(set(base_model.data.spec.countries_exclude) | {c})),
        )
        st = _refit_and_register(loco_spec, data_version, cfg, ModelCls, fcfg_cls,
                                 args.max_iter, args.tol, cdir, excl_country=c,
                                 side="runner", n_excluded=n_excl, base_J=fd.J,
                                 warm_source=loco_warm)
        lfd, fr = st["fd"], st["fr"]
        conv = "OK " if fr.converged else "MAX"
        print(f"    runner_{c}: {conv} drop {n_excl:,} runners -> "
              f"I={lfd.I:,} J={lfd.J:,} (-{fd.J - lfd.J}) N={lfd.N:,}  "
              f"loglik={fr.loglik_final:.2f} {st['wall_s']:.1f}s", flush=True)
        n_new += 1

    # ---- meta ------------------------------------------------------------
    done = sorted(d.name[len("runner_"):] for d in out_dir.glob("runner_*")
                  if (d / "manifest.json").is_file())
    meta = dict(
        source_fit_dir=str(fit_dir), source_solver=src_solver, solver=solver,
        slice_slug=slug, stem=stem, identity_key=man.get("identity_key", "?"),
        side="runner",
        source_I=fd.I, source_J=fd.J, source_N=fd.N,
        registered_loglik=float(base_model.log_lik()),
        baseline_cold_vs_registered_spearman=gap_spearman,
        baseline_cold_vs_registered_maxdv=gap_maxdv,
        max_iter=args.max_iter, tol=args.tol,
        countries_done=done,
        country_counts={str(k): int(v) for k, v in counts.head(30).items()},
        updated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )
    tmp = meta_path.with_name(meta_path.name + ".tmp")
    tmp.write_text(json.dumps(meta, indent=2))
    tmp.replace(meta_path)

    return (f"done: {slug}/{out_dir.name}  (+{n_new} country refits, "
            f"{len(done)} total)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True,
                    help="a fit dir, a slice dir, or the results/models root "
                         "(walked for fit dirs). A non-existent trailing component "
                         "is matched as a prefix against siblings.")
    ap.add_argument("--pattern", nargs="+", default=None, metavar="PAT",
                    help="select fits by leaf dir name (glob or substring); kept "
                         "if it matches ANY. Ignored when --fit is a single fit dir.")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--countries", nargs="+", default=["USA"],
                   help="ISO country codes to leave out, one refit each.")
    g.add_argument("--top-k", type=int, default=None,
                   help="instead of --countries, take the K highest-volume "
                        "countries (by athlete count) in the slice.")
    ap.add_argument("--redo", nargs="+", default=None, metavar="C",
                    help="re-run these countries even if already on disk.")
    ap.add_argument("--solver", default=None, choices=list(_SOLVERS),
                    help="refit solver (default: same as the source fit).")
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--overwrite", action="store_true",
                    help="wipe all country dirs (incl. baseline) and start fresh.")
    args = ap.parse_args()

    roots = resolve_roots(args.fit)
    fits: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        for fdir in discover_fits(root, args.pattern):
            if fdir not in seen:
                seen.add(fdir)
                fits.append(fdir)
    if not fits:
        hint = f" matching {args.pattern}" if args.pattern else ""
        where = roots[0] if len(roots) == 1 else f"{len(roots)} roots"
        raise SystemExit(f"no fits found under {where}{hint}")
    print(f"discovered {len(fits)} fit(s) under {len(roots)} root(s)")

    statuses = [run_fit(fdir, args) for fdir in fits]
    print("\n---- summary ----")
    for st in statuses:
        print(f"  {st}")


if __name__ == "__main__":
    main()
