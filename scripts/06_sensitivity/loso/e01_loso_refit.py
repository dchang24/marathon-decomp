"""e01 -- leave-one-series-out (LOSO) sensitivity: drop a venue, refit, save.

How much does the rest of the race-difficulty ranking move when one whole race
SERIES (all editions of a venue -- e.g. every Boston, every Valencia) is pulled
out of the network? This is the surgical, race-side companion to the runner-side
country LOCO: it isolates a single venue's leverage instead of bundling a whole
country's races. The removed races are gone, so the probe is the SPILLOVER onto
the SURVIVING races' v_j (plus u_i / d_i / aging). Post-fit qc/plot scripts are
deferred -- this script only produces the full fits they consume.

EXCLUSION = drop a series by ``SliceSpec.races_exclude=(series_key,)``.
``load_slice`` then re-runs the full cascade: athletes who only ran that venue
fall below min_race_count, the design is re-reduced to its giant component, and
A_n/A_e recompute. So a refit honestly includes the connectivity loss of pulling
a hub out (n_dropped_races can EXCEED the venue's own editions when removal
strands other races -- itself a leverage signal).

DETERMINISTIC. One refit per series, no replicates/seeds. Re-running a series
already on disk is a no-op unless ``--overwrite`` (whole study) or ``--redo s``.

APPENDABLE BY CONSTRUCTION. Each series is its OWN registry-style fit dir:

    results/sensitivity/loso/{slice_slug}/{stem}__{key}/
        baseline/          warm refit of the FULL slice (matched-estimator anchor)
        series_boston/     fit.pkl + manifest.json  (excl_series=boston)
        series_valencia/   ...
        loso_meta.json     source fit + config signature + series done + field table

Series labels drop a trailing ``_marathon`` for brevity: ``boston_marathon`` ->
dir ``series_boston``, resample_tag ``loso_boston``. The FULL series_key is
recorded natively in ``manifest.slice.races_exclude`` and as the explicit
``series_key`` column; the short label as ``excl_series``. ``--series`` accepts
either form. Adding a new ``--series`` value just creates a new subdir; nothing
already computed is touched.

Every refit is WARM-STARTED from the baseline (surviving v by race_id, surviving
u/d by athlete_id; aging carries over). The optimum is init-invariant, so this
reaches the same place a cold fit would but converges fast.

DISCOVERY (mirrors e10_bootstrap / loco e01). ``--fit`` may be a single fit dir,
a slice dir, or the ``results/models`` root; a non-existent trailing component is
matched as a prefix against siblings. ``--pattern`` filters by leaf dir name.

Run::

    # test three venues on ALL_B mrc5
    python scripts/06_sensitivity/loso/e01_loso_refit.py \
        --fit results/models/ALL_B_14-25_mrc5/full_nu8p00_best__02ba1954 \
        --series boston london valencia
    # or the K highest-field venues automatically
    python scripts/06_sensitivity/loso/e01_loso_refit.py \
        --fit results/models/ALL_B_14-25_mrc5/full_nu8p00_best__02ba1954 --top-k 8
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

OUT_ROOT = RESULTS_DIR / "sensitivity" / "loso"

_SOLVERS: dict[str, tuple[type, type]] = {
    "als": (Model, FitterConfig),
    "anderson": (ModelAnderson, AndersonFitterConfig),
}
_NAME_TO_SHORT = {Model.NAME: "als", ModelAnderson.NAME: "anderson"}

TOL = 1e-11        # match the production bootstrap tolerance
MAX_ITER = 2000
_SUFFIX = "_marathon"


# ---------------------------------------------------------------------------
# Discovery (mirrors loco e01)
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
        if "loso" in d.parts:          # skip our own output tree
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
# Series selection / naming
# ---------------------------------------------------------------------------

def _short_series(key: str) -> str:
    """Concise label for a series_key: drop a trailing '_marathon'."""
    return key[:-len(_SUFFIX)] if key.endswith(_SUFFIX) else key


def _series_field(fd) -> pd.Series:
    """Total field (sum of edition field sizes) per series_key, descending."""
    field = np.bincount(fd.col_idx, minlength=fd.J)
    return (pd.Series(field, index=pd.Series(fd.race_series, dtype=object))
            .groupby(level=0).sum().sort_values(ascending=False))


def _resolve_series(fd, args) -> list[str]:
    """Return a list of full series_keys to drop. `--top-k` ranks by total field;
    otherwise `--series` tokens are matched as a full key, a short label
    (`_marathon` stripped), or token+`_marathon`."""
    field = _series_field(fd)
    full_keys = list(field.index)
    if args.top_k:
        return full_keys[:args.top_k]
    short_to_full: dict[str, str] = {_short_series(k): k for k in full_keys}
    full_set = set(full_keys)
    out: list[str] = []
    for tok in args.series:
        if tok in full_set:
            out.append(tok)
        elif tok in short_to_full:
            out.append(short_to_full[tok])
        elif tok + _SUFFIX in full_set:
            out.append(tok + _SUFFIX)
        else:
            avail = ", ".join(sorted(short_to_full)[:40])
            raise SystemExit(f"--series '{tok}' not in slice. Available labels: {avail} ...")
    return out


# ---------------------------------------------------------------------------
# Warm start
# ---------------------------------------------------------------------------

def _build_warmstart(src_params: dict, src_ath_ids, src_race_ids, tgt_fd) -> dict:
    """Map a source fit's params onto `tgt_fd`'s athlete/race ordering.

    The target athlete/race sets are subsets of the source (LOSO drops a venue's
    races + the athletes/races stranded by it), so every target id is present in
    the source. u_i / d_i are gathered by athlete_id; v_j by race_id; the global
    aging coeffs carry over unchanged. init='warmstart' consumes u, v, d,
    theta_aging, gamma; sigma2/omega_d2/nu re-estimate from there."""
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
                        out_dir: Path, *, excl_series: str, series_key: str,
                        n_editions: int, base_I: int, base_J: int,
                        warm_source=None) -> dict:
    """Load the (LOSO or full) slice, refit, register a full fit dir.

    `warm_source`, if given, is `(src_params, src_ath_ids, src_race_ids)`: the
    fit warm-starts from it. Writes into a temp dir and atomically renames on
    success, so a killed refit never leaves a half-written `series_{label}/`."""
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
        resample_tag=f"loso_{excl_series}",
        warn_on_overwrite=False,
        excl_series=excl_series, series_key=series_key, loso_unit="series",
        n_editions=int(n_editions),
        n_dropped_races=int(base_J - fd.J),         # >= n_editions if removal strands races
        n_dropped_runners=int(base_I - fd.I),
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
    meta_path = out_dir / "loso_meta.json"

    if args.overwrite:
        for d in out_dir.iterdir():
            if d.is_dir():
                shutil.rmtree(d)

    field = _series_field(fd)
    series_keys = _resolve_series(fd, args)
    n_ed = pd.Series(fd.race_series, dtype=object).value_counts()

    print(f"\n=== LOSO refit: {slug}/{stem} (solver={solver}) ===", flush=True)
    print(f"    source I={fd.I:,} J={fd.J:,} N={fd.N:,}  "
          f"loglik={base_model.log_lik():.2f}", flush=True)
    print(f"    series={[_short_series(k) for k in series_keys]}", flush=True)

    # ---- baseline: warm refit of the FULL slice (matched anchor) ----------
    base_dir = out_dir / "baseline"
    gap_spearman = gap_maxdv = float("nan")
    base_warm = (base_model.params, fd.athlete_ids, fd.race_ids)
    if not (base_dir / "manifest.json").is_file() or args.overwrite:
        st = _refit_and_register(base_model.data.spec, data_version, cfg, ModelCls,
                                 fcfg_cls, args.max_iter, args.tol, base_dir,
                                 excl_series="(none)", series_key="(none)",
                                 n_editions=0, base_I=fd.I, base_J=fd.J,
                                 warm_source=base_warm)
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

    # Warm source for the LOSO refits = the on-disk baseline fit.
    base_cold_model, base_fd, *_ = _load_point_fit(base_dir)
    loso_warm = (base_cold_model.params, base_fd.athlete_ids, base_fd.race_ids)

    # ---- per-series LOSO refits ------------------------------------------
    n_new = 0
    for fullkey in series_keys:
        label = _short_series(fullkey)
        sdir = out_dir / f"series_{label}"
        redo = args.overwrite or (args.redo and label in args.redo)
        if (sdir / "manifest.json").is_file() and not redo:
            print(f"    series_{label}: already present (skip)", flush=True)
            continue
        n_editions = int(n_ed.get(fullkey, 0))
        if n_editions == 0:
            print(f"    series_{label}: 0 races for '{fullkey}' in slice (skip)", flush=True)
            continue
        loso_spec = replace(
            base_model.data.spec,
            races_exclude=tuple(sorted(set(base_model.data.spec.races_exclude) | {fullkey})),
        )
        st = _refit_and_register(loso_spec, data_version, cfg, ModelCls, fcfg_cls,
                                 args.max_iter, args.tol, sdir, excl_series=label,
                                 series_key=fullkey, n_editions=n_editions,
                                 base_I=fd.I, base_J=fd.J, warm_source=loso_warm)
        lfd, fr = st["fd"], st["fr"]
        conv = "OK " if fr.converged else "MAX"
        print(f"    series_{label}: {conv} drop {n_editions} editions -> "
              f"I={lfd.I:,} (-{fd.I - lfd.I:,}) J={lfd.J:,} (-{fd.J - lfd.J}) "
              f"N={lfd.N:,}  loglik={fr.loglik_final:.2f} {st['wall_s']:.1f}s", flush=True)
        n_new += 1

    # ---- meta ------------------------------------------------------------
    done = sorted(d.name[len("series_"):] for d in out_dir.glob("series_*")
                  if (d / "manifest.json").is_file())
    meta = dict(
        source_fit_dir=str(fit_dir), source_solver=src_solver, solver=solver,
        slice_slug=slug, stem=stem, identity_key=man.get("identity_key", "?"),
        unit="series",
        source_I=fd.I, source_J=fd.J, source_N=fd.N,
        registered_loglik=float(base_model.log_lik()),
        baseline_vs_registered_spearman=gap_spearman,
        baseline_vs_registered_maxdv=gap_maxdv,
        max_iter=args.max_iter, tol=args.tol,
        series_done=done,
        series_field={_short_series(k): int(v) for k, v in field.head(40).items()},
        series_editions={_short_series(k): int(n_ed.get(k, 0)) for k in field.head(40).index},
        updated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )
    tmp = meta_path.with_name(meta_path.name + ".tmp")
    tmp.write_text(json.dumps(meta, indent=2))
    tmp.replace(meta_path)

    return (f"done: {slug}/{out_dir.name}  (+{n_new} series refits, {len(done)} total)")


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
    g.add_argument("--series", nargs="+", default=["boston"],
                   help="series to leave out, one refit each. Accepts the short "
                        "label (boston) or the full key (boston_marathon).")
    g.add_argument("--top-k", type=int, default=None,
                   help="instead of --series, take the K highest-field series "
                        "(sum of edition field sizes) in the slice.")
    ap.add_argument("--redo", nargs="+", default=None, metavar="S",
                    help="re-run these series (short labels) even if on disk.")
    ap.add_argument("--solver", default=None, choices=list(_SOLVERS),
                    help="refit solver (default: same as the source fit).")
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--overwrite", action="store_true",
                    help="wipe all series dirs (incl. baseline) and start fresh.")
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
