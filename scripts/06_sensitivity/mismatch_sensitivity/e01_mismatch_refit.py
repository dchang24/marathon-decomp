"""e01 -- identity-mismatch sensitivity: perturb, refit, record (S5.7).

Treats a registered production fit's slice as ground truth and asks how the fit
moves when additional identity error is injected at a controlled per-athlete
rate ``p``, for three operations:

    break  recall loss   -- split a fraction p of runners into two sub-runners
    join   precision loss -- fuse a fraction p of runners with a compatible partner
    both   join then break, each at rate p

For each (op, p, replicate) it builds a perturbed slice (``perturb.perturb``),
re-fits the SAME model config as the source fit (cold ``mean`` init, tight tol --
the athlete set changes so a warm start from the point fit is invalid), and
records the race-side (v) and global (theta_aging/gamma) factors per replicate.
``run_id=0`` (op="none", p=0) is a cold refit of the UNPERTURBED slice -- the
self-contained comparison baseline.

APPENDABLE (Monte-Carlo style). Re-running ADDS replicates: it counts the
existing replicates per (op, p) cell on disk and appends ``--N`` more, continuing
``run_id`` from that cell's current max. The replicate seed is ``seed0 + run_id``
and depends only on the replicate index, NOT on op/p -- so (a) re-runs never
repeat a draw, and (b) replicates at the same index share a seed across cells
(common random numbers), which de-noises the v(p) trend. ``--overwrite`` wipes
and restarts. A new ``--ops`` / ``--p`` value is simply a new cell that starts at
run_id 1; you never have to redo what's already there.

CRASH-SAFE. Each replicate is written the instant it finishes via a temp file +
atomic os.replace, so a kill (OOM, native abort, Ctrl-C) loses at most the
in-flight replicate and never truncates the live parquet. Just re-run the same
command to top up.

DISCOVERY (mirrors scripts/03_model_fit/e10_bootstrap.py). ``--fit`` may be a
single fit dir, a slice dir, or the whole ``results/models`` root; it is walked
for fit dirs (those with a manifest.json). A non-existent trailing component is
matched as a prefix against siblings, so ``results/models/Po10`` selects every
``Po10_*`` slice. ``--pattern`` filters by the leaf dir name: a token with glob
chars (``* ? [``) uses fnmatch, else substring; a fit is kept if it matches ANY.
Examples: 'full', 'AD_*', 'agingS4gv'.

Outputs under results/sensitivity/mismatch_sensitivity/{slice}/{model}/:
  runs.parquet          one row per replicate: op, p, run_id, seed, record_frac,
                        n_break_moves, n_join_moves, I/J/N, solver, tol, scalars,
                        convergence, wall_s.
  race_factors.parquet  long: op, p, run_id, race_id, v.
  global_coeffs.parquet long: op, p, run_id, block, k, value (theta_aging/gamma).
  meta.json             provenance + run signature (for append-compat checks).

Run::

    # one fit by partial path + name pattern, 50 reps per cell
    python scripts/06_sensitivity/mismatch_sensitivity/e01_mismatch_refit.py --fit results/models/Po10_M_14-25_mrc2 --pattern full --ops break join both --p 0.001 0.002 0.005 0.01 0.02 0.05 0.1 --N 50
    # later: top up to 100 reps per cell (just re-run with --N 50)
    # every full fit across all Po10 slices:
    python scripts/06_sensitivity/mismatch_sensitivity/e01_mismatch_refit.py \
        --fit results/models/Po10 --pattern full --N 50
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))          # this dir (perturb)
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

import perturb as P  # noqa: E402

OUT_ROOT = RESULTS_DIR / "sensitivity" / "mismatch_sensitivity"

# short solver name -> (model class, fitter-config class)
_SOLVERS: dict[str, tuple[type, type]] = {
    "als": (Model, FitterConfig),
    "anderson": (ModelAnderson, AndersonFitterConfig),
}
# manifest "solver" (== model NAME) -> short name
_NAME_TO_SHORT = {
    Model.NAME: "als",
    ModelAnderson.NAME: "anderson",
}

SEED0 = 90210
TOL = 1e-10
MAX_ITER = 2000


# ---------------------------------------------------------------------------
# Discovery (mirrors scripts/03_model_fit/e10_bootstrap.py)
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
    """Resolve `--fit` to a list of existing search roots.

    If the path exists, it is the single root. Otherwise the trailing component
    is treated as a prefix (or glob, if it has glob chars) matched against the
    sibling directories of its parent -- so `results/models/Po10` matches every
    `Po10_*` slice dir, and `results/models/ALL` matches every `ALL_*`."""
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
    """Every fit dir under `root` (or `root` itself if it is one) matching
    `patterns` on the leaf dir name. A fit dir is one holding a manifest.json."""
    if (root / "manifest.json").is_file():
        return [root]   # an explicit single fit dir is never pattern-filtered
    out = []
    for man in sorted(root.rglob("manifest.json")):
        d = man.parent
        if _name_matches(d.name, patterns):
            out.append(d)
    return out


# ---------------------------------------------------------------------------
# Loading the point fit
# ---------------------------------------------------------------------------

def _load_point_fit(fit_dir: Path):
    """Reload a registered point fit + its (cached) slice. Guards that the saved
    factors still match the slice shape -- a mismatch means the export changed."""
    import pickle
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
    return model, fd, short, man


# ---------------------------------------------------------------------------
# Row / factor collection
# ---------------------------------------------------------------------------

def _scalar_row(params: dict, fr) -> dict:
    return dict(
        sigma2=float(params["sigma2"]),
        omega_d2=float(params["omega_d2"]),
        nu=float(params["nu"]),
        loglik_final=(float(fr.loglik_final) if fr is not None else np.nan),
        rss_final=(float(fr.rss_final) if fr is not None else np.nan),
        n_iter=(int(fr.n_iter) if fr is not None else 0),
        converged=(bool(fr.converged) if fr is not None else True),
    )


def _collect_factors(op: str, p: float, run_id: int, model, fd):
    """Long rows for v_j (race) and the active global coeff vectors, one
    replicate. Per-athlete factors (u/d) are NOT collected: break/join change the
    athlete SET, so u_i/d_i are not comparable across perturbations."""
    params = model.params
    cfg = model.config
    v = np.asarray(params["v"], dtype=np.float64)
    race = [dict(op=op, p=p, run_id=run_id, race_id=int(fd.race_ids[j]), v=float(v[j]))
            for j in range(fd.J)]
    glob: list[dict] = []
    if getattr(cfg, "use_phi12", False):
        for k, val in enumerate(np.atleast_1d(np.asarray(params["theta_aging"], float))):
            glob.append(dict(op=op, p=p, run_id=run_id, block="theta_aging",
                             k=int(k), value=float(val)))
    if getattr(cfg, "use_gamma", False):
        for k, val in enumerate(np.atleast_1d(np.asarray(params["gamma"], float))):
            glob.append(dict(op=op, p=p, run_id=run_id, block="gamma",
                             k=int(k), value=float(val)))
    return race, glob


# ---------------------------------------------------------------------------
# Crash-safe persistence (atomic temp-file + os.replace, append in place)
# ---------------------------------------------------------------------------

def _atomic_write_parquet(path: Path, df: pd.DataFrame, *, retries: int = 5) -> None:
    """Write `df` to `path` atomically. Retries the replace a few times: OneDrive
    can briefly hold a handle on the target and raise PermissionError on Windows."""
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False)
    for i in range(retries):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if i == retries - 1:
                raise
            time.sleep(0.2 * (i + 1))


def _append_parquet(path: Path, new_rows: list[dict]) -> None:
    """Atomically concat new rows onto an existing parquet (if any) and rewrite."""
    if not new_rows:
        return
    df_new = pd.DataFrame(new_rows)
    if path.is_file():
        df_new = pd.concat([pd.read_parquet(path), df_new], ignore_index=True)
    _atomic_write_parquet(path, df_new)


# ---------------------------------------------------------------------------
# Append-compatibility signature
# ---------------------------------------------------------------------------

def _signature(identity_key: str, fd, args) -> dict:
    """Fields that define the perturbation *distribution*; they must match across
    appends into one output dir. ops/p are intentionally NOT here -- adding a new
    op or p value is just a new cell, never a conflict. solver/tol live per-row."""
    return dict(identity_key=identity_key, J=int(fd.J),
                split_mode=args.split_mode, tau=float(args.tau),
                seed0=int(args.seed0))


# ---------------------------------------------------------------------------
# Run one fit
# ---------------------------------------------------------------------------

def run_fit(fit_dir: Path, args) -> str:
    base_model, fd, src_solver, man = _load_point_fit(fit_dir)
    solver = args.solver or src_solver
    identity_key = man.get("identity_key", "?")
    key = fit_dir.name.rpartition("__")[2] or "nokey"
    stem = fit_dir.name.rpartition("__")[0] or fit_dir.name
    slug = fit_dir.parent.name
    out_dir = OUT_ROOT / slug / f"{stem}__{key}"
    runs_path = out_dir / "runs.parquet"
    race_path = out_dir / "race_factors.parquet"
    glob_path = out_dir / "global_coeffs.parquet"
    meta_path = out_dir / "meta.json"

    if args.overwrite and out_dir.is_dir():
        shutil.rmtree(out_dir)

    sig = _signature(identity_key, fd, args)

    # ---- append bookkeeping -----------------------------------------------
    fresh = not runs_path.is_file()
    existing = None
    if not fresh:
        prev = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        for fld, want in sig.items():
            got = prev.get(fld)
            if got is not None and got != want:
                return (f"REFUSE: {slug}/{out_dir.name} has {fld}={got} != {want}; "
                        f"a different perturbation distribution -- --overwrite to redo.")
        existing = pd.read_parquet(runs_path, columns=["op", "p", "run_id"])

    out_dir.mkdir(parents=True, exist_ok=True)
    ModelCls, CfgCls = _SOLVERS[solver]
    cfg = base_model.config
    fcfg = CfgCls(max_outer_iter=args.max_iter, tol=args.tol, stop_criterion="loglik",
                  init="mean", record_trace=False, verbose=0)

    print(f"\n=== mismatch sensitivity: {slug}/{stem} (solver={solver}) ===", flush=True)
    print(f"    registered I={fd.I:,} J={fd.J:,} N={fd.N:,}  "
          f"loglik={base_model.log_lik():.2f}", flush=True)
    print(f"    +{args.N} reps/cell  ops={args.ops} p={args.p} "
          f"split_mode={args.split_mode} seed0={args.seed0}", flush=True)

    def cell_max_run_id(op: str, p: float) -> int:
        if existing is None or len(existing) == 0:
            return 0
        m = existing[(existing.op == op) & (np.isclose(existing.p, p))]
        return int(m.run_id.max()) if len(m) else 0

    def persist_meta() -> None:
        n_done = int((pd.read_parquet(runs_path, columns=["op"])["op"] != "none").sum())
        meta = dict(
            **sig,                                   # carries identity_key, J, split_mode, tau, seed0
            source_fit_dir=str(fit_dir), source_solver=src_solver, solver=solver,
            slice_slug=slug, stem=stem, I=fd.I, N=fd.N,
            registered_loglik=float(base_model.log_lik()),
            max_iter=args.max_iter, tol=args.tol,
            n_replicates=n_done,
            updated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )
        tmp = meta_path.with_name(meta_path.name + ".tmp")
        tmp.write_text(json.dumps(meta, indent=2))
        os.replace(tmp, meta_path)

    def persist(run_row: dict, race_rows: list[dict], glob_rows: list[dict]) -> None:
        _append_parquet(runs_path, [run_row])
        _append_parquet(race_path, race_rows)
        _append_parquet(glob_path, glob_rows)

    # ---- run 0: cold refit of the UNPERTURBED slice (baseline) -------------
    # Written once on fresh creation. The registered fit was built by warm-start
    # continuation and may land at a marginally higher loglik; comparing perturbed
    # cold refits against it would mis-attribute that cold-vs-warm gap to the
    # perturbation, so we compare against this matched cold baseline and log the gap.
    if fresh:
        base_cold = ModelCls(fd, cfg, fcfg)
        base_fr = base_cold.fit()
        vw = np.asarray(base_model.params["v"]); vc = np.asarray(base_cold.params["v"])
        gap_spearman = float(spearmanr(vw, vc).statistic)
        gap_maxdv = float(np.abs(vw - vc).max())
        print(f"    cold baseline loglik={base_fr.loglik_final:.2f} "
              f"(vs registered d={base_model.log_lik() - base_fr.loglik_final:.2f}; "
              f"v spearman={gap_spearman:.5f} max|dv|={gap_maxdv:.4f})", flush=True)
        base_row = dict(op="none", p=0.0, run_id=0, seed=-1, record_frac=0.0,
                        n_break_moves=0, n_join_moves=0,
                        n_break_splits=0, n_join_merges=0, n_athletes_base=fd.I,
                        I=fd.I, J=fd.J, N=fd.N,
                        solver="point", tol=float("nan"), wall_s=0.0,
                        **_scalar_row(base_cold.params, base_fr))
        r0, g0 = _collect_factors("none", 0.0, 0, base_cold, fd)
        persist(base_row, r0, g0)
        persist_meta()

    # ---- the (op, p) x N sweep --------------------------------------------
    ps = [float(x) for x in args.p]
    t0 = time.perf_counter()
    n_new = 0
    for op in args.ops:
        for p in ps:
            base_n = cell_max_run_id(op, p)
            for k in range(1, args.N + 1):
                run_id = base_n + k
                seed = args.seed0 + run_id           # depends only on the index (CRN)
                rng = np.random.default_rng(seed)
                pert = P.perturb(fd, op, p, rng, split_mode=args.split_mode, tau=args.tau)
                pfd = pert.fd
                model = ModelCls(pfd, cfg, fcfg)
                ti = time.perf_counter()
                fr = model.fit()
                dt = time.perf_counter() - ti

                run_row = dict(op=op, p=p, run_id=run_id, seed=seed,
                               record_frac=pert.record_frac,
                               n_break_moves=pert.n_break_moves,
                               n_join_moves=pert.n_join_moves,
                               n_break_splits=pert.n_break_splits,
                               n_join_merges=pert.n_join_merges,
                               n_athletes_base=pert.n_athletes_base,
                               I=pfd.I, J=pfd.J, N=pfd.N,
                               solver=solver, tol=float(args.tol), wall_s=dt,
                               **_scalar_row(model.params, fr))
                rb, gb = _collect_factors(op, p, run_id, model, pfd)
                persist(run_row, rb, gb)
                persist_meta()
                n_new += 1

                # realised per-athlete rate: each split touches 1 athlete, each
                # merge touches 2 (absorber + absorbed), over the baseline count.
                p_real = ((pert.n_break_splits + 2 * pert.n_join_merges)
                          / max(pert.n_athletes_base, 1))
                conv = "OK " if fr.converged else "MAX"
                print(f"    {op:5s} p={p:<6g} run={run_id:3d} {conv} J={pfd.J:4d} "
                      f"p_real={p_real:.4f} rec_frac={pert.record_frac:.3f} "
                      f"loglik={fr.loglik_final:.2f} {dt:.1f}s", flush=True)

    persist_meta()
    total = int((pd.read_parquet(runs_path, columns=["op"])["op"] != "none").sum())
    return (f"done: {slug}/{out_dir.name}  (+{n_new} refits, {total} total, "
            f"{time.perf_counter() - t0:.1f}s)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True,
                    help="a fit dir, a slice dir, or the results/models root "
                         "(walked for fit dirs). A non-existent trailing component "
                         "is matched as a prefix against siblings, so "
                         "'results/models/Po10' selects every Po10_* slice.")
    ap.add_argument("--pattern", nargs="+", default=None, metavar="PAT",
                    help="select fits by leaf dir name: glob chars (* ? [) use "
                         "fnmatch, else substring. Kept if it matches ANY. "
                         "Ignored when --fit is a single fit dir.")
    ap.add_argument("--ops", nargs="+", default=["break", "join", "both"],
                    choices=["break", "join", "both"])
    ap.add_argument("--p", nargs="+",
                    default=["0.001", "0.002", "0.005", "0.01", "0.02", "0.05", "0.10"],
                    help="per-athlete perturbation rates.")
    ap.add_argument("--N", type=int, default=50,
                    help="replicates to ADD per (op, p) cell this run (appends).")
    ap.add_argument("--split-mode", default="datecut",
                    choices=["datecut", "iid", "singleton"])
    ap.add_argument("--tau", type=float, default=P.TAU_LOGTIME,
                    help="join finish-time gate (|delta median log t|); default ~log(1.2).")
    ap.add_argument("--solver", default=None, choices=list(_SOLVERS),
                    help="refit solver (default: same as the source fit).")
    ap.add_argument("--seed0", type=int, default=SEED0,
                    help="replicate run_id n uses seed seed0 + n (must match on append).")
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--overwrite", action="store_true",
                    help="wipe the output dir and start fresh.")
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
