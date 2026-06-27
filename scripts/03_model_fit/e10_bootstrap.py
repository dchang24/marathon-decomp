"""Bayesian (weighted) athlete bootstrap for uncertainty of the race-side and
global factors of an already-fitted model.

Quantifies how `v_j`, the aging curve (`theta_aging` / `gamma`) and the global
scalars (`sigma2`, `omega_d2`, `nu`) move as the athlete sample varies, via
strictly-positive Dirichlet weights over athletes folded into the fit as a
per-cell weight (`model.boot_w_cell`). See `marathon_decomp.resample` for why
this scheme (not the classic with-replacement athlete bootstrap) and how the
weight composes with the runner-reliability weight.

WHY WE DO NOT STORE u_i / d_i (the conceptual point). Both gauges the model
imposes every iteration are functions of `v` and the *fixed* race dates only,
never of `u`: `mean(v)=0` uses `v.mean()`, and the beta=0 APC gauge removes the
slope `c = cov(v, t_j)/var(t_j)` (model.apply_apc_gauge_beta0 reads only `v` and
`t_j`). So the saved `v_j`, `theta_aging`, `gamma` and scalars are already in the
final gauge at convergence and need nothing else to be reconstructed. Because the
athlete bootstrap keeps the race SET (and dates) fixed, every replicate shares
the same gauge reference, so the stored v_j / curve / scalars are directly
comparable across replicates. `u_i` / `d_i` would only be needed for fitted
values, athlete-ranking spread, or a *different* gauge defined on the u
distribution -- none of which is the goal here. (`--save-athlete` keeps them
anyway, for ranking-stability audits.)

To plot the aging curve from a replicate's `theta_aging` / `gamma`, use the
*point fit's* basis metadata (`fit.pkl` -> model_extra) with
`marathon_decomp.aging.aging_curve_from_payload` -- the basis is data-derived and
identical across replicates (athletes/races are fixed), so one basis serves all.

DISCOVERY. `--fit` may be a single fit dir, a slice dir, or the whole
`results/models` root; it is walked for fit dirs (those with a manifest.json).
`--pattern` selects which to bootstrap by the leaf dir name: a token with glob
chars (`* ? [`) is matched with fnmatch, otherwise as a substring. A fit is kept
if it matches ANY pattern. Examples: 'baseline_nu8*', 'baseline', 'AD_'.

APPENDABLE. Output lands in `<fit_dir>/bootstrap/`. Re-running ADDS replicates:
run N=5 for a dry run, then N=100 later and you have 105. Replicate `run_id`
continues from the existing max and each replicate's seed is `seed0 + run_id`, so
appends never repeat a draw and nothing is wasted. `--overwrite` wipes and
restarts. Run 0 (seed=-1) is the loaded point estimate, written once on first
creation for reference.

CRASH-SAFE. Each replicate is written the instant it finishes, through a temp
file + atomic os.replace -- so a kill (OOM, native abort, Ctrl-C) loses at most
the in-flight replicate and never truncates the live parquet. Just re-run the
same command to top up from where it stopped (run_id continues from the max on
disk). The optional `--save-athlete` table is streamed to per-replicate part
files under `athlete_parts/` (folded into athlete_factors.parquet at the end, or
on the next run after a kill) so the big I-row table is never rewritten per
replicate.

Outputs (under <fit_dir>/bootstrap/):
  runs.parquet           one row per replicate: run_id, seed, scheme, the
                         per-run solver + tol used, scalars (sigma2, omega_d2,
                         nu), loglik/rss, n_iter, converged, wall_s. solver/tol
                         are stored per row so appends with different settings
                         stay self-describing.
  race_factors.parquet   long: run_id, seed, race_idx, race_id, v.
  global_coeffs.parquet  long: run_id, seed, block, k, value (theta_aging/gamma);
                         only blocks active in the model are written.
  athlete_factors.parquet  (only --save-athlete) long: run_id, seed, u, d.
  meta.json              provenance + cumulative replicate count + run signature.

Run::

    # whole slice (default: anderson refits, warm-started). ALS is too slow.
    python scripts/03_model_fit/e10_bootstrap.py --fit results/models/Po10_M_14-25_mrc2 --N 100
    # pick models across all slices by name fragment
    python scripts/03_model_fit/e10_bootstrap.py --fit results/models --pattern baseline_nu8 --N 100
    # a single fit, dry run then top up
    python scripts/03_model_fit/e10_bootstrap.py --fit results/models/Po10_M_14-25_mrc2/baseline_nu8p00_best__cfc1baf4 --N 5
    python scripts/03_model_fit/e10_bootstrap.py --fit results/models/Po10_M_14-25_mrc2/baseline_nu8p00_best__cfc1baf4 --N 100
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pickle
import shutil
import time
from datetime import datetime, timezone
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
import pandas as pd

from marathon_decomp import (
    AndersonFitterConfig,
    FitterConfig,
    Model,
    ModelAnderson,
    bayesian_athlete_weights,
    boot_cell_weights,
    load_slice,
    registry,
)

# short solver name -> (model class, fitter-config class)
_SOLVERS: dict[str, tuple[type, type]] = {
    "als": (Model, FitterConfig),
    "anderson": (ModelAnderson, AndersonFitterConfig),
}

# ---- tunables (CLI flags override) ----------------------------------------
SEED0 = 12345          # replicate run_id n uses seed SEED0 + n
TOL = 1e-11            # tight: replicates must converge or warm-start biases low
MAX_ITER = 2000
CONCENTRATION = 1.0    # Dirichlet conc; 1.0 = standard Rubin Bayesian bootstrap
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Discovery
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
# Loading a point fit + warm-start config
# ---------------------------------------------------------------------------

def _load_point_fit(fit_dir: Path):
    """Read a saved fit, returning (model, payload, fd, identity_key). Guards
    that the slice reconstructed from the saved spec still matches the fit's
    shape -- a mismatch means the data export changed and the warm start is
    invalid."""
    with open(fit_dir / "fit.pkl", "rb") as f:
        payload = pickle.load(f)
    fd = load_slice(payload["spec"], payload["data_version"])
    model = registry.load_fit(fit_dir, fd)

    v = model.params.get("v")
    if v is not None and len(v) != fd.J:
        raise SystemExit(
            f"{fit_dir}: saved v has length {len(v)} but the slice now has "
            f"J={fd.J} races -- the data export changed since this fit. Re-fit "
            f"the point estimate before bootstrapping it.")

    man = json.loads((fit_dir / "manifest.json").read_text())
    return model, payload, fd, man.get("identity_key", "?")


def _warm_model_cfg(base_cfg, base_params: dict):
    """Copy the model config, seeding variance components from the point fit so
    the EB ridges start warm. Only seeds a component that is being estimated."""
    repl: dict = {"sigma2_init": float(base_params["sigma2"])}
    if base_cfg.omega_d2_fixed is None and not base_cfg.freeze_eb_prior:
        repl["omega_d2_init"] = float(base_params["omega_d2"])
    return dataclasses.replace(base_cfg, **repl)


def _make_fitter(solver: str, base_params: dict, *, init: str,
                 max_iter: int, tol: float):
    ModelCls, CfgCls = _SOLVERS[solver]
    kw = dict(max_outer_iter=max_iter, tol=tol, stop_criterion="loglik",
              record_trace=False, verbose=0)
    if init == "warmstart":
        kw.update(init="warmstart", warmstart=base_params)
    else:
        kw.update(init="mean")
    return ModelCls, CfgCls(**kw)


# ---------------------------------------------------------------------------
# Factor collection (race + global + scalars; athlete optional)
# ---------------------------------------------------------------------------

def _scalar_row(run_id: int, seed: int, scheme: str, params: dict, fr,
                wall_s: float, *, solver: str, tol: float) -> dict:
    return dict(
        run_id=run_id, seed=seed, scheme=scheme,
        solver=solver, tol=float(tol),     # per-run settings (may vary on append)
        sigma2=float(params["sigma2"]),
        omega_d2=float(params["omega_d2"]),
        nu=float(params["nu"]),
        loglik_final=(float(fr.loglik_final) if fr is not None else np.nan),
        rss_final=(float(fr.rss_final) if fr is not None else np.nan),
        n_iter=(int(fr.n_iter) if fr is not None else 0),
        converged=(bool(fr.converged) if fr is not None else True),
        wall_s=float(wall_s),
    )


def _collect_factors(run_id: int, seed: int, model, fd, *, save_athlete: bool):
    """Long rows for v_j (race), active global coeff vectors, and optionally
    u_i / d_i, for one replicate."""
    p = model.params
    cfg = model.config
    v = np.asarray(p["v"], dtype=np.float64)
    race = [dict(run_id=run_id, seed=seed, race_idx=j,
                 race_id=int(fd.race_ids[j]), v=float(v[j]))
            for j in range(fd.J)]

    glob: list[dict] = []
    if cfg.use_phi12:
        for k, val in enumerate(np.atleast_1d(np.asarray(p["theta_aging"], float))):
            glob.append(dict(run_id=run_id, seed=seed, block="theta_aging",
                             k=int(k), value=float(val)))
    if cfg.use_gamma:
        for k, val in enumerate(np.atleast_1d(np.asarray(p["gamma"], float))):
            glob.append(dict(run_id=run_id, seed=seed, block="gamma",
                             k=int(k), value=float(val)))

    ath: list[dict] = []
    if save_athlete:
        u = np.asarray(p["u"], dtype=np.float64)
        dd = np.asarray(p["d"], dtype=np.float64)
        for i in range(fd.I):
            ath.append(dict(run_id=run_id, seed=seed, athlete_idx=i,
                            athlete_id=int(fd.athlete_ids[i]),
                            u=float(u[i]), d=float(dd[i])))
    return race, glob, ath


# ---------------------------------------------------------------------------
# Crash-safe persistence helpers
#
# Each replicate is written the moment it finishes, so a kill (OOM, native
# abort, Ctrl-C) never loses completed work -- it can be topped up on the next
# run. Writes go through a temp file + os.replace so an interrupted write can
# never truncate/corrupt the live parquet (only the throwaway .tmp). The small
# tables (runs / race_factors / global_coeffs) are read-concat-rewritten per
# replicate (cheap: a few hundred rows each). The optional athlete table is
# huge (I rows per replicate), so it is streamed to per-replicate part files
# under athlete_parts/ and folded into athlete_factors.parquet only at the end
# (or, after a kill, on the next run's consolidation pass).
# ---------------------------------------------------------------------------

def _atomic_write_parquet(path: Path, df: pd.DataFrame, *, retries: int = 5) -> None:
    """Write `df` to `path` atomically (temp file + os.replace). Retries the
    replace a few times: OneDrive can briefly hold a handle on the target and
    raise PermissionError on Windows."""
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


def _write_athlete_part(out_dir: Path, run_id: int, ath_rows: list[dict]) -> None:
    """Stream one replicate's u_i / d_i to its own part file (no rewrite of the
    growing consolidated table)."""
    if not ath_rows:
        return
    parts = out_dir / "athlete_parts"
    parts.mkdir(exist_ok=True)
    _atomic_write_parquet(parts / f"{run_id}.parquet", pd.DataFrame(ath_rows))


def _consolidate_athlete(out_dir: Path) -> None:
    """Fold athlete_parts/*.parquet into athlete_factors.parquet. Idempotent and
    append-safe: rows already in the consolidated file whose run_id is not among
    the parts (e.g. a legacy table from an older run) are preserved."""
    parts = out_dir / "athlete_parts"
    files = sorted(parts.glob("*.parquet"), key=lambda p: int(p.stem)) if parts.is_dir() else []
    final = out_dir / "athlete_factors.parquet"
    frames = [pd.read_parquet(f) for f in files]
    part_ids = {int(f.stem) for f in files}
    if final.is_file():
        old = pd.read_parquet(final)
        old = old[~old["run_id"].isin(part_ids)]
        if len(old):
            frames.insert(0, old)
    if not frames:
        return
    _atomic_write_parquet(final, pd.concat(frames, ignore_index=True))


# ---------------------------------------------------------------------------
# Bootstrap one fit
# ---------------------------------------------------------------------------

def run_fit(fit_dir: Path, args) -> str:
    out_dir = fit_dir / "bootstrap"
    runs_path = out_dir / "runs.parquet"
    meta_path = out_dir / "meta.json"

    if args.overwrite and out_dir.is_dir():
        shutil.rmtree(out_dir)

    base_model, payload, fd, identity_key = _load_point_fit(fit_dir)
    base_cfg = base_model.config
    base_params = {k: (np.asarray(v).copy() if np.ndim(v) else v)
                   for k, v in base_model.params.items()}

    # ---- compatibility / append bookkeeping --------------------------------
    fresh = not runs_path.is_file()
    start_run_id = 0
    if not fresh:
        prev = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        if prev.get("identity_key") not in (None, identity_key):
            return (f"REFUSE: {fit_dir.name}/bootstrap has identity "
                    f"{prev.get('identity_key')} != fit {identity_key}; "
                    f"--overwrite to redo.")
        if prev.get("J") not in (None, fd.J):
            return (f"REFUSE: {fit_dir.name}/bootstrap was built on J="
                    f"{prev.get('J')} != current J={fd.J}; --overwrite to redo.")
        if prev.get("concentration") not in (None, float(args.concentration)):
            return (f"REFUSE: {fit_dir.name}/bootstrap used concentration="
                    f"{prev.get('concentration')} != {args.concentration}; a "
                    f"different Dirichlet draw is a different distribution -- "
                    f"--overwrite to redo or match --concentration.")
        start_run_id = int(pd.read_parquet(runs_path)["run_id"].max())

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {fit_dir.parent.name}/{fit_dir.name} ===", flush=True)
    print(f"    {payload['name']} v{payload['fitter_version']}  "
          f"I={fd.I:,} J={fd.J:,} N={fd.N:,}  base loglik={base_model.log_lik():.4f}",
          flush=True)
    print(f"    +{args.N} replicates (existing through run_id={start_run_id}), "
          f"solver={args.solver} init={args.init} tol={args.tol:g} "
          f"conc={args.concentration}", flush=True)

    ModelCls = _SOLVERS[args.solver][0]
    warm_cfg = _warm_model_cfg(base_cfg, base_params)
    base_loglik = float(base_model.log_lik())

    def persist_meta(total: int) -> None:
        """Write meta.json atomically. Refreshed after every replicate so the
        compatibility checks (identity/J/concentration) and the cumulative count
        are never stale, even if a later replicate is killed."""
        meta = dict(
            scheme="bayesian_athlete",
            identity_key=identity_key,
            source_name=payload["name"],
            source_fitter_version=payload["fitter_version"],
            source_fit_dir=str(fit_dir),
            solver=args.solver, init=args.init,
            concentration=float(args.concentration),
            seed0=int(args.seed0), tol=args.tol, max_iter=args.max_iter,
            I=fd.I, J=fd.J, N_obs=fd.N,
            n_replicates=total,             # cumulative (run 0 is the point fit)
            base_loglik=base_loglik,
            save_athlete=bool(args.save_athlete),
            updated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )
        tmp = meta_path.with_name(meta_path.name + ".tmp")
        tmp.write_text(json.dumps(meta, indent=2))
        os.replace(tmp, meta_path)

    def persist_replicate(run_row: dict, race_rows: list[dict],
                          glob_rows: list[dict], ath_rows: list[dict]) -> None:
        """Durably write one replicate the instant it finishes (atomic)."""
        _append_parquet(runs_path, [run_row])
        _append_parquet(out_dir / "race_factors.parquet", race_rows)
        _append_parquet(out_dir / "global_coeffs.parquet", glob_rows)
        _write_athlete_part(out_dir, run_row["run_id"], ath_rows)
        persist_meta(run_row["run_id"])

    # Run 0: the loaded point estimate (reference), written only on first creation.
    if fresh:
        base_row = _scalar_row(0, -1, "base", base_model.params,
                               base_model.fit_result, 0.0,
                               solver="point", tol=np.nan)
        r0, g0, a0 = _collect_factors(0, -1, base_model, fd,
                                      save_athlete=args.save_athlete)
        persist_replicate(base_row, r0, g0, a0)

    t0 = time.perf_counter()
    n_conv = 0
    for k in range(1, args.N + 1):
        run_id = start_run_id + k
        seed = args.seed0 + run_id
        rng = np.random.default_rng(seed)
        g = bayesian_athlete_weights(fd, rng, concentration=args.concentration)

        cfg = warm_cfg if args.init == "warmstart" else base_cfg
        _, fcfg = _make_fitter(args.solver, base_params, init=args.init,
                               max_iter=args.max_iter, tol=args.tol)
        model = ModelCls(fd, cfg, fcfg)
        model.boot_w_cell = boot_cell_weights(fd, g)

        ti = time.perf_counter()
        fr = model.fit()
        dt = time.perf_counter() - ti

        run_row = _scalar_row(run_id, seed, "bayesian_athlete",
                              model.params, fr, dt,
                              solver=args.solver, tol=args.tol)
        rb, gb, ab = _collect_factors(run_id, seed, model, fd,
                                      save_athlete=args.save_athlete)
        persist_replicate(run_row, rb, gb, ab)
        n_conv += int(fr.converged)

        conv = "OK " if fr.converged else "MAX"
        print(f"    run {run_id:4d} seed={seed} iters={fr.n_iter:4d} {conv} "
              f"loglik={fr.loglik_final:.4f} sigma2={model.params['sigma2']:.5f} "
              f"{dt:.2f}s", flush=True)

    if args.save_athlete:
        _consolidate_athlete(out_dir)

    total = int(pd.read_parquet(runs_path)["run_id"].max())
    persist_meta(total)
    return (f"done: {fit_dir.name}/bootstrap  (+{args.N}, {n_conv}/{args.N} "
            f"converged; {total} replicates total, {time.perf_counter()-t0:.1f}s)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True,
                    help="a fit dir, a slice dir, or the results/models root "
                         "(walked for fit dirs). A non-existent trailing "
                         "component is matched as a prefix against siblings, so "
                         "'results/models/Po10' selects every Po10_* slice.")
    ap.add_argument("--pattern", nargs="+", default=None, metavar="PAT",
                    help="select fits by leaf dir name: glob chars (* ? [) use "
                         "fnmatch, else substring. Kept if it matches ANY. "
                         "Ignored when --fit is a single fit dir.")
    ap.add_argument("--N", type=int, default=10,
                    help="number of replicates to ADD this run (appends).")
    ap.add_argument("--solver", choices=list(_SOLVERS), default="anderson",
                    help="refit solver (default anderson; ALS is too slow). The "
                         "fixed point is solver-invariant.")
    ap.add_argument("--init", choices=["warmstart", "mean"], default="warmstart",
                    help="'warmstart' from the point fit (default, fast); 'mean' "
                         "cold-starts each replicate (convergence audit).")
    ap.add_argument("--seed0", type=int, default=SEED0,
                    help="replicate run_id n uses seed seed0 + run_id.")
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--concentration", type=float, default=CONCENTRATION,
                    help="Dirichlet concentration (1.0 = Rubin bootstrap). Must "
                         "match across appends into one output dir.")
    ap.add_argument("--save-athlete", action="store_true",
                    help="also write u_i / d_i per replicate (large; ranking "
                         "stability only -- NOT needed for v/curve/global spread).")
    ap.add_argument("--overwrite", action="store_true",
                    help="wipe each fit's bootstrap/ dir and start fresh.")
    args = ap.parse_args()

    roots = resolve_roots(args.fit)
    fits: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        for fd in discover_fits(root, args.pattern):
            if fd not in seen:
                seen.add(fd)
                fits.append(fd)
    if not fits:
        hint = f" matching {args.pattern}" if args.pattern else ""
        where = roots[0] if len(roots) == 1 else f"{len(roots)} roots"
        raise SystemExit(f"no fits found under {where}{hint}")
    print(f"discovered {len(fits)} fit(s) under {len(roots)} root(s)")

    statuses = [run_fit(fd, args) for fd in fits]
    print("\n---- summary ----")
    for st in statuses:
        print(f"  {st}")


if __name__ == "__main__":
    main()
