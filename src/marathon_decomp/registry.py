"""Run registry for model-fit sweeps.

A *fit* is identified by (slice, model config, resample tag). The solver
(als/anderson) and the numeric fitter knobs are **provenance**, not identity —
als and anderson are the same model and converge to the same fixed point, so by
default they land in the same folder. To keep solver / init / seed variants
side-by-side (convergence studies), give them distinct `resample_tag`s.

Layout — you own every path segment; the registry only appends an 8-char
identity hash so reusing a stem can't silently clobber a different config:

    results/fits/{slice_slug}/{your_stem}__{key}/
        manifest.json     # canonical: full identity + metrics, flat columns
        fit.pkl           # BaseModel.save() payload (params by default)

`manifest.json` is the source of truth; `results/registry.parquet` is a
*derived* index you can delete and `rebuild_registry()` from the manifests.
New config knobs therefore need no migration — they show up as new columns
on the next rebuild, with NaN for older rows.

Typical sweep cell::

    from marathon_decomp import registry, load_slice, Model

    stem = f"{registry.terms_tag(cfg)}-poly{cfg.degree}"      # you name it
    parent = RESULTS_DIR / "fits" / registry.slice_slug(spec)
    fit_dir = registry.fit_path(parent, stem, spec, cfg)

    pf = registry.preflight(parent, stem, spec, cfg)
    if pf.already_run:
        continue                          # exact params already on disk
    pf.warn()                             # prints any same-name/diff-param clashes

    m = Model(load_slice(spec), cfg, fcfg); m.fit()
    registry.register_fit(m, fit_dir)
"""
from __future__ import annotations

import hashlib
import json
import math
import warnings
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseModel, SaveSpec
from .config import RESULTS_DIR
from .data import SliceSpec
from .models import Model, ModelAnderson

DEFAULT_FITS_ROOT = RESULTS_DIR / "fits"
DEFAULT_REGISTRY = RESULTS_DIR / "registry.parquet"

# solver NAME -> class, for load_fit()
_MODEL_BY_NAME: dict[str, type[BaseModel]] = {
    Model.NAME: Model,
    ModelAnderson.NAME: ModelAnderson,
}


# ---------------------------------------------------------------------------
# Identity: what makes two fits "the same"
# ---------------------------------------------------------------------------

def _identity_dict(spec: SliceSpec, cfg: Any, resample_tag: str) -> dict[str, Any]:
    """The fields that define fit identity: the data slice, the model config,
    and the resample tag. Fitter numerics are deliberately excluded — they are
    provenance (recorded in the manifest), not identity.
    """
    return {
        "slice": asdict(spec),
        "config": asdict(cfg) if is_dataclass(cfg) else dict(cfg),
        "resample_tag": resample_tag,
    }


def identity_key(
    spec: SliceSpec, cfg: Any, resample_tag: str = "base", *, n: int = 8,
) -> str:
    """Stable short hash of fit identity (slice + model config + resample tag)."""
    payload = json.dumps(_identity_dict(spec, cfg, resample_tag),
                         sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:n]


# ---------------------------------------------------------------------------
# Optional name helpers (you may ignore these and name folders by hand)
# ---------------------------------------------------------------------------

def slice_slug(spec: SliceSpec) -> str:
    """Readable slice stem: e.g. ``Po10_M_14-25_mrc5``.

    cohort_gender_YY-YY_mrcN, with gender B for both sexes. Append your own
    markers for special slices (e.g. ``+ "_LOCO-KEN"``).
    """
    gender = {"M": "M", "W": "W", "ALL": "B"}.get(spec.sex, spec.sex)
    ylo = spec.date_lo.year % 100 if spec.date_lo else 0
    yhi = spec.date_hi.year % 100 if spec.date_hi else 99
    return f"{spec.cohort}_{gender}_{ylo:02d}-{yhi:02d}_mrc{spec.min_race_count}"


def terms_tag(cfg: Any) -> str:
    """2-char term toggle tag: A=aging(+gamma), D=drift, x=off."""
    a = "A" if getattr(cfg, "use_phi12", False) else "x"
    d = "D" if getattr(cfg, "use_d", False) else "x"
    return f"{a}{d}"


def nu_tag(nu: float) -> str:
    """Noise tag for a fit stem: ``L2`` for Gaussian (nu=inf), else ``nu5p23``
    (2-dp, decimal point -> 'p'). Matches the Brent-refined 2-dp grid."""
    return "L2" if not math.isfinite(float(nu)) else f"nu{float(nu):.2f}".replace(".", "p")


def model_stem(cfg: Any, solver: str) -> str:
    """Readable fit-folder stem ``{terms}_{nu}_{solver}``.

    `terms` is the 2-char toggle tag (e.g. ``AD``, ``Ax``), collapsed to the
    word ``baseline`` when every optional term is off (the rank-1 u+v model).
    `nu` is :func:`nu_tag`; `solver` is the short solver name (``als`` /
    ``anderson``). Examples: ``baseline_L2_als``, ``baseline_nu5p23_anderson``,
    ``AD_L2_als``.
    """
    terms = terms_tag(cfg)
    if terms == "xx":
        terms = "baseline"
    nu = float(getattr(cfg, "nu", float("inf")))
    return f"{terms}_{nu_tag(nu)}_{solver}"


def fit_path(
    parent: Path | str, stem: str,
    spec: SliceSpec, cfg: Any, resample_tag: str = "base",
) -> Path:
    """Compose ``parent/{stem}__{key}``. The key is the identity hash."""
    return Path(parent) / f"{stem}__{identity_key(spec, cfg, resample_tag)}"


# ---------------------------------------------------------------------------
# Preflight: catch reruns and same-name/different-param clashes
# ---------------------------------------------------------------------------

@dataclass
class PreflightReport:
    key: str
    stem: str
    exact_matches: list[Path]                   # same identity already on disk
    name_clashes: list[tuple[Path, dict[str, dict[str, Any]]]]  # (path, field->{existing,current})

    @property
    def already_run(self) -> bool:
        return bool(self.exact_matches)

    def summary(self) -> str:
        lines: list[str] = []
        if self.exact_matches:
            lines.append(
                f"[already run] identity {self.key} exists at: "
                + ", ".join(str(p) for p in self.exact_matches)
            )
        for path, diff in self.name_clashes:
            fields = ", ".join(
                f"{k}: {v['existing']!r}->{v['current']!r}" for k, v in diff.items()
            )
            lines.append(
                f"[name clash] {path.name} shares stem '{self.stem}' but differs in: {fields}"
            )
        return "\n".join(lines) if lines else f"[clear] no clashes for {self.stem}__{self.key}"

    def warn(self) -> None:
        """Emit a warning for any clash (no-op if clear)."""
        if self.exact_matches or self.name_clashes:
            warnings.warn(self.summary(), stacklevel=2)


def preflight(
    parent: Path | str, stem: str,
    spec: SliceSpec, cfg: Any, resample_tag: str = "base",
) -> PreflightReport:
    """Check `parent` for fits clashing with this (stem, identity).

    - exact_matches: any ``*__{key}`` dir (same params, even under a different
      stem) — i.e. this exact fit has already been run.
    - name_clashes: ``{stem}__*`` dirs with a *different* key, annotated with
      which identity fields differ (read from their manifests).
    """
    parent = Path(parent)
    key = identity_key(spec, cfg, resample_tag)
    current = _flatten_identity(_identity_dict(spec, cfg, resample_tag))

    exact_matches: list[Path] = []
    name_clashes: list[tuple[Path, dict]] = []
    if parent.is_dir():
        for d in sorted(parent.iterdir()):
            if not d.is_dir() or "__" not in d.name:
                continue
            d_stem, _, d_key = d.name.rpartition("__")
            if d_key == key:
                exact_matches.append(d)
            elif d_stem == stem:
                name_clashes.append((d, _diff_identity(d, current)))
    return PreflightReport(key=key, stem=stem,
                           exact_matches=exact_matches, name_clashes=name_clashes)


def _diff_identity(fit_dir: Path, current: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Which identity fields differ between an on-disk fit and `current`."""
    man = fit_dir / "manifest.json"
    if not man.is_file():
        return {}
    try:
        existing = json.loads(man.read_text()).get("_identity_flat", {})
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, dict[str, Any]] = {}
    for k in set(existing) | set(current):
        ev, cv = existing.get(k), current.get(k)
        if ev != cv:
            out[k] = {"existing": ev, "current": cv}
    return out


# ---------------------------------------------------------------------------
# Recording a fit
# ---------------------------------------------------------------------------

def register_fit(
    model: BaseModel,
    fit_dir: Path | str,
    *,
    save: SaveSpec | None = None,
    resample_tag: str = "base",
    warn_on_overwrite: bool = True,
    **extra: Any,
) -> Path:
    """Write ``fit.pkl`` + ``manifest.json`` into `fit_dir`.

    The manifest denormalizes every SliceSpec / model-config / fitter field
    plus headline metrics into flat columns, so the parquet index can query
    any knob regardless of how the folder was named. `extra` kwargs become
    additional manifest columns.
    """
    fit_dir = Path(fit_dir)
    spec = model.data.spec
    cfg = model.config
    ident = _identity_dict(spec, cfg, resample_tag)
    key = identity_key(spec, cfg, resample_tag)

    man_path = fit_dir / "manifest.json"
    if warn_on_overwrite and man_path.is_file():
        diff = _diff_identity(fit_dir, _flatten_identity(ident))
        if diff:
            fields = ", ".join(f"{k}: {v['existing']!r}->{v['current']!r}"
                               for k, v in diff.items())
            warnings.warn(
                f"overwriting {fit_dir} which holds a DIFFERENT config "
                f"(differs in: {fields})", stacklevel=2,
            )

    fit_dir.mkdir(parents=True, exist_ok=True)
    save = save or SaveSpec.default()
    model.save(man_path.with_name("fit.pkl"), what=save)

    manifest = _build_manifest(model, fit_dir, key, resample_tag, ident, extra)
    man_path.write_text(json.dumps(manifest, indent=2, default=_json_default))
    return fit_dir


def _build_manifest(
    model: BaseModel, fit_dir: Path, key: str, resample_tag: str,
    ident: dict[str, Any], extra: dict[str, Any],
) -> dict[str, Any]:
    spec, cfg, fitter = model.data.spec, model.config, getattr(model, "fitter", None)
    fr = model.fit_result

    row: dict[str, Any] = {
        "identity_key": key,
        "resample_tag": resample_tag,
        "solver": model.NAME,
        "fitter_version": model.VERSION,
        "data_version": model.data.data_version,
        "slice_key": model.data.cache_key,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    try:
        row["fit_dir"] = str(fit_dir.relative_to(RESULTS_DIR)).replace("\\", "/")
    except ValueError:
        row["fit_dir"] = str(fit_dir)

    # denormalized knobs (flat columns) — names don't collide across the three
    _merge_no_clobber(row, asdict(spec))
    _merge_no_clobber(row, asdict(cfg))
    if fitter is not None:
        _merge_no_clobber(row, asdict(fitter), skip={"warmstart"})

    # headline metrics
    row.update(
        I=model.data.I, J=model.data.J, N=model.data.N,
        n_params_naive=model.n_params_naive(),
        loglik=model.log_lik(),
        aic=model.aic(), bic=model.bic(),
        aic_naive=model.aic(effective=False), bic_naive=model.bic(effective=False),
    )
    if fr is not None:
        row.update(converged=fr.converged, n_iter=fr.n_iter,
                   rss_final=fr.rss_final, loglik_final=fr.loglik_final)

    row.update(extra)
    # keep the raw identity (flattened) for fast, exact diffing in preflight
    row["_identity_flat"] = _flatten_identity(ident)
    return row


def _merge_no_clobber(dst: dict, src: dict, skip: set[str] | None = None) -> None:
    skip = skip or set()
    for k, v in src.items():
        if k in skip:
            continue
        if k in dst and dst[k] != v:
            warnings.warn(f"manifest column name collision on {k!r}; keeping first",
                          stacklevel=3)
            continue
        dst[k] = _jsonable(v)


def _flatten_identity(ident: dict[str, Any]) -> dict[str, Any]:
    """Flatten the nested identity dict into ``slice.<f>`` / ``config.<f>`` keys."""
    out: dict[str, Any] = {"resample_tag": ident["resample_tag"]}
    for grp in ("slice", "config"):
        for k, v in ident[grp].items():
            out[f"{grp}.{k}"] = _jsonable(v)
    return out


# ---------------------------------------------------------------------------
# Index build + query
# ---------------------------------------------------------------------------

def rebuild_registry(
    root: Path | str = DEFAULT_FITS_ROOT,
    out: Path | str = DEFAULT_REGISTRY,
) -> pd.DataFrame:
    """Glob ``root/**/manifest.json`` into one DataFrame and write it to `out`."""
    root, out = Path(root), Path(out)
    rows: list[dict] = []
    for man in root.rglob("manifest.json"):
        try:
            d = json.loads(man.read_text())
        except (OSError, json.JSONDecodeError):
            warnings.warn(f"skipping unreadable manifest {man}", stacklevel=2)
            continue
        d.pop("_identity_flat", None)   # internal; not an index column
        rows.append(d)
    df = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not df.empty:
        df.to_parquet(out, index=False)
    return df


def query(
    root: Path | str = DEFAULT_FITS_ROOT,
    out: Path | str = DEFAULT_REGISTRY,
    *,
    rebuild: bool = False,
    **filters: Any,
) -> pd.DataFrame:
    """Load the index (rebuilding if missing/stale-by-request) and apply
    equality filters. Scalar filters match by ==; list/tuple match by isin.
    """
    out = Path(out)
    if rebuild or not out.is_file():
        df = rebuild_registry(root, out)
    else:
        df = pd.read_parquet(out)
    for col, val in filters.items():
        if col not in df.columns:
            warnings.warn(f"filter column {col!r} not in registry", stacklevel=2)
            return df.iloc[0:0]
        if isinstance(val, (list, tuple, set)):
            df = df[df[col].isin(list(val))]
        else:
            df = df[df[col] == val]
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Loading a saved fit back into a model
# ---------------------------------------------------------------------------

def load_fit(fit_dir_or_row: Path | str | pd.Series | dict, data: Any) -> BaseModel:
    """Reload a saved fit as a live model bound to `data` (a FitData).

    Accepts a fit directory, a path to its ``fit.pkl``, or a registry row
    (Series/dict) carrying ``fit_dir`` + ``solver``.
    """
    if isinstance(fit_dir_or_row, (pd.Series, dict)):
        row = fit_dir_or_row
        pkl = (RESULTS_DIR / row["fit_dir"] / "fit.pkl")
        cls = _MODEL_BY_NAME[row["solver"]]
        return cls.load(pkl, data)

    p = Path(fit_dir_or_row)
    if p.is_dir():
        p = p / "fit.pkl"
    # solver name lives in the sibling manifest
    man = p.with_name("manifest.json")
    name = json.loads(man.read_text())["solver"] if man.is_file() else Model.NAME
    return _MODEL_BY_NAME[name].load(p, data)


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

def _jsonable(v: Any) -> Any:
    if isinstance(v, (date, datetime)):
        return v.isoformat()
    if isinstance(v, tuple):
        return list(v)
    return v


def _json_default(v: Any) -> Any:
    if isinstance(v, (date, datetime)):
        return v.isoformat()
    return str(v)
