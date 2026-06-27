from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # src/marathon_decomp/ -> src/ -> project root

try:
    from ._private_paths import (
        MASTER_DB_PATH,
        SPLINK_MATCHES_PATH,
        RULE_MATCHES_PATH,
    )
except ImportError:
    MASTER_DB_PATH = SPLINK_MATCHES_PATH = RULE_MATCHES_PATH = None

DATA_DIR      = ROOT / "data" / "race_results"        # real export (no version subdir)
FAKE_DATA_DIR = ROOT / "data" / "race_results_fake"   # synthetic, shareable copy
MISC_DIR      = ROOT / "data" / "misc"                # misc reference/QC tables
REFERENCE_DIR = ROOT / "data" / "reference"           # curated reference tables (race_catalogue, ...)

# covariate data-prep domains (weather + course elevation; see scripts/01_data_prep/)
WEATHER_DIR           = ROOT / "data" / "weather"                 # derived weather features + query points
WEATHER_OPENMETEO_DIR = WEATHER_DIR / "openmeteo"                # raw per-race hourly open-meteo CSVs
COURSE_PROFILE_DIR = ROOT / "data" / "course_profile"            # raw/ gpx, processed/ csv, derived elevation
COURSE_RAW_DIR     = COURSE_PROFILE_DIR / "raw"                   # source gpx (NOT shared)
COURSE_PROCESSED_DIR = COURSE_PROFILE_DIR / "processed"          # slimmed lat/lon/ele csv (shareable)
RESULTS_DIR = ROOT / "results"
SCRATCH_DIR = ROOT / "scratch"
DOCS_DIR    = ROOT / "docs"
CACHE_DIR   = ROOT / "cache" / "fitdata"   # on-disk FitData cache (see data.py)
PAPER_FIG_DIR = ROOT / "paper" / "v1_final" / "figures"


def display_path(p) -> str:
    """Render a path for human-readable provenance lines (e.g. the "CSV -> ..."
    footers in QC summaries / manifests).

    Returns the path *relative to the repo ROOT* in POSIX form, so summary
    outputs never embed an absolute machine path (username, OneDrive, etc.).
    Falls back to the bare filename for paths outside the repo. Use this instead
    of ``str(path)`` anywhere a path is written into a file a reader might see.
    """
    p = Path(p).resolve()
    try:
        return p.relative_to(ROOT).as_posix()
    except ValueError:
        return p.name