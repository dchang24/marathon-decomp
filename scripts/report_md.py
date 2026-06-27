"""Shared helper: render QC/summary numbers as a single readable Markdown file.

The ``q*`` / ``q99`` summary scripts collect their paper-quoted numbers as
``(section, label, value, source)`` rows. This module turns those rows into one
Markdown document -- an H1 title, optional italic subtitle lines, then one H2 +
table per section (columns: Quantity | Value | Source). The point is that the
paper can cite a single ``.md`` per analysis that is readable in any text editor
and cleanly git-trackable, instead of pointing footnotes at raw CSV/parquet.

Shared lib (no script prefix); importable as ``from report_md import ...`` because
the summary scripts put ``scripts/`` on ``sys.path``.

Markdown body is ASCII (table pipes, backticks, hyphens), so printing it to a
cp1252 Windows console is safe.
"""
from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Iterable, Sequence


Row = tuple[str, str, str, str]   # (section, label, value, source)


def _cell(s: object) -> str:
    """Make a value safe for a Markdown table cell (escape pipes, trim)."""
    return str(s).replace("|", "\\|").strip()


def render_markdown(
    title: str,
    rows: Iterable[Row],
    subtitle: str | Sequence[str] | None = None,
    notes: Sequence[str] | None = None,
) -> str:
    """Render ``(section, label, value, source)`` rows as a Markdown report.

    ``subtitle`` is one or more italic lines under the title (model / decision
    context). ``notes`` is an optional trailing bullet list. The ``source`` of
    each row is rendered as inline code so file/column references stay verbatim.
    """
    out: list[str] = [f"# {title}", ""]
    if subtitle:
        for line in ([subtitle] if isinstance(subtitle, str) else list(subtitle)):
            out.append(f"*{line}*")
            out.append("")
    out += [f"_Generated {date.today().isoformat()}._", ""]

    last_section = None
    for section, label, value, source in rows:
        if section != last_section:
            out += ["", f"## {section}", "",
                    "| Quantity | Value | Source |", "|---|---|---|"]
            last_section = section
        out.append(f"| {_cell(label)} | {_cell(value)} | `{_cell(source)}` |")

    if notes:
        out += ["", "## Notes", ""]
        out += [f"- {n}" for n in notes]
    out.append("")
    return "\n".join(out)


def write_markdown(path: Path, text: str) -> Path:
    """Write ``text`` to ``path`` (utf-8, parents created); returns the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8")
    return path
