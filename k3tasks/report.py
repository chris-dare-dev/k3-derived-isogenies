"""The report of a run directory: the scheduler's tables, then whether each computation was reproduced.

Each task stands under one word, its verdict or the scheduler's status (see `expectation.py`).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

from taskbatch import report as scheduler_report
from taskbatch.report import load_manifest
from taskbatch.store import ResultStore

from . import expectation
from .tasks import TaskPayload


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:g}"
    return str(value).replace("|", "\\|")


def _table(header: list[str], rows: list[list[Any]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(_cell(item) for item in row) + " |" for row in rows]
    return [*lines, ""]


def _differences(record: Mapping[str, Any]) -> str:
    parts = expectation.differences(record)
    work_dir = record.get("detail", {}).get("work_dir")
    if work_dir:
        parts.append(f"work directory {work_dir}")
    return "; ".join(parts)


def _peak(record: Mapping[str, Any] | None) -> float | None:
    """The task's peak memory: of all its processes together where its scope recorded it, else of its
    largest process."""
    if record is None:
        return None
    peak = record.get("detail", {}).get("scope_peak_mb") or record["peak_rss_mb"]
    return None if peak is None else float(peak)


def reproduction_lines(run_dir: Path) -> list[str]:
    latest = ResultStore(run_dir / "results.jsonl").latest()
    manifest = load_manifest(run_dir)
    ids = [task["id"] for task in manifest] or sorted(latest)
    rows = []
    counts: Counter[str] = Counter()
    for task_id in ids:
        record = latest.get(task_id)
        listed = record["task"] if record else next((t for t in manifest if t["id"] == task_id), None)
        task = cast("TaskPayload | None", listed)
        word = expectation.verdict(record)
        counts[word] += 1
        compared = ", ".join(expectation.compared_with(record))
        rows.append(
            [
                f"`{task_id}`",
                task["group"] if task else "",
                word,
                compared,
                record["wall_seconds"] if record else None,
                _peak(record),
            ]
        )
    summary = ", ".join(f"{counts[word]} {word}" for word in expectation.ORDER if counts[word])
    lines = ["## Reproduction", "", f"{summary}.", ""]
    lines += _table(["task", "group", "verdict", "compared with", "wall (s)", "peak (MB)"], rows)
    problems = [
        (task_id, latest[task_id])
        for task_id in ids
        if task_id in latest and expectation.verdict(latest[task_id]) == expectation.DIFFERS
    ]
    if problems:
        lines += ["## What differs", ""]
        lines += _table(["task", "difference"], [[f"`{t}`", _differences(r)] for t, r in problems])
    caveats = [
        (t, latest[t]["detail"]["caveat"])
        for t in ids
        if t in latest and latest[t].get("detail", {}).get("caveat")
    ]
    if caveats:
        lines += ["## Caveats", ""]
        lines += _table(["task", "caveat"], [[f"`{t}`", c] for t, c in caveats])
    return lines


def render(run_dir: Path, *, full: bool = False) -> str:
    return (
        scheduler_report.render(run_dir, full=full).rstrip()
        + "\n\n"
        + "\n".join(reproduction_lines(run_dir)).rstrip()
        + "\n"
    )
