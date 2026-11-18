"""The report of a run directory: the scheduler's tables, then whether each computation was reproduced.

Each task ends with one verdict:

    reproduced   everything it was compared with is identical to the stored files
    differs      something differs; the table says what, and the work directory was kept
    checked      no stored file to compare with, but the lines its output must contain are there
    ran          the command ran to the end and nothing was compared (its output was not kept)
    timeout, memory, error, crashed   the scheduler's status: no verdict
    not run      no record yet
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from taskbatch import report as scheduler_report
from taskbatch.store import ResultStore

ORDER = ["reproduced", "differs", "checked", "ran", "timeout", "memory", "error", "crashed", "unsupported", "not run"]


def verdict(record: dict[str, Any] | None) -> str:
    if record is None:
        return "not run"
    if record["status"] in ("ok", "mismatch"):
        return str(record.get("detail", {}).get("verdict") or ("differs" if record["status"] == "mismatch" else "ran"))
    return str(record["status"])


def _differences(record: dict[str, Any]) -> str:
    detail = record.get("detail", {})
    parts = []
    for label, note in (detail.get("differences") or {}).items():
        if (record.get("value") or {}).get(label) == "identical":
            continue
        shown = ", ".join(f"{k} {v}" for k, v in note.items() if k != "mode")
        parts.append(f"{label} ({note.get('mode')}): {shown}")
    for key, value in (record.get("value") or {}).items():
        expected = (record.get("predicted") or {}).get(key)
        if key.startswith(("expect ", "lines ")) and value != expected:
            parts.append(f"{key}: {value}, expected {expected}")
    if detail.get("work_dir"):
        parts.append(f"work directory {detail['work_dir']}")
    return "; ".join(parts)


def _peak(record: dict[str, Any] | None) -> float | None:
    """The task's peak memory: of all its processes together where its scope recorded it, else of its
    largest process."""
    if record is None:
        return None
    return record.get("detail", {}).get("scope_peak_mb") or record["peak_rss_mb"]


def reproduction_lines(run_dir: Path) -> list[str]:
    latest = ResultStore(run_dir / "results.jsonl").latest()
    manifest = scheduler_report.load_manifest(run_dir)
    ids = [task["id"] for task in manifest] or sorted(latest)
    rows = []
    counts: Counter[str] = Counter()
    for task_id in ids:
        record = latest.get(task_id)
        task = record["task"] if record else next((t for t in manifest if t["id"] == task_id), {})
        word = verdict(record)
        counts[word] += 1
        compared = ", ".join(k for k in (record or {}).get("value") or {} if not k.startswith(("expect ", "lines ")))
        rows.append(
            [
                f"`{task_id}`",
                task.get("group", ""),
                word,
                compared,
                record["wall_seconds"] if record else None,
                _peak(record),
            ]
        )
    summary = ", ".join(f"{counts[word]} {word}" for word in ORDER if counts[word])
    lines = ["## Reproduction", "", f"{summary}.", ""]
    lines += scheduler_report._table(["task", "group", "verdict", "compared with", "wall (s)", "peak (MB)"], rows)
    problems = [(task_id, latest[task_id]) for task_id in ids if task_id in latest and verdict(latest[task_id]) == "differs"]
    if problems:
        lines += ["## What differs", ""]
        lines += scheduler_report._table(["task", "difference"], [[f"`{t}`", _differences(r)] for t, r in problems])
    caveats = [(t, latest[t]["detail"]["caveat"]) for t in ids if t in latest and latest[t].get("detail", {}).get("caveat")]
    if caveats:
        lines += ["## Caveats", ""]
        lines += scheduler_report._table(["task", "caveat"], [[f"`{t}`", c] for t, c in caveats])
    return lines


def render(run_dir: Path, *, full: bool = False) -> str:
    return scheduler_report.render(run_dir, full=full).rstrip() + "\n\n" + "\n".join(reproduction_lines(run_dir)).rstrip() + "\n"
