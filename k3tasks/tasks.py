"""From a batch description to tasks, and from a task to its result (taskbatch/PROTOCOL.md).

A `[[tasks]]` table of a batch file selects computations of the registry:

    [[tasks]]
    kind = "computation"
    ids = ["n27/ledger", "n02/enumerate-*"]   # ids or shell patterns
    notes = [2, 3]                            # and/or notes
    groups = ["lattices"]                     # and/or claim families
    tags = ["heavy"]                          # and/or tags
    exclude = ["n03/enumerate-q49"]           # ids or patterns left out
    exclude_tags = ["heavy"]
    procs = 4                                 # at most this many processes per task
    keep_work = false                         # keep the work directory even when reproduced
    limits = { time_seconds = 600, memory_mb = 2000 }   # read by the scheduler only

A computation selected by two tables is taken once, with the first table's settings.
"""

from __future__ import annotations

from typing import Any

from . import registry, workspace

TASK_KEYS = {"kind", "ids", "notes", "groups", "tags", "exclude", "exclude_tags", "procs", "keep_work"}


def _list(table: dict[str, Any], key: str, kind: type) -> list[Any] | None:
    value = table.get(key)
    if value is None:
        return None
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        value = [value]
    if not isinstance(value, list) or not all(isinstance(v, kind) and not isinstance(v, bool) for v in value):
        raise ValueError(f"{key} must be a list of {kind.__name__}")
    return value


def expand_batch(description: dict[str, Any]) -> list[dict[str, Any]]:
    """The tasks of a batch description: one per selected computation."""
    tasks: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, table in enumerate(description.get("tasks", [])):
        if table.get("kind") != "computation":
            raise ValueError(f"tasks[{index}]: unknown kind {table.get('kind')!r}; the only kind is 'computation'")
        unknown = set(table) - TASK_KEYS
        if unknown:
            raise ValueError(f"tasks[{index}]: unknown key {sorted(unknown)[0]!r}")
        selectors = {key: _list(table, key, str) for key in ("ids", "groups", "tags", "exclude", "exclude_tags")}
        notes = _list(table, "notes", int)
        if not any(selectors[key] is not None for key in ("ids", "groups", "tags")) and notes is None:
            raise ValueError(f"tasks[{index}]: select computations with ids, notes, groups or tags")
        procs_cap = table.get("procs")
        if procs_cap is not None and (not isinstance(procs_cap, int) or isinstance(procs_cap, bool) or procs_cap < 1):
            raise ValueError(f"tasks[{index}]: procs must be a positive whole number")
        keep = table.get("keep_work", False)
        if not isinstance(keep, bool):
            raise ValueError(f"tasks[{index}]: keep_work must be true or false")
        chosen = registry.select(notes=notes, **selectors)
        if not chosen:
            raise ValueError(f"tasks[{index}] selects no computation")
        for computation in chosen:
            if computation.id in seen:
                continue
            seen.add(computation.id)
            procs = computation.procs if procs_cap is None else min(computation.procs, procs_cap)
            tasks.append(
                {
                    "id": computation.id,
                    "kind": "computation",
                    "set": index,
                    "group": computation.group,
                    "cost": round(computation.expected_seconds(procs), 1),
                    "procs": procs,
                    "expected_memory_mb": computation.expected_memory_mb(procs),
                    "keep_work": keep,
                    "command": computation.command(procs),
                }
            )
    return tasks


def run_task(task: dict[str, Any], places: workspace.Places | None = None) -> dict[str, Any]:
    """Runs one computation and compares what it produced with the stored files."""
    computation = registry.get(task["id"])
    places = places or workspace.Places.from_environment()
    procs = int(task.get("procs", computation.procs))
    execution = workspace.execute(computation, places, procs)
    if execution.exit_code != 0:
        message = workspace.tail(execution.stderr if not computation.merge_stderr else execution.stdout)
        raise workspace.WorkspaceError(
            f"{computation.script} exited with status {execution.exit_code}; work directory kept at "
            f"{execution.work_dir}\n{message}"
        )
    result = workspace.evaluate(computation, execution, places)
    return workspace.finish(execution, result, keep=bool(task.get("keep_work")))
