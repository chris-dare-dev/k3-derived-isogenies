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

A table with `control = { predicted = { ... } }` makes negative controls: the same computations under ids
ending in "!control", with the given parts of what is predicted replaced.  The replacement is wrong on
purpose, so each control must end `mismatch`; if one ends `ok`, the comparison is broken.

    [[tasks]]
    kind = "computation"
    ids = ["n27/ledger"]
    control = { predicted = { "sixlines_ledger.log" = "differs" } }

A second kind computes nothing.  It shows what the scheduler's limits look like in the records:

    [[tasks]]
    kind = "synthetic"
    items = [{ name = "sleeps", sleep_seconds = 300 }, { name = "allocates", allocate_mb = 1000 }]

Controls and synthetic tasks belong to no claim.
"""

from __future__ import annotations

import time
from typing import Any

from . import registry, workspace

TASK_KEYS = {"kind", "ids", "notes", "groups", "tags", "exclude", "exclude_tags", "procs", "keep_work", "control"}
SYNTHETIC_KEYS = {"name", "cost", "sleep_seconds", "allocate_mb", "fail"}
CONTROL_SUFFIX = "!control"
CONTROL_GROUP = "negative control"


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
        if table.get("kind") == "synthetic":
            for task in _synthetic(table, index):
                if task["id"] in seen:
                    raise ValueError(f"tasks[{index}]: two synthetic tasks are named {task['id']}")
                seen.add(task["id"])
                tasks.append(task)
            continue
        if table.get("kind") != "computation":
            raise ValueError(
                f"tasks[{index}]: unknown kind {table.get('kind')!r}; the kinds are 'computation' and 'synthetic'"
            )
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
        control = table.get("control")
        if control is not None and (
            not isinstance(control, dict) or set(control) != {"predicted"} or not isinstance(control["predicted"], dict)
        ):
            raise ValueError(f"tasks[{index}]: control must be a table with the one key predicted, itself a table")
        chosen = registry.select(notes=notes, **selectors)
        if not chosen:
            raise ValueError(f"tasks[{index}] selects no computation")
        for computation in chosen:
            task_id = computation.id + (CONTROL_SUFFIX if control else "")
            if task_id in seen:
                continue
            seen.add(task_id)
            procs = computation.procs if procs_cap is None else min(computation.procs, procs_cap)
            task = {
                "id": task_id,
                "kind": "computation",
                "set": index,
                "group": computation.group,
                "cost": round(computation.expected_seconds(procs), 1),
                "procs": procs,
                "expected_memory_mb": computation.expected_memory_mb(procs),
                "keep_work": keep,
                "command": computation.command(procs),
                "note": computation.note,
                "predicted": computation.predicted() or None,
                "compared": [{"file": check.stored, "mode": check.mode} for check in computation.checks],
                "caveat": computation.caveat or None,
            }
            if control:
                unknown = set(control["predicted"]) - set(computation.predicted())
                if unknown:
                    raise ValueError(
                        f"tasks[{index}]: the control predicts {sorted(unknown)[0]!r}, which {computation.id} does "
                        f"not give; it gives {', '.join(computation.predicted()) or 'nothing'}"
                    )
                task |= {
                    "group": CONTROL_GROUP,
                    "control": True,
                    "predicted": computation.predicted() | control["predicted"],
                    "wrong_on_purpose": control["predicted"],
                }
            tasks.append(task)
    return tasks


def _synthetic(table: dict[str, Any], index: int) -> list[dict[str, Any]]:
    """Tasks that compute nothing: they sleep, hold memory or fail, to show the scheduler's limits."""
    unknown = set(table) - {"kind", "items"}
    if unknown:
        raise ValueError(f"tasks[{index}]: unknown key {sorted(unknown)[0]!r}")
    items = table.get("items")
    if not isinstance(items, list) or not all(isinstance(item, dict) and "name" in item for item in items):
        raise ValueError(f"tasks[{index}]: items must be a list of tables with a name")
    tasks = []
    for item in items:
        unknown = set(item) - SYNTHETIC_KEYS
        if unknown:
            raise ValueError(f"tasks[{index}]: unknown key {sorted(unknown)[0]!r} in an item")
        tasks.append(
            {
                "id": f"synthetic/{item['name']}",
                "kind": "synthetic",
                "set": index,
                "group": "synthetic",
                "cost": item.get("cost", 0),
                "procs": 1,
                "expected_memory_mb": int(item.get("allocate_mb", 0)),
                "sleep_seconds": item.get("sleep_seconds", 0),
                "allocate_mb": item.get("allocate_mb", 0),
                "fail": bool(item.get("fail", False)),
            }
        )
    return tasks


def _run_synthetic(task: dict[str, Any]) -> dict[str, Any]:
    if task["fail"]:
        raise workspace.WorkspaceError("synthetic task asked to fail")
    held = bytearray(b"x") * (int(task["allocate_mb"]) << 20)
    time.sleep(float(task["sleep_seconds"]))
    return {
        "status": "ok",
        "value": len(held) >> 20,
        "predicted": None,
        "agrees": None,
        "verdict": "ran",
        "method": "synthetic",
    }


def run_task(task: dict[str, Any], places: workspace.Places | None = None) -> dict[str, Any]:
    """Runs one computation and compares what it produced with the stored files."""
    if task.get("kind") == "synthetic":
        return _run_synthetic(task)
    control = task.get("wrong_on_purpose") if task.get("control") else None
    computation = registry.get(str(task["id"]).removesuffix(CONTROL_SUFFIX))
    places = places or workspace.Places.from_environment()
    procs = int(task.get("procs", computation.procs))
    execution = workspace.execute(computation, places, procs)
    if execution.exit_code != 0:
        message = workspace.tail(execution.stderr if not computation.merge_stderr else execution.stdout)
        raise workspace.WorkspaceError(
            f"{computation.script} exited with status {execution.exit_code}; work directory kept at "
            f"{workspace.shown(execution.work_dir)}\n{message}"
        )
    result = workspace.evaluate(computation, execution, places, control)
    return workspace.finish(execution, result, keep=bool(task.get("keep_work")), discard=control is not None)
