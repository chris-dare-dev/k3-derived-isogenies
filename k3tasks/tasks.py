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

`expand_batch` gives the tasks of a batch description and `run_task` the result of one task: the two
commands of the runner.  `plan_batch` reads a batch file as the scheduler does and gives its tasks with
their limits, for everything else that needs the tasks of a batch.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, TypedDict

from taskbatch.config import BatchConfig, load_config
from taskbatch.model import Limits, Task
from taskbatch.settings import load_settings

from . import expectation, registry, workspace

TASK_KEYS = {
    "kind",
    "ids",
    "notes",
    "groups",
    "tags",
    "exclude",
    "exclude_tags",
    "procs",
    "keep_work",
    "control",
}
SYNTHETIC_KEYS = {"name", "cost", "sleep_seconds", "allocate_mb", "fail"}
CONTROL_SUFFIX = "!control"
CONTROL_GROUP = "negative control"


class Compared(TypedDict):
    """A stored file that a rerun is compared with, and how (compare.py)."""

    file: str
    mode: str


class _Placed(TypedDict):
    """What every task has: its place in the batch, and what it is expected to take."""

    id: str
    set: int  # index of its [[tasks]] table: the scheduler takes the limits of that table
    group: str
    cost: float  # expected seconds
    procs: int
    expected_memory_mb: int


class _Control(TypedDict, total=False):
    """What a negative control has besides: the parts of what is predicted that were replaced."""

    control: bool
    wrong_on_purpose: dict[str, Any]


class ComputationTask(_Placed, _Control):
    """One computation of the registry, or its negative control."""

    kind: Literal["computation"]
    keep_work: bool
    command: str
    note: int
    predicted: dict[str, Any] | None  # keys and values as in expectation.py; None when nothing is compared
    compared: list[Compared]
    caveat: str | None


class SyntheticTask(_Placed):
    """A task that computes nothing: it sleeps, holds memory or fails."""

    kind: Literal["synthetic"]
    sleep_seconds: float
    allocate_mb: float
    fail: bool


# A task as the scheduler hands it back to the runner and stores it in a record: the payload of
# taskbatch/PROTOCOL.md.  The committed records hold these keys, so none is added, dropped or renamed.
TaskPayload = ComputationTask | SyntheticTask


def _list(table: dict[str, Any], key: str, kind: type) -> list[Any] | None:
    value = table.get(key)
    if value is None:
        return None
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        value = [value]
    if not isinstance(value, list) or not all(isinstance(v, kind) and not isinstance(v, bool) for v in value):
        raise ValueError(f"{key} must be a list of {kind.__name__}")
    return value


def expand_batch(description: dict[str, Any]) -> list[TaskPayload]:
    """The tasks of a batch description: one per selected computation.

    Every way to the tasks of a batch passes here (the runner's `expand`, and `plan_batch` for a batch
    file), so each of them refuses a registry that is not consistent.
    """
    problems = registry.validate()
    if problems:
        raise ValueError("the registry is inconsistent: " + "; ".join(problems))
    tasks: list[TaskPayload] = []
    seen: set[str] = set()
    for index, table in enumerate(description.get("tasks", [])):
        if table.get("kind") == "synthetic":
            for synthetic in _synthetic(table, index):
                if synthetic["id"] in seen:
                    raise ValueError(f"tasks[{index}]: two synthetic tasks are named {synthetic['id']}")
                seen.add(synthetic["id"])
                tasks.append(synthetic)
            continue
        if table.get("kind") != "computation":
            kind = table.get("kind")
            raise ValueError(
                f"tasks[{index}]: unknown kind {kind!r}; the kinds are 'computation' and 'synthetic'"
            )
        unknown = set(table) - TASK_KEYS
        if unknown:
            raise ValueError(f"tasks[{index}]: unknown key {sorted(unknown)[0]!r}")
        selectors = {
            key: _list(table, key, str) for key in ("ids", "groups", "tags", "exclude", "exclude_tags")
        }
        notes = _list(table, "notes", int)
        if not any(selectors[key] is not None for key in ("ids", "groups", "tags")) and notes is None:
            raise ValueError(f"tasks[{index}]: select computations with ids, notes, groups or tags")
        procs_cap = table.get("procs")
        if procs_cap is not None and (
            not isinstance(procs_cap, int) or isinstance(procs_cap, bool) or procs_cap < 1
        ):
            raise ValueError(f"tasks[{index}]: procs must be a positive whole number")
        keep = table.get("keep_work", False)
        if not isinstance(keep, bool):
            raise ValueError(f"tasks[{index}]: keep_work must be true or false")
        control = table.get("control")
        if control is not None and (
            not isinstance(control, dict)
            or set(control) != {"predicted"}
            or not isinstance(control["predicted"], dict)
        ):
            raise ValueError(
                f"tasks[{index}]: control must be a table with the one key predicted, itself a table"
            )
        chosen = registry.select(notes=notes, **selectors)
        if not chosen:
            raise ValueError(f"tasks[{index}] selects no computation")
        for computation in chosen:
            task_id = computation.id + (CONTROL_SUFFIX if control else "")
            if task_id in seen:
                continue
            seen.add(task_id)
            procs = computation.procs if procs_cap is None else min(computation.procs, procs_cap)
            predicted = computation.predicted()
            task: ComputationTask = {
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
                "predicted": predicted or None,
                "compared": [{"file": check.stored, "mode": check.mode} for check in computation.checks],
                "caveat": computation.caveat or None,
            }
            if control:
                unknown = set(control["predicted"]) - set(predicted)
                if unknown:
                    gives = ", ".join(predicted) or "nothing"
                    raise ValueError(
                        f"tasks[{index}]: the control predicts {sorted(unknown)[0]!r}, which "
                        f"{computation.id} does not give; it gives {gives}"
                    )
                task["group"] = CONTROL_GROUP
                task["control"] = True
                task["predicted"] = predicted | control["predicted"]
                task["wrong_on_purpose"] = control["predicted"]
            tasks.append(task)
    return tasks


@dataclass(frozen=True)
class Plan:
    """A batch file and its tasks."""

    config: BatchConfig
    tasks: list[TaskPayload]

    def limits(self, task: TaskPayload) -> Limits:
        """The limits the scheduler gives a task of the batch."""
        return self.config.limits_for(Task.from_payload(dict(task)))


def plan_batch(path: Path) -> Plan:
    """The tasks of a batch file, read and expanded as the scheduler does it."""
    config = load_config(path, description_keys=load_settings().description_keys)
    return Plan(config, expand_batch(config.description))


def _synthetic(table: dict[str, Any], index: int) -> list[SyntheticTask]:
    """Tasks that compute nothing: they sleep, hold memory or fail, to show the scheduler's limits."""
    unknown = set(table) - {"kind", "items"}
    if unknown:
        raise ValueError(f"tasks[{index}]: unknown key {sorted(unknown)[0]!r}")
    items = table.get("items")
    if not isinstance(items, list) or not all(isinstance(item, dict) and "name" in item for item in items):
        raise ValueError(f"tasks[{index}]: items must be a list of tables with a name")
    tasks: list[SyntheticTask] = []
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


def _run_synthetic(task: SyntheticTask) -> dict[str, Any]:
    if task["fail"]:
        raise workspace.WorkspaceError("synthetic task asked to fail")
    held = bytearray(b"x") * (int(task["allocate_mb"]) << 20)
    time.sleep(float(task["sleep_seconds"]))
    return {
        "status": "ok",
        "value": len(held) >> 20,
        "predicted": None,
        "agrees": None,
        "verdict": expectation.RAN,
        "method": "synthetic",
    }


def run_task(task: TaskPayload, places: workspace.Places | None = None) -> dict[str, Any]:
    """Runs one computation and compares what it produced with the stored files."""
    if task["kind"] == "synthetic":
        return _run_synthetic(task)
    computation = registry.get(task["id"].removesuffix(CONTROL_SUFFIX))
    return workspace.reproduce(
        computation,
        places or workspace.Places.from_environment(),
        task["procs"],
        control=task.get("wrong_on_purpose") if task.get("control") else None,
        keep=task["keep_work"],
    )
