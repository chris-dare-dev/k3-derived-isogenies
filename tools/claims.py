#!/usr/bin/env python3
"""claims.py -- the claim tables, written from the claim registry, the batch files and the records.

Usage:  python3 tools/claims.py tables              write the generated parts of CLAIMS.md and GUIDE.md,
                                                     and claims/RECORDS.md
        python3 tools/claims.py check               fail if they are not what `tables` writes, if the claim
                                                     registry and the batch files do not fit together, or if
                                                     a stored line named as evidence is not what it says
        python3 tools/claims.py adopt RUN_DIR ...   copy the records of finished runs (results/<batch>)
                                                     into records/<batch>
        python3 tools/claims.py words               print what a committed record must not name (this
                                                     machine's name and directories), one per line; `make
                                                     report` hands them to the report as a warning

What is read:  claims/claims.toml (the claims and the tasks that bear on them), batches/*.toml (which tasks
exist), records/<batch>/results.jsonl (what was computed here), the import statements of experiments/*.py
(what the scripts share), and, when the data archive is unpacked, the stored lines named as evidence.
A task without a record reads "not re-run", whatever the stored results say.  A record of a task that
has changed since (another command, another predicted value) reads "stale" and supports nothing.

`adopt` is tools/adopt.py, which says what it copies and what it refuses.
"""

from __future__ import annotations

import ast
import fnmatch
import os
import re
import sys
import tomllib
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, TypedDict, TypeGuard, cast

REPOSITORY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPOSITORY))

import adopt  # noqa: E402

from k3tasks import expectation  # noqa: E402
from k3tasks.tasks import NOT_COMPUTED, ComputationTask, SyntheticTask, plan_batch  # noqa: E402
from taskbatch.record import answers  # noqa: E402
from taskbatch.store import ResultStore  # noqa: E402

REGISTRY = REPOSITORY / "claims" / "claims.toml"
RECORDS = REPOSITORY / "records"
DETAIL = REPOSITORY / "claims" / "RECORDS.md"
SCRIPTS = REPOSITORY / "experiments"
NOT_RERUN = "not re-run"
STALE = "stale"
STATUSES = ["ok", "mismatch", "unsupported", "timeout", "memory", "error", "crashed", STALE, NOT_RERUN]
STORED_NAME = re.compile(r"""["']((?:sixlines|ellk3)_[A-Za-z0-9_]+\.(?:jsonl|json|log))["']""")


class RegistryError(ValueError):
    """The registry, the batch files and the records do not fit together."""


class _Listed(TypedDict):
    batches: list[str]  # the names of the batch files that contain the task


class ListedComputation(ComputationTask, _Listed):
    pass


class ListedSynthetic(SyntheticTask, _Listed):
    pass


# every task of the batch files by id
Tasks = dict[str, ListedComputation | ListedSynthetic]


# ---------------------------------------------------------------------------------------------------
# reading


def load_registry() -> dict[str, Any]:
    with REGISTRY.open("rb") as handle:
        registry: dict[str, Any] = tomllib.load(handle)
    return registry


def batch_tasks() -> Tasks:
    """Every task of every batch file by id, each with the names of the batches that contain it."""
    tasks: Tasks = {}
    for path in sorted((REPOSITORY / "batches").glob("*.toml")):
        for task in plan_batch(path).tasks:
            listed = cast("ListedComputation | ListedSynthetic", {**task, "batches": []})
            tasks.setdefault(task["id"], listed)["batches"].append(path.stem)
    return tasks


def latest_records(root: Path) -> dict[str, dict[str, Any]]:
    """The most recent record of each task over the run directories under `root`."""
    latest: dict[str, dict[str, Any]] = {}
    if not root.is_dir():
        return latest
    for run_dir in sorted(path for path in root.iterdir() if (path / "results.jsonl").is_file()):
        for task_id, record in ResultStore(run_dir / "results.jsonl").latest().items():
            known = latest.get(task_id)
            if known is None or record["started_at"] > known["started_at"]:
                latest[task_id] = record | {"batch": run_dir.name}
    return latest


def standing(task: ListedComputation | ListedSynthetic, record: dict[str, Any] | None) -> str:
    """Where a task stands in the tables: the status of its latest record, if that record still answers it."""
    if record is None:
        return NOT_RERUN
    if not answers(
        record, {key: value for key, value in task.items() if key != "batches"}, ignore=NOT_COMPUTED
    ):
        return STALE
    return str(record["status"])


def _natural(text: str) -> list[Any]:
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text)]


def counted(task: ListedComputation | ListedSynthetic) -> TypeGuard[ListedComputation]:
    """Whether a task can count for a claim: negative controls and synthetic tasks never do."""
    return task["kind"] == "computation" and not task.get("control")


def tasks_of(claim: dict[str, Any], tasks: Tasks) -> list[ListedComputation]:
    """The tasks of a claim: those of its groups and those its patterns match."""
    groups = set(claim.get("groups", []))
    patterns = claim.get("tasks", [])
    chosen = [
        task
        for task_id, task in tasks.items()
        if counted(task)
        and (task["group"] in groups or any(fnmatch.fnmatchcase(task_id, pattern) for pattern in patterns))
    ]
    return sorted(chosen, key=lambda task: _natural(task["id"]))


def claim_tasks(claim: dict[str, Any], tasks: Tasks) -> list[str]:
    """The ids of the tasks of a claim."""
    return [task["id"] for task in tasks_of(claim, tasks)]


def validate(registry: dict[str, Any], tasks: Tasks) -> None:
    problems = []
    claims = {claim["id"]: claim for claim in registry["claim"]}
    ids = [claim["id"] for claim in registry["claim"]]
    problems += [f"two claims have the id {i}" for i in sorted({i for i in ids if ids.count(i) > 1})]
    groups = {task["group"] for task in tasks.values() if counted(task)}
    claimed: set[str] = set()
    for claim in registry["claim"]:
        name = claim["id"]
        table = registry["table"].get(claim["table"])
        if table is None:
            problems.append(f"{name}: no table {claim['table']!r}")
        elif len(claim["cells"]) != len(table["columns"]):
            problems.append(f"{name}: {len(claim['cells'])} cells for {len(table['columns'])} columns")
        selects = bool(claim.get("groups")) or bool(claim.get("tasks"))
        if selects == ("no_task" in claim):
            problems.append(f"{name}: a claim has either tasks (groups, tasks) or a reason under no_task")
        problems += [
            f"{name}: no task has the group {g!r}" for g in claim.get("groups", []) if g not in groups
        ]
        for pattern in claim.get("tasks", []):
            if re.search(r"\d\*", pattern):
                problems.append(
                    f"{name}: in {pattern!r} the * follows a digit and also matches longer numbers"
                )
            if not any(counted(tasks[task_id]) for task_id in fnmatch.filter(tasks, pattern)):
                problems.append(f"{name}: no task of any batch matches {pattern!r}")
        listed = claims.get(claim.get("records", name))
        if listed is None or "records" in listed:
            problems.append(f"{name}: records must name a claim that lists its own tasks")
        elif not set(claim_tasks(claim, tasks)) <= set(claim_tasks(listed, tasks)):
            problems.append(f"{name}: the claim {listed['id']} does not list all of its tasks")
        for item in claim.get("evidence", []):
            if set(item) - {"file", "line", "lines", "text"} or "file" not in item or "text" not in item:
                problems.append(f"{name}: evidence needs a file, a line or lines, and a text")
            elif ("line" in item) == ("lines" in item):
                problems.append(f"{name}: evidence in {item['file']} needs either a line or lines")
        claimed.update(claim_tasks(claim, tasks))
    for task_id, task in tasks.items():
        if counted(task) and task_id not in claimed:
            problems.append(f"the task {task_id} ({', '.join(task['batches'])}) belongs to no claim")
    if problems:
        raise RegistryError("\n".join(problems))


def evidence_problems(registry: dict[str, Any], data: Path) -> list[str] | None:
    """What is wrong with the stored lines named as evidence; None if the data archive is not unpacked."""
    files: dict[str, list[str]] = {}
    problems = []
    for claim in registry["claim"]:
        for item in claim.get("evidence", []):
            path = data / item["file"]
            if not path.is_file():
                return None
            if item["file"] not in files:
                files[item["file"]] = path.read_text(encoding="utf-8", errors="replace").splitlines()
            lines = files[item["file"]]
            first, last = (item["line"], item["line"]) if "line" in item else item["lines"]
            if not 1 <= first <= last <= len(lines):
                problems.append(f"{claim['id']}: {item['file']} has no lines {first} to {last}")
            elif item["text"] not in lines[first - 1]:
                problems.append(
                    f"{claim['id']}: line {first} of {item['file']} does not contain {item['text']!r}"
                )
    return problems


# ---------------------------------------------------------------------------------------------------
# what the scripts share


def imports(path: Path) -> set[str]:
    """The top-level modules a script imports, anywhere in its source."""
    found: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            found.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            found.add(node.module.split(".")[0])
    return found


def shared_rows() -> list[list[str]]:
    """Per script of experiments/: the other scripts it imports, other packages, stored files it names."""
    scripts = sorted(SCRIPTS.glob("*.py"))
    names = {path.stem for path in scripts}
    rows = []
    for path in scripts:
        found = imports(path)
        stored = sorted(set(STORED_NAME.findall(path.read_text(encoding="utf-8"))))
        rows.append([
            f"`{path.name}`",
            ", ".join(f"`{name}`" for name in sorted(found & names)) or "—",
            ", ".join(f"`{name}`" for name in sorted(found - names - sys.stdlib_module_names)) or "—",
            str(len(stored)) if stored else "—",
        ])  # fmt: skip
    return rows


# ---------------------------------------------------------------------------------------------------
# writing


def cell(value: Any, limit: int = 160) -> str:
    """A value as it appears in a table cell."""
    if value is None:
        text = ""
    elif isinstance(value, bool):
        text = "yes" if value else "no"
    elif isinstance(value, dict):
        text = ", ".join(f"{key} {cell(item, 10**9)}" for key, item in value.items())
    elif isinstance(value, list):
        text = "(" + ", ".join(cell(item, 10**9) for item in value) + ")"
    elif isinstance(value, int):
        text = f"{value:,}"
    else:
        text = str(value)
    text = text.replace("|", "\\|")
    return text if len(text) <= limit else text[: limit - 1] + "…"


def table(header: list[str], rows: list[list[str]]) -> list[str]:
    return [
        "| " + " | ".join(header) + " |",
        "|" + "---|" * len(header),
        *("| " + " | ".join(row) + " |" for row in rows),
    ]


def anchor(claim_id: str) -> str:
    return "claim-" + re.sub(r"[^a-z0-9]+", "-", claim_id.lower()).strip("-")


def stamp(records: Sequence[Mapping[str, Any]]) -> str:
    """When and from which commits a set of records comes."""
    days = sorted({record["started_at"][:10] for record in records})
    commits = sorted({(record.get("code") or {}).get("commit") or "unknown commit" for record in records})
    when = days[0] if len(days) == 1 else f"{days[0]} to {days[-1]}"
    return f"{when}, " + ", ".join(f"`{commit[:7]}`" for commit in commits)


def _plural(count: int, word: str) -> str:
    return word if count == 1 else word + "s"


def summary(
    claim: dict[str, Any],
    registry: dict[str, Any],
    tasks: Tasks,
    latest: dict[str, dict[str, Any]],
) -> str:
    """The last cell of a claim's row: what the scheduler's records say about it."""
    chosen = claim_tasks(claim, tasks)
    if not chosen:
        return f"no task: {claim['no_task']}"
    covers = claim.get("covers") or registry["table"][claim["table"]].get("covers")
    scope = f" Covers: {covers}." if covers else ""
    counts = Counter(standing(tasks[task_id], latest.get(task_id)) for task_id in chosen)
    recorded = [
        latest[task_id]
        for task_id in chosen
        if standing(tasks[task_id], latest.get(task_id)) not in (STALE, NOT_RERUN)
    ]
    link = f"[tasks](claims/RECORDS.md#{anchor(claim.get('records', claim['id']))})"
    if counts[NOT_RERUN] == len(chosen):
        batches = sorted({batch for task_id in chosen for batch in tasks[task_id]["batches"]})
        where = ", ".join(f"`{batch}`" for batch in batches)
        return f"**{NOT_RERUN}**: {len(chosen)} {_plural(len(chosen), 'task')} in {where}, no record; {link}.{scope}"
    if counts["ok"] == len(chosen):
        words = Counter(expectation.verdict(record) for record in recorded)
        how = ", ".join(
            f"{words[word]} {text}" for word, text in expectation.ENDED_WELL.items() if words[word]
        )
        lead = f"all {len(chosen)} tasks **ok** ({how})" if len(chosen) > 1 else f"its task is **ok** ({how})"
    else:
        parts = [f"{status} {counts[status]}" for status in STATUSES if counts[status]]
        lead = "tasks: " + "; ".join(parts) + f" (of {len(chosen)})"
    when = f"{stamp(recorded)}; " if recorded else ""
    return f"{lead}; {when}{link}.{scope}"


def evidence_cell(claim: dict[str, Any]) -> str:
    parts = []
    for item in claim.get("evidence", []):
        where = f"`{item['file']}`, "
        where += (
            f"line {item['line']}" if "line" in item else f"lines {item['lines'][0]} to {item['lines'][1]}"
        )
        parts.append(where)
    return "; ".join(parts)


def claims_table(
    slug: str, registry: dict[str, Any], tasks: Tasks, latest: dict[str, dict[str, Any]]
) -> list[str]:
    spec = registry["table"][slug]
    rows = [
        [*claim["cells"], evidence_cell(claim), summary(claim, registry, tasks, latest)]
        for claim in registry["claim"]
        if claim["table"] == slug
    ]
    return table([*spec["columns"], "Stored lines", "Re-run under the scheduler"], rows)


def plain_table(slug: str, registry: dict[str, Any]) -> list[str]:
    """A table of the registry as it was written by hand, without the generated columns."""
    rows = [claim["cells"] for claim in registry["claim"] if claim["table"] == slug]
    return table(registry["table"][slug]["columns"], rows)


def compared_cell(task: ListedComputation) -> str:
    """What a task's output is compared with: what the task predicts, in words."""
    modes = {item["file"]: item["mode"] for item in task["compared"]}
    parts = []
    for key, value in (task["predicted"] or {}).items():
        if key in modes:
            parts.append(f"`{key}` ({modes[key]}): {value}")
        else:
            parts.append(f"{cell(key)}: {cell(value)}")
    return "; ".join(parts) or "nothing: its output was not kept"


def record_cells(task: ListedComputation, record: dict[str, Any] | None) -> list[str]:
    if record is None:
        return [NOT_RERUN, ""]
    commit = ((record.get("code") or {}).get("commit") or "unknown")[:7]
    if standing(task, record) == STALE:
        return [
            f"{STALE}: recorded {record['started_at'][:10]} at `{commit}` for the task as it was then",
            "",
        ]
    state = f"{record['status']}, {record['started_at'][:10]}, `{commit}`, {record['wall_seconds']:g} s"
    differing = expectation.differing(record)
    found = expectation.verdict(record) + (f": {cell(differing)}" if differing else "")
    return [state, found]


def detail(registry: dict[str, Any], tasks: Tasks, latest: dict[str, dict[str, Any]]) -> str:
    """claims/RECORDS.md: every claim with its tasks, one row per task."""
    lines = [
        "# Tasks and records, claim by claim",
        "",
        "Generated by `tools/claims.py tables`; do not edit. For every claim of [CLAIMS.md](../CLAIMS.md) that has tasks: "
        "each task with what its output is compared with, and the latest record of the task under the scheduler, if "
        f'there is one. A task without a record reads "{NOT_RERUN}": the stored result is then the only support, and it '
        f'is the original output, not a run of this repository. A record reads "{STALE}" when its task has changed '
        "since it was computed (another command, or other stored files or lines to compare with): it supports nothing "
        "until the task is run again.",
        "",
        'A task is one command of a note\'s *Reproduce* section. "Compared with" names the stored files of the data '
        "archive ([DATA.md](../DATA.md)) and the mode of each comparison ([GUIDE.md](../GUIDE.md#what-a-task-does)), "
        "or the lines its output must contain. The verdict is **reproduced** when everything compared is identical, "
        "**checked** when there is no stored file but the output has the lines it must have, **ran** when the command "
        "ended and nothing was compared, and **differs** otherwise. A record names the day it was computed and the "
        "commit of the code.",
        "",
    ]
    for claim in registry["claim"]:
        chosen = tasks_of(claim, tasks)
        if not chosen or "records" in claim:
            continue
        lines += [f'<a name="{anchor(claim["id"])}"></a>', "", f"## {claim['id']}", "", claim["cells"][0], ""]
        rows = [
            [
                f"`{task['id']}`",
                ", ".join(task["batches"]),
                compared_cell(task),
                *record_cells(task, latest.get(task["id"])),
                cell(task["caveat"], 400),
            ]
            for task in chosen
        ]
        lines += table(["Task", "Batches", "Compared with", "Record", "Verdict", "Caveat"], rows)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def blocks(
    registry: dict[str, Any], tasks: Tasks, latest: dict[str, dict[str, Any]]
) -> dict[Path, dict[str, list[str]]]:
    """The generated parts of the hand-written files, by file and name."""
    shared = table(
        ["Script", "Imports from `experiments/`", "Other packages", "Stored files it names"], shared_rows()
    )
    return {
        REPOSITORY / "CLAIMS.md": {
            **{slug: claims_table(slug, registry, tasks, latest) for slug in registry["table"]},
            "shared": shared,
        },
        REPOSITORY / "GUIDE.md": {"main": plain_table("main", registry)},
    }


def splice(text: str, name: str, lines: list[str], path: Path) -> str:
    begin, end = f"<!-- BEGIN GENERATED {name}: tools/claims.py tables -->", f"<!-- END GENERATED {name} -->"
    if text.count(begin) != 1 or text.count(end) != 1:
        raise RegistryError(f"{path.name} needs exactly one pair of markers for the generated part {name!r}")
    head, rest = text.split(begin)
    _, tail = rest.split(end)
    return head + begin + "\n" + "\n".join(lines) + "\n" + end + tail


def outputs() -> dict[Path, str]:
    registry, tasks, latest = load_registry(), batch_tasks(), latest_records(RECORDS)
    validate(registry, tasks)
    written = {DETAIL: detail(registry, tasks, latest)}
    for path, parts in blocks(registry, tasks, latest).items():
        text = path.read_text(encoding="utf-8")
        for name, lines in parts.items():
            text = splice(text, name, lines, path)
        written[path] = text
    return written


# ---------------------------------------------------------------------------------------------------
# commands


def command_tables() -> int:
    for path, text in outputs().items():
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(REPOSITORY)}")
    return 0


def command_check() -> int:
    stale = [
        path
        for path, text in outputs().items()
        if not path.exists() or path.read_text(encoding="utf-8") != text
    ]
    for path in stale:
        print(
            f"{path.relative_to(REPOSITORY)} is not what the registry and the records give: run `make tables`",
            file=sys.stderr,
        )
    if not stale:
        print("the generated tables agree with the claim registry, the batch files and the records")
    data = Path(os.environ.get("K3_DATA_DIR", SCRIPTS))
    wrong = evidence_problems(load_registry(), data)
    if wrong is None:
        print(f"stored lines named as evidence: not checked, the data archive is not unpacked in {data}")
    elif wrong:
        print("\n".join(wrong), file=sys.stderr)
    else:
        print("stored lines named as evidence: each is where the registry says")
    return 1 if stale or wrong else 0


def command_adopt(run_dirs: list[Path]) -> int:
    for run_dir in run_dirs:
        done = adopt.adopt(run_dir, RECORDS, REPOSITORY)
        print(f"records/{run_dir.name}: {len(done.records)} records, {stamp(done.records)}")
        for (old, now), count in sorted(done.moved.items()):
            print(f"   {count} records: commit {old[:7]} is {now[:7]} in the present history (same tree)")
    return 0


def main(arguments: list[str]) -> int:
    try:
        if arguments == ["tables"]:
            return command_tables()
        if arguments == ["check"]:
            return command_check()
        if arguments == ["words"]:
            print(*adopt.private_words(REPOSITORY), sep="\n")
            return 0
        if len(arguments) >= 2 and arguments[0] == "adopt":
            return command_adopt([Path(argument) for argument in arguments[1:]])
    except (RegistryError, ValueError, OSError) as problem:
        print(f"claims.py: {problem}", file=sys.stderr)
        return 1
    print((__doc__ or "").split("\n\n")[1], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
