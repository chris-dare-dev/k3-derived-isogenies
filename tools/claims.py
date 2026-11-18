#!/usr/bin/env python3
"""claims.py -- the claim tables, written from the claim registry, the batch files and the records.

Usage:  python3 tools/claims.py tables              write the generated parts of CLAIMS.md and README.md,
                                                     and claims/RECORDS.md
        python3 tools/claims.py check               fail if they are not what `tables` writes, if the claim
                                                     registry and the batch files do not fit together, or if
                                                     a stored line named as evidence is not what it says
        python3 tools/claims.py adopt RUN_DIR ...   copy the records of finished runs (results/<batch>)
                                                     into records/<batch>

What is read:  claims/claims.toml (the claims and the tasks that bear on them), batches/*.toml (which tasks
exist), records/<batch>/results.jsonl (what was computed here), the import statements of experiments/*.py
(what the scripts share), and, when the data archive is unpacked, the stored lines named as evidence.
A task without a record reads "not re-run", whatever the stored results say.

`adopt` copies results.jsonl and manifest.json unchanged, and writes adopted.json next to them: the commit
of each record with the hash of its tree.  It refuses records that were computed from a working tree with
uncommitted changes, and records that name this machine or a directory outside the repository.  One field
can change on the way: if the history was rewritten after the run (messages only, every tree as it was),
a record's commit is replaced by the commit that now has the same tree, and `adopt` says so.
"""

from __future__ import annotations

import ast
import fnmatch
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tomllib
from collections import Counter
from pathlib import Path
from typing import Any

REPOSITORY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPOSITORY))

from k3tasks.report import verdict  # noqa: E402
from k3tasks.tasks import expand_batch  # noqa: E402
from taskbatch.config import load_config  # noqa: E402
from taskbatch.store import ResultStore  # noqa: E402

REGISTRY = REPOSITORY / "claims" / "claims.toml"
RECORDS = REPOSITORY / "records"
DETAIL = REPOSITORY / "claims" / "RECORDS.md"
SCRIPTS = REPOSITORY / "experiments"
NOT_RERUN = "not re-run"
STATUSES = ["ok", "mismatch", "unsupported", "timeout", "memory", "error", "crashed", NOT_RERUN]
VERDICTS = {
    "reproduced": "reproduced",
    "checked": "checked by the lines of their output",
    "ran": "ran with nothing to compare",
}
# kinds of task that never count for a claim
NOT_CLAIMED = {"synthetic"}
STORED_NAME = re.compile(r"""["']((?:sixlines|ellk3)_[A-Za-z0-9_]+\.(?:jsonl|json|log))["']""")


class RegistryError(ValueError):
    """The registry, the batch files and the records do not fit together."""


# ---------------------------------------------------------------------------------------------------
# reading


def load_registry() -> dict[str, Any]:
    with REGISTRY.open("rb") as handle:
        registry: dict[str, Any] = tomllib.load(handle)
    return registry


def batch_tasks() -> dict[str, dict[str, Any]]:
    """Every task of every batch file by id, each with the names of the batches that contain it."""
    tasks: dict[str, dict[str, Any]] = {}
    for path in sorted((REPOSITORY / "batches").glob("*.toml")):
        for task in expand_batch(load_config(path).description):
            tasks.setdefault(task["id"], task | {"batches": []})["batches"].append(path.stem)
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


def _natural(text: str) -> list[Any]:
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text)]


def counted(task: dict[str, Any]) -> bool:
    """Whether a task can count for a claim: negative controls and synthetic tasks never do."""
    return not task.get("control") and task["kind"] not in NOT_CLAIMED


def claim_tasks(claim: dict[str, Any], tasks: dict[str, dict[str, Any]]) -> list[str]:
    """The ids of the tasks of a claim: those of its groups and those its patterns match."""
    groups = set(claim.get("groups", []))
    patterns = claim.get("tasks", [])
    chosen = {
        task_id
        for task_id, task in tasks.items()
        if counted(task)
        and (task["group"] in groups or any(fnmatch.fnmatchcase(task_id, pattern) for pattern in patterns))
    }
    return sorted(chosen, key=_natural)


def validate(registry: dict[str, Any], tasks: dict[str, dict[str, Any]]) -> None:
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


def stamp(records: list[dict[str, Any]]) -> str:
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
    tasks: dict[str, dict[str, Any]],
    latest: dict[str, dict[str, Any]],
) -> str:
    """The last cell of a claim's row: what the scheduler's records say about it."""
    chosen = claim_tasks(claim, tasks)
    if not chosen:
        return f"no task: {claim['no_task']}"
    covers = claim.get("covers") or registry["table"][claim["table"]].get("covers")
    scope = f" Covers: {covers}." if covers else ""
    recorded = [latest[task_id] for task_id in chosen if task_id in latest]
    counts = Counter(record["status"] for record in recorded)
    counts[NOT_RERUN] = len(chosen) - len(recorded)
    link = f"[tasks](claims/RECORDS.md#{anchor(claim.get('records', claim['id']))})"
    if not recorded:
        batches = sorted({batch for task_id in chosen for batch in tasks[task_id]["batches"]})
        where = ", ".join(f"`{batch}`" for batch in batches)
        return f"**{NOT_RERUN}**: {len(chosen)} {_plural(len(chosen), 'task')} in {where}, no record; {link}.{scope}"
    if counts["ok"] == len(chosen):
        words = Counter(verdict(record) for record in recorded)
        how = ", ".join(f"{words[word]} {text}" for word, text in VERDICTS.items() if words[word])
        lead = f"all {len(chosen)} tasks **ok** ({how})" if len(chosen) > 1 else f"its task is **ok** ({how})"
    else:
        parts = [f"{status} {counts[status]}" for status in STATUSES if counts[status]]
        lead = "tasks: " + "; ".join(parts) + f" (of {len(chosen)})"
    return f"{lead}; {stamp(recorded)}; {link}.{scope}"


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
    slug: str, registry: dict[str, Any], tasks: dict[str, dict[str, Any]], latest: dict[str, dict[str, Any]]
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


def compared_cell(task: dict[str, Any]) -> str:
    """What a task's output is compared with: the predicted part of its record, in words."""
    modes = {item["file"]: item["mode"] for item in task.get("compared", [])}
    parts = []
    for key, value in (task.get("predicted") or {}).items():
        if key in modes:
            parts.append(f"`{key}` ({modes[key]}): {value}")
        else:
            parts.append(f"{cell(key)}: {cell(value)}")
    return "; ".join(parts) or "nothing: its output was not kept"


def record_cells(record: dict[str, Any] | None) -> list[str]:
    if record is None:
        return [NOT_RERUN, ""]
    commit = ((record.get("code") or {}).get("commit") or "unknown")[:7]
    state = f"{record['status']}, {record['started_at'][:10]}, `{commit}`, {record['wall_seconds']:g} s"
    differing = {
        key: value
        for key, value in (record.get("value") or {}).items()
        if value != (record.get("predicted") or {}).get(key)
    }
    found = verdict(record) + (f": {cell(differing)}" if differing else "")
    return [state, found]


def detail(
    registry: dict[str, Any], tasks: dict[str, dict[str, Any]], latest: dict[str, dict[str, Any]]
) -> str:
    """claims/RECORDS.md: every claim with its tasks, one row per task."""
    lines = [
        "# Tasks and records, claim by claim",
        "",
        "Generated by `tools/claims.py tables`; do not edit. For every claim of [CLAIMS.md](../CLAIMS.md) that has tasks: "
        "each task with what its output is compared with, and the latest record of the task under the scheduler, if "
        f'there is one. A task without a record reads "{NOT_RERUN}": the stored result is then the only support, and it '
        "is the original output, not a run of this repository.",
        "",
        'A task is one command of a note\'s *Reproduce* section. "Compared with" names the stored files of the data '
        "archive ([DATA.md](../DATA.md)) and the mode of each comparison ([README.md](../README.md#what-a-task-does)), "
        "or the lines its output must contain. The verdict is **reproduced** when everything compared is identical, "
        "**checked** when there is no stored file but the output has the lines it must have, **ran** when the command "
        "ended and nothing was compared, and **differs** otherwise. A record names the day it was computed and the "
        "commit of the code.",
        "",
    ]
    for claim in registry["claim"]:
        chosen = claim_tasks(claim, tasks)
        if not chosen or "records" in claim:
            continue
        lines += [f'<a name="{anchor(claim["id"])}"></a>', "", f"## {claim['id']}", "", claim["cells"][0], ""]
        rows = [
            [
                f"`{task_id}`",
                ", ".join(tasks[task_id]["batches"]),
                compared_cell(tasks[task_id]),
                *record_cells(latest.get(task_id)),
                cell(tasks[task_id].get("caveat"), 400),
            ]
            for task_id in chosen
        ]
        lines += table(["Task", "Batches", "Compared with", "Record", "Verdict", "Caveat"], rows)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def blocks(
    registry: dict[str, Any], tasks: dict[str, dict[str, Any]], latest: dict[str, dict[str, Any]]
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
        REPOSITORY / "README.md": {"main": plain_table("main", registry)},
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


def _git(*arguments: str) -> str | None:
    done = subprocess.run(
        ["git", "-C", str(REPOSITORY), *arguments], capture_output=True, text=True, check=False
    )
    return done.stdout.strip() if done.returncode == 0 else None


def renumbered(commit: str) -> str:
    """The commit of the present history that a record's commit stands for.

    A record names the commit of the code that computed it.  If the history was rewritten since (its
    messages changed, its trees not), that commit is no longer in it; the commit that now has the same
    tree is the same code under a new number.  Any other case is refused.
    """
    if _git("merge-base", "--is-ancestor", commit, "HEAD") is not None:
        return commit
    tree = _git("rev-parse", f"{commit}^{{tree}}")
    if tree is None:
        raise RegistryError(f"the commit {commit[:7]} of a record is not in this repository")
    history = [line.split() for line in (_git("log", "--format=%H %T", "HEAD") or "").splitlines()]
    same = [number for number, its_tree in history if its_tree == tree]
    if len(same) != 1:
        raise RegistryError(
            f"the commit {commit[:7]} of a record is not in the history, and {len(same)} commits of the history "
            "have its tree: the record cannot be placed"
        )
    return same[0]


def private_words() -> list[str]:
    """Text that a committed record must not contain: this machine's name and directories outside the tree."""
    words = [str(Path.home()), str(REPOSITORY)]
    host = socket.gethostname()
    return [word for word in [*words, host] if len(word) > 3]


def command_adopt(run_dirs: list[Path]) -> int:
    for run_dir in run_dirs:
        records = list(ResultStore(run_dir / "results.jsonl").latest().values())
        if not records:
            raise RegistryError(f"no records in {run_dir}")
        dirty = [r["task_id"] for r in records if (r.get("code") or {}).get("dirty") is not False]
        if dirty:
            raise RegistryError(
                f"{run_dir}: {len(dirty)} records, among them {dirty[0]}, were computed from a working tree with "
                "uncommitted changes or outside a checkout; commit the code and run the batch again"
            )
        text = (run_dir / "results.jsonl").read_text(encoding="utf-8")
        manifest = (run_dir / "manifest.json").read_text(encoding="utf-8")
        for word in private_words():
            if word.lower() in (text + manifest).lower():
                raise RegistryError(
                    f"{run_dir}: a record names this machine or a directory outside the repository; "
                    "records are committed, so it is not adopted"
                )
        # every line is copied as it is, unless its commit was renumbered by a rewrite of the history:
        # then that one field changes, and nothing else
        moved: Counter[tuple[str, str]] = Counter()
        trees: dict[str, str] = {}
        lines = []
        for line in text.splitlines():
            try:
                record = json.loads(line)
                commit = record["code"]["commit"]
            except (json.JSONDecodeError, KeyError, TypeError):
                lines.append(line)
                continue
            now = renumbered(commit)
            if now != commit:
                record["code"]["commit"] = now
                moved[commit, now] += 1
                line = json.dumps(record, sort_keys=True, ensure_ascii=False, allow_nan=False)
            trees[now] = _git("rev-parse", f"{now}^{{tree}}") or "unknown"
            lines.append(line)
        target = RECORDS / run_dir.name
        target.mkdir(parents=True, exist_ok=True)
        (target / "results.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")
        shutil.copyfile(run_dir / "manifest.json", target / "manifest.json")
        sidecar = {"trees": dict(sorted(trees.items()))}
        (target / "adopted.json").write_text(json.dumps(sidecar, indent=1) + "\n", encoding="utf-8")
        adopted = list(ResultStore(target / "results.jsonl").latest().values())
        print(f"records/{run_dir.name}: {len(adopted)} records, {stamp(adopted)}")
        for (old, now), count in sorted(moved.items()):
            print(f"   {count} records: commit {old[:7]} is {now[:7]} in the present history (same tree)")
    return 0


def main(arguments: list[str]) -> int:
    try:
        if arguments == ["tables"]:
            return command_tables()
        if arguments == ["check"]:
            return command_check()
        if len(arguments) >= 2 and arguments[0] == "adopt":
            return command_adopt([Path(argument) for argument in arguments[1:]])
    except (RegistryError, ValueError, OSError) as problem:
        print(f"claims.py: {problem}", file=sys.stderr)
        return 1
    print((__doc__ or "").split("\n\n")[1], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
