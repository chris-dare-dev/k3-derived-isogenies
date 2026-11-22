"""The registry of claims and the tables generated from it."""

from __future__ import annotations

import json
import re
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path
from typing import Any

from taskbatch.settings import REPOSITORY
from taskbatch.store import ResultStore

from . import helpers

sys.path.insert(0, str(REPOSITORY / "src" / "bookkeeping"))

import claims


def record(
    task_id: str,
    status: str = "ok",
    day: str = "2026-11-17",
    verdict: str = "reproduced",
    **payload: Any,
) -> dict[str, Any]:
    """The record of a task of TASKS as it is now; `payload` changes the task it was computed for."""
    known = {key: value for key, value in TASKS.get(task_id, {"kind": "k"}).items() if key != "batches"}
    return dict(
        helpers.record(
            known | {"id": task_id} | payload,
            status,
            value={"f.log": "identical" if status == "ok" else "differs"},
            predicted={"f.log": "identical"},
            agrees=status == "ok",
            wall_seconds=0.5,
            started_at=f"{day}T10:00:00+00:00",
            detail={"verdict": verdict if status == "ok" else "differs"},
            code={"commit": "a" * 40, "dirty": False},
        )
    )


def task(task_id: str, group: str, batch: str, **more: Any) -> dict[str, Any]:
    return {
        "id": task_id, "kind": "computation", "group": group, "batches": [batch],
        "predicted": {"f.log": "identical"}, "compared": [{"file": "f.log", "mode": "text"}], "caveat": None,
    } | more  # fmt: skip


TASKS = {
    "n02/a": task("n02/a", "lattices", "lattices"),
    "n02/b": task("n02/b", "lattices", "heavy", caveat="its output was not kept"),
    "n27/ledger": task("n27/ledger", "ledger", "smoke"),
    "n27/ledger!control": task("n27/ledger!control", "negative control", "negative-control", control=True),
    "synthetic/sleeps": task("synthetic/sleeps", "synthetic", "limits-demo") | {"kind": "synthetic"},
}
REGISTRY = {
    "table": {"t": {"columns": ["Result"]}, "l": {"columns": ["ID"], "covers": "the recount"}},
    "claim": [
        {"id": "lattices", "table": "t", "cells": ["Lattices"], "groups": ["lattices"]},
        {"id": "ledger", "table": "t", "cells": ["Ledger"], "groups": ["ledger"]},
        {"id": "H1", "table": "l", "cells": ["H1"], "tasks": ["n27/ledger"], "records": "ledger"},
        {"id": "H29", "table": "l", "cells": ["H29"], "no_task": "no section"},
    ],
}
LATTICES, LEDGER, H1, H29 = REGISTRY["claim"]


def summary(claim: dict[str, Any], latest: dict[str, dict[str, Any]]) -> str:
    return claims.summary(claim, REGISTRY, TASKS, latest)


class Registry(unittest.TestCase):
    def test_the_registry_fits_the_batch_files(self) -> None:
        registry, tasks = claims.load_registry(), claims.batch_tasks()
        claims.validate(registry, tasks)
        self.assertEqual(len(registry["claim"]), 58)
        # every computation of every batch belongs to a claim; controls and synthetic tasks to none
        claimed = {t for claim in registry["claim"] for t in claims.claim_tasks(claim, tasks)}
        self.assertEqual(len(claimed), 163)
        self.assertFalse([t for t in claimed if t.endswith("!control") or t.startswith("synthetic/")])

    def test_the_main_rows_are_the_groups_of_the_task_registry(self) -> None:
        registry, tasks = claims.load_registry(), claims.batch_tasks()
        rows = {c["id"]: c for c in registry["claim"] if c["table"] in ("main", "other")}
        self.assertEqual(len(rows), 11)
        self.assertEqual(sum(len(claims.claim_tasks(c, tasks)) for c in rows.values()), 163)
        self.assertIn("no_task", rows["torsors"])

    def test_the_ledger_rows_are_those_of_note_27_word_for_word(self) -> None:
        note = (REPOSITORY / "research" / "docs" / "27-hypothesis-ledger.md").read_text(encoding="utf-8")
        rows = {}
        for line in note.splitlines():
            if re.match(r"\| H\d+ \|", line):
                cells = [
                    cell.strip() for cell in re.split(r"(?<!\\) \| ", line.strip()[2:].removesuffix(" |"))
                ]
                rows[cells[0]] = cells[:4]
        self.assertEqual(len(rows), 47)
        ledger = {c["id"]: c for c in claims.load_registry()["claim"] if c["table"] == "ledger"}
        self.assertEqual(set(ledger), set(rows))
        for name, claim in ledger.items():
            with self.subTest(name):
                self.assertEqual(claim["cells"], rows[name])
                # H1 to H28 are the sections of the ledger script; the others have no computation
                self.assertEqual("no_task" in claim, int(name[1:]) > 28)

    def test_what_does_not_fit_is_reported(self) -> None:
        for change, message in (
            ({"groups": ["curves"]}, "no task has the group 'curves'"),
            ({"groups": [], "tasks": ["n99/*"]}, "no task of any batch matches 'n99/*'"),
            ({"cells": ["a", "b"]}, "2 cells for 1 columns"),
            ({"table": "missing"}, "no table 'missing'"),
            ({"no_task": "none"}, "either tasks (groups, tasks) or a reason"),
            ({"records": "H1"}, "records must name a claim that lists its own tasks"),
            ({"records": "ledger"}, "the claim ledger does not list all of its tasks"),
            ({"evidence": [{"file": "f.log", "text": "x"}]}, "needs either a line or lines"),
            ({"evidence": [{"file": "f.log", "line": 1}]}, "evidence needs a file"),
        ):
            registry = REGISTRY | {"claim": [LATTICES | change, LEDGER, H1, H29]}
            with self.subTest(change=change), self.assertRaises(claims.RegistryError) as caught:
                claims.validate(registry, TASKS)
            self.assertIn(message, str(caught.exception))
        with self.assertRaisesRegex(claims.RegistryError, "n27/ledger .smoke. belongs to no claim"):
            claims.validate(REGISTRY | {"claim": [LATTICES, H29]}, TASKS)
        with self.assertRaisesRegex(claims.RegistryError, "two claims have the id lattices"):
            claims.validate(REGISTRY | {"claim": [LATTICES, LATTICES, LEDGER, H1, H29]}, TASKS)
        with self.assertRaisesRegex(claims.RegistryError, "follows a digit"):
            claims.validate(REGISTRY | {"claim": [LATTICES, LEDGER | {"tasks": ["n2*"]}, H1, H29]}, TASKS)
        claims.validate(REGISTRY, TASKS)

    def test_controls_and_synthetic_tasks_never_count_for_a_claim(self) -> None:
        everything = {"id": "x", "table": "t", "cells": ["x"], "tasks": ["*"]}
        self.assertEqual(claims.claim_tasks(everything, TASKS), ["n02/a", "n02/b", "n27/ledger"])
        self.assertEqual(claims.claim_tasks(LEDGER, TASKS), ["n27/ledger"])

    def test_a_claim_without_a_record_reads_not_re_run(self) -> None:
        text = summary(LATTICES, {})
        self.assertTrue(text.startswith("**not re-run**: 2 tasks in `heavy`, `lattices`, no record"))
        self.assertNotIn("ok", text)

    def test_a_claim_is_ok_only_when_every_task_is(self) -> None:
        one = {"n02/a": record("n02/a")}
        self.assertIn("tasks: ok 1; not re-run 1 (of 2)", summary(LATTICES, one))
        self.assertNotIn("all", summary(LATTICES, one))
        both = one | {"n02/b": record("n02/b", verdict="ran")}
        self.assertIn(
            "all 2 tasks **ok** (1 reproduced, 1 ran with nothing to compare)", summary(LATTICES, both)
        )
        differs = one | {"n02/b": record("n02/b", status="mismatch")}
        self.assertIn("tasks: ok 1; mismatch 1 (of 2)", summary(LATTICES, differs))
        self.assertIn("2026-11-17, `aaaaaaa`", summary(LATTICES, both))

    def test_a_claim_without_tasks_says_why_and_a_table_says_what_its_tasks_cover(self) -> None:
        self.assertEqual(summary(H29, {}), "no task: no section")
        text = summary(H1, {"n27/ledger": record("n27/ledger")})
        self.assertIn("its task is **ok** (1 reproduced)", text)
        self.assertIn("[tasks](data/claims/RECORDS.md#claim-ledger)", text)
        self.assertTrue(text.endswith("Covers: the recount."))

    def test_the_detail_page_lists_comparison_record_and_caveat(self) -> None:
        page = claims.detail(REGISTRY, TASKS, {"n02/a": record("n02/a")})
        self.assertIn(
            "| `n02/a` | lattices | `f.log` (text): identical | ok, 2026-11-17, `aaaaaaa`, 0.5 s | reproduced |  |",
            page,
        )
        self.assertIn(
            "| `n02/b` | heavy | `f.log` (text): identical | not re-run |  | its output was not kept |", page
        )
        self.assertNotIn("## H1", page)  # listed under the claim it names
        self.assertNotIn("!control", page)
        page = claims.detail(REGISTRY, TASKS, {"n02/a": record("n02/a", status="mismatch")})
        self.assertIn("| differs: f.log differs |", page)

    def test_a_record_of_a_task_that_changed_since_is_stale_and_supports_nothing(self) -> None:
        old = record("n02/b", predicted={"f.log": "identical", "g.log": "identical"})
        self.assertEqual(claims.standing(TASKS["n02/b"], old), "stale")
        self.assertEqual(claims.standing(TASKS["n02/b"], record("n02/b", group="moved", set=3)), "ok")
        latest = {"n02/a": record("n02/a"), "n02/b": old}
        text = summary(LATTICES, latest)
        self.assertIn("tasks: ok 1; stale 1 (of 2)", text)
        self.assertNotIn("all", text)
        self.assertIn(
            "| `n02/b` | heavy | `f.log` (text): identical "
            "| stale: recorded 2026-11-17 at `aaaaaaa` for the task as it was then |  |",
            claims.detail(REGISTRY, TASKS, latest),
        )
        # a claim whose only records are stale names no day and no commit
        text = summary(LATTICES, {"n02/b": old})
        self.assertIn("tasks: stale 1; not re-run 1 (of 2); [tasks]", text)

    def test_a_changed_caveat_or_estimate_leaves_the_record_standing(self) -> None:
        task = TASKS["n02/b"]
        for changed in (
            {"caveat": "words written before the rerun"},
            {"note": 3},
            {"procs": 2},
            {"expected_memory_mb": 500},
            {"keep_work": True},
        ):
            with self.subTest(changed=changed):
                self.assertEqual(claims.standing(task, record("n02/b", **changed)), "ok")
        # what the task computes, and what it is compared with, still make a record stale
        for changed in (
            {"command": "python other.py"},
            {"compared": [{"file": "f.log", "mode": "exact"}]},
            {"predicted": {"f.log": "identical", "g.log": "identical"}},
        ):
            with self.subTest(changed=changed):
                self.assertEqual(claims.standing(task, record("n02/b", **changed)), "stale")

    def test_stored_lines_named_as_evidence_are_checked_against_the_data(self) -> None:
        item = {"file": "f.log", "line": 2, "text": "162 agree"}
        registry = REGISTRY | {"claim": [LEDGER | {"evidence": [item]}]}
        with tempfile.TemporaryDirectory() as scratch:
            data = Path(scratch)
            self.assertIsNone(claims.evidence_problems(registry, data))  # not unpacked: not checked
            (data / "f.log").write_text("ledger\n170 comparisons, 162 agree\n")
            self.assertEqual(claims.evidence_problems(registry, data), [])
            (data / "f.log").write_text("ledger\n170 comparisons, 161 agree\n")
            self.assertIn("line 2 of f.log does not contain", claims.evidence_problems(registry, data)[0])
            (data / "f.log").write_text("ledger\n")
            self.assertIn("has no lines 2 to 2", claims.evidence_problems(registry, data)[0])


class Files(unittest.TestCase):
    def test_the_committed_tables_are_what_the_generator_writes(self) -> None:
        stale = [
            path.name for path, text in claims.outputs().items() if path.read_text(encoding="utf-8") != text
        ]
        self.assertEqual(stale, [], "run `make tables`")

    def test_generated_parts_need_their_markers(self) -> None:
        begin, end = (
            "<!-- BEGIN GENERATED x: src/bookkeeping/claims.py tables -->",
            "<!-- END GENERATED x -->",
        )
        text = f"before\n{begin}\nold\n{end}\nafter\n"
        self.assertEqual(
            claims.splice(text, "x", ["new"], Path("f.md")), f"before\n{begin}\nnew\n{end}\nafter\n"
        )
        with self.assertRaises(claims.RegistryError):
            claims.splice("no markers", "x", ["new"], Path("f.md"))

    def test_the_adopted_controls_ended_as_they_must(self) -> None:
        expected = {
            "negative-control": {"mismatch": 3},
            "limits-demo": {"ok": 1, "timeout": 1, "memory": 1},
        }
        for batch, statuses in expected.items():
            path = claims.RECORDS / batch / "results.jsonl"
            if not path.exists():
                continue
            found: dict[str, int] = {}
            for row in ResultStore(path).latest().values():
                found[row["status"]] = found.get(row["status"], 0) + 1
            self.assertEqual(found, statuses, batch)

    def test_adopted_records_come_from_committed_code_and_name_no_machine(self) -> None:
        for path in sorted(claims.RECORDS.glob("*/results.jsonl")):
            text = path.read_text(encoding="utf-8")
            self.assertNotRegex(text, r"/home/|/Users/|/root/", path.parent.name)
            for line in text.splitlines():
                self.assertIs(json.loads(line)["code"]["dirty"], False, path.parent.name)
            trees = json.loads((path.parent / "adopted.json").read_text(encoding="utf-8"))["trees"]
            self.assertEqual(set(trees), {json.loads(line)["code"]["commit"] for line in text.splitlines()})


class Records(unittest.TestCase):
    def test_the_latest_record_of_a_task_wins_across_run_directories(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            for name, day, status in (("a", "2026-11-17", "timeout"), ("b", "2026-11-18", "ok")):
                (Path(scratch) / name).mkdir()
                (Path(scratch) / name / "results.jsonl").write_text(
                    json.dumps(record("t", status=status, day=day)) + "\n"
                )
            latest = claims.latest_records(Path(scratch))
            self.assertEqual((latest["t"]["status"], latest["t"]["batch"]), ("ok", "b"))


if __name__ == "__main__":
    unittest.main()
