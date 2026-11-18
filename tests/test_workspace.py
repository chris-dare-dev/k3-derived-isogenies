"""Running a computation in its own work directory, with toy scripts and toy stored results.

The toy script behaves like the scripts of experiments/: it finds its inputs next to itself, writes its
outputs there, and may append to a file that also exists among the stored results.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import textwrap
import unittest
from unittest import mock
from pathlib import Path

from k3tasks.registry import STDOUT, Check, Computation
from k3tasks.workspace import Places, WorkspaceError, evaluate, execute, finish, prepare

TOY_SCRIPT = textwrap.dedent(
    """
    import json, os, sys
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, HERE)
    import toy_helper
    mode = sys.argv[1]
    rows = [json.loads(line) for line in open(os.path.join(HERE, "toy_input.jsonl"))]
    if mode == "fail":
        print("toy.py: error: asked to fail", file=sys.stderr)
        sys.exit(3)
    if mode == "vandal":
        open(sys.argv[2], "a").write("vandalised\\n")
    if mode == "workers":
        print("workers", sys.argv[2])
    if mode == "start":
        import multiprocessing
        print("start", multiprocessing.get_start_method(), __name__)
    total = sum(r["n"] for r in rows) + (1 if mode == "wrong" else 0)
    print(toy_helper.banner(), "total", total, "(0.%d s)" % len(mode))
    with open(os.path.join(HERE, "toy_output.jsonl"), "a") as fh:
        for r in rows:
            fh.write(json.dumps({"n": r["n"] * 2, "seconds": len(mode)}) + "\\n")
    """
)

STORED = {
    "toy_input.jsonl": "".join(json.dumps({"n": n}) + "\n" for n in (1, 2, 3)),
    "toy_output.jsonl": "".join(json.dumps({"n": n, "seconds": 99}) + "\n" for n in (2, 4, 6)),
    "toy.log": "toy total 6 (5.1 s)\n",
    "sixlines_ledger.log": "marker\n",
}


def toy(mode: str = "ok", **changes: object) -> Computation:
    fields: dict = {
        "id": "n99/toy",
        "note": 99,
        "group": "ledger",
        "script": "toy.py",
        "args": (mode,),
        "seconds": 1,
        "checks": (
            Check(STDOUT, "toy.log", "text"),
            Check("toy_output.jsonl", "toy_output.jsonl", "records"),
        ),
        "expect": (r"total \d+",),
    }
    fields.update(changes)
    return Computation(**fields)


class Workspace(unittest.TestCase):
    def setUp(self) -> None:
        scratch = tempfile.TemporaryDirectory()
        self.addCleanup(scratch.cleanup)
        root = Path(scratch.name)
        scripts, data = root / "experiments", root / "data"
        scripts.mkdir()
        data.mkdir()
        (scripts / "toy.py").write_text(TOY_SCRIPT)
        (scripts / "toy_helper.py").write_text("def banner():\n    return 'toy'\n")
        for name, content in STORED.items():
            (data / name).write_text(content)
        self.places = Places(scripts=scripts, data=data, work=root / "work", python=sys.executable)
        self.stored_before = {p.name: p.read_bytes() for p in data.iterdir()}

    def tearDown(self) -> None:
        after = {p.name: p.read_bytes() for p in self.places.data.iterdir()}
        self.assertEqual(after, self.stored_before, "a stored file was changed")

    def run_toy(self, computation: Computation, keep: bool = False) -> dict:
        execution = execute(computation, self.places, 1)
        self.assertEqual(execution.exit_code, 0, execution.stderr.read_text())
        return finish(execution, evaluate(computation, execution, self.places), keep)

    def test_work_directory_links_scripts_and_copies_inputs_but_not_outputs(self) -> None:
        work = prepare(toy(), self.places)
        self.assertTrue((work / "toy.py").is_symlink())
        self.assertTrue((work / "toy_helper.py").is_symlink())
        self.assertFalse((work / "toy_input.jsonl").is_symlink())
        self.assertEqual((work / "toy_input.jsonl").read_text(), STORED["toy_input.jsonl"])
        # the computation writes toy_output.jsonl, so it starts absent even though a stored copy exists
        self.assertFalse((work / "toy_output.jsonl").exists())

    def test_reproduced(self) -> None:
        result = self.run_toy(toy())
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["verdict"], "reproduced")
        self.assertTrue(result["agrees"])
        self.assertEqual(result["value"], result["predicted"])
        self.assertIn("toy_output.jsonl", result["digests"])
        self.assertNotIn("work_dir", result)
        self.assertEqual(list(self.places.work.iterdir()), [])

    def test_keep_work(self) -> None:
        result = self.run_toy(toy(), keep=True)
        self.assertTrue(Path(result["work_dir"]).is_dir())

    def test_differs_and_keeps_the_work_directory(self) -> None:
        result = self.run_toy(toy("wrong"))
        self.assertEqual(result["status"], "mismatch")
        self.assertEqual(result["verdict"], "differs")
        self.assertEqual(result["value"]["toy.log"], "differs")
        self.assertEqual(result["value"]["toy_output.jsonl"], "identical")
        self.assertEqual(result["differences"]["toy.log"]["rerun"], "toy total 7 (#s)")
        self.assertTrue(Path(result["work_dir"]).is_dir())

    def test_a_control_ends_as_a_mismatch_and_leaves_no_work_directory(self) -> None:
        execution = execute(toy(), self.places, 1)
        result = evaluate(toy(), execution, self.places, control={"toy.log": "differs"})
        self.assertEqual((result["status"], result["agrees"], result["verdict"]), ("mismatch", False, "differs"))
        self.assertEqual(result["value"]["toy.log"], "identical")
        self.assertEqual(result["predicted"]["toy.log"], "differs")
        finish(execution, result, keep=False, discard=True)
        self.assertNotIn("work_dir", result)
        self.assertFalse(execution.work_dir.exists())

    def test_checked_by_output_lines(self) -> None:
        result = self.run_toy(toy(checks=(), counts=((r"^toy total", 1), (r"never", 0))))
        self.assertEqual(result["verdict"], "checked")
        self.assertEqual(result["status"], "ok")
        result = self.run_toy(toy(checks=(), expect=(r"total 99",)))
        self.assertEqual(result["verdict"], "differs")
        self.assertEqual(result["value"]["expect /total 99/"], "absent")

    def test_ran_when_nothing_is_compared(self) -> None:
        result = self.run_toy(toy(checks=(), expect=(), caveat="output not kept"))
        self.assertEqual((result["verdict"], result["agrees"], result["status"]), ("ran", None, "ok"))
        self.assertEqual(result["caveat"], "output not kept")

    def test_procs_reach_the_command(self) -> None:
        computation = toy("workers", args=("workers", "{procs}"), procs=4, checks=(), expect=(r"^workers 3$",))
        execution = execute(computation, self.places, 3)
        self.assertEqual(evaluate(computation, execution, self.places)["verdict"], "checked")

    def test_scripts_that_need_fork_get_it(self) -> None:
        computation = toy("start", args=("start",), checks=(), expect=(r"^start fork __main__$",))
        with mock.patch("k3tasks.registry.FORK_SCRIPTS", frozenset({"toy.py"})):
            self.assertTrue(computation.needs_fork)
            result = self.run_toy(computation)
        self.assertEqual(result["verdict"], "checked")

    def test_failure_is_reported(self) -> None:
        execution = execute(toy("fail"), self.places, 1)
        self.assertEqual(execution.exit_code, 3)
        self.assertIn("asked to fail", execution.stderr.read_text())

    def test_a_changed_stored_file_is_an_error(self) -> None:
        target = self.places.data / "toy.log"
        computation = toy("vandal", args=("vandal", str(target)))
        with self.assertRaisesRegex(WorkspaceError, "toy.log changed"):
            execute(computation, self.places, 1)
        target.write_text(STORED["toy.log"])  # undo, for tearDown
        os.utime(target)

    def test_missing_data_is_explained(self) -> None:
        (self.places.data / "sixlines_ledger.log").unlink()
        self.stored_before.pop("sixlines_ledger.log")
        with self.assertRaisesRegex(WorkspaceError, "not installed"):
            prepare(toy(), self.places)

    def test_missing_stored_reference(self) -> None:
        computation = toy(checks=(Check(STDOUT, "absent.log", "text"),))
        with self.assertRaisesRegex(WorkspaceError, "absent.log"):
            prepare(computation, self.places)
