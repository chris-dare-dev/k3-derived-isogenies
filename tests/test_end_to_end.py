"""The scheduler with the repository's runner, on a stand-in for the ledger script and its stored log.

K3_SCRIPTS_DIR and K3_DATA_DIR point the runner at a scratch copy, so this needs neither numpy, cypari2
nor the data archive.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from k3tasks.report import render

REPOSITORY = Path(__file__).resolve().parent.parent

BATCH = """
name = "e2e"
[scheduler]
workers = 2
memory_budget_mb = 1000
[limits]
time_seconds = 60
memory_mb = 400
[[tasks]]
kind = "computation"
ids = ["n27/ledger"]
"""

LEDGER = """import os
HERE = os.path.dirname(os.path.abspath(__file__))
print(open(os.path.join(HERE, "sixlines_ledger.log")).read(){suffix}, end="")
"""


class EndToEnd(unittest.TestCase):
    def setUp(self) -> None:
        scratch = tempfile.TemporaryDirectory()
        self.addCleanup(scratch.cleanup)
        self.root = Path(scratch.name)
        self.scripts = self.root / "experiments"
        self.scripts.mkdir()
        (self.scripts / "sixlines_ledger.log").write_text("170 comparisons, 162 agree\n")
        self.batch = self.root / "e2e.toml"
        self.batch.write_text(BATCH)
        self.run_dir = self.root / "run"
        self.env = os.environ | {
            "K3_SCRIPTS_DIR": str(self.scripts),
            "K3_WORK_DIR": str(self.root / "work"),
            "K3_PYTHON": sys.executable,
        }
        self.env.pop("K3_DATA_DIR", None)

    def taskbatch(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, "-m", "taskbatch", *arguments, str(self.batch), "--out", str(self.run_dir),
             "--memory-guard", "poll"] if arguments[0] == "run" else
            [sys.executable, "-m", "taskbatch", *arguments, str(self.batch), "--out", str(self.run_dir)],
            capture_output=True, text=True, cwd=REPOSITORY, env=self.env, timeout=120,
        )  # fmt: skip

    def records(self) -> list[dict]:
        return [json.loads(line) for line in (self.run_dir / "results.jsonl").read_text().splitlines()]

    def test_plan(self) -> None:
        done = self.taskbatch("plan")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("n27/ledger", done.stdout)
        self.assertIn("1 tasks with limits 60 s, 400 MB", done.stdout)

    def test_reproduced(self) -> None:
        (self.scripts / "sixlines_ledger.py").write_text(LEDGER.format(suffix=""))
        done = self.taskbatch("run")
        self.assertEqual(done.returncode, 0, done.stderr)
        (record,) = self.records()
        self.assertEqual(record["status"], "ok")
        self.assertEqual(record["detail"]["verdict"], "reproduced")
        self.assertEqual(record["value"], {"sixlines_ledger.log": "identical"})
        self.assertIn("stdout", record["detail"]["digests"])
        self.assertIn("data", record["runner"])
        self.assertEqual(list((self.root / "work").iterdir()), [])
        report = render(self.run_dir)
        self.assertIn("## Reproduction", report)
        self.assertIn("1 reproduced.", report)
        # a second run computes nothing
        again = self.taskbatch("run")
        self.assertEqual(again.returncode, 0)
        self.assertEqual(len(self.records()), 1)

    def test_differs(self) -> None:
        (self.scripts / "sixlines_ledger.py").write_text(LEDGER.format(suffix=".replace('162', '161')"))
        done = self.taskbatch("run")
        self.assertEqual(done.returncode, 1)
        (record,) = self.records()
        self.assertEqual(record["status"], "mismatch")
        self.assertEqual(record["detail"]["differences"]["sixlines_ledger.log"]["rerun"], "170 comparisons, 161 agree")
        report = render(self.run_dir)
        self.assertIn("1 differs.", report)
        self.assertIn("## What differs", report)
        self.assertTrue(Path(record["detail"]["work_dir"]).is_dir())

    def test_a_failing_script_is_an_error(self) -> None:
        (self.scripts / "sixlines_ledger.py").write_text("raise SystemExit('ledger: error: no stored files')\n")
        done = self.taskbatch("run")
        self.assertEqual(done.returncode, 1)
        (record,) = self.records()
        self.assertEqual(record["status"], "error")
        self.assertIn("exited with status 1", record["error"])

    def test_missing_data_is_an_error(self) -> None:
        (self.scripts / "sixlines_ledger.py").write_text(LEDGER.format(suffix=""))
        (self.scripts / "sixlines_ledger.log").unlink()
        self.taskbatch("run")
        (record,) = self.records()
        self.assertEqual(record["status"], "error")
        self.assertIn("not installed", record["error"])
