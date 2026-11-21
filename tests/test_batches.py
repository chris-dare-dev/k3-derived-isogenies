"""The batch files: they load, expand, and give every task limits that fit its expected time and memory."""

from __future__ import annotations

import contextlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from k3tasks import __main__ as command_line
from k3tasks import registry
from k3tasks.tasks import ComputationTask, SyntheticTask, expand_batch, plan_batch, run_task

REPOSITORY = Path(__file__).resolve().parent.parent
BATCHES = REPOSITORY / "batches"

sys.path.insert(0, str(REPOSITORY / "tools"))

import claims  # noqa: E402

# a task's time limit is at least this multiple of its expected time (and a minute more)
TIME_MARGIN = 1.5


class Batches(unittest.TestCase):
    def test_there_are_batches(self) -> None:
        names = {path.stem for path in BATCHES.glob("*.toml")}
        self.assertTrue({"smoke", "heavy", "lattices", "fourfold", "cm-method"} <= names)

    def test_limits_fit_the_expected_time_and_memory(self) -> None:
        for path in sorted(BATCHES.glob("*.toml")):
            plan = plan_batch(path)
            self.assertTrue(plan.tasks, path.name)
            for payload in plan.tasks:
                if payload["kind"] == "synthetic":
                    continue  # these are stopped at their limits on purpose
                limits = plan.limits(payload)
                with self.subTest(batch=path.stem, task=payload["id"]):
                    self.assertGreaterEqual(limits.time_seconds, TIME_MARGIN * payload["cost"] + 60)
                    self.assertGreaterEqual(limits.memory_mb, payload["expected_memory_mb"])
                    self.assertLessEqual(limits.memory_mb, plan.config.memory_budget_mb)

    def test_processes_stay_within_the_machine(self) -> None:
        # at most 16 processes at a time: workers times the largest procs of a task
        for path in sorted(BATCHES.glob("*.toml")):
            plan = plan_batch(path)
            with self.subTest(batch=path.stem):
                self.assertLessEqual(plan.config.workers * max(task["procs"] for task in plan.tasks), 16)

    def test_heavy_computations_only_in_the_heavy_batch(self) -> None:
        heavy = {c.id for c in registry.select(tags=["heavy"])}
        for path in sorted(BATCHES.glob("*.toml")):
            ids = {task["id"] for task in plan_batch(path).tasks}
            if path.stem == "heavy":
                self.assertEqual(heavy - ids, set())
            else:
                self.assertEqual(ids & heavy, set(), path.stem)

    def test_every_computation_is_in_a_batch(self) -> None:
        batched = set()
        for path in BATCHES.glob("*.toml"):
            batched |= {task["id"] for task in plan_batch(path).tasks}
        self.assertEqual(sorted({c.id for c in registry.all_computations()} - batched), [])

    def test_smoke_is_quick(self) -> None:
        tasks = plan_batch(BATCHES / "smoke.toml").tasks
        ids = {task["id"] for task in tasks}
        self.assertTrue({"n27/ledger", "n20/selftest", "n14/prym-check-17"} <= ids)
        self.assertLess(sum(task["cost"] for task in tasks), 60)

    def test_the_no_data_batch_is_what_needs_no_stored_result_and_is_quick(self) -> None:
        tasks = plan_batch(BATCHES / "no-data.toml").tasks
        self.assertEqual({task["id"] for task in tasks}, {c.id for c in registry.select(tags=["no-data"])})
        self.assertTrue(all(task["compared"] == [] and task["seeds"] is None for task in tasks))
        self.assertLess(sum(task["cost"] for task in tasks), 60)

    def test_the_command_line_says_in_one_line_whether_the_stored_results_are_installed(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            output, errors = io.StringIO(), io.StringIO()
            with (
                mock.patch.dict(os.environ, {"K3_DATA_DIR": scratch}),
                contextlib.redirect_stdout(output),
                contextlib.redirect_stderr(errors),
            ):
                self.assertEqual(command_line.main(["data"]), 1)
                (Path(scratch) / "sixlines_ledger.log").write_text("installed\n")
                self.assertEqual(command_line.main(["data"]), 0)
        (line,) = errors.getvalue().splitlines()
        self.assertIn("the stored results are not installed", line)
        self.assertIn("`make data`", line)
        self.assertIn("`make no-data` runs the 6 computations that need none", line)
        self.assertIn("stored results: installed in", output.getvalue())


class Expansion(unittest.TestCase):
    def test_selectors_and_settings(self) -> None:
        tasks = expand_batch(
            {
                "tasks": [
                    {"kind": "computation", "ids": ["n02/enumerate-p37"], "procs": 4, "keep_work": True},
                    {"kind": "computation", "notes": [2], "exclude": ["n02/enumerate-*"]},
                ]
            }
        )
        first = tasks[0]
        self.assertEqual(
            (first["id"], first["procs"], first["set"], first["keep_work"]), ("n02/enumerate-p37", 4, 0, True)
        )
        self.assertEqual(first["cost"], 2 * registry.get("n02/enumerate-p37").seconds)
        self.assertIn("--procs 4", first["command"])
        self.assertEqual({t["set"] for t in tasks[1:]}, {1})
        self.assertNotIn("n02/enumerate-p7", [t["id"] for t in tasks[1:]])

    def test_a_computation_is_taken_once(self) -> None:
        tasks = expand_batch(
            {
                "tasks": [
                    {"kind": "computation", "ids": ["n27/ledger"]},
                    {"kind": "computation", "groups": ["ledger"], "procs": 2},
                ]
            }
        )
        self.assertEqual([(t["id"], t["set"]) for t in tasks], [("n27/ledger", 0)])

    def test_procs_never_exceed_the_recorded_number(self) -> None:
        tasks = expand_batch({"tasks": [{"kind": "computation", "ids": ["n27/ledger"], "procs": 8}]})
        self.assertEqual(tasks[0]["procs"], 1)

    def test_a_control_is_the_same_computation_with_a_wrong_prediction(self) -> None:
        plain, control = expand_batch(
            {
                "tasks": [
                    {"kind": "computation", "ids": ["n27/ledger"]},
                    {
                        "kind": "computation",
                        "ids": ["n27/ledger"],
                        "control": {"predicted": {"sixlines_ledger.log": "differs"}},
                    },
                ]
            }
        )
        self.assertEqual(plain["predicted"], {"sixlines_ledger.log": "identical"})
        self.assertEqual(plain["compared"], [{"file": "sixlines_ledger.log", "mode": "exact"}])
        self.assertNotIn("control", plain)
        self.assertEqual(control["id"], "n27/ledger!control")
        self.assertEqual((control["group"], control["control"]), ("negative control", True))
        self.assertEqual(control["predicted"], {"sixlines_ledger.log": "differs"})
        self.assertEqual(control["command"], plain["command"])

    def test_the_negative_control_batch_holds_controls_only(self) -> None:
        tasks = plan_batch(BATCHES / "negative-control.toml").tasks
        self.assertEqual(len(tasks), 3)
        for task in tasks:
            self.assertTrue(task["control"] and task["id"].endswith("!control"), task["id"])
            self.assertNotEqual(
                task["predicted"], registry.get(task["id"].removesuffix("!control")).predicted()
            )

    def test_controls_and_synthetic_tasks_stay_in_their_own_batches(self) -> None:
        for path in sorted(BATCHES.glob("*.toml")):
            for task in plan_batch(path).tasks:
                if task.get("control"):
                    self.assertEqual(path.stem, "negative-control", task["id"])
                if task["kind"] == "synthetic":
                    self.assertEqual(path.stem, "limits-demo", task["id"])

    def test_synthetic_tasks(self) -> None:
        tasks = expand_batch({"tasks": [{"kind": "synthetic", "items": [
            {"name": "sleeps", "sleep_seconds": 5}, {"name": "holds", "allocate_mb": 3},
        ]}]})  # fmt: skip
        self.assertEqual([t["id"] for t in tasks], ["synthetic/sleeps", "synthetic/holds"])
        self.assertEqual({t["group"] for t in tasks}, {"synthetic"})
        result = run_task(tasks[1])
        self.assertEqual((result["status"], result["value"], result["agrees"]), ("ok", 3, None))

    def test_no_way_to_the_tasks_of_a_batch_passes_an_inconsistent_registry(self) -> None:
        smoke = BATCHES / "smoke.toml"
        with mock.patch.object(registry, "validate", return_value=["n00/x: unknown group g"]):
            with self.assertRaisesRegex(ValueError, "the registry is inconsistent: n00/x: unknown group g"):
                expand_batch({"tasks": [{"kind": "computation", "ids": ["n27/ledger"]}]})
            with self.assertRaisesRegex(ValueError, "the registry is inconsistent"):
                plan_batch(smoke)
            with self.assertRaisesRegex(ValueError, "the registry is inconsistent"):
                claims.batch_tasks()
            errors = io.StringIO()
            with contextlib.redirect_stderr(errors):
                self.assertEqual(command_line.main(["list", "--batch", str(smoke)]), 2)
            self.assertIn("the registry is inconsistent", errors.getvalue())

    def test_a_plan_gives_each_task_the_limits_of_its_table(self) -> None:
        plan = plan_batch(BATCHES / "limits-demo.toml")
        limits = {task["id"]: plan.limits(task).as_dict() for task in plan.tasks}
        self.assertEqual(
            limits,
            {
                "n27/ledger": {"time_seconds": 120, "memory_mb": 1000},
                "synthetic/sleeps-past-its-time-limit": {"time_seconds": 3, "memory_mb": 300},
                "synthetic/allocates-past-its-memory-limit": {"time_seconds": 120, "memory_mb": 300},
            },
        )

    def test_a_task_has_the_keys_of_its_type_and_no_other(self) -> None:
        # the keys are in the committed records: a control has two more, and nothing else varies
        kinds = set()
        for path in sorted(BATCHES.glob("*.toml")):
            for task in plan_batch(path).tasks:
                if task["kind"] == "synthetic":
                    expected = SyntheticTask.__required_keys__
                elif task.get("control"):
                    expected = ComputationTask.__required_keys__ | ComputationTask.__optional_keys__
                else:
                    expected = ComputationTask.__required_keys__
                self.assertEqual(set(task), expected, task["id"])
                kinds.add((task["kind"], bool(task.get("control"))))
        self.assertEqual(kinds, {("computation", False), ("computation", True), ("synthetic", False)})
        self.assertEqual(ComputationTask.__optional_keys__, {"control", "wrong_on_purpose"})
        self.assertEqual(SyntheticTask.__optional_keys__, frozenset())

    def test_refusals(self) -> None:
        for table in (
            {"kind": "synthetic", "ids": ["n27/ledger"]},
            {"kind": "computation"},
            {"kind": "computation", "ids": ["n27/ledger"], "colour": "red"},
            {"kind": "computation", "ids": ["n27/ledger"], "procs": 0},
            {"kind": "computation", "ids": ["n27/ledger"], "keep_work": "yes"},
            {"kind": "computation", "ids": ["n99/*"]},
            {"kind": "computation", "notes": [1]},
            {"kind": "computation", "notes": ["2"]},
            {"kind": "computation", "ids": ["n27/ledger"], "control": {"predicted": 3}},
            {
                "kind": "computation",
                "ids": ["n27/ledger"],
                "control": {"predicted": {"another.log": "differs"}},
            },
            {"kind": "computation", "ids": ["n27/ledger"], "control": {"value": 1}},
            {"kind": "census", "ids": ["n27/ledger"]},
            {"kind": "synthetic", "items": [{"sleep_seconds": 1}]},
            {"kind": "synthetic", "items": [{"name": "a", "colour": "red"}]},
            {"kind": "synthetic", "items": [{"name": "a"}, {"name": "a"}]},
        ):
            with self.subTest(table=table), self.assertRaises(ValueError):
                expand_batch({"tasks": [table]})
