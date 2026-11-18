"""Comparing a rerun's files with the stored ones: the modes, and what normalisation blanks out."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from k3tasks.compare import DIFFERS, IDENTICAL, MISSING, compare, normalise_text


class Normalisation(unittest.TestCase):
    def test_run_times_are_blanked(self) -> None:
        self.assertEqual(
            normalise_text("p = 7: 3 classes (8.1s)\np = 11: done (12 s)\nPOP; m = 5 (11 s)\ntook 3 minutes"),
            ["p = 7: 3 classes (#s)", "p = 11: done (#s)", "POP; m = 5 (#s)", "took #s"],
        )

    def test_numbers_that_are_not_times_stay(self) -> None:
        self.assertEqual(normalise_text("e = [6, -369, -20172] signs [1]"), ["e = [6, -369, -20172] signs [1]"])

    def test_shell_and_gnu_time_lines_are_dropped(self) -> None:
        text = "result 1\nreal\t0m20.574s\nuser\t0m21.068s\nsys\t0m0.055s\n\tCommand being timed: \"python x\"\n"
        text += "\tMaximum resident set size (kbytes): 4000\n"
        self.assertEqual(normalise_text(text), ["result 1"])

    def test_pari_stack_notice_is_dropped(self) -> None:
        text = "PARI stack size set to 400000000 bytes, maximum size set to 2000003072\nselftest passed\n"
        self.assertEqual(normalise_text(text), ["selftest passed"])

    def test_path_prefixes_are_removed(self) -> None:
        self.assertEqual(
            normalise_text("wrote /work/dir/out.json", ("/work/dir",)), ["wrote out.json"]
        )

    def test_run_times_inside_printed_records_are_blanked(self) -> None:
        self.assertEqual(
            normalise_text('  exception: {"p": 11, "seconds": 61.3, "setup_s": 0.4, "agrees": false}'),
            ['  exception: {"p": 11, "seconds": #, "setup_s": #, "agrees": false}'],
        )

    def test_dicts_of_counts_are_compared_in_sorted_order(self) -> None:
        first = normalise_text("status: {'ok': 296, 'not ordinary': 86, 'reducible': 272}")
        self.assertEqual(first, normalise_text("status: {'reducible': 272, 'ok': 296, 'not ordinary': 86}"))
        self.assertNotEqual(first, normalise_text("status: {'reducible': 271, 'ok': 296, 'not ordinary': 86}"))
        # other dicts keep their order
        self.assertEqual(normalise_text("{'b': [1], 'a': 2}"), ["{'b': [1], 'a': 2}"])

    def test_trailing_blank_lines_and_spaces(self) -> None:
        self.assertEqual(normalise_text("a  \nb\n\n\n"), ["a", "b"])


class Modes(unittest.TestCase):
    def setUp(self) -> None:
        scratch = tempfile.TemporaryDirectory()
        self.addCleanup(scratch.cleanup)
        self.dir = Path(scratch.name)

    def files(self, produced: str, stored: str) -> tuple[Path, Path]:
        a, b = self.dir / "produced", self.dir / "stored"
        a.write_text(produced)
        b.write_text(stored)
        return a, b

    def test_exact(self) -> None:
        self.assertEqual(compare(*self.files("x\n", "x\n"), "exact").verdict, IDENTICAL)
        outcome = compare(*self.files("x\ny\n", "x\nz\n"), "exact")
        self.assertEqual(outcome.verdict, DIFFERS)
        self.assertEqual(outcome.note["line"], 2)
        self.assertEqual((outcome.note["rerun"], outcome.note["stored"]), ("y", "z"))

    def test_exact_sees_run_times(self) -> None:
        self.assertEqual(compare(*self.files("done (3 s)\n", "done (4 s)\n"), "exact").verdict, DIFFERS)

    def test_text_ignores_run_times(self) -> None:
        self.assertEqual(compare(*self.files("done (3 s)\n", "done (4 s)\n"), "text").verdict, IDENTICAL)
        self.assertEqual(compare(*self.files("done 3\n", "done 4\n"), "text").verdict, DIFFERS)

    def test_text_contained(self) -> None:
        stored = "first command\nA\n\nsecond command\nB 1 (2 s)\nC\nthird\n"
        self.assertEqual(compare(*self.files("second command\nB 1 (9 s)\nC\n", stored), "text-contained").verdict, IDENTICAL)
        outcome = compare(*self.files("second command\nB 2\nC\n", stored), "text-contained")
        self.assertEqual(outcome.verdict, DIFFERS)
        self.assertEqual(outcome.note["line"], 2)
        self.assertEqual(outcome.note["stored"], "B 1 (#s)")

    def test_records_in_any_order_without_times(self) -> None:
        produced = json.dumps({"a": 1, "seconds": 3}) + "\n" + json.dumps({"a": 2, "seconds": 1}) + "\n"
        stored = json.dumps({"a": 2, "seconds": 9}) + "\n" + json.dumps({"seconds": 8, "a": 1}) + "\n"
        self.assertEqual(compare(*self.files(produced, stored), "records", ignore_keys=("seconds",)).verdict, IDENTICAL)
        self.assertEqual(compare(*self.files(produced, stored), "records").verdict, DIFFERS)

    def test_records_counts_multiplicity(self) -> None:
        line = json.dumps({"a": 1}) + "\n"
        outcome = compare(*self.files(line, line * 2), "records")
        self.assertEqual(outcome.verdict, DIFFERS)
        self.assertEqual((outcome.note["not_stored"], outcome.note["not_recomputed"]), (0, 1))

    def test_records_subset(self) -> None:
        one, two = json.dumps({"a": 1}) + "\n", json.dumps({"a": 2}) + "\n"
        self.assertEqual(compare(*self.files(one, one + two), "records-subset").verdict, IDENTICAL)
        outcome = compare(*self.files(one + json.dumps({"a": 3}) + "\n", one + two), "records-subset")
        self.assertEqual(outcome.verdict, DIFFERS)
        self.assertEqual(outcome.note["not_stored"], 1)
        self.assertEqual(compare(*self.files("", one), "records-subset").verdict, DIFFERS)

    def test_json(self) -> None:
        produced = json.dumps({"x": [1, 2], "seconds": 1})
        self.assertEqual(compare(*self.files(produced, json.dumps({"x": [1, 2], "seconds": 5})), "json",
                                 ignore_keys=("seconds",)).verdict, IDENTICAL)
        outcome = compare(*self.files(produced, json.dumps({"x": [1, 3]})), "json", ignore_keys=("seconds",))
        self.assertEqual(outcome.verdict, DIFFERS)
        self.assertEqual(outcome.note["keys"], ["x"])

    def test_missing_files(self) -> None:
        stored = self.dir / "stored"
        stored.write_text("x")
        self.assertEqual(compare(self.dir / "absent", stored, "exact").verdict, MISSING)
        self.assertEqual(compare(stored, self.dir / "absent", "exact").verdict, MISSING)

    def test_unknown_mode(self) -> None:
        with self.assertRaises(ValueError):
            compare(*self.files("", ""), "fuzzy")
