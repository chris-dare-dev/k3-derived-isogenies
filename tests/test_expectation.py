"""The words of a result: the keys of what is predicted and found, the verdict, and what differs.

The keys and the words are in the committed records, so they are spelled out here once.
"""

from __future__ import annotations

import unittest
from typing import Any

from k3tasks import expectation, registry
from k3tasks.expectation import compared_with, differences, differing, judge, observed, predicted, verdict

OUTPUT = ["selftest passed", "  pair 1 2 equal", "  pair 1 3 equal", "total 6"]


def record(status: str = "ok", **fields: Any) -> dict[str, Any]:
    return {"status": status, "value": None, "predicted": None, "detail": {}} | fields


class Keys(unittest.TestCase):
    def test_the_keys_and_the_words_are_those_of_the_records(self) -> None:
        self.assertEqual(
            predicted(["a.log", "b.jsonl"], ["selftest passed"], [(r"pair .* equal$", 2)]),
            {
                "a.log": "identical",
                "b.jsonl": "identical",
                "expect /selftest passed/": "found",
                "lines /pair .* equal$/": 2,
            },
        )
        self.assertEqual(
            observed(OUTPUT, ["selftest passed", "total 99"], [(r"pair .* equal$", 3)]),
            {"expect /selftest passed/": "found", "expect /total 99/": "absent", "lines /pair .* equal$/": 2},
        )
        self.assertEqual(
            expectation.ORDER,
            ["reproduced", "differs", "checked", "ran", "timeout", "memory", "error", "crashed",
             "unsupported", "not run"],
        )  # fmt: skip
        self.assertEqual(list(expectation.ENDED_WELL), ["reproduced", "checked", "ran"])

    def test_what_is_found_has_the_keys_of_what_is_predicted(self) -> None:
        for computation in registry.all_computations():
            found = observed([], computation.expect, computation.counts)
            files = {check.label for check in computation.checks}
            self.assertEqual(set(computation.predicted()), files | set(found), computation.id)

    def test_the_registry_predicts_every_comparison(self) -> None:
        self.assertEqual(registry.get("n27/ledger").predicted(), {"sixlines_ledger.log": "identical"})
        self.assertEqual(registry.get("n14/prym-check-17").predicted()[r"lines /^\s*pair .*\bequal$/"], 18)
        self.assertEqual(registry.get("n02/selftest").predicted(), {"expect /selftest passed/": "found"})
        self.assertEqual(registry.get("n13/groupsize-41").predicted(), {})

    def test_a_stored_file_is_told_from_a_line_of_the_output(self) -> None:
        value = {"a.log": "identical", "expect /x/": "found", "lines /y/": 3}
        self.assertEqual(compared_with(record(value=value)), ["a.log"])
        self.assertEqual(compared_with(None), [])
        self.assertEqual(compared_with(record()), [])


class Judging(unittest.TestCase):
    def test_everything_identical_is_reproduced(self) -> None:
        found = judge({"a.log": "identical"}, {"a.log": "identical"}, files_compared=True)
        self.assertEqual((found.status, found.agrees, found.verdict), ("ok", True, "reproduced"))

    def test_any_difference_is_a_mismatch(self) -> None:
        found = judge({"a.log": "differs"}, {"a.log": "identical"}, files_compared=True)
        self.assertEqual((found.status, found.agrees, found.verdict), ("mismatch", False, "differs"))
        found = judge({"lines /y/": 2}, {"lines /y/": 3}, files_compared=False)
        self.assertEqual((found.status, found.agrees, found.verdict), ("mismatch", False, "differs"))

    def test_lines_alone_are_checked_and_nothing_at_all_only_ran(self) -> None:
        found = judge({"expect /x/": "found"}, {"expect /x/": "found"}, files_compared=False)
        self.assertEqual((found.status, found.agrees, found.verdict), ("ok", True, "checked"))
        found = judge({}, {}, files_compared=False)
        self.assertEqual((found.status, found.agrees, found.verdict), ("ok", None, "ran"))


class Records(unittest.TestCase):
    def test_the_verdict_of_a_record(self) -> None:
        self.assertEqual(verdict(None), "not run")
        self.assertEqual(verdict(record(detail={"verdict": "checked"})), "checked")
        self.assertEqual(verdict(record()), "ran")
        self.assertEqual(verdict(record("mismatch")), "differs")
        for status in ("timeout", "memory", "error", "crashed", "unsupported"):
            self.assertEqual(verdict(record(status, detail={"verdict": "reproduced"})), status)

    def test_what_differs(self) -> None:
        found = record(
            "mismatch",
            value={"a.log": "differs", "b.log": "identical", "expect /x/": "absent", "lines /y/": 2},
            predicted={"a.log": "identical", "b.log": "identical", "expect /x/": "found", "lines /y/": 3},
            detail={
                "differences": {
                    "a.log": {"mode": "text", "line": 4, "produced": "p", "stored": "s"},
                    "b.log": {"mode": "text-contained", "lines_not_stored": 2},
                }
            },
        )
        self.assertEqual(differing(found), {"a.log": "differs", "expect /x/": "absent", "lines /y/": 2})
        self.assertEqual(
            differences(found),
            [
                "a.log (text): line 4, produced p, stored s",
                "expect /x/: absent, expected found",
                "lines /y/: 2, expected 3",
            ],
        )

    def test_a_control_that_was_reproduced_differs_from_its_wrong_prediction(self) -> None:
        found = record("mismatch", value={"a.log": "identical"}, predicted={"a.log": "differs"})
        self.assertEqual(differences(found), ["a.log: identical, predicted differs"])
        self.assertEqual(differing(record()), {})


if __name__ == "__main__":
    unittest.main()
