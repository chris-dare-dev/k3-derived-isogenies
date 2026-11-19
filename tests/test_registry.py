"""The registry of computations against the repository: scripts, stored files, notes, batches."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

from k3tasks import registry
from k3tasks.registry import STDERR, STDOUT

REPOSITORY = Path(__file__).resolve().parent.parent
EXPERIMENTS = REPOSITORY / "experiments"


def archive_manifest() -> set[str]:
    """The file names of the data archive, from the manifest table of DATA.md."""
    text = (REPOSITORY / "DATA.md").read_text(encoding="utf-8")
    manifest = text.split("## Manifest", 1)[1]
    return set(re.findall(r"^\| `([^`]+)` \|", manifest, flags=re.MULTILINE))


class Registry(unittest.TestCase):
    def test_consistent(self) -> None:
        self.assertEqual(registry.validate(), [])

    def test_every_script_exists(self) -> None:
        for c in registry.all_computations():
            with self.subTest(c.id):
                self.assertTrue((EXPERIMENTS / c.script).is_file(), c.script)

    def test_every_stored_file_is_in_the_archive(self) -> None:
        manifest = archive_manifest()
        self.assertEqual(len(manifest), 171)
        for c in registry.all_computations():
            for check in c.checks:
                with self.subTest(c.id, stored=check.stored):
                    self.assertIn(check.stored, manifest)

    def test_every_log_of_the_archive_is_compared_or_accounted_for(self) -> None:
        compared = {check.stored for c in registry.all_computations() for check in c.checks}
        # logs that record a timing or an enumeration rather than a command's output
        not_outputs = {
            "sixlines_enumerate_29_37.log",
            "sixlines_enumerate_41.log",
            "sixlines_enumerate_43.log",
            "sixlines_delta41_primary.log",
            "sixlines_delta41_more.log",
            "sixlines_ledger.log",
        }
        logs = {name for name in archive_manifest() if name.endswith(".log")}
        self.assertEqual(sorted(logs - compared - not_outputs), [])

    def test_lines_left_out_of_a_comparison_are_stated(self) -> None:
        left_out = [(c.id, check.unstored_tail) for c in registry.all_computations() for check in c.checks
                    if check.unstored_tail]  # fmt: skip
        self.assertEqual(left_out, [("n26/list", 1)])
        self.assertIn("the last line of the rerun is not", registry.get("n26/list").caveat)

    def test_the_key_added_after_the_early_coset_runs_is_compared_where_it_is_stored(self) -> None:
        # of the stored records of note 29's coset runs, only two have the key identified_by
        compared = {
            c.id for c in registry.select(ids=["n29/cosets-*"])
            if any("identified_by" not in check.ignore_keys for check in c.checks if check.mode == "records-subset")
        }  # fmt: skip
        self.assertEqual(compared, {"n29/cosets-43-2-coord", "n29/cosets-43-2-mix4500"})
        for c in registry.select(ids=["n29/cosets-*"]):
            if c.id not in compared:
                self.assertIn("identified_by", c.caveat, c.id)

    def test_note_matches_id(self) -> None:
        for c in registry.all_computations():
            self.assertTrue(c.id.startswith(f"n{c.note:02d}/"), c.id)
            self.assertTrue(any((REPOSITORY / "docs").glob(f"{c.note:02d}-*.md")), c.id)

    def test_scripts_named_in_a_reproduce_section_have_a_computation(self) -> None:
        named = set()
        for note in sorted((REPOSITORY / "docs").glob("[0-9][0-9]-*.md")):
            text = note.read_text(encoding="utf-8")
            if "Reproduce" in text:
                section = text.split("Reproduce", 1)[1]
                named |= set(re.findall(r"(?:python|\$PY) -I (\w+\.py)", section))
        scripts = {c.script for c in registry.all_computations()}
        self.assertEqual(sorted(named - scripts), [])

    def test_streams_and_written_files(self) -> None:
        c = registry.get("n02/cm-list")
        self.assertEqual({check.produced for check in c.checks}, {STDOUT, "sixlines_cm.jsonl"})
        self.assertEqual(c.written(), {"sixlines_cm.jsonl"})
        self.assertIn(STDERR, {check.produced for check in registry.get("n04/forms-enumerate-p7").checks})

    def test_procs_scale_the_expected_time(self) -> None:
        c = registry.get("n02/enumerate-p37")
        self.assertEqual(c.procs, 8)
        self.assertEqual(c.expected_seconds(8), c.seconds)
        self.assertEqual(c.expected_seconds(4), 2 * c.seconds)
        self.assertEqual(c.expected_seconds(16), c.seconds)
        self.assertIn("--procs 4", " ".join(c.argv(4)))
        self.assertGreater(c.expected_memory_mb(8), c.expected_memory_mb(2))
        single = registry.get("n27/ledger")
        self.assertEqual(single.expected_seconds(1), single.seconds)
        self.assertEqual(single.argv(4), [])

    def test_select(self) -> None:
        self.assertEqual([c.id for c in registry.select(ids=["n27/ledger"])], ["n27/ledger"])
        enumerations = registry.select(ids=["n02/enumerate-*"])
        self.assertEqual(len(enumerations), 9)
        self.assertTrue(all(c.note == 2 for c in registry.select(notes=[2])))
        self.assertTrue(all("heavy" in c.tags for c in registry.select(tags=["heavy"])))
        self.assertFalse(
            any("heavy" in c.tags for c in registry.select(groups=["lattices"], exclude_tags=["heavy"]))
        )
        self.assertNotIn("n02/enumerate-p37", [c.id for c in registry.select(notes=[2], exclude=["*-p37"])])
        with self.assertRaises(ValueError):
            registry.select(ids=["n99/*"])
        with self.assertRaises(ValueError):
            registry.select(groups=["nonsense"])
        with self.assertRaises(KeyError):
            registry.get("n99/nothing")
