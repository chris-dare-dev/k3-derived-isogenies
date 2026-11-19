"""What the harness shares with the scripts it runs, and the checks with what they check: the statements
of CLAIMS.md ("What the scripts share"), kept true."""

from __future__ import annotations

import sys
import unittest

from taskbatch.settings import REPOSITORY

sys.path.insert(0, str(REPOSITORY / "tools"))

import claims

EXPERIMENTS = REPOSITORY / "experiments"
HARNESS = {"k3tasks", "taskbatch", "runners", "task_runner", "claims", "adopt", "tools"}


class Independence(unittest.TestCase):
    def test_no_script_of_experiments_imports_the_harness(self) -> None:
        scripts = sorted(EXPERIMENTS.glob("*.py"))
        self.assertEqual(len(scripts), 30)
        for script in scripts:
            self.assertFalse(claims.imports(script) & HARNESS, script.name)

    def test_the_harness_imports_no_script_and_nothing_outside_the_standard_library(self) -> None:
        scripts = {path.stem for path in EXPERIMENTS.glob("*.py")}
        modules = [*(REPOSITORY / "k3tasks").glob("*.py"), *(REPOSITORY / "taskbatch").glob("*.py")]
        modules += [REPOSITORY / "runners" / "task_runner.py", *(REPOSITORY / "tools").glob("*.py")]
        for module in modules:
            found = claims.imports(module)
            self.assertFalse(found & scripts, module.name)
            self.assertLessEqual(
                found, sys.stdlib_module_names | {"k3tasks", "taskbatch", "adopt"}, module.name
            )

    def test_the_task_registry_uses_no_private_name_of_the_scheduler(self) -> None:
        for module in [*(REPOSITORY / "k3tasks").glob("*.py"), REPOSITORY / "runners" / "task_runner.py"]:
            text = module.read_text(encoding="utf-8")
            for line in text.splitlines():
                if "taskbatch" in line or "scheduler_report" in line:
                    self.assertNotRegex(line, r"(taskbatch[\w.]*|scheduler_report)\._\w", module.name)
                    self.assertNotRegex(line, r"from taskbatch[\w.]* import .*\b_\w", module.name)

    def test_the_ledger_shares_no_code_with_the_scripts_it_recounts(self) -> None:
        scripts = {path.stem for path in EXPERIMENTS.glob("*.py")}
        self.assertEqual(claims.imports(EXPERIMENTS / "sixlines_ledger.py") & scripts, set())

    def test_the_re_derivation_checks_import_the_scripts_they_check(self) -> None:
        scripts = {path.stem for path in EXPERIMENTS.glob("*.py")}
        self.assertEqual(
            claims.imports(EXPERIMENTS / "sixlines_review.py") & scripts,
            {"sixlines_cm", "sixlines_dictionary", "sixlines_prym"},
        )
        self.assertEqual(
            claims.imports(EXPERIMENTS / "sixlines_review2.py") & scripts,
            {"sixlines_cmmethod", "sixlines_delta_labels", "sixlines_jacobian"},
        )
        self.assertIn("sixlines_review2", claims.imports(EXPERIMENTS / "sixlines_delta41.py"))


if __name__ == "__main__":
    unittest.main()
