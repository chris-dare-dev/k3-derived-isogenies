"""Every computation behind the notes, as a task: one command of a note's *Reproduce* section.

A computation names the script in `experiments/` and its arguments, the note it belongs to, the claim
family it supports (the `group` of the reports), how long it took when the note was written and with how
many processes, and how to tell whether a rerun reproduced it: which stored files its output is compared
with, in which mode (see `compare.py`), and which lines its output must contain.

Times are those quoted in the notes, on a 16-core machine, or measured where a note gives none.  They are
the task's `cost`, so a batch runs the quick computations first and its `[[limits.tier]]` tables can give
the long ones more time.  `{procs}` in the arguments is the number of processes a task may start; a batch
can lower it, and the expected time then grows in proportion.

A computation runs in a work directory of its own, holding links to the scripts and copies of the stored
results; it writes there and never into `experiments/` (see `workspace.py`).
"""

from __future__ import annotations

import fnmatch
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from . import expectation
from .compare import MODES

STDOUT = "@stdout"
STDERR = "@stderr"
# the keys of a stored record that hold run times
DEFAULT_IGNORE = ("seconds", "setup_s")

GROUPS = {
    "ledger": "Ledger: every headline number recomputed from the stored files (note 27)",
    "lattices": "Six-line surfaces match Hermitian lattices (notes 2 to 5)",
    "ring-steps": "Lattice graphs of ring-changing steps (notes 6 to 12)",
    "fourfold": "Reduction to abelian fourfolds (notes 14 to 16, 18, 21)",
    "group-size": "The group G and its size (note 13)",
    "cm-method": "The CM method computes the action (notes 20, 24, 25)",
    "quadratic-character": "One bit leaks: delta is the character chi_-1 (notes 19, 23, 25, 26)",
    "leakage": "Nothing else leaks: the pre-registered audit (note 28)",
    "order-3": "No character of order 3 is read (note 29)",
    "below-16": "Below Picard number 16: two-neighbours of elliptic K3 surfaces (note 22)",
}


# Scripts whose worker pools read module state the workers inherit by fork, the default start method on
# Linux up to Python 3.13. From Python 3.14 the default is forkserver, under which these pools fail, so
# the harness starts them with the fork method set.
FORK_SCRIPTS = frozenset({"sixlines_fq.py", "sixlines_nonsplit.py"})


@dataclass(frozen=True)
class Check:
    """A file the rerun writes (or its standard output or error), compared with a stored file."""

    produced: str
    stored: str
    mode: str
    ignore_keys: tuple[str, ...] = DEFAULT_IGNORE
    # last lines of the command's output that the stored log is known to lack (mode text-contained)
    unstored_tail: int = 0

    @property
    def label(self) -> str:
        return self.stored


@dataclass(frozen=True)
class Seed:
    """Stored records that a written file starts with, for a command that reads the file it appends to.

    The records of `file` whose `field` is `value` are copied, as they are, into the work directory
    before the command runs.  Only what the command appends after them is compared with the stored file.
    """

    file: str
    field: str
    value: str


@dataclass(frozen=True)
class Computation:
    id: str
    note: int
    group: str
    script: str
    args: tuple[str, ...]
    seconds: float  # wall time when the note was written, with `procs` processes
    procs: int = 1  # processes the command starts at most; the value of {procs}
    memory_mb: int = 500  # expected peak of the main process
    worker_mb: int = 300  # expected peak of each worker process, when procs > 1
    checks: tuple[Check, ...] = ()
    expect: tuple[str, ...] = ()  # regular expressions that must each match a line of the output
    counts: tuple[tuple[str, int], ...] = ()  # regular expression, number of lines it must match
    writes: tuple[str, ...] = ()  # other files with stored names that the command writes
    merge_stderr: bool = False  # the stored log was written with 2>&1
    caveat: str = ""  # what is known to differ, or why nothing is compared
    tags: tuple[str, ...] = field(default=())
    seeds: tuple[Seed, ...] = ()  # stored records a written file starts with (see Seed)

    @property
    def needs_fork(self) -> bool:
        return self.script in FORK_SCRIPTS

    @property
    def parallel(self) -> bool:
        return any("{procs}" in argument for argument in self.args)

    def argv(self, procs: int) -> list[str]:
        return [argument.replace("{procs}", str(procs)) for argument in self.args]

    def expected_seconds(self, procs: int) -> float:
        """The expected wall time with `procs` processes, scaled from the time on record."""
        if not self.parallel or procs >= self.procs:
            return self.seconds
        return self.seconds * self.procs / max(procs, 1)

    def expected_memory_mb(self, procs: int) -> int:
        return self.memory_mb + (self.worker_mb * procs if self.parallel and procs > 1 else 0)

    def written(self) -> set[str]:
        """Files with stored names that the command writes: they are not copied into its work directory."""
        return {c.produced for c in self.checks if not c.produced.startswith("@")} | set(self.writes)

    def command(self, procs: int | None = None) -> str:
        return " ".join(["python", "-I", self.script, *self.argv(procs or self.procs)])

    def predicted(self) -> dict[str, Any]:
        """What a rerun is expected to give: each stored file identical, each line of the output there."""
        return expectation.predicted((check.label for check in self.checks), self.expect, self.counts)

    def as_dict(self) -> dict[str, Any]:
        return {
            "note": self.note,
            "group": self.group,
            "script": self.script,
            "args": list(self.args),
            "seconds": self.seconds,
            "procs": self.procs,
        }


# building the registry ----------------------------------------------------------------------------

_REGISTRY: dict[str, Computation] = {}


def _check(produced: str, stored: str | tuple[str, str], default_mode: str, unstored_tail: int = 0) -> Check:
    if isinstance(stored, tuple):
        return Check(produced, stored[0], stored[1], unstored_tail=unstored_tail)
    return Check(produced, stored, default_mode, unstored_tail=unstored_tail)


def C(
    id: str,
    script: str,
    args: str | list[str],
    seconds: float,
    *,
    procs: int = 1,
    memory_mb: int = 500,
    worker_mb: int = 300,
    stdout: str | tuple[str, str] | None = None,
    stdout_unstored_tail: int = 0,
    stderr: str | tuple[str, str] | None = None,
    outputs: Mapping[str, str | tuple[Any, ...]] | None = None,
    expect: tuple[str, ...] = (),
    counts: dict[str, int] | None = None,
    writes: tuple[str, ...] = (),
    merge_stderr: bool = False,
    caveat: str = "",
    tags: tuple[str, ...] = (),
    seeds: tuple[Seed, ...] = (),
    group: str,
) -> Computation:
    note = int(id.split("/")[0].lstrip("n"))
    checks = []
    if stdout is not None:
        checks.append(_check(STDOUT, stdout, "text", stdout_unstored_tail))
    if stderr is not None:
        checks.append(_check(STDERR, stderr, "text"))
    for produced, target in (outputs or {}).items():
        # a mode alone compares the file with the stored file of the same name
        checks.append(Check(produced, produced, target) if isinstance(target, str) else Check(produced, *target))
    computation = Computation(
        id=id,
        note=note,
        group=group,
        script=script,
        args=tuple(args.split() if isinstance(args, str) else args),
        seconds=seconds,
        procs=procs,
        memory_mb=memory_mb,
        worker_mb=worker_mb,
        checks=tuple(checks),
        expect=expect,
        counts=tuple((counts or {}).items()),
        writes=writes,
        merge_stderr=merge_stderr,
        caveat=caveat,
        tags=tags,
        seeds=seeds,
    )
    if id in _REGISTRY:
        raise ValueError(f"two computations with the id {id}")
    _REGISTRY[id] = computation
    return computation


PASSED = (r"selftest passed",)
HEAVY = ("heavy",)

# Note 27: the ledger ----------------------------------------------------------------------------

C("n27/ledger", "sixlines_ledger.py", [], 1, stdout=("sixlines_ledger.log", "exact"), group="ledger")

# Notes 2 to 5: six-line surfaces and Hermitian lattices ---------------------------------------------

G = "lattices"
C("n02/selftest", "sixlines.py", "selftest", 4, expect=PASSED, group=G)
_ENUMERATE = {7: 2, 11: 3, 13: 5, 17: 10, 19: 15, 23: 30, 29: 90, 31: 130, 37: 420}
for _p, _s in _ENUMERATE.items():
    C(f"n02/enumerate-p{_p}", "sixlines.py", f"enumerate --p {_p} --out sixlines_p{_p}.jsonl --procs {{procs}}",
      _s, procs=8, worker_mb=200, outputs={f"sixlines_p{_p}.jsonl": "records"}, group=G)
_SPLIT = " ".join(f"sixlines_p{p}.jsonl" for p in (7, 11, 13, 17, 19, 23, 29, 31, 37))
C("n02/cm-list", "sixlines_cm.py", f"--list --out sixlines_cm.jsonl {_SPLIT}", 60,
  stdout="sixlines_cm.log", outputs={"sixlines_cm.jsonl": "records"}, group=G)
C("n03/fq-selftest", "sixlines_fq.py", "selftest", 1, expect=PASSED, group=G)
C("n03/enumerate-q25", "sixlines_fq.py", "enumerate --p 5 --out sixlines_q25.jsonl --procs {procs}", 60,
  procs=8, worker_mb=200, outputs={"sixlines_q25.jsonl": "records"}, group=G)
_NEIGHBOURS = " ".join(f"sixlines_p{p}.jsonl" for p in (13, 17, 19, 23, 29, 31, 37)) + " sixlines_q25.jsonl"
C("n03/neighbours", "sixlines_cm.py", f"--list --neighbours --out sixlines_neighbours.jsonl {_NEIGHBOURS}", 180,
  stdout="sixlines_neighbours.log", outputs={"sixlines_neighbours.jsonl": "records"}, group=G)
C("n03/enumerate-q49", "sixlines_fq.py", "enumerate --p 7 --out sixlines_q49.jsonl --procs {procs} --reduced",
  85 * 60, procs=14, worker_mb=250, merge_stderr=True, stdout="sixlines_q49.log",
  outputs={"sixlines_q49.jsonl": "records"}, tags=HEAVY, group=G)
C("n03/q49-cm", "sixlines_cm.py", "--list --neighbours --out sixlines_q49_cm.jsonl sixlines_q49.jsonl", 300,
  stdout="sixlines_q49_cm.log", outputs={"sixlines_q49_cm.jsonl": "records"}, group=G)
C("n04/forms-selftest", "sixlines_forms.py", "selftest", 3, expect=PASSED, group=G)
C("n04/enumerate-p5", "sixlines.py", "enumerate --p 5 --out sixlines_p5.jsonl", 2,
  outputs={"sixlines_p5.jsonl": "records"}, group=G)
for _p, _q, _lat, _s in ((5, "q25", "sixlines_neighbours.jsonl", 20), (7, "q49", "sixlines_q49_cm.jsonl", 60)):
    C(f"n04/forms-enumerate-p{_p}", "sixlines_forms.py", f"enumerate --p {_p} --out sixlines_forms_p{_p}.jsonl",
      _s, stderr=f"sixlines_forms_p{_p}.log", outputs={f"sixlines_forms_p{_p}.jsonl": "records"}, group=G)
    C(f"n04/forms-compare-p{_p}", "sixlines_forms.py",
      f"compare sixlines_forms_p{_p}.jsonl sixlines_{_q}.jsonl --split sixlines_p{_p}.jsonl "
      f"--out sixlines_forms_p{_p}_compare.jsonl", 60, stdout=f"sixlines_forms_p{_p}_compare.log",
      outputs={f"sixlines_forms_p{_p}_compare.jsonl": "records"}, group=G)
    C(f"n04/types-{_q}", "sixlines_types.py", f"{_lat} sixlines_{_q}.jsonl --out sixlines_types_{_q}.jsonl", 60,
      stdout=f"sixlines_types_{_q}.log", outputs={f"sixlines_types_{_q}.jsonl": "records"}, group=G)
C("n05/nonsplit-selftest", "sixlines_nonsplit.py", "selftest", 1, expect=PASSED, group=G)
for _p, _s in {11: 5, 13: 10, 17: 60, 19: 120, 23: 660}.items():
    C(f"n05/nonsplit-enumerate-p{_p}", "sixlines_nonsplit.py",
      f"enumerate --p {_p} --out sixlines_nonsplit_p{_p}.jsonl --procs {{procs}}", _s, procs=12, worker_mb=200,
      stderr=f"sixlines_nonsplit_p{_p}.log", outputs={f"sixlines_nonsplit_p{_p}.jsonl": "records"}, group=G)
C("n05/nonsplit-lattices-small", "sixlines_nonsplit.py",
  "lattices --list --out sixlines_nonsplit_lat_small.jsonl sixlines_forms_p5.jsonl sixlines_forms_p7.jsonl "
  "sixlines_nonsplit_p11.jsonl sixlines_nonsplit_p13.jsonl", 30,
  stdout="sixlines_nonsplit_lat_small.log", outputs={"sixlines_nonsplit_lat_small.jsonl": "records"}, group=G)
for _p, _s in {17: 60, 19: 120, 23: 480}.items():
    C(f"n05/nonsplit-lattices-p{_p}", "sixlines_nonsplit.py",
      f"lattices --list --out sixlines_nonsplit_lat_p{_p}.jsonl sixlines_nonsplit_p{_p}.jsonl", _s,
      stdout=f"sixlines_nonsplit_lat_p{_p}.log", outputs={f"sixlines_nonsplit_lat_p{_p}.jsonl": "records"}, group=G)

# Notes 6 to 12: ring-changing steps, lattice side ---------------------------------------------------

G = "ring-steps"
for _id, _script, _extra, _name, _s, _tags in (
    ("n06/rational", "sixlines_rational.py", "", "sixlines_rational", 180, ()),
    ("n07/vertical", "sixlines_vertical.py", "", "sixlines_vertical", 180, ()),
    ("n08/levels", "sixlines_levels.py", "", "sixlines_levels", 5, ()),
    ("n09/ramified", "sixlines_ramified.py", "", "sixlines_ramified", 120, ()),
    ("n09/ramified-q7", "sixlines_ramified.py", "--search 7,3 --search 7,5", "sixlines_ramified_q7", 1800, HEAVY),
    ("n10/deep", "sixlines_deep.py", "", "sixlines_deep", 780, ()),
    ("n11/rank6", "sixlines_rank6.py", "", "sixlines_rank6", 240, ()),
    ("n12/hermitian", "sixlines_hermitian.py", "", "sixlines_hermitian", 60, ()),
    ("n12/hermitian-q7", "sixlines_hermitian.py", "--search 7,3", "sixlines_hermitian_q7", 420, ()),
):
    C(_id, _script, f"{_extra} --out {_name}.jsonl".strip(), _s, memory_mb=1000,
      stdout=f"{_name}.log", outputs={f"{_name}.jsonl": "records"}, tags=_tags, group=G)

# Note 13: the group G ------------------------------------------------------------------------------

C("n13/groupsize", "sixlines_groupsize.py",
  "--cap 300 --out sixlines_groupsize.jsonl --fourfold 41 61 101 --real --all 13 17 19 23 25 29 31 37 41 43 47 49 "
  "53 59 61 67 71 73 79 83 89 97 101 --sample 211:40 503:40 1009:40 2003:30 5003:30 10007:30 30011:20 100003:12 "
  "125:20 625:20 3125:20 343:20 2401:20", 17 * 60, memory_mb=1500,
  stdout="sixlines_groupsize.log", outputs={"sixlines_groupsize.jsonl": "records"}, group="group-size")
C("n13/groupsize-41", "sixlines_groupsize.py", "--cap 300 --fourfold 41 --real", 60, memory_mb=1000,
  caveat="a rerun quoted in note 18 section 2; its output was not kept", group="group-size")
C("n13/groupsize-all-small", "sixlines_groupsize.py", "--cap 300 --all 13 17 19 23", 30, memory_mb=1000,
  caveat="a rerun quoted in note 18 section 2; its output was not kept", group="group-size")

# Notes 14 to 16, 18 and 21: the fourfold -------------------------------------------------------------

G = "fourfold"
_PRYM_PRIMES = "5 7 11 13 17 19 23 29 31 37"
C("n14/prym-check-17", "sixlines_prym.py", "--check 17 2 3 4 11", 3,
  counts={r"^\s*pair .*\bequal$": 18, r"^\s*pair .*\btwist$": 2, r"^\s*pair .*\bdegenerate$": 10,
          r"mismatch": 0},
  caveat="no stored log; DATA.md records 18 equal, 2 twist, 10 degenerate, no mismatch", group=G)
C("n14/prym-run-11-13", "sixlines_prym.py", "--workers {procs} --run 11 13 --out sixlines_prym_11_13.jsonl", 20,
  procs=4, outputs={"sixlines_prym_11_13.jsonl": ("sixlines_prym.jsonl", "records-subset")},
  counts={r"^p = 1[13]: ": 2}, group=G)
C("n14/prym-run", "sixlines_prym.py", f"--workers {{procs}} --out sixlines_prym.jsonl --run {_PRYM_PRIMES}", 1800,
  procs=14, stdout="sixlines_prym.log", outputs={"sixlines_prym.jsonl": "records"}, tags=HEAVY, group=G)
C("n14/prym-summary", "sixlines_prym.py", "--summary sixlines_prym.jsonl", 5,
  caveat="prints the statistics of the stored file; no stored output to compare with", group=G)
C("n14/prym-other-half", "sixlines_prym.py", "--workers {procs} --other-half --run 13 17 19", 300, procs=8,
  caveat="the control of note 14 section 5.3; its output was not kept", group=G)
C("n15/dictionary", "sixlines_dictionary.py", "dictionary --workers {procs} --bound 2000 --out sixlines_dictionary.jsonl",
  240, procs=8, stdout="sixlines_dictionary.log", outputs={"sixlines_dictionary.jsonl": "records"}, group=G)
C("n15/association", "sixlines_dictionary.py", f"association {_PRYM_PRIMES}", 240,
  stdout="sixlines_dictionary_association.log", group=G)
C("n16/prym-fq-q25", "sixlines_prym_fq.py", "--workers {procs} --select all --pairs 0 --run 25 --out sixlines_prym_q25.jsonl",
  600, procs=14, stdout="sixlines_prym_q25.log", outputs={"sixlines_prym_q25.jsonl": "records"}, group=G)
C("n16/prym-fq-q49", "sixlines_prym_fq.py",
  "--workers {procs} --select irreducible --pairs 3 --run 49 --out sixlines_prym_q49.jsonl", 44 * 60, procs=14,
  stdout="sixlines_prym_q49.log", outputs={"sixlines_prym_q49.jsonl": "records"}, tags=HEAVY, group=G)
C("n16/prym-fq-13-17", "sixlines_prym_fq.py", "--workers {procs} --select all --pairs 0 --run 13 17", 20, procs=12,
  caveat="its output was not kept", group=G)
for _name, _args, _s, _procs in (
    ("sign-minus", "sign-minus sixlines_prym.jsonl", 5, 1),
    ("inverting-7-4", "inverting 7 4", 30, 1),
    ("step-8", "step 8", 60, 1),
    ("other-half", "other-half 13 17 19 --workers {procs}", 150, 10),
    ("kernel", "kernel --workers {procs}", 180, 10),
    ("inverting-11-1", "inverting 11 1", 480, 1),
):
    C(f"n18/review-{_name}", "sixlines_review.py", _args, _s, procs=_procs,
      stdout=("sixlines_review.log", "text-contained"), group=G)
_FP_LOG = ("sixlines_forms_fp.log", "text-contained")
C("n21/forms-fp-selftest", "sixlines_forms_fp.py", "selftest", 1, stdout=_FP_LOG, group=G)
C("n21/forms-fp-commuting", "sixlines_forms_fp.py", "commuting sixlines_prym.jsonl", 5, stdout=_FP_LOG, group=G)
C("n21/forms-fp-types", "sixlines_forms_fp.py", f"types {_PRYM_PRIMES} --workers {{procs}}", 10, procs=4,
  stdout=_FP_LOG, group=G)
C("n21/forms-fp-involution", "sixlines_forms_fp.py", "involution 11 13 17 19 23 29 31 37 --workers {procs}", 90,
  procs=4, caveat="uses random points: the counts of points vary between runs, so the log is not compared",
  group=G)
for _name, _args, _s in (
    ("k2-7-17", "7 11 13 17 --pairs 0 --k 2", 15 * 60),
    ("k2-19", "19 --pairs 6 --k 2", 10 * 60),
    ("k2-23", "23 --pairs 2 --limit 150 --k 2", 16 * 60),
    ("k2-29", "29 --pairs 2 --limit 60 --k 2", 19 * 60),
    ("k2-37", "37 --pairs 2 --limit 20 --k 2", 18 * 60),
    ("k3-7", "7 --pairs 0 --k 3", 2 * 60),
    ("k3-11", "11 --pairs 2 --k 3", 3 * 60),
    ("k3-13", "13 --pairs 1 --limit 12 --k 3", 13 * 60),
):
    C(f"n21/inverting-{_name}", "sixlines_forms_fp.py", f"inverting {_args} --workers {{procs}} --out sixlines_forms_fp.jsonl",
      _s, procs=4, stdout=_FP_LOG, outputs={"sixlines_forms_fp.jsonl": ("sixlines_forms_fp.jsonl", "records-subset")},
      tags=HEAVY if _s > 15 * 60 else (), group=G)

# Notes 20, 24 and 25: the CM method -------------------------------------------------------------------

G = "cm-method"
_CMJ = {"sixlines_cmmethod.jsonl": ("sixlines_cmmethod.jsonl", "records-subset")}
C("n20/selftest", "sixlines_cmmethod.py", "selftest --digits 60", 1,
  stdout=("sixlines_cmmethod_selftest.log", "exact"), group=G)
C("n20/sweep", "sixlines_cmmethod.py",
  "sweep 17 19 23 29 31 37 --digits-list 300 700 1200 --primes 3 --all-subsets --maxn 8 --workers {procs} "
  "--out sixlines_cmmethod.jsonl", 45 * 60, procs=4, worker_mb=1000, stdout="sixlines_cmmethod.log", outputs=_CMJ,
  caveat="the script was corrected after this sweep (note 20 section 7): the surfaces of a polynomial are now "
  "selected by e1 and e2, and the ratios with J3 were added, so the records of q = 31 and 37 and the log are "
  "expected to differ", tags=HEAVY, group=G)
C("n20/sweep-rerun", "sixlines_cmmethod.py",
  "sweep 29 31 37 --select 29,6 29,10 31,2 31,22 37,18 37,22 --digits-list 300 700 1200 --primes 6 --all-subsets "
  "--maxn 8 --workers {procs} --out sixlines_cmmethod.jsonl", 30 * 60, procs=4, worker_mb=1000,
  stdout="sixlines_cmmethod_rerun.log", outputs=_CMJ, tags=HEAVY, group=G)
C("n20/sweep-rerun2", "sixlines_cmmethod.py",
  "sweep 29 --select 29,10 --digits-list 1200 --primes 6 --all-subsets --workers 1 --out sixlines_cmmethod.jsonl",
  12 * 60, memory_mb=1500, stdout="sixlines_cmmethod_rerun2.log", outputs=_CMJ, group=G)
C("n20/points-maximal", "sixlines_cmmethod.py", "points 17 2 1 1 --digits 30 --maximal", 120,
  caveat="its output was not kept", group=G)
C("n20/points-nolevel", "sixlines_cmmethod.py", "points 17 2 1 1 --digits 30 --nolevel", 120,
  caveat="its output was not kept", group=G)
C("n20/action-23", "sixlines_cmmethod.py", "action 23 2 1 1 --digits 700 --primes 3 --all-subsets", 300,
  memory_mb=1500, caveat="its output was not kept", group=G)
C("n25/rerun-sweep-17-23", "sixlines_cmmethod.py",
  "sweep 17 23 --select 17,2 23,10 --digits-list 300 700 --primes 3 --all-subsets --maxn 8 --workers {procs} "
  "--out sixlines_cmmethod_17_23.jsonl", 300, procs=2, worker_mb=1000,
  caveat="a rerun quoted in note 25 section 2; its output was not kept", group=G)
for _args, _s in (("37 18 1 1 --digits 1200", 83), ("29 2 1 1 --digits 1200", 109), ("37 2 1 1 --digits 2000", 395)):
    _q, _e = _args.split()[:2]
    C(f"n24/delta-labels-{_q}-{_e}", "sixlines_delta_labels.py",
      f"run {_args} --primes 8 --workers {{procs}} --out sixlines_delta_labels.jsonl", _s, procs=14, worker_mb=500,
      stdout=("sixlines_delta_labels.log", "text-contained"),
      outputs={"sixlines_delta_labels.jsonl": ("sixlines_delta_labels.jsonl", "records-subset")}, group=G)
_R2_LOG = ("sixlines_review2.log", "text-contained")
_R2J = {"sixlines_review2.jsonl": ("sixlines_review2.jsonl", "records-subset")}
C("n25/review2-theta", "sixlines_review2.py", "theta", 15, stdout=_R2_LOG, group=G)
for _args, _s in (("37 18 1 1 --digits 1200", 60), ("29 2 1 1 --digits 1200", 90), ("37 2 1 1 --digits 2000", 360)):
    _q, _e = _args.split()[:2]
    C(f"n25/review2-galois-{_q}-{_e}", "sixlines_review2.py",
      f"galois {_args} --workers {{procs}} --out sixlines_review2.jsonl", _s, procs=8, worker_mb=500,
      stdout=_R2_LOG, outputs=_R2J, group=G)

# Notes 19, 23, 25 and 26: the quadratic character ------------------------------------------------------

G = "quadratic-character"
C("n19/characters-surfaces", "sixlines_characters.py", "surfaces --out sixlines_characters.jsonl", 20,
  stdout="sixlines_characters.log", outputs={"sixlines_characters.jsonl": "records"}, group=G)
C("n19/characters-genus", "sixlines_characters.py",
  "genus --workers {procs} --max-q 10007 --out sixlines_characters_genus.jsonl", 40, procs=8,
  stdout="sixlines_characters_genus.log", outputs={"sixlines_characters_genus.jsonl": "records"}, group=G)
C("n19/characters-dyadic", "sixlines_characters.py", "dyadic sixlines_characters_genus.jsonl", 5,
  stdout="sixlines_characters_dyadic.log", group=G)
C("n25/review2-stored", "sixlines_review2.py", "stored", 1, stdout=_R2_LOG, group=G)
C("n25/review2-levels", "sixlines_review2.py", "levels", 1, stdout=_R2_LOG, group=G)
C("n25/review2-pairing", "sixlines_review2.py", "pairing 29 2 3 12 27 1 3 --seed 3", 10, stdout=_R2_LOG, group=G)
C("n25/review2-delta", "sixlines_review2.py", "delta --tw 3 --seed 2026 --workers {procs} --out sixlines_review2.jsonl",
  60, procs=8, stdout=_R2_LOG, outputs=_R2J, group=G)
_JAC_RANDOM = ("the classes are random (seeded): the counts of generators tried can differ between runs, "
               "the group structures and verdicts do not")
C("n23/jacobian-validate", "sixlines_jacobian.py", "--validate-set --workers {procs} --out sixlines_jacobian_validate.jsonl",
  300, procs=4, worker_mb=600, stdout="sixlines_jacobian_validate.log",
  outputs={"sixlines_jacobian_validate.jsonl": "records"}, caveat=_JAC_RANDOM, group=G)
C("n23/jacobian-lattice", "sixlines_jacobian.py", "--lattice", 60, stdout="sixlines_jacobian_lattice.log", group=G)
C("n23/jacobian-levels", "sixlines_jacobian.py", "--levels --workers {procs} --out sixlines_jacobian_levels.jsonl",
  420, procs=4, worker_mb=600, stdout="sixlines_jacobian_levels.log",
  outputs={"sixlines_jacobian_levels.jsonl": "records"}, caveat=_JAC_RANDOM, group=G)
C("n23/jacobian-character", "sixlines_jacobian.py",
  "--character --pairs 3 --workers {procs} --out sixlines_jacobian_character.jsonl", 50 * 60, procs=4, worker_mb=600,
  stdout="sixlines_jacobian_character.log", outputs={"sixlines_jacobian_character.jsonl": "records"},
  caveat=_JAC_RANDOM, tags=HEAVY, group=G)
C("n23/jacobian-character-extra", "sixlines_jacobian.py",
  "--character-extra sixlines_jacobian_character.jsonl --workers {procs} --out sixlines_jacobian_character_extra.jsonl",
  35 * 60, procs=4, worker_mb=600, stdout="sixlines_jacobian_character_extra.log",
  outputs={"sixlines_jacobian_character_extra.jsonl": "records"},
  caveat=_JAC_RANDOM + "; the stored run was made in two pieces and joined, so its log differs from one run's",
  tags=HEAVY, group=G)
C("n23/jacobian-one-29", "sixlines_jacobian.py", "--one 29 2 3 5 20 1 2", 30,
  caveat="a rerun quoted in note 25 section 2; its output was not kept", group=G)
_D41_LOG = ("sixlines_delta41.log", "text-contained")
C("n26/enumerate-p41", "sixlines.py", "enumerate --p 41 --out sixlines_p41.jsonl --procs {procs}", 17 * 60, procs=14,
  worker_mb=200, outputs={"sixlines_p41.jsonl": "records"}, tags=HEAVY, group=G)
C("n26/cm-list-p41", "sixlines_cm.py", "--list --out sixlines_cm_p41.jsonl sixlines_p41.jsonl", 11,
  stdout="sixlines_cm_p41.log", outputs={"sixlines_cm_p41.jsonl": "records"}, group=G)
C("n26/dyadic", "sixlines_delta41.py", "dyadic", 1, caveat="no stored log", group=G)
C("n26/selfcheck", "sixlines_delta41.py", "selfcheck --n 3 --workers {procs}", 120, procs=6, stdout=_D41_LOG, group=G)
C("n26/list", "sixlines_delta41.py", "list", 1, stdout=_D41_LOG, stdout_unstored_tail=1, group=G,
  caveat="the stored log has the first 31 of the 32 lines that `list` prints: it lists surfaces 0 to 10 of the "
         "e1 = 22 torsor, which its own heading gives 12 surfaces, and lacks the line of surface 11.  The 31 "
         "lines are compared; the last line of the rerun is not")
C("n26/prym", "sixlines_delta41.py", "prym --pairs 0 --workers {procs}", 6 * 60, procs=14,
  outputs={"sixlines_prym_p41.jsonl": "records"}, group=G)
C("n26/validate", "sixlines_delta41.py", "validate --n 2", 120, stdout=_D41_LOG, group=G)
_DELTA = {"sixlines_delta41.jsonl": ("sixlines_delta41.jsonl", "records-subset")}
_DELTA_NOTE = ("delta appends to sixlines_delta41.jsonl and skips what is there; the rerun starts from an empty "
               "file, so it recomputes every record of its reading")
_DELTA_SEED = (Seed("sixlines_delta41.jsonl", "reading", "primary"),)
_DELTA_SEEDED_NOTE = ("delta appends to sixlines_delta41.jsonl and skips what is there; this reading chooses its pairs "
                      "from the records of the primary reading, so the rerun starts from the stored primary records "
                      "and only the records it appends are compared with the stored file")
C("n26/delta-primary", "sixlines_delta41.py", "delta --reading primary --pairs 6 --workers {procs}", 30 * 60,
  procs=12, worker_mb=500, outputs=_DELTA, caveat=_DELTA_NOTE, tags=HEAVY, group=G)
C("n26/delta-primary-surfaces", "sixlines_delta41.py",
  "delta --reading primary --pairs 0 --surfaces 6:7 6:15 6:2 6:3 22:3 --workers {procs}", 18 * 60, procs=14,
  worker_mb=500, outputs=_DELTA, caveat=_DELTA_NOTE, tags=HEAVY, group=G)
for _reading in ("twist1", "twist2", "seed", "n16", "r2"):
    C(f"n26/delta-{_reading}", "sixlines_delta41.py", f"delta --reading {_reading} --per-surface 1 --workers {{procs}}",
      160, procs=14, worker_mb=500, outputs=_DELTA, caveat=_DELTA_SEEDED_NOTE, seeds=_DELTA_SEED, tags=HEAVY,
      group=G)
C("n26/delta-cofactor", "sixlines_delta41.py", "delta --reading cofactor --pairs 2 --workers {procs}", 11 * 60,
  procs=14, worker_mb=500, outputs=_DELTA, caveat=_DELTA_SEEDED_NOTE, seeds=_DELTA_SEED, tags=HEAVY, group=G)
_COSETS41 = {"sixlines_delta41_cosets.jsonl": ("sixlines_delta41_cosets.jsonl", "records-subset")}
C("n26/cosets-22", "sixlines_delta41.py", "cosets 41 22 1 1 --digits 1200 2000 3000 --primes 8 --workers {procs} --no-delta",
  30 * 60, procs=4, memory_mb=1500, worker_mb=1000, stdout="sixlines_delta41_cosets_22.log", outputs=_COSETS41,
  tags=HEAVY, group=G)
C("n26/cosets-6", "sixlines_delta41.py", "cosets 41 6 1 1 --digits 3000 4200 --primes 8 --workers {procs} --no-delta",
  52 * 60, procs=8, memory_mb=1500, worker_mb=1000, stdout="sixlines_delta41_cosets_6.log", outputs=_COSETS41,
  tags=HEAVY, group=G)
C("n26/compare", "sixlines_delta41.py", "compare", 1, stdout=("sixlines_delta41_compare.log", "text-contained"),
  outputs={"sixlines_delta41_compare.json": ("sixlines_delta41_compare.json", "json")}, group=G)
for _args in ("37 18 1 1 --digits 1200", "29 2 1 1 --digits 1200", "37 2 1 1 --digits 2000"):
    _q, _e = _args.split()[:2]
    C(f"n26/signs-cosets-{_q}-{_e}", "sixlines_delta41.py",
      f"cosets {_args} --primes 3 --frob-primes 2 --workers {{procs}} --out sixlines_delta41_signs.jsonl", 180,
      procs=14, worker_mb=500, stdout=("sixlines_delta41_signs.log", "text-contained"),
      outputs={"sixlines_delta41_signs.jsonl": ("sixlines_delta41_signs.jsonl", "records-subset")}, group=G)
C("n26/signs", "sixlines_delta41.py", "signs", 1, stdout=("sixlines_delta41_compare.log", "text-contained"), group=G)

# Note 28: the leakage audit --------------------------------------------------------------------------

G = "leakage"
C("n28/features", "sixlines_leakage.py", "features", 4, memory_mb=800, stdout="sixlines_leakage_features.log",
  outputs={"sixlines_leakage_features.json": ("sixlines_leakage_features.json", "json")}, group=G)
C("n28/audit", "sixlines_leakage.py", "audit --workers {procs}", 129, procs=12, memory_mb=1000, worker_mb=400,
  stdout="sixlines_leakage.log", outputs={"sixlines_leakage.jsonl": "records"}, group=G)
C("n28/extra", "sixlines_leakage.py", "extra --workers {procs}", 12 * 60, procs=12, memory_mb=1000, worker_mb=400,
  stdout="sixlines_leakage_extra.log", outputs={"sixlines_leakage_extra.jsonl": "records"}, tags=HEAVY, group=G)

# Note 29: characters of order 3 ---------------------------------------------------------------------------

G = "order-3"
# The script writes the key identified_by in every record since the q = 43 reruns of section 7.4.  Of the
# 18 stored records two have it: the run with --identify-by-coordinate and the 4,500-digit run with --mix.
# For those two the key is compared; for every other run the stored record predates the key.
_COSETS3 = {"sixlines_order3_cosets.jsonl": ("sixlines_order3_cosets.jsonl", "records-subset")}
_NO_KEY = ("the script now writes the key identified_by, added for the q = 43 reruns of section 7.4, which the "
           "stored record of this run lacks, so that key is not compared")
C("n29/enumerate-p43", "sixlines.py", "enumerate --p 43 --out sixlines_p43.jsonl --procs {procs}", 25 * 60, procs=14,
  worker_mb=200, outputs={"sixlines_p43.jsonl": "records"}, tags=HEAVY, group=G)
C("n29/cm-list-p43", "sixlines_cm.py", "--list --out sixlines_cm_p43.jsonl sixlines_p43.jsonl", 15,
  stdout="sixlines_cm_p43.log", outputs={"sixlines_cm_p43.jsonl": "records"}, group=G)
C("n29/fields", "sixlines_order3.py", "fields", 10, stdout="sixlines_order3_fields.log", group=G)
_EARLY_COSETS = {"sixlines_order3_cosets.jsonl": ("sixlines_order3_cosets.jsonl", "records-subset",
                                                  ("seconds", "identified_by"))}
for _args, _s in (("23 2 1 1 --digits 700", 60), ("29 6 1 1 --digits 700", 90), ("29 2 1 1 --digits 1200", 300)):
    _q, _e = _args.split()[:2]
    C(f"n29/cosets-{_q}-{_e}", "sixlines_order3.py", f"cosets {_args}", _s, memory_mb=1000, outputs=_EARLY_COSETS,
      caveat="no log was kept of this run; " + _NO_KEY, group=G)
for _args, _log, _s, _procs in (
    ("37 2 1 1 --digits 1200 2000", "sixlines_order3_cosets_37.log", 600, 6),
    ("37 6 1 1 --digits 1200 2000", "sixlines_order3_cosets_37.log", 600, 6),
    ("41 22 1 1 --digits 2000", "sixlines_order3_cosets_41.log", 600, 6),
    ("41 2 1 1 --digits 2000", "sixlines_order3_cosets_41.log", 600, 6),
    ("43 14 1 1 --digits 2000 3000", "sixlines_order3_cosets_43.log", 900, 6),
    ("43 22 1 1 --digits 2000 3000", "sixlines_order3_cosets_43.log", 900, 6),
    ("43 14 1 1 --digits 3000 4500", "sixlines_order3_cosets_43_14.log", 1260, 6),
    ("43 2 1 1 --digits 3000 4500", "sixlines_order3_cosets_43_2.log", 1260, 7),
):
    _q, _e = _args.split()[:2]
    _digits = "-".join(_args.split("--digits ")[1].split())
    C(f"n29/cosets-{_q}-{_e}-{_digits}", "sixlines_order3.py", f"cosets {_args} --workers {{procs}}", _s,
      procs=_procs, memory_mb=1500, worker_mb=1000, stdout=(_log, "text-contained"), outputs=_EARLY_COSETS,
      caveat=_NO_KEY, tags=HEAVY, group=G)
for _name, _extra, _log in (
    ("mix", "--digits 3000 --mix 0 2 -1 1 0", "sixlines_order3_cosets_43_2_mix.log"),
    ("coord", "--digits 3000 --mix 0 2 -1 1 0 --identify-by-coordinate", "sixlines_order3_cosets_43_2_coord.log"),
    ("mix4500", "--digits 4500 --mix 0 2 -1 1 0", "sixlines_order3_cosets_43_2_mix4500.log"),
):
    C(f"n29/cosets-43-2-{_name}", "sixlines_order3.py", f"cosets 43 2 1 1 {_extra} --workers {{procs}}", 1260, procs=4,
      memory_mb=1500, worker_mb=1000, stdout=_log, outputs=_EARLY_COSETS if _name == "mix" else _COSETS3,
      caveat="the stored run read the period points of the 3000/4500-digit run from sixlines_order3_points/, "
      "which is not in the archive; the rerun computes them again, so it takes longer and its log may differ "
      "in the lines about the cache" + ("; " + _NO_KEY if _name == "mix" else ""), tags=HEAVY, group=G)
C("n29/labelcheck", "sixlines_order3.py", "labelcheck", 1, stdout="sixlines_order3_labelcheck.log", group=G)
C("n29/features", "sixlines_order3.py", "features", 12, memory_mb=800, stdout="sixlines_order3_features.log",
  outputs={"sixlines_order3_features.json": ("sixlines_order3_features.json", "json")}, group=G)
C("n29/selfcheck", "sixlines_order3.py", "selfcheck", 3, memory_mb=800, stdout="sixlines_order3_selfcheck.log",
  group=G)
C("n29/audit", "sixlines_order3.py", "audit --workers {procs}", 56, procs=10, memory_mb=1000, worker_mb=400,
  stdout="sixlines_order3_audit.log", outputs={"sixlines_order3.jsonl": "records"}, group=G)
C("n29/features-q43", "sixlines_order3.py", "features --with-43 --out sixlines_order3_features_q43.json", 12,
  memory_mb=800, stdout="sixlines_order3_features_q43.log",
  outputs={"sixlines_order3_features_q43.json": ("sixlines_order3_features_q43.json", "json")}, group=G)
C("n29/selfcheck-q43", "sixlines_order3.py", "selfcheck --features sixlines_order3_features_q43.json", 3,
  memory_mb=800, stdout="sixlines_order3_selfcheck_q43.log", group=G)
C("n29/audit-q43", "sixlines_order3.py",
  "audit --workers {procs} --features sixlines_order3_features_q43.json --out sixlines_order3_q43.jsonl", 70,
  procs=7, memory_mb=3000, worker_mb=400, stdout="sixlines_order3_q43.log",
  outputs={"sixlines_order3_q43.jsonl": "records"}, group=G)
C("n29/holdout43", "sixlines_order3_holdout43.py", [], 5, memory_mb=800, stdout="sixlines_order3_holdout43.log",
  group=G)

# Note 22: below Picard number 16 ---------------------------------------------------------------------------

G = "below-16"
for _name, _args, _core_hours, _procs in (
    ("e8a2-p5", "--family e8a2 --p 5 --count 16000 --check2 5", 2.3, 10),
    ("e8a2-p7", "--family e8a2 --p 7 --count 4000 --check2 5", 5.4, 10),
    ("e8-p5", "--family e8 --p 5 --count 1200", 0.8, 10),
    ("e8a2-p5-control", "--family e8a2 --p 5 --count 1500 --control", 0.1, 10),
    ("e8-p5-control", "--family e8 --p 5 --count 600 --control", 0.05, 10),
    ("generic-p5", "--family generic --p 5 --count 10 --check2 7 --bnf-limit 900", 0.9, 5),
    ("generic-p5-seed2", "--family generic --p 5 --count 30 --seed 2 --check2 7 --bnf-limit 900", 2.8, 10),
    ("generic-p5-control", "--family generic --p 5 --count 10 --control --bnf-limit 60", 2.5, 5),
):
    _file = "ellk3_" + _name.replace("-", "_")
    C(f"n22/{_name}", "ellk3_two_neighbours.py", f"{_args} --workers {{procs}} --out {_file}.jsonl",
      round(_core_hours * 3600 / _procs), procs=_procs, worker_mb=400, stdout=f"{_file}.log",
      outputs={f"{_file}.jsonl": "records"}, tags=HEAVY if _core_hours > 0.2 else (), group=G)
for _p in (5, 7):
    C(f"n22/graph-e8a2-p{_p}", "ellk3_two_neighbours.py", f"--graph ellk3_e8a2_p{_p}.jsonl", 10,
      stdout=("ellk3_e8a2_graph.log", "text-contained"), group=G)


# looking computations up ------------------------------------------------------------------------------------


def all_computations() -> list[Computation]:
    return list(_REGISTRY.values())


def get(computation_id: str) -> Computation:
    try:
        return _REGISTRY[computation_id]
    except KeyError:
        raise KeyError(f"no computation {computation_id!r}; `python3 -m k3tasks list` lists them") from None


def select(
    *,
    ids: list[str] | None = None,
    notes: list[int] | None = None,
    groups: list[str] | None = None,
    tags: list[str] | None = None,
    exclude: list[str] | None = None,
    exclude_tags: list[str] | None = None,
) -> list[Computation]:
    """The computations matching every selector given; an id may be a shell pattern ("n02/enumerate-*")."""
    chosen = all_computations()
    if ids is not None:
        unmatched = [pattern for pattern in ids if not any(fnmatch.fnmatchcase(c.id, pattern) for c in chosen)]
        if unmatched:
            raise ValueError(f"no computation matches {unmatched[0]!r}")
        chosen = [c for c in chosen if any(fnmatch.fnmatchcase(c.id, pattern) for pattern in ids)]
    if notes is not None:
        chosen = [c for c in chosen if c.note in notes]
    if groups is not None:
        unknown = set(groups) - set(GROUPS)
        if unknown:
            raise ValueError(f"unknown group {sorted(unknown)[0]!r}; the groups are {', '.join(GROUPS)}")
        chosen = [c for c in chosen if c.group in groups]
    if tags is not None:
        chosen = [c for c in chosen if set(tags) & set(c.tags)]
    if exclude_tags:
        chosen = [c for c in chosen if not set(exclude_tags) & set(c.tags)]
    if exclude:
        chosen = [c for c in chosen if not any(fnmatch.fnmatchcase(c.id, pattern) for pattern in exclude)]
    return chosen


def validate() -> list[str]:
    """Problems with the registry itself; empty when it is consistent."""
    problems = []
    for c in all_computations():
        if c.group not in GROUPS:
            problems.append(f"{c.id}: unknown group {c.group}")
        for check in c.checks:
            if check.mode not in MODES:
                problems.append(f"{c.id}: unknown mode {check.mode}")
            if check.produced == STDERR and c.merge_stderr:
                problems.append(f"{c.id}: compares the standard error, but it is merged into the output")
            if check.unstored_tail and (check.mode != "text-contained" or not c.caveat):
                problems.append(f"{c.id}: lines are left out of a comparison without the mode or the caveat for it")
        compared = {check.produced: check for check in c.checks}
        for seed in c.seeds:
            if seed.file not in compared or compared[seed.file].mode not in ("records", "records-subset"):
                problems.append(f"{c.id}: seeds {seed.file}, which it does not compare as records")
            if not c.caveat:
                problems.append(f"{c.id}: seeds {seed.file} without a caveat that says so")
        if c.parallel and c.procs < 2:
            problems.append(f"{c.id}: takes {{procs}} but records one process")
        if not c.parallel and c.procs != 1:
            problems.append(f"{c.id}: records {c.procs} processes but takes no {{procs}}")
        if not c.checks and not c.expect and not c.counts and not c.caveat:
            problems.append(f"{c.id}: nothing is checked and no caveat says why")
    return problems
