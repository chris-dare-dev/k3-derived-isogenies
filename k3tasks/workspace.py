"""Running one computation in a work directory of its own, away from the stored results.

The scripts find their inputs next to themselves (`os.path.dirname(os.path.abspath(__file__))`) and write
their outputs there or in the current directory, and several append to or overwrite stored files.  So a
computation never runs in `experiments/`.  Its work directory, under `results/work/`, holds

- a symbolic link to every script of `experiments/` (the scripts only read each other), and
- a copy of every stored result file, except the files the computation writes: those start absent, as
  they would on a fresh run.

Since `abspath` does not resolve links, a script run through its link takes the work directory as its
own: it reads the copies and writes next to them.  The stored files are never opened for writing.  As a
second guard, the size and modification time of every stored file are taken before and after the run,
and a change is an error.

The work directory is deleted when the rerun reproduced everything it was compared with, and kept for
inspection otherwise (or always, with `keep_work`).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .compare import IDENTICAL, compare, sha256
from .registry import STDERR, STDOUT, Computation

REPOSITORY = Path(__file__).resolve().parent.parent
DATA_SUFFIXES = (".jsonl", ".json", ".log")
SCRIPT_SUFFIXES = (".py", ".txt")
# a stored file that every installed copy of the data has; its absence means the archive is not unpacked
DATA_MARKER = "sixlines_ledger.log"


class WorkspaceError(RuntimeError):
    """The computation cannot be run here: missing data, missing script, or a stored file changed."""


@dataclass(frozen=True)
class Places:
    scripts: Path  # experiments/: the scripts
    data: Path  # where the stored results are unpacked (experiments/ by default)
    work: Path  # parent of the work directories
    python: str  # the interpreter that runs the scripts (needs numpy and cypari2)

    @classmethod
    def from_environment(cls) -> Places:
        """The default places, each of which an environment variable can move.

        K3_SCRIPTS_DIR (experiments/), K3_DATA_DIR (experiments/), K3_WORK_DIR (results/work),
        K3_PYTHON (the interpreter running this code).
        """
        scripts = Path(os.environ.get("K3_SCRIPTS_DIR", REPOSITORY / "experiments"))
        return cls(
            scripts=scripts,
            data=Path(os.environ.get("K3_DATA_DIR", scripts)),
            work=Path(os.environ.get("K3_WORK_DIR", REPOSITORY / "results" / "work")),
            python=os.environ.get("K3_PYTHON", sys.executable),
        )


def stored_files(data: Path) -> list[Path]:
    return sorted(p for p in data.iterdir() if p.is_file() and p.suffix in DATA_SUFFIXES)


def _signature(paths: list[Path]) -> dict[str, tuple[int, int]]:
    signature = {}
    for path in paths:
        status = path.stat()
        signature[path.name] = (status.st_size, status.st_mtime_ns)
    return signature


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")


def check_installed(computation: Computation, places: Places) -> None:
    script = places.scripts / computation.script
    if not script.is_file():
        raise WorkspaceError(f"no script {script}")
    if not (places.data / DATA_MARKER).is_file():
        raise WorkspaceError(
            f"the stored results are not installed in {places.data}: unpack the data archive as DATA.md "
            "describes (make data), or set K3_DATA_DIR"
        )
    missing = [c.stored for c in computation.checks if not (places.data / c.stored).is_file()]
    if missing:
        raise WorkspaceError(f"stored file {missing[0]} not found in {places.data}")


@dataclass
class Execution:
    work_dir: Path
    stdout: Path
    stderr: Path
    exit_code: int
    seconds: float


def prepare(computation: Computation, places: Places) -> Path:
    """A fresh work directory with links to the scripts and copies of the stored files it may read."""
    check_installed(computation, places)
    work_dir = places.work / f"{_slug(computation.id)}-{os.getpid()}"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    for script in sorted(places.scripts.iterdir()):
        if script.is_file() and script.suffix in SCRIPT_SUFFIXES:
            (work_dir / script.name).symlink_to(script.resolve())
    written = computation.written()
    for stored in stored_files(places.data):
        if stored.name not in written and not (work_dir / stored.name).exists():
            shutil.copy2(stored, work_dir / stored.name)
    return work_dir


# Runs `python -I script args` with the fork start method, which was the default on Linux up to Python 3.13.
FORK_LAUNCHER = (
    "import multiprocessing, runpy, sys; multiprocessing.set_start_method('fork'); "
    "sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name='__main__')"
)


def execute(computation: Computation, places: Places, procs: int) -> Execution:
    """Runs the computation's command in a fresh work directory; the stored files are not touched."""
    before = _signature(stored_files(places.data))
    work_dir = prepare(computation, places)
    stdout, stderr = work_dir / "_stdout.txt", work_dir / "_stderr.txt"
    script = [str(work_dir / computation.script), *computation.argv(procs)]
    command = [places.python, "-I", *(["-c", FORK_LAUNCHER] if computation.needs_fork else []), *script]
    started = time.monotonic()
    with stdout.open("wb") as out, stderr.open("wb") as err:
        done = subprocess.run(
            command,
            cwd=work_dir,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=subprocess.STDOUT if computation.merge_stderr else err,
            check=False,
        )
    seconds = time.monotonic() - started
    after = _signature(stored_files(places.data))
    changed = sorted(name for name in before if after.get(name) != before[name])
    if changed:
        raise WorkspaceError(f"stored file {changed[0]} changed while the computation ran ({len(changed)} in all)")
    return Execution(work_dir, stdout, stderr, done.returncode, seconds)


def produced_path(execution: Execution, produced: str) -> Path:
    if produced == STDOUT:
        return execution.stdout
    if produced == STDERR:
        return execution.stderr
    return execution.work_dir / produced


def evaluate(computation: Computation, execution: Execution, places: Places) -> dict:
    """The result line of the runner protocol, from what the command produced."""
    value: dict[str, object] = {}
    predicted: dict[str, object] = {}
    differences: dict[str, object] = {}
    digests: dict[str, str] = {}
    prefixes = (str(execution.work_dir), str(places.data.resolve()), str(places.scripts.resolve()))
    for check in computation.checks:
        outcome = compare(
            produced_path(execution, check.produced),
            places.data / check.stored,
            check.mode,
            ignore_keys=check.ignore_keys,
            strip_prefixes=prefixes,
            unstored_tail=check.unstored_tail,
        )
        value[check.label] = outcome.verdict
        predicted[check.label] = IDENTICAL
        if outcome.produced_sha256:
            digests[check.produced if not check.produced.startswith("@") else check.produced[1:]] = outcome.produced_sha256
        if outcome.verdict != IDENTICAL or outcome.note:
            differences[check.label] = {"mode": check.mode} | outcome.note
    output = execution.stdout.read_text(errors="replace")
    lines = output.splitlines()
    for pattern in computation.expect:
        key = f"expect /{pattern}/"
        value[key] = "found" if any(re.search(pattern, line) for line in lines) else "absent"
        predicted[key] = "found"
    for pattern, count in computation.counts:
        key = f"lines /{pattern}/"
        value[key] = sum(1 for line in lines if re.search(pattern, line))
        predicted[key] = count
    digests.setdefault("stdout", sha256(execution.stdout))

    agrees = value == predicted if predicted else None
    if agrees is False:
        verdict = "differs"
    elif computation.checks:
        verdict = "reproduced"
    elif predicted:
        verdict = "checked"
    else:
        verdict = "ran"
    return {
        "status": "mismatch" if agrees is False else "ok",
        "value": value or None,
        "predicted": predicted or None,
        "agrees": agrees,
        "verdict": verdict,
        "method": computation.command(),
        "note": computation.note,
        "script_seconds": round(execution.seconds, 3),
        "output_lines": len(lines),
        "digests": digests,
        "differences": differences,
        "caveat": computation.caveat or None,
    }


def finish(execution: Execution, result: dict, keep: bool) -> dict:
    """Deletes the work directory of a reproduced computation; otherwise names it in the result."""
    if keep or result["status"] != "ok":
        try:
            result["work_dir"] = str(execution.work_dir.relative_to(REPOSITORY))
        except ValueError:
            result["work_dir"] = str(execution.work_dir)
    else:
        shutil.rmtree(execution.work_dir, ignore_errors=True)
    return result


def tail(path: Path, lines: int = 15) -> str:
    try:
        return "\n".join(path.read_text(errors="replace").splitlines()[-lines:])
    except OSError:
        return ""
