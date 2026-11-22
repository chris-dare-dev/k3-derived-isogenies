#!/usr/bin/env python3
"""The entry point the scheduler calls (src/taskbatch/PROTOCOL.md).  One process expands a batch or runs one
task.

    python3 -I src/runners/task_runner.py expand '<batch description as JSON>'
    python3 -I src/runners/task_runner.py run '<task as JSON>'

"expand" prints one line "@task <JSON>" per task and a last line "@meta <JSON>" with the number of tasks
and the versions in use.  "run" runs the task's script in a work directory of its own and prints one line
"@result <JSON>".  An error ends the process with status 1 and the message on the standard error.

The scripts run under $K3_PYTHON (default: the interpreter running this file), which needs numpy and
cypari2; this file itself needs the standard library only.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from k3tasks import __version__
from k3tasks.tasks import expand_batch, run_task
from k3tasks.workspace import DATA_MARKER, Places, WorkspaceError


def emit(mark: str, content: object) -> None:
    print(mark, json.dumps(content, ensure_ascii=False), flush=True)


def versions(places: Places) -> dict[str, object]:
    """What a record should say about the software and the stored data it was compared with."""
    out: dict[str, object] = {"k3tasks": __version__, "runner_python": platform.python_version()}
    try:
        probe = subprocess.run(
            [places.python, "-I", "-c",
             "import platform, numpy, cypari2; print(platform.python_version(), numpy.__version__, "
             "cypari2.Pari()('version()'))"],
            capture_output=True, text=True, timeout=60, check=False,
        )  # fmt: skip
        if probe.returncode == 0:
            python, numpy, pari = probe.stdout.split(" ", 2)
            out |= {"python": python, "numpy": numpy, "pari": pari.strip()}
        else:
            out["scripts_python"] = "numpy or cypari2 missing"
    except (OSError, subprocess.TimeoutExpired):
        out["scripts_python"] = "cannot start"
    marker = places.data / DATA_MARKER
    out["data"] = hashlib.sha256(marker.read_bytes()).hexdigest()[:12] if marker.exists() else "not installed"
    return out


def scope_peak_mb() -> float | None:
    """Peak memory of the whole task, all its processes together, when it runs in a systemd scope of its own.

    The scheduler's `peak_rss_mb` is the largest resident size of one process; a script with workers uses
    more.  Under the cgroup guard the task's scope records the peak of everything in it.
    """
    try:
        line = Path("/proc/self/cgroup").read_text().splitlines()[-1]
        path = line.split("::", 1)[1]
        if "taskbatch-" not in path or not path.endswith(".scope"):
            return None
        return round(int(Path("/sys/fs/cgroup" + path, "memory.peak").read_text()) / 2**20, 1)
    except (OSError, IndexError, ValueError):
        return None


def main(arguments: list[str]) -> int:
    if len(arguments) != 2 or arguments[0] not in {"expand", "run"}:
        print("usage: task_runner.py (expand|run) '<JSON>'", file=sys.stderr)
        return 2
    try:
        request = json.loads(arguments[1])
        if arguments[0] == "run":
            result = run_task(request)
            peak = scope_peak_mb()
            if peak is not None:
                result["scope_peak_mb"] = peak
            emit("@result", result)
        else:
            tasks = expand_batch(request)
            for task in tasks:
                emit("@task", task)
            emit("@meta", {"tasks": len(tasks)} | versions(Places.from_environment()))
    except (WorkspaceError, ValueError, KeyError, OSError) as problem:
        message = str(problem) if not isinstance(problem, KeyError) else f"missing key {problem}"
        print(f"task_runner.py: error: {message}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
