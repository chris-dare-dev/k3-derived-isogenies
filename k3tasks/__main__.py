"""Command line of the task registry.

    python3 -m k3tasks list [--batch FILE] [--group G] [--note N]   the computations, with times and limits
    python3 -m k3tasks report RUN_DIR [--full]                     the scheduler's report and the reproduction table
    python3 -m k3tasks check                                       is the registry consistent, are the data installed?

Planning and running go through the scheduler: `python3 -m taskbatch plan|run batches/<name>.toml`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from taskbatch.config import ConfigError, load_config

from . import registry, report
from .tasks import expand_batch
from .workspace import DATA_MARKER, Places, stored_files


def _list(arguments: argparse.Namespace) -> int:
    if arguments.batch:
        config = load_config(arguments.batch)
        tasks = expand_batch(config.description)
        rows = []
        for task in tasks:
            limits = config.limits_for(_as_task(task))
            rows.append((task, limits))
        print(f"{'id':38} {'group':20} {'procs':>5} {'expect s':>9} {'limit s':>8} {'expect MB':>9} {'limit MB':>8}")
        for task, limits in rows:
            print(f"{task['id']:38} {task['group']:20} {task['procs']:5d} {task['cost']:9g} "
                  f"{limits.time_seconds:8g} {task['expected_memory_mb']:9d} {limits.memory_mb:8d}")
        return 0
    chosen = registry.select(
        groups=[arguments.group] if arguments.group else None, notes=[arguments.note] if arguments.note else None
    )
    print(f"{'id':38} {'group':20} {'procs':>5} {'expect s':>9} {'expect MB':>9}  checks")
    for c in chosen:
        checks = ", ".join(f"{k.stored} ({k.mode})" for k in c.checks) or ("output lines" if c.expect or c.counts else "-")
        tag = " [heavy]" if "heavy" in c.tags else ""
        print(f"{c.id:38} {c.group:20} {c.procs:5d} {c.seconds:9g} {c.expected_memory_mb(c.procs):9d}  {checks}{tag}")
    print(f"\n{len(chosen)} computations")
    return 0


def _as_task(payload: dict) -> object:
    from taskbatch.model import Task

    return Task.from_payload(payload)


def _report(arguments: argparse.Namespace) -> int:
    if not (arguments.run_dir / "results.jsonl").exists():
        print(f"k3tasks: no results.jsonl in {arguments.run_dir}", file=sys.stderr)
        return 2
    sys.stdout.write(report.render(arguments.run_dir, full=arguments.full))
    return 0


def _check(arguments: argparse.Namespace) -> int:
    problems = registry.validate()
    places = Places.from_environment()
    for problem in problems:
        print(f"registry: {problem}")
    print(f"{len(registry.all_computations())} computations in {len(registry.GROUPS)} groups")
    if (places.data / DATA_MARKER).exists():
        stored = {p.name for p in stored_files(places.data)}
        missing = sorted({c.stored for c in registry.all_computations() for c in c.checks} - stored)
        print(f"stored results: {len(stored)} files in {places.data}")
        for name in missing:
            print(f"stored file named in the registry but not installed: {name}")
            problems.append(name)
    else:
        print(f"stored results: not installed in {places.data} (see DATA.md, or make data)")
    print(f"scripts run with: {places.python}")
    return 1 if problems else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="k3tasks", description="The computations behind the notes, as tasks.")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="list the computations, or the tasks of a batch with their limits")
    listing.add_argument("--batch", type=Path)
    listing.add_argument("--group", choices=sorted(registry.GROUPS))
    listing.add_argument("--note", type=int)
    rep = commands.add_parser("report", help="the report of a run directory, with the reproduction table")
    rep.add_argument("run_dir", type=Path)
    rep.add_argument("--full", action="store_true")
    commands.add_parser("check", help="check the registry and the installed data")
    arguments = parser.parse_args(argv)
    try:
        return {"list": _list, "report": _report, "check": _check}[arguments.command](arguments)
    except (ConfigError, ValueError, json.JSONDecodeError) as problem:
        print(f"k3tasks: {problem}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
