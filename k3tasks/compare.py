"""Comparing what a rerun produced with the stored file it should reproduce.

A comparison has a mode, chosen per stored file in the registry:

    exact            the bytes are the same
    text             the same text once run times are blanked out (see `normalise_text`)
    text-contained   the normalised text of the rerun is a contiguous part of the stored text: the stored
                     log holds the output of several commands, and the rerun is one of them
    records          a file of JSON lines with the same records, in any order, once the keys that hold
                     run times are dropped
    records-subset   every record of the rerun is in the stored file (which holds more: several runs
                     were appended to it, or the rerun covers part of the stored run)
    json             one JSON document, equal once the keys that hold run times are dropped

Each comparison ends in a verdict, "identical" or "differs", with a short account of the first difference.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

MODES = ("exact", "text", "text-contained", "records", "records-subset", "json")
IDENTICAL = "identical"
DIFFERS = "differs"
MISSING = "missing"

# Run times as the scripts print them: "(3 s)", "(8.1s)", "; 11 s)", "12 minutes", "0m20.574s"
_DURATION = re.compile(r"\b\d+(?:\.\d+)?\s?(?:s|sec|secs|seconds|min|mins|minutes|hours)\b")
# a run time inside a JSON record that a script prints into its log: "seconds": 61.3
_JSON_SECONDS = re.compile(r'"(seconds|setup_s)": -?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?')
_SHELL_TIME = re.compile(r"^(?:real|user|sys)\s+\d+m[\d.]+s\s*$")
# the report of /usr/bin/time -v: an indented "Key: value" line
_GNU_TIME = re.compile(r"^\t[A-Z][^:]*(?:\([^)]*\))?[^:]*: .*$")
# PARI's notice when a script sets its stack; the maximum it reports depends on the machine's memory, and
# whether it reaches the log depends on the cypari2 build
_PARI_STACK = re.compile(r"^\s*PARI stack size set to ")
# a flat dict of counts as Python prints it, "{'ok': 296, 'reducible': 272}": its order is the order in
# which worker results arrived, so the items are sorted before lines are compared
_COUNT_DICT = re.compile(r"\{'[^'{}]*': -?\d+(?:, '[^'{}]*': -?\d+)*\}")
_MAX_SHOWN = 160


def _sorted_counts(match: re.Match[str]) -> str:
    return "{" + ", ".join(sorted(match.group(0)[1:-1].split(", "))) + "}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def normalise_text(text: str, strip_prefixes: tuple[str, ...] = ()) -> list[str]:
    """The lines of a log with run times blanked out and the given path prefixes removed.

    The prefixes are the absolute paths of the work directory and the stored data, so that a script
    that prints the path of a file it wrote prints the same line wherever it ran.  Lines written by the
    shell's `time` and by `/usr/bin/time -v` around a stored log are dropped, and so is PARI's notice of
    its stack size.  The items of a printed dict of counts are put in sorted order.
    """
    lines = []
    for line in text.splitlines():
        if _SHELL_TIME.match(line) or _GNU_TIME.match(line) or _PARI_STACK.match(line):
            continue
        for prefix in strip_prefixes:
            if prefix:
                line = line.replace(prefix.rstrip("/") + "/", "")
        line = _JSON_SECONDS.sub(r'"\1": #', _DURATION.sub("#s", line))
        lines.append(_COUNT_DICT.sub(_sorted_counts, line).rstrip())
    while lines and not lines[-1]:
        lines.pop()
    return lines


def _drop_keys(value: Any, keys: frozenset[str]) -> Any:
    if isinstance(value, dict):
        return {k: _drop_keys(v, keys) for k, v in value.items() if k not in keys}
    if isinstance(value, list):
        return [_drop_keys(v, keys) for v in value]
    return value


def _records(path: Path, keys: frozenset[str]) -> tuple[Counter[str], int]:
    """Canonical JSON of each record of a JSON-lines file, and the number of lines that are not JSON."""
    found: Counter[str] = Counter()
    unreadable = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                unreadable += 1
                continue
            found[json.dumps(_drop_keys(record, keys), sort_keys=True)] += 1
    return found, unreadable


def _short(text: str) -> str:
    return text if len(text) <= _MAX_SHOWN else text[: _MAX_SHOWN - 3] + "..."


@dataclass
class Comparison:
    verdict: str
    produced_sha256: str | None = None
    note: dict[str, Any] = field(default_factory=dict)


def compare(
    produced: Path,
    stored: Path,
    mode: str,
    *,
    ignore_keys: tuple[str, ...] = (),
    strip_prefixes: tuple[str, ...] = (),
) -> Comparison:
    """Compares a produced file with a stored one in the given mode."""
    if mode not in MODES:
        raise ValueError(f"unknown comparison mode {mode!r}")
    if not produced.exists():
        return Comparison(MISSING, None, {"why": f"the rerun wrote no {produced.name}"})
    digest = sha256(produced)
    if not stored.exists():
        return Comparison(MISSING, digest, {"why": f"no stored {stored.name}"})
    keys = frozenset(ignore_keys)

    if mode == "exact":
        same = digest == sha256(stored)
        note = {} if same else _first_difference(
            produced.read_text(errors="replace").splitlines(), stored.read_text(errors="replace").splitlines()
        )
        return Comparison(IDENTICAL if same else DIFFERS, digest, note)

    if mode in ("text", "text-contained"):
        mine = normalise_text(produced.read_text(errors="replace"), strip_prefixes)
        theirs = normalise_text(stored.read_text(errors="replace"), strip_prefixes)
        if mode == "text":
            same = mine == theirs
            return Comparison(IDENTICAL if same else DIFFERS, digest, {} if same else _first_difference(mine, theirs))
        start = _find_block(mine, theirs)
        if start is not None:
            return Comparison(IDENTICAL, digest, {"at_line": start + 1})
        return Comparison(DIFFERS, digest, _best_block_difference(mine, theirs))

    if mode == "json":
        mine_doc = _drop_keys(json.loads(produced.read_text(encoding="utf-8")), keys)
        theirs_doc = _drop_keys(json.loads(stored.read_text(encoding="utf-8")), keys)
        same = mine_doc == theirs_doc
        note = {} if same else {"why": "the documents differ", "keys": _differing_keys(mine_doc, theirs_doc)}
        return Comparison(IDENTICAL if same else DIFFERS, digest, note)

    mine_records, bad_mine = _records(produced, keys)
    theirs_records, _ = _records(stored, keys)
    extra = mine_records - theirs_records
    note: dict[str, Any] = {"records": mine_records.total(), "stored_records": theirs_records.total()}
    if bad_mine:
        note["unreadable_lines"] = bad_mine
    if mode == "records":
        missing = theirs_records - mine_records
        same = not extra and not missing and not bad_mine
        if not same:
            note |= {"not_stored": extra.total(), "not_recomputed": missing.total()}
            if extra:
                note["example_not_stored"] = _short(next(iter(extra)))
            if missing:
                note["example_not_recomputed"] = _short(next(iter(missing)))
    else:
        same = not extra and not bad_mine and mine_records.total() > 0
        if extra:
            note |= {"not_stored": extra.total(), "example_not_stored": _short(next(iter(extra)))}
        elif not mine_records:
            note["why"] = "the rerun wrote no records"
    return Comparison(IDENTICAL if same else DIFFERS, digest, note)


def _first_difference(mine: list[str], theirs: list[str]) -> dict[str, Any]:
    for index, (a, b) in enumerate(zip(mine, theirs, strict=False)):
        if a != b:
            return {"line": index + 1, "rerun": _short(a), "stored": _short(b)}
    return {"line": min(len(mine), len(theirs)) + 1, "rerun_lines": len(mine), "stored_lines": len(theirs)}


def _find_block(block: list[str], text: list[str]) -> int | None:
    if not block:
        return 0
    first = block[0]
    for start, line in enumerate(text):
        if line == first and text[start : start + len(block)] == block:
            return start
    return None


def _best_block_difference(block: list[str], text: list[str]) -> dict[str, Any]:
    """Where the rerun's lines agree longest with a stretch of the stored log, and the line after."""
    best_start, best_length = None, -1
    for start, line in enumerate(text):
        if block and line == block[0]:
            length = 0
            while length < len(block) and start + length < len(text) and text[start + length] == block[length]:
                length += 1
            if length > best_length:
                best_start, best_length = start, length
    if best_start is None:
        return {"why": "the first line of the rerun is not in the stored log", "rerun": _short(block[0]) if block else ""}
    stored_line = text[best_start + best_length] if best_start + best_length < len(text) else "<end of log>"
    return {
        "line": best_length + 1,
        "stored_line": best_start + best_length + 1,
        "rerun": _short(block[best_length]),
        "stored": _short(stored_line),
    }


def _differing_keys(a: Any, b: Any) -> list[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        return sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))[:10]
    return []
