"""What a rerun must give, what it gave, and the verdict: the words of a task's result.

A result compares up to three kinds of thing.  Each has a key of its own, the same in `value` (what the
rerun gave) and in `predicted` (what it must give):

    <name of a stored file>   how the comparison with that file ended; "identical" is predicted
    expect /<expression>/     "found" or "absent": whether a line of the output matches; "found" is predicted
    lines /<expression>/      the number of lines of the output that match

and the task ends with one verdict:

    reproduced   everything it was compared with is identical to the stored files
    differs      something differs; the report says what, and the work directory was kept
    checked      no stored file to compare with, but the lines its output must contain are there
    ran          the command ran to the end and nothing was compared (its output was not kept)
    timeout, memory, error, crashed, unsupported   the scheduler's status: no verdict
    not run      no record yet

These keys and words are written into the committed records, so none of them may change.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from taskbatch.model import Status

from .compare import IDENTICAL

FOUND = "found"
ABSENT = "absent"

REPRODUCED = "reproduced"
DIFFERS = "differs"
CHECKED = "checked"
RAN = "ran"
NOT_RUN = "not run"

# every word a task can stand under in a report, in the order of its summary
_NO_VERDICT = (Status.TIMEOUT, Status.MEMORY, Status.ERROR, Status.CRASHED, Status.UNSUPPORTED)
ORDER = [REPRODUCED, DIFFERS, CHECKED, RAN, *(str(status) for status in _NO_VERDICT), NOT_RUN]

# the verdicts of a task that ended well, each as the claim tables word it
ENDED_WELL = {
    REPRODUCED: "reproduced",
    CHECKED: "checked by the lines of their output",
    RAN: "ran with nothing to compare",
}

_EXPECT = "expect "
_LINES = "lines "


def _expect_key(expression: str) -> str:
    return f"{_EXPECT}/{expression}/"


def _count_key(expression: str) -> str:
    return f"{_LINES}/{expression}/"


def of_the_output(key: str) -> bool:
    """Whether a key of `value` or `predicted` is about the lines of the output, not a stored file."""
    return key.startswith((_EXPECT, _LINES))


def predicted(
    files: Iterable[str], expect: Iterable[str], counts: Iterable[tuple[str, int]]
) -> dict[str, Any]:
    """What a rerun must give: each stored file identical, each expected line there, each count met."""
    expected: dict[str, Any] = dict.fromkeys(files, IDENTICAL)
    expected |= {_expect_key(expression): FOUND for expression in expect}
    expected |= {_count_key(expression): count for expression, count in counts}
    return expected


def observed(
    lines: Sequence[str], expect: Iterable[str], counts: Iterable[tuple[str, int]]
) -> dict[str, Any]:
    """What the lines of an output give for the expressions that `predicted` names."""
    found: dict[str, Any] = {}
    for expression in expect:
        found[_expect_key(expression)] = (
            FOUND if any(re.search(expression, line) for line in lines) else ABSENT
        )
    for expression, _ in counts:
        found[_count_key(expression)] = sum(1 for line in lines if re.search(expression, line))
    return found


@dataclass(frozen=True)
class Judgement:
    status: str  # the status the runner reports: ok or mismatch
    agrees: bool | None  # None when nothing was predicted
    verdict: str


def judge(value: Mapping[str, Any], expected: Mapping[str, Any], *, files_compared: bool) -> Judgement:
    """The outcome of a command that ran to the end, from what it gave and what it must give."""
    agrees = value == expected if expected else None
    if agrees is False:
        return Judgement(str(Status.MISMATCH), False, DIFFERS)
    if files_compared:
        return Judgement(str(Status.OK), agrees, REPRODUCED)
    return Judgement(str(Status.OK), agrees, CHECKED if expected else RAN)


def verdict(record: Mapping[str, Any] | None) -> str:
    """The word a task stands under: its verdict if its command ran to the end, else its status."""
    if record is None:
        return NOT_RUN
    if record["status"] in (Status.OK, Status.MISMATCH):
        return str(
            record.get("detail", {}).get("verdict")
            or (DIFFERS if record["status"] == Status.MISMATCH else RAN)
        )
    return str(record["status"])


def compared_with(record: Mapping[str, Any] | None) -> list[str]:
    """The stored files a record was compared with."""
    value = (record["value"] if record else None) or {}
    return [key for key in value if not of_the_output(key)]


def differing(record: Mapping[str, Any]) -> dict[str, Any]:
    """The part of what a record gave that is not what was predicted."""
    expected = record.get("predicted") or {}
    return {key: value for key, value in (record.get("value") or {}).items() if value != expected.get(key)}


def differences(record: Mapping[str, Any]) -> list[str]:
    """What differs in a record, one account per stored file or expression."""
    value = record.get("value") or {}
    expected = record.get("predicted") or {}
    parts = []
    for label, note in (record.get("detail", {}).get("differences") or {}).items():
        if value.get(label) == IDENTICAL:
            continue
        shown = ", ".join(f"{k} {v}" for k, v in note.items() if k != "mode")
        parts.append(f"{label} ({note.get('mode')}): {shown}")
    for key, found in value.items():
        if of_the_output(key) and found != expected.get(key):
            parts.append(f"{key}: {found}, expected {expected.get(key)}")
        elif found == IDENTICAL and expected.get(key) != found:
            parts.append(f"{key}: identical, predicted {expected.get(key)}")
    return parts
