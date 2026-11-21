# Guide: derived isogenies of ordinary K3 surfaces

The full reference for [the README](README.md): results table, limits, layout, reproduction and batches.

*A negative result for isogeny-style cryptography, with notes, code and data.*

## The question

Supersingular elliptic curves with a class-group action underlie CSIDH. Ordinary K3 surfaces over a finite field carry a similar structure: those with the same Frobenius fall into torsors under a class group of a CM field, and the group acts by twisted derived equivalences ("derived isogenies"). Does this give a **new hard problem** for cryptography?

## The answer

**No.** In the one family where everything can be computed, double planes branched along six lines (Picard number 16), the problem reduces to a CSIDH-type class-group action on abelian fourfolds with multiplication by Z[i], and it is worse for cryptography than CSIDH: the group has about 0.01·q² elements, so a 256-bit group needs about 525-bit instances (512 for CSIDH-512), steps cost about ℓ⁴ instead of ℓ, and no honest party can compute a step at a useful size. The decision problem leaks one bit, a quadratic character, for 80 of the 477 examined polynomials that have a quadratic character. Below Picard number 16 there is no known reduction, but also no computable step, so there is nothing to build on. All of this is measured on small instances; at cryptographic size it is argument.

The full account is [SUMMARY.md](SUMMARY.md), a self-contained summary paper of about 12 pages. The detailed notes are in [docs/](docs/README.md).

## Main results

<!-- BEGIN GENERATED main: tools/claims.py tables -->
| Result | Numbers | Where | Evidence |
|---|---|---|---|
| Surfaces form torsors under a CM class group, acting by twisted derived equivalences (ρ̄ ≥ 12, p ≥ 5) | Assembly of Taelman, Valloni, Bayer-Fluckiger, Bragg–Yang | notes [1](docs/01-prior-art.md), [6](docs/06-torsors-under-cm-class-groups.md) | [P], [D] |
| Six-line surfaces match Hermitian lattices | 44 of 44 polynomials (p ≤ 37), 169 of 169 non-split, 32 of 32 for the step at p | notes [2](docs/02-six-line-count.md)–[5](docs/05-non-split-count.md) | [E] |
| **Reduction to abelian fourfolds.** P is a second exterior power of the Frobenius quartic of a fourfold with Z[i]-multiplication (Proposition Y); the fourfold is built over F_p as a Prym variety of a genus-5 curve (Paranjape) | 278 surfaces (534 configurations), 4,983 pairs of lines, no exception; also over F₂₅, F₄₉ | notes [13](docs/13-hardness-on-paper.md), [14](docs/14-fourfold-of-a-six-line-surface.md), [15](docs/15-fourfold-map-and-the-group-action.md), [16](docs/16-fourfold-over-f25-and-f49.md) | [D], [E]; equivariance a sketch |
| **The group** G = Cl(E)/⟨Cl(E₀), ramified primes⟩ | \|G\| ≈ 0.01·q² over 512 groups to q = 100,003; 2²⁵⁶ needs q ≈ 2¹³¹ | note [13](docs/13-hardness-on-paper.md) | [E] in range, [H] beyond |
| **The CM method computes the action on small groups** | 98 of 98 surfaces in 17 fields (150 of 150 in 21 fields with notes 26 and 29); ℓ-steps equal Frobenius in characteristic ℓ, 346 of 346; about 190 digits of precision per group element | notes [20](docs/20-cm-method-and-the-action.md), [24](docs/24-the-form-mod-4-is-the-character.md), [27](docs/27-hypothesis-ledger.md) | [E] |
| **One bit leaks.** A sign δ on the fourfold is the character χ₋₁: the genus-character attack of Castryck–Sotáková–Vercauteren carries over (it affects class-group actions with 2-torsion, not CSIDH itself) | 60 of 60 surfaces in 5 torsors at 3 primes; χ₋₁ exists for 80 of the 477 examined polynomials that have a quadratic character | notes [19](docs/19-quadratic-characters-from-a-surface.md), [23](docs/23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](docs/24-the-form-mod-4-is-the-character.md), [25](docs/25-checks-of-notes-20-23-24.md), [26](docs/26-delta-is-the-character-at-q-41.md) | [E], [D] |
| **Nothing else leaks** in a pre-registered audit | about 11,400 candidate invariants; only δ survives correction | note [28](docs/28-leakage-audit.md) | [E] |
| **No character of order 3 is read either**, in a second audit | 7 torsors with 3 ∣ \|G\| (72 surfaces, q = 23 to 41): no cubic residue symbol is a character on the two torsors where one exists, and no candidate separates the cosets beyond chance on any of the seven. At q = 37 the cube roots of unity are cubes, so a reading R with R³ = c³·d would have a constant symbol there and be invisible to these tests. Three more torsors at q = 43 (48 surfaces), where they are not cubes, were scored after the run: no cubic symbol is a character there either, and their one candidate past Holm does not replicate on a held-out torsor (note 29 §7.4) | note [29](docs/29-characters-of-order-3.md) | [E], [D] |
| **Below Picard 16**: no reduction known and no computable step | a reading of the literature; half a step at ℓ = 2 tested on 21,240 elliptic K3 pairs | notes [17](docs/17-below-picard-16.md), [22](docs/22-two-neighbours-of-elliptic-k3.md) | [P], [M], [E] |
| **Ledger**: every headline number recomputed from the stored files | 162 of 170 match the numbers as first written; the three slips behind the 8 others are corrected in the notes | note [27](docs/27-hypothesis-ledger.md) | [E] |
<!-- END GENERATED main -->

Evidence labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived in this work; **[E]** measured by the scripts in this repository; **[H]** heuristic; **[R]** re-derived in a separate re-derivation check, with its own code.

## Limits

- **Small instances only.** All surface data are over fields with q ≤ 49; the torsors have at most 36 surfaces (at q = 49), and the action was evaluated on at most 16. A cryptographic instance would need q ≈ 2¹³¹ and |G| ≈ 2²⁵⁶.
- **Argument, not measurement, at cryptographic size.** There the verdict rests on Proposition Y (valid for all q), on the equivariance of the surface-to-fourfold map (a sketch, never tested at a step ℓ ≠ p on the fourfold), on the CM points being canonical lifts (a sketch), on recalled costs of attacks, and on extrapolating |G| ≈ 0.01·q² ([note 27](docs/27-hypothesis-ledger.md) §1).
- **One family.** Only six-line double planes are treated explicitly.
- **No outside review.** The re-derivation checks (notes 18 and 25, and the check whose findings notes 26 to 29 state) and the ledger are part of this work; nothing here has been refereed.

## Repository layout

| Path | What it is |
|---|---|
| [SUMMARY.md](SUMMARY.md) | Self-contained summary paper |
| [docs/](docs/README.md) | 29 notes with the arguments, computations and checks, and an index |
| [experiments/](experiments/) | The Python scripts behind the notes, and `requirements.txt` |
| [CLAIMS.md](CLAIMS.md) | Every row of the main results and every hypothesis of the ledger, with the stored lines that carry its numbers, the tasks that bear on it and what their records say |
| [claims/](claims/) | `claims.toml`, the registry of claims from which the tables of CLAIMS.md and the table of main results above are written, and `RECORDS.md`, every task with its record |
| [records/](records/) | The records of the batches that were run under the scheduler, each from committed code |
| [batches/](batches/) | Batch files: which computations of the notes to rerun, with which time and memory limits |
| [k3tasks/](k3tasks/) | The registry of computations (one per command of a note's *Reproduce* section), the work directories they run in, and the comparison of a rerun with the stored files |
| [taskbatch/](taskbatch/) | The batch scheduler: one process per task, time and memory limits, one record per task, reports ([PROTOCOL.md](taskbatch/PROTOCOL.md)) |
| [runners/](runners/) | `task_runner.py`, the command through which the scheduler expands batches and runs tasks |
| [tools/](tools/) | `claims.py`: writes the generated tables and checks them; `adopt.py`: copies the records of a finished run into `records/` |
| [tests/](tests/) | Unit tests of the scheduler, the registry, the batch files and the claim tables (standard library only) |
| `Makefile`, `pyproject.toml` | `make test`, `make lint`, `make plan`, `make run`, `make report`, `make record`, `make tables`, `make check`; the settings of the scheduler, ruff and mypy |
| `Dockerfile`, `.python-version` | Python 3.12 with the packages of `experiments/requirements.txt`, for machines without them |
| [notebooks/](notebooks/README.md) | Two notebooks for readers without a command line: explore the stored results, run small instances live |
| `binder/` | Configuration for running the notebooks on a free Binder server |
| [DATA.md](DATA.md) | The stored results: what they are and how to install them |
| [LICENSE](LICENSE), [LICENSE-TEXT.md](LICENSE-TEXT.md) | MIT for code; CC BY 4.0 for the prose and data |
| [CITATION.cff](CITATION.cff) | How to cite |

## Reproducing

**Without a command line**: see [notebooks/](notebooks/README.md).

**From a terminal.** Requirements: Python 3 (tested with 3.12), `numpy` and `cypari2` (a binary wheel on Linux; see `experiments/requirements.txt` for the versions used).

```bash
python -m venv .venv
.venv/bin/pip install -r experiments/requirements.txt
```

Install the stored results as described in [DATA.md](DATA.md) (they unpack into `experiments/`). Then, from `experiments/`, with `python -I` as in the notes:

```bash
cd experiments
../.venv/bin/python -I sixlines_ledger.py > ledger-rerun.log && diff ledger-rerun.log sixlines_ledger.log
../.venv/bin/python -I sixlines_cmmethod.py selftest --digits 60
../.venv/bin/python -I sixlines_prym.py --check 17 2 3 4 11
```

The first recomputes the 170 comparisons of note 27 from the stored files (about a second) and should reproduce the stored log exactly. The second checks the theta-constant formulas of note 20 without CM (a second). The third builds the fourfold of one surface over F₁₇ for all 30 ordered pairs of lines (a few seconds). Every note that computes something ends with a *Reproduce* section giving the exact commands for its results; some take hours, and several overwrite or append to the stored files, so rerun them in a copy.

What was verified on a fresh clone with the archive unpacked is listed in [DATA.md](DATA.md#verification).

## Running a batch

The commands of the notes' *Reproduce* sections are also registered as tasks, so that they can be rerun in batches under time and memory limits, each compared with the stored file it should reproduce. One task is one command with one answer: for example the ledger, or the enumeration of six lines over F₃₇, or the coset run at q = 41 and e₁ = 22. `python3 -m k3tasks list` lists all 163, with the time each took when its note was written.

```bash
make venv
```

```bash
make data
```

```bash
make plan BATCH=batches/smoke.toml
```

```bash
make run BATCH=batches/smoke.toml
```

```bash
make report BATCH=batches/smoke.toml
```

`make venv` installs numpy and cypari2 in `~/.cache/k3-derived-isogenies/venv`, which runs the scripts (`K3_PYTHON` names another interpreter). `make data` checks the archive in `dist/` against its SHA-256 and unpacks it into `experiments/`, adding files only (see [DATA.md](DATA.md)). `plan` lists the tasks a run would compute, with their expected times and limits, and computes nothing. `run` computes every task that has no record yet; a run can be stopped with Ctrl-C and continued later. `report` prints the scheduler's tables and the reproduction table. Through make it also warns if the run names this machine's name or a directory outside the repository, which `make adopt` would refuse (`k3tasks report RUN --private-words FILE` takes the words from a file, `-` for standard input). The same commands without make are `python3 -m taskbatch plan|run batches/smoke.toml` and `python3 -m k3tasks report results/smoke`. The scheduler and the registry need the standard library only; `make test` runs their unit tests.

### What a task does

A task never runs in `experiments/` and never writes a stored file. It runs in a work directory of its own, `results/work/<task>-<pid>/`, which holds links to the scripts and copies of the stored results, except the files the task writes, which start absent as on a fresh run. The scripts take the directory they are run from as their own, so they read the copies and write next to them. The sizes and times of the stored files are checked before and after, and a change is an error. The command is the note's own, `python -I script …`. Two scripts, `sixlines_fq.py` and `sixlines_nonsplit.py`, hand tables to their worker processes by fork, the default start method on Linux up to Python 3.13; from Python 3.14 the default is forkserver, under which they fail, so the harness starts them with the fork method set and otherwise unchanged.

When the command ends, its output and the files it wrote are compared with the stored files named in the registry, in one of these modes:

| Mode | Identical when |
|---|---|
| `exact` | the bytes are the same (the ledger, the CM self-test) |
| `text` | the text is the same once run times (also those inside printed JSON records), the shell's `time` report and PARI's stack-size notice are blanked out, and the items of each printed dict of counts are sorted (worker results arrive in any order) |
| `text-contained` | the output is a contiguous part of a stored log that holds the output of several commands |
| `records` | the JSON records are the same in any order, without the keys that hold run times (`seconds`, `setup_s`) |
| `records-subset` | every record of the rerun is in the stored file, which holds more (appended runs) |
| `json` | the JSON document is the same, without run times |

Each task ends with a verdict: **reproduced** (identical to every stored file it was compared with), **differs** (the report says what, and the work directory is kept), **checked** (nothing stored to compare with, but its output has the lines it must have, such as `selftest passed`), **ran** (its output was not kept, so nothing is compared; the registry says why), or the scheduler's status when it did not finish: `timeout`, `memory`, `error`, `crashed`. Each record in `results/<batch>/results.jsonl` also has the wall and processor time, the peak memory (of the largest process, and under the cgroup guard of all the task's processes together, which the reproduction table shows), the SHA-256 of the output and of every file compared, the command line, and the versions of Python, numpy, PARI and the stored data.

### Limits

Each task has a time limit and a memory limit, and a batch has a number of workers and a memory budget: the memory limits of the tasks running at the same time never add up to more than the budget, and `run` refuses a budget above the memory available now. Where `systemd-run --user` works, each task runs in its own systemd scope with `MemoryMax` and no swap, so the kernel stops it at its limit together with all the processes it started; elsewhere the scheduler polls the task's resident memory. The scheduler stops a task at its time limit, and no task outlives the scheduler. Several scripts start their own processes (`--procs`, `--workers`): a batch sets how many with `procs`, and the expected time grows in proportion. Results and work directories are kept under `results/`, never in `/tmp`. Check `free -h` before a long batch.

### The batch file

```toml
name = "example"                 # results go to results/example

[scheduler]
workers = 2                      # tasks at a time
memory_budget_mb = 8000          # the memory limits of running tasks never add up to more

[limits]                         # limits of one task
time_seconds = 600
memory_mb = 1500

[[limits.tier]]                  # tasks expected to take 400 s or more get more time
min_cost = 400
time_seconds = 3600

[[tasks]]
kind = "computation"
groups = ["lattices"]            # also: ids (or patterns such as "n02/enumerate-*"), notes, tags
exclude_tags = ["heavy"]         # also: exclude (ids or patterns)
procs = 6                        # processes a task may start, at most
limits = { memory_mb = 4000 }    # this table's own limits
```

A task's cost is the time its command took when its note was written, scaled to the processes the batch allows; tasks run cheapest first, and tiers give the long ones more time. A task gets the limits of its table if there are any, otherwise those of the highest tier whose `min_cost` it reaches, otherwise the defaults. On the command line `--workers`, `--memory-budget`, `--time-limit` and `--memory-limit` replace the batch's settings. `taskbatch/PROTOCOL.md` is the contract between the scheduler and the runner.

| File | Content | Time |
|---|---|---|
| `smoke.toml` | the ledger (byte for byte), the CM self-test, the fourfold over F₁₇, the self-tests of the enumeration scripts, and small reruns of notes 8, 25 and 29 | about 10 seconds |
| `lattices.toml` | notes 2 to 5: enumerations over F_p (p ≤ 37) and F₂₅, lattice counts, neighbour steps, forms over the prime field, non-split counts | about 30 minutes |
| `ring-steps.toml` | notes 6 to 12: the lattice graphs of ring-changing steps | 9 minutes |
| `fourfold.toml` | notes 13 to 16, 18 and 21: the group sizes, the fourfold over F_p and F₂₅, the dictionary, the checks of note 18, forms over F_p | 32 minutes |
| `cm-method.toml` | notes 20, 24 and 25: the CM method at one field, the period points, the delta labels, the Galois checks | 19 minutes |
| `quadratic-character.toml` | notes 19, 23, 25 and 26: candidate characters, the Jacobian checks, the steps at q = 41 that read the stored enumeration | 18 minutes |
| `audits.toml` | notes 28 and 29: the features and the pre-registered audits, the q = 43 audit and its held-out test | 8 minutes |
| `below-16.toml` | note 22: the neighbour graphs and the two short controls | half a minute |
| `negative-control.toml` | three tasks whose predictions are wrong on purpose: each must end `mismatch`, and the run ends with status 1 | seconds |
| `limits-demo.toml` | the ledger next to a task that sleeps past its time limit and one that allocates past its memory limit: `timeout`, `memory`, and the run ends with status 3 | seconds |
| `heavy.toml` | everything that takes from a quarter of an hour to several hours (the enumerations over F₄₉ and at q = 41 and 43, the full Prym runs, the long CM sweeps, the delta and coset runs at q = 41 and 43, the elliptic K3 runs), and the audits again: not run by default | 12 h 22 min (the one run, two workers) |

The times are the wall times of one run on a 16-core machine, with each batch's own workers and processes. `heavy.toml` was run once as a batch, in 12 h 22 min with its two workers; the times its notes record add up to 19 hours. Of its 47 tasks, 29 are ok; the 18 that are not are listed in [CLAIMS.md](CLAIMS.md#known-differences-between-reruns-and-stored-files).

Every computation of the registry is in at least one batch, and the unit tests check that every task's limits are above its expected time and memory.

### Records and claims

`make run` writes into `results/`, which is not in git. A finished run becomes part of the repository through

```bash
make record BATCH=batches/smoke.toml
```

which runs the batch and copies its records into `records/smoke/`. Records are adopted only if they were computed from committed code: each names the commit of the code that computed it, and a record from a working tree with uncommitted changes is refused, as is one that names the machine it ran on. A file that git does not yet know counts as a change, so the records of one batch are committed before the next batch is recorded; or several batches are run with `make run` and adopted together with `make adopt BATCHES="smoke audits"`. Next to the records, `adopted.json` gives the tree of each commit they name. Nothing in `records/` is edited by hand.

[CLAIMS.md](CLAIMS.md) connects the records to the results: for each row of the table of main results and each hypothesis of the ledger, the stored lines that carry its numbers, the tasks that bear on it, and what their records say. A task without a record reads "not re-run". Its tables, [claims/RECORDS.md](claims/RECORDS.md) and the table of main results above are written by `make tables` from [claims/claims.toml](claims/claims.toml), the batch files and the records; `make check` fails if they are not what those give, or if a stored line named as evidence is not where the claim registry says.

The records in this repository were computed on Python 3.12 with the packages of `experiments/requirements.txt`, the interpreter of the stored results ([DATA.md](DATA.md#rerunning-against-the-stored-files)). `make lint` runs ruff and mypy; the `Dockerfile` gives an image with all of it:

```bash
docker build -t k3-derived-isogenies .
```

## Licence and citation

Code: [MIT](LICENSE). Text (the notes, the summary and this guide) and data: [CC BY 4.0](LICENSE-TEXT.md). Citation metadata: [CITATION.cff](CITATION.cff). The papers the notes cite are not included.

## A note on method

The idea for this project and its direction are my own. AI tools (Anthropic's Claude) were used extensively in carrying it out, including the derivations, code, computations and drafting. The results have not been independently reviewed. If you would like more detail on how the work was done, get in touch.
