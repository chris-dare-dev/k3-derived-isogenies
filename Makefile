# Targets run from the repository root.
#   make test                          the unit tests of the scheduler, the task registry and the claim tables
#                                      (standard library only)
#   make lint                          ruff and mypy (neither is needed to run anything)
#   make venv                          numpy and cypari2 in a virtual environment under ~/.cache
#   make data                          check the data archive in dist/ and unpack it into experiments/ (DATA.md)
#   make plan BATCH=batches/smoke.toml list what a run would compute, with limits, without computing
#   make run BATCH=batches/smoke.toml  compute the tasks of a batch that have no record yet
#   make report BATCH=batches/smoke.toml REPORT_ARGS=--full
#   make smoke                         the smoke batch, computed afresh, and its report
#   make record BATCH=batches/smoke.toml   run, then copy the records into records/ (committed code only);
#                                      EXPECT=1 for the negative control, EXPECT=3 for the limits demo.
#                                      Commit what a record wrote before the next one: a file that git does
#                                      not know makes the working tree differ from its commit.  To record
#                                      several batches in one go, `make run` each and then `make adopt`
#   make adopt BATCHES="smoke audits"  copy the records of finished runs into records/
#   make tables                        write the generated tables of CLAIMS.md, GUIDE.md and claims/RECORDS.md
#   make check                         is the registry consistent, are the data installed, and are those
#                                      tables what the claim registry and the records give?
#   make list                          every computation of the registry
# ARGS is passed to the scheduler, for example ARGS="--workers 2 --memory-budget 8000".
# The interpreter of the stored results and of the records is Python 3.12 (.python-version):
# `make venv PYTHON=python3.12`.
PYTHON  ?= python3
VENV    ?= $(HOME)/.cache/k3-derived-isogenies/venv
# the interpreter that runs the scripts of experiments/ (needs numpy and cypari2)
K3_PYTHON ?= $(if $(wildcard $(VENV)/bin/python),$(VENV)/bin/python,$(PYTHON))
export K3_PYTHON
BATCH   ?= batches/smoke.toml
ARGS    ?=
REPORT_ARGS ?=
ARCHIVE ?= dist/k3-derived-isogenies-data-v1.tar.xz
SHA256  := 7a16e7d759533558269147fef4b931c3cdaae66eae78fd41463cf2c9977b3714

RUN_DIR = results/$(basename $(notdir $(BATCH)))
# the status `run` must end with for `record` to adopt: 0, or 1 when every task is a negative control,
# or 3 when tasks are stopped at their limits on purpose
EXPECT  ?= 0
CHECKED   = taskbatch k3tasks runners tests tools
FORMATTED = taskbatch k3tasks runners tests tools

BATCHES ?= $(basename $(notdir $(BATCH)))

.PHONY: test lint format check list plan run report smoke record adopt tables venv data clean-work

test:
	$(PYTHON) -m unittest discover -s tests -t .

lint:
	ruff check $(CHECKED)
	ruff format --check $(FORMATTED)
	mypy

format:
	ruff format $(FORMATTED)

check:
	$(PYTHON) -m k3tasks check
	$(PYTHON) tools/claims.py check

list:
	$(PYTHON) -m k3tasks list

plan:
	$(PYTHON) -m taskbatch plan $(BATCH) $(ARGS)

run:
	$(PYTHON) -m taskbatch run $(BATCH) $(ARGS)

report:
	$(PYTHON) -m k3tasks report $(RUN_DIR) $(REPORT_ARGS)

smoke:
	rm -rf results/smoke
	$(PYTHON) -m taskbatch run batches/smoke.toml --quiet $(ARGS)
	$(PYTHON) -m k3tasks report results/smoke

record:
	status=0; $(PYTHON) -m taskbatch run $(BATCH) $(ARGS) || status=$$?; \
	if [ $$status -ne $(EXPECT) ]; then echo "run ended with status $$status, expected $(EXPECT)"; exit 1; fi
	$(PYTHON) tools/claims.py adopt $(RUN_DIR)

adopt:
	$(PYTHON) tools/claims.py adopt $(addprefix results/,$(BATCHES))

tables:
	$(PYTHON) tools/claims.py tables

venv:
	$(PYTHON) -m venv $(VENV)
	$(VENV)/bin/pip install -r experiments/requirements.txt

data:
	echo "$(SHA256)  $(ARCHIVE)" | sha256sum -c
	tar --skip-old-files -xJf $(ARCHIVE)

# work directories that were kept because a rerun differed or was stopped
clean-work:
	rm -rf results/work
