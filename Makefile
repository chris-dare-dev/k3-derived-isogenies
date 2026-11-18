# Targets run from the repository root.
#   make test                          the unit tests of the scheduler and the task registry (standard library only)
#   make venv                          numpy and cypari2 in a virtual environment under ~/.cache
#   make data                          check the data archive in dist/ and unpack it into experiments/ (DATA.md)
#   make plan BATCH=batches/smoke.toml list what a run would compute, with limits, without computing
#   make run BATCH=batches/smoke.toml  compute the tasks of a batch that have no record yet
#   make report BATCH=batches/smoke.toml REPORT_ARGS=--full
#   make list                          every computation of the registry
# ARGS is passed to the scheduler, for example ARGS="--workers 2 --memory-budget 8000".
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

.PHONY: test check list plan run report venv data clean-work

test:
	$(PYTHON) -m unittest discover -s tests -t .

check:
	$(PYTHON) -m k3tasks check

list:
	$(PYTHON) -m k3tasks list

plan:
	$(PYTHON) -m taskbatch plan $(BATCH) $(ARGS)

run:
	$(PYTHON) -m taskbatch run $(BATCH) $(ARGS)

report:
	$(PYTHON) -m k3tasks report $(RUN_DIR) $(REPORT_ARGS)

venv:
	$(PYTHON) -m venv $(VENV)
	$(VENV)/bin/pip install -r experiments/requirements.txt

data:
	echo "$(SHA256)  $(ARCHIVE)" | sha256sum -c
	tar --skip-old-files -xJf $(ARCHIVE)

# work directories that were kept because a rerun differed or was stopped
clean-work:
	rm -rf results/work
