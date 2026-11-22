# Python 3.12 with numpy and cypari2 at the versions of research/experiments/requirements.txt, make, git, xz
# and the lint tools, for machines without them: the interpreter of the stored results and of the records.
#   docker build -t k3-derived-isogenies .
#   docker run --rm --network none --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work -w /work k3-derived-isogenies make test lint check
# After `make data` (data/README.md), `make smoke` runs there too; without systemd the memory limits
# are enforced by polling (ARGS="--memory-guard poll" is the default the scheduler falls back to).
FROM python:3.12-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends make git xz-utils \
 && rm -rf /var/lib/apt/lists/*
COPY research/experiments/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --require-hashes -r /tmp/requirements.txt \
 && pip install --no-cache-dir ruff==0.17.5 mypy==2.4.1
