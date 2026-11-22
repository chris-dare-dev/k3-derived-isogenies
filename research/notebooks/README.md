# Notebooks

Two notebooks for readers who prefer not to use a command line.

| Notebook | What it does | Time |
|---|---|---|
| [01-explore-the-results.ipynb](01-explore-the-results.ipynb) | Recomputes the headline numbers from the stored results: the ledger, the size of the group, δ against the cosets of ker χ₋₁, the leakage audit | a few seconds |
| [02-small-live-instances.ipynb](02-small-live-instances.ipynb) | Runs small instances live: every six-line surface over F₁₃, the abelian fourfold of one surface over F₁₇, the CM method and the group action at q = 19, the self-test of the theta formulas | about 20 seconds |

**In the browser.** The Binder badge at the top of each notebook starts a free server with everything installed (`.binder/` holds the configuration: Python 3.12, numpy, cypari2 from its binary wheel, matplotlib; `postBuild` fetches the data).

**Locally.** Install the requirements of `.binder/requirements.txt` and Jupyter, put the data archive in `dist/` or the repository root (or set `K3DI_DATA_ARCHIVE`), and open the notebooks from this folder. `k3data.ensure_data()` checks the archive against the SHA-256 in [data/README.md](../../data/README.md) and unpacks it into `research/experiments/`.

The notebooks call the scripts in `research/experiments/` exactly as the notes do (`python -I script …`), so what they print is what the notes quote. They write nothing into `research/experiments/`.

**Longer reruns.** The notebooks keep to computations of seconds. To rerun the other commands of the notes, compared with the stored files and under time and memory limits, use the batches from a terminal: `make run BATCH=data/batches/smoke.toml` repeats the notebooks' checks of the ledger, the fourfold over F₁₇ and the theta formulas in half a minute, and the other files of `data/batches/` cover the claim families of the README, up to `heavy.toml` (see "Running a batch" in [GUIDE.md](../../GUIDE.md#running-a-batch)).
