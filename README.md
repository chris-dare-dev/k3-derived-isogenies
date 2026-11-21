# Derived isogenies of K3 surfaces

<!--
  Badges: the GitHub-dynamic ones (CI, last commit) render once the repo is
  public. The static ones always render.
-->
[![CI](https://github.com/chris-dare-dev/k3-derived-isogenies/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/chris-dare-dev/k3-derived-isogenies/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](pyproject.toml)
[![Code: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![Text and data: CC BY 4.0](https://img.shields.io/badge/text%20%26%20data-CC%20BY%204.0-green.svg)](LICENSE-TEXT.md)
[![Status](https://img.shields.io/badge/status-not%20peer%20reviewed-orange.svg)](#status)
[![Last commit](https://img.shields.io/github/last-commit/chris-dare-dev/k3-derived-isogenies.svg)](https://github.com/chris-dare-dev/k3-derived-isogenies/commits)

**Does a class-group action on ordinary K3 surfaces give a new hard problem for cryptography? No.**
Notes, code and stored data behind a negative result for isogeny-style cryptography.

Ordinary K3 surfaces over a finite field with the same Frobenius form torsors under a CM class group, which acts by twisted derived equivalences ("derived isogenies"), much as in CSIDH.

## Result

- **It reduces to abelian fourfolds.** For double planes branched along six lines (Picard number 16) the action is a CSIDH-type action on fourfolds with Z[i]-multiplication. Checked on 278 surfaces and 4,983 pairs of lines, with no exception.
- **It is worse than CSIDH.** The group has about 0.01·q² elements, so a 256-bit group needs about 525-bit instances (512 for CSIDH-512), and a step costs about ℓ⁴ instead of ℓ.
- **One bit leaks.** A quadratic character, for 80 of the 477 examined polynomials that have one. An audit of about 11,400 candidate invariants and a second one for order-3 characters found nothing else.
- **Below Picard number 16** there is no known reduction and no computable step.

## Read

| | |
|---|---|
| [SUMMARY.md](SUMMARY.md) | Self-contained summary paper, about 12 pages |
| [docs/](docs/README.md) | 29 notes: arguments, computations, checks |
| [CLAIMS.md](CLAIMS.md) | Each claim, the script behind it, and what its stored record says |
| [DATA.md](DATA.md) | The stored results and how to install them |
| [GUIDE.md](GUIDE.md) | Results table, limits, layout, running batches |

## Reproduce

```bash
make test                    # scheduler, registry and claim tables: standard library only, nothing to install
make venv && make no-data    # numpy and cypari2, then the six computations that need no stored results (a minute)
make data && make smoke      # the stored results (DATA.md), then the ledger of note 27 and a dozen quick checks
```

Or use the [notebooks](notebooks/README.md) with no command line, or `docker build -t k3-derived-isogenies .`.

## Status

All measurements are on small instances (q ≤ 49); at cryptographic size the verdict is argument, not measurement. The idea and direction are mine, and Claude was used extensively for the derivations, code, computations and drafting. Nothing here has been independently reviewed.

Cite via [CITATION.cff](CITATION.cff).
