# Data

The stored results of the scripts are **not in git**. They are distributed as one compressed archive, to be attached to a release of this repository.

| | |
|---|---|
| File | `k3-derived-isogenies-data-v1.tar.xz` |
| Size | 3,447,652 bytes (3.3 MiB) compressed; 109,790,993 bytes (about 110 MB) unpacked |
| Files | 171 |
| SHA-256 | `7a16e7d759533558269147fef4b931c3cdaae66eae78fd41463cf2c9977b3714` |
| Where to get it | the `v1` release of the repository: <https://github.com/chris-dare-dev/k3-derived-isogenies/releases/download/v1/k3-derived-isogenies-data-v1.tar.xz> |

## Installing it

From the root of a clone of this repository:

```bash
echo "7a16e7d759533558269147fef4b931c3cdaae66eae78fd41463cf2c9977b3714  k3-derived-isogenies-data-v1.tar.xz" | sha256sum -c
tar -xJf k3-derived-isogenies-data-v1.tar.xz
```

The archive contains one directory, `experiments/`, and unpacks its files next to the scripts, which is where the scripts look for them: most resolve the paths relative to their own location, and all of them work when run from `experiments/`, as the notes do. It adds files only; it does not overwrite any script. The `.gitignore` keeps the unpacked files out of git. The notebooks do the same through `notebooks/k3data.py`, which checks the SHA-256 first.

## What is in it

Every stored result file behind the notes: the `.jsonl`, `.json` and `.log` outputs of the scripts in `experiments/`, whose names begin with `sixlines_` or `ellk3_`. These are the files the notes quote and the files that the ledger (`sixlines_ledger.py`, [note 27](docs/27-hypothesis-ledger.md)) and the leakage audit (`sixlines_leakage.py`, [note 28](docs/28-leakage-audit.md)) read. Each note's *Reproduce* section names the files its commands write; [docs/README.md](docs/README.md) gives the script of each note.

| Files | Written by | Notes | Content |
|---|---|---|---|
| `sixlines_p{5,…,43}.jsonl`, `sixlines_q25.jsonl`, `sixlines_q49.jsonl` | `sixlines.py`, `sixlines_fq.py` | 2, 3, 26, 29 | Every configuration of six rational lines over F_q, with its character sums and associate |
| `sixlines_cm.jsonl`, `sixlines_cm_p41.jsonl`, `sixlines_cm_p43.jsonl`, `sixlines_q49_cm.jsonl`, `sixlines_neighbours.jsonl` and their `.log` files | `sixlines_cm.py` | 2, 3, 5, 26, 29 | Hermitian lattice counts and neighbour steps per Frobenius polynomial |
| `sixlines_enumerate_29_37.log`, `sixlines_enumerate_41.log`, `sixlines_enumerate_43.log`, `sixlines_q49.log` | `sixlines.py`, `sixlines_fq.py` | 2, 3, 26, 29 | Logs of the enumerations |
| `sixlines_forms*`, `sixlines_types*`, `sixlines_nonsplit*`, `sixlines_rational*` | the scripts of the same names | 4, 5, 6 | Forms over the prime field, non-split forms, rational steps |
| `sixlines_vertical*`, `sixlines_levels*`, `sixlines_ramified*`, `sixlines_deep*`, `sixlines_rank6*`, `sixlines_hermitian*` | the scripts of the same names | 7–12 | Lattice graphs for ring-changing steps (lattice side only) |
| `sixlines_groupsize.*` | `sixlines_groupsize.py` | 13 | The group G and its size for admissible polynomials up to q = 100,003 |
| `sixlines_prym*.jsonl` and their `.log` files | `sixlines_prym.py`, `sixlines_prym_fq.py`, `sixlines_delta41.py` | 14, 16, 26 | The fourfold's Frobenius quartic per pair of lines, against the surface |
| `sixlines_dictionary*`, `sixlines_review*` | `sixlines_dictionary.py`, `sixlines_review.py`, `sixlines_review2.py` | 15, 18, 25 | Class-group dictionary; the re-derivation checks |
| `sixlines_characters*` | `sixlines_characters.py` | 19 | Candidate quadratic characters on surfaces |
| `sixlines_cmmethod*` | `sixlines_cmmethod.py` | 20 | CM method: class polynomials, labels, action, self-test |
| `sixlines_forms_fp*` | `sixlines_forms_fp.py` | 21 | Fourfold forms for sign − and reducible P |
| `ellk3_*` | `ellk3_two_neighbours.py` | 22 | Elliptic K3 pairs below Picard number 16, with controls |
| `sixlines_jacobian*`, `sixlines_delta_labels*`, `sixlines_delta41*` | the scripts of the same names | 23, 24, 26 | Jacobian arithmetic, δ, and its comparison with the CM cosets |
| `sixlines_ledger.log` | `sixlines_ledger.py` | 27 | The ledger's output (the reference for the reproduction check) |
| `sixlines_leakage*` | `sixlines_leakage.py` | 28 | Features, audit records and logs of the leakage audit |
| `sixlines_order3*` | `sixlines_order3.py`, `sixlines_order3_holdout43.py` | 29 | Cosets of ker χ₃, features, audit records and logs of the order-3 audit; with `_q43`, the run with q = 43 (§7.4); the held-out test of its survivor |

Not in it: the texts of the papers the notes cite.

## Notes on the files

- `sixlines_prym.jsonl` and `sixlines_prym.log` were written by the current `sixlines_prym.py` (the command of [note 14](docs/14-fourfold-of-a-six-line-surface.md), 8 workers), after the classification fix described in [note 16](docs/16-fourfold-over-f25-and-f49.md): the field `degrees` lists the degrees of the factors of P with multiplicity. The files read by the ledger, and the counts quoted in the notes, depend only on whether P is irreducible, which the fix does not change.
- In `sixlines_leakage_features.log`, the line reporting the output file is the one the script prints (`wrote sixlines_leakage_features.json`).
- In `sixlines_enumerate_41.log`, a timing log written by `/usr/bin/time -v`, the absolute path of the Python interpreter is replaced by `python`. The same is done in `sixlines_enumerate_43.log`, `sixlines_order3_audit.log`, `sixlines_order3_q43.log` and `sixlines_order3_cosets_43_2.log`.
- In `sixlines_order3_cosets_43_2_mix.log`, `_43_2_coord.log` and `_43_2_mix4500.log`, the absolute path of the cache of period points is replaced by `sixlines_order3_points/`, where the script in this repository keeps it. That cache itself is not in the archive.
- The archive stores no user or group names (owner 0, group 0), and the files keep their modification times.

## How it was built

```bash
tar --sort=name --owner=0 --group=0 --numeric-owner --format=gnu -cf - experiments | xz -9 > k3-derived-isogenies-data-v1.tar.xz
sha256sum k3-derived-isogenies-data-v1.tar.xz > k3-derived-isogenies-data-v1.tar.xz.sha256
```

run in a staging directory holding `experiments/` with the 171 files of the manifest below.

## Verification

Run on 2026-11-15 (Linux, Python 3.12, numpy 2.5.4, cypari2 2.2.4) on a fresh clone of this repository with the archive unpacked:

| Run | Result |
|---|---|
| SHA-256 of the archive against the value above | OK |
| `sixlines_ledger.py`, compared with the stored `sixlines_ledger.log` | identical: 170 comparisons, 162 agree, 8 differ (the slips described in note 27 §1) |
| `sixlines_cmmethod.py selftest --digits 60` | `selftest: passed`; identical to the stored `sixlines_cmmethod_selftest.log` |
| `sixlines_prym.py --check 17 2 3 4 11` | 30 ordered pairs: 18 equal, 2 twist, 10 degenerate, no mismatch |
| `sixlines_prym.py --workers 4 --run 11 13` | 41 records, identical to the stored records for p = 11 and 13 |
| `notebooks/01-explore-the-results.ipynb` | runs in about 3 s; ledger identical; 512 groups, median \|G\|/q² 0.0092 |
| `notebooks/02-small-live-instances.ipynb` | runs in about 12 s; the enumeration over F₁₃ is identical to the stored file; CM action at q = 19: 18 of 18 steps equal Frobenius |

The archive was rebuilt on 2026-11-16 with the 26 files of note 29, 145 to 171 files. On a fresh copy of this repository with the new archive unpacked (same software):

| Run | Result |
|---|---|
| SHA-256 of the archive against the value above | OK |
| `sixlines_ledger.py`, compared with the stored `sixlines_ledger.log` | identical |
| `sixlines_order3.py features`, then `audit --workers 10` | `sixlines_order3_features.json` and `sixlines_order3.jsonl` identical to the stored files |
| `sixlines_order3.py features --with-43 --out sixlines_order3_features_q43.json`, then `audit --workers 10 --features sixlines_order3_features_q43.json --out sixlines_order3_q43.jsonl` | both identical to the stored files; the logs identical up to the times |
| `sixlines_order3_holdout43.py` | identical to the stored `sixlines_order3_holdout43.log` |
| `notebooks/01-explore-the-results.ipynb` and `02-small-live-instances.ipynb`, with the archive in `dist/` and no data in `experiments/` | both unpack the archive through `k3data.py`'s SHA-256 check and run (2.5 s and 12 s); ledger identical; the enumeration over F₁₃ identical to the stored file |

The other rows of the first table were not rerun with the new archive.

## Rerunning against the stored files

The batches of `batches/` (see "Running a batch" in the [README](README.md#running-a-batch)) rerun the commands of the notes and compare what they write with these files. They read the files from `experiments/`, or from the directory named by `K3_DATA_DIR`, copy them into a work directory under `results/work/` and never write to them; a rerun's own output stays under `results/`, which git ignores. `make data` runs the two commands above and skips files that are already there.

The records of the batches that were run are in `records/`, and [CLAIMS.md](CLAIMS.md) says what they show for each result. They were computed on Python 3.12.14 with numpy 2.5.4, cypari2 2.2.4 (PARI 2.17.2) and cysignals 1.12.6, the versions of `experiments/requirements.txt`; each record names its interpreter, numpy and PARI. The stored results are documented as runs on Python 3.12 with the same numpy and cypari2 (the verification above); the stored files do not record their interpreter. Two scripts, `sixlines_fq.py` and `sixlines_nonsplit.py`, rely on the fork start method, which is the default on Linux up to Python 3.13; on a later Python they run only through the harness, which sets it ([CLAIMS.md](CLAIMS.md#known-differences-between-reruns-and-stored-files)).

## Licence

The data are released under [CC BY 4.0](LICENSE-TEXT.md), like the text.

## Manifest

| File | Bytes | Date |
|---|---|---|
| `ellk3_e8_p5.jsonl` | 739981 | 2026-11-14 |
| `ellk3_e8_p5.log` | 1391 | 2026-11-14 |
| `ellk3_e8_p5_control.jsonl` | 207213 | 2026-11-14 |
| `ellk3_e8_p5_control.log` | 654 | 2026-11-14 |
| `ellk3_e8a2_graph.log` | 1017 | 2026-11-14 |
| `ellk3_e8a2_p5.jsonl` | 9848106 | 2026-11-14 |
| `ellk3_e8a2_p5.log` | 1464 | 2026-11-14 |
| `ellk3_e8a2_p5_control.jsonl` | 512972 | 2026-11-14 |
| `ellk3_e8a2_p5_control.log` | 698 | 2026-11-14 |
| `ellk3_e8a2_p7.jsonl` | 2519750 | 2026-11-14 |
| `ellk3_e8a2_p7.log` | 2018 | 2026-11-14 |
| `ellk3_generic_p5.jsonl` | 6929 | 2026-11-14 |
| `ellk3_generic_p5.log` | 637 | 2026-11-14 |
| `ellk3_generic_p5_control.jsonl` | 3845 | 2026-11-14 |
| `ellk3_generic_p5_control.log` | 362 | 2026-11-14 |
| `ellk3_generic_p5_seed2.jsonl` | 21522 | 2026-11-14 |
| `ellk3_generic_p5_seed2.log` | 678 | 2026-11-14 |
| `sixlines_characters.jsonl` | 392391 | 2026-11-14 |
| `sixlines_characters.log` | 9520 | 2026-11-14 |
| `sixlines_characters_dyadic.log` | 340 | 2026-11-14 |
| `sixlines_characters_genus.jsonl` | 163094 | 2026-11-14 |
| `sixlines_characters_genus.log` | 1209 | 2026-11-14 |
| `sixlines_cm.jsonl` | 30532 | 2026-11-10 |
| `sixlines_cm.log` | 7271 | 2026-11-10 |
| `sixlines_cm_p41.jsonl` | 20397 | 2026-11-15 |
| `sixlines_cm_p41.log` | 3103 | 2026-11-15 |
| `sixlines_cm_p43.jsonl` | 23748 | 2026-11-15 |
| `sixlines_cm_p43.log` | 3472 | 2026-11-15 |
| `sixlines_cmmethod.jsonl` | 37607 | 2026-11-14 |
| `sixlines_cmmethod.log` | 51162 | 2026-11-14 |
| `sixlines_cmmethod_rerun.log` | 33309 | 2026-11-14 |
| `sixlines_cmmethod_rerun2.log` | 3178 | 2026-11-14 |
| `sixlines_cmmethod_selftest.log` | 1556 | 2026-11-14 |
| `sixlines_deep.jsonl` | 87949 | 2026-11-11 |
| `sixlines_deep.log` | 47934 | 2026-11-11 |
| `sixlines_delta41.jsonl` | 388811 | 2026-11-15 |
| `sixlines_delta41.log` | 7248 | 2026-11-15 |
| `sixlines_delta41_compare.json` | 9086 | 2026-11-15 |
| `sixlines_delta41_compare.log` | 7329 | 2026-11-15 |
| `sixlines_delta41_cosets.jsonl` | 10889 | 2026-11-15 |
| `sixlines_delta41_cosets_22.log` | 5867 | 2026-11-15 |
| `sixlines_delta41_cosets_6.log` | 7938 | 2026-11-15 |
| `sixlines_delta41_more.log` | 959 | 2026-11-15 |
| `sixlines_delta41_primary.log` | 313 | 2026-11-15 |
| `sixlines_delta41_signs.jsonl` | 13240 | 2026-11-15 |
| `sixlines_delta41_signs.log` | 13741 | 2026-11-15 |
| `sixlines_delta_labels.jsonl` | 6024 | 2026-11-14 |
| `sixlines_delta_labels.log` | 10304 | 2026-11-14 |
| `sixlines_dictionary.jsonl` | 103519 | 2026-11-12 |
| `sixlines_dictionary.log` | 487 | 2026-11-12 |
| `sixlines_dictionary_association.log` | 1303 | 2026-11-12 |
| `sixlines_enumerate_29_37.log` | 369 | 2026-11-10 |
| `sixlines_enumerate_41.log` | 942 | 2026-11-15 |
| `sixlines_enumerate_43.log` | 942 | 2026-11-15 |
| `sixlines_forms_fp.jsonl` | 655504 | 2026-11-14 |
| `sixlines_forms_fp.log` | 44257 | 2026-11-14 |
| `sixlines_forms_p5.jsonl` | 10605 | 2026-11-10 |
| `sixlines_forms_p5.log` | 314 | 2026-11-10 |
| `sixlines_forms_p5_compare.jsonl` | 176 | 2026-11-10 |
| `sixlines_forms_p5_compare.log` | 1007 | 2026-11-10 |
| `sixlines_forms_p7.jsonl` | 41801 | 2026-11-10 |
| `sixlines_forms_p7.log` | 324 | 2026-11-10 |
| `sixlines_forms_p7_compare.jsonl` | 514 | 2026-11-10 |
| `sixlines_forms_p7_compare.log` | 3166 | 2026-11-10 |
| `sixlines_groupsize.jsonl` | 221989 | 2026-11-12 |
| `sixlines_groupsize.log` | 16357 | 2026-11-12 |
| `sixlines_hermitian.jsonl` | 49922 | 2026-11-11 |
| `sixlines_hermitian.log` | 3523 | 2026-11-11 |
| `sixlines_hermitian_q7.jsonl` | 506769 | 2026-11-11 |
| `sixlines_hermitian_q7.log` | 58511 | 2026-11-11 |
| `sixlines_jacobian_character.jsonl` | 58501 | 2026-11-14 |
| `sixlines_jacobian_character.log` | 11735 | 2026-11-14 |
| `sixlines_jacobian_character_extra.jsonl` | 31792 | 2026-11-14 |
| `sixlines_jacobian_character_extra.log` | 1123 | 2026-11-14 |
| `sixlines_jacobian_lattice.log` | 3567 | 2026-11-14 |
| `sixlines_jacobian_levels.jsonl` | 118354 | 2026-11-14 |
| `sixlines_jacobian_levels.log` | 3180 | 2026-11-14 |
| `sixlines_jacobian_validate.jsonl` | 16136 | 2026-11-14 |
| `sixlines_jacobian_validate.log` | 1013 | 2026-11-14 |
| `sixlines_leakage.jsonl` | 5895633 | 2026-11-15 |
| `sixlines_leakage.log` | 7266 | 2026-11-15 |
| `sixlines_leakage_extra.jsonl` | 90843 | 2026-11-15 |
| `sixlines_leakage_extra.log` | 1929 | 2026-11-15 |
| `sixlines_leakage_features.json` | 21480354 | 2026-11-15 |
| `sixlines_leakage_features.log` | 3257 | 2026-11-15 |
| `sixlines_ledger.log` | 42040 | 2026-11-15 |
| `sixlines_levels.jsonl` | 1364 | 2026-11-11 |
| `sixlines_levels.log` | 1244 | 2026-11-11 |
| `sixlines_neighbours.jsonl` | 38465 | 2026-11-10 |
| `sixlines_neighbours.log` | 13466 | 2026-11-10 |
| `sixlines_nonsplit_lat_p17.jsonl` | 6216 | 2026-11-10 |
| `sixlines_nonsplit_lat_p17.log` | 4561 | 2026-11-10 |
| `sixlines_nonsplit_lat_p19.jsonl` | 7703 | 2026-11-10 |
| `sixlines_nonsplit_lat_p19.log` | 5593 | 2026-11-10 |
| `sixlines_nonsplit_lat_p23.jsonl` | 12982 | 2026-11-10 |
| `sixlines_nonsplit_lat_p23.log` | 9337 | 2026-11-10 |
| `sixlines_nonsplit_lat_small.jsonl` | 5682 | 2026-11-10 |
| `sixlines_nonsplit_lat_small.log` | 4615 | 2026-11-10 |
| `sixlines_nonsplit_p11.jsonl` | 234042 | 2026-11-10 |
| `sixlines_nonsplit_p11.log` | 290 | 2026-11-10 |
| `sixlines_nonsplit_p13.jsonl` | 466611 | 2026-11-10 |
| `sixlines_nonsplit_p13.log` | 291 | 2026-11-10 |
| `sixlines_nonsplit_p17.jsonl` | 1407141 | 2026-11-10 |
| `sixlines_nonsplit_p17.log` | 294 | 2026-11-10 |
| `sixlines_nonsplit_p19.jsonl` | 2222839 | 2026-11-10 |
| `sixlines_nonsplit_p19.log` | 299 | 2026-11-10 |
| `sixlines_nonsplit_p23.jsonl` | 4862099 | 2026-11-10 |
| `sixlines_nonsplit_p23.log` | 302 | 2026-11-10 |
| `sixlines_order3.jsonl` | 5117077 | 2026-11-15 |
| `sixlines_order3_audit.log` | 13092 | 2026-11-15 |
| `sixlines_order3_cosets.jsonl` | 26855 | 2026-11-16 |
| `sixlines_order3_cosets_37.log` | 8774 | 2026-11-15 |
| `sixlines_order3_cosets_41.log` | 3777 | 2026-11-15 |
| `sixlines_order3_cosets_43.log` | 6303 | 2026-11-15 |
| `sixlines_order3_cosets_43_14.log` | 1819 | 2026-11-15 |
| `sixlines_order3_cosets_43_2.log` | 6859 | 2026-11-16 |
| `sixlines_order3_cosets_43_2_coord.log` | 1472 | 2026-11-15 |
| `sixlines_order3_cosets_43_2_mix.log` | 1754 | 2026-11-15 |
| `sixlines_order3_cosets_43_2_mix4500.log` | 1271 | 2026-11-16 |
| `sixlines_order3_features.json` | 11306949 | 2026-11-15 |
| `sixlines_order3_features.log` | 1771 | 2026-11-15 |
| `sixlines_order3_features_q43.json` | 12173369 | 2026-11-15 |
| `sixlines_order3_features_q43.log` | 2196 | 2026-11-15 |
| `sixlines_order3_fields.log` | 873 | 2026-11-15 |
| `sixlines_order3_holdout43.log` | 385 | 2026-11-16 |
| `sixlines_order3_labelcheck.log` | 244 | 2026-11-15 |
| `sixlines_order3_q43.jsonl` | 5718240 | 2026-11-15 |
| `sixlines_order3_q43.log` | 17180 | 2026-11-15 |
| `sixlines_order3_selfcheck.log` | 496 | 2026-11-15 |
| `sixlines_order3_selfcheck_q43.log` | 696 | 2026-11-15 |
| `sixlines_p11.jsonl` | 1348 | 2026-11-10 |
| `sixlines_p13.jsonl` | 2350 | 2026-11-10 |
| `sixlines_p17.jsonl` | 6830 | 2026-11-10 |
| `sixlines_p19.jsonl` | 10837 | 2026-11-10 |
| `sixlines_p23.jsonl` | 24032 | 2026-11-10 |
| `sixlines_p29.jsonl` | 64590 | 2026-11-10 |
| `sixlines_p31.jsonl` | 85885 | 2026-11-10 |
| `sixlines_p37.jsonl` | 183155 | 2026-11-10 |
| `sixlines_p41.jsonl` | 284262 | 2026-11-15 |
| `sixlines_p43.jsonl` | 347851 | 2026-11-15 |
| `sixlines_p5.jsonl` | 86 | 2026-11-10 |
| `sixlines_p7.jsonl` | 260 | 2026-11-10 |
| `sixlines_prym.jsonl` | 4177245 | 2026-11-15 |
| `sixlines_prym.log` | 10334 | 2026-11-15 |
| `sixlines_prym_p41.jsonl` | 285088 | 2026-11-15 |
| `sixlines_prym_q25.jsonl` | 678463 | 2026-11-12 |
| `sixlines_prym_q25.log` | 906 | 2026-11-12 |
| `sixlines_prym_q49.jsonl` | 741322 | 2026-11-12 |
| `sixlines_prym_q49.log` | 545 | 2026-11-12 |
| `sixlines_q25.jsonl` | 47319 | 2026-11-10 |
| `sixlines_q49.jsonl` | 836344 | 2026-11-10 |
| `sixlines_q49.log` | 205 | 2026-11-10 |
| `sixlines_q49_cm.jsonl` | 41744 | 2026-11-10 |
| `sixlines_q49_cm.log` | 9882 | 2026-11-10 |
| `sixlines_ramified.jsonl` | 1175637 | 2026-11-11 |
| `sixlines_ramified.log` | 35462 | 2026-11-11 |
| `sixlines_ramified_q7.jsonl` | 11526252 | 2026-11-11 |
| `sixlines_ramified_q7.log` | 352557 | 2026-11-11 |
| `sixlines_rank6.jsonl` | 84542 | 2026-11-11 |
| `sixlines_rank6.log` | 8241 | 2026-11-11 |
| `sixlines_rational.jsonl` | 31560 | 2026-11-11 |
| `sixlines_rational.log` | 6285 | 2026-11-11 |
| `sixlines_review.log` | 13113 | 2026-11-14 |
| `sixlines_review2.jsonl` | 5995 | 2026-11-15 |
| `sixlines_review2.log` | 13927 | 2026-11-15 |
| `sixlines_types_q25.jsonl` | 856 | 2026-11-10 |
| `sixlines_types_q25.log` | 987 | 2026-11-10 |
| `sixlines_types_q49.jsonl` | 4769 | 2026-11-10 |
| `sixlines_types_q49.log` | 5038 | 2026-11-10 |
| `sixlines_vertical.jsonl` | 180177 | 2026-11-11 |
| `sixlines_vertical.log` | 30042 | 2026-11-11 |
