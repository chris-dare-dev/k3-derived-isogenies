# 27. A ledger of the hardness hypotheses

*A ledger of every claim that notes [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md), [17](17-below-picard-16.md), [18](18-checks-of-notes-13-to-16.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [21](21-fourfold-forms-over-fp.md), [22](22-two-neighbours-of-elliptic-k3.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md) and [26](26-delta-is-the-character-at-q-41.md) make about the hardness question of note 13 §2, and of the evidence for each, recounted from the stored `.jsonl`, `.json` and `.log` files by `experiments/sixlines_ledger.py`. Nothing was recomputed beyond the stored files: no surface, curve, lattice or CM point. It prepares the leakage audit of surface invariants against the CM labels, which is [note 28](28-leakage-audit.md), and the rows touched by the test of characters of order 3 in [note 29](29-characters-of-order-3.md) are brought up to date. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `sixlines_ledger.py` from stored files; **[H]** heuristic; **[R]** re-derived in a separate re-derivation check, with its own code (not part of this repository). Labels of other notes are kept as they are there. Checked once, by a separate re-derivation check that recounted all 170 comparisons and about 50 further numbers with its own code from the stored files, before reading the ledger's script; every comparison reproduces. Its corrections are made in the text; what it did not cover is in §6. The classification itself is not re-checked.*

---

## 1. Result

**Mismatches between the notes and their files [E].** The script makes 170 comparisons between a number quoted in a note and the same number recounted from the files. 162 agree. The 8 that differ are 7 distinct numbers (the row for q = 17 is compared twice) and come down to three issues; none of them changes a conclusion. 8 of the 170 comparisons read the recounted number from a summary line of a log rather than from records [R].

1. **Note 13 §3.3, the statistics of |G|/q².** On the 512 groups that note 13 §8 defines (maximal orders, q in the table, without q = 125, 343, 625, 2,401, 3,125 and without the 24 duplicate "real" records), `sixlines_groupsize.jsonl` gives median 0.0092 and quartiles 0.0067 and 0.0128. The values that had been given for the 512, 0.0095, 0.0068 and 0.0129, are those of the 486 groups with q ≥ 41, the population of note 13's first fit: 0.00952, 0.00680 and 0.01289, under every one of 13 quantile rules, with the range 0.0003 to 0.024 [R]. (Over all 636 maximal records of the file the values are 0.0095, 0.0067 and 0.0129, with maximum 0.031.) On the 512, the range (0.0003 to 0.024), the two fitted exponents (2.01 and 1.98), the constant 0.011 and the group counts all match. The statement |G| ≈ 0.01·q² is not affected. Note 13 now states the recounted values on the 512, with those of the 486 beside them.
2. **Note 13 §3.3, the table row for q = 17.** The median of |G| in this row had been rounded to 2; the file gives 1.5 (the groups have orders 1 and 2), and the column |G|/q² = 0.0052 in the same row is computed from 1.5. The other 30 rows match the file in every column. Note 13 now states the median as 1.5.
3. **Configurations counted as surfaces.** Note 14 §1 item 2, note 13 §1, note 18 and note 21 had called the 534 surfaces. They are configurations: six lines up to projective equivalence and relabelling, as the table of note 14 §4 calls them. A surface is a configuration together with its associate (note 20 §2), and the 534 configurations make up 278 surfaces. In the same way, the 37 objects over F₂₅ and the 1,108 over F₄₉ of note 16 §1 item 2, which had been called surfaces, are 37 and 1,108 configurations, which make up 21 and 564 surfaces. The pair counts (4,983, 858, 3,324) are right. Note 14 §1 item 2 and note 16 §1 item 2 now state the counts as configurations, with the numbers of surfaces they make up, and notes 6, 13, 14 (§1, §5.2, §7), 16 (§1, §3, §4), 18 and 21 now call configurations what they count. The other counts of configurations in these notes make up these numbers of surfaces [R]: 1,145 configurations with irreducible P over F₂₅ and F₄₉, 585 surfaces; the 1,090 of them not defined over the prime field, 556; 1,329 with sign + and reducible P that have a computed pair, 1,124; all 1,564 with sign + and reducible P, 1,359; 1,892 with sign −, 1,447; over F₂₅, 131, 160 and 168 configurations (sign + and reducible P with a computed pair, all with sign + and reducible P, sign −), 119, 148 and 136; 1,080 ordinary configurations over F₄₉, 540.

**Quoted counts with no stored file.** Two counts cannot be recounted from any stored file. Each needs a rerun, which is named in the Reproduce section:

- note 13 §4.2: the 23,945 polynomials of the Proposition Y check (only the 265 of `--fourfold 41 61 101` are in the log);
- note 20 §5.1: the 234 odd lattices of the maximal order (a loop that is not a command).

Three further counts are not recounted by the ledger, but their files exist and an independent recount from them agrees [R]:

- note 23 §4: the counts of lattice classes behind the level test are in `sixlines_jacobian_lattice.log`, a log the ledger does not parse; they recount as 1/6, 5/14 and 1/2;
- note 24 §4: the statement that with three primes per field in place of eight the groups and orbits are the same. Note 26 §7.1 reran the three torsors with three primes into `sixlines_delta41_signs.jsonl`, and the partition equals that of the eight-prime runs in 3 of 3;
- note 21 §1 item 9 and §8: the 14,394 non-degenerate pairs are a stored line of `sixlines_forms_fp.log`; only their per-pair data are not stored.

**The ledger.** §2 lists 47 hypotheses. By class:

| class | hypotheses |
|---|---|
| proved on paper (sketches included, with their label) | H1, H6, H14, H40 |
| supported by data | H2–H5, H7–H13, H15–H25, H27, H28, H41 |
| heuristic | H26, H31–H33, H35–H37, H46, H47 |
| open | H29, H30, H34, H38, H42–H45 |
| contradicted (superseded by a later note) | H39, and the superseded forms of H5 (note 15 §3.1 "Why" and §3.3, corrected by note 18) and H15 (note 13 §6) |

**What carries the verdict.** The surfaces in the data have q ≤ 49 (and p = 41 at the last). The groups have q ≤ 100,003, and the largest torsor on which the action was evaluated has 16 elements. The verdict "a CSIDH-type problem on abelian fourfolds, one bit of leak where χ₋₁ exists" rests at cryptographic sizes on five main things. Four of them are arguments and one is an extrapolation:

- Proposition Y, an identity valid for all q (H1, [D], re-derived in note 18);
- the equivariance of surface ↦ fourfold (H29, a sketch with five gaps, never tested at a prime ℓ ≠ p on the Prym variety itself);
- the CM points being the canonical lifts (H30, a sketch);
- the generic attacks and the cost of zeta functions and steps (H35 to H37, [M], with [P] for Costa–Harvey–Kedlaya in H36);
- the law |G| ≈ 0.01·q², extrapolated from q ≤ 10⁵ to q ≈ 2¹³¹ (H10, [E] in range, [H] beyond).

It also uses:

- Deligne's description of ordinary abelian varieties with polarisations (H31, [M], not read), which H32 names as a condition;
- the torsor structure: Propositions A to C of note 6 and Lemma 9 of note 13, which have no row of their own, and H11 beyond its range, on which H14 and H40 rest;
- for "one bit of leak where χ₋₁ exists", H38, which is open beyond 5 of the 80 polynomials;
- for its being only one bit, H42 and H43, which are open.

No measurement in this repository is at a size where the problem could be hard. The data test the dictionary, the group and the leak on small instances, and nothing else (§4).

---

## 2. The ledger

**Columns.** "Source" is the note, by number, and the section where the claim is made. "Evidence" is the count recomputed from the files, agreeing of tested, and the population. "q" is the range covered. "Falsifier" is the result that would refute the hypothesis, and where it is in reach: **S** in the stored data, **R** by a rerun or extension of an existing script, **N** in reach of neither. "Audit" marks the hypotheses the leakage audit (§5, [note 28](28-leakage-audit.md)) bears on: **yes** if a correlation between a surface invariant and the CM labels would confirm or refute it; **uses** if the audit relies on it being true. Unless marked otherwise, every class-group count is under GRH (PARI without certification). The hypotheses fall into five groups, G1 to G5; `sixlines_ledger.py` prints them in that order, under the headings GROUP G1 to GROUP G5.

| ID | Statement | Source | Class | Evidence (recomputed) | q | Falsifier (reach) | Audit |
|---|---|---|---|---|---|---|---|
| **G1** | **Reduction to abelian fourfolds** | | | | | | |
| H1 | P(±q) = 64q⁵(μ², ν²) and P = Λ²χ for χ = X⁴ − 4(ν + μi)X³ + e₁X² − 4q(ν − μi)X + q² | 13 Prop. Y; 18 | proved on paper [D, re-derived by 18]; Honda–Tate [M] for the fourfold | 1,606 of 1,606 stored polynomials, of which 1,562 were built from (e₁, μ, ν) by the equations that Proposition Y restates and cannot fail; the 44 of enumerated surfaces test it; χχ̄ irreducible 32 of 44 | 13 to 100,003 | an admissible P that is not Λ²χ (S; excluded by the identity) | no |
| H2 | The curve of note 14 has Λ²χ = P or P(−t) for every surface of sign + with irreducible P | 14, 16, 26 | supported by data; an isogeny invariant only | 10,672 of 10,672 pairs on 1,735 configurations = 891 surfaces; reducible P: 12,873 + 2,829 pairs, 0 exceptions | p 13 to 41; 25, 49 | a pair giving another polynomial (R, any p) | no |
| H3 | χ of the curve is the quartic of H1 (c₃ = p·c̄₁, 4 ∣ c₁, {\|Re c₁\|, \|Im c₁\|} = {4μ, 4ν}) | 14 §4 | supported | 17,856 of 17,856 pairs (μ, ν from the stored e) | 11 to 37 | one pair failing (S: none) | no |
| H4 | Sign + and irreducible P: Frobenius commutes with J; no collinear triple, no anomalous pair | 14 §5.2, 18, 21 §8 | supported; [D] a collinear triple forces a root ±p | 0 inverting of 5,792 pairs seen; 0 of 534 configurations with a collinear triple (log) | 13 to 37 | an inverting pair on such a surface (R) | no |
| H5 | f : Cl(K) → Cl(E) is a homomorphism, f(Cl K⁺) ⊂ H, the odd parts of Cl⁻(K) and G are isomorphic. *Corrected:* f is induced by the square of the map of tori; the cokernel has an element of order 4 in 46 cases (18); the maximal-order lattices are not surfaces (20 §5.1) | 15 §3, 18 findings 1–2, 20 | supported, as corrected; [D] sketch of 18 for every P | 278 of 278 polynomials, 77,946 primes; kernel 6/205/64/3; cokernel 10/148/105/15 | 17 to 101 | a prime whose image is not the predicted class (S: none) | no |
| H6 | A configuration and its associate have the same curves | 15 §4 | proved [D] and measured | 52,207 of 52,207 pairs; 7,643 skipped | 5 to 37 | a mismatch (S: none) | no |
| H7 | The CM method's class polynomial mod 𝔭 has exactly the enumerated surfaces as roots | 20 §6, 24, 26, 29 | supported | 126 of 126 surfaces in 19 fields; under note 20's last-record rule 74 of 74; with the first records 74 of 87. Note 29 adds q = 37, e₁ = 6 and q = 41, e₁ = 2, 12 of 12 each: 150 surfaces in 21 fields | 17 to 41; \|G\| ≤ 16 | an extra or missing root (R at any small \|G\|) | uses |
| H8 | On even lattices of R the ideals act through G, freely, transitively; the step at 𝔩 is Frobenius conjugation in characteristic ℓ | 20 §7 | supported; Shimura–Taniyama [M]; orientation fitted on one field | 378 ideals, 0 conflicts; 390/390 compose; congruence 346/346; inverse fails 136/136; another pair fails 1,156/1,156; 594/594 squares | 17 to 37 | a conflict or failed congruence (R) | uses |
| H9 | The Prym variety separates the levels at ℓ as the lattices predict | 23 §4; weight: 25 finding 10 | supported, weakly | 422 pairs on 29 surfaces, 0 disagree within a surface; counts 1/6, 5/14, 1/2; 5 of 7 two-block surfaces self-associated | 29, 37 | a surface whose pairs disagree (S: none); the other three cases (R, needs F_{p⁸}, F_{p¹⁵}, F_{p²⁴}) | no |
| H29 | Surface ↦ fourfold is equivariant for the two actions | 15 §2, 18 finding 3 | **open**: a sketch with five gaps; indirect support from H23, H24 (an invariant of the Prym variety follows the CM action) | no direct test at ℓ ≠ p on a Prym variety | — | an isogeny A → A/A[𝔞] landing off the predicted surface (N: no method for isogenies of these Pryms, 15 §5.1) | uses |
| H30 | The CM points are the canonical lifts, so the action of note 20 is the action on the K3 surfaces | 20 §7.7, 18 finding 3 | **open** (sketch [D]) | none beyond H7 | — | a CM point reducing to a surface with a different lift (N) | uses |
| H31 | Self-dual hermitian Z[i]-lattices with π are the principally polarised ordinary fourfolds with Z[i] (Deligne, with polarisations); Honda–Tate gives the isogeny class | 13 §4.3, 15 §2, 18 §5 | heuristic [M, not read] | — | — | (N) | uses |
| H32 | Hence the K3 problem in this family is vectorisation for a class-group action on fourfolds with Z[i]: no harder, not a new assumption | 13 §1 item 4, 14 §6 | heuristic: conditional on H29, H30, Deligne with polarisations [M] | — | — | a failure of H29 or H30 (N) | no |
| **G2** | **Size and structure of the group** | | | | | | |
| H10 | \|G\| ≈ 0.01·q² for maximal orders; fitted exponent ≈ 2 | 13 §3.3 | supported in range; [H] the class-number heuristic beyond | 512 groups: median 0.0092 (recounted, see §1), range 0.0003 to 0.024; exponent 2.01 (q ≥ 41), 1.98 (q ≥ 211) | 13 to 100,003 | a drift of \|G\|/q² with q (R: larger q are cheap) | no |
| H11 | One prime of E ramified over E₀, above 2, the relative discriminant of norm 64; \|G\| even 509/512; odd part cyclic 506/512 | 13 §3.4 | supported | as stated, 512 of 512 for the ramification | 13 to 100,003 | a field with two ramified primes (R) | no |
| H12 | The polynomials of surfaces are exactly the admissible ones | 13 §3.2 | supported, **now out of sample** | 7 primes, 44 of 44; at q = 41, enumerated after the list of `sixlines_groupsize.jsonl` was written: 21 of 21, lattices = \|G\| in 6 of 6; at q = 43 (`sixlines_cm_p43.jsonl`, written for note 29): 24 of 24, lattices = \|G\| in 10 of 10 maximal orders [R] | 13 to 43 | an admissible P without surfaces (R: run at q = 41 and 43 and passed: nine primes, two out of sample) | no |
| H13 | The order is maximal only if μ = ν = 1 | 13 §3.2 | supported; no proof | 270 of 270 maximal have μ = ν = 1 (1,264 admissible) | 13 to 101 | a maximal order with μν > 1 (R; above 101 it was assumed when sampling) | no |
| H14 | One genus: the lattices of a maximal order number \|G\| | 13 §3.3, 19 Lemma 2 | proved on paper [D] given one ramified prime; measured | 215 of 215 records (191 polynomials) | 13 to 211 | a count ≠ \|G\| (S: none) | yes |
| H15 | Every character of order 2 of G is that of E(√d)/E, d ∈ E₀; at most 3 + 2-rank Cl(E₀). *Supersedes* 13 §6 "at most 3" | 19 §2, 18 finding 4 | supported; [D] Lemma 1 | 699 of 699 matched by an S-unit; 2-rank 0..4: 3/287/160/28/2 | 13 to 10,007 | a character with no d (S: none) | no |
| H16 | Surfaces per polynomial, all orders, ≈ 0.006·q² | 13 §3.3 | supported, small range | 7 of 7 averages; q = 41: 243 surfaces (469 configurations) on 21 polynomials, 0.0069·q² (`sixlines_cm_p41.jsonl`); q = 43: 304 surfaces on 24, 0.0069·q² (`sixlines_cm_p43.jsonl`) [R]. From q = 37 on the constant is nearer 0.007 | 13 to 43 | a drift of the average with q (S) | no |
| H33 | 2²⁵⁶ elements need q ≈ 2¹³¹; a surface is named by ~525 bits | 13 §5, 17 §6 | heuristic (arithmetic on H10) | — | — | H10 failing (R) | no |
| **G3** | **Computability for the honest party** | | | | | | |
| H17 | The CM method evaluates the action only for small \|G\|; the precision it needs grows with \|G\| | 20 §9, 24, 26, 29 | supported in range, as an upper bound; [H] beyond | first working digits 300 (\|G\| = 2), 1,200 (8), 1,200 to 2,000 (12), 3,000 (16), the only precisions tried at 12 (q = 29) and 16; 7 fields of note 20 skipped or failed. At \|G\| = 12 recognition needed 2,000 digits at q = 37 (1,200 failed) and worked at 2,000 at q = 41. The first q = 43 runs stalled in a diagnostic factorisation after H had been recognised, so they bear on neither (note 29 §4.2). On the reruns \|G\| = 12 needed 3,000 digits at q = 43 (2,000 failed); \|G\| = 24 needed 4,500, for a coordinate chosen to be injective modulo 43 (note 29 §7.4) | 17 to 43 | a run at \|G\| ~ 100 with low precision (R) | uses |
| H18 | δ costs about a minute per curve and needs F_{p⁴}; about half the pairs are unusable | 23, 25, 26 | supported in range | 106 of 242, 250 of 532 computed; median 74 s and 86 s | 29 to 41 | — | no |
| H34 | No step at ℓ ≠ p can be computed on the K3 side itself | 13 §2, 15 §1 | **open** (absence of a method) | — | — | a computed step (N) | no |
| H35 | On the fourfold a step costs ~ℓ⁴ by theta functions; the Prym variety has no theta coordinates | 13 §5, 15 §5.1, 17 §4.2 | heuristic [M: Lubicz–Robert, abstract only] | — | — | (N) | no |
| H36 | The zeta function at q ≈ 2¹³¹ prime is out of reach (~√p); instances need small characteristic | 13 §5, 17 §5.2 | heuristic [M; P for Costa–Harvey–Kedlaya in 17] | — | — | (N) | no |
| H37 | Given the action: classical attack ~√\|G\|, quantum Kuperberg 2^{O(√log\|G\|)} | 13 §5 | heuristic [M] | — | — | (N) | no |
| **G4** | **The decision problem** | | | | | | |
| H19 | The ten square classes (type) and the branch count are constant on a torsor; a function of (q mod 16, e₁ mod 16) | 19 §4.2 | supported; [D] by H14, conditional on the ten classes being local data at 2, which note 19 §7 does not prove | 36 of 36 torsors, 419 surfaces; 22 combinations, 0 with two types; the rule in (q mod 16, e₁ mod 16) tested on 10 of the 22 [R]; and constant on all 19 labelled torsors of note 28, including the two at q = 41 (28 surfaces) that note 19 did not have | 13 to 49 | a torsor with two types (S: none) | yes |
| H20 | Conic, j2, j4, quartic product: not characters (neither constant nor balanced) | 19 §4.3 | supported (negative) | neither on 15/23, 19/34, 24/34, 14/23 torsors; against labels (note 28), neither characters nor functions of the label difference on the 9 torsors with full labels (p 0.16 to 1), and no separation of the cosets of ker χ₋₁ on the 5 coset torsors (p 0.18 to 1) | 17 to 49 | — (S: a candidate matching the labels; the audit) | yes |
| H21 | Square classes of J₂, J₂J₆, J₃J₅ are not characters of the labels | 20 §7.6 | supported (negative) | J₂J₆: 6 of 8 torsors not, 1 character; J₃J₅: 5 not, 1 constant; J₂ not recomputable from the stored ratios. Computed in note 28 from the J-vector of each configuration, on its 9 torsors: J₂ a character on none, J₂J₆ on one (the Z/4 of q = 23, e₁ = 10, as note 20 found), J₃J₅ on none | 23 to 37 | (S) | yes |
| H22 | χ₋₁ is a character of G for 80 of 477 polynomials, exactly when K⁺ has a dyadic prime of odd degree; never for q ≡ 3 mod 4; d = ±2 never | 19 §2.4–2.5 | supported; [D] sketch | 80 of 477; 477 of 477 for the criterion; 53/126 (1 mod 8), 27/59 (5 mod 8), 0/292 | 13 to 10,007 | (S: none) | yes |
| H23 | δ is an invariant of the surface and balanced on each torsor | 23 §5.4, 26 §6 | supported | 60 of 60 surfaces, 356 pairs; as an invariant of the surface tested on 54 of the 60, since on 6 of note 23's 32 every computed pair comes from one configuration [R]; chance of agreement by coin tosses 2⁻⁷⁴ and 2⁻²²²; balance by chance 1.4 % (three torsors); reproduced as the positive control of note 28 | 29, 37, 41 | a surface whose pairs disagree (R: the exponent-64 curves need F_{p⁸}) | yes |
| H24 | δ(g·X) = χ₋₁(g)·δ(X) on the five torsors computed | 24, 26 | supported | 60 of 60; chance for a balanced function 1 in 7,470,540 and 1 in 2,972,970, together 1 in 2.2·10¹³ (for Z/4 × Z/2 only 1 in 3 against the two other characters); in note 28 δ is recovered as the positive control, first of 6,579 candidates on the coset torsors, exact p 4.5·10⁻¹⁴ | 29, 37, 41 | a surface on the wrong coset (R: q = 73, 89, 97, F₄₉) | yes |
| H25 | The cosets come out the same by the orbits of ideals, the factors of H over E(i), and ξ mod 4 on the lattices; Frobenius at ℓ ≡ 3 mod 4 has cycles of the order of its class | 25 §4, 26 §5 | supported | 5 of 5 torsors; 288 primes ℓ ≡ 3 mod 4 (160 + 128) and 200 with ℓ ≡ 1 (120 + 80), 488 in all, 0 exceptions | 29 to 41 | (S: none) | uses |
| H26 | δ = ζ·Γ(r)/H′(r): the sign of δ is fixed by the lattice | 26 §7.1 | heuristic: 5 of 5, rule formed after 4, so tested at 1 in 2 (1 in 16 before looking) | 60 of 60 surfaces | 29 to 41 | a sixth torsor against it (R) | yes |
| H38 | For every polynomial with χ₋₁ the decision problem leaks one bit (attack of Castryck–Sotáková–Vercauteren, *not* on CSIDH: 25 finding 5) | 24 §5, 25 | **open** beyond 5 of 80 polynomials; on the five, [E] (H24) with the inference of note 24 §2 | 5 polynomials of 80; none with q ≡ 1 mod 16 or p > 41 | 29 to 41 | δ failing on one torsor (R) | yes |
| H39 | The characters of order 2 of G are at most 3 genus characters, a number bounded whatever q, and none of them leaks | 13 §6 | **contradicted**: the bound by 18 finding 4 and H15; that none leaks by H24 | — | — | — | — |
| H40 | The method of Castryck–Sotáková–Vercauteren has nothing to read on the K3 lattice: no local invariant of the lattice is a character | 19 §2.2 | proved on paper [D] from H14 | consistent with H19 to H21; in note 28 no residue symbol of order 2 of the invariant monomials in J₂, …, J₆, or in the power sums of the ten products and the conic (up to 270 and 2,110 per field, orders 2 and 4), is a character; an order-4 reading could not reach significance there, since quartic symbols exist on only 3 of its 9 torsors | — | a local invariant of T that is a character (S/R: the audit) | yes |
| H41 | Over F₂₅, F₄₉ an invariant fixed by Frobenius conjugation cannot be a character that is −1 on 𝔭; such a character exists for 9 of 12 polynomials | 19 §2.3, §4.4 | supported; [D] criterion | 9 of 12; candidates fixed by conjugation 37/37, 1,080/1,080 (note 19's count; not recounted by the ledger, recounted from `sixlines_characters.jsonl` [R]) | 25, 49 | (S) | yes |
| H42 | The other characters of order 2 (the 397 of 477 polynomials with no χ₋₁, and the characters other than χ₋₁ of the 80 with it) have no reading from a surface | 19 §5, 24 §5 | **open** (absence) | the negatives of H20, H21; in note 28, within its families no invariant reads any of the characters of order 2 of the 9 fully labelled groups (none has χ₋₁): Holm-corrected p = 1 for all 4,838 candidates; planted exact readings on every torsor are found 200 times in 200, readings on 4 of the 9 torsors 37 times in 200; characters of order 4 were tested there only through their squares | 13 to 49 | a surface invariant matching one (S/R: the audit) | yes |
| H43 | The class of g cannot be computed from two surfaces in time polynomial in log q (search problem; "no easier" half) | 13 §7, 14 §8 | **open**: "the absence of one is not evidence" (13 §7, of an attack) | in note 28 no pair function f(X, Y) in its families is a function of the label difference beyond chance, on groups of order at most 8; in note 29 none separates the cosets of a character of order 3, on 7 torsors with \|G\| ≤ 12; the label remains a function of the complete invariants through the class polynomial, which such a test cannot see and which costs a class polynomial of degree \|G\| | — | an invariant computable in time polynomial in log q that determines labels beyond cosets (S/R: the audits) | yes |
| **G5** | **Below Picard number 16** | | | | | | |
| H27 | The two Jacobians of a (4, 4) double cover have the local structure at 2 of a 2-step: (t − 1)² ∣ P mod 2, a split pair when exactly twice, \|Br\| = 4·square | 22 §4.2 | supported; controls differ | 3,613 of 3,613; 1,349 of 1,349 split; 7,415 of 7,415 Artin–Tate numbers; controls: 381 of 492 without (t − 1), 17 of 49 inert, 298 of 843 four times a square | p = 5, 7; ρ̄ = 12, 10, 2 | a pair without it (S: none) | no |
| H28 | Under condition C, X₁ ≇ X₂ and the class of 𝔩/𝔩̄ is never trivial; \|G\| ≈ 0.08p⁴ (ρ̄ = 12), 0.065p⁵ (ρ̄ = 10). *Not tested:* that the step is that class | 22 §4.3 | supported (the classes), open (the step) | 972 of 972; Σ1/\|G\| = 23.5 expected isomorphic, 0 seen; ∏(1 − 1/\|G\|) = 2·10⁻¹⁰ for e8a2 at p = 5 | 5, 7 | an isomorphic pair (S: none); the step (R: a closed cycle, ~4·10⁵ surfaces) | no |
| H44 | No reduction to abelian varieties is known for any full family with ρ̄ ≤ 14 | 17 §2 | **open** (literature as of 2026-11-14 [P]) | — | — | a construction (N) | no |
| H45 | Below 16 there is no usable action: no step at an odd prime, half a step at ℓ = 2 | 17, 22 | **open** (absence) | H27, H28 | — | (N) | no |
| H46 | \|G\| ≈ q^{(20−ρ̄)/2}; an instance has ~512 bits at every ρ̄ | 13 §3.3, 17 §6 | heuristic; the exponent measured over a range of q only at ρ̄ = 16; at ρ̄ = 12 two primes give a slope of 3.8, at ρ̄ = 10 one prime gives only a constant [R] | H10, H28 | ≤ 7 below 16 | (R) | no |
| H47 | The Kuga–Satake variety exists over F_q (Madapusi Pera [P]); its simple factor has dimension 2^{(20−ρ̄)/2} | 17 §3–4 | heuristic [D, not checked against a source] | — | — | (N) | no |

H1–H28 are the IDs of the sections of `sixlines_ledger.py` and its log. H29–H47 are the hypotheses that have no section in the script, placed in their groups.

The rows H19 to H21, H23, H24, H40, H42 and H43 include the results of the audit of §5, carried out in [note 28](28-leakage-audit.md) with the tests fixed before the join; H42 and H43 stay open. H14, H22, H26 and H38 are not affected by it. The rows H7, H12, H16, H17 and H43 include the data and results of the test of characters of order 3 in [note 29](29-characters-of-order-3.md). The ledger has no row for characters of order above 2 as such; note 29 is the test of order 3, and H43 is narrowed by it again.

---

## 3. By group

### 3.1 G1: The reduction to abelian fourfolds

- **H1 [E].** The identity was checked by integer arithmetic on all 1,606 distinct polynomials stored in `sixlines_groupsize.jsonl`, q from 13 to 100,003, and on the 44 of `sixlines_cm.jsonl`. For each, μ and ν were recovered from P(±q) and the second exterior power of χ was formed. There were no exceptions. The 23,945 polynomials of note 13 §4.2 are not in a file. The reason the check adds little: note 18 re-derived the identity [R], so the data only guard against a slip in the formula; and 1,562 of the 1,606 were built from (e₁, μ, ν) by the same equations that Proposition Y restates, so they cannot fail, and only the 44 of enumerated surfaces test anything [R]. The identity holds for every q. It is one of the few statements in the ledger that reaches cryptographic size, and it says nothing about surfaces.
- **H2, H3 [E].** These are recounted from `sixlines_prym.jsonl`, `sixlines_prym_q25.jsonl`, `sixlines_prym_q49.jsonl` and `sixlines_prym_p41.jsonl`. Every column of the table of note 14 §4 matches. So do note 16's pair counts, note 26's 1,507/763/744/173, and note 18's 809 degenerate pairs and 0 inverting ones. The one mismatch is the noun "surfaces" (§1 item 3). Reading: the agreement is between Frobenius polynomials, which are isogeny invariants. It shows the curve has the right zeta function, not that its Prym variety is the right one inside the isogeny class.
- **H4.** The per-pair data of note 21's `types` run are in `sixlines_forms_fp.log` only. The log lines on sign + and irreducible P (534 configurations, 0 with a collinear triple, 14,394 pairs of the type of the sign) are printed by the ledger, not recounted.
- **H5 [E].** Every count of note 15 §3.2 matches `sixlines_dictionary.jsonl`, and the 46 of note 18 matches `sixlines_review.log`. The paragraph "Why" of note 15 is corrected by note 18, and note 20 §5 shows that the kernel and cokernel compare objects that are not surfaces. What stands is the isomorphism of the odd parts. Like H1, it concerns class groups only.
- **H7, H8 [E].** The rule of note 20 is applied: the last record of each field and command in `sixlines_cmmethod.jsonl` (51 records, 39 kept). Under it every total of note 20 §§6 and 7 matches. With the first record of each field instead, 74 of 87 enumerated points are matched. The difference is the first sweep, which selected surfaces by e₁ alone and so mixed in a second polynomial at q = 31 and 37 (note 20, Reproduce). So the rule matters, and the files bear it out. The 19 fields, with 126 surfaces, come from `sixlines_cmmethod.jsonl` (15), `sixlines_delta_labels.jsonl` (2 new fields; q = 37, e₁ = 18 again) and `sixlines_delta41_cosets.jsonl` (2).
- **H9 [E].** `sixlines_jacobian_levels.jsonl` gives 422 pairs and no disagreement inside a surface. With the association recorded in `sixlines_p29.jsonl` and `sixlines_p37.jsonl`, 5 of the 7 two-block surfaces are self-associated, as note 25 found. So in two of the three polynomials the count reduces to "the object with extra automorphisms is on the upper level".
- **H29, H30, H32.** These carry the reduction. No file bears on them directly. The indirect evidence is that δ, an invariant of the Prym variety of note 14 that is not an isogeny invariant, splits 60 surfaces exactly as the CM action does (H24). If the map were not equivariant on these torsors, that agreement would have no reason. It is evidence on five torsors with |G| ≤ 16, not a test of a step.

### 3.2 G2: Size and structure of the group

- **H10 [E].** See §1 items 1 and 2. 30 of the 31 rows of the table of note 13 §3.3 match in every column, and the row for q = 17 is now corrected there. The fits match. The median and quartiles that had been given for the 512 are those of the 486 groups with q ≥ 41, the population of the first fit; note 13 now states the values on the 512. The law in range is well supported. Beyond q = 100,003 it is the class-number heuristic of note 13 [H], which is standard for CM fields. q ≈ 2¹³¹ is 2¹¹⁴ times past the data.
- **H11 [E].** All counts of note 13 §3.4 match. Note 18's 212/212 and 258/258 also match, but the two populations add to 470 of 512. The other 42 groups have h(E₀) even and at most one even cyclic factor. Note 18's statement ("more than one only when h(E₀) is even") is true of those 42 as well, so the 42 are no exception.
- **H12 [E], new.** `sixlines_groupsize.jsonl` was written on 2026-11-12. It lists 21 admissible polynomials at q = 41 (its `all` run), and the enumeration at p = 41 of note 26 (`sixlines_cm_p41.jsonl`, 2026-11-15) finds exactly those 21. For the 6 with maximal order, the lattice counts equal the |G| of note 13 (12, 16, 26, 12, 10, 4). This is a prediction made three days before the test, by a different script, and no note had compared the two. At p = 43, enumerated for [note 29](29-characters-of-order-3.md) (`sixlines_p43.jsonl`, `sixlines_cm_p43.jsonl`), the 24 polynomials found are exactly the 24 admissible ones of the same list, and the lattice counts equal |G| for the 10 with maximal order [R]. Reading: the admissible polynomials are exactly those that occur, now at nine primes, two of them out of sample. The falsifier of H12 has been run and passed.
- **H13.** 270 of 270 for q ≤ 101. Above 101 the sampling drew only μ = ν = 1, so the data cannot test it there.
- **H14 [E].** 215 records, which are 191 distinct polynomials once the excluded q = 125 is left out (note 13's 191). Lemma 2 of note 19 makes this a consequence of H11 [D]; that lemma was derived once and has not been re-derived.
- **H15, H16 [E].** All counts match `sixlines_characters_genus.jsonl` and `sixlines_cm.jsonl`. For H16 the 469 at q = 41 that this row first quoted are configurations; `sixlines_cm_p41.jsonl` gives 243 surfaces on 21 polynomials, 0.0069·q², and `sixlines_cm_p43.jsonl` 304 on 24, also 0.0069·q² [R].

### 3.3 G3: Computability for the honest party

- **H17 [E].** These are the first working precisions at which the polynomials were recognised. They come from the first record of each field in `sixlines_cmmethod.jsonl` and from the records of notes 24 and 26:

  | \|G\| | 2 | 4 | 6 | 8 | 12 | 16 |
  |---|---|---|---|---|---|---|
  | digits | 300 (5 fields) | 300, 700, 700 | 700, 700 | 1,200 (5) | 1,200, 2,000, 2,000 | 3,000 |

  The 1,200 digits at |G| = 12 (q = 29) and the 3,000 at |G| = 16 were the only precisions tried there, so they are upper bounds [R]. The growth is roughly linear in |G| over this range: about 190 digits per element. The method needs a class polynomial of degree |G|. Seven fields of note 20 were skipped (10 to 26 lattices per type) or failed; note 24 did two of them (q = 29, e₁ = 2 and q = 37, e₁ = 2) with parallel period points, and a rerun of note 20 one. [Note 29](29-characters-of-order-3.md) adds two fields at |G| = 12, q = 37, e₁ = 6 (2,000 digits; 1,200 failed) and q = 41, e₁ = 2 (2,000), and the first runs at q = 43, which do not bear on this ledger entry: they stalled in the full factorisation of the lcm of the denominators of H for a progress line, after H had been recognised (note 29 §4.2). So they say nothing about the cost of recognition or about the precision. The reruns (note 29 §7.4) needed 3,000 digits for \|G\| = 12 and for Z/6 × Z/2 at q = 43, where 2,000 failed. For \|G\| = 24, 3,000 digits gave the class polynomial of one coordinate chosen to be injective modulo 43, and 4,500 the polynomials of the other invariants. So at \|G\| = 12 the first working precision is 1,200 to 3,000 digits across four primes, and the count of surfaces against q also limits which coordinates work. This is the honest party's only working tool for the action. At |G| ≈ 2²⁵⁶ it does not run, and that is an argument, not a measurement.
- **H18 [E].** δ takes 74 seconds (median, p = 29, 37) and 86 seconds (p = 41) per curve. 106 of 242 and 250 of 532 curves give a value; the rest need a pairing of higher order. There are no data above p = 41.
- **H34 to H37.** No file bears on these. They are what makes the problem a candidate at all: nobody can step on the K3 side, the fourfold side needs theta coordinates that the Prym variety does not have, the zeta function at a large prime is out of reach, and the attacks are generic. All are [M] or absence of method.

### 3.4 G4: The decision problem

- **H19, H20 [E].** Surfaces were rebuilt from configurations through the `assoc` field of `sixlines_p{p}.jsonl`, `sixlines_q25.jsonl` and `sixlines_q49.jsonl`. Then every number of note 19 §4.2 and §4.3 comes out: 36 torsors, 419 surfaces, |G| surfaces in each, the four rows of constant/balanced/neither, and c₁ taking two or more values on 263 configurations (137 with ±c₁, 131 with c₁ and its conjugate).
- **H21 [E].** The square classes of J₂J₆ and J₃J₅ equal those of the stored ratios J₆/J₂³ and J₃J₅/J₂⁴, since they differ by J₂⁴. Against the stored labels of note 20, with the field charted by J₃ left out as there, they give 8 torsors and the counts of note 20 §7.6 exactly. The row for J₂ alone cannot be recomputed from the stored files, because only ratios are stored; note 28 computes it from the J-vector of each configuration (H21 in §2).
- **H22 [E].** All counts match, including the dyadic criterion (`sixlines_characters_dyadic.log`, 477 of 477) and the 19,564 primes.
- **H23, H24 [E].** The ledger recounts from `sixlines_jacobian_character.jsonl`, `_extra.jsonl`, `sixlines_delta41.jsonl`, `sixlines_delta_labels.jsonl` and `sixlines_delta41_compare.json`:
  - 242/106/126/9/1 and 532/250/269/8/3/2;
  - δ constant on 32 of 32 and 28 of 28 surfaces;
  - u constant on 10 and on 0 of them;
  - note 23's first-run counts 11, 12, 15 and 29;
  - the controls, 36 of 36 and 83 of 83;
  - the agreement with the cosets, 32 of 32 and 28 of 28.

  The chances quoted are recomputed exactly: 1.4 %, 2⁻⁷⁴, 1/7,470,540 = (2/924)²·(2/70), and 1/2,972,970 = (2/12,870)·(2/924). For all five torsors together a balanced function agrees by chance with probability 4.5·10⁻¹⁴. Two things that probability does not cover. It assumes δ is balanced and otherwise arbitrary; for Z/4 × Z/2 another character is an alternative at odds of 1 in 3 (note 24 §6). And it says nothing about the 75 polynomials not tested.
- **H25 [E].** Recounted from `sixlines_review2.jsonl` and the `galois` records of `sixlines_delta41_cosets.jsonl`: 160 + 128 primes with ℓ ≡ 3 mod 4 and 120 + 80 with ℓ ≡ 1, no exception.
- **H26.** The log `sixlines_delta41_compare.log` gives 5 of 5. The rule was formed after four torsors had been seen, as note 26 says, so the fifth is the only out-of-sample point: one chance in two.
- **H38, H39, H40, H42, H43.** These are statements about what can be computed. The stored data reach them only through the negative results above and the positive one of H24. The audit of [note 28](28-leakage-audit.md) is the one test the stored data allow: any stored invariant of a configuration against the labels. §5 lists what can be joined. [Note 29](29-characters-of-order-3.md) extends it to characters of order 3.

### 3.5 G5: Below Picard number 16

- **H27, H28 [E].** Every count of note 22 §§4.1 to 4.4 that the files hold matches, including the controls, Σ1/|G| = 21.54 + 1.64 + 0.30 = 23.48, the product ∏(1 − 1/|G|) = 2.0·10⁻¹⁰ (a Poisson estimate e^{−21.5} would give 4.4·10⁻¹⁰), the medians 0.083, 0.078 and 0.065, and the partner counts of `ellk3_e8a2_graph.log`.
- **H44 to H47.** These are literature, absence of method, and heuristics. No file bears on them beyond the three small measurements of the group size.

---

## 4. What the data cannot reach

All the surfaces in this repository are over fields with at most 49 elements. All the groups have q ≤ 100,003. All the torsors on which the action was evaluated have |G| ≤ 16. A hard instance would have q ≈ 2¹³¹ and |G| ≈ 2²⁵⁶. For each hypothesis, how its evidence can bear on that size:

| How it bears at cryptographic size | Hypotheses |
|---|---|
| **An identity or proof valid for all q.** The data only check the algebra | H1, H6; H14 and H40 where exactly one prime ramifies (H11); the [D] parts of H4 and H15, and the [D] sketches of H5 and H22 |
| **A measured law, extrapolated.** Valid in range, the extrapolation is a heuristic | H10, H11, H13, H16, H17, H18, H46 |
| **Small instances of a structural statement.** Each passed test makes an error in the construction less likely. None measures difficulty | H2, H3, H7, H8, H9, H12, H19–H28, H38, H41 |
| **Argument only.** No file bears on it directly; H29 has the indirect evidence of H23 and H24 (§3.1) | H29, H30, H31, H32, H33, H35, H36, H37, H44, H47 |
| **Absence of a method.** Neither data nor argument can confirm it | H34, H42, H43, H45 |

H39 is in no row: it is contradicted (§2) and bears on nothing.

So the hardness verdict comes from the bottom three rows. What the data establish is that the objects are what the arguments say they are: the curve has the right zeta function, the group has the predicted order and structure, the CM points are the surfaces, and one character is computable. They establish nothing about the cost of any computation at a size where cost would matter. The one leak found (H24) is a structural fact. It would hold at every size if H38 holds, and that is not tested beyond q = 41.

---

## 5. For the leakage audit

Stored data that can be joined without a rerun (`AUDIT` section of the log):

- **Full labels in G.** 67 surfaces (134 configurations) in 15 fields of note 20 §7.5: 14 fields fully labelled, with 66 surfaces, and the one labelled surface of q = 37, e₁ = 18 [R]; with the five invariants J₄/J₂², J₆/J₂³, J₃²/J₂³, J₃J₅/J₂⁴, J₅²/J₂⁵ modulo q.
- **Cosets of ker χ₋₁ only.** The 60 surfaces of notes 24 and 26 (118 configurations not already labelled).
- **Joined to the invariants of `sixlines_characters.jsonl`** (type, conic, j2, j4, quartic, quartic type, orientation, branch, c₁ values): all 134 labelled configurations and 62 of those with cosets only. The q = 41 configurations have no record there; their invariants would need a rerun of `sixlines_characters.py surfaces` on `sixlines_p41.jsonl`.
- **Joined to δ, u, M and the kernel tallies**, per pair of lines: the files of notes 23 and 26.

The audit bears on H14, H19–H24, H26, H38, H40 to H43, and it uses H7, H8, H17, H25 and H29 to H31. An invariant that tracked a full label beyond the coset of ker χ₋₁ would refute H42. If it followed a label exactly, it would refute H43. One that tracked the coset where χ₋₁ does not exist would contradict H22. With torsors of at most 16 surfaces, a chance agreement is likely: one balanced function matches one given order-2 character, up to sign, with probability 2/C(n, n/2) on a torsor of n surfaces, 1/3, 1/35 and 1/6,435 for n = 4, 8 and 16, which is smaller than 1/|G| [R]. So the audit needs exact nulls and a correction for the number of candidates tried, as note 19 §4.3 did with its 4,194,303 functions; note 28 uses both.

---

## 6. What this does not give

- A check of anything beyond the stored files. The two counts that have no file (§1) are not confirmed, and the reruns that would confirm them were not run. The separate re-derivation check of this note recounted the stored files only: it did not cover the 23,945 polynomials and the 234 odd lattices, the per-pair data behind the 14,394 pairs, whether the stored files are right, any re-derivation of the [D] statements (Lemma 2 of note 19, Lemma 9 of note 13, the sketch of H22), or the [M] and [P] sources.
- A check of the scripts that wrote the files. A count that matches its file shows that the prose reports the file, not that the file is right.
- Any new evidence on hardness. The one new comparison, H12 at q = 41 and 43, is about which polynomials occur.
- A judgement of the [D] and [M] statements on which the classification rests. "Proved on paper" means a derivation exists in this repository and, for H1 and H6, that note 18 re-derived it. H14 and H40 rest on the one-genus step, Lemma 2 of note 19, which was derived once and not re-derived; note 18 re-derived Lemma 9 only. Beyond those re-derivations nobody has checked it.
- The audit itself. It is [note 28](28-leakage-audit.md): no further leak was found within the families tested, and δ is the only survivor on the coset torsors. Characters of order 3 are tested in [note 29](29-characters-of-order-3.md), with the same result.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_ledger.py > sixlines_ledger.log
```

It reads only stored files, takes about a second on one process, and is deterministic (two runs gave identical logs). Its last lines list every file read and every mismatch. The only computation beyond counting is integer arithmetic on stored polynomial coefficients, and one call of PARI's `polisirreducible` on the 44 stored degree-8 polynomials χχ̄.

Reruns that would supply the counts with no file:

```bash
python -I sixlines_groupsize.py --fourfold 25 43 47 49 53 59 67 71 73 79 83 89 97 103 107 109 113 121 125 127 131 137 139 149 151 157 163 167 169 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257
python -I sixlines_forms_fp.py types 5 7 11 13 17 19 23 29 31 37 --workers 4
```

The count of 234 odd lattices of note 20 §5.1 has no command. It is a loop over `CMField(q, e1, mu, nu).lattices()` for the 21 fields.

## Sources

No new source was read. In this repository: notes [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md), [17](17-below-picard-16.md), [18](18-checks-of-notes-13-to-16.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [21](21-fourfold-forms-over-fp.md), [22](22-two-neighbours-of-elliptic-k3.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [26](26-delta-is-the-character-at-q-41.md), [28](28-leakage-audit.md), [29](29-characters-of-order-3.md). Stored files: the 42 listed at the end of `experiments/sixlines_ledger.log`.
