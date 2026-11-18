# 26. δ is the character at q = 41: a test out of sample

*A test out of sample of the claim of [note 24](24-the-form-mod-4-is-the-character.md) that the sign δ of [note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) §5 is the character χ₋₁. Note 24 finds it on three torsors, at q = 29 and 37, both 5 mod 8; [note 25](25-checks-of-notes-20-23-24.md) §7 names the nearest test out of sample, q = 41 ≡ 1 mod 8. This note does that test: the enumeration at p = 41, the curves of [note 14](14-fourfold-of-a-six-line-surface.md), the CM method of [note 20](20-cm-method-and-the-action.md) at 12 and 16 lattices per type, and δ with a pairing of order 32. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_delta41.py` or by the script named; **[H]** heuristic; **[R]** re-derived in a separate re-derivation check. Statements of other notes keep the label they have there. Checked once by a separate re-derivation from the stored files, with its own code: δ and u recomputed from all 783 records, the class polynomial for e₁ = 22 rebuilt with independent recognition, G and χ₋₁ rebuilt. The claim of §1 item 1 stands and no number was found wrong; what the check confirmed is marked [R], its corrections are made in the text, and what it did not cover is in §8. Notes 20, 23 and 24, on which this note rests, have been checked once, in note 25.*

Notation as in notes 23 and 24: P the Frobenius polynomial of the surface, with coefficients e = (e₁, e₂, e₃); E = Q[F] of degree 6; G the group that acts simply transitively on the set S_P of surfaces with polynomial P; 𝔭 = (q, F − e₁); χ₋₁ the character of G of the extension E(i)/E; C the curve of genus 5 of note 14 for an ordered pair of the six lines, C_tw its quartic twists, J its automorphism of order 4, π Frobenius.

---

## 1. Result

1. **The prediction of note 24 §2 holds at q = 41, on both polynomials with the character [E].** The surfaces of each polynomial split into the two cosets of ker χ₋₁, and δ is constant on each coset with opposite values:

   | q | e | G | surfaces | with a δ | δ the same on all pairs of the surface | cosets | δ on coset A, on coset B | surfaces that agree |
   |---|---|---|---|---|---|---|---|---|
   | 41 | (6, −369, −20172) | Z/16 | 16 | 16 | 16 of 16 | 8 and 8 | −1, +1 | 16 of 16 |
   | 41 | (22, −369, −73964) | Z/12 | 12 | 12 | 12 of 12 | 6 and 6 | −1, +1 | 12 of 12 |
   | both | | | 28 | 28 | 28 of 28 | | | 28 of 28 |

   With note 24 that is 60 surfaces of 60 in five torsors at three primes, and the first with q ≡ 1 mod 8. Both sides were reproduced independently: δ from the stored records of all 28 surfaces, and the split for e₁ = 22 from a class polynomial rebuilt and recognised with independent code (§5) [R].
2. **The definition was fixed first and not changed [D].** δ is given at q ≡ 1 mod 8 by the formula of note 23 §5.2 without change. What changes is how its factor u is read: over F_{p⁴} the Tate pairing has order 32 and not 16, the eigenvalue π₁ has π₁⁴ − 1 of valuation 8 at 1 + i and not 7, and the number M of note 23 §5.3 is a class modulo 16 inside Z/32 where it was one number modulo 16 (§2). The definition and the reading rule were written into the docstring of the script before any curve at p = 41 was computed, and every run that produced a value of δ was started before any coset had been seen (§3). Two controls planned there return nothing, for a reason that could have been seen on paper (item 9).
3. **The surfaces at p = 41 [E].** 1,975,506 configurations of six lines, 2,948 classes, 1,234 of them self-associated; 17 minutes 20 seconds on 14 processes. The two polynomials have 32 and 24 configurations, none self-associated: 16 and 12 surfaces, the numbers of lattices and the orders of G that note 19 lists. The curve of note 14 was computed for all 30 ordered pairs of each of the 56 configurations: 1,507 computed, 173 degenerate, and in all 1,507 the second exterior power of the quartic is P or P(−t).
4. **The cosets, by the CM method and two further routes [E].** The CM method of note 20 works at 12 and at 16 lattices per type, at 2,000 and 3,000 digits. In both fields the class polynomial modulo 𝔭 has exactly the enumerated surfaces as its roots (12 of 12, 16 of 16). The orbits of the 48 ideals above eight primes ℓ ≡ 1 mod 4 (note 24), the two factors of the class polynomial over E(i) (note 25 §4.3) and the class modulo 4 of the polarisation ξ of each lattice (note 25 §4.4) give the same split, in both fields. The indicator polynomial of note 24 was recognised, so its route was used; the other two were run as well. The three routes check the code and the recognition; they are not three independent determinations. Once an indicator is recognised as i times a polynomial over E, with values that square to −1, its partition is necessarily the partition by the factors over E(i) (as note 25 §4.5 says of note 24) [D]. What the orbit route adds is that the 48 ideals act on the lattices with exactly those orbits.
5. **δ [E].** 532 curves tried over F_{41⁴}, 250 computed, on 5 to 23 pairs of lines for each surface of the first polynomial and 4 to 21 for each of the second. On every surface all computed pairs give the same δ. The factor u alone is the same on all pairs of no surface of the 28; with half of the normalisation it is constant on 3 and on none. δ is +1 on 8 surfaces and −1 on 8, and +1 on 6 and −1 on 6. The twists C₁ and C₂ and other random classes give the same u on 83 of 83 controls; over F_{41⁴} the twists are isomorphic to C₀, so these controls check bookkeeping only (§6).
6. **δ is the form modulo 4 of the lattice, here too [E].** The sign that note 25 reads on the lattice (L, ξ) of the CM lift, the unit part modulo 4 of −iξ at the prime of degree 1 above 2, splits the surfaces as δ does: 28 of 28. Its indicator polynomial is i times a polynomial over E, recognised in exactly one of the eight ways and with no sign between the types, in both fields.
7. **The sign as well, on all five torsors [E, found after the fact].** With ζ the fourth root of unity in F_q by which the script lets J act, δ(X) = ζ·Γ(r)/H′(r) for the indicator Γ of the lattice form: δ equals the sign of the lattice at the prime of E(i) above 𝔭 where i goes to ζ, with the same rule in the two torsors at q = 41 and the three of note 24 (§7.1). The rule was formed after four of the five torsors had been seen, so only the fifth tests it: one chance in 2 if the signs were coin tosses (one in 16 is the figure before looking).
8. **By chance? [D]** A function that takes each value on half of a torsor agrees with a given split into halves, up to sign, with probability 2/C(16, 8) = 1/6,435 on 16 surfaces and 2/C(12, 6) = 1/462 on 12: 1 in 2,972,970 for the two. For Z/16 and Z/12 there is no other character of order 2.
9. **What did not work [E].** (a) Over F_{41⁴}, on 269 of the 532 curves the twelve classes (1 − J)y_a are killed by 64 and not by 32 (§2.2; this is not the exponent of the part with J² = −1, since (1 − J)y_a need not lie in that part), and the method does not apply to them, as for the 126 curves of note 23 with exponent above 16. These are the curves called "of exponent 64" below. (b) The control with N = 16 over F_{41⁴} cannot be run as I defined it: the classes are killed by 32, not by 16 (28 runs of 28). (c) The extension to the curves of (a), multiplying the classes by 2, gives no reading on any of 111 curves: doubling the point multiplies M by 4, and M is then 0 or 8 modulo 16. (d) Over F_{41²} with N = 16 only vectors with M = 2·(odd) can be read (the rule would also read odd M, but no stored vector has one), and only 4 curves of 28 have one; those 4 agree with the primary reading. (e) The pairing over F_p gives nothing at all (§2.4). None of these touches the 250 curves of item 5.
10. **What this is not.** Not a proof: two more torsors and one more prime. One prime with q ≡ 1 mod 8, and it is 9 mod 16; of the 53 polynomials of note 19 with the character and q ≡ 1 mod 8, two are tested. No step at a prime ℓ ≡ 3 mod 4 was applied to a lattice or a surface, as in note 24. And the two sides are the two sides of note 24: if one of them is wrong in a way that survives the checks of notes 24 and 25, it is wrong here too.

So the statement of note 24 §5 stands and now covers the kind of prime that note 25 flagged as untested. Nothing here changes the conclusion of note 13 or the sizes of its §3.3.

---

## 2. What changes at q ≡ 1 mod 8 [D]

This section was worked out before any curve at p = 41 was computed, and its content is the docstring of `sixlines_delta41.py`. The separate check reproduced its valuations, k = min(n, ⌈s/2⌉) and N/ϖ₁ modulo 32 [R], and found two derivations imprecise, though they give the right numbers; they are corrected in §2.2 and §2.5. The rest of the section is not re-derived. The valuations quoted are those of the command `dyadic`, which uses the quartic χ and nothing else.

### 2.1 The invariant and its normalisation do not change

u is the class of h(m) = E(m, i·m) on the line V₁ of the 2-adic Tate module, in Q₂^× modulo norms from Q₂(i) and powers of 2 (note 23 §5.2). It is defined over the algebraic closure and does not see q modulo 8. The normalisation

> δ = (−1)^{(u−1)/2}·sign(Re c₁·Im c₁)·sign(c₂)

uses two facts only: u changes sign with the conjugate structure (c₁ ↦ c̄₁), and u is the same for the four twists (c₁ ↦ i^tw·c₁, c₂ ↦ (−1)^tw·c₂). Neither depends on q. It needs c₁ = 4(±1 ± i), that is (μ, ν) = (1, 1) in the quartic of note 13, Proposition Y; for both polynomials at q = 41 the formula of note 13 gives (μ, ν) = (1, 1), as at q = 29 and 37. So δ is defined by the same formula, with (c₁, c₂) the coefficients of the quartic of the curve C₀ as `sixlines_prym.py` gives them.

### 2.2 Which 2-power torsion is rational over which extension

Let π₁ in Z₂[i] be the root of χ on V₁ and s_r the valuation of π₁^r − 1 at 1 + i. If T₁ is the part of the Tate module in V₁, the points of V₁/T₁ that are rational over F_{p^r} form Z[i]/(1 + i)^{s_r}.

The prime of degree 1 above 2 is its own conjugate, so π₁·π̄₁ = p. A unit a + bi of Z₂[i] has norm 1 modulo 8 when its even coordinate is divisible by 4, and 5 modulo 8 otherwise. So:

- **p ≡ 5 mod 8.** π₁ = i^k·(a + 2bi) with a and b odd. For even k, π₁² − 1 ≡ 4(−1 + i) and π₁² + 1 ≡ −2 + 4i modulo 8; for odd k the two are exchanged, π₁² − 1 ≡ 2 + 4i and π₁² + 1 ≡ 4 + 4i (for example on twist 0 at q = 29, where π₁ ≡ 6 + 11i and s₂ = 2). Either way π₁⁴ − 1 has valuation 5 + 2 = 7 for every twist. The V₁-part of the rational torsion over F_{p⁴} is Z[i]/(1 + i)⁷, of exponent 16.
- **p ≡ 1 mod 8.** π₁ = i^k·(1 + 4β) with β in Z₂[i], so π₁⁴ − 1 = 16β·(a unit) and s₄ = 8 + v(β) ≥ 8. For p = 41 and both polynomials β is a unit [E], so s₄ = 8 and the V₁-part over F_{p⁴} is Z[i]/(1 + i)⁸ ≅ (Z/16)².

(s₁, s₂, s₄) for the four twists of the quartic with c₁ = 4(1 + i) and c₂ = e₁:

| q | e₁ | twist 0 | twist 1 | twist 2 | twist 3 | v₂(q^r − 1), r = 1, 2, 4 | v(χ′(π₁)), every twist |
|---|---|---|---|---|---|---|---|
| 29 | 2 | (1, 2, 7) | (2, 5, 7) | (1, 2, 7) | (3, 5, 7) | 2, 3, 4 | 6 |
| 37 | 2 | (1, 2, 7) | (3, 5, 7) | (1, 2, 7) | (2, 5, 7) | 2, 3, 4 | 6 |
| 37 | 18 | (1, 2, 7) | (2, 5, 7) | (1, 2, 7) | (3, 5, 7) | 2, 3, 4 | 6 |
| 41 | 6 | (2, 6, 8) | (1, 2, 8) | (4, 6, 8) | (1, 2, 8) | 3, 4, 5 | 6 |
| 41 | 22 | (4, 6, 8) | (1, 2, 8) | (2, 6, 8) | (1, 2, 8) | 3, 4, 5 | 6 |

At q = 41 the twist with π₁ ≡ 1 mod 4 has (4, 6, 8), the one with π₁ ≡ −1 mod 4 has (2, 6, 8), and the other two have (1, 2, 8). An independent root finder in Z₂[i] reproduces all 20 triples, v(χ′(π₁)) = 6, and the unit β at q = 41 [R].

The roots of unity in F_{41^r} have order 8, 16, 32 for r = 1, 2, 4. So over F_{41⁴} the Tate pairing T_N is available with N = 32. χ′(π₁) has valuation 6 as at q = 29 and 37, so π − π₁ is still 2 times a unit on the complement of V₁, which is what the argument of note 23 §5.3 for the complement uses.

For the group as a whole only a measurement is available [E]: the twelve classes (1 − J)y_a are killed by 32 on 261 of the 530 curves that the setup accepts and by 64 and not 32 on 269. All 82 pairs whose points a, b, c of note 14 are rational are of the first kind; the pairs with conjugate points are of both kinds (179 and 269). At q = 29 and 37 the two exponents were 16 and above 16, with the same pattern.

### 2.3 How u is read

The identity of note 23 §5.3 does not use q modulo 8. For X in the kernel of π − π₁ on the part with J² = −1, killed by N = 2ⁿ, with lift X̃ = γm/ϖ₁ + (a part in the complement), ϖ₁ = π₁^r − 1:

> Z := log T_N(X, iX) + i·log T_N(X, X) = M·conj(N/ϖ₁) mod N, M = N(γ)·h(m).

N/ϖ₁ has valuation 2n − s at 1 + i, so Z determines M modulo 2^k with k = min(n, ⌈s/2⌉):

| | s | N | conj(N/ϖ₁) | M is known modulo |
|---|---|---|---|---|
| q = 29, 37, r = 4 | 7 | 16 | (1 + i)·unit | 16: one number modulo N |
| q = 41, r = 4 | 8 | 32 | 2·unit | 16: two numbers modulo N |
| q = 41, r = 2, two of the twists | 6 | 16 | 2·unit | 8 |
| q = 41, r = 1, one twist | 4 | 8 | 2·unit | 4 |

> **Reading rule.** A vector of the kernel lattice gives u when its solutions M modulo N form exactly one class modulo 2^k, not the class of 0, with v₂(M) + 2 ≤ k; then u is the odd part of M modulo 4. Among the vectors that give u, the first with the least v₂(M) is taken, its point X is formed, and X is verified directly: its order, πX = π₁X, J²X = −X, and the two pairings of X computed directly equal the values assembled from the matrices. A curve counts as computed when all four hold.

For s = 7 and N = 16 this is the rule of `sixlines_jacobian.py`. Its code accepts a vector only when the solution modulo N is unique, which never happens at q = 41; that is the one place where the routine of note 23 had to be written again (`read_u`). The steps before it (twelve classes, four matrices of pairings, the kernel modulo N) are those of note 23 with N = 2^{v₂(p⁴ − 1)}, as its code already has it.

### 2.4 The complement, and why the pairing over F_p is not used

Note 23 §5.3: a kernel vector can have a component t/2 in the complement of V₁, t a lattice vector, and it adds (N/4)·E(t, ϖt) to the exponents. This is 0 modulo N when ϖ = π^r − 1 is divisible by 4 on the complement: for r = 4, where it is divisible by 8, and for r = 2 on the twists with π ≡ ±1 modulo 2. For r = 1 it is divisible by 2 only, and there M is known modulo 4 only, so a reading would need M odd. So over F_p nothing can be read, although the V₁-part of the rational torsion is (Z/4)² there. This is a sketch, as in note 23. What the data say about it is in §6.

### 2.5 The convention for Z

Note 23 fixed Z = M·conj(N/ϖ₁), against M·(N/ϖ₁), by the other choice having no solution. At q = 41 and r = 4, N/ϖ₁ = 2ε with ε a unit, and the two conventions have the same solutions for M exactly when M·Im ε ≡ 0 mod 8. N/ϖ₁ is 10 + 20i modulo 32 for e₁ = 6 and 18 modulo 32 for e₁ = 22, the same for all twists (automatically, since π₁⁴ does not change under the twists). For e₁ = 6, ε ≡ 5 + 10i modulo 16, so M·Im ε ≡ 0 mod 8 means 4 ∣ M: the conventions differ on every vector with M ≢ 0 mod 4, odd M included, and not only, as first written here, on those with M ≡ 2 mod 4 [D, R]. No stored vector has odd M (§6). For e₁ = 22 they coincide on every vector. Where they coincide the choice does not matter.

---

## 3. The order of events

The claim "fixed first and not changed" can be checked against the files only in part, so here is what was done, in order.

1. The valuations of §2.2 were computed from χ (no curve).
2. `sixlines_delta41.py` was written with §2 as its docstring: the primary reading (C₀, r = 4, N = 32, seed 1, twelve classes), five controls and one extension. `read_u` was run on stored curves of note 23 at p = 29 and 37 and returned the stored M, u, δ and order (6 of 6, command `selfcheck`).
3. The enumeration, the list of polynomials and the curves at p = 41. The checks of note 23 §3 on 12 curves (command `validate`).
4. One trial of the primary reading on 8 curves of two surfaces, not stored; the records of those curves in the file are from step 5.
5. The primary run: the first 6 stored pairs of each configuration, 336 curves. At the same time the CM computation for e₁ = 22 was started with `--no-delta`, which keeps δ out of it.
6. δ by surface was looked at. One surface of the first polynomial had no usable pair of 12 and four surfaces had 1 to 3. The counts were 8 and 7 with one missing, and 6 and 6. A second run was started for those five surfaces on all their stored pairs, with the controls and the extension after it, and the CM computation for e₁ = 6.
7. Only then the log of the cosets for e₁ = 22 was read and compared with the δ of step 5: 12 of 12.
8. The second run ended: the missing surface has δ = −1 on 21 pairs, the value that balance required. The CM computation for e₁ = 6 ended after that and was compared: 16 of 16.

What the files show: the docstring of the script; the records of `sixlines_delta41.jsonl` in the order they were written (the first 336 are step 5); the logs. The record order and the file times agree with steps 5 to 8 [R]: run 1 is exactly the first 6 stored pairs of each configuration, and run 2 exactly the remaining pairs of the five surfaces named; all δ of e₁ = 22 were on disk about 40 seconds before its coset split was printed, and all δ of e₁ = 6 about half an hour before its split. What the files do not show:

- **Step 2.** The only `selfcheck` output in the files ends `sixlines_delta41.log`, after the output of `validate`, and that file was last written after the last `delta` run. So the run of step 2 rests on this account.
- **That the docstring was not edited later.** Its block "THE DEFINITION" was edited once, after the comparison, to change a cross-reference, and in no other way. The list of commands under it grew: `selfcheck` and the choice of twist for the control over F_{p²} were added before any curve at p = 41, the option `--surfaces` in step 6, and after the comparison the command `signs` (§7.1), the option `--out` of `cosets` and a remark on `--no-delta`. These later additions sit below "THE DEFINITION", so the docstring's own statement that nothing below that line changed after the comparison except a cross-reference is too narrow; this paragraph is the complete account. The routine `read_u` and the table of readings were not touched after step 2, except for that choice of twist.

No value of δ was dropped or recomputed, and the script has no parameter that was tuned.

---

## 4. The surfaces at p = 41 [E]

`sixlines.py enumerate --p 41`: 1,975,506 configurations in general position, 2,948 classes under PGL₃ × S₆, 1,234 self-associated. 17 minutes 20 seconds of wall time on 14 processes, 2.7 hours of processor time; the first 4.4 minutes are the canonical forms, on one process. `sixlines_cm.py --list` (11 seconds): 469 classes with ordinary irreducible P in 21 polynomials, and the number of surfaces equals the number of lattices in all 21 (6 with a maximal order, 15 without).

| e | configurations with P, with P(−t) | self-associated | surfaces | lattices per type (note 2) | class number of E, of E₀ | index of R_min |
|---|---|---|---|---|---|---|
| (6, −369, −20172) | 18, 14 | 0 | 16 | 16 | 64, 2 | 1 |
| (22, −369, −73964) | 15, 9 | 0 | 12 | 12 | 48, 2 | 1 |

The curves, command `prym`: all 30 ordered pairs of each of the 56 configurations, 336 seconds on 14 processes. 1,507 computed and 173 degenerate; the second exterior power of the quartic of the curve is P in 763 and P(−t) in 744, and nothing else occurs. The points a, b, c of note 14 are rational in 245 pairs.

`validate`: the checks of note 23 §3 on the first pair of four configurations, on C₀ and C₁ over F_p and C₀ over F_{p²}, two random classes each: 12 curves, 24 classes. The order kills the class, J⁴ = 1, Frobenius^r = 1 and commutes with J, χ(Frobenius) kills the Prym part with i acting as J³ on 24 of 24 and with i acting as J on 0 of 24. So the identification of note 23 holds at p = 41.

---

## 5. The cosets [E]

Command `cosets q e1 1 1 --no-delta`. It runs `sixlines_delta_labels.run_one` (note 24) and `sixlines_review2.cmd_galois` (note 25) without changing them; their fixed input files are replaced in memory by the p = 41 files, and the lattices and period points are computed once for both.

| q | e₁ | digits | period points: wall time, processor time | chart, coordinate | roots of H mod 𝔭 that are enumerated surfaces | primes ℓ ≡ 1 mod 4 | classes of G seen | group generated, by type | orbit indicator recognised with | values of Γ/H′ |
|---|---|---|---|---|---|---|---|---|---|---|
| 41 | 22 | 2,000 | 23 minutes on 4 processes, 88 minutes | J₂, (0, 1, 0, 0, 0) | 12 of 12 | 389, 541, 601, 977, 2,017, 2,089, 2,557, 2,897 | 6 of 12 | 6, 6, 6 | divisor i, signs (−, +), and no other | 9, 32 |
| 41 | 6 | 3,000 | 50 minutes on 8 processes, 6.4 hours | J₃, (2, −2, 7, −7, −8, 3, 5) | 16 of 16 | 1,069, 1,753, 2,137, 2,377, 2,557, 2,801, 4,013, 4,937 | 8 of 16 | 8, 8, 8 | divisor i, signs (+, −), and no other | 9, 32 |

- **e₁ = 22.** At 1,200 digits the class polynomial is recognised but the polynomials of the invariants are not; at 2,000 both are.
- **e₁ = 6.** Run at 3,000 digits only; whether fewer suffice was not tried. One surface has J₂ ≡ 0 modulo 41, so the chart J₂ is not integral at 𝔭; in the chart J₃ the eight coordinates of note 20 give 12 to 15 distinct roots, and the seeded search of note 24 is used, with its acceptance rule (16 simple roots modulo 𝔭). The first coordinate that passes is the one in the table.
- 9 and 32 are the two square roots of −1 modulo 41, as the check of note 24 §3 requires.

The other two routes, in the same command:

| q | e₁ | H over E | H over E(i) | split by the factors equals the orbits | form modulo 4: lattices per type with s = +1, −1 | form indicator recognised with | split by the form equals the orbits |
|---|---|---|---|---|---|---|---|
| 41 | 22 | 12 | 6, 6 | yes | 6, 6 in each type | divisor i, signs (+, +), and no other | yes |
| 41 | 6 | 16 | 8, 8 | yes | 8, 8 in each type | divisor i, signs (+, +), and no other | yes |

And the classes of primes, as in note 25 §4.3:

| q | e₁ | G | primes of degree 1 with ℓ ≡ 1 mod 4, classes | with ℓ ≡ 3 mod 4, classes | the two sets of classes | degrees of the factors of H mod 𝔩 equal the order of the class |
|---|---|---|---|---|---|---|
| 41 | 22 | Z/12 | 40, 6 | 58, 6 | disjoint; the first a subgroup | 98 of 98 |
| 41 | 6 | Z/16 | 40, 8 | 70, 8 | disjoint; the first a subgroup | 110 of 110 |

"Primes" in this table and in §8 are prime ideals of E of degree 1: the 128 with ℓ ≡ 3 mod 4 lie over 46 rational primes (22 and 24 in the two fields), and the 80 with ℓ ≡ 1 mod 4 over 38 [R].

So in both fields the class polynomial is irreducible over E, the group that acts is G, the subgroup found is ker χ₋₁, and Frobenius at a prime of degree 1 with ℓ ≡ 3 mod 4 exchanges the two cosets (the argument of note 25, finding 2).

**Checked independently [R].** The tables of this section reproduce from the stored records. For e₁ = 6, H was not recomputed (it needs 3,000 digits).

- **e₁ = 22, rebuilt.** H was rebuilt from the period points at 1,200 and 2,000 digits, reusing the code of note 20 for the lattices and period points only, and recognised with independent code. Its 12 roots modulo 𝔭 are the stored ones, the invariants at each root are those of exactly two associated configurations, and its two factors over E(i) reduce to the stored split at both primes above 𝔭. The cubic class field of E gives the three cosets of the kernel of a character of order 3 ([note 29](29-characters-of-order-3.md)), and each coset of ker χ₋₁ meets each of them in 2 surfaces. The quartic class field gives four blocks of 3, and each coset of ker χ₋₁ is a union of two of them, which a random split does with probability 3/462.
- **G and χ₋₁, rebuilt.** For e₁ = 22, Cl(E) = Z/24 × Z/2, and the kernel of Cl(E) → G is generated by the ramified prime above 2 and the image of Cl(E₀), so G = Z/12; for e₁ = 6, G = Z/16. E(i)/E is unramified (its relative discriminant has norm 1), and χ₋₁(𝔩) = (−1)^{(N𝔩−1)/2} holds on about 1,200 primes for each field.
- **GRH.** For e₁ = 22 the split itself does not depend on GRH: H modulo 𝔩 is irreducible of degree 12 at 146 primes, so the Galois group of its splitting field over E is cyclic and has one subgroup of index 2. The structure of Cl(E) and G still rests on GRH.
- **Complex conjugation.** For e₁ = 22, H has its coefficients in E₀, not only in E: all 13 are fixed by F ↦ q²/F. It is irreducible over E and over E₀(i), and it already splits into the same two factors over the totally real field E₀(√(4q² − s²)), s = F + q²/F. So complex conjugation preserves each coset and is not what makes the split. Frobenius at 𝔭 is trivial, since 𝔭 is principal, and association keeps a surface on its root.

One thing to know when reading the stored record and log: the code of notes 24 and 25 prints "equals the split by delta: False" and "delta is constant on each orbit …: False". That is because δ was withheld from it here (`--no-delta`). The comparison is made by `compare` (§7).

---

## 6. δ [E]

Command `delta --reading primary`, over F_{41⁴} with N = 32, on C₀. First the first 6 stored pairs of each configuration, then all stored pairs of the five surfaces that had fewer than 4 computed pairs. All counts of this section were recounted independently from the records and agree [R].

| e₁ | surfaces | configurations | curves tried | computed | exponent 64 | no kernel vector gives u | not verified | degenerate | pairs with a δ, by surface | surfaces with δ = +1, −1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | 16 | 32 | 352 | 164 | 182 | 4 | 2 | 0 | 5 to 23 | 8, 8 |
| 22 | 12 | 24 | 180 | 86 | 87 | 4 | 1 | 2 | 4 to 21 | 6, 6 |
| both | 28 | 56 | 532 | 250 | 269 | 8 | 3 | 2 | | 14, 14 |

- **Every curve computed passed the four direct checks** of the reading rule. In the 3 marked "not verified" the point is not in the kernel of π − π₁ (the twelve classes did not span, as in the one case of note 23); they are not counted, and their values would have agreed with their surfaces.
- **The point.** X has order 16 on all 250 curves. M is 12 on 119, 4 on 121, and twice an odd number on 10. So u = 3 on 125 curves and u = 1 on 125.
- **Inside a curve.** On each of the 250 curves every kernel vector that gives u gives the same u, and no kernel vector is without a solution M. Vectors with M ≡ 8 or 0 modulo 16 give nothing, as the rule says. The 8 curves with no usable vector have only such vectors.
- **The complement (§2.4).** On 240 of the 250 curves the kernel tallies fit M = N(γ)h(m) exactly: order 16 with v₂(M) = 2, order 8 with v₂(M) ≥ 3 [R]. On the 10 curves whose chosen M is 2·(odd), there are order-16 vectors with M = 2·(odd) and with M = 4·(odd) together, and 9 order-8 vectors with M ≡ 0 mod 16, which the sketch does not predict [R]. Throughout, over each kernel the values of M that are not 0 modulo 8 have one odd part modulo 4, as at q = 29 and 37, so u and δ are unchanged. What the data confirm is that statement, not the sketch.
- **The convention (§2.5).** Over the 211 computed records for e₁ = 6 with N = 32 (primary and controls), 1,855 kernel vectors have M not 0 modulo 16; 40 of them, on 7 records, have M ≡ 2 mod 4, and exactly those 40 have no solution with the other convention. No kernel vector in the 783 records has odd M, so these are all the vectors on which the conventions differ, and "no solution with the other convention ⇔ M ≡ 2 mod 4" holds on 211 of 211 records [R]. For e₁ = 22 the two conventions have the same solutions on all 1,073 such vectors of 122 records, as §2.5 says they must. So the convention of note 23 is confirmed at q = 41 on 40 vectors and cannot be tested on the rest.
- **u depends on the pair, δ does not.** u is the same on all computed pairs of 0 of the 28 surfaces. u·sign(Re c₁·Im c₁) is constant on 3, u·sign(c₂) on 0, and the product of the two signs without u on 0. δ is constant on 28 of 28, over 250 pairs.
- **Balance.** After the first run the first polynomial stood at 8 and 7 with one surface without a pair. The second run gave it δ = −1 on 21 of 21 pairs.
- **Time.** The stored times are 41 to 106 seconds per curve in the first run (12 processes) and 66 to 128 seconds in the second (14 processes) [R].

Controls (note 23 §5.5), one computed curve of each surface unless said otherwise, each compared with the primary u of the same curve:

| control | runs | computed | same u | |
|---|---|---|---|---|
| other pairs of lines, of the configuration and of its associate | 250 | 250 | δ the same on 28 surfaces of 28 | the primary run itself |
| the twist C₁, F_{41⁴}, N = 32 | 28 | 27 | 27 of 27 | one curve with no usable vector |
| the twist C₂, F_{41⁴}, N = 32 | 28 | 28 | 28 of 28 | |
| other random classes (seed 7) | 28 | 28 | 28 of 28 | |
| F_{41²}, N = 16, the twist with s = 6 | 28 | 4 | 4 of 4 | 24 curves have no vector with M = 2·(odd) |
| F_{41⁴}, N = 16 | 28 | 0 | | the classes are killed by 32 and not by 16: not a usable control as defined |
| extension: the classes doubled, on curves of exponent 64 | 111 | 0 | | every vector has M ≡ 0 or 8 mod 16 |

The twist and seed rows check less than their names suggest. Over F_{41⁴} every element of F₄₁^× is a fourth power, so the twists C₁ and C₂ are isomorphic to C₀ there, and these controls check only the bookkeeping of twists and random classes on the same group [D]. The evidence that δ does not depend on the curve is the 250 pairs of the primary run.

The last two rows are my own design and should not have been planned: a class killed by 32 cannot be paired with N = 16 without first cutting down to the 16-torsion, and a point 2X′ has M = 4·M(X′) with M(X′) even on every curve seen, which leaves no digit to read. They cost 28 and 111 runs and say nothing about δ either way. The curves of exponent 64 stay out of reach of this method; they would need the pairing of order 64, over F_{41⁸}.

---

## 7. The comparison [E]

Command `compare`. A configuration is (b, c, e, f) of note 2; a surface is a configuration and its associate, on both sides. "Orbit" is the orbit of the ideals above primes ℓ ≡ 1 mod 4, named A for the smaller value of Γ/H′ of its indicator. "Factor" is the factor of H over E(i) of which the surface is a root at one of the two primes above 𝔭. "Form" is the value of Γ/H′ for the indicator of the sign s read on the lattices; equal values mean equal s.

**q = 41, e = (6, −369, −20172), G = Z/16**

| orbit | root of the class polynomial mod 𝔭 | configurations (b, c, e, f) | factor | form | δ | pairs with a δ |
|---|---|---|---|---|---|---|
| A | 3 | (2, 5, 17, 8), (2, 17, 5, 8) | 2 | 32 | −1 | 7 |
| A | 14 | (2, 5, 27, 19), (2, 7, 30, 34) | 2 | 32 | −1 | 21 |
| A | 25 | (2, 5, 7, 39), (2, 7, 5, 39) | 2 | 32 | −1 | 21 |
| A | 33 | (2, 3, 16, 10), (2, 7, 39, 30) | 2 | 32 | −1 | 6 |
| A | 34 | (3, 7, 10, 32), (3, 8, 20, 9) | 2 | 32 | −1 | 23 |
| A | 38 | (2, 6, 7, 9), (2, 6, 35, 26) | 2 | 32 | −1 | 6 |
| A | 39 | (2, 17, 27, 20), (2, 7, 24, 23) | 2 | 32 | −1 | 8 |
| A | 40 | (2, 3, 19, 23), (2, 3, 8, 19) | 2 | 32 | −1 | 7 |
| B | 2 | (2, 5, 26, 18), (2, 6, 12, 13) | 1 | 9 | +1 | 7 |
| B | 8 | (3, 8, 36, 9), (3, 4, 24, 27) | 1 | 9 | +1 | 6 |
| B | 10 | (2, 7, 36, 12), (2, 7, 8, 6) | 1 | 9 | +1 | 7 |
| B | 12 | (2, 5, 31, 15), (2, 7, 18, 27) | 1 | 9 | +1 | 6 |
| B | 15 | (2, 5, 35, 24), (2, 6, 7, 15) | 1 | 9 | +1 | 5 |
| B | 17 | (2, 5, 17, 16), (2, 6, 8, 11) | 1 | 9 | +1 | 7 |
| B | 21 | (2, 3, 19, 36), (2, 5, 27, 24) | 1 | 9 | +1 | 21 |
| B | 29 | (2, 5, 35, 16), (2, 6, 8, 38) | 1 | 9 | +1 | 6 |

**q = 41, e = (22, −369, −73964), G = Z/12**

| orbit | root of the class polynomial mod 𝔭 | configurations (b, c, e, f) | factor | form | δ | pairs with a δ |
|---|---|---|---|---|---|---|
| A | 0 | (2, 5, 27, 39), (2, 5, 35, 22) | 1 | 32 | −1 | 4 |
| A | 5 | (2, 5, 16, 7), (2, 3, 27, 9) | 1 | 32 | −1 | 6 |
| A | 22 | (2, 5, 17, 15), (2, 5, 26, 16) | 1 | 32 | −1 | 21 |
| A | 32 | (2, 6, 24, 23), (2, 6, 31, 38) | 1 | 32 | −1 | 6 |
| A | 33 | (2, 5, 19, 14), (2, 17, 13, 25) | 1 | 32 | −1 | 5 |
| A | 38 | (2, 5, 16, 17), (2, 6, 31, 24) | 1 | 32 | −1 | 7 |
| B | 1 | (2, 7, 36, 33), (2, 7, 8, 23) | 2 | 9 | +1 | 5 |
| B | 11 | (2, 5, 16, 30), (2, 6, 27, 11) | 2 | 9 | +1 | 5 |
| B | 17 | (2, 5, 24, 16), (2, 17, 11, 20) | 2 | 9 | +1 | 8 |
| B | 25 | (2, 6, 24, 4), (2, 6, 31, 23) | 2 | 9 | +1 | 5 |
| B | 28 | (2, 17, 27, 15), (2, 7, 24, 11) | 2 | 9 | +1 | 8 |
| B | 39 | (2, 6, 27, 29), (2, 7, 23, 25) | 2 | 9 | +1 | 6 |

In both tables the four columns orbit, factor, form and δ are the same split, and δ is constant on A, constant on B, and the two differ. The comparison is between partitions: the name of an orbit, the number of a factor and the sign of δ for a polynomial are conventions.

### 7.1 The sign between δ and the form, on all five torsors [E]

This was not planned and was looked for after the comparison above. Note 25 §4.4 says δ is the form of the lattice "up to one sign for each polynomial", and note 23 calls the overall sign of δ for a polynomial a convention. In the two tables the form value 32 goes with δ = −1 in both polynomials, and 32 is ζ = n₀^{(q−1)/4}, the fourth root of unity by which `sixlines_jacobian.py` lets J act (n₀ = 6 the least primitive root of 41). The form value at a root r is Γ(r)/H′(r) = s/i, with s the sign of the lattice, so the rule would be

> δ(X) = ζ·Γ(r)/H′(r) modulo q: δ is the sign s of the lattice of the CM lift, read at the prime of E(i) above 𝔭 where i goes to ζ.

To see whether it holds beyond q = 41 the command `cosets` was run on the three torsors of note 24 (with δ of note 23), which stores the form value of each root; `signs` then tests the rule:

| q | e₁ | ζ | surfaces with (form value, δ) | δ = ζ·(form value) on every surface |
|---|---|---|---|---|
| 29 | 2 | 12 | (12, −1): 6; (17, +1): 6 | yes |
| 37 | 2 | 31 | (31, −1): 6; (6, +1): 6 | yes |
| 37 | 18 | 31 | (31, −1): 4; (6, +1): 4 | yes |
| 41 | 6 | 32 | (32, −1): 8; (9, +1): 8 | yes |
| 41 | 22 | 32 | (32, −1): 6; (9, +1): 6 | yes |

Five torsors of five, 60 surfaces of 60. If the sign of each polynomial were a coin toss, one rule with one free sign would fit five torsors once in 16, but that is the figure before looking. I formed the rule after seeing four of the five rows; the fifth (q = 37, e₁ = 2) came afterwards and is the only test of it: one chance in 2. So this is weak evidence, and of a different standing from §1 item 1: it says that the sign of δ for a polynomial is possibly not a convention but fixed by the lattice once ζ is fixed. The table reproduces from the stored files, with δ at q = 29 and 37 recomputed from the records of note 23 [R]. The values of the orbit indicator follow no single rule: δ = +ζ·(value) at (29, 2) and (37, 2), but −ζ·(value) at (37, 18) and at q = 41 [R]. That indicator's sign is not canonical, which is why the form indicator is the one used here. Nothing in the test of §1 depends on it. It rests on conventions I did not trace: which of J and J³ is i, the orientation of note 20, and the order of the types.

The rerun also gives, for the three earlier torsors, the same orbits, factors and form split as notes 24 and 25 (12, 12 and 8 surfaces).

---

## 8. Check and limits

- **One separate check.** The check was made with its own code from the stored files, by the same kind of method as this note; it is not an outside review. It recomputed δ, u, Z and the reading rule from the stored M and c of all 783 records (250 primary readings, 3 not verified, 87 controls), recounted every count of §1, §6 and §7, rebuilt H for e₁ = 22 with its split, and rebuilt G and χ₋₁ for both fields. It tried to explain δ by something else: c₁, the two signs, the compare flag, the configuration's S₁, rational points and the pair itself, and all 2⁷ products of sign features with and without u. δ is a function of none of them, and the only non-trivial product constant on all 28 surfaces is the defined one [R].
- **What the check did not cover.** The pairing values themselves, which are stored as booleans only (their internal consistency, Z against M, parity and kernel tallies, was checked); c₁ and c₂ from the curves; H for e₁ = 6; the orbit and form routes, whose permutations are not stored; the identification of the Galois action on the roots with the action of G (note 20 §7.3); the enumeration at p = 41; a certification of the class groups; the trial of step 4 of §3; and whether the docstring was edited, which nothing records. The scripts of notes 14, 20, 23, 24 and 25 are used unchanged and carry their own limits; note 25 read the Python side of three of them and not the GP code of the group law, of Miller's algorithm or of the enumeration of lattices.
- **The new code** is `read_u`, a copy of `character_pair` of note 23 with the reading rule of §2.3 in place of that of note 23 §5.3, which reads M from Z = M·conj(N/ϖ₁) mod N. It reproduces the stored M, u, δ and order on 6 stored curves of note 23. The agreement of §1 would be hard to obtain from a wrong reading: u is not constant on any surface, and the two signs that turn it into δ are functions of the quartic alone.
- **The derivations of §2 are mine.** The check reproduced their numbers and corrected two imprecisions (§2.2, §2.5); it did not re-derive the identity for Z or the argument for the complement, which are those of note 23, which calls them sketches. Their support at q = 41 is the direct verification of each point and the consistency over each kernel (§6).
- **"Fixed first" rests on my account** of §3 for the part the files do not show (step 2, and the docstring). The file times support "δ before cosets" for both polynomials [R].
- **One prime.** q = 41 ≡ 9 mod 16, π₁⁴ − 1 of valuation exactly 8. For q ≡ 1 mod 16, or when β of §2.2 is not a unit, the valuation is larger and M is known to more digits; that case was not met. q = 73, 89, 97 would be the next primes ≡ 1 mod 8; none is enumerated.
- **Half of the curves.** 269 of 532 have exponent 64 and are not computed. If δ on those disagreed with the rest, δ would depend on the pair. Every surface has at least 4 computed pairs, and computed pairs from both of its configurations (28 of 28), so which surfaces got a δ is not biased. But exponent 64 is not spread evenly: in the first run it occurs on 96 of the 168 curves of coset A (δ = −1) and on 77 of the 166 of coset B, Fisher p = 0.06 with both polynomials together and 0.38 and 0.36 for each polynomial over all runs; 4 of the 5 surfaces of the second run are in A [R]. The data cannot exclude that the curves of exponent 64 would carry the other sign.
- **The conjugation of c₁.** On 169 of the 250 computed curves π₁ ≡ π̄₁ modulo 16, so the direct check πX = π₁X cannot tell c₁ from c̄₁, and the check "P or P(−t)" holds for every c₁ = 4(±1 ± i) [R]. On those curves sign(Re c₁·Im c₁) rests on the point counts of note 14 and on `validate` (12 curves). An error on single curves would show as a δ that varies within a surface, which does not occur; an error on all curves would only flip the overall sign, which is a convention.
- **No step outside the kernel** was applied to a lattice or a surface. That δ changes sign under such a step is the inference of note 24 §2, supported here as in note 25 by the irreducibility of H over E and the cycle lengths of Frobenius at 128 prime ideals of degree 1 with ℓ ≡ 3 mod 4, over 46 rational primes.
- **"The action" is that of note 20**, carried to F_q by reduction of CM points; that it is the action of the class group on the K3 surfaces themselves is the sketch of note 20 §7.7.
- **Recognition** of H, of the invariants and of the indicators is by the rule of note 20. For e₁ = 6 the coordinate comes from the seeded search, whose acceptance rule does not see δ.
- **Class groups** are PARI's without certification (GRH); the factorisations over E and E(i) are PARI's `nffactor`. For e₁ = 22 the split itself does not depend on GRH (§5).
- **The enumeration at p = 41** was run once. Its support is that the two polynomials have the predicted numbers of surfaces, that all 21 polynomials match their lattice counts, and that all 1,507 curves give P or P(−t).

---

## 9. What this does not give

- A proof that δ(g·X) = χ₋₁(g)·δ(X), at q = 41 or in general.
- The statement for the other 75 polynomials of note 19 with the character, for any prime above 41, or for q ≡ 1 mod 16.
- δ on the pairs of lines with exponent 64, or any pairing of order above 32.
- A reading of u over F_p or, for most curves, over F_{p²}: the computation still needs F_{p⁴}.
- A step at a prime ℓ ≡ 3 mod 4, on lattices or on surfaces; full labels in G for the 28 surfaces.
- Why δ equals the lattice form. It is observed on 60 surfaces; the Prym variety has a polarisation of type (1, 1, 1, 2) and the lattice a principal one (note 25 §7).
- Anything about the other characters of order 2, the polynomials without χ₋₁, non-maximal orders, or the search problem.
- A costed algorithm for δ at cryptographic size. The cost paragraph of note 24 §5 stands as it is, with one more data point: 40 to 130 seconds for one curve at p = 41.

---

## Reproduce

From `experiments/`, with `PY=python`:

```bash
$PY -I sixlines.py enumerate --p 41 --out sixlines_p41.jsonl --procs 14
$PY -I sixlines_cm.py --list --out sixlines_cm_p41.jsonl sixlines_p41.jsonl > sixlines_cm_p41.log
$PY -I sixlines_delta41.py dyadic
$PY -I sixlines_delta41.py selfcheck --n 3 --workers 6
$PY -I sixlines_delta41.py list
$PY -I sixlines_delta41.py prym --pairs 0 --workers 14
$PY -I sixlines_delta41.py validate --n 2
$PY -I sixlines_delta41.py delta --reading primary --pairs 6 --workers 12
$PY -I sixlines_delta41.py cosets 41 22 1 1 --digits 1200 2000 3000 --primes 8 --workers 4 --no-delta
$PY -I sixlines_delta41.py delta --reading primary --pairs 0 --surfaces 6:7 6:15 6:2 6:3 22:3 --workers 14
for R in twist1 twist2 seed n16 r2; do $PY -I sixlines_delta41.py delta --reading $R --per-surface 1 --workers 14; done
$PY -I sixlines_delta41.py delta --reading cofactor --pairs 2 --workers 14
$PY -I sixlines_delta41.py cosets 41 6 1 1 --digits 3000 4200 --primes 8 --workers 8 --no-delta
$PY -I sixlines_delta41.py compare
for a in "37 18 1 1 --digits 1200" "29 2 1 1 --digits 1200" "37 2 1 1 --digits 2000"; do $PY -I sixlines_delta41.py cosets $a --primes 3 --frob-primes 2 --workers 14 --out sixlines_delta41_signs.jsonl; done
$PY -I sixlines_delta41.py signs
```

Times on 16 cores: 17 minutes; 11 seconds; a second; 2 minutes; a second; 6 minutes; 2 minutes; 30 minutes; 30 minutes (6 at 1,200 digits, 24 at 2,000, while the run before it was using 12 processes); 18 minutes; 13 minutes for the five controls; 11 minutes; 52 minutes; a second; 9 minutes for the three torsors of §7.1; a second. Several were run at the same time, so these are upper bounds. `delta` appends to `sixlines_delta41.jsonl` and skips what is already there, so the order of the two primary runs matters for the choice of pairs and not for the result.

Files written, all new: `sixlines_p41.jsonl`, `sixlines_enumerate_41.log`, `sixlines_cm_p41.jsonl`, `sixlines_cm_p41.log`, `sixlines_prym_p41.jsonl`, `sixlines_delta41.jsonl` (783 records), `sixlines_delta41_cosets.jsonl` (two records for each polynomial), `sixlines_delta41_compare.json`, `sixlines_delta41_signs.jsonl` (§7.1; the three reruns at q = 29 and 37, which read δ from the files of note 23), and the logs `sixlines_delta41.log` (`list`, `validate`, `selfcheck`), `sixlines_delta41_primary.log`, `sixlines_delta41_more.log`, `sixlines_delta41_cosets_22.log`, `sixlines_delta41_cosets_6.log`, `sixlines_delta41_signs.log`, `sixlines_delta41_compare.log` (`compare`, then `signs`). No file of an earlier note was changed, and none of their scripts. The script imports `sixlines_prym.py`, `sixlines_jacobian.py`, `sixlines_cmmethod.py`, `sixlines_delta_labels.py` and `sixlines_review2.py`; `selfcheck` reads `sixlines_jacobian_character.jsonl`, and the reruns of §7.1 read it and `sixlines_jacobian_character_extra.jsonl`, `sixlines_cm.jsonl`, `sixlines_p29.jsonl` and `sixlines_p37.jsonl`. The random classes are seeded. No package was installed.

## Sources

No new source was read. In this repository: notes [2](02-six-line-count.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [29](29-characters-of-order-3.md). The paper of Castryck, Sotáková and Vercauteren (CRYPTO 2020, https://eprint.iacr.org/2020/151) is cited as read in notes 23 and 25.
