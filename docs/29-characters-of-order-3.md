# 29. Characters of order 3 of G against invariants of six-line surfaces

*[Note 28](28-leakage-audit.md) could not see characters of order 3 of G: among its labelled torsors only two have 3 ∣ \|G\|, both with q ≡ 2 mod 3, where no residue symbol of order 3 exists, and its statistics could not reach significance on them in any case. This note tests them. It computes the three cosets of the kernel of a character χ₃ of order 3 on seven torsors, by factoring the class polynomial of [note 20](20-cm-method-and-the-action.md) over the cubic class field, and tests residue symbols of order 3 and the families of note 28 against those cosets with exact permutation tests fixed before the join. Notes [19](19-quadratic-characters-from-a-surface.md), 20, [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [26](26-delta-is-the-character-at-q-41.md) and 28 supply the data and the code. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_order3.py`; **[H]** heuristic; **[R]** re-derived in a separate re-derivation check, with its own code (not part of this repository). Statements of other notes keep the label they have there. Checked once by a separate re-derivation from the stored files, with its own code: the cubic class fields of all seven torsors, the class polynomial over L for the two Z/6 torsors, 8,880 cubic feature values from the raw configurations, and every statistic, plant and shuffle of §5 to §7. What it confirmed is marked [R], its corrections are made in the text, and what it did not cover is listed in §9.*

---

## 1. Result

**No invariant in the families tested reads a character of order 3, or separates the cosets of its kernel, beyond chance [E].** Characters of order 3 can be tested as characters only on the two torsors with q ≡ 1 mod 3 (q = 37, G = Z/12 twice, 24 surfaces). There are 184 residue symbols of order 3 on them. None is a character of order 3 on either torsor, and none is a function of the coset difference beyond chance (best p 0.030, Holm 1); at q = 37 the cube roots of unity are themselves cubes, which limits what these symbols can see (§2). Among all 2,839 candidates of that population the smallest p is 0.0046, Holm 1; it is the exact indicator of one coset on one torsor, at chance level (§5.1). On the five torsors with q ≡ 2 mod 3 (48 surfaces) there is no symbol of order 3 (§2), so the tests there ask only whether some invariant separates the cosets. None does: the smallest p of 2,448 is 1.4·10⁻⁴, Holm 0.34, and it does not replicate on the held-out torsors at q = 37. The per-pair fourfold quantities do not separate them either: the best is u, p 0.0061, Holm 0.080. The positive control comes out as predicted. The populations are small, and the power has the limits that item 7 states. In detail:

1. **The cosets [E, R].** There are seven torsors with 3 ∣ \|G\| and 72 surfaces: Z/12 twice at q = 37; Z/6 at q = 23 and at q = 29; Z/12 at q = 29; and Z/12 twice at q = 41. On each:
   - the class polynomial H modulo 𝔭 has exactly the enumerated surfaces as its roots;
   - over the cubic class field L of ker χ₃, H splits into three factors of degree \|G\|/3, whose roots at one prime of L above 𝔭 partition the surfaces;
   - the partition is the same at every prime of L above 𝔭;
   - at 60 or 62 primes of E of degree 1 per torsor, 424 in all, the relative cubic splits exactly when the class of the prime lies in ker χ₃ (424 of 424).

   Two of the seven class polynomials are new: q = 37, e₁ = 6, which note 20 skipped, and q = 41, e₁ = 2. Two checks tie the cosets to other computations, and neither is a test of an invariant (§7.1). On the two Z/6 torsors that have full labels from note 20, the cosets are the labels modulo 3. On the three torsors that also have the cosets of ker χ₋₁ from notes 24 and 26, each coset of ker χ₃ meets each coset of ker χ₋₁ in exactly 2 surfaces, as Z/12 ≅ Z/4 × Z/3 requires. The separate check confirmed the stored cubic fields on all seven torsors at 4,569 primes and refactored H over L for the two Z/6 torsors [R].
2. **q = 43 did not finish in this run [E].** None of its three torsors had a class polynomial when the tests were run: Z/24 at e₁ = 2, Z/12 at e₁ = 14, Z/6 × Z/2 at e₁ = 22. So in this run the population with q ≡ 1 mod 3 is q = 37 alone. Where each run stopped, and why, is in §4.2. **After the run (§7.4)**, reruns gave cosets on all three, and A₄₃ was scored. No cubic symbol is a character of order 3 there either: 0 of 552 cells, best Holm 0.28, at a prime where the cube roots of unity are not cubes. On A, now five torsors, the smallest p is 3.3·10⁻⁵, Holm 0.095. A₄₃ has one candidate past Holm, the quadratic symbol of R2(X) + R2(Y) (Holm 0.0057). It reverses direction at q = 37 and does not replicate on the held-out Z/24 torsor (p 0.22).
3. **Order 3 at q = 37 [E, R].** The 184 symbols are:
   - 40 J-monomials;
   - 121 monomials in the power sums and the conic;
   - 22 fixed quantities of the lines;
   - the cubic symbol of the CM coordinate.

   Each is defined on both torsors and agrees on the two configurations of every surface. Of the 368 (candidate, torsor) cells, 0 are exact characters of order 3. 259 of the 368 contain a surface where the quantity vanishes, and such a cell cannot be a character, so only 109 could have been one [R]. The diff test has best p 0.030 (E3, p₁p₂p₃p₄²C) and the pair test 0.012 (J3, J₄J₅J₆²), both Holm 1. The cubic type of the ten products has p 0.46.
4. **The families of note 28 [E].** On population A, the smallest p is 0.0046, for the quadratic symbol of J₃J₅ as a pair of values, Holm 1, Benjamini–Hochberg (BH) 1. It is the exact indicator of one coset on q = 37, e₁ = 2 (§5.1). 57 of the 2,832 defined candidates have p ≤ 0.05. On population B the smallest is 1.4·10⁻⁴, Holm 0.34, BH 0.17, for the quartic symbol of p₁²p₂p₃²p₄C². All its weight is on one torsor (q = 41, e₁ = 2). On the two torsors at q = 37, which were not used to find it, it has p = 1 (§7.2). 67 of 2,443 have p ≤ 0.05.
5. **The fourfold side per pair of lines [E].** On three torsors with 354 records: u has p 0.0061 (Holm 0.080), the status of the computation 0.015, M mod 16 0.030, and the rest p ≥ 0.32. δ has p 0.999 against the cosets of order 3. That is expected but not forced: δ is the character χ₋₁ (notes 24 and 26), and on Z/12 the cosets of ker χ₋₁ and of ker χ₃ are independent (item 1), but the statistic is over records, 6 to 48 per surface, not over surfaces.
6. **The positive control [E].** δ against the cosets of ker χ₋₁ on the three torsors that have them, scored with the same code: p = 1.01·10⁻⁸ = (1/462)³, the number fixed beforehand. It is three orders of magnitude below the thresholds of A (1.8·10⁻⁵) and B (2.0·10⁻⁵).
7. **Power [E, R].** Planted signals were scored against the same nulls. "Found" means p ≤ 0.05/m (§6.2). Counts are out of 200 plants each. A plant "replaced with probability r" replaces each value by a uniformly random element of Z/3 with probability r, so the value is wrong with probability 2r/3.

   | plant | A (2 torsors) | B (5 torsors) |
   |---|---|---|
   | exact order-3 character on every torsor | 200 | 200 |
   | each value replaced with probability 1/8 (wrong with probability 1/12) | 173 | 196 |
   | each value replaced with probability 1/4 (wrong with probability 1/6) | 76 | 161 |
   | exact on 1 torsor | 4 | 10 |
   | indicator of one coset | 0 | 200 |
   | a random function of c(Y) − c(X) | 200 | 200 |
   | a random function (null) | 0 | 0 |

   With values actually wrong with probability 1/8, replays find 127 and 131 of 200 on A and 172 and 178 on B; with probability 1/4, 30 and 31 on A and 74 on B [R]. On A the indicator of a coset cannot be found at all: its best possible p on two torsors of 12 is (35/5,775)² = 3.7·10⁻⁵, above the threshold.
8. **Shuffled labels [E].** Population A was rerun with the cosets permuted at random, seeds 1 to 10. The smallest p ranged from 3.4·10⁻⁵ to 0.0039, against 0.0046 for the true cosets, which are less extreme than all ten.

**What this means [D, H].** A CSIDH-style attack on the decision problem reads a character from each object (note 24 found one reading, δ for χ₋₁). For order 3 such a reading would need a cubic quantity. Within the families of note 28 and the cubic families added here, nothing reads χ₃ on the two torsors where a cubic symbol exists. That holds for quantities whose cubic class moves by a non-cube between cosets: at q = 37 the cube roots of unity are cubes, so a reading of the kind §2 describes would have a constant cubic symbol there. For q ≡ 2 mod 3, every quantity defined over F_q is a cube (§2), so a reading of χ₃ there would need a quantity outside F_q, and none was tried. The test would have seen a cubic symbol that is a character on both torsors at q = 37, even with one value in twelve wrong. It would not have seen a reading on one torsor only, an indicator of one coset at q = 37, or a weak bias. The groups have order 6 and 12, and only two torsors have q ≡ 1 mod 3. So this is not evidence about cost at any size where the problem could be hard. It extends the negative result of note 28 from characters of order 2 to order 3, tested as characters on the two torsors at q = 37 and, on the other five, only through separation tests together with the argument of §2. The leak that certainly exists is unchanged: the label as a function of the complete invariants, through the class polynomial (note 28 §7). Here that leak is the factorisation of H over L.

---

## 2. On paper [D]

**The characters and their kernel.** On every torsor treated the 3-part of G is cyclic (Z/3 in Z/6, Z/12, Z/24, Z/6 × Z/2), so G has two characters of order 3, χ₃ and χ̄₃ = χ₃², with one kernel of index 3. Its class field is a cyclic cubic extension L/E, unramified, inside the field generated over E by the roots of the class polynomial H (note 19 §2.3, and note 25 §4.3 for the quadratic case).

**A partition is a labelling.** Every bijection of a set of three elements onto Z/3 is x ↦ ±x + c, because the affine group of Z/3 has order 6 = 3!. So the partition of the surfaces into the three cosets of ker χ₃ fixes the map c : surfaces → Z/3 up to translation and sign, which is exactly the ambiguity of a label (note 20 §7.5) composed with χ₃ or χ̄₃. Nothing beyond the partition is needed for a test of order 3.

**The cosets from the class polynomial.** Over L, H splits into three factors of degree N/3, permuted by Gal(L/E) ≅ Z/3, whose roots over C are the three cosets. Over a prime field the class of 𝔭 = (q, F − e₁) in G is trivial (note 19 §2.4), so 𝔭 splits completely in L; at one prime 𝔓 of L above 𝔭 the roots of each factor are the surfaces of one coset. This is note 25 §4.3 with E(i) replaced by L, and L is computed by PARI's `bnrclassfield` from the subgroup of Cl(E) that is the preimage of ker χ₃.

**The character test.** For a function f from the surfaces to Z/3: f(Y) − f(X) is a function F of c(Y) − c(X) for all X, Y if and only if f = e·c + constant with e ∈ Z/3. (Take X, Y, Z in the three cosets: F(a + b) = F(a) + F(b), so F is additive.) If f is not constant, e = ±1 and f is χ₃ or χ̄₃ of the label, up to a constant.

**Why the residue field matters.** L is not a Kummer extension of E: ζ₃ is not in E. On the ten polynomials of §3, E has no quadratic subfield at all and x² + x + 1 stays irreducible over it [E, `sixlines_order3.py fields`]. It becomes one over E(ζ₃): L(ζ₃) = E(ζ₃)(∛d). A reading of χ₃ by a cubic residue symbol is a function R of the six lines with R(X)³ = (cube)·d, the analogue of note 19 §2.3, and its symbol lives in the residue field of a prime of E(ζ₃) above 𝔭: F_q when q ≡ 1 mod 3, F_{q²} when q ≡ 2 mod 3.

**At q = 37 such an R has a constant cubic symbol.** 37 ≡ 1 mod 9, so the cube roots of unity are cubes in F₃₇: they are 10 and 26, and 26 = 2¹². R takes values r₀·ζ^{k(X)} with R³ = c³·d, so its values on the three cosets differ by cube roots of unity, and its cubic symbol is the same on all three. Such an R would enter the J3, E3, N3 and X3 tests of §3 as a constant, and they cannot see it at q = 37. At q = 43 (≡ 7 mod 9) they could; on B, ζ₃ is not a cube in F_{q²} for q = 23, 29 and 41. For the eleven absolute invariants of B2, such a reading would still show, as a zero of v(Y) − v(X) exactly when X and Y share a coset. The statements of §1 and §8 about cubic quantities are limited accordingly: to quantities whose cubic class moves by a non-cube between cosets.

**For q ≡ 2 mod 3 there is no cubic symbol to test.** Then 3 ∤ q − 1, so x ↦ x³ is a bijection of F_q* and every element of F_q is a cube. In F_{q²}, for x ∈ F_q*, x^{(q²−1)/3} = (x^{q−1})^{(q+1)/3} = 1, so every element of F_q is a cube in F_{q²} as well. All the invariants of note 28 and of §3 are defined over F_q (the lines are rational), so their cubic symbols are trivial over F_q and over F_{q²}. A cubic symbol over F_{q²} can only be non-trivial on a quantity that is not in F_q; for x ∈ F_{q²} it is the cubic character of x̄/x on the norm-one torus. Square roots of F_q-quantities give nothing either: for x² ∈ F_q, x̄/x = ±1, and −1 is a cube in F_{q²} because (q + 1)/3 is even when q ≡ 5 mod 6, so the cubic symbol of x is trivial [D, R]. A non-trivial symbol needs x̄/x of order divisible by 3 on the norm-one torus. I know of no natural quantity of six rational lines that has this property; the candidates need an auxiliary choice (a root of a quadratic formed from the invariants, a point of an auxiliary curve). So on B only the non-homomorphic tests of §3 are run: whether any invariant, through its quadratic or quartic symbol or its value, separates the cosets of order 3.

**What CSIDH suggests [M].** In the analogy of note 19 §2.2, G corresponds to Cl(O)², and squaring is a bijection on 3-parts, so a character of order 3 of G corresponds to one of Cl(O). For imaginary quadratic orders the characters read from curves by pairings are genus characters, of order 2; I do not know a method that reads a character of order 3 from two curves. The test here is empirical and says nothing about that question for curves.

---

## 3. The tests, as fixed before the join [D]

These tests were fixed before any invariant was tabulated against a coset of order 3 (the order of events, and what "before" can be checked against, are in §4.3). What had been seen then: the cosets of ker χ₃ on q = 23, e₁ = 2; q = 29, e₁ = 6; q = 29, e₁ = 2 (roots of H and their partition, no invariant other than the CM coordinate), and the agreement of the first two with note 20's labels modulo 3. The text below follows the tests as then fixed, with the three amendments made before the join.

**Labels.** On each torsor, c(X) ∈ Z/3 is the coset of ker χ₃ containing the surface X, χ₃ a character of order 3 of G (all target groups have cyclic 3-part, so there is one kernel). The partition determines c up to translation and sign, and every statistic below is invariant under both.

**Populations.** A: q ≡ 1 mod 3 (q = 37, e₁ = 2, 6; q = 43, e₁ = 2, 14, 22, those that converge). B: q ≡ 2 mod 3 (q = 23, e₁ = 2; q = 29, e₁ = 2, 6; q = 41, e₁ = 2, 22, those that converge). F: the per-pair-of-lines fourfold records, on the torsors of A ∪ B that have them (q = 37, e₁ = 2 from note 23 and q = 41, e₁ = 22 from note 26). Holm is applied over all candidates of each population separately.

**Surface value.** As in note 28 §2.2: for a function v of a configuration, the surface value is the sorted pair over the configuration and its associate; the *diff* variant is used only on torsors where v agrees on the two configurations of every surface.

**Candidates.**

- *A only, residue symbols of order 3* (discrete logarithm modulo 3 to the least primitive root, 'z' for 0):
  - **J3**: monomials J₂^{a₂}⋯J₆^{a₆}, a ∈ (Z/3)⁵ \ {0}, weight Σ k·a_k ≡ 0 mod 3, one of each conjugate pair {a, 2a} (40 candidates);
  - **E3**: monomials p₁^{a₁}⋯p₅^{a₅}C^b, weights 2, 4, 6, 8, 10, 2 in units of t, the same rule (121 candidates). −1 is a cube, so the sign of C, which changes under odd relabellings of the six lines (note 28 §2.2), does not matter;
  - **N3**, fixed natural quantities of the lines: σ₃, σ₆, σ₉, the elementary symmetric functions of the ten t_s²; the products σ_kσ_l, 1 ≤ k ≤ l ≤ 10, k + l ≡ 0 mod 3 (18); V = ∏_{s<s′}(t_s² − t_{s′}²), weight 90 (its sign does not matter). 22 candidates;
  - **X3**: the cubic symbol of the CM coordinate (the root of H mod 𝔭), and the *cubic type*: the vector of the cubic symbols of the ten t_s, in normal form under S₆, a common shift and the sign (2 candidates).
  - Every monomial and quantity is checked for invariance under random A ↦ gAD on three configurations of each field, as in note 28, and dropped if it fails.
- *A and B, the families of note 28* with the orders available at q: S (note 19's bits and type, the quadratic and, for q ≡ 1 mod 4, quartic symbol of the CM coordinate), J and E residue symbols of orders 2 and 4 (no order 8), P (histograms of the stored curve data, where stored).
- **B2, relational**: for the eleven absolute invariants v of note 28 (I₁…I₅, R₂…R₅, C²/p₁², the CM coordinate) the residue symbols of v(Y) − v(X) and v(Y) + v(X), of each order d ∈ {2, 3, 4} with d | q − 1.
- **F, per pair of lines** (population F): the 13 quantities of note 28 §2.2 (δ, u, M mod 16, v₂(M), order of the point, kernel shape, status, sign(Re c₁·Im c₁), sign(c₂), c₁, twist, r, s).

**Statistics.**

- *diff* (order-3 symbols, A): N_T·MI(f(Y) − f(X); c(Y) − c(X)) over the N_T = n(n − 1) ordered pairs X ≠ Y, in bits; a pair with a zero value gets its own value ('z', f(X), f(Y)). f(Y) − f(X) is a function of c(Y) − c(X) exactly when f = ±c + constant, that is when f is a character of order 3 (§2).
- *pair* (every surface candidate, A and B): N_T·MI((f(X), f(Y)); c(Y) − c(X)).
- *B2*: N_T·MI(h(X, Y); c(Y) − c(X)).
- *F*: R_T·MI(x(record); c(surface of the record)) over the R_T records of a torsor.
- Summed over the torsors where the candidate is defined and not constant.

**Null.** λ replaced by a uniformly random bijection onto G, independently on each torsor; c is then a uniformly random balanced labelled partition into three cosets. Exact: for n ≤ 12 all partitions are enumerated (15 for n = 6, 5,775 for n = 12, unordered; the statistics do not depend on the order of the cosets); for n = 24 (q = 43, e₁ = 2, if it converges) the surface statistics, which depend only on the contingency table of f against c, are computed exactly over all tables with their multivariate hypergeometric probabilities; B2 and F are not evaluated on a torsor of 24 surfaces, because their null cannot be enumerated (1.6·10⁹ partitions). A candidate whose table enumeration exceeds 5·10⁶ tables on a torsor is left out on that torsor. Combination over torsors by exact convolution on a grid of 0.01 bit-pairs, rounded to the conservative side, as in note 28. Candidates constant or undefined on every torsor enter with p = 1.

**Correction.** Holm over all candidates of a population, α = 0.05; Bonferroni and BH reported.

**Controls.**

- *Positive*: δ against the cosets of ker χ₋₁ on q = 37, e₁ = 2 and q = 41, e₁ = 22, with the same code (two cosets, all balanced 2-partitions, statistic n·MI(f; coset)). Expected p = (1/462)² = 4.7·10⁻⁶. Reported against the thresholds of populations A and B.
- *Shuffled labels*: population A rerun 10 times with the cosets of every torsor permuted at random (seeds 1 to 10).
- *Plants* on A and on B, 200 each, against the same nulls and the Bonferroni threshold 0.05/m of the population: an exact order-3 character on every torsor (diff); the same with each value replaced at random with probability 1/8 and 1/4; exact on 2 torsors only and on 1 only, random elsewhere; a two-valued function of the coset, f = [c = c₀] (pair); a random function h(c(Y) − c(X)) (B2, on the torsors with n ≤ 12); a random function (null).
- *Survivors*, if any: other torsors held out, the configurations separately, conjugate characters (the symbol 2f in place of f, which gives the same statistic by construction), relabelling (another primitive root), and other pairs of lines for F.

**Amendments, before the join.**

1. q = 29, e₁ = 2 also has the fourfold records of note 23 and the cosets of ker χ₋₁, which the list above forgot. Population F and the positive control take every torsor of A ∪ B with the data, as the rule says, so both include it; the expected p of the control becomes (1/462)³ = 1.0·10⁻⁸.
2. The q = 43 torsors are reported separately from q = 37. Besides population A as above, its two parts A₃₇ (q = 37) and A₄₃ (q = 43) are each scored with Holm over their own candidates, and each gets its own plants (the same eight kinds, 200 each) and its own power figures. Population A stays the primary test. If a q = 43 torsor does not converge, the note says at which precision and step it stopped.
3. **A debugging run.** The audit code was run once on a copy of the features of the three torsors of B then available, in which the cosets of every torsor (of order 3 and of ker χ₋₁) were permuted at random (seed 999), with 10 plants and 2 shuffles, to find errors in the code. That output says nothing about the true cosets and was not used. The table method was checked against the enumeration of partitions on q = 29, e₁ = 2 (25 candidates; largest difference in p, observed statistic and mean 3.7·10⁻¹³).

---

## 4. The run [E]

### 4.1 The torsors

Command `cosets`, one record per torsor in `sixlines_order3_cosets.jsonl`. "Primes" are the primes of E of degree 1 used for the class-field check: those whose class lies in ker χ₃ and where the cubic splits, and those where neither holds. No other case occurs.

| population | q | e₁ | G | surfaces | digits | chart, coordinate | primes: in ker χ₃ and split, neither | cosets of ker χ₋₁ | fourfold records | time |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 37 | 2 | Z/12 | 12 | 2,000 (1,200 failed) | J₂, (0, 1, 0, 0, 0) | 22, 40 | yes | 82 | 21 min |
| A | 37 | 6 | Z/12 | 12 | 2,000 (1,200 failed) | J₃, (0, 1, 0, 0, 0) | 10, 52 | no | — | 18 min |
| B | 23 | 2 | Z/6 | 6 | 700 | J₂, (0, 0, 1, 0, 0) | 12, 48 | no | — | 28 s |
| B | 29 | 2 | Z/12 | 12 | 1,200 | J₃, (3, −5, −7, −2, −7, 1, 6) | 18, 42 | yes | 92 | 6 min |
| B | 29 | 6 | Z/6 | 6 | 700 | J₃, (0, 1, 0, 0, 0) | 16, 44 | no | — | 29 s |
| B | 41 | 2 | Z/12 | 12 | 2,000 | J₂, (3, 1, 2, 5, 7) | 16, 44 | no | — | 19 min |
| B | 41 | 22 | Z/12 | 12 | 2,000 | J₂, (0, 1, 0, 0, 0) | 16, 44 | yes | 180 | 16 min |

- On every torsor the roots of H modulo 𝔭 are the enumerated surfaces. H has three factors of degree \|G\|/3 over L, with the same partition at each of the three primes of L above 𝔭, and the class-field check passes.
- **What the class-field check checks.** The test at 424 primes compares `bnrclassfield` with the subgroup of Cl(E) built from the Smith form that `sixlines_cmmethod.e_setup` returns; it does not check that subgroup. The splitting of H into three factors over L does. The separate check confirmed, with its own class-field code, that the stored cubic fields are cyclic and unramified with norm group the preimage of the unique subgroup of index 3 of G on all seven torsors (at 4,569 primes), that 𝔭 is principal in all ten fields, and that the subgroup is right [R].
- At q = 29, e₁ = 2 the coordinate is a seeded combination of the invariants (note 24's search, with its acceptance rule). The CM coordinate of X3 and of B2 is a different function on the two torsors of A: J₆/J₂³ at e₁ = 2 and J₂J₄/J₃² at e₁ = 6.
- At q = 37 both runs at 1,200 digits failed in the same way. One coordinate gave a class polynomial with 12 simple roots modulo 𝔭, but the polynomials of the invariants were not recognised. Both runs at 2,000 digits succeeded.
- **New class polynomials.** In `sixlines_cmmethod.jsonl`, note 20 marks q = 37, e₁ = 6 as "skipped: 12 lattices per type"; q = 41, e₁ = 2 appears in no earlier record. So for these two fields the roots of H are the enumerated surfaces for the first time: 12 of 12 each ([note 27](27-hypothesis-ledger.md), H7).
- **No quadratic subfield in E.** For all ten polynomials of §3, including q = 43, E has no quadratic subfield and x² + x + 1 is irreducible over E (command `fields`). This is the fact §2 uses.

### 4.2 q = 43

| e₁ | G | lattices per type | where it stopped |
|---|---|---|---|
| 14 | Z/12 | 12 | 36 period points at 2,000 digits in 617 s. H for the first coordinate (chart J₂, (1, 0, 0, 0, 0)) was recognised; the run then sat in the full factorisation of the lcm of its denominators for the progress line (below) until it was stopped, about 5 h 50 min after its start. Next in the queue was 3,000 digits. |
| 22 | Z/6 × Z/2 | — | queued behind e₁ = 14 in the same job, never started |
| 2 | Z/24 | 24 | 72 period points at 3,000 digits in 3,521 s of wall time (7.0 h of processor time). No line was written after the points, so the run was in the first coordinate's recognition or in the same factorisation, presumably the latter, when it was stopped, about 5 h 20 min after its start. Next in the queue were 4,500 and 6,000 digits. |

Both jobs ran at 80–95 % of a core and were stopped by hand; neither wrote a record. **The stall was not in recognition.** In `Torsor._setup_chart` of `sixlines_cmmethod.py` the progress line factors the lcm of the denominators of H completely, and that line runs only after `rec.poly` has returned H. A rerun with that line factoring only over primes below 10⁶ (`patch_setup_chart` in `sixlines_order3.py`, §9; the acceptance rule is unchanged) finds a composite cofactor of 97 digits in that lcm at e₁ = 14, and of 95 and 113 digits at e₁ = 22, and takes 0.1 to 0.2 seconds per coordinate. So these runs say nothing about the cost of recognition or about how the precision grows with \|G\| (note 27, H17). The reruns and their results are in §7.4.

**The enumeration at p = 43** that these runs need was made for this note with the commands of note 26 §4 (`sixlines.py enumerate`, then `sixlines_cm.py --list`): 2,432,120 configurations of six lines in general position, 3,603 classes, 1,437 of them self-associated (25 minutes on 14 processes); 608 classes with ordinary irreducible P, in 24 polynomials, and in all 24 the number of surfaces equals the number of lattices (10 with a maximal order, 14 without). Note 27 (H12, H16) uses it.

### 4.3 The order of events

"Fixed before the join" can be checked against the files only in part, so here is what was done, in order. The times are those of the records, the logs and the files.

1. **About 14:33.** The coset job for q = 37 was started, before §3 was written. It wrote its 1,200-digit failure record at about 14:40. The cosets of q = 23, e₁ = 2; q = 29, e₁ = 6; and q = 29, e₁ = 2 were computed around this time; the coset record of q = 23, e₁ = 2 holds I₁ to I₅, five of the eleven scalars of B2, next to each coset.
2. **14:35 to 14:44.** §3 and its amendments were written.
3. **14:39.** A features file was written with the three torsors of B then available, for the debugging run of §3 (amendment 3). By the definition of item 7 below it was a join: it held the true cosets next to every candidate value on those three torsors, whatever it was used for. That the joined values were not looked at rests on my account. The file was overwritten at the join, so its time rests on this account too.
4. **14:43:00.** The coset job for q = 41 was started, 45 seconds before the last change to the script (14:43:45). So the runs at q = 37, e₁ = 2 and q = 41, e₁ = 22 used the earlier version of the script, and the other two runs at q = 37 and 41 the final one. The test code was last changed before any record with cosets at q = 37 or 41 existed.
5. **From about 14:47.** The q = 43 jobs were started (§4.2).
6. **About 14:59 to 15:24.** The coset records of q = 41 and 37 were written: the first at about 14:59 and 15:01, the last at 15:18 and 15:24.
7. **19:09.** The command `features` was rerun on all cosets. This was the join.
8. **19:10.** The command `audit` was run: 56 seconds on 10 processes, 8 minutes of processor time. This was the run.
9. **19:13 to 19:15.** `labelcheck` and `selfcheck` were run.

**Deviations from §3: none in the tests.** Three things to know:

1. **Population A is q = 37 only.** In this run the rule "those that converge" leaves A = A₃₇ (§4.2). The candidate scores of A and A₃₇ are therefore identical. Their plants differ only by seed (2103 and 2137), and A₄₃ has neither scores nor plants.
2. **Processes.** `audit` ran on 10 processes instead of its default 12. The output is deterministic, so only the time changes.
3. **One self-test was uninformative.** The script checks that the cubic invariance filter drops quantities of the wrong weight. On q = 37, e₁ = 6 that self-test used a configuration with J₂ = p₁ = 0, so it could not show a change (0 of 6 trials). On q = 37, e₁ = 2 it dropped J₂ and p₁ in 4 of 6 trials each (2 in 3 is expected) and kept J₃ in 6 of 6. The filter itself dropped none of the 184 cubic quantities on either torsor, and no symbol of order 2 or 4 on any of the seven.

---

## 5. Results [E]

### 5.1 Population A: q = 37, two torsors of 12, m = 2,839

Monomials are written out: a J-exponent a = (a₂, …, a₆), an E-exponent a = (a₁, …, a₅, b) for p₁^{a₁}⋯p₅^{a₅}C^b.

| family | variant | candidates (constant or undefined) | best candidate | statistic | exact p | Holm |
|---|---|---|---|---|---|---|
| J3, cubic symbols of J-monomials | diff | 40 (0) | J₄J₅J₆² | 46.7 | 0.039 | 1 |
| J3 | pair | 40 (0) | J₄J₅J₆² | 99.3 | 0.012 | 1 |
| E3, cubic symbols of power sums and conic | diff | 121 (0) | p₁p₂p₃p₄²C | 48.3 | 0.030 | 1 |
| E3 | pair | 121 (0) | p₁p₂p₃p₄²C | 78.5 | 0.032 | 1 |
| N3, natural quantities | diff | 22 (0) | V = ∏(t_s² − t_{s′}²) | 27.9 | 0.031 | 1 |
| N3 | pair | 22 (0) | σ₅σ₇ | 36.7 | 0.049 | 1 |
| X3, cubic symbol of the CM coordinate | diff | 1 (0) | | 14.6 | 0.78 | 1 |
| X3, the same and the cubic type | pair | 2 (0) | cubic type | 245.1 | 0.46 | 1 |
| S, note 19 | pair | 13 (3) | j4 | 26.4 | 0.042 | 1 |
| J, d = 2 | pair | 15 (1) | J₃J₅ | 57.4 | 0.0046 | 1 |
| J, d = 4 | pair | 255 (1) | J₃²J₅² (the same function) | 57.4 | 0.0046 | 1 |
| E, d = 2 | pair | 63 (0) | p₂ | 26.4 | 0.042 | 1 |
| E, d = 4 | pair | 2,047 (0) | p₂p₃²C² | 121.0 | 0.012 | 1 |
| P, curve histograms | pair | 11 (2) | parity of the count of c₁ = −4 + 4i | 17.5 | 0.17 | 1 |
| B2, d = 2 | | 22 (0) | I₃(Y) + I₃(X) | 9.4 | 0.10 | 1 |
| B2, d = 3 | | 22 (0) | I₅(Y) + I₅(X) | 16.4 | 0.092 | 1 |
| B2, d = 4 | | 22 (0) | CM coordinate, sum | 19.0 | 0.17 | 1 |

- **The count.** 2,839 = 80 + 242 + 44 + 3 (cubic, both variants) + 13 + 15 + 255 + 63 + 2,047 + 11 (note 28's families, pair) + 66 (B2), as §3 lists them. Seven are constant on both torsors and enter with p = 1.
- **Calibration.** 57 of the 2,832 defined candidates have p ≤ 0.05 (2.0 %) and 4 have p ≤ 0.01. Discrete exact p-values are conservative.
- **The four smallest p are one function, and it is the exact indicator of one coset [E, R].** They are J₃J₅ (d = 2), J₃J₄J₅ (d = 2) and their d = 4 squares, all with per-torsor p 0.0061 and 0.76. J₄ is a square on these surfaces, so J₃J₄J₅ has the same square class as J₃J₅; on these two torsors the four candidates give one partition of the surfaces. On q = 37, e₁ = 2 the quadratic symbol of J₃J₅ is 0 on the four surfaces of one coset of ker χ₃ and 1 on the other eight; its per-torsor p is 35/5,775, the floor for a two-valued function. On e₁ = 6 it is unrelated (p 0.76). This is a case §8 lists as one the test would not see, and it occurred. It is chance, as Holm says: exact separations of one coset on one torsor number 32 distinct functions on A and 41 on B, against 45 and 55 expected under the null [R].
- **No exact character of order 3.** None of the 368 (candidate, torsor) cells of the diff variant is one: f(Y) − f(X) is a function of c(Y) − c(X) on no torsor for any cubic symbol. 259 of the cells contain a surface where the quantity vanishes, 110 at e₁ = 2 and 149 at e₁ = 6, and such a cell cannot be a character; so only 109 cells could have been. The best of these agrees with an affine function of the coset on 9 of 12 surfaces (N3, σ₁σ₈, e₁ = 2) [R]. The plants of §6.2 have no zeros, so their power applies to zero-free symbols.

### 5.2 Population B: q ≡ 2 mod 3, five torsors (6, 12, 6, 12, 12), m = 2,448

| family | variant | candidates (constant) | best candidate | statistic | exact p | Holm |
|---|---|---|---|---|---|---|
| S, note 19 | pair | 13 (3) | j4 | 30.4 | 0.073 | 1 |
| J, d = 2 | pair | 15 (1) | J₆ | 81.5 | 0.0067 | 1 |
| J, d = 4 | pair | 255 (1) | J₂²J₄J₆² | 86.1 | 0.0037 | 1 |
| E, d = 2 | pair | 63 (0) | p₂p₄ | 56.2 | 0.013 | 1 |
| E, d = 4 | pair | 2,047 (0) | p₁²p₂p₃²p₄C² | 296.3 | 1.4·10⁻⁴ | 0.34 |
| P, curve histograms | pair | 11 (0) | histogram of c₂ | 384.9 | 0.12 | 1 |
| B2, d = 2 | | 22 (0) | R₅(Y) − R₅(X) | 21.4 | 0.10 | 1 |
| B2, d = 4 | | 22 (0) | R₄(Y) + R₄(X) | 35.0 | 0.12 | 1 |

No residue symbol of order 3, and no B2 symbol of order 3, exists on B, because 3 ∤ q − 1. Order 4 is defined at q = 29 and 41 only. 67 of 2,443 defined candidates have p ≤ 0.05 (2.7 %) and 19 have p ≤ 0.01. The smallest p and its conjugate (p₁²p₂p₃²p₄C², a = 232302, the same statistic) are taken up in §7.2. BH gives q ≥ 0.17 for every candidate.

### 5.3 Population F: per pair of lines, m = 13

The records are 92 at q = 29, e₁ = 2 (note 23), 82 at q = 37, e₁ = 2 (note 23) and 180 at q = 41, e₁ = 22 (note 26): 354 records. On most surfaces there are 6 records (12 at q = 41, e₁ = 22; 16 on one or two surfaces, 48 on one), and about half have status ok: 39 of 82, 37 of 92 and 86 of 180 [R]. The statistic is R_T·MI(x; coset of order 3 of the record's surface).

| quantity | exact p | Holm | per torsor (29/2, 37/2, 41/22) |
|---|---|---|---|
| u | 0.0061 | 0.080 | 0.23, 0.0064, 0.23 |
| status of the computation | 0.015 | 0.19 | 0.0088, 0.33, 0.33 |
| M mod 16 | 0.030 | 0.33 | 0.72, 0.0055, 0.26 |
| order of the point | 0.32 | 1 | |
| v₂(M) | 0.45 | 1 | |
| sign(c₂) | 0.47 | 1 | |
| kernel shape | 0.56 | 1 | |
| sign(Re c₁·Im c₁) | 0.89 | 1 | |
| c₁ | 0.97 | 1 | |
| δ | 0.999 | 1 | |
| twist, r, s | constant | 1 | |

No quantity survives. The two smallest come from different torsors: u from q = 37 and the status from q = 29. They are not a pattern across torsors.

---

## 6. Controls and power [E]

### 6.1 The positive control

δ against the cosets of ker χ₋₁, on q = 29, e₁ = 2; q = 37, e₁ = 2; and q = 41, e₁ = 22. The code is the same, with two cosets, all 462 balanced splits of 12, and the statistic n·MI(δ; coset). Per torsor p = 1/462 = 0.0022, the statistic is 12 bits, its maximum, and the combined p is 1.01·10⁻⁸. That is (1/462)³, as amendment 1 of §3 expected.

### 6.2 Plants

200 of each kind. "Found" means p ≤ 0.05/m: 1.76·10⁻⁵ on A and A₃₇ (m = 2,839), and 2.04·10⁻⁵ on B (m = 2,448). "Replaced at random with probability r" replaces a value by a uniformly random element of Z/3, so it is wrong with probability 2r/3: 1/12 and 1/6 in the two rows.

| plant | A: found, p ≤ 0.05, median p | A₃₇ (seed 2137) | B |
|---|---|---|---|
| an exact order-3 character on every torsor (diff) | 200, 200, 3.0·10⁻⁸ | 200, 200, 3.0·10⁻⁸ | 200, 200, 2.3·10⁻¹⁴ |
| the same, each value replaced at random w.p. 1/8 | 173, 200, 7.5·10⁻⁷ | 174, 200, 7.5·10⁻⁷ | 196, 200, 8.7·10⁻¹² |
| w.p. 1/4 | 76, 196, 5.8·10⁻⁵ | 88, 197, 3.4·10⁻⁵ | 161, 200, 2.1·10⁻⁸ |
| exact on 2 torsors, random elsewhere | 200, 200, 3.0·10⁻⁸ | 200, 200, 3.0·10⁻⁸ | 172, 200, 5.5·10⁻⁶ |
| exact on 1 torsor, random elsewhere | 4, 200, 1.1·10⁻⁴ | 4, 200, 1.2·10⁻⁴ | 10, 159, 1.7·10⁻⁴ |
| f = [c = c₀] on every torsor (pair) | 0, 200, 3.7·10⁻⁵ | 0, 200, 3.7·10⁻⁵ | 200, 200, 8.9·10⁻⁹ |
| h(c(Y) − c(X)), h random (B2 form) | 200, 200, 3.0·10⁻⁸ | 200, 200, 3.0·10⁻⁸ | 200, 200, 2.3·10⁻¹⁴ |
| a random function to Z/3 (null) | 0, 8, 0.53 | 0, 8, 0.60 | 0, 6, 0.55 |

- **The exact plants check the convolution.** On A the median is (1/5,775)² = 3.0·10⁻⁸, for two torsors with 5,775 balanced partitions each. On B it is (1/15)²·(1/5,775)³ = 2.3·10⁻¹⁴. On A "exact on 2 torsors" is the same plant as "on every torsor".
- **The indicator of a coset on A.** The partitions of 12 surfaces with a given block of four number 35. So the best possible p of a two-valued function on one torsor is 35/5,775, and on two torsors (35/5,775)² = 3.7·10⁻⁵, above the threshold. A two-valued invariant that picked out one coset exactly on both torsors at q = 37 would come out with p ≈ 4·10⁻⁵ and Holm about 0.1, and would not count as found. On B, with five torsors, it would be found every time.
- **A reading on one torsor** is out of reach on both populations (4 and 10 of 200), as in note 28.
- **Seeds [R].** Each figure in the table comes from one seed, and a replay plant by plant with exact p-values gives the same counts in all 24 rows. With further seeds of the same generator, "exact on 1 torsor" on A is found 6 to 15 times in 200 (8 seeds, mean 11); the 1/4 plant on A 81 to 97 times; and the 1/4 plant on B 141 to 156 times, with median p 1.4–2.6·10⁻⁷. At the true error rates of the plants (§1 item 7) the counts are lower. The conclusions of this section do not change.

### 6.3 Shuffled labels

Population A, cosets permuted at random on each torsor, seeds 1 to 10, all 2,832 defined candidates each time:

| seed | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | true cosets |
|---|---|---|---|---|---|---|---|---|---|---|---|
| smallest p | 0.0039 | 0.0013 | 0.0013 | 3.4·10⁻⁵ | 0.0010 | 7.2·10⁻⁵ | 0.0027 | 0.0022 | 0.0033 | 0.00053 | 0.0046 |
| p ≤ 0.05 | 68 | 153 | 102 | 345 | 218 | 354 | 109 | 108 | 37 | 232 | 57 |
| p ≤ 0.01 | 6 | 18 | 15 | 86 | 42 | 101 | 16 | 26 | 7 | 65 | 4 |

The true cosets are less extreme than all ten shuffles on the smallest p and on the count at 0.01, and than nine of ten at 0.05. Two shuffles (seeds 4 and 6) came within a factor of 2 to 4 of the threshold: with a few thousand correlated candidates on two torsors, the minimum of the null p-values is not far from 0.05/m. A survivor at q = 37 would therefore have needed the held-out checks of §3 before it meant anything.

---

## 7. After the fact [E]

None of this was in the pre-registration, and none of it changes §5 or §6.

### 7.1 Labels modulo 3, and the grid with χ₋₁

- **Labels modulo 3.** On the two Z/6 torsors with full labels from note 20 (q = 23, e₁ = 2 and q = 29, e₁ = 6), the partition by the label modulo 3 is the partition by the cosets from L: 2, 2, 2 surfaces on both (`labelcheck`). This had been seen before §3 was written, as §3 says; it is now in a log.
- **The grid with χ₋₁.** On q = 29, e₁ = 2; q = 37, e₁ = 2; and q = 41, e₁ = 22, the cross-table of the coset of ker χ₋₁ (notes 24 and 26, from the orbits of ideals above primes ℓ ≡ 1 mod 4) against the coset of ker χ₃ (here, from the factors over L) has 2 surfaces in each of its 6 cells. In Z/12 ≅ Z/4 × Z/3 this must hold. For two independent balanced partitions it would happen with probability 8,100/34,650 = 0.23 on one torsor and 0.013 on three. The two computations share the code for the lattices, the period points and the recognition of H, and nothing after that: one uses ideals acting on lattices, the other factors H over L. So this is a consistency check of the cubic route, not an independent derivation of it. For q = 41, e₁ = 22 the separate check of note 26 rebuilt H with its own recognition, and its own cubic class field gives the same three cosets [R].

### 7.2 The smallest p of B

The quartic symbol of p₁²p₂p₃²p₄C², on q = 41, e₁ = 2, takes these values on the three cosets. The values are discrete logarithms modulo 4 to the least primitive root, and 'z' marks a surface where the monomial is 0 in F_q.

| coset | values on its 4 surfaces |
|---|---|
| first | 1, 1, 1, 1 |
| second | 2, 2, 2, 2 |
| third | 0, 3, 3, z |

The coset is a function of the value, so the per-torsor p is 1/5,775, the floor. On the other torsors where it is defined it is not:
- q = 29, e₁ = 2: p 0.33;
- q = 29, e₁ = 6: p 1;
- q = 41, e₁ = 22: p 0.15;
- the two torsors at q = 37, which were not used to find it: p = 1 on both, the smallest possible statistic.

How often does this happen by chance? Of the 7,315 (candidate, torsor) cells of B on torsors of 12, 1,759 have a unique maximising partition. So about 1,759/5,775 = 0.30 cells are expected at the floor under the null. Candidates of order 4 come in conjugate pairs, so that is about 0.15 independent events. One event, seen as two cells, is what the "exact on 1 torsor" plants look like, and the plants show it is not found after correction (10 of 200). It is not a character of order 3: it takes five values (0, 1, 2, 3 and a zero), not three, and it does not repeat.

### 7.3 The table method

`selfcheck` compares the exact table distribution with the full enumeration of partitions, on 25 random candidates of each of the five torsors of 12. The largest difference in p, observed statistic or mean is 3.8·10⁻¹³. It had been run on q = 29, e₁ = 2 before the join. No torsor of 24 converged, so the table method was not needed in the run.

---

### 7.4 q = 43, after the run [E]

The reruns of §4.2 gave cosets on all three torsors at q = 43, and A₄₃ was then scored as Amendment 2 asks. No test was changed. This subsection came after §5 to §7.3 and changes none of their numbers.

**The cosets.** Logs `sixlines_order3_cosets_43.log`, `_43_14.log`, `_43_2.log`, `_43_2_mix.log` and `_43_2_coord.log`; the records are appended to `sixlines_order3_cosets.jsonl`.

| e₁ | G | status | digits | chart, coordinate | primes: in ker χ₃ and split, neither | time |
|---|---|---|---|---|---|---|
| 14 | Z/12 | ok | 3,000 (2,000: check failed) | J₂, (1, 1, 0, 0, 0) | 12, 50 | 37 min |
| 22 | Z/6 × Z/2 | ok | 3,000 (2,000: not obtained) | J₃, (1, 0, 0, 0, 0) | 16, 44 | 35 min |
| 2 | Z/24 | ok | 4,500 (3,000: H only, surfaces identified by the coordinate) | J₂, (0, 2, −1, 1, 0) | 12, 50 | 2 h 33 min of points, then 31 s |

- **e₁ = 14 at 2,000 digits.** H and four of the five invariant polynomials were right. The fifth, J₅²/J₂⁵, was not: on 11 of the 12 roots it differs from (J₃J₅/J₂⁴)²/(J₃²/J₂³), which the enumerated surfaces satisfy. So the check "the roots of H are the enumerated surfaces" failed. At 3,000 digits it passes, with the same partition.
- **e₁ = 2: no coordinate.** With 24 surfaces in F₄₃, a coordinate is injective modulo 𝔭 with probability 3.1·10⁻⁴ (for a random map). None of the eight fixed coordinates, and none of the 300 random ones that note 26's `setup_more` adds, gave 24 simple roots.
- **e₁ = 2: the coordinate used.** (0, 2, −1, 1, 0) was chosen because it is injective on the 24 enumerated surfaces; up to sign it is one of two with coefficients at most 2 (`cosets --mix`). Its H is recognised at 3,000 digits with 24 simple roots modulo 𝔭, but the polynomials of the other invariants are not.
- **e₁ = 2: the weaker check.** The surfaces are identified by the value of the coordinate (`--identify-by-coordinate`). The check becomes: the 24 roots of H modulo 𝔭 are exactly the coordinate values of the 24 enumerated surfaces. They are, and over L, H has three factors of degree 8. That is weaker than the check on the other nine torsors, which compares all five invariants. At 4,500 digits, from stored period points, the polynomials of all five invariants are recognised for this coordinate, and the full check passes with the same partition (`sixlines_order3_cosets_43_2_mix4500.log`); the fixed coordinates still fail there. `features --with-43` reads that record and writes a file identical byte for byte to the one audited. The choice of coordinate used the enumerated surfaces; H, L and the partition do not.
- **All three.** The partition is the same at every prime of L above 𝔭, and the class-field check passes at all 60 or 62 primes.

**The run.** No test was changed.
- **Commands.** `features --with-43 --out sixlines_order3_features_q43.json`; `selfcheck --features sixlines_order3_features_q43.json` (the torsors of 12, including e₁ = 14 and 22: largest difference 3.9·10⁻¹³); `audit --features sixlines_order3_features_q43.json --out sixlines_order3_q43.jsonl` on 7 processes, 88 s, log `sixlines_order3_q43.log`.
- **Reproduction.** B and F come out identical to §5.2 and §5.3 (2,448 and 13 p-values, equal to 10⁻¹²), and A₃₇ identical to §5.1.
- **Candidates at q = 43.** There are 271 per surface, because 4 ∤ 42 leaves no symbols of order 4. The cubic filter dropped nothing.

| population | torsors (surfaces) | m (defined) | smallest p | candidate | Holm | p ≤ 0.05 |
|---|---|---|---|---|---|---|
| A | 37: 12, 12; 43: 24, 12, 12 | 2,839 (2,832) | 3.3·10⁻⁵ | B2 sym d=2 of R2(Y) + R2(X) | 0.095 | 56 |
| A₄₃ | 43: 24, 12, 12 | 499 (495) | 1.15·10⁻⁵ | the same | **0.0057** | 21 |

- **Order 3 at q = 43.** There are 184 cubic symbols on 3 torsors, 552 cells, and 0 exact characters. The best diff p is 5.6·10⁻⁴ (E3, a = 011020, all of it on e₁ = 2), Holm 0.28; the best pair p is 0.0018, and the cubic type has p 0.61. Here 43 ≡ 7 mod 9, so the cube roots of unity are not cubes, and the reading R³ = c³·d of §2 would show.
- **Shuffled labels on A**, seeds 1 to 10: the smallest p ranges from 7.2·10⁻⁵ to 0.0039. The true cosets give 3.3·10⁻⁵, below all ten; that is the candidate below.
- **Positive control:** 1.01·10⁻⁸, as before.
- **Plants, 200 each.** A plant counts as found when p ≤ 0.05/m (1.8·10⁻⁵ on A, 1.0·10⁻⁴ on A₄₃).

  | plant | A | A₄₃ |
  |---|---|---|
  | exact order-3 character on every torsor | 200 | 200 |
  | replaced with probability 1/8 | 200 | 200 |
  | replaced with probability 1/4 | 198 | 188 |
  | exact on 2 torsors | 200 | 200 |
  | exact on 1 torsor | 62 | 132 |
  | indicator of one coset | 200 | 200 |
  | a random function of c(Y) − c(X) (torsors of 12) | 200 | 200 |
  | a random function (null) | 0 | 0 |

**The survivor of A₄₃, and how it breaks.**
- **What it tests.** B2 sym d=2 takes the quadratic symbol of R2(X) + R2(Y), R2 = p₂/p₁², over ordered pairs, against c(Y) − c(X). For a symmetric quantity that amounts to whether X and Y share a coset.
- **Where it is defined.** B2 is defined on the two torsors of 12 at q = 43; it is not scored at n = 24. It has p 0.018 on e₁ = 14 and 5.2·10⁻⁴ = 3/5,775 on e₁ = 22.
- **On e₁ = 22** the statistic takes its largest value over all partitions. Of the 13 defined pairs inside a coset with non-zero sum, none is a square; 2 sums are 0 and 3 are undefined. Across cosets the split is 20 squares and 20 non-squares.

The checks:
1. **Direction.** A reading would have one direction. The share of non-squares inside cosets is 13/13 (43/22), 14/18 (43/14), 7/14 (37/6) and 5/18 (37/2), against 0.42 to 0.57 across cosets. At q = 37 the per-torsor p-values are 0.076 and 1; on B the candidate has p 0.74.
2. **The held-out torsor: Z/24 at q = 43**, which the audit did not score for B2. A test was written down before the values were looked at (`sixlines_order3_holdout43_prereg.txt`; script `sixlines_order3_holdout43.py`):
   - T1: the statistic over pairs i < j, against 200,000 random balanced labellings: p = 0.018.
   - T2: the share of non-squares inside cosets, one-sided: p = 0.34 (0.451 against 0.437 under the null).

   T1 as written is not the survivor's statistic. Over i < j, how c(Y) − c(X) splits between +1 and −1 depends on the numbering of the surfaces (by root), and that split carries T1: 11 of the 15 zero sums fall at +1. With the survivor's own definition, over ordered pairs, p = 0.22.
3. **No other held-out torsor.** Below q = 43 the only torsor with q ≡ 1 mod 3 and 3 ∣ \|G\| is q = 37, e₁ = 30 (Z/3, three surfaces), where every pair lies across cosets. The next primes with such torsors (61, 67, 73, 79) have no enumeration.

So it does not replicate in any form that could be a reading of χ₃. It is of the kind §8 lists, like the smallest p of §5.1: a near-exact separation on one torsor of 12, and the plant "exact on 1 torsor" is found 132 times in 200 on A₄₃. It passes Holm in A₄₃ (0.0057) and not in A (0.095), the primary population (Amendment 2).


---

## 8. What this says [D, H]

**For the decision problem.** δ, the one reading found, is a sign attached to a 2-adic quantity of the fourfold and read over F_{q⁴} (notes 23, 26). An analogue for χ₃ would presumably have to read a 3-adic or cubic quantity [H]. On the K3 side, at q ≡ 1 mod 3, the obvious ones are cubic residue symbols of invariants of the six lines. They are tested here, and on the two torsors where they exist none reads χ₃, exactly or approximately, among the quantities whose cubic class moves by a non-cube between cosets (§2). At q ≡ 2 mod 3 they are all trivial (§2). So a reading of χ₃ there would need a quantity in F_{q²} \ F_q whose ratio to its conjugate has order divisible by 3, or a pairing of order 3 on the Prym variety over an extension where μ₃ is rational. Neither was attempted.

**What would have been seen.**
- An order-3 character read exactly on both torsors at q = 37, or with one value in twelve wrong, by a symbol that never vanishes there.
- Any function of the coset difference that holds on all torsors.
- On B, an invariant that singles out one coset.

**What would not.**
- A reading on one torsor only.
- An invariant that singles out one coset at q = 37 (§6.2); one such function occurred, at chance level (§5.1).
- A reading at q = 37 by a quantity R with R³ = c³·d, whose cubic symbol is constant there (§2).
- Weak biases.
- Anything outside the families.

**For the ledger** ([note 27](27-hypothesis-ledger.md)).
- **H43 (the search problem)** is narrowed again, now for characters of order 3: as characters on the two torsors at q = 37, and by separation of the cosets on all seven, with \|G\| ≤ 12. It stays open.
- **H42** concerns characters of order 2 and is not affected.
- **H7** gains two fields, 24 surfaces, with the roots of H equal to the enumerated surfaces.
- **H17**: the first q = 43 runs stalled in the progress line's factorisation of the denominators of H, after H had been recognised (§4.2), so they bear neither on the cost of recognition nor on the precision. On the reruns (§7.4), \|G\| = 12 and Z/6 × Z/2 needed 3,000 digits (2,000 failed), and for \|G\| = 24 the class polynomial of a coordinate chosen to be injective modulo 43 was recognised at 3,000 digits, and the other invariants at 4,500.

**For note 28 §7.** Its remark that characters of order 3 were out of its reach is answered here within the families: as characters on the two torsors at q = 37, and by separation of the cosets on all seven torsors with 3 ∣ \|G\| at q ≤ 41.

---

## 9. Check and limits

- **The check.** One separate re-derivation from the stored files, with its own code. It confirmed the algebra of §2; the stored cubic fields on all seven torsors at 4,569 primes, with 𝔭 principal on all ten fields and the subgroup right; H over L for the two Z/6 torsors; 8,880 cubic feature values, recomputed from the raw configurations; every statistic, plant median and shuffle of §5 to §7, by its own exact computation; the 424 class-field checks, the labels modulo 3, the grid with χ₋₁ and its probability 8,100/34,650; the two new fields; and that the run matches §3. Its corrections are in the text: the error rates of the plants (§1 item 7, §6.2), the order of events (§4.3), the cube roots of unity at q = 37 (§2), and the indicator of one coset (§5.1). No test result changed.
- **What the check did not cover:** H on the five Z/12 torsors; GRH; how §2 and §3 were written before the join; the script as it stood before its last change; the order-2 and order-4 features, which it scored from stored values and did not recompute; and the table method for n = 24, which was not used.
- **"Fixed first" rests on the files and on my account** (§4.3). The test code was last changed before any record with cosets at q = 37 or 41 existed. The features file of 14:39 was a join on the three torsors of B then available; that its values were not looked at, and that no invariant was tabulated against a coset between 14:44 and the join, rests on my account.
- **The cosets** rest on PARI's `bnrclassfield` and `nffactor` and on class groups computed under GRH (PARI without certification). The subgroup handed to `bnrclassfield` is the preimage of ker χ₃ in Cl(E). It contains the kernel of Cl(E) → G, the subgroup of note 13, Lemma 9, and it is built from the Smith form that `sixlines_cmmethod.e_setup` returns. The checks are the class-field test at 424 primes (which checks `bnrclassfield` against that subgroup), the splitting of H into three factors over L (which checks the subgroup), the agreement with note 20's labels on the two Z/6 torsors, and the grid with χ₋₁ on three Z/12 torsors.
- **The class polynomials** are recognised by note 20's rule, and the lattices and period points come from its code, unchanged (only the progress line of `_setup_chart` was later replaced, below).
- **The generator.** The cubic symbols use the least primitive root of F_q. Another generator maps the values by x ↦ 2x on Z/3, which leaves the diff test and the mutual information unchanged.
- **Two torsors with q ≡ 1 mod 3.** Population A is two torsors of 12 surfaces at one prime, where the cube roots of unity are cubes. With q = 43 it has five, three of them at a prime where they are not (§7.4, after the run). The power figures of §6.2 are for the two torsors at q = 37; §7.4 gives them for five.
- **Fourfold data** exist on three torsors, with 6 records on most surfaces and about half of them with status ok (§5.3).
- **The script in this repository** differs from the one that ran in its docstring and help texts (cross-references, the list of commands), in where `cosets` keeps its period points (`sixlines_order3_points/` next to the script), in `patch_setup_chart`, and in the option `--with-43` of `features`, added after the run. `patch_setup_chart` replaces `Torsor._setup_chart` of `sixlines_cmmethod.py` so that its progress line factors the denominators of H only over primes below 10⁶ and prints the sizes of the cofactors; it changes only that log line, and the acceptance rule is unchanged. It and the point cache were added at about 22:15, after the q = 43 jobs had started and after the last change of §4.3 item 4, and no record used by §5 to §7 was made with them. By default `features` leaves out q = 43, so it reads the same records as the run. On a fresh copy with the stored files, `fields`, `labelcheck`, `features`, `selfcheck` and `audit` reproduce `sixlines_order3_features.json` and `sixlines_order3.jsonl` byte for byte, and their logs up to the times [E].
- **The q = 43 runs** of this run were stopped (§4.2); the reruns are in §7.4. They append records to `sixlines_order3_cosets.jsonl`; `features` leaves q = 43 out unless given `--with-43`, so a rerun reproduces §5 to §7 from the stored records either way.

---

## 10. What this does not give

- A proof that no invariant of a six-line surface reads a character of order 3, or the class.
- Any test with q ≡ 1 mod 3 beyond q = 37, any torsor of 24, or any group with 9 ∣ \|G\|.
- A reading of χ₃ from quantities in F_{q²} \ F_q, or a pairing of order 3 on the Prym variety.
- Anything about cost at cryptographic size.
- A check of the cosets of ker χ₋₁ or of the labels of note 20, which are taken as given.

---

## Reproduce

From `experiments/`, with `PY=python`:

```bash
$PY -I sixlines.py enumerate --p 43 --out sixlines_p43.jsonl --procs 14
$PY -I sixlines_cm.py --list --out sixlines_cm_p43.jsonl sixlines_p43.jsonl > sixlines_cm_p43.log
$PY -I sixlines_order3.py fields > sixlines_order3_fields.log
$PY -I sixlines_order3.py cosets 23 2 1 1 --digits 700
$PY -I sixlines_order3.py cosets 29 6 1 1 --digits 700
$PY -I sixlines_order3.py cosets 29 2 1 1 --digits 1200
for a in "37 2 1 1 --digits 1200 2000" "37 6 1 1 --digits 1200 2000"; do $PY -I sixlines_order3.py cosets $a --workers 6; done > sixlines_order3_cosets_37.log
for a in "41 22 1 1 --digits 2000" "41 2 1 1 --digits 2000"; do $PY -I sixlines_order3.py cosets $a --workers 6; done > sixlines_order3_cosets_41.log
$PY -I sixlines_order3.py labelcheck > sixlines_order3_labelcheck.log
$PY -I sixlines_order3.py features > sixlines_order3_features.log
$PY -I sixlines_order3.py selfcheck > sixlines_order3_selfcheck.log
$PY -I sixlines_order3.py audit --workers 10 > sixlines_order3_audit.log
```

- **The `cosets` lines** are as the records show them (q, e₁ and the precisions). No log was kept of the three runs at q = 23 and 29, and the worker counts of those runs are not recorded. `cosets` appends to `sixlines_order3_cosets.jsonl` and keeps the period points it computes in `sixlines_order3_points/`, so that a second run at the same precision reuses them; `features` uses the last record with status ok for each torsor.
- **The enumeration at p = 43** (the first two lines) takes 25 minutes on 14 processes and is needed only for the q = 43 runs and for note 27 (H12, H16); its timing log is `sixlines_enumerate_43.log`.
- **q = 43.** The first coset runs were `cosets 43 14 1 1 --digits 2000 3000` and `cosets 43 22 1 1 --digits 2000 3000` in one job, and `cosets 43 2 1 1 --digits 3000 4500 6000 --workers 8`. They need the enumeration at p = 43 (`sixlines_p43.jsonl`, `sixlines_cm_p43.jsonl`). The reruns of §7.4 were:
  ```bash
  for a in "43 14 1 1 --digits 2000 3000" "43 22 1 1 --digits 2000 3000"; do $PY -I sixlines_order3.py cosets $a --workers 6; done > sixlines_order3_cosets_43.log
  $PY -I sixlines_order3.py cosets 43 14 1 1 --digits 3000 4500 --workers 6 > sixlines_order3_cosets_43_14.log
  $PY -I sixlines_order3.py cosets 43 2 1 1 --digits 3000 4500 --workers 7 > sixlines_order3_cosets_43_2.log
  $PY -I sixlines_order3.py cosets 43 2 1 1 --digits 3000 --mix 0 2 -1 1 0 --workers 4 > sixlines_order3_cosets_43_2_mix.log
  $PY -I sixlines_order3.py cosets 43 2 1 1 --digits 3000 --mix 0 2 -1 1 0 --identify-by-coordinate --workers 4 > sixlines_order3_cosets_43_2_coord.log
  $PY -I sixlines_order3.py cosets 43 2 1 1 --digits 4500 --mix 0 2 -1 1 0 --workers 4 > sixlines_order3_cosets_43_2_mix4500.log
  $PY -I sixlines_order3.py features --with-43 --out sixlines_order3_features_q43.json > sixlines_order3_features_q43.log
  $PY -I sixlines_order3.py selfcheck --features sixlines_order3_features_q43.json > sixlines_order3_selfcheck_q43.log
  $PY -I sixlines_order3.py audit --workers 7 --features sixlines_order3_features_q43.json --out sixlines_order3_q43.jsonl > sixlines_order3_q43.log
  $PY -I sixlines_order3_holdout43.py > sixlines_order3_holdout43.log
  ```
  The last three `cosets` runs read the period points that the third one stored in `sixlines_order3_points/`. The record of e₁ = 2 that `features --with-43` reads is the last one, at 4,500 digits.
- **Times.** `features` takes 12 seconds on one process. `audit` takes 56 seconds on 10 processes and is deterministic, with fixed seeds. `labelcheck` and `selfcheck` take a few seconds. The coset runs take from 30 seconds to 21 minutes each (§4.1).
- **Files written.** `sixlines_order3_cosets.jsonl`, `sixlines_order3_features.json`, `sixlines_order3.jsonl` (12,950 records: candidates, shuffles, plants, the control), and the logs named above.
- **Files read.** `sixlines_cm.jsonl`, `sixlines_cm_p41.jsonl`, `sixlines_cm_p43.jsonl`, `sixlines_groupsize.jsonl`, `sixlines_cmmethod.jsonl`, `sixlines_p{23,29,37,41}.jsonl`, and everything that note 28's `build_population`, `prym_tables` and `fourfold_records` read.
- **Code imported, not changed.** `sixlines_cmmethod.py`, `sixlines_cm.py`, `sixlines_delta_labels.py`, `sixlines_characters.py` and `sixlines_leakage.py`.

## Sources

No new source was read. In this repository: notes [13](13-hardness-on-paper.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [26](26-delta-is-the-character-at-q-41.md), [27](27-hypothesis-ledger.md), [28](28-leakage-audit.md). PARI/GP's `bnrclassfield` as documented in the PARI manual [M].
