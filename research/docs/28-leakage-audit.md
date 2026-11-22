# 28. A leakage audit of surface invariants against the CM labels

*The leakage audit that [note 27](27-hypothesis-ledger.md) §5 prepares: whether any invariant computable from a six-line surface carries information about its label in G beyond the known bit χ₋₁. It joins the labels of [note 20](20-cm-method-and-the-action.md) §7.5 and the cosets of [notes 24](24-the-form-mod-4-is-the-character.md) and [26](26-delta-is-the-character-at-q-41.md) to the invariants of [note 19](19-quadratic-characters-from-a-surface.md), the J-invariants of note 20, the curve data of [notes 14](14-fourfold-of-a-six-line-surface.md) and 26, the fourfold data of [note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) and 26, and two families of residue symbols added here, and tests them with exact permutation tests fixed before the join. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_leakage.py`; **[H]** heuristic; **[R]** re-derived in a separate re-derivation check. Statements of other notes keep the label they have there. Checked once by a separate re-derivation: a re-implementation of the whole audit with its own code, from the raw configurations, gives the stored p-value, statistic and per-torsor values of every one of the 11,417 candidates, and "no further leak found" stands. What the check confirmed is marked [R]; its corrections concern what the audit could see, not what it found, and are made in the text (§1, §2.2, §3, §5, §7); what it did not cover is in §8. Notes 19, 20, 23 and 24, which supply the data, have no re-check beyond note 25's check of notes 20, 23 and 24; note 26 has its own check.*

---

## 1. Result

**No further leak was found within the families tested [E].** On the 9 torsors with full labels in G (56 surfaces, groups Z/4, Z/6, Z/4 × Z/2 and Z/8, none with the character χ₋₁), none of 4,838 pair functions built from surface invariants is a function of the label difference beyond chance: the smallest exact p-value is 0.0038, and after the Holm correction every candidate has p = 1. On the 5 torsors where only the cosets of ker χ₋₁ are known (60 surfaces), the one candidate that survives the correction among 6,579 is the positive control δ (p = 4.5·10⁻¹⁴, Holm 3.0·10⁻¹⁰); the next is at p = 0.0010, Holm 1. In detail:

1. **The labelled set [E].** Note 20 labels 67 surfaces in 15 fields: 14 torsors are fully labelled, with 66 surfaces, and the fifteenth (q = 37, e₁ = 18) has 1 of its 8 surfaces labelled. 5 of the 14 have 2 surfaces and carry no information for any test here, which leaves 9 torsors and 56 surfaces (§2). The 5 torsors of notes 24 and 26 have cosets only. Every stored value that could be compared was recomputed from the configurations and agrees: 1,488 values of note 19 on 196 configurations, the c₁ lists of 196, the CM invariants of 252, δ of all 60 surfaces from the per-pair records. A separate recomputation of 308,648 feature values from the raw configurations agrees as well, and no label leaks into the features [R].
2. **The tests were fixed before the join [D].** Six families of candidates, the pair functions formed from them, the statistic (mutual information between the pair value and the label difference, summed over torsors), the null (all relabellings of each torsor, enumerated exactly), the correction (Holm over all candidates of a population) and the controls were written down before any invariant was tabulated against a label (§3). One amendment before the join: no residue symbols of order 8 of the new family E (65,535 monomials would have dominated the correction).
3. **Characters of any order that a symbol can see [E].** For a residue symbol f of order d, "f(Y) − f(X) is a function of g" is the same as "f = χ∘λ + constant for a homomorphism χ : G → Z/d". 1,360 such symbols of order 2 and 4 (J-monomials, power sums and the conic, note 19's bits, the CM coordinate) were tested this way on the 9 torsors: 34 are defined on all 9, and 1,278, all of order 4, on only 2 or 3 [R]. 1,030 more were formed and are constant, or differ between a configuration and its associate (1,026 of them symbols, and the four parity bits of P). None is a character beyond chance (best p 0.0038, Holm 1). A symbol of order 4 of an F_q-quantity exists only for q ≡ 1 mod 4, which is 3 of the 9 torsors (29/10, 37/22, 29/6), and they include only 2 of the 7 torsors with a character of order 4. A quartic symbol that read an order-4 character exactly on both would get p = (1/3)(1/315) = 1.1·10⁻³, or 1.1·10⁻⁴ if it also read the order-2 character of 29/6, far above the threshold 1.0·10⁻⁵; every one of the 3,578 quartic J- and E-candidates has a smallest attainable p above the threshold [R]. So characters of order 4 were tested only through their squares. Note 20 §7.6's three square classes come out as there: J₂J₆ is a character on one torsor (the Z/4 of q = 23, e₁ = 10), J₃J₅ on none; the square class of J₂, which note 27 could not recompute, is a character on none of the 9.
4. **Relational tests on pairs of surfaces [E].** Legendre and quartic symbols of v(X) − v(Y) and v(X) + v(Y), for eleven absolute invariants v: 44 candidates, best p 0.064. The ordered pair (f(X), f(Y)) of every surface invariant, which catches a leak that is not a homomorphism: best p 0.0062, Holm 1.
5. **The fourfold side per pair of lines [E].** On the coset torsors, over 774 records of pairs of lines: u (p 0.70), M mod 16 (0.43), the valuation of M (0.44), the order of the point (0.77), the shape of the 2-primary kernel (0.77), the status of the computation (0.31), c₁ (0.015), the signs of Re c₁·Im c₁ (0.12) and of c₂ (0.023). None survives the correction. Only δ, which combines u with the two signs, separates the cosets.
6. **Nothing on the K3 side reads χ₋₁ [E].** On the coset torsors no candidate other than δ separates the cosets: note 19's conic, j2, j4 and quartic bits (p 0.17 to 1), the quartic type (0.096), the J symbols of orders 2, 4 and 8 and the E symbols of orders 2 and 4 (best 0.0010), the curve histograms (best 0.18).
7. **The positive control [E].** δ comes out first of 6,579 candidates on the coset torsors, with exact p = 4.5·10⁻¹⁴, the number note 27 §3.4 computed by hand (see §5.1), as a surface invariant and again per pair of lines. A test that missed it would have been too weak to mean anything; this one does not miss it.
8. **Power [E].** Planted signals scored against the same exact null and the same correction (m = 4,838, threshold 1.0·10⁻⁵) on the 9 torsors: an order-2 character read exactly on every torsor is found 200 times in 200 (median p 7.4·10⁻¹⁰), with each value replaced at random with probability 1/8 (about one value in 16 wrong) 165 in 200, with probability 1/4 (about one in 8 wrong) 61 in 200, and with values actually wrong one time in 8 or in 4, 63 and 1 in 200 [R]; read on 4 of the 9 torsors only, 37 in 200; on 2 of 9, 2 in 200. An order-4 character on the 7 torsors that have one: 200 in 200, but no candidate can produce that plant (item 3). A random function of g: 200 in 200. Random functions: 0 in 200 after correction, 11 in 200 at p ≤ 0.05 before (§5.2). Most candidates could never reject: 3,670 of the 3,802 defined candidates cannot reach the threshold even with their statistic at its maximum on every torsor, so the figures are for a correction about 20 times stricter than needed (§3) [R].
9. **Shuffled labels [E].** With the labels of every torsor permuted at random, the same 3,802 candidates give a smallest p between 0.00082 and 0.0044 in 10 runs, against 0.0038 for the true labels, which is less extreme than 8 of the 10 (§6).

**What this means for the search problem at small sizes [D, H].** Within these families, two surfaces of a torsor tell nothing about the class that relates them beyond what chance gives, on groups of order at most 8, except χ₋₁ through δ where it exists. The audit would have seen any order-2 character that one of the candidates reads on most torsors, and any exact relation f(X, Y) = h(g) among the relational candidates; characters of order 4 it saw only through their squares (item 3). It cannot see: a leak through a function outside the families; a leak on one or two torsors only; a character of order 3 (only the two Z/6 torsors have one, no residue symbol of order 3 exists there because 3 ∤ q − 1, and no statistic here could reach the threshold with one, §7; [note 29](29-characters-of-order-3.md) is the test of order 3); and the one leak that certainly exists, the label itself as a function of the complete invariants (J₂ : … : J₆), computed through the class polynomial of degree |G| (§7). The torsors are tiny, so none of this is evidence about the cost of the search problem at cryptographic size. It removes the cheapest kind of attack, a fixed residue symbol, at the sizes where it could be tested, and narrows H42 and H43 of the ledger without settling them.

---

## 2. The labelled set [E]

### 2.1 Torsors

A surface is a configuration and its associate (note 20 §2); every labelled surface here has two configurations, none is self-associated.

| population | q | e₁ | G | surfaces | used | χ₋₁ a character |
|---|---|---|---|---|---|---|
| L | 23 | 10 | Z/4 | 4 | yes | no (q ≡ 3 mod 4) |
| L | 29 | 10 | Z/4 | 4 | yes | no |
| L | 23 | 2 | Z/6 | 6 | yes | no |
| L | 31 | 18 | Z/4 | 4 | yes | no |
| L | 31 | 2 | Z/4 × Z/2 | 8 | yes | no |
| L | 31 | 14 | Z/8 | 8 | yes | no |
| L | 31 | 6 | Z/8 | 8 | yes | no |
| L | 37 | 22 | Z/8 | 8 | yes | no |
| L | 29 | 6 | Z/6 | 6 | yes | no |
| L | 17, 19, 19, 23, 31 | 2, 2, 6, 6, 22 | Z/2 | 2 each | no: one value of g | — |
| L | 37 | 18 | Z/4 × Z/2 | 8, one labelled | no | (in C) |
| C | 37 | 18 | Z/4 × Z/2 | 8 | yes | yes |
| C | 29 | 2 | Z/12 | 12 | yes | yes |
| C | 37 | 2 | Z/12 | 12 | yes | yes |
| C | 41 | 22 | Z/12 | 12 | yes | yes |
| C | 41 | 6 | Z/16 | 16 | yes | yes |

- **L**: 9 torsors, 56 surfaces, 112 configurations, from the last `action` record of each field in `sixlines_cmmethod.jsonl` (note 27's rule). With the 10 surfaces of the five torsors of order 2 and the one labelled surface of q = 37, e₁ = 18 that is note 20's 67.
- **C**: 5 torsors, 60 surfaces, 120 configurations; the label is the coset of ker χ₋₁ (the orbit of the ideals above primes ℓ ≡ 1 mod 4).
- No torsor of L has χ₋₁: q ≡ 3 mod 4 for 23 and 31, and for the three with q ≡ 1 mod 4 the character does not exist (note 19 §4.5 lists the prime-field torsors that have it: 29/2, 37/2, 37/18). So on L every leak would be "beyond χ₋₁", and the positive control has to come from C and from planted signals.

### 2.2 Which invariants exist where

| family | what | on | candidates |
|---|---|---|---|
| S | note 19's type, branch, branch_zero, conic, j2, j4, twist, c₁ values; the Legendre symbol of the CM coordinate | all | |
| S | quartic, quartic type, orientation; the quartic symbol of the CM coordinate | q ≡ 1 mod 4 (29, 37, 41) | |
| J | residue symbols of J₂^{a₂}⋯J₆^{a₆} with weight Σ k·a_k ≡ 0 mod d: d = 2 (15) | all | 15 |
| J | d = 4 (255) | q ≡ 1 mod 4 | 255 |
| J | d = 8 (4,095) | q = 41 (C only) | 4,095 |
| E | residue symbols of p₁^{a₁}⋯p₅^{a₅}C^b, p_k = Σ_s t_s^{2k} over the ten splits, C the conic invariant, weight condition as for J: d = 2 (63) | all | 63 |
| E | d = 4 (2,047) | q ≡ 1 mod 4 | 2,047 |
| P | histograms over the stored pairs of lines of c₁, c₂, `compare`, `rational_point`, c₄ rational; parity of the count of each of the four values ±4 ± 4i of c₁ | all (`sixlines_prym.jsonl`, `sixlines_prym_p41.jsonl`) | 11 |
| B2 | Legendre (and, for q ≡ 1 mod 4, quartic; for q = 41 on C, octic) symbols of v(Y) − v(X) and v(Y) + v(X), v in J₄/J₂², J₆/J₂³, J₃²/J₂³, J₃J₅/J₂⁴, J₅²/J₂⁵, p_k/p₁^k (k = 2..5), C²/p₁², and the CM coordinate of note 20 (which on 23/10, 31/18 and 37/22 equals J₄/J₂²) | all | 44 on L, 66 on C |
| F | per pair of lines: δ (control), u, M mod 16, v₂(M), order of the point, kernel shape, status, sign(Re c₁·Im c₁), sign(c₂), c₁, twist, r, s | C only: 774 records (242 of note 23, 532 of note 26), 356 with a value of δ | 13 |

- **The surface value.** For a function v of the canonical configuration (b, c, e, f), the surface value is the sorted pair (v(C), v(C′)) over the configuration and its associate. This is an invariant of the surface whether or not v is intrinsic. A symbol is used as a group element (the "diff" variant below) only on torsors where it agrees on the two configurations of every surface; the 1,024 E-monomials with an odd power of the conic never qualify. The reason is not association, under which the conic does not change sign, as was first written here: C(b, e, c, f) = λ²·C(b, c, e, f) on 126 of 126 surfaces, with λ the scaling of the J-vector [R]. The conic changes sign under odd relabellings of the six lines (the sign character of S₆), while J and the p_k do not. So where χ_d(−1) ≠ 1 (d = 2 at q ≡ 3 mod 4; d = 4 at q ≡ 5 mod 8), the symbol of a monomial with an odd power of C depends on the ordering that the canonical form picks [R]. As functions of the surface through the canonical forms these are still valid candidates for the pair variant.
- **Weights [D, checked].** Under l_i ↦ a_i·l_i and a change of coordinates g, each t_s = [ijk][lmn] is multiplied by μ = (∏a_i)·det(g)², so J_k by μ^k and p_k by μ^{2k}, and the conic by μ². A d-th power residue symbol of a monomial is therefore invariant under GL₃ × torus exactly when its weight is divisible by d; for a monomial with an odd power of C it is not invariant under relabelling (bullet above). The script checks every monomial kept on three configurations of every field under three random rescalings (0 dropped). That the filter does drop monomials of the wrong weight is not checked in the script or its logs; rerun independently, it holds: the symbols of J₃ modulo squares, and of J₂, p₁ and C modulo fourth powers, change under about half of random GL₃ × torus transformations [R].
- **The symbols** are discrete logarithms modulo d to the least primitive root of F_q, with a separate value for 0. The choice of generator is arbitrary but fixed for a field; it does not depend on the labels.
- **The J-vector** is computed from the configuration with note 20's `cm_tmat` and `cm_J`, and the stored invariants of all 252 configurations (chart J₂, or J₃ for q = 41, e₁ = 6 and for q = 29, e₁ = 6 and e₁ = 2) come out of it.
- **Note 19's invariants** are recomputed by `sixlines_characters.candidates` for every labelled configuration, including the 56 at q = 41 that have no record in `sixlines_characters.jsonl`; where a record exists (196 configurations) the 1,488 stored values and the 196 c₁ lists agree.

**A by-product for H19 [E].** The type and the branch count are constant on all 19 torsors, including the two at q = 41 that note 19 did not have (28 surfaces).

---

## 3. The tests, as fixed before the join [D]

These tests were fixed before the join, that is before any invariant was tabulated against a label. The file in which they were first written down is not part of this repository, so the pre-registration can be checked only against this text. What had been seen of the data then: the structure of the records, note 20 §7.5's printed table for q = 23, e₁ = 2 and note 26 §7's tables for q = 41. The text below follows the tests as then fixed. The counts of candidates were filled in after the run; two rules settled while writing the script, still before the join, are marked as such (the agreement rule for the diff variant in §2.2, and that a candidate is left out on a torsor where it is constant, which cannot change its p-value). The amendment at the end of this section is the only change to the list.

**Pair functions (population L).** For a candidate and two surfaces X ≠ Y of a torsor, with g = λ(Y) − λ(X):

- *diff*, for symbols and ±1 bits: f(Y) − f(X) in Z/d. It is a function of g if and only if f = χ∘λ + c for a homomorphism χ : G → Z/d [D: φ(x + g) − φ(x) = F(g) for all x forces F to be additive]. So this variant is the character test, for every character of order dividing d, not only of order 2.
- *pair*, for every surface candidate: the ordered pair (f(X), f(Y)). This catches a leak that is not a homomorphism.
- *B2*: the relational symbols.

**Statistic.** N_T·MI_T(pair value; g) in bits on each torsor, N_T = n(n − 1) ordered pairs, summed over the torsors where the candidate is defined and not constant.

**How the tests handle the ambiguity of labels.** A label is defined up to translation and up to the orientation x ↦ −x (note 20 §3). Translation leaves every g unchanged. Any automorphism of G, the orientation included, is a bijection of the values of g and leaves the mutual information unchanged. The null below is invariant under both. The character test is invariant because the set of homomorphisms is.

**Null.** λ is replaced by a uniformly random bijection from the torsor to G, independently on each torsor. The distribution is exact: all n! bijections are enumerated (n ≤ 8), and the torsors are combined by an exact convolution on a grid of 0.01 bit-pairs, rounded to the conservative side. A secondary statistic, the number of torsors on which the pair value is an exact non-constant function of g, has its p-value from the Poisson-binomial tail of the exact per-torsor probabilities.

**Population C.** The label is a bit. The statistic is n_T·MI(f(X); coset(X)) for surface candidates, N_T·MI(f(X, Y); same coset) for B2, and R_T·MI(x(pair of lines); coset) over the R_T records of a torsor for F. The null is a uniformly random balanced split of the surfaces, all of them enumerated (12,870 for 16 surfaces).

**Correction.** Holm over all candidates of a population together, α = 0.05; Bonferroni and Benjamini–Hochberg q-values reported. Candidates that were formed but are constant on every torsor, or never defined, enter with p = 1: 1,036 on L, 9 on C. So m = 4,838 on L and 6,579 on C. Most of these can never reject [R]. Of the 3,802 defined candidates of L, 3,670 cannot reach the threshold 1.0·10⁻⁵ even with their statistic at its maximum on every torsor, among them all 3,578 quartic J- and E-symbols; 225 cannot reach 0.01. On C, 3,380 of the 4,094 order-8 J-symbols cannot reach 7.6·10⁻⁶, because with 8 values on 12 or 16 surfaces they are nearly injective. Such candidates enlarge m without any chance of rejecting. Tarone's label-free count gives an effective m of 237 on L and 3,631 on C, so on L the correction is about 20 times stricter than needed (§5.2). This is the analogue of note 19 §4.3's count of 4,194,303 functions; the difference is that the labels allow a test of each function rather than a census of those that happen to be constant or balanced.

**Controls.** Positive: δ on C must come out with Holm-corrected p < 0.05, or the audit is declared too weak. Power: planted signals on L scored against the same null and m (§5.2). Null: random ±1 functions. Survivors, if any, to be attacked with held-out torsors, other pairs of lines, shuffled labels and the configurations separately.

**Amendment, before the join.** As first fixed, the tests of §3 also listed residue symbols of order 8 at q = 41 for the family E. There are 65,535 such monomials on two torsors; they would have multiplied m by ten and with it every threshold. They were dropped before the join. Order 8 is kept for J (4,095 monomials).

---

## 4. Results [E]

### 4.1 Population L: full labels, 9 torsors, m = 4,838

| family | variant | candidates (of which constant or undefined) | best candidate | statistic | exact p | Holm |
|---|---|---|---|---|---|---|
| S note 19 and the CM coordinate | diff | 6 (0) | j4 | 40.8 | 0.18 | 1 |
| S | pair | 13 (3) | quartic symbol of the CM coordinate | 84.8 | 0.14 | 1 |
| J residue symbols of J-monomials | diff | 270 (2) | J₂²J₃²J₅², d = 4 | 29.8 | 0.014 | 1 |
| J | pair | 270 (2) | the same | 43.2 | 0.014 | 1 |
| E power sums and conic | diff | 2,110 (1,024) | p₂p₄, d = 4 | 57.5 | 0.0038 | 1 |
| E | pair | 2,110 (0) | p₂p₄, d = 2 | 96.3 | 0.0062 | 1 |
| P curve histograms | diff | 4 (4) | — | — | 1 | 1 |
| P | pair | 11 (1) | parity of the count of c₁ = −4 + 4i | 169.4 | 0.086 | 1 |
| B2 relational symbols | | 44 (0) | quartic symbol of J₆/J₂³(Y) + J₆/J₂³(X) | 40.5 | 0.064 | 1 |

- **Calibration.** Of the 3,802 candidates that are defined, 200 have p ≤ 0.05 (5.3 %) and 4 have p ≤ 0.01 (0.1 %). The discrete exact p-values are conservative, and many candidates cannot reach 0.01 at all.
- **The best candidate**, p₂p₄ of order 4, is defined on the three torsors with q ≡ 1 mod 4, with per-torsor p 1, 0.038 and 0.10. Its conjugate p₂³p₄³ gives the same value, as it must.
- **Secondary statistic.** The candidate that is an exact character on the most torsors relative to chance is E d = 2, p₁p₂p₃p₅, on 1 of 9 torsors, p 0.086; Bonferroni 1. Over all 4,309 (candidate, torsor) cells of the diff variant, 60 are exact characters: 57 on torsors of 4 surfaces, where one balanced ±1 function in three is the character, 2 on torsors of 6 and 1 on a torsor of 8.
- **Note 20 §7.6 again, with J₂ itself.** On the 9 torsors: the square class of J₂ is a character on none (p 0.47); of J₂J₆ on one, the Z/4 of q = 23, e₁ = 10, as note 20 found (p 0.51); of J₃J₅ on none (p 0.96). J₄ is a square on every surface (constant).

### 4.2 Population C: cosets, 5 torsors, m = 6,579

| family | candidates (constant) | best candidate | statistic | exact p | Holm |
|---|---|---|---|---|---|
| S note 19, CM coordinate, δ | 14 (3) | **δ (positive control)** | 60.0 | **4.5·10⁻¹⁴** | **3.0·10⁻¹⁰** |
| S without δ | 13 (3) | quartic type | 60.0 | 0.096 | 1 |
| J, d = 2, 4, 8 | 4,365 (3) | J₂²J₃²J₄³J₅², d = 8 | 18.4 | 0.0013 | 1 |
| E, d = 2, 4 | 2,110 (0) | p₂p₃²p₄p₅²C², d = 4 | 35.5 | 0.0010 | 1 |
| P curve histograms | 11 (0) | parity of the count of c₁ = −4 − 4i | 13.5 | 0.18 | 1 |
| B2 | 66 (0) | octic symbol of R₅(Y) + R₅(X) | 61.4 | 0.012 | 1 |
| F per pair of lines | 13 (3) | **δ per pair (control)** | 348.4 | **4.5·10⁻¹⁴** | **3.0·10⁻¹⁰** |
| F without δ | 12 (3) | c₁ | 23.2 | 0.015 | 1 |

The candidates of note 19 on the coset torsors, which note 19 §4.5 looked at on two of them without labels: conic 0.17, j2 0.85, j4 1, quartic 0.69, quartic type 0.096, orientation 0.26, the Legendre symbol of the coordinate 1. The fourfold quantities: in §1 item 5. 353 of 6,570 defined candidates have p ≤ 0.05 (5.4 %), 40 have p ≤ 0.01 (two of them δ).

---

## 5. The positive control and the power [E]

### 5.1 δ

δ is +1 on one coset and −1 on the other on all five torsors (notes 24, 26). Under the null of a balanced split, the chance of a perfect agreement up to sign is 2/C(n, n/2) per torsor: 1/35, 1/462, 1/462, 1/462, 1/6,435. The product is 4.5·10⁻¹⁴, the number of note 27 §3.4, and the convolution returns it exactly, as a surface invariant and again over the 774 records of pairs of lines (the per-pair statistic is 348.4 bit-records rather than 356 because the records are not balanced between the cosets). It is first of 6,579, eight orders of magnitude below the first Holm threshold, 0.05/6,579 = 7.6·10⁻⁶. So the test does not miss the one leak that is known.

### 5.2 Power on population L

Each plant is a pair function on all 9 torsors, scored against the exact null; "found" means p ≤ 0.05/4,838 = 1.0·10⁻⁵, the Bonferroni threshold, which is the first step of Holm. 200 plants each, seed 2026.

| plant | found after correction | p ≤ 0.05 before | median p |
|---|---|---|---|
| an order-2 character, exact, on every torsor | 200 | 200 | 7.4·10⁻¹⁰ |
| the same, each value replaced at random with probability 1/8 (about one value in 16 wrong) | 165 | 199 | 2.4·10⁻⁷ |
| the same, probability 1/4 (about one value in 8 wrong) | 61 | 187 | 1.3·10⁻⁴ |
| an exact order-2 character on 4 of the 9 torsors, random elsewhere | 37 | 198 | 2.6·10⁻⁴ |
| on 2 of the 9 | 2 | 154 | 0.012 |
| an exact order-4 character, on the 7 torsors that have one | 200 | 200 | 7.5·10⁻¹² |
| f(X, Y) = h(g) for a random h : G → {0, 1} | 200 | 200 | 5.9·10⁻¹⁶ |
| a random ±1 function (null) | 0 | 11 | 0.56 |

The table is correct for the plants as defined, and an independent replay matches it [R]. Two rows need care:

- **The rows with errors.** A value "replaced at random with probability 1/8 (1/4)" is wrong with probability 1/16 (1/8). With values actually wrong one time in 8 and one in 4, the counts are 63 and 1 of 200, with median p 1.4·10⁻⁴ and 0.093 [R].
- **The order-4 row.** No candidate of L can produce this plant: quartic symbols exist only on the 3 torsors with q ≡ 1 mod 4 (§1 item 3), where an exact quartic reading gives p = 1.1·10⁻³ and is not found (0 of 20 independent plants) [R].

At the effective m = 237 of §3 (threshold 0.05/237), the plants of the second, third and fourth rows would be found 184, 110 and 97 times in 200, and the best candidate of L (p 0.0038) would still fail [R].

Read: a candidate that reads a character of order 2 on every torsor, or on most of them with few errors, would have been found. A leak present on fewer than half the torsors, or one that is only a bias, would mostly not. The one-torsor case is hopeless: on Z/8 a balanced function is a given character, up to sign, with probability 1/35.

### 5.3 Power on population C

The same on the 5 coset torsors, against m = 6,579 (threshold 7.6·10⁻⁶); 200 plants each, seed 2027, scored as surface candidates. These plants were added after the audit had run; they were not in the pre-registration.

| plant | found after correction | p ≤ 0.05 before | median p |
|---|---|---|---|
| the coset bit itself | 200 | 200 | 4.5·10⁻¹⁴ |
| each value replaced at random with probability 1/8 (about one value in 16 wrong) | 199 | 200 | 4.9·10⁻¹¹ |
| probability 1/4 (about one value in 8 wrong) | 160 | 199 | 1.2·10⁻⁷ |
| the coset bit on 2 of the 5 torsors, random elsewhere | 86 | 200 | 1.3·10⁻⁵ |
| on 1 of the 5 | 2 | 186 | 0.0045 |
| a random ±1 function (null) | 0 | 12 | 0.48 |

So a K3-side quantity that agreed with χ₋₁ on seven surfaces in eight would have been found on these torsors most of the time (173 times in 200 with values actually wrong one time in 8); one that agreed on three surfaces in four would mostly not (21 times in 200, median p 9.9·10⁻⁴) [R]; and one that read it on a single torsor would not.

---

## 6. Trying to break the best candidates [E]

Nothing survived the correction, so there is no survivor to break. The shuffled-label control was pre-registered for survivors; it was run anyway, as a check that the best p-values are what chance gives:

- **Shuffled labels.** The whole family of L was rerun 10 times with the labels of every torsor permuted at random (seeds 1 to 10). Each run repeats all 3,802 defined candidates against exact nulls:

  | seed | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | true labels |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | smallest p | 0.00082 | 0.0013 | 0.0044 | 0.0021 | 0.0019 | 0.00095 | 0.0016 | 0.0019 | 0.0038 | 0.0029 | 0.0038 |
  | candidates with p ≤ 0.05 | 180 | 123 | 66 | 189 | 111 | 300 | 153 | 155 | 210 | 92 | 200 |
  | with p ≤ 0.01 | 46 | 21 | 4 | 15 | 24 | 72 | 17 | 36 | 31 | 6 | 4 |

  Two of the ten runs were reproduced independently [R]. The true labels sit inside the spread of random ones on all three rows. The counts vary by a factor of five between runs because the candidates are strongly correlated (a monomial, its powers and its products with constants on a torsor move together), which is why the correction is Holm and not a count of small p-values.
- **The best candidates by torsor.** The best candidate of L, p₂p₄ of order 4, has its weight on one torsor (q = 37, e₁ = 22, per-torsor p 0.038; 0.10 on q = 29, e₁ = 6); on q = 29, e₁ = 10 its observed statistic is the smallest possible (p 1). The best of C without δ, an E-symbol of order 4, has per-torsor p 0.49, 0.18, 0.015, 0.38 and 0.0051. Neither pattern is that of a function that reads something uniformly; both are what the smallest of a few thousand p-values looks like.

The other checks planned for survivors (other pairs of lines, the two configurations separately) were not needed. The J symbols, and the E symbols with an even power of C, take the same value on both configurations of every surface, so the second would change nothing for them; all 23,312 E values that differ between the two configurations of a surface have an odd power of C (§2.2) [R].

---

## 7. What the audit says about the search problem at small sizes [D, H]

**What an attacker would try, and what is excluded.** The CSIDH-style attack on the decision problem reads a character from each of two objects; note 24 found one such reading, δ for χ₋₁. A cheaper attacker would look for a fixed formula: a power residue symbol of an invariant, or of a combination of two surfaces' invariants. On the 9 torsors with full labels, none of 4,838 formulas of that kind reads any character of order 2, or any function of g, beyond chance; characters of order 4 were tested only through their squares (§1 item 3). On groups of these sizes the audit would have found a formula that worked on most torsors. That is the content of the negative result.

**What it cannot exclude.**

- **The label is a function of the invariants.** (J₂ : J₃ : J₄ : J₅ : J₆) determines the surface, so the label up to translation is a function of it, through the class polynomial of degree |G| over E and its reduction modulo 𝔭, which is how note 20 computed the labels. Every near-injective candidate has the null distribution of an injective one, and the statistic has no power on it (three of the five histogram candidates of P have p = 1 on L for this reason; the histograms of c₂ and of `compare` have 0.55 and 0.48 [R]). The audit tests simple functions, not all functions.
- **Leaks on few torsors, or small biases** (§5.2).
- **Characters of order 3.** Only the two Z/6 torsors have one, and with q = 23 and 29 no symbol of order 3 exists. Nor can the pair variant catch a reading: for fixed g, X runs over every surface, so (f(X), f(Y)) is never an exact non-constant function of g; the random function of g in §5.2 is planted directly as the pair value, which only the relational candidates (B2) can produce. And every statistic of this audit is invariant under affine relabellings, so on a Z/6 torsor the per-torsor p is at least 12/720 = 1/60: a reading confined to the two Z/6 torsors has p ≥ 2.8·10⁻⁴, and an exact order-3 character there gives p = 0.0044, in either variant [D, R]. So a character of order 3 was out of reach of this test whatever the family. [Note 29](29-characters-of-order-3.md) is the test of it.
- **Characters of order 8, and of order 4 at q ≡ 3 mod 4.** No residue symbol of an F_q-quantity has that order there, so characters of order 8 on the three Z/8 torsors, and of order 4 on the torsors with q ≡ 3 mod 4, are not testable by symbols (§1 item 3).
- **Other families.** Pairings on the Prym variety other than δ, higher-degree invariants not in E, point counts of auxiliary curves beyond those stored, and anything over extensions of F_q were not tried. δ itself needs F_{p⁴}; an analogue for another character might need a larger field, and nothing here tests that.
- **Larger groups.** At |G| ≤ 8 the search problem is solved by listing the torsor. The audit is about structure, not cost. A formula that leaked one bit at these sizes would be expected to leak at every size, as δ appears to (H38 is open beyond five polynomials); a formula that failed here might work only on larger groups, for instance one that sees an element of order 16. Groups with |G| ≥ 9 and full labels need the CM method at larger precision (note 27, H17).

**For the decision problem.** On the coset torsors, nothing on the K3 side and nothing on the fourfold side other than δ sees χ₋₁. Note 19's conclusion that "nothing computed on the K3 lattice or from square classes is a character" extends to the J- and E-symbols and to full labels.

**For the ledger** (note 27): H19 holds on two more torsors; H20 and H21 are confirmed against labels, and the J₂ row of H21 is now computed; H23 and H24 are reproduced as the control; H40 is consistent; H42 and H43 stay open, with a negative result within stated families and stated power. H22, H26 and H38 are not affected. Note 27 records these results in its table (§2) and in its §6.

---

## 8. Check and limits

- **One separate check.** Besides the agreement with stored values (§2.2), the positive control, which the script reproduces to the digit, and the plants, the audit was re-implemented independently, with its own code, from the raw configurations: features, exact nulls, convolution, Holm, plants and shuffles. It gives the stored p-value, statistic and per-torsor values of all 11,417 candidates, the δ product 4.5025·10⁻¹⁴, the product 7.404·10⁻¹⁰ of the first plant, the tables of §4, §5.2 and §5.3 and two shuffled-label runs; 308,648 feature values recomputed from the raw configurations agree, and no label leaks into the features [R]. It is a check by the same kind of method, not an outside review.
- **What the check did not cover.** Note 19's own S invariants, which were compared with note 19's file by the audit only; the labels and cosets themselves; eight of the ten shuffle seeds; the byte-identical rerun of `features`; and the pre-registration, which has no file in this repository.
- **The exact null** is exact for the pre-registered statistic. The combination over torsors rounds to a grid of 0.01 bit-pairs and shifts the threshold by one bin per torsor to the conservative side. Planted exact characters give the product of the per-torsor probabilities (7.4·10⁻¹⁰ for the median plant of the first row of the table in §5.2 (an exact order-2 character on every torsor) is (1/3)³·(1/10)²·(3/35)·(1/35)³, the product over the nine groups), which checks the convolution.
- **The correction is Holm**, which is valid under any dependence and so conservative here: the candidates are strongly correlated (each d = 4 symbol and its conjugate give the same statistic, for one). Benjamini–Hochberg gives q ≥ 0.85 for every candidate except δ.
- **The population.** 9 torsors with full labels, all from note 20 and the single rule of note 27 for choosing records; none has χ₋₁ and none has more than 8 surfaces. The labels come from one computation per field.
- **The fourfold data** exist on the coset torsors only, and there on the stored pairs of lines of each surface, of which 356 of 774 records have a value of δ (most of the others have exponent 64, note 26). Nothing on the fourfold side was tested against full labels.
- **Canonical forms.** The P candidates depend on the stored computation of the curve on the canonical configuration. They are invariants of the configuration as functions, not as geometry.
- **Generator.** Symbols of order 4 and 8 depend on the choice of primitive root, which is fixed per field; a different choice applies an automorphism of Z/d to every value, which leaves the diff test and the mutual information unchanged.
- **Not tested**: joint functions of several bits (note 19 §4.3 did all functions of four), quartic or octic symbols of the ratios t_s/t_s′ one by one (only their normal form under relabelling, note 19's quartic type, is relabelling-invariant), and the two configurations separately for the S and P candidates.
- **The script grew after the audit ran.** The command `extra` was added after `audit` had written its files; `features` and `audit` were not changed (beyond a cross-reference in the docstring), and `features` rerun afterwards gives a byte-identical `sixlines_leakage_features.json`. The stored `sixlines_leakage_features.log` is from that rerun, not from the run that wrote the JSON.

---

## 9. What this does not give

- A proof that no surface invariant reads a character of G other than χ₋₁, or the class.
- Any test at |G| > 8 with full labels, or of the 75 polynomials with χ₋₁ not yet computed, or of non-maximal orders.
- Any test of fourfold quantities other than those stored by notes 23 and 26, or of the fourfold side against full labels.
- Anything about cost at cryptographic size.
- A check of the labels themselves: the audit takes notes 20, 24 and 26 as given.

---

## Reproduce

From `experiments/`, with `PY=python`:

```bash
$PY -I sixlines_leakage.py features > sixlines_leakage_features.log
$PY -I sixlines_leakage.py audit --workers 12 > sixlines_leakage.log
$PY -I sixlines_leakage.py extra --workers 12 > sixlines_leakage_extra.log
```

`features` reads the stored files and writes `sixlines_leakage_features.json` (4 seconds, one process, PARI). `audit` writes `sixlines_leakage.jsonl`, one record per candidate with its per-torsor statistics and p-values, and the plants; 2 minutes 9 seconds on 12 processes (19 minutes of processor time), numpy only. `extra` writes `sixlines_leakage_extra.jsonl`: the 10 shuffled-label reruns and the plants on C; 12 minutes on 12 processes (1.7 hours of processor time). All three are deterministic (fixed seeds). Files read: `sixlines_cmmethod.jsonl`, `sixlines_delta_labels.jsonl`, `sixlines_delta41_cosets.jsonl`, `sixlines_delta41_compare.json`, `sixlines_characters.jsonl`, `sixlines_p{17,19,23,29,31,37,41}.jsonl`, `sixlines_prym.jsonl`, `sixlines_prym_p41.jsonl`, `sixlines_jacobian_character.jsonl`, `sixlines_jacobian_character_extra.jsonl`, `sixlines_delta41.jsonl`; and the code of `sixlines_cmmethod.py` (`cm_tmat`, `cm_J`) and `sixlines_characters.py` (`candidates`, `forms`), which are imported and not changed.

## Sources

No new source was read. In this repository: notes [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [26](26-delta-is-the-character-at-q-41.md), [27](27-hypothesis-ledger.md), [29](29-characters-of-order-3.md). The method of Castryck, Sotáková and Vercauteren is as described in notes 19, 23 and 25.
