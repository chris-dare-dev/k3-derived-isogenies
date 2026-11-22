# 18. Re-derivation checks of notes 13 to 16

*Re-derivation checks of notes [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md) and [16](16-fourfold-over-f25-and-f49.md). It re-derives the statements listed in §2, reads the GP code of `sixlines_prym.py` line by line, reads four of the sources and the abstract of a fifth, reruns pieces of the four scripts, and adds five checks (`experiments/sixlines_review.py`). The corrections it finds are stated in their corrected form in notes 13 to 15. Its own [D] statements are not re-checked (§7).*

*Labels, all for work done in these checks: **[R]** a statement of the checked notes, re-derived here. **[E]** measured here, by `sixlines_review.py` or by a rerun of the script of a checked note where that is stated. **[P]** read here in the cited source. **[D]** derived here and not in the checked notes; not independently verified.*

---

## 1. Outcome

| Question | Answer |
|---|---|
| Any error that changes the conclusion (a CSIDH-type problem on abelian fourfolds with multiplication by Z[i], not a new source of hardness)? | **None found** |
| Note 13: Lemma 9, the shape P(±q) = 64q⁵·(square), Proposition Y | **Correct** [R]. The Artin–Tate formula and the squareness, quoted there from memory, were read at a source [P]. Two imprecise statements (findings 4, 5) |
| Note 14: translation of Paranjape, Lemmas 1 to 4, trace formulas, GP code | **Correct** [R]. The code computes what the text says. I found no way for the comparison to pass for a wrong reason. Both choices that note 14 calls "fitted" now have independent support (findings 8, 9) |
| Note 14 §5.2, sign − and reducible P | **Does not threaten the main claim** (§4). New: in every anomalous pair tested, 2,524 of one kind and 17 of the other, the curve gives the polynomial of the surface over F_{p²} [E] |
| Note 15 §2, why the map respects the action | A sketch with stated gaps. No error; one more gap (finding 3) |
| Note 15 §3, the dictionary between the two class groups | **Two errors** (findings 1, 2). The homomorphism the script checks is the square of the map of tori of §2; the formula of §3.3 was off by a square; "killed by 2" was false for the cokernel in 46 of 278 cases (all corrected in note 15). The orders, and the isomorphism of the odd parts, stand, and get a proof sketch that covers every polynomial |
| Note 16 | No error found. The regression on p = 13, 17 was rerun and agrees |
| Statements quoted from memory | Of eleven, six were read at a source or re-derived. Of the five that were not, two enter the reduction: Honda–Tate, and Deligne's description of ordinary abelian varieties with polarisations (§5) |

The standing of the conclusion after these checks is the standing the notes give it: the transfer to fourfolds is tested on Frobenius polynomials and not proved, and the equivariance of the map for the two group actions is a sketch that no computation tests.

---

## 2. What was checked, and how

| Note | Statement | What I did | Result |
|---|---|---|---|
| 13 | Lemma 9 (1) | Valuations at a ramified prime, as in the proof | holds [R]; only for odd primes, as stated |
| 13 | Lemma 9 (2) | Rewrote the double coset through s ↦ s/s̄, place by place | holds [R] |
| 13 | Lemma 9 (3) | Order of O(q_N) from Lombardo–Peters–Schütt, Remarks 1.5: Sp(4, F₂) × Z/2, 1,440 | holds [R, P] |
| 13 | P(q) = 64q⁵·\|Br\|, twist P(−t) | Source read; normalisation worked out; hypotheses checked for the surface and its twist | holds [P, R], §5 |
| 13 | e₂, e₃ in terms of (e₁, μ, ν); h(±2q) | Recomputed | holds [R] |
| 13 | Proposition Y (1), (2) | Recomputed the three coefficients of the second exterior power; the argument for Q(π) ∋ i | holds [R] |
| 13 | Proposition Y (3) | Recomputed | imprecise, finding 5 |
| 13 | §6, number of genus characters | Compared with the 2-ranks in `sixlines_groupsize.jsonl` | incomplete, finding 4 |
| 13 | `--real`, `--fourfold 41`, `--all 13 17 19 23` | Rerun (5 seconds) | same output: \|G\| = N in all 24; 21 polynomials at q = 41 |
| 14 | §2.1, what Paranjape proves | Read §§2 to 4 of the same scan | as reported [P], with two remarks (§5) |
| 14 | §2.2 steps 1 to 5 | Re-derived from Paranjape's §3 (below) | holds [R] |
| 14 | Lemmas 1 to 4 | Re-derived, including the existence of Q and the descent datum | hold [R] |
| 14 | §2.4 trace formulas | Re-derived from the fixed points of J^m·Frobenius^k, including branch points and poles | hold [R] |
| 14 | GP functions `pr_lambda`, `pr_curve`, `pr_half`, `pr_chain`, `pr_point`, `pr_findq`, `pr_trace` | Read line by line against §2 | compute what the text says |
| 14 | `--check 17 2 3 4 11`; `--summary sixlines_prym.jsonl` | Rerun; counts of §4 recomputed from the file | same output; counts agree except one (finding 6) |
| 14 | §5.2 | New measurements | §4 |
| 14 | §5.3, the other half | New check on p = 13, 17, 19 | the expectation holds, finding 8 |
| 15 | §2 | Read the lattice statement in Lombardo–Peters–Schütt; went through the steps | sketch stands; finding 3 |
| 15 | §3.1 claim and its justification | Recomputed the map on characters and the two composites; recorded the structure of kernel and cokernel | findings 1, 2 |
| 15 | §3.3 formula | Computed Ψ(𝔞) in PARI | finding 2 |
| 15 | §4 association | Re-derived the [D] statement; reran `association 5 7 11 13 17` | holds [R]; same counts |
| 16 | whole note | Read `sixlines_prym_fq.py` against `sixlines_prym.py`; reran `--run 13 17` | 140 and 618 pairs, 63/77 and 312/306, as stated |

**The translation of Paranjape (note 14 §2.2).** [R] The fibre of S → E over e is the line of pairs {x, e − x}, and the section D_i meets it at {p_i, e − p_i}. Two sections meet over e = p_i + p_j, so the cross-ratio of the four sections is a function of degree 2 on E, invariant under e ↦ −e because Σp_i = 0; this identifies E/± with the λ-line and gives step 3. For an isogeny E′ → E with kernel {0, κ}, the map E′/± → E/± is ramified at the points with 2x = κ, which lie over the two points of order 2 of E that are not the image of E′[2]; this gives step 2, with τ the image of E′[2]. Paranjape's text places the two fixed fibres over o′ and over "one of" the two points of order 2 of E′ that do not lie over o ("say q′"); so the pair (λ₅, λ₆) and the place of τ are as in step 2, and the paper does not say which of the two points it is. That is the choice between the two halves of τ (finding 8). Step 4: L² = O(Σp_i)·M and Σp_i = 0 give 2c₄ = τ. The value λ_c = N₁/D₁ is the image of the fixed points of the involution of the t-line that exchanges L5 and L6 and commutes with the involution of λ; on E′ that involution is the translation by a point of order 2 over τ, whose fixed points x have 2x over τ. So λ_c is a half of τ, which the script also asserts.

**Lemma 2.** [R] Changing the sign of one of a, b, c replaces {p_i} by −{p_i} up to a point of order 2 (with h′ = h − a the new points are −p₂, −p₁, −p₄, −p₃ in the labelling of the proof). So only on the twist of the lemma can the branch divisor be translated to a rational one, and the choice of twist is forced, not a convention.

**Lemma 3.** [R] 2Q = a + b + c + ρ in both cases; the sum of B is 2ρ and 4c″ = 4ρ, so 2B − 4(c″) − 4(o) is principal. Q exists because H¹ of the Galois group of F_{p²}/F_p with values in E(F_{p²}) is zero.

**Lemma 4.** [R] The map (x, v) ↦ (Frob x, G(Frob x)/v^p) lands on the curve because F·Frob(F) = G⁴, and its square is the Frobenius of v⁴ = F over F_{p²} because G has coefficients in F_p; so it is a descent datum. For a curve with no automorphism beyond the powers of J the type, commuting or inverting, is the same for every form over F_p [D]: the only isomorphism between (E, B) and its conjugate is the translation by the point ε of Lemma 2, and it fixes the class L; so the conjugate of L is L or L·M according to Frob(c₄) = ±c₄, whatever the form.

**Trace formulas.** [R] With J : v ↦ iv, the map J^m·Frobenius^k fixes the points over x with χ₄(F(x)) = i^{−m} when k is even, and those with ψ(x) = i^m when k is odd in the twisted case; in the second case v^{q+1} = F(x)^{(q+1)/4} does not depend on v, so all four points over x are fixed or none. Summing i^{−mj} times the number of fixed points over m gives T_k = −Σχ₄(F(x)) for even k and −Σψ(x)⁻¹ for odd k on the same eigenspace. At a branch point the contribution to the eigenspace of i is 0; at a pole of order 4 it is the character of the leading coefficient; both are as in the script.

**Could the comparison pass for a wrong reason?** I looked for these and found none.

- *Leakage.* The curve side receives only (p, configuration, pair). The numbers S₁, S₂, S₃ of the surface are read only in `compare`.
- *Selection.* A pair is computed or skipped by its type and by degeneracy, both decided before any point is counted; no pair is retried with another convention. For the 534 configurations with irreducible P no pair is of the skipped type: the stored records show 0 inverting pairs, 4,983 computed and 809 degenerate.
- *Loose comparison.* `compare` asks for three real integers equal to (e₁, e₂, e₃) or to (−e₁, e₂, −e₃). The constant term p² is an input of the comparison; a wrong one would give a mismatch, not a pass.
- *Conventions across fields.* The square roots are tied to one element (i for p ≡ 3 mod 4, a square root of the least primitive root otherwise). Replacing it by its negative applies a field automorphism to the whole computation and leaves the four counts unchanged [R]; so the three traces are consistent even though PARI chooses the root separately in each field.
- *Chance.* The other half of τ gives a quartic in Z[i][X] with the right absolute values and fails in 106 of 116 pairs (note 14 §5.3); so the test separates two curves that share E, the branch points and E′.
- *What the comparison cannot see.* χ against χ̄, π against −π, and the quartic twist. So it fixes the fourfold up to a quadratic twist and the choice of i, which is all the second exterior power can fix.
- *Coverage.* Mode B (p ≡ 3 mod 4, twisted Frobenius) is exercised only by note 14: 152 configurations with irreducible P and 1,846 pairs at p = 19, 23, 31. Note 16 has q ≡ 1 mod 4 and adds nothing for mode B.

---

## 3. Findings

Severity: **error** (a statement of the text is false), **gap** (a step is not justified), **imprecise**, **confirmation** (an open or fitted item of the checked notes that is now settled).

### 3.1 Finding 1. Note 15 §3.1: the map that is checked is not the map of tori of §2 — error

Write x_j for the image of x ∈ K under the embedding π ↦ π_j over Q(i), and ρ(x) = (x_j x_k)_{j<k}.

- [R] ρ does not map K^× into E^×. An automorphism that fixes π₁π₂ and moves i sends x₁x₂ to the conjugate of x₃x₄. So ρ(x) lies in E only if conj(x₃x₄) = x₁x₂, which for x·x̄ = m rational means N_{K/Q(i)}(x) = m². Note 15 §2 names the condition "x·x̄ rational"; the second condition is needed too, and it is the "determinant p²" of note 14 §5.2. In general ρ(x) lies in E(i).
- [R] The ideal Ψ(𝔓) = (ℓ, P_a(t)) of §3.1, which is what `sixlines_dictionary.py` computes for the primes over both square roots of −1, is induced by N_{E(i)/E}∘ρ : x ↦ x_j x_k·conj(x_l x_m), {l, m} the complement of {j, k}. On characters, the embedding of E at the pair {j, k} pulls back to σ_j + σ_k + σ̄_l + σ̄_m. On the subtorus where ρ has values in E, Ψ = ρ².
- [R] With the map in the other direction, Θ(y)_j = ∏_{k≠j} y_{jk}, the composites are

  > Θ∘Ψ(x) = (x/x̄)²·N·N̄², N = N_{K/Q(i)}(x), and Ψ∘Θ(y) = (y/ȳ)²·N_{E/Q}(y)².

  The composites x²·N and (y/ȳ)·N(y) are those of ρ, not of Ψ. So both composites are the *fourth* power on Cl⁻(K) and on G, not the square.
- [D] Consequences. The cokernel of Cl⁻(K) → G is killed by 4. If f(a) lies in H then a⁴ lies in the group generated by Cl(K⁺) and the images under Θ of the primes ramified over E₀, whose classes in Cl⁻(K) have order at most 2; so the kernel is killed by 8. Both are 2-groups, hence the odd parts of Cl⁻(K) and G are isomorphic, for every P with χ·χ̄ irreducible and maximal orders. This is more than the run shows, and it is a sketch.
- [E] `sixlines_review.py kernel`, the run of note 15 with the structure recorded:

  | | trivial | Z/2 | (Z/2)² | (Z/2)³ | Z/4 | Z/2 × Z/4 |
  |---|---|---|---|---|---|---|
  | kernel | 6 | 205 | 64 | 3 | 0 | 0 |
  | cokernel | 10 | 148 | 69 | 5 | 36 | 10 |

  The kernel is killed by 2 in all 278 cases. The cokernel is not in 46.

**What stands:** the Claim of §3.1 as stated (orders that are powers of 2), the two tables of orders, the isomorphism of the odd parts, and item 3 of §1. **What was corrected:** the paragraph "Why" of note 15 §3.1, where the composites are fourth powers and not squares, and where the cokernel is not killed by 2; note 15 §3.1 now states the corrected version.

### 3.2 Finding 2. Note 15 §3.3, §1 item 4 and §5.1: the image of 𝔞 was off by a square — error

- [E] For 𝔞 = 𝔓₁𝔓₂·conj(𝔓₃)·conj(𝔓₄) and Ψ as defined in §3.1, Ψ(𝔞) = (ℓ·𝔩₁₂/𝔩₃₄)², not ℓ·𝔩₁₂/𝔩₃₄: 48 primes out of 48, on 8 polynomials (`sixlines_review.py step`). 𝔞 has norm ℓ⁴ and Ψ(𝔞) has norm ℓ¹²; ℓ·𝔩₁₂/𝔩₃₄ has norm ℓ⁶.
- [R] By hand: Ψ(𝔓₁)Ψ(𝔓₂) = 𝔩₁₂²𝔩₁₃𝔩₁₄𝔩₂₃𝔩₂₄, and the two conjugate primes give the same again, since conj(𝔩_{jk}) is 𝔩 at the complementary pair.
- The displayed formula is true for ρ: N_{K/Q(i)}(𝔞) = (ℓ²) and 𝔞·conj(𝔞) = (ℓ), so 𝔞 is in the domain of ρ, and ρ(𝔞) = ℓ·𝔩₁₂/𝔩₃₄. So the sentence that matters for the reduction is right: the isogeny with kernel A[𝔞] corresponds to the step 𝔩₁₂/𝔩₃₄, the class of 𝔩₁₂ in G.
- What changes: f sends the class of 𝔞 in Cl⁻(K) to the *fourth* power of the class of 𝔩₁₂, not to its square. So the cokernel of f is explained only in part by the multiplier ℓ of 𝔞; a further factor 2 comes from f being induced by ρ² and not by ρ. In "For the reduction", the element of G attached to an isogeny with kernel A[𝔟] is ρ(𝔟), the square root of the ideal that the formula on primes gives. In §5.1, ρ(π) = F and Ψ((π)) = (F)².
- [D] Of the 16 ideals ∏_{j∈S} 𝔓_j·∏_{j∉S} conj(𝔓_j), all with 𝔞·conj(𝔞) = (ℓ), only the 6 with two primes of each kind are in the domain of ρ, one for each prime 𝔩_{jk}. The other 10 are isogenies of the same kind between principally polarised fourfolds with multiplication by Z[i]. I took them not to be steps of the surface and did not work out what they correspond to; [note 20](20-cm-method-and-the-action.md) finds that they are steps too. The group that note 15 §3.3 says it did not compute should be defined with this condition.

### 3.3 Finding 3. Note 15 §2, "From C to F_q" — gap

"The fourfold of the canonical lift has complex multiplication … A lift on which Frobenius acts is the canonical one." A lift with complex multiplication need not be one on which Frobenius acts: the quasi-canonical lifts have complex multiplication too. The step closes if the fourfold of the canonical lift is defined over the unramified ring W(F_q), where for odd p the only lift with complex multiplication is the canonical one [recalled, not read]; and that needs the six-line model of the canonical lift over W(F_q), which is the first gap note 15 lists. The four gaps that note 15 §2 listed when this check was made are real, and no other was found; note 15 §2 now lists this one as a fifth. Nothing in these checks tests §2.

### 3.4 Finding 4. Note 13 §6: the count of genus characters — imprecise

The bound "at most 3" counts characters attached to the ramified places. [D] A quadratic character ψ₀ of Cl(E₀) that is trivial on the norm of the ramified prime gives a character ψ₀∘N_{E/E₀} of G: it is trivial on the classes extended from E₀, whose norms are squares. These are of the same nature and their number is the 2-rank of Cl(E₀), up to one; it is not bounded by the number of ramified places. [E] In the 512 groups of note 13, G has more than one cyclic factor of even order only when h(E₀) is even (212 of 212), and at most one whenever h(E₀) is odd (258 of 258). So "their number is bounded whatever q" is not shown. The observed number is at most 4, so the conclusion that these characters give a few bits at most, none of them shown to leak, is unaffected in the data. Note 13 §6 now states the characters ψ₀∘N_{E/E₀} and this bound.

### 3.5 Finding 5. Note 13, Proposition Y (2) and (3) — imprecise

- (3): changing the sign of ν replaces χ(X) by χ̄(−X): π by −π *and* χ by χ̄.
- (2): the hypothesis "if its middle coefficient is prime to p" is automatic. [R] P ordinary has one unit root, so two roots of χ are units and the other two have the valuation of q, and χ·χ̄ has four unit roots.
- [R] A consistency that note 13 does not state: Proposition Y gives P(q) = 4q⁵·(Im a)² and P(−q) = 4q⁵·(Re a)² for a = the coefficient of X³. So the squareness of P(±q)/q⁵ follows from the existence of χ, and Artin–Tate adds that |Br| = (Im a/4)².

### 3.6 Finding 6. Note 14 §5.2: one count — imprecise

Note 14 §5.2 counted 673 of the 1,891 configurations with sign − that have pairs where Frobenius inverts J as also having pairs of commuting type: it is 672. 673 configurations with sign − have pairs of commuting type, and one of these has no inverting pair.

### 3.7 Finding 7. Note 14 §5.2: the account is incomplete, as note 14 says — see §4

### 3.8 Finding 8. Note 14 §5.3 and §7: the other half of τ — confirmation

Note 14 §5.3 expected, without checking it, that the other half gives the fourfold of the surface X′ in which L6 is replaced by its conjugate under the involution of λ. [E] `sixlines_review.py other-half 13 17 19`: for every configuration X and ordered pair, X′ was identified in the list, and the quartic computed with the other half was compared with the polynomial of X′.

| p | pairs computed with the other half, X′ of sign + | giving P or P(−t) of X′ |
|---|---|---|
| 13 | 178 | 178 |
| 17 | 688 | 688 |
| 19 | 1,034 | 1,034 |

So the half of τ is no longer a choice fitted on one surface: with one half the curve is that of X, with the other that of X′, and this was tested on surfaces that played no part in fixing it. It agrees with Paranjape's text, where the two lifts of the involution fix the fibres over {o′, q′} or over {p′, r′}: the same λ, λ₅, λ₆ and two surfaces. Note 15 §4 is the other combination: both lines conjugated gives the same surface.

### 3.9 Finding 9. Note 14 §7: the convention ψ⁻¹ — confirmation

Note 14 §7 settled the convention "the eigenspace sees ψ⁻¹ for odd k" on one surface at p = 19. It follows from the fixed points (§2, trace formulas) [R]. So neither of the two fitted choices rests on a fit any more.

### 3.10 Finding 10. Note 16 — no finding

The port to F_q follows `sixlines_prym.py` step by step. The embedding of F_q into the larger fields is one of two; the other gives the conjugate curve and the same counts [R]. The correction to the test for "P irreducible" is in both scripts. The two long runs (10 and 44 minutes) were not rerun.

---

## 4. Note 14 §5.2: sign −, reducible P, and whether it matters

### 4.1 What note 14 says

Sign + should go with a Frobenius that commutes with J and has determinant p² on the i-eigenspace, sign − with one that inverts J. In the data this holds for irreducible P. For reducible P and for sign − there are pairs of the other type, for which note 14 §5.2 refers to note 21, and nothing was computed in the inverting case.

### 4.2 New measurements [E]

**Pairs of commuting type on surfaces with sign −** (`sixlines_review.py sign-minus sixlines_prym.jsonl`; these are the 2,524 pairs on 673 configurations that note 14 records and does not compare).

- In all 2,524, c₃ = −p·c̄₁ and c₂ is purely imaginary. So the constant term of χ is −p² in every one, not only in the example of note 14. (c₁ = 0 in 1,380 of them; c₁ and c₂ never vanish together.)
- With that constant term, let χ₂ be the quartic whose roots are the squares of those of χ, the polynomial of Frobenius² on the eigenspace. Its second exterior power is the polynomial of the *twisted* surface over F_{p²}, P₂(−t), in all 2,524.

**Pairs of inverting type** (`sixlines_review.py inverting`). For p ≡ 3 mod 4 these are the pairs with c₄ rational: the curve v⁴ = F is defined over F_p, Frobenius² commutes with J, and its traces are the sums of mode A over F_{p²}, F_{p⁴}, F_{p⁶}.

| p | configurations | pairs | χ₂ has rational coefficients | second exterior power of χ₂ = P₂ |
|---|---|---|---|---|
| 7 | 2, sign − | 8 | 8 | 8 |
| 11 | 7, sign −; 2, sign + and reducible P | 9 | 9 | 9 |

Rational coefficients are what an inverting Frobenius forces (it exchanges the two eigenspaces and commutes with its square). The comparison is the test: three integers of size up to p⁶, one accepted outcome in these 17 pairs.

**The type against the signs of X and X′** (from the run of finding 8, p = 13, 17, 19; X′ as there).

| sign of X, of X′ | half of X: commuting / inverting | half of X′: commuting / inverting |
|---|---|---|
| +, − | 936 / 0 | 0 / 984 |
| −, + | 0 / 1,048 | 994 / 0 |
| +, + | 906 / 278 | 906 / 278 |
| −, − | 234 / 1,074 | 234 / 1,074 |

When X and X′ have different signs, the type follows the sign without exception, in 3,962 pairs: the half that belongs to the surface of sign + is of commuting type and the half that belongs to the surface of sign − is of inverting type. Every anomalous pair (commuting on sign −, inverting on sign +) has X and X′ of the same sign; by Lemma 1 these are the pairs for which the branch values b₁, b₂ are rational. This follows from the collinear triple of §4.3 and is not a criterion, since 46,334 other pairs have rational b₁, b₂ ([note 21](21-fourfold-forms-over-fp.md)).

### 4.3 Reading

- [P] In Lombardo–Peters–Schütt the moduli space is the quotient by SO*(T) = SO⁺(T) ∪ τ̃·SO⁺(T), where τ̃ has determinant −1, and on the fourfold side by U*((2, 2)) = U((2, 2)) ⋊ ⟨τ⟩, where τ carries a fourfold to an isomorphic one with the conjugate action. So over C the surface determines the fourfold with the ring Z[i], and not which generator is i; and a Frobenius of determinant −1 on the transcendental part corresponds to one that inverts J. This supports the first paragraph of §5.2.
- [E] So does the table above, in the 3,962 pairs where it can be read off.
- [D] In the anomalous pairs the form over F_p that Lemma 4 writes is a curve whose Frobenius is right over F_{p²}, up to the twist, and wrong over F_p by a root of unity, of order 8 on the eigenspace and of order 4 on the second exterior power ([note 21](21-fourfold-forms-over-fp.md)) (for the commuting pairs on sign −, its square has the eigenvalues of −F²). No twist of the curve by a power of J supplies that: such a twist multiplies the second exterior power by ±1. The remark after Lemma 4 in §2 (a curve without further automorphisms has no form over F_p of the other type) does not apply to these pairs: by [note 21](21-fourfold-forms-over-fp.md) (not re-checked), in every anomalous pair the curve does have automorphisms beyond the powers of J: the pair lies in a collinear triple (three of the intersection points of the lines on a line), one of the branch values b₁, b₂ is 0, 1 or ∞, and the curve has an involution γ with γJ = J⁻¹γ (checked directly in 15,973 cases). So the fourfold of the surface is the Prym variety of a form of the curve, the twist by γ. It is simple over F_p and over F_{p²} in 822 of the 1,144 commuting pairs with c₁ ≠ 0, and splits over F_{p⁴}, where γ is defined; so the explanation of note 14 §5.2, that the fourfold is not simple over F_p and has other forms, is not the one that applies in those pairs.
- [D] Why irreducible P is not exposed to this. If the curve is right over some extension, the second exterior power of its Frobenius commutes with F and a power of it is ± a power of F; when P is irreducible, with six distinct roots for that power, it is ω·F for a root of unity ω of E(i). An inverting Frobenius is excluded, since χ₂ would then be rational and P₂ would have the root p² twice. So ω = ±1, the two accepted outcomes, or another root of unity of E(i): ±i always, and more only if E(i) contains √2 or √−3. Any of those would show as a second exterior power that is not real, which the script records as "no"; there were none in 4,983 pairs over F_p and 4,182 over F₂₅ and F₄₉. So this is excluded by the data, not by the argument.
- [E] Size of the phenomenon. The anomalous surfaces are not a thin set over F_p. Of the 3,990 configurations 3,456 (87%) have sign − or reducible P, and all of those have Picard number at least 18 over the algebraic closure.

  | p | 13 | 17 | 19 | 23 | 29 | 31 | 37 |
  |---|---|---|---|---|---|---|---|
  | sign +, P irreducible | 4% | 7% | 7% | 11% | 12% | 13% | 15% |
  | sign +, P reducible | 58% | 50% | 45% | 42% | 41% | 39% | 37% |
  | sign − | 38% | 43% | 48% | 47% | 47% | 48% | 48% |

  Sign − is half of all configurations, as it must be if Frobenius has determinant −1 half of the time. Why reducible P is so frequent among sign + I did not establish. One ingredient [D]: if Frobenius on the transcendental part were i times the second exterior power of a linear map with constant term −p², then P(p) = 2p⁵·(Re c₁ + Im c₁)², which Artin–Tate (64p⁵ times a square) allows only if it is 0. So a surface whose Frobenius is of that kind must have the root p or −p.

### 4.4 Verdict

It does not threaten the main claim. The claim concerns ordinary surfaces with irreducible P. For those: every pair tested is of commuting type, the determinant is p², and the second exterior power is P or P(−t), over F_p, F₂₅ and F₄₉. The unexplained behaviour needs Picard number at least 18, and in every such pair tested the curve still gives the polynomial of the surface over F_{p²}. What it does show is that Lemma 4 is not a rule that produces "the form of the surface" in general: it produces the only form the curve has when the curve has no automorphisms beyond the powers of J, and that is the form of the surface exactly when the surface has one of that type; in the anomalous pairs the curve also has the twist by γ (§4.3).

Left open in note 14 §5.2 and by these checks, and taken up in [note 21](21-fourfold-forms-over-fp.md) (not re-checked):

- *The form over F_p of the fourfold* for surfaces of sign − and for the pairs of inverting type on sign +. For sign − it is the Weil restriction of an abelian surface over F_{p²}, with Frobenius² given by the polynomial of the surface, and for an anomalous pair it is the twist of the curve by an involution (§4.3).
- *Inverting pairs for p ≡ 1 mod 4*, where c₄ is not rational (`pr_trace` sets ω = 0 in that case). The descent is that of mode B, but the traces of odd powers of an inverting Frobenius vanish, so only Frobenius² matters; it is computed with a copy of `pr_trace` in which ω is set.
- *Inverting pairs for p ≥ 19*, which cost p⁶ operations each with three traces. χ₂ is rational and self-reciprocal, so two traces suffice and a pair costs p⁴ operations; 1,981 inverting pairs were computed for p = 7 to 37 (all pairs for p ≤ 17, samples above), and the third trace on 35 of them.

---

## 5. Statements quoted from memory

| Statement | Where it is used | Status after these checks |
|---|---|---|
| Artin–Tate: P(q) = \|disc NS\|·q⁵·\|Br\|, and \|Br\| is a square | note 13 §3.2; the shape of P; Proposition Y | **Read** [P]: Kedlaya–Sutherland, Theorem 8 (after Elsenhans–Jahnel and Milne): L₁(1) = \|Δ\|·#Br(X), #Br(X) a perfect square, for a K3 surface that satisfies the Tate conjecture. With their normalisation L(0) = q this is P(q) = 64q⁵·\|Br\| when NS over F_q has rank 16 and discriminant 64 [R]. Milne and Liu–Lorenzini–Raynaud themselves were not read; nor was a source for the Tate conjecture in odd characteristic |
| The hypotheses of that formula for the surface and its twist | same | [R] The sixteen classes are rational on both. The covering involution is +1 on them and −1 on the transcendental part, so the twist has P(−t). Discriminant 64: Theorem 3.4 of Lombardo–Peters–Schütt [P] and the fact that specialisation of Néron–Severi groups has no cokernel torsion prime to p [recalled, not read] |
| Honda–Tate | note 13, Proposition Y (2) | Not read |
| Deligne's description of ordinary abelian varieties, with polarisations | note 13 §4.3; note 15 §2 | Not read. Taelman's introduction confirms only that the theorem exists in the form used for K3 surfaces. The statement for polarisations and for Z[i]-lattices is the part of the dictionary I could not check |
| Nygaard's canonical lift; Taelman's theorem, unconditional for ρ̄ ≥ 12 and p ≥ 5 | note 15 §2; note 13 §9 | **Read** [P]: Taelman, Theorem B and the list after it |
| Paranjape, §§2 to 4 | note 14 §2.1 | **Read** [P]. The five bullets of note 14 §2.1 are in the text. Two remarks: the scan of Lemma 2 says "rank 2" over Q[i] where the proof gives Hyp ⊕ Hyp, rank 4, as note 14 has it; and the choice between the two fixed fibres is left open in the paper ("say q′") |
| The reconstruction of the lost formulas of Paranjape §2 | note 14 §2.1, last paragraph | Consistent with what survives [P]: §3 has L ≅ O(2·p) and M ≅ O(2·p − 2·q), an algebra O ⊕ L⁻¹ ⊕ L⁻¹⊗M ⊕ M, and J acting by −i on L and −1 on M |
| Lombardo–Peters–Schütt §§2.1, 2.2, Theorem 3.4, Proposition 4.4 | note 13 §4.1; note 15 §2 | **Read** [P]; as quoted. The lattice statement that note 15 §2 needs is their (12) and the lines after it: T ⊥ iT = Λ²V_R over R = Z[i]. Note that their quotients are by groups that contain τ and τ̃ (§4.3) |
| Lubicz–Robert: isogenies with kernel (Z/ℓ)^g by theta functions | note 13 §4.3, §5; note 15 §5 | Abstract read [P]: isogenies with a prescribed kernel for abelian varieties in theta coordinates, any dimension, odd characteristic. The cost in dimension 4 was not checked. It needs theta coordinates, which the Prym variety with its polarisation does not come with (note 15 §5.1 says so) |
| Kuperberg; Castryck–Sotáková–Vercauteren | note 13 §5, §6 | Not read |
| Harvey; Costa–Harvey–Kedlaya; Matsumoto; Matsumoto–Sasaki–Yoshida | note 13 §5; note 14 §6 | Not read |

The reduction (notes 13 §4, 14, 15) depends on the first eight rows, of which Honda–Tate and Deligne's description were not read. The last three rows enter only the comparison of costs with CSIDH and the remark on the reverse map.

---

## 6. Corrections to the checked notes

| Where | Correction |
|---|---|
| Note 13 §3.2, Brauer groups | The formula and the squareness were read at a source; hypotheses checked; no error |
| Note 13, Proposition Y | Item 3: χ̄(−X). Item 2: the hypothesis is automatic |
| Note 13 §6 | The count of characters omits those from Cl(E₀) |
| Note 14 §5.2 | 672, not 673; the measurements of §4.2 |
| Note 14 §5.3 | The expectation is confirmed on 1,900 pairs |
| Note 14 §7 | Neither fitted choice rests on one surface any more |
| Note 15 §1 | Items 3 and 4 corrected as in §3.1 and §3.3 |
| Note 15 §2, gaps | One more gap |
| Note 15 §3.1 | The paragraph "Why" replaced by the corrected statement; structure of kernel and cokernel |
| Note 15 §3.3 | Ψ(𝔞) is the square; what holds for ρ |
| Note 15 §6 | The structure of kernel and cokernel recorded |

Note 16 needs none.

---

## 7. Still unverified

- **The equivariance of surface ↦ fourfold** (note 15 §2). It is the step on which "reduces to" rests, and it remains a sketch with five gaps. No computation in the four notes or here tests it at a prime other than p.
- **Deligne's theorem with polarisations**, and with it the sentence of note 13 §4.3 that identifies self-dual hermitian Z[i]-lattices with principally polarised fourfolds.
- **The Prym variety against the principally polarised fourfold.** Notes 14 and 16 and these checks see only the isogeny class with its Z[i]-structure.
- **My own [D] statements**: the composites and exponents of finding 1; the ten other kernels of finding 2; the remark that the type is the same for all forms of a curve without extra automorphisms; the argument of §4.3 for irreducible P.
- **Propositions A to C of note 6**, on which Lemma 9 rests; and the observation of note 13 §3.2, that the admissible polynomials are exactly those that occur, beyond the rerun for q = 13 to 23.
- **Note 13 §3.3 and §3.4**: the table for q above 23, the heuristic for other families, the fits. Class groups are PARI's, correct under GRH, here as there.
- **The two long runs of note 16**, and the full run of note 14 (I recomputed its statistics from the stored file and reran one configuration).
- **Paranjape's proof**, and that his construction, made over C, describes the same surface in characteristic p. The evidence for the second is the point counts.
- **Inverting pairs** beyond the 17 of §4.2 (note 21 computes 1,981; not checked here).
- Five of the eleven quoted statements (§5).

---

## 8. Reproduce

From `experiments/`:

```bash
python -I sixlines_review.py sign-minus sixlines_prym.jsonl
python -I sixlines_review.py inverting 7 4
python -I sixlines_review.py step 8
python -I sixlines_review.py other-half 13 17 19 --workers 10
python -I sixlines_review.py kernel --workers 10
python -I sixlines_review.py inverting 11 1
```

A few seconds, half a minute, one minute, two and a half minutes, three minutes and eight minutes. The output of the six commands, in this order, is `sixlines_review.log`. They need `sixlines_prym.jsonl`, `sixlines_cm.jsonl`, `sixlines_groupsize.jsonl` and the files `sixlines_p{p}.jsonl`. The random points used for a few exceptional values do not affect the counts.

The reruns of the scripts of the checked notes quoted in §2:

```bash
python -I sixlines_prym.py --check 17 2 3 4 11
python -I sixlines_prym.py --summary sixlines_prym.jsonl
python -I sixlines_groupsize.py --cap 300 --fourfold 41 --real
python -I sixlines_groupsize.py --cap 300 --all 13 17 19 23
python -I sixlines_prym_fq.py --workers 10 --select all --pairs 0 --run 13 17
python -I sixlines_dictionary.py association 5 7 11 13 17
```

## Sources

- Paranjape, *Abelian varieties associated to certain K3 surfaces*, Compositio Math. 68 (1988), 11–22: §§2 to 4 read in the scan [P]; displayed formulas not legible.
- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, https://arxiv.org/abs/1209.5997 : §§1.1 to 1.3, 2.1, 2.2, Theorems 3.4 and 3.5, Proposition 4.4 [P].
- Kedlaya, Sutherland, *A census of zeta functions of quartic K3 surfaces over F₂*, https://arxiv.org/abs/1511.06945 : Theorems 1, 2 and 8 [P].
- Taelman, *Ordinary K3 surfaces over a finite field*, https://arxiv.org/abs/1711.09225 : introduction, Theorems A to C [P].
- Lubicz, Robert, *Computing isogenies between abelian varieties*, https://arxiv.org/abs/1001.2016 : abstract [P].
- In this repository: notes [2](02-six-line-count.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md).

