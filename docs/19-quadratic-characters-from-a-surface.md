# 19. Can a quadratic character of G be computed from a surface?

*Takes up the question that [note 13](13-hardness-on-paper.md) §6 and §9 leave open, whether quadratic characters of G can be computed from a surface. This is the decision version of the hardness question, as genus characters are for CSIDH (Castryck, Sotáková, Vercauteren). Builds on notes 2, 3, 6 and 13 to 18; the one character it leaves open is computed in notes 23 and 24. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_characters.py`. The derivations are not re-checked.*

Notation as in note 13: P the Frobenius polynomial on the transcendental part, E = Q[F] of degree 6, E₀ = Q(x) its real cubic subfield with x = F + q²/F, and for a maximal order G = Cl(E)/⟨Cl(E₀), ramified primes⟩, acting simply transitively on the set S_P of surfaces with polynomial P (note 13, Lemma 9). A function φ from S_P to ±1 *is a character* if φ(g·X) = χ(g)·φ(X) for a character χ of G of order at most 2. Such a φ is constant on S_P, or takes each value on exactly half of it.

---

## 1. Result

**No candidate of this note computes a character.** On paper, the method of CSIDH has nothing to read on the K3 side; on the fourfold side it could reach at most one character, which exists for 80 of 477 polynomials. This note does not compute it; notes 23 and 24 do, on three torsors (item 7). In detail:

1. **Which characters exist [D, E].** Every character of order 2 of G is the character of an unramified extension E(√d)/E with d in the real field E₀ (Lemma 1). On a step at a prime of degree 1 above ℓ it is the Legendre symbol of d at the corresponding root of h mod ℓ. Their number is at most 3 + (2-rank of Cl(E₀)); note 13 §6 said at most 3, and [note 18](18-checks-of-notes-13-to-16.md), finding 4, already pointed to the characters that come from Cl(E₀). On the 480 fields of note 13 with q ≤ 10,007 it is 0, 1, 2, 3, 4 for 3, 287, 160, 28 and 2 of them. For 437 of the 699 characters d can be taken a unit of E₀; the other 262 need an ideal class of even order of E₀. No torsion subgroup is attached to such a d.
2. **They are not genus characters in the sense of CSIDH [D].** The surfaces of one polynomial form a single genus: their lattices with F are locally isometric at every prime (Lemma 2, which holds because exactly one finite prime of E₀ ramifies in E). The characters of Castryck, Sotáková and Vercauteren separate genera, by a pairing on torsion at a ramified prime. Here there is nothing for them to separate. G is the analogue of the group of squares Cl(O)², not of Cl(O); a character of order 2 of G is the analogue of a character of order 4 of Cl(O), for which no method is known [M].
3. **The ten square classes are the same for all surfaces of a torsor [E].** For the surface εw² = l₁⋯l₆ the ten classes of ε[ijk][lmn] in F_q*/F_q*², one for each split of the lines into two triples, are well defined up to relabelling. Their normal form is the same for all surfaces of each of the 36 polynomials with maximal order (419 surfaces, torsors of up to 36), and is a function of q and e₁ modulo 16. So λ₅λ₆, (1 − λ₅)(1 − λ₆), the choice of twist of E, and the discriminant of the branch quadratic of note 14, which are all products of these classes or determined by them, carry no information about the class. Whether the half c₄ is rational depends only on q mod 4.
4. **Finer invariants vary, and are not characters [E].** The square classes of the conic invariant and of two power sums of the ten products, and the product of the ten quartic characters, are invariants of the surface (for q ≡ 1 mod 4) that take both values inside a torsor. On the torsors with at least two surfaces each of them is neither constant nor balanced in 14 to 24 cases. Among the 4,194,303 functions of the four together, one is constant or balanced on all 23 torsors with q ≡ 1 mod 4, and it is an artefact (§4.3). Tested against labels in [note 28](28-leakage-audit.md): on the 9 torsors of note 20 with full labels (56 surfaces), none of the four is a character or a function of the label difference; neither are the quartic type, the orientation or the c₁ values (p 0.15 to 1). On the 5 coset torsors of notes 24 and 26 none of them separates the cosets of ker χ₋₁ (p 0.096 to 1). The same holds for the residue symbols of order 2 and 4 of monomials in J₂, …, J₆ and in the power sums of the ten products and the conic, and of order 8 at q = 41 for the J-monomials only, with a correction over all candidates. As tests of characters the order-4 symbols could not reach significance: quartic symbols exist only on the 3 labelled torsors with q ≡ 1 mod 4 (note 28 §1 item 3), so the negative result is for characters of order 2, and for order 4 only through their squares. The type and the branch count are constant on all 19 labelled torsors, including two at q = 41.
5. **The unit of c₁ is not an invariant of the surface [D, E].** The surface determines its fourfold with multiplication by Z[i] only up to the quadratic twist and the conjugate structure, because −1 and conjugation act trivially on the second exterior power. In the stored quartics, c₁ takes two or more values among the pairs of one configuration in 263 of the 343 configurations that have two such pairs.
6. **Frobenius conjugation [D, E].** A square class in F_q of a quantity defined over F_p, a count, or a point count of an auxiliary curve is unchanged by Frobenius conjugation; all candidates above are (1,117 configurations over F₂₅ and F₄₉, no exception). The lattice side has a character with value −1 on the class of 𝔭/𝔭̄ in 9 of the 12 polynomials with maximal order over F₂₅ and F₄₉. That character cannot be computed by any candidate of this kind. It needs a quantity that Frobenius moves, such as a quartic character over F₄₉ with value ±i; the one symmetric choice I found is 0 on 148 of 248 surfaces and is not an invariant.
7. **One character has the form CSIDH's method computes; this note does not compute it, notes 23 and 24 do [D, E].** For d = −1 the character is "the degree ℓ of the step modulo 4". It is a character of G for 80 of the 477 polynomials with a character (never for q ≡ 3 mod 4; d = 2 and d = −2 never give one), and exactly for those whose real quartic field K⁺ on the fourfold side has a prime above 2 of odd degree (477 of 477): there the polarised fourfolds should fall into two genera. On the fourfold it is the multiplier modulo 4 of a polarised isogeny, which for elliptic curves is read off a self-pairing of order 4 [M]. Doing the same here needs arithmetic in the Jacobian of the curve of genus 5 and the unproved equivariance of note 15. None of the candidates of items 3 to 5 matches it on the torsors where it exists. [Note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) §5 computes it, with arithmetic in the Jacobian of the curve of genus 5. The invariant δ is the class modulo 4 of the hermitian form of the polarisation on the line of the 2-adic Tate module at the prime of degree 1 above 2, normalised by the unit of c₁ and the sign of c₂ so that the conjugate structure and the twists drop out; it is read off a Tate pairing of order 16 over F_(p⁴), not off a point of order 4. On the three torsors over prime fields of §4.5 (32 surfaces, 106 curves computed of 242 tried) it is the same for all pairs of lines of a surface, on 32 surfaces of 32, and it takes each value on exactly half of each torsor: 6 and 6, 6 and 6, 4 and 4. It is none of the candidates of items 3 to 5. That it is the character, δ(g·X) = χ₋₁(g)·δ(X), is not shown in note 23, where no surface of these torsors is labelled by a class; [note 24](24-the-form-mod-4-is-the-character.md) shows it on those three torsors. Against the orbits of the steps at primes ℓ ≡ 1 mod 4, computed by the CM method of note 20, δ is constant on each of the two orbits and opposite between them: 32 surfaces of 32. So for these polynomials χ₋₁ of the class that relates two surfaces is computed from the two surfaces.
8. **The class group of the field of the fourfold sees no quadratic character [D, E].** Note 18 (findings 1 and 2) showed that the homomorphism Cl(K) → G of note 15 is induced by the square of the map of tori. A consequence: every character of order 2 of G is trivial on its image (19,564 primes of K over 477 polynomials, no exception), because χ_d(Ψ(𝔓)) = (N_{E₀/Q}(d)/ℓ). So the cokernel of note 15 maps onto G/G², and the cosets of its image are not sets of surfaces that the fourfolds tell apart: the steps between them are the isogenies that multiply the polarisation by ℓ.

For the decision problem: nothing computed in this note distinguishes g·X from a random surface of the torsor. This does not hold where χ₋₁ is a character of G (80 of the 477 polynomials): the invariant δ of note 23, read on the fourfold, distinguishes the two cosets of ker χ₋₁ on the three torsors tested ([note 24](24-the-form-mod-4-is-the-character.md)), so there the decision problem leaks one bit, as for the class group actions in the paper of Castryck, Sotáková and Vercauteren (that paper says CSIDH itself is not affected). The claim was re-checked in [note 25](25-checks-of-notes-20-23-24.md): it stands on the three torsors, and δ equals a number read off the lattice of the CM lift (§4.4 there). Items 2 to 6 stand: nothing computed on the K3 lattice or from square classes is a character. The bits that the analogy with CSIDH would leak are constant on a torsor, and the characters that do exist have no known geometric reading, with the exception of χ₋₁ of item 7. This is evidence from p ≤ 37 and two fields of p² elements, and a negative result about a list of candidates, not a proof that no character is computable.

---

## 2. On paper

### 2.1 The characters of order 2 of G [D]

**Lemma 1.** Let χ be a character of order 2 of G, seen as a character of Cl(E) that is trivial on the image of Cl(E₀) and on the primes ramified over E₀. Then the corresponding quadratic extension of E is E(√d) with d in E₀*. Conversely, for d in E₀* with E(√d)/E unramified, the character of E(√d)/E is trivial on the image of Cl(E₀), and it is a character of G if and only if the ramified primes split in E(√d).

*Proof.* For an ideal 𝔞 of E, 𝔞·𝔞̄ is extended from E₀, so χ(𝔞̄) = χ(𝔞). The class field L of χ is therefore stable under complex conjugation, and L/E₀ is Galois of degree 4. If it were cyclic, its only element of order 2 would fix E; but complex conjugation at a real place of E₀ is an element of order 2 that does not, because E is totally complex. So the group is (Z/2)² and L = E(√d) with d in E₀. Conversely, the character of E(√d)/E on an ideal 𝔟 extended from E₀ is that of E₀(√d)/E₀ on the norm 𝔟², which is trivial. ∎

- **On a step.** Let 𝔩 be a prime of E of degree 1 above ℓ, at the root λ of P mod ℓ, and x_ℓ = λ + q²/λ the root of h below it. Then χ_d(𝔩) = (d(x_ℓ)/ℓ). The step of note 3 at 𝔩 changes φ by this sign.
- **How many.** If E(√d)/E is unramified, d has even valuation at every prime of E₀ that is not ramified in E. Such d, modulo squares, form a group of 2-rank 4 + r, with r the 2-rank of Cl(E₀) modulo the dyadic ramified prime: −1, two units, a generator at the dyadic prime, and r classes. The element x² − 4q² gives E itself. So G has at most 3 + r characters of order 2, before the conditions at 2.
- **None is attached to a ramified prime in the way CSIDH's are.** The d that occur are units of E₀ and generators of squares of ideals (§2.5). They depend on the arithmetic of E₀, that is on P.

### 2.2 One genus [D]

By note 6 a surface of S_P is a lattice L with F; for a maximal order L is an ideal 𝔞 of E with the form Tr(κ·x·ȳ), κ in E₀, and two such pairs are isomorphic when they differ by (y𝔞, κ/(yȳ)).

**Lemma 2.** Suppose exactly one finite prime 𝔮 of E₀ ramifies in E. Then any two pairs (𝔞, κ), (𝔞′, κ′) with unimodular form and with κ, κ′ of the same sign at each real place are isometric, compatibly with E, over Z_ℓ for every ℓ.

*Proof.* Put 𝔟 = 𝔞′/𝔞 and β = κ/κ′, so that 𝔟·𝔟̄ = (β) and β is totally positive. At a prime v of E₀ write 𝔟 = b_v·O; the two lattices are isometric at v when the unit u_v = β/(b_v·b̄_v) is the norm of a unit. This holds when v is unramified. For the Hilbert symbol of β and x² − 4q²: it is 1 at the real places because β is positive, and at the finite v ≠ 𝔮 because v splits or β has even valuation. By the product formula it is 1 at 𝔮, so β, and with it u_𝔮, is a norm from E_𝔔; an element whose norm is a unit is a unit. ∎

- In all 512 fields of note 13 exactly one finite prime ramifies (note 13 §3.4), and note 13 found the number of lattices equal to |G| wherever it was counted. Lemma 2 is the reason.
- **Consequence.** Any invariant that depends only on the lattice with its form and F over Z₂, or over any Z_ℓ, is constant on S_P. That covers the action of F on T/4T and on any T/2ᵏT, the discriminant form with F, and the Brauer group with its pairings. For the fourfold of a maximal order it also covers the action of Frobenius on the points of order n, which is multiplication by π on O_K/n whatever the class. It does not cover pairings on the fourfold that use its polarisation (§2.4).
- **Comparison with CSIDH.** For an imaginary quadratic order, the lattices of the curves of one endomorphism ring, with their forms, fall into several genera, and Cl(O)/Cl(O)² permutes the genera simply transitively. The assigned characters are functions of the genus; that is why a pairing at a ramified prime computes them [M]. The group of norm-one idèles modulo units, which is what acts here (note 13, Lemma 9), is for a quadratic field Cl(O) modulo its ambiguous classes, isomorphic to Cl(O)² by squaring. It acts on one genus. A character of order 2 of Cl(O)² is the restriction of a character of order 4 of Cl(O). I know of no method to compute those from two curves [M].
- **Non-maximal orders.** There S_P is a union of torsors with different local structure at the primes dividing the index (notes 7 to 12), and local invariants can separate them (§4.2: several types on one polynomial occur, and only for even index).

### 2.3 What a character would be on surfaces [D, sketch]

Assume the dictionary of note 6. The canonical lifts of the surfaces of S_P are conjugate over E, under a group isomorphic to G, and √d lies in their field of moduli. So a character with class field E(√d) is given by a function R of the moduli of the canonical lift with R(X^σ) = σ(R(X)) and R(X)² = c²·d for a constant c: modulo the prime above p, R takes the two values ±r₀ on S_P, with r₀² = c²·d(x₀) and x₀ the unit root of h. In words: a character is a function of the six lines whose square is the same for all surfaces of the torsor, and which is not itself the same.

- For CSIDH such functions exist for all orders at once, as self-pairings [M]. Here R would have to produce a square root of a unit of E₀ that changes with P.
- **Frobenius conjugation.** R(X^(p)) = R(X)^p. So φ(X^(p)) = φ(X) exactly when r₀ lies in F_p, that is when d is a square at 𝔭: χ_d(𝔭) = 1, as note 3 requires. If χ_d(𝔭) = −1 then r₀ and −r₀ are conjugate, and both have the same square class in F_q.

**Criterion.** A candidate that is unchanged by Frobenius conjugation can only be a character that is trivial on the class of 𝔭. Square classes in F_q of quantities defined over F_p, numbers of pairs of lines with a property, and numbers of points of auxiliary curves are all of this kind. A character with χ(𝔭) = −1 needs a value in a set on which Frobenius acts: over F₄₉ a fourth root of unity ±i, over F₂₅ an eighth or a third root of unity.

### 2.4 The fourfold side [D]

**The class group of K sees no quadratic character.** Note 18, findings 1 and 2: the ideal map Ψ of note 15 §3.1 is induced by x ↦ x_j·x_k·conj(x_l·x_m), the square of the map of tori ρ(x) = (x_j·x_k) where ρ has values in E; the isogeny with kernel A[𝔞] of note 15 §3.3 is the step ρ(𝔞) = ℓ·𝔩₁₂/𝔩₃₄, and Ψ(𝔞) is its square.

- **[D]** Let 𝔓 be a prime of K of degree 1 above ℓ, at the root a of χ. Ψ(𝔓) is the product of the primes of E at the roots a·b of P, b the other three roots of χ, and these lie over the three roots of h mod ℓ, one for each way to split the four roots in two pairs. So χ_d(Ψ(𝔓)) = (N_{E₀/Q}(d)/ℓ). For a unit d the norm is ±1, and ℓ ≡ 1 mod 4 because K contains i. In general N(d) is ± a power of 2 times a square, and χ_d is trivial on the image of Cl(K) when it is ± a square.
- **[E]** For every polynomial of §2.5 with a character, every character is trivial on Ψ(𝔓) for the first 40 primes 𝔓 of degree 1.
- So the cokernel of note 15, of order 2 to 8, maps onto G/G². It does not describe sets of surfaces that the fourfolds separate. Ideal classes of K are isogenies that preserve the polarisation, and they give squares of steps; the steps themselves are isogenies that multiply the polarisation by ℓ, and under them the fourfolds of all surfaces of S_P are related.

**What the surface determines.** The automorphisms ±1, ±i of a fourfold A with its Z[i]-structure act on the second exterior power through ±1. So the forms of A over F_q with the same surface are A and its twist by −1, with Frobenius π and −π; and the fixed part of the second exterior power is the same for the structures i and −i. The surface therefore gives χ only up to χ(X) ↦ χ(−X) and χ ↦ χ̄. For a maximal order c₁ = 4(±1 ± i), and all four values belong to the same surface. A rule that picks one of the two signs by the action on torsion would need exactly one of π, −π to be 1 modulo 2(1 + i). At p = 17 with c₁ = 4 + 4i and c₂ = −2, χ(1) = 8(27 + 8i) and χ(−1) = 8(45 − 8i), with 27 + 8i and 45 − 8i of odd norm; a multiple of 2(1 + i) would give a multiple of 64. So neither is.

**The character of d = −1.** If E(i)/E is unramified and the ramified prime splits in it, χ₋₁ is a character of G, and χ₋₁(𝔩) = (−1/ℓ) for a prime 𝔩 of degree 1 above ℓ: the degree of the step modulo 4.

- Over a prime field 𝔭 is trivial in G, so χ₋₁ can only exist for p ≡ 1 mod 4, and χ₂ only for p ≡ ±1 mod 8.
- On the fourfold the step at 𝔩₁₂ is an isogeny φ with φ*λ′ = ℓ·λ, so Weil pairings satisfy e′(φP, φQ) = e(P, Q)^ℓ. A fourth root of unity attached canonically to A by its pairing of order 4 and its endomorphisms would be raised to the power ℓ: unchanged for ℓ ≡ 1 mod 4, conjugated for ℓ ≡ 3. This is how the character δ of an imaginary quadratic order is computed for elliptic curves [M]; in the paper of Castryck, Sotáková and Vercauteren, §3.2 (read for note 23), the invariant is a Tate pairing T₄(P, Q) of a point of order 4 with a generator of the cyclic 2-primary torsion on the floor of the 2-volcano.
- The local class of the hermitian lattice of A at a prime v of K⁺ above 2 changes by the Hilbert symbol (ℓ, −1)_v, which is (−1/ℓ) when K⁺_v has odd degree over Q₂ and 1 when the degree is even. So the polarised lattices of the fourfolds of S_P should fall into two genera, separated by the multiplier modulo 4, exactly when K⁺ has a prime above 2 of odd local degree. The second exterior power does not see this: scaling the hermitian form by a unit u scales the form on it by u².
- **[E]** `sixlines_characters.py dyadic`: χ₋₁ is a character of G exactly when K⁺ has a prime above 2 of odd local degree: 477 of the 477 polynomials of §2.5 with a character. In the 80 cases where it is, 2 is unramified in K⁺ with primes of degree 1 and 3; in the other 397 it is totally ramified. This agreement is between a computation in E and one in K, and it is what the sketch predicts.
- What it would take: a point of order 4, or a class killed by (1 + i)³, on the Prym variety over a small extension of F_q; its Weil pairing with its image under an element of O_K that does not change sign under the twist, for instance π²; and the compatibility of the polarisations of the Prym varieties of X and of g·X. None of this was done here; note 23 reads the invariant off a Tate pairing of order 16 over F_(p⁴) instead of a point of order 4. The Prym variety is not principally polarised, and the equivariance of X ↦ A_X is not proved (note 15).

So the question about pairings has this answer. A pairing with values ±1, on the 2-torsion or the (1 + i)-torsion, is unchanged by a step, whose multiplier ℓ is odd; and the square classes that such pairings produce from six lines are the same on a whole torsor in the data (§4.2). The first level at which a pairing could see a character is 4, for the one character χ₋₁, where it exists.

### 2.5 The class-group side in numbers [E]

`sixlines_characters.py genus`, on the 480 polynomials with maximal order of note 13 with q ≤ 10,007, q = p or p² (GRH). Note 13's 32 polynomials with q = 30,011 and 100,003 are left out (§6).

| | |
|---|---|
| 2-rank of G: 0, 1, 2, 3, 4 | 3, 287, 160, 28, 2 polynomials |
| characters of order 2 in all (a basis for each G) | 699 |
| matched on all test primes by an S-unit d of E₀ | 699 |
| d a unit of E₀ | 437 |
| d not a unit: needs an ideal class of E₀ | 262 |
| d of one sign at the three real places, up to the factor x² − 4q² | 184 |
| d of mixed signs | 515 |
| 2-rank 2 or more with h(E₀) odd | 0 of 248 |
| polynomials with a character | 477 |
| d = −1 gives a character: all; q ≡ 1 mod 8; q ≡ 5 mod 8; q ≡ 3 mod 4 | 80 of 477; 53 of 126; 27 of 59; 0 of 292 |
| d = 2 or d = −2 gives a character | 0 |
| all characters trivial on Ψ(𝔓) for the first 40 primes 𝔓 of K (§2.4) | 477 of 477 (19,564 primes) |
| over F₂₅ and F₄₉ (12 polynomials): the class of 𝔭 trivial; not trivial but a square; some character −1 on it | 2; 1; 9 |

- "Not a unit" means that no product of units and of generators at the primes above 2 matches the character; a generator of the square of a non-principal ideal does.
- When h(E₀) is odd G has at most one character of order 2, as note 18 found on the 512 fields. The two cases of 2-rank 4 have h(E₀) = 8 and 32.
- d = −1 never gives a character for q ≡ 3 mod 4: over a prime field that is §2.4, and no field of p² elements in this sample is 3 mod 4.

---

## 3. Candidates

Let l₁, …, l₆ be the linear forms and [ijk] the determinant of three of them. For the surface X : εw² = l₁⋯l₆:

- **The ten classes.** For a split {i, j, k} ∪ {l, m, n}, t = [ijk][lmn] is multiplied by a₁⋯a₆ when the forms are rescaled, and ε by the same factor, so the class of ε·t in F_q*/F_q*² depends only on X and the labelling. A relabelling permutes the ten classes, with signs when −1 is not a square. The quadratic twist changes all ten. **Type** = the normal form of the vector of ten classes under the 720 relabellings.
- **Quantities of note 14 that are products of them [D, asserted in the script].** λ₅ = [135][245]/([145][235]), so λ₅λ₆ is the product of the classes of the four splits 135|246, 136|245, 145|236, 146|235, and (1 − λ₅)(1 − λ₆) that of 125|346, 126|345, 145|236, 146|235. The choice of twist of E (note 14, Lemma 2) is determined by these two.
- **Branch.** The number of pairs of lines for which the two branch values of the cross-ratio map are rational: the discriminant of the branch quadratic of note 14.
- **Quartic (q ≡ 1 mod 4).** The quartic characters χ₄(ε·t) of the ten products are defined up to a common sign. Their product, which is χ₄(ε²·∏[ijk]) over all twenty brackets, is an invariant with values ±1; the twist changes its sign. **Quartic type** = the normal form of the ten values under relabelling, the common sign and the choice of χ₄. **Orientation** = the sign of (Σ over the classes that are squares)·(Σ over the others, divided by i), which changes with the choice of χ₄.
- **Conic.** The square class of [123][145][246][356] − [124][135][236][456], which vanishes when the six lines touch a conic. Association and odd relabellings change its sign.
- **j2, j4.** The square classes of Σt² and Σt⁴ over the ten splits.
- **c₁.** From the quartics stored by notes 14 and 16: the values of c₁ over the pairs of one configuration in which a, b, c are rational, where the curve is normalised without a choice of Q.

Not candidates: the trace of E, which is not constant on a surface (note 15 §5.1), and anything that depends on a labelling.

---

## 4. Tests [E]

`sixlines_characters.py surfaces`: the 76 irreducible ordinary polynomials of notes 2 and 3 (44 over the primes 13 to 37, 4 over F₂₅, 28 over F₄₉), 839 surfaces in 1,651 configurations. 36 polynomials have maximal order (419 surfaces); G has even order for 33 of them and 34 have at least two surfaces.

### 4.1 Well defined on a surface

Agree / differ, over all 1,651 configurations: the value on a configuration and on its associate, and the value after one random relabelling of the six lines.

| candidate | associate, q ≡ 1 mod 4 | associate, q ≡ 3 mod 4 | relabelling, q ≡ 1 | relabelling, q ≡ 3 |
|---|---|---|---|---|
| type | 1,499 / 0 | 152 / 0 | 1,499 / 0 | 152 / 0 |
| branch | 1,499 / 0 | 152 / 0 | 1,499 / 0 | 152 / 0 |
| j2, j4 | 1,499 / 0 | 152 / 0 | 1,499 / 0 | 152 / 0 |
| conic | 1,499 / 0 | 66 / 86 | 1,499 / 0 | 71 / 81 |
| quartic, quartic type | 1,499 / 0 | not defined | 1,499 / 0 | not defined |
| orientation | 1,317 / 182 | not defined | 1,297 / 202 | not defined |

- The conic class is an invariant of the surface only for q ≡ 1 mod 4; the orientation is not one.
- **c₁.** 343 configurations have two or more pairs in which a, b, c are rational. In 263 of them c₁ takes two or more values; c₁ and −c₁ both occur in 137, c₁ and its conjugate in 131.

### 4.2 The type is a function of the polynomial

- On each of the 36 polynomials with maximal order all surfaces have the same type.
- Over these, the type is a function of (q mod 16, e₁ mod 16), in the 22 combinations that occur. For q ≡ 1 mod 4 two types occur, with six and with four classes that are squares: the first exactly when (e₁ − 2q)/4 is 2 or 3 modulo 4.
- The number of pairs with rational branch values is a function of the type over all 839 surfaces: 7 for the types with q ≡ 1 mod 4 that occur for maximal orders, 5 or 9 for q ≡ 3 mod 4, and 15 for two types that occur only for other orders.
- Whether c₄ is rational: in every pair that is not degenerate, exactly when q ≡ 1 mod 4 (asserted).
- Of the 40 polynomials whose order is not maximal, 21 have one type; these include all 9 with odd index. The 19 with several types all have even index.

This is what §2.2 predicts if the ten classes are local data at 2. I expect that they describe the action of (F − 1)/2 on T/2T, which for a maximal order is multiplication by an element of O_E on O_E/2; I did not prove it.

### 4.3 Constant or balanced on a torsor

For the candidates with values ±1 that are invariants of the surface. Torsors of maximal orders with at least two surfaces: 34, of which 23 have q ≡ 1 mod 4.

| candidate | torsors | constant | balanced | neither |
|---|---|---|---|---|
| conic (q ≡ 1 mod 4) | 23 | 3 | 5 | 15 |
| j2 | 34 | 4 | 11 | 19 |
| j4 | 34 | 5 | 5 | 24 |
| quartic (q ≡ 1 mod 4) | 23 | 4 | 5 | 14 |

- Examples: on the 26 surfaces of q = 37, e₁ = 10 (G cyclic of order 26) the conic class is +1 on 15 and −1 on 11, the quartic product on 11 and 15. On the 36 surfaces of q = 49, e₁ = 6 the four candidates split 22/14, 22/14, 14/22 and 12/24.
- The balanced cases are mostly small torsors: for a random function a torsor of 2 is balanced half the time.
- **All functions of the four.** 23 combinations of values (conic, j2, j4, quartic) occur on the 23 torsors with q ≡ 1 mod 4; a power sum is 0 on a few surfaces. Of the 4,194,303 functions of the combination, up to sign, one is constant or balanced on every torsor: the indicator of a combination that occurs on a single surface in the whole data, on the torsor of two surfaces of q = 25, e₁ = 14.
- **The quartic type** takes 215 values on the 363 surfaces with q ≡ 1 mod 4, counted torsor by torsor: 11 on a torsor of 12, 24 on one of 26. Over F₄₉ it takes each value on a surface and on its Frobenius conjugate. It is close to a complete invariant, and so says nothing by itself.

### 4.4 Frobenius conjugation

| candidate | F₂₅, agree / differ | F₄₉, agree / differ |
|---|---|---|
| type, branch, conic, j2, j4, quartic, quartic type | 37 / 0 | 1,080 / 0 |
| quartic values without identifying χ₄ and its conjugate | 37 / 0 | 176 / 904 |
| orientation | 37 / 0 | 712 / 368 |

On the lattice side, for the polynomials with data:

| field | polynomials with maximal order | conjugation trivial in G | not trivial, but a square in G | some character is −1 on it |
|---|---|---|---|---|
| F₂₅ | 2 | 1 | 0 | 1 |
| F₄₉ | 10 | 1 | 1 | 8 |

The two trivial cases are the polynomials with maximal order on which note 3 found every surface fixed, (14, 175, −17500) over F₂₅ and (46, −833, −220892) over F₄₉.

So over F₄₉ a complete set of characters cannot be computed by candidates of the first row. The orientation is the only candidate here that Frobenius moves, and it is 0 on 148 of the 248 surfaces of maximal orders over F₄₉ (and not an invariant, §4.1). For a candidate that Frobenius conjugation negates, "balanced" holds automatically on a torsor where conjugation has no fixed surface, so the data cannot confirm such a candidate; they can only reject it.

### 4.5 The torsors where χ₋₁ is a character

Among the 36: q = 29, e₁ = 2 and q = 37, e₁ = 2 (G cyclic of order 12, so χ₋₁ is the only character), q = 37, e₁ = 18 (G = Z/4 × Z/2), and five over F₄₉ (e₁ = 6, 10, 22, 26, 38). On the two torsors of 12 surfaces: conic 7/5 and 7/5; j2 6/5 (one 0) and 7/5; j4 4/8 and 7/4 (one 0); quartic 5/7 and 6/6. None is balanced on both. Over F₄₉ the character is −1 on Frobenius conjugation, and the candidates of the first row of §4.4 cannot be it.

---

## 5. What this means for the decision problem

- **What CSIDH leaks does not leak here.** There the assigned characters are one bit each of g, computable from X and g·X. Their analogue here is the genus of the lattice, and S_P is one genus. The decision problem for the six-line family is therefore not broken by the method of Castryck, Sotáková and Vercauteren applied to the K3 lattice.
- **The characters that exist are of the other kind.** They correspond to characters of order 4 in the quadratic case, and by that analogy a method for them would also be a method for those, which is not known [M]. There is one exception:
- **χ₋₁.** Where it exists it is the multiplier of the fourfold isogeny modulo 4, and the fourfold has a pairing of order 4. This is the one place where the known method might transfer, and on the three torsors tested it does (notes 23 and 24). It concerns 80 of 477 polynomials in note 13's sample, and one bit.
- **The simple quantities are useless in both directions.** Square classes of products of brackets are functions of P. The finer invariants I tried behave like unstructured functions on the torsor. With labels this can be tested, and it holds within the families tried ([note 28](28-leakage-audit.md)): No residue symbol of order 2 of the invariants, and no relational symbol of two surfaces, is a character or a function of the class beyond chance, on groups of order at most 8. Characters of order 4 were tested only through their squares: quartic symbols exist on 3 of the 9 labelled torsors, and no candidate could reach significance with an order-4 reading. Note 28 states the power of the test. It is still a negative result about lists.
- **Note 13 §6** stands in its conclusion ("none is shown to") and is corrected in its count and its description of the characters.

---

## 6. Check and limits

- **Internal checks.** The script asserts: the number of surfaces of a polynomial with maximal order is |G|; G agrees with note 13 for every polynomial of the `genus` run; every character of order 2 found is matched by an S-unit d of E₀ on all test primes (60 plus the number of generators); each ideal Ψ(𝔓) has norm ℓ³; the product formulas of §3 for λ₅λ₆ and (1 − λ₅)(1 − λ₆) in every non-degenerate pair; the two descriptions of the branch values agree; the product of the twenty brackets is a square; over a prime field the unit-root prime is trivial in G. All pass.
- **Identification of d.** A character is matched to d on finitely many primes, 60 more than there are unknowns. Lemma 1 is the proof; the run supports it and supplies the signs and the kind of d.
- **The largest fields are left out of §2.5.** For note 13's 32 polynomials with q = 30,011 and 100,003 the classes of some 110 ideals each are slow to compute; a run on them was stopped unfinished. `genus --min-q 30011` runs them.
- **Small fields for the surfaces.** p ≤ 37 and q = 25, 49. Torsors have at most 36 surfaces. "Neither" is decisive for a candidate; "balanced" on small torsors is not evidence for one.
- **One relabelling per configuration** was tried, drawn with a fixed seed. The type is a normal form over all 720 by construction.
- **Maximal orders only** for the character tests. For the others S_P is not one torsor.
- **Candidates are a list.** Polynomial invariants of higher degree, quartic or octic characters of other expressions, and anything built from the curve C beyond its point counts were not tried. In particular no pairing on the Prym variety was computed here (note 23 computes one).
- **§2.3 and §2.4 are sketches.** They assume the dictionary of note 6 and, for the fourfold, the equivariance that note 15 could not test. The statement that pairings of order 2 give only the ten classes is an expectation supported by §4.2, not a computation of pairings.
- **Class groups** are PARI's without certification (GRH).
- The method of Castryck, Sotáková and Vercauteren is described from memory; only the abstract was read for this note (their §3.2 was read for note 23, §5.1).

---

## 7. What this does not give

- A proof that no quadratic character of G is computable from a surface.
- The character χ₋₁ on surfaces, or a test of whether a pairing of order 4 on the Prym variety computes it; notes 23 and 24 give χ₋₁ on three torsors.
- Any candidate that Frobenius conjugation moves and that is an invariant of the surface, for the characters with value −1 on 𝔭.
- A proof that the ten square classes are the action of (F − 1)/2 on T/2T.
- The class-group side for note 13's 32 polynomials with q = 30,011 and 100,003.
- Anything for non-maximal orders, where local invariants at the primes dividing the index do separate the torsors of one polynomial; that is the question of levels of notes 7 to 12, not of characters.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_characters.py surfaces --out sixlines_characters.jsonl > sixlines_characters.log
python -I sixlines_characters.py genus --workers 8 --max-q 10007 --out sixlines_characters_genus.jsonl > sixlines_characters_genus.log
python -I sixlines_characters.py dyadic sixlines_characters_genus.jsonl > sixlines_characters_dyadic.log
```

About 20 seconds, 40 seconds with 8 processes, and a few seconds. The first needs the files of notes 2, 3, 14 and 16 (`sixlines_p{p}.jsonl`, `sixlines_q25.jsonl`, `sixlines_q49.jsonl`, `sixlines_cm.jsonl`, `sixlines_q49_cm.jsonl`, `sixlines_neighbours.jsonl`, `sixlines_prym*.jsonl`); the second needs `sixlines_groupsize.jsonl` of note 13, the third the output of the second. All are deterministic.

## Sources

- Castryck, Sotáková, Vercauteren, *Breaking the decisional Diffie–Hellman problem for class group actions using genus theory*, https://eprint.iacr.org/2020/151 : abstract read [P]; the method (assigned characters, Tate pairings of order m dividing the discriminant, the characters δ and ε at 2) recalled [M].
- Recalled, not re-read [M]: Castryck, Houben, Merz, Mula, van Buuren, Vercauteren on self-pairings and class group actions; Matsumoto, Sasaki, Yoshida on the six-line family; the invariants of six points in the plane (Coble).
- In this repository: notes [2](02-six-line-count.md), [3](03-torsor-action.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md), [18](18-checks-of-notes-13-to-16.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [28](28-leakage-audit.md).

