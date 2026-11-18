# 15. Does the surface-to-fourfold map respect the group action?

*[Note 14](14-fourfold-of-a-six-line-surface.md) built the abelian fourfold of a six-line surface over F_p and left open whether the map carries the action of the class group on surfaces to the action on fourfolds. This note asks how far that can be tested. Builds on notes 2, 3, 6, 13 and 14. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_dictionary.py`. Status: re-checked in [note 18](18-checks-of-notes-13-to-16.md), which corrected §3.1 and §3.3 and found one more gap in §2; the corrected versions are stated here.*

---

## 1. Result

1. **The direct test was not done, and cannot be done with what exists.** It needs one step at a prime ℓ ≠ p computed for explicit objects on one side. On surfaces there is none ([note 3](03-torsor-action.md) §4). On fourfolds it means computing an isogeny with kernel (Z/ℓ)⁴ from the Prym variety of a curve of genus 5, or deciding whether two given Prym varieties are related by one. I know no invariant that decides this, and I did not attempt the isogeny (§5). A test was carried out later, by the CM method, in [note 20](20-cm-method-and-the-action.md) (§5.1).
2. **Why the map should respect the action [D, sketch].** The two actions are defined on lattices, and the map on lattices is M ↦ fixed part of the second exterior power of M over Z[i]. An isogeny of fourfolds changes M, and the K3 lattice changes with it. So the statement reduces to two inputs: the correspondence over C holds for lattices and not only up to isogeny (Lombardo–Peters–Schütt, read), and it is compatible with canonical lifts (not proved). §2 gives the argument and its gaps.
3. **The dictionary between the two groups is explicit and consistent [D, E].** For 278 polynomials there is one homomorphism from the class group of the field K of the fourfold (degree 8) to that of the field E of the surface (degree 6) that sends each of 77,946 primes of K to the ideal the dictionary predicts. It maps Cl(K⁺) into the subgroup that note 13 divides out. The homomorphism f is induced by x ↦ x_j x_k·conj(x_l x_m), which is the square of the map of tori of §2 where that map is defined ([note 18](18-checks-of-notes-13-to-16.md), finding 1). The induced map Cl⁻(K) → G has kernel of order 1 to 8 and cokernel of order 1 to 8, both powers of 2; the cokernel is not always killed by 2: it has an element of order 4 in 46 of the 278 cases. The odd parts of the two groups are isomorphic in all 278 cases. So "the same group up to 2-groups" of note 13 §4.3 holds, with the map written down (§3). These are the groups of the maximal order of K; the lattices of the surfaces are not lattices of the maximal order ([note 20](20-cm-method-and-the-action.md)), so this kernel and cokernel do not describe surfaces.
4. **What a step is on each side [D].** For a prime ℓ ≡ 1 mod 4 at which χ has roots a₁, …, a₄, the isogeny of the fourfold with kernel A[𝔓₁𝔓₂·conj(𝔓₃𝔓₄)], of order ℓ⁴ and stable under i, corresponds to the step of the surface at the prime of E at the root a₁a₂ of P, divided by its conjugate. The image of Cl⁻(K) in G misses a subgroup of index 2 to 8 because such an isogeny multiplies the polarisation by ℓ and is not an element of the norm-one torus; its square is (§3.3). f sends the class of 𝔞 in Cl⁻(K) to the fourth power of the class of 𝔩₁₂ (note 18, finding 2).
5. **Association is invisible to the curve [E].** Replacing both lines of a pair by their conjugates under the involution of the cross-ratio map gives the associated configuration: 52,207 pairs out of 52,207 that are not degenerate, over all 3,990 configurations. The curve of note 14 is unchanged by this replacement. So a configuration and its associate, which are the same K3 surface, have the same curves (§4).

So the question "does the map respect the action" has the answer "it should, for a reason that does not depend on computation, and the group-level consequences that can be checked all hold". A test on surfaces is not made here; note 20 makes one by the CM method. The conclusion of notes 13 and 14 for hardness is unchanged.

---

## 2. Why the map should respect the action [D, sketch]

**Lattices.** By [note 6](06-torsors-under-cm-class-groups.md), Propositions A to C, an ordinary surface X over F_q with polynomial P is a lattice L in a fixed space V of rank 6 with F acting, and the group is T¹_E(A_f) modulo the stabiliser and the global elements, with T¹_E the torus of elements of norm 1 in E = Q[F]. By Deligne's description [M], an ordinary fourfold with Frobenius polynomial χ·χ̄ and multiplication by Z[i] is a lattice M in a fixed space of rank 4 over Q(i), with π acting, and isogenies commuting with i are the lattices in the same space.

**The map on lattices.** Over C the correspondence sends a principally polarised fourfold A with H¹(A, Z) = M to the surface whose transcendental lattice is the fixed part of Λ²_{Z[i]}M under the antilinear involution, with the form multiplied by 2 [P: Lombardo–Peters–Schütt §2.2 and Theorem 3.4]. This is a statement about lattices. If A′ is A divided by a subgroup that is maximal isotropic for the polarisation and stable under i, with lattice M′ ⊃ M, then A′ is again principally polarised and its surface has the lattice built from M′.

**The map on tori.** x ↦ (x_j x_k)_{j<k} sends the torus GT_K of elements x of K with x·x̄ = m rational and N_{K/Q(i)}(x) = m² (note 18, finding 1) to the corresponding torus of E, and scalars to scalars. Proposition Y of note 13 says exactly that F on the K3 side is the image of π. So if the lattice of X is built from the lattice M of A, the lattice of g·X for g the image of x is built from x·M.

**From C to F_q.** Both descriptions over F_q come from canonical lifts: the lattice of X is that of its canonical lift (Nygaard, Taelman [M]) and the lattice of A that of its canonical lift (Deligne [M]). What is needed is that the fourfold of the canonical lift of X is the canonical lift of the fourfold of X. The fourfold of the canonical lift has complex multiplication, because the Mumford–Tate group of the surface is a torus and so is its preimage; and its reduction is the fourfold of X if the construction of note 14 commutes with reduction, which it does where it is defined by the formulas of note 14 §2. A lift on which Frobenius acts is the canonical one. A lift with complex multiplication need not be the canonical one, since quasi-canonical lifts have it too; what is needed is that Frobenius itself lifts (the last gap below, found in [note 18](18-checks-of-notes-13-to-16.md), finding 3).

**Gaps.**

- That the canonical lift of a six-line surface is a six-line surface [H].
- That the formulas of note 14, which pass through a Prym variety that is isogenous to the principally polarised fourfold and not equal to it, give the same lattice statement up to a fixed isogeny of degree a power of 2. I did not determine that isogeny.
- The prime 2, where the K3 form is the hermitian form times 2 and the fields are ramified.
- The identification of the two actions with one convention for the direction of a step (the sign left open in Proposition D.1 of note 6).
- That Frobenius lifts to the fourfold of the canonical lift. That follows if the six-line model of the canonical lift, and with it the curve of note 14, is defined over the unramified ring W(F_q).

None of this is a computation, and none of it was tested on surfaces.

---

## 3. The dictionary between the two groups

### 3.1 Statement [D]

Let P be as in note 13, χ its quartic over Z[i], E = Q[t]/P and K = Q[X]/(χ·χ̄), with E₀ and K⁺ the real subfields of index 2. Assume χ·χ̄ irreducible.

A prime 𝔓 of K of degree 1 over a prime ℓ ≡ 1 mod 4 is (ℓ, i − i_ℓ, π − a) with i_ℓ² = −1 and a a root of χ mod (ℓ, i − i_ℓ). The map of tori sends an element with valuation 1 at 𝔓 and 0 at the other primes over ℓ to an element with valuation 1 at each prime of E at a root a·b of P, b the other three roots of χ. So define

> Ψ(𝔓) = (ℓ, P_a(t)), P_a(t) = a³·(χ/(X − a))(t/a) mod ℓ,

an ideal of norm ℓ³ of E.

**Claim.** [𝔓] ↦ [Ψ(𝔓)] extends to a homomorphism f : Cl(K) → Cl(E); it maps the classes extended from K⁺ into H = ⟨classes extended from E₀, primes of E ramified over E₀⟩; and the induced map from Cl⁻(K) = Cl(K)/Cl(K⁺) to G = Cl(E)/H, the group of note 13, Lemma 9, has kernel and cokernel of order a power of 2.

*Why.* E embeds in the Galois closure of K as the field generated by a product of two conjugates of π over Q(i). x ↦ (x_j x_k) maps K^× to E(i)^×, and into E^× only on the subtorus where x·x̄ = m is rational and N_{K/Q(i)}(x) = m². The ideal Ψ(𝔓) defined above is the image of 𝔓 under the norm from E(i) to E of that map, x ↦ x_j x_k·conj(x_l x_m) with {l, m} the complement of {j, k}; on the subtorus it is the square of x ↦ x_j x_k. It sends principal ideals to principal ideals. With the map in the other direction, y ↦ (∏_{k≠j} y_{jk})_j, the two composites are x ↦ (x/x̄)²·N·N̄² with N = N_{K/Q(i)}(x), and y ↦ (y/ȳ)²·N_{E/Q}(y)². Q(i) and Q have class number 1, and y/ȳ has the class of y² modulo H, so both composites are fourth powers on Cl⁻(K) and on G (not squares; see [note 18](18-checks-of-notes-13-to-16.md), finding 1). So the cokernel is killed by 4 and the kernel by 8; both are 2-groups, and the odd parts are isomorphic for every P. This is a sketch; the run below is the check. Measured on the 278 polynomials (`sixlines_review.py kernel`, note 18): the kernel is (Z/2)^k with k ≤ 3 in all of them; the cokernel has an element of order 4 in 46.

### 3.2 Check [E]

`sixlines_dictionary.py dictionary` takes the 44 polynomials of the enumerated surfaces and the 246 further polynomials with maximal order and q ≤ 101 of note 13. For 12 of the 44, χ·χ̄ is reducible and the polynomial is skipped. For each of the other 278 it computes the classes of 𝔓 and of Ψ(𝔓) for all primes ℓ ≡ 1 mod 4 up to 2,000 and all roots a, between 210 and 340 primes 𝔓 per polynomial, and checks:

| check | holds for |
|---|---|
| the ideal (ℓ, i − i_ℓ, π − a) has norm ℓ and Ψ(𝔓) has norm ℓ³ (the script stops otherwise) | all |
| the classes of the 𝔓 generate Cl(K) | 278 of 278 |
| one homomorphism f with f([𝔓]) = [Ψ(𝔓)] for every 𝔓, compatible with the relations of Cl(K) | 278 of 278 |
| f(Cl(K⁺)) ⊂ H | 278 of 278 |
| the odd parts of Cl⁻(K) and G have the same elementary divisors | 278 of 278 (224 with a non-trivial odd part) |

The 2-part:

| order of the kernel of Cl⁻(K) → G | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| polynomials | 6 | 205 | 64 | 3 |

| order of the cokernel | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| polynomials | 10 | 148 | 105 | 15 |

The two groups have the same order in 209 cases and the same elementary divisors in these 209. The largest G is cyclic of order 246, with Cl⁻(K) cyclic of order 246 (q = 101, (e₁, e₂, e₃) = (10, −6969, −204020)).

These are the class groups of the maximal orders. The lattices of the surfaces are not lattices of the maximal order of K, which are all odd as hermitian Z[i]-lattices, so the kernel and cokernel of this section do not describe surfaces; on the lattices that do, G acts freely, and transitively in 14 of the 15 fields tested ([note 20](20-cm-method-and-the-action.md)).

### 3.3 What a step is on the fourfold [D]

Let ℓ ≡ 1 mod 4 with χ split mod (ℓ, i − i_ℓ), roots a₁, …, a₄, primes 𝔓₁, …, 𝔓₄, and let 𝔩_{jk} be the prime of E at a_j a_k. Then a₁a₂a₃a₄ = q², so 𝔩₃₄ is the conjugate of 𝔩₁₂, and, with ρ : x ↦ (x_j x_k) the map of tori of §2, which is defined on 𝔞,

> ρ(𝔞) = ℓ·𝔩₁₂/𝔩₃₄ and Ψ(𝔞) = ρ(𝔞)² = (ℓ·𝔩₁₂/𝔩₃₄)² for 𝔞 = 𝔓₁𝔓₂·conj(𝔓₃)·conj(𝔓₄).

The norm of Ψ(𝔞) is ℓ¹²; the formula for Ψ(𝔞) holds for 48 primes out of 48 on 8 polynomials (`sixlines_review.py step`, [note 18](18-checks-of-notes-13-to-16.md), finding 2, which corrected an earlier Ψ(𝔞) = ℓ·𝔩₁₂/𝔩₃₄). Of the 16 ideals made of one of 𝔓_j, conj(𝔓_j) for each j, only the 6 with two of each kind are in the domain of ρ; the other ten are steps too ([note 20](20-cm-method-and-the-action.md)).

𝔞·conj(𝔞) = (ℓ) and 𝔞 has norm ℓ⁴: the isogeny A → A/A[𝔞] has kernel (Z/ℓ)⁴, stable under i and maximal isotropic, as note 13 §4.3 expected. On the surface it is the step 𝔩₁₂/𝔩₃₄, which under the identification of Lemma 9 is the class of 𝔩₁₂ in G.

This also explains the cokernel. The element 𝔞 has multiplier ℓ and is not in the norm-one torus of K; 𝔞²/ℓ is, and corresponds to the class of 𝔞 in Cl⁻(K), which ρ sends to the class of 𝔩₁₂ squared, and f, which is ρ² there, to its fourth power (note 18, finding 2). So Cl⁻(K) sees the squares of the steps, and the group that acts on principally polarised fourfolds by isogenies of this kind is larger than Cl⁻(K) by elements of order 2. I did not compute that larger group.

Every character of order 2 of G is trivial on the image of f ([note 19](19-quadratic-characters-from-a-surface.md) §2.4): for the character of E(√d)/E, d in E₀, and a prime 𝔓 of K of degree 1 above ℓ, the value on Ψ(𝔓) is (N_{E₀/Q}(d)/ℓ), and it was 1 for 19,564 primes over 477 polynomials. So the cokernel of f maps onto G/G². It is not a set of surfaces that the fourfolds separate: the steps between the cosets are the isogenies with multiplier ℓ. The character of d = −1, where it is a character of G, is that multiplier modulo 4.

**For the reduction of note 13 §4.3.** Given two surfaces, an isogeny between their fourfolds with kernel A[𝔟] gives the element of G that relates the surfaces as the image ρ(𝔟) of 𝔟, the square root of the ideal that the formula for Ψ on primes gives (note 18, finding 2), provided §2 holds. The ambiguity is in a group of order at most 8 in these data.

---

## 4. Association [E, D]

For a pair (L5, L6), the cross-ratio map λ of note 14 has degree 2 on the pencil of lines through L5 ∩ L6. Let L5′, L6′ be the lines conjugate to L5, L6 under its involution.

- **[D]** The curve of note 14 for the pair (L5, L6) is built from λ, λ(L5), λ(L6) and the ratio N₁/D₁, which do not change when L5, L6 are replaced by L5′, L6′ together. So the configurations (L1, …, L4, L5, L6) and (L1, …, L4, L5′, L6′) have the same curve.
- **[E]** `sixlines_dictionary.py association`: for every configuration of the ten primes and every pair, the second configuration is the associated one (Gale transform), in 52,207 of the 52,207 pairs for which L5′, L6′ exist, are distinct from L5, L6 and from each other, and give lines in general position; 7,643 pairs are skipped. In 26,867 of them the configuration is its own associate.

A configuration and its associate are the same K3 surface (note 2). This is the statement of Lombardo–Peters–Schütt that the fourfolds parametrise six lines up to relabelling and association, seen over F_p. It gives no step: note 3 §4 already noted that the natural partner of degree 8 leads back to the associated configuration.

---

## 5. What a real test would take

Any of these would be one, in decreasing order of what I think is feasible.

1. **An isogeny of a Prym variety.** For a polynomial with several surfaces, say the six of p = 23, take ℓ = 5 or 13, form 𝔞 as in §3.3, compute A/A[𝔞] and recognise it as the fourfold of another surface. Isogenies with kernel (Z/ℓ)^g are computed with theta functions [M: Lubicz, Robert]; in dimension 4 with a polarisation that is not principal this is a project, not a script. The lattice side of note 3 predicts the permutation, so the outcome could be compared with a cycle type.
2. **An invariant of the fourfold that the lattice side can compute.** For polynomials with non-maximal orders the level of a surface is known on the lattice side (notes 7 to 12). On the fourfold it shows in the group structure of the rational ℓ-torsion, which needs arithmetic in the Jacobian of the curve of genus 5. This would test the map on levels, not on steps inside a level.
3. **Frobenius conjugation over F_{p²}.** Note 3 compared it with the lattice side for surfaces. The map of note 14 commutes with it by construction, so on the fourfold side it only tests Deligne's theory.

I did none of them.

### 5.1 Attempt at item 1

I looked for a route to the isogeny that fits the tools used in this repository (PARI, exact arithmetic over small fields) and found none. Nothing was computed beyond the two small observations below.

**What the problem is.** Take a surface X with its curve C_X over F_p, form the kernel A[𝔞] of §3.3 in the Prym variety, and decide which of the finitely many surfaces Y with the same polynomial has A_Y ≅ A_X/A[𝔞]. The kernel is the easy half. The hard half is to recognise the quotient.

**Routes examined, and where each stops. [D]**

- *Local data cannot recognise the quotient.* All fourfolds of a torsor have the same zeta function, the same groups of points as modules over the endomorphism ring, and isomorphic torsion with its pairing. Which Y is A_X/A[𝔞] is the class of the module Hom(A_X, A_Y), a global invariant.
- *Theta functions.* The known algorithm for an isogeny with kernel (Z/ℓ)^g needs theta coordinates of the source. For a Prym variety of a curve of genus 5 that is not hyperelliptic there is no formula for them from the curve, and the polarisation is not principal. This is a research problem in itself.
- *A correspondence between the curves.* An isogeny A_X → A_Y is a divisor class on C_X × C_Y. Its numerical class is known in advance, but finding an effective divisor in it, or showing there is none, needs the Picard group of a surface of general type.
- *An invariant with a modular equation.* The elliptic curve E of note 14 is not one. At p = 23, for the polynomial (2, 207, −2116), the curves E of the 15 pairs of a single surface already have traces 0, ±4 and ±8, so E is not constant up to isogeny even on one surface, and steps cannot be read from isogenies between the E's. [E, a one-off computation, not in `experiments/`]
- *Characteristic 0.* With numerical period matrices of the curves C one could detect homomorphisms between two Prym varieties by lattice reduction. But the surfaces to compare are canonical lifts, whose coordinates are not known; getting them needs the inverse of the period map of the surfaces (Matsumoto's theta constants, or periods of the K3 surfaces themselves), independently of the fourfold. Neither is available here.
- *Characteristic ℓ.* By the theory of complex multiplication the ℓ-step in characteristic p is Frobenius conjugation in characteristic ℓ on the reductions of the same lifts. Matching the surfaces over F_p with those in characteristic ℓ again needs the lifts.

**The one isogeny of the right kind that is explicit. [D]** Frobenius. Over F_q its kernel is A[(π)], and the map of tori ρ of §3.3 sends π to F, with Ψ((π)) = (F)² (note 18, finding 2), because the conjugates of F are the products of two conjugates of π (Proposition Y). At p, with a₁, a₂ the two roots of χ that are not units, (π) = 𝔓₁𝔓₂·conj(𝔓₃𝔓₄) and (F) = p·𝔭/𝔭̄ with 𝔭 the prime of E at a₁a₂: the formula of §3.3 with ℓ = p. So the dictionary sends the Frobenius isogeny of the fourfold to the step 𝔭/𝔭̄ of the surface, which note 3 identified with Frobenius conjugation and tested against the lattices over F₂₅ and F₄₉. The map of note 14 commutes with Frobenius conjugation because it is given by formulas over F_p. This is consistent, and it is automatic: it does not test the map at a prime other than p.

**What would be needed.** Either theta coordinates of these Prym varieties with their multiplication by i, or the moduli point of a canonical lift computed without the fourfold. With the second, the lattice side would label every surface, and homomorphisms between the fourfolds of neighbours could be checked numerically.

So the test of item 1 was not carried out here, and I did not see how to carry it out with the means used up to this note.

A test was carried out in [note 20](20-cm-method-and-the-action.md), by the last two routes of the list above taken together. The moduli points of the lifts are computed from lattices of the fourfold by Matsumoto's theta constants (the CM method); for 15 polynomials with p ≤ 37 they reduce to exactly the enumerated surfaces (74 of 74). The step at a prime ℓ ≡ 1 mod 4 that splits in K is then Frobenius conjugation in characteristic ℓ: the ideal 𝔓_j𝔓_k·conj(𝔓_l𝔓_m) acts on the roots of the class polynomial as x ↦ x^ℓ modulo the prime of E at a_j a_k in 346 of 346 cases tested, and its inverse in none of 136. It also changes two statements of this note, stated in their corrected form in §1 item 3, §3.2 and §3.3: the kernel and cokernel of §3.2 do not describe surfaces, and the ten ideals of §3.3 outside the domain of ρ are steps too. The Prym variety of note 14 is still not compared with these lattices, and one sign was fixed by the test itself.

---

## 6. Check and limits

- The dictionary run uses PARI's class groups without certification (GRH). Its checks are assertions on classes of ideals; a failure would have shown as a prime whose image is not the predicted class. There was none in 77,946.
- The claim of §3.1 is argued in a few lines and not written out in detail; the run is the evidence for the form in which it is stated here, in particular for the choice of H and of Cl⁻(K). The run records the orders of kernel and cokernel, not their structure; the structure was measured with `sixlines_review.py kernel` in note 18, finding 1 (§3.1).
- Primes ℓ ≡ 3 mod 4 and primes of K of higher degree are not used; they are not needed to generate the class groups in these cases.
- The 12 polynomials with χ·χ̄ reducible are skipped; there K is not a field.
- Only maximal orders of E and K. For the 20 polynomials of the enumerated surfaces whose order is not maximal, the groups compared are those of the maximal orders, not the groups that act.
- §2 is an argument with stated gaps, and §3.3 depends on it for its meaning. Nothing in this note touches a surface except §4.

---

## 7. What this does not give

- A test of the action on surfaces at a prime ℓ ≠ p. This was the question, and it stays open here, also after the attempt of §5.1; note 20 carries out a test by the CM method.
- The group that acts on principally polarised fourfolds with multiplication by Z[i], and its exact relation to G at 2.
- The reverse map from a fourfold to six lines, and any step.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_dictionary.py dictionary --workers 14 --bound 2000 --out sixlines_dictionary.jsonl > sixlines_dictionary.log
python -I sixlines_dictionary.py association 5 7 11 13 17 19 23 29 31 37 > sixlines_dictionary_association.log
```

A few minutes each. The first needs `sixlines_cm.jsonl` and `sixlines_groupsize.jsonl`, the second the files `sixlines_p{p}.jsonl` and `sixlines_prym.py`.

## Sources

- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, https://arxiv.org/abs/1209.5997 : §2.2 and Theorem 3.4, as in note 13 [P].
- Recalled, not re-read [M]: Deligne on ordinary abelian varieties; Nygaard and Taelman on canonical lifts of ordinary K3 surfaces; Lubicz and Robert on isogenies by theta functions.
- In this repository: notes [2](02-six-line-count.md), [3](03-torsor-action.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md).

