# 20. The CM method and the action

*Notes [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md) and [15](15-fourfold-map-and-the-group-action.md) reduce the hardness question for the class-group action on six-line surfaces to abelian fourfolds with multiplication by Z[i], but there no surface over F_p is labelled by an ideal class, so the action is not evaluated (note 15 §5.1, [note 18](18-checks-of-notes-13-to-16.md) §7). This note supplies the labels by the CM method: lattice, period point, theta constants, six lines over a number field, reduction modulo p. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_cmmethod.py`. Checked by re-derivation in [note 25](25-checks-of-notes-20-23-24.md), except for the points listed in §8. Where a derivation enters the computation, its only support is that the computation comes out; §8 says what that covers.*

---

## 1. Result

1. **The CM method works for this family [E].** For 15 polynomials of the enumerated surfaces (p from 17 to 37, maximal order on the K3 side, χ·χ̄ irreducible, at most 8 surfaces), the invariants of six lines computed from lattices over C, recognised as elements of the field E and reduced modulo the prime of E above p at which F is a unit, are the invariants of the surfaces found by enumeration: 74 surfaces of 74 (148 configurations), none missing and none extra (§6).
2. **The lattices are not those of the maximal order [E, D].** Note 15 worked with fourfolds with multiplication by the maximal order of K. For all 21 fields, every principally polarised lattice of the maximal order is odd as a hermitian Z[i]-lattice (234 lattices), and a six-line surface needs an even one. The lattices of the surfaces are modules over R = Z[i, π, q/π, (π − 1)/(1 + i)] and its conjugate, an order of index 8 in O_K, and they have index 2 in their O_K-span (64 and 16 in the one field treated where 2 is not totally ramified in K). The generator (π − 1)/(1 + i) is the condition that the six lines are rational. With it, the number of even lattices of each type equals the number of surfaces in every field treated (§5).
3. **Formulas [P].** The theta constants, the ten splittings, the invariants J₂, …, J₆ and the groups are taken from Matsumoto–Terasoma and Clingher–Malmendier–Shaska, who quote Matsumoto (1993). The ten squares see the fourfold with a level structure modulo 1 + i, which is an ordering of the six lines; S₆ acts on it. The labels use only the invariants of the full group, so no level structure is computed (§2).
4. **The action on lattices is the group G of note 13, through the dictionary of note 15 as corrected by note 18 [E].** For a prime ℓ ≡ 1 mod 4 that splits completely in K, the permutation of the lattices induced by 𝔓_j𝔓_k·conj(𝔓_l𝔓_m) depends only on the class in G of the prime of E at a_j a_k: 378 ideals, no conflict; products compose; a non-trivial class has no fixed point. In 13 of the 15 fields every class of G occurs among the ideals tried, in one more the classes that occur generate G, and there G acts simply transitively on the lattices of a type: no two lattices give one surface, and no surface is out of reach. In the remaining field (q = 37, e₁ = 18) the ideals above primes ℓ ≡ 1 mod 4 give a subgroup of index 2, the kernel of the character "ℓ mod 4" that note 19, item 7, found for this polynomial; the other half needs primes ℓ ≡ 3 mod 4, which the script does not treat (§7.2).
5. **The ℓ-step is Frobenius conjugation in characteristic ℓ [E].** Raising to the ℓ-th power on the roots of the class polynomial modulo the prime of E at a_j a_k is the permutation of item 4: 346 cases of 346 tested (32 more could not be tested at the precision reached). With the inverse permutation it fails in 136 of 136, and with the ideal of another pair that acts differently in 1,156 of 1,156 (§7.3). This is the test that note 15 §5.1 could not make. One sign, the orientation of the period point, was fixed on one field (q = 23, e₁ = 10: 12 cases at four primes, of which the 8 at the first three primes are in the table of §7.3) and then held in the 128 cases of the other fields that depend on it.
6. **The ten other ideals [D, E].** The ideals with one or three of the 𝔓_j, which note 18 (finding 2) left open, also permute the lattices of the surfaces, and the square of each is the class that the map of tori predicts: 594 of 594 tested (§7.4).
7. **Surfaces over F_p now carry labels [E].** 67 of the 74 surfaces are labelled by an element of G, up to a translation and the orientation, and each class acts on the labels as its translation (§7.5). In one field the polynomials that carry the action to F_p were not recognised at the precision reached.
8. **A first use of the labels [E].** The square classes in F_q of J₂, J₂J₆ and J₃J₅ are invariants of a surface. Against the labels, none of them is a character of G: on the 8 torsors with at least 4 surfaces where this could be tested, each is a character or constant on at most one (§7.6).

What this is not. It is not a test of the Prym construction of note 14: the lattices come from the correspondence over C, not from the curve. It does not touch hardness: the method needs a class polynomial of degree |G|, so it evaluates the action only where the group is small. And the formulas of Matsumoto's paper were read in two later papers, not in the paper itself.

---

## 2. The formulas [P]

**The domain and the groups.** D = {τ ∈ M₂(C) : (τ − τ\*)/2i positive definite}. Γ = U₂₂(Z[i]) = {g ∈ GL₄(Z[i]) : g·I·g\* = I}, I = [[0, −1₂], [1₂, 0]], acts by τ ↦ (g₁₁τ + g₁₂)(g₂₁τ + g₂₂)⁻¹. Γ(1 + i) is the subgroup of g ≡ 1 mod 1 + i, and the transposition τ ↦ ᵗτ extends both groups (Matsumoto–Terasoma §3.2, §3.3; Clingher–Malmendier–Shaska (2.11) to (2.15)).

**Theta constants.** For a, b in Q(i)²,

> Θ_{a,b}(τ) = Σ_{n ∈ Z[i]²} e[ ½ (n + a) τ (n + a)\* + Re((n + a) b\*) ],  e[x] = exp(2πix),

and Θ[ab] = Θ_{a/(1+i), b/(1+i)} for a, b in {0, 1}². The ten characteristics with a·b\* ≡ 0 mod 1 + i are matched with the ten splittings of {1, …, 6} into two triples (Matsumoto–Terasoma (19) and Table 1):

| ⟨J⟩ | 123 | 124 | 125 | 134 | 135 | 145 | 234 | 235 | 245 | 345 |
|---|---|---|---|---|---|---|---|---|---|---|
| [a₁a₂b₁b₂] | 1111 | 0011 | 0010 | 0001 | 0000 | 1100 | 1001 | 1000 | 0100 | 0110 |

Their Remark 2 says this table differs from Matsumoto's because the basis of the transcendental lattice differs.

**The inverse of the period map.** Let x be the 3 × 6 matrix of the lines and x⟨J⟩ = D(J)·D(Jᶜ) the product of two minors. Then Θ[ab]² is a form of weight 2 with character det for Γ(1 + i) extended by the transposition, Θ[ab]²(ᵗτ) = Θ[ab]²(τ), and

> [ Θ⟨J⟩(τ)² ]_J = [ x⟨J⟩ ]_J in P⁹

when τ is the period point of the surface of x (Matsumoto–Terasoma, Proposition 7, quoted from Matsumoto; their main theorem (4) gives the constant, which is not needed here). Association is the transposition (their Proposition 6 (4)).

**Level and invariants.** Γ/Γ(1 + i) ≅ Sp₄(F₂) ≅ S₆ acts on the ten squares as the relabelling of the lines (Clingher–Malmendier–Shaska, Theorem 2.8 and the text after (2.41); Lombardo–Peters–Schütt, Remarks 1.5). The ten products t₁, …, t₁₀ of (2.1) satisfy the fifteen linear relations (2.2) and span a space of dimension 5; with x₁, …, x₆ of (2.34) and s_j = Σx_k^j,

> J₂ = s₂/12, J₃ = s₃/12, J₄ = (4s₄ − s₂²)/64, J₅ = (5s₂s₃ − 12s₅)/240, J₆ = (3s₂³ − 4s₃² − 18s₂s₄ + 24s₆)/576

are invariant under S₆, J_k(τ) is a form of weight 2k for the full group, and [J₂ : … : J₆] in P(2, 3, 4, 5, 6) is the moduli point of six lines up to relabelling and association (their (2.32), Lemma 2.17, Corollary 2.25). In (2.32) the text has 20s₅ where (2.31) and (2.33) give 12s₅; the script uses 12, and the invariance under S₆ is checked (§4).

So **what the labels see**: the ten squares are functions of the fourfold with a symplectic basis of L/(1 + i)L, up to the conjugate structure; that is an ordering of the six lines up to association, and S₆ permutes the orderings. A label in this note is the point [J₂ : … : J₆], which is the K3 surface: no ordering and no choice between a configuration and its associate. The enumeration of note 2 lists configurations, so a surface that is not self-associated appears there twice.

**The fourfold side.** Principally polarised fourfolds with an automorphism of order 4 of type (2, 2) correspond to 2-planes in C⁴ = H₁ ⊗ R on which the hermitian form is positive, and the surface has transcendental lattice the fixed part of the second exterior power (Lombardo–Peters–Schütt §2.1, §2.2, as in note 13 §4.1). Their lattice is Z[i]⁴ with the form of matrix I.

**Not read.** Matsumoto (1993) and Matsumoto–Sasaki–Yoshida (1992) themselves. The theorem of Shimura and Taniyama used in §7.1 is recalled [M].

---

## 3. From a lattice to a period point [D]

K = Q(i)[X]/χ, with χ the quartic of note 13, Proposition Y, and φ₁, …, φ₄ its embeddings over Q(i). A **type** is a pair T = {j, k}; its CM type is {φ_j, φ_k, φ̄_l, φ̄_m}, so that i has eigenvalues i, i, −i, −i. E acts on the holomorphic 2-form through F ↦ π_jπ_k; a type and its complement give conjugate embeddings of E and the same surfaces, so three types are used: {1,2}, {1,3}, {1,4}.

A **polarised lattice** is a pair (L, ξ): L ⊂ K a lattice stable under Z[i], ξ in K with ξ̄ = −ξ and ξ·L = conj(L^∨) for the trace dual L^∨, so that Tr_{K/Q}(ξxȳ) is a principal polarisation. It has type T when Im φ(ξ) < 0 for φ = φ_j, φ_k and > 0 for φ_l, φ_m. (L, ξ) and (vL, ξ/(v v̄)) are the same object.

- **The hermitian lattice.** s(x, y) = 2·Tr_{K/Q(i)}(ξxȳ) is a unimodular skew-hermitian form on L with values in Z[i].
- **Parity.** s(x, x) lies in iZ. The lattice is **even** if it lies in 2iZ for all x. The form I of §2 is even, so the fourfolds of six-line surfaces are the even ones. Odd unimodular lattices of the same signature exist (the fourth power of the elliptic curve C/Z[i] with i acting by (i, i, −i, −i)); the sources do not mention them.
- **A standard basis.** For an even lattice the script finds b₁, …, b₄ with s(b₃, b₁) = s(b₄, b₂) = 1 and the other pairings 0: an isotropic vector from PARI's `qfsolve` on the form of rank 8, a partner by the Euclidean algorithm in Z[i], the orthogonal complement, and the same again in rank 2. It asserts the Gram matrix.
- **The period point.** With B₁ = (φ_r(b_k)) for r in {l, m} and k = 1, 2, and B₂ the same for k = 3, 4,

  > τ = conj(B₁⁻¹·B₂).

  The script asserts that τ lies in D. Another standard basis moves τ by Γ, and the J_k do not change.

**The orientation.** From the description of Lombardo–Peters–Schütt I first took τ = −B₁⁻¹B₂, the plane of H₁ on which i acts holomorphically. That is the complex conjugate of the formula above and gives the conjugate surface: the set of CM points is the same, and every test of §6 and of §7.2 passes with either. The test of §7.3 separates them, and I chose the formula by it, on the field of q = 23, e₁ = 10 (12 cases at four primes; the 8 at the first three primes are in the table of §7.3). I have no derivation of the sign that I trust; a plausible reason is that the period of the surface lives in cohomology and the plane of Lombardo–Peters–Schütt in homology. This is the sign that [note 6](06-torsors-under-cm-class-groups.md) leaves open in Proposition D.1. The comparison of [note 24](24-the-form-mod-4-is-the-character.md) does not depend on it: the other choice replaces a labelling of a torsor by its composite with x ↦ g₀ − x, which sends cosets of a subgroup to cosets ([note 25](25-checks-of-notes-20-23-24.md), finding 12).

---

## 4. Checks that do not use CM [E]

`sixlines_cmmethod.py selftest`, at 60 digits:

| check | result |
|---|---|
| the fifteen relations (2.2) on the products of minors of a random integer matrix | 15 of 15, exactly |
| J₂, …, J₆ under the 720 relabellings of that matrix | unchanged, exactly |
| [J₂ : … : J₆] of the associated matrix | equal, exactly |
| the relations (2.2) on the ten Θ² at three random points of D (not symmetric) | defect below 10⁻⁷⁰ |
| the six other Θ[ab] at these points | below 10⁻⁷⁷ |
| J_k(gτ) = det(g)^k·det(g₂₁τ + g₂₂)^{2k}·J_k(τ) for the five generators of Γ of Clingher–Malmendier–Shaska (2.13) | relative defect below 10⁻⁷¹ |
| J_k(ᵗτ) = J_k(τ) | below 10⁻⁷⁵ |
| a curve of genus 2, y² = ∏(x − p_i) with p = (−3, −2, −1, 1, 2, 3): its period matrix by numerical integration, against the six lines (p_i², −p_i, 1) tangent to a conic | [J(Θ²)] = [J(lines)] to 2·10⁻⁷², and J₄/J₂² below 10⁻⁷⁴ |

The first four rows test the table of §2 and the numbering of the ten products together: a wrong match breaks the linear relations. The last row is the inverse of the period map on the locus of Kummer surfaces, where the period matrix is that of the curve; its homology basis is chosen by trying sign patterns and keeping those that give a symmetric matrix in the Siegel space, and the comparison is between points of the weighted projective space, so the basis does not matter. Off that locus the inverse of the period map is tested only by §6.

---

## 5. Which lattices [D, E]

### 5.1 The maximal order gives no six-line surface [E]

For the 21 polynomials of the enumerated surfaces with maximal order on the K3 side and χ·χ̄ irreducible, the script lists the classes (𝔞, ξ) of the maximal order of K for the three types (`points --maximal`, and a loop over the fields that is not a command of the script; see Reproduce).

| | fields | lattices | even | odd |
|---|---|---|---|---|
| maximal order of K | 21 | 234 | 0 | 234 |

Their number per type is half the number of surfaces in each of the 21 fields. So the "kernel of order 2" of note 15 §3.2 counts objects that are not surfaces of this family.

### 5.2 The order [D]

Over F_q the fourfold has π and q/π as endomorphisms. If the six lines are rational, Frobenius fixes the ordering, which by §2 is a basis of L/(1 + i)L up to the conjugate structure; so π acts trivially on L/(1 + i)L and σ = (π − 1)/(1 + i) is an endomorphism. This is the fourfold form of the condition "(F − 1)/2 integral" of note 13 §3.2. The quadratic twist replaces π by iπ and gives the same condition. So

> R = Z[i, π, q/π, σ, σ̄].

In the fields where 2 is totally ramified, [O_K : R] = 8, against 512 for Z[i, π, q/π].

**Enumeration.** Every lattice of R is, up to a scalar of K, of the form 𝔞 ∩ (u·M′ + 8·O_K) with 𝔞 an ideal class representative prime to 2, M′ one of the lattices between R and O_K, and u a unit at 2 modulo the stabiliser of M′. The script lists the M′ (as submodules of O_K/R), the units (as the lattices uR), keeps the M′ that are isomorphic to their dual at 2, solves for ξ among the generators of the right ideal times units, and asserts ξL = conj(L^∨) for every lattice it returns. Equivalence under units is removed exactly.

### 5.3 Counts [E]

For q = 17, e₁ = 2, where there are 2 surfaces:

| order | index in O_K | lattices between it and O_K | polarised lattices, three types | even, per type | odd, per type |
|---|---|---|---|---|---|
| O_K | 1 | 1 | 3 | 0 | 1 |
| R | 8 | 4 | 21 | 2 | 5 |
| Z[i, π, q/π] | 512 | 64 | 597 | 78 | 121 |

Without the level condition there are 78 even lattices per type: they are the surfaces over F₁₇ with this polynomial whose lines are permuted by Frobenius [D, not checked against an enumeration]. The table of §6 gives, for each field, the number of even lattices per type next to the number of surfaces. They agree in every field. The even lattices have index 2 in their O_K-span, except for q = 37, e₁ = 18, where 2 is not totally ramified in K, [O_K : R] = 64 and the index is 16.

---

## 6. The comparison with the enumerated surfaces [E]

For each field: the three types give 3N points of D. For a coordinate x (J₄/J₂², or another of the ratios below when its class polynomial is not squarefree modulo 𝔭), the polynomials ∏(X − x) of the three types are the images under the three embeddings F ↦ π_jπ_k of one polynomial H with coefficients in E; the script solves for the coefficients on the integral basis of E and accepts them when they are rationals of small height at the working precision. The same is done for Σ_j I(j)·∏_{l≠j}(X − x_l) for each of the five invariants I = J₄/J₂², J₆/J₂³, J₃²/J₂³, J₃J₅/J₂⁴, J₅²/J₂⁵. These are reduced at the prime 𝔭 = (q, F − e₁) of E, where F is a unit, and the five invariants at each root of H mod 𝔭 are compared with those of the configurations of `sixlines_cm.jsonl`, computed from their lines over F_q.

| q | e₁ | Cl(K) | G | [O_K : R] | even lattices per type | odd per type | digits | surfaces enumerated (configurations) | CM points modulo 𝔭 | in both | on one side only |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 17 | 2 | Z/2 | Z/2 | 8 | 2 | 5 | 300 | 2 (4) | 2 | 2 | 0 |
| 19 | 2 | Z/2 | Z/2 | 8 | 2 | 9 | 300 | 2 (4) | 2 | 2 | 0 |
| 19 | 6 | Z/2 | Z/2 | 8 | 2 | 1 | 300 | 2 (4) | 2 | 2 | 0 |
| 23 | 2 | Z/6 | Z/6 | 8 | 6 | 3 | 700 | 6 (12) | 6 | 6 | 0 |
| 23 | 6 | Z/2 | Z/2 | 8 | 2 | 9 | 300 | 2 (4) | 2 | 2 | 0 |
| 23 | 10 | Z/4 | Z/4 | 8 | 4 | 18 | 300 | 4 (8) | 4 | 4 | 0 |
| 29 | 6* | Z/6 | Z/6 | 8 | 6 | 15 | 700 | 6 (12) | 6 | 6 | 0 |
| 29 | 10 | Z/4 | Z/4 | 8 | 4 | 10 | 1,200 | 4 (8) | 4 | 4 | 0 |
| 31 | 2 | Z/8 × Z/2 | Z/4 × Z/2 | 8 | 8 | 4 | 1,200 | 8 (16) | 8 | 8 | 0 |
| 31 | 6 | Z/8 | Z/8 | 8 | 8 | 36 | 1,200 | 8 (16) | 8 | 8 | 0 |
| 31 | 14 | Z/8 | Z/8 | 8 | 8 | 4 | 1,200 | 8 (16) | 8 | 8 | 0 |
| 31 | 18 | Z/4 | Z/4 | 8 | 4 | 2 | 700 | 4 (8) | 4 | 4 | 0 |
| 31 | 22 | Z/2 | Z/2 | 8 | 2 | 9 | 300 | 2 (4) | 2 | 2 | 0 |
| 37 | 18 | Z/4 | Z/4 × Z/2 | 64 | 8 | 20 | 1,200 | 8 (16) | 8 | 8 | 0 |
| 37 | 22 | Z/8 | Z/8 | 8 | 8 | 20 | 1,200 | 8 (16) | 8 | 8 | 0 |
| all | | | | | 74 | 165 | | 74 (148) | 74 | 74 | 0 |

- The row marked \* has a surface with J₂ ≡ 0 mod q. There the script uses seven ratios with J₃ in the denominator (J₂³/J₃², J₂J₄/J₃², J₆/J₃², J₄J₅/J₃³, J₄³/J₃⁴, J₅³/J₃⁵, J₂²J₅/J₃³) and compares those.
- "Digits" is the first working precision of 300, 700, 1,200 at which H and the five polynomials were recognised. For q = 29, e₁ = 10 this is 700; the record used is a run at 1,200, made so that the polynomials of §7.3 are recognised as well.
- The assignment of the three types to the three embeddings was not searched: it is the one of §3.
- The denominators of H on the integral basis of E contain 2, 5 and a few other primes; the coefficients are integral at 𝔭 in every field.
- An agreement by chance of one surface has probability about q⁻⁵ per candidate.

---

## 7. The action

### 7.1 The dictionary as a statement about Frobenius [D]

The reflex field of (K, type T) is E, embedded by F ↦ π_jπ_k, and its reflex type consists of the embeddings F ↦ π_1π_n, n ≠ 1, when K is embedded by φ₁. Let ℓ ≡ 1 mod 4 split completely in K, with roots a₁, …, a₄ of χ modulo (ℓ, i − i_ℓ), 𝔓_m = (ℓ, i − i_ℓ, π − a_m), and 𝔩_jk = (ℓ, F − a_j a_k) in E. Then the reflex norm of 𝔩_jk is

> 𝔟_jk = 𝔓_j·𝔓_k·conj(𝔓_l)·conj(𝔓_m),

the ideal 𝔞 of note 15 §3.3. By the theorem of Shimura and Taniyama [M], the Frobenius of a prime above 𝔩_jk sends the CM point of (L, ξ) to that of (𝔟_jk⁻¹L, ℓξ). So the statement of note 15 §3.3, with the correction of note 18, finding 2, is: *the step of the surface at 𝔩_jk is the isogeny with kernel A[𝔟_jk], and both are Frobenius conjugation in characteristic ℓ.* The map of tori ρ of note 18 is the inverse direction of the reflex norm.

For a lattice of R the ideal acts by L ↦ L + 8·𝔟⁻¹·(L·O_K), which changes L only at the primes of 𝔟.

### 7.2 On lattices [E]

For each field, three or six primes ℓ (the first that split completely and give six distinct products a_j a_k; they lie between 73 and 4,561). For each ℓ and each of the six pairs, the script computes the permutation of the N lattices of each type induced by 𝔟_jk, by locating the image among the known points numerically at 50 digits, and the class of 𝔩_jk in G = Cl(E)/H computed by PARI as in note 15.

| q | e₁ | G | primes ℓ | ideals 𝔟_jk | conflicts | classes of G seen | sums that compose | fixed points of non-trivial classes |
|---|---|---|---|---|---|---|---|---|
| 17 | 2 | Z/2 | 761, 1,753, 1,873 | 18 | 0 | 2 of 2 | 4 of 4 | 0 |
| 19 | 2 | Z/2 | 293, 373, 389 | 18 | 0 | 2 of 2 | 4 of 4 | 0 |
| 19 | 6 | Z/2 | 449, 877, 1,049 | 18 | 0 | 2 of 2 | 4 of 4 | 0 |
| 23 | 2 | Z/6 | 73, 613, 1,453 | 18 | 0 | 6 of 6 | 36 of 36 | 0 |
| 23 | 6 | Z/2 | 1,277, 1,297, 1,669 | 18 | 0 | 2 of 2 | 4 of 4 | 0 |
| 23 | 10 | Z/4 | 257, 317, 1,277 | 18 | 0 | 4 of 4 | 16 of 16 | 0 |
| 29 | 6 | Z/6 | 353, 409, 1,721, 1,901, 2,081, 2,137 | 36 | 0 | 6 of 6 | 36 of 36 | 0 |
| 29 | 10 | Z/4 | 1,741, 2,129, 2,557, 2,617, 2,713, 3,617 | 36 | 0 | 4 of 4 | 16 of 16 | 0 |
| 31 | 2 | Z/4 × Z/2 | 313, 617, 797, 1,093, 1,693, 2,273 | 36 | 0 | 7 of 8 | 42 of 42 | 0 |
| 31 | 6 | Z/8 | 337, 373, 661 | 18 | 0 | 8 of 8 | 64 of 64 | 0 |
| 31 | 14 | Z/8 | 181, 1,009, 1,373 | 18 | 0 | 8 of 8 | 64 of 64 | 0 |
| 31 | 18 | Z/4 | 229, 409, 613 | 18 | 0 | 4 of 4 | 16 of 16 | 0 |
| 31 | 22 | Z/2 | 761, 1,753, 1,873, 3,089, 3,253, 3,469 | 36 | 0 | 2 of 2 | 4 of 4 | 0 |
| 37 | 18 | Z/4 × Z/2 | 941, 2,713, 4,297, 4,357, 4,457, 4,561 | 36 | 0 | 4 of 8 | 16 of 16 | 0 |
| 37 | 22 | Z/8 | 281, 701, 1,021, 2,309, 3,257, 3,449 | 36 | 0 | 8 of 8 | 64 of 64 | 0 |
| all | | | | 378 | 0 | 69 of 74 | 390 of 390 | 0 |

- **No conflict**: two ideals whose primes of E have the same class in G induce the same permutation.
- **A free action**: the permutations of two classes compose to that of the sum; distinct classes give distinct permutations; a non-trivial class has no fixed point.
- **Transitive in 14 fields**: in 13 fields every class of G was seen, and for q = 31, e₁ = 2 the seven classes seen generate G. There G acts simply transitively on the N lattices of each type.
- **One field with half of G**: for q = 37, e₁ = 18, G = Z/4 × Z/2 and the 36 ideals give only the cyclic subgroup of order 4 generated by (1, 1). For this polynomial, and for no other of the 15, note 19 records that d = −1 gives a character of G (`sixlines_characters_genus.jsonl`, field `rational`), whose value on a step is ℓ mod 4; all ℓ used here are 1 mod 4. So the two orbits are what that character predicts, and joining them needs primes ℓ ≡ 3 mod 4, where the primes of K have degree 2 and the script has no code. In PARI's group, for both classes of ℓ, the classes of the primes of degree 1 of E with ℓ ≡ 1 mod 4 are this subgroup of order 4, and those with ℓ ≡ 3 mod 4 are its other coset (84 primes). The two orbits are also the roots of the two factors of the class polynomial over E(i), and the two values modulo 4 of the polarisation ξ at the prime of degree 1 above 2 (all three established in [note 25](25-checks-of-notes-20-23-24.md), finding 13, §4.3 and §4.4).
- This is the first check of Lemma 9 of note 13 on the fourfold side: the group is Cl(E) modulo the classes from E₀ and the ramified primes, and it has the order of the set.

**What becomes of the kernel and cokernel of note 15.** They compared Cl(K)/Cl(K⁺), for the maximal order, with G. The lattices of the surfaces are not lattices of the maximal order (§5), and the isogenies that move them have multiplier ℓ, so they are not elements of Cl(K)/Cl(K⁺) (note 15 §3.3 said this of the cokernel). On the lattices that do correspond to surfaces, no two give the same surface, and in 14 of the 15 fields every class of G is reached by ideals 𝔟_jk with ℓ ≡ 1 mod 4; in the fifteenth these reach the kernel of the one character that note 19 attached to the multiplier.

### 7.3 Frobenius in characteristic ℓ [E]

For a permutation s of the points, A_s = Σ_j x_{s(j)}·∏_{l≠j}(X − x_l) is recognised over E as in §6. If s is the action of the Frobenius of 𝔩, then A_s ≡ X^ℓ·H′ modulo (𝔩, H). The script tests this for s the permutation of 𝔟_jk and 𝔩 = 𝔩_jk, for all ℓ and pairs, when H is squarefree modulo 𝔩 and A_s was recognised. "Not tested" are the cases where one of the two fails; for q = 37, e₁ = 18, which has most of them, A_s, whose coefficients are larger than those of H, was not recognised at 1,200 digits.

| q | e₁ | G | congruence holds | not tested | inverse permutation fails | another pair with a different permutation fails |
|---|---|---|---|---|---|---|
| 17 | 2 | Z/2 | 18 of 18 | 0 | 0 of 0 | 48 of 48 |
| 19 | 2 | Z/2 | 18 of 18 | 0 | 0 of 0 | 48 of 48 |
| 19 | 6 | Z/2 | 18 of 18 | 0 | 0 of 0 | 48 of 48 |
| 23 | 2 | Z/6 | 14 of 14 | 4 | 8 of 8 | 52 of 52 |
| 23 | 6 | Z/2 | 18 of 18 | 0 | 0 of 0 | 48 of 48 |
| 23 | 10 | Z/4 | 18 of 18 | 0 | 8 of 8 | 64 of 64 |
| 29 | 6 | Z/6 | 34 of 34 | 2 | 24 of 24 | 136 of 136 |
| 29 | 10 | Z/4 | 36 of 36 | 0 | 20 of 20 | 120 of 120 |
| 31 | 2 | Z/4 × Z/2 | 36 of 36 | 0 | 16 of 16 | 148 of 148 |
| 31 | 6 | Z/8 | 18 of 18 | 0 | 12 of 12 | 76 of 76 |
| 31 | 14 | Z/8 | 18 of 18 | 0 | 12 of 12 | 72 of 72 |
| 31 | 18 | Z/4 | 18 of 18 | 0 | 8 of 8 | 48 of 48 |
| 31 | 22 | Z/2 | 36 of 36 | 0 | 0 of 0 | 96 of 96 |
| 37 | 18 | Z/4 × Z/2 | 10 of 10 | 26 | 0 of 0 | 0 of 0 |
| 37 | 22 | Z/8 | 36 of 36 | 0 | 28 of 28 | 152 of 152 |
| all | | | 346 of 346 | 32 | 136 of 136 | 1,156 of 1,156 |

- **Controls.** With s⁻¹ in place of s, where they differ, the congruence fails in every case. With the permutation of another pair at the same ℓ, it fails whenever that permutation differs and holds whenever it is the same.
- The left side, X^ℓ, is computed in characteristic ℓ from the class polynomial alone; the right side comes from the lattices and the ideal. They share H.
- This is the route "characteristic ℓ" of note 15 §5.1, and with §6 it is the first test of the action at a prime other than p.
- The orientation of §3 was chosen so that this holds on one field. On the other fields it is a prediction; the column "inverse fails" counts the cases that would have exposed a wrong sign.

### 7.4 The ten other ideals [D, E]

Note 18, finding 2: of the sixteen ideals 𝔟_S = ∏_{j∈S}𝔓_j·∏_{j∉S}conj(𝔓_j), only the six with |S| = 2 are in the domain of ρ, and the others were taken "not to be steps of the surface". They are steps [E]: each of them maps an even lattice of R to an even lattice of R of the same type, so it permutes the surfaces, and since it commutes with G it is an element g_S of G wherever G is transitive.

[D] 𝔟_S² times the power 2 − |S| of the prime of Q(i) above ℓ is in the domain of ρ, and a scalar does not move a lattice. So

> g_S² = ∏ [𝔩_jk]^{d_jk} over the three splittings {j,k} ∪ {l,m}, d_jk = |S ∩ {j,k}| − |S ∩ {l,m}|.

For |S| = 1 this is the product of the three classes [𝔩_jk] with j ∈ S, for |S| = 3 the product of the three with j, k ∈ S, and for |S| = 0 or 4 it is trivial.

| ideals with \|S\| ≠ 2 | permutation is that of a class seen in §7.2 | square tested | square is the predicted class |
|---|---|---|---|
| 630 | 610 | 594 | 594 |

The 20 permutations that are not those of a class seen in §7.2 are in the field q = 31, e₁ = 2, where one class was not seen. Which of the square roots g_S is, the formula does not say.

### 7.5 Labels over F_q [E]

The class polynomial modulo 𝔭 = (q, F − e₁) has N roots in F_q, the surfaces. For each class c seen in §7.2 the polynomial A_s of its permutation, reduced modulo 𝔭, sends a root to a root. Starting from one root with label 0, the script labels the others and then checks every class on every root. For q = 37, e₁ = 18 the polynomials A_s were not recognised and only the first surface has a label.

| fields | surfaces | labelled | a class moves a labelled surface to the surface with the translated label |
|---|---|---|---|
| 15 | 74 | 67 | 389 of 389 |

The six surfaces of q = 23, e₁ = 2, where G is cyclic of order 6:

| label in Z/6 | (J₄/J₂², J₆/J₂³, J₃²/J₂³, J₃J₅/J₂⁴, J₅²/J₂⁵) mod 23 | configurations (b, c, e, f) of note 2 |
|---|---|---|
| 0 | (1, 21, 1, 1, 1) | (2, 3, 9, 7), (2, 3, 9, 19) |
| 1 | (4, 21, 12, 19, 9) | (2, 3, 8, 7), (2, 6, 9, 20) |
| 2 | (2, 6, 5, 4, 17) | (2, 3, 5, 18), (2, 3, 7, 19) |
| 3 | (2, 1, 21, 10, 19) | (2, 3, 9, 5), (2, 3, 9, 21) |
| 4 | (8, 2, 17, 4, 5) | (2, 3, 17, 4), (2, 3, 20, 6) |
| 5 | (12, 17, 18, 8, 1) | (2, 6, 13, 21), (2, 3, 11, 9) |

A label is defined up to a translation (the choice of the first surface, or of the prime above 𝔭 in the class field) and up to the orientation of §3.

### 7.6 A first use: three square classes are not characters [E]

Note 19 could only ask whether a function on a torsor is constant or balanced. With labels one can ask whether it is a character: f(x) = ±ψ(label of x) for a character ψ of order 2 of G. J₂, J₂J₆ and J₃J₅ have even weight, so their square classes in F_q depend only on the surface; they are read off the invariants of §6. On the torsors with at least 4 surfaces, an even group, all surfaces labelled and J₂ ≠ 0:

| | torsors | not a character, not constant | a character | constant | has a zero |
|---|---|---|---|---|---|
| square class of J₂ | 8 | 6 | 0 | 0 | 2 |
| square class of J₂J₆ | 8 | 6 | 1 | 0 | 1 |
| square class of J₃J₅ | 8 | 5 | 0 | 1 | 2 |

The two hits are on groups of order 4, where one sign pattern in four is a character or constant. This agrees with note 19: nothing of this kind computes a character. It is three functions on eight torsors.

### 7.7 Why this is the action on the surfaces [D, sketch]

That this permutation is the class-group action on the surfaces of [note 13](13-hardness-on-paper.md) needs the CM points to be the canonical lifts of the surfaces: they have complex multiplication and reduce to them, and [D, sketch] they are defined over a field unramified at 𝔭 because the class polynomial is separable modulo 𝔭; note 18, finding 3, says this is what is needed.

---

## 8. Check and limits

- **Not a proof.** Every statement of §§5 to 7 is a count on the fields listed.
- **Sources.** The formulas were read in Matsumoto–Terasoma and Clingher–Malmendier–Shaska. Matsumoto's paper was not read. One misprint was met (§2).
- **One fitted sign.** The orientation of §3.
- **Recognition.** A coefficient is accepted when a rational of total size below 0.6 times the working digits fits it to within 10⁻⁽ᵈ⁻¹⁵⁾. A wrong acceptance would show as a mismatch in §6 or a failure in §7.3; there was none. Where a polynomial A_s was not recognised at the precision reached, the case is counted as skipped, not as passed.
- **Numerical location.** The permutation of §7.2 is found by comparing five invariants at 50 digits, with the nearest point closer than 10⁻¹² and the next farther than 10⁻⁶ (asserted).
- **Class groups** are PARI's without certification (GRH).
- **Coverage.** Of the 24 polynomials of the enumerated surfaces with maximal order, 3 have χ·χ̄ reducible and were not treated, 6 have 10, 12 or 26 surfaces and were skipped for time, and 15 were run. Only primes ℓ ≡ 1 mod 4 that split completely in K were used, three or six per field.
- **The enumeration of §5.2** rests on my argument that the list is complete. Its support is that the count equals the number of surfaces in every field.
- **Level and rationality.** That rational lines mean π ≡ 1 on L/(1 + i)L is derived, not read; its support is the same count.
- **My own [D] statements**: the reflex type and reflex norm of §7.1, the formula of §7.4, the description of the action on a lattice of a non-maximal order, the remark on canonical lifts.
- **Re-derivation check** ([note 25](25-checks-of-notes-20-23-24.md), findings 11 to 13). The formulas of §2 agree with the two sources, including the misprint; the sixteen theta constants agree with a direct sum; the self-test, the three lattice counts of §5.3 and two fields of §§6 and 7 were rerun with the same result; every total of §§6 and 7 was recomputed from the stored file. Not checked: that the list of lattices of §5.2 is complete, and the count of 234 in §5.1.

---

## 9. What this does not give

- **A test of note 14.** The lattice of the Prym variety of the curve of note 14 is not computed. The fourfold here is the one that Lombardo–Peters–Schütt attach to the periods of the surface; that Paranjape's Prym variety is isogenous to it is in the sources, and nothing here compares the two beyond what notes 14 and 16 did.
- **Anything about hardness.** The class polynomial has degree |G|, and the working precision needed went from 300 digits for |G| = 2 to 1,200 for |G| = 8, with about 17 minutes for the 24 period points of one field at that precision. At the sizes of note 13 §3.3 the method cannot be run. It is the analogue of computing a Hilbert class polynomial to evaluate the action of CSIDH.
- **A configuration, as opposed to its associate.** That needs the form of weight 4 with the odd character and the fourfold with its i, not up to the conjugate structure.
- **The six lines themselves over F_q.** Only the five invariants are computed; the lines were taken from the enumeration.
- **Polynomials with a non-maximal order on the K3 side** (20 of the 44), **with χ·χ̄ reducible** (3 of the 24 with maximal order), and **fields of p² elements**. Over F_{p²} Frobenius conjugation of the enumerated surfaces could be compared with the ideal above p directly.
- **Which element of G an ideal with |S| odd is**, beyond its square.
- **A derivation of the orientation**, and of the statement on canonical lifts.
- **The candidates of note 19 against the labels.** §7.6 tests three new functions; the conic invariant, the quartic characters and the other candidates of note 19 §3 were not rerun.
- **Primes ℓ ≡ 3 mod 4.** The second orbit for q = 37, e₁ = 18 is identified without them in [note 24](24-the-form-mod-4-is-the-character.md) (next item).
- **Four of the six fields with 10 to 26 surfaces.** Two of the six, q = 29 and q = 37 with e₁ = 2 (12 surfaces each), are done in [note 24](24-the-form-mod-4-is-the-character.md): with the period points computed in parallel they take 75 and 360 seconds, and the class polynomial modulo 𝔭 has exactly the 12 enumerated surfaces as roots. That makes 98 surfaces of 98 in 17 fields. In both, and for q = 37, e₁ = 18, the two orbits of the ideals above primes ℓ ≡ 1 mod 4 are separated on the surfaces over F_q by the invariant δ of [note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), so the second orbit is identified without a prime ℓ ≡ 3 mod 4. For q = 29, e₁ = 2 none of the eight coordinates of this script separates the 12 surfaces modulo 29; `sixlines_delta_labels.py` tries more.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_cmmethod.py selftest --digits 60 > sixlines_cmmethod_selftest.log
python -I sixlines_cmmethod.py sweep 17 19 23 29 31 37 --digits-list 300 700 1200 --primes 3 --all-subsets --maxn 8 --workers 4 --out sixlines_cmmethod.jsonl > sixlines_cmmethod.log
python -I sixlines_cmmethod.py sweep 29 31 37 --select 29,6 29,10 31,2 31,22 37,18 37,22 --digits-list 300 700 1200 --primes 6 --all-subsets --maxn 8 --workers 4 --out sixlines_cmmethod.jsonl > sixlines_cmmethod_rerun.log
python -I sixlines_cmmethod.py sweep 29 --select 29,10 --digits-list 1200 --primes 6 --all-subsets --workers 1 --out sixlines_cmmethod.jsonl > sixlines_cmmethod_rerun2.log
python -I sixlines_cmmethod.py points 17 2 1 1 --digits 30 --maximal
python -I sixlines_cmmethod.py points 17 2 1 1 --digits 30 --nolevel
python -I sixlines_cmmethod.py action 23 2 1 1 --digits 700 --primes 3 --all-subsets
```

A few seconds for the self-test; about 45 minutes for the first sweep and 30 for the second on a loaded machine with 4 workers, 12 minutes for the third; two minutes each for the two `points` commands. The second sweep repeats six fields after two corrections to the script (the surfaces of a polynomial are now selected by e₁ and e₂, where the first run selected by e₁ alone and counted the surfaces of a second polynomial at q = 31 and 37; and the ratios with J₃ were added) and with six primes; the third repeats one field at a fixed precision. A sweep appends one record per field and command to `sixlines_cmmethod.jsonl`, and the last record of a field is the one used here; the text goes to the log. It needs `sixlines_cm.jsonl` and the files `sixlines_p{p}.jsonl`. No package was installed: the script uses cypari2 and numpy only. The count of §5.1 was made by a loop over the fields calling `CMField(q, e1, mu, nu).lattices()`; it is not a command of the script.

## Sources

- Matsumoto, Terasoma, *Thomae type formula for K3 surfaces given by double covers of the projective plane branching along six lines*, https://arxiv.org/abs/1001.4865 : §2.1, §2.2, §2.5, §3.1 to §3.4, Table 1, Propositions 6 and 7, Remarks 1 and 2 [P].
- Clingher, Malmendier, Shaska, *Six line configurations and string dualities*, https://arxiv.org/abs/1806.07460 : §2 to Corollary 2.25 [P].
- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, https://arxiv.org/abs/1209.5997 : §1.1, §1.2, §2.1, §2.2 [P].
- Not read: Matsumoto, *Theta functions on the bounded symmetric domain of type I₂,₂ and the period map of a 4-parameter family of K3 surfaces*, Math. Ann. 295 (1993), 383–409; Matsumoto, Sasaki, Yoshida, *The monodromy of the period map of a 4-parameter family of K3 surfaces and the hypergeometric function of type (3, 6)*, Internat. J. Math. 3 (1992).
- Recalled, not re-read [M]: the main theorem of complex multiplication (Shimura, Taniyama).
- In this repository: notes [2](02-six-line-count.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [18](18-checks-of-notes-13-to-16.md), [19](19-quadratic-characters-from-a-surface.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md).

