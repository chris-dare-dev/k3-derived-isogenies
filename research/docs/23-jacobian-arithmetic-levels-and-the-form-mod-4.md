# 23. Arithmetic in the Jacobian of the genus-5 curve, levels, and the form mod 4

*Every test of the correspondence between six-line surfaces and abelian fourfolds in notes [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md) and [18](18-checks-of-notes-13-to-16.md) uses point counts, which are isogeny invariants. This note supplies the group law on the Jacobian of the curve of note 14 and uses it for two tests that point counts cannot make: the levels of notes 7 and 8, and the character χ₋₁ that [note 19](19-quadratic-characters-from-a-surface.md) §1 item 7 leaves uncomputed. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_jacobian.py`; **[H]** heuristic. Checked by re-derivation in [note 25](25-checks-of-notes-20-23-24.md) (findings 4 and 7 to 10); otherwise the derivations are supported only by the computations that come out.*

---

## 1. Result

1. **The arithmetic exists and passes its checks [E].** `sixlines_jacobian.py` computes in Jac(C) for the curve C : v⁴ = F of note 14 and its quartic twists, over F_p and F_{p^r}, p ≡ 1 mod 4: the group law, J, Frobenius, a canonical form, the structure of ℓ-primary parts, and the Tate pairing. Divisor classes are spaces of sections of a line bundle of degree 12 (Khuri-Makdisi's method), stored as power series at one rational point. On 35 curves over F_p, F_{p²}, F_{p³} and F_{29⁴}, the value at 1 of the Frobenius polynomial kills all 70 random classes, J and Frobenius satisfy their polynomials, and χ(Frobenius) vanishes on the Prym part with one of the two identifications of i with J and never with the other (§3). No computer algebra system and no package was needed.
2. **Levels: the prediction, and it holds [D, E].** The rule Λ² gives one statement for all six pairs (polynomial, odd ℓ) with rational lines and two levels: *a surface is on the lower level at ℓ exactly when Frobenius acts through a single Jordan block on the part of A[𝔩] at the multiple root of χ modulo 𝔩; on the upper level there are two blocks* (§4.1). For three of the six a quartic twist makes that root equal to 1, and the number of blocks is then read off the group Prym(F_p)[ℓ^∞]. Result: 3 polynomials, 29 surfaces, 53 configurations, 422 pairs of lines computed; on every surface all pairs agree; the surfaces with two blocks and with one block number 1 and 6, 5 and 14, 1 and 2, and the lattice classes on the upper and lower level number 1 and 6, 5 and 14, 1 and 2. No exception. The other three polynomials need the group over F_{p⁸}, F_{p¹⁵}, F_{p²⁴} and were not done (§4.4).
3. **So the fourfold of note 14 carries more than the zeta function of the surface [E].** All surfaces of one polynomial have the same point counts, on both sides. The endomorphism ring at ℓ of the fourfold separates them in the way the lattices of note 7 say. This is the first agreement between the two sides on an invariant that is not an isogeny invariant.
4. **χ₋₁: an invariant of the fourfold that should read it, and its normalisation [D].** Where 2 has a prime of degree 1 in K⁺, the 2-adic Tate module has a line V₁ over Q₂(i) on which Frobenius acts by a scalar π₁, and the class u of the hermitian form of the polarisation on V₁, modulo norms and powers of 2, is 1 or 3 modulo 4. An isogeny that multiplies the polarisation by ℓ multiplies u by ℓ. u changes sign with the conjugate structure; δ = (−1)^{(u−1)/2}·sign(Re c₁·Im c₁)·sign(c₂), with c₁, c₂ the coefficients of χ for the curve, does not, and is unchanged by twists (§5.2). The method of Castryck, Sotáková and Vercauteren was read, and it is of this kind (§5.1).
5. **It can be computed, over F_{p⁴} [D, E].** The projector that polynomials in Frobenius provide is divisible by 8 on V₁ and useless. Instead the kernel of π − π₁ on the rational 16-torsion is found by linear algebra on matrices of Tate pairings, and u is read off a self-pairing of an element of that kernel (§5.3). About 75 seconds for one curve.
6. **δ is an invariant of the surface, in the data [E].** On the three torsors over prime fields where note 19 says χ₋₁ is a character (q = 29, e₁ = 2; q = 37, e₁ = 2; q = 37, e₁ = 18): 32 surfaces, 64 configurations, 242 curves tried, 106 computed. On each of the 32 surfaces all computed pairs give the same δ (2 to 6 pairs each). u itself is the same for all pairs on only 10 of the 32, and each of the two partial normalisations on 12 and 15 of the first 29: the factor that §5.2 derives is the one that works. Twists 1 and 2 and other random classes give the same u on 36 of 36 controls.
7. **δ is balanced on all three torsors [E].** +1 on 6 surfaces and −1 on 6; 6 and 6; 4 and 4. The first run left one surface of each torsor without a usable pair, at 5 to 6, 5 to 6 and 3 to 4; balance predicted +1 for each of the three, and each came out +1 (2, 6 and 5 pairs). δ is none of the candidates of note 19 (this note, §5.4).
8. **What this is, and what it is not.** A function with values ±1 on surfaces that is well defined, is not a function of the polynomial, and takes each value on exactly half of each torsor, as a non-trivial character must. This note does **not** show that it is the character: that needs δ(g·X) = χ₋₁(g)·δ(X) for known g, and here no surface of these torsors is labelled by a class (note 20 labels surfaces by the CM method but skipped the torsors with 12 lattices). A function chosen at random would be balanced on all three torsors with probability about 1.4 %. [Note 24](24-the-form-mod-4-is-the-character.md) makes the test, and δ is the character on all three torsors: with the CM method of note 20 extended to the two torsors of 12 surfaces, the steps at primes ℓ ≡ 1 mod 4 (48 ideals per field) generate a subgroup of index 2 of G with two orbits on the surfaces, and δ is constant on each orbit with opposite values: 12 of 12, 12 of 12 and 8 of 8 surfaces. The function that separates the orbits is i times a function over E, so the character is χ₋₁. A balanced function would agree by chance once in 7,470,540. That is three torsors at two primes; [note 26](26-delta-is-the-character-at-q-41.md) adds two torsors at q = 41 (§5.5). So this is the first information about the class that can be computed from two surfaces in this repository: one bit, for the 80 of 477 polynomials of note 19, all with q ≡ 1 mod 4.

   [Note 25](25-checks-of-notes-20-23-24.md) (findings 4 and 8) recomputed δ on nine curves, three of each torsor, on the twist C₃ with other random classes and with the 2-adic arithmetic written again: 9 of 9 equal the stored value. And u can be read on the lattices of note 20: for an even lattice (L, ξ) it is the unit part modulo 4 of the image of −iξ in Q₂ at the prime of degree 1 above 2. On the 96 lattices of the three fields this splits the surfaces exactly as δ does. So δ is, up to one sign for each polynomial, that invariant of the lattice of the CM lift, for which "a step with multiplier ℓ multiplies u by ℓ" is immediate.
9. **Along the way [E].** The Prym varieties of different pairs of lines of one surface differ at 2: for 126 of the 242 curves, all of them pairs with conjugate branch points, the rational 2-primary torsion over F_{p⁴} has exponent above 16 and the method of §5.3 does not apply. Some twists have no usable point over F_p at all; an unguarded loop over random points hung a run on them.

---

## 2. The arithmetic

### 2.1 The curve and the line bundle [D]

Notation of note 14: E : Y² = W(X) over F_p, F the function on E with divisor 2B − 4(c″) − 4(O), B the four branch points, and C : v⁴ = F. `sixlines_prym.py` supplies W, the points and F as a chain of lines. This note works with

> C_tw : v⁴ = n₀^tw·F, tw = 0, 1, 2, 3, n₀ the least primitive root of p,

the four quartic twists, with J : v ↦ ζv, ζ = n₀^((p−1)/4). Only p ≡ 1 mod 4 is handled (c″ rational, Frobenius commuting with J); the points of B may be conjugate over F_{p²}.

- **Base point.** z is the first rational point of E, outside the support of F and with Y(z) ≠ 0, at which n₀^tw·F is a fourth power; P₀ is a point of C above z and t = X − X(z) a local parameter at P₀.
- **Line bundle.** L = O(3π\*(z)), of degree 12 ≥ 2g + 1 = 11. It is stable under J and under Frobenius, so both act on its sections.
- **Sections.** A function Σ f_j·v^j, j = 0, …, 3, with f_j on E, has poles bounded by m·3π\*(z) exactly when f_j lies in L_E(D_{m,j}), D_{m,j} = 3m(z) + ⌊j·div(F)/4⌋, the floor taken coefficient by coefficient. So

  > H⁰(L^m) = ⊕_j L_E(D_{m,j})·v^j, of dimension 12m − 4.

- **Riemann–Roch on E.** For a divisor D of degree n ≥ 1 with sum S in the group law, let g have divisor D − (n − 1)(O) − (S); `pr_chain` of note 14 writes g as a product of lines. Then L_E(D) = g⁻¹·L_E((n − 1)(O) + (S)), and L_E((n − 1)(O) + (S)) is spanned by the monomials X^a·Y^b of pole order at most n if S = O, and by (μ − μ(−S))/(X − X(S)), μ a monomial of pole order at most n + 1, otherwise.
- **Series.** X, Y and v are expanded at P₀ in t (square and fourth roots of series by Newton iteration). A section s of L^m is stored as the first 64 coefficients of t^{3m}·s; this is injective for m ≤ 5 (degree 60). The product of two sections is the product of the two truncated series. Chains through conjugate points are evaluated over F_{p²}, and the result is checked to have coefficients in F_p.

So the whole curve enters through five matrices of 64 rows with 8, 20, 32, 44 and 56 columns over F_p, and the grading by j. They serve over every F_{p^r}.

### 2.2 The group law [M for the method, D for the details]

This is the "small model" of Khuri-Makdisi. A class x of degree 0 is represented by an effective divisor D of degree 12 with x = [D − 3π\*(z)], stored as the subspace W_D = H⁰(L²(−D)) of H⁰(L²), of dimension 8. Two operations on subspaces are used: the span of the products of two spaces, and {s in H⁰(L^a) : s·W ⊂ U}. Both are linear algebra, and by Riemann–Roch the dimension of every result is known in advance; the script uses a few random elements in place of whole bases and accepts a result only when its dimension is the predicted one.

- **addflip(x₁, x₂) = −(x₁ + x₂).** H⁰(L⁴(−D₁ − D₂)) = W_{D₁}·W_{D₂} (dimension 20); H⁰(L³(−D₁ − D₂)) = {s in H⁰(L³) : s·H⁰(L) ⊂ H⁰(L⁴(−D₁ − D₂))} (dimension 8); for a section s of it, div(s) = D₁ + D₂ + D₃, and W_{D₃} = {u in H⁰(L²) : u·H⁰(L³(−D₁ − D₂)) ⊂ s·H⁰(L²)}.
- **Negation.** For s in W_D, div(s) = D + D′ and W_{D′} = {u : u·W_D ⊂ s·H⁰(L²)}.
- **Multiples** by a binary chain of addflips, keeping track of the sign.
- **Zero.** x = 0 exactly when {u in H⁰(L) : u·H⁰(L) ⊂ W_D} ≠ 0.
- **Canonical form.** Negate; in the result take the section that vanishes to the highest order at P₀ (it is unique up to a scalar); negate with respect to it. This gives the divisor of the class with the largest multiplicity at P₀, and its space in PARI's kernel normal form is a key for hashing. Equality of classes is tested either way.
- **J and Frobenius.** J^e multiplies the part of degree j by ζ^{ej}; Frobenius raises the coefficients of a basis to the p-th power.
- **Classes of points.** Let s₀ be a section of L with a simple zero at P₀ and D₁ = div(s₀) − P₀. H⁰(L²(−D₁)) is computed once (dimension 9); for a rational point R of C the sections in it that vanish at R form W_D for D = R + D₁, the class [R − P₀]. Values at R come from expressions (A(X) + Y·B(X))/Den(X)·v^j with coefficients in F_p, fitted to the series. A "random class" is a sum of three classes of random points over F_{p^r}.
- **Curves with few points.** Some twists have almost no points over F_p (one curve met has none outside the fibre of z, seven have four). There a random class is the trace of a random class over F_{p^{2r}}, brought down through its canonical form, which has entries in the smaller field.

### 2.3 Groups and pairings [D]

- **Structure of an ℓ-primary group** of known order ℓ^w, given by a projection of random classes: relations between the generators are found by discrete logarithms in the subgroup generated so far, by brute force on its elements of order ℓ (hashed by the canonical form), and the invariants are read off the Smith form of the relation matrix. The run stops when the subgroup has order ℓ^w, so a failure to generate is detected.
- **Tate pairing.** For nP = 0, Miller's algorithm in this model (following Bruin): each addition gives sections s of L³ and u of L² with D_a + D_b − D_{a+b} − 3π\*(z) = div(s/(u·1)), where 1 is the constant function as a section of L, and the function f_P with divisor n(D_P − 3π\*(z)) is a product of such ratios. A ratio of two sections of L^a is evaluated on an effective divisor Y of degree 12 as the quotient of the determinants of the two multiplication maps H⁰(L²)/W_Y → H⁰(L^{2+a})/H⁰(L^{2+a}(−Y)). Then T_n(P, Q) = (f_P(D_Q)/f_P(div s′))^{(q−1)/n} for a random section s′ of L. A vanishing determinant (supports that meet) is detected, and the representatives are changed.

### 2.4 Cost

One addflip takes of the order of 30 ms over F_p in PARI through cypari2 on a loaded machine, and a multiple by the group order about a second; over F_{p⁴} everything is several times slower, and one computation of §5 (twelve multiples by a cofactor of 75 bits, 48 Miller functions) takes about 75 seconds. No package was installed.

---

## 3. Validation [E]

`sixlines_jacobian.py --validate-set`: for p = 13, 17, 29, 37, the first two configurations with irreducible P, two stored pairs each, twists 0 and 1 over F_p, twist 0 over F_{p²}, one curve per prime with twist 2 over F_{p³}, and one curve over F_{29⁴}; two random classes x on each curve. 40 runs; 5 were refused by the setup (no base point z: four at p = 13, one at p = 17). On the other 35 curves (25 over F_p, 6 over F_{p²}, 3 over F_{p³}, 1 over F_{29⁴}), 70 classes:

| check | holds |
|---|---|
| x ≠ 0 | 70 of 70 |
| n·x = 0, n = #E(F_q)·\|χ_tw,r(1)\|² the value at 1 of the Frobenius polynomial of Jac(C) | 70 |
| (n + 1)·x = x | 70 |
| J⁴x = x | 70 |
| Frobenius^r x = x, and J∘Frobenius = Frobenius∘J on x | 70, 70 |
| the canonical form of x + 0 is that of x | 70 |
| (J² + 1)(1 − J²)x = 0 | 70 |
| χ(Frobenius) kills (1 − J²)x, with i acting as J³ | 70 |
| the same with i acting as J (a control) | 0 |
| (Frobenius² − a·Frobenius + p) kills (1 + J + J² + J³)x, a the trace of E | 70 |

- The order n is the product of two orders: that of E and that of the Prym part, |χ(1)|²; the quotient by J² is an elliptic curve isogenous to E and adds nothing (note 14 §2.4).
- Over F_p Frobenius is the identity on rational classes, and the check on χ says that the Gaussian integer χ(1) kills the Prym part; over F_{p²} and beyond it is a test of Frobenius itself.
- **The convention.** J : v ↦ ζv on functions acts on classes as the inverse of the element i for which the traces of note 14 are computed: χ holds with i ↦ J³ and fails with i ↦ J. Everywhere below "i" on classes means J³.
- **The pairing.** Bilinearity in both arguments, T_n(P, nQ) = 1, independence of the representatives, T(JP, JQ) = T(P, Q) and T(πP, πQ) = T(P, Q)^p were checked by hand on a few points (n = 4 over F₂₉, n = 8 over F_{29²}); these runs are not logged. The logged check is in the computation of §5: on each of the curves there, the two pairings of the final point computed directly agree with the values assembled from the matrices by bilinearity and J-invariance (§5.4). A logged check of the identities, for one curve over F_{29⁴} with N = 16, is in [note 25](25-checks-of-notes-20-23-24.md) (finding 9, `sixlines_review2.py pairing`): additivity in each argument, T(z, 16R) = 1, invariance under J, and T(πz, πy) = T(z, y)^p all hold.
- **One trap.** PARI's `ffprimroot` is random, so "the" root of unity of order n changed between calls until it was fixed once per curve; before that, logarithms of pairings differed by n/2 from one call to the next.
- Every dimension that Riemann–Roch predicts is asserted at each step of each operation; none failed in any run of this note.

---

## 4. Levels

### 4.1 What the correspondence predicts [D, E]

Let ℓ be an odd prime, ℓ ≠ p, M the lattice of the fourfold at ℓ, of rank 4 over Z_ℓ[i], with π acting, and T the lattice of the surface at ℓ with F acting. The correspondence on lattices (note 15 §2) is T = the fixed part of Λ²_{Z_ℓ[i]}M under the antilinear involution, with F = Λ²π.

- If ℓ ≡ 1 mod 4, Z_ℓ[i] = Z_ℓ × Z_ℓ, M = M_𝔩 ⊕ M_𝔩̄, and the fixed part is Λ²_{Z_ℓ}M_𝔩. So T/ℓT = Λ²(M_𝔩/ℓM_𝔩).
- If ℓ ≡ 3 mod 4, T ⊗ Z_ℓ[i] = Λ²M, and T/ℓT ⊗ F_{ℓ²} = Λ²(M/ℓM).

In both cases the Jordan type of F on T/ℓT is that of the second exterior power of π on a space of dimension 4 over k = Z[i]/𝔩. The primes of interest are those where χ has a multiple root a modulo 𝔩, of multiplicity m, the other roots being simple; the Jordan type of π at a is a partition of m, and it is the only unknown.

`sixlines_jacobian.py --lattice` computes both sides for the six pairs (polynomial, ℓ) with rational lines, ℓ odd, and two levels (the rows marked "rational" of note 7 §5 and two more without rational steps): for every lattice class of notes 2 and 7, the dimensions of ker g(F)^k on L/ℓL for each multiple factor g of P mod ℓ; and the same dimensions for Λ² of each Jordan type.

| q | (e₁, e₂, e₃) | ℓ | P mod ℓ: (degree, multiplicity) | m | classes with End(L) maximal at ℓ: ker dims, type of π | classes of index ℓ²: ker dims, type of π |
|---|---|---|---|---|---|---|
| 29 | (14, 87, −23548) | 13 | (1,1), (1,2), (1,2), (1,1) | 2 | 1: [2, 2], (1, 1) | 6: [1, 2], (2) |
| 29 | (18, 87, −30276) | 3 | (2,3) | 3 | 1: [2, 3, 3], (2, 1) | 2: [1, 2, 3], (3) |
| 37 | (14, −185, −38332) | 11 | (2,2), (2,1) | 2 | 1: [2, 2], (1, 1) | 6: [1, 2], (2) |
| 37 | (34, 1591, 38332) | 5 | (1,1), (1,1), (2,2) | 2 | 5: [2, 2], (1, 1) | 14: [1, 2], (2) |
| 37 | (38, 1591, 27380) | 7 | (1,2), (1,2), (2,1) | 2 | 2: [2, 2], (1, 1) | 8: [1, 2], (2) |
| 37 | (66, 4551, 169756) | 3 | (1,6) | 4 | 1: [2, 4, 6, 6, 6, 6], (3, 1) | 2: [2, 3, 4, 5, 6, 6], (4) |

In every case each class matches exactly one Jordan type, and the rule is the same:

> **Prediction.** A surface is on the lower level at ℓ exactly when π acts on the part of A[𝔩] at the multiple root through a single Jordan block; on the upper level there are two blocks.

One block means that M is cyclic over Z_ℓ[i][π] at the maximal ideal (𝔩, π − a), that is, the endomorphism ring of the fourfold is no larger than Z[i, π] there. This is the fourfold's version of note 7, Proposition I. In the last row the eigenspace of F has dimension 2 on both levels (note 8), and the levels differ only in ker (F + q)²; on the fourfold they differ already in the eigenspace.

### 4.2 What is measured [D]

For a root a of χ modulo 𝔩 and the twist C_tw, Frobenius has the root a·i^tw, with i the image of i in k. If a·i^tw = 1, the part of A_tw[𝔩] at this root that is fixed by Frobenius is the kernel of π − 1 on it, and its dimension over k is the number of Jordan blocks. So

> number of blocks = (rank of Prym(C_tw)(F_p)[ℓ] − contribution of the simple roots equal to 1)/[k : F_ℓ],

and the rank is read off the structure of the ℓ-primary part of the group of rational points, which §2.3 computes. A root of order r > 1 for every twist would need the group over F_{p^r}.

| q | e₁ | ℓ | least order of the multiple root over the four twists | |
|---|---|---|---|---|
| 29 | 14 | 13 | 1 | done |
| 37 | 34 | 5 | 1 | done |
| 37 | 66 | 3 | 1 | done |
| 29 | 18 | 3 | 8 | not done: needs Jac(C)(F_{29⁸}) |
| 37 | 14 | 11 | 15 | not done |
| 37 | 38 | 7 | 24 | not done |

### 4.3 Results [E]

`sixlines_jacobian.py --levels`: for each of the three polynomials, every configuration with polynomial P(t) or P(−t), grouped into surfaces by association, and for each configuration the pairs of lines stored by note 14 (at most 8). For each pair: the twist on which the multiple root is 1, and the structure of the ℓ-primary part of Prym(C_tw)(F_p).

| q | e₁ | ℓ | surfaces | configurations | pairs computed | lattice classes: upper, lower | surfaces with two blocks, one block | groups found, two blocks | groups found, one block |
|---|---|---|---|---|---|---|---|---|---|
| 29 | 14 | 13 | 7 | 13 | 104 | 1, 6 | 1, 6 | (Z/13)²: 8 pairs | Z/169: 96 |
| 37 | 34 | 5 | 19 | 35 | 278 | 5, 14 | 5, 14 | (Z/5)²: 28; Z/25 × Z/5: 27 | Z/25: 117; Z/125: 106 |
| 37 | 66 | 3 | 3 | 5 | 40 | 1, 2 | 1, 2 | (Z/3)⁴: 8 | (Z/9)²: 32 |
| all | | | 29 | 53 | 422 | 7, 22 | 7, 22 | | |

- **Pairs of one surface agree.** On each of the 29 surfaces every computed pair gives the same number of blocks: 422 pairs, no exception. The pairs use all four twists, and the order of the group changes with the pair (5² or 5³ in the second row), but not whether it is cyclic. Two pairs of the second row are degenerate for the script (no base point or a coincidence of points) and were not computed.
- **The counts agree** with the lattice side in all three polynomials: 7 surfaces on the upper level and 22 on the lower, as predicted. The counts weigh less than they look ([note 25](25-checks-of-notes-20-23-24.md), finding 10): the surfaces with two blocks number 1, 5 and 1, and 1, 3 and 1 of them are self-associated; no self-associated surface has one block. Apart from that, the agreement of counts is one number: 2 of the 16 other surfaces of the second polynomial have two blocks, against 2 classes with 2 automorphisms on the upper level. That all pairs of a surface agree (422 pairs) is not affected. The prediction itself is forced by the lattice computation and the rule Λ², not fitted.
- **Which surfaces.** The lattice side does not say which surface is on which level; notes 2 and 7 match surfaces and lattice classes only in number. One finer agreement: the lattice classes with 4 automorphisms are all on the upper level (1, 3 and 1 of them), and the self-associated configurations are 1, 3 and 1 in number and all have two blocks (that a self-associated configuration is a surface with an extra automorphism is my guess at why; I did not check it). In the second row the upper level has two more classes, with 2 automorphisms, and two more surfaces have two blocks.
- **What the agreement rests on.** The prediction uses the lattice classes of note 7 (computed there with no reference to equations) and the rule Λ². The measurement uses the curve of note 14 and the group of its rational points. Point counts cannot give either side: all surfaces of a row have the same zeta function.

This is the first test in this repository of the correspondence surface ↔ fourfold on an invariant that is not an isogeny invariant. It is passed in every case that could be computed.

### 4.4 What stopped the other three

The multiple root has order 8, 15 and 24 for every twist (§4.2), so its eigenspace is rational only over F_{p⁸}, F_{p¹⁵}, F_{p²⁴}. Over F_{29⁸} one group operation would take about a second and one random element of ℓ-power order about ten minutes; the case (29, 18, ℓ = 3) with 3 surfaces is within reach of a longer run and was not done. The other two are not within reach of this implementation.

---

## 5. The character χ₋₁

### 5.1 What Castryck, Sotáková and Vercauteren do [P]

Read in their paper, §3.2 and Theorem 10 (the conference version; note 19 had only the abstract).

- For an ordinary curve E over F_q with 4 | q − 1 they walk down the 2-volcano to a curve E₀ with E₀(F_q)[2^∞] cyclic of order 2^v, take P of order 4 and Q of order 2^v with 2^{v−2}Q = P, and form the reduced Tate pairing T₄(P, Q). It is a primitive fourth root of unity that does not depend on the choice of P and Q (another choice is (kP, kQ) with k odd, and k² ≡ 1 mod 4). For E₀′ = [𝔟]E₀ the compatibility T₄(φP, φQ) = T₄(P, Q)^{deg φ} gives N(𝔟) mod 4 = log_{T₄(P,Q)} T₄(P′, Q′), which is the character δ.
- Their Theorem 10 is a second form: for q ≡ 1 mod 4 and trace ≡ 0 mod 4, with E : y² = x³ + ax² + bx, δ([𝔞]) = (b′/b)^{(q−1)/4}.

So the method attaches to each curve a canonical primitive fourth root of unity in F_q, and an isogeny of degree d raises it to the power d. It uses the cyclic structure of the rational 2-primary torsion to make the root canonical. The analogy of note 19 §2.4 is right in this form: the canonical fourth root of unity that it seeks from a pairing of a point of order 4 is this T₄(P, Q), not a Weil pairing of P with itself.

### 5.2 The invariant on the fourfold [D]

Let 2 be unramified in K⁺ with primes v₁ of degree 1 and v₃ of degree 3, as in the 80 cases of note 19, and w₁ the prime of K above v₁, with K_{w₁} = Q₂(i). On the 2-adic Tate module, V₁ = the part where π acts through its root π₁ in Z₂[i]; it is a line over Q₂(i), orthogonal to the rest for the form E of any polarisation compatible with the multiplication. For a vector m of V₁ put h(m) = E(m, i·m). Then E(αm, βm) = Im(ᾱβ)·h(m) for α, β in Q₂(i), and h(αm) = N(α)·h(m).

> **u(A, λ, ι)** = the class of h(m) in Q₂^×/N(Q₂(i)^×) modulo the powers of 2: the odd part of h(m) modulo 4, in {1, 3}.

- It does not depend on m, nor on a lattice: it is an invariant of the polarised fourfold with its structure ι : Z[i] → End up to isogenies that respect λ up to a power of 2. An isogeny with φ\*λ′ = ℓ·λ multiplies u by ℓ. This is the character of note 19 §2.4, read on the line V₁ and not on a lattice; the Prym variety and its polarisation of type (1, 1, 1, 2) can be used as they are.
- **Conjugate structure.** Replacing ι by its conjugate replaces h by −h, so u by −u. A surface gives ι only up to conjugation (note 19 §1 item 5). The structure is pinned by χ: of the two, one has trace c₁ and the other c̄₁.
- **Twists.** The four twists C_tw have the same Jacobian over F_{p⁴} with the same J and the same line V₁, so the same u; their quartics have c₁·i^tw. The surface with polynomial P (not P(−t)) is the one with c₂ = e₁.
- **The candidate.** For a maximal order c₁ = 4(±1 ± i). So for a curve with data (c₁, c₂) and e₁ > 0 put

  > δ = (−1)^{(u−1)/2}·sign(Re c₁·Im c₁)·sign(c₂).

  The second factor is +1 when c₁ is ±4(1 + i) and −1 when it is ±4(1 − i); the third moves to the twist that belongs to the surface, where c₁ is multiplied by i. δ is unchanged by ι ↦ ῑ and by the twists. It is defined up to one sign for each polynomial (the choice of 1 + i against 1 − i, and of the fourth root of unity ζ in F_p, which enters J, the traces and the pairing together).
- **What is not controlled.** The Prym varieties of two pairs of lines of one surface are isogenous, and the isogeny need not carry one polarisation to a multiple of the other: the compatible polarisations of the isogeny class form a cone in K⁺, and a change by a totally positive β multiplies u by the Hilbert symbol (β, −1) at v₁. So u could depend on the pair. That is the first thing to test.
- **If the correspondence respects the actions** and the polarisations of the Prym varieties of the surfaces of a torsor are matched, δ(g·X) = χ₋₁(g)·δ(X), and δ takes each value on half of the torsor.

### 5.3 How it is computed [D]

**The pairing identity.** Write G = Jac(C)(F_{p^r})[2^∞] = ϖ⁻¹T/T with T the Tate module and ϖ = π^r − 1. For X in G with NX = 0 and any Y in G, with lifts X̃, Ỹ in ϖ⁻¹T,

> T_N(X, Y) = ζ_N^{E(N·X̃, ϖ·Ỹ)}.

If X lies in the kernel of π − π₁ and in the part with J² = −1, write X̃ = γm/ϖ₁ + (a part in the complement of V₁), ϖ₁ = π₁^r − 1. Then

> Z := log T_N(X, i·X) + i·log T_N(X, X) = M·conj(N/ϖ₁) mod N, M = N(γ)·h(m),

provided the complement contributes nothing (below). ϖ₁ is known from χ, so Z gives M modulo N/(1 + i), and u is the odd part of M modulo 4. Two real numbers are measured and one is unknown: that Z is a real multiple of conj(N/ϖ₁) is a check of each measurement. ζ_N is taken with ζ_N^{N/4} = ζ.

**Why the obvious projector fails [D, E].** χ/(X − π₁), evaluated at π, kills the complement of V₁ and maps into it. But π ≡ 1 mod 2 on the whole Tate module of the Prym variety for the right twist (on the one curve where I looked, Prym(C)(F_p)[2^∞] is killed by 2 and has order 2⁸), so the three differences π₁ − π_j are divisible by 2 and this operator is divisible by 8 on V₁: χ′(π₁) has valuation 6 at 1 + i for all three polynomials. On that curve, over F_p, F_{p²} and F_{p⁴}, its image after removing the part of E had exponent at most 2 (trials not logged). The eigenline is rational to order 16 there, but polynomials in π cannot reach it.

**What is done instead.** Over F_{p⁴}, N = 16:

1. twelve random classes y_a of 2-power order and z_a = (1 − J)y_a, which lie in the part with J² = −1 and are killed by 16;
2. the matrices log T₁₆(J^k z_a, y_b), k = 0, 1, 3, and log T₁₆(πz_a, y_b), by Miller's algorithm (48 functions, each evaluated on 13 divisors);
3. the lattice of integer vectors c with Σ c_a(π − π₁)z_a pairing trivially with every y_b, by linear algebra over Z/16; T₁₆ is perfect, so these are the combinations in the kernel of π − π₁;
4. for each basis vector, Z from the matrices (T(z_a, z_b) and T(z_a, i·z_b) are differences of the entries, because J preserves the pairing), and M;
5. for the vector with the least 2-adic valuation of M, the point X itself, with a direct check that X has the stated order, that πX = π₁X, that J²X = −X, and that the two pairings of X computed directly give the same Z.

**The complement [D].** The other elements of the kernel of π − π₁ come from the complement of V₁, where π − π₁ is 2 times an invertible map; their lifts are halves of lattice vectors, and their contribution to E(N·X̃, ϖ·X̃) is 4·E(t, i·ϖt) with ϖ divisible by 8, which is 0 modulo 16. This uses that the complement is orthogonal to V₁ and that π ≡ 1 mod 2; it is a sketch. The data agree with it: in all 142 stored records, the kernel vectors whose M is not divisible by 8 all have the same odd part modulo 4. (Kernel vectors also give M = 2, 6, 10 and 14, for elements of order 16, as §5.4 says; the statement is made precise in [note 25](25-checks-of-notes-20-23-24.md), finding 7.)

### 5.4 Results [E]

`sixlines_jacobian.py --character --pairs 3`, then `--character-extra` for the surfaces with fewer than two usable pairs (all their stored pairs) and for the controls. A surface is a configuration and its associate; the pairs are those stored by note 14.

| q | e₁ | G (note 19 §4.5) | surfaces | configurations | curves tried | computed | exponent above 16 | degenerate | not verified | surfaces with δ = +1, −1 | surfaces where all pairs agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 29 | 2 | Z/12 | 12 | 24 | 92 | 37 | 49 | 6 | 0 | 6, 6 | 12 of 12 |
| 37 | 2 | Z/12 | 12 | 24 | 82 | 39 | 41 | 1 | 1 | 6, 6 | 12 of 12 |
| 37 | 18 | Z/4 × Z/2 | 8 | 16 | 68 | 30 | 36 | 2 | 0 | 4, 4 | 8 of 8 |
| all | | | 32 | 64 | 242 | 106 | 126 | 9 | 1 | 16, 16 | 32 of 32 |

- **Every curve computed passed its own checks**: the point X has the stated order, πX = π₁X, J²X = −X, and the two pairings of X computed directly equal the values assembled from the matrices. In the one curve marked "not verified" the point was not in the kernel (the twelve classes did not span); it is not counted, and its value would have agreed with its surface.
- **The form of Z.** In 105 of the 106 curves every vector of the kernel lattice gives a Z that is a real multiple of conj(16/ϖ₁), as §5.3 requires. In one curve 9 vectors do not (their classes y_b did not span, so the lattice was too large); the vector used there passed the direct checks. The convention is fixed by the data: with N/ϖ₁ in place of its conjugate there is no solution.
- **M** is 12 in 71 computed curves and 4 in 60 (the element of the kernel has order 8), and 10, 14 or 2 in 11 (order 16); these counts include the controls. So M = 4·(odd) or 2·(odd), and u is that odd number modulo 4: 3 on 56 of the 106 curves and 1 on 50.
- **u depends on the pair, δ does not.** u is the same on all computed pairs of a surface for 10 of the 32 surfaces. On the 29 surfaces of the first run, u, u·sign(Re c₁·Im c₁) and u·sign(c₂) are constant on 11, 12 and 15; δ on 29. After the second run δ is constant on 32 of 32, with 2 to 6 pairs each (106 pairs; if the values were independent coin tosses, agreement on all surfaces would have probability 2⁻⁷⁴).
- **Balance, and a prediction.** After the first run the surfaces with δ = +1 and −1 numbered 5 and 6, 5 and 6, 3 and 4, with one surface of each torsor missing. A character takes each value on half of the torsor, so each of the three missing surfaces had to have δ = +1. The second run gave +1 on all their computed pairs: 2, 6 and 5 pairs.
- **Not one of the candidates of note 19.** Against the conic class, j2, j4, the quartic product and the orientation of `sixlines_characters.jsonl`, configuration by configuration, δ agrees and disagrees in mixed proportions on each torsor (for instance the conic class: 5 and 16, 15 and 5, 8 and 2); and it is not the flag that says whether the configuration or its twist has the polynomial P.

### 5.5 Attempts to break it [E]

| attempt | outcome |
|---|---|
| other pairs of lines of the same configuration, and of the associated configuration | same δ: 32 surfaces of 32, 106 pairs |
| the twists C₁ and C₂ in place of C₀ (the first usable pair of the first four surfaces of each torsor) | same u: 24 of 24 |
| other random classes y_a (another seed), same curves | same u: 12 of 12 |
| u without the normalisation, or with half of it | fails: constant on 10 to 15 surfaces only |
| the other identification of i with J in π₁ | would not give a kernel: χ fails with it (§3) |
| the other convention for Z | no solution M on any curve |
| a second prime | both primes available, 29 and 37, give the same picture |
| the three surfaces left out by the first run | the values that balance required |

What was **not** tried here, and could still break it:

- **A torsor with q ≡ 1 mod 8.** Both primes are 5 mod 8. The five torsors over F₄₉ of note 19 §4.5 are 1 mod 8, and there Frobenius conjugation acts with χ₋₁ = −1, which would be a test with known g; the script handles prime fields only. [Note 26](26-delta-is-the-character-at-q-41.md) tries q = 41 ≡ 1 mod 8, and δ is not broken: on the two torsors with the character there it is the same on all computed pairs of each of 28 surfaces (250 pairs) and is the character χ₋₁. The code of `character_pair` accepts a kernel vector only when M is unique modulo N, which cannot happen there (N = 32, M known modulo 16); the reading is redone in `sixlines_delta41.py`. The torsors over F₄₉ are not done.
- **The pairs with conjugate branch points and exponent above 16**, 126 of 242. If δ on those disagreed with the rest, it would depend on the pair. They need the pairing of order 32, that is F_{p⁸}.
- **Equivariance.** Balance is necessary, not sufficient. A function that changes sign under some fixed-point-free involution of the torsor is balanced too.

### 5.6 What it would mean [H]

If δ(g·X) = χ₋₁(g)·δ(X): given two surfaces of a torsor, the value of χ₋₁ on the class that relates them is δ(X)·δ(Y). The computation needs the six lines, one pair, the quartic χ of the curve, and arithmetic in Jac(C)(F_{p⁴}); apart from χ, which note 14 gets from point counts of C, nothing in it grows faster than a power of log p. That would break the decision problem of note 13 §6 for the polynomials where χ₋₁ is a character of G, by one bit, as the character δ does in the paper of §5.1 for class-group actions with even class number, such as those of supersingular curves over F_p with p ≡ 1 mod 4 (that paper's attack does not affect CSIDH, whose class group has no 2-torsion; [note 25](25-checks-of-notes-20-23-24.md), finding 5). It would not touch the 397 of 477 polynomials where 2 is totally ramified in K⁺, nor any q ≡ 3 mod 4, and it would say nothing about the other characters of note 19 §2.5, which have no such description.

It also bears on note 19 §2.2: there the K3 lattices of a torsor form one genus, so no local invariant of the surface lattice can be a character. δ is not local on that lattice: the second exterior power forgets the scaling of the hermitian form by a unit (it sees its square), and δ reads exactly that scaling at the prime above 2.

---

## 6. Check and limits

- **Checks.** One implementation. The checks made here are internal: §3, the dimension assertions, the direct verification of each point of §5, and the agreement of independent pairs. [Note 25](25-checks-of-notes-20-23-24.md) checks §§3 to 5 independently (findings 4 and 7 to 10).
- **Fitted choices: none after the fact, two fixed by a failing alternative.** The identification i ↦ J³ and the convention Z = M·conj(N/ϖ₁) were each decided by the other choice failing outright (§3, §5.4), on the first curve, before any comparison between surfaces. The formula for δ was written down before the first comparison of two pairs, and was not changed.
- **Prime fields with p ≡ 1 mod 4 only**: 13, 17, 29, 37 for §3; 29 and 37 for §4 and §5. Curves where Frobenius inverts J, p ≡ 3 mod 4, and the surfaces over F₂₅ and F₄₉ are not handled.
- **Pairs.** Only pairs stored by note 14 (at most 8 of the 30 per configuration for p ≥ 29). 2 pairs in §4, 5 runs in §3 and 9 curves in §5 are refused by the setup as degenerate although note 14 computes them (no base point z, as far as I looked; not examined further).
- **The level test of §4 covers three of six cases**, 29 of the 49 surfaces. In the three it is a count on every surface, not a sample. Which surface is on which level is not known on the lattice side, so the agreement is of numbers (and of the self-associated configurations); it could in principle be a coincidence of counts 1/6, 5/14, 1/2 against 1/6, 5/14, 1/2.
- **The rule of §4.1** assumes the lattice form of the correspondence of note 15 §2 at ℓ, which is what is being tested.
- **§5: the lattice statements are sketches.** That u is multiplied by ℓ under a step (tested on the lattices of note 20 in [note 25](25-checks-of-notes-20-23-24.md), finding 4; see §1 item 8), that the complement of V₁ does not contribute, and the relation between the Tate and Weil pairings (recalled) are not proved here. What is measured is: a number modulo 4 attached to a curve by a procedure that verifies its own point, and its behaviour over pairs and surfaces.

- **Random classes.** Sums of three points, or traces from F_{p^{2r}}. That they generate is checked where it matters: in §4 the subgroup found has the full order, in §5 the point is verified directly.
- **Sources recalled.** The method of §2.2 and the evaluation by determinants of §2.3 were implemented from memory of the papers named below, not from the texts.

---

## 7. What this does not give

- A proof of anything. Every statement of §§3 to 5 is a count on the curves listed.
- The equivariance δ(g·X) = χ₋₁(g)·δ(X). So this note does not show that a character of G is computable, or that the decision problem is broken; [note 24](24-the-form-mod-4-is-the-character.md) establishes the equivariance on the three torsors.
- The level test for (29, 18, ℓ = 3), (37, 14, ℓ = 11), (37, 38, ℓ = 7), and any test at ℓ = 2, where the Prym varieties of different pairs differ and note 20 finds that the lattices of the surfaces are not those of the maximal order.
- The endomorphism ring of the Prym variety itself; only the number of Jordan blocks at one maximal ideal.
- A Weil pairing as such (only the Tate pairing), and no pairing of order above 16.
- Arithmetic for p ≡ 3 mod 4, for the twisted Frobenius of note 14 Lemma 4, or over F_{p²} as a base field.
- A step of the action on either side, or the reverse map from a fourfold to six lines.
- Anything about the other quadratic characters of G, or about polynomials where χ₋₁ is not a character.

**The two computations that would settle item 8 of §1.** (a) Labels for one of the three torsors by the CM method of [note 20](20-cm-method-and-the-action.md), which skipped them at 12 lattices per type; with labels, δ of this note (`sixlines_jacobian_character.jsonl` and `_extra.jsonl`, field `delta`, by configuration) is tested against χ₋₁ directly. (b) The same invariant over F₄₉, where Frobenius conjugation is a known element with χ₋₁ = −1. (a) is done in [note 24](24-the-form-mod-4-is-the-character.md), for all three torsors: δ is constant on the two cosets of ker χ₋₁ and opposite between them, on 32 surfaces of 32. (b) is not done.

## Reproduce

From `experiments/`:

```bash
python -I sixlines_jacobian.py --validate-set --workers 4 --out sixlines_jacobian_validate.jsonl > sixlines_jacobian_validate.log
python -I sixlines_jacobian.py --lattice > sixlines_jacobian_lattice.log
python -I sixlines_jacobian.py --levels --workers 4 --out sixlines_jacobian_levels.jsonl > sixlines_jacobian_levels.log
python -I sixlines_jacobian.py --character --pairs 3 --workers 4 --out sixlines_jacobian_character.jsonl > sixlines_jacobian_character.log
python -I sixlines_jacobian.py --character-extra sixlines_jacobian_character.jsonl --workers 4 --out sixlines_jacobian_character_extra.jsonl > sixlines_jacobian_character_extra.log
```

About 5 minutes, 1 minute, 7 minutes, 50 minutes and 35 minutes with 4 processes on a loaded machine. The last command was run in two pieces (the first torsor, then `--torsor 37 2 --torsor 37 18` after a fault in the driver was repaired) and the outputs were joined.

One curve: `--validate p b c e f i j [--tw k] [--r r]` (all checks of §3 on the configuration (b, c, e, f) and the ordered pair (i, j)); `--one p b c e f i j [--tw k]` (the computation of §5.3 on one curve). The script needs `sixlines_prym.py`, `sixlines_prym.jsonl` and the files `sixlines_p{p}.jsonl` of note 2; `--lattice` also needs `sixlines_cm.py` and `sixlines_vertical.py`. The classes are random (PARI's generator, seeded), so the counts of generators tried can differ between runs; the group structures, the values of M and all verdicts do not depend on the random choices. No package was installed.

## Sources

- Castryck, Sotáková, Vercauteren, *Breaking the decisional Diffie–Hellman problem for class group actions using genus theory*, CRYPTO 2020, https://eprint.iacr.org/2020/151 : §2.5, §3 (Theorems 7 and 8, §3.1, §3.2) and Theorem 10 read in the copy of the authors' institute [P].
- Recalled, not re-read [M]: Khuri-Makdisi, *Linear algebra algorithms for divisors on an algebraic curve*, Math. Comp. 73 (2004), and *Asymptotically fast group operations on Jacobians of general curves*, Math. Comp. 76 (2007), for the representation of divisor classes by spaces of sections and the two operations on them; Bruin, *The Tate pairing for abelian varieties over finite fields*, J. Théor. Nombres Bordeaux 23 (2011), for the evaluation of a function on a divisor by determinants; the relation T_n(P, Q) = e_n(P, (π − 1)R) with nR = Q between the Tate and Weil pairings; that the kernel of the norm of a ramified double cover is connected (Mumford, *Prym varieties I*); Jacobowitz on hermitian lattices over ramified dyadic rings. None of these enters except through computations that are checked in §3 and §5.
- In this repository: notes [2](02-six-line-count.md), [7](07-ring-changing-steps.md), [8](08-below-the-top.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [18](18-checks-of-notes-13-to-16.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [24](24-the-form-mod-4-is-the-character.md), [25](25-checks-of-notes-20-23-24.md), [26](26-delta-is-the-character-at-q-41.md).
