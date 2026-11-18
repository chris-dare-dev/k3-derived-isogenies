# 11. Multiplicity 6 with ℓ⁴ dividing P(±q)

*Treats multiplicity 6 with ℓ⁴ dividing P(±q), the case for which [note 10](10-fourth-power.md), which describes the rational steps at eigenvalue ±q when ℓ⁴ divides P(±q) in multiplicity 2 and 4, gives computed examples only; builds on notes 8 to 10. Labels: **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_rank6.py`. Not re-checked.*

---

## 1. Result

Notation is that of notes 8, 9 and 10: ℓ odd, prime to q and to disc T; the eigenvalue is q (for −q replace F by −F); V is the part of T ⊗ Q_ℓ on which F ≡ q; τ = (F − q)/(F + q), g = F/q; 𝓛 is the set of unimodular τ-stable lattices in V. Here V has rank 6 and ℓ⁴ is the exact power of ℓ dividing P(q). No pair of the six-line data is of this kind; everything measured is on searched polynomials.

1. **Four patterns (Lemma 7).** The valuations of F − q on the three pairs of eigenvalues, with v(ℓ) = 1, are (1, ½, ½), (½, ¾, ¾), (3/2, ¼, ¼) or (⅔, ⅔, ⅔), in twelve configurations of primes.
2. **Three kinds of lattices (Proposition U).** By the Jordan blocks of τ on L/ℓL:
   - *leaves*, blocks (1, 5): exactly one step;
   - *middle* lattices, blocks (3, 3): ℓ + 1 steps;
   - *wide* lattices, blocks (1, 1, 1, 3) or (1, 1, 2, 2), where L/τL ≅ (Z/ℓ)⁴: ℓ² + ℓ + 1 steps, respectively ℓ + 1 or 2ℓ² + ℓ + 1.

   L is wide exactly when ℓ/τ preserves it, and is not a leaf exactly when τ³/ℓ preserves it. A leaf hangs on a lattice that is not a leaf, or is joined to one other leaf and to nothing else.
3. **The wide lattices are a case already solved (Theorem R).** They are the unimodular lattices stable under ψ = ℓ/τ, with the same steps, and ψ satisfies the hypothesis of notes 8 and 9 (its determinant is ℓ²). So they form the volcanoes of those parts: none in the pattern (3/2, ¼, ¼); two trees joined at their roots in (½, ¾, ¾); a volcano of depth 1 in (⅔, ⅔, ⅔); and in (1, ½, ½) the product of the lattices of a rank-2 part with a rank-4 volcano, a step moving in one factor only.
4. **In the pattern (3/2, ¼, ¼) there are exactly two lattices, exchanged by one step (Corollary 1).** This is 194 of the 235 polynomials found.
5. **A description by a quartic order (Theorem S).** When −P(q)P(−q) is a square in Q_ℓ, an explicit polynomial f of degree 4 over Z_ℓ is computed from P. The lattices of 𝓛 are the lattices in B = Q_ℓ[x]/f stable under x, up to scalars and of one parity; a step is a pair ℓΛ ⊂ Λ′ ⊂ Λ of index ℓ². Leaves are the lattices xS for S = Z_ℓ[x]/f; the others are the lattices stable under (x − λ₀)³/ℓ; the wide ones those stable under (x − λ₀)²/ℓ. When −P(q)P(−q) is not a square the same holds after the unramified quadratic extension, and 𝓛 is the part fixed by its Galois group.
6. **Check [E].** For q = 5, ℓ = 3, on 235 polynomials and 18,119 lattices: items 1, 2, 4 and the tree statements of item 3 hold wherever they could be tested; in all 102 cases where −P(q)P(−q) is a square and a starting lattice was found, the graph of lattices and the graph of the quartic order have the same invariant. Details and limits are in §7.

What this is not: a closed formula for the graph in each of the twelve configurations. Outside the pattern (3/2, ¼, ¼) the graph is given by Theorem R for its wide part and by Theorem S as a finite computation with 4 × 4 matrices; I have not turned the latter into a formula. §8 lists what is missing.

---

## 2. Patterns

Let x_j^{±1}, j = 1, 2, 3, be the eigenvalues of g on V in an algebraic closure, and w_j = v(x_j − 1) with v(ℓ) = 1; it is also the valuation of τ and of F − q at x_j. A prime 𝔩 of E above ℓ with ramification index e, residue degree f and v_𝔩(F − q) = v carries e·f eigenvalues, all with w = v/e.

**Lemma 7.** [D] w₁ + w₂ + w₃ = 2, and up to order (w₁, w₂, w₃) is one of

| pattern | primes with F ≡ q, as (e, f, v) |
|---|---|
| (1, ½, ½) | one of: a pair (1, 1, 1); an inert (1, 2, 1). With one of: a pair (2, 1, 1); a prime (2, 2, 1), inert or ramified over E₀; two primes (2, 1, 1) ramified over E₀ |
| (½, ¾, ¾) | (2, 1, 1) and (4, 1, 3), both ramified over E₀ |
| (3/2, ¼, ¼) | (2, 1, 3) and (4, 1, 1), both ramified over E₀ |
| (⅔, ⅔, ⅔) | a pair (3, 1, 2); or one prime (3, 2, 2), inert over E₀ |

"Pair" means two primes 𝔩 ≠ 𝔩̄; the others are their own conjugates.

*Proof.* Σ_𝔩 f·v = 4 is the sum of w over the six eigenvalues, and x, x⁻¹ have the same w. A prime equal to its conjugate carries e·f/2 pairs {x, x⁻¹}, a pair of primes e·f of them. A prime with f odd equal to its conjugate is ramified over E₀, and at a prime ramified over E₀ e is even and v is odd (τ is skew). Listing the ways to cover three pairs: one prime or pair covering all three forces w = ⅔ with (e, f) = (3, 1) for a pair or (3, 2) for one prime, the cases (6, 1), (2, 3), (1, 3), (1, 6) giving v = 4 at a prime ramified over E₀ or a non-integral v. A cover by one pair {x, x⁻¹} and two others with a common w_b < 1 gives w_b = ½, ¼ or ¾ and w_a = 1, 3/2, ½; w_b = ½ comes from e = 2 and v = 1, w_b = ¼ and ¾ from (4, 1) ramified over E₀ with v = 1 and 3; w_a = 1 from e = 1 and v = 1 (a prime (2, 1) would need v = 2), w_a = 3/2 and ½ from (2, 1) with v = 3 and 1. Three separate covers give (1, ½, ½) with two primes (2, 1, 1). ∎

In the search of §7: (3/2, ¼, ¼) 194 polynomials, (½, ¾, ¾) 13, (1, ½, ½) 19, (⅔, ⅔, ⅔) 9.

---

## 3. Three kinds of lattices

**Proposition U.** [D] Let L be in 𝓛.

1. The sizes of the Jordan blocks of τ on L/ℓL are (1, 5), (3, 3), (1, 1, 1, 3) or (1, 1, 2, 2). In the first two L/τL ≅ (Z/ℓ²)²; in the last two L/τL ≅ (Z/ℓ)⁴.
2. The numbers of rational steps are: 1 for (1, 5); ℓ + 1 for (3, 3); ℓ² + ℓ + 1 for (1, 1, 1, 3); for (1, 1, 2, 2), 2ℓ² + ℓ + 1 if −P(q)P(−q) is a square in Q_ℓ and ℓ + 1 if not.
3. L is wide, that is L/τL ≅ (Z/ℓ)⁴, exactly when ψ = ℓ/τ maps L into L. L is not a leaf, that is its blocks are not (1, 5), exactly when τ³/ℓ maps L into L.
4. A leaf has one neighbour. If the neighbour is a leaf, the two form a connected component of the graph.

*Proof.* (1) By note 10, Proposition O, the number μ of blocks is 2 or 4, with L/τL of the two shapes, and blocks of even size come in pairs; the partitions of 6 into 2 or 4 parts with this property are the four listed. (2) is the formula of Proposition O, with (ν, ρ) = (1, 1), (0, 2), (3, 1), (2, 2). In the last case the form on W(L) modulo its radical has rank 2, and L/ℓL is the orthogonal sum of that plane and of a hyperbolic space of rank 4 (the radical and a dual to it); so the plane is hyperbolic exactly when L/ℓL is, that is when det V = −1 up to squares, and det V = P(q)P(−q) up to squares (§5, Step 0). (3) L/τL is killed by ℓ exactly when ℓL ⊂ τL; and the blocks have size ≤ 3 exactly when τ³L ⊂ ℓL. (4) follows from (2). ∎

So the graph consists of the *core*, the lattices that are not leaves, and of leaves attached to it, unless it is a single edge.

---

## 4. The wide lattices

**Theorem R.** [D] The wide lattices are exactly the unimodular lattices of V stable under ψ = ℓ/τ, and two of them are joined by a rational step exactly when they are ℓ-neighbours. The element ψ is skew and det ψ = ℓ² up to a unit. Consequently:

- pattern (3/2, ¼, ¼): ψ is not integral, and there is no wide lattice;
- pattern (½, ¾, ¾): ψ has valuation 1 at both primes, which are ramified over E₀. By note 9, Theorem N, applied to ψ: the wide lattices form two rooted trees of depth n − 1 joined at their roots, where ℓ^{2n} is the index of Z_ℓ[ψ] in the maximal order; 2ℓ^j of them of type j;
- pattern (⅔, ⅔, ⅔): ψ has valuation 1 at primes with e = 3. By note 8, Theorem M, applied to ψ, with n = 1: for one prime inert over E₀, one top and ℓ + 1 lattices below it, or no wide lattice; for a pair, a line of tops, each with ℓ − 1 lattices below it;
- pattern (1, ½, ½): V = V₁ ⊥ V₂ with V₁ of rank 2 (the part with w = 1) and V₂ of rank 4. Every wide lattice is L₁ ⊕ L₂ with L₁ a unimodular lattice of V₁ stable under the maximal order (one lattice, a line of translates, or none: note 10, Proposition Q) and L₂ a unimodular τ-stable lattice of V₂, where ℓ² exactly divides the determinant of τ (notes 8 and 9, or note 10, Theorem P). A step between wide lattices changes exactly one of L₁, L₂, by a step of its own graph.

*Proof.* If L is wide, ψL ⊂ L. Conversely let L be unimodular with ψL ⊂ L. The form ⟨ψ⁻¹x, y⟩ mod Z_ℓ on L/ψL is alternating and perfect, as in note 8, Lemma 1, so L/ψL ≅ A × A; its order is ℓ^{6−4} = ℓ², so L/ψL ≅ (Z/ℓ)² and ℓL ⊂ ψL, that is τL ⊂ L. So L is in 𝓛 and is wide. A rational step between two lattices of 𝓛 is the relation of being ℓ-neighbours (on V all eigenvalues of F are ≡ q), which does not involve τ or ψ.

The valuations of ψ are 1 − w_j. In (3/2, ¼, ¼) one is negative. In (½, ¾, ¾) and (⅔, ⅔, ⅔) they are positive, ψ is topologically nilpotent, skew, generates A, and has determinant ℓ²; notes 8 and 9 use nothing else about τ, and the configurations are (d), and (a) or (b) with e = 3, of note 8, Lemma 3. In (a), if no lattice is stable under the maximal order, the argument of note 9, Step 1, gives a line of the residue field F_{ℓ²} equal to its orthogonal for the norm form, which does not exist; so there is no wide lattice at all. In (1, ½, ½) the characteristic polynomial of ψ is, modulo ℓ, a product of two coprime factors, so Z_ℓ[ψ] contains the idempotents of V₁ and V₂ and L = L₁ ⊕ L₂. L₁ is stable under ψ₁ and under ℓ/ψ₁ = τ₁, so Proposition Q applies; L₂ is stable under ψ₂, of determinant ℓ² on V₂, hence under τ₂ = ℓ/ψ₂ by the argument above, and conversely. Two sums L₁ ⊕ L₂ and L₁′ ⊕ L₂′ meet in index ℓ exactly when one summand is the same and the other meets in index ℓ. ∎

In (1, ½, ½) the wide part is a product in which a step moves one coordinate, where in rank 4 (note 10) a step moves both.

---

## 5. The quartic order

### 5.1 Statement

Write P(t) = t⁶ − e₁t⁵ + e₂t⁴ − e₃t³ + q²e₂t² − q⁴e₁t + q⁶.

**Theorem S.** [D] Suppose −P(q)P(−q) is a square in Q_ℓ. Let κ be 1 if qP(−q) is a square in Q_ℓ and a fixed non-residue otherwise, and

> s₊² = κ³P(−q)/q³,  s₋² = −κ³P(q)/q³,  a = (s₊ + s₋)/(2qκ),  d = (s₊ − s₋)/2,
>
> f(x) = x⁴ − a x³ + κe₁ x² − d x + q²κ².

Then f has coefficients in Z_ℓ and f ≡ (x − λ₀)⁴ mod ℓ. Put B = Q_ℓ[x]/f, S = Z_ℓ[x]/f and ε = x − λ₀. For every form for which 𝓛 is not empty:

1. 𝓛 is in bijection with the classes up to scalars of the lattices Λ in B stable under S, of one of the two parities of v(det Λ); the parity depends on the form.
2. The lattices of Λ and Λ′ are joined by a rational step exactly when, after scaling, ℓΛ ⊂ Λ′ ⊂ Λ with Λ/Λ′ of order ℓ².
3. The Jordan blocks of ε on Λ/ℓΛ and of τ on L/ℓL correspond: (4) to (1, 5); (3, 1) to (3, 3); (2, 2) to (1, 1, 1, 3); (2, 1, 1) to (1, 1, 2, 2).
4. The leaves are the lattices xS, x ∈ B^×. The lattices that are not leaves are those stable under ε³/ℓ, and the wide ones those stable under ε²/ℓ.
5. B^× acts on the graph of all classes of S-lattices, of both parities, and is transitive on the leaves.

Changing the signs of s₊, s₋ replaces x by −x or by qκ/x and does not change S.

**The other case.** [D] If −P(q)P(−q) is not a square in Q_ℓ, it is one in the unramified quadratic extension k, and the theorem holds for V ⊗ k over k, with O_k-lattices and index ℓ⁴. A unimodular lattice L of V gives L ⊗ O_k, and this identifies 𝓛 with the lattices of 𝓛(V ⊗ k) fixed by Gal(k/Q_ℓ), and the steps with the steps between fixed lattices. The action of Gal(k/Q_ℓ) is made explicit in [note 12](12-rank-six-galois-action.md): the order over O_k comes from an order Z_ℓ[X]/f₀ over Z_ℓ, the conjugation sends the class of a lattice to the class of its dual for a hermitian form, and 𝓛 is the set of the stable lattices that are self-dual, or self-dual up to ℓ. This case is tested there.

### 5.2 Consequences

**Corollary 1.** [D] In the pattern (3/2, ¼, ¼), 𝓛 is empty or consists of two lattices joined by a step. Both are leaves.

*Proof.* Let k′ be Q_ℓ or k as above, and μ_i the eigenvalues of x divided by a square root of qκ, so that the eigenvalues of g are the μ_iμ_j and μ₁² = (μ₁μ₂)(μ₁μ₃)(μ₁μ₄). Of the three factors one is ≡ 1 to valuation 3/2 and the other two to valuation ¼, and these two are conjugate under the tame group of a totally ramified quartic extension without being inverse to each other, so that their leading terms differ by a factor ±i; hence v(μ₁² − 1) = ¼. So B ⊗ k′ contains an element of valuation ¼ and is a totally ramified quartic field with uniformiser ε, and S ⊗ O_{k′} is its maximal order. Its lattices are the powers of the maximal ideal, four classes up to scalars, two of each parity, and the two of one parity are joined. So V ⊗ k′ has exactly two unimodular lattices fixed by g. If k′ = k and one of them is fixed by the Galois group, so is its only neighbour. ∎

**Corollary 2.** [D] In all patterns but (3/2, ¼, ¼), every leaf is joined to a lattice of the core.

*Proof.* Suppose first that −P(q)P(−q) is a square. By Theorem S (5) it is enough to look at the leaf S. Its neighbour is ε²S + ℓS. Write ε⁴ = −(c₃ε³ + c₂ε² + c₁ε + c₀) with ℓ dividing every c_i. Then (ε³/ℓ)·ε² = −((c₃/ℓ)ε⁴ + (c₂/ℓ)ε³ + (c₁/ℓ)ε² + (c₀/ℓ)ε) lies in ε²S + ℓS exactly when ℓ² divides c₀, and (ε³/ℓ)·ℓ lies in it. Now c₀ is the product of the four conjugates of ε, and v(c₀) = Σ v(μ_i − 1). In the three patterns this sum is at least 2: each μ_i² − 1 has valuation ≥ ½, as a product of three numbers ≡ 1 to valuation ≥ ½. If −P(q)P(−q) is not a square, the same over k shows that the neighbour of L ⊗ O_k is not a leaf, and it is L′ ⊗ O_k for the neighbour L′ of L; the Jordan blocks do not change with the extension. ∎

**Number of steps.** In the case of a square, the steps from Λ are the planes of Λ/ℓΛ stable under ε: 1, ℓ + 1, ℓ² + ℓ + 1 and 2ℓ² + ℓ + 1 for the four block types, in agreement with Proposition U.

### 5.3 What the graphs look like

From the run of §7, for ℓ = 3. "Square" refers to −P(q)P(−q). When no prime is a pair the whole connected component is listed, otherwise a ball.

| pattern | primes (e, f, v) | −P(q)P(−q) | polynomials | lattices: wide / middle / leaves |
|---|---|---|---|---|
| (3/2, ¼, ¼) | (2, 1, 3), (4, 1, 1) | square | 106 | 2, all: 0 / 0 / 2 |
| (3/2, ¼, ¼) | (2, 1, 3), (4, 1, 1) | not a square | 88 | 2, all: 0 / 0 / 2 |
| (½, ¾, ¾) | (2, 1, 1), (4, 1, 3) | square | 8 | 26, all: 2 / 6 / 18 |
| (½, ¾, ¾) | (2, 1, 1), (4, 1, 3) | not a square | 5 | 26, all: 2 / 6 / 18 |
| (1, ½, ½) | (1, 2, 1), (2, 1, 1), (2, 1, 1) | square | 2 | 26, all: 2 / 0 / 24; or 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (1, 2, 1), (2, 2, 1) | square | 5 | 26, all: 2 / 0 / 24; or 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (1, 2, 1), (2, 2, 1) | not a square | 3 | 53, all: 5 / 0 / 48; or 26, all: 2 / 0 / 24; or 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (2, 2, 1), pair (1, 1, 1) | square | 2 | ball of 269: 37 / 112 / 120 |
| (1, ½, ½) | (2, 2, 1), pair (1, 1, 1) | not a square | 4 | ball of 116: 16 / 28 / 72; or of 98: 14 / 24 / 60 |
| (1, ½, ½) | pair (1, 1, 1), pair (2, 1, 1) | square | 2 | ball of 575: 91 / 328 / 156; or of 404: 68 / 240 / 96 |
| (1, ½, ½) | (2, 1, 1), (2, 1, 1), pair (1, 1, 1) | not a square | 1 | no starting lattice found |
| (⅔, ⅔, ⅔) | (3, 2, 2) | not a square | 4 | 44, all: 5 / 3 / 36 |
| (⅔, ⅔, ⅔) | pair (3, 1, 2) | square | 5 | ball of 224: 23 / 111 / 90 |

Where two or three graphs are given they belong to different forms of the same polynomial. The balls have radius 4 around the starting lattice.

Three of these in words; in each, every leaf hangs on a wide lattice with ℓ² + ℓ + 1 = 13 steps.

- (½, ¾, ¾), 26 lattices: the two wide lattices are joined (Theorem R with n = 1); each has three middle lattices and nine leaves; each middle lattice is joined to one wide lattice and three middle ones. So the eight lattices of the core form the complete bipartite graph on 4 + 4 vertices.
- (1, ½, ½) with the primes (1, 2, 1) and (2, 2, 1) or two (2, 1, 1), 26 lattices: two wide lattices joined, each with twelve leaves, and no middle lattice. With 53 lattices: five wide lattices forming a star, whose centre has no other step, the four others having twelve leaves each. With 8 lattices: eight middle lattices forming the complete bipartite graph on 4 + 4 vertices; this occurs for the forms for which the rank-2 part has no unimodular lattice, so that there is no wide lattice.
- (⅔, ⅔, ⅔) with one prime, 44 lattices: five wide lattices forming a star (Theorem R), whose centre has ℓ + 1 = 4 steps and no other; each of the four others has also nine leaves and three middle lattices; each of the three middle lattices is joined to those four.

---

## 6. Proof of Theorem S

**Step 0: the quadratic space.** As in note 10, Step 0, with θ = g − g⁻¹: det V = N_{A₀/Q_ℓ}(−4α²θ²) = −Π(a_j² − 4) up to squares, a_j = x_j + x_j⁻¹, because A₀ now has degree 3. And Π(4 − a_j²) = det(g − 1)·det(g + 1) = P(q)P(−q)/q¹². A unimodular lattice of rank 6 over Z_ℓ is determined by its determinant [M]; the space Λ²(Q_ℓ⁴) with the form x ∧ y ∈ Λ⁴ ≅ Q_ℓ is the sum of three hyperbolic planes, of determinant −1. So V is isometric to it exactly when −P(q)P(−q) is a square.

**Step 1: lattices.** GL₄(Q_ℓ) acts on Λ²(Q_ℓ⁴) by Λ²h, multiplying the form by det h. The pairs (h, c) with c² = det h, acting by c⁻¹Λ²h, map onto SO(V) with kernel the (z, z²) [M]. For Λ = hZ_ℓ⁴ with v(det h) = 2j, L(Λ) = ℓ^{−j}Λ²Λ is unimodular and depends only on the class of Λ up to scalars. SO(V) is transitive on unimodular lattices, so all of them are of this form. If L(hZ_ℓ⁴) = L(Z_ℓ⁴), write h = k₁·diag(ℓ^{a_i})·k₂ with k₁, k₂ in GL₄(Z_ℓ); then a_i + a_j = j for all i < j, so the a_i are equal and h is a scalar times an element of GL₄(Z_ℓ). So the unimodular lattices correspond bijectively to the classes of lattices with v(det) even.

The lattice Λ′ spanned by e₁, e₂, ℓe₃, ℓe₄ gives ℓ⁻¹Λ²Λ′, which meets Λ²(Z_ℓ⁴) in index ℓ. GL₄(Z_ℓ) is transitive on the planes of F_ℓ⁴, and there are (ℓ² + 1)(ℓ² + ℓ + 1) of them, as many as isotropic lines in the hyperbolic space of rank 6 over F_ℓ. So the neighbours of L(Λ) are the L(Λ′) with ℓΛ ⊂ Λ′ ⊂ Λ of index ℓ².

**Step 2: fixed points.** g corresponds to (h, c). The eigenvalues of g are the λ_iλ_j/c, λ_i the eigenvalues of h; they are distinct, so the λ_i are, and B = Q_ℓ[h] is an étale algebra of degree 4, over which Q_ℓ⁴ is free of rank one. If g fixes L(Λ₀) then hΛ₀ = ℓ^kΛ₀, and replacing (h, c) by (ℓ^{−k}h, ℓ^{−2k}c) we may take k = 0. Then g fixes L(Λ) exactly when hΛ = Λ, that is when Λ is stable under Z_ℓ[h]. The lattices with v(det Λ) even relative to Z_ℓ⁴ are one parity class of the lattices in B; which one, relative to the lattice Z_ℓ[h] itself, depends on the position of that lattice, hence on the form. This gives (1) and (2), given f.

**Step 3: the polynomial.** Let x⁴ − a x³ + b x² − d x + c² be the characteristic polynomial of h. The polynomial with roots λ_iλ_j is x⁶ − b x⁵ + (ad − c²)x⁴ − (a²c² + d² − 2bc²)x³ + c²(ad − c²)x² − c⁴b x + c⁶, and P(t) = Π(t − qλ_iλ_j/c). Scaling h by a unit z multiplies c by z², so c = qκ with κ = 1 or the fixed non-residue. Comparing coefficients, e₁ = b/κ, e₂ = (ad − c²)/κ², e₃ = ((a² − 2b)c² + d²)/κ³, and substituting,

> P(−q) = q³(aqκ + d)²/κ³,  P(q) = −q³(aqκ − d)²/κ³.

So κ is determined by P(−q), and a, d are as stated, up to the signs of the two square roots. Changing both signs replaces h by −h, and changing one replaces h by c·h⁻¹; both generate the same order. All eigenvalues of g are ≡ 1, so λ_iλ_j ≡ c for all i < j, the λ_i are congruent to each other, and f ≡ (x − λ₀)⁴.

**Step 4: Jordan blocks.** On L/ℓL = Λ²(Λ/ℓΛ), g acts as Λ²(1 + N) with N the nilpotent ε/λ₀ on Λ/ℓΛ. For N of type (3, 1), (2, 2), (2, 1, 1) the decomposition of Λ² into Λ² of the blocks and tensor products of two blocks gives (3, 3), (3, 1, 1, 1), (2, 2, 1, 1), using that the tensor square of a block of size 2 is a block of size 3 and one of size 1 in odd characteristic. For one block of size 4, with Ne₁ = 0 and Ne_i = e_{i−1}, the operator D = Λ²(1 + N) − 1 has rank 4 and D² has rank 3 in every odd characteristic, by direct computation on the e_i ∧ e_j; so its blocks are (5, 1). This gives (3), and (4): one block of size 4 means that Λ/ℓΛ, hence Λ, is generated by one element over S; blocks of size ≤ 3, respectively ≤ 2, mean ε³Λ ⊂ ℓΛ, respectively ε²Λ ⊂ ℓΛ. (5): multiplication by an element of B^× maps S-lattices to S-lattices and preserves the relation of (2); it is transitive on the xS.

**Step 5: the other case.** Over k the same proof applies to V ⊗ k, whose determinant is −1 up to squares in k, with ℓ as uniformiser and residue field F_{ℓ²}. A unimodular O_k-lattice of V ⊗ k stable under the Galois group is L ⊗ O_k for a unique unimodular lattice L of V, since k is unramified, and L ∩ L′ has index ℓ in L exactly when (L ∩ L′) ⊗ O_k has O_k-colength one. ∎

---

## 7. Check [E]

`sixlines_rank6.py` runs over the polynomials t³h(t + q²/t), h a monic integer cubic with three real roots in (−2q, 2q), that are irreducible, have q as a root of multiplicity 6 mod ℓ with ℓ⁴ exactly dividing P(q), and no prime above ℓ ramified over E₀ at the eigenvalue −q. For each one it

- reads the pattern off the primes and asserts that it is one of the four (Lemma 7);
- for two forms in each class of valuations that have a unimodular lattice of the simplest kind, follows the steps from that lattice: the whole component if no prime is a pair, the ball of radius 4 otherwise; the local rule of note 10 is asserted on every lattice;
- asserts Proposition U: the block types are the four listed, L/τL has the stated shape, leaves have one step, and in the finite case two leaves are joined only when they are the whole graph;
- in the finite case, lists all lattices independently when the list is small enough and asserts that the steps reach all of them; asserts that the wide lattices form a tree (Theorem R); and in the pattern (3/2, ¼, ¼) asserts that there are exactly two lattices, with the list complete (Corollary 1);
- when −P(q)P(−q) is a square, computes f to 40 ℓ-adic digits, asserts f ≡ (x − λ₀)⁴, follows the steps between lattices stable under the companion matrix of f, and asserts: in the finite case, that the graph of 𝓛 has the same colour-refinement invariant as the graph of one of the two parities; otherwise, that the ball of radius 3 around a leaf, with degrees, has the same invariant as the ball around the lattice S (Theorem S).

For q = 5, ℓ = 3, eigenvalue +q: 235 polynomials, 130 with −P(q)P(−q) a square and 105 not. All assertions pass. A starting lattice was found for 188 of them (355 forms in all), and 18,119 lattices were met. In the pattern (3/2, ¼, ¼) the two lattices, with a complete list, were found for 245 forms. The comparison with the quartic order was made, and agrees, for all 102 polynomials of the square case that have a starting lattice (199 forms). For 47 polynomials no starting lattice of the simplest kind was found among the forms tried, and only the pattern was checked.

Limits: the comparison of graphs is by an invariant; in the infinite case it is one ball; when −P(q)P(−q) is not a square nothing is compared with a quartic order, and the statement about the extension k is not tested here (it is in [note 12](12-rank-six-galois-action.md)); the depth n of Theorem R is 1 in all examples, so the trees found are an edge or a star.

---

## 8. What this does not give

- A closed formula for the graph in each configuration outside (3/2, ¼, ¼). The middle lattices and the way leaves and middle lattices attach to the wide ones are given only through Theorem S, that is by a computation for each polynomial; the three descriptions of §5.3 are observations.
- The case where −P(q)P(−q) is not a square, beyond Proposition U, Theorem R and Corollary 1: the action of the Galois group on the lattices of the quartic order over k is not made explicit here. It is made explicit in [note 12](12-rank-six-galois-action.md), as the duality of a hermitian form (Theorem T there), with the consequences that leaves exist for one class of forms only and a formula for their number; the graph is still obtained by a computation for each polynomial, not by a closed formula.
- Multiplicity 6 with ℓ⁶ or more dividing P(±q). Theorem S and its proof do not use the exponent 4, but the patterns, Proposition U and Theorem R do.
- Ranks above 6, for other families.
- Which class of forms the transcendental lattice of a surface has; and surfaces at all: the polynomials are searched ones, not known to come from surfaces, and the six-line data has no such case.
- The realisation of the steps by twisted equivalences, as in note 10 §8.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_rank6.py --out sixlines_rank6.jsonl > sixlines_rank6.log
```

It takes about four minutes. An assertion fails if any statement of the list fails. Options: `--search q,l`, `--both-signs`, `--forms`, `--radius`.

## Sources

As in note [10](10-fourth-power.md). In addition, recalled and not re-read [M]: the isomorphism between SO of a rank-6 space of determinant −1 and SL₄ modulo ±1, through the second exterior power (Knus, Merkurjev, Rost, Tignol, *The Book of Involutions*, §15.D).
