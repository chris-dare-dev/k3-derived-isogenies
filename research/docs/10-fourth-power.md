# 10. Eigenvalue ±q when ℓ⁴ divides P(±q)

*Treats eigenvalue ±q when ℓ⁴ divides P(±q), the case that notes [8](08-below-the-top.md) and [9](09-ramified-over-e0.md), which assume that ℓ² is the exact power of ℓ dividing P(±q), leave open; builds on notes 6 to 9. Labels: **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_deep.py`. Not re-checked.*

---

## 1. Result

Notation is that of notes 6, 7, 8 and 9. ℓ is an odd prime not dividing q or disc T, the eigenvalue is q (for −q replace F by −F), V is the part of T ⊗ Q_ℓ on which F ≡ q, m its rank, and ℓ^c the exact power of ℓ dividing P(q). Notes 8 and 9 treat c = 2. Here c is arbitrary; c is always even when a lattice exists, so the next case is c = 4. No pair (polynomial, ℓ) of the six-line data has c ≥ 4, so everything measured here is on searched polynomials.

1. **Note 8 expected "products of volcanoes". That is right in rank 4 and wrong in rank 6.**
2. **A local rule for any c (Proposition O).** The number of rational steps from a lattice is read off the sizes of the Jordan blocks of F − q on L/ℓL. For c = 4 it takes the values 0, 1, 2, ℓ + 1 (as for c = 2) and also ℓ² + 1, (ℓ + 1)², ℓ² + ℓ + 1, 2ℓ² + ℓ + 1 and ℓ³ + ℓ² + ℓ + 1.
3. **The Brauer group is no longer constant along steps.** For c = 2, Br(X)[ℓ^∞] ≅ (Z/ℓ)² at every vertex. For c = 4 it is (Z/ℓ)⁴ at some vertices and (Z/ℓ²)² at others of the same graph.
4. **Rank 2: nothing changes (Proposition Q).** Every lattice is stable under the maximal order, whatever c is. There is one lattice and no step, or a line of translates. This is the most frequent way for ℓ⁴ to divide P(±q).
5. **Rank 4: a complete description for every c (Theorem P).** Two numbers t₁, t₂ are computed from the factor of P at ℓ. Either they lie in Q_ℓ, and the graph is the product of two ordinary volcanoes, of quadratic algebras K₁, K₂ and depths r₁, r₂ read off t₁ − 2 and t₂ − 2; or they are conjugate over the unramified quadratic extension of Q_ℓ, and the graph is one volcano with branching ℓ² in place of ℓ; or no lattice exists. In the product a step moves both coordinates at once. The case c = 2 is the case where one factor is a single edge, which is why notes 8 and 9 saw one volcano.
6. **Rank 6: no closed description.** The local rule holds, and the graphs are computed in each configuration found (§6). They are not products: a typical one has lattices with ℓ² + ℓ + 1 steps and contains cycles.
7. **Check [E].** For q = 5, ℓ = 3: Theorem P agrees with the lattices found on all 285 rank-4 polynomials with c ≥ 4 and on a sample of 117 with c = 2; Proposition O holds on all 14,803 lattices met, in ranks 4 and 6. Details and limits are in §7.

So for the six-line family, whose transcendental lattice has rank 6, the rational steps at eigenvalue ±q are now described except when the eigenvalue has multiplicity 6 and ℓ⁴ divides P(±q). That case is treated in [note 11](11-rank-six-fourth-power.md) for ℓ⁴ exactly: four patterns of valuations; in the most frequent one, two lattices exchanged by one step; in general three kinds of lattices, the wide ones forming a volcano of notes 8 and 9 for the element ℓ/τ, and the whole graph given by the lattices over an explicit quartic order. It is a reduction and not a closed formula in every configuration.

---

## 2. A local rule for any c

Let 𝓛 be the set of unimodular τ-stable lattices in V, with τ = (F − q)/(F + q) as in note 8. For L in 𝓛, τ is nilpotent on L/ℓL and skew for the form, which is non-degenerate on L/ℓL.

**Proposition O.** [D] Let L be in 𝓛, and let s₁, …, s_μ be the sizes of the Jordan blocks of τ on L/ℓL, so that s₁ + … + s_μ = m. Let ν be the number of blocks of size 1 and ρ = μ − ν.

1. L/τL ≅ A × A for a finite abelian group A of order ℓ^{c/2}; in particular c is even. It is dual to the part of Br(X)[ℓ^∞] on which Frobenius acts trivially (note 6, Proposition G).
2. μ is the number of generators of L/τL. It is even and 2 ≤ μ ≤ c.
3. Blocks of even size occur an even number of times.
4. The eigenspace W(L) = ker(τ on L/ℓL) has dimension μ. The radical of the form on W(L) has dimension ρ, and the form induced on the quotient is non-degenerate of rank ν.
5. The number of rational steps from L with eigenvalue q is

   > (ℓ^ρ − 1)/(ℓ − 1) + ℓ^ρ · i(ν),

   where i(ν) is the number of isotropic lines of a non-degenerate quadratic form of rank ν over F_ℓ: i(0) = i(1) = 0; i(2) = 0 or 2; i(3) = ℓ + 1; i(4) = ℓ² + 1 or (ℓ + 1)²; i(5) = ℓ³ + ℓ² + ℓ + 1.

*Proof.* (1) and the parity of μ are note 7, Proposition I (3): the form ⟨(g − 1)⁻¹x, y⟩ mod Z_ℓ on L/(g − 1)L, g = F/q, is alternating and perfect. (2) W(L) is the kernel of τ on L/ℓL, and its dimension is that of the cokernel L/(τL + ℓL). (3) is the standard constraint on nilpotent elements of an orthogonal Lie algebra in odd characteristic [M]. (4) The radical of the form on W = ker τ is W ∩ W^⊥ = W ∩ im τ, since τ is skew. In a Jordan basis, ker τ is spanned by the last vectors of the blocks and ker τ ∩ im τ by those of the blocks of size ≥ 2. (5) A line of W is isotropic when it lies in the radical, or when its image in the quotient is an isotropic line; each line of the quotient has ℓ^ρ lifts not in the radical. The steps are the isotropic eigenlines (note 7, Proposition H). ∎

For c = 2, μ = 2 and this is Proposition L of note 8. For c = 4:

| μ | L/τL | blocks of size 1 | steps |
|---|---|---|---|
| 2 | (Z/ℓ²)² | 0 | ℓ + 1 |
| 2 | (Z/ℓ²)² | 1 | 1 |
| 2 | (Z/ℓ²)² | 2 | 0 or 2 |
| 4 | (Z/ℓ)⁴ | 0 | ℓ³ + ℓ² + ℓ + 1 |
| 4 | (Z/ℓ)⁴ | 1 | ℓ² + ℓ + 1 |
| 4 | (Z/ℓ)⁴ | 2 | ℓ + 1 or 2ℓ² + ℓ + 1 |
| 4 | (Z/ℓ)⁴ | 3 | ℓ² + ℓ + 1 |
| 4 | (Z/ℓ)⁴ | 4 | ℓ² + 1 or (ℓ + 1)² |

μ = 4 exactly when ℓL ⊂ τL. Both values of μ occur in one graph (§4, §6), so the group structure of Br(X)[ℓ^∞] changes along rational steps although its order does not.

---

## 3. Rank 2

**Proposition Q.** [D] Let m = 2, so that A = Q_ℓ[F] on V is a quadratic étale algebra with A₀ = Q_ℓ. Then every lattice of 𝓛 is stable under the maximal order O of A, for every c:

- A unramified field: 𝓛 has one element and there is no step;
- A = Q_ℓ × Q_ℓ: 𝓛 is the line of translates 𝔩^j𝔩̄^{−j}L, and the steps are the two translations;
- A ramified: 𝓛 is empty (and c is odd).

*Proof.* Let L be a τ-stable lattice with End(L) = O_j = Z_ℓ + ℓ^jO; then L = xO_j [M]. The form is Tr(γ x z̄) with γ in Q_ℓ, and the trace dual of O_j is (ℓ^jθ)⁻¹O_j with θ = ω − ω̄, O = Z_ℓ[ω]. So L = L^∨ exactly when N(x)γℓ^jθ is a unit of O_j. It is a multiple sθ of θ with s in Q_ℓ. Its norm −s²θ² is a unit only if θ² is a unit, that is A is unramified over Q_ℓ, and s is a unit; and sθ = 2sω − s·Tr(ω) lies in O_j only if ℓ^j divides 2s. So j = 0. The O-lattices and their steps are those of note 7, Proposition J. ∎

Here τ = ℓ^{c/2} times a unit, Br(X)[ℓ^∞] ≅ (Z/ℓ^{c/2})² on this part, and no step changes the ring. In a separate count over the search range of §7 (q = 5, ℓ = 3, eigenvalue +q), 1,293 of the 1,793 polynomials with c = 4 have rank 2; the script skips them.

---

## 4. Rank 4

### 4.1 Volcanoes

For a quadratic étale algebra K over a local field k with residue field of Q elements, and r ≥ 0, let Vol(K, r) be the graph of homothety classes of lattices in K stable under the order O_k + ϖ^rO_K, two classes being joined when they have representatives Λ′ ⊂ Λ of index Q. The *level* of a class is the j ≤ r with End(Λ) = O_k + ϖ^jO_K. [M, Kohel for k = Q_ℓ]:

- level 0 is one vertex if K is an unramified field, a line (infinite in both directions) if K is split, two joined vertices if K is ramified;
- a vertex of level 0 has Q + 1 neighbours in all if r ≥ 1; those not of level 0 have level 1;
- a vertex of level 1 ≤ j < r has one neighbour of level j − 1 and Q of level j + 1; a vertex of level r ≥ 1 has one neighbour, of level r − 1.

Vertices have a *parity*, which changes along every edge (the tree is bipartite).

### 4.2 The numbers t₁, t₂

Let g = F/q on V. The characteristic polynomial of g + g⁻¹ = (F + V)/q on V is the square of a polynomial x² − Sx + Π with coefficients in Z_ℓ, congruent to (x − 2)² mod ℓ; it is the factor at ℓ, with roots ≡ 2, of the minimal polynomial of (F + V)/q. Put

> t₁, t₂ = the roots of z² − Πz + (S² − 2Π − 4),  D = Π² − 4(S² − 2Π − 4).

If x^{±1}, y^{±1} are the eigenvalues of g on V, then {t₁, t₂} = {xy + (xy)⁻¹, x/y + y/x}, and

> D = (t₁ − t₂)² = det(g − 1)·det(g + 1),  (t₁ − 2)(t₂ − 2) = S² − 4Π,  v_ℓ(D) = c.

### 4.3 Statement

**Theorem P.** [D] Let m = 4.

1. *No lattice.* If v_ℓ(D) is odd, 𝓛 is empty for every form.
2. *Product of two volcanoes.* Suppose D is a square in Q_ℓ, so t₁, t₂ ∈ Q_ℓ. For i = 1, 2 let w_i = v_ℓ(t_i − 2) ≥ 1 and
   - if w_i is odd: K_i ramified, r_i = (w_i − 1)/2;
   - if w_i is even: r_i = w_i/2, and K_i is split or the unramified field according as (t_i − 2)(t_i + 2)/ℓ^{w_i} is a square or not.

   Then for every form for which 𝓛 is not empty there is ε in {0, 1} such that 𝓛 is in bijection with the pairs (x₁, x₂), x_i a vertex of Vol(K_i, r_i) over Q_ℓ, whose parities differ by ε; and (x₁, x₂), (y₁, y₂) are joined by a rational step exactly when x₁, y₁ are neighbours and x₂, y₂ are neighbours. Changing ε gives an isomorphic graph unless K₁ and K₂ are both unramified fields; in that case both values occur, for different forms, and only ε = 0 has a lattice stable under the maximal order.
3. *One volcano over Q_{ℓ²}.* Suppose D is not a square and v_ℓ(D) is even. Let k be the unramified quadratic extension of Q_ℓ and w = min(v_ℓ(Π − 4), c/2), which is v_k(t₁ − 2). If w is odd let K be a ramified quadratic extension of k and r = (w − 1)/2. If w is even let r = w/2, and K = k × k or the unramified quadratic extension of k according as (S² − 4Π)/ℓ^{2w} is a square in Q_ℓ or not. Then for every form for which 𝓛 is not empty, 𝓛 with the rational steps is Vol(K, r) over k, with Q = ℓ².
4. c = 2·v_ℓ(t₁ − t₂) ≥ 2·min(w₁, w₂) in (2), with equality if w₁ ≠ w₂; and c ≥ 2w in (3).

**Consequences.**

- *c = 2.* Then min(w₁, w₂) = 1 in (2): one factor, say the second, is a ramified algebra with r₂ = 0, a single edge. The second coordinate is determined by the parity of the first, and the graph is Vol(K₁, r₁). This gives back notes 8 and 9 in rank 4: K₁ unramified with r₁ = n in configuration (a), split with r₁ = n in (b), ramified with r₁ = n − 1 in (c) and (d). In (3), c = 2 forces w = 1 and r = 0: two lattices joined by a step.
- *A product with both factors non-trivial needs c ≥ 4.* The smallest are unramified or split with r = 1 times ramified with r = 1 (w = 2 and 3), and two algebras with w₁ = w₂ = 2.
- *Number of steps.* It is the product of the numbers of neighbours of the two coordinates: (ℓ + 1)² when both have all their neighbours, ℓ + 1 when one of them is at its lowest level r_i ≥ 1, and one when both are.
- *Number of lattices* when no K_i is split. For K₁ unramified with r₁ = 1 and K₂ ramified with r₂ = 1: (ℓ + 1)(ℓ + 2), which is 20 for ℓ = 3. In (3) with K unramified and r = 1: ℓ² + 2, which is 11 for ℓ = 3.
- *In (2), when neither volcano is a single edge, the graph contains cycles of length 4*: from (x₁, x₂) move to (y₁, y₂), (x₁, z₂), (y₁′, y₂) and back. So it is not a tree, unlike all graphs of notes 8 and 9.

### 4.4 Which configurations of primes give what, for c = 4

From the search of §7 (q = 5, ℓ = 3, eigenvalue +q). The primes are those of E above ℓ with F ≡ q; e, f, v are the ramification index, the residue degree and the valuation of F − q.

| primes with F ≡ q | graph | polynomials |
|---|---|---|
| (e=2, f=1, v=1), (e=2, f=1, v=3) | product of two single edges: two lattices | 100 |
| (e=2, f=1, v=1), (e=2, f=1, v=3) | over Q_{ℓ²}, ramified, r = 0: two lattices | 91 |
| (e=1, f=4, v=1) | over Q_{ℓ²}, unramified, r = 1: ℓ² + 2 lattices | 31 |
| pair (e=1, f=2, v=1) | unramified r = 1 × split r = 1 | 22 |
| pair (e=1, f=2, v=1) | unramified r = 1 × split r = 2 | 1 |
| pair (e=2, f=1, v=2) | ramified r = 1 × split r = 1 | 10 |
| (e=2, f=2, v=2) | unramified r = 1 × ramified r = 1: (ℓ + 1)(ℓ + 2) lattices | 9 |
| (e=1, f=2, v=1), pair (e=1, f=1, v=1) | over Q_{ℓ²}, split, r = 1 | 1 |

"pair" is a pair 𝔩 ≠ 𝔩̄; the other primes are their own conjugates. So 191 of the 265 polynomials have two lattices and one step, as in note 9 with n = 1, although c = 4; the others have a graph that is not a tree or has branching ℓ².

---

## 5. Proof of Theorem P

**Step 0: the quadratic space.** Write θ = g − g⁻¹, so that A = A₀ ⊕ A₀θ and θ² = (g + g⁻¹)² − 4 ∈ A₀. The form is ⟨x, z⟩ = Tr_{A/Q_ℓ}(αx z̄) with α in A₀, and for x = x₀ + x₁θ it equals 2Tr_{A₀/Q_ℓ}(αx₀²) − 2Tr_{A₀/Q_ℓ}(αθ²x₁²). So det V = N(θ²) = (a² − 4)(b² − 4) = D up to squares, where a, b are the two values of g + g⁻¹. If a unimodular lattice exists, det V is a unit times a square: this gives (1). A quadratic space of rank 4 over Q_ℓ with a unimodular lattice is determined by its determinant [M]: it is (M₂(Q_ℓ), det) if the determinant is a square, and (H₂(k), det), the hermitian 2 × 2 matrices over k = Q_ℓ(√D), if not.

**Step 1: lattices and trees, D a square.** GL₂ × GL₂ acts on M₂ by x ↦ h₁xh₂⁻¹, multiplying det by det h₁/det h₂, and the pairs with det h₁ = det h₂, modulo scalars (z, z), map onto SO(V) [M]. For lattices Λ₁, Λ₂ in Q_ℓ² with v(det Λ₁) − v(det Λ₂) = 2j even, L(Λ₁, Λ₂) = ℓ^{−j}Hom(Λ₂, Λ₁) is a unimodular lattice of V, and depends only on the homothety classes. SO(V) acts transitively on unimodular lattices, so every one is of this form. If h₁M₂(Z_ℓ)h₂⁻¹ = M₂(Z_ℓ) with det h₁ = det h₂, put u = h₁h₂⁻¹. Since 1 lies in h₁M₂(Z_ℓ)h₁⁻¹, u lies in M₂(Z_ℓ); its determinant is 1, so u ∈ GL₂(Z_ℓ) and h₁M₂(Z_ℓ)h₁⁻¹ = M₂(Z_ℓ); hence both h_i are scalars times elements of GL₂(Z_ℓ). So the unimodular lattices of V correspond bijectively to the pairs of vertices of the tree of Q_ℓ² for which v(det Λ₁) ≡ v(det Λ₂) mod 2.

The neighbours of L(o, o) = M₂(Z_ℓ) are its (ℓ + 1)² images under GL₂(Z_ℓ) × GL₂(Z_ℓ) of one of them, since that group acts transitively on the isotropic lines of M₂(F_ℓ), the matrices of rank 1. The lattice ℓ·Hom(Z_ℓ ⊕ ℓZ_ℓ, Z_ℓ ⊕ ℓ⁻¹Z_ℓ) meets M₂(Z_ℓ) in index ℓ, so it is a neighbour, and it is L(y₁, y₂) with y₁, y₂ at distance 1 from o. So L(x₁, x₂) and L(y₁, y₂) are neighbours exactly when d(x₁, y₁) = d(x₂, y₂) = 1.

**Step 2: fixed points.** g corresponds to a pair (h₁, h₂), and L(x₁, x₂) is fixed by g exactly when h₁ fixes x₁ and h₂ fixes x₂. The eigenvalues of g are the λ_a/μ_b, for λ_± and μ_± the eigenvalues of h₁ and h₂; they are distinct, so λ₊ ≠ λ₋ and μ₊ ≠ μ₋, and K_i = Q_ℓ[h_i] is a quadratic étale algebra. If 𝓛 is not empty, each h_i fixes a vertex, h_iΛ_i = ℓ^{k_i}Λ_i, and det h₁ = det h₂ gives k₁ = k₂; dividing both by ℓ^{k₁} we may take h_i a unit of the maximal order of K_i. Then h_i fixes the class of Λ exactly when h_iΛ = Λ, that is when Λ, viewed in K_i, is stable under Z_ℓ[h_i] = Z_ℓ + ℓ^{r_i}O_{K_i}. So the fixed vertices of h_i are the vertices of Vol(K_i, r_i), embedded in the tree. The condition of Step 1 compares v(det Λ₁) and v(det Λ₂); in terms of a parity chosen on each volcano it says that the two parities differ by a constant ε, which depends on where the two volcanoes sit in the tree. Given the values of K_i and r_i, (2) follows from Step 1.

**Step 3: K_i and r_i.** A contains K₁ ⊗ K₂, acting by left and right multiplication, and both have dimension 4, so A = K₁ ⊗ K₂; the adjoint of left multiplication by h for the form det is left multiplication by the conjugate of h, so the involution of A is the product of the conjugations. Put u_i = h_i/h̄_i. With ρ the automorphism of A that conjugates the second factor only, g·ρ(g) = h₁²/N(h₂) = u₁ and ρ(g)/g = u₂. On eigenvalues ρ exchanges x and y, so u₁ and u₂ have the values (xy)^{±1} and (y/x)^{±1}, and t_i = u_i + u_i⁻¹. Write h_i = α + βω with O_{K_i} = Z_ℓ[ω]; then r_i = v_ℓ(β) and

> t_i − 2 = (u_i − 1)²/u_i = (h_i − h̄_i)²/N(h_i) = β²(ω − ω̄)²/N(h_i).

So v_ℓ(t_i − 2) is 2r_i if K_i is unramified or split and 2r_i + 1 if K_i is ramified, and K_i = Q_ℓ(u_i) is obtained by adjoining a square root of t_i² − 4 = (t_i − 2)(t_i + 2), where t_i + 2 ≡ 4 is a unit. The identities of §4.2 are direct: t₁ + t₂ = ab and t₁t₂ = a² + b² − 4, with a = x + x⁻¹ and b = y + y⁻¹, so S = a + b, Π = ab; D = (t₁ − t₂)² = (a² − 4)(b² − 4), and det(g ∓ 1) = (2 ∓ a)(2 ∓ b). This proves (2) and the first half of (4).

When K₁ and K₂ are both unramified fields, give a vertex of level j the parity of j; then for ε = 0 the pair of tops is in 𝓛 and for ε = 1 it is not. In the other cases one volcano has an automorphism that changes parity (a translation, or the exchange of the two tops), so the two graphs are isomorphic. Both values of ε occur: replacing h₂ by its conjugate under diag(ℓ, 1) gives an element of SO(V) with the same eigenvalues, hence the same A-module V with another form, and moves the second volcano to the other parity.

**Step 4: D not a square.** Then k = Q_ℓ(√D) is unramified, and GL₂(k) acts on H₂(k) by x ↦ hx h̄ᵗ, multiplying det by N_{k/Q_ℓ}(det h); the pairs (h, c) with c² = N(det h), modulo the (z, N(z)), map onto SO(V) [M]. For an O_k-lattice Λ = hO_k² in k², L(Λ) = ℓ^{−v(det h)}·hH₂(O_k)h̄ᵗ is unimodular and depends only on the homothety class; as in Step 1 this is a bijection between the vertices of the tree of k² and the unimodular lattices of V, and neighbours correspond to neighbours: there are ℓ² + 1 of each, GL₂(O_k) is transitive on both, and the lattice for h = diag(1, ℓ⁻¹) meets H₂(O_k) in index ℓ. The fixed points of g are the O_k-lattices in K = k[h] stable under O_k[h], that is Vol(K, r) over k. The eigenvalues of g are λ_aσ(λ_b)/c, with σ extending the conjugation of k, so xy = λ₊/λ₋ = u and t₁ = u + u⁻¹ ∈ k, t₂ = σ(t₁). As in Step 3, v_k(t₁ − 2) is 2r or 2r + 1, and K is split over k exactly when (t₁ − 2)(t₁ + 2)/ℓ^w is a square in k. A unit of k is a square exactly when its norm to Q_ℓ is one, and N((t₁ − 2)(t₁ + 2)) = (S² − 4Π)S². Finally t₁ − 2 = (Π − 4)/2 + √D/2. This proves (3) and the rest of (4). ∎

**Remarks.**

- In (3) the class of the form plays no role beyond existence: scaling the form by ℓ gives a space with no unimodular lattice.
- The proof uses that SO of a rank-4 space is, up to centre, SL₂ × SL₂ or SL₂ over a quadratic extension. In general 𝓛 is the set of fixed points of F/q among the unimodular lattices of V, which are vertices of the building of SO(V). For rank 4 that building is a product of two trees or one tree; for rank 6 it is the building of SL₄ or of a unitary group in four variables, of dimension 3, and the fixed points of a torus element are not described by trees. I have not worked that case out.

---

## 6. Rank 6

Proposition O applies. Beyond it I have only computations. For each configuration of primes found with c = 4 and m = 6 (q = 5, ℓ = 3, eigenvalue +q), the script takes one or two polynomials, chooses a form with a unimodular lattice of the simplest kind, and follows the rational steps from it: the whole connected component ("all") when no prime is a pair 𝔩 ≠ 𝔩̄, and the ball of radius 3 otherwise. The table gives the number of polynomials of each configuration and, for one of them, the lattices by number of steps.

| primes with F ≡ q | polynomials | lattices followed | lattices by number of steps |
|---|---|---|---|
| (e=2, f=1, v=3), (e=4, f=1, v=1) | 194 | 2, all | 2 with 1 |
| (e=2, f=1, v=1), (e=4, f=1, v=3) | 13 | 26, all | 2 with 13, 6 with 4, 18 with 1 |
| (e=1, f=2, v=1), (e=2, f=2, v=1) | 8 | 26, all | 2 with 13, 24 with 1 |
| (e=2, f=2, v=1), pair (e=1, f=1, v=1) | 6 | 80, ball | 12 with 13, 20 with 4, 48 with 1 |
| pair (e=3, f=1, v=2) | 5 | 152, ball | 7 with 22, 10 with 13, 81 with 4, 54 with 1 |
| (e=3, f=2, v=2) | 4 | 44, all | 4 with 13, 4 with 4, 36 with 1 |
| pair (e=1, f=1, v=1), pair (e=2, f=1, v=1) | 2 | 287, ball | 25 with 22, 26 with 13, 176 with 4, 60 with 1 |
| (e=1, f=2, v=1), (e=2, f=1, v=1), (e=2, f=1, v=1) | 2 | 8, all | 8 with 4 |
| (e=2, f=1, v=1), (e=2, f=1, v=1), pair (e=1, f=1, v=1) | 1 | 80, ball | 12 with 13, 20 with 4, 48 with 1 |

The script records e, f, v and whether a prime is its own conjugate; two rows merge a prime that is inert over E₀ with one that is ramified over E₀. For ℓ = 3: 13 = ℓ² + ℓ + 1, 22 = 2ℓ² + ℓ + 1 and 4 = ℓ + 1.

Observations [E]:

- The most frequent configuration, two primes ramified over E₀ with (e, v) = (2, 3) and (4, 1), gives two lattices exchanged by one step, as in note 9 with n = 1. It accounts for 194 of the 235 polynomials.
- In all the others there are lattices with L/τL ≅ (Z/ℓ)⁴ and ℓ² + ℓ + 1 = 13 steps (three blocks of size 1), and lattices with (Z/ℓ²)² and 1 or ℓ + 1 steps. Some have lattices with 4 or 22 = 2ℓ² + ℓ + 1 steps (two blocks of size 1).
- The graphs contain cycles. With (e, v) = (2, 1) and (4, 3): 26 lattices; two with 13 steps, joined to each other; six with 4 steps, each joined to one of the two and to three of the six; eighteen with one step.
- Only one class of forms per polynomial is followed, and only one connected component.

---

## 7. Check [E]

`sixlines_deep.py` searches the polynomials t³h(t + q²/t), h a monic integer cubic with three real roots in (−2q, 2q), with ℓ^c exactly dividing P(q), c even; keeps the irreducible ones with no prime above ℓ ramified over E₀ at the eigenvalue −q; and computes m from the factorisation of h over Q_ℓ.

**Rank 4.** For each polynomial with c ≥ 4, and for one in 25 of those with c = 2:

- it computes S, Π, D, t₁, t₂ to 60 ℓ-adic digits, checks v_ℓ(D) = c, and derives the prediction of Theorem P;
- for several forms Tr(αx ȳ) on O_E in each class of valuations of α, it lists all unimodular F-stable lattices in a box of O_E-lattices that must contain them (all of 𝓛 when no K is split);
- it follows the steps from one of them, testing Proposition O on every lattice met;
- when no K is split it asserts that the lattices reached are exactly those listed, and that the graph and the predicted graph have the same colour-refinement invariant (eight rounds from the degrees); when some K is split it compares in the same way the ball of radius 2 around the starting lattice, with degrees, with the balls of the predicted graph;
- when v_ℓ(D) is odd it asserts that no form has a lattice.

Colour refinement does not prove that two graphs are isomorphic; it separates all the predicted graphs that occur from each other.

| c | predicted graph | polynomials | lattices met |
|---|---|---|---|
| 2 | over Q_{ℓ²}: ramified r = 0 | 64 | 530 |
| 2 | ramified r = 0 × ramified r = 0 | 14 | 80 |
| 2 | ramified r = 0 × ramified r = 1 | 8 | 192 |
| 2 | ramified r = 0 × split r = 1 | 11 | 678 |
| 2 | ramified r = 0 × split r = 2 | 2 | 186 |
| 2 | ramified r = 0 × unramified r = 1 | 16 | 240 |
| 2 | ramified r = 0 × unramified r = 2 | 2 | 102 |
| 4 | over Q_{ℓ²}: ramified r = 0 | 91 | 1,062 |
| 4 | over Q_{ℓ²}: split r = 1 | 1 | 174 |
| 4 | over Q_{ℓ²}: unramified r = 1 | 31 | 1,023 |
| 4 | ramified r = 0 × ramified r = 0 | 100 | 1,224 |
| 4 | ramified r = 1 × split r = 1 | 10 | 2,640 |
| 4 | ramified r = 1 × unramified r = 1 | 9 | 540 |
| 4 | split r = 1 × unramified r = 1 | 22 | 3,186 |
| 4 | split r = 2 × unramified r = 1 | 1 | 255 |
| 6 | over Q_{ℓ²}: ramified r = 0 | 7 | 84 |
| 6 | ramified r = 0 × ramified r = 0 | 9 | 116 |
| 6 | ramified r = 1 × ramified r = 1 | 1 | 32 |
| 6 | ramified r = 1 × split r = 2 | 1 | 408 |
| 6 | split r = 1 × split r = 1 | 1 | 681 |
| 6 | unramified r = 1 × unramified r = 1 | 1 | 75 |

All agree: 402 polynomials (117 with c = 2, 265 with c = 4, 20 with c = 6) and 13,508 lattices, counted once per form tried. Every case with c = 2 has a factor that is a single edge, as the first consequence of Theorem P says. In the one case with two unramified fields, forms with ε = 0 and forms with ε = 1 both occur, and each graph matches only its own prediction.

Two limits of the run. When the box holds more than 1,500 lattices at some stage it is not listed; the steps are then followed from a lattice of the simplest kind (O_E plus a maximal isotropic subspace of its dual), and the assertion that the steps reach every lattice is not made. For three polynomials some forms had neither a listed box nor such a lattice and were skipped; other forms of the same polynomials were tested.

**Rank 6.** The polynomials of §6: Proposition O is tested on every lattice met: 17 polynomials, 1,295 lattices. For one of them no starting lattice of the simplest kind was found and nothing was followed.

**Rank 2** is not run: Proposition Q leaves nothing to compare beyond note 7.

---

## 8. What this does not give

- Rank 6 with c ≥ 4, beyond the local rule: here there is no description of the graph, no proof of the observations of §6, and only one class of forms followed. For the six-line family this is the remaining case. For c = 4 it is treated in [note 11](11-rank-six-fourth-power.md) to the extent stated there: the observations of §6 are explained (Proposition U, Theorem R, Corollary 1), and the graph is that of an explicit quartic order (Theorem S).
- Ranks above 6, which occur for other families with ρ̄ ≥ 12 (rank up to 10).
- Which ε, and in general which class of forms, the transcendental lattice of a surface has.
- Surfaces. All polynomials here come from a search and are not known to be Frobenius polynomials of surfaces of the family; the data of notes 2 to 5 has no case with c ≥ 4.
- Realisation of the steps by twisted equivalences: as in notes 8 and 9, through Propositions G and H and their unchecked steps. I have not re-read those proofs for the case where L/τL is not (Z/ℓ)².
- The comparison of graphs is by an invariant, and for infinite graphs by one ball of radius 2.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_deep.py --out sixlines_deep.jsonl > sixlines_deep.log
```

It takes about 13 minutes. An assertion fails if a graph differs from the prediction or if the local rule fails on a lattice. Options: `--search q,l`, `--both-signs`, `--cmin`, `--u-every`, `--rank6`, `--radius`.

## Sources

As in notes [6](06-torsors-under-cm-class-groups.md), [7](07-ring-changing-steps.md). In addition, recalled and not re-read [M]:

- the classification of unimodular quadratic lattices over Z_ℓ, ℓ odd, by rank and determinant (O'Meara, *Introduction to Quadratic Forms*, §92);
- the isomorphisms between SO of a rank-4 space and SL₂ × SL₂ or SL₂ over a quadratic extension, modulo centres (for instance Knus, Merkurjev, Rost, Tignol, *The Book of Involutions*, §15);
- the local structure of the lattices over a quadratic order (Kohel's thesis, as in note 7), and its analogue over an unramified extension.
