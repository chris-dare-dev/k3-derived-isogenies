# 12. Multiplicity 6 with ℓ⁴: the Galois action in the non-square case

*Makes explicit, for multiplicity 6 with ℓ⁴ dividing P(±q) and −P(q)P(−q) not a square in Q_ℓ, the Galois action that [note 11](11-rank-six-fourth-power.md) describes only as an action over the unramified quadratic extension; builds on notes 10 and 11. Labels: **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_hermitian.py`. Not re-checked.*

---

## 1. Result

Notation is that of note 11: ℓ odd, prime to q and to disc T; the eigenvalue is q (for −q replace F by −F); V is the part of T ⊗ Q_ℓ on which F ≡ q, of rank 6; g = F/q; ℓ⁴ is the exact power of ℓ dividing P(q); 𝓛 is the set of unimodular g-stable lattices of V. Here −P(q)P(−q) is not a square in Q_ℓ. As P(−q) ≡ (2q)⁶ mod ℓ is a square, this says that

> u = −P(q)/ℓ⁴

is a unit that is not a square. Let k = Q_ℓ(w), w² = u, be the unramified quadratic extension, with ring of integers O_k and conjugation σ.

1. **The Galois action is hermitian duality (Theorem T).** Over k the lattices are the lattices of a quartic order (note 11, Theorem S). That order comes from an order S₀ = Z_ℓ[X]/f₀ over Z_ℓ, with f₀ explicit in terms of P, and the conjugation of k sends the class of a lattice Λ to the class of its dual Λ^# for a hermitian form Tr(β v v̄′) on B₀ ⊗ k, B₀ = Q_ℓ[X]/f₀.
2. **So 𝓛 is a set of self-dual hermitian lattices.** The lattices of 𝓛 are the O_k-lattices M of B₀ ⊗ k stable under S₀ with M^# = M (type 0) or M^# = ℓ⁻¹M (type 1). A rational step joins a lattice M of type 0 and a lattice M′ of type 1 with ℓM ⊂ M′ ⊂ M. The kinds of note 11 are read off the Jordan blocks of X on M/ℓM.
3. **The form enters through β only, and decides whether there are leaves (Corollary 3).** β is determined by the form up to norms and scalars. There are leaves for exactly one class of β, that of 1/f₀′(X); for it S₀ ⊗ O_k itself is a leaf. For the other classes every lattice is middle or wide.
4. **Number of leaves (Corollary 4).** If no factor of f₀ over Q_ℓ has even residue degree the graph is finite, and for the class with leaves the number of leaves of one type is ℓ^c·Π_𝔭(1 + 1/N𝔭)/(1 + 1/ℓ), where 𝔭 runs over the primes of B₀ and ℓ^c is the index of S₀ in the maximal order.
5. **Check [E].** For ℓ = 3 and q = 5, 7, on the 1,130 searched polynomials with −P(q)P(−q) not a square (105 and 1,025): the identities of Theorem T (1) hold to 40 ℓ-adic digits; the graphs of all 1,200 classes of β were computed (14,334 lattices) and satisfy items 2 to 4; and for the 927 polynomials where the computation of note 11 has a starting lattice (1,779 forms, 55,464 lattices) the graph of 𝓛 has the same invariant as the graph of a class of β, the one with leaves being that of the leaf class. Details and limits are in §6.

What this is not: a closed formula for the whole graph in each configuration. As in note 11, the middle and wide lattices are obtained by a finite computation, now with lattices of rank 8 over Z_ℓ. And the rule that attaches β to a given form on O_E is not made explicit; §7 lists what is missing.

---

## 2. Statement

Write P(t) = t⁶ − e₁t⁵ + …, so that A₂ = e₁/q is the trace of g. Let ρ be the square root of P(−q)/q⁶ in Z_ℓ with ρ ≡ 8 mod ℓ, and

> f₀(X) = X⁴ + [ (2uℓ²/q³)·X³ + (12 − 2A₂)u·X² − (2u²ℓ²/q³)·X + (2 + A₂ − ρ)u² ] / (2 + A₂ + ρ).

Put B₀ = Q_ℓ[X]/f₀, S₀ = Z_ℓ[X]/f₀, θ the class of X, B = B₀ ⊗ k, S = S₀ ⊗ O_k, and v ↦ v̄ the map 1 ⊗ σ of B. For β in B₀^× let

> h_β(v, v′) = Tr_{B/k}(β·v·v̄′),

a hermitian form on the k-vector space B, and M^# = {v : h_β(v, M) ⊂ O_k} the dual of an O_k-lattice M.

**Theorem T.** [D] Suppose −P(q)P(−q) is not a square in Q_ℓ.

1. f₀ has coefficients in Z_ℓ and f₀ ≡ X⁴ mod ℓ. The element ξ = (w + θ)/(w − θ) of S satisfies ξ·ξ̄ = 1 and N_{B/k}(ξ) = 1, and the eigenvalues of g are the six products μ_iμ_j, i < j, of the four roots μ_i of the characteristic polynomial of ξ.
2. For every form for which 𝓛 is not empty there is a β in B₀^× such that 𝓛 is in bijection with the set of O_k-lattices M of B that are stable under S and satisfy M^# = M (*type 0*) or M^# = ℓ⁻¹M (*type 1*).
3. Two lattices of 𝓛 are joined by a rational step exactly when they are of different types and the corresponding M of type 0 and M′ of type 1 satisfy ℓM ⊂ M′ ⊂ M. Then M/M′ and M′/ℓM have dimension 2 over the residue field F_{ℓ²} of k.
4. The Jordan blocks of θ on M/ℓM over F_{ℓ²} and of τ on L/ℓL correspond: (4) to (1, 5), a leaf; (1, 3) to (3, 3), a middle lattice; (2, 2) to (1, 1, 1, 3) and (1, 1, 2) to (1, 1, 2, 2), the wide lattices. The leaves are the lattices xS, x in B^×; the lattices that are not leaves are those stable under θ³/ℓ, and the wide ones those stable under θ²/ℓ.
5. β is determined by the form up to multiplication by the elements b·b̄, b in B^×, and by Q_ℓ^×; forms that differ by a unit of Z_ℓ have the same β. Multiplying β by ℓ exchanges the two types. Necessarily v_ℓ(N_{B₀/Q_ℓ}(β·f₀′(θ))) is even.

**The Galois action.** [D] With λ₀ a square root of q in O_k, x = λ₀ξ generates the same order, S = O_k[x], and its characteristic polynomial is the polynomial f of note 11, Theorem S, for V ⊗ k, where κ = 1, for one choice of the signs of s₊ and s₋. So the classes of lattices of that theorem over k are the classes up to k^× of the S-stable O_k-lattices Λ of B, and the conjugation of k acts on them by

> [Λ] ↦ [Λ^#],  Λ^# = (β·f₀′(θ))⁻¹ · conj(S : Λ),  (S : Λ) = {b in B : bΛ ⊂ S},

conj denoting the image under v ↦ v̄.

On the generator, x̄ = ±q/x, with the sign + exactly when q is a square mod ℓ. The fixed classes are those of the lattices of types 0 and 1.

**The classes of β.** [D] Let B₀ = ΠC_i be the decomposition into fields, with ramification indices e_i and residue degrees f_i, and v_i the normalised valuations. C_i ⊗ k is a field exactly when f_i is odd. The class of β in the sense of (5) is given by the parities of v_i(β·f₀′(θ)) at the i with f_i odd, up to adding the vector (e_i mod 2); their sum is even. So there are 2^{n−1} classes if n ≥ 1 factors have odd f_i and all these e_i are even, 2^{n−2} if some e_i is odd, and one class if n = 0. I call *leaf class* the class of β = 1/f₀′(θ), where all parities are 0.

---

## 3. Consequences

Let T be the group of the x in B^× with x·x̄ = 1. Multiplication by an element of T maps S-stable lattices to S-stable lattices and commutes with the dual, so T acts on the graph of 𝓛, preserving types.

**Corollary 3 (leaves).** [D] 𝓛 contains a leaf exactly when β is in the leaf class. For β = 1/f₀′(θ): S is a leaf of type 0; the leaves of type 0 are the xS with x in T, a single orbit of T; and there are leaves of type 1 exactly when e_i is even for every i with f_i odd.

For β outside the leaf class every lattice of 𝓛 is middle or wide. In the runs of §6 these graphs are, in all 70 cases found, eight middle lattices forming the complete bipartite graph on 4 + 4 vertices.

**Corollary 4 (finiteness, number of leaves).** [D] If some f_i is even, 𝓛 is empty or infinite. If all f_i are odd, 𝓛 is finite, and in the leaf class the number of leaves of type 0 is

> ℓ^c · Π_i (1 + ℓ^{−f_i}) / (1 + ℓ⁻¹),  ℓ^c = [O_{B₀} : S₀],

and the number of leaves of type 1 is the same or zero, according to Corollary 3.

For example, in the pattern (3/2, ¼, ¼) B₀ is a totally ramified quartic field and S₀ its maximal order: one leaf of each type, which is Corollary 1 of note 11. In the pattern (½, ¾, ¾) B₀ is again a totally ramified quartic field, with c = 2 for ℓ = 3 in the examples: 9 + 9 leaves.

**Remark (levels; not tested).** [D] The argument of Corollary 3 applies to any order R with S ⊂ R ⊂ O_B in place of S. The ring of the b in B with bM ⊂ M, for M in 𝓛, is such an order, stable under conjugation, so R = R₀ ⊗ O_k with R₀ an order of B₀ containing S₀. If R₀ is a Gorenstein order, the lattices of 𝓛 with ring exactly R are the xR with x·x̄ in a fixed class, and those of one type are one orbit of T, with [T : T ∩ R^×] elements, or there are none. So 𝓛 is cut into levels by the orders R₀, as a volcano is, and the leaves are the level R₀ = S₀. In the pattern (½, ¾, ¾) with ℓ = 3 the numbers 9 + 3 + 1 of leaves, middle and wide lattices of one type are the indices ℓ², ℓ, 1 of a chain S₀ ⊂ S₀′ ⊂ O_{B₀}, which fits. This does not give the steps between levels, and it says nothing when R₀ is not Gorenstein, which happens for orders of degree 4; I have not computed the rings in the examples.

**What the graphs look like.** From the runs of §6, ℓ = 3, q = 5 and 7 together. For each configuration of primes of E: the factors of f₀ over Q_ℓ as (degree, e, f), and for each class of β the lattices found, as wide / middle / leaves. Where two graphs are given for the leaf class they belong to different polynomials.

| pattern | primes of E (e, f, v) | factors of f₀ (degree, e, f) | polynomials | leaf class | other class |
|---|---|---|---|---|---|
| (3/2, ¼, ¼) | (2, 1, 3), (4, 1, 1) | (4, 4, 1) | 816 | 2, all: 0 / 0 / 2; leaves 1 + 1, c = 0 | none |
| (½, ¾, ¾) | (2, 1, 1), (4, 1, 3) | (4, 4, 1) | 74 | 26, all: 2 / 6 / 18; leaves 9 + 9, c = 2 | none |
| (1, ½, ½) | (1, 2, 1), (2, 2, 1) | (1, 1, 1), (1, 1, 1), (2, 2, 1) | 31 | 53, all: 5 / 0 / 48; leaves 48 + 0, c = 3 (28 polynomials); or 161, all: 17 / 0 / 144; leaves 144 + 0, c = 4 (3) | 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (1, 2, 1), (2, 2, 1) | (2, 2, 1), (2, 2, 1) | 37 | 26, all: 2 / 0 / 24; leaves 12 + 12, c = 2 (32); or 80, all: 8 / 0 / 72; leaves 36 + 36, c = 3 (5) | 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (1, 2, 1), (2, 1, 1), (2, 1, 1) | (2, 2, 1), (2, 2, 1) | 2 | 80, all: 8 / 0 / 72; leaves 36 + 36, c = 3 | 8, all: 0 / 8 / 0 |
| (1, ½, ½) | (1, 2, 1), pair (2, 1, 1) | (2, 1, 2), (2, 2, 1) | 23 | infinite; ball of 17 around S: 5 / 0 / 12 | none |
| (1, ½, ½) | (2, 2, 1), pair (1, 1, 1) | (4, 2, 2) | 74 | infinite; ball of 44 around S: 8 / 12 / 24 | none |
| (1, ½, ½) | (2, 1, 1), (2, 1, 1), pair (1, 1, 1) | (4, 2, 2) | 25 | infinite; ball of 44 around S: 8 / 12 / 24 | none |
| (⅔, ⅔, ⅔) | (3, 2, 2) | (1, 1, 1), (3, 3, 1) | 48 | 44, all: 5 / 3 / 36; leaves 36 + 0, c = 3 | none |

"Leaves a + b" are the numbers of leaves of types 0 and 1; they are the numbers of Corollary 4, for instance 3³·(4/3)³/(4/3) = 48 with three primes of residue degree 1 and 3³·(4/3)²/(4/3) = 36 with two. In the two rows with the primes (1, 2, 1), (2, 2, 1), the prime (2, 2, 1) is inert over E₀ in the first and ramified over E₀ in the second (looked up for the three polynomials of the run with q = 5). In the graphs of 8 lattices each lattice of one type is joined to the four of the other type. The graphs with 17 and with 8 wide lattices appear only for q = 7; I have not checked that they are the trees of depth 2 that Theorem R of note 11 would give. The types also explain a feature of the graph of 44 lattices of note 11: its 4 lattices of type 1 are the outer wide lattices, and the centre, the middle lattices and the leaves are of type 0.

---

## 4. Proof of Theorem T

**Step 0: a hermitian model of V.** Let W = k⁴ with the hermitian form h(v, v′) = Σ v_i σ(v′_i), the lattice Λ₀ = O_k⁴ and vol = e₁ ∧ e₂ ∧ e₃ ∧ e₄. Let H = Λ²h, a hermitian form on Λ²W, and define ★ on Λ²W by ω ∧ ★ω′ = H(ω, ω′)·vol. Then ★ is σ-semilinear, ★(e₁ ∧ e₂) = e₃ ∧ e₄, ★(e₁ ∧ e₃) = −e₂ ∧ e₄, ★(e₁ ∧ e₄) = e₂ ∧ e₃, and ★² = 1. So V₀ = {ω : ★ω = ω} is a Q_ℓ-form of Λ²W, and on V₀ the form ω ∧ ω′/vol of note 11 equals H and takes values in Q_ℓ. A basis of V₀ is e₁₂ + e₃₄, e₁₃ − e₂₄, e₁₄ + e₂₃, w(e₁₂ − e₃₄), w(e₁₃ + e₂₄), w(e₁₄ − e₂₃), with Gram matrix diag(2, 2, 2, −2u, −2u, −2u): a unimodular lattice of determinant −u up to squares. By note 11, Step 0, det V = P(q)P(−q) = −u up to squares, and V has a unimodular lattice when 𝓛 is not empty. Unimodular lattices over Z_ℓ of the same rank and determinant are isometric [M], so V ≅ V₀, and V ⊗ k = Λ²W with the conjugation acting as ★.

**Step 1: the conjugate of a lattice.** For y in GL₄(k) with adjoint y† for h, H(Λ²y·ω, ω′) = H(ω, Λ²y†·ω′) and Λ²y·ω ∧ Λ²y·η = det(y)·ω ∧ η give

> ★ ∘ Λ²y ∘ ★⁻¹ = σ(det y) · Λ²(y†)⁻¹.

For Λ = yΛ₀ the dual lattice for h is Λ^# = (y†)⁻¹Λ₀. As ★ preserves Λ²Λ₀, ★(Λ²Λ) is a multiple of Λ²(Λ^#). In the notation of note 11, Step 1, where L(Λ) is the unimodular lattice ℓ^{−j}Λ²Λ of V ⊗ k: the conjugate of L(Λ) is L(Λ^#).

**Step 2: lattices of V.** A unimodular lattice L of V gives L ⊗ O_k, and this identifies the unimodular lattices of V with those of V ⊗ k fixed by the conjugation (note 11, Step 5). By note 11, Step 1, over k, the latter are the L(Λ), [Λ] a class with v(det) even relative to Λ₀, and by Step 1 above the fixed ones are those with Λ^# = cΛ, c in k^×. Replacing Λ by zΛ replaces c by c/(zσ(z)), so exactly one multiple of Λ, up to units, has Λ^# = Λ or Λ^# = ℓ⁻¹Λ. Conversely, if Λ^# = ℓ^{−a}Λ, comparing the indices of Λ and of Λ^# relative to Λ₀ = Λ₀^# gives v(det Λ) = 2a, which is even. So the unimodular lattices of V correspond to the lattices of types 0 and 1 of W.

Two lattices of V are ℓ-neighbours exactly when, after scaling, ℓΛ ⊂ Λ′ ⊂ Λ with quotients of dimension 2 over F_{ℓ²} (note 11, Steps 1 and 5). If Λ is of type 0 then Λ′^# contains Λ′ with a quotient of length 4, so Λ′ is of type 1 as it stands. If Λ is of type 1 the same count gives Λ′^# = ℓ⁻²Λ′, so ℓ⁻¹Λ′ is of type 0 and contains Λ with ℓ·(ℓ⁻¹Λ′) ⊂ Λ. So neighbours are of different types, as stated in (3).

**Step 3: g.** Over k, g = c⁻¹Λ²y with c² = det y (note 11, Step 1). By the formula of Step 1, g commutes with ★ exactly when Λ²(y·y†) is the scalar c·σ(c), that is when y·y† is a scalar m, which is then in Q_ℓ^×. Let L(Λ) be fixed by g, Λ of type 0 or 1. Then yΛ = zΛ, and taking duals m = zσ(z) up to a unit. Units of Z_ℓ are norms from k, so after dividing y by an element of k^× it is unitary and preserves Λ. On Λ/ℓΛ, with the form h or ℓ⁻¹h, it is a unitary transformation whose eigenvalues λ_i satisfy λ_iλ_j ≡ c for all i < j, because the eigenvalues of g are ≡ 1. So they are all congruent to one λ, and λσ(λ) ≡ 1 because the set of eigenvalues of a unitary transformation is stable under λ ↦ σ(λ)⁻¹. Dividing by a unit of norm 1 that lifts λ, and then by a fourth root, of norm 1 and ≡ 1, of the determinant, we get

> ξ in SU(h), unipotent on Λ/ℓΛ, with g = Λ²ξ.

So g fixes L(Λ) exactly when ξΛ = Λ.

**Step 4: the algebra.** The eigenvalues of g are the μ_iμ_j, μ_i those of ξ, and they are distinct; so the μ_i are distinct, B = k[ξ] is an étale algebra of degree 4 and W is free of rank one over it. The adjoint for h maps ξ to ξ⁻¹, so it preserves B and induces on it a σ-semilinear involution b ↦ b̄. Its fixed points B₀ form a Q_ℓ-algebra with B₀ ⊗ k = B. As ξ + 1 is invertible on Λ, t = (ξ − 1)/(ξ + 1) is defined, t̄ = −t, and θ = wt lies in B₀; k[θ] = k[ξ], so B₀ = Q_ℓ[θ], and O_k[ξ] = O_k[θ] = Z_ℓ[θ] ⊗ O_k. Choosing a generator of W over B, h becomes Tr_{B/k}(β v v̄′) for a unique β with β̄ = β; another generator multiplies β by some b·b̄. A lattice is stable under ξ exactly when it is stable under O_k[ξ]. With Steps 2 and 3 this gives (2) and (3), once the characteristic polynomial of θ is identified.

Multiplying the form of V by a unit does not change its unimodular lattices nor g; multiplying h by a scalar of Q_ℓ^× does not change the group nor the duality up to scalars, and ℓ·h exchanges the types. The determinant of h_β on the basis θ^i of S is N_{B₀/Q_ℓ}(β)·disc(f₀) up to sign, and N(f₀′(θ)) = disc(f₀) up to sign; the form h has a unimodular lattice, so its determinant has even valuation. This is (5).

**Step 5: the polynomial.** Let X⁴ − A₁X³ + A₂X² − A₃X + 1 be the characteristic polynomial of ξ. Then A₂ = Σμ_iμ_j is the trace of g, and the identities of note 11, Step 3, read with q = κ = 1,

> det(1 + g) = (A₁ + A₃)²,  det(1 − g) = −(A₁ − A₃)².

The left sides are P(−q)/q⁶ and P(q)/q⁶ = −uℓ⁴/q⁶. From ξ̄ = ξ⁻¹, σ(A₁) = A₃. So ρ = A₁ + A₃ lies in Z_ℓ and is ≡ 8 since the μ_i are ≡ 1; and A₁ − A₃ is a square root of uℓ⁴/q⁶ that σ changes into its opposite, so it is wℓ²/q³ for one of the two square roots w of u, which is the one we take. Now

> (1 − t)⁴·χ_ξ((1 + t)/(1 − t)) = (2 + A₂ + ρ)t⁴ + 2(A₁ − A₃)t³ + (12 − 2A₂)t² − 2(A₁ − A₃)t + (2 + A₂ − ρ),

and substituting t = θ/w and multiplying by u² gives (2 + A₂ + ρ)·f₀(θ) = 0, with 2 + A₂ + ρ ≡ 16 a unit. The roots θ_i = w(μ_i − 1)/(μ_i + 1) have positive valuation, so f₀ is integral and ≡ X⁴. And ξ = (1 + t)/(1 − t) = (w + θ)/(w − θ). This is (1).

**Step 6: Jordan blocks and leaves.** ξ − 1 = 2t/(1 − t) is a unit multiple of θ in S, so they have the same blocks on M/ℓM. Over F_{ℓ²}, L/ℓL becomes Λ²(M/ℓM) with g acting as Λ²ξ, and the computation of note 11, Step 4, applies; the blocks of g do not change with the extension of scalars. S is a local ring with residue field F_{ℓ²}, so M/ℓM has one block exactly when M is generated by one element. The rest of (4) is as in note 11.

**Step 7: the formula for the dual.** For an S-stable lattice Λ, Λ^# = β⁻¹·(Λ̄)^∨, where ∨ is the dual for the trace form of B over k. For the order S generated by one element, S^∨ = f₀′(θ)⁻¹S [M], and (Λ̄)^∨ = {b : bΛ̄ ⊂ S^∨} = f₀′(θ)⁻¹(S : Λ̄). With x = λ₀ξ: the eigenvalues of F are the qμ_iμ_j = (λ₀μ_i)(λ₀μ_j), which is the normalisation of note 11, Step 3, with κ = 1; and x̄ = σ(λ₀)/ξ = σ(λ₀)λ₀/x, where σ(λ₀) = λ₀ exactly when q is a square in Q_ℓ. ∎

---

## 5. Proofs of the corollaries

**Norms.** For a field C_i ⊗ k, unramified of degree 2 over C_i, the norms are the elements of even valuation; when C_i ⊗ k is a product of two fields exchanged by the conjugation, every element of C_i^× is a norm [M]. This gives the description of the classes of β in §2; the vector (e_i mod 2) is the class of ℓ, and the condition of Theorem T (5) is that Σf_iv_i(βf₀′(θ)) is even, where the terms with f_i even do not count. The norm map from S^× to S₀^× is surjective, because S is étale of degree 2 over the local ring S₀, with residue fields F_{ℓ²} and F_ℓ [M].

**Corollary 3.** By Step 7, (xS)^# = (βf₀′(θ))⁻¹·x̄⁻¹S. So xS is of type 0 exactly when x·x̄ is in (βf₀′(θ))⁻¹·S₀^×, and of type 1 exactly when it is in ℓ·(βf₀′(θ))⁻¹·S₀^×. As S₀^× consists of norms, a leaf of type 0 exists exactly when βf₀′(θ) is a norm, and one of type 1 exactly when ℓβf₀′(θ) is. For β = 1/f₀′(θ): x·x̄ = s·s̄ with s in S^×, so xS = (x/s)S with x/s in T; and ℓ is a norm exactly when its valuation e_i is even in each C_i with f_i odd. ∎

**Corollary 4.** If some f_i is even, C_i ⊗ k is a product of two copies of a field C′ and T contains a copy of C′^×, which is not compact; the stabiliser in T of a lattice is compact, so every orbit is infinite. If all f_i are odd, B is a product of fields and T is the group O¹ of units of norm 1 of the maximal order O_B. A self-dual lattice M generates an ideal 𝔞 = O_B·M with 𝔣𝔞 ⊂ M ⊂ 𝔞, 𝔣 the conductor of S, and M^# = M bounds the valuations of 𝔞 above and below; so 𝓛 is finite. In the leaf class the leaves of type 0 are T/(T ∩ S^×) = O¹/S¹. From the two exact sequences 1 → S¹ → S^× → S₀^× → 1 and 1 → O¹ → O_B^× → O_{B₀}^× → 1,

> [O¹ : S¹] = [O_B^× : S^×] / [O_{B₀}^× : S₀^×] = ( ℓ^{2c}·Π(1 − ℓ^{−2f_i})/(1 − ℓ⁻²) ) / ( ℓ^c·Π(1 − ℓ^{−f_i})/(1 − ℓ⁻¹) ),

using for an order S′ that is a local ring with residue field of q′ elements, of index N in its maximal order, [O^× : S′^×] = N·Π_𝔭(1 − 1/N𝔭)/(1 − 1/q′) [M]. The leaves of type 1, if any, are a translate of those of type 0 by an element of norm ℓ times a unit. ∎

---

## 6. Check [E]

`sixlines_hermitian.py` runs over the same polynomials as `sixlines_rank6.py` and keeps those for which −P(q)P(−q) is not a square. An O_k-lattice of B is represented as a Z_ℓ-lattice of rank 8 in B₀ ⊕ wB₀ stable under multiplication by w, with the symmetric form Tr_{B₀/Q_ℓ}(β(v₁v₁′ − u·v₂v₂′)); this form is half the trace from k to Q_ℓ of h_β, and it has the same dual lattices. f₀ is computed to 40 ℓ-adic digits and replaced by an integer polynomial. For each polynomial the script

- asserts f₀ ≡ X⁴ mod ℓ, and, by exact arithmetic in Q(w) on the integer polynomial, that the characteristic polynomial of Λ²ξ agrees with that of g, and det ξ with 1, to at least 34 digits (Theorem T (1); ξ·ξ̄ = 1 holds by construction);
- factors f₀ over Q_ℓ, finds for each class of β a representative π(θ)/f₀′(θ) with π a polynomial with small coefficients, and asserts that the number of classes is the one of §2;
- for the leaf class asserts that S is of type 0; for the other classes looks for a stable self-dual lattice between S ∩ S^# and S + S^#;
- follows the steps, found by going through all planes of M/ℓM over F_{ℓ²} and keeping those that are stable under θ and totally isotropic: the whole component if all f_i are odd, the ball of radius 3 otherwise; on every lattice met it asserts the type, that the blocks are of one of the four kinds, and that the number of steps is 1, ℓ + 1, ℓ² + ℓ + 1 and ℓ + 1 for (4), (1, 3), (2, 2) and (1, 1, 2) (Theorem T (4), with note 11, Proposition U);
- asserts that leaves occur exactly in the leaf class (Corollary 3) and, in the finite case, that the numbers of leaves of the two types are those of Corollary 4;
- computes the graph of 𝓛 as in note 11, from forms on O_E (two for each class of valuations) and the steps along isotropic eigenlines, and asserts: in the finite case, that it has the same number of lattices and the same colour-refinement invariant as the graph of one of the classes of β, which is the leaf class exactly when the graph has leaves; in the infinite case, that the ball of radius 3 around a leaf has the same invariant as the ball around S (Theorem T (2), (3)).

**Run for q = 5, ℓ = 3, eigenvalue +q.** Of the 235 polynomials of note 11, 105 have −P(q)P(−q) not a square: 88 of the pattern (3/2, ¼, ¼), 5 of (½, ¾, ¾), 8 of (1, ½, ½), 4 of (⅔, ⅔, ⅔). All assertions pass. The characteristic polynomials agree to all 40 digits. There are 108 classes of β in all (three polynomials have two classes, the others one); a starting lattice was found for every one, and 858 hermitian lattices were met. On the side of O_E a starting lattice was found for 86 polynomials, 156 forms and 1,842 lattices; the comparison was made for all of them, on the whole graph for 150 forms and on a ball for 6 (two polynomials). For 84 of the 100 polynomials with a finite graph, among them the three with two classes, every class of β was met as the graph of some form on O_E; for the other 16 no form on O_E had a starting lattice of the kind the script of note 11 looks for.

**Run for q = 7, ℓ = 3, eigenvalue +q.** 1,025 polynomials with −P(q)P(−q) not a square (923 others are squares): 728 of the pattern (3/2, ¼, ¼), 69 of (½, ¾, ¾), 184 of (1, ½, ½), 44 of (⅔, ⅔, ⅔). All assertions pass. There are 1,092 classes of β (67 polynomials have two), a starting lattice was found for every one, and 13,476 hermitian lattices were met. On the side of O_E a starting lattice was found for 841 polynomials, 1,623 forms and 53,622 lattices; the comparison was made for all of them, on the whole graph for 1,313 forms and on a ball for 310 (88 polynomials). For 722 of the 908 polynomials with a finite graph every class of β was met as the graph of some form on O_E, among them 36 of the 67 with two classes.

**Eigenvalue −q.** The polynomials with eigenvalue −q are the P(−t) for the polynomials above, so they give nothing new. I ran `--both-signs` once for q = 5 as a test of the signs in the script: same classes and same graphs, all assertions pass. Its output is not kept.

**Limits.**

- The graphs are compared by an invariant, and in the infinite case on one ball; forms on O_E are not matched with classes of β, only the graphs are.
- For 203 polynomials (19 and 184) the script of note 11 found no starting lattice, and only the statements on the hermitian side were tested.
- Classes other than the leaf class occur for 70 polynomials (3 and 67), always with the same graph of 8 lattices; for 10 of them no form on O_E with that graph was met, and for 21 none with the graph of the leaf class, so there one of the two graphs is a prediction only.
- On the hermitian side only the component of the starting lattice is listed, without an independent list of all self-dual stable lattices; in the finite case its size equals that of the graph of note 11, which was checked there to be complete when the list was small enough.
- The formula for the action on the classes that are not fixed (§2) is not tested by itself; only its fixed points are.
- ℓ = 3 only: the steps are found by going through all planes over F_{ℓ²}, which is slow for larger ℓ.
- As in note 11, the polynomials are searched ones. For q = 5 the wide lattices are an edge or a star, depth n = 1 in Theorem R of note 11; for q = 7 larger sets of wide lattices occur, which I have not compared with that theorem.

---

## 7. What this does not give

- A closed formula for the graph in each configuration. Outside the leaf count of Corollary 4 and the statements of note 11 (Proposition U, Theorem R), the middle and wide lattices and the way they are joined come from a computation for each polynomial; the descriptions of §3 are observations.
- The rule that attaches the class of β to a form given on O_E, that is to an element α of E₀. In the check the graphs are matched, not the forms. Nor is it proved here that different classes of β come from forms that differ by more than a unit; the counts agree in the examples.
- A test for ℓ other than 3, and of depth n > 1 in Theorem R, as in note 11.
- Multiplicity 6 with ℓ⁶ or more. Steps 0 to 4 and 6, 7 of the proof do not use the exponent 4; the formula for f₀ (through u) and the list of block types do.
- As in note 11: ranks above 6; which class of forms the lattice of a surface has; surfaces at all, the polynomials being searched ones and the six-line data having no such case; the realisation of the steps by twisted equivalences.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_hermitian.py --out sixlines_hermitian.jsonl > sixlines_hermitian.log
```

It takes about a minute. An assertion fails if any statement of the list of §6 fails. Options: `--search q,l`, `--both-signs`, `--forms`, `--radius`. The run for q = 7, which takes about seven minutes, is

```bash
python -I sixlines_hermitian.py --search 7,3 --out sixlines_hermitian_q7.jsonl > sixlines_hermitian_q7.log
```

## Sources

As in notes [10](10-fourth-power.md) and [11](11-rank-six-fourth-power.md). In addition, recalled and not re-read [M]: the classification of unimodular quadratic lattices over Z_ℓ, ℓ odd, by rank and determinant (O'Meara, *Introduction to Quadratic Forms*, §92); the isomorphism between the special unitary group of a hermitian space of dimension 4 and the spin group of a quadratic space of rank 6 whose discriminant is the quadratic extension (Knus, Merkurjev, Rost, Tignol, *The Book of Involutions*, §15.D); the dual of Z[θ] for the trace form (Serre, *Local Fields*, ch. III, §6); the norm group of an unramified extension (Serre, ch. V, §2).
