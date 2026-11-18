# 8. The structure below the top level at eigenvalue ±q

*Describes the graph of lattices below the top level when the eigenvalue is ±q and ℓ² exactly divides P(±q): an ℓ-volcano with explicit depth and rings (Theorem M), checked on the six-line data. Builds on [note 7](07-ring-changing-steps.md), and proves, in two of the four configurations, the statement left open there that a stable eigenline leads to a larger ring. Labels: **[D]** derived here; **[E]** measured by `experiments/sixlines_levels.py`. Nothing here has been independently verified.*

---

## 1. Result

Notation is that of notes 6 and 7. ℓ is an odd prime not dividing q or disc T, and the eigenvalue is q; for −q replace F by −F, that is X by its quadratic twist. Only the part of T ⊗ Z_ℓ on which F ≡ q matters, and steps with eigenvalue q change nothing else.

**Hypothesis (U).** ℓ² is the exact power of ℓ dividing P(q). Equivalently Br(X)[ℓ^∞] has order ℓ². This holds in all 48 pairs (polynomial, ℓ) of the data with odd torsion in the Brauer group.

Under (U):

1. **Every lattice has 0, 1, 2 or ℓ + 1 rational steps with eigenvalue q** (Proposition L). It has ℓ + 1 exactly when ℓ/(F − q)² preserves it.
2. **There are four configurations of primes** (Lemma 3). Two of them have a lattice that is maximal at ℓ: (a) one prime 𝔩 = 𝔩̄, inert over E₀, ramified over ℓ; (b) a pair 𝔩 ≠ 𝔩̄ of ramified primes of degree one.
3. **In (a) and (b) the graph is an ℓ-volcano** (Theorem M). Write ℓ^{2n+1} for the index of Z_ℓ[F] in the maximal order, on the part where F ≡ q; it is always of this shape and n ≥ 1. Then:
   - every lattice L has a *level* i with 0 ≤ i ≤ n, given by [O_E·L : L] = ℓ^{2i};
   - End(L) = Z_ℓ[F] + 𝔏^{2i}O_E at ℓ, with 𝔏 = 𝔩 in (a) and 𝔩𝔩̄ in (b); its index in the maximal order is ℓ^{2i}, and it is never Z_ℓ[F];
   - from a lattice of level i ≥ 1 exactly one step goes up, to (F − q)L + 𝔏^{i−1}O′ where O′ is its top; if i < n there are ℓ more, all going down; if i = n there are none;
   - from a top there are ℓ + 1 steps down in (a), and two translations and ℓ − 1 steps down in (b);
   - the lattices below one top form a tree, with (ℓ + 1)ℓ^{i−1} of them at level i in (a) and (ℓ − 1)ℓ^{i−1} in (b);
   - the eigenline going up is the only one stable under End(L).
4. **Every edge is a twisted equivalence defined over F_q**, with Brauer classes of order ℓ, and Br(X)[ℓ^∞] ≅ (Z/ℓ)² at every vertex (Propositions G and H).
5. **Check [E].** For the 9 pairs of note 7 in configurations (a) and (b), all the statements above were tested on each of the 221 lattice classes and their 410 neighbours, including the formula for the parent. All hold.

So the picture of the elliptic-curve volcano holds here as well, with one difference from the volcano at a double root (note 7, Proposition K): along this one the Brauer class of every edge is rational.

---

## 2. Setting

Let V be the part of T ⊗ Q_ℓ on which F ≡ q, of dimension m (the multiplicity of q as a root of P mod ℓ), A = Q_ℓ[F] acting on V, A₀ the subalgebra fixed by the involution, and O, O₀ their maximal orders. Put

> τ = (F − q)/(F + q).

Since F̄ = q²/F, τ̄ = −τ, so ⟨τx, y⟩ = −⟨x, τy⟩. On V, F + q is a unit of the local ring Z_ℓ[F] and F = q(1 + τ)/(1 − τ), so Z_ℓ[F] = Z_ℓ[τ] =: R. A lattice is F-stable exactly when it is τ-stable, and the eigenspace of F for q on L/ℓL is the kernel of τ. Let 𝓛 be the set of unimodular τ-stable lattices in V.

Hypothesis (U) says that det τ on V is ℓ² times a unit, so [L : τL] = ℓ² for every τ-stable lattice.

**Lemma 1.** [D] Under (U), for every L in 𝓛: L/τL ≅ (Z/ℓ)², so ℓL ⊂ τL ⊂ L ⊂ τ⁻¹L ⊂ ℓ⁻¹L, and the eigenspace W(L) = τ⁻¹(ℓL)/ℓL is a plane.

*Proof.* coker(F − q) = L/τL has order ℓ² and is of the shape A × A (note 7, Proposition I (3)). ∎

**Lemma 2 (edges).** [D] Assume (U). If L and L′ in 𝓛 are joined by a step, then M = L + L′ is τ-stable and M^∨ = L ∩ L′ = τM. Conversely, if M is τ-stable with M^∨ = τM and ℓM ⊂ τM, then M/τM is a non-degenerate quadratic plane over F_ℓ; the members of 𝓛 between τM and M are its isotropic lines, so there are two or none; and if there are two they are joined by a step.

*Proof.* Write L′ = (L ∩ L′) + Z_ℓy with y = x̃/ℓ. Then τy ∈ L because τx ∈ ℓL, so τL′ ⊂ L, and by symmetry τL ⊂ L′. So τM ⊂ L ∩ L′, and both have index ℓ² in M. Also M^∨ = L^∨ ∩ L′^∨. Conversely the self-dual lattices between M^∨ and M are the Lagrangian subgroups of M/M^∨; each contains τM, hence is τ-stable; two of them meet in τM, of index ℓ in each; and for y in L′, ℓy lies in L with τ(ℓy) ∈ ℓL. ∎

Call such an M with two isotropic lines an *edge*.

**Proposition L.** [D] Assume (U) and let L be in 𝓛. The number of rational steps from L with eigenvalue q is

- ℓ + 1 if ℓL ⊂ τ²L, that is if ℓ/τ² lies in End(L);
- otherwise 1 if m ≥ 4, and 0 or 2 if m = 2.

*Proof.* By note 7, Proposition I (3), the radical of the form on W is W ∩ im(τ), taken in L/ℓL. Here W = τ⁻¹(ℓL)/ℓL and im(τ) = τL/ℓL, so W ⊂ im(τ) exactly when ℓL ⊂ τ²L, and then all ℓ + 1 lines are isotropic. Otherwise look at τ on L/ℓL restricted to its generalised kernel, of dimension m: it is nilpotent with a kernel of dimension 2, so it has two Jordan blocks, and the radical is spanned by the kernel vectors of the blocks of size ≥ 2. If W is not in the image, one block has size 1. For m ≥ 4 the other has size ≥ 3, the radical is a line, the form on W has rank 1, and its only isotropic line is the radical. For m = 2 the form is non-degenerate. ∎

**Lemma 3 (configurations).** [D] Assume (U) and m ≥ 4. The primes of E above ℓ at which F ≡ q are one of:

- (a) one prime 𝔩 = 𝔩̄, inert over E₀, with f = 2 and e = m/2;
- (b) two primes 𝔩 ≠ 𝔩̄, with f = 1 and e = m/2;
- (c) one prime 𝔩 = 𝔩̄, ramified over E₀, with f = 2;
- (d) two primes, each equal to its conjugate and ramified over E₀, with f = 1.

In every case F − q has valuation 1 at these primes. A lattice maximal at ℓ can exist only in (a) and (b).

*Proof.* 2 = Σ f_𝔩·v_𝔩(τ) with every v_𝔩(τ) ≥ 1. At a prime ramified over E₀ a skew element has odd valuation, and a prime with f = 1 equal to its conjugate is ramified over E₀. So the sum is one prime with f = 2 and valuation 1, necessarily its own conjugate, or two primes with f = 1 and valuation 1, conjugate to each other or each its own conjugate. The last sentence is note 7, Proposition J (1). ∎

In the data: (a) 7 pairs, (b) 2 pairs, (c) 1 pair, (d) 3 pairs; the other 35 pairs with (U) have m = 2.

---

## 3. Theorem M

Assume (U), configuration (a) or (b) with e ≥ 2, and that 𝓛 contains a lattice maximal at ℓ. Call the unimodular O-lattices *tops*. Put 𝔏 = 𝔩 in (a) and 𝔏 = 𝔩𝔩̄ in (b), and R₀ = Z_ℓ[τ²] ⊂ O₀.

**Theorem M.** [D]

1. *Depth.* [O₀ : R₀] = ℓ^n and [O : R] = ℓ^{2n+1} for an integer n ≥ 1. If e is odd, n = (e − 1)/2. If e = 2, n is the valuation of the trace of τ/ω (notation of §4).
2. *Level.* For L in 𝓛, [O·L : L] = ℓ^{2i} with 0 ≤ i ≤ n. There is a unique top O′ = top(L) with 𝔏^iO′ ⊂ L ⊂ 𝔏^{−i}O′, and L = R·x + 𝔏^iO′ for every x in L that generates 𝔏^{−i}O′ over O. In (a) there is one top; in (b) the tops are the translates 𝔩^j𝔩̄^{−j}O′, j ∈ Z.
3. *Ring.* End(L) = R + 𝔏^{2i}O, of index ℓ^{2i} in O. In particular End(L) determines the level, and End(L) ≠ R.
4. *Steps from a top.* In (a): ℓ + 1 steps, all to level 1. In (b): two translations to the neighbouring tops and ℓ − 1 steps to level 1.
5. *Steps from level i ≥ 1.* Exactly one step goes to level i − 1; it leads to τL + 𝔏^{i−1}top(L). If i < n there are exactly ℓ other steps, to lattices of level i + 1 with the same top. If i = n there is no other step.
6. *Stable lines.* For i ≥ 1 the eigenline of the step going up is the only line of W(L) stable under End(L).
7. *Shape.* The lattices with a given top form a rooted tree under the steps, with (ℓ + 1)ℓ^{i−1} vertices at level i in (a) and (ℓ − 1)ℓ^{i−1} in (b), 1 ≤ i ≤ n.

**Consequences.**

- At a lattice of level n the plane W has one isotropic line and ℓ anisotropic ones; at lower levels all ℓ + 1 are isotropic. So a surface at the bottom has ℓ² − 1 rational Brauer classes of order ℓ, of which ℓ − 1 give a step.
- By Proposition H (2) and (3), every edge is realised by a twisted equivalence defined over F_q.
- In configurations (a) and (b) the statement left open in note 7 after Proposition I holds: the stable eigenline leads to a ring containing End(L), and the others do not.
- The exponent of ℓ in [O_E : Z[F, V]] is 2n + 1 when the other primes above ℓ contribute nothing. In the data it is 3 in eight pairs (n = 1) and 5 in one (n = 2).

---

## 4. Proof of Theorem M

**Normal form.** In (a) let W be the unramified quadratic extension of Q_ℓ, ω ∈ W with ω² = δ a non-square unit of Z_ℓ; then A = A₀ ⊗ W = A₀ ⊕ A₀ω with A₀ a totally ramified field of degree e. In (b), A = K × K with K totally ramified of degree e, the involution exchanges the factors, A₀ = K embedded diagonally, and ω = (1, −1), δ = 1. In both cases ω̄ = −ω, O = O₀ ⊕ O₀ω (ℓ is odd), and

> π = τ/ω

lies in A₀ and is a uniformiser of O₀, by Lemma 3. So O₀ = Z_ℓ[π], τ² = δπ², R₀ = Z_ℓ[π²] and R = R₀ ⊕ R₀τ. Let k′ = O/πO = F_ℓ[ω], a plane over F_ℓ with the quadratic form N(a) = aā; it is anisotropic in (a), and in (b) it is hyperbolic and its two isotropic lines are the eigenlines of ω. For any line F_ℓa of k′, its orthogonal is F_ℓωa. A line is isotropic exactly when it is an eigenline of ω.

Fix a top O′; replacing V by a translate we write O′ = O. The form is ⟨x, z⟩ = λ₀(Re(x z̄)), with Re(x) = (x + x̄)/2 and λ₀: A₀ → Q_ℓ linear; since O is unimodular, {c ∈ A₀ : λ₀(cO₀) ⊂ Z_ℓ} = O₀. Valuations on A₀ are normalised by v(π) = 1.

**Step 1: R₀ and the number n.** R₀ contains elements of every even valuation. Let κ be the least odd valuation of an element of R₀. An element c₀ + c₁π² + … with c₀ not a unit has valuation ≥ 2 because e ≥ 2, so κ ≥ 3; put n = (κ − 1)/2 ≥ 1. The valuations of elements of R₀ are the even numbers and the odd numbers ≥ κ, so O₀/R₀ has length n, and O₀/πR₀ has length n + 1. This gives (1). If e is odd, the terms c_jπ^{2j} (j < e) of an element of R₀ have pairwise distinct valuations e·v(c_j) + 2j, so the least odd valuation is that of ℓ, namely e. If e = 2 and π² + a₁π + a₀ = 0, then R₀ = Z_ℓ + Z_ℓa₁π.

*Parity.* For j ≤ n:

> O₀ = R₀ + πR₀ + π^{2j}O₀,  (R₀ + π^{2j}O₀) ∩ (πR₀ + π^{2j}O₀) = π^{2j}O₀.  (P_j)

Indeed below 2j the elements of R₀ have even valuation and those of πR₀ odd valuation, so a sum r + πr′ lies in π^{2j}O₀ only if both terms do; and the two images in O₀/π^{2j}O₀ have length j each.

**Step 2: the shape of a lattice relative to a top.** For L in 𝓛 let S_k(L) ⊂ k′ be the image of L ∩ π^kO under π^kO → π^kO/π^{k+1}O ≅ k′, and d_k its dimension.

- τ = πω maps L ∩ π^kO into L ∩ π^{k+1}O, so ωS_k ⊂ S_{k+1} and d_k is non-decreasing.
- ⟨ , ⟩ induces a perfect pairing between π^kO/π^{k+1}O and π^{−1−k}O/π^{−k}O, a non-zero multiple of Re(ab̄). Since L is self-dual, the dual of L ∩ π^kO + π^{k+1}O is L ∩ π^{−1−k}O + π^{−k}O, so S_{−1−k} = S_k^⊥ and d_k + d_{−1−k} = 2.

Hence there is i ≥ 0 with d_k = 0 for k < −i, d_k = 1 for −i ≤ k ≤ i − 1 and d_k = 2 for k ≥ i. So π^iO ⊂ L ⊂ π^{−i}O, [π^{−i}O : L] = ℓ^{2i}, and i = 0 exactly when L = O. For i ≥ 1 the lines satisfy S_{k+1} = ωS_k; call S_{−i} = F_ℓȳ the *direction* of L relative to O.

In the same way, for an edge M (Lemma 2), S_k(M)^⊥ = ωS_{−2−k}(M), so d_k(M) + d_{−2−k}(M) = 2 and d_k(M) = 1 exactly for −1−j ≤ k ≤ −1+j, for some j ≥ 0, the *type* of M.

**Step 3: the top of a lattice.** *Claim: for every L in 𝓛 there is a unique top relative to which L = O or the direction of L is not isotropic.* In (a) there is one top and no isotropic line. In (b) write V = Ke ⊕ Kf with both summands isotropic and τ acting by π and −π. L ∩ Ke is stable under π, hence under Z_ℓ[π] = O_K, so L ∩ Ke = π^aO_Ke and L ∩ Kf = π^bO_Kf; by self-duality the projections of L are π^{−b}O_Ke and π^{−a}O_Kf. So L modulo (L ∩ Ke) ⊕ (L ∩ Kf) is the graph of an isomorphism φ from D = π^{−b}O_K/π^aO_K to π^{−a}O_K/π^bO_K, and β(z, w) = λ₀(zφ(w)) mod Z_ℓ is a perfect pairing on D. Integrality of L says β(z, w) + β(w, z) = 0. As ℓ is odd, β is alternating and non-degenerate, so the order ℓ^{a+b} of D is a square: a + b = 2i. (Here β(z, w) = ⟨ze, φ(w)f⟩, which is λ₀(zφ(w)) up to a constant.) Then O″ = π^{(a−b)/2}O_Ke ⊕ π^{(b−a)/2}O_Kf is a top with π^iO″ ⊂ L ⊂ π^{−i}O″, and L maps onto both factors of π^{−i}O″/π^{1−i}O″, so its direction is not an eigenline of ω. Uniqueness: if the direction relative to a top O″ is not isotropic, L contains x = π^{−i}y with y a unit of O″, so O·L = π^{−i}O″, which determines O″ and i.

Call this top top(L) and i the level; then O·L = π^{−i}top(L) and [O·L : L] = ℓ^{2i}. From now on O = top(L) when L is given, and x = π^{−i}y ∈ L with y ∈ O^×.

**Step 4: generators, the bound i ≤ n, and the ring.** Let i ≥ 1. The elements τ^kx = π^{k−i}ω^ky, 0 ≤ k ≤ 2i − 1, lie in L and have non-zero image in π^{k−i}O/π^{k−i+1}O. So L₁ = R·x + π^iO ⊂ L has d_k(L₁) ≥ d_k(L) for all k, hence L₁ = L.

If r ∈ R₀ had odd valuation 2j + 1 < 2i, then rx ∈ L would have image in F_ℓȳ at k = 2j + 1 − i, where S_k = F_ℓω^{2j+1}ȳ = F_ℓωȳ. These lines differ because ȳ is not an eigenvector of ω. So κ > 2i − 1, that is i ≤ n. This proves (2).

π^{2i}O·L ⊂ π^iO ⊂ L, so R + π^{2i}O ⊂ End(L). Conversely if a ∈ End(L) then a ∈ O, and ax = rx + π^io with r ∈ R and o ∈ O, so (a − r)y ∈ π^{2i}O and a ∈ R + π^{2i}O. The index: R + π^{2i}O = (R₀ + π^{2i}O₀) ⊕ (πR₀ + π^{2i}O₀)ω, and for i ≤ n each summand has index ℓ^i in O₀ by (P_i). This proves (3).

**Step 5: edges at a lattice.** Let M be an edge containing L, of level i. Since [M : L] = ℓ and d_k(M) ≥ d_k(L), exactly one d_k increases by one, and the shapes of Step 2 leave two possibilities: type i (the change is at k = −1 − i) or type i − 1 (the change is at k = i − 1, only for i ≥ 1).

*Ends of an edge.* Let M have type j relative to O. A unimodular L″ between τM and M has d_k(L″) = d_k(M) except at one k, so its interval of ones is [−j, j − 1] or [−1 − j, j]: its level relative to O is j or j + 1, and it is j exactly when L″ ⊃ π^jO, that is when L″ = π^jO + τM. The lattice π^jO + τM is indeed unimodular: its dual is π^{−j}O ∩ M, which contains it, and both have index ℓ in M (the image of π^jO in M/τM is a line because d_{j−1}(M) = 1, and M/(π^{−j}O ∩ M) is S_{−1−j}(M)). So an edge of type j has one end of relative level j and one of relative level j + 1.

*The edge going up.* Let i ≥ 1. An edge of type i − 1 at L contains π^{i−1}O, so it is M_up = L + π^{i−1}O, which has index ℓ over L. It is an edge: its dual is L ∩ π^{1−i}O; τM_up = τL + π^iO is contained in that; [L : L ∩ π^{1−i}O] = ℓ because d_{−i} = 1; and [L : τL + π^iO] = ℓ because [L : τL] = ℓ² and π^iO ⊄ τL (otherwise π^{i−1}O ⊂ L). So M_up^∨ = τM_up, ℓM_up ⊂ τM_up, and M_up/τM_up contains the isotropic line L/τM_up, so it is hyperbolic. Its other end is π^{i−1}O + τM_up = τL + π^{i−1}O, of relative level i − 1. Its direction is ω times that of L, not isotropic, so its top is O and its level is i − 1.

*Edges going down.* Let M = L + Z_ℓw be an edge of type i. Its other end has relative level i + 1 and direction the image of w. If i ≥ 1, then τw ∈ L has image in S_{−i}(L) = F_ℓȳ, so the direction is F_ℓω⁻¹ȳ, not isotropic: the end has top O and level i + 1. If i = 0, the ℓ + 1 lines of τ⁻¹O/O = π⁻¹O/O ≅ k′ are all isotropic for the form of note 7 because e ≥ 2 (Proposition J (4)). Those with non-isotropic direction lead to level 1 with top O. In (b) the two with isotropic direction are the lines stable under O, and their steps are the translations (note 7, Proposition J (2)). This proves (4), and (5) apart from the number of edges going down.

**Step 6: how many edges go down.** By Proposition L, L has ℓ + 1 isotropic eigenlines exactly when ρ = ℓ/τ² lies in End(L), and one otherwise. ρ lies in A₀, so by (3) the condition is ρ ∈ R₀ + π^{2i}O₀.

*For i ≤ n − 1 it holds.* By (P_{i+1}) write ρ ≡ r + πr′ mod π^{2i+2} with r, r′ in R₀. Then ρπ² = ℓ/δ lies in R₀, and ℓ/δ − rπ² ≡ π·r′π² mod π^{2i+2}. The left side is in R₀ and the right side in πR₀, so by (P_{i+1}) both are ≡ 0, and v(r′) ≥ 2i − 1. Valuations in R₀ below κ are even, so v(r′) ≥ 2i and ρ ≡ r mod π^{2i}.

*For i = n it fails.* Otherwise a lattice of level n would have ℓ + 1 isotropic eigenlines, hence by Step 5 an edge going down to a lattice of level n + 1, contradicting Step 4. Lattices of level n exist, by induction from the top using the previous paragraph and (4).

So a lattice of level 1 ≤ i < n has ℓ + 1 edges, one going up and ℓ going down; one of level n has only the edge going up. This completes (5).

**Step 7: stable lines.** Let i ≥ 1 and identify W(L) with τ⁻¹L/L. The line of the edge going up is (L + π^{i−1}O)/L, stable under End(L) = R + π^{2i}O. The element a = π^{2i}ω of End(L) maps τ⁻¹L into π^{i−1}O and kills π^{i−1}O modulo L. It does not kill τ⁻¹L/L: for w = τ⁻¹x = π^{−i−1}ω⁻¹y, the image of aw in π^{i−1}O/π^iO is F_ℓȳ (ω² is a scalar), which is not S_{i−1}(L) = F_ℓωȳ. So a acts on the plane as a non-zero nilpotent, whose only invariant line is its image. This proves (6).

**Step 8: the tree.** Every lattice of level i ≥ 1 has a unique edge going up, and its other end has level i − 1 and the same top. So the lattices with a given top form a rooted tree. Distinct edges going down have distinct lower ends, and the counts in (4) and (5) give (7). ∎

**Remarks.**

- The proof uses (U) through Lemma 1 (every L/τL is (Z/ℓ)²) and through v(τ) = 1 (π = τ/ω is a uniformiser).
- Step 1 says that n is the number of odd valuations missing from Z_ℓ[τ²]. Equivalently, π ↦ −π is a well-defined ring automorphism of O₀/π^{2j} exactly for j ≤ n; this is the content of (P_j). So the depth measures how far A looks like it has the extra symmetry τ ↦ −τ over A₀.
- In (b) with e = 1 and in (a) with e = 1 one has m = 2 and there is nothing below the top (note 7, Proposition J).

---

## 5. Check on the data [E]

`sixlines_levels.py` takes the pairs (polynomial, ℓ) of note 7 with a step going down at eigenvalue ±q: 9 pairs, all with ℓ = 3. For every lattice class L it computes [O_E·L : L] and End(L), lists the isotropic eigenlines and builds the neighbours, computes their levels directly, tests whether ℓL ⊂ (F ∓ q)²L at ℓ, and builds the lattice (F ∓ q)L + (F ∓ q)^{2i−1}O_E·L, which is τL + 𝔏^{i−1}top(L). It asserts, for every class:

1. [O_E·L : L] = ℓ^{2i};
2. [O_E : End(L)] = ℓ^{2i} at ℓ;
3. the levels are 0, …, n with 2n + 1 the exponent of ℓ in [O_E : Z[F, V]];
4. ℓ + 1 isotropic eigenlines if i < n and one if i = n, and ℓ + 1 exactly when ℓL ⊂ (F ∓ q)²L;
5. for i ≥ 1 exactly one neighbour has level i − 1, it equals the lattice of the formula, and all others have level i + 1;
6. from level 0: all neighbours at level 1 in (a); two at level 0 and ℓ − 1 at level 1 in (b).

All assertions pass: 9 pairs, 221 lattice classes, 410 neighbours.

| q | (e₁, e₂, e₃) | configuration | e | n | classes at levels 0, …, n |
|---|---|---|---|---|---|
| 37 | (66, 4551, 169756) | (a) | 3 | 1 | 1, 2 |
| 49 | (78, 5439, 240100) | (a) | 3 | 1 | 2, 8 |
| 13 | (18, 351, 4732) | (a) | 3 | 1 | 3, 6 |
| 19 | (12, 627, 5776) | (a) | 3 | 1 | 4, 16 |
| 19 | (24, 399, 5776) | (a) | 3 | 1 | 4, 16 |
| 17 | (22, 391, 5780) | (a) | 2 | 2 | 3, 6, 18 |
| 23 | (28, 391, 4232) | (a) | 2 | 1 | 6, 24 |
| 23 | (16, 667, 4232) | (b) | 2 | 1 | 10, 20 |
| 23 | (22, 391, 10580) | (b) | 2 | 1 | 24, 48 |

For e = 3 the theorem gives n = 1 whatever the polynomial, as observed. The numbers of classes are global: a top class has ℓ + 1 or ℓ − 1 lattices below it at level 1, which fall into fewer classes when the top has extra automorphisms (ratios 2 and 4 in (a), 2 in (b), then 3 from level 1 to level 2).

---

## 6. What this does not give

- Configurations (c) and (d) of Lemma 3, where no lattice is maximal at ℓ. Proposition L applies to them. In the data every lattice there has exactly one step: in (c) the lattices are paired and the step is not a translation, in (d) it is a translation (note 7 §5). This is not always so: by [note 9](09-ramified-over-e0.md) (Theorem N) the index of Z_ℓ[F] at these primes is ℓ^{2n}; for n = 1, the case of the data, there are two lattices exchanged by one step; for n ≥ 2 the lattices form two trees of depth n − 1 joined at their roots.
- Cases where (U) fails, so that ℓ⁴ divides P(±q) and Br(X)[ℓ^∞] is larger: F − q of valuation ≥ 2 at a prime, three or more primes, or a prime with f = 4. None occurs in the data. There W can have dimension 4 or 6, and one might expect products of volcanoes, not a single one. [Note 10](10-fourth-power.md) works these cases out. In rank 2 nothing changes. In rank 4 the expectation is right, for every power of ℓ: the graph is the product of two volcanoes, a step moving in both at once, or one volcano with branching ℓ² (Theorem P); the case of this note is the one where a factor is a single edge. In rank 6 it is wrong: the graphs are not products, and only a local rule for the number of steps is proved (Proposition O).
- The general statement that an eigenline stable under End(T) leads to a ring containing End(T) when the eigenvalue is ±q. It is proved here only in configurations (a) and (b); [note 9](09-ramified-over-e0.md) proves it in (c) and (d) as well (Theorem N (6)), so it holds whenever (U) holds and m ≥ 4.
- Classes, as opposed to lattices: how many classes of the genus-like sets of Proposition A lie at each level is not computed in closed form; the theorem gives the local count and the data the global one.
- The realisation by surfaces. The theorem is about lattices. Through Propositions H and G it predicts twisted equivalences defined over F_q along every edge; those propositions carry the unchecked steps listed in note 6 §8. Nothing is tested on equations.

---

## Reproduce

From `experiments/`, after `sixlines_vertical.py` (note 7):

```bash
python -I sixlines_levels.py --out sixlines_levels.jsonl > sixlines_levels.log
```

It takes a few seconds. An assertion fails if any of the six statements fails on any lattice class.

## Sources

As in notes [6](06-torsors-under-cm-class-groups.md) and [7](07-ring-changing-steps.md). No new source is used; the argument is elementary lattice theory.
