# 9. Eigenvalue ±q at primes ramified over E₀

*Treats the two configurations of primes, with no lattice maximal at ℓ, that [note 8](08-below-the-top.md) left open; builds on notes 6, 7 and 8. Labels: **[D]** derived here; **[E]** measured by `experiments/sixlines_ramified.py`. Not re-checked.*

---

## 1. Result

Notation is that of notes 6, 7 and 8. ℓ is an odd prime not dividing q or disc T, the eigenvalue is q (for −q replace F by −F), and τ = (F − q)/(F + q). Hypothesis (U) of note 8 is assumed: ℓ² is the exact power of ℓ dividing P(q). The multiplicity m of q as a root of P mod ℓ is at least 4.

The two configurations are those of note 8, Lemma 3, in which the primes of E above ℓ with F ≡ q are ramified over E₀:

- (c) one prime 𝔩 = 𝔩̄ of degree two;
- (d) two primes 𝔩₁, 𝔩₂ of degree one, each its own conjugate.

Write 𝔏 for 𝔩 in (c) and 𝔩₁𝔩₂ in (d). Then (Theorem N):

1. **One base lattice.** There is exactly one O_E-lattice N whose dual is 𝔏⁻¹N at these primes. Every unimodular F-stable lattice L lies between 𝔏^jN and 𝔏^{−1−j}N for one j ≥ 0, its *type*, and [O_E·L : L] = ℓ^{2j+1}.
2. **Depth.** The index of Z_ℓ[F] in the maximal order, at these primes, is ℓ^{2n} with n ≥ 1. The types that occur are 0, …, n − 1.
3. **Ring.** End(L) = Z_ℓ[F] + 𝔏^{2j+1}O_E at ℓ, of index ℓ^{2j+1} in the maximal order. It is never maximal and never Z_ℓ[F].
4. **Two tops.** There are exactly two lattices of type 0. They are joined by a step. In (d) this step is the translation at 𝔩₁ (equally at 𝔩₂). In (c) it is not a translation, and no isometry commuting with F maps one of the two lattices to the other.
5. **Steps.** From a lattice of type j ≥ 1 exactly one step goes up, to (F − q)L + 𝔏^{j−1}N. If j < n − 1 there are ℓ more steps from L, all going down to type j + 1; if j = n − 1 there are none. This holds for j = 0 as well, with the step of item 4 in place of the step going up.
6. **Shape.** The lattices form two rooted trees of depth n − 1 whose roots are joined by one edge. There are 2ℓ^j lattices of type j and 2(ℓ^n − 1)/(ℓ − 1) in all.
7. **Stable line.** At every lattice the eigenline of the step going up (across, for type 0) is the only one stable under End(L).

So the observation of note 8 §6, that every lattice in these configurations has exactly one step, is the case n = 1: two lattices exchanged by one step. For n ≥ 2 the graph is again a volcano, with two tops in place of one.

**Check [E].** In the four pairs (polynomial, ℓ) of the data in these configurations n = 1, and the statements hold on all 76 lattice classes. To test n ≥ 2 the script searches Weil polynomials of the same shape with q = 5 and q = 7: 91,634 polynomials, of which 3,832 have n ≥ 2, up to n = 4. On every one of them the complete list of lattices, obtained in two independent ways, agrees with the theorem (212,516 lattices in all). These polynomials are not known to come from surfaces of the family.

With note 8 this describes the rational steps at eigenvalue ±q in all four configurations of Lemma 3, that is whenever (U) holds and m ≥ 4.

---

## 2. Setting

As in note 8 §2: V is the part of T ⊗ Q_ℓ on which F ≡ q, A = Q_ℓ[F] acting on V, A₀ the subalgebra fixed by the involution, O and O₀ the maximal orders, R = Z_ℓ[F] = Z_ℓ[τ], and 𝓛 the set of unimodular τ-stable lattices in V. Lemmas 1 and 2 and Proposition L of note 8 are used as stated there.

In configurations (c) and (d), τ is skew and has valuation 1 at each prime (Lemma 3), and each prime is ramified over A₀. So

> A = A₀ ⊕ A₀τ,  O = O₀ ⊕ O₀τ,  ϖ := τ² generates the radical of O₀,  τO = 𝔏.

In (c) A is a field and A₀ has residue field F_{ℓ²}. In (d) A = A₁ × A₂ and A₀ = A₁,₀ × A₂,₀ with residue fields F_ℓ. Put

> k′ = O/τO = O₀/ϖO₀,

an F_ℓ-algebra of dimension 2: the field F_{ℓ²} in (c), and F_ℓ × F_ℓ in (d). The involution acts trivially on k′. For an element or a lattice, "level k" refers to the filtration by the τ^kO, and τ^kO/τ^{k+1}O is identified with k′ by dividing by τ^k; multiplication by τ is then the identity of k′.

Write Re(a₀ + a₁τ) = a₀. V is free of rank one over A; after choosing a generator, V = A and ⟨x, z⟩ = λ(Re(x z̄)) for a Q_ℓ-linear map λ: A₀ → Q_ℓ.

**Lemma 4 (the base lattice).** [D]

1. There is exactly one O-lattice N in V with N^∨ = τ⁻¹N.
2. Take V = A with N = O. Then λ(O₀) ⊂ Z_ℓ, and t(c) = λ(ϖ⁻¹c) mod Z_ℓ defines a non-zero F_ℓ-linear map t: k′ → ℓ⁻¹Z_ℓ/Z_ℓ ≅ F_ℓ with no non-zero ideal of k′ in its kernel. The pairing between τ^kO/τ^{k+1}O and τ^{−2−k}O/τ^{−1−k}O induced by ⟨ , ⟩ is, up to sign, β(a, b) = t(ab) on k′.
3. β is a non-degenerate symmetric form on k′. An isotropic vector of β is a unit of k′, and for an isotropic line S the elements a of k′ with aS ⊂ S are those of F_ℓ.

*Proof.* (1) The dual of an O-lattice is an O-lattice. In (c) the O-lattices are the τ^kN₁ for one of them, N₁, and (τ^kN₁)^∨ = τ^{c−k}N₁ where N₁^∨ = τ^cN₁. If c were even, τ^{c/2}N₁ would be a unimodular O-lattice, which does not exist at a prime ramified over E₀ (note 7, Proposition J (1); the reason is that the different of A over A₀ has odd valuation). So c is odd and k = (c + 1)/2 is the only solution of c − k = k − 1. In (d) the same holds at each of the two primes.

(2) Re(O·τ⁻¹O) = ϖ⁻¹Re(τO) = O₀, and N^∨ = τ⁻¹N says that {a ∈ A₀ : λ(aO₀) ⊂ Z_ℓ} = O₀. So t is well defined, and t(c·k′) = 0 forces ϖ⁻¹c ∈ O₀, that is c = 0 in k′. For a, b in O, ⟨τ^ka, τ^{−2−k}b⟩ = ±λ(ϖ⁻¹Re(ab̄)), and Re(ab̄) ≡ ab mod ϖ because the involution is trivial on k′.

(3) Non-degenerate by (2). In (d) the two lines that are ideals of k′ are spanned by the idempotents e_i, and β(e_i, e_i) = t(e_i) ≠ 0 by (2). So an isotropic vector s is a unit, and as ∈ F_ℓs gives a ∈ F_ℓ. ∎

From now on V = A and N = O.

**Lemma 5 (the number n).** [D] Let R₀ = Z_ℓ[ϖ] ⊂ O₀, and let G_k ⊂ k′ be the set of leading coefficients at level k of R₀ ∩ ϖ^kO₀, that is its image in ϖ^kO₀/ϖ^{k+1}O₀ ≅ k′. Then:

1. F_ℓ ⊂ G_k ⊂ G_{k+1}, and G₀ = F_ℓ. Let n ≥ 1 be the least k with G_k = k′.
2. [O₀ : R₀ + ϖ^kO₀] = ℓ^{min(k, n)}, R₀ ⊃ ϖ^nO₀, [O₀ : R₀] = ℓ^n, R = R₀ ⊕ R₀τ and [O : R] = ℓ^{2n}.
3. The leading coefficients of R at level k (for the filtration by τ^kO) form F_ℓ for k < 2n and k′ for k ≥ 2n.
4. ρ = ℓ/ϖ lies in R₀ + ϖ^kO₀ exactly when k ≤ n − 1.

*Proof.* (1) ϖ^k lies in R₀, multiplication by ϖ maps R₀ ∩ ϖ^kO₀ into R₀ ∩ ϖ^{k+1}O₀, and R₀ maps onto F_ℓ in k′ because ϖ maps to 0. R₀ has finite index, so G_k = k′ for large k.

(2) The length of O₀/(R₀ + ϖ^kO₀) is the number of i < k with G_i = F_ℓ. The rest follows, with R = Z_ℓ[τ] = R₀ + R₀τ.

(3) τ^kO = ϖ^{⌈k/2⌉}O₀ ⊕ ϖ^{⌊k/2⌋}O₀τ, so the leading coefficients of R at levels 2a and 2a + 1 are both G_a.

(4) Let I_i = R₀ ∩ ϖ^iO₀. For i < n, G_i = F_ℓ gives I_i = Z_ℓϖ^i + I_{i+1}. So I₁ ⊂ ϖR₀ + I_n. As ℓ ∈ I₁ and I_n ⊂ ϖ^nO₀, ρ lies in R₀ + ϖ^{n−1}O₀. Conversely R₀ + ϖ^nO₀ = R₀, and ρ is not in R₀: R₀ ≅ Z_ℓ[X]/(h) with h the characteristic polynomial of ϖ on A₀, and h(0) is ℓ² times a unit by (U); if ℓ = ϖ·r(ϖ) then ℓ − X·r(X) is a multiple of h, and X = 0 gives ℓ ∈ ℓ²Z_ℓ. ∎

**Two examples of n.** If m = 4 in (c), A₀ is the unramified quadratic extension, ϖ = ℓu with u a unit, and n = 1 + s where ℓ^s is the exact power of ℓ for which u is congruent to an element of Z_ℓ. If both primes have e = 2 in (d), A₀ = Q_ℓ × Q_ℓ, ϖ = (ℓu₁, ℓu₂), and n = 1 + v_ℓ(u₁ − u₂). So n = 1 is the general case and n ≥ 2 needs a congruence.

---

## 3. Theorem N

Assume (U) and configuration (c) or (d). Let N, k′, β and n be as in Lemmas 4 and 5.

**Theorem N.** [D]

1. *Existence.* 𝓛 is not empty exactly when β has an isotropic vector. It then has two isotropic lines, S₊ and S₋.
2. *Type and direction.* For L in 𝓛 there is j with 0 ≤ j ≤ n − 1, the type, such that τ^jN ⊂ L ⊂ τ^{−1−j}N, both of index ℓ^{2j+1}, and O·L = τ^{−1−j}N. The image of L in τ^{−1−j}N/τ^{−j}N ≅ k′ is S₊ or S₋, the *direction* of L. For every x in L that generates τ^{−1−j}N over O, L = R·x + τ^jN.
3. *Ring.* End(L) = R + τ^{2j+1}O, of index ℓ^{2j+1} in O.
4. *Type 0.* The lattices of type 0 are L₊ and L₋, the inverse images of S₊ and S₋ in τ⁻¹N. They are joined by a step. In (d), L₋ = εL₊ for the isometry ε = (−1, 1) of V = A₁ × A₂, and this is the translation of note 7 at 𝔩₁. In (c) no element u of A with uū = 1 maps L₊ to L₋.
5. *Steps.* From L of type j ≥ 1 exactly one step leads to a lattice of type j − 1; it is τL + τ^{j−1}N. Every other step from a lattice of type j ≥ 0 leads to a lattice of type j + 1 with the same direction, except the step between L₊ and L₋. There are ℓ of them if j < n − 1 and none if j = n − 1.
6. *Stable lines.* The eigenline of the step to τL + τ^{j−1}N (to the other top, if j = 0) is the only line of W(L) stable under End(L).
7. *Shape.* The lattices with direction S₊ form a rooted tree under the steps, with root L₊ and ℓ^j vertices of type j, 0 ≤ j ≤ n − 1; the same for S₋; and the only other edge joins L₊ and L₋.

**Consequences.**

- For n = 1 the set 𝓛 has two elements, exchanged by a step, and each has exactly one rational step with eigenvalue q. This is what the data shows.
- A lattice of type n − 1 has one isotropic eigenline and ℓ anisotropic ones; a lattice of lower type has ℓ + 1 isotropic ones.
- In all four configurations of Lemma 3, End(L) is never Z_ℓ[F]. The indices in O of the rings that occur are 1, ℓ², …, ℓ^{2n} in (a) and (b), where [O : R] = ℓ^{2n+1}, and ℓ, ℓ³, …, ℓ^{2n−1} in (c) and (d), where [O : R] = ℓ^{2n}.
- The statement left open in note 7 after Proposition I now holds whenever (U) holds and m ≥ 4: the eigenline stable under End(L) leads to a lattice whose ring contains End(L), and the other isotropic eigenlines lead to smaller rings. In (c) and (d) at type 0 the ring is the same on both sides.
- By Proposition H (2) and (3), every edge is realised by a twisted equivalence defined over F_q, as in note 8.

---

## 4. Proof of Theorem N

**Step 1: the shape of a lattice.** For L in 𝓛 let S_k(L) ⊂ k′ be the set of leading coefficients of L ∩ τ^kO at level k, and d_k its dimension.

- τ maps L ∩ τ^kO into L ∩ τ^{k+1}O, so S_k ⊂ S_{k+1}.
- L is self-dual and (τ^kO)^∨ = τ^{−1−k}O, so the dual of L ∩ τ^kO + τ^{k+1}O is L ∩ τ^{−2−k}O + τ^{−1−k}O. By Lemma 4 (2), S_{−2−k} = S_k^⊥ for β, and d_k + d_{−2−k} = 2.

So d_{−1} = 1 and S_{−1} is its own orthogonal: an isotropic line S. As d_k is non-decreasing, there is j ≥ 0 with d_k = 0 for k < −1 − j, d_k = 1 for −1 − j ≤ k ≤ j − 1, and d_k = 2 for k ≥ j; and all the lines S_k are equal to S, because they are nested. Hence τ^jO ⊂ L ⊂ τ^{−1−j}O with both indices ℓ^{2j+1}. L contains an element x = τ^{−1−j}y with ȳ spanning S; by Lemma 4 (3) y is a unit of O, so O·L = τ^{−1−j}O. In particular 𝓛 is empty if β is anisotropic.

**Step 2: generators, the bound j ≤ n − 1, and the ring.** With x as above, τ^ix for 0 ≤ i ≤ 2j lies in L and has non-zero leading coefficient at level i − 1 − j. So L₁ = R·x + τ^jO ⊂ L has d_k(L₁) ≥ d_k(L) for all k, and L₁ = L.

By Lemma 5 there is r in R₀ of level 2n (for the filtration of O) with leading coefficient c outside F_ℓ. Then rx ∈ L has level 2n − 1 − j and leading coefficient cȳ, which is not in S by Lemma 4 (3). If j ≥ n, then −1 − j ≤ 2n − 1 − j ≤ j − 1 and S_{2n−1−j}(L) = S, a contradiction. So j ≤ n − 1. This proves (2), given Step 1.

τ^{2j+1}O·L ⊂ τ^jO ⊂ L, so R + τ^{2j+1}O ⊂ End(L). Conversely let a be in End(L) ⊂ O. Then ax = rx + τ^jo with r ∈ R and o ∈ O, so (a − r)y ∈ τ^{2j+1}O and a ∈ R + τ^{2j+1}O. Finally R + τ^{2j+1}O = (R₀ + ϖ^{j+1}O₀) ⊕ (R₀ + ϖ^jO₀)τ, of index ℓ^{j+1}·ℓ^j in O by Lemma 5 (2), since j + 1 ≤ n. This proves (3).

**Step 3: type 0.** A lattice L with O ⊂ L ⊂ τ⁻¹O is τ-stable, since τL ⊂ O. It is unimodular exactly when L/O is an isotropic line of τ⁻¹O/O = N^∨/N for the form with values in Q_ℓ/Z_ℓ, which by Lemma 4 (2) is β up to sign. So the lattices of type 0 correspond to the isotropic lines of β: none or two. This proves (1). If there are two, M = τ⁻¹O satisfies M^∨ = τM and ℓM ⊂ τM, so it is an edge (note 8, Lemma 2) and L₊, L₋ are joined by a step.

In (d), ε = (−1, 1) lies in O, εε̄ = 1, and ε commutes with τ; so ε is an isometry of V that preserves 𝓛 and N. It acts on k′ = F_ℓ × F_ℓ by (−1, 1). The form is β((a₁, a₂), (b₁, b₂)) = γ₁a₁b₁ + γ₂a₂b₂ with γ_i ≠ 0, its isotropic lines are spanned by (1, ±s) with s² = −γ₁/γ₂, and ε exchanges them. So εL₊ = L₋. If a is a uniformiser at 𝔩₁ and a unit at 𝔩₂, then a/ā ≡ ε mod τO, and 1 + τO preserves L₊ (it maps L₊ into L₊ + τ·τ⁻¹O). So (a/ā)L₊ = L₋: the step is the translation at 𝔩₁, and likewise at 𝔩₂ since −ε has the same effect.

In (c), an element u of A with uū = 1 is a unit, and its image in k′ = F_{ℓ²} has square 1 because the involution is trivial on k′. So u ∈ ±1 + τO and uL₊ = L₊. This proves (4). Since the automorphisms of V that commute with F and preserve the form are these u, L₊ and L₋ are not isomorphic as lattices with Frobenius and form.

**Step 4: edges.** Let M be an edge: τ-stable with M^∨ = τM. The dual of M ∩ τ^kO + τ^{k+1}O is τM ∩ τ^{−2−k}O + τ^{−1−k}O = τ(M ∩ τ^{−3−k}O + τ^{−2−k}O). So S_{−3−k}(M) = S_k(M)^⊥ and d_k(M) + d_{−3−k}(M) = 2: there is i ≥ 0, the *type* of M, with d_k(M) = 1 exactly for −1 − i ≤ k ≤ i − 2. The edge of type 0 is τ⁻¹O.

*Edges at a lattice.* Let M be an edge containing L of type j. Since [M : L] = ℓ and d_k(M) ≥ d_k(L), exactly one d_k increases by one. The shapes allow two cases: M has type j (the change is at k = j − 1), and then M ⊃ τ^{j−1}O, so M = L + τ^{j−1}O; or M has type j + 1 (the change is at k = −2 − j).

*Ends of an edge.* Let M be an edge of type i ≥ 1 and L″ a unimodular lattice between τM and M. Then d_k(L″) = d_k(M) except at one k, so L″ has type i or i − 1, and type i − 1 exactly when L″ ⊃ τ^{i−1}O, that is when L″ = τM + τ^{i−1}O. That lattice is unimodular: its dual is M ∩ τ^{−i}O, which contains it, and both have index ℓ in M (the first because d_{i−2}(M) = 1, the second because d_{−1−i}(M) = 1). The plane M/τM then has an isotropic line, so it has exactly two. So an edge of type i ≥ 1 joins one lattice of type i − 1 to one of type i. Both have the direction S_{−1}(M), a line because i ≥ 1.

*The edge going up.* Let L have type j ≥ 1 and M_up = L + τ^{j−1}O, of index ℓ over L. Its dual is L ∩ τ^{−j}O, of index ℓ in L because d_{−1−j}(L) = 1. τM_up = τL + τ^jO is contained in it, and has index ℓ in L: [L : τL] = ℓ² by (U), τL + τ^jO ≠ τL because τ^{j−1}O ⊄ L, and τL + τ^jO ≠ L because it lies in L ∩ τ^{−j}O. So M_up^∨ = τM_up; also ℓM_up ⊂ τM_up. So M_up is an edge of type j, and its other end is τM_up + τ^{j−1}O = τL + τ^{j−1}O, of type j − 1.

So the steps from L are: the step along the edge of type j, which exists and is unique (for j = 0 it is the step between L₊ and L₋); and steps along edges of type j + 1, leading to lattices of type j + 1 with the direction of L. This proves (5) apart from the number of edges of type j + 1.

**Step 5: how many edges go down.** By Proposition L, L has ℓ + 1 isotropic eigenlines exactly when ρ = ℓ/τ² lies in End(L), and one otherwise (m ≥ 4 in these configurations, since each prime is ramified). ρ lies in A₀, so by (3) the condition is ρ ∈ R₀ + ϖ^{j+1}O₀, which by Lemma 5 (4) holds exactly when j + 1 ≤ n − 1. So L has ℓ edges of type j + 1 if j < n − 1 and none if j = n − 1. This completes (5).

**Step 6: stable lines.** Identify W(L) with τ⁻¹L/L and let L have type j ≥ 0 and x = τ^{−1−j}y as in Step 1. The line of the edge of type j is (L + τ^{j−1}O)/L; it is stable under End(L) = R + τ^{2j+1}O because τ^{2j+1}O·τ^{j−1}O ⊂ L. Choose c in O with image in k′ outside F_ℓ and put a = τ^{2j+1}c ∈ End(L). Then a maps τ⁻¹L into τ^{j−1}O and τ^{j−1}O into τ^jO ⊂ L. It does not kill τ⁻¹L/L: a·τ⁻¹x = τ^{j−1}cy has leading coefficient cȳ at level j − 1, which is not in S_{j−1}(L) = S. So a acts on the plane W(L) as a non-zero nilpotent, and its only invariant line is its image. This proves (6).

**Step 7: the trees.** By Step 4 every lattice of type j ≥ 1 has exactly one neighbour of type j − 1, with the same direction, and distinct edges going down from a lattice have distinct lower ends. By Step 3 there is one lattice of type 0 for each direction, and by Step 5 a lattice of type j < n − 1 has ℓ neighbours of type j + 1. This gives (7) and the counts. ∎

**Remarks.**

- The argument is the one of note 8 with two changes. The filtration is by powers of τ itself, on which τ acts as the identity of k′, so the lines S_k of a lattice are all equal, where in note 8 they alternate between a line and its image under ω. And duality is centred at −1, which forces a middle line that is isotropic; this is why there are two tops and no lattice stable under O.
- (U) enters through Lemma 1 and Proposition L of note 8, through v(τ) = 1, and through h(0) in Lemma 5 (4).
- In (c) the two trees consist of lattices that are not isomorphic to each other at ℓ. Globally the two tops lie in different genera of lattices with Frobenius; note 7 recorded the step between them as a step between genera with the same ring.

---

## 5. Check [E]

`sixlines_ramified.py` runs in two steps.

**Step A, the data.** The pairs (polynomial, ℓ) of note 7 in configurations (c) and (d): four pairs, all with ℓ = 3, eigenvalue −q and n = 1. For every lattice class L the script asserts: [O_E·L : L] = ℓ; End(L) = Z_ℓ[F] + 𝔏O_E at ℓ; one isotropic eigenline among the ℓ + 1, and ℓL ⊄ (F + q)²L; that line is stable under End(L) and no other is; the neighbour L′ satisfies L + L′ = O_E·L and L ∩ L′ = 𝔏·O_E·L at ℓ; and L′ is a translate of L exactly in (d).

| q | (e₁, e₂, e₃) | configuration | e of the primes | n | classes |
|---|---|---|---|---|---|
| 23 | (34, 391, −2116) | (c) | 2 | 1 | 24 |
| 49 | (66, 5439, 297724) | (d) | 2, 4 | 1 | 4 |
| 19 | (18, 399, 10108) | (d) | 2, 4 | 1 | 36 |
| 19 | (30, 399, 1444) | (d) | 2, 4 | 1 | 12 |

All assertions pass on the 76 classes.

**Step B, deeper examples.** The data has no case with n ≥ 2, so the script looks for some. It lists the polynomials t³h(t + q²/t) with h a monic integer cubic with three real roots in (−2q, 2q), such that ±q is a root of multiplicity ≥ 4 mod ℓ and ℓ² exactly divides P(±q); keeps those that are irreducible, in configuration (c) or (d), with no other prime above ℓ ramified over E₀; and for each one:

- chooses a form ⟨x, y⟩ = Tr(αx ȳ) on O_E, α in E₀, for which N = O_E has dual 𝔏⁻¹ at ℓ and a unimodular lattice exists;
- lists the lattices by following the steps from one lattice of type 0;
- lists them again without using steps: every F-stable lattice between 𝔏^n and 𝔏^{−1−n} of the right index is reached by descending through maximal F-stable sublattices, and the unimodular ones are kept;
- asserts that the two lists are equal, and then statements (1) to (7) of the docstring on every lattice: position and index, the ring as an order, the numbers 2ℓ^j, the number of isotropic eigenlines and the criterion ℓL ⊂ (F ∓ q)²L, the formula for the lattice above, the types of all neighbours, the single stable line, which neighbours are translates, and that the graph is a tree.

| q | ℓ | configuration | n | polynomials | lattices for each |
|---|---|---|---|---|---|
| 5 | 3 | (c) | 1 | 3,304 | 2 |
| 5 | 3 | (c) | 2 | 238 | 8 |
| 5 | 3 | (c) | 3 | 12 | 26 |
| 5 | 3 | (d) | 1 | 3,974 | 2 |
| 5 | 3 | (d) | 2 | 90 | 8 |
| 5 | 7 | (c) | 1 | 538 | 2 |
| 5 | 7 | (c) | 2 | 6 | 16 |
| 5 | 7 | (d) | 1 | 366 | 2 |
| 7 | 3 | (c) | 1 | 26,390 | 2 |
| 7 | 3 | (c) | 2 | 1,844 | 8 |
| 7 | 3 | (c) | 3 | 150 | 26 |
| 7 | 3 | (c) | 4 | 6 | 80 |
| 7 | 3 | (d) | 1 | 33,688 | 2 |
| 7 | 3 | (d) | 2 | 1,100 | 8 |
| 7 | 3 | (d) | 3 | 46 | 26 |
| 7 | 3 | (d) | 4 | 4 | 80 |
| 7 | 5 | (c) | 1 | 10,254 | 2 |
| 7 | 5 | (c) | 2 | 230 | 12 |
| 7 | 5 | (c) | 3 | 8 | 62 |
| 7 | 5 | (d) | 1 | 9,288 | 2 |
| 7 | 5 | (d) | 2 | 98 | 12 |

In total 91,634 polynomials and 212,516 lattices: 8,528 and 19,396 for q = 5, 83,106 and 193,120 for q = 7.

All assertions pass. The second list goes one type beyond n − 1, so it also confirms that no lattice of type n exists. The polynomials come in pairs P(t), P(−t), which exchange the eigenvalues q and −q; both are tested.

---

## 6. What this does not give

- Cases where (U) fails (ℓ⁴ divides P(±q)), as in note 8 §6.
- The case m = 2 under (U), where note 8, Proposition L, gives 0 or 2 steps from each lattice; it is not revisited here. This case and the previous one are treated in [note 10](10-fourth-power.md): m = 2 for every power of ℓ (every lattice is stable under the maximal order), m = 4 for every power of ℓ (Theorem P, which contains the rank-4 cases of this note), and m = 6 with ℓ⁴ only through a local rule and computations.
- A closed formula for the number of classes, as opposed to lattices, of each type.
- Surfaces. The polynomials of Step B are Weil polynomials of the right shape, but I have not checked the remaining conditions for them to be Frobenius polynomials of K3 surfaces, and nothing is claimed about the six-line family for them. The four pairs of Step A are from the family, and there n = 1.
- The realisation of the steps by twisted equivalences rests on Propositions G and H, with the unchecked steps listed in note 6 §8. Nothing is tested on equations.

---

## Reproduce

From `experiments/`, after `sixlines_vertical.py` (note 7):

```bash
python -I sixlines_ramified.py --out sixlines_ramified.jsonl > sixlines_ramified.log
```

It takes about two minutes (Step A, then q = 5 with ℓ = 3 and ℓ = 7). The rows with q = 7 come from

```bash
python -I sixlines_ramified.py --search 7,3 --search 7,5 --out sixlines_ramified_q7.jsonl > sixlines_ramified_q7.log
```

which repeats Step A and takes about half an hour. An assertion fails if any statement fails on any lattice.

## Sources

As in notes [6](06-torsors-under-cm-class-groups.md), [7](07-ring-changing-steps.md) and [8](08-below-the-top.md). No new source is used.
