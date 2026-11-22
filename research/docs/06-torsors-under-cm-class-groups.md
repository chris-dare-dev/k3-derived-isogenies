# 6. Ordinary K3 surfaces over finite fields as torsors under CM class groups

*Assembles [notes 1](01-prior-art.md) to [5](05-non-split-count.md) into one statement with its evidence: ordinary K3 surfaces over a finite field form torsors under CM class groups (Proposition A), the steps are twisted derived equivalences over F̄_q (Propositions E and F), with their fields of definition (Proposition G). Proposition H and §6.7 summarise the ring-changing steps and fourfold computations of [notes 7](07-ring-changing-steps.md) to [16](16-fourfold-over-f25-and-f49.md). Its only new computation is the Artin–Tate check of §6.6. Labels: **[P]** read in the cited source's TeX; **[M]** recalled, not checked; **[D]** derived in this repository, with the argument given here or in the cited note; **[E]** measured by the scripts in `experiments/`. Nothing in it has been refereed. §8 lists what to check first.*

---

## 1. Summary

Taelman classified ordinary K3 surfaces over F_q by triples (M, F, 𝒦): a K3 lattice, a Frobenius endomorphism and an ample cone. This note reads that classification on the transcendental part.

1. **Structure (Proposition A).** For p ≥ 5 and geometric Picard rank ≥ 12, the surfaces with a fixed Néron–Severi part and a fixed rational Hermitian form on the transcendental part fall into genera. Each genus, with the Néron–Severi lattice marked, is a principal homogeneous space under a finite abelian group G_T, a quotient of the norm-one idèles of the CM field E = Q(F). Every element of G_T is the class of 𝔩/𝔩̄ for infinitely many degree-one primes 𝔩. Maximal and non-maximal orders are treated alike. The proof is Taelman's theorem and class field theory.
2. **A missing hypothesis (Remark B).** Conditions (M1)–(M5) of Taelman's Theorem B do not fix which real place carries the positive plane. The essential image is the triples that are positive at the place given by the embedding ι.
3. **A family where everything is computable (§5).** For double planes branched along six lines the transcendental lattice is L(2) with L unimodular of rank 6, the field E has degree 6, and the orders, the gluing and the count become explicit.
4. **Evidence (§6).** Four comparisons between enumerated surfaces and enumerated lattices, with no disagreement:

| comparison | notes | field sizes | polynomials | agree |
|---|---|---|---|---|
| number of surfaces, all lines rational | 2, 5 | p = 13 to 37 | 44 | 44 |
| the same, and Frobenius conjugation against the step at p | 3, 5 | 25, 49 | 32 | 32 |
| Frobenius action on NS of the surfaces fixed by conjugation | 4, 5 | 25, 49 | 9, with 29 surfaces | 29 of 29 |
| number of forms with lines permuted by an involution, by type | 5 | p = 5 to 23 | 169 | 169 |

5. **Derived equivalences (Propositions E and F).** Over F̄_q, any two surfaces of a genus are related by a single twisted derived equivalence, with Brauer classes of order ℓ, for infinitely many primes ℓ. The neighbour at a prime 𝔩 above ℓ is a moduli space of twisted sheaves of rank ℓ. This uses Bragg–Yang's construction for reflections, and the proof of their Torelli theorem, not only its statement.
6. **Rationality (Proposition G).** The partner and the isogeny are defined over F_q. The twisted equivalence is defined over F_{q^d}, with d the order of λ/q mod ℓ; over F_q itself exactly when 𝔩 divides F − q, which forces ℓ² to divide |Br(X)|. In the data Br(X) is usually trivial or a 2-group, and rational steps exist for 34 of the 490 pairs (polynomial, twist): 21 with a pair 𝔩 ≠ 𝔩̄ of degree-one primes dividing F ∓ q, and 13 more without one (§6.6).
7. **Steps that change the endomorphism ring (Proposition H, note 7).** All steps at an odd prime ℓ ≠ p are the isotropic lines of T/ℓT, and the rational ones are the eigenlines of F. The ring changes along a volcano at a double root of P mod ℓ, where the equivalence is never defined over F_q, and at ramified primes with eigenvalue ±q, where it is. In the data: 14 pairs (polynomial, ℓ) with ring-changing steps, 9 of them with a Brauer class of order 3 defined over F_q. When Br(X)[ℓ^∞] has order ℓ², the lattices below a maximal one form an ℓ-volcano whose depth and rings are explicit (note 8); when the primes are ramified over E₀ and no lattice is maximal, they form two such trees joined at their roots (note 9). When the order is ℓ⁴ or more and the eigenvalue has multiplicity 4, the graph is a product of two volcanoes or one volcano with branching ℓ² (note 10). With multiplicity 6 and order ℓ⁴ it is the graph of lattices over an explicit quartic order, and most often two lattices exchanged by one step (note 11); when −P(q)P(−q) is not a square in Q_ℓ the lattices are the self-dual ones for a hermitian form over the unramified quadratic extension (note 12).
8. **Not established (§7).** The primes 2 and p, the case where ℓ⁶ divides P(±q) and ±q has multiplicity 6 mod ℓ, a closed formula in every configuration when ℓ⁴ does, and any explicit step at a prime ℓ ≠ p on equations.

---

## 2. Setting

Fix a prime p ≥ 5, q = p^f, and an embedding ι: W(F_q) → C.

**Taelman's theorem.** [P] To an ordinary K3 surface X over F_q one attaches M = H²(X_can^ι, Z), with X_can the canonical lift, the endomorphism F of M matching geometric Frobenius on ℓ-adic cohomology, and the cone 𝒦 spanned by the ample classes. They satisfy:

- (M1) M is even unimodular of signature (3, 19);
- (M2) ⟨Fx, Fy⟩ = q²⟨x, y⟩;
- (M3) F is semisimple with eigenvalues of absolute value q;
- (M4) M ⊗ Z_p = M⁰ ⊕ M¹ ⊕ M² with F M^s = q^s M^s and ranks 1, 20, 1;
- (M5) 𝒦 is an F-stable connected component of the positive cone of NS(M, F) minus the walls of its (−2)-classes, where NS(M, F) = {x : F^d x = q^d x for some d}.

Theorem B of Taelman, *Ordinary K3 surfaces over a finite field*: X ↦ (M, F, 𝒦) is fully faithful, and essentially surjective if potential semistable reduction (⋆) holds. The theorem in his subsection "Unconditional results" gives an equivalence without (⋆) in three cases; the one used here is p ≥ 5 and rank NS ≥ 12.

**The transcendental part.** Write N = NS(M, F), T = N^⊥ and r = rank T = 22 − ρ̄. On N, F = q·g with g of finite order. In the proof of essential surjectivity Taelman shows [P]: T ⊗ Q is a one-dimensional vector space over E = Q[F|T], E is a CM field, complex conjugation sends F to V = q²/F, and the place v | p with v(F) > 0 has E_v = Q_p. So:

- the characteristic polynomial P of F on T is irreducible of degree r, and [E : Q] = r;
- choosing a generator, T becomes a lattice in E, stable under Z[F, V], with form ⟨x, y⟩ = Tr_{E/Q}(β x ȳ) for some β in the totally real subfield E₀;
- β is positive at exactly one real place of E₀ and negative at the others (signature (2, r − 2));
- by (M4), T ⊗ Z_p = T⁰ ⊕ T¹ ⊕ T², with T⁰ and T² of rank 1.

Let T¹ = {x ∈ E^× : x x̄ = 1}, as an algebraic torus over Q. It is the group of E-linear isometries of T ⊗ Q.

---

## 3. Proposition A: genera are torsors

**Marked surfaces.** Fix (N, g, 𝒦) as above. An *N-marked* ordinary K3 surface is a pair (X, m) with m: N → NS(X̄) an isometry carrying q·g to Frobenius and 𝒦 to the ample cone. Unmarked surfaces are marked ones modulo Aut(N, g, 𝒦).

**The action.** Let (M, F, 𝒦) be the triple of (X, m). For x in T¹(A_f), let x·M be the lattice in M ⊗ Q whose completion at each ℓ is (1 ⊕ x_ℓ)(M ⊗ Z_ℓ), with 1 on N and x_ℓ on T.

**The group.** Let U_T ⊂ T¹(A_f) be the elements with xT̂ = T̂ that act trivially on the discriminant group T^∨/T, and

> G_T = T¹(Q) \ T¹(A_f) / U_T.

**Proposition A.** [D, from Taelman's theorem] Let p ≥ 5 and let (X, m) be an N-marked ordinary K3 surface over F_q with ρ̄ ≥ 12.

1. For every x in T¹(A_f), (x·M, F, 𝒦) is the triple of an N-marked ordinary K3 surface x·X over F_q, unique up to isomorphism. It has the same zeta function as X.
2. x·X ≅ X as marked surfaces if and only if x ∈ T¹(Q)·U_T. So the orbit of (X, m) is a principal homogeneous space under the finite abelian group G_T.
3. Two marked surfaces (X, m), (X′, m′) lie in one orbit if and only if T(X) ⊗ Z_ℓ and T(X′) ⊗ Z_ℓ are isometric as Z_ℓ[F]-lattices for every ℓ, including p, compatibly with the gluing to N. (The rational Hermitian forms then agree up to a norm, by the Hasse norm theorem.) Call such a class a *genus*. The marked surfaces with given (N, g, 𝒦), given P and given rational Hermitian form are a finite disjoint union of genera.
4. G_T is a quotient of a ray class group of E, by s ↦ s/s̄. Hence every element of G_T is the class of the idèle that is a uniformiser at 𝔩, its inverse at 𝔩̄, and 1 elsewhere, for infinitely many degree-one primes 𝔩 ≠ 𝔩̄ of E. On lattices that element is T ↦ 𝔩 𝔩̄⁻¹ T when 𝔩 is prime to the conductor of End(T).

**Proof.** (1) x_ℓ is an isometry commuting with F, so x·M is even unimodular and F-stable, and (M2), (M3) hold. x_p commutes with the idempotents of Z_p[F] that cut out M⁰, M¹, M², so (M4) holds. N and 𝒦 are untouched, so (M5) holds and rank NS ≥ 12. The positive plane is the same subspace of M ⊗ R, so the positivity of Remark B holds. Taelman's unconditional theorem gives the surface; (M ⊗ Q, F) is unchanged, so is the zeta function. (2) An isomorphism of marked triples x·M → M is the identity on N and an E-linear isometry on T ⊗ Q, that is multiplication by some e in T¹(Q). It maps x·M to M exactly when ex preserves T̂ and the graph of the gluing, that is when ex ∈ U_T. The group is abelian, so every point of the orbit has the same stabiliser. T¹(Q) \ T¹(A_f) is compact and U_T is open, so G_T is finite. (3) If the local lattices are isometric compatibly with the gluing, the local isometries are elements x_ℓ of T¹(Q_ℓ), equal to units for almost all ℓ. The converse is clear. Finiteness: all such lattices lie between fixed bounds in T ⊗ Q. (4) By Hilbert 90 at each place, s ↦ s/s̄ maps the finite idèles of E onto T¹(A_f); it kills E^× modulo T¹(Q) and has open kernel modulo U_T. Chebotarev gives the primes. ∎

**Remarks.**

- *Unmarked surfaces.* Aut(N, g, 𝒦) acts on the marked surfaces through its image in the centraliser of g in O(q_N), by changing the gluing. In the six-line family the groups involved are explicit and give the counts of §5.
- *Orders.* End(T) contains Z[F, V]; at p it contains the three idempotents and F/q on T¹; and F acts on T^∨/T as g acts on N^∨/N. Different genera can have different endomorphism rings. A step at a prime not dividing the conductor preserves End(T). This is the "volcano" structure of the elliptic-curve case: the horizontal steps are in Proposition A, the vertical ones are passages between genera and are described in Proposition H.
- *Relation to the literature.* For End(T) = O_E and over C, the corresponding statements are Bayer-Fluckiger's theorem that K3 surfaces with CM by O_E and fixed discriminant ideal form a principal homogeneous space under her group 𝒞 (no dyadic place of E₀ ramified in E), and Valloni's K3 class group U_E(Q) \ U_E(A_f) / Ũ. [P] The statement over F_q, the congruence condition in U_T, and non-maximal orders are the additions. I found no source for them; that is weak evidence that there is none.
- *The prime p.* Elements with x_p = 1 are treated in Proposition E. The step with x_p = (p, 1, p⁻¹) is Frobenius conjugation (Proposition D); it is an isogeny and not a derived isogeny, by the slope condition in Bragg–Yang's existence theorem.

### Proposition E: the steps are derived isogenies over F̄_q

**Proposition E.** [D, from Bragg–Yang] Let p ≥ 5, X an ordinary K3 surface over F_q with ρ̄ ≥ 12, and x ∈ T¹(A_f) with x_p = 1. Then X ⊗ F̄_q and (x·X) ⊗ F̄_q are derived isogenous: there is a composition of twisted Fourier–Mukai equivalences between them whose action on ℓ-adic and crystalline cohomology is x. Since every class of G_T contains such an x (Proposition A (4)), any two surfaces of one genus are derived isogenous over F̄_q.

**Proof.** Write k = F̄_q, W = W(k), K = W[1/p], and fix an isomorphism K̄ ≅ C extending ι.

1. *A partner in characteristic p.* The lattice x·H²(X_k) has the same crystalline part as H²(X_k) and is unimodular at every ℓ. Bragg–Yang's existence theorem [P] gives a K3 surface X′ over k and a derived isogeny f: X_k ⇝ X′ with f*(H²(X′)) = x·H²(X_k); the slope condition holds because x_p = 1.
2. *Néron–Severi.* x is trivial on N ⊗ Q, and NS is saturated in ℓ-adic and in crystalline cohomology, so f* maps NS(X′) isomorphically onto N. Composing f with reflections in (−2)-curves of X′ and a sign, which are derived isogenies and integral, f* maps the ample cone of X′ to 𝒦.
3. *Lifting.* f is now a derived isogeny between surfaces of finite height that is integral on crystalline cohomology and matches Néron–Severi lattices and ample cones. These are the hypotheses used in the proof of Bragg–Yang's Torelli theorem up to its last step [P]: their theorem on lifting derived isogenies and their lemma on perfectly liftable isogenies give liftings over W of both surfaces, to which all line bundles lift, and an isogeny between the generic fibres that lifts f up to a change on Néron–Severi by an automorphism preserving the ample cone. In the proof of the lemma the lifting of the first surface is any lifting over W that lifts all line bundles; take the canonical lift X_can ⊗ W. The lifting X′_W of X′ is then defined by transporting the Hodge filtration through f. (Since f is an isomorphism of F-crystals it respects the slope decomposition, so X′_W is the canonical lift of X′; this is not needed.)
4. *Integrality.* Over C the lifted isogeny is a rational Hodge isometry H²(X′_W(C), Q) → H²(X_can(C), Q) = M ⊗ Q. By smooth and proper base change it carries H²(X′_W(C), Z_ℓ) onto x_ℓ(M ⊗ Z_ℓ) for ℓ ≠ p, and by the argument at the end of that proof (crystalline integrality gives p-adic étale integrality, Cais–Liu) onto M ⊗ Z_p. So it carries H²(X′_W(C), Z) onto x·M, up to an automorphism of (N, 𝒦), and the ample cone onto 𝒦.
5. *Conclusion.* (x·X)_can(C) has H² equal to x·M with the same Hodge structure and cone, by Taelman's construction. By the Torelli theorem X′_W(C) ≅ (x·X)_can(C). Both are generic fibres of smooth projective families over finite extensions of W, so by Matsusaka–Mumford the special fibres are isomorphic: X′ ≅ (x·X) ⊗ k. ∎

Step 3 uses a proof, with "Z-integral" weakened to "integral on crystalline cohomology and on Néron–Severi", and with a prescribed lifting of one surface. I read the proof and believe both changes are harmless, but this is the step to have checked (§8).

### Proposition F: one step is one twisted equivalence

**Proposition F.** [D, from Bragg–Yang] Let p ≥ 5 and X be ordinary over F_q with ρ̄ ≥ 12. Let ℓ ≠ p be a prime not dividing disc T nor the conductor of End(T), and 𝔩 ≠ 𝔩̄ a degree-one prime of E above ℓ. Let x be the idèle that is a uniformiser at 𝔩, its inverse at 𝔩̄ and 1 elsewhere, and Y = x·X. Then over k = F̄_q:

1. there are Brauer classes α on X_k and α′ on Y_k, both of order ℓ, and a Fourier–Mukai equivalence D^b(Y_k, α′) ≅ D^b(X_k, α);
2. α is the class whose kernel on T ⊗ Z_ℓ is 𝔩T;
3. Y_k is the moduli space of stable α-twisted sheaves on X_k with twisted Mukai vector (ℓ, 0, 0), taken with respect to the B-field b/ℓ of the proof and a generic polarisation.

By Proposition A (4) every class of G_T is the class of such an x, for infinitely many ℓ. So any two surfaces of one genus are twisted derived equivalent over F̄_q by a single equivalence.

**Proof.** Bragg–Yang's key geometric input for their existence theorem is a statement about reflections [P]: if b ∈ H²(X_k) = H²_ét(X_k, Ẑ^p) × H²_cris(X_k/W) is primitive, n = b²/2 is an integer, and b/n is a mixed B-field, then the moduli space X′ of stable α_{b/n}-twisted sheaves with Mukai vector (n, 0, 0) is a K3 surface, with a Brauer class α′ and an equivalence D^b(X′, α′) ≅ D^b(X_k, α_{b/n}) whose action f on cohomology satisfies f(H²(X′)) = s_b(H²(X_k)), where s_b is the reflection in b at every prime. It remains to choose b with s_b(H²(X_k)) = x·H²(X_k), and to identify X′.

*At ℓ.* Since ℓ does not divide the conductor, T ⊗ Z_ℓ is free of rank one over O_E ⊗ Z_ℓ and splits as T_𝔩 ⊕ T_𝔩̄ ⊕ T_rest according to the primes above ℓ. The summands T_𝔩 and T_𝔩̄ have rank one, each is isotropic, and both are orthogonal to T_rest, because ⟨ax, y⟩ = ⟨x, āy⟩. Since ℓ does not divide disc T, the pairing between them is perfect, and M ⊗ Z_ℓ = N_ℓ ⊕ T_ℓ. Choose generators e of T_𝔩 and f of T_𝔩̄ with ⟨e, f⟩ = 1 and put b_ℓ = ℓe + f. Then b_ℓ is primitive, b_ℓ² = 2ℓ, and

> s_b(e) = −f/ℓ,  s_b(f) = −ℓe,  s_b = 1 on the orthogonal complement of e and f.

So s_b(M ⊗ Z_ℓ) = ℓT_𝔩 ⊕ ℓ⁻¹T_𝔩̄ ⊕ T_rest ⊕ N_ℓ = x(M ⊗ Z_ℓ).

*Away from ℓ.* At each prime r ≠ ℓ, p, and in crystalline cohomology at p, take any primitive b_r with b_r² = 2ℓ; one exists because the lattice is unimodular and contains a hyperbolic plane. Since ℓ is a unit there, s_b preserves the lattice. So n = ℓ and s_b(H²(X_k)) = x·H²(X_k).

*B-field.* b/ℓ is integral away from ℓ. At p it lies in H²_cris(X_k/W), which consists of crystalline B-fields for surfaces of finite height [P]. So b/ℓ is a mixed B-field. Its Brauer class α is the ℓ-adic class of f/ℓ: it sends t ∈ T ⊗ Z_ℓ to ⟨f, t⟩/ℓ mod Z_ℓ, has order ℓ, and has kernel ℓT_𝔩 ⊕ T_𝔩̄ ⊕ T_rest = 𝔩T ⊗ Z_ℓ.

*Identification.* The equivalence is in particular a derived isogeny X′ ⇝ X_k with image x·H²(X_k). Steps 2 to 5 of the proof of Proposition E apply to it and give X′ ≅ Y_k. The class α′ is not trivial and has order ℓ: the equivalence identifies the ℓ-adic twisted transcendental lattices, the one of (X_k, α) is 𝔩T ⊗ Z_ℓ, of index ℓ in the unimodular lattice x(T ⊗ Z_ℓ) = T(X′) ⊗ Z_ℓ. ∎

**Remarks.**

- The proposition inherits the unchecked step of Proposition E (the identification), and nothing else beyond Bragg–Yang's statement for reflections.
- *Consistency over C.* T and 𝔩𝔩̄⁻¹T share the sublattice 𝔩T of index ℓ, which is the twisted transcendental lattice of (X_can, α) and of (Y_can, α′). For ρ ≥ 12 Huybrechts–Stellari's corollary [P] gives D^b(X_can(C), α) ≅ D^b(Y_can(C), α′). So the statement is the reduction mod p of a known one, as it should be.
- *The hypothesis on 𝔩* is the Elkies condition of note 1 §3.2: P has a simple root λ mod ℓ with λ ≢ q²/λ.
- *Rationality* is the subject of Proposition G.

### Proposition G: over which field the equivalence is defined

Keep the notation of Proposition F. Let λ ∈ F_ℓ^× be the image of F in O_E/𝔩, and d the order of λ/q in F_ℓ^×; so d divides ℓ − 1.

**Proposition G.** [D]

1. Y is defined over F_q, and so is the isogeny between X and Y as a correspondence with rational coefficients.
2. The subgroup generated by α in Br(X_k) is stable under Frobenius, which acts on it by λ/q or its inverse. So α is defined over F_{q^m} exactly when d divides m, and then it comes from a unique class in H²(X ⊗ F_{q^m}, μ_ℓ).
3. When d divides m, the equivalence of Proposition F is defined over F_{q^m}: there are classes α_m ∈ Br(X ⊗ F_{q^m}) and α′_m ∈ Br(Y ⊗ F_{q^m}) of order ℓ and an F_{q^m}-linear equivalence D^b(Y ⊗ F_{q^m}, α′_m) ≅ D^b(X ⊗ F_{q^m}, α_m). When d does not divide m there is none with this α.
4. In particular the equivalence is defined over F_q exactly when 𝔩 divides F − q. Then 𝔩̄ divides F − q as well, ℓ divides the index of Z[F, V] in O_E, and ℓ² divides |Br(X)|.
5. If ℓ does not divide the index of Z[F, V] in O_E, then d ≥ 3: the equivalence is not defined over F_q or F_{q²}.

**Proof.** (1) Y is defined over F_q by Proposition A. An isogeny is determined, up to equivalence, by its action on cohomology; here that action is x, which commutes with Frobenius. The cycle is defined over some F_{q^n}; the average of its n Frobenius conjugates is a Galois-invariant cycle with rational coefficients and the same action.

(2) α is the functional t ↦ ⟨f, t⟩/ℓ on T ⊗ Z_ℓ(1), with f a generator of T_𝔩̄(1). Geometric Frobenius acts on H²(X_k, Z_ℓ(1)) by F/q, and F acts on T_𝔩̄ by the image of V = q²/F at 𝔩, that is by q²/λ. So Frobenius multiplies f, and α, by q/λ. For the descent: Gal(F̄_q/F_q) has cohomological dimension 1 on torsion modules and H¹(X_k, μ_ℓ) = 0, so H²(X, μ_ℓ) → H²(X_k, μ_ℓ)^Γ is an isomorphism; and because ℓ does not divide disc N, H²(X_k, μ_ℓ) is the direct sum of NS/ℓ and T(1)/ℓ as a Galois module, so an invariant Brauer class has an invariant lift.

(3) Replace F_q by F_{q^m} and assume α invariant, with its class α_0 ∈ H²(X, μ_ℓ). The twisted Chern character, hence the condition "Mukai vector (ℓ, 0, 0)", depends only on the class in H²(μ_ℓ) [P], so the moduli problem of Proposition F is defined over F_q. Its coarse space M_0 is a K3 surface over F_q with M_0 ⊗ k = X′, the universal twisted sheaf is defined over F_q, and a Fourier–Mukai functor is an equivalence if it is one after extending the base field. So D^b(M_0, α′_0) ≅ D^b(X, α_0) over F_q. Its action f on cohomology commutes with Frobenius. The Hodge isometry of Proposition E, step 4, is f on ℓ-adic cohomology, so it commutes with F; after composing with the element of ±W(N) that moves the image of the ample cone of M_0 back to 𝒦, which commutes with g because both cones are Frobenius-stable, it is an isomorphism of triples. Taelman's functor is fully faithful over F_q, so M_0 ≅ Y over F_q. The last sentence of (3) holds because a class in Br(X ⊗ F_{q^m}) restricts to an invariant class.

(4) d = 1 means F ≡ q mod 𝔩. Conjugating, V ≡ q mod 𝔩̄, and F = q²/V ≡ q mod 𝔩̄. Then F and V have the same image (q, q) in O_E/𝔩 × O_E/𝔩̄, so Z[F, V] is not maximal at ℓ. Both functionals, for 𝔩 and for 𝔩̄, are invariant, so Br(X_k)^Γ contains (Z/ℓ)², and it is Br(X) at ℓ by (2). (5) If (λ/q)² = 1 then λ = q²/λ in F_ℓ, so F ≡ V mod 𝔩 and mod 𝔩̄, and the same argument gives that ℓ divides the index. ∎

**Remarks.**

- *The picture.* Over F̄_q every class of G_T is a single twisted equivalence, at infinitely many primes. Over F_q only finitely many primes can occur, those dividing F − q, and they all divide |Br(X)|. This is the situation of elliptic curves: an ℓ-isogeny with Galois-stable kernel is defined over F_q, but its kernel points are rational only over an extension of degree dividing ℓ − 1, and over F_q only if ℓ divides the number of points. Here the isogeny is the correspondence of (1), the kernel points are the Brauer class, and |Br(X)| plays the role of the number of points.
- *The Brauer group.* For ℓ prime to 2p and to disc N, (2) gives Br(X)[ℓ^∞] ≅ Hom(T ⊗ Z_ℓ / (F − q), Q_ℓ/Z_ℓ). Its order is the ℓ-part of P(q), in agreement with the Artin–Tate formula. If ℓ does not divide |Br(X)|, no twisted equivalence defined over F_q starts at X with a Brauer class of order divisible by ℓ.
- *Other invariant classes.* An invariant Brauer class of order ℓ whose kernel is not stable under O_E is treated in Proposition H. Such classes do not give partners with a smaller endomorphism ring when 𝔩 = 𝔩̄ is unramified with residue field F_{ℓ²}: then those classes are not isotropic and have no partner in the family (Proposition H (1); corrected in [note 7](07-ring-changing-steps.md) §6).

### Proposition H: all steps at ℓ, and the ones that change the ring

Proofs, the precise statements and the data are in [note 7](07-ring-changing-steps.md). Let ℓ be an odd prime not dividing q or disc T.

**Proposition H.** [D, note 7 Propositions H to K]

1. The surfaces of the family joined to X ⊗ F̄_q by one twisted equivalence with Brauer classes of order ℓ are the Kneser ℓ-neighbours of T, one for each isotropic line F_ℓ·x of T/ℓT; the Brauer class is that of x/ℓ. A class that is not isotropic has no partner in the family other than X.
2. The step is defined over F_q exactly when x is an eigenvector of F on T/ℓT. If the eigenvalue is λ, the equivalence is defined over F_{q^m} exactly when (λ/q)^m = 1.
3. The eigenlines with eigenvalue λ form a projective space of dimension μ − 1, with μ the number of generators of T over Z[F, V] at (ℓ, F − λ): one line where End(T) = Z[F, V], at least ℓ + 1 otherwise. For λ ≠ ±q all are isotropic. For λ = ±q, μ is even.
4. From a lattice maximal at ℓ, the steps that keep the ring are the translations 𝔩⁻¹𝔩̄T at degree-one primes; every other rational step goes down. At a double root λ ≠ ±q of P mod ℓ this is the ℓ-volcano of a quadratic order, and the equivalence is defined only over F_{q^d} with d ≥ 3. With eigenvalue ±q it happens at a pair 𝔩 ≠ 𝔩̄ of ramified primes (ℓ − 1 descents) and at a prime inert over E₀ and ramified over ℓ (ℓ + 1 descents); then the equivalence is defined over F_q for X or for its quadratic twist.
5. [D, note 8 Theorem M] In the two situations with eigenvalue ±q, suppose ℓ² is the exact power of ℓ dividing P(±q), and let ℓ^{2n+1} be the index of Z_ℓ[F] in the maximal order on the part where F ≡ ±q. Then n ≥ 1, and the lattices below a maximal one form a tree of depth n: a lattice L has a level i ≤ n with [O_E·L : L] = ℓ^{2i}; End(L) = Z_ℓ[F] + 𝔏^{2i}O_E, 𝔏 the product of the primes involved; from level i ≥ 1 one step goes up, along the only eigenline stable under End(L), and ℓ go down if i < n, none if i = n. Every edge is a twisted equivalence defined over F_q.
6. [D, note 9 Theorem N] Under the same hypothesis, suppose instead that the primes above ℓ with F ≡ ±q are ramified over E₀ (one of degree two, or two of degree one), with product 𝔏; then no lattice is maximal at ℓ. The index of Z_ℓ[F] at these primes is ℓ^{2n} with n ≥ 1. There is one O_E-lattice N with dual 𝔏⁻¹N, and every lattice L has a type j ≤ n − 1 with 𝔏^jN ⊂ L ⊂ 𝔏^{−1−j}N and End(L) = Z_ℓ[F] + 𝔏^{2j+1}O_E. There are two lattices of type 0, joined by a step, which is a translation exactly in the case of two primes. From type j ≥ 1 one step goes up, to (F ∓ q)L + 𝔏^{j−1}N; ℓ more go down if j < n − 1, none if j = n − 1. So the lattices form two trees of depth n − 1 joined at their roots; for n = 1 there are two lattices, exchanged by one step.
7. [D, note 10] Without the hypothesis on P(±q): the number of rational steps from a lattice L is determined by the sizes of the Jordan blocks of F ∓ q on L/ℓL (Proposition O), and L/(F ∓ q)L, hence the shape of Br(X)[ℓ^∞], can change along a step. If ±q has multiplicity 2 mod ℓ, every lattice is stable under the maximal order and the only steps are translations. If it has multiplicity 4, two numbers t₁, t₂ computed from the factor of P at ℓ decide: either they are in Q_ℓ and the graph is the product of two volcanoes of quadratic algebras, with depths read off t_i − 2, a step moving in both; or they are conjugate and the graph is one volcano with branching ℓ²; or no lattice exists (Theorem P). Items 5 and 6 in multiplicity 4 are the case where one factor is a single edge. For multiplicity 6 with ℓ⁴ dividing P(±q) there is no description.
8. [D, note 11] Multiplicity 6 with ℓ⁴ exactly dividing P(±q). The valuations of F ∓ q on the three pairs of eigenvalues are (1, ½, ½), (½, ¾, ¾), (3/2, ¼, ¼) or (⅔, ⅔, ⅔). In the pattern (3/2, ¼, ¼) there are exactly two lattices, exchanged by one step. In general a lattice is a leaf (one step), a middle lattice (ℓ + 1 steps) or a wide one (L/(F ∓ q)L ≅ (Z/ℓ)⁴); the wide lattices are the unimodular lattices stable under ℓ/τ, to which items 5 and 6 apply, and every leaf hangs on a lattice that is not a leaf. If −P(q)P(−q) is a square in Q_ℓ, the lattices are those of an explicit order Z_ℓ[x]/f of degree 4, up to scalars and of one parity, a step being a pair ℓΛ ⊂ Λ′ ⊂ Λ of index ℓ²; if not, the same holds over the unramified quadratic extension and the lattices are the ones fixed by its Galois group. This gives the graph by a computation for each polynomial, not by a formula.
9. [D, note 12 Theorem T] In item 8, suppose −P(q)P(−q) is not a square in Q_ℓ, and let k be the unramified quadratic extension. There is an explicit polynomial f₀ of degree 4 over Z_ℓ, f₀ ≡ X⁴ mod ℓ, such that with S₀ = Z_ℓ[X]/f₀ and B₀ = Q_ℓ[X]/f₀ the unimodular F-stable lattices are the O_k-lattices M of B₀ ⊗ k stable under S₀ with M^# = M or M^# = ℓ⁻¹M, the dual being taken for a hermitian form Tr(β v v̄′), β in B₀; a step joins M and M′ of the two kinds with ℓM ⊂ M′ ⊂ M. The form on T enters only through the class of β up to norms and scalars. Leaves exist for exactly one class, that of β = 1/f₀′(X); if no factor of f₀ has even residue degree the graph is finite and the number of leaves of one kind is ℓ^c·Π_𝔭(1 + 1/N𝔭)/(1 + 1/ℓ), ℓ^c the index of S₀ in the maximal order of B₀. The Galois action left open in item 8 is the map from the class of a lattice to the class of its dual.

So Proposition F is the case of one eigenline at a degree-one prime, and Proposition G (4) is the case of eigenvalue q.

---

## 4. Remark B: positivity

Let u be the eigenvalue of F on M⁰, a p-adic unit. Through ι it is a complex number, and the plane P_ι ⊂ M ⊗ R spanned by the eigenvectors of ι(u) and its conjugate is definite.

**Remark B.** [D] For a triple coming from a surface, P_ι is positive definite, because it is (H^{2,0} ⊕ H^{0,2}) ∩ H²(X_can^ι, R) by Taelman's Proposition on the Hodge decomposition. Conditions (M1)–(M5) do not imply this, because they do not mention ι. Suppose r > 2 and take the triple (M, F, 𝒦) of a surface with respect to ι. Choose another embedding ι′ = τ∘ι, with τ an automorphism of C that moves ι(u) to a conjugate over a different real place of E₀. The same triple satisfies (M1)–(M5), but T has only one positive plane, P_ι, so P_ι′ is negative definite, and the triple is not in the image of the functor defined with ι′. In the proof of essential surjectivity the condition is needed where the Torelli theorem is applied to the Hodge structure defined by (M4): the period map reaches only Hodge structures with ⟨ω, ω̄⟩ > 0. With the condition "P_ι is positive definite" added, the proof goes through unchanged.

This was read off arXiv:1711.09225v2, the latest arXiv version. I have not seen the published version. The data decide the matter independently: the number of surfaces equals the number of lattices positive at one place, not at any of the three (§6, first row; the three numbers are equal, so the sum would be three times too large).

---

## 5. The six-line family

**The lattices.** [D, note 2 §2; E, note 4 §4] Let X: w² = ℓ₁⋯ℓ₆ with the six lines in general position, resolved at its 15 nodes, over a field of odd characteristic. The 15 exceptional curves, the hyperplane class and the six line classes span a lattice N ⊂ NS(X̄) of rank 16 and discriminant group (Z/2)⁶. Its orthogonal complement is T₆ = L(2) with L odd unimodular of signature (2, 4). The group O(q_N) has order 1440; it is S₆ × Z/2, with S₆ the permutations of the lines and the central element z induced by association (the Gale transform), which gives the second six-line structure on the same surface.

If ρ̄ = 16, then E has degree 6 and E₀ is a totally real cubic field.

**Forms.** Over F_p consider six-line sets stable under Frobenius, with every line defined over F_{p²}. Frobenius permutes the lines by π with π² = 1: the identity (all lines rational, "split"), or one, two or three transpositions (types A, B, C).

**Proposition C.** [D, notes 2 §2 and 5 §2] Let X₀ be the double plane of a form of type π over F_p, ordinary with ρ̄ = 16, and F its Frobenius on L.

1. Carry the action of F on L/2L, with the quadratic form −x² mod 4, to (N^∨/N, q_N) by any isometry. The result lies in the conjugacy class of the permutation π in O(q_N). In particular (F − 1)/2 ∈ End(L) in the split case, and (F² − 1)/2, (F − V)/2 ∈ End(L) in general.
2. Conversely a unimodular lattice (L, F), with form Tr(α x ȳ) positive at the place of ι, with End(L) containing those elements and the part at p of §3, whose class in O(q_N) is that of π, comes from a form of type π.
3. The number of forms of type π with polynomial P(t) or P(−t) is the sum, over the isomorphism classes of such lattices, of the number of double cosets Z_{S₆}(π) \ Z_{O(q_N)}(π) / μ(End L). The summand is 2 when μ(End L) = ±1 (the form and its associate), and for π = 1 the count of surfaces is one per lattice.

The same holds over F_q with q² in place of p² for split configurations.

**Proposition D.** [D, notes 3 §2 and 4 §4] Let X be split over F_{p²}, and x the idèle that is p at the unit-root place v, p⁻¹ at v̄ and 1 elsewhere.

1. The Frobenius conjugate X^(p) is x·X. This follows from the main theorem of complex multiplication (Rizov; Valloni's form for ideal lattices) applied to the canonical lift.
2. X ≅ X^(p) if and only if F or −F has a square root ϖ in End(L), when μ(E) = ±1. Then ϖ is the Frobenius of a model over F_p, and its class in O(q_N) says how that Frobenius acts on NS: by a permutation of the lines, or by exchanging the two six-line structures.

---

## 6. Evidence [E]

All enumerations are exhaustive in the stated range. Surfaces are enumerated by normal forms of the configuration; Frobenius polynomials come from character sums over three extensions; lattices are enumerated with PARI as ideals with a Hermitian form together with self-dual sublattices for the non-maximal orders.

**6.1 Split count** (Proposition C for π = 1; notes 2 and 5 §6). Primes 7 ≤ p ≤ 37: 3,989 classes of configurations. For every ordinary irreducible polynomial the number of surfaces equals the number of lattices.

| p | 13 | 17 | 19 | 23 | 29 | 31 | 37 | total |
|---|---|---|---|---|---|---|---|---|
| polynomials | 1 | 2 | 2 | 4 | 9 | 10 | 16 | 44 |
| surfaces = lattices | 1 | 3 | 4 | 14 | 45 | 58 | 153 | 278 |

In 24 polynomials the order Z[F, V, (F − 1)/2] with its part at p is maximal. In the other 20 it is not, the lattices over all intermediate orders are needed, and in 4 of those there is no lattice over O_E at all. Also: (F − 1)/2 is integral in all 3,989 classes; the lattice count is the same for the three real places; counts per polynomial range from 1 to 26.

**6.2 Conjugation** (Proposition D.1; note 3). Over F₂₅ and F₄₉: 365 and 6,246 classes, 32 polynomials, 561 surfaces. The counts agree in all 32, and the cycle type of X ↦ X^(p) on the surfaces equals that of x on the lattices in all 32: 23 free involutions, 3 trivial, 6 mixed (for example 6 fixed and 4 swapped pairs).

**6.3 The Frobenius of the model over F_p** (Proposition D.2; note 4). A lattice is fixed by x exactly when ±F has a square root in End(L): all 561 lattice classes. For the 29 fixed surfaces the class of the root in O(q_N) equals the class read off from the configuration, polynomial by polynomial; five of the seven classes of involutions occur. The 217 forms over F₇ and 56 over F₅, enumerated separately, account for exactly the surfaces of "permutation" kind, and their zeta functions over F_p square to those over F_{p²}.

**6.4 Non-split count** (Proposition C for π ≠ 1; note 5). Primes 5 ≤ p ≤ 23: 51,804 classes of forms, 10,812 with ordinary irreducible polynomial, in 169 polynomials. The number of forms of each type A, B, C equals the prediction in all 169, on 5,468 lattices. In 16 polynomials Q(i) ⊂ E and the double-coset count is needed. In 9 polynomials the same lattice enumeration also contains lattices with F ≡ 1 mod 2; they predict the split configurations of 6.1 exactly (42 of 42).

**6.5 What the lattice side predicts beyond this** (note 3 §3). For maximal orders every step 𝔩/𝔩̄ at an Elkies prime ℓ ≤ 19 is a translation of the torsor; for example 26 surfaces at p = 37 form one 26-cycle under ℓ = 5. For non-maximal orders the steps preserve the level. None of this is tested on surfaces. By Proposition F each such step is, over F̄_p, one twisted derived equivalence with Brauer classes of order ℓ; for example the 12 surfaces with a common polynomial at p = 29 are joined in a cycle by equivalences twisted by classes of order 7.

**6.6 Artin–Tate and rational steps** (`sixlines_rational.py`). For a surface of the family with Frobenius q·g on N, the Artin–Tate formula reads |Br(X)| = 2^(16 − rank N^g)·P(q)/(q⁵·|disc N^g|), with P(−q) for the quadratic twist; N^g has rank 16, 12, 10, 10 and discriminant 64, 256, 1024, 1024 for g = 1, A, B, C. The right-hand side is the square of an integer in every case: 152 of 152 for split surfaces (76 polynomials, both twists) and 702 of 702 for the non-split forms (351 pairs of polynomial and type). This is an independent check of every Frobenius polynomial in notes 2, 3 and 5, since it constrains e₃ given e₁ and e₂. The values are 1 (506 times), 4 (222), 9 (89), 16 (27) and 25 (10).

So Br(X) is trivial or a 2-group in most cases, and then Proposition G allows no rational step at an odd prime. Of the 490 pairs (polynomial, twist), 48 have odd torsion in Br(X), always of order ℓ² with ℓ = 3 or 5. In 21 of them there is a pair 𝔩 ≠ 𝔩̄ of degree-one primes dividing F ∓ q (one split polynomial over F₄₉, twenty non-split ones), and the lattice classes are all maximal at ℓ. The step is then a permutation of the lattice classes that should be realised by twisted equivalences defined over F_q: for example, for the forms of type B over F₁₉ with polynomial (20, 399, 8664), twisted, the 38 lattice classes fall into two cycles of length 19 under the step at ℓ = 3. In the other 27 there is no such pair, but that does not rule out a rational step (6.7; [note 7](07-ring-changing-steps.md) §6): 14 of them have no rational step and 13 have one. As in 6.5, none of this is tested on surfaces.

**6.7 Rational steps and the endomorphism ring** (Proposition H; note 7, `sixlines_vertical.py`). For the 60 polynomials where an odd prime ℓ ≠ p divides the index of Z[F, V] in O_E (61 pairs of polynomial and ℓ), every eigenline of F on L/ℓL was listed for every lattice class, and the neighbour built for the isotropic ones: 5,398 eigenlines, 2,564 neighbours. Every neighbour is isomorphic to exactly one lattice class of the lists of 6.1 to 6.4, so those lists are closed under rational steps. 1,518 steps are translations, 610 go down, 412 go up, and 24 keep the ring without being translations.

| what happens at ℓ | pairs | ℓ |
|---|---|---|
| volcano at a double root λ ≠ ±q; equivalence over F_{q^d}, d from 3 to 12 | 5 | 7, 7, 7, 11, 13 |
| descent from a pair of ramified primes, Brauer class rational | 2 | 3 |
| descent from an inert ramified prime, Brauer class rational | 7 | 3 |
| a step between two genera with the same ring, Brauer class rational | 1 | 3 |
| translations only | 27 | 3, 5 |
| no rational step | 19 | 3, 5, 11 |

In each of the first three rows the counts per lattice are those of Proposition H (4): ℓ + 1 isotropic eigenlines at the top with 2 or 0 translations, and exactly one step going up from each lower lattice. One polynomial over F₁₇ has three levels, of index 1, 9 and 81 at ℓ = 3. In the nine pairs with a rational Brauer class, Br(X)[3^∞] is (Z/3)² at every level. None of this is tested on surfaces.

For those nine pairs `sixlines_levels.py` tests Proposition H (5) on every lattice class: the level, the ring, the depth, the number of isotropic eigenlines, the formula (F ∓ q)L + 𝔏^{i−1}top(L) for the lattice above, and the levels of all neighbours. All hold, on 221 classes and 410 neighbours.

Four pairs fall under Proposition H (6): the one of the row "a step between two genera" and three of those with translations only. All have n = 1. `sixlines_ramified.py` tests the statements on their 76 lattice classes. Since the data has no case with n ≥ 2, it also searches the Weil polynomials t³h(t + q²/t) with q = 5 and q = 7 for such configurations at ℓ = 3, 5 and 7, lists all lattices at ℓ for each in two independent ways, and tests the whole of Proposition H (6): 91,634 polynomials, 3,832 of them with n from 2 to 4, and 212,516 lattices. All hold. Those polynomials are not known to come from surfaces.

No pair of the data has ℓ⁴ dividing P(±q). `sixlines_deep.py` tests Proposition H (7) on searched polynomials with q = 5 and ℓ = 3: for multiplicity 4, the graph predicted from t₁, t₂ is compared with the lattices found, on 285 polynomials with ℓ⁴ or ℓ⁶ dividing P(q) and on 117 with ℓ² exactly; the rule for the number of steps is tested on every lattice met, 14,803 in all, including 1,295 in multiplicity 6. All agree. The graphs are compared by a colour-refinement invariant, and by one ball of radius 2 when they are infinite.

`sixlines_rank6.py` tests Proposition H (8) on the 235 searched polynomials of multiplicity 6 with ℓ⁴ exactly (q = 5, ℓ = 3): the patterns, the three kinds, the two lattices of the pattern (3/2, ¼, ¼), and, in the 102 cases where −P(q)P(−q) is a square and a starting lattice was found, the agreement with the quartic order, by the same kind of invariant. 18,119 lattices are met; all assertions pass. For 47 polynomials no starting lattice was found and only the pattern was checked.

`sixlines_hermitian.py` tests Proposition H (9) on the 105 of these polynomials for which −P(q)P(−q) is not a square, and on 1,025 more found for q = 7, ℓ = 3: the identities that define f₀, to 40 ℓ-adic digits; for each of the 1,200 classes of β the graph of self-dual lattices (14,334 lattices), with the kinds, the numbers of steps, the class that has leaves and the number of leaves; and, for the 927 polynomials where the method of `sixlines_rank6.py` has a starting lattice (1,779 forms), that the graph of unimodular F-stable lattices has the same invariant as the graph of a class of β. All assertions pass. The forms are not matched with the classes of β, only the graphs are.

`sixlines_groupsize.py` (note 13) works on polynomials only. For the seven primes of the first comparison of §1, the polynomials of the enumerated surfaces are exactly the ordinary irreducible Weil polynomials with P(q) = 64q⁵μ², P(−q) = 64q⁵ν² and (F − 1)/2 integral, 44 of 44. For a maximal order the group acting on the surfaces is Cl(E)/⟨Cl(E₀), ramified primes⟩; its order equals the number of lattices in all 191 cases counted, and is about 0.01·q² over 512 polynomials with q from 13 to 100,003 (class groups under GRH). And P is the second exterior power over Z[i] of the Frobenius polynomial of an abelian fourfold with multiplication by Z[i], on 23,989 polynomials.

`sixlines_prym.py` (note 14) builds that fourfold from the six lines over F_p, as the Prym variety of a curve of genus 5 that covers an elliptic curve cyclically of degree 4 (Paranjape's construction), and counts points of the curve over F_p, F_{p²}, F_{p³}. For the 534 enumerated configurations with sign + and irreducible P (p from 13 to 37; they make up 278 surfaces, [note 27](27-hypothesis-ledger.md) §1) and 4,983 ordered pairs of lines, the second exterior power of the characteristic polynomial of Frobenius on the i-eigenspace is P or P(−t), without exception; likewise for 12,873 pairs on 1,329 configurations with reducible P. For surfaces with sign − no comparison is made, and for them and for reducible P the form of the curve over F_p is not fully understood. `sixlines_prym_fq.py` (note 16) does the same over F₂₅ and F₄₉: 1,145 configurations with irreducible P (37 over F₂₅ and 1,108 over F₄₉, making up 21 and 564 surfaces, note 27 §1; 540 of the 564 are ordinary, and with the 21 over F₂₅ they are the 561 surfaces of §6.2, while the other 24 have an irreducible polynomial that is not ordinary), 4,182 pairs, no exception; 1,090 of these configurations are not defined over the prime field.

`sixlines_dictionary.py` (note 15) compares the class groups of the two fields, E of degree 6 for the surface and K of degree 8 for the fourfold, on 278 polynomials: the map of tori x ↦ (x_j x_k) gives a homomorphism Cl(K) → Cl(E) that matches the predicted image of each of 77,946 primes, and the induced map from Cl(K)/Cl(K⁺) to the group that acts on the surfaces has kernel and cokernel of order 1 to 8. It also checks, on 52,207 pairs of lines, that conjugating both lines of a pair under the involution of the cross-ratio map gives the associated configuration. Neither is a test of the action on surfaces at ℓ ≠ p; none exists.

---

## 7. What is not established

1. **Derived equivalences: what Propositions E and F do not give.** Bragg–Yang's existence theorem alone does not identify their partner X′ with the surface x·X of Taelman's theorem. Proposition E makes that identification by lifting, over F̄_q, and Proposition F shows that a step at an Elkies prime is a single equivalence. Still open: primes ℓ dividing disc T or the conductor (in the six-line family, ℓ = 2), and elements with x_p ≠ 1 on the slope-one part.
2. **Rationality: what Propositions G and H do not give.** They settle the steps of prime order ℓ, for ℓ odd and prime to p and to disc T, inside the family. Open: the 2-torsion and p-torsion of Br(X); steps at factors of P mod ℓ of degree 2, which change the ring through sublattices of index ℓ² and are not rational ℓ-steps; the graph when ℓ⁶ divides P(±q) and ±q has multiplicity 6 mod ℓ, and a closed description when ℓ⁴ does (ℓ² exactly is notes 8 and 9; multiplicity 2 and 4 are note 10; multiplicity 6 with ℓ⁴ is notes 11 and 12, as a reduction to a quartic order and, when −P(q)P(−q) is not a square, to self-dual hermitian lattices over it); the untwisted partners outside the family; and therefore the full class of surfaces joined to X by chains of twisted equivalences defined over F_q.
3. **An explicit step at ℓ ≠ p.** No twisted partner is exhibited. There is no modular polynomial for the family, explicit partners exist only at ℓ = 2, which divides the discriminant, and no invariant separates the surfaces of a torsor (note 3 §4). A possible route (note 13 §4): the polynomial P of a six-line surface is the second exterior power over Z[i] of the Frobenius polynomial of an abelian fourfold with multiplication by Z[i], the two moduli spaces coincide over C (Lombardo–Peters–Schütt), and an ℓ-step should be an isogeny of fourfolds with kernel (Z/ℓ)⁴, for which formulas exist. The correspondence over F_q is not established. Note 14 builds the fourfold from the six lines over F_p and finds the predicted Frobenius polynomial on every enumerated surface with irreducible P; the map from a fourfold back to six lines, which a step would need, is not made explicit.
4. **Geometry assumed in §5.** That a K3 surface over F_p with NS(X̄) ≅ N and the stated Frobenius is a six-line double plane; that N is primitive in H² for every member; that the nef cone is stable under the line permutations. The counts support all three.
5. **Scope of the data.** ρ̄ = 16 in one family; p ≤ 37; forms split over F_{p²} only; polynomials are taken from the surface side, so a polynomial with lattices and no surfaces would not be seen. The step tested on surfaces has order 2.
6. **Proposition D.1** is a derivation from the CM main theorem with the sign convention left open; for order 2 the sign does not matter.

---

## 8. What to check first

1. Remark B against the published version of Taelman's paper.
2. Proposition A (2): that the stabiliser needs triviality on T^∨/T and not only xT̂ = T̂. The six-line counts depend on it through the gluing (one surface per lattice when Aut(N) → O(q_N) is onto, two forms per lattice in the non-split case).
3. Proposition A (4) at primes dividing the conductor of End(T) or disc T, where the step is defined but is not multiplication by an ideal.
4. Proposition E, step 3: that the proof of Bragg–Yang's Torelli theorem applies with integrality assumed only on crystalline cohomology and on Néron–Severi, and with the canonical lift as the chosen lifting of the first surface. In particular their sentence that the modified isogeny differs from the given one on Néron–Severi only by an automorphism preserving the ample cone.
5. Proposition F: the choice of b, and that Bragg–Yang's statement for reflections applies to it (b/ℓ a mixed B-field with only an ℓ-adic Brauer class).
6. Proposition G (3): that the moduli space of twisted sheaves and its universal sheaf descend as stated, and the use of Taelman's full faithfulness to identify it with Y over F_q. Proposition H uses the same two steps for every isotropic eigenline; its part (1), second sentence, also uses that a twisted equivalence induces an isometry of ℓ-adic twisted transcendental lattices (note 7, Proposition H (4)).
7. The code path shared by every comparison in §6: `sixlines_cm.py`. Its sign test was wrong at low precision until note 5 §5; the outputs did not change when it was fixed, and the search for self-dual lattices was replaced and reproduces every earlier figure (note 5 §6).

---

## Reproduce

The figures of §6.1 to §6.5 come from the "Reproduce" sections of notes 2, 3, 4 and 5, and from the output files listed for them in [../DATA.md](../../data/README.md). For §6.6, from `experiments/`:

```bash
python -I sixlines_rational.py --out sixlines_rational.jsonl > sixlines_rational.log
```

It reads the result files of those notes and takes about three minutes. For §6.7:

```bash
python -I sixlines_vertical.py --out sixlines_vertical.jsonl > sixlines_vertical.log
python -I sixlines_levels.py --out sixlines_levels.jsonl > sixlines_levels.log
python -I sixlines_ramified.py --out sixlines_ramified.jsonl > sixlines_ramified.log
python -I sixlines_deep.py --out sixlines_deep.jsonl > sixlines_deep.log
python -I sixlines_rank6.py --out sixlines_rank6.jsonl > sixlines_rank6.log
python -I sixlines_hermitian.py --out sixlines_hermitian.jsonl > sixlines_hermitian.log
python -I sixlines_hermitian.py --search 7,3 --out sixlines_hermitian_q7.jsonl > sixlines_hermitian_q7.log
python -I sixlines_prym.py --workers 14 --out sixlines_prym.jsonl --run 5 7 11 13 17 19 23 29 31 37 > sixlines_prym.log
python -I sixlines_prym_fq.py --workers 14 --select all --pairs 0 --run 25 --out sixlines_prym_q25.jsonl > sixlines_prym_q25.log
python -I sixlines_prym_fq.py --workers 14 --select irreducible --pairs 3 --run 49 --out sixlines_prym_q49.jsonl > sixlines_prym_q49.log
python -I sixlines_dictionary.py dictionary --workers 14 --bound 2000 --out sixlines_dictionary.jsonl > sixlines_dictionary.log
python -I sixlines_dictionary.py association 5 7 11 13 17 19 23 29 31 37 > sixlines_dictionary_association.log
python -I sixlines_groupsize.py --cap 300 --out sixlines_groupsize.jsonl --fourfold 41 61 101 --real --all 13 17 19 23 25 29 31 37 41 43 47 49 53 59 61 67 71 73 79 83 89 97 101 --sample 211:40 503:40 1009:40 2003:30 5003:30 10007:30 30011:20 100003:12 125:20 625:20 3125:20 343:20 2401:20 > sixlines_groupsize.log
```

About three minutes for the first and a few seconds for the second.

## Sources

- Taelman, *Ordinary K3 surfaces over a finite field*, arXiv:1711.09225v2: https://arxiv.org/abs/1711.09225 (Theorems A and B; conditions (M1)–(M5); the proof of essential surjectivity; "Unconditional results") [P]
- Bayer-Fluckiger, *K3 surfaces with maximal complex multiplication*: https://arxiv.org/abs/2408.11643 (the group 𝒞; the corollary on surfaces with CM by O_E and discriminant ideal 𝔇) [P]
- Valloni, *Complex multiplication and Brauer groups of K3 surfaces*: https://arxiv.org/abs/1804.08763 (ideal lattices; K3 class group; main theorem in Shimura's form) [P]
- Rizov, *Complex multiplication for K3 surfaces*: https://arxiv.org/abs/math/0508018 [P]
- Bragg, Yang, *Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*: https://arxiv.org/abs/2102.01193 (introduction: existence theorem and Torelli theorem for derived isogenies; B-fields: ℓ-adic, crystalline and mixed, and the description of crystalline B-fields for finite height; §"Construction of Derived Isogenies": the lattice lemma and the proposition on reflections; §"Lifting derived isogenies to characteristic 0": perfect liftings, theorem on lifting; §"Twisted Derived Torelli Theorem": liftable isogenies, the lemma and the proof of the Torelli theorem) [P]
- Huybrechts, Stellari, *Equivalences of twisted K3 surfaces*: https://arxiv.org/abs/math/0409030 (corollary for ρ ≥ 12: equivalence exists if and only if the twisted transcendental lattices are Hodge isometric) [P]
- Huybrechts, *Motives of isogenous K3 surfaces*, Comment. Math. Helv. 94 (2019) [M]
- Matsumoto, Sasaki, Yoshida, Internat. J. Math. 3 (1992); Dolgachev, Ortland, *Point sets in projective spaces and theta functions* [M]
- Nikulin, *Integral symmetric bilinear forms and some of their applications* [M]
- Tate, *On the conjecture of Birch and Swinnerton-Dyer and a geometric analog*, Sém. Bourbaki 306 (1966); Milne, *On a conjecture of Artin and Tate*, Ann. of Math. 102 (1975) (the Artin–Tate formula); Liu, Lorenzini, Raynaud, *On the Brauer group of a surface*, Invent. Math. 159 (2005) (the order of Br is a square) [M]
