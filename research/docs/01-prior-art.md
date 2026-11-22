# 1. Prior art: what is already known

*Checks what the literature already gives for the class-group action on ordinary K3 surfaces over a finite field and its realisation by twisted derived equivalences, records which statements about the sources need care, and points to the notes where the new work and computation are done. Labels: **[P]** read in the cited source's TeX; **[M]** recalled or taken from an abstract, not checked; **[D]** derived here, with the argument given; **[H]** heuristic.*

---

## 1. Summary

The starting point is hypothesis (a) of [../SUMMARY.md](../../SUMMARY.md) §1.1: by analogy with CSIDH, twisted derived isogenies of ordinary K3 surfaces over F_q should form a class-group action [H]. Made precise, it consists of two statements and three computations:

- **Statement 1 (class-group action).** The ordinary K3 surfaces over F_q whose Taelman triple has the same Frobenius polynomial, the same NS-genus and the same local data as a given surface X form a torsor under a finite abelian group G, a quotient of a (norm-one) class group of the order O = End(T) of E = Q(F), and prime-to-p twisted derived equivalences act through ideals of O, compatibly with Rizov's CM action on the canonical lift.
- **Statement 2 (Kneser neighbours).** Two such surfaces whose lattices T are ℓ-neighbours are related by a twisted equivalence with Brauer classes of order ℓ defined over F_q, and conversely.
- **The explicit example**: two non-isomorphic surfaces with the same zeta function, the ideal class separating them, and a step realised by an explicit twisted Fourier–Mukai partner. **The count**: the number of surfaces with a given Frobenius polynomial against the class number of the Hermitian genus. **The character question**: which characters of G are computable from the surface, as genus characters are for elliptic curves.

**Statement 1 is an assembly of published results, not an open problem.** I found no source that states it over F_q, but every step of the proof is in print, including in two sources listed at the end of §2. The proof has four steps: (i) E = Q(F) is a CM field over which T_Q has dimension one; (ii) the lattices T are ideal lattices of E, and their classes form a torsor under a norm-one class group; (iii) every such lattice is the transcendental lattice of a surface; (iv) the surfaces so obtained are related by twisted derived equivalences. Where each step is in print:

| Step of the proof | Where it already is | Status |
|---|---|---|
| (i) T_Q is one-dimensional over E = Q(F), E is CM | Taelman, proof of essential surjectivity ("X has complex multiplication") | [P] |
| (ii) lattices are ideals, classes form a torsor under a norm-one class group | Valloni 2021, §§6–9 (ideal lattices, "K3 class group" G_K3(E)); Bayer-Fluckiger 2024, Classification section (group C⁺(E), principal homogeneous space) | [P], maximal order only |
| (iii)–(iv) every lattice in the isogeny class is realised and reached by twisted derived equivalences | Bragg–Yang, existence theorem, over an algebraically closed field, p ≥ 5 | [P] |
| Compatibility with the CM main theorem | Taelman's descent step uses exactly this reciprocity; Valloni §8 gives the ideal-lattice form | [P] |

So **what remains is a worked example of an essentially known theorem, plus a short assembly note** (the assembly is [note 6](06-torsors-under-cm-class-groups.md)). What is not in print, as far as two searches and keyword scans of the seven main sources show:

1. the statement for **non-maximal orders** Z[F, q²/F] ⊆ O ⊊ O_E, including the F-stable sublattices that are not O-stable (the "volcano" direction); for odd ℓ ≠ p these are classified in [note 7](07-ring-changing-steps.md) (Propositions H to K), with the field of definition of each twisted equivalence, and computed for the six-line data;
2. the **Néron–Severi side over F_q**: the gluing, the action of F on NS, and twists (§4.3);
3. **F_q-rationality** of the individual twisted equivalences and Brauer classes (Bragg–Yang work over k = k̄);
4. any **explicit example or count**.

A search that finds nothing is weak evidence. Items 1–3 are small; someone who knows this literature may regard them as exercises.

---

## 2. Statements about the literature that need care

Six points about these sources are easy to get wrong. The table gives, for each, the correct statement and its source.

| Point | Correct statement and source |
|---|---|
| Rizov's CM theorem: which paper | Two different papers. *Complex multiplication for K3 surfaces* is arXiv:math/0508018 (2005), and Valloni still cites it as an arXiv preprint. The Crelle paper (648 (2010), 13–67) is *Kuga–Satake abelian varieties of K3 surfaces in mixed characteristic*. The main theorem is the preprint's corollary titled "Complex multiplication for K3 surfaces", which follows its theorem that the period morphism is defined over Q: for σ ∈ Aut(C/E) and an idèle s with art(s) = σ on E^ab, there is a Hodge isometry η: P²(X, Q(1)) → P²(X^σ, Q(1)) with η(r_X(s)·t) = σ(t) on adelic cohomology. [P] |
| Which CM theorem Taelman uses | Not Rizov's but his own integral version: *Complex multiplication and Shimura stacks*, arXiv:1707.01236, Theorem 2 (models of a CM K3 over extensions of E classified by their Galois representations). [P] for the citation, [M] for its content |
| Zarhin's results on the transcendental part | Two separate results. Field CM or totally real: Zarhin, *Hodge groups of K3 surfaces*, Crelle 341 (1983). Irreducibility: Yu–Yui Prop. 3.2 gives R(X, T) = Q(X, T)^r with r dividing the height h; for ordinary X (h = 1) it is irreducible, which they attribute to Zarhin, *Transcendental cycles*, Thm 1.1. [P] |
| Exclusion of the totally real case | Done by Taelman. Any endomorphism of V_X commuting with F respects the Hodge structure, End_Hdg(V_X) is a field, so V_X is a cyclic Q[F]-module, E is generated by F, dim_E V_X = 1, and then E is CM. He also records σ(F) = q²/F, so **F/q is the norm-one element**, and that the place v ∣ p with v(F) > 0 has E_v = Q_p. [P] |
| When Taelman's equivalence is unconditional | Three cases: an ample L with L² < p − 4; Pic contains a hyperbolic plane and p ≥ 5; geometric Picard rank ≥ 12 and p ≥ 5. For degree 2 the first case needs **p ≥ 7**. [P] |
| Valloni's neighbour correspondence and F_q | Title is *Neighbors and arithmetic of isogenous K3 surfaces*, arXiv:2206.02560. It is over C and number fields only. Neighbours correspond to **isotropic** lines in Br(X)[d], for d prime to disc T(X) (a later subsection relaxes this). X_ℓ has NS(X_ℓ) ≅ NS(X); the surface Y_ℓ with T(Y_ℓ) the index-d sublattice exists when ρ ≥ 3. The map ℓ ↦ X_ℓ is injective only when End_Hdg is Q or totally real, so **the CM case is exactly the one the paper sets aside**. No class groups, no finite fields. [P] |

Checked in the sources [P]: Taelman's Theorem B and (M1)–(M5); Bragg–Yang (Algebra Number Theory 17 (2023), 1069–1126), p ≥ 5, k algebraically closed, slope-< 1 condition, prime-to-p derived isogeny classes matched with prime-to-p Hecke orbits, strong approximation giving singleton double quotients; Huybrechts–Stellari 2005, T(X, α) = ker(T(X) → Q/Z), and for ρ ≥ 12 every orientation-preserving Hodge isometry of twisted Mukai lattices lifts to an equivalence.

Yang's paper is arXiv:1810.08546. Its main theorem needs p > 18d + 4 and goes through Kisin's work on Langlands–Rapoport; Bragg–Yang's second theorem strengthens it.

**Further sources:**

- **Valloni, *Complex multiplication and Brauer groups of K3 surfaces*, arXiv:1804.08763** (Adv. Math. 2021 [M]). This, not the neighbours paper, is the relevant Valloni paper. It defines ideal lattices (I, α) with form Tr(α x ȳ), the type of a CM K3, the main theorem in Shimura's form (X^τ has type (s/s̄)·(I, α) when art(s) = τ), and the **K3 class group G_K3(E) = U_E(Q) \ U_E(A_f) / Ũ** with its class field. It notes that s ↦ s/s̄ maps the finite idèles onto G_K3(E) by Hilbert 90. This is the group G of Statement 1 for the maximal order.
- **Bayer-Fluckiger, *K3 surfaces with maximal complex multiplication*, arXiv:2408.11643** (to appear in Documenta Math.). O_E-lattices with the same discriminant ideal and signature form a principal homogeneous space under C⁺(E) = {(I, α) : α I Ī = O_E, α positive at the required places}/∼. If no dyadic place of E₀ ramifies in E and the discriminant module has length < 22 − [E:Q], the complex K3 surfaces with CM by O_E and discriminant ideal 𝔇 form a **principal homogeneous space under C⁺(E)**. This is Statement 1 over C.
- **Srivastava, *On derived equivalences of K3 surfaces in positive characteristic*, Doc. Math. 24 (2019)**: Fourier–Mukai partners of an ordinary K3 correspond to those of the geometric generic fibre of its canonical lift, and automorphisms lift. **Kurama, arXiv:2410.14065**: Fourier–Mukai equivalences lift to canonical lifts. Together these are the untwisted part of the dictionary. [M], abstracts only.

---

## 3. What the sources give for free [D]

**3.1 The group is a Galois group, and one step always suffices.** Let O = End_{Z[F]}(T) and G_O = T¹(Q) \ T¹(A_f) / (T¹(A_f) ∩ Ô^×). The map s ↦ s/s̄ from the finite idèles of E onto T¹(A_f) is surjective (Hilbert 90 place by place, with units at unramified places), and its composite to G_O has open kernel containing E^×. So G_O is a quotient of a ray class group of E, that is G_O ≅ Gal(H_O/E) for a finite abelian extension H_O (Valloni's K3 class field when O = O_E). By Chebotarev, every element of G_O is the Frobenius of infinitely many degree-one primes 𝔩 of E, and the uniformiser at 𝔩 acts on lattices by T ↦ 𝔩 𝔩̄⁻¹ T. Hence:

- G_O is not only generated by the classes [𝔩 𝔩̄⁻¹]; **every** element is such a class, for infinitely many ℓ prime to p and to the discriminant.
- Any two vertices of the torsor are ℓ-neighbours for infinitely many ℓ, and each vertex is its own ℓ-neighbour for infinitely many ℓ. Connectivity of the Cayley graph for *small* ℓ is a question about small generators of G_O, as for CSIDH.
- By Valloni's Shimura-form theorem, the ℓ-neighbour step on the canonical lift is Galois conjugation by Frob_𝔩. For the maximal order this is the precise "compatibility with Rizov" of Statement 1.

**3.2 Which ℓ give steps.** Let L ⊂ T be F-stable of index ℓ. Then F acts on T/L by a scalar λ, so L ⊇ (ℓ, F − λ)·T. If ℓ does not divide [O_E : Z[F]] and λ is a simple root of the Frobenius polynomial mod ℓ, then 𝔩 = (ℓ, F − λ) is a degree-one prime and L = 𝔩T. The matching line in T/ℓT is an eigenline of F, and ⟨Fx, Fx⟩ = q²⟨x, x⟩ shows it is isotropic unless λ² ≡ q² mod ℓ, that is unless 𝔩 = 𝔩̄. So:

- the analogue of an Elkies prime is: the Frobenius polynomial on T has a simple root λ mod ℓ with λ ≢ q²/λ;
- the Brauer class is then isotropic and the neighbour is 𝔩 𝔩̄⁻¹ T, in agreement with Valloni's Kneser correspondence;
- if ℓ divides [O_E : O], there are F-stable sublattices that are not O-stable. These change the order. Statement 1, as formulated, does not treat them.

**3.3 Realisation.** Take x ∈ T¹(A_f) with trivial component at p, acting on T ⊗ A_f and trivially on NS. The lattice x·H²(X) is F-stable and keeps the p-part, so Bragg–Yang's existence theorem gives a derived isogeny X ⇝ X′ over F̄_q (p ≥ 5) to a surface X′ with the right ℓ-adic and crystalline lattices. Identifying X′ with the surface that Taelman's equivalence, where unconditional, attaches to the lattice x·H²(X) needs the isogeny to lift to the canonical lift; [note 6](06-torsors-under-cm-class-groups.md) §3 does this over F̄_q, using the proof of Bragg–Yang's Torelli theorem (Proposition E), and shows that a step at an Elkies prime is one twisted equivalence, from their construction for reflections (Proposition F). Those arguments are unchecked. By 3.1 every class has such a representative. So "every element of G is realised by a composition of twisted equivalences" holds over F̄_q, given Proposition E of note 6. Whether each twisted Fourier–Mukai kernel and each Brauer class is defined over F_q is not addressed by either paper.

---

## 4. Refinements of Statements 1 and 2 [D]

**4.1 Rank.** T has rank 22 − ρ̄ with ρ̄ the geometric Picard rank, which can exceed the Picard rank ρ over F_q. Over a finite field ρ̄ is even, so [E:Q] is even, as it must be.

**4.2 The torsor is fixed by the Hermitian genus.** The Frobenius polynomial and the NS-genus alone do not fix the torsor. The form on T_Q ≅ E is Tr(β x ȳ) with β ∈ E₀^×. Two surfaces with the same zeta function have β, β′ whose ratio need not be a norm from E, even when the quadratic spaces are isometric. The "same local data" of Statement 1 is the **E-Hermitian genus**: T(Y) ⊗ Z_ℓ ≅ T(X) ⊗ Z_ℓ as Z_ℓ[F]-lattices with form, for every ℓ including p. Then β′/β is a local norm everywhere, hence a global norm (Hasse, E/E₀ cyclic), and the set is a torsor under G_O. Bayer-Fluckiger's formulation is the natural larger one: fixing only the discriminant ideal gives a torsor under C⁺(E), which contains G_{O_E} and counts several genera at once. For the count, the number to compare with is |C⁺|, not |G|.

**4.3 The NS side includes twists.** Taelman's F acts on all of M. On NS(M, F) it has finite order and preserves the cone. Surfaces with the same ideal lattice T differ by the lattice NS within its genus, the glue isomorphism, and the conjugacy class of F on NS, that is by **twists over F_q**. Statement 1 fixes (NS, F|NS, K) up to isometry, so that the count does not count twists as new torsor points. When T embeds uniquely in the K3 lattice (Bayer-Fluckiger's length condition; automatic for ρ̄ ≥ 12, since the discriminant group then needs at most rank T ≤ 10 generators), the complex surface is determined by T, and only the twists remain.

**4.4 The explicit prime is the bad prime.** The natural place for the explicit example is ℓ = 2 on double planes, because that is where explicit twisted partners exist. But disc T is even for every degree-2 K3 surface (it is 2 when ρ̄ = 1), so ℓ = 2 violates the coprimality hypothesis of Kneser's correspondence, and dyadic ramification is exactly what Bayer-Fluckiger excludes. Van Geemen's several types of 2-torsion classes on a double plane are this failure made visible. So for ℓ prime to the discriminant, the range covered by Kneser's correspondence, Statement 2 cannot be tested on that example: it needs explicit ℓ = 3 steps, and the ℓ ∣ disc case starts from Valloni's "More neighbors" subsection.

**4.5 The character question.** By the Artin–Tate formula [M], |Br(X)|·|disc NS(X)| is determined by the zeta function, so it is constant on the torsor. Characters of G must come from finer data, such as the discriminant form or the Galois structure of Br(X̄)[2], not from discriminants. [H]

---

## 5. Where the work is done

1. **Statement 1 as a proposition**, for O = O_E, ρ̄ ≥ 12, p ≥ 5, fixed (NS, F|NS, K) and fixed Hermitian genus, with G_{O_E} = Gal(H/E) and the one-step property of 3.1: [note 6](06-torsors-under-cm-class-groups.md), Proposition A, with the steps in Propositions E to G.
2. **Non-maximal orders**: O in terms of (M, F), the F-stable non-O-stable sublattices of 3.2, and the resulting graph, the part with no source: [note 7](07-ring-changing-steps.md), and the count over all orders in [note 2](02-six-line-count.md) §4.
3. **The explicit example at ρ̄ ≥ 16**, not on generic double sextics: [note 2](02-six-line-count.md). A generic ordinary double sextic has T of rank up to 20, so E has degree up to 20 and its class group is expensive. Double planes branched along six lines have ρ̄ ≥ 16, so E has degree ≤ 6, the transcendental Frobenius polynomial follows from point counts over F_p, F_{p²}, F_{p³}, Taelman is unconditional, and Huybrechts–Stellari's ρ ≥ 12 lifting applies. The Kummer sublocus (lines tangent to a conic) reduces to ordinary abelian surfaces and is a check against known theory; ρ̄ = 20 reduces to binary quadratic forms.
4. **The count**: ordinary surfaces in such a family over small F_p with a given transcendental polynomial, against |C⁺(E)| times the number of twists: [note 2](02-six-line-count.md).
5. **The character question**, in the form of 4.5: [note 13](13-hardness-on-paper.md) §6, and notes [19](19-quadratic-characters-from-a-surface.md) and [24](24-the-form-mod-4-is-the-character.md).

---

## 6. Not checked

- Bragg–Yang and Yang were read in their introductions and by keyword only. The claim in 3.3 uses the existence theorem as stated in the introduction.
- Ito–Ito–Koshikawa, Huybrechts–Stellari 2007 (Căldăraru's conjecture), Taelman's *Shimura stacks* paper and Kisin's paper were not opened.
- Whether Br(X) → Br(X̄)^F is surjective for a K3 surface over F_q, which bears on item 3 of §1. It is, on ℓ-primary torsion for ℓ ≠ p not dividing disc NS, by [note 6](06-torsors-under-cm-class-groups.md), Proposition G: H²(X, μ_ℓ) → H²(X̄, μ_ℓ)^F is an isomorphism, and NS/ℓ splits off. This answers item 3 of §1 for steps at Elkies primes.
- A further Valloni title, *Fields of definition of K3 surfaces with complex multiplication*, was not located.
- The statements in §3 and §4 are my derivations. 3.1 is a standard class-field-theory argument; 4.2 and 4.3 should be re-derived before they are relied on.

---

## Sources

- Taelman, *Ordinary K3 surfaces over a finite field*: https://arxiv.org/abs/1711.09225 ; *Complex multiplication and Shimura stacks*: https://arxiv.org/abs/1707.01236
- Valloni, *Complex multiplication and Brauer groups of K3 surfaces*: https://arxiv.org/abs/1804.08763 ; *Neighbors and arithmetic of isogenous K3 surfaces*: https://arxiv.org/abs/2206.02560
- Bayer-Fluckiger, *K3 surfaces with maximal complex multiplication*: https://arxiv.org/abs/2408.11643
- Bragg, Yang, *Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*: https://arxiv.org/abs/2102.01193
- Yang, *Isogenies between K3 surfaces over F̄_p*: https://arxiv.org/abs/1810.08546
- Rizov, *Complex multiplication for K3 surfaces*: https://arxiv.org/abs/math/0508018
- Yu, Yui, *K3 surfaces of finite height over finite fields*: https://arxiv.org/abs/0709.1979
- Huybrechts, Stellari, *Equivalences of twisted K3 surfaces*: https://arxiv.org/abs/math/0409030
- Srivastava, *On derived equivalences of K3 surfaces in positive characteristic*: https://arxiv.org/abs/1809.08970
- Kurama, *Fourier–Mukai partners of abelian varieties and K3 surfaces in positive and mixed characteristics*: https://arxiv.org/abs/2410.14065
- Ito, Ito, Koshikawa, *CM liftings of K3 surfaces over finite fields and their applications to the Tate conjecture*: https://arxiv.org/abs/1809.09604 (bibliography used for the Rizov references)
