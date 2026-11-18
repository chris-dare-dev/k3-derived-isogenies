# 7. Rational steps that change the endomorphism ring

*Classifies all rational steps at an odd prime ℓ ≠ p, including those that change the endomorphism ring, and computes them on the six-line data. Builds on [note 6](06-torsors-under-cm-class-groups.md), whose Propositions F and G treat a step at a prime 𝔩 ≠ 𝔩̄ where End(T) is maximal. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_vertical.py`; **[H]** heuristic. Nothing here has been independently verified.*

---

## 1. Result

Notation is that of note 6: X is an ordinary K3 surface over F_q, p ≥ 5, ρ̄ ≥ 12, with transcendental lattice T, Frobenius F on T, V = q²/F, E = Q(F), E₀ its real subfield, and P the characteristic polynomial of F on T. Throughout, ℓ is an odd prime that does not divide q or disc T. In the six-line family disc T = 2⁶, so this is every odd ℓ ≠ p.

1. **All steps at ℓ (Proposition H).** The surfaces of the family joined to X ⊗ F̄_q by one twisted equivalence with Brauer classes of order ℓ are the Kneser ℓ-neighbours of T: one for each isotropic line of T/ℓT. Such a step is defined over F_q exactly when the line is an eigenline of F. If the eigenvalue is λ, the equivalence is defined over F_{q^m} exactly when (λ/q)^m = 1. A Brauer class that is not isotropic gives no new surface of the family.
2. **How many, and which change the ring (Proposition I).** The eigenlines with eigenvalue λ form a projective space of dimension μ − 1, where μ is the number of generators that T needs over Z[F, V] at the maximal ideal (ℓ, F − λ). So there is one eigenline when End(T) = Z[F, V] there, and at least ℓ + 1 otherwise. An eigenline that is not stable under End(T) leads to a lattice whose ring does not contain End(T).
3. **From a lattice that is maximal at ℓ (Proposition J).** The steps that keep the ring are the translations 𝔩⁻¹𝔩̄ at degree-one primes. All other rational steps go down, and they occur in a short list of configurations. In two of them the eigenvalue is ±q, so the twisted equivalence is defined over F_q: a pair 𝔩 ≠ 𝔩̄ of ramified primes, and a prime that is inert over E₀ and ramified over ℓ.
4. **A double root λ ≠ ±q gives the volcano of the elliptic-curve case** (Proposition K), and its twisted equivalences are never defined over F_q or F_{q²}.
5. **Data (§5).** Among the 245 polynomials of notes 2, 3 and 5, an odd prime ℓ ≠ p divides the index of Z[F, V] in O_E for 60, giving 61 pairs (polynomial, ℓ). Rational steps that change the ring occur in 14 pairs. In 9 of them, all at ℓ = 3, the Brauer class is defined over F_q, for the member of each twist pair whose Frobenius is ≡ q at the prime.

| what happens at ℓ | pairs | ℓ |
|---|---|---|
| volcano at a double root λ ≠ ±q | 5 | 7, 7, 7, 11, 13 |
| descent from a pair of ramified primes, eigenvalue ±q | 2 | 3 |
| descent from an inert ramified prime, eigenvalue ±q | 7 | 3 |
| a step between two genera with the same ring | 1 | 3 |
| translations only | 27 | 3, 5 |
| no rational step | 19 | 3, 5, 11 |

Every one of the 2,564 neighbours built is unimodular, F-stable, and isomorphic to exactly one lattice class already in the lists of notes 2, 3 and 5. Nothing is tested on equations.

§6 corrects two statements of note 6, which now states them in corrected form.

---

## 2. Proposition H: steps are isotropic lines

For an isotropic line F_ℓ·x of T/ℓT, choose x̃ ≡ x mod ℓT with ⟨x̃, x̃⟩ ≡ 0 mod ℓ², and put

> T′_x = Z·x̃/ℓ + { y ∈ T : ⟨x, y⟩ ≡ 0 mod ℓ }.

It depends only on the line. It is unimodular at ℓ, equals T at the other primes, and T ∩ T′_x has index ℓ in both. Every lattice with these properties arises once (Kneser) [M]. For T of rank 6 and discriminant 2⁶ the number of isotropic lines is (ℓ² + 1)(ℓ² + ℓ + 1) if ℓ ≡ 1 mod 4 and (ℓ³ + 1)(ℓ + 1) if ℓ ≡ 3 mod 4; the script confirms 112, 806 and 2,752 for ℓ = 3, 5, 7 [E].

Let α_x ∈ Br(X_k)[ℓ], k = F̄_q, be the class of x̃/ℓ: the functional t ↦ ⟨x̃, t⟩/ℓ mod Z_ℓ on T ⊗ Z_ℓ. Its kernel is T ∩ T′_x.

**Proposition H.** [D, from Taelman and Bragg–Yang as in note 6]

1. For every isotropic line there is a K3 surface Y_k over k with the same N, cone and crystalline part as X_k and with T(Y_k) = T′_x, and a Fourier–Mukai equivalence D^b(Y_k, α′) ≅ D^b(X_k, α_x) with α′ of order ℓ.
2. The following are equivalent: F·x ∈ F_ℓ·x; T′_x is F-stable; Y_k has a model Y over F_q whose triple is (N ⊕ T′_x glued as for X, F, 𝒦), and the isogeny X ⇝ Y is defined over F_q as a correspondence with rational coefficients.
3. If F·x = λx, Frobenius multiplies α_x by λ/q. So α_x, and the equivalence of (1), are defined over F_{q^m} exactly when (λ/q)^m = 1; over F_q exactly when λ = q.
4. If α ∈ Br(X_k) has order ℓ and D^b(X_k, α) ≅ D^b(Y′, β) with Y′ a surface of the family (so disc T(Y′) = disc T), then T(Y′) is isomorphic to T, or α = c·α_x for an isotropic x and T(Y′) is isomorphic to T′_x. In particular a class that is not isotropic leads to no surface of the family other than X_k.

**Proof.** (1) T′_x ⊗ Z_ℓ is unimodular and nothing else changes, so M′ = (N ⊕ T′_x with the gluing of M) is an even unimodular lattice with NS part N. Some power F^n fixes the line, so (M′, F^n, 𝒦) satisfies (M1)–(M5) and the positivity of Remark B, and Taelman's theorem gives a surface over F_{q^n}. For the equivalence, lift x̃ to an isotropic vector f of T ⊗ Z_ℓ (Hensel), choose an isotropic e with ⟨e, f⟩ = 1, and put b_ℓ = ℓe + f. Then b_ℓ is primitive, b_ℓ² = 2ℓ, and the reflection s_b sends e ↦ −f/ℓ, f ↦ −ℓe and fixes the complement, so s_b(T ⊗ Z_ℓ) = T′_x ⊗ Z_ℓ. This is the computation of Proposition F with T_𝔩 and T_𝔩̄ replaced by Z_ℓ·e and Z_ℓ·f; the rest of that proof (Bragg–Yang's statement for reflections, the B-field b/ℓ, and steps 2 to 5 of Proposition E to identify the moduli space with Y_k) uses only the image lattice and applies unchanged.

(2) F is a similitude of T ⊗ Z_ℓ with multiplier q², a unit. So F·T′_x = T′_{Fx}, and T′_x is F-stable exactly when F fixes the line. Then n = 1 in (1). The correspondence is defined over F_q by the averaging argument of Proposition G (1).

(3) Geometric Frobenius acts on H²(X_k, Z_ℓ(1)) by F/q, hence on the class of x̃/ℓ by λ/q. The descent is that of Proposition G (2) and (3), which uses nothing about the line except that Frobenius fixes the class.

(4) The equivalence induces an E-linear rational isometry φ: T(X) ⊗ Q → T(Y′) ⊗ Q that carries the twisted transcendental lattice K = ker α onto that of (Y′, β), [P for the ℓ-adic realisations, Bragg–Yang; over C this is Huybrechts–Stellari]. The twisted lattice of (Y′, β) has index ord(β) in T(Y′), and comparing discriminants gives ord(β) = ℓ. So φ⁻¹T(Y′) is a lattice that is unimodular at ℓ, contains K with index ℓ, and equals T away from ℓ. Write α = ⟨x, ·⟩/ℓ. If ⟨x, x⟩ ≢ 0 mod ℓ, then T ⊗ Z_ℓ = Z_ℓx ⊕ x^⊥ and K ⊗ Z_ℓ = ℓZ_ℓx ⊕ x^⊥ has cyclic discriminant group Z/ℓ², with a single isotropic subgroup; the only unimodular overlattice is T. If x is isotropic, the discriminant group of K is a hyperbolic plane over F_ℓ with two isotropic lines, giving T and T′_x. ∎

**Remarks.**

- *Inverse.* If T′ is the step from T along x with eigenvalue λ, then T is the step from T′ along a line with eigenvalue q²/λ.
- *Outside the family.* Every class α of order ℓ, isotropic or not, should also have an untwisted partner: a K3 surface Y_α with T(Y_α) = ker α, of discriminant ℓ²·disc T ([M] over C, by Huybrechts–Stellari; [H] in characteristic p). It is not in the family, since its Néron–Severi lattice has discriminant divisible by ℓ². Its ring is End(ker α), which contains Z + ℓ·End(T). I have not computed these.
- *What is inherited.* The unchecked step of Proposition E, and the descent of Proposition G (3).

---

## 3. Proposition I: counting the eigenlines

Since F is invertible at ℓ, Z_ℓ[F] contains V, and Z_ℓ[F] = Z_ℓ[t]/P is a product of local rings, one for each irreducible factor g of P mod ℓ. Accordingly T ⊗ Z_ℓ = ⊕ T_g, and T_g is orthogonal to T_h unless h = g*, the factor whose roots are q²/(roots of g). For λ ∈ F_ℓ let W_λ ⊂ T/ℓT be the eigenspace of F, and μ_λ its dimension.

**Proposition I.** [D]

1. μ_λ is the minimal number of generators of T as a Z[F, V]-module at the maximal ideal 𝔪_λ = (ℓ, F − λ). It is 1 exactly when End(T) and Z[F, V] agree at 𝔪_λ.
2. If λ ≠ ±q, every line of W_λ is isotropic. So there are (ℓ^μ − 1)/(ℓ − 1) rational steps with eigenvalue λ, and as many with eigenvalue q²/λ.
3. If λ = ±q, the radical of the form on W_λ is W_λ ∩ (F − λ)(T/ℓT), and μ_λ is even. For λ = q, W_q is dual to Br(X)[ℓ].
4. Let x span an eigenline. If T′_x is stable under End(T), the line is stable under End(T). If λ ≠ ±q the converse holds, and then End(T′_x) ⊇ End(T) at ℓ.

**Proof.** (1) W_λ is the kernel of F − λ on T/ℓT, which has the dimension of the cokernel T/𝔪_λT; Nakayama. Z_ℓ[t]/P is generated by one element, hence Gorenstein, and over a Gorenstein order a lattice is locally free exactly when its ring of multipliers is the order (Bass) [M].

(2) For x, y ∈ W_λ, q²⟨x, y⟩ = ⟨Fx, Fy⟩ = λ²⟨x, y⟩.

(3) For x ∈ W_λ, Vx = (q²/λ)x = λx, so ⟨(F − λ)y, x⟩ = ⟨y, (V − λ)x⟩ = 0. So the image of F − λ lies in W_λ^⊥, and they have the same dimension. For the parity take λ = q and g = F/q, an isometry of T ⊗ Z_ℓ. On coker(g − 1) put b(x, y) = ⟨(g − 1)⁻¹x, y⟩ mod Z_ℓ. A direct computation gives b(x, y) + b(y, x) = −⟨x, y⟩ and b(x, x) = −⟨x, x⟩/2, both in Z_ℓ, and b is non-degenerate because T is unimodular at ℓ. So coker(F − q) carries a non-degenerate alternating form and is of the shape A × A. This is the lattice form of the fact that |Br(X)| is a square. The last sentence is the description of Br(X)[ℓ^∞] in note 6.

(4) Put O = End(T) ⊗ Z_ℓ; it is stable under the involution because T is unimodular at ℓ. If T′_x is O-stable, so is T ∩ T′_x = x^⊥, and for a ∈ O and y ∈ x^⊥, ⟨āx, y⟩ = ⟨x, ay⟩ ≡ 0, so āx lies in (x^⊥)^⊥ = F_ℓx. Conversely let λ ≠ ±q and g = t − λ. Then x̃ can be taken in T_g, which is isotropic, and T′_x = (T_g + Z_ℓx̃/ℓ) ⊕ (T_{g*} ∩ x^⊥) ⊕ (the rest). If the line is O-stable, both summands are. ∎

In every computed case with λ = ±q an O-stable eigenline also leads to a ring containing End(T) [E]; no proof is given here. It is proved when ℓ² is the exact power of ℓ in P(±q): in [note 8](08-below-the-top.md) when a lattice maximal at ℓ exists (Theorem M (3), (5), (6)), and in [note 9](09-ramified-over-e0.md) when none does (Theorem N (6)). So it holds whenever ℓ² is the exact power of ℓ in P(±q) and ±q has multiplicity ≥ 4 mod ℓ. Otherwise it is open.

---

## 4. Lattices maximal at ℓ, and volcanoes

**Proposition J.** [D] Suppose End(T) is maximal at ℓ. For a prime 𝔩 of E above ℓ, with ramification index e and residue degree f, let S_𝔩 = 𝔩^{e−1}T/ℓT, the socle of the 𝔩-part of T/ℓT; it has dimension f.

1. No prime above ℓ is ramified over E₀. So each prime is either different from its conjugate, or inert over E₀ with f even.
2. The eigenlines stable under O_E are the S_𝔩 with f = 1. Each is isotropic, and its step is the translation T′ = 𝔩⁻¹𝔩̄T. This extends Proposition F to ramified degree-one primes and to eigenvalue ±q.
3. Every other isotropic eigenline leads to a lattice that is not maximal at ℓ.
4. S_𝔩 is orthogonal to S_𝔩′ unless 𝔩′ = 𝔩̄ and e = 1. Hence, when F ≡ λ mod 𝔩 with λ ∈ F_ℓ:

| configuration | eigenvalue | isotropic eigenlines in the socles | of which translations | going down |
|---|---|---|---|---|
| 𝔩₁, 𝔩₂ of degree one, same λ, 𝔩₂ ≠ 𝔩̄₁ | any | ℓ + 1 | 2 | ℓ − 1 |
| 𝔩 ≠ 𝔩̄ of degree two | any | ℓ + 1 | 0 | ℓ + 1 |
| 𝔩 ≠ 𝔩̄ of degree one, e = 1, F ≡ λ mod both | ±q | 2 | 2 | 0 |
| 𝔩 ≠ 𝔩̄ of degree one, e ≥ 2, F ≡ λ mod both | ±q | ℓ + 1 | 2 | ℓ − 1 |
| 𝔩 = 𝔩̄, f = 2, e = 1 | ±q | 0 | 0 | 0 |
| 𝔩 = 𝔩̄, f = 2, e ≥ 2 | ±q | ℓ + 1 | 0 | ℓ + 1 |
| 𝔩 = 𝔩̄, f = 4, e = 1 | ±q | ℓ² + 1 | 0 | ℓ² + 1 |

In the rows with 𝔩 = 𝔩̄ the eigenvalue is forced to be ±q.

**Proof.** (1) A unimodular O_E-lattice at ℓ is 𝔞 with form Tr(α x ȳ), α ∈ E₀, and α𝔞𝔞̄ equal to the inverse different. At a prime 𝔓 ramified over E₀ the valuations of α and of 𝔞𝔞̄ are even, and that of the different of E is odd, because ℓ is odd and the different of E over E₀ has valuation 1 there. (2) T/ℓT is the sum of the O_E/𝔩^e, and a line stable under O_E is a simple submodule with residue field F_ℓ. By (1) such a prime has 𝔩 ≠ 𝔩̄, so the T_𝔩 are isotropic. With x in S_𝔩, x̃/ℓ generates 𝔩⁻¹T_𝔩 over T_𝔩 and x^⊥ meets T_𝔩̄ in 𝔩̄T_𝔩̄. (3) is Proposition I (4), since a lattice maximal at ℓ is stable under O_E. (4) T_𝔩 is orthogonal to T_𝔩′ unless 𝔩′ = 𝔩̄, because ⟨ax, y⟩ = ⟨x, āy⟩. For a ∈ 𝔩^{e−1} and b ∈ 𝔩̄^{e−1}, ⟨ax, by⟩ = ⟨b̄ax, y⟩ with b̄a ∈ 𝔩^{2e−2}, which lies in ℓO_E at 𝔩 exactly when e ≥ 2; for e = 1 the socle is the whole summand and the pairing is perfect. The table follows: a totally isotropic plane has ℓ + 1 lines; a hyperbolic plane has 2; for 𝔩 = 𝔩̄ and e = 1 the form on S_𝔩 = F_{ℓ^f} is the trace to F_ℓ of a Hermitian norm form, which is anisotropic for f = 2 and an elliptic quadric with ℓ² + 1 points for f = 4. If 𝔩 = 𝔩̄ the involution acts on the residue field as the Frobenius of F_{ℓ^f} over F_{ℓ^{f/2}} and fixes λ, so λ = q²/λ. ∎

The table lists the cases where W_λ is the sum of the socles. If F − λ vanishes to higher order at a ramified prime, W_λ is larger; that does not occur in the data.

**Proposition K (double root).** [D, M] Let λ ≠ ±q be a root of P mod ℓ of multiplicity exactly 2, h the factor of P over Z_ℓ lifting (t − λ)², A = Q_ℓ[t]/h, and ℓ^c the index of Z_ℓ[t]/h in the maximal order of A. Then the part T_h of T is homothetic to exactly one order Z_ℓ + ℓ^i O_A, 0 ≤ i ≤ c, its *level*, and the steps with eigenvalue λ are those of an ℓ-volcano:

- level c ≥ 1 (the floor): one step, to level c − 1;
- level 0 < i < c: ℓ + 1 steps, one to level i − 1 and ℓ to level i + 1;
- level 0 < c: ℓ + 1 steps, of which 2, 1 or 0 are translations as A is split, ramified or an unramified field, and the others go to level 1;
- c = 0: A is a ramified field, and the single step is a translation.

The index of End(T) in O_E at ℓ is ℓ^{2i}, the parts at h and h* contributing equally. The twisted equivalences are defined over F_{q^d} with d the order of λ/q, and d ≥ 3.

**Proof.** Steps with eigenvalue λ are the overlattices of T_h of index ℓ stable under F (the dual part T_{h*} follows). T_h is a rank-two lattice with multiplication by an order of the quadratic algebra A, and the classification of its neighbours is the standard one for lattices in a quadratic algebra (Kohel's thesis) [M]. d ≥ 3 because λ² ≠ q². ∎

So the equivalence realising a volcano step is never rational over F_q or F_{q²}. The steps with a rational Brauer class are the last rows of the table in Proposition J.

---

## 5. Computation [E]

`sixlines_vertical.py` takes each of the 245 polynomials of notes 2, 3 and 5 with its list of lattice classes, and each odd prime ℓ ≠ p dividing the index of Z[F, V] in O_E. For each class it lists the eigenlines of F on L/ℓL, builds the neighbour of each isotropic one, finds the class of the neighbour in the list, compares the rings at ℓ, and tests whether the neighbour is the translate (a/ā)L by a uniformiser a at one prime above ℓ.

**Totals.** 60 polynomials, 61 pairs (polynomial, ℓ), 1,545 pairs (lattice class, ℓ), 5,398 eigenlines: 2,834 anisotropic and 2,564 isotropic. Of the 2,564 steps, 1,518 are translations, 610 go down, 412 go up, and 24 keep the ring without being translations.

**Checks.** Every neighbour is unimodular and F-stable. Every neighbour is isomorphic to exactly one class of the list, so the lists of notes 2, 3 and 5 are closed under rational steps. The numbers of steps satisfy n(i → j, λ)·|Aut j| = n(j → i, q²/λ)·|Aut i| for every pair of classes. coker(F − q) and coker(F + q) are of the shape A × A at every vertex: trivial, or (Z/ℓ)² for exactly one of the two signs (48 pairs), in agreement with §6.6 of note 6.

**The 14 pairs with ring-changing steps.** Levels are given by the index of End(L) in O_E at ℓ; the entries are numbers of lattice classes.

| q | (e₁, e₂, e₃) | lines | ℓ | primes above ℓ involved | eigenvalue | classes at index 1, ℓ², ℓ⁴ | from the top | d |
|---|---|---|---|---|---|---|---|---|
| 29 | (14, 87, −23548) | rational | 13 | two of degree one | 5, 7 | 1, 6 | 2 translations, 12 down | 12 |
| 37 | (38, 1591, 27380) | rational | 7 | 𝔩 ≠ 𝔩̄ of degree two | 1, 4 | 2, 8 | 8 down | 3 or 6 |
| 13 | (6, −65, −2028) | permuted | 7 | 𝔩 ≠ 𝔩̄ of degree two | 3, 5 | 6, 24 | 8 down | 3 or 6 |
| 19 | (4, 19, 5776) | permuted | 7 | 𝔩 ≠ 𝔩̄ of degree two | 1, 4 | 4, 32 | 8 down | 3 or 6 |
| 17 | (10, −153, −5780) | permuted | 11 | 𝔩 ≠ 𝔩̄ of degree two | 2, 7 | 9, 54 | 12 down | 5 or 10 |
| 23 | (16, 667, 4232) | permuted | 3 | 𝔩 ≠ 𝔩̄, e = 2, f = 1 | −q | 10, 20 | 2 translations, 2 down | 1 |
| 23 | (22, 391, 10580) | permuted | 3 | 𝔩 ≠ 𝔩̄, e = 2, f = 1 | −q | 24, 48 | 2 translations, 2 down | 1 |
| 17 | (22, 391, 5780) | permuted | 3 | 𝔩 = 𝔩̄, e = 2, f = 2 | −q | 3, 6, 18 | 4 down | 1 |
| 23 | (28, 391, 4232) | permuted | 3 | 𝔩 = 𝔩̄, e = 2, f = 2 | −q | 6, 24 | 4 down | 1 |
| 37 | (66, 4551, 169756) | rational | 3 | 𝔩 = 𝔩̄, e = 3, f = 2 | −q | 1, 2 | 4 down | 1 |
| 49 | (78, 5439, 240100) | rational | 3 | 𝔩 = 𝔩̄, e = 3, f = 2 | −q | 2, 8 | 4 down | 1 |
| 13 | (18, 351, 4732) | permuted | 3 | 𝔩 = 𝔩̄, e = 3, f = 2 | −q | 3, 6 | 4 down | 1 |
| 19 | (12, 627, 5776) | permuted | 3 | 𝔩 = 𝔩̄, e = 3, f = 2 | −q | 4, 16 | 4 down | 1 |
| 19 | (24, 399, 5776) | permuted | 3 | 𝔩 = 𝔩̄, e = 3, f = 2 | −q | 4, 16 | 4 down | 1 |

Here d is the degree of the extension of F_q over which the twisted equivalence is defined. Polynomials are normalised with e₁ > 0, and eigenvalue −q means d = 1 for the member of each twist pair with Frobenius −F, the one whose polynomial is P(−t) and whose Brauer group has order 9 by the Artin–Tate formula; for the other member d = 2.

What the rows show:

- **Rows 1 to 5** are Proposition K with c = 1, top split in row 1 and an unramified field in rows 2 to 5. On the floor each class has one eigenline per eigenvalue, going up.
- **Rows 6 and 7** are the fourth row of the table in Proposition J. On the floor each class has four eigenlines, one isotropic and stable, going up, and three anisotropic.
- **Rows 8 to 14** are the sixth row of that table. In row 8 the descent continues: at index ℓ² each class has one eigenline going up and three going down to index ℓ⁴, and at index ℓ⁴ one going up and three anisotropic. For each type of form the nine classes form a tree with 1, 2 and 6 vertices at the three levels.
- In rows 6 to 14 the Brauer group at ℓ is (Z/3)² at every level. A surface at the top has 8 rational Brauer classes of order 3, all isotropic; one at the bottom has 2 isotropic ones and 6 that lead nowhere in the family.
- The lists stop above the level of Z[F, V] in rows 6 to 14: the index of Z[F, V] at 3 is 3⁵ in row 8 and 3³ in the others, and the lowest lattices have index 3⁴ and 3².

**The other 47 pairs.**

- *Translations only* (27 pairs). In 21 the eigenvalue is −q, W is a hyperbolic plane, and the two isotropic lines are the two translations of Proposition G; these are the 21 polynomials of note 6 §6.6. In 3, at ℓ = 5, the eigenvalues are simple roots and ℓ divides the index only through a factor of degree 2. In 3, at ℓ = 3, both primes above 3 are ramified over E₀, no lattice is maximal at 3, and each class has one isotropic eigenline; its step is translation by a/ā for a uniformiser a at either prime, a fixed-point-free involution of the classes with a rational Brauer class.
- *A step between genera with the same ring* (1 pair): q = 23, (34, 391, −2116), ℓ = 3. The prime above 3 has e = 2, f = 2 and is ramified over E₀. All 24 classes have ring of index 3. Each has one isotropic eigenline; its step pairs the classes into 12 pairs, and it is not a translation, since translation by a/ā fixes every lattice. Reason [D]: over the maximal order there is a lattice M with dual 𝔓M, and M/𝔓M = F_9 carries a form over F_3 that is a hyperbolic plane; its two isotropic lines give the two unimodular lattices, each is the neighbour of the other, and every unit of norm one is ≡ ±1 mod 𝔓 and fixes both. So the two lattices lie in different genera in the sense of Proposition A, with the same endomorphism ring, and they are joined by a twisted equivalence defined over F_q.
- *No rational step* (19 pairs). In 14 the eigenvalue is −q, the prime is inert over E₀ and unramified, and all eigenlines are anisotropic (fifth row of the table): the surface has ℓ² − 1 rational Brauer classes of order ℓ and none of them leads to another surface of the family. In 5 there is no eigenvalue in F_ℓ.

---

## 6. Corrections to note 6, and what this does not give

**Corrections.**

1. Note 6, remark after Proposition G: invariant Brauer classes with non-stable kernel do not give rational twisted partners with a smaller endomorphism ring when 𝔩 = 𝔩̄ has residue field F_{ℓ²}. For 𝔩 = 𝔩̄ unramified with f = 2 there is no such partner inside the family: those classes are anisotropic and have no partner in the family (Proposition H (4)). The partner outside the family, with T = ker α, may exist; see the remark in §2. The descents inside the family come from ramified primes.
2. Note 6 §6.6: the 27 pairs with odd torsion in Br(X) and no Elkies pair do not all have a prime above ℓ that is its own conjugate and no step. Of these 27, 14 have no rational step; 7 have descents from an inert ramified prime; 2 have a pair 𝔩 ≠ 𝔩̄ of ramified primes, with two translations and two descents; 3 have a translation at primes ramified over E₀; 1 has the step between genera described above.

**Not established.**

- Everything in §7 of note 6 that Propositions E, F and G inherit: the identification step of Proposition E and the descent of Proposition G (3) are used again here.
- Proposition H (4) relies on the equivalence inducing an isometry of ℓ-adic twisted transcendental lattices; I took this from Bragg–Yang's realisations without rereading the statement for this purpose.
- A general description of the graph below the top level when the eigenvalue is ±q, beyond the case where ℓ² is the exact power of ℓ dividing P(±q). The data show one eigenline going up from every lattice that is not maximal, and descents continuing in one case; Proposition K proves this only for double roots λ ≠ ±q. Under the hypothesis that ℓ² is the exact power of ℓ dividing P(±q), which holds throughout the data, [note 8](08-below-the-top.md) shows that the graph is an ℓ-volcano of depth n, with 2n + 1 the exponent of ℓ in the index of Z_ℓ[F], and End = Z_ℓ[F] + 𝔏^{2i}O_E at level i, checked on all 221 lattice classes concerned; and [note 9](09-ramified-over-e0.md) does the remaining configurations, with primes ramified over E₀ and no maximal lattice, under the same hypothesis: two trees of depth n − 1 joined at their roots, with ℓ^{2n} the index of Z_ℓ[F] and End = Z_ℓ[F] + 𝔏^{2j+1}O_E at type j. In the data n = 1: two lattices exchanged by one step.
- Steps at eigenvalues that are not in F_ℓ. A factor of degree 2 of P mod ℓ with multiplicity ≥ 2 changes the ring through sublattices of index ℓ², which are compositions of two steps over F_{q²} and are not rational ℓ-steps. This is why, in nine pairs of the data, the lattice classes have two different rings at ℓ and no rational ℓ-step joins them.
- ℓ = 2 and ℓ = p.
- The untwisted partners outside the family.
- Nothing is tested on surfaces. The predictions are concrete: for example, the 17 six-line sets of type A over F₁₇ with polynomial ±P, P given by (22, 391, 5780), each taken in the twist with |Br| = 9, lie on 9 lattice classes at three levels (1, 2 and 6 classes), joined by twisted equivalences with Brauer classes of order 3 defined over F₁₇.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_vertical.py --out sixlines_vertical.jsonl > sixlines_vertical.log
```

It reads the result files of notes 2, 3 and 5 and takes about three minutes. `--only q,e1,e2,e3` runs one polynomial.

## Sources

As in [note 6](06-torsors-under-cm-class-groups.md). In addition:

- Kneser, *Klassenzahlen definiter quadratischer Formen*, Arch. Math. 8 (1957) (neighbours) [M]
- Kohel, *Endomorphism rings of elliptic curves over finite fields*, thesis, Berkeley 1996 (the volcano) [M]
- Bass, *On the ubiquity of Gorenstein rings*, Math. Z. 82 (1963) (lattices over Gorenstein orders) [M]
