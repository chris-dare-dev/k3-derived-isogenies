# 22. A step at 2 below Picard number 16: two neighbours of elliptic K3 surfaces

*Tests the lead of [note 17](17-below-picard-16.md) §5.3: a double cover of P¹ × P¹ branched along a curve of bidegree (4, 4) has two genus-one fibrations, and their Jacobians should be two elliptic K3 surfaces joined by a step at 2 of the class-group action. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/ellk3_two_neighbours.py`; **[H]** heuristic. Not re-checked. Class groups are PARI's without certification, so they are correct under GRH.*

---

## 1. Result

The pair was built over F_5 and F_7 for 21,240 double covers in three families, with geometric Picard number 12, 10 and 2. Everything that the lattice side predicts and that can be read off the two surfaces holds, without exception. The one statement that would identify *which* element of the group the step is could not be exercised.

1. **The construction is right [E].** For every double cover Y and both fibrations, the genus-one curve of Y over each point of F_p and the fibre of the Jacobian have the same number of points (159,680 fibres). The two Jacobians X₁, X₂ have the same Frobenius polynomial in every sample. These two facts are forced by Lang's theorem and say only that the equations are right.
2. **The arithmetic at 2 is as predicted [E].** For a pair, with P irreducible and ordinary (3,613 cases):
   - t − 1 divides P mod 2 at least twice, always. For random surfaces of the same fibre types that do not come from a double cover it does not divide P at all in 77% of cases.
   - When it divides exactly twice (1,349 cases) the two primes above 2 are a conjugate pair of degree 1, always. In the controls 17 of 49 such cases are inert instead.
   - The order of the Brauer group given by the Artin–Tate formula is 4 times a square, always. In the controls it is even in 35% of cases.
3. **The second Jacobian is never the first [E].** Under a condition that forces the endomorphism ring to be maximal (condition C, 972 pairs) the prediction is: X₁ ≅ X₂ over F_p exactly when the class c of 𝔩/𝔩̄ is trivial, 𝔩 a prime of degree 1 above 2. In all 972 the class is not trivial and the surfaces are not isomorphic. If X₂ were a random member of the isogeny class, about 23 isomorphic pairs would be expected.
4. **A surface has at most two partners [E].** Under condition C a surface has exactly two classes of order 2 that give such a double cover, so it can occur with at most two different partners. Of the 1,594 surfaces of the family with a cusp none has three; 23 occur in three or four pairs. With several pairs of primes above 2, where more partners are allowed, 5 surfaces have three.
5. **Not tested: that the step is c and not another element.** The class c was never trivial, so "isomorphic exactly when trivial" was confirmed only in the direction "not trivial, not isomorphic". This is not bad luck: a trivial class needs an integer of E all of whose absolute values are 2, which only finitely many fields contain (§5). No cycle of the step was closed either; the sampling is too thin.
6. **Group sizes [E].** For the 972 fields, |G| has median 0.08·p⁴ at Picard number 12 (both primes) and 0.065·p⁵ at Picard number 10. This is the first measurement of the heuristic q^{(20−ρ̄)/2} of note 13 away from ρ̄ = 16, at two exponents and tiny p.
7. **Picard number 2 gave almost nothing [E].** Of 40 generic double covers over F_5, 18 gave a polynomial of degree 20, 15 of them reducible. The 3 irreducible ones satisfy item 2 (t − 1 divides P mod 2 four or six times) and none satisfies condition C. Each sample costs about 15 minutes.

So below Picard number 16 there is now one family in which a surface is named, its zeta function computed, and the pair (X, step of X) written down with equations, at ℓ = 2. It is half of a step in the sense of note 17: from a double cover to its two Jacobians. Going from a given surface to its partner still needs the double cover, which is a 2-descent that was not done. Nothing here changes the verdict of note 17.

| Family | ρ̄ | p | Pairs | With a polynomial | Irreducible, ordinary | Condition C | Class trivial | X₁ ≅ X₂ |
|---|---|---|---|---|---|---|---|---|
| E8 point and a cusp | 12 | 5 | 16,000 | 6,867 | 2,621 | 722 | 0 | 0 (2 outside the scope, §4.3) |
| E8 point and a cusp | 12 | 7 | 4,000 | 1,801 | 777 | 201 | 0 | 0 |
| E8 point | 10 | 5 | 1,200 | 533 | 212 | 49 | 0 | 0 |
| generic | 2 | 5 | 40 | 18 | 3 | 0 | — | 0 |

---

## 2. The construction

### 2.1 The pair

Let f(s, u) be of degree 4 in each variable over F_p, p ≥ 5, and Y: w² = f(s, u), resolved. Y → P¹_u has fibres w² = (quartic in s), and Y → P¹_s likewise. For a binary quartic a s⁴ + b s³ + c s² + d s + e put

> I = 12ae − 3bd + c²,  J = 72ace + 9bcd − 27ad² − 27eb² − 2c³.

The Jacobian of w² = quartic is y² = x³ − 27 I x − 27 J [M]. With the coefficients polynomials in u of degree 4, I and J have degrees 8 and 12, so

> X₁: y² = x³ − 27 I₁(u) x − 27 J₁(u),  X₂: y² = x³ − 27 I₂(s) x − 27 J₂(s)

are Weierstrass equations of elliptic K3 surfaces: the relative Jacobians of the two fibrations.

Over C [P: van Geemen, *Some remarks on Brauer groups of K3 surfaces*, 1.2, 4.1, 4.7, 4.8, 5.5]: for a general elliptic K3 surface X with a section and α in Br(X) of order 2, the genus-one fibration Y attached to α has T_Y = ker α; when NS(Y) is his lattice Λ_{2,0}, Y is a double cover of P¹ × P¹ branched along a (4, 4) curve, it has two genus-one fibrations, and their Jacobians are not isomorphic in general, the discriminant group (Z/2)² of T_Y having two maximal isotropic subgroups. So T(X₁) and T(X₂) are the two even unimodular overlattices of T(Y): they are 2-neighbours.

### 2.2 Three families

A generic f gives transcendental rank 20 and a Frobenius polynomial of degree 20, which costs point counts over F_{p^10}. To have small groups and many samples I imposed singular points on the branch curve. A rational double point of the branch curve lies in a fibre of both fibrations, so it adds the same root lattice to NS(Y), NS(X₁) and NS(X₂).

| Family | Branch curve | Fibres of X_i | NS(X_i) | rank T | Moduli |
|---|---|---|---|---|---|
| generic | smooth | 24 I₁ | U | 20 | 18 |
| e8 | a point of type E8 at (0, 0), tangent s = u | II* at 0, 14 I₁ | U ⊕ E8 | 12 | 10 |
| e8a2 | in addition a cusp at (∞, ∞), tangent 1/s = λ/u, λ random | II* at 0, I₃ at ∞, 11 I₁ | U ⊕ E8 ⊕ A2 | 10 | 8 |

- **How the singular points are imposed [D].** With t = s − u, a point of type E8 at the origin with tangent t = 0 is the vanishing of the coefficients of t^a u^b for a/3 + b/5 < 1: eleven linear conditions on the 25 coefficients of f. The tangent cannot be a line of a ruling: such a line would meet the curve five times and lie in it. A cusp at (∞, ∞) with tangent t′ = 1/s − λ/u = 0 is five more linear conditions, a/2 + b/3 < 1. The script draws a random vector of the solution space.
- **The fibre types are checked, not assumed [E].** A sample is kept only if the discriminants of both X₁ and X₂ have the expected shape: degree 24 and squarefree (generic); valuations of (A, B, Δ) at 0 equal to (≥ 4, 5, 10) and Δ/u¹⁰ squarefree of degree 14 (e8) or 11 with the fibre at ∞ multiplicative (e8a2). About one draw in eight passes for e8a2 at p = 5.
- **Why unimodular or nearly so [D].** NS(X_i) = U and U ⊕ E8 are unimodular, so T is unimodular and there is no gluing: one surface for each lattice (T, F). For U ⊕ E8 ⊕ A2 the discriminant group is Z/3 and both of its isometries come from automorphisms of the Néron–Severi lattice that preserve the ample cone (the identity and the inversion of the fibration), so again one surface for each lattice. And 2 does not divide the discriminant of T in any of the three.
- **Taelman's theorem is unconditional here**, by case (ii) of his Theorem 4.6: the Picard lattice contains a hyperbolic plane and p ≥ 5 [P, note 17 §5.1].

### 2.3 The Frobenius polynomial

For an elliptic surface with these fibres the L-function L(T) of the generic fibre has degree N = 20, 12, 10, with reciprocal roots of absolute value p: the eigenvalues of Frobenius on H² modulo the classes of the fibre types [M]. The script computes S_d, the sum over the points v of P¹(F_{p^d}) of the trace of the fibre (the trace of Frobenius for a smooth fibre, ±1 for a node, 0 for a cusp), for d ≤ N/2, from one point in each Frobenius orbit; then L = exp(Σ S_d T^d/d), and the upper half of the coefficients from the functional equation with sign +. P(t) = t^N L(1/t).

- When the middle coefficient of L vanishes the sign is not determined and the sample is set aside: 12,021 of 21,240. These include all the surfaces with sign −, whose polynomial has the roots p or −p and whose Picard number is larger.
- Checks on each polynomial [E]: the truncated series has integer coefficients; all roots have absolute value p (9,219 of 9,219 pairs, 975 of 975 controls).
- A second implementation, summing over all points of P¹(F_{p^d}), is kept in the script; the two agree on six surfaces, two of each family, for d ≤ 6. The second Jacobian is computed with it, for d ≤ 5 (e8a2), 6 (e8), 7 (generic).

---

## 3. The lattice prediction

Notation of note 6 §2 and note 13 §3.1: E = Q[F] of degree N, V = p²/F its conjugate, E₀ the real subfield, G = Cl(E)/⟨Cl(E₀), primes ramified over E₀⟩ (note 13, Lemma 9, whose proof does not use the rank).

### 3.1 When the endomorphism ring is maximal: condition C

**Lemma 1.** [D] Let P be irreducible and ordinary and (T, F) an even unimodular lattice with characteristic polynomial P coming from a surface. Then End(T) = O_E if:

1. *at p:* the polynomial c₁t^{N−2} + (c₂/p)t^{N−3} + … + c_{N−1}/p^{N−2} is squarefree mod p, where P = t^N + c₁t^{N−1} + …;
2. *at odd ℓ ≠ p:* ℓ does not divide [O_E : Z[F, V]];
3. *at 2:* P mod 2 = (t − 1)²·(squarefree), and above 2 there are two primes 𝔩 ≠ 𝔩̄ of degree 1.

*Proof.* (1) By (M4) of Taelman's theorem T ⊗ Z_p = T⁰ ⊕ T¹ ⊕ T² with ranks 1, N − 2, 1 and F = p·θ on T¹ with θ invertible. The polynomial displayed is, up to the unit c₁, the reduction of the characteristic polynomial of θ; if it is squarefree, Z_p[θ] is étale, hence maximal. (2) is clear. (3) By Hensel's lemma T ⊗ Z_2 splits according to the coprime factors of P mod 2, orthogonally except for conjugate factors. On the squarefree part Z_2[F] is maximal. The part belonging to (t − 1)² is a lattice Λ of rank 2 in E_𝔩 × E_𝔩̄ = Q_2 × Q_2 with form β(ab′ + a′b). Write Λ = ⟨(2^i, c), (0, 2^j)⟩. Its Gram matrix is β·[[2^{i+1}c, 2^{i+j}], [2^{i+j}, 0]]. Unimodular means v(β) + i + j = 0, and even then means v(c) ≥ j, so Λ = ⟨(2^i, 0), (0, 2^j)⟩: a module over Z_2 × Z_2. ∎

**Condition C** is (1), (2), (3) together with: P irreducible, ordinary, and the maximal order certified by PARI. At (3) the evenness of T does the work: without it the two primes could be glued.

### 3.2 The step

Assume condition C, so T(X₁) = 𝔞, a fractional ideal with the form Tr(β x ȳ).

- [D] F ≡ 1 mod 2 on the part of 𝔞/2𝔞 belonging to (t − 1)², of dimension 2, and F − 1 is invertible on the rest. So 𝔞 has exactly three F-stable sublattices of index 2: 𝔩𝔞, 𝔩̄𝔞, and the one glued along the diagonal. The first two have discriminant form the hyperbolic plane U(2), with two isotropic lines; the third is isometric at 2 to ⟨2⟩ ⊕ ⟨−2⟩ with one. In van Geemen's terms these are two surfaces of type Λ_{2,0} and one of type Λ_{2,1}. **So X₁ has exactly two rational Brauer classes of order 2 that give a double cover of P¹ × P¹.**
- [D] The even unimodular lattices between 𝔩𝔞 and its dual 𝔩̄⁻¹𝔞 are 𝔞 and 𝔩𝔩̄⁻¹𝔞. So if T(Y) = 𝔩𝔞 then **T(X₂) = 𝔩𝔩̄⁻¹·T(X₁)**: X₂ = c·X₁ with c the class of 𝔩/𝔩̄ in G, or its inverse if T(Y) = 𝔩̄𝔞.
- [D] (T(X₁), F) ≅ (T(X₂), F) if and only if 𝔩/𝔩̄ = (x) with x x̄ = 1. With §2.2 and the full faithfulness of Taelman's functor: **X₁ ≅ X₂ over F_p if and only if 𝔩/𝔩̄ is generated by an element of relative norm 1.**
- [D] An isomorphism of the K3 surfaces over F_p can be taken to respect the fibrations: for NS = U ⊕ E8 the elliptic fibration is unique, and for U ⊕ E8 ⊕ A2 two fibrations of the same type differ by an isometry of NS preserving the ample cone, which acts by ±1 on Z/3 and so comes from an automorphism over F_p. So the test on Weierstrass equations (§4.3) is the right one. This step I did not check against a source.

**What is assumed and not proved.** That the lattices of Taelman's theorem for Y, X₁, X₂ stand in the relation that van Geemen proves over C. This needs the relative Jacobian of the canonical lift of Y to be the canonical lift of X_i. Proposition E of [note 6](06-torsors-under-cm-class-groups.md) argues in this way for general steps; here ℓ = 2, which note 6 leaves out (its §7 items 1, 2).

### 3.3 Consequences that do not need condition C

For any pair with P irreducible of degree N:

- [D] T(Y) = an F-stable sublattice of index 2 of T(X₁), so F fixes a nonzero vector of T/2T, and its conjugate: **(t − 1)² divides P mod 2.** (The multiplicity of t − 1 is even for every P, because P mod 2 is self-reciprocal.)
- [D] ρ̄(Y) = ρ̄(X₁), so NS(Y) = U(2) ⊕ (the same root lattice), of discriminant 4 times that of NS(X₁). The Artin–Tate formula for Y and for X₁, which have the same polynomial, gives |Br(X₁)| = 4·|Br(Y)|, and the order of a Brauer group is a square [M]. **So |Br(X₁)| is 4 times a square.**
- For X with NS(X) the lattice of the fibre types the formula reads |Br(X)|·|disc| = P(p)/p^{N−1} [M], with disc = 1 for the generic and e8 families. For e8a2: disc = 3 if the two branches of the node at ∞ are rational; if not, the rational part of NS has discriminant 2 and the eigenvalue −p contributes a factor 2, so |Br(X)| = P(p)/p^{N−1} [D].

---

## 4. Results [E]

### 4.1 Validation

| | e8a2, p = 5 | e8a2, p = 7 | e8, p = 5 | generic, p = 5 |
|---|---|---|---|---|
| Pairs | 16,000 | 4,000 | 1,200 | 40 |
| Fibres of Y against fibres of X₁ and X₂, over F_p | 108,763 of 108,763 | 42,407 of 42,407 | 8,165 of 8,165 | 345 of 345 |
| Trace sums of X₁ and X₂ equal | all | all | all | all |
| With a polynomial (middle coefficient not 0) | 6,867 | 1,801 | 533 | 18 |
| Weil bound | all | all | all | all |
| P(p) = 0 | 1,319 | 344 | 130 | 11 |
| Reducible | 3,644 | 891 | 261 | 15 |
| Irreducible, not ordinary | 602 | 133 | 60 | 0 |
| Irreducible and ordinary | 2,621 | 777 | 212 | 3 |

At these primes more than half of the polynomials are reducible, most with the root p or −p: the Picard number has jumped. Those are outside the prediction.

### 4.2 Above 2

Pairs with P irreducible and ordinary, against controls: random Weierstrass equations with the same fibre types, kept by the same test on the discriminant, not coming from a double cover.

| | e8a2, p = 5 | e8a2, p = 7 | e8, p = 5 | generic, p = 5 | control e8a2, p = 5 | control e8, p = 5 |
|---|---|---|---|---|---|---|
| Irreducible and ordinary | 2,621 | 777 | 212 | 3 | 345 | 147 |
| t − 1 does not divide P mod 2 | 0 | 0 | 0 | 0 | 269 | 112 |
| divides exactly twice | 964 | 312 | 73 | 0 | 31 | 18 |
|  of which two primes of degree 1 | 964 | 312 | 73 | — | 22 | 10 |
|  of which one prime of degree 2 | 0 | 0 | 0 | — | 9 | 8 |
| divides 4 times or more | 1,657 | 465 | 139 | 3 | 45 | 17 |
| \|Br\| by Artin–Tate, P(p) ≠ 0: a square | 5,548 of 5,548 | 1,457 of 1,457 | 403 of 403 | 7 of 7 | 608 of 608 | 235 of 235 |
|  4 times a square | 5,548 | 1,457 | 403 | 7 | 214 | 84 |

- The case "one prime of degree 2 above (t − 1)²" is what Lemma 1 (3) excludes for a surface with a rational double cover of this kind: there F ≡ 1 on a part of T/2T of dimension 2, but the two isotropic lines are exchanged by Frobenius. It occurs in a third of the controls with multiplicity 2 and in none of the 1,349 pairs.
- In the controls the number is even exactly when it is 4 times a square, in 35% of the polynomials with P(p) ≠ 0. The 23% of the rows above (t − 1 divides P mod 2) is over irreducible ordinary polynomials only.
- The orders found for pairs with irreducible P are 4, 16, 36, 64, 100, with 4 in two thirds to three quarters of the cases.
- The Artin–Tate count was made for every pair with P(p) ≠ 0. For reducible P the Néron–Severi lattice is larger than the lattice of the fibre types, so there the number is not the order of the Brauer group; it was 4 times a square all the same.
- The generic control, 10 equations, gave one irreducible ordinary polynomial; t − 1 does not divide it mod 2.

### 4.3 The isomorphism test

Isomorphism over F_p is tested on the Weierstrass equations: (A₂, B₂)∘M = (μ⁴A₁, μ⁶B₁) for some M in PGL₂(F_p) and μ in F_p^×, all M tried.

| | e8a2, p = 5 | e8a2, p = 7 | e8, p = 5 |
|---|---|---|---|
| Condition C | 722 pairs, 164 polynomials | 201, 176 | 49, 48 |
| Class of 𝔩/𝔩̄ trivial | 0 | 0 | 0 |
| X₁ ≅ X₂ | 0 | 0 | 0 |
| \|G\|: least, median, greatest | 2, 52, 92 | 5, 188, 580 | 53, 202, 653 |
| Median of \|G\|/p⁴, resp. p⁵ (quartiles) | 0.083 (0.051, 0.115) | 0.078 (0.051, 0.104) | 0.065 (0.045, 0.093) |
| G cyclic | 705 | 167 | 46 |
| The class of 𝔩 generates G | 449 | 93 | 25 |
| Least order of the class | 2 | 5 | 13 |
| Σ 1/\|G\| over the pairs | 21.5 | 1.6 | 0.3 |

- **Two ways of computing the class agree in all 972 cases**: directly, 𝔩/𝔩̄ = (y) and y·ȳ a relative norm of a unit; and through G, the class of 𝔩 in Cl(E)/⟨Cl(E₀), ramified primes⟩. This is a check of Lemma 9 (2) of note 13 at degrees 10 and 12.
- **Σ 1/|G|** is the number of isomorphic pairs to expect if X₂ were a uniformly random member of the torsor of X₁: 23.5 in all, none observed. At p = 5 in the family e8a2 alone the probability of none is 2·10⁻¹⁰. So X₂ is not a random member; it is consistent with X₂ = c·X₁, c ≠ 1.
- **Outside condition C** no pair with irreducible ordinary P is isomorphic either (2,641 pairs). The only isomorphic pairs in all 21,240 are two in the family e8a2 at p = 5: one with reducible P, where the two equations coincide, and one with P not ordinary.
- Condition C fails for the remaining irreducible ordinary pairs because t − 1 divides P mod 2 four times or more (2,264), or, with multiplicity 2 (377): another factor of P mod 2 is repeated (216), the order is not maximal at p (81), an odd prime divides the index (69), or two of these (11).

### 4.4 Partners

For the family e8a2 an isomorphism over F_p fixes 0 and ∞, so the least image of (A, B) under u ↦ cu and (A, B) ↦ (μ⁴A, μ⁶B) is a canonical label. For each polynomial the script forms the graph on the labels with an edge for each pair.

| | Condition C, p = 5 | Condition C, p = 7 | Several pairs of primes above 2, p = 5 |
|---|---|---|---|
| Pairs, polynomials | 722, 164 | 201, 176 | 1,182, 288 |
| Surfaces | 1,194 | 400 | 1,949 |
| Surfaces occurring in 3 or 4 pairs | 23 | 0 | 51 |
| Distinct partners: 1, 2, 3 | 1,108, 86, 0 | 400, 0, 0 | 1,805, 139, 5 |
| Polynomials with more surfaces than \|G\| | 0 | 0 | — |
| Closed cycles | 0 | 0 | — |

- Under condition C no surface has three partners, as §3.2 requires; of the 23 that occur in at least three pairs, 15 have two partners and 8 have one. When t − 1 divides P mod 2 four times or more, three partners occur, and never more than the multiplicity.
- The number of surfaces seen for a polynomial never exceeds |G|, but no polynomial was sampled in more pairs than half of |G|, so this says little. For the one polynomial with |G| = 2 both surfaces were seen and they form the pair.
- The longest chain of steps seen has 6 surfaces. A cycle would have the length of the order of c, 7 or more in all but five cases.

---

## 5. Why the class is never trivial

[D] Suppose 𝔩/𝔩̄ = (x) with x x̄ = 1. Then z = 2x is an integer of E, since (2) = 𝔩𝔩̄𝔟 and (z) = 𝔩²𝔟, and z z̄ = 4: all conjugates of z have absolute value 2, and z is not in E₀. There are finitely many algebraic integers of degree at most 20 with all conjugates of absolute value 2, so only finitely many fields E can have a trivial class; if z generates a quadratic subfield, it is Q(√−7) or Q(√−15), which would be an extra endomorphism.

More generally, if cⁿ is trivial then E contains zₙ with all absolute values 2ⁿ, and if zₙ generates E, |d_E| ≤ 2^{(n+1)N(N−1)}. So the order of c grows at least like log₂|d_E| / N(N − 1) [D, a weak bound].

This is the K3 form of a familiar fact: in an imaginary quadratic field in which 2 splits, a prime above 2 is principal only for the discriminant −7. The classes of the primes used for steps are not uniformly distributed in the group; the small ones have large order. In the data c generates G in 58% of the cases and its least order is 2.

Consequence for the test: "X₁ ≅ X₂ exactly when c is trivial" cannot be exercised by drawing double covers at random. It would need a surface whose field contains such a z, found from the side of the polynomials.

---

## 6. What this does and does not show

**Shown, on surfaces with ρ̄ = 12 and 10, with equations:**

- Pairs of elliptic K3 surfaces in the same isogeny class, not isomorphic, produced by one explicit construction at ℓ = 2.
- The local structure at 2 that the lattice description requires of such a pair (§3.1, §3.3), which random surfaces of the same family do not have.
- The count of partners that the lattice description gives from P mod 2.
- The group G and the class c for 972 pairs; |G| ≈ 0.08·p⁴ and 0.065·p⁵.

**Not shown:**

- That X₂ = c^{±1}·X₁ rather than another element with the same properties. Everything measured is also true of any g ≠ 1 that depends only on (X₁, the Brauer class). No invariant separates the surfaces of a torsor (note 3 §4), so this needs a closed cycle or a count: all surfaces with one polynomial, and all the double covers joining them. In the family e8a2 at p = 5 that is about 4·10⁵ surfaces and 8·10⁶ forms f, which is within reach and was not done.
- A step from a given surface. The construction starts from Y. From X₁ and one of its two Brauer classes one needs the 2-covering, a descent over F_p(u) (note 17 §5.3).
- Anything at an odd prime ℓ.
- The relation of §3.2 between the canonical lifts, which is assumed.
- Picard number 2 beyond three polynomials. The generic family needs point counts over F_{5^10}: 10 to 18 minutes for each surface, and 15 of 18 polynomials were reducible. A larger prime would reduce that share and cost p¹⁰.
- p > 7. The families at p = 5 are small: 164 polynomials account for the 722 pairs of condition C.
- Certified class groups; and the two derivations marked "not checked against a source" (§3.2 on fibrations, §2.3 on the degree of L, which the Weil bound supports).

**For the hardness question** ([../SUMMARY.md](../../SUMMARY.md) §1.1, [note 13](13-hardness-on-paper.md), [note 17](17-below-picard-16.md)). This is the first step computed on K3 surfaces below Picard number 16 in this repository, and it is at one of the two primes, 2 and 3, where explicit geometry for a step is known ([M], [../SUMMARY.md](../../SUMMARY.md) §1.1). It makes the action at 2 testable in these families and does not make it usable: one prime, half of a step, and (§5) a prime whose class has large order but which alone would need walks of length |G|. The verdict of note 17, no reduction known and no usable action, stands.

---

## 7. Check [E]

`ellk3_two_neighbours.py`, with PARI through cypari2.

- **Sampling.** `--family generic|e8|e8a2`, `--p`, `--count`, `--seed`; sample i uses the seed 1000003·seed + i, so runs are reproducible. Forms are redrawn until both Jacobians have the fibre types of the family.
- **Per sample.** The fibre comparison over F_p for both fibrations; the trace sums S_d of X₁ for d ≤ N/2 by Frobenius orbits and of X₂ for d ≤ `--check2` by the direct sum; the polynomial; the Weil bound numerically; the Artin–Tate number; the isomorphism test over all of PGL₂(F_p); the local data of Lemma 1; under condition C the class of 𝔩/𝔩̄ both ways, G, and the order of the class. PARI errors if the series is not integral, if the two primes above 2 are not conjugate, or if a norm-one generator fails its own verification; none occurred.
- **`--control`.** Random Weierstrass equations with the fibre types of the family; the same data without the pair.
- **`--graph FILE`.** The partner graph of §4.4 from a result file of the family e8a2.
- **`--summary FILE`.** The summary of a run from its result file. The `.log` files were written with it after the runs.

**Limits.**

- Class groups under GRH. The maximal order is computed with trial division to 10⁶ and accepted only if `nfcertify` returns nothing; no pair was rejected for this, two control polynomials were.
- The fibre comparison and the equality of the trace sums are consequences of Lang's theorem; they test the formulas for I and J and the code, not the prediction.
- The sign of the functional equation is taken to be + whenever the middle coefficient is not 0. A polynomial with sign − and nonzero middle coefficient does not exist, so this is safe; the Weil bound holds for all.
- Samples whose polynomial has middle coefficient 0 are dropped whatever their sign, which loses a few of sign +.
- The batches at p = 5 of the families e8a2 and e8 and their controls were run with the direct trace sums for X₁; e8a2 at p = 7 and the generic runs with the sums by orbits, after the two were compared.
- The first generic run was stopped and restarted because the direct sums over F_{5^10} were too slow.

---

## Reproduce

From `experiments/`:

```bash
python -I ellk3_two_neighbours.py --family e8a2 --p 5 --count 16000 --workers 10 --check2 5 --out ellk3_e8a2_p5.jsonl
python -I ellk3_two_neighbours.py --family e8a2 --p 7 --count 4000 --workers 10 --check2 5 --out ellk3_e8a2_p7.jsonl
python -I ellk3_two_neighbours.py --family e8 --p 5 --count 1200 --workers 10 --out ellk3_e8_p5.jsonl
python -I ellk3_two_neighbours.py --family e8a2 --p 5 --count 1500 --workers 10 --control --out ellk3_e8a2_p5_control.jsonl
python -I ellk3_two_neighbours.py --family e8 --p 5 --count 600 --workers 10 --control --out ellk3_e8_p5_control.jsonl
python -I ellk3_two_neighbours.py --family generic --p 5 --count 10 --workers 5 --check2 7 --bnf-limit 900 --out ellk3_generic_p5.jsonl
python -I ellk3_two_neighbours.py --family generic --p 5 --count 30 --seed 2 --workers 10 --check2 7 --bnf-limit 900 --out ellk3_generic_p5_seed2.jsonl
python -I ellk3_two_neighbours.py --family generic --p 5 --count 10 --workers 5 --control --bnf-limit 60 --out ellk3_generic_p5_control.jsonl
python -I ellk3_two_neighbours.py --graph ellk3_e8a2_p5.jsonl
python -I ellk3_two_neighbours.py --graph ellk3_e8a2_p7.jsonl
```

Each run prints its summary; `--summary FILE` prints it again. About 2.3, 5.4 and 0.8 core-hours for the first three, 3.7 for the two generic runs together, and minutes for the controls except the generic one, about 2.5 core-hours. The two `--graph` outputs are in `ellk3_e8a2_graph.log`.

## Sources

- van Geemen, *Some remarks on Brauer groups of K3 surfaces*, https://arxiv.org/abs/math/0408006 : 1.2, 4.1, 4.7, 4.8, 5.3, 5.5, read for note 17 [P].
- Taelman, *Ordinary K3 surfaces over a finite field*, https://arxiv.org/abs/1711.09225 : Theorem 4.6, read for note 17 [P].
- Recalled, not re-read [M]: the invariants of a binary quartic and the Jacobian of y² = quartic (Weil; Cremona); Lang's theorem; the L-function of an elliptic curve over a function field and its degree (Grothendieck, Ogg, Shafarevich); Kodaira's fibre types and Tate's algorithm in characteristic at least 5; the Artin–Tate formula and the squareness of the order of the Brauer group (Milne; Liu, Lorenzini, Raynaud); the singularities E8 and A2 of plane curves by weights; the uniqueness of the elliptic fibration when NS = U ⊕ E8.
- In this repository: notes [3](03-torsor-action.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [17](17-below-picard-16.md); [../SUMMARY.md](../../SUMMARY.md) §1.1.

