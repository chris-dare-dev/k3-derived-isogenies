# 17. Picard number below 16

*A reading of the literature; nothing is computed beyond the arithmetic of the tables. It takes up the case that [note 13](13-hardness-on-paper.md) §1 item 5 leaves aside: Picard numbers below 16, where no transfer to abelian varieties is known. Builds on notes 2, 6 and 13 to 16. Labels: **[P]** read in the cited source for this note; **[M]** recalled, not checked; **[D]** derived here; **[H]** heuristic. Statements taken from other notes of this repository are cited by note and keep the label they have there. Not re-checked.*

Throughout, X is an ordinary K3 surface over F_q, ρ̄ its geometric Picard number, r = 22 − ρ̄ the rank of the transcendental part, E = Q[F] the CM field of degree r and m = r/2. Since E is a CM field, r is even, so **ρ̄ is even** [D, from Taelman as quoted in note 6 §2]: the values below 16 are 14, 12, 10, 8, 6, 4, 2.

---

## 1. Result

**The question.** For ρ̄ < 16, could the class-group action be a hardness assumption that does not reduce to abelian varieties?

**Answer.** No reduction to abelian varieties is known in any full family with ρ̄ ≤ 14. That is not a hardness assumption, for three reasons.

1. **Nobody can evaluate the action.** No step at any prime ℓ has been computed on a K3 surface with ρ̄ < 16, and no formula for one exists. The attacker and the honest party lack the same operation (note 13 §2). Explicit geometry for a step exists only at ℓ = 2 and 3 (starting hypothesis (d) of [../SUMMARY.md](../../SUMMARY.md) §1.1 [M]; see §5.3). The case closest to a step inside one family is that of elliptic K3 surfaces, over C, and it is half of a step (§5.3).
2. **The reduction exists in principle for every ρ̄.** The Kuga–Satake abelian variety of a K3 surface over a field of odd characteristic exists without any hypothesis on ρ̄, as a morphism of moduli spaces over Z[1/2] (Madapusi Pera) [P]. What is missing below 16 is equations for that morphism, not the morphism. "Does not reduce" can only mean "nobody has the formulas" (§3).
3. **The ceiling is CSIDH's whatever happens.** The group is abelian for every family (note 6, Proposition A), so Kuperberg's algorithm applies given the action [M], and by the square-root count an instance has about 512 bits for a group of 2²⁵⁶ elements at every ρ̄ (§6) [H]. A new assumption could at best equal CSIDH's security with steps that do not exist.

**Verdict by range.**

| ρ̄ | Verdict | Basis |
|---|---|---|
| 20, 18 | Reduces to elliptic curves and abelian surfaces | Starting hypothesis (b) in [../SUMMARY.md](../../SUMMARY.md) §1.1 [M]; not re-examined here |
| 16 | Reduces to abelian fourfolds in the six-line family, tested on point counts | Notes 13, 14, 15, 16. For other families with ρ = 16 the correspondence is algebraic over C in some (§2.1) and nothing was done over F_q |
| any, inside a family with an automorphism of order 4 or 3 acting on H^{2,0} | Reduces, over C, to abelian varieties of dimension m = (22 − ρ)/2 with multiplication by Q(i) or Q(√−3); an algebraic correspondence is known for cyclic quartics (ρ = 8) and for one family with ρ = 2 | §2.2 [P]. Not examined over F_q |
| 14 | No reduction known. No usable action. A transfer to abelian 8-folds is the one case below 16 where a computation on the abelian side is conceivable | §2.3, §4 |
| 12 | No reduction known. No usable action. The abelian varieties have dimension 16 | §4 |
| 10, 8, 6, 4, 2 (full families) | No reduction known, and the abelian varieties, of dimension 32 to 512, cannot be written down by theta functions. No usable action, so nothing to assume | §4, §5 |
| Open | Whether the steps at ℓ = 2, 3 on elliptic K3 surfaces can be computed over F_q (§5.3); whether the 8-folds at ρ̄ = 14 have a geometric construction (§2.3); whether the class can be read off a surface (note 13 §7, unchanged) | |

**Two facts used below.**

- Taelman's Theorem 4.6 has three unconditional cases [P]: an ample class of degree below p − 4; a hyperbolic plane in the Picard lattice and p ≥ 5; ρ̄ ≥ 12 and p ≥ 5. The first two hold at any Picard number. So the classification is unconditional for elliptic K3 surfaces with a section and for double planes when p ≥ 7 (§5.1). Note 13 §9 states all three.
- The Kuga–Satake variety of a transcendental lattice of rank r has dimension 2^{r−2} = 2^{20−ρ} [P: van Geemen's survey, 8.1]. What matters for a transfer is a simple factor, which for a surface over F_q has dimension 2^{m−1} = 2^{(20−ρ̄)/2}, or half of that (§4.1) [D].

---

## 2. Where the Kuga–Satake correspondence is known to be algebraic

The Kuga–Satake–Hodge conjecture asks for a cycle on X × A × A inducing T(X)_Q ⊂ H²(A × A, Q), with A a Kuga–Satake variety [P: van Geemen's survey, 10.2]. All statements in this section are over C.

### 2.1 Picard number 16 and above

- Floccari (2023): the conjecture "is known to hold in many cases for K3 surfaces of Picard number at least 17" (he cites Morrison) "but it is wide open otherwise. There are a couple of families of K3 surfaces of general Picard rank 16 for which the conjecture is known, namely, the family of double covers of the plane branched at 6 lines" (Paranjape) "and the family of K3 surfaces with 15 nodes in P⁴" (Ingalls, Logan, Patashnick) [P]. Note "many cases", not all, above 16; the statement "ρ ≥ 17 through abelian surfaces" of the starting hypotheses (../SUMMARY.md §1.1, [M]) was not re-examined.
- Ingalls, Logan, Patashnick: the 15-nodal surfaces are covered, up to maps of finite degree, by the square of a curve of genus 7; the correspondence is explicit in the sense used in this repository: curves and maps with equations [P: abstract].
- Voisin (2021), using Markman and O'Grady: the correspondence is algebraic for hyper-Kähler varieties of generalised Kummer type, whose Kuga–Satake varieties are powers of abelian fourfolds of Weil type [P: abstract and Theorem 0.1 as quoted]. Floccari deduces it for the K3 surfaces S with T(S)_Q inside Λ_{Kum₃}(2) ⊗ Q: "infinitely many new families of K3 surfaces of general Picard rank 16" [P: Theorem 1.3]. These cycles come from moduli of sheaves on hyper-Kähler sixfolds. They are algebraic and not explicit: no curve, no equations.
- Varesco: if the correspondence is algebraic for a family of generic Picard number 16, the Hodge conjecture holds for all powers of its members [P: abstract]. This is a consequence, not a new case.
- Ingalls, Logan, Patashnick, Remark 2.25: for a 5-dimensional family of generic Picard rank 15 in their construction they "do not know whether this is to be expected" [P].

So for families cut out by a Picard lattice alone, the known results stop at 16.

### 2.2 Below 16: only surfaces with a special automorphism

- **Cyclic quartics, ρ = 8.** X: t⁴ = F(x, y, z), with F = 0 a smooth plane quartic C. The transcendental part V has dimension 14 and the primitive algebraic part dimension 7, so ρ = 8 for general F; the automorphism t ↦ it makes V a vector space of dimension 7 over Q(i) [P: van Geemen–Izadi 4.2]. The family has dimension 6 and its periods fill a 6-ball (Kondo, as reported there).
  - The *half twist* V_{1/2} is a Hodge structure of weight one on the same space: an abelian variety A_C of dimension 7 with an action of Q(i) of type (1, 6), "Kondo's abelian variety", and V ⊂ V_{1/2} ⊗ H¹(E_i) with E_i the elliptic curve with multiplication by Z[i] [P: van Geemen–Izadi 4.2, 5.1; van Geemen, *Half twists*, introduction].
  - The correspondence is algebraic and explicit: the tangent lines of C give an isotrivial elliptic surface over C which, after a base change C′ → C, is a product E_i × C′ mapping onto X, and V_{1/2} sits in H¹(C′) [P: van Geemen's survey, Example 11.3; van Geemen–Izadi, Theorem 5.2 for d = 4, k = 2].
  - The whole Kuga–Satake variety, of dimension 2¹², is a sum of pieces S_0, …, S_7 with dim_Q S_i = 2·C(7, i); only S_0 = H¹(E_i) and S_1 = V_{1/2} are needed to carry V [P: survey, Theorem 11.2, stated there with proof "to appear"].
- **Triple covers of a quadric branched along a curve of genus 4, ρ = 2.** S: f₂(x₀..x₃) = 0, g₃(x₀..x₃) + x₄³ = 0 in P⁴. Bolognesi and Laterveer prove the Kuga–Satake–Hodge conjecture for these, through cyclic cubic fourfolds and van Geemen–Izadi; they cover two 9-dimensional components of the moduli of K3 surfaces with a non-symplectic automorphism of order 3, whose general member has Picard number 2 [P: abstract, Theorem 4.4, Remark 6.8]. Their Remark 4.3 lists the known cases as: six lines, 15 nodes, cyclic quartics, "some other K3 surfaces with Picard number 16" (Floccari) [P].
- Voisin's earlier example with an automorphism of order 3 is described in van Geemen's survey 11.1 [P]; her paper was not read.

In both families the transcendental lattice is a hermitian space over an imaginary quadratic field K and the surface moves in a family of dimension m − 1, half the dimension 20 − ρ of the family with the same Picard lattice. These are not the general surfaces of their Picard number.

### 2.3 Picard number 14: the abelian side is identified, the correspondence is not

- Poon (2025): for a Picard lattice P of rank 14 whose determinant is a square, the Kuga–Satake variety of a very general member is isogenous to A₊⁴ × A₋⁴ with A₊, A₋ non-isogenous simple abelian 8-folds whose endomorphism algebra is a totally definite quaternion algebra; X ↦ A₊ gives a map from the moduli space of P-polarised K3 surfaces to a moduli space of such 8-folds, induced by the isomorphism of the domains of type IV₆ and II₄ [P: abstract, Theorem 3.1.1]. The construction is on lattices and periods (worked in Magma for T = U ⊕ U(2) ⊕ D₄(−1)). No cycle, curve or equation is given, and none is claimed.
- Koike, Shiga, Takayama, Tsutsui (2000) study the family with P = D₄³ ⊕ ⟨−2⟩ ⊕ ⟨2⟩: double covers of P¹ × P¹ branched along four curves of bidegree (1, 1). They present it as the first step towards extending to rank 8 the theta-constant description that Matsumoto gave for six lines, and state: "we do not discuss the representation of the inverse of the period mapping" [P: introduction].
- I found no statement that the correspondence is algebraic for any family of generic Picard number 14, and Floccari's and Bolognesi–Laterveer's lists (2023, 2024) contain none [P for the lists; the absence is as of those lists].

### 2.4 Summary by Picard number

| Generic ρ of the family | Algebraic correspondence known | Explicit (curves, equations) |
|---|---|---|
| ≥ 17 | "many cases" (Morrison, via Floccari's summary) [P] | Shioda–Inose and Kummer structures [M] |
| 16 | six lines; 15 nodes; Floccari's families [P] | six lines (Paranjape; note 14), 15 nodes (squares of genus-7 curves) [P]; Floccari's: no |
| 15 | none; one family named as unknown [P] | — |
| 14 | none found; target 8-folds identified on periods [P] | no |
| 13 to 9 | none found | — |
| 8 | cyclic quartics only (6 moduli out of 12) [P] | yes: product of E_i and a cover of the quartic curve [P] |
| 7 to 3 | none found | — |
| 2 | genus-4 triple covers only (9 moduli out of 18) [P] | through cubic fourfolds; not read in detail |

---

## 3. What a transfer needs

### 3.1 The correspondence need not be explicit, and need not be a cycle

- **Existence over finite fields is unconditional.** "Given any field k of odd characteristic p and a polarized K3 surface (X, ξ) over k, there exists a finite separable extension k′/k and an abelian variety A over k′, the Kuga–Satake abelian variety, such that the Z_ℓ and crystalline realizations of the primitive cohomology PH²(X, ξ) embed naturally within those of H¹(A) ⊗ H¹(A)" [P: Madapusi Pera, Theorem 3]. The construction "gives rise, away from characteristic 2, to an open immersion from the moduli of primitively polarized K3 surfaces (of any fixed degree) to a certain regular integral model for a Shimura variety of orthogonal type", over which the Kuga–Satake abelian scheme lives [P: abstract, §3]. It is obtained by lifting, applying the complex construction and reducing (Deligne), an "indirect construction" in Rizov's words [P: both abstracts].
- **Isogenies match.** Yang: "isogenies between Kuga-Satake abelian varieties induce isogenies between K3 surfaces, in the context of integral models of Shimura varieties" [P: abstract]. Bragg–Yang: "the derived prime-to-p isogeny classes of K3's match up with the notion of prime-to-p Hecke orbit on the period domains of Kuga-Satake morphisms" [P: introduction].
- **Consequence.** [D, sketch] For a fixed family, X ↦ A(X) is a fixed morphism between quasi-projective schemes over Z[1/2] (with level structures), so it is given by finitely many polynomials that do not depend on q. Evaluating it at a point over F_q costs time polynomial in log q, with a constant depending on the family. Together with the matching of isogenies this says: vectorisation for the K3 surfaces of a fixed family reduces to vectorisation for abelian varieties with complex multiplication, of a fixed dimension, by an algorithm that exists and that nobody can write down. I did not check that the morphism is equivariant for the group of Proposition A at the level of the lattices of Taelman and Deligne; for six lines this is note 15, on lattices only.
- So the Hodge conjecture is not what stands between the K3 problem and abelian varieties. An algebraic cycle is not necessary (the morphism exists without it) and not sufficient (Floccari's cycles give nothing to compute with).

### 3.2 What an attacker would need

1. **The map as an algorithm.** From the equations of X over F_q, in time polynomial in log q, an abelian variety B(X) in the isogeny class predicted by P. Three kinds of source are known:
   - a geometric construction: a curve built from X whose Jacobian or Prym contains B. Exists for six lines (note 14), 15 nodes, cyclic quartics;
   - modular equations: coordinates of X and theta constants of B as functions on the same Shimura variety. Exists for six lines (Matsumoto, as reported by Koike et al. [P]); not for rank 14 [P];
   - lifting: the canonical lift, complex periods, the Clifford algebra, then reduction. This needs the embedding ι and a root of a class polynomial of degree |G| (note 13 §7), so it costs at least as much as the problem [H].
2. **The right point, not the right isogeny class.** The Frobenius polynomial of B follows from P (for six lines, Proposition Y of note 13; in general the eigenvalues are √q·∏(α_j/q)^{±1/2} over the m conjugate pairs [D]). So the isogeny class of B is free, and agreement of point counts tests only that. The secret is the class inside the torsor. Notes 14 and 16 are tests of this first kind.
3. **Equivariance.** A homomorphism from the group acting on surfaces to the class group acting on the B's, with small kernel, such that the map commutes with the actions. Expected from §3.1; checked on lattices for six lines (note 15).
4. **Arithmetic on the abelian side.** A representation of B and isogeny formulas for the primes used (§4.2).
5. **Not needed:** a cycle; a proof (a map verified on point counts and on one step would serve an attacker); definition over F_q rather than over an extension of bounded degree.

### 3.3 What the transfer is worth

Items 1 to 4 give the attacker the *action*: on the abelian side a step is an isogeny, which can be computed in principle. They do not give an attack. No algorithm for vectorisation on abelian varieties is known beyond the generic ones, which need only an oracle for the action and apply equally on the K3 side [M]. So a transfer is worth to the attacker exactly what it is worth to the honest party: it supplies the steps. This is why "reduces to abelian varieties" answered the *novelty* question at ρ̄ ≥ 16 and never the *security* one.

Possible extra value on the abelian side, none established: pairings and genus characters for the decision problem (note 13 §6), and whatever is known about the endomorphism ring.

---

## 4. Does the dimension make the transfer useless

### 4.1 The dimensions

- The Kuga–Satake variety of a lattice of rank n has dimension 2^{n−2} [P: survey 8.1]. Madapusi Pera's, built on the full Clifford algebra of the primitive lattice of rank 21 [P: §3], has dimension 2²⁰ [D] whatever ρ̄ is.
- Simple factors, when the Mumford–Tate group is as large as possible [P: survey, Theorem 7.7 and 8.1]: for n = 2m with d = (−1)^m·disc, the even Clifford algebra is M_{2^{m−2}}(D) × M_{2^{m−2}}(D) with D a quaternion algebra over Q when d is a square, and M_{2^{m−2}}(D) with D a quaternion algebra over Q(√d) otherwise; for n = 2m + 1 it is M_{2^{m−1}}(D). So a simple factor has dimension between 2^{m−2} and 2^m [D]. Known values: 4 for ρ = 16 (Lombardo–Peters–Schütt, note 13 §4.1); 8 for ρ = 14 with square determinant [P: Poon].
- **Over F_q.** [D] The surface has complex multiplication by E, of degree 2m. As a representation of the torus of E, H¹ of the Kuga–Satake variety is a sum of copies of the spin representation, whose 2^m weights (±½, …, ±½) each occur once. When the Galois group of E is the full group (Z/2)^m ⋊ S_m these weights form one orbit, so the simple factor B has complex multiplication by a field of degree 2^m and

  > dim B = 2^{m−1} = 2^{(20−ρ̄)/2};

  when the Galois group preserves the two half-spin representations there are two factors of dimension 2^{m−2}. For ρ̄ = 16 this gives 4 and the field of degree 8 of note 13, Proposition Y. Not checked against a source.
- So log₂ dim B equals log_q |G| under the heuristic of note 13: the abelian variety grows exponentially in the exponent of the group size.
- In a family with an automorphism (§2.2) the relevant abelian variety is the half twist, of dimension m = (22 − ρ̄)/2 and with multiplication by the same field E [D: the Frobenius eigenvalues on V are those of the half twist times the Frobenius of the elliptic curve with multiplication by K]. For cyclic quartics, 7 instead of 64. A cyclic quartic can be ordinary only for p ≡ 1 mod 4 [D, from the Hasse invariant of t⁴ = F; not checked]; I expect p ≡ 1 mod 3 for the family of order 3 and did not derive it.

### 4.2 What can be written down

Cost model: a point of level 4ℓ on an abelian variety of dimension g is represented by g(g + 1)/2·4^g coordinates after compression [P: Lubicz–Robert, abstract]; theta null points have 2^g coordinates at level 2 and 4^g at level 4, and an isogeny with kernel (Z/ℓ)^g costs of the order of ℓ^g operations [M]. The figures are arithmetic [D].

| ρ̄ | m | dim B | theta constants, level 2 / level 4 | 3^g |
|---|---|---|---|---|
| 16 | 3 | 4 | 16 / 256 | 81 |
| 14 | 4 | 8 | 256 / 65,536 | 6,561 |
| 12 | 5 | 16 | 65,536 / 4.3·10⁹ | 4.3·10⁷ |
| 10 | 6 | 32 | 4.3·10⁹ / 1.8·10¹⁹ | 1.9·10¹⁵ |
| 8 | 7 | 64 | 1.8·10¹⁹ / 3.4·10³⁸ | 3.4·10³⁰ |
| 6 | 8 | 128 | 3.4·10³⁸ / — | 10⁶¹ |
| 4 | 9 | 256 | 10⁷⁷ / — | 10¹²² |
| 2 | 10 | 512 | 10¹⁵⁴ / — | 10²⁴⁴ |

- ρ̄ = 14: an 8-fold is within reach of theta-function arithmetic. If a construction of A₊ from the surface were found, this case would go the way of ρ̄ = 16.
- ρ̄ = 12: storable, with steps of the order of 10⁷ to 10¹¹ field operations for ℓ = 3 to 5.
- ρ̄ ≤ 10: B cannot be represented by theta constants. A transfer would have to present B inside the Jacobian of a curve attached to X, as Paranjape's Prym does at 16 and E_i × C′ does for cyclic quartics. Nothing of the kind is known in a full family.

### 4.3 What this means for a hardness claim

- **Asymptotically the dimension is a constant.** ρ̄ is fixed by the family while q grows, so 4^g and ℓ^g with ℓ polynomial in log q are polynomial in log q. A claim "the transfer is infeasible because g = 64" is a statement about constants of the order of 10¹⁹ to 10¹⁵⁴, not about complexity.
- **The claim would be fragile.** It would rest on the absence of a compact description of a specific abelian variety that is known to exist and to be algebraic in X. Paranjape's fourfold is a Prym of a curve of genus 5 and not a point given by 256 theta constants; a similar finding at lower ρ̄ would remove the claim without touching any hard problem.
- **And it protects nothing.** By §3.3 the transfer supplies steps and no attack. If the steps existed on the K3 side the generic attacks would run there; while they do not exist, there is no instance to attack. In neither case does the size of B enter the security.
- So the honest formulation for ρ̄ ≤ 14 is: *the vectorisation problem is that of a commutative group action whose known relation to abelian varieties is not effective*. That is weaker than "new", and it has no content until someone can compute the action.

---

## 5. The honest party

Needs, as in note 13 §2: to name a surface, to know P, to apply a step.

### 5.1 The action on paper: Taelman and the propositions of note 6

- Taelman, Theorem 4.6 [P]: X ↦ (M, F, 𝒦) is an equivalence, without the hypothesis (⋆), between
  - (i) surfaces with an ample L, L² < p − 4, and triples with λ in M ∩ 𝒦, λ² < p − 4;
  - (ii) for p ≥ 5, surfaces whose geometric Picard lattice contains a hyperbolic lattice, and triples whose NS(M, F) does;
  - (iii) for p ≥ 5, surfaces with ρ̄ ≥ 12, and triples with rank NS(M, F) ≥ 12.
- Proposition A of note 6 uses (iii). Its proof uses only that the action leaves N and 𝒦 untouched (note 6, proof of A (1)), so it goes through under (i) or (ii) [D, not re-read line by line]. Hence the torsor structure holds, on the terms of note 6, for:
  - ρ̄ = 14, 12: every family, p ≥ 5, as written;
  - ρ̄ ≤ 10: elliptic K3 surfaces with a section, p ≥ 5; double planes for p ≥ 7; quartics for p ≥ 11.
- Proposition F (one step is one twisted equivalence) quotes Huybrechts–Stellari's corollary for ρ ≥ 12. Below 12 it is not available in that form. For elliptic surfaces the equivalence between a genus-one fibration and its Jacobian with a Brauer class is classical (Căldăraru; Donagi–Pantev) [M].
- Existence of surfaces with a given polynomial: Ito proves it for p ≥ 7 and every m from 1 to 10, that is ρ̄ from 20 down to 2, "up to finite extensions of the base fields", by constructing a complex surface with complex multiplication "which admits an elliptic fibration with a section, or an ample line bundle of low degree" [P: abstract, Theorem 1.3]; Taelman's earlier theorem assumed semistable reduction [P]. Kedlaya and Sutherland note that this leaves "an uncontrolled base field extension" [P: abstract]. Theorem 4.6 is the statement without extension, in its three cases.

So the status of Taelman's theorem is not what separates ρ̄ < 16 from ρ̄ = 16. On paper the action is there for ρ̄ = 14 and 12 exactly as for 16, and for elliptic surfaces down to ρ̄ = 2.

### 5.2 Naming a surface and computing its zeta function

- **Naming.** Equations are easy: quartics, double planes, Weierstrass equations y² = x³ + a₈(t)x + b₁₂(t), double covers of P¹ × P¹ branched along four (1, 1)-curves (ρ = 14, Koike et al. [P]).
- **But the torsor fixes N = NS(X̄), and ρ̄ is even.** [D] In a family whose generic Picard number is odd (quartics, double planes: 1) every surface over F_q has extra classes, and N changes from surface to surface. The honest party would have to find N for each surface. A family with even generic Picard number avoids this: elliptic surfaces with a section (N = U, ρ̄ = 2 for the general member [P: van Geemen, §3, over C]), or a family polarised by a lattice of even rank.
- **Zeta function.** For a nondegenerate hypersurface in a toric variety over F_p the controlled-reduction method is quasi-linear in p; for the Dwork quartic at p = 2²⁰ − 3 it took 22h7m on one core with the full simplex and 1m33s after reducing the polytope; examples go up to p ≈ 10⁶ [P: Costa–Harvey–Kedlaya]. Harvey's algorithm computes the zeta function of an arbitrary variety in time p^{1/2+o(1)}, with no implementation for surfaces reported [P: Harvey, Theorem 1.1; Costa–Harvey–Kedlaya, introduction]. For q = p^f with p small the p-adic methods are polynomial in f [M].

### 5.3 A step: the one lead, elliptic surfaces at ℓ = 2 and 3

Read in van Geemen, *Some remarks on Brauer groups of K3 surfaces* [P], over C:

- For an elliptic K3 surface X → P¹ with a section, Br(X) is the Tate–Shafarevich group. An element α of order b is a genus-one fibration Y → P¹ without a section, with a multisection of degree b and relative Jacobian X; T_Y = ker α (his 1.2, §3, 4.1).
- For general X, NS(Y) is a lattice Λ_{b,c} of rank 2, and "in general Y has two genus one fibrations and … the relative Jacobians of these fibrations are not always isomorphic" (4.1, 4.7).
- (b, c) = (2, 0): Y is the double cover of P¹ × P¹ branched along a curve of bidegree (4, 4); the two projections are the two fibrations; their Jacobians are not isomorphic in general, because the discriminant group (Z/2)² of T_Y has two maximal isotropic subgroups (4.8, 5.5).
- (b, c) = (3, 0): Y is the intersection of a quadric of rank 4 and a cubic in P⁴, and the two rulings of the quadric cut out the two fibrations in plane cubics (5.7).
- The 2-torsion of Br(X) is the 2-torsion of the Jacobian of a curve on X, by Recillas' trigonal construction (§7).

Reading [D]: for (b, 0) the two Jacobians X₁, X₂ have T(X₁) ∩ T(X₂) = T(Y) of index b in each. They are Kneser b-neighbours, which is a step of the action ([note 7](07-ring-changing-steps.md), Proposition H). So:

- **The upper half of a step is explicit at b = 2 and b = 3 in a family with ρ = 2.** From Y to X₁ and X₂ one takes the Jacobian of a binary quartic or of a ternary cubic over the function field, by classical invariant theory [M]. A random (4, 4) curve over F_q gives a pair of 2-neighbours with equations.
- **The lower half is a descent problem.** Given X = X₁ and α one must construct Y: a 2-covering or 3-covering of the generic fibre over F_q(t). For b = 2 the classes α over F_q correspond to rational 2-torsion points on the Jacobian of a trigonal curve, of genus 10 for a general surface [D, from his §7]. I did not find this carried out for K3 surfaces.
- **It does not make a usable action.** Two primes cannot generate the group efficiently: with k generators, walks of length N reach at most (2N + 1)^k elements, so k = 2 needs N of the order of |G|^{1/2} [D]. CSIDH uses 74 primes [M]. And ℓ = 2 is outside Propositions G and H of note 6, which assume ℓ odd.
- **It is the only place found where a step with ρ̄ < 16 could be computed and compared with the lattice prediction**: P from point counts, the genus and G from the lattices as in note 2, and one pair of 2-neighbours from one (4, 4) curve. The group has about q⁹ elements [H], so at ρ̄ = 2 no torsor could be enumerated (with an E8 point and a cusp on the branch curve the family at p = 5 has about 4·10⁵ surfaces, [note 22](22-two-neighbours-of-elliptic-k3.md)); a test would be of one step, not of the action.

For general b the model of Y_{b,0} lies in P^{b+1} on a quadric of rank 4 (his 5.3); equations for the Jacobians of genus-one curves of degree b ≥ 5 were not looked up.

The case b = 2 is tested over F_5 and F_7 in [note 22](22-two-neighbours-of-elliptic-k3.md), on 21,240 double covers with ρ̄ = 12, 10 and 2. The two Jacobians have the same polynomial and are never isomorphic when the polynomial is irreducible and ordinary; t − 1 divides P mod 2 at least twice, the primes above 2 are a split pair when it divides exactly twice, and the Brauer order is 4 times a square, in every case, none of which holds for random surfaces of the same families. Under a condition that forces the maximal order (972 pairs) the class of 𝔩/𝔩̄ is never trivial and no surface has more than two partners. Not tested: that the step is that class and not another element; no cycle was closed and no step was taken from a given surface. The group has about 0.08·p⁴ elements at ρ̄ = 12 and 0.065·p⁵ at ρ̄ = 10.

### 5.4 Answer to "is there any family below 16 where all three can be done"

| Family | ρ̄ | Name a surface with known N | P | A step |
|---|---|---|---|---|
| Lattice-polarised, rank 14 (e.g. four (1,1)-curves) | 14 | yes | yes, small p or p-adic in f | no |
| Cyclic quartics | 8 generic | yes | yes | no; on the abelian side, an isogeny of a 7-fold with multiplication by Z[i], not written |
| Elliptic with section | 2 generic | yes | yes | half of one, at ℓ = 2 and 3 [P over C] |
| Quartics, double planes | 2 or more, N varies | N must be found | yes | ℓ = 2 down-steps to other families (starting hypothesis (d), ../SUMMARY.md §1.1) [M] |

No family has all three.

---

## 6. Sizes and parameters

All under the heuristic |G| ≈ q^{(20−ρ̄)/2} of note 13 §3.3, which has been measured at ρ̄ = 16 only. [H] throughout; the arithmetic is [D].

| ρ̄ | r | \|G\| | q for \|G\| ≈ 2²⁵⁶ | Moduli, bits to name a surface | Zeta function over a prime field of that size |
|---|---|---|---|---|---|
| 16 | 6 | q² | 2¹²⁸ (2¹³¹ with the constant 0.01 of note 13) | 4, 512 | out of reach; small characteristic needed |
| 14 | 8 | q³ | 2^85.3 | 6, 512 | out of reach, also at p^{1/2} |
| 12 | 10 | q⁴ | 2⁶⁴ | 8, 512 | p^{1/2} = 2³² if such an algorithm were implemented |
| 10 | 12 | q⁵ | 2^51.2 | 10, 512 | p^{1/2} ≈ 2^25.6 |
| 8 | 14 | q⁶ | 2^42.7 | 12, 512 | p^{1/2} ≈ 2^21.3 |
| 6 | 16 | q⁷ | 2^36.6 | 14, 512 | linear in p: about 2·10⁶ core-hours |
| 4 | 18 | q⁸ | 2³² | 16, 512 | linear in p: about 9·10⁴ core-hours |
| 2 | 20 | q⁹ | 2^28.4 | 18, 512 | linear in p: about 8·10³ core-hours |

- **Instance length is always about 512 bits.** The family has 20 − ρ̄ moduli and the group is the square root of the number of points, so (20 − ρ̄)·log₂ q = 2·log₂|G|. This is CSIDH-512's ratio at every Picard number: lowering ρ̄ gains nothing in size.
- **The field shrinks as ρ̄ drops.** This is the one thing that gets easier. At ρ̄ = 16 a prime field of 2¹³¹ elements puts P out of reach (note 13 §5). At ρ̄ = 2 the field has about 2²⁸ elements and a quasi-linear method is within a few thousand core-hours, by linear extrapolation from 22h7m at p ≈ 2²⁰ for a generic-shape quartic; the extrapolation ignores the loss the authors report once the working modulus exceeds 2⁶² and assumes a Weierstrass model costs the same.
- **Attacks, given the action:** √|G| = 2¹²⁸ classically, Kuperberg quantumly, as in note 13 §5, at every ρ̄.
- **Automorphism families.** [H, D] There the family has m − 1 moduli and E = K·E₀, and the same square-root count gives |G| ≈ q^{(m−1)/2} = q^{(20−ρ̄)/4}, the square root of the general value: q³ for cyclic quartics, q^4.5 for the genus-4 family. Derived on paper for m = 2 only, by counting the polynomials X² − aX + uq with u a unit of K; not measured.
- **Per step**, on the abelian side: ℓ^g with g from the table of §4.2, against ℓ⁴ at ρ̄ = 16 and ℓ for CSIDH.

---

## 7. What this does not give

- Any computation. No surface was named here, no zeta function computed, no step taken. [Note 22](22-two-neighbours-of-elliptic-k3.md) does this at ℓ = 2 for pairs of elliptic K3 surfaces with ρ̄ = 12 and 10 (and three with ρ̄ = 2): surfaces named, zeta functions computed, and the pair (surface, its neighbour at 2) written down from a double cover. The verdict of §1 is unchanged by it.
- A check of §3.1's "in principle" reduction: equivariance of the Kuga–Satake morphism for the group of Proposition A, the level structures, the degree of the extension k′/k.
- The dimension formula 2^{m−1} of §4.1 from a source, or the half-twist statement over F_q.
- Proposition A under Taelman's cases (i) and (ii), written out; Propositions E to H below ρ̄ = 12.
- The descent of §5.3 over F_q, or any verification that the two Jacobians of a (4, 4) double cover are neighbours in the sense of note 6 at ℓ = 2, where note 6 is silent.
- The heuristic group size at any ρ̄ other than 16.
- A search of the literature after 2025 beyond the arXiv listing for "Kuga-Satake" (48 entries) taken on 2026-11-14.
- The case ρ̄ ≥ 17 over C, where Floccari says "many cases" and starting hypothesis (b) (../SUMMARY.md §1.1) says all [M].

---

## Sources

Read for this note (the passages named), all [P]:

- van Geemen, *Kuga-Satake varieties and the Hodge conjecture*, https://arxiv.org/abs/math/9903146 : §5.1, Theorem 7.7, §8.1, §10, §11.
- van Geemen, Izadi, *Half twists and the cohomology of hypersurfaces*, https://arxiv.org/abs/math/0008170 : abstract, introduction, 4.1–4.2, 5.1–5.4.
- van Geemen, *Half twists of Hodge structures of CM-type*, https://arxiv.org/abs/math/0008076 : abstract, introduction, 1.1–1.4.
- van Geemen, *Some remarks on Brauer groups of K3 surfaces*, https://arxiv.org/abs/math/0408006 : abstract, introduction, 1.2, 4.1–4.8, 5.1–5.9.
- Poon, *Kuga-Satake construction on families of K3 surfaces of Picard rank 14*, https://arxiv.org/abs/2504.08270 : abstract, introduction, 2.1, 3.1 to Theorem 3.1.1 and its proof, opening of §4.
- Koike, Shiga, Takayama, Tsutsui, *Study on the family of K3 surfaces induced from the lattice (D₄)³ ⊕ ⟨−2⟩ ⊕ ⟨2⟩*, https://arxiv.org/abs/math/0003145 : introduction.
- Floccari, *Sixfolds of generalized Kummer type and K3 surfaces*, https://arxiv.org/abs/2210.02948 : abstract, introduction to Theorem 1.3.
- Voisin, *Footnotes to papers of O'Grady and Markman*, https://arxiv.org/abs/2106.06979 : abstract, introduction to Theorem 0.1.
- Bolognesi, Laterveer, *A 9-dimensional family of K3 surfaces with finite-dimensional motive*, https://arxiv.org/abs/2310.11981 : abstract, introduction, Remarks 2.2, 4.3, 6.8, Theorem 4.4 (statement).
- Ingalls, Logan, Patashnick, *Explicit coverings of families of elliptic surfaces by squares of curves*, https://arxiv.org/abs/2009.07807 : abstract, §1 opening, Remark 2.25.
- Varesco, *The Hodge conjecture for powers of K3 surfaces of Picard number 16*, https://arxiv.org/abs/2203.09778 ; *Hodge similarities, algebraic classes, and Kuga-Satake varieties*, https://arxiv.org/abs/2304.02519 : abstracts.
- Schlickewei, *The Hodge conjecture for self-products of certain K3 surfaces*, https://arxiv.org/abs/0907.2503 : abstract, introduction.
- Madapusi Pera, *The Tate conjecture for K3 surfaces in odd characteristic*, https://arxiv.org/abs/1301.6326 : abstract, introduction to Theorem 3 and the paragraph on integral models.
- Rizov, *Kuga-Satake abelian varieties in positive characteristic*, https://arxiv.org/abs/math/0608497 : abstract, introduction opening.
- Yang, *Isogenies between K3 surfaces over F̄_p*, https://arxiv.org/abs/1810.08546 : abstract, §1 to Theorem 1.4 (statement begun).
- Bragg, Yang, *Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*, https://arxiv.org/abs/2102.01193 : the paragraph of the introduction on Hecke orbits.
- Taelman, *Ordinary K3 surfaces over a finite field*, https://arxiv.org/abs/1711.09225 : introduction, Theorem 4.6 and its proof.
- Ito, *Unconditional construction of K3 surfaces over finite fields with given L-function in large characteristic*, https://arxiv.org/abs/1612.05382 : abstract, Theorem 1.3, Remark 1.4.
- Taelman, *K3 surfaces over finite fields with given L-function*, https://arxiv.org/abs/1507.08547 : abstract, statements of Theorems 1 to 3.
- Ito, Ito, Koshikawa, *CM liftings of K3 surfaces over finite fields and their applications to the Tate conjecture*, https://arxiv.org/abs/1809.09604 : abstract.
- Kedlaya, Sutherland, *A census of zeta functions of quartic K3 surfaces over F₂*, https://arxiv.org/abs/1511.06945 : abstract.
- Costa, Harvey, Kedlaya, *Zeta functions of nondegenerate hypersurfaces in toric varieties via controlled reduction in p-adic cohomology*, https://arxiv.org/abs/1806.00368 : abstract, introduction, §6 to Example 6.2.
- Harvey, *Computing zeta functions of arithmetic schemes*, https://arxiv.org/abs/1402.3439 : abstract, Theorem 1.1.
- Lubicz, Robert, *Computing isogenies between abelian varieties*, https://arxiv.org/abs/1001.2016 : abstract.
- Costa, Elsenhans, Jahnel, Voight, *Explicit modularity of K3 surfaces with complex multiplication of large degree*, https://arxiv.org/abs/2502.15052 : abstract and first page; not used.
- Opened, abstract only, not used: van Geemen, *Real multiplication on K3 surfaces and Kuga Satake varieties*, https://arxiv.org/abs/math/0609839 ; van Geemen, *Fourfolds of Weil type and the spinor map*, https://arxiv.org/abs/2205.00483 ; Voisin, *A generalization of the Kuga-Satake construction*, https://arxiv.org/abs/math/0503122 ; Floccari, Varesco, *Algebraic cycles on hyper-Kähler varieties of generalized Kummer type*, https://arxiv.org/abs/2308.04865 .

Not read; where they appear above it is through the papers that quote them, or [M]:

- Paranjape, *Abelian varieties associated to certain K3 surfaces* (read for note 14, not again here); Lombardo, Peters, Schütt (read for note 13).
- Voisin, *Remarks on zero-cycles of self-products of varieties* (1996); Morrison on Picard number at least 17; Kondo on curves of genus 3 and 4; Markman; O'Grady; Matsumoto, and Matsumoto–Sasaki–Yoshida; Artebani–Sarti; Laterveer; Garbagnati–Penegini.
- Deligne, *La conjecture de Weil pour les surfaces K3*; André; Kisin, Lam, Shankar, Srinivasan, *Finiteness of reductions of Hecke orbits*.
- Căldăraru, and Donagi–Pantev, on genus-one fibrations and twisted sheaves; Huybrechts–Stellari; Recillas; the classical formulas for Jacobians of genus-one curves.
- Kuperberg; Castryck et al. (CSIDH); Lauder–Wan and the deformation method; the costs of theta-function isogenies beyond Lubicz–Robert's abstract.
- The bodies of all papers above beyond the passages named.

In this repository: notes [2](02-six-line-count.md), [7](07-ring-changing-steps.md), [6](06-torsors-under-cm-class-groups.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [16](16-fourfold-over-f25-and-f49.md), [22](22-two-neighbours-of-elliptic-k3.md); the starting hypotheses are stated in [../SUMMARY.md](../../SUMMARY.md) §1.1.
