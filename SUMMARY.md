# Derived isogenies of ordinary K3 surfaces are not a new hard problem: the six-line family reduces to CSIDH on abelian fourfolds

*Christopher Dare, November 2026.*

**Evidence labels**, as used in the notes in [research/docs/](research/docs/README.md): **[P]** read in the cited source; **[M]** recalled, not checked against a source; **[D]** derived in this work, with the argument given in the cited note; **[E]** measured by the scripts in this repository (`research/experiments/`); **[H]** heuristic; **[R]** re-derived in a separate re-derivation check, with its own code ([note 18](research/docs/18-checks-of-notes-13-to-16.md) checks notes 13 to 16, [note 25](research/docs/25-checks-of-notes-20-23-24.md) notes 20, 23 and 24, and notes 26 to 29 state the findings of their own check, whose code is not part of this repository). A label in this summary is the label of the cited note.

---

## Abstract

Ordinary K3 surfaces over a finite field with a fixed Frobenius fall into torsors under finite abelian class groups of CM fields, and the group acts through twisted derived equivalences ("derived isogenies"). This looks like the setting of CSIDH, one dimension higher, and raises the question of whether it gives a new hard problem for cryptography. In the one family where everything can be computed, double planes branched along six lines (geometric Picard number 16), it does not. The surface problem becomes a class-group action on ordinary abelian fourfolds with multiplication by Z[i], made explicit over F_p through Paranjape's construction. The group has about 0.01·q² elements, so a 256-bit group needs q ≈ 2¹³¹ and about 525-bit instances, against 512 bits for CSIDH-512, with steps that cost about ℓ⁴ rather than ℓ. The honest party has no way to compute steps at useful sizes; the only working tool, the CM method, needs a class polynomial of degree |G|. The decision problem leaks one bit for 80 of the 477 examined polynomials that have a quadratic character: a quadratic character χ₋₁ is read off a sign δ computed on the fourfold. This carries over the genus-character attack of Castryck, Sotáková and Vercauteren, which affects class-group actions with 2-torsion (not CSIDH itself). A leakage audit of about 11,400 further candidate invariants found nothing else, and a test of characters of order 3 (as characters on two torsors, by separation of their cosets on seven) found none either. Below Picard number 16 no reduction is known, but neither is any computable step, so there is nothing to build on. All surface data have q ≤ 49, the torsors in the data have at most 36 surfaces, and the action was evaluated on torsors of at most 16. At cryptographic size the verdict rests on argument, not on measurement.

---

## 1. The question

### 1.1 Why K3 surfaces look promising

Isogeny-based cryptography uses a group acting on a set of geometric objects. In CSIDH [M: Castryck, Lange, Martindale, Panny, Renes, 2018] the set is the supersingular elliptic curves over F_p with a given endomorphism ring, the group is the class group of an imaginary quadratic order, and a group element of small norm ℓ acts by an isogeny of degree ℓ that can be computed in time polynomial in ℓ. The hard problem is *vectorisation*: given two curves, find the group element that takes one to the other. The best quantum attack is Kuperberg's subexponential hidden-shift algorithm [M]; the best classical attack is a meet-in-the-middle search in about √|G| steps [M].

Ordinary K3 surfaces over F_q have the same shape on paper. Taelman's theorem [P: Taelman, *Ordinary K3 surfaces over a finite field*] describes them by lattice data: a K3 lattice M, an endomorphism F (Frobenius) and an ample cone. On the transcendental part T of M, F generates a CM field E, and T becomes a Hermitian lattice over an order of E. Lattices in one genus are permuted by ideal classes, as elliptic curves are in CSIDH. Over the algebraic closure the corresponding surfaces are related by twisted derived equivalences, which Bragg and Yang showed exist in positive characteristic [P: Bragg, Yang, *Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*].

Two features made this worth a careful look:

- polarised K3 surfaces have a 19-dimensional moduli space, so the group could be large compared with the description of an instance;
- derived isogenies have no obvious analogue of the Vélu formulas, so a step might be easy for someone holding extra structure and hard for everyone else.

**The starting hypotheses.** Four hypotheses, formed by analogy with CSIDH before any computation, set the questions of this work. The notes refer to them as "the starting hypotheses":

- (a) twisted derived isogenies of ordinary K3 surfaces form a CSIDH-type class-group action that honest parties cannot run [H];
- (b) for geometric Picard number ρ̄ ≥ 17 (in practice ρ̄ = 18 and 20) the surfaces and their action reduce to elliptic curves and abelian surfaces [M];
- (c) by analogy with elliptic curves, the class of a surface in its torsor cannot be read off the surface: the genus is local data, and the class depends on the embedding and on the canonical lift [H];
- (d) explicit geometry for a step is known only at ℓ = 2 and 3, and naming a general double plane canonically is hard [M].

The notes in [research/docs/](research/docs/README.md) make (a) precise in one family where every object can be computed, test (c) at small sizes, and use (b) and (d) as background ([note 13](research/docs/13-hardness-on-paper.md), [note 17](research/docs/17-below-picard-16.md)).

### 1.2 The problems

Fix q = p^f, a family of K3 surfaces, and a Frobenius polynomial P. Let S_P be the set of isomorphism classes of ordinary surfaces in the family over F_q whose Frobenius acts on the transcendental part with characteristic polynomial P. Following [note 13](research/docs/13-hardness-on-paper.md) §2:

- **Vectorisation.** Given X, Y in one torsor of S_P, find g with g·X = Y, as a product of classes of small primes.
- **Parallelisation.** Given X, a·X, b·X, find ab·X.
- **Decision.** Given X, a·X, b·X and Z, decide whether Z = ab·X.

A protocol, and every known attack, needs three things: to name a surface canonically, to know P, and to apply a group element. The third is the crux. The question is therefore twofold. Is the problem new, that is, not reducible to a known group action? And is it usable, that is, can an honest party compute steps while an adversary cannot invert them?

---

## 2. Ordinary K3 surfaces as torsors

[Note 1](research/docs/01-prior-art.md) found that the basic statement is an assembly of published results [P]: Taelman's proof already gives the CM structure; Valloni's K3 class group [P: Valloni, *Complex multiplication and Brauer groups of K3 surfaces*] and Bayer-Fluckiger's torsor theorem [P: Bayer-Fluckiger, *K3 surfaces with maximal complex multiplication*] give the group over C; and Bragg–Yang's existence theorem gives the realisation by derived isogenies. The note, [note 6](research/docs/06-torsors-under-cm-class-groups.md), writes it out for finite fields:

- **Proposition A [D, from Taelman's theorem and class field theory].** For p ≥ 5 and geometric Picard rank ρ̄ ≥ 12, marked ordinary K3 surfaces with a fixed Néron–Severi part and a fixed rational Hermitian form on T fall into genera. Each genus is a principal homogeneous space under a finite abelian group G_T, a quotient of the norm-one idèles of E. Every element of G_T is the class of 𝔩/𝔩̄ for infinitely many primes 𝔩 of degree one. Non-maximal orders are included.
- **Remark B [D].** Taelman's conditions (M1)–(M5) do not fix which real place of E carries the positive plane; the essential image is the triples positive at the place given by the chosen embedding.
- **Propositions E and F [D, from Bragg–Yang's construction and the proof of their Torelli theorem; unchecked].** Over F̄_q any two surfaces of a genus are related by one twisted derived equivalence, with Brauer classes of order ℓ, for infinitely many ℓ. The neighbour at a prime 𝔩 above ℓ is a moduli space of twisted sheaves of rank ℓ.
- **Proposition G [D].** The partner is defined over F_q. The equivalence is defined over F_{q^d}, d the order of λ/q mod ℓ, and over F_q exactly when 𝔩 divides F − q, which forces ℓ² to divide |Br(X)|. The Artin–Tate formula gives a perfect square for |Br(X)| in all 854 cases checked (note 6 §6.6) [E].
- **Proposition H and notes 7–12 [D, E on lattices].** All steps at an odd prime ℓ ≠ p are Kneser neighbours, the isotropic lines of T/ℓT. Steps that change the endomorphism ring form ℓ-volcanoes and related graphs, worked out case by case ([7](research/docs/07-ring-changing-steps.md), [8](research/docs/08-below-the-top.md), [9](research/docs/09-ramified-over-e0.md), [10](research/docs/10-fourth-power.md), [11](research/docs/11-rank-six-fourth-power.md), [12](research/docs/12-rank-six-galois-action.md)). These parts are checked on lattices only, on polynomials not all known to come from surfaces, and they are not needed for the verdict.

So the group action exists in every family with ρ̄ ≥ 12 and is abelian. What remains is to compute anything in it.

---

## 3. The six-line family

A six-line surface is the double cover of the projective plane branched along six lines in general position, w² = l₁l₂⋯l₆. Over F̄_q it has Picard number at least 16 (the 15 exceptional curves over the intersection points, and a line class), so T has rank 6, E has degree 6, and its totally real subfield E₀ is cubic. The moduli space has dimension 4. Six lines up to projective equivalence and relabelling are a *configuration*; a configuration and its *associate* (a classical involution of six points in the dual plane [P: Dolgachev, Ortland]) give the same surface, so surfaces are pairs of configurations (note 20 §2; note 27 §1 item 3 explains why the two counts must be kept apart).

What makes the family tractable [D, note 2 and note 6 §5]: T = L(2) with L unimodular of rank 6, the orders and the gluing to Néron–Severi are explicit, and surfaces can be enumerated by brute force for small q. The first four comparisons between enumerated surfaces and enumerated lattices found no disagreement [E]:

| Comparison | Notes | Fields | Polynomials | Agree |
|---|---|---|---|---|
| Number of surfaces with all lines rational, against Hermitian lattices | [2](research/docs/02-six-line-count.md), [5](research/docs/05-non-split-count.md) | p = 7 to 37 (3,989 classes) | 44 | 44 |
| Frobenius conjugation against the step at p | [3](research/docs/03-torsor-action.md), 5 | F₂₅, F₄₉ | 32 | 32 |
| Frobenius action on NS of surfaces fixed by conjugation | [4](research/docs/04-forms-over-the-prime-field.md), 5 | F₂₅, F₄₉ | 9 (29 surfaces) | 29 |
| Forms whose lines are permuted by an involution, by type | [5](research/docs/05-non-split-count.md) | p = 5 to 23 (51,804 classes) | 169 | 169 |

These tests check the classification. None of them tests a step at a prime ℓ ≠ p on surfaces: no such step can be written down on the K3 side (note 3 §4, note 13 §2).

---

## 4. The reduction to abelian fourfolds

### 4.1 Over the complex numbers

Six-line surfaces are classically related to abelian fourfolds of Weil type with multiplication by Q(i). Paranjape [P: Paranjape, *Abelian varieties associated to certain K3 surfaces*, Compositio Math. 68 (1988)] constructs the fourfold as a Prym variety, and Lombardo, Peters and Schütt [P: *Abelian fourfolds of Weil type and certain K3 double planes*] show that the two moduli spaces coincide over C. Matsumoto's theta functions on the domain of U(2, 2) [P via Matsumoto–Terasoma and Clingher–Malmendier–Shaska] give the period map explicitly.

### 4.2 Proposition Y: polynomials

[Note 13](research/docs/13-hardness-on-paper.md) §4 [D, re-derived in note 18 [R]]. The Frobenius polynomials of the enumerated surfaces are exactly the *admissible* polynomials with P(q) = 64q⁵μ² and P(−q) = 64q⁵ν², where μ² and ν² are the orders of the Brauer groups of the surface and of its quadratic twist [E for seven primes; later nine, two of them out of sample, note 27 H12]. P is the second exterior power over Z[i] of the quartic

  χ = X⁴ − 4(ν + μi)X³ + e₁X² − 4q(ν − μi)X + q²,

and χ·χ̄ is the Weil polynomial of an ordinary abelian fourfold with multiplication by Z[i]. The identity holds for all q; note 27 checked it by integer arithmetic on all 1,606 stored polynomials up to q = 100,003 [E]. It is one of the few statements here that reaches cryptographic size, and it says nothing about surfaces by itself: it matches isogeny classes, not objects.

### 4.3 The fourfold built over F_p

[Note 14](research/docs/14-fourfold-of-a-six-line-surface.md) makes Paranjape's construction explicit and rational [D, E]. From six lines and an ordered pair of them it builds an elliptic curve E (a double cover of the cross-ratio line), four points and a class of degree 2 on E, and a curve C of genus 5 that is a cyclic cover of E of degree 4. The fourfold is the Prym variety of C → E. The characteristic polynomial of Frobenius on the i-eigenspace of H¹(C), computed from point counts of C over F_p, F_{p²}, F_{p³}, has second exterior power P or P(−t), the polynomial of the surface or of its twist:

- 534 configurations with irreducible P (278 surfaces), 4,983 ordered pairs of lines, p ≤ 37: no exception [E];
- 12,873 more pairs on 1,329 configurations with reducible P: no exception [E];
- over F₂₅ and F₄₉ ([note 16](research/docs/16-fourfold-over-f25-and-f49.md)): 37 and 1,108 configurations (21 and 564 surfaces), 858 and 3,324 pairs, every ordinary polynomial of note 3: no exception [E].

Two choices are fixed by the surface: the twist of E, and which half of a point of order 2 is used; a control with the other half gives, in 1,900 of 1,900 pairs, the polynomial of the surface with one line replaced by its conjugate rather than that of the surface ([note 18](research/docs/18-checks-of-notes-13-to-16.md), finding 8) [E]. Surfaces of sign − (Frobenius anticommuting with the action of i) and surfaces with reducible P need a twist by an extra involution of the curve, worked out in [note 21](research/docs/21-fourfold-forms-over-fp.md) [D, E, unchecked]; the main case, sign + with irreducible P, contains no collinear triple and is unaffected.

### 4.4 Does the map carry the group action?

This is the step the reduction needs, and it is the weakest link.

- **On class groups [E].** [Note 15](research/docs/15-fourfold-map-and-the-group-action.md) finds, for 278 polynomials, a homomorphism from the class group of the degree-8 CM field K of the fourfold to the group G of the surface that matches the predicted images of 77,946 primes. Note 18 corrected it: the map checked is the square of the map of tori, so the composites are fourth powers; the orders and the isomorphism of the odd parts stand [R].
- **On lattices over C [E].** The CM method of note 20 (§6 below) acts on fourfold lattices, and the ideal 𝔓_j𝔓_k·conj(𝔓_l𝔓_m) of K acts on the surfaces as the class in G of the prime of E at a_j a_k: 378 ideals, no conflict.
- **On the Prym variety itself [E, indirect].** An invariant δ of the Prym variety that is not an isogeny invariant splits 60 surfaces exactly as the CM action does (§7). If the map were not equivariant on those torsors, that agreement would have no reason.
- **Not done.** Equivariance is a sketch with five gaps (note 15 §2; note 27 H29). No step at ℓ ≠ p was ever computed on a Prym variety: an isogeny with kernel (Z/ℓ)⁴ of the Prym variety of a genus-5 curve needs theta coordinates or canonical lifts, and note 15 §5.1 found no feasible route.

**Conclusion of §4.** If the map is equivariant over F_q, which the evidence supports and nobody has proved, the vectorisation problem for six-line K3 surfaces is the vectorisation problem for a class-group action on ordinary abelian fourfolds with multiplication by Z[i]: CSIDH-type, and not a new assumption.

---

## 5. The group and its size

For a polynomial whose order is maximal, the group acting on the surfaces is

  G = Cl(E) / ⟨Cl(E₀), ramified primes⟩

(note 13 Lemma 9 [D, re-derived in note 18]). Its size was computed (class groups under GRH, with PARI) for all admissible polynomials up to q = 101 and for samples up to q = 100,003: 512 groups in all [E]. Selected rows of note 13 §3.3:

| q | groups | median \|G\| | median \|G\|/q² |
|---|---|---|---|
| 37 | 7 | 10 | 0.0073 |
| 101 | 25 | 100 | 0.0098 |
| 1,009 | 40 | 11,658 | 0.0115 |
| 10,007 | 30 | 872,620 | 0.0087 |
| 100,003 | 12 | 83,306,780 | 0.0083 |

Over the 512 groups |G|/q² has median 0.0092 and quartiles 0.0067 and 0.0128 (recounted in note 27; the values first quoted in note 13, 0.0095, 0.0068 and 0.0129, are those of the 486 groups with q ≥ 41 used for the first fit), with range 0.0003 to 0.024. A least-squares fit of log|G| against log q gives exponent 1.98 over q ≥ 211 with constant 0.011 [E]. The class-number heuristic gives |G| of order q^{(20 − ρ̄)/2} for transcendental rank 22 − ρ̄, the square root of the number of points of the moduli space, as for elliptic curves [H].

Consequences [H]:

- A group of 2²⁵⁶ elements needs q ≈ 2¹³¹. A surface is named by four elements of F_q (its four moduli, in a canonical form that note 13 §2 gives [D]), about 525 bits. CSIDH-512 has a group of the same size with 512-bit instances. The law |G| ≈ 0.01·q² is extrapolated 2¹¹⁴ times beyond the data (note 27 H10).
- Structure [E]: |G| is even in 509 of 512; the odd part is cyclic in 506; exactly one prime of E is ramified over E₀ in every field, above 2. No Pohlig–Hellman reduction is known for a group action without computable invariants of subgroup orbits [H], so the 2-part matters only through characters (§7).

---

## 6. Computability: the honest party cannot run it

### 6.1 What would be needed

A step at ℓ on the K3 side is a moduli space of twisted sheaves of rank ℓ (Proposition F), with no equations (note 6 §7). On the fourfold side it is an isogeny with kernel (Z/ℓ)⁴ commuting with i, costing of the order of ℓ⁴ field operations once the kernel is found [M: Lubicz, Robert, *Computing isogenies between abelian varieties*], against of the order of ℓ for an elliptic curve. The zeta function of a K3 surface at a prime of cryptographic size costs at least about √p [M: Harvey; Costa, Harvey, Kedlaya], so instances would have small characteristic and q = p^f. Given the action at unit cost, the attacks are CSIDH's [M]: Kuperberg's algorithm, 2^{O(√log|G|)} evaluations, because G is abelian; and about √|G| ≈ 0.1·q classical evaluations.

### 6.2 The CM method works, on small groups

[Note 20](research/docs/20-cm-method-and-the-action.md) evaluates the action through complex multiplication [E; re-checked in note 25]:

1. Enumerate the principally polarised Hermitian lattices of the fourfold side. They are not those of the maximal order of K: all 234 such lattices in 21 fields are odd, and a six-line surface needs an even one. The right lattices are modules over R = Z[i, π, q/π, (π − 1)/(1 + i)], the last generator expressing that the six lines are rational; their number equals the number of surfaces in every field treated.
2. Compute their period points, the theta constants on the domain of U(2, 2) at 300 to 3,000 digits, and the invariants J₂, …, J₆ of six lines (formulas from Matsumoto–Terasoma and Clingher–Malmendier–Shaska [P]; Matsumoto's own paper not read).
3. Recognise the class polynomial over E and reduce it modulo the prime of E above p at which Frobenius is a unit.

Results:

- **Points.** The roots of the class polynomial are exactly the enumerated surfaces: 74 of 74 in 15 fields (note 20), 98 of 98 in 17 fields with the two torsors of 12 surfaces added in [note 24](research/docs/24-the-form-mod-4-is-the-character.md) [E]. With notes 26 and [29](research/docs/29-characters-of-order-3.md), 150 surfaces in 21 fields (note 27 H7).
- **Action.** The permutations induced by ideals depend only on their class in G, compose, are free, and are transitive in 14 of 15 fields (in the fifteenth, primes ℓ ≡ 1 mod 4 generate the kernel of χ₋₁, as note 19 predicted) [E].
- **Steps are Frobenius.** Raising to the ℓ-th power on the roots of the class polynomial modulo a prime above ℓ equals the permutation of the ℓ-step in 346 of 346 cases; with the inverse it fails in 136 of 136 [E]. One sign, the orientation of the period point, was fixed on one field and then held on the others.
- **Labels.** 67 surfaces over F_p are labelled by elements of G.

### 6.3 Why this does not help an honest party

The CM method needs a class polynomial of degree |G|. The first precisions at which the polynomials were recognised (note 27 H17) [E]:

| \|G\| | 2 | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| digits | 300 | 300–700 | 700 | 1,200 | 1,200–2,000 | 3,000 |

These are the first precisions tried that worked, so they are upper bounds. The first runs at q = 43 stalled in a diagnostic factorisation, after the class polynomial had been recognised, and bear on neither cost nor precision. On the reruns, |G| = 12 needed 3,000 digits there, and |G| = 24 needed 4,500 digits, for a coordinate chosen to be injective modulo 43 ([note 29](research/docs/29-characters-of-order-3.md) §4.2, §7.4). At |G| ≈ 2²⁵⁶ it cannot run. That is an argument, not a measurement. It is the same obstacle as for ordinary elliptic curves (the Couveignes–Rostovtsev–Stolbunov setting [M]), without the Vélu formulas that make CSIDH's steps cheap. No other way to evaluate a step at a prime ℓ ≠ p was found on either side.

---

## 7. The decision problem: one bit leaks

### 7.1 Which characters could be read

Genus theory gives the attack of Castryck, Sotáková and Vercauteren [P: *Breaking the decisional Diffie–Hellman problem for class group actions using genus theory*, CRYPTO 2020]: a character of order 2 of the class group that can be evaluated on single objects decides DDH with advantage 1/2. [Note 19](research/docs/19-quadratic-characters-from-a-surface.md) [D, E] showed that every character of order 2 of G is that of an unramified extension E(√d)/E with d in E₀, and that the CSV pairings, which separate genera, have nothing to read on the K3 lattice: the surfaces of one polynomial form a single genus, and G is the analogue of Cl(O)², so its characters of order 2 are analogues of characters of order 4. None of the natural invariants tested on 419 surfaces computes a character. One character stood out as within reach: that of d = −1, χ₋₁, the degree of a step modulo 4. It exists for 80 of the 477 examined polynomials that have a quadratic character ([note 19](research/docs/19-quadratic-characters-from-a-surface.md)), exactly when the real quartic field of the fourfold has a prime above 2 of odd degree, and only for q ≡ 1 mod 4 [E].

### 7.2 The sign δ is the character χ₋₁

[Note 23](research/docs/23-jacobian-arithmetic-levels-and-the-form-mod-4.md) implements arithmetic in the Jacobian of the genus-5 curve C (Khuri-Makdisi's linear algebra on spaces of sections), validated on 35 curves [E], and defines an invariant δ = ±1 from the Hermitian form modulo 4 at the prime of degree 1 above 2, read off a Tate pairing of order 16 over F_{p⁴} [D]. On the three torsors where χ₋₁ exists at p = 29 and 37, δ is independent of the pair of lines (32 surfaces, 106 curves) and balanced.

[Note 24](research/docs/24-the-form-mod-4-is-the-character.md) joined δ to the CM action [E]. The steps at primes ℓ ≡ 1 mod 4 generate a subgroup of index 2 with two orbits, and δ is constant on each orbit and opposite between them: 32 surfaces of 32 (G = Z/12, Z/12, Z/4 × Z/2). A balanced function would agree by chance with probability 1 in 7,470,540. The separating function is i times a function over E, so the extension is E(i) and the character is χ₋₁.

[Note 26](research/docs/26-delta-is-the-character-at-q-41.md) tested this out of sample at q = 41, the first prime ≡ 1 mod 8 in the tests, with the definition of δ written down before any curve was computed [E]: 16 of 16 surfaces (G = Z/16) and 12 of 12 (G = Z/12). In all, **60 of 60 surfaces in five torsors at three primes**; chance probability 4.5·10⁻¹⁴ for a balanced function (note 27 H24). Three routes to the cosets agree: orbits of ideals, the factors of the class polynomial over E(i), and the class of the polarisation modulo 4 on each lattice (note 25). They check the code and the recognition rather than giving three separate determinations: an indicator recognised as i times a polynomial over E must give the split over E(i) (note 26 §1). A separate check of note 26 recomputed δ from all 783 records at q = 41 and rebuilt the class polynomial for e₁ = 22 with its own recognition, and its own cubic and quartic class fields give the same split; for that field the split does not depend on GRH [R].

**Consequence [D, from the above; H for the cost].** For a polynomial with χ₋₁, the bit χ₋₁(g) of the class g with Y = g·X is δ(X)·δ(Y), computable from the two surfaces through their fourfolds. For X, aX, bX and a candidate Z = cX, the relation c = ab forces δ(X)δ(cX) = δ(aX)δ(bX), which a random c satisfies half of the time. That is the CSV attack on the decision problem, carried to the fourfold. The countermeasures are theirs: steps only in the kernel, or polynomials without the character. It touches 80 of the 477 examined polynomials that have a quadratic character (all with q ≡ 1 mod 4), and says nothing about the search problem.

**On the reference.** The attack of Castryck, Sotáková and Vercauteren affects class-group actions whose class group has even order. Their paper says it does not affect CSIDH, whose class group (p ≡ 3 mod 4) has no 2-torsion [P; [note 25](research/docs/25-checks-of-notes-20-23-24.md), finding 5]. So the analogue here is the attack on group actions with 2-torsion, not an attack on CSIDH.

### 7.3 The leakage audit: nothing else

[Note 28](research/docs/28-leakage-audit.md) [E; every p-value reproduced by an separate re-implementation, R] joined every stored surface invariant to the CM labels, with the candidates, statistic (mutual information with the label difference), exact permutation null and Holm correction fixed before the join:

- **Full labels**, 9 torsors (56 surfaces, |G| ≤ 8, none with χ₋₁): 4,838 candidates (residue symbols of order 2 and 4 of J-monomials and power sums, note 19's invariants, the CM coordinate, curve data, relational symbols of two surfaces). Smallest p 0.0038; Holm 1 for all.
- **Cosets of ker χ₋₁**, 5 torsors (60 surfaces): 6,579 candidates. δ comes first (p = 4.5·10⁻¹⁴, the positive control); the next is at p = 0.0010, Holm 1. Nothing on the K3 side reads χ₋₁.
- **Power.** A planted order-2 character on every torsor is found 200 of 200 times; on 4 of 9 torsors, 37 of 200; with values actually wrong one time in 8, 63 of 200. Characters of order 4 were tested only through their squares: a quartic symbol of a quantity over F_q exists only for q ≡ 1 mod 4, on 3 of the 9 torsors, and no candidate could reach the threshold with an order-4 reading. Most candidates could never reject at all (3,670 of 3,802 on the labelled torsors), so the correction is stricter than it needs to be. Shuffled labels give the same spread of p-values as the true ones.

About 11,400 candidates in all. The audit cannot see a leak outside its families, a leak on one or two torsors, characters of order 3 (§7.4), or anything at |G| > 8; and it cannot see the one leak that certainly exists, the label as a function of the complete invariants computed through the class polynomial.

### 7.4 Characters of order 3

The audit of note 28 could not see characters of order 3: only two of its labelled torsors have 3 ∣ |G|, both with q ≡ 2 mod 3. [Note 29](research/docs/29-characters-of-order-3.md) [E; re-checked, see §9] tests them, with the tests fixed before the join (one early join on three torsors, note 29 §4.3):

- **The cosets.** On the seven torsors with 3 ∣ |G| at q = 23 to 41 (72 surfaces; G = Z/6 twice and Z/12 five times), the class polynomial factors into three pieces over the cubic class field L of ker χ₃, and the roots of the pieces modulo 𝔭 are the three cosets of the kernel. Two of the class polynomials are new, and in both the roots are the enumerated surfaces. The cosets agree with the labels modulo 3 on the two Z/6 torsors, and meet the cosets of ker χ₋₁ in a regular grid on the three Z/12 torsors that have both. At q = 43 no class polynomial was obtained in the run itself (§6.3; note 29 §4.2). Reruns after it gave the cosets on all three q = 43 torsors (Z/24, Z/12, Z/6 × Z/2; 48 surfaces), the Z/24 one at 4,500 digits with a coordinate chosen to be injective modulo 43 (note 29 §7.4).
- **On paper [D].** For q ≡ 2 mod 3 every element of F_q is a cube, also in F_{q²}, so no cubic residue symbol of a quantity defined over F_q, or of its square root, can read χ₃. For q ≡ 1 mod 3 a cubic symbol can, but at q = 37 the cube roots of unity are themselves cubes, so a quantity R with R³ = c³·d, the natural reading, would have a constant symbol there.
- **Characters, at q = 37** (two torsors of 12): none of 184 cubic residue symbols (of J-monomials, of the power sums and the conic, of natural quantities of the lines, of the CM coordinate) is a character of order 3 on either torsor, and none is a function of the coset difference beyond chance (best p 0.030, Holm 1). Among 2,839 candidates the smallest p is 0.0046, Holm 1: the exact indicator of one coset on one torsor, which occurs at chance level.
- **Separation, at q ≡ 2 mod 3** (five torsors, 48 surfaces): no symbol of order 2 or 4, value or relational symbol separates the cosets (smallest p 1.4·10⁻⁴ among 2,448, Holm 0.34, not replicated on the held-out torsors at q = 37).
- **Characters, at q = 43** (three torsors, 48 surfaces, scored after the run): none of the 184 cubic symbols is a character of order 3 on any of them (best Holm 0.28), at a prime where the cube roots of unity are not cubes. With them, population A has smallest p 3.3·10⁻⁵, Holm 0.095. A₄₃ alone has one candidate past Holm, the quadratic symbol of R2(X) + R2(Y) (Holm 0.0057). It reverses direction at q = 37 and does not replicate on the held-out Z/24 torsor (p 0.22) (note 29 §7.4).
- **The fourfold.** On the three torsors with fourfold data (q = 29, 37 and 41, one each), no quantity of the fourfold per pair of lines separates the cosets either (best p 0.0061, Holm 0.080).
- **Power.** The positive control, δ against the cosets of ker χ₋₁ on those three torsors, gives p = (1/462)³ = 1.0·10⁻⁸, the predicted value. An exact order-3 character on both torsors at q = 37 is found 200 times in 200, and with one value in twelve wrong 173 times. A reading on one torsor only, the indicator of one coset at q = 37, or a weak bias would not be found.

So the negative result of §7.3 extends to characters of order 3, within stated families: as characters on the two torsors at q = 37, and on all seven as a separation of the cosets. A reading of χ₃ at q ≡ 2 mod 3 would need a quantity outside F_q or a pairing of order 3 on the Prym variety, and neither was tried.

---

## 8. Below Picard number 16

Could the action be a hardness assumption that does *not* pass through abelian varieties? [Note 17](research/docs/17-below-picard-16.md) [P, M, D; a reading of the literature, nothing computed] says no, for three reasons:

1. **Nobody can evaluate the action.** No step at any prime has been computed on a K3 surface with ρ̄ < 16, and no formula exists. The adversary and the honest party lack the same operation, which is not hardness.
2. **The reduction exists in principle.** The Kuga–Satake abelian variety exists over finite fields for every ρ̄, as a morphism of moduli spaces over Z[1/2] [P: Madapusi Pera]. Only formulas are missing. Its simple factor over F_q has dimension 2^{(20 − ρ̄)/2}: 8 at ρ̄ = 14, 16 at ρ̄ = 12, then 32 to 512, far beyond theta-function methods. At ρ̄ = 14 the abelian side is 8-folds with quaternionic multiplication, known on periods only [P: Poon].
3. **The ceiling is CSIDH's.** The group is abelian at every ρ̄, so Kuperberg applies given the action, and by the square-root count an instance has about 512 bits for a group of 2²⁵⁶ at every ρ̄ [H].

One lead was tested. On a double cover of P¹ × P¹ branched along a (4, 4) curve, the Jacobians X₁, X₂ of its two genus-one fibrations should be 2-neighbours [P over C: van Geemen]. [Note 22](research/docs/22-two-neighbours-of-elliptic-k3.md) built 21,240 such pairs over F₅ and F₇ at ρ̄ = 12, 10, 2 [E]: everything the lattice side predicts and that can be read off the pair holds, and random surfaces of the same families fail each check; under a condition forcing the maximal order (972 pairs) the class of 𝔩/𝔩̄ is never trivial. That the step is the predicted class, rather than another element, was not tested. It is half of a step at one prime, and |G| ≈ 0.08·p⁴ at ρ̄ = 12. The verdict stands: below 16 there is no known reduction, but also no computable action and so nothing to assume.

---

## 9. How the claims were checked

- **Three re-derivation checks.** [Note 18](research/docs/18-checks-of-notes-13-to-16.md) checks notes 13, 14, 15, 16: no error that changes the conclusion; Lemma 9 and Proposition Y re-derived; two errors in note 15 §3 corrected. [Note 25](research/docs/25-checks-of-notes-20-23-24.md) checks notes 20, 23, 24: the claim that δ is χ₋₁ does not fail; two further routes to the same split; the CSIDH reference corrected; two gaps closed or narrowed. A third check, with its own code (not part of this repository), re-derived notes 26 to 29 from the stored files, and each of those notes states its findings: δ = χ₋₁ at q = 41 reproduced from both sides, with no number wrong; all 170 comparisons of the ledger recounted, with five statements corrected; every p-value of the leakage audit reproduced, with its account of power corrected (order 4 was not testable, and the plants with errors had been misread); and the cosets and statistics of the order-3 test reproduced, with the error rates of its plants and its order of events corrected. No test result and no conclusion changed. None of the checks covers the GP code of the Jacobian group law or the completeness of the lattice enumeration of note 20.
- **A hypothesis ledger.** [Note 27](research/docs/27-hypothesis-ledger.md) lists 47 hypotheses in five groups, classed as proved on paper (4), supported by data (25), heuristic (9), open (8) or contradicted (1, plus superseded forms of two others), each with its falsifier. `research/experiments/sixlines_ledger.py` makes 170 comparisons between numbers as first quoted in the notes and numbers recomputed from the stored files: 162 agree. The 8 that differ (7 distinct numbers, one of them compared twice) come from three slips, all now corrected in the notes and none changing a conclusion: the median and quartiles of |G|/q² in note 13, one table row of note 13, and configurations counted as "surfaces" in notes 14 and 16.
- **Pre-registered tests.** Notes 26, 28 and 29 wrote their definitions, statistics and nulls before joining data to labels, and say so, with the order of events (note 26 §3, note 28 §3, note 29 §4.3). Note 29 records one early join, on three torsors, while its tests were being written.
- **In this repository.** On a fresh clone with the data archive unpacked, the ledger reproduces its stored log byte for byte, the self-test of the CM method passes and matches its stored log, `sixlines_prym.py` reproduces the stored records for p = 11 and 13, and both notebooks run ([data/README.md](data/README.md#verification)).

None of this is an outside review.

---

## 10. Limits

1. **Small instances only.** All surface data have q ≤ 49. The torsors in the data have at most 36 surfaces (at q = 49), those on which the action was evaluated at most 16, and the leakage audit's fully labelled torsors at most 8. Group sizes reach q = 100,003; a hard instance would need q ≈ 2¹³¹.
2. **Argument, not measurement, at cryptographic size.** Note 27 §1 names what carries the verdict there: Proposition Y (valid for all q); the equivariance of surface ↦ fourfold (a sketch with five gaps, never tested at a prime ℓ ≠ p on the Prym variety); the CM points being the canonical lifts (a sketch); the generic attacks and the costs of zeta functions and steps (recalled, [M]); and the extrapolation of |G| ≈ 0.01·q². No measurement here is at a size where the problem could be hard.
3. **One family.** The reduction is made explicit only for six-line double planes. Other families with ρ̄ = 16 were not treated over F_q.
4. **Mostly maximal orders.** G and its size are computed for maximal orders; non-maximal orders appear in the classification (notes 6–12) but not in the hardness analysis beyond counts.
5. **Sources partly recalled.** Several costs and the CSIDH background are [M]. Matsumoto's 1993 paper was not read; its formulas were taken from two later papers.
6. **No outside review.** Nothing here has been refereed; see [A note on method](#a-note-on-method).

---

## 11. Open questions

1. **Equivariance.** Prove that the surface ↦ Prym map of note 14 carries the class-group action on surfaces to that on fourfolds, over F_q, or compute one step at a prime ℓ ≠ p on the Prym variety (theta coordinates of a Prym variety of a genus-5 curve, or canonical lifts).
2. **Reading the class off a surface.** Is there any computation from a surface, beyond the class polynomial of degree |G|, that reveals information about its label? Notes 28 and 29 are negative results on groups of order at most 16, within stated families.
3. **χ₋₁ in general.** Prove that δ is χ₋₁ for every polynomial where the character exists, including q ≡ 1 mod 16 and the curves of exponent 64 not yet computed (note 26 §8).
4. **Characters of order 3.** Can χ₃ be read at all: from a quantity in F_{q²} \ F_q, by a pairing of order 3 on the Prym variety, or by a cubic symbol at a prime q ≡ 1 mod 3 where the cube roots of unity are not cubes (q = 43 is the next such prime in the data; on its three torsors, scored after the run of note 29, none of 184 cubic symbols is a character, note 29 §7.4)?
5. **Below 16.** Can the half-step at ℓ = 2 on elliptic K3 surfaces (note 22) be completed and closed into a cycle? Do the 8-folds at ρ̄ = 14 have a geometric construction?
6. **Non-maximal orders and ring-changing steps** as a source of structure for an adversary (the volcanoes of notes 7–12).
7. **A literature check** before any claim that something here is new: the classification is an assembly of known theorems (note 1), and the six-line/fourfold correspondence is classical over C.

---

## References

### Notes in this repository

All in [research/docs/](research/docs/README.md): [1](research/docs/01-prior-art.md) prior art; [2](research/docs/02-six-line-count.md), [3](research/docs/03-torsor-action.md), [4](research/docs/04-forms-over-the-prime-field.md), [5](research/docs/05-non-split-count.md) counts and the step at p; [6](research/docs/06-torsors-under-cm-class-groups.md) torsors under CM class groups; [7](research/docs/07-ring-changing-steps.md)–[12](research/docs/12-rank-six-galois-action.md) ring-changing steps; [13](research/docs/13-hardness-on-paper.md) hardness on paper; [14](research/docs/14-fourfold-of-a-six-line-surface.md), [15](research/docs/15-fourfold-map-and-the-group-action.md), [16](research/docs/16-fourfold-over-f25-and-f49.md), [21](research/docs/21-fourfold-forms-over-fp.md) the fourfold; [17](research/docs/17-below-picard-16.md), [22](research/docs/22-two-neighbours-of-elliptic-k3.md) below Picard 16; [18](research/docs/18-checks-of-notes-13-to-16.md), [25](research/docs/25-checks-of-notes-20-23-24.md) re-derivation checks; [19](research/docs/19-quadratic-characters-from-a-surface.md), [23](research/docs/23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](research/docs/24-the-form-mod-4-is-the-character.md), [26](research/docs/26-delta-is-the-character-at-q-41.md) characters and δ; [20](research/docs/20-cm-method-and-the-action.md) the CM method; [27](research/docs/27-hypothesis-ledger.md) the ledger; [28](research/docs/28-leakage-audit.md) the audit; [29](research/docs/29-characters-of-order-3.md) characters of order 3.

### Literature cited by the notes

Each note's *Sources* section says which passages were read [P] and which were recalled [M]. The main ones:

- Bayer-Fluckiger, *K3 surfaces with maximal complex multiplication*, arXiv:2408.11643.
- Bragg, Yang, *Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*, arXiv:2102.01193.
- Castryck, Lange, Martindale, Panny, Renes, *CSIDH: an efficient post-quantum commutative group action*, ASIACRYPT 2018 [M in the notes].
- Castryck, Sotáková, Vercauteren, *Breaking the decisional Diffie–Hellman problem for class group actions using genus theory*, CRYPTO 2020, IACR ePrint 2020/151.
- Clingher, Malmendier, Shaska, *Six line configurations and string dualities*, arXiv:1806.07460.
- Costa, Harvey, Kedlaya, *Zeta functions of nondegenerate hypersurfaces in toric varieties via controlled reduction in p-adic cohomology*, arXiv:1806.00368.
- Dolgachev, Ortland, *Point sets in projective spaces and theta functions*, Astérisque 165 (1988).
- van Geemen, *Kuga–Satake varieties and the Hodge conjecture*, arXiv:math/9903146; *Some remarks on Brauer groups of K3 surfaces*, arXiv:math/0408006.
- Harvey, *Computing zeta functions of arithmetic schemes*, arXiv:1402.3439.
- Huybrechts, Stellari, *Equivalences of twisted K3 surfaces*, arXiv:math/0409030.
- Ito, Ito, Koshikawa, *CM liftings of K3 surfaces over finite fields and their applications to the Tate conjecture*, arXiv:1809.09604.
- Kneser, *Klassenzahlen definiter quadratischer Formen*, Arch. Math. 8 (1957).
- Kohel, *Endomorphism rings of elliptic curves over finite fields*, thesis, Berkeley 1996.
- Kuperberg, *A subexponential-time quantum algorithm for the dihedral hidden subgroup problem*, SIAM J. Comput. 35 (2005) [M in the notes].
- Khuri-Makdisi, *Linear algebra algorithms for divisors on an algebraic curve*, Math. Comp. 73 (2004); *Asymptotically fast group operations on Jacobians of general curves*, Math. Comp. 76 (2007) [M in the notes].
- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, arXiv:1209.5997.
- Lubicz, Robert, *Computing isogenies between abelian varieties*, arXiv:1001.2016.
- Madapusi Pera, *The Tate conjecture for K3 surfaces in odd characteristic*, arXiv:1301.6326.
- Matsumoto, *Theta functions on the bounded symmetric domain of type I₂,₂ and the period map of a 4-parameter family of K3 surfaces*, Math. Ann. 295 (1993), 383–409 (not read; its formulas are taken from the two later papers below).
- Matsumoto, Terasoma, *Thomae type formula for K3 surfaces given by double covers of the projective plane branching along six lines*, arXiv:1001.4865.
- Paranjape, *Abelian varieties associated to certain K3 surfaces*, Compositio Math. 68 (1988), 11–22.
- Poon, *Kuga–Satake construction on families of K3 surfaces of Picard rank 14*, arXiv:2504.08270.
- Rizov, *Complex multiplication for K3 surfaces*, arXiv:math/0508018; *Kuga–Satake abelian varieties in positive characteristic*, arXiv:math/0608497.
- Taelman, *Ordinary K3 surfaces over a finite field*, arXiv:1711.09225; *K3 surfaces over finite fields with given L-function*, arXiv:1507.08547.
- Valloni, *Complex multiplication and Brauer groups of K3 surfaces*, arXiv:1804.08763; *Neighbors and arithmetic of isogenous K3 surfaces*, arXiv:2206.02560.
- Yang, *Isogenies between K3 surfaces over F̄_p*, arXiv:1810.08546.

The texts of these papers are not included in this repository.

---

## A note on method

The idea for this project and its direction are my own. AI tools (Anthropic's Claude) were used extensively in carrying it out, including the derivations, code, computations and drafting. The results have not been independently reviewed. If you would like more detail on how the work was done, get in touch.
