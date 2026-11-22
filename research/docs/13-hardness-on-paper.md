# 13. The hardness question on paper

*Asks, on paper, whether finding a path between two ordinary K3 surfaces with the same Frobenius is a hard problem, and whether it is a new one; starting hypothesis (a), stated in [../SUMMARY.md](../../SUMMARY.md) §1.1, expects a class-group action of CSIDH type, by analogy [H]. Uses [note 6](06-torsors-under-cm-class-groups.md) for the group action, and builds on notes 2, 3, 11 and 12. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_groupsize.py`; **[H]** heuristic. Checked in [note 18](18-checks-of-notes-13-to-16.md), whose corrections are stated here in corrected form; Propositions A to H of note 6, on which it rests, have unchecked steps (note 6 §8).*

---

## 1. Result

For the six-line family (ρ̄ = 16, transcendental rank 6), and only on paper:

1. **Size of the group [E, H].** For a polynomial whose order is maximal, the group acting on the surfaces is G = Cl(E)/⟨Cl(E₀), ramified primes⟩ (Lemma 9), and |G| ≈ 0.01·q² from q = 13 to q = 100,003. The heuristic for other Picard numbers is q^{(20 − ρ̄)/2}. A group of 2²⁵⁶ elements needs q near 2¹³¹.
2. **Which polynomials occur [E].** The Frobenius polynomials of the enumerated surfaces are exactly the admissible polynomials with P(q) = 64q⁵μ² and P(−q) = 64q⁵ν², for the seven primes with data (nine since, two of them out of sample: [note 27](27-hypothesis-ledger.md), H12). Here μ² is the order of the Brauer group and ν² that of the quadratic twist [M]. The order is maximal only when μ = ν = 1; in general an isogeny class is a stack of levels, and the ring-changing steps of notes 7 to 12 are the steps between them.
3. **The problem transfers to abelian fourfolds [P over C; D for polynomials].** P is the second exterior power over Z[i] of the quartic χ = X⁴ − 4(ν + μi)X³ + e₁X² − 4q(ν − μi)X + q², and χ·χ̄ is the Weil polynomial of an ordinary abelian fourfold with multiplication by Z[i] (Proposition Y). Over C the six-line surfaces and these fourfolds have the same moduli space (Lombardo–Peters–Schütt). A step of order ℓ on the K3 side should be one isogeny with kernel (Z/ℓ)⁴ commuting with i on the fourfold side.
4. **So at ρ̄ = 16 the problem is a CSIDH-type class-group action, as starting hypothesis (a) ([../SUMMARY.md](../../SUMMARY.md) §1.1) expects by analogy [H], and it reduces to abelian varieties, as hypothesis (b) states for ρ̄ ≥ 17 [M].** If the correspondence holds over F_q, the K3 problem in this family is the vectorisation problem for a class-group action on abelian fourfolds: no harder, and not a new assumption. The steps become computable in principle on the fourfold side, at a cost far above CSIDH's, for a group of the square root of the size one would expect from the instance length.
5. **What is not settled.** The correspondence is not proved over F_q. This note checks only its shadow on polynomials; [note 14](14-fourfold-of-a-six-line-surface.md) builds the fourfold over F_p from the six lines, as the Prym variety of a curve of genus 5 (Paranjape's construction, read and made explicit). The characteristic polynomial of Frobenius on its i-eigenspace is the quartic χ of item 3, and its second exterior power is the polynomial of the surface or of its twist, for all 534 enumerated configurations with irreducible P (278 surfaces, a configuration and its associate being one surface) and all 4,983 pairs of lines computed. So the correspondence is tested on point counts of surfaces and curves, and the map from a surface to its fourfold is explicit. Still open after note 14: a proof; that the map carries the group action over; the reverse map; any step (steps for groups of order at most 8 are computed by the CM method in [note 20](20-cm-method-and-the-action.md); see §9). Whether the class can be read off a surface, the "no easier" half, has no evidence either way in this note (see §7 for the test of note 28). This note runs nothing on surfaces. Picard numbers below 16, where no such transfer is known, are not touched here (see §9 and [note 17](17-below-picard-16.md)).

So the conclusion of this note on paper is that this is a CSIDH-type problem in dimension 4, not a new source of hardness.

---

## 2. The problems

Fix q = p^f, the six-line family, and a polynomial P. Let S_P be the set of isomorphism classes of ordinary split six-line surfaces over F_q with Frobenius polynomial P on the transcendental part.

- **Structure of S_P.** By Propositions A and C of note 6, S_P is a finite union of torsors. Each torsor is one genus, with one endomorphism ring; in the split case one surface corresponds to one lattice (L, F).
- **Vectorisation.** Given X and Y in one torsor, find g in the group with g·X = Y, as a product of classes of small primes.
- **Parallelisation.** Given X, a·X, b·X, find ab·X.
- **Decision.** Given X, a·X, b·X and Z, decide whether Z = ab·X.

What an algorithm, honest or not, needs:

1. *To name a surface.* Six lines up to projective equivalence, relabelling and association have four moduli; putting four lines in standard position and taking the least of the 1,440 images gives a canonical label [D]. This is easy here, unlike for general double planes, where naming a surface canonically is hard.
2. *To know P.* The zeta function of one surface.
3. *To apply a group element.* A step at a prime ℓ: on the K3 side this is the moduli space of twisted sheaves of Proposition F, with no equations (note 6 §7 item 3).

Without item 3 neither a protocol nor any of the known attacks can be run. That is not hardness: both sides lack the same operation.

---

## 3. The group

### 3.1 The maximal order

Notation of note 6 §2: E = Q[F] of degree 6, E₀ its real cubic subfield, T = L(2) with L unimodular, T¹ the torus of the x with x·x̄ = 1.

**Lemma 9.** [D] Suppose End(L) = O_E.

1. No odd prime of E₀ ramifies in E.
2. The isomorphism classes of lattices (L′, F) locally isometric to L form a torsor under
   > G = T¹(Q) \ T¹(A_f) / Stab(L̂) = Cl(E) / ⟨ image of Cl(E₀), classes of the primes of E ramified over E₀ ⟩.
3. The group G_T of Proposition A, which acts on marked surfaces, maps onto G with a kernel of order at most 1,440, the order of O(q_N).

*Proof.* (1) L is a fractional ideal I with α·I·Ī·𝔡_E = O_E, α in E₀. At a prime 𝔔 of E ramified over an odd prime of E₀, the valuations of α, of I·Ī and of 𝔡_{E₀} are even, and that of 𝔡_{E/E₀} is 1. (2) The map s ↦ s/s̄ from the finite idèles of E onto T¹(A_f) is surjective by Hilbert 90. Suppose s/s̄ = (e/ē)·u with e in E^× and u in the stabiliser, so u is a unit of norm 1 at every prime. At a prime inert or split over E₀, u = v/v̄ with v a unit. At a ramified prime every element of E_𝔔^× is sent to a unit. So s is in E^×·A_{E₀,f}^×·Ô_E^× times the E_𝔔^× at ramified 𝔔, and conversely. (3) The stabiliser of L̂ and the group U_T differ only at 2, where U_T asks for the identity on T^∨/T = L/2L; the quotient embeds in the group of isometries of that discriminant form, which is O(q_N) up to sign. ∎

In the data [E]: for the 24 polynomials of the enumerated surfaces with maximal order, the number of lattices, which equals the number of surfaces (note 2), is exactly |G|. The primes ramified over E₀ are described in §3.4.

### 3.2 Which polynomials occur

Write P(t) = t³h(t + q²/t), h(x) = x³ − e₁x² + (e₂ − 3q²)x − (e₃ − 2q²e₁). Ordinary means p ∤ e₁, q | e₂, q² | e₃.

**Brauer groups.** The Artin–Tate formula for a K3 surface over F_q whose Néron–Severi group, of rank 16 and discriminant 64, is defined over F_q gives P(q) = 64·q⁵·|Br(X)|, and the order of the Brauer group is a square [M]. Both were read at a source in [note 18](18-checks-of-notes-13-to-16.md) §5, with no error found: Kedlaya–Sutherland (https://arxiv.org/abs/1511.06945, Theorem 8, after Elsenhans–Jahnel and Milne) give L₁(1) = |Δ|·#Br(X) with #Br(X) a perfect square, for a K3 surface that satisfies the Tate conjecture; with their normalisation this is P(q) = 64q⁵·|Br(X)| when the Néron–Severi group over F_q has rank 16 and discriminant 64. That holds for the surface and for its twist: the sixteen classes are rational on both, and the covering involution acts by −1 on the transcendental part. The quadratic twist has polynomial P(−t). So

> P(q) = 64 q⁵ μ²,  P(−q) = 64 q⁵ ν²,  that is  h(2q) = 64q²μ²,  h(−2q) = −64q²ν²,

with μ² = |Br(X)| and ν² the same for the twist. Given (e₁, μ, ν) this determines P: e₂ = q(16μ² + 16ν² − q), e₃ = −2q²(16μ² − 16ν² + e₁).

**Observation.** [E] Call P *admissible* if it has this shape with e₁ > 0, is ordinary and irreducible, h has three real roots in (−2q, 2q), and (F − 1)/2 is integral (the split condition of Proposition C).

- For q = 13, 17, 19, 23, 29, 31, 37 the admissible polynomials are exactly the 44 polynomials of the enumerated surfaces (up to t ↦ −t).
- For every q up to 101 treated, an admissible polynomial has R_min = O_E only if μ = ν = 1.

The 23 values of q are the primes from 13 to 101 and 25, 49, with 1,264 admissible polynomials in all: from 1 at q = 13 (0.006·q²) to 188 at q = 101 (0.018·q²). Of these, 309 have μ = ν = 1 (32 at q = 101), and 270 of the 309 have maximal order. For the others the minimal order is not maximal: the isogeny class is a stack of levels, one for each order, and the steps between levels are the ring-changing steps of Proposition H at the primes dividing μν. This is the K3 form of t² − 4q = v²D for elliptic curves, with μν in the role of v.

### 3.3 Size

[E] For every q of the first block the script lists all admissible polynomials; for the larger q it draws polynomials with μ = ν = 1 at random. G is computed for those with maximal order.

| q | polynomials | with maximal order | \|G\|: min, median, max | median \|G\|/q² | surfaces per polynomial (data) |
|---|---|---|---|---|---|
| 13 | 1 | 1 | 1, 1, 1 | 0.0059 | 1.00 = 0.0059·q² |
| 17 | 2 | 2 | 1, 1.5, 2 | 0.0052 | 1.50 = 0.0052·q² |
| 19 | 2 | 2 | 2, 2, 2 | 0.0055 | 2.00 = 0.0055·q² |
| 23 | 4 | 3 | 2, 4, 6 | 0.0076 | 3.50 = 0.0066·q² |
| 25 | 4 | 2 | 2, 5, 8 | 0.0080 |  |
| 29 | 9 | 3 | 4, 6, 12 | 0.0071 | 5.00 = 0.0059·q² |
| 31 | 10 | 6 | 2, 8, 10 | 0.0083 | 5.80 = 0.0060·q² |
| 37 | 16 | 7 | 3, 10, 26 | 0.0073 | 9.56 = 0.0070·q² |
| 41 | 21 | 6 | 4, 12, 26 | 0.0071 |  |
| 43 | 24 | 10 | 2, 16, 28 | 0.0087 |  |
| 47 | 28 | 10 | 4, 20, 34 | 0.0091 |  |
| 49 | 28 | 10 | 4, 26, 36 | 0.0108 |  |
| 53 | 41 | 11 | 12, 24, 44 | 0.0085 |  |
| 59 | 52 | 13 | 10, 30, 52 | 0.0086 |  |
| 61 | 56 | 14 | 8, 35, 56 | 0.0094 |  |
| 67 | 70 | 18 | 16, 38, 84 | 0.0085 |  |
| 71 | 82 | 18 | 10, 48, 92 | 0.0095 |  |
| 73 | 89 | 17 | 6, 48, 82 | 0.0090 |  |
| 79 | 104 | 23 | 12, 56, 96 | 0.0090 |  |
| 83 | 118 | 19 | 18, 64, 116 | 0.0093 |  |
| 89 | 142 | 22 | 12, 74, 148 | 0.0093 |  |
| 97 | 173 | 28 | 20, 80, 184 | 0.0085 |  |
| 101 | 188 | 25 | 28, 100, 246 | 0.0098 |  |
| 211 | sampled | 40 | 12, 406, 968 | 0.0091 |  |
| 503 | sampled | 40 | 556, 2,845, 5,808 | 0.0112 |  |
| 1,009 | sampled | 40 | 1,668, 11,658, 19,360 | 0.0115 |  |
| 2,003 | sampled | 30 | 10,256, 38,498, 97,502 | 0.0096 |  |
| 5,003 | sampled | 30 | 55,176, 255,077, 546,516 | 0.0102 |  |
| 10,007 | sampled | 30 | 104,942, 872,620, 2,367,144 | 0.0087 |  |
| 30,011 | sampled | 20 | 1,410,472, 9,762,270, 16,730,664 | 0.0108 |  |
| 100,003 | sampled | 12 | 19,135,960, 83,306,780, 173,629,952 | 0.0083 |  |

- Over the 512 groups, |G|/q² has median 0.0092 and quartiles 0.0067 and 0.0128 (recounted from `sixlines_groupsize.jsonl` in [note 27](27-hypothesis-ledger.md) §1); it ranges from 0.0003 to 0.024. The values 0.0095, 0.0068 and 0.0129 first given here are those of the 486 groups with q ≥ 41, the population of the first fit below (0.00952, 0.00680, 0.01289, with the same range); over all 636 maximal records of the file, which include the 24 "real" duplicates and q = 125, 343, 625, 2,401, 3,125, they are 0.0095, 0.0067 and 0.0129. In the table, the median |G| for q = 17 is 1.5 (groups of orders 1 and 2), and the ratio 0.0052 in that row is computed from it (note 27 §1).
- A least-squares fit of log|G| against log q gives exponent 2.01 over q ≥ 41 (486 groups) and 1.98 over q ≥ 211 (242 groups), with constant 0.011 in the second.
- The last column is the average number of surfaces per polynomial among the enumerated surfaces, over all orders: about 0.006·q².
- Wherever the lattices were counted as in note 2, 191 polynomials including the 24 of the surfaces with maximal order, their number is exactly |G|: one genus.

**Heuristic for other families.** [H] The class number formula gives |G| of the order of √(d_E/d_{E₀}) up to logarithms. Here d_{E₀} is at most of the order of q⁴, because one root of h is a p-adic unit and the two others are divisible by q, and E/E₀ is ramified at 2 only; so |G| is of the order of q². The moduli space has dimension 4 and about q² polynomials occur, so q⁴ surfaces are spread over q² classes of size q², which agrees. The same count for transcendental rank r gives |G| ≈ q^{(r − 2)/2} = q^{(20 − ρ̄)/2}: the square root of the number of points of the moduli space, as for elliptic curves.

### 3.4 Structure

[E] Over the 512 groups:

- |G| is even in 509. The number of cyclic factors of even order is 1 in 297, 2 in 176, 3 in 34, 4 in 2, and 0 in 3.
- The odd part of G is cyclic in 506; G itself is cyclic in 296.
- For q ≥ 1,000 the largest prime factor of |G| is about |G|^0.55 in the median.
- In all 512 fields exactly one prime of E is ramified over E₀. It lies above 2, and the relative discriminant has norm 64.

No Pohlig–Hellman reduction is known for a group action without invariants of the orbits of a subgroup [H], so the factorisation of |G| matters only through §6.

---

## 4. The transfer to abelian fourfolds

### 4.1 Over the complex numbers

Read in Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes* [P]:

- The principally polarised abelian fourfolds A with an automorphism of order 4 of type (2, 2), that is multiplication by Z[i] preserving the polarisation, are parametrised by the domain H₂ of U(2, 2) (their §2.1).
- The second exterior power over C of H¹(A), for the complex structure given by i, carries an antilinear involution whose fixed part T(A) is a Hodge structure of rank 6 with Hodge numbers (1, 4, 1); and A ↦ T(A) induces an isomorphism of the moduli space of these fourfolds with D₄/SO*(T), T = U² ⊕ ⟨−1⟩² (their §2.2).
- For a generic six-line double plane X the transcendental lattice is T(2) (their Theorem 3.4), and the moduli space of six lines up to relabelling and correlation is that quotient (after Matsumoto, Sasaki, Yoshida).
- The Kuga–Satake variety of T(A) is isogenous to A⁴ (their Proposition 4.4, from Lombardo).

Not read here: Paranjape's paper, where the correspondence is shown to be algebraic (it is read and made explicit in [note 14](14-fourfold-of-a-six-line-surface.md) §2.1). According to the summaries in the two papers above and in Schlickewei's, A is the Prym variety of a cyclic cover of degree 4 of an elliptic curve by a curve of genus 5 built from the six lines [M].

This is the construction of [note 12](12-rank-six-galois-action.md), Step 0, with Q(i) in place of the unramified extension k: a quadratic space of rank 6 whose discriminant is a quadratic field is the fixed part of the second exterior power of a hermitian space of rank 4 over that field. For T = I_{2,4} the determinant is 1, so the field is Q(i).

### 4.2 Polynomials

**Proposition Y.** [D] Let P be ordinary with P(q) = 64q⁵μ² and P(−q) = 64q⁵ν², and

> χ(X) = X⁴ − 4(ν + μi)·X³ + e₁·X² − 4q(ν − μi)·X + q²  in Z[i][X].

1. The roots of P are the six products x_j·x_k, j < k, of the roots of χ.
2. The roots of χ have absolute value √q and χ·χ̄ has integer coefficients. If its middle coefficient is prime to p, it is the Weil polynomial of an isogeny class of ordinary abelian fourfolds over F_q (Honda–Tate [M]). If moreover χ·χ̄ is irreducible, Q(π) contains i, so these fourfolds have multiplication by Q(i), and Frobenius has characteristic polynomial χ over Q(i).
3. Changing the sign of μ replaces χ by χ̄, and changing the sign of ν replaces χ(X) by χ̄(−X), that is π by −π and χ by χ̄ together (corrected in [note 18](18-checks-of-notes-13-to-16.md), finding 5, from "replaces π by −π").

*Proof.* (1) For a quartic x⁴ − ax³ + bx² − dx + c² the polynomial of the products of two roots is t⁶ − bt⁵ + (ad − c²)t⁴ − (a²c² + d² − 2bc²)t³ + … (note 11, Step 3). With a = 4(ν + μi), d = q·ā, c = q, b = e₁ this gives e₂ = q·a·ā − q² = q(16μ² + 16ν² − q) and e₃ = q²(a² + ā²) − 2e₁q² = −2q²(16μ² − 16ν² + e₁), the values of §3.2. (2) The six products have absolute value q, so each root has absolute value √q; the coefficients of χ̄ are the conjugates of those of χ. The condition on the middle coefficient follows from P being ordinary: one root of P is a unit, so two roots of χ are units and the other two have the valuation of q, and χ·χ̄ has four unit roots (note 18, finding 5). If χ·χ̄ is irreducible, Q(π) and Q(i)[X]/χ are fields of degree 8, the first inside the second. ∎

Items 1 and 2 are re-derived in note 18, finding 5, with no error.

In the data [E]: (1), and the integrality and absolute values of (2), hold for the 44 polynomials of the enumerated surfaces and for 23,945 admissible polynomials, those of the primes from 41 to 257 and of q = 25, 49, 121, 125, 169, checked from the coefficients. In all of them the middle coefficient is prime to p; χ·χ̄ is irreducible for 32 of the 44 and 23,127 of the 23,945.

### 4.3 What this suggests, and what it does not prove

**Expected dictionary.** [D, sketch] The argument of note 12, Steps 0 to 4, with Z[i] for O_k, gives: unimodular lattices L with F correspond to self-dual hermitian Z[i]-lattices M of rank 4 with an endomorphism π, π·π̄ = q, up to scalars and up to conditions at 2. By Deligne's description of ordinary abelian varieties [M], such (M, π) are the principally polarised ordinary abelian fourfolds with multiplication by Z[i] in the isogeny class of Proposition Y. Two lattices are ℓ-neighbours exactly when ℓM ⊂ M′ ⊂ M with M′ self-dual for the form divided by ℓ: an isogeny A → A′ with kernel (Z/ℓ)⁴, maximal isotropic in A[ℓ] and stable under i. The group of Lemma 9 is then, up to 2-groups, a quotient of the corresponding unitary class group of the field Q(i)[π] of degree 8. This is made explicit and checked in [note 15](15-fourfold-map-and-the-group-action.md): a homomorphism Cl(K) → Cl(E), K = Q(i)[π], given on primes, induces Cl(K)/Cl(K⁺) → G with kernel and cokernel of order 1 to 8 (powers of 2) on 278 polynomials, and the odd parts are isomorphic. The isogeny with kernel A[𝔓₁𝔓₂·conj(𝔓₃𝔓₄)] corresponds to the step at the prime of E at the root a₁a₂.

**Not proved.** That the surface X over F_q and the fourfold are related by a correspondence defined over F_q, so that the lattices of Taelman's theorem and of Deligne's match; the conditions at 2; non-maximal orders. I checked polynomials only here, and built no fourfold. The fourfold is built in [note 14](14-fourfold-of-a-six-line-surface.md) for surfaces over F_p, p ≤ 37, and its Frobenius polynomial agrees with Proposition Y in every case with irreducible P. The matching of lattices is not proved or tested there.

**If it holds:**

- *For the attacker.* Vectorisation for six-line surfaces reduces to vectorisation for these fourfolds, once X ↦ A is computed; over C there are explicit formulas by theta constants (Matsumoto [M]).
- *For the honest party.* Isogenies of principally polarised abelian varieties with kernel (Z/ℓ)^g can be computed with theta functions, in time polynomial in ℓ^g [M: Lubicz, Robert]. So the step that has no equations on the K3 side has them, in principle, on the fourfold side.
- *For the question.* Nothing specific to K3 surfaces remains in this family. It is a class-group action on abelian varieties of dimension 4 with extra structure, as the starting hypotheses ([../SUMMARY.md](../../SUMMARY.md) §1.1, hypothesis (b)) state for ρ̄ ≥ 17 through abelian surfaces [M].

---

## 5. Attacks and costs, given the action

All of this assumes item 3 of §2 at unit cost.

- **Quantum.** The group is abelian for every family (Proposition A), so vectorisation is a hidden shift problem; Kuperberg's algorithm uses 2^{O(√log|G|)} evaluations [M]. This is the CSIDH situation, and it does not depend on §4.
- **Classical.** Meet in the middle or a random walk: about √|G| ≈ 0.1·q evaluations [M]. Listing the family costs q⁴.
- **Parameters.** [H] With |G| ≈ 0.01·q², a group of 2²⁵⁶ elements needs q ≈ 2¹³¹, and a surface is named by four elements of F_q, about 525 bits. CSIDH-512 has the same group size with 512 bits.
- **Knowing P.** For q = p of that size no method is known: computing the zeta function of a K3 surface costs at least about √p [M: Harvey; Costa, Harvey, Kedlaya]. For p small and q = p^f the p-adic methods are polynomial in f [M]. So instances of that size would have small characteristic.
- **One step.** On the fourfold side, of the order of ℓ⁴ operations in the field where the kernel is defined, after finding that kernel in A[ℓ]; against of the order of ℓ for an elliptic curve.

So even with everything granted, this is CSIDH with slower steps, the same generic attacks, and no gain.

---

## 6. The decision problem

Genus characters (Castryck, Sotáková, Vercauteren for elliptic curves) are the characters of order 2 of the group that can be computed from the object. Carried over to this family, that description is misleading ([note 19](19-quadratic-characters-from-a-surface.md)): the characters that Castryck, Sotáková and Vercauteren compute separate genera, and the surfaces of one polynomial are a single genus, because one finite prime ramifies; G is the analogue of Cl(O)², and its characters of order 2 are the analogue of characters of order 4. So that method has nothing to read on the K3 lattice.

- [D] By Lemma 9 (1), for a maximal order E/E₀ is unramified at odd primes. The characters that genus theory provides come from the places of E₀ ramified in E, and number at most one less than these places: the three real places and the dyadic primes that ramify. So their number is bounded whatever q. They are not all the characters of order 2 of G ([note 18](18-checks-of-notes-13-to-16.md), finding 4): a quadratic character ψ₀ of Cl(E₀) that is trivial on the norm of the ramified prime gives the character ψ₀∘N_{E/E₀} of G, and the number of these is not bounded by the number of ramified places. Every character of order 2 of G is the character of an unramified extension E(√d)/E with d in E₀, and their number is at most 3 + the 2-rank of Cl(E₀) (note 19); on the 480 fields with q ≤ 10,007, d is a unit of E₀ for 437 of 699 characters and needs an ideal class of E₀ for the other 262. So a bound on the number of characters of order 2 of G that does not depend on q is not shown.
- [E] In all 512 fields of §3 exactly one dyadic prime ramifies, which leaves at most 3 characters attached to ramified places. G has 1 or 2 cyclic factors of even order in 473 cases and at most 4. G has more than one cyclic factor of even order only when h(E₀) is even (212 of 212) and at most one when h(E₀) is odd (258 of 258) (note 18, finding 4).
- Whether any of them can be computed from a surface. One can: the character χ₋₁ of d = −1, the degree of the step modulo 4, exists for 80 of 477 polynomials in note 19, all with q ≡ 1 mod 4, exactly when the real quartic field of the fourfold has a prime above 2 of odd degree; it is computed from a surface on three torsors in [note 24](24-the-form-mod-4-is-the-character.md) (see §9 of this note). For the others, on the 419 surfaces of the 36 polynomials with maximal order (note 19): the square classes of the quantities of note 14 (λ₅λ₆, (1 − λ₅)(1 − λ₆), the twist of E, the discriminant of the branch quadratic) are the same for all surfaces of a torsor; the conic invariant, two power sums and a product of quartic characters vary and are neither constant nor balanced on 14 to 24 torsors; the unit of c₁ is not an invariant of the surface. None computes a character. [Note 28](28-leakage-audit.md) tests the characters other than χ₋₁ against full labels, on 9 torsors with groups Z/4, Z/6, Z/4 × Z/2 and Z/8, none of which has χ₋₁. The 4,838 candidate functions are built from note 19's invariants, from residue symbols of order 2 and 4 of the J-invariants and of power sums, from the curve data, and from relational symbols of two surfaces. None reads a character or any function of the class after correction (smallest exact p 0.0038, Holm 1). A character of order 2 read exactly on every torsor would have been found (200 of 200 plants). The negative result is for characters of order 2: a quartic residue symbol of a quantity over F_q exists only for q ≡ 1 mod 4, which is 3 of the 9 torsors, and no candidate could reach the threshold with an order-4 reading, so characters of order 4 were tested only through their squares (note 28 §1). Characters of order 3, which exist on two of the 9 torsors, are tested in [note 29](29-characters-of-order-3.md): none is read on the two torsors where a cubic symbol exists, and no candidate separates the cosets of its kernel on any of the seven torsors with 3 ∣ |G|.
- For non-maximal orders further characters come from the primes dividing μν, where the Brauer group is not trivial; notes 11 and 12 describe the lattices there, not characters.

So one bit of the class leaks this way, χ₋₁, for the polynomials where it exists; for the other characters none is shown to, now as a measured negative within the families of notes 28 and 29.

---

## 7. The other half: can the class be read off a surface

Starting hypothesis (c) ([../SUMMARY.md](../../SUMMARY.md) §1.1) is that it cannot, by analogy with elliptic curves: the genus is local data, and the class depends on the embedding ι and the canonical lift [H]. Nothing found since changes this, in either direction:

- The data one can compute from X over F_q, the zeta function, the Néron–Severi lattice with Frobenius, the Brauer group, is the same for all surfaces of a torsor (note 3 §4).
- Recovering the class through the canonical lift and its complex periods amounts to choosing a root of a class polynomial whose degree is the size of the torsor.
- The fourfold of §4 moves the question to abelian varieties with complex multiplication, where it is the usual vectorisation problem.

I looked for no attack beyond these remarks, and the absence of one is not evidence.

A first test of this half, at small sizes, is in [note 28](28-leakage-audit.md). On 9 torsors with |G| ≤ 8 and full labels, no pair function f(X, Y) of the surface invariants in the families of note 28 is a function of the class g with Y = g·X beyond chance. On the coset torsors the only bit read is χ₋₁, through δ. The test would have found an exact relation f(X, Y) = h(g) (200 of 200 plants). It cannot see the label as a function of the complete invariants through the class polynomial, which is the honest party's own computation, and it says nothing about cost at cryptographic size. [Note 29](29-characters-of-order-3.md) extends the negative result to characters of order 3: as characters on the two torsors at q = 37, and by separation of the cosets on seven torsors with |G| ≤ 12. Neither changes the conclusion above.

---

## 8. Check [E]

`sixlines_groupsize.py` works with PARI on the polynomials only.

- **real.** For the 44 polynomials of `sixlines_cm.jsonl`, which come from enumerated surfaces with p = 13 to 37: asserts that h(2q) and −h(−2q) are 64q² times squares; for the 24 with maximal order computes G and asserts that the count of lattices of note 2 is reproduced; it equals |G| in all 24.
- **all.** For q = 25, 49 and the primes from 13 to 101: runs through all (e₁, m, n) with h(2q) = 4q²m², h(−2q) = −4q²n², e₁ > 0, and keeps the admissible ones, 1,264 in all. Asserts that m and n are multiples of 4; that for the seven primes with data the admissible polynomials are exactly those of the surfaces; and that the order is maximal only if μ = ν = 1. For the 270 with maximal order computes G; when the class number of E is at most 300 it also counts the lattices and asserts that the count is a multiple of |G| (190 cases, equal in all).
- **sample.** For q = 211, 503, 1,009, 2,003, 5,003, 10,007, 30,011, 100,003: polynomials with μ = ν = 1 and e₁ drawn at random, kept if admissible with maximal order, 12 to 40 for each q; G for each.
- **fourfold.** For the polynomials of the first item and all admissible polynomials of the given q: asserts Proposition Y (1) by computing the characteristic polynomial of the second exterior power of the companion matrix of χ over Z[i], and that χ·χ̄ has integer coefficients and roots of absolute value √q.

All assertions pass.

**Limits.**

- Class groups are computed by PARI without certification, so they are correct under GRH.
- For q above 101 only polynomials with μ = ν = 1 are drawn, on the strength of the observation for q ≤ 101; and, the class numbers being above 300, the existence of a lattice with the right signature is checked for one of the 242 only.
- For q = p^f with f ≥ 3 the script inherits from `sixlines_cm.py` an approximation of the p-adic generator that is too coarse, so "maximal order" is not tested reliably there. The run includes q = 125, 343, 625, 2,401 and 3,125; I left them out of the table and the statistics. Their median |G|/q² is between 0.010 and 0.014.
- Maximal orders only. The constant 0.006 for whole isogeny classes comes from 44 polynomials with p ≤ 37.
- The fourfold check is an identity between polynomials. It says nothing about surfaces or fourfolds beyond their zeta functions.

---

## 9. What this does not give

- The correspondence of §4 over F_q. It is tested on surfaces, at the level of Frobenius polynomials, in [note 14](14-fourfold-of-a-six-line-surface.md), and not proved; until it is proved, items 3 and 4 of §1 are a conditional statement about this family.
- Anything for ρ̄ < 16, where no algebraic Kuga–Satake correspondence is known in general and the transfer argument is not available. This is examined on paper in [note 17](17-below-picard-16.md). (1) Taelman's theorem is unconditional for ρ̄ ≥ 12, and his Theorem 4.6 also when the Picard lattice contains a hyperbolic plane (p ≥ 5) or there is an ample class of degree below p − 4, at any Picard number; so the torsor exists on paper for ρ̄ = 14 and 12 as for 16, and for elliptic surfaces down to ρ̄ = 2. (2) Below 16 an algebraic correspondence is known only for surfaces with an automorphism of order 4 or 3 (cyclic quartics, ρ = 8; one family with ρ = 2); at ρ = 14 the abelian side is identified on periods as 8-folds with quaternionic multiplication (Poon), without a construction. (3) The simple factor of the Kuga–Satake variety of a surface over F_q has dimension 2^{(20−ρ̄)/2}: 8, 16, 32, …, 512 for ρ̄ = 14, 12, 10, …, 2. (4) No family below 16 has a named surface, a zeta function and a step together; the nearest is elliptic surfaces at ℓ = 2 and 3, where van Geemen gives half of a step over C. (5) With |G| ≈ q^{(20−ρ̄)/2} an instance has about 512 bits for 2²⁵⁶ at every ρ̄, and the field shrinks to about 2²⁸ at ρ̄ = 2, where the zeta function over a prime field comes within reach. Verdict: no reduction known for ρ̄ ≤ 14, and no usable action, so no assumption.
- The group for non-maximal orders, and the size of a whole isogeny class beyond the average from the enumerated surfaces.
- A proof of the shape P(±q) = 64q⁵·(square) from the surfaces, which is quoted, and of "maximal only if μ = ν = 1", which is observed.
- Whether quadratic characters of G other than χ₋₁ can be computed from a surface. The natural candidates were tested in [note 19](19-quadratic-characters-from-a-surface.md) and none is a character; the square classes among them are functions of the polynomial (§6, where the wider audit of note 28 is also reported). For the character of d = −1, which note 19 left open as within reach of a pairing of order 4 on the fourfold, [note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) gives a candidate computed from a surface: a number modulo 4 from a Tate pairing on the Jacobian of the curve of note 14 over F_(p⁴), normalised by the quartic of the curve. On the three torsors over F₂₉ and F₃₇ where the character exists it does not depend on the pair of lines (32 surfaces of 32) and is +1 on exactly half of each torsor. [Note 24](24-the-form-mod-4-is-the-character.md) shows that it is the character on the three torsors (32 surfaces of 32), against the action computed by the CM method of note 20, and [note 25](25-checks-of-notes-20-23-24.md) checked this and found that the claim stands. So for the polynomials where χ₋₁ is a character (80 of 477 in note 19, all with q ≡ 1 mod 4) the decision problem is broken by one bit: one bit of the class leaks there, by the method of Castryck, Sotáková and Vercauteren carried to the fourfold. The remedies are the countermeasures of §5.3 of their paper for class group actions that have such a character; they are not "the same as for CSIDH", since the paper says its attack does not affect CSIDH, which uses p ≡ 3 mod 4 and a class group without 2-torsion (note 25). The search problem is not touched. The evidence is three torsors at p = 29 and 37, and [note 26](26-delta-is-the-character-at-q-41.md) adds two at q = 41: 60 of 60 surfaces in five torsors at three primes. Note 23 also tests the correspondence of §4 beyond Frobenius polynomials: the endomorphism ring of the fourfold at ℓ separates the surfaces of one polynomial exactly as the levels of note 7 predict (29 surfaces, 3 polynomials, 422 pairs, no exception).
- Any experiment on hardness. No step is computed in this note. In another family a step at 2 is computed on the K3 side in [note 22](22-two-neighbours-of-elliptic-k3.md), as a pair: the two Jacobians of a double cover of P¹ × P¹ branched along a (4, 4) curve, for elliptic K3 surfaces with ρ̄ = 12 and 10 over F_5 and F_7. The local structure at 2 and the non-isomorphism agree with the lattice prediction in all cases; that the step is the class of 𝔩/𝔩̄ is not tested. The heuristic of §3.3 is measured there: |G| ≈ 0.08·p⁴ and 0.065·p⁵. In this family steps are computed in [note 20](20-cm-method-and-the-action.md), for groups of order at most 8: on lattices of the fourfold, and carried to the surfaces over F_p by class polynomials (the CM method, with Matsumoto's theta constants). In 15 fields the group that acts is G of Lemma 9, and 67 enumerated surfaces carry a label in G. This is not an experiment on hardness: the method needs a class polynomial of degree |G|. The lattices are even hermitian Z[i]-lattices over an order of index 8 in the maximal order of K (64 in one field), not over the maximal order.
- A check of the statements quoted from memory, in particular the cost of zeta functions in small characteristic and of isogenies in dimension 4.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_groupsize.py --cap 300 --out sixlines_groupsize.jsonl --fourfold 41 61 101 --real --all 13 17 19 23 25 29 31 37 41 43 47 49 53 59 61 67 71 73 79 83 89 97 101 --sample 211:40 503:40 1009:40 2003:30 5003:30 10007:30 30011:20 100003:12 125:20 625:20 3125:20 343:20 2401:20 > sixlines_groupsize.log
```

It takes about seventeen minutes, most of it for q from 60 to 101 and for the two largest primes. The sampling is seeded, so the output is the same on each run. The other values of q for the fourfold identity print one line of counts and take about a quarter of an hour:

```bash
python -I sixlines_groupsize.py --fourfold 25 43 47 49 53 59 67 71 73 79 83 89 97 103 107 109 113 121 125 127 131 137 139 149 151 157 163 167 169 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257
```

I ran that list in nine pieces and did not keep their output; together they give 23,680 polynomials, 22,902 of them with χ·χ̄ irreducible. The 23,945 and 23,127 of §4.2 are these plus the 265 and 225 of q = 41, 61, 101 in the first line of the stored `sixlines_groupsize.log` (the 44 and 32 there are the enumerated surfaces).

## Sources

- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, https://arxiv.org/abs/1209.5997 : read, §§1 to 2.2, Theorem 3.4, §4.3 [P].
- Lombardo, *Abelian varieties of Weil type and Kuga–Satake varieties*, https://arxiv.org/abs/math/0211224 : abstract and introduction read [P].
- Paranjape, *Abelian varieties associated to certain K3 surfaces*, Compositio Math. 68 (1988), https://numdam.org/item/CM_1988__68_1_11_0/ : not read here; described after the two papers above and Schlickewei, *The Hodge conjecture for self-products of certain K3 surfaces*, https://arxiv.org/abs/0907.2503 [M]. Read in [note 14](14-fourfold-of-a-six-line-surface.md) §2.1.
- Recalled, not re-read [M]: Matsumoto, Sasaki, Yoshida on the period map of the six-line family, and Matsumoto on its theta functions; the Artin–Tate formula and the squareness of the order of the Brauer group (Milne; Liu, Lorenzini, Raynaud); Honda–Tate theory and Deligne's description of ordinary abelian varieties; Kuperberg's algorithm; Castryck, Sotáková, Vercauteren on genus characters; Lubicz and Robert on isogenies by theta functions; Harvey, and Costa, Harvey, Kedlaya, on zeta functions.
- Kedlaya, Sutherland, https://arxiv.org/abs/1511.06945 , Theorem 8: the Artin–Tate formula for K3 surfaces, read in [note 18](18-checks-of-notes-13-to-16.md) §5 [P].
- In this repository: notes [2](02-six-line-count.md), [3](03-torsor-action.md), [6](06-torsors-under-cm-class-groups.md), [11](11-rank-six-fourth-power.md), [12](12-rank-six-galois-action.md); later results cited in the text are from notes 14, 15, 17 to 20, 22 to 25 and 27 to 29. The starting hypotheses are stated in [../SUMMARY.md](../../SUMMARY.md) §1.1.

