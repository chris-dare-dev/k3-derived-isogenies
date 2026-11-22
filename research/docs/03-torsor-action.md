# 3. The action on the torsor

*Tests the group action on the torsors of [note 2](02-six-line-count.md) (§6 item 1 there, "label the torsor") in the one direction where it can be computed on both sides, Frobenius conjugation over F₂₅ and F₄₉. Labels: **[P]** read in the cited source; **[D]** derived here; **[E]** measured by the scripts in `experiments/`.*

---

## 1. Result

Note 2 showed that surfaces and lattices are equally many. This note tests whether the **group action** agrees, in the one direction where it can be computed on both sides.

- **Lattice side, computed.** For every polynomial of note 2, the permutation of the lattice classes by the step 𝔩/𝔩̄ was computed for each Elkies prime ℓ ≤ 19, and by the step 𝔭/𝔭̄ at the unit-root prime 𝔭 above p. [E]
- **Surface side, one step.** An ℓ-neighbour step between two explicit surfaces cannot be tested with what exists today (§4). But the 𝔭/𝔭̄ step has a geometric meaning: it is **Frobenius conjugation** X ↦ X^(p), by the main theorem of complex multiplication applied to the canonical lift (§2). Over F_p it is the identity. Over F_{p²} it is an involution on the surfaces with a given Frobenius polynomial, and it can be computed.
- **Comparison over F₂₅ and F₄₉.** All six-line double planes over F₂₅ (365 classes) and over F₄₉ (6,246 classes) were enumerated. For all 32 ordinary irreducible Frobenius polynomials (4 over F₂₅, 28 over F₄₉, one of them computed with the faster lattice search of [note 5](05-non-split-count.md) §6), the number of surfaces equals the number of lattices, and **the cycle type of Frobenius conjugation on the surfaces equals the cycle type of the 𝔭/𝔭̄ step on the lattices in all 32**. [E] The F₂₅ cases are below; F₄₉ is in §6.

| polynomial over F₂₅ (e₁, e₂, e₃) | surfaces | lattices | order index | conjugation on surfaces | 𝔭/𝔭̄ on lattices |
|---|---|---|---|---|---|
| (2, 175, −2500) | 8 | 8 | 1 | 2, 2, 2, 2 | 2, 2, 2, 2 |
| (6, 175, −7500) | 6 | 6 | 25 | 2, 2, 2 | 2, 2, 2 |
| (14, 175, −17500) | 2 | 2 | 1 | 1, 1 | 1, 1 |
| (26, 1375, 27500) | 5 | 5 | 16 | 1, 1, 1, 2 | 1, 1, 1, 2 |

The four cases are different in kind: a free involution, a free involution on a torsor that lives entirely over non-maximal orders at 5, a trivial class, and a mixed case where three surfaces are defined over F₅ up to isomorphism and two are exchanged. The lattice computation reproduces each. This is the first check of the action itself, not only of the count. F₄₉ adds 28 polynomials and 540 surfaces, with the same outcome (§6).

---

## 2. Why Frobenius conjugation is the 𝔭/𝔭̄ step [D]

Let X be ordinary over F_q, q = p^f, with canonical lift X_can over W(F_q), and ι: W(F_q) → C. The conjugate surface X^(p) has canonical lift (X_can)^σ, with σ the Frobenius of W(F_q). So its lattice is that of the complex surface Z^τ, where Z = X_can ⊗_ι C and τ ∈ Aut(C) extends σ.

- E embeds in Q_p through the unit root, so τ fixes E, and on E^ab it lies in the decomposition group at the unit-root place v, over the Frobenius.
- By the main theorem (Rizov; Valloni's lattice form, note 1), the lattice of Z^τ is (s/s̄)·T(Z) for an idèle s with art(s) = τ. Here s is a uniformiser at v, up to units, and units act trivially because the order contains the idempotents at p (condition (M4)).
- So T(X^(p)) = x·T(X) with x = (p at v, p⁻¹ at v̄, 1 elsewhere). As a lattice, x·L = ξL + pL, where ξ ∈ R_min[1/p] approximates x p-adically.

Consistency: (F/q) is a norm-one element with divisor (𝔭̄/𝔭)^f times a unit at the slope-1 part, so the step has order dividing f, as conjugation does. In the data the step is the identity for all 40 polynomials over F_p of the first run of note 2, and for all 44 in the stored `sixlines_neighbours.log`. [E]

Frobenius conjugation is an isogeny in Yang's sense but not a derived isogeny: Bragg–Yang's slope condition excludes it [P, note 1]. So this tests the CM action at p, and says nothing yet about twisted Fourier–Mukai partners.

---

## 3. The neighbour graphs on the lattice side [E]

For a degree-one prime 𝔩 = (ℓ, F − λ) with λ a simple root of P mod ℓ and λ ≢ q²/λ (the Elkies condition of note 1 §3.2), the step is L ↦ 𝔩̄⁻¹(ℓL + (F − λ)L). Classes are compared by an exact isomorphism test: γL′ = L with γγ̄ = 1.

- **44 polynomials, 85 steps.** These are the 40 polynomials over F_p of the first run of note 2 and the 4 over F₂₅, and the counts of this section are of that run, whose output was not kept. The stored `sixlines_neighbours.log` and `.jsonl` are a later run that also has the other 4 polynomials over F_p: 48 polynomials, 88 steps. All four have non-maximal order, so the counts of maximal orders below are the same in it. Of the 297 pairs (polynomial, ℓ) with ℓ ∈ {3, …, 19}, ℓ ≠ p and ℓ prime to the index, 76 have an Elkies prime: about one in four (in the stored run, 79 of 321).
- **Maximal order (26 polynomials: 24 over F_p, 2 over F₂₅): every step is a translation.** All cycles of a given step have the same length. Examples at p = 29, 12 surfaces, class number 48: ℓ = 7 and ℓ = 19 give a 12-cycle, ℓ = 11 gives three 4-cycles, ℓ = 17 gives four 3-cycles. At p = 37, 26 surfaces: ℓ = 5, 7, 17 give a 26-cycle and ℓ = 11 two 13-cycles. So these groups are cyclic of order 12 and 26, and one small prime generates them.
- **Connectivity.** Steps with ℓ ≤ 19 connect the whole torsor in 20 of the 26 maximal cases. In the other 6 there is no Elkies prime ≤ 19, or too few.
- **Non-maximal orders (18 polynomials: 16 over F_p, 2 over F₂₅; 20 over F_p with the other four): the steps preserve the order.** The classes split into levels by endomorphism ring, and ℓ-steps with ℓ prime to the index act inside each level. Example at p = 37, 14 surfaces, index 4: two levels of 7, each a 7-cycle under ℓ = 11, 17, 19. No ℓ-step prime to the index connects two levels; that needs steps at the primes dividing the index (here 2), which are the vertical edges of the "volcano" and are not computed here. The steps that join the levels at an odd prime ℓ ≠ p are worked out in [note 7](07-ring-changing-steps.md): they are the eigenlines of F on T/ℓT that are not stable under End(T). In this note's data they occur at a double root of P mod ℓ, where the levels form an ℓ-volcano.

**What this labels.** The graph is a statement about lattices. Transferred to surfaces through the dictionary it says, for example, that the 12 surfaces at p = 29 form a single 12-cycle of twisted derived equivalences under ℓ = 7. It does not say which surface sits where. For a torsor of two elements there is nothing to choose, so the labelling is forced: that covers 10 of the 44 polynomials. For three or more, an independent geometric test of one step is needed.

---

## 4. Why ℓ-steps cannot be tested on surfaces yet

- **No modular polynomial.** For elliptic curves one tests ℓ-isogeny with Φ_ℓ(j, j′). The analogue here is a Hecke correspondence on the 4-dimensional moduli space of six lines. It is not known explicitly, and it cannot be interpolated from the data.
- **No explicit partner.** A step is two twisted equivalences, X ← Z → Y. Explicit Z exist for ℓ = 2 (nets of quadrics, van Geemen), and ℓ = 2 is excluded (note 1 §4.4). The diagonal net Σaᵢuᵢ² = Σbᵢuᵢ² = Σcᵢuᵢ² = 0 does give a canonical degree-8 partner Z_C of X_C, but its quotient by sign changes is the associated configuration, which is the same K3 surface. So it realises a 2-step down and no neighbour. [D] For an Elkies prime ℓ not dividing the discriminant or the conductor, [note 6](06-torsors-under-cm-class-groups.md), Proposition F, describes the partner, though not by equations: the ℓ-neighbour over F̄_p is the moduli space of stable α-twisted sheaves of twisted Mukai vector (ℓ, 0, 0), with α the Brauer class of order ℓ whose kernel is 𝔩T, and the two surfaces are related by one twisted derived equivalence. The argument rests on Bragg–Yang and has not been independently verified. It gives no test on the enumerated surfaces.
- **No invariant.** All surfaces of a torsor have the same NS, the same T as a Z-lattice (it is I_{2,4}(2) in every case, note 2 §2) and the same zeta function. Only the F-structure differs.
- **CM method.** Computing class polynomials for a moduli coordinate, by complex periods, p-adic canonical lifts or the Chinese remainder theorem over many p, would label everything. Each needs far more than the enumerations done here.

The honest position: the action is confirmed at p, by Frobenius conjugation, and is a prediction at ℓ ≠ p.

A possible route to a test goes through abelian fourfolds ([note 13](13-hardness-on-paper.md) §4): the Frobenius polynomial of a six-line surface is the second exterior power over Z[i] of the polynomial χ of an abelian fourfold with multiplication by Z[i] (Proposition Y there, checked on polynomials), and an ℓ-step should be an isogeny of such fourfolds with kernel (Z/ℓ)⁴. The fourfold has been built from the six lines over F_p and has the predicted Frobenius ([notes 14](14-fourfold-of-a-six-line-surface.md) and [15](15-fourfold-map-and-the-group-action.md)). An ℓ-step still cannot be tested: it would need an isogeny with kernel (Z/ℓ)⁴ of a Prym variety of a curve of genus 5, which was not computed. Note 15 §3.3 says which isogeny corresponds to which step, and §5 lists what a test would take.

---

## 5. Limits and next steps

- **What the test measures.** The cycle type of an involution is one number, the count of fixed points. It agreed in 32 of 32 polynomials (values 0 to 6 on torsors of 2 to 36 surfaces). It does not show that the bijection between surfaces and lattices is equivariant, only that the two involutions are conjugate.
- **One polynomial over F₄₉ not computed by the search used here** (discriminant group of order 4,096), as for the four of note 2. The faster search of note 5 §6 computes it (§6).
- **Order 2 only.** Over F_{p²} the step is an involution, so g and g⁻¹ are not distinguished. F_{p³} would test an element of order 3 with its orientation, but costs p⁶ times more per character sum.
- **Derived statements.** §2 is my derivation from the cited theorems; the sign convention (g or g⁻¹) is not fixed and does not matter for f = 2.
- **Next.** (a) F_{p³} for p = 5, with a faster character-sum kernel; and the non-split forms over F₇, to test the descent remark of §6 (done in [note 4](04-forms-over-the-prime-field.md)). (b) Vertical steps at 2, to connect the levels in non-maximal cases. (c) An ℓ = 3 step through elliptic fibrations: X and Y would need a common genus-one fibration of index 3; this is the route to a first explicit twisted partner.

---

## 6. F₄₉ [E]

**Enumeration.** 4,282,922 configurations over F₄₉ in general position, 6,246 classes, 85 minutes on 14 processes. Of these, 1,080 classes have an ordinary irreducible polynomial. They fall into 28 polynomials up to the twist and give 540 K3 surfaces. None is self-associated.

**Count.** The lattice count was computed here for 27 of the 28. **#K3 = N_total(τ) in all 27** (526 surfaces), for each of the three real places. 10 have maximal order and 17 do not. The indices are 2, 3, 4, 8, 9, 16 and 27, so the non-maximal part is now also seen at the prime 3 and not only at 2. In 5 of the 17 there is no lattice over the maximal order. The 28th polynomial, (2, 3871, −9604), is computed with the faster lattice search of [note 5](05-non-split-count.md) §6: 14 surfaces and 14 lattices. So the count agrees in 28 of 28 over F₄₉.

**Action at p.** **The cycle type of Frobenius conjugation on the surfaces equals that of the 𝔭/𝔭̄ step on the lattices in all 27.**

| kind | polynomials | detail |
|---|---|---|
| free involution (all surfaces in swapped pairs) | 21 | torsors of 4 to 36 surfaces |
| identity (every surface fixed) | 2 | (46, −833, −220892): 4 surfaces, maximal order, h(E) = 8; (118, 10927, 585844): 2 surfaces, index 8 |
| mixed | 4 | below |

| polynomial (e₁, e₂, e₃) | surfaces | index | fixed | swapped pairs | lattice side |
|---|---|---|---|---|---|
| (10, 1519, 182476) | 14 | 16 | 6 | 4 | 6 fixed, 4 pairs |
| (58, 1519, −48020) | 14 | 16 | 6 | 4 | 6 fixed, 4 pairs |
| (78, 5439, 240100) | 10 | 27 | 2 | 4 | 2 fixed, 4 pairs |
| (82, 5439, 220892) | 6 | 27 | 2 | 2 | 2 fixed, 2 pairs |

The 28th polynomial, (2, 3871, −9604), has 14 surfaces with 2 fixed and 6 pairs, and its 𝔭/𝔭̄ step, computed in note 5 §6, also has 2 fixed points and 6 swapped pairs. So the cycle type agrees in 28 of 28 over F₄₉, and in 32 of 32 with F₂₅.

**Neighbour steps.** 47 steps at Elkies primes ℓ ≤ 19, in 22 of the 27 polynomials. In the 10 maximal cases every step is again a translation: for example 36 surfaces, class number 144, where ℓ = 11 and ℓ = 19 give two 18-cycles and ℓ = 17 twelve 3-cycles. In the non-maximal cases 13 of 30 steps have cycles of unequal length, such as 8 + 24 or 4 + 4 + 12 + 12 on 32 surfaces. That is what levels of different sizes give: a step is a translation inside each level, by an element whose order depends on the level.

**A remark on the fixed surfaces.** [D] A surface fixed by conjugation has a model over F₇, up to a twist. Over F₇ no split six-line configuration has an irreducible ordinary polynomial (note 2: 3 classes, none). So the 22 fixed surfaces of the 27 polynomials above (24 with the 28th) would come from forms over F₇ on which Frobenius permutes the lines or exchanges the two polarisations. Those forms are not enumerated in note 2 (§5 there); counting them is an independent test of the fixed-point numbers above. [Note 4](04-forms-over-the-prime-field.md) does this: the forms over F₇ were enumerated (217 classes). Of the 24 fixed surfaces, 10 have a conjugation-stable configuration and exactly 20 forms; the other 14 are fixed because Frobenius exchanges the two polarisations. The lattice side predicts the Frobenius action on NS of every fixed surface, and it agrees in all 22 computed cases.

---

## Reproduce

From `experiments/`, with `python -I`:

```bash
python -I sixlines_fq.py selftest
```

```bash
python -I sixlines_fq.py enumerate --p 5 --out sixlines_q25.jsonl --procs 8
```

```bash
python -I sixlines_cm.py --list --neighbours --out sixlines_neighbours.jsonl sixlines_p{13,17,19,23,29,31,37}.jsonl sixlines_q25.jsonl > sixlines_neighbours.log
```

```bash
python -I sixlines_fq.py enumerate --p 7 --out sixlines_q49.jsonl --procs 14 --reduced > sixlines_q49.log 2>&1
```

```bash
python -I sixlines_cm.py --list --neighbours --out sixlines_q49_cm.jsonl sixlines_q49.jsonl > sixlines_q49_cm.log
```

The F₂₅ enumeration takes under a minute, the neighbour computation about three minutes. The F₄₉ enumeration takes 85 minutes on 14 processes; `--reduced` computes S₃ once per orbit of association and conjugation. Its lattice computation takes about five minutes.

## Sources

As in notes 1 and 2. For §2: Rizov, *Complex multiplication for K3 surfaces*, https://arxiv.org/abs/math/0508018 ; Valloni, https://arxiv.org/abs/1804.08763 ; Taelman, https://arxiv.org/abs/1711.09225 (the descent step of the essential-surjectivity proof uses the same reciprocity at p).
