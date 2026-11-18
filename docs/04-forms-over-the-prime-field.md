# 4. Forms over the prime field

*Finds the models over the prime field of the surfaces of [note 3](03-torsor-action.md) that are isomorphic to their Frobenius conjugate, and compares how Frobenius acts on their Néron–Severi lattice with the lattice prediction. Builds on the remark on fixed surfaces in note 3 §6 and on [note 2](02-six-line-count.md) §6 item 2 (non-split forms). Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by the scripts in `experiments/`.*

---

## 1. Result

Note 3 found 24 surfaces over F₄₉ and 5 over F₂₅ that are isomorphic to their Frobenius conjugate, and matched their number with the lattice side (22 of the 24 lie in polynomials whose lattices were computed there, the other 2 in the polynomial computed in [note 5](05-non-split-count.md) §6). This note looks at the models of those surfaces over the prime field.

- **Enumeration.** All sets of six lines over F₇ that are stable under Frobenius, with every line defined over F₄₉ and Frobenius moving at least one line, were enumerated up to PGL₃(F₇): **217 classes**. Over F₅: 56 classes. [E]
- **The fixed-surface count is confirmed without the lattice code.** Over F₄₉, 10 of the 24 fixed surfaces have a conjugation-stable configuration, and the F₇ enumeration finds exactly their 20 forms (two six-line polarisations each). The other 14 are fixed only because Frobenius exchanges the two polarisations, so they have no six-line model over F₇. [E]
- **The zeta functions agree.** For all 217 + 56 forms, the Frobenius polynomial over F_p, from point counts over F_p, F_{p²}, F_{p³}, has roots whose squares are the roots of the polynomial over F_{p²} from note 3. [E]
- **The lattices predict how Frobenius acts on NS.** [E] A lattice L is fixed by the 𝔭/𝔭̄ step exactly when ±F has a square root ϖ in End(L). Then ϖ is the Frobenius over F_p, and its action on L/2L is one of eight classes in O(q_N) = S₆ × Z/2. That class agrees with the one read off from the surface in **all 29 fixed surfaces** (24 over F₄₉, 5 over F₂₅, in 9 polynomials; the lattices of the two in (2, 3871, −9604) are computed with the faster search of note 5 §6), and no lattice has such a ϖ in the 23 polynomials without fixed surfaces.

| field | polynomials computed | with fixed surfaces | fixed surfaces | NS-action class agrees |
|---|---|---|---|---|
| F₂₅ | 4 | 2 | 5 | 5 |
| F₄₉ | 28 | 7 | 24 | 24 |

This is a stronger test than note 3's. There the invariant was one number per polynomial. Here each fixed surface carries a class with five observed values, and the multiset of classes agrees polynomial by polynomial. It is also the first check of the dictionary for surfaces where Frobenius acts non-trivially on NS.

---

## 2. Forms and their classification [D]

**Types.** A form is a Frobenius-stable set of six lines in general position, all defined over F_{p²}. Frobenius permutes the lines by an involution:

| type | cycle type | lines |
|---|---|---|
| A | (2, 1⁴) | four rational, one conjugate pair |
| B | (2², 1²) | two rational, two conjugate pairs |
| C | (2³) | three conjugate pairs |

**Normal form.** PGL₃(F_p) acts simply transitively on the frames (four labelled lines in general position) with a given Frobenius action. This is Lang's theorem for PGL₃; the script also counts the frames of types B and C and finds |PGL₃(F_p)| in each case (5,630,688 for p = 7). [E] So every form is a fixed frame of its type plus one conjugate pair, and the enumeration runs over the 2,394 non-rational points of the dual plane over F₄₉.

**Classes.** A labelled set of six lines has no automorphisms. So two labelled forms with the same Frobenius action are equivalent over F_p exactly when they are equivalent over F_{p²}. The class over F_p is the orbit of the split normal form (b, c, e, f) under the relabellings that commute with the involution; the class over F_{p²} is the orbit under all 720.

**Zeta function.** For a form with product of lines f₀ (coefficients in F_p), S_k = Σ χ(f₀) over P²(F_{p^k}) is the trace of Frobenius^k on T₆, as in the split case: the correction from the nodes cancels against the trace on N. For odd k a conjugate pair contributes the quadratic character of F_{p^{2k}} at one line of the pair.

---

## 3. The enumeration against the F_{p²} data [E]

| | p = 5 | p = 7 |
|---|---|---|
| classes of forms (A, B, C) | 56 (10, 26, 20) | 217 (37, 114, 66) |
| their classes over F_{p²} | 42 | 183 |
| conjugation-stable classes in the F_{p²} file | 45 | 190 |
| stable classes without automorphisms | 21 | 131 |
| forms per such class | exactly 1, in all 21 | exactly 1, in all 131 |
| stable classes with no form | 3 | 7 |
| forms whose class is not stable | 0 | 0 |
| polynomial over F_p squares to the one over F_{p²} | 56 of 56 | 217 of 217 |

- Every form lands on a conjugation-stable class, and every stable class without automorphisms has exactly one form. That is descent, and it confirms the conjugation map of note 3 by a separate enumeration.
- The 3 and 7 stable classes with no form all have an automorphism of order 2 (class size 360) and a reducible polynomial. Why they have no form was not examined.
- Association preserves the type, in all 42 and 183 classes: the associated configuration of a form of type A, B, C has forms of the same types.

**Forms with an irreducible ordinary polynomial over F₄₉: 20**, in three polynomials.

| polynomial over F₄₉ | fixed surfaces | forms over F₇ | types | polynomial over F₇ (up to t ↦ −t) |
|---|---|---|---|---|
| (10, 1519, 182476) | 6 | 12 | 4 A, 4 B, 4 C | t⁶ − 2t⁵ + 7t⁴ + 196t³ + 343t² − 4802t + 7⁶ |
| (2, 3871, −9604) | 2 | 4 | 4 B | t⁶ + 4t⁵ + 7t⁴ − 392t³ + 343t² + 9604t + 7⁶ |
| (118, 10927, 585844) | 2 | 4 | 2 B, 2 C | t⁶ − 8t⁵ + 91t⁴ − 392t³ + 4459t² − 19208t + 7⁶ |

So 10 fixed surfaces have 20 forms, two per surface. The other 14 fixed surfaces over F₄₉ (4 + 6 + 2 + 2, in four polynomials) have C^σ equivalent to the associated configuration and not to C.

In the first row, six surfaces over F₇ have the same Frobenius polynomial on T₆ and three different Frobenius actions on NS. Types B and C even have the same trace on NS (three fixed nodes each), hence the same zeta function over F₇, and differ only in the integral action.

---

## 4. The Frobenius action on NS, on both sides

**The eight classes.** [D] Let X₀ be a model over F_p of a fixed surface. Its Frobenius acts on N by p·g with g an involution, and g has an image in O(q_N). The script builds N from the 15 exceptional curves, h and the six line classes, and finds |O(q_N)| = 1440 with centre of order 2, generated by the line permutations and the central element z. [E] This replaces the [M] citation of note 2 §2. The involution classes are (π, 0) and (π, z), with π of type 1, A, B or C.

- **Surface side.** If C^σ is the relabelling of C by π, the class is (π, 0): "lines". If C^σ is the relabelling by π of the associated configuration with its natural labelling, the class is (π, z): "swap". The formula (b, c, e, f) ↦ (b, e, c, f) of note 2 is the Gale transform followed by the relabelling (14)(25)(36); the natural labelling is the one that commutes with relabelling, which the script asserts for every configuration it handles. That z is the image of association follows because association commutes with S₆ and is not induced by an automorphism.
- **Lattice side.** If L is fixed by the step x = 𝔭/𝔭̄, then xL = γ⁻¹L with γγ̄ = 1. Since x²L = (p²/F)L, the element γ²F/p² is a root of unity ζ, and ϖ = p/γ satisfies ϖ² = F/ζ and ϖL ⊂ L. When E has no fourth root of unity, ζ = ±1 and ±F = ϖ². Conversely a root ϖ of ±F in End(L) gives γ = p/ϖ. The action of ϖ on L/2L, with the quadratic form −(x·x) mod 4, is carried to (A_N, q_N) by any isometry and its class is read from the table of 1440 elements.

**Measured.** [E]

1. **Fixed ⇔ a square root of ±F in End(L)**, for every lattice class of the 31 polynomials (547 classes). In 7 polynomials over F₄₉, F or −F is a square in E. In one of them, (66, 5439, 297724) with index 9, the root lies in no End(L), and no surface is fixed.
2. **Classes agree** in all 8 polynomials with fixed surfaces whose lattices were computed here, and in the ninth, below the table:

| field, polynomial | index | fixed surfaces: classes (both sides) |
|---|---|---|
| F₄₉ (10, 1519, 182476) | 16 | 2 × (A, lines), 2 × (B, lines), 2 × (C, lines) |
| F₄₉ (118, 10927, 585844) | 8 | (B, lines), (C, lines) |
| F₄₉ (46, −833, −220892) | 1 | 4 × (B, swap) |
| F₄₉ (58, 1519, −48020) | 16 | 2 × (A, swap), 4 × (B, swap) |
| F₄₉ (78, 5439, 240100) | 27 | 2 × (B, swap) |
| F₄₉ (82, 5439, 220892) | 27 | 2 × (B, swap) |
| F₂₅ (14, 175, −17500) | 1 | 2 × (B, swap) |
| F₂₅ (26, 1375, 27500) | 16 | see below |

3. **The case Q(i) ⊂ E.** For (26, 1375, 27500) over F₂₅ both F and −F are squares, and the three fixed surfaces are self-associated. Each has two forms over F₅, with Frobenius ϖ and iϖ, one of each kind. Both sides give the same three pairs: {(A, lines), (B, swap)}, {(B, lines), (A, swap)}, {(C, lines), (B, swap)}.

The polynomial (2, 3871, −9604), whose lattices were not computed here, has 2 × (B, lines) on the surface side. With the faster lattice search of [note 5](05-non-split-count.md) §6, F is a square, 2 lattices are fixed, both of class (B, lines). So the classes agree for all 24 fixed surfaces over F₄₉, and for 29 of 29 with F₂₅.

---

## 5. What this adds, and limits

- **Non-split torsors.** A model over F₇ of a fixed surface is an ordinary K3 surface with ρ̄ = 16 on which Frobenius is not p on NS. Its lattice is (L, ϖ) with ϖ acting on L/2L by the class of g. The agreement in §4 says that the lattices with a given such action are as many as the surfaces, for these polynomials. This is the dependence on F|NS described in note 1 §4.3.
- **Not a full count over F₇.** Two models over F₇ of one surface over F₄₉ can differ by a cocycle in Aut(X), which is infinite. So the number of surfaces over F₇ up to isomorphism is not determined here; only the fixed surfaces over F₄₉ and their NS-action classes are. Forms with lines not all defined over F₄₉ (Frobenius of order 3, 4, 5 or 6 on the lines) are not enumerated.
- **Swap surfaces have no model in the family.** The 14 + 2 surfaces of kind "swap" are fixed only on the evidence of the configuration data; the F_p enumeration cannot see them. A model over F_p would need both polarisations, for example as a (2, 2)-type correspondence, and was not looked for.
- **Derived statements.** §2 and §4 are my derivations. The identification of z with association uses that the class of the association isometry is constant on the family; it is supported by the agreement and not proved here.
- **Small sample.** 27 fixed surfaces in 8 polynomials, classes (A, lines) 3, (B, lines) 4, (C, lines) 4, (A, swap) 3, (B, swap) 16, counting each form of the three self-associated surfaces; with the two (B, lines) surfaces of (2, 3871, −9604), 29 in 9 polynomials. The classes (1, swap) and (C, swap) did not occur.
- **Still no derived equivalence.** As in note 3, everything here is at p. The ℓ-steps remain a prediction.

**Next.** (a) Forms over F₁₁ and F₁₃ against an enumeration over F₁₂₁ and F₁₆₉: the F_p side is seconds, the F_{p²} side is the cost. Or drop the F_{p²} side: enumerate forms over F_p for larger p and compare with a direct count of lattices (L, ϖ) with prescribed action on L/2L, which is the non-split analogue of note 2. This second variant is [note 5](05-non-split-count.md), for p ≤ 23: the forms of each type are as many as the lattices predict, in all 136 polynomials of its first run and in all 169 with the faster search of its §6. (b) F_{p³}, as in note 3 §5.

---

## Reproduce

From `experiments/`, with `python -I`:

```bash
python -I sixlines_forms.py selftest
```

```bash
python -I sixlines_forms.py enumerate --p 7 --out sixlines_forms_p7.jsonl 2> sixlines_forms_p7.log
```

```bash
python -I sixlines_forms.py compare sixlines_forms_p7.jsonl sixlines_q49.jsonl --split sixlines_p7.jsonl --out sixlines_forms_p7_compare.jsonl > sixlines_forms_p7_compare.log
```

```bash
python -I sixlines_types.py sixlines_q49_cm.jsonl sixlines_q49.jsonl --out sixlines_types_q49.jsonl > sixlines_types_q49.log
```

For p = 5 replace 7 by 5, `q49` by `q25` and `sixlines_q49_cm.jsonl` by `sixlines_neighbours.jsonl`; the split file comes from `python -I sixlines.py enumerate --p 5 --out sixlines_p5.jsonl`. Each step takes under a minute. The inputs `sixlines_q25.jsonl`, `sixlines_q49.jsonl` and the two lattice files are those of note 3.

## Sources

As in notes 1, 2 and 3. For association as the Gale transform: Dolgachev, Ortland, *Point sets in projective spaces and theta functions*, ch. III [M]. For Lang's theorem: Lang, *Algebraic groups over finite fields*, Amer. J. Math. 78 (1956) [M].
