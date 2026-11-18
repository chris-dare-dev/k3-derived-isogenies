# 5. The non-split count

*Counts the six-line forms over F_p on which Frobenius permutes the lines by an involution, against a direct count of lattices: the non-split analogue of [note 2](02-six-line-count.md), and item (a) of [note 4](04-forms-over-the-prime-field.md) §5. Also gives the faster lattice search used in notes 2 to 4. Labels: **[P]** read in the cited source; **[D]** derived here; **[E]** measured by the scripts in `experiments/`.*

---

## 1. Result

Note 2 counted six-line double planes with all lines rational, where Frobenius is p on NS. This note counts the forms where Frobenius permutes the lines by an involution (types A, B, C of note 4), directly over F_p, and compares with a direct count of lattices. No data over F_{p²} is used.

| p | classes of forms | with ordinary irreducible polynomial | polynomials | lattice count computed | agree | forms in agreeing cases | lattices |
|---|---|---|---|---|---|---|---|
| 5 | 56 | 3 | 1 | 1 | 1 | 3 | 3 |
| 7 | 217 | 20 | 3 | 3 | 3 | 20 | 10 |
| 11 | 1,351 | 192 | 9 | 9 | 9 | 192 | 96 |
| 13 | 2,664 | 455 | 17 | 12 | 12 | 295 | 149 |
| 17 | 7,936 | 1,490 | 32 | 23 | 23 | 1,057 | 533 |
| 19 | 12,475 | 2,528 | 40 | 34 | 34 | 2,160 | 1,080 |
| 23 | 27,105 | 6,124 | 67 | 54 | 54 | 4,632 | 2,316 |
| **total** | **51,804** | **10,812** | **169** | **136** | **136** | **8,359** | **4,187** |

The table, and §3–§4 below, are a first run with the plain lattice search of note 2. With the faster lattice search of §6 all 169 polynomials are computed, and **all 169 agree**: 10,812 forms on 5,468 lattices; §6 has the complete figures.

**In all 136 polynomials of the first run the number of forms of each type A, B, C equals the lattice prediction.** [E] That is three numbers per polynomial, up to 88 + 88 + 88 forms on 132 lattices. The other 33 polynomials were not computed in that run, because their discriminant groups exceed the search limit; §6 computes them, and they agree.

So the dictionary of note 1 also holds, as a count, when Frobenius acts on NS by a non-trivial involution, and the lattice side separates the three types correctly.

**A sign bug in the lattice script was found and fixed** (§5). It changes no result of notes 2, 3 and 4.

---

## 2. The prediction [D]

**Lattices.** Let X₀ be a form of type π over F_p with ordinary irreducible polynomial P on T₆ = L(2), and E = Q(F). Frobenius is p·ρ(π) on N. Gluing N and T₆ forces F to act on L/2L, which is the discriminant group of T₆, by the class of ρ(π) in O(q_N) = S₆ × Z/2 (note 4 §4). In particular F² ≡ 1 on L/2L, so End(L) contains

> R′ = Z[F, V, (F² − 1)/2, (V² − 1)/2, (F − V)/2] and the unit-root part F/p,  V = p²/F.

R′ replaces the order R_min of note 2, which had (F − 1)/2. The class of F on L/2L is not a property of the order: it is computed for each lattice, as in note 4.

**Gluings.** For a lattice (L, F) with F in the class g = ρ(π), the F-equivariant gluings form a torsor under the centraliser Z(g) = Z_{S₆}(π) × ⟨z⟩ of g in O(q_N). Two gluings give the same form when they differ by a relabelling in Z_{S₆}(π) or by a root of unity of End(L). So

> #{forms of type π whose polynomial is P(t) or P(−t)} = Σ over lattices L of class (π, lines) of #( Z_{S₆}(π) \ Z(g) / μ(End L) ).

When E has only the roots of unity ±1 the summand is 2: the two forms are the configuration and its associate. A form is counted once; its two quadratic twists have polynomials P(t) and P(−t). For π = 1 this is the statement of note 2.

This counts pairs (surface, six-line structure). It does not need to know whether the form and its associate are the same surface over F_p.

**Surface side.** The forms are enumerated as in note 4. The polynomial needs S₁, S₂, S₃ over F_p, F_{p²}, F_{p³}. The product of the six lines is a product of linear and quadratic forms with coefficients in F_p, so for fixed (y, z) each factor is a translate of a fixed table in x, and the sum over F_{p³} never leaves F_{p³}. The sign of the functional equation, open when e₃ = 0, is fixed by requiring (F² − 1)/2 to be integral; exactly one sign is valid for every one of the 51,804 forms. [E]

---

## 3. What was measured [E]

1. **Per-type agreement in 136 of 136 polynomials.** Examples:

| p | polynomial (e₁, e₂, e₃) | index of R′ | lattices by class | forms A, B, C |
|---|---|---|---|---|
| 7 | (4, 7, −392) | 2 | 2 B | 0, 4, 0 |
| 11 | (12, 99, 0) | 8 | 3 B, 1 C | 0, 6, 2 |
| 11 | (10, 99, 484) | 16 | 8 A, 8 B | 16, 16, 0 |
| 23 | (8, −69, 4232) | 8 | 46 B, 46 C | 0, 92, 92 |
| 23 | (10, −345, −10580) | 4 | 44 A, 44 B, 44 C | 88, 88, 88 |

   Of the 136 polynomials, 39 have forms of one type only, 61 of two types and 36 of all three. The lattice side gives zero for the absent types in every case.
2. **Fourth roots of unity.** Four polynomials have Q(i) ⊂ E (p = 5, 13, 17, 17). There each lattice carries one form, not two, because i exchanges the two gluings; the double-coset formula gives the observed numbers (for example 1, 1, 1 forms on 3 lattices).
3. **Consistency with notes 3 and 4.** At p = 7 the three polynomials are those of note 4 §3, now with lattices computed over F₇ directly: 6, 2 and 2 lattices, equal to the numbers of fixed surfaces over F₄₉. This includes (4, 7, −392), whose lattices over F₄₉ could not be computed in note 3.
4. **Only "lines" classes occur.** All 4,187 lattices over R′ for these polynomials have a class (π, lines) with π of type A, B or C. None has F ≡ 1 on L/2L, and none has a class (π, swap). So for these polynomials every lattice is accounted for by a six-line form.

---

## 4. Limits

- **33 polynomials not computed in the first run** (2,453 forms; computed in §6), all with discriminant groups above 1,024. The index of R′ in the maximal order is up to 2,048 in these cases. One attempt at p = 13 with the limit raised to 4,096 did not finish the first such polynomial in 9 minutes. The subgroup search is the plain one of note 2.
- **Derived statements.** §2 is my derivation. It assumes that the nef cone is stable under the line permutations, so that every gluing gives a surface with a six-line structure, and it uses O(q_N) = S₆ × Z/2 as computed in note 4. The agreement supports it and does not prove it.
- **Only forms split over F_{p²}.** Frobenius of order 3, 4, 5 or 6 on the lines is not treated, nor are forms where Frobenius exchanges the two polarisations ("swap"). For those the order and the gluing count change.
- **Polynomials seen from the surface side only.** Lattices are computed for polynomials that occur among the forms. Polynomials with lattices of class A, B or C and no form would be missed; a Honda–Tate enumeration (note 2 §6 item 3) would close that.
- **Positivity.** Lattices are counted for one real place; independence of the place was checked on two polynomials only.
- **Still at p.** Nothing here concerns ℓ-neighbour steps or derived equivalences.

---

## 5. The sign bug in `sixlines_cm.py`

The lattice script decided the signs of α at the three real places by evaluating at 15 significant digits. The value of α at a place can be many orders of magnitude smaller than its coefficients (10⁻¹⁷ against 10¹⁸ was seen), and then the sign was noise. At p = 13, polynomial (8, −65, −2704), this produced 16 classes of mixed signatures (2, 4) and (4, 2) where there are 8 of signature (2, 4). It showed up because note 4's classification needs an isometry to q_N, which does not exist in signature (4, 2).

`sl_mask` now raises the precision until every value is far above the rounding error. After the fix:

- the outputs of note 2 (`sixlines_cm.jsonl`), note 3 (`sixlines_neighbours.jsonl`, `sixlines_q49_cm.jsonl`) and note 4 (`sixlines_types_q{25,49}.jsonl`), and their logs, are **byte-identical** to the earlier ones; [E]
- so no count or verdict in notes 2, 3 and 4 changes. The earlier runs were not affected, as far as their outputs show.

---

## 6. A faster search, and all 169 polynomials [E]

**The search.** The self-dual lattices between K and J were found by a plain search over subgroups of J/K, limited to groups of order 1,024. `sixlines_cm.py` uses `sl_lagr` instead:

- J/K splits into its ℓ-primary parts, and a self-dual lattice is a sum of one maximal isotropic R-submodule of each part;
- for each ℓ, the isotropic R-stable lattices H ⊇ K are explored from K. Every strictly larger one contains H + Rx for some x of order ℓ in H^⊥/H, so only the lines of that F_ℓ-space are tried;
- everything is exact lattice arithmetic in PARI (Hermite normal forms), with no enumeration of group elements.

The old search is kept behind `--slow`. On every polynomial computed before, in notes 2, 3, 4 and this one (40 + 31 + 136 = 207 polynomials), the new search gives the same counts, cycle types and classes; the output files differ only in the lines that were "not computed". Note 2's comparison takes 6 seconds, down from 66.

**All 169 polynomials agree.**

| p | polynomials | agree | forms | lattices | of which Q(i) ⊂ E |
|---|---|---|---|---|---|
| 5 | 1 | 1 | 3 | 3 | 1 |
| 7 | 3 | 3 | 20 | 10 | 0 |
| 11 | 9 | 9 | 192 | 96 | 0 |
| 13 | 17 | 17 | 455 | 242 | 6 |
| 17 | 32 | 32 | 1,490 | 773 | 9 |
| 19 | 40 | 40 | 2,528 | 1,268 | 0 |
| 23 | 67 | 67 | 6,124 | 3,076 | 0 |
| **total** | **169** | **169** | **10,812** | **5,468** | **16** |

The 33 new cases have index of R′ up to 2,048 and are the harder ones:

- **Unequal types.** For example p = 23, (2, 207, −2116), index 512: 138 lattices, with 72, 144 and 48 forms of types A, B, C, as predicted. In the first run the three types were equally many whenever all three occurred.
- **Sixteen polynomials with Q(i) ⊂ E**, against four before; the double-coset formula holds in each (for example p = 13, (2, 247, −676): 15 lattices carrying one split configuration and 6 + 10 + 6 forms).
- **Split lattices appear, and match note 2.** §3 item 4 is true of the first 136 polynomials only. In 9 of the new ones some lattices have F ≡ 1 on L/2L. These are exactly the 9 polynomials of note 2 for 13 ≤ p ≤ 23, and the same formula with π = 1 predicts the number of configurations with all lines rational: **42 predicted, 42 found in note 2's enumeration, polynomial by polynomial**. So one lattice enumeration over R′ gives the split count of note 2 and the three non-split counts together. For p = 23, (2, 207, −2116): 12 split configurations and 72 + 144 + 48 forms, on 6 + 36 + 72 + 24 lattices.
- No lattice has a "swap" class in any of the 169 polynomials.

**Also computed:** the four polynomials left open in note 2 (all agree: 7, 7, 19 and 10 surfaces, so 44 of 44) and the one left open over F₄₉ in notes 3 and 4 (14 surfaces and 14 lattices; conjugation and the 𝔭/𝔭̄ step both have 2 fixed points and 6 swapped pairs; both fixed surfaces have class (B, lines) on both sides).

**Limits that remain:** those of §4 other than the first. The search has a cap of 200,000 intermediate lattices, not reached here.

---

## Reproduce

From `experiments/`, with `python -I`:

```bash
python -I sixlines_nonsplit.py selftest
```

```bash
for p in 11 13 17 19 23; do python -I sixlines_nonsplit.py enumerate --p $p --out sixlines_nonsplit_p$p.jsonl --procs 12 2> sixlines_nonsplit_p$p.log; done
```

```bash
python -I sixlines_nonsplit.py lattices --list --out sixlines_nonsplit_lat_small.jsonl sixlines_forms_p5.jsonl sixlines_forms_p7.jsonl sixlines_nonsplit_p11.jsonl sixlines_nonsplit_p13.jsonl > sixlines_nonsplit_lat_small.log
```

```bash
for p in 17 19 23; do python -I sixlines_nonsplit.py lattices --list --out sixlines_nonsplit_lat_p$p.jsonl sixlines_nonsplit_p$p.jsonl > sixlines_nonsplit_lat_p$p.log; done
```

The self-test compares the fast sums with the stored values of note 4 for p = 5 and 7 (273 forms). Enumeration on 12 processes: seconds up to p = 13, 1 minute at 17, 2 at 19, 11 at 23. Lattices, with the search of §6 and all polynomials: under half a minute up to p = 13, then 1, 2 and 8 minutes. Add `--slow` to `sixlines_cm.py` for the original search; `sixlines_nonsplit.py lattices` always uses the new one.

## Sources

As in notes 1, 2 and 4.
