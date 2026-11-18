# 16. The fourfold over F₂₅ and F₄₉

*Extends the comparison of [note 14](14-fourfold-of-a-six-line-surface.md) from prime fields to the two fields of [note 3](03-torsor-action.md), where most surfaces are not defined over the prime field. Builds on notes 2, 3, 14 and 15. Labels: **[D]** derived here; **[E]** measured by `experiments/sixlines_prym_fq.py`. Checked in [note 18](18-checks-of-notes-13-to-16.md), finding 10 (no finding; the two long runs were not rerun).*

---

## 1. Result

1. **The construction of note 14 works over F_q, q = p² [D, E].** Nothing changes in it except that the six lines have coefficients in F_q. Since 25 and 49 are 1 mod 4, i is in F_q, and only the case "c₄ rational" of note 14, Lemma 4 can occur with Frobenius commuting with J.
2. **Every surface with irreducible P gives its polynomial [E].** Over F₂₅: 37 configurations, all 4 polynomials of note 3, 858 ordered pairs of lines. Over F₄₉: 1,108 configurations, all 28 polynomials of note 3, 3,324 pairs. The 37 and 1,108 configurations make up 21 and 564 surfaces, a configuration and its associate being one surface (recounted in [note 27](27-hypothesis-ledger.md) §1 from `sixlines_prym_q25.jsonl` and `sixlines_prym_q49.jsonl` with the association of `sixlines_q25.jsonl` and `sixlines_q49.jsonl`). In every pair the second exterior power of the quartic χ, computed from the points of the curve over F_q, F_{q²}, F_{q³}, is P or P(−t), and the coefficient of X in χ is q times the conjugate of that of X³.
3. **Most of these configurations are not defined over the prime field [E].** 30 of the 37 configurations and 1,060 of the 1,108 are not isomorphic to their Frobenius conjugate. So the comparison no longer rests on surfaces over F_p only, and it now covers the torsors on which note 3 tested the action of the class of 𝔭/𝔭̄.
4. **Reducible P over F₂₅ [E].** All 365 configurations over F₂₅ were run. For the 131 configurations with sign + and reducible P that have a computable pair, all 2,829 pairs give P or P(−t).

This removes the limit "prime fields only" of note 14 §7 for these two fields. It does not test a step: the map commutes with Frobenius conjugation by construction (note 15 §5.1).

---

## 2. What changes over F_q [D]

- **Square roots.** Note 14 fixed signs by the smaller representative in F_p. Here every square root is taken once, in F_q, and carried into the larger fields by an embedding. A point of E that is conjugate to its negative is written (X, ω·r) with r in F_q and ω a square root of a fixed non-square n₀ of F_q; changing the sign of ω changes both such points at once, which only translates the branch points by a point of order 2 (note 14, Lemma 2).
- **The auxiliary point Q** of note 14, Lemma 3 lies in E(F_{q²}) and satisfies Q − Frob_q(Q) = a + b. It is found by search over E(F_{q²}) and fixed by the least coordinates.
- **Fields.** The traces need F_q, F_{q²}, F_{q³}, and F_{q⁶} when two of a, b, c are conjugate to their negatives; for q = 49 that is F_{7¹²}.
- **Types.** With q ≡ 1 mod 4, Frobenius commutes with J exactly when c₄ is rational. The other case is recorded and not computed.

Over a prime field with p ≡ 1 mod 4 the script follows the same steps as `sixlines_prym.py` with different sign conventions. Run on p = 13 and 17 it gives the counts of note 14: 140 and 618 computed pairs, of which 63 and 312 give P and 77 and 306 give P(−t).

---

## 3. Results [E]

| | F₂₅ | F₄₉ |
|---|---|---|
| configurations run | 365 (all) | 1,108 (sign +, P irreducible) |
| ordered pairs tried per configuration | 30 | until 3 are computed |
| sign +, P irreducible | 37 | 1,108 |
| of these ordinary | 37 | 1,080 |
| with a computed pair | 37 | 1,108 |
| not isomorphic to the Frobenius conjugate | 30 | 1,060 |
| pairs computed | 858 | 3,324 |
| pairs giving P or P(−t) | 858 | 3,324 |
| pairs where Frobenius inverts J | 0 | 0 |
| ordinary polynomials up to the twist | 4 | 28 |
| sign +, P reducible: configurations with a computed pair, pairs, exceptions | 131, 2,829, 0 | not run |
| sign − | 168, not compared | not run |

- The 1,080 ordinary configurations over F₄₉ are the classes of note 3 §6, which make up its 540 surfaces. The other 28 have an irreducible polynomial that is not ordinary; they pass as well.
- In all 7,011 computed pairs of surfaces with sign +, c₃ = q·c̄₁.
- Over F₂₅, 29 configurations with sign + and reducible P have no computable pair, and 336 pairs of such configurations are of the type where Frobenius inverts J. Both phenomena are as in note 14 §5.2.
- 2,503 pairs over F₂₅ and 378 over F₄₉ were skipped as degenerate.

---

## 4. Check and limits

- **Regression.** The counts on p = 13, 17 quoted in §2 equal those of note 14, with a different implementation of the signs.
- **A correction to the classification in `sixlines_prym.py`.** It called P irreducible when the degrees of its factors, counted with multiplicity as one number each, were (6). That would also accept a sixth power such as (t − q)⁶, which occurs over F₂₅ (configuration (2, 3, 3, 4)). No polynomial of note 14 was affected: its 534 configurations are exactly those of the 44 polynomials of `sixlines_cm.jsonl`, which I rechecked. Both scripts now list the degrees with multiplicity.
- For F₄₉ only 3 pairs per configuration and only configurations with irreducible P; a full run would take about ten times longer.
- Internal checks as in note 14 §7: the scripts stop if a character value is not a fourth root of unity; none did.
- The surfaces over F₂₅ and F₄₉ with sign − are not compared, as in note 14.
- Fields F_{p²} with p = 5, 7 only.

---

## 5. What this does not give

- A test of the group action at a prime other than p. Frobenius conjugation is the step at p, and the map commutes with it by construction.
- Anything beyond Frobenius polynomials: the isogeny class of the fourfold with its Z[i]-structure, as in note 14.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_prym_fq.py --workers 14 --select all --pairs 0 --run 25 --out sixlines_prym_q25.jsonl > sixlines_prym_q25.log
python -I sixlines_prym_fq.py --workers 14 --select irreducible --pairs 3 --run 49 --out sixlines_prym_q49.jsonl > sixlines_prym_q49.log
python -I sixlines_prym_fq.py --workers 12 --select all --pairs 0 --run 13 17
```

About 10 minutes, 44 minutes and 20 seconds with these numbers of processes. They need `sixlines_q25.jsonl` and `sixlines_q49.jsonl` of note 3 and, for the third, the files of note 2.

## Sources

- In this repository: notes [2](02-six-line-count.md), [3](03-torsor-action.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md).
