# Notes

The detailed arguments, computations and checks behind [SUMMARY.md](../../SUMMARY.md). Each note states what it builds on, labels every statement by its evidence, says what it does not show, and ends with the commands that reproduce its numbers (from `experiments/`, with `python -I`) and its sources.

Evidence labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived in this work, with the argument given; **[E]** measured by the scripts in this repository (`experiments/`); **[H]** heuristic; **[R]** re-derived in a separate re-derivation check, with its own code (notes 18 and 25; notes 26 to 29 state the findings of their own check, whose code is not part of this repository).

**Reading order.** For the hardness question alone: notes 6, 13, 14, 20, 24, 27, then 17. Notes 7 to 12 (ring-changing steps) are lattice theory that the verdict does not need.

| Note | Topic | Scripts |
|---|---|---|
| **Torsors and the six-line family** | | |
| [1](01-prior-art.md) | What the literature already gives: the class-group action is an assembly of Taelman, Valloni, Bayer-Fluckiger and Bragg–Yang | |
| [2](02-six-line-count.md) | Every double plane branched along six rational lines over F_p, 7 ≤ p ≤ 37, counted against unimodular Hermitian lattices: all counts agree | `sixlines.py`, `sixlines_cm.py` |
| [3](03-torsor-action.md) | The action on the torsor: ℓ-steps on lattices; the step at p (Frobenius conjugation) on surfaces over F₂₅ and F₄₉ | `sixlines_fq.py`, `sixlines_cm.py --neighbours` |
| [4](04-forms-over-the-prime-field.md) | Forms over the prime field; the Frobenius action on Néron–Severi read off lattices and surfaces | `sixlines_forms.py`, `sixlines_types.py` |
| [5](05-non-split-count.md) | Forms whose lines are permuted by an involution, p ≤ 23: 169 polynomials, all counts as predicted | `sixlines_nonsplit.py` |
| [6](06-torsors-under-cm-class-groups.md) | **Ordinary K3 surfaces over F_q as torsors under CM class groups** (Propositions A–H), with the evidence of notes 2–5 and the Artin–Tate check | `sixlines_rational.py` |
| **Steps that change the endomorphism ring** | | |
| [7](07-ring-changing-steps.md) | All rational steps at odd ℓ ≠ p are Kneser neighbours; volcanoes and descents | `sixlines_vertical.py` |
| [8](08-below-the-top.md) | Below the top level at eigenvalue ±q: ℓ-volcanoes | `sixlines_levels.py` |
| [9](09-ramified-over-e0.md) | Eigenvalue ±q at primes ramified over E₀ | `sixlines_ramified.py` |
| [10](10-fourth-power.md) | Eigenvalue ±q when ℓ⁴ divides P(±q) | `sixlines_deep.py` |
| [11](11-rank-six-fourth-power.md) | Multiplicity 6 with ℓ⁴: lattices over an explicit quartic order | `sixlines_rank6.py` |
| [12](12-rank-six-galois-action.md) | Multiplicity 6 with ℓ⁴, the non-square case: self-dual hermitian lattices | `sixlines_hermitian.py` |
| **Hardness and the reduction to fourfolds** | | |
| [13](13-hardness-on-paper.md) | **The hardness question on paper**: the group G and \|G\| ≈ 0.01·q², Proposition Y, attacks and costs, the decision problem | `sixlines_groupsize.py` |
| [14](14-fourfold-of-a-six-line-surface.md) | **The abelian fourfold built over F_p** as the Prym variety of a genus-5 curve (Paranjape's construction made explicit) | `sixlines_prym.py` |
| [15](15-fourfold-map-and-the-group-action.md) | Does the surface-to-fourfold map respect the action? Checked on class groups, not on steps | `sixlines_dictionary.py` |
| [16](16-fourfold-over-f25-and-f49.md) | The fourfold construction over F₂₅ and F₄₉ | `sixlines_prym_fq.py` |
| [18](18-checks-of-notes-13-to-16.md) | Re-derivation checks of notes 13 to 16 | `sixlines_review.py` |
| [21](21-fourfold-forms-over-fp.md) | Which fourfold over F_p belongs to a surface of sign − or with reducible P | `sixlines_forms_fp.py` |
| **Computing the action, and what leaks** | | |
| [19](19-quadratic-characters-from-a-surface.md) | Can a quadratic character of G be computed from a surface? χ₋₁ singled out | `sixlines_characters.py` |
| [20](20-cm-method-and-the-action.md) | **The CM method**: the action evaluated on small groups; ℓ-steps equal Frobenius in characteristic ℓ | `sixlines_cmmethod.py` |
| [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) | Arithmetic in the Jacobian of the genus-5 curve; the level test; the sign δ from a Tate pairing | `sixlines_jacobian.py` |
| [24](24-the-form-mod-4-is-the-character.md) | **δ is the character χ₋₁** on three torsors: one bit of the class leaks | `sixlines_delta_labels.py` |
| [25](25-checks-of-notes-20-23-24.md) | Re-derivation checks of notes 20, 23 and 24 | `sixlines_review2.py` |
| [26](26-delta-is-the-character-at-q-41.md) | Out-of-sample test at q = 41: δ is χ₋₁ on two more torsors (60 of 60 in all) | `sixlines_delta41.py` |
| [28](28-leakage-audit.md) | **Leakage audit** of about 11,400 candidate invariants against the CM labels: nothing beyond δ | `sixlines_leakage.py` |
| [29](29-characters-of-order-3.md) | Characters of order 3 of G: none is read on the two torsors where a cubic symbol exists, and no candidate separates the cosets on any of the seven with 3 ∣ \|G\| | `sixlines_order3.py` |
| **Below Picard number 16** | | |
| [17](17-below-picard-16.md) | No reduction known and no computable step; a reading of the literature | |
| [22](22-two-neighbours-of-elliptic-k3.md) | A step at 2 on elliptic K3 surfaces, tested on 21,240 double covers | `ellk3_two_neighbours.py` |
| **Synthesis** | | |
| [27](27-hypothesis-ledger.md) | **The hypothesis ledger**: 47 hypotheses classed; 170 quoted numbers recomputed from the stored files | `sixlines_ledger.py` |

Propositions and theorems carry letters (A to U, and Y in note 13), and lemmas carry numbers or letters, assigned across the notes in the order the results were introduced; some letters and numbers were never used, and notes 14, 19, 21 and 22 label their lemmas locally. The notes are numbered in the order they were written, so a note may refer forward to a later one that settles a question it leaves open.
