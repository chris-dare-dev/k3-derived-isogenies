# 25. Re-derivation checks of notes 20, 23 and 24

*A re-derivation check of notes [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) and [24](24-the-form-mod-4-is-the-character.md), which together claim the first leak of information about the class in the group action on six-line surfaces: one bit, the character χ₋₁, for the polynomials where it exists. It reads the Python side of the three scripts, the two sources of note 20 and the paper of Castryck, Sotáková and Vercauteren, reruns pieces of all three scripts, and adds six checks (`experiments/sixlines_review2.py`). Not re-checked.*

*Labels, all for work done in this check: **[R]** a statement of the checked notes, re-derived here. **[E]** measured here, by `sixlines_review2.py` or by a rerun of the script of a checked note where that is stated. **[P]** read here in the cited source. **[D]** derived here and not in the checked notes; not independently verified.*

---

## 1. Outcome

**The claim of note 24 §1 does not fail.** I found no error in it, no shared input that could produce the agreement, and no simpler quantity that reproduces it. The split of the surfaces into two orbits was obtained again by two routes that do not use the action of ideals on lattices, and δ was recomputed on nine curves with the delicate arithmetic redone. What I did find is one wrong reference (to CSIDH), two gaps in the argument of note 24 §2 and §3 that are closed or narrowed here, and that the level test of note 23 §4 carries less weight than its counts suggest.

| Question | Answer |
|---|---|
| Does the claim of note 24 §1 (δ is χ₋₁ on the three torsors) fail anywhere? | **No.** Reproduced; reached by two further routes (§4.3, §4.4) [E] |
| Could the agreement have a cause other than the one claimed? | None found. The two computations share the list of configurations and nothing else [R]. Thirteen simpler candidates built from the stored data do not separate the orbits [E]. The cause is now visible: δ is, up to one sign for each polynomial, the form modulo 4 of the lattice that note 20 attaches to the surface (§4.4) [D, E] |
| Note 24 §2, from two orbits of a subgroup of index 2 to all of G | The logic is valid. One step, "the group that acts is G", rested on a count that does not suffice for Z/4 × Z/2. Supplied here: the class polynomial is irreducible over E, and Frobenius at 160 primes with ℓ ≡ 3 mod 4 has cycles of the right length (finding 2) |
| Note 24 §2, the subgroup is ker χ₋₁ | **Correct** [E]: in PARI's group the classes of the primes of degree 1 with ℓ ≡ 1 mod 4 are a subgroup of index 2 and those with ℓ ≡ 3 mod 4 are its other coset, in all three fields |
| Note 24 §3, the indicator polynomial and its reduction | **Correct** [R]. The code does what the text says. The rerun of one field gives the stored record exactly. One identification is thinner than the text says (finding 3) |
| Note 23 §5, the pairing and the kernel step | The code computes what the text says [R]. δ recomputed on nine curves, on a twist not used before, with the 2-adic arithmetic written again: 9 of 9. Seven identities of the pairing are now logged. One sentence is imprecise (finding 7). The GP code of the group law was not read |
| Note 20, formulas | **Agree with the two sources** [P]: the theta constants, the table of ten characteristics, the ten products, the fifteen relations, the six x_k, J₂ to J₆ and the misprint. The sixteen theta constants were recomputed by a direct sum [E] |
| Note 20, the order R, even lattices, the orientation, the Frobenius test | Reruns give the stated numbers; every total of §§6 and 7 was recomputed from the stored file. The orientation is fitted as note 20 says, and the result of note 24 does not depend on it [D]. The completeness of the enumeration of §5.2 was not checked |
| Note 23 §4, the level test | The prediction is forced by the lattice computation and the rule Λ², so it is not fitted [R]. The counts are right. But five of the seven surfaces with two blocks are the self-associated ones, so the count is informative in one polynomial only (finding 10) |
| Note 24 §5, decisional Diffie–Hellman and cost | The statement is correct as worded for the decision problem. The reference to CSIDH is not: the paper says its attack does not affect CSIDH (finding 5). The cost paragraph is fair and labelled as heuristic |

The standing of the claim after this check: on the three torsors it is as solid as a computation of this size can make it, and it now has a mechanism on the lattice side. Beyond them it is the conjecture that note 24 says it is. Of the 80 polynomials of note 19 with the character, 53 have q ≡ 1 mod 8 (note 19 §2.5; recounted from `sixlines_characters_genus.jsonl` [E]), and both primes tested are 5 mod 8.

---

## 2. What was checked, and how

| Note | Statement | What I did | Result |
|---|---|---|---|
| 24 | §1 items 1, 2, table of §4 | Reran `run 37 18 1 1 --digits 1200 --primes 8 --workers 8` (90 seconds) | the record equals the stored record, field by field [E] |
| 24 | "They share the list of surfaces and nothing else" | Read `sixlines_delta_labels.py` against the text; listed the files each side reads | holds [R], §4.1 |
| 24 | δ could be a function of simpler data | `stored`: 13 candidates against the stored orbits, by record and by surface | none separates the orbits [E], §4.2 |
| 24 | the two orbits | `galois`: H factored over E and E(i) by PARI, one factor reduced at the primes above 𝔭 | same split in all three fields [E], §4.3 |
| 24 | the two orbits, again | `galois`: the form modulo 4 read on the lattices | same split in all three fields [D, E], §4.4 |
| 24 | §2, χ₋₁ of a step is (−1/ℓ) and the subgroup is its kernel | `galois`: classes in G of 120 primes with ℓ ≡ 1 mod 4 and 160 with ℓ ≡ 3 mod 4 | holds [E] |
| 24 | §2, the inference to all of G | Went through the argument | valid; finding 2 |
| 24 | §3, I_T is i times a polynomial over E; reduction | Re-derived; read the code | holds [R]; finding 3 |
| 24 | §1 item 3, 1 in 7,470,540 | Recomputed: 462 · 462 · 35 | holds [R] |
| 24 | §5 | Read the paper of Castryck, Sotáková and Vercauteren, abstract, §5.3 and the conclusion | findings 5, 6 |
| 23 | §5.3, the formula for Z | Re-derived from the pairing identity; read `character_pair` line by line | holds [R] |
| 23 | §5.4 counts | Recounted from the two stored files | 242, 106, 126, 9, 1; M and u as stated; u constant on 10 of 32 [E] |
| 23 | δ | `delta`: nine curves, twist 3, another seed, own 2-adic arithmetic; and `--one 29 2 3 5 20 1 2` of the script of note 23 | 9 of 9 agree; the stored record of note 23 is reproduced exactly [E] |
| 23 | §3, the pairing | `pairing`: seven identities on one curve over F_{29⁴} | all hold [E] |
| 23 | §5.3, values of M over the kernel | Recounted from the stored `kernel` fields | finding 7 |
| 23 | §4.1, the prediction | Re-derived the Jordan types of Λ² for a double root; reran `--lattice` | holds [R]; same output |
| 23 | §4.3 | `levels`: recount of the stored file, with self-association | counts as stated; finding 10 |
| 23 | 1.4 %, 2⁻⁷⁴ | Recomputed | hold [R] |
| 20 | §2 | Read Matsumoto–Terasoma (19), Table 1, Remarks 1 and 2, Propositions 6 and 7, formula (4); Clingher–Malmendier–Shaska (2.1), (2.2), (2.12) to (2.16), Theorem 2.8, (2.30) to (2.34), Lemma 2.17, Corollary 2.25; compared with the text and with `cm_tsq`, `cm_J`, `cm_tmat`, `cm_rel` | agree [P]; finding 11 |
| 20 | theta constants | `theta`: direct fourfold sum for the sixteen characteristics at three points | agree to 10⁻³⁶ [E] |
| 20 | §4 | Reran `selftest --digits 60` | output identical to the stored log |
| 20 | §5.3 | Reran `points 17 2 1 1 --digits 30 --maximal` and `--nolevel` | 0 even and 1 odd per type; 512, 64, 597, 78 even per type [E] |
| 20 | §§6, 7 for two fields | Reran `sweep 17 23 --select 17,2 23,10` | same counts, including 8 of 8 inverse permutations failing [E] |
| 20 | §§6, 7 totals | Recomputed from `sixlines_cmmethod.jsonl`, last record of each field | 74, 165, 378, 0, 69, 390, 346, 32, 136, 1,156, 630, 610, 594, 67, 389: all as stated [E] |
| 20 | §3, the orientation | Went through what the other sign does | finding 12 |
| 20 | §7.1 | Checked that the type norm of π is F² and the norm of the ideal 𝔟_jk is ℓ⁴ | consistent [R]; the reflex norm itself not re-derived |

Not read: the GP code of the group law, the Miller functions and the series arithmetic in `sixlines_jacobian.py` (its first 520 lines), and the GP and Python code that enumerates the lattices of R in `sixlines_cmmethod.py`. Their support is the checks that come out.

---

## 3. Findings

Each is classed as **error**, **gap** (something asserted with less support than the text gives it), **imprecise**, or **confirmation**.

### 3.1 Finding 1. Note 24 §1, items 1 and 2: δ is χ₋₁ on the three torsors — confirmation

Reproduced, and reached by two other routes (§4.3 and §4.4). The rerun of q = 37, e₁ = 18 gives the stored record. The split of the 8, 12 and 12 surfaces is the same whether it comes from the orbits of the 48 ideals (note 24), from the factorisation of H over E(i), from the form modulo 4 of the lattices, or from δ.

### 3.2 Finding 2. Note 24 §2: "the group acting is G and not a smaller group with two orbits" — gap, closed here

The text says this is covered "by note 20's count (the lattices of each type number |G|) and by the order of the subgroup found, |G|/2". For G = Z/12 that suffices if G acts at all: an element outside the subgroup that preserved the orbits would give a retraction of Z/12 onto its subgroup of order 6, and there is none. For G = Z/4 × Z/2 it does not suffice [D]: the subgroup found is generated by (1, 1), the element (0, 1) could act trivially, and the picture would be the same, eight lattices, a subgroup of order 4 acting freely with two orbits of four.

What excludes it [E]: H is irreducible over E in all three fields (one factor of degree 8, 12, 12), so the Galois group is transitive on the roots; and at 160 primes of degree 1 of E with ℓ ≡ 3 mod 4 (44, 70 and 46, all those with ℓ below about 600 where H is integral and squarefree) every factor of H modulo the prime has degree equal to the order in G of the class of the prime, as at the 120 primes with ℓ ≡ 1 mod 4. A transitive abelian group on |G| points acts simply transitively. And since H = H₁·H₂ over E(i) with the two orbits as the roots of H₁ and H₂, an automorphism that moves i moves one orbit to the other: Frobenius at a prime of degree 1 with ℓ ≡ 3 mod 4 exchanges the orbits. That is the step outside the kernel that note 24 §1 item 7 and §6 say was not computed. It is a statement about the Galois action; that Frobenius at 𝔩 is the step of the class of 𝔩 is the theorem of Shimura and Taniyama, recalled in note 20 and tested there for ℓ ≡ 1 mod 4 only.

### 3.3 Finding 3. Note 24 §3: "the Galois action on the roots, which note 20 §7.3 identifies with the action of G" — gap, narrowed

Note 20 §7.3 was run on 15 fields. For q = 37, e₁ = 18 it tested 10 cases and skipped 26; the two torsors of 12 surfaces were not among the 15, and `sixlines_delta_labels.py` has no Frobenius test. So on the three torsors of note 24 the identification rests on 10 cases in one of them and on the other fields. The cycle lengths of finding 2 are a weaker test (the cycle type, not the permutation) but cover all three torsors and both cosets: 280 primes, no exception [E]. The comparison of permutations of note 20 §7.3 on the two torsors of 12 is not done.

### 3.4 Finding 4. Note 24 §1 item 4 and note 23 §5.2: what the agreement is evidence of — confirmation, sharpened

Note 23 §5.2 derives δ from the hermitian form on the line V₁ of the 2-adic Tate module and calls its lattice statements sketches. They can be tested on the lattices of note 20 directly, which had not been done (§4.4). For an even lattice (L, ξ) the form on V₁ = K_{w₁} is h(m) = 2·N(m)·t with t the image of −iξ in Q₂, so note 23's u is the unit part of t modulo 4 [D]. Computed for the 96 lattices of the three fields [E]: it takes each value on half of the lattices of each type; its indicator polynomial is i times a polynomial over E, recognised in exactly one of the eight ways and with no sign between the types; and modulo 𝔭 it splits the surfaces exactly as δ does, 32 of 32.

So the content of the agreement of note 24 is: **the number that note 23 computes in characteristic p from the curve of note 14 is, up to one sign for each polynomial, the class modulo 4 at w₁ of the polarisation of the lattice of the CM lift.** On the lattice side the character property is immediate: a step with multiplier ℓ sends (L, ξ) to (𝔟⁻¹L, ℓξ) and multiplies t by ℓ, for ℓ ≡ 3 mod 4 as well as ℓ ≡ 1. Note 24 §1 item 4 says the agreement is a test of the Prym variety of note 14; this says of what: its polarised line V₁ has the class of the lattice of Lombardo–Peters–Schütt, on every surface and every usable pair of lines.

### 3.5 Finding 5. Note 24 §1 item 6 and the last paragraph of §1; the remarks in note 13 §9 and note 19 §1: which attack on the decision problem carries over — imprecise

Read in the paper [P]: the abstract says the attack "works very efficiently for all supersingular elliptic curves over F_p with p ≡ 1 mod 4"; §5.3 says the simplest countermeasure is a setting where the class group has no 2-torsion, "e.g. supersingular elliptic curves over F_p with p ≡ 3 mod 4. This corresponds precisely to the CSIDH setting, so our attack does not impact CSIDH"; the conclusion repeats it. So the attack is on class group actions with a class group of even order (the scheme of Couveignes, Rostovtsev and Stolbunov, and supersingular curves with p ≡ 1 mod 4), and CSIDH is the case it does not touch. So what carries over is the attack of Castryck, Sotáková and Vercauteren on the decision problem for class group actions, not an attack on CSIDH's decision problem, and the remedies are the countermeasures of their §5.3 (a group without the character, or cofactors, which is what steps in the kernel only amounts to); the notes now say this. In the same way, supersingular curves with p ≡ 1 mod 4 are the paper's case and not a form of CSIDH, and note 24 §5 and note 23 §5.6 now treat them so. The mathematics of note 24 is not affected. One consequence of reading the paper's wording: its first countermeasure is a group with no character of order 2 at all, which is more than "a polynomial without the character". Note 19 finds a character of order 2 on 477 of the 480 groups of its list; for 397 of them none has a known description on a surface, which is a statement about what is known, not a countermeasure.

### 3.6 Finding 6. Note 24 §5: the statement and the cost — confirmation, with one omission

[R] If δ(g·X) = χ(g)·δ(X) then δ(X)·δ(abX) = δ(aX)·δ(bX), and a surface drawn at random satisfies it with probability 1/2 when χ is not trivial: an advantage of 1/2, as stated. Restricting a and b to the kernel removes it, as stated. The cost paragraph claims no more than a heuristic. My runs take 30 to 35 seconds for one curve on a machine that is not loaded [E], against the 75 it gives. The omission: the two primes tested are 5 mod 8, where the pairing has order 16; for q ≡ 1 mod 8 the order is 32 or more and the method has not been run, and 53 of the 80 polynomials are of that kind [E]. Note 24 §6 now gives the size of this limit.

### 3.7 Finding 7. Note 23 §5.3: the values of M over the kernel — imprecise

In the stored records the values of M over the kernel include M = 2, 6, 10 and 14 (twice an odd number, for elements of order 16), as §5.4 of note 23 says, and note 23 §5.3 now says so too. What holds, and is what the argument needs [E]: in each of the 142 records, all kernel vectors whose M is not divisible by 8 have the same odd part modulo 4. M = 8 says nothing, because M is known modulo 16. In one record 9 vectors give no M, as note 23 says.

### 3.8 Finding 8. Note 23 §5: δ — confirmation

[R] From T_N(X, Y) = ζ^{E(N·X̃, ϖ·Ỹ)} and X̃ = γm/ϖ₁, the exponents of T(X, X) and T(X, iX) are M·Im and M·Re of conj(N/ϖ₁), so Z = M·conj(N/ϖ₁) as stated; `zof` assembles exactly these two pairings from the matrices. The point X is verified by the group law, not by pairings, so the linear algebra is only a search.

[E] `delta`: three curves of each torsor (two surfaces of one orbit and one of the other), each on the twist C₃, which note 23 did not use, with seed 2026. The root of χ in Z₂[i] is found by lifting one binary digit at a time, not by `factorpadic`; N/ϖ₁ and M are computed from the two pairings of X evaluated directly; X is the last suitable kernel vector, not the first. Result: 9 of 9 computed, π X = π₁X, J²X = −X and the order hold on each, δ equals the stored value of the surface on each, and the other convention for Z has no solution on any. The valuation of χ′(π₁) at 1 + i is 6 on all nine, as §5.3 says, and that of ϖ₁ is 7, as in every stored record. The routine of note 23 itself, on one stored curve, returns the stored record including the kernel tallies.

What this does not make independent: the group law and Miller's algorithm are those of the script of note 23.

### 3.9 Finding 9. Note 23 §3: the pairing identities "checked by hand … not logged" — gap, closed for one curve

`pairing 29 2 3 12 27 1 3` logs, over F_{29⁴} with N = 16: additivity in each argument, T(3z, y) = T(z, y)³, T(z, 16R) = 1 and T(z, y + 16R) = T(z, y) for a random class R, T(Jz, Jy) = T(z, y), and T(πz, πy) = T(z, y)^p. All hold, on values that are not trivial (logarithms 13, 6, 4, 13).

### 3.10 Finding 10. Note 23 §4: the level test — imprecise as to its weight

[R] The prediction is not fitted. For a double root a of χ and simple roots b, c, the factor of P at ab is the tensor product of the block at a with a line, so F has kernel dimensions [1, 2] when π has one block and [2, 2] when it has two, and the lattice computation says which level is which. `--lattice` rerun gives the stored output.

[E] `levels`: the counts are as stated, 1 and 6, 5 and 14, 1 and 2, on 422 pairs with no disagreement inside a surface. But the surfaces with two blocks are, by polynomial, 1, 5 and 1, of which 1, 3 and 1 are self-associated; no self-associated surface has one block. On the lattice side the classes with 4 automorphisms number 1, 3 and 1 and are all on the upper level. Note 23 notes this as "one finer agreement" and guesses that a self-associated configuration is a surface with an extra automorphism; association is the transposition of the period point (Matsumoto–Terasoma, Proposition 6 (4)), which supports the guess, and I did not check it further. If it is right, the matter is the other way round: in the first and third polynomial the agreement of counts is the agreement "one object with an extra automorphism, and it is on the upper level", which could hold for reasons that have nothing to do with the rule Λ². What is left that is not explained that way is one number: in q = 37, e₁ = 34, exactly 2 of the 16 surfaces that are not self-associated have two blocks, and the lattice side has 2 classes with 2 automorphisms on the upper level. The solid part of §4 is the other statement, that every pair of lines of a surface gives the same number of blocks (422 pairs). The sentence "it is passed in every case that could be computed" is true; "the first agreement between the two sides on an invariant that is not an isogeny invariant" is better carried by note 24.

### 3.11 Finding 11. Note 20 §2: the formulas — confirmation

[P] The definition (19) of Θ_{a,b}, Table 1, Remark 2, Proposition 7 (1) and (2), and that association is the transposition (Proposition 6 (4), formula (16)) are in Matsumoto–Terasoma as quoted. The ten products (2.1), the fifteen relations (2.2), the x_k of (2.34) and J₂ to J₆ are in Clingher–Malmendier–Shaska as coded, and the misprint is real: (2.32) has 20s₅, (2.31) has 12s₅, and (2.33) gives s₅ = 60J₂J₃ − 20J₅, which is 12. The order of the ten characteristics in `cm_tsq` is Table 1 carried to the numbering of (2.1); I typed it again from the two sources and compared.

[E] `theta`: the sixteen constants by a direct sum over Z[i]² at three points of D, against `cm_theta16`: largest difference 3·10⁻³⁷ at 38 digits; the six characteristics outside the table vanish; the ten squares against `cm_tsq`, which sums over a quarter plane: 9·10⁻³⁷.

### 3.12 Finding 12. Note 20 §3: the fitted orientation — confirmation, and it does not enter note 24

Note 20 states the fit correctly: one field, then 128 of the 136 cases of the column "inverse fails" on other fields. [D] With the other sign every period point is replaced by its complex conjugate, the class polynomial by its image under the complex conjugation of E, and a labelling of the torsor by its composite with x ↦ g₀ − x. That map sends a coset of a subgroup to a coset. So the split into two orbits, which is all that note 24 uses, is the same for either orientation; so are the routes of §4.3 and §4.4.

### 3.13 Finding 13. Note 20 §§5 to 7: counts — confirmation, with one thing not checked

All totals recomputed from the stored file agree with the text, and the reruns of two fields and of the three lattice counts for q = 17 agree. For q = 37, e₁ = 18 note 20 says the 36 ideals give the kernel of the character "ℓ mod 4"; confirmed in PARI's group for ℓ of both classes [E]. Not checked: that the list of lattices of R in §5.2 is complete (note 20's own argument; its support is still that the count equals the number of surfaces, now in 17 fields), and the count of 234 odd lattices of the maximal order, which was a loop that is not a command.

### 3.14 Finding 14. Statements that cannot be checked — gap, minor

- Note 24 §4: "with three primes per field … the groups and orbits are the same (the first runs, not kept)".
- Note 23 §6: "the formula for δ was written down before the first comparison of two pairs".
- `sixlines_jacobian.py` and `sixlines_delta_labels.py` were last modified two and three minutes after their last stored output. Both reproduce a stored record exactly (findings 1 and 8), so the stored output is that of the current scripts as far as tested.

The second does not matter for note 24: the factor sign(Re c₁·Im c₁)·sign(c₂) depends only on the quartic of the curve, and the only thing it could have been fitted to is agreement between pairs of one surface, which uses no information about orbits. `stored` confirms that without it, or with half of it, nothing separates the orbits.

---

## 4. Could the agreement of note 24 have another cause?

### 4.1 Shared input [R]

`sixlines_jacobian.py --character` reads `sixlines_p{p}.jsonl` (configurations, their associates, the sums S) and `sixlines_prym.jsonl` (pairs of lines and the quartic of each curve), and computes in characteristic p. `sixlines_delta_labels.py` reads `sixlines_cm.jsonl` for the configurations of the polynomial, computes lattices, period points, H, the permutations of 48 ideals and the indicator Γ, and names the two orbits by the smaller value of Γ/H′. Only then does it open the two files of note 23. The acceptance rules (N simple roots; a rational of bounded size; values that are square roots of −1) do not see δ. The seeded search for a coordinate at q = 29 is accepted on the number of simple roots. The two sides meet in the configurations (b, c, e, f) and nowhere else. A surface is a configuration and its associate on both sides; `stored` asserts it.

### 4.2 Simpler candidates [E]

`stored`, over the 142 records with status "ok" (106 curves and 36 controls). For each candidate, the number of records on orbit A and on orbit B for each of its two values; a candidate separates the orbits when each value occurs on one orbit only.

| candidate | q = 29, e₁ = 2 | q = 37, e₁ = 2 | q = 37, e₁ = 18 | separates |
|---|---|---|---|---|
| δ | −: 31, 0; +: 0, 18 | −: 0, 24; +: 27, 0 | −: 23, 0; +: 0, 19 | yes |
| u alone | −: 21, 10; +: 10, 8 | −: 15, 11; +: 12, 13 | −: 9, 8; +: 14, 11 | no |
| sign(Re c₁·Im c₁) | −: 19, 10; +: 12, 8 | −: 11, 7; +: 16, 17 | −: 16, 15; +: 7, 4 | no |
| sign(c₂) | −: 25, 14; +: 6, 4 | −: 24, 14; +: 3, 10 | −: 16, 11; +: 7, 8 | no |
| the two signs together, without u | −: 10, 10; +: 21, 8 | −: 15, 13; +: 12, 11 | −: 14, 8; +: 9, 11 | no |
| u·sign(Re c₁·Im c₁) | −: 6, 14; +: 25, 4 | −: 24, 10; +: 3, 14 | −: 7, 11; +: 16, 8 | no |
| u·sign(c₂) | −: 12, 10; +: 19, 8 | −: 11, 17; +: 16, 7 | −: 7, 15; +: 16, 4 | no |
| the configuration, not its twist, has P (the sign of S₁) | −: 13, 10; +: 18, 8 | −: 8, 11; +: 19, 13 | −: 14, 5; +: 9, 14 | no |
| u times that flag | −: 12, 12; +: 19, 6 | −: 13, 16; +: 14, 8 | −: 13, 11; +: 10, 8 | no |
| comparison flag of note 14 (equal, twist) | 13, 10; 18, 8 | 11, 15; 16, 9 | 17, 11; 6, 8 | no |

Also not: the kind of branch points (four values), M, the order of X. No configuration of the three torsors is self-associated. By surface, the pair of flags "has P" of its two configurations takes each of its three values on both orbits in the first torsor and is no function of the orbit in the others. The factors of δ do their work only together: the signs alone are a function of the quartic, which is the same up to units and conjugation for all surfaces of a torsor, and u alone changes with the pair.

### 4.3 The orbits from the class polynomial alone [E]

`galois` computes H as note 20 does (lattices, period points, recognition over E) and then uses no ideal. PARI factors H over E and over E(i); one factor over E(i) is reduced at each of the two primes above 𝔭 = (q, F − e₁).

| q | e₁ | digits | H over E | H over E(i) | roots modulo q of one factor, at the two primes | equals the orbits of note 24 | equals the split by δ |
|---|---|---|---|---|---|---|---|
| 37 | 18 | 1,200 | 8 | 4, 4 | 11, 18, 22, 30 and 16, 17, 23, 31 | yes | yes |
| 29 | 2 | 1,200 | 12 | 6, 6 | 7, 10, 20, 21, 24, 26 and 0, 3, 6, 14, 18, 23 | yes | yes |
| 37 | 2 | 2,000 | 12 | 6, 6 | 0, 8, 11, 16, 27, 29 and 7, 13, 15, 22, 26, 33 | yes | yes |

This shares with note 24 the lattices, the period points and H, and replaces the rest of its computation of the orbits (permutations of ideals located numerically, orbits, and the recognition of Γ of note 24 §3) by an exact factorisation. It also shows directly that the quadratic extension is E(i), for Z/4 × Z/2 as well, where there are three candidates.

The classes, in the same command:

| q | e₁ | G | primes of degree 1 with ℓ ≡ 1 mod 4, classes | with ℓ ≡ 3 mod 4, classes | the two sets of classes | degrees of the factors of H mod 𝔩 equal the order of the class |
|---|---|---|---|---|---|---|
| 37 | 18 | Z/4 × Z/2 | 40, 4 | 44, 4 | disjoint; the first a subgroup | 84 of 84 |
| 29 | 2 | Z/12 | 40, 6 | 70, 6 | disjoint; the first a subgroup | 110 of 110 |
| 37 | 2 | Z/12 | 40, 6 | 46, 6 | disjoint; the first a subgroup | 86 of 86 |

### 4.4 The form modulo 4 on the lattices [D, E]

[D] Let (L, ξ) be a polarised lattice of note 20 §3, with E(x, y) = Tr_{K/Q}(ξxȳ). K ⊗ Q₂ is the product of K_{w₁} = Q₂(i) and a field of degree 6, the two factors are orthogonal for E because w₁ is its own conjugate, and Frobenius acts on K_{w₁} through the root π₁. So K_{w₁} is the line V₁ of note 23 §5.2. For m in it, h(m) = E(m, i·m) = Tr(−iξ·N(m)) = 2·N(m)·t, where t is the image in Q₂ of −iξ, an element of K⁺. The class of h modulo norms and powers of 2 is the unit part of t modulo 4. It does not change when (L, ξ) is replaced by (vL, ξ/(v v̄)).

[E] `galois` finds π₁ in Z₂[i] to 188 binary digits, evaluates −iξ there for each even lattice, checks that the result lies in Q₂, and takes the unit part modulo 4 as a sign s. Then, as in note 24 §3 but with s in place of the orbit: the indicator Σ_j s(j) ∏_{l≠j}(X − x_l) of each type, the eight recognitions, the reduction at 𝔭.

| q | e₁ | lattices per type with s = +1, −1 | recognised with | values of Γ/H′ at the roots | split of the surfaces equals the orbits of note 24 | equals the split by δ |
|---|---|---|---|---|---|---|
| 37 | 18 | 4, 4 in each type | divisor i, signs (+, +), and no other | 6, 31 | yes | yes |
| 29 | 2 | 6, 6 in each type | divisor i, signs (+, +), and no other | 12, 17 | yes | yes |
| 37 | 2 | 6, 6 in each type | divisor i, signs (+, +), and no other | 6, 31 | yes | yes |

Unlike the name of an orbit, s has a meaning in each type, and the three types fit one polynomial over E with no sign between them.

### 4.5 Reading

Four descriptions of one split of 32 surfaces agree: orbits of ideals above primes ℓ ≡ 1 mod 4 (note 24), factors of H over E(i), the form modulo 4 of the lattice, and δ from the Tate pairing in characteristic p. The first three are computed over C from the same lattices and are not independent of each other as computations; they are three different functions of those lattices. The fourth shares nothing with them but the configurations. I see no way for the fourth to match by an accident of code, and the chance figure of note 24 (1 in 7,470,540 for a balanced function) is computed correctly.

What is established is narrower than "δ is the character" may suggest and is worth stating once: on these three torsors, δ of a surface over F_q equals a number read off the lattice of its CM lift, for which the character property is a line of algebra. What is not established is in §7.

---

## 5. Sources read for this check

| Statement in the checked notes | Status |
|---|---|
| Note 20 §2, the formulas of the two papers | read [P], finding 11 |
| Note 23 §5.1 and note 24 §5, the method and scope of Castryck, Sotáková and Vercauteren | abstract, §5.3 and conclusion read [P]; finding 5. §3.2 and Theorem 10, quoted in note 23, not re-read |
| Note 20 §7.1, the theorem of Shimura and Taniyama | not read; still recalled |
| Note 20, Matsumoto (1993) | not read |
| Note 23 §2.2, §2.3, the methods of Khuri-Makdisi and of Bruin | not read; still recalled |

---

## 6. Statements affected

The findings bear on these statements, which the notes concerned state in the form the findings give:

- **Note 24**: §1 item 1 (checked, and the two further routes); §1 item 6 and the paragraph after item 7 (CSIDH, finding 5); §1 item 7 (the step outside the kernel, finding 2); §2 (finding 2); §3 (finding 3); §6, the bullet on two primes (53 of the 80 polynomials have q ≡ 1 mod 8).
- **Note 23**: §1 item 8 (checked; the form on the lattices); §3, the pairing (finding 9); §4.3 (finding 10); §5.3 (finding 7); §6, the bullet on the lattice statements (finding 4).
- **Note 20**: §1 (checked); §3 (finding 12); §7.2, the field with half of G (finding 13).
- **Note 13** §9 and **note 19** §1: the two places that say "as for CSIDH" (finding 5).

---

## 7. Still unverified

- **Any torsor beyond the three.** Two primes, both 5 mod 8; groups of order 8 and 12. The nearest test that is out of sample is q = 41, where note 19 lists two polynomials with the character (e₁ = 6 with G = Z/16, e₁ = 22 with G = Z/12) and 41 ≡ 1 mod 8, so the pairing has order 32. It needs the enumeration at p = 41, the curves of note 14 there, and the CM method at 16 lattices per type. With §4.4 the lattice side of such a test is cheap: s is read off ξ. This test is carried out in [note 26](26-delta-is-the-character-at-q-41.md), and the claim holds there: 16 surfaces of 16 for e₁ = 6 and 12 of 12 for e₁ = 22, with the cosets from the orbits of ideals, from the factors over E(i) and from ξ modulo 4, all three the same. The pairing over F_{41⁴} has order 32 as said here; the eigenvalue has π₁⁴ − 1 of valuation 8, so u is still read modulo 16. Half of the curves have exponent 64 and are not computed. The sign between δ and the form of the lattice, which finding 4 leaves free for each polynomial, follows one rule on all five torsors (note 26 §7.1; found after the fact). Note 26 is not re-checked.
- **That δ equals the lattice form in general**, and why. §4.4 is a measurement on 32 surfaces. The Prym variety of note 14 has a polarisation of type (1, 1, 1, 2) and the lattice a principal one; that the two give the same class at w₁ for every pair of lines is observed, not derived. My derivation of the lattice form is [D] and unchecked.
- **The 126 pairs with exponent above 16**, as note 23 says.
- **The group law and Miller's algorithm** of `sixlines_jacobian.py`: not read. Their support is note 23 §3, the seven identities of finding 9, and that everything built on them comes out.
- **The enumeration of the lattices of R** (note 20 §5.2) and the count of 234; the long sweeps of note 20, of which I reran two fields of fifteen; the runs of note 23, of which I reran ten curves of 242 and recounted the rest from the stored files.
- **The comparison of permutations (note 20 §7.3) on the two torsors of 12 surfaces**; for q = 37, e₁ = 18 it stands at 10 cases.
- **That Frobenius at 𝔩 is the step of the class of 𝔩 for ℓ ≡ 3 mod 4**: the theorem recalled in note 20; the cycle lengths agree, the permutation was not compared with anything, because no ideal above such a prime is applied to a lattice.
- **That the CM points are the canonical lifts** (note 20 §7.7), on which the identification of the action of note 20 with the action on the K3 surfaces rests.
- **The level test** for the three polynomials note 23 could not reach, and which surface is on which level on the lattice side.
- **Class groups** are PARI's without certification, here as there; the factorisations over E and E(i) are PARI's `nffactor`.
- Matsumoto (1993), Shimura–Taniyama, Khuri-Makdisi, Bruin: not read.

---

## 8. Reproduce

From `experiments/`:

```bash
python -I sixlines_review2.py stored
python -I sixlines_review2.py levels
python -I sixlines_review2.py theta
python -I sixlines_review2.py pairing 29 2 3 12 27 1 3 --seed 3
python -I sixlines_review2.py delta --tw 3 --seed 2026 --workers 8 --out sixlines_review2.jsonl
python -I sixlines_review2.py galois 37 18 1 1 --digits 1200 --workers 8 --out sixlines_review2.jsonl
python -I sixlines_review2.py galois 29 2 1 1 --digits 1200 --workers 8 --out sixlines_review2.jsonl
python -I sixlines_review2.py galois 37 2 1 1 --digits 2000 --workers 8 --out sixlines_review2.jsonl
```

A second each for the first two, 15 seconds, 10 seconds, one minute, one minute, a minute and a half, and six minutes, with 8 processes. The output of the eight commands, in this order, is `sixlines_review2.log`; the last four also append records to `sixlines_review2.jsonl`. They need the stored files of notes 2, 14, 20, 23 and 24 and import `sixlines_cmmethod.py`, `sixlines_delta_labels.py` and `sixlines_jacobian.py` without changing them. The random classes of `pairing` and `delta` are seeded.

The reruns of the scripts of the checked notes quoted in §2 (output not kept in this repository):

```bash
python -I sixlines_delta_labels.py run 37 18 1 1 --digits 1200 --primes 8 --workers 8
python -I sixlines_cmmethod.py selftest --digits 60
python -I sixlines_cmmethod.py sweep 17 23 --select 17,2 23,10 --digits-list 300 700 --primes 3 --all-subsets --maxn 8 --workers 2
python -I sixlines_cmmethod.py points 17 2 1 1 --digits 30 --maximal
python -I sixlines_cmmethod.py points 17 2 1 1 --digits 30 --nolevel
python -I sixlines_jacobian.py --lattice
python -I sixlines_jacobian.py --one 29 2 3 5 20 1 2
```

The `sweep` and the `run` take an `--out` file; I gave them one outside `experiments/`, so nothing was appended to the stored files.

## Sources

- Matsumoto, Terasoma, *Thomae type formula for K3 surfaces given by double covers of the projective plane branching along six lines*, https://arxiv.org/abs/1001.4865 : introduction with formula (4), §3.3 Proposition 6, §3.4 with (19), Remarks 1 and 2, Table 1, Proposition 7 [P].
- Clingher, Malmendier, Shaska, *Six line configurations and string dualities*, https://arxiv.org/abs/1806.07460 : (2.1), (2.2), (2.12) to (2.16), Theorem 2.8, (2.30) to (2.34), Lemma 2.17, Corollary 2.25 [P].
- Castryck, Sotáková, Vercauteren, *Breaking the decisional Diffie–Hellman problem for class group actions using genus theory*, CRYPTO 2020, https://eprint.iacr.org/2020/151 : abstract, §5.3, conclusion [P].
- In this repository: notes [2](02-six-line-count.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [18](18-checks-of-notes-13-to-16.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [24](24-the-form-mod-4-is-the-character.md).

