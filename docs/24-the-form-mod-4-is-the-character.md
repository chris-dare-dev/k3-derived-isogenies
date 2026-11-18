# 24. The form mod 4 is the character χ₋₁

*[Note 23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md) §5 computes a sign δ for each surface of three torsors, from a Tate pairing on the Prym variety of the curve of [note 14](14-fourfold-of-a-six-line-surface.md), and leaves open whether δ is the character χ₋₁ of [note 19](19-quadratic-characters-from-a-surface.md) §1 item 7: whether δ(g·X) = χ₋₁(g)·δ(X). Its §7 names the computation that would decide it: labels for those torsors by the CM method of [note 20](20-cm-method-and-the-action.md). This note does that computation. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_delta_labels.py`; **[H]** heuristic. Statements of other notes keep the label they have there. This note and notes 20 and 23, on which it rests, are re-checked in [note 25](25-checks-of-notes-20-23-24.md).*

Notation as in notes 13 and 19: P the Frobenius polynomial of the surface, E = Q[F] of degree 6, G the group that acts simply transitively on the set S_P of surfaces with polynomial P, 𝔭 = (q, F − e₁) the prime of E used by note 20. χ₋₁ is the character of G of the unramified extension E(i)/E, where it is one; its value on the step at a prime of degree 1 above ℓ is (−1/ℓ).

---

## 1. Result

1. **δ is the character on all three torsors [E].** The steps at primes ℓ ≡ 1 mod 4 lie in the kernel of χ₋₁. On the lattices of note 20 they generate a subgroup of index 2 of G in each of the three fields, with two orbits. Carried to the surfaces over F_q, the two orbits are exactly the set where δ = +1 and the set where δ = −1:

   | q | e₁ | G | surfaces | primes ℓ ≡ 1 mod 4 used | ideals | group they generate | orbits | δ on orbit A, on orbit B | surfaces that agree |
   |---|---|---|---|---|---|---|---|---|---|
   | 29 | 2 | Z/12 | 12 | 8 | 48 | order 6 | 6 and 6 | −1, +1 | 12 of 12 |
   | 37 | 2 | Z/12 | 12 | 8 | 48 | order 6 | 6 and 6 | +1, −1 | 12 of 12 |
   | 37 | 18 | Z/4 × Z/2 | 8 | 8 | 48 | order 4 | 4 and 4 | −1, +1 | 8 of 8 |
   | all | | | 32 | 24 | 144 | | | | 32 of 32 |

   Since δ is constant on each coset of the kernel and takes opposite values on the two, δ(g·X) = χ(g)·δ(X) for every g in G, with χ the character whose kernel is that subgroup. These are all the torsors on which note 23 computes δ; [note 26](26-delta-is-the-character-at-q-41.md) adds two at q = 41 (§6).

   [Note 25](25-checks-of-notes-20-23-24.md) (finding 1, §4) found no error here: the rerun of one field gives the stored record. The same split of the 32 surfaces comes out by two routes that apply no ideal to a lattice: the two factors of the class polynomial over E(i), and the class modulo 4 of the polarisation ξ of each lattice at the prime of degree 1 above 2, which is note 23's u read on the lattice. None of 13 simpler functions of the stored data separates the orbits.
2. **The character is χ₋₁ [E].** The function that is +1 on one orbit and −1 on the other is i times a function defined over E: of the eight ways to recognise it (divided by 1 or by i, with a sign for each type), exactly one succeeds in each field, and it is a division by i. Its values at the surfaces are the two square roots of −1 in F_q (6 and 31 modulo 37, 12 and 17 modulo 29). So the two orbits are separated by E(i)/E, as note 19 said they would be. For G = Z/12 this is the only character of order 2; for Z/4 × Z/2 it is one of three.
3. **By chance? [D]** Note 23 already knew that δ takes each value on half of each torsor. A function with that property agrees with a given split into two halves, up to sign, with probability 2/C(8, 4) = 1/35 on 8 surfaces and 2/C(12, 6) = 1/462 on 12. For the three together: 1 in 7,470,540.
4. **Two independent computations agree [E].** The orbits come from lattices over C: period points, theta constants at 1,200 or 2,000 digits, a class polynomial over E, reduction modulo 𝔭. δ comes from arithmetic in the Jacobian of a curve of genus 5 over F_{p⁴}. They share the list of surfaces and nothing else. Note 20 says of itself that it is not a test of the Prym variety of note 14. This is one: an invariant of that Prym variety which is not an isogeny invariant follows the action computed on the lattices.
5. **The CM method on two more fields [E].** Note 20 skipped the two torsors of 12 surfaces. With the period points computed in parallel they take 75 and 360 seconds, and the class polynomial modulo 𝔭 has exactly the 12 enumerated surfaces as roots in both. That brings note 20's comparison to 98 surfaces of 98 in 17 fields.
6. **What follows for the decision problem [D, from items 1 and 2; H for the cost].** For a polynomial P where χ₋₁ is a character of G, and two surfaces X and Y of S_P, the bit χ₋₁(g) of the class g with Y = g·X is δ(X)·δ(Y), and it is computed from the two surfaces without knowing g. For four surfaces X, aX, bX, cX, the relation c = ab forces δ(X)·δ(cX) = δ(aX)·δ(bX), which a random c satisfies half of the time. That is the attack of Castryck, Sotáková and Vercauteren on the decisional Diffie–Hellman problem for class-group actions, here through the fourfold. Their attack applies to class-group actions with even class number; it does not affect CSIDH, whose class group has no 2-torsion (their §5.3 and conclusion; [note 25](25-checks-of-notes-20-23-24.md), finding 5). It applies to the polynomials of note 19 item 7: 80 of 477 in its range, none with q ≡ 3 mod 4. The remedies are the countermeasures of their §5.3: use only steps in the kernel, or a polynomial without the character.
7. **What this is not.** Not a proof: three torsors, two primes, 32 surfaces (five torsors at three primes, 60 surfaces, with [note 26](26-delta-is-the-character-at-q-41.md)). Not an attack on the search problem: one bit of g, for a minority of polynomials. And no ideal above a prime with χ₋₁ = −1 was applied to a lattice; that δ changes sign under such a step is inferred from there being two orbits with opposite values (§2). The step outside the kernel is covered on the Galois side in [note 25](25-checks-of-notes-20-23-24.md) (finding 2): the class polynomial is irreducible over E and splits into two factors over E(i) whose roots are the two orbits, so Frobenius at a prime of degree 1 with ℓ ≡ 3 mod 4 exchanges the orbits; at 160 such primes its cycles have the length of the order of the class in G.

So the conclusion of note 13 stands and is sharpened in one place. The problem is of the type of CSIDH, and the attack of Castryck, Sotáková and Vercauteren on the decision problem for class-group actions with even class number carries over to it, for the polynomials where the character exists. That nothing computed in note 19 distinguishes g·X from a random surface of the torsor does not hold for those polynomials, as note 19 §1 records.

---

## 2. The prediction [D]

Note 20 §7.2 computes, for a prime ℓ ≡ 1 mod 4 that splits completely in K and each of the six pairs {j, k}, the permutation of the lattices of the surfaces induced by the ideal 𝔓_j𝔓_k·conj(𝔓_l𝔓_m), and finds that it depends only on the class in G of the prime of E at a_j a_k. That prime has degree 1 above ℓ, so χ₋₁ of its class is (−1/ℓ) = +1 (note 19, Lemma 1). Hence:

- every such ideal acts through an element of ker χ₋₁;
- if these elements generate ker χ₋₁, their orbits on S_P are the two cosets of ker χ₋₁;
- δ(g·X) = χ₋₁(g)·δ(X) for all g then says exactly: δ is constant on each orbit, and the two constants are opposite.

The converse is what makes the test complete. If δ is constant on each of two orbits of a subgroup H of index 2, with opposite values, then for g in H the surfaces X and g·X are in one orbit and δ agrees, and for g outside H they are in different orbits and δ differs. So δ(g·X) = χ_H(g)·δ(X) for every g in G, although no element outside H was ever applied. What a direct step at a prime ℓ ≡ 3 mod 4 would add is a check that the group acting is G and not a smaller group with two orbits. For G = Z/12 that is covered by note 20's count (the lattices of each type number |G|) and by the order of the subgroup found, |G|/2. For G = Z/4 × Z/2 the count does not cover it: there an element outside the subgroup could act trivially and the picture would be the same. It is covered, in all three fields, by the irreducibility of the class polynomial over E and by the cycle lengths of Frobenius at primes of both classes modulo 4 ([note 25](25-checks-of-notes-20-23-24.md), finding 2).

The comparison is between two partitions of a set of surfaces. No sign enters: the overall sign of δ for a polynomial is a convention (note 23 §5.2), and so is the name of an orbit.

---

## 3. From orbits of lattices to orbits of surfaces [D]

Note 20 carries the action to F_q through polynomials A_s over E that express the permutation of a class on the roots of the class polynomial H, and for q = 37, e₁ = 18 these were not recognised at the precision reached. Only the split into two orbits is needed here, and it has a smaller description.

For each of the three types T let x_{T,1}, …, x_{T,N} be the values of the coordinate at the lattices, so that ∏(X − x_{T,j}) is the image of H under the embedding of E belonging to T. Let s_T(j) = +1 on one orbit of the subgroup and −1 on the other, and

> I_T(X) = Σ_j s_T(j) ∏_{l ≠ j} (X − x_{T,l}).

The function s on the roots of H is fixed by the subgroup and negated by the rest of G. Under the Galois action on the roots, which note 20 §7.3 identifies with the action of G, s is therefore √d times an element of E at each root, with E(√d) the quadratic extension belonging to the character; for χ₋₁, d = −1. So I_T/i, with a sign ε_T for each type (the name of the orbit, and the value of i under the embedding), is the image of one polynomial Γ in E[X]. The script tries the divisors 1 and i and the four choices of (ε₂, ε₃), with the acceptance rule of note 20 (a rational of total size below 0.6 times the working digits within 10^−(digits − 15)).

On these three torsors the identification of note 20 §7.3 rests on 10 tested cases, all for q = 37, e₁ = 18; the two torsors of 12 surfaces were not in note 20, and this script has no Frobenius test. [Note 25](25-checks-of-notes-20-23-24.md) (finding 3 and §4.3) adds the cycle lengths at 280 primes, which agree with the orders of the classes, and a factorisation of H over E(i) by PARI, which gives the same two sets of roots without the recognition of Γ.

Modulo 𝔭, at a root r of H, the value Γ(r)/H′(r) is s(r)/i, one of the two square roots of −1 in F_q. Surfaces with the same value are in the same orbit. That the values are square roots of −1 is asserted by the script; it would fail for a wrong d or a wrongly recognised Γ.

The surface of a root is identified as in note 20 §6: the five invariants (seven in the chart with J₃ in the denominator) at the root are those of an enumerated configuration.

Two changes to the procedure of note 20, both in the new script and neither in `sixlines_cmmethod.py`:

- The period points are computed in forked worker processes. The functions called are the same.
- For q = 29, e₁ = 2 one surface has J₂ ≡ 0 modulo 29, and in the other chart each of the eight coordinates that note 20 tries gives two of the twelve surfaces the same value modulo 29 (8 to 11 distinct roots). The script then tries pseudo-random integer combinations of the seven invariants from a fixed seed, with the same acceptance rule: N simple roots modulo 𝔭. The first that passes is (3, −5, −7, −2, −7, 1, 6).

---

## 4. Results [E]

`sixlines_delta_labels.py run q e1 1 1 --primes 8`. δ is read from `sixlines_jacobian_character.jsonl` and `sixlines_jacobian_character_extra.jsonl` of note 23 (records with status "ok"; the count includes its controls with twists and other random classes).

| q | e₁ | digits | chart, coordinate | roots of H mod 𝔭 that are enumerated surfaces | primes ℓ | classes of G seen | order of the group generated, by type | Γ recognised with | values of Γ/H′ |
|---|---|---|---|---|---|---|---|---|---|
| 29 | 2 | 1,200 | J₃, (3, −5, −7, −2, −7, 1, 6) | 12 of 12 | 89, 733, 1,453, 1,493, 1,933, 2,153, 3,061, 3,581 | 6 of 12 | 6, 6, 6 | divisor i, signs (−, −) | 12, 17 |
| 37 | 2 | 2,000 | J₂, (0, 1, 0, 0, 0) | 12 of 12 | 241, 421, 661, 1,181, 1,381, 1,801, 2,521, 3,413 | 6 of 12 | 6, 6, 6 | divisor i, signs (+, +) | 6, 31 |
| 37 | 18 | 1,200 | J₂, (1, 2, 3, 0, 0) | 8 of 8 | 941, 2,713, 4,297, 4,357, 4,457, 4,561, 4,673, 4,817 | 4 of 8 | 4, 4, 4 | divisor i, signs (+, +) | 6, 31 |

- In each field exactly one of the eight recognitions succeeds.
- For q = 37, e₁ = 2 the polynomials of the invariants are not recognised at 1,200 digits; at 2,000 they are.
- With three primes per field in place of eight the groups and orbits are the same (the first runs, not kept). No prime ℓ ≡ 1 mod 4 joins the two orbits: 144 ideals.

The surfaces, by orbit. A configuration is (b, c, e, f) of note 2; a surface is a configuration and its associate.

**q = 29, e₁ = 2**

| orbit | root of the class polynomial mod 𝔭 | configurations (b, c, e, f) | δ | records of note 23 |
|---|---|---|---|---|
| A | 0 | (2, 3, 12, 27), (2, 3, 6, 17) | −1 | 7 |
| A | 3 | (2, 3, 7, 12), (2, 5, 15, 27) | −1 | 7 |
| A | 6 | (2, 3, 8, 7), (2, 3, 5, 20) | −1 | 6 |
| A | 14 | (2, 5, 19, 14), (2, 5, 19, 20) | −1 | 2 |
| A | 18 | (2, 5, 4, 28), (2, 3, 22, 24) | −1 | 4 |
| A | 23 | (2, 5, 6, 16), (2, 5, 24, 22) | −1 | 5 |
| B | 7 | (2, 5, 20, 4), (2, 3, 12, 14) | +1 | 3 |
| B | 10 | (2, 5, 14, 20), (2, 5, 13, 22) | +1 | 2 |
| B | 20 | (2, 3, 7, 17), (2, 5, 11, 7) | +1 | 6 |
| B | 21 | (2, 3, 9, 23), (2, 5, 9, 14) | +1 | 3 |
| B | 24 | (2, 3, 10, 25), (2, 5, 11, 28) | +1 | 2 |
| B | 26 | (2, 3, 25, 19), (2, 5, 6, 20) | +1 | 2 |

**q = 37, e₁ = 2**

| orbit | root of the class polynomial mod 𝔭 | configurations (b, c, e, f) | δ | records of note 23 |
|---|---|---|---|---|
| A | 0 | (2, 5, 23, 20), (2, 6, 10, 34) | +1 | 3 |
| A | 8 | (2, 3, 20, 6), (2, 3, 29, 15) | +1 | 4 |
| A | 11 | (2, 5, 18, 31), (2, 5, 31, 30) | +1 | 3 |
| A | 16 | (2, 3, 9, 7), (2, 3, 10, 18) | +1 | 5 |
| A | 27 | (2, 3, 6, 30), (2, 5, 7, 26) | +1 | 6 |
| A | 29 | (2, 3, 4, 18), (2, 3, 4, 26) | +1 | 6 |
| B | 7 | (2, 5, 32, 20), (2, 5, 21, 27) | −1 | 2 |
| B | 13 | (2, 3, 14, 33), (2, 5, 29, 20) | −1 | 5 |
| B | 15 | (2, 5, 24, 33), (2, 5, 6, 8) | −1 | 3 |
| B | 22 | (2, 3, 24, 19), (2, 3, 26, 29) | −1 | 4 |
| B | 26 | (2, 3, 8, 26), (2, 5, 9, 22) | −1 | 6 |
| B | 33 | (2, 3, 17, 13), (2, 5, 20, 17) | −1 | 4 |

**q = 37, e₁ = 18**

| orbit | root of the class polynomial mod 𝔭 | configurations (b, c, e, f) | δ | records of note 23 |
|---|---|---|---|---|
| A | 11 | (2, 5, 25, 20), (2, 6, 30, 24) | −1 | 2 |
| A | 18 | (2, 3, 10, 24), (2, 3, 20, 35) | −1 | 9 |
| A | 22 | (2, 3, 8, 29), (2, 3, 26, 21) | −1 | 7 |
| A | 30 | (2, 5, 35, 19), (2, 3, 11, 5) | −1 | 5 |
| B | 16 | (3, 4, 23, 32), (3, 7, 9, 8) | +1 | 5 |
| B | 17 | (2, 3, 20, 7), (2, 3, 16, 27) | +1 | 2 |
| B | 23 | (2, 5, 19, 33), (2, 3, 10, 17) | +1 | 7 |
| B | 31 | (2, 5, 15, 27), (2, 5, 15, 30) | +1 | 5 |

In every row δ is the same on all records of the surface (note 23's result), and in every table δ is constant on A, constant on B, and the two differ.

---

## 5. What it means for the decision problem

**The statement [D].** Let P be a polynomial for which χ₋₁ is a character of G: by note 19 item 7, q ≡ 1 mod 4 and the real quartic field K⁺ of the fourfold has a prime of degree 1 above 2. If δ(g·X) = χ₋₁(g)·δ(X) holds on S_P, then for any two surfaces of S_P

> χ₋₁([Y : X]) = δ(X)·δ(Y),

where [Y : X] is the class that carries X to Y. This is established on three torsors (§1) and is a conjecture beyond them.

**What it breaks.** The decisional problem: given X, aX, bX and Z, decide whether Z = abX. If so, δ(X)·δ(Z) = δ(aX)·δ(bX). A surface Z drawn at random from S_P satisfies this with probability 1/2. One bit, as in the paper of Castryck, Sotáková and Vercauteren for class-group actions with even class number, such as those of supersingular curves over F_p with p ≡ 1 mod 4 [P, in note 23 §5.1] (their attack does not affect CSIDH, whose class group has no 2-torsion; [note 25](25-checks-of-notes-20-23-24.md), finding 5). It does not help to find a, and it says nothing when a and b are restricted to ker χ₋₁, for instance when only steps at primes ℓ ≡ 1 mod 4 are used.

**How many polynomials.** 80 of the 477 polynomials of note 19 that have a character of order 2 (q ≤ 10,007), none with q ≡ 3 mod 4. The other quadratic characters of G (note 19 §2.5: up to 4 per group) have no description of this kind, and nothing here bears on them.

**Cost [H].** Note 23 §5.6: given the six lines and one usable pair, δ needs the quartic χ of the curve and arithmetic in Jac(C)(F_{p⁴}), about 75 seconds at p = 29 and 37. χ is known from P up to the unit of c₁ and the twist (note 13, Proposition Y; note 23 §5.2), a bounded number of candidates, and the right one is the one whose value at Frobenius kills a random class of the Prym part (note 23 §3 does this as a validation). So the computation should not need a point count of the curve of genus 5, and every step is polynomial in log p in principle. Nobody has run it beyond p = 37, about half of the pairs of lines are unusable (2-primary torsion of exponent above 16, note 23 item 9), and the linear algebra of the group law has not been costed.

**Why this does not contradict note 19.** Note 19 §2.2 shows that the K3 lattices of a torsor form one genus, so no local invariant of the lattice of the surface is a character, and concludes that the genus-theory method of Castryck, Sotáková and Vercauteren has nothing to read on the K3 side. That stands. δ is read on the fourfold, where the hermitian form over Z[i] has a class modulo 4 at the prime above 2 that the second exterior power forgets (note 23 §5.6). The surface determines the fourfold, by Paranjape's construction, so the bit is computable from the surface all the same.

---

## 6. Check and limits

- **Re-derivation check.** [Note 25](25-checks-of-notes-20-23-24.md) checks this note and notes 20 and 23. The agreement of §1 is evidence for both of them as well as for the statement tested: a wrong period map, a wrong class polynomial or a wrong pairing would have to produce the same split of 8, 12 and 12 surfaces.
- **Three torsors, two primes.** q = 29 and 37, both 5 modulo 8; groups of order 8 and 12. Note 19 lists 80 polynomials; 3 are tested. The five torsors over F₄₉, where q ≡ 1 mod 8 and Frobenius conjugation is a known element with χ₋₁ = −1, are not handled by the script of note 23. Of the 80 polynomials, 53 have q ≡ 1 mod 8 and 27 have q ≡ 5 mod 8 (note 19 §2.5), so the kind not tested here is the larger ([note 25](25-checks-of-notes-20-23-24.md), finding 6). [Note 26](26-delta-is-the-character-at-q-41.md) tests q = 41 ≡ 1 mod 8, on the two polynomials with the character there (e₁ = 6 with G = Z/16, e₁ = 22 with G = Z/12), and the statement holds: δ is constant on each coset of ker χ₋₁ and opposite between them on 16 surfaces of 16 and 12 of 12. That makes five torsors at three primes, 60 surfaces of 60. δ is given by the same formula; u is read from a pairing of order 32 over F_{41⁴}, by a rule written down before the curves were computed. Note 26 is not re-checked. The torsors over F₄₉ are still not done.
- **No step outside the kernel.** Only primes ℓ ≡ 1 mod 4 that split completely in K are used, as in note 20. The sign change of δ under the other coset is the inference of §2.
- **"The action" is that of note 20.** It is the action of ideals on the lattices of the fourfolds that Lombardo–Peters–Schütt attach to the surfaces, carried to F_q by reduction of CM points. That it is the class-group action on the K3 surfaces of note 13 is the sketch of note 20 §7.7.
- **Recognition.** Γ is accepted by the rule of note 20. A wrong acceptance would have to give values that are square roots of −1 at all N roots, and the right split.
- **Numerical location** of the permuted lattices is at 50 digits, with note 20's assertion (nearest point within 10⁻¹², next beyond 10⁻⁶).
- **The count of §1 item 3** assumes only that δ is balanced. It does not assume that δ is otherwise random; a δ that is a different character would also be balanced. For Z/12 there is no other character of order 2. For Z/4 × Z/2 there are two others, so there the chance is better put as 1 in 3 against that alternative.
- **Class groups** are PARI's without certification (GRH).
- **Coordinates for q = 29, e₁ = 2** come from a seeded pseudo-random search. The acceptance rule does not look at δ.

---

## 7. What this does not give

- A proof that δ(g·X) = χ₋₁(g)·δ(X), on these torsors or in general. Note 23 §5.2 gives the reason to expect it; its steps are sketches.
- The statement for the other 77 polynomials of note 19, for any q ≡ 1 mod 8, or for any p above 37. Two of the 77 are tested in [note 26](26-delta-is-the-character-at-q-41.md), both at q = 41 ≡ 1 mod 8, and the statement holds on them (28 surfaces of 28). Still open after that: the other 75, any q ≡ 1 mod 16, any p above 41.
- A step at a prime ℓ ≡ 3 mod 4, on lattices or on surfaces.
- Full labels in G for the surfaces of these three torsors. Only the coset of ker χ₋₁ is computed. The polynomials A_s of note 20 §7.5 were not attempted at the higher precision.
- Anything about the other characters of order 2 of G, the polynomials without χ₋₁, or non-maximal orders.
- An attack on the search problem, or any change to the sizes of note 13 §3.3.
- A costed algorithm for δ at cryptographic size.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_delta_labels.py run 37 18 1 1 --digits 1200 --primes 8 --workers 14 --out sixlines_delta_labels.jsonl
python -I sixlines_delta_labels.py run 29 2 1 1 --digits 1200 --primes 8 --workers 14 --out sixlines_delta_labels.jsonl
python -I sixlines_delta_labels.py run 37 2 1 1 --digits 2000 --primes 8 --workers 14 --out sixlines_delta_labels.jsonl
```

The three were run in this order with their output collected in `sixlines_delta_labels.log`: 83, 109 and 395 seconds on 16 cores. Each appends one record to `sixlines_delta_labels.jsonl`. The script imports `sixlines_cmmethod.py` and needs its inputs (`sixlines_cm.jsonl`, `sixlines_p29.jsonl`, `sixlines_p37.jsonl`) and the two files `sixlines_jacobian_character.jsonl` and `sixlines_jacobian_character_extra.jsonl` of note 23; it does not recompute δ. No package was installed.

## Sources

No new source was read. In this repository: notes [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [19](19-quadratic-characters-from-a-surface.md), [20](20-cm-method-and-the-action.md), [23](23-jacobian-arithmetic-levels-and-the-form-mod-4.md), [25](25-checks-of-notes-20-23-24.md), [26](26-delta-is-the-character-at-q-41.md). The paper of Castryck, Sotáková and Vercauteren (CRYPTO 2020, https://eprint.iacr.org/2020/151) is cited as read in note 23 §5.1.
