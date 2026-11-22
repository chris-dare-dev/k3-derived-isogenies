# 21. Which fourfold over F_p belongs to a surface of sign − or with reducible P

*Takes up the item that [note 14](14-fourfold-of-a-six-line-surface.md) §5.2 and [note 18](18-checks-of-notes-13-to-16.md) §4.4 leave open: for surfaces of sign −, and for some pairs of lines on surfaces with reducible P, the curve of note 14 is of the "wrong" type, and the form of the fourfold that belongs to the surface over F_p is not identified there. This note identifies it. Labels: **[P]** read in the cited source; **[M]** recalled, not checked; **[D]** derived here; **[E]** measured by `experiments/sixlines_forms_fp.py`; **[H]** heuristic. Not re-checked: the derivations are mine, and where they enter a test the only support is that the test comes out. It does not move the hardness conclusion of [note 13](13-hardness-on-paper.md), which concerns surfaces with irreducible P (§8).*

Terms. X is the double plane w² = l₁⋯l₆ over F_p, P the polynomial of Frobenius on its part of rank 6, and the *sign* of X is the sign of the functional equation of P. For an ordered pair of lines, C is the curve of note 14 with its automorphism J, χ the quartic of Frobenius on the i-eigenspace, χ₂ that of Frobenius². The pair is of *commuting* or *inverting* type according to what Frobenius does to J (note 14, Lemma 4). A pair is *anomalous* if its type is not the one the sign asks for: commuting on a surface of sign −, inverting on a surface of sign +. A *collinear triple* is a splitting of the six lines in three pairs {a, b}, {c, d}, {e, f} such that the points L_a ∩ L_b, L_c ∩ L_d, L_e ∩ L_f lie on a line.

---

## 1. Result

1. **The fourfold of a surface of sign − [D, E].** Write P = (t² − p²)(t² − g t + p²)(t² − g′ t + p²). On the fourfold that belongs to X over F_p, Frobenius inverts J, and Frobenius² has on the i-eigenspace the polynomial

   > h(X) = X⁴ − (g g′/p) X³ + (g² + g′² − 2p²) X² − p g g′ X + p⁴.

   The polynomial of Frobenius over F_p on the whole fourfold is h(X²). So the fourfold is, up to isogeny, the Weil restriction to F_p of an abelian surface over F_{p²} with polynomial h; J exists only over F_{p²}. The passage between h and P is the second exterior power for a map that inverts J (Lemma B). It gives P from h up to P(−t), over F_p: the sign of g g′ is not visible over F_{p²}.
2. **Test on pairs of inverting type [E].** χ₂ was computed for 1,641 such pairs on surfaces of sign − (p = 7 to 37). In all 1,499 pairs that lie in at most one collinear triple, χ₂ = h, sign included. In the 142 pairs that lie in two collinear triples, χ₂ = h in 127 and χ₂ = h(−X) in 15 (listed in §5.2); there the curve has at least 16 automorphisms and h is a square.
3. **Where the anomalous pairs come from [D, E].** If the pair lies in a collinear triple, the curve has an involution γ with γJ = J⁻¹γ, so its automorphism group contains a dihedral group of order 8 and it has forms over F_p of both types. All 5,360 anomalous pairs among the 101,006 non-degenerate ordered pairs of the 3,990 configurations lie in exactly one collinear triple; none of the 85,646 pairs outside a collinear triple is anomalous. For an anomalous pair the form that belongs to the surface is the Prym variety of the twist of the curve by γ. Note 18 §4.3 states this.
4. **Test on the 2,524 pairs of commuting type on sign − [D, E].** For these the rule is: χ(X) = −¼·r((1 − εi)X) with ε = ±1 and r ∈ Z[X] such that r(X)·r(−X) = 16·h(X²/2). The quartics that this allows were computed from the surface alone and χ is one of them in 2,524 of 2,524. For the 1,144 pairs with c₁ ≠ 0, χ determines h with its sign and hence P up to P(−t): this is a test over F_p, and it has no exception. For the 1,380 pairs with c₁ = 0, h is a square, χ is even and it fits h and h(−X) equally: no more is tested than in note 18.
5. **The algebra [D, E].** In 822 of the 1,144 pairs the Prym variety of the curve is simple over F_p and over F_{p²}, with a commutative endomorphism algebra of degree 8; it becomes isogenous to a square only over F_{p⁴}, and γ is defined over F_{p⁴}. Frobenius has to be twisted by γ, an element of order 2, and the sign of the surface is the class, inner or outer, of the action of Frobenius on the dihedral group ⟨J, γ⟩ (§4.2). The curve's Frobenius π and the right one φ satisfy π² = J^(±1)·φ²: they differ by a square root of ±i, a root of unity of order 8, on the eigenspace, and by ±i on the second exterior power.
6. **Inverting pairs on surfaces of sign + [D, E].** There P(±t) is the second exterior power of a quartic f ∈ Z[X], the right form is of commuting type with χ = f up to X ↦ iX, and the anomalous curve has χ₂ = the quartic of the squares of the roots of f. 340 of 340 pairs computed agree.
7. **p ≡ 1 mod 4 [D, E].** For an inverting Frobenius the traces of odd powers on the eigenspace vanish identically, so the descent of the curve to F_p (which exists, as in Lemma 4) carries nothing beyond Frobenius². χ₂ comes from the sums of mode A over F_{p²} and F_{p⁴}; two traces suffice because χ₂ is rational and self-reciprocal, so a pair costs p⁴ operations and not p⁶ (§6).
8. **Why reducible P is frequent for sign + [E, H].** Two causes. A collinear triple forces a root ±p of P; 1,070 of the 1,564 configurations with sign + and reducible P have one, and none of the 534 with irreducible P. Of the 494 others, all have a root ±p because one of the two integers μ, ν of note 13 is 0, and these integers are at most 3 in absolute value at p = 37 (§7).
9. **The main case is not touched.** No configuration with sign + and irreducible P has a collinear triple (0 of 534), none of its 14,394 non-degenerate pairs is anomalous, and none of the rules above applies to it. The statements of note 14 for that case stand as they are.

What stays open is in §9.

---

## 2. Frobenius on the fourfold and on the surface [D]

Notation. A is the fourfold (the Prym variety of note 14, or a form of it), V its first cohomology (ℓ-adic, or rational for a lift with complex multiplication), regarded as a vector space of dimension 4 over Q(i) through J. H is the hermitian form made from the polarisation, H(x, y) = E(x, Jy) + i·E(x, y). The second exterior power Λ²V over Q(i) has dimension 6; an antilinear involution ⋆ of it is defined by ω′ ∧ ⋆ω = H(ω, ω′)·vol, up to a constant of absolute value 1; T is the fixed part of ⋆, a rational form of Λ²V, and iT is the fixed part of −⋆. The rank-6 part of the surface is T (Paranjape, Lemma 3; Lombardo–Peters–Schütt, (12)). A Frobenius of a form of A over F_p that respects the ring Z[i] either commutes with J (*linear*) or inverts it (*antilinear*).

### 2.1 Linear

π is a similitude, H(πx, πy) = p·H(x, y), and Λ²π∘⋆ = (det π / p²)·⋆∘Λ²π, where det is the determinant over Q(i), the constant term of χ. So:

- det π = p²: Λ²π maps T to T. Its polynomial on T is the second exterior power of χ. This is the case of note 14 §3, sign +.
- det π = −p²: Λ²π exchanges T and iT. No multiple c·Λ²π with c a root of unity of Q(i) other than ±i maps T to T, and i·Λ²π has determinant +p⁶ on T. So a form with constant term −p² is not the form of a surface of sign −: determinant −1 on T needs a Frobenius that inverts J.

### 2.2 Antilinear: the second exterior power of an antilinear Frobenius

Let φ be antilinear. Then E(φx, Jφy) = −E(φx, φJy) = −p·E(x, Jy), so

> H(φx, φy) = −p · conj H(x, y).

Π = φ² is linear, and if Πv = av then Π(φv) = ā·φv. So the characteristic polynomial h of Π over Q(i) has rational coefficients, and its roots are a, ā, a′, ā′ with |a| = |a′| = p:

> h(X) = X⁴ − A X³ + B X² − p²A X + p⁴, A = 2 Re a + 2 Re a′, B = 2p² + 4 Re a · Re a′.

Since φ anticommutes with J, its characteristic polynomial on V over Q is even, and its square has polynomial h² there; so it is h(X²).

**Lemma B.** Assume a, ā, a′, ā′ distinct. Then Λ²φ commutes with exactly two of the involutions θ⋆, a pair ±⋆₀, with fixed parts T₀ and iT₀. On T₀ its characteristic polynomial is

> P(t) = (t² − p²)(t² − g t + p²)(t² − g′ t + p²), with g g′ = p·A and g² + g′² = B + 2p²,

and on iT₀ it is P(−t). Conversely P determines h:

> h(X) = X⁴ − (g g′/p) X³ + (g² + g′² − 2p²) X² − p g g′ X + p⁴.

In words: if γ, γ′ are roots of the two quadratic factors, the roots of h are γγ′/p, γγ̄′/p and their conjugates.

*Proof.* Take v with Πv = av and w = φv, so φw = av; the same with v′, w′ for a′. H is diagonal in the basis v, w, v′, w′ because Π is a similitude with multiplier p² and a·conj(b) = p² only for b = a among the four eigenvalues. Put H(v, v) = η, H(v′, v′) = η′; then H(w, w) = −pη and H(w′, w′) = −pη′, so the signature is (2, 2) whatever the signs of η, η′. The star maps e_j ∧ e_k to z·H(e_j, e_j)H(e_k, e_k)·(±e_l ∧ e_m), antilinearly, with |z|² = 1/(p²η²η′²). Λ²φ maps v∧w to −a·v∧w, v∧v′ to w∧w′, w∧w′ to aa′·v∧v′, v∧w′ to a′·w∧v′ and w∧v′ to a·v∧w′, antilinearly. Comparing Λ²φ∘⋆ and ⋆∘Λ²φ on v∧w gives z/z̄ = aa′/p², which fixes z up to sign: z = ±β/(p²|ηη′|) with β² = aa′. Then the fixed part of ⋆ meets ⟨v∧w, v′∧w′⟩ in a plane on which Λ²φ is x ↦ −a·x̄, with eigenvalues p and −p; it meets ⟨v∧v′, w∧w′⟩ in a plane on which Λ²φ is multiplication by ∓sβ (s the sign of ηη′), and ⟨v∧w′, w∧v′⟩ in a plane on which it is multiplication by ∓sβ′ with β′ = aβ̄/p, so β′² = aā′ and ββ′ = p·a. So γ = ∓sβ, γ′ = ∓sβ′ and γγ′ = p·a, γγ̄′ = p·a′. ∎

- The relation is a statement over F_p. Over F_{p²} one sees only the squares of the roots of P, which determine g² and g′² and so P up to the two changes g ↦ −g and g′ ↦ −g′ separately; h(−X) belongs to the polynomial with g′ replaced by −g′. The sign of the odd coefficients of h is the information that F_{p²} does not have.
- The proof is a computation with signs, so it is also checked numerically: `sixlines_forms_fp.py selftest` builds 200 random antilinear maps of this kind, computes the real form they preserve and the polynomial on it, and finds the relation in all, to 10⁻¹⁰ [E].
- X and its quadratic twist have the same h: φ and Jφ have the same square.

### 2.3 Which abelian variety

By Tate's theorem [M] an abelian variety over F_p is determined up to isogeny by the polynomial of Frobenius on V. For an antilinear Frobenius this is h(X²), the polynomial of the Weil restriction from F_{p²} to F_p of an abelian surface B over F_{p²} with polynomial h. Over F_{p²} the Weil restriction is B × B^(p), the two factors are isogenous because h is rational, and J is an automorphism that exchanges them up to isogeny; it is defined over F_{p²} and not over F_p. So:

> **a six-line surface of sign − has, as its fourfold over F_p, the Weil restriction of an abelian surface over F_{p²} (up to isogeny), with h given by Lemma B.** Over F_{p²} the fourfold is isogenous to B², and its endomorphism algebra contains M₂(Q(Π)).

That h is the polynomial of an abelian surface is Honda–Tate [M]; I did not check the condition at p for the non-ordinary cases.

### 2.4 Forms when there is an antilinear involution

Suppose the endomorphism algebra of A over the algebraic closure contains γ with γJ = J⁻¹γ and γ² = 1, so that Q(i) and γ generate M₂(Q), and suppose Frobenius normalises the group D = ⟨J, γ⟩, dihedral of order 8 (§3: this is the case when D is the automorphism group of the curve). By Skolem–Noether a Frobenius is n·π₀ with π₀ commuting with J and γ, and n in the normaliser of Q(i) in GL₂(Q), of finite order modulo scalars. Up to powers of J and rational scalars there are four cases:

| | n | type | what it satisfies | Frobenius acts on D by |
|---|---|---|---|---|
| L⁺ | 1 | linear | det = p²; χ has rational coefficients up to X ↦ iX | an inner automorphism |
| A⁺ | γ | antilinear | φ² = π₀² | an inner automorphism |
| L⁻ | 1 + εJ | linear | det = −p²; π² = 2εJ·π₀² | the outer class: γ ↦ γJ^(±1) |
| A⁻ | γ(1 + J) | antilinear | φ² = 2π₀² | the outer class |

L⁺ and A⁺ are twists of each other by γ (the cocycle that sends Frobenius to γ), and so are L⁻ and A⁻; no twist by an element of D passes from one of these pairs of rows to the other. In the cases L⁻, A⁻ put σ = 2π₀, with characteristic polynomial r ∈ Z[X] and roots of absolute value √(2p). Then:

- A⁻: h is the polynomial of σ²/2, that is r(X)·r(−X) = 16·h(X²/2).
- L⁻: π = (1 + εJ)σ/2, so χ(X) = −¼·r((1 − εi)X). Hence c₃ = −p·c̄₁, c₂ is purely imaginary, and c₁ = y(ε + i) with y real: the shape that note 18 §4.2 found, with one more relation (|Re c₁| = |Im c₁|).
- The ratio of the two Frobenius maps of one pair of rows is γ up to a power of J. For L⁻ and A⁻: π⁻¹γπ = γJ^(±1), so γ is defined over F_{p⁴} and not over F_{p²}, and φ = πγ has φ² = π²·J^(∓1). So π and φ differ by a square root of ±i on V, a root of unity of order 8, and by ±i on the second exterior power: this is what note 18 §4.3 saw.

Which rows can be the form of a surface is decided by the real form T. For a linear Frobenius this is §2.1: L⁺ can, L⁻ cannot. For an antilinear one it depends on γ: if Λ²γ∘⋆ = u·⋆∘Λ²γ, then for n = γw with w in Q(i) the unit is u·(w̄/w)², which is u for w = 1 and −u for w = 1 + i. So A⁺ maps T to T if u = 1 and A⁻ does if u = −1; §3, Lemma F, gives u = −1.

---

## 3. The involution: three collinear points [D, E]

Fix the pair (L5, L6) and the cross-ratio map λ = N/D of note 14 §2.2. The fibre of λ over 0 consists of the lines of the pencil through L1 ∩ L3 and through L2 ∩ L4; over ∞, through L1 ∩ L4 and L2 ∩ L3; over 1, through L1 ∩ L2 and L3 ∩ L4.

**Lemma D.** 0 is a branch value of λ exactly when L1 ∩ L3, L2 ∩ L4 and L5 ∩ L6 are collinear; likewise 1 and ∞ for the other two ways to split {1, 2, 3, 4}. λ₅ and λ₆ are never 0, 1, ∞ for lines in general position. So a branch value b₁ or b₂ is one of 0, 1, ∞ exactly when the pair lies in a *collinear triple*: a splitting of the six lines in three pairs whose three intersection points are collinear. Then the other branch value is rational too, since Δ has rational coefficients.

**Lemma E.** Suppose 0 is a branch value. Then the point a of E over λ = 0 has order 2, and ψ(x) = a + ρ − x maps the branch set B = {Q, a + ρ − Q, b + ρ − Q, c + ρ − Q} of note 14, Lemma 3, to itself (2Q = a + b + c + ρ). The divisor of F·(F∘ψ) is 4 times a principal divisor, F·(F∘ψ) = c₀·G⁴ with a constant c₀, and

> γ : (x, v) ↦ (ψ(x), c₀^(1/4)·G(x)/v)

is an automorphism of C over the algebraic closure with γJγ⁻¹ = J⁻¹ and γ² = 1. So the automorphism group of C contains the dihedral group D = ⟨J, γ⟩ of order 8. The same holds with b for the value 1 and c for ∞.

*Proof.* ψ exchanges Q with a + ρ − Q, and b + ρ − Q with a − b + Q = c + ρ − Q (use 2Q = a + b + c + ρ and 2a = 0). The divisor of F is 2B − 4(c₄ + ρ) − 4(o), that of F∘ψ is 2B − 4(a − c₄) − 4(a + ρ), and B − (c₄ + ρ) − (o) − (a − c₄) − (a + ρ) has degree 0 and sum −2a = 0. For γ² = 1: γ² lies over the identity of E, so it is a power of J. ψ has four fixed points, none of them in B when the four points of B are distinct. Over a fixed point x₀ the cover is unramified; if y lies over x₀ and γ(y) = J^m y, then γ²(y) = γ(J^m y) = J^(−m)γ(y) = y; a power of J with a fixed point over x₀ is the identity. ∎

So the curves of pairs in a collinear triple have forms over F_p of both types, and the remark after Lemma 4 in note 18 §2, that for a curve with no automorphism beyond the powers of J the type is the same for every form over F_p, does not apply to them.

[E] `sixlines_forms_fp.py types`, all 3,990 configurations, all 15 splittings of each and all 30 ordered pairs:

- The pair lies in a collinear triple ⟺ one of 0, 1, ∞ is a root of the cubic of E: 15,360 non-degenerate ordered pairs with both, 85,646 with neither, none with one of the two.
- A splitting is a collinear triple ⟺ exchanging the two lines in each of its three pairs gives the associated configuration (the Gale transform of note 2, with its labels): 3,385 splittings with both, 56,465 with neither, none with one of the two. This is the classical picture: association is an involution of the moduli space of six ordered lines, and the 15 loci where three intersection points are collinear are where it is undone by a relabelling [E here; I did not look for the statement in a source].

**Lemma F.** For the involution γ of Lemma E, Λ²γ exchanges T and iT: Λ²γ∘⋆ = −⋆∘Λ²γ.

*Proof (a sketch that uses Paranjape's description as quoted in notes 14 and 18).* In his construction W = (C × C)/G maps to Pic²(E′) by the sum, the surface is W/ι̃ for a lift ι̃ of the involution ι of S, and T is the part of Λ²V ⊕ conj(Λ²V) fixed by ι̃, so ⋆ is ι̃ followed by conjugation. The lifts of ι are ι̃j^k with j = (J, 1), which acts on Λ²V as i. ι̃ and ι̃j² lie over the same involution D ↦ D₀ − D of Pic²(E′), whose four fixed points lie over λ₅ and λ₆, and they are told apart by which two of the four fibres they fix: those over L5 and L6 for ι̃, those over the conjugate lines L5′ and L6′ for ι̃j². Now γ̂ = (γ, γ) is an automorphism of W. On S it commutes with ι (both lie over e ↦ 2ρ − e, and on a fibre their composite is the Möbius map that exchanges the four marked points in pairs, in either order). On Pic²(E′) it is D ↦ D₁ − D with D₁ the class of 2x′ + 2γ(x′) for any point x′ of E′; the class x′ + γ(x′) is a fixed point, and its norm is the point a + ρ, which lies over the branch value 0 and not over λ₅ or λ₆. So D₁ = D₀ + κ, with κ the point of order 2 that generates the kernel of E′ → E, and γ̂ moves the four fixed points of D ↦ D₀ − D by κ: it carries the fibres over L5, L6 to those over L5′, L6′. Hence γ̂ι̃γ̂⁻¹ = ι̃j², which acts on Λ²V ⊕ conj(Λ²V) as −ι̃. ∎

- A pair can lie in two collinear triples (both branch values among 0, 1, ∞). Then C has two such involutions, their product is an automorphism of order dividing 4 over the translation by τ that commutes with J, and the group has order at least 16 [D].

---

## 4. Which form belongs to the surface

### 4.1 Outside the collinear triples [E]

For a pair that lies in no collinear triple the type follows the sign: 85,646 non-degenerate ordered pairs on the 3,990 configurations, none anomalous. This is what the remark of note 18 §2 predicts for a curve whose only automorphisms are the powers of J, if the surface determines the curve with its forms; and it holds on every configuration, including the 2,052 that have a collinear triple elsewhere.

### 4.2 In a collinear triple [D, E]

By Lemma E the curve has the forms of §2.4, if D is its whole automorphism group, and by §2.1, §2.4 and Lemma F only two of the four rows can be the form of a surface: L⁺, for sign +, and A⁻, for sign −. The class of the action of Frobenius on D, inner or outer, is the same for all forms of one curve. So [D]:

> **the sign of the surface is the class of the action of Frobenius on D: inner for sign +, outer for sign −.** For sign + the forms of the curve are L⁺ (right) and A⁺ (anomalous); for sign − they are A⁻ (right) and L⁻ (anomalous). In both cases the right form is the anomalous one twisted by γ.

[E] The data agree with each part of this:

- **Sign +.** Where the curve of note 14 is of commuting type it gives P or P(−t) (note 14). Where it is of inverting type, χ₂ is the quartic of the squares of the roots of f, as A⁺ requires and A⁻ does not (§5.3).
- **Sign −.** In all 982 configurations of sign − with a collinear triple there is r ∈ Z[X] with r(X)·r(−X) = 16·h(X²/2): twice Frobenius² is a square in Q[X]/h, as A⁻ requires. In 477 of them Frobenius² itself is not a square (no r₀ ∈ Z[X] with r₀(X)·r₀(−X) = h(X²)), which excludes A⁺; in the other 505 both hold and the polynomial does not decide. Every pair of commuting type on these surfaces is L⁻ (§5.1). For comparison, among the 910 configurations of sign − without a collinear triple all four combinations occur: neither is a square in 420, only twice Frobenius² in 98, only Frobenius² in 346, both in 46.

Which of the two forms the recipe of note 14 writes down is decided by Lemma 4 there, that is by whether c₄ is rational; it is not a property of the surface. On one surface and one collinear triple the three pairs can be of different types: at p = 13, configuration (2, 3, 5, 6), sign −, the triple {1, 5}, {2, 4}, {3, 6} has the pair (1, 5) of commuting type and (3, 6) of inverting type [E]. All anomalous pairs lie in exactly one collinear triple; the 1,420 non-degenerate pairs that lie in two are all of the type of the sign.

---

## 5. The tests [E]

All from `sixlines_forms_fp.py`; the output is `experiments/sixlines_forms_fp.log`.

### 5.1 Pairs of commuting type on surfaces of sign −

The 2,524 pairs on 673 configurations stored by note 14 (`commuting`). For each, h is computed from (e₁, e₂) of the surface, then all r ∈ Z[X] of the shape X⁴ − aX³ + bX² − 2paX + 4p² with r(X)·r(−X) = 16·h(X²/2), then the quartics −¼·r((1 − εi)X) for ε = ±1. The curve side enters only in the comparison.

| | pairs | χ is one of the predicted quartics | r ∈ Z[X] for | h from χ gives |
|---|---|---|---|---|
| c₁ ≠ 0 | 1,144 | 1,144 | exactly one ε | P or P(−t), and no other polynomial: 1,144 |
| c₁ = 0 | 1,380 | 1,380 | both ε | P or P(−t) for one ε, another polynomial of sign − for the other: 1,380 |

- No exception. In all 2,524: c₃ = −p·c̄₁, c₂ purely imaginary (note 18), and also |Re c₁| = |Im c₁|, the pair lies in exactly one collinear triple, and a, b, c are not all rational.
- The number of predicted quartics is 4 for 952 of the pairs with c₁ ≠ 0 (the four χ(i^k X), which are the quartic twists of the curve), 6 or 8 for the other 192; for c₁ = 0 it is 2 or 6.
- What is tested beyond note 18. Note 18 compared the second exterior power of χ₂ with the polynomial of the twisted surface over F_{p²}; that fixes h up to h(−X). Here the sign is fixed too: for c₁ ≠ 0 the one ε for which r is integral gives h and not h(−X), 1,144 times out of 1,144. The polynomial with g′ replaced by −g′ is excluded each time.
- What is not tested: for c₁ = 0 (h a square, χ even) the curve does not distinguish h from h(−X). And the curve does not see which of X and its quadratic twist it belongs to: the sign δ with Σσ_jσ_k/2 = δ·e₁ (sum over the two pairs of roots of r that give γ and γ′) is +1 in 512 pairs, −1 in 502, and e₁ = 0 in 130. That is as it should be, since X and its twist have the same h.
- Structure, for the 1,144: χ·χ̄ is irreducible over Q in 822 and so is χ₂·χ̄₂, while χ₄ (Frobenius⁴) is rational and irreducible of degree 4 in 822. So there the Prym variety is simple over F_p and F_{p²} with endomorphism algebra the field Q[X]/(χ·χ̄) = K′(i), K′ = Q[X]/h, and over F_{p⁴} it is isogenous to the square of an abelian surface with endomorphism algebra M₂(K′) [D from Tate's theorem, M]. In the other 322, h is reducible (degrees 2, 2 in 192, and 1, 1, 2 in 130) and the fourfold already splits over F_p. So the explanation of note 14 §5.2, that the fourfold is not simple over F_p and has other forms, is not the one that applies in these pairs (as note 18 §4.3 says): over F_{p²} it is simple in most cases, and what matters is the automorphism γ over F_{p⁴}.

### 5.2 Pairs of inverting type on surfaces of sign −

`inverting`: χ₂ from two traces (§6), compared with h from the surface. One accepted outcome: c₁(χ₂) = g g′/p and c₂(χ₂) = g² + g′² − 2p², two integers of size up to 4p and 6p².

| p | pairs computed | in no collinear triple, χ₂ = h | in one, χ₂ = h | in two: χ₂ = h / χ₂ = h(−X) | sign +, χ₂ as predicted (§5.3) | pairs taken |
|---|---|---|---|---|---|---|
| 7 | 10 | 0 of 0 | 0 of 0 | 10 / 0 | 0 of 0 | all |
| 11 | 104 | 50 of 50 | 12 of 12 | 38 / 4 | 8 of 8 | all |
| 13 | 198 | 134 of 134 | 44 of 44 | 16 / 4 | 50 of 50 | all |
| 17 | 702 | 518 of 518 | 144 of 144 | 34 / 6 | 118 of 118 | all |
| 19 | 312 | 173 of 173 | 113 of 113 | 26 / 0 | 52 of 52 | at most 6 per configuration |
| 23 | 201 | 131 of 131 | 66 of 66 | 3 / 1 | 66 of 66 | 2 per configuration, first 150 configurations |
| 29 | 90 | 66 of 66 | 24 of 24 | 0 / 0 | 30 of 30 | 2 per configuration, first 60 |
| 37 | 24 | 14 of 14 | 10 of 10 | 0 / 0 | 16 of 16 | 2 per configuration, first 20 |
| all | 1,641 | 1,086 of 1,086 | 413 of 413 | 127 / 15 | 340 of 340 | |

"Pairs computed" counts the pairs on surfaces of sign −; the 17 pairs of note 18 §4.2 are among those of p = 7 and 11.

- **In at most one collinear triple: 1,499 pairs, all with χ₂ = h.** 1,086 of them lie in no collinear triple and are the ordinary case for sign −; 413 lie in one.
- **In two collinear triples: 142 pairs, 127 with χ₂ = h and 15 with χ₂ = h(−X).** The exceptions:

  | p | configuration | P | h | ordered pairs with χ₂ = h(−X) |
  |---|---|---|---|---|
  | 11 | (2, 3, 8, 9) | (t² − 121)(t + 11)²(t² − 2t + 121) | (X² + 2X + 121)² | (1, 3), (3, 1), (2, 4), (4, 2) |
  | 13 | (2, 3, 5, 6) | (t² − 169)(t − 13)²(t² − 2t + 169) | (X² − 2X + 169)² | (1, 3), (3, 1), (2, 4), (4, 2) |
  | 17 | (2, 3, 3, 14) | (t² − 289)(t − 17)²(t² + 14t + 289) | (X² + 14X + 289)² | (2, 5), (5, 2) |
  | 17 | (2, 6, 6, 10) | (t² − 289)(t − 17)²(t² − 10t + 289) | (X² − 10X + 289)² | (1, 4), (4, 1), (2, 5), (5, 2) |
  | 23 | (2, 3, 7, 18) | (t² − 529)(t + 23)²(t² − 26t + 529) | (X² + 26X + 529)² | (1, 2) |

  In each of them h is the square of a quadratic polynomial, and P has the root p or −p with multiplicity 3. With two involutions the curve has an automorphism u over the translation by τ that commutes with J (§3), so there are more forms than those of §2.4, and twisting by u can change Frobenius² by −1 [D]. I did not work out which sign the recipe of note 14 produces in this case.
- χ₂ has rational coefficients in every pair computed, as an inverting Frobenius forces.
- The sign is the content of the test over F_p: h(−X) is the polynomial for g′ ↦ −g′, which has the same squares of roots. 720 of the 1,499 pairs have χ₂ without a repeated factor; for a χ₂ that is the square of a quadratic polynomial the sign is still tested unless that polynomial is even.

### 5.3 Pairs of inverting type on surfaces of sign +

These are all anomalous and all lie in exactly one collinear triple. For each surface the quartics f = X⁴ − aX³ + bX² − paX + p² in Z[X] with second exterior power P(t) or P(−t) are found by search, and χ₂ is compared with the quartic of the squares of the roots of f.

- 340 pairs computed (table of §5.2); 340 agree. The prediction is a single quartic for 196 of them and one of two for 144 (when both P(t) and P(−t) are second exterior powers).
- On all 1,564 configurations of sign + with reducible P, with or without a collinear triple, such an f exists; on none of the 534 with irreducible P.

### 5.4 The involution, directly

`involution`: for every pair in a collinear triple and the point of order 2 among a, b, c, the value F(x)·F(ψ(x)) is computed at about 140 random points of E(F_{p²}); its class modulo fourth powers is constant in 15,973 of 15,973 cases (p = 11 to 37; a pair in two collinear triples counts twice). With one of the other two points in place of the point of order 2 it is constant in 0 of 28,034. So F∘ψ is a constant times G⁴/F, which is Lemma E.

---

## 6. Inverting pairs for p ≡ 1 mod 4 [D, E]

Note 18 computed inverting pairs only for p ≡ 3 mod 4, where c₄ is rational and the curve v⁴ = F is defined over F_p. For p ≡ 1 mod 4 an inverting pair has c₄ conjugate to −c₄.

- **The descent.** It is that of Lemma 4, mode B, unchanged: F·Frob(F) = G⁴ with G over F_p, and Φ : (x, v) ↦ (Frob x, G(Frob x)/v^p) is a Frobenius for a form of the curve over F_p. With i in F_p, J(v) = iv gives Φ∘J = J⁻¹∘Φ: inverting.
- **It carries no information on the eigenspace.** For x in E(F_p) and v⁴ = F(x): (v^(p+1))⁴ = F(x)·Frob(F)(x) = G(x)⁴, and (iv)^(p+1) = −v^(p+1) because p ≡ 1 mod 4. So v^(p+1)/G(x) takes two opposite values, each on two of the four points, and J^m∘Φ fixes two of the four points for two values of m that differ by 2, and none for the others. The contribution of x to the trace on the i-eigenspace, Σ_m i^(−m)·(number of fixed points of J^m Φ over x), is 0. The same holds for every odd power of Φ. This is as it must be: Φ exchanges the eigenspaces of J. So there is no analogue of the odd traces of mode B.
- **Frobenius².** Φ² is the Frobenius over F_{p²} of v⁴ = F, with F defined over F_{p²}. Its traces are the sums of mode A over F_{p²}, F_{p⁴}, F_{p⁶}. `pr_trace` computes them correctly once the square root ω of the least primitive root is available in the field; it sets ω = 0 when a, b, c are rational, because note 14 never needed c₄ outside F_p in that case. The new script uses a copy `fp_trace` in which ω is always set, and nothing else is changed.
- **Two traces suffice.** χ₂ has rational coefficients and its roots are closed under a ↦ p²/a, so χ₂ = X⁴ − c₁X³ + c₂X² − p²c₁X + p⁴ and the third trace is determined. A pair then costs p⁴ operations instead of p⁶. The third trace was computed as a check on 35 pairs (10 at p = 7, 13 at p = 11, 12 at p = 13, the last with p ≡ 1 mod 4): c₃ = p²·c₁ in all 35.
- **How far.** All inverting pairs for p = 7, 11, 13, 17; for p = 19, up to 6 pairs per configuration; for p = 23, 2 pairs for each of the first 150 configurations that have one; for p = 29 the first 60 and for p = 37 the first 20. Where the number of pairs is capped, half of them are taken among the pairs in exactly one collinear triple. p = 31 was not run. The runs took 95 minutes with 4 processes; a pair costs about 3.5 seconds at p = 17 and 100 seconds at p = 37 with two traces, and about 4 minutes at p = 13 with three.

---

## 7. Why P is so often reducible for sign + [E, H]

Note 18 §4.3 records that 37% to 58% of all configurations have sign + and reducible P and gives one ingredient. Two causes account for all 1,564.

**A collinear triple (1,070 configurations).** [D] With γ, the fourfold is isogenous over the algebraic closure to the square of an abelian surface, and for sign + (form L⁺) χ is rational up to X ↦ iX. The roots of a rational χ come in pairs with product p, so its second exterior power has the root p twice, and P has the root p or −p. [E] All 1,070 configurations with sign + and a collinear triple have reducible P, and none of the 534 with irreducible P has a collinear triple.

How many configurations have one [E]:

| p | configurations | with a collinear triple | collinear splittings (of 15 per configuration) | p × rate | 1 − (1 − rate)¹⁵ |
|---|---|---|---|---|---|
| 5 | 1 | 1 | 10 | 3.33 | 100% |
| 7 | 3 | 3 | 16 | 2.49 | 100% |
| 11 | 15 | 14 | 48 | 2.35 | 97% |
| 13 | 26 | 25 | 72 | 2.40 | 95% |
| 17 | 74 | 63 (85%) | 141 | 2.16 | 87% |
| 19 | 117 | 92 (79%) | 196 | 2.12 | 83% |
| 23 | 257 | 175 (68%) | 320 | 1.91 | 73% |
| 29 | 682 | 377 (55%) | 616 | 1.75 | 61% |
| 31 | 905 | 469 (52%) | 744 | 1.70 | 57% |
| 37 | 1,910 | 833 (44%) | 1,222 | 1.58 | 48% |
| all | 3,990 | 2,052 | 3,385 | | |

A collinear triple is one condition, so a given splitting is collinear at a rate of about 1/p; the measured rate is 1.6/p at p = 37 and falls with p (the count is of configurations in general position, which are a smaller share of all sextuples for small p). With 15 splittings the share of configurations with at least one is close to what independence gives, the last column. It tends to 15/p. So this cause is large for p ≤ 37 only because 15/p is.

**μν = 0 without a collinear triple (494 configurations).** For sign +, P(p) = 64p⁵μ² and P(−p) = 64p⁵ν² with integers μ, ν, and c₁ = 4(±ν ± μi) up to a unit (note 13, note 14 §1). [E] P is reducible exactly when μν = 0, on all 2,098 configurations of sign +. [H] |c₁| is of the size of √p, so μ and ν are integers of size √p/4, and one of them vanishes often. At p = 37 the 565 configurations of sign + without a collinear triple have (|μ|, |ν|) up to order:

| (0, 0) | (0, 1) | (0, 2) | (0, 3) | (1, 1) | (1, 2) | (1, 3) |
|---|---|---|---|---|---|---|
| 44 | 164 | 61 | 3 | 168 | 120 | 5 |

so μν = 0 in 272 of 565 (48%). The share should fall like a constant over √p; in this range it does not fall yet (2 of 7 at p = 17, 6 of 14 at p = 19, 17 of 45 at p = 23, 82 of 165 at p = 29, 114 of 230 at p = 31, 272 of 565 at p = 37). I have no model for the constant. The same small integers explain why c₁ = 0 in 1,380 of the 2,524 commuting pairs on sign −: there c₁ = y(ε + i) with y divisible by 4 and |y| ≤ 2√(2p), and the values found are 0, ±4, ±8, ±12.

**The anomalous pairs and the signs of X and X′.** Note 18 §4.2 observes that every anomalous pair has X and X′ of the same sign, that is, rational b₁, b₂. This is now a consequence [D, given §4.1]: an anomalous pair lies in a collinear triple, so one of b₁, b₂ is 0, 1 or ∞ and the other is rational (Lemma D). The converse is false: 46,334 non-degenerate pairs have rational b₁, b₂ and are not anomalous. The exact statement is: anomalous ⟹ in exactly one collinear triple, and among the 13,940 pairs that are in exactly one, the 5,360 anomalous ones are those for which Lemma 4 gives the type opposite to the sign.

---

## 8. The main case

Nothing found here changes a pair on a surface with sign + and irreducible P.

- [E] None of the 534 such configurations has a collinear triple, and none of their 14,394 non-degenerate pairs is anomalous.
- [D] A collinear triple forces a root ±p of P (§7), so this is not an accident of the list: for these surfaces the curve of every pair has, as far as the argument goes, no involution of the kind of Lemma E, and the forms of §2.4 do not exist.
- The rule of note 14 for this case (commuting type, determinant p², second exterior power P or P(−t)) is the row L⁺ of §2.4 without the other rows.

Nothing in note 14 is corrected, and what is found concerns only the anomalous pairs. The third bullet of note 18 §4.3 states it (the right form is the Prym variety of a form of the curve; the curve has more automorphisms than the powers of J; over F_{p²} the fourfold is mostly simple); what is corrected is the reading of note 18 §4.2 as a criterion (the same sign of X and X′ is necessary for an anomalous pair, not sufficient).

---

## 9. What this does not give

- **A construction of the correspondence over F_p.** The fourfold of a surface of sign − is identified by its Frobenius polynomial, h(X²), and by the type of its Frobenius. That the cycle of Paranjape is defined over F_p for the Weil restriction, or for the twist of the curve by γ, is not shown; nor is the surface W/ι̃ of his construction followed through the twist.
- **Lemma F** is a sketch on top of Paranjape's description of the two lifts as used in notes 14 and 18; the paper was not re-read for this note. Its consequence (sign − goes with the outer class) is observed on 982 configurations and forced by the polynomial in 477 of them.
- **That D is the whole automorphism group** of the curve of a pair in exactly one collinear triple, which §4.2 uses, is not shown.
- **Pairs in two collinear triples.** The rule holds up to X ↦ −X in χ₂ and I do not predict the sign (15 exceptions in 142 pairs).
- **Pairs with c₁ = 0** (1,380 of 2,524): no test over F_p beyond note 18.
- **Which type Lemma 4 gives for a pair in a collinear triple**, in terms of the surface. It is decided by the rationality of c₄; I did not look for a criterion on the six lines.
- **Honda–Tate and Tate's isogeny theorem**, used in §2.3 and §5.1, are quoted from memory. "Weil restriction of an abelian surface" is a statement up to isogeny; the polarisation and the isomorphism class were not examined.
- **The statement about association** in §3 is measured, not read in a source.
- **The 494 configurations of sign + with reducible P and no collinear triple**: only the heuristic of §7. Presumably they lie on other divisors of the moduli space where the fourfold acquires an antilinear endomorphism that is not an automorphism of the curve; I did not test this.
- **Inverting pairs beyond those of §6**: p = 31, most of p = 23, 29, 37, and the third trace beyond p = 13.
- Anything about the group action, or about Picard number below 16. The surfaces of this note all have Picard number at least 18 over the algebraic closure.

---

## Reproduce

From `experiments/`:

```bash
python -I sixlines_forms_fp.py selftest
python -I sixlines_forms_fp.py commuting sixlines_prym.jsonl
python -I sixlines_forms_fp.py types 5 7 11 13 17 19 23 29 31 37 --workers 4
python -I sixlines_forms_fp.py involution 11 13 17 19 23 29 31 37 --workers 4
python -I sixlines_forms_fp.py inverting 7 11 13 17 --pairs 0 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 19 --pairs 6 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 23 --pairs 2 --limit 150 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 29 --pairs 2 --limit 60 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 37 --pairs 2 --limit 20 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 7 --pairs 0 --k 3 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 11 --pairs 2 --k 3 --workers 4 --out sixlines_forms_fp.jsonl
python -I sixlines_forms_fp.py inverting 13 --pairs 1 --limit 12 --k 3 --workers 4 --out sixlines_forms_fp.jsonl
```

`sixlines_forms_fp.log` is the output of these commands in this order, and `sixlines_forms_fp.jsonl` has one record per run and inverting pair, 2,016 in all: 1,981 pairs with two traces, then 35 of them again with three (the file is appended to; delete it before a rerun). The first three take ten seconds together and `involution` a minute and a half; the inverting runs take 95 minutes with 4 processes (12 minutes for p = 17, 10 for p = 19, 16 for p = 23, 19 for p = 29, 18 for p = 37, 13 for the third trace at p = 13). They need `sixlines_prym.jsonl` (note 14) and `sixlines_prym.py`, which is imported and not changed. `involution` uses random points; its counts of points vary between runs, the counts of pairs do not. Everything else is deterministic.

## Sources

- Paranjape, *Abelian varieties associated to certain K3 surfaces*, Compositio Math. 68 (1988), 11–22: as read for notes 14 and 18 [P]; used here for the lattice T inside the second exterior power and for the construction that §3 and §4.2 refer to. Not re-read for this note.
- Lombardo, Peters, Schütt, *Abelian fourfolds of Weil type and certain K3 double planes*, https://arxiv.org/abs/1209.5997 : as quoted in note 18 §4.3 [P there]; not re-read.
- Recalled, not checked [M]: Tate's isogeny theorem; Honda–Tate; the Skolem–Noether theorem; the classical theory of association of six points in the plane (Coble).
- In this repository: notes [2](02-six-line-count.md), [13](13-hardness-on-paper.md), [14](14-fourfold-of-a-six-line-surface.md), [15](15-fourfold-map-and-the-group-action.md), [18](18-checks-of-notes-13-to-16.md).
