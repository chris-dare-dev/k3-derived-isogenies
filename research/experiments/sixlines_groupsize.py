#!/usr/bin/env python
"""
sixlines_groupsize.py -- size and structure of the group acting on six-line K3 surfaces with a
given Frobenius polynomial, as q grows (note 13).

For an ordinary six-line double plane over F_q with all lines rational and geometric Picard number 16,
the Frobenius polynomial on the transcendental part is

    P(t) = t^3 h(t + q^2/t),   h(x) = x^3 - e1 x^2 + (e2 - 3 q^2) x - (e3 - 2 q^2 e1),

with p not dividing e1, q | e2, q^2 | e3 (ordinary), E = Q[t]/P a CM field of degree 6 and E0 = Q[x]/h.
When the endomorphism ring of the transcendental lattice is the maximal order O_E, the surfaces with
that polynomial in one genus form a torsor under

    G = Cl(E) / < image of Cl(E0), classes of the primes of E ramified over E0 >

(note 13, Lemma 9; unmarked surfaces, split type).  The script measures |G| and its structure.

  real      the 44 polynomials of sixlines_cm.jsonl, which come from enumerated surfaces (p <= 37):
            number of surfaces per polynomial, and for the maximal-order ones G and the check that the
            number N of lattices is a multiple of |G|;
  all q     every polynomial with h(2q) = 4 q^2 m^2, h(-2q) = -4 q^2 n^2 (the shape found in all 44 real
            cases), e1 > 0, three real roots of h in (-2q, 2q), P irreducible, (F - 1)/2 integral.  All
            have m = 4 mu, n = 4 nu.  For those with R_min = O_E (always mu = nu = 1, asserted): G, and
            the number N of unimodular hermitian O_E-lattices of the right signature (the count of
            note 2), when the class number allows it;
  sample    random polynomials with mu = nu = 1 and R_min = O_E for larger q; G only, and N when the
            class number is at most --cap;
  fourfold  for the real polynomials and all polynomials of the given q: with a = 4 (nu + mu i),
                chi(X) = X^4 - a X^3 + e1 X^2 - q conj(a) X + q^2   over Z[i],
            asserts that P is the polynomial whose roots are the products of two roots of chi, and that
            chi conj(chi) is a Weil polynomial of weight 1 with integer coefficients (note 13, Proposition Y).

PARI's bnfinit is used without certification, so class groups are correct under GRH.

  python -I sixlines_groupsize.py [--out FILE] [--real] [--all q ...] [--sample q:count ...] [--fourfold q ...] [--cap H]

Run from experiments/.
"""
import argparse
import json
import math
import os
import random
import statistics
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402

pari(r"""
gs_group(bnf, q) = {
  \\ [h(E), h(E0), elementary divisors of G = Cl(E) / <Cl(E0), primes ramified over E0>, |d_E| / d_E0^2,
  \\  number of primes of E ramified over E0, cyc of Cl(E)]
  my(nf = bnf.nf, P = nf.pol, t = variable(P), F = Mod(t, P), s = F + q^2 / F, h = minpoly(s, 'x),
     bnf0 = bnfinit(h, 1), cyc = bnf.cyc, k = #cyc, ram = abs(nf.disc) / bnf0.disc^2, M, H, D, nram = 0,
     up = (I -> my(two = idealtwoelt(bnf0, I)); idealhnf(nf, two[1], subst(lift(nfbasistoalg(bnf0, two[2])), 'x, s))));
  M = matdiagonal(cyc);
  for(j = 1, #bnf0.gen, if(k, M = concat(M, Mat(bnfisprincipal(bnf, up(bnf0.gen[j]), 0)))));
  foreach(factor(ram)[, 1]~, l,
    foreach(idealprimedec(bnf0, l), pr,
      my(fa = idealfactor(nf, up(pr)));
      if(matsize(fa)[1] == 1 && fa[1, 2] == 2,
        nram++; if(k, M = concat(M, Mat(bnfisprincipal(bnf, fa[1, 1], 0)))))));
  if(k == 0, return([1, bnf0.no, [], ram, nram, []]));
  H = mathnf(M); D = select(d -> d > 1, matsnf(H));
  [bnf.no, bnf0.no, D, ram, nram, cyc];
}
""")
GP_group = pari("gs_group")
pari(r"""
gs_fourfold(P, q, e1, mu, nu) = {
  \\ chi = X^4 - a X^3 + e1 X^2 - q conj(a) X + q^2 over Z[i], a = 4 (nu + mu i).  Returns
  \\ [prod_{i<j} (t - x_i x_j) == P  for the roots x_i of chi,  chi * conj(chi) in Z[X],  its roots have absolute
  \\  value sqrt(q),  chi * conj(chi) is irreducible,  the middle coefficients show slopes 0 and 1 only]
  my(a = 4 * (nu + mu * I), ab = conj(a), chi = 'X^4 - a * 'X^3 + e1 * 'X^2 - q * ab * 'X + q^2,
     chib = 'X^4 - ab * 'X^3 + e1 * 'X^2 - q * a * 'X + q^2, Wy = chi * chib, W = real(Wy), t = variable(P),
     C = matrix(4, 4, r, c, if(r == c + 1, 1, if(c == 4, -polcoef(chi, r - 1), 0))), idx = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]], L, cp,
     pp = factor(q)[1, 1]);
  L = matrix(6, 6, r, c, C[idx[r][1], idx[c][1]] * C[idx[r][2], idx[c][2]] - C[idx[r][1], idx[c][2]] * C[idx[r][2], idx[c][1]]);
  cp = charpoly(L, t);
  [imag(cp) == 0 && real(cp) == P, imag(Wy) == 0, vecmax(apply(z -> abs(abs(z) - sqrt(q)), polroots(W))) < 1e-9 * sqrt(q),
   polisirreducible(W), polcoef(W, 4, 'X) % pp != 0];
}
""")
GP_fourfold = pari("gs_fourfold")


def poly_from(q, e1, m, n):
    """The coefficients e = (e1, e2, e3) with h(2q) = 4 q^2 m^2 and h(-2q) = -4 q^2 n^2."""
    return [e1, q * (m * m + n * n - q), -2 * q * q * (m * m - n * n + e1)]


def shape(q, e):
    """(h(2q) / q^2, -h(-2q) / q^2) if both are integers, else None."""
    e1, e2, e3 = e
    c2, c3 = e2 - 3 * q * q, e3 - 2 * q * q * e1
    a = 8 * q ** 3 - 4 * q * q * e1 + 2 * q * c2 - c3
    b = 8 * q ** 3 + 4 * q * q * e1 + 2 * q * c2 + c3
    if a % (q * q) or b % (q * q):
        return None
    return a // (q * q), b // (q * q)


def is_square(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def admissible(q, p, e):
    """Ordinary, three real roots of h in (-2q, 2q), (F - 1)/2 integral, P irreducible."""
    e1, e2, e3 = e
    if e1 % p == 0:
        return None
    co = cm.charpoly(q, e)
    if not cm.two_adic_ok(co):
        return None
    h = pari(f"x^3 - ({e1})*x^2 + ({e2 - 3 * q * q})*x - ({e3 - 2 * q * q * e1})")
    if int(pari.polsturm(h, [-2 * q, 2 * q])) != 3:
        return None
    P = cm.pol(co)
    if not pari.polisirreducible(P):
        return None
    return P


def fourfold(q, e):
    """The Weil polynomial of the abelian fourfold with multiplication by Z[i] attached to P (note 13, Proposition Y)."""
    a, b = shape(q, e)
    mu, nu = math.isqrt(a) // 8, math.isqrt(b) // 8
    res = [bool(int(x)) for x in GP_fourfold(cm.pol(cm.charpoly(q, e)), q, e[0], mu, nu)]
    assert res[0], "P is not the second exterior power over Z[i] of the quartic"
    assert res[1] and res[2], "the degree-8 polynomial is not a Weil polynomial of weight 1"
    return {"irreducible": res[3], "ordinary": res[4]}


def run_fourfold(qs):
    rows = [json.loads(line) for line in open("sixlines_cm.jsonl")]
    stats = Counter()
    for r in rows:
        f = fourfold(r["p"], r["e"])
        stats["real"] += 1
        stats["real, degree-8 polynomial irreducible"] += f["irreducible"]
        stats["real, fourfold ordinary"] += f["ordinary"]
    for q in qs:
        p = prime_of(q)
        for e in candidates_all(q, p):
            if admissible(q, p, e) is None:
                continue
            f = fourfold(q, e)
            stats["enumerated"] += 1
            stats["enumerated, degree-8 polynomial irreducible"] += f["irreducible"]
            stats["enumerated, fourfold ordinary"] += f["ordinary"]
    print("fourfold:", dict(stats), flush=True)


def analyse(q, p, e, cap, want_count=True):
    """G for one polynomial; None if R_min is not the maximal order."""
    P = admissible(q, p, e)
    if P is None:
        return None
    a, b = shape(q, e)
    mu_nu = [math.isqrt(a) // 8, math.isqrt(b) // 8]
    assert is_square(a) and is_square(b) and a % 64 == 0 and b % 64 == 0, "h(2q), h(-2q) are not 64 q^2 times squares"
    nf = pari.nfinit(P)
    gens = cm.GP_sl_gens(nf, q)
    if not all(bool(cm.GP_sl_integral(nf, x)) for x in gens):
        return {"q": q, "e": list(e), "maximal": False, "integral": False, "mu_nu": mu_nu}
    if int(cm.GP_sl_index(nf, gens)) != 1:
        return {"q": q, "e": list(e), "maximal": False, "mu_nu": mu_nu}
    assert mu_nu == [1, 1], "a maximal order with mu or nu > 1"
    bnf = None
    for prec in (38, 200):                              # PARI: "bug in small_norm (precision too low)" for some fields
        try:
            pari(f"default(realprecision, {prec})")
            bnf = pari.bnfinit(pari.nfinit(P), 1)
            hE, h0, D, ram, r2, cyc = GP_group(bnf, q)
            break
        except Exception:
            bnf = None
        finally:
            pari("default(realprecision, 38)")
    if bnf is None:
        return {"q": q, "e": list(e), "maximal": True, "mu_nu": mu_nu, "failed": True}
    D = [int(d) for d in D]
    G = math.prod(D)
    rec = {"q": q, "e": list(e), "maximal": True, "mu_nu": mu_nu, "h": int(hE), "h0": int(h0), "G": G, "G_cyc": D,
           "ramification": int(ram), "ramified_primes": int(r2), "log_q_G": round(math.log(G) / math.log(q), 4) if G > 1 else 0.0}
    if want_count and int(hE) <= cap:
        _, cnt, unit_index, w, _, _, _, _, _ = cm.GP_sl_count(bnf, q)
        N = [int(c) * int(unit_index) for c in cnt]
        rec["N"] = N
        assert all(x % G == 0 for x in N), "the number of lattices is not a multiple of |G|"
    return rec


def largest_prime(n):
    return max(int(x) for x in pari.factor(n)[0]) if n > 1 else 1


def summary(q, recs):
    failed = sum(1 for r in recs if r.get("failed"))
    mx = [r for r in recs if r.get("maximal") and not r.get("failed")]
    if not mx:
        return f"q={q}: no polynomial with maximal order"
    counted = [r for r in mx if "N" in r]
    occ = [r for r in counted if max(r["N"]) > 0]
    use = occ if len(counted) == len(mx) and occ else mx      # all maximal ones unless every one was counted
    Gs = sorted(r["G"] for r in use)
    cyclic = sum(len(r["G_cyc"]) <= 1 for r in use)
    oddcyc = sum(sum(1 for d in r["G_cyc"] if d >> (d & -d).bit_length() - 1 > 1) <= 1 for r in use)
    tworank = Counter(sum(d % 2 == 0 for d in r["G_cyc"]) for r in use)
    even = sum(r["G"] % 2 == 0 for r in use)
    unram = sum(r["ramification"] & (r["ramification"] - 1) == 0 for r in use)
    lp = statistics.median(math.log(largest_prime(r["G"])) / math.log(r["G"]) for r in use if r["G"] > 1) if any(
        r["G"] > 1 for r in use) else 0
    return (f"q={q}: {len(recs)} admissible, {sum(1 for r in recs if r.get("integral", True))} with integral generators, {len(mx)} maximal, "
            f"{len(occ)} with a lattice (of {len(counted)} counted); over the {len(use)} "
            f"{'with a lattice' if use is occ else 'maximal'}: |G| min {Gs[0]} median {statistics.median(Gs)} max {Gs[-1]}, "
            f"median log_q|G| {statistics.median(r['log_q_G'] for r in use):.3f}, median |G|/q^2 "
            f"{statistics.median(r['G'] / q ** 2 for r in use):.4f}, cyclic {cyclic}, odd part cyclic {oddcyc}, "
            f"2-ranks {dict(sorted(tworank.items()))}, bnfinit failed {failed}, |G| even {even}, "
            f"E/E0 unramified at odd primes {unram}, median log(largest prime)/log|G| {lp:.2f}")


def prime_of(q):
    return int(pari.factor(q)[0][0])


def run_real(fh, cap):
    rows = [json.loads(line) for line in open("sixlines_cm.jsonl")]
    print(f"real: {len(rows)} polynomials from enumerated surfaces")
    by = {}
    for r in rows:
        by.setdefault(r["p"], []).append(r)
    for q in sorted(by):
        rs = by[q]
        for r in rs:
            s = shape(q, r["e"])
            assert s and is_square(s[0]) and is_square(s[1]) and s[0] % 4 == 0 and s[1] % 4 == 0, "shape of h(+-2q)"
        k3 = [r["k3"] for r in rs]
        line = (f"real q={q}: {len(rs)} polynomials, {sum(k3):.0f} surfaces, per polynomial mean {sum(k3) / len(rs):.2f} "
                f"max {max(k3):.0f}, mean/q^2 {sum(k3) / len(rs) / q ** 2:.4f}")
        recs = []
        for r in rs:
            if r["index"] != 1:
                continue
            rec = analyse(q, prime_of(q), r["e"], cap)
            assert rec and rec["maximal"] and rec["N"] == r["N"], "disagreement with sixlines_cm.jsonl"
            rec["surfaces"] = r["k3"]
            rec["source"] = "real"
            recs.append(rec)
            if fh:
                fh.write(json.dumps(rec) + "\n")
        if recs:
            line += "; maximal order: " + ", ".join(f"|G|={r['G']} N={r['N'][0]}" for r in recs)
        print(line, flush=True)


def candidates_all(q, p):
    top = math.isqrt(16 * q) + 1
    for e1 in range(1, 6 * q):
        if e1 % p == 0:
            continue
        for m in range(1, top + 1):
            for n in range(1, top + 1):
                yield poly_from(q, e1, m, n)


def run_all(q, fh, cap):
    p = prime_of(q)
    recs = []
    for e in candidates_all(q, p):
        if admissible(q, p, e) is None:
            continue
        rec = analyse(q, p, e, cap)
        recs.append(rec)
        if rec and fh:
            rec["source"] = "all"
            fh.write(json.dumps(rec) + "\n")
    real = {tuple(json.loads(line)["e"]) for line in open("sixlines_cm.jsonl") if json.loads(line)["p"] == q}
    note = ""
    if real:
        assert real == {tuple(r["e"]) for r in recs}, "the admissible polynomials are not those of the enumerated surfaces"
        note = "  [same polynomials as the enumerated surfaces]"
    shapes = Counter(tuple(r["mu_nu"]) for r in recs)
    note += f"  (mu, nu): {dict(sorted(shapes.items()))}"
    print("all    " + summary(q, recs) + note, flush=True)


def run_sample(q, count, fh, cap, rng):
    """Random polynomials with mu = nu = 1 (m = n = 4), the only shape with a maximal order met in 'all'."""
    p = prime_of(q)
    recs, got = [], 0
    order = list(range(1, 6 * q))
    rng.shuffle(order)
    for e1 in order:
        if got >= count:
            break
        e = poly_from(q, e1, 4, 4)
        if e1 % p == 0 or admissible(q, p, e) is None:
            continue
        rec = analyse(q, p, e, cap)
        recs.append(rec)
        if rec and rec.get("maximal") and not rec.get("failed"):
            got += 1
            if fh:
                rec["source"] = "sample"
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
    print("sample " + summary(q, recs), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--all", type=int, nargs="*", default=[])
    ap.add_argument("--sample", nargs="*", default=[], help="q:count")
    ap.add_argument("--fourfold", type=int, nargs="*", help="check the fourfold identity on the real polynomials and on all of these q")
    ap.add_argument("--cap", type=int, default=3000, help="largest class number of E for which lattices are counted")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    fh = open(args.out, "w") if args.out else None
    rng = random.Random(args.seed)
    if args.fourfold is not None:
        run_fourfold(args.fourfold)
    if args.real:
        run_real(fh, args.cap)
    for q in args.all:
        run_all(q, fh, args.cap)
    for item in args.sample:
        q, count = (int(x) for x in item.split(":"))
        run_sample(q, count, fh, args.cap, rng)


if __name__ == "__main__":
    main()
