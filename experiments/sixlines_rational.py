#!/usr/bin/env python
"""
sixlines_rational.py -- rationality over F_q of the twisted equivalences (note 6, Proposition G).

1. Artin--Tate.  For a surface of the family over F_q whose Frobenius is q g on N (g = 1 or an
   involution of type A, B, C) and has polynomial P on the rank-6 part,

       |Br(X)| = 2^(16 - rank N^g) P(q) / (q^5 |disc N^g|),

   which must be the square of an integer.  The quadratic twist has P(-q) in place of P(q).

2. Rational steps.  The step at a degree-one prime l != lbar of E is a twisted equivalence defined
   over F_q exactly when l divides F - q (F + q for the twist).  For every polynomial the script
   factors F -+ q, and for each such prime computes the permutation of the lattice classes that are
   maximal at the prime below.

  python -I sixlines_rational.py [--out FILE]

Reads the result files of notes 2, 3 and 5.  Run from experiments/.
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from math import isqrt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402

pari(r"""
sr_primes(nf, q, s) = {
  \\ primes of E dividing F - s q that are of degree one, unramified, not self-conjugate, odd, prime to q
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = q^2 / F, sig = lift(V), fa, out = List(), pr, cpr);
  fa = idealfactor(nf, F - s * q);
  for(i = 1, matsize(fa)[1],
    pr = fa[i, 1];
    if(pr.p == 2 || q % pr.p == 0 || pr.f != 1 || pr.e != 1, next);
    cpr = idealhnf(nf, nfgaloisapply(nf, sig, pr));
    if(cpr == idealhnf(nf, pr), next);
    listput(out, [pr, fa[i, 2]]));
  Vec(out);
}
""")
pari(r"""
sr_maximal_at(nf, L, l) = {
  my(n = poldegree(nf.pol), w = nf.zk, J);
  J = sl_hnfq(matconcat(vector(n, j, sl_mulmat(nf, w[j]) * L)));
  (abs(matdet(J) / matdet(L)) % l) != 0;
}
""")
pari(r"""
sr_step(nf, q, L, pr) = {
  \\ the lattice equal to (a / abar) L at l and to L elsewhere, a a uniformiser at pr and a unit at
  \\ the other primes above l; L is assumed maximal at l
  my(P = nf.pol, l = pr.p, dec = idealprimedec(nf, l), fact, a, ab, M, K = 2);
  fact = matrix(#dec, 2, i, j, if(j == 1, dec[i], if(dec[i] == pr, 1, 0)));
  a = nfbasistoalg(nf, idealappr(nf, fact));
  ab = sl_conjel(nf, q, a);
  M = sl_hnfq(sl_mulmat(nf, a / ab) * L);
  sl_hnfq(concat(sl_latcap(M, L / l^K), l^K * L));
}
""")
GP_primes, GP_max, GP_step = pari("sr_primes"), pari("sr_maximal_at"), pari("sr_step")


def invariant_lattices():
    """Rank and |disc| of N^g for g = 1 and the three involution types."""
    duads = list(itertools.combinations(range(6), 2))
    half = Fraction(1, 2)
    gens = [[Fraction(int(i == j)) for i in range(16)] for j in range(16)]
    for i in range(6):
        gens.append([-half if i in d else Fraction(0) for d in duads] + [half])
    H = pari.mathnf(pari("[" + ";".join(",".join(str(2 * g[r]) for g in gens) for r in range(16)) + "]"))
    B = [[Fraction(int(H[c][r]), 2) for c in range(16)] for r in range(16)]
    D = [-2] * 15 + [2]
    Binv = cm.mat_inv(B)
    G = [[sum(B[r][i] * D[r] * B[r][j] for r in range(16)) for j in range(16)] for i in range(16)]
    out = {}
    for name, perm in (("1", (0, 1, 2, 3, 4, 5)), ("A", (1, 0, 2, 3, 4, 5)),
                       ("B", (1, 0, 3, 2, 4, 5)), ("C", (1, 0, 3, 2, 5, 4))):
        Pm = [[Fraction(0)] * 16 for _ in range(16)]
        for d, du in enumerate(duads):
            Pm[duads.index(tuple(sorted((perm[du[0]], perm[du[1]]))))][d] = Fraction(1)
        Pm[15][15] = Fraction(1)
        A = cm.int_mat(cm.mat_mul(cm.mat_mul(Binv, Pm), B))
        K = pari.matkerint(pari.matrix(16, 16, [A[r][c] - (r == c) for r in range(16) for c in range(16)]))
        rk = len(K)
        Kc = [[int(K[c][r]) for c in range(rk)] for r in range(16)]
        Gk = [[sum(Kc[a][i] * G[a][b] * Kc[b][j] for a in range(16) for b in range(16)) for j in range(rk)]
              for i in range(rk)]
        out[name] = (rk, int(abs(cm._det(Gk))))
    return out


def brauer_order(inv, q, e, sign, typ):
    co = cm.charpoly(q, e)
    P = sum(c * (sign * q) ** (6 - i) for i, c in enumerate(co))
    rk, disc = inv[typ]
    return Fraction(2 ** (16 - rk) * P, q ** 5 * disc)


def is_square(v):
    return v.denominator == 1 and v.numerator > 0 and isqrt(v.numerator) ** 2 == v.numerator


def rational_steps(q, e, split):
    """For both twists: the primes giving rational steps and their action on the lattice classes."""
    bnf = pari.bnfinit(cm.pol(cm.charpoly(q, e)), 1)
    nf = cm.GP_nf(bnf)
    found = {s: list(GP_primes(nf, q, s)) for s in (1, -1)}
    if not any(found.values()):
        return []
    gens = cm.GP_sl_gens(nf, q) if split else pari("sl_gens2")(nf, q)
    w = int(pari("(b)->b.tu[1]")(bnf))
    classes = cm.lattice_classes(bnf, q, gens, w)
    out = []
    for s, prs in found.items():
        for pr, mult in prs:
            l = int(pari("(p)->p.p")(pr))
            top = [i for i, (L, al) in enumerate(classes) if bool(GP_max(nf, L, l))]
            perm = {}
            for i in top:
                L2 = GP_step(nf, q, classes[i][0], pr)
                hits = [j for j in top if cm.GP_sl_iso(bnf, q, L2, classes[i][1], classes[j][0], classes[j][1])]
                assert len(hits) == 1, hits
                perm[i] = hits[0]
            order = sorted(perm)
            cyc = cm.cycle_type([order.index(perm[i]) for i in order]) if order else []
            out.append({"sign": s, "l": l, "multiplicity": int(mult), "classes": len(classes),
                        "maximal_at_l": len(top), "cycles": cyc})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    args = ap.parse_args()
    inv = invariant_lattices()
    print("rank and |disc| of the invariant lattices N^g:", inv)
    import sixlines_nonsplit as ns
    ns.setup_lattice()
    cases = []
    rows = [json.loads(line) for line in open("sixlines_cm.jsonl")]
    rows += [r for r in map(json.loads, open("sixlines_neighbours.jsonl")) if r["p"] == 25]
    rows += [json.loads(line) for line in open("sixlines_q49_cm.jsonl")]
    for r in rows:
        cases.append((r["p"], r["e"], True, {"1": {1: r["a_plus"], -1: r["a_minus"]}}))
    for name in ("small", "p17", "p19", "p23"):
        for r in map(json.loads, open(f"sixlines_nonsplit_lat_{name}.jsonl")):
            cases.append((r["p"], r["e"], False, {t: {1: n, -1: n} for t, n in r["forms"].items()}))
    stats, orders, fh = Counter(), Counter(), open(args.out, "w") if args.out else None
    for q, e, split, types in cases:
        rec = {"q": q, "e": e, "split": split, "brauer": {}}
        for typ, signs in types.items():
            for s in (1, -1):
                v = brauer_order(inv, q, e, s, typ)
                ok = is_square(v)
                stats[("split" if split else "non-split", "square" if ok else "NOT a square")] += 1
                if ok:
                    orders[("split" if split else "non-split", int(v))] += 1
                rec["brauer"][f"{typ}{'+' if s == 1 else '-'}"] = str(v)
        rec["steps"] = rational_steps(q, e, split)
        if rec["steps"]:
            print(f"  q={q} e={e} {'split' if split else 'non-split'} |Br|={rec['brauer']}")
            for st in rec["steps"]:
                print(f"      twist {'+' if st['sign'] == 1 else '-'}: l={st['l']} (multiplicity {st['multiplicity']}); "
                      f"{st['maximal_at_l']} of {st['classes']} lattice classes maximal at l; cycles {st['cycles']}",
                      flush=True)
        stats[("split" if split else "non-split", "with a rational step" if rec["steps"] else "no rational step")] += 1
        if fh:
            fh.write(json.dumps(rec) + "\n")
    print("Artin-Tate and steps:", {" ".join(k): v for k, v in sorted(stats.items())})
    print("orders of Br:", {f"{k[0]} {k[1]}": v for k, v in sorted(orders.items())})


if __name__ == "__main__":
    main()
