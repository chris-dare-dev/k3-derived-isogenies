#!/usr/bin/env python
"""
sixlines_levels.py -- the structure below the top level at eigenvalue +-q (note 8), checked on
every lattice.

Setting: an odd prime l prime to q, a polynomial P with +-q a root of P mod l of multiplicity >= 4,
exactly l^2 dividing P(+-q), and a lattice class that is maximal at l.  Note 8 proves, for a
unimodular F-stable lattice L, with lam = +-q the eigenvalue and g = F - lam:

  (1) level:    [O_E L : L] = l^(2i) for an integer i >= 0, the level; i = 0 exactly for maximal L;
  (2) ring:     [O_E : End(L)] = l^(2i) at l;
  (3) depth:    the levels are 0, ..., n with 2n + 1 the exponent of l in [O_E : Z[F, V]];
  (4) lines:    L has l + 1 isotropic eigenlines if i < n and exactly one if i = n;
                equivalently l + 1 exactly when l L is contained in g^2 L at l;
  (5) parent:   for i >= 1 exactly one step goes to level i - 1, and it is the lattice
                g L + g^(2i-1) O_E L (at l); every other step goes to level i + 1;
  (6) top:      from level 0 the steps go to level 1, except two translations to level 0 when the
                primes above l are a pair l != lbar.

The script takes the pairs (polynomial, l) of sixlines_vertical.jsonl that have a step going down at
an eigenvalue +-q, recomputes the lattice classes and the neighbours, and tests (1) to (6) for every
class.  Levels of neighbours are computed directly, without identifying their classes.

  python -I sixlines_levels.py [--out FILE]

Run from experiments/, after sixlines_vertical.py.
"""
import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402
import sixlines_vertical as sv  # noqa: E402

pari(r"""
lv_span(nf, L) = {
  \\ O_E L
  my(n = poldegree(nf.pol), w = nf.zk);
  sl_hnfq(matconcat(vector(n, j, sl_mulmat(nf, w[j]) * L)));
}
""")
pari(r"""
lv_level2(nf, L, l) = valuation(abs(matdet(L) / matdet(lv_span(nf, L))), l);
""")
pari(r"""
lv_at(X, L, l, K) = sl_hnfq(concat(sl_latcap(X, L / l^K), l^K * L));
""")
pari(r"""
lv_parent(nf, L, g, i, l, K) = {
  \\ the lattice equal to g L + g^(2i-1) O_E L at l and to L at the other primes
  my(X = sl_hnfq(concat(sl_mulmat(nf, g) * L, sl_mulmat(nf, g^(2 * i - 1)) * lv_span(nf, L))));
  lv_at(X, L, l, K);
}
""")
pari(r"""
lv_crit(nf, L, g, l, K) = {
  \\ is l L contained in g^2 L at l ?
  my(X = lv_at(sl_mulmat(nf, g^2) * L, L, l, K));
  denominator(X^(-1) * (l * L)) == 1;
}
""")
GP_level2, GP_parent, GP_crit = pari("lv_level2"), pari("lv_parent"), pari("lv_crit")


def check_pair(q, e, split, l, index_val):
    P = cm.pol(cm.charpoly(q, e))
    bnf = pari.bnfinit(P, 1)
    nf = cm.GP_nf(bnf)
    gens = cm.GP_sl_gens(nf, q) if split else pari("sl_gens2")(nf, q)
    classes = cm.lattice_classes(bnf, q, gens, int(pari("(b)->b.tu[1]")(bnf)))
    F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(nf)
    lam = next(s * q % l for s in (1, -1) if int(pari.subst(P, "t", s * q)) % l == 0)
    sign = 1 if (lam - q) % l == 0 else -1
    g = F - sign * q
    dec = [[int(x) for x in row] for row in sv.GP_dec(nf, q, l)]
    pair = any(not r[2] and r[0] >= 2 and r[4 if sign == 1 else 5] > 0 for r in dec)
    assert int(pari.valuation(pari.subst(P, "t", sign * q), l)) == 2, "l^2 does not divide P(+-q) exactly"
    K = index_val + 4
    stats = Counter()
    levels = []
    for L, al in classes:
        v = int(GP_level2(nf, L, l))
        assert v % 2 == 0, "(1) the index of L in O_E L is not a square"
        levels.append(v // 2)
    n = max(levels)
    assert index_val == 2 * n + 1, "(3) depth"
    assert set(levels) == set(range(n + 1)), "(3) a level is empty"
    for (L, al), i in zip(classes, levels):
        O = sv.GP_order(nf, L)
        assert int(pari.valuation(pari.matdet(O), l)) == 2 * i, "(2) ring"
        G = sv.GP_gram(nf, q, al)
        nb = [st[3] for st in sv.GP_steps(nf, q, L, G, O, l) if int(st[0]) == lam and int(st[2])]
        crit = bool(GP_crit(nf, L, g, l, K))
        assert crit == (i < n), "(4) criterion against level"
        assert len(nb) == (l + 1 if crit else 1), "(4) number of isotropic eigenlines"
        lv = [int(GP_level2(nf, M, l)) // 2 for M in nb]
        if i == 0:
            want = Counter({0: 2, 1: l - 1}) if pair else Counter({1: l + 1})
            assert Counter(lv) == want, "(6) steps from the top"
        else:
            par = GP_parent(nf, L, g, i, l, K)
            hits = [k for k, M in enumerate(nb) if M == par]
            assert len(hits) == 1, "(5) the parent is not among the neighbours"
            assert lv[hits[0]] == i - 1, "(5) level of the parent"
            assert all(x == i + 1 for k, x in enumerate(lv) if k != hits[0]), "(5) the other steps"
        stats[i] += 1
        stats["neighbours"] += len(nb)
    return {"q": q, "e": e, "split": split, "l": l, "eigenvalue": "+q" if sign == 1 else "-q",
            "primes": "pair" if pair else "inert", "n": n,
            "classes_per_level": [stats[i] for i in range(n + 1)], "neighbours": stats["neighbours"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    args = ap.parse_args()
    import sixlines_nonsplit as ns
    ns.setup_lattice()
    fh = open(args.out, "w") if args.out else None
    total = Counter()
    for r in map(json.loads, open("sixlines_vertical.jsonl")):
        if not any(s[0] == "down" and (x["d+"] == 1 or x["d-"] == 1) for x in r["rows"] for s in x["steps"]):
            continue
        rec = check_pair(r["q"], r["e"], r["split"], r["l"], r["index_val"])
        print(f"q={rec['q']} e={rec['e']} {'split' if rec['split'] else 'non-split'} l={rec['l']} "
              f"eigenvalue {rec['eigenvalue']} primes: {rec['primes']}; depth n={rec['n']}; "
              f"classes per level {rec['classes_per_level']}; {rec['neighbours']} neighbours; (1)-(6) hold",
              flush=True)
        total["pairs"] += 1
        total["classes"] += sum(rec["classes_per_level"])
        total["neighbours"] += rec["neighbours"]
        if fh:
            fh.write(json.dumps(rec) + "\n")
    print("summary:", dict(total))


if __name__ == "__main__":
    main()
