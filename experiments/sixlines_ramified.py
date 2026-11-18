#!/usr/bin/env python
"""
sixlines_ramified.py -- eigenvalue +-q at primes ramified over E_0 (note 9), checked on every lattice.

Setting: an odd prime l prime to q and a polynomial P, with exactly l^2 dividing P(lam), lam = +-q, such
that the primes of E above l at which F = lam are ramified over E_0: one prime of degree two
(configuration (c)) or two primes of degree one (configuration (d)).  Write LL for their product and
g = F - lam.  Note 9 proves, for the unimodular F-stable lattices L:

  (1) base:    there is one O_E-lattice N with dual LL^-1 N; every L has LL^j N < L < LL^(-1-j) N for a
               type j >= 0, and [O_E L : L] = l^(2j+1);
  (2) ring:    End(L) = Z_l[F] + LL^(2j+1) O_E at l, of index l^(2j+1) in O_E;
  (3) depth:   the types are 0, ..., n-1 with l^(2n) the index of Z_l[F] in O_E at these primes;
               there are 2 l^j lattices of type j, or no lattice at all;
  (4) lines:   L has l + 1 isotropic eigenlines if j < n-1 and exactly one if j = n-1; equivalently
               l + 1 exactly when l L is contained in g^2 L at l;
  (5) top:     the two lattices of type 0 are joined by a step; it is a translation in (d) and not in (c);
  (6) parent:  for j >= 1 exactly one step goes to type j-1, to g L + LL^(j-1) N; all other steps
               from any lattice go to type j+1;
  (7) stable:  the eigenline of the step of (5), (6) is the only one stable under End(L).

Part A takes the pairs (polynomial, l) of sixlines_vertical.jsonl in configurations (c), (d); all have
n = 1.  For every lattice class it tests (1), (2), (4), (5), (7).

Part B searches the Weil polynomials t^3 h(t + q^2/t), h a cubic with three real roots in (-2q, 2q),
for a given q and l, keeps the irreducible ones in configuration (c) or (d) with no other prime above l
ramified over E_0, and for each builds a form on O_E for which N = O_E, lists the lattices twice,

  - by following the steps from one lattice of type 0,
  - by listing every F-stable lattice between LL^n N and LL^(-1-n) N and keeping the unimodular ones,

compares the two lists, and tests (1) to (7) on every lattice.  These polynomials are not claimed to
come from surfaces of the family: the statements are about lattices.

  python -I sixlines_ramified.py [--out FILE] [--search q,l ...] [--every K]

Run from experiments/, after sixlines_vertical.py.  The default search is q = 5 with l = 3 and l = 7.
--every K tests only one polynomial in K among those with n = 1; all those with n >= 2 are tested.
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
import sixlines_levels as lv  # noqa: E402,F401  (defines lv_span, lv_level2, lv_at, lv_crit)

pari(r"""
rm_steps(nf, L, G, Rord, l, lam) = {
  \\ the lines of the eigenspace of F for lam on L / l L, for a lattice L that is F-stable and
  \\ unimodular at l.  Returns [mu, [[isotropic, neighbour, stable under Rord], ...]].
  my(P = nf.pol, n = poldegree(P), Fm = sl_mulmat(nf, Mod(variable(P), P)), Li = L^(-1),
     FL = Li * Fm * L, GL = L~ * G * L, D, W, mu, out = List(), om, f, x, gx, c, k, t, xt, K, M, st);
  if(denominator(FL) != 1, error("not F-stable"));
  D = denominator(GL);
  if(D % l == 0, error("not integral at l"));
  GL = D * GL;
  if(matdet(GL) % l == 0, error("not unimodular at l"));
  om = vector(n, j, Li * sl_mulmat(nf, nfbasistoalg(nf, Rord[,j])) * L);
  W = lift(matker(Mod(FL - lam * matid(n), l))); mu = #W;
  forvec(e = vector(mu, i, [0, l - 1]),
    f = 0; for(i = 1, mu, if(e[i] != 0, f = e[i]; break));
    if(f != 1, next);
    x = W * e~;
    st = 1;
    for(j = 1, n, if(matrank(Mod(concat(Mat(x), om[j] * x), l)) > 1, st = 0; break));
    gx = GL * x; c = x~ * gx;
    if(c % l, listput(out, [0, 0, st]); next);
    k = 0; for(j = 1, n, if(gx[j] % l, k = j; break));
    t = lift(Mod(-c / l, l) / Mod(2 * gx[k], l));
    xt = x; xt[k] += l * t;
    if((xt~ * GL * xt) % l^2, error("lift is not isotropic mod l^2"));
    K = matrix(n, n, a, b, if(b == k, if(a == k, l, 0),
                              if(a == b, 1, if(a == k, lift(-Mod(gx[b], l) / Mod(gx[k], l)), 0))));
    M = sl_hnfq(concat(K, xt / l));
    if(abs(matdet(M)) != 1, error("neighbour has the wrong index"));
    if(denominator(M~ * GL * M) != 1, error("neighbour is not integral"));
    if(denominator(M^(-1) * FL * M) != 1, error("neighbour is not F-stable"));
    listput(out, [1, sl_hnfq(L * M), st]));
  [mu, Vec(out)];
}
""")
pari(r"""
rm_maxsub(nf, X, g, l) = {
  \\ the sublattices of index l of X that contain l X + g X: the maximal F-stable sublattices of X at
  \\ the primes above l where g is not a unit
  my(n = #X, Y = sl_hnfq(concat(l * X, sl_mulmat(nf, g) * X)), C = X^(-1) * Y, Ph, d, out = List(), f, phi, k, K);
  Ph = lift(matker(Mod(C~, l))); d = #Ph;
  forvec(e = vector(d, i, [0, l - 1]),
    f = 0; for(i = 1, d, if(e[i] != 0, f = e[i]; break));
    if(f != 1, next);
    phi = Ph * e~;
    k = 0; for(j = 1, n, if(phi[j] % l, k = j; break));
    K = matrix(n, n, a, b, if(b == k, if(a == k, l, 0),
                              if(a == b, 1, if(a == k, lift(-Mod(phi[b], l) / Mod(phi[k], l)), 0))));
    listput(out, sl_hnfq(X * K)));
  Vec(out);
}
""")
pari(r"""
rm_unimodular(X, G, l) = {
  my(GL = X~ * G * X, D = denominator(GL));
  D % l != 0 && matdet(D * GL) % l != 0;
}
""")
pari(r"""
rm_setup(nf, q, l, sgn) = {
  \\ [LL, al0, units]: the product of the primes above l at which F = sgn q; an element al0 of E_0 that
  \\ generates the inverse different of E_0 at l; and some elements of E_0 that are units at l
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = q^2 / F, s = F + V, h = minpoly(s, 'x),
     nf0 = nfinit(h), dec = idealprimedec(nf, l), dec0 = idealprimedec(nf0, l), LL = idealhnf(nf, 1), a0, un = List());
  for(i = 1, #dec, if(idealval(nf, F - sgn * q, dec[i]) > 0, LL = idealmul(nf, LL, dec[i])));
  a0 = idealappr(nf0, matrix(#dec0, 2, i, j, if(j == 1, dec0[i], -idealval(nf0, nf0.diff, dec0[i]))));
  forvec(c = vector(3, i, [0, l - 1]),
    my(u = nfbasistoalg(nf0, c~));
    if(u != 0 && gcd(numerator(norm(u)), l) == 1, listput(un, subst(lift(u), 'x, s))));
  [LL, subst(lift(nfbasistoalg(nf0, a0)), 'x, s), Vec(un)];
}
""")
pari(r"""
rm_ring(nf, q, LL, k, l) = {
  \\ the order equal to Z[F, V] + LL^k O_E at l and to O_E elsewhere
  my(P = nf.pol, F = Mod(variable(P), P), B = sl_order(nf, [F, q^2 / F]));
  sv_local(sl_hnfq(concat(B, idealpow(nf, LL, k))), l);
}
""")
GP_steps, GP_maxsub, GP_unimodular, GP_setup, GP_ring = (
    pari(n) for n in ("rm_steps", "rm_maxsub", "rm_unimodular", "rm_setup", "rm_ring"))
GP_level2, GP_crit, GP_at, GP_span = pari("lv_level2"), pari("lv_crit"), pari("lv_at"), pari("lv_span")
HNF = pari("sl_hnfq")
MUL = pari("sl_mulmat")
INSIDE = pari("(A, B) -> denominator(A^(-1) * B) == 1")          # is B contained in A ?
CONCAT = pari("(A, B) -> sl_hnfq(concat(A, B))")
IDPOW = pari("(nf, I, k) -> sl_hnfq(idealhnf(nf, idealpow(nf, I, k)))")


def val(n, l):
    v = 0
    while n % l == 0:
        n //= l
        v += 1
    return v


def configuration(dec, sign):
    """'c', 'd' or None from the rows [e, f, self-conjugate, e over E_0, v(F - q), v(F + q)]."""
    col = 4 if sign == 1 else 5
    pr = [p for p in dec if p[col] > 0]
    if any(p[3] == 2 for p in dec if p[col] == 0):
        return None                                   # another prime above l is ramified over E_0
    if len(pr) == 1 and pr[0][1] == 2 and pr[0][2] and pr[0][3] == 2 and pr[0][col] == 1:
        return "c"
    if len(pr) == 2 and all(p[1] == 1 and p[2] and p[3] == 2 and p[col] == 1 for p in pr):
        return "d"
    return None


def local_depth(nf, q, LL, l):
    """n with l^(2n) the index of Z_l[F] in O_E at the primes of LL."""
    v = int(pari.valuation(pari.matdet(GP_ring(nf, q, LL, 40, l)), l))
    assert v % 2 == 0 and v < 40, "the index at these primes is not an even power"
    return v // 2


# ---------------------------------------------------------------------------------------------------
# Part A: the lattice classes of the data
# ---------------------------------------------------------------------------------------------------

def check_data_pair(q, e, split, l):
    P = cm.pol(cm.charpoly(q, e))
    bnf = pari.bnfinit(P, 1)
    nf = cm.GP_nf(bnf)
    F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(nf)
    dec = [[int(x) for x in row] for row in sv.GP_dec(nf, q, l)]
    found = [(s, configuration(dec, s)) for s in (1, -1)
             if int(pari.subst(P, "t", s * q)) % l == 0 and configuration(dec, s)]
    assert len(found) == 1
    sign, cfg = found[0]
    lam = sign * q % l
    g = F - sign * q
    assert val(int(pari.subst(P, "t", sign * q)), l) == 2
    LL = GP_setup(nf, q, l, sign)[0]
    n = local_depth(nf, q, LL, l)
    assert n == 1
    gens = cm.GP_sl_gens(nf, q) if split else pari("sl_gens2")(nf, q)
    classes = cm.lattice_classes(bnf, q, gens, int(pari("(b)->b.tu[1]")(bnf)))
    K = 6
    ring = GP_ring(nf, q, LL, 1, l)
    for L, al in classes:
        assert int(GP_level2(nf, L, l)) == 1, "(1) index of L in O_E L"
        O = sv.GP_order(nf, L)
        assert sv.GP_local(O, l) == ring, "(2) ring"
        G = sv.GP_gram(nf, q, al)
        lines = [st for st in sv.GP_steps(nf, q, L, G, O, l) if int(st[0]) == lam]
        assert len(lines) == l + 1 and all(int(st[1]) == 2 for st in lines)
        iso = [st for st in lines if int(st[2])]
        assert len(iso) == 1 and not GP_crit(nf, L, g, l, K), "(4) one isotropic eigenline"
        assert [int(st[4]) for st in lines if int(st[2])] == [1], "(7) the isotropic line is stable"
        assert not any(int(st[4]) for st in lines if not int(st[2])), "(7) no other line is stable"
        Lp = iso[0][3]
        span = GP_span(nf, L)
        assert int(GP_level2(nf, Lp, l)) == 1
        assert CONCAT(L, Lp) == GP_at(span, L, l, K), "(5) L + L' = O_E L"
        assert pari("sl_latcap")(L, Lp) == GP_at(MUL(nf, g) * span, L, l, K), "(5) L n L' = LL O_E L"
        tr = any(Lp == M for M in sv.GP_translates(nf, q, L, l, K))
        assert tr == (cfg == "d"), "(5) translation"
    return {"part": "A", "q": q, "e": e, "split": split, "l": l, "eigenvalue": "+q" if sign == 1 else "-q",
            "configuration": cfg, "n": n, "classes": len(classes)}


# ---------------------------------------------------------------------------------------------------
# Part B: polynomials found by a search, with all lattices at l listed
# ---------------------------------------------------------------------------------------------------

def candidates(q, l, sign):
    """Cubics h = x^3 - a x^2 + b x - c with three real roots in (-2q, 2q), a double root 2 sign q
    mod l, and exactly l^2 dividing h(2 sign q)."""
    x0, B = 2 * sign * q, 2 * q
    for a in range(-3 * B + 1, 3 * B):
        for b in range(-3 * B * B + 1, 3 * B * B):
            if (3 * x0 * x0 - 2 * a * x0 + b) % l or a * a - 3 * b <= 0:
                continue
            if not (3 * B * B - 2 * a * B + b > 0 and 3 * B * B + 2 * a * B + b > 0):
                continue
            base = x0 ** 3 - a * x0 ** 2 + b * x0
            for c in range(-B ** 3 + ((base + B ** 3) % (l * l)), B ** 3, l * l):
                hv = base - c
                if hv == 0 or val(hv, l) != 2:
                    continue
                if not (-B ** 3 - a * B * B - b * B - c < 0 < B ** 3 - a * B * B + b * B - c):
                    continue
                if a * a * b * b - 4 * b ** 3 - 4 * a ** 3 * c + 18 * a * b * c - 27 * c * c <= 0:
                    continue
                yield [a, b + 3 * q * q, c + 2 * a * q * q]


def tops(nf, q, G, LL, l):
    """The unimodular lattices between N = O_E and LL^-1."""
    n = 6
    N = pari.matid(n)
    M0 = IDPOW(nf, LL, -1)
    cols = [M0[i] for i in range(n)]
    vecs = [v for v in cols if not INSIDE(N, pari.Mat(v))]
    out = []
    seen = set()
    for i, u in enumerate(vecs):
        for w in vecs[i + 1:]:
            if CONCAT(CONCAT(N, pari.Mat(u)), pari.Mat(w)) != M0:
                continue
            for x in [u + c * w for c in range(l)] + [w]:
                X = CONCAT(N, pari.Mat(x))
                if str(X) not in seen and GP_unimodular(X, G, l):
                    seen.add(str(X))
                    out.append(X)
            return out
    raise AssertionError("LL^-1 / O_E is not a plane")


def all_lattices(nf, G, LL, g, l, J):
    """Every lattice between LL^J and LL^(-1-J) that is F-stable and unimodular at l."""
    low = IDPOW(nf, LL, J)
    layer = {str(X): X for X in [IDPOW(nf, LL, -1 - J)]}
    for _ in range(2 * J + 1):
        nxt = {}
        for X in layer.values():
            for Y in GP_maxsub(nf, X, g, l):
                if INSIDE(Y, low):
                    nxt.setdefault(str(Y), Y)
        layer = nxt
    return {k: X for k, X in layer.items() if GP_unimodular(X, G, l)}


def check_polynomial(q, e, l, sign, cfg, nf, n, LL, al0, units):
    F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(nf)
    lam = sign * q % l
    g = F - sign * q
    N = pari.matid(6)
    K = 2 * n + 6
    G = top = None
    for u in units:                                   # a form with N^dual = LL^-1 N and a lattice
        G = sv.GP_gram(nf, q, al0 * u)
        assert GP_at(G ** (-1), N, l, K) == IDPOW(nf, LL, -1), "the dual of O_E is not LL^-1"
        top = tops(nf, q, G, LL, l)
        assert len(top) in (0, 2)
        if top:
            break
    assert top, "no form with a unimodular lattice was found"
    rings = [GP_ring(nf, q, LL, k, l) for k in range(2 * n + 2)]
    assert rings[2 * n] == rings[2 * n + 1]
    # follow the steps from one top
    verts = {str(top[0]): top[0]}
    queue = [top[0]]
    nbrs = {}
    while queue:
        L = queue.pop()
        O = sv.GP_order(nf, L)
        mu, lines = GP_steps(nf, L, G, O, l, lam)
        assert int(mu) == 2 and len(lines) == l + 1
        nbrs[str(L)] = (O, [(bool(int(a)), M, bool(int(st))) for a, M, st in lines])
        for a, M, st in lines:
            if int(a) and str(M) not in verts:
                verts[str(M)] = M
                queue.append(M)
    # list all lattices independently
    full = all_lattices(nf, G, LL, g, l, n)
    assert set(full) == set(verts), "the two lists of lattices differ"
    typ = {}
    for k, L in verts.items():
        v = int(GP_level2(nf, L, l))
        assert v % 2 == 1, "(1) the index of L in O_E L is not an odd power"
        typ[k] = (v - 1) // 2
    count = Counter(typ.values())
    assert count == Counter({j: 2 * l ** j for j in range(n)}), "(3) number of lattices of each type"
    edges = 0
    for k, L in verts.items():
        j = typ[k]
        O, lines = nbrs[k]
        assert INSIDE(IDPOW(nf, LL, -1 - j), L) and INSIDE(L, IDPOW(nf, LL, j)), "(1) position"
        assert sv.GP_local(O, l) == rings[2 * j + 1], "(2) ring"
        assert int(pari.valuation(pari.matdet(O), l)) == 2 * j + 1, "(2) index of the ring"
        iso = [(M, st) for a, M, st in lines if a]
        crit = bool(GP_crit(nf, L, g, l, K))
        assert crit == (j < n - 1), "(4) criterion against type"
        assert len(iso) == (l + 1 if crit else 1), "(4) number of isotropic eigenlines"
        edges += len(iso)
        if j == 0:
            want = top[1] if k == str(top[0]) else top[0]
            assert k in (str(top[0]), str(top[1]))
        else:
            want = GP_at(CONCAT(MUL(nf, g) * L, IDPOW(nf, LL, j - 1)), N, l, K)
        hits = [i for i, (M, st) in enumerate(iso) if M == want]
        assert len(hits) == 1, "(5), (6) the step going up or across"
        assert all(typ[str(M)] == j + 1 for i, (M, st) in enumerate(iso) if i != hits[0]), "(6) the other steps"
        assert typ[str(want)] == max(j - 1, 0)
        assert iso[hits[0]][1], "(7) the line going up is stable"
        assert sum(st for a, M, st in lines) == 1, "(7) another line is stable"
        shifts = list(sv.GP_translates(nf, q, L, l, K))
        for i, (M, st) in enumerate(iso):
            assert any(M == T for T in shifts) == (cfg == "d" and j == 0 and i == hits[0]), "(5) translation"
    assert edges == 2 * (len(verts) - 1), "the graph is not a tree"
    return {"part": "B", "q": q, "e": e, "l": l, "eigenvalue": "+q" if sign == 1 else "-q",
            "configuration": cfg, "n": n, "lattices_per_type": [count[j] for j in range(n)], "steps": edges}


def search(q, l, every, fh, total):
    seen = Counter()
    for sign in (1, -1):
        for e in candidates(q, l, sign):
            P = cm.pol(cm.charpoly(q, e))
            if not pari.polisirreducible(P):
                continue
            nf = pari.nfinit(P)
            dec = [[int(x) for x in row] for row in sv.GP_dec(nf, q, l)]
            cfg = configuration(dec, sign)
            if not cfg:
                continue
            LL, al0, units = GP_setup(nf, q, l, sign)
            n = local_depth(nf, q, LL, l)
            seen[(cfg, n)] += 1
            if n == 1 and seen[(cfg, n)] % every != 1 and every != 1:
                continue
            rec = check_polynomial(q, e, l, sign, cfg, nf, n, LL, al0, list(units))
            total[("B", cfg, n)] += 1
            total["lattices"] += sum(rec["lattices_per_type"])
            total["steps"] += rec["steps"]
            if n >= 2:
                print(f"q={q} e={e} l={l} eigenvalue {rec['eigenvalue']} configuration ({cfg}) n={n}: "
                      f"lattices per type {rec['lattices_per_type']}; (1)-(7) hold", flush=True)
            if fh:
                fh.write(json.dumps(rec) + "\n")
    print(f"search q={q} l={l}: polynomials found by (configuration, n): "
          f"{ {f'({c}) n={n}': k for (c, n), k in sorted(seen.items())} }", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    ap.add_argument("--search", action="append", help="q,l (repeatable); default 5,3 and 5,7")
    ap.add_argument("--every", type=int, default=1)
    args = ap.parse_args()
    import sixlines_nonsplit as ns
    ns.setup_lattice()
    fh = open(args.out, "w") if args.out else None
    total = Counter()
    for r in map(json.loads, open("sixlines_vertical.jsonl")):
        P = cm.pol(cm.charpoly(r["q"], r["e"]))
        signs = [s for s in (1, -1) if int(pari.subst(P, "t", s * r["q"])) % r["l"] == 0
                 and sum(p[0] * p[1] for p in r["primes"] if p[4 if s == 1 else 5] > 0) >= 4
                 and configuration(r["primes"], s)]
        if not signs:
            continue
        rec = check_data_pair(r["q"], r["e"], r["split"], r["l"])
        print(f"data: q={rec['q']} e={rec['e']} {'split' if rec['split'] else 'non-split'} l={rec['l']} "
              f"eigenvalue {rec['eigenvalue']} configuration ({rec['configuration']}) n={rec['n']}: "
              f"{rec['classes']} classes; (1), (2), (4), (5), (7) hold", flush=True)
        total[("A", rec["configuration"])] += 1
        total["classes"] += rec["classes"]
        if fh:
            fh.write(json.dumps(rec) + "\n")
    for item in args.search or ["5,3", "5,7"]:
        q, l = (int(x) for x in item.split(","))
        search(q, l, args.every, fh, total)
    print("summary:", {str(k): v for k, v in sorted(total.items(), key=str)})


if __name__ == "__main__":
    main()
