#!/usr/bin/env python
"""
sixlines_hermitian.py -- multiplicity 6, l^4 exactly dividing P(+-q), -P(q)P(-q) not a square (note 12).

Setting of note 11: an odd prime l prime to q, lam = +-q a root of P mod l of multiplicity 6, exactly
l^4 dividing P(lam), g = F / lam.  Here -P(lam)P(-lam) is not a square in Q_l.  Write k for the
unramified quadratic extension of Q_l, w for a square root of

      u = -P(lam) / l^4                    (a unit of Z_l that is not a square),

A2 = trace(g), rho the square root of P(-lam) / q^6 that is 8 mod l, and

      f0 = X^4 + [ 2 u l^2/q^3  X^3 + (12 - 2 A2) u X^2 - 2 u^2 l^2/q^3  X + (2 + A2 - rho) u^2 ] / (2 + A2 + rho).

Note 12 states, with B0 = Q_l[X]/f0, S0 = Z_l[X]/f0, theta the class of X, B = B0 (x) k, S = S0 (x) O_k,
bar the conjugation of k on B, and h(v, v') = Tr_{B/k}(beta v bar(v')) for beta in B0:

  (1) f0 = X^4 mod l, and xi = (1 + theta/w)/(1 - theta/w) has  xi bar(xi) = 1,  norm 1, and the
      eigenvalues of g are the products of two eigenvalues of xi;
  (2) for every form with a unimodular F-stable lattice there is a beta, determined up to norms from B
      and up to Q_l^*, such that the unimodular F-stable lattices are the S-stable O_k-lattices M of B
      with M^# = M ("type 0") or M^# = M / l ("type 1"), # the dual for h; a rational step joins M of
      type 0 and M' of type 1 with l M < M' < M;
  (3) the Jordan blocks of theta on M / l M over the residue field of k are (4) [leaf, 1 step],
      (1,3) [middle, l + 1], (2,2) [wide, l^2 + l + 1] or (1,1,2) [wide, l + 1];
  (4) there are leaves exactly for the class of beta = 1/f0'(theta), for which S itself is a leaf;
  (5) if no factor of f0 over Q_l has even residue degree the graph is finite and, in the class of
      (4), the number of leaves of one type is  l^c prod_p (1 + 1/N p) / (1 + 1/l),  p the primes of
      B0 and l^c = [O_B0 : S0];  there are leaves of both types exactly when every prime of B0 has even
      ramification index.

The script works with the Z_l-lattices of rank 8 in B0 + w B0 stable under theta and w, and the form
Tr_{B0/Q_l}(beta (v1 v1' - u v2 v2')).  For every polynomial found it asserts (1), (3), (4), (5) on
the graphs of all classes of beta, and compares them with the graphs of unimodular F-stable lattices
computed as in note 11 from forms on O_E: equality of colour-refinement invariants (finite cases), or
of the balls around a leaf (infinite cases).

  python -I sixlines_hermitian.py [--out FILE] [--search q,l ...] [--forms K] [--both-signs]

Run from experiments/.  Default: q = 5 with l = 3, eigenvalue +q.  The enumeration of the steps is by
brute force over the planes of a space of dimension 4 over F_{l^2}; it is meant for l = 3.
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402
import sixlines_vertical as sv  # noqa: E402
import sixlines_deep as dp  # noqa: E402
import sixlines_rank6 as r6  # noqa: E402

PREC = 40

pari(r"""
hm_hnf(X) = my(d = denominator(X)); mathnf(d * X) / d;
""")
pari(r"""
hm_local(X, l) = {
  \\ the lattice equal to X at l and to Z^n elsewhere
  my(a = 0, X1, K, n = matsize(X)[1]);
  while(denominator(X * l^a) % l == 0, a++);
  X1 = X * l^a; K = valuation(matdet(X1), l) + 1;
  mathnf(concat(lift(Mod(X1, l^K)), l^K * matid(n))) / l^a;
}
""")
pari(r"""
hm_f0(P, q, sgn, l, prec) = {
  \\ [u, f0] with f0 lifted from Z_l to Z; 0 if -P(lam)P(-lam) is a square
  my(t = variable(P), G = subst(P, t, sgn * t), A2 = -polcoef(G, 5) / q, Pm = subst(G, t, -q), Pp = subst(G, t, q), u, rho, lead, c);
  if(valuation(Pp, l) != 4, error("l^4 does not exactly divide P(lam)"));
  u = -Pp / l^4;
  if(kronecker(u, l) != -1, return(0));
  rho = sqrt(Pm / q^6 + O(l^prec)); if(lift(Mod(rho - 8, l)) != 0, rho = -rho);
  if(lift(Mod(rho - 8, l)) != 0, error("the square root of P(-lam)/q^6 is not 8 mod l"));
  lead = 2 + A2 + rho;
  c = [2 * u * l^2 / q^3, (12 - 2 * A2) * u, -2 * u^2 * l^2 / q^3, (2 + A2 - rho) * u^2] / lead;
  [u, 'x^4 + lift(c[1] + O(l^prec)) * 'x^3 + lift(c[2] + O(l^prec)) * 'x^2 + lift(c[3] + O(l^prec)) * 'x + lift(c[4] + O(l^prec))];
}
""")
pari(r"""
hm_wedge(M) = {
  my(idx = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);
  matrix(6, 6, r, c, M[idx[r][1], idx[c][1]] * M[idx[r][2], idx[c][2]] - M[idx[r][1], idx[c][2]] * M[idx[r][2], idx[c][1]]);
}
""")
pari(r"""
hm_unitary(f0, u, P, q, sgn, l) = {
  \\ (1): l-adic agreement of  charpoly(Lambda^2 xi)  with the polynomial of g, and of  det(xi)  with 1.
  \\ (xi bar(xi) = 1 holds by construction: theta is fixed by the conjugation and w changes sign.)
  my(C = matcompanion(f0), w = Mod('w, 'w^2 - u), t = C / w, xi = (matid(4) + t) * (matid(4) - t)^(-1), tt = variable(P),
     d1 = charpoly(hm_wedge(xi), 'T) - subst(P, tt, sgn * q * 'T) / q^6, d3 = lift(matdet(xi) - 1), m = 10^6);
  for(i = 0, 5, my(c = lift(polcoef(d1, i, 'T)));
    for(j = 0, 1, my(z = polcoef(c, j, 'w)); if(z != 0, m = min(m, valuation(z, l)))));
  for(j = 0, 1, my(z = polcoef(d3, j, 'w)); if(z != 0, m = min(m, valuation(z, l))));
  m;
}
""")
pari(r"""
hm_setup(f0, u, bpol) = {
  \\ [Theta, W, Gram] on B0 + w B0, for beta = bpol(theta) / f0'(theta)
  my(C = matcompanion(f0), th = Mod('x, f0), be = subst(bpol, 'x, th) / subst(deriv(f0), 'x, th),
     T = matrix(4, 4, i, j, trace(be * th^(i + j - 2))), Z = matrix(4, 4));
  [matconcat([C, Z; Z, C]), matconcat([Z, u * matid(4); matid(4), Z]), matconcat([T, Z; Z, -u * T])];
}
""")
pari(r"""
hm_dual(M, G, l) = hm_local((M~ * G)^(-1), l);
""")
pari(r"""
hm_type(M, G, l) = {
  \\ a with M^# = l^-a M, or -1
  my(D = hm_dual(M, G, l));
  if(D == M, 0, if(D == M / l, 1, -1));
}
""")
pari(r"""
hm_ranks(M, Th, W, l) = {
  \\ M must be stable under Theta and W; the ranks over F_{l^2} of Theta^j on M / l M, j = 0..4
  my(Mi = M^(-1), A = Mi * Th * M, B = Mi * W * M, Tm, X = matid(8) * Mod(1, l), r = vector(5));
  if(denominator(A) != 1 || denominator(B) != 1, error("not stable"));
  Tm = Mod(A, l); r[1] = 4;
  for(j = 1, 4, X = X * Tm; r[j + 1] = matrank(X) / 2);
  r;
}
""")
pari(r"""
hm_planes(M, a, Th, W, G, l) = {
  \\ the neighbours of M, of type a: the M B + l M (divided by l if a = 1) for the planes B of M / l M over
  \\ F_{l^2} that are stable under Theta and totally isotropic for l^-a G
  my(Mi = M^(-1), Tm = Mod(Mi * Th * M, l), Wm = Mod(Mi * W * M, l), Gm = Mod(l^(-a) * M~ * G * M, l), E = matrix(8, 0), out = List());
  for(i = 1, 8, my(c = matid(8)[, i] * Mod(1, l), T = concat(E, concat(Mat(c), Mat(Wm * c)))); if(matrank(T) == #E + 2, E = T));
  if(#E != 8, error("no basis over F_{l^2}"));
  forsubset([4, 2], piv,
    my(i = piv[1], j = piv[2], free1 = select(k -> k > i && k != j, [1..4]), free2 = select(k -> k > j, [1..4]));
    forvec(c1 = vector(2 * #free1, k, [0, l - 1]),
      my(u1 = E[, 2 * i - 1], B1);
      for(k = 1, #free1, u1 += c1[2 * k - 1] * E[, 2 * free1[k] - 1] + c1[2 * k] * E[, 2 * free1[k]]);
      B1 = concat(Mat(u1), Mat(Wm * u1));
      if(B1~ * Gm * B1 == 0,
        forvec(c2 = vector(2 * #free2, k, [0, l - 1]),
          my(u2 = E[, 2 * j - 1], B, X);
          for(k = 1, #free2, u2 += c2[2 * k - 1] * E[, 2 * free2[k] - 1] + c2[2 * k] * E[, 2 * free2[k]]);
          B = concat(B1, concat(Mat(u2), Mat(Wm * u2)));
          if(B~ * Gm * B == 0 && matrank(concat(B, Tm * B)) == 4,
            X = hm_hnf(concat(M * lift(B), l * M));
            listput(out, if(a, X / l, X)))))));
  Vec(out);
}
""")
pari(r"""
hm_maxsub(X, Th, W, l) = {
  \\ the maximal sublattices of X stable under Theta and W
  my(Xi = X^(-1), Tm = Mod(Xi * Th * X, l), Wm = Mod(Xi * W * X, l), K = matker(Tm~), out = List());
  forvec(c = vector(#K, k, [0, l - 1]),
    if(c != vector(#K),
      my(phi = (K * c~)~, H = matker(matconcat([Mat(phi); Mat(phi * Wm)])));
      listput(out, hm_hnf(concat(X * lift(H), l * X)))));
  Set(Vec(out));
}
""")
pari(r"""
hm_factors(f0, l, prec) = {
  \\ the factors of f0 over Q_l, lifted to Z, with [e, f, valuation of the discriminant of the maximal order]
  my(F = factorpadic(f0, l, prec)[, 1]~);
  vector(#F, i, my(g = lift(F[i]));
    if(poldegree(g) == 1, [g, 1, 1, 0],
      my(nf = nfinit([g, [l]]), pr = idealprimedec(nf, l));
      if(#pr != 1, error("a factor is not irreducible over Q_l"));
      [g, pr[1].e, pr[1].f, valuation(nf.disc, l)]));
}
""")
GP_f0, GP_unitary, GP_setup, GP_dual, GP_type, GP_ranks, GP_planes, GP_maxsub, GP_factors, GP_hnf = (
    pari(n) for n in ("hm_f0", "hm_unitary", "hm_setup", "hm_dual", "hm_type", "hm_ranks", "hm_planes", "hm_maxsub",
                      "hm_factors", "hm_hnf"))
INSIDE = pari("(X, Y) -> denominator(Y^(-1) * X) == 1")       # X inside Y

KIND = {(4,): "leaf", (1, 3): "middle", (2, 2): "wide", (1, 1, 2): "wide"}


def blocks(M, Th, W, l):
    """Sizes of the Jordan blocks of theta on M / l M over F_{l^2}."""
    r = [int(x) for x in GP_ranks(M, Th, W, l)]
    atleast = [r[j] - r[j + 1] for j in range(4)]              # number of blocks of size >= j + 1
    out = []
    for j in range(4, 0, -1):
        out += [j] * (atleast[j - 1] - (atleast[j] if j < 4 else 0))
    return tuple(sorted(out))


def explore(Th, W, G, l, seed, radius=None, cap=6000):
    """Follow the steps from a lattice of type 0.  Returns (adjacency, degree, blocks), or None."""
    def key(M, a):
        return f"{a}:{M}"
    verts = {key(seed, 0): (seed, 0)}
    dist = {key(seed, 0): 0}
    nbrs, inv = {}, {}
    queue = [key(seed, 0)]
    for k in queue:
        M, a = verts[k]
        assert int(GP_type(M, G, l)) == a, "a lattice met is not of the expected type"
        inv[k] = blocks(M, Th, W, l)
        ns = GP_planes(M, a, Th, W, G, l)
        nbrs[k] = [key(X, 1 - a) for X in ns]
        assert len(set(nbrs[k])) == len(nbrs[k])
        if radius is not None and dist[k] >= radius:
            continue
        for X in ns:
            kk = key(X, 1 - a)
            if kk not in verts:
                verts[kk] = (X, 1 - a)
                dist[kk] = dist[k] + 1
                queue.append(kk)
        if len(verts) > cap:
            return None
    adj = {k: [x for x in nbrs[k] if x in nbrs] for k in nbrs}
    return adj, {k: len(nbrs[k]) for k in nbrs}, inv


def selfdual_inside(X, Th, W, G, l, cap=3000):
    """All stable M = M^# with X^# < M < X, by descent through maximal stable sublattices."""
    layer = {str(X): X}
    found = {}
    while layer:
        nxt = {}
        for Y in layer.values():
            D = GP_dual(Y, G, l)
            if D == Y:
                found[str(Y)] = Y
            elif INSIDE(D, Y):
                for Z in GP_maxsub(Y, Th, W, l):
                    nxt.setdefault(str(Z), Z)
        layer = nxt
        if len(layer) > cap:
            return None
    return found


def classes(fac, f0, l):
    """One beta = pi(theta) / f0'(theta) for each class: {parities of v(pi) at the primes of odd f: pi}.
    The classes are taken up to the class of l; the sum of the parities is even."""
    odd = [i for i, (g, e, f, d) in enumerate(fac) if f % 2]
    ell = tuple(fac[i][1] % 2 for i in odd)
    seen = {}
    x = pari("'x")
    for b in itertools.product([0, 1, -1, l, -l, 2, l * l], repeat=4):
        if not any(b):
            continue
        pol = sum(c * x ** i for i, c in enumerate(b))
        tot = pari.polresultant(f0, pol)
        if tot == 0:
            continue
        vals = []
        for i in odd:
            v = int(pari.valuation(pari.polresultant(fac[i][0], pol), l))
            assert v % fac[i][2] == 0
            vals.append(v // fac[i][2] % 2)
        if sum(vals) % 2:
            continue
        vals = tuple(vals)
        cls = min(vals, tuple((a + b) % 2 for a, b in zip(vals, ell)))
        tot = int(pari.valuation(tot, l))
        if cls not in seen or tot < seen[cls][0]:
            seen[cls] = (tot, pol)
    want = 2 ** max(len(odd) - 1, 0) // (2 if any(ell) and len(odd) > 1 else 1)
    assert len(seen) == want, "a class of forms has no representative among the elements tried"
    return {c: pol for c, (tot, pol) in seen.items()}, ell


def leaf_count(fac, f0, l):
    """(5): the number of leaves of one type in the finite case."""
    disc = int(pari.valuation(pari.poldisc(f0), l))
    c2 = disc - sum(d for g, e, f, d in fac)
    assert c2 % 2 == 0 and c2 >= 0
    n = Fraction(l) ** (c2 // 2) / (1 + Fraction(1, l))
    for g, e, f, d in fac:
        n *= 1 + Fraction(1, l ** f)
    assert n.denominator == 1
    return int(n), c2 // 2


def hermitian_side(P, q, sign, l, args):
    """The graphs of all classes of beta.  Returns None if -P(lam)P(-lam) is a square."""
    res = GP_f0(P, q, sign, l, PREC)
    if res == 0:
        return None
    u, f0 = res
    assert pari("(f, l) -> lift(Mod(f, l)) == 'x^4")(f0, l), "(1) f0 is not X^4 mod l"
    agree = int(GP_unitary(f0, u, P, q, sign, l))
    assert agree >= PREC - 6, "(1) xi is not unitary of norm 1 with the eigenvalues of g"
    fac = [(g, int(e), int(f), int(d)) for g, e, f, d in GP_factors(f0, l, PREC)]
    finite = all(f % 2 for g, e, f, d in fac)
    cls, ell = classes(fac, f0, l)
    out = {"u": int(u), "factors": sorted([int(pari.poldegree(g)), e, f] for g, e, f, d in fac), "finite": finite,
           "classes": [], "graphs": {}, "ball": None}
    for c, pol in sorted(cls.items()):
        Th, W, G = GP_setup(f0, u, pol)
        I = pari.matid(8)
        leafclass = not any(c)
        if leafclass:
            assert int(GP_type(I, G, l)) == 0, "(4) S is not self-dual for beta = 1/f0'"
            seed = I
        else:
            box = selfdual_inside(GP_hnf(pari.concat(I, GP_dual(I, G, l))), Th, W, G, l)
            seed = next(iter(box.values())) if box else None
        rec = {"class": list(c), "pi": str(pol), "leaf class": leafclass}
        out["classes"].append(rec)
        if seed is None:
            rec["lattices"] = "no starting lattice found"
            continue
        res = explore(Th, W, G, l, seed, radius=None if finite else args.radius, cap=args.cap)
        if res is None:
            rec["lattices"] = "too many"
            continue
        adj, degree, inv = res
        for k in adj:                                          # (3)
            want = {(4,): 1, (1, 3): l + 1, (2, 2): l * l + l + 1, (1, 1, 2): l + 1}
            assert inv[k] in want, "(3) a fifth type of blocks"
            assert degree[k] == want[inv[k]], "(3) number of steps"
        leaves = [k for k in adj if inv[k] == (4,)]
        assert bool(leaves) == leafclass, "(4) leaves in the wrong class"
        shape = Counter((KIND[inv[k]], inv[k], degree[k]) for k in adj)
        rec.update({"lattices": len(adj), "closed": finite, "types": [sum(k[0] == "0" for k in adj), sum(k[0] == "1" for k in adj)],
                    "shape": [[a, list(b), d, n] for (a, b, d), n in sorted(shape.items())]})
        if finite:
            if leafclass:                                      # (5)
                n, cc = leaf_count(fac, f0, l)
                l0 = sum(k[0] == "0" for k in leaves)
                l1 = len(leaves) - l0
                assert l0 == n and l1 == (0 if any(ell) else n), "(5) number of leaves"
                rec["leaves"] = [l0, l1]
                rec["index"] = cc
            out["graphs"][(len(adj), dp.whole_hash(adj))] = leafclass
        elif leafclass:
            out["ball"] = dp.ball_hash(adj, degree, "0:" + str(I), args.radius)
    return out


def check_polynomial(q, e, l, sign, args):
    P = cm.pol(cm.charpoly(q, e))
    herm = hermitian_side(P, q, sign, l, args)
    if herm is None:
        return None
    case = dp.Case(q, e, l, sign)
    pat = r6.pattern(case.info)
    finite = all(case.self)
    assert finite == herm["finite"], "finite on one side and not on the other"
    graphs, ball = herm.pop("graphs"), herm.pop("ball")
    direct = []
    seen = set()
    for G in case.forms(args.forms):
        seed = dp.lagrangian_seed(case, G)
        if seed is None:
            continue
        res = case.explore(G, seed, radius=None if finite else args.radius + 1, cap=args.cap)
        if res is None:
            direct.append({"lattices": "too many"})
            continue
        adj, degree, inv = res
        leaves = [k for k in adj if inv[k][0] == (1, 5)]
        rec = {"lattices": len(adj), "leaves": len(leaves)}
        if finite:
            got = (len(adj), dp.whole_hash(adj))
            assert got in graphs, "(2) the graph of lattices is not the graph of a class of beta"
            assert graphs[got] == bool(leaves), "(4) the graph with leaves is not that of the leaf class"
            seen.add(got)
            rec["compared"] = "whole"
        else:
            dist = dp.ball(adj, str(seed), args.radius)
            near = sorted((dist[k], k) for k in dist if inv[k][0] == (1, 5))
            if near and ball is not None:
                leaf = near[0][1]
                sub = case.explore(G, pari(leaf), radius=args.radius, cap=args.cap)
                if sub is not None:
                    assert dp.ball_hash(sub[0], sub[1], leaf, args.radius) == ball, "(2) the ball around a leaf differs"
                    rec["compared"] = "ball"
        direct.append(rec)
    herm.update({"q": q, "e": e, "l": l, "eigenvalue": "+q" if sign == 1 else "-q", "primes": case.info, "pattern": pat,
                 "direct": direct, "all classes met": finite and seen == set(graphs)})
    return herm


def search(q, l, args, fh, total):
    for sign in ((1, -1) if args.both_signs else (1,)):
        for e, c in dp.candidates(q, l, sign, 4):
            if c != 4:
                continue
            P = cm.pol(cm.charpoly(q, e))
            if not pari.polisirreducible(P):
                continue
            if dp.predict(P, q, sign, l)[0] != 6:
                continue
            nf0 = pari.nfinit(P)
            dec = [[int(x) for x in row] for row in sv.GP_dec(nf0, q, l)]
            if any(p[3] == 2 for p in dec if p[4 if sign == 1 else 5] == 0):
                continue
            try:
                rec = check_polynomial(q, e, l, sign, args)
            except AssertionError:
                print(f"FAILED: q={q} e={e} l={l} sign={sign}", flush=True)
                raise
            if rec is None:
                total["square (not treated here)"] += 1
                continue
            total["polynomials"] += 1
            total[rec["pattern"]] += 1
            total["classes of beta"] += len(rec["classes"])
            total["classes with lattices"] += sum(isinstance(r["lattices"], int) for r in rec["classes"])
            total["hermitian lattices"] += sum(r["lattices"] for r in rec["classes"] if isinstance(r["lattices"], int))
            total["forms compared"] += sum("compared" in r for r in rec["direct"])
            total["polynomials compared"] += any("compared" in r for r in rec["direct"])
            total["polynomials with every class met"] += bool(rec["all classes met"])
            if rec["pattern"] != "(3/2, 1/4, 1/4)" or len(rec["classes"]) != 1:
                print(f"q={q} e={e} l={l} {rec['eigenvalue']} primes {rec['primes']} pattern {rec['pattern']} f0 {rec['factors']}: "
                      + "; ".join(f"class {r['class']}: {r['lattices']}{'' if r.get('closed', True) else ' (ball)'}"
                                  + (f" types {r['types']}" if 'types' in r else "")
                                  + (f" leaves {r['leaves']}" if 'leaves' in r else "") for r in rec["classes"])
                      + " | direct: " + ", ".join(f"{r['lattices']}{' ' + r['compared'] if 'compared' in r else ''}" for r in rec["direct"]),
                      flush=True)
            if fh:
                fh.write(json.dumps(rec) + "\n")
    print(f"search q={q} l={l} done", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    ap.add_argument("--search", action="append", help="q,l (repeatable); default 5,3")
    ap.add_argument("--forms", type=int, default=2, help="forms tried per class of valuations on the side of O_E")
    ap.add_argument("--radius", type=int, default=3, help="radius of the balls when the graph is infinite")
    ap.add_argument("--cap", type=int, default=6000)
    ap.add_argument("--both-signs", action="store_true")
    args = ap.parse_args()
    fh = open(args.out, "w") if args.out else None
    total = Counter()
    for item in args.search or ["5,3"]:
        q, l = (int(x) for x in item.split(","))
        search(q, l, args, fh, total)
    print("summary:", {str(k): v for k, v in sorted(total.items(), key=str)})


if __name__ == "__main__":
    main()
