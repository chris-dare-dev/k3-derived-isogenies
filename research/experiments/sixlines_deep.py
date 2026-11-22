#!/usr/bin/env python
"""
sixlines_deep.py -- eigenvalue +-q when l^4 divides P(+-q) (note 10).

Setting: an odd prime l prime to q, lam = +-q, V the part of the lattice at l on which F = lam, of rank m,
and l^c the exact power of l dividing P(lam).  Note 10 states:

  (O) local rule.  Let s_1, ..., s_mu be the sizes of the Jordan blocks of F - lam on its generalised
      kernel in L / l L, nu the number of blocks of size 1 and rho = mu - nu.  Then coker(F - lam) is a
      group A x A with mu generators, even-sized blocks come in pairs, and the number of rational steps is
          (l^rho - 1)/(l - 1) + l^rho * i(nu),
      with i(nu) the number of isotropic lines of a non-degenerate form of rank nu over F_l.

  (P) rank 4.  Let x^2 - S x + Pd be the factor over Q_l of the characteristic polynomial of
      (F + V)/lam with roots = 2 mod l, and t_1, t_2 the roots of z^2 - Pd z + (S^2 - 2 Pd - 4), of
      discriminant D with v(D) = c.
      - D a square: the graph is the product of two volcanoes Vol(K_i, r_i), i = 1, 2: its vertices
        are the pairs (x_1, x_2) of vertices whose parities differ by a fixed eps in {0, 1}, and a
        step moves both coordinates to a neighbour.  Here w_i = v(t_i - 2); K_i is ramified with
        r_i = (w_i - 1)/2 if w_i is odd; otherwise r_i = w_i / 2 and K_i is split or unramified
        according as (t_i - 2)(t_i + 2) / l^w_i is a square or not.  (Vertices of a volcano have a
        parity that changes along every edge.)
      - D not a square, v(D) even: one volcano over the unramified quadratic extension K of Q_l, with
        branching l^2: w = min(v(Pd - 4), v(D)/2), ramified with r = (w - 1)/2 if w is odd, otherwise
        r = w/2 and split or unramified according as (S^2 - 4 Pd)/l^(2w) is a square or not.
      - v(D) odd: there is no unimodular lattice.
      Vol(K, r) with branching Q: the tops are one vertex (unramified), a line (split) or two joined
      vertices (ramified); a vertex of level 1 <= j < r has Q neighbours of level j + 1 and one of
      level j - 1; a top has neighbours of level 1 to make Q + 1 in all; levels stop at r.

The script searches the Weil polynomials t^3 h(t + q^2/t), h a cubic with three real roots in (-2q, 2q),
with l^c exactly dividing P(lam), c >= cmin, keeps the irreducible ones with no prime above l ramified
over E_0 at the other eigenvalue, and

  rank 4:  computes the prediction (P) from l-adic factorisations; lists unimodular lattices for several
           classes of forms; follows the steps from one of them (all of them when the graph is finite,
           a ball otherwise); compares the graph with the predicted one by colour refinement; and tests
           (O) on every lattice met.  When no lattice is predicted it checks that none is found.
  rank 6:  for --rank6 K polynomials of each configuration of primes, follows the steps from one
           lattice, tests (O) on every lattice, and prints the shape of what it finds.

  python -I sixlines_deep.py [--out FILE] [--search q,l ...] [--cmin C] [--u-every K] [--rank6 K]

Run from experiments/.  Default: q = 5 with l = 3, eigenvalue +q only (P(-t) has the same structure
at -q; --both-signs runs it too); cmin 4; one polynomial in 25 of those with c = 2 and rank 4 is tested
as well (--u-every), as a comparison with notes 8 and 9.
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402
import sixlines_vertical as sv  # noqa: E402
import sixlines_ramified as rm  # noqa: E402
from sixlines_ramified import (GP_steps, GP_maxsub, GP_unimodular, GP_setup, GP_ring, GP_at,  # noqa: E402
                               IDPOW, INSIDE, CONCAT, val)

pari(r"""
dp_primes(nf, q, l, sgn) = {
  \\ the primes above l at which F = sgn q, with [e, f, self-conjugate, v(F - sgn q)]
  my(P = nf.pol, F = Mod(variable(P), P), dec = idealprimedec(nf, l), prs = List(), info = List());
  for(i = 1, #dec, my(v = idealval(nf, F - sgn * q, dec[i]));
    if(v > 0,
      listput(prs, dec[i]);
      listput(info, [dec[i].e, dec[i].f,
        idealhnf(nf, nfgaloisapply(nf, lift(q^2 / F), dec[i])) == idealhnf(nf, dec[i]), v])));
  [Vec(prs), Vec(info)];
}
""")
pari(r"""
dp_box(nf, G, prs, self, R, l, w) = {
  \\ the largest O_E-lattice that can contain a unimodular lattice stable under R (0 if there is none);
  \\ for a pair of conjugate primes the two exponents are allowed to differ by at most w
  my(n = poldegree(nf.pol), D = lv_at(G^(-1), matid(n), l, 60), cond = sl_conductor(nf, R),
     B = idealhnf(nf, 1), d, f, a);
  for(i = 1, #prs,
    d = -idealval(nf, D, prs[i]); f = idealval(nf, cond, prs[i]);
    if(self[i],
      a = ceil((-f - d) / 2); if(2 * a > -d, return(0)),
      a = ceil((-d - f - w) / 2));
    B = idealmul(nf, B, idealpow(nf, prs[i], a)));
  sl_hnfq(idealhnf(nf, B));
}
""")
pari(r"""
dp_dual(X, G, l) = lv_at((X~ * G)^(-1), matid(#X), l, 60);
""")
pari(r"""
dp_elsewhere(nf, G, prs, l) = {
  \\ is O_E unimodular for G at the primes above l that are not in prs ?
  my(n = poldegree(nf.pol), D = lv_at(G^(-1), matid(n), l, 60), dec = idealprimedec(nf, l), ok = 1, hit);
  for(i = 1, #dec,
    hit = 0; for(j = 1, #prs, if(idealhnf(nf, dec[i]) == idealhnf(nf, prs[j]), hit = 1));
    if(!hit && idealval(nf, D, dec[i]) != 0, ok = 0));
  ok;
}
""")
pari(r"""
dp_jordan(nf, L, g, l) = {
  \\ sizes of the Jordan blocks of g on its generalised kernel in L / l L
  my(n = #L, M = Mod(L^(-1) * sl_mulmat(nf, g) * L, l), r = vector(n + 2, k, n - matrank(M^(k-1))), out = List());
  for(k = 1, n, my(c = 2 * r[k+1] - r[k] - r[k+2]); for(j = 1, c, listput(out, k)));
  Vec(out);
}
""")
pari(r"""
dp_forms(nf, q, l) = {
  \\ elements of E_0 of valuation 1 at one prime of E_0 above l and 0 at the others
  my(P = nf.pol, F = Mod(variable(P), P), s = F + q^2 / F, h = minpoly(s, 'x), nf0 = nfinit(h),
     dec0 = idealprimedec(nf0, l));
  vector(#dec0, i, subst(lift(nfbasistoalg(nf0, idealappr(nf0, matrix(#dec0, 2, a, b, if(b == 1, dec0[a], a == i))))), 'x, s));
}
""")
pari(r"""
dp_predict(P, q, sgn, l, prec) = {
  \\ the prediction (P).  Returns [m, c, kind, data...]
  my(t = variable(P), F = Mod(t, P), h = minpoly(F + q^2 / F, 'x), fa = factorpadic(h, l, prec), c = 1, d = 0,
     S, Pd, D, vD, w, out);
  for(i = 1, #fa~,
    if(valuation(subst(fa[i,1], 'x, 2 * sgn * q), l) > 0, c *= fa[i,1]^fa[i,2]; d += poldegree(fa[i,1]) * fa[i,2]));
  if(d != 2, return([2 * d]));
  S = -polcoef(c, 1) / (sgn * q); Pd = polcoef(c, 0) / q^2;
  D = Pd^2 - 4 * (S^2 - 2 * Pd - 4); vD = valuation(D, l);
  if(vD % 2, return([4, vD, "empty"]));
  if(issquare(D),
    my(sq = sqrt(D), tt = [(Pd + sq) / 2, (Pd - sq) / 2]);
    out = vector(2, i, w = valuation(tt[i] - 2, l);
      if(w % 2, ["ram", (w - 1) / 2], [if(issquare((tt[i] - 2) * (tt[i] + 2) / l^w), "split", "inert"), w / 2]));
    return([4, vD, "S", out[1], out[2]]));
  w = min(valuation(Pd - 4, l), vD / 2);
  if(w % 2, return([4, vD, "N", ["ram", (w - 1) / 2]]));
  [4, vD, "N", [if(issquare((S^2 - 4 * Pd) / l^(2 * w)), "split", "inert"), w / 2]];
}
""")
GP_primes, GP_box, GP_dual, GP_elsewhere, GP_jordan, GP_forms, GP_predict = (
    pari(n) for n in ("dp_primes", "dp_box", "dp_dual", "dp_elsewhere", "dp_jordan", "dp_forms", "dp_predict"))
INDEXVAL = pari("(A, B, l) -> valuation(abs(matdet(B) / matdet(A)), l)")


# ---------------------------------------------------------------------------------------------------
# graphs
# ---------------------------------------------------------------------------------------------------

def refine(adj, colour, rounds=8):
    """Colour refinement; returns the multiset of final colours."""
    for _ in range(rounds):
        colour = {v: hash((colour[v], tuple(sorted(colour[u] for u in adj[v])))) for v in adj}
    return tuple(sorted(Counter(colour.values()).items()))


def whole_hash(adj):
    return refine(adj, {v: len(adj[v]) for v in adj})


def ball(adj, root, radius):
    dist = {root: 0}
    queue = [root]
    for v in queue:
        if dist[v] < radius:
            for u in adj[v]:
                if u not in dist:
                    dist[u] = dist[v] + 1
                    queue.append(u)
    return dist


def ball_hash(adj, degree, root, radius):
    """Hash of the ball around root, each vertex coloured by its distance and its degree in the graph."""
    dist = ball(adj, root, radius)
    sub = {v: [u for u in adj[v] if u in dist] for v in dist}
    return refine(sub, {v: (dist[v], degree[v]) for v in dist})


def volcano(kind, r, Q, half):
    """Vol(K, r) with branching Q; for a split K the tops -half..half.  Returns (parity, level, adj)."""
    par, lev, adj = [], [], []

    def new(p, l):
        par.append(p)
        lev.append(l)
        adj.append([])
        return len(par) - 1

    def join(a, b):
        adj[a].append(b)
        adj[b].append(a)

    def grow(v, children):
        todo = [(v, children)]
        while todo:
            a, k = todo.pop()
            if lev[a] >= r:
                continue
            for _ in range(k):
                b = new(1 - par[a], lev[a] + 1)
                join(a, b)
                todo.append((b, Q))
    if kind == "inert":
        grow(new(0, 0), Q + 1)
    elif kind == "ram":
        a, b = new(0, 0), new(1, 0)
        join(a, b)
        grow(a, Q)
        grow(b, Q)
    else:
        tops = [new(k % 2, 0) for k in range(-half, half + 1)]
        for a, b in zip(tops, tops[1:]):
            join(a, b)
        for a in tops:
            grow(a, Q - 1)
    return par, lev, adj


def predicted_graph(pred, l, eps, half):
    """The graph of (P): adjacency, and the vertices whose two levels are 0 (centred, for split K)."""
    if pred[2] == "N":
        par, lev, adj = volcano(pred[3][0], pred[3][1], l * l, half)
        g = {i: list(a) for i, a in enumerate(adj)}
        centre = [i for i in g if lev[i] == 0]
        if pred[3][0] == "split":
            centre = [centre[len(centre) // 2]]
        return g, centre
    v1 = volcano(pred[3][0], pred[3][1], l, half)
    v2 = volcano(pred[4][0], pred[4][1], l, half)
    verts = [(i, j) for i in range(len(v1[0])) for j in range(len(v2[0])) if (v1[0][i] + v2[0][j]) % 2 == eps]
    g = {(i, j): [(a, b) for a in v1[2][i] for b in v2[2][j]] for i, j in verts}
    tops1 = [i for i in range(len(v1[0])) if v1[1][i] == 0]
    tops2 = [j for j in range(len(v2[0])) if v2[1][j] == 0]
    if pred[3][0] == "split":
        tops1 = tops1[len(tops1) // 2 - 1: len(tops1) // 2 + 1]
    if pred[4][0] == "split":
        tops2 = tops2[len(tops2) // 2 - 1: len(tops2) // 2 + 1]
    return g, [(i, j) for i in tops1 for j in tops2 if (i, j) in g]


def predict(P, q, sign, l):
    """(P) as a Python list: [m] or [4, c, kind, (type, r), ...]."""
    raw = GP_predict(P, q, sign, l, 60)
    if len(raw) == 1:
        return [int(raw[0])]
    return [int(raw[0]), int(raw[1]), str(raw[2])] + [(str(raw[i][0]), int(raw[i][1])) for i in range(3, len(raw))]


def is_infinite(pred):
    return any(x[0] == "split" for x in pred[3:])


# ---------------------------------------------------------------------------------------------------
# lattices
# ---------------------------------------------------------------------------------------------------

def iso_count(nu, l):
    """Possible numbers of isotropic lines of a non-degenerate form of rank nu over F_l."""
    return {0: {0}, 1: {0}, 2: {0, 2}, 3: {l + 1}, 4: {(l + 1) ** 2, l * l + 1},
            5: {l ** 3 + l * l + l + 1}, 6: {(l * l + l + 1) * (l * l + 1), (l ** 3 + 1) * (l + 1)}}[nu]


class Case:
    """One polynomial, one sign, with the l-adic data that does not depend on the form."""

    def __init__(self, q, e, l, sign):
        self.q, self.e, self.l, self.sign = q, e, l, sign
        self.P = cm.pol(cm.charpoly(q, e))
        self.nf = pari.nfinit(self.P)
        F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(self.nf)
        self.g = F - sign * q
        self.lam = sign * q % l
        self.LL, self.al0, units = GP_setup(self.nf, q, l, sign)
        self.units = list(units)
        self.prs, info = GP_primes(self.nf, q, l, sign)
        self.info = sorted([int(x) for x in row] for row in info)
        self.self = [int(row[2]) for row in info]
        self.R = GP_ring(self.nf, q, self.LL, 60, l)
        self.odd = list(GP_forms(self.nf, q, l))

    def forms(self, limit, skip=0):
        """Gram matrices of forms on O_E that are unimodular at the other primes above l: for each
        class of valuations, those of the units number skip to skip + limit - 1 that qualify."""
        twists = [1]
        for o in self.odd:
            twists += [t * o for t in twists]
        for o in twists:
            k = 0
            for u in self.units:
                G = sv.GP_gram(self.nf, self.q, self.al0 * u * o)
                if GP_elsewhere(self.nf, G, self.prs, self.l):
                    k += 1
                    if k > skip:
                        yield G
                    if k >= skip + limit:
                        break

    def box_lattices(self, G, w=1, cap=1500):
        """All unimodular F-stable lattices inside the box of dp_box; None if there are too many."""
        B = GP_box(self.nf, G, self.prs, self.self, self.R, self.l, w)
        if B == 0:
            return {}
        low = GP_dual(B, G, self.l)
        if not INSIDE(B, low):
            return {}
        idx = int(INDEXVAL(B, low, self.l))
        if idx % 2:
            return {}
        layer = {str(B): B}
        for _ in range(idx // 2):
            nxt = {}
            for X in layer.values():
                for Y in GP_maxsub(self.nf, X, self.g, self.l):
                    if INSIDE(Y, low):
                        nxt.setdefault(str(Y), Y)
            layer = nxt
            if len(layer) > cap:
                return None
        return {k: X for k, X in layer.items() if GP_unimodular(X, G, self.l)}

    def steps(self, L, G):
        """Neighbours of L, after testing the local rule (O) on L."""
        l = self.l
        O = sv.GP_order(self.nf, L)
        mu, lines = GP_steps(self.nf, L, G, O, l, self.lam)
        mu = int(mu)
        iso = [M for a, M, st in lines if int(a)]
        blocks = sorted(int(x) for x in GP_jordan(self.nf, L, self.g, l))
        cok = [int(x) for x in sv.GP_coker(self.nf, L, self.g, l)]
        assert len(blocks) == mu and len(cok) == mu and mu % 2 == 0, "(O) number of generators"
        assert cok[0::2] == cok[1::2], "(O) the cokernel is not of the shape A x A"
        assert all(blocks.count(s) % 2 == 0 for s in set(blocks) if s % 2 == 0), "(O) even blocks are not paired"
        nu = blocks.count(1)
        rho = mu - nu
        allowed = {(l ** rho - 1) // (l - 1) + l ** rho * i for i in iso_count(nu, l)}
        assert len(iso) in allowed, "(O) number of rational steps"
        return iso, (tuple(blocks), tuple(cok))

    def explore(self, G, seed, radius=None, cap=6000):
        """Follow the steps from seed: the whole component, or the ball of the given radius.
        Returns (adjacency, degree, invariants); adjacency is restricted to the vertices met."""
        verts = {str(seed): seed}
        dist = {str(seed): 0}
        nbrs, inv = {}, {}
        queue = [str(seed)]
        for k in queue:
            iso, inv[k] = self.steps(verts[k], G)
            nbrs[k] = [str(M) for M in iso]
            if radius is not None and dist[k] >= radius:
                continue
            for M in iso:
                if str(M) not in verts:
                    verts[str(M)] = M
                    dist[str(M)] = dist[k] + 1
                    queue.append(str(M))
            if len(verts) > cap:
                return None
        adj = {k: [x for x in nbrs[k] if x in nbrs] for k in nbrs}
        return adj, {k: len(nbrs[k]) for k in nbrs}, inv


def check_rank4(case, pred, forms_limit, radius):
    """Compare the lattices of one polynomial with the prediction (P)."""
    l = case.l
    kind = pred[2]
    found = Counter()
    lattices = 0
    first = list(case.forms(forms_limit))
    more = case.forms(10 ** 6, skip=forms_limit) if kind != "empty" else []
    for number, G in enumerate(itertools.chain(first, more)):
        if number >= len(first) and any(k.startswith("eps") for k in found):
            break                                      # further forms are tried only until a lattice is found
        box = case.box_lattices(G)
        if box is None:                                # too many lattices to list: take one, if any
            seed = lagrangian_seed(case, G)
            if seed is None:
                found["not listed"] += 1
                continue
        elif not box:
            found["no lattice"] += 1
            continue
        else:
            seed = next(iter(box.values()))
        assert kind != "empty", "a unimodular lattice exists although v(D) is odd"
        if not is_infinite(pred):
            res = case.explore(G, seed)
            assert res is not None
            adj, degree, inv = res
            lattices += len(adj)
            assert box is None or set(box) == set(adj), "the lattices reached by steps are not all the lattices"
            got = whole_hash(adj)
            ok = [eps for eps in ((0, 1) if kind == "S" else (0,))
                  if whole_hash(predicted_graph(pred, l, eps, 0)[0]) == got]
            assert ok, "the graph is not the predicted one"
            found[f"eps={ok}"] += 1
        else:
            res = case.explore(G, seed, radius=radius)
            assert res is not None
            adj, degree, inv = res
            lattices += len(adj)
            got = ball_hash(adj, degree, str(seed), radius)
            rootrad = pred[3][1] + (pred[4][1] if kind == "S" else 0) + 2
            half = rootrad + radius + 2
            ok = []
            for eps in ((0, 1) if kind == "S" else (0,)):
                pg, centre = predicted_graph(pred, l, eps, half)
                deg = {v: len(pg[v]) for v in pg}
                roots = set()
                for c in centre:
                    roots |= set(ball(pg, c, rootrad))
                if any(ball_hash(pg, deg, v, radius) == got for v in roots):
                    ok.append(eps)
            assert ok, "the ball is not a ball of the predicted graph"
            found[f"eps={ok}"] += 1
    if kind == "empty":
        assert set(found) <= {"no lattice", "not listed"}
    else:
        assert any(k != "no lattice" for k in found), "no form with a unimodular lattice was found"
    return dict(found), lattices


def lagrangian_seed(case, G):
    """A unimodular lattice N + (maximal isotropic subspace of N^dual / N), N = O_E; None if there is none."""
    l = case.l
    N = pari.matid(6)
    D = GP_at(G ** (-1), N, l, 60)
    if not INSIDE(D, N):
        return None
    k = int(INDEXVAL(D, N, l))
    if k == 0:
        return N
    if k % 2:
        return None
    cols = [D[i] for i in range(6)]
    gens = []
    X = N
    for v in cols:
        Y = CONCAT(X, pari.Mat(v))
        if Y != X:
            gens.append(v)
            X = Y
    if len(gens) != k:                                # N^dual / N is not killed by l
        return None
    L = N
    chosen = 0
    for coeffs in itertools.product(range(l), repeat=k):
        if not any(coeffs):
            continue
        x = sum((c * v for c, v in zip(coeffs, gens)), start=0 * gens[0])
        Y = CONCAT(L, pari.Mat(x))
        if Y == L:
            continue
        GY = pari("(Y, G, l) -> my(M = Y~ * G * Y, d = denominator(M)); d % l != 0")(Y, G, l)
        if GY:
            L = Y
            chosen += 1
            if chosen == k // 2:
                stable = INSIDE(L, rm.MUL(case.nf, case.g) * L)
                return L if stable and GP_unimodular(L, G, l) else None
    return None


def survey_rank6(case, forms_limit, cap):
    """Follow the steps from one lattice and describe what is found."""
    for G in case.forms(forms_limit):
        seed = lagrangian_seed(case, G)
        if seed is None:
            continue
        infinite = not all(case.self)
        res = case.explore(G, seed, radius=3 if infinite else None, cap=cap)
        if res is None:
            return {"lattices": f"more than {cap}"}
        adj, degree, inv = res
        shape = Counter((inv[k][0], inv[k][1], degree[k]) for k in adj)
        return {"lattices": len(adj), "closed": not infinite,
                "shape": [[list(b), list(c), d, n] for (b, c, d), n in sorted(shape.items())]}
    return {"lattices": 0}


def candidates(q, l, sign, cmin):
    """Cubics with three real roots in (-2q, 2q) and l^cmin dividing h(2 sign q); yields (e, c)."""
    x0, B = 2 * sign * q, 2 * q
    mod = l ** cmin
    for a in range(-3 * B + 1, 3 * B):
        for b in range(-3 * B * B + 1, 3 * B * B):
            if a * a - 3 * b <= 0:
                continue
            if not (3 * B * B - 2 * a * B + b > 0 and 3 * B * B + 2 * a * B + b > 0):
                continue
            base = x0 ** 3 - a * x0 ** 2 + b * x0
            for c in range(-B ** 3 + ((base + B ** 3) % mod), B ** 3, mod):
                hv = base - c
                if hv == 0:
                    continue
                if not (-B ** 3 - a * B * B - b * B - c < 0 < B ** 3 - a * B * B + b * B - c):
                    continue
                if a * a * b * b - 4 * b ** 3 - 4 * a ** 3 * c + 18 * a * b * c - 27 * c * c <= 0:
                    continue
                yield [a, b + 3 * q * q, c + 2 * a * q * q], val(hv, l)


def describe(pred):
    if pred[2] == "empty":
        return "none"
    parts = [f"{x[0]} r={x[1]}" for x in pred[3:]]
    return ("product: " if pred[2] == "S" else "over Q_l^2: ") + " x ".join(sorted(parts))


def search(q, l, args, fh, total):
    seen6 = Counter()
    for sign in ((1, -1) if args.both_signs else (1,)):
        low = min(args.cmin, 2) if args.u_every else args.cmin
        count2 = 0
        for e, c in candidates(q, l, sign, low):
            if c % 2 or (c < args.cmin and c != 2):
                continue
            P = cm.pol(cm.charpoly(q, e))
            if not pari.polisirreducible(P):
                continue
            pred = predict(P, q, sign, l)
            m = pred[0]
            if m == 2 or (c == 2 and m != 4):
                continue
            if c == 2:
                count2 += 1
                if count2 % args.u_every != 1 and args.u_every != 1:
                    continue
            nf0 = pari.nfinit(P)
            dec = [[int(x) for x in row] for row in sv.GP_dec(nf0, q, l)]
            if any(p[3] == 2 for p in dec if p[4 if sign == 1 else 5] == 0):
                continue
            case = Case(q, e, l, sign)
            if m == 4:
                assert pred[1] == c, "v(D) is not the exponent of l in P(lam)"
                try:
                    found, lattices = check_rank4(case, pred, args.forms, args.radius)
                except AssertionError:
                    print(f"FAILED: q={q} e={e} l={l} sign={sign} c={c} primes {case.info} predicted {pred}", flush=True)
                    raise
                key = f"rank 4, c={c}, {describe(pred)}"
                total[key] += 1
                total["lattices tested"] += lattices
                rec = {"q": q, "e": e, "l": l, "eigenvalue": "+q" if sign == 1 else "-q", "rank": 4, "c": c,
                       "primes": case.info, "predicted": describe(pred), "found": found, "lattices": lattices}
                if c >= 4:
                    print(f"q={q} e={e} l={l} {rec['eigenvalue']} c={c} primes {case.info}: {describe(pred)}; "
                          f"forms {found}; {lattices} lattices; agrees", flush=True)
            else:
                key = (c, tuple(map(tuple, case.info)))
                seen6[key] += 1
                if seen6[key] > args.rank6 or c != 4:
                    continue
                res = survey_rank6(case, args.forms, args.cap)
                total["rank 6 polynomials followed"] += 1
                rec = {"q": q, "e": e, "l": l, "eigenvalue": "+q" if sign == 1 else "-q", "rank": 6, "c": c,
                       "primes": case.info, **res}
                print(f"q={q} e={e} l={l} {rec['eigenvalue']} c={c} rank 6 primes {case.info}: {res}", flush=True)
            if fh:
                fh.write(json.dumps(rec) + "\n")
    print(f"search q={q} l={l}: rank 6 polynomials by (c, primes [e, f, self-conjugate, v]): "
          f"{ {str(k): v for k, v in sorted(seen6.items())} }", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    ap.add_argument("--search", action="append", help="q,l (repeatable); default 5,3")
    ap.add_argument("--cmin", type=int, default=4)
    ap.add_argument("--u-every", type=int, default=25, help="test one in K of the rank-4 polynomials with c = 2 (0: none)")
    ap.add_argument("--rank6", type=int, default=2, help="polynomials followed per configuration in rank 6")
    ap.add_argument("--forms", type=int, default=3, help="forms tried per class of valuations")
    ap.add_argument("--radius", type=int, default=2, help="radius of the ball compared when the graph is infinite")
    ap.add_argument("--cap", type=int, default=3000)
    ap.add_argument("--both-signs", action="store_true",
                    help="also eigenvalue -q (the polynomials P(-t), with the same structure)")
    args = ap.parse_args()
    fh = open(args.out, "w") if args.out else None
    total = Counter()
    for item in args.search or ["5,3"]:
        q, l = (int(x) for x in item.split(","))
        search(q, l, args, fh, total)
    print("summary:", {k: v for k, v in sorted(total.items())})


if __name__ == "__main__":
    main()
