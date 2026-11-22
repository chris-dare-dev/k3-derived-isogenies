#!/usr/bin/env python
"""
ellk3_two_neighbours.py -- the two Jacobians of a double cover of P^1 x P^1 branched along a curve of
bidegree (4, 4), over F_p, against the lattice prediction for a step at 2 (note 22; the lead of note 17, §5.3).

Y: w^2 = f(s, u), f of bidegree (4, 4), has two genus-one fibrations, over the u-line and over the s-line.
Their relative Jacobians are the elliptic K3 surfaces

    X1: y^2 = x^3 - 27 I1(u) x - 27 J1(u),      X2: y^2 = x^3 - 27 I2(s) x - 27 J2(s),

with I, J the invariants of f as a binary quartic in the other variable.  Over C the transcendental lattices
of X1 and X2 are the two unimodular overlattices of T(Y) (van Geemen), so X1 and X2 should be 2-neighbours.

Families of f:
  generic   all 25 coefficients random: X_i has 24 fibres I_1, NS = U, transcendental rank 20;
  e8        the branch curve has a singular point of type E8 at (0, 0) with tangent s = u: X_i has a fibre
            II* at 0 and 14 fibres I_1, NS = U + E8, transcendental rank 12 (geometric Picard number 10);
  e8a2      in addition a cusp at (infinity, infinity) with a random tangent: X_i has II* at 0, I_3 at
            infinity and 11 fibres I_1, NS = U + E8 + A2, transcendental rank 10 (geometric Picard number 12).

For each f the script
  1. builds X1, X2 and checks the fibre types from the discriminant;
  2. checks, fibre by fibre over F_p, that the genus-one curve of Y and the fibre of X_i have the same
     number of points;
  3. computes the L-function of X1 (degree N = 20, 12 or 10) from the traces of the fibres over F_{p^d},
     d <= N/2, and the functional equation; and the same sums for X2 up to --check2;
  4. tests whether X1 and X2 are isomorphic over F_p as Weierstrass fibrations (all of PGL_2(F_p));
  5. on the lattice side, from P alone: irreducibility, ordinarity, the order Z[F, V] at p, at 2 and at
     other primes, the primes above 2.  When there is exactly one pair l, conj(l) of primes of degree 1
     above 2 and the order is forced to be maximal (condition C), computes whether l / conj(l) is generated
     by an element of relative norm 1, the group G = Cl(E) / <Cl(E0), ramified primes> of note 13, Lemma 9,
     and the order of the class of l in G.

Prediction tested: under condition C,   X1 = X2 over F_p   if and only if   l / conj(l) = (x), x conj(x) = 1.

  control   random Weierstrass equations of the same fibre types that do not come from a double cover:
            the same lattice data, for comparison of what happens above 2.

PARI's bnfinit is used without certification, so class groups are correct under GRH.

  python -I ellk3_two_neighbours.py --family e8|e8a2|generic --p 5 --count 200 [--seed 1] [--workers 12] [--out FILE]
                                    [--control] [--check2 D] [--bnf-limit SECONDS]
  python -I ellk3_two_neighbours.py --graph FILE
  python -I ellk3_two_neighbours.py --summary FILE     (the summary of a run, from its result file)

  --graph   for a result file of the family e8a2: the graph of the pairs (X1, X2) for each polynomial.  If X2 is
            always the image of X1 under the class of l / conj(l) or its inverse, a surface has at most two
            distinct partners and every closed cycle has length the order of that class.

Run from experiments/.
"""
import argparse
import json
import random
import sys
import time
from collections import Counter
from math import comb
from multiprocessing import Pool

import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9, silent=True)
for _v in ("s", "u", "v", "U", "V", "t", "T", "x", "z"):
    pari(_v)

GP_SOURCE = r"""
ep_wei(f, p) = {
  \\ f in Z[s, u] of bidegree (4, 4).  [A1, B1, A2, B2] in F_p[u]: the Jacobians of the fibrations over the
  \\ u-line and over the s-line are y^2 = x^3 + A x + B.
  my(q1 = vector(5, i, polcoef(f, 5 - i, 's)), q2 = vector(5, i, subst(polcoef(f, 5 - i, 'u), 's, 'u)), r = List());
  foreach([q1, q2], q,
    my(a = q[1] * Mod(1, p), b = q[2] * Mod(1, p), c = q[3] * Mod(1, p), d = q[4] * Mod(1, p), e = q[5] * Mod(1, p), I, J);
    I = 12*a*e - 3*b*d + c^2; J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3;
    listput(r, -27 * I); listput(r, -27 * J));
  Vec(r);
}

ep_types(A, B, p, fam) = {
  \\ 1 if the fibres are as expected: 24 I_1 (fam 0); II* at 0 and 14 I_1, smooth at infinity (fam 1);
  \\ II* at 0, I_3 at infinity and 11 I_1 (fam 2).
  my(D = 4*A^3 + 27*B^2);
  if(D == 0 || poldegree(D) != if(fam == 2, 21, 24), return(0));
  if(fam == 2 && (poldegree(A) != 8 || poldegree(B) != 12), return(0));
  if(fam == 0, return(issquarefree(D)));
  if(A != 0 && valuation(A, 'u) < 4, return(0));
  if(B == 0 || valuation(B, 'u) != 5 || valuation(D, 'u) != 10, return(0));
  issquarefree(D / 'u^10);
}

ep_tr(a, b) = {
  \\ trace of Frobenius on the fibre y^2 = x^3 + a x + b over the field of a and b: good, node (+1 split, -1 not), cusp (0).
  if(4*a^3 + 27*b^2 != 0, return(ellap(ellinit([a, b]))));
  if(a == 0, 0, if(issquare(-2*a*b), 1, -1));
}

ep_ev(v, x, one) = { my(r = v[#v] * one); forstep(i = #v - 1, 1, -1, r = r * x + v[i]); r; }

ep_S(A, B, p, d) = {
  \\ sum over the points of P^1(F_{p^d}) of the trace of the fibre.
  my(g = ffgen([p, d], 'z), one = g^0, r = ffprimroot(g), q = p^d, x = one, S,
     va = concat(Vecrev(lift(A)), vector(8 - poldegree(A))), vb = concat(Vecrev(lift(B)), vector(12 - poldegree(B))));
  S = ep_tr(va[9] * one, vb[13] * one) + ep_tr(va[1] * one, vb[1] * one);
  for(k = 1, q - 1, S += ep_tr(ep_ev(va, x, one), ep_ev(vb, x, one)); x = x * r);
  S;
}

ep_Sall(A, B, p, m) = {
  \\ [S_1, ..., S_m] as ep_S, from one point in each Frobenius orbit: a closed point of degree e with trace a
  \\ contributes e (alpha^j + beta^j) to S_{e j}, alpha + beta = a, alpha beta = p^e (a node: e (+-1)^j; a cusp: 0).
  my(S = vector(m), va = concat(Vecrev(lift(A)), vector(8 - poldegree(A))), vb = concat(Vecrev(lift(B)), vector(12 - poldegree(B))), pts, a, b, tr, t0, t1, t2);
  for(e = 1, m,
    my(g = ffgen([p, e], 'z), one = g^0, r = ffprimroot(g), q = p^e, x = one, k2, mn);
    for(k = if(e == 1, -2, 0), q - 2,
      \\ k = -2, -1: the points infinity and 0.  Otherwise x = r^k, kept if k is the least element of its orbit
      \\ under multiplication by p mod q - 1 and the orbit has e elements.
      if(k == -2, a = va[9] * one; b = vb[13] * one, k == -1, a = va[1] * one; b = vb[1] * one,
        k2 = k; mn = 1;
        for(i = 1, e - 1, k2 = k2 * p % (q - 1); if(k2 <= k, mn = 0; break));
        if(mn, a = ep_ev(va, x, one); b = ep_ev(vb, x, one));
        x = x * r;
        if(!mn, next));
      tr = ep_tr(a, b);
      if(4*a^3 + 27*b^2 != 0,
        t0 = 2; t1 = tr; for(j = 1, m \ e, S[e*j] += e * t1; t2 = tr * t1 - q * t0; t0 = t1; t1 = t2),
        for(j = 1, m \ e, S[e*j] += e * tr^j))));
  S;
}

ep_fibres(f, A, B, p, which) = {
  \\ over F_p: for each u with a smooth fibre, the number of points of w^2 = f(., u) (which = 1) or w^2 = f(u, .)
  \\ (which = 2) on its smooth model must be p + 1 - trace.  Returns [number of fibres compared, number that agree].
  my(n = 0, ok = 0, chi = (c -> if(c == 0, 0, if(issquare(c), 1, -1))));
  for(u0 = 0, p - 1,
    my(a = subst(A, 'u, Mod(u0, p)), b = subst(B, 'u, Mod(u0, p)), g, cnt);
    if(4*a^3 + 27*b^2 == 0, next);
    g = if(which == 1, subst(f, 'u, u0), subst(subst(f, 's, u0), 'u, 's)) * Mod(1, p);
    cnt = 1 + chi(polcoef(g, 4, 's)) + sum(s0 = 0, p - 1, 1 + chi(subst(g, 's, Mod(s0, p))));
    n++; if(cnt == p + 1 - ep_tr(a, b), ok++));
  [n, ok];
}

ep_pol(Sv, p) = {
  \\ P(t) = prod (t - alpha_i) from S_1..S_m, N = 2m, with the functional equation of sign +.  0 if the middle
  \\ coefficient vanishes (sign not determined).
  my(m = #Sv, N = 2*m, L = exp(sum(d = 1, m, Sv[d] * 'T^d / d) + O('T^(m + 1))), c = vector(m + 1, i, polcoef(L, i - 1)), w);
  if(denominator(c) != 1, error("L-series not integral"));
  if(c[m + 1] == 0, return(0));
  w = vector(N + 1, i, if(i <= m + 1, c[i], p^(2*(i - 1) - N) * c[N + 2 - i]));
  sum(i = 1, N + 1, w[i] * 't^(N + 1 - i));
}

ep_weil(P, p) = vecmax(apply(z -> abs(abs(z) - p), polroots(P))) < 1e-9 * p;

ep_iso(A1, B1, A2, B2, p) = {
  \\ 1 if (A2, B2) o M = (mu^4 A1, mu^6 B1) for some M in PGL_2(F_p) and mu in F_p^*; 2 if only up to a
  \\ quadratic twist; 0 otherwise.
  my(h = (A, n) -> sum(i = 0, n, polcoef(A, i) * 'u^i * 'v^(n - i)), a1 = h(A1, 8), b1 = h(B1, 12), a2 = h(A2, 8), b2 = h(B2, 12),
     res = 0, one = Mod(1, p),
     try = (m11, m12, m21, m22) ->
       my(lu = (m11 * 'u + m12 * 'v) * one, lv = (m21 * 'u + m22 * 'v) * one, a = subst(subst(a2, 'u, 'U), 'v, 'V), b = subst(subst(b2, 'u, 'U), 'v, 'V), ka, kb, lam);
       a = subst(subst(a, 'U, lu), 'V, lv); b = subst(subst(b, 'U, lu), 'V, lv);
       if(b1 == 0 || b == 0, return(0));
       kb = 0; for(i = 0, 12, my(c1 = polcoef(polcoef(b1, i, 'u), 12 - i, 'v)); if(c1 != 0, kb = polcoef(polcoef(b, i, 'u), 12 - i, 'v) / c1; break));
       if(kb == 0 || b != kb * b1, return(0));
       if(a1 == 0, return(if(a == 0, -1, 0)));
       ka = 0; for(i = 0, 8, my(c1 = polcoef(polcoef(a1, i, 'u), 8 - i, 'v)); if(c1 != 0, ka = polcoef(polcoef(a, i, 'u), 8 - i, 'v) / c1; break));
       if(ka == 0 || a != ka * a1, return(0));
       lam = kb / ka; if(lam^2 != ka, return(0));
       if(issquare(lam), 1, 2));
  for(b = 0, p - 1, for(c = 0, p - 1, for(d = 0, p - 1,
    if((d - b*c) % p, my(r = try(1, b, c, d)); if(r == 1, return(1)); if(r == 2, res = 2)))));
  for(d = 0, p - 1, my(r = try(0, 1, 1, d)); if(r == 1, return(1)); if(r == 2, res = 2));
  res;
}

ep_local(P, p) = {
  \\ From P alone.  [status, m2, sq2, cat, idx, odd, nf]:
  \\  status  "ok", or why the polynomial is set aside;
  \\  m2      multiplicity of t - 1 in P mod 2;  sq2  1 if the other factors of P mod 2 are simple;
  \\  cat     how the factor (t - 1)^m2 decomposes: "none", "split" (two primes of degree 1), "ramified", "inert", "multi";
  \\  idx     [O_E : Z[F, V]];  odd  its part prime to 2p.
  my(N = poldegree(P), t = variable(P), c, hb, f2, m2 = 0, sq2 = 1, nf, F, V, M, idx, odd, dec, d1, cat);
  if(!polisirreducible(P), return(["reducible", 0, 0, "", 0, 0, 0]));
  c = vector(N, k, polcoef(P, N - k));
  if(c[1] % p == 0, return(["not ordinary", 0, 0, "", 0, 0, 0]));
  hb = sum(k = 1, N - 1, (c[k] / p^(k - 1)) * t^(N - 1 - k));
  if(denominator(content(hb)) != 1, error("Newton polygon"));
  f2 = factor(P * Mod(1, 2));
  for(i = 1, matsize(f2)[1], if(poldegree(f2[i, 1]) == 1, m2 = f2[i, 2], if(f2[i, 2] > 1, sq2 = 0)));
  nf = nfinit([P, 10^6]);
  if(#nfcertify(nf), return(["uncertified", m2, sq2, "", 0, 0, 0]));
  F = Mod(t, P); V = p^2 / F;
  M = Mat(concat(vector(N, i, nfalgtobasis(nf, F^(i - 1))), vector(N - 1, j, nfalgtobasis(nf, V^j))));
  idx = abs(matdet(mathnf(M)));
  odd = idx / (2^valuation(idx, 2) * p^valuation(idx, p));
  dec = idealprimedec(nf, 2); d1 = select(pr -> pr.f == 1, dec);
  cat = if(m2 == 0, "none", m2 > 2, "multi", #d1 == 2 && d1[1].e == 1 && d1[2].e == 1, "split", #d1 == 1, "ramified", #d1 == 0, "inert", "other");
  [if(!issquarefree(hb * Mod(1, p)), "p not maximal", "ok"), m2, sq2, cat, idx, odd, nf];
}

ep_global(nf, p) = {
  \\ Under condition C.  [l/conj(l) generated by an element of relative norm 1, |G|, elementary divisors of G,
  \\ order of the class of l in G, h(E), h(E0), order of the class of l/conj(l) in Cl(E)]
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = p^2 / F, sig = lift(V), cj = (z -> nfgaloisapply(nf, sig, z)),
     d1 = select(pr -> pr.f == 1, idealprimedec(nf, 2)), l = idealhnf(nf, d1[1]), lb = cj(l), bnf, I, pr, y, uu, r, Mn, vu, e, ctriv = 0,
     s = F + V, h, bnf0, up, cyc, k, Mm, ram, H, D, vl, w, ordG, Gsize, ordcl);
  if(lb != idealhnf(nf, d1[2]) || lb == l, error("the two primes above 2 are not conjugate"));
  bnf = bnfinit(nf, 1);
  I = idealdiv(nf, l, lb);
  pr = bnfisprincipal(bnf, I);
  ordcl = if(#bnf.cyc == 0, 1, lcm(vector(#bnf.cyc, i, bnf.cyc[i] / gcd(bnf.cyc[i], pr[1][i]))));
  if(pr[1] == 0,
    y = pr[2]; uu = nfeltmul(nf, y, cj(y)); r = #bnf.fu;
    vu = bnfisunit(bnf, uu); if(#vu == 0, error("y conj(y) is not a unit"));
    if(r == 0, ctriv = (nfalgtobasis(nf, uu) == nfalgtobasis(nf, 1)),
      Mn = matrix(r, r);
      for(i = 1, r, my(fi = nfalgtobasis(nf, bnf.fu[i]), ni = bnfisunit(bnf, nfeltmul(nf, fi, cj(fi)))); for(j = 1, r, Mn[j, i] = ni[j]));
      e = matsolve(Mn, vector(r, j, vu[j])~);
      ctriv = (denominator(e) == 1);
      if(ctriv,
        my(xx = nfeltdiv(nf, y, nffactorback(nf, bnf.fu, e)));
        if(nfalgtobasis(nf, nfeltmul(nf, xx, cj(xx))) != nfalgtobasis(nf, 1) || idealhnf(nf, xx) != I, error("norm-one generator")))));
  h = minpoly(s, 'x); bnf0 = bnfinit(h, 1);
  up = (J -> my(two = idealtwoelt(bnf0, J)); idealhnf(nf, two[1], subst(lift(nfbasistoalg(bnf0, two[2])), 'x, s)));
  cyc = bnf.cyc; k = #cyc;
  if(k == 0, return([ctriv, 1, [], 1, 1, bnf0.no, ordcl]));
  Mm = matdiagonal(cyc);
  for(j = 1, #bnf0.gen, Mm = concat(Mm, Mat(bnfisprincipal(bnf, up(bnf0.gen[j]), 0))));
  ram = abs(nf.disc) / bnf0.disc^2;
  foreach(factor(ram)[, 1]~, ll,
    foreach(idealprimedec(bnf0, ll), q0,
      my(fa = idealfactor(nf, up(q0)));
      if(matsize(fa)[1] == 1 && fa[1, 2] == 2, Mm = concat(Mm, Mat(bnfisprincipal(bnf, fa[1, 1], 0))))));
  H = mathnf(Mm); Gsize = matdet(H); D = select(d -> d > 1, matsnf(H));
  vl = bnfisprincipal(bnf, l, 0); w = H^(-1) * vl; ordG = denominator(w);
  [ctriv, Gsize, D, ordG, bnf.no, bnf0.no, ordcl];
}

ep_globalt(nf, p, lim) = alarm(lim, ep_global(nf, p));
"""
for _chunk in GP_SOURCE.strip().split("\n\n"):
    pari(_chunk)   # one definition per call

GP_S = pari("ep_S")
GP_Sall = pari("ep_Sall")
GP_fibres = pari("ep_fibres")
GP_tr = pari("ep_tr")
GP_globalt = pari("ep_globalt")
GP_iso = pari("ep_iso")
GP_local = pari("ep_local")
GP_pol = pari("ep_pol")
GP_types = pari("ep_types")
GP_wei = pari("ep_wei")
GP_weil = pari("ep_weil")

FAMILIES = {"generic": (0, 20), "e8": (1, 12), "e8a2": (2, 10)}
E8_LOW = [(0, j) for j in range(5)] + [(1, j) for j in range(4)] + [(2, 0), (2, 1)]


A2_LOW = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]


def e8_space(p, cusp=None):
    """Basis over F_p of the forms f(s, u) = sum c[i][j] s^i u^j, 0 <= i, j <= 4, such that g(t, u) = f(t + u, u)
    has no monomial t^a u^b with a/3 + b/5 < 1: a point of type E8 (or worse) at the origin, tangent s = u.
    With cusp = lam: in addition, in the coordinates 1/s, 1/u at (infinity, infinity), g'(t', u') = f'(t' + lam u', u')
    has no monomial t'^a u'^b with a/2 + b/3 < 1: a cusp (or worse) there, with tangent 1/s = lam / u."""
    rows = []
    for a, b in E8_LOW:
        row = [0] * 25
        for i in range(a, 5):
            j = b - (i - a)
            if 0 <= j <= 4:
                row[5 * i + j] = comb(i, a) % p
        rows.append(row)
    if cusp is not None:
        for a, b in A2_LOW:
            row = [0] * 25
            for i in range(a, 5):
                j = b - (i - a)
                if 0 <= j <= 4:
                    row[5 * (4 - i) + (4 - j)] = comb(i, a) * pow(cusp, i - a, p) % p
            rows.append(row)
    ker = pari.matker(pari.matrix(len(rows), 25, [x for r in rows for x in r]) * pari.Mod(1, p))
    return [[int(pari.lift(ker[c][r])) for r in range(25)] for c in range(len(ker))]


def sample_form(rng, p, family, basis):
    if family == "generic":
        return [rng.randrange(p) for _ in range(25)]
    c = [0] * 25
    for vec in basis:
        k = rng.randrange(p)
        c = [(x + k * y) % p for x, y in zip(c, vec)]
    return c


def form_pari(c):
    return pari("+".join(f"{c[5 * i + j]}*s^{i}*u^{j}" for i in range(5) for j in range(5)))


def coeffs(P, deg):
    return [int(pari.lift(pari.polcoef(P, k))) for k in range(deg + 1)]


def sample_weierstrass(rng, p, family):
    """A random equation with the fibre types of the family, not coming from a double cover.  For e8a2 the
    fibre I_3 at infinity is obtained by rejection in ep_types."""
    if family == "generic":
        a = [rng.randrange(p) for _ in range(9)]
        b = [rng.randrange(p) for _ in range(13)]
    else:
        a = [0] * 4 + [rng.randrange(p) for _ in range(5)]
        b = [0] * 5 + [rng.randrange(p) for _ in range(8)]
    A = pari("+".join(f"{x}*u^{i}" for i, x in enumerate(a))) * pari.Mod(1, p)
    B = pari("+".join(f"{x}*u^{i}" for i, x in enumerate(b))) * pari.Mod(1, p)
    return A, B


def label(A, B, p):
    """Least image of the coefficients of (A, B) under u -> c u and (A, B) -> (mu^4 A, mu^6 B): a canonical
    form for the family e8a2, where an isomorphism over F_p fixes 0 (the fibre II*) and infinity (the fibre I_3)."""
    best = None
    for c in range(1, p):
        for mu in range(1, p):
            m4, m6 = pow(mu, 4, p), pow(mu, 6, p)
            cand = (tuple(m4 * pow(c, i, p) * a % p for i, a in enumerate(A)), tuple(m6 * pow(c, i, p) * b % p for i, b in enumerate(B)))
            if best is None or cand < best:
                best = cand
    return best


def graph(path):
    """For the records of a run of the family e8a2: for each polynomial, the graph whose vertices are the surfaces
    X1, X2 up to isomorphism over F_p and whose edges are the pairs (X1, X2).  If X2 is the image of X1 under the
    class c of l / conj(l) or its inverse, the graph is a subgraph of the Cayley graph of c on the torsor: degree
    at most 2, and every cycle has length the order of c."""
    recs = [json.loads(line) for line in open(path)]
    recs = [r for r in recs if "P" in r and r["status"] == "ok" and r["sq2"] == 1 and r["index_odd"] == 1 and "A2" in r]
    out = []
    for name, sel in (("condition C", lambda r: r["cond_C"]), ("several pairs above 2 (multi)", lambda r: r["cat2"] == "multi")):
        byP = {}
        for r in recs:
            if sel(r):
                byP.setdefault(tuple(r["P"]), []).append(r)
        nsurf = nedges = loops = 0
        degs, deg_over, cyc_ok, cyc_bad, over_G, sat = Counter(), 0, 0, 0, 0, 0
        for P, rs in byP.items():
            p = rs[0]["p"]
            adj = {}
            for r in rs:
                a, b = label(r["A1"], r["B1"], p), label(r["A2"], r["B2"], p)
                if a == b:
                    loops += 1
                adj.setdefault(a, set()).add(b)
                adj.setdefault(b, set()).add(a)
            nsurf += len(adj)
            nedges += sum(len(v) for v in adj.values()) // 2
            bound = 2 if name == "condition C" else rs[0]["m2"]
            for v, nb in adj.items():
                degs[len(nb)] += 1
                if len(nb) > bound:
                    deg_over += 1
            if name == "condition C" and "G_size" in rs[0]:
                if len(adj) > rs[0]["G_size"]:
                    over_G += 1
                if len(adj) == rs[0]["G_size"]:
                    sat += 1
                # components that are cycles: every vertex of degree 2 and as many edges as vertices
                seen = set()
                for v in adj:
                    if v in seen:
                        continue
                    comp, stack = set(), [v]
                    while stack:
                        w = stack.pop()
                        if w not in comp:
                            comp.add(w)
                            stack.extend(adj[w])
                    seen |= comp
                    if len(comp) > 2 and all(len(adj[w]) == 2 for w in comp):
                        if len(comp) == rs[0]["ord_G"]:
                            cyc_ok += 1
                        else:
                            cyc_bad += 1
        out.append(f"{name}: {sum(len(v) for v in byP.values())} pairs, {len(byP)} polynomials, {nsurf} surfaces, {nedges} distinct edges, {loops} pairs with X1 = X2")
        out.append(f"  number of distinct partners of a surface: {dict(sorted(degs.items()))}; above the bound: {deg_over}")
        if name == "condition C":
            out.append(f"  polynomials with more surfaces than |G|: {over_G}; with exactly |G|: {sat}")
            out.append(f"  closed cycles: {cyc_ok} of length the order of the class, {cyc_bad} of another length")
    return out


def lattice_side(P, p, bnf_limit):
    loc = GP_local(P, p)
    rec = {"status": str(loc[0]), "m2": int(loc[1]), "sq2": int(loc[2]), "cat2": str(loc[3]),
           "index": int(loc[4]), "index_odd": int(loc[5])}
    rec["cond_C"] = bool(rec["status"] == "ok" and rec["cat2"] == "split" and rec["sq2"] == 1 and rec["index_odd"] == 1)
    if rec["cond_C"]:
        g = GP_globalt(loc[6], p, bnf_limit)
        if str(pari.type(g)) == "t_ERROR":
            rec["global_error"] = str(g)[:120]
        else:
            rec.update({"c_trivial": bool(g[0]), "G_size": int(g[1]), "G_div": [int(d) for d in g[2]],
                        "ord_G": int(g[3]), "hE": int(g[4]), "hE0": int(g[5]), "ord_cl": int(g[6])})
    return rec


def one(job):
    p, family, seed, control, check2, bnf_limit = job
    rng = random.Random(seed)
    fam = FAMILIES[family][0]
    N = FAMILIES[family][1]
    m = N // 2
    basis = e8_space(p) if fam == 1 else None
    rec = {"p": p, "family": family, "seed": seed, "control": control}
    t0 = time.time()
    for attempt in range(1, 1000001):
        if fam == 2 and not control:
            lam = rng.randrange(1, p)
            basis = e8_space(p, lam)
            rec["cusp_tangent"] = lam
        if control:
            A1, B1 = sample_weierstrass(rng, p, family)
            if GP_types(A1, B1, p, fam):
                break
        else:
            c = sample_form(rng, p, family, basis)
            f = form_pari(c)
            A1, B1, A2, B2 = GP_wei(f, p)
            if GP_types(A1, B1, p, fam) and GP_types(A2, B2, p, fam):
                break
    else:
        rec["error"] = "no valid sample"
        return rec
    rec["attempts"] = attempt
    rec["A1"], rec["B1"] = coeffs(A1, 8), coeffs(B1, 12)
    if not control:
        rec["f"] = c
        rec["A2"], rec["B2"] = coeffs(A2, 8), coeffs(B2, 12)
        rec["fibres1"] = [int(x) for x in GP_fibres(f, A1, B1, p, 1)]
        rec["fibres2"] = [int(x) for x in GP_fibres(f, A2, B2, p, 2)]
        rec["iso"] = int(GP_iso(A1, B1, A2, B2, p))
    S1 = [int(x) for x in GP_Sall(A1, B1, p, m)]
    rec["S1"] = S1
    if not control:
        S2 = [int(GP_S(A2, B2, p, d)) for d in range(1, min(m, check2) + 1)]
        rec["S2"] = S2
        rec["S_equal"] = S2 == S1[:len(S2)]
    P = GP_pol(S1, p)
    if P == 0:
        rec["error"] = "middle coefficient zero"
        return rec
    rec["P"] = [int(pari.polcoef(P, k)) for k in range(N, -1, -1)]
    rec["weil"] = bool(GP_weil(P, p))
    # Artin-Tate: P(p) / p^(N-1) = |Br(X)| |disc NS(X)|, with NS(X) the classes defined over F_p.  The lattice of
    # the fibre types has discriminant 1 (generic, e8); for e8a2 it is 3 if the fibre I_3 is split and the rational
    # part has discriminant 2 with a factor 2 from the eigenvalue -p, so P(p) / p^(N-1) = |Br(X)|, if it is not.
    br = pari.subst(P, "t", p) / pari(p) ** (N - 1)
    if fam == 2:
        rec["split_inf"] = int(GP_tr(pari.polcoef(A1, 8), pari.polcoef(B1, 12)))
        if rec["split_inf"] == 1:
            br = br / 3
    rec["br"] = str(br)
    rec["br_square"] = bool(pari.denominator(br) == 1 and pari.issquare(br))
    rec["br_4square"] = bool(pari.denominator(br / 4) == 1 and pari.issquare(br / 4))
    rec.update(lattice_side(P, p, bnf_limit))
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def summary(recs, control):
    out = []
    ok = [r for r in recs if "P" in r]
    out.append(f"samples {len(recs)}, with a polynomial {len(ok)}")
    if not control:
        nf1 = sum(r["fibres1"][0] for r in recs if "fibres1" in r)
        of1 = sum(r["fibres1"][1] for r in recs if "fibres1" in r)
        nf2 = sum(r["fibres2"][0] for r in recs if "fibres2" in r)
        of2 = sum(r["fibres2"][1] for r in recs if "fibres2" in r)
        out.append(f"fibres of Y against fibres of X1: {of1} of {nf1} agree; against X2: {of2} of {nf2}")
        out.append(f"trace sums of X1 and X2 equal: {sum(r.get('S_equal', False) for r in ok)} of {len(ok)}")
        out.append(f"X1 and X2 isomorphic over F_p: {dict(Counter(r['iso'] for r in ok))}  (0 no, 1 yes, 2 twist only)")
        nz = [r for r in ok if r["br"] != "0"]
        out.append(f"order of Br(X1) by Artin-Tate is 4 times a square: {sum(r['br_4square'] for r in nz)} of {len(nz)} with P(p) not 0")
    out.append(f"Weil bound holds: {sum(r['weil'] for r in ok)} of {len(ok)}")
    nz = [r for r in ok if r["br"] != "0"]
    out.append(f"order of Br by Artin-Tate is a square: {sum(r['br_square'] for r in nz)} of {len(nz)} with P(p) not 0")
    out.append(f"  and is even: {sum(r['br_square'] and int(r['br']) % 2 == 0 for r in nz)}")
    out.append(f"status: {dict(Counter(r['status'] for r in ok))}")
    irr = [r for r in ok if r["status"] not in ("reducible", "not ordinary", "uncertified")]
    out.append(f"irreducible, ordinary, certified: {len(irr)}")
    out.append(f"  multiplicity of t - 1 in P mod 2: {dict(sorted(Counter(r['m2'] for r in irr).items()))}")
    out.append(f"  above 2: {dict(Counter(r['cat2'] for r in irr))}")
    cc = [r for r in irr if r["cond_C"]]
    out.append(f"condition C: {len(cc)}")
    done = [r for r in cc if "c_trivial" in r]
    out.append(f"  class groups computed: {len(done)}")
    if done:
        out.append(f"  direct test and class in G agree: {sum((r['ord_G'] == 1) == r['c_trivial'] for r in done)} of {len(done)}")
        gs = sorted(r["G_size"] for r in done)
        out.append(f"  |G|: min {gs[0]}, median {gs[len(gs) // 2]}, max {gs[-1]}")
        out.append(f"  order of the class of l in G: {dict(sorted(Counter(r['ord_G'] for r in done).items()))}")
        if not control:
            tab = Counter((r["c_trivial"], r["iso"] == 1) for r in done)
            out.append("  prediction (class trivial <=> isomorphic):")
            out.append(f"    trivial and isomorphic      {tab[(True, True)]}")
            out.append(f"    trivial and not isomorphic  {tab[(True, False)]}")
            out.append(f"    not trivial and isomorphic  {tab[(False, True)]}")
            out.append(f"    not trivial, not isomorphic {tab[(False, False)]}")
    if not control:
        rest = [r for r in ok if not r.get("cond_C")]
        out.append(f"outside condition C: {len(rest)}, of which isomorphic {sum(r['iso'] == 1 for r in rest)}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=sorted(FAMILIES), default="e8")
    ap.add_argument("--p", type=int, default=5)
    ap.add_argument("--count", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--check2", type=int, default=6)
    ap.add_argument("--bnf-limit", type=int, default=600)
    ap.add_argument("--out")
    ap.add_argument("--graph", help="analyse the pairs of an e8a2 result file and exit")
    ap.add_argument("--summary", help="print again the summary of a result file and exit")
    args = ap.parse_args()
    if args.summary:
        recs = [json.loads(line) for line in open(args.summary)]
        print(f"family {recs[0]['family']}, p = {recs[0]['p']}, control = {recs[0]['control']}, seed = {min(r['seed'] for r in recs) // 1000003}")
        for line in summary(recs, recs[0]["control"]):
            print(line)
        return 0
    if args.graph:
        for line in graph(args.graph):
            print(line)
        return 0
    jobs = [(args.p, args.family, args.seed * 1000003 + i, args.control, args.check2, args.bnf_limit) for i in range(args.count)]
    recs = []
    out = open(args.out, "w") if args.out else None
    with Pool(args.workers) as pool:
        for rec in pool.imap_unordered(one, jobs):
            recs.append(rec)
            if out:
                out.write(json.dumps(rec) + "\n")
                out.flush()
    print(f"family {args.family}, p = {args.p}, control = {args.control}, seed = {args.seed}")
    for line in summary(recs, args.control):
        print(line)


if __name__ == "__main__":
    sys.exit(main())
