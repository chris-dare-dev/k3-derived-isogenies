#!/usr/bin/env python
"""
sixlines_cm.py -- Frobenius polynomials of six-line double planes and the lattice count they
should match (note 2).

Input: the rows written by `sixlines.py enumerate`.  For each class of configurations the
characteristic polynomial of Frobenius on the rank-6 part is

    Pc(t) = t^6 - e1 t^5 + e2 t^4 - e3 t^3 + p^2 e2 t^2 - p^4 e1 t + p^6,

with e1, e2, e3 the elementary symmetric functions obtained from S_1, S_2, S_3.  The quadratic
twist w^2 = c f replaces Pc(t) by Pc(-t); polynomials are grouped up to that sign.

For each ordinary, Q-irreducible Pc the script computes with PARI:

  E = Q[t]/Pc, its class group, units, roots of unity;
  the index in O_E of the smallest order R_min that the endomorphism ring of the transcendental
      lattice must contain: Z[F, p^2/F, (F-1)/2] together with the p-adic unit F/p on the slope-1 part;
  N[tau], the number of classes of unimodular hermitian O_E-lattices (I, alpha) of rank one
      (alpha I Ibar D_E = O_E, alpha in E_0) that are positive at the real place tau of E_0 and
      negative at the other two.

Prediction (note 2): if R_min = O_E, the number of K3 surfaces over F_p with polynomial Pc, which is
(self-associated classes) + (other classes)/2 over the classes with polynomial Pc(t) or Pc(-t),
equals N[tau], and N[tau] does not depend on tau.  If R_min is smaller, lattices over the
intermediate orders also contribute, so N[tau] is only a lower bound.

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import sys
from collections import defaultdict
from fractions import Fraction
from math import comb

import numpy as np
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9, silent=True)

pari(r"""
sl_gens(nf, p) = {
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = p^2 / F, pp = factor(p)[1, 1], M = 4,
     n = (pp - 1) * pp^(M - 1), e0, e1);
  e0 = lift(lift(Mod(Mod(1, pp^M) * t, P)^n));
  e1 = 1 - Mod(e0, P) - subst(e0, t, V);
  [(F - 1) / 2, (V - 1) / 2, F, V, e1 * F / p, e1 * V / p];
}
""")
pari(r"""
sl_integral(nf, g) = denominator(nfalgtobasis(nf, g)) == 1;
""")
pari(r"""
sl_index(nf, gens) = {
  my(n = poldegree(nf.pol), B, C, cols);
  B = Mat(nfalgtobasis(nf, 1));
  while(1,
    cols = B;
    for(i = 1, #B, for(j = 1, #gens, cols = concat(cols, nfeltmul(nf, gens[j], B[,i]))));
    C = mathnf(cols);
    if(C == B, break);
    B = C);
  if(#B < n, 0, abs(matdet(B)));
}
""")
pari(r"""
sl_mask(a, t, rts) = {
  \\ Signs of a real element at the three real places, as a bit mask of the negative ones.  The places
  \\ are ordered by the real part of the root of P, as in rts.  The value at a place can be many orders
  \\ of magnitude below the coefficients, so the precision is raised until every sign is certain.
  my(P, f, d, prec = 60, r, v, R, S, done = 0);
  if(type(a) != "t_POLMOD", return(if(a > 0, 0, 7)));
  P = a.mod; f = lift(a); d = poldegree(P);
  while(!done,
    prec *= 2;
    localprec(prec);
    r = select(z -> imag(z) > 0, Vec(polroots(P)));
    r = vecsort(r, (x, y) -> sign(real(x) - real(y)));
    v = vector(#r, j, real(subst(f, t, r[j])));
    R = vecmax(apply(abs, r));
    S = sum(i = 0, d - 1, abs(polcoef(f, i)) * R^i);
    done = (vecmin(apply(abs, v)) > S * 10.^(10 - prec)));
  sum(j = 1, 3, if(v[j] < 0, 2^(j - 1), 0));
}
""")
pari(r"""
sl_count(bnf, p) = {
  my(nf = bnf.nf, P = nf.pol, t = variable(P), F = Mod(t, P), V = p^2 / F, sig = lift(V),
     cyc = bnf.cyc, gen = bnf.gen, fu = bnf.fu, w = bnf.tu[1], zeta = bnf.tu[2],
     s = F + V, h, bnf0, eta, rts, Sigma, cnt = [0, 0, 0], nB, LA, Mm, cc, el, ex, idxA, idxB,
     I, J, pr, beta, alpha, found, sa, m, real_gen = 0, princ = 0, E);
  h = minpoly(s, 'x);
  bnf0 = bnfinit(h, 1);
  rts = select(r -> imag(r) > 0, Vec(polroots(P)));
  rts = vecsort(rts, (a, b) -> sign(real(a) - real(b)));
  eta = vector(2, k, subst(lift(bnf0.fu[k]), 'x, s));
  Sigma = Set(concat(vector(2, a, concat(vector(2, b, vector(2, c,
            bitxor(bitxor(if(a == 2, 7, 0), if(b == 2, sl_mask(eta[1], t, rts), 0)),
                   if(c == 2, sl_mask(eta[2], t, rts), 0))))))));
  nB = 0;
  for(a = 0, 1, for(b = 0, 1,
    m = sl_mask(eta[1]^a * eta[2]^b, t, rts);
    if(m == 0 || m == 7, nB++)));
  idxB = 4 / nB;
  Mm = Mat(vector(3, k, nfalgtobasis(nf, s^(k - 1))));
  LA = matrix(2, 2);
  for(i = 1, 2,
    el = nfbasistoalg(nf, fu[i]);
    el = el * subst(lift(el), t, V);
    cc = matinverseimage(Mm, nfalgtobasis(nf, el));
    ex = bnfisunit(bnf0, Mod(cc[1] + cc[2] * 'x + cc[3] * 'x^2, h));
    LA[i, 1] = ex[1]; LA[i, 2] = ex[2]);
  idxA = abs(matdet(LA));
  my(todo = List());
  if(#cyc == 0, listput(todo, []), forvec(e = vector(#cyc, i, [0, cyc[i] - 1]), listput(todo, e)));
  for(q = 1, #todo,
    I = if(#cyc == 0, idealhnf(nf, 1), idealfactorback(nf, gen, todo[q]));
    J = idealmul(nf, idealmul(nf, I, nfgaloisapply(nf, sig, I)), nf.diff);
    pr = bnfisprincipal(bnf, idealinv(nf, J), 3);
    if(pr[1] != 0, next);
    princ++;
    beta = nfbasistoalg(nf, pr[2]);
    found = 0;
    for(a = 0, w - 1, for(b1 = 0, 1, for(b2 = 0, 1,
      if(!found,
        alpha = beta * nfbasistoalg(nf, zeta)^a * nfbasistoalg(nf, fu[1])^b1 * nfbasistoalg(nf, fu[2])^b2;
        if(subst(lift(alpha), t, V) == alpha, found = 1)))));
    if(!found, next);
    real_gen++;
    sa = sl_mask(alpha, t, rts);
    for(tau = 1, 3,
      m = bitxor(sa, 7 - 2^(tau - 1));
      if(setsearch(Sigma, m), cnt[tau]++)));
  [bnf.no, cnt, idxA / idxB, w, bnf0.no, #Sigma, princ, real_gen, nf.disc];
}
""")

pari(r"""
sl_order(nf, gens) = {
  my(B, C, cols);
  B = Mat(nfalgtobasis(nf, 1));
  while(1,
    cols = B;
    for(i = 1, #B, for(j = 1, #gens, cols = concat(cols, nfalgtobasis(nf, nfbasistoalg(nf, B[,i]) * gens[j]))));
    C = mathnf(cols);
    if(C == B, break);
    B = C);
  B;
}
""")
pari(r"""
sl_conductor(nf, B) = {
  my(P = nf.pol, n = poldegree(P), w = nf.zk, T, D, A);
  T = matrix(n, n, i, j, nfelttrace(nf, Mod(w[i], P) * Mod(w[j], P)));
  D = (B~ * T)^(-1);
  A = idealhnf(nf, nfbasistoalg(nf, D[,1]));
  for(j = 2, n, A = idealadd(nf, A, nfbasistoalg(nf, D[,j])));
  idealinv(nf, idealmul(nf, A, nf.diff));
}
""")
pari(r"""
sl_ideals(nf, p, cond) = {
  my(P = nf.pol, t = variable(P), sig = lift(p^2 / Mod(t, P)), fa = idealfactor(nf, cond),
     k = matsize(fa)[1], L = List(), A);
  if(k == 0, return([idealhnf(nf, 1)]));
  forvec(e = vector(k, i, [0, fa[i, 2]]),
    A = idealhnf(nf, idealfactorback(nf, Vec(fa[, 1]), e));
    if(idealhnf(nf, nfgaloisapply(nf, sig, A)) == A && issquare(idealnorm(nf, A)), listput(L, A)));
  Vec(L);
}
""")
pari(r"""
sl_lattices(bnf, p, afr) = {
  my(nf = bnf.nf, P = nf.pol, t = variable(P), F = Mod(t, P), V = p^2 / F, sig = lift(V),
     cyc = bnf.cyc, gen = bnf.gen, fu = bnf.fu, w = bnf.tu[1], zeta = bnf.tu[2],
     s = F + V, h = minpoly(s, 'x), bnf0 = bnfinit(h, 1), eta, rts, I, J, pr, beta, alpha, found,
     todo = List(), res = List());
  rts = select(r -> imag(r) > 0, Vec(polroots(P)));
  rts = vecsort(rts, (a, b) -> sign(real(a) - real(b)));
  eta = vector(2, k, subst(lift(bnf0.fu[k]), 'x, s));
  if(#cyc == 0, listput(todo, []), forvec(e = vector(#cyc, i, [0, cyc[i] - 1]), listput(todo, e)));
  for(q = 1, #todo,
    I = if(#cyc == 0, idealhnf(nf, 1), idealhnf(nf, idealfactorback(nf, gen, todo[q])));
    J = idealmul(nf, idealmul(nf, idealmul(nf, I, nfgaloisapply(nf, sig, I)), nf.diff), afr);
    pr = bnfisprincipal(bnf, idealinv(nf, J), 3);
    if(pr[1] != 0, next);
    beta = nfbasistoalg(nf, pr[2]);
    found = 0;
    for(a = 0, w - 1, for(b1 = 0, 1, for(b2 = 0, 1,
      if(!found,
        alpha = beta * nfbasistoalg(nf, zeta)^a * nfbasistoalg(nf, fu[1])^b1 * nfbasistoalg(nf, fu[2])^b2;
        if(subst(lift(alpha), t, V) == alpha, found = 1)))));
    if(!found, next);
    listput(res, [I, alpha, sl_mask(alpha, t, rts)]));
  my(Mm = Mat(vector(3, k, nfalgtobasis(nf, s^(k - 1)))), LA = matrix(2, 2), el, cc, ex, span, reps = List(), cls = List(), u, m, v);
  for(i = 1, 2,
    el = nfbasistoalg(nf, fu[i]);
    el = el * subst(lift(el), t, V);
    cc = matinverseimage(Mm, nfalgtobasis(nf, el));
    ex = bnfisunit(bnf0, Mod(cc[1] + cc[2] * 'x + cc[3] * 'x^2, h));
    LA[i, 1] = ex[1] % 2; LA[i, 2] = ex[2] % 2);
  span = Set([[0, 0], [LA[1, 1], LA[1, 2]], [LA[2, 1], LA[2, 2]], [(LA[1, 1] + LA[2, 1]) % 2, (LA[1, 2] + LA[2, 2]) % 2]]);
  for(a = 0, 1, for(b = 0, 1,
    u = eta[1]^a * eta[2]^b;
    m = sl_mask(u, t, rts);
    if(m == 0 || m == 7,
      if(m == 7, u = -u);
      v = 0;
      for(q = 1, #cls, if(setsearch(span, [(a + cls[q][1]) % 2, (b + cls[q][2]) % 2]), v = 1));
      if(!v, listput(cls, [a, b]); listput(reps, u)))));
  [Vec(res), eta, [sl_mask(eta[1], t, rts), sl_mask(eta[2], t, rts)], nfbasistoalg(nf, zeta), Vec(reps)];
}
""")
pari(r"""
sl_latdata(nf, p, I, alpha, afr, gens, zeta) = {
  my(P = nf.pol, t = variable(P), V = p^2 / Mod(t, P), n = poldegree(P), w = nf.zk, cw, G, K, mm);
  cw = vector(n, j, subst(lift(Mod(w[j], P)), t, V));
  G = matrix(n, n, i, j, nfelttrace(nf, alpha * Mod(w[i], P) * cw[j]));
  K = idealhnf(nf, idealmul(nf, afr, I));
  mm = ((g) -> matconcat(vector(n, j, nfalgtobasis(nf, Mod(g, P) * Mod(w[j], P)))));
  [G, K, vector(#gens, i, mm(gens[i])), vector(n, i, mm(w[i])), mm(zeta)];
}
""")

pari(r"""
sl_hnfq(M) = { my(d = denominator(M)); mathnf(d * M) / d; }
""")
pari(r"""
sl_mulmat(nf, g) = {
  my(P = nf.pol, n = poldegree(P), w = nf.zk);
  matconcat(vector(n, j, nfalgtobasis(nf, Mod(g, P) * Mod(w[j], P))));
}
""")
pari(r"""
sl_conjel(nf, q, g) = { my(P = nf.pol, t = variable(P)); subst(lift(Mod(g, P)), t, q^2 / Mod(t, P)); }
""")
pari(r"""
sl_translate_p(nf, q, L) = {
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = q^2 / F, pp = factor(q)[1, 1], M = 4,
     n = (pp - 1) * pp^(M - 1), e0, e0m, e2m, xi);
  e0 = lift(lift(Mod(Mod(1, pp^M) * t, P)^n));
  e0m = Mod(e0, P); e2m = subst(e0, t, V);
  xi = pp * e0m + (1 - e0m - e2m) + e2m / pp;
  sl_hnfq(concat(sl_mulmat(nf, xi) * L, pp * L));
}
""")
pari(r"""
sl_translate_l(nf, q, L, l, lam) = {
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = q^2 / F, L1, A, B, D);
  L1 = sl_hnfq(concat(l * L, sl_mulmat(nf, F - lam) * L));
  if(abs(matdet(L1) / matdet(L)) != l, error("sublattice of index ", matdet(L1) / matdet(L)));
  A = L1 / l;
  B = sl_mulmat(nf, V - lam)^(-1) * L1;
  D = sl_hnfq(concat((A~)^(-1), (B~)^(-1)));
  sl_hnfq((D~)^(-1));
}
""")
pari(r"""
sl_iso(bnf, q, L1, a1, L2, a2) = {
  my(nf = bnf.nf, P = nf.pol, n = poldegree(P), w = nf.zk, wt = bnf.tu[1],
     zeta = nfbasistoalg(nf, bnf.tu[2]), f1 = nfbasistoalg(nf, bnf.fu[1]), f2 = nfbasistoalg(nf, bnf.fu[2]),
     J1, J2, pr, g0, nu, ex, M2, b, u, g1, x1, x2);
  J1 = sl_hnfq(matconcat(vector(n, j, sl_mulmat(nf, w[j]) * L1)));
  J2 = sl_hnfq(matconcat(vector(n, j, sl_mulmat(nf, w[j]) * L2)));
  pr = bnfisprincipal(bnf, idealmul(nf, J2, idealinv(nf, J1)), 3);
  if(pr[1] != 0, return(0));
  g0 = nfbasistoalg(nf, pr[2]);
  nu = (a1 / a2) / (g0 * sl_conjel(nf, q, g0));
  ex = bnfisunit(bnf, nu);
  if(#ex == 0, return(0));
  x1 = bnfisunit(bnf, f1 * sl_conjel(nf, q, f1)); x2 = bnfisunit(bnf, f2 * sl_conjel(nf, q, f2));
  M2 = [x1[1], x1[2]; x2[1], x2[2]];
  b = [ex[1], ex[2]] * M2^(-1);
  if(denominator(b) != 1, return(0));
  u = f1^b[1] * f2^b[2];
  if(u * sl_conjel(nf, q, u) != nu, return(0));
  g1 = g0 * u;
  for(a = 0, wt - 1, if(sl_hnfq(sl_mulmat(nf, g1 * zeta^a) * L1) == L2, return(1)));
  0;
}
""")
pari(r"""
sl_elkies(P, q, l) = {
  my(r = polrootsmod(P, l), out = List(), lam, mu);
  for(i = 1, #r,
    lam = lift(r[i]);
    if(lam == 0 || subst(P', variable(P), Mod(lam, l)) == 0, next);
    mu = lift(Mod(q, l)^2 / Mod(lam, l));
    if(mu == lam, next);
    if(lam < mu, listput(out, lam)));
  Vec(out);
}
""")

pari(r"""
sl_rbasis(nf, gens) = { my(B = sl_order(nf, gens)); vector(#B, i, nfbasistoalg(nf, B[,i])); }
""")
pari(r"""
sl_latcap(A, B) = { my(D = sl_hnfq(concat((A~)^(-1), (B~)^(-1)))); sl_hnfq((D~)^(-1)); }
""")
pari(r"""
sl_lagr(G, K, I, Rm, Om, Zm, w, maxstates) = {
  \\ Self-dual lattices L with K < L < J, stable under R, with O_E L = J, up to the roots of unity.
  \\ G: Gram matrix on the integral basis; I, K: Z-bases of J and of its dual K; Rm: multiplication
  \\ matrices of a Z-basis of R; Om: of the integral basis; Zm: of a generator of the roots of unity.
  \\ Everything is done in coordinates on J.  J/K splits into its l-primary parts, and a self-dual L
  \\ is a sum of one maximal isotropic R-submodule of each part.  For each prime l the isotropic
  \\ R-stable lattices H are explored from K: every strictly larger one contains H + R x for some x
  \\ of order l in H^perp / H, so the lines of that F_l-space are the only extensions to try.
  \\ Returns 0 if more than maxstates lattices H are met.
  my(n = #G, Ii = I^(-1), GJ = I~ * G * I, Kc = mathnf(Ii * K), d = abs(matdet(Kc)),
     ra = vector(#Rm, i, Ii * Rm[i] * I), oa = vector(#Om, i, Ii * Om[i] * I), za = Ii * Zm * I,
     fa, parts, nstates = 0, res);
  if(sl_hnfq(GJ^(-1)) != Kc, error("K is not the dual of J"));
  fa = factor(d);
  parts = vector(matsize(fa)[1]);
  for(k = 1, #parts,
    my(l = fa[k, 1], a = fa[k, 2], Jl, seen = Map(), queue, pos = 0, term = List(), H, Hd, T, U, Bi, r, x, M, f);
    if(a % 2, return([]));
    Jl = mathnf(concat(Kc, (d / l^a) * matid(n)));
    queue = List([Kc]); mapput(seen, Kc, 1);
    while(pos < #queue,
      pos++; H = queue[pos];
      Hd = sl_hnfq((H~ * GJ)^(-1));
      T = sl_latcap(Hd, H / l);
      U = l * H^(-1) * T;
      Bi = lift(matimage(Mod(U, l)));
      r = #Bi;
      if(r == 0,
        if(mathnf(matconcat(vector(#oa, j, oa[j] * H))) == Jl, listput(term, H));
        next);
      forvec(e = vector(r, i, [0, l - 1]),
        f = 0; for(i = 1, r, if(e[i] != 0, f = e[i]; break));
        if(f != 1, next);
        x = (H * (Bi * e~)) / l;
        M = mathnf(concat(concat(H, x), matconcat(vector(#ra, i, ra[i] * x))));
        if(denominator(M~ * GJ * M) != 1, next);
        if(!mapisdefined(seen, M),
          mapput(seen, M, 1); listput(queue, M);
          nstates++; if(nstates > maxstates, return(0)))));
    parts[k] = Vec(term));
  res = [Kc];
  for(k = 1, #parts,
    my(nr = List());
    for(i = 1, #res, for(j = 1, #parts[k], listput(nr, mathnf(concat(res[i], parts[k][j])))));
    res = Vec(nr));
  if(w > 2,
    my(reps = Map(), outl = List(), cur, best);
    for(i = 1, #res,
      best = res[i]; cur = mathnf(za * res[i]);
      while(cur != res[i], if(cmp(cur, best) < 0, best = cur); cur = mathnf(za * cur));
      if(!mapisdefined(reps, best), mapput(reps, best, 1); listput(outl, best)));
    res = Vec(outl));
  apply(L -> sl_hnfq(I * L), res);
}
""")

GP_nf = pari("(b)->b.nf")
GP_sl_lagr, GP_sl_rbasis = pari("sl_lagr"), pari("sl_rbasis")
GP_sl_gens, GP_sl_integral, GP_sl_index, GP_sl_count = (pari(n) for n in ("sl_gens", "sl_integral", "sl_index", "sl_count"))
GP_sl_ideals, GP_sl_lattices, GP_sl_latdata = (pari(n) for n in ("sl_ideals", "sl_lattices", "sl_latdata"))
GP_sl_order, GP_sl_conductor = pari("sl_order"), pari("sl_conductor")
GP_sl_hnfq, GP_sl_translate_p, GP_sl_translate_l, GP_sl_iso, GP_sl_elkies = (
    pari(n) for n in ("sl_hnfq", "sl_translate_p", "sl_translate_l", "sl_iso", "sl_elkies"))


# ---------------------------------------------------------------- lattices over non-maximal orders

def fmat(M):
    """PARI matrix -> list of rows of Fractions."""
    cols = [[Fraction(int(pari.numerator(x)), int(pari.denominator(x))) for x in col] for col in M]
    return [list(r) for r in zip(*cols)]


def mat_mul(A, B):
    return [[sum(a * b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


def mat_inv(A):
    n = len(A)
    M = [list(row) + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                M[r] = [x - M[r][c] * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def int_mat(A):
    assert all(x.denominator == 1 for row in A for x in row)
    return [[int(x) for x in row] for row in A]


LIMIT = 200000          # largest number of intermediate lattices in one search
NEIGHBOURS = False
SLOW = False            # use the original search (lagrangians), for comparison


def self_dual_lattices(nf, p, I, alpha, afr, gens, rbasis, zeta, w):
    """The unimodular R-lattices between K = afr I and I, as PARI matrices in integral-basis coordinates.
    None if the search is too large."""
    if SLOW:
        G, K, Rm, Om, Zm = GP_sl_latdata(nf, p, I, alpha, afr, gens, zeta)
        res = lagrangians(fmat(G), fmat(K), fmat(I), [fmat(m) for m in Rm], [fmat(m) for m in Om],
                          fmat(Zm), w, want_lattices=True)
        return None if res is None else [GP_sl_hnfq(pari_matrix(cols)) for cols in res]
    G, K, Rm, Om, Zm = GP_sl_latdata(nf, p, I, alpha, afr, rbasis, zeta)
    res = GP_sl_lagr(G, K, I, Rm, Om, Zm, w, LIMIT)
    return None if res.type() == "t_INT" else list(res)


def lagrangians(G, K, I, Rmats, Omats, Zmat, w, want_lattices=False):
    """The original search, kept for comparison (--slow); sl_lagr replaces it.
    Number of R-stable self-dual lattices L with K < L < J = I and O_E L = J, up to roots of unity.
    G: Gram matrix of the form on the integral basis; I, K: Z-bases (columns) of J and of K = a J;
    Rmats, Omats, Zmat: multiplication matrices of the generators of R, of the integral basis, of zeta."""
    n = len(I)
    Iinv = mat_inv(I)
    GJ = mat_mul(mat_mul([list(r) for r in zip(*I)], G), I)
    Kc = int_mat(mat_mul(Iinv, K))
    H = pari.mathnf(pari.matrix(n, n, [x for row in Kc for x in row]))
    Hc = [[int(H[j][i]) for j in range(n)] for i in range(n)]          # rows
    diag = [Hc[i][i] for i in range(n)]
    size = 1
    for d in diag:
        size *= d
    dual = mat_inv(GJ)                                                  # columns generate the dual of J
    assert abs(Fraction(1) / _det(GJ)) == size, "K is not the dual of J"
    if size > 1024:
        return None
    root = int(round(size ** 0.5))
    assert root * root == size

    def red(v):
        v = list(v)
        for i in range(n - 1, -1, -1):
            q = v[i] // diag[i]
            if q:
                for r in range(n):
                    v[r] -= q * Hc[r][i]
        return tuple(v)

    def act(A, v):
        return red([sum(A[r][c] * v[c] for c in range(n)) for r in range(n)])

    def add(u, v):
        return red([a + b for a, b in zip(u, v)])

    ract = [int_mat(mat_mul(mat_mul(Iinv, A), I)) for A in Rmats]
    oact = [int_mat(mat_mul(mat_mul(Iinv, A), I)) for A in Omats]
    zact = int_mat(mat_mul(mat_mul(Iinv, Zmat), I))
    elems = [red(v) for v in itertools.product(*[range(d) for d in diag])]
    zero = red([0] * n)

    def pairing_ok(u, v):
        return sum(u[i] * GJ[i][j] * v[j] for i in range(n) for j in range(n)).denominator == 1

    def group(gens):
        S = {zero}
        for g in gens:
            if g in S:
                continue
            layer, x = [], g
            while x != zero:
                layer.append(x)
                x = add(x, g)
            S = S | {add(s, y) for s in S for y in layer}
        return S

    def rmodule(x):
        orbit, todo = {x}, [x]
        while todo:
            y = todo.pop()
            for A in ract:
                z = act(A, y)
                if z not in orbit:
                    orbit.add(z)
                    todo.append(z)
        return frozenset(group(sorted(orbit)))

    cyc = {x: rmodule(x) for x in elems}
    iso = {x: all(pairing_ok(u, v) for u in cyc[x] for v in cyc[x]) and len(cyc[x]) <= root for x in elems}
    start = frozenset([zero])
    seen, frontier, found = {start}, [start], set()
    while frontier:
        Hs = frontier.pop()
        if len(Hs) == root:
            span = group(sorted({act(A, v) for A in oact for v in Hs}))
            if len(span) == size:
                found.add(Hs)
            continue
        for x in elems:
            if x in Hs or not iso[x]:
                continue
            C = cyc[x]
            if len(Hs) * len(C) > root * len(Hs & C):
                continue
            if not all(pairing_ok(u, v) for u in Hs for v in C):
                continue
            new = frozenset(add(u, v) for u in Hs for v in C)
            if new not in seen:
                seen.add(new)
                frontier.append(new)
    if w > 2:
        orbits = set()
        for Hs in found:
            orb, cur = [Hs], frozenset(act(zact, v) for v in Hs)
            while cur != Hs:
                orb.append(cur)
                cur = frozenset(act(zact, v) for v in cur)
            orbits.add(min(orb, key=sorted))
        found = orbits
    if want_lattices:
        # Z-generators of each lattice, as columns in integral-basis coordinates
        kcols = [list(c) for c in zip(*Hc)]
        return [[[sum(I[r][c] * v[c] for c in range(n)) for r in range(n)] for v in kcols + sorted(Hs)]
                for Hs in sorted(found, key=sorted)]
    return len(found)


def _det(A):
    n, M, d = len(A), [list(r) for r in A], Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d


def all_orders_count(bnf, p, gens, w):
    """N_total[tau]: unimodular hermitian lattices over every order between R_min and O_E.
    Returns (totals, per-ideal detail) or (None, reason)."""
    nf = GP_nf(bnf)
    totals, detail = [0, 0, 0], []
    cond = GP_sl_conductor(nf, GP_sl_order(nf, gens))
    rbasis = GP_sl_rbasis(nf, gens)
    for afr in GP_sl_ideals(nf, p, cond):
        lats, eta, emask, zeta, ureps = GP_sl_lattices(bnf, p, afr)
        emask = [int(m) for m in emask]
        norm = int(pari.idealnorm(nf, afr))
        per_tau = [0, 0, 0]
        for I, alpha, sa in lats:
            for tau in range(3):
                target = 7 - (1 << tau)
                unit = None
                for a in range(2):
                    for b in range(2):
                        for c in range(2):
                            if unit is None and int(sa) ^ (7 * a) ^ (emask[0] * b) ^ (emask[1] * c) == target:
                                unit = (-1) ** a * eta[0] ** b * eta[1] ** c
                if unit is None:
                    continue
                for eps in ureps:
                    res = self_dual_lattices(nf, p, I, alpha * unit * eps, afr, gens, rbasis, zeta, w)
                    if res is None:
                        return None, f"discriminant group of order {norm} too large"
                    per_tau[tau] += len(res)
        detail.append({"norm": norm, "count": per_tau})
        totals = [t + c for t, c in zip(totals, per_tau)]
    return totals, detail


def pari_matrix(cols):
    """Columns of Fractions -> PARI matrix."""
    n, m = len(cols[0]), len(cols)
    return pari("[" + ";".join(",".join(f"{cols[c][r].numerator}/{cols[c][r].denominator}" for c in range(m))
                               for r in range(n)) + "]")


def lattice_classes(bnf, q, gens, w, tau=0):
    """All classes (L, alpha) over the orders containing R_min, positive at the real place tau.
    L is a Z-basis in integral-basis coordinates (hermite normal form)."""
    nf = GP_nf(bnf)
    out = []
    cond = GP_sl_conductor(nf, GP_sl_order(nf, gens))
    rbasis = GP_sl_rbasis(nf, gens)
    for afr in GP_sl_ideals(nf, q, cond):
        lats, eta, emask, zeta, ureps = GP_sl_lattices(bnf, q, afr)
        emask = [int(m) for m in emask]
        target = 7 - (1 << tau)
        for I, alpha, sa in lats:
            unit = None
            for a in range(2):
                for b in range(2):
                    for c in range(2):
                        if unit is None and int(sa) ^ (7 * a) ^ (emask[0] * b) ^ (emask[1] * c) == target:
                            unit = (-1) ** a * eta[0] ** b * eta[1] ** c
            if unit is None:
                continue
            for eps in ureps:
                al = alpha * unit * eps
                res = self_dual_lattices(nf, q, I, al, afr, gens, rbasis, zeta, w)
                if res is None:
                    return None
                out.extend((L, al) for L in res)
    return out


def permutation(bnf, q, classes, translate):
    """The permutation of the lattice classes induced by a translation, as a list of images."""
    perm = []
    for L, al in classes:
        L2 = translate(L)
        hits = [j for j, (Lj, aj) in enumerate(classes) if GP_sl_iso(bnf, q, L2, al, Lj, aj)]
        assert len(hits) == 1, f"translate matches {len(hits)} classes"
        perm.append(hits[0])
    return perm


def cycle_type(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        n, j = 0, i
        while j not in seen:
            seen.add(j)
            j = perm[j]
            n += 1
        out.append(n)
    return sorted(out)


def neighbour_data(bnf, q, gens, w, index, primes=(3, 5, 7, 11, 13, 17, 19)):
    """Cycle types of translation by p / pbar and by l / lbar for the Elkies primes l."""
    nf = GP_nf(bnf)
    classes = lattice_classes(bnf, q, gens, w)
    if classes is None:
        return None
    P = pari("(b)->b.pol")(bnf)
    pp = int(pari.factor(q)[0][0])
    data = {"classes": len(classes),
            "p": cycle_type(permutation(bnf, q, classes, lambda L: GP_sl_translate_p(nf, q, L)))}
    perms = []
    for l in primes:
        if l == pp or index % l == 0:
            continue
        for lam in GP_sl_elkies(P, q, l):
            perm = permutation(bnf, q, classes, lambda L: GP_sl_translate_l(nf, q, L, l, lam))
            perms.append(perm)
            data.setdefault("l", []).append([l, int(lam), cycle_type(perm)])
    # orbits of the group generated by all the l-steps found
    comp = list(range(len(classes)))

    def find(x):
        while comp[x] != x:
            comp[x] = comp[comp[x]]
            x = comp[x]
        return x
    for perm in perms:
        for i, j in enumerate(perm):
            comp[find(i)] = find(j)
    data["components"] = len({find(i) for i in range(len(classes))})
    return data


def elementary(S):
    s1, s2, s3 = S
    e2n, e3n = s1 * s1 - s2, s1 ** 3 - 3 * s1 * s2 + 2 * s3
    return s1, e2n // 2, e3n // 6, (e2n % 2 == 0 and e3n % 6 == 0)


def charpoly(p, e, eps=1):
    """Monic characteristic polynomial for functional-equation sign eps (descending powers of t)."""
    e1, e2, e3 = e
    return [1, -e1, e2, -e3, eps * p * p * e2, -eps * p ** 4 * e1, eps * p ** 6]


def weil_ok(p, coeffs):
    """All roots have absolute value p (checked factor by factor, so repeated roots are harmless)."""
    for f in pari.factor(pol(coeffs))[0]:
        if not all(abs(abs(complex(z)) - p) < 1e-9 * p for z in pari.polroots(f)):
            return False
    return True


def choose_sign(p, e):
    """Signs eps for which the polynomial is a Weil polynomial with (F - 1)/2 integral.
    eps = -1 forces e3 = 0."""
    return [eps for eps in (1, -1) if (eps == 1 or e[2] == 0)
            and weil_ok(p, charpoly(p, e, eps)) and two_adic_ok(charpoly(p, e, eps))]


def two_adic_ok(coeffs):
    """Is Pc(1 + 2 s) divisible by 64, i.e. is (F - 1)/2 integral on the rank-6 part?"""
    n = len(coeffs) - 1
    shifted = [0] * (n + 1)
    for i, c in enumerate(coeffs):
        d = n - i
        for k in range(d + 1):
            shifted[k] += c * comb(d, k) * 2 ** k
    return all(v % 64 == 0 for v in shifted)


def pol(coeffs):
    n = len(coeffs) - 1
    return pari("+".join(f"({c})*t^{n - i}" for i, c in enumerate(coeffs)))


def analyse(path, do_cm=True):
    rows = [json.loads(line) for line in open(path)]
    p = rows[0]["p"]                       # the size q of the field
    prime = rows[0].get("prime", p)
    by_cfg = {tuple(r["cfg"]): r for r in rows}
    stats = defaultdict(int)
    groups = defaultdict(lambda: {"plus": [], "minus": []})
    for r in rows:
        e1, e2, e3, integral = elementary(r["S"])
        r["e"] = [e1, e2, e3]
        r["self_assoc"] = r["assoc"] == r["cfg"]
        stats["classes"] += 1
        stats["self_assoc"] += r["self_assoc"]
        if not integral:
            stats["non_integral_e"] += 1
        signs = choose_sign(p, (e1, e2, e3))
        if len(signs) != 1:
            stats["no_valid_sign" if not signs else "ambiguous_sign"] += 1
            continue
        eps = signs[0]
        co = charpoly(p, (e1, e2, e3), eps)
        if eps == -1:
            stats["sign_minus"] += 1
        a = by_cfg[tuple(r["assoc"])]
        if [abs(a["S"][0]), a["S"][1], abs(a["S"][2])] != [abs(r["S"][0]), r["S"][1], abs(r["S"][2])] \
                or (a["S"][0] * r["S"][2] != r["S"][0] * a["S"][2]):
            stats["assoc_zeta_mismatch"] += 1
        r["ordinary"] = e1 % prime != 0
        r["irreducible"] = bool(pari.polisirreducible(pol(co)))
        r["degrees"] = [int(pari.poldegree(f)) for f in pari.factor(pol(co))[0]]
        if not r["ordinary"]:
            stats["not_ordinary"] += 1
        elif not r["irreducible"]:
            stats["ordinary_reducible"] += 1
            if not r["self_assoc"]:
                stats["reducible_not_self_assoc"] += 1
        else:
            stats["ordinary_irreducible"] += 1
            if r["self_assoc"]:
                stats["irreducible_self_assoc"] += 1
            sign = 1 if e1 > 0 else -1
            groups[(sign * e1, e2, sign * e3)]["plus" if sign == 1 else "minus"].append(r)
    out = []
    for key in sorted(groups):
        g = groups[key]
        rec = {"p": p, "e": list(key), "a_plus": len(g["plus"]), "a_minus": len(g["minus"]),
               "n_self_assoc": sum(r["self_assoc"] for r in g["plus"] + g["minus"]),
               "configs": [r["cfg"] for r in g["plus"] + g["minus"]],
               "orbits": [r["orbit"] for r in g["plus"] + g["minus"]]}
        members = g["plus"] + g["minus"]
        if "frob" in members[0]:
            # K3 surfaces = association classes; Frobenius conjugation acts on them
            k3s = sorted({tuple(sorted([tuple(r["cfg"]), tuple(r["assoc"])])) for r in members})
            where = {c: i for i, pair in enumerate(k3s) for c in pair}
            rec["frob_cycles"] = cycle_type([where[tuple(by_cfg[pair[0]]["frob"])] for pair in k3s])
        n_all = rec["a_plus"] + rec["a_minus"]
        rec["k3"] = rec["n_self_assoc"] + (n_all - rec["n_self_assoc"]) / 2
        if do_cm:
            P = pol(charpoly(p, key))
            bnf = pari.bnfinit(P, 1)
            gens = GP_sl_gens(GP_nf(bnf), p)
            rec["gens_integral"] = all(bool(GP_sl_integral(GP_nf(bnf), x)) for x in gens)
            rec["index"] = int(GP_sl_index(GP_nf(bnf), gens)) if rec["gens_integral"] else None
            h, cnt, unit_index, w, h0, nsig, princ, real_gen, disc = GP_sl_count(bnf, p)
            rec.update({"h": int(h), "h0": int(h0), "w": int(w), "unit_index": int(unit_index),
                        "sign_group": int(nsig), "principal_classes": int(princ),
                        "real_generator_classes": int(real_gen), "disc": int(disc),
                        "N": [int(c) * int(unit_index) for c in cnt]})
            rec["N_total"], rec["N_detail"] = None, None
            if rec["gens_integral"]:
                rec["N_total"], rec["N_detail"] = all_orders_count(bnf, p, gens, int(w))
                if NEIGHBOURS and rec["N_total"] is not None:
                    rec["neighbours"] = neighbour_data(bnf, p, gens, int(w), rec["index"])
        out.append(rec)
    return p, stats, out


def verdict(rec):
    """Compare the number of K3 surfaces, k3 = (self-associated classes) + (other classes)/2, with the
    number of lattices: N (maximal order) and N_total (all orders containing R_min)."""
    N, T = rec["N"], rec["N_total"]
    if len(set(N)) != 1:
        return "N depends on tau"
    if T is None:
        return "not computed: " + str(rec["N_detail"]) if rec["index"] != 1 else \
            ("match" if rec["k3"] == N[0] else "MISMATCH")
    if len(set(T)) != 1:
        return "N_total depends on tau"
    if rec["index"] == 1:
        assert T == N, (T, N)
        return "match" if rec["k3"] == N[0] else "MISMATCH"
    return "match (non-maximal)" if rec["k3"] == T[0] else "MISMATCH (non-maximal)"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", help="write one JSON row per polynomial")
    ap.add_argument("--no-cm", action="store_true", help="classification only")
    ap.add_argument("--list", action="store_true", help="print every polynomial")
    ap.add_argument("--neighbours", action="store_true",
                    help="also compute the permutations of the lattice classes by p/pbar and l/lbar")
    ap.add_argument("--limit", type=int, default=200000,
                    help="largest number of intermediate lattices in one search for self-dual lattices")
    ap.add_argument("--slow", action="store_true",
                    help="use the original subgroup search (discriminant groups up to 1024)")
    args = ap.parse_args()
    global LIMIT, NEIGHBOURS, SLOW
    LIMIT = args.limit
    SLOW = args.slow
    NEIGHBOURS = args.neighbours
    fh = open(args.out, "w") if args.out else None
    for path in args.files:
        p, stats, out = analyse(path, do_cm=not args.no_cm)
        print(f"p={p}: " + ", ".join(f"{k}={v}" for k, v in sorted(stats.items())))
        if args.no_cm:
            continue
        tally = defaultdict(int)
        for rec in out:
            rec["verdict"] = verdict(rec)
            tally[rec["verdict"]] += 1
            if fh:
                fh.write(json.dumps(rec) + "\n")
            if args.list or not rec["verdict"].startswith("match"):
                print(f"  e={rec['e']} a=({rec['a_plus']},{rec['a_minus']}) self_assoc={rec['n_self_assoc']} "
                      f"k3={rec['k3']:g} N={rec['N']} N_total={rec['N_total']} h={rec['h']} "
                      f"h0={rec['h0']} unit_index={rec['unit_index']} index={rec['index']} w={rec['w']} "
                      f"-> {rec['verdict']}")
                if rec.get("neighbours"):
                    nb = rec["neighbours"]
                    print(f"      lattice side: p-step {nb['p']}; components under l-steps {nb['components']}; "
                          + "; ".join(f"l={l} lam={lam}: {ct}" for l, lam, ct in nb.get("l", [])))
                if "frob_cycles" in rec:
                    print(f"      surface side: Frobenius conjugation {rec['frob_cycles']}")
        print(f"p={p}: {len(out)} polynomials: " + ", ".join(f"{k}: {v}" for k, v in sorted(tally.items())))
    if fh:
        fh.close()


if __name__ == "__main__":
    sys.exit(main())
