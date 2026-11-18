#!/usr/bin/env python
"""
sixlines_prym.py -- the abelian fourfold of a six-line double plane over F_p, built from the six lines (note 14).

For six rational lines L1..L6 in general position over F_p and an ordered pair (L5, L6), Paranjape's construction
(Compositio Math. 68 (1988)) is carried out over F_p:

  1. t = L5 / L6 parametrises the lines m_t through L5 n L6.  lam(t) = N(t) / D(t) is the cross-ratio of the four
     points m_t n L1..L4.  lam5 = lam(0), lam6 = lam(oo), and b1, b2 are the two branch values of lam.
  2. E is the double cover of the lam-line branched at lam5, lam6, b1, b2, with origin o over lam5 and the point
     tau of order 2 over lam6.  In the coordinate x = 1 / (lam - lam5) it is Y^2 = W(X), a monic cubic.
  3. a, b, c are points of E over lam = 0, 1, oo.  E is taken in the quadratic twist on which Frobenius negates an
     even number of them: none, or two.
  4. c4 is the point with 2 c4 = tau whose lam-coordinate is NOT lam_c = N_1 / D_1 (the ratio of the coefficients
     of t); the other half belongs to another surface (--other-half).
  5. The branch points are B = {Q, a + rho - Q, b + rho - Q, c + rho - Q} and the pole is c4 + rho, where
     Q = O and rho = -(a + b + c) if a, b, c are rational, and otherwise Q is a point of E(F_{p^2}) with
     Q - Frob(Q) = the sum of the two that are not, and rho = Q + Frob(Q) - the rational one.
  6. F is the function on E with divisor 2 B - 4 (c4 + rho) - 4 (O), leading coefficient 1 at O, and
     C : v^4 = F is the curve of genus 5.  J : v -> i v, and the fourfold is the Prym variety of C -> E.

The traces T_k of Frobenius^k on the i-eigenspace of H^1(C) are character sums over E(F_{p^k}):

  mode A  (c4 rational, p = 1 mod 4):  T_k = - sum chi_4(F(x)).
  mode B  (c4 conjugate to -c4, p = 3 mod 4):  F Frob(F) = G^4 with G rational, the curve descends with
          Frobenius (x, v) -> (Frob x, G(Frob x) / v^p), and T_k = - sum G(x) / F(x)^((q+1)/4) for odd k,
          T_k = - sum chi_4(F(x)) for even k  (q = p^k).

In the two remaining cases Frobenius anticommutes with J; they are only recorded.

From T_1, T_2, T_3 the quartic chi = X^4 - c1 X^3 + c2 X^2 - c3 X + p^2 is formed, and its second exterior power
t^6 - c2 t^5 + (c1 c3 - p^2) t^4 - (c1^2 p^2 + c3^2 - 2 c2 p^2) t^3 + ... is compared with the polynomial of the
surface, taken from the character sums S_1, S_2, S_3 stored in sixlines_p{p}.jsonl (note 2).  "equal" means
they agree, "twist" that the second exterior power is P(-t), the polynomial of the quadratic twist.

  --check p b c e f     one configuration, all 30 ordered pairs, verbose
  --run p [p ...]       all configurations in sixlines_p{p}.jsonl; --out file for the records
  --pairs n             stop after n computed ordered pairs per configuration (default 8; 0 = all 30 are tried)
  --full-below P        try all 30 ordered pairs when p < P (default 25)
  --workers n           processes (default 8)
  --other-half          use the other half of tau in step 4 (a control: the result is then a different polynomial)
  --summary file        recompute the statistics from an output file

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cypari2 import Pari  # noqa: E402

pari = Pari()
pari.allocatemem(4 * 10 ** 8)
pari("'t; 'la; 'x")

pari(r"""
pr_lambda(L, idx, p) = {
  \\ [N, D]: the cross-ratio of the points of L1..L4 on the line L5 = t L6 is N(t) / D(t)
  my(l5 = L[idx[5],], l6 = L[idx[6],], B, Z, s = vector(4), co);
  for(j = 1, 3, Z = vector(3, k, k == j); B = Mod(matconcat([l5; l6; Z]), p); if(matdet(B) != 0, break));
  for(i = 1, 4, co = matsolve(B~, Mod(L[idx[i],]~, p)); s[i] = -(co[1] * 't + co[2]) / co[3]);
  [(s[1] - s[3]) * (s[2] - s[4]), (s[1] - s[4]) * (s[2] - s[3])];
}
""")

pari(r"""
pr_curve(N, D, p, other = 0) = {
  \\ [1, W, [X(a), X(b), X(c)], rat, Xc, is Yc rational, are b1 and b2 rational] or [0] (degenerate).
  \\ E : Y^2 = W(X) is the twist on which an even number of a, b, c are not rational: rat = 0 if all three are
  \\ rational, else the index of the rational one.
  my(l5, l6, lc, Dl, d0, d1, d2, C3, be, sq, c = Mod(1, p), rat = 0, m, W, e1, Xc);
  if(subst(D, 't, 0) == 0 || polcoef(D, 2) == 0 || polcoef(N, 2) == 0 || subst(N, 't, 0) == 0 || polcoef(D, 1) == 0, return([0]));
  l5 = subst(N, 't, 0) / subst(D, 't, 0); l6 = polcoef(N, 2) / polcoef(D, 2); lc = polcoef(N, 1) / polcoef(D, 1);
  if(l5 == l6 || l5 == 1 || l6 == 1 || lc == l5, return([0]));
  Dl = poldisc(N - 'la * D, 't);
  d0 = polcoef(Dl, 0, 'la); d1 = polcoef(Dl, 1, 'la); d2 = polcoef(Dl, 2, 'la);
  C3 = (1 + (l5 - l6) * 'x) * (d2 * (l5 * 'x + 1)^2 + d1 * 'x * (l5 * 'x + 1) + d0 * 'x^2);
  if(poldegree(C3) != 3 || poldisc(C3) == 0, return([0]));
  be = [-1 / l5, 1 / (1 - l5), Mod(0, p)];
  sq = vector(3, j, issquare(subst(C3, 'x, be[j])));
  if(!sq[3], error("the points over lam = oo are not rational"));
  if(!sq[1] && !sq[2], rat = 3);
  if(sq[1] + sq[2] == 1, c = Mod(znprimroot(p), p); rat = if(sq[1], 2, 1));
  m = c * pollead(C3);
  W = m^2 * c * subst(C3, 'x, 'x / m);
  e1 = m / (l6 - l5); Xc = m / (lc - l5);
  if(subst(W, 'x, e1) != 0 || (Xc - e1)^2 != subst(deriv(W, 'x), 'x, e1), error("lam_c is not a half of tau"));
  if(!other, Xc = 2 * e1 - Xc);
  [1, W, m * be, rat, Xc, issquare(subst(W, 'x, Xc)), issquare(poldisc(Dl, 'la))];
}
""")

pari(r"""
pr_line(E, A, B) = {
  \\ [u, v, w]: the line u Y + v X + w through A and B (tangent if A = B), u = 1 unless it is vertical
  my(s);
  if(A == [0], if(B == [0], return([0, 0, 1])); return([0, 1, -B[1]]));
  if(B == [0], return([0, 1, -A[1]]));
  if(A[1] == B[1],
    if(A[2] != B[2] || A[2] == 0, return([0, 1, -A[1]]));
    s = (3 * A[1]^2 + 2 * E.a2 * A[1] + E.a4) / (2 * A[2]),
    s = (B[2] - A[2]) / (B[1] - A[1]));
  [1, -s, s * A[1] - A[2]];
}
""")

pari(r"""
pr_half(E, pos, R) = {
  \\ lines [num, den] of the function with divisor (R) + sum(pos) + (-S) - n (O), S = R + sum(pos); also S
  my(num = List(), den = List(), A = R, S, n = #pos);
  if(n == 0, if(R != [0], listput(num, pr_line(E, R, ellneg(E, R)))); return([Vec(num), Vec(den), R]));
  for(j = 1, n,
    S = elladd(E, A, pos[j]);
    if(!(A == [0] && j < n),
      listput(num, pr_line(E, A, pos[j]));
      if(j < n && S != [0], listput(den, pr_line(E, S, ellneg(E, S)))));
    A = S);
  [Vec(num), Vec(den), A];
}
""")

pari(r"""
pr_chain(E, D, R = [0]) = {
  \\ D = [[P, n], ...] with sum n P = O (the coefficient of O is implied): lines of the function with divisor D
  \\ and leading coefficient 1 at O.  R shifts the intermediate points.
  my(pos = List(), neg = List(), a, b);
  for(i = 1, #D, for(j = 1, abs(D[i][2]), if(D[i][2] > 0, listput(pos, D[i][1]), listput(neg, D[i][1]))));
  a = pr_half(E, Vec(pos), R); b = pr_half(E, Vec(neg), R);
  if(a[3] != b[3], error("divisor not principal"));
  [concat(a[1], b[2]), concat(a[2], b[1])];
}
""")

pari(r"""
pr_eval(ch, X, Y) = {
  my(n = 1, d = 1);
  for(i = 1, #ch[1], n *= ch[1][i][1] * Y + ch[1][i][2] * X + ch[1][i][3]);
  for(i = 1, #ch[2], d *= ch[2][i][1] * Y + ch[2][i][2] * X + ch[2][i][3]);
  [n, d];
}
""")

pari(r"""
pr_value(E, D, X, Y, Rs) = {
  \\ value at (X, Y), a point outside the support of D, of the function with divisor D, by a shifted chain
  my(v);
  for(i = 1, #Rs, v = pr_eval(pr_chain(E, D, Rs[i]), X, Y); if(v[1] != 0 && v[2] != 0, return(v[1] / v[2])));
  error("no chain avoids the point");
}
""")

pari(r"""
pr_point(W, X, one, om, n0) = {
  \\ the point with this X in F_p: Y = r if W(X) = r^2, else Y = om r with W(X) = n0 r^2 and om^2 = n0; r is the
  \\ root with the smaller representative, so that the signs agree in all the fields used
  my(w = subst(W, 'x, X), p = w.mod, r);
  if(issquare(w), r = lift(sqrt(w)); if(2 * r > p, r = p - r); return([X * one, r * one]));
  r = lift(sqrt(w / n0)); if(2 * r > p, r = p - r);
  [X * one, om * r];
}
""")

pari(r"""
pr_findq(W, bx, rat, Xc, p) = {
  \\ two of the points a, b, c are conjugate to their negatives.  Returns [u, v, u', v'] in F_p with
  \\ Q = (u + v om, u' + v' om) in E(F_{p^2}) and Q - Frob(Q) = their sum.  The solutions differ by rational
  \\ points, which translate the curve.  The one with the least coordinates is taken among those for which no
  \\ branch point falls on the pole c4 + rho or on its conjugate, so that the choice is the same over every field
  \\ and in every run.
  my(g = ffgen(p^2, 'g), one = g^0, z = ffprimroot(g), n0 = if(p % 4 == 3, Mod(-1, p), Mod(znprimroot(p), p)), om = sqrt(n0 * one),
     a2 = polcoef(W, 2) * one, a4 = polcoef(W, 1) * one, a6 = polcoef(W, 0) * one, E = ellinit([0, a2, 0, a4, a6]),
     r = vector(3, j, pr_point(W, bx[j], one, om, n0)), c4 = pr_point(W, Xc, one, om, n0), tgt = [0], e2 = (p^2 - 1) / 2,
     X, rhs, Y, ys, co, Q, rho, B4, cs, csb, bad, key, best = 0, bestbad = 0);
  co = (t -> my(u = (t + t^p) / 2, v = (t - t^p) / (2 * om)); [lift(Mod(polcoef(u.pol, 0), p)), lift(Mod(polcoef(v.pol, 0), p))]);
  for(j = 1, 3, if(j != rat, tgt = elladd(E, tgt, r[j])));
  for(n = 0, p^2 - 1,
    X = if(n == 0, 0 * one, z^n);
    rhs = ((X + a2) * X + a4) * X + a6;
    if(rhs == 0, ys = [0 * one], if(rhs^e2 == one, Y = sqrt(rhs); ys = [Y, -Y], next));
    foreach(ys, y,
      Q = [X, y];
      if(ellsub(E, Q, [X^p, y^p]) != tgt, next);
      rho = ellsub(E, elladd(E, Q, [X^p, y^p]), r[rat]);
      B4 = concat([Q], vector(3, j, ellsub(E, elladd(E, r[j], rho), Q)));
      cs = elladd(E, c4, rho); csb = elladd(E, ellneg(E, c4), rho);
      bad = (cs == [0] || csb == [0]);
      for(i = 1, 4, if(B4[i] == cs || B4[i] == csb || B4[i] == [0], bad = 1));
      key = concat(co(X), co(y));
      if(bad,
        if(bestbad === 0 || lex(key, bestbad) < 0, bestbad = key),
        if(best === 0 || lex(key, best) < 0, best = key))));
  if(best === 0 && bestbad === 0, error("no point Q with Q - Frob(Q) = a + b"));
  if(best === 0, bestbad, best);
}
""")

pari(r"""
pr_trace(W, bx, rat, Qc, Xc, p, k, modeB) = {
  \\ [n0, n1, n2, n3]: the number of points x of E(F_{p^k}) at which the character takes the value i^j (at the
  \\ poles of order 4 the leading coefficient is used).  T_k = -(n0 - n2) - (n1 - n3) i.  0 if degenerate.
  my(q = p^k, big = if((modeB || rat) && k % 2, q^2, q), g = ffgen(big, 'g), one = g^0, z = ffprimroot(g), zs = z^((big - 1) / (q - 1)),
     a2 = polcoef(W, 2) * one, a4 = polcoef(W, 1) * one, a6 = polcoef(W, 0) * one, E = ellinit([0, a2, 0, a4, a6]),
     n0, om, r, i4, c4, Q = [0], rho, B4, cs, csb, D, Dc, DG, ch, chG, cnt = vector(4), tab,
     exq = (q - 1) / 4, exh = (q - 1) / 2, X, rhs, Y, ys, v, vG, w, exc = List(), Rs = List(), R, hit, ext, g2, emb, E2, U, V, tgt, zeros);
  \\ p = 3 mod 4: i4 is a square root of -1 in F_{p^2} and the points that are not rational have Y in i4 F_p.
  \\ p = 1 mod 4: i4 is in F_p and those points have Y in om F_p, om a square root of the least primitive root.
  if(p % 4 == 3,
    n0 = Mod(-1, p); i4 = sqrt(-one); om = i4,
    n0 = Mod(znprimroot(p), p); i4 = n0^((p - 1) / 4) * one; om = if(rat, sqrt(n0 * one), 0));
  r = vector(3, j, pr_point(W, bx[j], one, om, n0));
  c4 = pr_point(W, Xc, one, om, n0);
  tab = [one, i4, -one, -i4];
  if(rat == 0,
    rho = ellneg(E, elladd(E, elladd(E, r[1], r[2]), r[3])),
    \\ two of the three points are conjugate to their negatives: Q in E(F_{p^2}) with Q - Frob(Q) = their sum
    Q = [Qc[1] * one + Qc[2] * om, Qc[3] * one + Qc[4] * om];
    tgt = [0]; for(j = 1, 3, if(j != rat, tgt = elladd(E, tgt, r[j])));
    if(ellsub(E, Q, [Q[1]^p, Q[2]^p]) != tgt, error("Q - Frob(Q) is not a + b"));
    rho = ellsub(E, elladd(E, Q, [Q[1]^p, Q[2]^p]), r[rat]));
  \\ the branch points Q, a + rho - Q, b + rho - Q, c + rho - Q and the pole c4 + rho
  B4 = concat([Q], vector(3, j, ellsub(E, elladd(E, r[j], rho), Q)));
  cs = elladd(E, c4, rho); csb = elladd(E, ellneg(E, c4), rho);
  for(i = 1, 4, if(B4[i] == cs || B4[i] == csb, return(0)); for(j = i + 1, 4, if(B4[i] == B4[j], return(0))));
  if(cs == [0] || csb == [0], return(0));
  zeros = select(P -> P != [0], B4);
  D = concat(vector(#zeros, j, [zeros[j], 2]), [[cs, -4]]);
  DG = concat(vector(#zeros, j, [zeros[j], 1]), [[cs, -1], [csb, -1]]);
  ch = pr_chain(E, D); if(modeB && k % 2, chG = pr_chain(E, DG));
  \\ O is a pole of order 4 with leading coefficient 1 unless it is a branch point
  if(#zeros == 4, cnt[1]++);
  for(n = 0, q - 1,
    X = if(n == 0, 0 * one, zs^n);
    rhs = ((X + a2) * X + a4) * X + a6;
    if(rhs == 0, ys = [0 * one], if(rhs^exh == one, Y = sqrt(rhs); ys = [Y, -Y], next));
    foreach(ys, y,
      v = pr_eval(ch, X, y);
      if(v[1] == 0 || v[2] == 0, listput(exc, [X, y]); next);
      if(modeB && k % 2,
        vG = pr_eval(chG, X, y);
        if(vG[1] == 0 || vG[2] == 0, listput(exc, [X, y]); next);
        w = (v[1] / v[2])^((q + 1) / 4) * vG[2] / vG[1],
        w = if(big == q, (v[1] * v[2]^3)^exq, (v[1] / v[2])^exq));
      hit = 0; for(j = 1, 4, if(w == tab[j], cnt[j]++; hit = 1; break));
      if(!hit, error("the character value is not a fourth root of unity"))));
  \\ points on the auxiliary lines, and the pole cs: the chains are shifted by points over a larger field
  if(#exc,
    ext = 2; while(big^ext < 10^5 || ext % 2 == 0, ext++);
    g2 = ffgen(big^ext, 'h); emb = ffembed(g, g2);
    E2 = ellinit([0, ffmap(emb, a2), 0, ffmap(emb, a4), ffmap(emb, a6)]);
    for(i = 1, 12, R = random(E2); if(R != [0], listput(Rs, R)));
    D = ffmap(emb, D); DG = ffmap(emb, DG); tab = ffmap(emb, tab));
  foreach(Vec(exc), pt,
    if(pt == B4[1] || pt == B4[2] || pt == B4[3] || pt == B4[4], next);
    X = ffmap(emb, pt[1]); Y = ffmap(emb, pt[2]);
    if(pt == cs,
      \\ the class of the leading coefficient modulo fourth powers is the value at cs of F M^4, M the line through
      \\ cs and a point U of E(F_q): divisor D + 4 (cs) + 4 (U) + 4 (-cs - U) - 12 (O)
      until(U != cs && V != cs && V != [0],
        until(rhs != 0 && rhs^exh == one, U = zs^random(q - 1); rhs = ((U + a2) * U + a4) * U + a6);
        U = [U, sqrt(rhs)]; V = ellneg(E, elladd(E, cs, U)));
      Dc = ffmap(emb, concat(vector(#zeros, j, [zeros[j], 2]), [[U, 4], [V, 4]]));
      w = pr_value(E2, Dc, X, Y, Vec(Rs))^exq,
      if(modeB && k % 2,
        w = pr_value(E2, D, X, Y, Vec(Rs))^((q + 1) / 4) / pr_value(E2, DG, X, Y, Vec(Rs)),
        w = pr_value(E2, D, X, Y, Vec(Rs))^exq));
    hit = 0; for(j = 1, 4, if(w == tab[j], cnt[j]++; hit = 1; break));
    if(!hit, error("the character value at an exceptional point is not a fourth root of unity")));
  \\ for the twisted Frobenius, J^m Frobenius fixes the points with psi = i^m, not i^(-m): the eigenspace sees 1 / psi
  if(modeB && k % 2, cnt = [cnt[1], cnt[4], cnt[3], cnt[2]]);
  cnt;
}
""")

GP = {n: pari(n) for n in ["pr_lambda", "pr_curve", "pr_findq", "pr_trace"]}
HERE = os.path.dirname(os.path.abspath(__file__))
OTHER_HALF = 0


def forms(cfg):
    b, c, e, f = cfg
    return [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, b, c), (1, e, f)]


def gauss_mul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gauss_lin(*terms):
    """sum of coefficient * Gaussian integer"""
    return (sum(k * u[0] for k, u in terms), sum(k * u[1] for k, u in terms))


def quartic(T, p):
    """(c1, c2, c3) from the traces T_1, T_2, T_3, as exact Gaussian integers; None if a division is not exact."""
    t1, t2, t3 = T
    two_c2 = gauss_lin((1, gauss_mul(t1, t1)), (-1, t2))
    six_c3 = gauss_lin((1, gauss_mul(t1, gauss_mul(t1, t1))), (-3, gauss_mul(t1, t2)), (2, t3))
    if any(v % 2 for v in two_c2) or any(v % 6 for v in six_c3):
        return None
    return t1, (two_c2[0] // 2, two_c2[1] // 2), (six_c3[0] // 6, six_c3[1] // 6)


def wedge(c, p):
    """(e1, e2, e3) of the second exterior power of X^4 - c1 X^3 + c2 X^2 - c3 X + p^2, as Gaussian integers."""
    c1, c2, c3 = c
    c4 = (p * p, 0)
    e2 = gauss_lin((1, gauss_mul(c1, c3)), (-1, c4))
    e3 = gauss_lin((p * p, gauss_mul(c1, c1)), (1, gauss_mul(c3, c3)), (-2 * p * p, c2))
    return c2, e2, e3


def surface_e(S):
    s1, s2, s3 = S
    return s1, (s1 * s1 - s2) // 2, (s1 ** 3 - 3 * s1 * s2 + 2 * s3) // 6


def one_pair(p, cfg, pair, kmax=3):
    """Status and data for the ordered pair (L5, L6) = pair of the configuration cfg."""
    rest = [i for i in range(1, 7) if i not in pair]
    idx = rest + list(pair)
    L = pari.matrix(6, 3, [x for row in forms(cfg) for x in row])
    N, D = GP["pr_lambda"](L, idx, p)
    cur = GP["pr_curve"](N, D, p, OTHER_HALF)
    st = int(cur[0])
    if st == 0:
        return {"status": "degenerate"}
    W, bx, rat, Xc, yrat, brat = cur[1], cur[2], int(cur[3]), cur[4], bool(int(cur[5])), bool(int(cur[6]))
    out = {"c4_rational": yrat, "b_rational": brat, "rational_point": rat}
    if yrat != (p % 4 == 1):
        out["status"] = "anticommuting"
        return out
    modeB = 0 if yrat else 1
    Qc = GP["pr_findq"](W, bx, rat, Xc, p) if rat else pari("[0, 0, 0, 0]")
    T = []
    for k in range(1, kmax + 1):
        cnt = GP["pr_trace"](W, bx, rat, Qc, Xc, p, k, modeB)
        if cnt == 0:
            return {"status": "degenerate"}
        n = [int(v) for v in cnt]
        T.append((-(n[0] - n[2]), -(n[1] - n[3])))
    out["T"] = T
    c = quartic(T, p)
    if c is None:
        out["status"] = "not integral"
        return out
    out["c"] = c
    out["wedge"] = wedge(c, p)
    out["c3_is_p_conj_c1"] = c[2] == (p * c[0][0], -p * c[0][1])
    out["status"] = "computed"
    return out


def compare(res, e):
    w = res["wedge"]
    if any(v[1] for v in w):
        return "no"
    got = tuple(v[0] for v in w)
    if got == tuple(e):
        return "equal"
    if got == (-e[0], e[1], -e[2]):
        return "twist"
    return "no"


def sign_of(p, e):
    """The signs of the functional equation that are possible for the surface (note 2, item 2)."""
    import sixlines_cm as cm
    return cm.choose_sign(p, e)


def check(p, cfg):
    rows = [json.loads(line) for line in open(os.path.join(HERE, f"sixlines_p{p}.jsonl"))]
    row = [r for r in rows if r["cfg"] == cfg][0]
    e = surface_e(row["S"])
    print(f"p = {p}, configuration {cfg}: surface (e1, e2, e3) = {e}, signs {sign_of(p, e)}")
    for pair in itertools.permutations(range(1, 7), 2):
        res = one_pair(p, cfg, pair)
        line = f"  pair {pair}: {res['status']}"
        if res["status"] == "computed":
            line += f"  c = {res['c']}  wedge = {res['wedge']}  {compare(res, e)}"
        elif "c4_rational" in res:
            line += f"  (c4 rational: {res['c4_rational']}, b1 b2 rational: {res['b_rational']})"
        print(line, flush=True)


def config_record(job):
    p, r, max_pairs = job
    e = surface_e(r["S"])
    rec = {"p": p, "cfg": r["cfg"], "e": list(e), "signs": sign_of(p, e), "pairs": []}
    if rec["signs"] == [1]:
        import sixlines_cm as cm
        P = cm.pol(cm.charpoly(p, e))
        rec["degrees"] = sorted(int(pari.poldegree(f)) for f, m in zip(*pari.factor(P)) for _ in range(int(m)))
        rec["ordinary"] = e[0] % p != 0
    seen = Counter()
    used = 0
    for pair in itertools.permutations(range(1, 7), 2):
        res = one_pair(p, r["cfg"], pair)
        seen[res["status"]] += 1
        if res["status"] == "computed":
            rec["pairs"].append({"pair": list(pair), "T": res["T"], "c": res["c"], "compare": compare(res, e),
                                 "c3_is_p_conj_c1": res["c3_is_p_conj_c1"], "c4_rational": res["c4_rational"],
                                 "rational_point": res["rational_point"]})
            used += 1
            if max_pairs and used >= max_pairs:
                break
        elif res["status"] == "not integral":
            rec["pairs"].append({"pair": list(pair), "T": res["T"], "compare": "not integral"})
    rec["seen"] = dict(seen)
    return rec


def run(p, fh, max_pairs, workers):
    rows = [json.loads(line) for line in open(os.path.join(HERE, f"sixlines_p{p}.jsonl"))]
    t0 = time.time()
    jobs = [(p, r, max_pairs) for r in rows]
    with Pool(workers) as pool:
        recs = pool.map(config_record, jobs, chunksize=4)
    stats = tally(recs, p, fh)
    print(f"p = {p}: {dict(sorted(stats.items()))}  ({time.time() - t0:.0f} s)", flush=True)
    return stats


def mu_nu(p, e):
    """(mu, nu) with P(p) = 64 p^5 mu^2 and P(-p) = 64 p^5 nu^2 (note 13), or None."""
    e1, e2, e3 = e
    c2, c3 = e2 - 3 * p * p, e3 - 2 * p * p * e1
    a = 8 * p ** 3 - 4 * p * p * e1 + 2 * p * c2 - c3
    b = 8 * p ** 3 + 4 * p * p * e1 + 2 * p * c2 + c3
    if a % (64 * p * p) or b % (64 * p * p) or a < 0 or b < 0:
        return None
    a, b = a // (64 * p * p), b // (64 * p * p)
    m, n = math.isqrt(a), math.isqrt(b)
    return (m, n) if m * m == a and n * n == b else None


def tally(recs, p, fh=None):
    stats = Counter()
    for rec in recs:
        signs, seen = rec["signs"], rec["seen"]
        stats["configurations"] += 1
        kind = "sign +" if signs == [1] else ("sign -" if signs == [-1] else f"signs {signs}")
        if signs == [1]:
            kind += ", P irreducible" if rec["degrees"] == [6] else ", P reducible"
        stats[kind] += 1
        comps = [x["compare"] for x in rec["pairs"]]
        if comps and signs != [1]:
            # the second exterior power has sign +, so it cannot be the polynomial of this surface
            stats[f"{kind}: has pairs of commuting type (not comparable)"] += 1
        elif comps:
            stats[f"{kind}: computed"] += 1
            if all(x in ("equal", "twist") for x in comps):
                stats[f"{kind}: every pair gives P or its twist"] += 1
            else:
                stats[f"{kind}: MISMATCH"] += 1
                print(f"  mismatch: p = {p} cfg = {rec['cfg']} e = {rec['e']} {rec['pairs']}", flush=True)
        elif seen.get("anticommuting"):
            stats[f"{kind}: only anticommuting pairs"] += 1
        else:
            stats[f"{kind}: no usable pair"] += 1
        if seen.get("anticommuting") and seen.get("computed"):
            stats[f"{kind}: both computed and anticommuting pairs"] += 1
        good = [x for x in rec["pairs"] if x["compare"] in ("equal", "twist")]
        rational = [x for x in good if x["rational_point"] == 0]
        if signs == [1]:
            stats["pairs computed (sign +)"] += len(comps)
            stats["pairs equal"] += sum(x == "equal" for x in comps)
            stats["pairs twist"] += sum(x == "twist" for x in comps)
        else:
            stats["pairs computed (other signs)"] += len(comps)
        stats["pairs with c3 = p conj(c1)"] += sum(x["c3_is_p_conj_c1"] for x in good)
        stats["pairs with c1 = 0 mod 4"] += sum(x["c"][0][0] % 4 == 0 and x["c"][0][1] % 4 == 0 for x in good)
        stats["pairs with a, b, c rational"] += len(rational)
        stats["pairs with a, b, c rational and c2 = 6 mod 8"] += sum(x["c"][1][0] % 8 == 6 for x in rational)
        if rational:
            stats["configurations with a rational pair"] += 1
            stats["configurations whose rational pairs agree on the twist"] += len({x["compare"] for x in rational}) == 1
        for key in ("anticommuting", "degenerate"):
            stats[f"pairs {key}"] += seen.get(key, 0)
        # the quartic of note 13: c1 = 4 (nu + mu i) up to a unit and conjugation
        if signs == [1] and good:
            mn = mu_nu(p, rec["e"])
            stats["pairs with |Re c1|, |Im c1| = 4 mu, 4 nu in some order"] += sum(
                mn is not None and sorted(abs(v) for v in x["c"][0]) == sorted(4 * v for v in mn) for x in good)
        if fh:
            fh.write(json.dumps(rec) + "\n")
    if fh:
        fh.flush()
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", nargs=5, type=int, metavar=("p", "b", "c", "e", "f"))
    ap.add_argument("--run", nargs="*", type=int, default=[])
    ap.add_argument("--pairs", type=int, default=8)
    ap.add_argument("--full-below", type=int, default=25)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out")
    ap.add_argument("--other-half", action="store_true")
    ap.add_argument("--summary")
    args = ap.parse_args()
    global OTHER_HALF
    OTHER_HALF = 1 if args.other_half else 0
    if args.summary:
        recs = [json.loads(line) for line in open(args.summary)]
        total = Counter()
        for p in sorted({r["p"] for r in recs}):
            st = tally([r for r in recs if r["p"] == p], p)
            print(f"p = {p}: {dict(sorted(st.items()))}")
            total += st
        print(f"total: {dict(sorted(total.items()))}")
    if args.check:
        check(args.check[0], args.check[1:])
    fh = open(args.out, "w") if args.out else None
    total = Counter()
    for p in args.run:
        total += run(p, fh, 0 if p < args.full_below else args.pairs, args.workers)
    if args.run:
        print(f"total: {dict(sorted(total.items()))}", flush=True)


if __name__ == "__main__":
    main()
