#!/usr/bin/env python
"""
sixlines_jacobian.py -- arithmetic in the Jacobian of the genus-5 curve of note 14, over F_p and F_{p^r}, and the
two tests of note 23 that need it: levels, and the hermitian form modulo 4 (the character of d = -1 of note 19).

The curve.  For six rational lines and an ordered pair of them, sixlines_prym.py gives an elliptic curve
E : Y^2 = W(X) over F_p and a function F on E with divisor 2 B - 4 (c'') - 4 (O), B four points.  Here

    C_tw : v^4 = n0^tw F      (n0 the least primitive root of p, tw = 0..3; p = 1 mod 4, c'' rational)

is the curve of genus 5 and its quartic twists, J : v -> zeta v with zeta = n0^((p-1)/4), and the fourfold is the
part of Jac(C) on which J^2 = -1.  Only the case where Frobenius commutes with J over F_p is handled (p = 1 mod 4).

The method is Khuri-Makdisi's (Math. Comp. 73 (2004), "Linear algebra algorithms for divisors on an algebraic
curve", the small model; recalled, not re-read), with sections represented by power series at one rational point:

  * z is a rational point of E outside the support of F with n0^tw F(z) a fourth power, P0 a point of C above z,
    t = X - X(z) a local parameter at P0.  L = O(3 pi^*(z)) has degree 12 >= 2g + 1 and is stable under J and
    Frobenius.
  * H^0(L^m) = sum over j = 0..3 of L_E(3m (z) + floor(j div(F) / 4)) v^j.  The Riemann-Roch spaces on E are
    written down from chains of lines (pr_chain of sixlines_prym.py).  A section s of L^m is stored as the first
    N coefficients of t^(3m) s; products of sections are products of truncated series.
  * a class x of degree 0 is an effective divisor D of degree 12 with x = [D - 3 pi^*(z)], stored as the subspace
    W_D = H^0(L^2(-D)) of H^0(L^2) (dimension 8).
  * addflip(x1, x2) = -(x1 + x2):  H^0(L^4(-D1-D2)) = W_D1 W_D2;  H^0(L^3(-D1-D2)) = {s : s H^0(L) in it};
    pick s there, div(s) = D1 + D2 + D3;  W_D3 = {u in H^0(L^2) : u H^0(L^3(-D1-D2)) in s H^0(L^2)}.
    Every space has a dimension known from Riemann-Roch, and each step is checked against it.
  * x = 0 exactly when {u in H^0(L) : u H^0(L) in W_D} is not 0.  A canonical form (the divisor D ~ D_x with the
    largest multiplicity at P0) gives equality tests and hashing.
  * J multiplies the j-th graded piece by zeta^j; Frobenius raises coefficients to the p-th power.
  * the class [R - P0] of a rational point R: sections of L^2 that vanish on a fixed divisor D1 = div(s0) - P0
    and at R; values at R come from expressions N(x, y) / Den(x) v^j recovered from the series.
  * Tate pairing (Miller's algorithm in this model, after Bruin): a function is a ratio of two sections of the
    same L^a, and its value on an effective divisor Y is the ratio of the determinants of the two multiplication
    maps H^0(L^2) / W_Y -> H^0(L^(2+a)) / H^0(L^(2+a)(-Y)).

Validation (validate()): the order #E(F_q) |chi(1)|^2 kills random classes; J^4 = 1; Frobenius^r = 1 and commutes
with J; (J^2 + 1)(1 - J^2) = 0; chi(Frobenius) = 0 on (1 - J^2) Jac, with chi the quartic of sixlines_prym.py and
i acting as J^3 (it fails with i acting as J: the control); Frob^2 - a Frob + p = 0 on (1 + J + J^2 + J^3) Jac.

  --validate p b c e f i j [--tw k] [--r r]   one configuration (b, c, e, f) and ordered pair (i, j), all checks
  --validate-set                              the checks on a fixed set of curves (validate_set())
  --lattice                                   the lattice side of the level test (lattice_side())
  --levels [--case p e1 l] [--pairs n]        the level test on fourfolds (levels())
  --character [--torsor p e1] [--pairs n]     the form modulo 4 on the three torsors of note 19, section 4.5 (character())
  --one p b c e f i j [--tw k] [--r r]        the same on one curve
  --workers n, --out file

Run from experiments/ with python -I.  Needs sixlines_prym.py, sixlines_prym.jsonl and the files of note 2;
--lattice also needs sixlines_cm.py and sixlines_vertical.py.
"""
import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sixlines_prym as sp  # noqa: E402
from sixlines_prym import pari  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
K0 = 3          # L = O(K0 pi^*(z)), degree 4 K0 = 12
MMAX = 5        # spaces H^0(L^m) for m <= MMAX
NCOEF = 12 * MMAX + 4   # coefficients stored
GENUS = 5

pari("'t; 'X")

# ------------------------------------------------------------------------------------------------------------
# setup over F_p (the points of B may lie in F_{p^2}; everything returned has coefficients in F_p)
# ------------------------------------------------------------------------------------------------------------

pari(r"""
jc_tofp(c, p) = {
  \\ an element of the big field that lies in F_p, as an integer
  my(P);
  if(type(c) != "t_FFELT", return(lift(Mod(c, p))));
  P = c.pol; if(poldegree(P) > 0, error("coefficient not in F_p"));
  lift(Mod(polcoef(P, 0), p));
}
""")

pari(r"""
jc_newton(S, y0, e, NP) = {
  \\ the e-th root (e = 2 or 4) of the series S = y0^e + O(t), with constant term y0
  my(y = y0 + O('t^NP), it = 0);
  while(2^it < 2 * NP, y = y - (y^e - S) / (e * y^(e - 1)); it++);
  if((y^e - S) != 0, error("root of a series failed"));
  y;
}
""")

pari(r"""
jc_mons(M, Xs, Ys) = {
  \\ [exponents, series] of the monomials x^a y^b with 2 a + 3 b <= M, b <= 1: a basis of L_E(M (O))
  my(ex = List(), se = List());
  for(a = 0, M \ 2, listput(ex, [a, 0]); listput(se, Xs^a));
  for(a = 0, (M - 3) \ 2, listput(ex, [a, 1]); listput(se, Xs^a * Ys));
  [Vec(ex), Vec(se)];
}
""")

pari(r"""
jc_ellrr(E, Dfin, cO, Xs, Ys) = {
  \\ series of a basis of L_E(D), D = sum Dfin[i][2] (Dfin[i][1]) + cO (O), of degree n >= 1:
  \\ f in L(D) iff f g in L((n - 1)(O) + (S)), S the sum of D and div(g) = D - (n - 1)(O) - (S).
  my(n = cO, S = [0], L = List(), ch, g = 1, mm, xS, yS, out = List(), c);
  for(i = 1, #Dfin, n += Dfin[i][2]; S = elladd(E, S, ellmul(E, Dfin[i][1], Dfin[i][2])); if(Dfin[i][2], listput(L, Dfin[i])));
  if(n < 1, error("degree < 1"));
  if(S != [0], listput(L, [S, -1]));
  if(#L,
    ch = pr_chain(E, Vec(L));
    for(i = 1, #ch[1], g *= ch[1][i][1] * Ys + ch[1][i][2] * Xs + ch[1][i][3]);
    for(i = 1, #ch[2], g /= ch[2][i][1] * Ys + ch[2][i][2] * Xs + ch[2][i][3]));
  if(S == [0],
    mm = jc_mons(n, Xs, Ys);
    for(i = 1, #mm[2], listput(out, mm[2][i] / g)),
    xS = S[1]; yS = -S[2];
    mm = jc_mons(n + 1, Xs, Ys);
    for(i = 1, #mm[2],
      if(mm[1][i] == [0, 0], next);
      c = xS^mm[1][i][1] * yS^mm[1][i][2];
      listput(out, (mm[2][i] - c) / ((Xs - xS) * g))));
  if(#out != n, error("dimension of L_E(D)"));
  Vec(out);
}
""")

pari(r"""
jc_fit(S, Xp, Yp, M, p, NP) = {
  \\ the function with series S (coefficients in F_p, no pole at z) as A(x) + y B(x) in L_E(M (O)); [A, B]
  my(mm = jc_mons(M, Xp, Yp), n = #mm[2], Mat, rhs, sol, A = 0, B = 0);
  if(n >= NP - 30, error("not enough precision to fit"));
  Mat = matrix(NP - 30, n, i, j, polcoef(mm[2][j], i - 1));
  rhs = vectorv(NP - 30, i, polcoef(S, i - 1));
  sol = matinverseimage(Mat, rhs);
  if(#sol == 0, error("fit failed"));
  for(j = 1, n, if(mm[1][j][2], B += lift(sol[j]) * 'X^mm[1][j][1], A += lift(sol[j]) * 'X^mm[1][j][1]));
  [A, B];
}
""")

pari(r"""
jc_setup(W, bx, rat, Qc, Xc, p, tw, k0, M, N) = {
  \\ the curve C_tw : v^4 = n0^tw F of sixlines_prym.py, a base point, and H^0(L^m) for m <= M as series.
  \\ Returns [status, ...]; status 0 = degenerate.
  my(big = if(rat, p^2, p), g = ffgen(big, 'g), one = g^0, n0 = Mod(znprimroot(p), p), i4, om, r, c4, Q = [0], tgt, rho,
     a2 = polcoef(W, 2) * one, a4 = polcoef(W, 1) * one, a6 = polcoef(W, 0) * one, E = ellinit([0, a2, 0, a4, a6]),
     B4, cs, csb, zeros, DF, coefO, ch, gam, NP = N + 60, z = 0, x0, y0, w0, fz, v, rt, Xp, Yp, Xs, Ys, Fs, Fp, vs,
     V = vector(M), grad = vector(M), bas, Dfin, cO, s, col, beta, betap, ev = List(), den, fit, xc, Ffit, Wp = W * Mod(1, p), sup);
  if(p % 4 != 1, error("p = 1 mod 4 only"));
  i4 = n0^((p - 1) / 4) * one; om = if(rat, sqrt(n0 * one), 0);
  r = vector(3, j, pr_point(W, bx[j], one, om, n0));
  c4 = pr_point(W, Xc, one, om, n0);
  if(rat == 0,
    rho = ellneg(E, elladd(E, elladd(E, r[1], r[2]), r[3])),
    Q = [Qc[1] * one + Qc[2] * om, Qc[3] * one + Qc[4] * om];
    tgt = [0]; for(j = 1, 3, if(j != rat, tgt = elladd(E, tgt, r[j])));
    if(ellsub(E, Q, [Q[1]^p, Q[2]^p]) != tgt, error("Q - Frob(Q) is not a + b"));
    rho = ellsub(E, elladd(E, Q, [Q[1]^p, Q[2]^p]), r[rat]));
  B4 = concat([Q], vector(3, j, ellsub(E, elladd(E, r[j], rho), Q)));
  cs = elladd(E, c4, rho); csb = elladd(E, ellneg(E, c4), rho);
  for(i = 1, 4, if(B4[i] == cs || B4[i] == csb, return([0])); for(j = i + 1, 4, if(B4[i] == B4[j], return([0]))));
  if(cs == [0] || csb == [0], return([0]));
  if(cs[1]^p != cs[1] || cs[2]^p != cs[2], error("the pole is not rational"));
  zeros = select(P -> P != [0], B4);
  DF = concat(vector(#zeros, j, [zeros[j], 2]), [[cs, -4]]);
  coefO = 4 - 2 * #zeros;
  ch = pr_chain(E, DF);
  gam = n0^tw * one;
  \\ the base point: the first rational z = (x0, y0), y0 != 0, outside the support, with gam F(z) a fourth power
  sup = concat(zeros, [cs, ellneg(E, cs)]);
  for(xx = 0, p - 1,
    x0 = xx * one; w0 = subst(Wp, 'x, Mod(xx, p));
    if(w0 == 0 || !issquare(w0), next);
    y0 = lift(sqrt(w0)); if(2 * y0 > p, y0 = p - y0); y0 = y0 * one;
    forstep(sg = 1, -1, -2,
      if(z, break);
      if(#select(P -> P[1] == x0, sup), next);
      fz = pr_eval(ch, x0, sg * y0);
      if(fz[1] == 0 || fz[2] == 0, next);
      v = gam * fz[1] / fz[2];
      if(v^((p - 1) / 4) != one, next);
      z = [x0, sg * y0]; rt = sqrtn(Mod(jc_tofp(v, p), p), 4);
      rt = vecmin(vector(4, j, lift(rt * n0^(j * (p - 1) / 4)))));
    if(z, break));
  if(!z, return([0]));
  \\ series at z in t = X - x0
  Xp = Mod(jc_tofp(z[1], p), p) + 't + O('t^NP);
  Yp = jc_newton(subst(Wp, 'x, Xp), Mod(jc_tofp(z[2], p), p), 2, NP);
  Xs = Xp * one; Ys = Yp * one;
  Fs = gam; for(i = 1, #ch[1], Fs *= ch[1][i][1] * Ys + ch[1][i][2] * Xs + ch[1][i][3]);
  for(i = 1, #ch[2], Fs /= ch[2][i][1] * Ys + ch[2][i][2] * Xs + ch[2][i][3]);
  Fp = apply(c -> Mod(jc_tofp(c, p), p), Fs);
  vs = jc_newton(Fp, Mod(rt, p), 4, NP);
  \\ H^0(L^m) = sum_j L_E(m k0 (z) + floor(j div(F) / 4)) v^j
  beta = 1; for(i = 1, #zeros, beta *= ('X - zeros[i][1]));
  betap = apply(c -> jc_tofp(c, p), beta) * Mod(1, p);
  xc = Mod(jc_tofp(cs[1], p), p);
  for(m = 1, M,
    col = List(); grad[m] = List();
    for(j = 0, 3,
      Dfin = concat([[z, m * k0]], concat(vector(#zeros, i, [zeros[i], j \ 2]), [[cs, -j]]));
      cO = (j * coefO) \ 4;
      bas = jc_ellrr(E, Dfin, cO, Xs, Ys);
      for(i = 1, #bas,
        s = apply(c -> Mod(jc_tofp(c, p), p), bas[i]) * vs^j * 't^(m * k0);
        if(valuation(s, 't) < 0, error("a section has too large a pole"));
        listput(col, vectorv(N, n, lift(polcoef(s, n - 1))));
        listput(grad[m], j);
        if(m == 2,
          \\ the same function as (A(x) + y B(x)) / den(x) v^j, for evaluation at points
          den = ('X - jc_tofp(z[1], p))^(m * k0) * betap^(j \ 2);
          fit = jc_fit(apply(c -> Mod(jc_tofp(c, p), p), bas[i]) * subst(den, 'X, Xp), Xp, Yp, cO + 2 * poldegree(den), p, NP);
          listput(ev, [j, fit[1], fit[2], lift(den)]))));
    V[m] = Mat(Vec(col)); grad[m] = Vec(grad[m]);
    if(#V[m] != 4 * m * k0 - 4 || matrank(V[m] * Mod(1, p)) != #V[m], error("dimension of H^0(L^m)")));
  \\ F as (A + y B) / (x - xc)^4, without the factor gam
  Ffit = jc_fit(Fp / lift(Mod(jc_tofp(gam, p), p)) * (Xp - xc)^4, Xp, Yp, 12, p, NP);
  [1, V, grad, Vec(ev), [Ffit[1], Ffit[2], lift(xc)], lift(Mod(jc_tofp(z[1], p), p)), lift(betap), lift(n0), ellcard(ellinit([0, polcoef(W, 2), 0, polcoef(W, 1), polcoef(W, 0)] * Mod(1, p)))];
}
""")

# ------------------------------------------------------------------------------------------------------------
# arithmetic over F_{p^r}.  A context cx is a vector:
#   [one, N, V, grad, LI, zeta, frob, V1sub (sections of L vanishing at P0), s0, W0, WD1c, ev, Fdata, W, p, r, gam]
# ------------------------------------------------------------------------------------------------------------

pari(r"""
jc_toep(a, N) = matrix(N, N, i, j, if(i >= j, a[i - j + 1], 0 * a[1]));
""")

pari(r"""
jc_lker(U) = mattranspose(matker(mattranspose(U)));
""")

pari(r"""
jc_rows(M, ir) = mattranspose(vecextract(mattranspose(M), ir));
""")

pari(r"""
jc_rand(B) = {
  \\ a random element of the span of the columns of B
  my(n = #B, f = B[1, 1], q = f.p^f.f, c);
  c = vectorv(n, i, random(f));
  B * c;
}
""")

pari(r"""
jc_mul(A, B, N, target) = {
  \\ the span of the products a b, a in A, b in B; its dimension must be target
  my(M, T, k = 0, im);
  M = matrix(N, 0);
  while(1,
    k++;
    T = jc_toep(if(k <= 2, jc_rand(B), if(k - 2 <= #B, B[, k - 2], error("product space too small"))), N);
    M = matconcat([M, T * A]);
    if(#M < target, next);
    im = matimage(M);
    if(#im == target, return(im));
    if(#im > target, error("product space too large"));
    M = im);
}
""")

pari(r"""
jc_div(K, Wd, Vm, N, target) = {
  \\ {s in span(Vm) : s w in U for w in Wd}, U given by K with K U = 0 of full rank.  target = expected dimension
  \\ (-1: unknown, then every column of Wd is used).
  my(M = matrix(0, #Vm), k = 0, ke, w);
  while(1,
    k++;
    if(target < 0,
      if(k > #Wd, return(if(#ke, Vm * ke, matrix(N, 0))));
      w = Wd[, k],
      w = if(k <= 2, jc_rand(Wd), if(k - 2 <= #Wd, Wd[, k - 2], error("quotient space too large"))));
    M = if(k == 1, K * (jc_toep(w, N) * Vm), matconcat([M; K * (jc_toep(w, N) * Vm)]));
    ke = matker(M);
    if(target >= 0,
      if(#ke == target, return(Vm * ke));
      if(#ke < target, error("quotient space too small"))));
}
""")

pari(r"""
jc_addflip(cx, W1, W2, tr = 0) = {
  \\ W3 with [D3] = -([D1] + [D2]).  With tr: [W3, s], s the section of L^3 with div(s) = D1 + D2 + D3.
  my(N = cx[2], V = cx[3], U4, U3, s, W3);
  U4 = jc_mul(W1, W2, N, 20);
  U3 = jc_div(jc_lker(U4), V[1], V[3], N, 8);
  s = jc_rand(U3);
  W3 = jc_div(jc_lker(jc_toep(s, N) * V[2]), U3, V[2], N, 8);
  if(tr, [W3, s], W3);
}
""")

pari(r"""
jc_flipwith(cx, W, s) = {
  \\ s in W = H^0(L^2(-D)), div(s) = D + D'; returns H^0(L^2(-D'))
  my(N = cx[2], V = cx[3]);
  jc_div(jc_lker(jc_toep(s, N) * V[2]), W, V[2], N, 8);
}
""")

pari(r"""
jc_flip(cx, W, tr = 0) = { my(s = jc_rand(W), R = jc_flipwith(cx, W, s)); if(tr, [R, s], R); }
""")

pari(r"""
jc_iszero(cx, W) = {
  my(N = cx[2], V = cx[3]);
  #jc_div(jc_lker(W), V[1], V[1], N, -1) > 0;
}
""")

pari(r"""
jc_maxval(W) = {
  \\ the element of span(W) with the largest valuation (unique up to a scalar when it is returned)
  my(k = #W - 1, ke, last);
  ke = matker(W[1..k, ]);
  while(#ke, last = ke; k++; ke = matker(W[1..k, ]));
  if(#last != 1, error("maxval"));
  W * last[, 1];
}
""")

pari(r"""
jc_canon(cx, W) = {
  \\ canonical basis of the canonical representative: D'' ~ D with the largest multiplicity at P0
  my(W1 = jc_flip(cx, W), W2 = jc_flipwith(cx, W1, jc_maxval(W1)));
  matker(jc_lker(W2));
}
""")

pari(r"""
jc_J(cx, W, m, e) = {
  \\ the automorphism v -> zeta^e v on a subspace of H^0(L^m)
  my(V = cx[3], gr = cx[4][m], LI = cx[5][m], C = LI[2] * jc_rows(W, LI[1]), z = cx[6]);
  for(i = 1, #gr, C[i, ] = z^(e * gr[i]) * C[i, ]);
  V[m] * C;
}
""")

pari(r"""
jc_frob(cx, W, k = 1) = { for(i = 1, k, W = ffmap(cx[7], W)); W; }
""")

pari(r"""
jc_init(dat, W, p, r, tw, k0) = {
  \\ dat = the output of jc_setup
  my(ff = ffgen([p, r], 'u), one = ff^0, N = #dat[2][1][, 1], V = dat[2] * one, grad = dat[3], LI = vector(#V), n0 = Mod(dat[8], p),
     zeta = lift(n0^((p - 1) / 4)) * one, frob = fffrobenius(ff, 1), V1sub, s0, W0, WD1, cx, ir);
  for(m = 1, #V, ir = matindexrank(V[m])[1]; LI[m] = [ir, jc_rows(V[m], ir)^(-1)]);
  \\ sections of L vanishing at P0, and one of them, s0, with a simple zero there
  V1sub = V[1] * matker(V[1][1..1, ]);
  s0 = 0; for(i = 1, #V1sub, if(V1sub[2, i] != 0, s0 = V1sub[, i]; break));
  if(s0 === 0, error("no section with a simple zero at P0"));
  W0 = jc_toep(s0, N) * V[1];
  cx = [one, N, V, grad, LI, zeta, frob, V1sub, s0, W0, 0, dat[4], dat[5], W, p, r, lift(n0^tw) * one, dat[6], dat[7], k0];
  \\ H^0(L^2(-D1)), D1 = div(s0) - P0, in coordinates of the basis of H^0(L^2)
  WD1 = jc_div(jc_lker(jc_toep(s0, N) * V[2]), V1sub, V[2], N, 9);
  cx[11] = LI[2][2] * jc_rows(WD1, LI[2][1]);
  cx;
}
""")

pari(r"""
jc_curvepoint(cx, xx) = {
  \\ the points of C over F_{p^r} with X = xx: a list of [x, y, v], empty if none or if xx is special
  my(one = cx[1], W = cx[14], p = cx[15], r = cx[16], q = p^r, Fd = cx[13], w, y, f, v, out = List(), zeta = cx[6]);
  if(xx == Fd[3] * one || xx == cx[18] * one || subst(cx[19], 'X, xx) == 0, return([]));
  w = subst(W, 'x, xx);
  if(w == 0 || w^((q - 1) / 2) != one, return([]));
  y = sqrt(w);
  forstep(sg = 1, -1, -2,
    f = cx[17] * (subst(Fd[1], 'X, xx) + sg * y * subst(Fd[2], 'X, xx)) / (xx - Fd[3])^4;
    if(f == 0 || f^((q - 1) / 4) != one, next);
    v = sqrtn(f, 4);
    for(j = 0, 3, listput(out, [xx, sg * y, v * zeta^j])));
  Vec(out);
}
""")

pari(r"""
jc_pointclass(cx, R) = {
  \\ W_D for D = R + D1: the class [R - P0]
  my(ev = cx[12], vals, ke);
  vals = vector(#ev, i, (subst(ev[i][2], 'X, R[1]) + R[2] * subst(ev[i][3], 'X, R[1])) / subst(ev[i][4], 'X, R[1]) * R[3]^ev[i][1]);
  ke = matker(Mat(vals * cx[11]));
  if(#ke != 8, error("point on the base divisor"));
  cx[3][2] * (cx[11] * ke);
}
""")

pari(r"""
jc_randpoint(cx) = {
  \\ a random usable point of C over F_{p^r}; 0 if 400 abscissas give none
  my(one = cx[1], L);
  for(i = 1, 400, L = jc_curvepoint(cx, random(one)); if(#L, return(L[1 + random(#L)])));
  0;
}
""")

pari(r"""
jc_countpoints(cx) = {
  \\ r = 1: the number of usable points of C over F_p
  my(one = cx[1], p = cx[15]);
  sum(xx = 0, p - 1, #jc_curvepoint(cx, xx * one));
}
""")

pari(r"""
jc_descend(cxs, cxb, K) = {
  \\ a matrix over F_{p^(2r)} with entries in F_{p^r}, as a matrix over the field of the context cxs
  my(m = ffembed(cxs[1], cxb[1]), R = ffmap(ffinvmap(m), K));
  if(type(R) != "t_MAT", error("entries not in the subfield"));
  R * cxs[1];
}
""")



# --- Tate pairing -------------------------------------------------------------------------------------------

pari(r"""
jc_evalprep(cx, WY) = {
  \\ for an effective divisor Y of degree 12 given by WY = H^0(L^2(-Y)): [b, K3, K5] with b a basis of a complement
  \\ of WY in H^0(L^2), and K_m (12 rows) a linear map on series whose kernel on H^0(L^m) is H^0(L^m(-Y))
  my(N = cx[2], V = cx[3], LI = cx[5], b, U, K = vector(2), m, C);
  b = matimage(matconcat([WY, V[2]]));
  if(#b != 20, error("evalprep"));
  b = b[, 9..20];
  for(a = 1, 2,
    m = [3, 5][a];
    U = jc_mul(WY, V[m - 2], N, 12 * m - 16);
    C = LI[m][2] * jc_rows(U, LI[m][1]);
    K[a] = jc_lker(C) * LI[m][2];
    K[a] = [LI[m][1], K[a]]);
  [b, K[1], K[2]];
}
""")

pari(r"""
jc_det(cx, prep, Ts, m) = {
  \\ determinant of multiplication by the section s of L^(m-2): H^0(L^2)/W_Y -> H^0(L^m)/H^0(L^m(-Y)), m = 3 or 5;
  \\ Ts = jc_toep(s, N)
  my(K = prep[if(m == 3, 2, 3)]);
  matdet(K[2] * jc_rows(Ts * prep[1], K[1]));
}
""")

pari(r"""
jc_miller(cx, WP, n, preps) = {
  \\ P = [D_P - D_0] with n P = 0.  Returns the vector of values f(Y) / (a constant that depends only on deg Y),
  \\ for the divisors Y of preps, where div(f) = n D_P - n D_0; 0 if some factor vanishes on Y.
  my(N = cx[2], V = cx[3], k0 = cx[20], bits = binary(n), W = WP, r, r2, s, t, u, val = vector(#preps, i, cx[1]), h, sig, ke, d1, d2);
  sig = vectorv(N, i, if(i == k0 + 1, cx[1], 0 * cx[1]));     \\ the constant function 1, a section of L with divisor D_0
  h = ((s, t) -> my(Ts = jc_toep(s, N), Tt = jc_toep(jc_toep(t, N) * sig, N)); vector(#preps, i, d1 = jc_det(cx, preps[i], Ts, 5); d2 = jc_det(cx, preps[i], Tt, 5); if(d1 == 0 || d2 == 0, 0, d1 / d2)));
  for(i = 2, #bits,
    r = jc_addflip(cx, W, W, 1); r2 = jc_flip(cx, r[1], 1);
    u = h(r[2], r2[2]); for(j = 1, #preps, val[j] = val[j]^2 * u[j]);
    W = r2[1];
    if(bits[i],
      r = jc_addflip(cx, W, WP, 1); r2 = jc_flip(cx, r[1], 1);
      u = h(r[2], r2[2]); for(j = 1, #preps, val[j] = val[j] * u[j]);
      W = r2[1]));
  \\ now [W] = n P = 0: the section u of L with divisor D_(n)
  ke = jc_div(jc_lker(W), V[1], V[1], N, -1);
  if(#ke != 1, error("n P is not 0"));
  for(j = 1, #preps,
    d1 = jc_det(cx, preps[j], jc_toep(ke[, 1], N), 3); d2 = jc_det(cx, preps[j], jc_toep(sig, N), 3);
    val[j] = if(d1 == 0 || d2 == 0, 0, val[j] * d1 / d2));
  val;
}
""")

pari(r"""
jc_tate(cx, WP, WQ, n) = {
  \\ the reduced Tate pairing T_n(P, Q) = (f_P(D_Q) / f_P(D_0'))^((q - 1) / n), D_0' = div(s') for a random section s' of L
  my(N = cx[2], V = cx[3], p = cx[15], r = cx[16], sp, pr, v, WQ2 = WQ);
  if((p^r - 1) % n, error("n does not divide q - 1"));
  for(try = 1, 8,
    sp = jc_rand(V[1]);
    pr = [jc_evalprep(cx, WQ2), jc_evalprep(cx, jc_toep(sp, N) * V[1])];
    v = jc_miller(cx, WP, n, pr);
    if(v[1] != 0 && v[2] != 0, return((v[1] / v[2])^((p^r - 1) / n)));
    WQ2 = jc_flip(cx, jc_flip(cx, WQ2)));
  error("Tate pairing: supports not disjoint");
}
""")

pari(r"""
jc_tatemat(cx, WPs, WQs, n) = {
  \\ the matrix of T_n(P_a, Q_b); an entry is 0 where the supports met (recompute it with jc_tate)
  my(N = cx[2], V = cx[3], p = cx[15], r = cx[16], sp = jc_rand(V[1]), preps, v, M = matrix(#WPs, #WQs));
  preps = concat(vector(#WQs, b, jc_evalprep(cx, jc_flip(cx, jc_flip(cx, WQs[b])))), [jc_evalprep(cx, jc_toep(sp, N) * V[1])]);
  for(a = 1, #WPs,
    v = jc_miller(cx, WPs[a], n, preps);
    for(b = 1, #WQs, M[a, b] = if(v[b] == 0 || v[#preps] == 0, 0, (v[b] / v[#preps])^((p^r - 1) / n))));
  M;
}
""")

pari(r"""
jc_rootof1(cx, n) = {
  \\ a primitive n-th root of unity z in F_{p^r} with z^(n/4) = zeta (n a power of 2, n >= 4)
  my(g = ffprimroot(cx[1]), p = cx[15], r = cx[16], w = g^((p^r - 1) / n));
  forstep(t = 1, n, 2, if((w^t)^(n / 4) == cx[6], return(w^t)));
  error("no root of unity");
}
""")

pari(r"""
jc_dlog(x, zeta, n) = { for(k = 0, n - 1, if(zeta^k == x, return(k))); error("not a power of zeta"); }
""")

GP = {n: pari(n) for n in ["jc_setup", "jc_init", "jc_addflip", "jc_flip", "jc_iszero", "jc_canon", "jc_J", "jc_frob",
                           "jc_curvepoint", "jc_pointclass", "jc_randpoint", "jc_rand", "jc_tate", "jc_dlog", "jc_tatemat", "jc_rootof1", "jc_countpoints", "jc_descend"]}


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gpow_i(k):
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][k % 4]


def curve_input(p, cfg, pair):
    """(W, bx, rat, Qc, Xc) of sixlines_prym.py for this configuration and ordered pair, or None."""
    rest = [i for i in range(1, 7) if i not in pair]
    idx = rest + list(pair)
    L = pari.matrix(6, 3, [x for row in sp.forms(cfg) for x in row])
    N, D = sp.GP["pr_lambda"](L, idx, p)
    cur = sp.GP["pr_curve"](N, D, p, 0)
    if int(cur[0]) == 0:
        return None
    W, bx, rat, Xc, yrat = cur[1], cur[2], int(cur[3]), cur[4], bool(int(cur[5]))
    if yrat != (p % 4 == 1) or p % 4 != 1:
        return None
    Qc = sp.GP["pr_findq"](W, bx, rat, Xc, p) if rat else pari("[0, 0, 0, 0]")
    return W, bx, rat, Qc, Xc


def chi_twist(c, p, tw):
    """coefficients (1, -c1, c2, -c3, p^2) -> list [a0..a4] of Gaussian integers of chi for C_tw, low degree first."""
    c1, c2, c3 = [tuple(v) for v in c]
    c1, c2, c3 = gmul(gpow_i(tw), c1), gmul(gpow_i(2 * tw), c2), gmul(gpow_i(3 * tw), c3)
    return [(p * p, 0), (-c3[0], -c3[1]), c2, (-c1[0], -c1[1]), (1, 0)]


class Jac:
    """Jac(C_tw)(F_{p^r}) for one configuration and pair.  Classes are GP matrices (bases of W_D)."""

    def __init__(self, p, cfg, pair, tw=0, r=1, seed=1):
        self.p, self.cfg, self.pair, self.tw, self.r = p, list(cfg), tuple(pair), tw, r
        inp = curve_input(p, cfg, pair)
        if inp is None:
            raise ValueError("degenerate or not of commuting type")
        self.W, self.rat = inp[0], inp[2]
        dat = GP["jc_setup"](inp[0], inp[1], inp[2], inp[3], inp[4], p, tw, K0, MMAX, NCOEF)
        if int(dat[0]) == 0:
            raise ValueError("degenerate")
        self.dat = dat
        self.cardE = int(dat[8])
        pari.setrand(seed)
        self.cx = GP["jc_init"](dat, self.W, p, r, tw, K0)
        self.zero = self.cx[9]

    def with_r(self, r):
        new = object.__new__(Jac)
        new.__dict__.update(self.__dict__)
        new.r = r
        new.cx = GP["jc_init"](self.dat, self.W, self.p, r, self.tw, K0)
        new._roots = {}
        for name in ("_few", "_big"):
            new.__dict__.pop(name, None)
        new.zero = new.cx[9]
        return new

    # group law
    def addflip(self, a, b):
        return GP["jc_addflip"](self.cx, a, b)

    def neg(self, a):
        return GP["jc_flip"](self.cx, a)

    def add(self, a, b):
        return self.neg(self.addflip(a, b))

    def sub(self, a, b):
        return self.addflip(self.neg(a), b)

    def is_zero(self, a):
        return bool(GP["jc_iszero"](self.cx, a))

    def eq(self, a, b):
        return self.is_zero(self.sub(a, b))

    def key(self, a):
        return str(GP["jc_canon"](self.cx, a))

    def mul(self, n, a):
        if n == 0:
            return self.zero
        if n < 0:
            n, a = -n, self.neg(a)
        na = None
        W, sg = a, 1
        for bit in bin(n)[3:]:
            W, sg = self.addflip(W, W), -sg
            if bit == "1":
                if sg == 1:
                    W, sg = self.addflip(W, a), -1
                else:
                    if na is None:
                        na = self.neg(a)
                    W, sg = self.addflip(W, na), 1
        return W if sg == 1 else self.neg(W)

    # endomorphisms
    def J(self, a, e=1):
        return GP["jc_J"](self.cx, a, 2, e % 4)

    def frob(self, a, k=1):
        return GP["jc_frob"](self.cx, a, k)

    def gauss(self, g, a, s=1):
        """(g0 + g1 J^s) a"""
        out = self.mul(g[0], a)
        if g[1]:
            out = self.add(out, self.mul(g[1], self.J(a, s)))
        return out

    def poly(self, coeffs, a, s=1):
        """sum_k coeffs[k](J^s) Frob^k a, coeffs Gaussian integers, low degree first (Horner)."""
        out = self.zero
        for g in reversed(coeffs):
            out = self.add(self.frob(out), self.gauss(g, a, s))
        return out

    # points
    def point(self, R):
        return GP["jc_pointclass"](self.cx, R)

    def random(self, n=3):
        """A random class: a sum of n classes [R - P0], R a point of C over F_{p^r}.  If C has too few such points
        (fewer than 12 over F_p, or none found over a larger field), the trace of such a sum over F_{p^(2r)}."""
        if not hasattr(self, "_few"):
            self._few = int(GP["jc_countpoints"](self.cx)) < 12 if self.r == 1 else GP["jc_randpoint"](self.cx) == 0
        if self._few:
            if not hasattr(self, "_big"):
                self._big = self.with_r(2 * self.r)
            big = self._big
            x = big.random(n)
            y = big.add(x, big.frob(x, self.r))
            return GP["jc_descend"](self.cx, big.cx, GP["jc_canon"](big.cx, y))
        out = self.zero
        for _ in range(n):
            for attempt in range(50):
                R = GP["jc_randpoint"](self.cx)
                if R == 0:
                    continue
                try:
                    cl = self.point(R)
                    break
                except Exception:
                    continue
            else:
                raise RuntimeError("no usable point")
            out = self.add(out, cl)
        return out

    def tate(self, P, Q, n):
        """the reduced Tate pairing T_n(P, Q), an n-th root of unity in F_{p^r}"""
        return GP["jc_tate"](self.cx, P, Q, n)

    def root_of_unity(self, n):
        """one primitive n-th root of unity z with z^(n/4) = zeta, the same in every call (ffprimroot is random)"""
        if not hasattr(self, "_roots"):
            self._roots = {}
        if n not in self._roots:
            self._roots[n] = GP["jc_rootof1"](self.cx, n)
        return self._roots[n]

    def tate_matrix(self, Ps, Qs, n):
        """logarithms of T_n(P_a, Q_b) to the base of the n-th root of unity z with z^(n/4) = zeta"""
        M = GP["jc_tatemat"](self.cx, Ps, Qs, n)
        z = self.root_of_unity(n)
        out = []
        for a in range(len(Ps)):
            row = []
            for b in range(len(Qs)):
                v = M[a, b]
                if v == 0:
                    v = self.tate(Ps[a], Qs[b], n)
                row.append(int(GP["jc_dlog"](v, z, n)))
            out.append(row)
        return out

    # orders
    def card_E(self):
        a = self.p + 1 - self.cardE
        return int(pari(f"polresultant(x^2 - ({a})*x + {self.p}, x^{self.r} - 1)"))

    def chi(self, c):
        return chi_twist(c, self.p, self.tw)

    def card_prym_gauss(self, c):
        ch = self.chi(c)
        pol = " + ".join(f"(({a}) + ({b})*y) * x^{k}" for k, (a, b) in enumerate(ch))
        res = pari(f"lift(polresultant(Mod(1, y^2 + 1) * ({pol}), x^{self.r} - 1, x))")
        return int(pari.polcoef(res, 0)), int(pari.polcoef(res, 1))

    def card(self, c):
        g = self.card_prym_gauss(c)
        return self.card_E() * (g[0] ** 2 + g[1] ** 2)



# ------------------------------------------------------------------------------------------------------------
# structure of a finite abelian l-group given by a projection of random classes
# ------------------------------------------------------------------------------------------------------------

class LGroup:
    """The l-group H of known order l^w generated by project(random class).  Relations between the generators
    are found by discrete logarithms in the subgroup generated so far (brute force on its elements of order l),
    and the invariants are those of the relation matrix."""

    def __init__(self, A, ell, w, project, max_gens=14):
        self.A, self.ell, self.w, self.project = A, ell, w, project
        self.ys, self.R = [], []          # generators kept, relation rows over them
        self.basis = []                   # (element, exponent): independent, orders l^exponent
        self.T = []                       # basis[i] = sum T[i][k] ys[k]
        self.table = None
        self.complete = False
        self.tried = 0
        for _ in range(max_gens):
            if self.size() == ell ** w:
                self.complete = True
                break
            self.add(project(A.random()))
            self.tried += 1
        self.complete = self.size() == ell ** w
        if self.size() > ell ** w:
            raise RuntimeError("subgroup larger than the expected order")

    def size(self):
        return self.ell ** sum(e for _, e in self.basis)

    def invariants(self):
        return sorted((e for _, e in self.basis), reverse=True)

    def order_exp(self, z):
        k = 0
        while not self.A.is_zero(z):
            z = self.A.mul(self.ell, z)
            k += 1
            if k > self.w:
                raise RuntimeError("element of order larger than the group")
        return k

    def socle_table(self):
        if self.table is None:
            A, ell = self.A, self.ell
            tab = {A.key(A.zero): ()}
            cur = [((), A.zero)]
            for b, e in self.basis:
                s = A.mul(ell ** (e - 1), b)
                nxt = []
                for co, el in cur:
                    x = el
                    for c in range(ell):
                        if c:
                            x = A.add(x, s)
                        nxt.append((co + (c,), x))
                cur = nxt
            self.table = {A.key(el): co for co, el in cur}
        return self.table

    def dlog(self, z):
        """coefficients of z in the basis, or None if z is not in the subgroup"""
        A, ell = self.A, self.ell
        co = [0] * len(self.basis)
        cur = z
        while True:
            k = self.order_exp(cur)
            if k == 0:
                return co
            if not self.basis:
                return None
            c = self.socle_table().get(A.key(A.mul(ell ** (k - 1), cur)))
            if c is None:
                return None
            for i, (b, e) in enumerate(self.basis):
                if c[i]:
                    if e < k:
                        return None
                    d = c[i] * ell ** (e - k)
                    co[i] += d
                    cur = A.sub(cur, A.mul(d, b))

    def add(self, y):
        A, ell = self.A, self.ell
        z, j = y, 0
        while True:
            co = self.dlog(z)
            if co is not None:
                break
            z = A.mul(ell, z)
            j += 1
        if j == 0:
            return
        k = len(self.ys)
        row = [-sum(co[i] * self.T[i][l] for i in range(len(self.basis))) for l in range(k)] + [ell ** j]
        self.R = [r + [0] for r in self.R] + [row]
        self.ys.append(y)
        k += 1
        U, V, D = pari.matsnf(pari.matrix(k, k, [x for r in self.R for x in r]), 1)
        Vi = V ** (-1)
        mod = ell ** self.w
        self.basis, self.T = [], []
        for i in range(k):
            d = int(D[i, i])
            if abs(d) == 1:
                continue
            e = 0
            while d % ell == 0:
                d //= ell
                e += 1
            assert abs(d) == 1
            Trow = [int(Vi[i, l]) % mod for l in range(k)]
            b = A.zero
            for l in range(k):
                if Trow[l]:
                    b = A.add(b, A.mul(Trow[l], self.ys[l]))
            self.basis.append((b, e))
            self.T.append(Trow)
        self.table = None


def val(n, ell):
    k = 0
    while n and n % ell == 0:
        n //= ell
        k += 1
    return k


def prym_ell_group(A, c, ell):
    """Structure of the l-primary part of the group of rational points of the part with J^2 = -1 (l odd)."""
    n = A.card(c)
    g = A.card_prym_gauss(c)
    w = val(g[0] ** 2 + g[1] ** 2, ell)
    m = n // ell ** val(n, ell)

    def project(x):
        y = A.mul(m, x)
        return A.sub(y, A.J(y, 2))

    return LGroup(A, ell, w, project)


# ------------------------------------------------------------------------------------------------------------
# step 2: levels
# ------------------------------------------------------------------------------------------------------------

LEVEL_CASES = [(29, [14, 87, -23548], 13), (37, [34, 1591, 38332], 5), (37, [66, 4551, 169756], 3)]

pari(r"""
jc_rootsat(c1, c2, c3, p, l, tw) = {
  \\ chi_tw modulo each prime above l: [[residue degree, [root, multiplicity, order of the root], ...], ...]
  my(out = List(), T, iis, chi0 = lift('X^4 - c1 * 'X^3 + c2 * 'X^2 - c3 * 'X + p^2), chi, fa, ii, ro, row);
  if(l % 4 == 1,
    iis = [lift(sqrt(Mod(-1, l))) * Mod(1, l), -lift(sqrt(Mod(-1, l))) * Mod(1, l)],
    T = ffgen(('y^2 + 1) * Mod(1, l), 'w); iis = [T]);
  for(k = 1, #iis,
    ii = iis[k];
    chi = subst(chi0, 'y, ii) * ii^0;
    fa = factor(chi);
    row = List([if(l % 4 == 1, 1, 2)]);
    for(j = 1, #fa~,
      if(poldegree(fa[j, 1]) == 1,
        ro = -polcoef(fa[j, 1], 0) * ii^tw;
        listput(row, [fa[j, 2], if(l % 4 == 1, znorder(ro), fforder(ro))])));
    listput(out, Vec(row)));
  Vec(out);
}
""")


def gstr(v):
    return f"Mod(({v[0]}) + ({v[1]}) * y, y^2 + 1)"


def roots_at(c, p, ell, tw):
    r = pari("jc_rootsat")(pari(gstr(c[0])), pari(gstr(c[1])), pari(gstr(c[2])), p, ell, tw)
    return [(int(row[0]), [(int(x[0]), int(x[1])) for x in list(row)[1:]]) for row in r]


def level_pair(job):
    """One configuration and pair: the twist on which a multiple root of chi modulo a prime above l becomes 1,
    and the structure of the l-part of the rational points of the fourfold of that twist over F_p."""
    p, cfg, pair, c, ell = job
    out = {"p": p, "cfg": cfg, "pair": list(pair), "c": c, "l": ell}
    best = None
    for tw in range(4):
        ra = roots_at(c, p, ell, tw)
        mult = [m for f, rows in ra for m, o in rows if o == 1 and m >= 2]
        if mult:
            best = (tw, ra)
            break
    if best is None:
        out["status"] = "no twist with a multiple root equal to 1"
        return out
    tw, ra = best
    try:
        A = Jac(p, cfg, tuple(pair), tw, 1, seed=1)
    except ValueError as err:
        out["status"] = str(err)
        return out
    G = prym_ell_group(A, c, ell)
    # roots equal to 1: (residue degree, multiplicity)
    ones = [(f, m) for f, rows in ra for m, o in rows if o == 1]
    simple = sum(f for f, m in ones if m == 1)
    fm = [(f, m) for f, m in ones if m >= 2]
    inv = G.invariants()
    out.update({"status": "ok" if G.complete else "group not generated", "tw": tw, "w": G.w, "invariants": inv,
                "rank": len(inv), "roots_equal_1": ones, "generators_tried": G.tried})
    if len(fm) == 1 and G.complete:
        f, m = fm[0]
        blocks, rem = divmod(len(inv) - simple, f)
        out["multiplicity"] = m
        out["blocks"] = blocks if rem == 0 else None
    return out


def surfaces_of(p, e):
    """the configurations with polynomial P(t) or P(-t), grouped into surfaces by association"""
    rows = {tuple(r["cfg"]): r for r in map(json.loads, open(os.path.join(HERE, f"sixlines_p{p}.jsonl")))}
    cfgs = []
    for cfg, r in rows.items():
        ee = sp.surface_e(r["S"])
        if list(ee) == list(e) or list(ee) == [-e[0], e[1], -e[2]]:
            cfgs.append(cfg)
    surf, seen = [], set()
    for cfg in sorted(cfgs):
        if cfg in seen:
            continue
        a = tuple(rows[cfg]["assoc"])
        seen.update([cfg, a])
        surf.append(sorted({cfg, a}))
    return surf


def stored_pairs(p):
    out = defaultdict(list)
    for r in map(json.loads, open(os.path.join(HERE, "sixlines_prym.jsonl"))):
        if r["p"] == p:
            for x in r["pairs"]:
                if x.get("compare") in ("equal", "twist"):
                    out[tuple(r["cfg"])].append((tuple(x["pair"]), [list(v) for v in x["c"]], x["compare"], x["rational_point"]))
    return out


def levels(workers, out_path, cases=None, max_pairs=8):
    fh = open(out_path, "w") if out_path else None
    for p, e, ell in LEVEL_CASES:
        if cases and (p, e[0], ell) not in cases:
            continue
        surf = surfaces_of(p, e)
        sp_ = stored_pairs(p)
        jobs, owner = [], []
        for si, cfgs in enumerate(surf):
            for cfg in cfgs:
                for pair, c, comp, rat in sp_[tuple(cfg)][:max_pairs]:
                    jobs.append((p, list(cfg), pair, c, ell))
                    owner.append(si)
        t0 = time.time()
        with Pool(workers) as pool:
            recs = pool.map(level_pair, jobs, chunksize=1)
        per = defaultdict(Counter)
        status = Counter()
        for si, rec in zip(owner, recs):
            rec["surface"] = si
            rec["e"] = e
            status[rec["status"]] += 1
            if rec["status"] == "ok":
                per[si][(rec.get("blocks"), tuple(rec["invariants"]))] += 1
            if fh:
                fh.write(json.dumps(rec) + "\n")
        if fh:
            fh.flush()
        print(f"p = {p}, e = {e}, l = {ell}: {len(surf)} surfaces, {sum(len(x) for x in surf)} configurations, "
              f"{len(jobs)} pairs, status {dict(status)}  ({time.time() - t0:.0f} s)")
        dist = Counter()
        for si, cfgs in enumerate(surf):
            vals = sorted(per[si].items(), key=str)
            blocks = {k[0] for k, n in vals}
            dist[tuple(sorted(blocks, key=str))] += 1
            print(f"  surface {si} {['-'.join(map(str, c)) for c in cfgs]}: "
                  + "; ".join(f"blocks {k[0]}, group {list(k[1])}: {n} pairs" for k, n in vals), flush=True)
        print(f"  surfaces by number of blocks at the multiple root: {dict(dist)}", flush=True)


# ------------------------------------------------------------------------------------------------------------
# step 2, lattice side: Jordan types of F on L / l L, and of the second exterior power of pi on M / l M
# ------------------------------------------------------------------------------------------------------------

ALL_LEVEL_CASES = [(29, [14, 87, -23548], 13), (29, [18, 87, -30276], 3), (37, [14, -185, -38332], 11),
                   (37, [34, 1591, 38332], 5), (37, [38, 1591, 27380], 7), (37, [66, 4551, 169756], 3)]

pari(r"""
jl_kdims(M, P, l) = {
  \\ for each irreducible factor g of P mod l with multiplicity >= 2: [g, deg, multiplicity, [dim ker g(M)^k / deg]]
  my(fa = factor(P * Mod(1, l)), out = List(), g, G, n = #M, d);
  for(j = 1, #fa~,
    if(fa[j, 2] < 2, next);
    g = fa[j, 1]; d = poldegree(g);
    G = subst(lift(g), variable(P), M);
    listput(out, [lift(g), d, fa[j, 2], vector(fa[j, 2], k, (n - matrank(G^k)) / d)]));
  Vec(out);
}
""")

pari(r"""
jl_wedge(A) = {
  my(idx = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]);
  matrix(6, 6, r, c, A[idx[r][1], idx[c][1]] * A[idx[r][2], idx[c][2]] - A[idx[r][1], idx[c][2]] * A[idx[r][2], idx[c][1]]);
}
""")

pari(r"""
jl_fourfold(c1, c2, c3, p, l, P) = {
  \\ chi modulo a prime above l.  For every Jordan type of pi at its multiple root (the other roots semisimple):
  \\ [partition, jl_kdims of the second exterior power, order of the root for the twists 0..3]
  my(chi0 = lift('X^4 - c1 * 'X^3 + c2 * 'X^2 - c3 * 'X + p^2), T, ii, chi, fa, out = List(), blocks, a, m = 0, rest = List(), A, parts, k, one, ord);
  if(l % 4 == 1, ii = lift(sqrt(Mod(-1, l))) * Mod(1, l); one = Mod(1, l), T = ffgen(('y^2 + 1) * Mod(1, l), 'w); ii = T; one = T^0);
  chi = subst(chi0, 'y, ii) * one;
  fa = factor(chi);
  for(j = 1, #fa~,
    if(fa[j, 2] >= 2, a = -polcoef(fa[j, 1], 0); m = fa[j, 2], for(e = 1, fa[j, 2], listput(rest, matcompanion(fa[j, 1])))));
  if(m == 0, return([]));
  ord = vector(4, tw, if(l % 4 == 1, znorder(a * ii^(tw - 1)), fforder(a * ii^(tw - 1))));
  parts = partitions(m);
  for(i = 1, #parts,
    blocks = List(Vec(rest));
    for(b = 1, #parts[i], k = parts[i][b]; listput(blocks, matrix(k, k, r, s, if(r == s, a, if(s == r + 1, one, 0 * one)))));
    A = matconcat(matdiagonal(Vec(blocks)));
    listput(out, [Vec(parts[i]), jl_kdims(jl_wedge(A), P, l), ord]));
  Vec(out);
}
""")


def lattice_side():
    """For each polynomial whose surfaces lie on two levels at an odd prime l: the lattice classes of notes 2 and
    7 with the index of End(L) at l and the dimensions of ker g(F)^k on L / l L; and the same dimensions for the
    second exterior power of pi with each Jordan type at the multiple root of chi.  A class and a Jordan type
    belong together when the dimensions agree."""
    import sixlines_cm as cm
    import sixlines_vertical as sv
    prym = [json.loads(line) for line in open(os.path.join(HERE, "sixlines_prym.jsonl"))]
    for q, e, ell in ALL_LEVEL_CASES:
        P = cm.pol(cm.charpoly(q, e))
        bnf = pari.bnfinit(P, 1)
        nf = cm.GP_nf(bnf)
        w = int(pari("(b)->b.tu[1]")(bnf))
        zeta = pari("(b)->nfbasistoalg(b.nf, b.tu[2])")(bnf)
        classes = cm.lattice_classes(bnf, q, cm.GP_sl_gens(nf, q), w)
        F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(nf)
        c = next(x["c"] for r in prym if r["p"] == q and r["e"] == e for x in r["pairs"] if x["compare"] == "equal")
        four = pari("jl_fourfold")(pari(gstr(c[0])), pari(gstr(c[1])), pari(gstr(c[2])), q, ell, P)
        types = {str(kd): [int(v) for v in part] for part, kd, od in four}
        orders = [int(v) for v in four[0][2]]
        tally = Counter()
        for L, al in classes:
            O = sv.GP_order(nf, L)
            lev = int(pari.valuation(pari.matdet(sv.GP_local(O, ell)), ell))
            aut = int(sv.GP_aut(nf, L, zeta, w))
            M = pari("(nf, L, F, l)->(L^(-1) * sl_mulmat(nf, F) * L) * Mod(1, l)")(nf, L, F, ell)
            kd = pari("jl_kdims")(M, P, ell)
            tally[(lev, aut, str([[int(v) for v in row[3]] for row in kd]), str(types.get(str(kd))))] += 1
        print(f"q = {q}, e = {e}, l = {ell}: {len(classes)} lattice classes; "
              f"P mod l has factors of degree and multiplicity {[(int(pari.poldegree(g)), int(m)) for g, m in zip(*pari.factor(P * pari.Mod(1, ell)))]}; "
              f"order of the multiple root of chi_tw, tw = 0..3: {orders}")
        for (lev, aut, kd, part), n in sorted(tally.items()):
            print(f"   {n} classes: index of End(L) at l = l^{lev}, {aut} automorphisms, dim ker g(F)^k = {kd}, "
                  f"Jordan type of pi at the multiple root: {part}")
        print(f"   second exterior power, by Jordan type: " + "; ".join(f"{[int(v) for v in part]} -> {[[int(v) for v in row[3]] for row in kd]}" for part, kd, od in four), flush=True)


# ------------------------------------------------------------------------------------------------------------
# step 3: the hermitian form at the prime of degree 1 above 2, modulo 4
# ------------------------------------------------------------------------------------------------------------

CHARACTER_CASES = [(29, [2, 87, -3364]), (37, [2, -185, -5476]), (37, [18, -185, -49284])]


def gadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def gmod(u, m):
    return (u[0] % m, u[1] % m)


def gval(x):
    """valuation at 1 + i"""
    if x[0] == 0 and x[1] == 0:
        return 999
    v = 0
    while (x[0] + x[1]) % 2 == 0:
        x = ((x[0] + x[1]) // 2, (x[1] - x[0]) // 2)
        v += 1
    return v


def gdiv1pi(x, k):
    for _ in range(k):
        x = ((x[0] + x[1]) // 2, (x[1] - x[0]) // 2)
    return x


def gpow(x, n, m):
    r = (1, 0)
    for _ in range(n):
        r = gmod(gmul(r, x), m)
    return r


def ginv(u, mod):
    n = (u[0] * u[0] + u[1] * u[1]) % mod
    ni = pow(n, -1, mod)
    return ((u[0] * ni) % mod, (-u[1] * ni) % mod)


def geval(co, x, m):
    r = (0, 0)
    for c in reversed(co):
        r = gmod(gadd(gmul(r, x), c), m)
    return r


def dyadic_root(co, B=40):
    """The roots of chi (Gaussian coefficients, low degree first) in Z_2[i], modulo 2^B, and the degrees of the
    factors of chi chibar over Q_2.  One root when the real quartic field has one prime of degree 1 above 2."""
    pol = " + ".join(f"(({a}) + ({b})*y) * X^{k}" for k, (a, b) in enumerate(co))
    N = pari(f"my(f = {pol}); lift(Mod(f, y^2+1) * Mod(subst(f, y, -y), y^2+1))")
    N = pari.polcoef(N, 0, "y")
    facs = list(pari.factorpadic(N, 2, 3 * B)[0])
    out = []
    for g in facs:
        if int(pari.poldegree(g)) != 2:
            continue
        a = -pari.polcoef(g, 1) / 2
        d = pari.polcoef(g, 0) - a * a
        if not pari.issquare(d):
            continue
        b = pari.sqrt(d)
        a, b = int(pari.lift(pari.Mod(a, 2 ** (2 * B)))), int(pari.lift(pari.Mod(b, 2 ** (2 * B))))
        for sg in (1, -1):
            x = (a % 2 ** B, (sg * b) % 2 ** B)
            if gval(geval(co, x, 2 ** B)) >= B:
                out.append(x)
    return out, [int(pari.poldegree(g)) for g in facs]


def character_pair(job):
    """One configuration and pair.  Over F_{p^r} with N = 2^v | p^r - 1: P-side elements z_a = (1 - J) y_a of the
    2-primary group, the kernel of Frobenius - pi_1 on their span by linear algebra on Tate pairings, and for X in
    it Z = log T_N(X, iX) + i log T_N(X, X), which is M conj(N / (pi_1^r - 1)) with M = N(gamma) h.
    The class of h modulo norms from Z_2[i], that is the odd part of M modulo 4, is the invariant."""
    p, cfg, pair, c, tw, r, mm, seed = job
    out = {"p": p, "cfg": cfg, "pair": list(pair), "c": c, "tw": tw, "r": r}
    co = chi_twist(c, p, tw)
    B = 40
    rts, degs = dyadic_root(co, B)
    out["dyadic_degrees"] = degs
    if len(rts) != 1:
        out["status"] = "no single root in Z_2[i]"
        return out
    x = rts[0]
    mod = 2 ** B
    varpi = gadd(gpow(x, r, mod), (-1, 0))
    s = gval(varpi)
    nn = val(p ** r - 1, 2)
    N = 2 ** nn
    if s > 2 * nn:
        out["status"] = "N / (pi_1^r - 1) is not integral"
        return out
    try:
        A = Jac(p, cfg, tuple(pair), tw, r, seed=seed)
    except ValueError as err:
        out["status"] = str(err)
        return out
    n = A.card(c)
    m = n >> val(n, 2)
    t0 = time.time()
    ys = [A.mul(m, A.random()) for _ in range(mm)]
    zs = [A.sub(y, A.J(y, 1)) for y in ys]
    if not all(A.is_zero(A.mul(N, z)) for z in zs):
        out["status"] = "exponent of (1 - J) G larger than N"
        return out
    G0 = A.tate_matrix(zs, ys, N)
    G1 = A.tate_matrix([A.J(z, 1) for z in zs], ys, N)
    G3 = A.tate_matrix([A.J(z, 3) for z in zs], ys, N)
    Gp = A.tate_matrix([A.frob(z) for z in zs], ys, N)
    al, be = x[0] % N, x[1] % N
    K = [[(Gp[a][b] - al * G0[a][b] - be * G3[a][b]) % N for b in range(mm)] for a in range(mm)]
    ker = pari("(K, N) -> my(m = #K); matkerint(matconcat([mattranspose(K), N * matid(m)]))[1..m, ]")(
        pari.matrix(mm, mm, [v for row in K for v in row]), N)
    vecs = [[int(ker[i, j]) % N for i in range(mm)] for j in range(len(ker))]

    def order(cv):
        img = [sum(cv[a] * G0[a][b] for a in range(mm)) % N for b in range(mm)]
        o = N
        while o > 1 and all((o // 2) * v % N == 0 for v in img):
            o //= 2
        return o

    def zof(cv):
        re = sum(cv[a] * cv[b] * (G1[a][b] - G0[a][b]) for a in range(mm) for b in range(mm)) % N
        im = sum(cv[a] * cv[b] * (G0[a][b] - G3[a][b]) for a in range(mm) for b in range(mm)) % N
        return (re, im)

    ups = gdiv1pi(varpi, s)
    w = gmul(gdiv1pi((N, 0), s), ginv(gmod(ups, mod), mod))      # N / varpi
    wb = (w[0] % N, (-w[1]) % N)

    def solve(Z):
        return sorted({M for M in range(N) if all((M * wb[k] - Z[k]) % N == 0 for k in (0, 1))})

    vecs = [v for v in vecs if order(v) > 1]
    tally = Counter()
    best = None
    for cv in vecs:
        Z = zof(cv)
        sol = solve(Z)
        tally[(order(cv), tuple(sol))] += 1
        if len(sol) == 1 and sol[0] and (best is None or val(sol[0], 2) < val(best[1], 2)):
            best = (cv, sol[0], Z)
    out.update({"s": s, "N": N, "kernel": sorted(([o, list(sol), k] for (o, sol), k in tally.items()), key=str),
                "unsolvable": sum(k for (o, sol), k in tally.items() if not sol)})
    if best is None:
        out["status"] = "no kernel element with a non-zero self-pairing"
        return out
    cv, M, Z = best
    X = A.zero
    for a in range(mm):
        if cv[a]:
            X = A.add(X, A.mul(cv[a], zs[a]))
    o = order(cv)
    d = A.sub(A.frob(X), A.add(A.mul(al, X), A.mul(be, A.J(X, 3))))
    Zd = (A.tate_matrix([X], [A.J(X, 3)], N)[0][0], A.tate_matrix([X], [X], N)[0][0])
    out["verified"] = {"order": A.is_zero(A.mul(o, X)) and not A.is_zero(A.mul(o // 2, X)), "kernel": A.is_zero(d),
                       "J^2 = -1": A.is_zero(A.add(X, A.J(X, 2))), "Z": list(Zd) == list(Z)}
    v2 = val(M, 2)
    out.update({"M": M, "v2_M": v2, "order_X": o, "seconds": round(time.time() - t0)})
    if v2 + 2 > nn:
        out["status"] = "M known to too few digits"
        return out
    u = (M >> v2) % 4
    c1, c2 = c[0], c[1]
    sg = (1 if c1[0] * c1[1] > 0 else -1) * (1 if c2[0] > 0 else -1)
    out.update({"u": u, "delta": (1 if u == 1 else -1) * sg, "status": "ok" if all(out["verified"].values()) else "not verified"})
    return out


def _guarded_character(job):
    try:
        return character_pair(job)
    except Exception as err:          # reported as a status, counted in the summary
        return {"p": job[0], "cfg": job[1], "pair": list(job[2]), "c": job[3], "tw": job[4], "r": job[5],
                "status": "error: " + str(err)[:200]}


def character(workers, out_path, cases=None, max_pairs=2, tw=0, r=4, mm=12):
    fh = open(out_path, "w") if out_path else None
    for p, e in CHARACTER_CASES:
        if cases and (p, e[0]) not in cases:
            continue
        surf = surfaces_of(p, e)
        sp_ = stored_pairs(p)
        jobs, owner = [], []
        for si, cfgs in enumerate(surf):
            for cfg in cfgs:
                for pair, c, comp, rat in sp_[tuple(cfg)][:max_pairs]:
                    jobs.append((p, list(cfg), pair, c, tw, r, mm, 1))
                    owner.append(si)
        t0 = time.time()
        with Pool(workers) as pool:
            recs = pool.map(_guarded_character, jobs, chunksize=1)
        per = defaultdict(list)
        status = Counter()
        for si, rec in zip(owner, recs):
            rec["surface"] = si
            rec["e"] = e
            status[rec["status"]] += 1
            per[si].append(rec)
            if fh:
                fh.write(json.dumps(rec) + "\n")
        if fh:
            fh.flush()
        print(f"p = {p}, e = {e}: {len(surf)} surfaces, {sum(len(x) for x in surf)} configurations, {len(jobs)} pairs, "
              f"status {dict(status)}  ({time.time() - t0:.0f} s)")
        dist = Counter()
        for si, cfgs in enumerate(surf):
            good = [r_ for r_ in per[si] if r_["status"] == "ok"]
            us = Counter(r_["u"] for r_ in good)
            ds = Counter(r_["delta"] for r_ in good)
            dist[tuple(sorted(ds))] += 1
            print(f"  surface {si} {['-'.join(map(str, c)) for c in cfgs]}: u {dict(us)}, delta {dict(ds)}; "
                  + "; ".join(f"{'-'.join(map(str, r_['cfg']))} {tuple(r_['pair'])} c1={tuple(r_['c'][0])} c2={r_['c'][1][0]} M={r_.get('M')} u={r_.get('u')}" for r_ in per[si]), flush=True)
        print(f"  surfaces by the values of delta over their pairs: {dict(dist)}", flush=True)




def validate(p, cfg, pair, tw=0, r=1, npts=3, verbose=True):
    res = sp.one_pair(p, cfg, pair)
    if res["status"] != "computed":
        print("pair not computed by sixlines_prym:", res["status"])
        return None
    t0 = time.time()
    A = Jac(p, cfg, pair, tw, r)
    t1 = time.time()
    c = res["c"]
    n = A.card(c)
    out = {"p": p, "cfg": cfg, "pair": list(pair), "tw": tw, "r": r, "order": n, "setup_s": round(t1 - t0, 2)}
    ok = Counter()
    for _ in range(npts):
        x = A.random()
        ok["x != 0"] += not A.is_zero(x)
        ok["order kills x"] += A.is_zero(A.mul(n, x))
        ok["(order + 1) x = x"] += A.eq(A.mul(n + 1, x), x)
        ok["J^4 x = x"] += A.eq(A.J(A.J(A.J(A.J(x)))), x)
        ok["Frob^r x = x"] += A.eq(A.frob(x, r), x)
        ok["J Frob = Frob J"] += A.eq(A.J(A.frob(x)), A.frob(A.J(x)))
        ok["key(x + 0) = key(x)"] += A.key(A.add(x, A.zero)) == A.key(x)
        y = A.sub(x, A.J(x, 2))          # (1 - J^2) x: in the part with J^2 = -1
        ok["(J^2 + 1)(1 - J^2) x = 0"] += A.is_zero(A.add(y, A.J(y, 2)))
        for s in (1, 3):
            ok[f"chi(Frob) kills (1 - J^2) x, i -> J^{s}"] += A.is_zero(A.poly(A.chi(c), y, s))
        w = A.add(A.add(x, A.J(x, 2)), A.add(A.J(x, 1), A.J(x, 3)))    # (1 + J + J^2 + J^3) x: the part from E
        a = p + 1 - A.cardE
        ok["(Frob^2 - a Frob + p) kills (1 + J + J^2 + J^3) x"] += A.is_zero(A.poly([(p, 0), (-a, 0), (1, 0)], w))
    out["checks"] = dict(ok)
    out["seconds"] = round(time.time() - t1, 2)
    if verbose:
        print(json.dumps(out))
    return out


def character_extra(workers, main_path, out_path, cases=None, r=4, mm=12):
    """Attempts to break the result of character(): (a) the surfaces with fewer than two usable pairs, on all their
    stored pairs; (b) the first usable pair of the first four surfaces of each torsor on the twists 1 and 2;
    (c) the same pairs with other random classes."""
    main = [json.loads(line) for line in open(main_path)]
    fh = open(out_path, "w") if out_path else None
    for p, e in CHARACTER_CASES:
        if cases and (p, e[0]) not in cases:
            continue
        rows = [x for x in main if x["p"] == p and x["e"] == e]
        ok = defaultdict(list)
        for x in rows:
            if x["status"] == "ok":
                ok[x["surface"]].append(x)
        surf = surfaces_of(p, e)
        sp_ = stored_pairs(p)
        done = {(tuple(x["cfg"]), tuple(x["pair"])) for x in rows}
        jobs, tag = [], []
        for si, cfgs in enumerate(surf):
            if len(ok[si]) < 2:
                for cfg in cfgs:
                    for pair, c, comp, rat in sp_[tuple(cfg)]:
                        if (tuple(cfg), tuple(pair)) not in done:
                            jobs.append((p, list(cfg), pair, c, 0, r, mm, 1))
                            tag.append(("more pairs", si))
        for si in sorted(k for k, v in ok.items() if v)[:4]:
            x = ok[si][0]
            for tw in (1, 2):
                jobs.append((p, x["cfg"], tuple(x["pair"]), x["c"], tw, r, mm, 1))
                tag.append((f"twist {tw}", si))
            jobs.append((p, x["cfg"], tuple(x["pair"]), x["c"], 0, r, mm, 7))
            tag.append(("other classes", si))
        with Pool(workers) as pool:
            recs = pool.map(_guarded_character, jobs, chunksize=1)
        ref = {si: {x["delta"] for x in v} for si, v in ok.items()}
        refu = {(tuple(x["cfg"]), tuple(x["pair"])): x["u"] for v in ok.values() for x in v}
        tally = Counter()
        newd = defaultdict(Counter)
        for (kind, si), rec in zip(tag, recs):
            rec.update({"surface": si, "e": e, "kind": kind})
            if fh:
                fh.write(json.dumps(rec) + "\n")
            if rec["status"] != "ok":
                tally[(kind, rec["status"][:40])] += 1
                continue
            if kind == "more pairs":
                newd[si][rec["delta"]] += 1
                tally[(kind, "agrees with the surface" if not ref.get(si) or ref[si] == {rec["delta"]} else "DIFFERS")] += 1
            else:
                tally[(kind, "same u" if refu[(tuple(rec["cfg"]), tuple(rec["pair"]))] == rec["u"] else "DIFFERENT u")] += 1
        if fh:
            fh.flush()
        print(f"p = {p}, e = {e}: {len(jobs)} runs: {dict(sorted(tally.items()))}")
        for si in sorted(newd):
            print(f"  surface {si} {['-'.join(map(str, c)) for c in surf[si]]}: before {sorted(ref.get(si, []))}, new pairs delta {dict(newd[si])}", flush=True)


def _validate_job(job):
    try:
        return validate(*job, npts=2, verbose=False)
    except Exception as err:          # a failure is a result to report, not to hide
        return {"p": job[0], "cfg": job[1], "pair": list(job[2]), "tw": job[3], "r": job[4], "error": str(err)}


def validate_set(workers, out_path):
    """The checks of validate() on a fixed set: for p = 13, 17, 29, 37 the first two configurations with irreducible
    P that have stored pairs, two pairs each, twists 0 and 1 over F_p; twist 0 over F_{p^2}; one curve per prime over
    F_{p^3}; one curve over F_{29^4}."""
    jobs = []
    for p in (13, 17, 29, 37):
        irr = {tuple(r["cfg"]) for r in map(json.loads, open(os.path.join(HERE, "sixlines_prym.jsonl")))
               if r["p"] == p and r.get("degrees") == [6]}
        st = stored_pairs(p)
        cfgs = sorted(c for c in st if c in irr)[:2]
        for ci, cfg in enumerate(cfgs):
            for pi, (pair, c, comp, rat) in enumerate(st[cfg][:2]):
                jobs += [(p, list(cfg), pair, tw, 1) for tw in (0, 1)]
                if pi == 0:
                    jobs.append((p, list(cfg), pair, 0, 2))
                if pi == 0 and ci == 0:
                    jobs.append((p, list(cfg), pair, 2, 3))
                if p == 29 and pi == 0 and ci == 0:
                    jobs.append((p, list(cfg), pair, 0, 4))
    with Pool(workers) as pool:
        recs = pool.map(_validate_job, jobs, chunksize=1)
    total, bad = Counter(), 0
    fh = open(out_path, "w") if out_path else None
    for rec in recs:
        if fh:
            fh.write(json.dumps(rec) + "\n")
        if "error" in rec:
            bad += 1
            print("ERROR", rec)
            continue
        for k, v in rec["checks"].items():
            total[k] += v
        total["classes"] += 2
        total["curves"] += 1
        total[f"curves over F_p^{rec['r']}"] += 1
    print(f"validate-set: {len(jobs)} runs, {bad} with an error")
    for k, v in total.items():
        print(f"  {k}: {v}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", nargs=7, type=int, metavar=("p", "b", "c", "e", "f", "i", "j"))
    ap.add_argument("--tw", type=int, default=0)
    ap.add_argument("--r", type=int, default=1)
    ap.add_argument("--validate-set", action="store_true")
    ap.add_argument("--levels", action="store_true")
    ap.add_argument("--lattice", action="store_true")
    ap.add_argument("--character", action="store_true")
    ap.add_argument("--character-extra", metavar="FILE", help="the output of --character; runs the controls of character_extra()")
    ap.add_argument("--torsor", nargs=2, type=int, action="append", metavar=("p", "e1"))
    ap.add_argument("--one", nargs=7, type=int, metavar=("p", "b", "c", "e", "f", "i", "j"), help="step 3 on one configuration and pair")
    ap.add_argument("--case", nargs=3, type=int, action="append", metavar=("p", "e1", "l"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pairs", type=int, default=8)
    ap.add_argument("--out")
    args = ap.parse_args()
    if args.validate_set:
        validate_set(args.workers, args.out)
    if args.levels:
        levels(args.workers, args.out, [tuple(x) for x in args.case] if args.case else None, args.pairs)
    if args.lattice:
        lattice_side()
    if args.character:
        character(args.workers, args.out, [tuple(x) for x in args.torsor] if args.torsor else None, args.pairs, args.tw, args.r if args.r > 1 else 4)
    if args.character_extra:
        character_extra(args.workers, args.character_extra, args.out, [tuple(x) for x in args.torsor] if args.torsor else None)
    if args.one:
        v = args.one
        res = sp.one_pair(v[0], v[1:5], tuple(v[5:7]))
        print(json.dumps(character_pair((v[0], v[1:5], tuple(v[5:7]), res["c"], args.tw, args.r if args.r > 1 else 4, 12, 1))))
    if args.validate:
        v = args.validate
        validate(v[0], v[1:5], tuple(v[5:7]), args.tw, args.r)


if __name__ == "__main__":
    main()
