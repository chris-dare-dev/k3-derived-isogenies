#!/usr/bin/env python
"""
sixlines_cmmethod.py -- the CM method for six-line double planes, and the action of the class group (note 20).

A six-line double plane w^2 = l1 ... l6 over C has a period point tau in the hermitian upper half space

    D = { tau in M_2(C) : (tau - tau^*) / 2i > 0 },

and the inverse of the period map is given by theta constants (Matsumoto 1993; formulas as stated in
Matsumoto-Terasoma, arXiv:1001.4865, section 3.4, and Clingher-Malmendier-Shaska, arXiv:1806.07460, section 2):

    Theta[a b](tau) = sum over n in Z[i]^2 of e[ (1/2) m tau m^* + Re(m b'^*) ],   m = n + a',
                      a' = a / (1 + i),  b' = b / (1 + i),  a, b in {0, 1}^2,  e[x] = exp(2 pi i x).

The ten characteristics with a b^* = 0 mod (1 + i) are indexed by the splittings of {1..6} into two triples
(table 1 of Matsumoto-Terasoma), and

    [ Theta_<J>(tau)^2 ]_J  =  [ D(J) D(J^c) ]_J   in P^9,

D(J) the 3 x 3 minor of the 3 x 6 matrix of the lines.  The squares are modular forms of weight 2 for the
congruence subgroup of level 1 + i of U(2,2; Z[i]), extended by tau -> transpose(tau); the quotient of the
full group by it is S_6 (relabelling the lines), and the transpose is association.  The polynomials
J2, J3, J4, J5, J6 in the ten squares (Clingher-Malmendier-Shaska (2.30)-(2.34)) are invariant under S_6, and
[J2 : J3 : J4 : J5 : J6] in P(2,3,4,5,6) is the moduli point of the six lines up to relabelling and association,
that is of the K3 surface.

The lattices.  K = Q(i)[X]/chi, chi = X^4 - 4(nu + mu i)X^3 + e1 X^2 - 4q(nu - mu i)X + q^2 (note 13, Proposition Y),
of degree 8 over Q, with embeddings phi_1..phi_4 over Q(i).  A type is a pair T = {j, k}.  A polarised lattice is
(L, xi): L a Z[i]-lattice in K, xi totally imaginary with xi L = conj(trace dual of L), Im phi(xi) < 0 exactly for
phi_j, phi_k.  The form 2 Tr_{K/Q(i)}(xi x conj(y)) is unimodular skew-hermitian over Z[i]; when it is even the
lattice has a basis b_1..b_4 with Gram matrix [[0, -1], [1, 0]], and the period point is

    tau = conj(B1^-1 B2),   B1 = (phi_r(b_k)) for r not in T and k = 1, 2,   B2 the same for k = 3, 4.

The surfaces with six rational lines over F_q are the even lattices stable under the order
R = Z[i, pi, q/pi, sigma, conj(sigma)], sigma = (pi - 1)/(1 + i); the lattices of the maximal order are odd.

  selftest                 checks that do not involve CM: the linear relations among the ten squares, the
                           S_6-invariance of J2..J6 on matrices, the transformation of J2..J6 under generators of
                           the full modular group, and a curve of genus 2 (six lines tangent to a conic) through its
                           numerically computed period matrix
  points q e1 mu nu        the polarised lattices of the field, their parity, period points and invariants
                           (--maximal: the maximal order only; --nolevel: the order Z[i, pi, q/pi])
  compare q e1 mu nu       class polynomials over the reflex field E from the even lattices, reduction at the prime
                           (q, F - e1) of E, comparison with the surfaces of sixlines_cm.jsonl
  action q e1 mu nu        [4a] the permutations of the lattices induced by the ideals above primes l that split
                           in K, against the classes in G = Cl(E)/H of the primes of E (note 13, Lemma 9; note 15
                           section 3.3 with the correction of note 18);  [4b] x -> x^l on the roots of the class
                           polynomial in characteristic l against these permutations;  [4c] labels in G for the
                           surfaces over F_q
  sweep p [p ...]          compare and action for every polynomial of sixlines_cm.jsonl with maximal order at the
                           given primes, raising the precision until the class polynomials are recognised

Uses cypari2 and numpy.  PARI's class groups are not certified (GRH).  Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import os
import re
import sys
import time
from fractions import Fraction

import cypari2


pari = cypari2.Pari()
pari.allocatemem(4 * 10 ** 9, silent=True)
pari("'x; 'y; 'z")

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ GP side

GPSRC = r"""
\\ ---------- theta constants on D ----------
\\ One-dimensional sums  sum_{x in Z} s^x exp(pi i t X^2 + al X),  X = x + sh,  for s = 1 and s = -1.
cm_th1(t, al, sh, R2) = {
  my(yt = imag(t), c = real(al) / (2 * Pi * yt), w = sqrt(R2 / yt) + 2, x0 = round(c - sh), X, v, r, e1, e2,
     s0 = 0, s1 = 0, q2 = exp(2 * Pi * I * t), sg);
  \\ upwards from x0
  X = x0 + sh; v = exp(Pi * I * t * X^2 + al * X); r = exp(Pi * I * t * (2 * X + 1) + al);
  sg = if(x0 % 2, -1, 1);
  for(k = 0, ceil(w),
    s0 += v; s1 += sg * v; v *= r; r *= q2; sg = -sg);
  \\ downwards from x0 - 1
  X = x0 - 1 + sh; v = exp(Pi * I * t * X^2 + al * X); r = exp(Pi * I * t * (-2 * X + 1) - al);
  sg = if((x0 - 1) % 2, -1, 1);
  for(k = 0, ceil(w),
    s0 += v; s1 += sg * v; v *= r; r *= q2; sg = -sg);
  [s0, s1];
}

\\ All sixteen Theta[a1 a2 b1 b2](tau); returns a vector indexed by 8 a1 + 4 a2 + 2 b1 + b2 + 1.
\\ With quarter = 1 the sum over m1 runs over a quarter plane (the terms are invariant under m -> i m for the ten
\\ even characteristics); the six other entries are then meaningless.
cm_theta16(tau, R2, quarter = 0) = {
  my(t11 = tau[1, 1], t12 = tau[1, 2], t21 = tau[2, 1], t22 = tau[2, 2], Y, Yi, B1, res = vector(16), a1, a2,
     X1, Y1, m1, al, be, FX, FY, c, sg, bx, by, idx, cx, cy, wt, bd);
  Y = (tau - conj(tau~)) / (2 * I);
  Yi = Y^(-1);
  bd = R2 * real(Yi[1, 1]) + 1;
  B1 = ceil(sqrt(bd)) + 1;
  for(a1 = 0, 1, for(a2 = 0, 1,
    for(x1 = if(quarter, 0, -B1 - 1), B1 + 1, for(y1 = if(quarter, 0, -B1 - 1), B1 + 1,
      X1 = x1 + a1 / 2; Y1 = y1 - a1 / 2;
      wt = 1;
      if(quarter,
        if(a1, Y1 = y1 + 1 / 2);
        if(X1 == 0 && Y1 == 0, wt = 1, if(X1 == 0, next, wt = 4)));
      if((X1^2 + Y1^2) > bd, next);
      m1 = X1 + I * Y1;
      al = Pi * I * (t12 * m1 + t21 * conj(m1));
      be = Pi * (t12 * m1 - t21 * conj(m1));
      FX = cm_th1(t22, al, a2 / 2, R2);
      FY = cm_th1(t22, be, -a2 / 2, R2);
      c = wt * exp(Pi * I * t11 * (X1^2 + Y1^2));
      sg = if((round(X1 - a1 / 2) - round(Y1 + a1 / 2) + a1) % 2, -1, 1);
      for(b1 = 0, 1, for(b2 = 0, 1,
        idx = 8 * a1 + 4 * a2 + 2 * b1 + b2 + 1;
        res[idx] += if(b1, sg, 1) * if(b2 && a2, -1, 1) * c * FX[b2 + 1] * FY[b2 + 1]));
    ));
  ));
  res;
}

\\ The ten squares t1..t10 in the order of Clingher-Malmendier-Shaska (2.1):
\\ t1 = <135>, t2 = <145>, t3 = <235>, t4 = <245>, t5 = <125>, t6 = <345>, t7 = <134>, t8 = <124>, t9 = <234>, t10 = <123>
\\ with the characteristics [a1 a2 b1 b2] of Matsumoto-Terasoma, table 1.
cm_tsq(tau, R2) = {
  my(th = cm_theta16(tau, R2, 1),
     ch = [[0,0,0,0], [1,1,0,0], [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,1], [0,0,1,1], [1,0,0,1], [1,1,1,1]]);
  vector(10, k, th[8 * ch[k][1] + 4 * ch[k][2] + 2 * ch[k][3] + ch[k][4] + 1]^2);
}

\\ J2..J6 from t1..t10  (Clingher-Malmendier-Shaska (2.34), (2.32))
cm_J(t) = {
  my(x = vector(6), s);
  x[1] = 2 * t[1] + 2 * t[5] - 3 * t[6] - t[7] - t[8];
  x[2] = -t[1] - t[5] - t[7] + 2 * t[8];
  x[3] = -t[1] + 2 * t[5] - t[7] - t[8];
  x[4] = -t[1] - t[5] + 3 * t[6] + 2 * t[7] + 2 * t[8];
  x[5] = -t[1] - t[5] + 2 * t[7] - t[8];
  x[6] = 2 * t[1] - t[5] - t[7] - t[8];
  s = vector(6, j, sum(k = 1, 6, x[k]^j));
  [s[2] / 12, s[3] / 12, (4 * s[4] - s[2]^2) / 64, (5 * s[2] * s[3] - 12 * s[5]) / 240,
   (3 * s[2]^3 - 4 * s[3]^2 - 18 * s[2] * s[4] + 24 * s[6]) / 576];
}

\\ the ten products of minors of a 3 x 6 matrix, same order
cm_tmat(A) = {
  my(D = (i, j, k) -> matdet(Mat([A[, i], A[, j], A[, k]])));
  [D(1,3,5) * D(2,4,6), D(1,4,5) * D(2,3,6), D(1,4,6) * D(2,3,5), D(1,3,6) * D(2,4,5), D(1,2,5) * D(3,4,6),
   D(1,2,6) * D(3,4,5), D(1,3,4) * D(2,5,6), D(1,2,4) * D(3,5,6), D(1,5,6) * D(2,3,4), D(1,2,3) * D(4,5,6)];
}

\\ the fifteen linear relations of (2.2); each row is a vector of coefficients on t1..t10
cm_rel() = {
  my(L = [[1,-2,-5,-9], [1,-2,-6,-7], [1,-3,-5,-10], [1,-3,-6,-8], [1,-4,-7,-10], [1,-4,-8,-9], [2,-3,7,-8],
          [2,-3,9,-10], [2,-4,5,-8], [2,-4,6,-10], [3,-4,5,-7], [3,-4,6,-9], [5,-6,-7,9], [5,-6,-8,10], [7,-8,-9,10]]);
  vector(#L, r, my(v = vector(10)); for(k = 1, 4, v[abs(L[r][k])] = sign(L[r][k])); v);
}

\\ action of g in U(2,2) on tau, and the automorphy factor det(C tau + D)
cm_act(g, tau) = {
  my(A = g[1..2, 1..2], B = g[1..2, 3..4], C = g[3..4, 1..2], Dd = g[3..4, 3..4]);
  [(A * tau + B) * (C * tau + Dd)^(-1), matdet(C * tau + Dd)];
}
"""


def load_gp(src):
    """One pari() call per function; comment lines are removed (the library parser does not skip them)."""
    lines = [l for l in src.split("\n") if not l.lstrip().startswith("\\\\")]
    lines = [re.sub(r"\s*\\\\.*$", "", l) for l in lines]
    cur = []
    for l in lines:
        cur.append(l)
        if l == "}":
            pari("\n".join(cur))
            cur = []
    assert not "".join(cur).strip(), cur


load_gp(GPSRC)


def setprec(digits):
    pari.default("realprecision", digits)


def R2_for(digits):
    """n Y n^* <= R2 keeps the terms above 10^-(digits + 10)."""
    return (digits + 10) * 2.302585093 / 3.141592653589793


# ------------------------------------------------------------------------------------------------ selftest

def random_tau(seed):
    """A point of D that is not symmetric, entries of moderate size."""
    pari.setrand(seed)
    return pari(r"""
      my(Y, X, v);
      v = vector(8, k, random(1.) - 1/2);
      Y = [1.1 + v[1]/3, (v[2] + I*v[3])/3; (v[2] - I*v[3])/3, 0.9 + v[4]/3];
      X = [v[5], v[6] + I*v[7]; v[6] - I*v[7], v[8]];
      X + I*Y
    """)


GENERATORS = r"""
[ [I,0,0,0; 0,1,0,0; 0,0,I,0; 0,0,0,1],
  [1,0,1,0; 0,1,0,0; 0,0,1,0; 0,0,0,1],
  [1,1,0,0; 0,1,0,0; 0,0,1,0; 0,0,-1,1],
  [0,1,0,0; 1,0,0,0; 0,0,0,1; 0,0,1,0],
  [0,0,1,0; 0,0,0,1; -1,0,0,0; 0,-1,0,0] ]
"""


def selftest(digits):
    setprec(digits)
    R2 = R2_for(digits)
    ok = True
    tol = pari(f"10.^(-{digits - 15})")

    # 1. S_6-invariance of J2..J6 and the relations (2.2) on matrices, exactly
    pari.setrand(7)
    A = pari("matrix(3, 6, i, j, random(41) - 20)")
    t = pari("cm_tmat")(A)
    rel = pari("cm_rel")()
    bad = [r for r in range(15) if pari("(v, t) -> v * t~")(rel[r], t) != 0]
    print(f"relations (2.2) on a random integer matrix: {15 - len(bad)} of 15 hold")
    ok &= not bad
    J0 = pari("cm_J")(t)
    ninv = 0
    for perm in itertools.permutations(range(1, 7)):
        Ap = pari("(A, p) -> matrix(3, 6, i, j, A[i, p[j]])")(A, list(perm))
        ninv += pari("cm_J")(pari("cm_tmat")(Ap)) == J0
    print(f"J2..J6 unchanged under relabelling the lines: {ninv} of 720 permutations")
    ok &= ninv == 720
    # association: the Gale transform of (I | M) is (I | M^-t) up to scaling the columns
    Aa = pari("A -> my(g = A[, 1..3]^(-1) * A); concat(matid(3), (g[, 4..6]~)^(-1))")(A)
    Ja = pari("cm_J")(pari("cm_tmat")(Aa))
    w = [2, 3, 4, 5, 6]
    lam = None
    same = True
    # weighted projective equality: Ja[k] = c^w[k] J0[k]
    r2 = Ja[0] / J0[0]
    for k in range(5):
        for l in range(5):
            if Ja[k] ** w[l] * J0[l] ** w[k] != Ja[l] ** w[k] * J0[k] ** w[l]:
                same = False
    print(f"[J2:...:J6] unchanged under association (weighted projective, up to sign conventions): {same}")
    ok &= same

    # 2. relations (2.2) among the theta squares at three random points
    for seed in (1, 2, 3):
        tau = random_tau(seed)
        t = pari("cm_tsq")(tau, R2)
        scale = max(abs(complex(x)) for x in t)
        worst = max(abs(complex(pari("(v, t) -> v * t~")(rel[r], t))) for r in range(15)) / scale
        th = pari("cm_theta16")(tau, R2)
        odd = [abs(complex(th[k])) for k in range(16)
               if k not in (0, 12, 8, 4, 2, 6, 1, 3, 9, 15)]
        print(f"tau #{seed}: relations (2.2) on theta squares, largest defect {worst:.1e}; "
              f"largest of the six other theta constants {max(odd):.1e}")
        ok &= worst < 1e-12 and max(odd) < 1e-12

    # 3. transformation of J_k under generators of the full group
    tau = random_tau(5)
    J = pari("cm_J")(pari("cm_tsq")(tau, R2))
    gens = pari(GENERATORS)
    for gi in range(5):
        g = gens[gi]
        chk = pari("g -> my(N = [0,0,-1,0; 0,0,0,-1; 1,0,0,0; 0,1,0,0]); g * N * conj(g~) == N")(g)
        tau2, f = pari("cm_act")(g, tau)
        J2 = pari("cm_J")(pari("cm_tsq")(tau2, R2))
        dg = pari.matdet(g)
        worst = 0.0
        for k in range(5):
            wk = k + 2
            pred = dg ** wk * f ** (2 * wk) * J[k]
            worst = max(worst, abs(complex((J2[k] - pred) / pred)))
        print(f"generator G{gi + 1} (unitary: {bool(chk)}, det {dg}): J_k(g tau) = det(g)^k det(C tau + D)^2k J_k(tau), "
              f"largest relative defect {worst:.1e}")
        ok &= bool(chk) and worst < 1e-12
    tauT = pari("t -> t~")(tau)
    JT = pari("cm_J")(pari("cm_tsq")(tauT, R2))
    worst = max(abs(complex((JT[k] - J[k]) / J[k])) for k in range(5))
    print(f"transpose: J_k(tau^t) = J_k(tau), largest relative defect {worst:.1e}")
    ok &= worst < 1e-12

    # 4. genus 2: six lines tangent to a conic, period matrix by numerical integration
    ok &= genus2_check(digits)
    print("selftest:", "passed" if ok else "FAILED")
    return ok


def genus2_check(digits):
    """y^2 = prod (x - p_i) with real p_1 < ... < p_6.  The lines (p^2, -p, 1) are tangent to a conic.

    Periods of dx/y and x dx/y over the cycles a1 = around [p1, p2], a2 = around [p3, p4], b1 = around [p2, p3] u
    ... are computed from the five real integrals between consecutive branch points; a symplectic basis is
    a1 = I1, a2 = I3, b1 = I2 + I4, b2 = I4 up to signs, and any error in it shows as a period matrix that is not
    symmetric with positive imaginary part.  The comparison is between weighted projective points, so the choice
    of basis and of labelling does not matter.
    """
    setprec(digits)
    R2 = R2_for(digits)
    res = pari(r"""
    (pts) -> {
      my(f = x -> prod(k = 1, 6, x - pts[k]), W, P, best = 0);
      \\ W[j][k] = int_{p_k}^{p_{k+1}} x^(j-1) dx / sqrt|f|, k = 1..5
      W = matrix(2, 5, j, k, intnum(x = [pts[k], -1/2], [pts[k + 1], -1/2], x^(j - 1) / sqrt(abs(f(x)))));
      W;
    }
    """)([-3, -2, -1, 1, 2, 3])
    # On [p_k, p_{k+1}] the sign of f is (-1)^(6-k): f < 0 on (p1,p2), > 0 on (p2,p3), < 0 on (p3,p4), ...
    # so the integrals over (p1,p2), (p3,p4), (p5,p6) are imaginary periods and the others real, up to factors 2.
    # Try the sign patterns of a standard homology basis and keep those giving a symmetric tau in the Siegel space.
    found = None
    cands = []
    for sa in itertools.product((1, -1), repeat=4):
        cand = pari(r"""
        (W, s) -> {
          my(A, B, tau);
          \\ a1 = 2 i W1, a2 = 2 i W3 ; b1 = 2 (W2 + W4)-type combinations
          A = Mat([2 * I * W[, 1], 2 * I * s[1] * W[, 3]]);
          B = Mat([2 * s[2] * W[, 2] + 2 * s[3] * W[, 4], 2 * s[4] * W[, 4]]);
          tau = A^(-1) * B;
          tau;
        }
        """)(res, list(sa))
        sym = abs(complex(cand[0, 1] - cand[1, 0]))
        y = pari("t -> imag(t)")(cand)
        pos = float(y[0, 0]) > 0 and float(pari.matdet(y)) > 0
        if sym < 1e-20 and pos:
            cands.append((sa, cand))
    if not cands:
        # other orientation: periods in the other order
        for sa in itertools.product((1, -1), repeat=4):
            cand = pari(r"""
            (W, s) -> {
              my(A, B, tau);
              A = Mat([2 * W[, 2] + 2 * s[3] * W[, 4], 2 * s[4] * W[, 4]]);
              B = Mat([2 * I * s[1] * W[, 1], 2 * I * s[2] * W[, 3]]);
              tau = A^(-1) * B;
              tau;
            }
            """)(res, list(sa))
            sym = abs(complex(cand[0, 1] - cand[1, 0]))
            y = pari("t -> imag(t)")(cand)
            pos = float(y[0, 0]) > 0 and float(pari.matdet(y)) > 0
            if sym < 1e-20 and pos:
                cands.append((sa, cand))
    print(f"genus 2: {len(cands)} sign patterns give a symmetric period matrix in the Siegel space")
    if not cands:
        return False
    A = pari("p -> matrix(3, 6, i, j, [p[j]^2, -p[j], 1][i])")([-3, -2, -1, 1, 2, 3])
    Jl = pari("cm_J")(pari("cm_tmat")(A))
    good = 0
    for sa, tau in cands:
        Jt = pari("cm_J")(pari("cm_tsq")(tau, R2))
        d = wp_distance(Jt, Jl)
        print(f"  signs {sa}: distance between [J(theta)] and [J(lines)] in P(2,3,4,5,6): {d:.1e};"
              f" J4(theta)/J2^2 = {abs(complex(Jt[2] / Jt[0] ** 2)):.1e}")
        good += d < 1e-12
    return good == len(cands)


def wp_distance(Ja, Jb):
    """0 when [Ja] = [Jb] in P(2,3,4,5,6); uses the degree-0 ratios with J2 and J3."""
    a, b = invariants(Ja), invariants(Jb)
    worst = 0.0
    for x, y in zip(a, b):
        size = float(abs(x) + abs(y))
        if size > 1e-25:
            worst = max(worst, float(abs(x - y)) / size)
    return worst


# ------------------------------------------------------------------------------------------------ CM points

GPCM = r"""
cm_field(q, e1, mu, nu) = {
  my(a = 4 * (nu + mu * 'y), ab = 4 * (nu - mu * 'y), chi, W, nf, bnf, iK, cK, rts, emb, pl, mi);
  chi = 'z^4 - a * 'z^3 + e1 * 'z^2 - q * ab * 'z + q^2;
  W = polresultant(chi, 'y^2 + 1, 'y);
  if(!polisirreducible(W), return(0));
  bnf = bnfinit(W, 1); nf = bnf.nf;
  iK = nfroots(nf, 'x^2 + 1);
  iK = if(subst(subst(chi, 'y, iK[1]), 'z, Mod('z, W)) == 0, iK[1], iK[2]);
  if(subst(subst(chi, 'y, iK), 'z, Mod('z, W)) != 0, error("no i"));
  cK = Mod(q, W) / Mod('z, W);
  if(subst(lift(iK), 'z, cK) != -iK, error("conjugation"));
  [bnf, iK, cK, chi, W];
}

cm_embeddings(F) = {
  my(W = F[5], iK = lift(F[2]), rts = polroots(W), pl = List());
  for(k = 1, #rts, if(imag(subst(iK, 'z, rts[k])) > 0, listput(pl, rts[k])));
  if(#pl != 4, error("embeddings"));
  vecsort(Vec(pl), (u, v) -> if(abs(real(u) - real(v)) > 1e-20, sign(real(u) - real(v)), sign(imag(u) - imag(v))));
}

cm_conj(F, z) = {
  subst(lift(z), 'z, F[3]);
}

cm_classreps(F) = {
  my(bnf = F[1], cyc = bnf.cyc, gen = bnf.gen, L = List(), Id);
  forvec(e = vector(#cyc, k, [0, cyc[k] - 1]),
    Id = idealhnf(bnf, 1);
    for(k = 1, #cyc, Id = idealmul(bnf, Id, idealpow(bnf, gen[k], e[k])));
    Id = idealred(bnf, Id);
    listput(L, [Vec(e), Id]));
  Vec(L);
}

cm_polarisations(F, Id) = {
  my(bnf = F[1], nf = bnf.nf, W = F[5], Ib, C, v, xi0, tu, fu, w, out = List(), u, xi);
  Ib = nfgaloisapply(nf, lift(F[3]), Id);
  C = idealinv(nf, idealmul(nf, idealmul(nf, Id, Ib), nf.diff));
  v = bnfisprincipal(bnf, C);
  if(v[1] != 0, return([]));
  xi0 = nfbasistoalg(nf, v[2]);
  tu = Mod(lift(bnf.tu[2]), W); w = bnf.tu[1]; fu = bnf.fu;
  for(a = 0, w - 1, forvec(e = vector(#fu, k, [0, 1]),
    u = tu^a * prod(k = 1, #fu, Mod(lift(fu[k]), W)^e[k]);
    xi = u * xi0;
    if(cm_conj(F, xi) == -xi, listput(out, xi))));
  Vec(out);
}

cm_isnorm(F, r) = {
  my(bnf = F[1], nf = bnf.nf, W = F[5], tu = Mod(lift(bnf.tu[2]), W), w = bnf.tu[1], y);
  for(a = 0, w - 1,
    if(nfeltissquare(nf, r * tu^a, &y),
      y = Mod(lift(nfbasistoalg(nf, y)), W);
      if(y * cm_conj(F, y) == r, return(1))));
  0;
}

cm_qicoords(F, Id) = {
  my(nf = F[1].nf, W = F[5], iK = lift(F[2]), M, zk = nf.zk, H = idealhnf(nf, Id), out = List(), al, c);
  M = matrix(8, 8, r, k, polcoef(if(k <= 4, 'z^(k - 1), lift(Mod(iK * 'z^(k - 5), W))), r - 1, 'z));
  for(m = 1, 8,
    al = lift(Mod(zk * H[, m], W));
    c = matsolve(M, vectorv(8, r, polcoef(al, r - 1, 'z)));
    listput(out, vector(4, k, c[k] + I * c[k + 4])));
  Vec(out);
}

cm_elt(F, v) = {
  my(W = F[5], iK = lift(F[2]));
  Mod(sum(k = 1, 4, (real(v[k]) + imag(v[k]) * iK) * 'z^(k - 1)), W);
}

cm_form(F, xi, bj, bk) = {
  my(z = xi * bj * cm_conj(F, bk));
  trace(z) - I * trace(F[2] * z);
}

cm_qf8(G) = {
  my(n = #G, M = matrix(2 * n, 2 * n));
  for(a = 1, n, for(b = 1, n,
    M[a, b] = real(I * G[a, b]); M[a, n + b] = real(G[a, b]);
    M[n + a, b] = -real(G[a, b]); M[n + a, n + b] = real(I * G[a, b])));
  M;
}

cm_tau(F, pl, S, B4) = {
  my(Sc = setminus([1, 2, 3, 4], S), Bm);
  Bm = matrix(4, 4, j, k, subst(lift(B4[k]), 'z, pl[j]));
  conj(matrix(2, 2, r, k, Bm[Sc[r], k])^(-1) * matrix(2, 2, r, k, Bm[Sc[r], k + 2]));
}

cm_reduce(tau) = {
  my(Y, U, H, Bh, t = tau, n = 0, ch = 1, k, c);
  while(ch && n < 200,
    ch = 0; n++;
    Y = (t - conj(t~)) / (2 * I);
    if(real(Y[1, 1]) > real(Y[2, 2]) * (1 + 1e-9), U = [0, 1; 1, 0]; t = U * t * conj(U~); ch = 1; next);
    c = round(Y[2, 1] / Y[1, 1]);
    if(c != 0, U = [1, 0; -c, 1]; t = U * t * conj(U~); ch = 1; next);
    H = (t + conj(t~)) / 2;
    Bh = [round(real(H[1, 1])), round(H[1, 2]); conj(round(H[1, 2])), round(real(H[2, 2]))];
    if(Bh != 0, t = t - Bh);
    if(norm(t[1, 1]) < 0.98,
      t = ([0, 0; 0, 1] * t + [-1, 0; 0, 0]) * ([1, 0; 0, 0] * t + [0, 0; 0, 1])^(-1); ch = 1; next);
    if(norm(matdet(t)) < 0.98, t = -t^(-1); ch = 1; next));
  t;
}
"""

load_gp(GPCM)


def G(name):
    return pari(name)


def gauss_div(a, b):
    """Nearest Gaussian integer to a / b."""
    return pari("(a, b) -> round(a / b)")(a, b)


def gauss_rowreduce(rows):
    """Row echelon form over Z[i] of a list of rows of Gaussian integers; returns the non-zero rows."""
    rows = [list(r) for r in rows]
    ncol = len(rows[0])
    out = []
    for col in range(ncol):
        while True:
            nz = [r for r in rows if r[col] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: int(pari.norm(r[col])))
            piv = nz[0]
            for r in nz[1:]:
                c = gauss_div(r[col], piv[col])
                for k in range(ncol):
                    r[k] = r[k] - c * piv[k]
        nz = [r for r in rows if r[col] != 0]
        if nz:
            out.append(nz[0])
            rows = [r for r in rows if r is not nz[0]]
    assert all(all(x == 0 for x in r) for r in rows)
    return out


def gauss_bezout(w):
    """x with sum x_m w_m = gcd, by row operations on (w | identity)."""
    n = len(w)
    rows = [[w[m]] + [pari(1) if k == m else pari(0) for k in range(n)] for m in range(n)]
    while True:
        nz = [r for r in rows if r[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda r: int(pari.norm(r[0])))
        piv = nz[0]
        for r in nz[1:]:
            c = gauss_div(r[0], piv[0])
            for k in range(n + 1):
                r[k] = r[k] - c * piv[k]
    piv = [r for r in rows if r[0] != 0][0]
    return piv[0], piv[1:]


def herm(G_, x, y):
    """x G y^* for rows x, y."""
    n = len(x)
    s = pari(0)
    for a in range(n):
        if x[a] == 0:
            continue
        for b in range(n):
            if y[b] != 0:
                s = s + x[a] * G_[a][b] * pari.conj(y[b])
    return s


def hyperbolic_pair(Gm):
    """For a unimodular even skew-hermitian Gram matrix Gm (list of lists) of rank n over Z[i]: rows e, f with
    s(e, e) = s(f, f) = 0, s(f, e) = 1, and a basis of the orthogonal complement.  None if no isotropic vector."""
    n = len(Gm)
    M8 = G("cm_qf8")(pari.matrix(n, n, [x for r in Gm for x in r]))
    sol = pari.qfsolve(M8)
    if pari.type(sol) == "t_INT":
        return None
    if pari.type(sol) == "t_MAT":
        sol = sol[0]
    e = [sol[a] + pari("I") * sol[n + a] for a in range(n)]
    g, _ = gauss_bezout(e)
    e = [x / g for x in e]
    assert herm(Gm, e, e) == 0
    w = [sum((Gm[m][k] * pari.conj(e[k]) for k in range(n)), pari(0)) for m in range(n)]   # s(b_m, e)
    g, x = gauss_bezout(w)
    assert pari.norm(g) == 1, "not unimodular"
    f = [xx / g for xx in x]
    assert herm(Gm, f, e) == 1
    sff = herm(Gm, f, f)
    assert pari.real(sff) == 0 and pari.imag(sff) % 2 == 0, "odd"
    c = pari("I") * (pari.imag(sff) / 2)
    f = [f[k] + c * e[k] for k in range(n)]
    assert herm(Gm, f, f) == 0 and herm(Gm, f, e) == 1 and herm(Gm, e, f) == -1
    comp = []
    for m in range(n):
        x = [pari(1) if k == m else pari(0) for k in range(n)]
        al = herm(Gm, x, f)
        be = -herm(Gm, x, e)
        comp.append([x[k] + al * e[k] + be * f[k] for k in range(n)])
    comp = gauss_rowreduce(comp)
    assert len(comp) == n - 2
    return e, f, comp


def standard_basis(Gm):
    """T (rows) with T G T^* = [[0, -1], [1, 0]] (blocks), for a 4 x 4 unimodular even skew-hermitian Gram matrix."""
    r = hyperbolic_pair(Gm)
    if r is None:
        return None
    e1, f1, comp = r
    G2 = [[herm(Gm, comp[a], comp[b]) for b in range(2)] for a in range(2)]
    r2 = hyperbolic_pair(G2)
    if r2 is None:
        return None
    e2c, f2c, rest = r2
    e2 = [e2c[0] * comp[0][k] + e2c[1] * comp[1][k] for k in range(4)]
    f2 = [f2c[0] * comp[0][k] + f2c[1] * comp[1][k] for k in range(4)]
    T = [e1, e2, f1, f2]
    chk = [[herm(Gm, T[a], T[b]) for b in range(4)] for a in range(4)]
    want = [[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]]
    assert all(chk[a][b] == want[a][b] for a in range(4) for b in range(4)), chk
    return T


# ------------------------------------------------------------------------------------------------ lattices of Z[i, pi]

GPLAT = r"""
lat_hnf(M) = {
  my(d = denominator(M));
  mathnf(d * M) / d;
}

lat_mulmat(nf, x) = {
  matconcat(vector(8, j, nfeltmul(nf, x, matid(8)[, j])));
}

lat_setup(F, q, lev) = {
  my(bnf = F[1], nf = bnf.nf, W = F[5], iK = lift(F[2]), gens, HR, H2, S, Mi, Mp, Cm, T, zk = nf.zk, sg);
  gens = [Mod(iK, W), Mod('z, W), Mod(q, W) / Mod('z, W)];
  if(lev,
    sg = (Mod('z, W) - 1) / (1 + Mod(iK, W));
    gens = concat(gens, [sg, Mod(cm_conj(F, sg), W)]));
  gens = apply(g -> lat_mulmat(nf, nfalgtobasis(nf, g)), gens);
  HR = matid(8)[, 1..1];
  while(1,
    H2 = lat_hnf(matconcat(concat([HR], vector(#gens, k, gens[k] * HR))));
    if(H2 == HR, break);
    HR = H2);
  if(denominator(HR) != 1, error("order not integral"));
  S = matsnf(HR, 1);
  Mi = lat_mulmat(nf, nfalgtobasis(nf, Mod(iK, W)));
  Mp = gens[2];
  Cm = matconcat(vector(8, j, nfalgtobasis(nf, Mod(cm_conj(F, Mod(zk[j], W)), W))));
  T = matrix(8, 8, j, k, nfelttrace(nf, nfeltmul(nf, matid(8)[, j], matid(8)[, k])));
  [HR, S, Mi, Mp, Cm, T, if(lev, gens[4], 0)];
}

lat_dualconj(LS, L) = {
  lat_hnf(LS[5] * (L~ * LS[6])^(-1));
}

lat_prod(nf, A, B) = {
  my(v = List());
  for(j = 1, #A, for(k = 1, #B, listput(v, nfeltmul(nf, A[, j], B[, k]))));
  lat_hnf(matconcat(Vec(v)));
}

lat_mod8(nf, x, L) = {
  mathnf(concat(lift(Mod(lat_mulmat(nf, x), 8)) * L, 8 * matid(8)));
}

lat_ideal(nf, L) = {
  my(Id = idealhnf(nf, L[, 1]));
  for(j = 2, #L, Id = idealadd(nf, Id, L[, j]));
  Id;
}

lat_units8(nf) = {
  my(pr = idealprimedec(nf, 2), rad, out = List(), B, c);
  rad = idealhnf(nf, 1);
  for(k = 1, #pr, rad = idealmul(nf, rad, pr[k]));
  for(j = 1, 24,
    B = idealhnf(nf, idealpow(nf, rad, j));
    for(k = 1, 8, listput(out, (matid(8)[, 1] + B[, k]) % 8)));
  if(#pr > 1 || pr[1].f > 1,
    forvec(v = vector(8, k, [0, 1]),
      c = vectorv(8, k, v[k]);
      if(c != 0 && nfeltnorm(nf, c) % 2 == 1, listput(out, c))));
  Set(Vec(out));
}

lat_qi(F, L) = {
  my(nf = F[1].nf, W = F[5], iK = lift(F[2]), M, zk = nf.zk, out = List(), al, c);
  M = matrix(8, 8, r, k, polcoef(if(k <= 4, 'z^(k - 1), lift(Mod(iK * 'z^(k - 5), W))), r - 1, 'z));
  for(m = 1, 8,
    al = lift(Mod(zk * L[, m], W));
    c = matsolve(M, vectorv(8, r, polcoef(al, r - 1, 'z)));
    listput(out, vector(4, k, c[k] + I * c[k + 4])));
  Vec(out);
}

lat_unitcands(F) = {
  my(bnf = F[1], W = F[5], tu = Mod(lift(bnf.tu[2]), W), w = bnf.tu[1], fu = bnf.fu, out = List());
  for(a = 0, w - 1, forvec(e = vector(#fu, k, [0, 1]),
    listput(out, tu^a * prod(k = 1, #fu, Mod(lift(fu[k]), W)^e[k]))));
  Vec(out);
}

lat_coprime(nf, Id) = {
  idealmul(nf, idealcoprime(nf, Id, 2), Id);
}

lat_sqrtnorm(F, r, L1, L2) = {
  my(bnf = F[1], nf = bnf.nf, W = F[5], tu = Mod(lift(bnf.tu[2]), W), w = bnf.tu[1], y, v);
  for(a = 0, w - 1,
    if(nfeltissquare(nf, r * tu^a, &y),
      v = Mod(lift(nfbasistoalg(nf, y)), W);
      for(b = 0, w - 1,
        if(v * tu^b * cm_conj(F, v * tu^b) == r && lat_hnf(lat_mulmat(nf, nfalgtobasis(nf, v * tu^b)) * L1) == L2, return(1)))));
  0;
}
"""

load_gp(GPLAT)

GPE = r"""
e_up(bnf, bnf0, sE, Id) = {
  my(two = idealtwoelt(bnf0, Id));
  idealhnf(bnf, two[1], subst(lift(nfbasistoalg(bnf0, two[2])), 'x, sE));
}

e_setup(q, e1, mu, nu) = {
  my(e2 = q * (16 * mu^2 + 16 * nu^2 - q), e3 = -2 * q^2 * (16 * mu^2 - 16 * nu^2 + e1), P, bnf, nf, sE, h, bnf0,
     cyc, k, M, rm, S);
  P = 'w^6 - e1 * 'w^5 + e2 * 'w^4 - e3 * 'w^3 + q^2 * e2 * 'w^2 - q^4 * e1 * 'w + q^6;
  if(!polisirreducible(P), return(0));
  bnf = bnfinit(P, 1); nf = bnf.nf;
  sE = Mod('w, P) + q^2 / Mod('w, P);
  h = minpoly(sE, 'x); bnf0 = bnfinit(h, 1);
  cyc = bnf.cyc; k = #cyc;
  if(k == 0, return([bnf, matrix(0, 0), [matrix(0, 0), matrix(0, 0), matrix(0, 0)], P]));
  M = matdiagonal(cyc);
  for(j = 1, #bnf0.gen, M = concat(M, Mat(bnfisprincipal(bnf, e_up(bnf, bnf0, sE, bnf0.gen[j]), 0))));
  rm = abs(nf.disc) / bnf0.disc^2;
  foreach(factor(rm)[, 1]~, l,
    foreach(idealprimedec(bnf0, l), pr,
      my(fa = idealfactor(nf, e_up(bnf, bnf0, sE, pr)));
      if(matsize(fa)[1] == 1 && fa[1, 2] == 2, M = concat(M, Mat(bnfisprincipal(bnf, fa[1, 1], 0))))));
  M = mathnf(M);
  S = matsnf(M, 1);
  [bnf, M, S, P];
}

e_class(ES, b, l) = {
  my(bnf = ES[1], P = ES[4], S = ES[3], Id, v, g, D);
  if(#bnf.cyc == 0, return([]));
  Id = idealadd(bnf, l, Mod('w - b, P));
  if(idealnorm(bnf, Id) != l, error("prime of E not of degree 1"));
  v = bnfisprincipal(bnf, Id, 0);
  g = S[1] * v; D = S[3];
  select(x -> x != -1, vector(#g, r, if(D[r, r] > 1, g[r] % D[r, r], -1)));
}

k_primes(F, q, mu, nu, e1, l, il) = {
  my(nf = F[1].nf, W = F[5], iK = Mod(lift(F[2]), W), chil, rts, out = List(), a, P1, P2);
  chil = Mod(1, l) * ('x^4 - 4 * (nu + mu * il) * 'x^3 + e1 * 'x^2 - 4 * q * (nu - mu * il) * 'x + q^2);
  rts = polrootsmod(chil, l);
  if(#rts != 4, return(0));
  for(j = 1, 4,
    a = lift(rts[j]);
    P1 = idealadd(nf, idealadd(nf, l, iK - il), Mod('z - a, W));
    P2 = idealadd(nf, idealadd(nf, l, iK + il), Mod('z - lift(Mod(q, l) / rts[j]), W));
    if(idealnorm(nf, P1) != l || idealnorm(nf, P2) != l, error("prime of K not of degree 1"));
    listput(out, [a, P1, P2]));
  Vec(out);
}

k_ideal(nf, prs, sub) = {
  my(B = idealhnf(nf, 1));
  for(j = 1, 4, B = idealmul(nf, B, if(bittest(sub, j - 1), prs[j][2], prs[j][3])));
  B;
}

lat_act(nf, L, B) = {
  lat_hnf(concat(L, 8 * idealmul(nf, idealinv(nf, B), lat_ideal(nf, L))));
}

e_embmat(ES, psis) = {
  my(zk = ES[1].nf.zk, M = matrix(6, 6), z);
  for(t = 1, 3, for(n = 1, 6,
    z = subst(lift(zk[n]), 'w, psis[t]);
    M[2 * t - 1, n] = real(z); M[2 * t, n] = imag(z)));
  M^(-1);
}

e_solve(Minv, vals) = {
  my(v = vectorv(6));
  for(t = 1, 3, v[2 * t - 1] = real(vals[t]); v[2 * t] = imag(vals[t]));
  Minv * v;
}

e_prime(ES, l, b) = {
  my(nf = ES[1].nf, P = ES[4], prs = idealprimedec(nf, l), pr = 0);
  for(k = 1, #prs, if(prs[k].f == 1 && nfeltval(nf, Mod('w - b, P), prs[k]) > 0, pr = prs[k]));
  if(pr === 0 || pr.e != 1, return(0));
  [nf, nfmodprinit(nf, pr), l, pr];
}

e_red(EP, c) = {
  my(y);
  if(nfeltval(EP[1], c, EP[4]) < 0, return(-1));
  y = nfmodpr(EP[1], c, EP[2]).pol;
  if(type(y) == "t_POL", y = polcoef(y, 0));
  lift(y) % EP[3];
}
"""

load_gp(GPE)


def order_submodules(LS):
    """The lattices M' with Z[i, pi, q/pi] = R in M' in O_K, as integer HNF matrices."""
    import numpy as np
    HR, S, Mi, Mp = LS[0], LS[1], LS[2], LS[3]
    U, D = S[0], S[2]
    d = [int(D[j, j]) for j in range(8)]
    idx = [j for j in range(8) if d[j] > 1]
    dd = [d[j] for j in idx]
    Ui = U ** (-1)

    def tomat(A):
        return np.array([[int(A[r, c]) for c in range(8)] for r in range(8)], dtype=object)
    Un, Uin = tomat(U), tomat(Ui)
    Ai, Ap = Un.dot(tomat(Mi)).dot(Uin), Un.dot(tomat(Mp)).dot(Uin)
    As = Un.dot(tomat(LS[6])).dot(Uin) if LS[6] != 0 else None
    n = len(idx)
    size = 1
    for x in dd:
        size *= x
    elts = list(itertools.product(*[range(x) for x in dd]))
    eidx = {e: k for k, e in enumerate(elts)}

    def app(A, a):
        return tuple(sum(int(A[idx[r]][idx[c]]) * a[c] for c in range(n)) % dd[r] for r in range(n))
    addt = np.array([[eidx[tuple((x + y) % m_ for x, y, m_ in zip(a, b, dd))] for b in elts] for a in elts],
                    dtype=np.int32)
    it = np.array([eidx[app(Ai, a)] for a in elts])
    pt = np.array([eidx[app(Ap, a)] for a in elts])
    st = np.array([eidx[app(As, a)] for a in elts]) if As is not None else pt

    def close(mask):
        while True:
            ids = np.nonzero(mask)[0]
            new = mask.copy()
            new[it[ids]] = True
            new[pt[ids]] = True
            new[st[ids]] = True
            ids2 = np.nonzero(new)[0]
            new[addt[np.ix_(ids2, ids2)].ravel()] = True
            if new.sum() == mask.sum():
                return new
            mask = new
    zero = np.zeros(size, dtype=bool)
    zero[0] = True
    cyc = []
    for k in range(size):
        mk = zero.copy()
        mk[k] = True
        cyc.append(close(mk))
    subs = {zero.tobytes(): zero}
    front = [zero]
    while front:
        nxt = []
        for U_ in front:
            ids1 = np.nonzero(U_)[0]
            for k in range(size):
                if U_[k]:
                    continue
                new = np.zeros(size, dtype=bool)
                new[addt[np.ix_(ids1, np.nonzero(cyc[k])[0])].ravel()] = True
                key = new.tobytes()
                if key not in subs:
                    subs[key] = new
                    nxt.append(new)
        front = nxt
    out = []
    for mask in subs.values():
        cols = [HR]
        lifts = []
        for k in np.nonzero(mask)[0]:
            a = elts[k]
            lifts.append([sum(int(Uin[r][idx[c]]) * a[c] for c in range(n)) for r in range(8)])
        Lm = pari.matconcat([HR, pari.matrix(8, len(lifts), [lifts[c][r] for r in range(8) for c in range(len(lifts))])])
        out.append((int(mask.sum()), pari.mathnf(Lm)))
    out.sort(key=lambda t: (t[0], str(t[1])))
    return out


TYPES = [(1, 2), (1, 3), (1, 4)]


class CMField:
    def __init__(self, q, e1, mu, nu):
        self.q, self.e1, self.mu, self.nu = q, e1, mu, nu
        F = G("cm_field")(q, e1, mu, nu)
        if F == 0:
            raise ValueError("chi conj(chi) is reducible")
        self.F = F
        self.bnf = F[0]
        self.W = F[4]
        self.cyc = [int(c) for c in pari("b -> b.cyc")(self.bnf)]
        self.aut = int(pari("b -> #nfgaloisconj(b.nf)")(self.bnf))
        self.reps = G("cm_classreps")(F)

    def embeddings(self):
        return G("cm_embeddings")(self.F)

    def lattices(self):
        """All (class exponents, ideal, xi, type, parity, T) with xi totally imaginary of the sign pattern of the
        type, up to xi -> xi v conj(v).  T is the standard basis (rows of K-elements) when the lattice is even."""
        F = self.F
        pl = self.embeddings()
        out = []
        for e, I in self.reps:
            xis = G("cm_polarisations")(F, I)
            if len(xis) == 0:
                continue
            rows = G("cm_qicoords")(F, I)
            den = pari.denominator(pari.Vec(pari.concat(rows)))
            rows = gauss_rowreduce([[x * den for x in r] for r in rows])
            assert len(rows) == 4
            basis = [G("cm_elt")(F, [x / den for x in r]) for r in rows]
            for S in TYPES:
                kept = []
                for xi in xis:
                    im = [float(pari.imag(pari.subst(pari.lift(xi), "z", pl[j]))) for j in range(4)]
                    if all((im[j] < 0) == ((j + 1) in S) for j in range(4)):
                        if not any(G("cm_isnorm")(F, xi / x2) for x2 in kept):
                            kept.append(xi)
                for xi in kept:
                    Gm = [[G("cm_form")(F, xi, basis[a], basis[b]) for b in range(4)] for a in range(4)]
                    det = pari.matdet(pari.matrix(4, 4, [x for r in Gm for x in r]))
                    assert det == 1, det
                    even = all(pari.imag(Gm[a][a]) % 2 == 0 for a in range(4))
                    T = None
                    if even:
                        Tc = standard_basis(Gm)
                        if Tc is not None:
                            T = [sum((Tc[a][k] * basis[k] for k in range(4)), pari(0)) for a in range(4)]
                    out.append({"e": [int(x) for x in e], "ideal": I, "xi": xi, "type": S, "even": bool(even), "T": T})
        return out

    def order_lattices(self, verbose=True, level=True):
        """All (L, xi) with L a lattice stable under R, xi L = conj(trace dual of L), xi totally imaginary with the
        sign pattern of a type S, up to (L, xi) -> (v L, xi / (v conj v)).

        R = Z[i, pi, q/pi, sigma, conj(sigma)] with sigma = (pi - 1)/(1 + i) when level is set (Frobenius is trivial
        on L/(1 + i)L: the six lines are rational), and R = Z[i, pi, q/pi] otherwise.

        Every R-lattice is, up to scaling, aa . (u M' + 8 O_K) with aa a class representative prime to 2, M' one of
        the lattices between R and O_K, and u a unit at 2 modulo the stabiliser of M'.
        """
        if getattr(self, "_lat", None) is not None and self._lat[0] == level:
            return [dict(r) for r in self._lat[1]]
        F, bnf = self.F, self.bnf
        nf = pari("b -> b.nf")(bnf)
        pl = self.embeddings()
        LS = G("lat_setup")(F, self.q, 1 if level else 0)
        HR, Cm = LS[0], LS[4]
        id8 = pari("8 * matid(8)")
        dual = pari("(LS, L) -> lat_hnf((L~ * LS[6])^(-1))")
        mulmat = G("lat_mulmat")
        mod8 = G("lat_mod8")
        subs = order_submodules(LS)
        gens = G("lat_units8")(nf)
        one = pari("matid(8)[, 1]")
        seen = {str(mod8(nf, one, HR)): one}
        front = [one]
        while front:
            nxt = []
            for u in front:
                for g in gens:
                    v = pari("(nf, a, b) -> nfeltmul(nf, a, b) % 8")(nf, u, g)
                    key = str(mod8(nf, v, HR))
                    if key not in seen:
                        seen[key] = v
                        nxt.append(v)
            front = nxt
        units = list(seen.values())
        assert pari("H -> denominator(8 * H^(-1)) == 1")(HR), "8 O_K is not contained in R"
        self.info = {"index_R": int(pari.matdet(HR)), "between": len(subs), "principal_mod_R": len(units)}
        if verbose:
            print(f"  [O_K : R] = {self.info['index_R']}, lattices between R and O_K: {len(subs)}, "
                  f"units at 2 modulo units of R: {len(units)}", flush=True)
        cands = G("lat_unitcands")(F)
        reps = [(e, G("lat_coprime")(nf, I)) for e, I in self.reps]
        out = []
        nself = 0
        for mi, (size, Mh) in enumerate(subs):
            N = G("lat_dualconj")(LS, Mh)
            C = dual(LS, G("lat_prod")(nf, Mh, dual(LS, N)))
            if G("lat_prod")(nf, C, Mh) != N:
                continue
            cgen = None
            # a generator of N over M' at 2, if there is one, is among the columns of a basis of the colon
            for b in [C[k] for k in range(8)]:
                ratio = abs(pari.matdet(mulmat(nf, b))) * pari.matdet(Mh) / pari.matdet(N)
                if ratio != 0 and pari.valuation(ratio, 2) == 0:
                    cgen = b
                    break
            if cgen is None:
                continue
            nself += 1
            M8 = pari.mathnf(pari.concat(Mh, id8))
            Q = {}
            for u in units:
                key = str(mod8(nf, u, Mh))
                if key not in Q:
                    ub = pari("(C, u) -> (C * u) % 8")(Cm, u)
                    uu = pari("(nf, a, b) -> nfeltmul(nf, a, b) % 8")(nf, u, ub)
                    Q[key] = (u, mod8(nf, uu, Mh))
            NI = G("lat_ideal")(nf, N)
            cel = pari("(nf, b) -> nfbasistoalg(nf, b)")(nf, cgen)
            for e, I in reps:
                Ib = pari("(nf, c, I) -> nfgaloisapply(nf, c, I)")(nf, pari.lift(F[2]), I)
                Cid = pari("(nf, NI, I, Ib) -> idealmul(nf, NI, idealinv(nf, idealmul(nf, I, Ib)))")(nf, NI, I, Ib)
                v = pari.bnfisprincipal(bnf, Cid)
                if v[0] != 0:
                    continue
                xi0 = pari("(nf, c) -> nfbasistoalg(nf, c)")(nf, v[1])
                for eps in cands:
                    xi = pari("(a, b, W) -> Mod(lift(a) * lift(b), W)")(xi0, eps, self.W)
                    if G("cm_conj")(F, xi) != -pari.lift(xi) and pari.Mod(G("cm_conj")(F, xi), self.W) != -xi:
                        continue
                    im = [float(pari.imag(pari.subst(pari.lift(xi), "z", pl[j]))) for j in range(4)]
                    S = tuple(j + 1 for j in range(4) if im[j] < 0)
                    if S not in TYPES:
                        continue
                    x = pari("(nf, a, b) -> nfalgtobasis(nf, a / b)")(nf, xi, cel)
                    assert int(pari.denominator(x)) % 2 == 1
                    x8 = pari("x -> lift(Mod(x, 8))")(x)
                    for key, (u, nu) in Q.items():
                        if mod8(nf, x8, nu) != M8:
                            continue
                        Lu = mod8(nf, u, Mh)
                        L = dual(LS, G("lat_hnf")(pari.concat(dual(LS, pari.idealhnf(nf, I)), dual(LS, Lu))))   # aa cap Lu
                        xic = pari("(nf, a) -> nfalgtobasis(nf, a)")(nf, xi)
                        assert G("lat_hnf")(mulmat(nf, xic) * L) == G("lat_dualconj")(LS, L), "not self-dual"
                        out.append({"e": [int(t) for t in e], "L": L, "xi": xi, "type": S, "M": mi,
                                    "index": int(pari.matdet(Mh))})
        self.info["self_dual_types"] = nself
        # equivalence under units
        keep = []
        for r in out:
            dup = False
            for r2 in keep:
                if r2["type"] != r["type"] or r2["e"] != r["e"] or r2["index"] != r["index"]:
                    continue
                if G("lat_sqrtnorm")(F, r2["xi"] / r["xi"], r2["L"], r["L"]):
                    dup = True
                    break
            if not dup:
                keep.append(r)
        if verbose:
            print(f"  self-dual local types: {nself}; unimodular lattices found: {len(out)}, up to units: {len(keep)}",
                  flush=True)
        for r in keep:
            r["even"], r["T"] = self.standardise(r["L"], r["xi"])
        self._lat = (level, [dict(r) for r in keep])
        return keep

    def standardise(self, L, xi):
        """(even, T): parity of the Z[i]-lattice L with the form 2 Tr_{K/Q(i)}(xi x conj y), and a standard basis."""
        F = self.F
        rows = G("lat_qi")(F, L)
        den = pari.denominator(pari.Vec(pari.concat(rows)))
        rows = gauss_rowreduce([[x * den for x in rw] for rw in rows])
        assert len(rows) == 4
        basis = [G("cm_elt")(F, [x / den for x in rw]) for rw in rows]
        Gm = [[G("cm_form")(F, xi, basis[a], basis[b]) for b in range(4)] for a in range(4)]
        det = pari.matdet(pari.matrix(4, 4, [x for rw in Gm for x in rw]))
        assert det == 1, det
        even = all(pari.imag(Gm[a][a]) % 2 == 0 for a in range(4))
        T = None
        if even:
            Tc = standard_basis(Gm)
            if Tc is not None:
                T = [sum((Tc[a][k] * basis[k] for k in range(4)), pari(0)) for a in range(4)]
        return bool(even), T

    def tau(self, lat, pl=None):
        pl = pl or self.embeddings()
        t = G("cm_tau")(self.F, pl, list(lat["type"]), lat["T"])
        return t


def invariants(J, chart="J2"):
    """Degree-0 functions on P(2,3,4,5,6).  Chart J2: J4/J2^2, J6/J2^3, J3^2/J2^3, J3 J5/J2^4, J5^2/J2^5.
    Chart J3 (for surfaces with J2 = 0 modulo q): J2^3/J3^2, J2 J4/J3^2, J6/J3^2, J4 J5/J3^3, J4^3/J3^4, J5^3/J3^5,
    J2^2 J5/J3^3."""
    J2, J3, J4, J5, J6 = J[0], J[1], J[2], J[3], J[4]
    if chart == "J2":
        return [J4 / J2 ** 2, J6 / J2 ** 3, J3 ** 2 / J2 ** 3, J3 * J5 / J2 ** 4, J5 ** 2 / J2 ** 5]
    return [J2 ** 3 / J3 ** 2, J2 * J4 / J3 ** 2, J6 / J3 ** 2, J4 * J5 / J3 ** 3, J4 ** 3 / J3 ** 4, J5 ** 3 / J3 ** 5,
            J2 ** 2 * J5 / J3 ** 3]


def cm_points(fld, digits, verbose=True, order=True, level=True):
    setprec(digits + 20)
    R2 = R2_for(digits)
    pl = fld.embeddings()
    lats = fld.order_lattices(verbose, level) if order else fld.lattices()
    pts = []
    for lat in lats:
        if lat["T"] is None:
            lat["J"] = None
            continue
        t0 = time.time()
        tau = fld.tau(lat, pl)
        Y = pari("t -> (t - conj(t~)) / (2 * I)")(tau)
        assert float(pari.real(Y[0, 0])) > 0 and float(pari.real(pari.matdet(Y))) > 0, "tau not in D"
        taur = G("cm_reduce")(tau)
        Yr = pari("t -> (t - conj(t~)) / (2 * I)")(taur)
        lat["tau"] = taur
        lat["detY"] = float(pari.real(pari.matdet(Yr)))
        lat["J"] = G("cm_J")(G("cm_tsq")(taur, R2))
        if verbose:
            print(f"  class {lat['e']} index {lat.get('index', 1)} type {lat['type']}: det Y = {lat['detY']:.3f}, {time.time() - t0:.1f} s",
                  flush=True)
    return lats


def cmd_points(a):
    fld = CMField(a.q, a.e1, a.mu, a.nu)
    print(f"K = Q[x]/({fld.W}), class group {fld.cyc}, automorphisms {fld.aut}, class representatives {len(fld.reps)}")
    lats = cm_points(fld, a.digits, order=not a.maximal, level=not a.nolevel)
    for S in TYPES:
        ev = [l for l in lats if l["type"] == S and l["even"]]
        od = [l for l in lats if l["type"] == S and not l["even"]]
        print(f"type {S}: {len(ev)} even lattices, {len(od)} odd; with a standard basis: {sum(l['T'] is not None for l in ev)}")
        for l in ev:
            if l["J"] is None:
                continue
            inv = invariants(l["J"])
            print("   class", l["e"], "index", l.get("index", 1), " ".join(f"{complex(x):.12g}" for x in inv))


# ------------------------------------------------------------------------------------------------ comparison

INV_NAMES = ["J4/J2^2", "J6/J2^3", "J3^2/J2^3", "J3 J5/J2^4", "J5^2/J2^5"]


def rational_poly(vals, digits):
    """prod (X - v) with rational coefficients, or None.  The test is that each coefficient is within
    10^-(digits - 15) relative of a rational whose numerator and denominator together have fewer than 0.6 * digits
    digits."""
    pol = pari("v -> my(P = 1); for(k = 1, #v, P *= ('X - v[k])); P")(vals)
    return rationalise(pol, digits)


def rationalise(pol, digits):
    n = int(pari.poldegree(pol))
    coeffs = []
    for k in range(n + 1):
        c = pari.polcoef(pol, k)
        if float(abs(pari.imag(c))) > 10.0 ** (-(digits - 15)) * max(1.0, float(abs(c))):
            return None
        c = pari.real(c)
        r = pari.bestappr(c, pari(10) ** int(0.3 * digits))
        err = abs(c - r)
        size = len(str(pari.numerator(r))) + len(str(pari.denominator(r)))
        if float(err) > 10.0 ** (-(digits - 15)) * max(1.0, float(abs(c))) or size > 0.6 * digits:
            return None
        coeffs.append(r)
    return pari.Pol(list(reversed(coeffs)), "X")


def interpolation_poly(xs, ys):
    """sum_k y_k prod_{j != k} (X - x_j)."""
    return pari("(x, y) -> sum(k = 1, #x, y[k] * prod(j = 1, #x, if(j == k, 1, 'X - x[j])))")(xs, ys)


def config_invariants(cfg, p, chart="J2"):
    b, c, e, f = cfg
    A = pari.matrix(3, 6, [1, 0, 0, 1, 1, 1, 0, 1, 0, 1, b, e, 0, 0, 1, 1, c, f])
    J = G("cm_J")(G("cm_tmat")(A * pari.Mod(1, p)))
    if J[0 if chart == "J2" else 1] == 0:
        return None
    return tuple(int(pari.lift(x)) for x in invariants(J, chart))


def load_surfaces(q, e1, e2):
    """Configurations of sixlines_cm.jsonl with p = q and coefficients e1, e2 (e1 up to the sign of the twist)."""
    rows = []
    with open(os.path.join(HERE, "sixlines_cm.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            if r["p"] == q and abs(r["e"][0]) == abs(e1) and r["e"][1] == e2:
                rows.append(r)
    return rows


def cm_data(fld, digits, verbose=True):
    """The even lattices with their invariants, grouped by type."""
    lats = cm_points(fld, digits, verbose=verbose)
    pts = [l for l in lats if l.get("J") is not None]
    return lats, pts


def rat(c, digits):
    """A rational within 10^-(digits - 15) (relative) of the real number c, of total size below 0.6 * digits."""
    r = pari.bestappr(c, pari(10) ** int(0.3 * digits))
    err = abs(c - r)
    size = len(str(pari.numerator(r))) + len(str(pari.denominator(r)))
    if float(err) > 10.0 ** (-(digits - 15)) * max(1.0, float(abs(c))) or size > 0.6 * digits:
        return None
    return r


class Recogniser:
    """Elements of E = Q[w]/P, as coordinates on the integral basis, from their values at the three complex
    embeddings w -> pi_j pi_k, {j, k} a type."""

    def __init__(self, fld, ES, digits, conj=False):
        pl = fld.embeddings()
        psis = [pl[T[0] - 1] * pl[T[1] - 1] for T in TYPES]
        if conj:
            psis = [pari.conj(x) for x in psis]
        self.Minv = G("e_embmat")(ES, psis)
        self.digits = digits

    def element(self, vals):
        r = G("e_solve")(self.Minv, vals)
        out = []
        for k in range(6):
            x = rat(r[k], self.digits)
            if x is None:
                return None
            out.append(x)
        return pari.Col(out)

    def poly(self, polys):
        n = int(pari.poldegree(polys[0]))
        coeffs = []
        for k in range(n + 1):
            c = self.element([pari.polcoef(p, k) for p in polys])
            if c is None:
                return None
            coeffs.append(c)
        return coeffs


def reduce_poly(coeffs, EP):
    """The polynomial with coefficients in E reduced at a prime of degree 1; None if a coefficient is not integral
    there."""
    out = []
    for c in coeffs:
        y = int(G("e_red")(EP, c))
        if y < 0:
            return None
        out.append(pari.Mod(y, EP[2]))
    return pari.Pol(list(reversed(out)), "X")


def denominators(coeffs):
    d = pari(1)
    for c in coeffs:
        d = pari.lcm(d, pari.denominator(c))
    return d


class Torsor:
    """The even lattices of a field, by type, with the class polynomials over E."""

    def __init__(self, fld, digits, verbose=True):
        self.fld, self.digits = fld, digits
        self.lats, pts = cm_data(fld, digits, verbose)
        self.pts = pts
        self.set_chart("J2")
        self.by = {T: [p for p in pts if p["type"] == T] for T in TYPES}
        self.e2 = fld.q * (16 * fld.mu ** 2 + 16 * fld.nu ** 2 - fld.q)
        self.N = len(self.by[TYPES[0]])
        assert all(len(self.by[T]) == self.N for T in TYPES)
        self.ES = G("e_setup")(fld.q, fld.e1, fld.mu, fld.nu)
        assert self.ES != 0, "P is reducible"
        self.rec = Recogniser(fld, self.ES, digits)
        self.rec_conj = Recogniser(fld, self.ES, digits, conj=True)
        self.mix = None

    def set_chart(self, chart):
        self.chart = chart
        for p in self.pts:
            p["inv"] = invariants(p["J"], chart)

    def prime(self, l, b):
        """The prime (l, w - b) of E, of degree 1 and unramified; 0 otherwise."""
        return G("e_prime")(self.ES, l, b)

    def coordinate(self, mix):
        return {T: [sum((m_ * v for m_, v in zip(mix, p["inv"])), pari(0)) for p in self.by[T]] for T in TYPES}

    def setup(self, q, bq, verbose=True):
        """Choose a coordinate whose class polynomial is recognised over E, integral at the prime (q, w - bq) and has
        N simple roots modulo it."""
        EP = self.prime(q, bq)
        assert EP != 0
        for chart in ("J2", "J3"):
            self.set_chart(chart)
            if self._setup_chart(q, EP, verbose):
                return True
        return False

    def _setup_chart(self, q, EP, verbose):
        for mix in ([1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 1, 0, 0, 0], [1, 2, 3, 0, 0], [3, 1, 2, 5, 7],
                    [2, -1, 3, 1, -2], [5, 3, -1, 2, 1]):
            xs = self.coordinate(mix)
            H = self.rec.poly([pari("v -> my(P = 1); for(k = 1, #v, P *= ('X - v[k])); P")(xs[T]) for T in TYPES])
            if H is None:
                if verbose:
                    print(f"  chart {self.chart}, coordinate {mix}: not recognised at {self.digits} digits")
                continue
            Hq = reduce_poly(H, EP)
            nroots = -1 if Hq is None else len(pari.polrootsmod(Hq, q))
            sqf = Hq is not None and pari.poldegree(pari.gcd(Hq, pari.deriv(Hq))) == 0
            if verbose:
                print(f"  chart {self.chart}, coordinate {mix}: denominators {pari.factor(denominators(H))[0]}, "
                      f"roots mod {q}: {nroots} (-1: not integral at the prime), squarefree: {bool(sqf)}")
            if sqf and nroots == self.N:
                self.mix, self.xs, self.H = mix, xs, H
                self.Gs = self.invariant_polys()
                if any(g is None for g in self.Gs) or any(reduce_poly(g, EP) is None for g in self.Gs):
                    if verbose:
                        print(f"  the polynomials of the invariants are not recognised at {self.digits} digits")
                    self.mix = None
                    return False
                return True
        return False

    def interpolate(self, ys, conj=False):
        """The polynomial of E[X] whose value at the type T is sum_j ys[T][j] prod_{l != j} (X - xs[T][l])."""
        rec = self.rec_conj if conj else self.rec
        return rec.poly([interpolation_poly(self.xs[T], ys[T]) for T in TYPES])

    def invariant_polys(self):
        n = len(self.pts[0]["inv"])
        return [self.interpolate({T: [p["inv"][k] for p in self.by[T]] for T in TYPES}) for k in range(n)]

    def perm_poly(self, perms, conj=False):
        return self.interpolate({T: [self.xs[T][perms[T][j]] for j in range(self.N)] for T in TYPES}, conj)

    def locate(self, L, xi, T, digits=50):
        """Index among the lattices of type T of the point of (L, xi)."""
        old = int(pari.default("realprecision"))
        try:
            setprec(digits + 10)
            even, Tb = self.fld.standardise(L, xi)
            assert even and Tb is not None, "the image is not an even lattice"
            tau = G("cm_tau")(self.fld.F, self.fld.embeddings(), list(T), Tb)
            inv = invariants(G("cm_J")(G("cm_tsq")(G("cm_reduce")(tau), R2_for(digits))), self.chart)
            best = []
            for j, p in enumerate(self.by[T]):
                d = max(float(abs(inv[k] - p["inv"][k]) / (abs(inv[k]) + abs(p["inv"][k]))) for k in range(len(inv)))
                best.append((d, j))
            best.sort()
            assert best[0][0] < 1e-12 and (len(best) == 1 or best[1][0] > 1e-6), best[:2]
            return best[0][1]
        finally:
            setprec(old)

    def act(self, B, beta):
        """The permutation, for each type, induced by (L, xi) -> (L + 8 B^-1 (L O_K), beta xi)."""
        nf = pari("b -> b.nf")(self.fld.bnf)
        perms = {}
        for T in TYPES:
            perm = []
            for p in self.by[T]:
                L2 = G("lat_act")(nf, p["L"], B)
                perm.append(self.locate(L2, p["xi"] * beta, T))
            assert sorted(perm) == list(range(self.N))
            perms[T] = tuple(perm)
        return perms


def enumerated_points(q, e1, e2, chart="J2"):
    enum_set = {}
    for r in load_surfaces(q, e1, e2):
        for cfg in r["configs"]:
            enum_set.setdefault(config_invariants(cfg, q, chart), []).append(tuple(cfg))
    return enum_set


def reduced_points(tor, q, bq):
    """{root of H mod q: tuple of the five invariants}."""
    EP = tor.prime(q, bq)
    Hq = reduce_poly(tor.H, EP)
    dH = pari.deriv(Hq)
    Gs = [reduce_poly(g, EP) for g in tor.Gs]
    out = {}
    for r in pari.polrootsmod(Hq, q):
        d = pari.subst(dH, "X", r)
        out[int(pari.lift(r))] = tuple(int(pari.lift(pari.subst(g, "X", r) / d)) for g in Gs)
    return out


def cmd_compare(a, tor=None):
    q = a.q
    if tor is None:
        fld = CMField(a.q, a.e1, a.mu, a.nu)
        print(f"K = Q[z]/({fld.W}), class group {fld.cyc}, automorphisms {fld.aut}")
        tor = Torsor(fld, a.digits)
    fld = tor.fld
    for S in TYPES:
        ev = [l for l in tor.lats if l["type"] == S and l["even"]]
        od = [l for l in tor.lats if l["type"] == S and not l["even"]]
        print(f"type {S}: {len(ev)} even lattices (indices in O_K {sorted(set(l['index'] for l in ev))}), {len(od)} odd")
    rec = {"cmd": "compare", "q": q, "e1": a.e1, "mu": a.mu, "nu": a.nu, "digits": a.digits, "class_group_K": fld.cyc,
           "index_R": fld.info["index_R"], "even_per_type": tor.N,
           "odd_per_type": sum(1 for l in tor.lats if l["type"] == TYPES[0] and not l["even"])}
    if tor.mix is None and not tor.setup(q, a.e1 % q):
        print("no usable coordinate")
        rec["status"] = "class polynomial not obtained"
        return rec
    enum_set = enumerated_points(q, a.e1, tor.e2, tor.chart)
    rec["enumerated_configurations"] = sum(len(v) for v in enum_set.values())
    rec["enumerated_points"] = len(enum_set)
    rec["chart"] = tor.chart
    red = reduced_points(tor, q, a.e1 % q)
    cm_set = set(red.values())
    common = cm_set & set(enum_set)
    rec.update({"status": "ok", "coordinate": tor.mix, "denominator_primes": [int(x) for x in pari.factor(denominators(tor.H))[0]],
                "cm_points_mod_q": len(cm_set), "in_both": len(common),
                "enumerated_only": len(set(enum_set) - cm_set), "cm_only": len(cm_set - set(enum_set))})
    print(f"enumerated: {rec['enumerated_configurations']} configurations, {len(enum_set)} moduli points")
    print(f"CM method: {len(cm_set)} points over F_{q} at the prime (q, F - e1) of E")
    print(f"in both: {len(common)}; enumerated only: {rec['enumerated_only']}; CM only: {rec['cm_only']}")
    for t in sorted(cm_set):
        print("   ", t, "enumerated as", enum_set.get(t, "-"))
    return rec


# ------------------------------------------------------------------------------------------------ the action

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def compose(p1, p2):
    """First p2, then p1."""
    return tuple(p1[x] for x in p2)


def inverse(p):
    out = [0] * len(p)
    for a_, b_ in enumerate(p):
        out[b_] = a_
    return tuple(out)


def split_primes(fld, tor, count, bound=20000):
    """Primes l = 1 mod 4, l != q, at which chi has four roots a_1..a_4 modulo (l, i - i_l) with six distinct
    products, the primes of E at these products being of degree 1."""
    out = []
    l = 4
    while len(out) < count and l < bound:
        l = int(pari.nextprime(l + 1))
        if l % 4 != 1 or l == fld.q or fld.q % l == 0:
            continue
        il = int(pari("l -> lift(sqrt(Mod(-1, l)))")(l))
        prs = G("k_primes")(fld.F, fld.q, fld.mu, fld.nu, fld.e1, l, il)
        if prs == 0:
            continue
        a = [int(prs[j][0]) for j in range(4)]
        if len(set(a)) < 4:
            continue
        b = {pr: a[pr[0]] * a[pr[1]] % l for pr in PAIRS}
        if len(set(b.values())) < 6:
            continue
        EPs = {pr: tor.prime(l, b[pr]) for pr in PAIRS}
        if any(EPs[pr] == 0 for pr in PAIRS):
            continue
        out.append({"l": l, "prs": prs, "b": b, "EP": EPs})
    return out


def cmd_action(a, tor=None):
    q = a.q
    if tor is None:
        fld = CMField(a.q, a.e1, a.mu, a.nu)
        print(f"K = Q[z]/({fld.W}), class group {fld.cyc}")
        tor = Torsor(fld, a.digits, verbose=False)
    fld = tor.fld
    nf = pari("b -> b.nf")(fld.bnf)
    N = tor.N
    cycG = [int(x) for x in pari("E -> my(D = E[3][3]); select(x -> x > 1, vector(#D, r, D[r, r]))")(tor.ES)] \
        if len(pari("E -> E[1].cyc")(tor.ES)) else []
    orderG = 1
    for x in cycG:
        orderG *= x
    print(f"even lattices per type: {N};  Cl(E) = {pari('E -> E[1].cyc')(tor.ES)},  G = Cl(E)/H = {cycG} of order {orderG}")
    rec = {"cmd": "action", "q": q, "e1": a.e1, "digits": a.digits, "N": N, "G": cycG, "class_group_K": fld.cyc}
    if tor.mix is None and not tor.setup(q, a.e1 % q):
        print("no usable coordinate")
        rec["status"] = "class polynomial not obtained"
        return rec
    primes = split_primes(fld, tor, a.primes)
    print(f"primes l used: {[P['l'] for P in primes]}")
    rec["primes"] = [P["l"] for P in primes]

    def add(c1, c2):
        return tuple((x + y) % m_ for x, y, m_ in zip(c1, c2, cycG))

    def mul(c, k):
        return tuple((x * k) % m_ for x, m_ in zip(c, cycG))

    # ---- characteristic 0: the permutations
    t0 = time.time()
    table = {}          # class in G -> permutations by type
    conflicts = 0
    for P in primes:
        l = P["l"]
        P["cls"] = {pr: tuple(int(x) for x in G("e_class")(tor.ES, P["b"][pr], l)) for pr in PAIRS}
        P["perm"] = {}
        subs = range(16) if a.all_subsets else [sum(1 << j for j in pr) for pr in PAIRS]
        for sub in subs:
            B = G("k_ideal")(nf, P["prs"], sub)
            P["perm"][sub] = tor.act(B, l)
        for pr in PAIRS:
            sub = (1 << pr[0]) | (1 << pr[1])
            c = P["cls"][pr]
            if c in table and table[c] != P["perm"][sub]:
                conflicts += 1
            table.setdefault(c, P["perm"][sub])
        print(f"  l = {l}: classes of the six primes of E {[P['cls'][pr] for pr in PAIRS]}, {time.time() - t0:.0f} s", flush=True)
    npairs = 6 * len(primes)
    print(f"[4a] ideals P_j P_k conj(P_l P_m) against the class of the prime of E at a_j a_k: {npairs} ideals, "
          f"{len(table)} of {orderG} classes of G seen, conflicts (same class, different permutation): {conflicts}")
    zero = tuple(0 for _ in cycG)
    ident = tuple(range(N))
    hom_ok = hom_n = 0
    for c1 in table:
        for c2 in table:
            c3 = add(c1, c2)
            if c3 in table:
                hom_n += 1
                hom_ok += all(table[c3][T] == compose(table[c1][T], table[c2][T]) for T in TYPES)
    distinct = len({tuple(table[c][T] for T in TYPES) for c in table})
    fixed = sum(1 for c in table if c != zero for T in TYPES for j in range(N) if table[c][T][j] == j)
    triv = (zero not in table) or all(table[zero][T] == ident for T in TYPES)
    print(f"     products: {hom_ok} of {hom_n} pairs of classes compose correctly; distinct permutations: {distinct} "
          f"for {len(table)} classes; trivial class acts trivially: {triv}; fixed points of non-trivial classes: {fixed}")
    rec.update({"ideals": npairs, "classes_seen": len(table), "conflicts": conflicts, "products_ok": hom_ok,
                "products": hom_n, "distinct_permutations": distinct, "fixed_points": fixed, "trivial_ok": bool(triv)})
    if a.all_subsets:
        sq_ok = sq_n = trans = other_n = 0
        for P in primes:
            for sub in range(16):
                S = [j for j in range(4) if sub >> j & 1]
                if len(S) == 2:
                    continue
                other_n += 1
                perm = P["perm"][sub]
                if any(all(table[c][T] == perm[T] for T in TYPES) for c in table):
                    trans += 1
                c = zero
                for pr, co in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
                    d = sum(1 for j in pr if j in S) - sum(1 for j in co if j in S)
                    c = add(c, mul(P["cls"][pr], d % (max(cycG) if cycG else 1) if cycG else 0))
                if c in table:
                    sq_n += 1
                    sq_ok += all(compose(perm[T], perm[T]) == table[c][T] for T in TYPES)
        print(f"     the other ten ideals at each prime: {other_n} ideals; act as an element of G: {trans}; "
              f"square equals the predicted class: {sq_ok} of {sq_n}")
        rec.update({"other_ideals": other_n, "other_in_G": trans, "other_square_ok": sq_ok, "other_square_tested": sq_n})

    # ---- Frobenius in characteristic l
    cache = {}

    def poly_of(perms, inv=False):
        key = (tuple(perms[T] for T in TYPES), inv)
        if key not in cache:
            p2 = {T: inverse(perms[T]) for T in TYPES} if inv else perms
            cache[key] = tor.perm_poly(p2)
        return cache[key]

    def frob_test(A, EP, l):
        Hl = reduce_poly(tor.H, EP)
        if A is None or Hl is None:
            return None
        Al = reduce_poly(A, EP)
        if Al is None or pari.poldegree(pari.gcd(Hl, pari.deriv(Hl))) != 0:
            return None
        return bool(pari("(H, A, l) -> (lift(Mod('X, H)^l) * deriv(H) - A) % H == 0")(Hl, Al, l))

    ok = n = skipped = 0
    inv_fail = inv_n = mis_fail = mis_n = mis_same_ok = mis_same = 0
    for P in primes:
        l = P["l"]
        for pr in PAIRS:
            sub = (1 << pr[0]) | (1 << pr[1])
            perms = P["perm"][sub]
            r = frob_test(poly_of(perms), P["EP"][pr], l)
            if r is None:
                skipped += 1
                continue
            n += 1
            ok += r
            if any(inverse(perms[T]) != perms[T] for T in TYPES):
                r2 = frob_test(poly_of(perms, inv=True), P["EP"][pr], l)
                if r2 is not None:
                    inv_n += 1
                    inv_fail += (not r2)
            for pr2 in PAIRS:
                if pr2 == pr:
                    continue
                sub2 = (1 << pr2[0]) | (1 << pr2[1])
                r3 = frob_test(poly_of(P["perm"][sub2]), P["EP"][pr], l)
                if r3 is None:
                    continue
                if P["perm"][sub2] == perms:
                    mis_same += 1
                    mis_same_ok += r3
                else:
                    mis_n += 1
                    mis_fail += (not r3)
    print(f"[4b] Frobenius in characteristic l: x -> x^l on the roots of the class polynomial modulo the prime of E at "
          f"a_j a_k equals the action of P_j P_k conj(P_l P_m): {ok} of {n} (skipped {skipped})")
    print(f"     controls: the inverse permutation fails in {inv_fail} of {inv_n}; the ideal of another pair with a "
          f"different permutation fails in {mis_fail} of {mis_n} (same permutation: holds in {mis_same_ok} of {mis_same})")
    rec.update({"frobenius_ok": ok, "frobenius_tested": n, "frobenius_skipped": skipped, "inverse_fails": inv_fail,
                "inverse_tested": inv_n, "mismatch_fails": mis_fail, "mismatch_tested": mis_n})

    # ---- labels over F_q
    EPq = tor.prime(q, a.e1 % q)
    red = reduced_points(tor, q, a.e1 % q)
    Hq = reduce_poly(tor.H, EPq)
    dH = pari.deriv(Hq)
    roots = sorted(red)
    moves = {}
    for c in table:
        A = poly_of(table[c])
        Aq = None if A is None else reduce_poly(A, EPq)
        if Aq is None:
            continue
        moves[c] = {r: int(pari.lift(pari.subst(Aq, "X", pari.Mod(r, q)) / pari.subst(dH, "X", pari.Mod(r, q)))) for r in roots}
    label = {roots[0]: zero}
    todo = [roots[0]]
    while todo:
        r = todo.pop()
        for c, mv in moves.items():
            r2 = mv[r]
            if r2 in red and r2 not in label:
                label[r2] = add(label[r], c)
                todo.append(r2)
    cons = cons_n = 0
    for c, mv in moves.items():
        for r in roots:
            if r in label and mv[r] in label:
                cons_n += 1
                cons += label[mv[r]] == add(label[r], c)
    enum_set = enumerated_points(q, a.e1, tor.e2, tor.chart)
    signs = {}
    try:
        with open(os.path.join(HERE, f"sixlines_p{q}.jsonl")) as fh:
            for line in fh:
                r_ = json.loads(line)
                signs[tuple(r_["cfg"])] = r_["S"][0]
    except OSError:
        pass
    print(f"[4c] labels over F_{q}: {len(label)} of {len(roots)} surfaces labelled; every class acts on the labels as "
          f"its translation: {cons} of {cons_n}")
    rows = []
    for r in roots:
        cfgs = enum_set.get(red[r], [])
        rows.append({"root": r, "invariants": list(red[r]), "label": list(label.get(r, ())),
                     "configurations": [list(c) for c in cfgs], "S1": [signs.get(c) for c in cfgs]})
        print(f"     label {label.get(r)}  invariants {red[r]}  configurations {cfgs}  S_1 {[signs.get(c) for c in cfgs]}")
    rec.update({"status": "ok", "chart": tor.chart, "labelled": len(label), "label_consistency": [cons, cons_n],
                "surfaces": rows})
    return rec


# ------------------------------------------------------------------------------------------------ sweep

def mu_nu(q, e):
    import math
    e1, e2, e3 = e
    c2, c3 = e2 - 3 * q * q, e3 - 2 * q * q * e1
    A = 8 * q ** 3 - 4 * q * q * e1 + 2 * q * c2 - c3
    B = 8 * q ** 3 + 4 * q * q * e1 + 2 * q * c2 + c3
    return math.isqrt(A // (64 * q * q)), math.isqrt(B // (64 * q * q))


def sweep_one(job):
    q, e, digits_list, nprimes, all_subsets, maxn = job
    import io
    import contextlib
    mu, nu = mu_nu(q, e)
    base = {"q": q, "e": e, "mu_nu": [mu, nu]}
    buf = io.StringIO()
    out = []
    try:
        with contextlib.redirect_stdout(buf):
            try:
                fld = CMField(q, e[0], mu, nu)
            except ValueError as ex:
                return [dict(base, cmd="sweep", status=str(ex))], buf.getvalue()
            ns = argparse.Namespace(q=q, e1=e[0], mu=mu, nu=nu, primes=nprimes, all_subsets=all_subsets)
            tor = None
            for digits in digits_list:
                ns.digits = digits
                t0 = time.time()
                try:
                    tor = Torsor(fld, digits, verbose=False)
                except ValueError as ex:
                    return [dict(base, cmd="sweep", status=str(ex))], buf.getvalue()
                if tor.N > maxn:
                    return [dict(base, cmd="sweep", status=f"skipped: {tor.N} lattices per type", N=tor.N)], buf.getvalue()
                ok = tor.setup(q, e[0] % q)
                print(f"digits {digits}: {tor.N} lattices per type, class polynomial {'found' if ok else 'not found'}, "
                      f"{time.time() - t0:.0f} s", flush=True)
                if ok:
                    break
            if tor is None or tor.mix is None:
                return [dict(base, cmd="sweep", status="class polynomial not obtained", N=tor.N if tor else None,
                             digits=digits_list[-1])], buf.getvalue()
            r1 = cmd_compare(ns, tor)
            r2 = cmd_action(ns, tor)
            out = [dict(base, **r1), dict(base, **r2)]
    except Exception as ex:   # noqa: BLE001  (a failure of one field must not stop the sweep)
        import traceback
        out.append(dict(base, cmd="sweep", status="error: " + repr(ex)))
        buf.write(traceback.format_exc())
    return out, buf.getvalue()


def cmd_sweep(a):
    from multiprocessing import Pool
    jobs = []
    with open(os.path.join(HERE, "sixlines_cm.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            if r["p"] in a.p and r["index"] == 1 and (not a.select or f"{r['p']},{r['e'][0]}" in a.select):
                jobs.append((r["p"], r["e"], a.digits_list, a.primes, a.all_subsets, a.maxn))
    print(f"{len(jobs)} polynomials with maximal order", flush=True)
    with Pool(a.workers) as pool:
        for recs, log in pool.imap_unordered(sweep_one, jobs):
            print("=" * 100)
            print(f"q = {recs[0]['q']}, e = {recs[0]['e']}")
            print(log, end="")
            for r in recs:
                if r.get("cmd") == "sweep":
                    print("status:", r["status"])
            sys.stdout.flush()
            if a.out:
                with open(a.out, "a") as fh:
                    for r in recs:
                        fh.write(json.dumps(r) + "\n")


# ------------------------------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("selftest")
    s.add_argument("--digits", type=int, default=60)
    for name in ("points", "compare", "action"):
        s = sub.add_parser(name)
        s.add_argument("q", type=int)
        s.add_argument("e1", type=int)
        s.add_argument("mu", type=int)
        s.add_argument("nu", type=int)
        s.add_argument("--digits", type=int, default=100)
        s.add_argument("--maximal", action="store_true", help="lattices of the maximal order only")
        s.add_argument("--nolevel", action="store_true", help="do not ask for (pi - 1)/(1 + i) in the order")
        s.add_argument("--out", help="append a JSON record to this file")
        s.add_argument("--primes", type=int, default=6, help="number of primes l (action)")
        s.add_argument("--all-subsets", action="store_true", help="all sixteen ideals above each l (action)")
    s = sub.add_parser("sweep")
    s.add_argument("p", type=int, nargs="+")
    s.add_argument("--digits-list", type=int, nargs="+", default=[200, 500, 1000])
    s.add_argument("--primes", type=int, default=3)
    s.add_argument("--all-subsets", action="store_true")
    s.add_argument("--maxn", type=int, default=12, help="skip fields with more lattices per type")
    s.add_argument("--workers", type=int, default=4)
    s.add_argument("--select", nargs="*", help="only the polynomials q,e1 listed")
    s.add_argument("--out")
    a = ap.parse_args()
    if a.cmd == "sweep":
        cmd_sweep(a)
        return
    if a.cmd == "selftest":
        sys.exit(0 if selftest(a.digits) else 1)
    if a.cmd == "points":
        cmd_points(a)
    if a.cmd in ("compare", "action"):
        rec = cmd_compare(a) if a.cmd == "compare" else cmd_action(a)
        if a.out:
            with open(a.out, "a") as fh:
                fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
