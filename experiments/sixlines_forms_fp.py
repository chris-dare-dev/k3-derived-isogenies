#!/usr/bin/env python
"""
sixlines_forms_fp.py -- which fourfold over F_p belongs to a six-line surface of sign - or with reducible P (note 21).

It uses the functions of sixlines_prym.py unchanged.  Notation: X is the double plane w^2 = l1 ... l6 over F_p,
P(t) = t^6 - e1 t^5 + e2 t^4 - e3 t^3 + ... the polynomial of Frobenius on its part of rank 6, chi the quartic of
Frobenius on the i-eigenspace of the curve of note 14 for an ordered pair of lines, chi_2 the quartic of Frobenius^2.

The rules that are tested (derived in note 21):

  sign -.   P = (t^2 - p^2)(t^2 - g t + p^2)(t^2 - g' t + p^2), and
            h(X) = X^4 - A X^3 + B X^2 - p^2 A X + p^4,   A = g g' / p,   B = g^2 + g'^2 - 2 p^2,
            is the quartic of Frobenius^2 on the i-eigenspace of the fourfold of X; the fourfold over F_p is
            isogenous to the Weil restriction of an abelian surface over F_{p^2} with polynomial h, and its
            polynomial over F_p is h(X^2).
            (a) pair of inverting type:   chi_2 = h.
            (b) pair of commuting type:   chi(X) = -(1/4) r((1 - eps i) X) with eps = +-1 and r in Z[X],
                r(X) r(-X) = 16 h(X^2 / 2).
  sign +, pair of inverting type.   P(t) or P(-t) is the second exterior power of a quartic
            f = X^4 - a X^3 + b X^2 - p a X + p^2 in Z[X], and chi_2 is the quartic f2 of the squares of the
            roots of f.
  A pair is called anomalous if its type is not the type of the sign (commuting on sign -, inverting on sign +).
  Anomalous pairs occur only for a pair {e, f} of lines such that, for a splitting of the other four in two
  pairs {a, b}, {c, d}, the points La n Lb, Lc n Ld, Le n Lf are collinear ("the pair lies in a collinear triple").

Subcommands:

  commuting FILE      The stored pairs of commuting type on surfaces with sign - (FILE is an output file of
                      sixlines_prym.py).  Rule (b), from the surface alone: the set of r with
                      r(X) r(-X) = 16 h(X^2/2) is computed from (e1, e2) and chi is compared with the predicted
                      quartics.  Also the reverse: h from chi, then P from h, compared with P and P(-t).

  types p [p ...]     For every configuration of FILE (default sixlines_prym.jsonl) with these p and every ordered
                      pair: the type by Lemma 4 of note 14, whether b1, b2 are rational, whether the pair lies
                      in a collinear triple, whether one of 0, 1, oo is a branch value of the cross-ratio map.
                      For every configuration and each of the 15 splittings in three pairs: collinear or not,
                      and whether exchanging the lines of each pair gives the associated configuration.
                      Tests of the rules on the polynomials of all surfaces with a collinear triple.

  inverting p [p ...] Pairs of inverting type: chi_2 from the character sums over F_{p^2}, F_{p^4} (and F_{p^6}
                      with --k 3) of the curve v^4 = F, which is defined over F_{p^2}.  For p = 1 mod 4 the
                      point c4 is not rational; pr_trace of sixlines_prym.py sets the square root om to 0 when
                      a, b, c are rational, so a copy of it with om always set is used.  --pairs n: at most n
                      pairs per configuration (0 = all), --limit m: the first m configurations with such a pair.
                      Cost: p^4 points per pair with --k 2, p^6 with --k 3.

  involution p [p ...]  Lemma E of note 21, directly.  For every pair that lies in a collinear triple: psi(x) = t - x
                      with t the sum of rho and the point of order 2 among a, b, c; F(x) F(psi(x)) is evaluated at
                      random points of E(F_{p^2}) and its class modulo fourth powers must be constant.  Control: the
                      same with one of the other two points in place of the point of order 2.

  selftest            The linear algebra behind the rule for sign -, numerically (numpy): for random a, a' of
                      absolute value p and a random hermitian form of signature (2, 2), an antilinear map phi of
                      C^4 with H(phi x, phi y) = -p conj H(x, y) and phi^2 with eigenvalues a, conj a, a', conj a'
                      is built, the real form T of the second exterior power that it preserves is computed, and
                      the polynomial of phi on T is compared with the rule.

  python -I sixlines_forms_fp.py selftest
  python -I sixlines_forms_fp.py involution 11 13 17 19 23 --workers 4
  python -I sixlines_forms_fp.py commuting sixlines_prym.jsonl
  python -I sixlines_forms_fp.py types 5 7 11 13 17 19 23 29 31 37 --workers 4
  python -I sixlines_forms_fp.py inverting 7 11 13 17 19 --pairs 0 --k 2 --workers 4 --out sixlines_forms_fp.jsonl
  python -I sixlines_forms_fp.py inverting 13 --pairs 1 --k 3 --limit 12 --workers 4

Run from experiments/.
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
import sixlines_prym as sp  # noqa: E402

pari = sp.pari
HERE = os.path.dirname(os.path.abspath(__file__))

# the quartic whose roots are the squares of the roots of f, and the second exterior power of a quartic
pari("ff_squares(f) = my(g = f * subst(f, 'X, -'X)); subst(substpol(g, 'X^2, 'Y), 'Y, 'X)")
pari("ff_wedge(f) = my(M = matcompanion(f), idx = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]], W); "
     "W = matrix(6, 6, r, s, M[idx[r][1], idx[s][1]] * M[idx[r][2], idx[s][2]] - M[idx[r][1], idx[s][2]] * M[idx[r][2], idx[s][1]]); "
     "charpoly(W, 'X)")
pari("ff_pattern(f) = my(F = factor(f)); vecsort(concat(vector(#F~, j, vector(F[j, 2], k, poldegree(F[j, 1])))))")


def setup_trace():
    """fp_trace: pr_trace of sixlines_prym.py with the square root om set also when a, b, c are rational."""
    src = open(os.path.join(HERE, "sixlines_prym.py")).read()
    a = src.index("pr_trace(W, bx, rat, Qc, Xc, p, k, modeB) = {")
    b = src.index('""")', a)
    body = src[a:b].replace("pr_trace(", "fp_trace(", 1)
    old = "om = if(rat, sqrt(n0 * one), 0));"
    assert body.count(old) == 1
    pari(body.replace(old, "om = sqrt(n0 * one));"))


# ---------------------------------------------------------------- polynomials of the surface

def surface_pol(p, e, sign):
    e1, e2, e3 = e
    return pari(f"X^6-({e1})*X^5+({e2})*X^4-({e3})*X^3+({sign})*(({p * p * e2})*X^2-({p ** 4 * e1})*X+{p ** 6})")


def h_of_surface(p, e):
    """(A, B) of h for a surface with sign -: g g' = e2 - p^2, g + g' = e1.  None if p does not divide g g'."""
    e1, e2, e3 = e
    assert e3 == 0
    gg = e2 - p * p
    if gg % p:
        return None
    return gg // p, e1 * e1 - 2 * gg - 2 * p * p


def hpol(p, A, B):
    return pari(f"X^4-({A})*X^3+({B})*X^2-({p * p * A})*X+{p ** 4}")


def surface_of_h(p, A, B):
    """(e1, e2) up to the sign of e1 for the polynomial of sign - that belongs to h, or None."""
    s = B + 2 * p * A + 2 * p * p
    if s < 0 or math.isqrt(s) ** 2 != s:
        return None
    return math.isqrt(s), p * A + p * p


def roots_r(p, A, B, scale):
    """All r = X^4 - a X^3 + b X^2 - scale p a X + scale^2 p^2 in Z[X] with r(X) r(-X) = scale^4 h(X^2 / scale).

    scale = 1: h is the quartic of the squares of the roots of r (Frobenius^2 is a square in the algebra of h).
    scale = 2: the roots of h are the squares of the roots of r, divided by 2."""
    out = []
    h = hpol(p, A, B)
    target = pari(f"{scale ** 4}") * pari.subst(h, "X", pari(f"X^2/{scale}"))
    bound = int(4 * math.sqrt(scale * p)) + 1
    for a in range(-bound, bound + 1):
        if (a * a - scale * A) % 2:
            continue
        b = (a * a - scale * A) // 2
        r = pari(f"X^4-({a})*X^3+({b})*X^2-({scale * p * a})*X+{scale * scale * p * p}")
        if r * pari.subst(r, "X", pari("-X")) == target:
            out.append(r)
    return out


def weil_quartics(p, e):
    """All (kappa, f): f = X^4 - a X^3 + b X^2 - p a X + p^2 in Z[X] whose second exterior power is P(kappa t)."""
    out = []
    P = surface_pol(p, e, 1)
    bound = int(4 * math.sqrt(p)) + 1
    for kappa in (1, -1):
        Pk = pari.subst(P, "X", pari(f"{kappa}*X"))
        b = kappa * e[0]                      # the coefficient of t^5 of the second exterior power is -b
        for a in range(-bound, bound + 1):
            f = pari(f"X^4-({a})*X^3+({b})*X^2-({p * a})*X+{p * p}")
            if pari("ff_wedge")(f) == Pk:
                out.append((kappa, f))
    return out


# ---------------------------------------------------------------- commuting

def run_commuting(path):
    st = Counter()
    exceptions = []
    for line in open(path):
        rec = json.loads(line)
        if rec["signs"] != [-1]:
            continue
        p, e = rec["p"], rec["e"]
        ab = h_of_surface(p, e)
        had = False
        L = sp.forms(rec["cfg"])
        duads = Counter()
        for syn in SYNTHEMES:
            pts = [cross(L[i - 1], L[j - 1], p) for i, j in syn]
            if det3(pts[0], pts[1], pts[2], p) == 0:
                duads.update(syn)
        for x in rec["pairs"]:
            if "c" not in x:
                continue
            had = True
            st["pairs"] += 1
            if ab is None:
                st["p does not divide g g'"] += 1
                exceptions.append((p, rec["cfg"], x["pair"], "p does not divide g g'"))
                continue
            A, B = ab
            h = hpol(p, A, B)
            c1, c2, c3 = x["c"]
            st["c3 = -p conj(c1)"] += c3 == [-p * c1[0], p * c1[1]]
            st["c2 purely imaginary"] += c2[0] == 0
            st["|Re c1| = |Im c1|"] += abs(c1[0]) == abs(c1[1])
            st["a, b, c not all rational"] += x["rational_point"] != 0
            zero = c1 == [0, 0]
            st["c1 = 0"] += zero
            st[f"|Re c1| = {abs(c1[0])}"] += 1
            chi = pari(f"X^4-(({c1[0]})+I*({c1[1]}))*X^3+(({c2[0]})+I*({c2[1]}))*X^2-(({c3[0]})+I*({c3[1]}))*X-{p * p}")
            # r_eps = -4 chi(X / (1 - eps i)); it has to be in Z[X]
            good = []
            for eps in (1, -1):
                r = pari(f"-4") * pari.subst(chi, "X", pari(f"X/(1-({eps})*I)"))
                if pari.imag(r) != 0 or pari.denominator(pari.content(r)) != 1:
                    continue
                st["r_eps in Z[X] (pairs x eps)"] += 1
                ok = r * pari.subst(r, "X", pari("-X")) == 16 * pari.subst(h, "X", pari("X^2/2"))
                good.append((eps, r, bool(ok)))
            agree = [g for g in good if g[2]]
            key = "c1 = 0" if zero else "c1 != 0"
            st[f"{key}: pairs"] += 1
            st[f"{key}: pairs in {duads[tuple(sorted(x['pair']))]} collinear triple(s)"] += 1
            st[f"{key}: r in Z[X] for {len(good)} of the two eps"] += 1
            st[f"{key}: r(X) r(-X) = 16 h(X^2/2) for {len(agree)} eps"] += 1
            if not agree:
                exceptions.append((p, rec["cfg"], x["pair"], "no eps agrees"))
                continue
            # the prediction from the surface alone: chi is one of -(1/4) r((1 - eps i) X), r from (e1, e2)
            pred = set()
            for r in roots_r(p, A, B, 2):
                for eps in (1, -1):
                    pred.add(str(pari("-1/4") * pari.subst(r, "X", pari(f"(1-({eps})*I)*X"))))
            st[f"{key}: chi is one of the predicted quartics"] += str(chi) in pred
            st[f"{key}: number of predicted quartics = {len(pred)}"] += 1
            if str(chi) not in pred:
                exceptions.append((p, rec["cfg"], x["pair"], "chi not among the predicted quartics"))
            # the reverse: h from the curve, then the polynomial of sign - that belongs to it
            res = set()
            for eps, r, ok in good:
                hc = pari.subst(pari("ff_squares")(r), "X", pari("2*X")) / 16
                Ac, Bc = -int(pari.polcoef(hc, 3)), int(pari.polcoef(hc, 2))
                s = surface_of_h(p, Ac, Bc)
                if s is None:
                    res.add("no")
                elif (s[0], s[1]) == (e[0], e[1]) or (-s[0], s[1]) == (e[0], e[1]):
                    res.add("P or P(-t)")
                else:
                    res.add("another polynomial")
            st[f"{key}: from the curve, polynomial of sign -: {sorted(res)}"] += 1
            # structure
            if not zero:
                eps, r, ok = agree[0]
                st[f"c1 != 0: h factors with degrees {pari('ff_pattern')(h)}"] += 1
                st[f"c1 != 0: chi conj(chi) factors over Q with degrees {pari('ff_pattern')(pari.real(chi * pari.conj(chi)))}"] += 1
                chi2 = pari("ff_squares")(chi)
                st[f"c1 != 0: chi_2 conj(chi_2) factors over Q with degrees {pari('ff_pattern')(pari.real(chi2 * pari.conj(chi2)))}"] += 1
                chi4 = pari("ff_squares")(chi2)
                st[f"c1 != 0: chi_4 in Q[X] = {pari.imag(chi4) == 0}, with degrees {pari('ff_pattern')(pari.real(chi4))}"] += 1
                # which of the surface and its twist: delta = (sum of sigma_j sigma_k / 2 over the two pairs) / e1
                m2 = int(pari.polcoef(r, 2)) - 4 * p
                if e[0] == 0:
                    st["c1 != 0: e1 = 0"] += 1
                else:
                    st[f"c1 != 0: delta = {pari(m2) / (2 * e[0])}"] += 1
        st["surfaces with such pairs"] += had
    print("commuting pairs on surfaces with sign -:")
    for k in sorted(st):
        print(f"  {k}: {st[k]}")
    print(f"  exceptions: {len(exceptions)}")
    for x in exceptions:
        print("   ", x)


# ---------------------------------------------------------------- types

def det3(a, b, c, p):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) % p


def cross(a, b, p):
    return ((a[1] * b[2] - a[2] * b[1]) % p, (a[2] * b[0] - a[0] * b[2]) % p, (a[0] * b[1] - a[1] * b[0]) % p)


def normal_form(L, p):
    """(b, c, e, f) of six labelled lines: the first three to the coordinate lines, the fourth to (1, 1, 1)."""
    if det3(L[0], L[1], L[2], p) == 0:
        return None
    co = lambda v: (det3(v, L[1], L[2], p), det3(L[0], v, L[2], p), det3(L[0], L[1], v, p))  # noqa: E731
    c4 = co(L[3])
    if 0 in c4:
        return None
    inv = [pow(x, p - 2, p) for x in c4]
    out = []
    for j in (4, 5):
        c = [x * i % p for x, i in zip(co(L[j]), inv)]
        if c[0] == 0:
            return None
        s = pow(c[0], p - 2, p)
        out += [c[1] * s % p, c[2] * s % p]
    return tuple(out)


def synthemes():
    out = []
    for a in range(2, 7):
        rest = [x for x in range(2, 7) if x != a]
        for b in rest[1:]:
            cd = [x for x in rest[1:] if x != b]
            out.append(((1, a), (rest[0], b), tuple(cd)))
    return out


SYNTHEMES = synthemes()


def types_one(job):
    p, cfg, sign, degrees, e = job
    L = sp.forms(cfg)
    b, c, ee, f = cfg
    # the associated configuration with its labels (Gale transform of the matrix (1 | A))
    gale = normal_form([(1, 1, 1), (1, b, ee), (1, c, f), (1, 0, 0), (0, 1, 0), (0, 0, 1)], p)
    coll, assoc = [], []
    for syn in SYNTHEMES:
        pts = [cross(L[i - 1], L[j - 1], p) for i, j in syn]
        coll.append(det3(pts[0], pts[1], pts[2], p) == 0)
        perm = list(range(6))
        for i, j in syn:
            perm[i - 1], perm[j - 1] = j - 1, i - 1
        assoc.append(normal_form([L[k] for k in perm], p) == gale)
    duads = Counter(d for syn, cl in zip(SYNTHEMES, coll) if cl for d in syn)
    Lm = pari.matrix(6, 3, [x for row in L for x in row])
    pairs = []
    for pair in itertools.permutations(range(1, 7), 2):
        rest = [i for i in range(1, 7) if i not in pair]
        N, D = sp.GP["pr_lambda"](Lm, rest + list(pair), p)
        cur = sp.GP["pr_curve"](N, D, p, 0)
        induad = duads[tuple(sorted(pair))]          # the number of collinear triples that contain the pair: 0, 1 or 2
        if int(cur[0]) == 0:
            pairs.append((pair, None, None, induad, None))
            continue
        W, bx = cur[1], cur[2]
        yrat, brat = bool(int(cur[5])), bool(int(cur[6]))
        special = any(pari.subst(W, "x", bx[j]) == 0 for j in range(3))
        pairs.append((pair, yrat == (p % 4 == 1), brat, induad, special))
    return p, cfg, sign, degrees, e, coll, assoc, pairs


def classify(sign, degrees):
    if sign == -1:
        return "sign -"
    return "sign +, P irreducible" if degrees == [6] else "sign +, P reducible"


def run_types(primes, path, workers):
    jobs = []
    for line in open(path):
        r = json.loads(line)
        if r["p"] in primes:
            assert len(r["signs"]) == 1
            jobs.append((r["p"], r["cfg"], r["signs"][0], r.get("degrees"), r["e"]))
    with Pool(workers) as pool:
        out = pool.map(types_one, jobs, chunksize=16)
    st = Counter()
    perp = {}
    exceptions = []
    for p, cfg, sign, degrees, e, coll, assoc, pairs in out:
        kind = classify(sign, degrees)
        pp = perp.setdefault(p, Counter())
        st["configurations"] += 1
        pp["configurations"] += 1
        pp[kind] += 1
        for cl, asc in zip(coll, assoc):
            st[f"splittings: collinear = {cl}, exchange gives the associate = {asc}"] += 1
            pp["collinear splittings"] += cl
        has = any(coll)
        st[f"{kind}: configurations"] += 1
        st[f"{kind}: with a collinear triple"] += has
        pp[f"{kind}: with a collinear triple"] += has
        pp["with a collinear triple"] += has
        ncomm = ninv = 0
        for pair, comm, brat, induad, special in pairs:
            if comm is None:
                st[f"{kind}: pairs degenerate, in a collinear triple = {induad > 0}"] += 1
                continue
            anomalous = comm != (sign == 1)
            ncomm += comm
            ninv += not comm
            st[f"{kind}: pairs {'anomalous' if anomalous else 'of the type of the sign'}, in a collinear triple = {induad > 0}"] += 1
            st[f"{kind}: pairs {'anomalous' if anomalous else 'of the type of the sign'}, in {induad} collinear triples"] += 1
            st[f"pairs in a collinear triple = {induad > 0}: b1, b2 rational = {brat}"] += 1
            st[f"pairs in a collinear triple = {induad > 0}: 0, 1 or oo is a branch value = {special}"] += 1
            st[f"pairs {'anomalous' if anomalous else 'of the type of the sign'}: b1, b2 rational = {brat}"] += 1
            if anomalous and not induad:
                exceptions.append((p, cfg, pair, "anomalous pair outside a collinear triple"))
        # the rules on the polynomial of the surface
        if has and sign == -1:
            ab = h_of_surface(p, e)
            if ab is None:
                st["sign -, collinear triple: p does not divide g g'"] += 1
                exceptions.append((p, cfg, None, "p does not divide g g'"))
            else:
                inner, outer = bool(roots_r(p, ab[0], ab[1], 1)), bool(roots_r(p, ab[0], ab[1], 2))
                st[f"sign -, collinear triple: Frobenius^2 a square = {inner}, 2 Frobenius^2 a square = {outer}"] += 1
                anom = sum(1 for pair, comm, brat, induad, special in pairs if comm)
                st[f"sign -, collinear triple, 2 Frobenius^2 a square = {outer}: configurations with a commuting pair"] += anom > 0
                if not inner and not outer:
                    exceptions.append((p, cfg, None, "sign -, collinear triple: neither Frobenius^2 nor 2 Frobenius^2 is a square"))
                if anom and not outer:
                    exceptions.append((p, cfg, None, "commuting pair but 2 Frobenius^2 is not a square"))
        if sign == -1 and not has:
            ab = h_of_surface(p, e)
            if ab is not None:
                inner, outer = bool(roots_r(p, ab[0], ab[1], 1)), bool(roots_r(p, ab[0], ab[1], 2))
                st[f"sign -, no collinear triple: Frobenius^2 a square = {inner}, 2 Frobenius^2 a square = {outer}"] += 1
        if sign == 1:
            mn = sp.mu_nu(p, e)
            st[f"sign +: P reducible = {degrees != [6]}, mu nu = 0: {None if mn is None else mn[0] * mn[1] == 0}"] += 1
            if not has:
                pp["sign +, no collinear triple"] += 1
                pp["sign +, no collinear triple, P reducible"] += degrees != [6]
                if p == max(primes) and mn is not None:
                    st[f"sign +, no collinear triple, p = {p}: (|mu|, |nu|) up to order = {tuple(sorted(mn))}"] += 1
            fs = weil_quartics(p, e)
            st[f"{kind}, collinear triple = {has}: P(+-t) is the second exterior power of a quartic in Z[X] = {bool(fs)}"] += 1
            if has and not fs:
                exceptions.append((p, cfg, None, "sign +, collinear triple: P is not a second exterior power over Z"))
    print("types:")
    for k in sorted(st):
        print(f"  {k}: {st[k]}")
    print("  by p (configurations; with a collinear triple; 1 - (1 - 1/p)^15; then by class: all / with a collinear triple):")
    for p in sorted(perp):
        pp = perp[p]
        n = pp["configurations"]
        rate = pp["collinear splittings"] / (15 * n)
        row = (f"    p = {p}: {n}; {pp['with a collinear triple']} = {100 * pp['with a collinear triple'] / n:.0f}%; "
               f"{100 * (1 - (1 - 1 / p) ** 15):.0f}%; collinear splittings {pp['collinear splittings']}, p x rate = {p * rate:.2f}, "
               f"1 - (1 - rate)^15 = {100 * (1 - (1 - rate) ** 15):.0f}%;")
        for kind in ("sign +, P irreducible", "sign +, P reducible", "sign -"):
            row += f"  {kind}: {pp[kind]} / {pp[kind + ': with a collinear triple']}"
        row += f"  sign +, no collinear triple: {pp['sign +, no collinear triple']}, of which P reducible {pp['sign +, no collinear triple, P reducible']}"
        print(row)
    print(f"  exceptions: {len(exceptions)}")
    for x in exceptions:
        print("   ", x)
    sys.stdout.flush()


# ---------------------------------------------------------------- inverting

def inverting_job(job):
    p, cfg, sign, e, pair, induad, kmax = job
    t0 = time.time()
    rest = [i for i in range(1, 7) if i not in pair]
    L = pari.matrix(6, 3, [x for row in sp.forms(cfg) for x in row])
    N, D = sp.GP["pr_lambda"](L, rest + list(pair), p)
    cur = sp.GP["pr_curve"](N, D, p, 0)
    W, bx, rat, Xc = cur[1], cur[2], int(cur[3]), cur[4]
    Qc = sp.GP["pr_findq"](W, bx, rat, Xc, p) if rat else pari("[0, 0, 0, 0]")
    T = []
    for k in range(1, kmax + 1):
        cnt = pari("fp_trace")(W, bx, rat, Qc, Xc, p, 2 * k, 0)
        if cnt == 0:
            return None
        n = [int(v) for v in cnt]
        T.append((-(n[0] - n[2]), -(n[1] - n[3])))
    rec = {"p": p, "cfg": cfg, "sign": sign, "e": e, "pair": list(pair), "in_collinear_triple": int(induad),
           "rational_point": rat, "T": T, "seconds": round(time.time() - t0, 1)}
    rec["rational"] = all(t[1] == 0 for t in T)
    if rec["rational"]:
        t1, t2 = T[0][0], T[1][0]
        c1, c2 = t1, (t1 * t1 - t2) // 2
        rec["c"] = [c1, c2]
        if kmax >= 3:
            six = t1 ** 3 - 3 * t1 * t2 + 2 * T[2][0]
            rec["c3_is_p2_c1"] = six == 6 * p * p * c1
        if sign == -1:
            ab = h_of_surface(p, e)
            rec["predicted"] = [list(ab)] if ab else []
        else:
            pred = set()
            for kappa, f in weil_quartics(p, e):
                f2 = pari("ff_squares")(f)
                pred.add((-int(pari.polcoef(f2, 3)), int(pari.polcoef(f2, 2))))
            rec["predicted"] = sorted(list(x) for x in pred)
        # the factorisation of chi_2 over Q: degrees of the irreducible factors, with multiplicity
        chi2 = hpol(p, c1, c2)
        rec["chi2_degrees"] = [int(d) for d in pari("ff_pattern")(chi2)]
        rec["chi2_squarefree"] = bool(pari.issquarefree(chi2))
        rec["agrees"] = [c1, c2] in rec["predicted"]
        rec["agrees_up_to_sign"] = any(abs(c1) == abs(x[0]) and c2 == x[1] for x in rec["predicted"])
    return rec


def inverting_list(job):
    """The pairs of inverting type of one configuration that are not degenerate before any count."""
    p, cfg, sign, degrees, e = job
    p, cfg, sign, degrees, e, coll, assoc, pairs = types_one(job)
    return [(p, cfg, sign, e, pair, induad) for pair, comm, brat, induad, special in pairs if comm is False]


def run_inverting(primes, path, max_pairs, limit, kmax, workers, out):
    fh = open(os.path.join(HERE, out), "a") if out else None
    total = Counter()
    for p in primes:
        t0 = time.time()
        cfgs = []
        for line in open(path):
            r = json.loads(line)
            if r["p"] == p:
                cfgs.append((p, r["cfg"], r["signs"][0], r.get("degrees"), r["e"]))
        with Pool(workers, initializer=setup_trace) as pool:
            lists = [x for x in pool.map(inverting_list, cfgs, chunksize=16) if x]
            if limit:
                lists = lists[:limit]
            jobs = []
            for lst in lists:
                # pairs in a collinear triple first, so that both kinds are present when max_pairs is small
                lst = sorted(lst, key=lambda x: (x[5] != 1, x[4]))
                if max_pairs:
                    half = max(1, max_pairs // 2)
                    a = [x for x in lst if x[5] == 1][:half]
                    b = [x for x in lst if x[5] != 1][:max_pairs - len(a)]
                    lst = a + b if len(a) + len(b) >= min(max_pairs, len(lst)) else lst[:max_pairs]
                jobs += [x + (kmax,) for x in lst]
            recs = pool.map(inverting_job, jobs, chunksize=1)
        st = Counter()
        st["configurations with a pair of inverting type"] = len(lists)
        for rec in recs:
            if rec is None:
                st["pairs degenerate"] += 1
                continue
            kind = "sign -" if rec["sign"] == -1 else "sign +"
            kind += f", in {rec['in_collinear_triple']} collinear triple(s)"
            st[f"{kind}: pairs"] += 1
            st[f"{kind}: chi_2 rational"] += rec["rational"]
            if rec["rational"]:
                sq = "chi_2 squarefree" if rec["chi2_squarefree"] else "chi_2 with a repeated factor"
                st[f"{kind}, {sq}: pairs"] += 1
                st[f"{kind}, {sq}: chi_2 as predicted"] += rec["agrees"]
                st[f"{kind}: chi_2 as predicted"] += rec["agrees"]
                st[f"{kind}: as predicted up to the sign of the odd coefficients only"] += rec["agrees_up_to_sign"] and not rec["agrees"]
                st[f"{kind}: number of predicted quartics = {len(rec['predicted'])}"] += 1
                if "c3_is_p2_c1" in rec:
                    st[f"{kind}: c3 = p^2 c1 (third trace)"] += rec["c3_is_p2_c1"]
            if not rec.get("agrees"):
                print(f"  exception: {json.dumps(rec)}", flush=True)
            if fh:
                fh.write(json.dumps(rec) + "\n")
        if fh:
            fh.flush()
        print(f"inverting, p = {p}, k = {kmax}, pairs per configuration = {max_pairs or 'all'}"
              f"{', first ' + str(limit) + ' configurations' if limit else ''}: {dict(sorted(st.items()))}  ({time.time() - t0:.0f} s)", flush=True)
        total += st
    if len(primes) > 1:
        print(f"inverting, total: {dict(sorted(total.items()))}", flush=True)


# ---------------------------------------------------------------- involution

pari(r"""
fp_invol(W, bx, rat, Qc, Xc, p, j) = {
  \\ [n, m]: at n random points x of E(F_{p^2}) the value F(x) F(t - x), t = (point number j of a, b, c) + rho, was
  \\ compared with its value at a first point; m of the ratios are fourth powers.  F as in pr_trace.
  my(q = p^2, g = ffgen(q, 'g), one = g^0, z = ffprimroot(g),
     a2 = polcoef(W, 2) * one, a4 = polcoef(W, 1) * one, a6 = polcoef(W, 0) * one, E = ellinit([0, a2, 0, a4, a6]),
     n0, om, r, c4, Q = [0], rho, B4, cs, zeros, D, ch, t, val, v0 = 0, ok = 0, tot = 0, X, rhs, P1, P2, f1, f2);
  if(p % 4 == 3, n0 = Mod(-1, p); om = sqrt(-one), n0 = Mod(znprimroot(p), p); om = sqrt(n0 * one));
  r = vector(3, jj, pr_point(W, bx[jj], one, om, n0));
  c4 = pr_point(W, Xc, one, om, n0);
  if(rat == 0,
    rho = ellneg(E, elladd(E, elladd(E, r[1], r[2]), r[3])),
    Q = [Qc[1] * one + Qc[2] * om, Qc[3] * one + Qc[4] * om];
    rho = ellsub(E, elladd(E, Q, [Q[1]^p, Q[2]^p]), r[rat]));
  B4 = concat([Q], vector(3, jj, ellsub(E, elladd(E, r[jj], rho), Q)));
  cs = elladd(E, c4, rho);
  for(i = 1, 4, if(B4[i] == cs || B4[i] == elladd(E, ellneg(E, c4), rho), return([0, 0])); for(k = i + 1, 4, if(B4[i] == B4[k], return([0, 0]))));
  if(cs == [0], return([0, 0]));
  zeros = select(P -> P != [0], B4);
  D = concat(vector(#zeros, jj, [zeros[jj], 2]), [[cs, -4]]);
  ch = pr_chain(E, D);
  t = elladd(E, r[j], rho);
  for(n = 1, 300,
    X = z^random(q - 1); rhs = ((X + a2) * X + a4) * X + a6;
    if(rhs == 0 || rhs^((q - 1) / 2) != one, next);
    P1 = [X, sqrt(rhs)]; P2 = ellsub(E, t, P1);
    if(P2 == [0], next);
    f1 = pr_eval(ch, P1[1], P1[2]); f2 = pr_eval(ch, P2[1], P2[2]);
    if(f1[1] == 0 || f1[2] == 0 || f2[1] == 0 || f2[2] == 0, next);
    val = f1[1] / f1[2] * f2[1] / f2[2];
    if(v0 == 0, v0 = val; next);
    tot++; if((val / v0)^((q - 1) / 4) == one, ok++));
  [tot, ok];
}
""")


def involution_job(job):
    p, cfg, sign, degrees, e = job
    L = sp.forms(cfg)
    Lm = pari.matrix(6, 3, [x for row in L for x in row])
    out = []
    for pair in itertools.permutations(range(1, 7), 2):
        rest = [i for i in range(1, 7) if i not in pair]
        N, D = sp.GP["pr_lambda"](Lm, rest + list(pair), p)
        cur = sp.GP["pr_curve"](N, D, p, 0)
        if int(cur[0]) == 0:
            continue
        W, bx, rat, Xc = cur[1], cur[2], int(cur[3]), cur[4]
        special = [j + 1 for j in range(3) if pari.subst(W, "x", bx[j]) == 0]
        if not special:
            continue
        Qc = sp.GP["pr_findq"](W, bx, rat, Xc, p) if rat else pari("[0, 0, 0, 0]")
        for j in (1, 2, 3):
            n, m = (int(v) for v in pari("fp_invol")(W, bx, rat, Qc, Xc, p, j))
            out.append((j in special, n, m))
    return out


def run_involution(primes, path, workers):
    for p in primes:
        jobs = []
        for line in open(path):
            r = json.loads(line)
            if r["p"] == p:
                jobs.append((p, r["cfg"], r["signs"][0], r.get("degrees"), r["e"]))
        with Pool(workers) as pool:
            res = [x for lst in pool.map(involution_job, jobs, chunksize=8) for x in lst]
        st = Counter()
        for special, n, m in res:
            if n == 0:
                st["degenerate"] += 1
                continue
            key = "point of order 2" if special else "control (another of a, b, c)"
            st[f"{key}: pairs"] += 1
            st[f"{key}: pairs with a constant class modulo fourth powers"] += n == m
            st[f"{key}: points"] += n
            st[f"{key}: points in the class of the first"] += m
        print(f"involution, p = {p}: {dict(sorted(st.items()))}", flush=True)


# ---------------------------------------------------------------- selftest

def run_selftest(trials=200, seed=1):
    import numpy as np
    rng = np.random.default_rng(seed)
    idx = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    def wedge(M):
        return np.array([[M[i, k] * M[j, l] - M[i, l] * M[j, k] for (k, l) in idx] for (i, j) in idx])

    def sgn(perm):
        s = 1
        for x in range(4):
            for y in range(x + 1, 4):
                if perm[x] > perm[y]:
                    s = -s
        return s

    worst = 0.0
    for _ in range(trials):
        p = float(rng.choice([5, 7, 11, 13, 17]))
        a, a2 = p * np.exp(1j * rng.uniform(0.1, 3.0)), p * np.exp(1j * rng.uniform(0.1, 3.0))
        h1, h2 = rng.uniform(0.5, 2) * rng.choice([-1, 1]), rng.uniform(0.5, 2) * rng.choice([-1, 1])
        # basis v, w = phi v, v', w' = phi v';  phi(x) = M conj(x);  H = diag(h1, -p h1, h2, -p h2)
        M = np.zeros((4, 4), dtype=complex)
        M[1, 0], M[0, 1], M[3, 2], M[2, 3] = 1, a, 1, a2
        h = [h1, -p * h1, h2, -p * h2]
        # the star: omega' ^ star(omega) = H(omega, omega') vol;  star(omega) = S conj(omega)
        S = np.zeros((6, 6), dtype=complex)
        for c, (i, j) in enumerate(idx):
            k, l = [x for x in range(4) if x not in (i, j)]
            S[idx.index((k, l)), c] = h[i] * h[j] * sgn((i, j, k, l)) / (p * abs(h1 * h2))
        W = wedge(M)
        A1, A2 = W @ S.conj(), S @ W.conj()
        nz = np.argmax(abs(A2))
        u = A1.flat[nz] / A2.flat[nz]
        assert abs(abs(u) - 1) < 1e-9 and np.allclose(A1, u * A2)
        S = np.sqrt(u) * S                    # now phi commutes with the star; -S is the other choice (T and i T)
        assert np.allclose(S @ S.conj(), np.eye(6)) and np.allclose(W @ S.conj(), S @ W.conj())
        gens = []
        for c in range(6):
            for z in (1, 1j):
                b = np.zeros(6, dtype=complex)
                b[c] = z
                gens.append(b + S @ b.conj())
        G = np.array([np.concatenate([g.real, g.imag]) for g in gens]).T          # 12 x 12, rank 6
        U, sv, _ = np.linalg.svd(G)
        assert sv[5] > 1e-6 and sv[6] < 1e-9
        Bas = U[:, :6]
        img = []
        for c in range(6):
            t = Bas[:6, c] + 1j * Bas[6:, c]
            y = W @ t.conj()
            img.append(np.concatenate([y.real, y.imag]))
        R = np.linalg.lstsq(Bas, np.array(img).T, rcond=None)[0]
        co = np.poly(R).real                  # t^6 - e1 t^5 + e2 t^4 - e3 t^3 + ...
        e1, e2, e3 = -co[1], co[2], -co[3]
        A, B = (e2 - p * p) / p, e1 * e1 - 2 * (e2 - p * p) - 2 * p * p
        At, Bt = 2 * a.real + 2 * a2.real, 2 * p * p + 4 * a.real * a2.real
        err = max(abs(e3), abs(co[6] + p ** 6) / p ** 5, abs(A - At), abs(B - Bt) / p)
        worst = max(worst, err)
        assert err < 1e-6 * p ** 2, (p, a, a2, e1, e2, e3, A, At, B, Bt)
    print(f"selftest: {trials} random antilinear maps; in each the polynomial on T is (t^2 - p^2)(t^4 - e1 t^3 + ...) with "
          f"determinant -p^6, and A = g g'/p, B = g^2 + g'^2 - 2 p^2 are the coefficients of the polynomial of phi^2; "
          f"largest error {worst:.1e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    a = sub.add_parser("involution")
    a.add_argument("primes", type=int, nargs="+")
    a.add_argument("--file", default="sixlines_prym.jsonl")
    a.add_argument("--workers", type=int, default=4)
    a = sub.add_parser("commuting")
    a.add_argument("file")
    a = sub.add_parser("types")
    a.add_argument("primes", type=int, nargs="+")
    a.add_argument("--file", default="sixlines_prym.jsonl")
    a.add_argument("--workers", type=int, default=4)
    a = sub.add_parser("inverting")
    a.add_argument("primes", type=int, nargs="+")
    a.add_argument("--file", default="sixlines_prym.jsonl")
    a.add_argument("--pairs", type=int, default=2)
    a.add_argument("--limit", type=int, default=0)
    a.add_argument("--k", type=int, default=2, choices=(2, 3))
    a.add_argument("--workers", type=int, default=4)
    a.add_argument("--out")
    args = ap.parse_args()
    if args.cmd == "selftest":
        run_selftest()
    elif args.cmd == "involution":
        run_involution(args.primes, args.file, args.workers)
    elif args.cmd == "commuting":
        run_commuting(args.file)
    elif args.cmd == "types":
        run_types(args.primes, args.file, args.workers)
    else:
        run_inverting(args.primes, args.file, args.pairs, args.limit, args.k, args.workers, args.out)


if __name__ == "__main__":
    main()
