#!/usr/bin/env python
"""
sixlines_review.py -- the re-derivation checks of notes 13, 14, 15 and 16 (note 18).

It uses the functions of sixlines_prym.py and sixlines_dictionary.py unchanged and adds five checks.

  sign-minus FILE     Surfaces with sign - in an output file of sixlines_prym.py, pairs of commuting type.  The
                      script of note 14 records (c1, c2, c3) for them and does not compare.  Here: whether
                      c3 = -p conj(c1) and c2 is purely imaginary (so that the constant term of chi is -p^2), and
                      whether, with that constant term, the second exterior power of the quartic chi_2 whose roots
                      are the squares of those of chi is the polynomial of the surface over F_{p^2}, or its twist.

  inverting p n       Pairs of the type "Frobenius inverts J" with c4 rational (so p = 3 mod 4), which note 14 does
                      not compute.  The curve v^4 = F is defined over F_p, and Frobenius^2 commutes with J.  Traces
                      of Frobenius^2, ^4, ^6 on the i-eigenspace by pr_trace (mode A over F_{p^2}); the quartic
                      chi_2 must have rational coefficients, and its second exterior power is compared with the
                      polynomial of the surface over F_{p^2}.  At most n pairs per configuration.  p = 7: seconds
                      per pair; p = 11: about a minute per pair.

  other-half p ...    Note 14, section 5.3.  For every configuration X and ordered pair (L5, L6): X' is the
                      configuration with L6 replaced by its conjugate under the involution of the cross-ratio map.
                      The quartic computed with the OTHER half of tau is compared with the polynomial of X' (and
                      of X).  Prediction of note 14: it is that of X'.

  step [n]            Note 15, section 3.3.  For the first n polynomials of sixlines_cm.jsonl with chi conj(chi)
                      irreducible and up to 6 primes l = 1 mod 4 at which chi splits: aa = P1 P2 conj(P3) conj(P4),
                      and Psi as defined in note 15 section 3.1 (the ideal (l, P_a(t)) for each prime).  Prints
                      whether Psi(aa) is l * l12 / l34 or its square.

  kernel              Note 15, section 3.1.  The run of `sixlines_dictionary.py dictionary` with the structure, and
                      not only the order, of the kernel and the cokernel of Cl^-(K) -> G.

  python -I sixlines_review.py sign-minus sixlines_prym.jsonl
  python -I sixlines_review.py inverting 7 4
  python -I sixlines_review.py inverting 11 1
  python -I sixlines_review.py other-half 13 17 19 [--workers n]
  python -I sixlines_review.py step 8
  python -I sixlines_review.py kernel [--workers n]

Run from experiments/.  PARI's bnfinit is used without certification (GRH) in step and kernel.
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

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
import sixlines_dictionary as sd  # noqa: E402
import sixlines_prym as sp  # noqa: E402

pari = sp.pari
HERE = os.path.dirname(os.path.abspath(__file__))

# the polynomial whose roots are the squares of the roots of f
pari("rv_squares(f) = my(g = f * subst(f, 'X, -'X)); subst(substpol(g, 'X^2, 'Y), 'Y, 'X)")
# the characteristic polynomial of the second exterior power of the companion matrix of a quartic
pari("rv_wedge(f) = my(M = matcompanion(f), idx = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]], W); "
     "W = matrix(6, 6, r, s, M[idx[r][1], idx[s][1]] * M[idx[r][2], idx[s][2]] - M[idx[r][1], idx[s][2]] * M[idx[r][2], idx[s][1]]); "
     "charpoly(W, 'X)")


def surface_pol(p, e, sign):
    """P(X) = X^6 - e1 X^5 + e2 X^4 - e3 X^3 + sign (p^2 e2 X^2 - p^4 e1 X + p^6)."""
    e1, e2, e3 = e
    return pari(f"X^6-({e1})*X^5+({e2})*X^4-({e3})*X^3+({sign})*(({p * p * e2})*X^2-({p ** 4 * e1})*X+{p ** 6})")


def over_fp2(chi2, P):
    """'equal', 'twist' or 'no': the second exterior power of chi2 against the polynomial of the squares of the roots of P."""
    w, P2 = pari("rv_wedge")(chi2), pari("rv_squares")(P)
    if w == P2:
        return "equal"
    if w == pari.subst(P2, "X", pari("-X")):
        return "twist"
    return "no"


# ---------------------------------------------------------------- sign-minus

def run_sign_minus(path):
    st = Counter()
    for line in open(path):
        r = json.loads(line)
        if r["signs"] != [-1]:
            continue
        p = r["p"]
        P = surface_pol(p, r["e"], -1)
        had = False
        for x in r["pairs"]:
            if "c" not in x:
                continue
            had = True
            c1, c2, c3 = x["c"]
            st["pairs"] += 1
            st["c3 = -p conj(c1)"] += c3 == [-p * c1[0], p * c1[1]]
            st["c2 purely imaginary"] += c2[0] == 0
            st["c1 = 0"] += c1 == [0, 0]
            chi = pari(f"X^4-(({c1[0]})+I*({c1[1]}))*X^3+(({c2[0]})+I*({c2[1]}))*X^2-(({c3[0]})+I*({c3[1]}))*X-{p * p}")
            st["over F_p^2: " + over_fp2(pari("rv_squares")(chi), P)] += 1
        st["surfaces with such pairs"] += had
    print("sign -, pairs of commuting type:", dict(sorted(st.items())))


# ---------------------------------------------------------------- inverting

def run_inverting(p, max_pairs):
    assert p % 4 == 3, "c4 rational and Frobenius inverting J needs p = 3 mod 4"
    rows = [json.loads(line) for line in open(os.path.join(HERE, f"sixlines_p{p}.jsonl"))]
    st = Counter()
    for r in rows:
        cfg = r["cfg"]
        e = sp.surface_e(r["S"])
        signs = sp.sign_of(p, e)
        done = 0
        for pair in itertools.permutations(range(1, 7), 2):
            if done >= max_pairs:
                break
            rest = [i for i in range(1, 7) if i not in pair]
            L = pari.matrix(6, 3, [x for row in sp.forms(cfg) for x in row])
            N, D = sp.GP["pr_lambda"](L, rest + list(pair), p)
            cur = sp.GP["pr_curve"](N, D, p, 0)
            if int(cur[0]) == 0:
                continue
            W, bx, rat, Xc, yrat = cur[1], cur[2], int(cur[3]), cur[4], bool(int(cur[5]))
            if not yrat:                      # with p = 3 mod 4 this is the commuting type of note 14
                continue
            Qc = sp.GP["pr_findq"](W, bx, rat, Xc, p) if rat else pari("[0, 0, 0, 0]")
            T, t0 = [], time.time()
            for k in (2, 4, 6):
                cnt = sp.GP["pr_trace"](W, bx, rat, Qc, Xc, p, k, 0)
                if cnt == 0:
                    T = None
                    break
                n = [int(v) for v in cnt]
                T.append((-(n[0] - n[2]), -(n[1] - n[3])))
            if T is None:
                continue
            done += 1
            st["pairs"] += 1
            c = sp.quartic(T, p)
            out = f"p = {p} cfg = {cfg} e = {e} signs = {signs} pair = {pair}: T = {T}"
            if c is None or any(v[1] for v in c):
                st["chi_2 not rational"] += 1
                print(out + "  chi_2 NOT RATIONAL", flush=True)
                continue
            c1, c2, c3 = (v[0] for v in c)
            st["chi_2 rational"] += 1
            st["c3 = p^2 c1"] += c3 == p * p * c1
            chi2 = pari(f"X^4-({c1})*X^3+({c2})*X^2-({c3})*X+{p ** 4}")
            res = sorted({over_fp2(chi2, surface_pol(p, e, s)) for s in signs})
            st[f"sign {signs}: over F_p^2 {res}"] += 1
            print(out + f"  chi_2 = {chi2}  over F_p^2: {res}  ({time.time() - t0:.0f} s)", flush=True)
    print(f"inverting, p = {p}:", dict(sorted(st.items())), flush=True)


# ---------------------------------------------------------------- other-half

def half_job(job):
    p, cfg, pair, half = job
    sp.OTHER_HALF = half
    res = sp.one_pair(p, cfg, pair)
    return tuple(cfg), pair, half, res.get("status"), res.get("wedge")


def against(w, e):
    if any(v[1] for v in w):
        return "no"
    got = tuple(v[0] for v in w)
    return "equal" if got == tuple(e) else ("twist" if got == (-e[0], e[1], -e[2]) else "no")


def run_other_half(p, workers):
    sd.CONJ = pari(r"""(N, D, p) -> my(l5, l6, f5, f6);
        if(subst(D, 't, 0) == 0 || polcoef(D, 2) == 0, return(0));
        l5 = subst(N, 't, 0) / subst(D, 't, 0); l6 = polcoef(N, 2) / polcoef(D, 2);
        f5 = N - l5 * D; f6 = N - l6 * D;
        if(polcoef(f5, 2) == 0 || polcoef(f6, 1) == 0, return(0));
        [lift(-polcoef(f5, 1) / polcoef(f5, 2)), lift(-polcoef(f6, 0) / polcoef(f6, 1))]""")
    rows = [json.loads(line) for line in open(os.path.join(HERE, f"sixlines_p{p}.jsonl"))]
    bykey = {sd.encode(r["cfg"], p): r for r in rows}
    byassoc = {sd.encode(r["assoc"], p): r for r in rows}
    pts, meta = [], []
    for r in rows:
        forms = sp.forms(r["cfg"])
        for pair in itertools.permutations(range(1, 7), 2):
            cl = sd.conjugate_lines(p, forms, pair)
            if cl is None:
                continue
            pts.append([forms[i - 1] for i in range(1, 7) if i not in pair] + [forms[pair[0] - 1], cl[1]])
            meta.append((tuple(r["cfg"]), pair))
    keys = sd.canonical_keys(np.array(pts, dtype=np.int64), p)
    target = {}
    for k, m in zip(keys, meta):
        if k != -1:
            target[m] = bykey.get(int(k)) or byassoc.get(int(k))
    jobs = [(p, list(cfg), pair, h) for (cfg, pair) in target for h in (0, 1)]
    with Pool(workers) as pool:
        out = pool.map(half_job, jobs, chunksize=8)
    rowof = {tuple(r["cfg"]): r for r in rows}
    st = Counter()
    for cfg, pair, half, status, w in out:
        X, Y = rowof[cfg], target[(cfg, pair)]
        if Y is None:
            st[f"half {half}: X' not in the list"] += 1
            continue
        eX, eY = sp.surface_e(X["S"]), sp.surface_e(Y["S"])
        sX, sY = sp.sign_of(p, eX), sp.sign_of(p, eY)
        name = "the half of note 14" if half == 0 else "the other half"
        key = f"{name}; sign of X {sX}, of X' {sY}: "
        if status != "computed":
            st[key + status] += 1
        else:
            st[key + f"computed, against X: {against(w, eX)}, against X': {against(w, eY)}"] += 1
            if half == 1 and sY == [1]:
                st["the other half, X' of sign +: pairs computed"] += 1
                st["the other half, X' of sign +: P or P(-t) of X'"] += against(w, eY) != "no"
            if half == 0 and sX == [1]:
                st["the half of note 14, X of sign +: pairs computed"] += 1
                st["the half of note 14, X of sign +: P or P(-t) of X"] += against(w, eX) != "no"
    print(f"other-half, p = {p}:")
    for k in sorted(st):
        print(f"  {k}: {st[k]}")
    sys.stdout.flush()


# ---------------------------------------------------------------- step

pari(r"""
rv_step(P, q, e1, mu, nu, B) = {
  \\ rows [l, norm of aa, norm of Psi(aa), Psi(aa) == l l12 / l34, Psi(aa) == (l l12 / l34)^2, aa conj(aa) == (l)]
  my(a = 4 * (nu + mu * I), chi = 'X^4 - a * 'X^3 + e1 * 'X^2 - q * conj(a) * 'X + q^2, W8 = real(chi * conj(chi)),
     t = variable(P), bnfK, bnfE, io = 0, piK, res = List(), il, chil, cb, rts, rb, PP, PPc, psi, lj, aa, img, step);
  if(!polisirreducible(W8), return(0));
  bnfK = bnfinit(W8, 1); bnfE = bnfinit(P, 1); piK = Mod('X, W8);
  foreach(nfroots(bnfK, 'x^2 + 1), z,
    if(piK^4 - 4 * (nu + mu * z) * piK^3 + e1 * piK^2 - 4 * q * (nu - mu * z) * piK + q^2 == 0, io = z));
  forprime(l = 5, B,
    if(l % 4 != 1 || (bnfE.nf.disc * bnfK.nf.disc * poldisc(P) * poldisc(W8) * q) % l == 0, next);
    il = lift(sqrt(Mod(-1, l)));
    chil = Mod(1, l) * 'x^4 - 4 * (nu + mu * il) * 'x^3 + e1 * 'x^2 - 4 * q * (nu - mu * il) * 'x + q^2;
    cb = Mod(1, l) * 'x^4 - 4 * (nu - mu * il) * 'x^3 + e1 * 'x^2 - 4 * q * (nu + mu * il) * 'x + q^2;
    rts = polrootsmod(chil, l); if(#rts != 4, next);
    rb = vector(4, j, q / rts[j]);
    \\ Psi of note 15: the prime (l, i -+ i_l, pi - a0) goes to (l, P_a0(t)), P_a0 = a0^3 (cl / (x - a0))(t / a0)
    psi = ((cl, a0) -> idealadd(bnfE, l, Mod(subst(lift(a0^3 * subst(cl / ('x - a0), 'x, 'x / a0)), 'x, t), P)));
    lj = ((u, v) -> idealadd(bnfE, l, Mod(t - lift(rts[u] * rts[v]), P)));
    PP = vector(4, j, idealadd(bnfK, idealadd(bnfK, l, io - il), piK - lift(rts[j])));
    PPc = vector(4, j, idealadd(bnfK, idealadd(bnfK, l, io + il), piK - lift(rb[j])));
    aa = idealmul(bnfK, idealmul(bnfK, PP[1], PP[2]), idealmul(bnfK, PPc[3], PPc[4]));
    img = idealmul(bnfE, idealmul(bnfE, psi(chil, rts[1]), psi(chil, rts[2])), idealmul(bnfE, psi(cb, rb[3]), psi(cb, rb[4])));
    step = idealmul(bnfE, l, idealdiv(bnfE, lj(1, 2), lj(3, 4)));
    listput(res, [l, idealnorm(bnfK, aa), idealnorm(bnfE, img), img == step, img == idealpow(bnfE, step, 2),
                  idealmul(bnfK, aa, idealmul(bnfK, idealmul(bnfK, PPc[1], PPc[2]), idealmul(bnfK, PP[3], PP[4]))) == idealhnf(bnfK, l)]);
    if(#res >= 6, break));
  Vec(res);
}
""")


def run_step(count):
    st = Counter()
    for line in open(os.path.join(HERE, "sixlines_cm.jsonl")):
        if st["polynomials"] >= count:
            break
        r = json.loads(line)
        q, e = r["p"], r["e"]
        mu, nu = sd.mu_nu(q, e)
        out = pari("rv_step")(cm.pol(cm.charpoly(q, e)), q, e[0], mu, nu, 5000)
        if out == 0:
            continue
        st["polynomials"] += 1
        for row in out:
            l, na, ni, first, second, selfdual = (int(v) for v in row)
            assert na == l ** 4 and selfdual == 1
            st["primes"] += 1
            st["norm of Psi(aa) is l^12"] += ni == l ** 12
            st["Psi(aa) = l l12 / l34"] += first
            st["Psi(aa) = (l l12 / l34)^2"] += second
        print(f"  q = {q} e = {e}: l = {[int(row[0]) for row in out]}", flush=True)
    print("step:", dict(sorted(st.items())))


# ---------------------------------------------------------------- kernel

def kernel_setup():
    """dc_one of sixlines_dictionary.py with two more outputs: elementary divisors of kernel and cokernel."""
    src = open(os.path.join(HERE, "sixlines_dictionary.py")).read()
    a = src.index("dc_one(P, q, e1, mu, nu, B) = {")
    b = src.index('""")', a)
    body = src[a:b].replace("dc_one(", "rv_two(", 1)
    old = "  [1, n, lmax, cK, cE, ok_hom, ok_rel, ok_sub, dom, tgt, coker,"
    assert body.count(old) == 1
    new = """  my(kerS = [], cokS = [], KI, A, Lam);
  if(kK && kE,
    cokS = select(d -> d > 1, matsnf(mathnf(concat(HE, F))));
    KI = matkerint(concat(F, HE)); A = matrix(kK, matsize(KI)[2], r, c, KI[r, c]);
    Lam = mathnf(concat(A, LK));
    kerS = select(d -> d > 1, matsnf(Lam^(-1) * LK)));
  if(kK && !kE, kerS = select(d -> d > 1, matsnf(LK)));
  if(!kK && kE, cokS = select(d -> d > 1, matsnf(HE)));
  [1, n, lmax, cK, cE, ok_hom, ok_rel, ok_sub, dom, tgt, coker, kerS, cokS,"""
    pari(body.replace(old, new))


def kernel_one(job):
    q, e, bound = job
    mu, nu = sd.mu_nu(q, e)
    res = pari("rv_two")(cm.pol(cm.charpoly(q, e)), q, e[0], mu, nu, bound)
    if res == 0 or int(res[0]) == 0:
        return None
    ker = sorted(int(x) for x in res[11])
    cok = sorted(int(x) for x in res[12])
    assert math.prod(ker) * int(res[9]) == int(res[8]) * int(res[10]) and math.prod(cok) == int(res[10])
    assert int(res[5]) and int(res[6]) and int(res[7])
    return ker, cok


def run_kernel(workers):
    jobs, seen = [], set()
    for line in open(os.path.join(HERE, "sixlines_cm.jsonl")):
        r = json.loads(line)
        jobs.append((r["p"], r["e"], 2000))
        seen.add((r["p"], tuple(r["e"])))
    for line in open(os.path.join(HERE, "sixlines_groupsize.jsonl")):
        r = json.loads(line)
        if r.get("source") == "all" and r.get("maximal") and (r["q"], tuple(r["e"])) not in seen:
            jobs.append((r["q"], r["e"], 2000))
            seen.add((r["q"], tuple(r["e"])))
    with Pool(workers, initializer=kernel_setup) as pool:
        recs = [r for r in pool.map(kernel_one, jobs, chunksize=2) if r]
    K, C = Counter(tuple(r[0]) for r in recs), Counter(tuple(r[1]) for r in recs)
    print(f"kernel: {len(recs)} polynomials")
    print("  elementary divisors of the kernel of Cl^-(K) -> G:", dict(sorted(K.items())))
    print("  elementary divisors of the cokernel:", dict(sorted(C.items())))
    print("  kernel killed by 2:", sum(n for k, n in K.items() if all(d == 2 for d in k)),
          " cokernel killed by 2:", sum(n for k, n in C.items() if all(d == 2 for d in k)),
          " cokernel killed by 4:", sum(n for k, n in C.items() if all(4 % d == 0 for d in k)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("sign-minus")
    a.add_argument("file")
    a = sub.add_parser("inverting")
    a.add_argument("p", type=int)
    a.add_argument("pairs", type=int)
    a = sub.add_parser("other-half")
    a.add_argument("primes", type=int, nargs="+")
    a.add_argument("--workers", type=int, default=8)
    a = sub.add_parser("step")
    a.add_argument("count", type=int, nargs="?", default=8)
    a = sub.add_parser("kernel")
    a.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    if args.cmd == "sign-minus":
        run_sign_minus(args.file)
    elif args.cmd == "inverting":
        run_inverting(args.p, args.pairs)
    elif args.cmd == "other-half":
        for p in args.primes:
            run_other_half(p, args.workers)
    elif args.cmd == "step":
        run_step(args.count)
    else:
        run_kernel(args.workers)


if __name__ == "__main__":
    main()
