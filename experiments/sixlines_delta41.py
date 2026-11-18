#!/usr/bin/env python
"""
sixlines_delta41.py -- the sign delta of note 23 against the character chi_{-1} at q = 41 = 1 mod 8 (note 26).

Notes 23 and 24 compute delta on three torsors at q = 29 and 37, both 5 mod 8, and find that it is the character
chi_{-1} of G.  Note 25 section 7 names the nearest test out of sample: q = 41, where note 19 lists two polynomials
with the character, e = [6, -369, -20172] (G = Z/16) and e = [22, -369, -73964] (G = Z/12).  This script does that
test.  It reads no file of the earlier notes except the scripts, and writes only files of its own:

    sixlines_p41.jsonl          written by `sixlines.py enumerate --p 41`      (input here)
    sixlines_cm_p41.jsonl       written by `sixlines_cm.py --list ... sixlines_p41.jsonl`   (input here)
    sixlines_prym_p41.jsonl     `prym`:   the curves of note 14 for the configurations of the two polynomials
    sixlines_delta41.jsonl      `delta`:  one record for each curve and reading
    sixlines_delta41_cosets.jsonl   `cosets`: the two cosets of ker chi_{-1} on the surfaces, three routes
    sixlines_delta41_compare.json   `compare`

THE DEFINITION AT q = 1 mod 8.  Written on 2026-11-15 before any curve at p = 41 was computed (the docstring was
in the file before the first run of `delta`; nothing below this line was changed after a comparison with the cosets
-- if it had been, note 26 would say so first; the one edit made later is a cross-reference, here and
  in the first line).

  Notation of note 23 section 5: T the 2-adic Tate module of the Prym variety, V_1 the line over Q_2(i) on which
  Frobenius pi acts by its root pi_1 in Z_2[i], h(m) = E(m, i m) on V_1, u = the odd part of h(m) modulo 4 (the
  class of h modulo norms from Q_2(i) and powers of 2), delta = (-1)^((u-1)/2) sign(Re c_1 Im c_1) sign(c_2), with
  (c_1, c_2) the coefficients of the quartic chi of the curve C_0 : v^4 = F as sixlines_prym.py gives them.

  What does not change.  u is defined on the Tate module over the algebraic closure and does not see q modulo 8.
  The normalisation does not either: it uses only that u changes sign with the conjugate structure (c_1 -> conj
  c_1) and is the same for the four twists (c_1 -> i^tw c_1, c_2 -> (-1)^tw c_2), and for both polynomials at
  q = 41 (mu, nu) = (1, 1), so c_1 = 4 (+-1 +- i) as at q = 29 and 37.  So delta is DEFINED BY THE SAME FORMULA.

  What changes is how u is read.
    * N(pi_1) = p.  For p = 5 mod 8, pi_1 = (unit) (1 + 2 i (odd)) and pi_1^4 - 1 has valuation s = 7 at 1 + i for
      every twist.  For p = 1 mod 8, pi_1 = i^k (1 + 4 beta), and for p = 41 and both polynomials beta is a unit:
      v(pi_1^r - 1) for r = 1, 2, 4 is (4, 6, 8) on the twist with pi_1 = 1 mod 4, (2, 6, 8) on the twist with
      pi_1 = -1 mod 4, and (1, 2, 8) on the other two (computed from chi alone; chi'(pi_1) has valuation 6 as at
      q = 29 and 37, so pi - pi_1 is still 2 times a unit on the complement of V_1).
    * roots of unity: v_2(p^r - 1) = 3, 4, 5 for r = 1, 2, 4 (against 2, 3, 4 for p = 5 mod 8).  So over F_{p^4}
      the Tate pairing has order N = 32; the V_1-part of the rational 2-power torsion is Z[i]/(1 + i)^8, of
      exponent 16.
    * the identity of note 23 section 5.3 is unchanged: for X in the kernel of pi - pi_1 on the part with
      J^2 = -1, killed by N,   Z := log T_N(X, iX) + i log T_N(X, X) = M conj(N / varpi_1)  mod N,
      varpi_1 = pi_1^r - 1, M = N(gamma) h(m).  It determines M modulo 2^k, k = min(n, ceil(s / 2)), N = 2^n,
      s = v(varpi_1).  At q = 29, 37: s = 7, n = 4, k = 4, M is one number modulo 16.  At q = 41, r = 4: s = 8,
      n = 5, k = 4: M is a class modulo 16, two numbers modulo 32.
    * READING RULE.  A kernel vector gives u when its solutions M modulo N form exactly one class modulo 2^k, not
      the class of 0, and v_2(M) + 2 <= k; then u = the odd part of M modulo 4.  (For s = 7 this is the rule of
      sixlines_jacobian.py.)  Among the kernel vectors that give u the one with the least v_2(M) is taken (the
      first such, as in note 23), its point X is formed, and X is verified directly: its order, pi X = pi_1 X,
      J^2 X = -X, and the two pairings of X computed directly equal the values assembled from the matrices.
      Status "ok" needs all four.
    * the complement of V_1 contributes (N / 4) E(t, varpi t), which is 0 modulo N when varpi = pi^r - 1 is
      divisible by 4 on the complement: for r = 4 (divisible by 8) and for r = 2 on the twists with pi = +-1
      modulo 2 there.  For r = 1 it is not, and there M is known modulo 4 only: THE PAIRING OVER F_p IS NOT USED.

  PRIMARY READING:  C_0, r = 4, N = 32, seed 1, twelve classes -- the computation of note 23 with N = 2^v_2(p^4-1)
  as its code already has it, and the reading rule above.  delta of a surface = the common value over its
  configurations and pairs with status ok.
  CONTROLS (note 23 section 5.5), each compared with the primary u of the same curve:
    (a) other pairs of lines of the configuration and of its associate (this is the primary run itself);
    (b) the twists C_1, C_2 over F_{p^4}, N = 32;
    (c) another seed (other random classes);
    (d) N = 16 over F_{p^4} (s = 8 = 2n: N / varpi_1 is a unit, M modulo 16);
    (e) F_{p^2} with N = 16 on the twists where s = 6 (M modulo 8, so only vectors with M = 2 (odd)).
  If the exponent of (1 - J) G is above N the curve is "not computed", as in note 23.  One extension is made and
  reported separately, never mixed into the primary count: with 2^c the least power for which 2^c N kills the
  twelve classes z_a, the classes 2^c z_a are used in their place (reading "cofactor").

Commands (run from experiments/ with python -I):

  list                         the configurations and surfaces of the two polynomials
  prym [--pairs n]             sixlines_prym.config_record on those configurations (n = 0: all 30 ordered pairs)
  dyadic                       the table of v(pi_1^r - 1) by twist, for the two polynomials and those of note 23
  validate [--n k]             the checks of note 23 section 3 (sixlines_jacobian.validate) on a few curves at p = 41
  selfcheck [--n k]            the routine read_u of this script on stored curves of note 23 at p = 29 and 37 (read
                               from sixlines_jacobian_character.jsonl): same M, u and delta as stored?
  delta --reading R            R in primary, twist1, twist2, seed, n16, r2, cofactor; appends to sixlines_delta41.jsonl
  cosets q e1 mu nu            the CM method at q: class polynomial, its roots against the enumerated surfaces, the
                               orbits of the ideals above primes l = 1 mod 4 (sixlines_delta_labels.run_one), the
                               factors over E(i) and the form modulo 4 of xi (sixlines_review2.cmd_galois).  With
                               --no-delta the imported code is given no delta, and its lines "equals the split by
                               delta" and "delta is constant on each orbit" then say False: the comparison is `compare`
  compare                      delta against the cosets
  signs                        the sign between delta and the form on the lattice, for the two polynomials at q = 41
                               and the three torsors of note 24: is delta = zeta * (value of the form indicator at
                               the root), zeta = n0^((q-1)/4) the fourth root of unity of sixlines_jacobian.py?  Reads
                               sixlines_delta41_signs.jsonl, written by `cosets 37 18 1 1 --digits 1200`,
                               `cosets 29 2 1 1 --digits 1200` and `cosets 37 2 1 1 --digits 2000` with
                               `--primes 3 --frob-primes 2 --out sixlines_delta41_signs.jsonl` (delta of note 23)

Uses cypari2 through the scripts of notes 14, 20, 23, 24 and 25, which are imported and not changed.  Their fixed
input names (sixlines_cm.jsonl, sixlines_prym.jsonl, sixlines_delta_labels.jsonl, the files of note 23) are
redirected here to the files above, in memory.  Class groups are PARI's without certification (GRH).
"""
import argparse
import collections
import json
import math
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

P = 41
CASES = [[6, -369, -20172], [22, -369, -73964]]
F_ENUM = "sixlines_p41.jsonl"
F_CM = "sixlines_cm_p41.jsonl"
F_PRYM = "sixlines_prym_p41.jsonl"
F_DELTA = "sixlines_delta41.jsonl"
F_COSETS = "sixlines_delta41_cosets.jsonl"
F_COMPARE = "sixlines_delta41_compare.json"


def jl(name):
    with open(os.path.join(HERE, name)) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def surface_e(S):
    s1, s2, s3 = S
    return [s1, (s1 * s1 - s2) // 2, (s1 ** 3 - 3 * s1 * s2 + 2 * s3) // 6]


def surfaces_of(e):
    """The configurations with polynomial P(t) or P(-t), grouped into surfaces by association (as
    sixlines_jacobian.surfaces_of, from sixlines_p41.jsonl)."""
    rows = {tuple(r["cfg"]): r for r in jl(F_ENUM)}
    cfgs = [cfg for cfg, r in rows.items()
            if surface_e(r["S"]) in (list(e), [-e[0], e[1], -e[2]])]
    surf, seen = [], set()
    for cfg in sorted(cfgs):
        if cfg in seen:
            continue
        a = tuple(rows[cfg]["assoc"])
        assert a in rows and surface_e(rows[a]["S"]) in (list(e), [-e[0], e[1], -e[2]])
        seen.update([cfg, a])
        surf.append(sorted({cfg, a}))
    return surf, rows


def cmd_list(a):
    cmrows = {tuple(r["e"]): r for r in jl(F_CM)}
    for e in CASES:
        surf, rows = surfaces_of(e)
        r = cmrows.get(tuple(e))
        print(f"q = {P}, e = {e}: {sum(len(s) for s in surf)} configurations, {len(surf)} surfaces "
              f"({sum(1 for s in surf if len(s) == 1)} self-associated)")
        if r:
            print(f"  sixlines_cm: a_plus {r['a_plus']}, a_minus {r['a_minus']}, self-associated {r['n_self_assoc']}, "
                  f"k3 {r['k3']:g}, N {r['N']}, N_total {r['N_total']}, index {r['index']}, h {r['h']}, h0 {r['h0']}, "
                  f"verdict {r['verdict']}")
        for i, s in enumerate(surf):
            print(f"  surface {i}: " + ", ".join(f"{c} S1 = {rows[c]['S'][0]}" for c in s))


# ------------------------------------------------------------------------------------------------ prym

def cmd_prym(a):
    import sixlines_prym as sp
    jobs = []
    for e in CASES:
        surf, rows = surfaces_of(e)
        for s in surf:
            for cfg in s:
                jobs.append((P, rows[cfg], a.pairs))
    t0 = time.time()
    with Pool(a.workers) as pool:
        recs = pool.map(sp.config_record, jobs, chunksize=1)
    with open(os.path.join(HERE, F_PRYM), "w") as fh:
        for rec in recs:
            fh.write(json.dumps(rec) + "\n")
    tally = collections.Counter()
    for rec in recs:
        for x in rec["pairs"]:
            tally[(x["compare"], x.get("rational_point"))] += 1
        for k, v in rec["seen"].items():
            tally[("seen", k)] += v
    print(f"{len(recs)} configurations, {time.time() - t0:.0f} s; pairs by (compare, rational_point) and statuses seen: "
          f"{dict(sorted(tally.items(), key=str))}")


def stored_pairs():
    out = collections.defaultdict(list)
    for r in jl(F_PRYM):
        for x in r["pairs"]:
            if x.get("compare") in ("equal", "twist"):
                out[tuple(r["cfg"])].append((tuple(x["pair"]), [list(v) for v in x["c"]], x["compare"], x["rational_point"]))
    return out


# ------------------------------------------------------------------------------------------------ dyadic

def cmd_dyadic(a):
    import sixlines_jacobian as sj
    print("valuation at 1 + i of pi_1^r - 1, r = 1, 2, 4, 8, and of chi'(pi_1), by twist; chi with c1 = 4(1 + i), c2 = e1")
    for q, e1 in [(29, 2), (37, 2), (37, 18), (41, 6), (41, 22)]:
        c = [[4, 4], [e1, 0], [4 * q, -4 * q]]
        for tw in range(4):
            co = sj.chi_twist(c, q, tw)
            rts, degs = sj.dyadic_root(co, 40)
            if len(rts) != 1:
                print(f"  q = {q}, e1 = {e1}, twist {tw}: {len(rts)} roots in Z_2[i], degrees {degs}")
                continue
            x = rts[0]
            mod = 2 ** 40
            dco = [(k * cf[0], k * cf[1]) for k, cf in enumerate(co)][1:]
            vd = sj.gval(sj.geval(dco, x, mod))
            ss = [sj.gval(sj.gadd(sj.gpow(x, r, mod), (-1, 0))) for r in (1, 2, 4, 8)]
            print(f"  q = {q}, e1 = {e1}, twist {tw}: pi_1 = {x[0] % 16} + {x[1] % 16} i mod 16, v(pi_1^r - 1) = {ss}, "
                  f"v(chi'(pi_1)) = {vd}, v_2(q^r - 1) = {[sj.val(q ** r - 1, 2) for r in (1, 2, 4, 8)]}")


# ------------------------------------------------------------------------------------------------ delta

def read_u(job):
    """One curve and one reading.  The computation of sixlines_jacobian.character_pair, with N a parameter, the
    reading rule of the docstring, and an optional cofactor."""
    import sixlines_jacobian as sj
    pari = sj.pari
    p, cfg, pair, c, tw, r, nn, mm, seed, cofactor, reading = job
    out = {"p": p, "cfg": list(cfg), "pair": list(pair), "c": c, "tw": tw, "r": r, "seed": seed, "reading": reading}
    t0 = time.time()
    co = sj.chi_twist(c, p, tw)
    B = 40
    rts, degs = sj.dyadic_root(co, B)
    out["dyadic_degrees"] = degs
    if len(rts) != 1:
        out["status"] = "no single root in Z_2[i]"
        return out
    x = rts[0]
    mod = 2 ** B
    varpi = sj.gadd(sj.gpow(x, r, mod), (-1, 0))
    s = sj.gval(varpi)
    assert (p ** r - 1) % 2 ** nn == 0
    N = 2 ** nn
    k = min(nn, (s + 1) // 2)
    out.update({"s": s, "N": N, "k": k})
    if s > 2 * nn:
        out["status"] = "N / (pi_1^r - 1) is not integral"
        return out
    try:
        A = sj.Jac(p, list(cfg), tuple(pair), tw, r, seed=seed)
    except ValueError as err:
        out["status"] = str(err)
        return out
    n = A.card(c)
    m = n >> sj.val(n, 2)
    ys = [A.mul(m, A.random()) for _ in range(mm)]
    zs = [A.sub(y, A.J(y, 1)) for y in ys]
    cof = 0
    while not all(A.is_zero(A.mul(N << cof, z)) for z in zs):
        cof += 1
        assert cof < 12
    out["exponent_over_N"] = cof                       # the twelve classes z_a are killed by 2^cof N and not less
    if cof and not cofactor:
        out["status"] = "exponent of (1 - J) G larger than N"
        out["seconds"] = round(time.time() - t0)
        return out
    if cof:
        zs = [A.mul(2 ** cof, z) for z in zs]
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

    ups = sj.gdiv1pi(varpi, s)
    w = sj.gmul(sj.gdiv1pi((N, 0), s), sj.ginv(sj.gmod(ups, mod), mod))      # N / varpi
    wb = (w[0] % N, (-w[1]) % N)                                             # its conjugate
    wo = (w[0] % N, w[1] % N)                                                # the other convention (a control)

    def solve(Z, ww):
        return sorted({M for M in range(N) if all((M * ww[j] - Z[j]) % N == 0 for j in (0, 1))})

    def reading_of(sol):
        """(M modulo 2^k, v_2(M), u) when the solutions are one class modulo 2^k, not 0, with v_2 + 2 <= k."""
        cls = {M % 2 ** k for M in sol}
        if len(cls) != 1 or len(sol) != N >> k:
            return None
        M = cls.pop()
        if M == 0:
            return None
        v2 = sj.val(M, 2)
        if v2 + 2 > k:
            return None
        return M, v2, (M >> v2) % 4

    vecs = [v for v in vecs if order(v) > 1]
    tally = collections.Counter()
    other = 0
    best = None
    us = collections.Counter()
    for cv in vecs:
        Z = zof(cv)
        sol = solve(Z, wb)
        rd = reading_of(sol)
        other += bool(solve(Z, wo)) and any(Z)
        tally[(order(cv), tuple(M % 2 ** k for M in sol[:1]), len(sol), None if rd is None else rd[2])] += 1
        if rd is not None:
            us[rd[2]] += 1
            if best is None or rd[1] < best[1][1]:
                best = (cv, rd, Z)
    out.update({"kernel_vectors": len(vecs),
                "kernel": sorted(([o, list(M), ns, u_, cnt] for (o, M, ns, u_), cnt in tally.items()), key=str),
                "unsolvable": sum(cnt for (o, M, ns, u_), cnt in tally.items() if ns == 0),
                "solvable_with_other_convention": other,
                "u_over_kernel": {str(key): v for key, v in sorted(us.items())}})
    if best is None:
        out["status"] = "no kernel vector that gives u"
        out["seconds"] = round(time.time() - t0)
        return out
    cv, (M, v2, u), Z = best
    X = A.zero
    for a in range(mm):
        if cv[a]:
            X = A.add(X, A.mul(cv[a], zs[a]))
    o = order(cv)
    d = A.sub(A.frob(X), A.add(A.mul(al, X), A.mul(be, A.J(X, 3))))
    Zd = (A.tate_matrix([X], [A.J(X, 3)], N)[0][0], A.tate_matrix([X], [X], N)[0][0])
    out["verified"] = {"order": A.is_zero(A.mul(o, X)) and not A.is_zero(A.mul(o // 2, X)), "kernel": A.is_zero(d),
                       "J^2 = -1": A.is_zero(A.add(X, A.J(X, 2))), "Z": list(Zd) == list(Z)}
    c1, c2 = c[0], c[1]
    sg = (1 if c1[0] * c1[1] > 0 else -1) * (1 if c2[0] > 0 else -1)
    out.update({"M": M, "v2_M": v2, "order_X": o, "Z": list(Z), "u": u, "sign": sg, "delta": (1 if u == 1 else -1) * sg,
                "u_constant_over_kernel": len(us) == 1,
                "status": "ok" if all(out["verified"].values()) else "not verified",
                "seconds": round(time.time() - t0)})
    return out


def _guarded(job):
    try:
        return read_u(job)
    except Exception as err:                     # reported as a status and counted
        return {"p": job[0], "cfg": list(job[1]), "pair": list(job[2]), "c": job[3], "tw": job[4], "r": job[5],
                "seed": job[8], "reading": job[10], "status": "error: " + str(err)[:300]}


READINGS = {            # name: (twists, r, n, seed, cofactor)
    "primary": ([0], 4, 5, 1, False),
    "twist1": ([1], 4, 5, 1, False),
    "twist2": ([2], 4, 5, 1, False),
    "seed": ([0], 4, 5, 7, False),
    "n16": ([0], 4, 4, 1, False),
    "r2": ([0, 1], 2, 4, 1, False),          # of these two twists, the one with v(pi_1^2 - 1) = 6
    "cofactor": ([0], 4, 5, 1, True),
}


def s_of(c, tw, r):
    import sixlines_jacobian as sj
    rts, degs = sj.dyadic_root(sj.chi_twist(c, P, tw), 40)
    return sj.gval(sj.gadd(sj.gpow(rts[0], r, 2 ** 40), (-1, 0))) if len(rts) == 1 else None


def cmd_delta(a):
    tws, r, nn, seed, cofactor = READINGS[a.reading]
    done = jl(F_DELTA) if os.path.exists(os.path.join(HERE, F_DELTA)) else []
    prim_ok = {(tuple(x["cfg"]), tuple(x["pair"])) for x in done if x["reading"] == "primary" and x["status"] == "ok"}
    prim_exp = {(tuple(x["cfg"]), tuple(x["pair"])) for x in done if x["reading"] == "primary"
                and x["status"].startswith("exponent")}
    have = {(x["reading"], tuple(x["cfg"]), tuple(x["pair"]), x["tw"]) for x in done}
    sp_ = stored_pairs()
    jobs, owner = [], []
    for e in CASES:
        if a.e1 and e[0] not in a.e1:
            continue
        surf, rows = surfaces_of(e)
        for si, cfgs in enumerate(surf):
            if a.surfaces and f"{e[0]}:{si}" not in a.surfaces:
                continue
            per_surface = 0
            for cfg in cfgs:
                pairs = sp_[tuple(cfg)]
                if a.reading == "primary":
                    pairs = pairs[:a.pairs] if a.pairs else pairs
                elif a.reading == "cofactor":
                    pairs = [x for x in pairs if (tuple(cfg), x[0]) in prim_exp]
                    pairs = pairs[:a.pairs] if a.pairs else pairs
                else:
                    pairs = [x for x in pairs if (tuple(cfg), x[0]) in prim_ok]
                for pair, c, comp, rat in pairs:
                    if a.reading not in ("primary", "cofactor") and a.per_surface and per_surface >= a.per_surface:
                        break
                    used = False
                    for tw in tws:
                        if a.reading == "r2" and s_of(c, tw, r) != 6:
                            continue
                        if (a.reading, tuple(cfg), tuple(pair), tw) in have:
                            continue
                        jobs.append((P, list(cfg), pair, c, tw, r, nn, 12, seed, cofactor, a.reading))
                        owner.append((e, si))
                        used = True
                    per_surface += used
    print(f"reading {a.reading}: {len(jobs)} runs", flush=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        recs = pool.map(_guarded, jobs, chunksize=1)
    with open(os.path.join(HERE, F_DELTA), "a") as fh:
        for (e, si), rec in zip(owner, recs):
            rec["e"] = e
            rec["surface"] = si
            fh.write(json.dumps(rec) + "\n")
    st = collections.Counter((rec["e"][0], rec["status"][:60]) for rec in recs)
    print(f"  {time.time() - t0:.0f} s; by (e1, status): {dict(sorted(st.items()))}", flush=True)


def cmd_validate(a):
    import sixlines_jacobian as sj
    sp_ = stored_pairs()
    k = 0
    for e in CASES:
        surf, rows = surfaces_of(e)
        for cfgs in surf[:a.n]:
            cfg = cfgs[0]
            for pair, c, comp, rat in sp_[tuple(cfg)][:1]:
                for tw, r in ((0, 1), (1, 1), (0, 2)):
                    try:
                        out = sj.validate(P, list(cfg), tuple(pair), tw, r, npts=2, verbose=False)
                    except Exception as err:
                        print(f"  {cfg} {pair} tw {tw} r {r}: refused: {str(err)[:80]}")
                        continue
                    print(f"  {cfg} {pair} tw {tw} r {r}: {json.dumps(out['checks']) if out else None}", flush=True)
                    k += 1


def cmd_selfcheck(a):
    recs = [r for r in jl("sixlines_jacobian_character.jsonl") if r["status"] == "ok"]
    jobs, ref = [], []
    for p in (29, 37):
        for r in [x for x in recs if x["p"] == p][:a.n]:
            jobs.append((p, r["cfg"], tuple(r["pair"]), r["c"], r["tw"], 4, 4, 12, 1, False, "selfcheck"))
            ref.append(r)
    with Pool(a.workers) as pool:
        out = pool.map(_guarded, jobs, chunksize=1)
    same = 0
    for r, o in zip(ref, out):
        ok = o["status"] == "ok" and (o["M"], o["u"], o["delta"], o["order_X"]) == (r["M"], r["u"], r["delta"], r["order_X"])
        same += ok
        print(f"  p = {r['p']} {r['cfg']} {r['pair']}: stored M {r['M']} u {r['u']} delta {r['delta']} order {r['order_X']}; "
              f"here status {o['status']}, M {o.get('M')} u {o.get('u')} delta {o.get('delta')} order {o.get('order_X')}, "
              f"s {o.get('s')} N {o.get('N')} k {o.get('k')}, u over the kernel {o.get('u_over_kernel')}: "
              f"{'same' if ok else 'DIFFERENT'}", flush=True)
    print(f"read_u reproduces the stored record of note 23: {same} of {len(ref)}")


# ------------------------------------------------------------------------------------------------ cosets

def cmd_cosets(a):
    """The CM method at q with the unchanged code of notes 20, 24 and 25; their fixed input files are replaced in
    memory by the p = 41 files, and the lattices and period points are computed once."""
    import sixlines_delta_labels as dl
    import sixlines_review2 as rv
    cm = dl.cm
    dl.WORKERS = a.workers
    q, e1 = a.q, a.e1

    load_surfaces0, load_delta0 = cm.load_surfaces, dl.load_delta
    f_out = a.out or os.path.join(HERE, F_COSETS)

    def load_surfaces(q_, e1_, e2_):
        if q_ != P:                          # another prime (a test of this command): the file of note 2
            return load_surfaces0(q_, e1_, e2_)
        return [r for r in jl(F_CM) if r["p"] == q_ and abs(r["e"][0]) == abs(e1_) and r["e"][1] == e2_]

    def load_delta(q_):
        out = {}
        path = os.path.join(HERE, F_DELTA)
        if q_ != P:
            return load_delta0(q_)
        if a.no_delta or not os.path.exists(path):
            return out
        for r in jl(F_DELTA):
            if r["p"] == q_ and r["reading"] == "primary" and r.get("status") == "ok":
                out.setdefault(tuple(r["cfg"]), []).append(r["delta"])
        return out

    cache = {}
    Field0, Torsor0 = cm.CMField, cm.Torsor

    def field(*args):
        if ("F",) + args not in cache:
            cache[("F",) + args] = Field0(*args)
        return cache[("F",) + args]

    def torsor(fld, digits, verbose=True):
        key = ("T", id(fld), digits)
        if key not in cache:
            cache[key] = Torsor0(fld, digits, verbose)
        return cache[key]

    cm.load_surfaces, dl.load_delta, cm.CMField, cm.Torsor = load_surfaces, load_delta, field, torsor
    rec = None
    for digits in a.digits:
        t0 = time.time()
        rec = dl.run_one(q, e1, a.mu, a.nu, digits, a.primes)
        rec["seconds"] = round(time.time() - t0)
        if rec.get("status") in ("ok", "indicator not recognised", "indicator recognised in more than one way"):
            break
        print(f"  status at {digits} digits: {rec.get('status')}", flush=True)
    rec["route"] = "orbits of ideals above primes l = 1 mod 4 (note 24)"
    with open(f_out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    if "surfaces" not in rec and rec.get("status") == "class polynomial not obtained":
        return
    # the routes of note 25: its code compares with "the stored orbits of note 24"; give it the record just made,
    # or, if the indicator was not recognised, a record with no surfaces (its comparisons then say False)
    stub = dict(rec)
    stub.setdefault("surfaces", [])
    jl0 = rv.jl
    rv.jl = lambda name: [stub] if name == "sixlines_delta_labels.jsonl" else jl0(name)
    ns = argparse.Namespace(q=q, e1=e1, mu=a.mu, nu=a.nu, digits=rec["digits"], workers=a.workers, primes=a.frob_primes,
                            out=f_out + ".galois")
    rv.cmd_galois(ns)
    # the split of the surfaces by each route, written out by configuration (the note 25 code prints roots only)
    fld = cm.CMField(q, e1, a.mu, a.nu)
    tor = cm.Torsor(fld, rec["digits"], verbose=False)
    red = cm.reduced_points(tor, q, e1 % q)
    enum_set = cm.enumerated_points(q, e1, tor.e2, tor.chart)
    with open(f_out + ".galois") as fh:
        gal = [json.loads(line) for line in fh][-1]
    out = {"q": q, "e1": e1, "route": "factors over E(i) and the form modulo 4 of xi (note 25)", "digits": rec["digits"],
           "chart": tor.chart, "coordinate": tor.mix, "galois": gal,
           "roots": {str(r): {"invariants": list(red[r]), "configurations": [list(c) for c in enum_set.get(red[r], [])]}
                     for r in sorted(red)}}
    # the form modulo 4 again, root by root (cmd_galois keeps only the values)
    if gal.get("form_indicator_recognised") and len(gal["form_indicator_recognised"]) == 1:
        pari = cm.pari
        s = {T: gal["form_mod_4_on_lattices"][str(T)] for T in cm.TYPES}
        I = {T: cm.interpolation_poly(tor.xs[T], s[T]) for T in cm.TYPES}
        name, e2, e3 = gal["form_indicator_recognised"][0]
        div = pari(1) if name == "1" else pari("I")
        Gam = tor.rec.poly([I[cm.TYPES[0]] / div, e2 * I[cm.TYPES[1]] / div, e3 * I[cm.TYPES[2]] / div])
        EPq = tor.prime(q, e1 % q)
        Gq, Hq = cm.reduce_poly(Gam, EPq), cm.reduce_poly(tor.H, EPq)
        dH = pari.deriv(Hq)
        for r in red:
            out["roots"][str(r)]["form_value"] = int(pari.lift(pari.subst(Gq, "X", pari.Mod(r, q)) / pari.subst(dH, "X", pari.Mod(r, q))))
    if "roots_of_one_factor" in gal:
        for side, part in enumerate(gal["roots_of_one_factor"]):
            for r in part:
                out["roots"][str(r)]["factor_side"] = side
    with open(f_out, "a") as fh:
        fh.write(json.dumps(out) + "\n")
    os.remove(f_out + ".galois")


# ------------------------------------------------------------------------------------------------ compare

def cmd_compare(a):
    recs = jl(F_DELTA)
    cos = jl(F_COSETS)
    result = {}
    for e in CASES:
        e1 = e[0]
        surf, rows = surfaces_of(e)
        where = {cfg: si for si, s in enumerate(surf) for cfg in s}
        mine = [r for r in recs if r["e"] == e]
        prim = [r for r in mine if r["reading"] == "primary"]
        st = collections.Counter(r["status"][:50] for r in prim)
        print(f"\nq = {P}, e = {e}: {len(surf)} surfaces, {sum(len(s) for s in surf)} configurations")
        print(f"  primary reading: {len(prim)} curves; by status: {dict(sorted(st.items()))}")
        ok = [r for r in prim if r["status"] == "ok"]
        dl_ = collections.defaultdict(list)
        ul_ = collections.defaultdict(list)
        for r in ok:
            dl_[r["surface"]].append(r["delta"])
            ul_[r["surface"]].append(r["u"])
        with_delta = sum(1 for si in range(len(surf)) if dl_[si])
        const = sum(1 for si in range(len(surf)) if len(set(dl_[si])) == 1)
        uconst = sum(1 for si in range(len(surf)) if len(set(ul_[si])) == 1)
        print(f"  surfaces with a delta: {with_delta}; delta the same on all pairs: {const}; u the same on all pairs: {uconst}; "
              f"pairs per surface: {sorted(len(dl_[si]) for si in range(len(surf)))}")
        print(f"  M over the curves: {dict(collections.Counter(r['M'] for r in ok))}; order of X: "
              f"{dict(collections.Counter(r['order_X'] for r in ok))}; u constant over the kernel vectors: "
              f"{sum(r['u_constant_over_kernel'] for r in ok)} of {len(ok)}; curves with unsolvable kernel vectors: "
              f"{sum(1 for r in ok if r['unsolvable'])}")
        delta = {si: (dl_[si][0] if len(set(dl_[si])) == 1 else None) for si in range(len(surf))}
        print(f"  surfaces with delta = +1, -1, mixed or none: {sum(1 for v in delta.values() if v == 1)}, "
              f"{sum(1 for v in delta.values() if v == -1)}, {sum(1 for v in delta.values() if v is None)}")
        # controls
        pu = {(tuple(r["cfg"]), tuple(r["pair"])): r["u"] for r in ok}
        for name in ("twist1", "twist2", "seed", "n16", "r2"):
            cs = [r for r in mine if r["reading"] == name]
            if not cs:
                continue
            good = [r for r in cs if r["status"] == "ok"]
            same = sum(1 for r in good if pu.get((tuple(r["cfg"]), tuple(r["pair"]))) == r["u"])
            print(f"  control {name}: {len(cs)} runs, {len(good)} ok, same u as the primary reading: {same} of {len(good)}; "
                  f"other statuses: {dict(collections.Counter(r['status'][:50] for r in cs if r['status'] != 'ok'))}")
        cf = [r for r in mine if r["reading"] == "cofactor"]
        if cf:
            good = [r for r in cf if r["status"] == "ok"]
            agree = sum(1 for r in good if delta.get(r["surface"]) == r["delta"])
            known = sum(1 for r in good if delta.get(r["surface"]) is not None)
            print(f"  extension (cofactor): {len(cf)} runs, {len(good)} ok; delta equals that of the surface: {agree} of {known}; "
                  f"statuses: {dict(collections.Counter(r['status'][:50] for r in cf))}")
        # the cosets
        res = {"surfaces": len(surf), "with_delta": with_delta, "delta_constant": const}
        lab = [r for r in cos if r.get("q") == P and r.get("e1") == e1 and "ideals" in r.get("route", "") and r.get("surfaces")]
        gal = [r for r in cos if r.get("q") == P and r.get("e1") == e1 and "E(i)" in r.get("route", "")]
        routes = {}
        if lab:
            routes["orbits of ideals"] = {frozenset(tuple(c) for c in s["configurations"]): s["orbit"] for s in lab[-1]["surfaces"]}
        if gal:
            g = gal[-1]["roots"]
            if all("factor_side" in v for v in g.values()):
                routes["factors over E(i)"] = {frozenset(tuple(c) for c in v["configurations"]): v["factor_side"] for v in g.values()}
            if all("form_value" in v for v in g.values()):
                v0 = min(v["form_value"] for v in g.values())
                routes["xi modulo 4"] = {frozenset(tuple(c) for c in v["configurations"]): int(v["form_value"] != v0) for v in g.values()}
        table = []
        for si, s in enumerate(surf):
            row = {"surface": si, "configurations": [list(c) for c in s], "delta": delta[si], "pairs": len(dl_[si])}
            for name, mp in routes.items():
                row[name] = mp.get(frozenset(s))
            table.append(row)
        for name, mp in routes.items():
            assert set(mp) == {frozenset(s) for s in surf}, f"{name}: the roots are not the enumerated surfaces"
            by = {0: [], 1: []}
            for row in table:
                by[row[name]].append(row["delta"])
            sizes = [len(by[0]), len(by[1])]
            vals = [sorted(set(x for x in by[k] if x is not None)) for k in (0, 1)]
            agrees = len(vals[0]) == 1 and len(vals[1]) == 1 and vals[0] != vals[1] and None not in by[0] + by[1]
            # the number of surfaces that agree with the better of the two assignments of signs to cosets
            best = max(sum(1 for row in table if row["delta"] == (sg if row[name] == 0 else -sg)) for sg in (1, -1))
            print(f"  route '{name}': cosets of sizes {sizes}; delta on coset 0: {dict(collections.Counter(by[0]))}, "
                  f"on coset 1: {dict(collections.Counter(by[1]))}; surfaces that agree: {best} of {len(surf)}; "
                  f"delta is constant on each coset and opposite: {agrees}")
            res[name] = {"sizes": sizes, "delta_by_coset": [dict(collections.Counter(str(x) for x in by[k])) for k in (0, 1)],
                         "agree": best, "holds": agrees}
        names = list(routes)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                same = {frozenset(frozenset(s) for s in surf if routes[names[i]][frozenset(s)] == k) for k in (0, 1)} == \
                       {frozenset(frozenset(s) for s in surf if routes[names[j]][frozenset(s)] == k) for k in (0, 1)}
                print(f"  routes '{names[i]}' and '{names[j]}' give the same split: {same}")
                res[f"{names[i]} = {names[j]}"] = same
        print("  surface  configurations                       delta  pairs  " + "  ".join(names))
        for row in table:
            print(f"  {row['surface']:5d}    {str(row['configurations']):38s} {str(row['delta']):>5s}  {row['pairs']:3d}    "
                  + "  ".join(str(row.get(nm)) for nm in names))
        res["table"] = table
        result[str(e1)] = res
    with open(os.path.join(HERE, F_COMPARE), "w") as fh:
        json.dump(result, fh, indent=1)


def cmd_signs(a):
    import cypari2
    pari = cypari2.Pari()

    def zeta(p):
        return pow(int(pari(f"lift(znprimroot({p}))")), (p - 1) // 4, p)

    cases = []
    rows = jl("sixlines_delta41_signs.jsonl")
    for i in range(0, len(rows), 2):
        lab, gal = rows[i], rows[i + 1]
        d = {s["root"]: s["delta"] for s in lab["surfaces"]}
        cases.append((lab["q"], lab["e1"], [(v["form_value"], d[int(r)][0]) for r, v in gal["roots"].items()]))
    cos = jl(F_COSETS)
    with open(os.path.join(HERE, F_COMPARE)) as fh:
        cmp_ = json.load(fh)
    for e in CASES:
        gal = [r for r in cos if r.get("e1") == e[0] and "E(i)" in r.get("route", "")][-1]
        tab = {frozenset(tuple(c) for c in row["configurations"]): row["delta"] for row in cmp_[str(e[0])]["table"]}
        cases.append((P, e[0], [(v["form_value"], tab[frozenset(tuple(c) for c in v["configurations"])])
                                for v in gal["roots"].values()]))
    good = 0
    for q, e1, pairs in sorted(cases):
        z = zeta(q)
        t = collections.Counter(pairs)
        rule = all((z * v - dlt) % q == 0 for v, dlt in pairs)
        anti = all((z * v + dlt) % q == 0 for v, dlt in pairs)
        good += rule
        print(f"q = {q}, e1 = {e1}: zeta = {z}; (value of the form indicator, delta): surfaces = {dict(sorted(t.items()))}; "
              f"delta = zeta * value on every surface: {rule}; delta = -zeta * value on every surface: {anti}")
    print(f"delta = zeta * value: {good} torsors of {len(cases)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    s = sub.add_parser("prym")
    s.add_argument("--pairs", type=int, default=0)
    s.add_argument("--workers", type=int, default=14)
    sub.add_parser("dyadic")
    s = sub.add_parser("validate")
    s.add_argument("--n", type=int, default=2)
    s = sub.add_parser("selfcheck")
    s.add_argument("--n", type=int, default=2)
    s.add_argument("--workers", type=int, default=4)
    s = sub.add_parser("delta")
    s.add_argument("--reading", choices=sorted(READINGS), default="primary")
    s.add_argument("--pairs", type=int, default=0, help="primary, cofactor: pairs per configuration (0 = all stored)")
    s.add_argument("--per-surface", type=int, default=1, help="controls: curves per surface")
    s.add_argument("--e1", type=int, nargs="*")
    s.add_argument("--surfaces", nargs="*", help="only these surfaces, as e1:index (the index of `list`)")
    s.add_argument("--workers", type=int, default=14)
    s = sub.add_parser("cosets")
    for name in ("q", "e1", "mu", "nu"):
        s.add_argument(name, type=int)
    s.add_argument("--digits", type=int, nargs="+", default=[2000])
    s.add_argument("--primes", type=int, default=8)
    s.add_argument("--frob-primes", type=int, default=40)
    s.add_argument("--workers", type=int, default=14)
    s.add_argument("--no-delta", action="store_true", help="do not read delta (the comparison is then left to `compare`)")
    s.add_argument("--out", help="write the records here and not to sixlines_delta41_cosets.jsonl (for a test at another q)")
    sub.add_parser("compare")
    sub.add_parser("signs")
    a = ap.parse_args()
    {"list": cmd_list, "prym": cmd_prym, "dyadic": cmd_dyadic, "validate": cmd_validate, "selfcheck": cmd_selfcheck,
     "delta": cmd_delta,
     "cosets": cmd_cosets, "compare": cmd_compare, "signs": cmd_signs}[a.cmd](a)


if __name__ == "__main__":
    main()
