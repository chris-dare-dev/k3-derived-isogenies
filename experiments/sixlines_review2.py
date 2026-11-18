#!/usr/bin/env python
"""
sixlines_review2.py -- the re-derivation checks of notes 20, 23 and 24 (note 25).

Nothing here changes the output of the checked scripts; their functions are imported and their stored files are
read.

  stored                      [note 24 section 1, note 23 section 5]  the stored values of delta
                              (sixlines_jacobian_character*.jsonl) against the stored orbits
                              (sixlines_delta_labels.jsonl): is the split into two orbits reproduced by anything
                              simpler than delta -- one factor of delta, a pair of factors, the flag "the
                              configuration, not its twist, has the polynomial P", the comparison flag of note 14,
                              self-association, the kind of branch points, the value M?  Also: surfaces are a
                              configuration and its associate on both sides; the values of M over the kernel of one
                              curve never contain both 4 and 12.
  galois q e1 mu nu           [note 24 sections 2, 3]  the two orbits without the action on lattices: the class
                              polynomial H over E (computed as in note 20) is factored over E and over E(i) by
                              PARI; over E(i) it should split into two factors of degree N/2, and the roots of one
                              factor modulo a prime above (q, F - e1) should be one of the stored orbits.  Second,
                              the form modulo 4 read on the lattices themselves: for each even lattice (L, xi) the
                              unit part modulo 4 of the image of -i xi in Q_2 at the prime of degree 1 above 2; its
                              indicator polynomial is recognised over E and reduced as in note 24 section 3, with
                              no ideal acting.  Then, for
                              primes l of both classes modulo 4 and every prime of degree 1 of E above l: the class
                              in G (PARI) and the degrees of the factors of H modulo that prime.  Checks that the
                              classes with l = 1 mod 4 form a subgroup of index 2 that the classes with l = 3 mod 4
                              avoid, and that Frobenius has cycles of the length of the order of the class.
  theta                       [note 20 section 2]  the sixteen theta constants of Matsumoto-Terasoma (19) by a direct
                              fourfold sum, against cm_theta16 and cm_tsq of sixlines_cmmethod.py.
  pairing p b c e f i j       [note 23 section 3]  logged checks of the Tate pairing of sixlines_jacobian.py over
                              F_{p^4}: bilinearity on both sides, T(P, N R) = 1, invariance under J, and
                              T(pi P, pi Q) = T(P, Q)^p.
  delta [--tw k] [--seed s]   [note 23 section 5]  delta recomputed on three curves of each torsor (one surface of
                              each orbit and one more), on another twist and with other random classes, and with the
                              2-adic arithmetic redone here: the root of chi in Z_2[i] by lifting digit by digit
                              (not by factorpadic), N / (pi_1^4 - 1) and the solution M from the two pairings of the
                              point computed directly.  Compared with the stored value of the surface.
  levels                      [note 23 section 4]  recount of sixlines_jacobian_levels.jsonl: surfaces by number of
                              blocks, and how many of them are self-associated.

Run from experiments/ with python -I.  cypari2 only.
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


def jl(name):
    with open(os.path.join(HERE, name)) as fh:
        return [json.loads(line) for line in fh]


def sgn(x):
    return 1 if x > 0 else -1


def surface_e(S):
    s1, s2, s3 = S
    return [s1, (s1 * s1 - s2) // 2, (s1 ** 3 - 3 * s1 * s2 + 2 * s3) // 6]


TORSORS = [(29, 2), (37, 2), (37, 18)]


def stored_tables():
    recs = []
    for name in ("sixlines_jacobian_character.jsonl", "sixlines_jacobian_character_extra.jsonl"):
        for r in jl(name):
            r["file"] = name
            recs.append(r)
    labels = {}
    for r in jl("sixlines_delta_labels.jsonl"):
        labels[(r["q"], r["e1"])] = r
    rows = {}
    for p in (29, 37):
        for r in jl(f"sixlines_p{p}.jsonl"):
            rows[(p, tuple(r["cfg"]))] = r
    prym = {}
    for r in jl("sixlines_prym.jsonl"):
        if r["p"] in (29, 37):
            for x in r["pairs"]:
                prym[(r["p"], tuple(r["cfg"]), tuple(x["pair"]))] = x
    return recs, labels, rows, prym


# ------------------------------------------------------------------------------------------------ stored

def cmd_stored(a):
    recs, labels, rows, prym = stored_tables()
    print(f"records of note 23: {len(recs)}; by status: {dict(collections.Counter(r['status'] for r in recs))}")
    main_ok = [r for r in recs if r["status"] == "ok" and r.get("kind") in (None, "more pairs")]
    print(f"computed curves that are not controls: {len(main_ok)}; controls with status ok: "
          f"{sum(1 for r in recs if r['status'] == 'ok') - len(main_ok)}")
    total_bad = 0
    for q, e1 in TORSORS:
        lab = labels[(q, e1)]
        orbit_of, root_of = {}, {}
        for s in lab["surfaces"]:
            cfgs = [tuple(c) for c in s["configurations"]]
            # a surface is a configuration and its associate, on both sides
            assoc = {tuple(rows[(q, c)]["assoc"]) for c in cfgs}
            assert assoc == set(cfgs), (cfgs, assoc)
            for c in cfgs:
                orbit_of[c] = s["orbit"]
                root_of[c] = s["root"]
        ok = [r for r in recs if r["status"] == "ok" and r["p"] == q and r["e"][0] == e1]
        print(f"\nq = {q}, e1 = {e1}: {len(lab['surfaces'])} surfaces, {len(ok)} records with status ok "
              f"({sum(1 for r in ok if r.get('kind') in (None, 'more pairs'))} curves, the rest controls)")
        curves = []
        for r in ok:
            cfg = tuple(r["cfg"])
            row = rows[(q, cfg)]
            c1, c2 = r["c"][0], r["c"][1]
            x = prym[(q, cfg, tuple(r["pair"]))]
            curves.append({
                "root": root_of[cfg], "orbit": orbit_of[cfg], "delta": r["delta"], "u": 1 if r["u"] == 1 else -1,
                "s1": sgn(c1[0] * c1[1]), "s2": sgn(c2[0]), "f": 1 if surface_e(row["S"]) == r["e"] else -1,
                "compare": x["compare"], "rat": x["rational_point"], "M": r["M"], "order_X": r["order_X"],
                "self": tuple(row["assoc"]) == cfg, "kernel": r["kernel"], "cfg": cfg})
        # delta is constant on each surface
        per = collections.defaultdict(set)
        for c in curves:
            per[c["root"]].add(c["delta"])
        print(f"  surfaces on which all records give one delta: {sum(1 for v in per.values() if len(v) == 1)} of {len(per)}")
        # candidates, curve by curve
        cands = [("delta", lambda c: c["delta"]), ("u alone", lambda c: c["u"]),
                 ("sign(Re c1 Im c1)", lambda c: c["s1"]), ("sign(c2)", lambda c: c["s2"]),
                 ("sign(Re c1 Im c1) sign(c2)", lambda c: c["s1"] * c["s2"]),
                 ("u sign(Re c1 Im c1)", lambda c: c["u"] * c["s1"]), ("u sign(c2)", lambda c: c["u"] * c["s2"]),
                 ("configuration has P", lambda c: c["f"]), ("comparison flag of note 14", lambda c: c["compare"]),
                 ("u * (configuration has P)", lambda c: c["u"] * c["f"]),
                 ("self-associated", lambda c: c["self"]), ("branch points (rational_point)", lambda c: c["rat"]),
                 ("M", lambda c: c["M"]), ("order of X", lambda c: c["order_X"])]
        print("  candidate                          value: (records on orbit 0, on orbit 1)        separates the orbits")
        for name, fn in cands:
            tab = collections.defaultdict(lambda: [0, 0])
            for c in curves:
                tab[fn(c)][c["orbit"]] += 1
            sep = all(v[0] == 0 or v[1] == 0 for v in tab.values()) and len(tab) > 1
            print(f"  {name:34s} {'; '.join(f'{k}: {tuple(v)}' for k, v in sorted(tab.items(), key=str)):46s} {'YES' if sep else 'no'}")
            if sep and name != "delta":
                total_bad += 1
        # candidates on surfaces: the flag "has P" of the two configurations of the surface
        tab = collections.defaultdict(lambda: [0, 0])
        for s in lab["surfaces"]:
            fs = tuple(sorted(1 if surface_e(rows[(q, tuple(c))]["S"])[0] == e1 else -1 for c in s["configurations"]))
            tab[fs][s["orbit"]] += 1
        print(f"  by surface, the flags 'has P' of its two configurations: "
              f"{'; '.join(f'{k}: {tuple(v)}' for k, v in sorted(tab.items()))}")
        # the values of M over the kernel of one curve
        mixed = 0
        for c in curves:
            # M is known modulo 16, so its odd part modulo 4 is known only when M is not divisible by 8
            odd = {(M >> (M & -M).bit_length() - 1) % 4 for o, sol, k in c["kernel"] for M in sol
                   if M % 8 and len(sol) == 1}
            mixed += len(odd) > 1
        print(f"  records whose kernel vectors give both classes of odd part of M modulo 4: {mixed} of {len(curves)}")
        total_bad += mixed
    print(f"\nsomething other than delta separates the orbits, or a curve is inconsistent: {total_bad} cases")


# ------------------------------------------------------------------------------------------------ galois

def cmd_galois(a):
    import sixlines_delta_labels as dl
    cm = dl.cm
    pari = cm.pari
    dl.WORKERS = a.workers
    q, e1 = a.q, a.e1
    t0 = time.time()
    fld = cm.CMField(q, e1, a.mu, a.nu)
    tor = cm.Torsor(fld, a.digits, verbose=False)
    N = tor.N
    if not tor.setup(q, e1 % q, verbose=False) and not dl.setup_more(tor, q, e1 % q):
        print("class polynomial not obtained")
        return
    print(f"q = {q}, e1 = {e1}: {N} even lattices per type, class polynomial recognised at {a.digits} digits, "
          f"chart {tor.chart}, coordinate {tor.mix}  ({time.time() - t0:.0f} s)", flush=True)
    ES = tor.ES
    bnfE, P = ES[0], ES[3]
    # H as a polynomial in x over Q[y]/P(y)
    Py = pari.subst(P, "w", "y")
    coeffs = [pari.subst(pari.lift(pari.nfbasistoalg(bnfE, c)), "w", "y") for c in tor.H]
    Hx = pari("(c, Py) -> sum(k = 1, #c, Mod(c[k], Py) * 'x^(k - 1))")(coeffs, Py)
    nfE = pari.nfinit(Py)
    fa = pari.nffactor(nfE, Hx)
    degs = [int(pari.poldegree(g)) for g in fa[0]]
    print(f"  H over E: factors of degree {degs}  ({time.time() - t0:.0f} s)", flush=True)
    # E(i)
    eq = pari("(nf) -> rnfequation(nf, 'x^2 + 1, 1)")(nfE)
    Qy = pari.subst(eq[0], "x", "y")
    ay = pari.subst(pari.lift(eq[1]), "x", "y")
    k = int(eq[2])
    Habs = pari("(c, ay, Qy) -> sum(k = 1, #c, subst(c[k], 'y, Mod(ay, Qy)) * 'x^(k - 1))")(coeffs, ay, Qy)
    nfA = pari.nfinit(Qy)
    fa2 = pari.nffactor(nfA, Habs)
    degs2 = [int(pari.poldegree(g)) for g in fa2[0]]
    print(f"  H over E(i): factors of degree {degs2}  ({time.time() - t0:.0f} s)", flush=True)
    rec = {"q": q, "e1": e1, "N": N, "digits": a.digits, "factors_over_E": degs, "factors_over_E_i": degs2}
    if sorted(degs2) == [N // 2, N // 2]:
        # reduce one factor at a prime above (q, F - e1): y_abs = i + k F  ->  i0 + k e1
        red = cm.reduced_points(tor, q, e1 % q)
        enum_set = cm.enumerated_points(q, e1, tor.e2, tor.chart)
        i0 = int(pari(f"lift(sqrt(Mod(-1, {q})))"))
        parts = []
        redfac = pari(r"""(nfA, g, Fa, ia, e1, i0, q) -> {
            my(prs = idealprimedec(nfA, q), pr = 0, md, v = Vec(g), y);
            for(j = 1, #prs, if(prs[j].f == 1 && prs[j].e == 1 && nfeltval(nfA, Fa - e1, prs[j]) > 0
                                && nfeltval(nfA, ia - i0, prs[j]) > 0, pr = prs[j]));
            if(pr === 0, error("no prime of degree 1 above (q, F - e1, i - i0)"));
            md = nfmodprinit(nfA, pr);
            Pol(vector(#v, j, if(nfeltval(nfA, v[j], pr) < 0, error("not integral"));
                              y = nfmodpr(nfA, v[j], md); if(type(y) == "t_FFELT", y = y.pol);
                              Mod(lift(polcoef(y + 0 * 'x, 0)), q)), 'X) }""")
        Fa = pari("(ay, Qy) -> Mod(ay, Qy)")(ay, Qy)
        ia = pari("(ay, Qy, k) -> Mod('y, Qy) - k * Mod(ay, Qy)")(ay, Qy, k)
        assert ia * ia == -1
        for ii in (i0, q - i0):
            h1 = redfac(nfA, fa2[0][0], Fa, ia, e1, ii, q)
            rts = sorted(int(pari.lift(r)) for r in pari.polrootsmod(h1, q))
            parts.append(rts)
        assert sorted(parts[0] + parts[1]) == sorted(red), (parts, sorted(red))
        mine = {frozenset(frozenset(tuple(c) for c in enum_set[red[r]]) for r in part) for part in parts}
        lab = [r for r in jl("sixlines_delta_labels.jsonl") if r["q"] == q and r["e1"] == e1][-1]
        stored = {frozenset(frozenset(tuple(c) for c in s["configurations"]) for s in lab["surfaces"] if s["orbit"] == o)
                  for o in (0, 1)}
        dlt = {frozenset(frozenset(tuple(c) for c in s["configurations"]) for s in lab["surfaces"] if s["delta"] == [d])
               for d in (1, -1)}
        rec["roots_of_one_factor"] = parts
        rec["equals_stored_orbits"] = mine == stored
        rec["equals_delta_split"] = mine == dlt
        print(f"  roots modulo {q} of one factor over E(i), at the two primes above (q, F - e1): {parts}")
        print(f"  this split of the surfaces equals the stored orbits of note 24: {mine == stored}; "
              f"equals the split by delta: {mine == dlt}", flush=True)
    # the form modulo 4 read on the lattices: for (L, xi) the polarisation is Tr(xi x conj(y)), and on the line
    # K_{w1} = Q_2(i) of K (x) Q_2 it is h(m) = 2 N(m) t with t the image of -i xi in Q_2.  Its unit part modulo 4
    # is note 23's u for the lattice; a step with multiplier l multiplies it by l.
    from fractions import Fraction
    co = [(q * q, 0), (-4 * q * a.nu, 4 * q * a.mu), (e1, 0), (-4 * a.nu, -4 * a.mu), (1, 0)]
    B = 200
    cand = {(x[0] % 2 ** (B - 12), x[1] % 2 ** (B - 12)) for x in roots_by_lifting(co, B)}
    rec["roots_of_chi_in_Z2_i"] = len(cand)
    if len(cand) == 1:
        pi1 = next(iter(cand))
        iK = fld.F[1]

        def unit_mod4(xi):
            tp = pari.lift(-iK * xi)
            re, im, pw = Fraction(0), Fraction(0), (1, 0)
            for kk in range(int(pari.poldegree(tp)) + 1):
                c = pari.polcoef(tp, kk)
                c = Fraction(int(pari.numerator(c)), int(pari.denominator(c)))
                re += c * pw[0]
                im += c * pw[1]
                pw = gmul(pw, pi1)
            v2 = lambda fr: ((fr.numerator & -fr.numerator).bit_length() - 1) - ((fr.denominator & -fr.denominator).bit_length() - 1)
            vr = v2(re)
            vi = v2(im) if im else 10 ** 6
            assert vi - vr > B - 60, ("the image of -i xi is not in Q_2", vr, vi)
            num, den = re.numerator, re.denominator
            num //= (num & -num)
            den //= (den & -den)
            return (num * pow(den, -1, 4)) % 4, vr

        s = {}
        for T in cm.TYPES:
            us = [unit_mod4(p_["xi"]) for p_ in tor.by[T]]
            s[T] = [1 if u_[0] == 1 else -1 for u_ in us]
        rec["form_mod_4_on_lattices"] = {str(T): s[T] for T in cm.TYPES}
        print(f"  unit part modulo 4 of -i xi at the prime of degree 1 above 2, as a sign, on the lattices of each type: "
              f"{[s[T] for T in cm.TYPES]}", flush=True)
        if all(sum(s[T]) == 0 for T in cm.TYPES):
            I = {T: cm.interpolation_poly(tor.xs[T], s[T]) for T in cm.TYPES}
            found = []
            for name, div in (("1", pari(1)), ("i", pari("I"))):
                for e2 in (1, -1):
                    for e3 in (1, -1):
                        Gam = tor.rec.poly([I[cm.TYPES[0]] / div, e2 * I[cm.TYPES[1]] / div, e3 * I[cm.TYPES[2]] / div])
                        if Gam is not None:
                            found.append((name, e2, e3, Gam))
            rec["form_indicator_recognised"] = [[f[0], f[1], f[2]] for f in found]
            print(f"  its indicator polynomial is recognised over E for (divisor, signs) = {rec['form_indicator_recognised']}")
            if len(found) == 1:
                EPq = tor.prime(q, e1 % q)
                Gq, Hq = cm.reduce_poly(found[0][3], EPq), cm.reduce_poly(tor.H, EPq)
                dH = pari.deriv(Hq)
                red = cm.reduced_points(tor, q, e1 % q)
                enum_set = cm.enumerated_points(q, e1, tor.e2, tor.chart)
                val = {r: int(pari.lift(pari.subst(Gq, "X", pari.Mod(r, q)) / pari.subst(dH, "X", pari.Mod(r, q)))) for r in red}
                v0 = min(val.values())
                mine = {frozenset(frozenset(tuple(c) for c in enum_set[red[r]]) for r in red if (val[r] == v0) == side)
                        for side in (True, False)}
                lab = [r for r in jl("sixlines_delta_labels.jsonl") if r["q"] == q and r["e1"] == e1][-1]
                stored = {frozenset(frozenset(tuple(c) for c in s_["configurations"]) for s_ in lab["surfaces"] if s_["orbit"] == o)
                          for o in (0, 1)}
                dlt = {frozenset(frozenset(tuple(c) for c in s_["configurations"]) for s_ in lab["surfaces"] if s_["delta"] == [d])
                       for d in (1, -1)}
                rec["form_values_at_roots"] = sorted(set(val.values()))
                rec["form_split_equals_stored_orbits"] = mine == stored
                rec["form_split_equals_delta_split"] = mine == dlt
                print(f"  values at the roots modulo {q}: {sorted(set(val.values()))}; the split of the surfaces by the form "
                      f"on the lattice equals the stored orbits: {mine == stored}; equals the split by delta: {mine == dlt}",
                      flush=True)
    # classes of primes of degree 1 and Frobenius on the roots of H
    D = ES[2][2]
    cyc = [int(D[r, r]) for r in range(len(D)) if int(D[r, r]) > 1]
    order_G = math.prod(cyc)
    cls = {1: set(), 3: set()}
    tally = collections.Counter()
    bad = []
    used = {1: 0, 3: 0}
    l = 2
    while min(used.values()) < a.primes and l < 5000:
        l = int(pari.nextprime(l + 1))
        if l == q or l == 2:
            continue
        for b in pari.polrootsmod(P, l):
            b = int(pari.lift(b))
            EP = tor.prime(l, b)
            if EP == 0:
                continue
            Hl = cm.reduce_poly(tor.H, EP)
            if Hl is None or int(pari.poldegree(pari.gcd(Hl, pari.deriv(Hl)))) != 0:
                continue
            g = tuple(int(x) for x in cm.G("e_class")(ES, b, l))
            order = 1
            for gr, d in zip(g, cyc):
                order = math.lcm(order, d // math.gcd(gr, d))
            fdeg = sorted({int(pari.poldegree(f)) for f in pari.factormod(Hl, l)[0]})
            cls[l % 4].add(g)
            used[l % 4] += 1
            good = fdeg == [order]
            tally[(l % 4, good)] += 1
            if not good:
                bad.append((l, b, g, order, fdeg))
    sub1 = cls[1]
    closed = all(tuple((x + y) % d for x, y, d in zip(g1, g2, cyc)) in sub1 for g1 in sub1 for g2 in sub1)
    rec.update({"G": cyc, "classes_l_1_mod_4": len(cls[1]), "classes_l_3_mod_4": len(cls[3]),
                "disjoint": not (cls[1] & cls[3]), "classes_1_mod_4_closed": closed,
                "frobenius_cycles_equal_order": {f"l = {k[0]} mod 4, {'yes' if k[1] else 'no'}": v for k, v in tally.items()}})
    print(f"  G = {cyc} of order {order_G}; primes of degree 1 of E, l < {l + 1}: {used[1]} with l = 1 mod 4 in "
          f"{len(cls[1])} classes, {used[3]} with l = 3 mod 4 in {len(cls[3])} classes; the two sets of classes are "
          f"disjoint: {not (cls[1] & cls[3])}; the first is closed under addition: {closed}")
    print(f"  all factors of H modulo the prime have degree the order of its class in G: "
          f"{ {('l = %d mod 4' % k[0], k[1]): v for k, v in sorted(tally.items())} }")
    if bad:
        print(f"  exceptions (l, root, class, order, degrees): {bad[:10]}")
    print(f"  ({time.time() - t0:.0f} s)")
    if a.out:
        with open(a.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


# ------------------------------------------------------------------------------------------------ theta

def cmd_theta(a):
    import sixlines_cmmethod as cm
    pari = cm.pari
    cm.setprec(38)
    pari(r"""
    rv_theta(tau, a, b, B) = {
      my(s = 0, m, al = [a[1], a[2]] / (1 + I), be = [b[1], b[2]] / (1 + I), ex);
      for(x1 = -B, B, for(y1 = -B, B, for(x2 = -B, B, for(y2 = -B, B,
        m = [x1 + I * y1 + al[1], x2 + I * y2 + al[2]];
        ex = (m * tau * conj(m)~) / 2 + real(m * conj(be)~);
        s += exp(2 * Pi * I * ex)))));
      s;
    }
    """)
    worst16, worst10 = 0.0, 0.0
    for seed in (1, 2, 3):
        tau = cm.random_tau(seed)
        R2 = cm.R2_for(30)
        th = cm.G("cm_theta16")(tau, R2, 0)
        direct = {}
        for a1 in (0, 1):
            for a2 in (0, 1):
                for b1 in (0, 1):
                    for b2 in (0, 1):
                        v = pari("rv_theta")(tau, [a1, a2], [b1, b2], 7)
                        direct[(a1, a2, b1, b2)] = v
                        d = float(abs(v - th[8 * a1 + 4 * a2 + 2 * b1 + b2]))
                        worst16 = max(worst16, d)
        tsq = cm.G("cm_tsq")(tau, R2)
        ch = [(0, 0, 0, 0), (1, 1, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 1, 1, 0), (0, 0, 0, 1), (0, 0, 1, 1),
              (1, 0, 0, 1), (1, 1, 1, 1)]          # Matsumoto-Terasoma table 1 in the order t1..t10 of (2.1), typed here
        for kk in range(10):
            worst10 = max(worst10, float(abs(direct[ch[kk]] ** 2 - tsq[kk])))
        odd = max(float(abs(v)) for key, v in direct.items() if key not in ch)
        print(f"  point {seed}: the six characteristics outside the table vanish in the direct sum to {odd:.1e}")
    print(f"the sixteen theta constants, direct fourfold sum against cm_theta16, three points: largest difference {worst16:.1e}")
    print(f"the ten squares in the order t1..t10, direct sum against cm_tsq (quarter-plane sum): largest difference {worst10:.1e}")


# ------------------------------------------------------------------------------------------------ pairing

def cmd_pairing(a):
    import sixlines_jacobian as sj
    p, cfg, pair = a.p, [a.b, a.c, a.e, a.f], (a.i, a.j)
    c = None
    for x in sj.stored_pairs(p)[tuple(cfg)]:
        if tuple(x[0]) == pair:
            c = x[1]
    assert c is not None, "pair not stored by note 14"
    t0 = time.time()
    A = sj.Jac(p, cfg, pair, 0, 4, seed=a.seed)
    n = A.card(c)
    m = n >> sj.val(n, 2)
    N = 16
    ys = [A.mul(m, A.random()) for _ in range(3)]
    zs = [A.sub(y, A.J(y, 1)) for y in ys]
    assert all(A.is_zero(A.mul(N, z)) for z in zs), "exponent above 16: choose another pair"
    R = A.random()
    lg = lambda P_, Q_: A.tate_matrix([P_], [Q_], N)[0][0]
    out = []
    t = {(i, j): lg(zs[i], ys[j]) for i in range(2) for j in range(2)}
    out.append(("values T(z_a, y_b), a, b = 1, 2 (logarithms modulo 16)", [t[(0, 0)], t[(0, 1)], t[(1, 0)], t[(1, 1)]]))
    out.append(("T(z1 + z2, y1) = T(z1, y1) T(z2, y1)", lg(A.add(zs[0], zs[1]), ys[0]) == (t[(0, 0)] + t[(1, 0)]) % N))
    out.append(("T(z1, y1 + y2) = T(z1, y1) T(z1, y2)", lg(zs[0], A.add(ys[0], ys[1])) == (t[(0, 0)] + t[(0, 1)]) % N))
    out.append(("T(3 z1, y2) = T(z1, y2)^3", lg(A.mul(3, zs[0]), ys[1]) == 3 * t[(0, 1)] % N))
    out.append(("T(z1, 16 R) = 1 for a random class R", lg(zs[0], A.mul(N, R)) == 0))
    out.append(("T(z1, y2 + 16 R) = T(z1, y2)", lg(zs[0], A.add(ys[1], A.mul(N, R))) == t[(0, 1)]))
    out.append(("T(J z1, J y2) = T(z1, y2)", lg(A.J(zs[0], 1), A.J(ys[1], 1)) == t[(0, 1)]))
    out.append(("T(pi z1, pi y2) = T(z1, y2)^p", lg(A.frob(zs[0]), A.frob(ys[1])) == p * t[(0, 1)] % N))
    print(f"p = {p}, configuration {cfg}, pair {pair}, over F_(p^4), N = 16, seed {a.seed}  ({time.time() - t0:.0f} s)")
    for name, v in out:
        print(f"  {name}: {v}")


# ------------------------------------------------------------------------------------------------ delta

def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def geval(co, x, mod):
    r = (0, 0)
    for cf in reversed(co):
        r = gmul(r, x)
        r = ((r[0] + cf[0]) % mod, (r[1] + cf[1]) % mod)
    return r


def roots_by_lifting(co, B):
    """All x in (Z/2^B)[i] with chi(x) = 0 modulo 2^B, by lifting one binary digit at a time."""
    cur = [(x, y) for x in (0, 1) for y in (0, 1) if geval(co, (x, y), 2) == (0, 0)]
    for k in range(1, B):
        nxt = []
        for x in cur:
            for d0 in (0, 1):
                for d1 in (0, 1):
                    y = (x[0] + d0 * 2 ** k, x[1] + d1 * 2 ** k)
                    if geval(co, y, 2 ** (k + 1)) == (0, 0):
                        nxt.append(y)
        cur = nxt
        assert len(cur) < 5000, "too many approximate roots"
    return cur


def v1pi(x):
    """valuation at 1 + i of a non-zero Gaussian integer"""
    n = x[0] * x[0] + x[1] * x[1]
    return (n & -n).bit_length() - 1


def delta_job(job):
    import sixlines_jacobian as sj
    pari = sj.pari
    p, cfg, pair, c, tw, seed, want = job
    out = {"p": p, "cfg": list(cfg), "pair": list(pair), "c": c, "tw": tw, "seed": seed, "stored_delta": want}
    t0 = time.time()
    try:
        # chi of the twist, low degree first, written out here
        iu = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        c1, c2, c3 = [tuple(v) for v in c]
        c1t, c2t, c3t = gmul(iu[tw % 4], c1), gmul(iu[2 * tw % 4], c2), gmul(iu[3 * tw % 4], c3)
        co = [(p * p, 0), (-c3t[0], -c3t[1]), c2t, (-c1t[0], -c1t[1]), (1, 0)]
        B = 26
        rts = roots_by_lifting(co, B)
        classes = {(x[0] % 2 ** 16, x[1] % 2 ** 16) for x in rts}
        out["roots_mod_2^16"] = len(classes)
        if len(classes) != 1:
            out["status"] = "not one root in Z_2[i]"
            return out
        pi1 = next(iter(classes))
        mod = 2 ** 16
        dco = [(k * cf[0], k * cf[1]) for k, cf in enumerate(co)][1:]
        out["v_chi_prime"] = v1pi(geval(dco, pi1, mod))            # valuation at 1 + i of chi'(pi_1)
        out["pi1_minus_1_mod_2"] = ((pi1[0] - 1) % 2, pi1[1] % 2)
        x = (1, 0)
        for _ in range(4):
            x = gmul(x, pi1)
            x = (x[0] % mod, x[1] % mod)
        varpi = ((x[0] - 1) % mod, x[1])
        s = v1pi(varpi)
        N = 16
        assert (p ** 4 - 1) % 16 == 0 and (p ** 4 - 1) % 32 != 0
        if s > 8:
            out["status"] = "N / varpi not integral"
            return out
        # conj(N / varpi) = N varpi / Norm(varpi)
        nrm = varpi[0] * varpi[0] + varpi[1] * varpi[1]
        odd = nrm >> s
        assert odd % 2 == 1
        num = (N * varpi[0], N * varpi[1])
        assert num[0] % 2 ** s == 0 and num[1] % 2 ** s == 0
        inv = pow(odd % N, -1, N)
        wb = ((num[0] >> s) * inv % N, (num[1] >> s) * inv % N)
        al, be = pi1[0] % N, pi1[1] % N

        A = sj.Jac(p, list(cfg), tuple(pair), tw, 4, seed=seed)
        n = A.card(c)
        m = n >> sj.val(n, 2)
        mm = 12
        ys = [A.mul(m, A.random()) for _ in range(mm)]
        zs = [A.sub(y, A.J(y, 1)) for y in ys]
        if not all(A.is_zero(A.mul(N, z)) for z in zs):
            out["status"] = "exponent above 16"
            return out
        G0 = A.tate_matrix(zs, ys, N)
        G3 = A.tate_matrix([A.J(z, 3) for z in zs], ys, N)
        Gp = A.tate_matrix([A.frob(z) for z in zs], ys, N)
        K = [[(Gp[a_][b_] - al * G0[a_][b_] - be * G3[a_][b_]) % N for b_ in range(mm)] for a_ in range(mm)]
        ker = pari("(K, N) -> my(m = #K); matkerint(matconcat([mattranspose(K), N * matid(m)]))[1..m, ]")(
            pari.matrix(mm, mm, [v for row in K for v in row]), N)
        vecs = [[int(ker[i, j]) % N for i in range(mm)] for j in range(len(ker))]
        # candidates: T(X, X) from the matrices, used only to choose; the last one with an odd multiple of 4 or 2
        def txx(cv):
            return sum(cv[a_] * cv[b_] * (G0[a_][b_] - G3[a_][b_]) for a_ in range(mm) for b_ in range(mm)) % N
        cand = [cv for cv in vecs if any(cv) and txx(cv) % 8 != 0]
        if not cand:
            cand = [cv for cv in vecs if any(cv) and txx(cv) != 0]
        if not cand:
            out["status"] = "no kernel vector with a non-zero self-pairing"
            return out
        cv = cand[-1]
        X = A.zero
        for a_ in range(mm):
            if cv[a_]:
                X = A.add(X, A.mul(cv[a_], zs[a_]))
        iX = A.J(X, 3)
        chk = {"pi X = pi_1 X": A.is_zero(A.sub(A.frob(X), A.add(A.mul(al, X), A.mul(be, iX)))),
               "J^2 X = -X": A.is_zero(A.add(X, A.J(X, 2))), "16 X = 0": A.is_zero(A.mul(16, X)),
               "4 X != 0": not A.is_zero(A.mul(4, X))}
        Z = (A.tate_matrix([X], [iX], N)[0][0], A.tate_matrix([X], [X], N)[0][0])
        sol = [M for M in range(N) if (M * wb[0] - Z[0]) % N == 0 and (M * wb[1] - Z[1]) % N == 0]
        sol_other = [M for M in range(N) if (M * wb[0] - Z[0]) % N == 0 and (-M * wb[1] - Z[1]) % N == 0]
        out.update({"s": s, "checks": chk, "Z": list(Z), "M": sol, "M_with_the_other_convention": sol_other,
                    "seconds": round(time.time() - t0)})
        if len(sol) != 1 or sol[0] == 0 or not all(chk.values()):
            out["status"] = "no single non-zero M, or a check failed"
            return out
        M = sol[0]
        v2 = (M & -M).bit_length() - 1
        if v2 > 2:
            out["status"] = "M known to too few digits"
            return out
        u = (M >> v2) % 4
        sg = sgn(c1[0] * c1[1]) * sgn(c2[0])
        out.update({"u": u, "delta": (1 if u == 1 else -1) * sg, "status": "ok"})
    except Exception as err:
        out["status"] = "error: " + str(err)[:200]
    return out


def cmd_delta(a):
    recs, labels, rows, prym = stored_tables()
    jobs = []
    for q, e1 in TORSORS:
        lab = labels[(q, e1)]
        picked, seen_orbit = [], collections.Counter()
        # the surfaces in the order of their roots; take the first of each orbit and then one more
        for s in sorted(lab["surfaces"], key=lambda s_: s_["root"]):
            if seen_orbit[s["orbit"]] >= (2 if s["orbit"] == 0 else 1):
                continue
            cfgs = {tuple(c) for c in s["configurations"]}
            ok = [r for r in recs if r["status"] == "ok" and r["p"] == q and tuple(r["cfg"]) in cfgs
                  and r.get("kind") in (None, "more pairs")]
            if not ok:
                continue
            r = ok[-1]                                    # the last computed pair of the surface
            seen_orbit[s["orbit"]] += 1
            picked.append((q, tuple(r["cfg"]), tuple(r["pair"]), r["c"], a.tw, a.seed, s["delta"][0]))
        jobs += picked
    print(f"{len(jobs)} curves, twist {a.tw}, seed {a.seed}", flush=True)
    with Pool(a.workers) as pool:
        res = pool.map(delta_job, jobs, chunksize=1)
    agree = 0
    for r in res:
        same = r.get("delta") == r["stored_delta"]
        agree += bool(same)
        print(f"  p = {r['p']}  {tuple(r['cfg'])} pair {tuple(r['pair'])}: status {r['status']}; "
              f"Z = {r.get('Z')}, M = {r.get('M')} (other convention: {r.get('M_with_the_other_convention')}), "
              f"v(chi'(pi_1)) = {r.get('v_chi_prime')}, s = {r.get('s')}, "
              f"u = {r.get('u')}, delta = {r.get('delta')}, stored for the surface {r['stored_delta']}: "
              f"{'agrees' if same else 'DOES NOT AGREE' if r['status'] == 'ok' else 'not computed'}; "
              f"checks {r.get('checks')}; {r.get('seconds')} s", flush=True)
    print(f"computed {sum(1 for r in res if r['status'] == 'ok')} of {len(res)}; delta agrees with the stored value of the surface: {agree}")
    if a.out:
        with open(a.out, "a") as fh:
            for r in res:
                fh.write(json.dumps(r) + "\n")


# ------------------------------------------------------------------------------------------------ levels

def cmd_levels(a):
    recs = jl("sixlines_jacobian_levels.jsonl")
    rows = {}
    for p in (29, 37):
        for r in jl(f"sixlines_p{p}.jsonl"):
            rows[(p, tuple(r["cfg"]))] = r
    cases = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in recs:
        cases[(r["p"], r["e"][0], r["l"])][r["surface"]].append(r)
    for (p, e1, l), surf in cases.items():
        two = one = mixed = 0
        self2 = self1 = 0
        pairs = 0
        for si, rs in surf.items():
            ok = [r for r in rs if r["status"] == "ok"]
            pairs += len(ok)
            bl = {r["blocks"] for r in ok}
            cfgs = {tuple(r["cfg"]) for r in rs}
            selfa = all(tuple(rows[(p, c)]["assoc"]) == c for c in cfgs)
            if bl == {2}:
                two += 1
                self2 += selfa
            elif bl == {1}:
                one += 1
                self1 += selfa
            else:
                mixed += 1
        print(f"q = {p}, e1 = {e1}, l = {l}: {len(surf)} surfaces, {pairs} pairs computed; two blocks {two} "
              f"(self-associated {self2}), one block {one} (self-associated {self1}), pairs disagree on {mixed}; "
              f"two blocks and not self-associated: {two - self2} of {len(surf) - self2 - self1}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("stored")
    s = sub.add_parser("galois")
    for name in ("q", "e1", "mu", "nu"):
        s.add_argument(name, type=int)
    s.add_argument("--digits", type=int, default=1200)
    s.add_argument("--workers", type=int, default=8)
    s.add_argument("--primes", type=int, default=40)
    s.add_argument("--out")
    sub.add_parser("theta")
    s = sub.add_parser("pairing")
    for name in ("p", "b", "c", "e", "f", "i", "j"):
        s.add_argument(name, type=int)
    s.add_argument("--seed", type=int, default=3)
    s = sub.add_parser("delta")
    s.add_argument("--tw", type=int, default=3)
    s.add_argument("--seed", type=int, default=2026)
    s.add_argument("--workers", type=int, default=8)
    s.add_argument("--out")
    sub.add_parser("levels")
    a = ap.parse_args()
    {"stored": cmd_stored, "galois": cmd_galois, "theta": cmd_theta, "pairing": cmd_pairing, "delta": cmd_delta,
     "levels": cmd_levels}[a.cmd](a)


if __name__ == "__main__":
    main()
