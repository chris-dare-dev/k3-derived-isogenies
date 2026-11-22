#!/usr/bin/env python
"""
sixlines_order3.py -- characters of order 3 of G against invariants of six-line surfaces (note 29).

Note 28 could not see characters of order 3: among its labelled torsors only two have 3 | |G|, both with
q = 2 mod 3.  This script

  cosets q e1 mu nu [--digits D ...]
        the three cosets of ker chi, chi a character of order 3 of G (G has cyclic 3-part on every torsor treated,
        so there is one such kernel).  The class polynomial H over E is computed with the code of notes 20
        and 24 (lattices, period points, recognition over E), unchanged except for the progress line of
        Torsor._setup_chart (patch_setup_chart) and a cache of the period points; PARI's bnrclassfield gives the cyclic cubic extension
        L of E cut out by ker chi (as a subgroup of Cl(E) containing the subgroup of note 13, Lemma 9, as
        sixlines_cmmethod.e_setup builds it); H is factored over L; each factor of degree N/3 is reduced at one prime
        of L of degree 1 above (q, F - e1).  Its roots are one coset.  Checks: the roots of H mod (q, F - e1) are
        exactly the enumerated surfaces; three factors of degree N/3 whose roots partition them; for the first
        primes l of E of degree 1, the relative cubic splits modulo l exactly when the class of the prime in G lies
        in ker chi (class field theory against e_class).  Appends to sixlines_order3_cosets.jsonl.
  labelcheck
        on the torsors with full labels in G (sixlines_cmmethod.jsonl, note 20 section 7.5) and 3 | |G|: the labels
        modulo the 3-part against the cosets of `cosets`.
  features
        the candidate invariants of every surface of every torsor with cosets (see note 29 section 3), written to
        sixlines_order3_features.json.
  audit
        the pre-registered tests of note 29 section 3: exact nulls, Holm over each population, controls, plants,
        shuffled labels (all in this one command).  Writes sixlines_order3.jsonl.
  fields
        for the ten polynomials of note 29 section 3: E has no quadratic subfield and x^2 + x + 1 is irreducible over E.
  selfcheck
        the exact table distribution against the full enumeration of partitions (note 29 section 7.3).

A coset labelling by Z/3 is determined up to translation and sign by the partition alone: every bijection of a set of
three elements onto Z/3 is affine.  So the partition is all a test of an order-3 character needs.

Fixed input names of the imported scripts (sixlines_cm.jsonl) are redirected in memory to sixlines_cm_p41.jsonl and
sixlines_cm_p43.jsonl at q = 41, 43.  No existing file is written.  Uses cypari2 and numpy.  Class groups are
PARI's without certification (GRH).  Run from experiments/ with python -I.
"""
import argparse
import collections
import itertools
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

F_COSETS = "sixlines_order3_cosets.jsonl"
F_FEATURES = "sixlines_order3_features.json"
F_OUT = "sixlines_order3.jsonl"


def jl(name):
    with open(os.path.join(HERE, name)) as fh:
        return [json.loads(line) for line in fh if line.strip()]


GPO3 = r"""
o3_subgroup(ES, r) = {
  my(bnf = ES[1], U = ES[3][1], k = #bnf.cyc, Lam = matid(k));
  Lam[r, r] = 3;
  mathnf(concat(U^(-1) * Lam, matdiagonal(bnf.cyc)));
}
o3_prime_above(nfA, Fa, e1, q) = {
  my(prs = idealprimedec(nfA, q), out = List());
  for(j = 1, #prs, if(prs[j].f == 1 && prs[j].e == 1 && nfeltval(nfA, Fa - e1, prs[j]) > 0, listput(out, prs[j])));
  Vec(out);
}
o3_reduce(nfA, g, pr, q) = {
  my(md = nfmodprinit(nfA, pr), v = Vec(g), y);
  Pol(vector(#v, j, if(nfeltval(nfA, v[j], pr) < 0, error("not integral"));
                    y = nfmodpr(nfA, v[j], md); if(type(y) == "t_FFELT", y = y.pol);
                    Mod(lift(polcoef(y + 0 * 'x, 0)), q)), 'X);
}
"""


def setup_cm(q):
    import sixlines_delta_labels as dl
    cm = dl.cm
    load0 = cm.load_surfaces

    def load_surfaces(q_, e1_, e2_):
        if q_ in (41, 43):
            return [r for r in jl(f"sixlines_cm_p{q_}.jsonl")
                    if r["p"] == q_ and abs(r["e"][0]) == abs(e1_) and r["e"][1] == e2_]
        return load0(q_, e1_, e2_)
    cm.load_surfaces = load_surfaces
    cm.load_gp(GPO3)
    patch_setup_chart(cm)
    cache_points(cm)
    return dl, cm


def patch_setup_chart(cm):
    """Torsor._setup_chart of sixlines_cmmethod with one change: the verbose line factors the lcm of the
    denominators of H only over primes below 10^6 (factor(d, 10^6)), and prints the cofactor's size.  The full
    factorisation in the original runs only after H has been recognised, and can run for hours when the lcm has a
    large composite factor: at q = 43, e1 = 14 it has one of 97 digits, and the first runs at q = 43 stalled there.
    Only that log line changes; the acceptance rule is unchanged."""
    pari = cm.pari

    def _setup_chart(self, q, EP, verbose):
        extra = [m for m in EXTRA_MIXES if len(m) == len(self.pts[0]["inv"])]
        for mix in extra + [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 1, 0, 0, 0], [1, 2, 3, 0, 0],
                            [3, 1, 2, 5, 7], [2, -1, 3, 1, -2], [5, 3, -1, 2, 1]]:
            t1 = time.time()
            xs = self.coordinate(mix)
            H = self.rec.poly([pari("v -> my(P = 1); for(k = 1, #v, P *= ('X - v[k])); P")(xs[T])
                               for T in cm.TYPES])
            if H is None:
                if verbose:
                    print(f"  chart {self.chart}, coordinate {mix}: not recognised at {self.digits} digits", flush=True)
                continue
            Hq = cm.reduce_poly(H, EP)
            nroots = -1 if Hq is None else len(pari.polrootsmod(Hq, q))
            sqf = Hq is not None and pari.poldegree(pari.gcd(Hq, pari.deriv(Hq))) == 0
            if verbose:
                fa = pari.factor(cm.denominators(H), 10 ** 6)
                big = [p for p in fa[0] if p > 10 ** 6]
                small = [p for p in fa[0] if p <= 10 ** 6]
                print(f"  chart {self.chart}, coordinate {mix}: denominators: primes below 10^6 {small}, "
                      f"cofactors of {[len(str(p)) for p in big]} digits; roots mod {q}: {nroots} (-1: not integral "
                      f"at the prime), squarefree: {bool(sqf)} ({time.time() - t1:.1f} s)", flush=True)
            if sqf and nroots == self.N:
                self.mix, self.xs, self.H = mix, xs, H
                self.Gs = self.invariant_polys()
                if any(g is None for g in self.Gs) or any(cm.reduce_poly(g, EP) is None for g in self.Gs):
                    if verbose:
                        print(f"  the polynomials of the invariants are not recognised at {self.digits} digits",
                              flush=True)
                    if IDENTIFY_BY_COORDINATE and mix in extra:
                        self.Gs = None
                        print("  surfaces identified by the coordinate instead (--identify-by-coordinate)", flush=True)
                        return True
                    self.mix = None
                    return False
                return True
        return False
    cm.Torsor._setup_chart = _setup_chart


# Coordinates tried before the eight of sixlines_cmmethod (cosets --mix).  With 24 surfaces in F_43 a coordinate is
# injective modulo q with probability 3e-4, so at q = 43, e1 = 2 the eight fixed and the 300 random coordinates all
# have repeated roots; --mix lets a coordinate found injective on the enumerated surfaces be tried first.  The
# acceptance rule is the same, and the check that the roots of H are the enumerated surfaces still applies.
EXTRA_MIXES = []
# cosets --identify-by-coordinate: if the class polynomial of a --mix coordinate is recognised with N simple roots
# modulo the prime but the polynomials of the other invariants are not, identify the roots with the enumerated
# surfaces by the value of the coordinate itself.  The coordinate is required to be injective on the enumerated
# surfaces, and the check becomes: the roots of H modulo (q, F - e1) are exactly the coordinate values of the
# enumerated surfaces.  The record says "identified_by": "coordinate".
IDENTIFY_BY_COORDINATE = False

POINT_CACHE = os.path.join(HERE, "sixlines_order3_points")


def cache_points(cm):
    """cm.cm_points (the parallel version of sixlines_delta_labels) with the period points kept on disk in
    POINT_CACHE, keyed by the field and the precision, so that a second run at the same precision does not recompute
    them.  Only tau, J and det Y are stored, as GP strings."""
    pari = cm.pari
    inner = cm.cm_points

    def cm_points(fld, digits, verbose=True, order=True, level=True):
        key = f"q{fld.q}_e{fld.e1}_mu{fld.mu}_nu{fld.nu}_d{digits}_o{int(order)}_l{int(level)}.json"
        path = os.path.join(POINT_CACHE, key)
        if os.path.exists(path):
            cm.setprec(digits + 20)
            lats = fld.order_lattices(False, level) if order else fld.lattices()
            with open(path) as fh:
                st = json.load(fh)
            assert len(st) == len(lats)
            for lat, s in zip(lats, st):
                lat["J"] = None
                if s is not None:
                    assert lat["T"] is not None
                    lat["tau"], lat["J"], lat["detY"] = pari(s[0]), pari(s[1]), s[2]
            print(f"  {sum(s is not None for s in st)} period points at {digits} digits read from {path}", flush=True)
            return lats
        lats = inner(fld, digits, verbose, order, level)
        os.makedirs(POINT_CACHE, exist_ok=True)
        st = [None if lat.get("J") is None else [str(lat["tau"]), str(lat["J"]), lat["detY"]] for lat in lats]
        with open(path + ".part", "w") as fh:
            json.dump(st, fh)
        os.replace(path + ".part", path)
        return lats
    cm.cm_points = cm_points


def cmd_cosets(a):
    dl, cm = setup_cm(a.q)
    pari = cm.pari
    dl.WORKERS = a.workers
    global IDENTIFY_BY_COORDINATE
    if a.mix:
        EXTRA_MIXES.append(list(a.mix))
    IDENTIFY_BY_COORDINATE = a.identify_by_coordinate
    q, e1 = a.q, a.e1
    for digits in a.digits:
        t0 = time.time()
        rec = {"q": q, "e1": e1, "mu": a.mu, "nu": a.nu, "digits": digits}
        fld = cm.CMField(q, e1, a.mu, a.nu)
        tor = cm.Torsor(fld, digits, verbose=False)
        N = tor.N
        rec["N"] = N
        print(f"q = {q}, e1 = {e1}: {N} even lattices per type, {digits} digits ({time.time() - t0:.0f} s)", flush=True)
        if not tor.setup(q, e1 % q, verbose=True) and not dl.setup_more(tor, q, e1 % q):
            rec["status"] = "class polynomial not obtained"
            rec["seconds"] = round(time.time() - t0)
            print(f"  class polynomial not obtained at {digits} digits", flush=True)
            with open(os.path.join(HERE, F_COSETS), "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            continue
        rec.update({"chart": tor.chart, "coordinate": tor.mix})
        print(f"  class polynomial recognised: chart {tor.chart}, coordinate {tor.mix} ({time.time() - t0:.0f} s)",
              flush=True)
        enum_set = cm.enumerated_points(q, e1, tor.e2, tor.chart)
        if tor.Gs is None:
            # identification by the coordinate: {coordinate value mod q: invariants of the enumerated surface}
            EPq = tor.prime(q, e1 % q)
            roots = {int(pari.lift(r)) for r in pari.polrootsmod(cm.reduce_poly(tor.H, EPq), q)}
            xval = {k: sum(m_ * v for m_, v in zip(tor.mix, k)) % q for k in enum_set if k is not None}
            injective = len(set(xval.values())) == len(xval) == len(enum_set)
            red = {x: k for k, x in xval.items()} if injective else {}
            rec["identified_by"] = "coordinate"
            rec["coordinate_injective_on_enumerated"] = injective
            rec["roots_mod_p"] = len(roots)
            rec["enumerated_surfaces"] = len(enum_set)
            rec["roots_equal_enumerated"] = injective and roots == set(red) and len(roots) == N
            print(f"  coordinate injective on the enumerated surfaces: {injective}", flush=True)
        else:
            red = cm.reduced_points(tor, q, e1 % q)
            rec["identified_by"] = "invariant polynomials"
            rec["roots_mod_p"] = len(red)
            rec["enumerated_surfaces"] = len(enum_set)
            rec["roots_equal_enumerated"] = set(red.values()) == set(enum_set) and len(red) == N
        print(f"  roots of H mod (q, F - e1): {len(red)}; enumerated surfaces: {len(enum_set)}; the same set: "
              f"{rec['roots_equal_enumerated']}", flush=True)
        ES = tor.ES
        D = ES[2][2]
        k = len(D)
        diag = [int(D[r, r]) for r in range(k)]
        cyc = [d for d in diag if d > 1]
        rec["G"] = cyc
        r3 = [r for r in range(k) if diag[r] % 3 == 0]
        if len(r3) != 1:
            rec["status"] = f"3-part of G not cyclic or trivial: {cyc}"
            print("  " + rec["status"])
            break
        r3 = r3[0]
        # G coordinates as e_class returns them: rows with D > 1, in order
        gpos = [r for r in range(k) if diag[r] > 1].index(r3)
        bnf = ES[0]
        Hs = pari("o3_subgroup")(ES, r3 + 1)
        pol = pari.bnrclassfield(bnf, Hs, 1)
        rec["relative_cubic"] = str(pol)
        assert int(pari.poldegree(pol)) == 3, pol
        print(f"  cubic extension of E for ker chi: {pol}  ({time.time() - t0:.0f} s)", flush=True)
        nfE = pari("b -> b.nf")(bnf)
        eq = pari.rnfequation(nfE, pol, 1)
        Qy = pari.subst(eq[0], "x", "y")
        ay = pari.subst(pari.lift(eq[1]), "x", "y")
        red_ = pari.polredbest(Qy, 1)
        Q2 = pari.subst(red_[0], pari.variable(red_[0]), "y")
        bmap = pari.subst(pari.lift(red_[1]), pari.variable(pari.lift(red_[1])), "y")   # root of Qy in terms of Q2
        ay2 = pari("(ay, b, Q2) -> lift(subst(ay, 'y, Mod(b, Q2)))")(ay, bmap, Q2)       # w in terms of Q2
        nfA = pari.nfinit(Q2)
        print(f"  absolute field of degree {int(pari.poldegree(Q2))}: nfinit done ({time.time() - t0:.0f} s)", flush=True)
        coeffs = [pari.lift(pari.nfbasistoalg(bnf, c)) for c in tor.H]
        Habs = pari("(c, ay, Q2) -> sum(k = 1, #c, subst(c[k], 'w, Mod(ay, Q2)) * 'x^(k - 1))")(coeffs, ay2, Q2)
        fa = pari.nffactor(nfA, Habs)
        degs = [int(pari.poldegree(g)) for g in fa[0]]
        rec["factors_over_L"] = degs
        print(f"  H over L: factors of degree {degs}  ({time.time() - t0:.0f} s)", flush=True)
        if sorted(degs) != [N // 3] * 3:
            rec["status"] = "H does not split into three factors of degree N/3 over L"
            print("  " + rec["status"])
            break
        Fa = pari("(ay, Q2) -> Mod(ay, Q2)")(ay2, Q2)
        prs = pari("o3_prime_above")(nfA, Fa, e1, q)
        rec["primes_of_L_above_p"] = len(prs)
        assert len(prs) >= 1
        sides = []
        for g in fa[0]:
            h1 = pari("o3_reduce")(nfA, g, prs[0], q)
            sides.append(sorted(int(pari.lift(x)) for x in pari.polrootsmod(h1, q)))
        allr = sorted(x for s in sides for x in s)
        rec["partition_of_roots"] = allr == sorted(red)
        print(f"  roots of the three factors at one prime of L above (q, F - e1): {sides}; a partition of the roots of "
              f"H: {rec['partition_of_roots']}", flush=True)
        assert rec["partition_of_roots"]
        # the same partition at the other primes above p, up to a permutation of the factors
        same = True
        for pr in prs[1:]:
            s2 = []
            for g in fa[0]:
                h1 = pari("o3_reduce")(nfA, g, pr, q)
                s2.append(frozenset(int(pari.lift(x)) for x in pari.polrootsmod(h1, q)))
            same &= set(s2) == {frozenset(s) for s in sides}
        rec["same_partition_at_other_primes"] = same
        # class field theory against e_class: the cubic splits mod a prime l of E of degree 1 iff its class is in ker chi
        P = ES[3]
        agree, tested, l = collections.Counter(), 0, 2
        while tested < a.check_primes and l < 20000:
            l = int(pari.nextprime(l + 1))
            if l == q:
                continue
            for b in pari.polrootsmod(P, l):
                b = int(pari.lift(b))
                EP = tor.prime(l, b)
                if EP == 0:
                    continue
                cf = [pari.polcoef(pol, j) for j in range(4)]
                try:
                    vals = [int(cm.G("e_red")(EP, pari.lift(c) if pari.type(c) == "t_POLMOD" else c)) for c in cf]
                except Exception:
                    continue
                if min(vals) < 0:
                    continue
                cub = pari.Pol(list(reversed([pari.Mod(v, l) for v in vals])), "X")
                if int(pari.poldegree(cub)) != 3 or int(pari.poldegree(pari.gcd(cub, pari.deriv(cub)))) != 0:
                    continue
                splits = len(pari.polrootsmod(cub, l)) == 3
                g = [int(x) for x in cm.G("e_class")(ES, b, l)]
                inker = g[gpos] % 3 == 0
                agree[(inker, splits)] += 1
                tested += 1
        rec["class_field_check"] = {f"in ker chi {k_[0]}, cubic splits {k_[1]}": v for k_, v in sorted(agree.items())}
        rec["class_field_check_ok"] = all(k_[0] == k_[1] for k_ in agree)
        print(f"  primes of E of degree 1 (l < {l + 1}): {rec['class_field_check']}; consistent: "
              f"{rec['class_field_check_ok']}", flush=True)
        surfaces = []
        for side, s in enumerate(sides):
            for x in s:
                cfgs = enum_set.get(red[x], [])
                surfaces.append({"root": x, "invariants": list(red[x]), "configurations": [list(c) for c in cfgs],
                                 "coset3": side})
        rec["surfaces"] = sorted(surfaces, key=lambda s: s["root"])
        rec["status"] = "ok" if rec["roots_equal_enumerated"] and rec["class_field_check_ok"] else "check failed"
        rec["seconds"] = round(time.time() - t0)
        print(f"  status {rec['status']} ({rec['seconds']} s)", flush=True)
        with open(os.path.join(HERE, F_COSETS), "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        break


def cmd_labelcheck(a):
    last = collections.OrderedDict()
    for r in jl("sixlines_cmmethod.jsonl"):
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r
    cos = {}
    for r in jl(F_COSETS):
        if r.get("status") == "ok":
            cos[(r["q"], r["e1"])] = r
    for (q, e, cmd), r in last.items():
        if cmd != "action" or r["status"] != "ok" or math.prod(r["G"]) % 3:
            continue
        if any(s["label"] is None for s in r["surfaces"]):
            print(f"q = {q}, e1 = {e[0]}: not all labelled")
            continue
        G = r["G"]
        i3 = [i for i, d in enumerate(G) if d % 3 == 0]
        assert len(i3) == 1
        part_lab = collections.defaultdict(set)
        for s in r["surfaces"]:
            part_lab[s["label"][i3[0]] % 3].add(frozenset(map(tuple, s["configurations"])))
        c = cos.get((q, e[0]))
        if c is None:
            print(f"q = {q}, e1 = {e[0]}, G = {G}: no cosets computed")
            continue
        part_cos = collections.defaultdict(set)
        for s in c["surfaces"]:
            part_cos[s["coset3"]].add(frozenset(map(tuple, s["configurations"])))
        same = set(map(frozenset, part_lab.values())) == set(map(frozenset, part_cos.values()))
        print(f"q = {q}, e1 = {e[0]}, G = {G}: labels mod 3 give cosets of sizes "
              f"{sorted(len(v) for v in part_lab.values())}; the class-field cosets are the same partition: {same}")


# ================================================================================================= features

POP_A = [(37, 2), (37, 6), (43, 2), (43, 14), (43, 22)]
POP_B = [(23, 2), (29, 2), (29, 6), (41, 2), (41, 22)]
SPLITS = [(J, tuple(i for i in range(6) if i not in J)) for J in itertools.combinations(range(6), 3) if 0 in J]


def split_perms():
    idx = {frozenset(map(frozenset, s)): k for k, s in enumerate(SPLITS)}
    out = []
    for sg in itertools.permutations(range(6)):
        out.append([idx[frozenset((frozenset(sg[i] for i in J), frozenset(sg[i] for i in Jc)))] for J, Jc in SPLITS])
    return out


SPERMS = None


def cubic_type(v):
    """Normal form of the vector of cubic symbols (0, 1, 2 or 3 for zero) of the ten t_s under S_6 (v -> v o sigma),
    a common shift and the sign."""
    import numpy as np
    global SPERMS
    if SPERMS is None:
        SPERMS = np.array(split_perms())
    v = np.array(v)
    best = None
    for e in (1, 2):
        for c in range(3):
            w = np.where(v == 3, 3, (e * v + c) % 3)
            W = w[SPERMS]                                  # (720, 10): W[k, s] = w[sigma_k(s)]
            m = min(map(tuple, W.tolist()))
            if best is None or m < best:
                best = m
    return "".join(map(str, best))


def tvals(pari, A, p):
    """t_s = D(J) D(J^c) over the ten splits, in the order of SPLITS."""
    D = lambda J: int(pari.lift(pari.matdet(pari.matrix(3, 3, [A[r, c] for r in range(3) for c in J]))))
    return [D(J) * D(Jc) % p for J, Jc in SPLITS]


def esym(xs, p):
    e = [1] + [0] * len(xs)
    for x in xs:
        for k in range(len(xs), 0, -1):
            e[k] = (e[k] + e[k - 1] * x) % p
    return e


def n3_quantities(t, p):
    x = [a * a % p for a in t]
    e = esym(x, p)
    out = {"N3 sigma3": e[3], "N3 sigma6": e[6], "N3 sigma9": e[9]}
    for k in range(1, 11):
        for l in range(k, 11):
            if (k + l) % 3 == 0:
                out[f"N3 sigma{k} sigma{l}"] = e[k] * e[l] % p
    V = 1
    for i in range(10):
        for j in range(i + 1, 10):
            V = V * (x[i] - x[j]) % p
    out["N3 V"] = V
    return out


def cubic_families():
    from sixlines_leakage import monomials, J_WEIGHTS, E_WEIGHTS

    def classes(ws):
        seen, out = set(), []
        for a in monomials(3, ws, lambda w: w % 3 == 0):
            b = tuple((2 * x) % 3 for x in a)
            if b in seen:
                continue
            seen.add(a)
            out.append(a)
        return out
    return ([(f"J3 a={''.join(map(str, a))}", a) for a in classes(J_WEIGHTS)],
            [(f"E3 a={''.join(map(str, a))}", a) for a in classes(E_WEIGHTS)])


def order3_values(pari, A, p, j3, e3, F):
    import sixlines_leakage as lk
    Q = lk.config_quantities(pari, A, p)
    t = tvals(pari, A, p)
    v = {}
    bj = [F.sym(x, 3) for x in Q["J"]]
    be = [F.sym(x, 3) for x in Q["p"] + [Q["C"]]]
    for name, a in j3:
        v[name] = lk.mono_symbol(bj, a, 3)
    for name, a in e3:
        v[name] = lk.mono_symbol(be, a, 3)
    for name, x in n3_quantities(t, p).items():
        v[name] = F.sym(x, 3)
    v["X3 cubic type"] = cubic_type([3 if x == 0 else F.sym(x, 3) for x in t])
    return v


def cubic_filter(pari, cfgs, p, rng, j3, e3, trials=3):
    """Names whose value changes under A -> g A D on one of the configurations."""
    import sixlines_leakage as lk
    F = lk.Fp(p)
    bad = set()
    for cfg in cfgs:
        A = lk.matrix_of(pari, cfg, p)
        v0 = order3_values(pari, A, p, j3, e3, F)
        for _ in range(trials):
            while True:
                g = pari.matrix(3, 3, [rng.randrange(p) for _ in range(9)]) * pari.Mod(1, p)
                if pari.matdet(g) != 0:
                    break
            Dg = pari.matdiagonal([rng.randrange(1, p) for _ in range(6)]) * pari.Mod(1, p)
            v1 = order3_values(pari, g * A * Dg, p, j3, e3, F)
            bad |= {k for k in v0 if v0[k] != v1[k]}
    return bad


def cubic_filter_selftest(pari, cfg, p, rng, j3, e3):
    """The filter must drop quantities of the wrong weight: J2, p1 and sigma1 (weight 2) mod cubes."""
    import sixlines_leakage as lk
    F = lk.Fp(p)
    A = lk.matrix_of(pari, cfg, p)
    drops = collections.Counter()
    for _ in range(6):
        while True:
            g = pari.matrix(3, 3, [rng.randrange(p) for _ in range(9)]) * pari.Mod(1, p)
            if pari.matdet(g) != 0:
                break
        Dg = pari.matdiagonal([rng.randrange(1, p) for _ in range(6)]) * pari.Mod(1, p)
        Q0, Q1 = lk.config_quantities(pari, A, p), lk.config_quantities(pari, g * A * Dg, p)
        drops["J2"] += F.sym(Q0["J"][0], 3) != F.sym(Q1["J"][0], 3)
        drops["p1"] += F.sym(Q0["p"][0], 3) != F.sym(Q1["p"][0], 3)
        drops["J3 (weight 3, kept)"] += F.sym(Q0["J"][1], 3) != F.sym(Q1["J"][1], 3)
    return dict(drops)


def last_cosets(with_43=False):
    """The last record with status ok per torsor.  q = 43 is left out unless with_43: the run of note 29 was made
    without it, and its records were added to the file afterwards (features --with-43)."""
    out = {}
    for r in jl(F_COSETS):
        if r.get("status") == "ok" and (with_43 or r["q"] != 43):
            out[(r["q"], r["e1"])] = r
    return out


def cmd_features(a):
    import sixlines_cmmethod as cmm  # noqa: F401  (cm_tmat, cm_J)
    import sixlines_characters as ch
    import sixlines_leakage as lk
    pari = cmm.pari
    rng = random.Random(103)
    cos = last_cosets(a.with_43)
    j3, e3 = cubic_families()
    print(f"cubic monomial classes: J3 {len(j3)}, E3 {len(e3)}")
    # chi_-1 cosets and delta where note 28 has them
    _, _, _, C = lk.build_population()
    chim1 = {}
    for t in C:
        for s in t["surfaces"]:
            chim1[(t["q"], t["e1"], frozenset(s["configs"]))] = (s["coset"], s["delta"])
    quart = lk.prym_tables()
    ffr = collections.defaultdict(list)
    for x in lk.fourfold_records():
        ffr[(x["p"], x["cfg"])].append(x)
    checks = collections.Counter()
    tors_out = []
    for pop, lst in (("A", POP_A), ("B", POP_B)):
        for (q, e1) in lst:
            if (q, e1) not in cos:
                print(f"  {pop} q={q} e1={e1}: no cosets; left out", flush=True)
                continue
            r = cos[(q, e1)]
            rows = {tuple(x["cfg"]): x for x in jl(f"sixlines_p{q}.jsonl")}
            F = lk.Fp(q)
            jf, ef = lk.j_family(q), lk.e_family(q)
            jf = [x for x in jf if x[1] in (2, 4)]
            some = [tuple(c) for s in r["surfaces"] for c in s["configurations"]][:3]
            bad = lk.invariance_filter(pari, jf + ef, some, q, rng)
            jf = [x for x in jf if x[0] not in bad]
            ef = [x for x in ef if x[0] not in bad]
            cubic = (q - 1) % 3 == 0
            bad3 = cubic_filter(pari, some, q, rng, j3, e3) if cubic else set()
            if cubic:
                checks[("cubic filter self-test (changes in 6 trials)",
                        json.dumps(cubic_filter_selftest(pari, some[0], q, rng, j3, e3)))] += 1
            checks[(f"q={q} e1={e1}: monomials dropped (orders 2, 4; order 3)", (len(bad), len(bad3)))] += 1
            tsurf = []
            for s in r["surfaces"]:
                cfgs = [tuple(c) for c in s["configurations"]]
                assert len(cfgs) == 2, (q, e1, cfgs)
                cvals = []
                for cfg in cfgs:
                    row = rows[cfg]
                    assert row["S"][0] in (e1, -e1)
                    twist = row["S"][0] != e1
                    c = ch.candidates(q, ch.forms(cfg, q), twist)
                    v = {"S twist": int(twist)}
                    for k in ("conic", "j2", "j4", "quartic", "orientation", "branch", "branch_zero"):
                        if k in c:
                            v["S " + k] = c[k]
                    v["S type"] = ch.show(c["type"])
                    if "quartic_type" in c:
                        v["S quartic_type"] = "".join(map(str, c["quartic_type"]))
                    qr = quart.get((q, cfg))
                    if qr:
                        v["S c1_values"] = str(sorted({tuple(x["c"][0]) for x in qr["pairs"]
                                                       if x.get("rational_point") == 0 and "c" in x}))
                        v.update(lk.prym_features(qr))
                    A = lk.matrix_of(pari, cfg, q)
                    Q = lk.config_quantities(pari, A, q)
                    ci = lk.chart_invariants(Q["J"], q, r["chart"])
                    checks[("stored invariants of the root reproduced from the configuration", ci == s["invariants"])] += 1
                    for name, d, aa in jf:
                        v[name] = lk.mono_symbol([F.sym(x, d) for x in Q["J"]], aa, d)
                    for name, d, aa in ef:
                        v[name] = lk.mono_symbol([F.sym(x, d) for x in Q["p"] + [Q["C"]]], aa, d)
                    for d in (2, 4):
                        if (q - 1) % d == 0:
                            v[f"S root symbol d={d}"] = F.sym(s["root"], d)
                    if cubic:
                        v.update({k: x for k, x in order3_values(pari, A, q, j3, e3, F).items() if k not in bad3})
                        v["X3 root symbol d=3"] = F.sym(s["root"], 3)
                    sc = lk.scalars(Q, q)
                    sc["root"] = s["root"]
                    v["_scalars"] = sc
                    cvals.append(v)
                keys = sorted(set().union(*[set(x) for x in cvals]) - {"_scalars"})
                sv = {}
                for k in keys:
                    vals = [json.dumps(x.get(k)) for x in cvals]
                    sv[k] = {"pair": sorted(vals), "agree": len(set(vals)) == 1, "value": cvals[0].get(k)}
                sc = {}
                for k in cvals[0]["_scalars"]:
                    vals = [x["_scalars"][k] for x in cvals]
                    sc[k] = vals[0] if len(set(vals)) == 1 else "differs"
                out = {"configs": [list(c) for c in cfgs], "coset3": s["coset3"], "root": s["root"], "values": sv,
                       "scalars": sc}
                cm1 = chim1.get((q, e1, frozenset(cfgs)))
                if cm1 is not None:
                    out["chi_m1_coset"], out["delta"] = cm1
                recs = [x["F"] for c in cfgs for x in ffr.get((q, c), [])]
                if recs:
                    out["fourfold"] = recs
                tsurf.append(out)
            hc = sum(1 for s in tsurf if "chi_m1_coset" in s)
            tors_out.append({"pop": pop, "q": q, "e1": e1, "G": r["G"], "chart": r["chart"], "n": len(tsurf),
                             "surfaces": tsurf, "chi_m1": hc == len(tsurf),
                             "fourfold": sum(1 for s in tsurf if "fourfold" in s)})
            print(f"  {pop} q={q} e1={e1} G={r['G']} n={len(tsurf)}: candidates per surface "
                  f"{len(tsurf[0]['values'])}; chi_-1 cosets {hc}; surfaces with fourfold records "
                  f"{tors_out[-1]['fourfold']}", flush=True)
    print("CHECKS")
    for k, v in sorted(checks.items(), key=str):
        print(f"  {k}: {v}")
    with open(os.path.join(HERE, a.out), "w") as fh:
        json.dump({"torsors": tors_out}, fh)
    print(f"wrote {a.out}")


# ================================================================================================= audit

GRID = 0.01


def partitions3(n):
    """All partitions of range(n) into three blocks of n/3, once each (canonical labels)."""
    import numpy as np
    m = n // 3
    out = []
    for A in itertools.combinations(range(1, n), m - 1):
        A = (0,) + A
        rest = [i for i in range(n) if i not in A]
        for B in itertools.combinations(rest[1:], m - 1):
            B = (rest[0],) + B
            lab = np.full(n, 2, dtype=np.int64)
            lab[list(A)] = 0
            lab[list(B)] = 1
            out.append(lab)
    return np.array(out)


def canon(lab):
    m, out, nxt = {}, [], 0
    for x in lab:
        if x not in m:
            m[x] = nxt
            nxt += 1
        out.append(m[x])
    return tuple(out)


def partitions2(n):
    import numpy as np
    out = []
    for A in itertools.combinations(range(1, n), n // 2 - 1):
        lab = np.ones(n, dtype=np.int64)
        lab[[0] + list(A)] = 0
        out.append(lab)
    return np.array(out)


class Tor:
    def __init__(self, t, labels=None, k=3):
        import numpy as np
        self.t = t
        self.n = n = t["n"]
        self.k = k
        self.obs = np.array(labels if labels is not None else [s["coset3"] for s in t["surfaces"]], dtype=np.int64)
        self.pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        self.I = np.array([i for i, j in self.pairs])
        self.J = np.array([j for i, j in self.pairs])
        self.enum = n <= 12
        if self.enum:
            self.rows = partitions3(n) if k == 3 else partitions2(n)
            self.index = {canon(r.tolist()): i for i, r in enumerate(self.rows)}
            self.D = (self.rows[:, self.J] - self.rows[:, self.I]) % k
            self.obs_row = self.index[canon(self.obs.tolist())]

    def row_of(self, lab):
        return self.index[canon(list(lab))]


def mi_counts(cnt, N):
    """N * MI from a contingency array (..., A, B) of counts summing to N over the last two axes."""
    import numpy as np
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = cnt.sum(axis=-1, keepdims=True)
        rb = cnt.sum(axis=-2, keepdims=True)
        t = np.where(cnt > 0, cnt * np.log2(cnt * N / (ra * rb)), 0.0)
    return t.sum(axis=(-1, -2))


TABLE_CAP = 5_000_000
_TABLES = {}


def tables(r, m):
    """All V x 3 tables with row sums r and column sums m, with log of the labelled count."""
    import numpy as np
    key = (tuple(r), m)
    if key in _TABLES:
        return _TABLES[key]
    out, lw = [], []
    lf = [math.lgamma(x + 1) for x in range(sum(r) + 1)]
    V = len(r)

    def rec(v, rem, cur, w):
        if len(out) > TABLE_CAP:
            return
        if v == V:
            if rem == (0, 0, 0):
                out.append(list(cur))
                lw.append(w)
            return
        tail = sum(r[v + 1:])
        rv = r[v]
        for x in range(min(rv, rem[0]) + 1):
            for y in range(min(rv - x, rem[1]) + 1):
                z = rv - x - y
                if z > rem[2]:
                    continue
                nrem = (rem[0] - x, rem[1] - y, rem[2] - z)
                if sum(nrem) != tail:
                    continue
                rec(v + 1, nrem, cur + [(x, y, z)], w + lf[rv] - lf[x] - lf[y] - lf[z])
    rec(0, (m, m, m), [], 0.0)
    if len(out) > TABLE_CAP:
        _TABLES[key] = None
        return None
    T = np.array(out, dtype=np.float64).reshape(len(out), V, 3)
    lw = np.array(lw)
    w = np.exp(lw - lw.max())
    _TABLES[key] = (T, w / w.sum())
    return _TABLES[key]


def table_stats(T, phi, ncat, rsum):
    """Pair statistic N * MI(phi(u(X), u(Y)); c(Y) - c(X)) for each table T (K, V, 3)."""
    import numpy as np
    K, V, _ = T.shape
    n = int(rsum.sum())
    N = n * (n - 1)
    O = np.zeros((V * V, ncat))
    O[np.arange(V * V), phi.ravel()] = 1.0
    out = np.empty(K)
    step = max(1, 200_000 // (V * V))
    for s in range(0, K, step):
        Ts = T[s:s + step]
        C = np.empty((Ts.shape[0], ncat, 3))
        for d in range(3):
            M = np.einsum("kaj,kbj->kab", Ts, np.roll(Ts, -d, axis=2))
            if d == 0:
                M[:, np.arange(V), np.arange(V)] -= rsum[None, :]
            C[:, :, d] = M.reshape(Ts.shape[0], V * V) @ O
        out[s:s + step] = mi_counts(C, N)
    return out


def surface_pair_codes(u, z3, variant, T):
    """Pair codes for a surface candidate: u codes over surfaces; z3[code] the value in Z/3 or None (zero)."""
    import numpy as np
    from sixlines_leakage import codes
    if variant == "pair":
        V = int(max(u)) + 1
        return np.array([u[i] * V + u[j] for i, j in T.pairs], dtype=np.int64)
    vals = []
    for i, j in T.pairs:
        a, b = z3[u[i]], z3[u[j]]
        vals.append(("z", a, b) if a is None or b is None else (b - a) % 3)
    return codes(vals)


def phi_matrix(V, z3, variant):
    import numpy as np
    cats, phi = {}, np.zeros((V, V), dtype=np.int64)
    for a in range(V):
        for b in range(V):
            if variant == "pair":
                key = (a, b)
            else:
                key = ("z", z3[a], z3[b]) if z3[a] is None or z3[b] is None else (z3[b] - z3[a]) % 3
            phi[a, b] = cats.setdefault(key, len(cats))
    return phi, len(cats)


def torsor_dist(T, spec, extra_labels=()):
    """(values, weights, observed, [observed under each extra labelling]) for one candidate on one torsor, or None."""
    import numpy as np
    from sixlines_leakage import mi_rows
    kind = spec[0]
    if kind == "surface":
        _, u, z3, variant = spec
        u = np.asarray(u, dtype=np.int64)
        if T.enum:
            f = surface_pair_codes(u, z3, variant, T)
            f = np.unique(f, return_inverse=True)[1]
            vals = mi_rows(T.D, f, T.k)
            w = np.full(len(vals), 1.0 / len(vals))
            return vals, w, vals[T.obs_row], [vals[T.row_of(l)] for l in extra_labels]
        V = int(u.max()) + 1
        r = np.bincount(u, minlength=V)
        tb = tables(r.tolist(), T.n // 3)
        if tb is None:
            return None
        Tt, w = tb
        phi, ncat = phi_matrix(V, z3, variant)
        vals = table_stats(Tt, phi, ncat, r.astype(float))

        def obs_of(lab):
            Ob = np.zeros((1, V, 3))
            for x, c in zip(u, lab):
                Ob[0, x, c] += 1
            return float(table_stats(Ob, phi, ncat, r.astype(float))[0])
        return vals, w, obs_of(T.obs), [obs_of(l) for l in extra_labels]
    if not T.enum:
        return None
    if kind == "pairs":
        f = np.unique(np.asarray(spec[1], dtype=np.int64), return_inverse=True)[1]
        vals = mi_rows(T.D, f, T.k)
    else:                                                       # records
        f = np.unique(np.asarray(spec[1], dtype=np.int64), return_inverse=True)[1]
        G = T.rows[:, np.asarray(spec[2])]
        vals = mi_rows(G, f, T.k)
    w = np.full(len(vals), 1.0 / len(vals))
    return vals, w, vals[T.obs_row], [vals[T.row_of(l)] for l in extra_labels]


def convolve_null(dists):
    """The exact distribution of the sum of the per-torsor statistics on the grid GRID (each value rounded down);
    dists = [(values, weights)]."""
    import numpy as np
    acc = np.array([1.0])
    for vals, w in dists:
        b = np.floor(vals / GRID + 1e-7).astype(np.int64)
        hist = np.bincount(b, weights=w)
        new = np.zeros(len(acc) + len(hist))
        nh, na = np.nonzero(hist)[0], np.nonzero(acc)[0]
        if len(nh) <= len(na):
            for ui in nh:
                new[ui:ui + len(acc)] += hist[ui] * acc
        else:
            for ai in na:
                new[ai:ai + len(hist)] += acc[ai] * hist
        acc = new
    return acc


def tail_from(acc, obs_list):
    """P(sum >= observed sum), conservative: the observed bins are summed and lowered by one per torsor less one."""
    obs_bin = sum(int(math.floor(o / GRID + 1e-7)) for o in obs_list)
    thr = max(obs_bin - len(obs_list) + 1, 0)
    return float(acc[thr:].sum()) if thr < len(acc) else 0.0


def tail_pw(dists):
    """Exact P(sum >= observed sum) on the grid GRID, conservative; dists = [(values, weights, observed)]."""
    return tail_from(convolve_null([(v, w) for v, w, _ in dists]), [o for _, _, o in dists])


_AT = None
_SH = None


def _init(tors, shuffles, k=3, lab_key="coset3"):
    global _AT, _SH
    _AT = [Tor(t, [s[lab_key] for s in t["surfaces"]], k) for t in tors]
    _SH = shuffles


def _score(job):
    name, per = job
    dists, detail, shd = [], {}, [[] for _ in range(len(_SH))]
    for ti, spec in per.items():
        T = _AT[ti]
        res = torsor_dist(T, spec, [sh[ti] for sh in _SH])
        if res is None:
            detail[ti] = {"left out": True}
            continue
        vals, w, obs, ext = res
        dists.append((vals, w, obs))
        for k, o in enumerate(ext):
            shd[k].append((vals, w, o))
        exact = None
        if spec[0] == "surface" and spec[3] == "diff":
            import numpy as np
            f = surface_pair_codes(np.asarray(spec[1], dtype=np.int64), spec[2], "diff", T)
            cnt = np.bincount(np.unique(f, return_inverse=True)[1]).astype(float)
            hN = float(-(cnt * np.log2(cnt / cnt.sum())).sum())        # N * H(f(Y) - f(X))
            exact = bool(hN > 0 and abs(obs - hN) < 1e-7)             # f(Y) - f(X) is a function of c(Y) - c(X)
        detail[ti] = {"stat": float(obs), "p_torsor": float(w[vals >= obs - 1e-9].sum()), "max": float(vals.max()),
                      "exact_order3_character": exact}
    if not dists:
        return {"name": name, "p": 1.0, "stat": 0.0, "torsors": 0, "detail": detail, "shuffled_p": [1.0] * len(_SH)}
    acc = convolve_null([(v, w) for v, w, _ in dists])
    return {"name": name, "p": tail_from(acc, [o for _, _, o in dists]), "stat": sum(o for _, _, o in dists),
            "torsors": len(dists), "detail": detail,
            "shuffled_p": [tail_from(acc, [o for _, _, o in d]) for d in shd]}


def holm(ps):
    from sixlines_leakage import holm as h
    return h(ps)


def candidates(tors, pop, attempted):
    """name -> {ti: spec}.  spec: ('surface', codes, z3 values per code, variant) | ('pairs', codes) |
    ('records', codes, surface index)."""
    from sixlines_leakage import codes, Fp
    cands = collections.defaultdict(dict)
    for ti, t in enumerate(tors):
        surf = t["surfaces"]
        q, n = t["q"], t["n"]
        keys = sorted(set().union(*[set(s["values"]) for s in surf]))
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        for k in keys:
            order3 = k.startswith(("J3 ", "E3 ", "N3 ", "X3 root"))
            attempted.add(k + " | pair")
            if order3:
                attempted.add(k + " | diff")
            vals = [s["values"].get(k) for s in surf]
            if any(v is None for v in vals):
                continue
            pv = [json.dumps(v["pair"]) for v in vals]
            if len(set(pv)) > 1:
                u = codes(pv).tolist()
                cands[k + " | pair"][ti] = ("surface", u, [None] * (max(u) + 1), "pair")
            if order3 and all(v["agree"] for v in vals):
                x = [v["value"] for v in vals]
                if len(set(map(str, x))) > 1:
                    u = codes(x).tolist()
                    z3 = [None] * (max(u) + 1)
                    for ui, xi in zip(u, x):
                        z3[ui] = None if xi == "z" else int(xi)
                    cands[k + " | diff"][ti] = ("surface", u, z3, "diff")
        if n > 12:
            continue
        F = Fp(q)
        for k in surf[0]["scalars"]:
            sv = [s["scalars"][k] for s in surf]
            if any(v == "differs" for v in sv):
                continue
            for d in (2, 3, 4):
                if (q - 1) % d:
                    continue
                for op, fn in (("-", lambda a, b: a - b), ("+", lambda a, b: a + b)):
                    nm = f"B2 sym d={d} of {k}(Y) {op} {k}(X)"
                    attempted.add(nm)
                    vals = [("u" if sv[i] is None or sv[j] is None else F.sym(fn(sv[j], sv[i]), d)) for i, j in pairs]
                    if len(set(map(str, vals))) > 1:
                        cands[nm][ti] = ("pairs", codes(vals).tolist())
    return cands


def f_candidates(tors, attempted):
    from sixlines_leakage import codes
    cands = collections.defaultdict(dict)
    for ti, t in enumerate(tors):
        recs = [(si, r) for si, s in enumerate(t["surfaces"]) for r in s.get("fourfold", [])]
        if not recs:
            continue
        for k in recs[0][1]:
            attempted.add(k)
            sub_ = [(si, r[k]) for si, r in recs if r.get(k) is not None]
            if len({json.dumps(v) for _, v in sub_}) > 1:
                cands[k][ti] = ("records", codes([json.dumps(v) for _, v in sub_]).tolist(), [si for si, _ in sub_])
    return cands


def run_population(pool_ctx, tors, cand, attempted, workers, shuffles, k=3, lab_key="coset3"):
    jobs = sorted(cand.items())
    with pool_ctx.Pool(workers, initializer=_init, initargs=(tors, shuffles, k, lab_key)) as pool:
        res = pool.map(_score, jobs, chunksize=2)
    have = {r["name"] for r in res}
    for nm in sorted(attempted - have):
        res.append({"name": nm, "p": 1.0, "stat": 0.0, "torsors": 0, "detail": {}, "constant": True,
                    "shuffled_p": [1.0] * len(shuffles)})
    ps = [r["p"] for r in res]
    from sixlines_leakage import bh
    for r, h_, b_ in zip(res, holm(ps), bh(ps)):
        r["holm"], r["bh"], r["bonferroni"] = h_, b_, min(1.0, len(res) * r["p"])
    return res


def shuffled_labels(tors, nseeds):
    out = []
    for seed in range(1, nseeds + 1):
        rng = random.Random(seed)
        lab = []
        for t in tors:
            c = [s["coset3"] for s in t["surfaces"]]
            rng.shuffle(c)
            lab.append(c)
        out.append(lab)
    return out


def plant_specs(tors, kind, rng, flip=0.0, nexact=None):
    """A planted candidate on the torsors of a population (true cosets)."""
    per = {}
    ids = list(range(len(tors)))
    chosen = set(ids) if nexact is None else set(rng.sample(ids, nexact))
    for ti, t in enumerate(tors):
        c = [s["coset3"] for s in t["surfaces"]]
        n = t["n"]
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        if kind == "char3":
            if ti in chosen:
                e, s0 = rng.choice((1, 2)), rng.randrange(3)
                f = [(e * x + s0) % 3 if rng.random() >= flip else rng.randrange(3) for x in c]
            else:
                f = [rng.randrange(3) for _ in c]
            if len(set(f)) > 1:
                per[ti] = ("surface", f, [0, 1, 2], "diff")
        elif kind == "indicator":
            c0 = rng.randrange(3)
            f = [int(x == c0) for x in c]
            per[ti] = ("surface", f, [None, None], "pair")
        elif kind == "h_of_d":
            if n > 12:
                continue
            h = [rng.randrange(2) for _ in range(3)]
            if len(set(h)) == 1:
                h[rng.randrange(3)] ^= 1
            per[ti] = ("pairs", [h[(c[j] - c[i]) % 3] for i, j in pairs])
        elif kind == "random":
            f = [rng.randrange(3) for _ in c]
            if len(set(f)) > 1:
                per[ti] = ("surface", f, [0, 1, 2], "diff")
    return per


PLANTS = [("char3", 0.0, None, "an exact order-3 character on every torsor"),
          ("char3", 0.125, None, "the same, each value replaced at random w.p. 1/8"),
          ("char3", 0.25, None, "the same, w.p. 1/4"),
          ("char3", 0.0, 2, "exact on 2 torsors, random elsewhere"),
          ("char3", 0.0, 1, "exact on 1 torsor, random elsewhere"),
          ("indicator", 0.0, None, "f = [c = c0] on every torsor (pair)"),
          ("h_of_d", 0.0, None, "h(c(Y) - c(X)), h random, torsors with n <= 12"),
          ("random", 0.0, None, "a random function to Z/3 (null)")]


def cmd_audit(a):
    from multiprocessing import get_context
    with open(os.path.join(HERE, a.features)) as fh:
        feat = json.load(fh)
    ctx = get_context("fork")
    out = open(os.path.join(HERE, a.out), "w")
    summary = {}
    for pop in ("A", "A37", "A43", "B"):
        tors = [t for t in feat["torsors"] if t["pop"] == pop[0] and (len(pop) == 1 or t["q"] == int(pop[1:]))]
        if not tors:
            continue
        attempted = set()
        cand = candidates(tors, pop, attempted)
        shuffles = shuffled_labels(tors, a.shuffles) if pop == "A" else []
        t0 = time.time()
        res = run_population(ctx, tors, cand, attempted, a.workers, shuffles)
        m = len(res)
        summary[pop] = m
        print(f"\nPOPULATION {pop}: " + ", ".join(f"q={t['q']} e1={t['e1']} G={t['G']} n={t['n']}" for t in tors)
              + f"; m = {m} candidates ({len(cand)} defined somewhere; {time.time() - t0:.0f} s)", flush=True)
        for r in res:
            r["pop"] = pop
            out.write(json.dumps(r) + "\n")
        report(res, m, tors)
        if shuffles:
            print("  shuffled labels (seed: smallest p, candidates with p <= 0.05, p <= 0.01):")
            for k in range(len(shuffles)):
                ps = sorted(r["shuffled_p"][k] for r in res)
                best = min(res, key=lambda r: r["shuffled_p"][k])["name"]
                print(f"    seed {k + 1}: {ps[0]:.3g} ({best}), {sum(p <= 0.05 for p in ps)}, {sum(p <= 0.01 for p in ps)}")
                out.write(json.dumps({"pop": pop, "shuffle_seed": k + 1, "min_p": ps[0], "best": best,
                                      "p<=0.05": sum(p <= 0.05 for p in ps), "p<=0.01": sum(p <= 0.01 for p in ps)}) + "\n")
            ps = sorted(r["p"] for r in res)
            print(f"    true labels: {ps[0]:.3g}, {sum(p <= 0.05 for p in ps)}, {sum(p <= 0.01 for p in ps)}")
        # plants
        rng = random.Random({"A": 2103, "B": 2104, "A37": 2137, "A43": 2143}[pop])
        jobs = []
        for kind, flip, nex, tag in PLANTS:
            for k in range(a.plants):
                per = plant_specs(tors, kind, rng, flip, nex)
                if per:
                    jobs.append((f"{tag} #{k}", per))
        with ctx.Pool(a.workers, initializer=_init, initargs=(tors, [])) as pool:
            pres = pool.map(_score, jobs, chunksize=4)
        print(f"  PLANTS on {pop} ({a.plants} each; found = p <= 0.05/m = {0.05 / m:.2g})")
        grp = collections.defaultdict(list)
        for r in pres:
            grp[r["name"].split(" #")[0]].append(r["p"])
            out.write(json.dumps({"pop": pop, "plant": r["name"], "p": r["p"]}) + "\n")
        for k, ps in grp.items():
            ps.sort()
            print(f"    {k}: found {sum(p <= 0.05 / m for p in ps)} of {len(ps)}, p <= 0.05 {sum(p <= 0.05 for p in ps)}, "
                  f"median p {ps[len(ps) // 2]:.2g}", flush=True)
    # population F: fourfold records against the cosets of order 3
    tors = [t for t in feat["torsors"] if t["fourfold"] and t["n"] <= 12]
    attempted = set()
    cand = f_candidates(tors, attempted)
    res = run_population(ctx, tors, cand, attempted, a.workers, [])
    print(f"\nPOPULATION F: " + ", ".join(f"q={t['q']} e1={t['e1']} records "
                                         f"{sum(len(s.get('fourfold', [])) for s in t['surfaces'])}" for t in tors)
          + f"; m = {len(res)}")
    for r in res:
        r["pop"] = "F"
        out.write(json.dumps(r) + "\n")
    report(res, len(res), tors, top=13)
    # positive control: delta against the cosets of ker chi_-1, the same code with two cosets
    ctl = [t for t in feat["torsors"] if t["chi_m1"]]
    from sixlines_leakage import codes
    cand = {}
    for ti, t in enumerate(ctl):
        u = codes([json.dumps(s["delta"]) for s in t["surfaces"]]).tolist()
        cand[ti] = ("records", u, list(range(t["n"])))          # n * MI(delta(X); coset(X)) over surfaces
    with ctx.Pool(1, initializer=_init, initargs=(ctl, [], 2, "chi_m1_coset")) as pool:
        r = pool.map(_score, [("delta against chi_-1", cand)])[0]
    print(f"\nPOSITIVE CONTROL delta against the cosets of ker chi_-1 on "
          + ", ".join(f"q={t['q']} e1={t['e1']}" for t in ctl)
          + f": p = {r['p']:.3g} (per torsor " + ", ".join(f"{d['p_torsor']:.3g}" for d in r["detail"].values())
          + "); thresholds 0.05/m: " + ", ".join(f"{k} {0.05 / v:.2g}" for k, v in summary.items()))
    out.write(json.dumps({"pop": "control", **r}) + "\n")
    out.close()


def family_of(name):
    if name.startswith("B2"):
        return "B2 sym d=" + name.split("d=")[1].split()[0]
    head = name.split()[0]
    if head in ("J", "E"):
        return head + " d=" + name.split("d=")[1].split()[0]
    return head


def report(res, m, tors, top=12):
    fam = collections.defaultdict(list)
    for r in res:
        v = r["name"].split(" | ")[-1] if " | " in r["name"] else ""
        fam[(family_of(r["name"]), v)].append(r)
    print("  family | variant | candidates (constant or undefined) | best | stat | exact p | Holm")
    for (f, v), rs in sorted(fam.items()):
        b = min(rs, key=lambda r: r["p"])
        nc = sum(1 for r in rs if r.get("constant"))
        print(f"  {f} | {v} | {len(rs)} ({nc}) | {b['name']} | {b['stat']:.2f} | {b['p']:.3g} | {b['holm']:.3g}")
    print(f"  the {top} smallest p:")
    for r in sorted(res, key=lambda r: r["p"])[:top]:
        det = "; ".join(f"q={tors[int(ti)]['q']}/{tors[int(ti)]['e1']}: "
                        + ("left out" if d.get("left out") else f"p {d['p_torsor']:.3g}"
                           + (" exact" if d.get("exact_order3_character") else ""))
                        for ti, d in r["detail"].items())
        print(f"    {r['name']}: stat {r['stat']:.2f}, p {r['p']:.3g}, Holm {r['holm']:.3g}, BH {r['bh']:.3g} [{det}]")
    ex = collections.Counter()
    for r in res:
        for ti, d in r["detail"].items():
            if d.get("exact_order3_character"):
                ex[(tors[int(ti)]["q"], tors[int(ti)]["e1"])] += 1
    if ex:
        print(f"  (candidate, torsor) cells that are exact order-3 characters: {dict(ex)}")
    ps = [r["p"] for r in res if not r.get("constant")]
    print(f"  calibration: {len(ps)} defined, p <= 0.05: {sum(p <= 0.05 for p in ps)}, p <= 0.01: {sum(p <= 0.01 for p in ps)}")


def cmd_selfcheck(a):
    """The table method against the enumeration of partitions on torsors of 12 surfaces."""
    import numpy as np
    with open(os.path.join(HERE, a.features)) as fh:
        feat = json.load(fh)
    rng = random.Random(5)
    for t in feat["torsors"]:
        if t["n"] != 12:
            continue
        T = Tor(t)
        attempted = set()
        cand = candidates([t], t["pop"], attempted)
        names = [k for k in sorted(cand) if 0 in cand[k] and cand[k][0][0] == "surface"]
        worst = 0.0
        for k in rng.sample(names, min(25, len(names))):
            spec = cand[k][0]
            v1, w1, o1, _ = torsor_dist(T, spec)
            T.enum = False
            v2, w2, o2, _ = torsor_dist(T, spec)
            T.enum = True
            p1, p2 = float(w1[v1 >= o1 - 1e-9].sum()), float(w2[v2 >= o2 - 1e-9].sum())
            m1, m2 = float((v1 * w1).sum()), float((v2 * w2).sum())
            worst = max(worst, abs(p1 - p2), abs(o1 - o2), abs(m1 - m2))
        print(f"q={t['q']} e1={t['e1']}: 25 candidates, largest difference between the two methods (p, observed, "
              f"mean): {worst:.2e}")


def cmd_fields(a):
    """zeta_3 is not in E: x^2 + x + 1 over E, and the quadratic subfields of E, for the ten polynomials."""
    import cypari2
    pari = cypari2.Pari()
    for (q, e1) in POP_A + POP_B:
        r = [x for x in jl("sixlines_groupsize.jsonl") if x["q"] == q and x["e"][0] == e1 and x.get("maximal")][0]
        e1_, e2, e3 = r["e"]
        P = pari(f"w^6 - ({e1_})*w^5 + ({e2})*w^4 - ({e3})*w^3 + ({q * q * e2})*w^2 - ({q ** 4 * e1_})*w + {q ** 6}")
        nf = pari.nfinit(P)
        nfac = len(pari.nffactor(nf, pari("x^2 + x + 1"))[0])
        print(f"q = {q}, e1 = {e1}, G = {r['G_cyc']}: factors of x^2 + x + 1 over E {nfac}; quadratic subfields of E "
              f"{[str(x[0]) for x in pari.nfsubfields(P, 2)]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("cosets")
    for name in ("q", "e1", "mu", "nu"):
        s.add_argument(name, type=int)
    s.add_argument("--digits", type=int, nargs="+", default=[1200])
    s.add_argument("--workers", type=int, default=8)
    s.add_argument("--check-primes", type=int, default=60)
    s.add_argument("--mix", type=int, nargs="+", default=None,
                   help="a coordinate (integer combination of the invariants of the chart) to try first")
    s.add_argument("--identify-by-coordinate", action="store_true",
                   help="for the --mix coordinate, identify surfaces by its value if the invariant polynomials fail")
    sub.add_parser("labelcheck")
    s = sub.add_parser("features")
    s.add_argument("--with-43", action="store_true", help="include the q = 43 records (added after the run of note 29)")
    s.add_argument("--out", default=F_FEATURES)
    sub.add_parser("fields")
    s = sub.add_parser("selfcheck")
    s.add_argument("--features", default=F_FEATURES)
    s = sub.add_parser("audit")
    s.add_argument("--workers", type=int, default=12)
    s.add_argument("--shuffles", type=int, default=10)
    s.add_argument("--plants", type=int, default=200)
    s.add_argument("--out", default=F_OUT)
    s.add_argument("--features", default=F_FEATURES)
    a = ap.parse_args()
    {"cosets": cmd_cosets, "labelcheck": cmd_labelcheck, "features": cmd_features, "audit": cmd_audit,
     "selfcheck": cmd_selfcheck, "fields": cmd_fields}[a.cmd](a)


if __name__ == "__main__":
    main()
