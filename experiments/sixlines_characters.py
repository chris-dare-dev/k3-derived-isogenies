#!/usr/bin/env python
"""
sixlines_characters.py -- quadratic characters of the group acting on six-line K3 surfaces, and candidates for
computing one from a surface (note 19).

G = Cl(E) / <Cl(E0), primes of E ramified over E0> acts simply transitively on the six-line surfaces with a given
Frobenius polynomial P when the order is maximal (note 13, Lemma 9).  A function phi from these surfaces to +-1
"is a character" if phi(g X) = chi(g) phi(X) for a character chi of order at most 2 of G.  Then phi is constant on
the torsor (chi trivial) or takes each value on exactly half of it.

  genus     the class-group side, on polynomials only (also: whether d = -1, 2, -2 give characters of G, and whether
            the classes of the ideals Psi(P) of note 15 lie in the kernel of every character of order 2).  For every polynomial with maximal order in
            sixlines_groupsize.jsonl (note 13): the characters of order 2 of G; each is identified with the
            character of E(sqrt d) / E for an S-unit d of the real cubic field E0, by its values on primes of
            degree 1; the signs of d at the three real places; whether d can be a unit of E0, or a unit away
            from 2.  Also the class of the unit-root prime above p, which is Frobenius conjugation (note 3).

  surfaces  the surface side, on the enumerated surfaces of the 76 irreducible ordinary polynomials of notes 2
            and 3 (prime fields 13..37, F_25, F_49).  For every configuration:
              * the ten square classes eps [ijk][lmn], one for each split of the six lines in two triples, of the
                surface eps w^2 = l1 ... l6, reduced to a normal form under relabelling ("type");
              * for q = 1 mod 4, the quartic characters of the ten eps [ijk][lmn], defined up to a common sign:
                their product ("quartic"), their normal form under relabelling ("quartic type"), and a sign that
                changes with the choice of the quartic character ("orientation");
              * the square classes of three polynomial invariants: the conic invariant
                [123][145][246][356] - [124][135][236][456], and the power sums of the ten products
                t = [ijk][lmn] of degree 2 and 4 ("conic", "j2", "j4");
              * the number of pairs of lines for which the two branch values of the cross-ratio map of note 14
                are rational ("branch");
            and from the stored quartics of notes 14 and 16 the values of c1 over the pairs of one surface.
            Each candidate is tested for: the same value on a configuration and its associate; the same value
            after a relabelling of the lines; constant or balanced on each torsor of a maximal order; invariance
            under Frobenius conjugation over F_25 and F_49.

Class groups are PARI's, without certification (GRH).

  python -I sixlines_characters.py genus [--workers n] [--min-q Q] [--max-q Q] [--out FILE] | --summary FILE [FILE ...]
  python -I sixlines_characters.py surfaces [--out FILE]
  python -I sixlines_characters.py dyadic FILE          FILE = the output of genus: compares "d = -1 gives a character
                                                        of G" with the primes above 2 of the real quartic field K+

Run from experiments/.
"""
import argparse
import itertools
import json
import os
import random
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402
from sixlines_fq import FieldExt  # noqa: E402
import sixlines_prym as sp  # noqa: E402  (pr_* helpers are not used; the module sets up PARI variables)
import sixlines_prym_fq as spq  # noqa: E402  (pq_lambda, pq_curve)

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ class groups

pari(r"""
ch_group(bnf, bnf0, q) = {
  \\ [H, s]: H = HNF of the relations of G = Cl(E) / <Cl(E0), ramified primes> on the generators of Cl(E)
  \\ (a k x k matrix, k = #bnf.cyc), s = F + q^2 / F as an element of E
  my(nf = bnf.nf, P = nf.pol, t = variable(P), F = Mod(t, P), s = F + q^2 / F, cyc = bnf.cyc, k = #cyc,
     ram = abs(nf.disc) / bnf0.disc^2, M,
     up = (I -> my(two = idealtwoelt(bnf0, I)); idealhnf(nf, two[1], subst(lift(nfbasistoalg(bnf0, two[2])), 'x, s))));
  if(k == 0, return([matrix(0, 0), s]));
  M = matdiagonal(cyc);
  for(j = 1, #bnf0.gen, M = concat(M, Mat(bnfisprincipal(bnf, up(bnf0.gen[j]), 0))));
  foreach(factor(ram)[, 1]~, l,
    foreach(idealprimedec(bnf0, l), pr,
      my(fa = idealfactor(nf, up(pr)));
      if(matsize(fa)[1] == 1 && fa[1, 2] == 2, M = concat(M, Mat(bnfisprincipal(bnf, fa[1, 1], 0))))));
  [mathnf(M), s];
}
""")

pari(r"""
ch_characters(H) = {
  \\ a basis of the characters of order 2 of Z^k / H: the rows v mod 2 with v H = 0 mod 2
  my(k = matsize(H)[1]);
  if(k == 0, return([]));
  apply(v -> lift(v~), Vec(matker(Mod(H~, 2))));
}
""")

pari(r"""
ch_inlattice(H, v) = {
  \\ is the class with exponent vector v trivial in Z^k / H ?
  if(matsize(H)[1] == 0, return(1));
  my(w = matsolve(H, v)); denominator(w) == 1;
}
""")

pari(r"""
ch_unitroot(bnf, q) = {
  \\ the class of the prime of E above p at which F is a unit and F-bar is not (the unit-root prime)
  my(nf = bnf.nf, p = factor(q)[1, 1], F = Mod(variable(nf.pol), nf.pol), out = 0);
  foreach(idealprimedec(nf, p), pr, if(nfeltval(nf, F, pr) == 0, out = pr));
  if(out == 0, error("no unit-root prime"));
  if(out.f != 1 || out.e != 1, error("the unit-root prime is not of degree 1"));
  bnfisprincipal(bnf, out, 0);
}
""")

pari(r"""
ch_sunits(bnf0, q) = {
  \\ generators of a group of S-units of E0 that contains every d with (d) = (square) (power of dyadic primes):
  \\ S = the primes above 2 and the primes dividing the generators of Cl(E0).
  \\ [list of elements as polynomials in x, number of them that are units, number that are units or dyadic]
  my(dy = idealprimedec(bnf0, 2), S = dy, T = dy, units, out);
  for(j = 1, #bnf0.gen, foreach(idealfactor(bnf0, bnf0.gen[j])[, 1]~, pr, if(!setsearch(Set(S), pr), S = concat(S, [pr]))));
  units = concat([bnf0.tu[2]], bnf0.fu);
  out = concat(apply(lift, units), apply(lift, bnfsunit(bnf0, dy)[1]));
  my(n1 = #units, n2 = #out);
  if(#S > #dy, out = concat(out, apply(lift, bnfsunit(bnf0, S)[1])));
  [out, n1, n2];
}
""")

pari(r"""
ch_primes(bnf, H, s, q, ds, count) = {
  \\ for `count` primes L of E of degree 1 above odd l not dividing the discriminant, q, or the norms of the d:
  \\ [exponent vector of the class of L, vector of (1 - Legendre(d(x_L) / l)) / 2 for the d in ds]
  \\ with x_L the image of s = F + q^2 / F in the residue field
  my(nf = bnf.nf, P = nf.pol, bad = 2 * q * nf.disc, out = List(), l = 2, rts, x, leg, ok, sl = lift(s));
  while(#out < count,
    l = nextprime(l + 1);
    if(bad % l == 0 || denominator(content(sl)) % l == 0, next);
    rts = polrootsmod(P, l);
    foreach(rts, r,
      x = subst(sl, variable(P), r);
      ok = 1; leg = vector(#ds, j, my(v = subst(ds[j] * Mod(1, l), 'x, x)); if(v == 0, ok = 0); if(v != 0 && issquare(v), 0, 1));
      if(!ok, next);
      listput(out, [bnfisprincipal(bnf, idealprimedec(nf, l)[select(pr -> nfeltval(nf, Mod(variable(P), P) - lift(r), pr) > 0, idealprimedec(nf, l), 1)[1]], 0), leg, l])));
  Vec(out);
}
""")

pari(r"""
ch_solve(A, b) = {
  \\ a solution c over F_2 of A c = b, or 0
  my(c = matinverseimage(Mod(A, 2), Mod(b, 2)));
  if(#c == 0, 0, lift(c));
}
""")

pari(r"""
ch_signs(bnf0, d) = {
  \\ number of real places of E0 at which d is negative
  my(r = polrootsreal(bnf0.pol));
  sum(j = 1, #r, subst(d, 'x, r[j]) < 0);
}
""")

pari(r"""
ch_psi(bnf, q, e1, count) = {
  \\ maximal orders (mu = nu = 1): chi = X^4 - 4 (1 + i) X^3 + e1 X^2 - 4 q (1 - i) X + q^2.  For `count` pairs
  \\ (l, root a of chi mod (l, i - i_l)), l = 1 mod 4: the exponent vector of the class of the ideal
  \\ Psi(P) = (l, a^3 (chi / (X - a))(t / a)) of E, of norm l^3 (note 15, section 3.1)
  my(nf = bnf.nf, P = nf.pol, t = variable(P), out = List(), l = 3, il, chi, rts, ca, Pa, I);
  while(#out < count,
    l = nextprime(l + 1);
    if(l % 4 != 1 || (2 * q * nf.disc) % l == 0, next);
    il = sqrt(Mod(-1, l));
    foreach([il, -il], j,
      chi = 'X^4 - 4 * (1 + j) * 'X^3 + e1 * 'X^2 - 4 * q * (1 - j) * 'X + q^2;
      rts = polrootsmod(lift(chi), l);
      foreach(rts, a,
        ca = (chi * Mod(1, l)) / ('X - a);
        Pa = lift(a^3 * subst(ca, 'X, t / a));
        I = idealadd(nf, l, Mod(Pa, P));
        if(idealnorm(nf, I) != l^3, error("Psi(P) is not of norm l^3"));
        listput(out, bnfisprincipal(bnf, I, 0)))));
  Vec(out);
}
""")

pari(r"""
ch_dyadic(q, e1) = {
  \\ maximal orders: the real quartic field K+ = Q(pi + q / pi) of the fourfold, pi a root of
  \\ chi = X^4 - 4 (1 + i) X^3 + e1 X^2 - 4 q (1 - i) X + q^2.  [e, f] for its primes above 2; 0 if chi conj(chi) is reducible
  my(chi = 'X^4 - 4 * (1 + I) * 'X^3 + e1 * 'X^2 - 4 * q * (1 - I) * 'X + q^2, W = real(chi * conj(chi)), m);
  if(!polisirreducible(W), return(0));
  m = minpoly(Mod('X, W) + q / Mod('X, W), 'y);
  apply(pr -> [pr.e, pr.f], idealprimedec(nfinit([m, [2]]), 2));
}
""")

GPC = {n: pari(n) for n in ["ch_dyadic", "ch_psi", "ch_group", "ch_characters", "ch_inlattice", "ch_unitroot", "ch_sunits", "ch_primes",
                            "ch_solve", "ch_signs"]}


def fields_of(q, e):
    """bnf of E and of E0 for the polynomial (q, e); PARI sometimes needs a higher precision."""
    P = cm.pol(cm.charpoly(q, e))
    e1, e2, e3 = e
    h = pari(f"x^3 - ({e1})*x^2 + ({e2 - 3 * q * q})*x - ({e3 - 2 * q * q * e1})")
    for prec in (38, 100, 250):
        try:
            pari(f"default(realprecision, {prec})")
            return pari.bnfinit(pari.nfinit(P), 1), pari.bnfinit(h, 1)
        except Exception:
            continue
        finally:
            pari("default(realprecision, 38)")
    raise RuntimeError(f"bnfinit failed for q = {q}, e = {e}")


def group_record(q, e, identify=True, nprimes=60):
    """The characters of order 2 of G for one polynomial with maximal order."""
    bnf, bnf0 = fields_of(q, e)
    H, s = GPC["ch_group"](bnf, bnf0, q)
    k = int(pari.matsize(H)[0])
    cyc = [int(d) for d in pari.matsnf(H) if int(d) > 1] if k else []
    chars = [[int(x) for x in v] for v in GPC["ch_characters"](H)] if k else []
    r2 = len(chars)
    assert r2 == sum(d % 2 == 0 for d in cyc)
    rec = {"q": q, "e": list(e), "G": 1, "G_cyc": cyc, "rank2": r2, "h0": int(pari("(b) -> b.no")(bnf0))}
    for d in cyc:
        rec["G"] *= d
    if k:
        up = [int(x) for x in GPC["ch_unitroot"](bnf, q)]
        rec["frobenius_trivial"] = bool(GPC["ch_inlattice"](H, pari.Col(up)))
        rec["frobenius_on_characters"] = [sum(a * b for a, b in zip(v, up)) % 2 for v in chars]
    else:
        rec["frobenius_trivial"], rec["frobenius_on_characters"] = True, []
    if not identify or r2 == 0:
        rec["characters"] = []
        return rec
    ds, n_units, n_dyadic = GPC["ch_sunits"](bnf0, q)
    n_units, n_dyadic = int(n_units), int(n_dyadic)
    data = GPC["ch_primes"](bnf, H, s, q, ds, nprimes + len(ds))
    cls = [[int(x) for x in row[0]] for row in data]
    A = pari.matrix(len(data), len(ds), [int(x) for row in data for x in row[1]])
    # the characters of the rational d = -1, 2, -2: "the degree of the step mod 4, mod 8"
    ells = [int(x) for x in pari("(v) -> apply(r -> r[3], v)")(data)]
    B = pari.matrix(len(data), r2, [sum(a * c for a, c in zip(v, cl)) % 2 for cl in cls for v in chars])
    rec["rational"] = {}
    for name, d in (("-1", -1), ("2", 2), ("-2", -2)):
        leg = pari.Col([0 if int(pari.kronecker(d, l)) == 1 else 1 for l in ells])
        rec["rational"][name] = GPC["ch_solve"](B, leg) != 0
    # note 15: the images of primes of K
    psi = GPC["ch_psi"](bnf, q, e[0], 40)
    rec["psi_tested"] = len(psi)
    rec["psi_in_kernel_of_all_characters"] = all(
        sum(a * int(c) for a, c in zip(v, cl)) % 2 == 0 for cl in psi for v in chars)
    rec["characters"] = []
    for v in chars:
        b = pari.Col([sum(a * c for a, c in zip(v, cl)) % 2 for cl in cls])
        level = None
        for name, n in (("unit", n_units), ("dyadic", n_dyadic), ("class group", len(ds))):
            if n == 0:
                continue
            sub = pari.matrix(len(data), n, [int(x) for row in data for x in list(row[1])[:n]])
            c = GPC["ch_solve"](sub, b)
            if c != 0:
                level = name
                break
        assert level is not None, f"a character of G is not of the form E(sqrt d), d in E0: q = {q}, e = {e}"
        d = pari(1)
        for cj, dj in zip(c, ds):
            if int(cj):
                d = d * dj
        d = pari("(d, b) -> lift(Mod(d, b.pol))")(d, bnf0)
        neg = int(GPC["ch_signs"](bnf0, d))
        rec["characters"].append({"level": level, "negative_places": min(neg, 3 - neg)})
    return rec


def _genus_job(job):
    q, e = job
    try:
        return group_record(q, e)
    except Exception as exc:                                       # recorded, not hidden
        return {"q": q, "e": list(e), "failed": str(exc)[:200]}


def cmd_genus(args):
    if args.summary:
        genus_summary([json.loads(line) for name in args.summary for line in open(name)])
        return
    rows = [json.loads(line) for line in open(os.path.join(HERE, "sixlines_groupsize.jsonl"))]
    jobs = []
    for r in rows:
        if not (r.get("maximal") and "G" in r) or r["source"] == "real":
            continue
        p = int(pari.factor(r["q"])[0][0])
        if r["q"] not in (p, p * p) or not args.min_q <= r["q"] <= args.max_q:   # note 13 leaves out q = p^f, f >= 3
            continue
        jobs.append((r["q"], tuple(r["e"])))
    known = {(r["q"], tuple(r["e"])): r for r in rows if "G" in r}
    fh = open(args.out, "w") if args.out else None
    recs = []
    with Pool(args.workers) as pool:
        for rec in pool.imap(_genus_job, jobs):                     # written as they come, in the order of the input
            recs.append(rec)
            assert "failed" in rec or rec["G"] == known[(rec["q"], tuple(rec["e"]))]["G"], "G differs from note 13"
            if fh:
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
    genus_summary(recs)


def genus_summary(recs):
    st = Counter()
    rank, level, places, frob, rat = Counter(), Counter(), Counter(), Counter(), Counter()
    bound = Counter()
    for rec in recs:
        st["polynomials"] += 1
        if "failed" in rec:
            st["failed"] += 1
            print("  failed:", rec, flush=True)
            continue
        rank[rec["rank2"]] += 1
        st["characters"] += rec["rank2"]
        for c in rec["characters"]:
            level[c["level"]] += 1
            places[c["negative_places"]] += 1
        if rec["rank2"]:
            st["polynomials with a character"] += 1
            st["  the image of Cl(K) (note 15) is in the kernel of all of them"] += rec["psi_in_kernel_of_all_characters"]
            st["  primes of K tested"] += rec["psi_tested"]
            for name, yes in rec["rational"].items():
                rat[(name, rec["q"] % 8)] += yes
                rat[("polynomials", rec["q"] % 8)] += name == "-1"
        bound["rank2 > 3"] += rec["rank2"] > 3
        bound["rank2 > 2 and h(E0) odd"] += rec["rank2"] > 2 and rec["h0"] % 2 == 1
        bound["h(E0) odd"] += rec["h0"] % 2 == 1
        p = int(pari.factor(rec["q"])[0][0])
        if rec["q"] != p:
            frob["torsors over F_{p^2}"] += 1
            frob["Frobenius conjugation trivial in G"] += rec["frobenius_trivial"]
            frob["Frobenius conjugation not trivial, a square (no character sees it)"] += (
                not rec["frobenius_trivial"] and not any(rec["frobenius_on_characters"]))
            frob["some character is -1 on Frobenius conjugation"] += any(rec["frobenius_on_characters"])
        else:
            assert rec["frobenius_trivial"], "the unit-root prime is not trivial in G over a prime field"
    print("genus:", dict(st))
    print("  2-rank of G:", dict(sorted(rank.items())))
    print("  every character of order 2 is E(sqrt d)/E with d in E0; the smallest group that contains d:", dict(level))
    print("  real places of E0 where d is negative (or positive, whichever is fewer):", dict(sorted(places.items())))
    print("  ", dict(bound))
    print("  ", dict(frob))
    print("  the character of E(sqrt d)/E for rational d is a character of G (by q mod 8; polynomials with a character):")
    for r in (1, 3, 5, 7):
        print(f"    q = {r} mod 8: {rat[('polynomials', r)]} polynomials; d = -1: {rat[('-1', r)]}, d = 2: {rat[('2', r)]}, d = -2: {rat[('-2', r)]}")
    sys.stdout.flush()


# ------------------------------------------------------------------------------------------------ surfaces

pari(r"""
ch_setup(q, modulus) = {
  \\ [generator w of F_q with the modulus of sixlines_fq.py, the non-square n0 with the least code, a fixed
  \\  element of order 4 if q = 1 mod 4]
  my(p = factor(q)[1, 1], w = if(q == p, ffgen(p, 'w), ffgen(Mod(1, p) * modulus, 'w)), one = w^0, el, n0 = 0);
  for(c = 1, q - 1, el = if(q == p, c * one, (c % p) * one + (c \ p) * w); if(!issquare(el), n0 = el; break));
  [w, n0, if(q % 4 == 1, ffprimroot(w)^((q - 1) / 4), 0)];
}
""")

pari(r"""
ch_branch(N, D) = {
  \\ the discriminant of the quadratic N' D - N D', whose roots are the two lines of the pencil at which the
  \\ cross-ratio map N / D ramifies
  poldisc(deriv(N, 't) * D - N * deriv(D, 't), 't);
}
""")

_FIELD = {}
TRIPLES = list(itertools.combinations(range(6), 3))
SPLITS = [(t, tuple(i for i in range(6) if i not in t)) for t in TRIPLES if 0 in t]
SPLIT_INDEX = {s[0]: i for i, s in enumerate(SPLITS)}
PERMS = list(itertools.permutations(range(6)))


def field(q):
    if q not in _FIELD:
        p = int(pari.factor(q)[0][0])
        if q == p:
            w, n0, i4 = pari("ch_setup")(q, 0)
        else:
            mod = FieldExt(p, 2)._irreducible()
            w, n0, i4 = pari("ch_setup")(q, pari(f"x^2 + {mod[1]}*x + {mod[0]}"))
        _FIELD[q] = (p, w, n0, i4)
    return _FIELD[q]


def element(code, q):
    p, w, n0, i4 = field(q)
    one = w ** 0
    return code * one if q == p else (code % p) * one + (code // p) * w


def forms(cfg, q):
    b, c, e, f = (element(v, q) for v in cfg)
    o, z = element(1, q), element(0, q)
    return [(o, z, z), (z, o, z), (z, z, o), (o, o, o), (o, b, c), (o, e, f)]


def eta(x):
    return 1 if pari.issquare(x) else -1


def det3(u, v, w):
    return u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])


def _sorted_sign(t):
    s = 1
    for i in range(3):
        for j in range(i + 1, 3):
            if t[i] > t[j]:
                s = -s
    return tuple(sorted(t)), s


def _split_of(t):
    return SPLIT_INDEX[t] if 0 in t else SPLIT_INDEX[tuple(i for i in range(6) if i not in t)]


def _relabel_table():
    """For each relabelling l'_i = l_{pm[i]}: the split and the sign with t'_j = sign * t_{split}."""
    table = []
    for pm in PERMS:
        row = []
        for A, B in SPLITS:
            a, sa = _sorted_sign([pm[i] for i in A])
            b, sb = _sorted_sign([pm[i] for i in B])
            row.append((_split_of(a), sa * sb))
        table.append(row)
    return table


RELABEL = _relabel_table()


def normal_form(beta, minus_one):
    """The largest image of the vector of ten square classes under the 720 relabellings; minus_one = eta(-1)."""
    best = None
    for row in RELABEL:
        v = tuple(beta[j] * (minus_one if s < 0 else 1) for j, s in row)
        if best is None or v > best:
            best = v
    return best


def quartic_form(v, minus_one, conj=True):
    """Normal form of the exponents (mod 4) of the quartic characters of the ten products: under relabelling, under
    v -> v + 2 (the products are defined up to a common square), and, if conj, under v -> -v (the choice of the
    quartic character).  minus_one = the exponent of chi_4(-1)."""
    best = None
    for row in RELABEL:
        u = tuple((v[j] + (minus_one if s < 0 else 0)) % 4 for j, s in row)
        for u2 in (u, tuple((x + 2) % 4 for x in u)):
            for u3 in ((u2, tuple((-x) % 4 for x in u2)) if conj else (u2,)):
                if best is None or u3 < best:
                    best = u3
    return best


def split_products(L):
    br = {t: det3(L[t[0]], L[t[1]], L[t[2]]) for t in TRIPLES}
    return br, [br[A] * br[B] for A, B in SPLITS]


def candidates(q, L, twist=False, with_pairs=True):
    """Candidates of a labelled configuration given by six linear forms L, for the surface eps w^2 = l1 ... l6
    with eps = 1, or eps = the fixed non-square if twist."""
    p, w, n0, i4 = field(q)
    one = w ** 0
    eps = n0 if twist else one
    br, t = split_products(L)
    beta = tuple(eta(eps * x) for x in t)
    out = {"beta": beta, "type": normal_form(beta, eta(-one))}
    conic = br[(0, 1, 2)] * br[(0, 3, 4)] * br[(1, 3, 5)] * br[(2, 4, 5)] - br[(0, 1, 3)] * br[(0, 2, 4)] * br[(1, 2, 5)] * br[(3, 4, 5)]
    j2 = sum((x * x for x in t), 0 * one)
    j4 = sum((x ** 4 for x in t), 0 * one)
    out["conic"] = 0 if conic == 0 else eta(conic)
    out["j2"] = 0 if j2 == 0 else eta(j2)
    out["j4"] = 0 if j4 == 0 else eta(j4)
    if q % 4 == 1:
        ex = (q - 1) // 4
        tab = {str(one): 0, str(i4): 1, str(-one): 2, str(-i4): 3}
        v = tuple(tab[str((eps * x) ** ex)] for x in t)
        m1 = tab[str((-one) ** ex)]
        assert all((x % 2 == 0) == (b == 1) for x, b in zip(v, beta))
        # the product of the ten quartic characters does not change under v -> v + 2 or under relabelling
        assert sum(v) % 2 == 0, "the product of the twenty brackets is not a square"
        out["quartic"] = 1 if sum(v) % 4 == 0 else -1
        out["quartic_type"] = quartic_form(v, m1)
        out["quartic_oriented"] = quartic_form(v, m1, conj=False)
        a = sum(1 if x == 0 else -1 for x in v if x % 2 == 0)
        b = sum(1 if x == 1 else -1 for x in v if x % 2 == 1)
        out["orientation"] = (a * b > 0) - (a * b < 0)
    if not with_pairs:
        return out
    Lm = pari.matrix(6, 3, [x for row in L for x in row])
    branch, c4 = [], set()
    for pair in itertools.combinations(range(1, 7), 2):
        rest = [i for i in range(1, 7) if i not in pair]
        N, D = spq.GP["pq_lambda"](Lm, rest + list(pair))
        d = pari("ch_branch")(N, D)
        branch.append(0 if d == 0 else eta(d))
        cur = spq.GP["pq_curve"](N, D, n0)
        if int(cur[0]) == 0:
            continue
        l5 = pari.subst(N, "t", 0) / pari.subst(D, "t", 0)
        l6 = pari.polcoef(N, 2) / pari.polcoef(D, 2)
        # [D] lam5 lam6 and (1 - lam5)(1 - lam6) are products of four of the ten classes (eps cancels)
        i1, i2, i3, i4_ = (x - 1 for x in rest)
        i5, i6 = pair[0] - 1, pair[1] - 1

        def cls(a, b, c):
            return beta[_split_of(tuple(sorted((a, b, c))))]
        pred1 = cls(i1, i3, i5) * cls(i1, i3, i6) * cls(i1, i4_, i5) * cls(i1, i4_, i6)
        pred2 = cls(i1, i2, i5) * cls(i1, i2, i6) * cls(i1, i4_, i5) * cls(i1, i4_, i6)
        assert eta(l5 * l6) == pred1 and eta((1 - l5) * (1 - l6)) == pred2, "lam5 lam6 is not the product of four classes"
        assert bool(int(cur[6])) == (branch[-1] == 1), "the two descriptions of the branch values differ"
        c4.add(int(cur[5]))
    out["branch"] = sum(1 for x in branch if x == 1)
    out["branch_zero"] = sum(1 for x in branch if x == 0)
    out["c4_rational"] = sorted(c4)
    return out


def load_torsors():
    tors = []
    for line in open(os.path.join(HERE, "sixlines_cm.jsonl")):
        t = json.loads(line)
        t["q"] = t["p"]
        tors.append(t)
    for line in open(os.path.join(HERE, "sixlines_neighbours.jsonl")):
        t = json.loads(line)
        if t["p"] == 25:
            t["q"] = 25
            tors.append(t)
    for line in open(os.path.join(HERE, "sixlines_q49_cm.jsonl")):
        t = json.loads(line)
        t["q"] = 49
        tors.append(t)
    return tors


def load_rows():
    rows = {}
    for q in (13, 17, 19, 23, 29, 31, 37):
        for line in open(os.path.join(HERE, f"sixlines_p{q}.jsonl")):
            r = json.loads(line)
            rows[(q, tuple(r["cfg"]))] = r
    for q in (25, 49):
        for line in open(os.path.join(HERE, f"sixlines_q{q}.jsonl")):
            r = json.loads(line)
            rows[(q, tuple(r["cfg"]))] = r
    return rows


def load_quartics():
    out = {}
    for name, key in (("sixlines_prym.jsonl", "p"), ("sixlines_prym_q25.jsonl", "q"), ("sixlines_prym_q49.jsonl", "q")):
        for line in open(os.path.join(HERE, name)):
            r = json.loads(line)
            out[(r[key], tuple(r["cfg"]))] = r
    return out


def show(v):
    return "".join("+" if x > 0 else "-" for x in v)


def shape(counter, n):
    """constant / balanced / neither, for the values of a function on a torsor of n surfaces that should be +-1"""
    if len(counter) == 1:
        return "constant"
    if set(counter) == {1, -1} and counter[1] == counter[-1]:
        return "balanced"
    return "neither"


INVARIANTS = ("type", "conic", "j2", "j4", "branch", "quartic", "quartic_type", "orientation")
BITS = ("conic", "j2", "j4", "quartic")


def cmd_surfaces(args):
    rng = random.Random(75)
    tors, rows, quart = load_torsors(), load_rows(), load_quartics()
    fh = open(args.out, "w") if args.out else None
    st = Counter()
    assoc_ok, relabel_ok, frob_ok = Counter(), Counter(), Counter()
    per_torsor, joint = [], {}
    type_rule, branch_rule = defaultdict(set), defaultdict(set)
    c1_stats, orient = Counter(), Counter()
    for t in tors:
        q, e = t["q"], t["e"]
        maximal = t["index"] == 1
        grp = group_record(q, tuple(e), identify=False) if maximal else None
        val = {}
        for cfg in map(tuple, t["configs"]):
            r = rows[(q, cfg)]
            L = forms(cfg, q)
            # the surface of the torsor is w^2 = l1 ... l6 if its trace is e1, and the quadratic twist otherwise
            twist = r["S"][0] != e[0]
            assert not twist or r["S"][0] == -e[0]
            c = candidates(q, L, twist)
            val[cfg] = c
            branch_rule[show(c["type"])].add((c["branch"], c["branch_zero"]))
            assert c["c4_rational"] in ([], [1 if q % 4 == 1 else 0]), "c4 rational does not follow q mod 4"
            # relabelling: the same candidates from the forms in another order (one random order each)
            c2 = candidates(q, [L[i] for i in rng.choice(PERMS[1:])], twist)
            for name in INVARIANTS:
                if name in c:
                    relabel_ok[(name, q % 4, c[name] == c2[name])] += 1
            # the stored quartics: values of c1 over the pairs in which a, b, c are rational
            qr = quart.get((q, cfg))
            c1 = sorted({tuple(x["c"][0]) for x in qr["pairs"] if x.get("rational_point") == 0 and "c" in x}) if qr else []
            if qr and sum(1 for x in qr["pairs"] if x.get("rational_point") == 0 and "c" in x) >= 2:
                c1_stats["configurations with two or more pairs in which a, b, c are rational"] += 1
                c1_stats["of these: c1 takes two or more values"] += len(c1) >= 2
                c1_stats["of these: c1 and -c1 both occur"] += any((-a, -b) in c1 for a, b in c1)
                c1_stats["of these: c1 and its conjugate both occur"] += any((a, -b) in c1 and b != 0 for a, b in c1)
            if fh:
                rec = {"q": q, "e": e, "cfg": list(cfg), "twist": twist, "type": show(c["type"]), "c1_values": c1}
                rec.update({k: c[k] for k in ("conic", "j2", "j4", "branch", "branch_zero", "quartic", "orientation") if k in c})
                if "quartic_type" in c:
                    rec["quartic_type"] = "".join(map(str, c["quartic_type"]))
                fh.write(json.dumps(rec) + "\n")
        # association and Frobenius conjugation
        k3, seen = [], set()
        for cfg in val:
            a = tuple(rows[(q, cfg)]["assoc"])
            for name in INVARIANTS:
                if name in val[cfg]:
                    assoc_ok[(name, q % 4, val[cfg][name] == val[a][name])] += 1
            if "frob" in rows[(q, cfg)]:
                fcfg = tuple(rows[(q, cfg)]["frob"])
                for name in INVARIANTS + ("quartic_oriented",):
                    if name in val[cfg]:
                        frob_ok[(name, q, val[cfg][name] == val[fcfg][name])] += 1
            key = min(cfg, a)
            if key not in seen:
                seen.add(key)
                k3.append(key)
        n = len(k3)
        assert n == round(t["k3"])
        st["polynomials"] += 1
        st["surfaces"] += n
        types = Counter(val[c]["type"] for c in k3)
        line = {"q": q, "e": e, "index": t["index"], "surfaces": n, "types": {show(k): v for k, v in types.items()}}
        if not maximal:
            st["not maximal: polynomials"] += 1
            st["not maximal: one type on the polynomial"] += len(types) == 1
            st["not maximal, odd index: polynomials"] += t["index"] % 2
            st["not maximal, odd index: one type on the polynomial"] += t["index"] % 2 == 1 and len(types) == 1
            per_torsor.append(line)
            continue
        st["maximal: polynomials"] += 1
        st["maximal: surfaces"] += n
        st["maximal: polynomials with |G| even"] += n % 2 == 0
        assert grp["G"] == n, "the number of surfaces is not |G|"
        line.update({"G_cyc": grp["G_cyc"], "rank2": grp["rank2"], "frobenius_trivial": grp["frobenius_trivial"],
                     "frobenius_on_characters": grp["frobenius_on_characters"]})
        st["maximal: one type on the torsor"] += len(types) == 1
        type_rule[(q % 16, e[0] % 16)].add(show(next(iter(types))) if len(types) == 1 else "several")
        line["bits"] = {}
        for name in BITS:
            if name not in val[k3[0]] or (name == "conic" and q % 4 == 3):
                continue                                        # the conic class is not an invariant of the surface there
            cnt = Counter(val[c][name] for c in k3)
            line["bits"][name] = [cnt.get(1, 0), cnt.get(-1, 0), cnt.get(0, 0), shape(cnt, n)]
        if q % 4 == 1:
            line["quartic_types"] = len({val[c]["quartic_type"] for c in k3})
            st["maximal, q = 1 mod 4: surfaces"] += n
            st["maximal, q = 1 mod 4: quartic types, summed over the torsors"] += line["quartic_types"]
            if n >= 2:
                joint[(q, tuple(e))] = [tuple(val[c][b] for b in BITS) for c in k3]
            if q == 49:
                for c in k3:
                    orient[val[c]["orientation"]] += 1
        per_torsor.append(line)

    print("surfaces:", dict(st))
    print("\nthe same value on a configuration and its associate (candidate, q mod 4): agree / differ")
    for name in INVARIANTS:
        for r in (1, 3):
            a, b = assoc_ok[(name, r, True)], assoc_ok[(name, r, False)]
            if a + b:
                print(f"  {name:13s} q = {r} mod 4: {a} / {b}")
    print("the same value after a random relabelling of the six lines: agree / differ")
    for name in INVARIANTS:
        for r in (1, 3):
            a, b = relabel_ok[(name, r, True)], relabel_ok[(name, r, False)]
            if a + b:
                print(f"  {name:13s} q = {r} mod 4: {a} / {b}")
    print("the same value on the Frobenius conjugate: agree / differ")
    for name in INVARIANTS + ("quartic_oriented",):
        print(f"  {name:16s} F_25: {frob_ok[(name, 25, True)]} / {frob_ok[(name, 25, False)]}   "
              f"F_49: {frob_ok[(name, 49, True)]} / {frob_ok[(name, 49, False)]}")
    print("orientation over F_49, surfaces of maximal orders:", dict(sorted(orient.items())))
    print("\nc1 of the stored quartics (notes 14, 16):", dict(c1_stats))

    print("\ntype of the ten square classes as a function of (q mod 16, e1 mod 16), maximal orders:")
    for key in sorted(type_rule):
        print(f"  q = {key[0]:2d}, e1 = {key[1]:2d} mod 16: {sorted(type_rule[key])}")
    assert all(len(v) == 1 for v in type_rule.values()), "the type is not a function of (q mod 16, e1 mod 16)"
    print("number of pairs of lines with rational branch values (and with a double branch value), by type, all 76 polynomials:")
    for key in sorted(branch_rule):
        print(f"  {key}: {sorted(branch_rule[key])}")
    assert all(len(v) == 1 for v in branch_rule.values()), "the branch count is not a function of the type"

    print("\nmaximal orders: q, e, G, [values on the characters of Frobenius conjugation], type, quartic types; "
          "then for each bit +1 / -1 / 0 and the shape")
    summary = defaultdict(Counter)
    for line in per_torsor:
        if line["index"] != 1:
            continue
        bits = "  ".join(f"{k} {v[0]}/{v[1]}/{v[2]} {v[3]}" for k, v in line["bits"].items())
        fro = "" if line["q"] in (13, 17, 19, 23, 29, 31, 37) else f" frob {line['frobenius_on_characters']}"
        print(f"  {line['q']:2d} {str(line['e']):22s} G {str(line['G_cyc']):10s}{fro} {list(line['types'])[0]} "
              f"{line.get('quartic_types', '-')}  {bits}")
        for k, v in line["bits"].items():
            if line["surfaces"] >= 2:
                summary[k][v[3] + (", |G| even" if line["surfaces"] % 2 == 0 else ", |G| odd")] += 1
    print("\nshape of each bit over the torsors of maximal orders with at least two surfaces:")
    for k in BITS:
        print(f"  {k:8s} {dict(sorted(summary[k].items()))}")

    # is any function of (conic, j2, j4, quartic) constant or balanced on every torsor ?  (q = 1 mod 4 only, where
    # all four are invariants of the surface)
    keys = sorted({v for vals in joint.values() for v in vals})
    found = []
    for mask in range(1, 2 ** (len(keys) - 1)):                 # f and -f are the same test
        f = {k: (1 if mask >> i & 1 else -1) for i, k in enumerate(keys)}
        moved = []
        for key, vals in joint.items():
            s = sum(f[v] for v in vals)
            if abs(s) == len(vals):
                continue
            if s != 0:
                break
            moved.append(key)
        else:
            found.append(([k for k in keys if f[k] == 1], moved))
    print(f"\nfunctions of (conic, j2, j4, quartic), up to sign, on the {len(keys)} value combinations that occur over the "
          f"{len(joint)} torsors with q = 1 mod 4 and at least two surfaces: {2 ** (len(keys) - 1) - 1} non-constant")
    print(f"  constant or balanced on every torsor: {len(found)}")
    for plus, moved in found:
        print(f"    +1 on {plus}; balanced on {len(moved)} torsors: {moved}")
    if fh:
        fh.close()


def cmd_dyadic(args):
    """Is "d = -1 gives a character of G" the same as "K+ has a prime above 2 of odd local degree" ?  (note 19, 2.4)"""
    table = Counter()
    for line in open(args.records):
        rec = json.loads(line)
        if "failed" in rec or not rec["rank2"]:
            continue
        dy = GPC["ch_dyadic"](rec["q"], rec["e"][0])
        if dy == 0:
            table["chi conj(chi) reducible"] += 1
            continue
        shape_ = tuple(sorted((int(x[0]), int(x[1])) for x in dy))
        odd = any(e * f % 2 for e, f in shape_)
        table[(rec["rational"]["-1"], odd, shape_)] += 1
    print("dyadic: polynomials with a character of order 2, by (d = -1 gives a character, K+ has a dyadic prime of odd "
          "local degree, [e, f] of the dyadic primes of K+):")
    for key, n in sorted(table.items(), key=str):
        print(f"  {key}: {n}")
    agree = sum(n for key, n in table.items() if isinstance(key, tuple) and key[0] == key[1])
    total = sum(n for key, n in table.items() if isinstance(key, tuple))
    print(f"  the two conditions agree in {agree} of {total}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("genus")
    g.add_argument("--workers", type=int, default=8)
    g.add_argument("--max-q", type=int, default=10 ** 9)
    g.add_argument("--min-q", type=int, default=0)
    g.add_argument("--summary", nargs="+", metavar="FILE", help="recompute the statistics from output files")
    g.add_argument("--out")
    g.set_defaults(func=cmd_genus)
    d = sub.add_parser("dyadic")
    d.add_argument("records", help="the output file of genus")
    d.set_defaults(func=cmd_dyadic)
    s = sub.add_parser("surfaces")
    s.add_argument("--out")
    s.set_defaults(func=cmd_surfaces)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
