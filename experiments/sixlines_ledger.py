"""Hypothesis ledger for the hardness question for six-line K3 surfaces (note 27).

Reads only stored .jsonl / .json / .log files of earlier notes and recomputes,
for each hypothesis with stored evidence, the population, the count that
agrees, the exceptions and the range of q. Where a note quoted a number the
script compares it with the recomputed one and prints MISMATCH when they
differ. It computes no surface, curve, lattice or CM point; the only
arithmetic is on stored polynomial coefficients and stored counts (and one
call of PARI's polisirreducible on stored polynomials, for H1).

"quoted" is the number as a note first stated it. Where a check fails, the note
has since been corrected to the recomputed value (note 27 §1 lists the cases).

Usage, from experiments/:
    python -I sixlines_ledger.py > sixlines_ledger.log
"""

import collections
import json
import math
import os
import re
import statistics
import sys
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
FILES_USED = collections.OrderedDict()
CHECKS = []          # (hid, what, recomputed, quoted, ok)
CURRENT = [None]


def path(name):
    return os.path.join(HERE, name)


def load(name):
    """All records of a .jsonl file, in file order."""
    out = []
    with open(path(name)) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    FILES_USED[name] = len(out)
    return out


def load_json(name):
    with open(path(name)) as fh:
        d = json.load(fh)
    FILES_USED[name] = 1
    return d


def read_log(name):
    with open(path(name)) as fh:
        s = fh.read()
    FILES_USED[name] = s.count("\n")
    return s


def head(hid, title):
    CURRENT[0] = hid
    print()
    print("=" * 78)
    print(f"{hid}  {title}")
    print("=" * 78)


def info(*a):
    print("   ", *a)


def check(what, recomputed, quoted, note=""):
    ok = recomputed == quoted
    CHECKS.append((CURRENT[0], what, recomputed, quoted, ok, note))
    tag = "ok      " if ok else "MISMATCH"
    extra = f"  ({note})" if note else ""
    print(f"    [{tag}] {what}: recomputed {recomputed!r}, quoted {quoted!r}{extra}")
    return ok


def check_close(what, recomputed, quoted, tol, note=""):
    ok = abs(recomputed - quoted) <= tol
    CHECKS.append((CURRENT[0], what, round(recomputed, 6), quoted, ok, note))
    tag = "ok      " if ok else "MISMATCH"
    extra = f"  ({note})" if note else ""
    print(f"    [{tag}] {what}: recomputed {recomputed:.4g}, quoted {quoted}{extra}")
    return ok


def qrange(qs):
    qs = sorted(set(qs))
    return f"q in {{{', '.join(map(str, qs))}}}" if len(qs) <= 12 else f"q from {qs[0]} to {qs[-1]} ({len(qs)} values)"


# ---------------------------------------------------------------- arithmetic

def P_at(e, q, t):
    e1, e2, e3 = e
    return (t**6 - e1 * t**5 + e2 * t**4 - e3 * t**3
            + q * q * e2 * t**2 - q**4 * e1 * t + q**6)


def isqrt_exact(n):
    if n < 0:
        return None
    r = math.isqrt(n)
    return r if r * r == n else None


def mu_nu(e, q):
    """(mu, nu) >= 0 with P(q) = 64 q^5 mu^2, P(-q) = 64 q^5 nu^2, or None."""
    a, b = P_at(e, q, q), P_at(e, q, -q)
    d = 64 * q**5
    if a % d or b % d:
        return None
    m, n = isqrt_exact(a // d), isqrt_exact(b // d)
    if m is None or n is None:
        return None
    return m, n


def wedge2_of_chi(e1, mu, nu, q):
    """e2, e3 of the second exterior power of Proposition Y's quartic (note 13)."""
    # chi = X^4 - a X^3 + b X^2 - d X + c^2, a = 4(nu + mu i), d = q conj(a), b = e1, c = q
    a = complex(4 * nu, 4 * mu)
    d = q * a.conjugate()
    c2 = q * q
    e2 = (a * d - c2)
    e3 = (a * a * c2 + d * d - 2 * e1 * c2)
    return e2, e3


def wedge2_quartic(c, q):
    """(e1, e2, e3) of the polynomial of products of two roots of
    chi = X^4 - c1 X^3 + c2 X^2 - c3 X + q^2 with Gaussian-integer c's
    (note 14 §3); returned as complex numbers."""
    c1, c2, c3 = (complex(*x) for x in c)
    return c2, c1 * c3 - q * q, c1 * c1 * q * q + c3 * c3 - 2 * c2 * q * q


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def largest_prime_factor(n):
    best, d = 1, 2
    while d * d <= n:
        while n % d == 0:
            best, n = d, n // d
        d += 1
    return max(best, n) if n > 1 else best


def quantiles(xs):
    xs = sorted(xs)
    q = statistics.quantiles(xs, n=4, method="inclusive")
    return q[0], statistics.median(xs), q[2]


def lsq_slope(xs, ys):
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    return b, my - b * mx


def balanced_chance(sizes):
    """Probability that a function balanced on torsors of these sizes agrees,
    up to sign on each torsor, with a given split into halves."""
    pr = Fraction(1)
    for n in sizes:
        pr *= Fraction(2, comb(n, n // 2))
    return pr


def balanced_prob(sizes):
    """Probability that a uniformly random +-1 function is balanced on all."""
    pr = Fraction(1)
    for n in sizes:
        pr *= Fraction(comb(n, n // 2), 2**n)
    return pr


def canon_e(e):
    """Polynomial up to t -> -t (the quadratic twist)."""
    e1, e2, e3 = e
    return (abs(e1), e2, e3 if e1 > 0 else -e3) if e1 != 0 else (0, e2, abs(e3))


# =================================================================== group a

def group_a():
    print("\n\n######## GROUP G1: the reduction to abelian fourfolds ########")

    # ---------------------------------------------------------------- H1
    head("H1", "Proposition Y: P(+-q) = 64 q^5 (mu^2, nu^2) and P = wedge^2 of chi (note 13 §4.2)")
    gs = load("sixlines_groupsize.jsonl")
    cm = load("sixlines_cm.jsonl")
    n = ok_shape = ok_wedge = 0
    qs = []
    bad = []
    seen = set()
    for r in gs:
        key = (r["q"], tuple(r["e"]))
        if key in seen:
            continue
        seen.add(key)
        n += 1
        q, e = r["q"], r["e"]
        mn = mu_nu(e, q)
        if mn is None:
            bad.append((q, e, "shape"))
            continue
        ok_shape += 1
        mu, nu = mn
        e2, e3 = wedge2_of_chi(e[0], mu, nu, q)
        if e2 == e[1] and e3 == e[2]:
            ok_wedge += 1
            qs.append(q)
        else:
            bad.append((q, e, "wedge"))
    info(f"files: sixlines_groupsize.jsonl; population: {n} distinct stored polynomials "
         f"(admissible ones of the 'all' runs, the 'sample' runs and the 24 'real'), {qrange([k[0] for k in seen])}")
    info(f"P(q) = 64q^5 mu^2 and P(-q) = 64q^5 nu^2 with integers: {ok_shape} of {n}")
    info(f"second exterior power of chi equals P: {ok_wedge} of {n}; exceptions: {bad[:5]}")
    n44 = sum(1 for r in cm if mu_nu(r["e"], r["p"]) is not None)
    w44 = 0
    for r in cm:
        mn = mu_nu(r["e"], r["p"])
        if mn and wedge2_of_chi(r["e"][0], mn[0], mn[1], r["p"]) == (r["e"][1], r["e"][2]):
            w44 += 1
    check("enumerated polynomials (sixlines_cm.jsonl) with the shape and P = wedge^2 chi", (n44, w44), (44, 44))
    log = read_log("sixlines_groupsize.log")
    m = re.search(r"fourfold: (\{.*?\})", log)
    four = eval(m.group(1)) if m else {}
    check("stored fourfold line: real polynomials with chi conj(chi) irreducible", four.get("real, degree-8 polynomial irreducible"), 32)
    info(f"stored fourfold line of the log: {four}")
    info("NOT IN STORED FILES: the 23,945 polynomials of note 13 §4.2 (q = 25, 49, 121, 125, 169 and the primes 41..257) "
         "include the 265 of --fourfold 41 61 101, which are in the stored log; the other 23,680 come from nine runs "
         "whose output was not kept. "
         "Rerun: sixlines_groupsize.py --fourfold 25 43 47 ... 257 (note 13, Reproduce).")
    try:
        import cypari2
        pari = cypari2.Pari()
        irr = 0
        for r in cm:
            mn = mu_nu(r["e"], r["p"])
            mu, nu, q = mn[0], mn[1], r["p"]
            chi = f"x^4 - 4*({nu} + {mu}*I)*x^3 + {r['e'][0]}*x^2 - 4*{q}*({nu} - {mu}*I)*x + {q*q}"
            chibar = f"x^4 - 4*({nu} - {mu}*I)*x^3 + {r['e'][0]}*x^2 - 4*{q}*({nu} + {mu}*I)*x + {q*q}"
            f8 = pari(f"real(({chi})*({chibar}))")
            irr += int(pari.polisirreducible(f8))
        check("chi conj(chi) irreducible among the 44 enumerated polynomials (PARI on stored e)", irr, 32)
    except Exception as exc:  # pragma: no cover
        info(f"PARI not available: {exc}")

    # ---------------------------------------------------------------- H2
    head("H2", "The curve of note 14 gives P or P(-t) for every surface with irreducible P (notes 14, 16, 26)")
    pr = load("sixlines_prym.jsonl")
    tab = collections.defaultdict(collections.Counter)
    pairs_signplus = []
    for r in pr:
        p = r["p"]
        t = tab[p]
        t["configurations"] += 1
        sign = r["signs"][0] if r["signs"] else None
        irr = r.get("degrees") == [6]
        comp = [x for x in r["pairs"]]
        if sign == 1 and irr:
            t["sign+ irr"] += 1
            t["sign+ irr with pair"] += bool(comp)
            t["sign+ irr pairs"] += len(comp)
            t["sign+ irr inverting seen"] += r["seen"].get("anticommuting", 0)
            t["sign+ irr degenerate seen"] += r["seen"].get("degenerate", 0)
        elif sign == 1:
            t["sign+ red"] += 1
            t["sign+ red with pair"] += bool(comp)
            t["sign+ red pairs"] += len(comp)
        else:
            t["sign-"] += 1
        if sign == 1:
            for x in comp:
                t["P"] += x["compare"] == "equal"
                t["P(-t)"] += x["compare"] == "twist"
                t["other"] += x["compare"] not in ("equal", "twist")
                pairs_signplus.append((r, x))
    tot = collections.Counter()
    for p in sorted(tab):
        tot.update(tab[p])
        info(f"p={p:>2}: " + ", ".join(f"{k} {v}" for k, v in tab[p].items()))
    quoted = {"configurations": 3990, "sign+ irr": 534, "sign+ irr with pair": 534, "sign+ irr pairs": 4983,
              "sign+ red": 1564, "sign+ red with pair": 1329, "sign+ red pairs": 12873, "sign-": 1892,
              "P": 9048, "P(-t)": 8808, "other": 0}
    for k, v in quoted.items():
        check(f"note 14 §4 table, column '{k}' (all p)", tot[k], v)
    check("note 18 §2: inverting pairs seen on the 534 configurations with irreducible P", tot["sign+ irr inverting seen"], 0)
    check("note 18 §2: degenerate pairs seen on the 534", tot["sign+ irr degenerate seen"], 809)
    # q25 / q49 / p41
    for name, q, quoted_ in (("sixlines_prym_q25.jsonl", 25, (37, 858, 131, 2829, 30)),
                             ("sixlines_prym_q49.jsonl", 49, (1108, 3324, None, None, 1060))):
        rr = load(name)
        irr = [r for r in rr if r["signs"] and r["signs"][0] == 1 and r.get("degrees") == [6]]
        red = [r for r in rr if r["signs"] and r["signs"][0] == 1 and r.get("degrees") not in (None, [6]) and r["pairs"]]
        npairs = sum(len(r["pairs"]) for r in irr)
        good = sum(1 for r in irr for x in r["pairs"] if x["compare"] in ("equal", "twist"))
        nred = sum(len(r["pairs"]) for r in red)
        goodred = sum(1 for r in red for x in r["pairs"] if x["compare"] in ("equal", "twist"))
        notprime = sum(1 for r in irr if not r["defined_over_prime_field"])
        info(f"{name}: sign + irreducible {len(irr)} configurations, {npairs} pairs, {good} give P or P(-t); "
             f"reducible with a pair {len(red)}, {nred} pairs, {goodred} good; not defined over F_p {notprime}")
        check(f"note 16, q={q}: configurations, pairs", (len(irr), npairs), quoted_[:2])
        check(f"note 16, q={q}: pairs giving P or P(-t)", good, npairs)
        if quoted_[2] is not None:
            check(f"note 16, q={q}: reducible P, surfaces with a pair, pairs", (len(red), nred), quoted_[2:4])
        check(f"note 16, q={q}: not isomorphic to the Frobenius conjugate", notprime, quoted_[4])
    r41 = load("sixlines_prym_p41.jsonl")
    n41 = sum(len(r["pairs"]) for r in r41)
    eq41 = sum(1 for r in r41 for x in r["pairs"] if x["compare"] == "equal")
    tw41 = sum(1 for r in r41 for x in r["pairs"] if x["compare"] == "twist")
    deg41 = sum(r["seen"].get("degenerate", 0) for r in r41)
    check("note 26 §4, p=41: configurations, computed, P, P(-t), degenerate", (len(r41), n41, eq41, tw41, deg41),
          (56, 1507, 763, 744, 173))
    surf = {}
    for p in (13, 17, 19, 23, 29, 31, 37):
        for r in load(f"sixlines_p{p}.jsonl"):
            surf[(p, tuple(r["cfg"]))] = min(tuple(r["cfg"]), tuple(r["assoc"]))
    irr534 = [r for r in pr if r["signs"] and r["signs"][0] == 1 and r.get("degrees") == [6]]
    nsurf = len({(r["p"], surf[(r["p"], tuple(r["cfg"]))]) for r in irr534})
    check("note 14 §1 item 2 says '534 surfaces': surfaces (a configuration and its associate) among the 534",
          nsurf, 534, "534 counts configurations, as the table of §4 says; the surfaces number 278")
    for q, f_, pf, quoted_s in ((25, "sixlines_q25.jsonl", "sixlines_prym_q25.jsonl", 37), (49, "sixlines_q49.jsonl", "sixlines_prym_q49.jsonl", 1108)):
        sq_ = {tuple(r["cfg"]): min(tuple(r["cfg"]), tuple(r["assoc"])) for r in load(f_)}
        irr_ = [r for r in load(pf) if r["signs"] and r["signs"][0] == 1 and r.get("degrees") == [6]]
        check(f"note 16 §1 item 2, q={q}: 'surfaces' with irreducible P, counted as surfaces",
              len({sq_[tuple(r['cfg'])] for r in irr_}), quoted_s, "the quoted number counts configurations")
    total_pairs = tot["sign+ irr pairs"] + 858 + 3324 + n41
    info(f"TOTAL for sign + and irreducible P: {534 + 37 + 1108 + len(r41)} configurations, {nsurf} + 21 + 564 + 28 surfaces "
         f"(p 13..37, q 25, 49, p 41), {total_pairs} pairs of lines, 0 exceptions")

    # ---------------------------------------------------------------- H3
    head("H3", "chi of the curve is Proposition Y's quartic: c3 = p conj(c1), c1 = 0 mod 4, {|Re c1|,|Im c1|} = {4mu,4nu}")
    n = a = b = c = 0
    excs = []
    for r, x in pairs_signplus:
        n += 1
        c1 = x["c"][0]
        a += bool(x.get("c3_is_p_conj_c1"))
        b += (c1[0] % 4 == 0 and c1[1] % 4 == 0)
        mn = mu_nu(r["e"], r["p"])
        if mn is not None and sorted([abs(c1[0]), abs(c1[1])]) == sorted([4 * mn[0], 4 * mn[1]]):
            c += 1
        elif len(excs) < 3:
            excs.append((r["p"], r["cfg"], x["pair"], c1, mn))
    info(f"files: sixlines_prym.jsonl; population: {n} computed pairs on surfaces of sign + (p = 11..37)")
    check("pairs with c3 = p conj(c1)", a, 17856)
    check("pairs with c1 = 0 mod 4", b, 17856)
    check("pairs with {|Re c1|, |Im c1|} = {4mu, 4nu}, mu, nu from P(+-p) of the stored e", c, 17856,
          f"exceptions {excs}" if excs else "")
    check("population", n, 17856)

    # ---------------------------------------------------------------- H4
    head("H4", "Sign + and irreducible P: no anomalous (inverting) pair, no collinear triple (notes 14 §5.2, 21 §8)")
    info("files: sixlines_prym.jsonl (seen counts), sixlines_forms_fp.log (types run)")
    flog = read_log("sixlines_forms_fp.log")
    m = re.findall(r"(?im)^.*(?:irreducible).*$", flog)
    for line in m[:6]:
        info("log: " + line.strip()[:160])
    info("note 21 quotes 0 of 534 with a collinear triple and 14,394 non-degenerate pairs, none anomalous;"
         " the per-pair data of the 'types' run are in the log only (lines above), not in a .jsonl.")

    # ---------------------------------------------------------------- H5
    head("H5", "Dictionary Cl(K) -> Cl(E) (note 15 §3, corrected by note 18 findings 1, 2)")
    dic = load("sixlines_dictionary.jsonl")
    ok = [r for r in dic if r["status"] == "ok"]
    check("polynomials treated / skipped as chi conj(chi) reducible", (len(ok), len(dic) - len(ok)), (278, 12))
    check("homomorphism found", sum(r["homomorphism"] for r in ok), 278)
    check("f(Cl(K+)) in H", sum(r["subgroups_respected"] for r in ok), 278)
    check("odd parts equal", sum(r["odd_parts_equal"] for r in ok), 278)
    check("odd part non-trivial", sum(1 for r in ok if odd_part(r["G"]) > 1), 224)
    check("primes of K checked", sum(r["primes"] for r in ok), 77946)
    kc = collections.Counter(r["kernel"] for r in ok)
    cc = collections.Counter(r["cokernel"] for r in ok)
    check("kernel orders 1, 2, 4, 8", tuple(kc[i] for i in (1, 2, 4, 8)), (6, 205, 64, 3))
    check("cokernel orders 1, 2, 4, 8", tuple(cc[i] for i in (1, 2, 4, 8)), (10, 148, 105, 15))
    check("Cl^-(K) and G of the same order", sum(1 for r in ok if r["minus_K"] == r["G"]), 209)
    check("same order implies same elementary divisors",
          sum(1 for r in ok if r["minus_K"] == r["G"] and r["minus_K_cyc"] == r["G_cyc"]), 209)
    big = max(ok, key=lambda r: r["G"])
    check("largest G (order, q, e)", (big["G"], big["q"], big["e"]), (246, 101, [10, -6969, -204020]))
    rlog = read_log("sixlines_review.log")
    m = re.search(r"elementary divisors of the cokernel: (\{.*?\})", rlog)
    cok = eval(m.group(1))
    order4 = sum(v for k, v in cok.items() if 4 in k)
    check("note 18 finding 1: cokernel with an element of order 4 (sixlines_review.log)", order4, 46)
    info(f"{qrange([r['q'] for r in ok])}; all checks are on class groups (GRH), none on a surface")

    # ---------------------------------------------------------------- H6
    head("H6", "Association is invisible to the curve (note 15 §4)")
    alog = read_log("sixlines_dictionary_association.log")
    m = re.search(r"total: (\{.*?\})", alog)
    t = eval(m.group(1))
    check("pairs giving the associated configuration / skipped / self-associated",
          (t["associated configuration"], t["no conjugate lines"], t["of which self-associated"]), (52207, 7643, 26867))
    check("exceptions (pairs - associated - skipped)", t["pairs"] - t["associated configuration"] - t["no conjugate lines"], 0)

    # ---------------------------------------------------------------- H7
    head("H7", "CM method: the class polynomial reduces to exactly the enumerated surfaces (notes 20, 24, 26)")
    cmm = load("sixlines_cmmethod.jsonl")
    last = collections.OrderedDict()
    for r in cmm:
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r   # rule of note 20: the last record of a field and command
    info("rule applied (note 20, Reproduce): for sixlines_cmmethod.jsonl the last record of each (field, command) is used;"
         f" {len(cmm)} records reduce to {len(last)}")
    comp = [r for k, r in last.items() if k[2] == "compare" and r["status"] == "ok"]
    act = [r for k, r in last.items() if k[2] == "action" and r["status"] == "ok"]
    first = collections.OrderedDict()
    for r in cmm:
        first.setdefault((r["q"], tuple(r["e"]), r["cmd"]), r)
    fcomp = [r for k, r in first.items() if k[2] == "compare" and r["status"] == "ok"]
    info(f"compare records used: {len(comp)} fields; with the FIRST record instead: {len(fcomp)} fields, "
         f"in_both {sum(r['in_both'] for r in fcomp)} of {sum(r['enumerated_points'] for r in fcomp)}")
    check("fields compared", len(comp), 15)
    check("enumerated surfaces, CM points mod p, in both, one side only",
          (sum(r["enumerated_points"] for r in comp), sum(r["cm_points_mod_q"] for r in comp),
           sum(r["in_both"] for r in comp), sum(r["enumerated_only"] + r["cm_only"] for r in comp)), (74, 74, 74, 0))
    check("configurations", sum(r["enumerated_configurations"] for r in comp), 148)
    check("even lattices per type = surfaces, in every field", sum(r["even_per_type"] == r["enumerated_points"] for r in comp), 15)
    check("even / odd lattices per type, totals", (sum(r["even_per_type"] for r in comp), sum(r["odd_per_type"] for r in comp)), (74, 165))
    dl = load("sixlines_delta_labels.jsonl")
    d41 = load("sixlines_delta41_cosets.jsonl")
    d41last = collections.OrderedDict(((r["q"], r["e1"]), r) for r in d41 if "surfaces_in_both" in r)
    more = [(r["q"], r["e1"], r["surfaces_in_both"], r["enumerated_points"], r["N"]) for r in dl] + \
           [(r["q"], r["e1"], r["surfaces_in_both"], r["enumerated_points"], r["N"]) for r in d41last.values()]
    info(f"further fields (sixlines_delta_labels.jsonl, sixlines_delta41_cosets.jsonl, last record per field): {more}")
    tot_in = sum(r["in_both"] for r in comp) + sum(x[2] for x in more)
    tot_en = sum(r["enumerated_points"] for r in comp) + sum(x[3] for x in more)
    nfields = len(comp) + len(more)
    info(f"records summed (q=37, e1=18 appears in note 20 and again in note 24): {tot_in} of {tot_en} in {nfields} records")
    tot91 = sum(r["in_both"] for r in comp) + sum(x[2] for x in more if (x[0], x[1]) in ((29, 2), (37, 2)))
    check("note 24: surfaces of the CM method in 17 fields", tot91, 98)
    distinct = {(r["q"], r["e"][0]) for r in comp} | {(x[0], x[1]) for x in more}
    tot_in_d = sum(r["in_both"] for r in comp if (r["q"], r["e"][0]) != (37, 18)) + sum(x[2] for x in more)
    check("distinct fields and surfaces", (len(distinct), tot_in_d), (19, 126))
    gsizes = sorted({x[4] for x in more} | {r["enumerated_points"] for r in comp})
    info(f"largest torsor reached: |G| = {max(gsizes)} ; q range {qrange([x[0] for x in distinct])}")

    # ---------------------------------------------------------------- H8
    head("H8", "The action on lattices is G, freely and transitively; the step is Frobenius in characteristic l (note 20 §7)")
    keys = ["ideals", "conflicts", "classes_seen", "products_ok", "products", "fixed_points",
            "frobenius_ok", "frobenius_tested", "frobenius_skipped", "inverse_fails", "inverse_tested",
            "mismatch_fails", "mismatch_tested", "other_ideals", "other_in_G", "other_square_tested",
            "other_square_ok", "labelled"]
    s = {k: sum(r[k] for r in act) for k in keys}
    gtot = sum(math.prod(r["G"]) for r in act)
    check("fields with an action record", len(act), 15)
    check("ideals, conflicts, fixed points", (s["ideals"], s["conflicts"], s["fixed_points"]), (378, 0, 0))
    check("classes seen of |G| summed", (s["classes_seen"], gtot), (69, 74))
    check("sums that compose", (s["products_ok"], s["products"]), (390, 390))
    check("Frobenius congruence holds / tested / not tested", (s["frobenius_ok"], s["frobenius_tested"], s["frobenius_skipped"]), (346, 346, 32))
    check("inverse permutation fails / tested", (s["inverse_fails"], s["inverse_tested"]), (136, 136))
    check("another pair fails / tested", (s["mismatch_fails"], s["mismatch_tested"]), (1156, 1156))
    check("other ideals: total, permutation of a class seen, square tested, square ok",
          (s["other_ideals"], s["other_in_G"], s["other_square_tested"], s["other_square_ok"]), (630, 610, 594, 594))
    check("surfaces labelled", s["labelled"], 67)
    lc = [r["label_consistency"] for r in act]
    check("label consistency (ok, tested) summed", (sum(x[0] for x in lc), sum(x[1] for x in lc)), (389, 389))
    trans = sum(1 for r in act if r["classes_seen"] == math.prod(r["G"]))
    info(f"fields where every class of G was seen: {trans} of {len(act)} (note 20: 13, plus q=31 e1=2 generated)")
    check("fields with every class seen", trans, 13)
    pstats = [(r["q"], r["e1"], min(r["primes"]), max(r["primes"])) for r in act]
    info(f"primes l used: from {min(x[2] for x in pstats)} to {max(x[3] for x in pstats)}")

    # ---------------------------------------------------------------- H9
    head("H9", "Levels: the fourfold separates the surfaces as the lattice levels predict (note 23 §4; weight: note 25 finding 10)")
    lev = load("sixlines_jacobian_levels.jsonl")
    okl = [r for r in lev if r["status"] == "ok"]
    by = collections.defaultdict(set)
    for r in okl:
        by[(tuple(r["e"]), r["surface"])].add(r["blocks"])
    incons = [k for k, v in by.items() if len(v) > 1]
    check("pairs computed / degenerate", (len(okl), len(lev) - len(okl)), (422, 2))
    check("surfaces / surfaces whose pairs disagree", (len(by), len(incons)), (29, 0))
    cfgs = {(tuple(r["e"]), tuple(r["cfg"])) for r in okl}
    check("configurations", len(cfgs), 53)
    counts = collections.defaultdict(collections.Counter)
    for (e, s_), v in by.items():
        counts[e][min(v)] += 1
    out = {e: (c[2], c[1]) for e, c in counts.items()}
    info(f"surfaces with two blocks, one block, by polynomial: {out}")
    check("two/one blocks per polynomial (as a sorted list)", sorted(out.values()), sorted([(1, 6), (5, 14), (1, 2)]))
    # self-association from the configuration files
    assoc = {}
    for p in (29, 37):
        for r in load(f"sixlines_p{p}.jsonl"):
            assoc[(p, tuple(r["cfg"]))] = tuple(r["assoc"])
    selfa = collections.Counter()
    for r in okl:
        if r["blocks"] == 2 and assoc.get((r["p"], tuple(r["cfg"]))) == tuple(r["cfg"]):
            selfa[(tuple(r["e"]), r["surface"])] += 1
    sa_by = collections.Counter(e for (e, s_) in selfa)
    info(f"two-block surfaces that are self-associated, by polynomial: {dict(sa_by)} (note 25: 1, 3 and 1)")
    check("self-associated two-block surfaces, total", sum(sa_by.values()), 5)
    info("what is left after note 25: one polynomial (q=37, e1=34) where 2 of 16 non-self-associated surfaces "
         "have two blocks against 2 lattice classes with 2 automorphisms; that all pairs of a surface agree "
         "(422 pairs) is the robust part. Lattice-class counts are in sixlines_jacobian_lattice.log, not re-read here.")


# =================================================================== group b

def group_b():
    print("\n\n######## GROUP G2: size and structure of the group ########")
    gs = load("sixlines_groupsize.jsonl")
    first_block = [13, 17, 19, 23, 25, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]
    sample_q = [211, 503, 1009, 2003, 5003, 10007, 30011, 100003]
    groups = [r for r in gs if r["maximal"] and "G" in r and
              ((r["source"] == "all" and r["q"] in first_block) or (r["source"] == "sample" and r["q"] in sample_q))]

    head("H10", "|G| ~ 0.01 q^2 for maximal orders, q = 13 .. 100,003 (note 13 §3.3)")
    info(f"files: sixlines_groupsize.jsonl; population: {len(groups)} groups, {qrange([r['q'] for r in groups])}; "
         "excluded as note 13 §8 says: q = 125, 343, 625, 2,401, 3,125 and the 'real' duplicates")
    check("number of groups", len(groups), 512)
    ratios = [r["G"] / r["q"] ** 2 for r in groups]
    q1, med, q3 = quantiles(ratios)
    check_close("median |G|/q^2", med, 0.0095, 0.00005)
    check_close("lower quartile", q1, 0.0068, 0.00005)
    check_close("upper quartile", q3, 0.0129, 0.00005)
    for lab, pop in (("all 636 maximal records (with the 24 'real' duplicates and q = 125 .. 3,125)",
                      [r for r in gs if r["maximal"] and "G" in r]),
                     ("612 records without the 'real' duplicates", [r for r in gs if r["maximal"] and "G" in r and r["source"] != "real"])):
        rr_ = [r["G"] / r["q"] ** 2 for r in pop]
        a_, b_, c_ = quantiles(rr_)
        info(f"same statistics over {lab}: quartiles {a_:.4f}, median {b_:.4f}, {c_:.4f}; range {min(rr_):.4f} .. {max(rr_):.4f}")
    check_close("min |G|/q^2", min(ratios), 0.0003, 0.00005)
    check_close("max |G|/q^2", max(ratios), 0.024, 0.0005)
    g41 = [r for r in groups if r["q"] >= 41]
    b, a = lsq_slope([math.log(r["q"]) for r in g41], [math.log(r["G"]) for r in g41])
    check("groups with q >= 41", len(g41), 486)
    check_close("fit exponent, q >= 41", b, 2.01, 0.005)
    g211 = [r for r in groups if r["q"] >= 211]
    b2, a2 = lsq_slope([math.log(r["q"]) for r in g211], [math.log(r["G"]) for r in g211])
    check("groups with q >= 211", len(g211), 242)
    check_close("fit exponent, q >= 211", b2, 1.98, 0.005)
    check_close("fit constant, q >= 211", math.exp(a2), 0.011, 0.0005)
    info("per q (polynomials, maximal, |G| min/median/max, median |G|/q^2) against note 13 §3.3:")
    table = {13: (1, 1, 1, 1, 1), 17: (2, 2, 1, 2, 2), 19: (2, 2, 2, 2, 2), 23: (4, 3, 2, 4, 6), 25: (4, 2, 2, 5, 8),
             29: (9, 3, 4, 6, 12), 31: (10, 6, 2, 8, 10), 37: (16, 7, 3, 10, 26), 41: (21, 6, 4, 12, 26),
             43: (24, 10, 2, 16, 28), 47: (28, 10, 4, 20, 34), 49: (28, 10, 4, 26, 36), 53: (41, 11, 12, 24, 44),
             59: (52, 13, 10, 30, 52), 61: (56, 14, 8, 35, 56), 67: (70, 18, 16, 38, 84), 71: (82, 18, 10, 48, 92),
             73: (89, 17, 6, 48, 82), 79: (104, 23, 12, 56, 96), 83: (118, 19, 18, 64, 116), 89: (142, 22, 12, 74, 148),
             97: (173, 28, 20, 80, 184), 101: (188, 25, 28, 100, 246),
             211: (None, 40, 12, 406, 968), 503: (None, 40, 556, 2845, 5808), 1009: (None, 40, 1668, 11658, 19360),
             2003: (None, 30, 10256, 38498, 97502), 5003: (None, 30, 55176, 255077, 546516),
             10007: (None, 30, 104942, 872620, 2367144), 30011: (None, 20, 1410472, 9762270, 16730664),
             100003: (None, 12, 19135960, 83306780, 173629952)}
    nbad = 0
    for q, (npol, nmax, lo, md, hi) in table.items():
        src = "all" if q in first_block else "sample"
        allq = [r for r in gs if r["q"] == q and r["source"] == src]
        mx = [r["G"] for r in allq if r["maximal"] and "G" in r]
        md_ = statistics.median(mx)
        md_ = int(md_) if md_ == int(md_) else md_
        got = (len(allq) if npol is not None else None, len(mx), min(mx), md_, max(mx))
        want = (npol, nmax, lo, md, hi)
        if got != want:
            nbad += 1
            check(f"table row q={q}", got, want)
    check("rows of the table of note 13 §3.3 that differ from the file", nbad, 0)

    head("H11", "Structure of G: one ramified prime, above 2, norm 64; parity, cyclicity (note 13 §3.4)")
    check("fields with exactly one prime of E ramified over E0 / relative discriminant of norm 64",
          (sum(r["ramified_primes"] == 1 for r in groups), sum(r["ramification"] == 64 for r in groups)), (512, 512))
    check("|G| even", sum(r["G"] % 2 == 0 for r in groups), 509)
    ev = collections.Counter(sum(1 for c in r["G_cyc"] if c % 2 == 0) for r in groups)
    check("number of even cyclic factors 1, 2, 3, 4, 0", tuple(ev[i] for i in (1, 2, 3, 4, 0)), (297, 176, 34, 2, 3))
    oddcyc = sum(1 for r in groups if sum(1 for c in r["G_cyc"] if odd_part(c) > 1) <= 1
                 and all(math.gcd(odd_part(a_), odd_part(b_)) == 1 for i, a_ in enumerate(r["G_cyc"]) for b_ in r["G_cyc"][i + 1:]))
    check("odd part cyclic", oddcyc, 506)
    check("G cyclic", sum(len(r["G_cyc"]) <= 1 for r in groups), 296)
    big = [r for r in groups if r["q"] >= 1000 and r["G"] > 1]
    ex = statistics.median(math.log(largest_prime_factor(r["G"])) / math.log(r["G"]) for r in big)
    check_close("median log(largest prime factor)/log|G|, q >= 1000", ex, 0.55, 0.01)
    # note 18 finding 4
    hk = [r for r in groups if "h0" in r]
    even_h0 = [r for r in hk if r["h0"] % 2 == 0]
    odd_h0 = [r for r in hk if r["h0"] % 2 == 1]
    a1 = sum(1 for r in groups if ev_count(r) > 1)
    a2 = sum(1 for r in even_h0 if ev_count(r) > 1)
    b1 = sum(1 for r in odd_h0 if ev_count(r) <= 1)
    info(f"note 18 finding 4: G with more than one even cyclic factor: {a1}, of which h(E0) even: {a2}; "
         f"h(E0) odd: {len(odd_h0)}, of which at most one even factor: {b1}; h(E0) even in all: {len(even_h0)}")
    check("note 18: (more than one even factor, all with h(E0) even)", (a1, a2), (212, 212))
    check("note 18: (h(E0) odd, all with at most one even factor)", (len(odd_h0), b1), (258, 258),
          f"the two populations quoted add to 470 of 512; h(E0) even with at most one even factor: "
          f"{len(even_h0) - a2}")

    head("H12", "Which polynomials occur: admissible = enumerated (note 13 §3.2), now also at q = 41 (note 26)")
    cm = load("sixlines_cm.jsonl")
    agree = []
    for q in (13, 17, 19, 23, 29, 31, 37):
        adm = {canon_e(r["e"]) for r in gs if r["q"] == q and r["source"] == "all"}
        enu = {canon_e(r["e"]) for r in cm if r["p"] == q}
        agree.append((q, len(adm), len(enu), adm == enu))
    info(f"(q, admissible, enumerated, equal): {agree}")
    check("primes where admissible = enumerated / polynomials", (sum(x[3] for x in agree), sum(x[2] for x in agree)), (7, 44))
    cm41 = load("sixlines_cm_p41.jsonl")
    adm41 = {canon_e(r["e"]) for r in gs if r["q"] == 41 and r["source"] == "all"}
    enu41 = {canon_e(r["e"]) for r in cm41}
    info(f"q = 41 (out of sample: the 'all' run of note 13 predates the enumeration of note 26): admissible {len(adm41)}, "
         f"enumerated {len(enu41)}, equal {adm41 == enu41}")
    check("q = 41: admissible set equals enumerated set", adm41 == enu41, True, "the prediction recorded in note 27, H12")
    m41 = {canon_e(r["e"]): r for r in gs if r["q"] == 41 and r["source"] == "all"}
    lat = [(r["N"][0], m41[canon_e(r["e"])].get("G")) for r in cm41 if r["verdict"] == "match" and canon_e(r["e"]) in m41]
    check("q = 41, maximal orders: lattices per type = |G| of note 13", sum(a_ == b_ for a_, b_ in lat), len(lat),
          f"{lat}")

    head("H13", "The order is maximal only when mu = nu = 1 (note 13 §3.2, q <= 101)")
    allr = [r for r in gs if r["source"] == "all"]
    m11 = [r for r in allr if tuple(r["mu_nu"]) == (1, 1)]
    mx = [r for r in allr if r["maximal"]]
    check("admissible / mu = nu = 1 / maximal", (len(allr), len(m11), len(mx)), (1264, 309, 270))
    check("maximal with (mu, nu) != (1, 1)", sum(1 for r in mx if tuple(r["mu_nu"]) != (1, 1)), 0)
    info(f"{qrange([r['q'] for r in allr])}; above 101 only mu = nu = 1 was sampled, so the observation is untested there")

    head("H14", "One genus: the number of lattices equals |G| (note 13 §3.3, note 19 Lemma 2)")
    withN = [r for r in gs if "N" in r and "G" in r and r["maximal"] and r["q"] not in (125, 343, 625, 2401, 3125)]
    eq = [r for r in withN if all(n_ == r["G"] for n_ in r["N"])]
    by_src = collections.Counter(r["source"] for r in withN)
    info(f"records with a lattice count and a group: {len(withN)} {dict(by_src)}")
    check("records where every type has |G| lattices", (len(eq), len(withN)), (len(withN), len(withN)))
    distinct = {(r["q"], tuple(r["e"])) for r in withN}
    check("distinct polynomials with lattices counted", len(distinct), 191, "note 13 §3.3; §8 says 190 cases for 'all'")
    check("'all' records with lattices counted", by_src["all"], 190)

    head("H15", "Characters of order 2 of G (note 19 §2.5)")
    gen = load("sixlines_characters_genus.jsonl")
    rk = collections.Counter(r["rank2"] for r in gen)
    check("polynomials / q <= 10,007", (len(gen), max(r["q"] for r in gen)), (480, 10007))
    check("2-rank 0..4", tuple(rk[i] for i in range(5)), (3, 287, 160, 28, 2))
    chars = [c for r in gen for c in r["characters"]]
    check("characters, unit, not unit", (len(chars), sum(c["level"] == "unit" for c in chars),
                                          sum(c["level"] != "unit" for c in chars)), (699, 437, 262))
    onesign = sum(1 for c in chars if c["negative_places"] in (0, 3))
    check("d of one sign at the three real places / mixed", (onesign, len(chars) - onesign), (184, 515))
    check("2-rank >= 2 with h(E0) odd / h(E0) odd",
          (sum(1 for r in gen if r["rank2"] >= 2 and r["h0"] % 2 == 1),
           sum(1 for r in gen if r["h0"] % 2 == 1)), (0, 248))
    check("polynomials with a character", sum(1 for r in gen if r["rank2"] >= 1), 477)

    head("H16", "Surfaces per polynomial over all orders ~ 0.006 q^2 (note 13 §3.3, last column)")
    per = collections.defaultdict(list)
    for r in cm:
        per[r["p"]].append(r["k3"])
    quoted = {13: 1.00, 17: 1.50, 19: 2.00, 23: 3.50, 29: 5.00, 31: 5.80, 37: 9.56}
    for p, v in quoted.items():
        check_close(f"average surfaces per polynomial, p={p}", statistics.fmean(per[p]), v, 0.006)
    allk = [(p, statistics.fmean(v) / p**2) for p, v in per.items()]
    info(f"average / q^2 by p: {[(p, round(x, 4)) for p, x in allk]}")


def ev_count(r):
    return sum(1 for c in r["G_cyc"] if c % 2 == 0)


# =================================================================== group c

def group_c():
    print("\n\n######## GROUP G3: computability for the honest party ########")
    head("H17", "The CM method needs precision growing with |G| (notes 20 §9, 24, 26)")
    cmm = load("sixlines_cmmethod.jsonl")
    last = collections.OrderedDict()
    for r in cmm:
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r
    rows = collections.defaultdict(list)
    first_digits = {}
    for r in cmm:
        if r["cmd"] == "compare" and r["status"] == "ok":
            k = (r["q"], r["e"][0])
            first_digits.setdefault(k, r["digits"])
    for (q, e, cmd), r in last.items():
        if cmd == "compare" and r["status"] == "ok":
            rows[r["enumerated_points"]].append(first_digits[(q, e[0])])
    for r in load("sixlines_delta_labels.jsonl"):
        rows[r["N"]].append(r["digits"])
    for r in collections.OrderedDict(((r["q"], r["e1"]), r) for r in load("sixlines_delta41_cosets.jsonl") if "surfaces_in_both" in r).values():
        rows[r["N"]].append(r["digits"])
    for n in sorted(rows):
        info(f"|G| = {n:>2}: first working precision that recognised the polynomials (digits): {sorted(rows[n])}")
    skipped = [(r["q"], r["e"], r["status"]) for r in last.values() if r["status"].startswith("skipped") or "not obtained" in r["status"]]
    info(f"fields skipped or failed in sixlines_cmmethod.jsonl (last record): {skipped}")
    info("largest |G| ever run: 16 (q = 41, e1 = 6, 3,000 digits, 50 min on 8 processes, note 26 §5).")
    info("Extrapolation [H]: the method builds a class polynomial of degree |G|; at |G| ~ 0.01 q^2 ~ 2^256 it cannot run.")

    head("H18", "delta costs about a minute per curve at p <= 41, and about half of the pairs are unusable (notes 23, 25, 26)")
    jc = load("sixlines_jacobian_character.jsonl") + load("sixlines_jacobian_character_extra.jsonl")
    d41 = load("sixlines_delta41.jsonl")
    pr = [r for r in jc if r.get("kind", "primary") in ("primary", "more pairs")]
    sec = [r["seconds"] for r in pr if "seconds" in r]
    sec41 = [r["seconds"] for r in d41 if r["reading"] == "primary" and r["status"] == "ok"]
    info(f"p = 29, 37: seconds per computed curve, median {statistics.median(sec)}, range {min(sec)}..{max(sec)}")
    info(f"p = 41: seconds per computed curve, median {statistics.median(sec41)}, range {min(sec41)}..{max(sec41)}")
    usable = (sum(r["status"] == "ok" for r in pr), len(pr))
    usable41 = (sum(r["status"] == "ok" for r in d41 if r["reading"] == "primary"),
                sum(1 for r in d41 if r["reading"] == "primary"))
    check("p = 29, 37: computed / tried", usable, (106, 242))
    check("p = 41: computed / tried", usable41, (250, 532))
    info("no data above p = 41; cost at cryptographic size is the heuristic of note 24 §5 [H].")


# =================================================================== group d

def surfaces_of(p, cfgs_file):
    """configuration -> canonical surface key (the smaller of cfg and assoc)."""
    out = {}
    for r in load(cfgs_file):
        c, a = tuple(r["cfg"]), tuple(r["assoc"])
        out[(p, c)] = min(c, a)
    return out


def group_d():
    print("\n\n######## GROUP G4: the decision problem ########")

    head("H19", "Square classes (type), branch count: the same on every surface of a torsor (note 19 §4.2)")
    ch = load("sixlines_characters.jsonl")
    gs = load("sixlines_groupsize.jsonl")
    maxi = {}
    for r in gs:
        if r["maximal"] and "G" in r and (r["source"] == "real" or (r["source"] == "all" and r["q"] in (25, 49))):
            maxi[(r["q"], canon_e(r["e"]))] = r["G"]
    surf = {}
    for p in (13, 17, 19, 23, 29, 31, 37):
        surf.update(surfaces_of(p, f"sixlines_p{p}.jsonl"))
    surf.update(surfaces_of(25, "sixlines_q25.jsonl"))
    surf.update(surfaces_of(49, "sixlines_q49.jsonl"))
    tors = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in ch:
        k = (r["q"], canon_e(r["e"]))
        if k in maxi:
            tors[k][surf[(r["q"], tuple(r["cfg"]))]].append(r)
    nsurf = sum(len(v) for v in tors.values())
    check("polynomials with maximal order / surfaces", (len(tors), nsurf), (36, 419))
    check("torsors whose number of surfaces is |G|", sum(len(v) == maxi[k] for k, v in tors.items()), 36)
    types_const = sum(1 for k, v in tors.items() if len({x["type"] for s in v.values() for x in s}) == 1)
    check("torsors on which the type is constant", types_const, 36)
    tmap = collections.defaultdict(set)
    for k, v in tors.items():
        for s in v.values():
            for x in s:
                tmap[(k[0] % 16, x["e"][0] % 16 if x["e"][0] > 0 else (-x["e"][0]) % 16)].add(x["type"])
    info(f"combinations (q mod 16, e1 mod 16) that occur: {len(tmap)}; with more than one type: "
         f"{sum(len(v) > 1 for v in tmap.values())}")
    br_const = sum(1 for k, v in tors.items() if len({x["branch"] for s in v.values() for x in s}) == 1)
    check("torsors on which the branch count is constant", br_const, 36)

    c1m = [r for r in ch if len({tuple(x) for x in r["c1_values"]}) >= 2]
    neg = sum(1 for r in c1m if any([-a, -b] in r["c1_values"] for a, b in r["c1_values"]))
    conj = sum(1 for r in c1m if any([a, -b] in r["c1_values"] and b != 0 for a, b in r["c1_values"]))
    check("note 19 §4.1: configurations where c1 takes two or more values / c1 and -c1 / c1 and its conjugate",
          (len(c1m), neg, conj), (263, 137, 131),
          "the file stores the set of values, not the number of pairs, so the 343 cannot be recounted")

    head("H20", "Conic, j2, j4, quartic product: neither constant nor balanced (note 19 §4.3)")
    res = {}
    for cand, cond in (("conic", lambda q: q % 4 == 1), ("j2", lambda q: True), ("j4", lambda q: True),
                       ("quartic", lambda q: q % 4 == 1)):
        c = collections.Counter()
        for (q, e), v in tors.items():
            if len(v) < 2 or not cond(q):
                continue
            vals = []
            for s, recs in v.items():
                sv = {x[cand] for x in recs if x.get(cand) is not None}
                vals.append(next(iter(sv)) if len(sv) == 1 else None)
            plus, minus = vals.count(1), vals.count(-1)
            if plus == len(vals) or minus == len(vals):
                c["constant"] += 1
            elif plus == minus and plus + minus == len(vals):
                c["balanced"] += 1
            else:
                c["neither"] += 1
            c["torsors"] += 1
        res[cand] = (c["torsors"], c["constant"], c["balanced"], c["neither"])
    check("conic (torsors, constant, balanced, neither)", res["conic"], (23, 3, 5, 15))
    check("j2", res["j2"], (34, 4, 11, 19))
    check("j4", res["j4"], (34, 5, 5, 24))
    check("quartic", res["quartic"], (23, 4, 5, 14))
    info("a value 0 of a power sum is counted as 'neither'; note 19 says a power sum is 0 on a few surfaces")

    head("H21", "J2 J6 and J3 J5 square classes are not characters of the CM labels (note 20 §7.6)")
    cmm = load("sixlines_cmmethod.jsonl")
    last = collections.OrderedDict()
    for r in cmm:
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r
    act = [r for k, r in last.items() if k[2] == "action" and r["status"] == "ok"]
    out = collections.defaultdict(collections.Counter)
    for r in act:
        q = r["q"]
        G = r["G"]
        if r["labelled"] != r["N"] or r["N"] < 4 or math.prod(G) % 2 or r.get("chart") == "J3":
            continue
        surfs = r["surfaces"]
        if any(s["invariants"][0] is None for s in surfs):
            continue
        for name, idx in (("J2J6 (via J6/J2^3)", 1), ("J3J5 (via J3J5/J2^4)", 3)):
            vals = []
            for s in surfs:
                x = s["invariants"][idx] % q
                vals.append(0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1))
            labels = [tuple(s["label"]) for s in surfs]
            verdict = classify_character(vals, labels, G)
            out[name][verdict] += 1
            out[name]["torsors"] += 1
    for name, c in out.items():
        info(f"{name}: {dict(c)}")
    if "J2J6 (via J6/J2^3)" in out:
        c = out["J2J6 (via J6/J2^3)"]
        check("J2J6 (torsors, not a character, a character, constant, zero)",
              (c["torsors"], c["neither"], c["character"], c["constant"], c["zero"]), (8, 6, 1, 0, 1),
              "J6/J2^3 has the square class of J2 J6; J2 itself is not stored, so that row is not recomputed")
        c = out["J3J5 (via J3J5/J2^4)"]
        check("J3J5 (torsors, not a character, a character, constant, zero)",
              (c["torsors"], c["neither"], c["character"], c["constant"], c["zero"]), (8, 5, 0, 1, 2))

    head("H22", "chi_{-1} is a character of G exactly for 80 of 477 polynomials, never for q = 3 mod 4 (note 19 §2.4, §2.5)")
    gen = load("sixlines_characters_genus.jsonl")
    wc = [r for r in gen if r["rank2"] >= 1]
    m1 = [r for r in wc if r["rational"]["-1"]]
    check("d = -1 gives a character / with a character", (len(m1), len(wc)), (80, 477))
    c8 = collections.Counter(r["q"] % 8 for r in wc)
    m8 = collections.Counter(r["q"] % 8 for r in m1)
    check("q = 1 mod 8: (with -1, with a character)", (m8[1], c8[1]), (53, 126))
    check("q = 5 mod 8", (m8[5], c8[5]), (27, 59))
    check("q = 3 mod 4", (m8[3] + m8[7], c8[3] + c8[7]), (0, 292))
    check("d = 2 or -2 gives a character", sum(1 for r in wc if r["rational"]["2"] or r["rational"]["-2"]), 0)
    dl = read_log("sixlines_characters_dyadic.log")
    m = re.search(r"the two conditions agree in (\d+) of (\d+)", dl)
    check("dyadic criterion agrees (sixlines_characters_dyadic.log)", (int(m.group(1)), int(m.group(2))), (477, 477))
    check("all characters trivial on Psi(P): polynomials, primes",
          (sum(r["psi_in_kernel_of_all_characters"] for r in wc), sum(r["psi_tested"] for r in wc)), (477, 19564))
    fq = [r for r in gen if r["q"] in (25, 49)]
    triv = sum(r["frobenius_trivial"] for r in fq)
    minus = sum(1 for r in fq if not r["frobenius_trivial"] and any(v for v in r["frobenius_on_characters"]))
    check("F25, F49: polynomials, Frobenius trivial, some character -1 on it", (len(fq), triv, minus), (12, 2, 9))
    info(f"{qrange([r['q'] for r in gen])}; share of polynomials with chi_-1: {len(m1)}/{len(wc)} = {len(m1)/len(wc):.3f}")

    head("H23", "delta is an invariant of the surface and balanced on each torsor (notes 23 §5.4, 26 §6)")
    jc = load("sixlines_jacobian_character.jsonl")
    jx = load("sixlines_jacobian_character_extra.jsonl")
    prim = jc + [r for r in jx if r["kind"] == "more pairs"]
    st = collections.Counter(r["status"] for r in prim)
    check("p=29,37 curves tried / ok / exponent above 16 / degenerate / not verified",
          (len(prim), st["ok"], st["exponent of (1 - J) G larger than N"], st["degenerate"], st["not verified"]),
          (242, 106, 126, 9, 1))
    ds = collections.defaultdict(list)
    us = collections.defaultdict(list)
    for r in prim:
        if r["status"] == "ok":
            ds[(r["p"], r["e"][0], r["surface"])].append(r["delta"])
            us[(r["p"], r["e"][0], r["surface"])].append(r["u"])
    const = sum(len(set(v)) == 1 for v in ds.values())
    check("surfaces / delta constant on all pairs", (len(ds), const), (32, 32))
    check("u alone constant on all pairs", sum(len(set(v)) == 1 for v in us.values()), 10)
    fr = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in jc:
        if r["status"] == "ok":
            c1, c2 = r["c"][0], r["c"][1]
            s1 = 1 if c1[0] * c1[1] > 0 else -1
            s2 = 1 if c2[0] > 0 else -1
            k = (r["p"], r["e"][0], r["surface"])
            uu = 1 if r["u"] % 4 == 1 else -1
            fr["u"][k].append(uu)
            fr["u*sgn(Re c1 Im c1)"][k].append(uu * s1)
            fr["u*sgn(c2)"][k].append(uu * s2)
            fr["delta"][k].append(r["delta"])
    cnt = {k: sum(len(set(v)) == 1 for v in d.values()) for k, d in fr.items()}
    info(f"first run (sixlines_jacobian_character.jsonl alone): surfaces {len(fr['delta'])}, constant: {cnt}")
    check("note 23 §5.4: first run, constant on (u, u*sgn(Re c1 Im c1), u*sgn(c2), delta) of 29 surfaces",
          (len(fr["delta"]), cnt["u"], cnt["u*sgn(Re c1 Im c1)"], cnt["u*sgn(c2)"], cnt["delta"]), (29, 11, 12, 15, 29),
          "§1 item 6 of note 23 quotes 'on 12 and 15 of the first 29' for the two partial normalisations")
    bal = collections.defaultdict(collections.Counter)
    for (p, e1, s), v in ds.items():
        bal[(p, e1)][v[0]] += 1
    info(f"delta +1 / -1 by torsor: {[(k, (c[1], c[-1])) for k, c in sorted(bal.items())]}")
    check("balanced on each of the three torsors", all(c[1] == c[-1] for c in bal.values()), True)
    npairs = sum(len(v) for v in ds.values())
    pr_agree = Fraction(1, 2 ** (npairs - len(ds)))
    check("chance that all pairs agree if independent coin tosses: exponent of 2", npairs - len(ds), 74)
    pb = balanced_prob([12, 12, 8])
    check_close("chance that a random function is balanced on 12, 12, 8 (note 23: 1.4 %)", float(pb), 0.014, 0.0005)
    ctrl = [r for r in jx if r["kind"] != "more pairs"]
    same = 0
    tot = 0
    prim_u = {(r["p"], tuple(r["cfg"]), tuple(r["pair"])): r["u"] for r in prim if r["status"] == "ok"}
    for r in ctrl:
        if r["status"] != "ok":
            continue
        k = (r["p"], tuple(r["cfg"]), tuple(r["pair"]))
        if k in prim_u:
            tot += 1
            same += prim_u[k] == r["u"]
    check("controls (twists 1, 2, other classes): same u as the primary curve", (same, tot), (36, 36))
    # p = 41
    d41 = load("sixlines_delta41.jsonl")
    p41 = [r for r in d41 if r["reading"] == "primary"]
    st = collections.Counter(r["status"] for r in p41)
    check("p=41 primary: tried / ok / exponent 64 / no vector / not verified / degenerate",
          (len(p41), st["ok"], st["exponent of (1 - J) G larger than N"], st["no kernel vector that gives u"],
           st["not verified"], st["degenerate"]), (532, 250, 269, 8, 3, 2))
    ds41 = collections.defaultdict(list)
    us41 = collections.defaultdict(list)
    for r in p41:
        if r["status"] == "ok":
            ds41[(r["e"][0], r["surface"])].append(r["delta"])
            us41[(r["e"][0], r["surface"])].append(r["u"])
    check("p=41 surfaces / delta constant / u constant",
          (len(ds41), sum(len(set(v)) == 1 for v in ds41.values()), sum(len(set(v)) == 1 for v in us41.values())),
          (28, 28, 0))
    b41 = collections.defaultdict(collections.Counter)
    for (e1, s), v in ds41.items():
        b41[e1][v[0]] += 1
    info(f"p=41 delta +1 / -1 by torsor: {[(k, (c[1], c[-1])) for k, c in sorted(b41.items())]}")
    check("p=41 balanced on both torsors", all(c[1] == c[-1] for c in b41.values()), True)
    np41 = sum(len(v) for v in ds41.values())
    info(f"p=41: pairs with a delta {np41}; chance that all pairs agree if independent: 2^-{np41 - len(ds41)}")
    ctrl41 = collections.Counter()
    prim41 = {(tuple(r["cfg"]), tuple(r["pair"])): r["u"] for r in p41 if r["status"] == "ok"}
    for r in d41:
        if r["reading"] in ("twist1", "twist2", "seed", "r2") and r["status"] == "ok":
            k = (tuple(r["cfg"]), tuple(r["pair"]))
            ctrl41[(r["reading"], "tested")] += k in prim41
            ctrl41[(r["reading"], "same")] += prim41.get(k) == r["u"]
    info(f"p=41 controls: {dict(ctrl41)}")
    check("p=41 controls twist1, twist2, seed: same u / tested",
          (sum(ctrl41[(x, 'same')] for x in ('twist1', 'twist2', 'seed')),
           sum(ctrl41[(x, 'tested')] for x in ('twist1', 'twist2', 'seed'))), (83, 83))

    head("H24", "delta is chi_{-1}: constant on the cosets of ker chi_{-1}, opposite between them (notes 24, 26)")
    dl = load("sixlines_delta_labels.jsonl")
    sizes = []
    agree = tot = 0
    for r in dl:
        by = collections.defaultdict(set)
        for s in r["surfaces"]:
            for d in s["delta"]:
                by[s["orbit"]].add(d)
        ok_ = len(by) == 2 and all(len(v) == 1 for v in by.values()) and by[0] != by[1]
        n_ = len(r["surfaces"])
        agree += n_ if ok_ else 0
        tot += n_
        sizes.append(n_)
        info(f"q={r['q']} e1={r['e1']}: surfaces {n_}, ideals' orbits {r['orbits'][0]}, delta by orbit "
             f"{ {k: sorted(v) for k, v in by.items()} }, primes {len(r['primes'])}")
    check("notes 23, 24: surfaces that agree / surfaces", (agree, tot), (32, 32))
    comp = load_json("sixlines_delta41_compare.json")
    agree41 = 0
    for e1, d in comp.items():
        for route in ("orbits of ideals", "factors over E(i)", "xi modulo 4"):
            info(f"q=41 e1={e1} {route}: agree {d[route]['agree']} of {d['surfaces']}, holds {d[route]['holds']}")
        agree41 += d["orbits of ideals"]["agree"]
        sizes.append(d["surfaces"])
    check("note 26: surfaces that agree (orbits of ideals)", agree41, 28)
    c91 = balanced_chance([12, 12, 8])
    c97 = balanced_chance([16, 12])
    check("note 24: chance 1 in", c91.denominator // c91.numerator, 7470540)
    check("note 26: chance 1 in", c97.denominator // c97.numerator, 2972970)
    call = balanced_chance(sizes)
    info(f"all five torsors {sizes}: chance for a balanced function 1 in {call.denominator // call.numerator:,} "
         f"(~{float(call):.2e}); for G = Z/4 x Z/2 note 24 §6 puts it at 1 in 3 against the two other characters")
    info("population: 60 surfaces, 5 torsors, q in {29, 37, 41}; 5 of the 80 polynomials with chi_-1 in note 19's list "
         "(q = 41 is inside the range of that list)")

    head("H25", "The cosets by three routes, and Frobenius at l = 3 mod 4 (note 25 §4.3, §4.4; note 26 §5)")
    rv = [r for r in load("sixlines_review2.jsonl") if "q" in r]
    for r in rv:
        info(f"q={r['q']} e1={r['e1']}: factors over E {r['factors_over_E']}, over E(i) {r['factors_over_E_i']}, "
             f"= orbits {r['equals_stored_orbits']}, = delta {r['equals_delta_split']}, form = delta "
             f"{r['form_split_equals_delta_split']}, cycles {r['frobenius_cycles_equal_order']}")
    check("three routes agree on the three torsors of note 24",
          sum(r["equals_stored_orbits"] and r["equals_delta_split"] and r["form_split_equals_delta_split"] for r in rv), 3)
    g41 = [r["galois"] for r in load("sixlines_delta41_cosets.jsonl") if "galois" in r]
    for g in g41:
        info(f"q=41 e1={g['e1']}: factors over E {g['factors_over_E']}, over E(i) {g['factors_over_E_i']}, "
             f"= orbits {g['equals_stored_orbits']}, form = orbits {g['form_split_equals_stored_orbits']}, "
             f"cycles {g['frobenius_cycles_equal_order']}")
    check("q=41: factors and form both equal the orbits", sum(g["equals_stored_orbits"] and g["form_split_equals_stored_orbits"] for g in g41), 2)
    l3 = [sum(v for k, v in r["frobenius_cycles_equal_order"].items() if k.startswith("l = 3")) for r in rv]
    l3_41 = [sum(v for k, v in g["frobenius_cycles_equal_order"].items() if k.startswith("l = 3")) for g in g41]
    nl3 = [sum(v for k, v in r["frobenius_cycles_equal_order"].items() if k.startswith("l = 3") and "yes" not in k) for r in rv + g41]
    check("primes l = 3 mod 4 with cycle lengths = order of the class: note 25, note 26", (sum(l3), sum(l3_41)), (160, 128))
    check("primes l = 3 mod 4 where they differ", sum(nl3), 0)
    check("primes l = 1 mod 4 with cycle lengths = order (note 25: 120)", sum(sum(v for k, v in r["frobenius_cycles_equal_order"].items() if k.startswith("l = 1")) for r in rv), 120)
    d9 = load("sixlines_review2.jsonl")
    dd = [r for r in d9 if "stored_delta" in r]
    check("delta recomputed on nine curves equal to stored", (sum(r["delta"] == r["stored_delta"] for r in dd), len(dd)), (9, 9))

    head("H26", "Sign rule delta = zeta * Gamma/H' on all five torsors (note 26 §7.1, found after the fact)")
    slog = read_log("sixlines_delta41_compare.log")
    for line in slog.splitlines():
        if "zeta" in line or "rule" in line.lower():
            info("log: " + line.strip()[:150])
    info("chance if the five signs were coin tosses: one rule fits with probability 2 * 2^-5 = 1/16; formed after 4 of 5")


def classify_character(vals, labels, G):
    """'zero' if a value is 0, 'constant', 'character' if vals = +-psi(label) for a
    non-trivial character psi of order 2 of G = prod Z/G_i, else 'neither'."""
    if 0 in vals:
        return "zero"
    if len(set(vals)) == 1:
        return "constant"
    import itertools
    evens = [i for i, g in enumerate(G) if g % 2 == 0]
    for bits in itertools.product((0, 1), repeat=len(evens)):
        if not any(bits):
            continue
        def psi(lab):
            s = 0
            for b, i in zip(bits, evens):
                s += b * lab[i]
            return 1 if s % 2 == 0 else -1
        ps = [psi(l) for l in labels]
        if ps == vals or [-x for x in ps] == vals:
            return "character"
    return "neither"


# =================================================================== group e

def group_e():
    print("\n\n######## GROUP G5: below Picard number 16 ########")
    head("H27", "Two Jacobians of a (4,4) double cover: the local structure at 2 that a 2-step needs (note 22 §4)")
    fam = {"e8a2 p5": "ellk3_e8a2_p5.jsonl", "e8a2 p7": "ellk3_e8a2_p7.jsonl", "e8 p5": "ellk3_e8_p5.jsonl"}
    gen = load("ellk3_generic_p5.jsonl") + load("ellk3_generic_p5_seed2.jsonl")
    data = {k: load(v) for k, v in fam.items()}
    data["generic p5"] = gen
    quoted = {"e8a2 p5": (16000, 6867, 2621, 722), "e8a2 p7": (4000, 1801, 777, 201), "e8 p5": (1200, 533, 212, 49),
              "generic p5": (40, 18, 3, 0)}
    tot_iro = 0
    for k, rr in data.items():
        withP = [r for r in rr if "P" in r]
        iro = [r for r in withP if r["status"] in ("ok", "p not maximal")]
        cc = sum(r["cond_C"] for r in withP)
        check(f"{k}: pairs, with a polynomial, irreducible ordinary, condition C", (len(rr), len(withP), len(iro), cc), quoted[k])
        none = sum(1 for r in iro if r["cat2"] in ("", "none"))
        split = sum(1 for r in iro if r["cat2"] == "split")
        multi = sum(1 for r in iro if r["cat2"] == "multi")
        info(f"{k}: irreducible ordinary: (t-1) does not divide P mod 2 exactly twice or more: {none}, exactly twice (split pair): {split}, 4+ times: {multi}")
        tot_iro += len(iro)
        br = [r for r in withP if int(r["br"]) != 0]
        info(f"{k}: Artin-Tate number with P(p) != 0: {len(br)}, a square {sum(r['br_square'] for r in br)}, "
             f"4 times a square {sum(r['br_4square'] for r in br)}")
    check("irreducible ordinary pairs, all families", tot_iro, 3613)
    ctl = []
    for name in ("ellk3_e8a2_p5_control.jsonl", "ellk3_e8_p5_control.jsonl"):
        rr = load(name)
        iro = [r for r in rr if r.get("status") in ("ok", "p not maximal")]
        c = collections.Counter(r["cat2"] for r in iro)
        br = [r for r in rr if "P" in r and int(r["br"]) != 0]
        info(f"control {name}: irreducible ordinary {len(iro)}, cat2 {dict(c)}; P(p) != 0: {len(br)}, 4*square "
             f"{sum(r['br_4square'] for r in br)}")
        ctl.append((len(iro), c["none"], c["split"] + c["inert"], c["inert"], c["multi"], len(br), sum(r["br_4square"] for r in br)))
    check("control e8a2: (irr. ord., t-1 does not divide, exactly twice, of which inert, 4+, P(p)!=0, 4*square)", ctl[0], (345, 269, 31, 9, 45, 608, 214))
    check("control e8", ctl[1], (147, 112, 18, 8, 17, 235, 84))
    pairs_ = [r for k in ("e8a2 p5", "e8a2 p7", "e8 p5", "generic p5") for r in data[k] if r.get("status") in ("ok", "p not maximal")]
    check("pairs: t-1 does not divide P mod 2 / exactly twice and split / inert",
          (sum(r["cat2"] in ("", "none") for r in pairs_), sum(r["cat2"] == "split" for r in pairs_), sum(r["cat2"] == "inert" for r in pairs_)), (0, 1349, 0))
    brs = [r for k in ("e8a2 p5", "e8a2 p7", "e8 p5", "generic p5") for r in data[k] if "P" in r and int(r["br"]) != 0]
    check("pairs with P(p) != 0: Artin-Tate number 4 times a square", (sum(r["br_4square"] for r in brs), len(brs)), (5548 + 1457 + 403 + 7,) * 2)

    head("H28", "Under condition C the two Jacobians are never isomorphic and the class of l/lbar is never trivial (note 22 §4.3)")
    allC = [r for k in ("e8a2 p5", "e8a2 p7", "e8 p5") for r in data[k] if r.get("cond_C")]
    check("condition C pairs / class trivial / X1 = X2", (len(allC), sum(r["c_trivial"] for r in allC), sum(bool(r["iso"]) for r in allC)), (972, 0, 0))
    s = {k: sum(1 / r["G_size"] for r in data[k] if r.get("cond_C")) for k in ("e8a2 p5", "e8a2 p7", "e8 p5")}
    for k, v in s.items():
        info(f"{k}: sum 1/|G| = {v:.2f}")
    check_close("sum 1/|G|, e8a2 p5", s["e8a2 p5"], 21.5, 0.05)
    check_close("sum 1/|G|, all (note 22 says 23.5)", sum(s.values()), 23.5, 0.05)
    pnone = math.prod(1 - 1 / r["G_size"] for r in data["e8a2 p5"] if r.get("cond_C"))
    info(f"probability of no isomorphic pair in e8a2 p5 if X2 were uniform in the torsor: prod(1 - 1/|G|) = {pnone:.2e}; "
         f"Poisson exp(-sum) = {math.exp(-s['e8a2 p5']):.2e}")
    check_close("note 22: 2e-10", pnone, 2e-10, 1e-10)
    for k, expo, qv in (("e8a2 p5", 4, 0.083), ("e8a2 p7", 4, 0.078), ("e8 p5", 5, 0.065)):
        p = 5 if "p5" in k else 7
        md = statistics.median(r["G_size"] / p**expo for r in data[k] if r.get("cond_C"))
        check_close(f"{k}: median |G|/p^{expo}", md, qv, 0.0006)
    glog = read_log("ellk3_e8a2_graph.log")
    m = re.findall(r"number of distinct partners of a surface: (\{.*?\})", glog)
    info(f"partner counts in ellk3_e8a2_graph.log: {m}")
    check("condition C, p=5: partners 1, 2, 3", tuple(eval(m[0]).get(i, 0) for i in (1, 2, 3)), (1108, 86, 0))


# ======================================================== inventory for the audit

def audit_inventory():
    print("\n\n######## INVENTORY FOR THE LEAKAGE AUDIT (note 28) ########")
    head("AUDIT", "Stored CM labels and stored per-configuration invariants that can be joined")
    cmm = load("sixlines_cmmethod.jsonl")
    last = collections.OrderedDict()
    for r in cmm:
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r
    lab = {}
    for (q, e, cmd), r in last.items():
        if cmd == "action" and r["status"] == "ok":
            for s_ in r["surfaces"]:
                if s_.get("label"):
                    for c in s_["configurations"]:
                        lab[(q, tuple(c))] = ("label", tuple(r["G"]), tuple(s_["label"]))
    for r in load("sixlines_delta_labels.jsonl"):
        for s_ in r["surfaces"]:
            for c in s_["configurations"]:
                lab.setdefault((r["q"], tuple(c)), ("coset", r["N"], s_["orbit"]))
    for e1, d in load_json("sixlines_delta41_compare.json").items():
        for s_ in d["table"]:
            for c in s_["configurations"]:
                lab.setdefault((41, tuple(c)), ("coset", d["surfaces"], s_["orbits of ideals"]))
    kinds = collections.Counter(v[0] for v in lab.values())
    info(f"configurations with a stored label in G: {kinds['label']}; with only the coset of ker chi_-1: {kinds['coset']}")
    ch = {(r["q"], tuple(r["cfg"])): r for r in load("sixlines_characters.jsonl")}
    joined = collections.Counter(v[0] for k, v in lab.items() if k in ch)
    info(f"of these, with the invariants of sixlines_characters.jsonl (type, conic, j2, j4, quartic, branch): {dict(joined)}")
    info("q = 41 configurations have labels (cosets) but no record in sixlines_characters.jsonl; note 19's invariants "
         "for them would have to be recomputed from sixlines_p41.jsonl (a rerun, not a read).")
    info("torsors with full labels: those of note 20 §7.5 (67 surfaces in 15 fields); with cosets only: the five of notes 24, 26.")


# ====================================================================== main

def main():
    print("sixlines_ledger.py: hypothesis ledger for six-line K3 surfaces (note 27). Reads stored files only.")
    print("'quoted' is the number as a note first stated it; where a check fails, the note now states the recomputed value.")
    group_a()
    group_b()
    group_c()
    group_d()
    group_e()
    audit_inventory()
    print("\n\n######## FILES USED ########")
    for k, v in FILES_USED.items():
        print(f"    {k}: {v} records/lines")
    print("\n######## MISMATCHES (quoted = as first stated; the notes now state the recomputed values) ########")
    bad = [c for c in CHECKS if not c[4]]
    for hid, what, rec, quo, ok, note in bad:
        print(f"    {hid}: {what}: recomputed {rec!r}, quoted {quo!r}" + (f"  ({note})" if note else ""))
    print(f"\n    {len(CHECKS)} comparisons, {len(CHECKS) - len(bad)} agree, {len(bad)} differ")


if __name__ == "__main__":
    main()
