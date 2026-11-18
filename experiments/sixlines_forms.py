#!/usr/bin/env python
"""
sixlines_forms.py -- non-split forms over F_p of six-line double planes (note 4).

A form is a set of six lines in general position that is stable under the p-Frobenius, with every
line defined over F_{p^2}.  Frobenius then permutes the lines by an involution, of type

  A  (2,1,1,1,1)   four rational lines and one conjugate pair
  B  (2,2,1,1)     two rational lines and two conjugate pairs
  C  (2,2,2)       three conjugate pairs

PGL_3(F_p) acts simply transitively on the frames (four labelled lines in general position) with a
given Frobenius action; `enumerate` checks this by counting the frames.  So every form of a type is
equivalent to a fixed frame of that type plus one conjugate pair, and the enumeration runs over the
non-rational points of the dual plane.

Two labelled forms with the same Frobenius action are equivalent over F_p exactly when they are
equivalent over F_{p^2} (a labelled set of six lines has no automorphisms).  So the class of a form
over F_p is the orbit of its split normal form (b, c, e, f) under the relabellings that commute with
the Frobenius involution, and its class over F_{p^2} is the orbit under all 720.

  enumerate --p P --out FILE    every class of forms, with its class over F_{p^2} and the character
                                sums S_1, S_2, S_3 of w^2 = l_1 ... l_6 over F_p, F_{p^2}, F_{p^3}
  compare FORMS Q_FILE          against the conjugation-stable classes of sixlines_fq.py, and the
                                Frobenius polynomial over F_p against the one over F_{p^2}
  selftest

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sixlines import Field, char_sum as split_char_sum  # noqa: E402
from sixlines_fq import FieldExt, Tables, Tower, decode  # noqa: E402

FRAME = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]
INVOLUTION = {"A": (0, 1, 2, 3, 5, 4), "B": (0, 1, 3, 2, 5, 4), "C": (1, 0, 3, 2, 5, 4)}
ALL_PERMS = list(itertools.permutations(range(6)))


def centraliser(pi):
    return [s for s in ALL_PERMS if all(s[pi[i]] == pi[s[i]] for i in range(6))]


# ---------------------------------------------------------------- geometry over F_{p^2}

def projective_points(q):
    """Normalised representatives of P^2(F_q), as an (n, 3) array of codes."""
    v = np.arange(q, dtype=np.int32)
    y, z = (g.ravel() for g in np.meshgrid(v, v, indexing="ij"))
    return np.concatenate([np.stack([np.ones_like(y), y, z], axis=1),
                           np.stack([np.zeros(q, dtype=np.int32), np.ones(q, dtype=np.int32), v], axis=1),
                           np.array([[0, 0, 1]], dtype=np.int32)])


def det3(T, u, v, w):
    return T.dot(T.cross(u, v), w)


def general_position(T, pts):
    ok = np.ones(np.broadcast(*[c for pt in pts for c in pt]).shape, dtype=bool)
    for i, j, k in itertools.combinations(range(len(pts)), 3):
        ok &= det3(T, pts[i], pts[j], pts[k]) != 0
    return ok


def normal_form(T, pts):
    """(b, c, e, f) of six labelled points in general position: the first four go to the frame."""
    cr = [T.cross(pts[1], pts[2]), T.cross(pts[2], pts[0]), T.cross(pts[0], pts[1])]
    lam = [T.dot(cr[i], pts[3]) for i in range(3)]
    mu = [T.mul(lam[1], lam[2]), T.mul(lam[0], lam[2]), T.mul(lam[0], lam[1])]
    out = []
    for j in (4, 5):
        x = [T.mul(mu[i], T.dot(cr[i], pts[j])) for i in range(3)]
        s = T.INV[x[0]]
        out += [T.mul(x[1], s), T.mul(x[2], s)]
    return np.stack(out, axis=1).astype(np.int32)


def const(pt, n):
    return [np.full(n, c, dtype=np.int32) for c in pt]


def cols(arr):
    return [arr[:, 0], arr[:, 1], arr[:, 2]]


def standard_frames(T, p):
    """A frame of each type, and the number of frames of each type (to compare with |PGL_3(F_p)|)."""
    pts = projective_points(T.q)
    rational = pts[(pts < p).all(axis=1)]
    nonrat = pts[~(pts < p).all(axis=1)]
    conj = T.FROB[nonrat]
    counts, frames = {}, {}
    # B: two rational points and a conjugate pair
    total, first = 0, None
    for r1, r2 in itertools.permutations(range(len(rational)), 2):
        n = len(nonrat)
        ok = general_position(T, [const(rational[r1], n), const(rational[r2], n), cols(nonrat), cols(conj)])
        total += int(ok.sum())
        if first is None and ok.any() and tuple(rational[r1]) == FRAME[0] and tuple(rational[r2]) == FRAME[1]:
            first = int(np.argmax(ok))
    counts["B"] = total
    frames["B"] = [FRAME[0], FRAME[1], tuple(int(v) for v in nonrat[first]), tuple(int(v) for v in conj[first])]
    # C: two conjugate pairs
    total, first = 0, None
    n = len(nonrat)
    for i in range(n):
        ok = general_position(T, [const(nonrat[i], n), const(conj[i], n), cols(nonrat), cols(conj)])
        total += int(ok.sum())
        if first is None and ok.any():
            first = (i, int(np.argmax(ok)))
    counts["C"] = total
    i, j = first
    frames["C"] = [tuple(int(v) for v in x) for x in (nonrat[i], conj[i], nonrat[j], conj[j])]
    frames["A"] = list(FRAME)
    return frames, counts, nonrat, conj


# ---------------------------------------------------------------- character sums

class Big:
    """F_{p^n} with vectorised arithmetic, containing F_{p^2} (codes of `base`) and F_{p^k} for k | n."""

    def __init__(self, base, n):
        p = self.p = base.p
        self.F = F = base if n == 2 else FieldExt(p, n)
        q = self.q = F.q
        primes = [r for r in range(2, 2000) if (q - 1) % r == 0 and all(r % s for s in range(2, r))]
        g = next(a for a in range(2, q) if all(F.power(a, (q - 1) // r) != 1 for r in primes))
        exp = np.zeros(q - 1, dtype=np.int64)
        x = 1
        for i in range(q - 1):
            exp[i] = x
            x = F.mul(x, g)
        assert x == 1 and len(set(exp.tolist())) == q - 1
        self.exp = exp
        self.log = np.zeros(q, dtype=np.int64)
        self.log[exp] = np.arange(q - 1)
        self.digits = F.digits.astype(np.int64)
        self.pw = p ** np.arange(n, dtype=np.int64)
        if n == 2:
            self.emb2 = np.arange(q, dtype=np.int64)
        else:
            c0, c1 = base.mod[0], base.mod[1]
            beta = next(x for x in range(q) if F.add(F.add(F.mul(x, x), F.scal(c1, x)), c0 % p) == 0)
            self.emb2 = np.array([F.add(code % p, F.scal(code // p, beta)) for code in range(p * p)], dtype=np.int64)

    def add(self, a, b):
        return ((self.digits[a] + self.digits[b]) % self.p) @ self.pw

    def mul(self, a, b):
        a, b = np.broadcast_arrays(np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64))
        out = self.exp[(self.log[a] + self.log[b]) % (self.q - 1)]
        return np.where((a == 0) | (b == 0), 0, out)

    def subfield(self, k):
        m = (self.q - 1) // (self.p ** k - 1)
        return np.concatenate([[0], self.exp[::m]])

    def chi_sub(self, k):
        """Quadratic character of F_{p^k}, as a table on the codes of the big field (0 off the subfield)."""
        m = (self.q - 1) // (self.p ** k - 1)
        chi = np.zeros(self.q, dtype=np.int64)
        idx = np.arange(0, self.q - 1, m)
        chi[self.exp[idx]] = np.where((idx // m) % 2 == 0, 1, -1)
        return chi

    def line_values(self, line, k):
        """Values of a x + b y + c z on the points (1,y,z), (0,1,z), (0,0,1) of P^2(F_{p^k})."""
        a, b, c = (int(v) for v in line)
        sub = self.subfield(k)
        by, cz = self.mul(b, sub), self.mul(c, sub)
        affine = self.add(self.add(by[:, None], cz[None, :]), a).ravel()
        return np.concatenate([affine, self.add(cz, b), [c]])


def form_sums(bigs, T, rational, pairs):
    """S_1, S_2, S_3 for the product of the rational lines and of the conjugate pairs."""
    out = []
    for k in (1, 2, 3):
        big = bigs[k]
        if k % 2:
            chi_k, chi_big = big.chi_sub(k), big.F.chi.astype(np.int64)
            prod = 1
            for ln in rational:
                prod = prod * chi_k[big.line_values(ln, k)]
            for ln in pairs:
                # chi_k(l lbar) = chi_k(norm of l) = character of the quadratic extension at l
                prod = prod * chi_big[big.line_values([big.emb2[v] for v in ln], k)]
        else:
            chi = big.F.chi.astype(np.int64)
            lines = list(rational) + list(pairs) + [tuple(int(T.FROB[v]) for v in ln) for ln in pairs]
            prod = 1
            for ln in lines:
                prod = prod * chi[big.line_values([big.emb2[v] for v in ln], k)]
        out.append(int(np.sum(prod)))
    return out


def split_form(kind, pts):
    """Rational lines and one line of each conjugate pair, from the labelled six."""
    if kind == "A":
        return pts[:4], [pts[4]]
    if kind == "B":
        return pts[:2], [pts[2], pts[4]]
    return [], [pts[0], pts[2], pts[4]]


# ---------------------------------------------------------------- commands

def cmd_enumerate(args):
    p, t0 = args.p, time.time()
    T = Tables(p)
    q = T.q
    frames, counts, nonrat, conj = standard_frames(T, p)
    order = p ** 3 * (p ** 3 - 1) * (p ** 2 - 1)
    counts["A"] = order                                  # the frame itself, by definition
    print(f"p={p}: |PGL_3(F_p)| = {order}; frames of type B: {counts['B']}, type C: {counts['C']} "
          f"({time.time() - t0:.1f}s)", file=sys.stderr)
    assert counts["B"] == order and counts["C"] == order
    bigs = {1: Big(T.base, 2), 3: Big(T.base, 6)}
    bigs[2] = bigs[1]
    rows = []
    for kind in "ABC":
        n = len(nonrat)
        pts = [const(x, n) for x in frames[kind]] + [cols(nonrat), cols(conj)]
        ok = general_position(T, pts)
        idx = np.nonzero(ok)[0]
        pts = [[c[idx] for c in pt] for pt in pts]
        cfg = normal_form(T, pts)
        key_p = T.keys(cfg, centraliser(INVOLUTION[kind]))
        key_q = T.keys(cfg, ALL_PERMS)
        classes = {}
        for i, k in enumerate(key_p.tolist()):
            classes.setdefault(k, []).append(i)
        for k in sorted(classes):
            members = classes[k]
            i = members[0]
            six = [tuple(int(c[i]) for c in pt) for pt in pts]
            assert len({int(key_q[m]) for m in members}) == 1
            rational, pairs = split_form(kind, six)
            rows.append({"p": p, "type": kind, "lines": [list(x) for x in six], "cfg": [int(v) for v in cfg[i]],
                         "candidates": len(members), "class_q": decode(int(key_q[i]), q),
                         "S": form_sums(bigs, T, rational, pairs)})
        print(f"p={p}: type {kind}: {len(idx)} forms on the standard frame, {len(classes)} classes "
              f"({time.time() - t0:.1f}s)", file=sys.stderr)
    with open(args.out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"p={p}: wrote {len(rows)} rows to {args.out}", file=sys.stderr)


def graeffe(co):
    """Monic polynomial whose roots are the squares of the roots of co (descending coefficients)."""
    n = len(co) - 1
    minus = [c * (-1) ** i for i, c in enumerate(co)]
    prod = [0] * (2 * n + 1)
    for i, a in enumerate(co):
        for j, b in enumerate(minus):
            prod[i + j] += a * b
    assert all(v == 0 for v in prod[1::2])
    return prod[0::2]


def twist(co):
    return [c * (-1) ** i for i, c in enumerate(co)]


def cmd_compare(args):
    import sixlines_cm as cm
    forms = [json.loads(line) for line in open(args.forms)]
    rows = [json.loads(line) for line in open(args.qfile)]
    p, q = forms[0]["p"], rows[0]["p"]
    by = {tuple(r["cfg"]): r for r in rows}
    T = Tables(p)
    stable = {tuple(r["cfg"]) for r in rows if r["frob"] == r["cfg"]}
    swapped = {tuple(r["cfg"]) for r in rows if r["frob"] == r["assoc"] != r["cfg"]}
    # split configurations over F_p, base-changed
    split = set()
    if args.split:
        srows = [json.loads(line) for line in open(args.split)]
        cfg = np.array([r["cfg"] for r in srows], dtype=np.int32)
        split = {tuple(decode(int(k), q)) for k in T.keys(cfg, ALL_PERMS)}
    hit = defaultdict(list)
    for f in forms:
        hit[tuple(f["class_q"])].append(f)
    from_forms = set(hit)
    print(f"p={p}: {len(forms)} classes of non-split forms "
          f"({dict(sorted(Counter(f['type'] for f in forms).items()))}), in {len(from_forms)} classes over F_{q}")
    print(f"  conjugation-stable classes over F_{q}: {len(stable)}; from split configurations over F_{p}: "
          f"{len(split)}; from non-split forms: {len(from_forms)}; from both: {len(split & from_forms)}")
    print(f"  stable classes with no form over F_{p}: {len(stable - split - from_forms)}; "
          f"forms whose class is not stable: {len(from_forms - stable)}")
    trivial = {c for c in stable if by[c]["orbit"] == 720}
    multi = Counter(len(hit[c]) for c in trivial)
    print(f"  stable classes without automorphisms: {len(trivial)}; forms per such class: {dict(sorted(multi.items()))}")
    print(f"  classes fixed only up to association: {len(swapped)}")

    # Frobenius polynomials: over F_p from the form, over F_q from the q-file
    stats, per_poly = Counter(), defaultdict(list)
    for f in forms:
        r = by[tuple(f["class_q"])]
        e1, e2, e3, integral = cm.elementary(r["S"])
        signs = cm.choose_sign(q, (e1, e2, e3))
        if not integral or len(signs) != 1:
            stats["F_q polynomial not determined"] += 1
            continue
        big = cm.charpoly(q, (e1, e2, e3), signs[0])
        a1, a2, a3, integral = cm.elementary(f["S"])
        if not integral:
            stats["F_p sums not integral"] += 1
            continue
        found = None
        for eps in (1, -1):
            if eps == -1 and a3 != 0:
                continue
            small = cm.charpoly(p, (a1, a2, a3), eps)
            if cm.weil_ok(p, small) and graeffe(small) in (big, twist(big)):
                found = small
        f["poly_p"], f["poly_q"] = found, big
        stats["roots of the F_p polynomial square to those of the F_q polynomial" if found
              else "MISMATCH between the F_p and F_q polynomials"] += 1
        f["S2_agrees"] = abs(f["S"][1]) == abs(r["S"][0])
        stats["S_2 of the form = +-S_1 of the class" if f["S2_agrees"] else "S_2 MISMATCH"] += 1
        ordinary = e1 % p != 0
        irreducible = bool(cm.pari.polisirreducible(cm.pol(big)))
        if ordinary and irreducible:
            sign = 1 if e1 > 0 else -1
            per_poly[(sign * e1, e2, sign * e3)].append(f)
    for k, v in sorted(stats.items()):
        print(f"  {k}: {v}")
    print(f"  forms with ordinary irreducible polynomial over F_{q}: {sum(len(v) for v in per_poly.values())}")
    out = []
    for key in sorted(per_poly):
        fs = per_poly[key]
        k3 = {tuple(sorted([tuple(f["class_q"]), tuple(by[tuple(f["class_q"])]["assoc"])])) for f in fs}
        types = dict(sorted(Counter(f["type"] for f in fs).items()))
        polys = sorted({tuple(f["poly_p"]) if f["poly_p"] else None for f in fs}, key=str)
        print(f"  e={list(key)}: {len(fs)} forms {types}, {len(k3)} K3 surfaces over F_{q}")
        for f in fs:
            print(f"      type {f['type']} class over F_{q} {f['class_q']} S={f['S']} polynomial over F_{p} {f['poly_p']}")
        out.append({"q": q, "e": list(key), "forms": len(fs), "types": types, "k3": len(k3),
                    "polys_p": [list(x) if x else None for x in polys]})
    if args.out:
        with open(args.out, "w") as fh:
            for rec in out:
                fh.write(json.dumps(rec) + "\n")


def cmd_selftest(args):
    ok = True
    p = 7
    T = Tables(p)
    bigs = {1: Big(T.base, 2), 3: Big(T.base, 6)}
    bigs[2] = bigs[1]
    # 1. a split configuration over F_p: the three sums against sixlines.py
    cfg = (2, 3, 4, 6)
    lines = FRAME + [(1, cfg[0], cfg[1]), (1, cfg[2], cfg[3])]
    mine = form_sums(bigs, T, lines, [])
    ref = [split_char_sum(Field(p, k), cfg) for k in (1, 2, 3)]
    ok &= mine == ref
    print(f"split {cfg}: {mine} against sixlines.py {ref} {'ok' if mine == ref else 'MISMATCH'}")
    # 2. a form of type A: S_2 against the split sum over F_{p^2}; S_1 by direct count over F_p
    b, c = 9, 23
    pair = (1, b, c)
    mine = form_sums(bigs, T, FRAME, [pair])
    ref2 = Tower(T.base, 1).char_sum((b, c, int(T.FROB[b]), int(T.FROB[c])))
    B = T.base
    direct = 0
    for pt in [(1, y, z) for y in range(p) for z in range(p)] + [(0, 1, z) for z in range(p)] + [(0, 0, 1)]:
        val = 1
        for ln in FRAME + [pair, tuple(int(T.FROB[v]) for v in pair)]:
            v = 0
            for coef, coord in zip(ln, pt):
                v = B.add(v, B.mul(coef, coord))
            val = B.mul(val, v)
        assert val < p
        direct += 0 if val == 0 else (1 if pow(val, (p - 1) // 2, p) == 1 else -1)
    ok &= mine[1] == ref2 and mine[0] == direct
    print(f"type A (1,{b},{c}): S_2 {mine[1]} against sixlines_fq.py {ref2}; S_1 {mine[0]} against direct {direct} "
          f"{'ok' if mine[1] == ref2 and mine[0] == direct else 'MISMATCH'}")
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("enumerate")
    e.add_argument("--p", type=int, required=True)
    e.add_argument("--out", required=True)
    e.set_defaults(func=cmd_enumerate)
    c = sub.add_parser("compare")
    c.add_argument("forms")
    c.add_argument("qfile")
    c.add_argument("--split", help="sixlines.py output over F_p, to account for split configurations")
    c.add_argument("--out")
    c.set_defaults(func=cmd_compare)
    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    args = ap.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()
