#!/usr/bin/env python
"""
sixlines_fq.py -- six-line double planes over F_q, q = p^2, with the action of Frobenius conjugation (note 3).

Same family and conventions as sixlines.py, but the last two lines (1,b,c), (1,e,f) have
coefficients in F_q.  Each row carries, besides the class size, the associated class and S_1, S_2, S_3
(character sums over F_q, F_{q^2}, F_{q^3}), the class of the conjugate configuration
(b^p, c^p, e^p, f^p).  Conjugation X -> X^(p) permutes the surfaces over F_q with a given Frobenius
polynomial; the dictionary predicts that it is translation by the class of p/pbar for the unit-root
prime p of E.  sixlines_cm.py --neighbours compares the two.

  enumerate --p P --out FILE    (q = P^2)
  selftest

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sixlines import Field  # noqa: E402


def poly_divides(d, f, p):
    """Does the monic polynomial d divide f over F_p?  (coefficient lists, lowest degree first)"""
    r = list(f)
    for i in range(len(r) - 1, len(d) - 2, -1):
        c = r[i] % p
        if c:
            for j, dj in enumerate(d):
                r[i - len(d) + 1 + j] -= c * dj
    return all(x % p == 0 for x in r[:len(d) - 1])


class FieldExt(Field):
    """F_{p^k} for any k, with a genuinely irreducible modulus."""

    def _irreducible(self):
        p, k = self.p, self.k
        if k == 1:
            return [0, 1]
        divisors = [list(c) + [1] for d in range(1, k // 2 + 1) for c in itertools.product(range(p), repeat=d)]
        for coeffs in itertools.product(range(p), repeat=k):
            poly = list(coeffs) + [1]
            if poly[0] and not any(poly_divides(d, poly, p) for d in divisors):
                return poly
        raise RuntimeError("no irreducible polynomial")


class Tower:
    """F_q inside F_{q^k}: embedding of the codes of F_q and the orbits of x -> x^q."""

    def __init__(self, base, k):
        p = base.p
        self.F = F = FieldExt(p, 2 * k)
        if k == 1:
            beta = p                                  # same modulus search order, so same field
            assert F.mod == base.mod
        else:
            c0, c1 = base.mod[0], base.mod[1]
            beta = next(x for x in range(F.q)
                        if F.add(F.add(F.mul(x, x), F.scal(c1, x)), c0 % p) == 0)
        self.emb = [F.add(d0, F.scal(d1, beta)) for d1 in range(p) for d0 in range(p)]
        self.emb = [self.emb[(code % p) * 1 + (code // p) * p] for code in range(p * p)]
        seen = np.zeros(F.q, dtype=bool)
        reps = []
        for a in range(1, F.q):
            if seen[a]:
                continue
            n, x = 0, a
            while not seen[x]:
                seen[x] = True
                x = F.power(x, p * p)
                n += 1
            reps.append((a, n))
        self.orbits = reps

    def char_sum(self, cfg):
        F, p = self.F, self.F.p
        b, c, e, f = (self.emb[v] for v in cfg)
        X, axes, total = F.X, tuple(range(F.k)), 0
        for y, weight in self.orbits:
            ts = [F.add(y, 1), F.add(F.mul(b, y), c), F.add(F.mul(e, y), f)]
            prod = X
            for t in ts:
                prod = prod * np.roll(X, [-int(d) for d in F.digits[t]], axes)
            total += weight * int(F.chi[y]) * int(prod.sum(dtype=np.int64))
        return total


class Tables:
    """Addition and multiplication tables of F_q, q = p^2."""

    def __init__(self, p):
        self.base = B = FieldExt(p, 2)
        q = self.q = B.q
        self.MT = np.array([[B.mul(a, b) for b in range(q)] for a in range(q)], dtype=np.int32)
        self.AT = np.array([[B.add(a, b) for b in range(q)] for a in range(q)], dtype=np.int32)
        self.NEG = np.array([B.scal(p - 1, a) for a in range(q)], dtype=np.int32)
        self.INV = np.array([0] + [B.power(a, q - 2) for a in range(1, q)], dtype=np.int32)
        self.FROB = np.array([B.power(a, p) for a in range(q)], dtype=np.int32)

    def mul(self, a, b):
        return self.MT[a, b]

    def sub(self, a, b):
        return self.AT[a, self.NEG[b]]

    def dot(self, u, v):
        return self.AT[self.AT[self.MT[u[0], v[0]], self.MT[u[1], v[1]]], self.MT[u[2], v[2]]]

    def cross(self, u, v):
        return [self.sub(self.mul(u[1], v[2]), self.mul(u[2], v[1])),
                self.sub(self.mul(u[2], v[0]), self.mul(u[0], v[2])),
                self.sub(self.mul(u[0], v[1]), self.mul(u[1], v[0]))]

    def general_configs(self):
        q = self.q
        v = np.arange(2, q, dtype=np.int32)
        b, c, e, f = (g.ravel() for g in np.meshgrid(v, v, v, v, indexing="ij"))
        m = self.sub(self.mul(b, f), self.mul(c, e))
        ok = (b != c) & (e != f) & (e != b) & (f != c) & (m != 0)
        ok &= self.AT[self.sub(m, self.sub(f, c)), self.sub(e, b)] != 0
        return np.stack([b[ok], c[ok], e[ok], f[ok]], axis=1)

    def keys(self, cfg, perms):
        """Smallest key of each configuration over the given relabellings."""
        q, n = self.q, len(cfg)
        one, zero = np.ones(n, dtype=np.int32), np.zeros(n, dtype=np.int32)
        pts = [[one, zero, zero], [zero, one, zero], [zero, zero, one], [one, one, one],
               [one, cfg[:, 0], cfg[:, 1]], [one, cfg[:, 2], cfg[:, 3]]]
        best = np.full(n, q ** 4, dtype=np.int64)
        for perm in perms:
            Q = [pts[i] for i in perm]
            cr = [self.cross(Q[1], Q[2]), self.cross(Q[2], Q[0]), self.cross(Q[0], Q[1])]
            lam = [self.dot(cr[i], Q[3]) for i in range(3)]
            mu = [self.mul(lam[1], lam[2]), self.mul(lam[0], lam[2]), self.mul(lam[0], lam[1])]
            img = []
            for j in (4, 5):
                x = [self.mul(mu[i], self.dot(cr[i], Q[j])) for i in range(3)]
                s = self.INV[x[0]]
                img.append((self.mul(x[1], s).astype(np.int64), self.mul(x[2], s).astype(np.int64)))
            key = ((img[0][0] * q + img[0][1]) * q + img[1][0]) * q + img[1][1]
            np.minimum(best, key, out=best)
        return best


_T, _CFG, _TOWERS = None, None, {}


def _keys_chunk(perms):
    return _T.keys(_CFG, perms)


def canonical(T, cfg, procs):
    global _T, _CFG
    _T, _CFG = T, cfg
    perms = list(itertools.permutations(range(6)))
    if len(cfg) < 200000 or procs == 1:
        return T.keys(cfg, perms)
    chunks = [perms[i::procs] for i in range(procs)]
    with Pool(procs) as pool:                       # forked: the workers see _T and _CFG
        parts = pool.map(_keys_chunk, chunks)
    return np.minimum.reduce(parts)


def decode(key, q):
    return [int(key // q ** 3), int(key // q ** 2 % q), int(key // q % q), int(key % q)]


def _init(p, kmax):
    base = FieldExt(p, 2)
    for k in range(1, kmax + 1):
        _TOWERS[k] = Tower(base, k)


def _sums(cfg):
    return [_TOWERS[k].char_sum(cfg) for k in sorted(_TOWERS)]


def _sums12(cfg):
    return [_TOWERS[k].char_sum(cfg) for k in (1, 2)]


def _sum3(cfg):
    return _TOWERS[3].char_sum(cfg)


def reduced_sums(pool, todo, uniq, assoc, frob):
    """S_1, S_2 for every class; S_3 once per orbit of association and conjugation.
    Conjugation preserves all S_k.  Association preserves them up to the twist sign, which S_1 fixes;
    classes with S_1 = 0 are computed directly."""
    s12 = pool.map(_sums12, todo, chunksize=8)
    index = {int(k): i for i, k in enumerate(uniq)}
    parent = list(range(len(todo)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(todo)):
        for other in (index[int(assoc[i])], index[int(frob[i])]):
            parent[find(i)] = find(other)
    direct = sorted({find(i) for i in range(len(todo))} | {i for i in range(len(todo)) if s12[i][0] == 0})
    s3 = dict(zip(direct, pool.map(_sum3, [todo[i] for i in direct], chunksize=2)))
    out = []
    for i in range(len(todo)):
        if i in s3:
            out.append(s12[i] + [s3[i]])
        else:
            r = find(i)
            assert abs(s12[i][0]) == abs(s12[r][0]) and s12[i][1] == s12[r][1]
            out.append(s12[i] + [s3[r] * (1 if s12[i][0] == s12[r][0] else -1)])
    return out, len(direct)


def cmd_enumerate(args):
    p, t0 = args.p, time.time()
    T = Tables(p)
    q = T.q
    cfg = T.general_configs()
    keys = canonical(T, cfg, args.procs)
    uniq, counts = np.unique(keys, return_counts=True)
    reps = np.array([decode(k, q) for k in uniq], dtype=np.int32)
    assoc = canonical(T, reps[:, [0, 2, 1, 3]], 1)
    frob = canonical(T, T.FROB[reps], 1)
    print(f"q={q}: {len(cfg)} configurations, {len(uniq)} classes, {int((assoc == uniq).sum())} self-associated, "
          f"{int((frob == uniq).sum())} conjugation-stable ({time.time() - t0:.1f}s)", file=sys.stderr)
    todo = [tuple(int(v) for v in r) for r in reps]
    with Pool(args.procs, initializer=_init, initargs=(p, 3)) as pool:
        if args.reduced:
            sums, n3 = reduced_sums(pool, todo, uniq, assoc, frob)
            print(f"q={q}: S_3 computed directly for {n3} of {len(todo)} classes", file=sys.stderr)
        else:
            sums = pool.map(_sums, todo, chunksize=2)
    with open(args.out, "w") as fh:
        for r, n, a, fr, s in zip(reps, counts, assoc, frob, sums):
            fh.write(json.dumps({"p": q, "prime": p, "cfg": [int(v) for v in r], "orbit": int(n),
                                 "assoc": decode(int(a), q), "frob": decode(int(fr), q), "S": s}) + "\n")
    print(f"q={q}: wrote {len(reps)} rows to {args.out} ({time.time() - t0:.1f}s)", file=sys.stderr)


def cmd_selftest(args):
    """Fast character sums against direct sums over the points of P^2."""
    ok = True
    for p, kmax in ((3, 2), (5, 1)):
        T = Tables(p)
        base, q = T.base, T.q
        cfg = T.general_configs()
        rng = np.random.default_rng(p)
        for k in range(1, kmax + 1):
            tw = Tower(base, k)
            F = tw.F
            for i in rng.choice(len(cfg), size=2, replace=False):
                c = tuple(int(v) for v in cfg[i])
                b, cc, e, f = (tw.emb[v] for v in c)
                lines = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, b, cc), (1, e, f)]
                pts = [(1, y, z) for y in range(F.q) for z in range(F.q)] + [(0, 1, z) for z in range(F.q)] + [(0, 0, 1)]
                slow = 0
                for pt in pts:
                    val = 1
                    for ln in lines:
                        v = 0
                        for coef, coord in zip(ln, pt):
                            v = F.add(v, F.mul(coef, coord))
                        val = F.mul(val, v)
                    slow += int(F.chi[val])
                fast = tw.char_sum(c)
                ok &= fast == slow
                print(f"p={p} k={k} cfg={c}: fast {fast} direct {slow} {'ok' if fast == slow else 'MISMATCH'}")
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("enumerate")
    e.add_argument("--p", type=int, required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--procs", type=int, default=8)
    e.add_argument("--reduced", action="store_true",
                   help="compute S_3 once per orbit of association and conjugation")
    e.set_defaults(func=cmd_enumerate)
    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    args = ap.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()
