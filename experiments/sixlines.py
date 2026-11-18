#!/usr/bin/env python
"""
sixlines.py -- double planes branched along six rational lines over F_p (note 2).

X_C : w^2 = l1 l2 l3 l4 l5 l6, with the six lines in general position (no three concurrent).  Its
resolution is a K3 surface of geometric Picard rank >= 16 on which Frobenius acts by p on the 16
classes of the 15 exceptional curves and the hyperplane class.  The trace of Frobenius^k on the
remaining rank-6 part of H^2 is the character sum

    S_k(C) = sum over P in P^2(F_{p^k}) of chi(l1 l2 ... l6 (P)).

Lines are written as points of the dual plane.  Up to PGL_3 the first four are the standard frame
(1,0,0), (0,1,0), (0,0,1), (1,1,1) and the last two are (1,b,c), (1,e,f), so a configuration is the
tuple (b, c, e, f).  The associated configuration (Gale transform) is (b, e, c, f).

  enumerate   all configurations over F_p, one representative per PGL_3(F_p) x S_6 class, with the
              class size, the class of the associated configuration, and S_1, S_2, S_3
  selftest    compare the fast character sums with direct point counts

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool

import numpy as np

FRAME = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]


# ---------------------------------------------------------------- configurations

def general_configs(p):
    """All (b, c, e, f) for which the six points are in general position (no three collinear)."""
    out = []
    vals = range(2, p)
    for b in vals:
        for c in vals:
            if b == c:
                continue
            for e in vals:
                if e == b:
                    continue
                for f in vals:
                    if f == e or f == c:
                        continue
                    m = (b * f - c * e) % p
                    if m == 0:
                        continue
                    if (m - (f - c) + (e - b)) % p == 0:
                        continue
                    out.append((b, c, e, f))
    return np.array(out, dtype=np.int64).reshape(-1, 4)


def points_of(cfg):
    """(n, 6, 3) array of the six points of each configuration."""
    n = len(cfg)
    pts = np.zeros((n, 6, 3), dtype=np.int64)
    for i, q in enumerate(FRAME):
        pts[:, i, :] = q
    pts[:, 4, 0] = 1
    pts[:, 4, 1] = cfg[:, 0]
    pts[:, 4, 2] = cfg[:, 1]
    pts[:, 5, 0] = 1
    pts[:, 5, 1] = cfg[:, 2]
    pts[:, 5, 2] = cfg[:, 3]
    return pts


def canonical_keys(cfg, p):
    """Smallest code b + p(c + p(e + p f))-style key over the 720 relabellings of the six points."""
    inv = np.array([0] + [pow(a, p - 2, p) for a in range(1, p)], dtype=np.int64)
    pts = points_of(cfg)
    best = np.full(len(cfg), p ** 4, dtype=np.int64)
    for perm in itertools.permutations(range(6)):
        q = pts[:, perm, :]
        cr = [np.cross(q[:, 1], q[:, 2]) % p, np.cross(q[:, 2], q[:, 0]) % p, np.cross(q[:, 0], q[:, 1]) % p]
        lam = [np.einsum("ij,ij->i", cr[i], q[:, 3]) % p for i in range(3)]
        mu = [lam[1] * lam[2] % p, lam[0] * lam[2] % p, lam[0] * lam[1] % p]
        img = []
        for j in (4, 5):
            x = [mu[i] * (np.einsum("ij,ij->i", cr[i], q[:, j]) % p) % p for i in range(3)]
            s = inv[x[0]]
            img.append((x[1] * s % p, x[2] * s % p))
        key = ((img[0][0] * p + img[0][1]) * p + img[1][0]) * p + img[1][1]
        np.minimum(best, key, out=best)
    return best


def decode(key, p):
    f = key % p
    e = key // p % p
    c = key // p ** 2 % p
    b = key // p ** 3
    return int(b), int(c), int(e), int(f)


# ---------------------------------------------------------------- finite fields F_{p^k}, k <= 3

class Field:
    def __init__(self, p, k):
        self.p, self.k, self.q = p, k, p ** k
        self.mod = self._irreducible()
        codes = np.arange(self.q)
        self.digits = np.stack([codes // p ** j % p for j in range(k)], axis=1)
        sq = np.zeros(self.q, dtype=bool)
        for a in range(1, self.q):
            sq[self.mul(a, a)] = True
        chi = np.where(sq, 1, -1).astype(np.int8)
        chi[0] = 0
        self.chi = chi
        self.X = np.zeros((p,) * k, dtype=np.int8)
        self.X[tuple(self.digits[:, j] for j in range(k))] = chi
        # Frobenius orbits on codes
        theta_pows = [self.power(p if k > 1 else 1, j * p) if k > 1 else 1 for j in range(k)]
        frob = np.zeros(self.q, dtype=np.int64)
        for a in range(self.q):
            v = 0
            for j in range(k):
                v = self.add(v, self.scal(int(self.digits[a, j]), theta_pows[j]))
            frob[a] = v
        seen = np.zeros(self.q, dtype=bool)
        reps = []
        for a in range(1, self.q):
            if seen[a]:
                continue
            n, x = 0, a
            while not seen[x]:
                seen[x] = True
                x = frob[x]
                n += 1
            reps.append((a, n))
        self.orbits = reps

    def _irreducible(self):
        p, k = self.p, self.k
        if k == 1:
            return [0, 1]
        for coeffs in itertools.product(range(p), repeat=k):
            poly = list(coeffs) + [1]
            if all(sum(cf * pow(x, i, p) for i, cf in enumerate(poly)) % p for x in range(p)):
                return poly
        raise RuntimeError("no irreducible polynomial")

    def vec(self, a):
        return [a // self.p ** j % self.p for j in range(self.k)]

    def code(self, v):
        return sum(int(d) * self.p ** j for j, d in enumerate(v))

    def add(self, a, b):
        return self.code([(x + y) % self.p for x, y in zip(self.vec(a), self.vec(b))])

    def scal(self, s, a):
        return self.code([s * x % self.p for x in self.vec(a)])

    def mul(self, a, b):
        p, k = self.p, self.k
        va, vb = self.vec(a), self.vec(b)
        prod = [0] * (2 * k - 1)
        for i, x in enumerate(va):
            if x:
                for j, y in enumerate(vb):
                    prod[i + j] += x * y
        for d in range(2 * k - 2, k - 1, -1):
            t = prod[d] % p
            if t:
                for i in range(k):
                    prod[d - k + i] -= t * self.mod[i]
            prod[d] = 0
        return self.code([x % p for x in prod[:k]])

    def power(self, a, n):
        r = 1
        while n:
            if n & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            n >>= 1
        return r


def char_sum(F, cfg):
    """S = sum over x, y in F_q^* of chi(x y (x + y + 1)(x + b y + c)(x + e y + f))."""
    b, c, e, f = cfg
    p, k, X = F.p, F.k, F.X
    axes = tuple(range(k))
    total = 0
    for y, weight in F.orbits:
        d = F.digits[y]
        t1 = [int(d[0]) + 1] + [int(v) for v in d[1:]]
        t2 = [b * int(d[0]) + c] + [b * int(v) for v in d[1:]]
        t3 = [e * int(d[0]) + f] + [e * int(v) for v in d[1:]]
        prod = X * np.roll(X, [-(v % p) for v in t1], axes)
        prod = prod * np.roll(X, [-(v % p) for v in t2], axes)
        prod = prod * np.roll(X, [-(v % p) for v in t3], axes)
        total += weight * int(F.chi[y]) * int(prod.sum(dtype=np.int64))
    return total


def direct_sum(F, cfg):
    """The same sum by direct field arithmetic over all of P^2(F_q) (slow; for the self-test)."""
    b, c, e, f = cfg
    lines = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, b, c), (1, e, f)]
    q = F.q
    pts = [(1, y, z) for y in range(q) for z in range(q)] + [(0, 1, z) for z in range(q)] + [(0, 0, 1)]
    total = 0
    for pt in pts:
        val = 1
        for ln in lines:
            v = 0
            for coef, coord in zip(ln, pt):
                v = F.add(v, F.scal(coef, coord))
            val = F.mul(val, v)
            if val == 0:
                break
        total += int(F.chi[val])
    return total


_FIELDS = {}


def _init(p, kmax):
    for k in range(1, kmax + 1):
        _FIELDS[k] = Field(p, k)


def _sums(cfg):
    return [char_sum(_FIELDS[k], cfg) for k in sorted(_FIELDS)]


# ---------------------------------------------------------------- commands

def cmd_enumerate(args):
    p = args.p
    t0 = time.time()
    cfg = general_configs(p)
    keys = canonical_keys(cfg, p)
    own = ((cfg[:, 0] * p + cfg[:, 1]) * p + cfg[:, 2]) * p + cfg[:, 3]
    uniq, counts = np.unique(keys, return_counts=True)
    assert set(uniq.tolist()) <= set(own.tolist())
    reps = np.array([decode(k, p) for k in uniq], dtype=np.int64)
    assoc = canonical_keys(reps[:, [0, 2, 1, 3]], p)
    print(f"p={p}: {len(cfg)} configurations, {len(uniq)} classes, "
          f"{int((assoc == uniq).sum())} self-associated ({time.time() - t0:.1f}s)", file=sys.stderr)
    with Pool(args.procs, initializer=_init, initargs=(p, 3)) as pool:
        sums = pool.map(_sums, [tuple(int(v) for v in r) for r in reps], chunksize=4)
    with open(args.out, "w") as fh:
        for r, n, a, s in zip(reps, counts, assoc, sums):
            fh.write(json.dumps({"p": p, "cfg": [int(v) for v in r], "orbit": int(n),
                                 "assoc": list(decode(int(a), p)), "S": s}) + "\n")
    print(f"p={p}: wrote {len(reps)} rows to {args.out} ({time.time() - t0:.1f}s)", file=sys.stderr)


def cmd_selftest(args):
    ok = True
    for p, kmax in ((5, 3), (7, 2), (11, 2)):
        cfg = general_configs(p)
        rng = np.random.default_rng(p)
        for k in range(1, kmax + 1):
            F = Field(p, k)
            for i in rng.choice(len(cfg), size=3, replace=False):
                c = tuple(int(v) for v in cfg[i])
                fast, slow = char_sum(F, c), direct_sum(F, c)
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
    e.add_argument("--procs", type=int, default=6)
    e.set_defaults(func=cmd_enumerate)
    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    args = ap.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()
