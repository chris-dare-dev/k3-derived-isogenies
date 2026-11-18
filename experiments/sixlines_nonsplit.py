#!/usr/bin/env python
"""
sixlines_nonsplit.py -- the count for non-split forms over F_p (note 5).

Surface side: the forms of sixlines_forms.py (six lines over F_{p^2}, Frobenius an involution of
type A, B or C on them), enumerated for larger p with a faster character sum.  No data over
F_{p^2} is used.

Lattice side: unimodular hermitian lattices L over the orders containing

    R' = Z[F, V, (F^2 - 1)/2, (V^2 - 1)/2, (F - V)/2] + the unit-root idempotent part,   V = p^2/F,

each with the class of F on L/2L in O(q_N) = S_6 x Z/2.  The prediction is

    #{forms of type pi with polynomial P(t) or P(-t)} = sum over lattices of class (pi, lines) of
        #( double cosets: S_6-centraliser, centraliser in O(q_N), roots of unity of End(L) ),

which is 2 per lattice when E has only the roots of unity +-1.

  enumerate --p P --out FILE [--procs N]
  lattices FORMS [--out FILE] [--limit N]
  selftest

Run from experiments/ with python -I.
"""
import argparse
import json
import os
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sixlines import Field  # noqa: E402
from sixlines_forms import (FRAME, INVOLUTION, centraliser, cols, const, general_position,  # noqa: E402
                            normal_form, projective_points, split_form)
from sixlines_fq import Tables  # noqa: E402


# ---------------------------------------------------------------- character sums with F_p-coefficients

class FormSum:
    """Sum over P^2(F_{p^k}) of chi(product of linear and quadratic forms with coefficients in F_p).
    For fixed (y, z) a linear form is a (x + t) and a quadratic form is A ((x + h)^2 - delta), so its
    character is a translate of a fixed table."""

    def __init__(self, p, k):
        self.p, self.k = p, k
        F = Field(p, k)
        q = self.q = F.q
        self.chi = F.chi.astype(np.int64)
        self.X = F.X
        self.axes = tuple(range(k))
        self.pw = p ** np.arange(k, dtype=np.int64)
        dig = F.digits.astype(np.int64)
        sq = np.array([F.mul(a, a) for a in range(q)], dtype=np.int64)
        self.R = np.empty((q,) + (p,) * k, dtype=np.int8)
        idx = tuple(dig[:, j] for j in range(k))
        for delta in range(q):
            self.R[delta][idx] = F.chi[((dig[sq] - dig[delta]) % p) @ self.pw]
        ys = [0] + [a for a, n in F.orbits]
        wts = [1] + [n for a, n in F.orbits]
        zero = np.zeros((1, k), dtype=np.int64)
        self.y = np.concatenate([dig[ys], dig[[1]]])                 # rows (y, 1), then the row (1, 0)
        self.z = np.concatenate([np.tile(dig[[1]], (len(ys), 1)), zero])
        self.yy = np.concatenate([dig[sq[ys]], dig[[1]]])
        self.yz = np.concatenate([dig[ys], zero])
        self.zz = self.z
        self.w = wts + [1]
        self.inv = [0] + [pow(a, p - 2, p) for a in range(1, p)]

    def sum(self, linear, quadratic):
        p, inv = self.p, self.inv
        tabs, shifts, const_ = [], [], 1
        for a, b, c in linear:
            const_ *= int(self.chi[a])
            shifts.append((-(b * self.y + c * self.z) * inv[a]) % p)
            tabs.append(None)
        for A, By, Bz, Cyy, Cyz, Czz in quadratic:
            const_ *= int(self.chi[A])
            i2 = inv[2 * A % p]
            shifts.append((-(By * self.y + Bz * self.z) * i2) % p)
            d = ((By * By * self.yy + 2 * By * Bz * self.yz + Bz * Bz * self.zz) * (i2 * i2 % p)
                 - (Cyy * self.yy + Cyz * self.yz + Czz * self.zz) * inv[A]) % p
            tabs.append((d @ self.pw).tolist())
        shifts = [s.tolist() for s in shifts]
        total = 0
        for r, wt in enumerate(self.w):
            prod = None
            for tab, sh in zip(tabs, shifts):
                src = self.X if tab is None else self.R[tab[r]]
                rolled = np.roll(src, sh[r], self.axes)
                prod = rolled if prod is None else prod * rolled
            total += wt * int(prod.sum(dtype=np.int64))
        return const_ * (total + 1)


_T, _SUMS = None, {}


def form_data(T, p, rational, pairs):
    """Linear and quadratic forms over F_p for the six lines, after a change of coordinates that makes
    every coefficient of x non-zero."""
    lines = list(rational) + list(pairs)
    for s in range(p):
        for t in range(p):
            xs = [int(T.AT[T.AT[ln[0], T.MT[s, ln[1]]], T.MT[t, ln[2]]]) for ln in lines]
            if all(xs):
                moved = [(x, int(ln[1]), int(ln[2])) for x, ln in zip(xs, lines)]
                lin = moved[:len(rational)]
                quad = []
                for a, b, c in moved[len(rational):]:
                    ab, bb, cb = (int(T.FROB[v]) for v in (a, b, c))
                    tr = lambda u, vb, ub, v: int(T.AT[T.MT[u, vb], T.MT[ub, v]])     # noqa: E731
                    quad.append((int(T.MT[a, ab]), tr(a, bb, ab, b), tr(a, cb, ab, c),
                                 int(T.MT[b, bb]), tr(b, cb, bb, c), int(T.MT[c, cb])))
                assert all(v < p for f in lin + quad for v in f)
                return lin, quad
    raise RuntimeError("no change of coordinates found")


def _sums(task):
    kind, six = task
    rational, pairs = split_form(kind, six)
    lin, quad = form_data(_T, _T.base.p, rational, pairs)
    return [_SUMS[k].sum(lin, quad) for k in (1, 2, 3)]


def first_frames(T, p, nonrat, conj):
    n = len(nonrat)
    ok = general_position(T, [const(FRAME[0], n), const(FRAME[1], n), cols(nonrat), cols(conj)])
    i = int(np.argmax(ok))
    frames = {"A": list(FRAME),
              "B": [FRAME[0], FRAME[1], tuple(int(v) for v in nonrat[i]), tuple(int(v) for v in conj[i])]}
    for i in range(n):
        ok = general_position(T, [const(nonrat[i], n), const(conj[i], n), cols(nonrat), cols(conj)])
        if ok.any():
            j = int(np.argmax(ok))
            frames["C"] = [tuple(int(v) for v in x) for x in (nonrat[i], conj[i], nonrat[j], conj[j])]
            break
    return frames


def cmd_enumerate(args):
    global _T
    p, t0 = args.p, time.time()
    _T = T = Tables(p)
    pts = projective_points(T.q)
    nonrat = pts[~(pts < p).all(axis=1)]
    conj = T.FROB[nonrat]
    frames = first_frames(T, p, nonrat, conj)
    tasks, meta = [], []
    for kind in "ABC":
        n = len(nonrat)
        six = [const(x, n) for x in frames[kind]] + [cols(nonrat), cols(conj)]
        idx = np.nonzero(general_position(T, six))[0]
        six = [[c[idx] for c in pt] for pt in six]
        cfg = normal_form(T, six)
        key = T.keys(cfg, centraliser(INVOLUTION[kind]))
        uniq, first, counts = np.unique(key, return_index=True, return_counts=True)
        for i, cnt in zip(first.tolist(), counts.tolist()):
            lines = [tuple(int(c[i]) for c in pt) for pt in six]
            tasks.append((kind, lines))
            meta.append({"p": p, "type": kind, "lines": [list(x) for x in lines],
                         "cfg": [int(v) for v in cfg[i]], "candidates": cnt})
        print(f"p={p}: type {kind}: {len(idx)} forms on the standard frame, {len(uniq)} classes "
              f"({time.time() - t0:.1f}s)", file=sys.stderr)
    for k in (1, 2, 3):
        _SUMS[k] = FormSum(p, k)
    print(f"p={p}: tables ready ({time.time() - t0:.1f}s)", file=sys.stderr)
    with Pool(args.procs) as pool:                       # forked: the workers share the tables
        sums = pool.map(_sums, tasks, chunksize=8)
    with open(args.out, "w") as fh:
        for m, s in zip(meta, sums):
            m["S"] = s
            fh.write(json.dumps(m) + "\n")
    print(f"p={p}: wrote {len(meta)} rows to {args.out} ({time.time() - t0:.1f}s)", file=sys.stderr)


# ---------------------------------------------------------------- lattice side

def setup_lattice():
    import sixlines_cm as cm
    import sixlines_types as ty
    pari = cm.pari
    pari(r"""
sl_gens2(nf, p) = {
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = p^2 / F, pp = factor(p)[1, 1], M = 4,
     n = (pp - 1) * pp^(M - 1), e0, e1);
  e0 = lift(lift(Mod(Mod(1, pp^M) * t, P)^n));
  e1 = 1 - Mod(e0, P) - subst(e0, t, V);
  [(F^2 - 1) / 2, (V^2 - 1) / 2, (F - V) / 2, F, V, e1 * F / p, e1 * V / p];
}
""")
    pari(r"""
sl_frob(nf) = Mod(variable(nf.pol), nf.pol);
""")
    pari(r"""
sl_zeta(bnf) = [bnf.tu[1], nfbasistoalg(bnf.nf, bnf.tu[2])];
""")
    return cm, ty, pari


def double_cosets(Z, H, K, compose):
    """Number of orbits of H x K on Z by x -> h x k."""
    seen, n = set(), 0
    for x in Z:
        if x in seen:
            continue
        n += 1
        todo = [x]
        seen.add(x)
        while todo:
            y = todo.pop()
            for h in H:
                for k in K:
                    v = compose(compose(h, y), k)
                    if v not in seen:
                        seen.add(v)
                        todo.append(v)
    return n


def lattice_counts(cm, ty, pari, bnf, p, QN, labels):
    """For every lattice class: the class of F on L/2L and the number of forms it carries."""
    nf = cm.GP_nf(bnf)
    gens = pari("sl_gens2")(nf, p)
    if not all(bool(cm.GP_sl_integral(nf, x)) for x in gens):
        return {"error": "the generators of R' are not integral"}
    w, zeta = pari("sl_zeta")(bnf)
    w = int(w)
    index = int(cm.GP_sl_index(nf, gens))
    classes = cm.lattice_classes(bnf, p, gens, w)
    if classes is None:
        return {"error": "lattice search too large", "index": index}
    F = pari("sl_frob")(nf)
    ident = tuple(range(64))
    out = {"index": index, "w": w, "classes": len(classes), "labels": Counter(), "forms": Counter()}
    for L, al in classes:
        G = ty.GP_gram(nf, p, L, al)
        Gi = [[int(G[c][r]) for c in range(6)] for r in range(6)]
        assert abs(int(pari.matdet(G))) == 1
        QT = [0] * 64
        for x in range(64):
            v = [x >> i & 1 for i in range(6)]
            QT[x] = (-sum(v[i] * Gi[i][j] * v[j] for i in range(6) for j in range(6))) % 4
        phi = ty.isometries(QT, QN, first=True)[0]

        def carry(elt):
            W = ty.GP_action(nf, L, elt)
            if int(pari.denominator(W)) != 1:
                return None
            wt = ty.table([sum((int(W[c][r]) % 2) << r for r in range(6)) for c in range(6)])
            g = [0] * 64
            for x in range(64):
                g[phi[x]] = phi[wt[x]]
            return tuple(g)
        g = carry(F)
        ct, eps = labels[g]
        out["labels"][(ct, eps)] += 1
        if eps:
            continue
        Z = [x for x in labels if ty.compose(x, g) == ty.compose(g, x)]
        H = [x for x in Z if labels[x][1] == 0]
        K = {ident}
        for a in range(1, w):
            u = carry(zeta ** a)
            if u is not None:
                K.add(u)
        assert all(k in Z for k in K)
        out["forms"][ct] += double_cosets(Z, H, sorted(K), ty.compose)
    return out


NAMES = {(1, 1, 1, 1, 1, 1): "1", (1, 1, 1, 1, 2): "A", (1, 1, 2, 2): "B", (2, 2, 2): "C"}


def cmd_lattices(args):
    cm, ty, pari = setup_lattice()
    cm.LIMIT = args.limit
    QN, labels = ty.lattice_N()
    from sixlines_forms import graeffe
    fh = open(args.out, "w") if args.out else None
    for path in args.forms:
        forms = [json.loads(line) for line in open(path)]
        p = forms[0]["p"]
        stats, groups = Counter(), defaultdict(Counter)
        for f in forms:
            e1, e2, e3, integral = cm.elementary(f["S"])
            if not integral:
                stats["sums not integral"] += 1
                continue
            # the sign of the functional equation is open when e3 = 0; F^2 is the Frobenius of a split
            # configuration, so (F^2 - 1)/2 is integral, and that decides it
            signs = [eps for eps in (1, -1) if (eps == 1 or e3 == 0)
                     and cm.weil_ok(p, cm.charpoly(p, (e1, e2, e3), eps))
                     and cm.two_adic_ok(graeffe(cm.charpoly(p, (e1, e2, e3), eps)))]
            if len(signs) != 1:
                stats["NO valid sign" if not signs else "sign ambiguous"] += 1
                continue
            stats["sign +" if signs[0] == 1 else "sign -"] += 1
            co = cm.charpoly(p, (e1, e2, e3), signs[0])
            if e1 % p == 0:
                stats["not ordinary"] += 1
            elif not bool(pari.polisirreducible(cm.pol(co))):
                stats["ordinary, reducible"] += 1
            else:
                stats["ordinary, irreducible"] += 1
                sign = 1 if e1 > 0 else -1
                groups[(sign * e1, e2, sign * e3)][f["type"]] += 1
        print(f"p={p}: {len(forms)} classes of forms: " + ", ".join(f"{k}: {v}" for k, v in sorted(stats.items())))
        tally, t0 = Counter(), time.time()
        for key in sorted(groups):
            seen = groups[key]
            bnf = pari.bnfinit(cm.pol(cm.charpoly(p, key)), 1)
            res = lattice_counts(cm, ty, pari, bnf, p, QN, labels)
            rec = {"p": p, "e": list(key), "forms": dict(seen)}
            if "error" in res:
                verdict = "not computed: " + res["error"]
                rec.update({"verdict": verdict, "index": res.get("index")})
            else:
                pred = {NAMES[ct]: n for ct, n in res["forms"].items() if NAMES[ct] != "1"}
                ok = all(seen.get(t, 0) == pred.get(t, 0) for t in "ABC")
                verdict = "match" if ok else "MISMATCH"
                rec.update({"verdict": verdict, "index": res["index"], "w": res["w"], "classes": res["classes"],
                            "predicted": pred, "split": res["forms"].get((1, 1, 1, 1, 1, 1), 0),
                            "labels": {f"{NAMES[ct]}{'s' if eps else ''}": n for (ct, eps), n in res["labels"].items()}})
            tally[verdict] += 1
            if args.list or not verdict.startswith("match"):
                print(f"  e={list(key)} forms={dict(sorted(seen.items()))} "
                      + (f"index={rec['index']} w={rec['w']} lattices={rec['classes']} classes={rec['labels']} "
                         f"predicted={dict(sorted(rec['predicted'].items()))} " if "predicted" in rec else "")
                      + f"-> {verdict}", flush=True)
            if fh:
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
        print(f"p={p}: {len(groups)} polynomials: " + ", ".join(f"{k}: {v}" for k, v in sorted(tally.items()))
              + f" ({time.time() - t0:.0f}s)", flush=True)
    if fh:
        fh.close()


def cmd_selftest(args):
    """The fast sums against the values stored by sixlines_forms.py for p = 5 and 7."""
    global _T
    ok = True
    for p in (5, 7):
        path = f"sixlines_forms_p{p}.jsonl"
        if not os.path.exists(path):
            print(f"{path} missing, skipped")
            continue
        _T = Tables(p)
        for k in (1, 2, 3):
            _SUMS[k] = FormSum(p, k)
        rows = [json.loads(line) for line in open(path)]
        bad = sum(_sums((r["type"], [tuple(x) for x in r["lines"]])) != r["S"] for r in rows)
        ok &= bad == 0
        print(f"p={p}: {len(rows)} forms, {bad} with different sums")
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("enumerate")
    e.add_argument("--p", type=int, required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--procs", type=int, default=8)
    e.set_defaults(func=cmd_enumerate)
    c = sub.add_parser("lattices")
    c.add_argument("forms", nargs="+")
    c.add_argument("--out")
    c.add_argument("--list", action="store_true", help="print every polynomial")
    c.add_argument("--limit", type=int, default=200000,
                   help="largest number of intermediate lattices in one search for self-dual lattices")
    c.set_defaults(func=cmd_lattices)
    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    args = ap.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()
