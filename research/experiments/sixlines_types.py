#!/usr/bin/env python
"""
sixlines_types.py -- how Frobenius over F_p acts on NS for the surfaces fixed by conjugation (note 4).

A six-line double plane X over F_q, q = p^2, that is isomorphic to its conjugate has a model over
F_p (up to the twist).  The p-Frobenius of the model acts on the rank-16 lattice N of NS by p g,
with g an involution in O(N).  Its image in O(q_N) = S_6 x Z/2 is

    (pi, 0)   Frobenius permutes the six lines by pi           (a form of type A, B or C)
    (pi, 1)   Frobenius exchanges the two six-line polarisations, and the labels by pi

A self-associated surface (E contains Q(i)) has two forms over F_p, one of each kind, and gets both labels.

Surface side: read off from the configuration (b, c, e, f) and its conjugate.
Lattice side: the surface is fixed exactly when +-F = w^2 with w in End(L); then w is the
p-Frobenius on the rank-6 part (both w and i w when i is in End(L)), and its action on L/2L, carried to the discriminant group of N,
must be in the same conjugacy class.

  python -I sixlines_types.py CM_FILE Q_FILE [--out FILE]

CM_FILE is the output of sixlines_cm.py --neighbours for Q_FILE.  Run from experiments/.
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402
from sixlines_forms import ALL_PERMS, normal_form  # noqa: E402
from sixlines_fq import Tables  # noqa: E402

pari(r"""
sl_gram(nf, q, L, al) = {
  my(n = #L, v = vector(n, i, nfbasistoalg(nf, L[,i])));
  matrix(n, n, i, j, trace(al * v[i] * sl_conjel(nf, q, v[j])));
}
""")
pari(r"""
sl_sqrt(nf) = {
  my(P = nf.pol, t = variable(P), r);
  my(out = List());
  r = nfroots(nf, x^2 - Mod(t, P)); if(#r, listput(out, [r[1], 1]));
  r = nfroots(nf, x^2 + Mod(t, P)); if(#r, listput(out, [r[1], -1]));
  Vec(out);
}
""")
pari(r"""
sl_action(nf, L, g) = L^(-1) * sl_mulmat(nf, g) * L;
""")
GP_gram, GP_sqrt, GP_action = pari("sl_gram"), pari("sl_sqrt"), pari("sl_action")


# ---------------------------------------------------------------- quadratic forms on (Z/2)^6

def table(images):
    """Values of the linear map with the given images of the basis bits, on all of (Z/2)^n."""
    out = [0]
    for v in images:
        out += [x ^ v for x in out]
    return out


def isometries(QA, QB, first=False):
    """Linear maps g of (Z/2)^6 with QB[g x] = QA[x] (Q with values in Z/4).  Tables of values."""
    found = []

    def extend(img):
        k = len(img).bit_length() - 1
        if k == 6:
            found.append(tuple(img))
            return
        for v in range(64):
            if QB[v] != QA[1 << k]:
                continue
            new = img + [x ^ v for x in img]
            if all(QB[new[x]] == QA[x] for x in range(len(img), len(new))):
                extend(new)
                if first and found:
                    return
    extend([0])
    return found


def compose(g, h):
    return tuple(g[h[x]] for x in range(64))


# ---------------------------------------------------------------- the lattice N and its discriminant form

def cycle_class(perm):
    return tuple(cm.cycle_type(list(perm)))


def lattice_N():
    """O(q_N) with every element labelled (cycle type of pi, eps); the values of q_N."""
    duads = list(itertools.combinations(range(6), 2))
    half = Fraction(1, 2)
    gens = [[Fraction(int(i == j)) for i in range(16)] for j in range(16)]
    for i in range(6):
        gens.append([-half if i in d else Fraction(0) for d in duads] + [half])
    M = pari("[" + ";".join(",".join(str(2 * g[r]) for g in gens) for r in range(16)) + "]")
    H = pari.mathnf(M)
    B = [[Fraction(int(H[c][r]), 2) for c in range(16)] for r in range(16)]          # columns: a basis of N
    D = [-2] * 15 + [2]
    G = [[sum(B[r][i] * D[r] * B[r][j] for r in range(16)) for j in range(16)] for i in range(16)]
    assert abs(cm._det(G)) == 64
    Ginv = cm.mat_inv(G)
    dual2 = cm.int_mat([[2 * x for x in row] for row in Ginv])                       # columns: 2 x, x in the dual

    def qval(y):
        v = sum(y[i] * G[i][j] * y[j] for i in range(16) for j in range(16)) / 2
        assert v.denominator == 1
        return int(v) % 4
    basis, span = [], {tuple([0] * 16): 0}
    for c in range(16):
        y = [dual2[r][c] for r in range(16)]
        key = tuple(v % 2 for v in y)
        if key in span:
            continue
        bit = 1 << len(basis)
        basis.append(y)
        for k, code in list(span.items()):
            span[tuple((a + b) % 2 for a, b in zip(k, key))] = code | bit
    assert len(basis) == 6 and len(span) == 64
    QN = [0] * 64
    for code in range(64):
        y = [sum(basis[i][r] for i in range(6) if code >> i & 1) for r in range(16)]
        QN[code] = qval(y)
    Binv = cm.mat_inv(B)
    group = isometries(QN, QN)
    assert len(group) == 1440
    centre = [g for g in group if all(compose(g, h) == compose(h, g) for h in group)]
    assert len(centre) == 2
    z = next(g for g in centre if g != tuple(range(64)))
    labels = {}
    for perm in ALL_PERMS:
        images = []
        for y in basis:
            v = [sum(B[r][c] * y[c] for c in range(16)) for r in range(16)]          # orthogonal coordinates
            w = [Fraction(0)] * 16
            for d, duad in enumerate(duads):
                w[duads.index(tuple(sorted((perm[duad[0]], perm[duad[1]]))))] = v[d]
            w[15] = v[15]
            y2 = [sum(Binv[r][c] * w[c] for c in range(16)) for r in range(16)]
            assert all(x.denominator == 1 for x in y2)
            images.append(span[tuple(int(x) % 2 for x in y2)])
        g = tuple(table(images))
        assert g in group
        for eps, h in ((0, g), (1, compose(g, z))):
            assert h not in labels
            labels[h] = (cycle_class(perm), eps)
    assert len(labels) == 1440
    return QN, labels


# ---------------------------------------------------------------- surface side

def relabel(T, cfg, perm):
    n = len(cfg)
    one, zero = np.ones(n, dtype=np.int32), np.zeros(n, dtype=np.int32)
    pts = [[one, zero, zero], [zero, one, zero], [zero, zero, one], [one, one, one],
           [one, cfg[:, 0], cfg[:, 1]], [one, cfg[:, 2], cfg[:, 3]]]
    return normal_form(T, [pts[i] for i in perm])


def gale(T, cfg):
    """The associated configuration with its natural labelling.  For the points e1, e2, e3 and the
    columns of A = (1,1,1), (1,b,c), (1,e,f), the associated points are the rows of A followed by
    e1, e2, e3; in normal form that is (b, e, c, f) relabelled by (14)(25)(36)."""
    return relabel(T, cfg[:, [0, 2, 1, 3]], (3, 4, 5, 0, 1, 2))


def surface_labels(T, by, configs):
    """For each surface of one polynomial that is fixed by conjugation, its labels (cycle type, eps).
    A surface has one label, or two when it is self-associated (then it has two forms over F_p)."""
    pairs = sorted({tuple(sorted([tuple(c), tuple(by[tuple(c)]["assoc"])])) for c in configs})
    out = []
    for pair in pairs:
        r = by[pair[0]]
        frob = tuple(r["frob"])
        if frob not in pair:
            continue
        cfg = np.array([r["cfg"]], dtype=np.int32)
        target = T.FROB[cfg][0].tolist()
        mine = []
        for eps, start in ((0, cfg), (1, gale(T, cfg))):
            hits = set()
            for perm in ALL_PERMS:
                # the labelled Gale transform commutes with relabelling
                assert gale(T, relabel(T, cfg, perm))[0].tolist() == relabel(T, gale(T, cfg), perm)[0].tolist()
                if relabel(T, start, perm)[0].tolist() == target:
                    hits.add(cycle_class(perm))
            assert len(hits) <= 1, hits
            mine += [(ct, eps) for ct in hits]
        assert mine
        out.append(tuple(sorted(mine)))
    return out


# ---------------------------------------------------------------- lattice side

def lattice_labels(bnf, q, gens, w, QN, labels):
    nf = cm.GP_nf(bnf)
    classes = cm.lattice_classes(bnf, q, gens, w)
    roots = [(r, int(sg)) for r, sg in GP_sqrt(nf)]
    info = {"classes": len(classes), "sqrt_signs": [sg for r, sg in roots], "fixed": 0, "stable": 0,
            "agree": True, "labels": []}
    for L, al in classes:
        fixed = bool(cm.GP_sl_iso(bnf, q, cm.GP_sl_translate_p(nf, q, L), al, L, al))
        G = GP_gram(nf, q, L, al)
        Gi = [[int(G[c][r]) for c in range(6)] for r in range(6)]
        assert abs(int(pari.matdet(G))) == 1
        QT = [0] * 64
        for x in range(64):
            v = [x >> i & 1 for i in range(6)]
            QT[x] = (-sum(v[i] * Gi[i][j] * v[j] for i in range(6) for j in range(6))) % 4
        phi = None
        mine = []
        for root, sg in roots:
            W = GP_action(nf, L, root)
            if int(pari.denominator(W)) != 1:
                continue
            wt = table([sum((int(W[c][r]) % 2) << r for r in range(6)) for c in range(6)])
            assert all(QT[wt[x]] == QT[x] for x in range(64))
            phi = phi or isometries(QT, QN, first=True)[0]
            g = [0] * 64
            for x in range(64):
                g[phi[x]] = phi[wt[x]]
            mine.append(labels[tuple(g)])
        info["fixed"] += fixed
        info["stable"] += bool(mine)
        info["agree"] &= fixed == bool(mine)
        if mine:
            info["labels"].append(tuple(sorted(mine)))
    return info


def show(lbls):
    names = {(1, 1, 1, 1, 1, 1): "1", (1, 1, 1, 1, 2): "A", (1, 1, 2, 2): "B", (2, 2, 2): "C"}
    c = Counter(" or ".join(f"({names.get(ct, str(ct))}, {'swap' if eps else 'lines'})" for ct, eps in one)
                for one in lbls)
    return ", ".join(f"{n} x {t}" for t, n in sorted(c.items())) or "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmfile")
    ap.add_argument("qfile")
    ap.add_argument("--out")
    args = ap.parse_args()
    QN, labels = lattice_N()
    print(f"O(q_N): 1440 elements = S_6 x Z/2, centre of order 2; values of q_N: {dict(sorted(Counter(QN).items()))}")
    rows = [json.loads(line) for line in open(args.qfile)]
    by = {tuple(r["cfg"]): r for r in rows}
    q, prime = rows[0]["p"], rows[0]["prime"]
    T = Tables(prime)
    tally, out = Counter(), []
    for rec in (json.loads(line) for line in open(args.cmfile)):
        if rec["p"] != q:
            continue
        surf = surface_labels(T, by, rec["configs"])
        line = f"  e={rec['e']} k3={rec['k3']:g} index={rec['index']}: surfaces fixed {len(surf)}: {show(surf)}"
        res = {"q": q, "e": rec["e"], "k3": rec["k3"],
               "surface": [[[list(ct), eps] for ct, eps in one] for one in surf]}
        if rec.get("N_total") is None:
            print(line + "; lattices not computed")
            tally["lattices not computed"] += 1
            out.append(res)
            continue
        P = cm.pol(cm.charpoly(q, rec["e"]))
        bnf = pari.bnfinit(P, 1)
        gens = cm.GP_sl_gens(cm.GP_nf(bnf), q)
        info = lattice_labels(bnf, q, gens, rec["w"], QN, labels)
        ok = info["agree"] and sorted(info["labels"]) == sorted(surf)
        sq = {(): "neither F nor -F is a square", (1,): "F is a square", (-1,): "-F is a square",
              (1, -1): "F and -F are squares"}[tuple(info["sqrt_signs"])]
        print(line + f"\n      lattices: {sq}; fixed by the p-step {info['fixed']}, stable under the root "
              f"{info['stable']}: {show(info['labels'])} -> {'match' if ok else 'MISMATCH'}")
        tally["match" if ok else "MISMATCH"] += 1
        tally["with fixed surfaces" if surf else "without fixed surfaces"] += 1
        res.update({"sqrt_signs": info["sqrt_signs"], "fixed": info["fixed"], "stable": info["stable"],
                    "lattice": [[[list(ct), eps] for ct, eps in one] for one in info["labels"]], "match": ok})
        out.append(res)
    print(f"q={q}: " + ", ".join(f"{k}: {v}" for k, v in sorted(tally.items())))
    if args.out:
        with open(args.out, "w") as fh:
            for res in out:
                fh.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
