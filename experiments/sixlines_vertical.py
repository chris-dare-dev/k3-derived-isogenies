#!/usr/bin/env python
"""
sixlines_vertical.py -- rational l-steps that change the endomorphism ring (note 7).

For a unimodular lattice L with Frobenius F and an odd prime l prime to q, an l-step is a neighbour

    L' = Z x~/l + { y in L : <x, y> = 0 mod l },      x~ = x mod l L,  <x~, x~> = 0 mod l^2,

one for every isotropic line F_l x of L / l L.  The step is defined over F_q exactly when the line is
stable under F, that is when x is an eigenvector of F on L / l L, with eigenvalue lam say; the twisted
equivalence that realises it is defined over F_q exactly when lam = q (lam = -q for the quadratic twist).

For every polynomial of notes 2, 3 and 5 and every odd prime l prime to q that divides the index of
Z[F, V] in O_E, the script takes each class of lattices, lists the F-stable lines of L / l L, builds
the neighbour for the isotropic ones, finds its class among the lattices of the polynomial, compares
the endomorphism rings at l, and tests whether the neighbour is a translate (a / abar) L of L by a
uniformiser a at one prime above l ("translation": a step inside the torsor of Proposition A).  It also computes the l-part of coker(F -+ q), which is dual to Br(X)[l^oo].

Checks made on the way: every neighbour is unimodular, F-stable and isomorphic to exactly one class of
the list; the numbers of steps satisfy  n(i -> j, lam) |Aut j| = n(j -> i, q^2/lam) |Aut i|.

  python -I sixlines_vertical.py [--out FILE] [--only q,e1,e2,e3]

Reads the result files of notes 2, 3 and 5.  Run from experiments/.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sixlines_cm as cm  # noqa: E402
from sixlines_cm import pari  # noqa: E402

pari(r"""
sv_order(nf, L) = {
  \\ the ring { a in E : a L < L }, as a Z-basis in integral-basis coordinates: it is the trace dual
  \\ of the product of L with its trace dual
  my(P = nf.pol, n = poldegree(P), w = nf.zk, T, Ld, M);
  T = matrix(n, n, i, j, nfelttrace(nf, Mod(w[i], P) * Mod(w[j], P)));
  Ld = (L~ * T)^(-1);
  M = sl_hnfq(matconcat(vector(n, j, sl_mulmat(nf, nfbasistoalg(nf, Ld[,j])) * L)));
  sl_hnfq((M~ * T)^(-1));
}
""")
pari(r"""
sv_gram(nf, q, al) = {
  my(P = nf.pol, n = poldegree(P), w = nf.zk, cw = vector(n, j, sl_conjel(nf, q, w[j])));
  matrix(n, n, i, j, nfelttrace(nf, al * Mod(w[i], P) * cw[j]));
}
""")
pari(r"""
sv_aut(nf, L, zeta, w) = sum(a = 0, w - 1, sl_hnfq(sl_mulmat(nf, zeta^a) * L) == L);
""")
pari(r"""
sv_coker(nf, L, g, l) = {
  \\ exponents of l in the elementary divisors of g on L
  my(M = L^(-1) * sl_mulmat(nf, g) * L, d);
  if(denominator(M) != 1, error("not stable"));
  d = matsnf(M);
  vecsort(select(e -> e > 0, apply(x -> valuation(x, l), Vec(d))));
}
""")
pari(r"""
sv_steps(nf, q, L, G, Rord, l) = {
  \\ the F-stable lines of L / l L.  One entry [lam, mu, kind, L', stable] per line: lam the eigenvalue,
  \\ mu the dimension of its eigenspace, kind = 1 if the line is isotropic (then L' is the neighbour)
  \\ and 0 if not, stable = 1 if the line is stable under the ring Rord = End(L)
  my(P = nf.pol, n = poldegree(P), Fm = sl_mulmat(nf, Mod(variable(P), P)), Li = L^(-1),
     FL = Li * Fm * L, GL = L~ * G * L, rts = polrootsmod(P, l), out = List(), om);
  if(denominator(FL) != 1 || denominator(GL) != 1 || abs(matdet(GL)) != 1, error("not a unimodular F-lattice"));
  om = vector(n, j, Li * sl_mulmat(nf, nfbasistoalg(nf, Rord[,j])) * L);
  for(i = 1, #rts,
    my(lam = lift(rts[i]), W = lift(matker(Mod(FL - lam * matid(n), l))), mu = #W, f, x, gx, c, k, t, xt, K, M, Lp, st);
    forvec(e = vector(mu, i, [0, l - 1]),
      f = 0; for(i = 1, mu, if(e[i] != 0, f = e[i]; break));
      if(f != 1, next);
      x = W * e~;
      st = 1;
      for(j = 1, n, if(matrank(Mod(concat(Mat(x), om[j] * x), l)) > 1, st = 0; break));
      gx = GL * x; c = x~ * gx;
      if(c % l, listput(out, [lam, mu, 0, 0, st]); next);
      k = 0; for(j = 1, n, if(gx[j] % l, k = j; break));
      t = lift(Mod(-c / l, l) / Mod(2 * gx[k], l));
      xt = x; xt[k] += l * t;
      if((xt~ * GL * xt) % l^2, error("lift is not isotropic mod l^2"));
      K = matrix(n, n, a, b, if(b == k, if(a == k, l, 0),
                                if(a == b, 1, if(a == k, lift(-Mod(gx[b], l) / Mod(gx[k], l)), 0))));
      M = sl_hnfq(concat(K, xt / l));
      if(abs(matdet(M)) != 1, error("neighbour has the wrong index"));
      if(denominator(M~ * GL * M) != 1, error("neighbour is not integral"));
      if(denominator(M^(-1) * FL * M) != 1, error("neighbour is not F-stable"));
      Lp = sl_hnfq(L * M);
      listput(out, [lam, mu, 1, Lp, st])));
  Vec(out);
}
""")
pari(r"""
sv_dec(nf, q, l) = {
  \\ the primes of E above l: [e, f, self-conjugate, ramification index over E_0, v(F - q), v(F + q)]
  my(P = nf.pol, t = variable(P), F = Mod(t, P), V = q^2 / F, s = F + V, h = minpoly(s, 'x),
     nf0 = nfinit(h), dec = idealprimedec(nf, l), dec0 = idealprimedec(nf0, l), out = List());
  for(i = 1, #dec,
    my(pr = dec[i], self = (idealhnf(nf, nfgaloisapply(nf, lift(V), pr)) == idealhnf(nf, pr)), er = 0, up);
    for(j = 1, #dec0,
      up = idealhnf(nf, l, subst(lift(nfbasistoalg(nf0, dec0[j].gen[2])), 'x, s));
      er = max(er, idealval(nf, up, pr)));
    listput(out, [pr.e, pr.f, self, er, idealval(nf, F - q, pr), idealval(nf, F + q, pr)]));
  vecsort(Vec(out));
}
""")
pari(r"""
sv_translates(nf, q, L, l, K) = {
  \\ for every prime pr above l: the lattice equal to (a / abar) L at l and to L elsewhere, with a a
  \\ uniformiser at pr and a unit at the other primes above l.  K bounds the exponent of l involved.
  my(dec = idealprimedec(nf, l), fact, a, M);
  vector(#dec, k,
    fact = matrix(#dec, 2, i, j, if(j == 1, dec[i], i == k));
    a = nfbasistoalg(nf, idealappr(nf, fact));
    M = sl_hnfq(sl_mulmat(nf, a / sl_conjel(nf, q, a)) * L);
    sl_hnfq(concat(sl_latcap(M, L / l^K), l^K * L)));
}
""")
pari(r"""
sv_local(O, l) = { my(d = abs(matdet(O))); mathnf(concat(O, l^valuation(d, l) * matid(#O))); }
""")
GP_order, GP_gram, GP_aut, GP_coker, GP_steps, GP_local, GP_dec, GP_translates = (
    pari(n) for n in ("sv_order", "sv_gram", "sv_aut", "sv_coker", "sv_steps", "sv_local", "sv_dec", "sv_translates"))


def val(n, l):
    v = 0
    while n % l == 0:
        n //= l
        v += 1
    return v


def order_mod(a, l):
    a %= l
    d, x = 1, a
    while x != 1:
        x = x * a % l
        d += 1
    return d


def pattern(P, q, l):
    """How P factors mod l: [(kind, multiplicity)] with kind '+q', '-q', 'lam=..' or 'deg n'."""
    fa = pari.factormod(P, l)
    out = []
    for i in range(len(fa[0])):
        g, m = fa[0][i], int(fa[1][i])
        if int(pari.poldegree(g)) == 1:
            lam = int(pari("(g)->lift(-polcoef(g,0))")(g))
            kind = "+q" if (lam - q) % l == 0 else "-q" if (lam + q) % l == 0 else f"lam={lam}"
        else:
            kind = f"deg {int(pari.poldegree(g))}"
        out.append([kind, m])
    return sorted(out)


def graph(bnf, q, classes, l):
    """The l-steps between the lattice classes of one polynomial."""
    nf = cm.GP_nf(bnf)
    w, zeta = (pari("(b)->b.tu[1]")(bnf), pari("(b)->nfbasistoalg(b.nf, b.tu[2])")(bnf))
    F = pari("(nf)->Mod(variable(nf.pol), nf.pol)")(nf)
    orders = [GP_order(nf, L) for L, al in classes]
    local = [GP_local(o, l) for o in orders]          # the ring at l: End(L) + l^v O_E
    names = {}
    for o in local:
        names.setdefault(str(o), (val(int(pari.matdet(o)), l), len(names)))
    level = [names[str(o)][0] for o in local]
    oid = [names[str(o)][1] for o in local]
    aut = [int(GP_aut(nf, L, zeta, w)) for L, al in classes]
    verts = []
    for i, (L, al) in enumerate(classes):
        verts.append({"level": level[i], "order": oid[i], "aut": aut[i],
                      "coker+": [int(x) for x in GP_coker(nf, L, F - q, l)],
                      "coker-": [int(x) for x in GP_coker(nf, L, F + q, l)]})
    edges = Counter()          # (i, j, lam, is a translation) -> number of lines
    lines = defaultdict(Counter)
    for i, (L, al) in enumerate(classes):
        G = GP_gram(nf, q, al)
        shifts = list(GP_translates(nf, q, L, l, max(level) + 3))
        for lam, mu, kind, Lp, st in GP_steps(nf, q, L, G, orders[i], l):
            lam, mu, st = int(lam), int(mu), int(st)
            lines[i][(lam, mu, "isotropic" if int(kind) else "anisotropic", "stable" if st else "not stable")] += 1
            if not int(kind):
                continue
            op = GP_order(nf, Lp)
            hits = [j for j in range(len(classes)) if orders[j] == op
                    and cm.GP_sl_iso(bnf, q, Lp, al, classes[j][0], classes[j][1])]
            assert len(hits) == 1, f"neighbour matches {len(hits)} classes"
            tr = any(Lp == M for M in shifts)
            assert not tr or oid[hits[0]] == oid[i]
            edges[(i, hits[0], lam, tr)] += 1
    total = Counter()
    for (i, j, lam, tr), n in edges.items():
        total[(i, j, lam)] += n
    for (i, j, lam), n in total.items():
        mu = int(pari("(q,lam,l)->lift(Mod(q,l)^2/Mod(lam,l))")(q, lam, l))
        assert n * aut[j] == total[(j, i, mu)] * aut[i], "the step and its inverse do not balance"
    return verts, edges, lines, local


def summarise(q, l, verts, edges, lines, orders):
    """Per (ring, eigenvalue): the lines and where the steps go."""
    nord = len({v["order"] for v in verts})
    rel = {}
    reps = {}
    for i, v in enumerate(verts):
        reps.setdefault(v["order"], orders[i])

    def contains(a, b):        # ring a contains ring b
        return int(pari.denominator(pari.matsolve(reps[a], reps[b]))) == 1
    rows = defaultdict(lambda: {"classes": 0, "lines": Counter(), "to": Counter()})
    for i, v in enumerate(verts):
        for lam in sorted({k[0] for k in lines[i]}):
            key = (v["order"], lam)
            rows[key]["classes"] += 1
            for k, n in lines[i].items():
                if k[0] == lam:
                    rows[key]["lines"][k[1:]] += n
    for (i, j, lam, tr), n in edges.items():
        a, b = verts[i]["order"], verts[j]["order"]
        if tr:
            kind = "translation"
        elif a == b:
            kind = "same ring"
        elif contains(b, a):
            kind = "up"
        elif contains(a, b):
            kind = "down"
        else:
            kind = "across"
        rel[(a, b)] = kind
        rows[(a, lam)]["to"][(kind, b, verts[j]["level"])] += n
    out = []
    for (o, lam), r in sorted(rows.items()):
        members = [i for i, v in enumerate(verts) if v["order"] == o]
        image = {i: [j for (a, j, m, tr), n in edges.items() if a == i and m == lam and verts[j]["order"] == o
                     for _ in range(n)] for i in members}
        cycles = None
        if all(len(t) == 1 for t in image.values()) and sorted(t[0] for t in image.values()) == members:
            cycles = cm.cycle_type([members.index(image[i][0]) for i in members])
        lev = next(v["level"] for v in verts if v["order"] == o)
        cok = sorted({(tuple(v["coker+"]), tuple(v["coker-"])) for v in verts if v["order"] == o})
        out.append({"order": o, "level": lev, "lam": lam,
                    "d+": order_mod(lam * pow(q, -1, l), l), "d-": order_mod(-lam * pow(q, -1, l), l),
                    "classes": r["classes"],
                    "lines": [[*k, n] for k, n in sorted(r["lines"].items())],
                    "steps": [[*k, n] for k, n in sorted(r["to"].items())],
                    "coker": [[list(a), list(b)] for a, b in cok], "cycles": cycles})
    return nord, out


def load_cases():
    cases = []
    rows = [json.loads(line) for line in open("sixlines_cm.jsonl")]
    rows += [r for r in map(json.loads, open("sixlines_neighbours.jsonl")) if r["p"] == 25]
    rows += [json.loads(line) for line in open("sixlines_q49_cm.jsonl")]
    cases += [(r["p"], r["e"], True) for r in rows]
    for name in ("small", "p17", "p19", "p23"):
        cases += [(r["p"], r["e"], False) for r in map(json.loads, open(f"sixlines_nonsplit_lat_{name}.jsonl"))]
    return cases


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    ap.add_argument("--only", help="q,e1,e2,e3: one polynomial")
    args = ap.parse_args()
    import sixlines_nonsplit as ns
    ns.setup_lattice()
    only = [int(x) for x in args.only.split(",")] if args.only else None
    fh = open(args.out, "w") if args.out else None
    tally = Counter()
    for q, e, split in load_cases():
        if only and [q, *e] != only:
            continue
        P = cm.pol(cm.charpoly(q, e))
        nf0 = pari.nfinit(P)
        index = int(pari("(nf,q)->my(P=nf.pol,F=Mod(variable(P),P)); sl_index(nf,[F,q^2/F])")(nf0, q))
        pp = int(pari.factor(q)[0][0])
        primes = [int(l) for l in pari.factor(index)[0] if int(l) not in (2, pp)]
        tally["polynomials"] += 1
        if not primes:
            continue
        tally["with an odd prime in the index"] += 1
        bnf = pari.bnfinit(P, 1)
        nf = cm.GP_nf(bnf)
        gens = cm.GP_sl_gens(nf, q) if split else pari("sl_gens2")(nf, q)
        classes = cm.lattice_classes(bnf, q, gens, int(pari("(b)->b.tu[1]")(bnf)))
        for l in primes:
            pat = pattern(P, q, l)
            dec = [[int(x) for x in row] for row in GP_dec(nf, q, l)]
            verts, edges, lines, orders = graph(bnf, q, classes, l)
            nord, rows = summarise(q, l, verts, edges, lines, orders)
            keep = ("translation", "same ring")
            change = sorted({s[0] for r in rows for s in r["steps"] if s[0] not in keep})
            rational = sorted({s[0] for r in rows for s in r["steps"]
                               if s[0] not in keep and (r["d+"] == 1 or r["d-"] == 1)})
            if any(s[0] == "same ring" for r in rows for s in r["steps"]):
                tally["pairs with a step that keeps the ring and is not a translation"] += 1
            tally["pairs (polynomial, l)"] += 1
            if any(r["steps"] for r in rows):
                tally["pairs with a rational step"] += 1
            if change:
                tally["pairs with a rational step changing the ring"] += 1
            if rational:
                tally["pairs with such a step and a rational Brauer class"] += 1
            print(f"q={q} e={e} {'split' if split else 'non-split'} l={l} index l^{val(index, l)} "
                  f"P mod l: {pat}; primes [e, f, self-conjugate, e over E_0, v(F-q), v(F+q)]: {dec}; "
                  f"{len(classes)} classes, {nord} rings at l", flush=True)
            for r in rows:
                print(f"    ring {r['order']} (index l^{r['level']}, {r['classes']} classes) lam={r['lam']} "
                      f"d=({r['d+']},{r['d-']}) coker(F-q),(F+q)={r['coker']} lines={r['lines']} steps={r['steps']}"
                      + (f" cycles={r['cycles']}" if r["cycles"] else ""))
            if fh:
                fh.write(json.dumps({
                    "q": q, "e": e, "split": split, "l": l, "index_val": val(index, l), "pattern": pat, "primes": dec,
                    "classes": len(classes), "rings": nord, "vertices": verts, "rows": rows,
                    "edges": [[i, j, lam, int(tr), n] for (i, j, lam, tr), n in sorted(edges.items())]}) + "\n")
    print("summary:", dict(tally))


if __name__ == "__main__":
    main()
