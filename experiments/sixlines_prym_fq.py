#!/usr/bin/env python
"""
sixlines_prym_fq.py -- the fourfold of a six-line double plane over F_q, q = p^2 (note 16).

The construction and the comparison of sixlines_prym.py (note 14), with F_q in place of F_p as the field of
definition of the six lines.  q = 25 and 49 are 1 mod 4, so i is in F_q and only one of the two types of note 14,
Lemma 4 can be computed: the half c4 of tau is rational and T_k = - sum chi_4(F(x)) over E(F_{q^k}), k = 1, 2, 3.
If c4 is conjugate to -c4, Frobenius inverts J and the pair is only recorded.

All choices of square roots are made once, in F_q, and carried to the larger fields by embeddings; a point that is
conjugate to its negative is written (X, om r) with r in F_q and om a square root of a fixed non-square n0 of F_q.

With a prime q = 1 mod 4 the script runs on the files of note 2 and must reproduce the comparisons of note 14
for the pairs of that type (--run 13 17 29 37).

  --check q b c e f     one configuration (codes of sixlines_fq.py), all 30 ordered pairs, verbose
  --run q [q ...]       configurations of sixlines_q{q}.jsonl (q = p^2) or sixlines_p{q}.jsonl (q prime)
  --select all|irreducible   which configurations (default irreducible: sign + and P irreducible)
  --pairs n             stop after n computed ordered pairs per configuration (default 4; 0 = all 30 are tried)
  --workers n           processes (default 8)
  --out FILE            records
  --summary FILE        recompute the statistics from an output file

Run from experiments/ with python -I.
"""
import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sixlines_cm as cm  # noqa: E402
import sixlines_prym as sp  # noqa: E402  (defines pr_line, pr_half, pr_chain, pr_eval, pr_value)
from sixlines_fq import FieldExt  # noqa: E402

pari = sp.pari
HERE = os.path.dirname(os.path.abspath(__file__))

pari(r"""
pq_code(x) = {
  \\ the integer code of an element of F_q: digits of its polynomial in the generator, as in sixlines_fq.py
  my(f = x.pol, p = x.p);
  if(type(f) != "t_POL", return(lift(f) % p));
  sum(j = 0, poldegree(f), (lift(polcoef(f, j)) % p) * p^j);
}
""")

pari(r"""
pq_sqrt(x) = {
  \\ the square root of x in its field with the smaller code
  my(r = sqrt(x));
  if(pq_code(-r) < pq_code(r), -r, r);
}
""")

pari(r"""
pq_lambda(L, idx) = {
  \\ [N, D]: the cross-ratio of the points of L1..L4 on the line L5 = t L6 is N(t) / D(t)
  my(l5 = L[idx[5],], l6 = L[idx[6],], one = L[1, 1]^0, B, Z, s = vector(4), co);
  for(j = 1, 3, Z = vector(3, k, (k == j) * one); B = matconcat([l5; l6; Z]); if(matdet(B) != 0, break));
  for(i = 1, 4, co = matsolve(B~, L[idx[i],]~); s[i] = -(co[1] * 't + co[2]) / co[3]);
  [(s[1] - s[3]) * (s[2] - s[4]), (s[1] - s[4]) * (s[2] - s[3])];
}
""")

pari(r"""
pq_curve(N, D, n0) = {
  \\ as pr_curve of sixlines_prym.py, over F_q; n0 is the fixed non-square.
  \\ [1, W, [X(a), X(b), X(c)], rat, Xc, is Yc rational, are b1 and b2 rational] or [0] (degenerate)
  my(one = n0^0, l5, l6, lc, Dl, d0, d1, d2, C3, be, sq, c = one, rat = 0, m, W, e1, Xc);
  if(subst(D, 't, 0) == 0 || polcoef(D, 2) == 0 || polcoef(N, 2) == 0 || subst(N, 't, 0) == 0 || polcoef(D, 1) == 0, return([0]));
  l5 = subst(N, 't, 0) / subst(D, 't, 0); l6 = polcoef(N, 2) / polcoef(D, 2); lc = polcoef(N, 1) / polcoef(D, 1);
  if(l5 == l6 || l5 == 1 || l6 == 1 || lc == l5, return([0]));
  Dl = poldisc(N - 'la * D, 't);
  d0 = polcoef(Dl, 0, 'la); d1 = polcoef(Dl, 1, 'la); d2 = polcoef(Dl, 2, 'la);
  C3 = (1 + (l5 - l6) * 'x) * (d2 * (l5 * 'x + 1)^2 + d1 * 'x * (l5 * 'x + 1) + d0 * 'x^2);
  if(poldegree(C3) != 3 || poldisc(C3) == 0, return([0]));
  be = [-1 / l5, 1 / (1 - l5), 0 * one];
  sq = vector(3, j, issquare(subst(C3, 'x, be[j])));
  if(!sq[3], error("the points over lam = oo are not rational"));
  if(!sq[1] && !sq[2], rat = 3);
  if(sq[1] + sq[2] == 1, c = n0; rat = if(sq[1], 2, 1));
  m = c * pollead(C3);
  W = m^2 * c * subst(C3, 'x, 'x / m);
  e1 = m / (l6 - l5); Xc = m / (lc - l5);
  if(subst(W, 'x, e1) != 0 || (Xc - e1)^2 != subst(deriv(W, 'x), 'x, e1), error("lam_c is not a half of tau"));
  Xc = 2 * e1 - Xc;
  [1, W, m * be, rat, Xc, issquare(subst(W, 'x, Xc)), issquare(poldisc(Dl, 'la))];
}
""")

pari(r"""
pq_point(W, X, n0) = {
  \\ [X, r, flag]: the point (X, r) if flag = 1, the point (X, om r) with om^2 = n0 if flag = 0; r in F_q
  my(w = subst(W, 'x, X));
  if(w == 0 || issquare(w), [X, if(w == 0, w, pq_sqrt(w)), 1], [X, pq_sqrt(w / n0), 0]);
}
""")

pari(r"""
pq_lift(P, emb, om) = [ffmap(emb, P[1]), if(P[3], ffmap(emb, P[2]), om * ffmap(emb, P[2]))];
""")

pari(r"""
pq_findq(W, r, rat, c4, n0) = {
  \\ as pr_findq: Q = (u + v om, u' + v' om) in E(F_{q^2}) with Q - Frob_q(Q) = the sum of the two points of r that
  \\ are conjugate to their negatives; [u, v, u', v'] in F_q, least in the order of the codes among the solutions
  \\ for which no branch point falls on the pole or its conjugate
  my(q = n0.p^n0.f, g = ffgen(q^2, 'h), one = g^0, emb = ffembed(n0, g), inv = ffinvmap(emb), z = ffprimroot(g), om = sqrt(ffmap(emb, n0)),
     a2 = ffmap(emb, polcoef(W, 2)), a4 = ffmap(emb, polcoef(W, 1)), a6 = ffmap(emb, polcoef(W, 0)), E = ellinit([0, a2, 0, a4, a6]),
     pt = vector(3, j, pq_lift(r[j], emb, om)), c4b = pq_lift(c4, emb, om), tgt = [0], e2 = (q^2 - 1) / 2,
     X, rhs, Y, ys, co, Q, rho, B4, cs, csb, bad, key, best = 0, bestbad = 0, bq = 0, bqbad = 0);
  co = (t -> my(u = (t + t^q) / 2, v = (t - t^q) / (2 * om)); [ffmap(inv, u), ffmap(inv, v)]);
  for(j = 1, 3, if(j != rat, tgt = elladd(E, tgt, pt[j])));
  for(n = 0, q^2 - 1,
    X = if(n == 0, 0 * one, z^n);
    rhs = ((X + a2) * X + a4) * X + a6;
    if(rhs == 0, ys = [0 * one], if(rhs^e2 == one, Y = sqrt(rhs); ys = [Y, -Y], next));
    foreach(ys, y,
      Q = [X, y];
      if(ellsub(E, Q, [X^q, y^q]) != tgt, next);
      rho = ellsub(E, elladd(E, Q, [X^q, y^q]), pt[rat]);
      B4 = concat([Q], vector(3, j, ellsub(E, elladd(E, pt[j], rho), Q)));
      cs = elladd(E, c4b, rho); csb = elladd(E, ellneg(E, c4b), rho);
      bad = (cs == [0] || csb == [0]);
      for(i = 1, 4, if(B4[i] == cs || B4[i] == csb || B4[i] == [0], bad = 1));
      co = concat((t -> my(u = (t + t^q) / 2, v = (t - t^q) / (2 * om)); [ffmap(inv, u), ffmap(inv, v)])(X),
                  (t -> my(u = (t + t^q) / 2, v = (t - t^q) / (2 * om)); [ffmap(inv, u), ffmap(inv, v)])(y));
      key = apply(pq_code, co);
      if(bad,
        if(bestbad === 0 || lex(key, bestbad) < 0, bestbad = key; bqbad = co),
        if(best === 0 || lex(key, best) < 0, best = key; bq = co))));
  if(best === 0 && bestbad === 0, error("no point Q with Q - Frob(Q) = a + b"));
  if(best === 0, bqbad, bq);
}
""")

pari(r"""
pq_trace(W, r, rat, Qc, c4, i4b, n0, k) = {
  \\ [n0, n1, n2, n3]: the number of points x of E(F_{q^k}) at which chi_4(F(x)) = i^j (at the poles of order 4
  \\ the leading coefficient is used).  T_k = -(n0 - n2) - (n1 - n3) i.  0 if degenerate.
  my(q0 = n0.p^n0.f, q = q0^k, big = if(rat && k % 2, q^2, q), g = ffgen(big, 'g), one = g^0, emb0 = ffembed(n0, g),
     z = ffprimroot(g), zs = z^((big - 1) / (q - 1)),
     a2 = ffmap(emb0, polcoef(W, 2)), a4 = ffmap(emb0, polcoef(W, 1)), a6 = ffmap(emb0, polcoef(W, 0)), E = ellinit([0, a2, 0, a4, a6]),
     om = if(rat, sqrt(ffmap(emb0, n0)), 0), pt = vector(3, j, pq_lift(r[j], emb0, om)), c4b = pq_lift(c4, emb0, om),
     i4 = ffmap(emb0, i4b), tab = [one, i4, -one, -i4], Q = [0], rho, B4, cs, csb, D, Dc, ch, cnt = vector(4),
     exq = (q - 1) / 4, exh = (q - 1) / 2, X, rhs, Y, ys, v, w, exc = List(), Rs = List(), R, hit, ext, g2, emb, E2, U, V, tgt, zeros);
  if(rat == 0,
    rho = ellneg(E, elladd(E, elladd(E, pt[1], pt[2]), pt[3])),
    Q = [ffmap(emb0, Qc[1]) + ffmap(emb0, Qc[2]) * om, ffmap(emb0, Qc[3]) + ffmap(emb0, Qc[4]) * om];
    tgt = [0]; for(j = 1, 3, if(j != rat, tgt = elladd(E, tgt, pt[j])));
    if(ellsub(E, Q, [Q[1]^q0, Q[2]^q0]) != tgt, error("Q - Frob(Q) is not a + b"));
    rho = ellsub(E, elladd(E, Q, [Q[1]^q0, Q[2]^q0]), pt[rat]));
  B4 = concat([Q], vector(3, j, ellsub(E, elladd(E, pt[j], rho), Q)));
  cs = elladd(E, c4b, rho); csb = elladd(E, ellneg(E, c4b), rho);
  for(i = 1, 4, if(B4[i] == cs || B4[i] == csb, return(0)); for(j = i + 1, 4, if(B4[i] == B4[j], return(0))));
  if(cs == [0] || csb == [0], return(0));
  zeros = select(P -> P != [0], B4);
  D = concat(vector(#zeros, j, [zeros[j], 2]), [[cs, -4]]);
  ch = pr_chain(E, D);
  if(#zeros == 4, cnt[1]++);
  for(n = 0, q - 1,
    X = if(n == 0, 0 * one, zs^n);
    rhs = ((X + a2) * X + a4) * X + a6;
    if(rhs == 0, ys = [0 * one], if(rhs^exh == one, Y = sqrt(rhs); ys = [Y, -Y], next));
    foreach(ys, y,
      v = pr_eval(ch, X, y);
      if(v[1] == 0 || v[2] == 0, listput(exc, [X, y]); next);
      w = if(big == q, (v[1] * v[2]^3)^exq, (v[1] / v[2])^exq);
      hit = 0; for(j = 1, 4, if(w == tab[j], cnt[j]++; hit = 1; break));
      if(!hit, error("the character value is not a fourth root of unity"))));
  if(#exc,
    ext = 2; while(big^ext < 10^5 || ext % 2 == 0, ext++);
    g2 = ffgen(big^ext, 'h); emb = ffembed(g, g2);
    E2 = ellinit([0, ffmap(emb, a2), 0, ffmap(emb, a4), ffmap(emb, a6)]);
    for(i = 1, 12, R = random(E2); if(R != [0], listput(Rs, R)));
    D = ffmap(emb, D); tab = ffmap(emb, tab));
  foreach(Vec(exc), pp,
    if(pp == B4[1] || pp == B4[2] || pp == B4[3] || pp == B4[4], next);
    X = ffmap(emb, pp[1]); Y = ffmap(emb, pp[2]);
    if(pp == cs,
      until(U != cs && V != cs && V != [0],
        until(rhs != 0 && rhs^exh == one, U = zs^random(q - 1); rhs = ((U + a2) * U + a4) * U + a6);
        U = [U, sqrt(rhs)]; V = ellneg(E, elladd(E, cs, U)));
      Dc = ffmap(emb, concat(vector(#zeros, j, [zeros[j], 2]), [[U, 4], [V, 4]]));
      w = pr_value(E2, Dc, X, Y, Vec(Rs))^exq,
      w = pr_value(E2, D, X, Y, Vec(Rs))^exq);
    hit = 0; for(j = 1, 4, if(w == tab[j], cnt[j]++; hit = 1; break));
    if(!hit, error("the character value at an exceptional point is not a fourth root of unity")));
  cnt;
}
""")

pari(r"""
pq_setup(q, modulus) = {
  \\ [generator w of F_q with the modulus of sixlines_fq.py, the non-square n0 with the least code, a square root i4 of -1]
  my(p = factor(q)[1, 1], w = if(q == p, ffgen(p, 'w), ffgen(Mod(1, p) * modulus, 'w)), one = w^0, el, n0 = 0);
  for(c = 1, q - 1, el = if(q == p, c * one, (c % p) * one + (c \ p) * w); if(!issquare(el), n0 = el; break));
  [w, n0, pq_sqrt(-one)];
}
""")

GP = {n: pari(n) for n in ["pq_lambda", "pq_curve", "pq_point", "pq_findq", "pq_trace", "pq_setup"]}
FIELD = {}


def field(q):
    """(prime, generator, n0, i4) for F_q, with the encoding of sixlines_fq.py when q = p^2."""
    if q not in FIELD:
        p = int(pari.factor(q)[0][0])
        if q == p:
            w, n0, i4 = GP["pq_setup"](q, 0)
        else:
            mod = FieldExt(p, 2)._irreducible()
            w, n0, i4 = GP["pq_setup"](q, pari(f"x^2 + {mod[1]}*x + {mod[0]}"))
        FIELD[q] = (p, w, n0, i4)
    return FIELD[q]


def element(code, q):
    p, w, n0, i4 = field(q)
    one = w ** 0
    return code * one if q == p else (code % p) * one + (code // p) * w


def forms(cfg, q):
    b, c, e, f = (element(v, q) for v in cfg)
    o, z = element(1, q), element(0, q)
    return [(o, z, z), (z, o, z), (z, z, o), (o, o, o), (o, b, c), (o, e, f)]


def one_pair(q, cfg, pair):
    p, w, n0, i4 = field(q)
    rest = [i for i in range(1, 7) if i not in pair]
    L = pari.matrix(6, 3, [x for row in forms(cfg, q) for x in row])
    N, D = GP["pq_lambda"](L, rest + list(pair))
    cur = GP["pq_curve"](N, D, n0)
    if int(cur[0]) == 0:
        return {"status": "degenerate"}
    W, bx, rat, Xc, yrat, brat = cur[1], cur[2], int(cur[3]), cur[4], bool(int(cur[5])), bool(int(cur[6]))
    out = {"c4_rational": yrat, "b_rational": brat, "rational_point": rat}
    if not yrat:
        out["status"] = "anticommuting"
        return out
    r = pari.Vec([GP["pq_point"](W, bx[j], n0) for j in range(3)])
    c4 = GP["pq_point"](W, Xc, n0)
    Qc = GP["pq_findq"](W, r, rat, c4, n0) if rat else pari.Vec([n0, n0, n0, n0])
    T = []
    for k in (1, 2, 3):
        cnt = GP["pq_trace"](W, r, rat, Qc, c4, i4, n0, k)
        if cnt == 0:
            return {"status": "degenerate"}
        n = [int(v) for v in cnt]
        T.append((-(n[0] - n[2]), -(n[1] - n[3])))
    out["T"] = T
    c = sp.quartic(T, q)
    if c is None:
        out["status"] = "not integral"
        return out
    out["c"] = c
    out["wedge"] = sp.wedge(c, q)
    out["c3_is_q_conj_c1"] = c[2] == (q * c[0][0], -q * c[0][1])
    out["status"] = "computed"
    return out


def describe(q, r):
    e = sp.surface_e(r["S"])
    p = field(q)[0]
    signs = cm.choose_sign(q, e)
    rec = {"q": q, "cfg": r["cfg"], "e": list(e), "signs": signs}
    if signs == [1]:
        P = cm.pol(cm.charpoly(q, e))
        rec["degrees"] = sorted(int(pari.poldegree(f)) for f, m in zip(*pari.factor(P)) for _ in range(int(m)))
        rec["ordinary"] = e[0] % p != 0
    if "frob" in r:
        rec["defined_over_prime_field"] = r["frob"] == r["cfg"] or r["frob"] == r["assoc"]
    return rec


def config_record(job):
    q, r, max_pairs = job
    rec = describe(q, r)
    e = tuple(rec["e"])
    rec["pairs"] = []
    seen = Counter()
    used = 0
    for pair in itertools.permutations(range(1, 7), 2):
        res = one_pair(q, r["cfg"], pair)
        seen[res["status"]] += 1
        if res["status"] == "computed":
            rec["pairs"].append({"pair": list(pair), "T": res["T"], "c": res["c"], "compare": sp.compare(res, e),
                                 "c3_is_q_conj_c1": res["c3_is_q_conj_c1"], "rational_point": res["rational_point"]})
            used += 1
            if max_pairs and used >= max_pairs:
                break
        elif res["status"] == "not integral":
            rec["pairs"].append({"pair": list(pair), "T": res["T"], "compare": "not integral"})
    rec["seen"] = dict(seen)
    return rec


def tally(recs):
    stats = Counter()
    for rec in recs:
        signs, seen = rec["signs"], rec["seen"]
        kind = "sign +" if signs == [1] else ("sign -" if signs == [-1] else f"signs {signs}")
        if signs == [1]:
            kind += ", P irreducible" if rec["degrees"] == [6] else ", P reducible"
        stats["configurations"] += 1
        stats[kind] += 1
        comps = [x["compare"] for x in rec["pairs"]]
        if comps and signs != [1]:
            stats[f"{kind}: has pairs of commuting type (not comparable)"] += 1
        elif comps:
            stats[f"{kind}: computed"] += 1
            if all(x in ("equal", "twist") for x in comps):
                stats[f"{kind}: every pair gives P or its twist"] += 1
            else:
                stats[f"{kind}: MISMATCH"] += 1
                print(f"  mismatch: q = {rec['q']} cfg = {rec['cfg']} e = {rec['e']} {rec['pairs']}", flush=True)
            if rec.get("defined_over_prime_field") is False:
                stats[f"{kind}: computed, not defined over the prime field"] += 1
        elif seen.get("anticommuting"):
            stats[f"{kind}: only pairs where Frobenius inverts J"] += 1
        else:
            stats[f"{kind}: no usable pair"] += 1
        if signs == [1]:
            good = [x for x in rec["pairs"] if x["compare"] in ("equal", "twist")]
            stats["pairs computed (sign +)"] += len(comps)
            stats["pairs equal"] += sum(x == "equal" for x in comps)
            stats["pairs twist"] += sum(x == "twist" for x in comps)
            stats["pairs with c3 = q conj(c1)"] += sum(x["c3_is_q_conj_c1"] for x in good)
            stats["pairs with conjugate points"] += sum(x["rational_point"] != 0 for x in good)
            stats["pairs where Frobenius inverts J (sign +)"] += seen.get("anticommuting", 0)
        stats["pairs degenerate"] += seen.get("degenerate", 0)
    return stats


def load(q, select):
    p = int(pari.factor(q)[0][0])
    name = f"sixlines_p{q}.jsonl" if q == p else f"sixlines_q{q}.jsonl"
    rows = [json.loads(line) for line in open(os.path.join(HERE, name))]
    if select == "irreducible":
        keep = []
        for r in rows:
            d = describe(q, r)
            if d["signs"] == [1] and d["degrees"] == [6]:
                keep.append(r)
        rows = keep
    return rows


def check(q, cfg):
    name = f"sixlines_p{q}.jsonl" if field(q)[0] == q else f"sixlines_q{q}.jsonl"
    row = [r for r in map(json.loads, open(os.path.join(HERE, name))) if r["cfg"] == cfg][0]
    e = sp.surface_e(row["S"])
    print(f"q = {q}, configuration {cfg}: surface (e1, e2, e3) = {e}, signs {cm.choose_sign(q, e)}")
    for pair in itertools.permutations(range(1, 7), 2):
        t0 = time.time()
        res = one_pair(q, cfg, pair)
        line = f"  pair {pair}: {res['status']}"
        if res["status"] == "computed":
            line += f"  c = {res['c']}  wedge = {res['wedge']}  {sp.compare(res, e)}  ({time.time() - t0:.1f} s)"
        print(line, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", nargs=5, type=int, metavar=("q", "b", "c", "e", "f"))
    ap.add_argument("--run", nargs="*", type=int, default=[])
    ap.add_argument("--select", choices=["all", "irreducible"], default="irreducible")
    ap.add_argument("--pairs", type=int, default=4)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out")
    ap.add_argument("--summary")
    args = ap.parse_args()
    if args.summary:
        recs = [json.loads(line) for line in open(args.summary)]
        for q in sorted({r["q"] for r in recs}):
            print(f"q = {q}: {dict(sorted(tally([r for r in recs if r['q'] == q]).items()))}")
    if args.check:
        check(args.check[0], args.check[1:])
    fh = open(args.out, "w") if args.out else None
    for q in args.run:
        t0 = time.time()
        rows = load(q, args.select)
        with Pool(args.workers) as pool:
            recs = pool.map(config_record, [(q, r, args.pairs) for r in rows], chunksize=2)
        if fh:
            for rec in recs:
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
        print(f"q = {q}: {dict(sorted(tally(recs).items()))}  ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
