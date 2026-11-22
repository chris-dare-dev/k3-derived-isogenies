#!/usr/bin/env python
"""
sixlines_leakage.py -- leakage audit of six-line surface invariants against the CM labels (note 28).

Question: does any invariant computable from a six-line surface carry information about its label in G beyond the
known bit chi_-1?  The tests were fixed before any invariant was joined to a label (the pre-registration is
reproduced in note 28).

  features   build the labelled set and compute every candidate invariant (PARI; one process).
             L = full labels in G (note 20 sec. 7.5, last 'action' record of each field in sixlines_cmmethod.jsonl);
             C = cosets of ker chi_-1 (notes 24, 26: sixlines_delta_labels.jsonl, sixlines_delta41_cosets.jsonl,
             sixlines_delta41_compare.json).  For each configuration of a labelled surface:
               S  note 19's invariants recomputed with sixlines_characters.candidates (checked against
                  sixlines_characters.jsonl where a record exists), c1_values from the stored quartics, the twist;
               J  J2..J6 of the six lines (sixlines_cmmethod's cm_tmat, cm_J; the stored ratios are checked);
                  d-th power residue symbols of monomials J^a of weight = 0 mod d;
               E  t_s = [ijk][lmn], p_k = sum t_s^(2k) (k = 1..5), the conic C: residue symbols of monomials;
               P  histograms of the stored per-pair curve data (sixlines_prym*.jsonl);
               B2 scalars (absolute invariants) for the relational symbols of v(X) - v(Y), v(X) + v(Y);
             and, for C, the per-pair fourfold data (sixlines_jacobian_character*.jsonl, sixlines_delta41.jsonl).
             Monomial candidates are kept only if invariant under random GL3 x diagonal rescalings (checked on three
             configurations of every field, before any label is read for the statistics).
             Writes sixlines_leakage_features.json.
  audit      exact permutation tests (all n! relabellings of each torsor of L; all balanced splits of each torsor of
             C), combined over torsors by exact convolution, Holm / Bonferroni / BH over all candidates of a
             population; the positive control (delta on C), planted signals for power, random functions for the
             null.  numpy only; --workers processes.  Writes sixlines_leakage.jsonl and prints the tables.

  python -I sixlines_leakage.py features [--out sixlines_leakage_features.json]
  python -I sixlines_leakage.py audit [--features FILE] [--out sixlines_leakage.jsonl] [--workers 12]
  python -I sixlines_leakage.py extra [--shuffles 10] [--plants 200] [--workers 12]
             shuffled-label reruns of population L and planted signals on population C.

Run from experiments/.  Deterministic (fixed seeds).
"""
import argparse
import collections
import itertools
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

H = 0.01            # grid of the convolution, in bit-pairs


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return [json.loads(l) for l in fh if l.strip()]


# ================================================================================================= features

def build_population():
    """The labelled torsors.  A torsor: q, e1, G (cyclic factors), surfaces [{configs, label | coset, root,
    invariants, delta}]."""
    last = collections.OrderedDict()
    for r in load("sixlines_cmmethod.jsonl"):
        last[(r["q"], tuple(r["e"]), r["cmd"])] = r
    L, L_small, L_partial = [], [], []
    for (q, e, cmd), r in last.items():
        if cmd != "action" or r["status"] != "ok":
            continue
        surf = [{"configs": [tuple(c) for c in s["configurations"]], "label": tuple(s["label"]) if s["label"] else None,
                 "root": s["root"], "invariants": s["invariants"]} for s in r["surfaces"]]
        tor = {"pop": "L", "q": q, "e": list(e), "G": list(r["G"]), "chart": r.get("chart", "J2"), "surfaces": surf}
        if any(s["label"] is None for s in surf):
            L_partial.append(tor)
        elif len(surf) <= 2:
            L_small.append(tor)
        else:
            L.append(tor)
    C = []
    for r in load("sixlines_delta_labels.jsonl"):
        surf = [{"configs": [tuple(c) for c in s["configurations"]], "coset": s["orbit"], "root": s["root"],
                 "invariants": s["invariants"], "delta": s["delta"][0] if len(set(s["delta"])) == 1 else None}
                for s in r["surfaces"]]
        C.append({"pop": "C", "q": r["q"], "e1": r["e1"], "G": None, "chart": r["chart"], "surfaces": surf})
    with open(os.path.join(HERE, "sixlines_delta41_compare.json")) as fh:
        cmp41 = json.load(fh)
    for r in load("sixlines_delta41_cosets.jsonl"):
        if "surfaces" not in r:
            continue
        tab = {tuple(sorted(map(tuple, s["configurations"]))): s for s in cmp41[str(r["e1"])]["table"]}
        surf = []
        for s in r["surfaces"]:
            t = tab[tuple(sorted(map(tuple, s["configurations"])))]
            assert t["orbits of ideals"] == s["orbit"]
            surf.append({"configs": [tuple(c) for c in s["configurations"]], "coset": s["orbit"], "root": s["root"],
                         "invariants": s["invariants"], "delta": t["delta"]})
        C.append({"pop": "C", "q": r["q"], "e1": r["e1"], "G": None, "chart": r["chart"], "surfaces": surf})
    return L, L_small, L_partial, C


class Fp:
    """Arithmetic helpers mod a prime p: discrete logarithms to the least primitive root."""

    def __init__(self, p):
        self.p = p
        g = 2
        while True:
            if all(pow(g, (p - 1) // f, p) != 1 for f in prime_factors(p - 1)):
                break
            g += 1
        self.g = g
        self.dlog = {}
        x = 1
        for k in range(p - 1):
            self.dlog[x] = k
            x = x * g % p

    def sym(self, x, d):
        """d-th power residue symbol as dlog mod d, or 'z' for 0."""
        x %= self.p
        return "z" if x == 0 else self.dlog[x] % d


def prime_factors(n):
    out, f = [], 2
    while f * f <= n:
        while n % f == 0:
            out.append(f)
            n //= f
        f += 1
    if n > 1:
        out.append(n)
    return sorted(set(out))


def monomials(d, weights, cond):
    """Exponent vectors a in (Z/d)^k, a != 0, with cond(sum w_k a_k)."""
    out = []
    for a in itertools.product(range(d), repeat=len(weights)):
        if any(a) and cond(sum(w * x for w, x in zip(weights, a))):
            out.append(a)
    return out


def mono_symbol(base_syms, a, d):
    """Residue symbol of a monomial from the symbols of its bases (dlogs mod d, or 'z')."""
    s = 0
    for x, k in zip(base_syms, a):
        if k % d == 0:
            continue
        if x == "z":
            return "z"
        s += k * x
    return s % d


J_WEIGHTS = (2, 3, 4, 5, 6)
E_WEIGHTS = (2, 4, 6, 8, 10, 2)        # p1..p5, C in units of the weight of t_s


def j_family(q):
    """(name, d, a) for the J family: weight = 0 mod d."""
    out = []
    for d in (2, 4, 8):
        if (q - 1) % d or (d == 8 and q != 41):
            continue
        for a in monomials(d, J_WEIGHTS, lambda w: w % d == 0):
            out.append((f"J d={d} a={''.join(map(str, a))}", d, a))
    return out


def e_family(q):
    """(name, d, a) for the E family: weight (in units of t) = 0 mod d; d in {2, 4} only (amendment, see note 28)."""
    out = []
    for d in (2, 4):
        if (q - 1) % d:
            continue
        for a in monomials(d, E_WEIGHTS, lambda w: w % d == 0):
            out.append((f"E d={d} a={''.join(map(str, a))}", d, a))
    return out


def config_quantities(pari, A, p):
    """J2..J6, the ten t_s, p_1..p_5, C, for a 3 x 6 matrix A over F_p (pari Mod entries)."""
    t = pari("cm_tmat")(A)
    J = [int(pari.lift(x)) for x in pari("cm_J")(t)]
    tv = [int(pari.lift(x)) for x in t]
    pk = [sum(pow(x, 2 * k, p) for x in tv) % p for k in range(1, 6)]
    D = lambda i, j, k: int(pari.lift(pari.matdet(pari.matrix(3, 3, [A[r, c] for r in range(3) for c in (i, j, k)]))))
    conic = (D(0, 1, 2) * D(0, 3, 4) * D(1, 3, 5) * D(2, 4, 5) - D(0, 1, 3) * D(0, 2, 4) * D(1, 2, 5) * D(3, 4, 5)) % p
    return {"J": J, "t": tv, "p": pk, "C": conic}


def matrix_of(pari, cfg, p):
    b, c, e, f = cfg
    return pari.matrix(3, 6, [1, 0, 0, 1, 1, 1, 0, 1, 0, 1, b, e, 0, 0, 1, 1, c, f]) * pari.Mod(1, p)


def inv_ratio(J, p):
    """Chart-J2 absolute invariants (the stored I1..I5), or None if J2 = 0."""
    J2, J3, J4, J5, J6 = J
    if J2 % p == 0:
        return None
    i = pow(J2, -1, p)
    return [J4 * i ** 2 % p, J6 * i ** 3 % p, J3 * J3 * i ** 3 % p, J3 * J5 * i ** 4 % p, J5 * J5 * i ** 5 % p]


def chart_invariants(J, p, chart):
    J2, J3, J4, J5, J6 = J
    if chart == "J2":
        return inv_ratio(J, p)
    if J3 % p == 0:
        return None
    i = pow(J3, -1, p)
    return [J2 ** 3 * i ** 2 % p, J2 * J4 * i ** 2 % p, J6 * i ** 2 % p, J4 * J5 * i ** 3 % p, J4 ** 3 * i ** 4 % p,
            J5 ** 3 * i ** 5 % p, J2 ** 2 * J5 * i ** 3 % p]


def scalars(Q, p):
    """Absolute invariants for B2: I1..I5 (chart J2), p_k / p_1^k (k = 2..5), C^2 / p_1^2.  None where undefined."""
    out = {}
    I = inv_ratio(Q["J"], p)
    for k in range(5):
        out[f"I{k + 1}"] = None if I is None else I[k]
    p1 = Q["p"][0]
    for k in range(2, 6):
        out[f"R{k}"] = None if p1 == 0 else Q["p"][k - 1] * pow(p1, -k, p) % p
    out["C2/p1^2"] = None if p1 == 0 else Q["C"] ** 2 * pow(p1, -2, p) % p
    return out


def invariance_filter(pari, fams, cfgs, p, rng, trials=3):
    """Drop monomial candidates whose symbol changes under A -> g A D (g in GL3, D diagonal)."""
    F = Fp(p)
    bad = set()
    for cfg in cfgs:
        A = matrix_of(pari, cfg, p)
        Q0 = config_quantities(pari, A, p)
        for _ in range(trials):
            while True:
                g = pari.matrix(3, 3, [rng.randrange(p) for _ in range(9)]) * pari.Mod(1, p)
                if pari.matdet(g) != 0:
                    break
            Dg = pari.matdiagonal([rng.randrange(1, p) for _ in range(6)]) * pari.Mod(1, p)
            Q1 = config_quantities(pari, g * A * Dg, p)
            for name, d, a in fams:
                if name.startswith("J"):
                    b0, b1 = [F.sym(x, d) for x in Q0["J"]], [F.sym(x, d) for x in Q1["J"]]
                else:
                    b0 = [F.sym(x, d) for x in Q0["p"] + [Q0["C"]]]
                    b1 = [F.sym(x, d) for x in Q1["p"] + [Q1["C"]]]
                if mono_symbol(b0, a, d) != mono_symbol(b1, a, d):
                    bad.add(name)
    return bad


def prym_tables():
    out = {}
    for name, key in (("sixlines_prym.jsonl", "p"), ("sixlines_prym_p41.jsonl", "p")):
        for r in load(name):
            out[(r[key], tuple(r["cfg"]))] = r
    return out


def prym_features(r):
    """Per configuration: histograms over its stored pairs (P family)."""
    pairs = [x for x in r["pairs"] if "c" in x]
    h = lambda key: tuple(sorted(collections.Counter(key(x) for x in pairs).items()))
    out = {
        "P hist c1": h(lambda x: tuple(x["c"][0])),
        "P hist c1 (rational a,b,c)": tuple(sorted(collections.Counter(tuple(x["c"][0]) for x in pairs
                                                                        if x.get("rational_point") == 0).items())),
        "P hist c2": h(lambda x: tuple(x["c"][1])),
        "P hist compare": h(lambda x: x.get("compare")),
        "P hist rational_point": h(lambda x: x.get("rational_point")),
        "P count c4_rational": sum(1 for x in pairs if x.get("c4_rational")),
        "P number of pairs": len(pairs),
    }
    c1c = collections.Counter(tuple(x["c"][0]) for x in pairs)
    for v in ((4, 4), (4, -4), (-4, 4), (-4, -4)):
        out[f"P parity of c1 = {v}"] = c1c.get(v, 0) % 2
    return out


def fourfold_records():
    """Per pair of lines, the fourfold data of notes 23 and 26 (primary readings only)."""
    recs = []
    for r in load("sixlines_jacobian_character.jsonl"):
        recs.append(r)
    for r in load("sixlines_jacobian_character_extra.jsonl"):
        if r.get("kind") == "more pairs":
            recs.append(r)
    for r in load("sixlines_delta41.jsonl"):
        if r.get("reading") == "primary":
            recs.append(r)
    out = []
    for r in recs:
        c1, c2 = r["c"][0], r["c"][1]
        ok = r["status"] == "ok"
        kern = tuple((k[0], tuple(k[1]), k[2]) for k in r.get("kernel", [])) if r.get("kernel") else None
        out.append({"p": r["p"], "cfg": tuple(r["cfg"]), "pair": tuple(r["pair"]), "F": {
            "F delta (control)": r.get("delta") if ok else None,
            "F u": r.get("u") if ok else None,
            "F M mod 16": (r["M"] % 16) if ok and r.get("M") is not None else None,
            "F v2_M": r.get("v2_M") if ok else None,
            "F order_X": r.get("order_X") if ok else None,
            "F kernel shape": str(kern) if kern is not None else None,
            "F status": r["status"],
            "F sign(Re c1 Im c1)": (c1[0] * c1[1] > 0) - (c1[0] * c1[1] < 0),
            "F sign(c2)": (c2[0] > 0) - (c2[0] < 0),
            "F tw": r.get("tw"), "F r": r.get("r"), "F s": r.get("s"),
            "F c1": str(tuple(c1)),
        }})
    return out


def cmd_features(args):
    import sixlines_cmmethod as cmm  # noqa: F401  (loads cm_tmat, cm_J into the PARI instance)
    import sixlines_characters as ch
    pari = cmm.pari
    rng = random.Random(100)
    L, L_small, L_partial, C = build_population()
    print("POPULATION")
    for tag, lst in (("L (full labels, |G| >= 3)", L), ("L, |G| = 2 (no information; excluded)", L_small),
                     ("L, only partly labelled (excluded)", L_partial), ("C (cosets of ker chi_-1)", C)):
        ns = [len(t["surfaces"]) for t in lst]
        nc = sum(len(s["configs"]) for t in lst for s in t["surfaces"])
        print(f"  {tag}: {len(lst)} torsors, {sum(ns)} surfaces, {nc} configurations; "
              + ", ".join(f"q={t['q']} e1={t['e'][0] if 'e' in t else t['e1']} G={t['G'] or 'coset'} n={len(t['surfaces'])}"
                          for t in lst))
    stored79 = {(r["q"], tuple(r["cfg"])): r for r in load("sixlines_characters.jsonl")}
    rows = {}
    for q in sorted({t["q"] for t in L + L_small + C}):
        for r in load(f"sixlines_p{q}.jsonl"):
            rows[(q, tuple(r["cfg"]))] = r
    quart = prym_tables()
    ffr = collections.defaultdict(list)
    for x in fourfold_records():
        ffr[(x["p"], x["cfg"])].append(x)
    checks = collections.Counter()
    tors_out = []
    for tor in L + L_small + C:
        q = tor["q"]
        e1 = tor["e"][0] if "e" in tor else tor["e1"]
        F = Fp(q)
        jf, ef = j_family(q), e_family(q)
        some = [c for s in tor["surfaces"] for c in s["configs"]][:3]
        bad = invariance_filter(pari, jf + ef, some, q, rng)
        checks[("monomials dropped by the invariance filter", len(bad))] += 1
        jf = [x for x in jf if x[0] not in bad]
        ef = [x for x in ef if x[0] not in bad]
        tsurf = []
        for s in tor["surfaces"]:
            cvals = []
            for cfg in s["configs"]:
                row = rows[(q, cfg)]
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
                v["S c1_values"] = str(sorted({tuple(x["c"][0]) for x in qr["pairs"]
                                               if x.get("rational_point") == 0 and "c" in x}) if qr else None)
                st = stored79.get((q, cfg))
                if st is not None:
                    for k in ("conic", "j2", "j4", "quartic", "orientation", "branch", "branch_zero", "type", "quartic_type"):
                        if k in st:
                            checks[("note 19 stored value reproduced", v["S " + k] == st[k])] += 1
                    checks[("note 19 c1_values reproduced", v["S c1_values"] == str(sorted(map(tuple, st["c1_values"]))))] += 1
                A = matrix_of(pari, cfg, q)
                Q = config_quantities(pari, A, q)
                ci = chart_invariants(Q["J"], q, tor["chart"])
                checks[("stored invariants of the surface reproduced from the configuration", ci == s["invariants"])] += 1
                for name, d, a in jf:
                    v[name] = mono_symbol([F.sym(x, d) for x in Q["J"]], a, d)
                for name, d, a in ef:
                    v[name] = mono_symbol([F.sym(x, d) for x in Q["p"] + [Q["C"]]], a, d)
                for d in (2, 4):
                    if (q - 1) % d == 0:
                        v[f"S root symbol d={d}"] = F.sym(s["root"], d)
                if qr:
                    v.update(prym_features(qr))
                v["_scalars"] = scalars(Q, q)
                cvals.append(v)
            # the surface value of each candidate: sorted pair over the configurations; 'agree' if they agree
            keys = sorted(set().union(*[set(x) for x in cvals]) - {"_scalars"})
            sv = {}
            for k in keys:
                vals = [json.dumps(x.get(k)) for x in cvals]
                sv[k] = {"pair": sorted(vals), "agree": len(set(vals)) == 1, "value": cvals[0].get(k)}
            sc = {}
            for k in cvals[0]["_scalars"]:
                vals = [x["_scalars"][k] for x in cvals]
                sc[k] = vals[0] if len(set(vals)) == 1 else "differs"
                checks[(f"scalar {k} the same on the configurations of a surface", len(set(vals)) == 1)] += 1
            sc["root"] = s["root"]
            out = {"configs": [list(c) for c in s["configs"]], "values": sv, "scalars": sc}
            for k in ("label", "coset", "delta"):
                if k in s:
                    out[k] = s[k]
            if tor["pop"] == "C":
                recs = [x for c in s["configs"] for x in ffr.get((q, c), [])]
                out["fourfold"] = [x["F"] for x in recs]
                checks[("C surfaces with fourfold records", bool(recs))] += 1
                dl = {x["F"]["F delta (control)"] for x in recs} - {None}
                checks[("C: delta of the records = delta of the surface", dl == {s["delta"]})] += 1
            tsurf.append(out)
        tors_out.append({"pop": tor["pop"], "q": q, "e1": e1, "G": tor["G"], "chart": tor["chart"],
                         "small": tor in L_small, "fam_sizes": {"J": len(jf), "E": len(ef)}, "surfaces": tsurf})
        print(f"  features: q={q} e1={e1} {tor['pop']} n={len(tsurf)}: J {len(jf)}, E {len(ef)} monomials kept "
              f"({len(bad)} dropped)", flush=True)
    print("CHECKS")
    for k, v in sorted(checks.items(), key=str):
        print(f"  {k}: {v}")
    with open(args.out, "w") as fh:
        json.dump({"torsors": tors_out}, fh)
    print(f"wrote {args.out}")


# ================================================================================================= audit

def group_elements(cyc):
    return list(itertools.product(*[range(n) for n in cyc]))


def sub(a, b, cyc):
    return tuple((x - y) % n for x, y, n in zip(a, b, cyc))


class TorsorL:
    """All n! relabellings of a fully labelled torsor, with the difference index of every ordered pair."""

    def __init__(self, tor):
        import numpy as np
        self.tor = tor
        self.cyc = tor["G"]
        self.elts = group_elements(self.cyc)
        self.n = len(self.elts)
        idx = {g: i for i, g in enumerate(self.elts)}
        self.lab = [idx[tuple(s["label"])] for s in tor["surfaces"]]
        diff = np.array([[idx[sub(b, a, self.cyc)] for b in self.elts] for a in self.elts], dtype=np.int16)  # diff[a, b] = b - a
        self.pairs = [(i, j) for i in range(self.n) for j in range(self.n) if i != j]
        perms = np.array(list(itertools.permutations(range(self.n))), dtype=np.int16)   # surface i -> label perms[k, i]
        I = np.array([i for i, j in self.pairs]), np.array([j for i, j in self.pairs])
        self.G = diff[perms[:, I[0]], perms[:, I[1]]]          # (K, N), values 1..n-1 (0 never: i != j)
        obs = np.array(self.lab, dtype=np.int16)
        self.obs_row = int(np.nonzero((perms == obs).all(axis=1))[0][0])
        self.K, self.N = self.G.shape


def mi_rows(G, f, ng):
    """N * MI(f; g) in bits for each row of G (K x N group indices < ng), f an int code vector of length N."""
    import numpy as np
    K, N = G.shape
    V = int(f.max()) + 1
    nf = np.bincount(f, minlength=V).astype(float)
    hf = -sum(c * math.log2(c / N) for c in nf if c)            # N * H(f)
    out = np.empty(K)
    step = max(1, 2_000_000 // max(1, N * 1))
    for s in range(0, K, step):
        Gs = G[s:s + step].astype(np.int64)
        k = Gs.shape[0]
        code = (np.arange(k)[:, None] * ng + Gs) * V + f[None, :]
        cnt = np.bincount(code.ravel(), minlength=k * ng * V).reshape(k, ng, V).astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(cnt > 0, cnt * np.log2(cnt), 0.0).sum(axis=2)          # sum_f c log c, per g
            ngc = cnt.sum(axis=2)
            u = np.where(ngc > 0, ngc * np.log2(ngc), 0.0)
        hfg = (u - t).sum(axis=1)                                 # N * H(f | g)
        out[s:s + k] = hf - hfg
    return out


def codes(values):
    """Integer codes for a list of hashable values."""
    import numpy as np
    m = {}
    return np.array([m.setdefault(json.dumps(v, sort_keys=True), len(m)) for v in values], dtype=np.int64)


def tail_p(dists):
    """dists: list of (array of per-row statistics, observed). Exact P(sum >= sum observed) on the grid H,
    conservative.  Each row of a torsor has the same probability."""
    import numpy as np
    acc = np.array([1.0])
    obs_bin = 0
    for vals, obs in dists:
        b = np.floor(vals / H + 1e-7).astype(np.int64)
        ob = int(math.floor(obs / H + 1e-7))
        obs_bin += ob
        u, c = np.unique(b, return_counts=True)
        pr = c / c.sum()
        new = np.zeros(len(acc) + int(u.max()))
        for ui, pi in zip(u, pr):
            new[ui:ui + len(acc)] += pi * acc
        acc = new
    thr = obs_bin - len(dists) + 1
    thr = max(thr, 0)
    return float(acc[thr:].sum()) if thr < len(acc) else 0.0


def poisson_binomial_tail(ps, k):
    dp = [1.0]
    for p in ps:
        nd = [0.0] * (len(dp) + 1)
        for i, x in enumerate(dp):
            nd[i] += x * (1 - p)
            nd[i + 1] += x * p
        dp = nd
    return sum(dp[k:])


# -------------------------------------------------------------- candidate construction, population L

ATTEMPTED = set()


def pad_constant(res, pop):
    """Candidates that were formed but are constant on every torsor enter the correction with p = 1."""
    have = {r["name"] for r in res}
    extra = sorted(ATTEMPTED - have)
    for k in extra:
        res.append({"name": k, "p": 1.0, "stat": 0.0, "torsors": 0, "exact_count": 0, "p_exact_count": 1.0,
                    "detail": {}, "constant": True})
    ATTEMPTED.clear()
    return len(extra)


def l_candidates(tors):
    """For population L: dict name -> {torsor index: pair code vector (length N)}.  Variants: 'diff' for symbol
    candidates (f(Y) - f(X) in Z/d), 'pair' (f(X), f(Y)) for all surface candidates, and the B2 relational ones."""
    import numpy as np
    cands = collections.defaultdict(dict)
    for ti, T in enumerate(tors):
        surf = T.tor["surfaces"]
        q = T.tor["q"]
        keys = sorted(set().union(*[set(s["values"]) for s in surf]))
        for k in keys:
            ATTEMPTED.add(k + " | pair")
            if symbol_order(k, None):
                ATTEMPTED.add(k + " | diff")
            vals = [s["values"].get(k) for s in surf]
            if any(v is None for v in vals):
                continue
            pairv = [json.dumps(v["pair"]) for v in vals]
            if len(set(pairv)) > 1:
                cands[k + " | pair"][ti] = codes([[pairv[i], pairv[j]] for i, j in T.pairs])
            d = symbol_order(k, vals)
            if d and all(v["agree"] for v in vals):
                x = [v["value"] for v in vals]
                if len(set(map(str, x))) > 1:
                    cands[k + " | diff"][ti] = codes([("z" if "z" in (x[i], x[j]) else (to_int(x[j]) - to_int(x[i])) % d)
                                                      for i, j in T.pairs])
        Fq = FpSmall(q)
        for k in surf[0]["scalars"]:
            sv = [s["scalars"][k] for s in surf]
            if any(v == "differs" for v in sv):
                continue
            for d in (2, 4):
                if (q - 1) % d:
                    continue
                for op, fn in (("-", lambda a, b: a - b), ("+", lambda a, b: a + b)):
                    ATTEMPTED.add(f"B2 sym d={d} of {k}(Y) {op} {k}(X)")
                    vals = [("u" if sv[i] is None or sv[j] is None else Fq.sym(fn(sv[j], sv[i]), d)) for i, j in T.pairs]
                    if len(set(map(str, vals))) > 1:
                        cands[f"B2 sym d={d} of {k}(Y) {op} {k}(X)"][ti] = codes(vals)
    return cands


def to_int(x):
    return int(x)


def symbol_order(k, vals):
    """The order d if candidate k is group-valued (residue symbols, ±1 bits), else 0."""
    if k.startswith("J d=") or k.startswith("E d="):
        return int(k.split("d=")[1].split()[0])
    if k.startswith("S root symbol d="):
        return int(k.split("d=")[1])
    if k in ("S conic", "S j2", "S j4", "S quartic") or k.startswith("P parity"):
        return 2
    return 0


class FpSmall:
    def __init__(self, p):
        self.F = Fp(p)

    def sym(self, x, d):
        return self.F.sym(x, d)


def bits_to_z2(v):
    """±1 -> 0/1; 0 -> 'z'."""
    return "z" if v == 0 else (0 if v == 1 else 1)


def prepare_l(tors):
    """Bits stored as ±1 are mapped to Z/2 before forming differences."""
    for T in tors:
        for s in T.tor["surfaces"]:
            for k in ("S conic", "S j2", "S j4", "S quartic"):
                if k in s["values"]:
                    v = s["values"][k]
                    v["value"] = bits_to_z2(v["value"])


# -------------------------------------------------------------- workers

_TORS = None


def _init_l(tors_json):
    global _TORS
    _TORS = [TorsorL(t) for t in tors_json]


def _score_l(job):
    """job: (name, {ti: code list}) -> result dict."""
    import numpy as np
    name, per = job
    dists, ex_p, ex_obs, detail = [], [], 0, {}
    for ti, f in per.items():
        T = _TORS[ti]
        f = np.asarray(f, dtype=np.int64)
        vals = mi_rows(T.G, f, T.n)
        obs = vals[T.obs_row]
        dists.append((vals, obs))
        # exact function of g (non-constant): N * MI = N * H(f)
        Nf = np.bincount(f).astype(float)
        hfull = -sum(c * math.log2(c / len(f)) for c in Nf if c)
        isfun = np.abs(vals - hfull) < 1e-9
        ex_p.append(float(isfun.mean()))
        ex_obs += bool(isfun[T.obs_row])
        detail[ti] = {"stat": float(obs), "p_torsor": float((vals >= obs - 1e-9).mean()), "function_of_g": bool(isfun[T.obs_row])}
    p = tail_p(dists)
    pf = poisson_binomial_tail(ex_p, ex_obs) if ex_obs else 1.0
    return {"name": name, "p": p, "stat": sum(o for _, o in dists), "torsors": len(per), "exact_count": ex_obs,
            "p_exact_count": pf, "detail": detail}


# -------------------------------------------------------------- population C

class TorsorC:
    """All balanced splits of a coset-labelled torsor (each split once, with its complement)."""

    def __init__(self, tor):
        import numpy as np
        self.tor = tor
        self.n = len(tor["surfaces"])
        self.obs = np.array([s["coset"] for s in tor["surfaces"]], dtype=np.int64)
        rows = []
        for S in itertools.combinations(range(self.n), self.n // 2):
            r = np.zeros(self.n, dtype=np.int64)
            r[list(S)] = 1
            rows.append(r)
        self.splits = np.array(rows)                  # (K, n): coset bit of each surface
        self.obs_row = int(np.nonzero((self.splits == self.obs).all(axis=1) | (self.splits == 1 - self.obs).all(axis=1))[0][0])
        self.pairs = [(i, j) for i in range(self.n) for j in range(self.n) if i != j]


def c_candidates(tors):
    """name -> {ti: ('surface', codes over surfaces) | ('pair', codes over ordered pairs) | ('record', codes, surface index)}"""
    cands = collections.defaultdict(dict)
    for ti, T in enumerate(tors):
        surf = T.tor["surfaces"]
        q = T.tor["q"]
        keys = sorted(set().union(*[set(s["values"]) for s in surf]))
        for k in keys:
            ATTEMPTED.add(k)
            vals = [s["values"].get(k) for s in surf]
            if any(v is None for v in vals):
                continue
            pv = [json.dumps(v["pair"]) for v in vals]
            if len(set(pv)) > 1:
                cands[k][ti] = ("surface", list(codes(pv)))
        dv = [s["delta"] for s in surf]
        cands["S delta (positive control)"][ti] = ("surface", list(codes(dv)))
        Fq = FpSmall(q)
        for k in surf[0]["scalars"]:
            sv = [s["scalars"][k] for s in surf]
            if any(v == "differs" for v in sv):
                continue
            for d in (2, 4, 8):
                if (q - 1) % d:
                    continue
                for op, fn in (("-", lambda a, b: a - b), ("+", lambda a, b: a + b)):
                    ATTEMPTED.add(f"B2 sym d={d} of {k}(Y) {op} {k}(X)")
                    vals = [("u" if sv[i] is None or sv[j] is None else Fq.sym(fn(sv[j], sv[i]), d)) for i, j in T.pairs]
                    if len(set(map(str, vals))) > 1:
                        cands[f"B2 sym d={d} of {k}(Y) {op} {k}(X)"][ti] = ("pair", list(codes(vals)))
        recs = [(si, r) for si, s in enumerate(surf) for r in s["fourfold"]]
        if recs:
            for k in recs[0][1]:
                ATTEMPTED.add(k)
                sub_ = [(si, r[k]) for si, r in recs if r.get(k) is not None]
                if len({json.dumps(v) for _, v in sub_}) > 1:
                    cands[k][ti] = ("record", list(codes([v for _, v in sub_])), [si for si, _ in sub_])
    return cands


_TC = None


def _init_c(tors_json):
    global _TC
    _TC = [TorsorC(t) for t in tors_json]


def _score_c(job):
    import numpy as np
    name, per = job
    dists, detail = [], {}
    for ti, spec in per.items():
        T = _TC[ti]
        kind = spec[0]
        if kind == "surface":
            f = np.asarray(spec[1], dtype=np.int64)
            Gm = T.splits                                         # (K, n)
        elif kind == "pair":
            f = np.asarray(spec[1], dtype=np.int64)
            I = np.array([i for i, j in T.pairs]), np.array([j for i, j in T.pairs])
            Gm = (T.splits[:, I[0]] != T.splits[:, I[1]]).astype(np.int64)
        else:
            f = np.asarray(spec[1], dtype=np.int64)
            Gm = T.splits[:, np.asarray(spec[2])]
        vals = mi_rows(Gm, f, 2)
        obs = vals[T.obs_row]
        dists.append((vals, obs))
        detail[ti] = {"stat": float(obs), "p_torsor": float((vals >= obs - 1e-9).mean())}
    return {"name": name, "p": tail_p(dists), "stat": sum(o for _, o in dists), "torsors": len(per), "detail": detail}


# -------------------------------------------------------------- corrections

def holm(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj = [0.0] * m
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * ps[i]))
        adj[i] = run
    return adj


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i], reverse=True)
    adj = [0.0] * m
    run = 1.0
    for r, i in enumerate(order):
        rank = m - r
        run = min(run, ps[i] * m / rank)
        adj[i] = run
    return adj


# -------------------------------------------------------------- planted signals (power), population L

def plants(tors, kind, rng, flip=0.0, frac=1.0):
    """A planted candidate: {ti: pair codes}."""
    import numpy as np
    per = {}
    chosen = list(range(len(tors)))
    if frac < 1.0:
        chosen = rng.sample(chosen, max(1, round(frac * len(tors))))
    for ti, T in enumerate(tors):
        cyc, elts = T.cyc, T.elts
        lab = [elts[i] for i in T.lab]
        if kind in ("chi2", "chi4"):
            d = 2 if kind == "chi2" else 4
            # characters of order exactly d: chi(x) = sum c_k x_k (n/ cyc_k) ... as a hom to Z/d
            # homomorphisms G -> Z/d of order exactly d: x -> sum a_k x_k with n_k a_k = 0 mod d
            homs = []
            for a in itertools.product(range(d), repeat=len(cyc)):
                if any((nk * ak) % d for nk, ak in zip(cyc, a)):
                    continue
                vals = {x: sum(ak * xk for ak, xk in zip(a, x)) % d for x in elts}
                if len(set(vals.values())) == d:
                    homs.append(vals)
            if not homs:
                continue
            if ti not in chosen:
                f = [rng.randrange(d) for _ in lab]
            else:
                chi = rng.choice(homs)
                f = [chi[x] if rng.random() >= flip else rng.randrange(d) for x in lab]
            per[ti] = list(codes([(f[j] - f[i]) % d for i, j in T.pairs]))
        elif kind == "random_pm1":
            f = [rng.randrange(2) for _ in lab]
            if len(set(f)) == 1:
                continue
            per[ti] = list(codes([(f[j] - f[i]) % 2 for i, j in T.pairs]))
        elif kind == "h_of_g":
            h = {g: rng.randrange(2) for g in elts}
            per[ti] = list(codes([h[sub(lab[j], lab[i], cyc)] for i, j in T.pairs]))
    return per


# -------------------------------------------------------------- main audit

def cmd_audit(args):
    from multiprocessing import get_context
    with open(args.features) as fh:
        feat = json.load(fh)
    torsL = [t for t in feat["torsors"] if t["pop"] == "L" and not t["small"]]
    torsC = [t for t in feat["torsors"] if t["pop"] == "C"]
    out = open(args.out, "w")
    print(f"population L: {len(torsL)} torsors, {sum(len(t['surfaces']) for t in torsL)} surfaces; "
          f"population C: {len(torsC)} torsors, {sum(len(t['surfaces']) for t in torsC)} surfaces")
    objsL = [TorsorL(t) for t in torsL]
    prepare_l(objsL)
    candL = l_candidates(objsL)
    jobsL = [(k, {ti: list(map(int, v)) for ti, v in per.items()}) for k, per in sorted(candL.items())]
    print(f"L: {len(jobsL)} candidates (pair functions)")
    ctx = get_context("fork")
    with ctx.Pool(args.workers, initializer=_init_l, initargs=(torsL,)) as pool:
        resL = pool.map(_score_l, jobsL, chunksize=4)
        # power and null controls, scored against the same exact null; m = number of real candidates
        rng = random.Random(2026)
        ctl_jobs = []
        for kind, flip, frac, tag in (("chi2", 0.0, 1.0, "P1 order-2 character, every torsor"),
                                      ("chi2", 0.125, 1.0, "P2 order-2 character, values flipped w.p. 1/8"),
                                      ("chi2", 0.25, 1.0, "P2' order-2 character, values flipped w.p. 1/4"),
                                      ("chi2", 0.0, 0.5, "P3 order-2 character on half the torsors"),
                                      ("chi2", 0.0, 0.25, "P3' order-2 character on a quarter of the torsors"),
                                      ("chi4", 0.0, 1.0, "P4 order-4 character where G has one"),
                                      ("h_of_g", 0.0, 1.0, "P5 a random function of g"),
                                      ("random_pm1", 0.0, 1.0, "N random ±1 function (null)")):
            for r in range(args.plants):
                per = plants(objsL, kind, rng, flip, frac)
                if per:
                    ctl_jobs.append((f"{tag} #{r}", per))
        resCtl = pool.map(_score_l, ctl_jobs, chunksize=4)
    nconstL = pad_constant(resL, "L")
    print(f"L: {nconstL} further candidates constant on every torsor, entered with p = 1")
    mL = len(resL)
    ps = [r["p"] for r in resL]
    for r, a, b in zip(resL, holm(ps), bh(ps)):
        r["holm"], r["bh"], r["bonferroni"] = a, b, min(1.0, mL * r["p"])
        r["pop"] = "L"
        out.write(json.dumps(r) + "\n")
    # C
    objsC = [TorsorC(t) for t in torsC]
    ATTEMPTED.clear()
    candC = c_candidates(objsC)
    jobsC = [(k, per) for k, per in sorted(candC.items())]
    print(f"C: {len(jobsC)} candidates")
    with ctx.Pool(args.workers, initializer=_init_c, initargs=(torsC,)) as pool:
        resC = pool.map(_score_c, jobsC, chunksize=4)
    nconstC = pad_constant(resC, "C")
    print(f"C: {nconstC} further candidates constant on every torsor, entered with p = 1")
    mC = len(resC)
    ps = [r["p"] for r in resC]
    for r, a, b in zip(resC, holm(ps), bh(ps)):
        r["holm"], r["bh"], r["bonferroni"] = a, b, min(1.0, mC * r["p"])
        r["pop"] = "C"
        out.write(json.dumps(r) + "\n")
    for r in resCtl:
        r["pop"] = "control L"
        r["bonferroni"] = min(1.0, mL * r["p"])
        out.write(json.dumps(r) + "\n")
    out.close()
    report(torsL, torsC, resL, resC, resCtl, mL, mC)


def family_of(name):
    if name.startswith("B2"):
        return "B2 relational symbols"
    if name.startswith("F "):
        return "F fourfold per pair of lines"
    head = name.split()[0]
    return {"S": "S note 19 invariants and root", "J": "J residue symbols of J-monomials",
            "E": "E residue symbols of power sums and conic", "P": "P curve histograms"}.get(head, head)


def report(torsL, torsC, resL, resC, resCtl, mL, mC):
    print("\nTORSORS")
    for t in torsL:
        print(f"  L q={t['q']} e1={t['e1']} G={t['G']} n={len(t['surfaces'])}")
    for t in torsC:
        print(f"  C q={t['q']} e1={t['e1']} n={len(t['surfaces'])}")
    for tag, res, m in (("L", resL, mL), ("C", resC, mC)):
        print(f"\nPOPULATION {tag}: m = {m} candidates")
        fam = collections.defaultdict(list)
        for r in res:
            v = r["name"].split(" | ")[-1] if " | " in r["name"] else ""
            fam[(family_of(r["name"]), v)].append(r)
        print("  family | variant | candidates | best candidate | stat | exact p | Holm | BH q")
        for (f, v), rs in sorted(fam.items()):
            b = min(rs, key=lambda r: r["p"])
            print(f"  {f} | {v} | {len(rs)} | {b['name']} | {b['stat']:.2f} | {b['p']:.3g} | {b['holm']:.3g} | {b['bh']:.3g}")
        print("  the 12 smallest p-values:")
        for r in sorted(res, key=lambda r: r["p"])[:12]:
            extra = f" exact-function torsors {r['exact_count']}/{r['torsors']} (p {r['p_exact_count']:.3g})" if "exact_count" in r else ""
            print(f"    {r['name']}: stat {r['stat']:.2f} on {r['torsors']} torsors, p {r['p']:.3g}, Holm {r['holm']:.3g}, "
                  f"BH {r['bh']:.3g};{extra}; per torsor p "
                  + ", ".join(f"{d['p_torsor']:.3g}" for d in r["detail"].values()))
        if tag == "L":
            ex = sorted((r for r in res if r["exact_count"]), key=lambda r: r["p_exact_count"])[:8]
            print("  secondary statistic, exact non-constant functions of g; smallest p:")
            for r in ex:
                print(f"    {r['name']}: {r['exact_count']} of {r['torsors']} torsors, p {r['p_exact_count']:.3g}, "
                      f"Bonferroni {min(1, m * r['p_exact_count']):.3g}; torsors "
                      + ", ".join(str(ti) for ti, d in r["detail"].items() if d["function_of_g"]))
    print(f"\nCONTROLS on L (Bonferroni at m = {mL}, alpha = 0.05)")
    grp = collections.defaultdict(list)
    for r in resCtl:
        grp[r["name"].split(" #")[0]].append(r)
    for k, rs in grp.items():
        det = sum(1 for r in rs if r["p"] * mL <= 0.05)
        unc = sum(1 for r in rs if r["p"] <= 0.05)
        med = sorted(r["p"] for r in rs)[len(rs) // 2]
        print(f"  {k}: {len(rs)} plants, detected after correction {det}, p <= 0.05 uncorrected {unc}, median p {med:.3g}")


def plants_c(tors, rng, flip=0.0, frac=1.0, kind="coset"):
    """A planted surface function on C: the coset bit with flips, on a fraction of the torsors; or random."""
    per = {}
    chosen = set(rng.sample(range(len(tors)), max(1, round(frac * len(tors)))))
    for ti, T in enumerate(tors):
        if kind == "random":
            f = [rng.randrange(2) for _ in range(T.n)]
        elif ti in chosen:
            f = [int(b) if rng.random() >= flip else rng.randrange(2) for b in T.obs]
        else:
            f = [rng.randrange(2) for _ in range(T.n)]
        if len(set(f)) > 1:
            per[ti] = ("surface", f)
    return per


def cmd_extra(args):
    """Shuffled-label reruns of population L (the distribution of the smallest p over all candidates when the
    labels carry no information), and planted signals on population C."""
    from multiprocessing import get_context
    with open(args.features) as fh:
        feat = json.load(fh)
    torsL = [t for t in feat["torsors"] if t["pop"] == "L" and not t["small"]]
    torsC = [t for t in feat["torsors"] if t["pop"] == "C"]
    ctx = get_context("fork")
    out = open(args.out, "w")
    print("SHUFFLED LABELS, population L: each run permutes the labels of every torsor at random and repeats the "
          "whole family of candidates")
    for seed in range(1, args.shuffles + 1):
        rng = random.Random(seed)
        sh = json.loads(json.dumps(torsL))
        for t in sh:
            labs = [s["label"] for s in t["surfaces"]]
            rng.shuffle(labs)
            for s, l in zip(t["surfaces"], labs):
                s["label"] = l
        objs = [TorsorL(t) for t in sh]
        prepare_l(objs)
        ATTEMPTED.clear()
        cand = l_candidates(objs)
        jobs = [(k, {ti: list(map(int, v)) for ti, v in per.items()}) for k, per in sorted(cand.items())]
        with ctx.Pool(args.workers, initializer=_init_l, initargs=(sh,)) as pool:
            res = pool.map(_score_l, jobs, chunksize=4)
        ps = sorted(r["p"] for r in res)
        rec = {"seed": seed, "candidates": len(res), "min_p": ps[0], "p<=0.05": sum(p <= 0.05 for p in ps),
               "p<=0.01": sum(p <= 0.01 for p in ps), "best": min(res, key=lambda r: r["p"])["name"]}
        out.write(json.dumps(rec) + "\n")
        print(f"  seed {seed}: {rec['candidates']} candidates, smallest p {rec['min_p']:.3g} ({rec['best']}), "
              f"p <= 0.05: {rec['p<=0.05']}, p <= 0.01: {rec['p<=0.01']}", flush=True)
    print("\nPLANTED SIGNALS, population C (Bonferroni at the m of the audit, " + str(args.mC) + ")")
    objsC = [TorsorC(t) for t in torsC]
    rng = random.Random(2027)
    jobs = []
    for flip, frac, tag in ((0.0, 1.0, "the coset bit itself"), (0.125, 1.0, "coset bit, values flipped w.p. 1/8"),
                            (0.25, 1.0, "coset bit, values flipped w.p. 1/4"), (0.0, 0.4, "coset bit on 2 of 5 torsors"),
                            (0.0, 0.2, "coset bit on 1 of 5 torsors"), (None, 1.0, "random ±1 function (null)")):
        for r in range(args.plants):
            per = plants_c(objsC, rng, flip or 0.0, frac, "random" if flip is None else "coset")
            jobs.append((f"{tag} #{r}", per))
    with ctx.Pool(args.workers, initializer=_init_c, initargs=(torsC,)) as pool:
        res = pool.map(_score_c, jobs, chunksize=4)
    grp = collections.defaultdict(list)
    for r in res:
        grp[r["name"].split(" #")[0]].append(r)
        out.write(json.dumps({"plant": r["name"], "p": r["p"]}) + "\n")
    for k, rs in grp.items():
        det = sum(1 for r in rs if r["p"] * args.mC <= 0.05)
        unc = sum(1 for r in rs if r["p"] <= 0.05)
        med = sorted(r["p"] for r in rs)[len(rs) // 2]
        print(f"  {k}: {len(rs)} plants, detected after correction {det}, p <= 0.05 uncorrected {unc}, median p {med:.3g}")
    out.close()


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("features")
    a.add_argument("--out", default=os.path.join(HERE, "sixlines_leakage_features.json"))
    b = sp.add_parser("audit")
    b.add_argument("--features", default=os.path.join(HERE, "sixlines_leakage_features.json"))
    b.add_argument("--out", default=os.path.join(HERE, "sixlines_leakage.jsonl"))
    b.add_argument("--workers", type=int, default=12)
    b.add_argument("--plants", type=int, default=200)
    c = sp.add_parser("extra")
    c.add_argument("--features", default=os.path.join(HERE, "sixlines_leakage_features.json"))
    c.add_argument("--out", default=os.path.join(HERE, "sixlines_leakage_extra.jsonl"))
    c.add_argument("--workers", type=int, default=12)
    c.add_argument("--shuffles", type=int, default=10)
    c.add_argument("--plants", type=int, default=200)
    c.add_argument("--mC", type=int, default=6579)
    args = ap.parse_args()
    if args.cmd == "features":
        cmd_features(args)
    elif args.cmd == "audit":
        cmd_audit(args)
    else:
        cmd_extra(args)


if __name__ == "__main__":
    main()
