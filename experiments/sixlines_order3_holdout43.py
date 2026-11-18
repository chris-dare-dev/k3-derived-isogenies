"""sixlines_order3_holdout43.py -- held-out test (note 29 section 7.4) of the A43 survivor
"B2 sym d=2 of R2(Y) + R2(X)" on q = 43, e1 = 2 (Z/24), where the audit does not score B2.  The test as written
down beforehand is in sixlines_order3_holdout43_prereg.txt (T1 over pairs i < j, T2); the last block, added after
that run, is the survivor's own definition over ordered pairs.  Run from experiments/ with
python -I sixlines_order3_holdout43.py."""
import json, itertools
import numpy as np
feat = json.load(open("sixlines_order3_features_q43.json"))
t = [t for t in feat["torsors"] if t["q"] == 43 and t["e1"] == 2][0]
q, surf = 43, t["surfaces"]
n = len(surf); assert n == 24
c = np.array([s["coset3"] for s in surf])
R2 = [s["scalars"]["R2"] for s in surf]
P = list(itertools.combinations(range(n), 2))
def sym(x):
    if x % q == 0: return 2
    return 0 if pow(x, (q - 1) // 2, q) == 1 else 1
f = np.array([3 if R2[i] is None or R2[j] is None else sym(R2[i] + R2[j]) for i, j in P])
I = np.array([i for i, j in P]); J = np.array([j for i, j in P])
def stats(lab):                      # lab: (B, n)
    g = (lab[:, J] - lab[:, I]) % 3      # (B, pairs)
    B = lab.shape[0]
    cnt = np.zeros((B, 4, 3))
    for a in range(4):
        for b in range(3):
            cnt[:, a, b] = ((f == a)[None, :] & (g == b)).sum(1)
    N = cnt.sum((1, 2), keepdims=True)
    pab = cnt / N; pa = pab.sum(2, keepdims=True); pb = pab.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        mi = np.nansum(np.where(pab > 0, pab * np.log2(pab / (pa * pb)), 0), (1, 2))
    same = (g == 0)
    sq = ((f == 0)[None, :] & same).sum(1); ns = ((f == 1)[None, :] & same).sum(1)
    return mi, ns / np.maximum(sq + ns, 1), sq, ns
obs = stats(c[None, :])
print("symbol counts over pairs (sq, nonsq, zero, undef):", np.bincount(f, minlength=4).tolist())
print("observed: MI %.5f bits; same-coset non-square share %.3f (squares %d, non-squares %d)" % (obs[0][0], obs[1][0], obs[2][0], obs[3][0]))
rng = np.random.default_rng(4302)
base = np.repeat([0, 1, 2], 8)
D = 200000; ge1 = ge2 = 0; m2 = []
for k in range(0, D, 20000):
    lab = np.array([rng.permutation(base) for _ in range(20000)])
    mi, sh, _, _ = stats(lab)
    ge1 += int((mi >= obs[0][0] - 1e-12).sum()); ge2 += int((sh >= obs[1][0] - 1e-12).sum()); m2.append(sh.mean())
print("T1 p = %.4g; T2 p = %.4g (null mean share %.3f)" % ((1 + ge1) / (1 + D), (1 + ge2) / (1 + D), np.mean(m2)))

# --- after the pre-registered run: the table, and the survivor's own definition (ordered pairs i != j)
g = (c[J] - c[I]) % 3
tab = np.zeros((4, 3), int)
for a, b in zip(f, g): tab[a, b] += 1
print("table f (sq, nonsq, zero, undef) x g (0, 1, 2) over i < j:\n", tab)
P2 = [(i, j) for i in range(n) for j in range(n) if i != j]
I = np.array([i for i, j in P2]); J = np.array([j for i, j in P2])
f = np.array([3 if R2[i] is None or R2[j] is None else sym(R2[i] + R2[j]) for i, j in P2])
obs = stats(c[None, :])
rng = np.random.default_rng(4302); ge = 0
for k in range(0, D, 20000):
    lab = np.array([rng.permutation(base) for _ in range(20000)])
    ge += int((stats(lab)[0] >= obs[0][0] - 1e-12).sum())
print("ordered pairs (the B2 definition): MI %.5f bits, p = %.4g" % (obs[0][0], (1 + ge) / (1 + D)))
# same-vs-different only, i<j
