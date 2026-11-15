"""Section 8 (Theorem 8.2, Corollary 8.3, Remark 8.4) and Remark 10.2 (open questions)."""
from common import check, finish

import math
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

from starforest import adjacency, roots, energy_secular, partitions


def rho(n, w, part, h=1.0):
    """Spectral radius of Gamma^h_pi via the largest root of Phi_h = 1 (Theorem 4.1(ii))."""
    return roots(n - w - sum(part) - len(part), w, part, h)[0]


def admissible(n, w, M):
    return [p for p in partitions(M) if n - w - M - len(p) >= 0]


def dominates(p, q):
    a, b = list(p) + [0] * 64, list(q) + [0] * 64
    return all(sum(a[:k]) >= sum(b[:k]) for k in range(1, 65))


print("Theorem 8.2: rho strictly increasing in the dominance order")
pairs = viol = 0
for h in (0.2, 0.5, 1.0, 2.5, 6.0):
    for M in range(2, 8):
        for w in range(0, 4):
            for n in range(w + M + 2, w + 2 * M + 3):
                P = admissible(n, w, M)
                val = {p: rho(n, w, p, h) for p in P}
                for p in P:
                    for q in P:
                        if p != q and dominates(p, q):
                            pairs += 1
                            viol += not val[p] > val[q] + 1e-10
print(f"  ({pairs} comparable pairs, M <= 7, w <= 3, h in {{0.2,0.5,1,2.5,6}})")
check("violations", viol, 0)

print("Corollary 8.3: rho(Gamma_(1^M)) = (n - 1 + sqrt((n+1)^2 - 4(2M + w)))/2")
worst = 0.0
for M in range(1, 6):
    for w in range(0, 4):
        for n in range(w + 2 * M, w + 2 * M + 4):
            exact = (n - 1 + math.sqrt((n + 1) ** 2 - 4 * (2 * M + w))) / 2
            A, _ = adjacency(n - w - 2 * M, w, (1,) * M)
            worst = max(worst, abs(np.linalg.eigvalsh(A)[-1] - exact))
check("largest deviation", f"{worst:.1e}", "0", ok=worst < 1e-10)

print("Remark 8.4: walk counts of the switched-off forest")


def mwc(ms, k):
    return sum((m + 1) * m ** (k // 2) if k % 2 == 0 else 2 * m ** ((k + 1) // 2) for m in ms)


ok = True
for ms in [(3, 2), (4, 1, 1), (5,), (2, 2, 2)]:
    A, _ = adjacency(0, 0, ms, h=0.0)
    AF = np.ones_like(A) - np.eye(len(A)) - A          # adjacency matrix of F itself
    one = np.ones(len(A))
    for k in range(0, 7):
        direct = one @ np.linalg.matrix_power(AF, k) @ one
        formula = len(A) if k == 0 else mwc(ms, k)
        ok &= abs(direct - formula) < 1e-9
check("mwc_k(F) = sum_l f_k(m_l) for k <= 6", ok, True)
inc = True
for M in range(2, 9):
    for p in partitions(M):
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                a, b = p[i], p[j]
                q = list(p)
                q[i], q[j] = a + 1, b - 1
                q = [v for v in q if v > 0]
                for k in range(2, 8):
                    inc &= mwc(q, k) > mwc(p, k)
check("a transfer (a,b) -> (a+1,b-1) strictly increases mwc_k for 2 <= k <= 7, M <= 8",
      inc, True)

print("Remark 10.2 (i)")
inc_pairs = {}
for M in range(1, 8):
    P = list(partitions(M))
    inc_pairs[M] = [(p, q) for i, p in enumerate(P) for q in P[i + 1:]
                    if not dominates(p, q) and not dominates(q, p)]
first = min(M for M in inc_pairs if inc_pairs[M])
check("smallest M with an incomparable pair", first, 6)
check("incomparable pairs at M = 6", inc_pairs[6],
      "[(4,1,1),(3,3)] and [(3,1,1,1),(2,2,2)]",
      ok=sorted(inc_pairs[6]) == sorted([((4, 1, 1), (3, 3)), ((3, 1, 1, 1), (2, 2, 2))]))
check("rho at n=10, w=0: (3,1,1,1) vs (2,2,2)",
      (round(rho(10, 0, (3, 1, 1, 1)), 4), round(rho(10, 0, (2, 2, 2)), 4)), (8.8446, 8.8437))
print(f"  (for comparison, (4,1,1) vs (3,3) at n=10, w=0: "
      f"{rho(10, 0, (4, 1, 1)):.4f} vs {rho(10, 0, (3, 3)):.4f})")
gap, ties = math.inf, 0
for M in range(2, 9):
    for w in range(0, 4):
        for n in range(w + M + 2, w + 2 * M + 3):
            vals = sorted(rho(n, w, p) for p in admissible(n, w, M))
            for a, b in zip(vals, vals[1:]):
                gap = min(gap, b - a)
                ties += b - a < 1e-9
print(f"  (M <= 8, w <= 3, w+M+2 <= n <= w+2M+2: smallest gap between radii {gap:.2e})")
check("partitions with equal spectral radius", ties, 0)

print("Remark 10.2 (ii): energy is not monotone in the dominance order")
for (n, w, M, claim) in [(14, 3, 6, ("(2,2,2)", "(2,1,1,1,1)")), (11, 0, 4, ("(4)", None))]:
    vals = {p: energy_secular(n - w - M - len(p), w, p) for p in admissible(n, w, M)}
    hi = max(vals, key=vals.get)
    lo = min(vals, key=vals.get)
    fmt = lambda p: "(" + ",".join(map(str, p)) + ")"
    if claim[1] is None:
        check(f"n={n}, w={w}, M={M}: largest at", fmt(hi), claim[0])
    else:
        check(f"n={n}, w={w}, M={M}: largest, smallest at", (fmt(hi), fmt(lo)), claim)

print("Remark 10.2 (iii): connected C with coronal_{-A(C)}(c) >= 1 on (0,1]")
grid = np.linspace(1e-4, 1.0, 4001)


def coronal_ok(G):
    B = nx.to_numpy_array(G)
    lam, V = np.linalg.eigh(-B)
    a = (V.T @ np.ones(len(B))) ** 2
    for cval in grid:
        d = cval - lam
        if np.any(np.abs(d) < 1e-9):
            continue
        if (a / d).sum() < 1 - 1e-9:
            return False
    return True


graphs = [G for G in graph_atlas_g()[1:] if 2 <= len(G) <= 6 and nx.is_connected(G)]
good = [G for G in graphs if coronal_ok(G)]
stars = [G for G in good if nx.is_tree(G) and max(dict(G.degree()).values()) == len(G) - 1]
bad_trees = [G for G in graphs if nx.is_tree(G) and G not in good]
print(f"  connected graphs on 2..6 vertices: {len(graphs)}, with the property: {len(good)}, "
      f"of which stars: {len(stars)}")
check("P_4 fails", coronal_ok(nx.path_graph(4)), False)
check("all stars K_{1,m} (m <= 5) have the property",
      all(coronal_ok(nx.star_graph(m)) for m in range(1, 6)), True)
check("the class is larger than the stars", len(good) > len(stars), True)
check("not all trees have it (smallest failing tree has 4 vertices)",
      min(len(G) for G in bad_trees), 4)

finish()
