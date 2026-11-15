"""Exhaustive checks over all members of small order (h = 1).

Theorem 4.2 (inertia), Corollary 4.3 (HL-index), Theorem 5.2 (forbidden interval),
Corollary 5.3 (energy as a positive part),
Theorem 6.1 with Remark 6.2 and Example 6.3, the dense-versus-secular comparison of
Section 5.1, and Remark 10.1 (no cospectral pair for n <= 9).
"""
from common import check, finish

import numpy as np

from starforest import adjacency, members, roots, energy_secular, energy_dense
from starforest.bounds import (invariants, lower_bound, km_bound, mcclelland_analogue,
                               mcclelland_gutman)

TOL = 1e-9

# ---------------------------------------------------------------- n <= 13
fam13 = list(members(3, 13))
check("number of members with 3 <= n <= 13", len(fam13), 1157)

bad = dict(inertia=0, gap=0, Epos=0, HL=0, bounds=0, km_hyp=0)
km_eq, km_vs_mcc, analog_gt, ratios_km, ratios_an = [], [], 0, [], []
for (n0, w, ms) in fam13:
    A, _ = adjacency(n0, w, ms)
    lam = np.sort(np.linalg.eigvalsh(A))[::-1]
    p = invariants(n0, w, ms)
    n, r, c = p["n"], p["r"], p["c"]
    E = float(np.abs(lam - c).sum())

    # Theorem 4.2
    npos, nneg = int((lam > TOL).sum()), int((lam < -TOL).sum())
    if (npos, nneg) != (r + 1, r + w - 1 + p["iota"]):
        bad["inertia"] += 1
    # Theorem 5.2
    if np.any((lam > TOL) & (lam < c - TOL)):
        bad["gap"] += 1
    # Corollary 5.3
    if abs(E - 2 * (lam[:r + 1].sum() - (r + 1) * c)) > 1e-8:
        bad["Epos"] += 1
    # Corollary 4.3
    H, L = (n + 1) // 2, -(-(n + 1) // 2)
    R0 = max(abs(lam[H - 1]), abs(lam[L - 1])) < 1e-8
    if R0 != (max(r + 1, r + w - 1 + p["iota"]) <= (n - 1) // 2):
        bad["HL"] += 1
    # Theorem 6.1
    lo, km = lower_bound(n0, w, ms), km_bound(n0, w, ms)
    if n >= 4 and p["rho_star"] ** 2 < p["Theta"] / (r + 1):
        bad["km_hyp"] += 1
    if not (lo - 1e-9 <= E <= km + 1e-9):
        bad["bounds"] += 1
    if abs(E - km) < 1e-9:
        km_eq.append((n0, w, ms))
    mg = mcclelland_gutman(n0, w, ms)
    if km >= mg:
        km_vs_mcc.append(((n0, w, ms), km, mg, E))
    an = mcclelland_analogue(n0, w, ms)
    if an > mg:
        analog_gt += 1
    if w >= 1:
        ratios_km.append(E / km)
        ratios_an.append(E / an)

print("Theorems and corollaries over all 1157 members")
for key, label in [("inertia", "Theorem 4.2: n+ = r+1, n- = r+w-1+iota"),
                   ("HL", "Corollary 4.3: HL-index criterion"),
                   ("gap", "Theorem 5.2: no eigenvalue in (0, sigma/n)"),
                   ("Epos", "Corollary 5.3: E = 2(s+ - (r+1)c)"),
                   ("km_hyp", "Theorem 6.1: rho_*^2 >= Theta/(r+1) for n >= 4"),
                   ("bounds", "Theorem 6.1: (6.2) <= E <= (6.3)")]:
    check(f"{label}, violations", bad[key], 0)

check("equality in (6.3) exactly for Gamma(0,w;1,...,1)",
      all(n0 == 0 and set(ms) == {1} for (n0, w, ms) in km_eq), True)
check("equality cases of (6.3) = all members Gamma(0,w;1^r) with 3<=n<=13",
      len(km_eq), sum(1 for (n0, w, ms) in fam13 if n0 == 0 and set(ms) == {1}),
      source="expected")

print("Remark 6.2")
check("members where (6.3) is not below the McClelland-type bound of Gutman et al.",
      [m[0] for m in km_vs_mcc], [(0, 1, (2,))])
(_, km0, mg0, E0), = km_vs_mcc
check("Gamma(0,1;2): (6.3), McClelland, E", (round(km0, 3), round(mg0, 3), round(E0, 3)),
      (6.208, 5.916, 5.230))
check("members where the McClelland-type analogue exceeds that bound", analog_gt, 310)

print("Example 6.3 (members with w >= 1)")
check("number of members with 3<=n<=13 and w>=1", len(ratios_km), 799)
check("E/(6.3): mean, min", (round(float(np.mean(ratios_km)), 3), round(float(np.min(ratios_km)), 3)),
      (0.921, 0.819))
check("E/(McClelland-type analogue): mean", round(float(np.mean(ratios_an)), 3), 0.686)

# ---------------------------------------------------------------- n <= 15
print("Section 5.1: dense eigensolver versus secular equation")
fam15 = list(members(3, 15))
check("number of members with 3 <= n <= 15", len(fam15), 2318)
worst = 0.0
for (n0, w, ms) in fam15:
    A, _ = adjacency(n0, w, ms)
    worst = max(worst, abs(energy_dense(A) - energy_secular(n0, w, ms)))
check("largest |E_dense - E_secular|", f"{worst:.1e}", "below 3e-14", ok=worst < 3e-14)

# ---------------------------------------------------------------- Remark 10.1
print("Remark 10.1: spectral determination, all members (r >= 0) with n <= 9")
seen, pairs = {}, []
for (n0, w, ms) in members(1, 9, r_min=0):
    A, _ = adjacency(n0, w, ms)
    key = tuple(np.round(np.sort(np.linalg.eigvalsh(A)), 6))
    if key in seen:
        pairs.append((seen[key], (n0, w, ms)))
    seen[key] = (n0, w, ms)
check("cospectral pairs", pairs, [])

finish()
