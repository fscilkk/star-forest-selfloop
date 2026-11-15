"""Section 7: loops of weight h.  Theorem 7.1, its sharpness example, the remarks after it,
Remark 7.2, and the claim of Section 5.1 that the secular procedure applies to Gamma^h."""
from common import check, finish

import numpy as np

from starforest import adjacency, members, spectrum, energy_secular, energy_dense

TOL = 1e-9


def statements(n0, w, ms, h):
    """Truth values of Theorem 7.1 (i), (ii), (iii) for Gamma^h(n0, w; ms)."""
    A, _ = adjacency(n0, w, ms, h)
    n, r = A.shape[0], len(ms)
    lam = np.sort(np.linalg.eigvalsh(A))[::-1]
    ch = h * (n - w) / n
    i = int((lam > TOL).sum()) == r + 1
    ii = not np.any((lam > TOL) & (lam < ch - TOL))
    E = float(np.abs(lam - ch).sum())
    iii = abs(E - 2 * (lam[:r + 1].sum() - (r + 1) * ch)) < 1e-8
    return i, ii, iii


fam11 = list(members(3, 11))
check("number of members with 3 <= n <= 11", len(fam11), 539)

print("Theorem 7.1 for 0 < h <= 1 (h = 0.1, 0.2, ..., 1.0), all 539 members")
viol = 0
for h in np.round(np.arange(0.1, 1.01, 0.1), 10):
    for m in fam11:
        viol += not all(statements(*m, h))
check("violations of (i), (ii) or (iii)", viol, 0)

print("Sharpness: Gamma^h(0,0;2) has spectrum {h+1, h, h-1}")
for h in (0.5, 1.0, 1.5, 3.0):
    s = spectrum(0, 0, (2,), h)
    check(f"h = {h}", [round(v, 9) for v in s], [h + 1, h, h - 1])
i, ii, _ = statements(0, 0, (2,), 1.5)
check("(i) and (ii) fail for Gamma^h(0,0;2) at h = 1.5", (i, ii), (False, False))
check("Gamma^h(0,0;1) = h I_2 satisfies (i), (ii) at h = 5",
      statements(0, 0, (1,), 5.0)[:2], (True, True))

print("After Theorem 7.1")
fam10 = list(members(2, 10))
check("members with n <= 10 violating (iii) at h = 1.001",
      sum(not statements(*m, 1.001)[2] for m in fam10), 0)
check("Gamma^h(0,4;2) violates (iii) at h = 2", not statements(0, 4, (2,), 2.0)[2], True)

print("Remark 7.2 (539 members with 3 <= n <= 11)")
check("(i) still holds at h = 1.001", sum(statements(*m, 1.001)[0] for m in fam11), 54)
check("(i) still holds at h = 1.5", sum(statements(*m, 1.5)[0] for m in fam11), 29)

print("Section 5.1: the secular procedure for Gamma^h (with mu_0 = n - 1 + max(h, 1))")
worst = 0.0
for h in (0.3, 0.7, 1.3, 2.0, 4.0):
    for m in members(2, 11):
        A, _ = adjacency(*m, h)
        worst = max(worst, abs(energy_dense(A) - energy_secular(*m, h)))
check("largest |E_dense - E_secular| over n <= 11 and h in {0.3,0.7,1.3,2,4}",
      f"{worst:.1e}", "< 1e-12", ok=worst < 1e-12)

finish()
