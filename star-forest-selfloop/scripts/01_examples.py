"""Worked examples of the paper: Example 3.7, the energy after Corollary 5.5,
Example 6.3, Example 9.2, Remark 9.3, Remark 4.4 and two small spectra quoted in the text."""
from common import check, close, finish

import numpy as np
import sympy as sp

from starforest import adjacency, spectrum, roots, energy_secular
from starforest.bounds import invariants, lower_bound, km_bound, mcclelland_analogue

x = sp.symbols("x")


def charpoly(n0, w, ms):
    A, _ = adjacency(n0, w, ms)
    return sp.factor(sp.Matrix(A.astype(int)).charpoly(x).as_expr())


print("Example 3.7")
check("chi of Gamma(1,2;2,2)", charpoly(1, 2, (2, 2)),
      sp.factor(x**2 * (x + 1) * (x**2 - 2) * (x**4 - 8*x**3 - x**2 + 12*x + 2)),
      ok=sp.expand(charpoly(1, 2, (2, 2)) - x**2*(x+1)*(x**2-2)*(x**4-8*x**3-x**2+12*x+2)) == 0)
q = np.roots([1, -8, -1, 12, 2])
q = sorted(q.real, reverse=True)
check("roots of the quartic", [round(v, 6) for v in q], [7.931312, 1.350650, -0.167528, -1.114435])
check("chi of Gamma(0,1;3,2)", sp.expand(charpoly(0, 1, (3, 2))),
      sp.expand(x**3 * (x**5 - 7*x**4 - 2*x**3 + 27*x**2 - x - 24)))

print("Energy after Corollary 5.5")
E = energy_secular(1, 2, (2, 2))
check("E(Gamma(1,2;2,2))", round(E, 3), 16.726)

print("Example 6.3")
p = invariants(1, 2, (2, 2))
check("T, Theta", (p["T"], p["Theta"]), (71, 68))
check("lower <= E <= KM", (round(lower_bound(1, 2, (2, 2)), 2), round(E, 2),
                           round(km_bound(1, 2, (2, 2)), 2)), (13.94, 16.73, 17.90))
check("McClelland-type analogue", round(mcclelland_analogue(1, 2, (2, 2)), 2), 23.90)

print("Example 9.2 (Gamma(1,1;3,2))")
th = [t for t in roots(1, 1, (3, 2)) if t > 0]
check("positive roots of Phi = 1", [round(t, 4) for t in th], [7.8638, 1.6239, 1.3633])


def charges(n0, w, ms):
    A, cl = adjacency(n0, w, ms)
    lam, V = np.linalg.eigh(A)
    Vp = V[:, lam > 1e-9]
    return 2 * np.einsum("ij,ij->i", Vp, Vp), cl


q2, cl = charges(1, 1, (3, 2))
got = [round(q2[cl["centers"][0]], 4), round(q2[cl["leaves"][0][0]], 4),
       round(q2[cl["centers"][1]], 4), round(q2[cl["leaves"][1][0]], 4),
       round(q2[cl["W"][0]], 4), round(q2[cl["N"][0]], 4)]
check("2P_vv (3-centre, 3-leaf, 2-centre, 2-leaf, W, N)", got,
      [1.3385, 0.4745, 1.3068, 0.6688, 0.2473, 0.3461])
check("sum of 2P_vv", round(q2.sum(), 10), 6.0)

print("Remark 9.3 (Gamma(0,1;4))")
q3, cl = charges(0, 1, (4,))
check("2P_vv on W and on a leaf", (round(q3[cl["W"][0]], 3), round(q3[cl["leaves"][0][0]], 3)),
      (0.519, 0.446))

print("Remark 4.4 (spanning star, n = 7)")
A, _ = adjacency(0, 0, (6,), h=0.0)                 # K_7 minus K_{1,6}, no loops
check("positive eigenvalues of K_7 - K_{1,6}", int((np.linalg.eigvalsh(A) > 1e-9).sum()), 1)
A, _ = adjacency(0, 0, (6,))
check("positive eigenvalues of Gamma(0,0;6)", int((np.linalg.eigvalsh(A) > 1e-9).sum()), 2)

print("Spectra quoted in the text")
for m in range(2, 6):
    s = [round(v, 9) for v in spectrum(0, 0, (m,))]
    check(f"Spec Gamma(0,0;{m}) = {{m, 1, 0^(m-1)}}", s, [float(m), 1.0] + [0.0] * (m - 1),
          ok=np.allclose(s, [m, 1] + [0] * (m - 1)))
check("Spec Gamma(0,0;1) = {1, 1}", [round(v, 9) for v in spectrum(0, 0, (1,))], [1.0, 1.0])

finish()
