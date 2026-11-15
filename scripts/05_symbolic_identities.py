"""Symbolic identities: Lemma 2.6, Remark 2.7, Theorem 3.2, Lemma 3.5, Lemma 5.1,
Lemma 8.1 and (8.2), together with the check of [Hampiholi et al., Thm 2.11] in Remark 3.4."""
from common import check, finish

import numpy as np
import sympy as sp

from starforest import adjacency, members

x, c, a, b, m = sp.symbols("x c a b m", real=True)
mm = sp.symbols("m", positive=True)


def h_(k, t):
    """Star correction h_m of Lemma 3.5 (h_0 = 0)."""
    return 0 if k == 0 else k * (k + 1 - 2 * t) / (t * (t ** 2 - k))


print("Lemma 2.6 and Remark 2.7 (stars K_{1,k}, k = 1..6)")
ok_det = ok_cor = ok_pf = True
for k in range(1, 7):
    B = sp.zeros(k + 1, k + 1)
    for j in range(1, k + 1):
        B[0, j] = B[j, 0] = 1
    ok_det &= sp.expand((x * sp.eye(k + 1) + B).det() - x ** (k - 1) * (x ** 2 - k)) == 0
    one = sp.ones(k + 1, 1)
    cor = sp.simplify((one.T * (x * sp.eye(k + 1) + B).inv() * one)[0])
    ok_cor &= sp.simplify(cor - ((k + 1) * x - 2 * k) / (x ** 2 - k)) == 0
    pf = (sp.sqrt(k) - 1) ** 2 / 2 / (x - sp.sqrt(k)) + (sp.sqrt(k) + 1) ** 2 / 2 / (x + sp.sqrt(k))
    ok_pf &= sp.simplify(pf - ((k + 1) * x - 2 * k) / (x ** 2 - k)) == 0
check("det(xI + A(K_{1,m})) = x^(m-1)(x^2 - m)", ok_det, True)
check("coronal_{-A(K_{1,m})} = ((m+1)x - 2m)/(x^2 - m)", ok_cor, True)
check("partial fractions of Remark 2.7", ok_pf, True)

print("Lemma 3.5 and Lemma 5.1")
check("((m+1)x - 2m)/(x^2-m) - (m+1)/x = h_m(x)",
      sp.simplify(((mm + 1) * x - 2 * mm) / (x ** 2 - mm) - (mm + 1) / x - h_(mm, x)), 0)
check("h_m(c) = g_m(c) - (m+1)/c", sp.simplify(h_(mm, c) + (mm + 1) / c
                                               - (2 * mm - (mm + 1) * c) / (mm - c ** 2)), 0)
check("g_1 in reduced form is 2/(1+c)", sp.simplify((2 - 2 * c) / (1 - c ** 2) - 2 / (1 + c)), 0)
check("g_m(1) = 1 for m >= 2", sp.simplify((2 * mm - (mm + 1)) / (mm - 1)), 1)

print("Lemma 8.1 (transfer) and (8.2) (merge)")
lhs = h_(a + 1, x) + h_(b - 1, x) - h_(a, x) - h_(b, x)
rhs = (a - b + 1) * x * (x - 1) ** 2 * (2 * x ** 2 - a - b) / (
    (x ** 2 - a) * (x ** 2 - b) * (x ** 2 - a - 1) * (x ** 2 - b + 1))
check("transfer identity (b >= 2)", sp.simplify(lhs - rhs), 0)
ok_b1 = all(sp.simplify(h_(A_ + 1, x) - h_(A_, x) - h_(1, x)
                        - rhs.subs({a: A_, b: 1})) == 0 for A_ in range(1, 7))
check("transfer identity (b = 1, h_0 = 0), a = 1..6", ok_b1, True)
merge = ((a * b * (x - 1) ** 2 * (2 * x ** 2 - a - b))
         / (x * (x ** 2 - a) * (x ** 2 - b) * (x ** 2 - a - b)))
check("merge identity", sp.simplify(h_(a + b, x) - h_(a, x) - h_(b, x) - merge), 0)

print("Theorem 3.2: chi = D(1 - Phi), all members with n <= 9")


def chi_formula(n0, w, ms):
    D = x ** n0 * (x + 1) ** w
    Phi = sp.Integer(n0) / x + sp.Integer(w) / (x + 1)
    for k in ms:
        D *= x ** (k - 1) * (x ** 2 - k)
        Phi += ((k + 1) * x - 2 * k) / (x ** 2 - k)
    return sp.Poly(sp.cancel(D * (1 - Phi)), x)


bad = 0
for mem in members(1, 9, r_min=0):
    A, _ = adjacency(*mem)
    true = sp.Matrix(A.astype(int)).charpoly(x)
    bad += (true - chi_formula(*mem)).is_zero is False
check("members violating Theorem 3.2", bad, 0)

print("Remark 3.4: [Hampiholi et al., Thm 2.11], K_p minus k stars K_{1,m-1}")


def cl_charpoly(p, mstar, k):
    """Integer coefficients of the characteristic polynomial of K_p minus k stars of order mstar."""
    A, _ = adjacency(p - mstar * k, 0, (mstar - 1,) * k, h=0.0)
    return np.rint(np.poly(np.linalg.eigvalsh(A))).astype(np.int64)


def coeffs(expr):
    return np.array(sp.Poly(sp.expand(expr), x).all_coeffs(), dtype=np.int64)


wrong, fixed_wrong, cases, failing = 0, 0, 0, set()
for p in range(3, 12):
    for ms_ in range(2, p + 1):
        for k in range(1, p // ms_ + 1):
            if p - 2 * k - 1 < 0:
                continue
            cases += 1
            base = (x + 1) ** (p - 2 * k - 1) * ((x + 1) ** 2 - (ms_ - 1)) ** (k - 1)
            cub = x ** 3 - (p - 3) * x ** 2 + (2 * ms_ * k - 2 * p - 2 * k - ms_ + 4) * x
            true = cl_charpoly(p, ms_, k)
            if not np.array_equal(true, coeffs(base * (cub + (p - ms_ * k) * (ms_ - 2)))):
                wrong += 1
                failing.add((k >= 2, ms_ >= 3))
            if not np.array_equal(true, coeffs(base * (cub + (p - ms_ * k + k - 1) * (ms_ - 2)))):
                fixed_wrong += 1
print(f"  (cases checked, p <= 11: {cases}; printed formula wrong in {wrong})")
check("...and wrong exactly when k >= 2 and m >= 3", failing, {(True, True)})
check("corrected constant (p - mk + k - 1)(m - 2) wrong in", fixed_wrong, 0)
check("chi of K_9 - 2K_{1,2}", sp.factor(sp.Poly(cl_charpoly(9, 3, 2).tolist(), x).as_expr()),
      sp.factor((x + 1) ** 4 * (x ** 2 + 2 * x - 1) * (x ** 3 - 6 * x ** 2 - 9 * x + 4)),
      ok=sp.expand(sp.Poly(cl_charpoly(9, 3, 2).tolist(), x).as_expr()
                   - (x + 1) ** 4 * (x ** 2 + 2 * x - 1) * (x ** 3 - 6 * x ** 2 - 9 * x + 4)) == 0)

finish()
