"""Secular-equation computation of the spectrum and energy (Sections 3-5, 7).

With y = x + 1 - h the secular function of Gamma^h(n0, w; m_1..m_r) is (eq. (7.2))

    Phi_h(x) = n0/y + w/(x+1) + sum_l ((m_l+1) y - 2 m_l) / (y^2 - m_l),

a star of size 1 being entered in reduced form 2/(y+1).  The eigenvalues are the roots of
Phi_h = 1 (one in each gap between consecutive poles, plus one above the largest pole,
Theorem 4.1) together with the structural eigenvalues, the zeros of det(xI + Z_h) that are
not cancelled by a pole (Corollary 3.6).  Nothing of size n is ever formed, so n may be huge.
"""
import math
from collections import Counter

try:
    import mpmath
except ImportError:  # only needed for the high-precision check
    mpmath = None

_TOL = 1e-12


def multiplicities(ms):
    """nu_m = number of stars of size m.  `ms` may be a sequence of star sizes or a mapping
    {m: nu_m}; the second form keeps the cost independent of the number of stars."""
    return Counter(dict(ms)) if hasattr(ms, "items") else Counter(ms)


def _order(n0, w, nu):
    return n0 + w + sum((m + 1) * v for m, v in nu.items())


def _num(mp):
    if mp:
        return mpmath.mpf, mpmath.sqrt
    return float, math.sqrt


def secular(x, n0, w, nu, h=1.0, mp=False):
    """Phi_h(x), grouping equal stars (Section 5.1)."""
    F, _ = _num(mp)
    y = x + 1 - F(h)
    val = F(0)
    if n0:
        val += n0 / y
    if w:
        val += w / (x + 1)
    for m, v in nu.items():
        if m == 1:
            val += v * 2 / (y + 1)
        else:
            val += v * ((m + 1) * y - 2 * m) / (y * y - m)
    return val


def _merge(values):
    """Merge numerically equal reals (e.g. x = -1 from W and from stars of size 1 at h = 1)."""
    out = []
    for v in sorted(values, reverse=True):
        if not out or abs(out[-1] - v) > _TOL:
            out.append(v)
    return out


def poles(n0, w, nu, h=1.0):
    """Distinct poles of Phi_h in decreasing order (eq. (3.3) for h = 1)."""
    P = []
    if n0:
        P.append(h - 1)
    if w:
        P.append(-1.0)
    for m in nu:
        if m == 1:
            P.append(h - 2)                    # y = -1; y = +1 is not a pole (Remark 2.7)
        else:
            P += [h - 1 + math.sqrt(m), h - 1 - math.sqrt(m)]
    return _merge(P)


def roots(n0, w, ms, h=1.0, mp=False, dps=40, maxit=400, above=None):
    """All roots theta_1 > ... > theta_t of Phi_h = 1, by bisection in the brackets of
    Theorem 4.1, the top bracket being closed by mu_0 = n - 1 + max(h, 1)
    (Proposition 5.7 and the paragraph after it).  With `above` set, only the brackets
    reaching above that value are searched (the energy needs only theta_1..theta_{d+1})."""
    nu = multiplicities(ms)
    n = _order(n0, w, nu)
    P = poles(n0, w, nu, h)
    if mp:
        mpmath.mp.dps = dps
    F, _ = _num(mp)
    brackets = [n - 1 + max(h, 1.0)] + P
    out = []
    for j in range(len(P)):
        if above is not None and brackets[j] <= above:
            break
        lo, hi = F(brackets[j + 1]), F(brackets[j])
        for _ in range(maxit):
            mid = (lo + hi) / 2
            if not (lo < mid < hi):
                break
            if secular(mid, n0, w, nu, h, mp) > 1:
                lo = mid
            else:
                hi = mid
        out.append((lo + hi) / 2)
    return out


def structural_spectrum(n0, w, ms, h=1.0):
    """Eigenvalues forced by the block structure (Corollary 3.6), as {value: multiplicity}.

    Zeros of det(xI + Z_h):  h-1 (mult. n0 + sum(m-1)),  -1 (mult. w),
    h-1 +- sqrt(m) (mult. nu_m) for m >= 2,  h and h-2 (mult. nu_1);
    one copy of each pole of Phi_h is removed.
    """
    nu = multiplicities(ms)
    zeros = Counter()

    def add(v, k):
        if k <= 0:
            return
        for key in zeros:
            if abs(key - v) <= _TOL:
                zeros[key] += k
                return
        zeros[v] += k

    add(h - 1, n0 + sum((m - 1) * v for m, v in nu.items()))
    add(-1.0, w)
    for m, v in nu.items():
        if m == 1:
            add(h, v)
            add(h - 2, v)
        else:
            add(h - 1 + math.sqrt(m), v)
            add(h - 1 - math.sqrt(m), v)
    for p in poles(n0, w, nu, h):
        for key in zeros:
            if abs(key - p) <= _TOL:
                zeros[key] -= 1
                break
    return {k: v for k, v in zeros.items() if v > 0}


def spectrum(n0, w, ms, h=1.0):
    """Full spectrum, decreasing (only sensible for moderate n)."""
    ev = list(roots(n0, w, ms, h))
    for v, k in structural_spectrum(n0, w, ms, h).items():
        ev += [v] * k
    return sorted(ev, reverse=True)


def energy_secular(n0, w, ms, h=1.0, mp=False, dps=40):
    """E(Gamma^h) = sum_i |lambda_i - c| with c = h*sigma/n, computed as
    2 * sum_{lambda > c} (lambda - c), which holds because the deviations sum to zero.
    For 0 < h <= 1 this is the formula of Corollaries 5.3 and 5.5 / Theorem 7.1(iii)."""
    nu = multiplicities(ms)
    n = _order(n0, w, nu)
    F, _ = _num(mp)
    if mp:
        mpmath.mp.dps = dps
    c = F(h) * (n - w) / n
    total = F(0)
    for t in roots(n0, w, nu, h, mp=mp, dps=dps, above=float(c)):
        if t > c:
            total += t - c
    for v, k in structural_spectrum(n0, w, nu, h).items():
        if v > c:
            total += k * (F(v) - c)
    return 2 * total


def energy_dense(A):
    """Energy from a dense eigensolver: sum |lambda_i - tr(A)/n|."""
    import numpy as np
    lam = np.linalg.eigvalsh(A)
    return float(np.abs(lam - np.trace(A) / A.shape[0]).sum())
