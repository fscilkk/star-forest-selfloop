"""Closed-form quantities of Section 6 (all at h = 1)."""
import math
from collections import Counter

from .family import order


def invariants(n0, w, ms):
    """Parameters used in Corollary 5.5 and Theorem 6.1."""
    nu = Counter(ms)
    n = order(n0, w, ms)
    r = len(ms)
    M = sum(ms)
    sigma = n - w
    c = sigma / n
    nu1 = nu.get(1, 0)
    eps = 1 if (w >= 1 or nu1 >= 1) else 0
    iota = 1 if n0 >= 1 else 0
    d = sum(1 for m in nu if m >= 2)
    T = n * n - 2 * M - w                                                # eq. (6.1)
    Theta = T - (w + nu1 - eps) - sum((v - 1) * m for m, v in nu.items() if m >= 2)
    Sigma_plus = sum((v - 1) * math.sqrt(m) for m, v in nu.items() if m >= 2) + nu1
    return dict(n=n, r=r, M=M, sigma=sigma, c=c, nu=nu, eps=eps, iota=iota, d=d,
                T=T, Theta=Theta, rho_star=T / n, Sigma_plus=Sigma_plus)


def lower_bound(n0, w, ms):
    """(6.2):  E >= 2(rho_* + Sigma^+ - (r+1)c)."""
    p = invariants(n0, w, ms)
    return 2 * (p["rho_star"] + p["Sigma_plus"] - (p["r"] + 1) * p["c"])


def km_bound(n0, w, ms):
    """(6.3):  Koolen-Moulton-type upper bound."""
    p = invariants(n0, w, ms)
    rs, r = p["rho_star"], p["r"]
    return 2 * rs + 2 * math.sqrt(r * (p["Theta"] - rs ** 2)) - 2 * (r + 1) * p["c"]


def mcclelland_analogue(n0, w, ms):
    """The McClelland-type analogue displayed after Theorem 6.1."""
    p = invariants(n0, w, ms)
    return 2 * math.sqrt((p["r"] + 1) * p["Theta"]) - 2 * (p["r"] + 1) * p["c"]


def mcclelland_gutman(n0, w, ms):
    """McClelland-type bound of Gutman et al. [Thm 6]:  sqrt(n (T - sigma^2/n))."""
    p = invariants(n0, w, ms)
    return math.sqrt(p["n"] * (p["T"] - p["sigma"] ** 2 / p["n"]))
