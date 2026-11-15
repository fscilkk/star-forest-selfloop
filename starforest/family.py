"""The family Gamma^h(n0, w; m_1, ..., m_r) of Definition 2.1 (and Section 7 for h != 1).

Gamma is K_n with the edges of r disjoint stars K_{1,m_l} deleted, a loop of weight h at
every vertex except the w vertices of W, and n0 further looped vertices N that meet no star.
Order: n = n0 + w + sum_l (m_l + 1).
"""
import numpy as np


def order(n0, w, ms):
    return n0 + w + sum(m + 1 for m in ms)


def adjacency(n0, w, ms, h=1.0):
    """Adjacency matrix of Gamma^h(n0, w; ms) and its vertex classes.

    Vertices are ordered star by star (centre first, then its leaves), then W, then N.
    Returns (A, classes) with classes = dict(centers, leaves, W, N).
    """
    n = order(n0, w, ms)
    A = np.ones((n, n)) - np.eye(n)
    centers, leaves = [], []
    i = 0
    for m in ms:
        c = i
        L = list(range(i + 1, i + 1 + m))
        A[c, L] = 0.0
        A[L, c] = 0.0
        centers.append(c)
        leaves.append(L)
        i += m + 1
    W = list(range(i, i + w))
    N = list(range(i + w, i + w + n0))
    diag = np.full(n, float(h))
    diag[W] = 0.0
    A[np.arange(n), np.arange(n)] = diag
    return A, dict(centers=centers, leaves=leaves, W=W, N=N)


def partitions(total, min_part=1, max_part=None):
    """Partitions of `total` into parts >= min_part, as nonincreasing tuples."""
    if max_part is None:
        max_part = total
    if total == 0:
        yield ()
        return
    for p in range(min(total, max_part), min_part - 1, -1):
        for rest in partitions(total - p, min_part, p):
            yield (p,) + rest


def members(nmin, nmax, r_min=1):
    """All (n0, w, ms) with nmin <= n <= nmax and at least r_min stars.

    ms is the nonincreasing tuple of star sizes m_l >= 1 (stars K_{1,m_l} have m_l+1
    vertices).  Each member of the family is listed exactly once.
    """
    for n in range(nmin, nmax + 1):
        for w in range(0, n + 1):
            for k in range(0, n - w + 1):          # k = number of vertices in the stars
                if k == 0:
                    if r_min == 0:
                        yield (n - w, w, ())
                    continue
                for orders in partitions(k, min_part=2):
                    ms = tuple(o - 1 for o in orders)
                    if len(ms) >= r_min:
                        yield (n - w - k, w, ms)
