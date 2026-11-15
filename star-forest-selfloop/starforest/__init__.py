"""Reproducibility code for "Inertia and energy of star-forest self-loop graphs".

Modules
-------
family   : the graphs Gamma^h(n0, w; m_1, ..., m_r), adjacency matrices, enumeration.
secular  : the secular function Phi_h, its poles, bracketed root finding, and the
           spectrum and energy computed without forming the adjacency matrix.
"""
from .family import order, adjacency, members, partitions
from .secular import (multiplicities, secular, poles, roots, structural_spectrum,
                      spectrum, energy_secular, energy_dense)

__all__ = ["order", "adjacency", "members", "partitions", "multiplicities", "secular",
           "poles", "roots", "structural_spectrum", "spectrum", "energy_secular",
           "energy_dense"]
