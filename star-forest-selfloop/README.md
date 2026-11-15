# Inertia and energy of star-forest self-loop graphs — reproducibility code

Code accompanying

> L. Kittiratanawasin and K. Rodtes, *Inertia and energy of star-forest self-loop graphs* (submitted).

It contains

* `starforest/` — a small package implementing the family Γ^h(n₀, w; m₁, …, m_r), its
  secular function Φ_h, and the bracketed root-finding algorithm of Proposition 5.7, which
  computes the spectrum and the energy **without forming the adjacency matrix**;
* `scripts/` — one script per group of claims in the paper. Each script prints every number
  it checks next to the value stated in the paper and marks it `[OK]` or `[FAIL]`;
* `results/` — the outputs of the last run; `table1_paper.tex` holds the rows of Table 1 exactly
  as printed in the paper;
* `figure/` — the TikZ source and PDF of Figure 1.

## Quick start

```bash
pip install -r requirements.txt
python run_all.py          # scripts 01–05 and a quick Table 1 (about 1 minute;
                           #   writes results/table1_quick.tex, leaving table1.tex untouched)
python run_all.py --full   # full Table 1 as well (a few more minutes, depending on the machine)
```

Every script can also be run on its own from `scripts/`, e.g. `python 02_exhaustive_checks.py`.
Tested with Python 3.8 (Windows 10) and Python 3.12 (Linux).

### Using the algorithm

```python
from starforest import energy_secular, spectrum

# Gamma(n0, w; m_1, ..., m_r) with star sizes given as a list ...
energy_secular(1, 2, (2, 2))                      # 16.7256...   (after Corollary 5.5)
spectrum(1, 2, (2, 2))                            # all 9 eigenvalues

# ... or, for very many stars, as {size: number of stars}; the cost then depends on n
# only through a logarithmic factor and the number of distinct star sizes (Section 5.1)
energy_secular(56_666_671, 10_000_000, {1: 1_666_667, 2: 1_666_667, 3: 1_666_667,
                                        4: 1_666_666, 5: 1_666_666})   # n = 10^8, E = 212941102.407
```

All functions accept a loop weight `h` (Section 7), e.g. `energy_secular(1, 2, (2, 2), h=0.5)`.

## Reproducing the paper

1. Install Python 3.8 or later and the packages: `pip install -r requirements.txt`.
2. From this folder run `python run_all.py --full`.
3. The last line should read `All scripts passed.` Every value the scripts check is then equal
   to the value printed in the paper (see the table below for where each one appears).

What to expect:

* Scripts 01–05 and the energies of Table 1 give exactly the values in the paper on any
  machine.
* The **timings** of Table 1, and the absolute dense-versus-secular discrepancies, depend on
  the machine and change from run to run. The table printed in the paper is kept unchanged in
  `results/table1_paper.tex`, with the full output of that run in
  `results/06_timing_table1_paper.txt` (Intel Core i7-4710HQ, 16 GB, Windows 10, Python 3.8.5,
  NumPy 1.24.4). A new run writes `results/table1.tex` and `results/06_timing_table1.txt` for
  your machine and never touches the `_paper` files.
* The checks that do not depend on the machine are the energies, the relative agreement of
  the two methods (below 10⁻¹⁴) and the 40-digit comparison (below 10⁻¹⁵).

## Where each claim of the paper is checked

| Paper | Claim | Script |
|---|---|---|
| Example 3.7 | χ_Γ of Γ(1,2;2,2) and Γ(0,1;3,2); roots of the quartic | `01_examples.py` |
| after Cor. 5.5 | E(Γ(1,2;2,2)) = 16.726 | `01_examples.py` |
| Example 6.3 | T = 71, Θ = 68; 13.94 ≤ 16.73 ≤ 17.90; 23.90 | `01_examples.py` |
| Example 6.3 | 799 members with w ≥ 1: mean 0.921, min 0.819; 0.686 | `02_exhaustive_checks.py` |
| Example 9.2, Remark 9.3 | positive roots; the values of 2P_vv | `01_examples.py` |
| Remark 4.4 | spanning star: 1 versus 2 positive eigenvalues | `01_examples.py` |
| Thm 4.2, Cor 4.3, Thm 5.2, Cor 5.3 | inertia, HL-index, forbidden interval, energy as a positive part, all 1157 members with 3 ≤ n ≤ 13 | `02_exhaustive_checks.py` |
| Theorem 6.1 | bounds, the hypothesis for n ≥ 4, the equality case | `02_exhaustive_checks.py` |
| Remark 6.2 | (6.3) below the McClelland-type bound except Γ(0,1;2) (6.208, 5.916, 5.230); 310 | `02_exhaustive_checks.py` |
| Section 5.1 | 2318 members with 3 ≤ n ≤ 15: dense and secular energies agree to below 3·10⁻¹⁴ | `02_exhaustive_checks.py` |
| Remark 10.1 | no cospectral pair with n ≤ 9 | `02_exhaustive_checks.py` |
| Theorem 7.1 | (i)–(iii) for 0 < h ≤ 1; sharpness Γ^h(0,0;2); Γ^h(0,0;1) | `03_weighted_loops.py` |
| after Thm 7.1 | (iii) at h = 1.001 for n ≤ 10; Γ^h(0,4;2) at h = 2 | `03_weighted_loops.py` |
| Remark 7.2 | (i) holds for 54 members at h = 1.001 and 29 at h = 1.5 (539 members) | `03_weighted_loops.py` |
| Section 5.1 | the secular procedure for Γ^h with μ₀ = n − 1 + max(h, 1) | `03_weighted_loops.py` |
| Theorem 8.2, Cor. 8.3 | dominance monotonicity for several h; ρ(Γ_(1^M)) | `04_dominance_and_open_questions.py` |
| Remark 8.4 | walk counts mwc_k(F), strict increase under transfers | `04_dominance_and_open_questions.py` |
| Remark 10.2 | (i) incomparable pairs at M = 6, 8.8446 vs 8.8437, no equal radii; (ii) energy not monotone; (iii) coronal condition, P₄ fails | `04_dominance_and_open_questions.py` |
| Lemmas 2.6, 3.5, 5.1, 8.1, (8.2), Thm 3.2 | symbolic identities; χ = D(1 − Φ) for all members with n ≤ 9 | `05_symbolic_identities.py` |
| Remark 3.4 | [Hampiholi et al., Thm 2.11] is wrong exactly when k ≥ 2 and m ≥ 3; corrected constant term | `05_symbolic_identities.py` |
| Table 1 | energies, dense versus secular timings, 40-digit check | `06_timing_table1.py` |

## Notes

* **Reading the output.** Each line `[OK] label: computed (paper: claimed)` compares a
  computed value with the value printed in the paper. A line marked `(expected: …)` compares
  with a value derived from a statement of the paper rather than printed in it.
* **Timings** (Table 1) depend on the machine and on the LAPACK build; see *Reproducing the
  paper* above. The absolute dense-versus-secular discrepancies also depend on the machine,
  since they are dominated by the rounding error of the dense eigensolver.
* The dense computations are single-threaded (the timing script sets `OMP_NUM_THREADS=1`
  and the analogous variables).
* High-precision evaluation uses `mpmath` at 40 significant digits.

## Licence

MIT, see `LICENSE`. If you use the code, please cite the paper (see `CITATION.cff`).
