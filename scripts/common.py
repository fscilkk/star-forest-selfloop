"""Shared helpers: import path and a uniform way of reporting each claim of the paper."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

_failures = []


def _plain(v):
    """Turn numpy scalars (also inside tuples/lists) into plain Python numbers for printing."""
    if isinstance(v, (list, tuple)):
        return type(v)(_plain(u) for u in v)
    if hasattr(v, "item") and not hasattr(v, "__len__"):
        return v.item()
    return v


def check(label, computed, claimed=None, ok=None, source="paper"):
    """Print 'label: computed (paper: claimed) [OK|FAIL]' and record failures.

    `source` names where the claimed value comes from: "paper" for a value printed in the
    paper, "expected" for a value derived from a statement of the paper rather than printed."""
    computed = _plain(computed)
    if ok is None:
        ok = computed == claimed
    tag = "OK" if ok else "FAIL"
    if claimed is None:
        print(f"  [{tag}] {label}: {computed}")
    else:
        print(f"  [{tag}] {label}: {computed}   ({source}: {claimed})")
    if not ok:
        _failures.append(label)


def close(a, b, digits):
    """True if a and b agree when rounded to the number of decimals printed in the paper."""
    return round(a, digits) == round(b, digits)


def finish():
    if _failures:
        print(f"\n{len(_failures)} check(s) FAILED:", *_failures, sep="\n  ")
        sys.exit(1)
    print("\nall checks passed")
