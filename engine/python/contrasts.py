from __future__ import annotations
from itertools import permutations
from .config import Contrast


def expand_contrasts(explicit, all_vs_all, conditions):
    out = list(explicit)
    if all_vs_all:
        have = {(c.numerator, c.denominator) for c in out}
        for a, b in permutations(sorted(set(conditions)), 2):
            if (a, b) not in have:
                out.append(Contrast(name=f"{a}_vs_{b}", numerator=a, denominator=b))
                have.add((a, b))
    return out
