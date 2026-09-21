#!/usr/bin/env python3
"""
metrics.py — canonical, dependency-light quality-indicator functions
shared across this package's experiment scripts. Round 8: extracted
from four separately-duplicated copies (run_baselines_comparison_v3.py,
run_ablation_v3.py, binary_baseline_v3.py,
build_binary_vs_decoder_shared_ref_v3.py) into one canonical
definition, per an external review's correct observation that a
pytest test for spacing had to import a full pymoo-dependent
experiment script just to reach this one function. Only depends on
numpy, so `tests/test_reproducibility.py`'s spacing test (and anyone
else) can import it without pymoo installed.
"""
import numpy as np


def compute_spacing(F):
    """Canonical Schott (1995) spacing: Manhattan (L1) nearest-neighbour
    distance, sample std (ddof=1). FIX (external review round 6-followup,
    confirmed real): this function previously used Euclidean (L2) distance
    and ddof=0 (population std) -- a different, undisclosed-as-different
    variant from the metric this paper cites by name and reference
    (Schott 1995)."""
    if len(F) < 2:
        return 0.0
    n = len(F)
    min_dists = np.full(n, np.inf)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = np.abs(F[i] - F[j]).sum()
            if d < min_dists[i]:
                min_dists[i] = d
    min_dists = min_dists[np.isfinite(min_dists)]
    return float(np.std(min_dists, ddof=1)) if len(min_dists) > 1 else 0.0
