#!/usr/bin/env python3
"""
coverage_threshold_sensitivity_v3.py — sensitivity of the coverage
objective C_tau to the fidelity threshold (canonical: 0.70) used to
decide whether a city pair counts as "covered". An external review
asked whether this threshold, used nowhere else in the oracle, is
arbitrary enough to threaten the coverage-dependent findings (the
fidelity/rate/coverage correlation structure, the effective-
dimensionality result, and the classical-baseline dominance check).

Re-decodes every point on the pooled evolutionary front (from its
stored 11-variable parameter vector, X_pareto in global_pareto_v3.npz)
and recomputes ONLY coverage at five threshold values, reusing the
canonical oracle's own path-search machinery (not a proxy) -- this
does not require re-running any evolutionary search, since coverage at
a different threshold is a deterministic function of the same
computed per-pair fidelities the oracle already derives.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from bdcz_oracle_v3 import N_NODES, decode_topology, _build_best_medium_adjacency, _werner_optimal_paths_from_source

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT_JSON = RESULTS / 'coverage_threshold_sensitivity_v3.json'
THRESHOLDS = [0.50, 0.60, 0.70, 0.80, 0.90]


def coverage_at_thresholds(eu, ev, et, thresholds):
    N = N_NODES
    total_pairs = N * (N - 1) // 2
    adj = _build_best_medium_adjacency(eu, ev, et)
    counts = {t: 0 for t in thresholds}
    for src in range(N):
        best_F, pred = _werner_optimal_paths_from_source(adj, src)
        for tgt in range(src + 1, N):
            if best_F[tgt] <= 0.0 or pred[tgt] < 0:
                continue
            for t in thresholds:
                if best_F[tgt] >= t:
                    counts[t] += 1
    return {t: counts[t] / total_pairs for t in thresholds}


def main():
    gp = np.load(RESULTS / 'global_pareto_v3.npz')
    X = gp['X_pareto']
    F = gp['F_pareto']
    n = len(X)
    print(f'Recomputing coverage at {len(THRESHOLDS)} thresholds for {n} pooled-front points...')

    cov_by_threshold = {t: [] for t in THRESHOLDS}
    for i, params in enumerate(X):
        result = decode_topology(params)
        if result is None:
            for t in THRESHOLDS:
                cov_by_threshold[t].append(0.0)
            continue
        eu, ev, et = result
        covs = coverage_at_thresholds(eu, ev, et, THRESHOLDS)
        for t in THRESHOLDS:
            cov_by_threshold[t].append(covs[t])
        if (i + 1) % 100 == 0:
            print(f'  {i+1}/{n}')

    cov_baseline = np.array(cov_by_threshold[0.70])
    correlations = {}
    for t in THRESHOLDS:
        if t == 0.70:
            continue
        rho, p = spearmanr(cov_baseline, cov_by_threshold[t])
        correlations[str(t)] = {'spearman_rho_vs_0.70': float(rho), 'p': float(p)}

    result = {
        'description': ('Coverage (C_tau) recomputed at 5 fidelity thresholds for every '
                         'point on the pooled evolutionary front, by re-decoding each '
                         'stored parameter vector and re-running the canonical exact '
                         'path search (not a proxy). Tests whether the canonical '
                         'threshold (0.70) is a load-bearing, arbitrary choice.'),
        'thresholds': THRESHOLDS,
        'n_points': n,
        'coverage_mean_by_threshold': {str(t): float(np.mean(cov_by_threshold[t])) for t in THRESHOLDS},
        'coverage_std_by_threshold': {str(t): float(np.std(cov_by_threshold[t], ddof=1)) for t in THRESHOLDS},
        'spearman_rank_correlation_vs_baseline_0.70': correlations,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    for t in THRESHOLDS:
        print(f'  threshold={t}: mean C_tau={result["coverage_mean_by_threshold"][str(t)]:.4f}')
    for t, c in correlations.items():
        print(f'  Spearman(0.70, {t}) = {c["spearman_rho_vs_0.70"]:.4f}')


if __name__ == '__main__':
    main()
