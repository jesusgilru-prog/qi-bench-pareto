#!/usr/bin/env python3
"""
connectivity_diagnostics_v3.py — structural connectivity of the
Pareto-front topologies.

External review (finding: "the formulation favours disconnected and
potentially degenerate topologies") noted that mean fidelity, rate and
latency are averaged only over CONNECTED city pairs, so hard pairs
vanish from those three objectives and are penalised only through
coverage. A solution can therefore sit on the front with very few
edges and near-zero coverage. The manuscript reported edge counts and
coverage but never reported how CONNECTED the front's topologies
actually are.

This script computes, for each of the 659 front points:
  * number of connected components (over the 100 cities),
  * size of the giant component,
  * fraction of the 4,950 city pairs that are reachable at all,
  * fraction of pairs that are reachable AND meet the fidelity
    threshold (i.e. the coverage objective itself),
so the manuscript can state plainly what share of the released front
consists of operationally global networks and what share does not.

Connectivity is computed on the same best-medium adjacency the oracle
uses for path search, so "reachable" here means exactly "reachable by
the oracle", not "has an edge in the stored topology".
"""
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import bdcz_oracle_v3 as O

RESULTS = Path(__file__).resolve().parent.parent / 'results'
FRONT_NPZ = RESULTS / 'global_pareto_v3.npz'
OUT_JSON = RESULTS / 'connectivity_diagnostics_v3.json'
OUT_NPZ = RESULTS / 'connectivity_diagnostics_v3.npz'

FID_THRESHOLD = 0.70   # same as oracle_evaluate()


def diagnostics_for(x):
    res = O.decode_topology(np.asarray(x, dtype=float))
    if res is None:
        return None
    eu, ev, et = res
    adj = O._build_best_medium_adjacency(eu, ev, et)
    n = O.N_NODES
    rows, cols = [], []
    for u in range(n):
        for v, _f, _r, _l in adj[u]:
            rows.append(u)
            cols.append(v)
    if not rows:
        return None
    g = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    n_comp, labels = connected_components(g, directed=False)
    _sizes, counts = np.unique(labels, return_counts=True)
    giant = int(counts.max())

    total_pairs = n * (n - 1) // 2
    reachable = 0
    covered = 0
    for s in range(n):
        best_F, pred = O._werner_optimal_paths_from_source(adj, s)
        for t in range(s + 1, n):
            if best_F[t] <= 0.0 or pred[t] < 0:
                continue
            reachable += 1
            if best_F[t] >= FID_THRESHOLD:
                covered += 1
    return {
        'n_components': int(n_comp),
        'giant_component_size': giant,
        'n_isolated_nodes': int((counts == 1).sum()),
        'reachable_pairs': reachable,
        'reachable_fraction': reachable / total_pairs,
        'covered_pairs': covered,
        'covered_fraction': covered / total_pairs,
    }


def main():
    d = np.load(FRONT_NPZ, allow_pickle=True)
    X = d['X_pareto']
    F = d['F_pareto']
    n_edges = d['n_edges']

    recs = []
    for i in range(X.shape[0]):
        r = diagnostics_for(X[i])
        if r is None:
            r = {'n_components': O.N_NODES, 'giant_component_size': 1,
                 'n_isolated_nodes': O.N_NODES, 'reachable_pairs': 0,
                 'reachable_fraction': 0.0, 'covered_pairs': 0, 'covered_fraction': 0.0}
        r['n_edges'] = int(n_edges[i])
        recs.append(r)

    comps = np.array([r['n_components'] for r in recs])
    giant = np.array([r['giant_component_size'] for r in recs])
    reach = np.array([r['reachable_fraction'] for r in recs])
    cover = np.array([r['covered_fraction'] for r in recs])
    edges = np.array([r['n_edges'] for r in recs])

    np.savez_compressed(OUT_NPZ, n_components=comps, giant_component_size=giant,
                        reachable_fraction=reach, covered_fraction=cover, n_edges=edges)

    def pct(mask):
        return float(mask.sum()) / len(recs)

    out = {
        'description': ('Structural connectivity of the 659 deduplicated Pareto-front '
                        'points, computed on the same best-medium adjacency the oracle '
                        'uses for path search. "Reachable" = some path exists (the set '
                        'the fidelity/rate/latency means are taken over); "covered" = '
                        'reachable AND end-to-end fidelity >= %.2f (the coverage '
                        'objective).' % FID_THRESHOLD),
        'n_points': len(recs),
        'n_cities': int(O.N_NODES),
        'total_city_pairs': int(O.N_NODES * (O.N_NODES - 1) // 2),
        'fidelity_threshold': FID_THRESHOLD,
        'connected_components': {
            'min': int(comps.min()), 'median': float(np.median(comps)),
            'mean': float(comps.mean()), 'max': int(comps.max()),
            'fraction_fully_connected': pct(comps == 1),
        },
        'giant_component_size': {
            'min': int(giant.min()), 'median': float(np.median(giant)),
            'mean': float(giant.mean()), 'max': int(giant.max()),
        },
        'reachable_fraction': {
            'min': float(reach.min()), 'median': float(np.median(reach)),
            'mean': float(reach.mean()), 'max': float(reach.max()),
        },
        'covered_fraction': {
            'min': float(cover.min()), 'median': float(np.median(cover)),
            'mean': float(cover.mean()), 'max': float(cover.max()),
        },
        'operational_bands': {
            'fraction_reachable_ge_0.99': pct(reach >= 0.99),
            'fraction_reachable_ge_0.90': pct(reach >= 0.90),
            'fraction_reachable_ge_0.50': pct(reach >= 0.50),
            'fraction_reachable_lt_0.10': pct(reach < 0.10),
            'fraction_coverage_lt_0.01': pct(cover < 0.01),
        },
        'sparse_solutions': {
            'note': ('front points with fewer than 20 stored edges -- the regime the '
                     'review flagged as "17-18 edges, coverage near zero"'),
            'count': int((edges < 20).sum()),
            'median_reachable_fraction': (float(np.median(reach[edges < 20]))
                                          if (edges < 20).any() else None),
            'median_components': (float(np.median(comps[edges < 20]))
                                  if (edges < 20).any() else None),
            'max_covered_fraction': (float(cover[edges < 20].max())
                                     if (edges < 20).any() else None),
        },
        'edges_vs_reachable_spearman': None,
    }
    from scipy.stats import spearmanr
    rho, p = spearmanr(edges, reach)
    out['edges_vs_reachable_spearman'] = {'rho': float(rho), 'p': float(p)}

    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    print(json.dumps({k: v for k, v in out.items() if k != 'description'}, indent=2))


if __name__ == '__main__':
    main()
