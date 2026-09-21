#!/usr/bin/env python3
"""
werner_path_validation_v3.py — validation of the oracle's maximum
Werner-composed-fidelity path search.

External review (finding: "the mathematical explanation of the
Werner-optimal route is misleading and incomplete") observed that under
the Werner-parameter substitution

    w = (4F - 1) / 3        <=>        F = (3w + 1) / 4

the swap composition F12 = F1 F2 + (1 - F1)(1 - F2)/3 becomes an exact
PRODUCT, w12 = w1 w2. Consequently:

  * an EXACT additive shortest-path formulation exists, with edge
    weights -log w_e (valid whenever 0 < w_e <= 1);
  * what is only a proxy is -log F_e, which the manuscript previously
    conflated with "no additive formulation exists";
  * optimality of the Dijkstra-like relaxation follows from w being a
    positive product with w_e <= 1 (path fidelity is non-increasing in
    additional hops and the sub-path of an optimal path is optimal),
    not merely from "monotonically non-increasing composition".

This script proves the algebra numerically and validates the
implementation three ways:

  1. the substitution identity w12 = w1 w2 holds to machine precision;
  2. on the paper's own 100-city instance, the oracle's
     _werner_optimal_paths_from_source() returns exactly the same
     end-to-end fidelities as a shortest-path search with -log w
     weights (scipy dijkstra);
  3. on small random graphs, both agree with exhaustive enumeration of
     all simple paths.

It also quantifies how often the -log F proxy (which the manuscript
correctly declines to use) selects a Werner-suboptimal path.
"""
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

import bdcz_oracle_v3 as O

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT_JSON = RESULTS / 'werner_path_validation_v3.json'

RNG_SEED = 20260730
N_RANDOM_GRAPHS = 200
RANDOM_GRAPH_NODES = 7
TOL = 1e-9


def w_of(F):
    return (4.0 * np.asarray(F) - 1.0) / 3.0


def F_of(w):
    return (3.0 * np.asarray(w) + 1.0) / 4.0


def check_substitution_identity(rng, n=200000):
    """w12 = w1 * w2 exactly, for the Werner composition used by the oracle."""
    F1 = rng.uniform(0.25, 1.0, n)
    F2 = rng.uniform(0.25, 1.0, n)
    F12 = O._werner_compose(F1, F2)
    lhs = w_of(F12)
    rhs = w_of(F1) * w_of(F2)
    return float(np.abs(lhs - rhs).max())


def adjacency_to_matrix(adj, n):
    """Dense F matrix (0 where no edge) from the oracle's adjacency list."""
    M = np.zeros((n, n))
    for u in range(n):
        for v, f, _r, _l in adj[u]:
            M[u, v] = max(M[u, v], f)
    return M


def neglogw_shortest_path(Fmat):
    """End-to-end fidelity from an exact additive shortest path on -log w."""
    n = Fmat.shape[0]
    W = w_of(Fmat)
    valid = (Fmat > 0) & (W > 0)
    cost = np.zeros((n, n))
    cost[valid] = -np.log(W[valid])
    sp = dijkstra(csr_matrix(np.where(valid, np.maximum(cost, 0.0), 0.0) * valid),
                  directed=False, unweighted=False)
    return F_of(np.exp(-sp))


def neglogF_shortest_path(Fmat):
    """The (approximate) -log F proxy, for the suboptimality count only."""
    n = Fmat.shape[0]
    valid = Fmat > 0
    cost = np.zeros((n, n))
    cost[valid] = -np.log(Fmat[valid])
    sp = dijkstra(csr_matrix(np.where(valid, np.maximum(cost, 0.0), 0.0) * valid),
                  directed=False, unweighted=False)
    return sp, valid


def exhaustive_best_F(Fmat, s, t):
    """Best Werner-composed fidelity over all simple s-t paths."""
    n = Fmat.shape[0]
    others = [x for x in range(n) if x not in (s, t)]
    best = 0.0
    for r in range(n - 1):
        for mid in itertools.permutations(others, r):
            path = (s,) + mid + (t,)
            f = 1.0
            ok = True
            for a, b in zip(path[:-1], path[1:]):
                if Fmat[a, b] <= 0:
                    ok = False
                    break
                f = O._werner_compose(f, Fmat[a, b])
            if ok:
                best = max(best, f)
    return best


def oracle_paths(adj, n):
    out = np.zeros((n, n))
    for s in range(n):
        bf, _pred = O._werner_optimal_paths_from_source(adj, s)
        out[s] = bf
    return out


def validate_on_paper_instance(rng):
    """Decode a few topologies from the paper's own bounds and compare."""
    n = O.N_NODES
    max_abs_diff = 0.0
    proxy_suboptimal_pairs = 0
    pairs_checked = 0
    topologies = 0
    for trial in range(8):
        params = np.array([
            rng.integers(2, 9), rng.integers(1, 5), rng.integers(0, 3),
            rng.uniform(1500, 6000), rng.uniform(2000, 12000), rng.uniform(200, 2000),
            rng.uniform(0.3, 0.9), rng.uniform(0.3, 0.9), rng.uniform(0.3, 0.9),
            rng.uniform(0.0, 1.0), trial,
        ], dtype=float)
        res = O.decode_topology(params, rng=np.random.default_rng(5000 + trial))
        if res is None:
            continue
        topologies += 1
        eu, ev, et = res
        adj = O._build_best_medium_adjacency(eu, ev, et)
        Fmat = adjacency_to_matrix(adj, n)

        F_oracle = oracle_paths(adj, n)
        F_additive = neglogw_shortest_path(Fmat)

        reach = (F_oracle > 0) & ~np.eye(n, dtype=bool)
        max_abs_diff = max(max_abs_diff, float(np.abs(F_oracle[reach] - F_additive[reach]).max()))

        # how often does the -log F proxy pick a strictly worse path?
        sp_proxy, valid = neglogF_shortest_path(Fmat)
        for s in range(n):
            _bf, _pred = O._werner_optimal_paths_from_source(adj, s)
            for t in range(s + 1, n):
                if not reach[s, t]:
                    continue
                pairs_checked += 1
                if not np.isfinite(sp_proxy[s, t]):
                    continue
                # re-compose the proxy's chosen path fidelity
                path = proxy_path(Fmat, sp_proxy, s, t)
                if path is None:
                    continue
                f = 1.0
                for a, b in zip(path[:-1], path[1:]):
                    f = O._werner_compose(f, Fmat[a, b])
                if f < F_oracle[s, t] - 1e-12:
                    proxy_suboptimal_pairs += 1
    return {
        'topologies_checked': topologies,
        'reachable_pairs_checked': pairs_checked,
        'max_abs_fidelity_difference_oracle_vs_neglogw': max_abs_diff,
        'proxy_neglogF_suboptimal_pairs': proxy_suboptimal_pairs,
        'proxy_neglogF_suboptimal_fraction': (proxy_suboptimal_pairs / pairs_checked
                                              if pairs_checked else 0.0),
    }


def proxy_path(Fmat, sp_proxy, s, t):
    """Reconstruct one -log F shortest path by greedy backward stepping."""
    n = Fmat.shape[0]
    path = [t]
    cur = t
    guard = 0
    while cur != s and guard < n + 2:
        guard += 1
        best = None
        for v in range(n):
            if Fmat[cur, v] <= 0:
                continue
            if abs(sp_proxy[s, v] + (-np.log(Fmat[cur, v])) - sp_proxy[s, cur]) < 1e-9:
                best = v
                break
        if best is None:
            return None
        path.append(best)
        cur = best
    if cur != s:
        return None
    return path[::-1]


def validate_on_random_graphs(rng):
    """Exhaustive enumeration vs oracle vs -log w on small graphs."""
    n = RANDOM_GRAPH_NODES
    max_diff_oracle = 0.0
    max_diff_additive = 0.0
    checked = 0
    for _ in range(N_RANDOM_GRAPHS):
        Fmat = np.zeros((n, n))
        for u in range(n):
            for v in range(u + 1, n):
                if rng.random() < 0.45:
                    f = rng.uniform(0.6, 0.99)
                    Fmat[u, v] = Fmat[v, u] = f
        adj = [[] for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if Fmat[u, v] > 0:
                    adj[u].append((v, Fmat[u, v], 1.0, 1.0))

        saved_n, saved_adjfun = O.N_NODES, None
        O.N_NODES = n
        try:
            F_oracle = np.zeros((n, n))
            for s in range(n):
                bf, _ = O._werner_optimal_paths_from_source(adj, s)
                F_oracle[s] = bf
        finally:
            O.N_NODES = saved_n
        F_additive = neglogw_shortest_path(Fmat)

        for s in range(n):
            for t in range(s + 1, n):
                if F_oracle[s, t] <= 0:
                    continue
                brute = exhaustive_best_F(Fmat, s, t)
                checked += 1
                max_diff_oracle = max(max_diff_oracle, abs(F_oracle[s, t] - brute))
                max_diff_additive = max(max_diff_additive, abs(F_additive[s, t] - brute))
    return {
        'n_graphs': N_RANDOM_GRAPHS,
        'nodes_per_graph': n,
        'pairs_checked': checked,
        'max_abs_diff_oracle_vs_exhaustive': max_diff_oracle,
        'max_abs_diff_neglogw_vs_exhaustive': max_diff_additive,
    }


def main():
    rng = np.random.default_rng(RNG_SEED)
    identity_err = check_substitution_identity(rng)
    paper = validate_on_paper_instance(rng)
    small = validate_on_random_graphs(rng)

    out = {
        'description': ('Validation of the maximum Werner-composed-fidelity path search. '
                        'Under w=(4F-1)/3 the swap composition is an exact product '
                        '(w12=w1*w2), so an exact additive shortest-path formulation with '
                        'weights -log w exists; -log F is the proxy that does NOT preserve '
                        'optimality. Checks: (1) the substitution identity; (2) agreement '
                        'between the oracle\'s Dijkstra-like relaxation and a -log w '
                        'shortest path on the paper\'s own 100-city instance; (3) agreement '
                        'of both with exhaustive enumeration on small random graphs.'),
        'rng_seed': RNG_SEED,
        'substitution_identity_max_abs_error': identity_err,
        'paper_instance': paper,
        'small_random_graphs': small,
        'tolerance': TOL,
        'all_checks_pass': bool(
            identity_err < 1e-12
            and paper['max_abs_fidelity_difference_oracle_vs_neglogw'] < TOL
            and small['max_abs_diff_oracle_vs_exhaustive'] < TOL
            and small['max_abs_diff_neglogw_vs_exhaustive'] < TOL
        ),
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    print(json.dumps({k: v for k, v in out.items() if k != 'description'}, indent=2))


if __name__ == '__main__':
    main()
