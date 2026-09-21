#!/usr/bin/env python3
"""
routing_policy_sensitivity_v3.py — is the pooled front's low effective
dimensionality (Section: Discussion, PC1 ~96% of variance) an artefact
of always evaluating fidelity, rate AND latency along the SAME
fidelity-optimal path, rather than each objective's own optimal path?
An external review raised this directly: if rate/latency were instead
measured along their own best route, the strong cross-objective
correlations that drive the PCA result might not hold.

Re-decodes every point on the pooled evolutionary front (from its
stored parameter vector) and recomputes the four aggregate objectives
under two ADDITIONAL routing policies, reusing the canonical oracle's
own adjacency construction (medium-collapse-to-highest-fidelity) so
only the path-selection rule differs:

  - fidelity-optimal (canonical, already used throughout this paper):
    maximise the exact Werner-composed end-to-end fidelity hop-by-hop.
  - rate-optimal: maximise the bottleneck (minimum) per-hop rate along
    the path (classic widest-path problem under a different edge
    weight), a Dijkstra-like relaxation exactly analogous to the
    canonical search but composing via min() instead of Werner
    composition.
  - latency-optimal: minimise the sum of per-hop latencies (standard
    additive shortest path).

For each policy, fidelity/rate/latency are all read off along THAT
policy's own selected path (not the fidelity-optimal one), and the
correlation matrix + PCA are recomputed exactly as in
analyze_effective_dimensionality.py, so the three policies' results
are directly comparable.
"""
import heapq
import json
from pathlib import Path

import numpy as np

from bdcz_oracle_v3 import (
    N_NODES, decode_topology, _build_best_medium_adjacency, _werner_compose,
)

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT_JSON = RESULTS / 'routing_policy_sensitivity_v3.json'
FID_THRESHOLD = 0.70


def rate_optimal_paths_from_source(adj, source):
    """Widest-path Dijkstra maximising the bottleneck (min) rate along the path."""
    N = N_NODES
    best_R = np.zeros(N)
    best_R[source] = np.inf
    pred = -np.ones(N, dtype=np.int64)
    visited = np.zeros(N, dtype=bool)
    heap = [(-np.inf, source)]
    while heap:
        neg_r, u = heapq.heappop(heap)
        if visited[u]:
            continue
        visited[u] = True
        r_u = -neg_r
        for v, _f, r_uv, _l in adj[u]:
            if visited[v]:
                continue
            cand = min(r_u, r_uv)
            if cand > best_R[v]:
                best_R[v] = cand
                pred[v] = u
                heapq.heappush(heap, (-cand, v))
    return best_R, pred


def latency_optimal_paths_from_source(adj, source):
    """Standard additive-shortest-path Dijkstra minimising summed latency."""
    N = N_NODES
    best_L = np.full(N, np.inf)
    best_L[source] = 0.0
    pred = -np.ones(N, dtype=np.int64)
    visited = np.zeros(N, dtype=bool)
    heap = [(0.0, source)]
    while heap:
        l_u, u = heapq.heappop(heap)
        if visited[u]:
            continue
        visited[u] = True
        for v, _f, _r, l_uv in adj[u]:
            if visited[v]:
                continue
            cand = l_u + l_uv
            if cand < best_L[v]:
                best_L[v] = cand
                pred[v] = u
                heapq.heappush(heap, (cand, v))
    return best_L, pred


def path_nodes_from_pred(pred, src, tgt):
    nodes = []
    cur = tgt
    while cur != src and cur >= 0:
        nodes.append(cur)
        cur = pred[cur]
    nodes.append(src)
    nodes.reverse()
    return nodes


def evaluate_along_path(adj_lookup, path_nodes, fiber_swap_success=0.50):
    """Given a path's node sequence, compute (F_e2e, rate_e2e, lat_e2e)
    by walking the SAME edges that path uses, regardless of which
    policy selected it."""
    hop_fids, hop_rates, hop_lats = [], [], []
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        f, r, l = adj_lookup.get((a, b), adj_lookup.get((b, a)))
        hop_fids.append(f)
        hop_rates.append(r)
        hop_lats.append(l)
    F_e2e = hop_fids[0]
    for f in hop_fids[1:]:
        F_e2e = _werner_compose(F_e2e, f)
    n_hops = len(hop_rates)
    rate_e2e = (min(hop_rates) if hop_rates else 1e-3) * (fiber_swap_success ** max(n_hops - 1, 0))
    lat_e2e = sum(hop_lats)
    return F_e2e, rate_e2e, lat_e2e


def build_adj_lookup(adj):
    lookup = {}
    for u, neighbours in enumerate(adj):
        for v, f, r, l in neighbours:
            lookup[(u, v)] = (f, r, l)
    return lookup


def evaluate_policy(eu, ev, et, policy):
    N = N_NODES
    total_pairs = N * (N - 1) // 2
    adj = _build_best_medium_adjacency(eu, ev, et)
    adj_lookup = build_adj_lookup(adj)

    n_pairs = 0
    sum_fid = sum_ratelog = sum_lat = 0.0
    n_covered = 0

    for src in range(N):
        if policy == 'rate':
            best_val, pred = rate_optimal_paths_from_source(adj, src)
            reachable = lambda tgt: best_val[tgt] > 0 and pred[tgt] >= 0
        else:
            best_val, pred = latency_optimal_paths_from_source(adj, src)
            reachable = lambda tgt: np.isfinite(best_val[tgt]) and pred[tgt] >= 0

        for tgt in range(src + 1, N):
            if not reachable(tgt):
                continue
            path_nodes = path_nodes_from_pred(pred, src, tgt)
            F_e2e, rate_e2e, lat_e2e = evaluate_along_path(adj_lookup, path_nodes)
            n_pairs += 1
            sum_fid += F_e2e
            sum_ratelog += np.log1p(rate_e2e)
            sum_lat += lat_e2e
            if F_e2e >= FID_THRESHOLD:
                n_covered += 1

    if n_pairs == 0:
        return 0.5, 0.0, 1.0, 0.0
    return sum_fid / n_pairs, sum_ratelog / n_pairs, sum_lat / n_pairs, n_covered / total_pairs


def pca_summary(F):
    Fs = (F - F.mean(axis=0)) / F.std(axis=0)
    cov = np.cov(Fs, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    ratio = eigvals / eigvals.sum()
    corr = np.corrcoef(F, rowvar=False)
    return ratio.tolist(), corr.tolist()


def main():
    gp = np.load(RESULTS / 'global_pareto_v3.npz')
    X = gp['X_pareto']
    n = len(X)
    print(f'Recomputing objectives under 2 alternative routing policies for {n} points...')

    results = {'fidelity_optimal_canonical': gp['F_pareto'].tolist()}
    for policy in ['rate', 'latency']:
        F_policy = []
        for i, params in enumerate(X):
            decoded = decode_topology(params)
            if decoded is None:
                F_policy.append([-0.5, 0.0, 1.0, 0.0])
                continue
            eu, ev, et = decoded
            f, r, l, c = evaluate_policy(eu, ev, et, policy)
            F_policy.append([-f, -r, l, -c])
            if (i + 1) % 100 == 0:
                print(f'  [{policy}] {i+1}/{n}')
        results[f'{policy}_optimal'] = F_policy

    out = {'description': ('Fidelity/rate/latency/coverage recomputed under 3 routing '
                            'policies (fidelity-optimal/canonical, rate-optimal, '
                            'latency-optimal) for every pooled-front point, to test '
                            'whether the low effective dimensionality (PCA) is an '
                            'artefact of always using the fidelity-optimal path.'),
           'n_points': n}
    for policy_key in ['fidelity_optimal_canonical', 'rate_optimal', 'latency_optimal']:
        F = np.array(results[policy_key])
        # F is stored in minimisation convention (-fid, -ratelog, lat, -cov); flip for PCA/corr readability
        F_readable = np.column_stack([-F[:, 0], -F[:, 1], F[:, 2], -F[:, 3]])
        ratio, corr = pca_summary(F_readable)
        out[policy_key] = {'pca_explained_variance_ratio': ratio, 'correlation_matrix': corr}

    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    for policy_key in ['fidelity_optimal_canonical', 'rate_optimal', 'latency_optimal']:
        pc1 = out[policy_key]['pca_explained_variance_ratio'][0]
        print(f'  {policy_key}: PC1 = {pc1*100:.1f}%')


if __name__ == '__main__':
    main()
