#!/usr/bin/env python3
"""
classical_baselines_v3.py — Canonical-oracle rerun of classical baselines
==============================================================================
Same 14 constructors as v2 (MST, k-NN, RGG, hub-and-spoke, region-aware
MST) with the same heterogeneous BEST_MEDIUM per-edge assignment, but
using bdcz_oracle_v3 (canonical, exact Werner-optimal search, real
SAT/FSO constants matching Paper 1).
"""
import numpy as np
import networkx as nx
import json
from pathlib import Path

from bdcz_oracle_v3 import oracle_evaluate, DIST, N_NODES, city_region, BEST_MEDIUM


def eval_topo(edges):
    if len(edges) < 5:
        return None
    eu = np.array([e[0] for e in edges], dtype=np.int64)
    ev = np.array([e[1] for e in edges], dtype=np.int64)
    et = BEST_MEDIUM[eu, ev]
    try:
        fid, rate_log, lat, cov = oracle_evaluate(eu, ev, et)
        return {
            'fidelity': float(fid), 'rate_log': float(rate_log),
            'latency': float(lat), 'coverage': float(cov),
            'F_minimization': [float(-fid), float(-rate_log), float(lat), float(-cov)],
        }
    except Exception as e:
        print(f"  Oracle error: {e}")
        return None


print("=" * 60); print("1. MST Geografico (Kruskal)"); print("=" * 60)
G_full = nx.Graph()
for i in range(N_NODES):
    for j in range(i + 1, N_NODES):
        G_full.add_edge(i, j, weight=DIST[i, j])
mst = nx.minimum_spanning_tree(G_full)
mst_edges = [(min(u, v), max(u, v)) for u, v in mst.edges()]
mst_result = eval_topo(mst_edges)
print(f"MST: |E|={len(mst_edges)}, metrics={mst_result}")

print("\n" + "=" * 60); print("2. k-NN Determinista"); print("=" * 60)
knn_results = {}
for k in [3, 5, 7, 10]:
    edges = set()
    for i in range(N_NODES):
        nbrs = np.argsort(DIST[i])[1:k + 1]
        for j in nbrs:
            edges.add((min(i, j), max(i, j)))
    edges = sorted(edges)
    result = eval_topo(edges)
    knn_results[f'k{k}'] = {'n_edges': len(edges), **(result or {})}
    print(f"k-NN k={k}: |E|={len(edges)}, metrics={result}")

print("\n" + "=" * 60); print("3. Random Geometric Graph"); print("=" * 60)
rgg_results = {}
for threshold_km in [500, 1000, 2000, 3000, 5000]:
    edges = [(i, j) for i in range(N_NODES) for j in range(i + 1, N_NODES)
             if DIST[i, j] < threshold_km]
    result = eval_topo(edges)
    rgg_results[f'thr{threshold_km}'] = {'n_edges': len(edges), **(result or {})}
    print(f"RGG thr={threshold_km}km: |E|={len(edges)}, metrics={result}")

print("\n" + "=" * 60); print("4. Hub-and-Spoke"); print("=" * 60)
hub_results = {}
for n_hubs in [5, 10, 15]:
    hubs = sorted(range(N_NODES), key=lambda i: mst.degree(i), reverse=True)[:n_hubs]
    hs_edges = set()
    for h in hubs:
        for i in range(N_NODES):
            if i != h:
                hs_edges.add((min(h, i), max(h, i)))
    for i, h1 in enumerate(hubs):
        for h2 in hubs[i + 1:]:
            hs_edges.add((min(h1, h2), max(h1, h2)))
    hs_edges = sorted(hs_edges)
    result = eval_topo(hs_edges)
    hub_results[f'hubs{n_hubs}'] = {'n_edges': len(hs_edges), 'hub_nodes': hubs, **(result or {})}
    print(f"Hub-and-spoke ({n_hubs} hubs): |E|={len(hs_edges)}, metrics={result}")

print("\n" + "=" * 60); print("5. Region-aware MST"); print("=" * 60)
regions = np.unique(city_region)
region_edges = set()
for r in regions:
    nodes_r = np.where(city_region == r)[0]
    if len(nodes_r) < 2:
        continue
    G_r = nx.Graph()
    for i in nodes_r:
        for j in nodes_r:
            if i < j:
                G_r.add_edge(i, j, weight=DIST[i, j])
    mst_r = nx.minimum_spanning_tree(G_r)
    for u, v in mst_r.edges():
        region_edges.add((min(u, v), max(u, v)))
for ri in range(len(regions)):
    for rj in range(ri + 1, len(regions)):
        nodes_i = np.where(city_region == regions[ri])[0]
        nodes_j = np.where(city_region == regions[rj])[0]
        best_d, best_pair = np.inf, None
        for ni in nodes_i:
            for nj in nodes_j:
                if DIST[ni, nj] < best_d:
                    best_d, best_pair = DIST[ni, nj], (min(ni, nj), max(ni, nj))
        if best_pair:
            region_edges.add(best_pair)
region_edges = sorted(region_edges)
region_mst_result = eval_topo(region_edges)
print(f"Region-aware MST: |E|={len(region_edges)}, metrics={region_mst_result}")

out = {
    'mst': {'n_edges': len(mst_edges), **(mst_result or {})},
    'knn': knn_results, 'rgg': rgg_results, 'hub_and_spoke': hub_results,
    'region_mst': {'n_edges': len(region_edges), **(region_mst_result or {})},
    'note': 'v3: canonical Paper-1-unified oracle (exact Werner-optimal search); BEST_MEDIUM per edge',
}
out_path = Path(__file__).resolve().parent.parent / 'results' / 'classical_baselines_v3_results.json'
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\nSaved: {out_path}")
