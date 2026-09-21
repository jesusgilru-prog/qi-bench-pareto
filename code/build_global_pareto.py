#!/usr/bin/env python3
"""
build_global_pareto.py — Global non-dominated Pareto front, all 7 algorithms
================================================================================
Global non-dominated union of NSGA-II + NSGA-III + MOEA/D + SMS-EMOA + RVEA +
SPEA2 + AGE-MOEA-II across all 30 seeds each (canonical v3 oracle), tagged
with source_algorithm. Also
reports per-solution resource use (edge count and medium composition),
addressing the reviewer's point that Pareto quality claims should not be
read without knowing whether solutions simply use more infrastructure
than the classical baselines they are compared against.

FIX (Round 3 review, confirmed real bug): the non-dominated union of 1001
raw records contained massive duplication across seeds/algorithms -- only
494 UNIQUE objective vectors, and MOEA/D's 278 raw records collapsed to
just 19 unique objective vectors (independently verified against the npz
before this fix: np.unique(F.round(6), axis=0)). Reporting "1,001 global
non-dominated topologies" with "MOEA/D contributes 28%" is therefore
misleading -- it counts the same handful of MOEA/D solutions, rediscovered
across different seeds, as if they were 278 distinct contributions. Fixed
by deduplicating on the OBJECTIVE VECTOR (rounded to 1e-6, the natural unit
for a Pareto front, which is defined in objective space) as the primary
reported front, while separately reporting the number of structurally
unique DECODED TOPOLOGIES (distinct edge sets) among the raw pooled
records, since two different topologies can occasionally map to
numerically identical objective vectors. Per-point provenance (every
(algorithm, seed) pair that produced a given unique objective vector, not
just the first one found) is preserved in global_pareto_v3_provenance.json
rather than being silently discarded by the dedup.
"""
import numpy as np
import pickle
import json
from pathlib import Path
from collections import Counter, defaultdict

from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from bdcz_oracle_v3 import decode_topology, _build_best_medium_adjacency

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step1_baselines_v3'

#  round 6 (external review B1/B2): 10->30 seeds, 3->5 algorithms
#  (added SMS-EMOA, RVEA), matching run_baselines_comparison_v3.py exactly.
#  round 14 (external review, ChatGPT Round 13 finding #1, confirmed real
#  bug): SPEA2 and AGE-MOEA-II were added as a seventh/eighth-algorithm
#  comparison in unified_seven_algorithm_comparison_v3.py (561-point
#  reference front, config in that script's own JSON) but this script kept
#  pooling only the original 5, so every downstream consumer of
#  global_pareto_v3.npz -- connectivity, PCA, threshold/routing
#  sensitivity, classical-baseline dominance, Figure 4, the released
#  package -- was still analysing the OLD 513-point front even after the
#  manuscript's headline comparison moved to 7 algorithms. Now pools all 7,
#  matching unified_seven_algorithm_comparison_v3.py's ALGOS/CKPT routing
#  exactly (MODERN_CKPT for spea2/age_moea2, MAIN_CKPT for the rest).
MAIN_ALGORITHMS = ['nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea']
MODERN_ALGORITHMS = ['spea2', 'age_moea2']
ALGORITHMS = MAIN_ALGORITHMS + MODERN_ALGORITHMS
MODERN_CKPT_DIR = RESULTS / 'step9_modern_algorithms_v3'
SEEDS = list(range(42, 72))


def ckpt_dir_for(algo):
    return MODERN_CKPT_DIR if algo in MODERN_ALGORITHMS else CKPT_DIR


all_F, all_X, all_algo, all_seed = [], [], [], []
for algo in ALGORITHMS:
    for seed in SEEDS:
        p = ckpt_dir_for(algo) / f'{algo}_seed{seed}.pkl'
        with open(p, 'rb') as fh:
            run = pickle.load(fh)
        F, X = run['F_pareto'], run['X_pareto']
        if len(F) == 0:
            continue
        all_F.append(F)
        all_X.append(X)
        all_algo.extend([algo] * len(F))
        all_seed.extend([seed] * len(F))

merged_F = np.vstack(all_F)
merged_X = np.vstack(all_X)
algo_arr = np.array(all_algo)
seed_arr = np.array(all_seed)

nd_mask = find_non_dominated(merged_F)
raw_F = merged_F[nd_mask]
raw_X = merged_X[nd_mask]
raw_algo = algo_arr[nd_mask]
raw_seed = seed_arr[nd_mask]

print(f'Total points pooled: {len(merged_F)}')
print(f'Global non-dominated RAW records (pre-dedup): {len(raw_F)}')

# ── Deduplicate by objective vector (the Pareto front's natural unit) ──────
key_to_indices = defaultdict(list)
for i, f in enumerate(raw_F):
    key_to_indices[tuple(np.round(f, 6))].append(i)

unique_keys = list(key_to_indices.keys())
rep_idx = np.array([key_to_indices[k][0] for k in unique_keys])  # first occurrence per unique vector
global_F = raw_F[rep_idx]
global_X = raw_X[rep_idx]
global_algo = raw_algo[rep_idx]        # primary (first-seen) source, kept for backward-compat plotting
global_seed = raw_seed[rep_idx]

# Full provenance: every (algorithm, seed) that produced each unique vector
provenance = []
exclusive_counts = Counter()
shared_count = 0
for k in unique_keys:
    idxs = key_to_indices[k]
    sources = sorted(set(zip(raw_algo[idxs].tolist(), raw_seed[idxs].tolist())))
    algos_present = sorted(set(a for a, s in sources))
    provenance.append({'objective_vector': list(k),
                        'sources': [{'algorithm': a, 'seed': int(s)} for a, s in sources]})
    if len(algos_present) == 1:
        exclusive_counts[algos_present[0]] += 1
    else:
        shared_count += 1

print(f'Unique objective vectors (deduplicated Pareto front): {len(unique_keys)}')
print(f'  exclusive to one algorithm: { {a: exclusive_counts.get(a, 0) for a in ALGORITHMS} }')
print(f'  shared across >=2 algorithms: {shared_count}')

with open(RESULTS / 'global_pareto_v3_provenance.json', 'w') as f:
    json.dump({'n_unique_objective_vectors': len(unique_keys),
                'n_raw_records_pre_dedup': int(len(raw_F)),
                'exclusive_by_algorithm': {a: exclusive_counts.get(a, 0) for a in ALGORITHMS},
                'shared_across_2plus_algorithms': shared_count,
                'per_point_provenance': provenance}, f, indent=2)
print(f'Saved: global_pareto_v3_provenance.json')

# ── Structurally unique DECODED TOPOLOGIES among the raw (pre-dedup) pool,
#    reported separately since a Pareto FRONT (objective space) and a
#    Pareto SET (design space) are conceptually distinct quantities.
#    FIX (Round 4 review, confirmed real): the package previously only
#    released 452 representative X vectors (one per unique objective
#    vector), while separately CLAIMING 599 (now: whatever this rerun's
#    count is) structurally distinct topologies -- those topologies were
#    never actually written to disk, only counted. Every unique topology's
#    full structure (edge list, medium, decoder params, objective vector,
#    and every (algorithm, seed) that (re)discovered it) is now saved to
#    unique_topologies_v3.npz/.json, so the released package matches what
#    the manuscript claims it releases. ─────────────────────────────────
topo_first_idx = {}          # key -> first raw index that produced it
topo_by_algo = defaultdict(set)
topo_provenance = defaultdict(list)
topo_eu, topo_ev, topo_et = {}, {}, {}
for i in range(len(raw_X)):
    result = decode_topology(raw_X[i])
    if result is None:
        continue
    eu, ev, et = result
    key = tuple(sorted(zip(eu.tolist(), ev.tolist(), et.tolist())))
    if key not in topo_first_idx:
        topo_first_idx[key] = i
        topo_eu[key], topo_ev[key], topo_et[key] = eu, ev, et
    topo_by_algo[raw_algo[i]].add(key)
    topo_provenance[key].append({'algorithm': str(raw_algo[i]), 'seed': int(raw_seed[i])})

topo_keys = list(topo_first_idx.keys())
n_unique_topologies = len(topo_keys)
topo_exclusive = Counter()
topo_shared = 0
for k in topo_keys:
    owners = [a for a in ALGORITHMS if k in topo_by_algo[a]]
    if len(owners) == 1:
        topo_exclusive[owners[0]] += 1
    else:
        topo_shared += 1
print(f'Structurally unique decoded topologies (raw pool, pre-dominance-pooling context): '
      f'{n_unique_topologies}')
print(f'  by algorithm (count of unique topologies each algorithm contributes to): '
      f'{ {a: len(topo_by_algo[a]) for a in ALGORITHMS} }')
print(f'  exclusive: { {a: topo_exclusive.get(a,0) for a in ALGORITHMS} }, shared: {topo_shared}')

# Save every unique topology's structure, not just a count of them.
# FIX (external review round 7, author-approved full scope, "safe
# serialisation"): previously stored eu/ev/et as dtype=object arrays
# (one variable-length array per topology), which requires
# allow_pickle=True to load -- a real portability/security concern for
# publicly-released data. Flattened into a single (edge_u, edge_v,
# edge_type) triple plus a topology_offsets index array (CSR-style: edge
# k of topology i is at flat index topology_offsets[i]+k, i's edges span
# [topology_offsets[i], topology_offsets[i+1])), all plain numeric
# dtypes. Loads with allow_pickle=False; the trusted-internal
# .pkl checkpoints elsewhere in results/ are unaffected by this change.
topo_idx_order = [topo_first_idx[k] for k in topo_keys]
_edge_u_list = [topo_eu[k] for k in topo_keys]
_edge_v_list = [topo_ev[k] for k in topo_keys]
_edge_t_list = [topo_et[k] for k in topo_keys]
_lengths = np.array([len(x) for x in _edge_u_list], dtype=np.int64)
topology_offsets = np.concatenate([[0], np.cumsum(_lengths)]).astype(np.int64)
edge_u = np.concatenate(_edge_u_list).astype(np.int64) if _edge_u_list else np.array([], dtype=np.int64)
edge_v = np.concatenate(_edge_v_list).astype(np.int64) if _edge_v_list else np.array([], dtype=np.int64)
edge_type = np.concatenate(_edge_t_list).astype(np.int64) if _edge_t_list else np.array([], dtype=np.int64)
np.savez(str(RESULTS / 'unique_topologies_v3.npz'),
         X=raw_X[topo_idx_order],
         F=raw_F[topo_idx_order],
         topology_offsets=topology_offsets,
         edge_u=edge_u, edge_v=edge_v, edge_type=edge_type)
with open(RESULTS / 'unique_topologies_v3_provenance.json', 'w') as f:
    json.dump({
        'n_unique_topologies': n_unique_topologies,
        'exclusive_by_algorithm': {a: topo_exclusive.get(a, 0) for a in ALGORITHMS},
        'shared_across_2plus_algorithms': topo_shared,
        'note': ('Row i of unique_topologies_v3.npz (X, F) corresponds to '
                 'provenance[i] below -- every (algorithm, seed) that (re)discovered this '
                 'exact topology. Topology i\'s edges are '
                 'edge_u/edge_v/edge_type[topology_offsets[i]:topology_offsets[i+1]] '
                 '(flat CSR-style storage, allow_pickle=False).'),
        'provenance': [topo_provenance[k] for k in topo_keys],
    }, f, indent=2)
print(f'Saved: unique_topologies_v3.npz, unique_topologies_v3_provenance.json '
      f'({n_unique_topologies} topologies, full structure)')

counts = Counter(global_algo.tolist())
print(f'\nUsing DEDUPLICATED objective-vector front for all downstream figures/tables/dominance checks:')
for a in ALGORITHMS:
    print(f'  source={a:8s}: {counts.get(a, 0)} points '
          f'({100*counts.get(a,0)/len(global_F):.1f}%)')

# ── Resource use per Pareto solution (edge count, medium composition) ──────
# FIX (Round 4 review): |E| is ambiguous across three different quantities
# a topology can report -- (a) raw medium-specific edges (a city pair with
# both a fibre AND a satellite link counts as 2), (b) unique city pairs
# connected by at least one medium, (c) the "effective" edges the oracle
# actually routes over, which collapses parallel media on the same pair to
# the single highest-fidelity one (_build_best_medium_adjacency). All three
# are now reported rather than only (a), since classical baselines (MST/
# kNN/RGG/hub-and-spoke) have exactly one medium per pair by construction,
# so only (b)/(c) are directly comparable to their |E|.
n_edges_list, n_fibre_list, n_sat_list, n_fso_list = [], [], [], []
n_unique_pairs_list, n_effective_list = [], []
for params in global_X:
    result = decode_topology(params)
    if result is None:
        n_edges_list.append(0); n_fibre_list.append(0)
        n_sat_list.append(0); n_fso_list.append(0)
        n_unique_pairs_list.append(0); n_effective_list.append(0)
        continue
    eu, ev, et = result
    n_edges_list.append(len(eu))
    n_fibre_list.append(int((et == 0).sum()))
    n_sat_list.append(int((et == 1).sum()))
    n_fso_list.append(int((et == 2).sum()))
    n_unique_pairs_list.append(len(set(zip(eu.tolist(), ev.tolist()))))
    adj = _build_best_medium_adjacency(eu, ev, et)
    n_effective_list.append(sum(len(nbrs) for nbrs in adj) // 2)
n_edges_arr = np.array(n_edges_list)

np.savez(str(RESULTS / 'global_pareto_v3.npz'),
         F_pareto=global_F, X_pareto=global_X,
         source_algorithm=global_algo, source_seed=global_seed,
         n_edges=n_edges_arr, n_fibre=np.array(n_fibre_list),
         n_satellite=np.array(n_sat_list), n_fso=np.array(n_fso_list),
         n_unique_pairs=np.array(n_unique_pairs_list),
         n_effective_edges=np.array(n_effective_list))

# Human-readable summary
raw_fid = -global_F[:, 0]
raw_ratelog = -global_F[:, 1]
raw_lat = global_F[:, 2]
raw_cov = -global_F[:, 3]
summary = {
    'n_points': int(len(global_F)),
    'n_points_description': 'unique objective vectors (deduplicated Pareto front, objective space)',
    'n_raw_records_pre_dedup': int(len(raw_F)),
    'n_pooled_before_dominance_filter': int(len(merged_F)),
    'n_unique_topologies': int(n_unique_topologies),
    'n_unique_topologies_description': ('structurally distinct decoded edge sets among the '
                                        'raw non-dominated pool (design space, before '
                                        'objective-vector dedup) -- see '
                                        'global_pareto_v3_provenance.json for full per-point '
                                        'algorithm/seed provenance'),
    'by_source_algorithm': {a: int(counts.get(a, 0)) for a in ALGORITHMS},
    'by_source_algorithm_description': ('count of unique objective vectors whose FIRST-SEEN '
                                        '(algo,seed) occurrence was this algorithm, in a fixed '
                                        'iteration order -- kept for backward-compatible '
                                        'plotting only; exclusive_objective_vectors_by_algorithm '
                                        'below is the scientifically meaningful figure (does '
                                        'not depend on iteration order)'),
    'exclusive_objective_vectors_by_algorithm': {a: exclusive_counts.get(a, 0) for a in ALGORITHMS},
    'shared_objective_vectors_2plus_algorithms': shared_count,
    'ranges': {
        'fidelity': [float(raw_fid.min()), float(raw_fid.max())],
        'log_rate': [float(raw_ratelog.min()), float(raw_ratelog.max())],
        'latency_s': [float(raw_lat.min()), float(raw_lat.max())],
        'coverage': [float(raw_cov.min()), float(raw_cov.max())],
    },
    'correlations': {
        'rate_vs_latency': float(np.corrcoef(raw_ratelog, raw_lat)[0, 1]),
        'rate_vs_coverage': float(np.corrcoef(raw_ratelog, raw_cov)[0, 1]),
        'fidelity_vs_rate': float(np.corrcoef(raw_fid, raw_ratelog)[0, 1]),
    },
    'resource_use': {
        'n_edges_raw_medium_specific': {'min': int(n_edges_arr.min()), 'median': float(np.median(n_edges_arr)),
                    'mean': float(n_edges_arr.mean()), 'max': int(n_edges_arr.max()),
                    'description': 'counts each (pair, medium) separately -- a pair with 2 media counts as 2'},
        'n_edges_unique_city_pairs': {'min': int(np.min(n_unique_pairs_list)),
                    'median': float(np.median(n_unique_pairs_list)),
                    'mean': float(np.mean(n_unique_pairs_list)), 'max': int(np.max(n_unique_pairs_list)),
                    'description': 'unique (u,v) pairs connected by at least one medium -- directly comparable to classical baselines, which have exactly one medium per pair'},
        'n_edges_effective_oracle': {'min': int(np.min(n_effective_list)),
                    'median': float(np.median(n_effective_list)),
                    'mean': float(np.mean(n_effective_list)), 'max': int(np.max(n_effective_list)),
                    'description': 'edges actually used for routing after collapsing parallel media per pair to the highest-fidelity one (equals n_edges_unique_city_pairs by construction)'},
        # kept for backward compatibility; equivalent to n_edges_raw_medium_specific
        'n_edges': {'min': int(n_edges_arr.min()), 'median': float(np.median(n_edges_arr)),
                    'mean': float(n_edges_arr.mean()), 'max': int(n_edges_arr.max())},
        # FIX (Round 4 review, confirmed real): mean_fraction_* below is
        # mean(n_medium)/mean(n_edges) -- a POOLED ratio across all edges in
        # all topologies, which implicitly weights larger topologies more
        # heavily and is NOT the same quantity as "the average per-topology
        # medium fraction". Verified the two differ materially (fibre: 5.4%
        # pooled vs 12.8% mean-of-per-topology-ratios). Both are now
        # reported, clearly labelled, rather than only the pooled one under
        # an ambiguous "mean fraction of edges" label.
        'pooled_fraction_fibre': float(np.mean(n_fibre_list) / max(np.mean(n_edges_arr), 1e-9)),
        'pooled_fraction_satellite': float(np.mean(n_sat_list) / max(np.mean(n_edges_arr), 1e-9)),
        'pooled_fraction_fso': float(np.mean(n_fso_list) / max(np.mean(n_edges_arr), 1e-9)),
        'pooled_fraction_description': ('sum(edges of medium)/sum(all edges) across the whole '
                                        'Pareto front, equivalently mean(n_medium)/mean(n_edges) '
                                        '-- implicitly weights larger topologies more heavily'),
        'per_topology_fraction_mean_fibre': float(np.mean(np.array(n_fibre_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_std_fibre': float(np.std(np.array(n_fibre_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_mean_satellite': float(np.mean(np.array(n_sat_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_std_satellite': float(np.std(np.array(n_sat_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_mean_fso': float(np.mean(np.array(n_fso_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_std_fso': float(np.std(np.array(n_fso_list) / np.maximum(n_edges_arr, 1))),
        'per_topology_fraction_description': ('mean +/- std of each topology\'s OWN medium '
                                              'fraction, then averaged across topologies -- '
                                              'every topology weighted equally regardless of size'),
    },
}
with open(RESULTS / 'global_pareto_v3_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
print(f'\nSaved: global_pareto_v3.npz, global_pareto_v3_summary.json')
