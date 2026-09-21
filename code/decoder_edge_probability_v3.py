#!/usr/bin/env python3
"""
decoder_edge_probability_v3.py — effective per-pair edge-inclusion
probability of the topology decoder.

External review (finding: "Algorithm 1 undirected-edge sampling is
ambiguous") pointed out that the decoder loops over every node u and
its top-k neighbours v, drawing one Bernoulli(p) per (u, v) visit,
while storing edges in a set keyed on the canonical unordered triple
(min(u,v), max(u,v), medium). A pair that appears in BOTH endpoints'
k-nearest-neighbour lists is therefore visited twice and receives two
independent draws, so its effective inclusion probability is
1 - (1 - p)^2, not p. A pair that appears in only one endpoint's list
receives a single draw and keeps probability p.

This script quantifies both facts on the actual 100-city instance:
(a) what fraction of candidate pairs are mutual k-NN, and
(b) a Monte-Carlo check that the realised inclusion frequencies match
    1 - (1 - p)^2 and p respectively.

The behaviour is a property of the decoder as implemented and used for
every result in the paper; this artefact exists so that Algorithm 1's
description in the manuscript can state it exactly rather than leave it
to the reader to infer from the pseudocode.
"""
import json
from pathlib import Path

import numpy as np

import bdcz_oracle_v3 as O

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT_JSON = RESULTS / 'decoder_edge_probability_v3.json'

N = O.N_NODES
DIST = O.DIST
MC_TRIALS = 400
MC_RNG_BASE = 1000


def knn_sets(k, max_d):
    nb = {}
    for u in range(N):
        row = DIST[u]
        cand = sorted((row[v], v) for v in range(N) if v != u and row[v] <= max_d)[:k]
        nb[u] = {v for _, v in cand}
    return nb


def pair_census(k, max_d):
    nb = knn_sets(k, max_d)
    pairs = sorted({(min(u, v), max(u, v)) for u in range(N) for v in nb[u]})
    mutual = {(a, b) for (a, b) in pairs if b in nb[a] and a in nb[b]}
    return pairs, mutual


def monte_carlo(k, max_d, p, trials=MC_TRIALS):
    pairs, mutual = pair_census(k, max_d)
    counts = {pr: 0 for pr in pairs}
    for t in range(trials):
        params = np.array([k, 1, 0, max_d, 0.0, 0.0, p, 0.05, 0.05, 0.0, t], dtype=float)
        res = O.decode_topology(params, enable_satellite=False, enable_fso=False,
                                rng=np.random.default_rng(MC_RNG_BASE + t))
        if res is None:
            continue
        eu, ev, et = res
        for a, b, medium in zip(eu, ev, et):
            if medium != 0:
                continue
            key = (int(min(a, b)), int(max(a, b)))
            if key in counts:
                counts[key] += 1
    mut = [counts[pr] / trials for pr in pairs if pr in mutual]
    one = [counts[pr] / trials for pr in pairs if pr not in mutual]
    return pairs, mutual, mut, one


def main():
    census = {}
    for k, max_d in [(3, 3000.0), (5, 5000.0), (8, 8000.0)]:
        pairs, mutual = pair_census(k, max_d)
        census[f'k{k}_maxd{int(max_d)}'] = {
            'k_fibre': k,
            'max_d_fibre_km': max_d,
            'distinct_candidate_pairs': len(pairs),
            'mutual_knn_pairs': len(mutual),
            'mutual_fraction': len(mutual) / len(pairs),
        }

    mc = {}
    for p in [0.25, 0.5, 0.75]:
        pairs, mutual, mut, one = monte_carlo(3, 3000.0, p)
        mc[f'p{p}'] = {
            'p_nominal': p,
            'expected_mutual': 1 - (1 - p) ** 2,
            'observed_mutual_mean': float(np.mean(mut)),
            'n_mutual_pairs': len(mut),
            'expected_one_sided': p,
            'observed_one_sided_mean': float(np.mean(one)),
            'n_one_sided_pairs': len(one),
        }

    out = {
        'description': ('Effective per-pair edge-inclusion probability of the topology '
                        'decoder. Pairs that are mutual k-NN receive two independent '
                        'Bernoulli(p) draws (one per endpoint visit) and are stored in a '
                        'canonicalised set, so their effective inclusion probability is '
                        '1-(1-p)^2; pairs appearing in only one endpoint\'s list keep p. '
                        'Monte-Carlo run against the actual decode_topology() used for '
                        'every result in the paper.'),
        'monte_carlo_config': {'trials': MC_TRIALS, 'rng_base': MC_RNG_BASE,
                               'k_fibre': 3, 'max_d_fibre_km': 3000.0,
                               'satellite_and_fso_disabled': True},
        'candidate_pair_census': census,
        'monte_carlo': mc,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    for key, d in census.items():
        print(f"  {key}: {d['mutual_knn_pairs']}/{d['distinct_candidate_pairs']} pairs mutual "
              f"({100 * d['mutual_fraction']:.1f}%)")
    for key, d in mc.items():
        print(f"  {key}: mutual obs={d['observed_mutual_mean']:.3f} "
              f"(expected {d['expected_mutual']:.3f}); one-sided obs="
              f"{d['observed_one_sided_mean']:.3f} (expected {d['expected_one_sided']:.3f})")


if __name__ == '__main__':
    main()
