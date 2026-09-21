#!/usr/bin/env python3
"""
sigma_ablation_stats_v3.py — Kruskal-Wallis omnibus + Holm-corrected
pairwise Mann-Whitney + Cliff's delta for the three sigma-ablation
variants (A_sigma_optimised, B_sigma_fixed, C_sigma_averaged), computed
from the same per-seed HV/IGD values `run_sigma_ablation_v3.py`'s
aggregate() derives from its checkpoints (same normalisation/reference
front), so the reported p-values match the checkpoints exactly rather
than being recomputed under a different convention.

Run after `run_sigma_ablation_v3.py --aggregate-only` has produced
sigma_ablation_v3_results.json (this script re-reads the raw .pkl
checkpoints directly, not that JSON, since the JSON only stores
mean/std, not the per-seed values needed for Mann-Whitney).
"""
import json
import pickle
from itertools import combinations
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from scipy.stats import kruskal, mannwhitneyu

from run_sigma_ablation_v3 import VARIANTS, SEEDS, HV_REF, CKPT_DIR, ckpt_path

OUT_JSON = Path(__file__).resolve().parent.parent / 'results' / 'sigma_ablation_v3_stats.json'


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running_max = 0.0
    for rank, idx in enumerate(order):
        corrected = (m - rank) * pvals[idx]
        running_max = max(running_max, corrected)
        adj[idx] = min(running_max, 1.0)
    return adj


def main():
    variant_fronts = {}
    for v in VARIANTS:
        variant_fronts[v] = []
        for s in SEEDS:
            with open(ckpt_path(v, s), 'rb') as fh:
                variant_fronts[v].append(pickle.load(fh))

    all_F = np.vstack([r['F_pareto'] for v in VARIANTS for r in variant_fronts[v] if len(r['F_pareto']) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_all = find_non_dominated(all_F)
    ref_front = normalise(np.unique(np.round(all_F[nd_all], 6), axis=0))
    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

    per_seed_hv = {}
    for v in VARIANTS:
        hv_l = []
        for r in variant_fronts[v]:
            F = r['F_pareto']
            hv_l.append(float(hv_calc(normalise(F))) if len(F) > 0 else 0.0)
        per_seed_hv[v] = hv_l

    variants = list(VARIANTS.keys())
    kw_stat, kw_p = kruskal(*[per_seed_hv[v] for v in variants])

    pairs = list(combinations(variants, 2))
    raw_p = []
    mw_results = []
    for a, b in pairs:
        u_stat, p = mannwhitneyu(per_seed_hv[a], per_seed_hv[b], alternative='two-sided')
        delta = cliffs_delta(per_seed_hv[a], per_seed_hv[b])
        raw_p.append(p)
        mw_results.append({'pair': [a, b], 'u_stat': float(u_stat), 'p_raw': float(p),
                            'cliffs_delta': float(delta)})
    adj_p = holm(np.array(raw_p))
    for i, r in enumerate(mw_results):
        r['p_holm'] = float(adj_p[i])

    out = {
        'kruskal_wallis': {'H': float(kw_stat), 'p': float(kw_p), 'metric': 'HV'},
        'pairwise_mannwhitney_holm': mw_results,
        'n_seeds': len(SEEDS),
        'per_seed_hv': per_seed_hv,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    print(f"Kruskal-Wallis H={kw_stat:.3f} p={kw_p:.4g}")
    for r in mw_results:
        print(f"  {r['pair'][0]} vs {r['pair'][1]}: p_holm={r['p_holm']:.4g} delta={r['cliffs_delta']:+.2f}")


if __name__ == '__main__':
    main()
