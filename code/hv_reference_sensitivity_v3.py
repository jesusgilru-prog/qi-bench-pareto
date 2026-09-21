#!/usr/bin/env python3
"""
hv_reference_sensitivity_v3.py — HV reference-point sensitivity check
(external review round 6, ChatGPT item 13). Recomputes HV for every
algorithm x seed at three normalised reference points (1.05, 1.10, 1.20 --
the main text uses 1.10) WITHOUT re-running any optimisation, using the
already-computed and already-normalised per-run fronts. Checks whether the
algorithm ranking (by mean HV) is stable across reference-point choices.
"""
import json
import pickle
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

from bdcz_oracle_v3 import dedup_front
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.indicators.hv import Hypervolume

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step1_baselines_v3'

with open(RESULTS / 'baselines_comparison_v3.json') as f:
    base = json.load(f)
algos = base['config']['algorithms']
seeds = base['config']['seeds']
g_min = np.array(base['normalisation']['g_min'])
g_range = np.array(base['normalisation']['g_range'])


def normalise(F):
    return (F - g_min) / g_range


REF_POINTS = [1.05, 1.10, 1.20]
hv_by_ref = {rp: {a: [] for a in algos} for rp in REF_POINTS}

for algo in algos:
    for seed in seeds:
        p = CKPT_DIR / f'{algo}_seed{seed}.pkl'
        if not p.exists():
            continue
        with open(p, 'rb') as fh:
            run = pickle.load(fh)
        F_raw = run['F_pareto']
        if len(F_raw) == 0:
            continue
        F_norm = normalise(F_raw)
        for rp in REF_POINTS:
            hv_calc = Hypervolume(ref_point=np.array([rp] * 4))
            hv_by_ref[rp][algo].append(float(hv_calc(F_norm)))

means = {rp: {a: float(np.mean(v)) if v else None for a, v in hv_by_ref[rp].items()}
         for rp in REF_POINTS}
rankings = {rp: sorted(algos, key=lambda a: means[rp][a] or 0.0, reverse=True) for rp in REF_POINTS}

print('Mean HV by algorithm at each reference point:')
for rp in REF_POINTS:
    print(f'  ref={rp}: ' + ', '.join(f'{a}={means[rp][a]:.4f}' for a in algos))
print()
print('Ranking by mean HV at each reference point:')
for rp in REF_POINTS:
    print(f'  ref={rp}: {" > ".join(rankings[rp])}')

# Per-seed Spearman rank correlation between ref=1.10 (main text) and the
# two alternatives, to check the ranking is not an artefact of one choice.
base_rp = 1.10
spearman_out = {}
for rp in REF_POINTS:
    if rp == base_rp:
        continue
    corrs = []
    for seed_idx in range(len(seeds)):
        v1 = [hv_by_ref[base_rp][a][seed_idx] for a in algos if seed_idx < len(hv_by_ref[base_rp][a])]
        v2 = [hv_by_ref[rp][a][seed_idx] for a in algos if seed_idx < len(hv_by_ref[rp][a])]
        if len(v1) == len(algos) and len(v2) == len(algos):
            rho, _ = spearmanr(v1, v2)
            corrs.append(rho)
    spearman_out[rp] = {'mean_spearman_vs_ref1.10': float(np.mean(corrs)) if corrs else None,
                         'n_seeds_compared': len(corrs)}
    print(f'\nPer-seed Spearman rank correlation, ref={rp} vs. ref={base_rp}: '
          f'mean={spearman_out[rp]["mean_spearman_vs_ref1.10"]:.3f} (n={len(corrs)} seeds)')

out = {
    'reference_points': REF_POINTS,
    'mean_hv_by_ref': means,
    'ranking_by_ref': rankings,
    'ranking_stable': len(set(tuple(r) for r in rankings.values())) == 1,
    'spearman_vs_ref1.10': spearman_out,
}
out_path = RESULTS / 'hv_reference_sensitivity_v3.json'
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2)
print(f'\nSaved: {out_path}')
print(f'Ranking stable across all 3 reference points: {out["ranking_stable"]}')
