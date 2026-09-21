#!/usr/bin/env python3
"""
categorical_operators_shared_ref_v3.py — shared-reference HV comparison
between the category-native-operator baseline (categorical_operators_baseline_v3.py)
and the decoder's own NSGA-II runs, using ONE common reference front/
normalisation so the two HV numbers are directly comparable (their own,
separately-normalised HV values in categorical_operators_v3_results.json
and baselines_comparison_v3.json are NOT comparable to each other, by
construction of that per-block normalisation -- this script builds the
one number that IS).

Uses the same 10 seeds (42-51) for both sides, reading directly from
each side's cached checkpoints. Run after both
`categorical_operators_baseline_v3.py --aggregate-only` and
`run_baselines_comparison_v3.py --aggregate-only` have completed.
"""
import json
import pickle
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from pymoo.util.nds.non_dominated_sorting import find_non_dominated

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CATOPS_CKPT = RESULTS / 'step5_catops_v3'
BASELINE_CKPT = RESULTS / 'step1_baselines_v3'
OUT_JSON = RESULTS / 'categorical_operators_vs_decoder_shared_ref_v3.json'

SEEDS = list(range(42, 52))
DENSITIES = ['d03', 'd10', 'd11_decoder_median', 'd15']
HV_REF = np.array([1.1] * 4)


def load_catops(density, seed):
    with open(CATOPS_CKPT / f'{density}_seed{seed}.pkl') as f:
        entry = json.load(f)
    return {'F_pareto': np.array(entry['pareto_front']).reshape(-1, 4) if entry['pareto_front'] else np.empty((0, 4))}


def load_decoder_nsga2(seed):
    with open(BASELINE_CKPT / f'nsga2_seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)


def main():
    catops_runs = {d: [load_catops(d, s) for s in SEEDS] for d in DENSITIES}
    decoder_runs = [load_decoder_nsga2(s) for s in SEEDS]

    all_F = np.vstack(
        [r['F_pareto'] for d in DENSITIES for r in catops_runs[d] if len(r['F_pareto']) > 0]
        + [r['F_pareto'] for r in decoder_runs if len(r['F_pareto']) > 0]
    )
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    hv_calc = Hypervolume(ref_point=HV_REF)

    def summarise(runs):
        hv_l = [float(hv_calc(normalise(r['F_pareto']))) for r in runs if len(r['F_pareto']) > 0]
        arr = np.array(hv_l)
        std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        return {'mean': float(np.mean(arr)), 'std': std, 'n_seeds': len(arr)}

    out = {
        'description': ('HV of the category-native-operator baseline (4 densities) and the '
                         'decoder\'s own NSGA-II, both computed under ONE shared reference '
                         'front/normalisation (unlike their block-specific own numbers '
                         'elsewhere), same 10 seeds (42-51) both sides.'),
        'decoder_nsga2': summarise(decoder_runs),
        'categorical_operators': {d: summarise(catops_runs[d]) for d in DENSITIES},
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    print(f"decoder_nsga2 HV={out['decoder_nsga2']['mean']:.4f}+/-{out['decoder_nsga2']['std']:.4f}")
    for d in DENSITIES:
        m = out['categorical_operators'][d]
        print(f"  {d:<22} HV={m['mean']:.4f}+/-{m['std']:.4f}")


if __name__ == '__main__':
    main()
