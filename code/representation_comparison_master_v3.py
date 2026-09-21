#!/usr/bin/env python3
"""
representation_comparison_master_v3.py — Round 8: a single master
comparison of every representation/operator combination tested in this
paper for the 4,950-candidate-edge topology-design problem, all scored
against ONE shared, deduplicated reference front and ONE shared
min-max normalisation, so their HV/IGD numbers are directly comparable
in one table -- unlike the paper's other categorical-baseline
comparisons, each of which uses its own block-specific reference front
(a source of real reviewer confusion flagged across multiple rounds:
the same decoder+NSGA-II run's HV reads differently in each block).

Representations compared (all NSGA-II, same evolutionary budget):
  - decoder            : the 11-variable parametric decoder (30 seeds)
  - relaxed_categorical : 4,950-variable real-valued relaxation,
                          SBX/PM operators (30 seeds)
  - category_native     : 4,950-variable genuinely discrete
                          representation, category-native operators
                          (uniform crossover, add/remove/change-medium
                          mutation), best density (d03); 30 seeds.

  FIX (Round 13, post satellite-visibility correction): the best
  category-native density used to be d15 (15% initial edge density),
  and was extended to 30 seeds on that basis. Under the corrected
  satellite model the density-vs-HV relationship inverted: d03 (3%
  initial density) now dominates d10/d11/d15 by roughly an order of
  magnitude (denser initial populations waste more of their edges on
  now-infeasible long-range satellite links than the category-native
  mutation operators can productively fix within budget). d03 has
  therefore been extended to 30 seeds instead, and is now the primary
  row; d10, d11_decoder_median and d15 are the secondary, 10-seed
  exploratory sweep.
"""
import hashlib
import json
import pickle
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from scipy.stats import kruskal, mannwhitneyu

from metrics import compute_spacing

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT_JSON = RESULTS / 'representation_comparison_master_v3.json'
HV_REF = np.array([1.1] * 4)

DECODER_SEEDS = list(range(42, 72))       # 30 seeds
CATNATIVE_PRIMARY_DENSITY = 'd03'          # best-performing density, Round 13
CATNATIVE_PRIMARY_SEEDS = list(range(42, 72))  # 30 seeds after Round 13 extension
CATNATIVE_OTHER_SEEDS = list(range(42, 52))  # 10 seeds, exploratory


def load_decoder(seed):
    with open(RESULTS / 'step1_baselines_v3' / f'nsga2_seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)['F_pareto']


def load_relaxed_categorical(seed):
    with open(RESULTS / 'step2_binary_v3' / f'nsga2_categorical_seed{seed}.pkl') as f:
        entry = json.load(f)
    pf = entry['pareto_front']
    return np.array(pf).reshape(-1, 4) if pf else np.empty((0, 4))


def load_category_native(density, seed):
    with open(RESULTS / 'step5_catops_v3' / f'{density}_seed{seed}.pkl') as f:
        entry = json.load(f)
    pf = entry['pareto_front']
    return np.array(pf).reshape(-1, 4) if pf else np.empty((0, 4))


def sha256_of_array(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).round(6).tobytes()).hexdigest()[:16]


def main():
    decoder_fronts = [load_decoder(s) for s in DECODER_SEEDS]
    relaxed_fronts = [load_relaxed_categorical(s) for s in DECODER_SEEDS]
    primary_native_fronts = [load_category_native(CATNATIVE_PRIMARY_DENSITY, s)
                             for s in CATNATIVE_PRIMARY_SEEDS]
    other_density_fronts = {
        d: [load_category_native(d, s) for s in CATNATIVE_OTHER_SEEDS]
        for d in ['d10', 'd11_decoder_median', 'd15']
    }

    all_F = np.vstack(
        [F for F in decoder_fronts if len(F) > 0]
        + [F for F in relaxed_fronts if len(F) > 0]
        + [F for F in primary_native_fronts if len(F) > 0]
    )
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_all = find_non_dominated(all_F)
    ref_front = normalise(np.unique(np.round(all_F[nd_all], 6), axis=0))
    hv_calc = Hypervolume(ref_point=HV_REF)

    def summarise(fronts):
        hv_l, spacing_l, npar_l = [], [], []
        for F in fronts:
            if len(F) == 0:
                hv_l.append(0.0); spacing_l.append(0.0); npar_l.append(0)
                continue
            Fn = normalise(F)
            hv_l.append(float(hv_calc(Fn)))
            spacing_l.append(compute_spacing(Fn))
            npar_l.append(len(F))
        arr = lambda l: np.array(l, dtype=float)
        def ms(a):
            std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
            return {'mean': float(np.mean(a)), 'std': std}
        return {'hv': ms(arr(hv_l)), 'spacing': ms(arr(spacing_l)),
                'n_pareto': ms(arr(npar_l)), 'hv_raw': hv_l, 'n_seeds': len(fronts)}

    primary = {
        'decoder': summarise(decoder_fronts),
        'relaxed_categorical_sbx_pm': summarise(relaxed_fronts),
        'category_native_d03': summarise(primary_native_fronts),
    }

    kw_stat, kw_p = kruskal(primary['decoder']['hv_raw'],
                             primary['relaxed_categorical_sbx_pm']['hv_raw'],
                             primary['category_native_d03']['hv_raw'])

    pairwise = {}
    names = list(primary.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            u, p = mannwhitneyu(primary[a]['hv_raw'], primary[b]['hv_raw'], alternative='two-sided')
            pairwise[f'{a}_vs_{b}'] = {'u_stat': float(u), 'p': float(p)}

    secondary = {d: summarise(fronts) for d, fronts in other_density_fronts.items()}

    out = {
        'description': ('Single master comparison: decoder, relaxed-categorical SBX/PM, '
                         'and category-native (best density d03, Round 13 -- was d15 '
                         'before the satellite-visibility oracle correction inverted the '
                         'density-vs-HV relationship) all scored against ONE '
                         'shared, deduplicated reference front and normalisation. Other '
                         'category-native densities (d10, d11_decoder_median, d15) '
                         'reported separately as a secondary, '
                         '10-seed exploratory sweep (not matched to the primary rows\' '
                         '30-seed statistical power).'),
        'reference_front_hash': sha256_of_array(ref_front),
        'reference_front_n_points': len(ref_front),
        'normalisation_bounds_hash': sha256_of_array(np.concatenate([g_min, g_range])),
        'reference_point': HV_REF.tolist(),
        'deduplicated': True,
        'feasible_only': True,
        'primary_30_seed_comparison': primary,
        'kruskal_wallis_hv': {'H': float(kw_stat), 'p': float(kw_p)},
        'pairwise_mannwhitney_hv': pairwise,
        'secondary_10_seed_exploratory_densities': secondary,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    for name, s in primary.items():
        print(f"  {name:28s} HV={s['hv']['mean']:.4f}+/-{s['hv']['std']:.4f} "
              f"spacing={s['spacing']['mean']:.4f} n_par={s['n_pareto']['mean']:.1f} (n={s['n_seeds']})")
    print(f"Kruskal-Wallis H={kw_stat:.3f} p={kw_p:.4g}")
    for k, v in pairwise.items():
        print(f"  {k}: p={v['p']:.4g}")


if __name__ == '__main__':
    main()
