#!/usr/bin/env python3
"""
run_ablation_v3.py — Canonical-oracle ablation, 8 variants, 30 seeds
========================================================================
Adds a fibre/satellite/FSO-only completeness sweep the reviewer asked
for (FSO-only, variant I, was missing -- only fibre-anchored and
satellite-anchored combinations existed before). Seed count history:
5 -> 10 (to match the main comparison's then-10-seed power) -> 30
(external review round 6-followup: every significant p-value at n=10
sat exactly at the Wilcoxon resolution floor; raised to 30 to match
the main comparison, which had the same fix applied in Round 6).
Uses bdcz_oracle_v3 (canonical, Paper-1-unified physics).

FIX (Round 3 review, confirmed real bug): the variants below previously
tried to disable a medium by fixing its k/p decoder parameters to 0.0 via
`build_bounds`'s `fixed` dict. This did NOT work: decode_topology floors
k_fibre/k_sat at 1 and p_fibre/p_sat/p_fso at 0.05 regardless of what the
parameter vector says, so e.g. 'satellite-only' could still emit fibre
edges at p=0.05 per candidate. Each variant now ALSO passes explicit
enable_fibre/enable_satellite/enable_fso flags to eval_single, which
hard-disable a medium's entire edge-generation block in decode_topology --
this is what actually, unambiguously isolates a medium, independent of
the (now merely search-space-narrowing, not disabling) 'fixed' dict.

Variants
--------
A  full_decoder     — 11 free parameters (baseline), all 3 media enabled
B  fiber_only       — satellite, FSO hard-disabled
C  fiber_satellite  — FSO hard-disabled
D  fiber_fso        — satellite hard-disabled
E  no_interregion   — interregion_bonus=0, all 3 media enabled
G  satellite_only   — fibre, FSO hard-disabled
H  satellite_fso    — fibre hard-disabled
I  fso_only         — fibre, satellite hard-disabled           [NEW]
"""

import numpy as np
from metrics import compute_spacing
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

from bdcz_oracle_v3 import XL as XL_BASE, XU as XU_BASE, eval_single, N_VAR, dedup_front

from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize
from pymoo.operators.sampling.lhs import LHS
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD

POP_SIZE = 104
N_GEN    = 80
SEEDS    = list(range(42, 72))   # FIX (external review round 6-followup, MAJOR
                                  # finding): 10 seeds made every one-vs-A Wilcoxon
                                  # p-value that WAS significant sit exactly at the
                                  # n=10 resolution floor (p_raw=2/2^10=0.001953125),
                                  # with Cliff's delta saturated at +-1.00 -- the same
                                  # pathology the manuscript itself explains (and
                                  # resolves) for the n=30 main comparison, but had
                                  # left uncorrected here. Raised to 30 seeds to match.
HV_REF   = np.array([1.1] * 4)

OUT_DIR  = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = OUT_DIR / 'step3_ablation_v3'
OUT_JSON = OUT_DIR / 'ablation_v3_results.json'
OUT_NPZ  = OUT_DIR / 'ablation_v3_pareto.npz'

VARIANTS = {
    'A_full_decoder': {'description': '11 free parameters (baseline)', 'fixed': {},
                       'media': (True, True, True)},
    'B_fiber_only': {'description': 'satellite and FSO hard-disabled',
                     'fixed': {1: 0.0, 2: 0.0, 7: 0.0, 8: 0.0},
                     'media': (True, False, False)},
    'C_fiber_satellite': {'description': 'FSO hard-disabled',
                          'fixed': {2: 0.0, 8: 0.0},
                          'media': (True, True, False)},
    'D_fiber_fso': {'description': 'satellite hard-disabled',
                    'fixed': {1: 0.0, 7: 0.0},
                    'media': (True, False, True)},
    'E_no_interregion': {'description': 'interregion_bonus=0',
                        'fixed': {9: 0.0},
                        'media': (True, True, True)},
    'G_satellite_only': {'description': 'fibre and FSO hard-disabled',
                         'fixed': {0: 0.0, 2: 0.0, 6: 0.0, 8: 0.0},
                         'media': (False, True, False)},
    'H_satellite_fso': {'description': 'fibre hard-disabled',
                        'fixed': {0: 0.0, 6: 0.0},
                        'media': (False, True, True)},
    'I_fso_only': {'description': 'fibre and satellite hard-disabled',
                  'fixed': {0: 0.0, 1: 0.0, 6: 0.0, 7: 0.0},
                  'media': (False, False, True)},
}


def build_bounds(fixed_params):
    xl, xu = XL_BASE.copy(), XU_BASE.copy()
    for idx, val in fixed_params.items():
        xl[idx] = val
        xu[idx] = val
    return xl, xu


class QuantumTopoVariantProblem(Problem):
    def __init__(self, xl, xu, media=(True, True, True)):
        super().__init__(n_var=N_VAR, n_obj=4, n_constr=0, xl=xl, xu=xu)
        self.enable_fibre, self.enable_satellite, self.enable_fso = media

    def _evaluate(self, X, out, *args, **kwargs):
        results = np.zeros((len(X), 4))
        for i, params in enumerate(X):
            results[i] = eval_single(params, enable_fibre=self.enable_fibre,
                                      enable_satellite=self.enable_satellite,
                                      enable_fso=self.enable_fso)
        out['F'] = np.column_stack([
            -results[:, 0], -results[:, 1], results[:, 2], -results[:, 3],
        ])




def ckpt_path(v, s): return CKPT_DIR / f'{v}_seed{s}.pkl'
def done_path(v, s): return CKPT_DIR / f'{v}_seed{s}_completed.txt'
def is_done(v, s): return done_path(v, s).exists()

def save_run(v, s, X, F, elapsed_min, description):
    payload = {'X_pareto': X, 'F_pareto': F, 'n_pareto': len(F),
               'elapsed_min': elapsed_min, 'seed': s, 'algorithm': 'nsga3',
               'variant': v, 'description': description}
    with open(ckpt_path(v, s), 'wb') as fh:
        pickle.dump(payload, fh)
    with open(done_path(v, s), 'w') as fh:
        fh.write(f'Completed at {datetime.now().isoformat()}\n')

def load_run(v, s):
    with open(ckpt_path(v, s), 'rb') as fh:
        return pickle.load(fh)


def run_worker(variants, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    ref_dirs = get_reference_directions('das-dennis', 4, n_partitions=6)
    print(f'Ablation v3: {len(variants)} variants x {len(seeds)} seeds (this shard), canonical oracle')

    for variant_name in variants:
        vcfg = VARIANTS[variant_name]
        xl, xu = build_bounds(vcfg['fixed'])
        for seed in seeds:
            if is_done(variant_name, seed):
                print(f'  [SKIP] {variant_name} seed={seed}')
                continue
            print(f'  [RUN]  {variant_name} seed={seed} ...', end=' ', flush=True)
            t0 = time.time()

            problem = QuantumTopoVariantProblem(xl, xu, media=vcfg['media'])
            algo = NSGA3(pop_size=POP_SIZE, ref_dirs=ref_dirs, sampling=LHS(),
                        crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
            res = minimize(problem, algo,
                            termination=MaximumGenerationTermination(N_GEN),
                            seed=seed, verbose=False, save_history=False)
            elapsed_min = (time.time() - t0) / 60.0

            if res.F is not None and len(res.F) > 0:
                pop_F, pop_X = res.pop.get('F'), res.pop.get('X')
                nd_mask = find_non_dominated(pop_F)
                X_par, F_par = pop_X[nd_mask], pop_F[nd_mask]
                F_par, X_par = dedup_front(F_par, X_par)
            else:
                X_par, F_par = np.empty((0, N_VAR)), np.empty((0, 4))

            save_run(variant_name, seed, X_par, F_par, elapsed_min, vcfg['description'])
            print(f'n_pareto={len(F_par)}, {elapsed_min:.1f} min')


def aggregate():
    missing = [(v, s) for v in VARIANTS for s in SEEDS if not is_done(v, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} (variant, seed) combos missing, e.g. {missing[:5]}')
        return

    variant_fronts, variant_X_dict, variant_elapsed = {}, {}, {}
    for variant_name in VARIANTS:
        variant_fronts[variant_name], variant_X_dict[variant_name] = [], []
        variant_elapsed[variant_name] = []
        for seed in SEEDS:
            run = load_run(variant_name, seed)
            variant_fronts[variant_name].append(run['F_pareto'])
            variant_X_dict[variant_name].append(run['X_pareto'])
            variant_elapsed[variant_name].append(run['elapsed_min'])

    all_F_stack = np.vstack([F for fs in variant_fronts.values() for F in fs if len(F) > 0])
    g_min, g_max = all_F_stack.min(axis=0), all_F_stack.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    all_A_fronts = [f for f in variant_fronts['A_full_decoder'] if len(f) > 0]
    all_A_X = [x for x in variant_X_dict['A_full_decoder'] if len(x) > 0]
    if all_A_fronts:
        merged_A, merged_A_X = np.vstack(all_A_fronts), np.vstack(all_A_X)
        nd_A = find_non_dominated(merged_A)
        ref_raw, _ = dedup_front(merged_A[nd_A], merged_A_X[nd_A])
        ref_front = normalise(ref_raw)
    else:
        nd_all = find_non_dominated(all_F_stack)
        # all_X_stack mirrors all_F_stack's construction order for dedup's X arg
        all_X_stack = np.vstack([X for xs in variant_X_dict.values() for X in xs if len(X) > 0])
        ref_raw, _ = dedup_front(all_F_stack[nd_all], all_X_stack[nd_all])
        ref_front = normalise(ref_raw)

    # FIX (external review round 6, item 15): `ref_front` above (IGD vs.
    # A_full_decoder only) structurally favours variant A -- its own points
    # are trivially at zero distance from a reference built from itself,
    # which is not a symmetric comparison across variants. Added a SECOND,
    # symmetric reference front built from the deduplicated union of ALL
    # eight variants' non-dominated points, so every variant (including A)
    # is scored against the same shared target. Both are reported:
    # `igd` (vs. A_full_decoder only, kept for continuity with prior
    # rounds -- a meaningful "distance-from-the-best-known-variant" reading)
    # and `igd_union` (vs. the shared 8-variant front, the primary symmetric
    # comparison metric going forward).
    all_X_stack_u = np.vstack([X for xs in variant_X_dict.values() for X in xs if len(X) > 0])
    nd_union = find_non_dominated(all_F_stack)
    ref_raw_union, _ = dedup_front(all_F_stack[nd_union], all_X_stack_u[nd_union])
    ref_front_union = normalise(ref_raw_union)

    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)
    igd_calc_union = IGD(ref_front_union)

    summary, per_run_details, combined_fronts, combined_X_dict = {}, {}, {}, {}
    for variant_name in VARIANTS:
        per_run_details[variant_name] = []
        hv_list, igd_list, sp_list, np_list, el_list = [], [], [], [], []

        igd_union_list = []
        for idx, seed in enumerate(SEEDS):
            F_raw = variant_fronts[variant_name][idx]
            elapsed_min = variant_elapsed[variant_name][idx]
            if len(F_raw) == 0:
                per_run_details[variant_name].append(dict(seed=seed, n_pareto=0,
                    hv=0.0, igd=None, igd_union=None, spacing=0.0, elapsed_min=elapsed_min))
                continue
            F_norm = normalise(F_raw)
            hv_val, igd_val = float(hv_calc(F_norm)), float(igd_calc(F_norm))
            igd_union_val = float(igd_calc_union(F_norm))
            sp_val = compute_spacing(F_norm)
            per_run_details[variant_name].append(dict(seed=seed, n_pareto=len(F_raw),
                hv=hv_val, igd=igd_val, igd_union=igd_union_val, spacing=sp_val, elapsed_min=elapsed_min))
            hv_list.append(hv_val); igd_list.append(igd_val); sp_list.append(sp_val)
            igd_union_list.append(igd_union_val)
            np_list.append(len(F_raw)); el_list.append(elapsed_min)

        valid_F = [f for f in variant_fronts[variant_name] if len(f) > 0]
        valid_X = [x for x in variant_X_dict[variant_name] if len(x) > 0]
        if valid_F:
            merged_F, merged_X = np.vstack(valid_F), np.vstack(valid_X)
            nd_mask = find_non_dominated(merged_F)
            F_dedup, X_dedup = dedup_front(merged_F[nd_mask], merged_X[nd_mask])
            combined_fronts[variant_name] = normalise(F_dedup)
            combined_X_dict[variant_name] = X_dedup
        else:
            combined_fronts[variant_name] = np.empty((0, 4))
            combined_X_dict[variant_name] = np.empty((0, N_VAR))

        def ms(lst):
            a = np.array([x for x in lst if x is not None and np.isfinite(x)])
            if len(a) == 0:
                return {'mean': None, 'std': None}
            # ddof=1 (external review round 6, item 10): sample, not
            # population, standard deviation for a sample of stochastic runs.
            std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
            return {'mean': float(np.mean(a)), 'std': std}

        summary[variant_name] = {
            'description': VARIANTS[variant_name]['description'],
            'hv': ms(hv_list), 'igd': ms(igd_list), 'igd_union': ms(igd_union_list),
            'spacing': ms(sp_list),
            'n_pareto': ms(np_list), 'elapsed_min': ms(el_list),
            'hv_raw': hv_list,
        }

    npz_dict = {'g_min': g_min, 'g_max': g_max, 'g_range': g_range}
    for variant_name in VARIANTS:
        npz_dict[f'F_{variant_name}'] = combined_fronts[variant_name]
        npz_dict[f'X_{variant_name}'] = combined_X_dict[variant_name]
    np.savez(str(OUT_NPZ), **npz_dict)

    output = {
        'config': {'algorithm': 'nsga3', 'pop_size': POP_SIZE, 'n_gen': N_GEN,
                   'seeds': SEEDS, 'n_var': N_VAR, 'variants': list(VARIANTS.keys())},
        'normalisation': {'g_min': g_min.tolist(), 'g_max': g_max.tolist(),
                          'g_range': g_range.tolist()},
        'summary': summary,
        'per_run_details': per_run_details,
    }
    with open(OUT_JSON, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    print(f'Saved: {OUT_NPZ}')
    for variant_name in VARIANTS:
        m = summary[variant_name]
        print(f"  {variant_name:<20} HV={m['hv']['mean']:.4f}+/-{m['hv']['std']:.4f} "
              f"n_par={m['n_pareto']['mean']:.1f}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--variants', type=str, default=None,
                    help='comma-separated subset of VARIANTS keys, default all')
    p.add_argument('--seeds', type=str, default=None,
                    help='comma-separated seeds or start:end range, default all SEEDS')
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()

    variants = args.variants.split(',') if args.variants else list(VARIANTS.keys())
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS

    if not args.aggregate_only:
        run_worker(variants, seeds)
    else:
        aggregate()


if __name__ == '__main__':
    main()
