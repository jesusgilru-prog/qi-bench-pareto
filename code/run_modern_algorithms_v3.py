#!/usr/bin/env python3
"""
run_modern_algorithms_v3.py — Round 8: an external review correctly
noted the main comparison's newest algorithm (RVEA, 2016) is dated for
a 2026 submission. Adds two post-2020 many-objective algorithms under
the EXACT same conditions as the main comparison (same 11-variable
decoder, same oracle, same bounds, same SBX/PM operators, same
LHS sampling, same 8,320-evaluation budget, 30 seeds):

  - SPEA2 (Zitzler, Laumanns & Thiele 2001) -- included as a
    long-standing, non-dominance-and-density baseline distinct from
    NSGA-II's crowding distance; genuinely predates 2020, included
    for completeness of the "classical Pareto-based" family rather
    than as the "modern" entry.
  - AGE-MOEA-II (Panichella 2022) -- a genuinely 2020s many-objective
    algorithm using adaptive geometry estimation instead of a fixed
    reference-direction set (unlike NSGA-III/RVEA/MOEA-D, all of
    which need ref_dirs), included specifically to test whether a
    contemporary many-objective method fares any better than
    NSGA-III/RVEA on this benchmark's low effective dimensionality
    (Section: Discussion). Requires `numba` (optional dependency,
    listed in requirements.txt as such).

Reuses QuantumTopoProblem, make_algorithm-equivalent construction and
n_gen_for's matched-budget logic directly from
run_baselines_comparison_v3.py rather than duplicating the oracle
wrapper.
"""
import argparse
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.algorithms.moo.age2 import AGEMOEA2
from pymoo.algorithms.moo.spea2 import SPEA2
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from scipy.stats import mannwhitneyu

from run_baselines_comparison_v3 import QuantumTopoProblem, POP_SIZE, N_GEN
from bdcz_oracle_v3 import dedup_front
from metrics import compute_spacing

SEEDS = list(range(42, 72))  # 30 seeds, matching the main comparison
HV_REF = np.array([1.1] * 4)
ALGORITHMS = ['spea2', 'age_moea2']

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step9_modern_algorithms_v3'
OUT_JSON = RESULTS / 'modern_algorithms_v3_results.json'


def make_algorithm(name):
    if name == 'spea2':
        return SPEA2(pop_size=POP_SIZE, sampling=LHS(),
                     crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'age_moea2':
        return AGEMOEA2(pop_size=POP_SIZE, sampling=LHS(),
                        crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    raise ValueError(name)


def ckpt_path(algo, seed): return CKPT_DIR / f'{algo}_seed{seed}.pkl'
def done_path(algo, seed): return CKPT_DIR / f'{algo}_seed{seed}_completed.txt'
def is_done(algo, seed): return done_path(algo, seed).exists()


def run_worker(algos, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    for algo in algos:
        for seed in seeds:
            if is_done(algo, seed):
                print(f'  [SKIP] {algo} seed={seed}')
                continue
            print(f'  [RUN]  {algo} seed={seed} ...', end=' ', flush=True)
            t0 = time.time()
            problem = QuantumTopoProblem()
            algorithm = make_algorithm(algo)
            res = minimize(problem, algorithm, termination=MaximumGenerationTermination(N_GEN),
                            seed=seed, verbose=False, save_history=False)
            elapsed_min = (time.time() - t0) / 60.0
            n_calls = problem.n_calls if hasattr(problem, 'n_calls') else POP_SIZE * N_GEN
            if res.F is not None and len(res.F) > 0:
                pop_F, pop_X = res.pop.get('F'), res.pop.get('X')
                nd_mask = find_non_dominated(pop_F)
                X_par, F_par = pop_X[nd_mask], pop_F[nd_mask]
                F_par, X_par = dedup_front(F_par, X_par)
            else:
                X_par, F_par = np.empty((0, 11)), np.empty((0, 4))
            payload = {'X_pareto': X_par, 'F_pareto': F_par, 'n_pareto': len(F_par),
                       'elapsed_min': elapsed_min, 'seed': seed, 'algo': algo, 'n_calls': n_calls}
            with open(ckpt_path(algo, seed), 'wb') as fh:
                pickle.dump(payload, fh)
            with open(done_path(algo, seed), 'w') as fh:
                fh.write(f'Completed at {datetime.now().isoformat()}\n')
            print(f'n_pareto={len(F_par)}, n_calls={n_calls}, {elapsed_min:.1f} min')


def aggregate():
    missing = [(a, s) for a in ALGORITHMS for s in SEEDS if not is_done(a, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} combos missing, e.g. {missing[:5]}')
        return

    runs = {a: [] for a in ALGORITHMS}
    for a in ALGORITHMS:
        for s in SEEDS:
            with open(ckpt_path(a, s), 'rb') as fh:
                runs[a].append(pickle.load(fh))

    # Pool together with the main comparison's own RAW per-seed fronts
    # (not the already-globally-filtered 659-point front) so the
    # normalisation bounds match Table 3's exact convention: Table 3's
    # own g_min/g_max come from the union of all 150 raw per-seed
    # fronts, which is WIDER than the globally-non-dominated 659-point
    # front (some individually-non-dominated points get filtered out of
    # the global front by other runs' points, but still legitimately
    # widen the normalisation bounding box) -- using the narrower
    # front here would make HV numbers look comparable to Table 3 while
    # actually sitting on a different scale.
    MAIN_ALGOS = ['nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea']
    MAIN_SEEDS = list(range(42, 72))
    main_raw_fronts = []
    for a in MAIN_ALGOS:
        for s in MAIN_SEEDS:
            with open(RESULTS / 'step1_baselines_v3' / f'{a}_seed{s}.pkl', 'rb') as f:
                main_raw_fronts.append(pickle.load(f)['F_pareto'])

    all_F = np.vstack([F for F in main_raw_fronts if len(F) > 0]
                       + [r['F_pareto'] for a in ALGORITHMS for r in runs[a] if len(r['F_pareto']) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_all = find_non_dominated(all_F)
    ref_front = normalise(np.unique(np.round(all_F[nd_all], 6), axis=0))
    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

    summary = {}
    for a in ALGORITHMS:
        hv_l, igd_l, spacing_l, npar_l, elapsed_l = [], [], [], [], []
        for r in runs[a]:
            F = r['F_pareto']
            if len(F) == 0:
                continue
            Fn = normalise(F)
            hv_l.append(float(hv_calc(Fn)))
            igd_l.append(float(igd_calc(Fn)))
            spacing_l.append(compute_spacing(Fn))
            npar_l.append(len(F))
            elapsed_l.append(r['elapsed_min'])

        def ms(lst):
            arr = np.array(lst)
            std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
            return {'mean': float(np.mean(arr)), 'std': std}

        summary[a] = {'hv': ms(hv_l), 'igd': ms(igd_l), 'spacing': ms(spacing_l),
                      'n_pareto': ms(npar_l), 'elapsed_min': ms(elapsed_l)}
        summary[a]['hv_raw'] = hv_l

    # Persist significance vs. NSGA-II/NSGA-III under this SAME reference,
    # so the manuscript's p-values are traceable to this JSON directly
    # rather than only reproducible by re-running an ad hoc script.
    main_hv_raw = {}
    for a in ['nsga2', 'nsga3']:
        hv_l = []
        for s in MAIN_SEEDS:
            with open(RESULTS / 'step1_baselines_v3' / f'{a}_seed{s}.pkl', 'rb') as f:
                F = pickle.load(f)['F_pareto']
            hv_l.append(float(hv_calc(normalise(F))) if len(F) > 0 else 0.0)
        main_hv_raw[a] = hv_l

    significance = {}
    for a in ALGORITHMS:
        for main_a in ['nsga2', 'nsga3']:
            u, p = mannwhitneyu(summary[a]['hv_raw'], main_hv_raw[main_a], alternative='two-sided')
            significance[f'{a}_vs_{main_a}'] = {'u_stat': float(u), 'p': float(p)}

    output = {
        'description': ('SPEA2 and AGE-MOEA-II under the exact same conditions as the '
                         'main 5-algorithm comparison (Table 3), re-scored against a '
                         'reference front that pools their own fronts WITH the main '
                         "comparison's own RAW per-seed fronts (step1_baselines_v3/*.pkl, "
                         'the same 150-front union Table 3 itself normalises against -- '
                         'NOT the already-globally-deduplicated global_pareto_v3.npz, '
                         'which has narrower bounds and would silently put HV on a '
                         'different scale despite looking comparable), so HV '
                         'here is directly comparable to Table 3, not a separately '
                         'normalised quantity (unlike this paper\'s other categorical-'
                         'baseline comparisons).'),
        'config': {'pop_size': POP_SIZE, 'n_gen': N_GEN, 'seeds': SEEDS, 'algorithms': ALGORITHMS},
        'summary': summary,
        'significance_vs_main': significance,
    }
    with open(OUT_JSON, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    for a in ALGORITHMS:
        m = summary[a]
        print(f"  {a:<12} HV={m['hv']['mean']:.4f}+/-{m['hv']['std']:.4f} "
              f"IGD={m['igd']['mean']:.4f}+/-{m['igd']['std']:.4f} n_par={m['n_pareto']['mean']:.1f}")
    for k, v in significance.items():
        print(f"  {k}: p={v['p']:.4g}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--algos', type=str, default=None)
    p.add_argument('--seeds', type=str, default=None)
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()
    algos = args.algos.split(',') if args.algos else ALGORITHMS
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS
    if not args.aggregate_only:
        run_worker(algos, seeds)
    missing = [(a, s) for a in ALGORITHMS for s in SEEDS if not is_done(a, s)]
    if missing:
        print(f'\n{len(missing)} (algo, seed) combos still missing, not aggregating yet.')
        return
    aggregate()


if __name__ == '__main__':
    main()
