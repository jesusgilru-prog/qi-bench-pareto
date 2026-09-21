#!/usr/bin/env python3
"""
run_multi_instance_v3.py — external review round 7 (author-approved full
scope): re-run the 5-algorithm comparison on 3 additional, smaller/denser
REAL city-subset instances (see multi_instance_oracle_v3.py for exactly
which real cities and why), to test whether the algorithm ranking is an
artefact of the specific 100-city benchmark or holds more generally.

Reduced seed count (10 per instance, not 30) since this is a
generalisation check across instances, not a new headline statistical
comparison for the main benchmark.
"""
import argparse
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.algorithms.moo.moead import MOEAD
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.moo.rvea import RVEA
from pymoo.algorithms.moo.sms import SMSEMOA
from pymoo.core.problem import Problem
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.util.ref_dirs import get_reference_directions

from bdcz_oracle_v3 import N_VAR, dedup_front
from multi_instance_oracle_v3 import Instance, INSTANCES

POP_SIZE, N_GEN = 104, 80
SEEDS = list(range(42, 52))  # 10 seeds per instance
HV_REF = np.array([1.1] * 4)
ALGORITHMS = ['nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea']

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step7_multi_instance_v3'
OUT_JSON = RESULTS / 'multi_instance_v3_results.json'


class InstanceProblem(Problem):
    """Decoder search-space bounds (XL/XU) are decoder-parameter bounds,
    not city-count-dependent, so they are kept IDENTICAL to the main
    100-city benchmark for a fair, matched-search-space comparison
    across instances."""
    def __init__(self, instance):
        from bdcz_oracle_v3 import XL, XU
        super().__init__(n_var=N_VAR, n_obj=4, n_constr=0, xl=XL, xu=XU)
        self.instance = instance
        self.n_calls = 0

    def _evaluate(self, X, out, *args, **kwargs):
        F = np.zeros((len(X), 4))
        for i, params in enumerate(X):
            fid, rate_log, lat, cov = self.instance.eval_single(params)
            F[i] = [-fid, -rate_log, lat, -cov]
        self.n_calls += len(X)
        out['F'] = F


def make_algorithm(name, ref_dirs):
    if name == 'nsga2':
        return NSGA2(pop_size=POP_SIZE, sampling=LHS(), crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'nsga3':
        return NSGA3(pop_size=POP_SIZE, ref_dirs=ref_dirs, sampling=LHS(), crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'moead':
        return MOEAD(ref_dirs=ref_dirs, n_neighbors=15, sampling=LHS(), crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'sms_emoa':
        return SMSEMOA(pop_size=POP_SIZE, sampling=LHS(), crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'rvea':
        return RVEA(ref_dirs=ref_dirs, pop_size=POP_SIZE, sampling=LHS(), crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    raise ValueError(name)


def n_gen_for(algo, ref_dirs):
    if algo == 'moead':
        target_evals = POP_SIZE * (N_GEN + 1)
        return max(1, round(target_evals / len(ref_dirs)) - 1)
    return N_GEN


def ckpt_path(inst, algo, s): return CKPT_DIR / f'{inst}_{algo}_seed{s}.pkl'
def done_path(inst, algo, s): return CKPT_DIR / f'{inst}_{algo}_seed{s}_completed.txt'
def is_done(inst, algo, s): return done_path(inst, algo, s).exists()


def run_worker(instances, algos, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    ref_dirs = get_reference_directions('das-dennis', 4, n_partitions=6)
    for inst_name in instances:
        instance = Instance(inst_name)
        for algo in algos:
            n_gen_algo = n_gen_for(algo, ref_dirs)
            for seed in seeds:
                if is_done(inst_name, algo, seed):
                    print(f'  [SKIP] {inst_name} {algo} seed={seed}')
                    continue
                print(f'  [RUN]  {inst_name} (N={instance.n_nodes}) {algo} seed={seed} ...', end=' ', flush=True)
                t0 = time.time()
                problem = InstanceProblem(instance)
                algorithm = make_algorithm(algo, ref_dirs)
                res = minimize(problem, algorithm, termination=MaximumGenerationTermination(n_gen_algo),
                                seed=seed, verbose=False, save_history=False)
                elapsed_min = (time.time() - t0) / 60.0
                if res.F is not None and len(res.F) > 0:
                    pop_F, pop_X = res.pop.get('F'), res.pop.get('X')
                    nd_mask = find_non_dominated(pop_F)
                    X_par, F_par = pop_X[nd_mask], pop_F[nd_mask]
                    F_par, X_par = dedup_front(F_par, X_par)
                else:
                    X_par, F_par = np.empty((0, N_VAR)), np.empty((0, 4))
                payload = {'X_pareto': X_par, 'F_pareto': F_par, 'n_pareto': len(F_par),
                           'elapsed_min': elapsed_min, 'seed': seed, 'algorithm': algo,
                           'instance': inst_name, 'n_calls': problem.n_calls}
                with open(ckpt_path(inst_name, algo, seed), 'wb') as fh:
                    pickle.dump(payload, fh)
                with open(done_path(inst_name, algo, seed), 'w') as fh:
                    fh.write(f'Completed at {datetime.now().isoformat()}\n')
                print(f'n_pareto={len(F_par)}, n_calls={problem.n_calls}, {elapsed_min:.1f} min')


def aggregate():
    instances = list(INSTANCES.keys())
    missing = [(i, a, s) for i in instances for a in ALGORITHMS for s in SEEDS if not is_done(i, a, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} combos missing, e.g. {missing[:5]}')
        return

    out = {'config': {'instances': INSTANCES, 'algorithms': ALGORITHMS, 'seeds': SEEDS, 'pop_size': POP_SIZE},
           'per_instance': {}}
    for inst_name in instances:
        runs = {a: [] for a in ALGORITHMS}
        for a in ALGORITHMS:
            for s in SEEDS:
                with open(ckpt_path(inst_name, a, s), 'rb') as fh:
                    runs[a].append(pickle.load(fh))

        all_F = np.vstack([r['F_pareto'] for a in ALGORITHMS for r in runs[a] if len(r['F_pareto']) > 0])
        g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
        g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

        def normalise(F):
            return (F - g_min) / g_range

        nd = find_non_dominated(all_F)
        ref_front = normalise(np.unique(np.round(all_F[nd], 6), axis=0))
        hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

        summary = {}
        for a in ALGORITHMS:
            hv_l, igd_l, npar_l = [], [], []
            for r in runs[a]:
                F = r['F_pareto']
                if len(F) == 0:
                    continue
                Fn = normalise(F)
                hv_l.append(float(hv_calc(Fn)))
                igd_l.append(float(igd_calc(Fn)))
                npar_l.append(len(F))
            def ms(lst):
                arr = np.array(lst)
                std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                return {'mean': float(np.mean(arr)), 'std': std}
            summary[a] = {'hv': ms(hv_l), 'igd': ms(igd_l), 'n_pareto': ms(npar_l)}
        ranking = sorted(ALGORITHMS, key=lambda a: summary[a]['hv']['mean'], reverse=True)
        out['per_instance'][inst_name] = {
            'n_nodes': len(INSTANCES[inst_name]['indices']),
            'n_evaluations': runs[ALGORITHMS[0]][0].get('n_calls'),
            'summary': summary, 'hv_ranking': ranking}

    with open(OUT_JSON, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    for inst_name in instances:
        print(f"  {inst_name}: ranking = {' > '.join(out['per_instance'][inst_name]['hv_ranking'])}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--instances', type=str, default=None)
    p.add_argument('--algos', type=str, default=None)
    p.add_argument('--seeds', type=str, default=None)
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()
    instances = args.instances.split(',') if args.instances else list(INSTANCES.keys())
    algos = args.algos.split(',') if args.algos else ALGORITHMS
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS
    if not args.aggregate_only:
        run_worker(instances, algos, seeds)
    else:
        aggregate()


if __name__ == '__main__':
    main()
