#!/usr/bin/env python3
"""
categorical_operators_baseline_v3.py — external review round 7
(author-approved full scope): a SECOND categorical baseline, using
operators native to a categorical representation (uniform categorical
crossover; add-edge, remove-edge, change-medium mutation) instead of the
SBX/polynomial-mutation real-valued relaxation used in binary_baseline_v3.py.

Purpose: binary_baseline_v3.py's poor performance shows a *relaxed*
categorical encoding under *continuous* operators struggles at this scale;
it does not by itself show that *no* categorical encoding/operator
combination could compete with the compact parametric decoder. This
script tests that directly with operators designed for discrete categories.

Representation: integer array of length 4,950 (one per candidate city
pair), values in {0=no edge, 1=fibre, 2=satellite, 3=FSO} -- genuinely
discrete throughout (no real-valued relaxation/rounding).

Tested at 4 initial densities: 3%, 10%, 15%, and 10.9% (the decoder's own
median edge density across the global Pareto front, verified directly
from global_pareto_v3.npz's n_edges field: median 540.5/4950).

NSGA-II only (matching the faster of the two algorithms in the original
categorical baseline), matched budget (pop=104, 80 generations).
"""
import argparse
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.problem import Problem
from pymoo.core.sampling import Sampling
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.optimize import minimize
from pymoo.util.nds.non_dominated_sorting import find_non_dominated

from bdcz_oracle_v3 import oracle_evaluate, N_NODES, dedup_front

N_EDGES = N_NODES * (N_NODES - 1) // 2  # 4950
POP_SIZE, N_GEN = 104, 80
HV_REF = np.array([1.1] * 4)
SEEDS = list(range(42, 52))  # 10 seeds, matched to the original categorical
                             # baseline's Round-4 seed count, per the reviewer's
                             # suggested 10-15 seeds for this secondary baseline.
DENSITIES = {'d03': 0.03, 'd10': 0.10, 'd11_decoder_median': 0.109, 'd15': 0.15}

OUT_DIR = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = OUT_DIR / 'step5_catops_v3'
OUT_JSON = OUT_DIR / 'categorical_operators_v3_results.json'

EDGE_MAP = np.zeros((N_EDGES, 2), dtype=np.int64)
_idx = 0
for _i in range(N_NODES):
    for _j in range(_i + 1, N_NODES):
        EDGE_MAP[_idx] = [_i, _j]
        _idx += 1


class CategoricalDensitySampling(Sampling):
    """Genuinely discrete initial population at a fixed target density,
    medium drawn uniformly among {fibre, satellite, FSO} for active edges."""
    def __init__(self, density):
        super().__init__()
        self.density = density

    def _do(self, problem, n_samples, **kwargs):
        X = np.zeros((n_samples, problem.n_var), dtype=np.int64)
        active = np.random.random((n_samples, problem.n_var)) < self.density
        cat = np.random.randint(1, 4, size=(n_samples, problem.n_var))
        X[active] = cat[active]
        return X


class UniformCategoricalCrossover(Crossover):
    """Standard uniform crossover directly on category labels: each gene
    independently inherited from parent 1 or 2 with probability 0.5."""
    def __init__(self):
        super().__init__(2, 2)

    def _do(self, problem, X, **kwargs):
        _, n_matings, n_var = X.shape
        mask = np.random.random((n_matings, n_var)) < 0.5
        Y = np.zeros_like(X)
        Y[0][mask] = X[0][mask]
        Y[0][~mask] = X[1][~mask]
        Y[1][mask] = X[1][mask]
        Y[1][~mask] = X[0][~mask]
        return Y


class EdgeMutation(Mutation):
    """Three category-native mutation operators, each with its own
    per-gene probability: add-edge (0 -> random medium), remove-edge
    (nonzero -> 0), change-medium (nonzero -> a DIFFERENT nonzero
    category). Probabilities chosen so the expected number of mutated
    genes per individual is comparable to a standard 1/n_var PM rate."""
    def __init__(self, p_add=1.0 / N_EDGES, p_remove=1.0 / N_EDGES, p_change=1.0 / N_EDGES):
        super().__init__()
        self.p_add, self.p_remove, self.p_change = p_add, p_remove, p_change

    def _do(self, problem, X, **kwargs):
        X = X.copy()
        is_zero = X == 0
        is_nonzero = ~is_zero

        add_mask = is_zero & (np.random.random(X.shape) < self.p_add)
        X[add_mask] = np.random.randint(1, 4, size=add_mask.sum())

        remove_mask = is_nonzero & (np.random.random(X.shape) < self.p_remove)
        X[remove_mask] = 0

        change_mask = is_nonzero & ~remove_mask & (np.random.random(X.shape) < self.p_change)
        n_change = change_mask.sum()
        if n_change > 0:
            old_vals = X[change_mask]
            new_vals = np.random.randint(1, 3, size=n_change)  # 1 or 2 -> remapped below
            new_vals = np.where(new_vals >= old_vals, new_vals + 1, new_vals)  # skip old value, stay in {1,2,3}
            X[change_mask] = new_vals
        return X


class CategoricalTopologyProblemV2(Problem):
    def __init__(self):
        super().__init__(n_var=N_EDGES, n_obj=4, n_constr=0, xl=0, xu=3)

    def _evaluate(self, X, out, *args, **kwargs):
        F = np.zeros((len(X), 4))
        for row_idx, row in enumerate(X.astype(np.int64)):
            active_idx = np.where(row > 0)[0]
            if len(active_idx) < 5:
                F[row_idx] = [0.0, 0.0, 1e6, 0.0]
                continue
            eu = EDGE_MAP[active_idx, 0].copy()
            ev = EDGE_MAP[active_idx, 1].copy()
            et = (row[active_idx] - 1).astype(np.int64)
            try:
                fid, rate_log, lat, cov = oracle_evaluate(eu, ev, et)
                F[row_idx] = [-fid, -rate_log, lat, -cov]
            except Exception:
                F[row_idx] = [0.0, 0.0, 1e6, 0.0]
        out['F'] = F


def ckpt_path(d, s): return CKPT_DIR / f'{d}_seed{s}.pkl'
def done_path(d, s): return CKPT_DIR / f'{d}_seed{s}_completed.txt'
def is_done(d, s): return done_path(d, s).exists()


def run_worker(densities, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    for dname in densities:
        density = DENSITIES[dname]
        for seed in seeds:
            if is_done(dname, seed):
                print(f'  [SKIP] {dname} seed={seed}')
                continue
            print(f'  [RUN]  {dname} (density={density}) seed={seed} ...', end=' ', flush=True)
            t0 = time.time()
            np.random.seed(seed)
            problem = CategoricalTopologyProblemV2()
            algo = NSGA2(pop_size=POP_SIZE, sampling=CategoricalDensitySampling(density),
                        crossover=UniformCategoricalCrossover(), mutation=EdgeMutation())
            res = minimize(problem, algo, ('n_gen', N_GEN), seed=seed, verbose=False)
            elapsed_min = (time.time() - t0) / 60.0
            if res.F is not None and len(res.F) > 0:
                nd_mask = find_non_dominated(res.F)
                F_par = res.F[nd_mask]
                F_par = np.unique(np.round(F_par, 6), axis=0)
            else:
                F_par = np.empty((0, 4))
            entry = {'seed': seed, 'density_name': dname, 'density': density,
                     'n_pareto': len(F_par), 'elapsed_min': elapsed_min,
                     'pareto_front': F_par.tolist()}
            with open(ckpt_path(dname, seed), 'w') as fh:
                json.dump(entry, fh)
            with open(done_path(dname, seed), 'w') as fh:
                fh.write(f'Completed at {datetime.now().isoformat()}\n')
            print(f'n_pareto={len(F_par)}, {elapsed_min:.1f} min')


def aggregate():
    densities = list(DENSITIES.keys())
    missing = [(d, s) for d in densities for s in SEEDS if not is_done(d, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} (density, seed) combos missing, e.g. {missing[:5]}')
        return

    all_fronts = {}
    for d in densities:
        all_fronts[d] = []
        for s in SEEDS:
            with open(ckpt_path(d, s)) as fh:
                entry = json.load(fh)
            all_fronts[d].append(np.array(entry['pareto_front']).reshape(-1, 4))

    combined = np.vstack([F for fs in all_fronts.values() for F in fs if len(F) > 0])
    g_min, g_max = combined.min(axis=0), combined.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_global = find_non_dominated(combined)
    ref_front = normalise(np.unique(np.round(combined[nd_global], 6), axis=0))
    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

    summary = {}
    for d in densities:
        hv_l, igd_l, npar_l, el_l = [], [], [], []
        for F in all_fronts[d]:
            if len(F) == 0:
                continue
            Fn = normalise(F)
            hv_l.append(float(hv_calc(Fn)))
            igd_l.append(float(igd_calc(Fn)))
            npar_l.append(len(F))
        def ms(lst):
            a = np.array(lst)
            std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
            return {'mean': float(np.mean(a)), 'std': std} if len(a) else {'mean': None, 'std': None}
        summary[d] = {'hv': ms(hv_l), 'igd': ms(igd_l), 'n_pareto': ms(npar_l)}

    output = {'config': {'algorithm': 'nsga2', 'pop_size': POP_SIZE, 'n_gen': N_GEN,
                         'seeds': SEEDS, 'densities': DENSITIES},
              'summary': summary}
    with open(OUT_JSON, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    for d in densities:
        m = summary[d]
        print(f"  {d:<24} HV={m['hv']['mean']:.4f}+/-{m['hv']['std']:.4f} n_par={m['n_pareto']['mean']:.1f}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--densities', type=str, default=None)
    p.add_argument('--seeds', type=str, default=None)
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()
    densities = args.densities.split(',') if args.densities else list(DENSITIES.keys())
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS
    if not args.aggregate_only:
        run_worker(densities, seeds)
    else:
        aggregate()


if __name__ == '__main__':
    main()
