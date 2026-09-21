#!/usr/bin/env python3
"""
run_sigma_ablation_v3.py — external review round 7 (author-approved full
scope): does the decoder's 11th variable (seed_offset / "sigma") act as a
genuine design parameter, or mostly as a discontinuous selector among
pseudorandom realisations of the same physical decoder settings?

Variants, matched budget (NSGA-III, pop=104, 80 generations), matched
seeds:

A_sigma_optimised   — the actual 11-variable decoder (params[10]=seed_offset
    is a free variable, exactly as used everywhere else in this paper).
B_sigma_fixed_{0,42,500} — Round 8: genuinely $n_{\mathrm{var}}=10$: the
    pymoo Problem itself has 10 decision variables (SBX/PM crossover and
    mutation only ever see and manipulate 10 reals), and seed_offset is
    padded on afterwards, right before calling the oracle, as one of
    three fixed constants (0, 42, 500) tested independently rather than
    a single arbitrary choice -- Round 7's version kept n_var=11 with
    the 11th dimension's bounds collapsed to a single point, which is
    behaviourally equivalent (SBX/PM cannot move a variable whose
    xl==xu) but was correctly flagged as not what "10-variable" should
    mean at the level of the search-space definition itself.
C_sigma_averaged    — genuinely 10 free variables (same fix as B);
    fitness is the mean objective vector over K=3 independent
    realisations (seed_offset = 0, 1, 2) of the SAME physical
    parameters, decoupling the search from any single realisation's
    luck. Costs 3x the oracle calls of A/B per individual (tracked
    explicitly in the results, not just in evolutionary-evaluation
    count) -- reported so a reader can judge fairness under either an
    equal-evaluations or equal-oracle-calls framing.

If A >> B/C, part of the reported performance comes from optimising over
realisations, not physical topology parameters. If B/C ~= A, sigma can be
dropped for a cleaner 10-variable representation. If C > B, averaging
recovers most of A's advantage without needing sigma as a search variable.
"""
import argparse
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.algorithms.moo.nsga3 import NSGA3
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

from bdcz_oracle_v3 import XL as XL11, XU as XU11, N_VAR as N_VAR11, dedup_front
from bdcz_oracle_v3 import decode_topology, oracle_evaluate

POP_SIZE, N_GEN = 104, 80
SEEDS = list(range(42, 57))  # 15 seeds, matched budget, matched algorithm (NSGA-III)
HV_REF = np.array([1.1] * 4)
K_REALISATIONS = 3

OUT_DIR = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = OUT_DIR / 'step4_sigma_v3'
OUT_JSON = OUT_DIR / 'sigma_ablation_v3_results.json'


N_VAR10 = 10
XL10, XU10 = XL11[:10], XU11[:10]
N_ORACLE_CALLS_MULTIPLIER = {'A_sigma_optimised': 1, 'C_sigma_averaged': K_REALISATIONS}


def eval_optimised(params):
    """Variant A: standard 11-var decoder, seed_offset is params[10]."""
    result = decode_topology(params)
    if result is None:
        return np.array([0.5, 0.0, 1.0, 0.0])
    eu, ev, et = result
    return np.array(oracle_evaluate(eu, ev, et))


def make_eval_fixed(sigma_value):
    """Variant B family: genuinely 10 free variables (the pymoo Problem
    below has n_var=10; SBX/PM never see an 11th dimension at all).
    seed_offset is padded on as a fixed constant right before the oracle
    call, using an explicit rng so decode_topology's own params[10]
    read (always evaluated, but unused when rng is given) never fails
    on a padded value."""
    def eval_fixed(params_10):
        params_11 = np.append(params_10, sigma_value)
        result = decode_topology(params_11, rng=np.random.default_rng(int(sigma_value)))
        if result is None:
            return np.array([0.5, 0.0, 1.0, 0.0])
        eu, ev, et = result
        return np.array(oracle_evaluate(eu, ev, et))
    return eval_fixed


def eval_averaged(params_10):
    """Variant C: genuinely 10 free variables; fitness = mean objective
    vector over K_REALISATIONS independent decoder realisations of the
    same 10 physical parameters. Costs K_REALISATIONS oracle calls per
    individual (tracked via N_ORACLE_CALLS_MULTIPLIER), not 1."""
    vals = []
    for k in range(K_REALISATIONS):
        params_11 = np.append(params_10, 0.0)
        result = decode_topology(params_11, rng=np.random.default_rng(k))
        if result is None:
            vals.append(np.array([0.5, 0.0, 1.0, 0.0]))
        else:
            eu, ev, et = result
            vals.append(np.array(oracle_evaluate(eu, ev, et)))
    return np.mean(vals, axis=0)


FIXED_SIGMA_VALUES = [0, 42, 500]

VARIANTS = {
    'A_sigma_optimised': {'eval_fn': eval_optimised, 'n_var': N_VAR11,
                          'xl': XL11, 'xu': XU11},
    **{f'B_sigma_fixed_{s}': {'eval_fn': make_eval_fixed(s), 'n_var': N_VAR10,
                              'xl': XL10, 'xu': XU10}
       for s in FIXED_SIGMA_VALUES},
    'C_sigma_averaged': {'eval_fn': eval_averaged, 'n_var': N_VAR10,
                         'xl': XL10, 'xu': XU10},
}


class SigmaVariantProblem(Problem):
    def __init__(self, eval_fn, xl, xu, n_var):
        super().__init__(n_var=n_var, n_obj=4, n_constr=0, xl=xl, xu=xu)
        self.eval_fn = eval_fn

    def _evaluate(self, X, out, *args, **kwargs):
        F = np.zeros((len(X), 4))
        for i, params in enumerate(X):
            fid, rate_log, lat, cov = self.eval_fn(params)
            F[i] = [-fid, -rate_log, lat, -cov]
        out['F'] = F


def ckpt_path(v, s): return CKPT_DIR / f'{v}_seed{s}.pkl'
def done_path(v, s): return CKPT_DIR / f'{v}_seed{s}_completed.txt'
def is_done(v, s): return done_path(v, s).exists()


def run_worker(variants, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    ref_dirs = get_reference_directions('das-dennis', 4, n_partitions=6)
    for v in variants:
        cfg = VARIANTS[v]
        for seed in seeds:
            if is_done(v, seed):
                print(f'  [SKIP] {v} seed={seed}')
                continue
            print(f'  [RUN]  {v} seed={seed} ...', end=' ', flush=True)
            t0 = time.time()
            problem = SigmaVariantProblem(cfg['eval_fn'], cfg['xl'], cfg['xu'], cfg['n_var'])
            algo = NSGA3(pop_size=POP_SIZE, ref_dirs=ref_dirs, sampling=LHS(),
                        crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
            res = minimize(problem, algo, termination=MaximumGenerationTermination(N_GEN),
                            seed=seed, verbose=False, save_history=False)
            elapsed_min = (time.time() - t0) / 60.0
            if res.F is not None and len(res.F) > 0:
                pop_F, pop_X = res.pop.get('F'), res.pop.get('X')
                nd_mask = find_non_dominated(pop_F)
                X_par, F_par = pop_X[nd_mask], pop_F[nd_mask]
                F_par, X_par = dedup_front(F_par, X_par)
            else:
                X_par, F_par = np.empty((0, cfg['n_var'])), np.empty((0, 4))
            payload = {'X_pareto': X_par, 'F_pareto': F_par, 'n_pareto': len(F_par),
                       'elapsed_min': elapsed_min, 'seed': seed, 'variant': v}
            with open(ckpt_path(v, seed), 'wb') as fh:
                pickle.dump(payload, fh)
            with open(done_path(v, seed), 'w') as fh:
                fh.write(f'Completed at {datetime.now().isoformat()}\n')
            print(f'n_pareto={len(F_par)}, {elapsed_min:.1f} min')


def aggregate():
    missing = [(v, s) for v in VARIANTS for s in SEEDS if not is_done(v, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} (variant, seed) combos missing, e.g. {missing[:5]}')
        return

    variant_fronts = {}
    for v in VARIANTS:
        variant_fronts[v] = []
        for s in SEEDS:
            with open(ckpt_path(v, s), 'rb') as fh:
                run = pickle.load(fh)
            variant_fronts[v].append(run)

    all_F = np.vstack([r['F_pareto'] for v in VARIANTS for r in variant_fronts[v] if len(r['F_pareto']) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_all = find_non_dominated(all_F)
    ref_front = normalise(np.unique(np.round(all_F[nd_all], 6), axis=0))
    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

    summary = {}
    for v in VARIANTS:
        hv_l, igd_l, npar_l, elapsed_l = [], [], [], []
        for r in variant_fronts[v]:
            F = r['F_pareto']
            if len(F) == 0:
                continue
            Fn = normalise(F)
            hv_l.append(float(hv_calc(Fn)))
            igd_l.append(float(igd_calc(Fn)))
            npar_l.append(len(F))
            elapsed_l.append(r['elapsed_min'])

        def ms(lst):
            a = np.array(lst)
            std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
            return {'mean': float(np.mean(a)), 'std': std}

        oracle_call_multiplier = N_ORACLE_CALLS_MULTIPLIER.get(v, 1)
        n_evolutionary_evals = POP_SIZE * N_GEN
        summary[v] = {'hv': ms(hv_l), 'igd': ms(igd_l), 'n_pareto': ms(npar_l),
                      'elapsed_min': ms(elapsed_l),
                      'n_evolutionary_evaluations': n_evolutionary_evals,
                      'oracle_calls_per_evaluation': oracle_call_multiplier,
                      'total_oracle_calls': n_evolutionary_evals * oracle_call_multiplier}

    output = {
        'config': {'algorithm': 'nsga3', 'pop_size': POP_SIZE, 'n_gen': N_GEN,
                   'seeds': SEEDS, 'k_realisations': K_REALISATIONS,
                   'fixed_sigma_values_tested': FIXED_SIGMA_VALUES,
                   'variants': list(VARIANTS.keys())},
        'summary': summary,
    }
    with open(OUT_JSON, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    for v in VARIANTS:
        m = summary[v]
        print(f"  {v:<20} HV={m['hv']['mean']:.4f}+/-{m['hv']['std']:.4f} "
              f"IGD={m['igd']['mean']:.4f}+/-{m['igd']['std']:.4f} n_par={m['n_pareto']['mean']:.1f}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--variants', type=str, default=None)
    p.add_argument('--seeds', type=str, default=None)
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
