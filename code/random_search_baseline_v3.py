#!/usr/bin/env python3
"""
random_search_baseline_v3.py — Round 8: the sanity-check baseline an
external review correctly noted was missing for an 11-dimensional
mixed search space: does the evolutionary search actually do better
than i.i.d. uniform random sampling at the SAME evaluation budget
(8,320 oracle calls), or could a naive search reach a comparable
front?

Uses the exact same 11-variable decoder, oracle, bounds and evaluation
budget as the main comparison (run_baselines_comparison_v3.py), but
replaces NSGA-II/NSGA-III/etc.'s selection+crossover+mutation entirely
with i.i.d. uniform sampling over [XL, XU] -- no population dynamics,
no selection pressure, just 8,320 independent random draws per seed,
non-dominated-filtered at the end. 30 seeds, matching the main
comparison's statistical power.
"""
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from pymoo.util.nds.non_dominated_sorting import find_non_dominated

from bdcz_oracle_v3 import XL, XU, N_VAR, decode_topology, oracle_evaluate, dedup_front

SEEDS = list(range(42, 72))  # 30 seeds, matching the main comparison
N_EVALS = 8320  # matching the main comparison's per-run budget (NSGA-II/III/SMS-EMOA/RVEA)
HV_REF = np.array([1.1] * 4)

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step8_random_search_v3'
OUT_JSON = RESULTS / 'random_search_baseline_v3_results.json'


def eval_single(params):
    result = decode_topology(params)
    if result is None:
        return np.array([0.5, 0.0, 1.0, 0.0])
    eu, ev, et = result
    fid, rate_log, lat, cov = oracle_evaluate(eu, ev, et)
    return np.array([fid, rate_log, lat, cov])


def ckpt_path(s): return CKPT_DIR / f'seed{s}.pkl'
def done_path(s): return CKPT_DIR / f'seed{s}_completed.txt'
def is_done(s): return done_path(s).exists()


def run_worker(seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        if is_done(seed):
            print(f'  [SKIP] seed={seed}')
            continue
        print(f'  [RUN]  seed={seed} ...', end=' ', flush=True)
        t0 = time.time()
        rng = np.random.default_rng(seed)
        X = rng.uniform(XL, XU, size=(N_EVALS, N_VAR))
        F_raw = np.array([eval_single(x) for x in X])
        # minimisation convention matching the rest of this paper: -fid, -ratelog, lat, -cov
        F_min = np.column_stack([-F_raw[:, 0], -F_raw[:, 1], F_raw[:, 2], -F_raw[:, 3]])
        nd_mask = find_non_dominated(F_min)
        F_par, X_par = F_min[nd_mask], X[nd_mask]
        F_par, X_par = dedup_front(F_par, X_par)
        elapsed_min = (time.time() - t0) / 60.0
        payload = {'X_pareto': X_par, 'F_pareto': F_par, 'n_pareto': len(F_par),
                   'elapsed_min': elapsed_min, 'seed': seed, 'n_evals': N_EVALS}
        with open(ckpt_path(seed), 'wb') as fh:
            pickle.dump(payload, fh)
        with open(done_path(seed), 'w') as fh:
            fh.write(f'Completed at {datetime.now().isoformat()}\n')
        print(f'n_pareto={len(F_par)}, {elapsed_min:.1f} min')


def aggregate():
    missing = [s for s in SEEDS if not is_done(s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} seeds missing, e.g. {missing[:5]}')
        return

    runs = []
    for s in SEEDS:
        with open(ckpt_path(s), 'rb') as fh:
            runs.append(pickle.load(fh))

    all_F = np.vstack([r['F_pareto'] for r in runs if len(r['F_pareto']) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    hv_calc = Hypervolume(ref_point=HV_REF)
    hv_l, npar_l = [], []
    for r in runs:
        F = r['F_pareto']
        if len(F) == 0:
            hv_l.append(0.0); npar_l.append(0)
            continue
        hv_l.append(float(hv_calc(normalise(F))))
        npar_l.append(len(F))

    def ms(lst):
        a = np.array(lst)
        std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
        return {'mean': float(np.mean(a)), 'std': std}

    output = {
        'description': ('i.i.d. uniform random sampling over the same 11-variable '
                         'bounds, same 8,320-evaluation budget, same oracle, 30 seeds. '
                         'HV computed against this baseline\'s OWN pooled reference '
                         'front (not directly comparable to Table 3\'s HV without a '
                         'shared-reference re-score, same caveat as every other '
                         'representation-comparison block in this paper).'),
        'config': {'n_evals': N_EVALS, 'seeds': SEEDS},
        'summary': {'hv': ms(hv_l), 'n_pareto': ms(npar_l)},
    }
    with open(OUT_JSON, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_JSON}')
    print(f"  random_search HV={output['summary']['hv']['mean']:.4f}+/-{output['summary']['hv']['std']:.4f} "
          f"n_par={output['summary']['n_pareto']['mean']:.1f}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--seeds', type=str, default=None)
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS
    if not args.aggregate_only:
        run_worker(seeds)
    else:
        print('  [aggregate-only] skipping runs')
    aggregate()


if __name__ == '__main__':
    main()
