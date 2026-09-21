#!/usr/bin/env python3
"""
run_convergence_v3.py — external review round 7 (author-approved full
scope): per-generation HV/IGD tracking for all 5 algorithms, addressing
the recurring reviewer request for convergence curves (the main
comparison used save_history=False, so only the final-generation front
was ever recorded).

Uses a REDUCED seed count (10, seeds 42-51) rather than the main
comparison's 30, since this is a secondary characterisation of
convergence behaviour, not a new headline statistical comparison
(matching the reviewer's own suggestion of 10-15 seeds for this kind of
supplementary analysis). Same algorithms, same matched-budget generation
counts, same population size as the main comparison
(run_baselines_comparison_v3.py, reused directly here for make_algorithm
and n_gen_for).

HV/IGD at each generation are computed against the SAME fixed reference
(HV_REF=1.1 in each objective; IGD reference front = the released 561-point
global Pareto front / pooled evolutionary front, numerically identical to
the main comparison's own reference front) and the SAME min-max
normalisation bounds already established in
unified_seven_algorithm_comparison_v3.json, so intermediate-generation
values are directly comparable to the main comparison's final-generation
numbers -- not a separately-normalised, non-comparable quantity.

FIX (round 14, ChatGPT Round 13 review finding #4, confirmed real): this
script previously sourced its normalisation from baselines_comparison_v3.json
and its IGD reference from global_pareto_v3.npz while that file still only
pooled 5 algorithms -- both stale relative to the manuscript's actual "main
comparison" (Tables 3-7), which since Round 13 is the 7-algorithm
unified_seven_algorithm_comparison_v3.json (g_min/g_max built from the union
of all 210 per-seed fronts) and the rebuilt 561-point global_pareto_v3.npz.
Figure 2's reported IGD values were therefore on a different scale than the
text describing them, and did not match the union used everywhere else. This
requires RE-RUNNING the archived optimisations (not just reprocessing): the
per-generation HV/IGD are computed live inside ConvergenceCallback.notify()
and only the resulting numbers -- not the raw per-generation archive -- are
checkpointed, so correcting the normalisation constants requires regenerating
step6_convergence_v3/*.pkl. The population dynamics themselves are seeded and
unaffected by this change (only the logged metric scale differs), so this is
a mechanical re-run, not a new experiment.
"""
import argparse
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from pymoo.core.callback import Callback
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.util.ref_dirs import get_reference_directions

from run_baselines_comparison_v3 import make_algorithm, n_gen_for, ALGORITHMS, POP_SIZE, N_GEN, QuantumTopoProblem
from bdcz_oracle_v3 import dedup_front

SEEDS = list(range(42, 52))  # 10 seeds, reduced from the main comparison's 30
                             # -- this is a convergence-behaviour characterisation,
                             # not a new headline statistical comparison.
HV_REF = np.array([1.1] * 4)

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = RESULTS / 'step6_convergence_v3'
OUT_JSON = RESULTS / 'convergence_v3_results.json'

with open(RESULTS / 'unified_seven_algorithm_comparison_v3.json') as f:
    _base = json.load(f)
G_MIN = np.array(_base['normalisation']['g_min'])
_G_MAX = np.array(_base['normalisation']['g_max'])
G_RANGE = np.where((_G_MAX - G_MIN) < 1e-12, 1e-9, _G_MAX - G_MIN)

_gp = np.load(RESULTS / 'global_pareto_v3.npz', allow_pickle=True)
_REF_RAW = _gp['F_pareto']


def normalise(F):
    return (F - G_MIN) / G_RANGE


_REF_FRONT_NORM = normalise(_REF_RAW)
_hv_calc = Hypervolume(ref_point=HV_REF)
_igd_calc = IGD(_REF_FRONT_NORM)


class ConvergenceCallback(Callback):
    """Tracks a genuine best-so-far archive: the non-dominated frontier
    of EVERY solution evaluated up to and including the current
    generation, not just the current population's own front
    (algorithm.opt/pop). Round 8 fix: the previous version read
    algorithm.opt each generation, which for a finite population can
    occasionally lose a previously-found non-dominated point to
    crowding/truncation, producing small non-monotonic dips in HV/IGD
    that are a population-dynamics artefact, not a property of the
    search's actual best-found front. Maintaining an explicit archive
    (union of every generation's population, re-filtered to
    non-dominated, every generation) guarantees HV is non-decreasing
    and archive-IGD is non-increasing by construction, since the
    non-dominated frontier of a growing candidate pool can only match
    or improve on the previous frontier."""

    def __init__(self):
        super().__init__()
        self.trace = []  # list of (n_eval, hv, igd, n_unique)
        self.archive_F = np.empty((0, 4))

    def notify(self, algorithm):
        F = algorithm.pop.get('F')
        combined = np.vstack([self.archive_F, F]) if len(self.archive_F) > 0 else F
        nd_mask = find_non_dominated(combined)
        F_nd = combined[nd_mask]
        F_dedup, _ = dedup_front(F_nd, F_nd)
        self.archive_F = F_dedup

        Fn = normalise(F_dedup)
        hv = float(_hv_calc(Fn)) if len(Fn) > 0 else 0.0
        igd = float(_igd_calc(Fn)) if len(Fn) > 0 else float('inf')
        self.trace.append({
            'n_eval': int(algorithm.evaluator.n_eval),
            'hv': hv, 'igd': igd, 'n_unique': int(len(F_dedup)),
        })


def ckpt_path(algo, s): return CKPT_DIR / f'{algo}_seed{s}.pkl'
def done_path(algo, s): return CKPT_DIR / f'{algo}_seed{s}_completed.txt'
def is_done(algo, s): return done_path(algo, s).exists()


def run_worker(algos, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    ref_dirs = get_reference_directions('das-dennis', 4, n_partitions=6)
    for algo in algos:
        n_gen_algo = n_gen_for(algo, ref_dirs)
        for seed in seeds:
            if is_done(algo, seed):
                print(f'  [SKIP] {algo} seed={seed}')
                continue
            print(f'  [RUN]  {algo} seed={seed} (n_gen={n_gen_algo}) ...', end=' ', flush=True)
            t0 = time.time()
            algorithm = make_algorithm(algo, ref_dirs)
            cb = ConvergenceCallback()
            problem = QuantumTopoProblem()
            res = minimize(problem, algorithm,
                            termination=MaximumGenerationTermination(n_gen_algo),
                            seed=seed, verbose=False, save_history=False, callback=cb)
            elapsed_min = (time.time() - t0) / 60.0
            payload = {'algo': algo, 'seed': seed, 'trace': cb.trace, 'elapsed_min': elapsed_min}
            with open(ckpt_path(algo, seed), 'wb') as fh:
                pickle.dump(payload, fh)
            with open(done_path(algo, seed), 'w') as fh:
                fh.write(f'Completed at {datetime.now().isoformat()}\n')
            print(f'{len(cb.trace)} generations recorded, {elapsed_min:.1f} min')


def _bootstrap_ci(values, n_boot=2000, alpha=0.05, rng=None):
    rng = rng or np.random.default_rng(0)
    values = np.asarray(values)
    boots = np.array([np.median(rng.choice(values, size=len(values), replace=True))
                       for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def aggregate():
    missing = [(a, s) for a in ALGORITHMS for s in SEEDS if not is_done(a, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} (algo, seed) combos missing, e.g. {missing[:5]}')
        return

    out = {'config': {'algorithms': ALGORITHMS, 'seeds': SEEDS, 'pop_size': POP_SIZE},
           'per_algorithm': {}}
    for algo in ALGORITHMS:
        traces = []
        for s in SEEDS:
            with open(ckpt_path(algo, s), 'rb') as fh:
                traces.append(pickle.load(fh)['trace'])
        # Align on n_eval (all seeds of a given algorithm share the same
        # generation schedule/pop size, so n_eval values coincide exactly).
        n_gens = min(len(t) for t in traces)
        curve = []
        for g in range(n_gens):
            n_eval = traces[0][g]['n_eval']
            hv_vals = [t[g]['hv'] for t in traces]
            igd_vals = [t[g]['igd'] for t in traces if np.isfinite(t[g]['igd'])]
            hv_lo, hv_hi = _bootstrap_ci(hv_vals)
            entry = {'n_eval': n_eval, 'hv_median': float(np.median(hv_vals)),
                     'hv_ci_lo': hv_lo, 'hv_ci_hi': hv_hi}
            if igd_vals:
                igd_lo, igd_hi = _bootstrap_ci(igd_vals)
                entry.update({'igd_median': float(np.median(igd_vals)),
                              'igd_ci_lo': igd_lo, 'igd_ci_hi': igd_hi})
            curve.append(entry)
        out['per_algorithm'][algo] = curve
    with open(OUT_JSON, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'Saved: {OUT_JSON}')


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
    else:
        aggregate()


if __name__ == '__main__':
    main()
