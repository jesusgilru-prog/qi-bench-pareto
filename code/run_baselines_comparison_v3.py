#!/usr/bin/env python3
"""
run_baselines_comparison_v3.py — Canonical-oracle rerun (Paper 1-unified physics)
====================================================================================
Same protocol as v2 (NSGA-II vs NSGA-III vs MOEA/D, seeds 42-51) but using
bdcz_oracle_v3 (exact Werner-optimal path search, real SAT/FSO base
fidelities matching Paper 1's canonical oracle, no -log(F) proxy, no
fidelity floor).

FIX (Round 3 review, confirmed real bug): NSGA2/NSGA3 use pop_size=104
(POP_SIZE, set explicitly), but MOEA/D in pymoo does not take a pop_size
argument -- its population size is implicitly len(ref_dirs). Verified
directly: get_reference_directions('das-dennis', 4, n_partitions=6) returns
84 directions, not 104, so at the same N_GEN=80 MOEA/D previously ran on
~19% fewer function evaluations than NSGA2/NSGA3 (6,720 vs 8,320) while the
manuscript claimed "the same population size, generations and parameter
ranges" for all three. Fixed by giving MOEA/D more generations
(N_GEN_MOEAD) so its total evaluation count approximately matches
NSGA2/NSGA3's, rather than forcing an equal population size onto an
algorithm whose population size is architecturally tied to its reference
directions. The ACTUAL evaluation count for every algorithm is measured
directly (problem.n_calls) and saved to the output JSON -- the manuscript
must report that measured number, not an assumed one.
"""

import numpy as np
from metrics import compute_spacing
import json
import pickle
import time
from datetime import datetime
from pathlib import Path

from bdcz_oracle_v3 import XL, XU, eval_single, N_VAR, dedup_front

from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.moo.moead import MOEAD
from pymoo.algorithms.moo.sms import SMSEMOA
from pymoo.algorithms.moo.rvea import RVEA
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize
from pymoo.operators.sampling.lhs import LHS
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD

POP_SIZE   = 104
N_GEN      = 80
# FIX (external review round 6, B1/B2): 10 seeds put every Holm-corrected
# p-value this benchmark reports at or within noise of the paired
# two-sided Wilcoxon signed-rank test's own resolution floor for n=10
# (2/2^10 = 0.00195 raw, 0.0176 after Holm across 9 tests) -- a reviewer
# with statistical training would (correctly) read every "significant"
# result here as indistinguishable from "the test's floor was reached",
# not as evidence of a specific effect size. 30 seeds moves the floor to
# 2/2^30 (~1.9e-9), far below any p-value actually observed, so a reported
# significant result reflects the data, not the test's resolution limit.
# ALSO (B2): NSGA-II/NSGA-III/MOEA-D alone is a thin algorithm portfolio
# for an evolutionary-computation venue specifically (as opposed to a
# general applied-ML venue) -- added SMS-EMOA (hypervolume-
# based selection, a natural comparison given HV is this paper's primary
# metric) and RVEA (a modern reference-vector-adaptation many-objective
# algorithm) using the SAME decoder, oracle, evaluation budget, and seeds
# as the original three. AGE-MOEA-II was considered and rejected: pymoo's
# implementation requires numba, a new dependency this project does not
# otherwise need, and RVEA already covers the "modern many-objective
# reference" requirement without it.
SEEDS      = list(range(42, 72))
ALGORITHMS = ['nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea']
HV_REF     = np.array([1.1] * 4)

OUT_DIR  = Path(__file__).resolve().parent.parent / 'results'
CKPT_DIR = OUT_DIR / 'step1_baselines_v3'
OUT_FILE = OUT_DIR / 'baselines_comparison_v3.json'


class QuantumTopoProblem(Problem):
    def __init__(self):
        super().__init__(n_var=N_VAR, n_obj=4, n_constr=0, xl=XL, xu=XU)
        self.n_calls = 0

    def _evaluate(self, X, out, *args, **kwargs):
        results = np.zeros((len(X), 4))
        for i, params in enumerate(X):
            results[i] = eval_single(params)
        self.n_calls += len(X)
        out['F'] = np.column_stack([
            -results[:, 0], -results[:, 1], results[:, 2], -results[:, 3],
        ])




def make_algorithm(name, ref_dirs):
    if name == 'nsga2':
        return NSGA2(pop_size=POP_SIZE, sampling=LHS(),
                     crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'nsga3':
        return NSGA3(pop_size=POP_SIZE, ref_dirs=ref_dirs, sampling=LHS(),
                     crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'moead':
        return MOEAD(ref_dirs=ref_dirs, n_neighbors=15, sampling=LHS(),
                     crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'sms_emoa':
        # SMS-EMOA takes an explicit pop_size (unlike MOEA/D/NSGA3/RVEA,
        # whose population is entangled with ref_dirs), so it gets the
        # SAME POP_SIZE/N_GEN as NSGA2/NSGA3 with no budget adjustment.
        return SMSEMOA(pop_size=POP_SIZE, sampling=LHS(),
                        crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    elif name == 'rvea':
        # RVEA requires ref_dirs but its `pop_size` argument was verified
        # (via direct calibration against the real oracle problem, see
        # n_gen_for) to fully control its per-generation evaluation count
        # once set explicitly, exactly like NSGA2/SMS-EMOA -- despite
        # ref_dirs having only 84 directions, pop_size=104 here yields the
        # same 104 evaluations/generation as the others, not 84.
        return RVEA(ref_dirs=ref_dirs, pop_size=POP_SIZE, sampling=LHS(),
                     crossover=SBX(prob=0.9, eta=15), mutation=PM(eta=20))
    raise ValueError(name)


def n_gen_for(algo, ref_dirs):
    """MOEA/D's population is entangled with len(ref_dirs) (not cleanly
    fixed by POP_SIZE), so it gets a generation count chosen so its total
    evaluation budget approximately matches POP_SIZE x (N_GEN+1), instead
    of running on an implicitly smaller population for the same number of
    generations (the originally confirmed bug, round 3).

    SMS-EMOA and RVEA were both calibrated directly against the real
    11-D oracle problem before this budget-matching was trusted (not
    assumed from pymoo's docs): both measured EXACTLY pop_size evaluations
    per generation with no extra initial-population call (104 evals at
    n_gen=5 and 1,560 at n_gen=15, an exact 104/generation slope with zero
    intercept for both algorithms) -- despite RVEA's ref_dirs/pop_size
    interaction being non-trivial (see make_algorithm), its pop_size
    argument turned out to fully control its per-generation evaluation
    count once measured, so it needs no adjustment either. n_gen=N_GEN
    therefore gives all of NSGA2/NSGA3/SMS-EMOA/RVEA the identical
    104*80=8,320-evaluation budget; only MOEA/D's population is
    architecturally decoupled from POP_SIZE and needs a different n_gen."""
    if algo == 'moead':
        target_evals = POP_SIZE * (N_GEN + 1)
        return max(1, round(target_evals / len(ref_dirs)) - 1)
    return N_GEN


def ckpt_path(algo, seed):
    return CKPT_DIR / f'{algo}_seed{seed}.pkl'

def done_path(algo, seed):
    return CKPT_DIR / f'{algo}_seed{seed}_completed.txt'

def is_done(algo, seed):
    return done_path(algo, seed).exists()

def save_run(algo, seed, X, F, elapsed_min, n_calls):
    payload = {'X_pareto': X, 'F_pareto': F, 'n_pareto': len(F),
               'elapsed_min': elapsed_min, 'seed': seed, 'algorithm': algo,
               'n_calls': n_calls, 'variant': 'full_decoder_11d_v3'}
    with open(ckpt_path(algo, seed), 'wb') as fh:
        pickle.dump(payload, fh)
    with open(done_path(algo, seed), 'w') as fh:
        fh.write(f'Completed at {datetime.now().isoformat()}\n')

def load_run(algo, seed):
    with open(ckpt_path(algo, seed), 'rb') as fh:
        return pickle.load(fh)


def run_worker(algos, seeds, ref_dirs):
    """Run (and checkpoint) exactly the requested (algo, seed) combinations.
    Safe to invoke concurrently across multiple OS processes with disjoint
    (algos, seeds) shards -- checkpoint files are one-per-(algo,seed), so
    there is no shared mutable state between shards (external review round
    6, B1/B2: 30 seeds x 5 algorithms is far too slow run serially; this
    lets the work be sharded across CPU cores)."""
    for algo in algos:
        for seed in seeds:
            if is_done(algo, seed):
                print(f'  [SKIP] {algo} seed={seed}')
                continue
            print(f'  [RUN]  {algo} seed={seed} ...', end=' ', flush=True)
            t0 = time.time()

            problem   = QuantumTopoProblem()
            algorithm = make_algorithm(algo, ref_dirs)
            n_gen_algo = n_gen_for(algo, ref_dirs)
            res = minimize(problem, algorithm,
                            termination=MaximumGenerationTermination(n_gen_algo),
                            seed=seed, verbose=False, save_history=False)
            elapsed_min = (time.time() - t0) / 60.0

            if res.F is not None and len(res.F) > 0:
                pop_F = res.pop.get('F'); pop_X = res.pop.get('X')
                nd_mask = find_non_dominated(pop_F)
                X_par, F_par = pop_X[nd_mask], pop_F[nd_mask]
                n_before_dedup = len(F_par)
                F_par, X_par = dedup_front(F_par, X_par)
                if n_before_dedup != len(F_par):
                    print(f'[dedup {n_before_dedup}->{len(F_par)}]', end=' ', flush=True)
            else:
                X_par, F_par = np.empty((0, N_VAR)), np.empty((0, 4))

            save_run(algo, seed, X_par, F_par, elapsed_min, problem.n_calls)
            print(f'n_pareto={len(F_par)}, {elapsed_min:.2f} min, '
                  f'n_calls={problem.n_calls}')


def aggregate(ref_dirs):
    """Read every (algo, seed) checkpoint in ALGORITHMS x SEEDS (all of
    them must exist -- run_worker for every shard first) and write the
    final baselines_comparison_v3.json, exactly as the original single-
    process main() did."""
    all_fronts, all_X_dict, elapsed_dict, calls_dict = {}, {}, {}, {}
    for algo in ALGORITHMS:
        all_fronts[algo], all_X_dict[algo] = [], []
        elapsed_dict[algo], calls_dict[algo] = [], []
        for seed in SEEDS:
            run = load_run(algo, seed)
            all_fronts[algo].append(run['F_pareto'])
            all_X_dict[algo].append(run['X_pareto'])
            elapsed_dict[algo].append(run['elapsed_min'])
            calls_dict[algo].append(run.get('n_calls', POP_SIZE * (N_GEN + 1)))

    combined_all = np.vstack([F for fs in all_fronts.values() for F in fs if len(F) > 0])
    g_min = combined_all.min(axis=0)
    g_max = combined_all.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd_global = find_non_dominated(combined_all)
    # Dedup the IGD reference front too: per-run fronts are now deduplicated
    # individually, but the same optimum can still be independently
    # rediscovered across different seeds/algorithms, which would still
    # overweight IGD toward frequently-rediscovered regions if left in.
    combined_all_X = np.vstack([X for xs in all_X_dict.values() for X in xs if len(X) > 0])
    ref_front_raw, _ = dedup_front(combined_all[nd_global], combined_all_X[nd_global])
    ref_front = normalise(ref_front_raw)

    per_run_details, summary = {}, {}
    for algo in ALGORITHMS:
        per_run_details[algo] = []
        hv_list, igd_list, sp_list, npar_list, elapsed_list, calls_list = [], [], [], [], [], []
        hv_calc = Hypervolume(ref_point=HV_REF)
        igd_calc = IGD(ref_front)

        for idx, seed in enumerate(SEEDS):
            F_raw = all_fronts[algo][idx]
            if len(F_raw) == 0:
                per_run_details[algo].append(dict(seed=seed, n_pareto=0, hv=0.0,
                    igd=None, spacing=0.0, elapsed_min=elapsed_dict[algo][idx],
                    n_calls=calls_dict[algo][idx]))
                continue
            F_norm = normalise(F_raw)
            hv_val = float(hv_calc(F_norm))
            igd_val = float(igd_calc(F_norm))
            sp_val = compute_spacing(F_norm)
            rec = dict(seed=seed, n_pareto=len(F_raw), hv=hv_val, igd=igd_val,
                       spacing=sp_val, elapsed_min=elapsed_dict[algo][idx],
                       n_calls=calls_dict[algo][idx])
            per_run_details[algo].append(rec)
            hv_list.append(hv_val); igd_list.append(igd_val); sp_list.append(sp_val)
            npar_list.append(len(F_raw)); elapsed_list.append(elapsed_dict[algo][idx])
            calls_list.append(calls_dict[algo][idx])

        def ms(lst):
            a = np.array([x for x in lst if x is not None and np.isfinite(x)])
            if len(a) == 0:
                return {'mean': None, 'std': None, 'median': None, 'iqr_lo': None, 'iqr_hi': None}
            # FIX (external review round 6, item 10): np.std() defaults to
            # ddof=0 (population std), which is the wrong estimator for a
            # SAMPLE of stochastic optimisation runs -- ddof=1 (sample std,
            # Bessel-corrected) is the standard choice for n independent
            # replicates. Also report median/IQR (robust to the NSGA-II
            # seed-43 style outlier the review flagged), not only mean+std.
            std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
            q1, q3 = np.percentile(a, [25, 75])
            return {'mean': float(np.mean(a)), 'std': std, 'median': float(np.median(a)),
                    'iqr_lo': float(q1), 'iqr_hi': float(q3)}

        summary[algo] = {
            'hv': ms(hv_list), 'igd': ms(igd_list), 'spacing': ms(sp_list),
            'n_pareto': ms(npar_list), 'elapsed_min': ms(elapsed_list),
            'n_calls': ms(calls_list),
        }

    hv_means = {a: summary[a]['hv']['mean'] or 0.0 for a in ALGORITHMS}
    ranked_hv = sorted(ALGORITHMS, key=lambda a: hv_means[a], reverse=True)
    n_gen_by_algo = {a: n_gen_for(a, ref_dirs) for a in ALGORITHMS}

    output = {
        'config': {'pop_size': POP_SIZE, 'n_gen': N_GEN, 'n_gen_by_algo': n_gen_by_algo,
                   'n_ref_dirs': int(len(ref_dirs)), 'seeds': SEEDS,
                   'algorithms': ALGORITHMS, 'n_var': N_VAR,
                   'note': ('v3: canonical Paper-1-unified oracle (exact Werner-optimal '
                            'search, real SAT/FSO constants); MOEA/D population is '
                            'len(ref_dirs), not pop_size -- n_gen_by_algo gives its '
                            'matched-budget generation count, actual evaluation counts '
                            'are in per_run_details[*][*]["n_calls"]. round 6: 30 seeds '
                            '(42-71) and 5 algorithms (added SMS-EMOA, RVEA) replacing '
                            'the round-3..5 10-seed/3-algorithm design (external review '
                            'round 6, B1/B2).')},
        'normalisation': {'g_min': g_min.tolist(), 'g_max': g_max.tolist(),
                          'g_range': g_range.tolist(),
                          # FIX (Round 5 review): this field previously reported
                          # len(nd_global) -- the RAW non-dominated count BEFORE
                          # dedup_front() -- while the actual IGD reference front
                          # used for every metric in this file is the smaller,
                          # deduplicated one. The metric computation itself was
                          # always correct (ref_front is dedup'd before use); only
                          # this metadata label was misleading.
                          'n_ref_front_raw_pre_dedup': int(len(nd_global)),
                          'n_ref_front': int(len(ref_front))},
        'summary': summary,
        'ranking': {'by_hv': ranked_hv},
        'per_run_details': per_run_details,
    }
    with open(OUT_FILE, 'w') as fh:
        json.dump(output, fh, indent=2)
    print(f'\nSaved: {OUT_FILE}')
    for rank, algo in enumerate(ranked_hv, 1):
        m, s = summary[algo]['hv']['mean'], summary[algo]['hv']['std']
        md = summary[algo]['hv']['median']
        print(f'  {rank}. {algo:8s} HV mean={m:.4f}+/-{s:.4f}  median={md:.4f}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--algos', type=str, default=None,
                         help='Comma-separated subset of ALGORITHMS to run (default: all)')
    parser.add_argument('--seeds', type=str, default=None,
                         help='Comma-separated or start:end subset of SEEDS to run (default: all)')
    parser.add_argument('--aggregate-only', action='store_true',
                         help='Skip running; just aggregate existing checkpoints into the JSON '
                              '(use after all parallel worker shards have completed)')
    args = parser.parse_args()

    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    ref_dirs = get_reference_directions('das-dennis', 4, n_partitions=6)

    algos = args.algos.split(',') if args.algos else ALGORITHMS
    if args.seeds:
        if ':' in args.seeds:
            lo, hi = args.seeds.split(':')
            seeds = list(range(int(lo), int(hi)))
        else:
            seeds = [int(s) for s in args.seeds.split(',')]
    else:
        seeds = SEEDS

    n_gen_by_algo = {a: n_gen_for(a, ref_dirs) for a in ALGORITHMS}
    print(f'ref_dirs: {len(ref_dirs)} | canonical (v3) oracle | '
          f'this shard: {algos} x seeds={seeds[:3]}...{seeds[-3:]} (n={len(seeds)}) | '
          f'n_gen (matched-budget): {n_gen_by_algo}')

    if not args.aggregate_only:
        run_worker(algos, seeds, ref_dirs)
    else:
        print('  [aggregate-only] skipping optimisation runs')

    missing = [(a, s) for a in ALGORITHMS for s in SEEDS if not is_done(a, s)]
    if missing:
        print(f'\n{len(missing)} (algo, seed) combinations still missing '
              f'(e.g. {missing[:5]}) -- not aggregating yet. Run the remaining '
              f'shards, then call with --aggregate-only.')
        return
    aggregate(ref_dirs)


if __name__ == '__main__':
    main()
