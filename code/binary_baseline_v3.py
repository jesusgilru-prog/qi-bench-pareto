#!/usr/bin/env python3
"""
binary_baseline_v3.py — True 4-way categorical binary/medium baseline
==========================================================================
Fixes a fairness gap the reviewer correctly identified in v2: assigning
each active edge its argmax-fidelity medium (BEST_MEDIUM) is NOT the same
representational freedom as the parametric decoder, because the decoder
can choose a lower-fidelity medium when it trades off better on rate,
latency or coverage. Here each of the 4,950 candidate city pairs is a
genuinely free categorical choice in {no edge, fibre, satellite, FSO},
searched via a real-valued relaxation (each variable in [-0.499, 3.499],
rounded to the nearest integer category at evaluation time) -- a standard
pragmatic encoding for categorical variables under continuous operators
(SBX/PM), used here instead of pymoo's stricter integer-native operators
for robustness. This is documented in the manuscript as a modelling
choice, not hidden.

Also reports HV/IGD/spacing against the SAME reference front used for
the main algorithm comparison (not just |F*|, which the reviewer
correctly noted does not by itself indicate front quality), and uses
bdcz_oracle_v3 (canonical, exact Werner-optimal oracle).
"""
import numpy as np
from metrics import compute_spacing
import json
import time
from pathlib import Path

from bdcz_oracle_v3 import oracle_evaluate, N_NODES

from pymoo.core.problem import Problem
from pymoo.core.sampling import Sampling
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD

N_EDGES = N_NODES * (N_NODES - 1) // 2  # 4950
SEEDS = list(range(42, 72))  # FIX (external review round 6-followup, confirmed
                              # real): the manuscript claimed 10 seeds here
                              # "match the statistical power of the main
                              # comparison", which became false once the main
                              # comparison was raised to 30 seeds in Round 6.
                              # Raised to 30 (42-71) to make that claim true
                              # again, rather than just rewording it. Seed
                              # history: 3 (orig.) -> 10 (Round 4) -> 30 (this).
POP_SIZE, N_GEN = 104, 80
HV_REF = np.array([1.1] * 4)

CKPT_DIR = Path(__file__).resolve().parent.parent / 'results' / 'step2_binary_v3'

EDGE_MAP = np.zeros((N_EDGES, 2), dtype=np.int64)
idx = 0
for i in range(N_NODES):
    for j in range(i + 1, N_NODES):
        EDGE_MAP[idx] = [i, j]
        idx += 1


class SparseCategoricalSampling(Sampling):
    """~3% of edges active initially (category != 0), medium drawn uniformly
    among {fibre, satellite, FSO} for those active edges."""
    def _do(self, problem, n_samples, **kwargs):
        X = np.zeros((n_samples, problem.n_var))
        active = np.random.random((n_samples, problem.n_var)) < 0.03
        cat = np.random.randint(1, 4, size=(n_samples, problem.n_var))
        X[active] = cat[active]
        # jitter within the category's real-valued band so SBX has room to work
        X += np.random.uniform(-0.3, 0.3, size=X.shape)
        return np.clip(X, -0.499, 3.499)


class CategoricalTopologyProblem(Problem):
    def __init__(self):
        super().__init__(n_var=N_EDGES, n_obj=4, n_constr=0, xl=-0.499, xu=3.499)

    def _evaluate(self, X, out, *args, **kwargs):
        cats = np.rint(X).astype(np.int64)  # nearest-integer rounding -> {0,1,2,3}
        F = np.zeros((len(X), 4))
        for row_idx, row in enumerate(cats):
            active_idx = np.where(row > 0)[0]
            if len(active_idx) < 5:
                F[row_idx] = [0.0, 0.0, 1e6, 0.0]
                continue
            eu = EDGE_MAP[active_idx, 0].copy()
            ev = EDGE_MAP[active_idx, 1].copy()
            et = (row[active_idx] - 1).astype(np.int64)  # 1,2,3 -> 0,1,2 (fibre,sat,fso)
            try:
                fid, rate_log, lat, cov = oracle_evaluate(eu, ev, et)
                F[row_idx] = [-fid, -rate_log, lat, -cov]
            except Exception:
                F[row_idx] = [0.0, 0.0, 1e6, 0.0]
        out["F"] = F




def ckpt_path(algo_name, seed): return CKPT_DIR / f'{algo_name}_seed{seed}.pkl'
def done_path(algo_name, seed): return CKPT_DIR / f'{algo_name}_seed{seed}_completed.txt'
def is_done(algo_name, seed): return done_path(algo_name, seed).exists()


def run_worker(algo_names, seeds):
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    for algo_name in algo_names:
        for seed in seeds:
            if is_done(algo_name, seed):
                print(f'  [SKIP] {algo_name} seed={seed}')
                continue
            print(f"\n{'='*60}\nRunning {algo_name} seed={seed}...\n{'='*60}")
            np.random.seed(seed)
            t0 = time.time()

            if 'nsga3' in algo_name:
                ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=6)
                algo = NSGA3(ref_dirs=ref_dirs, pop_size=POP_SIZE,
                            sampling=SparseCategoricalSampling(),
                            crossover=SBX(prob=0.9, eta=15),
                            mutation=PM(eta=20))
            else:
                algo = NSGA2(pop_size=POP_SIZE, sampling=SparseCategoricalSampling(),
                            crossover=SBX(prob=0.9, eta=15),
                            mutation=PM(eta=20))

            problem = CategoricalTopologyProblem()
            res = minimize(problem, algo, ('n_gen', N_GEN), seed=seed, verbose=True)
            elapsed = time.time() - t0

            if res.F is not None and len(res.F) > 0:
                nd_mask = find_non_dominated(res.F)
                F_par = res.F[nd_mask]
                # FIX (Round 4 review, same confirmed issue as the decoder
                # experiments): dedup by objective vector before reporting
                # |F*|/HV/IGD/spacing -- see dedup_front's docstring in
                # bdcz_oracle_v3.py for why the raw non-dominated population
                # routinely contains a large fraction of exact duplicates.
                F_par = np.unique(np.round(F_par, 6), axis=0)
            else:
                F_par = np.empty((0, 4))

            result_entry = {
                'seed': seed, 'n_pareto': len(F_par),
                'elapsed_min': elapsed / 60,
                'pareto_front': F_par.tolist(),
            }
            with open(ckpt_path(algo_name, seed), 'w') as fh:
                json.dump(result_entry, fh)
            with open(done_path(algo_name, seed), 'w') as fh:
                fh.write('done\n')
            print(f"\n  Seed {seed}: |F*|={len(F_par)}, time={elapsed/60:.1f} min")


def aggregate():
    algo_names = ['nsga2_categorical', 'nsga3_categorical']
    missing = [(a, s) for a in algo_names for s in SEEDS if not is_done(a, s)]
    if missing:
        print(f'Cannot aggregate: {len(missing)} (algo, seed) combos missing, e.g. {missing[:5]}')
        return

    results = {'nsga2_categorical': [], 'nsga3_categorical': []}
    all_fronts = {}
    for algo_name in algo_names:
        all_fronts[algo_name] = []
        for seed in SEEDS:
            with open(ckpt_path(algo_name, seed)) as fh:
                entry = json.load(fh)
            results[algo_name].append(entry)
            all_fronts[algo_name].append(np.array(entry['pareto_front']).reshape(-1, 4))

    _finish(results, all_fronts)


def _finish(results, all_fronts):
    # ── HV/IGD/spacing against a reference front shared ONLY between the two
    # categorical conditions (nsga2_categorical, nsga3_categorical) below --
    # NOT shared with the decoder. FIX (external review round 6-followup): this
    # block's output key was previously named 'metrics_vs_shared_reference',
    # which reads as if it were the decoder-vs-categorical shared reference
    # computed in build_binary_vs_decoder_shared_ref_v3.py (the numbers the
    # manuscript actually cites) -- it is a different, self-contained
    # computation and was never referenced by the manuscript, but the
    # misleading name could make a reader think the two files disagree on the
    # same quantity. Renamed for clarity; values are unchanged.
    combined = np.vstack([F for fs in all_fronts.values() for F in fs if len(F) > 0])
    if len(combined) > 0:
        g_min, g_max = combined.min(axis=0), combined.max(axis=0)
        g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)
        def normalise(F): return (F - g_min) / g_range
        nd_global = find_non_dominated(combined)
        ref_front = normalise(np.unique(np.round(combined[nd_global], 6), axis=0))
        hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

        metrics_summary = {}
        for algo_name, fronts in all_fronts.items():
            hv_list, igd_list, sp_list = [], [], []
            for F in fronts:
                if len(F) == 0:
                    continue
                F_norm = normalise(F)
                hv_list.append(float(hv_calc(F_norm)))
                igd_list.append(float(igd_calc(F_norm)))
                sp_list.append(compute_spacing(F_norm))
            def ms(lst):
                a = np.array(lst)
                std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
                return {'mean': float(a.mean()), 'std': std} if len(a) else {'mean': None, 'std': None}
            metrics_summary[algo_name] = {'hv': ms(hv_list), 'igd': ms(igd_list), 'spacing': ms(sp_list)}
        results['metrics_vs_categorical_only_reference'] = metrics_summary
        print("\n=== HV/IGD/spacing vs categorical-only reference front "
              "(NOT the decoder-shared reference -- see "
              "build_binary_vs_decoder_shared_ref_v3.py for that) ===")
        print(json.dumps(metrics_summary, indent=2))

    out_path = Path(__file__).resolve().parent.parent / 'results' / 'binary_baseline_v3_results.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


def parse_seeds(spec):
    if ':' in spec:
        a, b = spec.split(':')
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(',')]


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--algos', type=str, default=None,
                    help='comma-separated subset of nsga2_categorical,nsga3_categorical')
    p.add_argument('--seeds', type=str, default=None,
                    help='comma-separated seeds or start:end range, default all SEEDS')
    p.add_argument('--aggregate-only', action='store_true')
    args = p.parse_args()

    algo_names = args.algos.split(',') if args.algos else ['nsga2_categorical', 'nsga3_categorical']
    seeds = parse_seeds(args.seeds) if args.seeds else SEEDS

    if not args.aggregate_only:
        run_worker(algo_names, seeds)
    else:
        aggregate()


if __name__ == '__main__':
    main()
