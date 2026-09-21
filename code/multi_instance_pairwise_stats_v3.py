#!/usr/bin/env python3
"""
multi_instance_pairwise_stats_v3.py — NSGA-II vs. NSGA-III pairwise HV
comparison per additional instance, using each instance's OWN
two-algorithm reference front/normalisation (not the instance's
all-five-algorithm reference used in multi_instance_v3_results.json's
'summary' block). The manuscript (Section: Multi-instance validation)
reports both this pairwise comparison and the all-five-algorithm
summary; this script exists so the pairwise numbers are traceable to a
committed artefact rather than only appearing in prose (regression:
these numbers were previously computed ad hoc and not persisted).

All FOUR additional instances are treated as one family of hypothesis
tests and Holm-corrected jointly (external review, finding #11: the
per-instance p-values were previously reported uncorrected, which
overstates the evidence for the smallest of them). Cliff's delta and a
percentile bootstrap CI on the mean HV difference are reported per
instance so that non-significant results can be read as effect-size
estimates rather than as evidence of no effect.

Run after `run_multi_instance_v3.py --aggregate-only` has produced
checkpoints in results/step7_multi_instance_v3/.
"""
import json
import pickle
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from scipy.stats import mannwhitneyu

RESULTS = Path(__file__).resolve().parent.parent / 'results'
CKPT = RESULTS / 'step7_multi_instance_v3'
OUT_JSON = RESULTS / 'multi_instance_pairwise_stats_v3.json'

# Todas las instancias usan ahora 10 semillas. global_n200 se corría antes
# con 5 por su coste; la recomputación posterior al fix satelital la extendió
# a 10, que era justo lo que el manuscrito listaba como trabajo futuro.
INSTANCES = {
    'global_n50': list(range(42, 52)),
    'regional_eu_na': list(range(42, 52)),
    'metro_europe': list(range(42, 52)),
    'global_n200': list(range(42, 52)),
}
HV_REF = np.array([1.1] * 4)
BOOT_N = 20000
BOOT_SEED = 20260730


def load(inst, algo, seed):
    with open(CKPT / f'{inst}_{algo}_seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)


def cliffs_delta(a, b):
    """Cliff's delta for a vs b: P(a>b) - P(a<b)."""
    a = np.asarray(a, dtype=float)[:, None]
    b = np.asarray(b, dtype=float)[None, :]
    return float(((a > b).sum() - (a < b).sum()) / (a.size * b.size))


def bootstrap_mean_diff_ci(a, b, n=BOOT_N, seed=BOOT_SEED, alpha=0.05):
    """Percentile bootstrap CI on mean(a) - mean(b) (independent samples)."""
    rng = np.random.default_rng(seed)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    diffs = np.empty(n)
    for i in range(n):
        diffs[i] = (rng.choice(a, a.size, replace=True).mean()
                    - rng.choice(b, b.size, replace=True).mean())
    lo, hi = np.percentile(diffs, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def holm(pvals):
    """Holm-Bonferroni step-down adjusted p-values, order preserved."""
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for i, idx in enumerate(order):
        running = max(running, min(1.0, (m - i) * pvals[idx]))
        adj[idx] = running
    return adj


def main():
    out = {}
    for inst, seeds in INSTANCES.items():
        runs = {a: [load(inst, a, s) for s in seeds] for a in ['nsga2', 'nsga3']}
        all_F = np.vstack([r['F_pareto'] for a in ['nsga2', 'nsga3']
                            for r in runs[a] if len(r['F_pareto']) > 0])
        g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
        g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

        def normalise(F):
            return (F - g_min) / g_range

        hv_calc = Hypervolume(ref_point=HV_REF)
        hv = {}
        for a in ['nsga2', 'nsga3']:
            hv[a] = [float(hv_calc(normalise(r['F_pareto']))) if len(r['F_pareto']) > 0 else 0.0
                      for r in runs[a]]
        u_stat, p = mannwhitneyu(hv['nsga2'], hv['nsga3'], alternative='two-sided')
        ci_lo, ci_hi = bootstrap_mean_diff_ci(hv['nsga2'], hv['nsga3'])
        out[inst] = {
            'nsga2_hv_mean': float(np.mean(hv['nsga2'])),
            'nsga2_hv_std': float(np.std(hv['nsga2'], ddof=1)),
            'nsga3_hv_mean': float(np.mean(hv['nsga3'])),
            'nsga3_hv_std': float(np.std(hv['nsga3'], ddof=1)),
            'mean_diff_nsga2_minus_nsga3': float(np.mean(hv['nsga2']) - np.mean(hv['nsga3'])),
            'mean_diff_bootstrap_ci95': [ci_lo, ci_hi],
            'cliffs_delta_nsga2_vs_nsga3': cliffs_delta(hv['nsga2'], hv['nsga3']),
            'mannwhitney_u': float(u_stat),
            'mannwhitney_p': float(p),
            'n_seeds': len(seeds),
            'hv_raw_nsga2': hv['nsga2'],
            'hv_raw_nsga3': hv['nsga3'],
        }

    names = list(out)
    adj = holm([out[k]['mannwhitney_p'] for k in names])
    for k, a in zip(names, adj):
        out[k]['mannwhitney_p_holm'] = float(a)
        out[k]['significant_holm_005'] = bool(a < 0.05)

    result = {
        'description': ('NSGA-II vs. NSGA-III pairwise HV comparison per additional '
                         'instance, using a reference front/normalisation built from '
                         'ONLY these two algorithms\' pooled checkpoints for that '
                         'instance (distinct from multi_instance_v3_results.json\'s '
                         'all-five-algorithm reference front). The four instances are '
                         'Holm-corrected jointly as one family; uncorrected p-values '
                         'are retained for transparency. Cliff\'s delta and a 20,000-'
                         'resample percentile bootstrap CI on the mean HV difference '
                         'are reported so non-significant instances can be read as '
                         'effect-size estimates rather than as evidence of no effect.'),
        'family_correction': 'Holm-Bonferroni over the 4 per-instance tests',
        'bootstrap': {'n_resamples': BOOT_N, 'seed': BOOT_SEED, 'level': 0.95},
        'per_instance': out,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    for inst, d in out.items():
        print(f"  {inst:15s} nsga2={d['nsga2_hv_mean']:.3f} nsga3={d['nsga3_hv_mean']:.3f} "
              f"p={d['mannwhitney_p']:.4g} p_holm={d['mannwhitney_p_holm']:.4g} "
              f"delta={d['cliffs_delta_nsga2_vs_nsga3']:+.2f} "
              f"CI95=[{d['mean_diff_bootstrap_ci95'][0]:+.4f},{d['mean_diff_bootstrap_ci95'][1]:+.4f}]")


if __name__ == '__main__':
    main()
