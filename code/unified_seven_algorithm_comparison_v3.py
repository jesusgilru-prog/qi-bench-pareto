#!/usr/bin/env python3
"""
unified_seven_algorithm_comparison_v3.py — one comparison, seven
algorithms, one reference front.

External review raised two connected objections to the way SPEA2 and
AGE-MOEA-II were reported:

  (#6) they were added as a side paragraph rather than integrated, so
       the abstract, the main table and the conclusion still presented
       NSGA-II as the sole winner among five, with no IGD, spacing,
       cardinality or runtime for the two extra algorithms and no joint
       family of corrected tests;
  (#7) "same normalisation convention" does not mean "same numbers":
       run_modern_algorithms_v3.py normalises against the union of the
       150 main per-seed fronts PLUS its own 60, which is a strictly
       wider bounding box than Table 3's 150-front union, so the two
       sets of HV values were not on identical scales despite the
       claim that they were.

This script settles both by re-scoring ALL SEVEN algorithms against a
single reference: bounds and reference front are built once from the
union of all 210 raw per-seed fronts (5 main algorithms x 30 seeds +
SPEA2 and AGE-MOEA-II x 30 seeds). Every reported HV, IGD and spacing
value therefore lives on exactly the same scale. It also runs:

  * a joint Kruskal-Wallis omnibus per metric across all seven;
  * Holm-corrected pairwise Mann-Whitney over the whole family of
    3 metrics x C(7,2)=21 pairs = 63 tests (independent samples,
    the manuscript's primary analysis);
  * a paired Friedman + Wilcoxon sensitivity analysis blocked on seed,
    Holm-corrected over the same family size.

Nothing is re-optimised: this is a re-scoring of already-committed
checkpoints, so it adds no new stochastic results.
"""
import json
import pickle
from itertools import combinations
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from pymoo.indicators.igd import IGD
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from scipy.stats import friedmanchisquare, kruskal, mannwhitneyu, rankdata, wilcoxon

from metrics import compute_spacing

RESULTS = Path(__file__).resolve().parent.parent / 'results'
MAIN_CKPT = RESULTS / 'step1_baselines_v3'
MODERN_CKPT = RESULTS / 'step9_modern_algorithms_v3'
OUT_JSON = RESULTS / 'unified_seven_algorithm_comparison_v3.json'

MAIN_ALGOS = ['nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea']
MODERN_ALGOS = ['spea2', 'age_moea2']
ALGOS = MAIN_ALGOS + MODERN_ALGOS
SEEDS = list(range(42, 72))
HV_REF = np.array([1.1] * 4)
METRICS = ['hv', 'igd', 'spacing']
HIGHER_IS_BETTER = {'hv': True, 'igd': False, 'spacing': False}


def ckpt_dir(algo):
    return MODERN_CKPT if algo in MODERN_ALGOS else MAIN_CKPT


def load(algo, seed):
    with open(ckpt_dir(algo) / f'{algo}_seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)


def cliffs_delta(a, b):
    a = np.asarray(a, dtype=float)[:, None]
    b = np.asarray(b, dtype=float)[None, :]
    return float(((a > b).sum() - (a < b).sum()) / (a.size * b.size))


def rank_biserial_paired(a, b):
    diffs = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    diffs = diffs[diffs != 0]
    if len(diffs) == 0:
        return 0.0
    r = rankdata(np.abs(diffs), method='average')
    pos, neg = r[diffs > 0].sum(), r[diffs < 0].sum()
    total = pos + neg
    return float((pos - neg) / total) if total > 0 else 0.0


def holm(pvals):
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for i, idx in enumerate(order):
        running = max(running, min(1.0, (m - i) * pvals[idx]))
        adj[idx] = running
    return adj


def main():
    runs = {a: [load(a, s) for s in SEEDS] for a in ALGOS}

    # ---- one reference for all seven -------------------------------------
    all_F = np.vstack([r['F_pareto'] for a in ALGOS for r in runs[a]
                       if len(r['F_pareto']) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def normalise(F):
        return (F - g_min) / g_range

    nd = find_non_dominated(all_F)
    ref_front = normalise(np.unique(np.round(all_F[nd], 6), axis=0))
    hv_calc, igd_calc = Hypervolume(ref_point=HV_REF), IGD(ref_front)

    per_seed = {a: {m: [] for m in METRICS} for a in ALGOS}
    for a in ALGOS:
        per_seed[a]['n_pareto'] = []
        per_seed[a]['elapsed_min'] = []
        for r in runs[a]:
            F = r['F_pareto']
            if len(F) == 0:
                per_seed[a]['hv'].append(0.0)
                per_seed[a]['igd'].append(float('nan'))
                per_seed[a]['spacing'].append(float('nan'))
                per_seed[a]['n_pareto'].append(0)
                per_seed[a]['elapsed_min'].append(float(r['elapsed_min']))
                continue
            Fn = normalise(F)
            per_seed[a]['hv'].append(float(hv_calc(Fn)))
            per_seed[a]['igd'].append(float(igd_calc(Fn)))
            per_seed[a]['spacing'].append(float(compute_spacing(Fn)))
            per_seed[a]['n_pareto'].append(int(len(F)))
            per_seed[a]['elapsed_min'].append(float(r['elapsed_min']))

    def ms(v):
        arr = np.asarray(v, dtype=float)
        return {'mean': float(np.mean(arr)),
                'std': float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
                'median': float(np.median(arr)),
                'iqr': [float(np.percentile(arr, 25)), float(np.percentile(arr, 75))]}

    summary = {a: {k: ms(per_seed[a][k])
                   for k in METRICS + ['n_pareto', 'elapsed_min']} for a in ALGOS}
    for a in ALGOS:
        # Per-seed values for every metric, so figures and any downstream
        # re-analysis are driven by this artefact rather than by a
        # differently-normalised one.
        for k in METRICS + ['n_pareto', 'elapsed_min']:
            summary[a][f'{k}_raw'] = per_seed[a][k]

    # ---- omnibus ---------------------------------------------------------
    omnibus = {}
    for m in METRICS:
        H, p = kruskal(*[per_seed[a][m] for a in ALGOS])
        omnibus[f'kruskal_{m}'] = {'H': float(H), 'p': float(p), 'k': len(ALGOS)}
    for m in METRICS:
        chi2, p = friedmanchisquare(*[per_seed[a][m] for a in ALGOS])
        omnibus[f'friedman_{m}'] = {'chi2': float(chi2), 'p': float(p), 'k': len(ALGOS)}

    # ---- pairwise, one joint Holm family ---------------------------------
    keys, raw_p, records = [], [], []
    for m in METRICS:
        for a, b in combinations(ALGOS, 2):
            u, p = mannwhitneyu(per_seed[a][m], per_seed[b][m], alternative='two-sided')
            keys.append(f'{m}:{a}_vs_{b}')
            raw_p.append(float(p))
            records.append({'metric': m, 'a': a, 'b': b, 'u_stat': float(u),
                            'p_raw': float(p),
                            'cliffs_delta': cliffs_delta(per_seed[a][m], per_seed[b][m])})
    adj = holm(raw_p)
    unpaired = {}
    for k, rec, ap in zip(keys, records, adj):
        rec['p_holm'] = float(ap)
        rec['significant_holm_005'] = bool(ap < 0.05)
        unpaired[k] = rec

    keys_p, raw_pp, records_p = [], [], []
    for m in METRICS:
        for a, b in combinations(ALGOS, 2):
            stat, p = wilcoxon(per_seed[a][m], per_seed[b][m])
            keys_p.append(f'{m}:{a}_vs_{b}')
            raw_pp.append(float(p))
            records_p.append({'metric': m, 'a': a, 'b': b, 'w_stat': float(stat),
                              'p_raw': float(p),
                              'rank_biserial': rank_biserial_paired(per_seed[a][m],
                                                                   per_seed[b][m])})
    adj_p = holm(raw_pp)
    paired = {}
    for k, rec, ap in zip(keys_p, records_p, adj_p):
        rec['p_holm'] = float(ap)
        rec['significant_holm_005'] = bool(ap < 0.05)
        paired[k] = rec

    # ---- top group by HV -------------------------------------------------
    order = sorted(ALGOS, key=lambda a: -summary[a]['hv']['mean'])
    best = order[0]
    tied_with_best = [a for a in ALGOS
                      if a == best
                      or not unpaired[f'hv:{min(a, best, key=ALGOS.index)}_vs_'
                                      f'{max(a, best, key=ALGOS.index)}']['significant_holm_005']]

    out = {
        'description': ('Unified comparison of all SEVEN algorithms (NSGA-II, NSGA-III, '
                        'MOEA/D, SMS-EMOA, RVEA, SPEA2, AGE-MOEA-II) re-scored against a '
                        'single reference front and single min-max normalisation built '
                        'from the union of all 210 raw per-seed fronts, so every HV, IGD '
                        'and spacing value is on the same scale. Supersedes the split '
                        'presentation in which SPEA2/AGE-MOEA-II were scored against a '
                        'wider bounding box than Table 3 while being described as '
                        'directly comparable to it. No new optimisation runs: this is a '
                        're-scoring of committed checkpoints.'),
        'config': {'algorithms': ALGOS, 'seeds': SEEDS, 'hv_ref_point': HV_REF.tolist(),
                   'n_fronts_pooled': int(len(ALGOS) * len(SEEDS)),
                   'holm_family_size': len(raw_p)},
        'normalisation': {'g_min': g_min.tolist(), 'g_max': g_max.tolist(),
                          'ref_front_size': int(ref_front.shape[0])},
        'summary': summary,
        'ranking_by_mean_hv': order,
        'top_group_indistinguishable_from_best_hv': sorted(tied_with_best),
        'omnibus': omnibus,
        'pairwise_unpaired_mannwhitney_holm': unpaired,
        'pairwise_paired_wilcoxon_holm': paired,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)

    print(f'Saved: {OUT_JSON}')
    print(f'{"algorithm":12s} {"HV":>18s} {"IGD":>18s} {"spacing":>18s} {"|F*|":>12s} {"min":>8s}')
    for a in order:
        s = summary[a]
        print(f'{a:12s} {s["hv"]["mean"]:.3f}+-{s["hv"]["std"]:.3f}      '
              f'{s["igd"]["mean"]:.3f}+-{s["igd"]["std"]:.3f}      '
              f'{s["spacing"]["mean"]:.3f}+-{s["spacing"]["std"]:.3f}      '
              f'{s["n_pareto"]["mean"]:6.1f}   {s["elapsed_min"]["mean"]:6.2f}')
    print('\nomnibus:', {k: round(v['p'], 6) for k, v in omnibus.items()})
    print('\nHV pairs significant after joint Holm:')
    for k, v in unpaired.items():
        if v['metric'] == 'hv' and v['significant_holm_005']:
            print(f'  {k}: p_holm={v["p_holm"]:.4g} delta={v["cliffs_delta"]:+.2f}')
    print('\ntop group indistinguishable from best on HV:',
          out['top_group_indistinguishable_from_best_hv'])


if __name__ == '__main__':
    main()
