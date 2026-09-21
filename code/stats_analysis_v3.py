#!/usr/bin/env python3
"""
stats_analysis_v3.py — Global Holm correction + rank-biserial for paired data
==================================================================================
Two statistical corrections requested by the second review round:
1. Holm-Bonferroni correction applied JOINTLY across all tests (3 metrics
   x C(k,2) algorithm pairs, k=5 algorithms as of round 6 -> 30 tests, was
   9 with the original 3 algorithms), not per-metric-family separately as
   in v2's script.
2. Matched-pairs rank-biserial correlation reported alongside (not instead
   of) Cliff's delta for the paired (by-seed) Wilcoxon tests -- rank-biserial
   is the standard paired-data effect size; Cliff's delta is more typically
   used for independent samples.

Round 6 (external review, item 9): a review questioned whether pairing by
seed INDEX is valid across algorithms with different population sizes and
architectures. It is: pairing by replicate index (not by literal shared
randomness) is the standard design in the EC statistics literature (Derrac,
Garcia & Herrera 2011, the reference both external reviews this round
themselves cite) for comparing algorithms across independent runs, and is
MORE powerful than treating the runs as independent -- switching to
Kruskal-Wallis/Mann-Whitney as primary would be a regression, not a fix.
We instead (a) increased seeds from 10 to 30, which was review round 6's
actual statistical concern (with 10 seeds, the paired two-sided Wilcoxon's
own minimum-achievable p-value, 2/2^10=0.00195, sits uncomfortably close to
several of this paper's previously-reported "significant" p-values, making
it impossible to tell a real large effect from the test's own resolution
floor; 30 seeds moves that floor to 2/2^30~1.86e-9), and (b) added an
independent-samples sensitivity check (Kruskal-Wallis + Mann-Whitney,
below) alongside the primary paired test, so a skeptical reader can confirm
the qualitative ranking does not depend on the pairing assumption.

Also runs the ablation's statistics. FIX (Round 3 review): the ablation's
8 variants share the SAME seeds (42-71 as of the round 6-followup 30-seed
expansion, applied for the same resolution-floor reason as the main
comparison above) -- every variant at seed s uses
the same RNG seed, so these are PAIRED, not independent, samples. The
previous version used Kruskal-Wallis (an independent-samples omnibus test)
+ Mann-Whitney U (independent-samples pairwise), which discards the pairing
structure and is the wrong test family for this design. Replaced with
Friedman (paired omnibus across all 8 variants, blocking on seed) +
Wilcoxon signed-rank vs the full decoder (paired pairwise) + Holm
correction + matched-pairs rank-biserial, mirroring exactly the method
already used above for the main algorithm comparison.
"""
import json
import numpy as np
from pathlib import Path
from scipy.stats import friedmanchisquare, wilcoxon, rankdata, kruskal, mannwhitneyu
from itertools import combinations

RESULTS = Path(__file__).resolve().parent.parent / 'results'


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    more = sum((x > y) for x in a for y in b)
    less = sum((x < y) for x in a for y in b)
    return (more - less) / (n1 * n2)


def rank_biserial_paired(a, b):
    """Matched-pairs rank-biserial correlation for the Wilcoxon signed-rank
    test: r = (sum of ranks for positive differences - sum for negative)
    / total sum of ranks. Range [-1, 1].

    FIX (Round 5 review): the previous rank computation
    (`np.argsort(np.argsort(x))+1`) assigns arbitrary, distinct integer
    ranks to tied values (breaking ties by array position) instead of the
    standard averaged rank Wilcoxon's own tie-correction assumes. Switched
    to `scipy.stats.rankdata(method='average')`. Verified this made no
    difference to any reported value in this manuscript (no exact ties
    occur in the HV/IGD/spacing differences actually computed here, since
    they are continuous floating-point outputs of independent stochastic
    runs), but is the methodologically correct implementation regardless
    of whether it happened to matter for this particular dataset."""
    diffs = np.asarray(a) - np.asarray(b)
    diffs = diffs[diffs != 0]
    if len(diffs) == 0:
        return 0.0
    abs_ranks = rankdata(np.abs(diffs), method='average')
    pos_sum = abs_ranks[diffs > 0].sum()
    neg_sum = abs_ranks[diffs < 0].sum()
    total = pos_sum + neg_sum
    return float((pos_sum - neg_sum) / total) if total > 0 else 0.0


def holm_correct(pvals):
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = min((m - rank) * pvals[idx], 1.0)
        running_max = max(running_max, val)
        adj[idx] = running_max
    return adj


def main():
    with open(RESULTS / 'baselines_comparison_v3.json') as f:
        base = json.load(f)

    algos = base['config']['algorithms']
    per_seed = base['per_run_details']
    seeds = base['config']['seeds']
    metrics = ['hv', 'igd', 'spacing']

    data = {m: {a: [] for a in algos} for m in metrics}
    for a in algos:
        recs = {r['seed']: r for r in per_seed[a]}
        for s in seeds:
            r = recs[s]
            for m in metrics:
                data[m][a].append(r[m] if r[m] is not None else np.nan)

    print('=' * 70)
    print(f'FRIEDMAN OMNIBUS TEST (paired by seed, {len(algos)} algorithms)')
    print('=' * 70)
    friedman_results = {}
    for m in metrics:
        arrs = [np.array(data[m][a]) for a in algos]
        stat, p = friedmanchisquare(*arrs)
        friedman_results[m] = {'statistic': float(stat), 'p': float(p)}
        print(f'  {m:10s}: chi2={stat:.3f}, p={p:.5f}')

    # ── Pairwise Wilcoxon, RAW p-values collected across ALL 9 tests ───────
    pairs = list(combinations(algos, 2))
    all_raw_p = []
    all_labels = []
    all_deltas = []
    all_rb = []
    for m in metrics:
        for a1, a2 in pairs:
            x, y = np.array(data[m][a1]), np.array(data[m][a2])
            try:
                stat, p = wilcoxon(x, y)
            except ValueError:
                p = 1.0
            all_raw_p.append(p)
            all_labels.append((m, f'{a1} vs {a2}'))
            all_deltas.append(cliffs_delta(x, y))
            all_rb.append(rank_biserial_paired(x, y))

    adj_p_global = holm_correct(np.array(all_raw_p))

    print()
    print('=' * 70)
    print(f'HOLM CORRECTION APPLIED JOINTLY ACROSS ALL {len(all_raw_p)} TESTS')
    print('=' * 70)
    pairwise_results = {m: {} for m in metrics}
    for (m, lbl), rp, ap, d, rb in zip(all_labels, all_raw_p, adj_p_global, all_deltas, all_rb):
        pairwise_results[m][lbl] = {
            'p_raw': float(rp), 'p_holm_global': float(ap),
            'cliffs_delta': float(d), 'rank_biserial': float(rb),
        }
        print(f'  {m:10s} {lbl:25s} p_raw={rp:.4f} p_holm({len(all_raw_p)})={ap:.4f} '
              f'delta={d:+.3f} r_rb={rb:+.3f}')

    hv_means = {a: np.mean(data['hv'][a]) for a in algos}
    ranked = sorted(algos, key=lambda a: hv_means[a], reverse=True)
    print(f'\nHV ranking: {" > ".join(ranked)}  (means: '
          f'{", ".join(f"{a}={hv_means[a]:.4f}" for a in ranked)})')

    # ── Independent-samples robustness check (external review round 6,
    #    item 9): the primary test above pairs runs by seed INDEX across
    #    algorithms with different population sizes/architectures -- a
    #    standard design in the EC statistics literature (Derrac, Garcia,
    #    Herrera 2011) where "seed index" is the blocking factor (matched
    #    replicate number), not a claim that algorithms share literal
    #    random state. We keep the paired Friedman/Wilcoxon design above as
    #    PRIMARY (it is the more powerful, standard choice for this design,
    #    and switching to an independent-samples test would be a regression,
    #    not an improvement) but additionally report Kruskal-Wallis +
    #    Mann-Whitney U (treating the same seed's per-algorithm runs as
    #    independent, discarding the pairing) as a sensitivity check: if the
    #    qualitative ranking and significance pattern hold under BOTH
    #    designs, the conclusion does not depend on the pairing assumption.
    print()
    print('=' * 70)
    print('SENSITIVITY CHECK: independent-samples tests (same data, no pairing)')
    print('=' * 70)
    independent_out = {'kruskal': {}, 'pairwise_mannwhitney_holm_global': {}}
    all_raw_p_ind, all_labels_ind = [], []
    for m in metrics:
        arrs = [np.array(data[m][a]) for a in algos]
        stat, p = kruskal(*arrs)
        independent_out['kruskal'][m] = {'statistic': float(stat), 'p': float(p)}
        print(f'  Kruskal-Wallis {m:10s}: H={stat:.3f}, p={p:.5f}')
        for a1, a2 in pairs:
            x, y = np.array(data[m][a1]), np.array(data[m][a2])
            try:
                stat_mw, p_mw = mannwhitneyu(x, y, alternative='two-sided')
            except ValueError:
                p_mw = 1.0
            all_raw_p_ind.append(p_mw)
            all_labels_ind.append((m, f'{a1} vs {a2}', cliffs_delta(x, y)))
    adj_p_ind = holm_correct(np.array(all_raw_p_ind))
    for (m, lbl, delta), rp, ap in zip(all_labels_ind, all_raw_p_ind, adj_p_ind):
        # Vargha-Delaney A12 = (delta+1)/2, reported alongside Cliff's delta
        # (external review round 6-followup: an independent-samples pairwise
        # test needs its own effect size, analogous to the paired design's
        # matched-pairs rank-biserial correlation).
        independent_out['pairwise_mannwhitney_holm_global'].setdefault(m, {})[lbl] = {
            'p_raw': float(rp), 'p_holm_global': float(ap),
            'cliffs_delta': float(delta), 'vargha_delaney_a12': float((delta + 1) / 2)}
        print(f'  Mann-Whitney {m:10s} {lbl:25s} p_raw={rp:.4f} p_holm={ap:.4f} delta={delta:+.3f}')

    # ── Ablation stats: Friedman (paired omnibus) + Wilcoxon signed-rank
    #    vs full decoder (paired pairwise) + Holm + rank-biserial. FIX
    #    (Round 3 review): the 8 variants share the same seeds, so this
    #    is paired data -- Kruskal-Wallis/Mann-Whitney (independent-samples
    #    tests) discarded that structure. See module docstring.
    print()
    print('=' * 70)
    print('ABLATION: FRIEDMAN OMNIBUS (paired by seed) + PAIRED WILCOXON vs FULL (HV)')
    print('=' * 70)
    ablation_out = {}
    abl_path = RESULTS / 'ablation_v3_results.json'
    if abl_path.exists():
        with open(abl_path) as f:
            abl = json.load(f)
        variants = abl['config']['variants']
        hv_raw = {v: abl['summary'][v]['hv_raw'] for v in variants}
        n_seeds_per_variant = {v: len(hv_raw[v]) for v in variants}
        if len(set(n_seeds_per_variant.values())) != 1:
            raise RuntimeError(
                f'Ablation variants do not all have the same number of valid '
                f'seed replicates ({n_seeds_per_variant}) -- the Friedman/'
                f'Wilcoxon tests below require perfectly paired (same-seed) '
                f'data across all variants. Fix the failing seed(s) before '
                f'trusting these statistics.')

        arrs = [np.array(hv_raw[v]) for v in variants]
        stat, p = friedmanchisquare(*arrs)
        print(f'  Friedman (all {len(variants)} variants, n={n_seeds_per_variant[variants[0]]} '
              f'paired seeds): chi2={stat:.3f}, p={p:.6f}')
        ablation_out['friedman'] = {'statistic': float(stat), 'p': float(p)}

        ref = 'A_full_decoder'
        others = [v for v in variants if v != ref]
        raw_p, deltas, rbs = [], [], []
        for v in others:
            try:
                stat_w, p_w = wilcoxon(hv_raw[ref], hv_raw[v])
            except ValueError:
                p_w = 1.0
            raw_p.append(p_w)
            deltas.append(cliffs_delta(hv_raw[ref], hv_raw[v]))
            rbs.append(rank_biserial_paired(hv_raw[ref], hv_raw[v]))
        adj_p = holm_correct(np.array(raw_p))
        ablation_out['pairwise_vs_full'] = {}
        for v, rp, ap, d, rb in zip(others, raw_p, adj_p, deltas, rbs):
            ablation_out['pairwise_vs_full'][v] = {
                'p_raw': float(rp), 'p_holm': float(ap),
                'cliffs_delta': float(d), 'rank_biserial': float(rb),
            }
            sig = 'DIFFERS' if ap < 0.05 else 'no sig. diff.'
            print(f'  A_full vs {v:20s} p_raw={rp:.4f} p_holm={ap:.4f} '
                  f'delta={d:+.3f} r_rb={rb:+.3f}  [{sig}]')
    else:
        print('  (ablation_v3_results.json not yet available)')

    output = {
        'friedman': friedman_results,
        'pairwise_wilcoxon_holm_global': pairwise_results,
        'hv_ranking': ranked,
        'hv_means': hv_means,
        'ablation': ablation_out,
        'independent_samples_sensitivity_check': independent_out,
    }
    with open(RESULTS / 'stats_analysis_v3.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f'\nSaved: {RESULTS / "stats_analysis_v3.json"}')


if __name__ == '__main__':
    main()
