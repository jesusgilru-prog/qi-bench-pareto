#!/usr/bin/env python3
"""
generate_tables_v3.py — Emit final LaTeX table/paragraph content from v3 JSON
==================================================================================
Same purpose as v2's generator: copy numbers programmatically from the
actual canonical-oracle result files, never by hand, into a single .tex
snippet ready to paste into the manuscript.
"""
import json
import numpy as np
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT = RESULTS / 'GENERATED_TABLES_v3.tex'

lines = []


def add(s=''):
    lines.append(s)


with open(RESULTS / 'baselines_comparison_v3.json') as f:
    base = json.load(f)
algos = base['config']['algorithms']
n_seeds = len(base['config']['seeds'])
summ = base['summary']
label = {'nsga2': 'NSGA-II', 'nsga3': 'NSGA-III', 'moead': 'MOEA/D',
         'sms_emoa': 'SMS-EMOA', 'rvea': 'RVEA'}

add('%% === TABLE: algo-comparison (v3, canonical oracle) ===')
add(r'\begin{table}[tbp]')
add(r'\centering')
add(f'\\caption{{Algorithm comparison on the quantum network topology '
    f'optimisation problem (${n_seeds}$ seeds, mean$\\,\\pm\\,$std, canonical oracle).}}')
add(r'\label{tab:algo-comparison}')
add(r'\small')
add(r'\setlength{\tabcolsep}{4pt}')
add(r'\resizebox{\columnwidth}{!}{%')
add(r'\begin{tabular}{lccccc}')
add(r'\toprule')
add(r'Algorithm & HV $\uparrow$ & IGD $\downarrow$ & Spacing $\downarrow$ & $|F^\star|$ & Time (min) \\')
add(r'\midrule')
hv_means = {a: summ[a]['hv']['mean'] for a in algos}
best_algo = max(hv_means, key=hv_means.get)
for a in algos:
    s = summ[a]
    bold = a == best_algo
    def fmt(m, dec=3):
        v, sd = s[m]['mean'], s[m]['std']
        txt = f'{v:.{dec}f} \\pm {sd:.{dec}f}'
        return f'\\mathbf{{{txt}}}' if bold else txt
    npar = f"{s['n_pareto']['mean']:.1f} \\pm {s['n_pareto']['std']:.1f}"
    time_s = f"{s['elapsed_min']['mean']:.2f} \\pm {s['elapsed_min']['std']:.2f}"
    add(f"{label[a]:8s} & ${fmt('hv')}$ & ${fmt('igd')}$ & ${fmt('spacing')}$ & ${npar}$ & ${time_s}$ \\\\")
add(r'\bottomrule')
add(r'\end{tabular}%')
add(r'}')
add(r'\end{table}')
add()

stats_path = RESULTS / 'stats_analysis_v3.json'
if stats_path.exists():
    with open(stats_path) as f:
        stats = json.load(f)
    add('%% === TABLE: friedman ===')
    add(r'\begin{table}[tbp]')
    add(r'\centering')
    add(f'\\caption{{Friedman omnibus test across the {len(algos)} algorithms '
        f'(paired by seed, $n={n_seeds}$).}}')
    add(r'\label{tab:friedman}')
    add(r'\small')
    add(r'\begin{tabular}{lcc}')
    add(r'\toprule')
    add(r'Metric & $\chi^2$ & $p$ \\')
    add(r'\midrule')
    for m, key in [('HV', 'hv'), ('IGD', 'igd'), ('Spacing', 'spacing')]:
        fr = stats['friedman'][key]
        add(f"{m} & {fr['statistic']:.3f} & {fr['p']:.5f} \\\\")
    add(r'\bottomrule')
    add(r'\end{tabular}')
    add(r'\end{table}')
    add()

    pw_n_tests = sum(len(v) for v in stats['pairwise_wilcoxon_holm_global'].values())
    n_pairs = len(list(stats['pairwise_wilcoxon_holm_global']['hv'].keys()))
    add(f'%% === TABLE: holm-wilcoxon (corrected jointly across all {pw_n_tests} tests) ===')
    add(r'\begin{table}[tbp]')
    add(r'\centering')
    add(f"\\caption{{Pairwise Wilcoxon signed-rank tests ($n={n_seeds}$, two-sided), "
        f"Holm-corrected jointly across all {pw_n_tests} tests (3 metrics $\\times$ {n_pairs} "
        f"pairs), with matched-pairs rank-biserial correlation $r_{{rb}}$.}}")
    add(r'\label{tab:wilcoxon}')
    add(r'\small')
    add(r'\setlength{\tabcolsep}{4pt}')
    add(r'\begin{tabular}{lccc}')
    add(r'\toprule')
    add(r'Comparison & HV $(p_{\mathrm{Holm}}, r_{rb})$ & IGD $(p_{\mathrm{Holm}}, r_{rb})$ & Spacing $(p_{\mathrm{Holm}}, r_{rb})$ \\')
    add(r'\midrule')
    pw = stats['pairwise_wilcoxon_holm_global']
    pair_labels = list(pw['hv'].keys())

    def fmt_p(p):
        # FIX (external review round 6-followup, confirmed real): a p-value
        # must never be printed as literally "0.000" (implies impossible
        # certainty) -- anything rounding to zero at 3 decimals is now
        # shown as "<0.001" instead.
        return '<0.001' if p < 0.0005 else f'{p:.3f}'

    for lbl in pair_labels:
        cells = []
        for m in ['hv', 'igd', 'spacing']:
            r = pw[m][lbl]
            sig = '^{**}' if r['p_holm_global'] < 0.01 else ('^{*}' if r['p_holm_global'] < 0.05 else '')
            cells.append(f"({fmt_p(r['p_holm_global'])}{sig}, {r['rank_biserial']:+.2f})")
        nice_lbl = (lbl.replace('nsga2', 'NSGA-II').replace('nsga3', 'NSGA-III')
                       .replace('moead', 'MOEA/D').replace('sms_emoa', 'SMS-EMOA')
                       .replace('rvea', 'RVEA'))
        add(f"{nice_lbl} & ${cells[0]}$ & ${cells[1]}$ & ${cells[2]}$ \\\\")
    add(r'\bottomrule')
    add(r'\end{tabular}')
    add(r'\end{table}')
    add()

    if 'independent_samples_sensitivity_check' in stats:
        ind = stats['independent_samples_sensitivity_check']
        add('%% === Independent-samples sensitivity check (inline text data,'
            ' external review round 6, item 9) ===')
        for m, key in [('HV', 'hv'), ('IGD', 'igd'), ('Spacing', 'spacing')]:
            kw = ind['kruskal'][key]
            add(f"%% Kruskal-Wallis {m}: H={kw['statistic']:.3f}, p={kw['p']:.5f}")
        n_agree = 0
        n_total = 0
        for m in ['hv', 'igd', 'spacing']:
            for lbl in pair_labels:
                n_total += 1
                paired_sig = pw[m][lbl]['p_holm_global'] < 0.05
                ind_sig = ind['pairwise_mannwhitney_holm_global'][m][lbl]['p_holm_global'] < 0.05
                if paired_sig == ind_sig:
                    n_agree += 1
        add(f"%% Paired (Wilcoxon) vs. independent-samples (Mann-Whitney) significance "
            f"pattern agrees on {n_agree}/{n_total} comparisons at alpha=0.05 "
            f"(both Holm-corrected jointly).")
        add()

abl_path = RESULTS / 'ablation_v3_results.json'
if abl_path.exists():
    with open(abl_path) as f:
        abl = json.load(f)
    variants = abl['config']['variants']
    abl_n_seeds = len(abl['config']['seeds'])
    add(f'%% === TABLE: ablation ({len(variants)} variants, {abl_n_seeds} seeds) ===')
    add(r'\begin{table}[tbp]')
    add(r'\centering')
    add(r'\caption{Decoder ablation: eight variants evaluated with NSGA-III '
        f'(${abl_n_seeds}$ seeds, mean$\\,\\pm\\,$std, canonical oracle).}}')
    add(r'\label{tab:ablation}')
    add(r'\small')
    add(r'\setlength{\tabcolsep}{3pt}')
    add(r'\resizebox{\columnwidth}{!}{%')
    add(r'\begin{tabular}{lcccc}')
    add(r'\toprule')
    add(r'Variant & HV $\uparrow$ & IGD $\downarrow$ & Spacing $\downarrow$ & $|F^\star|$ \\')
    add(r'\midrule')
    s = abl['summary']
    best_v = max(variants, key=lambda v: s[v]['hv']['mean'] or 0)
    names = {
        'A_full_decoder': 'A. Full decoder', 'B_fiber_only': 'B. Fibre only',
        'C_fiber_satellite': 'C. Fibre + satellite', 'D_fiber_fso': 'D. Fibre + FSO',
        'E_no_interregion': 'E. No inter-region bonus',
        'G_satellite_only': 'G. Satellite only', 'H_satellite_fso': 'H. Satellite + FSO',
        'I_fso_only': 'I. FSO only',
    }
    for v in variants:
        m = s[v]
        bold = v == best_v
        def fmt(key, dec=4):
            val, sd = m[key]['mean'], m[key]['std']
            txt = f'{val:.{dec}f} \\pm {sd:.{dec}f}'
            return f'\\mathbf{{{txt}}}' if bold else txt
        npar = f"{m['n_pareto']['mean']:.0f} \\pm {m['n_pareto']['std']:.0f}"
        add(f"{names.get(v,v):<28s} & ${fmt('hv')}$ & ${fmt('igd')}$ & ${fmt('spacing')}$ & ${npar}$ \\\\")
    add(r'\bottomrule')
    add(r'\end{tabular}%')
    add(r'}')
    add(r'\end{table}')
    add()

    if stats_path.exists() and 'pairwise_vs_full' in stats.get('ablation', {}):
        add('%% === Ablation pairwise-vs-full stats (inline text data) ===')
        add(f"%% Friedman omnibus (paired by seed): {stats['ablation'].get('friedman')}")
        for v, r in stats['ablation']['pairwise_vs_full'].items():
            sig = 'differs significantly' if r['p_holm'] < 0.05 else 'not significantly different'
            add(f"%% {v}: p_holm={r['p_holm']:.4f} delta={r['cliffs_delta']:+.3f} "
                f"r_rb={r['rank_biserial']:+.3f} -> {sig}")
        add()

cb_path = RESULTS / 'classical_baselines_v3_results.json'
if cb_path.exists():
    with open(cb_path) as f:
        cb = json.load(f)
    add('%% === Classical baselines (v3, raw data) ===')
    def row(name, d):
        if not d or 'fidelity' not in d:
            return None
        return (name, d['n_edges'], d['fidelity'], d['rate_log'], d['latency'], d['coverage'])
    all_rows = []
    if 'mst' in cb: all_rows.append(row('MST', cb['mst']))
    for k, v in cb.get('knn', {}).items(): all_rows.append(row(f'k-NN {k}', v))
    for k, v in cb.get('rgg', {}).items(): all_rows.append(row(f'RGG {k}', v))
    for k, v in cb.get('hub_and_spoke', {}).items(): all_rows.append(row(f'Hub {k}', v))
    if 'region_mst' in cb: all_rows.append(row('Region-MST', cb['region_mst']))
    for r in all_rows:
        if r:
            add(f"%% {r[0]:15s} |E|={r[1]:5d} F={r[2]:.3f} Rlog={r[3]:.3f} L={r[4]:.4f} C={r[5]:.3f}")
    add()

bb_path = RESULTS / 'binary_baseline_v3_results.json'
if bb_path.exists():
    with open(bb_path) as f:
        bb = json.load(f)
    add('%% === Binary/categorical baseline (v3) summary ===')
    if 'metrics_vs_shared_reference' in bb:
        add(f"%% metrics_vs_shared_reference={json.dumps(bb['metrics_vs_shared_reference'])}")
    for k in ['nsga2_categorical', 'nsga3_categorical']:
        if k in bb:
            v = bb[k]
            sizes = [r['n_pareto'] for r in v]
            times = [r['elapsed_min'] for r in v]
            add(f"%% {k}: |F*|={np.mean(sizes):.1f}+/-{np.std(sizes):.1f}  "
                f"time={np.mean(times):.2f}+/-{np.std(times):.2f} min")
    add()

gp_path = RESULTS / 'global_pareto_v3_summary.json'
if gp_path.exists():
    with open(gp_path) as f:
        gp = json.load(f)
    add('%% === Global Pareto front summary (v3) ===')
    add(f"%% n_points={gp['n_points']}  by_algo={gp['by_source_algorithm']}")
    add(f"%% ranges={gp['ranges']}")
    add(f"%% correlations={gp['correlations']}")
    add(f"%% resource_use={gp['resource_use']}")
    add()

dom_path = RESULTS / 'dominance_check_v3.json'
if dom_path.exists():
    with open(dom_path) as f:
        dom = json.load(f)
    add('%% === Dominance check (v3): classical baselines vs global front ===')
    add(f"%% n_dominated={dom['n_dominated']}/{dom['n_classical']}")
    for row in dom['rows']:
        add(f"%%   {row['name']:20s} dominated={row['dominated']}")
    add()

with open(OUT, 'w') as f:
    f.write('\n'.join(lines))
print(f'Saved: {OUT}')
print('\n'.join(lines[:40]))
