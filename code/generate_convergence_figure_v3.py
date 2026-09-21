#!/usr/bin/env python3
"""
generate_convergence_figure_v3.py — external review round 7 (author-
approved full scope): HV/IGD vs. evaluations for all 5 algorithms,
median + bootstrap 95% CI (10 seeds, see run_convergence_v3.py).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# FIX (external review: figure labels too small once the figure is
# downscaled into the 6.3 in text block). Same rcParams bump as
# generate_figures_v3.py.
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 14,
})


RESULTS = Path(__file__).resolve().parent.parent / 'results'
PAPER = Path(__file__).resolve().parent.parent.parent / 'paper2'
ARTIFACT_FIGURES = Path(__file__).resolve().parent.parent / 'figures'

LABEL = {'nsga2': 'NSGA-II', 'nsga3': 'NSGA-III', 'moead': 'MOEA/D',
         'sms_emoa': 'SMS-EMOA', 'rvea': 'RVEA'}
COLOR = {'nsga2': '#4C72B0', 'nsga3': '#DD8452', 'moead': '#55A868',
         'sms_emoa': '#C44E52', 'rvea': '#8172B2'}

d = json.load(open(RESULTS / 'convergence_v3_results.json'))
algos = d['config']['algorithms']

# FIX (MDPI proofreading, 2026-09-21): at figsize height 4.5in the rotated
# y-axis label of panel (a) was clipped by the tight bounding box
# ("...95% bootstrap C"). Taller figure so both labels fit in full.
fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
for algo in algos:
    curve = d['per_algorithm'][algo]
    n_eval = [c['n_eval'] for c in curve]
    hv = [c['hv_median'] for c in curve]
    hv_lo = [c['hv_ci_lo'] for c in curve]
    hv_hi = [c['hv_ci_hi'] for c in curve]
    axes[0].plot(n_eval, hv, color=COLOR[algo], label=LABEL[algo], linewidth=1.3)
    axes[0].fill_between(n_eval, hv_lo, hv_hi, color=COLOR[algo], alpha=0.15)

    igd = [c.get('igd_median') for c in curve if c.get('igd_median') is not None]
    n_eval_igd = [c['n_eval'] for c in curve if c.get('igd_median') is not None]
    igd_lo = [c['igd_ci_lo'] for c in curve if c.get('igd_ci_lo') is not None]
    igd_hi = [c['igd_ci_hi'] for c in curve if c.get('igd_ci_hi') is not None]
    axes[1].plot(n_eval_igd, igd, color=COLOR[algo], label=LABEL[algo], linewidth=1.3)
    axes[1].fill_between(n_eval_igd, igd_lo, igd_hi, color=COLOR[algo], alpha=0.15)

axes[0].set_xlabel('Objective evaluations')
axes[0].set_ylabel('Hypervolume (median, 95% bootstrap CI)')
axes[0].set_title('(a) HV vs. evaluations')
axes[1].set_xlabel('Objective evaluations')
axes[1].set_ylabel('IGD (median, 95% bootstrap CI)')
axes[1].set_title('(b) IGD vs. evaluations')
axes[1].legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig(PAPER / 'convergence_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(ARTIFACT_FIGURES / 'convergence_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.close()
print('Saved convergence_v3.pdf')

# Print key summary numbers for the manuscript text
for algo in algos:
    curve = d['per_algorithm'][algo]
    half_evals = curve[len(curve) // 2]['n_eval']
    hv_half = curve[len(curve) // 2]['hv_median']
    hv_final = curve[-1]['hv_median']
    pct_of_final_at_half = 100 * hv_half / hv_final if hv_final else float('nan')
    print(f'{algo:10s}: final n_eval={curve[-1]["n_eval"]}, HV_final={hv_final:.4f}, '
          f'HV_at_half_budget(n_eval={half_evals})={hv_half:.4f} '
          f'({pct_of_final_at_half:.1f}% of final)')
