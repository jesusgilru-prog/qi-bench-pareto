#!/usr/bin/env python3
"""
generate_figures.py — Regenerate all 3 paper figures from v2 results
==========================================================================
Replaces the earlier draft's NSGA-III-only figures with versions built
from the corrected 11-D, heterogeneous-baseline, global-front data.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / 'results'
# FIX (Round 3 rerun): this script lives at .../QI-Bench-Pareto-v1/code/, so
# parent.parent is QI-Bench-Pareto-v1, not paper2_v4 -- the real manuscript
# directory (where main.tex's \includegraphics resolves filenames) is one
# level further up. Confirmed by the actual figure files' prior location.
PAPER = Path(__file__).resolve().parent.parent.parent / 'paper2'
# FIX (external review round 6-followup, confirmed real): this script only
# ever wrote into PAPER, never into QI-Bench-Pareto-v1/figures/ -- the
# artefact's own figures/ directory silently went stale every time this
# script was rerun (confirmed: figures/ still held pre-Round-6, 3-algorithm
# plots while paper2/ had the current 5-algorithm ones, shipped together in
# the same review package). Now writes to both locations from one save call
# each, so there is exactly one generation step and no way for them to
# drift apart again.
ARTIFACT_FIGURES = Path(__file__).resolve().parent.parent / 'figures'


def save_both(fig_or_plt, filename):
    plt.savefig(PAPER / filename)
    plt.savefig(ARTIFACT_FIGURES / filename)


# FIX (external review: figures "occupy a small fraction of the page and
# their labels and points are too small to read at normal size"). Every
# figure here is drawn several inches wider than the 6.3 in text block it
# is typeset into, so all text in the PDF is downscaled by roughly a half.
# Sizing the defaults up compensates; figure widths in main.tex are set to
# \textwidth so nothing is scaled down further.
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 14,
})

# ── Figure 1: algorithm comparison (HV, IGD bar/error plots) ────────────────
# FIX (external review: SPEA2 and AGE-MOEA-II were reported in a side
# paragraph rather than integrated). This figure is now driven by
# unified_seven_algorithm_comparison_v3.json, in which all seven
# algorithms are scored against ONE shared reference front and
# normalisation, so the panels are directly comparable across all seven
# and match Table 3 exactly.
with open(RESULTS / 'unified_seven_algorithm_comparison_v3.json') as f:
    base = json.load(f)
# Ordered by mean HV, as in Table 3.
algos = base['ranking_by_mean_hv']
label = {'nsga2': 'NSGA-II', 'nsga3': 'NSGA-III', 'moead': 'MOEA/D',
         'sms_emoa': 'SMS-EMOA', 'rvea': 'RVEA', 'spea2': 'SPEA2',
         'age_moea2': 'AGE-MOEA-II'}
palette = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2',
           '#937860', '#DA8BC3']
algo_colors = [palette[i % len(palette)] for i in range(len(algos))]
summ = base['summary']

# FIX (external review round 6, item 8): a bar+errorbar plot hides
# individual-run structure -- a single high-HV outlier run can visually
# look identical to a genuinely higher, consistent mean. Switched to a
# boxplot (median/IQR, robust to outliers) with raw per-seed points
# jittered on top, so any outlier is visible as a specific labelled point
# rather than folded invisibly into a mean+std bar.
rng = np.random.RandomState(0)

def raw_vals(metric):
    return [[v for v in summ[a][f'{metric}_raw'] if v is not None and np.isfinite(v)]
            for a in algos]

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
x = np.arange(len(algos)) + 1
for ax, metric, title, ylab in [
    (axes[0], 'hv', '(a) Hypervolume', 'Hypervolume (higher is better)'),
    (axes[1], 'igd', '(b) IGD', 'IGD (lower is better)'),
]:
    vals = raw_vals(metric)
    bp = ax.boxplot(vals, positions=x, widths=0.5, showfliers=False, patch_artist=True)
    for patch, c in zip(bp['boxes'], algo_colors):
        patch.set_facecolor(c); patch.set_alpha(0.35)
    for i, v in enumerate(vals):
        jitter = rng.uniform(-0.12, 0.12, size=len(v))
        ax.scatter(np.full(len(v), x[i]) + jitter, v, color=algo_colors[i],
                   s=14, alpha=0.8, zorder=3, edgecolors='black', linewidths=0.3)
    # FIX (round 14, ChatGPT Round 13 finding #15, confirmed real): at
    # rotation=20 with no horizontal anchor, adjacent long labels
    # (AGE-MOEA-II / NSGA-III) visually touch. A steeper rotation with
    # right-anchored labels gives each label more vertical clearance.
    ax.set_xticks(x)
    ax.set_xticklabels([label[a] for a in algos], rotation=35, ha='right',
                        rotation_mode='anchor')
    ax.set_ylabel(ylab)
    ax.set_title(title)
plt.tight_layout()
save_both(plt, 'algo_comparison_v3.pdf')
plt.close()
print('Saved algo_comparison_v3.pdf')

# ── Figure 2: ablation bar chart ────────────────────────────────────────────
abl_path = RESULTS / 'ablation_v3_results.json'
if abl_path.exists():
    with open(abl_path) as f:
        abl = json.load(f)
    variants = abl['config']['variants']
    abl_n_seeds = len(abl['config']['seeds'])
    # FIX (external review round 6-followup, confirmed real): this title was
    # a hardcoded "10 seeds" string, missed when the ablation was raised to
    # 30 seeds -- the compiled manuscript showed this figure's own title
    # contradicting its caption on the same page. Made dynamic, and (same
    # rationale as Figure 1's Round-6 fix) switched from a bar+errorbar
    # plot to a boxplot with jittered per-seed points now that n=30 makes
    # the per-seed distribution worth showing directly.
    per_run = abl['per_run_details']

    def abl_raw_vals(metric):
        return [[r[metric] for r in per_run[v] if r.get(metric) is not None] for v in variants]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    xv = np.arange(len(variants)) + 1
    vals = abl_raw_vals('hv')
    bp = ax.boxplot(vals, positions=xv, widths=0.5, showfliers=False, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('#4C72B0'); patch.set_alpha(0.35)
    for i, v in enumerate(vals):
        jitter = rng.uniform(-0.12, 0.12, size=len(v))
        ax.scatter(np.full(len(v), xv[i]) + jitter, v, color='#4C72B0',
                   s=12, alpha=0.8, zorder=3, edgecolors='black', linewidths=0.3)
    ax.set_xticks(xv)
    ax.set_xticklabels([v.split('_')[0] for v in variants])
    ax.set_ylabel('Hypervolume')
    ax.set_title(f'Decoder ablation (8 variants, {abl_n_seeds} seeds)')
    plt.tight_layout()
    save_both(plt, 'ablation_v3.pdf')
    plt.close()
    print('Saved ablation_v3.pdf')

# ── Figure 3: pairwise projections of the global Pareto front ──────────────
# FIX (Round 4 review, confirmed real): `source_algorithm` in the npz is
# only the FIRST-SEEN algorithm for each unique point (kept for
# backward-compat plotting, see build_global_pareto.py), so a point found
# by multiple algorithms was silently coloured as if exclusive to one --
# ~5% of points are actually shared. Now loads the real per-point
# provenance and adds an explicit fourth "Shared" category (distinct
# marker) instead of arbitrarily attributing shared points to whichever
# algorithm happened to run first in the pooling order.
gp = np.load(RESULTS / 'global_pareto_v3.npz', allow_pickle=True)
F = gp['F_pareto']  # -fid, -ratelog, +lat, -cov
src = gp['source_algorithm']
try:
    with open(RESULTS / 'global_pareto_v3_provenance.json') as f:
        _prov = json.load(f)
    _n_algos_per_point = np.array([
        len(set(s['algorithm'] for s in p['sources'])) for p in _prov['per_point_provenance']
    ])
    src_category = np.array([
        'shared' if n > 1 else s for n, s in zip(_n_algos_per_point, src)
    ])
except FileNotFoundError:
    src_category = src  # fall back to first-seen-only if provenance file is absent
fid, ratelog, lat, cov = -F[:, 0], -F[:, 1], F[:, 2], -F[:, 3]
data = {'Fidelity': fid, 'Log-rate': ratelog, 'Latency (s)': lat, 'Coverage': cov}
names = list(data.keys())
# FIX (round 14, ChatGPT Round 13 finding #1, confirmed real bug): the
# released pooled front now unions all SEVEN algorithms
# (build_global_pareto.py), but this line still hardcoded the original
# five, so points whose first-seen/exclusive source was SPEA2 or
# AGE-MOEA-II matched no entry in the `for algo in front_algos + ['shared']`
# scatter loop below and were silently dropped from the figure entirely
# (not merely missing from the legend). Derives its list from whichever
# algorithms actually appear in the data instead of a fixed five.
front_algos = [a for a in algos if (src_category == a).any()]
front_colors = [palette[i % len(palette)] for i in range(len(front_algos))]
colors = dict(zip(front_algos, front_colors))
colors['shared'] = '#000000'
marker_shapes = ['o', '^', 's', 'D', 'v', 'P', 'X']
markers = {a: marker_shapes[i % len(marker_shapes)] for i, a in enumerate(front_algos)}
markers['shared'] = '*'

# FIX (external review: "the 16-panel matrix is hard to read at normal
# size; labels and points are too small"). The figure is drawn at 12 in
# wide and typeset at ~6.3 in, so everything is downscaled by ~0.5 in the
# PDF. Axis labels, ticks and markers are sized here to survive that
# downscale rather than to look right in the standalone PDF.
LBL_FS, TICK_FS = 15, 12
fig, axes = plt.subplots(4, 4, figsize=(12, 12.8))
for i in range(4):
    for j in range(4):
        ax = axes[i, j]
        if i == j:
            ax.hist(data[names[i]], bins=20, color='gray')
        else:
            for algo in list(front_algos) + ['shared']:
                lbl = label.get(algo, 'Shared')
                mask = src_category == algo
                sz = 70 if algo == 'shared' else 18
                ax.scatter(data[names[j]][mask], data[names[i]][mask],
                           c=colors[algo], marker=markers[algo], s=sz, alpha=0.8,
                           label=lbl if (i, j) == (1, 0) else None)
        ax.tick_params(axis='both', labelsize=TICK_FS)
        if i == 3:
            ax.set_xlabel(names[j], fontsize=LBL_FS)
        if j == 0:
            ax.set_ylabel(names[i], fontsize=LBL_FS)
handles, labels = axes[1, 0].get_legend_handles_labels()
# FIX (external review round 6-followup, confirmed real): with 5 algorithms
# + 'Shared' = 6 legend entries, ncol=3 wraps to 2 rows, and the reserved
# top margin (rect=[...,0.97] plus this bbox_to_anchor) only had room for
# one row before the figure was saved without bbox_inches='tight' -- the
# rendered PDF's legend was visibly cut down to 3 of 6 categories (NSGA-III,
# SMS-EMOA, Shared only). Fixed by using one row (ncol = number of entries)
# and saving with bbox_inches='tight' so nothing outside the axes is ever
# silently clipped regardless of legend size.
fig.legend(handles, labels, loc='upper center', ncol=len(handles),
           bbox_to_anchor=(0.5, 0.995), markerscale=2, fontsize=15)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(PAPER / 'pareto_pairwise_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(ARTIFACT_FIGURES / 'pareto_pairwise_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.close()
print('Saved pareto_pairwise_v3.pdf')

# Print correlations for the caption
print('rate_vs_latency:', np.corrcoef(ratelog, lat)[0, 1])
print('rate_vs_coverage:', np.corrcoef(ratelog, cov)[0, 1])
print('fidelity_vs_rate:', np.corrcoef(fid, ratelog)[0, 1])
