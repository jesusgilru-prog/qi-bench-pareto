#!/usr/bin/env python3
"""
generate_topology_maps_v3.py — external review round 7 (author-approved
full scope): geographic maps of 4 representative topologies, since the
paper otherwise never shows an actual designed network.

1. Max-fidelity solution (from the pooled evolutionary front).
2. Max-coverage solution (from the pooled evolutionary front).
3. Knee-point solution (closest to the ideal corner in normalised,
   minimisation-convention objective space, L2 distance).
4. The non-dominated classical constructor, rebuilt deterministically
   via the exact same construction as classical_baselines_v3.py.
   FIX (post satellite-visibility fix, Round 13): under the corrected
   oracle every hub-and-spoke variant is now dominated by the pooled
   evolutionary front (their long hub-to-spoke satellite links lose
   most of their reach once a single satellite can no longer bridge
   arbitrary distances); the Random Geometric Graph at a 5,000 km
   threshold is the classical constructor that is now non-dominated
   (confirmed in Table~tab:classical-baselines), so panel (d) plots
   that instead of the 10-hub hub-and-spoke network.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

from bdcz_oracle_v3 import DIST, N_NODES, BEST_MEDIUM, _build_best_medium_adjacency


def n_effective_edges(eu, ev, et):
    """Unique-pair / effective edge count (external review, ChatGPT Round 13
    finding #14, confirmed real): |E| as previously reported for panels
    (a)-(c) was the RAW medium-specific count (a pair with both fibre AND
    satellite counts twice), while panel (d)'s RGG has exactly one medium
    per pair by construction, so its raw count already equals its
    effective count. Reporting only the raw count for (a)-(c) therefore
    made all four panels look directly comparable on network density when
    they were not. Mirrors build_global_pareto.py's n_effective_edges."""
    if len(eu) == 0:
        return 0
    adj = _build_best_medium_adjacency(eu, ev, et)
    return sum(len(nbrs) for nbrs in adj) // 2

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

MEDIUM_COLOR = {0: '#4C72B0', 1: '#DD8452', 2: '#55A868'}  # fibre, satellite, FSO
MEDIUM_LABEL = {0: 'Fibre', 1: 'Satellite', 2: 'FSO'}

# ── Load city coordinates ───────────────────────────────────────────────────
# cities.csv now carries all 200 rows (the original 100 plus the 100
# GeoNames-sourced cities appended for the global_n200 instance,
# Section: Multi-instance validation), while these maps depict the main
# 100-city benchmark. Rows with index >= N_NODES belong to the scale-up
# instance and are skipped here rather than overflowing the arrays.
lats, lons = np.zeros(N_NODES), np.zeros(N_NODES)
with open(Path(__file__).resolve().parent.parent / 'cities.csv') as f:
    for row in csv.DictReader(f):
        i = int(row['index'])
        if i >= N_NODES:
            continue
        lats[i] = float(row['latitude'])
        lons[i] = float(row['longitude'])


def plot_topology(ax, eu, ev, et, title):
    ax.set_global()
    ax.coastlines(linewidth=0.5, color='gray')
    for medium in [0, 1, 2]:
        mask = et == medium
        for u, v in zip(eu[mask], ev[mask]):
            ax.plot([lons[u], lons[v]], [lats[u], lats[v]],
                    color=MEDIUM_COLOR[medium], linewidth=0.9, alpha=0.75,
                    transform=ccrs.Geodetic())
    active_nodes = np.unique(np.concatenate([eu, ev])) if len(eu) > 0 else np.array([], dtype=int)
    ax.scatter(lons[active_nodes], lats[active_nodes], s=6, color='black',
               zorder=5, transform=ccrs.PlateCarree())
    ax.set_title(title, fontsize=9)


# ── 1-3: evolutionary solutions from the pooled front ──────────────────────
# FIX (safe serialisation, external review round 7): unique_topologies_v3.npz
# now stores edges as a flat (edge_u, edge_v, edge_type) triple plus a
# topology_offsets CSR-style index, not per-topology object arrays --
# loads with allow_pickle=False. See build_global_pareto.py.
gp = np.load(RESULTS / 'unique_topologies_v3.npz')
F = gp['F']
offsets = gp['topology_offsets']


def topo_edges(idx):
    lo, hi = offsets[idx], offsets[idx + 1]
    return gp['edge_u'][lo:hi], gp['edge_v'][lo:hi], gp['edge_type'][lo:hi]


fid, ratelog, lat, cov = -F[:, 0], -F[:, 1], F[:, 2], -F[:, 3]

idx_maxfid = int(np.argmax(fid))
idx_maxcov = int(np.argmax(cov))

g_min, g_max = F.min(axis=0), F.max(axis=0)
g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)
Fn = (F - g_min) / g_range  # minimisation convention, normalised to [0,1]
dist_to_ideal = np.linalg.norm(Fn, axis=1)
idx_knee = int(np.argmin(dist_to_ideal))

eu_mf, ev_mf, et_mf = topo_edges(idx_maxfid)
eu_mc, ev_mc, et_mc = topo_edges(idx_maxcov)
eu_kn, ev_kn, et_kn = topo_edges(idx_knee)
eff_mf = n_effective_edges(eu_mf, ev_mf, et_mf)
eff_mc = n_effective_edges(eu_mc, ev_mc, et_mc)
eff_kn = n_effective_edges(eu_kn, ev_kn, et_kn)
print(f'max-fidelity idx={idx_maxfid}  F=fid{fid[idx_maxfid]:.3f} rate{ratelog[idx_maxfid]:.2f} '
      f'lat{lat[idx_maxfid]:.4f} cov{cov[idx_maxfid]:.3f}  |E|raw={len(eu_mf)} |E|eff={eff_mf}')
print(f'max-coverage idx={idx_maxcov}  F=fid{fid[idx_maxcov]:.3f} rate{ratelog[idx_maxcov]:.2f} '
      f'lat{lat[idx_maxcov]:.4f} cov{cov[idx_maxcov]:.3f}  |E|raw={len(eu_mc)} |E|eff={eff_mc}')
print(f'knee idx={idx_knee}  F=fid{fid[idx_knee]:.3f} rate{ratelog[idx_knee]:.2f} '
      f'lat{lat[idx_knee]:.4f} cov{cov[idx_knee]:.3f}  |E|raw={len(eu_kn)} |E|eff={eff_kn}')

# ── 4: non-dominated classical constructor -- RGG at 5,000 km, rebuilt
# deterministically exactly as classical_baselines_v3.py builds it ────────
rgg_edges = sorted((i, j) for i in range(N_NODES) for j in range(i + 1, N_NODES)
                   if DIST[i, j] < 5000)
rgg_eu = np.array([e[0] for e in rgg_edges])
rgg_ev = np.array([e[1] for e in rgg_edges])
rgg_et = BEST_MEDIUM[rgg_eu, rgg_ev]
print(f'RGG (5,000 km): |E|={len(rgg_edges)}')

# ── Plot all four ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 7), subplot_kw={'projection': ccrs.Robinson()})

plot_topology(axes[0, 0], eu_mf, ev_mf, et_mf,
              f'(a) Max-fidelity ($F={fid[idx_maxfid]:.3f}$, $|E|_{{\\mathrm{{eff}}}}={eff_mf}$)')
plot_topology(axes[0, 1], eu_mc, ev_mc, et_mc,
              f'(b) Max-coverage ($C_\\tau={cov[idx_maxcov]:.3f}$, $|E|_{{\\mathrm{{eff}}}}={eff_mc}$)')
plot_topology(axes[1, 0], eu_kn, ev_kn, et_kn,
              f'(c) Knee point ($|E|_{{\\mathrm{{eff}}}}={eff_kn}$)')
plot_topology(axes[1, 1], rgg_eu, rgg_ev, rgg_et,
              f'(d) RGG (5,000 km), non-dominated ($|E|_{{\\mathrm{{eff}}}}={len(rgg_edges)}$)')

handles = [plt.Line2D([0], [0], color=MEDIUM_COLOR[m], lw=2, label=MEDIUM_LABEL[m]) for m in [0, 1, 2]]
fig.legend(handles=handles, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.0), fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(PAPER / 'topology_maps_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(ARTIFACT_FIGURES / 'topology_maps_v3.pdf', bbox_inches='tight', pad_inches=0.05)
plt.close()
print('Saved topology_maps_v3.pdf')

# Save selection metadata for the manuscript/reproducibility
out = {
    'max_fidelity': {'idx': idx_maxfid, 'fidelity': float(fid[idx_maxfid]), 'rate_log': float(ratelog[idx_maxfid]),
                     'latency': float(lat[idx_maxfid]), 'coverage': float(cov[idx_maxfid]),
                     'n_edges_raw': int(len(eu_mf)), 'n_edges_effective': int(eff_mf)},
    'max_coverage': {'idx': idx_maxcov, 'fidelity': float(fid[idx_maxcov]), 'rate_log': float(ratelog[idx_maxcov]),
                     'latency': float(lat[idx_maxcov]), 'coverage': float(cov[idx_maxcov]),
                     'n_edges_raw': int(len(eu_mc)), 'n_edges_effective': int(eff_mc)},
    'knee_point': {'idx': idx_knee, 'fidelity': float(fid[idx_knee]), 'rate_log': float(ratelog[idx_knee]),
                  'latency': float(lat[idx_knee]), 'coverage': float(cov[idx_knee]),
                  'n_edges_raw': int(len(eu_kn)), 'n_edges_effective': int(eff_kn)},
    'rgg_5000km': {'n_edges_raw': len(rgg_edges), 'n_edges_effective': len(rgg_edges)},
    'note': ('n_edges_raw counts each (city-pair, medium) separately -- a pair with both '
             'fibre and satellite counts twice. n_edges_effective collapses parallel media '
             'per pair to the single highest-fidelity one (unique city pairs actually '
             'routed over); this is the figure captions\' |E| and the only one directly '
             'comparable to classical baselines, which have one medium per pair by '
             'construction. External review (ChatGPT Round 13, finding #14).'),
}
with open(RESULTS / 'topology_maps_v3_selection.json', 'w') as f:
    json.dump(out, f, indent=2)
print('Saved topology_maps_v3_selection.json')
