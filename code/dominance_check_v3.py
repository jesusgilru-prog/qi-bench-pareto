#!/usr/bin/env python3
"""
dominance_check.py — Fair classical/binary-baseline dominance check
========================================================================
Checks whether each classical-constructor point (now using heterogeneous
BEST_MEDIUM, not hardcoded fibre) is Pareto-dominated by the corrected
global non-dominated front (union of NSGA-II/NSGA-III/MOEA-D, 11-D
oracle). Reports the count/fraction dominated, honestly, rather than
asserting "every classical baseline is dominated" without rechecking
after the fairness fix — the heterogeneous-medium classical baselines
are substantially stronger than the old fibre-only ones, so this needs
to be re-verified rather than assumed.
"""
import json
import numpy as np
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / 'results'

global_pareto = np.load(RESULTS / 'global_pareto_v3.npz')
global_F = global_pareto['F_pareto']  # minimisation space: -fid,-ratelog,+lat,-cov

with open(RESULTS / 'classical_baselines_v3_results.json') as f:
    cb = json.load(f)

rows = []
if 'mst' in cb and 'F_minimization' in cb['mst']:
    rows.append(('MST', cb['mst']['n_edges'], cb['mst']['F_minimization']))
for k, v in cb.get('knn', {}).items():
    if 'F_minimization' in v:
        rows.append((f'k-NN {k}', v['n_edges'], v['F_minimization']))
for k, v in cb.get('rgg', {}).items():
    if 'F_minimization' in v:
        rows.append((f'RGG {k}', v['n_edges'], v['F_minimization']))
for k, v in cb.get('hub_and_spoke', {}).items():
    if 'F_minimization' in v:
        rows.append((f'Hub-and-spoke {k}', v['n_edges'], v['F_minimization']))
if 'region_mst' in cb and 'F_minimization' in cb['region_mst']:
    rows.append(('Region-MST', cb['region_mst']['n_edges'], cb['region_mst']['F_minimization']))


def is_dominated(point, front):
    """True if `front` contains a point that dominates `point` (minimisation)."""
    p = np.array(point)
    return bool(np.any(np.all(front <= p, axis=1) & np.any(front < p, axis=1)))


print(f"{'Method':<22s} {'|E|':>6s} {'Fid':>7s} {'Rlog':>7s} {'Lat':>8s} {'Cov':>7s}  Dominated?")
n_dominated = 0
for name, n_edges, f_min in rows:
    dom = is_dominated(f_min, global_F)
    n_dominated += dom
    fid, rlog, lat, cov = -f_min[0], -f_min[1], f_min[2], -f_min[3]
    print(f"{name:<22s} {n_edges:>6d} {fid:>7.3f} {rlog:>7.3f} {lat:>8.4f} {cov:>7.3f}  "
          f"{'YES' if dom else 'NO (non-dominated!)'}")

print(f"\n{n_dominated}/{len(rows)} classical constructors dominated by the "
      f"global {len(global_F)}-point front.")

with open(RESULTS / 'dominance_check_v3.json', 'w') as f:
    json.dump({
        'n_classical': len(rows),
        'n_dominated': int(n_dominated),
        'rows': [{'name': r[0], 'n_edges': r[1], 'dominated': is_dominated(r[2], global_F)} for r in rows],
    }, f, indent=2)
print(f"Saved: {RESULTS / 'dominance_check_v3.json'}")
