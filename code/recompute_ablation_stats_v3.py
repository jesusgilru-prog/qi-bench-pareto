#!/usr/bin/env python3
"""
SUPERSEDED as of the round 6-followup 30-seed ablation rerun: run_ablation_v3.py's
own aggregate() now computes ddof=1 stats AND a genuine per-seed igd_union
directly (the raw per-seed F arrays this script's workaround needed are no
longer missing, since the ablation was actually rerun rather than patched
post-hoc). Do NOT run this script against the current, 30-seed
ablation_v3_results.json -- it expects the old 10-seed shape and will
overwrite the real per-seed igd_union with a coarser pooled approximation.
Kept only for historical reference; not part of the reproduction steps in
README.md any more.

recompute_ablation_stats_v3.py — external review round 6, items 10 and 15,
applied to the EXISTING ablation_v3_results.json / ablation_v3_pareto.npz
without re-running the 80 (8 variant x 10 seed) optimisations.

Item 10 (sample std): `per_run_details` already stores per-seed hv/igd/
spacing, so mean/std can be recomputed with ddof=1 (sample std) directly
from that -- no rerun needed.

Item 15 (IGD_union): a genuinely symmetric, non-A-biased reference front,
built from the union of all eight variants' ALREADY-POOLED, deduplicated
combined fronts in ablation_v3_pareto.npz. This gives one IGD_union value
per variant (computed against the variant's own combined/deduplicated
front, not per-seed mean+/-std -- the raw per-seed F arrays needed for a
per-seed IGD_union were not retained after the original run and are not
worth a fresh 80-run optimisation campaign just for this secondary
statistic this round). Disclosed as a scope decision, not hidden.
"""
import json
import numpy as np
from pathlib import Path
from pymoo.util.nds.non_dominated_sorting import find_non_dominated
from pymoo.indicators.igd import IGD

from bdcz_oracle_v3 import dedup_front

RESULTS = Path(__file__).resolve().parent.parent / 'results'

with open(RESULTS / 'ablation_v3_results.json') as f:
    d = json.load(f)
variants = d['config']['variants']

# ── Item 10: ddof=1 recompute from stored per-seed values ──────────────────
for v in variants:
    for metric in ['hv', 'igd', 'spacing']:
        vals = [r[metric] for r in d['per_run_details'][v] if r.get(metric) is not None]
        if not vals:
            continue
        a = np.array(vals)
        std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
        old_std = d['summary'][v][metric]['std']
        d['summary'][v][metric]['std'] = std
        d['summary'][v][metric]['mean'] = float(np.mean(a))  # unchanged, recomputed for consistency
        print(f'{v:22s} {metric:8s}: ddof=0 std was {old_std:.5f} -> ddof=1 std {std:.5f}')

# ── Item 15: pooled IGD_union from the combined fronts in the npz ──────────
npz = np.load(RESULTS / 'ablation_v3_pareto.npz', allow_pickle=True)
g_min, g_max = npz['g_min'], npz['g_max']
g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

def normalise(F):
    return (F - g_min) / g_range

all_F = np.vstack([npz[f'F_{v}'] for v in variants if len(npz[f'F_{v}']) > 0])
all_X = np.vstack([npz[f'X_{v}'] for v in variants if len(npz[f'X_{v}']) > 0])
nd = find_non_dominated(all_F)
ref_raw, _ = dedup_front(all_F[nd], all_X[nd])
ref_front_union = normalise(ref_raw)
igd_union_calc = IGD(ref_front_union)

print(f'\nShared union reference front: {len(ref_front_union)} deduplicated points '
      f'(pooled across all {len(variants)} variants).')
print('Pooled IGD_union per variant (against the shared front, one value per variant, '
      'not mean+/-std over seeds -- see module docstring):')
igd_union_out = {}
for v in variants:
    Fv = normalise(npz[f'F_{v}'])
    if len(Fv) == 0:
        continue
    val = float(igd_union_calc(Fv))
    igd_union_out[v] = val
    print(f'  {v:22s}: IGD_union = {val:.4f}  (vs.\\ existing A-biased IGD = '
          f'{d["summary"][v]["igd"]["mean"]:.4f})')
    d['summary'][v]['igd_union_pooled'] = val

with open(RESULTS / 'ablation_v3_results.json', 'w') as f:
    json.dump(d, f, indent=2)
print(f'\nUpdated {RESULTS / "ablation_v3_results.json"} in place '
      f'(ddof=1 stds + igd_union_pooled field added).')
