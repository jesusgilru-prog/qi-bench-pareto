#!/usr/bin/env python3
"""
analyze_effective_dimensionality.py — PCA on the global Pareto front's
objective vectors (external review round 6, ChatGPT item 11 / ESWA-style
"effective objective dimensionality" request).

The four objectives (fidelity, log-rate, latency, coverage) are not four
independent conflicts: this project's own review history already found
strong pairwise correlations among fidelity/rate/latency (all mutually
"aligned" once minimize/maximize directions are accounted for), with
coverage genuinely in tension with the other three. This script quantifies
that directly via a standardised PCA on the deduplicated global front,
rather than asserting a specific variance-explained number without having
computed it (a specific numeric claim from an external review was
deliberately NOT asserted in an earlier round for exactly this reason).
"""
import json
import numpy as np
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / 'results'

gp = np.load(RESULTS / 'global_pareto_v3.npz', allow_pickle=True)
F = gp['F_pareto']  # columns: -fidelity, -log_rate, +latency, -coverage (minimised)
# Convert to a "higher is better except latency" natural sign for interpretability,
# then standardise (zero mean, unit variance) before PCA -- PCA on raw objectives
# would be dominated by whichever objective happens to have the largest numeric
# range, not by which is most informative.
fid, ratelog, lat, cov = -F[:, 0], -F[:, 1], F[:, 2], -F[:, 3]
X = np.column_stack([fid, ratelog, lat, cov])
names = ['fidelity', 'log_rate', 'latency', 'coverage']

Xc = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
corr = np.corrcoef(Xc, rowvar=False)

# PCA via SVD on the standardised, centred data (equivalent to eigendecomposition
# of the correlation matrix for standardised inputs).
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
explained_var = (S ** 2) / (len(Xc) - 1)
explained_ratio = explained_var / explained_var.sum()
cum_ratio = np.cumsum(explained_ratio)

print('Correlation matrix (fidelity, log_rate, latency, coverage):')
print(np.array2string(corr, precision=3, suppress_small=True))
print()
print('PCA explained variance ratio (standardised inputs):')
for i, (r, c) in enumerate(zip(explained_ratio, cum_ratio), 1):
    print(f'  PC{i}: {r*100:5.2f}%  (cumulative {c*100:5.2f}%)')
print()
print('Loadings (rows=PC, cols=objective):')
print('       ' + '  '.join(f'{n:>10s}' for n in names))
for i, row in enumerate(Vt, 1):
    print(f'  PC{i}: ' + '  '.join(f'{v:+10.3f}' for v in row))

out = {
    'n_points': int(len(X)),
    'objective_names': names,
    'correlation_matrix': corr.tolist(),
    'pca_explained_variance_ratio': explained_ratio.tolist(),
    'pca_cumulative_variance_ratio': cum_ratio.tolist(),
    'pca_loadings': Vt.tolist(),
}
out_path = RESULTS / 'effective_dimensionality_v3.json'
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2)
print(f'\nSaved: {out_path}')
