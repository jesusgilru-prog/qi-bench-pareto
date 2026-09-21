#!/usr/bin/env python3
"""
random_search_vs_decoder_shared_ref_v3.py — NSGA-II sobre el decodificador
frente a búsqueda aleatoria i.i.d., ambos puntuados contra el MISMO frente de
referencia y la misma normalización.

Este fichero existe porque `random_search_vs_decoder_shared_ref_v3.json` era
el único resultado citado por el manuscrito que ningún script regeneraba: se
calculó a mano en su día y sobrevivió a la recomputación posterior al fix
satelital con la fecha y los números antiguos, sin que nada avisara. Es
exactamente el mismo defecto de trazabilidad que ya se corrigió para las
comparaciones multi-instancia.

El frente de referencia se construye con la unión de los frentes crudos por
semilla de AMBOS competidores, de modo que los dos HV viven en la misma
escala. Esa escala NO es la de la tabla de comparación de algoritmos (que
normaliza contra los siete algoritmos), así que los valores de aquí solo son
comparables entre sí.
"""
import json
import pickle
from pathlib import Path

import numpy as np
from pymoo.indicators.hv import Hypervolume
from scipy.stats import mannwhitneyu

RESULTS = Path(__file__).resolve().parent.parent / 'results'
DECODER_CKPT = RESULTS / 'step1_baselines_v3'
RANDOM_CKPT = RESULTS / 'step8_random_search_v3'
OUT_JSON = RESULTS / 'random_search_vs_decoder_shared_ref_v3.json'

SEEDS = list(range(42, 72))
HV_REF = np.array([1.1] * 4)


def load_decoder(seed):
    with open(DECODER_CKPT / f'nsga2_seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)['F_pareto']


def load_random(seed):
    with open(RANDOM_CKPT / f'seed{seed}.pkl', 'rb') as f:
        return pickle.load(f)['F_pareto']


def cliffs_delta(a, b):
    a = np.asarray(a, float)[:, None]
    b = np.asarray(b, float)[None, :]
    return float(((a > b).sum() - (a < b).sum()) / (a.size * b.size))


def main():
    dec = [load_decoder(s) for s in SEEDS]
    rnd = [load_random(s) for s in SEEDS]

    all_F = np.vstack([F for F in dec + rnd if len(F) > 0])
    g_min, g_max = all_F.min(axis=0), all_F.max(axis=0)
    g_range = np.where((g_max - g_min) < 1e-12, 1e-9, g_max - g_min)

    def norm(F):
        return (F - g_min) / g_range

    hv = Hypervolume(ref_point=HV_REF)
    hv_dec = [float(hv(norm(F))) if len(F) else 0.0 for F in dec]
    hv_rnd = [float(hv(norm(F))) if len(F) else 0.0 for F in rnd]

    u, p = mannwhitneyu(hv_dec, hv_rnd, alternative='two-sided')
    md, mr = float(np.mean(hv_dec)), float(np.mean(hv_rnd))

    out = {
        'description': ('NSGA-II (decoder) vs. i.i.d. random search, SAME shared '
                        'reference front and normalisation built from the union of '
                        'both competitors\' raw per-seed fronts, 30 seeds each. These '
                        'HV values are NOT on the scale of the seven-algorithm '
                        'comparison table; only the within-file comparison is '
                        'meaningful.'),
        'n_seeds': len(SEEDS),
        'normalisation': {'g_min': g_min.tolist(), 'g_max': g_max.tolist()},
        'decoder_hv_mean': md,
        'decoder_hv_std': float(np.std(hv_dec, ddof=1)),
        'random_search_hv_mean': mr,
        'random_search_hv_std': float(np.std(hv_rnd, ddof=1)),
        'relative_advantage': (md - mr) / mr if mr else None,
        'cliffs_delta_decoder_vs_random': cliffs_delta(hv_dec, hv_rnd),
        'mannwhitney_u': float(u),
        'mannwhitney_p': float(p),
        'decoder_hv_raw': hv_dec,
        'random_search_hv_raw': hv_rnd,
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'Saved: {OUT_JSON}')
    print(f"  decoder       HV = {md:.4f} +- {out['decoder_hv_std']:.4f}")
    print(f"  random search HV = {mr:.4f} +- {out['random_search_hv_std']:.4f}")
    print(f"  ventaja relativa = {100 * out['relative_advantage']:.1f}%  "
          f"delta = {out['cliffs_delta_decoder_vs_random']:+.2f}  p = {p:.3g}")


if __name__ == '__main__':
    main()
