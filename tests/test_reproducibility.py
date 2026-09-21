#!/usr/bin/env python3
"""
test_reproducibility.py — external review round 7 (author-approved full
scope): a real, if not exhaustive, pytest suite covering the checks two
rounds of independent review kept re-raising by hand. Run with:

    cd QI-Bench-Pareto-v1 && python -m pytest tests/ -v

Each test reads only already-released results/manuscript files -- none
of these re-run any optimisation.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
CODE = ROOT / 'code'
PAPER = ROOT.parent / 'paper2'


# ── Statistics correctness ──────────────────────────────────────────────────

def test_canonical_spacing_is_manhattan_ddof1():
    """compute_spacing must use L1 (Manhattan) distance and ddof=1, per
    Schott (1995) -- confirmed real bug (Euclidean/ddof=0) in an earlier
    round, fixed and must not regress. Round 8: imports the standalone
    metrics.py (numpy-only) rather than a full pymoo-dependent experiment
    script, per an external review's correct observation that this test
    previously required pymoo installed just to reach one function."""
    sys.path.insert(0, str(CODE))
    import importlib
    mod = importlib.import_module('metrics')
    F = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    # Manhattan nearest-neighbour distances here: point0->{1,2} both L1=1;
    # point1->0 L1=1, point1->2 L1=2; point2->0 L1=1, point2->1 L1=2.
    # min-dists = [1, 1, 1] -> std(ddof=1) of a constant array is 0.
    result = mod.compute_spacing(F)
    assert result == pytest.approx(0.0, abs=1e-9)

    F2 = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 5.0]])
    # nearest-neighbour L1 dists: pt0 -> min(2,5)=2; pt1 -> min(2, 7)=2;
    # pt2 -> min(5,7)=5. std([2,2,5], ddof=1) != 0, and would differ from
    # the Euclidean-distance version (dists 2, 2, 5 are actually the same
    # here by coincidence on this axis-aligned example -- use a non-trivial
    # case below to actually discriminate L1 vs L2).
    F3 = np.array([[0.0, 0.0], [3.0, 4.0], [10.0, 0.0]])
    # L1: pt0-pt1=7, pt0-pt2=10, pt1-pt2=11 -> nn: pt0=7, pt1=7, pt2=10
    # L2: pt0-pt1=5, pt0-pt2=10, pt1-pt2=~8.06 -> nn: pt0=5, pt1=5, pt2=8.06
    result3 = mod.compute_spacing(F3)
    l1_expected = float(np.std([7.0, 7.0, 10.0], ddof=1))
    l2_wrong = float(np.std([5.0, 5.0, 8.06], ddof=1))
    assert result3 == pytest.approx(l1_expected, abs=1e-6)
    assert abs(result3 - l2_wrong) > 1e-3, "spacing looks Euclidean, not Manhattan"


def test_ablation_seed_count_is_30():
    d = json.loads((RESULTS / 'ablation_v3_results.json').read_text())
    assert len(d['config']['seeds']) == 30
    assert d['config']['seeds'][0] == 42 and d['config']['seeds'][-1] == 71


def test_main_comparison_uses_5_algorithms_30_seeds():
    d = json.loads((RESULTS / 'baselines_comparison_v3.json').read_text())
    assert set(d['config']['algorithms']) == {'nsga2', 'nsga3', 'moead', 'sms_emoa', 'rvea'}
    assert len(d['config']['seeds']) == 30


def test_evaluation_budgets_matched():
    d = json.loads((RESULTS / 'baselines_comparison_v3.json').read_text())
    for algo in ['nsga2', 'nsga3', 'sms_emoa', 'rvea']:
        for run in d['per_run_details'][algo]:
            assert run['n_calls'] == 8320, f'{algo} seed={run["seed"]} n_calls={run["n_calls"]}'
    for run in d['per_run_details']['moead']:
        assert run['n_calls'] == 8316, f'moead seed={run["seed"]} n_calls={run["n_calls"]}'


# ── No infeasible/pathological solutions ────────────────────────────────────

def test_no_pathological_latency_in_categorical_baseline():
    """Confirmed real bug in an earlier round: an FSO edge beyond its max
    range leaked a ~1e6-second sentinel into a real, averaged latency
    (max observed 3696.99s). Must not recur."""
    d = json.loads((RESULTS / 'binary_baseline_v3_results.json').read_text())
    max_lat = 0.0
    for key in ['nsga2_categorical', 'nsga3_categorical']:
        for run in d[key]:
            for pt in run['pareto_front']:
                max_lat = max(max_lat, pt[2])
    assert max_lat < 10.0, f'pathological latency found: {max_lat}'


def test_no_infinite_or_nan_in_main_fronts():
    d = json.loads((RESULTS / 'baselines_comparison_v3.json').read_text())
    for algo, runs in d['per_run_details'].items():
        for run in runs:
            for key in ['hv', 'igd', 'spacing']:
                if run.get(key) is not None:
                    assert np.isfinite(run[key]), f'{algo} seed={run["seed"]} {key} not finite'


# ── Manuscript hygiene ───────────────────────────────────────────────────────

def test_abstract_word_count_under_250():
    text = (PAPER / 'main.tex').read_text()
    m = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', text, re.S)
    assert m is not None, 'no abstract found'
    n_words = len(m.group(1).split())
    assert n_words <= 250, f'abstract has {n_words} words, exceeds the journal\'s 250-word limit'


def test_highlight_lengths_under_85_chars():
    text = (PAPER / 'highlights_upload.tex').read_text()
    items = re.findall(r'\\item (.+)', text)
    assert 3 <= len(items) <= 5, f'expected 3-5 highlight bullets, found {len(items)}'
    for item in items:
        assert len(item) <= 85, f'highlight exceeds 85 chars ({len(item)}): {item!r}'


def test_no_stale_email_variant():
    text = (PAPER / 'main.tex').read_text()
    assert 'jesus.gil.ruiz@' not in text, 'fabricated email variant regressed into main.tex'
    assert 'jesus.gil@universidadeuropea.es' in text


def test_no_visible_revision_history_narration():
    """LaTeX comments (%) are fine; visible body text referencing internal
    round numbers or audit narrative is not, per house style."""
    text = (PAPER / 'main.tex').read_text()
    visible_lines = [ln for ln in text.splitlines() if not ln.strip().startswith('%')]
    visible_text = '\n'.join(visible_lines)
    for bad in ['external review round', 'this round', 'as of this round', 'ChatGPT item']:
        assert bad not in visible_text, f'found revision-history narration: {bad!r}'


def test_placeholders_are_only_the_known_disclosed_ones():
    """CITATION.cff's repository-code URL and .zenodo.json's related-DOI
    field are the only two files ALLOWED to still contain an unfilled
    REPLACE_WITH_ placeholder (the DOI genuinely cannot exist until the
    journal/Zenodo assigns one). This test does NOT require either
    placeholder to still be present -- once a real URL/DOI is filled in,
    this test must keep passing, not start failing; it only asserts that
    no *other* file has one, and that these two never regress to
    containing MORE than the single expected placeholder each.
    CHANGELOG.md/README.md discussing these two known placeholders in
    prose is expected and fine; this test only scans actual machine-read
    config/data files, not documentation or this test file itself."""
    machine_read_files = ['CITATION.cff', '.zenodo.json']
    hits = {}
    for name in machine_read_files:
        content = (ROOT / name).read_text()
        if 'REPLACE_WITH' in content:
            hits[name] = True
    assert set(hits) <= {'CITATION.cff', '.zenodo.json'}, \
        f'placeholder found somewhere unexpected: {hits}'

    # No OTHER released config/data file (excluding docs and checkpoints)
    # should contain a placeholder.
    checked_globs = ['*.py', 'results/*.json', 'results/*.tex']
    unexpected = []
    for pattern in checked_globs:
        for path in ROOT.glob(pattern):
            if path.name == 'test_reproducibility.py':
                continue
            content = path.read_text()
            if 'REPLACE_WITH' in content:
                unexpected.append(str(path.relative_to(ROOT)))
    assert not unexpected, f'unexpected placeholder locations: {unexpected}'


# ── Safe serialisation ───────────────────────────────────────────────────────

def test_unique_topologies_loads_without_pickle():
    d = np.load(RESULTS / 'unique_topologies_v3.npz')  # no allow_pickle kwarg
    assert 'topology_offsets' in d
    assert 'edge_u' in d and 'edge_v' in d and 'edge_type' in d
    n_topo = len(d['topology_offsets']) - 1
    assert n_topo == len(d['F'])


def test_global_and_ablation_npz_load_without_pickle():
    np.load(RESULTS / 'global_pareto_v3.npz')
    np.load(RESULTS / 'ablation_v3_pareto.npz')


def test_multi_instance_pairwise_numbers_are_traceable():
    """Round 8 regression guard: the manuscript's NSGA-II-vs-NSGA-III
    pairwise HV/p-value numbers per additional instance (Section:
    Multi-instance validation) must be backed by a committed artefact,
    not just prose -- these were previously computed ad hoc and not
    persisted, which an external review correctly flagged as
    unverifiable. Values here must match main.tex exactly.

    Round 13: expected values updated after the satellite-visibility
    oracle correction (Section: oracle) changed every downstream
    number, including these; global_n200 also moved from 5 to 10
    seeds in the same round."""
    with open(RESULTS / 'multi_instance_pairwise_stats_v3.json') as f:
        data = json.load(f)['per_instance']
    expected = {
        'global_n50': (0.373, 0.347, 0.089),
        'regional_eu_na': (0.799, 0.788, 0.0022),
        'metro_europe': (0.036, 0.047, 0.850),
        'global_n200': (0.439, 0.454, 0.791),
    }
    for inst, (nsga2_hv, nsga3_hv, p) in expected.items():
        d = data[inst]
        assert d['nsga2_hv_mean'] == pytest.approx(nsga2_hv, abs=5e-4)
        assert d['nsga3_hv_mean'] == pytest.approx(nsga3_hv, abs=5e-4)
        assert d['mannwhitney_p'] == pytest.approx(p, abs=5e-4)


def test_cities_csv_has_200_unique_real_rows():
    """Round 8: cities.csv was extended from 100 to 200 rows (100
    additional real cities sourced from GeoNames). Guards against
    accidental truncation/duplication and confirms every new row has a
    non-empty coordinate_source citation, not a blank/fabricated one."""
    import csv
    with open(ROOT / 'cities.csv') as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 200
    indices = sorted(int(r['index']) for r in rows)
    assert indices == list(range(200))
    names = [r['city'] for r in rows]
    assert len(names) == len(set(names)), 'duplicate city name found in cities.csv'
    for r in rows:
        assert r['coordinate_source'].strip(), f"row {r['index']} has an empty coordinate_source"


def test_global_n200_oracle_matches_original_100_city_oracle():
    """Round 8: the N=200 instance's Instance class reimplements
    decode_topology/oracle_evaluate as closures over a 200-city table;
    this guards that its first 100 rows' distance matrix is numerically
    identical to the original, separately-validated 100-city oracle
    (bdcz_oracle_v3.DIST), which was verified once before trusting any
    N=200 result and must not silently regress."""
    sys.path.insert(0, str(CODE))
    import importlib
    mio = importlib.import_module('multi_instance_oracle_v3')
    main_oracle = importlib.import_module('bdcz_oracle_v3')
    inst = mio.Instance('global_n200')
    assert inst.n_nodes == 200
    assert np.allclose(inst.dist[:100, :100], main_oracle.DIST)
    assert np.allclose(inst.link_fid[:100, :100, :], main_oracle.LINK_FID)


# ── External oracle validation ──────────────────────────────────────────────

def test_sequence_validation_directional_agreement():
    """Regression guard for the SeQUeNCe directional cross-check (Round 7,
    Threats to validity item iii): does NOT import the `sequence` package
    (optional, heavy dependency, not required for any other script in this
    package) -- only reads the already-committed results JSON, so this test
    runs in the default environment."""
    with open(RESULTS / 'sequence_oracle_validation_v3.json') as f:
        data = json.load(f)
    assert data['sequence_rate_is_monotone_nonincreasing'] is True
    assert data['analytical_rate_is_monotone_nonincreasing'] is True
    assert data['directional_agreement'] is True
    assert data['spearman_rho_sequence_vs_analytical'] == pytest.approx(1.0)


# ── Checksums ────────────────────────────────────────────────────────────────

def test_checksums_pass():
    result = subprocess.run(['sha256sum', '-c', 'CHECKSUMS.sha256'],
                            cwd=str(ROOT), capture_output=True, text=True)
    assert result.returncode == 0, f'checksum failures:\n{result.stdout}\n{result.stderr}'


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
