# QI-Bench-Pareto-v1

Reproducibility package for *"Many-Objective Optimisation of Global Quantum
Internet Topologies: A Compact Parametric Representation with an Analytical
BDCZ Oracle"* (Jesús Gil Ruiz, Diego Rubén Rodríguez Regadera and Rafael Muñoz Gil;
accepted for publication in *Entropy*, MDPI, 2026).

## Contents

```
QI-Bench-Pareto-v1/
├── README.md                            (this file)
├── CHANGELOG.md                         itemised round-by-round audit trail
├── Makefile                             test/figures/tables/paper/verify/package targets
├── LICENSE-CODE                         MIT (code/)
├── LICENSE-DATA                         CC-BY-4.0 (cities.csv, results/)
├── CITATION.cff                         citation metadata
├── .zenodo.json                         Zenodo deposit metadata
├── cities.csv                           200 cities (100 original + 100 added Round 8, GeoNames-sourced):
│                                          index, city, country, lat, lon, region, coordinate_source
├── requirements.txt                     exact package versions used
├── tests/
│   └── test_reproducibility.py          Round 8: 18 pytest checks (spacing formula, budgets,
│                                          no-pathological-latency, manuscript hygiene, safe
│                                          serialisation, checksums) -- run via `make test`
├── code/
│   ├── bdcz_oracle_v3.py                       canonical BDCZ oracle + 11-var decoder
│   ├── run_baselines_comparison_v3.py          NSGA-II/NSGA-III/MOEA-D/SMS-EMOA/RVEA x 30 seeds
│   ├── run_ablation_v3.py                      8-variant decoder ablation x 30 seeds
│   ├── binary_baseline_v3.py                   fully categorical 4,950-var baseline x 30 seeds
│   ├── classical_baselines_v3.py               MST/kNN/RGG/hub-and-spoke baselines
│   ├── build_global_pareto.py                  global non-dominated front + resource stats
│   ├── build_binary_vs_decoder_shared_ref_v3.py  decoder vs. categorical, shared ref (Round 4)
│   ├── stats_analysis_v3.py                    Kruskal-Wallis + globally-Holm-corrected
│   │                                             Mann-Whitney (primary), plus a common-seed
│   │                                             paired (Friedman/Wilcoxon) sensitivity check
│   ├── dominance_check_v3.py                   classical-baseline vs global-front dominance
│   ├── generate_tables_v3.py                   LaTeX table generator from results JSON
│   ├── generate_figures_v3.py                  regenerates algo_comparison/ablation/pareto_pairwise figures
│   ├── generate_topology_maps_v3.py            Round 7: 4 real topologies on a world map (cartopy)
│   ├── generate_convergence_figure_v3.py       Round 7: HV/IGD vs. evaluations figure
│   ├── recompute_ablation_stats_v3.py          SUPERSEDED (kept for history only, see its
│   │                                             docstring) -- do not run against the current
│   │                                             30-seed ablation_v3_results.json
│   ├── hv_reference_sensitivity_v3.py          Round 6: HV re-computed at 3 reference points,
│   │                                             from cached per-run fronts (no rerun)
│   ├── analyze_effective_dimensionality.py     Round 6: PCA/correlation analysis of the
│   │                                             global front's 4 objectives
│   ├── run_sigma_ablation_v3.py                Round 8: sigma-optimised vs. 3 independently-tested
│   │                                             fixed values vs. K-averaged, genuinely n_var=10
│   │                                             for the fixed/averaged variants
│   ├── categorical_operators_baseline_v3.py    Round 7: category-native operators baseline
│   ├── multi_instance_oracle_v3.py             Round 7/8: parameterised oracle for city subsets
│   │                                             AND the genuine N=200 scale-up (verified
│   │                                             numerically identical to bdcz_oracle_v3.py on
│   │                                             the original 100-city set)
│   ├── run_multi_instance_v3.py                Round 7/8: 5-algorithm comparison on 3 subset
│   │                                             instances plus a genuinely independent N=200
│   ├── sequence_oracle_validation_v3.py        Round 7: directional cross-check of the fibre
│   │                                             transmittance-vs-distance law against SeQUeNCe
│   │                                             (external discrete-event simulator); optional
│   │                                             dependency, see the script's docstring
│   ├── sigma_ablation_stats_v3.py              Round 8: Kruskal-Wallis/Mann-Whitney/Holm/Cliff's
│   │                                             delta for the sigma ablation, from checkpoints
│   ├── categorical_operators_shared_ref_v3.py  Round 8: shared-reference HV comparison,
│   │                                             category-native-operator baseline vs. decoder
│   ├── multi_instance_pairwise_stats_v3.py     Round 8: NSGA-II-vs-NSGA-III pairwise HV/p-value
│   │                                             per additional instance, own 2-algorithm
│   │                                             reference front (persists numbers the
│   │                                             manuscript cites, previously ad hoc)
│   ├── metrics.py                              Round 8: canonical compute_spacing (numpy-only,
│   │                                             no pymoo), extracted from 4 duplicated copies
│   ├── coverage_threshold_sensitivity_v3.py    Round 8: C_tau recomputed at 5 fidelity
│   │                                             thresholds (0.50-0.90) for the pooled front
│   ├── routing_policy_sensitivity_v3.py        Round 8: objectives recomputed under
│   │                                             rate-optimal/latency-optimal routing (tests
│   │                                             whether low effective dimensionality is a
│   │                                             routing-policy artefact)
│   ├── random_search_baseline_v3.py            Round 8: i.i.d. uniform random sampling
│   │                                             sanity-check baseline, same bounds/budget/oracle
│   ├── run_modern_algorithms_v3.py             Round 8: SPEA2 + AGE-MOEA-II (2022) under the
│   │                                             exact main-comparison conditions
│   ├── representation_comparison_master_v3.py Round 8: decoder + relaxed-categorical +
│   │                                             category-native (d15, 30 seeds), ONE shared
│   │                                             reference front/normalisation
│   ├── unified_seven_algorithm_comparison_v3.py  Round 9: ALL SEVEN algorithms re-scored
│   │                                             against ONE reference front built from all 210
│   │                                             raw per-seed fronts; joint Kruskal-Wallis +
│   │                                             Holm over 3 metrics x 21 pairs = 63 tests,
│   │                                             plus the paired Friedman/Wilcoxon sensitivity
│   ├── generate_seven_algorithm_tables_v3.py  Round 9: emits the LaTeX table bodies for the
│   │                                             above, so the manuscript's numbers are
│   │                                             generated rather than transcribed
│   ├── connectivity_diagnostics_v3.py          Round 9: components, giant-component size and
│   │                                             reachable-pair fraction for all 659 front
│   │                                             points, on the adjacency the oracle routes over
│   ├── werner_path_validation_v3.py            Round 9: proves w=(4F-1)/3 makes the swap
│   │                                             composition an exact product, and validates the
│   │                                             oracle's path search against a -log w shortest
│   │                                             path and against exhaustive enumeration
│   └── decoder_edge_probability_v3.py          Round 9: effective per-pair edge-inclusion
│                                                 probability (mutual k-NN pairs get two draws,
│                                                 so 1-(1-p)^2), with a Monte-Carlo check
├── figures/
│   ├── algo_comparison_v3.pdf           boxplot + jittered per-seed points (7 algorithms, 30 seeds)
│   ├── ablation_v3.pdf                  boxplot + jittered per-seed points (8 variants, 30 seeds)
│   ├── pareto_pairwise_v3.pdf
│   ├── convergence_v3.pdf               Round 7: HV/IGD vs. evaluations, 5 algorithms
│   └── topology_maps_v3.pdf             Round 7: 4 representative topologies on a world map
├── results/
│   ├── baselines_comparison_v3.json         per-seed + aggregate metrics (5 algorithms, 30 seeds)
│   ├── ablation_v3_results.json             per-variant + per-seed metrics (8 variants, 30 seeds), incl. per-seed igd_union
│   ├── binary_baseline_v3_results.json      categorical baseline (own reference front, 30 seeds)
│   ├── binary_vs_decoder_shared_ref_v3.json decoder vs.\ categorical, shared reference front
│   ├── classical_baselines_v3_results.json  classical constructor metrics
│   ├── global_pareto_v3.npz                 global front: F, X, source_algorithm, source_seed, edge/medium counts
│   │                                          (F/X deduplicated by objective vector -- see Round 3 below)
│   ├── global_pareto_v3_summary.json        ranges, correlations, resource use, dedup counts
│   ├── global_pareto_v3_provenance.json     every (algorithm, seed) that produced each unique point
│   ├── unique_topologies_v3.npz              ALL structurally-distinct topologies (X, F, plus edge_u/edge_v/edge_type
│   │                                            + topology_offsets, flat CSR-style, allow_pickle=False) -- Round 4
│   ├── unique_topologies_v3_provenance.json  per-topology (algorithm, seed) provenance -- Round 4
│   ├── dominance_check_v3.json               classical-baseline dominance results
│   ├── stats_analysis_v3.json               Kruskal-Wallis/Mann-Whitney (primary) + Friedman/Wilcoxon (sensitivity)
│   ├── effective_dimensionality_v3.json      Round 6: PCA loadings + variance ratios
│   ├── hv_reference_sensitivity_v3.json      Round 6: HV ranking stability across reference points
│   ├── sigma_ablation_v3_results.json        Round 7: sigma-optimised vs. fixed vs. averaged
│   ├── sigma_ablation_v3_stats.json          Round 7: Kruskal-Wallis/Mann-Whitney for the sigma ablation
│   ├── categorical_operators_v3_results.json Round 7: category-native-operator baseline, 4 densities
│   ├── categorical_operators_vs_decoder_shared_ref_v3.json  Round 7: the manuscript's own
│   │                                          10-seed category-native-vs-decoder comparison
│   │                                          (Section: categorical baseline, ~6.4x gap)
│   │                                          still legitimately cites this file's d15 row as-is;
│   │                                          representation_comparison_master_v3.json (Round 8)
│   │                                          is the SEPARATE, 30-seed, master version of the
│   │                                          same comparison (~3.5x gap) -- both are cited for
│   │                                          different, clearly-scoped purposes, neither is stale
│   ├── convergence_v3_results.json           Round 7: per-generation HV/IGD, 5 algorithms
│   ├── multi_instance_v3_results.json        Round 7: 5-algorithm comparison, 3 additional real-city instances
│   ├── multi_instance_pairwise_stats_v3.json Round 8: NSGA-II-vs-NSGA-III pairwise HV/p-value per instance
│   │                                        (Round 9: Holm-corrected across all four instances,
│   │                                        with Cliff's delta and bootstrap CIs)
│   ├── unified_seven_algorithm_comparison_v3.json  Round 9: all 7 algorithms, ONE reference
│   │                                        front, 63-test joint Holm family, both designs
│   ├── connectivity_diagnostics_v3.json     Round 9: components/reach for all 659 front points
│   ├── connectivity_diagnostics_v3.npz      Round 9: the same, per point
│   ├── werner_path_validation_v3.json       Round 9: Werner-product identity + path-search checks
│   ├── decoder_edge_probability_v3.json     Round 9: effective per-pair edge probability
│   ├── GENERATED_SEVEN_ALGO_TABLES_v3.tex   Round 9: generated LaTeX table bodies
│   ├── topology_maps_v3_selection.json       Round 7: selection indices/objectives for the 4 mapped topologies
│   ├── sequence_oracle_validation_v3.json    Round 7: SeQUeNCe directional cross-check (rate vs. distance)
│   ├── coverage_threshold_sensitivity_v3.json  Round 8: C_tau ranking robustness, 5 thresholds
│   ├── routing_policy_sensitivity_v3.json    Round 8: PCA/correlations under 3 routing policies
│   ├── random_search_baseline_v3_results.json  Round 8: random-search baseline (own reference)
│   ├── random_search_vs_decoder_shared_ref_v3.json  Round 8: shared-reference vs. NSGA-II
│   ├── modern_algorithms_v3_results.json     Round 8: SPEA2 + AGE-MOEA-II, shared reference with Table 3
│   └── representation_comparison_master_v3.json  Round 8: decoder/relaxed-categorical/category-native, one shared reference
└── CHECKSUMS.sha256
```

## Reproducing the results

Requirements are pinned in `requirements.txt` (`pymoo==0.6.1.6`, `numpy`,
`scipy`, `networkx`, `matplotlib`, `cartopy`, `pytest`).

```bash
cd code/
python run_baselines_comparison_v3.py   # ~1-2 hours single-threaded, 150 runs (5 algorithms x 30
                                         # seeds; canonical oracle is exact, not a fast proxy).
                                         # Supports --algos/--seeds/--aggregate-only for sharding
                                         # across processes -- see the script's argparse help.
python run_ablation_v3.py               # ~3 hours single-threaded, 240 runs (8 variants x 30
                                         # seeds as of the round 6-followup expansion). Supports
                                         # --variants/--seeds/--aggregate-only for sharding.
python classical_baselines_v3.py        # <1 min, deterministic
python binary_baseline_v3.py            # ~3 hours single-threaded, 60 runs (2 conditions x 30
                                         # seeds as of the round 6-followup expansion). Supports
                                         # --algos/--seeds/--aggregate-only for sharding.
python build_global_pareto.py           # seconds, reads checkpoints above
python build_binary_vs_decoder_shared_ref_v3.py  # seconds, reads decoder checkpoints + binary results above
python stats_analysis_v3.py             # seconds, reads results JSON
python dominance_check_v3.py            # seconds
python generate_tables_v3.py            # seconds, emits GENERATED_TABLES_v3.tex
python generate_figures_v3.py           # seconds, regenerates algo_comparison/ablation/pareto_pairwise
python generate_topology_maps_v3.py     # seconds, reads unique_topologies_v3.npz (requires cartopy)
python generate_convergence_figure_v3.py  # seconds, reads convergence_v3_results.json
# recompute_ablation_stats_v3.py is SUPERSEDED -- do not run, see its docstring
python hv_reference_sensitivity_v3.py   # seconds, reads cached checkpoints
python analyze_effective_dimensionality.py  # seconds, reads global_pareto_v3.npz
python run_sigma_ablation_v3.py         # ~3.5 hours single-threaded, 75 runs (5 variants -- 1
                                         # optimised, 3 independently-tested fixed values, 1
                                         # K=3-averaged -- x 15 seeds). Genuinely n_var=10 for
                                         # the fixed/averaged variants (Round 8). Supports
                                         # --variants/--seeds/--aggregate-only.
python categorical_operators_baseline_v3.py  # ~5 hours single-threaded, 40 runs (4 densities x
                                         # 10 seeds) + 20 more seeds for d15 only (Round 8, to
                                         # reach 30 seeds for the master comparison below).
                                         # Supports --densities/--seeds/--aggregate-only.
python run_convergence_v3.py            # ~1 hour single-threaded, 50 runs (5 algorithms x 10
                                         # seeds, tracks best-so-far archive HV/IGD per
                                         # generation, Round 8). Supports
                                         # --algos/--seeds/--aggregate-only.
python run_multi_instance_v3.py         # ~1 hour for the 3 subset instances (150 runs, 3
                                         # instances x 5 algorithms x 10 seeds) plus ~11-23
                                         # hours single-threaded for global_n200 (25 runs, 5
                                         # algorithms x 5 seeds -- N=200 is much slower per
                                         # oracle call than N=100, Round 8). Supports
                                         # --instances/--algos/--seeds/--aggregate-only.
python sequence_oracle_validation_v3.py # ~10 sec, requires the separate `sequence` package
                                         # (pip install sequence==0.8.5); NOT required for any
                                         # other script in this package. Runs SeQUeNCe's own
                                         # discrete-event entanglement-generation protocol at 7
                                         # fibre distances x 10 seeds; see the script's module
                                         # docstring for what is and isn't cross-validated.
python sigma_ablation_stats_v3.py       # seconds, reads step4_sigma_v3/ checkpoints
python categorical_operators_shared_ref_v3.py  # seconds, reads step1/step5 checkpoints
python multi_instance_pairwise_stats_v3.py     # seconds, reads step7 checkpoints
python coverage_threshold_sensitivity_v3.py    # ~1 min, re-decodes the pooled front at 5 thresholds
python routing_policy_sensitivity_v3.py        # ~2 min, re-decodes the pooled front under 2 extra routing policies
python random_search_baseline_v3.py     # ~5.5 hours single-threaded, 30 seeds x 8,320 i.i.d.
                                         # random evaluations. Supports --seeds/--aggregate-only.
python run_modern_algorithms_v3.py      # ~1.5 hours single-threaded, 60 runs (SPEA2 + AGE-MOEA-II
                                         # x 30 seeds); AGE-MOEA-II requires the optional `numba`
                                         # dependency. Supports --algos/--seeds/--aggregate-only.
python representation_comparison_master_v3.py  # ~1 min, reads step1/step2/step5 checkpoints
python unified_seven_algorithm_comparison_v3.py  # ~2 min, re-scores all 7 algorithms on ONE
                                         # reference front; no new optimisation runs
python generate_seven_algorithm_tables_v3.py     # seconds, emits the manuscript's table bodies
python connectivity_diagnostics_v3.py    # ~15 min, re-decodes all 659 front points
python werner_path_validation_v3.py      # ~3 min, identity + oracle + exhaustive-enumeration checks
python decoder_edge_probability_v3.py    # ~1 min, Monte-Carlo over the real decode_topology()
```

Or run the whole pipeline's tests/figures/paper via the `Makefile`
(see `make test`/`make figures`/`make paper`/`make verify`).

The algorithm-comparison, ablation and categorical-baseline runs are all
seeded `42..71` (30 seeds each, as of the round 6-followup expansion). The
analytical oracle is fully deterministic given a topology, so the only
source of run-to-run variation is the evolutionary algorithm's own RNG
stream, fixed per seed. The oracle performs an exact Werner-optimal
widest-path search (`O(N)` Dijkstra-like searches per topology) rather
than a faster `-log F` shortest-path proxy, which is why runs are
markedly slower than a naive `floyd_warshall`-based implementation would
be — this is a deliberate choice, not an oversight (see
`bdcz_oracle_v3.py`'s module docstring).

## Provenance

This package has been through eight rounds of internal pre-submission
validation and reproducibility checks. See `CHANGELOG.md` for the full,
itemized audit trail (real code/data fixes verified against results
before being written down, not cosmetic text edits) -- kept out of
this README so it stays a normal reproducibility-package README
rather than a debugging narrative. Headlines: Round 4 found and fixed
a confirmed bug where HV/IGD/spacing/cardinality were computed on
non-deduplicated fronts, which reversed a reported finding (MOEA/D's
apparent spacing advantage over NSGA-II/NSGA-III was an artefact of
duplicate objective vectors and does not survive deduplication).
Round 6 expanded the algorithm comparison from 3 algorithms x 10 seeds
to 5 algorithms x 30 seeds (adding SMS-EMOA and RVEA), which both
resolved the original Wilcoxon p-value resolution floor and revealed
a genuine, seed-dependent bimodality in NSGA-II's hypervolume (7 of
30 seeds reach a qualitatively higher-HV region) that was invisible
at n=10. Round 7 fixed a real oracle bug (an infeasible-medium
sentinel leaking into a real, averaged latency in the categorical
baseline), corrected the spacing metric to Schott's canonical
Manhattan-distance/ddof=1 formula, made independent-samples testing
(Kruskal-Wallis/Mann-Whitney) the primary statistical analysis, and
added convergence curves, a sigma ablation, a category-native-operator
baseline, topology maps, multi-instance validation across three
additional real-city benchmarks, and a directional cross-check of the
oracle's fibre transmittance-vs-distance law against SeQUeNCe, an
independent discrete-event quantum network simulator (Spearman
$\rho=1.0$ across seven distances; see `sequence_oracle_validation_v3.py`
for the honest scope of what this does and does not validate).
Round 8, prompted by an external review requesting a formal physics
audit before any further work, found and fixed a real, confirmed bug
in the oracle's fibre latency formula (`PHYSICS_ORACLE_AUDIT.md` has
the full derivation and citations): a spurious factor made latency
scale quadratically rather than linearly with distance. Since latency
is a directly-optimised objective, every evolutionary campaign in
this package was fully recomputed under the corrected oracle (735
runs). The most consequential downstream change: the ablation's
headline finding reverses from "fibre and satellite are both
required" to "fibre is necessary, but satellite's contribution beyond
fibre alone is not statistically detectable" -- reported directly,
not smoothed over (`CHANGELOG.md` item 75 has the full list of
consequences, including a NSGA-II-vs-NSGA-III instance-dependent
reversal on the smallest additional multi-instance benchmark and a
new disagreement between the primary and paired statistical designs
on that same comparison). The identical latency bug was confirmed
present in the companion QI-Bench-100 paper's own canonical oracle
copy, but was explicitly left untouched there, out of scope for this
package.

## License

Code (`code/`) is released under the MIT license (`LICENSE-CODE`) --
an independent pre-submission technical audit flagged that CC-BY-4.0,
a content/data license, is not an appropriate license for software (it
lacks the patent and liability terms a software license needs). Data (`cities.csv`, and the
result CSV/JSON/NPZ files under `results/`) remains released under
CC-BY-4.0 (`LICENSE-DATA`), for which CC-BY is the standard choice.
