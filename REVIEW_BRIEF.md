# Review Brief — QI-Bench-Pareto-v1

Short orientation for a reviewer or editor, distinct from `README.md`
(code/reproducibility documentation) and `CHANGELOG.md` (full itemised
round-by-round audit trail). Read this first; follow the links below for
depth.

## What this package supports

*"Many-Objective Optimisation of Global Quantum Internet Topologies: A
Compact Parametric Representation with an Analytical BDCZ Oracle"*
(Jesús Gil Ruiz, Diego Rubén Rodríguez Regadera and Rafael Muñoz Gil; accepted for publication in *Entropy*, MDPI, 2026).

The paper designs global quantum-internet topologies (fibre,
low-Earth-orbit satellite, free-space optical) over 100 cities
(C(100,2) = 4,950 candidate edges), jointly optimising fidelity,
entanglement rate, latency and coverage. Direct per-edge categorical
optimisation over 4,950 edges is shown to degenerate; the paper's core
contribution is an 11-variable parametric representation that
deterministically decodes into a topology via $k$-nearest-neighbour
rules, plus a benchmark (seven many-objective algorithms, 30 seeds each;
an eight-variant decoder ablation; classical network-design baselines).

## Headline, verifiable results

All numbers below are read directly from files in `results/`, not
transcribed from memory — see the "Verify a number" column. Current as
of the pooled-front unification fix (CHANGELOG Round 11, package Round
14); every number below post-dates that fix.

| Claim | Value | Verify against |
|---|---|---|
| Top group (mean HV, 30 seeds) | NSGA-II $0.353\pm0.040$, SPEA2 $0.352$, AGE-MOEA-II $0.347$ — mutually indistinguishable, each significantly ahead of NSGA-III/SMS-EMOA/RVEA/MOEA-D | `unified_seven_algorithm_comparison_v3.json` |
| NSGA-II bimodality | 6/30 seeds in a higher-HV cluster (NSGA-III: same, 6/30, 4 shared seeds) | `unified_seven_algorithm_comparison_v3.json` → `summary.nsga2.hv_raw` |
| Fibre necessary under NSGA-III, satellite not detectably so beyond it | fibre-free variants collapse HV significantly (p_Holm<1.4e-8); fibre-only/fibre+satellite indistinguishable from the full decoder (p_Holm=0.088/1.000). Run with NSGA-III only -- now the algorithm separated from the top group -- not confirmed for NSGA-II/SPEA2/AGE-MOEA-II | `ablation_v3_results.json`, `stats_analysis_v3.json` |
| Pooled evolutionary front size | 561 unique objective vectors (all 7 algorithms), 944 topologies | `global_pareto_v3_summary.json`, `unique_topologies_v3.npz` |
| Front connectivity | 0% of 561 points connect all 100 cities; median reaches 50.9% of pairs | `connectivity_diagnostics_v3.json` |
| Classical constructors non-dominated | 1 of 14 (RGG at 5,000 km) — was 3 (hub-and-spoke) before the satellite fix; re-checked against the 561-point front | `dominance_check_v3.json` |
| Effective dimensionality | PC1 = 96.1%, PC1+PC2 = 99.1% | `effective_dimensionality_v3.json` |
| HV ranking stable across reference points | Spearman 0.983 (1.05), 0.987 (1.20) vs. 1.10 | `hv_reference_sensitivity_v3.json` |
| SeQUeNCe directional cross-check (rate vs. distance) | Spearman $\rho=1.0$, 7 distances (fibre model bit-identical pre/post satellite fix) | `sequence_oracle_validation_v3.json` |
| Random search vs. NSGA-II | NSGA-II wins but modestly (+8.8%, p=5.0e-4) | `random_search_vs_decoder_shared_ref_v3.json` |
| Genuinely independent N=200 instance (100 new real cities, GeoNames) | NSGA-III leads NSGA-II (point estimate), not significant, now at n=10 seeds; tested with the original 5 algorithms only | `multi_instance_v3_results.json` → `per_instance.global_n200` |
| Satellite line-of-sight bug: found, fixed, everything recomputed | 14.6% of satellite edges had no line of sight before the fix; 925 evolutionary runs recomputed | `CHANGELOG.md` items 91-105 |
| Pooled-front unification bug: found, fixed, every consumer re-run | Released front was still 5-algorithm/513-point after the main comparison moved to 7 algorithms/561 points; root-caused to `build_global_pareto.py` and fixed | `CHANGELOG.md` items 106-121 |

## Package map

```
QI-Bench-Pareto-v1/
├── code/       all analysis/experiment scripts (see README.md for the full list + run order)
├── results/    every JSON/NPZ underlying every number in the paper
├── figures/    the 3 figures as shipped in the manuscript
├── cities.csv  the 100-city benchmark (reverse-geocoded, city/country/region)
└── CHECKSUMS.sha256  sha256 of every file above, for integrity verification
```

Reproduction commands and expected runtimes are in `README.md`. All
evolutionary runs are seeded and the analytical oracle is deterministic
given a topology, so results are exactly reproducible given the same
`pymoo==0.6.1.6` version pinned in `requirements.txt`.

## What changed most recently (Round 11 in this file's numbering,
"Round 14" in the package's download-numbering)

An independent external review of the Round 13 PDF found that the
satellite-visibility fix below (Round 10/13) had been applied to the
manuscript's headline algorithm comparison (which moved to all seven
algorithms) but not propagated to the released pooled Pareto front:
`build_global_pareto.py`, the source of `global_pareto_v3.npz` and
everything computed from it, still hardcoded the original five
algorithms. Root-caused to that single script, fixed, and every
downstream consumer (connectivity, PCA, threshold/routing sensitivity,
classical-baseline dominance, the topology-map figure, the released
package itself) re-run against the corrected 561-point, seven-algorithm
front. Headline conclusions survive (RGG at 5,000 km remains the sole
non-dominated classical constructor; 0% of front points fully connect
all 100 cities), but every number computed from the front shifted
slightly (PC1 95.6%→96.1%, median reachable-pair fraction 54.8%→50.9%,
etc. -- see the table above). The same review also found and the fixes
included: a second, independent bug in the pairwise-projection figure
that silently dropped SPEA2/AGE-MOEA-II points from the legend; a
raw-vs-effective edge-count mismatch across Figure 4's panels; a
mis-scaled convergence figure (Figure 2) still normalised against the
old five-algorithm comparison; an inverted statistical claim (NSGA-III
vs. SMS-EMOA spacing); a stale threats-to-validity paragraph still
describing the pre-recomputation statistics; and a dozen smaller
terminology/cross-reference/number-drift issues throughout the prose.
Full itemised list: `CHANGELOG.md`, "Round 11" section, items 106-121.

## What changed before that (Round 10 in this file's numbering,
"Round 13" in the package's download-numbering)

An independent external review found a critical, confirmed-real physics
bug: the satellite link model computed slant range with a flat-Earth
approximation and had NO line-of-sight visibility limit, while the
decoder's own search bounds allowed satellite links out to 15,000 km. A
single satellite at the modelled 500 km altitude can see both ground
stations only up to 4,891 km apart; beyond that the Earth blocks the
line of sight. On the previously released front, 14.6% of all satellite
edges (69,800 of 476,827) were physically impossible. Two smaller,
related bugs in the same link-budget calculation were found and fixed
alongside it (an inconsistent two-arm efficiency, and a downlink
latency that summed both arms' transit times instead of using the
single simultaneous transit). All three are fixed, and **every result
in this package was recomputed from scratch**: 925 evolutionary runs
across all nine experimental blocks, plus every derived statistic and
figure. The headline result changed in the process: NSGA-II vs
NSGA-III, not significant in the previous round, is now significant
under both a 5- and a 7-algorithm Holm family; the classical-baseline
non-dominated exemplar changed from hub-and-spoke to a random geometric
graph; and no point on the released front now connects all 100 cities
(was 28.1% before the fix). Full itemised list: `CHANGELOG.md`,
"Round 10" section, items 91-105.

## What changed before that (Round 8)

An external review of the Round-10 PDF produced findings across the
oracle's mathematics, the decoder's specification, the statistical
framing and the bibliography. Every one was verified against the real
files before being acted on; those that turned out to be false
positives are listed as discarded at the end of `CHANGELOG.md`'s
Round 8 section. Four changes are substantive enough that a reader of
an earlier round needs to know about them:

1. **The paper's headline algorithm result changed.** SPEA2 and
   AGE-MOEA-II were previously reported in a side paragraph, scored
   against a wider bounding box than Table 3 while being described as
   directly comparable to it. All seven algorithms are now re-scored
   against one reference front built from all 210 raw per-seed fronts
   (`unified_seven_algorithm_comparison_v3.py`), with a single joint
   Holm family of 3 metrics x 21 pairs = 63 tests. The raw p-values
   did not change, but at 63 tests instead of 30, NSGA-II vs NSGA-III
   is no longer significant on any metric. The result is now stated
   as a top group -- NSGA-II, SPEA2 and AGE-MOEA-II mutually
   indistinguishable -- rather than as a single winner.
2. **The Werner path search's justification was wrong and is now a
   proof.** Under `w=(4F-1)/3` the swap composition is an exact
   product, so an exact additive shortest-path formulation exists and
   optimality follows from the standard multiplicative-shortest-path
   argument, not from the "monotonically non-increasing" hand-wave the
   manuscript used. `werner_path_validation_v3.py` verifies the
   identity, the oracle against an independent `-log w` shortest path
   on 18,674 real pairs, and both against exhaustive enumeration. It
   also shows the `-log F` proxy the paper declines to use is wrong
   for 18.0% of reachable pairs, not "a small fraction".
3. **The decoder's effective edge probability is not its nominal
   one.** A pair in both endpoints' k-NN lists is drawn twice, so its
   inclusion probability is `1-(1-p)^2`; 48.7-55.4% of candidate
   fibre pairs are mutual k-NN. Verified by Monte-Carlo against the
   real `decode_topology()` (`decoder_edge_probability_v3.py`). No
   result changes; Algorithm 1 now states it.
4. **Most of the released front is not a connected network.** Only
   28.1% of the 659 points connect all 100 cities; the median point
   has 10 components and reaches 27.1% of city pairs
   (`connectivity_diagnostics_v3.py`). Reported in a new subsection,
   and the front should be filtered by connectivity or minimum
   coverage before use.

Also: the four multi-instance NSGA-II/NSGA-III tests are now
Holm-corrected as one family (which demotes `regional_eu_na` from
significant to suggestive), a complete link-level oracle appendix was
added so the physics can be audited without the code, four reference
metadata errors were fixed, and float placement and figure font sizes
were corrected. Manuscript is now 85 pages.

## What changed in Round 8

An external review asked for a formal physics audit of the oracle's
latency formulas, with citations, before any further changes. That
audit (`PHYSICS_ORACLE_AUDIT.md`) found a real, confirmed bug: fibre
latency scaled quadratically rather than linearly with distance
(verified two independent ways -- a naive sequential-repeater model,
and Sangouard et al., *Rev. Mod. Phys.* 83, 33 (2011), Eq. 12 --
neither supports latency increasing with the number of repeater
spans at fixed total distance; the manuscript's own "sum of per-link
propagation times" definition also contradicted the quadratic code).
Fixed, and **all 735 evolutionary runs across every campaign were
fully recomputed** under the corrected oracle, since latency is a
directly-optimised objective. The most consequential result: the
ablation's headline finding reverses from "fibre and satellite are
both required" to "fibre is necessary, satellite's contribution
beyond fibre alone is not statistically detectable" -- reported
directly, not smoothed over. Other real, verified changes: the
common-seed paired statistical design now disagrees with the primary
independent-samples design on the NSGA-II-vs-NSGA-III comparison
itself (7 of 30 total comparisons disagree, up from 0); NSGA-III's
own seed-bimodality grew from 3/30 to 10/30 high-HV seeds, with
majority overlap with NSGA-II's; the `metro_europe` multi-instance
result now reverses in NSGA-III's favour (not significant either
direction). See `CHANGELOG.md` item 75 for the complete, itemised
list of consequences. The identical bug is confirmed present in the
companion QI-Bench-100 paper's own canonical oracle but was
deliberately left untouched there.

## What changed in Round 7

Full external-review author-approved scope: corrected the spacing metric
to Schott's canonical Manhattan-distance/ddof=1 formula (`ablation_v3_results.json`,
`baselines_comparison_v3.json`); made independent-samples testing
(Kruskal-Wallis/Mann-Whitney, `stats_analysis_v3.json`) the primary
statistical analysis, with the prior Friedman/Wilcoxon paired design
retained as a common-seed sensitivity check; added a sigma ablation
(`sigma_ablation_v3_results.json`), a category-native-operator baseline
(`categorical_operators_v3_results.json`), convergence curves
(`convergence_v3_results.json`), topology maps (`topology_maps_v3_selection.json`),
multi-instance validation across 3 additional real-city benchmarks
(`multi_instance_v3_results.json`), and a directional cross-check of the
oracle's fibre transmittance-vs-distance law against SeQUeNCe, an
independent discrete-event quantum network simulator (Spearman
$\rho=1.0$ across 7 distances; `sequence_oracle_validation_v3.json` --
validates the rate dimension only, explicitly does not validate
fidelity, see the script's docstring). Also converted `unique_topologies_v3.npz`
to a pickle-free flat/offsets format and added a 15-test pytest suite +
`Makefile`. Full itemised list: `CHANGELOG.md`, Round 7 section.

## What changed in Round 6

Two independent pre-submission technical audits (received together,
AI-assisted and reviewed by the author -- not formal journal peer
review) both proposed
expanding the original 3-algorithm/10-seed design; the combined,
most ambitious scope was adopted: 5 algorithms (added SMS-EMOA, RVEA)
x 30 seeds, plus a symmetric IGD variant, an HV reference-point
sensitivity check, and a real (not deferred) effective-dimensionality
analysis. Full itemised list: `CHANGELOG.md`, Round 6 section.

## Known limitations (not hidden — see the manuscript's Threats to
## validity, Section 6, and `CHANGELOG.md` item 28 for the full list)

None of the five algorithms had their operators/hyperparameters tuned
for this specific problem (matched, untuned settings throughout); the
NSGA-II bimodality's root cause is reported but not investigated
further; the ablation's per-seed symmetric IGD would need a fresh
80-run campaign not performed this round (a pooled, not per-seed,
value is reported instead); FSO's marginal contribution is scoped to
this benchmark's sparse, intercontinental city distribution and its
150 km range, not a general claim about FSO's importance at every
network scale.
