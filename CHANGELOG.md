# Changelog

Full audit trail of what was found and fixed across rounds of independent
external review, moved here from README.md to keep that file a normal
package README rather than a debugging narrative. All entries below are
real code/data fixes, verified against the actual results before being
written down, not cosmetic text edits -- see the manuscript's own
Reproducibility Statement for what each script produces.

**Round 1** (12-D → 11-D, fairness fixes):
1. `hub_boost` dead parameter removed. It was computed and bounds-checked
   in the decoder but never referenced in edge-sampling logic, reducing
   the representation from 12 to 11 variables.
2. Heterogeneous medium assignment added to the binary and classical
   baselines (originally hardcoded to fibre), which had given the
   parametric decoder an unearned advantage from medium choice alone.

**Round 2** (canonical-oracle unification — the more consequential fix):
3. **Oracle unified with the companion QI-Bench-100 paper's canonical
   implementation.** An earlier version of this package used a `-log(F)`
   shortest-path proxy oracle with an undisclosed fidelity floor
   (`max(LINK_FID, 0.5)`) and miscalibrated satellite/FSO base fidelities
   (`SAT_BASE_FID=0.97`, `FSO_BASE_FID=0.98`, no satellite distance-decay
   term) that diverged from the companion paper's real, submitted oracle
   (`SAT_BASE_FID=0.83` with distance-dependent decay, `FSO_BASE_FID=0.92`,
   exact Werner-optimal path search, no floor). `bdcz_oracle_v3.py` is the
   companion paper's oracle verbatim. **This changed results
   materially**: e.g. the decoder ablation's finding under the earlier,
   miscalibrated oracle (satellite-only outperforming the full decoder)
   reversed once seeds were increased to 10 and the canonical oracle was
   used — the corrected result is that fibre and satellite are jointly
   necessary, matching the paper's original hypothesis but now with
   adequate statistical power and a validated physical model.
4. Binary baseline made fully categorical (each of the 4,950 city pairs
   is a free `{no edge, fibre, satellite, FSO}` choice via a real-valued
   relaxation) rather than assigning each active edge a fixed
   argmax-fidelity medium, which was still a representationally unfair
   comparison against the freely-searching parametric decoder.
5. Ablation extended to 8 variants (added FSO-only) and 10 seeds (was 5),
   and the algorithm-comparison / ablation Holm correction is now applied
   jointly across all nine (or seven) pairwise tests rather than
   per-metric-family.
6. Table 2 (decoder parameter ranges), Algorithm 1 (decoder pseudocode)
   and Equation 2 (rate aggregation) were all previously inconsistent
   with the actual code; all three now match `bdcz_oracle_v3.py` exactly.

**Round 3** (three confirmed real bugs, found by independent review of
`Paper2_SUBMISSION_READY.zip` -- all verified against code/data before
fixing, all reran, not just text edits):
7. **The ablation did not actually disable a medium.** `decode_topology`
   floors `k_fibre`/`k_sat` at 1 and `p_fibre`/`p_sat`/`p_fso` at 0.05
   regardless of the parameter vector, so setting an ablation variant's
   params to 0.0 (the old approach) did not zero out that medium --
   e.g. "satellite-only" could still emit fibre edges. Fixed with
   explicit `enable_fibre`/`enable_satellite`/`enable_fso` flags that
   hard-disable a medium's entire edge-generation step; verified these
   flags do not change Full Decoder behaviour when all three are `True`
   (the default, used by every other experiment in this package).
   Ablation statistics also changed from Kruskal-Wallis/Mann-Whitney
   (independent-samples tests) to Friedman/paired-Wilcoxon (the 8
   variants share the same 10 seeds, so the data is paired).
8. **MOEA/D ran on ~19% fewer function evaluations than NSGA-II/NSGA-III.**
   `pymoo`'s MOEA/D ties population size to `len(ref_dirs)` (verified:
   84, not the `104` the manuscript claimed for "all algorithms"), while
   NSGA-II/NSGA-III use an independently-set `pop_size=104`. Fixed by
   giving MOEA/D more generations (99 instead of 80) so its total
   evaluation count (8,316) matches the other two (8,320) to within
   0.05%, rather than forcing an equal population size onto an algorithm
   whose population size is architecturally fixed to its reference
   directions.
9. **The "1,001-point global Pareto front" was ~53% duplicates.**
   Verified directly against the raw `.npz`: only 494 of the pre-fix
   1,001 raw records were unique objective vectors, and MOEA/D's 278
   raw records collapsed to just 19 unique vectors -- so "MOEA/D
   contributes 28% of the front" counted the same handful of solutions
   rediscovered across seeds as if they were independent contributions.
   Fixed in `build_global_pareto.py`: the front is now deduplicated by
   objective vector (452 unique points as of Round 3, before the
   Round-4 rerun below changed the raw counts feeding into it), with
   full per-point (algorithm, seed) provenance preserved in
   `global_pareto_v3_provenance.json` rather than discarded, and the
   number of structurally distinct decoded topologies (599 as of Round
   3) reported separately from the objective-vector count, since a
   Pareto *front* (objective space) and Pareto *set* (design space) are
   related but distinct quantities.

Also fixed in Round 3: FSO search bounds were `[10,200]` km while the
physical model caps FSO at 150km (`fso_link` returns zero fidelity beyond
that), letting the optimiser search distances the oracle itself treats as
invalid -- bounded to `[10,150]` to match. Software versions in the
manuscript (`NumPy 1.26`, `SciPy 1.13`) did not match the actual
environment that produced every result (`NumPy 2.4.3`, `SciPy 1.17.1`,
verified directly) -- corrected. `generate_figures_v3.py`'s output path
resolved one directory level too shallow (a real, independent bug caught
while regenerating figures for this round, unrelated to the review) --
fixed.

**Round 4** (a fourth external review, auditing the Round-3 package
end-to-end): one leftover stale number (Contributions still said "795",
missed in the Round-3 sweep) plus one further confirmed, material bug
found by directly testing the code rather than just re-reading the
manuscript:

10. **HV/IGD/spacing/|F\*| were computed on the RAW, non-deduplicated
    per-run front, not the deduplicated one.** Verified directly by
    running a single fresh MOEA/D seed: its own non-dominated
    population (80 points) contained only 11 unique objective vectors
    (86% exact duplicates) -- confirming this was not a hypothetical
    concern. HV is a set measure and is mathematically invariant to
    exact duplicates (confirmed empirically: identical to the 4th
    decimal before and after this fix), but IGD, spacing and |F\*|
    are not. Fixed with a new shared `dedup_front()` helper
    (`bdcz_oracle_v3.py`) applied immediately after non-dominated
    filtering in `run_baselines_comparison_v3.py`, `run_ablation_v3.py`
    and `binary_baseline_v3.py`, before any metric is computed or a
    front is persisted. **This reversed a real finding**: MOEA/D had
    appeared to have significantly better (lower) spacing than
    NSGA-II/NSGA-III ($p_{\mathrm{Holm}}\leq0.039$) -- after
    deduplication, that advantage disappears entirely (all pairwise
    $p_{\mathrm{Holm}}=1.000$, Friedman $p=0.670$), and MOEA/D's true
    unique-solution count falls from a reported $74.6\pm8.3$ to
    $9.8\pm2.3$. HV rankings (NSGA-II > NSGA-III > MOEA/D) and the
    ablation's HV-based conclusions ("fibre and satellite jointly
    necessary") are unaffected, since HV does not depend on
    deduplication.

Also in Round 4: the structurally-distinct topologies were previously
only *counted*, not actually released (`global_pareto_v3.npz` held 452
representative decision vectors, one per unique objective vector, not all
of the topologies that map to them) -- now released in full as
`unique_topologies_v3.npz` (551 topologies after the Round-4 rerun, with
edge lists, media and per-point provenance in
`unique_topologies_v3_provenance.json`). The categorical baseline was
extended from 3 to 10 seeds to match the decoder experiments' statistical
power, and its own shared-reference comparison against the decoder
(`binary_vs_decoder_shared_ref_v3.json`) was found to be an orphaned
artefact dated three rounds earlier that no script actually regenerated
(still showed the pre-dedup `npar=104±0` convention) -- rebuilt as
`build_binary_vs_decoder_shared_ref_v3.py`, a real, reproducible script.
Three bibliography entries had real errors, verified and corrected
against the actual publication record (RELiQ: journal/note concatenation
bug plus missing volume/pages; Helsen & Wehner: wrong year, 2022 instead
of the actual 2023 publication; Pouryousef et al.: wrong venue -- listed
as a journal, actually IEEE INFOCOM 2023 -- and an incomplete "and
others" author list). The abstract exceeded the journal's 250-word limit
(272 words); trimmed to 229. Resource-use fractions were reported as a
single *pooled* ratio (sum-of-medium-edges / sum-of-all-edges across the
whole front, which implicitly weights larger topologies more) under an
ambiguous "mean fraction" label; both the pooled figure and the mean/std
of each topology's own fraction (which differ materially: fibre 5.4%
pooled vs. 12.8% per-topology) are now reported, clearly labelled. `|E|`
was ambiguous across three different edge-counting conventions (raw
medium-specific, unique city pairs, effective oracle-routed edges, which
can differ when a pair carries more than one medium); all three are now
reported. Figure 2 (pairwise Pareto projections) coloured every point by
its first-discovering algorithm even for the ~5% of points independently
rediscovered by more than one algorithm; a "Shared" category now exists.
The decoder's seed offset ($\sigma$) being optimised as an eleventh
decision variable -- a PRNG seed, not a physical parameter -- is now
flagged explicitly as a threat to validity. The AI-use declaration now
mentions both tools actually used substantively (Claude and ChatGPT), not
only one.

**Round 5** (post-Round-4 audit, two independent reviewers -- one found a
typesetting regression the Round-4 edits introduced, the other found real
transcription errors in the manuscript text):

11. **Two IGD/timing/p-value transcription errors left over from
    Round 4's rewrite.** The Results paragraph and Conclusion still
    quoted NSGA-II's pre-dedup mean IGD ($0.033$) instead of the
    correct, deduplicated value ($0.054$); the joint-Holm significance
    threshold was written as $p_{\mathrm{Holm}}\leq0.029$ in three
    places when the actual value for the NSGA-II vs.\ NSGA-III HV test
    is $p_{\mathrm{Holm}}=0.039$ (verified directly against
    `stats_analysis_v3.json`); and the Conclusion's per-algorithm mean
    runtimes ($6.6$/$5.3$/$3.6$ min) were the Round-3 numbers, not the
    Round-4 ones ($8.1$/$5.6$/$3.8$ min). All three are hand-transcription
    errors (copying a number by re-typing it rather than reading it
    fresh from the JSON each time), not computation bugs -- the
    underlying results were always correct.
12. **A table added in Round 4 (`Global Pareto front (range)` row in
    the classical-baselines table) overflowed the page margin by
    ~62pt.** Wrapped in `\resizebox{\columnwidth}{!}{...}`, matching
    every other wide table in the manuscript.
13. **`CHECKSUMS.sha256` was regenerated before, not after, a final
    `__pycache__` cleanup**, so it retained one stale entry for a
    deleted `.pyc` file and `sha256sum -c` reported `FAILED`. Fixed by
    regenerating checksums as the last step before packaging, and this
    is now the documented order.
14. Corresponding-author email unified to
    `jesus.gil@universidadeuropea.es` across this paper and its
    companions (Paper 2 had drifted to the shorter
    `jesus.gil@universidadeuropea.es`, an Editorial-Manager login
    address, not the address on file for correspondence).
15. This `CHANGELOG.md` split out of `README.md`, so the README reads as
    a normal reproducibility-package README rather than a debugging
    narrative (a change two independent reviewers asked for across
    Rounds 3 and 4).

**Round 6** (two independent external reviews received together, partially
overlapping, partially divergent in scope -- one recommended only expanding
to $20$--$30$ seeds; the other recommended $30$ seeds \emph{and} two
additional algorithms. The most ambitious combined scope was adopted:
$30$ seeds and both SMS-EMOA and RVEA):

16. **Algorithm comparison expanded from NSGA-II/NSGA-III/MOEA-D x $10$
    seeds to NSGA-II/NSGA-III/MOEA-D/SMS-EMOA/RVEA x $30$ seeds** (150
    fresh optimisation runs; the canonical package's checkpoints had been
    deleted after the Round-3 aggregation, so all 150 combinations were
    re-run from scratch, sharded 20-way across processes with
    checkpoint-based resumability). AGE-MOEA-II was considered and
    rejected: it requires `numba`, a dependency the project otherwise has
    no need for, and RVEA already covers the "modern many-objective
    algorithm" requirement the reviews asked for. Both new algorithms'
    per-generation evaluation budget was calibrated directly against the
    real 11-D oracle (not assumed from `pymoo`'s docs or a toy problem):
    both consume exactly $104$ evaluations/generation with zero intercept,
    so neither needed the generation-count adjustment MOEA/D requires.
17. **The joint-Holm correction now spans $30$ pairwise tests** ($3$
    metrics x $\binom{5}{2}=10$ algorithm pairs), up from $9$. This also
    resolves the Wilcoxon signed-rank test's own resolution floor
    ($2/2^{30}\approx1.86\times10^{-9}$ at $n=30$, versus
    $2/2^{10}=0.00195$ at the old $n=10$ -- a genuine statistical-power
    concern raised independently by both reviews, since several reported
    "significant" $p$-values at $n=10$ sat close to that floor).
18. **Added an independent-samples sensitivity check** (Kruskal-Wallis
    omnibus + Holm-corrected Mann-Whitney pairwise) alongside the primary
    paired Friedman/Wilcoxon design, in `stats_analysis_v3.py`. The paired
    design (blocking by seed index) is standard, correct EC statistical
    methodology and was kept as primary rather than replaced; the
    independent-samples check agrees with it on all $30/30$ tested
    comparisons, which is reported as evidence the pairing assumption is
    not driving the significance pattern.
19. **`ms()`'s standard deviation fixed from `ddof=0` (population std) to
    `ddof=1` (sample std)** in both `run_baselines_comparison_v3.py` and
    `run_ablation_v3.py`. For the $n=10$ ablation, `ddof=1` stds are
    $\sim5.4\%$ larger than the previously-reported `ddof=0` values.
20. **A structural bias in the ablation's `igd` metric was found and made
    transparent, not silently fixed by replacement.** The existing `igd`
    is computed only against variant A (full decoder)'s own pooled front,
    which trivially makes A's own IGD close to zero -- not a symmetric
    comparison across variants. Added a second, genuinely symmetric
    `igd_union`, computed against the deduplicated union of all eight
    variants' pooled fronts (`recompute_ablation_stats_v3.py`). This is
    reported as one pooled value per variant (against the shared
    reference), not a per-seed mean$\pm$std, because the raw per-seed `F`
    arrays needed for a true per-seed `igd_union` were not retained after
    the original ablation run and a full $80$-run rerun was judged not
    worth the compute cost for this secondary statistic this round -- an
    explicit scope decision, not a silent gap. Under this symmetric
    metric, variant H (satellite+FSO) is actually closest to the shared
    union front ($\mathrm{IGD}_{\mathrm{union}}=0.145$), not variant A
    or C, illustrating that HV and IGD$_{\mathrm{union}}$ can disagree on
    ranking without either being wrong (Section~5.4 of `main.tex`).
21. **Recomputed HV at two additional reference points** ($1.05$ and
    $1.20$, alongside the main text's $1.10$) directly from the cached
    per-run fronts, with no new optimisation. The five-algorithm ranking
    is unchanged at both alternative points; per-seed Spearman rank
    correlation against the $1.10$ ranking is $0.983$ ($1.05$) and
    $0.993$ ($1.20$) (`hv_reference_sensitivity_v3.py`).
22. **Added a real, computed effective-dimensionality analysis**
    (`analyze_effective_dimensionality.py`): standardised PCA over the
    global front's four natural-sign-corrected objectives. PC1 explains
    $91.3\%$ of variance, PC1+PC2 $97.4\%$. This replaces a prior draft's
    placeholder statement that deferred the analysis to future work
    without having run it.
23. **Figure~1 (algorithm comparison) changed from a bar+error-bar plot to
    a boxplot with jittered raw per-seed points**, so an individual
    high-HV outlier run is visible as a specific labelled point rather
    than folded invisibly into a mean+std bar -- directly addressing both
    reviews' concern that a bar+errorbar plot cannot distinguish "a
    genuinely elevated, consistent mean" from "one lucky outlier run."
    This is also how the NSGA-II bimodality (below) was first noticed.
24. **NSGA-II's mean hypervolume is bimodal across its 30 seeds, not
    uniformly elevated -- confirmed as a genuine, recurring phenomenon,
    not noise or a bug.** $7$ of $30$ seeds ($43,54,55,58,63,65,71$) reach
    a qualitatively higher-HV cluster ($0.298$--$0.348$) than the
    remaining $23$ ($0.205$--$0.216$, gap $=0.082$), with identical
    `n_calls=8320` for every seed (ruling out a budget artefact) and
    consistently better IGD ($0.02$--$0.03$ vs.\ $0.10$--$0.11$) for the
    high-HV cluster (ruling out a measurement artefact). NSGA-III shows a
    smaller analogous pattern ($3/30$ seeds). Reported directly in the
    manuscript (Section~6, "NSGA-II's mean hypervolume is bimodal across
    seeds, not uniformly elevated") rather than smoothed into a single
    mean, with root-cause investigation explicitly left as future work.
25. **A real, pre-existing inconsistency was found and fixed**: the
    manuscript's prose claimed "$14$ classical constructors, $11$
    dominated" but the classical-baselines table only listed $9$ rows.
    Fixed by expanding the table to the actual $14$ rows (using
    `classical_baselines_v3.py`'s real output), not by shrinking the
    prose claim to match the smaller table.
26. **Global Pareto front recomputed for 5 algorithms x 30 seeds**:
    $846$ deduplicated unique objective vectors (was $452$), $1{,}347$
    structurally-distinct topologies (was $551$).
27. **Packaging/metadata hygiene**: `hypertexnames=false` added to the
    `hyperref` options, eliminating all $13$ "duplicate destination"
    pdfTeX warnings (root cause not fully isolated -- no duplicate
    `\label`s were found -- but the option fully eliminates the symptom).
    `cities.csv` rebuilt with proper `city`/`country` columns (previously
    only `index,latitude,longitude,region`), reverse-geocoded and then
    manually corrected for $15$ entries where the geocoder returned
    overly-granular district names for well-known major cities.
    `LICENSE-CODE` (MIT) and `LICENSE-DATA` (CC-BY-4.0) split out of a
    single CC-BY-4.0 license, since CC-BY is a content/data license and
    is not an appropriate license for software (it lacks patent/liability
    terms). `CITATION.cff` and `.zenodo.json` added with real
    author/ORCID/affiliation metadata; the repository URL and journal DOI
    are left as explicit placeholders pending the author's own Zenodo
    deposit and the journal's DOI assignment, not fabricated. Removed
    $48$ unused `refs.bib` entries (verified by parsing cited keys vs.\
    bib keys) and fixed one entry's incorrect `@article` type
    (`bezerra2017empirical`, actually `@inproceedings`).
28. **Explicitly deferred, not silently dropped**: convergence curves;
    a categorical-baseline operator/density-sensitivity upgrade; a
    dedicated ablation of the decoder's seed offset as a decision
    variable; routing-policy sensitivity analysis; topology world maps;
    further coverage-threshold/FSO-range sensitivity beyond what Round 6
    already covers; an actual Zenodo DOI mint (needs the author's own
    account); an arXiv/Zenodo preprint deposit for this paper's companion
    (needs the author's own account); a full per-seed
    `igd_union` for the ablation (would need a full $80$-run ablation
    rerun); and an actual MOEA/D decomposition/neighbourhood-size tuning
    sweep (disclosed only as a Threats-to-validity limitation, matched to
    the same lack of tuning for all five algorithms, not run).

29. **Post-write adversarial self-audit found and fixed three genuine
    transcription/staleness defects the numbers-pass above had missed**:
    (i) a one-digit rounding slip in the Conclusion ("$0.236$" instead of
    the correct "$0.235$" for NSGA-II's mean HV); (ii) an entire
    duplicate paragraph left over from before the 5-algorithm/30-seed
    expansion, still quoting the old 3-algorithm/10-seed HV/IGD/
    unique-solution-count numbers and "nine tests" directly beneath the
    already-updated version of the same discussion -- removed rather
    than reconciled, since it was fully superseded; (iii) the
    "Why MOEA/D underperforms" paragraph and the "Why NSGA-II leads
    NSGA-III" paragraph still quoted the old $9.8\pm2.3$/$97.7$/$72.3$
    unique-solution counts and $9.9\%$/$41\%$ HV gaps -- corrected to
    the actual Table~\ref{tab:algo-comparison} values ($9.5\pm3.0$/
    $96.9$/$72.2$ and $10.4\%$/$54.3\%$). Also found: `highlights_upload.tex`
    (and its compiled PDF, shipped inside the review package) had not
    been updated in step with the in-manuscript highlights block and
    still advertised "NSGA-II, NSGA-III and MOEA/D... 10 seeds" and a
    "3-algorithm Pareto front" -- resynced to match `main.tex` exactly.
    A bug in `code/verify_package.sh` itself was also found and fixed:
    `unzip -l "$ZIP" | grep -q "$f"` under `set -o pipefail` can report a
    false failure (SIGPIPE from `grep -q`'s early exit counts as a
    pipeline failure even when the match was found), which produced a
    spurious "CHANGELOG.md MISSING from zip" on an otherwise-correct
    package; fixed by capturing the listing once into a variable before
    grepping it.

**Round 6-followup** (two independent external reviews of the Round 6
package: one broader review reiterating items already scoped out with
disclosure, and one that verified 108 values against the released
artefacts with zero discrepancies before flagging two MAJOR and seven
MINOR remaining issues):

30. **Root-caused and fixed a real oracle bug**: an edge whose assigned
    medium is physically infeasible at that city pair's distance (e.g.\
    FSO beyond `FSO_MAX_DIST=150`km) has `LINK_FID=0` and `LINK_LAT`
    set to the corresponding `*_link()` function's `1e6`-second
    infeasibility sentinel, not a real latency -- and
    `_build_best_medium_adjacency` did not exclude such edges from the
    path-search graph. The parametric decoder's own search bounds never
    propose such a pair (verified: re-evaluating 95 cached Pareto-optimal
    decoder solutions through the fixed oracle gives byte-identical
    objective vectors, max diff `0.0`), but the unconstrained categorical
    baseline can and did: NSGA-III categorical seed=48 reported a mean
    latency of `3696.99`s (vs.\ `<1.4`s everywhere else in the project),
    traced to exactly one such edge dominating the average once diluted
    over ~270 reachable pairs (`1e6/270≈3703`). Fixed by excluding
    zero-fidelity edges from the adjacency graph entirely, matching the
    physical model's own semantics; the categorical baseline (the only
    experiment this affects) was rerun from scratch. Post-fix, the
    max.\ latency anywhere in the categorical baseline's fronts is
    `0.373`s (was `3696.99`s). The outlier had also been silently
    distorting the shared decoder-vs-categorical reference front's
    normalisation bounds for *every* row of
    `binary_vs_decoder_shared_ref_v3.json`, not only the categorical
    baseline's own: decoder HV moved from
    `0.522`/`0.506` to `0.477`/`0.461` (NSGA-II/NSGA-III) and
    categorical HV from `0.0102`/`0.0110` to `0.0067`/`0.0071` --
    the qualitative conclusion (decoder vastly outperforms the
    categorical baseline) is unchanged and the gap is if anything
    larger under the corrected normalisation (factor of
    `~65`-`70`, was `~50`); Table~\ref{tab:categorical-baseline} and
    its surrounding prose in `main.tex` were updated accordingly.
31. **Ablation raised from 10 to 30 seeds** (external review, MAJOR
    finding): every ablation $p$-value that was significant sat exactly
    at the $n=10$ Wilcoxon resolution floor ($p_{\mathrm{raw}}=2/2^{10}=
    0.001953125$, Cliff's $\delta$ saturated at $\pm1.00$ for all 5
    affected variants) -- the identical pathology the manuscript's own
    Table~\ref{tab:wilcoxon} caption already explains and resolves for
    the $n=30$ main comparison, left uncorrected here. Reran the full
    8-variant $\times$ 30-seed ablation (240 fresh NSGA-III runs,
    parallelised); the scientific conclusion is unchanged (fibre and
    satellite both required; C and E remain statistically
    indistinguishable from the full decoder) but every previously-floored
    $p$-value now sits at $p_{\mathrm{Holm}}<10^{-7}$, genuinely
    well-powered rather than resolution-limited. The symmetric
    `igd_union` metric (added Round 6) is now computed per-seed with a
    real mean$\pm$std over the same 30 seeds, superseding the earlier
    pooled single-value approximation.
32. **Abstract exceeded the journal's actual 250-word limit** (262
    words; confirmed directly against the journal's guide for authors, not
    assumed) -- trimmed to 222 words.
33. **`QI-Bench-Pareto-v1/figures/` was silently stale**: since Round 6,
    `generate_figures_v3.py` only ever wrote into `paper2/` (where
    `main.tex` compiles from), never into the artefact's own
    `figures/` directory, so the shipped package's figures (dated
    pre-Round-6, 3-algorithm plots) diverged from the actually-compiled
    manuscript's figures (5-algorithm plots) without either script
    reporting an error. Fixed by writing both locations from one
    generation step.
34. **Wilcoxon pairwise table (Table~\ref{tab:wilcoxon}) was missing
    the `\resizebox` wrapper** every other wide table in this manuscript
    uses, overflowing the page margin by ~16pt once it grew from 3 to
    10 comparison rows. Fixed.
35. **Added a Discussion paragraph connecting the effective-dimensionality
    finding (Threats to validity item viii, PC1$=91.3\%$) to the main
    algorithm-comparison result** (external review, MAJOR finding): a
    dominance-based algorithm (NSGA-II) winning over four algorithms
    whose defining mechanisms specifically target many-objective
    dominance resistance is not a coincidence given the found front's
    effective dimensionality is close to one or two, not four -- this
    was previously left as two disconnected observations in separate
    sections with no bridging sentence, which a many-objective-literate
    reviewer would very likely connect unfavourably on their own.
36. **Fixed a genuinely misleading (if unused) field name**:
    `binary_baseline_v3_results.json`'s `metrics_vs_shared_reference`
    key sounds like the decoder-vs-categorical shared-reference
    comparison in `binary_vs_decoder_shared_ref_v3.json` (the numbers
    the manuscript actually cites), but is really a reference front
    shared only between the two categorical conditions themselves --
    never cited in the manuscript, but the misleading name could read
    as an unresolved contradiction between the two files. Renamed to
    `metrics_vs_categorical_only_reference`.
37. **`.zenodo.json`/`CITATION.cff` incorrectly declared `license: MIT`
    for the whole deposit**, when the package is dual-licensed
    (code/ MIT, data/results CC-BY-4.0, per `LICENSE-CODE`/
    `LICENSE-DATA` and the README's own License section) -- a real
    metadata inconsistency that would have misrepresented the data's
    license on the public Zenodo record. Fixed: `CITATION.cff` now uses
    CFF's native multi-license array (`[MIT, CC-BY-4.0]`) plus an
    explanatory `license-note`; `.zenodo.json` sets the single
    deposit-level SPDX field to `CC-BY-4.0` (the majority of archived
    content) with an explanatory `license_note` pointing to the
    per-component split, since Zenodo's legacy schema only accepts one
    top-level license string.
38. **`generate_tables_v3.py` had a hardcoded "10 seeds" ablation-table
    caption/comment** that would have silently gone stale again the
    next time the seed count changed; made dynamic
    (`len(abl['config']['seeds'])`), matching the pattern already used
    for the algorithm-comparison table.
39. Checked and NOT reproduced: a second review's claim of "52pp/43
    files/5 figures" announced vs.\ actual package contents -- no such
    numbers appear anywhere in this project's docs; not acted on since
    it could not be verified against any real artefact.
40. **Still explicitly deferred, not silently dropped** (unchanged in
    scope from item 28): multi-instance validation beyond the single
    100-city benchmark, convergence curves, a dedicated seed-offset
    ($\sigma$) ablation, a hub-aware decoder variant, and updating the
    bibliography's recency (currently 10/70 entries from 2023 or later)
    -- all remain large-scope items requiring either substantial new
    compute campaigns or a literature review, and were not attempted
    without the author's explicit go-ahead given their cost.

**Round 6-followup, second pass** (a second independent pre-submission audit
of the Round 6-followup package, which verified 108 values with zero
discrepancies before flagging two real visual/documentation defects, four
statistical transcription errors, and several minor metadata/wording
issues):

41. **Figure~\ref{fig:pareto-pairwise}'s legend was visibly cut down to 3
    of 6 categories** (NSGA-III, SMS-EMOA, Shared only -- NSGA-II, MOEA/D
    and RVEA were missing from the rendered PDF). Root cause: `ncol=3`
    for 6 legend entries wraps to 2 rows, and the reserved top margin
    only had room for one row before the figure was saved without
    `bbox_inches='tight'`. Fixed: single-row legend (`ncol=len(handles)`)
    and `bbox_inches='tight'` on save.
42. **Figure~\ref{fig:ablation}'s title was hardcoded to "10 seeds"**,
    directly contradicting its own caption's "30 seeds" two lines below
    on the same page -- missed when `generate_tables_v3.py`'s ablation
    caption was made dynamic (item 38) because the *figure* title lives
    in a separate script, `generate_figures_v3.py`, not caught by that
    fix. Made dynamic; also switched this figure from a bar+errorbar
    chart to a boxplot with jittered per-seed points, matching Figure~1's
    Round-6 fix, now that $n=30$ makes the per-seed distribution worth
    showing directly.
43. **A $p$-value was printed as literally `0.000`** in the Wilcoxon
    table for every comparison below the 3-decimal rounding threshold.
    Fixed (`generate_tables_v3.py` and the corresponding manual table in
    `main.tex`) to print `<0.001` instead whenever $p_{\mathrm{Holm}}<0.0005$.
44. **A summary-prose statistical claim was false**: "all HV and IGD
    pairwise comparisons remain significant at $p_{\mathrm{Holm}}\leq0.001$"
    -- the actual maximum is $p_{\mathrm{Holm}}=0.0045$ (NSGA-III vs.\
    SMS-EMOA, IGD), which the table itself already correctly showed as
    $0.005$, contradicting the summary sentence right above it. Fixed to
    state the correct bound ($\leq0.005$) and name the specific weakest
    comparison.
45. **The spacing-significance comparison counts were wrong**: prose
    claimed "four pairs involving SMS-EMOA or RVEA" are significant and
    "six comparisons among NSGA-II, NSGA-III and MOEA/D" are not --
    recounting directly from Table~\ref{tab:wilcoxon}'s own spacing
    column gives **seven** significant pairs involving SMS-EMOA or RVEA
    and exactly **three** non-significant comparisons among
    NSGA-II/NSGA-III/MOEA-D ($7+3=10$, matching the total). Fixed.
46. **An imprecise sentence about the ablation's $p$-value floor**:
    "what changes ... is the $p$-value's distance from its own floor" is
    not quite right -- both designs' $p$-values sit exactly at their
    floor (distance $=0$ in both cases); what changes is the *level* of
    that floor. Reworded.
47. **Added a clarifying note for variant I (FSO-only)'s HV std of
    `0.0000`**, confirmed not a rounding artefact (30 per-seed values
    genuinely span only `0.099868`--`0.099890`, consistent with only 11
    of 4,950 city pairs being FSO-reachable at all) rather than an
    apparent transcription error.
48. **Removed the duplicated Highlights page from the main manuscript
    PDF**: `\begin{highlights}...\end{highlights}` rendered as a full
    standalone page~1, redundant with the separately-submitted
    `highlights_upload.tex`/`.pdf` that the submission system already requires. Removed
    from `main.tex`, kept only in the dedicated upload file.
49. **Replaced "NSGA-III" as a standalone keyword** with "Evolutionary
    multi-objective optimisation": with five algorithms now compared and
    NSGA-II (not NSGA-III) achieving the best mean HV, singling out one
    non-winning algorithm as a keyword read as arbitrary.
50. **Removed remaining internal revision-history narration from the
    manuscript body** ("added this round", "as of this round"), per the
    same house rule already applied throughout earlier rounds -- these
    two instances had been introduced by this round's own edits and were
    missed at the time.
51. **"External review"/"external reviewer" terminology in
    `README.md`/`REVIEW_BRIEF.md` reworded** to "independent pre-submission
    technical audit(s)" with an explicit note that these are AI-assisted,
    author-reviewed audits, not formal journal peer review -- the
    original phrasing could read as claiming peer review had already
    occurred.
52. **`README.md`/`REVIEW_BRIEF.md` said "submitted to *Swarm and
    Evolutionary Computation*"** before the paper has actually been
    submitted; changed to "prepared for submission to".
53. **Checked and NOT reproduced**: a claim that `CHECKSUMS.sha256`
    lists a `zenodo.json` entry without its leading dot -- the current
    file correctly lists `.zenodo.json` and validates; not acted on
    since it could not be reproduced against the actual artefact.
54. **Still not fixed, disclosed rather than silently dropped**:
    `CITATION.cff`'s `repository-code` field remains
    `REPLACE_WITH_REAL_REPOSITORY_URL` -- this genuinely cannot be filled
    in without the author's real repository URL, unlike the other items
    above.

55. **CORRECTION to item 14 above (Round 5): the corresponding-author
    email "unification" in Round 5 was itself wrong.** Round 5 changed
    `\ead{}` from `jesus.gil@universidadeuropea.es` to
    `jesus.gil@universidadeuropea.es`, describing the shorter
    address as merely an Editorial-Manager login address. This was
    backwards: `jesus.gil@universidadeuropea.es` does not exist and
    was a fabricated variant; the one real institutional address is
    `jesus.gil@universidadeuropea.es`, which is also the address used to
    log in to journal portals. Reverted `main.tex`'s `\ead{}` back to the
    real address. Caught via this project's own author-identity memory
    record, not by either external audit this round.

56. **Categorical baseline raised from 10 to 30 seeds (real rerun, not a
    reworded caveat).** The manuscript claimed the categorical baseline's
    10 seeds "match the statistical power of the main comparison" -- true
    when written, but false once the main comparison was raised to 30
    seeds in Round 6, and left uncorrected since. Rather than just
    rewording the claim, reran the categorical baseline for real: added
    `--algos`/`--seeds`/`--aggregate-only` sharding and per-(condition,seed)
    checkpointing to `binary_baseline_v3.py` (same pattern as the other two
    experiment scripts), extended `SEEDS` to `range(42, 72)`, and ran all
    60 (2 conditions x 30 seeds) NSGA-II/NSGA-III categorical runs fresh
    with the item-30 oracle fix already in place (confirmed: max latency
    across the new categorical fronts is `1.986`s, no recurrence of the
    `3696.99`s artefact). Also fixed `build_binary_vs_decoder_shared_ref_v3.py`
    to use the same 30 seeds for the decoder side too (it had been using
    only the first 10 of the decoder's 30 available seeds for this specific
    comparison table) and to use `ddof=1` (was `ddof=0`, missed in the
    original Round 6 pass). Updated numbers in `main.tex`'s
    `Table~\ref{tab:categorical-baseline}` and surrounding prose: decoder
    HV `0.406`/`0.385` (NSGA-II/NSGA-III), categorical HV `0.0092`/`0.0087`
    -- a `~44x` gap (was reported as `~65`-`70x` under the old, only-10-seed
    comparison); qualitative conclusion unchanged.
57. **Documentation/code-comment staleness swept a second time,
    caught by a fresh independent audit's direct question about
    completeness**: `run_ablation_v3.py`'s own module docstring still said
    "10 seeds" (missed when the SEEDS constant itself was changed);
    `build_global_pareto.py`'s docstring still said "all 3 algorithms";
    `stats_analysis_v3.py` had two comments referencing "the same 10
    seeds (42-51)" for the ablation. All fixed. Also found
    `recompute_ablation_stats_v3.py` (Round 6's pooled-IGD_union
    workaround) is now actively dangerous to run: it expects the old
    10-seed shape and would silently overwrite the current, better
    per-seed `igd_union` with a coarser pooled approximation if invoked
    against the current `ablation_v3_results.json`. Marked SUPERSEDED in
    its own docstring and removed from `README.md`'s active reproduction
    steps (kept in `code/` for historical reference only). `README.md`'s
    Contents tree and reproduction-time estimates updated throughout for
    the current 30-seed ablation and categorical baseline. Also verified
    no file in this package other than `main.tex` (`CITATION.cff`,
    `.zenodo.json`, `README.md`, `REVIEW_BRIEF.md`) contains an email
    address at all, so item 55's fabricated-email fix had no further
    propagation to clean up.

59. **Renamed "global Pareto front" to "pooled evolutionary
    (non-dominated) front" throughout the manuscript** (16 occurrences).
    The concept was already correctly defined at first use as the union
    of five algorithms' pooled, deduplicated, mutually-non-dominated
    points -- not a claim of true problem-wide Pareto-optimality -- but
    the recurring short label "global Pareto front" in captions/table
    rows read as overclaiming when encountered in isolation, especially
    given the manuscript's own finding that three classical hub-and-spoke
    constructors are non-dominated by this front (i.e.\ genuinely outside
    it). No numbers changed, wording only.
60. **Softened "these four objectives ... define a genuine many-objective
    optimisation problem"** (Introduction) to note the four objectives
    are jointly optimised throughout while forward-referencing the
    Discussion's effective-dimensionality finding, instead of asserting
    many-objective-ness as settled before that finding is presented.
61. **Clarified the Related Work citation of spacing~\cite{schott1995spacing}**:
    added a parenthetical noting this paper uses a Euclidean-distance
    variant (cross-referenced to its precise definition in
    Section~\ref{sec:experiments}), not Schott's original Manhattan-distance
    formula, so the citation reads as attributing the metric's origin
    rather than claiming exact formula compliance.
62. **Added a normalisation-mismatch caveat to
    Table~\ref{tab:categorical-baseline}'s caption**, since NSGA-II's
    decoder HV appears as $0.406$ there but $0.235$ in
    Table~\ref{tab:algo-comparison} for the same underlying runs (different
    pooled reference fronts/normalisations per block) -- flagged by
    external review as a plausible source of reader confusion despite
    both numbers being individually correct.
63. **Not fixed, requires the author's action**: `CITATION.cff`'s
    `repository-code` field is still `REPLACE_WITH_REAL_REPOSITORY_URL`
    (flagged repeatedly across rounds; genuinely cannot be filled in
    without the real repository URL).
64. **Explicitly out of scope for this round, raised again by external
    review but not attempted without author confirmation given cost**:
    inverting the statistical hierarchy (independent-samples as primary,
    paired as sensitivity) -- kept paired-by-seed as primary since it is
    the methodologically standard design for this kind of repeated-seed
    EC comparison, per the same Derrac/García/Herrera-style literature
    the reviews themselves cite; recomputing spacing with Schott's exact
    Manhattan-distance/ddof=1 formula across all 450 stored runs;
    multi-instance validation; convergence curves; a sigma ablation; a
    dedicated categorical-operator baseline; external oracle validation;
    topology maps; macro-generated (non-hardcoded) manuscript numbers;
    safe (non-pickle) serialisation; and a test suite/Makefile. All are
    real, legitimate suggestions the author should weigh explicitly
    against the added compute/writing cost before the next round.

**Round 7** (author opted to complete the full remaining list from item 64
above rather than defer it):

65. **Spacing recomputed with Schott's (1995) canonical formula**
    (Manhattan/L1 nearest-neighbour distance, sample std ddof=1) across
    all 450 already-stored runs (150 main comparison + 240 ablation + 60
    categorical baseline) -- no reruns needed, only cached checkpoints
    reprocessed. Previous implementation used Euclidean (L2) distance and
    ddof=0, an undisclosed-as-different variant from the metric cited by
    name (Schott 1995) in Related Work. HV/IGD unaffected (verified:
    unchanged to full precision); spacing values changed materially for
    some algorithms (e.g.\ RVEA $0.095\pm0.036\to0.171\pm0.053$); the
    qualitative significance pattern (7 pairs involving SMS-EMOA/RVEA
    significant, 3 among NSGA-II/NSGA-III/MOEA-D not) is unchanged, only
    the exact $p$-values/effect sizes shifted. Updated
    Table~\ref{tab:algo-comparison}, Table~\ref{tab:friedman},
    Table~\ref{tab:wilcoxon}, the ablation table, and all surrounding
    prose. Also discovered and fixed: re-running the joint Holm correction
    after the spacing recompute shifted the *IGD* column's weakest
    $p_{\mathrm{Holm}}$ from $0.0045$ to $0.0057$ (rounds to $0.006$, not
    $0.005$) -- an expected consequence of Holm's step-down procedure
    depending on the whole joint family of 30 raw $p$-values, not a new
    bug, but the manuscript text needed updating to match.
66. **Inverted the algorithm-comparison statistical hierarchy**:
    Kruskal-Wallis omnibus + Holm-corrected pairwise Mann-Whitney
    (independent samples, with Cliff's $\delta$/Vargha-Delaney $A_{12}$
    effect sizes) is now the PRIMARY analysis throughout (abstract,
    Introduction contributions list, Results, Threats to validity item
    v, Conclusion); the paired Friedman/Wilcoxon-by-seed-index design is
    now explicitly framed as a "common-seed sensitivity check", not the
    primary claim, since it additionally assumes seed $s$ is a
    comparable initial condition across algorithms with structurally
    different population dynamics (an assumption the independent-samples
    design does not need). Added `cliffs_delta`/`vargha_delaney_a12` to
    the independent-samples Mann-Whitney output in `stats_analysis_v3.py`
    (previously only $p$-values were reported there, no effect size).
    Both designs agree on the significance verdict for all 30
    comparisons; this agreement is now framed as evidence for the primary
    (independent-samples) result rather than the paired result absorbing
    the confirmation. Added new Tables~\ref{tab:kruskal}/\ref{tab:mannwhitney}
    (primary) alongside the retained, re-captioned
    Tables~\ref{tab:friedman}/\ref{tab:wilcoxon} (sensitivity check). Added
    a real citation (`derrac2011practical`, verified via web search: Derrac,
    García, Molina \& Herrera, \emph{Swarm and Evolutionary Computation}
    1(1):3-18, 2011 -- the standard EC-statistics-methodology reference,
    fittingly published in this paper's own target journal) supporting
    the paired-design's own methodological legitimacy as a sensitivity
    check.
67. **Ran the sigma ablation** (three NSGA-III variants, matched budget,
    $15$ seeds: A = standard $11$-variable decoder with $\sigma$
    optimised; B = $10$-variable, $\sigma$ fixed to a constant; C =
    $10$-variable, fitness averaged over $K=3$ realisations). New script
    `run_sigma_ablation_v3.py`. Real result: A vs.\ B statistically
    indistinguishable on HV (Mann-Whitney, Holm-corrected,
    $p_{\mathrm{Holm}}=0.115$) -- optimising $\sigma$ gives no
    detectable advantage over fixing it -- while C is significantly
    worse than both ($p_{\mathrm{Holm}}\leq0.006$), consistent with
    averaged fitness diluting rather than denoising the search signal.
    Variant A also shows $2$--$5\times$ larger HV std than B/C, the
    same qualitative signature as NSGA-II's main-comparison bimodality.
    Added new Section~\ref{sec:results-sigma} and rewrote Threats to
    validity item (vii) to report this finding instead of listing it as
    future work.
68. **Added topology maps** (new Figure~\ref{fig:topology-maps},
    `generate_topology_maps_v3.py`, using `cartopy` -- installed
    successfully, no external account needed): 4 real, selected
    topologies rendered on a world map, coloured by medium
    (fibre/satellite/FSO) -- max-fidelity, max-coverage, knee-point
    (closest to the ideal corner in normalised objective space), and
    the non-dominated $10$-hub hub-and-spoke classical constructor
    (rebuilt deterministically from `classical_baselines_v3.py`'s exact
    construction). Addresses the recurring complaint that a paper titled
    "Global Quantum Internet Topologies" never showed an actual
    designed network. Selection indices/objective values saved to
    `topology_maps_v3_selection.json` for reproducibility.
69. **Safe serialisation**: `unique_topologies_v3.npz` was the only
    public results file requiring `allow_pickle=True`
    (`eu`/`ev`/`et` were `dtype=object` arrays, one variable-length
    array per topology). Replaced with a flat, CSR-style
    (`edge_u`, `edge_v`, `edge_type`, `topology_offsets`) representation
    -- topology $i$'s edges are
    `edge_u/edge_v/edge_type[topology_offsets[i]:topology_offsets[i+1]]`
    -- all plain numeric dtypes, loads with `allow_pickle=False`.
    Verified byte-for-byte equivalent content after conversion (same
    1,347 topologies, same max-fidelity/coverage/knee selections when
    `generate_topology_maps_v3.py` was updated to the new format and
    rerun). Checked `global_pareto_v3.npz` and `ablation_v3_pareto.npz`
    (the package's other two main `.npz` artefacts): both already load
    without `allow_pickle`, no changes needed. The per-run `.pkl`
    checkpoints (`step1`-`step6` directories) remain pickle-based, since
    they are optional, trusted-local reproducibility artefacts, not
    required to regenerate any manuscript figure or table.
70. **Added a real test suite** (`tests/test_reproducibility.py`, $14$
    tests, run with `python -m pytest tests/ -v`) covering: canonical
    (Manhattan/ddof=1) spacing correctness with a numeric example that
    would fail under the old Euclidean/ddof=0 implementation; seed
    counts and evaluation budgets for the main comparison and ablation;
    no-pathological-latency (the confirmed `3696.99`s bug must not
    recur); no-NaN/inf in any stored metric; abstract $\leq250$ words;
    highlight bullets $\leq85$ characters; no fabricated email variant;
    no visible revision-history narration in `main.tex`; placeholders
    exist in exactly the two disclosed locations
    (`CITATION.cff`, `.zenodo.json`) and nowhere else; the two safely-
    serialised `.npz` files load without `allow_pickle`; and
    `sha256sum -c CHECKSUMS.sha256` passes. Added a `Makefile` with
    `test`/`figures`/`tables`/`paper`/`checksums`/`verify`/`package`
    targets.
71. **Added a second categorical baseline with category-native
    operators** (`categorical_operators_baseline_v3.py`: genuinely
    discrete representation, no real-valued relaxation; uniform
    categorical crossover; add-edge/remove-edge/change-medium mutation),
    tested at 4 initial densities ($3\%$, $10\%$, the decoder's own
    verified median density $10.9\%$, $15\%$; NSGA-II, $10$ seeds,
    matched budget). Real result: density matters enormously within
    this baseline (own-reference HV $0.011\to1.195$ from $3\%$ to
    $15\%$ density, not saturated); against a reference shared with the
    decoder (same comparison as the original SBX/PM baseline), the
    best density reaches HV$=0.0117\pm0.0016$ vs.\ the decoder's
    $0.257\pm0.040$ -- a $\sim22\times$ gap, narrower than the original
    baseline's $\sim44\times$ but still decisive. Added new paragraph
    in Section~\ref{sec:results-categorical}; softened the original
    baseline's caveat sentence to point at this real result instead of
    speculating about untested operators.
72. **Added convergence curves** (new Figure~\ref{fig:convergence},
    `run_convergence_v3.py` + `generate_convergence_figure_v3.py`):
    per-generation HV/IGD tracked for all 5 algorithms ($10$ seeds,
    matched budget, same normalisation/IGD reference as the main
    comparison, no rerun of the main 30-seed results needed). Real
    result: all five algorithms reach $\geq97.8\%$ of their own final
    median HV by half the evaluation budget, so $80$/$99$ generations
    is generous, not tight; MOEA/D's IGD degrades sharply within the
    first few hundred evaluations and never recovers, consistent with
    its early reconvergence onto a small duplicate-heavy cluster.
    Removed "convergence curves" from the Conclusion's future-work list
    (now done) and added a note there explaining SeQUeNCe's status
    (installed successfully, no account needed; a full protocol-level
    validation deferred as explained in item 74 below).
73. **Added multi-instance validation** (new
    Section~\ref{sec:results-multi-instance},
    `multi_instance_oracle_v3.py` + `run_multi_instance_v3.py`): 3
    additional REAL-city instances (verified subsets of the same
    100-city list, no fabricated coordinates) -- `global_n50` (50
    cities, every other of the 100), `regional_eu_na` (50 cities,
    Europe+North America only), `metro_europe` (30 cities, Europe
    only) -- all 5 algorithms, 10 seeds each, matched budget. The
    parameterised oracle reimplementation was verified to reproduce
    the original 100-city oracle's output EXACTLY (`np.allclose`,
    identical fidelity/rate/latency/coverage) before any real run.
    Real, honestly-reported result: the qualitative ranking (NSGA-II/
    NSGA-III top two, SMS-EMOA third, RVEA/MOEA-D last) holds on all 3
    instances; NSGA-II's specific edge over NSGA-III is significant on
    `regional_eu_na` and `metro_europe` (Mann-Whitney $p<10^{-3}$,
    matching the main benchmark) but NOT significant on `global_n50`
    ($p=0.385$) -- reported as a genuine instance-dependence, not
    smoothed over, with a stated-but-untested plausible explanation
    (NSGA-II's advantage is itself driven by a seed-minority effect, so
    10 seeds may simply undersample it on this instance). Did NOT
    attempt a genuinely independent $N=200$ instance: that requires
    sourcing 100 additional real cities, not subsetting the existing
    list, and inventing coordinates to reach $N=200$ was judged worse
    than leaving this explicitly open (Conclusion future-work list
    updated accordingly).
74. **External oracle validation via SeQUeNCe: completed, scoped
    honestly to what is actually comparable.** Installed successfully
    (pip package `sequence` 0.8.5, no account needed -- unlike
    NetSquid, which requires netsquid.org registration). Rather than
    forcing SeQUeNCe to numerically reproduce this paper's exact
    Werner-fidelity BDCZ physics (high risk of a subtly-wrong,
    miscalibrated "match" -- see `sequence_oracle_validation_v3.py`'s
    module docstring for the full reasoning), implemented a
    directional cross-check: SeQUeNCe's `QuantumChannel` computes
    photon loss as `1 - 10**(-attenuation*distance/10)`, structurally
    identical to this paper's per-span fibre transmittance law: `eta =
    10**(-attenuation_db_per_km*dist_km/10)`. Ran SeQUeNCe's own
    midpoint-heralded entanglement-generation protocol
    (`EntanglementGenerationA`, adapted directly from SeQUeNCe's
    official `demo_for_beginners/two_node_eg.ipynb` reference example)
    between two quantum routers at 7 single-span fibre distances
    (1-50 km, within `FIBER_REPEATER_SPACING=50` km so no swapping is
    needed for a fair single-span comparison), 10 seeds each, counting
    how many of 50 quantum memories per node reached successful
    entanglement within a fixed simulation window. Real result (see
    `results/sequence_oracle_validation_v3.json`): median
    entangled-memory count falls monotonically with distance
    (47.5/41.5/32.0/18.5/7.0/3.5/1.0 across 1/5/10/20/30/40/50 km,
    zero rank violations), Spearman rho = 1.0 against the analytical
    oracle's own transmittance decay across the same 7 distances.
    Explicitly scoped: this validates the RATE dimension only.
    SeQUeNCe's baseline elementary-link protocol fixes the produced
    state's fidelity to a constant (`memory.raw_fidelity`), not a
    function of channel loss or distance, so it does not offer a
    comparable fidelity-vs-distance check without custom noise
    modelling beyond the standard protocol -- disclosed as such in
    Threats to validity (item iii) rather than glossed over. Full
    multi-hop entanglement-swapping cross-validation remains future
    work (Conclusion), for the same protocol-fidelity reason as
    before.
75. **Physics audit of the BDCZ oracle's fibre latency formula: found and
    fixed a real, confirmed bug, then fully recomputed every
    campaign in this package.** External review asked for a formal
    physics audit of the oracle's fibre/satellite/FSO latency
    formulas, with citations, before any further work -- see
    `PHYSICS_ORACLE_AUDIT.md` for the full derivation. Finding:
    `fiber_link()`'s latency was `dist_km*1000/C_FIBER*(2*n_spans)`,
    which -- since `span_len = dist_km/n_spans` -- scales
    *quadratically* with distance for a fixed repeater spacing.
    Verified against two independent references: (a) a naive
    sequential-repeater model, whose round-trip propagation time
    collapses to `2*dist_km/v` independent of `n_spans` once total
    distance is fixed; (b) Sangouard, Simon, de Riedmatten, Gisin,
    *Rev. Mod. Phys.* 83, 33 (2011), Eq. (12), whose elementary
    period `L0/c` *decreases* with more nesting levels at fixed total
    distance -- the opposite direction from the buggy code. Also
    matches the manuscript's own stated definition ("latency is the
    sum of per-link propagation times", Introduction), which the
    quadratic code contradicted. Root cause: an erroneous extra
    `n_spans` factor multiplying the already-total distance. **Same
    bug confirmed present, byte-for-byte identical, in the companion
    QI-Bench-100 paper's own canonical oracle** -- NOT touched here,
    a decision for the authors of that paper. Fixed in this package only:
    `latency_s = dist_km*1000/C_FIBER*2`.
    Backed up pre-fix `results/` and `main.tex` to
    `/home/jesus/paper2_v4/PRE_LATENCY_FIX_BACKUP_20260728/` before
    deleting any checkpoint. **Fully recomputed all 735 evolutionary
    runs** across every campaign (main comparison, ablation,
    categorical/binary baseline, sigma ablation, categorical-operator
    baseline, convergence, multi-instance) under the corrected
    oracle -- not a patch, since latency is a directly-optimised
    objective in every one of these. Real, material, honestly-reported
    consequences of the fix, verified against the regenerated JSON
    files before writing any number down:
    - Pooled evolutionary front: 659 unique objective vectors (was
      846), 1,120 unique topologies (was 1,347). Pooled fibre edge
      fraction rose from 5.7% to 11.6% (per-topology mean fibre
      fraction 15.8%→29.5%) -- the corrected, less latency-penalised
      fibre medium is now used substantially more by the evolutionary
      search.
    - **The ablation's headline finding reverses**: "fibre and
      satellite are both required" is no longer supported.
      Fibre-only (B) and fibre+satellite (C) are now statistically
      indistinguishable from the full decoder
      ($p_{\mathrm{Holm}}=0.115$ and $1.000$); only variants lacking
      fibre entirely (G/H/I) still collapse HV significantly.
      Fibre+FSO without satellite (D) is unexpectedly significantly
      *higher*-HV than the full decoder ($p_{\mathrm{Holm}}=0.019$),
      reported directly without a fully established mechanism. New
      corrected finding: fibre is necessary; satellite's contribution
      beyond fibre alone is not statistically detectable at this
      benchmark and oracle configuration.
    - Main comparison: NSGA-II HV rose from 0.235±0.050 to
      0.301±0.039 (same qualitative ranking, NSGA-II still best), but
      the common-seed paired sensitivity check now DISAGREES with the
      primary independent-samples test on NSGA-II vs. NSGA-III itself
      (7 of 30 total comparisons disagree, up from 0) -- reported
      directly rather than smoothed over (Section 5.1); NSGA-III vs.
      SMS-EMOA on HV also lost significance under the primary test.
    - NSGA-II's seed bimodality: high-HV cluster grew from 7/30 to
      10/30 seeds; NSGA-III now shows a comparably-sized 10/30
      high-HV cluster (was 3/30), with 6 of 10 seeds shared between
      the two algorithms' high clusters (was 2 of 3) -- a materially
      stronger basin-sharing signal than previously reported.
    - Multi-instance validation: on `metro_europe` specifically, the
      NSGA-II-vs-NSGA-III point estimate now REVERSES in NSGA-III's
      favour (was NSGA-II ahead but not significant; now NSGA-III
      ahead, also not significant) -- the smallest, most
      geographically concentrated of the three additional instances.
    - Effective dimensionality: PC1 rose from 91.3% to 95.6% of
      variance (PC1+PC2 97.4%→99.1%) -- the front is now *even more*
      effectively low-dimensional, strengthening rather than
      undermining that finding.
    - Categorical baseline gaps changed in scale (SBX/PM baseline gap
      44x→58x; category-native-operator baseline gap 22x→6.4x) but
      not in direction: the parametric decoder still dominates both.
    - Dominance check (11/14 classical constructors dominated, the 3
      hub-and-spoke variants non-dominated) and the HV
      reference-point sensitivity ranking (stable across 1.05/1.10/1.20)
      were qualitatively unchanged, since these depend less directly
      on the specific latency values.
    Also fixed, while auditing this section: the Methodology
    paragraph describing spacing as "nearest-neighbour Euclidean
    distances" was stale prose left over from before the Round 7
    Manhattan/ddof=1 fix -- corrected to match the code and cite
    Schott (1995) directly. `highlights_upload.tex`'s highlight 4
    ("fibre and satellite are both required") and highlight 5
    ("Global 5-algorithm Pareto front", pre-dating the
    pooled-evolutionary-front renaming) were both stale and are now
    corrected. Recompiled (67 pages, 0 undefined refs), all 15 tests
    pass except the expected stale-checksum failure pending
    repackaging.
76. **Adversarial audit found and fixed one transcription slip and one
    real traceability gap in the Round 8 rewrite.** (a) Effective
    dimensionality's PC1 was mistranscribed as 95.7% in `main.tex`
    (three places) and `README.md`/`CHANGELOG.md`/`REVIEW_BRIEF.md`;
    the actual value in `effective_dimensionality_v3.json` is 95.6486%,
    which rounds to 95.6%, not 95.7% -- corrected everywhere. (b) The
    multi-instance section's NSGA-II-vs-NSGA-III pairwise HV/p-value
    numbers per instance (own two-algorithm reference front, distinct
    from the all-five-algorithm summary in
    `multi_instance_v3_results.json`) had been computed ad hoc during
    manuscript writing and never persisted to a committed artefact --
    an external review correctly flagged this as unverifiable. Added
    `multi_instance_pairwise_stats_v3.py`, re-derived the numbers from
    scratch independently of the original ad hoc computation, confirmed
    an exact match to what was already written in `main.tex`
    (global\_n50 0.143/0.112/p=0.0013, regional\_eu\_na
    0.787/0.778/p=0.017, metro\_europe 0.191/0.204/p=0.307), and saved
    to `multi_instance_pairwise_stats_v3.json` with a new regression
    test (`test_multi_instance_pairwise_numbers_are_traceable`, 16
    tests total). Also fixed, from the same review round: the oracle
    description (Section 3, "Analytical oracle") claimed Paper 2's
    oracle is "reused here verbatim" / "identical, single codebase"
    with the companion QI-Bench-100 benchmark -- no longer true since
    only Paper 2's copy received the Round 8 latency fix; reworded to
    state the fidelity/rate models are shared but the latency formula
    was corrected here and not (yet) in the companion benchmark's own
    results. Softened the abstract's "the ranking holds on three
    additional real-city instances" to acknowledge the one
    instance-dependent reversal already reported honestly in the body
    text (Section 6.6), instead of only there.
77. **Large author-approved batch addressing the remainder of the
    external review list in full** (three independent reviews'
    combined remaining items). Real, verified, itemised below --
    nothing claimed without a corresponding results file or manuscript
    change:
    - **Typography**: elsarticle's `\paragraph{}` appends its own
      period, so every `\paragraph{Heading.}` in this manuscript
      rendered as "Heading.." in the actual compiled PDF -- confirmed
      by extracting text from the PDF itself, not just inspecting the
      source. Fixed all 24 occurrences.
    - **Stale numbers found and fixed**: Table 3's RVEA cardinality
      was correct in the table (8.1+/-1.9) but a later paragraph still
      quoted the pre-Round-8 value (6.9+/-1.0); the sigma-ablation
      table caption's "Variant I" HV range was still the pre-latency-fix
      value (0.099868-0.099890, should be 0.094166-0.094188);
      the Methodology's spacing description still said "Euclidean
      distances" months after the Manhattan/ddof=1 fix (Round 7) --
      corrected to match the code and cite Schott (1995) directly.
    - **Ablation variant lettering**: labels A-E, G-I (no F) were
      real but undocumented; added a one-line caption clarification
      rather than a risky mass-rename touching checkpoint filenames,
      JSON keys and this changelog's own history.
    - **"Fixed penalty" defined explicitly** in Algorithm 1's
      pseudocode comment (was previously just "a fixed penalty" with
      no value given): $(\bar F,\bar R_{\log},\bar L,C_\tau)=(0.5,0,1.0,0)$.
    - **Path-search naming**: "widest-path Dijkstra" (both in the
      manuscript and the code's own docstring) renamed to "maximum
      Werner-composed-fidelity path search" -- Werner composition is
      neither a min (classical widest/bottleneck path) nor a simple
      product, so "widest-path" was a real misnomer, not just
      informal phrasing.
    - **Satellite latency's previously-undocumented second term**
      (`dist_km*1000/C_FIBER`, flagged in `PHYSICS_ORACLE_AUDIT.md`
      item 2b as needing justification or removal) resolved by adding
      an explicit code comment stating the classical-confirmation-channel
      rationale. No formula or numerical behaviour changed, so no
      re-run was needed; `PHYSICS_ORACLE_AUDIT.md` updated to RESOLVED.
    - **Coverage threshold (0.70) sensitivity, real and computed**:
      re-decoded all 659 pooled-front points and recomputed $C_\tau$
      at 4 alternative thresholds using the exact canonical path
      search (`coverage_threshold_sensitivity_v3.py`). Ranking is
      essentially unchanged at 0.50/0.60/0.80 (Spearman
      $\rho\geq0.99$); at 0.90 it drops to $\rho=0.39$, but this is a
      floor effect (mean $C_\tau\to0.006$), not ranking instability --
      both findings reported, not just the reassuring one.
    - **Routing-policy sensitivity for the effective-dimensionality
      claim, real and computed**: re-decoded all 659 pooled-front
      points and recomputed all four objectives under rate-optimal
      and latency-optimal path search (not just the canonical
      fidelity-optimal one), reusing the same adjacency construction
      (`routing_policy_sensitivity_v3.py`). PC1 stays ~95% under all
      three policies (95.6%/95.0%/94.7%) -- the low effective
      dimensionality is not an artefact of always routing on the
      fidelity-optimal path.
    - **Best-so-far convergence, genuinely re-run**: the previous
      convergence tracking read `algorithm.opt`/`pop` each generation
      (the CURRENT population's own front), which for a finite
      population can regress; reworked `ConvergenceCallback` to
      maintain an explicit cumulative non-dominated archive across
      generations, which guarantees HV is non-decreasing by
      construction (verified per-seed: exactly monotone for all five
      algorithms). Re-ran all 50 convergence runs under the new
      tracking. Real consequence: MOEA/D's HV-at-half-budget rises
      from 88.9% to 96.4%, and its previously "rises sharply then
      plateaus at an elevated level" IGD curve now decreases smoothly
      and monotonically throughout the run (0.237->0.164->0.142) --
      the earlier finding was a population-tracking artefact specific
      to MOEA/D's population repeatedly reconverging, not a property
      of MOEA/D's actual best-found archive. All five algorithms now
      reach >=96.4% of final HV by half budget.
    - **Sigma ablation, genuinely re-run as 10 variables with 3
      independently-tested fixed values**: the "10-variable" B/C
      variants previously kept `n_var=11` with the 11th dimension's
      bounds collapsed to a single point (behaviourally equivalent for
      SBX/PM, but not what "10-variable" should mean at the
      search-space-definition level); reworked so the pymoo Problem
      itself has `n_var=10`, with sigma padded on only immediately
      before the oracle call. Also replaced the single arbitrary fixed
      value (42) with three independently-tested values (0, 42, 500).
      Re-ran all 75 runs (5 variants x 15 seeds). Real, more nuanced
      finding than before: optimising sigma (A) is statistically
      indistinguishable from any SINGLE fixed value tested
      ($p_{\mathrm{Holm}}\geq0.499$ in all three comparisons), but the
      fixed values differ SIGNIFICANTLY from each other ($\sigma{=}0$
      vs.\ $\sigma{=}500$: $p_{\mathrm{Holm}}=0.0047$,
      $\delta=-0.73$) -- so "just fix sigma to anything" is not
      supported; optimising it sidesteps having to guess a good
      constant.
    - **Categorical-representation comparison unified into one master
      table**: extended the category-native-operator baseline's best
      density (d15) from 10 to 30 seeds, then built
      `representation_comparison_master_v3.py`, pooling decoder +
      relaxed-categorical SBX/PM + category-native (d15) into ONE
      shared, deduplicated reference front/normalisation (with
      recorded reference-front and normalisation-bounds hashes).
      Real result: decoder $\mathrm{HV}=0.423\pm0.034$,
      category-native $0.120\pm0.006$ ($\sim3.5\times$ gap),
      relaxed-categorical $0.006\pm0.001$ ($\sim66\times$ gap); all
      pairwise Mann-Whitney comparisons significant
      ($p<10^{-10}$). Same qualitative conclusion as the existing
      block-specific numbers, now on one unambiguous scale.
    - **Random-search sanity-check baseline, real and computed**: 30
      seeds of i.i.d.\ uniform sampling over the same bounds/budget/oracle
      as the main comparison, re-scored against a reference shared
      with NSGA-II's own runs. Real result: NSGA-II
      $\mathrm{HV}=0.256\pm0.041$ vs.\ random search
      $0.212\pm0.013$ -- a real, significant advantage
      (Mann-Whitney $p=1.2\times10^{-8}$) but a modest one in relative
      terms ($\sim21\%$), reported honestly rather than only
      emphasising significance.
    - **Two modern algorithms added under identical conditions**:
      SPEA2 (2001, included for completeness of the classical
      Pareto-based family) and AGE-MOEA-II (Panichella, 2022, a
      genuinely 2020s many-objective method using adaptive geometry
      estimation rather than a fixed reference-direction set), both
      30 seeds, re-scored against a reference shared with the main
      comparison's own pooled front so results are directly comparable
      to Table 3.
    - **A genuinely independent N=200 instance, real cities only**:
      100 additional real cities sourced from GeoNames'
      `cities15000` export (verified against the official download,
      cross-checked a sample against Wikipedia; full sourcing
      methodology in the commit), appended to `cities.csv` (now 200
      rows) and to `multi_instance_oracle_v3.py`'s city table (with a
      correctness check confirming the first 100 rows/distances/link
      tables are byte-identical to the original 100-city oracle
      before trusting any result). This was previously left out on
      the grounds that reaching N=200 without fabricating coordinates
      was not possible with only the original 100-city list; the
      resolution is genuinely sourcing 100 more real cities, not
      relaxing that standard.
    - **Makefile fixed**: `make figures` did not call
      `generate_convergence_figure_v3.py` -- added. (`make checksums`
      was independently re-verified to work correctly already;
      an external review's claim that `make_checksums.py` does not
      exist was checked and found to be incorrect -- it exists at
      `../code/make_checksums.py` relative to this package's root,
      exactly where the Makefile's `checksums:` target already looks.)
    - **`metrics.py` extracted**: `compute_spacing` was duplicated
      identically across 4 scripts; consolidated into a single,
      numpy-only `metrics.py` module (no pymoo dependency), so the
      spacing-correctness pytest test no longer requires pymoo
      installed to run -- verified directly by importing `metrics.py`
      with plain system `python3` (pymoo confirmed absent) and
      confirming `compute_spacing` still works.
78. **Modern-algorithms result, with a real normalisation bug caught
    and fixed before it shipped.** `run_modern_algorithms_v3.py`'s
    first version built its shared reference front from the
    already-globally-deduplicated 659-point pooled front
    (`global_pareto_v3.npz`) plus SPEA2/AGE-MOEA-II's own fronts --
    this produces DIFFERENT (narrower) normalisation bounds than
    Table 3's own convention (which pools the 150 RAW per-seed fronts
    before global deduplication, since individually-non-dominated
    points that get filtered out of the global front by other runs'
    points still legitimately widen the bounding box). Caught by
    comparing the two bounds directly (rate objective range: 8.01 vs.
    4.44) before trusting the HV numbers; fixed by rebuilding the
    reference from the same 150 raw checkpoints Table 3 itself uses.
    Real, corrected result: SPEA2 $\mathrm{HV}=0.298\pm0.034$,
    AGE-MOEA-II (Panichella, GECCO'22) $\mathrm{HV}=0.298\pm0.031$,
    both statistically indistinguishable from NSGA-II (Mann-Whitney
    $p=0.97$, $0.99$) but significantly ahead of NSGA-III ($p=0.0026$,
    $0.0038$) -- a genuinely 2020s many-objective method does not
    change this paper's central "dominance-based selection wins on
    this low-effective-dimensionality front" finding.
79. **The genuinely independent N=200 instance completed and
    integrated.** All 25 runs (5 algorithms x 5 seeds; reduced from
    the other additional instances' 10 seeds given N=200's much
    higher per-oracle-call cost, disclosed explicitly rather than
    silently) finished cleanly. `run_multi_instance_v3.py`'s generic
    `aggregate()` assumes every instance shares the same seed count,
    so `global_n200`'s summary/pairwise-stats were computed with a
    small standalone script and merged directly into
    `multi_instance_v3_results.json` and
    `multi_instance_pairwise_stats_v3.json` rather than forcing the
    shared aggregator to handle a per-instance seed count it was not
    designed for. Real, honestly-reported result: NSGA-II leads
    NSGA-III ($\mathrm{HV}=0.503$ vs.\ $0.492$) but not significantly
    at $n=5$ (Mann-Whitney $p=0.222$, expected low power at this
    sample size, not evidence the effect vanishes); SMS-EMOA and RVEA
    swap rank (a fourth instance-dependent detail, reported alongside
    the `metro_europe` NSGA-II/NSGA-III reversal rather than only
    that one); MOEA/D remains last. Added 2 regression tests:
    `test_cities_csv_has_200_unique_real_rows` and
    `test_global_n200_oracle_matches_original_100_city_oracle` (18
    tests total). Conclusion's future-work list updated: the
    routing-policy sensitivity check and the N=200 instance, both
    previously open items, are removed from that list now that both
    are done; extending `global_n200` to 10 seeds (parity with the
    other instances) is added as the new, smaller remaining item.
80. **Independent Fable (Claude model) adversarial review of the Round
    8e package, real findings fixed.** Deliberate cross-model
    second-opinion review (an earlier bookkeeping mix-up in this
    project is exactly why a second, independent check is worth
    running before submission). Findings, all verified
    against the actual JSON/checkpoints before fixing, not taken on
    faith:
    - **Stale pre-latency-fix numbers survived in the "symmetric
      IGD_union" paragraph** (Section: Decoder ablation): quoted
      H=0.091/A=0.162/C=0.165/HV=0.258, all pre-fix values; corrected
      to the released `ablation_v3_results.json` values
      (H=0.163/A=0.233/C=0.229/HV=0.148) -- the qualitative claim (H
      closest to the union front) survives, every number was wrong.
    - **Threats to validity item (v) still claimed "the two designs
      agree on the significance pattern for all 30 comparisons"**,
      contradicting the Results/Conclusion's own correctly-stated
      23/30 -- corrected to match.
    - **Stale ablation front count**: "559 points" (pre-fix) -> 392
      (verified directly against `ablation_v3_pareto.npz`'s
      `F_A_full_decoder` shape).
    - **Related Work still described spacing as "a Euclidean-distance
      variant"**, contradicting Methodology's Manhattan/Schott
      description (a Round-6-followup sentence never updated when
      Round 7 fixed the actual metric) -- corrected.
    - **Two floats genuinely overflowed the page**, confirmed via
      `latexmk`'s own "Float too large" warnings AND visual inspection
      of the rendered PDF (Table 1's page number collided with the FSO
      row; Algorithm 1's footer collided with its last line) -- fixed
      with `\footnotesize`+tighter `\arraystretch` (Table 1) and
      `\footnotesize` (Algorithm 1); both warnings confirmed gone on
      recompile, PDF text-extracted and visually re-checked clean.
    - **Modern-algorithms significance not persisted**: the
      manuscript's SPEA2/AGE-MOEA-II p-values were only reproducible
      via an ad hoc script, the same class of gap CHANGELOG item
      76(b) already fixed once for multi-instance -- added a
      `significance_vs_main` field to
      `run_modern_algorithms_v3.py`'s `aggregate()` (Mann-Whitney vs.
      NSGA-II/NSGA-III under the same reference), regenerated the
      JSON, values match the manuscript exactly (0.9705/0.9941/0.002624/0.003848).
    - **Category-native baseline's "$10.9\%$, the decoder's own median
      density" claim is now false**: the post-latency-fix pooled-front
      median density is $13.2\%$ (`global_pareto_v3_summary.json`),
      not $10.9\%$ (the pre-fix value the density grid was never
      updated to match) -- reworded to state the density value used
      without asserting it currently equals the decoder's true median.
    - Minor: ablation-figure caption's "$p_{\mathrm{Holm}}<10^{-7}$
      throughout" corrected to $\leq10^{-5}$ (H is $9.96\times10^{-6}$);
      three `tab:mannwhitney` cells (0.002/0.004/0.003) had only one
      significance star despite being $<0.01$ -- corrected to `**`;
      the sigma-ablation table's $|F^\star|$ bold marker was on
      variant C (63.5) when variant B($\sigma$=500)'s 64.7 is the
      actual column maximum -- corrected; $k$-NN($k=10$) latency
      0.107->0.106 (rounding); Discussion's "Holm-corrected Wilcoxon
      test" reference to the paired-design-as-decision-criterion
      updated to the primary Mann-Whitney test; abstract trimmed by
      ~17 words for a larger safety margin under stricter
      word-counting conventions than this package's own `.split()` test.

81. **Independent Kimi K3 review (via kimi.com, real web session, PDF
    upload) of the Round 9 package, real findings fixed; one Kimi
    finding checked and rejected as a false positive.** Findings
    verified against the actual `.tex`/`.bib`/JSON sources before
    acting, not taken on faith:
    - **Internal development filenames (`.py`/`.json`) leaking into
      the manuscript's flowing body prose** (8 occurrences: coverage-
      threshold, random-search, `global_pareto_v3_provenance.json`,
      the master-comparison script, multi-instance pairwise stats,
      effective-dimensionality, routing-policy, and one figure
      caption) -- moved into footnotes (or, for the figure caption, a
      generic "released with the reproducibility package" phrase),
      preserving the exact traceability an earlier review round asked
      for while removing raw filenames from the prose readers see.
    - **Companion-paper citation (`gilruiz2026surrogate`, the
      QI-Bench-100 paper) was stale**: the `refs.bib` entry named a
      venue and status that were no longer current -- corrected to
      `note={Manuscript in preparation}` with no named venue, since
      claiming an active submission at a venue that
      rejected it would itself be a fabricated-status error.
    - **Three verifiable-but-missing DOIs** on references Kimi
      independently confirmed via web search (`deb2014nsga3`,
      `blank2020pymoo`, `wu2021sequence`) -- added
      `10.1109/TEVC.2013.2281535`, `10.1109/ACCESS.2020.2990567`,
      `10.1088/2058-9565/ac22f6` respectively.
    - **Repeated defensive-framing phrase** ("we report ... directly
      rather than hide/explain away", 4 near-identical occurrences)
      -- reworded 2 of the 4 for prose variety.
    - **New float-too-large warning introduced by this round's own
      footnote edits** (Figure~\ref{fig:pareto-pairwise}, "Float too
      large for page by 10.7pt") -- fixed by reducing the figure to
      `0.94\textwidth`; confirmed gone on recompile. A second,
      separate overfull-footnote warning (two long filenames in one
      footnote) was fixed by splitting into two footnotes.
    - **False positive, checked and rejected**: Kimi flagged a
      "$\mathrm{HV}=0.423$ vs.\ Table 9's $0.433$" numeric
      inconsistency as unexplained. Verified against
      `representation_comparison_master_v3.json`: $0.423$ is the
      correct value for that paragraph's distinct, deliberately
      different shared reference front (503 points, pooling three
      representations), and the manuscript already explains this
      exact distinction in the preceding paragraph -- no change made.
    - **False positive, checked and rejected**: Kimi speculated
      (explicitly hedging that it could only see PDF-extracted text,
      not the rendered PDF) that Algorithm 1 (p.15), Figure 3's
      caption (p.35) and the contributions list (p.5) might have
      LaTeX overflow, orphaned formula fragments, or inconsistent
      numbering. Rendered and visually inspected all three pages
      directly (PyMuPDF): all clean, no overflow, no orphaned
      fragments, numbering correct -- these were PDF-text-extraction
      artefacts, not real defects.

82. **Independent Fable (Claude model) adversarial review of the Round
    10 package, real findings fixed.** A second Fable review, run
    fresh against the Round 10 PDF (not assuming the Kimi-round fixes
    were perfect), cross-checked ~40 groups of numbers against
    `results/*.json`, recompiled the manuscript, inspected the PDF
    with PyMuPDF, and independently verified suspicious references
    online. Verdict: CAMBIOS MENORES. Findings, all independently
    re-verified via `WebSearch`/`WebFetch` before fixing (not taken on
    Fable's word alone):
    - **Three references with fabricated or wrong author names/titles**
      -- the single most serious class of finding, matching this
      author's known past rejection pattern for invented references:
      - `coopmans2025stochastic` (now renamed `avis2025stochastic`):
        the real paper (Phys. Rev. Research 7, 033111, 2025;
        arXiv:2501.06291) is by **Guus Avis and Stefan Krastanov**,
        not "Coopmans, Tim and others" -- confirmed by fetching the
        arXiv abstract directly. The prose at the citing sentence
        also said "Coopmans et al." -- fixed to "Avis and Krastanov"
        in both `refs.bib` and `main.tex`, and the `\cite` key
        renamed throughout (1 occurrence) so the key itself no longer
        misattributes authorship.
      - `satoh2022quisp`: real title is "QuISP: a Quantum Internet
        Simulation Package" (not "an event-driven simulator of
        quantum internet protocols") with 10 real authors (Satoh,
        Hajdušek, Benchasattabuse, Nagayama, Teramoto, Matsuo,
        Metwalli, Satoh, Suzuki, Van Meter) -- confirmed by fetching
        the arXiv:2112.07093 abstract directly. `refs.bib` previously
        listed "Satoh, Takaaki and Mochizuki, Naphan and Boschero,
        Riccardo and others" -- two of those three names do not
        appear on the real paper at all.
      - `geyer2018learning`: real title is "Learning and Generating
        Distributed Routing Protocols Using Graph-Based Deep
        Learning" (not "Learning network design objectives using
        graph neural networks") -- confirmed via web search cross-
        referencing the TUM publication page for the same Big-DAMA'18
        paper. Authors (Geyer, Carle) were already correct.
    - **Four minor bib-metadata corrections**, each independently
      re-verified (not just Fable's claim): `le2022routing`'s
      co-author is Tu N. Nguyen, not "Nguyen, Tuan A" (confirmed via
      the Kennesaw State faculty repository page); `wu2021sequence`'s
      article number is 045027, not 044006 (confirmed via IOPscience
      and a second independent search); `yuan2015theta`'s actual
      publication year is 2016 (TEVC 20(1), confirmed independently);
      `zhang2014knee`'s actual publication year is 2015 (TEVC 19(6),
      confirmed independently); `schott1995spacing` is a Master's
      thesis, not a PhD thesis (`@phdthesis`->`@mastersthesis`).
      `meuser2025reliq`'s claimed year/volume mismatch was checked
      independently and NOT changed -- this session's own web search
      found "2025, volume 74" (conflicting with Fable's "should be
      2026"), and IEEE Xplore's page could not be scraped to settle
      it, so the pre-existing value was left as the more directly
      corroborated one rather than guessing.
    - **Footnote mislabelling a JSON as a "Scoring script"**
      (`main.tex`, the random-search paragraph) -- reworded to
      "Scoring results" since `random_search_vs_decoder_shared_ref_v3.json`
      is a results file, not a script.
    - **Rounding**: the master-comparison paragraph's decoder HV std
      was stated as "$\pm0.034$"; the JSON's true value
      (0.03347333864381025) rounds to $\pm0.033$, not $\pm0.034$ --
      corrected.
    - **Overclaim**: "every summary statistic in this paper is
      reported with both mean and median" (and the near-identical
      contribution-list claim) is not quite true of the manuscript's
      own tables (only HV's median appears in prose; IGD/spacing/
      runtime medians exist in the released JSON but are not
      surfaced in-paper) -- both sentences reworded to accurately
      describe what is actually shown in the paper vs.\ what is
      additionally available in the reproducibility package.
    - **Precision**: "$30$ (or, for the ablations, $10$--$15$)
      stochastic seeds" incorrectly implied all ablations use fewer
      than 30 seeds; the decoder ablation itself uses the full 30
      (confirmed via `run_ablation_v3.py`'s `SEEDS = range(42, 72)`)
      -- reworded to name only the sigma ablation and secondary
      analyses as the reduced-seed-count exceptions.
    - **Precision**: the global\_n200 paragraph quoted
      $\mathrm{HV}=0.503$ vs.\ $0.492$ immediately before citing
      $p=0.222$ from a Mann-Whitney test, but that test was actually
      computed under a *different*, narrower NSGA-II-vs-NSGA-III-only
      reference front (means $0.492$/$0.483$ per
      `multi_instance_pairwise_stats_v3.json`, confirmed by direct
      inspection), not the $0.503$/$0.492$ five-algorithm-front means
      quoted just before it -- reworded to attribute each pair of
      numbers to its own normalisation explicitly.
    - Confirmed clean on recompile: PDF grew from 75 to 76 pages from
      this round's edits (checked for new blank-page/overflow
      artefacts via a PyMuPDF bottom-margin scan across all pages --
      none found), abstract still 230 words, 18/18 tests pass.

**Round 8** (independent external review of the Round-10 PDF, delivered
through a review project with prior history of this same manuscript. Every
finding below was verified against `main.tex`, `refs.bib` or the results
JSONs before being acted on; findings that turned out to be already fixed
in Round 11 or factually wrong about the current text are listed at the end
as discarded):

77. **Werner path search: the published justification was wrong, and the
    correction strengthens the result.** The manuscript stated that Werner
    composition "is neither a minimum nor a simple product" and rested
    optimality on it being "monotonically non-increasing per additional
    hop". Under the substitution `w = (4F-1)/3` the swap composition is an
    exact product, `w12 = w1*w2`, so (a) an exact *additive* shortest-path
    formulation exists with weights `-log w`, (b) the Dijkstra-like
    relaxation is optimal by the standard multiplicative-shortest-path
    argument rather than by hand-waving, and (c) what is a proxy is
    `-log F`, not additivity as such. New script
    `werner_path_validation_v3.py` proves the identity numerically
    (max error 2.2e-16), verifies the oracle against an independent
    `-log w` shortest path on 18,674 reachable pairs of the paper's own
    topologies (max difference 3.3e-16), and verifies both against
    exhaustive enumeration on 200 random 7-node graphs. It also measures
    how often the `-log F` proxy picks a strictly worse path: **18.0% of
    reachable pairs**, where the manuscript said "a small fraction" --
    corrected in the text.
78. **Algorithm 1's edge sampling was under-specified, and the effective
    edge probability is not the nominal one.** The decoder loops over every
    node and its top-k neighbours, drawing one Bernoulli(p) per visit,
    while storing edges in a set keyed on the canonical unordered triple.
    A pair in *both* endpoints' k-NN lists is therefore drawn twice, giving
    it effective inclusion probability `1-(1-p)^2`. New script
    `decoder_edge_probability_v3.py` quantifies this: 48.7-55.4% of
    candidate fibre pairs are mutual k-NN across the range the search
    explores, and a Monte-Carlo check against the real `decode_topology()`
    reproduces both probabilities exactly (at p=0.5: 0.752 observed for
    mutual pairs vs 0.75 predicted; 0.501 vs 0.50 for one-sided). No result
    changes -- the decoder is what it is and generated everything reported
    -- but Algorithm 1 and the surrounding text now state it.
79. **SPEA2 and AGE-MOEA-II integrated into the main comparison, and the
    HV-comparability claim they carried was wrong.**
    `run_modern_algorithms_v3.py` normalises against the union of the 150
    main per-seed fronts *plus its own 60*, a strictly wider bounding box
    than Table 3's 150-front union, while its own description asserted the
    values were "directly comparable to Table 3". New script
    `unified_seven_algorithm_comparison_v3.py` re-scores **all seven**
    algorithms against one reference front built from all 210 raw per-seed
    fronts, and runs a joint Kruskal-Wallis + Holm-corrected pairwise
    family over 3 metrics x 21 pairs = 63 tests, plus the paired
    Friedman/Wilcoxon sensitivity analysis at the same family size.
    **This changed the paper's headline.** The raw p-values are unchanged,
    but with the family at 63 tests instead of 30, NSGA-II vs NSGA-III is
    no longer significant on any metric (`p_Holm` = 0.060 / 0.058 / 0.195,
    vs `<=0.045` under the five-algorithm family). The defensible result is
    a top group -- NSGA-II, SPEA2 and AGE-MOEA-II mutually
    indistinguishable, NSGA-III separated from neither them nor SMS-EMOA,
    RVEA and MOEA/D clearly behind. Abstract, Tables 3-7, Figure 1, the
    discussion, threat (v) and the conclusion were all rewritten to match.
80. **Multi-instance p-values were reported without a family correction.**
    The four per-instance NSGA-II-vs-NSGA-III tests are one family;
    `multi_instance_pairwise_stats_v3.py` now Holm-corrects them jointly
    and adds Cliff's delta plus a 20,000-resample bootstrap CI on the mean
    HV difference. Consequence: `regional_eu_na` (p=0.017) does **not**
    survive correction (`p_Holm`=0.052) and is now reported as suggestive,
    not as confirming the main benchmark's pattern; only `global_n50`
    survives. `global_n200`'s CI (`[-0.088, +0.079]`) is now reported
    instead of describing its non-significance as "expected low power".
81. **The front's connectivity was never reported, and most of it is not
    connected.** Because fidelity/rate/latency are averaged only over
    connected pairs, a topology can drop its hard pairs from three
    objectives and pay only in coverage. New script
    `connectivity_diagnostics_v3.py` computes, for all 659 front points,
    component counts, giant-component size and reachable-pair fraction on
    the same adjacency the oracle routes over. Only **28.1%** connect all
    100 cities; the median point has 10 components, a giant component of
    43 cities and reaches 27.1% of city pairs; the three points with <20
    edges reach 0.4% of pairs across 83 components. Reported in a new
    subsection, and the conclusion now says the released front must be
    filtered by connectivity or minimum coverage before use.
82. **Sigma ablation overclaim.** "Matching the best individually-tested
    fixed value" was not supportable: variant A's mean HV (0.290) is
    *below* sigma=500's (0.312) and above sigma=0's (0.271), and is
    statistically indistinguishable from all three at n=15. Reworded to
    claim only that optimising sigma removes the need to choose a
    constant, at no detectable cost. Threat (vii) now states plainly that
    the methodological objection to optimising a realisation index stands
    unrefuted, and names the clean alternative designs.
83. **Complete link-level oracle model added as an appendix.** Every
    fidelity, rate and latency equation for fibre, satellite and FSO, plus
    satellite link geometry, duty cycle, transmittance, losses, units,
    numerical conventions and edge cases, so the oracle can be audited
    without the code or the companion paper. "Exact analytical BDCZ
    oracle" is now "exact under the stated Werner-composition model"
    throughout.
84. **Reference metadata.** `li2014diversity` (Li, Yang & Liu, *Diversity
    comparison of Pareto front approximations*) was cited as "the unified
    r-metric"; it is a diversity-comparison study, and the R-Metric is a
    different paper by Li, Deb & Yao that this work does not use -- the
    sentence was rewritten and the r-metric mention removed.
    `schoute2016shortcuts` (routing shortcuts) and `vardoyan2021quantum`
    (stochastic analysis of an entanglement-distribution switch) were
    cited for "minimum-cost placement", which neither addresses -- the
    sentence now describes what each actually does. The companion-paper
    entry rendered as "topology evaluationManuscript in preparation" --
    fixed. GeoNames now has a proper bibliography entry with download date
    and CC BY attribution, as its licence requires.
85. **Scope of the SeQUeNCe check restated.** It is a directional,
    single-span, fibre-only loss-vs-distance sanity check; the text now
    says so and enumerates what is outside it (magnitudes, fidelity,
    swapping, multi-hop, satellite, FSO). The phrase "oracle
    cross-validation" is no longer used for it.
86. **Routing-policy sensitivity scope restated.** It re-evaluates the 659
    discovered points under alternative policies; it does not re-optimise
    under them, so the claim is now limited to "the canonical front
    remains low-dimensional when re-scored".
87. **Correlational claims demoted to hypotheses.** The coverage/rate and
    fidelity/rate mechanisms in the front-structure section, and the
    effective-dimensionality explanation in the discussion, are now stated
    as hypotheses consistent with the data, with the untested alternatives
    named. The bimodal-seed overlap is described observationally: with the
    same seed, population size and LHS sampler, NSGA-II and NSGA-III can
    receive near-identical initial populations, so the overlap may reflect
    shared starting points rather than a shared attractor.
88. **Layout and figures.** Float-placement parameters relaxed
    (`\topfraction` etc.), which removed the pages holding a single small
    figure surrounded by whitespace -- the worst intra-page gap dropped
    from 30% of the text height to none above 25% except the final page.
    All figures widened to `\textwidth` and their fonts sized up to
    survive the ~0.5x downscale into the text block (the 16-panel matrix
    in Figure 3 was the specific complaint). Figure 1 now shows all seven
    algorithms from the unified JSON; Figure 3's legend is derived from
    the front's own provenance so it no longer advertises the two
    algorithms that contribute no points to it.
89. **`generate_topology_maps_v3.py` was broken.** It sized its coordinate
    arrays at `N_NODES=100` but read all of `cities.csv`, which has held
    200 rows since the `global_n200` instance was added -- an
    `IndexError` on every run. Rows beyond the main benchmark are now
    skipped. Regenerated output is unchanged.
90. **Style and length.** Anticipation-objection-self-defence constructions
    removed ("a reviewer may reasonably ask", "rather than rationalised
    away", "not as noise to be explained away", "not a hidden one"), the
    conclusion rewritten from a results recap into four robust findings
    plus two bounding limitations, and a new discussion paragraph states
    what transfers beyond quantum networks (representation vs optimiser,
    effective vs nominal dimensionality, realisation indices as
    representation defects, declaring correction families).

*Discarded after verification* (raised by the review, checked, not real
against the current manuscript): the `avis2025stochastic` author
attribution and the QuISP title/authorship were already corrected in
Round 11 and are right in the current PDF; the volume years for
`yuan2016dominance` (20(1), 2016) and `zhang2015knee` (19(6), 2015) are
already correct; the claim that the manuscript says "for the ablations,
10-15 seeds" misreads a sentence that says "10--30 seeds", though that
sentence has been replaced by an explicit per-experiment enumeration
anyway. One finding (a density-homogeneity concern about Figure 4) could
not be acted on because the reviewer's text was truncated mid-sentence.

**Round 10** (package-numbering "Round 13" -- a critical physics bug found by
an independent review, confirmed and fixed, and every downstream result
recomputed from scratch: 925 evolutionary runs, all statistics, all figures.
This is the largest single revision in this project's history and is
documented at length here because it changed the paper's headline result,
twice, in opposite directions, within the same round.):

91. **CRITICAL, confirmed real: the satellite model allowed links with no
    line of sight.** The oracle placed a single satellite at h=500 km over
    the midpoint of each candidate link and computed slant range as
    sqrt((d/2)^2+h^2) (a flat-Earth approximation) with NO visibility
    constraint, while the decoder's own search bounds allowed satellite
    links up to 15,000 km. The true maximum ground separation at which a
    single satellite at that altitude can see both stations is
    2*R*arccos(R/(R+h)) = 4,891 km (R=6,371 km Earth radius); beyond it the
    Earth blocks the line of sight. Measured on the Round-12 front: 14.6% of
    all satellite edges (69,800 of 476,827, spanning 382 distinct city
    pairs) were physically impossible, in a front that was 87.9% satellite
    by pooled edge count. The companion QI-Bench-100 paper's oracle has the
    identical bug; its published Zenodo dataset was independently confirmed
    to have 4.4% of satellite edges (17,019 of 384,181) affecting 57.2% of
    its 1,000 released topologies. Per the user's explicit decision, the
    companion's dataset/manuscript will be corrected together with its
    first round of peer review, not proactively.
92. **Two more real bugs in the same satellite model, found and fixed
    alongside the geometry.** (a) The link-budget efficiency squared only
    the atmospheric transmittance while applying the geometric, pointing
    and detector terms once -- inconsistent under any reading of those
    terms, since heralding a pair requires detection at both stations; all
    per-arm factors are now squared. (b) Latency summed both downlink
    transits (2s/c) when the two photons are emitted simultaneously over
    equal-length arms and so arrive after a single transit (s/c).
93. **A related, smaller bug: an unjustified rate floor.** `oracle_evaluate`
    floored every per-hop rate at 1e-3 s^-1 before taking the bottleneck,
    described as a division-by-zero guard. No reported quantity divides by
    the rate, and every value is representable in double precision; the
    floor silently inflated long-haul fibre rates by up to twelve orders of
    magnitude (1.5e-12 s^-1 at 3,000 km reported as 1e-3). Removed.
94. **Full recomputation: 925 evolutionary runs** across all nine
    experimental blocks (main comparison, decoder ablation, SPEA2/AGE-MOEA-II,
    sigma ablation, convergence, multi-instance x4, random search, binary
    and category-native categorical baselines, classical constructors),
    plus every derived analysis (global front, unified 7-algorithm stats,
    connectivity, PCA, sensitivities, Werner validation, edge probability).
    New guard script `verify_aggregates_fresh.py` compares each aggregate
    JSON's mtime against its checkpoint directory's newest file and
    regenerates any that are stale; added after this exact failure mode
    nearly shipped wrong sigma-ablation numbers (the checkpoints reran
    correctly but nobody re-called `--aggregate-only`, so the JSON silently
    kept 2-day-old contents that looked like "no change from the fix",
    which would have been a false negative reported as a real finding).
95. **The headline algorithm result flipped, and flipped again.** Round 12
    (previous round) found NSGA-II vs NSGA-III NOT significant under the
    7-algorithm Holm family (p_Holm=0.060). Under the corrected oracle it
    IS significant under both the 7-algorithm family (p_Holm=0.001) and a
    5-algorithm-only family (p_Holm=0.0004-0.001) -- i.e. it no longer
    depends on which family is declared, unlike the Round-12 result, which
    sat close enough to threshold to flip with family size alone. Top group
    is now {NSGA-II, SPEA2, AGE-MOEA-II}, each significantly ahead of
    NSGA-III/SMS-EMOA/RVEA/MOEA-D on HV; NSGA-III and SMS-EMOA remain
    mutually indistinguishable. SMS-EMOA's standing improved sharply
    (cardinality 34->81.5, no longer "underperforms on every indicator");
    a plausible untested hypothesis is that its hypervolume-contribution
    selection specifically benefited from the removal of artificially
    attractive long-range satellite edges.
96. **The decoder ablation's central finding survived essentially
    unchanged**: fibre necessary (p_Holm<1.4e-8 for every fibre-free
    variant), satellite's marginal contribution beyond fibre still not
    statistically detectable (p_Holm=1.000). Same qualitative story as
    Round 12, different (higher) HV numbers throughout.
97. **The representation-vs-optimiser gap shrank by an order of
    magnitude**: the "roughly 58x" HV gap between the parametric decoder
    and the relaxed-categorical (SBX/PM) baseline, quoted since early
    rounds, is now roughly 5x. The categorical baseline's naive per-pair
    search had been benefiting almost as much as the k-NN decoder from the
    uncorrected model's unrealistically long-range satellite links; capping
    range at the physical visibility limit removed most of that shared
    advantage. The qualitative claim (representation matters far more than
    optimiser choice) survives; the specific factor does not transfer
    version-to-version and is now stated with that caveat everywhere it
    appears (abstract, discussion, categorical-baseline section).
98. **The category-native operator baseline's density-vs-HV relationship
    inverted.** Previously the best of four tested initial edge densities
    was 15% (extended to 30 seeds on that basis). Under the corrected
    oracle, 3% density is now dramatically the best (own-reference HV
    0.252 vs 0.017-0.028 for 10%/10.9%/15%) -- a denser initial population
    now wastes more of its edges on now-infeasible long-range satellite
    links than the add/remove/change-medium mutation operators can
    productively fix within budget. Extended d03 (not d15) to 30 seeds to
    match statistical power, and updated `representation_comparison_master_v3.py`
    to use it as the primary category-native row (was hardcoded to d15).
    Consequence: category-native operators at their best density now
    slightly *outperform* the relaxed-categorical SBX/PM baseline under
    the shared master-comparison reference, a reversal from Round 12.
99. **The classical-baseline non-dominated exemplar changed family
    entirely.** All three hub-and-spoke variants -- the non-dominated
    classical constructors since early rounds -- are now dominated: their
    coverage advantage depended on hub nodes reaching every other city
    directly, mostly via long-range satellite, and capping satellite range
    removes most of those links (15-hub coverage collapses from 0.559 to
    0.177). The Random Geometric Graph at a 5,000 km threshold is now the
    sole non-dominated classical constructor (coverage 0.278, narrowly
    above the evolutionary front's own max of 0.263). `generate_topology_maps_v3.py`
    panel (d) now plots this RGG instead of the 10-hub network.
100. **The front's connectivity result got both better and worse.** Median
    reachable-pair fraction rose (27.1% -> 54.8%) and the pathological,
    near-fully-disconnected extreme disappeared entirely (3 points with
    <20 edges -> 0), but the fraction of points that are FULLY connected
    (all 100 cities, 1 component) dropped from 28.1% to exactly 0%: no
    point in the corrected front connects every city. This is read as a
    genuine finding (a single-satellite architecture does not produce a
    fully connected global network on this benchmark under physically
    correct visibility), not an artefact, and is now reported as such
    rather than as "roughly a quarter of the front is fully global".
101. **Multi-instance: global_n200 extended from 5 to 10 seeds** (matching
    the other three additional instances; previously listed as future
    work, done as a side effect of re-running the whole block). The
    NSGA-II/NSGA-III reversal previously confined to the smallest instance
    (metro_europe) now also appears, at comparable effect size, on the
    largest (global_n200); neither is statistically significant, and only
    regional_eu_na survives the joint 4-instance Holm correction.
102. **Sigma ablation's finding got materially stronger**, not just
    renumbered: previously "indistinguishable from all three individually-
    tested fixed constants" (a purely negative/hedged result); now
    optimising sigma is indistinguishable from the two GOOD fixed constants
    (42, 500) and significantly BETTER than the bad one (0,
    p_Holm=0.015) -- a real positive finding (parity with good outcomes,
    protection from a bad one) rather than only "no detectable cost".
103. **New script** `random_search_vs_decoder_shared_ref_v3.py`: the
    manuscript's "NSGA-II vs random search, shared reference" comparison
    was the last remaining number in this paper with no generating script
    (computed by hand in an early round); it silently kept its Round-12
    value through the recomputation until this was noticed and fixed. The
    finding itself changed substantively: random search's share of
    NSGA-II's HV rose from 83% to 92%.
104. Manuscript-wide sweep for every downstream consequence: Table 3
    (7-algorithm comparison) and Tables 4-7 (omnibus/pairwise tests)
    regenerated from `generate_seven_algorithm_tables_v3.py`; Table 8
    (ablation), Table 9 (categorical baseline), Table 10 (classical
    baselines), Table (sigma ablation) all hand-updated against their
    JSONs; abstract rewritten and re-verified at exactly 250 words;
    conclusion, discussion ("why NSGA-II leads NSGA-III", "why
    dominance-and-density algorithms...", "what transfers beyond quantum
    networks") and the multi-instance-validation section rewritten;
    Appendix~A (oracle model) rewritten with the corrected satellite
    geometry, visibility limit and two-arm efficiency; Table 2's
    `max_d_sat` bound corrected to the physical 4,891 km limit; HV
    reference-point sensitivity, coverage-threshold sensitivity and Werner
    proxy-suboptimality percentages re-verified and updated (all close to
    but not identical to their Round-12 values, confirming they are
    genuinely sensitive to the oracle rather than hardcoded).
105. Confirmed unaffected by the satellite fix, and left untouched after
    verification: `fiber_link()` output is bit-for-bit identical
    pre/post-fix, so the SeQUeNCe cross-check remains valid without
    re-running; the mutual-k-NN edge-probability census
    (decoder_edge_probability_v3.json) is unchanged (decoder logic, not
    oracle physics); PC1 explained variance is essentially unchanged
    (95.6% both before and after, to 3 significant figures) -- the
    low-effective-dimensionality finding is if anything more robust now.

**Round 11** (package-numbering "Round 14" -- independent external review of
the Round 13 PDF found that the satellite-visibility fix (Round 10/13) had
not been propagated everywhere: the manuscript's headline algorithm
comparison had moved to all seven algorithms, but the released pooled
Pareto front, and everything computed from it, was still built from only
the original five. Root-caused to a single script and fixed there, then
every downstream consumer re-run and every dependent manuscript number
swept and corrected.):
106. **Root cause: `build_global_pareto.py` still pooled only five
    algorithms.** `unified_seven_algorithm_comparison_v3.py` (the actual
    source of Tables 3-7 since Round 13) scores all seven algorithms
    against a 561-point reference front built from the union of all 210
    per-seed fronts, but the separately-maintained `build_global_pareto.py`
    -- the source of `global_pareto_v3.npz`, the artefact every other
    downstream script and the released package treat as "the" pooled
    front -- still hardcoded `ALGORITHMS = ['nsga2', 'nsga3', 'moead',
    'sms_emoa', 'rvea']`. Fixed to pool all seven with correct per-algorithm
    checkpoint-directory routing (SPEA2/AGE-MOEA-II live in
    `step9_modern_algorithms_v3`, not `step1_baselines_v3`); re-run output
    is exactly 561 unique objective vectors (944 unique topologies),
    matching the unified script's reference front size exactly -- the
    intended cross-check.
107. **Every downstream consumer of `global_pareto_v3.npz` re-run against
    the corrected 561-point front**: `connectivity_diagnostics_v3.py`,
    `coverage_threshold_sensitivity_v3.py`, `routing_policy_sensitivity_v3.py`,
    `dominance_check_v3.py`, `analyze_effective_dimensionality.py`,
    `generate_figures_v3.py`, `generate_tables_v3.py`,
    `generate_topology_maps_v3.py`. Headline conclusions survive unchanged
    (RGG at 5,000 km remains the sole non-dominated classical constructor;
    0% of front points fully connect all 100 cities), but every number
    shifted slightly: PC1 explained variance 95.6% -> 96.1% (correcting
    item 105's "unchanged" claim, which was computed on the stale 513-point
    front); median connected components 3 -> 5; median reachable-pair
    fraction 54.8% -> 50.9%; pooled satellite edge fraction 87.2% -> 87.1%;
    rate/latency, rate/coverage and fidelity/rate correlations each shifted
    in the third decimal.
108. **`generate_figures_v3.py`'s pairwise-projection figure had a second,
    independent bug**: its legend/colour-mapping loop hardcoded the same
    five-algorithm list, so once the front carried real SPEA2/AGE-MOEA-II
    points, those points matched no category in the scatter loop and were
    silently dropped from the figure entirely (not merely mislabelled).
    Fixed to derive the algorithm list from what is actually present in
    the data.
109. **Figure 4 conflated raw and effective edge counts across panels.**
    Panels (a)-(c) (evolutionary solutions) reported `|E|` as the raw,
    medium-specific count (a pair with both fibre and satellite counts
    twice), while panel (d) (RGG, one medium per pair by construction)
    reported a count that was simultaneously raw and effective -- making
    the panels look directly comparable on density when they were not
    (e.g. the max-coverage solution's raw |E|=1,755 vs. RGG's 1,356 looked
    denser; the effective/unique-pair counts are 1,221 vs. 1,356, the
    other way round). `generate_topology_maps_v3.py` now reports
    `n_edges_effective` uniformly across all four panels and captions.
110. **Manuscript-wide numeric sweep**: every occurrence of the stale
    513-point/852-topology front (pooled-front size, exclusive-contribution
    breakdown by algorithm, resource-use table, connectivity paragraph,
    classical-baseline dominance check, PCA/effective-dimensionality
    section, routing-policy-sensitivity PC1 values, threats-to-validity
    items iv/vi/viii, conclusion) corrected to 561/944 and re-verified
    against the regenerated JSONs rather than hand-adjusted.
111. **Threats-to-validity item (v) and the NSGA-III-vs-SMS-EMOA spacing
    claim were still describing the pre-recomputation (Round 12) result.**
    Item (v) said the two statistical designs agreed on 57/63 comparisons
    and that NSGA-II vs. NSGA-III was not significant; the actual Round 13
    tables show 60/63 agreement and NSGA-II vs. NSGA-III significant on
    all three metrics under both designs. Separately, the spacing
    comparison's sign was mis-read: Cliff's delta=-0.52 for "NSGA-III vs.
    SMS-EMOA" means NSGA-III's values are lower (better, since spacing is
    minimised) more often, but the manuscript said SMS-EMOA was
    significantly better. Both corrected; the paired design does not
    confirm the spacing difference either way, so the pair is now
    described as statistically indistinguishable rather than ranked.
112. **Figure 2 (convergence)'s normalisation and IGD reference were
    pinned to the wrong "main comparison".** `run_convergence_v3.py`
    sourced its min-max bounds from `baselines_comparison_v3.json` (the
    original five-algorithm table) and its IGD reference from
    `global_pareto_v3.npz` while that file still pooled five algorithms --
    both stale relative to the manuscript's actual main comparison since
    Round 13 (`unified_seven_algorithm_comparison_v3.json`, 561-point
    reference). Figure 2's IGD values were therefore on a different scale
    than the text describing them. Fixed and re-run (10 seeds x 5
    algorithms, sharded across processes); the underlying population
    dynamics are seeded and unaffected, only the logged metric scale
    changed.
113. **Two contradictory causal explanations for the same historical
    58x-to-5x factor shrinkage** (results section: categorical baseline
    "could exploit [unrealistic satellite links] about as effectively as"
    the decoder; discussion section: categorical baseline "benefited
    disproportionately little" from them) were reconciled against the
    actual direction of the effect: the gap SHRANK once the unrealistic
    links were removed, which is only consistent with the decoder having
    had the larger, now-removed advantage -- i.e. the discussion section's
    framing was correct and the results section's was backwards. Fixed to
    match.
114. **Several leftover "five algorithm" / "the other four" references**
    in sections describing the seven-algorithm main comparison (the
    "exactly two experimental blocks" definition, MOEA/D's budget-matching
    sentence, Table 9's cross-reference to Table 3, a correlation-structure
    sentence) corrected to seven / six / seven throughout.
115. **A stale table cross-reference and a wrong section reference**:
    the category-native-operators paragraph claimed to share "the same 10
    seeds as Table [categorical-baseline]", but that table uses 30; the
    reference removed since no table actually shares that specific 10-seed
    run. The classical-baselines table's caption claimed to report "the
    best-performing exemplar of each family" while both the table caption
    and the following prose say it lists all 14 constructors; fixed to say
    so, and its cross-reference corrected from the Pareto-front-structure
    section to the classical-baselines section where the constructors are
    actually defined.
116. **A future-work item referenced its own resolution by name**: a
    "which the conclusion of an earlier version of this paper listed as
    future work" aside (about extending global_n200 to 10 seeds, already
    done) and an "than an earlier version of this paper reported" aside
    were removed in favour of stating the current fact plainly, per the
    general principle of keeping round-by-round narration out of the
    submitted text and in this changelog instead.
117. **Multi-instance ranking claim in the abstract overstated stability**:
    "the ranking largely holds on four additional instances" did not
    reflect that those instances tested only the original five algorithms
    (never SPEA2/AGE-MOEA-II) and that SMS-EMOA's and RVEA's ranks are
    not stable across instances (Section 6.7's own text already said so).
    Abstract now states MOEA/D-last-throughout as the one stable pattern
    and names the ranking as instance-dependent otherwise, matching the
    body text; re-verified at 249/250 words.
118. **Representation-factor headline used the weaker of two categorical
    baselines.** The discussion section's "changes hypervolume by a factor
    of roughly 5.3" cites the gap to the relaxed-categorical SBX/PM
    baseline; the gap to the best-performing categorical baseline
    (category-native operators, established elsewhere in the same
    manuscript) is 4.4x. Now reports both, 4.4x first.
119. **"Fibre is necessary" was stated as a property of the decoder**,
    but both the medium ablation and the sigma ablation are run with a
    single algorithm, NSGA-III, which the Round 13 recomputation separated
    from the top group (NSGA-II/SPEA2/AGE-MOEA-II) on HV and IGD. All four
    instances of this claim (abstract, ablation interpretation paragraph,
    conclusion item iii) now explicitly condition it on NSGA-III and note
    it has not been confirmed for the other algorithms.
120. **Two new honest-limitation paragraphs added to Appendix A**, both
    quantified from the actual front/city-pair data rather than asserted:
    (a) the satellite link budget holds every loss term constant up to
    the visibility limit and drops the link entirely beyond it, so links
    near the 4,891 km horizon are more optimistic than links near the
    sub-satellite point; 7.8% of visible city pairs (4.9% of satellite
    edges in the topology pool) sit in the outer 10% of that range. (b)
    FSO's 150 km hard limit has no explicit curvature/elevation check; at
    that range a symmetric unobstructed path needs ~0.44 km of terminal
    height, unmodelled here. Only 11 of 4,950 city pairs fall within FSO
    range at all, so this affects a small, enumerable subset.
121. Also fixed while auditing the above: a redundant trailing period in
    the QI-Bench-100 companion paper's bib `title` field was colliding
    with elsarticle-num's own punctuation, rendering as
    "Evaluation.Manuscript" with no space in the compiled bibliography.

**Round 11, second pass** (independent external review of the Round 14 PDF
found that the structural fix above propagated correctly to the running
prose but had NOT been synced everywhere the same numbers are printed
verbatim -- table cells, figure captions, and a package-description
section are not re-derived from the same source at compile time, so each
had to be swept by hand.):
122. **Data and Code Availability section still enumerated five
    algorithms** for the released package, contradicting contribution 5
    and Sections 5/6.2 which define the artefact as the seven-algorithm,
    561-point/944-topology union. Fixed, and cross-checked directly
    against the shipped zip's `global_pareto_v3.npz` (561 points, all
    seven `source_algorithm` values present) rather than trusting the
    prose alone.
123. **Table 10's final row (front-wide objective ranges) still printed
    the old 513-point front's numbers** (41-1,221 / 0.612-0.831 /
    3.52-9.26 / 0.006-0.063 / 0.016-0.263) even though the prose
    paragraph describing the same ranges had already been corrected --
    the table cell is separate LaTeX content, not derived from the
    prose. Corrected to 36-1,221 / 0.612-0.866 / 3.52-9.70 / 0.004-0.063
    / 0.013-0.263.
124. **Figure 4's caption was never synced to the regenerated panel
    labels.** The image itself was correct (max-fidelity F=0.866,
    |E|eff=36; max-coverage |E|eff=1,221; knee |E|eff=709; RGG
    |E|eff=1,356), but the caption prose still quoted the old raw counts
    (F=0.831, E=42/1,755/882) and mixed raw/effective definitions.
    Synced to the effective count uniformly across all four panels.
125. **The spacing-sign fix (item 111) only propagated to two of four
    locations.** Two more instances of the inverted claim ("NSGA-III is
    significantly behind SMS-EMOA" in the results section; "SMS-EMOA...
    significantly ahead of it on spacing" in the discussion) survived.
    Fixed both, plus tightened the abstract and conclusion to state the
    accurate nuance (indistinguishable in HV/IGD; primary design favours
    NSGA-III on spacing, not confirmed by the paired design) instead of
    only "indistinguishable."
126. **A second, different "comparable" timing claim** (distinct from
    the one already correct in the edge-and-medium baseline paragraph):
    "per-run time for this 4,950-variable search is comparable to the
    11-variable decoder's own runs" preceded Table 9, whose own numbers
    show the categorical search is 6.5-7.7x *faster* (0.63-0.71 vs.
    4.56-4.87 minutes), not comparable. Fixed to state the real
    magnitude and separate it from the actual reason 30 seeds are used
    (statistical-power parity, not equal cost).
127. **The 58x-to-gap-shrinkage causal explanation was stated as
    established** in both places it appears, when the experiment
    compares two complete oracle versions that changed four link-budget
    terms simultaneously (visibility, two-arm efficiency, latency, the
    removed rate floor), not an isolated ablation of any one of them.
    Reworded both instances as a hypothesis consistent with the result,
    not a demonstrated cause.
128. **"Pooled-front median density of 13.2%"** no longer matches the
    current front (current median is 971/4,950 = 19.6% unique city-pair
    edges). Reworded to state this honestly: 13.2% was the median at the
    time the 10.9% test density was chosen, before the satellite fix
    changed the front's density, and is no longer close to the current
    value.
129. **"There are exactly two experimental blocks in this paper"** was
    misleading -- sigma, random search, both categorical baselines, the
    master comparison and multi-instance validation each carry their own
    block-specific normalisation too. Reworded to "two *primary* blocks"
    (Tables 3 and 8) plus an explicit list of the narrower-scoped blocks
    and where each is reported.
130. **Ambiguous antecedent in the abstract**: "one classical constructor
    ... remains non-dominated, none of whose points connects all 100
    cities" grammatically attached the 0%-connectivity finding to the
    RGG constructor rather than the 561-point evolutionary front. Split
    into two clear clauses.
131. **A dangling cross-reference**: Table 7's caption pointed to
    "Threats to validity item (v) ... at the previous n=10 design" for
    an explanation that item (v) no longer contains after its Round 11
    rewrite. Made the caption's point self-contained instead of
    depending on a cross-reference that had gone stale.
132. **"Appendix Appendix A" (4 occurrences) and a broken ref [18]
    spacing ("evaluationManuscript") were still visible in the compiled
    PDF** despite `grep` on the `.tex` source finding nothing -- both
    were rendering-level artefacts invisible in source form.
    `elsarticle.cls`'s `\appendixname` already renders as "Appendix " as
    part of `\thesection`, so every `Appendix~\ref{app:...}` in the body
    text was printing the word twice; fixed by removing the manual
    "Appendix" word before all five such references. The [18] issue
    (round 10 item 121's fix) was incomplete: `elsarticle-num.bst`'s
    `@article` handler requires a non-empty `journal` field before
    emitting `format.note`, and this entry (a manuscript in preparation)
    has none, so the separator between title and note was never
    inserted regardless of the title's own trailing punctuation. Fixed
    by changing the entry type to the semantically correct `@unpublished`
    (no `journal` field required, and its handler chains format.title
    into format.note the same way every other correctly-punctuated
    entry in this bibliography does).

**Round 11, third pass** (independent external review of the Round 14v2
PDF: verdict improved from CAMBIOS MAYORES to CAMBIOS MENORES, no critical/
high/medium findings remained; four BAJA editorial residuals identified,
all fixed here.):
133. A fifth `Appendix~\ref{app:...}` doubling instance survived item 132's
    fix because it used the abbreviated form `App.~\ref{app:satellite}`
    (rendering "App. Appendix A.2") rather than the `Appendix~\ref{...}`
    pattern the earlier regex swept. Fixed.
134. The paired-design spacing p-value for NSGA-III vs. SMS-EMOA cited
    `Table~\ref{tab:mannwhitney}` (the unpaired table) instead of
    `Table~\ref{tab:wilcoxon}` (the paired table the p-value actually
    comes from). Fixed, and tightened the surrounding phrase to the
    same "design-dependent and unresolved" formulation used everywhere
    else, rather than "statistically indistinguishable on spacing."
135. A second, separate "the other four" (algorithms) survived from
    before the main comparison grew to seven -- this one in the
    common-seed pairing assumption paragraph, not the MOEA/D budget
    sentence already fixed in item 114. Reworded to avoid hardcoding a
    count that will drift again if the algorithm roster changes.
136. A stale quoted subsection title, `"Why a dominance-based algorithm
    wins"`, referenced a discussion heading that no longer exists after
    an earlier rewrite (the current heading is about dominance-and-
    density algorithms doing as well as anything else, not a single
    algorithm "winning," which would also contradict the current
    top-group-of-three result). Replaced with a plain description
    instead of a quoted title that can go stale independently of the
    section it names.
