# Scientific follow-up of S4, S6 and S6b

## Decision

The current evidence does **not** establish that L1-RBF/Apricot is a generally better representative sampler than uniform random sampling. It establishes selective redistribution of cell-type representation, including useful DA-eligibility gains for some extremely rare populations and substantial losses for others. Random sampling has better composition fidelity and typical-cell 60D proximity in these benchmarks. A narrower diversity/rare-population claim remains a hypothesis until matched random-DA and established-sketching comparisons are available.

This audit only reads existing scientific inputs. New files in this report directory are audit derivatives. No PeakVI, Apricot, neighbor fitting, DA testing or publication-figure generation was run. Both audit programs ran as compute-node batch jobs 11137396 and 11137397. Figure candidates and production results were left unchanged.

## 1. Benchmark validation and limitations

Validation outputs: [benchmark_validation.tsv](benchmark_validation.tsv), [benchmark_protocol_versions.tsv](benchmark_protocol_versions.tsv), [benchmark_comparisons.tsv](benchmark_comparisons.tsv), [all_saved_metrics.tsv](all_saved_metrics.tsv). Scripts: `publication_figures/scripts/audit_saved_benchmarks.py`, `audit_cached_diagnostics.py`. Existing implementation: `publication_figures/scripts/deferred_benchmarks.py`.

- Complete: 315 runs = 15 production references and 300 random runs; each has all-query and nonselected-query summaries. Exactly seeds 2001–2020 at every atlas/budget. No duplicate metric keys.
- Budgets match exactly within each atlas/budget. Current production IDs map to unique valid atlas positions. All 300 random memberships were reconstructed in memory, without replacement, using `SeedSequence([seed, atlas_index, pct])`. Every reconstructed membership hash and population-count vector matches the stored run. These are distinct RNG streams, not nested budget prefixes. No random ID lists were written by this audit.
- Same canonical `results/peakvi/{atlas}/peakvi_latent.npy`, raw 60D, no scaling/standardization, float32 index buffers, Euclidean NNDescent, 30 graph neighbors, index seed 1705, one thread, exact float64 recomputation to returned candidates. Provenance versions agree: NumPy 2.4.6, pandas 2.3.3, SciPy 1.17.1, PyNNDescent 0.6.0. Paths/shapes are verified; historical byte identity is not provable because saved provenance lacks latent-content hashes.
- Same adaptive accuracy rule, **not identical final search effort**: production uses epsilon 0.5 in 7 runs and 1.0 in 8; random uses 0.5 in 166 and 1.0 in 134. Initial requested epsilon in provenance is not the effective effort; use metrics `effective_epsilon` and `accuracy_attempts`.
- All runs pass >=99% of 128 exact validation queries within 1% distance error plus absolute tolerance 1e-7. Minimum pass fraction is 127/128 = 0.9921875. **Important unresolved tail issue:** maximum recorded sampled relative error is 5.4400 (544% excess distance); the rule permits one badly missed query. This is finite-sample calibration, not proof of globally accurate q99 or every radius classification. Distances remain approximate upper bounds; ANN errors can be method-dependent. Keep the broad median/q90 result qualified and prioritize an independent exact sensitivity check before strong tail claims.
- Full-query rows evaluate every atlas cell, set selected cells to true zero, and include them for both methods. Nonselected rows exclude their respective selected cells. The latter query identities differ between methods, so they are a useful sensitivity analysis, not the same fixed test cohort. A common held-out query cohort would be an additional control. The observed disadvantage persists in nonselected summaries.
- One fixed radius per atlas, reused for all methods/budgets/seeds: Kidney 2.966494; HECA 2.142815; Fetal 2.176918. These are raw within-atlas distance units, not cross-atlas comparable or biological boundaries. Radius is the 90th percentile of nearest-nonself distances at 50,000 fixed sampled atlas queries (seed 1703), using a separately calibrated full-atlas index. Radius-estimation uncertainty is not in seed bands. Source: `publication_figures/benchmarks/radius/*_radius.tsv`.
- Shading is the empirical linearly interpolated 2.5th–97.5th percentiles of 20 random-seed values; not confidence intervals, not proposed-method uncertainty. The caption explicitly says this, but still carries a historical “caption template” heading. Plot legends themselves do not spell out the band meaning. No captions or figures were changed in this audit.
- Every median and band endpoint in `tables/S4_S6_measured_benchmark_summary.tsv` was reproduced from saved run metrics. S6b source counts were checked through reconstruction of all population-count vectors. Original sources for every panel are `benchmarks/full/*_metrics.tsv` and `*_populations.tsv`; plotting definitions in `scripts/plot_random_benchmark_when_available.py`. `tables/S6b_measured_rare_population_counts.tsv` exposes the example-level source rows.

**Report these qualifications before extending the analysis:** unequal adaptive search effort; occasional large calibrated ANN misses; missing latent/feature content fingerprints; historical diagnostics for different selections; and the absence of random DA outputs. None invalidates the exact composition/retention comparisons. No evidence was found of a budget, membership, seed-count or annotation-leakage error in random sampling.

## 2. What S4/S6 quantitatively show

Each cell below is proposed / random-seed median. Coverage uses the fixed within-atlas radius and all atlas queries. All five budgets and both query cohorts are in [benchmark_comparisons.tsv](benchmark_comparisons.tsv).

| Atlas | Budget | Composition TV | q90 nearest distance | Fraction within radius |
|---|---:|---:|---:|---:|
| Kidney | 2% | .1150 / .0394 | 4.653 / 4.308 | .144 / .438 |
| Kidney | 10% | .1088 / .0192 | 3.984 / 3.565 | .511 / .721 |
| HECA | 2% | .1420 / .0180 | 3.135 / 2.679 | .395 / .666 |
| HECA | 10% | .1278 / .0075 | 2.743 / 2.397 | .631 / .803 |
| Fetal | 2% | .1191 / .0215 | 3.569 / 2.948 | .230 / .561 |
| Fetal | 10% | .1006 / .0088 | 3.067 / 2.539 | .517 / .767 |

Across all 15 atlas/budget combinations, proposed TV, JS, median distance and q90 distance are worse than **every** measured random run. q99 is worse than every random run in 13/15 combinations, so “random wins every coverage metric everywhere” would be false. With selected cells excluded, proposed q90 remains worse than all random runs and fixed-radius coverage remains lower than all random runs in 15/15 comparisons. These are descriptive comparisons, not replicate-based hypothesis tests of the proposed algorithm. Source: [all_saved_metrics.tsv](all_saved_metrics.tsv), `cached_audit_11137397.log`.

Proposed retains all 24/101/78 annotated types at every budget. At 2%, random retains every type in 20/20 Kidney runs, 5/20 HECA runs, and 15/20 Fetal runs. At 10%, all three atlases have complete retention in all 20 random runs. Thus 2% complete retention is a meaningful observed reliability improvement for HECA, modest for Fetal, and not distinguishing in Kidney; it is not a universal novel capability.

## 3. Rare-population gains and losses

Random entries below are median selected counts; parentheses give the number of random seeds meeting >=10 target and >=20 rest cells. Current proposed counts/enrichments and random percentile/min/max counts for every budget: [rare_population_comparisons.tsv](rare_population_comparisons.tsv).

| Atlas/population | Full count | 2% proposed / random | Proposed enrichment, 2% | 10% proposed / random |
|---|---:|---:|---:|---:|
| HECA stem/progenitor | 56 | 14 / 1 (0/20 eligible) | 12.50× | 47 / 5 (1/20 eligible) |
| HECA naive CD4 T | 219 | 1 / 5 (1/20 eligible) | .228× | 4 / 21.5 (20/20 eligible) |
| Fetal Megakaryocytes? | 60 | 11 / 1 (0/20 eligible) | 9.17× | 20 / 6 (0/20 eligible) |
| Fetal ganglion | 678 | 1 / 11.5 (17/20 eligible) | .0737× | 28 / 68 (20/20 eligible) |
| Kidney Mes | 118 | 4 / 2 (0/20 eligible) | 1.695× | 17 / 12.5 (15/20 eligible) |
| Kidney Macula_Densa | 161 | 3 / 2.5 (0/20 eligible) | .932× | 10 / 14.5 (20/20 eligible) |

The two extremely rare enriched examples become DA-eligible at 2%; random never does in these measured seeds. Conversely, HECA naive CD4 remains ineligible even at 10% whereas every random run is eligible, and Fetal ganglion is ineligible at 2% whereas most random runs are eligible. These are demonstrated eligibility gains/losses, **not demonstrated effect fidelity or discovery gains**. Keep the uncertain Megakaryocytes? annotation exactly as supplied; biological validation of that label remains separate.

## 4. Mechanism: evidence versus conjecture

Current enrichment correlates negatively with full population abundance. At 2%, descriptive across-type Spearman rho is −.624 Kidney, −.175 HECA, −.479 Fetal (full series in [descriptive_associations.tsv](descriptive_associations.tsv)). Abundance alone cannot explain which rare populations gain or lose, as the examples show; enrichment also mathematically contains abundance in its denominator. These are not independent causal effects.

Existing k=30 nonself Euclidean neighbor metrics are saved for Kidney/Fetal. Their cell-order checksums match current atlas names; we re-aggregated them using **current** selected IDs, without building a neighbor index. Latent-byte provenance is not archived, and original kNN diagnostics lack the later exact-calibration gate. Treat the following as descriptive cached-geometry evidence, not exact current-latent validation. Sources: `diagnostics/peakvi_celltype_separability/{kidney_atac,fetal_atac}/per_cell_knn_sparsity.npz`; [cached_cell_metric_selection_summary.tsv](cached_cell_metric_selection_summary.tsv).

| Atlas | Budget | Median cell mean-30NN distance: selected / full atlas |
|---|---:|---:|
| Kidney | 2% | 3.608 / 2.795 |
| Kidney | 10% | 3.229 / 2.795 |
| Fetal | 2% | 2.469 / 1.930 |
| Fetal | 10% | 2.204 / 1.930 |

Current selected cells occupy lower-density regions by this proxy. Yet the across-type correlation of median neighbor distance with enrichment at 2% is −.100 in Kidney and +.488 in Fetal. A uniform “all rare enriched types are sparse” explanation is unsupported. Mes has mean-30NN-distance population median 2.711 versus Macula_Densa 2.024; Fetal Megakaryocytes? 2.166 versus ganglion 1.905. No corresponding HECA cached kNN file was found. kNN purity describes annotation mixing, not the actual distance to the nearest other population. A validated between-population-distance analysis is missing.

Raw activation caches align to the current cell-name order file but lack archived feature-content checksums. Re-aggregation shows current 2% selected-cell median raw RBF mass below atlas median in all three atlases: Kidney 5.702 versus 8.937; HECA 3.352 versus 6.232; Fetal 4.321 versus 8.457. This is consistent with normalization amplifying weak total-similarity profiles, not proof that low mass itself is advantageous. At 2%, enrichment versus population-median raw mass has rho −.602, −.421, −.507 respectively. [cached_cell_metric_selection_summary.tsv](cached_cell_metric_selection_summary.tsv), [cached_geometry_associations.tsv](cached_geometry_associations.tsv).

Own-anchor counts also do not explain enrichment: saved nearest-cell-assigned landmark counts are zero for enriched HECA stem/progenitor and Fetal Megakaryocytes?, but one for depleted naive CD4 and ganglion. Anchors are descriptive nearest-cell labels, not annotation-guided centroids or features exclusive to a type. Conditional normalized residual feature gains for the current subsets are **not** established by these counts.

Cached feature support differs in the focused pairs: effective support 3.721 versus 2.132 for HECA stem/progenitor versus naive CD4; 29.018 versus 30.509 for Fetal Megakaryocytes? versus ganglion. Derived L1 initial gains are 9.027 versus 5.767 for the HECA pair, but 24.487 versus 26.511 for the Fetal pair. Thus neither initial gain nor breadth provides a universal per-population explanation. The derived gain uses saved per-cell raw singleton gain × sqrt(current target mass/raw row mass); it is a float32-derived diagnostic, not a rerun of Apricot. [cached_focus_geometry.tsv](cached_focus_geometry.tsv).

**Historical diagnostic trap:** `diagnostics/rbf_activation_marginal_gain/*` computes residual gains on the **raw** RBF feature matrix, not the current row-normalized matrix. `featurebased_selection_behavior` and related separability tables include historical selected counts: e.g. Fetal Megakaryocytes? 1 at 2%, current 11. Their old enrichment/gain associations are not evidence about current completed selections. Intrinsic geometry/profile summaries can be cautiously reused after alignment; historical selection-conditioned residual gains cannot. [historical_diagnostic_focus.tsv](historical_diagnostic_focus.tsv); diagnostic source `diagnostics/diagnose_rbf_activation_and_marginal_gain.py`.

## 5. What the optimized objective actually rewards

Production `scripts/03_subset_selection_on_peakvi.py` sets x_ij = c r_ij / sum_j r_ij, with c the atlas-median raw RBF row mass, then maximizes F(S) = sum_j sqrt(sum_{i in S} x_ij), using unweighted sqrt FeatureBasedSelection and its default optimizer. Each selected cell contributes the same total feature mass c. There is no explicit annotation, abundance, local-density, nearest-cell-distance or DA term.

Writing C_j = sum_selected x_ij, the marginal gain is sum_j [sqrt(C_j+x_ij)−sqrt(C_j)]. This favors additions to under-covered feature dimensions, but it need not cover cell-density modes or minimize distances to typical cells. At the first step, diffuse normalized rows can have higher sum sqrt(x_ij) than sharply concentrated ones. Normalization removes total-similarity magnitude and can amplify very weak raw profiles; in HECA the recorded 2% row scales range .240868–316.529681 (`results/subset_selection/peakvi/heca_atac/subset_selection_summary_peakvi_02pct.tsv`). Scaling all normalized rows by common c only scales the objective by sqrt(c), leaving mathematical rankings unchanged apart from numerical/optimizer details.

The observed coverage loss is a real disadvantage under a **cell-mass-weighted nearest-distance definition of representativeness**. It is compatible with a feature-diversity objective, but compatibility is not validation of biological usefulness. Conversely, mean/q90 cell-weighted coverage cannot by itself rule out gains in under-represented regions. Those gains require independently measured conditional/local coverage and downstream DA evidence. Optimized feature coverage is a mechanism diagnostic, not independent biological validation.

## 6. Highest-priority staged DA experiment — proposal only

No random DA outputs are present under `publication_figures/benchmarks/`; current S3 contains production versus parent-atlas comparisons only. Existing `scripts/09_1_run_signac_da_all_cell_types.R` uses Seurat objects with BPCells `marker_features(method='wilcoxon')`, not a hidden different Signac test. Preserve the exact implementation: >=10 target, >=20 rest; BH within type across tested matrix peaks before min-detection filtering; detection >=.01 in either group; adjusted P<=.05; |average log2FC|>=.25. The effect uses the existing group-size-dependent +1 pseudocount formula. Do not replace the test with a library default FindMarkers call.

**Stage A: six jobs**, seed 2001, 2% and 10%, all three atlases. This directly covers both enriched/depleted focal examples and both severe and less severe sampling budgets. Use saved full-atlas significant references; no full-atlas DA rerun. Reconstruct the **same** already-benchmarked random memberships, verify their recorded hashes, then write only approved DA-input ID files. Existing full-benchmark runs did not write these ID lists; the current `deferred_random_DA.py` defaults to absent `benchmarks/composition` memberships. The existing pilot batch file therefore must not be submitted unchanged. Prepare verified inputs and pass an explicit memberships directory first; validate raw IDs against the Seurat/BPCells cell names too.

**Stage B: 12 additional jobs**, seeds 2002 and 2003 at the same budgets/atlases, if Stage A demonstrates useful signals and acceptable cost (18 total). Three seeds characterize limited repeatability, not a 20-seed uncertainty distribution. Only then consider all budgets: 45 random DA runs for three seeds. Retain the 20-seed exact composition/eligibility evidence independently; it does not imply 20-seed DA evidence.

**Resources:** initial request per pilot: four CPUs, 128 GB memory, one-hour wall limit, one global concurrent job to serialize BPCells reads. This is a conservative request consistent with the existing pilot design, not measured random-run memory. The six corresponding successful production jobs took 1,343.49 seconds total (22.39 min; .3732 job-hours), with CPU RSS 1.305–8.865 GiB. Their four-CPU allocation corresponds to 1.493 allocated core-hours, not measured CPU utilization. For random jobs, provision 3–5× that total as a planning allowance (~1.1–1.9 serial hours / 4.5–7.5 allocated core-hours), **not a prediction**. Scheduling and postprocessing add time. The initial six one-hour requests cap requested allocation at 24 core-hours if all use their full limit. Eighteen comparable jobs imply 67.17 min at production-equivalent cost, with 3–5× allowance ~3.4–5.6 serial hours; revise from the first actual random job. Source: `tables/S6_DA_cost_reference.tsv`, existing `run_deferred_DA_pilot_AFTER_APPROVAL.sbatch`. Scratch/reference-output size may dominate; record peak RSS, elapsed time, read failures and output sizes. Do not silently lower memory based only on production RSS.

### Predeclare analysis and denominators

- Fixed reference universe R: unique full-atlas-significant (cell_type, peak) pairs, identical across methods and seeds. Freeze direction by full effect >0/open or <0/closed. No refitting or reference choice based on subset results.
- Eligibility: record every annotated type, target/rest counts, success/failure and missingness. Available cohort A: R pairs with finite subset effects. Report |A|/|R| and eligible-type counts, both pooled and by type/direction. Distinguish eligibility from computational failure.
- Spearman/Pearson and median absolute error on A, using unmodified effects. Separately pool signed effects across all pairs and correlate within each type and direction; never substitute an average of type correlations for pooled rho. Predeclare minimum finite/variable pairs (the existing safe correlation requires >=3 and nonzero variation), and report n for every estimate.
- Direction agreement: matching nonzero signs / |A|. Switch rate: opposite nonzero signs / |A|; zeros count in neither numerator. Strict recovery: subset significant AND matching direction / |R|. Report open/closed separately as well as overall.
- Show original analysis and a **fixed within-atlas common-eligibility cohort**, defined by the intersection across both methods, all included seeds and both budgets; preserve the same reference pairs for comparisons. Also use pairwise-common finite effect pairs for directly matched effect correlations, reporting any further missingness. A cohort restricted to production alone is insufficient.
- Report pooled peak-pair weighted values, per-type values/distributions and clearly defined macro summaries. For rare types, include both the gains and losses in this report, including types without computable correlations. Preserve biological annotation uncertainty.
- Seed bands are empirical seed variability. Peaks and cells within a reference population are not independent biological replicates. Parent/subset overlap, one-versus-rest composition and the pseudocount influence effects. Do not attribute differences solely to statistical power, and do not claim external truth from self-recapitulation.

## 7. Prioritized decisions and further approvals

1. **Approve Stage A random DA first**, after preparing hash-verified memberships and an output-only wrapper. Stop on input mismatch or analysis errors. Continue to two additional seeds only after reviewing costs/results. No expensive jobs have been launched by this audit.
2. **Coverage sensitivity control:** independent shared validation queries for representative 2%/10% production/random runs, stratified across dense/sparse regions, with bounded exact searches; review severe sampled misses and radius uncertainty. Current ANN tail gate allows errors too large for unqualified tail claims. Resource/time should be piloted, not extrapolated from a nonexistent exact full-atlas run. Leave figures unchanged pending review.
3. **Current-objective mechanism audit:** stream saved raw RBF features once, use the actual current L1 transform, and calculate current selected feature totals and candidate marginal gains for focal populations versus abundance-matched controls. Reuse existing density caches; obtain missing HECA density and nearest-other-population distances with a separately calibrated bounded analysis. Compare full distributions/quantiles and abundance-adjusted descriptive associations, not merely UMAP positions. Residual gain after the completed selection is not the candidate's gain when it competed during greedy selection; conditional counterfactual additions/removals may be needed. Proposed allocation: 1–2 CPUs, 16–32 GB, initially 1–2 hours for streaming features; HECA neighbor work needs its own measured pilot. These are resource ceilings, not measured runtimes. No such full-feature or neighbor job was run here.
4. **Established sketching comparator:** Geosketch is warranted because it directly targets geometry/diversity, the claimed niche here. Use the same saved 60D matrices, exact budgets and three explicit seeds initially; compare eligibility, composition, local/global distances and the same DA endpoints. No new PeakVI. Treat donor handling and input transforms consistently. Original method: Hie et al., Cell Systems 2019, https://pmc.ncbi.nlm.nih.gov/articles/PMC6597305/ ; implementation https://github.com/brianhie/geosketch . Resources are not measured in this environment; pilot before expansion.

### Manuscript claim assessment

**Demonstrated improvements:** complete retention versus some 2% HECA/Fetal random runs; selective enrichment/eligibility of certain extremely rare annotations; historical L1 ablation improves retention versus the historical unnormalized implementation. These are distinct comparisons, not an overall win over random or established sketchers.

**Random preserves better:** annotated proportions and typical-cell 60D nearest-selected coverage in the tested runs. Exact composition evidence is robust to ANN concerns. **Biological advantage:** plausible for selected enriched populations, unproven until matched DA; depleted populations supply clear counterexamples to a universal rare-type claim.

**Application Note readiness:** this is an editorial/scientific judgment, not a prediction of acceptance. A broad representativeness/biological-fidelity superiority claim is currently unsupported. A competitive narrower feature-diversity/targeted rare-population application is possible only if downstream gains are demonstrated and relevant baseline tradeoffs are reported. Algorithm changes are not automatically required; first test whether this tradeoff is useful. If matched random DA is equal or better across both pooled and per-type endpoints, do not rescue the claim using the objective value or cherry-picked rare examples. Reconsider the objective or manuscript scope. Freeze publication-figure decisions until these tests are reviewed.

## Reproduce this read-only audit

`sbatch publication_figures/audit_saved_benchmarks.sbatch` reads existing results and writes only audit derivatives. The cached diagnostic audit is `publication_figures/scripts/audit_cached_diagnostics.py` and must also run through a compute-node batch job. Neither command runs new DA or creates figures. No need to rerun the atlas-scale benchmark for this report.
