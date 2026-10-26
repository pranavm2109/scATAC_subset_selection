# Deferred analysis plan — NOT EXECUTED

## Scope and available software
Saved production `results/peakvi/{atlas}/peakvi_latent.npy` matrices contain raw 60D Gaussian PeakVI vectors; `cell_names.txt` supplies row identities. No training or Apricot recomputation is needed. Installed in scvi: pynndescent 0.6.0, scikit-learn 1.9.0, scipy 1.17.1. faiss and hnswlib were not found. These are local inventory findings, not assumptions about other environments.

Do not use UMAP coordinates for S4/S6. Do not allocate N-by-k full distance matrices. The supplied script uses a selected-reference PyNNDescent index, full-atlas query chunks of 10,000, and exact distances to returned candidates. Atlas latent arrays are memory-mapped; there is no persistent ANN/data cache. Index and query parameters, versions, selected-ID hash, and coordinate transform (none) are recorded.

## Measurements
1. Primary S4: every full-atlas cell queried against each of the five saved production subsets. Report all-cell and nonselected-cell distance median/q90/q99 separately. The all-cell measure includes true zeros for selected cells, an increasingly important fraction at larger budgets.
2. Accuracy control for each index: 128 reproducibly sampled nonselected query cells, compared to exact brute-force nearest distances using bounded 4,096-reference blocks. If the initial search effort misses the threshold, retry query epsilon at 0.5, 1.0, then 2.0 as needed, reusing the same exact reference distances; record effective epsilon and every accuracy attempt. Abort if the maximum effort still gives fewer than 99% of queries within 1% of exact distance (absolute tolerance 1e-7). This diagnostic is finite-sample evidence, not a guarantee for every query; inspect q99/max error and increase ANN effort if needed. Exact nearest-neighbor identity can differ for tied or almost-tied candidates; distance accuracy is the criterion.
3. Fixed-radius coverage: first calibrate ONE atlas-specific radius, independently of method/budget. Proposed primary radius is the 90th percentile of full-atlas nearest-nonself distance at 50,000 uniformly sampled query cells, using a full-atlas ANN index and separate exact calibration. Query four candidates and exclude self. Use the same radius for all selected methods/budgets/seeds. This defines proximity at a local atlas sampling scale, not a biological boundary. Duplicate/zero-distance behavior must be reviewed. Radius estimation uncertainty is not included in seed-variation intervals. Since ANN candidate distances are recomputed exactly, distances are upper bounds on the true nearest distance; coverage is conservative conditional on this estimated fixed radius. Calibration script adds three full-atlas index tasks. Omitting radius from the pilot remains valid; do not silently substitute a different threshold per subset.
4. S6: independent uniform sampling WITHOUT replacement, exact production cell budgets, seeds 2001–2020. Membership generation uses only N, budget and RNG, never annotations. Independent streams per atlas and budget; random selections are not nested. Existing annotations enter only after membership is fixed. Retain complete per-population counts/proportions/enrichment, a declared rare threshold (<0.1% full abundance), TV and base-2 JS. DA eligibility remains >=10 target, >=20 rest.
5. Compare proposed selection to the empirical random distribution. Shading represents across-seed 2.5th–97.5th percentiles (20 seeds), NOT confidence intervals. No uncertainty band is invented for the single production run. Report all seeds and missing/failed tasks. A measured-results-only plotter refuses incomplete 20-seed comparisons.

## Staging and resources (estimates, not measured NN runtimes)
| Stage | Tasks | Proposed resources | Planning expectation / decision |
|---|---:|---|---|
| Random composition only | 300 random + 15 production references | 2 CPUs, 8–16 GiB, 30–60 min/job allowance | Small count calculations after ID/metadata reads. No latent search. Write IDs only when needed for approved DA. |
| Radius calibration | 3 | 2 CPUs, 32 GiB, 8 h ceiling each | Full-atlas NN index, largest is 1.45M vectors. Unmeasured; pilot before accepting its output. |
| Coverage pilot | 12 | 2 CPUs, 32 GiB, 8 h ceiling; max 2 concurrent | Each atlas at 2% and 10%, production + random seed 2001; shared 50,000 query sample. Planning range: Kidney minutes to ~1 h; Fetal/HECA ~0.5–4 h/task, potentially longer with JIT/index calibration. These ranges are unvalidated. Measure actual wall time/RSS/accuracy before expanding. |
| Full coverage comparison | 315 | Resource limits revised after pilot | 15 production + 300 random indices. Every atlas cell queried. Full-query work is ~1.6x/14.4x/29x the 50k query portion for Kidney/Fetal/HECA; index build time does not scale by that factor. Estimate expansion as sum(index build + full query + calibration), not pilot total times those factors. Could require hundreds of aggregate job-hours; do not launch all tasks before pilot review. |

Raw latent float32 storage is approximately 18.5 MiB Kidney, 165 MiB Fetal, 332 MiB HECA; these exclude metadata, graphs, ANN candidates, JIT and temporary arrays. 32 GiB is a conservative initial request, not a claim that it will be used. Default one ANN thread favors reproducibility; resource/timing comparisons must preserve thread configuration. No saved ANN index cache is introduced.

`tables/S4_S6_DEFERRED_task_manifest.tsv` is metadata only: generating it did not sample cells. Historical unnormalized summary tables exist, but original rules referenced production selected-ID paths that now contain L1 selections. No verified preserved unnormalized ID list was found. New 60D comparisons against that variant require restored verified memberships or separately approved recomputation; the `historical` mode demands an explicit path.

## Proposed approval commands — do not execute until approved
Run from project root. First cheap composition, and separately radius/pilot; inspect metrics before full expansion.

```bash
# Counts-only benchmark, no ANN; existing-output plotting job never runs this.
BENCHMARK_PHASE=composition sbatch --export=ALL --mem=16000M --time=01:00:00 \
  --array=0-314%2 publication_figures/run_deferred_coverage_AFTER_APPROVAL.sbatch

# Radius calibration, if fixed-radius coverage is approved:
sbatch publication_figures/run_deferred_radius_AFTER_APPROVAL.sbatch

# After radius jobs succeed, 12-task NN pilot (50k shared queries):
BENCHMARK_PHASE=pilot sbatch --export=ALL \
  --array=0,1,84,85,105,106,189,190,210,211,294,295%2 \
  publication_figures/run_deferred_coverage_AFTER_APPROVAL.sbatch

# Full-query 315-task expansion requires a SECOND review of pilot cost/accuracy.
# BENCHMARK_PHASE=full sbatch --export=ALL --array=0-314%2 \
#   publication_figures/run_deferred_coverage_AFTER_APPROVAL.sbatch
```

Pilot and full outputs use separate directories. Scripts refuse existing task outputs to avoid accidental overwrites. Use a fresh directory for retries with changed NN settings. The random plotter is not run until all required measured files exist.

## Staged random DA comparison
Reuse the measured full-atlas reference; do not rerun full-atlas DA. First run seed 2001 at 2% and 10% for each atlas (six DA jobs). If useful and affordable, expand to seeds 2001–2003 at all budgets (45 subset DA jobs). Do not use six pilot jobs to imply a 20-seed DA distribution. A second expansion may be warranted after measuring variability.

Production subset DA wall times sum to 0.783 h for one seed-equivalent across all 15 atlas/budget combinations; the six pilot-equivalent jobs sum to 0.373 h. Three seed-equivalents would be 2.35 aggregate job-hours, or 9.40 core-hours at four CPUs, before extra scheduling/I/O/postprocessing. These are EXTRAPOLATIONS from current proposed-method job logs, not measured random-run costs. Budget at least 3–5x timing headroom until random jobs are measured; reference output size and BPCells conversion also matter. `tables/S6_DA_cost_reference.tsv` lists each observed wall time and CPU RSS with exact source paths. Production observed subset RSS spans approximately 0.9–8.9 GiB; an initial 128 GB request preserves the existing production rule's conservative memory budget. No new DA job has been launched.

`deferred_random_DA.py` prepares the exact production DA invocation, original thresholds, existing full-atlas significant reference, and per-task scratch location. It requires previously generated random IDs and --execute on SLURM. Output is exclusively under publication_figures/benchmarks/random_DA. Serialize jobs reading the same atlas BPCells object; initially use global concurrency one or separate serial chains per atlas. Capture Snakemake/sacct peak memory and wall time for the pilot. Request four CPUs, 128 GB, one-hour initial wall limit per pilot job, then revise from evidence. Existing tests use no donor covariate; preserve these definitions for a fair comparison, while documenting cell-level and overlapping-reference limitations. Postprocess random references using the same eligibility/count definitions as S3; do not fabricate pooled random Spearman from aggregate type correlations.

No S4/S6 publication result or comparative advantage claim is generated until the measurements exist.

An unsubmitted six-job pilot file is provided: `run_deferred_DA_pilot_AFTER_APPROVAL.sbatch`. Its global array concurrency is one. It requires the approved composition stage to have written random IDs; it is not a dependency or automatic continuation of any existing-output plotting job. Submit only after explicit DA-pilot approval. Seed-expansion and random DA plotting should be reviewed after pilot outputs exist.

## Array-input regression fix
Radius array 11136447 failed before calibration because a read-only memory-map view reached PyNNDescent/Numba. Index and query inputs now use private writable C-contiguous float32 copies; saved latent files remain read-only. Compute-node regression job 11136611 passed radius/nonself and subset-reference checks on 512 existing Kidney vectors, each with 128 exact validation queries. Radius passed at epsilon 0.2; subset required epsilon 0.5. This is a software check, not an atlas-scale accuracy or benchmark result. Dependent pilot array 11136448 must be cancelled and resubmitted against a new radius job ID.
