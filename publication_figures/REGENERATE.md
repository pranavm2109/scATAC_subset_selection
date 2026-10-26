# One figure location

All publication figure exports now go to `publication_figures/figures/`. Sources are in `scripts/`, measured summary tables in `tables/`, captions in `captions/`. PDF proofs are inside `figures/proofs/`. The duplicate `polished/` and `final_minimal/` directories have been removed after migration and successful regeneration.

From the project directory:

```bash
sbatch publication_figures/regenerate_all_figures.sbatch
```

This regenerates S1, S2, S3, S3b, S5, S7 and S8 in PDF, SVG and 600-dpi PNG, validates cohort denominators and historical selected-ID membership, checks PDF width/fonts, and renders proofs. Final caption qualifications are restored from `captions/templates/`. No training, selection or DA testing is run. Legacy existing-figure batch entrypoints delegate to this same script. S4/S6 are additionally plotted only when all 315 full-benchmark metric files exist and the measured-results plotter validates completeness. An incomplete benchmark does not prevent regeneration of the seven existing figures.

## S4/S6 next stage

Radius calibration and the 12-task coverage pilot completed. All pilot runs passed the declared finite-sample ANN distance-accuracy gate (>=99% of 128 validation queries within 1% of exact distance). This is sampled validation, not a guarantee for every query. Review `tables/S4_S6_pilot_review.tsv` for measured pilot timings. Pilot uses only 50,000 queries, 2% and 10%, and one random seed; it cannot supply the final 20-seed comparison.

The next command is an explicit full-benchmark launch, not a figure-only command. It uses saved raw 60D PeakVI and saved production memberships, with 300 matched-budget random selections (20 seeds × 5 budgets × 3 atlases), plus 15 production references. It queries every atlas cell. No PeakVI, Apricot or DA rerun. Two tasks run concurrently, with existing accuracy checks and fixed atlas-specific radii. Do not rerun radius or pilot jobs.

```bash
BENCHMARK_PHASE=full sbatch --export=ALL,BENCHMARK_PHASE=full \
  --array=0-314%2 publication_figures/run_deferred_coverage_AFTER_APPROVAL.sbatch
```

This writes measured data to `benchmarks/full/`, not figure duplicates. It has not been launched by the consolidation job. Its full-query runtime is not established by the pilot. Existing partial outputs are protected against overwriting; after a failure, inspect the failed task before retrying, rather than resubmitting the entire array.

After all 315 tasks complete successfully:

```bash
sbatch publication_figures/regenerate_all_figures.sbatch
```

That writes S4_measured_60D_coverage, S6_measured_random_comparison and S6b_measured_rare_population_counts into the same `figures/` directory, with the existing figures. S6 covers composition/retention versus random sampling; random-sampling DA comparisons remain a separate, unexecuted analysis and are not implied by S6.
