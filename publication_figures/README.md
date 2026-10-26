# Quantitative publication figure collection

**Current entrypoint:** [REGENERATE.md](REGENERATE.md). All current supplementary figures are in `figures/`, regenerated with `sbatch publication_figures/regenerate_all_figures.sbatch` as PDF, SVG and 600-dpi PNG. Duplicate polishing directories have been consolidated. S4/S6 require the complete measured benchmark described in that guide; plotting does not launch it.

New, question-led figures for the scATAC representative-subset project. Original production results and original 3x5 figures are unchanged. All generated composites target **178 mm double-column width**. Current supplementary exports are PDF, SVG, and 600-DPI PNG. Dense grid points/densities are rasterized inside otherwise vector PDF/SVG containers; editable text is retained.

Start with [the main candidate](figures/MAIN_candidate_quantitative.pdf), [S1](figures/S1_normalization.pdf), and [S3](figures/S3_DA_fidelity.pdf). Use the BioRender workflow separately if journal space allows; do not overload the quantitative main composite with the full pipeline diagram.

## Inventory and figure guide

The initial [figure inventory](FIGURE_INVENTORY.md) separates existing-output figures from unexecuted new analyses. [QA_REPORT.md](QA_REPORT.md) and [FIGURE_QA.tsv](tables/FIGURE_QA.tsv) record final checks.

| Figure / question and intended finding | Outputs: PDF / SVG / PNG | Script | Figure-specific data and caption |
|---|---|---|---|
| **S1. Does L1 normalization improve selection?** Historical matched-budget retention/composition comparison; distance evidence explicitly restricted to historical 2D UMAP. | [PDF](figures/S1_normalization.pdf) · [SVG](figures/S1_normalization.svg) · [PNG](figures/S1_normalization.png) | [plot_existing_summaries.py](scripts/plot_existing_summaries.py), `s1()` | [2% methods](tables/S1_normalization_2pct.tsv), [budget series](tables/S1_normalization_budgets.tsv), [identity checks](tables/S1_historical_subset_identity.tsv), [caption](captions/S1_normalization.md) |
| **S2. Presence versus proportional representation.** Retention, TV/JS, all-population enrichment/abundance, deliberately selected rare-cell examples and DA ineligibility. | [PDF](figures/S2_population_representation.pdf) · [SVG](figures/S2_population_representation.svg) · [PNG](figures/S2_population_representation.png) | same script, `s2()` | [composition](tables/S2_composition.tsv), [all populations](tables/S2_all_population_metrics.tsv), [rare examples](tables/S2_rare_examples.tsv), [caption](captions/S2_population_representation.md) |
| **S3. DA fidelity and eligibility.** Pooled effect/significance curves, availability, exact common-eligibility count metrics, per-type open/closed correlation distributions. | [PDF](figures/S3_DA_fidelity.pdf) · [SVG](figures/S3_DA_fidelity.svg) · [PNG](figures/S3_DA_fidelity.png) | same script, `s3()` | [pooled](tables/S3_pooled.tsv), [per-type](tables/S3_per_cell_type_direction.tsv), [fixed cohort](tables/S3_fixed_eligibility.tsv), [cohort members](tables/S3_common_eligible_types.tsv), [caption](captions/S3_DA_fidelity.md) |
| **S3b. Direction-specific significant recovery.** Open/closed peak recovery for all references and fixed-common-eligible types. | [PDF](figures/S3b_direction_specific_recovery.pdf) · [SVG](figures/S3b_direction_specific_recovery.svg) · [PNG](figures/S3b_direction_specific_recovery.png) | same script, `s3()` | [direction summaries](tables/S3b_direction_specific_recovery.tsv), [caption](captions/S3b_direction_specific_recovery.md) |
| **S4. Original 60D coverage.** No measured result generated; plan and script only. | **Deferred** | [deferred_benchmarks.py](scripts/deferred_benchmarks.py) | [approval/resource plan](plans/S4_S6_APPROVAL_PLAN.md) |
| **S5. Computational costs.** Selection fit/runtime memory separated from conversion, training and features. | [PDF](figures/S5_computational_cost.pdf) · [SVG](figures/S5_computational_cost.svg) · [PNG](figures/S5_computational_cost.png) | same summary script, `s5()` | [selection](tables/S5_selection.tsv), [stages](tables/S5_stage_costs.tsv), [caption](captions/S5_computational_cost.md) |
| **S6. Random sampling distribution.** No random selection or DA baseline executed. | **Deferred** | [benchmark script](scripts/deferred_benchmarks.py), [future plotter](scripts/plot_random_benchmark_when_available.py), [DA wrapper](scripts/deferred_random_DA.py) | [315 metadata-only tasks](tables/S4_S6_DEFERRED_task_manifest.tsv), [observed DA cost references](tables/S6_DA_cost_reference.tsv), [approval plan](plans/S4_S6_APPROVAL_PLAN.md) |
| **S7. Revised descriptive UMAP grid.** Legible, shared atlas coordinates, no quantitative coverage claim. | [PDF](figures/S7_descriptive_UMAP_grid.pdf) · [SVG](figures/S7_descriptive_UMAP_grid.svg) · [PNG](figures/S7_descriptive_UMAP_grid.png) | [revised_grids.py](scripts/revised_grids.py), `umap_grid()` | [counts](tables/S7_descriptive_UMAP_grid.tsv), [exact sources](tables/S7_S8_source_manifest.tsv), [caption](captions/S7_descriptive_UMAP_grid.md) |
| **S8. Revised DA effect-density grid.** Signed axes, shared density legend, minimal repeated labels; all other metrics in S3. | [PDF](figures/S8_DA_effect_density_grid.pdf) · [SVG](figures/S8_DA_effect_density_grid.svg) · [PNG](figures/S8_DA_effect_density_grid.png) | same grid script, `da_grid()` | [density bins](tables/S8_density_bins.tsv), [panel metrics](tables/S8_panel_metrics.tsv), [exact sources](tables/S7_S8_source_manifest.tsv), [caption](captions/S8_DA_effect_density_grid.md) |
| **Main candidate.** Four-panel methodological/biological story: normalization retention, composition TV, DA direction agreement, strict recovery. | [PDF](figures/MAIN_candidate_quantitative.pdf) · [SVG](figures/MAIN_candidate_quantitative.svg) · [PNG](figures/MAIN_candidate_quantitative.png) | summary script, `main_candidate()` | [normalization](tables/MAIN_normalization.tsv), [composition](tables/MAIN_composition.tsv), [DA](tables/MAIN_DA.tsv), [caption](captions/MAIN_candidate_quantitative.md) |

## Legend/layout revision

All summary-panel legends now occupy dedicated bands outside the data axes, with 7-pt text and clear dashed-line keys. Rendering rejects legend clipping, legend/panel overlap, and axis-label collisions. See [LEGEND_QA.tsv](tables/LEGEND_QA.tsv). To regenerate this revision without rereading large grids:

```bash
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
sbatch publication_figures/rework_legends.sbatch
```

## Exact original sources and scientific limitations

[Source manifest](tables/SOURCE_MANIFEST.tsv) enumerates the original small tables with hashes and figure associations. Grid source paths are enumerated in [S7/S8 manifest](tables/S7_S8_source_manifest.tsv). Every caption identifies its original inputs and derivative TSVs. Captions are **technical working notes for author review**, not manuscript-ready prose.

- [Historical provenance](plans/HISTORICAL_PROVENANCE.md): common historical UMAP geometry is verified within each ablation by implementation, not by archived embedding hashes across regenerated runs. Do not call these 60D distances. No preserved unnormalized membership file was verified; summaries are still valid historical observations.
- Fixed-common-eligible direction/recovery metrics are reconstructed EXACTLY from per-type counts. Pooled common-cohort Spearman cannot be reconstructed from per-type correlations; no such value is invented.
- S3 correlation boxplots represent population heterogeneity, not sampling confidence intervals. Some directions have very few reference peaks. Nonfinite correlations are excluded from boxes; the full table preserves all rows.
- All selected budgets remain independent selections. Every population being present does not establish preserved proportions or sufficient counts for DA.
- S8 reads existing peak effects and aggregates density only; no new DA tests. Square-bin density is explicitly distinguished from the original hexbin rendering. Display clipping is documented.
- No random baseline advantage, 60D coverage result, repeated-timing uncertainty, or independent-donor validation is claimed.

## Reproduce only the authorized existing-output figures

```bash
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
sbatch publication_figures/run_existing_figures.sbatch
```

This writes only `publication_figures/` derivatives. It does not invoke the production Snakemake pipeline, train PeakVI, compute subsets, or perform DA testing. It reads small summaries plus existing UMAP/detail tables for grids; run on compute nodes, not login nodes. Original figures remain in `figures/publication/`.

The resource/command plan for new S4/S6 benchmarks is separate and **awaits approval**. It proposes 20 random seeds, staged NN accuracy/resource calibration, and a smaller staged DA comparison with measured-production cost references. Do not submit those jobs as part of figure regeneration.
