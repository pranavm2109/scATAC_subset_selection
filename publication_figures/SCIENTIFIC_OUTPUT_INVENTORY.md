# Scientific-output inventory and preservation guide

The complete pre-reorganization file inventory, sizes and SHA-256 hashes are in `archive/pre_reorganization_20261018_235114/inventory_before.tsv`. Its 1,571 publication files were archived before packaging; `BACKUP_VERIFIED.json` records verification of each archived file. The archive also includes the original manuscript-preparation materials. No permanent deletion was performed.

| Material | Exact existing area / identifying files | Preservation decision |
|---|---|---|
| Full benchmark measurements | `benchmarks/full/*_metrics.tsv`, `*_populations.tsv` | Retain all measurements, including all-query and nonselected-query records. |
| Benchmark seed and membership provenance | `benchmarks/full/*_provenance.json` | Retain every seed, membership hash and query/search protocol record. These records do not necessarily contain explicit membership lists for all 20 seeds. |
| Pilot and radius calibration | `benchmarks/pilot/`, `benchmarks/radius/` | Retain separately from full measurements; pilot files are not presumed obsolete duplicates. |
| Explicit random DA memberships | `benchmarks/random_DA_stage_A_2001/memberships/*_selected_cells.txt` | Preserve all six reconstructed seed-2001 memberships at 2%/10% for three atlases. |
| Membership/Seurat validation and DA preflight | `benchmarks/random_DA_stage_A_2001/preflight/membership_manifest.tsv`, `seurat_validation.tsv`, `production_script.sha256`, `R_sessionInfo.txt`, smoke-test inputs/outputs and `PASS` | Preserve all provenance and validation, including small smoke-test artifacts. |
| Saved DA pilot results | `benchmarks/random_DA_stage_A_2001/runs/` | Retain per-atlas/per-budget DA outputs, status and effects. No DA recomputation is part of packaging. |
| DA comparison summaries and resources | `reports/random_DA_stage_A_2001/` | Preserve pooled/per-type/original/common-cohort metrics, rare-population metrics, source definitions, failures and actual resource records. Random DA has one seed, not 20. |
| Scientific audits | `reports/scientific_followup/` | Preserve REPORT.md and all supporting comparison, protocol, validation, geometry and historical-diagnostic tables. |
| Derived plotting inputs / source hashes | `tables/`, especially `SOURCE_MANIFEST.tsv`, `S7_S8_source_manifest.tsv`, identity checks and preservation tables | Retain original table paths. `figure_tables/` exposes existing tables through links and contains the final Figure 1/S9 input copies. |
| Plotting and analysis code | `scripts/`, root `*.sbatch`, `plans/` | Preserve every original script and plan. `figure_scripts/` links original dependencies and adds final-only Figure 1/S9 wrappers; historical analysis launchers are not invoked. |
| Execution and accounting records | `logs/`, DA resource tables and SLURM accounting exports | Preserve successful and failed run records; they are scientific/computational provenance. |
| Accepted supplementary exports/captions | `figures/S*.*`, `captions/S*.md`, caption templates and proofs | Original exports and captions remain unchanged. `supplementary_figures/` contains byte-identical accepted exports and the new S9 overview. |
| Final main Figure 1 | `main_figures/Figure_1_rare_population_representation_and_DA_fidelity.*` | Final-named copies of the reviewed manuscript-preparation exports; exact values and sources in `figure_tables/Figure_1_sources/`. |
| Former main overview | `figures/MAIN_candidate_quantitative.*`, `tables/MAIN_*.tsv`, its original caption and function | Preserve originals; regenerate the same supported four panels as S9, using byte-identical table copies. |
| Historical drafts / duplicate exports | `DELETION_MANIFEST.tsv`; `reports/publication_reorganization_20261018_235114/exact_duplicate_inventory.tsv` | Proposed export cleanup only. Duplicate contents do not justify deleting scientific provenance. No deletions authorized or executed. |

Upstream production objects, PeakVI representations, selected-cell lists and production DA outputs outside `publication_figures/` are referenced by the existing source manifests. They remain in their original project directories and are not relocated, recomputed or represented as included in this directory-only archive.

The current manuscript structure is in `FIGURE_INDEX.md`; earlier inventories/plans/README files are preserved historical records and may describe stages before the benchmarks were completed.
