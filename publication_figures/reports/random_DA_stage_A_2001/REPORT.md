# Stage A matched-budget random DA pilot — seed 2001 only

## Scope and validation

Exactly six random DA runs: Kidney/HECA/Fetal × 2%/10%, seed 2001. Memberships were reconstructed from the completed benchmark, their recorded SHA-256 hashes verified, and all IDs validated against both Seurat and BPCells columns. The unchanged production script passed an 80-cell, 100-peak real-input smoke test before pilot submission. Pilot array 11137405 used global concurrency one, four CPUs, 128000 MB requested RAM and one hour per task. Preflight artifacts are in `publication_figures/benchmarks/random_DA_stage_A_2001/preflight/`; actual random outputs in `runs/{atlas}/seed2001/{fraction}/`. No production result or publication figure was overwritten.

Production implementation: `scripts/09_1_run_signac_da_all_cell_types.R` unchanged; BPCells Wilcoxon, >=10 target and >=20 rest cells, BH within type, adjusted P<=.05, |average log2FC|>=.25, detection >=.01 in either group, original pseudocount effect formula. Production subsets come from `results/signac_da/{atlas}/{fraction}/`; full references from `results/signac_da/{atlas}/100pct/signac_da_peaks.tsv`. Comparison verifies identical reference keys/directions/effects, threshold definitions and Boolean flags; see `definition_QA.tsv` and `sources.tsv`.

## Recorded computational outcomes

The table records actual elapsed seconds, Slurm batch peak RSS and CPU accounting, not requested memory as measured use. `slurm_accounting.psv` retains raw task/step records. TotalCPU is measured aggregate CPU time; four times wall time is allocated core time and differs from CPU utilization.

| atlas | subset_size | job_id | state | exit_code | CPUs | requested_memory | elapsed_seconds | peak_batch_RSS_GiB | TotalCPU |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | 11137405_0 | COMPLETED | 0:0 | 4 | 131072000K | 77 | 2.0174 | 00:32.025 |
| kidney_atac | 10pct | 11137405_1 | COMPLETED | 0:0 | 4 | 131072000K | 136 | 2.0188 | 00:42.292 |
| heca_atac | 02pct | 11137405_2 | COMPLETED | 0:0 | 4 | 131072000K | 513 | 7.8994 | 01:54.618 |
| heca_atac | 10pct | 11137405_3 | COMPLETED | 0:0 | 4 | 131072000K | 744 | 7.9719 | 02:42.486 |
| fetal_atac | 02pct | 11137405_4 | COMPLETED | 0:0 | 4 | 131072000K | 485 | 4.9043 | 01:08.561 |
| fetal_atac | 10pct | 11137405_5 | COMPLETED | 0:0 | 4 | 131072000K | 535 | 3.0628 | 01:32.199 |

## Original full-reference comparison

Pooled Spearman correlates signed effects over all available full-atlas significant peak–cell-type pairs. Availability denominator is all reference pairs. Direction agreement and switches divide by available finite effects; zeros enter neither numerator. Strict recovery divides by ALL full-reference pairs, including ineligible types. All rates are fractions, not percentages. Pooled estimates weight types with more reference peaks more heavily; do not infer per-type performance from them.

| atlas | subset_size | method | n_eligible_types | n_reference | n_available | spearman | direction_agreement | sign_switch_rate | strict_recovery |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | production | 22 | 1703705 | 1589649 | 0.6584 | 0.7995 | 0.2004 | 0.0892 |
| kidney_atac | 02pct | random | 21 | 1703705 | 1550224 | 0.6517 | 0.7946 | 0.2047 | 0.0976 |
| kidney_atac | 10pct | production | 24 | 1703705 | 1703705 | 0.8600 | 0.9395 | 0.0604 | 0.3866 |
| kidney_atac | 10pct | random | 24 | 1703705 | 1703705 | 0.8698 | 0.9520 | 0.0480 | 0.4162 |
| heca_atac | 02pct | production | 89 | 8118523 | 7359128 | 0.7963 | 0.8802 | 0.1198 | 0.3616 |
| heca_atac | 02pct | random | 89 | 8118523 | 7590718 | 0.7535 | 0.8674 | 0.1326 | 0.3265 |
| heca_atac | 10pct | production | 99 | 8118523 | 8079938 | 0.9309 | 0.9636 | 0.0364 | 0.6816 |
| heca_atac | 10pct | random | 97 | 8118523 | 8036234 | 0.9289 | 0.9835 | 0.0165 | 0.6589 |
| fetal_atac | 02pct | production | 73 | 3655720 | 3474643 | 0.7568 | 0.8593 | 0.1407 | 0.2745 |
| fetal_atac | 02pct | random | 68 | 3655720 | 3303898 | 0.7878 | 0.8827 | 0.1173 | 0.1933 |
| fetal_atac | 10pct | production | 78 | 3655720 | 3655720 | 0.8911 | 0.9354 | 0.0646 | 0.5978 |
| fetal_atac | 10pct | random | 77 | 3655720 | 3643283 | 0.9008 | 0.9548 | 0.0452 | 0.4699 |

## Fixed common-eligibility comparison

The cohort is the within-atlas intersection of types eligible in BOTH methods at BOTH budgets, for this single random seed. It is fixed across the four analyses, not separately redefined per budget or method. `common_eligible_types.tsv` lists membership. This deliberately excludes some of the rare populations whose eligibility differs, so the original-reference and eligibility tables remain essential. Reference denominators match across methods and budgets within this cohort.

| atlas | subset_size | method | n_cohort_types | n_reference | n_available | spearman | direction_agreement | strict_recovery |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | production | 21 | 1550224 | 1550224 | 0.6669 | 0.8019 | 0.0975 |
| kidney_atac | 02pct | random | 21 | 1550224 | 1550224 | 0.6517 | 0.7946 | 0.1073 |
| kidney_atac | 10pct | production | 21 | 1550224 | 1550224 | 0.8726 | 0.9446 | 0.4043 |
| kidney_atac | 10pct | random | 21 | 1550224 | 1550224 | 0.8844 | 0.9558 | 0.4300 |
| heca_atac | 02pct | production | 83 | 7077576 | 7077576 | 0.7871 | 0.8781 | 0.3972 |
| heca_atac | 02pct | random | 83 | 7077576 | 7077576 | 0.7949 | 0.8940 | 0.3715 |
| heca_atac | 10pct | production | 83 | 7077576 | 7077576 | 0.9455 | 0.9629 | 0.7324 |
| heca_atac | 10pct | random | 83 | 7077576 | 7077576 | 0.9452 | 0.9863 | 0.7144 |
| fetal_atac | 02pct | production | 64 | 3136090 | 3136090 | 0.7491 | 0.8495 | 0.2889 |
| fetal_atac | 02pct | random | 64 | 3136090 | 3136090 | 0.8051 | 0.8808 | 0.2223 |
| fetal_atac | 10pct | production | 64 | 3136090 | 3136090 | 0.8967 | 0.9299 | 0.5989 |
| fetal_atac | 10pct | random | 64 | 3136090 | 3136090 | 0.9173 | 0.9543 | 0.4907 |

## Correlations on exactly the same finite reference pairs

Original-reference correlations above may use different available pairs. Here both methods correlate on the identical pairwise finite intersection. Per-type and open/closed versions are in `paired_finite_effect_correlations.tsv`. These correlations do not replace original strict-recovery denominators.

| atlas | subset_size | cohort | method | n_identical_finite_pairs | spearman | median_absolute_error |
| --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | original_reference | production | 1550224 | 0.6669 | 0.6424 |
| kidney_atac | 02pct | original_reference | random | 1550224 | 0.6517 | 0.6596 |
| kidney_atac | 02pct | common_eligible | production | 1550224 | 0.6669 | 0.6424 |
| kidney_atac | 02pct | common_eligible | random | 1550224 | 0.6517 | 0.6596 |
| kidney_atac | 10pct | original_reference | production | 1703705 | 0.8600 | 0.3822 |
| kidney_atac | 10pct | original_reference | random | 1703705 | 0.8698 | 0.3331 |
| kidney_atac | 10pct | common_eligible | production | 1550224 | 0.8726 | 0.3610 |
| kidney_atac | 10pct | common_eligible | random | 1550224 | 0.8844 | 0.3091 |
| heca_atac | 02pct | original_reference | production | 7077576 | 0.7871 | 1.1022 |
| heca_atac | 02pct | original_reference | random | 7077576 | 0.7949 | 0.9025 |
| heca_atac | 02pct | common_eligible | production | 7077576 | 0.7871 | 1.1022 |
| heca_atac | 02pct | common_eligible | random | 7077576 | 0.7949 | 0.9025 |
| heca_atac | 10pct | original_reference | production | 8012074 | 0.9295 | 0.7034 |
| heca_atac | 10pct | original_reference | random | 8012074 | 0.9300 | 0.4428 |
| heca_atac | 10pct | common_eligible | production | 7077576 | 0.9455 | 0.6487 |
| heca_atac | 10pct | common_eligible | random | 7077576 | 0.9452 | 0.3796 |
| fetal_atac | 02pct | original_reference | production | 3136090 | 0.7491 | 0.6539 |
| fetal_atac | 02pct | original_reference | random | 3136090 | 0.8051 | 0.6051 |
| fetal_atac | 02pct | common_eligible | production | 3136090 | 0.7491 | 0.6539 |
| fetal_atac | 02pct | common_eligible | random | 3136090 | 0.8051 | 0.6051 |
| fetal_atac | 10pct | original_reference | production | 3643283 | 0.8906 | 0.3835 |
| fetal_atac | 10pct | original_reference | random | 3643283 | 0.9008 | 0.3229 |
| fetal_atac | 10pct | common_eligible | production | 3136090 | 0.8967 | 0.3791 |
| fetal_atac | 10pct | common_eligible | random | 3136090 | 0.9173 | 0.2853 |

## Open and closed reference peaks

Direction is fixed by the full-atlas reference effect. Correlations within open/closed groups restrict the effect range and answer a different question from pooled signed correlations. Per-type counterparts are in `per_type_metrics.tsv`.

| atlas | subset_size | method | direction | n_reference | n_available | spearman | direction_agreement | sign_switch_rate | strict_recovery |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | production | open | 589853 | 502544 | 0.6142 | 0.8760 | 0.1239 | 0.2031 |
| kidney_atac | 02pct | production | closed | 1113852 | 1087105 | 0.3236 | 0.7641 | 0.2358 | 0.0288 |
| kidney_atac | 02pct | random | open | 589853 | 493913 | 0.6217 | 0.8808 | 0.1189 | 0.2130 |
| kidney_atac | 02pct | random | closed | 1113852 | 1056311 | 0.2962 | 0.7543 | 0.2448 | 0.0365 |
| kidney_atac | 10pct | production | open | 589853 | 589853 | 0.8095 | 0.9542 | 0.0458 | 0.5543 |
| kidney_atac | 10pct | production | closed | 1113852 | 1113852 | 0.5828 | 0.9318 | 0.0682 | 0.2978 |
| kidney_atac | 10pct | random | open | 589853 | 589853 | 0.8301 | 0.9592 | 0.0408 | 0.5853 |
| kidney_atac | 10pct | random | closed | 1113852 | 1113852 | 0.5960 | 0.9481 | 0.0519 | 0.3266 |
| heca_atac | 02pct | production | open | 2705101 | 2457871 | 0.8021 | 0.8841 | 0.1159 | 0.5396 |
| heca_atac | 02pct | production | closed | 5413422 | 4901257 | 0.4785 | 0.8783 | 0.1217 | 0.2726 |
| heca_atac | 02pct | random | open | 2705101 | 2320101 | 0.8513 | 0.9520 | 0.0480 | 0.5358 |
| heca_atac | 02pct | random | closed | 5413422 | 5270617 | 0.4013 | 0.8301 | 0.1699 | 0.2219 |
| heca_atac | 10pct | production | open | 2705101 | 2688567 | 0.8879 | 0.9272 | 0.0728 | 0.7724 |
| heca_atac | 10pct | production | closed | 5413422 | 5391371 | 0.7903 | 0.9817 | 0.0183 | 0.6363 |
| heca_atac | 10pct | random | open | 2705101 | 2622884 | 0.9434 | 0.9818 | 0.0182 | 0.8092 |
| heca_atac | 10pct | random | closed | 5413422 | 5413350 | 0.7773 | 0.9843 | 0.0157 | 0.5838 |
| fetal_atac | 02pct | production | open | 2618809 | 2440295 | 0.5341 | 0.8679 | 0.1321 | 0.3136 |
| fetal_atac | 02pct | production | closed | 1036911 | 1034348 | 0.4440 | 0.8392 | 0.1608 | 0.1758 |
| fetal_atac | 02pct | random | open | 2618809 | 2269329 | 0.5883 | 0.9211 | 0.0789 | 0.2158 |
| fetal_atac | 02pct | random | closed | 1036911 | 1034569 | 0.3449 | 0.7984 | 0.2016 | 0.1365 |
| fetal_atac | 10pct | production | open | 2618809 | 2618809 | 0.7452 | 0.9250 | 0.0750 | 0.6155 |
| fetal_atac | 10pct | production | closed | 1036911 | 1036911 | 0.7497 | 0.9617 | 0.0383 | 0.5532 |
| fetal_atac | 10pct | random | open | 2618809 | 2606372 | 0.7659 | 0.9485 | 0.0515 | 0.4800 |
| fetal_atac | 10pct | random | closed | 1036911 | 1036911 | 0.7377 | 0.9706 | 0.0294 | 0.4445 |

## Per-cell-type summaries

`per_type_metrics.tsv` contains every reference-bearing type, original/common cohorts and open/closed/all directions. `eligibility.tsv` retains ALL annotated types, including types with no reference peaks. Macro summaries below use equally weighted type medians; correlations with fewer than three finite pairs or zero variance are undefined. Undefined correlations are excluded from the median and their finite count is explicit. Types with no full-reference peaks have undefined recovery, not zero. These distributions are not biological replicates or confidence intervals.

| atlas | subset_size | method | cohort | direction | n_types | n_finite_spearman | median_per_type_spearman | median_per_type_strict_recovery | median_per_type_direction_agreement |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| fetal_atac | 02pct | production | common_eligible | all | 64 | 64 | 0.6549 | 0.1603 | 0.8517 |
| fetal_atac | 02pct | production | original_reference | all | 78 | 73 | 0.6480 | 0.1361 | 0.8611 |
| fetal_atac | 02pct | random | common_eligible | all | 64 | 64 | 0.6061 | 0.0866 | 0.8848 |
| fetal_atac | 02pct | random | original_reference | all | 78 | 68 | 0.6061 | 0.0593 | 0.8860 |
| fetal_atac | 10pct | production | common_eligible | all | 64 | 64 | 0.8850 | 0.5110 | 0.9683 |
| fetal_atac | 10pct | production | original_reference | all | 78 | 78 | 0.8720 | 0.5110 | 0.9666 |
| fetal_atac | 10pct | random | common_eligible | all | 64 | 64 | 0.8149 | 0.3605 | 0.9507 |
| fetal_atac | 10pct | random | original_reference | all | 78 | 77 | 0.7766 | 0.2847 | 0.9545 |
| heca_atac | 02pct | production | common_eligible | all | 83 | 83 | 0.8122 | 0.3687 | 0.9274 |
| heca_atac | 02pct | production | original_reference | all | 101 | 89 | 0.8183 | 0.3107 | 0.9275 |
| heca_atac | 02pct | random | common_eligible | all | 83 | 83 | 0.8202 | 0.3110 | 0.9562 |
| heca_atac | 02pct | random | original_reference | all | 101 | 89 | 0.7955 | 0.2031 | 0.9485 |
| heca_atac | 10pct | production | common_eligible | all | 83 | 83 | 0.9404 | 0.8433 | 0.9797 |
| heca_atac | 10pct | production | original_reference | all | 101 | 99 | 0.9399 | 0.8140 | 0.9783 |
| heca_atac | 10pct | random | common_eligible | all | 83 | 83 | 0.9528 | 0.8414 | 0.9936 |
| heca_atac | 10pct | random | original_reference | all | 101 | 97 | 0.9265 | 0.7576 | 0.9923 |
| kidney_atac | 02pct | production | common_eligible | all | 21 | 21 | 0.5796 | 0.0342 | 0.8039 |
| kidney_atac | 02pct | production | original_reference | all | 24 | 22 | 0.5640 | 0.0289 | 0.8002 |
| kidney_atac | 02pct | random | common_eligible | all | 21 | 21 | 0.6295 | 0.0395 | 0.8360 |
| kidney_atac | 02pct | random | original_reference | all | 24 | 21 | 0.6295 | 0.0269 | 0.8360 |
| kidney_atac | 10pct | production | common_eligible | all | 21 | 21 | 0.8384 | 0.3703 | 0.9453 |
| kidney_atac | 10pct | production | original_reference | all | 24 | 24 | 0.8254 | 0.3135 | 0.9443 |
| kidney_atac | 10pct | random | common_eligible | all | 21 | 21 | 0.8482 | 0.4270 | 0.9620 |
| kidney_atac | 10pct | random | original_reference | all | 24 | 24 | 0.8307 | 0.4075 | 0.9583 |

## Enriched AND depleted focal populations

Target/rest counts and eligibility:

| atlas | subset_size | method | cell_type | n_target_cells | n_rest_cells | status |
| --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | production | Macula_Densa | 3 | 1614 | insufficient_target_cells |
| kidney_atac | 02pct | production | Mes | 4 | 1613 | insufficient_target_cells |
| kidney_atac | 02pct | random | Macula_Densa | 1 | 1616 | insufficient_target_cells |
| kidney_atac | 02pct | random | Mes | 1 | 1616 | insufficient_target_cells |
| kidney_atac | 10pct | production | Macula_Densa | 10 | 8074 | success |
| kidney_atac | 10pct | production | Mes | 17 | 8067 | success |
| kidney_atac | 10pct | random | Macula_Densa | 15 | 8069 | success |
| kidney_atac | 10pct | random | Mes | 16 | 8068 | success |
| heca_atac | 02pct | production | Haematopoietic stem and progenitor cell | 14 | 28996 | success |
| heca_atac | 02pct | production | Naive CD4 T cell | 1 | 29009 | insufficient_target_cells |
| heca_atac | 02pct | random | Haematopoietic stem and progenitor cell | 2 | 29008 | insufficient_target_cells |
| heca_atac | 02pct | random | Naive CD4 T cell | 5 | 29005 | insufficient_target_cells |
| heca_atac | 10pct | production | Haematopoietic stem and progenitor cell | 47 | 145004 | success |
| heca_atac | 10pct | production | Naive CD4 T cell | 4 | 145047 | insufficient_target_cells |
| heca_atac | 10pct | random | Haematopoietic stem and progenitor cell | 4 | 145047 | insufficient_target_cells |
| heca_atac | 10pct | random | Naive CD4 T cell | 21 | 145030 | success |
| fetal_atac | 02pct | production | Ganglion cells | 1 | 14411 | insufficient_target_cells |
| fetal_atac | 02pct | production | Megakaryocytes? | 11 | 14401 | success |
| fetal_atac | 02pct | random | Ganglion cells | 10 | 14402 | success |
| fetal_atac | 02pct | random | Megakaryocytes? | 1 | 14411 | insufficient_target_cells |
| fetal_atac | 10pct | production | Ganglion cells | 28 | 72033 | success |
| fetal_atac | 10pct | production | Megakaryocytes? | 20 | 72041 | success |
| fetal_atac | 10pct | random | Ganglion cells | 84 | 71977 | success |
| fetal_atac | 10pct | random | Megakaryocytes? | 4 | 72057 | insufficient_target_cells |

Full-reference metrics (NA correlation means insufficient or invariant available effects; missing availability must not be interpreted as concordance):

| atlas | subset_size | method | cell_type | n_reference | n_available | spearman | direction_agreement | sign_switch_rate | strict_recovery |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_atac | 02pct | production | Macula_Densa | 71750 | 0 | NA | NA | NA | 0.0000 |
| kidney_atac | 02pct | production | Mes | 42306 | 0 | NA | NA | NA | 0.0000 |
| kidney_atac | 02pct | random | Macula_Densa | 71750 | 0 | NA | NA | NA | 0.0000 |
| kidney_atac | 02pct | random | Mes | 42306 | 0 | NA | NA | NA | 0.0000 |
| kidney_atac | 10pct | production | Macula_Densa | 71750 | 71750 | 0.7743 | 0.8246 | 0.1754 | 0.3208 |
| kidney_atac | 10pct | production | Mes | 42306 | 42306 | 0.8236 | 0.9270 | 0.0730 | 0.1637 |
| kidney_atac | 10pct | random | Macula_Densa | 71750 | 71750 | 0.8280 | 0.8910 | 0.1090 | 0.3367 |
| kidney_atac | 10pct | random | Mes | 42306 | 42306 | 0.8335 | 0.9749 | 0.0251 | 0.4182 |
| heca_atac | 02pct | production | Haematopoietic stem and progenitor cell | 18488 | 18488 | 0.7774 | 0.9916 | 0.0084 | 0.7076 |
| heca_atac | 02pct | production | Naive CD4 T cell | 24160 | 0 | NA | NA | NA | 0.0000 |
| heca_atac | 02pct | random | Haematopoietic stem and progenitor cell | 18488 | 0 | NA | NA | NA | 0.0000 |
| heca_atac | 02pct | random | Naive CD4 T cell | 24160 | 0 | NA | NA | NA | 0.0000 |
| heca_atac | 10pct | production | Haematopoietic stem and progenitor cell | 18488 | 18488 | 0.9609 | 0.9941 | 0.0059 | 0.7226 |
| heca_atac | 10pct | production | Naive CD4 T cell | 24160 | 0 | NA | NA | NA | 0.0000 |
| heca_atac | 10pct | random | Haematopoietic stem and progenitor cell | 18488 | 0 | NA | NA | NA | 0.0000 |
| heca_atac | 10pct | random | Naive CD4 T cell | 24160 | 24160 | 0.8032 | 0.9922 | 0.0078 | 0.0230 |
| fetal_atac | 02pct | production | Megakaryocytes? | 12437 | 12437 | 0.6456 | 0.9675 | 0.0325 | 0.1344 |
| fetal_atac | 02pct | production | Ganglion cells | 51433 | 0 | NA | NA | NA | 0.0000 |
| fetal_atac | 02pct | random | Megakaryocytes? | 12437 | 0 | NA | NA | NA | 0.0000 |
| fetal_atac | 02pct | random | Ganglion cells | 51433 | 51433 | 0.5518 | 0.9292 | 0.0708 | 0.1227 |
| fetal_atac | 10pct | production | Megakaryocytes? | 12437 | 12437 | 0.8430 | 0.9993 | 0.0007 | 0.6888 |
| fetal_atac | 10pct | production | Ganglion cells | 51433 | 51433 | 0.6067 | 0.9506 | 0.0494 | 0.3709 |
| fetal_atac | 10pct | random | Megakaryocytes? | 12437 | 0 | NA | NA | NA | 0.0000 |
| fetal_atac | 10pct | random | Ganglion cells | 51433 | 51433 | 0.6729 | 0.9362 | 0.0638 | 0.4093 |

## Interpretation limits and stopping rule

This is ONE random seed per atlas/budget, not a random-seed distribution. Do not infer statistical superiority, uncertainty bands or robustness across seeds. These results evaluate reproduction of a parent-atlas estimate, not external biological truth; parent/subset overlap affects concordance. Differences may involve eligibility, effect missingness, sample size, the group-size-dependent pseudocount and changes in one-versus-rest composition. The pilot cannot identify a causal biological or power-only explanation. An eligibility gain is a real evaluability gain; it is not by itself validated biological discovery. Assess pooled, macro/per-type and common-pair results together, and retain both enriched and depleted examples.

The original effect formula is preserved even for peaks undetected in both subset groups: with zero counts, its +1 count pseudocount can yield log2(n_rest/n_target), a nonzero effect driven by group sizes. These pairs fail the detection requirement for strict recovery but can enter available-effect direction/correlation metrics under the original definitions. The frequency of this scenario was not separately tabulated in this completed pilot; its impact remains unquantified. This is a definition qualification, not a new biological effect or a reason to change the approved analysis silently.

No additional seeds, Geosketch, coverage sensitivity jobs or publication figures were launched. Stop after reporting this pilot. Any further experiment requires explicit user approval.

## Descriptive cross-case summary

Counts below summarize the six atlas/budget point comparisons for this one random seed. Higher is favorable for these three metrics. Counts are not significance tests, independent biological replications or evidence of robustness across random seeds. Original-reference correlations can use different available pairs; consult the identical-finite-pair table as well.

| cohort | n_atlas_budget_cases | spearman_production_higher | spearman_random_higher | strict_recovery_production_higher | strict_recovery_random_higher | direction_agreement_production_higher | direction_agreement_random_higher |
| --- | --- | --- | --- | --- | --- | --- | --- |
| original_reference | 6 | 3 | 3 | 4 | 2 | 2 | 4 |
| common_eligible | 6 | 2 | 4 | 4 | 2 | 1 | 5 |

## Scientific interpretation of the completed pilot

1. **Strict recovery improves in HECA and Fetal, not Kidney.** Original-reference recovery (proposed/random) is Kidney 2% 8.92/9.76%, Kidney 10% 38.66/41.62%, HECA 2% 36.16/32.65%, HECA 10% 68.16/65.89%, Fetal 2% 27.45/19.33%, and Fetal 10% 59.78/46.99%. These are single-seed point comparisons, not method-level significance tests. Source: `pooled_metrics.tsv`, original_reference/all.

2. **The HECA/Fetal recovery gains are not solely an eligibility-denominator artifact.** With exactly the same fixed common-eligible reference pairs and complete effect availability in both methods, proposed strict recovery remains higher in both HECA and Fetal budgets. Fetal gains are 28.89 vs 22.23% at 2%, and 59.89 vs 49.07% at 10%. This rules out eligibility-related missing effects as the sole explanation for those cohort-specific differences, but does not isolate cell-number power, selected-cell characteristics, one-versus-rest composition or pseudocount effects. Source: `pooled_metrics.tsv`, common_eligible/all.

3. **Effect fidelity is not uniformly improved.** Within fixed common cohorts, random has higher pooled signed Spearman in four of six cases and higher pooled direction agreement in five of six. In HECA 2%, the apparent original pooled correlation advantage reverses on identical finite pairs (proposed .7871, random .7949). For Fetal 10%, the equal-type median signed correlation favors proposed (.8850 vs .8149), whereas pooled correlation favors random (.8967 vs .9173). Pooling weights peak-rich types and includes between-type structure; equal-type medians answer another question. Neither should substitute for the other. Sources: `paired_finite_effect_correlations.tsv`, `macro_type_summaries.tsv`, `pooled_metrics.tsv`.

4. **Direction-specific tradeoffs matter.** HECA 10% proposed strict recovery is higher for closed reference peaks (63.63 vs 58.38%) but lower for open peaks (77.24 vs 80.92%). Fetal strict-recovery gains occur for both open and closed reference peaks at both budgets. Report these alongside the signed pooled totals. Source: `pooled_metrics.tsv`, original_reference/open and closed.

5. **There are concrete rare-type gains and counterexamples.** Proposed HECA stem/progenitor strict recovery is 70.76% at 2% and 72.26% at 10%; random is ineligible at both budgets. Proposed Megakaryocytes? recovery is 13.44% and 68.88%; random is again ineligible. Conversely, random alone is eligible for naive CD4 T at 10% and ganglion at 2%. At 10%, Mes has almost equal target counts (17 proposed, 16 random) but markedly different strict recovery (16.37 vs 41.82%); target-cell count alone cannot explain every difference. This does not establish a biological mechanism or an external truth advantage. Sources: `eligibility.tsv`, `rare_population_metrics.tsv`.

**Conclusion:** this pilot supplies evidence of selective DA-recovery advantages that the global nearest-distance benchmark could not establish, alongside real effect-fidelity and depleted-population disadvantages. It supports a tradeoff-based scientific interpretation, not a claim of universal biological superiority. Robustness across random seeds remains unknown. Stop here as approved; no additional experiments have been launched.
