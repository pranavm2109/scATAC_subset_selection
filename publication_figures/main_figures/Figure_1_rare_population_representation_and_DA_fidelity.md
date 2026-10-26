# Figure 1: selective rare-population representation and DA evaluability

All panels use the **2% global sampling budget**: HECA 29,010 cells, Fetal 14,412 cells. The four rows retain original annotation names: Haematopoietic stem and progenitor cell; Megakaryocytes?; Naive CD4 T cell; Ganglion cells. Their full-atlas populations contain 56, 60, 219 and 678 cells, respectively. Question marks indicate supplied annotation uncertainty, not a newly validated identity.

## Panel A — rare-population counts

The colored diamond shows the one production selection. Each grey dot is one saved, matched-budget uniform random selection, seeds 2001–2020. Overlapping integer counts are stacked deterministically by seed order; no memberships are generated and no data points are discarded. The central strip spans the interquartile range and the thin line spans the observed minimum–maximum. These describe seed variability, not confidence intervals. Seed 2001 is outlined to identify the random selection used in panel C. The dotted vertical line is 10 target cells.

Source: `publication_figures/tables/S6b_measured_rare_population_counts.tsv`, filtered to subset_size=02pct, the four specified original annotations, all 20 random seeds and the production selection. Production counts are 14/11/1/1 in the stated row order. The corresponding enrichment ratios are 12.50×, 9.17×, .228× and .0737×, relative to each population's original proportion. This is descriptive selective enrichment/depletion, not preservation of proportions.

## Panel B — DA eligibility

“Yes”/“No” reports production eligibility. The random column reports the number of eligible runs out of 20; the small grey strip encodes that same fraction without an uncertainty interval. Eligibility uses >=10 target AND >=20 rest cells. All rest-cell requirements are satisfied for these examples, so the target-cell threshold determines eligibility. Production eligibility is Yes/Yes/No/No; random eligibility counts are 0/20, 0/20, 1/20 and 17/20. The naive CD4 example shows depletion despite both methods usually being ineligible at this budget; it is not an additional production advantage.

Sources: saved `DA_eligible` values in the S6b table, cross-checked for the two displayed DA runs against `publication_figures/reports/random_DA_stage_A_2001/eligibility.tsv`. No new DA or statistical tests are performed. Fraction of random runs eligible is not a validated probability for all future samples.

## Panel C — strict reference-peak recovery

P = one production selection; R = **random seed 2001 only**. Bars use the existing original-reference, all-direction strict-recovery metrics: full-atlas significant peak–cell-type pairs that are also significant in the subset with matching direction, divided by all full-atlas significant reference pairs for that type. No reference, threshold or denominator is redefined for display.

The measured eligible cases at 2% are production stem/progenitor (70.76%), production Megakaryocytes? (13.44%) and random ganglion (12.27%). All other displayed method/population cases are DA-ineligible, shown as text with no bar. The source count-based strict-recovery statistic is zero for these ineligible cases, but their effects were unavailable; presenting those zeros as measured failure would be misleading. Source statistics are retained separately from display NA in `figure_tables/Figure_1_sources/figure1_panel_values.tsv`. A genuinely eligible run with zero recovery would instead receive a 0% label at the axis; none occurs among these displayed eligible cases.

Source: `publication_figures/reports/random_DA_stage_A_2001/rare_population_metrics.tsv`, cohort=original_reference, direction=all, subset_size=02pct; eligibility from the matching status table. Per-type full-reference denominators are 18,488 / 12,437 / 24,160 / 51,433 in row order. Displayed random seed 2001 has 2 / 1 / 5 / 10 target cells; production has 14 / 11 / 1 / 1. Do not substitute the 20-seed median count for the DA pilot's actual count.

## Technical caption

**Figure 1. Annotation-blind selection selectively enriches exceptionally rare populations and can enable DA testing at a small sampling budget.** All panels use independent selections at 2% of each atlas. (A) Selected-cell counts for four deliberately chosen enriched and depleted annotations: colored diamonds, production selection; grey dots, 20 uniform random seeds; outline, random seed 2001. Strips show interquartile and observed min–max ranges, not confidence intervals. Full-atlas counts appear at left. (B) Eligibility for the unchanged >=10 target / >=20 rest-cell DA rule; random values count eligible runs among 20 seeds. (C) Strict same-direction significant full-reference peak recovery, measured for one production and one seed-2001 random subset per atlas. DA-ineligible cases have no measured bar and are distinct from eligible 0% recovery. Random DA was not evaluated over 20 seeds. Significance retains the production BPCells Wilcoxon implementation, BH within cell type, adjusted P<=.05, |average log2FC|>=.25 and detection >=1% in either group. Examples are descriptive and not an unbiased sample of populations. Parent/subset overlap, cell-type-specific sample sizes, one-versus-rest composition, the original pseudocount formula and annotation uncertainty limit biological and causal interpretation. Overall composition and typical-cell 60D proximity favor random sampling (S4/S6); DA fidelity is not uniformly improved (Stage A tables). BioRender workflow illustration is a separate graphical abstract.

## Export and preservation

Final exports are byte-identical copies of the reviewed 178-mm PDF, SVG and 600-DPI PNG, with the definitive manuscript filename. The original reviewed files and all source values remain preserved. `figure_scripts/plot_figure_1_final.py` reproduces the same plotting logic with updated output names and paths; it was not necessary to change or regenerate any artist for this rename.
