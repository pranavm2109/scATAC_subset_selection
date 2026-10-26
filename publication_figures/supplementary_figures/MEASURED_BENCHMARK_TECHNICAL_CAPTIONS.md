# Measured benchmark figures — current technical reading guide

These notes accompany the preserved exports, and do not replace or modify their historical shared caption template.

**S4:** Rows show median, 90th-percentile and 99th-percentile nearest-selected distances and the fraction within a fixed radius, using saved 60-dimensional PeakVI coordinates. Columns are Kidney, HECA and Fetal. Distances use calibrated approximate neighbor search; query protocols and accuracy limitations are recorded in the saved benchmark metrics and `reports/scientific_followup/REPORT.md`. All-query distances include selected cells at zero. Calibration does not guarantee exact tail distances.

**S6:** Rows show annotated types retained, total variation distance and Jensen–Shannon divergence in bits. Columns are Kidney, HECA and Fetal. These are composition/representation comparisons, not DA tests.

**S6b:** Six deliberately chosen population examples show counts across budgets: Kidney Mes and Macula_Densa; HECA Haematopoietic stem and progenitor cell and Naive CD4 T cell; Fetal Megakaryocytes? and Ganglion cells. The dotted line marks 10 target cells; DA also requires at least 20 rest cells. Original annotation uncertainty is retained. These examples do not establish universal rare-population enrichment.

For all three: one production selection per budget is compared with 20 saved uniform-random seeds (2001–2020). Shading is the empirical 2.5th–97.5th percentile range across seeds, not a confidence interval. Inputs: `tables/S4_S6_measured_benchmark_summary.tsv`, `tables/S6b_measured_rare_population_counts.tsv`, and `benchmarks/full/*_metrics.tsv` / `*_populations.tsv`; exact protocols, seed/membership hashes and validation records remain in benchmarks. The random DA pilot is separate and has seed 2001 only. Neither S4 nor S6/S6b should be described as a 20-seed DA benchmark.

**Line key:** colored solid lines/points represent the proposed production selection; grey dashed lines show the random median. Grey shading shows the random seed range defined above. This key is particularly necessary for S6b, whose accepted artwork has no inline method legend; its appearance is intentionally unchanged.
