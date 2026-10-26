# Figure QA

All eight figure sets have PDF, SVG and PNG exports. PDF media boxes are 178 mm wide; PNG previews are 450 DPI; SVG text remains editable. Dense grid artists are intentionally rasterized within vector containers. See tables/FIGURE_QA.tsv for physical dimensions and file checks.

Numerical checks: composition/eligibility counts align with production status tables; count-derived DA metrics reproduce the production panel summary; fixed-common-eligibility denominators remain constant within atlas; streamed density-bin totals equal all 15 production available-effect denominators. Existing small source-file hashes remain unchanged. No sampling or ANN baseline was executed.

Layout checks: summary-panel legends occupy dedicated non-data bands at 7 pt. Rendering-time assertions reject legends clipped by the page or overlapping any panel, title, or axis label; axis-label pairs are also checked for collisions. See tables/LEGEND_QA.tsv. Unused legend bands are collapsed.

Visual review: all figure previews inspected at their declared print sizes; population-retention annotation crowding was corrected, unrelated atlas UMAP distances are not connected by lines, and density-grid whitespace reduced. Fonts are specified at print size (typically 7–9 pt; compact grid annotations/ticks 6–6.5 pt). A final author/journal proof is still appropriate. Descriptive UMAP marker saturation is not a quantitative coverage result.

Deferred work: S4/S6 are scripts/plans only. No result placeholders or simulated baseline points appear in the figure collection. Historical embedding identity across regenerations is not established; within-run common embeddings are verified from the analysis implementation.
