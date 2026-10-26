# Figure inventory — existing-output work first

All figures are new derivatives. Production scripts/results and original grids are preserved.

| Figure | Question / scope | Existing sources | Status / new analysis needed |
|---|---|---|---|
| S1 | RBF normalization: retention, composition, historical UMAP distances | diagnostics/rbf_normalization_ablation/*/comparison/method_summary.tsv; diagnostics/rbf_normalization_generalization/*/fraction_method_summary.tsv | Generate now. Historical within-run methods share a single UMAP input by script design; no claim of 60D coverage or identity with regenerated production UMAPs. |
| S2 | Presence versus proportional composition, rare populations | results/cell_type_proportions/peakvi/*/cell_type_proportions_peakvi.tsv; results/signac_da/*/*/signac_da_status.tsv | Generate now; TV/JS calculated from existing small count tables, no new cell selection. |
| S3 + S3b | DA effect fidelity, strict recovery, eligibility, open/closed heterogeneity | results/publication_panels/signac_da_recapitulation_3x5_summary.tsv; results/signac_da_capture/*/da_peak_capture_summary.tsv | Generate now. Fixed-common-eligibility count metrics exact from summaries; pooled Spearman under fixed eligibility not reconstructible from per-type correlations. |
| S4 | Quantitative original 60D latent coverage | saved latent.npy, cell_names.txt, selected IDs | Scripts/plan only. New nearest-neighbor processing needs approval. No S4 result figure until measurements exist. |
| S5 | Selection and upstream computational costs | production selection summaries and logs/benchmarking | Generate now. Individual successful runs; CPU RSS, not GPU memory. |
| S6 | Matched-budget random reference distributions | saved cell IDs/annotations and latent matrices | Scripts/plan only; at least 20 seeds. No random benchmark or DA jobs launched. |
| S7 | Revised descriptive 3x5 UMAP grid | saved production UMAP coordinates and selected IDs | Generate on compute node; no retraining or coverage inference. |
| S8 | Revised 3x5 DA effect-density grid | existing DA detail TSVs and publication metrics | Generate on compute node by streaming existing records; no new tests; histogram density rather than hexbin. |
| Main candidate | Normalization, composition, DA direction agreement and strict recovery | small tables used by S1–S3 | Generate now, four panels, no unsupported random-baseline or 60D coverage claim. |

## Design
All composites target 178 mm double-column width. PDF and SVG exports contain editable text and vector axes; dense scatter/density artists may be embedded raster layers. PNG previews are 450 DPI. Fonts target 7–9 pt at actual print size; captions are separate technical working notes, not manuscript prose.

## Scope safeguards
No PeakVI fitting, Apricot selection, production DA, baseline generation, or latent coverage jobs are part of existing-output figure generation. S4/S6 scripts require an explicit --execute argument. Historical UMAP summary provenance is weaker than an archived embedding hash: scripts establish within-run comparability, not bytewise equivalence to current embeddings. All current derived metrics retain original populations and independent subset-size memberships.
