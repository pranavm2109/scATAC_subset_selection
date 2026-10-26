# Manuscript-ready publication materials

Use [FIGURE_INDEX.md](FIGURE_INDEX.md) for the current manuscript structure. Existing README/REGENERATE files are historical records; their older main-candidate/deferred-benchmark statements are superseded by this index. No historical file was overwritten.

```
publication_figures/
  main_figures/             # definitive Figure 1, three formats and technical caption
  supplementary_figures/    # ten accepted figures unchanged + S9 overview; captions
  figure_scripts/           # final-only plotting scripts + links to existing scripts/batch files
  figure_tables/            # final source tables + links to existing derived/QA tables
  benchmarks/               # full, pilot, radius, random_DA_stage_A_2001 (untouched)
  reports/                  # existing scientific audits + publication_reorganization_20261018_235114/
  archive/pre_reorganization_20261018_235114/  # verified complete pre-change tar + inventory
  figures/, scripts/, tables/, captions/, logs/, plans/  # compatibility paths, unchanged
  FIGURE_INDEX.md / .tsv
  DELETION_MANIFEST.tsv      # proposals only; nothing deleted
  PUBLICATION_README.md
```

The graphical abstract remains the separate BioRender file outside this directory. No experiments, PeakVI, Apricot or DA were run. To reproduce only S9 on a compute node, use the final plotting script through the documented plotting batch entrypoint; do not use historical regenerate-all scripts for this finalization.

Archive and inventory: `archive/pre_reorganization_20261018_235114`. QC and PDF proofs: `reports/publication_reorganization_20261018_235114`. Supplementary sequential numbering is proposed in the index, not applied to preserved exports. Main Figure 1 exports are exact reviewed copies with the final filename. S9 retains the original overview panels, metrics and technical qualifications.

## Reproduce Figure 1 and S9 only (optional; no submission needed now)

```bash
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
sbatch publication_figures/figure_scripts/render_final_figures.sbatch
```

This plotting-only entrypoint reads saved tables and writes the final named Figure 1/S9 exports. It does not regenerate S1–S8 or run analyses. The archival/reorganization driver is a one-time operation, not a figure-regeneration command. Nothing needs rerunning for the delivered exports.

See `SCIENTIFIC_OUTPUT_INVENTORY.md` for the scientific-data preservation guide and the latest report's `VISUAL_REVIEW.md` for the individual PDF review.

## All final manuscript figures (current complete entrypoint)

```bash
sbatch publication_figures/figure_scripts/render_all_publication_figures.sbatch
```

This renders Figure 1 plus S1–S9, including S3b/S6b, using saved inputs, and checks all 12 export sets. It stages legacy plotting functions so original `figures/`, `tables/` and captions remain unchanged. Final-named exports in `main_figures/` and `supplementary_figures/` are intentionally refreshed when this command is run. Reviewed technical captions remain in place. It does not launch benchmarks or DA. It requires completed, compatible historical diagnostic/benchmark tables; a fresh Snakemake run alone does not supply these. See the updated project `pseudocode.txt` for the complete workflow and run commands.
