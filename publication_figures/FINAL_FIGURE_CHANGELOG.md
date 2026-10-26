# Final minimal corrections

Only S1 and S3 are revised in `figures/`. Their accepted predecessors remain in `../polished/figures/`.

- S1 C: requested title, “Nearest-selected-cell distances in 2D UMAP (2%),” wrapped over two lines to fit the existing panel. All other panels retain their data and styling. Caption explicitly separates historical 2D distances from 60-dimensional PeakVI coverage, records the repeated historical membership checks, and preserves the coordinate-provenance limitation.
- S3 D: reference-pair availability uses circles; cell-type availability uses squares. Both use solid lines, with an explicit marker legend. E/F retain solid original/all-reference and dashed fixed-common-cohort curves. No nine-panel layout, atlas colors, numerical values, axis limits or statistical definitions changed.
- Updated S1/S3 exports: PDF with embedded TrueType font programs, SVG with vector text, and 600-dpi PNG. Source retains the existing font hierarchy and 178-mm figure width.
- S7 already explicitly states the 3,000-selected-cell display cap, seed, background cap, actual n and displayed shown counts. It is unchanged. Kidney 2% has fewer than 3,000 selected cells, so its displayed and actual counts correctly coincide.
- S2, S3b, S5 and S8 remain unchanged. Their accepted versions and S7 remain in `../polished/figures/`; they are not re-exported or resampled. S4/S6 are excluded entirely.

Compute-only QC rechecks exact plotted x/y and bar values against the accepted code, common-cohort membership, reference/available denominators, and all 14 available historical L1 cell-set identities. Accepted figure, caption and input-table hashes must remain unchanged. PDF proofs for all seven figures are rendered into `proofs/` at native page dimensions (178 mm wide).

Unresolved provenance: no archived historical UMAP coordinate hash proves equality to current regenerated coordinates. Shared coordinates within each historical comparison are documented by the original implementation. Historical unnormalized memberships cannot be recovered from reused paths; HECA 10% has no saved historical L1 list. These limitations remain disclosed; no new coverage claim or analysis is introduced.

Regenerate the two revised figures and all seven PDF proofs, on a compute node:

```bash
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
sbatch publication_figures/final_minimal/run_final.sbatch
```
