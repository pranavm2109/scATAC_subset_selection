# Publication cleanup report

- Full pre-change inventory: 1571 files, 7,381,798,470 bytes. Every publication file was archived and its archived contents verified against SHA-256; manuscript-preparation source materials were included too. Archive: `archive/pre_reorganization_20261018_235114`.
- Definitive Figure 1: renamed exports are byte-identical to reviewed originals, with unchanged three-panel layout and measurements. Original candidate files retained.
- S9: regenerated only the original four-panel overview from exact copied input tables; no analysis or metrics recomputed. Panel code unchanged except function/data/output designations. Original overview exports, script and inputs retained.
- Ten finalized supplementary figures copied byte-identically. Their original exports, captions and all scientific inputs/scripts/reports were verified unchanged after packaging.
- Figure index maps scripts, source tables and captions, with an optional sequential-label proposal. Main and S9 are regular final exports; accepted supplementary copies retain their current names and artwork.
- No deletion executed. Exact-content duplicate inventory is descriptive only; DELETION_MANIFEST.tsv lists proposed obsolete export locations and replacements. Scientific data, memberships, seed records, scripts, reports and provenance are not deletion candidates.
- Existing relative paths are retained. New script/table folders link to those dependencies instead of moving them.
- Export QC: all 12 PDF/SVG/PNG sets checked for 178-mm width, embedded TrueType PDF fonts, parseable SVG and 600-DPI PNG. All PDFs rendered to proofs at native page dimensions. Human visual review is recorded separately after the job.

## Scientific qualifications retained

20 random seeds inform representation/eligibility; random DA seed 2001 only. DA-ineligible is distinct from measured zero recovery. Rare examples are selective, not universal gains. Historical 2D UMAP normalization evidence is not current 60D coverage. Benchmark distances use approximate search with recorded accuracy limits. S9 panel C is direction agreement, not effect-size correlation. Full-reference recovery denominators and changing effect availability are preserved. No new biological, causal or superiority claim follows from packaging.

## Open editorial matters

Optional sequential numbering is not applied; authors must choose one scheme consistently. The historical S4/S6 caption is still a template and is preserved; a current measured-benchmark technical reading guide accompanies the exports and links the scientific audit. Final journal prose and the BioRender license/version need author review. Historical README/plans remain accessible and clearly superseded by PUBLICATION_README.md/FIGURE_INDEX.md.

## Completion and execution record

All twelve PDF proofs were individually inspected; see `VISUAL_REVIEW.md`. S6b relies on the explicit companion caption line key because its accepted artwork has no inline method legend. No accepted artwork was changed to address this.

Job 11138961 completed the inventory and verified archive, then stopped on a missing `grouped` plotting-helper import while rendering S9. Job 11138965 resumed from that same verified archive after correcting the wrapper and completed packaging/QC. Both logs are preserved under `manuscript_preparation/logs/`. Neither job ran analyses or scientific experiments. The final standalone S9 wrapper routes the original save/layout helper into temporary staging and then the final supplementary folder; it does not depend on the one-time packaging work directory. The optional compute-node reproduction command is documented in PUBLICATION_README.md. No further jobs were submitted.

Job: 11138965. Timezone for version stamp: America/New_York.
