# Manuscript evidence inventory — frozen existing results

Read-only synthesis. No new experiments, analysis reruns, algorithm changes or supplementary redesign. Evidence supports selective benefits and explicit tradeoffs; it does not establish universal representativeness or biological superiority. Stage A comprises one random DA seed, not replicated DA evidence.

## Theme 1 — L1 normalization versus unnormalized landmark features

| Supported finding | Existing evidence | Appropriate reference | Qualification |
|---|---|---|---|
| At 2%, L1 retains 24/24 Kidney, 101/101 HECA and 78/78 Fetal types, versus historical unnormalized counts 22/24, 90/101 and 74/78 | `publication_figures/tables/S1_normalization_2pct.tsv`; original `diagnostics/rbf_normalization_ablation/{atlas}/comparison/method_summary.tsv` | S1 A | One recorded selection per configuration; unnormalized Apricot is not uniform random sampling |
| Historical composition TV decreases from .2436 to .1150 Kidney, .2231 to .1420 HECA, .2000 to .1191 Fetal under L1 | Same S1 source table | S1 B | Improvement relative to this unnormalized implementation, not proof of random-baseline superiority |
| Kidney/Fetal multibudget normalization summaries exist | `S1_normalization_budgets.tsv` | S1 D | HECA multibudget normalization comparison is absent |

Historical 2D UMAP distances in S1 C are not 60D coverage. Within historical comparisons the geometry was shared, but no archived coordinate-content hash establishes equality to current regenerated coordinates. Fourteen available historical L1 ID lists match production; HECA 10% historical list is absent. Original unnormalized membership paths were reused and cannot establish current baseline membership. Cite the existing S1 caption and `S1_current_identity_check.tsv`; do not inflate this into independent replication.

## Theme 2 — complete presence and selective exceptional-rarity enrichment

| Supported finding | Existing evidence | Appropriate reference | Qualification |
|---|---|---|---|
| All 24/101/78 annotated types are present at every tested 2/4/6/8/10% budget | `publication_figures/tables/S2_composition.tsv`, `S2_all_population_metrics.tsv` | S2 caption; S7 descriptive overlays | Presence is not sufficient cell count for DA, proportional fidelity or quantitative coverage |
| At 2%, HECA stem/progenitor: 14 of 56 cells; Fetal Megakaryocytes?: 11 of 60; enrichment 12.50× and 9.17× | S6b source table and `reports/scientific_followup/rare_population_comparisons.tsv` | Candidate Figure 1 A/B; frozen S6b | Deliberately chosen examples, not universal rare-population protection; uncertain annotation retained |
| No measured random seed meets eligibility for those two populations at 2%; naive CD4 is eligible in 1/20 seeds and ganglion in 17/20 | Same 20-seed source table | Candidate Figure 1 B; S6b | Empirical frequency among the measured seeds, not population-level probability or CI |
| Counterexamples: at 2%, proposed naive CD4 count 1 versus random median 5; ganglion 1 versus median 11.5 | Same source table | Candidate Figure 1 A/B; S6b | Both methods are ineligible for naive CD4 in the displayed DA seed; ganglion is a direct opposing eligibility example |
| Complete 2% retention occurs in random 20/20 Kidney, 5/20 HECA, 15/20 Fetal runs | `reports/scientific_followup/benchmark_comparisons.tsv` | S6; evidence tables | Complete presence is most distinguishing for HECA here; one proposed run, not replicated method reliability |
| Random sampling better preserves composition and typical-cell nearest-selected proximity | `tables/S4_S6_measured_benchmark_summary.tsv`; `reports/scientific_followup/benchmark_comparisons.tsv` | Frozen S4/S6 | Exact composition contrasts; approximate NN coverage with finite-sample calibration and occasional large misses; fixed-radius uncertainty not included in seed bands |

The source files above are under `publication_figures/`. Evaluating annotations after selection does not make annotations selection inputs. Retain independence of the five production selections; do not depict nested prefixes. Source implementations and summary metadata: `scripts/03_subset_selection_on_peakvi.py`, `results/subset_selection/peakvi/{atlas}/subset_selection_summary_peakvi_{fraction}.tsv`.

## Theme 3 — downstream DA gains and losses

| Supported finding | Existing evidence | Appropriate reference | Qualification |
|---|---|---|---|
| Original-reference strict recovery at 2%: proposed/random HECA 36.16/32.65%, Fetal 27.45/19.33%, Kidney 8.92/9.76% | `reports/random_DA_stage_A_2001/pooled_metrics.tsv` | Stage A evidence tables, candidate C for focal populations only | One random seed; these are whole-atlas pooled reference-pair rates, not the four-population bars |
| At 10%: HECA 68.16/65.89%, Fetal 59.78/46.99%, Kidney 38.66/41.62% | Same table | Stage A evidence tables | No inference across untested random seeds or intermediate DA budgets |
| HECA/Fetal strict-recovery gains persist in common-eligible cohorts; fixed membership counts 21 Kidney / 83 HECA / 64 Fetal | `common_eligible_types.tsv`, `pooled_metrics.tsv`, `paired_finite_effect_correlations.tsv` | Stage A tables; S3 explains production metric concepts | Cohort is the intersection of BOTH methods at BOTH pilot budgets. Earlier production-only S3 cohort is different |
| Effect fidelity is mixed: common-cohort random pooled Spearman higher in 4/6 cases, direction agreement higher in 5/6 | `descriptive_case_counts.tsv`, pooled and paired-finite tables | Stage A tables; frozen S3/S3b are production-only context | Pooled signed, within-type signed and within-direction correlations differ; no method-level significance test |
| At 2%, proposed HSPC strict recovery 70.76%, Megakaryocytes? 13.44%; seed-2001 random is ineligible. Random ganglion recovers 12.27%; proposed is ineligible | `rare_population_metrics.tsv`, `eligibility.tsv` | Candidate Figure 1 C | Ineligible is unavailable DA, not measured zero biological fidelity. Recovery is relative to the parent-atlas significant reference |
| Other rare-type disadvantages remain: naive CD4 random eligible only at 10%; Mes 10% strict recovery proposed/random 16.37/41.82% despite target counts 17/16 | Same tables | Stage A rare-population evidence tables; S2/S6b representation context | Cell count alone is not a sufficient mechanism; does not isolate power, rest composition or selected-cell differences |

## Required disclosures

- No universal superiority, proportional-representativeness or uniformly improved 60D coverage claim. The optimized feature objective is not independent biological validation.
- Twenty random seeds support representation/eligibility variability; random DA uses seed 2001 only at 2%/10%. Proposed selections have one production run per budget. No invented uncertainty intervals or across-seed DA claims.
- Examples were selected for illustration, including gains and losses. Most evidence is relative to annotated populations and parent-atlas reference estimates; there is no external biological truth validation here.
- DA uses the unchanged BPCells Wilcoxon implementation within Seurat/Signac objects: >=10 target/>=20 rest, BH within type across matrix peaks, adjusted P<=.05, |average log2FC|>=.25, detection >=.01 in either group. Full references are fixed. Strict recovery includes ineligible types in the original denominator; direction agreement/correlation require available effects.
- Pooled peak-pair weighting can disagree with equally weighted per-type summaries. Shared finite pairs and fixed common eligibility mitigate availability confounding but do not equalize per-type sample size or isolate biological mechanisms.
- Parent/subset overlap, changing one-versus-rest composition and the group-size-dependent +1 pseudocount affect interpretation. Reference recovery is recall-like and is not external specificity or validation of newly discovered peaks. Do not attribute differences solely to power.
- Historical residual-gain diagnostics were calculated on earlier/raw-feature selections and cannot establish the current L1 mechanism. Saved density diagnostics have provenance limitations. Do not assert a proven mechanism of rare-population preference.
- S4 distances are approximate NN upper bounds with 128-query validation; the audit found occasional severe sampled misses despite the gate. Radii are atlas-specific and fixed, not biological thresholds. Nonselected-cohort query identities differ between methods. Composition measurements do not depend on this ANN caveat.
- S7 is descriptive 2D UMAP, with deterministic display sampling capped at 3,000 selected and 75,000 background points. Actual n and shown are distinct. S8 uses display-only clipping; correlations use unmodified effects.
- Runtimes are individual successful production/pilot runs, with different inputs/load/hardware; no controlled scaling, speed superiority or repeated timing uncertainty is established.

## Defensible manuscript framing

“Annotation-blind L1-normalized submodular feature selection selectively increases representation of some exceptionally rare annotated populations, enabling DA testing at low global sampling budgets. These gains coexist with depleted populations, worse global composition/proximity fidelity and dataset-/population-dependent DA tradeoffs.”

This framing is descriptive of current measurements. It does not replace the need to disclose the single-seed DA pilot or establish generalizable superiority. BioRender supplies the separate graphical abstract; candidate Figure 1 supplies the focused quantitative illustration. No further experiments or figure redesign are authorized by this preparation task.
