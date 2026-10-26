from common import *
import os,subprocess,json,hashlib
assert os.environ.get('SLURM_JOB_ID')
dest=OUT/'reports/random_DA_stage_A_2001';dest.mkdir(exist_ok=True,parents=True)
job='11137405'
cmd=['sacct','-j',job,'--units=K','--format=JobID,State,ExitCode,AllocCPUS,ReqMem,ElapsedRaw,MaxRSS,TotalCPU,Start,End,NodeList','-n','-P']
result=subprocess.run(cmd,capture_output=True,text=True,check=True)
(dest/'slurm_accounting.psv').write_text('JobID|State|ExitCode|AllocCPUS|ReqMem|ElapsedRaw|MaxRSS|TotalCPU|Start|End|NodeList\n'+result.stdout)
raw=pd.read_csv(dest/'slurm_accounting.psv',sep='|');resources=[]
for i in range(6):
    key=f'{job}_{i}';r=raw[raw.JobID==key]
    if len(r)!=1:continue
    r=r.iloc[0];batch=raw[raw.JobID==key+'.batch'];mem=batch.MaxRSS.iloc[0] if len(batch) else np.nan
    rss=float(str(mem).rstrip('K'))/1024**2 if pd.notna(mem) else np.nan
    resources.append(dict(task=i,atlas=ATLASES[i//2],subset_size=['02pct','10pct'][i%2],job_id=key,state=r.State,exit_code=r.ExitCode,CPUs=r.AllocCPUS,requested_memory=r.ReqMem,elapsed_seconds=r.ElapsedRaw,peak_batch_RSS_GiB=rss,TotalCPU=r.TotalCPU,start=r.Start,end=r.End,node=r.NodeList))
cost=pd.DataFrame(resources);cost.to_csv(dest/'resources.tsv',sep='\t',index=False)
p=pd.read_csv(dest/'pooled_metrics.tsv',sep='\t');elig=pd.read_csv(dest/'eligibility.tsv',sep='\t');typ=pd.read_csv(dest/'per_type_metrics.tsv',sep='\t');pair=pd.read_csv(dest/'paired_finite_effect_correlations.tsv',sep='\t')
previous=load('S3_pooled')
for rec in p[(p.method=='production')&(p.cohort=='original_reference')&(p.direction=='all')].itertuples():
    old=previous[(previous.atlas==rec.atlas)&(previous.subset_size==rec.subset_size)].iloc[0]
    assert rec.n_reference==old.n_atlas_reference_peak_cell_type_pairs and rec.n_available==old.n_effects_available
    assert np.isclose(rec.spearman,old.spearman_rho,rtol=0,atol=1e-8)
    assert np.isclose(rec.direction_agreement,old.direction_concordance,rtol=0,atol=1e-12)
    assert np.isclose(rec.strict_recovery,old.strict_recovery_fraction,rtol=0,atol=1e-12)
rare=pd.read_csv(dest/'rare_population_metrics.tsv',sep='\t');common=pd.read_csv(dest/'common_eligible_types.tsv',sep='\t');macro=pd.read_csv(dest/'macro_type_summaries.tsv',sep='\t')
def md(df):
    def fmt(v):return 'NA' if pd.isna(v) else f'{v:.4f}' if isinstance(v,(float,np.floating)) else str(v)
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
q=p[(p.cohort=='original_reference')&(p.direction=='all')][['atlas','subset_size','method','n_eligible_types','n_reference','n_available','spearman','direction_agreement','sign_switch_rate','strict_recovery']]
g=p[(p.cohort=='common_eligible')&(p.direction=='all')][['atlas','subset_size','method','n_cohort_types','n_reference','n_available','spearman','direction_agreement','strict_recovery']]
r=rare[(rare.cohort=='original_reference')&(rare.direction=='all')][['atlas','subset_size','method','cell_type','n_reference','n_available','spearman','direction_agreement','sign_switch_rate','strict_recovery']]
status=elig[['atlas','subset_size','method','cell_type','n_target_cells','n_rest_cells','status']]
focus={'kidney_atac':['Mes','Macula_Densa'],'heca_atac':['Haematopoietic stem and progenitor cell','Naive CD4 T cell'],'fetal_atac':['Megakaryocytes?','Ganglion cells']}
rs=pd.concat([status[(status.atlas==a)&status.cell_type.isin(ts)] for a,ts in focus.items()])
paired_pooled=pair[(pair.cell_type=='POOLED')&(pair.direction=='all')][['atlas','subset_size','cohort','method','n_identical_finite_pairs','spearman','median_absolute_error']]
direction=p[(p.cohort=='original_reference')&(p.direction!='all')][['atlas','subset_size','method','direction','n_reference','n_available','spearman','direction_agreement','sign_switch_rate','strict_recovery']]
text='''# Stage A matched-budget random DA pilot — seed 2001 only

## Scope and validation

Exactly six random DA runs: Kidney/HECA/Fetal × 2%/10%, seed 2001. Memberships were reconstructed from the completed benchmark, their recorded SHA-256 hashes verified, and all IDs validated against both Seurat and BPCells columns. The unchanged production script passed an 80-cell, 100-peak real-input smoke test before pilot submission. Pilot array 11137405 used global concurrency one, four CPUs, 128000 MB requested RAM and one hour per task. Preflight artifacts are in `publication_figures/benchmarks/random_DA_stage_A_2001/preflight/`; actual random outputs in `runs/{atlas}/seed2001/{fraction}/`. No production result or publication figure was overwritten.

Production implementation: `scripts/09_1_run_signac_da_all_cell_types.R` unchanged; BPCells Wilcoxon, >=10 target and >=20 rest cells, BH within type, adjusted P<=.05, |average log2FC|>=.25, detection >=.01 in either group, original pseudocount effect formula. Production subsets come from `results/signac_da/{atlas}/{fraction}/`; full references from `results/signac_da/{atlas}/100pct/signac_da_peaks.tsv`. Comparison verifies identical reference keys/directions/effects, threshold definitions and Boolean flags; see `definition_QA.tsv` and `sources.tsv`.

## Recorded computational outcomes

The table records actual elapsed seconds, Slurm batch peak RSS and CPU accounting, not requested memory as measured use. `slurm_accounting.psv` retains raw task/step records. TotalCPU is measured aggregate CPU time; four times wall time is allocated core time and differs from CPU utilization.

'''+md(cost[['atlas','subset_size','job_id','state','exit_code','CPUs','requested_memory','elapsed_seconds','peak_batch_RSS_GiB','TotalCPU']])+'''
## Original full-reference comparison

Pooled Spearman correlates signed effects over all available full-atlas significant peak–cell-type pairs. Availability denominator is all reference pairs. Direction agreement and switches divide by available finite effects; zeros enter neither numerator. Strict recovery divides by ALL full-reference pairs, including ineligible types. All rates are fractions, not percentages. Pooled estimates weight types with more reference peaks more heavily; do not infer per-type performance from them.

'''+md(q)+'''
## Fixed common-eligibility comparison

The cohort is the within-atlas intersection of types eligible in BOTH methods at BOTH budgets, for this single random seed. It is fixed across the four analyses, not separately redefined per budget or method. `common_eligible_types.tsv` lists membership. This deliberately excludes some of the rare populations whose eligibility differs, so the original-reference and eligibility tables remain essential. Reference denominators match across methods and budgets within this cohort.

'''+md(g)+'''
## Correlations on exactly the same finite reference pairs

Original-reference correlations above may use different available pairs. Here both methods correlate on the identical pairwise finite intersection. Per-type and open/closed versions are in `paired_finite_effect_correlations.tsv`. These correlations do not replace original strict-recovery denominators.

'''+md(paired_pooled)+'''
## Open and closed reference peaks

Direction is fixed by the full-atlas reference effect. Correlations within open/closed groups restrict the effect range and answer a different question from pooled signed correlations. Per-type counterparts are in `per_type_metrics.tsv`.

'''+md(direction)+'''
## Per-cell-type summaries

`per_type_metrics.tsv` contains every reference-bearing type, original/common cohorts and open/closed/all directions. `eligibility.tsv` retains ALL annotated types, including types with no reference peaks. Macro summaries below use equally weighted type medians; correlations with fewer than three finite pairs or zero variance are undefined. Undefined correlations are excluded from the median and their finite count is explicit. Types with no full-reference peaks have undefined recovery, not zero. These distributions are not biological replicates or confidence intervals.

'''+md(macro[macro.direction=='all'])+'''
## Enriched AND depleted focal populations

Target/rest counts and eligibility:

'''+md(rs)+'''
Full-reference metrics (NA correlation means insufficient or invariant available effects; missing availability must not be interpreted as concordance):

'''+md(r)+'''
## Interpretation limits and stopping rule

This is ONE random seed per atlas/budget, not a random-seed distribution. Do not infer statistical superiority, uncertainty bands or robustness across seeds. These results evaluate reproduction of a parent-atlas estimate, not external biological truth; parent/subset overlap affects concordance. Differences may involve eligibility, effect missingness, sample size, the group-size-dependent pseudocount and changes in one-versus-rest composition. The pilot cannot identify a causal biological or power-only explanation. An eligibility gain is a real evaluability gain; it is not by itself validated biological discovery. Assess pooled, macro/per-type and common-pair results together, and retain both enriched and depleted examples.

The original effect formula is preserved even for peaks undetected in both subset groups: with zero counts, its +1 count pseudocount can yield log2(n_rest/n_target), a nonzero effect driven by group sizes. These pairs fail the detection requirement for strict recovery but can enter available-effect direction/correlation metrics under the original definitions. `n_available_undetected_both` in pooled/per-type tables counts pairs with both recorded detection fractions zero. This is a definition qualification, not a new biological effect or a reason to change the approved analysis silently.

No additional seeds, Geosketch, coverage sensitivity jobs or publication figures were launched. Stop after reporting this pilot. Any further experiment requires explicit user approval.
'''
overview=[]
for cohort in ['original_reference','common_eligible']:
    w=p[(p.cohort==cohort)&(p.direction=='all')].pivot(index=['atlas','subset_size'],columns='method',values=['spearman','strict_recovery','direction_agreement'])
    rec=dict(cohort=cohort,n_atlas_budget_cases=len(w))
    for col in ['spearman','strict_recovery','direction_agreement']:
        rec[col+'_production_higher']=int((w[col]['production']>w[col]['random']).sum())
        rec[col+'_random_higher']=int((w[col]['random']>w[col]['production']).sum())
    overview.append(rec)
overview=pd.DataFrame(overview);overview.to_csv(dest/'descriptive_case_counts.tsv',sep='\t',index=False)
text+='\n## Descriptive cross-case summary\n\nCounts below summarize the six atlas/budget point comparisons for this one random seed. Higher is favorable for these three metrics. Counts are not significance tests, independent biological replications or evidence of robustness across random seeds. Original-reference correlations can use different available pairs; consult the identical-finite-pair table as well.\n\n'+md(overview)
if 'n_available_undetected_both' not in p.columns:
    text=text.replace('`n_available_undetected_both` in pooled/per-type tables counts pairs with both recorded detection fractions zero.','The frequency of this scenario was not separately tabulated in this completed pilot; its impact remains unquantified.')
text+='''
## Scientific interpretation of the completed pilot

1. **Strict recovery improves in HECA and Fetal, not Kidney.** Original-reference recovery (proposed/random) is Kidney 2% 8.92/9.76%, Kidney 10% 38.66/41.62%, HECA 2% 36.16/32.65%, HECA 10% 68.16/65.89%, Fetal 2% 27.45/19.33%, and Fetal 10% 59.78/46.99%. These are single-seed point comparisons, not method-level significance tests. Source: `pooled_metrics.tsv`, original_reference/all.

2. **The HECA/Fetal recovery gains are not solely an eligibility-denominator artifact.** With exactly the same fixed common-eligible reference pairs and complete effect availability in both methods, proposed strict recovery remains higher in both HECA and Fetal budgets. Fetal gains are 28.89 vs 22.23% at 2%, and 59.89 vs 49.07% at 10%. This rules out eligibility-related missing effects as the sole explanation for those cohort-specific differences, but does not isolate cell-number power, selected-cell characteristics, one-versus-rest composition or pseudocount effects. Source: `pooled_metrics.tsv`, common_eligible/all.

3. **Effect fidelity is not uniformly improved.** Within fixed common cohorts, random has higher pooled signed Spearman in four of six cases and higher pooled direction agreement in five of six. In HECA 2%, the apparent original pooled correlation advantage reverses on identical finite pairs (proposed .7871, random .7949). For Fetal 10%, the equal-type median signed correlation favors proposed (.8850 vs .8149), whereas pooled correlation favors random (.8967 vs .9173). Pooling weights peak-rich types and includes between-type structure; equal-type medians answer another question. Neither should substitute for the other. Sources: `paired_finite_effect_correlations.tsv`, `macro_type_summaries.tsv`, `pooled_metrics.tsv`.

4. **Direction-specific tradeoffs matter.** HECA 10% proposed strict recovery is higher for closed reference peaks (63.63 vs 58.38%) but lower for open peaks (77.24 vs 80.92%). Fetal strict-recovery gains occur for both open and closed reference peaks at both budgets. Report these alongside the signed pooled totals. Source: `pooled_metrics.tsv`, original_reference/open and closed.

5. **There are concrete rare-type gains and counterexamples.** Proposed HECA stem/progenitor strict recovery is 70.76% at 2% and 72.26% at 10%; random is ineligible at both budgets. Proposed Megakaryocytes? recovery is 13.44% and 68.88%; random is again ineligible. Conversely, random alone is eligible for naive CD4 T at 10% and ganglion at 2%. At 10%, Mes has almost equal target counts (17 proposed, 16 random) but markedly different strict recovery (16.37 vs 41.82%); target-cell count alone cannot explain every difference. This does not establish a biological mechanism or an external truth advantage. Sources: `eligibility.tsv`, `rare_population_metrics.tsv`.

**Conclusion:** this pilot supplies evidence of selective DA-recovery advantages that the global nearest-distance benchmark could not establish, alongside real effect-fidelity and depleted-population disadvantages. It supports a tradeoff-based scientific interpretation, not a claim of universal biological superiority. Robustness across random seeds remains unknown. Stop here as approved; no additional experiments have been launched.
'''
(dest/'REPORT.md').write_text(text)
qa=pd.read_csv(dest/'definition_QA.tsv',sep='\t')
assert len(qa)==12 and len(p)==72 and len(cost)==6
assert cost.state.eq('COMPLETED').all() and cost.exit_code.eq('0:0').all()
assert len(pd.read_csv(dest/'failures.tsv',sep='\t'))==0
executions=[json.loads(path.read_text()) for path in (OUT/'benchmarks/random_DA_stage_A_2001/runs').glob('*/seed2001/*/execution.json')]
assert len(executions)==6 and all(x['returncode']==0 for x in executions)
started=min(x['started_unix'] for x in executions)
figure_checks=[]
for path in sorted((OUT/'figures').glob('*')):
    if path.suffix not in ['.pdf','.svg','.png']:continue
    digest=hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda:fh.read(1048576),b''):digest.update(block)
    figure_checks.append(dict(path=str(path),mtime_precedes_pilot=path.stat().st_mtime<started,sha256=digest.hexdigest()))
pd.DataFrame(figure_checks).to_csv(dest/'figure_preservation_check.tsv',sep='\t',index=False)
(dest/'QA_summary.json').write_text(json.dumps(dict(DA_tasks=6,successful=6,failed=0,seed=2001,definition_checks=12,pooled_metric_rows=72,all_figure_mtimes_precede_pilot=all(r['mtime_precedes_pilot'] for r in figure_checks),further_experiments_submitted=False),indent=2))
print('Report written:',dest/'REPORT.md',flush=True)
print(q.to_string(index=False),flush=True)
print(g.to_string(index=False),flush=True)
print(rs.to_string(index=False),flush=True)
