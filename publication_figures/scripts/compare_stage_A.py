"""Summarize only the approved six DA runs; no figures or additional testing."""
from common import *
import os,json,gc,hashlib,subprocess
from scipy.stats import spearmanr,pearsonr
assert os.environ.get('SLURM_JOB_ID')
base=OUT/'benchmarks/random_DA_stage_A_2001';dest=OUT/'reports/random_DA_stage_A_2001';dest.mkdir(parents=True,exist_ok=True)
assert hashlib.sha256((ROOT/'scripts/09_1_run_signac_da_all_cell_types.R').read_bytes()).hexdigest()==(base/'preflight/production_script.sha256').read_text().strip()
def put(name,frame):frame.to_csv(dest/f'{name}.tsv',sep='\t',index=False)
def cor(x,y,method='spearman'):
    if len(x)<3 or np.std(x)==0 or np.std(y)==0:return np.nan
    return float((spearmanr if method=='spearman' else pearsonr)(x,y).statistic)
def summarize(t):
    available=t.effect_available.to_numpy();x=t.full_avg_log2FC.to_numpy()[available];y=t.subset_avg_log2FC.to_numpy()[available];nr=len(t);na=len(x)
    return dict(n_reference=nr,n_available=na,availability=na/nr if nr else np.nan,n_available_undetected_both=int((available&t.subset_pct_target.eq(0)&t.subset_pct_rest.eq(0)).sum()),n_same=int(t.same_direction_effect.sum()),n_switch=int(t.direction_switched_effect.sum()),n_strict=int(t.strict_same_direction_capture.sum()),direction_agreement=t.same_direction_effect.sum()/na if na else np.nan,sign_switch_rate=t.direction_switched_effect.sum()/na if na else np.nan,strict_recovery=t.strict_same_direction_capture.sum()/nr if nr else np.nan,spearman=cor(x,y),pearson=cor(x,y,'pearson'),median_absolute_error=float(np.median(np.abs(x-y))) if na else np.nan)
eligibility=[];cohorts=[];metrics=[];pertype=[];paired=[];checks=[];sources=[];failures=[]
for a in ATLASES:
    statuses={};eligible_sets=[]
    for f in ['02pct','10pct']:
        for method in ['production','random']:
            folder=ROOT/f'results/signac_da/{a}/{f}' if method=='production' else base/f'runs/{a}/seed2001/{f}'
            if not (folder/'signac_da_status.tsv').exists():failures.append(dict(atlas=a,subset_size=f,method=method,error='Missing status output'));continue
            s=pd.read_csv(folder/'signac_da_status.tsv',sep='\t');statuses[method,f]=s
            assert not s.cell_type.duplicated().any()
            assert s.fdr_threshold.eq(.05).all() and s.abs_log2fc_threshold.eq(.25).all() and s.min_pct.eq(.01).all()
            assert s.test.eq('wilcox').all() and s.multiple_testing.eq('BH_within_cell_type').all()
            eligible=(s.n_target_cells>=10)&(s.n_rest_cells>=20)
            assert np.array_equal(s.status.eq('success'),eligible)
            eligible_sets.append(set(s.loc[eligible,'cell_type']))
            eligibility.extend(s.assign(method=method,subset_size=f,seed=123 if method=='production' else 2001).to_dict('records'))
            sources.append(dict(atlas=a,subset_size=f,method=method,status_source=str(folder/'signac_da_status.tsv'),effects_source=str(folder/'signac_da_reference_effects.tsv')))
    if len(statuses)!=4:continue
    common=set.intersection(*eligible_sets)
    cohorts.extend(dict(atlas=a,cell_type=t,definition='eligible in both methods at both 2% and 10%; seed2001 only') for t in sorted(common))
    reference=pd.read_csv(ROOT/f'results/signac_da/{a}/100pct/signac_da_peaks.tsv',sep='\t',usecols=['cell_type','peak','avg_log2FC','direction','adjusted_p_value','pct_target','pct_rest'])
    assert not reference.duplicated(['cell_type','peak']).any()
    assert reference.adjusted_p_value.le(.05).all() and reference.avg_log2FC.abs().ge(.25).all()
    assert np.maximum(reference.pct_target,reference.pct_rest).ge(.01).all()
    for f in ['02pct','10pct']:
        pair={}
        for method in ['production','random']:
            folder=ROOT/f'results/signac_da/{a}/{f}' if method=='production' else base/f'runs/{a}/seed2001/{f}'
            path=folder/'signac_da_reference_effects.tsv'
            if not path.exists():failures.append(dict(atlas=a,subset_size=f,method=method,error='Missing reference effects'));continue
            cols=['cell_type','peak','full_avg_log2FC','subset_avg_log2FC','full_direction','effect_available','same_direction_effect','direction_switched_effect','strict_same_direction_capture','subset_is_significant','subset_adjusted_p_value','subset_pct_target','subset_pct_rest','subset_analysis_status']
            t=pd.read_csv(path,sep='\t',usecols=cols,dtype={'cell_type':'category','full_direction':'category','subset_analysis_status':'category'})
            assert len(t)==len(reference)
            assert np.array_equal(t.cell_type.astype(str).to_numpy(),reference.cell_type.to_numpy()) and np.array_equal(t.peak.to_numpy(),reference.peak.to_numpy())
            assert np.allclose(t.full_avg_log2FC,reference.avg_log2FC,rtol=1e-12,atol=1e-12)
            assert np.array_equal(t.full_direction.astype(str),reference.direction)
            x=t.full_avg_log2FC.to_numpy();y=t.subset_avg_log2FC.to_numpy();available=np.isfinite(y)&t.subset_analysis_status.eq('success').to_numpy()
            assert np.array_equal(t.effect_available,available)
            same=available&(np.sign(x)==np.sign(y));switch=available&(np.sign(x)!=np.sign(y))&(y!=0)
            sig=available&(np.maximum(t.subset_pct_target,t.subset_pct_rest)>=.01)&np.isfinite(t.subset_adjusted_p_value)&(t.subset_adjusted_p_value<=.05)&(np.abs(y)>=.25)
            assert np.array_equal(t.same_direction_effect,same) and np.array_equal(t.direction_switched_effect,switch)
            assert np.array_equal(t.subset_is_significant,sig) and np.array_equal(t.strict_same_direction_capture,sig&same)
            for cohort in ['original_reference','common_eligible']:
                u=t if cohort=='original_reference' else t[t.cell_type.isin(common)]
                for direction in ['all','open','closed']:
                    v=u if direction=='all' else u[u.full_direction==direction]
                    metrics.append(dict(atlas=a,subset_size=f,method=method,cohort=cohort,direction=direction,n_eligible_types=len(statuses[method,f].query("status=='success'")),n_cohort_types=reference.cell_type.nunique() if cohort=='original_reference' else len(common),**summarize(v)))
                for typ,v in u.groupby('cell_type',observed=True,sort=False):
                    for direction in ['all','open','closed']:
                        w=v if direction=='all' else v[v.full_direction==direction]
                        pertype.append(dict(atlas=a,subset_size=f,method=method,cohort=cohort,cell_type=typ,direction=direction,**summarize(w)))
            pair[method]=t[['cell_type','full_direction','full_avg_log2FC','subset_avg_log2FC','effect_available']].copy()
            checks.append(dict(atlas=a,subset_size=f,method=method,n_reference=len(t),reference_keys_effects_identical=True,flags_thresholds_verified=True))
            del t;gc.collect()
        if len(pair)==2:
            p=pair['production'];r=pair['random'];both=p.effect_available.to_numpy()&r.effect_available.to_numpy();cm=p.cell_type.isin(common).to_numpy()
            for cohort,mask in [('original_reference',both),('common_eligible',both&cm)]:
                for typ in ['POOLED']+list(reference.cell_type.unique()):
                    tm=mask if typ=='POOLED' else mask&p.cell_type.eq(typ).to_numpy()
                    for direction in ['all','open','closed']:
                        dm=tm if direction=='all' else tm&p.full_direction.eq(direction).to_numpy()
                        x=p.full_avg_log2FC.to_numpy()[dm]
                        for method,t in pair.items():
                            y=t.subset_avg_log2FC.to_numpy()[dm]
                            paired.append(dict(atlas=a,subset_size=f,cohort=cohort,cell_type=typ,direction=direction,method=method,n_identical_finite_pairs=len(x),spearman=cor(x,y),pearson=cor(x,y,'pearson'),median_absolute_error=float(np.median(np.abs(x-y))) if len(x) else np.nan))
            del pair,p,r;gc.collect()
        print('Compared',a,f,'common cohort',len(common),flush=True)
    del reference;gc.collect()
put('eligibility',pd.DataFrame(eligibility));put('common_eligible_types',pd.DataFrame(cohorts));put('pooled_metrics',pd.DataFrame(metrics));put('per_type_metrics',pd.DataFrame(pertype));put('paired_finite_effect_correlations',pd.DataFrame(paired));put('definition_QA',pd.DataFrame(checks));put('sources',pd.DataFrame(sources));put('failures',pd.DataFrame(failures,columns=['atlas','subset_size','method','error']))
per=pd.DataFrame(pertype)
macro=per.groupby(['atlas','subset_size','method','cohort','direction'],observed=True).agg(n_types=('cell_type','nunique'),n_finite_spearman=('spearman','count'),median_per_type_spearman=('spearman','median'),median_per_type_strict_recovery=('strict_recovery','median'),median_per_type_direction_agreement=('direction_agreement','median')).reset_index()
put('macro_type_summaries',macro)
focus={'kidney_atac':['Mes','Macula_Densa'],'heca_atac':['Haematopoietic stem and progenitor cell','Naive CD4 T cell'],'fetal_atac':['Megakaryocytes?','Ganglion cells']}
rare=pd.concat([per[(per.atlas==a)&per.cell_type.isin(types)] for a,types in focus.items()]);put('rare_population_metrics',rare)
print('Failures:',failures,flush=True)
wide=pd.DataFrame(metrics)
keys=['atlas','subset_size','cohort','direction']
delta=wide[wide.method=='production'].merge(wide[wide.method=='random'],on=keys,suffixes=('_production','_random'),validate='one_to_one')
assert np.array_equal(delta.n_reference_production,delta.n_reference_random)
for col in ['availability','spearman','pearson','direction_agreement','sign_switch_rate','strict_recovery','median_absolute_error']:
    delta[col+'_production_minus_random']=delta[col+'_production']-delta[col+'_random']
put('pooled_method_differences',delta)
print(pd.DataFrame(metrics).query("direction=='all'").to_string(index=False),flush=True)
