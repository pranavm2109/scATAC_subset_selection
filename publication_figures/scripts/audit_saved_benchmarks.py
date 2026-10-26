"""Read saved results and verify memberships; no NN fitting, selection or DA."""
from common import *
from pathlib import Path
import os,json,hashlib
from scipy.stats import spearmanr
assert os.environ.get('SLURM_JOB_ID'),'Compute-node batch only'
dest=OUT/'reports'/'scientific_followup';dest.mkdir(parents=True,exist_ok=True)
def put(name,d):d.to_csv(dest/f'{name}.tsv',sep='\t',index=False)
base=OUT/'benchmarks'/'full'
paths=sorted(base.glob('*_metrics.tsv'));assert len(paths)==315
d=pd.concat([pd.read_csv(p,sep='\t') for p in paths],ignore_index=True)
assert len(d)==630 and not d.duplicated(['atlas','subset_size','method','seed','cohort']).any()
assert d.query_mode.eq('full_atlas').all() and d.fraction_within_one_percent.ge(.99).all()
assert d.n_exact_validation.eq(128).all()
pops=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(base.glob('*_populations.tsv'))],ignore_index=True)
rows=[];identity=[];provenances=[]
for atlas in ATLASES:
    names=pd.read_csv(ROOT/f'results/peakvi/{atlas}/cell_names.txt',header=None,dtype=str)[0].to_numpy()
    assert pd.Index(names).is_unique
    latent=np.load(ROOT/f'results/peakvi/{atlas}/peakvi_latent.npy',mmap_mode='r');assert latent.shape==(len(names),60)
    meta=read(f'results/cell_metadata/{atlas}/{atlas}_cell_metadata.tsv',dtype=str).set_index('cell_id')
    labels=meta.loc[names,'cell_type'];assert labels.notna().all()
    codes,types=pd.factorize(labels,sort=True);fullcounts=np.bincount(codes,minlength=len(types))
    radius=pd.read_csv(OUT/f'benchmarks/radius/{atlas}_radius_production_02pct_seed2001_radius.tsv',sep='\t').iloc[0]
    for f in FRACS:
        q=d[(d.atlas==atlas)&(d.subset_size==f)&(d.cohort=='all_queries')]
        assert len(q)==21 and set(q[q.method=='random'].seed)==set(range(2001,2021))
        assert q.n_selected.nunique()==1 and q.n_atlas.eq(len(names)).all() and q.n_query.eq(len(names)).all()
        assert np.allclose(q.fixed_radius,radius.radius) and q.ann_index_seed.eq(1705).all() and q.threads.eq(1).all()
        k=int(q.n_selected.iloc[0]);pct=int(f[:2])
        for r in q.itertuples():
            if r.method=='production':
                ids=pd.read_csv(ROOT/f'results/subset_selection/peakvi/{atlas}/selected_cells_peakvi_{f}.txt',header=None,dtype=str)[0]
                assert not ids.duplicated().any();chosen=pd.Index(names).get_indexer(ids);assert (chosen>=0).all();chosen=np.sort(chosen)
            else:
                chosen=np.sort(np.random.default_rng(np.random.SeedSequence([int(r.seed),ATLASES.index(atlas),pct])).choice(len(names),k,replace=False))
            assert len(chosen)==len(np.unique(chosen))==k
            token=hashlib.sha256('\n'.join(names[chosen]).encode()).hexdigest();assert token==r.selected_ID_sha256
            pop=pops[(pops.atlas==atlas)&(pops.subset_size==f)&(pops.method==r.method)&(pops.seed==r.seed)].set_index('cell_type').loc[types]
            count=np.bincount(codes[chosen],minlength=len(types));assert np.array_equal(count,pop.subset_count)
            assert np.array_equal(fullcounts,pop.full_count)
            non=d[(d.atlas==atlas)&(d.subset_size==f)&(d.method==r.method)&(d.seed==r.seed)&(d.cohort=='nonselected_queries')].iloc[0]
            assert non.n_query==len(names)-k
            prov=json.loads((base/f'{atlas}_evaluate_{r.method}_{f}_seed{r.seed}_provenance.json').read_text())
            assert prov['latent_path']==str(ROOT/f'results/peakvi/{atlas}/peakvi_latent.npy')
            assert prov['coordinate_transform']=='none; saved raw 60D Euclidean PeakVI'
            provenances.append(dict(atlas=atlas,subset_size=f,method=r.method,seed=r.seed,versions=json.dumps(prov['versions'],sort_keys=True)))
            identity.append(dict(atlas=atlas,subset_size=f,method=r.method,seed=r.seed,n_selected=k,hash_matches=True,counts_match=True,IDs_valid=True))
        for cohort in ['all_queries','nonselected_queries']:
            u=d[(d.atlas==atlas)&(d.subset_size==f)&(d.cohort==cohort)];pro=u[u.method=='production'].iloc[0];rnd=u[u.method=='random']
            record=dict(atlas=atlas,subset_size=f,cohort=cohort,n_selected=k,radius=radius.radius,n_random_seeds=20)
            for metric in ['n_retained','n_DA_eligible','TV','JS_bits','median_distance','q90_distance','q99_distance','coverage_within_radius']:
                record[metric+'_proposed']=pro[metric];record[metric+'_random_median']=rnd[metric].median();record[metric+'_random_q025']=rnd[metric].quantile(.025);record[metric+'_random_q975']=rnd[metric].quantile(.975)
                record[metric+'_random_min']=rnd[metric].min();record[metric+'_random_max']=rnd[metric].max()
            record['random_complete_retention_runs']=int(rnd.n_retained.eq(pro.n_types).sum())
            rows.append(record)
    print('Validated memberships and counts:',atlas,flush=True)
put('benchmark_validation',pd.DataFrame(identity));put('benchmark_protocol_versions',pd.DataFrame(provenances));put('benchmark_comparisons',pd.DataFrame(rows))
# Verify every numerical S4/S6 summary against underlying metrics.
summary=load('S4_S6_measured_benchmark_summary')
for r in summary.itertuples():
    q=d[(d.atlas==r.atlas)&(d.subset_size==r.subset_size)&(d.method==r.method)&(d.cohort=='all_queries')][r.metric]
    assert len(q)==r.n_seeds and np.isclose(q.median(),r.median) and np.isclose(q.quantile(.025),r.seed_q025) and np.isclose(q.quantile(.975),r.seed_q975)
print('PASS: all panel summary values reproduced',flush=True)
ct=load('S2_all_population_metrics');focus={'heca_atac':['Haematopoietic stem and progenitor cell','Naive CD4 T cell'],'fetal_atac':['Megakaryocytes?','Ganglion cells'],'kidney_atac':['Mes','Macula_Densa']}
rare=[];historical=[];associations=[]
for atlas in ATLASES:
    for f in FRACS:
        current=ct[(ct.atlas==atlas)&(ct.subset_size==f)]
        for typ in focus[atlas]:
            pro=current[current.cell_type==typ].iloc[0];rnd=pops[(pops.atlas==atlas)&(pops.subset_size==f)&(pops.method=='random')&(pops.cell_type==typ)]
            rare.append(dict(atlas=atlas,subset_size=f,cell_type=typ,full_count=pro.full_count,full_percent=100*pro.full_proportion,proposed_count=pro.subset_count,enrichment=pro.enrichment,DA_eligible=pro.DA_eligible,random_count_median=rnd.subset_count.median(),random_count_q025=rnd.subset_count.quantile(.025),random_count_q975=rnd.subset_count.quantile(.975),random_count_min=rnd.subset_count.min(),random_count_max=rnd.subset_count.max(),random_eligible_seeds=int(rnd.DA_eligible.sum())))
        covariates={'log10_full_count':np.log10(current.full_count.to_numpy())}
        geometry=ROOT/f'diagnostics/peakvi_celltype_separability/{atlas}/cell_type_separability.tsv'
        if geometry.exists():
            g=pd.read_csv(geometry,sep='\t').set_index('cell_type').reindex(current.cell_type)
            for col in ['median_knn_distance','median_knn_radius','knn_purity','silhouette']:covariates[col]=g[col].to_numpy()
        features=pd.read_csv(ROOT/f'diagnostics/rbf_activation_marginal_gain/{atlas}/cell_type_activation_marginal_gain.tsv',sep='\t').set_index('cell_type').reindex(current.cell_type)
        for col in ['median_total_rbf_mass','median_effective_landmark_support','median_initial_marginal_gain','n_landmarks']:covariates['historical_geometry_'+col]=features[col].to_numpy()
        for name,x in covariates.items():
            y=current.log2_enrichment.to_numpy();ok=np.isfinite(x)&np.isfinite(y)
            associations.append(dict(atlas=atlas,subset_size=f,covariate=name,n_types=int(ok.sum()),spearman_rho=spearmanr(x[ok],y[ok]).statistic,interpretation='descriptive across types; historical geometry provenance not hash verified' if name!='log10_full_count' else 'current counts; descriptive'))
    for kind,path in [('raw_activation',ROOT/f'diagnostics/rbf_activation_marginal_gain/{atlas}/cell_type_activation_marginal_gain.tsv'),('feature_behavior',ROOT/f'diagnostics/featurebased_selection_behavior/{atlas}/cell_type_featurebased_diagnostics.tsv')]:
        if not path.exists():continue
        g=pd.read_csv(path,sep='\t')
        for typ in focus[atlas]:
            t=g[g.cell_type==typ]
            if len(t):
                rec=t.iloc[0].to_dict();rec.update(diagnostic_kind=kind,source=str(path.relative_to(ROOT)))
                for f in FRACS:
                    if f'selected_n_{f}' in rec:
                        n=int(ct[(ct.atlas==atlas)&(ct.cell_type==typ)&(ct.subset_size==f)].subset_count.iloc[0]);rec[f'current_n_{f}']=n;rec[f'historical_count_matches_{f}']=rec[f'selected_n_{f}']==n
                historical.append(rec)
put('rare_population_comparisons',pd.DataFrame(rare));put('descriptive_associations',pd.DataFrame(associations));put('historical_diagnostic_focus',pd.DataFrame(historical))
put('all_saved_metrics',d)
print('Saved read-only audit derivatives in',dest,flush=True)
print(pd.DataFrame(rows).query("cohort=='all_queries'")[["atlas","subset_size","TV_proposed","TV_random_median","q90_distance_proposed","q90_distance_random_median","coverage_within_radius_proposed","coverage_within_radius_random_median","random_complete_retention_runs"]].to_string(index=False),flush=True)
print(pd.DataFrame(rare).query("subset_size in ['02pct','10pct']").to_string(index=False),flush=True)
