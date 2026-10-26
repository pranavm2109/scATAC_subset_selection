from common import *
import os,hashlib
from scipy.stats import spearmanr
assert os.environ.get('SLURM_JOB_ID')
dest=OUT/'reports/scientific_followup'
def put(n,d):d.to_csv(dest/f'{n}.tsv',sep='\t',index=False)
ct=load('S2_all_population_metrics');focus={'heca_atac':['Haematopoietic stem and progenitor cell','Naive CD4 T cell'],'fetal_atac':['Megakaryocytes?','Ganglion cells'],'kidney_atac':['Mes','Macula_Densa']}
cells=[];groups=[];assoc=[]
for atlas in ATLASES:
    names=pd.read_csv(ROOT/f'results/peakvi/{atlas}/cell_names.txt',header=None,dtype=str)[0].to_numpy()
    labels=read(f'results/cell_metadata/{atlas}/{atlas}_cell_metadata.tsv',dtype=str).set_index('cell_id').loc[names,'cell_type'].to_numpy()
    raw=np.load(ROOT/f'diagnostics/rbf_activation_marginal_gain/{atlas}/per_cell_activation_marginal_gain.npz')
    ordered=pd.read_csv(ROOT/str(raw['cell_name_order_file'][0]),header=None,dtype=str)[0].to_numpy();assert np.array_equal(names,ordered)
    mass=raw['total_rbf_mass'];gain=raw['initial_marginal_gain'];support=raw['effective_landmark_support'];assert len(mass)==len(names)
    prod=read(f'results/subset_selection/peakvi/{atlas}/subset_selection_summary_peakvi_02pct.tsv').iloc[0]
    values={'raw_RBF_mass':mass,'effective_landmark_support':support,'L1_initial_gain_from_saved_raw':gain*np.sqrt(float(prod.rbf_row_sum_target)/mass)}
    nn=ROOT/f'diagnostics/peakvi_celltype_separability/{atlas}/per_cell_knn_sparsity.npz'
    if nn.exists():
        z=np.load(nn);h=hashlib.sha256()
        for n in names:h.update(str(n).encode());h.update(b'\n')
        assert str(z['cell_name_sha256'])==h.hexdigest()
        values['mean_30NN_distance']=z['mean_knn_distance'];values['30NN_radius']=z['knn_radius']
    for f in ['02pct','10pct']:
        ids=pd.read_csv(ROOT/f'results/subset_selection/peakvi/{atlas}/selected_cells_peakvi_{f}.txt',header=None,dtype=str)[0]
        chosen=pd.Index(names).get_indexer(ids);assert (chosen>=0).all();mask=np.zeros(len(names),dtype=bool);mask[chosen]=True
        for key,x in values.items():
            cells.append(dict(atlas=atlas,subset_size=f,metric=key,full_median=float(np.median(x)),selected_median=float(np.median(x[mask])),nonselected_median=float(np.median(x[~mask])),provenance='Saved geometry; cell order verified; latent/feature content hashes unavailable'))
        for typ in focus[atlas]:
            sub=labels==typ;r=ct[(ct.atlas==atlas)&(ct.subset_size==f)&(ct.cell_type==typ)].iloc[0]
            rec=dict(atlas=atlas,subset_size=f,cell_type=typ,current_count=r.subset_count,current_enrichment=r.enrichment)
            for key,x in values.items():rec[key+'_population_median']=float(np.median(x[sub]))
            lm=read(f'diagnostics/landmark_celltype_coverage/{atlas}/cell_type_landmark_counts.tsv');rec['assigned_anchor_landmarks']=int(lm[lm.cell_type==typ].n_landmarks.iloc[0]);groups.append(rec)
        current=ct[(ct.atlas==atlas)&(ct.subset_size==f)]
        for key,x in values.items():
            means=pd.Series(x).groupby(labels).median().reindex(current.cell_type).to_numpy();y=current.log2_enrichment.to_numpy()
            assoc.append(dict(atlas=atlas,subset_size=f,metric=key,n_types=len(y),spearman_rho=spearmanr(means,y).statistic))
    print('Summarized saved per-cell diagnostics:',atlas,flush=True)
put('cached_cell_metric_selection_summary',pd.DataFrame(cells));put('cached_focus_geometry',pd.DataFrame(groups));put('cached_geometry_associations',pd.DataFrame(assoc))
print(pd.DataFrame(cells).to_string(index=False),flush=True);print(pd.DataFrame(groups).to_string(index=False),flush=True)
d=pd.read_csv(dest/'all_saved_metrics.tsv',sep='\t')
print('Effective epsilon by method:',d[d.cohort=='all_queries'].groupby(['method','effective_epsilon']).size().to_string(),flush=True)
print('Min calibration rate',d.fraction_within_one_percent.min(),'max sampled relative error',d.relative_error_max.max(),flush=True)
q=pd.read_csv(dest/'benchmark_comparisons.tsv',sep='\t');q=q[q.cohort=='nonselected_queries']
print('Nonselected q90: proposed > worst random?',(q.q90_distance_proposed>q.q90_distance_random_max).sum(),'of',len(q),flush=True)
print('Nonselected fixed radius: proposed < worst random?',(q.coverage_within_radius_proposed<q.coverage_within_radius_random_min).sum(),'of',len(q),flush=True)
print('Radii:',q.groupby('atlas').radius.first().to_string(),flush=True)
for metric in ['TV','JS_bits','median_distance','q90_distance','q99_distance']:
    a=pd.read_csv(dest/'benchmark_comparisons.tsv',sep='\t');a=a[a.cohort=='all_queries'];print(metric,'proposed > all 20 random in',int((a[metric+'_proposed']>a[metric+'_random_max']).sum()),'of15',flush=True)
print('Counts and timing:',load('S6_DA_cost_reference').to_string(index=False),flush=True)
