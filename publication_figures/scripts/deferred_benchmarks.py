"""S4/S6 approval-gated benchmark; NEVER called by run_existing_figures.sbatch.

Default prints a plan. --execute is required for loading/processing saved latents.
No training, Apricot fit, or DA testing. Random memberships use no annotations.
No persistent NN cache. Each atlas/budget/method/seed is an independent task.
"""
from pathlib import Path
import argparse,json,time,hashlib,os
import numpy as np
import pandas as pd
from common import ROOT,OUT,ATLASES,FRACS,composition


def select_positions(method,n,k,atlas,pct,seed,cell_ids,explicit_ids=None):
    if method=='random':
        # Budgets are independent, not nested prefixes. Atlas and budget streams differ.
        rng=np.random.default_rng(np.random.SeedSequence([seed,ATLASES.index(atlas),pct]))
        return np.sort(rng.choice(n,k,replace=False))
    path=Path(explicit_ids) if explicit_ids else ROOT/f'results/subset_selection/peakvi/{atlas}/selected_cells_peakvi_{pct:02d}pct.txt'
    ids=pd.read_csv(path,header=None,dtype=str)[0]
    if ids.isna().any() or ids.duplicated().any():raise ValueError('Invalid selected IDs')
    pos=pd.Index(cell_ids).get_indexer(ids)
    if len(pos)!=k or np.any(pos<0):raise ValueError('Budget or ID alignment mismatch')
    return np.sort(pos)


def exact_minima(query,reference,exclude_reference_positions=None,block=4096):
    from scipy.spatial.distance import cdist
    mins=np.full(len(query),np.inf)
    for start in range(0,len(reference),block):
        dist=cdist(query,reference[start:start+block],metric='euclidean')
        if exclude_reference_positions is not None:
            local=np.asarray(exclude_reference_positions)-start;ok=(local>=0)&(local<dist.shape[1]);dist[np.flatnonzero(ok),local[ok]]=np.inf
        mins=np.minimum(mins,dist.min(axis=1))
    return mins


def nearest(index,reference,query,k=1,epsilon=.2,exclude_positions=None):
    # Numba/PyNNDescent kernels require writable C-contiguous float32 input.
    query_buffer=np.array(query,dtype=np.float32,order='C',copy=True)
    idx,_=index.query(query_buffer,k=k,epsilon=epsilon)
    if np.any(idx<0):raise RuntimeError('ANN query did not return valid references')
    distance=np.linalg.norm(query.astype(np.float64)[:,None,:]-reference[idx].astype(np.float64),axis=2)
    if exclude_positions is not None:distance[idx==np.asarray(exclude_positions)[:,None]]=np.inf
    result=distance.min(axis=1)
    if not np.isfinite(result).all():raise RuntimeError('No nonself candidate')
    return result


def validate_index(index,reference,full,chosen,seed,epsilon,radius_mode=False):
    rng=np.random.default_rng(seed)
    if radius_mode:
        qpos=rng.choice(len(full),min(128,len(full)),replace=False)
        exact=exact_minima(full[qpos],reference,exclude_reference_positions=qpos)
    else:
        mask=np.ones(len(full),dtype=bool);mask[chosen]=False;pool=np.flatnonzero(mask);qpos=rng.choice(pool,min(128,len(pool)),replace=False)
        exact=exact_minima(full[qpos],reference)
    attempts=[]
    for effort in sorted(set([epsilon]+[v for v in [.5,1.,2.] if v>epsilon])):
        approx=nearest(index,reference,full[qpos],k=4 if radius_mode else 1,epsilon=effort,exclude_positions=qpos if radius_mode else None)
        rel=(approx-exact)/np.maximum(exact,1e-12);within=float(np.mean(approx<=exact*1.01+1e-7));attempts.append(f'{effort}:{within:.6f}')
        stats=dict(n_exact_validation=len(exact),fraction_within_one_percent=within,relative_error_q99=float(np.quantile(rel,.99)),relative_error_max=float(rel.max()),effective_epsilon=effort,accuracy_attempts=';'.join(attempts))
        if within>=.99:return stats
    raise RuntimeError(f'ANN accuracy gate failed even after increased search effort: {stats}; increase index quality before reporting.')


def radius_calibration(args,latent,prefix):
    from pynndescent import NNDescent
    # A view of mmap_mode='r' stays read-only; allocate a private writable buffer.
    reference=np.array(latent,dtype=np.float32,order='C',copy=True)
    assert reference.flags.writeable and reference.flags.c_contiguous
    started=time.time();index=NNDescent(reference,metric='euclidean',n_neighbors=30,random_state=1701,n_jobs=args.threads,low_memory=True);index.prepare()
    validation=validate_index(index,latent,latent,None,1702,args.epsilon,True)
    epsilon=validation['effective_epsilon']
    q=np.random.default_rng(1703).choice(len(latent),min(50000,len(latent)),replace=False)
    distances=[]
    for start in range(0,len(q),args.chunk):
        ii=q[start:start+args.chunk];distances.append(nearest(index,latent,latent[ii],k=4,epsilon=epsilon,exclude_positions=ii))
    distances=np.concatenate(distances);r=float(np.quantile(distances,.90))
    if not r>0:raise RuntimeError('Nonpositive radius: duplicates require separate scientific review')
    record=dict(atlas=args.atlas,radius=r,definition='90th percentile of nonself nearest-neighbor distance in full-atlas raw 60D PeakVI; fixed sampled query set',n_radius_queries=len(q),query_seed=1703,index_seed=1701,epsilon=epsilon,threads=args.threads,seconds=time.time()-started,**validation)
    pd.DataFrame([record]).to_csv(str(prefix)+'_radius.tsv',sep='\t',index=False)


def evaluate(args,latent,cell_ids,prefix):
    from pynndescent import NNDescent
    p=pd.read_csv(ROOT/'results/publication_panels/peakvi_subset_overlay_3x5_summary.tsv',sep='\t');f=f'{args.pct:02d}pct';row=p[(p.atlas==args.atlas)&(p.fraction==f)].iloc[0]
    k=int(row.n_selected_cells);n=len(latent);assert n==row.n_atlas_cells
    chosen=select_positions(args.method,n,k,args.atlas,args.pct,args.seed,cell_ids,args.selected_ids)
    selected_hash=hashlib.sha256('\n'.join(cell_ids[chosen]).encode()).hexdigest()
    # Annotations are loaded only AFTER membership has been fixed.
    meta=pd.read_csv(ROOT/f'results/cell_metadata/{args.atlas}/{args.atlas}_cell_metadata.tsv',sep='\t',usecols=['cell_id','cell_type'],dtype=str).set_index('cell_id')
    if meta.index.has_duplicates:raise ValueError('Duplicated metadata IDs')
    labels=meta.loc[cell_ids,'cell_type'];codes,types=pd.factorize(labels,sort=True)
    if (codes<0).any():raise ValueError('Missing annotations: define treatment before benchmarking')
    fullcounts=np.bincount(codes,minlength=len(types));subcounts=np.bincount(codes[chosen],minlength=len(types));tv,js=composition(fullcounts/n,subcounts/k)
    population=pd.DataFrame(dict(cell_type=types,full_count=fullcounts,full_proportion=fullcounts/n,subset_count=subcounts,subset_proportion=subcounts/k,enrichment=(subcounts/k)/(fullcounts/n),rare_full_below_0p1pct=fullcounts/n<.001,DA_eligible=(subcounts>=10)&(k-subcounts>=20)))
    population=population.assign(atlas=args.atlas,subset_size=f,method=args.method,seed=args.seed if args.method=='random' else 123)
    population.to_csv(str(prefix)+'_populations.tsv',sep='\t',index=False)
    if args.write_selected_ids:pd.Series(cell_ids[chosen]).to_csv(str(prefix)+'_selected_cells.txt',index=False,header=False)
    if args.composition_only:
        pd.DataFrame([dict(atlas=args.atlas,subset_size=f,method=args.method,seed=args.seed if args.method=='random' else 123,n_atlas=n,n_selected=k,n_retained=int((subcounts>0).sum()),n_types=len(types),TV=tv,JS_bits=js,selected_ID_sha256=selected_hash)]).to_csv(str(prefix)+'_composition_metrics.tsv',sep='\t',index=False)
        return
    query=np.arange(n) if args.query_sample==0 else np.sort(np.random.default_rng(1704).choice(n,min(args.query_sample,n),replace=False))
    # Same query stream for every budget, method and seed within an atlas.
    reference=np.array(latent[chosen],dtype=np.float32,order='C',copy=True);started=time.time()
    index=NNDescent(reference,metric='euclidean',n_neighbors=30,random_state=1705,n_jobs=args.threads,low_memory=True);index.prepare()
    validation=validate_index(index,reference,latent,chosen,1706,args.epsilon)
    epsilon=validation['effective_epsilon']
    selected_mask=np.zeros(n,dtype=bool);selected_mask[chosen]=True;distances=np.empty(len(query),dtype=float)
    for start in range(0,len(query),args.chunk):
        ii=query[start:start+args.chunk];distances[start:start+len(ii)]=nearest(index,reference,latent[ii],epsilon=epsilon)
    distances[selected_mask[query]]=0.0
    radius=None
    if args.radius_table:
        rr=pd.read_csv(args.radius_table,sep='\t');rr=rr[rr.atlas==args.atlas]
        if len(rr)!=1:raise ValueError('Expected one fixed radius per atlas')
        radius=float(rr.iloc[0].radius)
    rows=[]
    for cohort,keep in [('all_queries',np.ones(len(query),dtype=bool)),('nonselected_queries',~selected_mask[query])]:
        dd=distances[keep]
        rows.append(dict(atlas=args.atlas,subset_size=f,method=args.method,seed=args.seed if args.method=='random' else 123,n_atlas=n,n_selected=k,n_query=len(dd),query_mode='full_atlas' if args.query_sample==0 else 'fixed_uniform_query_sample',query_seed=1704,cohort=cohort,median_distance=float(np.median(dd)),q90_distance=float(np.quantile(dd,.90)),q99_distance=float(np.quantile(dd,.99)),mean_distance=float(dd.mean()),fixed_radius=radius,coverage_within_radius=float(np.mean(dd<=radius)) if radius is not None else np.nan,TV=tv,JS_bits=js,n_retained=int((subcounts>0).sum()),n_types=len(types),n_DA_eligible=int(((subcounts>=10)&(k-subcounts>=20)).sum()),selected_ID_sha256=selected_hash,ann_index_seed=1705,epsilon=epsilon,threads=args.threads,coverage_processing_seconds=time.time()-started,**validation))
    pd.DataFrame(rows).to_csv(str(prefix)+'_metrics.tsv',sep='\t',index=False)
    from importlib.metadata import version
    (Path(str(prefix)+'_provenance.json')).write_text(json.dumps(dict(arguments=vars(args),latent_path=str(ROOT/f'results/peakvi/{args.atlas}/peakvi_latent.npy'),latent_shape=list(latent.shape),latent_bytes=int(latent.nbytes),selected_ID_sha256=selected_hash,coordinate_transform='none; saved raw 60D Euclidean PeakVI',versions={x:version(x) for x in ['numpy','scipy','pynndescent','pandas']},interpretation='Approximate nearest distances recomputed exactly to returned candidates; upper bounds on true nearest distance. Radius coverage is conservative conditional on the fixed radius, which is itself ANN-estimated.'),indent=2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--atlas',choices=ATLASES,required=True);ap.add_argument('--mode',choices=['evaluate','radius'],default='evaluate');ap.add_argument('--pct',type=int,choices=[2,4,6,8,10],default=2);ap.add_argument('--method',choices=['production','random','historical'],default='production');ap.add_argument('--seed',type=int,default=2001);ap.add_argument('--selected-ids');ap.add_argument('--radius-table');ap.add_argument('--query-sample',type=int,default=0,help='0=every atlas cell; otherwise shared fixed uniform query sample');ap.add_argument('--threads',type=int,default=1,help='1 favors reproducibility; parallel ANN variability must be recorded');ap.add_argument('--chunk',type=int,default=10000);ap.add_argument('--epsilon',type=float,default=.2);ap.add_argument('--outdir',default=str(OUT/'benchmarks'));ap.add_argument('--composition-only',action='store_true',help='Only inexpensive count metrics after membership generation; no ANN index');ap.add_argument('--write-selected-ids',action='store_true');ap.add_argument('--execute',action='store_true');args=ap.parse_args()
    if args.method=='historical' and not args.selected_ids:ap.error('Historical unnormalized IDs must be explicitly supplied and provenance verified; production IDs now represent L1.')
    if args.query_sample<0 or args.threads<1 or args.chunk<1 or args.epsilon<=0:ap.error('Invalid query/compute parameters')
    if not args.execute:
        print(json.dumps(dict(status='PLAN ONLY; no latent loaded and no sampling executed',arguments=vars(args),proposed_seeds=list(range(2001,2021))),indent=2));return
    if 'SLURM_JOB_ID' not in os.environ:ap.error('Execute on a SLURM compute node, not a login node')
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);tag=f'{args.atlas}_{args.mode}_{args.method}_{args.pct:02d}pct_seed{args.seed}';prefix=out/tag
    if list(out.glob(tag+'_*')):raise FileExistsError('Refusing to overwrite existing benchmark artifacts; choose a new --outdir')
    latent=np.load(ROOT/f'results/peakvi/{args.atlas}/peakvi_latent.npy',mmap_mode='r');cell_ids=pd.read_csv(ROOT/f'results/peakvi/{args.atlas}/cell_names.txt',header=None,dtype=str)[0].to_numpy()
    assert latent.shape==(len(cell_ids),60) and len(set(cell_ids))==len(cell_ids)
    if args.mode=='radius':radius_calibration(args,latent,prefix)
    else:evaluate(args,latent,cell_ids,prefix)

if __name__=='__main__':main()
