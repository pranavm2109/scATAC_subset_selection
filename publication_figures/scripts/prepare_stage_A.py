from common import *
import os,hashlib,json
assert os.environ.get('SLURM_JOB_ID')
base=OUT/'benchmarks/random_DA_stage_A_2001';members=base/'memberships';members.mkdir(parents=True,exist_ok=True);(base/'preflight').mkdir(exist_ok=True)
rows=[]
for a in ATLASES:
    names=pd.read_csv(ROOT/f'results/peakvi/{a}/cell_names.txt',header=None,dtype=str)[0].to_numpy()
    for pct in [2,10]:
        f=f'{pct:02d}pct';m=read(f'publication_figures/benchmarks/full/{a}_evaluate_random_{f}_seed2001_metrics.tsv');r=m[m.cohort=='all_queries'].iloc[0]
        chosen=np.sort(np.random.default_rng(np.random.SeedSequence([2001,ATLASES.index(a),pct])).choice(len(names),int(r.n_selected),replace=False))
        ids=names[chosen];assert len(ids)==len(set(ids))==r.n_selected
        token=hashlib.sha256('\n'.join(ids).encode()).hexdigest();assert token==r.selected_ID_sha256
        path=members/f'{a}_evaluate_random_{f}_seed2001_selected_cells.txt'
        content='\n'.join(ids)+'\n'
        if path.exists():assert path.read_text()==content
        else:path.write_text(content)
        rows.append(dict(atlas=a,subset_size=f,n_selected=len(ids),seed=2001,sha256=token,ids_file=str(path)))
pd.DataFrame(rows).to_csv(base/'preflight/membership_manifest.tsv',sep='\t',index=False)
source=ROOT/'scripts/09_1_run_signac_da_all_cell_types.R'
(base/'preflight/production_script.sha256').write_text(hashlib.sha256(source.read_bytes()).hexdigest()+'\n')
print('Verified and wrote six exact existing random memberships',flush=True)
