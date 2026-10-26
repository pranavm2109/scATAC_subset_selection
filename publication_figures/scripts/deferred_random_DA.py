"""Approval-gated DA on pre-existing random memberships; no automatic sampling."""
from pathlib import Path
import argparse,os,subprocess,json
from common import ROOT,OUT,ATLASES,FRACS,read

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--atlas',choices=ATLASES,required=True);ap.add_argument('--pct',type=int,choices=[2,4,6,8,10],required=True);ap.add_argument('--seed',type=int,default=2001);ap.add_argument('--memberships-dir',default=str(OUT/'benchmarks'/'random_DA_stage_A_2001'/'memberships'));ap.add_argument('--outdir',default=str(OUT/'benchmarks'/'random_DA_stage_A_2001'/'runs'));ap.add_argument('--execute',action='store_true');args=ap.parse_args()
    f=f'{args.pct:02d}pct';ids=Path(args.memberships_dir)/f'{args.atlas}_evaluate_random_{f}_seed{args.seed}_selected_cells.txt';dest=Path(args.outdir)/args.atlas/f'seed{args.seed}'/f
    cmd=['/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/r_env2/bin/Rscript',str(ROOT/'scripts/09_1_run_signac_da_all_cell_types.R'),str(ROOT/f'data/processed/{args.atlas}/{args.atlas}_barebones_seurat.rds'),str(ids),args.atlas,f,str(dest/'signac_da_peaks.tsv'),str(dest/'signac_da_status.tsv'),str(dest/'signac_da_reference_effects.tsv'),str(ROOT/f'results/signac_da/{args.atlas}/100pct/signac_da_peaks.tsv'),'10','20','0.05','0.25','0.01']
    if not args.execute:print(json.dumps(dict(status='PLAN ONLY; no DA executed',command=cmd,serialization='Never overlap jobs reading the same atlas BPCells object'),indent=2));return
    if 'SLURM_JOB_ID' not in os.environ:ap.error('Use a SLURM compute node')
    if dest.exists():raise FileExistsError('Refusing to overwrite previous random DA outputs')
    if not ids.exists():raise FileNotFoundError('Generate approved random memberships before DA')
    expected=read('results/publication_panels/peakvi_subset_overlay_3x5_summary.tsv');expected=int(expected[(expected.atlas==args.atlas)&(expected.fraction==f)].iloc[0].n_selected_cells)
    values=ids.read_text().splitlines();assert len(values)==len(set(values))==expected
    import hashlib,time
    expected_hash=read(f'publication_figures/benchmarks/full/{args.atlas}_evaluate_random_{f}_seed{args.seed}_metrics.tsv').iloc[0].selected_ID_sha256
    actual_hash=hashlib.sha256('\n'.join(values).encode()).hexdigest()
    if actual_hash!=expected_hash:raise ValueError('Membership hash differs from completed random benchmark')
    validation=Path(args.memberships_dir).parent/'preflight'/'seurat_validation.tsv'
    checked=__import__('pandas').read_csv(validation,sep='\t')
    assert len(checked[(checked.atlas==args.atlas)&(checked.subset_size==f)&(checked.n_selected==expected)])==1
    dest.mkdir(parents=True);env=os.environ.copy();env['SIGNAC_DA_TMPDIR']=str(dest/'scratch');env.update(OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4')
    record=dict(command=cmd,slurm_job_id=os.environ['SLURM_JOB_ID'],array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'),membership_sha256=actual_hash,production_script_sha256=hashlib.sha256((ROOT/'scripts/09_1_run_signac_da_all_cell_types.R').read_bytes()).hexdigest(),full_reference=str(ROOT/f'results/signac_da/{args.atlas}/100pct/signac_da_peaks.tsv'),started_unix=time.time())
    (dest/'execution.json').write_text(json.dumps(record,indent=2))
    with (dest/'DA.log').open('w') as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
    record.update(returncode=result.returncode,finished_unix=time.time());(dest/'execution.json').write_text(json.dumps(record,indent=2))
    result.check_returncode()
if __name__=='__main__':main()
