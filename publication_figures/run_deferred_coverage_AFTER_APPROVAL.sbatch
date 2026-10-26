#!/bin/bash
#SBATCH --job-name=scatac_coverage
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=32000M
#SBATCH --time=08:00:00
#SBATCH --output=publication_figures/logs/deferred_coverage_%A_%a.log
# No array default: reviewer chooses pilot or approved expansion indices explicitly.
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python - <<'PY'
import os,subprocess,sys
import pandas as pd
idx=int(os.environ['SLURM_ARRAY_TASK_ID']);r=pd.read_csv('publication_figures/tables/S4_S6_DEFERRED_task_manifest.tsv',sep='\t').set_index('task_index').loc[idx]
phase=os.environ.get('BENCHMARK_PHASE','pilot')
if phase not in ['pilot','full','composition']:raise ValueError('Invalid benchmark phase')
cmd=[sys.executable,'publication_figures/scripts/deferred_benchmarks.py','--atlas',r.atlas,'--pct',str(r.pct),'--method',r.method,'--seed',str(r.seed),'--threads','1','--outdir',f'publication_figures/benchmarks/{phase}','--execute']
if phase=='composition':cmd+=['--composition-only','--write-selected-ids']
else:
 radius=f'publication_figures/benchmarks/radius/{r.atlas}_radius_production_02pct_seed2001_radius.tsv'
 cmd+=['--query-sample','50000' if phase=='pilot' else '0','--radius-table',radius]
subprocess.run(cmd,check=True)
PY
