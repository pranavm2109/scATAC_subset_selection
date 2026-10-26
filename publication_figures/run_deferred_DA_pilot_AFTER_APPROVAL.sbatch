#!/bin/bash
#SBATCH --job-name=scatac_random_DA
#SBATCH --account=cedar
#SBATCH --cpus-per-task=4
#SBATCH --mem=128000M
#SBATCH --time=01:00:00
#SBATCH --array=0-5%1
#SBATCH --output=publication_figures/logs/deferred_DA_%A_%a.log
# Global concurrency one serializes access to shared BPCells atlas objects.
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python - <<'PY'
import os,sys,subprocess
from pathlib import Path
import hashlib
root=Path('publication_figures/benchmarks/random_DA_stage_A_2001')
assert (root/'preflight/PASS').exists(),'Preflight and small-input test must pass first'
assert hashlib.sha256(Path('scripts/09_1_run_signac_da_all_cell_types.R').read_bytes()).hexdigest()==(root/'preflight/production_script.sha256').read_text().strip(),'Production script changed after preflight'
atlases=['kidney_atac','heca_atac','fetal_atac'];i=int(os.environ['SLURM_ARRAY_TASK_ID'])
assert 0<=i<6
subprocess.run([sys.executable,'publication_figures/scripts/deferred_random_DA.py','--atlas',atlases[i//2],'--pct',str([2,10][i%2]),'--seed','2001','--memberships-dir',str(root/'memberships'),'--outdir',str(root/'runs'),'--execute'],check=True)
PY
