#!/bin/bash
#SBATCH --job-name=scatac_radius
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=32000M
#SBATCH --time=08:00:00
#SBATCH --array=0-2%1
#SBATCH --output=publication_figures/logs/deferred_radius_%A_%a.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
ATLASES=(kidney_atac heca_atac fetal_atac)
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/scripts/deferred_benchmarks.py --atlas "${ATLASES[$SLURM_ARRAY_TASK_ID]}" --mode radius --threads 1 --outdir publication_figures/benchmarks/radius --execute
