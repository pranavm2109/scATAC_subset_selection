#!/bin/bash
#SBATCH --job-name=scatac_figures
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=16000M
#SBATCH --time=02:00:00
#SBATCH --output=publication_figures/logs/existing_figures_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export MPLCONFIGDIR="${TMPDIR:-/tmp}/scatac-publication-mpl-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PYTHON=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python
bash publication_figures/regenerate_all_figures.sbatch
