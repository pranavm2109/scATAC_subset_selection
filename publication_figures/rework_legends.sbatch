#!/bin/bash
#SBATCH --job-name=scatac_legend_QA
#SBATCH --account=cedar
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M
#SBATCH --time=00:20:00
#SBATCH --output=publication_figures/logs/legend_QA_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export MPLCONFIGDIR="${TMPDIR:-/tmp}/scatac-legend-mpl-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PYTHON=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python
bash publication_figures/regenerate_all_figures.sbatch
