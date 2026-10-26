#!/bin/bash
#SBATCH --job-name=scatac_all_publication_figures
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=publication_figures/logs/all_publication_figures_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export MPLCONFIGDIR="${TMPDIR:-/tmp}/all-publication-figures-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/figure_scripts/render_all_publication_figures.py
