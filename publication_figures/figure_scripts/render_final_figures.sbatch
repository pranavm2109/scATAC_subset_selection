#!/bin/bash
#SBATCH --job-name=scatac_final_figures_only
#SBATCH --account=cedar
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=publication_figures/logs/final_figures_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export MPLCONFIGDIR="${TMPDIR:-/tmp}/final-figures-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/figure_scripts/plot_figure_1_final.py
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/figure_scripts/plot_s9_overview.py
