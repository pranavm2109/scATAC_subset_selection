#!/bin/bash
#SBATCH --job-name=scatac_DA_compare
#SBATCH --account=cedar
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000M
#SBATCH --time=01:00:00
#SBATCH --output=publication_figures/logs/DA_compare_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/scripts/compare_stage_A.py
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/scripts/report_stage_A.py
