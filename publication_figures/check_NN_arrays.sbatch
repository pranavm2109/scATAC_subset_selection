#!/bin/bash
#SBATCH --job-name=scatac_NN_check
#SBATCH --account=cedar
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M
#SBATCH --time=00:15:00
#SBATCH --output=publication_figures/logs/NN_array_check_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/scripts/check_NN_array_inputs.py
