#!/bin/bash
#SBATCH --job-name=scatac_audit
#SBATCH --account=cedar
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=publication_figures/logs/scientific_audit_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python publication_figures/scripts/audit_saved_benchmarks.py
