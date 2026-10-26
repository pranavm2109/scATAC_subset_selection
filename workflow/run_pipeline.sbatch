#!/bin/bash
#SBATCH --job-name=scatac_pipeline_controller
#SBATCH --account=cedar
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/pipeline_controller_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
set +u
source /home/exacloud/gscratch/CEDAR/mishrap/miniconda/etc/profile.d/conda.sh
conda activate snakemake
set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Controller runs on a compute node; individual rules request their own SLURM resources.
# Optional extra arguments, e.g. --dry-run, are forwarded to Snakemake.
snakemake --snakefile workflow/Snakefile --executor slurm \
    --jobs 4 --local-cores 1 --default-resources slurm_account=cedar \
    --rerun-incomplete --keep-going --latency-wait 120 --printshellcmds \
    "$@" all
