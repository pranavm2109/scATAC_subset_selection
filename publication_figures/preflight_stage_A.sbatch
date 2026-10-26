#!/bin/bash
#SBATCH --job-name=scatac_DA_preflight
#SBATCH --account=cedar
#SBATCH --cpus-per-task=4
#SBATCH --mem=128000M
#SBATCH --time=01:00:00
#SBATCH --output=publication_figures/logs/DA_preflight_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PYTHON=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python
RSCRIPT=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/r_env2/bin/Rscript
"$PYTHON" publication_figures/scripts/prepare_stage_A.py
"$RSCRIPT" publication_figures/scripts/validate_stage_A.R
BASE=publication_figures/benchmarks/random_DA_stage_A_2001/preflight
export SIGNAC_DA_TMPDIR="$PWD/$BASE/scratch"
"$RSCRIPT" scripts/09_1_run_signac_da_all_cell_types.R "$BASE/smoke_input.rds" "$BASE/smoke_ids.txt" kidney_atac 02pct "$BASE/smoke_da.tsv" "$BASE/smoke_status.tsv" "$BASE/smoke_effects.tsv" "$BASE/smoke_reference.tsv" 10 20 0.05 0.25 0.01
touch "$BASE/PASS"
