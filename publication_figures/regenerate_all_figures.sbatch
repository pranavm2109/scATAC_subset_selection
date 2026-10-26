#!/bin/bash
#SBATCH --job-name=scatac_figures
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=publication_figures/logs/regenerate_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export MPLCONFIGDIR="${TMPDIR:-/tmp}/scatac-figures-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PYTHON=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python
"$PYTHON" publication_figures/scripts/validate_publication_inputs.py
"$PYTHON" publication_figures/scripts/regenerate_publication.py
"$PYTHON" publication_figures/scripts/check_pdf_exports.py
mkdir -p publication_figures/figures/proofs
for name in S1_normalization S2_population_representation S3_DA_fidelity S3b_direction_specific_recovery S5_computational_cost S7_descriptive_UMAP_grid S8_DA_effect_density_grid; do
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -dFirstPage=1 -dLastPage=1 -sOutputFile="publication_figures/figures/proofs/${name}.png" "publication_figures/figures/${name}.pdf"
done
# Render S4/S6 only when their complete measured benchmark exists.
shopt -s nullglob
benchmark_files=(publication_figures/benchmarks/full/*_metrics.tsv)
if [[ ${#benchmark_files[@]} -eq 315 ]]; then
    "$PYTHON" publication_figures/scripts/plot_random_benchmark_when_available.py --benchmark-dir publication_figures/benchmarks/full
else
    echo 'S4/S6 deferred: complete full benchmark not yet available. No analyses launched.'
fi
