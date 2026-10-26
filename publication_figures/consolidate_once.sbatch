#!/bin/bash
#SBATCH --job-name=scatac_consolidate
#SBATCH --account=cedar
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=publication_figures/logs/consolidate_%j.log
set -euo pipefail
cd /home/groups/CEDAR/mishrap/scATAC_subset_selection
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON=/home/exacloud/gscratch/CEDAR/mishrap/miniconda/envs/scvi/bin/python
"$PYTHON" publication_figures/scripts/consolidate_figures.py
bash publication_figures/regenerate_all_figures.sbatch
"$PYTHON" - <<'PY'
from pathlib import Path
import shutil,pandas as pd
pub=Path('publication_figures')
for name in ['S1_normalization','S2_population_representation','S3_DA_fidelity','S3b_direction_specific_recovery','S5_computational_cost','S7_descriptive_UMAP_grid','S8_DA_effect_density_grid']:
    for ext in ['pdf','svg','png']:assert (pub/'figures'/f'{name}.{ext}').stat().st_size>0
for name in ['polished','final_minimal']:shutil.rmtree(pub/name)
files=list((pub/'benchmarks/pilot').glob('*_metrics.tsv'))
d=pd.concat([pd.read_csv(p,sep='\t') for p in files]);d=d[d.cohort=='all_queries']
assert len(d)==12 and d.fraction_within_one_percent.ge(.99).all()
d[['atlas','subset_size','method','seed','n_query','fraction_within_one_percent','effective_epsilon','coverage_processing_seconds']].to_csv(pub/'tables/S4_S6_pilot_review.tsv',sep='\t',index=False)
print('PASS: 12 pilot accuracy gates. Pilot processing seconds:',d.coverage_processing_seconds.min(),d.coverage_processing_seconds.max())
print('Consolidation complete; duplicate directories deleted.')
PY
