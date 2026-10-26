from common import *
import os
import plot_existing_summaries as summaries
import revised_grids as grids
assert os.environ.get('SLURM_JOB_ID'),'Submit the batch script; do not run on login node'
for fun in [summaries.s1,summaries.s2,summaries.s3,summaries.s5,grids.umap_grid,grids.da_grid]:
    print('Generating',fun.__name__,flush=True);fun()
# Keep final reviewed captions, including provenance qualifications.
for p in (OUT/'captions'/'templates').glob('*.md'):(OUT/'captions'/p.name).write_text(p.read_text())
print('Generated seven supplementary figures in publication_figures/figures',flush=True)
