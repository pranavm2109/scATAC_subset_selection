"""One-time, user-authorized consolidation; execute on a compute node."""
from pathlib import Path
import os,shutil
assert os.environ.get('SLURM_JOB_ID')
root=Path.cwd();pub=root/'publication_figures';old=pub/'polished';final=pub/'final_minimal';scripts=pub/'scripts'
style=(final/'scripts/style.py').read_text().replace('parents[3]','parents[2]').replace("OUT=ROOT/'publication_figures'/'final_minimal'","OUT=ROOT/'publication_figures'")
(scripts/'common.py').write_text(style)
source=(old/'scripts/summaries.py').read_text().replace('from style import *','from common import *')
source=source.replace("panel(c,'C','Historical 2D distances (2%)')","panel(c,'C','Nearest-selected-cell distances\\nin 2D UMAP (2%)')")
source=source.replace("100*t.n_evaluable_cell_types/t.n_reference_cell_types,'--'","100*t.n_evaluable_cell_types/t.n_reference_cell_types,'-s'")
source=source.replace("color='#555555',label='Reference pairs'","color='#555555',marker='o',label='Reference pairs'")
source=source.replace("color='#555555',ls='--',label='Cell types'","color='#555555',marker='s',ls='-',label='Cell types'")
(scripts/'plot_existing_summaries.py').write_text(source)
(scripts/'revised_grids.py').write_text((old/'scripts/grids.py').read_text().replace('from style import *','from common import *'))
validation=(old/'scripts/run.py').read_text().split("d=load('S3_per_cell_type_direction')",1)[1].split('for fun in [',1)[0]
validation="from common import *\nimport os\nassert os.environ.get('SLURM_JOB_ID')\nd=load('S3_per_cell_type_direction')"+validation
validation+='\nidentity=load("S1_historical_subset_identity")\nfor r in identity.itertuples():\n    assert set(pd.read_csv(ROOT/r.historical,header=None,dtype=str)[0])==set(pd.read_csv(ROOT/r.production,header=None,dtype=str)[0])\nprint("PASS: cohort, denominator and historical membership checks")\n'
(scripts/'validate_publication_inputs.py').write_text(validation)
(scripts/'check_pdf_exports.py').write_text((final/'scripts/pdf_qc.py').read_text().replace('from style import *','from common import *').replace("folder=OUT if name in ['S1_normalization','S3_DA_fidelity'] else ROOT/'publication_figures'/'polished'","folder=OUT"))
(pub/'captions'/'templates').mkdir(exist_ok=True)
for srcdir in [old/'captions',final/'captions']:
    for p in srcdir.glob('*.md'):
        text=p.read_text().replace('polished/tables/','tables/')
        (pub/'captions'/p.name).write_text(text)
        (pub/'captions'/'templates'/p.name).write_text(text)
for srcdir in [old/'figures',final/'figures']:
    for p in srcdir.iterdir():
        if p.is_file():shutil.copy2(p,pub/'figures'/p.name)
for srcdir in [old/'tables',final/'tables']:
    for p in srcdir.glob('*.tsv'):shutil.copy2(p,pub/'tables'/p.name)
# Preserve the review record without retaining duplicate output directories.
(pub/'FINAL_FIGURE_CHANGELOG.md').write_text((final/'CHANGELOG.md').read_text())
for name in ['common.py','plot_existing_summaries.py','revised_grids.py','validate_publication_inputs.py','check_pdf_exports.py']:
    compile((scripts/name).read_text(),name,'exec')
print('Canonical sources, captions and latest outputs transferred; regenerating before deletion.',flush=True)
