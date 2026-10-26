from common import *
import re,os
assert os.environ.get('SLURM_JOB_ID')
rows=[]
for name in ['S1_normalization','S2_population_representation','S3_DA_fidelity','S3b_direction_specific_recovery','S5_computational_cost','S7_descriptive_UMAP_grid','S8_DA_effect_density_grid']:
    folder=OUT
    raw=(folder/'figures'/f'{name}.pdf').read_bytes()
    m=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',raw);assert m,name
    width=(float(m[3])-float(m[1]))/72*25.4
    assert abs(width-178)<.01,(name,width)
    assert b'/FontFile2' in raw and b'/Subtype /Type3' not in raw,name
    assert b'/FontDescriptor' in raw,name
    rows.append(dict(figure=name,width_mm=width,embedded_TrueType=True,type3_fonts=False))
table('PDF_export_QA',pd.DataFrame(rows))
print('PASS: all seven PDFs are 178 mm wide with embedded TrueType font programs and no Type 3 fonts')
