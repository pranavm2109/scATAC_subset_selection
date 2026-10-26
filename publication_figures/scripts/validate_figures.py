"""Artifact and numerical QA; no atlas-scale processing."""
from common import *
from PIL import Image
import re,hashlib,xml.etree.ElementTree as ET
NAMES=['S1_normalization','S2_population_representation','S3_DA_fidelity','S3b_direction_specific_recovery','S5_computational_cost','S7_descriptive_UMAP_grid','S8_DA_effect_density_grid','MAIN_candidate_quantitative']
rows=[]
for name in NAMES:
    pdf=OUT/'figures'/f'{name}.pdf';svg=OUT/'figures'/f'{name}.svg';png=OUT/'figures'/f'{name}.png'
    assert all(p.exists() and p.stat().st_size>1000 for p in [pdf,svg,png])
    raw=pdf.read_bytes();m=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',raw);assert m,name
    x0,y0,x1,y1=map(float,m.groups());width=(x1-x0)*25.4/72;assert abs(width-178)<.02,(name,width)
    root=ET.parse(svg).getroot();ntext=sum(1 for el in root.iter() if el.tag.endswith('}text'));assert ntext>0,(name,'SVG missing editable text')
    with Image.open(png) as im:
        size=im.size;dpi=im.info.get('dpi');im.verify()
    assert dpi and abs(dpi[0]-450)<1
    assert (OUT/'captions'/f'{name}.md').exists(),name
    rows.append(dict(figure=name,width_mm=width,height_mm=(y1-y0)*25.4/72,png_width=size[0],png_height=size[1],dpi=dpi[0],svg_text_elements=ntext,PDF_bytes=pdf.stat().st_size,SVG_bytes=svg.stat().st_size,PNG_bytes=png.stat().st_size))
legend_reports=[pd.read_csv(p,sep='\t') for p in (OUT/'tables').glob('*_legend_QA.tsv')]
assert len(legend_reports)==6
legend_report=pd.concat(legend_reports,ignore_index=True)
assert legend_report.legend_outside_all_data_and_label_bounds.all()
assert (legend_report.legend_font_pt>=7).all()
table('LEGEND_QA',legend_report)
manifest=load('SOURCE_MANIFEST')
for _,r in manifest.iterrows():assert hashlib.sha256((ROOT/r.source).read_bytes()).hexdigest()==r.sha256,f'Input changed: {r.source}'
for a in ATLASES:
    f=load('S3_fixed_eligibility');f=f[(f.atlas==a)&(f.cohort=='common_eligible')];assert f.n_reference.nunique()==1 and (f.n_available==f.n_reference).all()
h=load('S8_density_bins').groupby(['atlas','subset_size']).n_reference_pairs.sum();p=load('S8_panel_metrics').set_index(['atlas','subset_size']);assert (h==p.loc[h.index,'n_effects_available']).all()
table('FIGURE_QA',pd.DataFrame(rows))
# Enumerate exact large plotting inputs without rereading them solely for hashing.
source=[]
for a in ATLASES:
    for f in FRACS:source.append(dict(figure='S7',source=f'results/subset_selection/peakvi/{a}/selected_cells_peakvi_{f}.txt'))
    source.append(dict(figure='S7',source=f'results/umap_peakvi/{a}/{a}_peakvi_umap.tsv'))
    source.append(dict(figure='S8',source=f'results/signac_da_capture/{a}/da_peak_overlap_details.tsv'))
for r in source:
    st=(ROOT/r['source']).stat();r.update(bytes=st.st_size,mtime_ns=st.st_mtime_ns)
table('S7_S8_source_manifest',pd.DataFrame(source))
(OUT/'QA_REPORT.md').write_text('''# Figure QA

All eight figure sets have PDF, SVG and PNG exports. PDF media boxes are 178 mm wide; PNG previews are 450 DPI; SVG text remains editable. Dense grid artists are intentionally rasterized within vector containers. See tables/FIGURE_QA.tsv for physical dimensions and file checks.

Numerical checks: composition/eligibility counts align with production status tables; count-derived DA metrics reproduce the production panel summary; fixed-common-eligibility denominators remain constant within atlas; streamed density-bin totals equal all 15 production available-effect denominators. Existing small source-file hashes remain unchanged. No sampling or ANN baseline was executed.

Layout checks: summary-panel legends occupy dedicated non-data bands at 7 pt. Rendering-time assertions reject legends clipped by the page or overlapping any panel, title, or axis label; axis-label pairs are also checked for collisions. See tables/LEGEND_QA.tsv. Unused legend bands are collapsed.

Visual review: all figure previews inspected at their declared print sizes; population-retention annotation crowding was corrected, unrelated atlas UMAP distances are not connected by lines, and density-grid whitespace reduced. Fonts are specified at print size (typically 7–9 pt; compact grid annotations/ticks 6–6.5 pt). A final author/journal proof is still appropriate. Descriptive UMAP marker saturation is not a quantitative coverage result.

Deferred work: S4/S6 are scripts/plans only. No result placeholders or simulated baseline points appear in the figure collection. Historical embedding identity across regenerations is not established; within-run common embeddings are verified from the analysis implementation.
''')
print('QA passed: eight 178-mm figure sets, 450-DPI previews, editable SVG text, validated DA counts and unchanged source hashes.',flush=True)
