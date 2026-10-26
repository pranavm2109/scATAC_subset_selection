"""Plot saved scientific outputs only; never train, select cells, search neighbors or run DA."""
from pathlib import Path
import os,sys,runpy,tempfile,shutil,subprocess,re,csv
assert os.environ.get('SLURM_JOB_ID'), 'Submit render_all_publication_figures.sbatch'
ROOT=Path(__file__).resolve().parents[2];PUB=ROOT/'publication_figures'
sys.path.insert(0,str(PUB/'scripts'))
STEMS=['S1_normalization','S2_population_representation','S3_DA_fidelity','S3b_direction_specific_recovery','S4_measured_60D_coverage','S5_computational_cost','S6_measured_random_comparison','S6b_measured_rare_population_counts','S7_descriptive_UMAP_grid','S8_DA_effect_density_grid']
# Historical diagnostics and measured benchmarks must already exist. A core pipeline run
# alone does not regenerate them. Validation must fail rather than invent missing results.
metrics=list((PUB/'benchmarks/full').glob('*_metrics.tsv'))
assert len(metrics)==315, 'Need the completed 315-run saved benchmark; no benchmark will be launched'
for folder in ['main_figures','supplementary_figures','figure_tables','reports']:(PUB/folder).mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix='publication_plot_',dir=PUB/'reports') as folder:
    stage=Path(folder)
    for sub in ['figures','tables','captions']:(stage/sub).mkdir()
    import common
    common.OUT=stage  # Legacy plotting helpers write into staging, not historical files.
    runpy.run_path(str(PUB/'scripts/validate_publication_inputs.py'),run_name='__main__')
    import plot_existing_summaries as summary
    import revised_grids as grids
    for fun in [summary.s1,summary.s2,summary.s3,summary.s5,grids.umap_grid,grids.da_grid]:
        print('Rendering',fun.__name__,flush=True);fun()
    import plot_random_benchmark_when_available as benchmark
    sys.argv=['plot_random_benchmark_when_available.py','--benchmark-dir',str(PUB/'benchmarks/full')]
    benchmark.main()  # Summarize saved small benchmark files; no NN/selection/DA work.
    for stem in STEMS:
        for ext in ['pdf','svg','png']:
            src=stage/'figures'/f'{stem}.{ext}';assert src.exists(),src
            shutil.copy2(src,PUB/'supplementary_figures'/src.name)
    # Current reviewed captions remain in place; do not copy old generated caption text.
    for p in (stage/'tables').glob('*.tsv'):
        shutil.copy2(p,PUB/'figure_tables'/('render_QA_'+p.name))
# Separate processes preserve each final figure's established matplotlib settings.
for script in ['plot_figure_1_final.py','plot_s9_overview.py']:
    subprocess.run([sys.executable,str(PUB/'figure_scripts'/script)],check=True)
main=PUB/'main_figures/Figure_1_rare_population_representation_and_DA_fidelity'
s9=PUB/'supplementary_figures/S9_quantitative_selection_and_biological_fidelity_overview'
proofs=PUB/'reports/final_figure_proofs';proofs.mkdir(exist_ok=True)
from PIL import Image
import xml.etree.ElementTree as ET
rows=[]
for base in [main]+[PUB/'supplementary_figures'/s for s in STEMS]+[s9]:
    raw=base.with_suffix('.pdf').read_bytes()
    m=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',raw);assert m,base
    width=(float(m[3])-float(m[1]))/72*25.4
    assert abs(width-178)<.02 and b'/FontFile2' in raw and b'/Subtype /Type3' not in raw,base
    ET.parse(base.with_suffix('.svg'))
    with Image.open(base.with_suffix('.png')) as im:
        dpi=im.info.get('dpi',(0,0))[0];assert abs(dpi-600)<1,(base,dpi)
    subprocess.run(['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-r150','-dFirstPage=1','-dLastPage=1','-sOutputFile='+str(proofs/(base.name+'.png')),str(base.with_suffix('.pdf'))],check=True)
    rows.append(dict(figure=base.name,width_mm=width,png_dpi=dpi,embedded_TrueType=True,svg_parsed=True,pdf_proof_rendered=True))
with (PUB/'figure_tables/final_render_QA.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('PASS: Figure 1 and eleven supplementary figures exported; inspect PDF proofs before submission.',flush=True)
