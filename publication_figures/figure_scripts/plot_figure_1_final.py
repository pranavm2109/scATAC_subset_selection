"""Figure preparation only: read existing summary tables; never run analyses."""
from pathlib import Path
import os,json,hashlib,textwrap,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

assert os.environ.get('SLURM_JOB_ID'),'Use the authorized plotting-only compute job'
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'publication_figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.6,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
FOCUS=[('heca_atac','Haematopoietic stem and progenitor cell'),('fetal_atac','Megakaryocytes?'),('heca_atac','Naive CD4 T cell'),('fetal_atac','Ganglion cells')]
COLORS={'kidney_atac':'#0072B2','heca_atac':'#D55E00','fetal_atac':'#009E73'};LABELS={'heca_atac':'HECA','fetal_atac':'Fetal'}
SOURCES={'representation':'publication_figures/tables/S6b_measured_rare_population_counts.tsv','eligibility':'publication_figures/reports/random_DA_stage_A_2001/eligibility.tsv','DA':'publication_figures/reports/random_DA_stage_A_2001/rare_population_metrics.tsv','comparison_QA':'publication_figures/reports/random_DA_stage_A_2001/definition_QA.tsv'}
frozen=list((ROOT/'publication_figures/figures').glob('S*.*'))+list((ROOT/'publication_figures/captions').glob('S*.md'))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as fh:
        for b in iter(lambda:fh.read(1048576),b''):h.update(b)
    return h.hexdigest()
before={str(p.relative_to(ROOT)):sha(p) for p in frozen if p.is_file()}
tables={k:pd.read_csv(ROOT/v,sep='\t') for k,v in SOURCES.items()}
representation=tables['representation'];elig=tables['eligibility'];effects=tables['DA']
plotdata=[];points=[]
for a,typ in FOCUS:
    r=representation[(representation.atlas==a)&(representation.cell_type==typ)&(representation.subset_size=='02pct')]
    random=r[r.method=='random'].sort_values('seed');proposed=r[r.method=='production']
    assert len(random)==20 and set(random.seed)==set(range(2001,2021)) and len(proposed)==1
    p=proposed.iloc[0];pilot=random[random.seed==2001].iloc[0]
    assert r.full_count.nunique()==1 and r.full_proportion.nunique()==1
    budget=int(elig[(elig.atlas==a)&(elig.subset_size=='02pct')&(elig.method=='production')].n_dataset_cells.iloc[0])
    assert np.array_equal(r.DA_eligible,(r.subset_count>=10)&(budget-r.subset_count>=20))
    for method in ['production','random']:
        status=elig[(elig.atlas==a)&(elig.cell_type==typ)&(elig.subset_size=='02pct')&(elig.method==method)]
        da=effects[(effects.atlas==a)&(effects.cell_type==typ)&(effects.subset_size=='02pct')&(effects.method==method)&(effects.cohort=='original_reference')&(effects.direction=='all')]
        assert len(status)==len(da)==1
        s=status.iloc[0];d=da.iloc[0];expected=p if method=='production' else pilot
        iseligible=bool(s.n_target_cells>=10 and s.n_rest_cells>=20)
        assert s.n_target_cells==expected.subset_count and iseligible==bool(expected.DA_eligible)
        assert (s.status=='success')==iseligible
        assert d.n_strict<=d.n_available<=d.n_reference and np.isclose(d.strict_recovery,d.n_strict/d.n_reference)
        if not iseligible:assert d.n_available==0 and d.n_strict==0
        plotdata.append(dict(atlas=a,cell_type=typ,budget='02pct',method=method,DA_seed=123 if method=='production' else 2001,full_population_count=int(p.full_count),full_population_percent=100*p.full_proportion,selected_count=int(s.n_target_cells),rest_count=int(s.n_rest_cells),proposed_enrichment=float(p.enrichment),random_eligible_seeds=int(random.DA_eligible.sum()),random_total_seeds=20,DA_eligible=iseligible,n_reference=int(d.n_reference),n_available=int(d.n_available),n_strict=int(d.n_strict),strict_recovery_source=float(d.strict_recovery),strict_recovery_display_percent=100*d.strict_recovery if iseligible else np.nan,display_status='measured' if iseligible else 'DA-ineligible; not a measured zero'))
    points.extend(random.to_dict('records'))
data=pd.DataFrame(plotdata);pd.DataFrame(points).to_csv(OUT/'figure_tables/Figure_1_sources/panel_A_random_seed_counts.tsv',sep='\t',index=False)
data.to_csv(OUT/'figure_tables/Figure_1_sources/figure1_panel_values.tsv',sep='\t',index=False)
pd.DataFrame([dict(role=k,path=v,sha256=sha(ROOT/v)) for k,v in SOURCES.items()]).to_csv(OUT/'figure_tables/Figure_1_sources/source_manifest.tsv',sep='\t',index=False)

fig=plt.figure(figsize=(178/25.4,142/25.4))
labelax=fig.add_axes([.018,.19,.225,.62]);axA=fig.add_axes([.278,.19,.217,.62]);axB=fig.add_axes([.545,.19,.17,.62]);axC=fig.add_axes([.77,.19,.21,.62])
for ax in [labelax,axA,axB,axC]:ax.set_ylim(-.55,3.65);ax.set_yticks([])
labelax.set_axis_off();axB.set_axis_off()
for ax in [axA,axC]:
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='x',length=3,width=.6)
axA.set_xlim(-.8,23);axA.set_xticks([0,10,20]);axA.set_xlabel('Selected cells')
axA.axvline(10,color='#A0A6AD',ls=':',lw=.8,zorder=0)
axC.set_xlim(0,100);axC.set_xticks([0,50,100]);axC.set_xlabel('Strict recovery (%)')
axB.set_xlim(0,1)
for ax,letter,title,sub in [(axA,'A','Rare-population counts','20 random seeds; 2% budget'),(axB,'B','DA eligibility','≥10 target; ≥20 rest'),(axC,'C','DA reference recovery','Random: seed 2001 only')]:
    ax.text(-.02,1.16,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,va='bottom')
    ax.text(0,1.10,title,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
    ax.text(0,1.035,sub,transform=ax.transAxes,fontsize=6.8,va='bottom')
labelax.text(0,1.10,'Original annotation',transform=labelax.transAxes,fontweight='bold',fontsize=8)
axB.text(.23,3.44,'Proposed',ha='center',fontsize=7)
axB.text(.77,3.44,'Random\n(20 seeds)',ha='center',va='center',fontsize=7)

for i,(a,typ) in enumerate(FOCUS):
    y=3-i;color=COLORS[a];r=representation[(representation.atlas==a)&(representation.cell_type==typ)&(representation.subset_size=='02pct')]
    rand=r[r.method=='random'].sort_values('seed');p=data[(data.atlas==a)&(data.cell_type==typ)&(data.method=='production')].iloc[0]
    labelax.text(0,y+.32,LABELS[a],color=color,fontweight='bold',fontsize=7.5,va='center')
    labelax.text(0,y+.02,textwrap.fill(typ,width=25),fontsize=7.5,va='center',linespacing=1.2)
    labelax.text(0,y-.32,f'Full atlas: n = {p.full_population_count:,}',fontsize=6.8,color='#555555',va='center')
    vals=rand.subset_count.to_numpy();baseline=y-.14
    q1,med,q3=np.quantile(vals,[.25,.5,.75])
    axA.plot([vals.min(),vals.max()],[baseline]*2,color='#B2B7BE',lw=.7,zorder=1)
    axA.plot([q1,q3],[baseline]*2,color='#D1D5DA',lw=5,zorder=1,solid_capstyle='butt')
    yy=np.full(20,baseline,dtype=float)
    for value in np.unique(vals):
        ix=np.flatnonzero(vals==value);yy[ix]+=(np.arange(len(ix))-(len(ix)-1)/2)*.055
    axA.scatter(vals,yy,s=10,color='#929AA4',linewidths=0,zorder=3)
    ix=int(np.flatnonzero(rand.seed.to_numpy()==2001)[0]);axA.scatter([vals[ix]],[yy[ix]],s=25,facecolors='none',edgecolors='#242A31',linewidths=.85,zorder=4)
    axA.plot([p.selected_count],[y+.30],marker='D',ms=4.5,color=color,ls='',zorder=5)
    axA.text(p.selected_count+.8,y+.30,f'{p.selected_count}',fontsize=7,va='center',color=color)
    yes=bool(p.DA_eligible)
    axB.add_patch(Rectangle((.055,y-.125),.35,.25,facecolor=color if yes else 'white',edgecolor=color if yes else '#B2B7BE',lw=.7))
    axB.text(.23,y,'Yes' if yes else 'No',ha='center',va='center',fontsize=7.5,color='white' if yes else '#555555',fontweight='bold')
    count=int(p.random_eligible_seeds)
    axB.text(.77,y+.04,f'{count}/20',ha='center',va='center',fontsize=8)
    axB.add_patch(Rectangle((.60,y-.17),.34,.035,facecolor='#ECEEF1',edgecolor='none'))
    if count:axB.add_patch(Rectangle((.60,y-.17),.34*count/20,.035,facecolor='#929AA4',edgecolor='none'))
    for method,offset,col,letter in [('production',.17,color,'P'),('random',-.17,'#929AA4','R')]:
        s=data[(data.atlas==a)&(data.cell_type==typ)&(data.method==method)].iloc[0]
        axC.text(-.065,y+offset,letter,transform=axC.get_yaxis_transform(),ha='right',va='center',fontsize=7,color=col)
        if s.DA_eligible:
            value=s.strict_recovery_display_percent
            axC.barh(y+offset,value,height=.17,color=col,linewidth=0)
            if value==0:axC.plot(0,y+offset,'|',color=col,ms=7,clip_on=False)
            axC.text(value+2.5,y+offset,f'{value:.1f}%',fontsize=7,va='center',color='#30343A')
        else:
            axC.text(3,y+offset,'DA-ineligible',fontsize=7,va='center',color='#666666',fontstyle='italic')
    if i<3:
        for ax in [labelax,axA,axB,axC]:ax.axhline(y-.50,color='#EBEDF0',lw=.55,zorder=0)

fig.legend(handles=[Line2D([],[],marker='D',ls='',color='#555555',ms=4,label='Proposed (atlas color; P)'),Line2D([],[],marker='o',ls='',color='#929AA4',ms=3,label='Uniform random (R)'),Line2D([],[],marker='o',ls='',mfc='none',mec='#242A31',ms=4,label='Seed 2001 highlighted in A')],loc='lower center',bbox_to_anchor=(.51,.075),ncol=3,frameon=False,fontsize=7,handletextpad=.4,columnspacing=1)
fig.text(.5,.043,'Representation and eligibility: 20 random seeds. DA recovery: one random seed (2001).',ha='center',fontsize=7)
fig.text(.5,.017,'DA-ineligible is not a measured 0% recovery; population examples are descriptive.',ha='center',fontsize=7,color='#555555')
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[]
for ax in fig.axes:
    texts.extend(t for t in ax.texts if t.get_visible() and t.get_text())
    if ax.axison:
        texts.extend(t for t in [ax.xaxis.label,*ax.get_xticklabels()] if t.get_visible() and t.get_text())
texts+=list(fig.texts)
for legend in fig.legends:texts.extend(legend.get_texts())
for t in texts:
    bbox=t.get_window_extent(renderer)
    assert bbox.x0>=0 and bbox.y0>=0 and bbox.x1<=fig.bbox.width and bbox.y1<=fig.bbox.height,('Clipped text',t.get_text())
for i,t in enumerate(texts):
    b=t.get_window_extent(renderer)
    for u in texts[i+1:]:
        assert not b.overlaps(u.get_window_extent(renderer)),('Text overlap',t.get_text(),u.get_text())
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'main_figures/Figure_1_rare_population_representation_and_DA_fidelity.{ext}',dpi=600)
plt.close(fig)
pdf=OUT/'main_figures/Figure_1_rare_population_representation_and_DA_fidelity.pdf';raw=pdf.read_bytes()
m=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',raw)
assert m and abs((float(m[3])-float(m[1]))/72*25.4-178)<.01
assert b'/FontFile2' in raw and b'/Subtype /Type3' not in raw
after={str(p.relative_to(ROOT)):sha(p) for p in frozen if p.is_file()};assert before==after
pd.DataFrame([dict(path=p,sha256=h,unchanged=True) for p,h in before.items()]).to_csv(OUT/'figure_tables/Figure_1_sources/frozen_supplementary_manifest.tsv',sep='\t',index=False)
(OUT/'figure_tables/Figure_1_sources/figure_QA.json').write_text(json.dumps(dict(width_mm=178,PNG_dpi=600,random_representation_seeds=20,random_DA_seed=2001,budget_percent=2,ineligible_displayed_as_zero=False,frozen_supplementary_unchanged=True,embedded_TrueType=True,text_bounds_and_overlaps_checked=True),indent=2))
print('PASS: final figure exported; source values checked; finalized supplementary figures/captions unchanged',flush=True)
