from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'publication_figures'
ATLASES=['kidney_atac','heca_atac','fetal_atac']
LABELS={'kidney_atac':'Kidney','heca_atac':'HECA','fetal_atac':'Fetal'}
COLORS={'kidney_atac':'#0072B2','heca_atac':'#D55E00','fetal_atac':'#009E73'}
FRACS=['02pct','04pct','06pct','08pct','10pct']
PCTS=[2,4,6,8,10]
METHODS=['baseline','l1','singleton_capped']
MCOL={'baseline':'#6D7785','l1':'#0072B2','singleton_capped':'#CC79A7'}
MLAB={'baseline':'Unnormalized','l1':'L1 normalized','singleton_capped':'Singleton capped'}
DIRCOL={'open':'#0072B2','closed':'#D55E00'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.25,'lines.markersize':3,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def read(rel,**kw):return pd.read_csv(ROOT/rel,sep='\t',**kw)
def table(name,df):
    p=OUT/'tables'/f'{name}.tsv';df.to_csv(p,sep='\t',index=False);return p

def load(name):return pd.read_csv(ROOT/'publication_figures'/'tables'/f'{name}.tsv',sep='\t')
def publication_grid(nrows,ncols,figsize):
    """Reserve a separate legend band under each plot row; never cover data."""
    fig=plt.figure(figsize=figsize,layout='constrained')
    gs=fig.add_gridspec(2*nrows,ncols,height_ratios=[1,.22]*nrows)
    fig._publication_grid=(gs,nrows,ncols)
    axs=np.empty((nrows,ncols),dtype=object)
    for i in range(nrows):
        for j in range(ncols):
            ax=fig.add_subplot(gs[2*i,j]);band=fig.add_subplot(gs[2*i+1,j]);band.set_axis_off()
            ax._publication_legend_band=band;axs[i,j]=ax
    return fig,axs.squeeze()

def save(fig,name):
    # Migrate every panel legend to its allocated, non-data legend band.
    moved=[]
    for ax in list(fig.axes):
        if not ax.axison:continue  # Dedicated legend bands are not source panels.
        old=ax.get_legend()
        if old is None:continue
        target=getattr(ax,'_publication_legend_band',None)
        if target is None:raise RuntimeError(f'{name}: panel legend lacks a dedicated band')
        handles=list(old.legend_handles);labels=[t.get_text() for t in old.get_texts()];old.remove()
        width_mm=fig.get_figwidth()*25.4*ax.get_position().width
        ncol=min(2,len(labels)) if width_mm>=60 else 1
        legend=target.legend(handles,labels,loc='center',ncol=ncol,frameon=False,fontsize=7,
                             borderaxespad=0,borderpad=.2,columnspacing=1.1,
                             handlelength=2.6,handletextpad=.5,labelspacing=.35)
        moved.append((ax,target,legend))
    grid=getattr(fig,'_publication_grid',None)
    if grid is not None:
        gs,nrows,ncols=grid;ratios=[1,.22]*nrows
        for row in range(nrows):
            if not any(ax.get_subplotspec().rowspan.start==2*row for ax,_,_ in moved):ratios[2*row+1]=.015
        gs.set_height_ratios(ratios)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();bounds=fig.bbox
    report=[]
    # Axis-label checks catch crowding that legend-only tests miss.
    data_axes=[ax for ax in fig.axes if ax.axison]
    for i,ax in enumerate(data_axes):
        for label in [ax.xaxis.label,ax.yaxis.label]:
            if not label.get_text():continue
            box=label.get_window_extent(renderer)
            for other in data_axes[i+1:]:
                for other_label in [other.xaxis.label,other.yaxis.label]:
                    if other_label.get_text():assert not box.overlaps(other_label.get_window_extent(renderer)),(name,'axis labels overlap')
    for ax,target,legend in moved:
        box=legend.get_window_extent(renderer)
        assert bounds.contains(box.x0,box.y0) and bounds.contains(box.x1,box.y1),(name,'legend clipped')
        for other in fig.axes:
            if other is target or not other.axison:continue
            # Include tick labels, axis labels and panel titles, not just data rectangles.
            other_box=other.get_tightbbox(renderer)
            assert not box.overlaps(other_box),(name,'legend overlaps plotted panel or labels')
        report.append(dict(figure=name,panel_title=ax.get_title(loc='left'),n_entries=len(legend.get_texts()),legend_font_pt=7,legend_outside_all_data_and_label_bounds=True))
    if report:table(name+'_legend_QA',pd.DataFrame(report))
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/'figures'/f'{name}.{ext}',dpi=600)
    plt.close(fig)
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=8)
    ax.text(-.14,1.075,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
    ax.tick_params(width=.6,length=3)
def budgets(ax):ax.set_xticks(PCTS);ax.set_xlabel('Subset (%)')
def atlas_lines(ax,df,y,percent=False,**kwargs):
    for a in ATLASES:
        t=df[df.atlas==a].sort_values('pct')
        ax.plot(t.pct,t[y]*(100 if percent else 1),'-o',color=COLORS[a],label=LABELS[a],**kwargs)
    budgets(ax)
def caption(name,text): (OUT/'captions'/f'{name}.md').write_text(text.strip()+'\n')
def composition(p,q):
    p=np.asarray(p,float);q=np.asarray(q,float);m=(p+q)/2
    js=0.
    for x in [p,q]:
        keep=x>0;js+=.5*np.sum(x[keep]*np.log2(x[keep]/m[keep]))
    return float(.5*np.abs(p-q).sum()),float(js)
