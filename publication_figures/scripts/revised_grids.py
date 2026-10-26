"""Display-only changes using saved coordinates, memberships and density bins."""
from common import *
from matplotlib.colors import LogNorm

def umap_grid():
    fig=plt.figure(figsize=(178/25.4,155/25.4),layout='constrained')
    gs=fig.add_gridspec(6,5,height_ratios=[1,.21]*3)
    summary=read('results/publication_panels/peakvi_subset_overlay_3x5_summary.tsv');rows=[]
    for i,a in enumerate(ATLASES):
        frame=read(f'results/umap_peakvi/{a}/{a}_peakvi_umap.tsv',dtype={'cell_name':str}).set_index('cell_name')
        assert not frame.index.has_duplicates
        xy=frame[['UMAP_1','UMAP_2']].to_numpy(float);assert np.isfinite(xy).all()
        bg=np.random.default_rng(123).choice(len(frame),min(75000,len(frame)),replace=False)
        lo=xy.min(axis=0);hi=xy.max(axis=0);pad=(hi-lo)*.035
        for j,f in enumerate(FRACS):
            ids=pd.read_csv(ROOT/f'results/subset_selection/peakvi/{a}/selected_cells_peakvi_{f}.txt',header=None,dtype=str)[0]
            assert not ids.duplicated().any()
            chosen=frame.loc[ids,['UMAP_1','UMAP_2']].to_numpy()
            r=summary[(summary.atlas==a)&(summary.fraction==f)].iloc[0]
            assert len(chosen)==r.n_selected_cells and len(frame)==r.n_atlas_cells
            show=np.random.default_rng(123).choice(len(chosen),min(3000,len(chosen)),replace=False)
            ax=fig.add_subplot(gs[2*i,j]);band=fig.add_subplot(gs[2*i+1,j]);band.axis('off')
            ax.scatter(xy[bg,0],xy[bg,1],s=.24,color='#B8BEC5',alpha=.48,linewidths=0,rasterized=True)
            ax.scatter(chosen[show,0],chosen[show,1],s=.65,color=COLORS[a],alpha=.62,linewidths=0,rasterized=True)
            ax.set_xlim(lo[0]-pad[0],hi[0]+pad[0]);ax.set_ylim(lo[1]-pad[1],hi[1]+pad[1]);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
            for spine in ax.spines.values():spine.set_visible(False)
            ax.set_title(f'{chr(65+i*5+j)}   {PCTS[j]}%',loc='left',fontsize=8,fontweight='bold',pad=3)
            if j==0:ax.set_ylabel(LABELS[a],color=COLORS[a],fontweight='bold',fontsize=8)
            band.text(.5,.5,f'n = {len(chosen):,}\nshown = {len(show):,}',ha='center',va='center',fontsize=7)
            rows.append(dict(atlas=a,subset_size=f,n_atlas=len(frame),n_selected=len(chosen),n_selected_shown=len(show),n_background_shown=len(bg),display_seed=123))
    fig.supxlabel('Grey: full-atlas background (up to 75,000 shown); color: selected cells (up to 3,000 shown)',fontsize=7)
    save(fig,'S7_descriptive_UMAP_grid');table('S7_display_counts',pd.DataFrame(rows))
    caption('S7_descriptive_UMAP_grid','''# S7 — saved production UMAPs
Rows: Kidney (blue), HECA (orange), Fetal (green). Columns: independently selected 2–10% budgets. Existing embedding coordinates, independent selection memberships and shared within-atlas limits are unchanged. Each panel reports the actual selected count n and display count shown. For legibility only, background and selected points are uniformly sampled without replacement, capped at 75,000 and 3,000 respectively using NumPy default_rng(123), separately for background and each selection in saved row order. This display sampling is not an experimental baseline or a change to selection. Point density is not a quantitative coverage statistic. Smaller markers and a fixed selected display cap prevent large selections from hiding the atlas background. Rasterized points; vector labels.
Sources: results/umap_peakvi/{atlas}/{atlas}_peakvi_umap.tsv; results/subset_selection/peakvi/{atlas}/selected_cells_peakvi_{fraction}.txt. Actual counts checked against results/publication_panels/peakvi_subset_overlay_3x5_summary.tsv. Display counts: polished/tables/S7_display_counts.tsv. Width 178 mm.''')

def da_grid():
    p=load('S3_pooled');b=load('S8_density_bins');hist={}
    for a in ATLASES:
        for f in FRACS:
            t=b[(b.atlas==a)&(b.subset_size==f)];h=np.zeros((80,80),dtype=np.int64)
            assert not t.duplicated(['x_bin','y_bin']).any()
            h[t.x_bin.astype(int),t.y_bin.astype(int)]=t.n_reference_pairs
            assert h.sum()==p[(p.atlas==a)&(p.subset_size==f)].iloc[0].n_effects_available
            hist[a,f]=h
    fig=plt.figure(figsize=(178/25.4,151/25.4),layout='constrained')
    gs=fig.add_gridspec(8,5,height_ratios=[1,.13,1,.13,1,.13,.10,.13])
    norm=LogNorm(vmin=1,vmax=max(h.max() for h in hist.values()))
    for i,a in enumerate(ATLASES):
        lim=float(p[p.atlas==a].display_axis_limit.iloc[0])
        for j,f in enumerate(FRACS):
            ax=fig.add_subplot(gs[2*i,j]);band=fig.add_subplot(gs[2*i+1,j]);band.axis('off')
            r=p[(p.atlas==a)&(p.subset_size==f)].iloc[0]
            im=ax.imshow(np.ma.masked_equal(hist[a,f].T,0),origin='lower',extent=[-lim,lim,-lim,lim],cmap='cividis',norm=norm,aspect='equal',rasterized=True)
            ax.plot([-lim,lim],[-lim,lim],color='#9D563C',ls='--',lw=.7)
            ax.set_xticks([-round(lim*.7,1),0,round(lim*.7,1)]);ax.set_yticks([-round(lim*.7,1),0,round(lim*.7,1)]);ax.tick_params(labelsize=7,length=2,pad=2)
            ax.set_title(f'{chr(65+i*5+j)}   {PCTS[j]}%',loc='left',fontsize=8,fontweight='bold',pad=3)
            if j==0:ax.set_ylabel(LABELS[a],color=COLORS[a],fontweight='bold',fontsize=8)
            else:ax.tick_params(labelleft=False)
            if i<2:ax.tick_params(labelbottom=False)
            band.text(.5,.5,f'rho = {r.spearman_rho:.2f}',ha='center',va='center',fontsize=7)
    cax=fig.add_subplot(gs[6,1:4]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.ax.tick_params(labelsize=7,length=2,pad=1)
    cb.set_label('Reference peak–cell-type pairs / bin (log scale)',fontsize=7,labelpad=2)
    note=fig.add_subplot(gs[7,:]);note.axis('off')
    note.text(.5,.3,'Reference filter: |full-atlas log2FC| ≥ 0.25; display-only axis clipping',ha='center',va='center',fontsize=7)
    fig.supxlabel('Full-atlas average log2 fold change',fontsize=8)
    fig.supylabel('Subset average log2 fold change',fontsize=8)
    save(fig,'S8_DA_effect_density_grid')
    caption('S8_DA_effect_density_grid','''# S8 — existing DA effect densities
Each 80 × 80 square-bin panel counts significant full-atlas reference peak–cell-type pairs with finite subset effects. Saved density bins are reused exactly; no testing or peak-level recomputation. Shared logarithmic color scale across all panels. Dashed diagonal: exact agreement. Correlations below panels are the existing pooled Spearman correlations of unmodified signed effects across cell types, not correlations computed from bins. They are distinct from per-cell-type, within-direction correlations in S3.
The gap near zero on the full-atlas axis reflects the reference eligibility filter |full-atlas average log2FC| ≥ 0.25; it is not a biological discontinuity or missing-data artifact. Bin width can partially obscure this gap. Within-atlas axis limits are shared across budgets and copied from production. Values beyond these limits were clipped to boundary bins for display only; boundary accumulation is not an effect mode. Correlations use unclipped effects. Each histogram sum is checked against the available-effect denominator. Availability changes with cell-type DA eligibility (S3). Full-reference significance and overlapping parent/subset limitations remain as documented in S3.
Sources: publication_figures/tables/S8_density_bins.tsv and S3_pooled.tsv, originally derived from results/signac_da_capture/{atlas}/da_peak_overlap_details.tsv and results/publication_panels/signac_da_recapitulation_3x5_summary.tsv. Width 178 mm. Density rasterized, text vector.''')

if __name__=='__main__':
    umap_grid();da_grid()
