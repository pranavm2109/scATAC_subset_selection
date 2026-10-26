"""Question-led figures from measured summaries. No model fitting or cell selection."""
from common import *
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap


def grouped(ax,h,column,methods):
    for j,m in enumerate(methods):
        vals=[h[(h.atlas==a)&(h.method==m)].iloc[0][column] for a in ATLASES]
        ax.bar(np.arange(3)+(j-(len(methods)-1)/2)*.24,vals,width=.22,color=MCOL[m],label=MLAB[m])
    ax.set_xticks(range(3),[LABELS[a] for a in ATLASES])


def s1():
    h=load('S1_normalization_2pct');g=load('S1_normalization_budgets')
    fig,axs=publication_grid(2,2,figsize=(178/25.4,160/25.4));a,b,c,d=axs.ravel()
    hh=h.copy();hh['retained_pct']=100*hh.n_cell_types_retained/hh.n_cell_types_full
    grouped(a,hh,'retained_pct',METHODS);a.set_ylim(0,108);a.set_ylabel('Annotated types retained (%)');panel(a,'A','Population retention at 2%')
    for i,atlas in enumerate(ATLASES):
        for j,m in enumerate(METHODS):
            r=h[(h.atlas==atlas)&(h.method==m)].iloc[0]
            a.text(i+(j-1)*.24,100*r.n_cell_types_retained/r.n_cell_types_full+1,f'{r.n_cell_types_retained:.0f}',ha='center',fontsize=7)
    a.legend(ncol=1,loc='lower right',frameon=False)
    grouped(b,h,'composition_total_variation',['baseline','l1']);b.set_ylabel('Composition TV distance');b.set_ylim(0,.28);panel(b,'B','Composition deviation at 2%')
    for j,m in enumerate(['baseline','l1']):
        for quant,marker,ls in [('q90','o','-'),('q99','s','--')]:
            vals=[h[(h.atlas==a)&(h.method==m)].iloc[0][f'umap_nearest_selected_{quant}'] for a in ATLASES]
            c.plot(np.arange(3)+(j-.5)*.1,vals,marker=marker,ls='None',color=MCOL[m],label=f'{MLAB[m]}, {quant[1:]}th')
    c.set_xticks(range(3),[LABELS[a] for a in ATLASES]);c.set_yscale('log');c.set_ylabel('Historical 2D UMAP distance');panel(c,'C','Nearest-selected-cell distances\nin 2D UMAP (2%)');c.legend(frameon=False,fontsize=7,loc='upper right')
    for atlas in ['kidney_atac','fetal_atac']:
        for m,ls in [('baseline','--'),('l1','-')]:
            t=g[(g.atlas==atlas)&(g.method==m)].sort_values('pct')
            d.plot(t.pct,t.composition_total_variation,ls+'o',color=COLORS[atlas],label=f'{LABELS[atlas]}, {MLAB[m]}')
    budgets(d);d.set_ylabel('Composition TV distance');panel(d,'D','Budget series: Kidney and Fetal');d.legend(frameon=False,fontsize=7)
    save(fig,'S1_normalization')
    caption('S1_normalization','''# S1 — technical working caption
Question: how does equalizing total RBF feature mass change selection?
178 mm width. A: historical 2% selections, annotated-type retention; integers above bars are types retained, with atlas totals 24/101/78. B: total variation (TV) of annotated-type proportions. C: 90th/99th percentiles of Euclidean distance from every cell to its nearest selected cell in the historical 2D UMAP embedding; selected cells contribute zero. All methods within each historical comparison use one shared embedding, verified in the analysis implementation. These are not 60D PeakVI distances and are not asserted to use byte-identical regenerated production UMAP coordinates. D: existing Kidney/Fetal budget series; HECA multi-budget comparison summary absent.
Sources: diagnostics/rbf_normalization_ablation/{atlas}/comparison/method_summary.tsv and diagnostics/rbf_normalization_generalization/{kidney_atac,fetal_atac}/fraction_method_summary.tsv. Derivatives: tables/S1_*.tsv. Source implementations: diagnostics/summarize_rbf_normalization_ablation.py and summarize_l1_normalization_generalization.py. One recorded selection per method/budget, no error bars or random-baseline inference. Same identities verified against current production wherever historical L1 lists exist; see S1_historical_subset_identity.tsv. Raw unnormalized objective values are not compared as independent validation.''')


def retention(ax,c):
    mat=np.stack([c[c.atlas==a].sort_values('pct').retention_fraction for a in ATLASES])
    ax.imshow(mat,vmin=0,vmax=1,cmap=ListedColormap(['#F1F3F5','#D9EAF2']),aspect='auto')
    for i,a in enumerate(ATLASES):
        t=c[c.atlas==a].sort_values('pct')
        for j,(_,r) in enumerate(t.iterrows()):ax.text(j,i,f'{100*r.retention_fraction:.0f}',ha='center',va='center',fontsize=7)
    ax.set_xticks(range(5),PCTS);ax.set_yticks(range(3),[f'{LABELS[a]} ({int(c[c.atlas==a].n_types.iloc[0])})' for a in ATLASES]);ax.set_xlabel('Selected cells (%)');ax.tick_params(length=0)
    for spine in ax.spines.values():spine.set_visible(False)


def s2():
    c=load('S2_composition');ct=load('S2_all_population_metrics')
    fig,axs=publication_grid(3,3,figsize=(178/25.4,218/25.4))
    eligibility=ct.groupby(['atlas','pct']).DA_eligible.agg(['sum','count']).reset_index()
    ax=axs[0,0]
    for i,atlas in enumerate(ATLASES):
        for j,pct in enumerate(PCTS):
            r=eligibility[(eligibility.atlas==atlas)&(eligibility.pct==pct)].iloc[0]
            ax.scatter(j,i,s=610,marker='s',color=COLORS[atlas],alpha=.13)
            ax.text(j,i,f'{100*r["sum"]/r["count"]:.1f}',ha='center',va='center',fontsize=7)
    ax.set_xlim(-.5,4.5);ax.set_ylim(2.5,-.5);ax.set_xticks(range(5),PCTS);ax.set_yticks(range(3),[LABELS[a] for a in ATLASES]);ax.set_xlabel('Subset (%)');ax.tick_params(length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    panel(ax,'A','DA-eligible types (%)');table('S2_DA_eligibility',eligibility)
    for ax,col,letter,title,yl in [(axs[0,1],'TV','B','Proportion deviation','TV distance'),(axs[0,2],'JS_bits','C','Proportion divergence','JS divergence (bits)')]:
        atlas_lines(ax,c,col);ax.set_ylabel(yl);panel(ax,letter,title)
    axs[0,1].legend(frameon=False,fontsize=6)
    for j,a in enumerate(ATLASES):
        ax=axs[1,j]
        for pct,col,mark in [(2,COLORS[a],'o'),(10,COLORS[a],'^')]:
            t=ct[(ct.atlas==a)&(ct.pct==pct)]
            for eligible in [True,False]:
                u=t[t.DA_eligible==eligible]
                ax.scatter(100*u.full_proportion,u.log2_enrichment,s=13,marker=mark,facecolors=col if eligible else 'none',edgecolors=col,linewidths=.65,alpha=.8)
        ax.axhline(0,color='#777777',lw=.65);ax.set_xscale('log');ax.set_xlabel('Full-atlas abundance (%)');ax.set_ylabel('Log2 proportion enrichment');ax.set_ylim(-4.2,4.2);panel(ax,chr(68+j),LABELS[a])
        if j==0:
            ax.legend(handles=[Line2D([],[],marker='o',color='#555555',ls='',label='2%'),Line2D([],[],marker='^',color='#555555',ls='',label='10%'),Line2D([],[],marker='o',mfc='none',color='#666666',ls='',label='DA-ineligible')],frameon=False,fontsize=7,loc='lower right')
    examples={'kidney_atac':[('Mes','Mes'),('Macula_Densa','Macula densa')],'heca_atac':[('Haematopoietic stem and progenitor cell','Stem/progenitor'),('Naive CD4 T cell','Naive CD4 T')],'fetal_atac':[('Megakaryocytes?','Megakaryocytes?'),('Ganglion cells','Ganglion')]}
    rare=[]
    for j,a in enumerate(ATLASES):
        ax=axs[2,j]
        for k,(name,label) in enumerate(examples[a]):
            t=ct[(ct.atlas==a)&(ct.cell_type==name)].sort_values('pct');rare.append(t)
            ax.plot(t.pct,t.subset_count,ls=['-','--'][k],marker=['o','s'][k],color=COLORS[a],label=label)
            # Hollow rings explicitly indicate failure of the complete target/rest eligibility rule.
            u=t[~t.DA_eligible];ax.scatter(u.pct,u.subset_count,s=35,facecolors='white',edgecolors=COLORS[a],zorder=4)
        ax.text(.98,10,'10 cells',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=7,color='#555555')
        ax.axhline(10,color='#777777',ls=':',lw=.8);ax.set_yscale('log');ax.set_ylim(.7,100);budgets(ax);ax.set_ylabel('Selected population cells');panel(ax,chr(71+j),f'{LABELS[a]}: rare examples');ax.legend(frameon=False,fontsize=6,loc='upper left')
    table('S2_rare_examples',pd.concat(rare));save(fig,'S2_population_representation')
    caption('S2_population_representation','''# S2 — technical working caption
Question: are populations present, proportionally represented, and sufficiently sampled for DA?
178 mm width. A: percentage of annotated cell types retained at each independent subset size; row labels give total types in parentheses. B/C: TV and Jensen–Shannon divergence of subset versus full-atlas proportions, computed from existing counts; JS uses base-2 logarithms (bits), not square-root JS distance. D–F: every population at 2% and 10%; enrichment = subset proportion / full-atlas proportion. Hollow markers mean <10 target cells or <20 rest cells, cross-checked against saved DA status tables. The zero line is equal proportions. G–I: deliberately selected rare-population examples, not an unbiased population summary; hollow rings mark DA ineligibility and the dotted line is 10 target cells. Full names/counts are in S2_rare_examples.tsv. Eosinophilic granulocytes and all other types remain in S2_all_population_metrics.tsv.
Sources: results/cell_type_proportions/peakvi/{atlas}/cell_type_proportions_peakvi.tsv and results/signac_da/{atlas}/{fraction}/signac_da_status.tsv. New summary calculations only; no new selections. Presence does not imply composition fidelity, DA eligibility, or superiority over random sampling. Examples are descriptive and should not be presented as a random sample.''')


def s3():
    p=load('S3_pooled');d=load('S3_per_cell_type_direction');f=load('S3_fixed_eligibility');dr=load('S3_direction_summary')
    fig,axs=publication_grid(3,3,figsize=(178/25.4,218/25.4))
    for ax,y,letter,title,yl,pct in [(axs[0,0],'direction_concordance','A','Effect direction','Agreement (%)',True),(axs[0,1],'strict_recovery_fraction','B','Significant-peak recovery','Strict recovery (%)',True),(axs[0,2],'spearman_rho','C','Pooled signed effects','Spearman rho',False)]:
        atlas_lines(ax,p,y,pct);panel(ax,letter,title);ax.set_ylabel(yl)
    axs[0,0].set_ylim(70,101);axs[0,1].set_ylim(0,80);axs[0,2].set_ylim(.5,1);axs[0,0].legend(frameon=False,fontsize=6)
    ax=axs[1,0]
    for a in ATLASES:
        t=p[p.atlas==a].sort_values('pct');ax.plot(t.pct,100*t.fraction_effects_available,'-o',color=COLORS[a]);ax.plot(t.pct,100*t.n_evaluable_cell_types/t.n_reference_cell_types,'-s',color=COLORS[a])
    ax.set_ylim(80,101);budgets(ax);ax.set_ylabel('Available / total (%)');panel(ax,'D','Changing eligibility')
    ax.legend(handles=[Line2D([],[],color='#555555',marker='o',label='Reference pairs'),Line2D([],[],color='#555555',marker='s',ls='-',label='Cell types')],frameon=False,fontsize=6)
    for ax,y,letter,title,yl in [(axs[1,1],'direction_agreement','E','Common eligible types','Direction agreement (%)'),(axs[1,2],'strict_recovery','F','Common eligible types','Strict recovery (%)')]:
        for a in ATLASES:
            for cohort,ls in [('common_eligible','--'),('all_reference','-')]:
                t=f[(f.atlas==a)&(f.cohort==cohort)].sort_values('pct');ax.plot(t.pct,100*t[y],ls+'o',color=COLORS[a],markersize=2.5)
        budgets(ax);ax.set_ylabel(yl);panel(ax,letter,title)
    axs[1,1].legend(handles=[Line2D([],[],color='#555555',ls='--',label='Common cohort'),Line2D([],[],color='#555555',label='All reference')],frameon=False,fontsize=7)
    for j,a in enumerate(ATLASES):
        ax=axs[2,j]
        for drn,offset in [('open',-.22),('closed',.22)]:
            vals=[d[(d.atlas==a)&(d.pct==pct)&(d.direction==drn)].effect_spearman_rho.dropna().to_numpy() for pct in PCTS]
            b=ax.boxplot(vals,positions=np.array(PCTS)+offset,widths=.36,patch_artist=True,showfliers=False,manage_ticks=False,medianprops={'color':'black','lw':.8},boxprops={'lw':.65},whiskerprops={'lw':.65},capprops={'lw':.65})
            for box in b['boxes']:box.set_facecolor(DIRCOL[drn]);box.set_alpha(.65)
        ax.axhline(0,color='#888888',lw=.6);ax.set_ylim(-.35,1.05);budgets(ax);ax.set_ylabel('Per-type Spearman rho');panel(ax,chr(71+j),LABELS[a])
    axs[2,0].legend(handles=[Line2D([],[],color=DIRCOL[x],lw=4,label=x.capitalize()) for x in ['open','closed']],frameon=False,fontsize=6,loc='lower right')
    for ax in axs[:2,:].ravel():ax.set_xlabel('')
    save(fig,'S3_DA_fidelity')
    fig,axs=publication_grid(1,3,figsize=(178/25.4,80/25.4))
    for j,a in enumerate(ATLASES):
        ax=axs[j]
        for drn in ['open','closed']:
            for cohort,ls in [('all_reference','-'),('common_eligible','--')]:
                t=dr[(dr.atlas==a)&(dr.direction==drn)&(dr.cohort==cohort)].sort_values('pct');ax.plot(t.pct,100*t.strict_recovery,ls+'o',color=DIRCOL[drn],markersize=2.5)
        ax.set_ylim(0,100);budgets(ax);ax.set_ylabel('Strict recovery (%)');panel(ax,chr(65+j),LABELS[a])
    axs[0].legend(handles=[Line2D([],[],color=DIRCOL[x],label=x.capitalize()) for x in ['open','closed']]+[Line2D([],[],color='#555555',label='All reference'),Line2D([],[],color='#555555',ls='--',label='Common cohort')],frameon=False,fontsize=7)
    save(fig,'S3b_direction_specific_recovery')
    # Explicit figure-specific source table for the companion.
    table('S3b_direction_specific_recovery',dr)
    caption('S3_DA_fidelity','''# S3 — technical working caption
Question: how do DA effect fidelity, statistical recovery, and eligibility change with budget?
178 mm width. A: matching nonzero effect signs / available finite subset effects at full-atlas significant peak–cell-type pairs. B: pairs significant in the subset with matching direction / all significant full-atlas reference pairs; ineligible types remain in this denominator. C: existing pooled Spearman across signed effects and cell types. D: fraction of reference pairs with available effects (solid) and fraction of reference cell types eligible (dashed). E/F: exact count-based metrics restricted to types eligible at every budget, alongside all-reference curves; fixed cohorts are listed in S3_common_eligible_types.tsv. G–I: per-type correlations separately within full-atlas open and closed reference peaks. Boxes show median, quartiles and 1.5-IQR whiskers; outliers are not drawn, no populations are removed from tables, and nonfinite correlations are excluded from boxes. These are across-population distributions, not confidence intervals or seed replicates. Small reference-peak sets can produce unstable correlations.
Fixed-common-eligibility pooled Spearman cannot be recovered from per-type correlations; it is intentionally not fabricated. Raw peak-level processing would be needed. Count sums were checked against the production publication summary. Zero effects are neither matching nor switched. Full-atlas significance: BH adjusted P <=0.05, |average log2FC| >=0.25 and detection >=1% in either group; subset strict significance uses the same thresholds. DA eligibility: >=10 target and >=20 rest cells. BH correction occurs within cell type. Effect reproduction does not require subset significance.
Sources: results/publication_panels/signac_da_recapitulation_3x5_summary.tsv; results/signac_da_capture/{atlas}/da_peak_capture_summary.tsv. Derived sources S3_*.tsv; original definitions scripts/09_1_run_signac_da_all_cell_types.R and 09_2_compute_signac_da_peak_capture.R. Subsets overlap their parent atlas. Changes can reflect sample size, eligibility, and altered one-versus-rest composition; power alone is not established as the cause.''')
    caption('S3b_direction_specific_recovery','''# S3b — technical working caption
178 mm width companion to S3. Same-direction statistically significant recovery separately for full-atlas open and closed reference peaks. Solid: all reference types; dashed: types eligible at all five budgets. Denominator is the direction-specific full-atlas reference-pair count in the relevant cohort, including pairs lacking subset effects for the all-reference cohort. This plot aggregates existing counts and does not rerun DA. See S3 caption for thresholds and limitations. Source: tables/S3b_direction_specific_recovery.tsv, derived from existing per-cell-type capture summaries.''')


def s5():
    p=load('S5_selection');s=load('S5_stage_costs')
    fig=plt.figure(figsize=(178/25.4,145/25.4),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1.4,.23,1.25]);a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[2,:]);legend_band=fig.add_subplot(gs[1,0]);legend_band.set_axis_off();a._publication_legend_band=legend_band
    for atlas in ATLASES:
        t=p[p.atlas==atlas].sort_values('pct');a.plot(t.n_selected,t.fit_hours,'-o',color=COLORS[atlas],label=LABELS[atlas]);b.plot(t.pct,t.peak_CPU_RSS_GiB,'-o',color=COLORS[atlas],label=LABELS[atlas])
    a.set_xscale('log');a.set_yscale('log');a.set_xlabel('Selected cells (log scale)');a.set_ylabel('Apricot fit time (h, log scale)');panel(a,'A','Subset-selection cost');a.legend(frameon=False)
    budgets(b);b.set_ylabel('Selection-job peak CPU RSS (GiB)');b.set_ylim(0,29);panel(b,'B','Measured selection memory')
    c.axis('off');panel(c,'C','Pipeline stages: single successful production runs');c.texts[0].set_x(-.03)
    stages=['Conversion','PeakVI','Landmark features','Selection (10%)'];rows=[]
    for atlas in ATLASES:
        rr=[]
        for st in stages:
            v=s[(s.atlas==atlas)&(s.stage==st)].iloc[0].job_seconds
            rr.append(f'{v/60:.1f} min' if st in ['Conversion','Landmark features'] else f'{v/3600:.2f} h')
        rows.append(rr)
    tab=c.table(cellText=rows,rowLabels=[LABELS[x] for x in ATLASES],colLabels=['Conversion /\nexport job','PeakVI\ntraining job','Landmark\nfeature job','Selection job\n(10%)'],loc='center',cellLoc='center',bbox=[.12,.16,.87,.75]);tab.auto_set_font_size(False);tab.set_fontsize(8)
    c.text(.12,.02,'Table: job wall time. Panel A: Apricot fit() time only.',transform=c.transAxes,fontsize=7)
    for (r,col),cell in tab.get_celld().items():
        cell.set_edgecolor('#D8DDE2');cell.set_linewidth(.5)
        if r==0:cell.set_facecolor('#EAF0F4')
    save(fig,'S5_computational_cost')
    caption('S5_computational_cost','''# S5 — technical working caption
178 mm width. A: recorded Apricot fit() seconds converted to hours, excluding input loading/normalization and all upstream training. B: measured successful selection-job peak CPU RSS in GiB (benchmark MiB /1024), not requested memory or GPU memory. C: successful Snakemake job wall times for conversion/export, PeakVI training, landmark-feature construction, and the independent 10% selection. Stage times are not combined into a misleading selection-only total. Conversion is not the earlier biological preprocessing. All five production selections are independent jobs, not prefix extraction.
Sources: results/subset_selection/peakvi/{atlas}/subset_selection_summary_peakvi_{fraction}.tsv; logs/benchmarking/subset_selection/peakvi/{atlas}_{fraction}.tsv; logs/benchmarking/{anndata,peakvi}/{atlas}.tsv; logs/benchmarking/subset_selection/peakvi_landmarks/{atlas}.tsv. Figure data S5_selection.tsv and S5_stage_costs.tsv. Individual successful runs, one seed (123); no timing uncertainty, controlled scaling comparison, or comparative speed claim. Different training epoch budgets and hardware/load can affect cross-atlas times. GPU memory, historical failed-run expenditure and donor preprocessing costs are not measured here.''')


def main_candidate():
    h=load('MAIN_normalization');c=load('MAIN_composition');p=load('MAIN_DA')
    fig,axs=publication_grid(2,2,figsize=(178/25.4,140/25.4));a,b,d,e=axs.ravel()
    hh=h.copy();hh['retained_pct']=100*hh.n_cell_types_retained/hh.n_cell_types_full
    grouped(a,hh,'retained_pct',['baseline','l1']);a.set_ylim(0,105);a.set_ylabel('Annotated types retained (%)');panel(a,'A','Normalization improves retention (2%)');a.legend(frameon=False,fontsize=6.5,loc='lower right')
    atlas_lines(b,c,'TV');b.set_ylabel('Composition TV distance');panel(b,'B','Population proportions remain altered');b.legend(frameon=False,fontsize=6.5)
    atlas_lines(d,p,'direction_concordance',True);d.set_ylabel('Direction agreement (%)');d.set_ylim(70,100);panel(d,'C','DA effect reproduction')
    atlas_lines(e,p,'strict_recovery_fraction',True);e.set_ylabel('Strict recovery (%)');e.set_ylim(0,80);panel(e,'D','DA significance and direction')
    save(fig,'MAIN_candidate_quantitative')
    caption('MAIN_candidate_quantitative','''# Main quantitative candidate — technical working caption
178 mm, four panels. A: historical matched-budget unnormalized versus L1 RBF ablation at 2%, retaining 22/24 versus 24/24 Kidney types, 90/101 versus 101/101 HECA, and 74/78 versus 78/78 Fetal. B: current production TV of cell-type proportions. C: matching effect direction among available subset effects at full-atlas significant peak–cell-type pairs. D: subset significance with matching direction divided by all full-atlas reference pairs. Budgets independently selected; no annotations guide production selection. Curves are measured individual runs, not replicate averages. C has budget-dependent eligibility; D retains the full reference denominator. Counts, thresholds, and eligibility limitations are documented in S1–S3. No random-sampling superiority or quantitative 60D coverage claim is made. Sources: tables/MAIN_*.tsv, with original inputs recorded in SOURCE_MANIFEST.tsv. Recommend as a main-figure candidate only if quantitative results take priority over the separate workflow illustration.''')

if __name__=='__main__':
    for fun in [s1,s2,s3,s5]:
        print('Plotting',fun.__name__,flush=True);fun()
