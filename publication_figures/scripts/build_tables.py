"""Read small existing summaries; write only new publication_figures derivatives."""
from common import *
import hashlib, importlib.metadata

def main():
    source_rows=[]
    def src(path,figure):
        p=ROOT/path
        source_rows.append(dict(figure=figure,source=str(path),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    hist=[];general=[];cells=[];comp=[];da=[];perf=[];stage=[];identity=[]
    for a in ATLASES:
        path=f'diagnostics/rbf_normalization_ablation/{a}/comparison/method_summary.tsv';src(path,'S1');h=read(path);hist.append(h)
        if a!='heca_atac':
            path=f'diagnostics/rbf_normalization_generalization/{a}/fraction_method_summary.tsv';src(path,'S1');general.append(read(path))
        path=f'results/cell_type_proportions/peakvi/{a}/cell_type_proportions_peakvi.tsv';src(path,'S2');d=read(path)
        N=int(d.atlas_count.sum());p=d.atlas_count/N
        for f,pct in zip(FRACS,PCTS):
            path=f'results/signac_da/{a}/{f}/signac_da_status.tsv';src(path,'S2,S3');st=read(path).set_index('cell_type')
            c=d[f'subset_{f}_count'];n=int(c.sum());q=c/n;r=q/p;tv,js=composition(p,q)
            comp.append(dict(atlas=a,subset_size=f,pct=pct,n_atlas=N,n_selected=n,n_types=len(d),n_retained=int((c>0).sum()),retention_fraction=float((c>0).mean()),TV=tv,JS_bits=js,full_count_min=int(d.atlas_count.min()),full_count_median=d.atlas_count.median(),full_count_max=int(d.atlas_count.max()),subset_count_min=int(c.min()),subset_count_median=c.median(),subset_count_max=int(c.max())))
            for i,row in d.iterrows():
                t=row.cell_type;eligible=bool(st.loc[t,'status']=='success');expected=bool(c.iloc[i]>=10 and n-c.iloc[i]>=20)
                assert eligible==expected and int(st.loc[t,'n_target_cells'])==int(c.iloc[i])
                cells.append(dict(atlas=a,cell_type=t,subset_size=f,pct=pct,full_count=int(row.atlas_count),full_proportion=float(p.iloc[i]),subset_count=int(c.iloc[i]),subset_proportion=float(q.iloc[i]),enrichment=float(r.iloc[i]),log2_enrichment=float(np.log2(r.iloc[i])) if r.iloc[i]>0 else np.nan,DA_eligible=eligible))
            path=f'results/subset_selection/peakvi/{a}/subset_selection_summary_peakvi_{f}.tsv';src(path,'S5');s=read(path).iloc[0]
            path=f'logs/benchmarking/subset_selection/peakvi/{a}_{f}.tsv';src(path,'S5');b=read(path).iloc[-1]
            perf.append(dict(atlas=a,subset_size=f,pct=pct,n_selected=int(s.n_selected_unique),fit_seconds=float(s.time_seconds),fit_hours=float(s.time_seconds)/3600,job_seconds=float(b['s']),peak_CPU_RSS_GiB=float(b.max_rss)/1024,seed=int(s.random_seed_selection),optimizer=s.selection_optimizer))
            old=ROOT/f'diagnostics/rbf_normalization_ablation/{a}/l1/selected_cells_{f}.txt';new=ROOT/f'results/subset_selection/peakvi/{a}/selected_cells_peakvi_{f}.txt'
            if old.exists():
                u=set(old.read_text().splitlines());v=set(new.read_text().splitlines());identity.append(dict(atlas=a,subset_size=f,identical_cell_sets=u==v,jaccard=len(u&v)/len(u|v),historical=str(old.relative_to(ROOT)),production=str(new.relative_to(ROOT))))
        path=f'results/signac_da_capture/{a}/da_peak_capture_summary.tsv';src(path,'S3');da.append(read(path))
        for name,rel in [('Conversion',f'logs/benchmarking/anndata/{a}.tsv'),('PeakVI',f'logs/benchmarking/peakvi/{a}.tsv'),('Landmark features',f'logs/benchmarking/subset_selection/peakvi_landmarks/{a}.tsv'),('Selection (10%)',f'logs/benchmarking/subset_selection/peakvi/{a}_10pct.tsv')]:
            src(rel,'S5');b=read(rel).iloc[-1];stage.append(dict(atlas=a,stage=name,job_seconds=float(b['s']),job_hours=float(b['s'])/3600,peak_CPU_RSS_GiB=float(b.max_rss)/1024))
    h=pd.concat(hist);g=pd.concat(general);g['pct']=g.fraction_label.str[:2].astype(int)
    c=pd.DataFrame(comp);ct=pd.DataFrame(cells);d=pd.concat(da,ignore_index=True);d['pct']=d.subset_size.str[:2].astype(int)
    assert len(c)==15 and len(ct)==1015 and len(d)==2030
    path='results/publication_panels/signac_da_recapitulation_3x5_summary.tsv';src(path,'S3');pooled=read(path);pooled['pct']=pooled.subset_size.str[:2].astype(int)
    fixed=[];direction=[];commonrows=[]
    for a in ATLASES:
        ad=d[d.atlas==a];elig=ad[['cell_type','subset_size','subset_analysis_status']].drop_duplicates();eligible=elig[elig.subset_analysis_status=='success'].groupby('cell_type').subset_size.nunique();common=set(eligible[eligible==5].index)
        for t in sorted(common):commonrows.append(dict(atlas=a,cell_type=t))
        for f,pct in zip(FRACS,PCTS):
            allrows=ad[ad.subset_size==f]
            for cohort,t in [('all_reference',allrows),('common_eligible',allrows[allrows.cell_type.isin(common)])]:
                nr=int(t.n_da_peaks_atlas.sum());na=int(t.n_reference_effects_available.sum());same=int(t.n_same_direction_effect.sum());strict=int(t.n_same_direction_overlap.sum());sw=int(t.n_sign_switched.sum())
                fixed.append(dict(atlas=a,subset_size=f,pct=pct,cohort=cohort,n_cell_types=t.cell_type.nunique(),n_reference=nr,n_available=na,direction_agreement=same/na,strict_recovery=strict/nr,n_sign_switched=sw,availability=na/nr))
                if cohort=='all_reference':
                    original=pooled[(pooled.atlas==a)&(pooled.subset_size==f)].iloc[0]
                    assert nr==original.n_atlas_reference_peak_cell_type_pairs and na==original.n_effects_available and sw==original.n_sign_switched
                    assert np.isclose(same/na,original.direction_concordance) and np.isclose(strict/nr,original.strict_recovery_fraction)
                else:assert nr==na
                for dr in ['open','closed']:
                    u=t[t.direction==dr];nr=int(u.n_da_peaks_atlas.sum());na=int(u.n_reference_effects_available.sum())
                    direction.append(dict(atlas=a,subset_size=f,pct=pct,cohort=cohort,direction=dr,n_reference=nr,n_available=na,direction_agreement=u.n_same_direction_effect.sum()/na,strict_recovery=u.n_same_direction_overlap.sum()/nr,median_per_type_spearman=u.effect_spearman_rho.median(),n_finite_correlations=u.effect_spearman_rho.notna().sum()))
    da_cost=[]
    for atlas in ATLASES:
        for pct in PCTS:
            rel=f'logs/benchmarking/signac_da/{atlas}/{pct:02d}pct.tsv';src(rel,'S6 plan');b=read(rel).iloc[-1]
            da_cost.append(dict(atlas=atlas,pct=pct,observed_production_wall_seconds=float(b['s']),observed_peak_CPU_RSS_GiB=float(b.max_rss)/1024,source=rel))
    table('S6_DA_cost_reference',pd.DataFrame(da_cost))
    table('S1_normalization_2pct',h);table('S1_normalization_budgets',g);table('S1_historical_subset_identity',pd.DataFrame(identity))
    table('S2_composition',c);table('S2_all_population_metrics',ct)
    table('S3_per_cell_type_direction',d);table('S3_pooled',pooled);table('S3_fixed_eligibility',pd.DataFrame(fixed));table('S3_direction_summary',pd.DataFrame(direction));table('S3_common_eligible_types',pd.DataFrame(commonrows))
    table('S5_selection',pd.DataFrame(perf));table('S5_stage_costs',pd.DataFrame(stage));table('SOURCE_MANIFEST',pd.DataFrame(source_rows))
    # Figure-specific tables for the recommended composite.
    table('MAIN_normalization',h[h.method.isin(['baseline','l1'])]);table('MAIN_composition',c);table('MAIN_DA',pooled)
    (OUT/'tables'/'software_versions.json').write_text(json.dumps({p:importlib.metadata.version(p) for p in ['numpy','pandas','matplotlib','scipy']},indent=2))
    print('Validated: 15 budgets, 1,015 population rows, 2,030 directional DA rows; count aggregates reproduce publication metrics.',flush=True)
if __name__=='__main__':main()
