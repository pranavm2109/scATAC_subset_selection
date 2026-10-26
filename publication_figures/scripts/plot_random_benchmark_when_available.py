"""Deferred S4/S6 plotter. Refuses incomplete measured benchmarks; no placeholders."""
from common import *
import argparse


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark-dir',required=True);args=ap.parse_args()
    paths=list(Path(args.benchmark_dir).glob('*_metrics.tsv'))
    if not paths:raise SystemExit('No measured benchmark files; no S4/S6 figure generated.')
    d=pd.concat([pd.read_csv(p,sep='\t') for p in paths],ignore_index=True);d=d[d.cohort=='all_queries']
    if d.duplicated(['atlas','subset_size','method','seed']).any():raise ValueError('Duplicate benchmark runs')
    for a in ATLASES:
        for f in FRACS:
            r=d[(d.atlas==a)&(d.subset_size==f)&(d.method=='random')];p=d[(d.atlas==a)&(d.subset_size==f)&(d.method=='production')]
            if r.seed.nunique()<20 or len(p)!=1:raise ValueError(f'Incomplete measured comparison: {a} {f}')
            if r.query_mode.nunique()!=1 or p.query_mode.iloc[0]!=r.query_mode.iloc[0] or r.query_seed.nunique()!=1 or p.query_seed.iloc[0]!=r.query_seed.iloc[0]:raise ValueError('Incomparable query protocols')
            if any(r.n_query!=p.n_query.iloc[0]):raise ValueError('Incomparable query sizes')
            if r.fixed_radius.isna().any() or p.fixed_radius.isna().any() or not np.allclose(r.fixed_radius,p.fixed_radius.iloc[0]):raise ValueError('Missing or incomparable fixed radius')
    if not d.groupby('atlas').fixed_radius.nunique().eq(1).all():raise ValueError('Radius must remain fixed across budgets within atlas')
    summaries=[]
    for metric in ['n_retained','TV','JS_bits','median_distance','q90_distance','q99_distance','coverage_within_radius']:
        for (a,f,m),u in d.groupby(['atlas','subset_size','method']):
            x=u[metric].dropna();summaries.append(dict(atlas=a,subset_size=f,method=m,metric=metric,n_seeds=len(x),median=x.median(),seed_q025=x.quantile(.025),seed_q975=x.quantile(.975)))
    s=pd.DataFrame(summaries);s['pct']=s.subset_size.str[:2].astype(int);table('S4_S6_measured_benchmark_summary',s)
    for name,metrics in [('S4_measured_60D_coverage',['median_distance','q90_distance','q99_distance','coverage_within_radius']),('S6_measured_random_comparison',['n_retained','TV','JS_bits'])]:
        fig,axs=publication_grid(len(metrics),3,figsize=(178/25.4,(72*len(metrics))/25.4))
        for i,metric in enumerate(metrics):
            for j,a in enumerate(ATLASES):
                ax=axs[i,j];r=s[(s.atlas==a)&(s.metric==metric)&(s.method=='random')].sort_values('pct');p=s[(s.atlas==a)&(s.metric==metric)&(s.method=='production')].sort_values('pct')
                ax.fill_between(r.pct,r.seed_q025,r.seed_q975,color='#B9C1C9',alpha=.55);ax.plot(r.pct,r['median'],'--',color='#666666',label='Random median');ax.plot(p.pct,p['median'],'-o',color=COLORS[a],label='Proposed');budgets(ax);ax.set_ylabel({'median_distance':'Median nearest distance','q90_distance':'90th-percentile distance','q99_distance':'99th-percentile distance','coverage_within_radius':'Fraction within fixed radius','n_retained':'Annotated types retained','TV':'Composition TV distance','JS_bits':'JS divergence (bits)'}[metric]);panel(ax,chr(65+i*3+j),LABELS[a])
        axs[0,0].legend(frameon=False);save(fig,name)
    population_paths=list(Path(args.benchmark_dir).glob('*_populations.tsv'))
    if not population_paths:raise ValueError('Missing per-population random benchmark outputs')
    pops=pd.concat([pd.read_csv(p,sep='\t') for p in population_paths],ignore_index=True)
    examples={'kidney_atac':['Mes','Macula_Densa'],'heca_atac':['Haematopoietic stem and progenitor cell','Naive CD4 T cell'],'fetal_atac':['Megakaryocytes?','Ganglion cells']}
    short={'Haematopoietic stem and progenitor cell':'Stem/progenitor','Naive CD4 T cell':'Naive CD4 T','Macula_Densa':'Macula densa'}
    selected=[];fig,axs=publication_grid(2,3,figsize=(178/25.4,144/25.4))
    for j,a in enumerate(ATLASES):
        for i,t in enumerate(examples[a]):
            u=pops[(pops.atlas==a)&(pops.cell_type==t)].copy();u['pct']=u.subset_size.str[:2].astype(int);selected.append(u)
            if u[u.method=='random'].groupby('pct').seed.nunique().min()<20:raise ValueError('Incomplete rare-cell replicates')
            r=u[u.method=='random'].groupby('pct').subset_count.agg(median='median',q025=lambda v:v.quantile(.025),q975=lambda v:v.quantile(.975));p=u[u.method=='production'].sort_values('pct')
            ax=axs[i,j];ax.fill_between(r.index,r.q025,r.q975,color='#B9C1C9',alpha=.55);ax.plot(r.index,r['median'],'--',color='#666666');ax.plot(p.pct,p.subset_count,'-o',color=COLORS[a]);ax.axhline(10,ls=':',lw=.7,color='#777777');budgets(ax);ax.set_ylabel('Selected population cells');ax.set_ylim(bottom=0);panel(ax,chr(65+i*3+j),f'{LABELS[a]}: {short.get(t,t)}')
    save(fig,'S6b_measured_rare_population_counts');table('S6b_measured_rare_population_counts',pd.concat(selected))
    caption('S4_S6_measured_benchmarks','''# Deferred measured benchmark caption template — review after execution
Distances use original saved 60D Euclidean PeakVI vectors with calibrated approximate neighbor search; inspect validation and query_mode fields before interpretation. Shaded bands are empirical 2.5th–97.5th percentiles across independently generated random selection seeds, NOT confidence intervals. Proposed method has one production selection per budget and no seed uncertainty estimate. All methods share the declared evaluation-query protocol. All-query distances include selected cells at zero; nonselected-query summaries are retained separately in benchmark metrics. Rare-cell panels use deliberately selected examples; dotted lines mark 10 target cells, with full rest-cell eligibility in the population tables. No DA fidelity claims follow from these count/coverage plots. Actual results must be reviewed before drafting conclusions.''')
if __name__=='__main__':main()
