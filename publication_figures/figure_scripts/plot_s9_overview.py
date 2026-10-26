"""S9: unchanged four-panel quantitative overview; saved summary tables only."""
from pathlib import Path
import os,sys,shutil,tempfile
assert os.environ.get('SLURM_JOB_ID'), 'Submit the plotting-only batch entrypoint'
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'publication_figures/scripts'))
import common
from common import *
from plot_existing_summaries import grouped
def load(name):return pd.read_csv(ROOT/'publication_figures/figure_tables'/f'{name}.tsv',sep='\t')
_legacy_save=save
def save(fig,name):
    """Preserve original layout checks; export to the final folders only."""
    previous=common.OUT
    with tempfile.TemporaryDirectory(prefix='s9_plot_',dir=ROOT/'publication_figures/reports') as folder:
        common.OUT=Path(folder)
        try:
            for sub in ['figures','tables']:(common.OUT/sub).mkdir()
            _legacy_save(fig,name)
            for ext in ['pdf','svg','png']:
                shutil.copy2(common.OUT/'figures'/f'{name}.{ext}',ROOT/'publication_figures/supplementary_figures'/f'{name}.{ext}')
            for p in (common.OUT/'tables').glob('*.tsv'):
                shutil.copy2(p,ROOT/'publication_figures/figure_tables'/p.name)
        finally:common.OUT=previous

def quantitative_overview():
    h=load('S9_normalization');c=load('S9_composition');p=load('S9_DA')
    fig,axs=publication_grid(2,2,figsize=(178/25.4,140/25.4));a,b,d,e=axs.ravel()
    hh=h.copy();hh['retained_pct']=100*hh.n_cell_types_retained/hh.n_cell_types_full
    grouped(a,hh,'retained_pct',['baseline','l1']);a.set_ylim(0,105);a.set_ylabel('Annotated types retained (%)');panel(a,'A','Normalization improves retention (2%)');a.legend(frameon=False,fontsize=6.5,loc='lower right')
    atlas_lines(b,c,'TV');b.set_ylabel('Composition TV distance');panel(b,'B','Population proportions remain altered');b.legend(frameon=False,fontsize=6.5)
    atlas_lines(d,p,'direction_concordance',True);d.set_ylabel('Direction agreement (%)');d.set_ylim(70,100);panel(d,'C','DA effect reproduction')
    atlas_lines(e,p,'strict_recovery_fraction',True);e.set_ylabel('Strict recovery (%)');e.set_ylim(0,80);panel(e,'D','DA significance and direction')
    save(fig,'S9_quantitative_selection_and_biological_fidelity_overview')

if __name__=='__main__':quantitative_overview()
