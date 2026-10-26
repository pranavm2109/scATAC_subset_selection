from common import *
import os
assert os.environ.get('SLURM_JOB_ID')
d=load('S3_per_cell_type_direction');p=load('S3_pooled');f=load('S3_fixed_eligibility');ct=load('S2_all_population_metrics');members=load('S3_common_eligible_types');checks=[]
for a in ATLASES:
    cohorts=[]
    for frac,pct in zip(FRACS,PCTS):
        status=read(f'results/signac_da/{a}/{frac}/signac_da_status.tsv')
        eligible=set(status.loc[(status.n_target_cells>=10)&(status.n_rest_cells>=20),'cell_type'])
        cc=ct[(ct.atlas==a)&(ct.pct==pct)]
        assert set(cc.loc[cc.DA_eligible,'cell_type'])==eligible
        cohorts.append(eligible)
        dd=d[(d.atlas==a)&(d.subset_size==frac)];pp=p[(p.atlas==a)&(p.subset_size==frac)].iloc[0]
        assert dd.n_da_peaks_atlas.sum()==pp.n_atlas_reference_peak_cell_type_pairs
        assert dd.n_reference_effects_available.sum()==pp.n_effects_available
        assert dd.n_sign_switched.sum()==pp.n_sign_switched
        assert dd.n_same_direction_overlap.sum()==pp.n_strict_same_direction_recovered
        assert np.isclose(dd.n_same_direction_effect.sum()/pp.n_effects_available,pp.direction_concordance)
        assert np.isclose(dd.n_sign_switched.sum()/pp.n_effects_available,pp.sign_switch_fraction)
        checks.append(dict(atlas=a,subset_size=frac,n_eligible=len(eligible),n_types=len(status),n_reference=int(pp.n_atlas_reference_peak_cell_type_pairs),n_available=int(pp.n_effects_available),denominators_verified=True))
    common=set.intersection(*cohorts)
    assert common==set(members.loc[members.atlas==a,'cell_type'])
    for frac in FRACS:
        dd=d[(d.atlas==a)&(d.subset_size==frac)&d.cell_type.isin(common)]
        ff=f[(f.atlas==a)&(f.subset_size==frac)&(f.cohort=='common_eligible')].iloc[0]
        assert ff.n_cell_types==len(common)
        assert dd.n_da_peaks_atlas.sum()==ff.n_reference
        assert dd.n_reference_effects_available.sum()==ff.n_available
        assert np.isclose(dd.n_same_direction_effect.sum()/ff.n_available,ff.direction_agreement)
        assert np.isclose(dd.n_same_direction_overlap.sum()/ff.n_reference,ff.strict_recovery)
table('eligibility_denominator_QA',pd.DataFrame(checks))

identity=load("S1_historical_subset_identity")
for r in identity.itertuples():
    assert set(pd.read_csv(ROOT/r.historical,header=None,dtype=str)[0])==set(pd.read_csv(ROOT/r.production,header=None,dtype=str)[0])
print("PASS: cohort, denominator and historical membership checks")
