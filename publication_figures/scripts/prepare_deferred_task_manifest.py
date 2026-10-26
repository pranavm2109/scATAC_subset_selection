"""Writes task metadata only: never selects random cells or reads latent matrices."""
from common import *
rows=[]
for a in ATLASES:
    for pct in PCTS:
        for method,seed in [('production',123)]+[('random',s) for s in range(2001,2021)]:
            rows.append(dict(task_index=len(rows),atlas=a,pct=pct,method=method,seed=seed))
d=pd.DataFrame(rows);table('S4_S6_DEFERRED_task_manifest',d)
pilot=d[(d.pct.isin([2,10]))&((d.method=='production')|((d.method=='random')&(d.seed==2001)))]
print('Full task count:',len(d),'(15 production + 300 random). Pilot array indices:',','.join(pilot.task_index.astype(str)))
