"""Small regression check on saved real vectors; no benchmark output generated."""
from common import ROOT
import numpy as np
from pynndescent import NNDescent
from deferred_benchmarks import nearest,validate_index
m=np.load(ROOT/'results/peakvi/kidney_atac/peakvi_latent.npy',mmap_mode='r')
x=m[:512];assert not x.flags.writeable
before=x.copy();buf=np.array(x,dtype=np.float32,order='C',copy=True)
assert buf.flags.writeable and buf.flags.c_contiguous
index=NNDescent(buf,metric='euclidean',n_neighbors=30,random_state=1701,n_jobs=1,low_memory=True);index.prepare()
d=nearest(index,buf,x[:16],k=4,exclude_positions=np.arange(16))
assert np.isfinite(d).all() and (d>=0).all()
print('Full-reference/nonself check:',validate_index(index,buf,x,None,1702,.2,True),flush=True)
chosen=np.arange(0,512,4);ref=np.array(x[chosen],dtype=np.float32,order='C',copy=True)
sub=NNDescent(ref,metric='euclidean',n_neighbors=30,random_state=1705,n_jobs=1,low_memory=True);sub.prepare()
print('Subset-reference check:',validate_index(sub,ref,x,chosen,1706,.2),flush=True)
assert np.array_equal(before,x) and not x.flags.writeable
print('PASS: writable index copies and read-only query inputs work; original vectors unchanged.',flush=True)
