from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import cho_factor, cho_solve
r=np.load(Path(__file__).with_name('result_arrays.npz'))
idx=r['selected_idx']; d=np.genfromtxt(str(Path(os.environ.get('PANTHEON_DATA_DIR',Path(__file__).resolve().parents[1]/'inputs'))/'Pantheon+SH0ES.dat'),names=True,dtype=None,encoding='utf8')
C=r['C']; n=len(idx); zz=r['z']; residual=r['residual']; yc=r['y']; zmax=max(zz)
cf=cho_factor(C,lower=True,check_finite=False)
ones=np.ones(n); ci1=cho_solve(cf,ones,check_finite=False)
inter=(ci1@residual)/(ci1@ones)
r0=residual-inter
cid=d['CID'][idx]
print('DATA SHAPE',n,'Unique',len(set(cid)), 'zminmax',zz.min(),zz.max(),'intercept',inter,'rresid low/high',np.quantile(r0,[0,.01,.5,.99,1]));
bins=np.array([.0233,.05,.10,.20,.30,.40,.50,.65,.8,1.0,1.3,1.8,2.5])
ix=np.digitize(zz,bins,right=True)-1
X=np.column_stack([ix==j for j in range(len(bins)-1)]).astype(float)
mat=X.T@cho_solve(cf,X,check_finite=False);cov=np.linalg.inv(mat);means=cov@(X.T@cho_solve(cf,r0,check_finite=False));errs=np.sqrt(np.diag(cov));counts=np.sum(X,axis=0).astype(int)
for i in range(len(bins)-1): print('bin',i,bins[i],bins[i+1],counts[i],round(means[i],5),round(errs[i],5))
print('SURVEYS',[(int(k),int(v)) for k,v in zip(*np.unique(d['IDSURVEY'][idx],return_counts=True))])
