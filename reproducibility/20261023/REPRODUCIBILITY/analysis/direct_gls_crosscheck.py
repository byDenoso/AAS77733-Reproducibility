"""Direct original-data Cholesky GLS crosscheck of the 20-survey reduced model."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np, json
from pathlib import Path
from scipy.linalg import cholesky, solve_triangular
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'analysis'
a=np.load(A/'result_arrays.npz')
expected=json.loads((A/'survey20_summary.json').read_text())
d=np.genfromtxt(ROOT/'inputs'/'Pantheon+SH0ES.dat',names=True,dtype=None,encoding='utf-8')
ids=np.asarray(d['IDSURVEY'][a['selected_idx']],dtype=int); vals=np.unique(ids)
L=cholesky(a['C'],lower=True,check_finite=False)
w=lambda x:solve_triangular(L,x,lower=True,check_finite=False)
G=np.column_stack([ids==v for v in vals]).astype(float)
X=w(G);basis=w(a['residual'])
B=CubicSpline(a['om_grid'],a['B'],axis=1)
best_om=float(expected['best_omega_m'])
basis=basis-a['U']@B(best_om)
fit=lambda M:np.linalg.lstsq(M,basis,rcond=None)[1][0]
no=fit(X)
step_pivot=float(expected['step_search']['pivot'])
step=w((a['z']>step_pivot).astype(float));withstep=fit(np.column_stack([X,step]))
lo,hi=expected['bin_family']['max_bin']; high=w(((a['z']>=lo)&(a['z']<hi)).astype(float));withbin=fit(np.column_stack([X,high]))
obs_step=float(no-withstep); obs_bin=float(no-withbin)
res={'n_surveys':int(len(vals)),'null_chi2_direct':float(no),'step_pivot':step_pivot,
     'step_local_delta_chi2_direct':obs_step,'high_bin': [float(lo),float(hi)],
     'high_bin_local_delta_chi2_direct':obs_bin,
     'step_vs_reduced_difference':float(obs_step-float(expected['step_search']['T'])),
     'bin_vs_reduced_difference':float(obs_bin-float(expected['bin_family']['T']))}
(A/'direct_gls_check.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
