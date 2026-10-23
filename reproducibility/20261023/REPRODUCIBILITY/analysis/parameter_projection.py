import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,solve_triangular
from scipy.integrate import cumulative_trapezoid
r=np.load(Path(__file__).with_name('result_arrays.npz'))
C=r['C']; z=r['z'];zh=r['zh'];om=float(r['om_hat']); L=cholesky(C,lower=True)
w=solve_triangular(L,np.ones(len(z)),lower=True); a=w/np.linalg.norm(w)
def white(v):
 w=solve_triangular(L,v,lower=True)
 return w-a*(a@w)
def mu(omega):
 g=np.linspace(0,max(z)*1.001,8192);v=1/np.sqrt(omega*(1+g)**3+1-omega);i=cumulative_trapezoid(v,g,initial=0)
 return 5*np.log10((299792.458/70)*(1+zh)*np.interp(z,g,i))+25
u=white((mu(om+1e-5)-mu(om-1e-5))/2e-5)
s=white((z>.36).astype(float))
corr=float((u@s)/np.sqrt((u@u)*(s@s)))
print('corr',corr,'se inflation',1/np.sqrt(1-corr*corr))
