from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]
a=np.load(R/'analysis/result_arrays.npz')
z=a['z'];zh=a['zh'];y=a['y'];C=a['C'];n=len(z)
d=np.genfromtxt(R/'inputs/Pantheon+SH0ES.dat',names=True,dtype=None,encoding='utf-8')
ids=d['IDSURVEY'][a['selected_idx']]
cf=cho_factor(C,lower=True,check_finite=False)
K=cho_solve(cf,np.eye(n),check_finite=False)
b=((z>=1)&(z<1.3)).astype(float)
X1=np.ones((n,1)); vals=np.unique(ids); X20=np.column_stack([ids==v for v in vals]).astype(float)
def mu(om):
 g=np.linspace(0,max(z)*1.001,4096)
 integ=cumulative_trapezoid(1/np.sqrt(om*(1+g)**3+1-om),g,initial=0)
 return 5*np.log10((299792.458/70)*(1+zh)*np.interp(z,g,integ))+25
for label,X in [('one',X1),('survey20',X20)]:
 for sv in (100,101):
  dr=np.flatnonzero(ids==sv);kee=np.flatnonzero(ids!=sv)
  ker=K[np.ix_(kee,kee)] - K[np.ix_(kee,dr)]@np.linalg.solve(K[np.ix_(dr,dr)],K[np.ix_(dr,kee)])
  xx=X[kee]; keepcol=(xx.sum(axis=0)>0);xx=xx[:,keepcol]
  bb=b[kee];xxk=xx.T@ker; S=xxk@xx;si=np.linalg.inv(S)
  Pb=ker@bb-(ker@xx)@(si@(xxk@bb));v=bb@Pb
  def objective(om):
   rr=(y-mu(om))[kee];Kr=ker@rr
   return rr@Kr-(xx.T@Kr)@(si@(xx.T@Kr))
  fit=minimize_scalar(objective,bounds=(.15,.5),method='bounded',options={'xatol':5e-8})
  T=((Pb@(y-mu(fit.x))[kee])**2/v)
  print(label,'drop survey',sv,'removed',len(dr),'remain in high bin',int(bb.sum()),'T',round(T,5),'om',round(fit.x,5))
