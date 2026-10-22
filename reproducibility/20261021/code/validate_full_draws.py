"""Direct raw data-space generation and noninterpolated GL96 optimizer checks.
This is an executor validation, distinct from the reviewer's independent code.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,time
from scipy.optimize import minimize_scalar
from scipy.linalg import solve_triangular
from numpy.polynomial.legendre import leggauss
from independent_step import Engine,ROOT,save_json
start=time.monotonic();e=Engine(all_partitions=True)
a=np.load(ROOT/'outputs/baseline/draws.npz');n=len(a['omega'])
selected=set(range(8))
for key in ['primary_T','joint_T','all_partitions_T','omega','simultaneous_bin_T']:
    ix=np.argsort(a[key]);selected.update(ix[:3].tolist());selected.update(ix[-3:].tolist())
for i in np.linspace(0,n-1,40,dtype=int):
    if len(selected)>=40:break
    selected.add(int(i))
selected=sorted(selected)[:40];needed=set(selected)
rng=np.random.default_rng(2026102101);eps={}
for i in range(n):
    z=rng.standard_normal(e.n)
    if i in needed:eps[i]=z
q,w=leggauss(96);zz=e.z[:,None]*(q+1)/2

def mu(om):
    integral=e.z/2*np.sum(w/np.sqrt(om*(1+zz)**3+1-om),axis=1)
    return 5*np.log10(299792.458/70*(1+e.zh)*integral)+25

def fit(fun):
    opt=minimize_scalar(fun,bounds=e.bounds,method='bounded',options={'xatol':1e-11})
    if not opt.success:raise RuntimeError(opt.message)
    options=[(opt.fun,opt.x)]+[(fun(x),x) for x in e.bounds]
    f,x=min(options);return x,f
results=[]
for ix in selected:
    mock=e.fit_mean+e.L@eps[ix] # materialized full correlated vector
    yw=e.project(e.white(mock))
    def r(om):return yw-e.project(e.white(mu(om)))
    def fun(om):res=r(om);return float(res@res)
    om,c=fit(fun);rr=r(om);sc=rr@e.Q[:,:e.nstep];T=float(np.max(sc*sc));sa=rr@e.all_Q
    joints=[]
    for k in range(e.nstep):
        def jf(o):
            rr=r(o);score=rr@e.Q[:,k];return float(rr@rr-score**2)
        oj,cj=fit(jf);joints.append((c-cj,oj))
    jm=max(x[0] for x in joints);piv=e.pivots[np.argmax([x[0] for x in joints])]
    bins=e.bin_operator@rr;bt=float(np.max((bins/e.bin_sigma)**2))
    x={'index':ix,'Omega_error':float(abs(om-a['omega'][ix])),'chi2_error':float(abs(c-a['chi2'][ix])),'primary_T_error':float(abs(T-a['primary_T'][ix])),'all_partitions_T_error':float(abs(np.max(sa*sa)-a['all_partitions_T'][ix])),'joint_T_error':float(abs(jm-a['joint_T'][ix])),'joint_pivot':float(piv),'stored_joint_pivot':float(a['joint_pivot'][ix]),'simultaneous_bin_T_error':float(abs(bt-a['simultaneous_bin_T'][ix]))}
    results.append(x);save_json(ROOT/'outputs/direct40_validation.json',{'status':'RUNNING','completed':len(results),'selected':selected,'checks':results,'seconds':time.monotonic()-start})
    print(ix,round(time.monotonic()-start,2),max(x[k] for k in x if k.endswith('_error')),flush=True)
maxima={k:max(x[k] for x in results) for k in results[0] if k.endswith('_error')}
save_json(ROOT/'outputs/direct40_validation.json',{'status':'COMPLETE','N':len(results),'selected':selected,'checks':results,'max_errors':maxima,'joint_pivot_mismatches':sum(x['joint_pivot']!=x['stored_joint_pivot'] for x in results),'pass':max(maxima.values())<1e-3,'seconds':time.monotonic()-start,'generator':'Materialized y*=fitted_mean+L epsilon; full correlated data-space model, GL96 distances and direct nuisance projection; no Chebyshev calls in objective'})
